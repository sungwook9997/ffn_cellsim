#!/usr/bin/env python3
"""TAG ``gen`` stage: authoritative context, BM25 excerpts, and grounded answer."""
from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from typing import Any

import duckdb
from tag_backend import llm
from tag_exec import DB_PATH

try:
    from supersession import match_topic, resolve_authoritative
except Exception:  # pragma: no cover - optional when imported outside this directory
    match_topic = resolve_authoritative = None

BM25_FLOOR = 6.0

GEN_SYS = (
    "You answer questions over the ffn_cellsim knowledge base. Use the SQL "
    "result rows as the structured ground truth, and the PDF evidence excerpts "
    "(BM25 over the reference corpus) only when they are relevant to the "
    "question. Be concise and precise. Cite the specific identifiers present "
    "(citation_key / KB-x / MC-* / VG-* / DOI) and, for excerpt-based claims, "
    "the [citation_key pN] tag. If the rows are empty and no excerpt applies, "
    "say no matching records were found. Do not invent data.\n"
    "AUTHORITATIVE-CHAIN RULE: when an 'AUTHORITATIVE RECORD' block is provided, "
    "it is the current, latest conclusion for that topic and OVERRIDES any older "
    "or superseded-looking statement in the SQL rows or PDF excerpts. Report it "
    "as the current answer; reference superseded records ONLY as history and say "
    "explicitly that they were superseded. Never present a superseded conclusion "
    "as the current status.\n"
    "RESULTS-AUDIT RULE: the run_audit table is the disk-grounded integrity verdict "
    "for headline result-claims (twin of source_audit for citations). A result is "
    "'done / validated' ONLY if its run_audit verdict is VERIFIED. If a result has "
    "verdict RETRACT (forbidden basal-footprint metric or value drift), NEEDS_REGEN "
    "(claimed artifact absent — number only in prose/PNG), or GPU_UNREPRODUCED "
    "(GPU-only, no committed build/CI trace), state that explicitly and do NOT report "
    "it as established — cite the verdict and its note.\n"
    "PARAMETER-PROVENANCE RULE: the param_audit table is the disk+citation-grounded "
    "integrity verdict for the simulation constants (constant -> KU -> citation -> "
    "verdict). A constant is 'literature-anchored / no magic number' ONLY if its "
    "param_audit verdict is VERIFIED. If verdict is VALUE_DRIFT (config value != its "
    "sourced value), UNSOURCED (no committed KU->source link), SOURCE_SUSPECT "
    "(fabrication-risk citation) or SOURCE_UNVERIFIED (citation only CHECK/NO_DOI), "
    "say so explicitly and do NOT present the constant as fully sourced — cite the "
    "verdict, the KU, and the citation_key."
)


def _load_kb_records(con: Any) -> list[dict[str, Any]]:
    """Load supersession records when the optional derived layer exists."""
    have = con.execute(
        "SELECT count(*) FROM information_schema.tables "
        "WHERE table_name='kb_record'"
    ).fetchone()[0]
    if not have:
        return []
    rows = con.execute(
        "SELECT topic, record_id, claim, status, is_primary, "
        "authoritative_as_of, conclusion, aliases FROM kb_record"
    ).fetchall()
    records = []
    for topic, record_id, claim, status, is_primary, as_of, conclusion, aliases in rows:
        try:
            parsed_aliases = json.loads(aliases) if aliases else []
        except Exception:
            parsed_aliases = []
        records.append(
            {
                "topic": topic,
                "record_id": record_id,
                "claim": claim,
                "status": status,
                "is_primary": is_primary == "true",
                "authoritative_as_of": as_of,
                "conclusion": conclusion,
                "aliases": parsed_aliases,
                "supersedes": [],
                "superseded_by": [],
            }
        )
    return records


def authoritative_context(con: Any, question: str) -> str:
    """Return the newest authoritative record for a matched question topic."""
    if match_topic is None or resolve_authoritative is None:
        return ""
    records = _load_kb_records(con)
    if not records:
        return ""
    topic = match_topic(records, question)
    if not topic:
        return ""
    authoritative = resolve_authoritative(records, topic)
    if not authoritative or authoritative["status"] != "authoritative":
        return ""
    superseded = [
        record
        for record in records
        if record["topic"] == topic and record["status"] == "superseded"
    ]
    lines = [
        f"AUTHORITATIVE RECORD for topic '{topic}' "
        f"(as_of {authoritative['authoritative_as_of']}): "
        f"{authoritative['record_id']}",
        f"  Current conclusion: {authoritative['conclusion']}",
    ]
    if superseded:
        lines.append(
            "  SUPERSEDED (historical — do NOT report as current): "
            + ", ".join(sorted(record["record_id"] for record in superseded))
        )
    lines.append(
        "  When answering about current/authoritative status, use the "
        "AUTHORITATIVE RECORD above; mention superseded findings only as "
        "history and say they were superseded."
    )
    return "\n".join(lines)


def content_search(
    question: str, k: int = 4, floor: float = BM25_FLOOR
) -> list[tuple[Any, ...]]:
    """Return BM25 paper excerpts at or above the measured relevance floor."""
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        have = con.execute(
            "SELECT count(*) FROM information_schema.tables "
            "WHERE table_name='paper_chunks'"
        ).fetchone()[0]
        if not have:
            return []
        return con.execute(
            "SELECT citation_key, page, substr(text,1,400) AS snippet, "
            "fts_main_paper_chunks.match_bm25(chunk_id, ?) AS score "
            "FROM paper_chunks WHERE score IS NOT NULL AND score >= ? "
            "ORDER BY score DESC LIMIT ?",
            [question, floor, k],
        ).fetchall()
    except Exception:
        return []
    finally:
        con.close()


def excerpts_to_text(excerpts: Sequence[Sequence[Any]]) -> str:
    """Format BM25 rows for the grounded-answer prompt."""
    if not excerpts:
        return "(none)"
    return "\n".join(
        f"- [{citation_key} p{page}] {str(snippet).strip()[:300]}"
        for citation_key, page, snippet, _ in excerpts
    )


def rows_to_text(
    cols: Sequence[str], rows: Sequence[Sequence[Any]], cap: int = 60
) -> str:
    """Format SQL rows as a compact pipe-delimited block."""
    head = " | ".join(cols)
    body = "\n".join(
        " | ".join("" if value is None else str(value) for value in row)
        for row in rows[:cap]
    )
    extra = f"\n... ({len(rows) - cap} more rows)" if len(rows) > cap else ""
    return f"{head}\n{body}{extra}"


def gen(
    question: str,
    sql: str,
    cols: Sequence[str],
    rows: Sequence[Sequence[Any]],
    excerpts: Sequence[Sequence[Any]],
    model: str | None,
    auth_context: str = "",
    *,
    llm_fn: Callable[..., str] = llm,
) -> str:
    """Generate a concise answer grounded in rows and relevant excerpts."""
    table = rows_to_text(cols, rows)
    auth_block = (
        "AUTHORITATIVE-CHAIN CONTEXT (overrides older statements):\n"
        f"{auth_context}\n\n"
        if auth_context
        else ""
    )
    prompt = (
        f"Question: {question}\n\n{auth_block}SQL used:\n{sql}\n\n"
        f"Result rows ({len(rows)}):\n{table}\n\n"
        "PDF evidence excerpts (use only if relevant):\n"
        f"{excerpts_to_text(excerpts)}\n\nAnswer:"
    )
    return llm_fn(prompt, system=GEN_SYS, model=model, max_tokens=1200)
