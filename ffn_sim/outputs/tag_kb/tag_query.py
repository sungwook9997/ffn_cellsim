#!/usr/bin/env python3
"""TAG (Table-Augmented Generation) query engine over the ffn_cellsim KB.

Canonical TAG loop (Biswal et al. 2024):
    syn  : NL question R  ->  DuckDB SQL  (Claude, given the live schema + the
           actual relation vocabulary from the edges table)
    exec : SQL  ->  rows  (DuckDB, read-only, SELECT-only guard, 1x self-repair)
    gen  : (R, rows)  ->  NL answer with SE/KB/MC/VG id + DOI citations  (Claude)

Backend tables live in `kb.duckdb` (built by `notion_to_duckdb.py`). The PDF
content layer (`paper_chunks`, FTS) is added by `references_ingest.py`; when
present, `gen` can pull evidence text for the semantic step.

LLM backend (dual, no new secret required):
    • anthropic SDK if ANTHROPIC_API_KEY is set
    • else `claude -p` (Claude Code headless) via subprocess

USAGE:
    python tag_query.py "Which Parameters feed a failing ValidationGate?"
    python tag_query.py --sql-only "How many SourceEvidence papers have no DOI?"
    python tag_query.py --model claude-opus-4-8 "..."
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import subprocess
import sys
import textwrap

import duckdb

# supersession / authoritative-chain helpers (request: latest authoritative wins)
try:
    from supersession import match_topic, resolve_authoritative
except Exception:                                            # pragma: no cover
    match_topic = resolve_authoritative = None

HERE = pathlib.Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"

# columns worth exposing distinct values for (helps syn ground filters)
CATEGORICAL = {"source_type", "status", "confidence", "unit", "subtopic", "kind"}
WRITE_RE = re.compile(r"\b(insert|update|delete|drop|create|alter|attach|copy|"
                      r"pragma|replace|truncate|install|load)\b", re.I)


# --------------------------------------------------------------------------- #
# LLM backend
# --------------------------------------------------------------------------- #
def llm(prompt: str, system: str = "", model: str | None = None,
        max_tokens: int = 1500, temperature: float = 0.0) -> str:
    """One-shot completion. anthropic SDK if key present, else `claude -p`.

    temperature defaults to 0.0: this is a factual KB QA backend, so the NL->SQL
    `syn` step and the `gen` step must be deterministic. Non-deterministic
    sampling (the SDK default 1.0) was the dominant source of TAG flakiness —
    the same question produced different SQL run-to-run. (The `claude -p`
    fallback has no temperature knob; it stays as-is.)
    """
    key = os.environ.get("ANTHROPIC_API_KEY", "").strip()
    if key:
        import anthropic
        client = anthropic.Anthropic(api_key=key)
        msg = client.messages.create(
            model=model or "claude-opus-4-8",
            max_tokens=max_tokens,
            temperature=temperature,
            system=system or None,
            messages=[{"role": "user", "content": prompt}],
        )
        return "".join(b.text for b in msg.content if b.type == "text")
    # fallback: Claude Code headless. Prompt via stdin (robust to length/escaping)
    full = f"{system}\n\n{prompt}" if system else prompt
    cmd = ["claude", "-p"]
    if model:
        cmd += ["--model", model]
    r = subprocess.run(cmd, input=full, text=True, capture_output=True, timeout=180)
    if r.returncode != 0:
        raise RuntimeError(f"claude -p failed: {r.stderr[:400]}")
    return r.stdout.strip()


# --------------------------------------------------------------------------- #
# schema introspection — give syn an accurate world model
# --------------------------------------------------------------------------- #
def schema_text(con) -> str:
    out: list[str] = [
        "TABLES (all node columns are TEXT; CAST when comparing numbers):",
        "KEY IDS — every node's `id` is an OPAQUE Notion-page UUID (e.g. "
        "'372120da…'), used only for JOINs via the edges table. The HUMAN-READABLE "
        "identifier is a separate column: knowledge_claim.kb_id ('KB-3.19'), "
        "source_evidence.citation_key ('Bell1978_Science') and .uid ('SE119'). "
        "To filter by a 'KB-x.y' value use kb_id, NEVER id. Prefer the direct "
        "source_evidence→edges→knowledge_claim path over joining through paper_refs.",
    ]
    tables = [r[0] for r in con.execute(
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_schema='main' ORDER BY 1").fetchall()]
    for t in tables:
        if t in ("_meta",):
            continue
        cols = [r[0] for r in con.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name = ? ORDER BY ordinal_position", [t]).fetchall()]
        out.append(f'  {t}({", ".join(cols)})')
        for c in cols:
            if c in CATEGORICAL:
                vals = [str(v[0]) for v in con.execute(
                    f'SELECT DISTINCT "{c}" FROM "{t}" WHERE "{c}" IS NOT NULL '
                    f'LIMIT 20').fetchall()]
                if vals and len(vals) <= 20:
                    out.append(f'      {c} ∈ {{{", ".join(vals)}}}')
    # the actual relation vocabulary — the key to correct joins
    out.append("\nRELATIONS — to traverse, JOIN a node's id to edges.src_id and "
               "edges.dst_id to the other node's id. Available (src_type, rel, dst_type):")
    for s, rel, d in con.execute(
            "SELECT DISTINCT src_type, rel, dst_type FROM edges ORDER BY 1,2,3").fetchall():
        out.append(f"  ({s}) -[{rel}]-> ({d})")
    out.append("\nNote: each node's relation column also stores a JSON array of "
               "target ids; prefer the edges table for joins.")

    # supersession / authoritative-chain layer (supersession.py) — only if built
    has_sup = con.execute(
        "SELECT count(*) FROM information_schema.tables "
        "WHERE table_name='kb_record'").fetchone()[0]
    if has_sup:
        out.append(
            "\nAUTHORITATIVE-CHAIN LAYER (conclusions that were overturned over "
            "time). A 'topic' may have several records; older ones were SUPERSEDED "
            "and MUST NOT be reported as the current answer.\n"
            "  kb_record(topic, record_id, claim, gate, status, is_primary, "
            "authoritative_as_of, conclusion, aliases, source_path, source_kind)\n"
            "      status ∈ {authoritative, superseded, open}\n"
            "  supersession(topic, newer_id, older_id, source_record, source_kind) "
            "— newer_id supersedes older_id\n"
            "  authoritative_record(topic, record_id, claim, gate, status, "
            "authoritative_as_of, conclusion, source_path) — VIEW: the SINGLE "
            "current authoritative record per topic.\n"
            "  RULE: for any 'current / latest / authoritative status' question, "
            "SELECT FROM authoritative_record (or filter kb_record to "
            "status='authoritative' and MAX(authoritative_as_of)); treat "
            "status='superseded' rows as historical only.")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# authoritative-chain context — guarantees the newest authoritative conclusion
# is in front of `gen`, regardless of what SQL `syn` happened to write.
# --------------------------------------------------------------------------- #
def _load_kb_records(con) -> list[dict]:
    have = con.execute("SELECT count(*) FROM information_schema.tables "
                       "WHERE table_name='kb_record'").fetchone()[0]
    if not have:
        return []
    rows = con.execute(
        "SELECT topic, record_id, claim, status, is_primary, "
        "authoritative_as_of, conclusion, aliases FROM kb_record").fetchall()
    recs = []
    for topic, rid, claim, status, isprim, asof, concl, aliases in rows:
        try:
            al = json.loads(aliases) if aliases else []
        except Exception:
            al = []
        recs.append({"topic": topic, "record_id": rid, "claim": claim,
                     "status": status, "is_primary": (isprim == "true"),
                     "authoritative_as_of": asof, "conclusion": concl,
                     "aliases": al, "supersedes": [], "superseded_by": []})
    return recs


def authoritative_context(con, question: str) -> str:
    """If the question targets a topic with a supersession chain, return a block
    naming the current authoritative record + the superseded ones, so `gen`
    cannot answer with an overturned conclusion. Empty string when N/A."""
    if match_topic is None or resolve_authoritative is None:
        return ""
    recs = _load_kb_records(con)
    if not recs:
        return ""
    topic = match_topic(recs, question)
    if not topic:
        return ""
    auth = resolve_authoritative(recs, topic)
    if not auth or auth["status"] != "authoritative":
        return ""
    superseded = [r for r in recs
                  if r["topic"] == topic and r["status"] == "superseded"]
    lines = [
        f"AUTHORITATIVE RECORD for topic '{topic}' "
        f"(as_of {auth['authoritative_as_of']}): {auth['record_id']}",
        f"  Current conclusion: {auth['conclusion']}",
    ]
    if superseded:
        lines.append("  SUPERSEDED (historical — do NOT report as current): "
                     + ", ".join(sorted(r["record_id"] for r in superseded)))
    lines.append("  When answering about current/authoritative status, use the "
                 "AUTHORITATIVE RECORD above; mention superseded findings only as "
                 "history and say they were superseded.")
    return "\n".join(lines)


SYN_SYS = (
    "You are a DuckDB SQL expert for the ffn_cellsim knowledge base, a "
    "Table-Augmented-Generation (TAG) backend. Given the schema, write ONE "
    "DuckDB SELECT query that answers the user's question. Rules: SELECT/WITH "
    "only (never write). All node columns are TEXT — CAST(... AS DOUBLE) for "
    "numeric comparisons. Traverse relations via the edges table. Prefer EXACT id "
    "equality (vg_id = 'VG-...', mc_id = 'MC-...', run_id = 'RUN-...') over broad "
    "LIKE '%...%' matching — a fuzzy LIKE on an id often over-matches a whole "
    "family (e.g. '%h7%' catches both Gate-A and Gate-B). Use LIKE only when the "
    "user explicitly asks for a family/prefix. Return ONLY a single ```sql fenced "
    "code block, no prose."
)

FEWSHOT = textwrap.dedent("""\
    Q: How many SourceEvidence papers have no DOI?
    ```sql
    SELECT count(*) AS n_missing_doi FROM source_evidence WHERE doi IS NULL;
    ```

    Q: List ValidationGates with status 'failing' and the ModelContract each tests.
    ```sql
    SELECT vg.title AS gate, mc.title AS contract
    FROM validation_gate vg
    JOIN edges e ON e.src_id = vg.id AND e.dst_type = 'model_contract'
    JOIN model_contract mc ON mc.id = e.dst_id
    WHERE lower(vg.status) = 'failing';
    ```

    Q: What is the current authoritative status / conclusion for KU-3.5 cortical tension?
    ```sql
    SELECT record_id, status, authoritative_as_of, conclusion
    FROM authoritative_record
    WHERE lower(topic) LIKE '%cortical%' OR lower(topic) LIKE '%3.5%';
    ```

    Q: Which RunResults are linked to the H.7 Gate-A gate, with commit and outcome?
    ```sql
    SELECT r.run_id, r.commit, r.outcome
    FROM run_result r
    JOIN edges e ON e.src_id = r.id AND e.dst_type = 'validation_gate'
    JOIN validation_gate g ON g.id = e.dst_id
    WHERE g.vg_id = 'VG-H7-gate-a';
    ```

    Q: Which headline result-claims are NOT verified (retracted / needs-regen / GPU-unreproduced)?
    ```sql
    SELECT claim_id, verdict, metric, note FROM run_audit
    WHERE verdict <> 'VERIFIED' ORDER BY verdict;
    ```

    Q: Which simulation constants are NOT fully sourced (value drift / unsourced / unverified citation)?
    ```sql
    SELECT claim_id, location, ku, citation_key, verdict, note FROM param_audit
    WHERE verdict <> 'VERIFIED' ORDER BY verdict;
    ```

    Q: Which code modules implement the focal-adhesion motor-clutch contract?
    ```sql
    SELECT cm.path, cm.status
    FROM code_mapping cm
    JOIN edges e ON e.src_id = cm.id AND e.dst_type = 'model_contract'
    JOIN model_contract mc ON mc.id = e.dst_id
    WHERE mc.mc_id = 'MC-U2-fa-motor-clutch';
    ```
    """)


def syn(con, question: str, model: str | None, error: str = "", prev: str = "") -> str:
    schema = schema_text(con)
    repair = ""
    if error:
        repair = (f"\nYour previous SQL failed:\n{prev}\nDuckDB error:\n{error}\n"
                  "Fix it and return only the corrected ```sql block.")
    prompt = f"{schema}\n\nEXAMPLES:\n{FEWSHOT}\nQ: {question}{repair}\n```sql\n"
    raw = llm(prompt, system=SYN_SYS, model=model, max_tokens=800)
    m = re.search(r"```sql\s*(.*?)```", raw, re.S | re.I)
    sql = (m.group(1) if m else raw).strip().rstrip(";").strip()
    return sql


# --------------------------------------------------------------------------- #
# exec — SELECT-only guard + run
# --------------------------------------------------------------------------- #
def exec_sql(sql: str):
    low = sql.lstrip().lower()
    if not (low.startswith("select") or low.startswith("with")):
        raise ValueError("refused: not a SELECT/WITH query")
    if ";" in sql.rstrip(";"):
        raise ValueError("refused: multiple statements")
    if WRITE_RE.search(sql):
        raise ValueError("refused: write keyword in query")
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        cur = con.execute(sql)
        cols = [d[0] for d in cur.description]
        rows = cur.fetchall()
        return cols, rows
    finally:
        con.close()


# BM25 relevance floor for the content layer. Measured 2026-06-08 across several
# probes: purely RELATIONAL questions (answerable from the SQL rows alone) top out
# at ~4.3-5.5 on paper_chunks via incidental word overlap, while genuine CONTENT
# questions score ~6.0 (natural phrasing) up to ~10.5 (keyword-dense). A floor of
# 6.0 drops the incidental-overlap noise the `gen` step otherwise had to explain
# away (e.g. unrelated RHOA/Arslan excerpts surfacing on a "which runs link to
# Gate-A" query) while keeping real content hits. Absolute BM25 is query-length
# dependent, so this is a pragmatic separator, not a hard guarantee.
BM25_FLOOR = 6.0


def content_search(question: str, k: int = 4, floor: float = BM25_FLOOR):
    """BM25 over paper_chunks (the PDF content layer) — the TAG semantic step.

    Only excerpts scoring >= `floor` are returned, so relational queries (whose
    top BM25 hits are incidental-overlap noise) get NO excerpts instead of
    polluting the answer. Returns [] if the content layer isn't ingested."""
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        have = con.execute("SELECT count(*) FROM information_schema.tables "
                           "WHERE table_name='paper_chunks'").fetchone()[0]
        if not have:
            return []
        return con.execute(
            "SELECT citation_key, page, substr(text,1,400) AS snippet, "
            "fts_main_paper_chunks.match_bm25(chunk_id, ?) AS score "
            "FROM paper_chunks WHERE score IS NOT NULL AND score >= ? "
            "ORDER BY score DESC LIMIT ?", [question, floor, k]).fetchall()
    except Exception:
        return []
    finally:
        con.close()


def excerpts_to_text(ex) -> str:
    if not ex:
        return "(none)"
    return "\n".join(f"- [{ck} p{pg}] {sn.strip()[:300]}" for ck, pg, sn, _ in ex)


def rows_to_text(cols, rows, cap: int = 60) -> str:
    head = " | ".join(cols)
    body = "\n".join(" | ".join("" if v is None else str(v) for v in r)
                     for r in rows[:cap])
    extra = f"\n... ({len(rows) - cap} more rows)" if len(rows) > cap else ""
    return f"{head}\n{body}{extra}"


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


def gen(question: str, sql: str, cols, rows, excerpts, model: str | None,
        auth_context: str = "") -> str:
    table = rows_to_text(cols, rows)
    auth_block = (f"AUTHORITATIVE-CHAIN CONTEXT (overrides older statements):\n"
                  f"{auth_context}\n\n") if auth_context else ""
    prompt = (f"Question: {question}\n\n{auth_block}SQL used:\n{sql}\n\n"
              f"Result rows ({len(rows)}):\n{table}\n\n"
              f"PDF evidence excerpts (use only if relevant):\n"
              f"{excerpts_to_text(excerpts)}\n\nAnswer:")
    return llm(prompt, system=GEN_SYS, model=model, max_tokens=1200)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="TAG query over the ffn_cellsim KB")
    ap.add_argument("question")
    ap.add_argument("--model", default=None, help="override LLM model id")
    ap.add_argument("--sql-only", action="store_true", help="print SQL, skip exec/gen")
    ap.add_argument("--no-gen", action="store_true", help="exec but skip NL answer")
    args = ap.parse_args()

    if not DB_PATH.exists():
        sys.exit(f"{DB_PATH} not found — run notion_to_duckdb.py first.")
    con = duckdb.connect(str(DB_PATH), read_only=True)

    sql = syn(con, args.question, args.model)
    print(f"\n\033[36m── syn (SQL) ──\033[0m\n{sql}\n")
    if args.sql_only:
        return

    cols = rows = None
    for attempt in range(2):                                 # up to 2 self-repairs
        try:
            cols, rows = exec_sql(sql)
            break
        except Exception as e:
            print(f"\033[33m── exec error, repairing ({attempt + 1}/2) ──\033[0m\n{e}\n")
            sql = syn(con, args.question, args.model, error=str(e), prev=sql)
            print(f"\033[36m── syn (repaired) ──\033[0m\n{sql}\n")
    if rows is None:                                         # final attempt, let it raise
        cols, rows = exec_sql(sql)

    print(f"\033[36m── exec ({len(rows)} rows) ──\033[0m\n{rows_to_text(cols, rows, cap=20)}\n")
    if args.no_gen:
        return
    auth_context = authoritative_context(con, args.question)
    if auth_context:
        print(f"\033[35m── authoritative chain (latest wins) ──\033[0m\n"
              f"{auth_context}\n")
    excerpts = content_search(args.question)
    if excerpts:
        print(f"\033[36m── content (BM25 over paper_chunks) ──\033[0m\n"
              f"{excerpts_to_text(excerpts)}\n")
    print(f"\033[36m── gen (answer) ──\033[0m\n"
          f"{gen(args.question, sql, cols, rows, excerpts, args.model, auth_context)}")


if __name__ == "__main__":
    main()
