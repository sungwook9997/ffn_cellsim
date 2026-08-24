#!/usr/bin/env python3
"""TAG ``syn`` stage: live DuckDB schema plus examples to NL-to-SQL."""
from __future__ import annotations

import hashlib
import re
import textwrap
from collections.abc import Callable
from typing import Any

from tag_backend import DEFAULT_CACHE_TTL, llm

CATEGORICAL = {"source_type", "status", "confidence", "unit", "subtopic", "kind"}

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


def schema_text(con: Any) -> str:
    """Render the live table, categorical-value, and relation vocabulary."""
    out: list[str] = [
        "TABLES (all node columns are TEXT; CAST when comparing numbers):",
        "KEY IDS — every node's `id` is an OPAQUE Notion-page UUID (e.g. "
        "'372120da…'), used only for JOINs via the edges table. The HUMAN-READABLE "
        "identifier is a separate column: knowledge_claim.kb_id ('KB-3.19'), "
        "source_evidence.citation_key ('Bell1978_Science') and .uid ('SE119'). "
        "To filter by a 'KB-x.y' value use kb_id, NEVER id. Prefer the direct "
        "source_evidence→edges→knowledge_claim path over joining through paper_refs.",
    ]
    tables = [
        row[0]
        for row in con.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema='main' ORDER BY 1"
        ).fetchall()
    ]
    for table in tables:
        if table == "_meta":
            continue
        cols = [
            row[0]
            for row in con.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name = ? ORDER BY ordinal_position",
                [table],
            ).fetchall()
        ]
        out.append(f'  {table}({", ".join(cols)})')
        for column in cols:
            if column in CATEGORICAL:
                values = [
                    str(row[0])
                    for row in con.execute(
                        f'SELECT DISTINCT "{column}" FROM "{table}" '
                        f'WHERE "{column}" IS NOT NULL ORDER BY 1 LIMIT 20'
                    ).fetchall()
                ]
                if values and len(values) <= 20:
                    out.append(f'      {column} ∈ {{{", ".join(values)}}}')
    out.append(
        "\nRELATIONS — to traverse, JOIN a node's id to edges.src_id and "
        "edges.dst_id to the other node's id. Available (src_type, rel, dst_type):"
    )
    for source, relation, destination in con.execute(
        "SELECT DISTINCT src_type, rel, dst_type FROM edges ORDER BY 1,2,3"
    ).fetchall():
        out.append(f"  ({source}) -[{relation}]-> ({destination})")
    out.append(
        "\nNote: each node's relation column also stores a JSON array of "
        "target ids; prefer the edges table for joins."
    )

    has_supersession = con.execute(
        "SELECT count(*) FROM information_schema.tables "
        "WHERE table_name='kb_record'"
    ).fetchone()[0]
    if has_supersession:
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
            "status='superseded' rows as historical only."
        )
    return "\n".join(out)


def schema_hash(schema: str) -> str:
    """Return the cache-key component that changes with the rendered schema."""
    return hashlib.sha256(schema.encode("utf-8")).hexdigest()


def static_cache_prefix(schema: str) -> tuple[str, str]:
    """Build the exact cached block and its visible schema-derived key."""
    digest = schema_hash(schema)
    prefix = f"SCHEMA_SHA256: {digest}\n\n{schema}\n\nEXAMPLES:\n{FEWSHOT}"
    return prefix, digest


def syn(
    con: Any,
    question: str,
    model: str | None,
    error: str = "",
    prev: str = "",
    *,
    cache_ttl: str = DEFAULT_CACHE_TTL,
    llm_fn: Callable[..., str] = llm,
) -> str:
    """Synthesize one guarded DuckDB query from a natural-language question."""
    schema = schema_text(con)
    cache_prefix, cache_key = static_cache_prefix(schema)
    repair = ""
    if error:
        repair = (
            f"\nYour previous SQL failed:\n{prev}\nDuckDB error:\n{error}\n"
            "Fix it and return only the corrected ```sql block."
        )
    dynamic_prompt = f"Q: {question}{repair}\n```sql\n"
    raw = llm_fn(
        dynamic_prompt,
        system=SYN_SYS,
        model=model,
        max_tokens=800,
        cache_prefix=cache_prefix,
        cache_key=cache_key,
        cache_ttl=cache_ttl,
    )
    match = re.search(r"```sql\s*(.*?)```", raw, re.S | re.I)
    return (match.group(1) if match else raw).strip().rstrip(";").strip()
