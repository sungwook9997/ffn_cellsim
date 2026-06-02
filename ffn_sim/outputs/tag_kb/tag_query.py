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
import os
import pathlib
import re
import subprocess
import sys
import textwrap

import duckdb

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
        max_tokens: int = 1500) -> str:
    """One-shot completion. anthropic SDK if key present, else `claude -p`."""
    key = os.environ.get("ANTHROPIC_API_KEY", "").strip()
    if key:
        import anthropic
        client = anthropic.Anthropic(api_key=key)
        msg = client.messages.create(
            model=model or "claude-opus-4-8",
            max_tokens=max_tokens,
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
    out: list[str] = ["TABLES (all node columns are TEXT; CAST when comparing numbers):"]
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
    return "\n".join(out)


SYN_SYS = (
    "You are a DuckDB SQL expert for the ffn_cellsim knowledge base, a "
    "Table-Augmented-Generation (TAG) backend. Given the schema, write ONE "
    "DuckDB SELECT query that answers the user's question. Rules: SELECT/WITH "
    "only (never write). All node columns are TEXT — CAST(... AS DOUBLE) for "
    "numeric comparisons. Traverse relations via the edges table. Return ONLY a "
    "single ```sql fenced code block, no prose."
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


def content_search(question: str, k: int = 4):
    """BM25 over paper_chunks (the PDF content layer) — the TAG semantic step.
    Returns [] if the content layer hasn't been ingested yet."""
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        have = con.execute("SELECT count(*) FROM information_schema.tables "
                           "WHERE table_name='paper_chunks'").fetchone()[0]
        if not have:
            return []
        return con.execute(
            "SELECT citation_key, page, substr(text,1,400) AS snippet, "
            "fts_main_paper_chunks.match_bm25(chunk_id, ?) AS score "
            "FROM paper_chunks WHERE score IS NOT NULL "
            "ORDER BY score DESC LIMIT ?", [question, k]).fetchall()
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
    "say no matching records were found. Do not invent data."
)


def gen(question: str, sql: str, cols, rows, excerpts, model: str | None) -> str:
    table = rows_to_text(cols, rows)
    prompt = (f"Question: {question}\n\nSQL used:\n{sql}\n\n"
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

    try:
        cols, rows = exec_sql(sql)
    except Exception as e:                                   # one self-repair pass
        print(f"\033[33m── exec error, repairing ──\033[0m\n{e}\n")
        sql = syn(con, args.question, args.model, error=str(e), prev=sql)
        print(f"\033[36m── syn (repaired) ──\033[0m\n{sql}\n")
        cols, rows = exec_sql(sql)

    print(f"\033[36m── exec ({len(rows)} rows) ──\033[0m\n{rows_to_text(cols, rows, cap=20)}\n")
    if args.no_gen:
        return
    excerpts = content_search(args.question)
    if excerpts:
        print(f"\033[36m── content (BM25 over paper_chunks) ──\033[0m\n"
              f"{excerpts_to_text(excerpts)}\n")
    print(f"\033[36m── gen (answer) ──\033[0m\n"
          f"{gen(args.question, sql, cols, rows, excerpts, args.model)}")


if __name__ == "__main__":
    main()
