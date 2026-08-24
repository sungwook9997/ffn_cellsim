#!/usr/bin/env python3
"""CLI orchestration for the modular ffn_cellsim TAG query pipeline.

Pipeline modules:
    ``tag_syn``      NL question -> DuckDB SQL, with schema-hashed prompt cache
    ``tag_exec``     SELECT-only guard -> read-only DuckDB rows
    ``tag_gen``      rows + relevant excerpts -> grounded natural-language answer
    ``tag_backend``  Anthropic SDK cache or explicit uncached ``claude -p`` fallback

The command-line interface remains the stable entry point:

    python tag_query.py "Which Parameters feed a failing ValidationGate?"
    python tag_query.py --sql-only "How many SourceEvidence papers have no DOI?"
    python tag_query.py --model claude-opus-5 "..."
"""
from __future__ import annotations

import argparse
import logging
import sys
from typing import Any

import duckdb
from tag_backend import (
    DEFAULT_CACHE_TTL,
    DEFAULT_MODEL,
    SUPPORTED_CACHE_TTLS,
    llm,
)
from tag_exec import DB_PATH, WRITE_RE, exec_sql
from tag_gen import (
    BM25_FLOOR,
    GEN_SYS,
    authoritative_context,
    content_search,
    excerpts_to_text,
    rows_to_text,
)
from tag_gen import (
    gen as _gen,
)
from tag_syn import (
    CATEGORICAL,
    FEWSHOT,
    SYN_SYS,
    schema_hash,
    schema_text,
    static_cache_prefix,
)
from tag_syn import (
    syn as _syn,
)

__all__ = [
    "BM25_FLOOR",
    "CATEGORICAL",
    "DB_PATH",
    "DEFAULT_CACHE_TTL",
    "DEFAULT_MODEL",
    "FEWSHOT",
    "GEN_SYS",
    "SYN_SYS",
    "WRITE_RE",
    "authoritative_context",
    "content_search",
    "exec_sql",
    "excerpts_to_text",
    "gen",
    "llm",
    "rows_to_text",
    "schema_hash",
    "schema_text",
    "static_cache_prefix",
    "syn",
]


def syn(
    con: Any,
    question: str,
    model: str | None,
    error: str = "",
    prev: str = "",
    *,
    cache_ttl: str = DEFAULT_CACHE_TTL,
) -> str:
    """Compatibility wrapper preserving the historical ``tag_query.syn`` API."""
    return _syn(
        con,
        question,
        model,
        error,
        prev,
        cache_ttl=cache_ttl,
        llm_fn=llm,
    )


def gen(
    question: str,
    sql: str,
    cols: list[str],
    rows: list[tuple[Any, ...]],
    excerpts: list[tuple[Any, ...]],
    model: str | None,
    auth_context: str = "",
) -> str:
    """Compatibility wrapper preserving the historical ``tag_query.gen`` API."""
    return _gen(
        question,
        sql,
        cols,
        rows,
        excerpts,
        model,
        auth_context,
        llm_fn=llm,
    )


def main() -> None:
    """Run the stable syn -> exec -> gen command-line pipeline."""
    parser = argparse.ArgumentParser(description="TAG query over the ffn_cellsim KB")
    parser.add_argument("question")
    parser.add_argument(
        "--model",
        default=None,
        help=f"override LLM model id (default: {DEFAULT_MODEL})",
    )
    parser.add_argument(
        "--cache-ttl",
        choices=sorted(SUPPORTED_CACHE_TTLS),
        default=DEFAULT_CACHE_TTL,
        help=(
            "Anthropic prompt-cache TTL; ignored with an explicit warning on "
            f"the claude CLI fallback (default: {DEFAULT_CACHE_TTL})"
        ),
    )
    parser.add_argument(
        "--sql-only", action="store_true", help="print SQL, skip exec/gen"
    )
    parser.add_argument(
        "--no-gen", action="store_true", help="exec but skip NL answer"
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[tag-kb] %(message)s")
    if not DB_PATH.exists():
        sys.exit(f"{DB_PATH} not found — run notion_to_duckdb.py first.")
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        sql = syn(con, args.question, args.model, cache_ttl=args.cache_ttl)
        print(f"\n\033[36m── syn (SQL) ──\033[0m\n{sql}\n")
        if args.sql_only:
            return

        cols = rows = None
        for attempt in range(2):
            try:
                cols, rows = exec_sql(sql)
                break
            except Exception as exc:
                print(
                    "\033[33m── exec error, repairing "
                    f"({attempt + 1}/2) ──\033[0m\n{exc}\n"
                )
                sql = syn(
                    con,
                    args.question,
                    args.model,
                    error=str(exc),
                    prev=sql,
                    cache_ttl=args.cache_ttl,
                )
                print(f"\033[36m── syn (repaired) ──\033[0m\n{sql}\n")
        if rows is None:
            cols, rows = exec_sql(sql)

        print(
            f"\033[36m── exec ({len(rows)} rows) ──\033[0m\n"
            f"{rows_to_text(cols, rows, cap=20)}\n"
        )
        if args.no_gen:
            return
        auth_context = authoritative_context(con, args.question)
        if auth_context:
            print(
                "\033[35m── authoritative chain (latest wins) ──\033[0m\n"
                f"{auth_context}\n"
            )
        excerpts = content_search(args.question)
        if excerpts:
            print(
                "\033[36m── content (BM25 over paper_chunks) ──\033[0m\n"
                f"{excerpts_to_text(excerpts)}\n"
            )
        print(
            "\033[36m── gen (answer) ──\033[0m\n"
            f"{gen(args.question, sql, cols, rows, excerpts, args.model, auth_context)}"
        )
    finally:
        con.close()


if __name__ == "__main__":
    main()
