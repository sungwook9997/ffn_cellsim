#!/usr/bin/env python3
"""TAG ``exec`` stage: SELECT-only validation and read-only DuckDB execution."""
from __future__ import annotations

import pathlib
import re
from typing import Any

import duckdb

HERE = pathlib.Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"

WRITE_RE = re.compile(
    r"\b(insert|update|delete|drop|create|alter|attach|copy|"
    r"pragma|replace|truncate|install|load)\b",
    re.I,
)


def validate_select(sql: str) -> None:
    """Reject writes, multiple statements, and non-query SQL."""
    low = sql.lstrip().lower()
    if not (low.startswith("select") or low.startswith("with")):
        raise ValueError("refused: not a SELECT/WITH query")
    if ";" in sql.rstrip(";"):
        raise ValueError("refused: multiple statements")
    if WRITE_RE.search(sql):
        raise ValueError("refused: write keyword in query")


def exec_sql(
    sql: str, *, db_path: pathlib.Path = DB_PATH
) -> tuple[list[str], list[tuple[Any, ...]]]:
    """Execute one validated query against the read-only materialised mirror."""
    validate_select(sql)
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        cursor = con.execute(sql)
        cols = [description[0] for description in cursor.description]
        return cols, cursor.fetchall()
    finally:
        con.close()
