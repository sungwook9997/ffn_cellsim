#!/usr/bin/env python3
"""Measured common mechanics for manifest-backed KB integrity verifiers.

Only ``verify_runs`` and ``verify_params`` share this skeleton. The source,
force, and gate-contract verifiers have distinct inputs and verdict workflows
and intentionally do not depend on this module.
"""
from __future__ import annotations

import re
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import yaml

try:
    import duckdb
except Exception:  # duckdb is optional: both CI gates must run without it
    duckdb = None

try:
    import kb_coverage
except Exception as exc:  # pragma: no cover - exercised by negative controls
    kb_coverage = None
    COVERAGE_IMPORT_ERROR: Exception | None = exc
else:
    COVERAGE_IMPORT_ERROR = None

_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def load_claims(manifest: Path) -> list[dict[str, Any]]:
    """Load the conventional top-level ``claims`` list from a YAML manifest."""
    if not manifest.exists():
        return []
    return (yaml.safe_load(manifest.read_text()) or {}).get("claims", [])


def verdict_counts(
    rows: Sequence[Sequence[Any]], verdict_index: int
) -> Counter[str]:
    """Count verifier rows by their verdict column."""
    return Counter(str(row[verdict_index]) for row in rows)


def audit_manifest(
    manifest: Path,
    *,
    audit_claim: Callable[[dict[str, Any]], tuple[str, str]],
    declared_aliases: Mapping[str, str],
    suspicion: Mapping[str, int],
    row_builder: Callable[[dict[str, Any], str, str, str, str], Sequence[Any]],
) -> list[tuple[Any, ...]]:
    """Apply shared declared-vs-observed drift semantics to manifest claims."""
    rows = []
    for claim in load_claims(manifest):
        verdict, note = audit_claim(claim)
        declared = declared_aliases.get(
            (claim.get("declared") or "").strip().lower(), ""
        )
        drift = (
            "DRIFT"
            if declared
            and suspicion.get(verdict, 9) < suspicion.get(declared, 9)
            else ""
        )
        rows.append(tuple(row_builder(claim, verdict, declared, drift, note)))
    return rows


def write_audit_table(
    db_path: Path,
    table_name: str,
    columns: Sequence[str],
    rows: Sequence[Sequence[Any]],
) -> bool:
    """Replace one derived audit table when DuckDB is available."""
    identifiers = [table_name, *columns]
    if not all(_IDENTIFIER.fullmatch(name) for name in identifiers):
        raise ValueError(f"unsafe DuckDB identifier in {identifiers!r}")
    if duckdb is None:
        return False
    con = duckdb.connect(str(db_path))
    try:
        con.execute(f'DROP TABLE IF EXISTS "{table_name}"')
        ddl = ", ".join(f'"{column}" TEXT' for column in columns)
        con.execute(f'CREATE TABLE "{table_name}" ({ddl})')
        placeholders = ",".join("?" for _ in columns)
        con.executemany(
            f'INSERT INTO "{table_name}" VALUES ({placeholders})', rows
        )
    finally:
        con.close()
    return True


def run_coverage_gate(
    coverage: ModuleType | None,
    import_error: Exception | None,
    *,
    kind: str,
    label: str,
    manifest_name: str,
) -> int:
    """Run the shared omission-side ratchet, failing if it is unavailable."""
    if coverage is None:
        print(
            f"[{label}] GATE FAILED — coverage ratchet unavailable: "
            f"{import_error!r}"
        )
        print(
            "  kb_coverage.py is PART of this gate (it is the only check on what "
            f"{manifest_name} OMITS). A missing coverage check is a FAILURE, "
            "not a skip — restore aleph/outputs/tag_kb/kb_coverage.py."
        )
        return 1
    return coverage.gate(kind, label)


def check_coverage(
    coverage: ModuleType | None,
    import_error: Exception | None,
    *,
    kind: str,
    label: str,
) -> None:
    """Print the shared non-blocking omission-side summary."""
    if coverage is None:
        print(f"           ⚠️  coverage ratchet unavailable: {import_error!r}")
    else:
        coverage.check(kind, label)
