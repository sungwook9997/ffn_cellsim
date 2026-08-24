#!/usr/bin/env python3
"""Audit the read-only boundary between external assays and Aleph observables."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


SCHEMA = "aleph.outer_library.observation_operator_registry.v1"
REPORT_SCHEMA = "aleph.outer_library.observation_operator_audit.v1"
ALLOWED_STATUSES = {
    "unavailable", "candidate_native", "native_validated", "surrogate_validated"
}
PROMOTABLE_STATUSES = {"native_validated", "surrogate_validated"}
VALIDATION_FIELDS = {
    "operator_id", "external_unit", "aleph_unit", "unit_conversion_verified",
    "geometry_matched", "validation_report_sha256",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(registry: dict[str, Any], repo_root: Path) -> dict[str, Any]:
    entries = registry.get("entries", [])
    duplicate_names = sorted({
        name for name in (entry.get("external_observable") for entry in entries)
        if sum(item.get("external_observable") == name for item in entries) > 1
    })
    audited = []
    for entry in entries:
        status = entry.get("status")
        source_path = entry.get("source_path")
        source_symbol = entry.get("source_symbol")
        source_exists = None
        source_symbol_present = None
        source_sha256 = None
        if source_path:
            source = repo_root / source_path
            source_exists = source.is_file()
            if source_exists:
                source_text = source.read_text(encoding="utf-8")
                source_symbol_present = bool(source_symbol and source_symbol in source_text)
                source_sha256 = _sha256(source)
        required_validation_fields_present = all(
            entry.get(field) is not None for field in VALIDATION_FIELDS
        )
        promotable = bool(
            status in PROMOTABLE_STATUSES
            and required_validation_fields_present
            and entry.get("unit_conversion_verified") is True
            and entry.get("geometry_matched") is True
        )
        audited.append({
            "external_observable": entry.get("external_observable"),
            "status": status,
            "status_allowed": status in ALLOWED_STATUSES,
            "operator_id": entry.get("operator_id"),
            "source_path": source_path,
            "source_exists": source_exists,
            "source_symbol_present": source_symbol_present,
            "source_sha256": source_sha256,
            "required_validation_fields_present": required_validation_fields_present,
            "promotable_to_numeric_sweep": promotable,
            "forbidden_shortcut_recorded": bool(entry.get("forbidden_shortcut")),
        })
    checks = {
        "schema_supported": registry.get("schema") == SCHEMA,
        "authority_is_conservative_inventory": (
            registry.get("authority") == "conservative_pre_sweep_inventory"
        ),
        "entries_present": bool(entries),
        "external_observables_unique": not duplicate_names,
        "all_statuses_allowed": all(item["status_allowed"] for item in audited),
        "all_declared_source_evidence_resolves": all(
            item["source_exists"] is not False
            and item["source_symbol_present"] is not False
            for item in audited
        ),
        "all_shortcuts_explicitly_forbidden": all(
            item["forbidden_shortcut_recorded"] for item in audited
        ),
        "unavailable_entries_have_no_operator": all(
            item["status"] != "unavailable" or item["operator_id"] is None
            for item in audited
        ),
        "validated_statuses_have_complete_validation": all(
            item["status"] not in PROMOTABLE_STATUSES
            or item["promotable_to_numeric_sweep"]
            for item in audited
        ),
        "no_unvalidated_entry_is_promotable": all(
            not item["promotable_to_numeric_sweep"]
            or item["status"] in PROMOTABLE_STATUSES
            for item in audited
        ),
    }
    passed = all(checks.values())
    return {
        "schema": REPORT_SCHEMA,
        "registry_schema": registry.get("schema"),
        "checks": checks,
        "passed": passed,
        "failed_checks": [name for name, value in checks.items() if not value],
        "duplicate_external_observables": duplicate_names,
        "entries": audited,
        "promotable_external_observables": sorted(
            item["external_observable"] for item in audited
            if item["promotable_to_numeric_sweep"]
        ),
        "unavailable_external_observables": sorted(
            item["external_observable"] for item in audited
            if item["status"] == "unavailable"
        ),
        "candidate_native_external_observables": sorted(
            item["external_observable"] for item in audited
            if item["status"] == "candidate_native"
        ),
        "aleph_source_files_mutated": False,
        "may_emit_numeric_parameter_range": False,
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(
        json.loads(args.registry.read_text(encoding="utf-8")),
        args.repo_root.resolve(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
