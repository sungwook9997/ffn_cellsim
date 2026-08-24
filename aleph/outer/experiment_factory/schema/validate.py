#!/usr/bin/env python
"""Dependency-free semantic validation for proposed experiment-factory records."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterable


SHA256 = re.compile(r"^[0-9a-f]{64}$")
EXPERIMENT_ID = re.compile(r"^aleph:experiment:[0-9a-f]{16}$")
OBSERVATION_ID = re.compile(r"^aleph:observation:[0-9a-f]{16}$")
PMC_ID = re.compile(r"^PMC[0-9]+$")

EXPERIMENT_REQUIRED = {
    "schema", "authority_status", "review_state", "experiment_record_id", "source",
    "biological_systems", "conditions", "protocols", "observations", "comparisons",
    "provenance",
}
CANDIDATE_REQUIRED = {
    "schema", "authority_status", "review_state", "selection_rank", "source_family_id",
    "pmcid", "title", "journal", "acquisition_pass", "publication_year", "year_stratum", "primary_domain",
    "discovery_domain", "modality_tags", "cell_state_tags", "mechanics_tags",
    "priority_score", "observation_locator_count", "leakage_component_id",
    "eligible_component_size", "related_out_of_pool_count", "selection_reason",
    "source_record_sha256", "source_payload_sha256",
}


def _require_keys(value: dict[str, Any], keys: set[str], path: str, errors: list[str]) -> None:
    for key in sorted(keys - value.keys()):
        errors.append(f"{path}: missing required key {key!r}")


def _sha(value: Any, path: str, errors: list[str], *, nullable: bool = False) -> None:
    if nullable and value is None:
        return
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        errors.append(f"{path}: expected lowercase SHA-256")


def _evidence(value: Any, path: str, errors: list[str]) -> None:
    required = {
        "document_sha256", "section_path", "locator_type", "locator_id",
        "text_or_asset_sha256", "char_start", "char_end",
    }
    if not isinstance(value, dict):
        errors.append(f"{path}: expected object")
        return
    _require_keys(value, required, path, errors)
    _sha(value.get("document_sha256"), f"{path}.document_sha256", errors)
    _sha(value.get("text_or_asset_sha256"), f"{path}.text_or_asset_sha256", errors)
    if not isinstance(value.get("section_path"), list):
        errors.append(f"{path}.section_path: expected array")
    start, end = value.get("char_start"), value.get("char_end")
    if start is not None and (not isinstance(start, int) or start < 0):
        errors.append(f"{path}.char_start: expected non-negative integer or null")
    if end is not None and (not isinstance(end, int) or end < 0):
        errors.append(f"{path}.char_end: expected non-negative integer or null")
    if isinstance(start, int) and isinstance(end, int) and end < start:
        errors.append(f"{path}: char_end precedes char_start")


def _evidence_array(value: Any, path: str, errors: list[str], *, nonempty: bool = True) -> None:
    if not isinstance(value, list):
        errors.append(f"{path}: expected array")
        return
    if nonempty and not value:
        errors.append(f"{path}: at least one exact evidence locator is required")
    for index, locator in enumerate(value):
        _evidence(locator, f"{path}[{index}]", errors)


def validate_candidate(record: Any) -> list[str]:
    errors: list[str] = []
    if not isinstance(record, dict):
        return ["$: expected object"]
    _require_keys(record, CANDIDATE_REQUIRED, "$", errors)
    if record.get("schema") != "aleph.external_training.annotation_candidate.v1":
        errors.append("$.schema: unsupported candidate schema")
    if record.get("authority_status") != "proposed":
        errors.append("$.authority_status: candidates must remain proposed")
    if record.get("review_state") != "gold_candidate":
        errors.append("$.review_state: candidate set is not human-verified gold")
    if record.get("acquisition_pass") not in {"first", "second"}:
        errors.append("$.acquisition_pass: expected first or second")
    year = record.get("publication_year")
    if not isinstance(year, int) or not 2010 <= year <= 2026:
        errors.append("$.publication_year: expected integer in 2010..2026")
    if not isinstance(record.get("selection_rank"), int) or record.get("selection_rank", 0) < 1:
        errors.append("$.selection_rank: expected positive integer")
    pmcid = record.get("pmcid")
    if pmcid is not None and (not isinstance(pmcid, str) or PMC_ID.fullmatch(pmcid) is None):
        errors.append("$.pmcid: expected PMCID or null")
    component = record.get("leakage_component_id")
    if not isinstance(component, str) or re.fullmatch(r"lc:[0-9a-f]{16}", component) is None:
        errors.append("$.leakage_component_id: invalid component ID")
    _sha(record.get("source_record_sha256"), "$.source_record_sha256", errors)
    _sha(record.get("source_payload_sha256"), "$.source_payload_sha256", errors)
    return errors


def validate_experiment(record: Any) -> list[str]:
    errors: list[str] = []
    if not isinstance(record, dict):
        return ["$: expected object"]
    _require_keys(record, EXPERIMENT_REQUIRED, "$", errors)
    if record.get("schema") != "aleph.external_training.experiment_record.v1":
        errors.append("$.schema: unsupported experiment schema")
    if record.get("authority_status") != "proposed":
        errors.append("$.authority_status: external annotations must remain proposed")
    if not isinstance(record.get("experiment_record_id"), str) or EXPERIMENT_ID.fullmatch(record["experiment_record_id"]) is None:
        errors.append("$.experiment_record_id: invalid deterministic ID")

    source = record.get("source")
    if not isinstance(source, dict):
        errors.append("$.source: expected object")
    else:
        _require_keys(source, {"source_family_id", "pmcid", "acquisition_pass", "publication_year", "primary_domain"}, "$.source", errors)

    systems = record.get("biological_systems")
    if not isinstance(systems, list) or not systems:
        errors.append("$.biological_systems: at least one system is required")
        systems = []
    system_ids = {x.get("biological_system_id") for x in systems if isinstance(x, dict)}

    conditions = record.get("conditions")
    if not isinstance(conditions, list) or not conditions:
        errors.append("$.conditions: at least one condition is required")
        conditions = []
    condition_ids = {x.get("condition_id") for x in conditions if isinstance(x, dict)}
    if len(condition_ids) != len(conditions):
        errors.append("$.conditions: condition IDs must be present and unique")
    for index, condition in enumerate(conditions):
        if isinstance(condition, dict) and condition.get("biological_system_id") not in system_ids:
            errors.append(f"$.conditions[{index}].biological_system_id: unresolved reference")

    protocols = record.get("protocols")
    if not isinstance(protocols, list) or not protocols:
        errors.append("$.protocols: at least one protocol is required")
        protocols = []
    protocol_ids = {x.get("protocol_id") for x in protocols if isinstance(x, dict)}
    if len(protocol_ids) != len(protocols):
        errors.append("$.protocols: protocol IDs must be present and unique")
    for index, protocol in enumerate(protocols):
        if isinstance(protocol, dict):
            _evidence_array(protocol.get("evidence"), f"$.protocols[{index}].evidence", errors)

    observations = record.get("observations")
    if not isinstance(observations, list):
        errors.append("$.observations: expected array")
        observations = []
    observation_ids: set[str] = set()
    for index, observation in enumerate(observations):
        path = f"$.observations[{index}]"
        if not isinstance(observation, dict):
            errors.append(f"{path}: expected object")
            continue
        oid = observation.get("observation_id")
        if not isinstance(oid, str) or OBSERVATION_ID.fullmatch(oid) is None:
            errors.append(f"{path}.observation_id: invalid deterministic ID")
        elif oid in observation_ids:
            errors.append(f"{path}.observation_id: duplicate ID")
        observation_ids.add(oid)
        if observation.get("protocol_id") not in protocol_ids:
            errors.append(f"{path}.protocol_id: unresolved reference")
        for cid in observation.get("condition_ids", []):
            if cid not in condition_ids:
                errors.append(f"{path}.condition_ids: unresolved condition {cid!r}")
        state = observation.get("value_state")
        value = observation.get("value")
        if state == "missing" and value is not None:
            errors.append(f"{path}.value: missing state requires null value")
        if state in {"reported", "extracted", "derived"} and value is None:
            errors.append(f"{path}.value: populated state requires a value")
        unit = observation.get("unit")
        if not isinstance(unit, dict):
            errors.append(f"{path}.unit: expected explicit unit object")
        elif isinstance(value, (int, float)) and state in {"reported", "extracted", "derived"}:
            if unit.get("reported") is None and unit.get("normalization_status") not in {"dimensionless", "unknown"}:
                errors.append(f"{path}.unit: numeric value requires reported unit, dimensionless, or explicit unknown")
            if unit.get("normalization_status") == "converted" and (not unit.get("conversion_formula") or not unit.get("conversion_evidence")):
                errors.append(f"{path}.unit: conversion requires formula and evidence")
        extraction = observation.get("extraction")
        if not isinstance(extraction, dict):
            errors.append(f"{path}.extraction: expected object")
        elif extraction.get("status") == "ambiguous" and not (extraction.get("ambiguity_reason") and extraction.get("candidates")):
            errors.append(f"{path}.extraction: ambiguous status requires reason and candidates")
        _evidence_array(observation.get("evidence"), f"{path}.evidence", errors)

    comparisons = record.get("comparisons")
    if not isinstance(comparisons, list):
        errors.append("$.comparisons: expected array")
        comparisons = []
    for index, comparison in enumerate(comparisons):
        path = f"$.comparisons[{index}]"
        if not isinstance(comparison, dict):
            errors.append(f"{path}: expected object")
            continue
        for cid in comparison.get("left_condition_ids", []) + comparison.get("right_condition_ids", []):
            if cid not in condition_ids:
                errors.append(f"{path}: unresolved condition {cid!r}")
        for oid in comparison.get("observation_ids", []):
            if oid not in observation_ids:
                errors.append(f"{path}: unresolved observation {oid!r}")
        if comparison.get("status") == "ambiguous" and not comparison.get("ambiguity_reason"):
            errors.append(f"{path}: ambiguous comparison requires a reason")
        _evidence_array(comparison.get("evidence"), f"{path}.evidence", errors)

    provenance = record.get("provenance")
    if not isinstance(provenance, dict):
        errors.append("$.provenance: expected object")
    else:
        _sha(provenance.get("source_receipt_sha256"), "$.provenance.source_receipt_sha256", errors)
        _sha(provenance.get("candidate_record_sha256"), "$.provenance.candidate_record_sha256", errors)
        _sha(provenance.get("record_sha256"), "$.provenance.record_sha256", errors, nullable=True)
        review_state = record.get("review_state")
        annotators = provenance.get("annotator_ids")
        if review_state == "adjudicated" and (not isinstance(annotators, list) or len(set(annotators)) < 2 or provenance.get("adjudicator_id") is None):
            errors.append("$.provenance: adjudicated requires two annotators and an adjudicator")
    return errors


def iter_records(path: Path) -> Iterable[tuple[int, Any]]:
    if path.suffix == ".jsonl":
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if line.strip():
                yield line_number, json.loads(line)
    else:
        yield 1, json.loads(path.read_text(encoding="utf-8"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--kind", choices=("auto", "experiment", "candidate"), default="auto")
    args = parser.parse_args(argv)
    failures = 0
    records = 0
    for raw_path in args.paths:
        path = Path(raw_path)
        try:
            rows = iter_records(path)
            for line_number, record in rows:
                records += 1
                kind = args.kind
                if kind == "auto":
                    kind = "candidate" if isinstance(record, dict) and record.get("schema") == "aleph.external_training.annotation_candidate.v1" else "experiment"
                errors = validate_candidate(record) if kind == "candidate" else validate_experiment(record)
                for error in errors:
                    failures += 1
                    print(f"{path}:{line_number}: {error}", file=sys.stderr)
        except (OSError, json.JSONDecodeError) as exc:
            failures += 1
            print(f"{path}: {exc}", file=sys.stderr)
    if failures:
        print(json.dumps({"records": records, "errors": failures, "status": "invalid"}, sort_keys=True))
        return 1
    print(json.dumps({"records": records, "errors": 0, "status": "valid"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
