#!/usr/bin/env python3
"""Audit external cortical-tension evidence and compare eligible engine records.

The comparator is deliberately stricter than a plotting script.  It redacts all
engine magnitudes unless the run records establish full-native CUDA execution,
all compartments, physical accepted steps, time-step independence, stationarity,
a non-blocked quantitative claim, protocol identity, a validated observation
operator, and measured between-seed scatter.  It never turns overlap with an
external interval into an Aleph validation-gate verdict.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median, stdev
from typing import Any, Iterable

import duckdb


CATALOG_SCHEMA = "aleph.outer_library.cortical_tension_reference_catalog.v1"
REPORT_SCHEMA = "aleph.outer_library.cortical_tension_comparison.v1"
REFERENCE_ID = "hosseini2020_mcf7_suspended_interphase_afm"
OBSERVABLE = "effective_laplace_surface_tension"
ENGINE_OBSERVABLE = "gamma_total_pn_per_um"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _record(catalog: dict[str, Any], reference_id: str = REFERENCE_ID) -> dict[str, Any]:
    matches = [item for item in catalog.get("records", []) if item.get("reference_id") == reference_id]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one reference {reference_id!r}; found {len(matches)}")
    return matches[0]


def audit_catalog(catalog: dict[str, Any], kb_path: Path) -> dict[str, Any]:
    """Cross-check the structured reference against the read-only KB snapshot."""
    reference = _record(catalog)
    evidence = reference["evidence"]
    connection = duckdb.connect(str(kb_path), read_only=True)
    try:
        rows = connection.execute(
            """
            SELECT k.kb_id, k.status, k.confidence, k.value_range_si,
                   se.uid, se.citation_key, se.doi, sa.verdict, vg.status, vg.band_expected
            FROM knowledge_claim k
            JOIN edges e ON e.src_id = k.id AND e.dst_type = 'source_evidence'
            JOIN source_evidence se ON se.id = e.dst_id
            LEFT JOIN source_audit sa ON sa.citation_key = se.citation_key
            LEFT JOIN validation_gate vg ON vg.vg_id = ?
            WHERE k.kb_id = ? AND se.citation_key = ?
            """,
            [
                reference["gate_context"]["validation_gate_id"],
                evidence["knowledge_claim_id"],
                evidence["citation_key"],
            ],
        ).fetchall()
    finally:
        connection.close()
    row = rows[0] if len(rows) == 1 else None
    value_text = row[3] if row else ""
    checks = {
        "catalog_schema_supported": catalog.get("schema") == CATALOG_SCHEMA,
        "single_primary_reference": len(catalog.get("records", [])) == 1,
        "outer_has_no_aleph_authority": (
            catalog.get("aleph_authority") == "none"
            and catalog.get("may_emit_aleph_parameter") is False
            and catalog.get("may_mutate_physics") is False
        ),
        "kb_join_unique": len(rows) == 1,
        "knowledge_claim_verified": bool(row and row[1] == "verified"),
        "source_evidence_identity_matches": bool(
            row
            and row[4] == evidence["source_evidence_uid"]
            and row[5] == evidence["citation_key"]
            and str(row[6]).lower().rstrip("/").endswith(evidence["doi"].lower())
        ),
        "source_audit_ok": bool(row and row[7] == evidence["source_audit_required_verdict"] == "OK"),
        "numeric_summary_matches_claim": bool(
            row
            and all(token in value_text for token in ("0.27 mN/m", "0.18-0.40", "n=27"))
            and reference["summary"]["median_pn_per_um"] == 270.0
            and reference["summary"]["q1_pn_per_um"] == 180.0
            and reference["summary"]["q3_pn_per_um"] == 400.0
            and reference["summary"]["n_cells"] == 27
        ),
        "draft_gate_not_promoted": bool(
            row
            and row[8] == reference["gate_context"]["status"] == "draft"
            and reference["gate_context"]["may_score_pass_fail"] is False
        ),
        "protocol_is_total_effective_laplace_tension": (
            reference["measurement_operator"]["external_observable"] == OBSERVABLE
            and "total effective Laplace" in reference["measurement_operator"]["definition"]
        ),
        "coverage_gaps_are_explicit": len(catalog.get("coverage_gaps", [])) >= 3,
        "forbidden_shortcuts_recorded": len(catalog.get("forbidden_shortcuts", [])) >= 5,
    }
    return {
        "schema": "aleph.outer_library.cortical_tension_catalog_audit.v1",
        "passed": all(checks.values()),
        "checks": checks,
        "failed_checks": [name for name, passed in checks.items() if not passed],
        "kb_snapshot": {"path": str(kb_path), "sha256": _sha256(kb_path)},
        "resolved_evidence": None if row is None else {
            "knowledge_claim_id": row[0],
            "knowledge_claim_status": row[1],
            "knowledge_claim_confidence": row[2],
            "source_evidence_uid": row[4],
            "citation_key": row[5],
            "doi": row[6],
            "source_audit_verdict": row[7],
            "validation_gate_status": row[8],
            "validation_gate_band": row[9],
        },
        "aleph_authority": "none",
        "may_score_validation_gate": False,
    }


def infer_reference(reference: dict[str, Any]) -> dict[str, Any]:
    """Expose the audited group summary as a protocol-conditioned Outer head."""
    summary = reference["summary"]
    return {
        "schema": "aleph.outer_library.cortical_tension_reference_head.v1",
        "head_id": "mcf7_suspended_interphase_effective_cortical_tension",
        "kind": "protocol_conditioned_external_group_summary",
        "refused": False,
        "state_posterior": [{
            "latent_id": "effective_laplace_surface_tension",
            "estimate": summary["median_pn_per_um"],
            "interval": [summary["q1_pn_per_um"], summary["q3_pn_per_um"]],
            "interval_kind": summary["interval_kind"],
            "unit": summary["unit"],
            "conditioning": reference["context"],
            "authority": "external_observation_summary",
        }],
        "uncertainty": {
            "kind": "reported_empirical_IQR_plus_figure_digitization_precision",
            "n_cells": summary["n_cells"],
            "digitization_precision_pn_per_um": summary["digitization_precision_pn_per_um"],
            "calibrated_predictive_distribution": False,
        },
        "applicability": {
            "status": "exact_context_only",
            "context": reference["context"],
            "measurement_operator": reference["measurement_operator"],
        },
        "limitations": [
            "one study and one cell-state/protocol group cannot train a general cortical-tension predictor",
            "the empirical IQR is not a calibrated predictive interval for a new laboratory",
            "the draft validation gate is context only and is not scored here",
        ],
        "aleph_dependency": False,
        "may_emit_Aleph_parameter": False,
        "may_score_validation_gate": False,
    }


def _operator_entry(registry: dict[str, Any]) -> dict[str, Any] | None:
    for item in registry.get("entries", []):
        if item.get("external_observable") == OBSERVABLE:
            return item
    return None


def _engine_readiness(record: dict[str, Any], reference: dict[str, Any]) -> dict[str, Any]:
    context = record.get("comparison_context", {})
    census = record.get("census", {})
    stationarity = record.get("measurements", {}).get("stationarity", {}).get(ENGINE_OBSERVABLE, {})
    value = stationarity.get("mean")
    checks = {
        "run_record_schema_supported": record.get("schema") == "ac.engine.observe/run-record@2",
        "self_stamped_build_present": bool(record.get("build", {}).get("commit")),
        "cuda_runtime": str(record.get("device", "")).startswith("cuda"),
        "full_native_population": (
            census.get("fraction_of_native") == 1.0
            or census.get("is_native_population") is True
        ),
        "full_native_compartments": context.get("full_native_compartments") is True,
        "physical_acceptance_gate_passed": context.get("physical_acceptance_gate_passed") is True,
        "time_step_independence_passed": context.get("time_step_independence_passed") is True,
        "quantitative_claim_not_blocked": bool(
            record.get("quantitative_claim_status")
            and record.get("quantitative_claim_status") != "BLOCKED"
        ),
        "target_reference_matches": context.get("target_reference_id") == reference["reference_id"],
        "cell_line_matches": context.get("cell_line") == reference["context"]["cell_line"],
        "surface_state_matches": context.get("surface_state") == reference["context"]["surface_state"],
        "cell_cycle_state_matches": context.get("cell_cycle_state") == reference["context"]["cell_cycle_state"],
        "temperature_matches": context.get("temperature_c") == reference["context"]["temperature_c"],
        "engine_operator_is_total_method_of_planes": context.get("measurement_operator") == "gamma_total_method_of_planes",
        "stationary_total_tension_mean_present": bool(
            stationarity.get("verdict") == "STATIONARY"
            and isinstance(value, (int, float))
            and math.isfinite(float(value))
        ),
        "seed_present": isinstance(record.get("config", {}).get("seed"), int),
    }
    return {
        "record_id": record.get("run_label") or record.get("build", {}).get("commit") or "unknown",
        "seed": record.get("config", {}).get("seed"),
        "checks": checks,
        "ready_except_operator_and_seed_ensemble": all(checks.values()),
        "failed_checks": [name for name, passed in checks.items() if not passed],
        "_value": float(value) if checks["stationary_total_tension_mean_present"] else None,
    }


def compare_engine_records(
    records: Iterable[dict[str, Any]],
    catalog: dict[str, Any],
    catalog_audit: dict[str, Any],
    operator_registry: dict[str, Any],
) -> dict[str, Any]:
    """Return a descriptive comparison or a value-redacted readiness report."""
    reference = _record(catalog)
    readiness = [_engine_readiness(item, reference) for item in records]
    operator = _operator_entry(operator_registry)
    operator_ready = bool(
        operator
        and operator.get("status") == "native_validated"
        and operator.get("geometry_matched") is True
        and operator.get("unit_conversion_verified") is True
        and isinstance(operator.get("validation_report_sha256"), str)
        and len(operator["validation_report_sha256"]) == 64
    )
    seeds = [item["seed"] for item in readiness if item["seed"] is not None]
    seed_scatter_measured = len(readiness) >= 3 and len(set(seeds)) == len(readiness)
    comparison_permitted = bool(
        catalog_audit.get("passed")
        and readiness
        and all(item["ready_except_operator_and_seed_ensemble"] for item in readiness)
        and operator_ready
        and seed_scatter_measured
    )
    public_readiness = [
        {key: value for key, value in item.items() if key != "_value"}
        for item in readiness
    ]
    blockers = []
    if not catalog_audit.get("passed"):
        blockers.append("external_reference_catalog_audit_failed")
    blockers.extend(
        f"engine_record[{index}]:{name}"
        for index, item in enumerate(readiness)
        for name in item["failed_checks"]
    )
    if not operator_ready:
        blockers.append("observation_operator_not_native_validated_and_geometry_matched")
    if not seed_scatter_measured:
        blockers.append("between_seed_scatter_not_measured_with_three_distinct_seed_records")
    result: dict[str, Any] = {
        "schema": REPORT_SCHEMA,
        "reference_id": reference["reference_id"],
        "comparison_permitted": comparison_permitted,
        "comparison_kind": "descriptive_protocol_matched_observable_comparison",
        "validation_gate_verdict": None,
        "engine_record_readiness": public_readiness,
        "operator_readiness": {
            "registry_entry": operator,
            "native_validated_and_geometry_matched": operator_ready,
        },
        "seed_ensemble": {
            "record_count": len(readiness),
            "distinct_seed_count": len(set(seeds)),
            "between_seed_scatter_measured": seed_scatter_measured,
        },
        "blockers": sorted(set(blockers)),
        "engine_summary": None,
        "external_reference": None,
        "descriptive_position": None,
        "redaction": None,
        "aleph_authority": "none",
        "may_emit_Aleph_parameter": False,
        "may_mutate_physics": False,
        "may_score_validation_gate": False,
    }
    if not comparison_permitted:
        result["redaction"] = {
            "engine_magnitudes_and_ratios_redacted": True,
            "reason": "At least one evidence, protocol, runtime, acceptance, convergence, operator, or seed-scatter precondition is unresolved.",
        }
        return result
    values = [item["_value"] for item in readiness]
    assert all(value is not None for value in values)
    numeric_values = [float(value) for value in values if value is not None]
    engine_median = median(numeric_values)
    engine_sd = stdev(numeric_values)
    summary = reference["summary"]
    result["engine_summary"] = {
        "median_pn_per_um": engine_median,
        "between_seed_sd_pn_per_um": engine_sd,
        "minimum_pn_per_um": min(numeric_values),
        "maximum_pn_per_um": max(numeric_values),
        "seed_count": len(numeric_values),
        "within_run_sem_used": False,
    }
    result["external_reference"] = {
        "median_pn_per_um": summary["median_pn_per_um"],
        "q1_pn_per_um": summary["q1_pn_per_um"],
        "q3_pn_per_um": summary["q3_pn_per_um"],
        "interval_kind": summary["interval_kind"],
        "n_cells": summary["n_cells"],
    }
    result["descriptive_position"] = {
        "engine_minus_reference_median_pn_per_um": engine_median - summary["median_pn_per_um"],
        "engine_over_reference_median": engine_median / summary["median_pn_per_um"],
        "engine_median_inside_reference_iqr": summary["q1_pn_per_um"] <= engine_median <= summary["q3_pn_per_um"],
        "not_a_gate_verdict": True,
    }
    return result


def write_reference_csv(catalog: dict[str, Any], output: Path) -> None:
    reference = _record(catalog)
    summary = reference["summary"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "reference_id", "cell_line", "surface_state", "cell_cycle_state",
            "temperature_c", "measurement_operator", "unit", "median",
            "q1", "q3", "interval_kind", "n_cells", "citation_key", "doi",
            "knowledge_claim_id", "source_audit_required_verdict",
        ], lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            "reference_id": reference["reference_id"],
            **reference["context"],
            "measurement_operator": reference["measurement_operator"]["method"],
            "unit": summary["unit"],
            "median": summary["median_pn_per_um"],
            "q1": summary["q1_pn_per_um"],
            "q3": summary["q3_pn_per_um"],
            "interval_kind": summary["interval_kind"],
            "n_cells": summary["n_cells"],
            "citation_key": reference["evidence"]["citation_key"],
            "doi": reference["evidence"]["doi"],
            "knowledge_claim_id": reference["evidence"]["knowledge_claim_id"],
            "source_audit_required_verdict": reference["evidence"]["source_audit_required_verdict"],
        })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--kb", type=Path, required=True)
    parser.add_argument("--operator-registry", type=Path, required=True)
    parser.add_argument("--engine-record", type=Path, action="append", default=[])
    parser.add_argument("--reference-csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    catalog = json.loads(args.catalog.read_text(encoding="utf-8"))
    audit = audit_catalog(catalog, args.kb)
    registry = json.loads(args.operator_registry.read_text(encoding="utf-8"))
    records = [json.loads(path.read_text(encoding="utf-8")) for path in args.engine_record]
    report = compare_engine_records(records, catalog, audit, registry)
    report["catalog_audit"] = audit
    report["reference_head"] = infer_reference(_record(catalog))
    report["inputs"] = {
        "catalog": {"path": str(args.catalog), "sha256": _sha256(args.catalog)},
        "operator_registry": {"path": str(args.operator_registry), "sha256": _sha256(args.operator_registry)},
        "engine_records": [
            {"path": str(path), "sha256": _sha256(path)} for path in args.engine_record
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.reference_csv:
        write_reference_csv(catalog, args.reference_csv)
    print(json.dumps({
        "catalog_audit_passed": audit["passed"],
        "comparison_permitted": report["comparison_permitted"],
        "blockers": report["blockers"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
