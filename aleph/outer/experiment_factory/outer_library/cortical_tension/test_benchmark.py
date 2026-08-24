from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from aleph.outer.experiment_factory.outer_library.cortical_tension.benchmark import (
    audit_catalog,
    compare_engine_records,
    infer_reference,
)


HERE = Path(__file__).resolve().parent
OUTER = HERE.parent
REPO = HERE.parents[4]
CATALOG = json.loads((HERE / "reference_catalog.json").read_text(encoding="utf-8"))
REGISTRY = json.loads((OUTER / "observation_operator_registry.json").read_text(encoding="utf-8"))
KB = REPO / "aleph/outputs/tag_kb/kb.duckdb"


def test_catalog_resolves_only_source_audit_ok_evidence():
    result = audit_catalog(CATALOG, KB)
    assert result["passed"] is True
    assert result["resolved_evidence"]["source_audit_verdict"] == "OK"
    assert result["resolved_evidence"]["validation_gate_status"] == "draft"
    assert result["may_score_validation_gate"] is False


def test_reference_head_preserves_context_and_non_authority():
    result = infer_reference(CATALOG["records"][0])
    posterior = result["state_posterior"][0]
    assert posterior["estimate"] == 270.0
    assert posterior["interval"] == [180.0, 400.0]
    assert posterior["conditioning"]["cell_line"] == "MCF-7"
    assert result["uncertainty"]["calibrated_predictive_distribution"] is False
    assert result["may_emit_Aleph_parameter"] is False
    assert result["may_score_validation_gate"] is False


def test_existing_force_accepted_record_is_value_redacted():
    record = json.loads((
        REPO / "aleph/outputs/ac/engine_cortex_ab/settled_on/record.json"
    ).read_text(encoding="utf-8"))
    result = compare_engine_records(
        [record], CATALOG, audit_catalog(CATALOG, KB), REGISTRY
    )
    assert result["comparison_permitted"] is False
    assert result["engine_summary"] is None
    assert result["external_reference"] is None
    assert result["redaction"]["engine_magnitudes_and_ratios_redacted"] is True
    blockers = set(result["blockers"])
    assert "engine_record[0]:physical_acceptance_gate_passed" in blockers
    assert "engine_record[0]:quantitative_claim_not_blocked" in blockers
    assert "observation_operator_not_native_validated_and_geometry_matched" in blockers


def _synthetic_record(seed: int, value: float) -> dict:
    return {
        "schema": "ac.engine.observe/run-record@2",
        "run_label": f"synthetic comparator control seed {seed}",
        "build": {"commit": "0123456789abcdef"},
        "device": "cuda:0",
        "quantitative_claim_status": "PERMITTED_BY_SYNTHETIC_TEST_ONLY",
        "config": {"seed": seed},
        "census": {"fraction_of_native": 1.0},
        "comparison_context": {
            "target_reference_id": "hosseini2020_mcf7_suspended_interphase_afm",
            "cell_line": "MCF-7",
            "surface_state": "suspended_rounded",
            "cell_cycle_state": "interphase",
            "temperature_c": 37.0,
            "measurement_operator": "gamma_total_method_of_planes",
            "full_native_compartments": True,
            "physical_acceptance_gate_passed": True,
            "time_step_independence_passed": True,
        },
        "measurements": {"stationarity": {"gamma_total_pn_per_um": {
            "verdict": "STATIONARY", "mean": value,
        }}},
    }


def test_three_seed_ready_control_computes_descriptive_metrics_not_a_gate():
    registry = copy.deepcopy(REGISTRY)
    entry = next(item for item in registry["entries"] if item["external_observable"] == "effective_laplace_surface_tension")
    entry.update({
        "status": "native_validated",
        "geometry_matched": True,
        "unit_conversion_verified": True,
        "validation_report_sha256": "a" * 64,
    })
    records = [_synthetic_record(0, 240.0), _synthetic_record(1, 270.0), _synthetic_record(2, 300.0)]
    result = compare_engine_records(
        records, CATALOG, audit_catalog(CATALOG, KB), registry
    )
    assert result["comparison_permitted"] is True
    assert result["engine_summary"]["median_pn_per_um"] == 270.0
    assert result["engine_summary"]["between_seed_sd_pn_per_um"] == pytest.approx(30.0)
    assert result["descriptive_position"]["engine_median_inside_reference_iqr"] is True
    assert result["validation_gate_verdict"] is None
    assert result["may_emit_Aleph_parameter"] is False


def test_protocol_mismatch_remains_redacted_even_with_three_seeds():
    records = [_synthetic_record(seed, 270.0) for seed in range(3)]
    records[1]["comparison_context"]["cell_line"] = "MDA-MB-231"
    result = compare_engine_records(
        records, CATALOG, audit_catalog(CATALOG, KB), REGISTRY
    )
    assert result["comparison_permitted"] is False
    assert "engine_record[1]:cell_line_matches" in result["blockers"]
    assert result["engine_summary"] is None
