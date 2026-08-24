from __future__ import annotations

import gzip
import json
from copy import deepcopy
from pathlib import Path

import pytest

from .build_pilot import DOMAIN_TARGETS, PASS_TARGETS, YEAR_TARGETS, build, read_jsonl
from .audit_pilot import audit
from .evaluate_annotations import agreement, evaluate_machine, prepare_adjudication, validate_submission


ROOT = Path(__file__).resolve().parents[4]
CANDIDATES = ROOT / "aleph/outer/experiment_factory/goldset/candidates_300.jsonl"
RECORDS = ROOT / "data/external_training/experiment_factory/extraction/full/records"
RESULTS = Path(__file__).parent / "results"


def submission(source: str, slot: str, annotator: str, value: str = "fibroblast", adjudicated: bool = False) -> dict:
    return {
        "schema": "aleph.external_training.annotation_submission.v1",
        "authority_status": "proposed",
        "review_state": "adjudicated" if adjudicated else "single_annotator",
        "assignment_id": f"fixture:{source}:{slot}",
        "annotation_slot": "ADJUDICATION" if adjudicated else slot,
        "annotator_id": annotator,
        "source_family_id": source,
        "claims": [{
            "claim_type": "cell_type",
            "experiment_anchor": "paragraph:1",
            "value": value,
            "evidence_key": "a" * 64 + "|paragraph|p1|0|10|" + "b" * 64,
            "status": "adjudicated" if adjudicated else "reported",
            "notes": None,
        }],
        "completed_at": "2026-08-05T19:00:00+09:00",
    }


def test_frozen_results_have_required_coverage_and_no_component_split():
    summary = json.loads((RESULTS / "summary.json").read_text())
    assert summary["selected_sources"] == 50
    assert summary["assignments"] == 100
    assert summary["domain_counts"] == DOMAIN_TARGETS
    assert summary["year_stratum_counts"] == YEAR_TARGETS
    assert summary["acquisition_pass_counts"] == PASS_TARGETS
    assert summary["modality_coverage"]["missing"] == []
    assert summary["selection_audit"]["split_components"] == []
    assert summary["human_annotations_completed"] == 0
    assert summary["machine_preannotation_is_ground_truth"] is False


def test_blinded_forms_are_empty_and_machine_hidden():
    forms_a = read_jsonl(RESULTS / "forms_A.jsonl")
    forms_b = read_jsonl(RESULTS / "forms_B.jsonl")
    assert len(forms_a) == len(forms_b) == 50
    assert {row["source"]["source_family_id"] for row in forms_a} == {row["source"]["source_family_id"] for row in forms_b}
    for row in forms_a + forms_b:
        assert row["annotator_id"] is None
        assert row["claims"] == []
        assert row["review_state"] == "unassigned"
        assert row["blinding"]["machine_preannotation_hidden"] is True


def test_machine_packets_are_explicit_candidates_and_have_no_raw_text():
    with gzip.open(RESULTS / "machine_preannotations.jsonl.gz", "rt") as handle:
        packets = [json.loads(line) for line in handle]
    assert len(packets) == 50
    assert sum(len(row["experiment_records"]) for row in packets) == 996
    encoded = json.dumps(packets)
    assert '"human_verified": true' not in encoded.lower()
    assert '"review_state": "adjudicated"' not in encoded.lower()
    assert '"evidence_text"' not in encoded.lower()
    assert '"caption"' not in encoded.lower()
    assert all(row["review_state"] == "machine_candidate" for row in packets)


def test_generator_replay_is_byte_identical(tmp_path):
    args = type("Args", (), {"candidates": CANDIDATES, "records_root": RECORDS, "output": tmp_path})()
    build(args)
    expected = json.loads((RESULTS / "SHA256SUMS.json").read_text())
    actual = json.loads((tmp_path / "SHA256SUMS.json").read_text())
    assert actual == expected
    for name in expected:
        assert (tmp_path / name).read_bytes() == (RESULTS / name).read_bytes()


def test_submission_validation_and_pair_agreement():
    a = submission("doi:test", "A", "ann-a")
    b = submission("doi:test", "B", "ann-b")
    assert validate_submission(a) == []
    report = agreement([a], [b])
    assert report["sources"] == 1
    assert report["exact_source_claim_sets"] == 1
    assert report["per_field_exact_claim_agreement"]["cell_type"]["f1"] == 1.0


def test_same_annotator_is_rejected_and_no_auto_adjudication():
    a = submission("doi:test", "A", "ann-a")
    b = submission("doi:test", "B", "ann-a")
    with pytest.raises(ValueError, match="same annotator"):
        prepare_adjudication([a], [b])
    b["annotator_id"] = "ann-b"
    b["claims"][0]["value"] = "epithelial"
    packets = prepare_adjudication([a], [b])
    assert packets[0]["review_state"] == "double_adjudication_pending"
    assert len(packets[0]["disagreements"]) == 2
    assert packets[0]["adjudicated_claims"] == []
    assert packets[0]["automatic_adjudication_performed"] is False


def test_machine_benchmark_refuses_non_adjudicated_reference():
    machine = {"doi:test": {"atomic_claims": []}}
    not_gold = submission("doi:test", "A", "ann-a")
    with pytest.raises(ValueError, match="review_state"):
        evaluate_machine(machine, [not_gold])


def test_machine_benchmark_masks_ambiguous_and_reports_exact_metrics():
    gold = submission("doi:test", "ADJUDICATION", "adj-1", adjudicated=True)
    evidence = gold["claims"][0]["evidence_key"]
    machine = {"doi:test": {"atomic_claims": [
        {"claim_type": "cell_type", "experiment_anchor": "paragraph:1", "value": "fibroblast", "evidence_key": evidence, "ambiguous": False},
        {"claim_type": "cell_type", "experiment_anchor": "paragraph:1", "value": "wrong", "evidence_key": evidence, "ambiguous": True},
    ]}}
    report = evaluate_machine(machine, [gold])
    assert report["adjudicated_sources"] == 1
    assert report["ambiguous_machine_claims_masked"] == 1
    assert report["per_field_exact_metrics"]["cell_type"]["f1"] == 1.0


def test_adjudicated_claim_status_is_required():
    gold = submission("doi:test", "ADJUDICATION", "adj-1", adjudicated=True)
    broken = deepcopy(gold)
    broken["claims"][0]["status"] = "reported"
    assert "unadjudicated_claim_in_gold" in validate_submission(broken, adjudicated=True)


def test_committed_pilot_audit_passes():
    report = audit(RESULTS, CANDIDATES)
    assert report["status"] == "pass"
    assert report["manifest_mismatches"] == []
    assert report["forbidden_packet_keys"] == []
    assert report["split_leakage_components"] == []
    assert report["benchmark_metrics_available"] is False
