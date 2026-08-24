import json
from pathlib import Path

from validate import validate_candidate, validate_experiment


HERE = Path(__file__).parent


def _load(name: str):
    return json.loads((HERE / "controls" / name).read_text(encoding="utf-8"))


def test_the_positive_control_is_valid():
    assert validate_experiment(_load("positive_minimal.json")) == []


def test_the_broken_evidence_hash_is_killed():
    errors = validate_experiment(_load("negative_missing_evidence_hash.json"))
    assert any("expected lowercase SHA-256" in error for error in errors)
    assert any("char_end precedes char_start" in error for error in errors)


def test_ambiguity_without_reason_and_candidates_is_killed():
    errors = validate_experiment(_load("negative_ambiguous_without_candidates.json"))
    assert any("ambiguous status requires reason and candidates" in error for error in errors)


def test_candidate_authority_cannot_be_promoted():
    record = {
        "schema": "aleph.external_training.annotation_candidate.v1",
        "authority_status": "accepted",
        "review_state": "gold",
    }
    errors = validate_candidate(record)
    assert any("must remain proposed" in error for error in errors)
    assert any("not human-verified gold" in error for error in errors)
