from __future__ import annotations

import hashlib

import pytest

from aleph.virtual_cell.step_receipt import (
    StepReceipt,
    predicate_vector_sha256,
    state_coverage_sha256,
    validate_step_receipt_chain,
)

_COVERAGE = (
    "component_state",
    "connector_binding_soa",
    "topology_free_list",
    "field_state",
    "accepted_clock",
)


def _hash(label: str) -> str:
    return hashlib.sha256(label.encode()).hexdigest()


def _receipt(**overrides: object) -> StepReceipt:
    predicates = overrides.get(
        "predicate_results",
        {"inner_converged": True, "balance_ok": True},
    )
    coverage = overrides.get("state_coverage", _COVERAGE)
    fields: dict[str, object] = {
        "run_id": "run-1",
        "decision_contract_sha256": _hash("external-runtime-contract"),
        "source_run_manifest_sha256": _hash("run-manifest"),
        "build_stamp_sha256": _hash("build"),
        "attempted_step_index": 9,
        "accepted_step_before": 4,
        "accepted_step_after": 5,
        "accepted_time_before_s": 0.375,
        "accepted_time_after_s": 0.5,
        "attempted_dt_s": 0.125,
        "accepted": True,
        "predicate_results": predicates,
        "predicate_results_sha256": predicate_vector_sha256(predicates),  # type: ignore[arg-type]
        "state_coverage": coverage,
        "state_coverage_sha256": state_coverage_sha256(coverage),  # type: ignore[arg-type]
        "state_hash_before": _hash("state-before"),
        "state_hash_after": _hash("state-after"),
        "rng_hash_before": _hash("rng-before"),
        "rng_hash_after": _hash("rng-after"),
        "rejected_state_restored": None,
    }
    fields.update(overrides)
    if "predicate_results" in overrides and "predicate_results_sha256" not in overrides:
        fields["predicate_results_sha256"] = predicate_vector_sha256(
            overrides["predicate_results"]  # type: ignore[arg-type]
        )
    if "state_coverage" in overrides and "state_coverage_sha256" not in overrides:
        fields["state_coverage_sha256"] = state_coverage_sha256(
            overrides["state_coverage"]  # type: ignore[arg-type]
        )
    return StepReceipt(**fields)  # type: ignore[arg-type]


def test_accepted_receipt_advances_accepted_clock_once() -> None:
    receipt = _receipt()
    assert receipt.accepted_step_after == 5
    assert receipt.accepted_time_after_s == pytest.approx(0.5)
    assert receipt.authority == "unverified-observer"


def test_rejected_receipt_requires_claimed_digest_rollback() -> None:
    receipt = _receipt(
        accepted=False,
        accepted_step_after=4,
        accepted_time_after_s=0.375,
        state_hash_after=_hash("state-before"),
        rng_hash_after=_hash("rng-before"),
        rejected_state_restored=True,
        predicate_results={"inner_converged": False, "balance_ok": True},
    )
    assert not receipt.accepted


@pytest.mark.parametrize(
    ("field", "value", "match"),
    [
        ("accepted_step_after", 5, "must not advance accepted_step"),
        ("accepted_time_after_s", 0.5, "must not advance accepted time"),
        ("rejected_state_restored", False, "restored-state"),
        ("state_hash_after", _hash("changed"), "state digest"),
        ("rng_hash_after", _hash("changed"), "RNG digest"),
    ],
)
def test_rejected_receipt_rejects_inconsistent_transaction_claim(
    field: str, value: object, match: str
) -> None:
    overrides: dict[str, object] = {
        "accepted": False,
        "accepted_step_after": 4,
        "accepted_time_after_s": 0.375,
        "state_hash_after": _hash("state-before"),
        "rng_hash_after": _hash("rng-before"),
        "rejected_state_restored": True,
    }
    overrides[field] = value
    with pytest.raises(ValueError, match=match):
        _receipt(**overrides)


def test_receipt_does_not_infer_acceptance_from_predicate_names() -> None:
    receipt = _receipt(predicate_results={"arbitrary_future_predicate": False})
    assert receipt.accepted
    assert receipt.authority == "unverified-observer"


def test_receipt_rejects_non_sha256_digests() -> None:
    with pytest.raises(ValueError, match="SHA-256"):
        _receipt(state_hash_before="state-before")


def test_predicate_vector_is_bound_to_receipt() -> None:
    with pytest.raises(ValueError, match="does not bind"):
        _receipt(predicate_results_sha256=_hash("wrong"))


def test_predicate_names_cannot_collide_after_normalization() -> None:
    with pytest.raises(ValueError, match="collide"):
        _receipt(predicate_results={"balance_ok": True, " balance_ok ": False})


def test_state_coverage_manifest_is_bound_to_receipt() -> None:
    with pytest.raises(ValueError, match="does not bind"):
        _receipt(state_coverage_sha256=_hash("incomplete"))


def test_accepted_receipt_cannot_claim_rejected_state_restoration() -> None:
    with pytest.raises(ValueError, match="must not claim"):
        _receipt(rejected_state_restored=True)


def test_receipt_chain_checks_attempt_clock_state_and_rng_continuity() -> None:
    first = _receipt()
    second = _receipt(
        attempted_step_index=10,
        accepted_step_before=5,
        accepted_step_after=5,
        accepted_time_before_s=0.5,
        accepted_time_after_s=0.5,
        accepted=False,
        state_hash_before=_hash("state-after"),
        state_hash_after=_hash("state-after"),
        rng_hash_before=_hash("rng-after"),
        rng_hash_after=_hash("rng-after"),
        rejected_state_restored=True,
    )
    assert validate_step_receipt_chain((first, second)) == (first, second)


def test_receipt_chain_detects_missing_attempt() -> None:
    first = _receipt()
    second = _receipt(
        attempted_step_index=11,
        accepted_step_before=5,
        accepted_step_after=6,
        accepted_time_before_s=0.5,
        accepted_time_after_s=0.625,
        state_hash_before=_hash("state-after"),
        rng_hash_before=_hash("rng-after"),
    )
    with pytest.raises(ValueError, match="missing, duplicate, or reordered"):
        validate_step_receipt_chain((first, second))
