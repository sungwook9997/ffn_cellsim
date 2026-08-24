"""Controls for the accepted-step admissibility oracle.

The three PI-named controls are ``test_control_1_*`` (balanced-but-stretched), ``test_control_2_*``
(broken adjoint) and ``test_control_3_*`` (nonfinite).  Everything else guards the oracle itself.
"""

from __future__ import annotations

import hashlib
import inspect
import math
import subprocess
import sys

import numpy as np
import pytest

from aleph.virtual_cell.admissibility import (
    ADMISSIBILITY_CLAUSE_ORDER,
    EXPANSION_ROUTES,
    AdmissibilityClause,
    AdmissibilityVerdict,
    AttemptedStepRecord,
    ClauseStatus,
    ConnectorForcePair,
    RollbackChannel,
    StateDigestSet,
    adjoint_closure_only,
    digest_array,
    digest_named_arrays,
    evaluate_admissibility,
    verify_rollback,
)
from aleph.virtual_cell.contracts import EvidenceSource
from aleph.virtual_cell.step_receipt import StepReceipt, state_coverage_sha256

# --- tolerances -------------------------------------------------------------------------------
# These are TEST FIXTURE values, not a ratified gate contract.  The oracle refuses to supply a
# default precisely because choosing them is a PI decision (see the module docstring).
REL_TOL = 1.0e-9
ABS_FLOOR_PN = 1.0e-12

# The inner-solve numbers below reproduce the magnitudes recorded in the archived full-native
# artifact (`STATE.md` tier-(a): `inner_converged: false`, `residual_start = residual_end = 2634.12`).
# They are used here as FIXTURE magnitudes for a CPU control; nothing here re-measures the engine.
NATIVE_STALLED_RESIDUAL_PN = 2634.12
INNER_TOLERANCE_PN = 1.0e-3


def _channel(
    prefix: str,
    *,
    clock_step: int = 7,
    clock_time_s: float = 0.875,
    n_edges: int = 3,
) -> dict[str, object]:
    return {
        "state": {
            "pos_um": np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float64) + len(prefix),
            "binding_soa": np.array([1, 0, 1, 1], dtype=np.int32),
        },
        "rng": {"counter": np.array([len(prefix)], dtype=np.uint64)},
        "topology": {
            "edges": np.array([[i, i + 1] for i in range(n_edges)], dtype=np.int32),
        },
        "clock": {
            "accepted_step": np.array([clock_step], dtype=np.int64),
            "accepted_time_s": np.array([clock_time_s], dtype=np.float64),
        },
    }


def _digests(
    prefix: str,
    *,
    clock_step: int = 7,
    clock_time_s: float = 0.875,
    n_edges: int = 3,
) -> StateDigestSet:
    channels = _channel(prefix, clock_step=clock_step, clock_time_s=clock_time_s, n_edges=n_edges)
    return StateDigestSet.from_channels(**channels)  # type: ignore[arg-type]


_BEFORE = _digests("before")
# The candidate mutated all four channels: it moved nodes, drew randoms, bound one more crosslink,
# and provisionally advanced the clock.  A rollback audit over a candidate that mutated nothing
# would be vacuous, so the fixture deliberately mutates every channel.
_AFTER = _digests("after-candidate", clock_step=8, clock_time_s=1.0, n_edges=4)


def _balanced_pair(connector_id: str, force: np.ndarray) -> ConnectorForcePair:
    """A correctly wired connector: the B side is the exact negation, perturbed by one ULP."""
    reaction = -force.copy()
    reaction[0] = np.nextafter(reaction[0], math.inf)
    return ConnectorForcePair(connector_id=connector_id, force_a_pn=force, force_b_pn=reaction)


def _record(**overrides: object) -> AttemptedStepRecord:
    fields: dict[str, object] = {
        "attempted_step_index": 41,
        "attempted_dt_s": 1.0e-4,
        "connector_pairs": (
            _balanced_pair("nmii_sf_motor", np.array([1234.5, -678.9, 42.0])),
            _balanced_pair("erm_cortex", np.array([3.5, 0.25, -1.75])),
        ),
        "inner_residual_pn": 1.0e-6,
        "inner_tolerance_pn": INNER_TOLERANCE_PN,
        "inner_iterations": 120,
        "inner_iteration_budget": 400,
        "candidate_state": {
            "pos_um": np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]], dtype=np.float64),
            "force_pn": np.array([[1.0, -1.0, 0.5]], dtype=np.float64),
        },
        "digests_before": _BEFORE,
        "digests_after": _AFTER,
        "rollback_digests": None,
    }
    fields.update(overrides)
    return AttemptedStepRecord(**fields)  # type: ignore[arg-type]


def _verdict(record: AttemptedStepRecord) -> AdmissibilityVerdict:
    return evaluate_admissibility(
        record, adjoint_rel_tol=REL_TOL, adjoint_abs_floor_pn=ABS_FLOOR_PN
    )


def _incumbent(record: AttemptedStepRecord) -> bool:
    return adjoint_closure_only(record, adjoint_rel_tol=REL_TOL, adjoint_abs_floor_pn=ABS_FLOOR_PN)


# --- baseline ---------------------------------------------------------------------------------


def test_a_wired_and_converged_step_passes_every_clause() -> None:
    verdict = _verdict(_record())
    assert verdict.admissible
    assert verdict.rejected_by == ()
    assert tuple(c.clause for c in verdict.clauses) == ADMISSIBILITY_CLAUSE_ORDER
    assert all(c.status is ClauseStatus.PASS for c in verdict.clauses)
    assert _incumbent(_record())


# --- CONTROL 1: balanced-but-stretched --------------------------------------------------------


def test_control_1_balanced_but_stretched_is_accepted_by_adjoint_only_and_rejected_here() -> None:
    """Force pairs cancel to round-off while the mechanical solve is 6 decades off tolerance."""
    record = _record(
        inner_residual_pn=NATIVE_STALLED_RESIDUAL_PN,
        inner_iterations=400,
    )

    # The incumbent force-only predicate ACCEPTS: Newton's third law holds by construction.
    assert _incumbent(record) is True

    verdict = _verdict(record)
    assert verdict.admissible is False
    assert verdict.rejected_by == (AdmissibilityClause.INNER_CONVERGED,)

    adjoint = verdict.outcome(AdmissibilityClause.ADJOINT_CLOSURE)
    convergence = verdict.outcome(AdmissibilityClause.INNER_CONVERGED)
    assert adjoint.status is ClauseStatus.PASS
    assert convergence.status is ClauseStatus.FAIL
    assert convergence.measured == pytest.approx(NATIVE_STALLED_RESIDUAL_PN)
    assert convergence.threshold == pytest.approx(INNER_TOLERANCE_PN)
    # The adjoint residual really is round-off small against the pair's own force traffic.
    assert adjoint.measured is not None and adjoint.threshold is not None
    assert adjoint.measured <= adjoint.threshold
    assert adjoint.measured < 1.0e-9
    assert "native cost failure" in verdict.expansions[0]


def test_control_1_numbers_are_the_ones_the_decision_card_quotes() -> None:
    """Pin the measured magnitudes the PI decision card cites, so the card cannot drift."""
    record = _record(inner_residual_pn=NATIVE_STALLED_RESIDUAL_PN, inner_iterations=400)
    worst = max(record.connector_pairs, key=lambda pair: pair.residual_pn)
    assert worst.connector_id == "nmii_sf_motor"
    assert worst.residual_pn == pytest.approx(2.2737367544323206e-13, rel=1e-9)
    assert worst.magnitude_scale_pn == pytest.approx(2818.9781552896075, rel=1e-9)
    ratio = worst.residual_pn / worst.magnitude_scale_pn
    assert ratio < 1.0e-15
    assert record.inner_residual_pn / record.inner_tolerance_pn == pytest.approx(2634120.0)


# --- CONTROL 2: broken adjoint ----------------------------------------------------------------


def test_control_2_broken_adjoint_one_sided_injection_rejects() -> None:
    """One connector scatters onto A without the equal-and-opposite reaction on B."""
    force = np.array([1234.5, -678.9, 42.0])
    one_sided = ConnectorForcePair(
        connector_id="nmii_sf_motor",
        force_a_pn=force,
        force_b_pn=-force + np.array([12.0, 0.0, 0.0]),  # 12 pN of the pair never lands
    )
    record = _record(
        connector_pairs=(one_sided, _balanced_pair("erm_cortex", np.array([3.5, 0.25, -1.75])))
    )

    verdict = _verdict(record)
    assert verdict.admissible is False
    assert verdict.rejected_by == (AdmissibilityClause.ADJOINT_CLOSURE,)
    outcome = verdict.outcome(AdmissibilityClause.ADJOINT_CLOSURE)
    assert outcome.measured == pytest.approx(12.0)
    assert "nmii_sf_motor" in outcome.detail
    # The incumbent global gate also catches this one — it is the positive control the SF lane ran.
    assert _incumbent(record) is False
    assert "missing component" in verdict.expansions[0]


def test_control_2b_global_gate_misses_what_the_per_connector_clause_catches() -> None:
    """Two connectors with equal and opposite one-sided errors cancel in a GLOBAL resultant."""
    leak = np.array([12.0, 0.0, 0.0])
    a_force = np.array([1234.5, -678.9, 42.0])
    b_force = np.array([-800.0, 15.0, 7.5])
    record = _record(
        connector_pairs=(
            ConnectorForcePair("nmii_sf_motor", a_force, -a_force + leak),
            ConnectorForcePair("erm_cortex", b_force, -b_force - leak),
        )
    )

    assert _incumbent(record) is True  # global |Σreaction + Σtraction| == 0
    verdict = _verdict(record)
    assert verdict.admissible is False
    assert AdmissibilityClause.ADJOINT_CLOSURE in verdict.rejected_by


# --- CONTROL 3: nonfinite ---------------------------------------------------------------------


def test_control_3_nan_residual_rejects_before_any_ordered_comparison() -> None:
    """A NaN comparison is False both ways; finiteness must therefore short-circuit."""
    # The trap this clause exists to close, stated as an executable fact:
    assert not (float("nan") > INNER_TOLERANCE_PN)
    assert not (float("nan") <= INNER_TOLERANCE_PN)

    record = _record(inner_residual_pn=float("nan"))
    verdict = _verdict(record)
    assert verdict.admissible is False
    assert verdict.rejected_by[0] is AdmissibilityClause.STATE_FINITE
    for clause in ADMISSIBILITY_CLAUSE_ORDER[1:]:
        outcome = verdict.outcome(clause)
        assert outcome.status is ClauseStatus.NOT_EVALUATED
        assert outcome.measured is None  # nothing was compared at all
        assert outcome.passed is False
    assert "inner_residual_pn" in verdict.outcome(AdmissibilityClause.STATE_FINITE).detail


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_control_3_nonfinite_state_is_rejected_and_is_invisible_to_the_incumbent(
    bad: float,
) -> None:
    record = _record(
        candidate_state={
            "pos_um": np.array([[0.1, 0.2, bad], [0.4, 0.5, 0.6]], dtype=np.float64),
            "force_pn": np.array([[1.0, -1.0, 0.5]], dtype=np.float64),
        }
    )
    assert _incumbent(record) is True  # the force channel is clean; the gate never sees the state
    verdict = _verdict(record)
    assert verdict.admissible is False
    assert verdict.rejected_by[0] is AdmissibilityClause.STATE_FINITE
    assert "candidate_state.pos_um" in verdict.outcome(AdmissibilityClause.STATE_FINITE).detail


def test_control_3_nonfinite_connector_force_is_named() -> None:
    record = _record(
        connector_pairs=(
            ConnectorForcePair(
                "nmii_sf_motor",
                np.array([1.0, 2.0, float("nan")]),
                np.array([-1.0, -2.0, 0.0]),
            ),
        )
    )
    verdict = _verdict(record)
    assert verdict.rejected_by[0] is AdmissibilityClause.STATE_FINITE
    assert "connector.nmii_sf_motor" in verdict.outcome(AdmissibilityClause.STATE_FINITE).detail


# --- scale safety -----------------------------------------------------------------------------


def test_adjoint_closure_is_scale_safe_where_an_absolute_threshold_is_not() -> None:
    """A single absolute pN threshold misjudges connectors whose traffic differs by decades."""
    heavy = ConnectorForcePair(
        "sf_motor_heavy",
        np.array([1.0e6, 0.0, 0.0]),
        np.array([-1.0e6 + 1.0e-4, 0.0, 0.0]),  # 1e-4 pN residual on 2e6 pN of traffic
    )
    light = ConnectorForcePair(
        "linc_light",
        np.array([1.0e-3, 0.0, 0.0]),
        np.array([-1.0e-3 + 1.0e-5, 0.0, 0.0]),  # 1e-5 pN residual on 2e-3 pN of traffic
    )
    assert heavy.residual_pn == pytest.approx(1.0e-4)
    assert light.residual_pn == pytest.approx(1.0e-5)

    # An absolute threshold of 1e-4.5 pN would pass the LIGHT pair (1e-5 < 3.2e-5) and fail the
    # HEAVY one, which is backwards: heavy is closed to 5e-11 relative, light only to 5.0e-3.
    absolute_threshold = 10.0**-4.5
    assert light.residual_pn < absolute_threshold < heavy.residual_pn
    assert heavy.residual_pn / heavy.magnitude_scale_pn == pytest.approx(5.0e-11, rel=1e-6)
    assert light.residual_pn / light.magnitude_scale_pn == pytest.approx(
        5.025125628140704e-3, rel=1e-9
    )

    heavy_verdict = _verdict(_record(connector_pairs=(heavy,)))
    light_verdict = _verdict(_record(connector_pairs=(light,)))
    assert heavy_verdict.admissible is True
    assert light_verdict.admissible is False


def test_absolute_floor_covers_the_near_zero_pair() -> None:
    """A pair with no force traffic has no relative scale; the floor decides."""
    idle = ConnectorForcePair("idle", np.zeros(3), np.zeros(3))
    assert idle.magnitude_scale_pn == 0.0
    assert _verdict(_record(connector_pairs=(idle,))).admissible is True


def test_tolerances_are_required_and_never_defaulted() -> None:
    for function in (evaluate_admissibility, adjoint_closure_only):
        signature = inspect.signature(function)
        for name in ("adjoint_rel_tol", "adjoint_abs_floor_pn"):
            parameter = signature.parameters[name]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
            assert parameter.default is inspect.Parameter.empty
    with pytest.raises(TypeError):
        evaluate_admissibility(_record())  # type: ignore[call-arg]
    with pytest.raises(ValueError, match="nonzero relative tolerance or absolute floor"):
        evaluate_admissibility(_record(), adjoint_rel_tol=0.0, adjoint_abs_floor_pn=0.0)


# --- iteration budget -------------------------------------------------------------------------


def test_iteration_budget_overrun_is_its_own_clause() -> None:
    record = _record(inner_iterations=401, inner_iteration_budget=400)
    verdict = _verdict(record)
    assert verdict.rejected_by == (AdmissibilityClause.ITERATION_BUDGET,)
    assert verdict.outcome(AdmissibilityClause.ITERATION_BUDGET).measured == pytest.approx(401.0)


# --- acceptance is not stationarity ------------------------------------------------------------


def test_acceptance_is_not_stationarity_a_far_from_rest_step_is_admissible() -> None:
    """A large, non-decaying, wholly changed candidate is admissible if every clause holds."""
    record = _record(
        candidate_state={
            "pos_um": np.array([[900.0, -900.0, 900.0]], dtype=np.float64),
            "force_pn": np.array([[5.0e4, -5.0e4, 5.0e4]], dtype=np.float64),
        },
        digests_after=_digests("limit-cycle", clock_step=8, clock_time_s=1.0),
    )
    verdict = _verdict(record)
    assert verdict.admissible is True
    assert record.digests_after.state != record.digests_before.state
    # No clause references the previous state at all.
    assert {c.clause for c in verdict.clauses} == set(ADMISSIBILITY_CLAUSE_ORDER)


# --- rollback ---------------------------------------------------------------------------------


def test_bit_exact_rollback_of_all_four_channels_verifies() -> None:
    record = _record(
        inner_residual_pn=NATIVE_STALLED_RESIDUAL_PN,
        rollback_digests=_BEFORE,
    )
    audit = verify_rollback(record)
    assert audit.bit_exact
    assert audit.leaked == ()
    assert set(audit.restored) == set(RollbackChannel)
    assert set(audit.mutated_by_candidate) == set(RollbackChannel)
    assert audit.vacuous is False

    verdict = _verdict(record)
    assert verdict.admissible is False
    assert verdict.transaction_sound is True
    assert verdict.outcome(AdmissibilityClause.ROLLBACK_BIT_EXACT).status is ClauseStatus.PASS


@pytest.mark.parametrize(
    ("leaked_channel", "rolled_back"),
    [
        # rollback restores state/rng/topology but silently leaves the clock advanced
        (
            RollbackChannel.CLOCK,
            StateDigestSet(
                state=_BEFORE.state,
                rng=_BEFORE.rng,
                topology=_BEFORE.topology,
                clock=_AFTER.clock,
            ),
        ),
        # rollback drops a topology edge
        (
            RollbackChannel.TOPOLOGY,
            StateDigestSet(
                state=_BEFORE.state,
                rng=_BEFORE.rng,
                topology=digest_named_arrays(
                    {"edges": np.array([[0, 1], [1, 2]], dtype=np.int32)}  # the [2, 3] edge is gone
                ),
                clock=_BEFORE.clock,
            ),
        ),
        # rollback leaves the RNG stream advanced
        (
            RollbackChannel.RNG,
            StateDigestSet(
                state=_BEFORE.state,
                rng=_AFTER.rng,
                topology=_BEFORE.topology,
                clock=_BEFORE.clock,
            ),
        ),
        # rollback misses a state array
        (
            RollbackChannel.STATE,
            StateDigestSet(
                state=_AFTER.state,
                rng=_BEFORE.rng,
                topology=_BEFORE.topology,
                clock=_BEFORE.clock,
            ),
        ),
    ],
)
def test_rollback_restoring_three_of_four_is_a_failure_and_names_the_leak(
    leaked_channel: RollbackChannel, rolled_back: StateDigestSet
) -> None:
    record = _record(inner_residual_pn=NATIVE_STALLED_RESIDUAL_PN, rollback_digests=rolled_back)
    audit = verify_rollback(record)
    assert audit.bit_exact is False
    assert audit.leaked == (leaked_channel,)
    assert len(audit.restored) == 3
    assert leaked_channel.value in audit.detail
    assert "3/4" not in audit.detail and "1/4" in audit.detail

    verdict = _verdict(record)
    assert verdict.transaction_sound is False
    clause = verdict.outcome(AdmissibilityClause.ROLLBACK_BIT_EXACT)
    assert clause.status is ClauseStatus.FAIL
    assert leaked_channel.value in clause.detail
    assert clause.measured == pytest.approx(3.0)
    assert any("missing component" in expansion for expansion in verdict.expansions)


def test_vacuous_rollback_is_flagged_because_it_proves_nothing() -> None:
    record = _record(
        inner_residual_pn=NATIVE_STALLED_RESIDUAL_PN,
        digests_after=_BEFORE,  # the candidate never mutated anything
        rollback_digests=_BEFORE,
    )
    audit = verify_rollback(record)
    assert audit.bit_exact is True
    assert audit.vacuous is True
    assert "vacuous" in audit.detail


def test_rollback_verification_requires_post_rollback_digests() -> None:
    with pytest.raises(ValueError, match="post-rollback digests"):
        verify_rollback(_record())


# --- byte-level digest semantics ---------------------------------------------------------------


def test_digest_is_byte_exact_where_float_equality_lies() -> None:
    positive_zero = np.array([0.0], dtype=np.float64)
    negative_zero = np.array([-0.0], dtype=np.float64)
    assert bool(positive_zero[0] == negative_zero[0])  # float compare says "restored"
    assert digest_array(positive_zero) != digest_array(negative_zero)  # the bytes say otherwise

    nan_a = np.array([np.nan], dtype=np.float64)
    nan_b = nan_a.copy()
    assert not bool(nan_a[0] == nan_b[0])  # float compare says "leaked"
    assert digest_array(nan_a) == digest_array(nan_b)  # the bytes say otherwise


def test_digest_separates_dtype_and_shape() -> None:
    assert digest_array(np.array([1, 2, 3], dtype=np.int32)) != digest_array(
        np.array([1, 2, 3], dtype=np.int64)
    )
    assert digest_array(np.array([[1, 2], [3, 4]], dtype=np.int32)) != digest_array(
        np.array([1, 2, 3, 4], dtype=np.int32)
    )


def test_digest_rejects_object_arrays() -> None:
    with pytest.raises(TypeError, match="object arrays"):
        digest_array(np.array([object()], dtype=object))


def test_named_array_digest_is_order_independent_but_name_sensitive() -> None:
    left = digest_named_arrays({"a": np.zeros(2), "b": np.ones(2)})
    right = digest_named_arrays({"b": np.ones(2), "a": np.zeros(2)})
    renamed = digest_named_arrays({"a": np.zeros(2), "c": np.ones(2)})
    assert left == right
    assert left != renamed


def test_named_array_digest_rejects_normalized_name_collisions() -> None:
    with pytest.raises(ValueError, match="collide"):
        digest_named_arrays({"pos": np.zeros(2), " pos ": np.ones(2)})


# --- authority and composition ------------------------------------------------------------------


def test_verdict_is_an_unverified_observer_and_cannot_claim_native_evidence() -> None:
    verdict = _verdict(_record())
    assert verdict.authority == "unverified-observer"
    assert verdict.evidence_source is EvidenceSource.ANALYTIC_ORACLE
    assert verdict.evidence_source is not EvidenceSource.NATIVE_ACCEPTED
    with pytest.raises(ValueError, match="unverified-observer"):
        AdmissibilityVerdict(
            attempted_step_index=0,
            admissible=True,
            clauses=verdict.clauses,
            rejected_by=(),
            rollback=None,
            expansions=(),
            authority="native-accepted",
        )


def test_verdict_clause_vector_composes_with_a_step_receipt() -> None:
    record = _record(inner_residual_pn=NATIVE_STALLED_RESIDUAL_PN, rollback_digests=_BEFORE)
    verdict = _verdict(record)
    coverage = ("component_state", "rng_stream", "topology_free_list", "accepted_clock")
    receipt = StepReceipt(
        run_id="admissibility-oracle",
        decision_contract_sha256=hashlib.sha256(b"pending-pi-signature").hexdigest(),
        source_run_manifest_sha256=hashlib.sha256(b"cpu-sandbox").hexdigest(),
        build_stamp_sha256=hashlib.sha256(b"build").hexdigest(),
        attempted_step_index=record.attempted_step_index,
        accepted_step_before=12,
        accepted_step_after=12,
        accepted_time_before_s=0.5,
        accepted_time_after_s=0.5,
        attempted_dt_s=record.attempted_dt_s,
        accepted=verdict.admissible,
        predicate_results=dict(verdict.predicate_results()),
        predicate_results_sha256=verdict.predicate_results_sha256(),
        state_coverage=coverage,
        state_coverage_sha256=state_coverage_sha256(coverage),
        state_hash_before=_BEFORE.state,
        state_hash_after=_BEFORE.state,
        rng_hash_before=_BEFORE.rng,
        rng_hash_after=_BEFORE.rng,
        rejected_state_restored=True,
    )
    assert receipt.accepted is False
    assert receipt.predicate_results["inner_converged"] is False
    assert receipt.predicate_results["adjoint_closure"] is True
    assert receipt.authority == "unverified-observer"


def test_skipped_clauses_are_false_in_the_predicate_vector() -> None:
    verdict = _verdict(_record(inner_residual_pn=float("nan")))
    results = verdict.predicate_results()
    assert results["state_finite"] is False
    assert results["adjoint_closure"] is False
    assert results["inner_converged"] is False


def test_every_clause_has_a_declared_expansion_route() -> None:
    for clause in AdmissibilityClause:
        assert clause in EXPANSION_ROUTES
        assert EXPANSION_ROUTES[clause].strip()


# --- record hygiene ------------------------------------------------------------------------------


def test_record_is_immutable_and_copies_its_arrays() -> None:
    positions = np.array([[0.1, 0.2, 0.3]], dtype=np.float64)
    record = _record(candidate_state={"pos_um": positions})
    positions[0, 0] = 99.0
    assert record.candidate_state["pos_um"][0, 0] == pytest.approx(0.1)
    with pytest.raises(ValueError):
        record.candidate_state["pos_um"][0, 0] = 5.0
    with pytest.raises((AttributeError, TypeError)):
        record.inner_residual_pn = 0.0  # type: ignore[misc]


def test_record_rejects_duplicate_connector_ids() -> None:
    pair = _balanced_pair("nmii_sf_motor", np.array([1.0, 0.0, 0.0]))
    with pytest.raises(ValueError, match="unique"):
        _record(connector_pairs=(pair, pair))


def test_record_requires_at_least_one_connector_pair() -> None:
    with pytest.raises(ValueError, match="at least one connector"):
        _record(connector_pairs=())


def test_record_accepts_a_nonfinite_residual_so_the_clause_can_catch_it() -> None:
    record = _record(inner_residual_pn=float("nan"))
    assert math.isnan(record.inner_residual_pn)
    with pytest.raises(ValueError, match="inner_tolerance_pn"):
        _record(inner_tolerance_pn=float("nan"))


def test_connector_pair_accepts_per_node_force_arrays() -> None:
    per_node = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.5, 0.5, 0.0]])
    pair = ConnectorForcePair("split_scatter", per_node, -per_node)
    assert pair.residual_pn == pytest.approx(0.0)
    # Σ|f_i| is the Higham scale and is NOT |Σf_i| (which would be 0.5·√2 here).
    assert pair.magnitude_scale_pn == pytest.approx(2.0 * (1.0 + 1.0 + math.sqrt(0.5)))


# --- CPU-only ------------------------------------------------------------------------------------


def test_importing_the_oracle_does_not_import_or_initialise_warp() -> None:
    probe = (
        "import sys;"
        "import aleph.virtual_cell.admissibility as m;"
        "assert 'warp' not in sys.modules, sorted(k for k in sys.modules if 'warp' in k);"
        "print(m.__name__)"
    )
    completed = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "aleph.virtual_cell.admissibility" in completed.stdout
