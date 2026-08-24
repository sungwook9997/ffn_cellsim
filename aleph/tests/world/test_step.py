"""PHASE 4 gates — the step runs, and nothing invents a verdict for it.

These are the properties that stop the two failures this engine has already shipped: a gate that
cannot fail (`balance_ok`, a Higham rounding bound) and one that cannot succeed (`descent_ratio`,
read post-rollback). Both were one-sided, and both passed review, so the gates below are written
against the SHAPE of that mistake rather than against either instance.
"""

from __future__ import annotations

import pytest

from aleph.world.step import (
    ACCEPTANCE_UNDEFINED,
    AdjointClosurePredicate,
    StepContext,
    UndefinedAcceptance,
    WorldStep,
)


class _Arena:
    def __init__(self) -> None:
        self.checks = 0

    def assert_partitioned(self) -> None:
        self.checks += 1


def test_the_self_check_runs() -> None:
    """The module's own Sanity Gate, executed by the suite as well as by __main__."""
    from aleph.world.step import _demo

    _demo()


def test_an_unjudged_step_is_not_accepted() -> None:
    """The single most important property here.

    `is_accepted` must be False for ACCEPTANCE_UNDEFINED. Treating an unjudged step as accepted is
    exactly the promotion the charter forbids, and it is the one a careless `if verdict:` would make.
    """
    ws = WorldStep(arena=_Arena(), predicate=UndefinedAcceptance())
    v = ws.step(dt_phys_s=0.05, solve=lambda: None)
    assert v.verdict == ACCEPTANCE_UNDEFINED
    assert v.is_accepted is False
    assert v.is_judged is False
    assert ws.as_record()["n_accepted"] == 0


def test_one_unjudged_step_keeps_accepted_step_evidence_null() -> None:
    """`CellStatePosterior.accepted_step_evidence` is nullable for this reason.

    A run where any step went unjudged may not offer itself as admissible evidence, however many of
    its other steps were judged.
    """
    ws = WorldStep(arena=_Arena(), predicate=AdjointClosurePredicate())
    ws.step(dt_phys_s=0.05, solve=lambda: None, assemble_residual=lambda: (1e-9, 1e-8))
    assert ws.as_record()["accepted_step_evidence"] == "SEE_PREDICATE"
    ws.step(dt_phys_s=0.05, solve=lambda: None)          # no residual -> unjudged
    assert ws.as_record()["accepted_step_evidence"] is None


def test_the_verdict_carries_what_judged_it() -> None:
    """A verdict without its judge is how a rounding bound got quoted as convergence for months."""
    ws = WorldStep(arena=_Arena(), predicate=AdjointClosurePredicate())
    v = ws.step(dt_phys_s=0.05, solve=lambda: None, assemble_residual=lambda: (1e-9, 1e-8))
    assert v.verdict == "ACCEPTED"
    assert "NOT convergence" in v.predicate_tests
    row = v.as_record()
    assert row["predicate"] == "adjoint_closure" and "ADJOINT CLOSURE ONLY" in row["predicate_tests"]


def test_a_predicate_must_explain_itself() -> None:
    class _Mute:
        name = "mute"
        tests = "  "

        def decide(self, ctx: StepContext) -> str:
            return "ACCEPTED"

    with pytest.raises(ValueError, match="empty"):
        WorldStep(arena=_Arena(), predicate=_Mute())

    with pytest.raises(TypeError, match="AcceptancePredicate"):
        WorldStep(arena=_Arena(), predicate=object())


def test_a_third_verdict_state_is_refused() -> None:
    """Anything that is neither a decision nor UNDEFINED reads as approval downstream."""
    class _Sneaky:
        name = "sneaky"
        tests = "invents a state"

        def decide(self, ctx: StepContext) -> str:
            return "PROBABLY_FINE"

    with pytest.raises(ValueError, match="third state"):
        WorldStep(arena=_Arena(), predicate=_Sneaky()).step(dt_phys_s=0.05, solve=lambda: None)


def test_the_partition_invariant_runs_every_step() -> None:
    """Promoted from build time: a step that corrupts a range is caught by the step that did it."""
    arena = _Arena()
    ws = WorldStep(arena=arena, predicate=UndefinedAcceptance())
    for _ in range(3):
        ws.step(dt_phys_s=0.05, solve=lambda: None)
    assert arena.checks == 3


def test_a_failing_solve_records_nothing() -> None:
    """A step that did not complete must not leave a verdict behind to be counted later."""
    ws = WorldStep(arena=_Arena(), predicate=UndefinedAcceptance())

    def _boom() -> None:
        raise RuntimeError("inner solve failed")

    with pytest.raises(RuntimeError):
        ws.step(dt_phys_s=0.05, solve=_boom)
    assert ws.verdicts == []


@pytest.mark.parametrize("dt", [0.0, -1.0, float("nan"), float("inf")])
def test_dt_is_checked(dt: float) -> None:
    ws = WorldStep(arena=_Arena(), predicate=UndefinedAcceptance())
    with pytest.raises(ValueError):
        ws.step(dt_phys_s=dt, solve=lambda: None)


def test_the_observables_hook_is_the_seat_e1_needs() -> None:
    """The seat (e) 1 needs — and gamma has now arrived through it.

    The point of the hook is that a stationarity predicate becomes writable WITHOUT editing the step
    loop: the day gamma is emitted, only the predicate passed in changes. ⚠ That day was 2026-08-21
    (`world/observe_gamma.py`, `984c1292`), and this docstring said "cortical tension is not emitted"
    for a day after it was. The hook held its side of the bargain — the step loop did not change.
    """
    class _NeedsGamma:
        name = "needs_gamma"
        tests = "would test stationarity of cortical tension; UNDEFINED while it is absent"

        def decide(self, ctx: StepContext) -> str:
            if "gamma_pn_per_um" not in ctx.observables:
                return ACCEPTANCE_UNDEFINED
            return "ACCEPTED"

    ws = WorldStep(arena=_Arena(), predicate=_NeedsGamma())
    assert ws.step(dt_phys_s=0.05, solve=lambda: None).verdict == ACCEPTANCE_UNDEFINED
    v = ws.step(dt_phys_s=0.05, solve=lambda: None, observables=lambda: {"gamma_pn_per_um": 3.7})
    assert v.verdict == "ACCEPTED"


def test_a_run_of_zero_steps_offers_no_evidence() -> None:
    """`len(judged) == len(verdicts)` is 0 == 0, so an empty run claimed its evidence was available.

    Found on a fresh read rather than by a failure: the driver refuses a zero-step run so it never
    reached an artifact. But `as_record` IS the artifact interface, and "no unjudged steps" is not
    the same claim as "every step was judged" — which is the vacuous pass this module exists to
    refuse, appearing inside the module itself.
    """
    ws = WorldStep(arena=_Arena(), predicate=AdjointClosurePredicate())
    r = ws.as_record()
    assert r["n_steps"] == 0
    assert r["accepted_step_evidence"] is None, "a run that took no step offers no evidence"

    # ...and one judged step is not zero steps: the distinction has to survive the fix.
    ws.step(dt_phys_s=0.05, solve=lambda: None, assemble_residual=lambda: (1e-9, 1e-8))
    assert ws.as_record()["accepted_step_evidence"] == "SEE_PREDICATE"
