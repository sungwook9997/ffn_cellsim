"""Classification + hook-swap guards for the per-connector device-run census.

The census's GPU-side observations are checked by its own positive/negative controls at run time.
What CANNOT be checked there is the classification table itself and the hook swap/restore — a bug in
either would make every row wrong in the same direction, which is exactly the failure a control does
not catch. Both are CPU-testable, so they are tested here.
"""

from __future__ import annotations

from dataclasses import dataclass

from aleph.scripts.ac_connector_devicerun_census import _Census, _force_channel, _verdict


def _rec(calls=(), launches=(), l1=0.0, ledger=0.0) -> dict:
    return {"calls": list(calls), "launches": list(launches), "l1_delta_pN": l1,
            "l1_before_pN": 0.0, "l1_after_pN": 0.0,
            "ledger_pN": ledger, "instrumented_hooks": [], "errors": []}


def test_verdict_never_folds_not_called_into_stub() -> None:
    """A connector the pipeline never reached is NOT a stub finding — the two must stay distinct."""
    assert _verdict(_rec()) == "NOT_CALLED"
    assert _verdict(_rec(calls=["accumulate"])) == "STUB"


def test_the_verdict_is_the_launch_observation_and_nothing_else() -> None:
    """A kernel that ran and moved no force still RAN. Folding force in is how 2.7e-21 pN read as force."""
    assert _verdict(_rec(calls=["accumulate"], launches=["k"])) == "DEVICE_RUN"
    assert _verdict(_rec(calls=["accumulate"], launches=["k"], l1=1e-9)) == "DEVICE_RUN"


def test_verdict_requires_a_launch_for_device_run() -> None:
    """Force moving with no kernel launched is not a device-run claim — it is somebody else's force."""
    assert _verdict(_rec(calls=["accumulate"], l1=42.0)) == "STUB"


def test_force_channel_states_the_fact_and_declines_the_judgement() -> None:
    """It must not grade 'physical' vs 'noise' — that needs a chosen scale the charter forbids in a gate."""
    assert _force_channel(_rec(), 551434) == "ZERO"
    assert _force_channel(_rec(l1=2.7e-21), 551434) == "NONZERO"
    assert _force_channel(_rec(l1=42.0), 551434) == "NONZERO"


@dataclass(frozen=True, slots=True)
class _Frozen:
    """The shape 7 of the 13 real connector runtime classes have: frozen + slots, refuses ``setattr``."""

    tag: int = 0

    def accumulate(self) -> None: ...
    def accumulate_ledger(self, _l) -> None: ...


def test_instrumentation_works_on_a_frozen_slots_runtime_and_restores_the_class() -> None:
    """Instance ``setattr`` is impossible on these; class-level wrapping is why the census works at all."""
    original = _Frozen.accumulate
    census = _Census(probe=None)
    census.instrument("c", _Frozen())
    assert _Frozen.accumulate is not original, "class hook was not wrapped"
    assert set(census.records["c"]["instrumented_hooks"]) == {"accumulate", "accumulate_ledger"}
    census.restore()
    assert _Frozen.accumulate is original, "restore did not put the class back as found"


def test_two_edges_sharing_one_runtime_object_share_one_record() -> None:
    """``ImmersedPorousTransfer`` serves 7 declared edges. One object cannot report per-edge, and the
    census must not pretend otherwise — both names must resolve to the SAME record."""
    census = _Census(probe=None)
    runtime = _Frozen()
    census.instrument("edge_a", runtime)
    census.instrument("edge_b", runtime)
    try:
        assert census.records["edge_a"] is census.records["edge_b"]
        assert len(census.by_object) == 1
    finally:
        census.restore()


def test_distinct_objects_of_one_class_are_recorded_separately() -> None:
    """The proxy dispatches on ``id(self)``, so one wrapped class still attributes per object."""
    census = _Census(probe=None)
    census.instrument("a", _Frozen(tag=1))
    census.instrument("b", _Frozen(tag=2))
    try:
        assert census.records["a"] is not census.records["b"]
        assert len(census.by_object) == 2
    finally:
        census.restore()


def test_negative_control_nulls_only_mechanics_hooks() -> None:
    """The nulled connector keeps its ledger hook, or its STUB verdict would be a harness artifact."""
    census = _Census(probe=None)
    nulled = census.nullify(_Frozen())
    assert nulled == ("accumulate",)
    assert "accumulate_ledger" not in nulled


def test_a_runtime_that_refuses_instrumentation_is_recorded_not_silently_skipped() -> None:
    """If wrapping ever fails, the row must say so — never read as NOT_CALLED, which is a code finding."""
    rec = _rec()
    rec["errors"].append("accumulate: not instrumentable on X")
    assert _verdict(rec) == "NOT_INSTRUMENTABLE"
