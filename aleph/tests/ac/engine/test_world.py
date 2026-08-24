"""Accepted physical-clock transaction gates for the component cell engine."""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from aleph.engine import CellActor, CellWorldTransaction, reference_cell_architecture


@dataclass(slots=True)
class _Participant:
    calls: list[tuple[object, ...]] = field(default_factory=list)

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted: object) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append(("ledger", ledger))


class _PartialParticipant:
    def snapshot_candidate(self) -> None:
        pass


def test_world_uses_one_predicate_and_one_physical_step_for_all_runtimes() -> None:
    actor = CellActor(reference_cell_architecture())
    membrane = _Participant()
    composite_fa = _Participant()
    actor.bind_component("membrane", membrane)
    actor.bind_connector("fa_actin_anchor", composite_fa)
    actor.bind_connector("integrin_collagen_clutch", composite_fa)
    world = CellWorldTransaction(actor)

    accepted_d = object()
    ledger = object()
    world.begin_candidate()
    world.accumulate_ledgers(ledger)
    world.finalize_candidate(accepted_d, dt_phys=0.05, rng_seed=9)

    expected = [
        ("snapshot",),
        ("ledger", ledger),
        ("rollback", accepted_d),
        ("commit", accepted_d, 0.05, 9),
    ]
    assert membrane.calls == expected
    assert composite_fa.calls == expected
    assert world.participants == (membrane, composite_fa)


def test_world_rejects_incomplete_transaction_api() -> None:
    actor = CellActor(reference_cell_architecture())
    actor.bind_component("membrane", _PartialParticipant())
    with pytest.raises(TypeError, match="incomplete transaction API"):
        CellWorldTransaction(actor)


def test_production_world_rejects_missing_component_or_connector_bindings() -> None:
    actor = CellActor(reference_cell_architecture())
    actor.bind_component("membrane", _Participant())
    with pytest.raises(ValueError, match="not fully bound"):
        CellWorldTransaction(actor, require_complete=True)


def test_world_rejects_registry_mutation_after_construction() -> None:
    actor = CellActor(reference_cell_architecture())
    actor.bind_component("membrane", _Participant())
    world = CellWorldTransaction(actor)
    actor.bind_component("cortex", _Participant())
    with pytest.raises(RuntimeError, match="bindings changed"):
        world.begin_candidate()


@pytest.mark.parametrize(
    ("dt_phys", "rng_seed"),
    [(0.0, 0), (-0.1, 0), (float("nan"), 0), (0.1, -1), (0.1, True), (0.1, 1.5)],
)
def test_world_rejects_invalid_physical_step_inputs(dt_phys: float, rng_seed: object) -> None:
    actor = CellActor(reference_cell_architecture())
    actor.bind_component("membrane", _Participant())
    world = CellWorldTransaction(actor)
    with pytest.raises(ValueError):
        world.finalize_candidate(object(), dt_phys=dt_phys, rng_seed=rng_seed)
