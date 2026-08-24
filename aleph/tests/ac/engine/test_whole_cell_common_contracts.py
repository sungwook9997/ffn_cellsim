"""Structural gates for the whole-cell common-contracts spine (Step 1, pre-biology).

CUDA-free gates for the four Step-1 deliverables of ``WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md``:

* ``propose_events`` as the enforced 5th scheduler method (component ``has_events`` + connector ``kinetics``),
* :class:`PopulationLedger` free-list discipline + the disjoint-ID no-double-count assert,
* :class:`CellState` schema, sourced-only provenance, and the ``CellState -/-> force`` hard boundary,
* :class:`CellTransaction.step` whole-cell orchestration order under one device predicate.

They inject spies (the CUDA lane injects real Warp-resident runtimes); no biology, no device state.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from aleph.engine import (
    CellActor,
    CellState,
    CellTransaction,
    CellWorldTransaction,
    Distribution,
    PopulationLedger,
    Rate,
    assert_cellstate_writes_no_force,
    assert_disjoint_populations,
    mcf7_reference_state,
    reference_cell_architecture,
)
from aleph.engine.contracts import (
    CellArchitecture,
    ComponentContract,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
    ConnectorScope,
)


# --------------------------------------------------------------------------------------------------
# Shared spies
# --------------------------------------------------------------------------------------------------
@dataclass(slots=True)
class _FullSpy:
    """A state-owner / connector double with the full transaction + ledger + mechanics + events API."""

    name: str | None = None
    component_a: str | None = None
    component_b: str | None = None
    mechanical_group: str | None = None
    with_events: bool = True
    calls: list[tuple] = field(default_factory=list)

    def accumulate(self, *args) -> None:
        self.calls.append(("accumulate", args))

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted, dt_phys, rng_seed) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger) -> None:
        self.calls.append(("ledger", ledger))

    def __getattr__(self, item):
        # propose_events exists only when with_events; keeps a "no events" spy genuinely event-less.
        if item == "propose_events":
            try:
                enabled = object.__getattribute__(self, "with_events")
            except AttributeError:
                raise AttributeError(item) from None
            if enabled:
                def _propose(rates, dt_phys, rng_seed, neighbors):
                    self.calls.append(("propose_events", dt_phys, rng_seed))
                return _propose
        raise AttributeError(item)


@dataclass(slots=True)
class _SpyClock:
    base_seed: int = 11
    calls: list[tuple] = field(default_factory=list)

    def advance(self, accepted, dt_phys) -> None:
        self.calls.append(("advance", accepted, dt_phys))


# --------------------------------------------------------------------------------------------------
# 1. propose_events enforcement
# --------------------------------------------------------------------------------------------------
def _mini_architecture(*, component_has_events: bool, connector_kinetic: bool) -> CellArchitecture:
    return CellArchitecture(
        components=(
            ComponentContract(
                "alpha", ComponentRole.ACTIVE_LOAD_PATH, "rep", "solve",
                has_events=component_has_events,
            ),
            ComponentContract("beta", ComponentRole.ACTIVE_LOAD_PATH, "rep", "solve"),
        ),
        connectors=(
            ConnectorContract(
                "alpha_beta_link", ConnectorFamily.TRANSIENT_ACTIN, "alpha", "beta",
                connector_kinetic, connector_kinetic,
            ),
        ),
    )


def _bind_all(architecture: CellArchitecture, *, events_everywhere: bool) -> CellActor:
    actor = CellActor(architecture)
    for component in architecture.components:
        actor.bind_component(component.name, _FullSpy(with_events=events_everywhere))
    for connector in architecture.connectors:
        actor.bind_connector(
            connector.name,
            _FullSpy(
                name=connector.name,
                component_a=connector.component_a,
                component_b=connector.component_b,
                with_events=events_everywhere,
            ),
        )
    return actor


def test_kinetic_connector_must_expose_propose_events() -> None:
    architecture = _mini_architecture(component_has_events=False, connector_kinetic=True)
    actor = _bind_all(architecture, events_everywhere=False)
    with pytest.raises(TypeError, match="propose_events"):
        actor.assert_fully_bound()
    ok = _bind_all(architecture, events_everywhere=True)
    ok.assert_fully_bound()  # does not raise


def test_has_events_component_must_expose_propose_events() -> None:
    architecture = _mini_architecture(component_has_events=True, connector_kinetic=False)
    actor = _bind_all(architecture, events_everywhere=False)
    with pytest.raises(TypeError, match="propose_events"):
        actor.assert_fully_bound()
    ok = _bind_all(architecture, events_everywhere=True)
    ok.assert_fully_bound()


def test_non_kinetic_connector_does_not_require_propose_events() -> None:
    architecture = _mini_architecture(component_has_events=False, connector_kinetic=False)
    actor = _bind_all(architecture, events_everywhere=False)
    actor.assert_fully_bound()  # no kinetics, no has_events -> propose_events not required


def test_reference_architecture_leaves_component_events_for_biology() -> None:
    # The Step-1 contract flips no component's has_events; that is a per-component biology decision.
    for component in reference_cell_architecture().components:
        assert component.has_events is False


# --------------------------------------------------------------------------------------------------
# 2. PopulationLedger — free-list discipline + disjoint-ID assert
# --------------------------------------------------------------------------------------------------
def test_population_allocate_release_conserves_fixed_capacity() -> None:
    pop = PopulationLedger("cortex", id_base=0, capacity=4)
    assert pop.dormant_count == 4 and pop.active_count == 0
    a = pop.allocate()
    b = pop.allocate()
    assert {a, b} <= set(range(0, 4))
    assert pop.active_count == 2 and pop.dormant_count == 2
    pop.release(a)
    assert pop.active_count == 1 and pop.dormant_count == 3
    pop.assert_invariants()  # active + dormant == capacity throughout


def test_population_sever_pops_one_more_slot() -> None:
    pop = PopulationLedger("cortex", id_base=10, capacity=3)
    parent = pop.allocate()
    child = pop.sever(parent)
    assert child != parent and pop.owns(child)
    assert pop.active_count == 2
    pop.assert_invariants()


def test_population_never_resamples_past_capacity() -> None:
    pop = PopulationLedger("filopodium", id_base=0, capacity=1)
    pop.allocate()
    with pytest.raises(RuntimeError, match="free-list exhausted"):
        pop.allocate()  # no growing to a target density


def test_population_seed_active_must_be_in_block() -> None:
    pop = PopulationLedger("sf_arc", id_base=100, capacity=5)
    pop.seed_active([100, 101])
    assert pop.active_count == 2
    with pytest.raises(ValueError, match="outside its block"):
        pop.seed_active([200])


def test_disjoint_populations_reject_overlapping_blocks() -> None:
    cortex = PopulationLedger("cortex", id_base=0, capacity=100)
    overlap = PopulationLedger("sf_arc", id_base=50, capacity=100)
    with pytest.raises(ValueError, match="overlapping global"):
        assert_disjoint_populations([cortex, overlap])


def test_disjoint_populations_accept_separate_blocks() -> None:
    cortex = PopulationLedger("cortex", id_base=0, capacity=100)
    sf = PopulationLedger("sf_arc", id_base=100, capacity=50)
    cortex.seed_active([0, 1, 2])
    sf.seed_active([100, 149])
    assert_disjoint_populations([cortex, sf])  # SF is a separate, disjoint population


def test_disjoint_populations_reject_duplicate_component() -> None:
    a = PopulationLedger("cortex", id_base=0, capacity=10)
    b = PopulationLedger("cortex", id_base=10, capacity=10)
    with pytest.raises(ValueError, match="distinct component"):
        assert_disjoint_populations([a, b])


# --------------------------------------------------------------------------------------------------
# 3. CellState — schema, sourced-only, force boundary
# --------------------------------------------------------------------------------------------------
def test_mcf7_reference_axes_are_the_ratified_list() -> None:
    state = mcf7_reference_state()
    assert state.axes == (
        "MCF7", "hybrid_E_M", "G1", "collagen_spread", "polarized", "physiological_resting",
    )


def test_cellstate_rate_and_inventory_carry_provenance() -> None:
    state = CellState(
        lineage="MCF7", emt="hybrid_E_M", cycle="G1", adhesion="collagen_spread",
        geometry="polarized", osmotic_state="physiological_resting",
        rates={("cortex", "nucleation"): Rate("nucleation", 0.5, "KB-3.18")},
        inventory={"cortex": Distribution("length", {"mean": 3.0}, "KB-3.18")},
    )
    assert state.rate("cortex", "nucleation").kb_source == "KB-3.18"
    assert state.initial_inventory("cortex").kind == "length"


def test_cellstate_unsourced_rate_raises_not_defaults() -> None:
    state = mcf7_reference_state()
    with pytest.raises(KeyError, match="no sourced rate"):
        state.rate("cortex", "nucleation")
    with pytest.raises(KeyError, match="no sourced initial inventory"):
        state.initial_inventory("cortex")


def test_rate_without_kb_source_is_rejected() -> None:
    with pytest.raises(ValueError, match="provenance"):
        Rate("nucleation", 0.5, "")


def test_cellstate_cannot_emit_force_or_count() -> None:
    with pytest.raises(ValueError, match="EMERGE"):
        Rate("force", 1.0, "KB-x")
    with pytest.raises(ValueError, match="EMERGE"):
        Distribution("tension", {"v": 1.0}, "KB-x")
    with pytest.raises(ValueError, match="EMERGE"):
        Distribution("filament_count", {"n": 70686.0}, "KB-x")


def test_cellstate_writes_no_force_boundary_holds() -> None:
    assert_cellstate_writes_no_force(mcf7_reference_state())


def test_cellstate_with_force_writer_is_rejected() -> None:
    class _Leaky(CellState):
        __slots__ = ()

        def set_tension(self, value):  # a forbidden mechanical writer
            return value

    leaky = _Leaky(
        lineage="MCF7", emt="hybrid_E_M", cycle="G1", adhesion="collagen_spread",
        geometry="polarized", osmotic_state="physiological_resting",
    )
    with pytest.raises(AssertionError, match="forbidden mechanical/topological writer"):
        assert_cellstate_writes_no_force(leaky)


# --------------------------------------------------------------------------------------------------
# 4. CellTransaction.step — whole-cell orchestration order
# --------------------------------------------------------------------------------------------------
def _two_participant_world() -> tuple[CellWorldTransaction, _FullSpy, _FullSpy]:
    actor = CellActor(reference_cell_architecture())
    cortex = _FullSpy(with_events=True)
    fa = _FullSpy(with_events=True)
    actor.bind_component("cortex", cortex)
    actor.bind_connector("fa_actin_anchor", fa)
    actor.bind_connector("integrin_collagen_clutch", fa)  # composite: one participant
    return CellWorldTransaction(actor), cortex, fa


def test_step_runs_snapshot_propose_solve_finalize_advance_in_order() -> None:
    world, cortex, fa = _two_participant_world()
    clock = _SpyClock(base_seed=11)
    txn = CellTransaction(world, clock)

    order: list[str] = []
    accepted = object()

    def solve() -> None:
        order.append("solve")

    # Wrap participant hooks to record global order via the clock/solve interleave.
    txn.step(dt_phys=0.05, solve=solve, accepted_d=accepted)

    # composite FA joint is one proposer, cortex the other -> exactly two proposers.
    assert len(txn.proposers) == 2
    # every participant saw snapshot -> propose_events -> ... -> rollback -> commit with the base seed.
    for spy in (cortex, fa):
        kinds = [c[0] for c in spy.calls]
        assert kinds == ["snapshot", "propose_events", "rollback", "commit"]
        assert spy.calls[1] == ("propose_events", 0.05, 11)
        assert spy.calls[-1] == ("commit", accepted, 0.05, 11)
    # solve ran between propose and finalize; clock advanced once with the SAME predicate.
    assert order == ["solve"]
    assert clock.calls == [("advance", accepted, 0.05)]


def test_step_requires_a_predicate_or_a_ledger_gate() -> None:
    world, _, _ = _two_participant_world()
    txn = CellTransaction(world, _SpyClock())
    with pytest.raises(ValueError, match="no acceptance predicate"):
        txn.step(dt_phys=0.05, solve=lambda: None)


def test_step_rejects_invalid_dt() -> None:
    world, _, _ = _two_participant_world()
    txn = CellTransaction(world, _SpyClock())
    with pytest.raises(ValueError, match="dt_phys"):
        txn.step(dt_phys=0.0, solve=lambda: None, accepted_d=object())


def test_transaction_asserts_disjoint_populations_at_build_and_after_step() -> None:
    world, _, _ = _two_participant_world()
    cortex = PopulationLedger("cortex", id_base=0, capacity=10)
    sf = PopulationLedger("sf_arc", id_base=10, capacity=10)
    cortex.seed_active([0, 1])
    sf.seed_active([10, 11])
    txn = CellTransaction(world, _SpyClock(), population_ledgers=(cortex, sf))
    txn.step(dt_phys=0.05, solve=lambda: None, accepted_d=object())  # disjoint assert passes each step

    # overlapping populations are a rejected build (blocks would double-count a filament).
    bad = PopulationLedger("sf_arc", id_base=0, capacity=10)
    bad.seed_active([0])
    with pytest.raises(ValueError, match="double-counted"):
        CellTransaction(world, _SpyClock(), population_ledgers=(cortex, bad))


def test_clock_must_expose_advance_and_base_seed() -> None:
    world, _, _ = _two_participant_world()
    with pytest.raises(TypeError, match="advance"):
        CellTransaction(world, object())
