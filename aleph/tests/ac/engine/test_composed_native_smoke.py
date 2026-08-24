"""CPU-safe smoke gate: the REAL private-array owners compose into ONE native cell, disjointly.

This is the structural gate for :func:`aleph.engine.composed_native.compose_native_cell_world` — the
real-owner analog of the census-double bridge.  It proves that ``sf_arc`` + ECM + NMII bind into one
``build_composed_cell_world`` triple and one top-level ``CellTransaction``, that the owners are bound
VERBATIM (never replaced by a census double), that the cortex is a PORT (never a participant), that the four
populations are DISJOINT (using the REAL ``sf_arc`` host ledger), and that the interior fluid column is
honestly EXCLUDED.

It is CUDA-free by construction (the numpy-reference / source pattern the project mandates for the dev Mac):
the real Warp-resident owners cannot be built here, so — exactly as the CUDA lane injects the reals — this
gate injects owner doubles that duck-type the real owner API and drives the composition + transaction
invariants that do NOT need a device.  No ``wp.launch`` runs.  ``test_native_composition_smoke.py`` and the
foundation static-contract gate are untouched and stay green.  A separate SOURCE-level test proves the CUDA
lane (:func:`build_native_composed_cell_world`) constructs the REAL owner builders, not doubles.
"""

from __future__ import annotations

import inspect

import pytest

from aleph.engine import composed_native
from aleph.engine.composed_native import (
    ECM_COMPONENT,
    INTERIOR_COLUMN_FLUID_CONNECTORS,
    NMII_COMPONENT,
    NMII_CORTEX_MOTOR,
    SF_COMPONENT,
    compose_native_cell_world,
)
from aleph.engine.population import PopulationLedger, assert_disjoint_populations
from aleph.engine.sf_population import build_sf_arc_population


# --------------------------------------------------------------------------------------------------
# Owner doubles: the CUDA lane injects the real Warp owners; here we duck-type their host API.
# They are NOT census doubles substituted BY the builder — they stand in for the objects the caller
# builds and passes, and the gate proves the builder binds these exact objects verbatim.
# --------------------------------------------------------------------------------------------------
class _OwnerDouble:
    """A real-owner stand-in with the full transaction/ledger/mechanics API + a call log."""

    def __init__(self, *, kind: str, n_minifilaments: int = 0) -> None:
        self.kind = kind
        self.n_minifilaments = n_minifilaments
        self.force_d = object()  # a stand-in for the owned device force array (never launched here)
        self.calls: list[tuple] = []

    # mechanics hooks (owner-specific; only invoked by the pipeline, which the gate does not run)
    def accumulate(self) -> None:
        self.calls.append(("accumulate",))

    def accumulate_internal(self) -> None:
        self.calls.append(("accumulate_internal",))

    def geometry(self):
        return self

    # transaction participant API
    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted, dt_phys, rng_seed) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger) -> None:
        self.calls.append(("ledger", ledger))


class _ConnectorDouble(_OwnerDouble):
    """The MOTOR connector double — a participant that also proposes events (kinetic connector)."""

    def __init__(self) -> None:
        super().__init__(kind="connector")
        self.name = NMII_CORTEX_MOTOR
        self.component_a = NMII_COMPONENT
        self.component_b = "cortex"

    def propose_events(self, rates, dt_phys, rng_seed, neighbors) -> None:
        self.calls.append(("propose_events", dt_phys, rng_seed))

    def accumulate_candidate(self, actuator, port) -> None:
        self.calls.append(("accumulate_candidate",))

    def compute_loads(self, actuator, port) -> None:
        self.calls.append(("compute_loads",))


class _PortDouble:
    """The cortex bind-target PORT — owns a force array but NO transaction API (never a participant)."""

    def __init__(self) -> None:
        self.component = "cortex"
        self.force_d = object()


class _SpyClock:
    """A device-free event clock stand-in (the CUDA lane uses a WarpEventClock)."""

    def __init__(self, base_seed: int = 7) -> None:
        self.base_seed = base_seed
        self.advances: list[tuple] = []

    def advance(self, accepted, dt_phys) -> None:
        self.advances.append((accepted, dt_phys))


def _real_sf_population():
    """The REAL host ``sf_arc`` disjoint population (NumPy only) + its id_base=1e6 ledger."""
    pop = build_sf_arc_population(n_ventral=8, n_dorsal=4, n_arc=4, n_cap=4, n_per_fiber=9)
    pop.assert_partitioned()
    return pop


def _ledgers(pop):
    """cortex [0,70686) · sf_arc [1e6,…) · ecm [2e6,…) · nmii [3e6,…) — every filament in exactly one."""
    cortex = PopulationLedger("cortex", id_base=0, capacity=70686)
    cortex.seed_active(range(0, 128))
    ecm = PopulationLedger(ECM_COMPONENT, id_base=2_000_000, capacity=64)
    ecm.seed_active(range(2_000_000, 2_000_064))
    nmii = PopulationLedger(NMII_COMPONENT, id_base=3_000_000, capacity=32)
    nmii.seed_active(range(3_000_000, 3_000_032))
    return (cortex, pop.ledger, ecm, nmii)


def _compose(pop):
    """Drive the REAL wiring over owner doubles (the CUDA lane injects the real Warp owners)."""
    sf = _OwnerDouble(kind="sf")
    ecm = _OwnerDouble(kind="ecm")
    nmii = _OwnerDouble(kind="nmii", n_minifilaments=32)
    connector = _ConnectorDouble()
    port = _PortDouble()
    clock = _SpyClock()
    cell = compose_native_cell_world(
        sf_owner=sf, ecm_owner=ecm, nmii_actuator=nmii, nmii_cortex_connector=connector,
        cortex_port=port, clock=clock, population_ledgers=_ledgers(pop),
    )
    return cell, {"sf": sf, "ecm": ecm, "nmii": nmii, "connector": connector, "port": port, "clock": clock}


# --------------------------------------------------------------------------------------------------
# 1. the composed world binds the REAL owners verbatim at their reference-architecture slots
# --------------------------------------------------------------------------------------------------
def test_composed_world_binds_real_owners_not_doubles() -> None:
    pop = _real_sf_population()
    cell, parts = _compose(pop)

    actor = cell.world.actor
    # each owner is bound VERBATIM (identity) — the builder never wraps it in a census double.
    assert actor.component_runtime(SF_COMPONENT) is parts["sf"]
    assert actor.component_runtime(ECM_COMPONENT) is parts["ecm"]
    assert actor.component_runtime(NMII_COMPONENT) is parts["nmii"]
    assert actor.connector_runtime(NMII_CORTEX_MOTOR) is parts["connector"]

    registered = set(cell.world.registered_components())
    assert {SF_COMPONENT, ECM_COMPONENT, NMII_COMPONENT} <= registered
    assert NMII_CORTEX_MOTOR in set(cell.world.registered_connectors())


def test_the_internal_arc_edge_is_registered_when_it_is_supplied() -> None:
    """`dorsal_arc_crosslink` ran for months without being REGISTERED, so nothing counted it.

    `build_sf_arc_state_owner` binds the connector by default and `SFArcStateOwner.accumulate` launches
    its `link_spring_kernel`, but the composer was never handed the object — so the actor did not know
    the edge existed and the 2026-08-11 device-run census reported it `NOT_BOUND` while its kernel ran.
    This pins the registration: an edge whose force is in the composed run must be COUNTED as bound.
    """
    from aleph.engine.composed_native import DORSAL_ARC_CROSSLINK

    class _ArcDouble:
        n_joints = 4
        name, component_a, component_b = DORSAL_ARC_CROSSLINK, SF_COMPONENT, SF_COMPONENT

    pop = _real_sf_population()
    arc = _ArcDouble()
    cell = compose_native_cell_world(
        sf_owner=_OwnerDouble(kind="sf"), ecm_owner=_OwnerDouble(kind="ecm"),
        nmii_actuator=_OwnerDouble(kind="nmii", n_minifilaments=32),
        nmii_cortex_connector=_ConnectorDouble(), cortex_port=_PortDouble(), clock=_SpyClock(),
        population_ledgers=_ledgers(pop), sf_internal_arc_connector=arc,
    )
    assert DORSAL_ARC_CROSSLINK in set(cell.world.registered_connectors())
    assert cell.world.actor.connector_runtime(DORSAL_ARC_CROSSLINK) is arc


def test_an_arc_connector_with_no_joints_is_not_registered() -> None:
    """A zero-joint connector launches nothing, so registering it would count a runtime that cannot run."""
    from aleph.engine.composed_native import DORSAL_ARC_CROSSLINK

    class _EmptyArc:
        n_joints = 0
        name, component_a, component_b = DORSAL_ARC_CROSSLINK, SF_COMPONENT, SF_COMPONENT

    pop = _real_sf_population()
    cell = compose_native_cell_world(
        sf_owner=_OwnerDouble(kind="sf"), ecm_owner=_OwnerDouble(kind="ecm"),
        nmii_actuator=_OwnerDouble(kind="nmii", n_minifilaments=32),
        nmii_cortex_connector=_ConnectorDouble(), cortex_port=_PortDouble(), clock=_SpyClock(),
        population_ledgers=_ledgers(pop), sf_internal_arc_connector=_EmptyArc(),
    )
    assert DORSAL_ARC_CROSSLINK not in set(cell.world.registered_connectors())


# --------------------------------------------------------------------------------------------------
# 2. cortex is a PORT — never bound as a component / never a transaction participant
# --------------------------------------------------------------------------------------------------
def test_cortex_is_a_port_never_a_participant() -> None:
    pop = _real_sf_population()
    cell, parts = _compose(pop)

    # cortex is not a bound component...
    assert "cortex" not in set(cell.world.registered_components())
    # ...and the port object is not among the transaction participants.
    assert parts["port"] not in cell.participants
    assert cell.cortex_port is parts["port"]
    # the three real owners + the motor connector ARE the deduplicated participants.
    for key in ("sf", "ecm", "nmii", "connector"):
        assert parts[key] in cell.participants


# --------------------------------------------------------------------------------------------------
# 3. the interior fluid column is honestly EXCLUDED (PI-gated Card-5 driver seam)
# --------------------------------------------------------------------------------------------------
def test_interior_fluid_column_is_excluded() -> None:
    pop = _real_sf_population()
    cell, _ = _compose(pop)
    registered_connectors = set(cell.world.registered_connectors())
    for connector in INTERIOR_COLUMN_FLUID_CONNECTORS:
        assert connector not in registered_connectors, f"{connector!r} must NOT be composed (PI-gated)"
    # the interior compartments are not bound components either.
    for component in ("membrane", "cytosol", "nucleus"):
        assert component not in set(cell.world.registered_components())


# --------------------------------------------------------------------------------------------------
# 4. the four populations are DISJOINT (REAL sf_arc ledger), asserted at build
# --------------------------------------------------------------------------------------------------
def test_populations_are_disjoint_across_all_owners() -> None:
    pop = _real_sf_population()
    cell, _ = _compose(pop)

    names = [led.component for led in cell.population_ledgers]
    assert set(names) == {"cortex", SF_COMPONENT, ECM_COMPONENT, NMII_COMPONENT}
    assert_disjoint_populations(cell.population_ledgers)  # no block overlap, no shared active id
    # the sf_arc ledger in the composed cell is the REAL host population's ledger (id_base 1e6).
    sf_ledger = next(led for led in cell.population_ledgers if led.component == SF_COMPONENT)
    assert sf_ledger is pop.ledger
    assert sf_ledger.block[0] == 1_000_000

    # a colliding ECM block (poaching the sf_arc namespace) is a REJECTED compose, not a silent double-count.
    bad = _ledgers(pop)
    collide = PopulationLedger(ECM_COMPONENT, id_base=pop.ledger.block[0], capacity=8)
    with pytest.raises(ValueError, match="overlapping global"):
        assert_disjoint_populations((bad[0], pop.ledger, collide, bad[3]))


# --------------------------------------------------------------------------------------------------
# 5. ONE CellTransaction wraps them: one predicate → every rollback+commit; clock advances; disjoint re-run
# --------------------------------------------------------------------------------------------------
def test_one_transaction_accepted_step_reaches_every_owner() -> None:
    pop = _real_sf_population()
    cell, parts = _compose(pop)

    accepted = object()
    solved = []
    cell.transaction.step(dt_phys=0.05, solve=lambda: solved.append(True), accepted_d=accepted)

    assert solved == [True]  # the caller-owned inner solve ran exactly once
    seed = parts["clock"].base_seed  # rng_seed threaded through commit == the clock's base_seed
    # every real owner saw snapshot, then the SAME predicate for rollback + commit.
    for key in ("sf", "ecm", "nmii", "connector"):
        calls = parts[key].calls
        assert ("snapshot",) in calls
        assert ("rollback", accepted) in calls
        assert ("commit", accepted, 0.05, seed) in calls
    # the connector proposed events; the clock advanced once on the same predicate.
    assert any(c[0] == "propose_events" for c in parts["connector"].calls)
    assert parts["clock"].advances == [(accepted, 0.05)]
    # the participants are deduplicated (no owner appears twice).
    assert len(cell.participants) == len(set(id(p) for p in cell.participants))


# --------------------------------------------------------------------------------------------------
# 6. SOURCE proof: the CUDA lane constructs the REAL owner builders, not census doubles
# --------------------------------------------------------------------------------------------------
def test_cuda_lane_constructs_real_owner_builders() -> None:
    src = inspect.getsource(composed_native.build_native_composed_cell_world)
    # the native builder must call the REAL owner builders (not _CensusRuntime / _FacadeSpy).
    assert "build_sf_arc_state_owner(" in src
    assert "build_ecm_state_owner(" in src
    assert "make_event_clock(" in src
    assert "compose_native_cell_world(" in src
    for forbidden in ("_CensusRuntime", "_FacadeSpy"):
        assert forbidden not in src

    # build_sf_arc_state_owner itself wires the real SF kernel-bound mechanics.
    sf_src = inspect.getsource(composed_native.build_sf_arc_state_owner)
    assert "build_sf_filament_mechanics(" in sf_src
    assert "build_sf_mechanics_topology(" in sf_src
