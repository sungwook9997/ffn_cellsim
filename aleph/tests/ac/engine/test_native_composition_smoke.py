"""CPU-safe smoke gate: the session's advanced components compose into ONE cell world, disjointly.

This session landed four components as separate native slices — cortex (GATE-A/B), ``sf_arc``
(KERNEL_BOUND over a disjoint population), native ECM (device-SoA topology + constitutive force),
and the interior-column ``CellTransaction`` scaffold (membrane/cortex/cytosol/nucleus).  Each has its
OWN per-slice native gate (``scripts/ac_gate_b_*``); none of them is yet assembled into the single
composed-world path (``build_composed_cell_world`` + a whole-cell ``CellTransaction``).  See
``docs/v2_audit/cell_engine/NATIVE_COMPOSITION_SCOPE_2026-07-25.md`` for the full gap analysis.

This gate proves the piece that IS ready now: the declarative ``reference_cell_architecture`` composed
through ``build_composed_cell_world`` binds cortex + sf_arc + ecm + the interior-column compartments
into ONE ``(actor, pipeline, transaction)`` triple, that the exact-once force-assembly pipeline
references each of them, and that their populations are DISJOINT — using the REAL ``sf_arc`` population
ledger (host NumPy, no device), not a placeholder count.

It is CUDA-free by construction (the numpy-reference / source pattern the project mandates for the dev
Mac): the CUDA lane injects the real Warp-resident owners; here we inject census/transaction doubles for
the CUDA-only owners (ECM SoA, the aliased surface owners, the Biot field) and drive the composition +
population invariants that do NOT need a device.  No ``wp.launch`` runs.  The foundation static-contract
gate (``test_whole_cell_common_contracts.py``) is untouched and stays green.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from aleph.engine import (
    ComposedCellRuntimes,
    FacadeDispatch,
    PopulationLedger,
    assert_disjoint_populations,
    build_composed_cell_world,
    dump_composed_world_state,
    reference_cell_architecture,
)
from aleph.engine.dispatch import canonical_facade_claims
from aleph.engine.ecm_world import ECMWorld
from aleph.engine.fluid_core import FluidCore
from aleph.engine.intermediate_filament_rig import IntermediateFilamentRig
from aleph.engine.microtubule_rig import MicrotubuleRig
from aleph.engine.nmii_actuator import NMIIActuator
from aleph.engine.protrusion import ProtrusionActors
from aleph.engine.sf_population import build_sf_arc_population
from aleph.engine.stress_fiber import StressFiberActor
from aleph.engine.surface_body import SurfaceBody

# The four components this session advanced (their canonical component names in the reference graph).
_ADVANCED_COMPONENTS = ("cortex", "sf_arc", "ecm")
# The interior-column vertical slice couples these four compartments through three fluid connectors.
_INTERIOR_COLUMN_COMPONENTS = ("membrane", "cortex", "cytosol", "nucleus")
_INTERIOR_COLUMN_CONNECTORS = (
    "surface_porous_transfer",       # cortex <-> cytosol (Peskin pressure traction)
    "membrane_cytosol_boundary",     # membrane <-> cytosol (Kedem-Katchalsky semipermeable flux)
    "nucleus_cytosol_boundary",      # nucleus  <-> cytosol (impermeable no-flux envelope)
)

_FACADE_METHOD = {
    SurfaceBody: "accumulate_mechanics",
    FluidCore: "candidate_iteration",
    StressFiberActor: "accumulate_mechanics",
    MicrotubuleRig: "accumulate_mechanics",
    IntermediateFilamentRig: "accumulate_mechanics",
    ProtrusionActors: "accumulate_mechanics",
    NMIIActuator: "accumulate_candidate",
    ECMWorld: "accumulate_mechanics",
}

_FA_GROUP = "alpha2beta1_collagen_series"
_FA_EDGES = ("fa_actin_anchor", "integrin_collagen_clutch")


@dataclass(slots=True)
class _FacadeSpy:
    """Dispatch-facade double exposing every canonical orchestration method the manifest may claim."""

    calls: list[tuple] = field(default_factory=list)

    def accumulate_mechanics(self, *args, **kwargs) -> None:
        self.calls.append(("accumulate_mechanics", args, kwargs))

    def accumulate_candidate(self, *args, **kwargs) -> None:
        self.calls.append(("accumulate_candidate", args, kwargs))

    def candidate_iteration(self, *args, **kwargs) -> None:
        self.calls.append(("candidate_iteration", args, kwargs))


class _CensusRuntime:
    """A state-owner / connector double carrying a census + the full transaction/ledger/mechanics API.

    The dev Mac has no CUDA, so the real Warp-resident owners (ECM device SoA, the surface owners that
    alias ``cell.pos_d``/``cell.f_d``, the Biot field) cannot be built here.  This double exposes the
    host-side count attributes ``dump_composed_world_state`` reads plus the snapshot/rollback/commit +
    ledger hooks the composed-world transaction validates.  It owns no geometry and launches no kernel.
    """

    def __init__(self, *, name=None, component_a=None, component_b=None, mechanical_group=None,
                 n_filaments=0, n_nodes=0, n_faces=0, n_endpoints=0) -> None:
        if name is not None:
            self.name = name
        if component_a is not None:
            self.component_a = component_a
        if component_b is not None:
            self.component_b = component_b
        if mechanical_group is not None:
            self.mechanical_group = mechanical_group
        self.n_filaments = n_filaments
        self.n_nodes = n_nodes
        self.n_faces = n_faces
        self.n_endpoints = n_endpoints
        self.calls: list[tuple] = []

    def accumulate(self, *args) -> None:
        self.calls.append(("accumulate", args))

    def propose_events(self, rates, dt_phys, rng_seed, neighbors) -> None:
        self.calls.append(("propose_events", dt_phys, rng_seed))

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted, dt_phys, rng_seed) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger) -> None:
        self.calls.append(("ledger", ledger))


def _real_sf_population():
    """Build the REAL host ``sf_arc`` disjoint population (NumPy only, no device) + its ledger."""
    pop = build_sf_arc_population(n_ventral=8, n_dorsal=4, n_arc=4, n_cap=4, n_per_fiber=9)
    pop.assert_partitioned()
    return pop


def _census_for(pop) -> dict[str, dict[str, int]]:
    """Per-component census; sf_arc's filament count is the REAL population's active count."""
    return {
        "cortex": {"n_filaments": 70686, "n_nodes": 494802, "n_faces": 12},
        "membrane": {"n_nodes": 40962, "n_faces": 81920},
        "cytosol": {"n_nodes": 4096},
        "nucleus": {"n_nodes": 2562, "n_faces": 5120},
        "sf_arc": {"n_filaments": pop.ledger.active_count, "n_nodes": int(pop.pos.shape[0])},
        "nmii": {"n_filaments": 2000, "n_nodes": 30000, "n_endpoints": 8000},
        "microtubule": {"n_filaments": 300, "n_nodes": 6000},
        "intermediate_filament": {"n_filaments": 1500, "n_nodes": 22000},
        "lamellipodium": {"n_filaments": 5000, "n_nodes": 40000},
        "filopodium": {"n_filaments": 60, "n_nodes": 1200},
        "ecm": {"n_filaments": 300, "n_nodes": 6000},
    }


def _runtimes(pop) -> ComposedCellRuntimes:
    """Drive the REAL composed-world builder over census doubles (CUDA lane injects the real owners)."""
    architecture = reference_cell_architecture()
    census = _census_for(pop)

    component_owners: dict[str, object] = {}
    for component in architecture.components:
        if not component.dynamically_evolving and not component.owns_geometry:
            component_owners[component.name] = object()  # world_boundary: inert reference frame
            continue
        component_owners[component.name] = _CensusRuntime(**census.get(component.name, {}))

    # ONE runtime per declared mechanical group, DERIVED from the architecture — a composite series is
    # one object however many series the composition declares.
    group_runtimes: dict[str, object] = {
        group: _CensusRuntime(mechanical_group=group, n_endpoints=256)
        for group in sorted({c.mechanical_group for c in architecture.connectors if c.mechanical_group})
    }
    connector_runtimes: dict[str, object] = {}
    for connector in architecture.connectors:
        if connector.mechanical_group is not None:
            connector_runtimes[connector.name] = group_runtimes[connector.mechanical_group]
            continue
        connector_runtimes[connector.name] = _CensusRuntime(
            name=connector.name, component_a=connector.component_a, component_b=connector.component_b,
        )

    facade_dispatch = {}
    for claim in canonical_facade_claims():
        method = _FACADE_METHOD[claim.owner_type]
        args = (object(), 0.05) if method == "candidate_iteration" else (object(),)
        facade_dispatch[claim.owner_type] = FacadeDispatch(_FacadeSpy(), args=args)

    return ComposedCellRuntimes(component_owners, connector_runtimes, facade_dispatch)


# --------------------------------------------------------------------------------------------------
# 1. one composed world binds every advanced component + the interior-column compartments
# --------------------------------------------------------------------------------------------------
def test_advanced_components_and_interior_column_bind_into_one_world() -> None:
    pop = _real_sf_population()
    world = build_composed_cell_world(_runtimes(pop), require_complete=True)
    registered = set(world.registered_components())

    for name in _ADVANCED_COMPONENTS:
        assert name in registered, f"advanced component {name!r} not bound into the composed world"
    for name in _INTERIOR_COLUMN_COMPONENTS:
        assert name in registered, f"interior-column compartment {name!r} not bound"
    for name in _INTERIOR_COLUMN_CONNECTORS:
        assert name in set(world.registered_connectors()), f"interior-column connector {name!r} not bound"

    # require_complete=True proves the whole graph (all 13 components / 32 connectors) is bindable at once.
    architecture = reference_cell_architecture()
    assert registered == {c.name for c in architecture.components}
    assert set(world.registered_connectors()) == {c.name for c in architecture.connectors}


# --------------------------------------------------------------------------------------------------
# 2. the exact-once force-assembly pipeline references each advanced component + interior connector
# --------------------------------------------------------------------------------------------------
def test_force_assembly_pipeline_references_each_advanced_component() -> None:
    pop = _real_sf_population()
    world = build_composed_cell_world(_runtimes(pop))
    facade_task_names = {task.name for task in world.facade_tasks}

    # each advanced component's mechanics is claimed by exactly one canonical facade task.
    expected_owner = {"cortex": "SurfaceBody", "sf_arc": "StressFiberActor", "ecm": "ECMWorld"}
    for component, facade in expected_owner.items():
        owner = world.pipeline.component_owner(component)
        assert owner == facade, f"{component!r} mechanics owned by {owner!r}, expected {facade!r}"
        assert owner in facade_task_names

    # the interior-column compartments + their three fluid connectors are each dispatched exactly once.
    for component in _INTERIOR_COLUMN_COMPONENTS:
        assert world.pipeline.component_owner(component) in facade_task_names
    for connector in _INTERIOR_COLUMN_CONNECTORS:
        assert world.pipeline.connector_owner(connector) in facade_task_names

    # the sf_arc + ecm connector edges are dispatched too (the load path each advanced component couples on).
    for connector in ("sf_cortex_transient", "dorsal_arc_crosslink", "ecm_crosslink", "fa_actin_anchor"):
        assert world.pipeline.connector_owner(connector) in facade_task_names


# --------------------------------------------------------------------------------------------------
# 3. the advanced components own DISJOINT populations (REAL sf_arc ledger, no shared-node weld)
# --------------------------------------------------------------------------------------------------
def test_advanced_component_populations_are_disjoint() -> None:
    pop = _real_sf_population()

    cortex = PopulationLedger("cortex", id_base=0, capacity=70686)
    cortex.seed_active(range(0, 128))                       # a representative active slice of the cortex block
    sf_arc = pop.ledger                                     # the REAL sf_arc ledger (id_base = 1_000_000)
    ecm = PopulationLedger("ecm", id_base=2_000_000, capacity=512)
    ecm.seed_active(range(2_000_000, 2_000_040))

    # cortex [0, 70686) · sf_arc [1e6, 1e6+40) · ecm [2e6, 2e6+512) — every filament in exactly one component.
    assert cortex.block[1] <= sf_arc.block[0] <= sf_arc.block[1] <= ecm.block[0]
    assert_disjoint_populations([cortex, sf_arc, ecm])     # no block overlap, no shared active id

    # a colliding ECM block (poaching the sf_arc namespace) is a REJECTED build, not a silent double-count.
    collide = PopulationLedger("ecm", id_base=sf_arc.block[0], capacity=8)
    with pytest.raises(ValueError, match="overlapping global"):
        assert_disjoint_populations([cortex, sf_arc, collide])


# --------------------------------------------------------------------------------------------------
# 4. the composed-world dump gives each advanced component a disjoint global filament-ID block
# --------------------------------------------------------------------------------------------------
def test_composed_world_dump_disjoint_blocks_for_advanced_components() -> None:
    pop = _real_sf_population()
    world = build_composed_cell_world(_runtimes(pop))
    dump = dump_composed_world_state(world)

    assert dump.actor_ids_are_unique()
    assert dump.filament_id_blocks_are_disjoint()

    cortex, sf_arc, ecm = dump.by_name("cortex"), dump.by_name("sf_arc"), dump.by_name("ecm")
    # three distinct actors — no shared node, no weld.
    assert len({cortex.actor_id, sf_arc.actor_id, ecm.actor_id}) == 3
    # the sf_arc census count is the REAL host population's active fiber count.
    assert sf_arc.n_filaments == pop.ledger.active_count
    # pairwise non-overlapping global filament blocks.
    for a, b in ((cortex, sf_arc), (cortex, ecm), (sf_arc, ecm)):
        (lo_a, hi_a), (lo_b, hi_b) = a.filament_id_range, b.filament_id_range
        assert hi_a <= lo_b or hi_b <= lo_a


# --------------------------------------------------------------------------------------------------
# 5. all advanced + interior-column owners are participants under ONE accepted-step predicate
# --------------------------------------------------------------------------------------------------
def test_advanced_owners_share_one_accepted_step_predicate() -> None:
    pop = _real_sf_population()
    world = build_composed_cell_world(_runtimes(pop))
    accepted, ledger = object(), object()

    # one snapshot / one ledger accumulate / one predicate to every rollback + commit.
    world.transaction.begin_candidate()
    world.transaction.accumulate_ledgers(ledger)
    world.transaction.finalize_candidate(accepted, dt_phys=0.05, rng_seed=7)

    expected = [("snapshot",), ("ledger", ledger), ("rollback", accepted), ("commit", accepted, 0.05, 7)]
    participants = world.transaction.participants
    # the advanced-component owners are among the deduplicated participants and all saw the same predicate.
    advanced_owners = [
        world.actor.component_runtime(name) for name in (*_ADVANCED_COMPONENTS, *_INTERIOR_COLUMN_COMPONENTS)
    ]
    for owner in advanced_owners:
        assert owner in participants
        assert owner.calls == expected
