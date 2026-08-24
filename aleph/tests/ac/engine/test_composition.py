"""Composed-world wiring gates: concrete-runtime registration + exact-once dispatch + dump ID namespace.

These are CUDA-free structural gates.  They inject runtime doubles (the CUDA lane injects the real
Warp-resident objects) and prove the composed-world builder binds every declared component/connector into one
discoverable ``(actor, pipeline, transaction)`` triple, that the α2β1–collagen FA series joint is one actor and
not two, and that the composed-world dump assigns a disjoint global actor/filament ID namespace so a
per-compartment renderer never double-draws.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from aleph.engine import (
    ComposedCellRuntimes,
    FacadeDispatch,
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
from aleph.engine.stress_fiber import StressFiberActor
from aleph.engine.surface_body import SurfaceBody

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


class _FullRuntime:
    """State-owner / connector double with the full transaction + ledger + mechanics API."""

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


_FA_GROUP = "alpha2beta1_collagen_series"
_FA_EDGES = ("fa_actin_anchor", "integrin_collagen_clutch")

# Per-component census used to exercise the dump ID namespace (values are illustrative, not physiological).
_COMPONENT_CENSUS = {
    "cortex": {"n_filaments": 70686, "n_nodes": 494802, "n_faces": 12},
    "membrane": {"n_nodes": 5000, "n_faces": 9996},
    "sf_arc": {"n_filaments": 400, "n_nodes": 8000},
    "nmii": {"n_filaments": 2000, "n_nodes": 30000, "n_endpoints": 8000},
    "microtubule": {"n_filaments": 300, "n_nodes": 6000},
    "intermediate_filament": {"n_filaments": 1500, "n_nodes": 22000},
    "lamellipodium": {"n_filaments": 5000, "n_nodes": 40000},
    "filopodium": {"n_filaments": 60, "n_nodes": 1200},
}


def _runtimes(*, share_fa_component: bool = False) -> ComposedCellRuntimes:
    architecture = reference_cell_architecture()

    component_owners: dict[str, object] = {}
    for component in architecture.components:
        if not component.dynamically_evolving and not component.owns_geometry:
            component_owners[component.name] = object()  # world_boundary: inert reference frame
            continue
        component_owners[component.name] = _FullRuntime(**_COMPONENT_CENSUS.get(component.name, {}))

    # ONE runtime per declared mechanical group, DERIVED from the architecture.  This used to name the
    # single `alpha2beta1_collagen_series` pair, which meant every new series (the two nascent adhesions,
    # PI 2026-08-09 option (a)) silently bound two runtimes and tripped the composite guard.
    group_runtimes: dict[str, object] = {
        group: _FullRuntime(mechanical_group=group, n_endpoints=256)
        for group in sorted({c.mechanical_group for c in architecture.connectors if c.mechanical_group})
    }
    fa_joint = group_runtimes[_FA_GROUP]
    if share_fa_component:
        component_owners["focal_adhesion"] = fa_joint

    connector_runtimes: dict[str, object] = {}
    for connector in architecture.connectors:
        if connector.mechanical_group is not None:
            connector_runtimes[connector.name] = group_runtimes[connector.mechanical_group]
            continue
        connector_runtimes[connector.name] = _FullRuntime(
            name=connector.name,
            component_a=connector.component_a,
            component_b=connector.component_b,
        )

    facade_dispatch = {}
    for claim in canonical_facade_claims():
        method = _FACADE_METHOD[claim.owner_type]
        args = (object(), 0.05) if method == "candidate_iteration" else (object(),)
        facade_dispatch[claim.owner_type] = FacadeDispatch(_FacadeSpy(), args=args)

    return ComposedCellRuntimes(
        component_owners=component_owners,
        connector_runtimes=connector_runtimes,
        facade_dispatch=facade_dispatch,
    )


def test_build_registers_every_component_and_connector_in_one_world() -> None:
    world = build_composed_cell_world(_runtimes())
    architecture = reference_cell_architecture()

    assert set(world.registered_components()) == {c.name for c in architecture.components}
    assert set(world.registered_connectors()) == {c.name for c in architecture.connectors}
    # every geometry-owning component and every connector is dispatched exactly once by a canonical facade.
    for component in architecture.components:
        if component.owns_geometry:
            assert world.pipeline.component_owner(component.name) in {
                task.name for task in world.facade_tasks
            }
    for connector in architecture.connectors:
        assert world.pipeline.connector_owner(connector.name) in {
            task.name for task in world.facade_tasks
        }
    assert len(world.facade_tasks) == 8


def test_world_transaction_runs_one_predicate_across_all_participants() -> None:
    world = build_composed_cell_world(_runtimes())
    accepted, ledger = object(), object()

    world.transaction.begin_candidate()
    world.transaction.accumulate_ledgers(ledger)
    world.transaction.finalize_candidate(accepted, dt_phys=0.05, rng_seed=7)

    expected = [("snapshot",), ("ledger", ledger), ("rollback", accepted), ("commit", accepted, 0.05, 7)]
    participants = world.transaction.participants
    assert participants  # non-empty
    for participant in participants:
        assert participant.calls == expected


def test_pipeline_dispatch_calls_each_facade_once_with_its_bindings() -> None:
    world = build_composed_cell_world(_runtimes())
    world.pipeline.run_candidate()
    for task in world.facade_tasks:
        facade = task.runtime
        assert len(facade.calls) == 1  # each canonical facade orchestration method fires exactly once


def test_composite_fa_joint_is_one_participant_not_two() -> None:
    world = build_composed_cell_world(_runtimes())
    fa_runtime = world.actor.connector_runtime("fa_actin_anchor")
    assert world.actor.connector_runtime("integrin_collagen_clutch") is fa_runtime
    # the shared joint appears exactly once among the deduplicated transaction participants.
    assert sum(1 for participant in world.transaction.participants if participant is fa_runtime) == 1


def test_missing_facade_dispatch_is_rejected() -> None:
    runtimes = _runtimes()
    broken = dict(runtimes.facade_dispatch)
    del broken[NMIIActuator]
    with pytest.raises(ValueError, match="missing the dispatch facade"):
        build_composed_cell_world(
            ComposedCellRuntimes(runtimes.component_owners, runtimes.connector_runtimes, broken)
        )


def test_incomplete_world_requires_complete_flag() -> None:
    runtimes = _runtimes()
    partial_components = {k: v for k, v in runtimes.component_owners.items() if k != "microtubule"}
    partial = ComposedCellRuntimes(partial_components, runtimes.connector_runtimes, runtimes.facade_dispatch)
    # the exact-once dispatch coverage still holds, so a bring-up world builds without require_complete...
    world = build_composed_cell_world(partial, require_complete=False)
    assert "microtubule" not in world.registered_components()
    # ...but a production world must have every declared binding.
    with pytest.raises(ValueError, match="not fully bound"):
        build_composed_cell_world(partial, require_complete=True)


def test_dump_assigns_disjoint_global_actor_and_filament_ids() -> None:
    world = build_composed_cell_world(_runtimes())
    dump = dump_composed_world_state(world)

    assert dump.actor_ids_are_unique()
    assert dump.filament_id_blocks_are_disjoint()
    # SF is a separate actor from the cortex — different ids and non-overlapping global filament blocks.
    cortex = dump.by_name("cortex")
    sf = dump.by_name("sf_arc")
    assert cortex.actor_id != sf.actor_id
    assert cortex.n_filaments == 70686 and sf.n_filaments == 400
    lo_c, hi_c = cortex.filament_id_range
    lo_s, hi_s = sf.filament_id_range
    assert hi_c <= lo_s or hi_s <= lo_c
    # per-compartment census is exposed for the viz.
    assert cortex.n_nodes == 494802 and cortex.n_faces == 12
    assert dump.n_filaments_total == sum(
        e.n_filaments for e in dump.entries
    )


def test_dump_counts_composite_fa_joint_once_with_alias() -> None:
    world = build_composed_cell_world(_runtimes())
    dump = dump_composed_world_state(world)

    fa_entries = [e for e in dump.entries if "fa_actin_anchor" in (e.name, *e.aliases)]
    assert len(fa_entries) == 1
    entry = fa_entries[0]
    assert entry.name == "fa_actin_anchor"
    assert "integrin_collagen_clutch" in entry.aliases  # the second edge is an alias, not a second actor


def test_dump_shared_fa_component_and_connector_collapse_to_one_actor() -> None:
    world = build_composed_cell_world(_runtimes(share_fa_component=True))
    dump = dump_composed_world_state(world)
    # binding the FA joint as both the focal_adhesion component and its two connector edges is ONE actor.
    fa_entry = dump.by_name("focal_adhesion")
    assert "fa_actin_anchor" in fa_entry.aliases
    assert "integrin_collagen_clutch" in fa_entry.aliases
    assert dump.actor_ids_are_unique()
