"""Composed-world assembly for the component-first Active Cell engine.

The declarative :mod:`~aleph.engine.contracts` architecture, the exact-once
:func:`~aleph.engine.dispatch.canonical_facade_claims` manifest, the :class:`CellActor` binding
registry, and the :class:`CellWorldTransaction` accepted-step coordinator are each independently validated.
Nothing, however, *binds concrete runtime instances* into one discoverable composed world — every engine
test builds spies inline.  This module is that missing integration glue: it takes the concrete facade,
component state-owner, and connector runtime objects a build assembles (on the CUDA lane) and returns a
validated ``(actor, pipeline, transaction)`` triple.

Design invariants preserved (never re-decided here):

* **Components own state; connectors are the only mechanical connection.**  Component state-owners bind as
  ``CellActor`` components; connector runtimes bind as connectors.  Dispatch facades orchestrate the candidate
  schedule but are *not* transaction participants — the state-owners and connector runtimes are.
* **Exact-once dispatch.**  The candidate task list is generated directly from ``canonical_facade_claims()``,
  so each geometry-owning component and each connector edge is claimed by exactly one facade method, once.
* **Composite joints bind once.**  The α2β1–collagen FA series joint is one runtime object bound under both
  ``fa_actin_anchor`` and ``integrin_collagen_clutch``; :class:`CellActor` and :class:`CellWorldTransaction`
  deduplicate it, so it is never a spring-in-series double count.
* **One physical clock.**  The returned :class:`CellWorldTransaction` is the sole snapshot/rollback/commit
  coordinator; this module advances no time and reads no device predicate on the host.

This module allocates no device memory and launches no kernel; it is CPU-importable and the CUDA lane injects
the real Warp-resident runtimes.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from aleph.engine.actor import CellActor
from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.dispatch import (
    CandidatePhase,
    CandidateTask,
    CellCandidatePipeline,
    canonical_facade_claims,
    validate_canonical_facade_claims,
)
from aleph.engine.fluid_core import FluidCore
from aleph.engine.nmii_actuator import NMIIActuator
from aleph.engine.surface_body import SurfaceBody
from aleph.engine.world import CellWorldTransaction

__all__ = [
    "FacadeDispatch",
    "ComposedCellRuntimes",
    "ComposedCellWorld",
    "canonical_facade_phase",
    "build_composed_cell_world",
]


# Deterministic candidate-phase for each dispatch facade's single canonical orchestration method.  The
# exact-once gate keys on ``(runtime, phase)`` and on global claim uniqueness, so these choices set ordering,
# not correctness; a production build may override them.  Fields solve after solid mechanics assembles; the
# head-resolved NMII motor graph assembles with the other connectors it drives.
_CANONICAL_FACADE_PHASE: dict[str, CandidatePhase] = {
    "aleph.engine.surface_body.SurfaceBody": CandidatePhase.ASSEMBLE_INTERNAL,
    "aleph.engine.stress_fiber.StressFiberActor": CandidatePhase.ASSEMBLE_INTERNAL,
    "aleph.engine.microtubule_rig.MicrotubuleRig": CandidatePhase.ASSEMBLE_INTERNAL,
    "aleph.engine.intermediate_filament_rig.IntermediateFilamentRig": CandidatePhase.ASSEMBLE_INTERNAL,
    "aleph.engine.protrusion.ProtrusionActors": CandidatePhase.ASSEMBLE_INTERNAL,
    "aleph.engine.ecm_world.ECMWorld": CandidatePhase.ASSEMBLE_INTERNAL,
    "aleph.engine.nmii_actuator.NMIIActuator": CandidatePhase.ASSEMBLE_CONNECTORS,
    "aleph.engine.fluid_core.FluidCore": CandidatePhase.ASSEMBLE_FIELDS,
}


def canonical_facade_phase(facade_type: type[object]) -> CandidatePhase:
    """Return the default candidate phase for one dispatch facade type.

    Keyed by FULLY-QUALIFIED NAME rather than by the class object (2026-07-28).  Keying by the object
    forced this module to import every rig at import time purely to build a lookup table, and because
    ``ac/engine/__init__`` imports this module, a cortex-motor run that touches no microtubule, no ECM,
    no protrusion and no intermediate filament still loaded ~5,800 lines of them.  A qualified name is
    exactly as unambiguous as the object for this purpose, and an unregistered facade still raises.
    """
    key = f"{facade_type.__module__}.{facade_type.__qualname__}"
    try:
        return _CANONICAL_FACADE_PHASE[key]
    except KeyError as error:
        raise KeyError(
            f"no canonical candidate phase registered for facade {facade_type.__name__!r}"
        ) from error


@dataclass(frozen=True, slots=True)
class FacadeDispatch:
    """One concrete dispatch facade instance plus the forwarded step arguments for its canonical method.

    ``facade`` is the object whose ``canonical_facade_claims`` method (e.g. ``accumulate_mechanics``) claims a
    set of component/connector edges.  ``args``/``kwargs`` are the per-step bindings forwarded to that method
    (a ``*StepBindings`` and, for the fluid field iteration, the physical ``dt``).  ``phase`` defaults to the
    canonical ordering for the facade type.
    """

    facade: object
    args: tuple[Any, ...] = ()
    kwargs: tuple[tuple[str, Any], ...] = ()
    phase: CandidatePhase | None = None


@dataclass(frozen=True, slots=True)
class ComposedCellRuntimes:
    """The concrete runtime objects a build assembles, keyed by declared architecture name / facade type.

    ``component_owners`` maps every declared component name to its state-owning runtime; ``connector_runtimes``
    maps every declared connector name to its connector runtime (a composite joint appears under each of its
    semantic edge names as the *same* object); ``facade_dispatch`` maps each of the eight canonical dispatch
    facade types to its :class:`FacadeDispatch`.  All three are supplied by the CUDA assembly lane; this module
    only wires and validates them.
    """

    component_owners: Mapping[str, object]
    connector_runtimes: Mapping[str, object]
    facade_dispatch: Mapping[type[object], FacadeDispatch]


@dataclass(frozen=True, slots=True)
class ComposedCellWorld:
    """A fully wired composed cell: binding registry, exact-once candidate schedule, and step transaction."""

    architecture: CellArchitecture
    actor: CellActor
    pipeline: CellCandidatePipeline
    transaction: CellWorldTransaction
    _facade_tasks: tuple[CandidateTask, ...] = field(default_factory=tuple)

    @property
    def participants(self) -> tuple[object, ...]:
        """Return the deduplicated transaction participants (state-owners + connector runtimes)."""
        return self.transaction.participants

    @property
    def facade_tasks(self) -> tuple[CandidateTask, ...]:
        """Return the per-facade canonical candidate tasks in registration order."""
        return self._facade_tasks

    def registered_components(self) -> tuple[str, ...]:
        """Return declared components that are bound to a concrete runtime."""
        missing = set(self.actor.missing_bindings()["components"])
        return tuple(c.name for c in self.architecture.components if c.name not in missing)

    def registered_connectors(self) -> tuple[str, ...]:
        """Return declared connectors that are bound to a concrete runtime."""
        missing = set(self.actor.missing_bindings()["connectors"])
        return tuple(c.name for c in self.architecture.connectors if c.name not in missing)


def build_composed_cell_world(
    runtimes: ComposedCellRuntimes,
    architecture: CellArchitecture | None = None,
    *,
    require_complete: bool = True,
) -> ComposedCellWorld:
    """Bind concrete runtimes into a validated ``(actor, pipeline, transaction)`` composed cell world.

    Args:
        runtimes: the concrete facade / component / connector runtime objects to register.
        architecture: the declarative composition (defaults to :func:`reference_cell_architecture`).
        require_complete: when ``True`` (production), assert every declared component and connector is bound
            with a full transaction/ledger/mechanics API before returning; when ``False``, allow an incomplete
            bring-up world (the exact-once dispatch coverage gate is still enforced).

    Returns:
        A :class:`ComposedCellWorld` whose ``pipeline`` dispatches each geometry-owning component and each
        connector edge exactly once (from ``canonical_facade_claims``) and whose ``transaction`` snapshots,
        rolls back, and commits the deduplicated state-owners/connector runtimes under one accepted predicate.

    Raises:
        ValueError: if a facade dispatch spec is missing, or the exact-once/coverage gates reject the wiring.
    """
    architecture = architecture or reference_cell_architecture()
    # The static facade manifest must match the architecture before we instantiate anything against it.
    validate_canonical_facade_claims(architecture)

    actor = CellActor(architecture)
    for component in architecture.components:
        runtime = runtimes.component_owners.get(component.name)
        if runtime is not None:
            actor.bind_component(component.name, runtime)
    for connector in architecture.connectors:
        runtime = runtimes.connector_runtimes.get(connector.name)
        if runtime is not None:
            actor.bind_connector(connector.name, runtime)

    tasks = _facade_candidate_tasks(runtimes.facade_dispatch)
    pipeline = CellCandidatePipeline(architecture, tasks)
    transaction = CellWorldTransaction(actor, require_complete=require_complete)
    return ComposedCellWorld(
        architecture=architecture,
        actor=actor,
        pipeline=pipeline,
        transaction=transaction,
        _facade_tasks=tasks,
    )


def _facade_candidate_tasks(
    facade_dispatch: Mapping[type[object], FacadeDispatch],
) -> tuple[CandidateTask, ...]:
    """Generate the exact-once candidate task list from the canonical facade manifest + concrete facades."""
    tasks: list[CandidateTask] = []
    for claim in canonical_facade_claims():
        spec = facade_dispatch.get(claim.owner_type)
        if spec is None:
            raise ValueError(
                f"composed world is missing the dispatch facade for {claim.owner_type.__name__!r}"
            )
        # The facade must expose the manifest's canonical method (production passes a real facade instance;
        # a CUDA-free wiring harness passes a double that implements the same method).  The CandidateTask
        # constructor re-checks callability; this gives a clearer composed-world error first.
        if not callable(getattr(spec.facade, claim.method_name, None)):
            raise ValueError(
                f"facade dispatch for {claim.owner_type.__name__!r} has no callable "
                f"{claim.method_name!r} method"
            )
        phase = spec.phase if spec.phase is not None else canonical_facade_phase(claim.owner_type)
        tasks.append(
            CandidateTask(
                name=claim.owner_type.__name__,
                phase=phase,
                runtime=spec.facade,
                method_name=claim.method_name,
                args=spec.args,
                kwargs=spec.kwargs,
                component_claims=claim.component_claims,
                connector_claims=claim.connector_claims,
            )
        )
    return tuple(tasks)
