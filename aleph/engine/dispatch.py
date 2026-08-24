"""Exact-once candidate-phase dispatch for the component-first cell engine.

The component/connector graph is only declarative until a composed world assigns every geometry-owning
component and every connector edge to one callable candidate task.  This module makes that assignment a
build-time invariant.  It owns no physical state, performs no device reads, and does not advance time.

A composite joint may claim multiple semantic connector edges in one task.  Conversely, the same runtime
object cannot be scheduled twice in one phase; a facade that needs several calls in a phase must expose one
canonical orchestration method.  This prevents the FA series joint and other shared runtimes from silently
double-counting force or energy.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

from aleph.engine.contracts import CellArchitecture


class CandidatePhase(StrEnum):
    """Deterministic phases inside one uncommitted physical-step candidate."""

    DISCOVER_GEOMETRY = "discover_geometry"
    ASSEMBLE_INTERNAL = "assemble_internal"
    ASSEMBLE_CONNECTORS = "assemble_connectors"
    ASSEMBLE_FIELDS = "assemble_fields"
    SOLVE_COUPLED = "solve_coupled"
    PROPOSE_KINETICS = "propose_kinetics"
    EVALUATE_GATES = "evaluate_gates"


@dataclass(frozen=True, slots=True)
class CandidateTask:
    """One exact runtime method call plus its authoritative graph coverage claims.

    ``component_claims`` identify component mechanics owned by this task. ``connector_claims`` identify
    semantic graph edges whose mechanics/field transfer this method actually dispatches. Extra helper tasks
    in another phase may leave both tuples empty, but every required claim must occur exactly once globally.
    """

    name: str
    phase: CandidatePhase
    runtime: object
    method_name: str
    args: tuple[Any, ...] = ()
    kwargs: tuple[tuple[str, Any], ...] = ()
    component_claims: tuple[str, ...] = ()
    connector_claims: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("candidate task name must be non-empty")
        if not isinstance(self.phase, CandidatePhase):
            raise TypeError("candidate task phase must be a CandidatePhase")
        if not self.method_name.strip():
            raise ValueError(f"candidate task {self.name!r} needs a method name")
        if not callable(getattr(self.runtime, self.method_name, None)):
            raise TypeError(
                f"candidate task {self.name!r} runtime has no callable {self.method_name!r}"
            )
        if len(set(self.component_claims)) != len(self.component_claims):
            raise ValueError(f"candidate task {self.name!r} repeats a component claim")
        if len(set(self.connector_claims)) != len(self.connector_claims):
            raise ValueError(f"candidate task {self.name!r} repeats a connector claim")
        keyword_names = tuple(name for name, _ in self.kwargs)
        if len(set(keyword_names)) != len(keyword_names):
            raise ValueError(f"candidate task {self.name!r} repeats a keyword argument")

    def dispatch(self) -> None:
        """Invoke the registered candidate method without inspecting device state."""
        method = getattr(self.runtime, self.method_name)
        method(*self.args, **dict(self.kwargs))


@dataclass(frozen=True, slots=True)
class CanonicalFacadeClaim:
    """Static proof that one real facade method owns named component/connector dispatch slots."""

    owner_type: type[object]
    method_name: str
    component_claims: tuple[str, ...]
    connector_claims: tuple[str, ...]

    def __post_init__(self) -> None:
        if not callable(getattr(self.owner_type, self.method_name, None)):
            raise TypeError(
                f"canonical facade {self.owner_type.__name__} has no callable "
                f"{self.method_name!r}"
            )


def canonical_facade_claims() -> tuple[CanonicalFacadeClaim, ...]:
    """Return the exact common-graph ownership manifest backed by real facade methods."""
    from aleph.engine.ecm_world import ECMWorld
    from aleph.engine.fluid_core import FluidCore
    from aleph.engine.intermediate_filament_rig import IntermediateFilamentRig
    from aleph.engine.microtubule_rig import MicrotubuleRig
    from aleph.engine.nmii_actuator import NMIIActuator
    from aleph.engine.protrusion import ProtrusionActors
    from aleph.engine.stress_fiber import StressFiberActor
    from aleph.engine.surface_body import SurfaceBody

    return (
        CanonicalFacadeClaim(
            SurfaceBody,
            "accumulate_mechanics",
            ("membrane", "cortex"),
            (
                "membrane_erm_cortex",
                "surface_porous_transfer",
                "membrane_cytosol_boundary",
                # PI-ratified 2026-07-25. Both are surface-surface non-penetration couplers whose shared
                # endpoint is the cortex, which this facade owns. membrane_cortex_contact is the compressive
                # channel the unilateral ERM tether structurally cannot supply (SURFACE_BODY_PLAN.md:257
                # names this facade as its owner); nucleus_cortex_contact is the only compressive path
                # between the shell and the envelope, every other nucleus solid edge being a tether.
                "membrane_cortex_contact",
                "nucleus_cortex_contact",
                # PI D5-A, 2026-07-28. The membrane's OUTWARD face: the free face of the asymmetric
                # world boundary (basal = 2D ECM + far-field anchor, free = media). This facade owns
                # the membrane surface quadrature the traction is applied over, so it owns the claim.
                # The exterior solve landed with T10 (`medium_exterior.ExteriorStokesMedium`) and its
                # connector identity with `MembraneMediumTraction` (2026-08-09), so this claim now does
                # assert that a medium traction is evaluated — the sentence that stood here said it did
                # not. What it still does NOT assert is a MAGNITUDE: `mu_medium` is a PI-GAP.
                "membrane_medium_traction",
            ),
        ),
        CanonicalFacadeClaim(
            FluidCore,
            "candidate_iteration",
            ("cytosol", "nucleus"),
            ("nucleus_cytosol_boundary",),
        ),
        CanonicalFacadeClaim(
            StressFiberActor,
            "accumulate_mechanics",
            ("sf_arc",),
            (
                "sf_cortex_transient",
                "actin_cap_linc",
                "sf_cytosol_transfer",
                "dorsal_arc_crosslink",
            ),
        ),
        CanonicalFacadeClaim(
            MicrotubuleRig,
            "accumulate_mechanics",
            ("microtubule",),
            (
                "mt_nucleus_linc",
                "mt_cortex_capture",
                "mt_sf_spectraplakin",
                "mt_cytosol_transfer",
            ),
        ),
        CanonicalFacadeClaim(
            IntermediateFilamentRig,
            "accumulate_mechanics",
            ("intermediate_filament",),
            ("if_nucleus_linc", "if_sf_plectin", "if_cytosol_transfer"),
        ),
        CanonicalFacadeClaim(
            ProtrusionActors,
            "accumulate_mechanics",
            ("lamellipodium", "filopodium"),
            (
                "lamellipodium_membrane_contact",
                "lamellipodium_cortex_seam",
                "lamellipodium_cytosol_transfer",
                "lamellipodium_nascent_fa",
                "filopodium_membrane_tip",
                "filopodium_cortex_root",
                "filopodium_cytosol_transfer",
                "filopodium_nascent_fa",
            ),
        ),
        CanonicalFacadeClaim(
            NMIIActuator,
            "accumulate_candidate",
            ("nmii",),
            (
                "nmii_sf_motor",
                "nmii_cortex_motor",
                "nmii_lamellipodium_motor",
                "nmii_filopodium_motor",
                # PI-ratified 2026-07-25: nmii was the only geometry-owning dynamic component with no fluid
                # coupling, so an all-heads-unbound minifilament was a free rigid body.
                "nmii_cytosol_transfer",
            ),
        ),
        CanonicalFacadeClaim(
            ECMWorld,
            "accumulate_mechanics",
            ("ecm",),
            (
                "fa_actin_anchor",
                "integrin_collagen_clutch",
                # PI 2026-08-09 option (a): each nascent adhesion is its own SERIES.  Only the
                # LIGAND-side halves are claimed here — the actin-side anchors already belong to the
                # protrusion facade that owns their filament end, and a connector may have exactly one
                # canonical owner.  The composite spring runs actin -> collagen either way, and ECM owns
                # the collagen end of all three series.
                "lamellipodium_nascent_clutch",
                "filopodium_nascent_clutch",
                "ecm_crosslink",
                "ecm_far_field_anchor",
                "membrane_ecm_contact",
            ),
        ),
    )


def validate_canonical_facade_claims(architecture: CellArchitecture) -> None:
    """Reject stale, missing, or duplicate facade ownership in the common dispatch manifest."""
    component_owners: dict[str, str] = {}
    connector_owners: dict[str, str] = {}
    for claim in canonical_facade_claims():
        owner = claim.owner_type.__name__
        for component in claim.component_claims:
            if component in component_owners:
                raise ValueError(f"component {component!r} has duplicate canonical facade owners")
            component_owners[component] = owner
        for connector in claim.connector_claims:
            if connector in connector_owners:
                raise ValueError(f"connector {connector!r} has duplicate canonical facade owners")
            connector_owners[connector] = owner
    expected_components = {
        component.name for component in architecture.components if component.owns_geometry
    }
    expected_connectors = {connector.name for connector in architecture.connectors}
    if set(component_owners) != expected_components or set(connector_owners) != expected_connectors:
        raise ValueError(
            "canonical facade manifest does not match the architecture; "
            f"missing_components={sorted(expected_components - set(component_owners))}, "
            f"extra_components={sorted(set(component_owners) - expected_components)}, "
            f"missing_connectors={sorted(expected_connectors - set(connector_owners))}, "
            f"extra_connectors={sorted(set(connector_owners) - expected_connectors)}"
        )


@dataclass(frozen=True, slots=True)
class CellCandidatePipeline:
    """Validated exact-once mechanics/proposal schedule for one cell architecture."""

    architecture: CellArchitecture
    tasks: tuple[CandidateTask, ...]
    _ordered_tasks: tuple[CandidateTask, ...] = field(init=False, repr=False)
    _component_owner: dict[str, str] = field(init=False, repr=False)
    _connector_owner: dict[str, str] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not self.tasks:
            raise ValueError("cell candidate pipeline requires at least one task")
        task_names = tuple(task.name for task in self.tasks)
        if len(set(task_names)) != len(task_names):
            raise ValueError("candidate task names must be unique")

        known_components = {component.name for component in self.architecture.components}
        known_connectors = {connector.name for connector in self.architecture.connectors}
        required_components = {
            component.name for component in self.architecture.components if component.owns_geometry
        }
        component_owner: dict[str, str] = {}
        connector_owner: dict[str, str] = {}
        runtime_phase_owner: dict[tuple[int, CandidatePhase], str] = {}

        for task in self.tasks:
            unknown_components = set(task.component_claims) - known_components
            if unknown_components:
                raise ValueError(
                    f"candidate task {task.name!r} claims unknown components "
                    f"{sorted(unknown_components)}"
                )
            unknown_connectors = set(task.connector_claims) - known_connectors
            if unknown_connectors:
                raise ValueError(
                    f"candidate task {task.name!r} claims unknown connectors "
                    f"{sorted(unknown_connectors)}"
                )
            for component in task.component_claims:
                previous = component_owner.setdefault(component, task.name)
                if previous != task.name:
                    raise ValueError(
                        f"component {component!r} mechanics is claimed by both "
                        f"{previous!r} and {task.name!r}"
                    )
            for connector in task.connector_claims:
                previous = connector_owner.setdefault(connector, task.name)
                if previous != task.name:
                    raise ValueError(
                        f"connector {connector!r} dispatch is claimed by both "
                        f"{previous!r} and {task.name!r}"
                    )
            runtime_phase = (id(task.runtime), task.phase)
            previous_task = runtime_phase_owner.setdefault(runtime_phase, task.name)
            if previous_task != task.name:
                raise ValueError(
                    f"runtime {type(task.runtime).__name__} is scheduled more than once in "
                    f"phase {task.phase.value!r}: {previous_task!r}, {task.name!r}"
                )

        missing_components = required_components - set(component_owner)
        missing_connectors = known_connectors - set(connector_owner)
        if missing_components or missing_connectors:
            raise ValueError(
                "candidate dispatch coverage is incomplete; "
                f"missing_components={sorted(missing_components)}, "
                f"missing_connectors={sorted(missing_connectors)}"
            )

        phase_rank = {phase: index for index, phase in enumerate(CandidatePhase)}
        ordered = tuple(
            task
            for _, task in sorted(
                enumerate(self.tasks),
                key=lambda item: (phase_rank[item[1].phase], item[0]),
            )
        )
        object.__setattr__(self, "_ordered_tasks", ordered)
        object.__setattr__(self, "_component_owner", component_owner)
        object.__setattr__(self, "_connector_owner", connector_owner)

    @property
    def ordered_tasks(self) -> tuple[CandidateTask, ...]:
        """Return tasks in deterministic phase then registration order."""
        return self._ordered_tasks

    def component_owner(self, name: str) -> str:
        """Return the task that claims one component's mechanics."""
        return self._component_owner[name]

    def connector_owner(self, name: str) -> str:
        """Return the task that dispatches one semantic connector edge."""
        return self._connector_owner[name]

    def run_phase(self, phase: CandidatePhase) -> None:
        """Dispatch one candidate phase once in deterministic order."""
        if not isinstance(phase, CandidatePhase):
            raise TypeError("phase must be a CandidatePhase")
        for task in self._ordered_tasks:
            if task.phase is phase:
                task.dispatch()

    def run_candidate(self) -> None:
        """Dispatch one uncommitted pass through all candidate phases.

        The enclosing nonlinear/world solver may invoke individual phases again as part of convergence. This
        method does not decide acceptance, commit state, or advance physical time.
        """
        for task in self._ordered_tasks:
            task.dispatch()


@dataclass(frozen=True, slots=True)
class ConnectorRuntimeBinding:
    """Which concrete class supplies one declared connector's mechanics.

    Distinct from :class:`CanonicalFacadeClaim`, and the two answer different questions.  A facade
    claim says WHO DISPATCHES an edge inside a candidate phase, and the manifest above has been
    complete for every declared edge since long before every edge had an implementation — because a
    facade can claim an edge whose runtime does not exist.  This says WHAT IS BOUND at the edge, and
    it is the only place the repo can answer "how many connectors actually have a runtime" without a
    session counting by hand.  It was counted by hand until 2026-08-09, and the number in `STATE.md`
    was therefore an assertion rather than a measurement.

    ``symbol`` is verified by IMPORT in ``test_dispatch``, and verified to be a concrete class rather
    than a ``typing.Protocol``: several connector seams in this tree are declared as Protocols with no
    implementation behind them (``ECMMembraneContact`` was one until this census was written), and a
    census that could not tell those apart would report the seam as a runtime.
    """

    connector: str
    module: str
    symbol: str
    #: Set when several semantic edges resolve to ONE composite runtime under a ``mechanical_group``.
    mechanical_group: str | None = None


def connector_runtime_census() -> tuple[ConnectorRuntimeBinding, ...]:
    """Return the concrete runtime bound at every declared connector, in architecture order.

    Hand-maintained, and deliberately in ONE place: the runtime modules do not declare which edges
    they serve — three of them now do (``immersed_transfer_edges``, ``filament_crosslink_edges``,
    ``connector_joint_edges``) and the rest predate the idea.  A single table that
    :func:`validate_connector_runtime_census` checks against the architecture, and whose symbols the
    gate imports, fails loudly on drift in either direction; the same list spread across ten modules
    would not.
    """
    engine = "aleph.engine."
    return tuple(
        ConnectorRuntimeBinding(connector=name, module=engine + module, symbol=symbol, mechanical_group=group)
        for name, module, symbol, group in (
            # --- surface / field couplers -----------------------------------------------------------
            ("membrane_erm_cortex", "erm_cortex_connector", "ErmCortexConnector", None),
            ("membrane_cytosol_boundary", "cytosol_connected", "MovingSemipermeableFluidBoundaryFacade", None),
            ("nucleus_cytosol_boundary", "fluid_core", "MovingImpermeableFluidBoundaryFacade", None),
            # --- immersed transfer: one runtime, seven edges (2026-08-09) ----------------------------
            ("surface_porous_transfer", "immersed_transfer", "ImmersedPorousTransfer", None),
            ("sf_cytosol_transfer", "immersed_transfer", "ImmersedPorousTransfer", None),
            ("mt_cytosol_transfer", "immersed_transfer", "ImmersedPorousTransfer", None),
            ("if_cytosol_transfer", "immersed_transfer", "ImmersedPorousTransfer", None),
            ("lamellipodium_cytosol_transfer", "immersed_transfer", "ImmersedPorousTransfer", None),
            ("filopodium_cytosol_transfer", "immersed_transfer", "ImmersedPorousTransfer", None),
            ("nmii_cytosol_transfer", "immersed_transfer", "ImmersedPorousTransfer", None),
            # --- filament crosslinks: one runtime, three families, five edges (2026-08-09) -----------
            ("sf_cortex_transient", "filament_crosslink", "FilamentCrosslinkConnector", None),
            ("lamellipodium_cortex_seam", "filament_crosslink", "FilamentCrosslinkConnector", None),
            ("filopodium_cortex_root", "filament_crosslink", "FilamentCrosslinkConnector", None),
            ("if_sf_plectin", "filament_crosslink", "FilamentCrosslinkConnector", None),
            ("mt_sf_spectraplakin", "filament_crosslink", "FilamentCrosslinkConnector", None),
            # --- contact + non-NMII motor: one runtime, two force laws, six edges (2026-08-09) -------
            ("membrane_cortex_contact", "connector_joints", "ContractJointConnector", None),
            ("membrane_ecm_contact", "connector_joints", "ContractJointConnector", None),
            ("nucleus_cortex_contact", "connector_joints", "ContractJointConnector", None),
            ("lamellipodium_membrane_contact", "connector_joints", "ContractJointConnector", None),
            ("filopodium_membrane_tip", "connector_joints", "ContractJointConnector", None),
            ("mt_cortex_capture", "connector_joints", "ContractJointConnector", None),
            # --- LINC ------------------------------------------------------------------------------
            ("actin_cap_linc", "linc_connector", "LincRigConnectorAdapter", None),
            ("mt_nucleus_linc", "linc_connector", "LincRigConnectorAdapter", None),
            ("if_nucleus_linc", "linc_connector", "LincRigConnectorAdapter", None),
            # --- NMII head-resolved crossbridges -----------------------------------------------------
            ("nmii_sf_motor", "cortex_motor_slice", "FilamentMotorConnector", None),
            ("nmii_cortex_motor", "cortex_motor_slice", "FilamentMotorConnector", None),
            ("nmii_lamellipodium_motor", "cortex_motor_slice", "FilamentMotorConnector", None),
            ("nmii_filopodium_motor", "cortex_motor_slice", "FilamentMotorConnector", None),
            # --- internal joints ---------------------------------------------------------------------
            ("dorsal_arc_crosslink", "sf_mechanics", "SFInternalArcJointConnector", None),
            ("ecm_crosslink", "collagen_crosslink", "CollagenCrosslinkConnector", None),
            # --- adhesion series: two semantic edges, ONE composite spring each ------------------------
            ("fa_actin_anchor", "ecm_world", "CompositeLoadPathECMClutchAdapter", "alpha2beta1_collagen_series"),
            ("integrin_collagen_clutch", "ecm_world", "CompositeLoadPathECMClutchAdapter", "alpha2beta1_collagen_series"),
            ("lamellipodium_nascent_fa", "connector_joints", "ContractJointConnector", "lamellipodium_nascent_series"),
            ("lamellipodium_nascent_clutch", "connector_joints", "ContractJointConnector", "lamellipodium_nascent_series"),
            ("filopodium_nascent_fa", "connector_joints", "ContractJointConnector", "filopodium_nascent_series"),
            ("filopodium_nascent_clutch", "connector_joints", "ContractJointConnector", "filopodium_nascent_series"),
            # --- environment boundaries ---------------------------------------------------------------
            ("ecm_far_field_anchor", "ecm_world", "ECMBoundaryAnchorFacade", None),
            ("membrane_medium_traction", "medium_exterior", "MembraneMediumTraction", None),
        )
    )


def validate_connector_runtime_census(architecture: CellArchitecture) -> None:
    """Reject a census that does not name a runtime for every declared connector, exactly once.

    Raises:
        ValueError: on a connector with no entry, an entry naming an undeclared connector, a duplicate
            entry, or an entry whose ``mechanical_group`` disagrees with the declaration it is for.
    """
    census = connector_runtime_census()
    seen: dict[str, ConnectorRuntimeBinding] = {}
    for binding in census:
        if binding.connector in seen:
            raise ValueError(f"connector {binding.connector!r} has duplicate runtime census entries")
        seen[binding.connector] = binding
    declared = {connector.name: connector for connector in architecture.connectors}
    if set(seen) != set(declared):
        raise ValueError(
            "connector runtime census does not match the architecture; "
            f"without_a_runtime={sorted(set(declared) - set(seen))}, "
            f"not_declared={sorted(set(seen) - set(declared))}"
        )
    for name, binding in seen.items():
        if binding.mechanical_group != declared[name].mechanical_group:
            raise ValueError(
                f"connector {name!r} census entry claims mechanical group "
                f"{binding.mechanical_group!r} but the contract declares "
                f"{declared[name].mechanical_group!r}"
            )


__all__ = [
    "CandidatePhase",
    "CandidateTask",
    "CanonicalFacadeClaim",
    "CellCandidatePipeline",
    "ConnectorRuntimeBinding",
    "canonical_facade_claims",
    "connector_runtime_census",
    "validate_canonical_facade_claims",
    "validate_connector_runtime_census",
]
