"""Two-array adjoint filament crosslinks: one force law, three families, five declared edges.

WHAT THIS CLOSES.  Five declared edges across three families are the SAME mechanics — a state-gated
Hookean central force between two interpolated segment points, in two components that own DISJOINT
arrays::

    sf_cortex_transient        transient_actin  sf_arc                -> cortex
    lamellipodium_cortex_seam  transient_actin  lamellipodium         -> cortex
    filopodium_cortex_root     transient_actin  filopodium            -> cortex
    if_sf_plectin              plectin          intermediate_filament -> sf_arc
    mt_sf_spectraplakin        spectraplakin    microtubule           -> sf_arc

The family label is biology; the force law is one.  Writing three near-copies of it would be three
places for the adjoint scatter to drift apart, and the adjoint is the property the balance gate exists
to check — so it is written once and the ``JointKind`` carries which biology a joint is, for the
diagnostic rather than for a branch in the kernel.

The ``transient_actin`` family also holds a FOURTH edge this module deliberately refuses:
``dorsal_arc_crosslink`` (``sf_arc -> sf_arc``, INTERNAL) has both endpoints in ONE array and binds
``link_spring_kernel`` directly over a ``(J, 2)`` pair table.  That is other physics and it already has
a runtime (:class:`~aleph.engine.sf_mechanics.SFInternalArcJointConnector`).

An inter-component joint must gather from two position arrays and scatter the equal-and-opposite pair
back into two force arrays.  ``sf_mechanics`` says exactly this in its own docstring and calls its
inter-component edge SEAMED for that reason.

**What differs across the five is CHEMISTRY, not mechanics.**  Five edges declare four chemistry cards
(``transient_actin_crosslink`` twice, then ``formin_fascin_root_coupling``, ``plectin_if_actin``,
``spectraplakin_mt_actin``).  Chemistry sets when a bond forms or breaks — when ``candidate_state``
flips — and does not change the force law.  So the rate law is injected, not baked in.

**Almost none of this is new code.**  :class:`~aleph.engine.load_path.LoadPathJointRuntime` is already a
fixed-capacity two-array joint SoA with the full transaction API, and
``_segment_joint_force_kernel`` is already generic: two positions, two forces, two segment tables,
barycentric weights, and a state gate that zeroes the force below ``JointState.ACTIN_ENGAGED``.  Nothing
in it is focal-adhesion-specific.  What was missing was a spec, a resolver and a ``JointKind`` — so that
is what this module adds, plus the thin edge-identity wrapper a connector slot needs.

WHY THE KINETICS DELEGATE IS REQUIRED RATHER THAN DEFAULTED.  All three edges declare ``kinetics=True``
and ``commit_on_accept=True``.  A transient crosslink with no rate law is not a transient crosslink; it is
a permanent weld, which is the one thing the charter says a connector may never be.  So
:meth:`FilamentCrosslinkConnector.propose_events` RAISES when no delegate is bound rather than
silently proposing nothing — binding the mechanical edge is allowed and useful on its own, but it may
never be mistaken for a running kinetic one.

KNOWN GAP, recorded rather than invented.  :class:`~aleph.engine.load_path.EndpointRole` declares roles
for stress-fibre and collagen endpoints only; there is no member for a lamellipodial seam or a filopodial
root.  Adding one is a contract change and therefore a PI decision, so this resolver validates the
endpoint COMPONENTS against the declared contract — which is the part that can be checked today — and
carries the caller's role through into the SoA without asserting a role it has no declaration for.

Sanity Gate (before first execution; :mod:`aleph.tests.ac.engine.test_transient_actin`):
  * dimensions: stiffness pN/um, rest length um, ``u`` a dimensionless barycentric in [0, 1].  This module
    introduces no constant; every magnitude arrives on the spec.
  * boundary cases: a joint list of length zero is refused by the joint runtime (an empty connector is a
    binding mistake, not a configuration); an edge name outside the family, or of INTERNAL scope, raises.
  * conservation invariant: the force law is central and equal-and-opposite by construction — the kernel
    scatters ``+f`` through A's barycentric weights and ``-f`` through B's, so the pair sums to zero
    exactly, at any configuration.  That is the adjoint property the balance gate tests.
  * sign sense: a stretched joint pulls its endpoints together (rest length is the zero-force separation);
    a joint below ``ACTIN_ENGAGED`` carries no force at all, so an unbound crosslink cannot pull.
  * precision / measurement protocol: no host readback; the connector produces no magnitude of its own.

CPU-importable: it allocates nothing itself and defers every allocation and launch to the joint runtime,
so the structural gate drives it with recording doubles.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorFamily,
    ConnectorScope,
    reference_cell_architecture,
)
from aleph.engine.load_path import (
    ActorRegistry,
    JointKind,
    JointState,
    LoadPathJointRuntime,
    PortRef,
    ResolvedJoint,
)

__all__ = [
    "CROSSLINK_FAMILIES",
    "FILAMENT_CROSSLINK_EDGES",
    "FilamentCrosslinkConnector",
    "FilamentCrosslinkJointSpec",
    "build_filament_crosslink_connector",
    "filament_crosslink_edges",
    "resolve_filament_crosslink_joint",
]

#: Families whose inter-component edges are this one force law.  A family is admitted here only when its
#: edges join two GEOMETRIC material points; ``actin_anchor`` is absent for that reason — its nascent
#: edges end on ``focal_adhesion``, which declares ``owns_geometry=False``, so there is nothing to pull
#: against until the endpoint shape is decided (an open PI item, not an omission).
CROSSLINK_FAMILIES = (
    ConnectorFamily.TRANSIENT_ACTIN,
    ConnectorFamily.PLECTIN,
    ConnectorFamily.SPECTRAPLAKIN,
)

#: Which ``JointKind`` the SoA records per family, so a diagnostic can tell the biologies apart.
_FAMILY_JOINT_KIND = {
    ConnectorFamily.TRANSIENT_ACTIN: JointKind.TRANSIENT_ACTIN_CROSSLINK,
    ConnectorFamily.PLECTIN: JointKind.PLECTIN_IF_ACTIN,
    ConnectorFamily.SPECTRAPLAKIN: JointKind.SPECTRAPLAKIN_MT_ACTIN,
}


def filament_crosslink_edges(architecture: CellArchitecture | None = None) -> tuple[str, ...]:
    """Return every inter-component edge this runtime serves, in architecture order.

    Read from the contract, so a new component that declares a crosslink in one of
    :data:`CROSSLINK_FAMILIES` is covered without editing this module.  INTERNAL edges are excluded —
    both their endpoints share one array, which is other physics with its own runtime.
    """
    arch = architecture or reference_cell_architecture()
    return tuple(
        c.name for c in arch.connectors
        if c.family in CROSSLINK_FAMILIES and c.scope is ConnectorScope.INTER_COMPONENT
    )


#: Convenience snapshot at import time.
FILAMENT_CROSSLINK_EDGES = filament_crosslink_edges()


@dataclass(frozen=True, slots=True)
class FilamentCrosslinkJointSpec:
    """One transient crosslink between an actin body's material point and a cortical material point."""

    joint_id: int
    edge: str
    port_a: PortRef
    port_b: PortRef
    stiffness_pn_per_um: float
    rest_um: float
    initial_state: JointState = JointState.ACTIN_ENGAGED

    def __post_init__(self) -> None:
        if self.joint_id < 0:
            raise ValueError("joint ID must be nonnegative")
        if not np.isfinite(self.stiffness_pn_per_um) or self.stiffness_pn_per_um <= 0.0:
            raise ValueError("joint stiffness must be finite and positive")
        if not np.isfinite(self.rest_um) or self.rest_um < 0.0:
            raise ValueError("joint rest length must be finite and nonnegative")


def _contract(edge: str, architecture: CellArchitecture):
    contract = next((c for c in architecture.connectors if c.name == edge), None)
    if contract is None:
        raise ValueError(f"connector {edge!r} is not declared in the architecture")
    if contract.family not in CROSSLINK_FAMILIES:
        admitted = ", ".join(sorted(f.value for f in CROSSLINK_FAMILIES))
        raise ValueError(
            f"connector {edge!r} is family {contract.family.value!r}; this runtime serves {admitted}"
        )
    if contract.scope is not ConnectorScope.INTER_COMPONENT:
        raise ValueError(
            f"connector {edge!r} is INTERNAL — both endpoints share one array, so it binds a "
            "single-array joint runtime, not the two-array adjoint one"
        )
    return contract


def resolve_filament_crosslink_joint(
    spec: FilamentCrosslinkJointSpec,
    registry: ActorRegistry,
    architecture: CellArchitecture,
) -> ResolvedJoint:
    """Validate one transient body<->cortex crosslink against the declared edge and resolve it.

    Endpoint COMPONENTS are checked against the contract (which is what can be checked today); endpoint
    ROLES are carried through unasserted — see the module docstring's recorded gap.

    Raises:
        ValueError: if the edge is not an inter-component ``transient_actin`` contract, or if either port
            belongs to a component the contract does not name.
    """
    contract = _contract(spec.edge, architecture)
    registry.validate_port(spec.port_a)
    registry.validate_port(spec.port_b)
    if spec.port_a.component != contract.component_a:
        raise ValueError(
            f"{spec.edge!r} endpoint A must belong to {contract.component_a!r}, "
            f"got {spec.port_a.component!r}"
        )
    if spec.port_b.component != contract.component_b:
        raise ValueError(
            f"{spec.edge!r} endpoint B must belong to {contract.component_b!r}, "
            f"got {spec.port_b.component!r}"
        )
    return ResolvedJoint(
        joint_id=spec.joint_id,
        kind=_FAMILY_JOINT_KIND[contract.family],
        port_a=spec.port_a,
        port_b=spec.port_b,
        stiffness_pn_per_um=spec.stiffness_pn_per_um,
        rest_um=spec.rest_um,
        initial_state=spec.initial_state,
        semantic_edges=(contract.name,),
    )


class FilamentCrosslinkConnector:
    """One inter-component ``transient_actin`` edge over a two-array adjoint joint runtime.

    Delegates mechanics and the whole transaction API to
    :class:`~aleph.engine.load_path.LoadPathJointRuntime`; what this adds is the edge identity a connector
    slot binds under, the contract-derived guard, and the kinetics seam.
    """

    #: fidelity marker (a resolved per-joint central force, never a lumped body-to-body spring).
    per_joint_central_force = True
    aggregate_or_lumped = False

    def __init__(
        self,
        *,
        name: str,
        joint_runtime: LoadPathJointRuntime,
        kinetics: object | None = None,
        architecture: CellArchitecture | None = None,
    ) -> None:
        arch = architecture or reference_cell_architecture()
        contract = _contract(name, arch)
        self.name = contract.name
        # Endpoints come from the contract, never from the caller.
        self.component_a = contract.component_a
        self.component_b = contract.component_b
        self.bidirectional = contract.bidirectional
        self.adjoint_transfer_required = contract.adjoint_transfer_required
        self.chemistry_card = contract.chemistry_card
        self._joints = joint_runtime
        self._kinetics = kinetics
        self.ledger_calls: list[object] = []

    @property
    def n_joints(self) -> int:
        """Number of crosslink slots; fixed at build, so a bond count never changes the population."""
        return self._joints.capacity

    def accumulate(self) -> None:
        """Scatter the equal-and-opposite crosslink force into BOTH owners' force arrays."""
        self._joints.accumulate()

    def propose_events(self, *args: object, **kwargs: object) -> object:
        """Propose bind/unbind onto the joint runtime's candidate state via the injected rate law.

        Raises:
            RuntimeError: if no kinetics delegate is bound.  All three edges declare ``kinetics=True``;
                a transient crosslink with no rate law is a permanent weld, so this refuses rather than
                proposing nothing and looking like a running kinetic edge.
        """
        if self._kinetics is None:
            raise RuntimeError(
                f"{self.name!r} declares kinetics=True but no rate law is bound: its chemistry card is "
                f"{self.chemistry_card!r}.  Bind a delegate, or do not call propose_events — a transient "
                "crosslink without kinetics is a permanent weld."
            )
        return self._kinetics.propose_events(self._joints, *args, **kwargs)

    def snapshot_candidate(self) -> None:
        self._joints.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        self._joints.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self._joints.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Record the ledger handle (the adjoint-work term reduces device-side on the native lane).

        ADDED 2026-08-09 by the connector runtime census.  Without it this class was missing a hook
        :meth:`~aleph.engine.actor.CellActor.assert_fully_bound` REQUIRES of every connector, so all
        five edges it serves were counted as having a runtime while a ``require_complete=True`` world
        would have rejected them — the census's first catch, and exactly the class of "green for the
        wrong reason" it was written to find.
        """
        self.ledger_calls.append(ledger)


def build_filament_crosslink_connector(
    *,
    name: str,
    specs: tuple[FilamentCrosslinkJointSpec, ...],
    registry: ActorRegistry,
    actor_a: object,
    actor_b: object,
    kinetics: object | None = None,
    architecture: CellArchitecture | None = None,
) -> FilamentCrosslinkConnector:
    """Resolve the specs and bind one inter-component transient-actin edge (CUDA lane).

    Args:
        name: the declared edge, e.g. ``"sf_cortex_transient"``.
        specs: one :class:`FilamentCrosslinkJointSpec` per crosslink slot, all on the same edge.
        registry: the actor registry the ports are validated against.
        actor_a: endpoint A's segment-actor runtime.
        actor_b: endpoint B's segment-actor runtime or port view.
        kinetics: the bind/unbind rate law.  Optional here; ``propose_events`` refuses without it.
        architecture: optional composition (defaults to :func:`reference_cell_architecture`).

    Returns:
        A :class:`FilamentCrosslinkConnector` ready to bind at the ``name`` connector slot.

    Raises:
        ValueError: if any spec names a different edge, or the edge is not an inter-component
            ``transient_actin`` contract.
    """
    arch = architecture or reference_cell_architecture()
    off_edge = tuple(sorted({s.edge for s in specs} - {name}))
    if off_edge:
        raise ValueError(f"specs for {name!r} also name {off_edge}; one connector binds ONE edge")
    joints = tuple(resolve_filament_crosslink_joint(spec, registry, arch) for spec in specs)
    return FilamentCrosslinkConnector(
        name=name,
        joint_runtime=LoadPathJointRuntime(actor_a, actor_b, joints),
        kinetics=kinetics,
        architecture=arch,
    )
