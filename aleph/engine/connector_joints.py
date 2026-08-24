"""The remaining two-array connector joints: contact, cortical capture, and nascent adhesion series.

WHAT THIS CLOSES.  Eleven declared edges of
:func:`~aleph.engine.contracts.reference_cell_architecture` had no runtime.  They arrive here as
**two** force laws, not eleven and not five::

    UNILATERAL (may push, may never pull)
      membrane_cortex_contact         contact  membrane      <-> cortex     kinetics=False
      membrane_ecm_contact            contact  membrane      <-> ecm        kinetics=False
      nucleus_cortex_contact          contact  nucleus       <-> cortex     kinetics=False
      lamellipodium_membrane_contact  contact  lamellipodium <-> membrane   kinetics=True
      filopodium_membrane_tip         contact  filopodium    <-> membrane   kinetics=True

    BILATERAL (a bound crossbridge or clutch resists in both directions)
      mt_cortex_capture               motor    microtubule   <-> cortex     kinetics=True
      lamellipodium_nascent_fa + _nascent_clutch   ONE composite series spring, lamellipodium -> ecm
      filopodium_nascent_fa    + _nascent_clutch   ONE composite series spring, filopodium    -> ecm

THE FAMILY LABEL IS NOT THE PHYSICS BOUNDARY, IN BOTH DIRECTIONS.  ``contact`` is ONE declared family
and TWO force behaviours only in the KINETIC sense — all five edges are the same unilateral gap law;
what the two ``brownian_ratchet`` edges add is a rate law that decides when a barbed end is positioned
to push at all.  Conversely ``motor``, ``actin_anchor`` and ``fa_clutch`` are three families and ONE
mechanics: a state-gated bilateral Hookean central force, which
:func:`~aleph.engine.load_path._segment_joint_force_kernel` already is.  So this module adds exactly
ONE kernel — the unilateral one — and binds everything else to machinery that was already generic.

WHY UNILATERAL IS A DIFFERENT LAW AND NOT A PARAMETER.  A bilateral spring at the physiological gap
would let the membrane PULL the cortex back when they separate, and let a lamellipodial barbed end tow
the membrane inward.  Neither is a contact; a contact that pulls is an adhesion nobody declared.  The
sign is therefore structural, gated on the gap, and cannot be turned off by a stiffness value.

WHY THE NASCENT ADHESIONS ARE SERIES AND NOT TWO SPRINGS.  ``focal_adhesion`` declares
``owns_geometry=False``, so an actin anchor ending on it has ONE geometric endpoint and nothing to pull
against.  PI 2026-08-09 option (a) resolved this the way the mature adhesion already worked: both
semantic edges register under one ``mechanical_group``, the composite spring runs actin -> collagen, and
the geometry-less FA is the middle.  :meth:`~aleph.engine.actor.CellActor.bind_connector` is what
enforces one runtime per group; :func:`~aleph.engine.load_path.resolve_fa_series_group` is what proves
the two edges actually chain through a single shared component.

WHAT THE CALLER MUST STILL SUPPLY, AND WHY IT IS NOT DEFAULTED.  Every kinetic edge here declares
``kinetics=True`` and ``commit_on_accept=True``.  :meth:`ContractJointConnector.propose_events` RAISES
without an injected rate law rather than proposing nothing: a kinetic connector that silently proposes
nothing is indistinguishable from a permanent weld, which is the one thing the charter says a connector
may never be.  No rate constant, stiffness or contact distance is defined in this module — every
magnitude arrives on a spec, so nothing here can become an undeclared constant.

KNOWN GAP, recorded rather than invented.  All five contact edges and both nascent clutches declare
``remap_on_accept=True``: their pair set is meant to be re-discovered by a broad phase on each accepted
step.  This runtime holds a FIXED-CAPACITY pair list resolved at build, exactly as the mature
``alpha2beta1_collagen_series`` composite joint has since it landed.  That is a real limitation and it
is shared with the incumbent path rather than newly introduced here; closing it is a broad-phase
landing, not a connector landing, and inventing one now would put a contact-search heuristic in the
tree that no gate asked for.

Sanity Gate (before first execution; :mod:`aleph.tests.ac.engine.test_connector_joints`):
  * dimensions: stiffness pN/um, ``rest_um`` um, ``u`` a dimensionless barycentric in [0, 1].  For a
    contact joint ``rest_um`` is the CONTACT DISTANCE — the separation at which the force vanishes —
    because it packs into the same SoA slot the bilateral law reads as a rest length.
  * boundary cases: an empty joint list is refused by the joint runtime (an empty connector is a
    binding mistake, not a configuration); an edge outside the served set, an INTERNAL edge, or a
    series group that does not chain all raise at construction rather than binding wrong.
  * conservation invariant: both laws scatter ``+f`` through A's barycentric weights and ``-f``
    through B's, so the pair sums to zero exactly at any configuration.  That adjoint is what the
    balance gate tests, and it is why one kernel may be substituted for the other at all.
  * sign sense: unilateral force is non-zero only when ``length < rest``, and then points A away from
    B — separated bodies feel nothing and touching bodies are pushed apart, never drawn together.
    Bilateral joints keep the crosslink convention: a stretched joint pulls its endpoints together,
    and a joint below ``JointState.ACTIN_ENGAGED`` carries no force at all.
  * precision / measurement protocol: no host readback; this module produces no magnitude of its own.

CPU-importable: it allocates nothing and launches nothing itself, deferring every allocation and launch
to :class:`~aleph.engine.load_path.LoadPathJointRuntime`, so the structural gate drives it with
recording doubles.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorContract,
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
    resolve_fa_series_group,
)

__all__ = [
    "AdhesionSeriesJointSpec",
    "ConnectorJointSpec",
    "ContractJointConnector",
    "adhesion_series_groups",
    "build_adhesion_series_connector",
    "build_connector_joint",
    "connector_joint_edges",
    "contact_pair_reference",
    "force_kernel_for",
    "resolve_adhesion_series_joint",
    "resolve_connector_joint",
]

#: The NMII actuator owns every MOTOR edge that touches the ``nmii`` component (head-resolved
#: Stam-Hocky crossbridges).  A MOTOR edge that does NOT touch it is a different motor on a different
#: geometry — today only cortical dynein capturing a microtubule — and is served here.  Stated as a
#: structural rule rather than a name list so a second non-NMII motor is covered without an edit.
_NMII_COMPONENT = "nmii"

#: The pre-existing series, whose composite runtime is ``ecm_world.CompositeLoadPathECMClutchAdapter``.
#: Named here only so :func:`resolve_adhesion_series_joint` records the ``JointKind`` that is already in
#: device SoAs for it, rather than relabelling a landed population.
_MATURE_FA_SERIES_GROUP = "alpha2beta1_collagen_series"

_ACTIVE_STATE_MIN = wp.constant(int(JointState.ACTIN_ENGAGED))


@wp.kernel
def _unilateral_contact_force_kernel(
    pos_a: wp.array(dtype=wp.vec3d),
    force_a: wp.array(dtype=wp.vec3d),
    segments_a: wp.array(dtype=wp.int32, ndim=2),
    pos_b: wp.array(dtype=wp.vec3d),
    force_b: wp.array(dtype=wp.vec3d),
    segments_b: wp.array(dtype=wp.int32, ndim=2),
    active: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    element_a: wp.array(dtype=wp.int32),
    element_b: wp.array(dtype=wp.int32),
    u_a: wp.array(dtype=wp.float64),
    u_b: wp.array(dtype=wp.float64),
    stiffness: wp.array(dtype=wp.float64),
    rest: wp.array(dtype=wp.float64),
    load: wp.array(dtype=wp.float64),
    energy: wp.array(dtype=wp.float64),
    joint_force_on_a: wp.array(dtype=wp.vec3d),
) -> None:
    """One thread per contact pair; force-free at and beyond ``rest``, repulsive inside it.

    Argument list and adjoint scatter are identical to
    :func:`~aleph.engine.load_path._segment_joint_force_kernel` — the same barycentric weights gather
    and scatter, so the pair force sums to zero exactly.  What differs is the single ``gap < 0`` branch:
    ``rest`` is read as a CONTACT DISTANCE, and a separated pair contributes nothing at all.
    """
    t = wp.tid()
    if active[t] == 0 or state[t] < _ACTIVE_STATE_MIN:
        load[t] = wp.float64(0.0)
        energy[t] = wp.float64(0.0)
        joint_force_on_a[t] = wp.vec3d(0.0, 0.0, 0.0)
        return
    ea = element_a[t]
    eb = element_b[t]
    ia0 = segments_a[ea, 0]
    ia1 = segments_a[ea, 1]
    ib0 = segments_b[eb, 0]
    ib1 = segments_b[eb, 1]
    wa1 = u_a[t]
    wa0 = wp.float64(1.0) - wa1
    wb1 = u_b[t]
    wb0 = wp.float64(1.0) - wb1
    xa = wa0 * pos_a[ia0] + wa1 * pos_a[ia1]
    xb = wb0 * pos_b[ib0] + wb1 * pos_b[ib1]
    delta = xb - xa
    length = wp.length(delta)
    gap = length - rest[t]
    pair_force = wp.vec3d(0.0, 0.0, 0.0)
    overlap = wp.float64(0.0)
    if gap < wp.float64(0.0) and length > wp.float64(1.0e-15):
        # gap < 0, so this vector points A AWAY from B and its negation pushes B away from A.
        pair_force = stiffness[t] * gap * delta / length
        overlap = -gap
    wp.atomic_add(force_a, ia0, wa0 * pair_force)
    wp.atomic_add(force_a, ia1, wa1 * pair_force)
    wp.atomic_add(force_b, ib0, -wb0 * pair_force)
    wp.atomic_add(force_b, ib1, -wb1 * pair_force)
    joint_force_on_a[t] = pair_force
    load[t] = stiffness[t] * overlap
    energy[t] = wp.float64(0.5) * stiffness[t] * overlap * overlap


def contact_pair_reference(
    separation_um: float,
    contact_distance_um: float,
    stiffness_pn_per_um: float,
) -> float:
    """Host oracle for the unilateral law: the SIGNED pair magnitude the kernel computes.

    Negative means the pair force vector points endpoint A away from endpoint B — a repulsion, the only
    sign a contact may ever have.  Zero at and beyond the contact distance.  This mirrors
    :func:`~aleph.engine.load_path.segment_pair_reference` for the bilateral law and exists for the same
    reason: the sign convention is the part of a contact that can be silently inverted, and a shape
    assertion downstream would not notice a connector that pulls.

    Args:
        separation_um: distance between the two interpolated endpoints [um].
        contact_distance_um: the gap at and beyond which the pair is force-free [um].
        stiffness_pn_per_um: contact stiffness [pN/um].

    Returns:
        ``k * (separation - contact_distance)`` when overlapping, else ``0.0``.
    """
    gap = float(separation_um) - float(contact_distance_um)
    return float(stiffness_pn_per_um) * gap if gap < 0.0 else 0.0


def _contract(architecture: CellArchitecture, name: str) -> ConnectorContract:
    contract = next((c for c in architecture.connectors if c.name == name), None)
    if contract is None:
        raise ValueError(f"connector {name!r} is not declared in the architecture")
    return contract


def _serves(contract: ConnectorContract) -> bool:
    """Whether this runtime is the one that binds ``contract``, decided structurally."""
    if contract.scope is not ConnectorScope.INTER_COMPONENT:
        return False
    if contract.family is ConnectorFamily.CONTACT:
        return True
    if contract.family is ConnectorFamily.MOTOR:
        return _NMII_COMPONENT not in (contract.component_a, contract.component_b)
    return False


def connector_joint_edges(architecture: CellArchitecture | None = None) -> tuple[str, ...]:
    """Return every SINGLE-edge connector this runtime serves, in architecture order.

    Read from the contract, so a new component declaring a contact — or a second non-NMII motor — is
    covered without editing this module.  The adhesion SERIES edges are excluded: they bind two at a
    time through :func:`build_adhesion_series_connector`, and :func:`adhesion_series_groups` lists them.
    """
    arch = architecture or reference_cell_architecture()
    return tuple(c.name for c in arch.connectors if _serves(c))


def adhesion_series_groups(architecture: CellArchitecture | None = None) -> tuple[str, ...]:
    """Return every declared adhesion-series ``mechanical_group`` name, in architecture order.

    Includes the mature ``alpha2beta1_collagen_series``: it is the same construction, and listing it
    from the contract is what keeps a fourth series from being invisible to this module.
    """
    arch = architecture or reference_cell_architecture()
    seen: list[str] = []
    for connector in arch.connectors:
        group = connector.mechanical_group
        if group is not None and group not in seen:
            seen.append(group)
    return tuple(seen)


def force_kernel_for(contract: ConnectorContract) -> object | None:
    """Return the pair law an edge is launched with, decided by its FAMILY and never by the caller.

    ``None`` means the bilateral default,
    :func:`~aleph.engine.load_path._segment_joint_force_kernel`.  Exposed rather than inlined into
    :func:`build_connector_joint` so the choice can be gated without allocating device memory: a
    contact edge silently getting the bilateral law is a connector that PULLS, and that is a wrong
    answer no shape assertion downstream would notice.
    """
    return _unilateral_contact_force_kernel if contract.family is ConnectorFamily.CONTACT else None


def _joint_kind(contract: ConnectorContract) -> JointKind:
    """Return which biology the SoA records for a single-edge joint, for diagnostics only."""
    if contract.family is ConnectorFamily.CONTACT:
        # The two kinetic contacts are the Brownian ratchet; the three non-kinetic ones are steric.
        # Read from `kinetics` rather than from the chemistry card so a new card cannot mislabel it.
        return (
            JointKind.BROWNIAN_RATCHET_CONTACT
            if contract.kinetics
            else JointKind.STERIC_CONTACT
        )
    return JointKind.CORTICAL_DYNEIN_CAPTURE


@dataclass(frozen=True, slots=True)
class ConnectorJointSpec:
    """One pair slot on a single-edge two-array connector.

    Attributes:
        joint_id: unique within the connector.
        edge: the declared connector name this slot belongs to.
        port_a: endpoint on the contract's ``component_a``.
        port_b: endpoint on the contract's ``component_b``.
        stiffness_pn_per_um: pair stiffness [pN/um].
        rest_um: zero-force separation [um].  On a CONTACT edge this is the contact distance — the
            physiological gap at and beyond which the pair is force-free.
        initial_state: the state the slot starts in; below ``ACTIN_ENGAGED`` it carries no force.
    """

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


@dataclass(frozen=True, slots=True)
class AdhesionSeriesJointSpec:
    """One composite adhesion: an actin material point sprung to a collagen ligand.

    There is no FA-side port because ``focal_adhesion`` owns no geometry — that is the whole reason the
    two semantic edges are one mechanical joint.  ``group`` names the ``mechanical_group`` the two edges
    declare, and the two endpoint COMPONENTS are derived from it rather than supplied.
    """

    joint_id: int
    group: str
    fa_cluster_id: int
    actin_port: PortRef
    ligand_port: PortRef
    stiffness_pn_per_um: float
    rest_um: float
    initial_state: JointState = JointState.ACTIN_ENGAGED

    def __post_init__(self) -> None:
        if self.joint_id < 0 or self.fa_cluster_id < 0:
            raise ValueError("joint and FA-cluster IDs must be nonnegative")
        if not np.isfinite(self.stiffness_pn_per_um) or self.stiffness_pn_per_um <= 0.0:
            raise ValueError("joint stiffness must be finite and positive")
        if not np.isfinite(self.rest_um) or self.rest_um < 0.0:
            raise ValueError("joint rest length must be finite and nonnegative")


def resolve_connector_joint(
    spec: ConnectorJointSpec,
    registry: ActorRegistry,
    architecture: CellArchitecture,
) -> ResolvedJoint:
    """Validate one single-edge pair against its declared contract and resolve it.

    Endpoint COMPONENTS are checked against the contract, which is what can be checked today.  Endpoint
    ROLES are carried through unasserted: :class:`~aleph.engine.load_path.EndpointRole` declares members
    for stress-fibre and collagen endpoints only, and adding one for a membrane quadrature or a cortical
    dynein site is a contract change and therefore a PI decision.

    Raises:
        ValueError: if the edge is not one this runtime serves, or either port belongs to a component
            the contract does not name.
    """
    contract = _contract(architecture, spec.edge)
    if not _serves(contract):
        raise ValueError(
            f"connector {spec.edge!r} (family {contract.family.value!r}) is not served by this "
            "runtime; it binds inter-component CONTACT edges and non-NMII MOTOR edges"
        )
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
        kind=_joint_kind(contract),
        port_a=spec.port_a,
        port_b=spec.port_b,
        stiffness_pn_per_um=spec.stiffness_pn_per_um,
        rest_um=spec.rest_um,
        initial_state=spec.initial_state,
        semantic_edges=(contract.name,),
    )


def resolve_adhesion_series_joint(
    spec: AdhesionSeriesJointSpec,
    registry: ActorRegistry,
    architecture: CellArchitecture,
) -> ResolvedJoint:
    """Validate and create exactly ONE mechanical spring for both of a series' semantic edges.

    Raises:
        RuntimeError: if the group resolves to more than one mechanical joint, which would count the
            adhesion's stiffness twice.
        ValueError: if either port belongs to a component the resolved series does not end on.
    """
    resolution = resolve_fa_series_group(architecture, spec.group)
    registry.validate_port(spec.actin_port)
    registry.validate_port(spec.ligand_port)
    if spec.actin_port.component != resolution.actin_component:
        raise ValueError(
            f"{spec.group!r} actin-side port must belong to {resolution.actin_component!r}, "
            f"got {spec.actin_port.component!r}"
        )
    if spec.ligand_port.component != resolution.ligand_component:
        raise ValueError(
            f"{spec.group!r} ligand-side port must belong to {resolution.ligand_component!r}, "
            f"got {spec.ligand_port.component!r}"
        )
    if resolution.mechanical_joint_count != 1:
        raise RuntimeError(
            f"{spec.group!r} semantic edges would double-count mechanical stiffness"
        )
    # The MATURE series predates this resolver and its kind is already written into device SoAs as 0;
    # only the nascent series are new.  Keyed on the declared group name, not on the actin component,
    # so a fourth series on `sf_arc` would still be recorded as nascent rather than silently mature.
    kind = (
        JointKind.FA_COMPOSITE_SERIES
        if spec.group == _MATURE_FA_SERIES_GROUP
        else JointKind.NASCENT_FA_SERIES
    )
    return ResolvedJoint(
        joint_id=spec.joint_id,
        kind=kind,
        port_a=spec.actin_port,
        port_b=spec.ligand_port,
        stiffness_pn_per_um=spec.stiffness_pn_per_um,
        rest_um=spec.rest_um,
        initial_state=spec.initial_state,
        semantic_edges=resolution.semantic_edges,
        fa_cluster_id=spec.fa_cluster_id,
    )


class ContractJointConnector:
    """One declared connector — or one adhesion series — over a two-array adjoint joint runtime.

    Mechanics and the whole transaction API are delegated to
    :class:`~aleph.engine.load_path.LoadPathJointRuntime`.  What this adds is the identity a connector
    slot binds under, read from the contract so a caller cannot mis-state it, and the kinetics seam.

    The several ``accumulate_*`` aliases exist because the owning facades declare different mechanics
    protocols for the same call — ``accumulate_contact`` on ``ECMWorld``, ``accumulate_actor`` on
    ``ProtrusionActors``, ``accumulate_motor`` on ``MicrotubuleRig``.  All of them dispatch one force
    launch and differ only in which non-owning view the facade hands in for signature parity; the live
    arrays the joint scatters into were borrowed at build.
    """

    #: fidelity marker (a resolved per-pair central force, never a lumped body-to-body spring).
    per_joint_central_force = True
    aggregate_or_lumped = False

    def __init__(
        self,
        *,
        name: str,
        component_a: str,
        component_b: str,
        joint_runtime: LoadPathJointRuntime,
        kinetics_declared: bool,
        chemistry_card: str | None,
        bidirectional: bool,
        adjoint_transfer_required: bool,
        mechanical_group: str | None = None,
        kinetics: object | None = None,
    ) -> None:
        self.name = name
        self.component_a = component_a
        self.component_b = component_b
        self.mechanical_group = mechanical_group
        self.chemistry_card = chemistry_card
        self.bidirectional = bidirectional
        self.adjoint_transfer_required = adjoint_transfer_required
        self._kinetics_declared = kinetics_declared
        self._joints = joint_runtime
        self._kinetics = kinetics
        self.mechanics_calls: list[object] = []
        self.ledger_calls: list[object] = []

    @property
    def n_joints(self) -> int:
        """Number of pair slots; fixed at build (see the module docstring's recorded broad-phase gap)."""
        return self._joints.capacity

    def accumulate(self) -> None:
        """Scatter the equal-and-opposite pair force into BOTH owners' force arrays."""
        self._joints.accumulate()

    def accumulate_contact(self, *views: object) -> None:
        """``ECMWorld`` mechanics-protocol alias for :meth:`accumulate`."""
        self.mechanics_calls.append(views)
        self._joints.accumulate()

    def accumulate_actor(self, *views: object) -> None:
        """``ProtrusionActors`` mechanics-protocol alias for :meth:`accumulate`."""
        self.mechanics_calls.append(views)
        self._joints.accumulate()

    def accumulate_motor(self, *views: object) -> None:
        """``MicrotubuleRig`` mechanics-protocol alias for :meth:`accumulate`."""
        self.mechanics_calls.append(views)
        self._joints.accumulate()

    def propose_events(self, *args: object, **kwargs: object) -> object:
        """Propose bind/unbind onto the joint runtime's candidate state via the injected rate law.

        Raises:
            RuntimeError: if the contract declares ``kinetics=False`` — a non-kinetic connector has no
                events to propose, and answering the call would invent a transition the contract
                forbids — or if it declares ``kinetics=True`` and no rate law is bound, because a
                kinetic connector that silently proposes nothing is a permanent weld.
        """
        if not self._kinetics_declared:
            raise RuntimeError(
                f"{self.name!r} declares kinetics=False: it has no events.  Its pair set changes only "
                "through a broad-phase remap, which is not a kinetic proposal."
            )
        if self._kinetics is None:
            raise RuntimeError(
                f"{self.name!r} declares kinetics=True but no rate law is bound: its chemistry card is "
                f"{self.chemistry_card!r}.  Bind a delegate, or do not call propose_events — a kinetic "
                "connector without a rate law is a permanent weld."
            )
        return self._kinetics.propose_events(self._joints, *args, **kwargs)

    def snapshot_candidate(self) -> None:
        self._joints.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        self._joints.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self._joints.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Record the ledger handle (the adjoint-work term reduces device-side on the native lane)."""
        self.ledger_calls.append(ledger)


def build_connector_joint(
    *,
    name: str,
    specs: tuple[ConnectorJointSpec, ...],
    registry: ActorRegistry,
    actor_a: object,
    actor_b: object,
    kinetics: object | None = None,
    architecture: CellArchitecture | None = None,
) -> ContractJointConnector:
    """Resolve the specs and bind ONE declared single-edge connector (CUDA lane).

    The force law is chosen from the contract's FAMILY, never from the caller: a ``CONTACT`` edge gets
    the unilateral kernel and everything else gets the bilateral one.

    Args:
        name: the declared edge, e.g. ``"membrane_cortex_contact"`` or ``"mt_cortex_capture"``.
        specs: one :class:`ConnectorJointSpec` per pair slot, all naming ``name``.
        registry: the actor registry the ports are validated against.
        actor_a: endpoint A's segment-actor runtime or borrowed view.
        actor_b: endpoint B's segment-actor runtime or borrowed view.
        kinetics: the rate law, required only before ``propose_events`` is called.
        architecture: optional composition (defaults to :func:`reference_cell_architecture`).

    Returns:
        A :class:`ContractJointConnector` ready to bind at the ``name`` connector slot.

    Raises:
        ValueError: if any spec names a different edge, or ``name`` is not an edge this runtime serves.
    """
    arch = architecture or reference_cell_architecture()
    off_edge = tuple(sorted({s.edge for s in specs} - {name}))
    if off_edge:
        raise ValueError(f"specs for {name!r} also name {off_edge}; one connector binds ONE edge")
    contract = _contract(arch, name)
    joints = tuple(resolve_connector_joint(spec, registry, arch) for spec in specs)
    kernel = force_kernel_for(contract)
    return ContractJointConnector(
        name=contract.name,
        component_a=contract.component_a,
        component_b=contract.component_b,
        joint_runtime=LoadPathJointRuntime(actor_a, actor_b, joints, force_kernel=kernel),
        kinetics_declared=contract.kinetics,
        chemistry_card=contract.chemistry_card,
        bidirectional=contract.bidirectional,
        adjoint_transfer_required=contract.adjoint_transfer_required,
        kinetics=kinetics,
    )


def build_adhesion_series_connector(
    *,
    group: str,
    specs: tuple[AdhesionSeriesJointSpec, ...],
    registry: ActorRegistry,
    actin_actor: object,
    ligand_actor: object,
    kinetics: object | None = None,
    architecture: CellArchitecture | None = None,
) -> ContractJointConnector:
    """Resolve the specs and build the ONE composite runtime both of a series' edges bind to.

    The returned object must be bound under BOTH semantic edge names — that is what makes the two edges
    one spring rather than two in parallel, and
    :meth:`~aleph.engine.actor.CellActor.bind_connector` refuses any other pairing.

    Args:
        group: the declared ``mechanical_group``, e.g. ``"lamellipodium_nascent_series"``.
        specs: one :class:`AdhesionSeriesJointSpec` per adhesion, all naming ``group``.
        registry: the actor registry the ports are validated against.
        actin_actor: the actin-side component's segment-actor runtime or borrowed view.
        ligand_actor: the collagen-side runtime or borrowed view.
        kinetics: the clutch bind/unbind rate law, required only before ``propose_events`` is called.
        architecture: optional composition (defaults to :func:`reference_cell_architecture`).

    Returns:
        A :class:`ContractJointConnector` carrying ``mechanical_group == group``.

    Raises:
        ValueError: if any spec names a different group, or the group does not resolve to a chained
            one-anchor/one-clutch series.
    """
    arch = architecture or reference_cell_architecture()
    off_group = tuple(sorted({s.group for s in specs} - {group}))
    if off_group:
        raise ValueError(
            f"specs for {group!r} also name {off_group}; one composite runtime binds ONE series"
        )
    resolution = resolve_fa_series_group(arch, group)
    joints = tuple(resolve_adhesion_series_joint(spec, registry, arch) for spec in specs)
    edges = {c.name: c for c in arch.mechanical_group(group)}
    kinetic = any(c.kinetics for c in edges.values())
    cards = " + ".join(resolution.chemistry_cards)
    return ContractJointConnector(
        name=group,
        component_a=resolution.actin_component,
        component_b=resolution.ligand_component,
        joint_runtime=LoadPathJointRuntime(actin_actor, ligand_actor, joints),
        kinetics_declared=kinetic,
        chemistry_card=cards,
        bidirectional=all(c.bidirectional for c in edges.values()),
        adjoint_transfer_required=all(c.adjoint_transfer_required for c in edges.values()),
        mechanical_group=group,
        kinetics=kinetics,
    )
