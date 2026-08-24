"""Runtime contract and minimal mechanics slice for SF--FA--ECM load paths.

The top-level :mod:`aleph.engine.contracts` graph declares that ``sf_arc`` reaches
``ecm`` through the geometry-less ``focal_adhesion`` component.  This module resolves
those two semantic edges to one composite series mechanical joint.  It also supplies
the typed internal dorsal-SF--transverse-arc connector that cannot be represented by
co-location in an actin array.

The device runtime is intentionally small: actors own separate position/force arrays,
while a fixed-capacity structure-of-arrays (SoA) owns interpolated segment ports, joint
state, loads, energy, and accepted-step kinetics state.  It is the VS-0 seam, not a
complete FA maturation or partner-search implementation.

Units are micrometres, piconewtons, and seconds.

Sanity Gate:
    * topology: only ventral ends or a dorsal basal end may enter an FA series joint;
      a dorsal free end may connect only to a transverse-arc material point;
    * generation: stale actor/entity generations fail before device upload;
    * double-count: the ACTIN_ANCHOR and FA_CLUTCH semantic edges resolve to exactly
      one spring energy and one force scatter;
    * conservation: the central pair force is equal-and-opposite, moment-free, and its
      interpolation/scatter are adjoints, so nodal virtual work equals ``-dU``;
    * transaction: committed state, age, and RNG epoch change only when the scheduler's
      device acceptance scalar is non-zero.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorFamily,
    ConnectorScope,
)

__all__ = [
    "ActorRecord",
    "ActorRegistry",
    "ActorSegmentRuntime",
    "BorrowedSegmentActorView",
    "CompositeFAJointSpec",
    "DorsalArcJointSpec",
    "ElementKind",
    "EndpointRole",
    "JointKind",
    "JointState",
    "LoadPathJointRuntime",
    "MechanicalGroupResolution",
    "PortRef",
    "ResolvedJoint",
    "SegmentPairReference",
    "resolve_dorsal_arc_joint",
    "resolve_fa_series_group",
    "resolve_fa_series_joint",
    "segment_pair_reference",
]

_FA_SERIES_GROUP = "alpha2beta1_collagen_series"
_FA_SEMANTIC_EDGES = frozenset({"fa_actin_anchor", "integrin_collagen_clutch"})
_FA_ACTIN_ROLES = frozenset(
    {
        "ventral_end_0",
        "ventral_end_1",
        "dorsal_basal_end",
    }
)


class ElementKind(IntEnum):
    """Geometry element addressed by a connector port."""

    NODE = 0
    SEGMENT = 1
    TRIANGLE = 2
    MATERIAL_POINT = 3


class EndpointRole(IntEnum):
    """Biological role of one load-path endpoint."""

    VENTRAL_END_0 = 0
    VENTRAL_END_1 = 1
    DORSAL_BASAL_END = 2
    DORSAL_FREE_END = 3
    TRANSVERSE_ARC_MATERIAL = 4
    COLLAGEN_LIGAND = 5

    @property
    def contract_name(self) -> str:
        """Stable snake-case name used by endpoint-role validation."""
        return {
            EndpointRole.VENTRAL_END_0: "ventral_end_0",
            EndpointRole.VENTRAL_END_1: "ventral_end_1",
            EndpointRole.DORSAL_BASAL_END: "dorsal_basal_end",
            EndpointRole.DORSAL_FREE_END: "dorsal_free_end",
            EndpointRole.TRANSVERSE_ARC_MATERIAL: "transverse_arc_material",
            EndpointRole.COLLAGEN_LIGAND: "collagen_ligand",
        }[self]


class JointKind(IntEnum):
    """Mechanical joint kinds in the first load-path slice."""

    FA_COMPOSITE_SERIES = 0
    DORSAL_ARC_INTERNAL = 1
    # Two-array adjoint filament crosslinks.  One force law, three biologies — the kind is what the SoA
    # carries so a diagnostic can tell them apart, not a branch in the kernel.
    # Additive: the value is written into device SoA, so members are appended, never renumbered.
    #: `sf_arc`/`lamellipodium`/`filopodium` -> `cortex` transient actin crosslink.
    TRANSIENT_ACTIN_CROSSLINK = 2
    #: `intermediate_filament` -> `sf_arc`, plectin.
    PLECTIN_IF_ACTIN = 3
    #: `microtubule` -> `sf_arc`, spectraplakin.
    SPECTRAPLAKIN_MT_ACTIN = 4
    #: Unilateral steric non-penetration; force-free at and beyond the physiological gap.
    #: `membrane_cortex_contact` / `membrane_ecm_contact` / `nucleus_cortex_contact`.
    STERIC_CONTACT = 5
    #: Unilateral polymerisation ratchet: a growing barbed end may PUSH a membrane, never pull it.
    #: `lamellipodium_membrane_contact` / `filopodium_membrane_tip`.
    BROWNIAN_RATCHET_CONTACT = 6
    #: A bound cortical-dynein crossbridge on a microtubule lattice site (`mt_cortex_capture`).
    CORTICAL_DYNEIN_CAPTURE = 7
    #: A nascent adhesion's composite protrusion-actin -> collagen series spring
    #: (`{lamellipodium,filopodium}_nascent_fa` + `_nascent_clutch`).
    NASCENT_FA_SERIES = 8


class JointState(IntEnum):
    """Common committed joint states; only engaged/reinforced states carry force."""

    FREE = 0
    LIGAND_BOUND = 1
    ACTIN_ENGAGED = 2
    REINFORCED = 3


_ACTIVE_STATE_MIN = wp.constant(int(JointState.ACTIN_ENGAGED))


def _device_storage_key(array: object) -> tuple[object, ...]:
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


@dataclass(frozen=True, slots=True)
class BorrowedSegmentActorView:
    """Non-owning segment arrays used to build a joint over real component state owners."""

    record: ActorRecord
    pos_d: wp.array
    force_d: wp.array
    segments_d: wp.array

    def __post_init__(self) -> None:
        arrays = (self.pos_d, self.force_d, self.segments_d)
        if not all(bool(getattr(getattr(array, "device", None), "is_cuda", False)) for array in arrays):
            raise ValueError("borrowed load-path actor arrays must reside on CUDA")
        if getattr(self.pos_d, "dtype", None) != wp.vec3d:
            raise TypeError("borrowed actor position must have dtype wp.vec3d")
        if getattr(self.force_d, "dtype", None) != wp.vec3d:
            raise TypeError("borrowed actor force must have dtype wp.vec3d")
        if getattr(self.segments_d, "dtype", None) != wp.int32:
            raise TypeError("borrowed actor segments must have dtype wp.int32")
        if self.pos_d.shape != self.force_d.shape or len(self.pos_d.shape) != 1:
            raise ValueError("borrowed actor position and force must have matching rank-1 shapes")
        if len(self.segments_d.shape) != 2 or int(self.segments_d.shape[1]) != 2:
            raise ValueError("borrowed actor segments must have shape (n_segments, 2)")
        if int(self.segments_d.shape[0]) != self.record.n_elements:
            raise ValueError("borrowed actor record does not match segment count")
        if len({str(array.device) for array in arrays}) != 1:
            raise ValueError("borrowed actor arrays must share one CUDA device")
        if _device_storage_key(self.pos_d) == _device_storage_key(self.force_d):
            raise ValueError("borrowed actor position and force must not alias")

    @property
    def device(self) -> str:
        """CUDA device shared with the owning component."""
        return str(self.pos_d.device)


@dataclass(frozen=True, slots=True)
class ActorRecord:
    """Build-time identity and generation of one runtime component actor."""

    component: str
    actor_id: int
    actor_generation: int
    entity_id: int
    entity_generation: int
    n_elements: int

    def __post_init__(self) -> None:
        if not self.component.strip():
            raise ValueError("actor component must be non-empty")
        if min(self.actor_id, self.actor_generation, self.entity_id, self.entity_generation) < 0:
            raise ValueError("actor/entity IDs and generations must be nonnegative")
        if self.n_elements <= 0:
            raise ValueError("actor must expose at least one geometry element")


@dataclass(frozen=True, slots=True)
class ActorRegistry:
    """Immutable build-time registry used to reject stale or cross-component ports."""

    actors: tuple[ActorRecord, ...]

    def __post_init__(self) -> None:
        ids = tuple(actor.actor_id for actor in self.actors)
        if len(ids) != len(set(ids)):
            raise ValueError("actor IDs must be unique")

    def actor(self, actor_id: int) -> ActorRecord:
        """Return an actor by stable ID."""
        for actor in self.actors:
            if actor.actor_id == actor_id:
                return actor
        raise KeyError(actor_id)

    def validate_port(self, port: PortRef) -> ActorRecord:
        """Validate ownership, generation, entity identity, and element bounds."""
        actor = self.actor(port.actor_id)
        if actor.component != port.component:
            raise ValueError(
                f"port component {port.component!r} does not own actor {port.actor_id} "
                f"({actor.component!r})"
            )
        if actor.actor_generation != port.actor_generation:
            raise ValueError(f"stale actor generation for actor {port.actor_id}")
        if actor.entity_id != port.entity_id or actor.entity_generation != port.entity_generation:
            raise ValueError(f"stale entity generation for actor {port.actor_id}")
        if not 0 <= port.element_id < actor.n_elements:
            raise ValueError(f"element {port.element_id} is out of bounds for actor {port.actor_id}")
        return actor


@dataclass(frozen=True, slots=True)
class PortRef:
    """Stable connector endpoint on an actor geometry element.

    ``local_coordinates[0]`` is the segment barycentric coordinate ``u`` in VS-0.
    The four-wide payload leaves room for triangle barycentrics without changing the
    SoA schema.
    """

    component: str
    actor_id: int
    actor_generation: int
    entity_id: int
    entity_generation: int
    element_kind: ElementKind
    element_id: int
    local_coordinates: tuple[float, float, float, float]
    role: EndpointRole

    def __post_init__(self) -> None:
        if not self.component.strip():
            raise ValueError("port component must be non-empty")
        if self.element_id < 0:
            raise ValueError("port element ID must be nonnegative")
        if len(self.local_coordinates) != 4 or not np.all(np.isfinite(self.local_coordinates)):
            raise ValueError("port local coordinates must be four finite values")
        if self.element_kind is ElementKind.SEGMENT:
            u = self.local_coordinates[0]
            if not 0.0 <= u <= 1.0:
                raise ValueError("segment barycentric coordinate must be in [0, 1]")

    @property
    def segment_u(self) -> float:
        """Return the segment interpolation coordinate or reject a non-segment port."""
        if self.element_kind is not ElementKind.SEGMENT:
            raise ValueError("VS-0 mechanics requires a SEGMENT port")
        return float(self.local_coordinates[0])


@dataclass(frozen=True, slots=True)
class MechanicalGroupResolution:
    """Two FA semantic graph edges resolved to one mechanical contribution."""

    name: str
    semantic_edges: tuple[str, str]
    chemistry_cards: tuple[str, str]
    #: How many mechanical joints these semantic edges actually resolve to.  Carries NO default on
    #: purpose: it defaulted to 1 while :func:`resolve_fa_series_group` never set it, so the
    #: ``!= 1`` double-count guard downstream could not fire under any input — a guard that was
    #: structurally unable to fail.  It must be measured from the declared connectors.
    mechanical_joint_count: int
    #: The geometry-less component both edges pass through (``focal_adhesion`` in every declared series).
    middle_component: str = ""
    #: The ACTIN_ANCHOR edge's far end — the actin population the composite spring starts on.
    actin_component: str = ""
    #: The FA_CLUTCH edge's far end — the ligand population the composite spring ends on.
    ligand_component: str = ""


def resolve_fa_series_group(
    architecture: CellArchitecture,
    group_name: str = _FA_SERIES_GROUP,
) -> MechanicalGroupResolution:
    """Resolve one geometry-less adhesion series and enforce its double-stiffness guard.

    GENERALISED 2026-08-09.  This used to assert the group's two edge NAMES and two chemistry CARDS
    against module constants, so it described the mature ``alpha2beta1_collagen_series`` and nothing
    else — while PI option (a) had just created two more series of the same shape
    (``lamellipodium_nascent_series``, ``filopodium_nascent_series``).  A name list is also the weaker
    check: it cannot tell whether the two edges actually CHAIN.  What makes a series a series is
    structural and is now what is verified — one ACTIN_ANCHOR edge, one FA_CLUTCH edge, both
    inter-component, both carrying a chemistry card, and sharing EXACTLY ONE component, which is the
    geometry-less middle the composite spring passes through.  Two edges that shared no component would
    be two unrelated springs filed under one group; two that shared both would be one edge written twice.
    """
    connectors = architecture.mechanical_group(group_name)
    names = frozenset(connector.name for connector in connectors)
    if len(connectors) != 2:
        raise ValueError(f"adhesion series {group_name!r} must contain exactly 2 edges; got {names}")
    families = frozenset(connector.family for connector in connectors)
    if families != frozenset({ConnectorFamily.ACTIN_ANCHOR, ConnectorFamily.FA_CLUTCH}):
        raise ValueError("FA series group needs exactly one ACTIN_ANCHOR and one FA_CLUTCH edge")
    if any(connector.scope is not ConnectorScope.INTER_COMPONENT for connector in connectors):
        raise ValueError("FA series semantic edges must be inter-component")
    cards = tuple(sorted(connector.chemistry_card or "" for connector in connectors))
    if not all(cards):
        raise ValueError(f"adhesion series {group_name!r} has an edge with no chemistry card")
    first, second = connectors
    middle = {first.component_a, first.component_b} & {second.component_a, second.component_b}
    if len(middle) != 1:
        raise ValueError(
            f"adhesion series {group_name!r} edges must share exactly one (geometry-less) component; "
            f"shared={sorted(middle)}"
        )
    # MEASURE the joint count rather than assuming it. The group itself is ONE composite joint; any
    # OTHER connector spanning the same component pair is a second independent load path across that
    # pair — the stiffness counted twice, which is exactly the condition the downstream `!= 1` guard
    # exists to catch and could never previously reach.
    #
    # Counting by distinct mechanical_group value would NOT work: `architecture.mechanical_group()`
    # selects by that very name, so every connector it returns shares it and the count is 1 by
    # construction — a second dead guard. The component PAIR is the observable that can actually
    # differ. Note `lamellipodium_nascent_fa` / `filopodium_nascent_fa` also touch focal_adhesion
    # but span different pairs (their own actin populations), so they are correctly not rivals.
    group_pairs = {
        frozenset((connector.component_a, connector.component_b)) for connector in connectors
    }
    #
    # REFINED 2026-08-09, when option (a) added two more FA series.  Spanning the same pair is not by
    # itself a double count: `focal_adhesion` is ONE component holding every adhesion cluster, so a
    # nascent adhesion's clutch spans `focal_adhesion <-> ecm` exactly as the mature one does while
    # being a different physical adhesion with its own composite spring.  What makes a connector a
    # RIVAL is spanning the pair WITHOUT a mechanical group of its own — an ungrouped second path
    # across the pair is an independent stiffness in parallel, which is the condition this guard was
    # written for.  A connector in its own declared group is a separate composite joint by
    # construction, and `bind_connector` already forces each such group to one runtime.
    rivals = tuple(
        connector
        for connector in architecture.connectors
        if connector.name not in names
        and frozenset((connector.component_a, connector.component_b)) in group_pairs
        and connector.mechanical_group is None
    )
    joint_count = 1 + len(rivals)
    hub = next(iter(middle))
    anchor = next(c for c in connectors if c.family is ConnectorFamily.ACTIN_ANCHOR)
    clutch = next(c for c in connectors if c.family is ConnectorFamily.FA_CLUTCH)
    return MechanicalGroupResolution(
        name=group_name,
        semantic_edges=tuple(sorted(names)),  # type: ignore[arg-type]
        chemistry_cards=cards,  # type: ignore[arg-type]
        mechanical_joint_count=joint_count,
        middle_component=hub,
        actin_component=anchor.component_b if anchor.component_a == hub else anchor.component_a,
        ligand_component=clutch.component_b if clutch.component_a == hub else clutch.component_a,
    )


@dataclass(frozen=True, slots=True)
class CompositeFAJointSpec:
    """One mechanical SF--ECM spring carrying both geometry-less FA semantic edges."""

    joint_id: int
    fa_cluster_id: int
    sf_port: PortRef
    ecm_port: PortRef
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


@dataclass(frozen=True, slots=True)
class DorsalArcJointSpec:
    """One dynamic internal connector from a dorsal free end to transverse-arc actin."""

    joint_id: int
    dorsal_port: PortRef
    arc_port: PortRef
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
class ResolvedJoint:
    """Validated mechanical joint ready for SoA packing."""

    joint_id: int
    kind: JointKind
    port_a: PortRef
    port_b: PortRef
    stiffness_pn_per_um: float
    rest_um: float
    initial_state: JointState
    semantic_edges: tuple[str, ...]
    fa_cluster_id: int = -1


def resolve_fa_series_joint(
    spec: CompositeFAJointSpec,
    registry: ActorRegistry,
    architecture: CellArchitecture,
) -> ResolvedJoint:
    """Validate endpoint roles and create exactly one SF--ECM mechanical joint."""
    resolution = resolve_fa_series_group(architecture)
    registry.validate_port(spec.sf_port)
    registry.validate_port(spec.ecm_port)
    if spec.sf_port.component != "sf_arc":
        raise ValueError("FA actin-side port must belong to sf_arc")
    if spec.sf_port.role.contract_name not in _FA_ACTIN_ROLES:
        raise ValueError("FA actin-side role must be a ventral end or dorsal basal end")
    if spec.ecm_port.component != "ecm" or spec.ecm_port.role is not EndpointRole.COLLAGEN_LIGAND:
        raise ValueError("FA ligand-side port must be a collagen ligand on ecm")
    _ = spec.sf_port.segment_u
    _ = spec.ecm_port.segment_u
    if resolution.mechanical_joint_count != 1:
        raise RuntimeError("FA semantic edges would double-count mechanical stiffness")
    return ResolvedJoint(
        joint_id=spec.joint_id,
        kind=JointKind.FA_COMPOSITE_SERIES,
        port_a=spec.sf_port,
        port_b=spec.ecm_port,
        stiffness_pn_per_um=spec.stiffness_pn_per_um,
        rest_um=spec.rest_um,
        initial_state=spec.initial_state,
        semantic_edges=resolution.semantic_edges,
        fa_cluster_id=spec.fa_cluster_id,
    )


def resolve_dorsal_arc_joint(
    spec: DorsalArcJointSpec,
    registry: ActorRegistry,
    architecture: CellArchitecture,
) -> ResolvedJoint:
    """Validate the typed same-component dorsal-free--arc internal connector."""
    registry.validate_port(spec.dorsal_port)
    registry.validate_port(spec.arc_port)
    if spec.dorsal_port.component != "sf_arc" or spec.arc_port.component != "sf_arc":
        raise ValueError("dorsal--arc connector must stay inside sf_arc")
    if spec.dorsal_port.role is not EndpointRole.DORSAL_FREE_END:
        raise ValueError("dorsal--arc endpoint A must be the dorsal free end")
    if spec.arc_port.role is not EndpointRole.TRANSVERSE_ARC_MATERIAL:
        raise ValueError("dorsal--arc endpoint B must be transverse-arc material")
    _ = spec.dorsal_port.segment_u
    _ = spec.arc_port.segment_u
    connector = next(
        (item for item in architecture.connectors if item.name == "dorsal_arc_crosslink"),
        None,
    )
    if connector is None or connector.scope is not ConnectorScope.INTERNAL:
        raise ValueError("architecture has no internal dorsal--arc connector contract")
    return ResolvedJoint(
        joint_id=spec.joint_id,
        kind=JointKind.DORSAL_ARC_INTERNAL,
        port_a=spec.dorsal_port,
        port_b=spec.arc_port,
        stiffness_pn_per_um=spec.stiffness_pn_per_um,
        rest_um=spec.rest_um,
        initial_state=spec.initial_state,
        semantic_edges=(connector.name,),
    )


@dataclass(frozen=True, slots=True)
class SegmentPairReference:
    """Host virtual-work oracle for one interpolated central-force segment joint."""

    force_a: npt.NDArray[np.float64]
    force_b: npt.NDArray[np.float64]
    joint_force_on_a: npt.NDArray[np.float64]
    load_pn: float
    energy_pn_um: float
    nodal_virtual_work_pn_um: float | None


def segment_pair_reference(
    pos_a: npt.ArrayLike,
    segment_a: tuple[int, int],
    u_a: float,
    pos_b: npt.ArrayLike,
    segment_b: tuple[int, int],
    u_b: float,
    stiffness_pn_per_um: float,
    rest_um: float,
    *,
    displacement_a: npt.ArrayLike | None = None,
    displacement_b: npt.ArrayLike | None = None,
) -> SegmentPairReference:
    """Scatter one Hookean pair and optionally evaluate nodal virtual work.

    ``joint_force_on_a`` points from endpoint A toward B in extension.  Nodal virtual
    work is ``sum(F_i dot delta_x_i) = -dU`` for the supplied displacement direction.
    """
    pa = np.asarray(pos_a, dtype=np.float64)
    pb = np.asarray(pos_b, dtype=np.float64)
    if pa.ndim != 2 or pb.ndim != 2 or pa.shape[1:] != (3,) or pb.shape[1:] != (3,):
        raise ValueError("actor positions must be (N, 3)")
    if not 0.0 <= u_a <= 1.0 or not 0.0 <= u_b <= 1.0:
        raise ValueError("segment coordinates must be in [0, 1]")
    ia0, ia1 = segment_a
    ib0, ib1 = segment_b
    xa = (1.0 - u_a) * pa[ia0] + u_a * pa[ia1]
    xb = (1.0 - u_b) * pb[ib0] + u_b * pb[ib1]
    delta = xb - xa
    length = float(np.linalg.norm(delta))
    force = np.zeros(3, dtype=np.float64)
    extension = length - rest_um
    if length > 1.0e-15:
        force = stiffness_pn_per_um * extension * delta / length
    fa = np.zeros_like(pa)
    fb = np.zeros_like(pb)
    fa[ia0] += (1.0 - u_a) * force
    fa[ia1] += u_a * force
    fb[ib0] -= (1.0 - u_b) * force
    fb[ib1] -= u_b * force
    work = None
    if displacement_a is not None or displacement_b is not None:
        if displacement_a is None or displacement_b is None:
            raise ValueError("both actor displacement arrays are required for virtual work")
        da = np.asarray(displacement_a, dtype=np.float64)
        db = np.asarray(displacement_b, dtype=np.float64)
        if da.shape != pa.shape or db.shape != pb.shape:
            raise ValueError("displacement arrays must match actor positions")
        work = float(np.sum(fa * da) + np.sum(fb * db))
    return SegmentPairReference(
        force_a=fa,
        force_b=fb,
        joint_force_on_a=force,
        load_pn=float(abs(stiffness_pn_per_um * extension)),
        energy_pn_um=float(0.5 * stiffness_pn_per_um * extension**2),
        nodal_virtual_work_pn_um=work,
    )


@wp.kernel
def _zero_vec3_kernel(values: wp.array(dtype=wp.vec3d)) -> None:
    values[wp.tid()] = wp.vec3d(0.0, 0.0, 0.0)


@wp.kernel
def _segment_joint_force_kernel(
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
    """One thread per joint; gather and scatter use the same barycentric weights."""
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
    pair_force = wp.vec3d(0.0, 0.0, 0.0)
    extension = length - rest[t]
    if length > wp.float64(1.0e-15):
        pair_force = stiffness[t] * extension * delta / length
    wp.atomic_add(force_a, ia0, wa0 * pair_force)
    wp.atomic_add(force_a, ia1, wa1 * pair_force)
    wp.atomic_add(force_b, ib0, -wb0 * pair_force)
    wp.atomic_add(force_b, ib1, -wb1 * pair_force)
    joint_force_on_a[t] = pair_force
    load[t] = wp.abs(stiffness[t] * extension)
    energy[t] = wp.float64(0.5) * stiffness[t] * extension * extension


@wp.kernel
def _refresh_rest_on_bind_kernel(
    pos_a: wp.array(dtype=wp.vec3d),
    segments_a: wp.array(dtype=wp.int32, ndim=2),
    pos_b: wp.array(dtype=wp.vec3d),
    segments_b: wp.array(dtype=wp.int32, ndim=2),
    active: wp.array(dtype=wp.int32),
    accepted: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    candidate_state: wp.array(dtype=wp.int32),
    element_a: wp.array(dtype=wp.int32),
    element_b: wp.array(dtype=wp.int32),
    u_a: wp.array(dtype=wp.float64),
    u_b: wp.array(dtype=wp.float64),
    rest: wp.array(dtype=wp.float64),
) -> None:
    """On a committed FREE→ENGAGED, set the rest length to the separation AT THAT MOMENT.

    The ``r0_bind`` convention (`bc3fef26`) — every joint starts force-free.  It held trivially while
    nothing moved, because the build-time separation WAS the separation at bind.  The moment a
    component integrates its own positions the two diverge, and a rebind at the stale rest injects
    ``k·Δ`` from nowhere: at α-actinin's 4.6e5 pN/µm a 10 nm drift is 4,600 pN, two orders above the
    bond's own rupture force.  Nothing would raise; the cell would simply be pulled by an artifact.

    Geometry is computed with the SAME barycentric weights as
    :func:`_segment_joint_force_kernel` — a different convention here would be its own defect.

    Ordering: ``accepted`` gates the write, so a rejected step leaves ``rest`` untouched and rollback
    stays correct by construction.  This must launch BEFORE ``_commit_joint_state_kernel``, which
    overwrites ``state`` and so destroys the evidence of the transition.
    """
    t = wp.tid()
    if active[t] == 0 or accepted[0] == 0:
        return
    if state[t] >= _ACTIVE_STATE_MIN or candidate_state[t] < _ACTIVE_STATE_MIN:
        return  # not a FREE→ENGAGED this step; an already-engaged joint keeps its own rest
    ea = element_a[t]
    eb = element_b[t]
    wa1 = u_a[t]
    wa0 = wp.float64(1.0) - wa1
    wb1 = u_b[t]
    wb0 = wp.float64(1.0) - wb1
    xa = wa0 * pos_a[segments_a[ea, 0]] + wa1 * pos_a[segments_a[ea, 1]]
    xb = wb0 * pos_b[segments_b[eb, 0]] + wb1 * pos_b[segments_b[eb, 1]]
    rest[t] = wp.length(xb - xa)


@wp.kernel
def _commit_joint_state_kernel(
    active: wp.array(dtype=wp.int32),
    accepted: wp.array(dtype=wp.int32),
    candidate_state: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    age: wp.array(dtype=wp.float64),
    rng_epoch: wp.array(dtype=wp.int32),
    dt_phys: wp.float64,
) -> None:
    """Commit candidate topology only under the scheduler-owned acceptance scalar."""
    t = wp.tid()
    if active[t] != 0 and accepted[0] != 0:
        previous = state[t]
        state[t] = candidate_state[t]
        if candidate_state[t] == previous:
            age[t] = age[t] + dt_phys
        else:
            age[t] = wp.float64(0.0)
        rng_epoch[t] = rng_epoch[t] + wp.int32(1)
    elif accepted[0] == 0:
        candidate_state[t] = state[t]


@wp.kernel
def _rollback_candidate_state_kernel(
    accepted: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    candidate_state: wp.array(dtype=wp.int32),
) -> None:
    """Reset scratch candidates after rejection; committed state is never mutated."""
    t = wp.tid()
    if accepted[0] == 0:
        candidate_state[t] = state[t]


class ActorSegmentRuntime:
    """Minimal CUDA-resident segment actor used by the VS-0 load-path slice."""

    def __init__(
        self,
        record: ActorRecord,
        positions: npt.ArrayLike,
        segments: npt.ArrayLike,
        *,
        device: str | None = None,
    ) -> None:
        dev = wp.get_device(device)
        if not dev.is_cuda:
            raise RuntimeError("ActorSegmentRuntime requires a CUDA device")
        pos = np.ascontiguousarray(positions, dtype=np.float64)
        seg = np.ascontiguousarray(segments, dtype=np.int32)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError("positions must be (N, 3)")
        if seg.ndim != 2 or seg.shape[1] != 2:
            raise ValueError("segments must be (M, 2)")
        if seg.size and (seg.min() < 0 or seg.max() >= pos.shape[0]):
            raise ValueError("segment topology references an invalid node")
        if record.n_elements != seg.shape[0]:
            raise ValueError("actor record element count does not match segment topology")
        self.record = record
        self.device = str(dev)
        self.n_nodes = int(pos.shape[0])
        self.segments_host = seg
        with wp.ScopedDevice(self.device):
            self.pos_d = wp.array(pos, dtype=wp.vec3d, device=self.device)
            self.force_d = wp.zeros(self.n_nodes, dtype=wp.vec3d, device=self.device)
            self.segments_d = wp.array(seg, dtype=wp.int32, ndim=2, device=self.device)

    def zero_force(self) -> None:
        """Clear this actor's force accumulator on device."""
        wp.launch(_zero_vec3_kernel, dim=self.n_nodes, inputs=[self.force_d], device=self.device)


class LoadPathJointRuntime:
    """Fixed-capacity CUDA SoA for one actor-pair's resolved segment joints.

    ``force_kernel`` selects the pair law the SoA is launched with.  It exists because the SoA, the
    generation-checked ports, the barycentric gather/scatter adjoint and the whole accepted-step
    transaction API are identical for every two-array joint in the graph, while the pair law is not:
    a crosslink is a bilateral Hookean spring, a contact is a UNILATERAL one that may only push.
    Every admitted kernel must take the exact argument list of
    :func:`_segment_joint_force_kernel` and must scatter ``+f`` through A's weights and ``-f`` through
    B's, because that adjoint is what the balance gate tests; what it may change is when and in which
    direction ``f`` is non-zero.
    """

    def __init__(
        self,
        actor_a: ActorSegmentRuntime | BorrowedSegmentActorView,
        actor_b: ActorSegmentRuntime | BorrowedSegmentActorView,
        joints: tuple[ResolvedJoint, ...],
        *,
        force_kernel: object | None = None,
    ) -> None:
        if actor_a.device != actor_b.device:
            raise ValueError("joint actors must reside on the same CUDA device")
        if not joints:
            raise ValueError("joint runtime needs at least one resolved joint")
        ids = tuple(joint.joint_id for joint in joints)
        if len(ids) != len(set(ids)):
            raise ValueError("joint IDs must be unique")
        for joint in joints:
            if joint.port_a.actor_id != actor_a.record.actor_id:
                raise ValueError("joint endpoint A does not belong to actor A")
            if joint.port_b.actor_id != actor_b.record.actor_id:
                raise ValueError("joint endpoint B does not belong to actor B")
        self.actor_a = actor_a
        self.actor_b = actor_b
        self.device = actor_a.device
        self.capacity = len(joints)
        self.force_kernel = force_kernel if force_kernel is not None else _segment_joint_force_kernel

        def i32(values: list[int]) -> wp.array:
            return wp.array(np.asarray(values, dtype=np.int32), dtype=wp.int32, device=self.device)

        def f64(values: list[float]) -> wp.array:
            return wp.array(np.asarray(values, dtype=np.float64), dtype=wp.float64, device=self.device)

        with wp.ScopedDevice(self.device):
            self.active_d = i32([1] * self.capacity)
            self.joint_id_d = i32([joint.joint_id for joint in joints])
            self.joint_kind_d = i32([int(joint.kind) for joint in joints])
            self.state_d = i32([int(joint.initial_state) for joint in joints])
            self.candidate_state_d = i32([int(joint.initial_state) for joint in joints])
            self.fa_cluster_id_d = i32([joint.fa_cluster_id for joint in joints])
            self.actor_a_id_d = i32([joint.port_a.actor_id for joint in joints])
            self.actor_a_generation_d = i32(
                [joint.port_a.actor_generation for joint in joints]
            )
            self.entity_a_id_d = i32([joint.port_a.entity_id for joint in joints])
            self.entity_a_generation_d = i32(
                [joint.port_a.entity_generation for joint in joints]
            )
            self.element_a_kind_d = i32([int(joint.port_a.element_kind) for joint in joints])
            self.element_a_id_d = i32([joint.port_a.element_id for joint in joints])
            self.role_a_d = i32([int(joint.port_a.role) for joint in joints])
            self.u_a_d = f64([joint.port_a.segment_u for joint in joints])
            self.actor_b_id_d = i32([joint.port_b.actor_id for joint in joints])
            self.actor_b_generation_d = i32(
                [joint.port_b.actor_generation for joint in joints]
            )
            self.entity_b_id_d = i32([joint.port_b.entity_id for joint in joints])
            self.entity_b_generation_d = i32(
                [joint.port_b.entity_generation for joint in joints]
            )
            self.element_b_kind_d = i32([int(joint.port_b.element_kind) for joint in joints])
            self.element_b_id_d = i32([joint.port_b.element_id for joint in joints])
            self.role_b_d = i32([int(joint.port_b.role) for joint in joints])
            self.u_b_d = f64([joint.port_b.segment_u for joint in joints])
            self.stiffness_d = f64([joint.stiffness_pn_per_um for joint in joints])
            self.rest_d = f64([joint.rest_um for joint in joints])
            self.age_d = f64([0.0] * self.capacity)
            self.rng_epoch_d = i32([0] * self.capacity)
            self.load_d = f64([0.0] * self.capacity)
            self.energy_d = f64([0.0] * self.capacity)
            self.joint_force_on_a_d = wp.zeros(
                self.capacity, dtype=wp.vec3d, device=self.device
            )

    def accumulate(self) -> None:
        """Add one mechanical contribution per joint to the two actor force arrays."""
        wp.launch(
            self.force_kernel,
            dim=self.capacity,
            inputs=[
                self.actor_a.pos_d,
                self.actor_a.force_d,
                self.actor_a.segments_d,
                self.actor_b.pos_d,
                self.actor_b.force_d,
                self.actor_b.segments_d,
                self.active_d,
                self.state_d,
                self.element_a_id_d,
                self.element_b_id_d,
                self.u_a_d,
                self.u_b_d,
                self.stiffness_d,
                self.rest_d,
            ],
            outputs=[self.load_d, self.energy_d, self.joint_force_on_a_d],
            device=self.device,
        )

    def snapshot_candidate(self) -> None:
        """Reset candidate topology to the current committed state by a D2D copy."""
        wp.copy(self.candidate_state_d, self.state_d)

    def rollback(self, accepted: wp.array) -> None:
        """Discard scratch transitions under a rejected outer-step predicate."""
        wp.launch(
            _rollback_candidate_state_kernel,
            dim=self.capacity,
            inputs=[accepted, self.state_d, self.candidate_state_d],
            device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit graph state, age, and RNG epoch only when ``accepted[0] != 0``.

        ``rng_seed`` is part of the common transaction hook.  VS-0 receives candidate
        transitions from a caller; the future GPU hazard proposal consumes the seed.
        """
        if not np.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if rng_seed < 0:
            raise ValueError("rng_seed must be nonnegative")
        # BEFORE the state commit, which overwrites `state` and destroys the FREE→ENGAGED evidence.
        wp.launch(
            _refresh_rest_on_bind_kernel,
            dim=self.capacity,
            inputs=[
                self.actor_a.pos_d, self.actor_a.segments_d,
                self.actor_b.pos_d, self.actor_b.segments_d,
                self.active_d, accepted, self.state_d, self.candidate_state_d,
                self.element_a_id_d, self.element_b_id_d, self.u_a_d, self.u_b_d,
            ],
            outputs=[self.rest_d],
            device=self.device,
        )
        wp.launch(
            _commit_joint_state_kernel,
            dim=self.capacity,
            inputs=[
                self.active_d,
                accepted,
                self.candidate_state_d,
                self.state_d,
                self.age_d,
                self.rng_epoch_d,
                wp.float64(dt_phys),
            ],
            device=self.device,
        )
