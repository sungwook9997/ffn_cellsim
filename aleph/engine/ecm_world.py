"""Component-first runtime seam for the collagen ECM world.

This module fixes ECM ownership and coupling before any new constitutive or kinetic law is
introduced.  ``ecm`` owns live collagen geometry plus damage, sleep, refinement, and topology
metadata on CUDA.  Collagen crosslinks and the alpha2beta1--collagen clutch remain graph-owned
connectors, while the far-field anchor is an explicit finite-domain boundary participant.

The collagen endpoint is a segment material point, not a permanently captured node.  Its
barycentric gather and scatter are therefore compatible with the adjoint segment-pair mechanics
in :mod:`aleph.engine.load_path`, and its topology epoch lets the accepted-step transaction
invalidate stale views after remeshing.

This is an ownership/orchestration slice.  Injected delegates launch the existing Warp kernels;
the facade never reads authoritative device values, performs a host neighbour query, advances a
private physical clock, or commits topology outside the global acceptance predicate.

Sanity Gate:
    * ownership: every authoritative ECM array is CUDA-resident, non-aliased, and disjoint from
      connector state;
    * physiological baseline: an unanchored finite ECM representative volume is rejected;
    * load path: the two FA semantic edges resolve to one alpha2beta1--collagen series runtime;
    * virtual work: a clutch endpoint is a segment barycentric material point and requires an
      equal-and-opposite, adjoint force scatter plus a work ledger;
    * topology: sleeping/refinement/damage state and collagen crosslinks participate in the same
      accepted-step transaction; and
    * accounting: component, crosslink, clutch, and boundary-reaction ledgers are all visited.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum
from typing import TYPE_CHECKING, Callable, Protocol, runtime_checkable

import warp as wp

if TYPE_CHECKING:
    from aleph.components.ecm.device_schema import ECMTopologyState
    from aleph.engine.ecm_mechanics import ECMConstitutiveForce
    from aleph.laws.ecm_library import ECMSpec

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentRole,
    ConnectorFamily,
    ConnectorScope,
)
from aleph.engine.ledger import GlobalCellLedger
from aleph.engine.load_path import (
    ActorRecord,
    ElementKind,
    EndpointRole,
    PortRef,
    resolve_fa_series_group,
)
from aleph.engine.runtime import (
    LedgerContributor,
    MechanicsContributor,
    TransactionParticipant,
)

ECM_COMPONENT = "ecm"
ECM_CROSSLINK_CONNECTOR = "ecm_crosslink"
FA_ECM_CLUTCH_CONNECTOR = "integrin_collagen_clutch"
FA_SERIES_GROUP = "alpha2beta1_collagen_series"
ECM_FAR_FIELD_ANCHOR = "ecm_far_field_anchor"
MEMBRANE_COMPONENT = "membrane"
MEMBRANE_ECM_CONTACT = "membrane_ecm_contact"


def _device_is_cuda(array: object) -> bool:
    """Return whether an array-like object declares a CUDA device."""
    device = getattr(array, "device", None)
    return bool(getattr(device, "is_cuda", False))


def _storage_key(array: object) -> tuple[object, ...]:
    """Return a conservative storage identity without copying device data."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


def _validate_device_array(
    array: object,
    *,
    label: str,
    dtype: object,
    ndim: int,
    nonempty: bool = True,
) -> None:
    """Validate CUDA array metadata without reading authoritative state."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != dtype:
        raise TypeError(f"{label} must have dtype {dtype}")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != ndim:
        raise ValueError(f"{label} must be a {ndim}-dimensional array")
    if nonempty and any(int(size) <= 0 for size in shape):
        raise ValueError(f"{label} must have no empty dimension")


def _require_same_device(arrays: tuple[object, ...], *, label: str) -> None:
    """Require all arrays to occupy one CUDA device."""
    devices = {str(getattr(array, "device", None)) for array in arrays}
    if len(devices) != 1:
        raise ValueError(f"{label} arrays must share one CUDA device")


def _require_distinct_storage(arrays: tuple[object, ...], *, label: str) -> None:
    """Reject authoritative arrays that alias device storage."""
    keys = tuple(_storage_key(array) for array in arrays)
    if len(set(keys)) != len(keys):
        raise ValueError(f"{label} arrays must own distinct storage")


class BoundaryAnchorMode(StrEnum):
    """Physiological closures allowed for a finite ECM representative volume."""

    FAR_FIELD_DIRICHLET = "far_field_dirichlet"
    CONTINUUM_IMPEDANCE = "continuum_impedance"


@dataclass(frozen=True, slots=True)
class ECMWorldSettings:
    """Mechanism switches for the first ECM world.

    No numerical material value lives here: stiffness, damage, remodelling, and kinetic parameters
    stay in literature-backed per-unit configuration.  These switches only reject structurally
    non-physiological or non-transactional production builds.
    """

    anchor_mode: BoundaryAnchorMode
    physiological_anchor_required: bool = True
    crosslink_kinetics: bool = True
    remodelling_enabled: bool = True
    damage_enabled: bool = True
    sleeping_enabled: bool = True
    refinement_enabled: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.anchor_mode, BoundaryAnchorMode):
            raise TypeError("anchor_mode must be a BoundaryAnchorMode")
        if not self.physiological_anchor_required:
            raise ValueError("production ECM requires a physiological far-field anchor")
        if not self.crosslink_kinetics:
            raise ValueError("collagen crosslink kinetics cannot be disabled in production ECM")
        if not self.remodelling_enabled or not self.damage_enabled:
            raise ValueError("production ECM requires remodelling and damage state")
        if not self.sleeping_enabled or not self.refinement_enabled:
            raise ValueError("ECM world requires conservative sleeping and refinement state")


@runtime_checkable
class ECMInternalCrosslink(Protocol):
    """Graph-owned collagen--collagen connector population."""

    name: str
    component_a: str
    component_b: str
    adjoint_transfer_required: bool

    def accumulate_internal(self, pos: wp.array, force: wp.array) -> None:
        """Add crosslink forces to the ECM-owned force array."""

    def snapshot_candidate(self) -> None:
        """Snapshot committed crosslink topology before candidate proposals."""

    def rollback(self, accepted: wp.array) -> None:
        """Discard rejected crosslink proposals under the scheduler predicate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit crosslink association/dissociation only after global acceptance."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute force/work/energy and topology counts."""


@runtime_checkable
class CompositeECMClutch(Protocol):
    """One graph-owned SF--FA--collagen composite series runtime."""

    mechanical_group: str
    ecm_endpoint: CollagenMaterialPointEndpointView
    equal_opposite_required: bool
    adjoint_transfer_required: bool
    work_ledger_required: bool

    def accumulate(self) -> None:
        """Add the single series-joint mechanics to both actor force arrays."""

    def snapshot_candidate(self) -> None:
        """Snapshot clutch state before candidate kinetics."""

    def rollback(self, accepted: wp.array) -> None:
        """Discard rejected clutch proposals."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit clutch kinetics only after global acceptance."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute equal/opposite force, moment, energy, and work closure terms."""


@runtime_checkable
class CompositeJointDelegate(Protocol):
    """Load-path joint mechanics/transaction seam wrapped by the ECM facade."""

    actor_a: object
    actor_b: object

    def accumulate(self) -> None:
        """Scatter the single composite series load to SF and live collagen arrays."""

    def snapshot_candidate(self) -> None:
        """Snapshot candidate clutch state."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore a rejected candidate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit clutch kinetics only after acceptance."""


@dataclass(frozen=True, slots=True)
class CompositeLoadPathECMClutchAdapter:
    """Make one load-path joint the sole SF--FA--ECM composite runtime."""

    joint: CompositeJointDelegate
    ecm_endpoint: CollagenMaterialPointEndpointView
    ledger: LedgerContributor
    mechanical_group: str = FA_SERIES_GROUP
    equal_opposite_required: bool = True
    adjoint_transfer_required: bool = True
    work_ledger_required: bool = True

    def __post_init__(self) -> None:
        actors = (self.joint.actor_a, self.joint.actor_b)
        components = {actor.record.component for actor in actors}
        if components != {"sf_arc", ECM_COMPONENT}:
            raise ValueError("composite FA joint must join sf_arc and ecm segment actors")
        ecm_actor = next(actor for actor in actors if actor.record.component == ECM_COMPONENT)
        if ecm_actor.record != self.ecm_endpoint.owner.actor:
            raise ValueError("composite FA joint ECM actor identity does not match ECM owner")
        if (
            _storage_key(ecm_actor.pos_d) != _storage_key(self.ecm_endpoint.position_d)
            or _storage_key(ecm_actor.force_d) != _storage_key(self.ecm_endpoint.force_d)
            or _storage_key(ecm_actor.segments_d) != _storage_key(self.ecm_endpoint.segments_d)
        ):
            raise ValueError("composite FA joint must borrow the live ECM owner arrays")

    def accumulate(self) -> None:
        """Dispatch the mechanical group exactly once."""
        self.joint.accumulate()

    def snapshot_candidate(self) -> None:
        """Forward clutch snapshot."""
        self.joint.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Forward rejected-candidate restoration."""
        self.joint.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Forward accepted clutch commit."""
        self.joint.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward the joint force/work/moment ledger."""
        self.ledger.accumulate_ledger(ledger)


@runtime_checkable
class BoundaryAnchorDelegate(Protocol):
    """Injected far-field boundary mechanics and accepted-state cache."""

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        """Apply the finite-domain far-field closure to ECM force/residual arrays."""

    def snapshot_candidate(self) -> None:
        """Snapshot prescribed-target or impedance history before a candidate."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore boundary history under the rejected predicate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted boundary history without advancing a private clock."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute far-field reaction, moment, and boundary-work terms."""


@dataclass(frozen=True, slots=True)
class MembraneContactSurfaceView:
    """Non-owning live membrane surface arrays used by collagen contact."""

    component: str
    position_d: wp.array
    force_d: wp.array

    def __post_init__(self) -> None:
        if self.component != MEMBRANE_COMPONENT:
            raise ValueError("ECM contact surface must be the registered membrane")
        for label, array in (("position_d", self.position_d), ("force_d", self.force_d)):
            _validate_device_array(
                array,
                label=f"membrane.{label}",
                dtype=wp.vec3d,
                ndim=1,
            )
        if self.position_d.shape != self.force_d.shape:
            raise ValueError("membrane contact position and force shapes must match")
        _require_same_device((self.position_d, self.force_d), label="membrane contact")
        _require_distinct_storage((self.position_d, self.force_d), label="membrane contact")


@runtime_checkable
class ECMMembraneContact(Protocol):
    """Graph-owned reversible segment--surface contact and mapping cache."""

    name: str
    component_a: str
    component_b: str
    adjoint_transfer_required: bool

    def accumulate_contact(
        self,
        membrane: MembraneContactSurfaceView,
        ecm: CollagenMaterialPointEndpointView,
    ) -> None:
        """Apply equal-and-opposite live membrane--collagen segment contact."""

    def snapshot_candidate(self) -> None:
        """Snapshot reversible contact candidate and broad-phase cache state."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore contact cache state under rejection."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted remap/cache state without contact kinetics."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute contact force, moment, work, and stale-generation counts."""


@dataclass(frozen=True, slots=True)
class ECMBoundaryAnchorFacade:
    """Explicit physiological closure for the truncated ECM domain."""

    delegate: BoundaryAnchorDelegate
    mode: BoundaryAnchorMode
    name: str = ECM_FAR_FIELD_ANCHOR
    component: str = ECM_COMPONENT
    component_a: str = ECM_COMPONENT
    component_b: str = "world_boundary"
    reaction_ledger_required: bool = True
    work_ledger_required: bool = True

    def __post_init__(self) -> None:
        if self.name != ECM_FAR_FIELD_ANCHOR:
            raise ValueError(f"expected ECM boundary {ECM_FAR_FIELD_ANCHOR!r}")
        if self.component != ECM_COMPONENT:
            raise ValueError("ECM far-field anchor must address the ecm component")
        if {self.component_a, self.component_b} != {ECM_COMPONENT, "world_boundary"}:
            raise ValueError("ECM far-field anchor endpoints must be ecm and world_boundary")
        if not isinstance(self.mode, BoundaryAnchorMode):
            raise TypeError("ECM boundary mode must be a BoundaryAnchorMode")
        if not self.reaction_ledger_required or not self.work_ledger_required:
            raise ValueError("ECM boundary must ledger reaction and boundary work")

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Apply boundary mechanics without hiding reaction forces."""
        self.delegate.accumulate_boundary(pos, force)

    def snapshot_candidate(self) -> None:
        """Forward the transaction snapshot hook."""
        self.delegate.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Forward the scheduler-owned predicate without host inspection."""
        self.delegate.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted boundary history only."""
        self.delegate.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward boundary reaction/work accounting."""
        self.delegate.accumulate_ledger(ledger)


@wp.kernel
def _far_field_reaction_ledger_kernel(
    force: wp.array(dtype=wp.vec3d),
    pinned_index: wp.array(dtype=wp.int32),
    reaction: wp.array(dtype=wp.vec3d),
) -> None:
    """Record then remove the Dirichlet constraint reaction on each pinned ECM node.

    One thread per pinned node.  ``reaction[t]`` stores the constraint reaction ``R = -f``
    (the internal force the far-field pin resists) *before* the node's force is zeroed, so the
    cell traction that reaches the truncated boundary is preserved in the ledger rather than
    silently discarded.  This is finite-domain boundary bookkeeping, not a collagen material law.
    """
    t = wp.tid()
    k = pinned_index[t]
    f = force[k]
    reaction[t] = -f
    force[k] = wp.vec3d(0.0, 0.0, 0.0)


@wp.kernel
def _far_field_reaction_rollback_kernel(
    accepted: wp.array(dtype=wp.int32),
    committed_reaction: wp.array(dtype=wp.vec3d),
    reaction: wp.array(dtype=wp.vec3d),
) -> None:
    """Discard the rejected candidate reaction; committed accepted reaction is never mutated."""
    t = wp.tid()
    if accepted[0] == 0:
        reaction[t] = committed_reaction[t]


@wp.kernel
def _far_field_reaction_commit_kernel(
    accepted: wp.array(dtype=wp.int32),
    reaction: wp.array(dtype=wp.vec3d),
    target: wp.array(dtype=wp.vec3d),
    committed_reaction: wp.array(dtype=wp.vec3d),
    committed_target: wp.array(dtype=wp.vec3d),
    boundary_work: wp.array(dtype=wp.float64),
) -> None:
    """Finalize accepted reaction and accumulate reaction-conjugate boundary work.

    For a static far-field target (the first MCF7 baseline) ``target == committed_target`` so the
    boundary-work increment is exactly zero; a prescribed moving anchor contributes
    ``sum(reaction . delta_target)``.  No private clock is advanced.
    """
    t = wp.tid()
    if accepted[0] != 0:
        dwork = wp.dot(reaction[t], target[t] - committed_target[t])
        wp.atomic_add(boundary_work, 0, dwork)
        committed_reaction[t] = reaction[t]
        committed_target[t] = target[t]


@dataclass(frozen=True, slots=True)
class FarFieldDirichletBoundaryRuntime:
    """Warp-kernel-backed reaction-preserving far-field Dirichlet closure.

    This concrete :class:`BoundaryAnchorDelegate` binds real ported Warp kernels behind the ECM
    boundary seam (evidence state KERNEL_BOUND): it owns the pinned-node index set, the physiological
    bulk-ECM target positions, and the reaction/boundary-work ledger arrays.  It never reads an
    authoritative device value on the host, advances no private physical clock, and mutates
    committed state only under the scheduler-owned acceptance scalar.

    A non-empty pinned set is required: an unanchored finite representative volume is rejected here as
    well as in :class:`ECMWorldSettings`.
    """

    pinned_index_d: wp.array
    target_position_d: wp.array
    reaction_d: wp.array
    committed_reaction_d: wp.array
    committed_target_d: wp.array
    boundary_work_d: wp.array

    def __post_init__(self) -> None:
        _validate_device_array(
            self.pinned_index_d,
            label="boundary.pinned_index_d",
            dtype=wp.int32,
            ndim=1,
        )
        n_pinned = int(self.pinned_index_d.shape[0])
        vec_arrays = (
            ("target_position_d", self.target_position_d),
            ("reaction_d", self.reaction_d),
            ("committed_reaction_d", self.committed_reaction_d),
            ("committed_target_d", self.committed_target_d),
        )
        for label, array in vec_arrays:
            _validate_device_array(array, label=f"boundary.{label}", dtype=wp.vec3d, ndim=1)
            if array.shape != (n_pinned,):
                raise ValueError(f"boundary.{label} must have one entry per pinned far-field node")
        _validate_device_array(
            self.boundary_work_d,
            label="boundary.boundary_work_d",
            dtype=wp.float64,
            ndim=1,
        )
        if self.boundary_work_d.shape != (1,):
            raise ValueError("boundary.boundary_work_d must hold exactly one accepted work scalar")
        arrays = (
            self.pinned_index_d,
            self.target_position_d,
            self.reaction_d,
            self.committed_reaction_d,
            self.committed_target_d,
            self.boundary_work_d,
        )
        _require_same_device(arrays, label="ECM far-field boundary")
        _require_distinct_storage(arrays, label="ECM far-field boundary")

    @property
    def device(self) -> str:
        """CUDA device shared by every boundary ledger array."""
        return str(self.pinned_index_d.device)

    @property
    def n_pinned(self) -> int:
        """Number of far-field pinned collagen nodes from host-visible metadata."""
        return int(self.pinned_index_d.shape[0])

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        """Record and remove the pinned-node reaction on the ECM force accumulator."""
        wp.launch(
            _far_field_reaction_ledger_kernel,
            dim=self.n_pinned,
            inputs=[force, self.pinned_index_d, self.reaction_d],
            device=self.device,
        )

    def snapshot_candidate(self) -> None:
        """Reset the candidate reaction to the last accepted reaction by a device copy."""
        wp.copy(self.reaction_d, self.committed_reaction_d)

    def rollback(self, accepted: wp.array) -> None:
        """Restore committed reaction under a rejected outer-step predicate."""
        wp.launch(
            _far_field_reaction_rollback_kernel,
            dim=self.n_pinned,
            inputs=[accepted, self.committed_reaction_d, self.reaction_d],
            device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted reaction and boundary work; no private clock is advanced."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if isinstance(rng_seed, bool) or not isinstance(rng_seed, int) or rng_seed < 0:
            raise ValueError("rng_seed must be a nonnegative integer")
        wp.launch(
            _far_field_reaction_commit_kernel,
            dim=self.n_pinned,
            inputs=[
                accepted,
                self.reaction_d,
                self.target_position_d,
                self.committed_reaction_d,
                self.committed_target_d,
                self.boundary_work_d,
            ],
            device=self.device,
        )

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward the committed far-field reaction resultant and boundary work to the global cell ledger.

        The typed :class:`~aleph.engine.ledger.GlobalCellLedger` now owns the ``add_far_field_reaction``
        slot: it reduces the committed per-pinned-node reaction into the boundary-reaction channel of the
        global force-balance gate and accumulates the boundary work — on the device, never a host readback.
        A ledger sink lacking the hook (a partial bring-up sink) is still tolerated as a schema gap, and no
        authoritative device value is read either way.
        """
        if isinstance(ledger, GlobalCellLedger):
            ledger.add_far_field_reaction(self.committed_reaction_d, self.boundary_work_d)
            return
        recorder = getattr(ledger, "add_far_field_reaction", None)
        if callable(recorder):
            recorder(self.committed_reaction_d, self.boundary_work_d)


@dataclass(frozen=True, slots=True)
class ECMStateOwner:
    """The ECM component and the CUDA arrays it exclusively owns.

    Connector bond states are deliberately absent.  ``endpoint_refcount_d`` is an ECM-owned
    conservative lock/cache used to prevent sleeping or remeshing underneath a live graph endpoint;
    the clutch remains authoritative for bond identity and chemistry.
    """

    name: str
    settings: ECMWorldSettings
    actor: ActorRecord
    position_d: wp.array
    force_d: wp.array
    segments_d: wp.array
    fiber_offsets_d: wp.array
    fiber_sleep_state_d: wp.array
    segment_refinement_level_d: wp.array
    segment_damage_d: wp.array
    endpoint_refcount_d: wp.array
    topology_epoch_d: wp.array
    mechanics: MechanicsContributor
    transaction: TransactionParticipant
    ledger: LedgerContributor

    def __post_init__(self) -> None:
        if self.name != ECM_COMPONENT or self.actor.component != ECM_COMPONENT:
            raise ValueError("ECM state owner and actor component must both be 'ecm'")
        _validate_device_array(
            self.position_d,
            label="ecm.position_d",
            dtype=wp.vec3d,
            ndim=1,
        )
        _validate_device_array(
            self.force_d,
            label="ecm.force_d",
            dtype=wp.vec3d,
            ndim=1,
        )
        if self.position_d.shape != self.force_d.shape:
            raise ValueError("ECM position and force arrays must have identical shape")
        _validate_device_array(
            self.segments_d,
            label="ecm.segments_d",
            dtype=wp.int32,
            ndim=2,
        )
        if int(self.segments_d.shape[1]) != 2:
            raise ValueError("ECM segments_d must have shape (n_segments, 2)")
        _validate_device_array(
            self.fiber_offsets_d,
            label="ecm.fiber_offsets_d",
            dtype=wp.int32,
            ndim=1,
        )
        if int(self.fiber_offsets_d.shape[0]) < 2:
            raise ValueError("ECM fiber_offsets_d must describe at least one collagen fiber")
        n_fibers = int(self.fiber_offsets_d.shape[0]) - 1
        n_segments = int(self.segments_d.shape[0])
        per_fiber = (("fiber_sleep_state_d", self.fiber_sleep_state_d, wp.int32),)
        per_segment = (
            ("segment_refinement_level_d", self.segment_refinement_level_d, wp.int32),
            ("segment_damage_d", self.segment_damage_d, wp.float64),
            ("endpoint_refcount_d", self.endpoint_refcount_d, wp.int32),
        )
        for label, array, dtype in per_fiber:
            _validate_device_array(array, label=f"ecm.{label}", dtype=dtype, ndim=1)
            if array.shape != (n_fibers,):
                raise ValueError(f"ecm.{label} must have one entry per collagen fiber")
        for label, array, dtype in per_segment:
            _validate_device_array(array, label=f"ecm.{label}", dtype=dtype, ndim=1)
            if array.shape != (n_segments,):
                raise ValueError(f"ecm.{label} must have one entry per collagen segment")
        _validate_device_array(
            self.topology_epoch_d,
            label="ecm.topology_epoch_d",
            dtype=wp.int32,
            ndim=1,
        )
        if self.topology_epoch_d.shape != (1,):
            raise ValueError("ecm.topology_epoch_d must contain exactly one accepted topology epoch")
        if self.actor.n_elements != n_segments:
            raise ValueError("ECM actor element count must equal allocated collagen segments")
        arrays = (
            self.position_d,
            self.force_d,
            self.segments_d,
            self.fiber_offsets_d,
            self.fiber_sleep_state_d,
            self.segment_refinement_level_d,
            self.segment_damage_d,
            self.endpoint_refcount_d,
            self.topology_epoch_d,
        )
        _require_same_device(arrays, label="ECM state")
        _require_distinct_storage(arrays, label="ECM state")
        n_nodes = getattr(self.mechanics, "n_nodes", None)
        if n_nodes is not None and int(n_nodes) != int(self.position_d.shape[0]):
            raise ValueError(
                f"ECM mechanics expects {int(n_nodes)} nodes but state owns "
                f"{int(self.position_d.shape[0])}"
            )

    @property
    def device(self) -> str:
        """CUDA device string used by this state owner."""
        return str(self.position_d.device)

    @property
    def n_nodes(self) -> int:
        """Allocated collagen nodes from host-visible array metadata."""
        return int(self.position_d.shape[0])

    @property
    def n_segments(self) -> int:
        """Allocated collagen segments from host-visible array metadata."""
        return int(self.segments_d.shape[0])

    @property
    def n_fibers(self) -> int:
        """Allocated collagen fiber identities from host-visible offset metadata."""
        return int(self.fiber_offsets_d.shape[0]) - 1

    def accumulate(self) -> None:
        """Launch collagen stretch/bend/excluded-volume mechanics on owned arrays."""
        self.mechanics.accumulate(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        """Snapshot damage, sleep, refinement, and topology state."""
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Restore all rejected ECM candidate mutations."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit remodelling/damage/refinement only under the final predicate.

        After the accepted-step transaction commit — which is the only place a remesh can change the
        committed segment rest lengths / active set and bump ``topology_epoch_d`` — the mechanics
        delegate's derived per-element parameters are re-derived from the *new* committed topology via
        :meth:`refresh_mechanics_parameters`.  This never happens inside :meth:`accumulate` (the
        accepted-step contract keeps candidate force pure) and reads no authoritative device value on
        the host.
        """
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)
        self.refresh_mechanics_parameters()

    def refresh_mechanics_parameters(self) -> None:
        """Re-derive the constitutive delegate's committed-topology parameters (``k_seg`` / ``α``).

        The ECM mechanics delegate (:class:`aleph.engine.ecm_mechanics.ECMConstitutiveForce`) owns
        derived per-element parameters — ``k_seg = EA/seg_rest`` and ``α = κ/seg³`` — keyed off the
        committed segment rest lengths.  A committed remesh bumps ``topology_epoch_d`` and changes those
        rest lengths, so the parameters must be re-derived from the new committed topology.  Invoked
        AFTER an accepted-step commit (never inside :meth:`accumulate`), fully device-resident (no host
        epoch readback), and a value-idempotent no-op when the committed topology did not actually
        change.  Delegates that carry no ``refresh_parameters`` hook (the build-once cortex/SF style, or
        a test spy) are left untouched — this is additive by construction and never a new requirement on
        the :class:`~aleph.engine.runtime.MechanicsContributor` protocol.
        """
        refresh = getattr(self.mechanics, "refresh_parameters", None)
        if callable(refresh):
            refresh()

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward ECM elastic/damage/population/topology accounting."""
        self.ledger.accumulate_ledger(ledger)

    def collagen_endpoint_view(self) -> CollagenMaterialPointEndpointView:
        """Return a non-owning, generation-aware view used by the FA clutch."""
        return CollagenMaterialPointEndpointView(self)


@dataclass(frozen=True, slots=True)
class CollagenMaterialPointEndpointView:
    """Non-owning collagen-segment endpoint view for alpha2beta1 ligation."""

    owner: ECMStateOwner

    @property
    def component(self) -> str:
        """Component name required by the connector graph."""
        return ECM_COMPONENT

    @property
    def role(self) -> EndpointRole:
        """Typed load-path role for collagen ligation."""
        return EndpointRole.COLLAGEN_LIGAND

    @property
    def position_d(self) -> wp.array:
        """Borrow the live collagen positions without acquiring ownership."""
        return self.owner.position_d

    @property
    def force_d(self) -> wp.array:
        """Borrow the ECM force accumulator used by adjoint scatter."""
        return self.owner.force_d

    @property
    def segments_d(self) -> wp.array:
        """Borrow the live segment topology."""
        return self.owner.segments_d

    @property
    def topology_epoch_d(self) -> wp.array:
        """Device epoch against which connector endpoint generations are validated."""
        return self.owner.topology_epoch_d

    @property
    def endpoint_refcount_d(self) -> wp.array:
        """Borrow the conservative sleep/refinement locks for live connector endpoints."""
        return self.owner.endpoint_refcount_d

    def port_ref(self, segment_id: int, u: float) -> PortRef:
        """Create a build-time typed port without reading device geometry.

        The returned barycentric coordinate is uploaded by the graph-owned joint runtime.  Any
        remesh must update the actor/entity generation only in the same accepted transaction that
        increments ``topology_epoch_d`` and remaps live connector endpoints.
        """
        if isinstance(segment_id, bool) or not isinstance(segment_id, int):
            raise TypeError("segment_id must be an integer")
        if not 0 <= segment_id < self.owner.n_segments:
            raise ValueError("segment_id is outside the allocated collagen topology")
        if isinstance(u, bool) or not isinstance(u, int | float) or not math.isfinite(u):
            raise ValueError("collagen material coordinate must be finite")
        if not 0.0 <= float(u) <= 1.0:
            raise ValueError("collagen material coordinate must be in [0, 1]")
        actor = self.owner.actor
        return PortRef(
            component=ECM_COMPONENT,
            actor_id=actor.actor_id,
            actor_generation=actor.actor_generation,
            entity_id=actor.entity_id,
            entity_generation=actor.entity_generation,
            element_kind=ElementKind.SEGMENT,
            element_id=segment_id,
            local_coordinates=(float(u), 0.0, 0.0, 0.0),
            role=EndpointRole.COLLAGEN_LIGAND,
        )


@dataclass(frozen=True, slots=True)
class ECMWorldStepBindings:
    """Non-owning connector and boundary bindings for one ECM candidate."""

    alpha2beta1_clutch: CompositeECMClutch
    crosslink_connector: ECMInternalCrosslink
    boundary_anchor: ECMBoundaryAnchorFacade
    membrane: MembraneContactSurfaceView
    membrane_contact: ECMMembraneContact

    def __post_init__(self) -> None:
        clutch = self.alpha2beta1_clutch
        if clutch.mechanical_group != FA_SERIES_GROUP:
            raise ValueError(f"ECM clutch must resolve mechanical group {FA_SERIES_GROUP!r}")
        if not (
            clutch.equal_opposite_required
            and clutch.adjoint_transfer_required
            and clutch.work_ledger_required
        ):
            raise ValueError("ECM clutch requires equal/opposite adjoint mechanics and work ledger")
        crosslink = self.crosslink_connector
        if crosslink.name != ECM_CROSSLINK_CONNECTOR:
            raise ValueError(f"expected graph connector {ECM_CROSSLINK_CONNECTOR!r}")
        if {crosslink.component_a, crosslink.component_b} != {ECM_COMPONENT}:
            raise ValueError("ECM crosslink connector must remain internal to ecm")
        if not crosslink.adjoint_transfer_required:
            raise ValueError("ECM crosslink connector requires an adjoint force scatter")
        if self.boundary_anchor.component != ECM_COMPONENT:
            raise ValueError("ECM boundary anchor must address ecm")
        if self.membrane_contact.name != MEMBRANE_ECM_CONTACT:
            raise ValueError(f"expected graph connector {MEMBRANE_ECM_CONTACT!r}")
        if {
            self.membrane_contact.component_a,
            self.membrane_contact.component_b,
        } != {MEMBRANE_COMPONENT, ECM_COMPONENT}:
            raise ValueError("membrane ECM contact endpoints must be membrane and ecm")
        if not self.membrane_contact.adjoint_transfer_required:
            raise ValueError("membrane ECM contact requires an adjoint force scatter")


@dataclass(frozen=True, slots=True)
class ECMWorldDofLedger:
    """Allocation metadata; active/sleeping populations remain device-ledger values."""

    allocated_nodes: int
    allocated_segments: int
    allocated_fibers: int
    sleep_state_slots: int
    refinement_state_slots: int
    damage_state_slots: int
    endpoint_refcount_slots: int
    topology_epoch_slots: int


def _find_connector(architecture: CellArchitecture, name: str) -> object:
    """Return a named connector or raise a targeted build error."""
    for connector in architecture.connectors:
        if connector.name == name:
            return connector
    raise ValueError(f"architecture must register {name!r}")


def _validate_architecture(architecture: CellArchitecture) -> None:
    """Validate the ECM component and the two graph-owned connector contracts."""
    ecm = architecture.component(ECM_COMPONENT)
    if ecm.role is not ComponentRole.ENVIRONMENT or not ecm.owns_geometry:
        raise ValueError("registered ecm must be a geometry-owning ENVIRONMENT component")

    resolution = resolve_fa_series_group(architecture)
    if resolution.name != FA_SERIES_GROUP or resolution.mechanical_joint_count != 1:
        raise ValueError("FA semantic edges must resolve to one alpha2beta1 collagen joint")
    clutch = _find_connector(architecture, FA_ECM_CLUTCH_CONNECTOR)
    if clutch.family is not ConnectorFamily.FA_CLUTCH:
        raise ValueError("integrin_collagen_clutch must use ConnectorFamily.FA_CLUTCH")
    if {clutch.component_a, clutch.component_b} != {"focal_adhesion", ECM_COMPONENT}:
        raise ValueError("integrin_collagen_clutch must join focal_adhesion and ecm")
    if clutch.chemistry_card != "alpha2beta1_collagen":
        raise ValueError("ECM clutch must use alpha2beta1_collagen chemistry")
    if not clutch.kinetics or not clutch.commit_on_accept:
        raise ValueError("ECM clutch kinetics must commit only on accepted steps")
    if not clutch.bidirectional or not clutch.adjoint_transfer_required:
        raise ValueError("ECM clutch must be bidirectional with adjoint transfer")

    crosslink = _find_connector(architecture, ECM_CROSSLINK_CONNECTOR)
    if crosslink.family is not ConnectorFamily.FIBER_CROSSLINK:
        raise ValueError("ecm_crosslink must use ConnectorFamily.FIBER_CROSSLINK")
    if crosslink.scope is not ConnectorScope.INTERNAL:
        raise ValueError("ecm_crosslink must be an internal connector")
    if {crosslink.component_a, crosslink.component_b} != {ECM_COMPONENT}:
        raise ValueError("ecm_crosslink must stay inside ecm")
    if crosslink.chemistry_card != "collagen_crosslink":
        raise ValueError("ecm_crosslink must use collagen_crosslink chemistry")
    if not crosslink.kinetics or not crosslink.commit_on_accept:
        raise ValueError("ECM crosslink kinetics must commit only on accepted steps")
    if not crosslink.bidirectional or not crosslink.adjoint_transfer_required:
        raise ValueError("ECM crosslinks require bidirectional adjoint mechanics")
    if not (
        crosslink.generation_required
        and crosslink.remap_on_accept
        and crosslink.blocks_sleep_refine
    ):
        raise ValueError("ECM crosslinks must lock and remap live topology generations")

    boundary = _find_connector(architecture, ECM_FAR_FIELD_ANCHOR)
    if boundary.family is not ConnectorFamily.ENVIRONMENT_BOUNDARY:
        raise ValueError("ecm_far_field_anchor must use ENVIRONMENT_BOUNDARY")
    if {boundary.component_a, boundary.component_b} != {ECM_COMPONENT, "world_boundary"}:
        raise ValueError("ECM far-field anchor must join ecm and world_boundary")
    if boundary.kinetics or boundary.commit_on_accept:
        raise ValueError("ECM far-field anchor is a non-kinetic boundary constraint")
    if not boundary.bidirectional or not boundary.adjoint_transfer_required:
        raise ValueError("ECM far-field anchor must expose reaction and adjoint work")

    contact = _find_connector(architecture, MEMBRANE_ECM_CONTACT)
    if contact.family is not ConnectorFamily.CONTACT:
        raise ValueError("membrane_ecm_contact must use ConnectorFamily.CONTACT")
    if {contact.component_a, contact.component_b} != {MEMBRANE_COMPONENT, ECM_COMPONENT}:
        raise ValueError("membrane_ecm_contact must join membrane and ecm")
    if contact.kinetics or contact.commit_on_accept:
        raise ValueError("membrane_ecm_contact must be a non-kinetic contact")
    if not contact.bidirectional or not contact.adjoint_transfer_required:
        raise ValueError("membrane_ecm_contact requires bidirectional adjoint mechanics")
    if not (
        contact.generation_required
        and contact.remap_on_accept
        and contact.blocks_sleep_refine
    ):
        raise ValueError("membrane_ecm_contact must lock and remap ECM topology generations")


@dataclass(frozen=True, slots=True)
class ECMWorld:
    """Coupled facade over live collagen state, graph connectors, and far-field closure."""

    architecture: CellArchitecture
    ecm: ECMStateOwner

    def __post_init__(self) -> None:
        _validate_architecture(self.architecture)
        if self.ecm.name != ECM_COMPONENT:
            raise ValueError("ECMWorld.ecm must be the registered ECM state owner")

    def _validate_bindings(self, bindings: ECMWorldStepBindings) -> None:
        """Reject a clutch endpoint borrowed from a different ECM actor."""
        if bindings.alpha2beta1_clutch.ecm_endpoint.owner is not self.ecm:
            raise ValueError("alpha2beta1 clutch endpoint must borrow this ECM state owner")
        if bindings.boundary_anchor.mode is not self.ecm.settings.anchor_mode:
            raise ValueError("ECM boundary mode must match the physiological world settings")

    def accumulate_mechanics(self, bindings: ECMWorldStepBindings) -> None:
        """Assemble collagen, crosslink, contact, and clutch, then close the far-field boundary LAST.

        The caller owns force zeroing and ensures this facade is the sole dispatch site for these
        bindings.  The clutch runtime scatters to both its SF and ECM actor arrays exactly once.

        Ordering is load-bearing, not stylistic.  ``boundary_anchor.accumulate`` runs
        :func:`_far_field_reaction_ledger_kernel`, which records ``reaction = -f`` at each pinned node
        and then ZEROES that node's force.  So it must run after every force contributor, and it used
        to run third — before ``membrane_contact`` and ``alpha2beta1_clutch``.  That had two distinct
        consequences, both fixed by moving it last:

        1. The recorded far-field reaction structurally EXCLUDED the clutch and membrane-contact
           contributions, so any traction read off this ledger was wrong by exactly the two channels
           that carry the cell's grip on the ECM — the ones a traction measurement is about.
        2. Worse, and separately: because the kernel zeroes the pinned node's force, the two channels
           that ran afterwards re-loaded those nodes and left them with a NET FORCE. The Dirichlet
           far-field pin was silently violated for precisely those two channels.
        """
        self._validate_bindings(bindings)
        self.ecm.accumulate()
        bindings.crosslink_connector.accumulate_internal(
            self.ecm.position_d,
            self.ecm.force_d,
        )
        bindings.membrane_contact.accumulate_contact(
            bindings.membrane,
            self.ecm.collagen_endpoint_view(),
        )
        bindings.alpha2beta1_clutch.accumulate()
        # LAST — records the reaction and zeroes the pin; anything after it would escape both.
        bindings.boundary_anchor.accumulate(self.ecm.position_d, self.ecm.force_d)

    def transaction_participants(
        self,
        bindings: ECMWorldStepBindings,
    ) -> tuple[
        ECMStateOwner,
        ECMInternalCrosslink,
        ECMBoundaryAnchorFacade,
        ECMMembraneContact,
        CompositeECMClutch,
    ]:
        """Return every owner/cache participating in accepted-step atomicity."""
        self._validate_bindings(bindings)
        return (
            self.ecm,
            bindings.crosslink_connector,
            bindings.boundary_anchor,
            bindings.membrane_contact,
            bindings.alpha2beta1_clutch,
        )

    def snapshot_candidate(self, bindings: ECMWorldStepBindings) -> None:
        """Snapshot ECM, crosslink, boundary, and clutch state as one candidate."""
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: ECMWorldStepBindings, accepted: wp.array) -> None:
        """Forward one device acceptance predicate to every ECM participant."""
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self,
        bindings: ECMWorldStepBindings,
        accepted: wp.array,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        """Atomically finalize topology, damage, remodelling, boundary, and clutch state."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if isinstance(rng_seed, bool) or not isinstance(rng_seed, int) or rng_seed < 0:
            raise ValueError("rng_seed must be a nonnegative integer")
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, bindings: ECMWorldStepBindings, ledger: object) -> None:
        """Collect component, connector, and boundary ledgers without deciding acceptance."""
        for participant in self.transaction_participants(bindings):
            participant.accumulate_ledger(ledger)

    def dof_ledger(self) -> ECMWorldDofLedger:
        """Return immutable allocation counts without reading active device state."""
        return ECMWorldDofLedger(
            allocated_nodes=self.ecm.n_nodes,
            allocated_segments=self.ecm.n_segments,
            allocated_fibers=self.ecm.n_fibers,
            sleep_state_slots=int(self.ecm.fiber_sleep_state_d.shape[0]),
            refinement_state_slots=int(self.ecm.segment_refinement_level_d.shape[0]),
            damage_state_slots=int(self.ecm.segment_damage_d.shape[0]),
            endpoint_refcount_slots=int(self.ecm.endpoint_refcount_d.shape[0]),
            topology_epoch_slots=int(self.ecm.topology_epoch_d.shape[0]),
        )


class _LegacyLocalMikadoMechanics(Protocol):
    """Narrow subset allowed through the legacy mechanics adapter."""

    node_off: int
    n_nodes: int
    owns_crosslinks: bool
    owns_boundary_anchor: bool

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch local collagen mechanics on caller-owned CUDA arrays."""


@dataclass(frozen=True, slots=True)
class LegacyLocalMikadoMechanicsAdapter:
    """Reuse only local collagen mechanics from a legacy compartment.

    The wrapped object must address a local array from offset zero and must not retain crosslink or
    boundary state.  Host KD-tree construction, host relaxation loops, fixed-node clutch endpoints,
    and hidden pinning are intentionally outside the adapter boundary.
    """

    compartment: _LegacyLocalMikadoMechanics

    def __post_init__(self) -> None:
        if int(self.compartment.node_off) != 0:
            raise ValueError("legacy ECM adapter requires a local array with node_off == 0")
        if int(self.compartment.n_nodes) <= 0:
            raise ValueError("legacy ECM compartment must contain collagen nodes")
        if self.compartment.owns_crosslinks:
            raise ValueError("legacy ECM adapter cannot own graph crosslink state")
        if self.compartment.owns_boundary_anchor:
            raise ValueError("legacy ECM adapter cannot hide the far-field boundary anchor")

    @property
    def n_nodes(self) -> int:
        """Number of collagen nodes required by the wrapped mechanics."""
        return int(self.compartment.n_nodes)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Delegate existing Warp collagen force kernels."""
        self.compartment.accumulate(pos, force)


def build_ecm_state_owner(
    topology: "ECMTopologyState",
    *,
    settings: ECMWorldSettings,
    actor: ActorRecord,
    transaction: TransactionParticipant,
    ledger: LedgerContributor,
    spec: "ECMSpec | str | None" = None,
    device: str | None = None,
    launch: Callable[..., object] = wp.launch,
    force_builder: Callable[..., "ECMConstitutiveForce"] | None = None,
) -> ECMStateOwner:
    """Assemble the ``ecm`` state owner with a REAL collagen constitutive force pass as its mechanics.

    This is the engine-side construction site the ECM owner was missing: it builds the concrete
    :class:`aleph.engine.ecm_mechanics.ECMConstitutiveForce` from the owner's on-GPU
    :class:`~aleph.components.ecm.device_schema.ECMTopologyState` SoA + the SOURCED collagen card
    (``build_collagen_constitutive_force`` — which allocates the owned ``k_seg_d`` / ``alpha_d`` derived
    parameter arrays, binds the SoA connectivity directly, and derives the parameters once via
    ``refresh_parameters``), then wires that force pass into :attr:`ECMStateOwner.mechanics`.  The owner's
    per-step force hook (:meth:`ECMStateOwner.accumulate`) therefore drives the collagen stretch/bend
    kernels under OWNER control, and :meth:`ECMStateOwner.commit_irreversible` re-derives the parameters
    after an accepted remesh — moving ``ecm`` from KERNEL_BOUND toward owner-driven CONNECTED.

    The ECM is engine-only (no incumbent driver computes collagen constitutive forces), so this is a
    GENUINE binding of net-new mechanistic physics, not a facade over a lumped incumbent.  The
    non-mechanics delegates (``actor`` / ``transaction`` / ``ledger``) stay separately owned exactly as
    the ECM transaction contract requires — this builder only fixes the mechanics seam.

    Args:
        topology: the on-GPU collagen Mikado SoA whose node/segment arrays the owner will own.
        settings: the physiological ECM world switches (rejects a non-physiological build).
        actor: the ``ecm`` actor record (``n_elements`` must equal the topology segment capacity).
        transaction: the owner's accepted-step transaction delegate (damage/sleep/refinement/remesh).
        ledger: the owner's ledger delegate.
        spec: an :class:`~aleph.laws.ecm_library.ECMSpec` or a registry key; defaults to ``collagen_I``.
        device: the CUDA device the topology lives on (defaults to the topology's device).
        launch: injectable Warp launcher (defaults to ``wp.launch``).
        force_builder: injectable force-pass factory (defaults to
            ``ecm_mechanics.build_collagen_constitutive_force``); the CUDA-free structural gate injects a
            recording double so the array wiring can be asserted off-GPU.

    Returns:
        The assembled :class:`ECMStateOwner` whose ``mechanics`` is the bound constitutive force pass.
    """
    from aleph.engine.ecm_mechanics import (
        COLLAGEN_MATERIAL_KEY,
        build_collagen_constitutive_force,
    )

    builder = force_builder or build_collagen_constitutive_force
    card = COLLAGEN_MATERIAL_KEY if spec is None else spec
    force_pass = builder(topology, card, device=device, launch=launch)
    return ECMStateOwner(
        name=ECM_COMPONENT,
        settings=settings,
        actor=actor,
        position_d=topology.position_d,
        force_d=topology.force_d,
        segments_d=topology.segments_d,
        fiber_offsets_d=topology.fiber_offsets_d,
        fiber_sleep_state_d=topology.fiber_sleep_state_d,
        segment_refinement_level_d=topology.segment_refinement_level_d,
        segment_damage_d=topology.segment_damage_d,
        endpoint_refcount_d=topology.endpoint_refcount_d,
        topology_epoch_d=topology.topology_epoch_d,
        mechanics=force_pass,
        transaction=transaction,
        ledger=ledger,
    )


__all__ = [
    "BoundaryAnchorDelegate",
    "BoundaryAnchorMode",
    "CollagenMaterialPointEndpointView",
    "CompositeLoadPathECMClutchAdapter",
    "CompositeECMClutch",
    "ECMBoundaryAnchorFacade",
    "ECM_COMPONENT",
    "ECM_CROSSLINK_CONNECTOR",
    "ECM_FAR_FIELD_ANCHOR",
    "ECMInternalCrosslink",
    "ECMMembraneContact",
    "ECMStateOwner",
    "ECMWorld",
    "ECMWorldDofLedger",
    "ECMWorldSettings",
    "ECMWorldStepBindings",
    "FarFieldDirichletBoundaryRuntime",
    "FA_ECM_CLUTCH_CONNECTOR",
    "FA_SERIES_GROUP",
    "LegacyLocalMikadoMechanicsAdapter",
    "MEMBRANE_COMPONENT",
    "MEMBRANE_ECM_CONTACT",
    "MembraneContactSurfaceView",
    "build_ecm_state_owner",
]
