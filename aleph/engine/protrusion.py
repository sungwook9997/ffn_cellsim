"""Non-invasive runtime seam for explicit lamellipodium and filopodium actors.

The two protrusions are distinct state owners.  A lamellipodium owns an adaptive,
explicit Arp2/3 mother/daughter graph; a filopodium owns an adaptive, explicit bundled
actin graph.  Neither may be replaced by a lumped active-stress patch.

This first slice fixes CUDA ownership, connector bindings, accepted-step event coverage,
and refinement ledgers.  Existing weave builders remain initialization/oracle assets;
no legacy runtime is mutated here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentContract,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
)
from aleph.engine.runtime import (
    LedgerContributor,
    MechanicsContributor,
    TransactionParticipant,
)
from aleph.components.weave.branch_angle import ARP23_K_THETA, ARP23_THETA0_RAD
from aleph.components.weave.branch_angle_warp import branch_angle_kernel

LAMELLIPODIUM_COMPONENT = "lamellipodium"
FILOPODIUM_COMPONENT = "filopodium"
MEMBRANE_COMPONENT = "membrane"
CORTEX_COMPONENT = "cortex"
CYTOSOL_COMPONENT = "cytosol"
FA_COMPONENT = "focal_adhesion"

LAMELLIPODIUM_REPRESENTATION = "local adaptive explicit Arp2/3 branched F-actin"
FILOPODIUM_REPRESENTATION = "local adaptive explicit bundled F-actin"

LAM_MEMBRANE = "lamellipodium_membrane_contact"
LAM_CORTEX = "lamellipodium_cortex_seam"
LAM_CYTOSOL = "lamellipodium_cytosol_transfer"
LAM_NASCENT_FA = "lamellipodium_nascent_fa"
FILO_MEMBRANE = "filopodium_membrane_tip"
FILO_CORTEX = "filopodium_cortex_root"
FILO_CYTOSOL = "filopodium_cytosol_transfer"
FILO_NASCENT_FA = "filopodium_nascent_fa"

LAMELLIPODIUM_EVENT_CHANNELS = frozenset(
    {"branching", "polymerization", "capping", "severing"}
)
FILOPODIUM_EVENT_CHANNELS = frozenset(
    {"bundling", "polymerization", "capping", "severing"}
)


def _storage_key(array: object) -> tuple[object, ...]:
    """Return device-storage identity without reading authoritative data."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


def _validate_device_array(
    array: object,
    *,
    label: str,
    dtype: object,
    rank: int = 1,
    shape0: int | None = None,
    shape1: int | None = None,
) -> None:
    """Validate CUDA residence, dtype, rank and selected extents without a readback."""
    device = getattr(array, "device", None)
    if not bool(getattr(device, "is_cuda", False)):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != dtype:
        raise TypeError(f"{label} must have dtype {dtype}")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != rank:
        raise ValueError(f"{label} must have rank {rank}")
    if any(int(extent) <= 0 for extent in shape):
        raise ValueError(f"{label} must be non-empty")
    if shape0 is not None and int(shape[0]) != shape0:
        raise ValueError(f"{label} must have leading extent {shape0}")
    if shape1 is not None and int(shape[1]) != shape1:
        raise ValueError(f"{label} must have trailing extent {shape1}")


@runtime_checkable
class ProtrusionTransaction(TransactionParticipant, Protocol):
    """Transaction delegate covering explicit protrusion topology and event channels."""

    component_name: str
    event_channels: frozenset[str]

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return every authoritative mutable array covered by snapshot/rollback."""


@dataclass(frozen=True, slots=True)
class ProtrusionActorView:
    """Non-owning geometry/topology view exposed to graph-owned connectors."""

    component: str
    representation: str
    position_d: wp.array
    force_d: wp.array
    fiber_offset_d: wp.array
    active_filament_d: wp.array
    active_node_d: wp.array
    persistent_filament_id_d: wp.array
    refinement_parent_id_d: wp.array
    refinement_level_d: wp.array
    actor_generation_d: wp.array
    topology_epoch_d: wp.array
    n_filament_capacity: int


@runtime_checkable
class ProtrusionGraphConnector(Protocol):
    """Graph-owned coupling from one protrusion actor to an external component."""

    name: str
    component_a: str
    component_b: str

    def accumulate_actor(self, actor: ProtrusionActorView) -> None:
        """Gather live actor geometry and scatter an adjoint external reaction."""

    def snapshot_candidate(self) -> None:
        """Snapshot connector state/caches for the candidate physical step."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore connector state under a rejected device predicate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit binding/remodelling only under final acceptance."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Add connector force/work/mass/topology terms without deciding acceptance."""


@dataclass(frozen=True, slots=True)
class RefinementClosureReport:
    """A posteriori fine<->coarse mapping closure for one protrusion actor."""

    identity_mismatch_count: int
    topology_mismatch_count: int
    force_residual_pn: float
    work_residual_pn_um: float

    def assert_closed(self, *, force_tolerance_pn: float, work_tolerance_pn_um: float) -> None:
        """Reject a mapping that changes identity, topology, force, or work."""
        if force_tolerance_pn < 0.0 or work_tolerance_pn_um < 0.0:
            raise ValueError("mapping tolerances must be nonnegative")
        if self.identity_mismatch_count != 0:
            raise ValueError("refinement mapping changed persistent filament identity")
        if self.topology_mismatch_count != 0:
            raise ValueError("refinement mapping changed active graph topology")
        if abs(self.force_residual_pn) > force_tolerance_pn:
            raise ValueError("refinement mapping does not close resultant force")
        if abs(self.work_residual_pn_um) > work_tolerance_pn_um:
            raise ValueError("refinement mapping does not close virtual work")


def _validate_common_actor_arrays(
    *,
    component: str,
    n_filament_capacity: int,
    position_d: object,
    force_d: object,
    fiber_offset_d: object,
    active_filament_d: object,
    active_node_d: object,
    persistent_filament_id_d: object,
    refinement_parent_id_d: object,
    refinement_level_d: object,
    actor_generation_d: object,
    topology_epoch_d: object,
) -> None:
    """Validate the shared CUDA ownership metadata of either protrusion actor."""
    if n_filament_capacity <= 0:
        raise ValueError(f"{component} filament capacity must be positive")
    _validate_device_array(position_d, label=f"{component}.position_d", dtype=wp.vec3d)
    n_nodes = int(position_d.shape[0])
    _validate_device_array(
        force_d,
        label=f"{component}.force_d",
        dtype=wp.vec3d,
        shape0=n_nodes,
    )
    _validate_device_array(
        fiber_offset_d,
        label=f"{component}.fiber_offset_d",
        dtype=wp.int32,
        shape0=n_filament_capacity + 1,
    )
    for label, array, dtype, length in (
        ("active_filament_d", active_filament_d, wp.int32, n_filament_capacity),
        ("active_node_d", active_node_d, wp.int32, n_nodes),
        ("persistent_filament_id_d", persistent_filament_id_d, wp.int64, n_filament_capacity),
        ("refinement_parent_id_d", refinement_parent_id_d, wp.int64, n_filament_capacity),
        ("refinement_level_d", refinement_level_d, wp.int32, n_filament_capacity),
        ("actor_generation_d", actor_generation_d, wp.int32, 1),
        ("topology_epoch_d", topology_epoch_d, wp.int32, 1),
    ):
        _validate_device_array(
            array,
            label=f"{component}.{label}",
            dtype=dtype,
            shape0=length,
        )
    arrays = (
        position_d,
        force_d,
        fiber_offset_d,
        active_filament_d,
        active_node_d,
        persistent_filament_id_d,
        refinement_parent_id_d,
        refinement_level_d,
        actor_generation_d,
        topology_epoch_d,
    )
    if len({str(array.device) for array in arrays}) != 1:
        raise ValueError(f"all {component} arrays must share one CUDA device")
    if _storage_key(position_d) == _storage_key(force_d):
        raise ValueError(f"{component} position and force arrays must not alias")


def _validate_transaction_coverage(
    transaction: ProtrusionTransaction,
    *,
    component: str,
    required_channels: frozenset[str],
    mutable_arrays: tuple[object, ...],
) -> None:
    """Ensure rollback covers geometry/topology and all required accepted events."""
    if transaction.component_name != component:
        raise ValueError(f"{component} transaction has the wrong component owner")
    missing_channels = required_channels - transaction.event_channels
    if missing_channels:
        raise ValueError(f"{component} transaction is missing event channels {sorted(missing_channels)}")
    covered = {_storage_key(array) for array in transaction.owned_arrays()}
    missing_arrays = {_storage_key(array) for array in mutable_arrays} - covered
    if missing_arrays:
        raise ValueError(f"{component} transaction does not cover all mutable topology arrays")


# --- Faessler 2020 EMBO J 39:e104254 sourced Arp2/3 branch anchors (KB weave/params_i0b4) ---
# theta0 = 70 deg rest branch angle (in-cell cryo-ET 68 +/- 9 deg); k_theta derived by equipartition
# k_theta = kT/sigma_theta^2 = 0.173 pN*um/rad^2. These are the ONLY admissible branch-angle constants;
# any other value is a magic-number substitution and is rejected at construction.
LAM_BRANCH_THETA0_RAD = ARP23_THETA0_RAD
LAM_BRANCH_K_THETA_PN_UM = ARP23_K_THETA

# --- Filopodium fascin-bundle geometry (FOUND) + crosslink stiffness (NOT FOUND, PI-GAP) ---
# PI_GAP_LITERATURE_SOURCING_2026-07-23 §4.  The bundle GEOMETRY is sourced (Courson & Rock 2010 JBC,
# 10.1074/jbc.M110.123117): strictly parallel filaments at ~8 nm inter-filament spacing, ~36-38 nm axial
# repeat.  No published SINGLE-MOLECULE fascin crosslink stiffness (pN/µm) exists anywhere; do NOT inherit the
# filamin 8.2e5 value.  So the rest spacing defaults to the sourced geometry, but ``k_fascin_pn_per_um`` is an
# UNSET slot with NO default — the bundle mechanism binds now and activates only when PI supplies a stiffness.
FASCIN_INTERFILAMENT_SPACING_UM = 0.008   # ~8 nm, Courson & Rock 2010 JBC 10.1074/jbc.M110.123117
FASCIN_AXIAL_REPEAT_UM = 0.036            # ~36-38 nm crossover repeat, same
FASCIN_K_CROSSLINK_PN_PER_UM: float | None = None   # NOT FOUND — PI-GAP slot, no default


@runtime_checkable
class KernelLauncher(Protocol):
    """Injectable Warp launch seam (``wp.launch`` in production; a capturing double in CPU tests)."""

    def __call__(
        self,
        kernel: object,
        *,
        dim: int,
        inputs: list[object],
        outputs: list[object],
        device: object,
    ) -> None:
        """Enqueue one kernel launch over ``dim`` threads against injected device arrays."""


def _warp_launch(
    kernel: object,
    *,
    dim: int,
    inputs: list[object],
    outputs: list[object],
    device: object,
) -> None:
    """Default launcher delegating to ``wp.launch`` (only ever executed on a CUDA device)."""
    wp.launch(kernel, dim=dim, inputs=inputs, outputs=outputs, device=device)


@dataclass(frozen=True, slots=True)
class LamellipodiumBranchAngleMechanics:
    """KERNEL_BOUND delegate: launches the real Arp2/3 angle-harmonic Warp kernel through the seam.

    This binds ``aleph.components.weave.branch_angle_warp.branch_angle_kernel`` (a Warp port of the ff
    ``network_warp`` branch-angle force whose CPU oracle is ``ac.weave.branch_angle``).  The delegate holds
    NO authoritative force/topology storage of its own: it reads the actor-owned branch triples/activity and
    scatters the harmonic-angle force directly into the actor-owned ``force`` array passed by the seam, so
    there is no private state path (the KERNEL_BOUND requirement).  The three nodal forces of each triple sum
    to zero inside the kernel, so momentum is conserved and Newton's third law holds per branch junction.

    Constants are the Faessler 2020 anchors only; an off-source theta0/k_theta is rejected (no magic numbers).
    The ``launch`` seam defaults to ``wp.launch`` and is injected with a capturing double under CPU tests
    (there is no CUDA on the dev Mac; the force/sign/work CUDA_UNIT gate runs on the A5000 in the GPU lane).
    """

    n_branch: int
    branch_triples_d: wp.array
    branch_active_d: wp.array
    theta0_rad: float = LAM_BRANCH_THETA0_RAD
    k_theta_pn_um: float = LAM_BRANCH_K_THETA_PN_UM
    kernel: object = branch_angle_kernel
    launch: KernelLauncher = _warp_launch

    def __post_init__(self) -> None:
        if self.n_branch <= 0:
            raise ValueError("lamellipodium branch mechanics needs a positive branch count")
        if abs(self.theta0_rad - LAM_BRANCH_THETA0_RAD) > 1.0e-12:
            raise ValueError("branch theta0 must be the sourced Faessler 2020 anchor (no magic numbers)")
        if abs(self.k_theta_pn_um - LAM_BRANCH_K_THETA_PN_UM) > 1.0e-12:
            raise ValueError(
                "branch k_theta must be the equipartition-derived Faessler 2020 anchor (no magic numbers)"
            )
        if self.kernel is not branch_angle_kernel:
            raise ValueError("lamellipodium branch mechanics must bind the real branch_angle_kernel")
        _validate_device_array(
            self.branch_triples_d,
            label="branch_mechanics.branch_triples_d",
            dtype=wp.int32,
            rank=2,
            shape0=self.n_branch,
            shape1=3,
        )
        _validate_device_array(
            self.branch_active_d,
            label="branch_mechanics.branch_active_d",
            dtype=wp.int32,
            shape0=self.n_branch,
        )

    def assert_bound_to(
        self,
        *,
        branch_triples_d: wp.array,
        branch_active_d: wp.array,
        n_branch: int,
    ) -> None:
        """Reject a delegate that carries a private copy of the actor's branch topology arrays."""
        if self.n_branch != n_branch:
            raise ValueError("branch mechanics capacity disagrees with the lamellipodium owner")
        if _storage_key(self.branch_triples_d) != _storage_key(branch_triples_d):
            raise ValueError("branch mechanics must launch against the actor-owned branch_triples_d")
        if _storage_key(self.branch_active_d) != _storage_key(branch_active_d):
            raise ValueError("branch mechanics must launch against the actor-owned branch_active_d")

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the real angle-harmonic kernel, scattering force into the actor-owned accumulator."""
        self.launch(
            self.kernel,
            dim=self.n_branch,
            inputs=[
                pos,
                self.branch_triples_d,
                self.branch_active_d,
                self.theta0_rad,
                self.k_theta_pn_um,
            ],
            outputs=[force],
            device=getattr(pos, "device", None),
        )


@dataclass(frozen=True, slots=True)
class LamellipodiumStateOwner:
    """CUDA owner of a local adaptive explicit Arp2/3 branched network."""

    n_filament_capacity: int
    n_branch_capacity: int
    position_d: wp.array
    force_d: wp.array
    fiber_offset_d: wp.array
    active_filament_d: wp.array
    active_node_d: wp.array
    persistent_filament_id_d: wp.array
    refinement_parent_id_d: wp.array
    refinement_level_d: wp.array
    branch_triples_d: wp.array
    branch_active_d: wp.array
    barbed_state_d: wp.array
    actor_generation_d: wp.array
    topology_epoch_d: wp.array
    mechanics: MechanicsContributor
    transaction: ProtrusionTransaction
    ledger: LedgerContributor
    component: str = LAMELLIPODIUM_COMPONENT
    representation: str = LAMELLIPODIUM_REPRESENTATION

    def __post_init__(self) -> None:
        if self.component != LAMELLIPODIUM_COMPONENT:
            raise ValueError("lamellipodium owner cannot be relabelled")
        if self.representation != LAMELLIPODIUM_REPRESENTATION:
            raise ValueError("lamellipodium must remain explicit Arp2/3 branched F-actin")
        if self.n_branch_capacity <= 0:
            raise ValueError("lamellipodium branch capacity must be positive")
        _validate_common_actor_arrays(
            component=self.component,
            n_filament_capacity=self.n_filament_capacity,
            position_d=self.position_d,
            force_d=self.force_d,
            fiber_offset_d=self.fiber_offset_d,
            active_filament_d=self.active_filament_d,
            active_node_d=self.active_node_d,
            persistent_filament_id_d=self.persistent_filament_id_d,
            refinement_parent_id_d=self.refinement_parent_id_d,
            refinement_level_d=self.refinement_level_d,
            actor_generation_d=self.actor_generation_d,
            topology_epoch_d=self.topology_epoch_d,
        )
        _validate_device_array(
            self.branch_triples_d,
            label="lamellipodium.branch_triples_d",
            dtype=wp.int32,
            rank=2,
            shape0=self.n_branch_capacity,
            shape1=3,
        )
        _validate_device_array(
            self.branch_active_d,
            label="lamellipodium.branch_active_d",
            dtype=wp.int32,
            shape0=self.n_branch_capacity,
        )
        _validate_device_array(
            self.barbed_state_d,
            label="lamellipodium.barbed_state_d",
            dtype=wp.int32,
            shape0=self.n_filament_capacity,
        )
        mutable = (
            self.position_d,
            self.fiber_offset_d,
            self.active_filament_d,
            self.active_node_d,
            self.persistent_filament_id_d,
            self.refinement_parent_id_d,
            self.refinement_level_d,
            self.branch_triples_d,
            self.branch_active_d,
            self.barbed_state_d,
            self.actor_generation_d,
            self.topology_epoch_d,
        )
        if len({_storage_key(array) for array in mutable}) != len(mutable):
            raise ValueError("lamellipodium authoritative arrays must not alias")
        _validate_transaction_coverage(
            self.transaction,
            component=self.component,
            required_channels=LAMELLIPODIUM_EVENT_CHANNELS,
            mutable_arrays=mutable,
        )
        if isinstance(self.mechanics, LamellipodiumBranchAngleMechanics):
            self.mechanics.assert_bound_to(
                branch_triples_d=self.branch_triples_d,
                branch_active_d=self.branch_active_d,
                n_branch=self.n_branch_capacity,
            )

    def geometry(self) -> ProtrusionActorView:
        """Expose non-owning geometry and persistent refinement identity."""
        return ProtrusionActorView(
            component=self.component,
            representation=self.representation,
            position_d=self.position_d,
            force_d=self.force_d,
            fiber_offset_d=self.fiber_offset_d,
            active_filament_d=self.active_filament_d,
            active_node_d=self.active_node_d,
            persistent_filament_id_d=self.persistent_filament_id_d,
            refinement_parent_id_d=self.refinement_parent_id_d,
            refinement_level_d=self.refinement_level_d,
            actor_generation_d=self.actor_generation_d,
            topology_epoch_d=self.topology_epoch_d,
            n_filament_capacity=self.n_filament_capacity,
        )

    def accumulate(self) -> None:
        """Launch explicit filament/branch mechanics against actor-owned arrays."""
        self.mechanics.accumulate(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger.accumulate_ledger(ledger)


@wp.kernel
def _fascin_bundle_kernel(
    pos: wp.array(dtype=wp.vec3d),               # (N,) node positions
    pairs: wp.array(dtype=wp.int32, ndim=2),     # (C, 2) fascin-crosslinked node pairs
    active: wp.array(dtype=wp.int32),            # (C,) 1 = engaged fascin crosslink, 0 = free slot
    k_fascin: wp.float64,                        # crosslink stiffness [pN/µm] — PI-GAP slot value
    rest_spacing: wp.float64,                    # inter-filament rest spacing [µm] (sourced geometry)
    force: wp.array(dtype=wp.vec3d),             # (N,) out (atomic accumulate)
) -> None:
    """Hookean fascin crosslink holding two bundled filaments at ``rest_spacing``; only engaged slots pull.

    One thread per crosslink slot.  Equal-and-opposite ``k·(L−rest)·û`` on the pair (Newton's third law), so a
    stretched bundle is pulled tight and momentum is conserved.  This is the actin-crosslink spring form
    (``ff.network_warp.link_spring_kernel``) with an activity gate; the stiffness is the injected PI-GAP value.
    """
    t = wp.tid()
    if active[t] == 0:
        return
    i = pairs[t, 0]
    j = pairs[t, 1]
    d = pos[j] - pos[i]
    length = wp.length(d)
    if length > wp.float64(1.0e-12):
        f = (k_fascin * (length - rest_spacing) / length) * d
        wp.atomic_add(force, i, f)
        wp.atomic_add(force, j, -f)


@dataclass(frozen=True, slots=True)
class FilopodiumFascinBundleMechanics:
    """PI-GAP-gated delegate: launches the real fascin bundle-crosslink Warp kernel through the seam.

    Mirrors :class:`LamellipodiumBranchAngleMechanics` for the filopodium.  It holds NO authoritative
    force/topology storage: it reads the actor-owned crosslink pairs/activity and scatters the crosslink spring
    force into the actor-owned ``force`` array passed by the owner, so there is no private state path.  The
    bundle GEOMETRY (``rest_spacing_um``) defaults to the sourced ~8 nm Courson & Rock 2010 value, but
    ``k_fascin_pn_per_um`` is an UNSET PI-GAP slot with NO default (no single-molecule fascin stiffness exists)
    — the delegate binds the mechanism now and only a caller-supplied, sourced stiffness activates it.
    """

    n_crosslink: int
    crosslink_pairs_d: wp.array
    crosslink_active_d: wp.array
    k_fascin_pn_per_um: float
    rest_spacing_um: float = FASCIN_INTERFILAMENT_SPACING_UM
    kernel: object = _fascin_bundle_kernel
    launch: KernelLauncher = _warp_launch

    def __post_init__(self) -> None:
        if self.n_crosslink <= 0:
            raise ValueError("filopodium fascin mechanics needs a positive crosslink count")
        if not (self.k_fascin_pn_per_um > 0.0) or self.k_fascin_pn_per_um != self.k_fascin_pn_per_um:
            raise ValueError(
                "k_fascin_pn_per_um is an UNSET PI-GAP slot: supply a finite positive sourced stiffness "
                "(no single-molecule fascin value exists; do not inherit filamin 8.2e5)"
            )
        if not (self.rest_spacing_um > 0.0):
            raise ValueError("fascin rest spacing must be positive [µm]")
        if self.kernel is not _fascin_bundle_kernel:
            raise ValueError("filopodium fascin mechanics must bind the real _fascin_bundle_kernel")
        _validate_device_array(
            self.crosslink_pairs_d,
            label="fascin_mechanics.crosslink_pairs_d",
            dtype=wp.int32,
            rank=2,
            shape0=self.n_crosslink,
            shape1=2,
        )
        _validate_device_array(
            self.crosslink_active_d,
            label="fascin_mechanics.crosslink_active_d",
            dtype=wp.int32,
            shape0=self.n_crosslink,
        )

    def assert_bound_to(
        self,
        *,
        crosslink_pairs_d: wp.array,
        crosslink_active_d: wp.array,
        n_crosslink: int,
    ) -> None:
        """Reject a delegate that carries a private copy of the actor's fascin crosslink arrays."""
        if self.n_crosslink != n_crosslink:
            raise ValueError("fascin mechanics capacity disagrees with the filopodium owner")
        if _storage_key(self.crosslink_pairs_d) != _storage_key(crosslink_pairs_d):
            raise ValueError("fascin mechanics must launch against the actor-owned crosslink_pairs_d")
        if _storage_key(self.crosslink_active_d) != _storage_key(crosslink_active_d):
            raise ValueError("fascin mechanics must launch against the actor-owned crosslink_active_d")

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the real fascin crosslink kernel, scattering force into the actor-owned accumulator."""
        self.launch(
            self.kernel,
            dim=self.n_crosslink,
            inputs=[
                pos,
                self.crosslink_pairs_d,
                self.crosslink_active_d,
                self.k_fascin_pn_per_um,
                self.rest_spacing_um,
            ],
            outputs=[force],
            device=getattr(pos, "device", None),
        )


def fascin_bundle_pair_reference(
    pos_i: object,
    pos_j: object,
    *,
    k_fascin_pn_per_um: float,
    rest_spacing_um: float,
) -> tuple[list[float], list[float], float]:
    """Pure-Python mirror of one engaged ``_fascin_bundle_kernel`` slot (force on i, force on j, load).

    CPU-green ground truth for :class:`FilopodiumFascinBundleMechanics`.  ``force_on_i`` points from i toward j
    when the pair is stretched beyond ``rest_spacing``; ``force_on_j`` is its negative (Newton's third law), so
    the pair is moment-free about the connecting line and ``load = |k·(L−rest)|``.
    """
    if not (k_fascin_pn_per_um > 0.0) or k_fascin_pn_per_um != k_fascin_pn_per_um:
        raise ValueError("k_fascin_pn_per_um must be a finite positive PI-supplied stiffness")
    if not (rest_spacing_um > 0.0):
        raise ValueError("rest_spacing_um must be positive")
    d = [float(pos_j[k]) - float(pos_i[k]) for k in range(3)]
    length = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) ** 0.5
    if length <= 1.0e-12:
        return [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0
    scale = k_fascin_pn_per_um * (length - rest_spacing_um) / length
    force_on_i = [scale * d[k] for k in range(3)]
    return force_on_i, [-f for f in force_on_i], abs(k_fascin_pn_per_um * (length - rest_spacing_um))


@dataclass(frozen=True, slots=True)
class FilopodiumStateOwner:
    """CUDA owner of a local adaptive explicit bundled actin protrusion."""

    n_filament_capacity: int
    n_crosslink_capacity: int
    position_d: wp.array
    force_d: wp.array
    fiber_offset_d: wp.array
    active_filament_d: wp.array
    active_node_d: wp.array
    persistent_filament_id_d: wp.array
    refinement_parent_id_d: wp.array
    refinement_level_d: wp.array
    bundle_id_d: wp.array
    polarity_d: wp.array
    crosslink_pairs_d: wp.array
    crosslink_active_d: wp.array
    barbed_state_d: wp.array
    actor_generation_d: wp.array
    topology_epoch_d: wp.array
    mechanics: MechanicsContributor
    transaction: ProtrusionTransaction
    ledger: LedgerContributor
    component: str = FILOPODIUM_COMPONENT
    representation: str = FILOPODIUM_REPRESENTATION

    def __post_init__(self) -> None:
        if self.component != FILOPODIUM_COMPONENT:
            raise ValueError("filopodium owner cannot be relabelled")
        if self.representation != FILOPODIUM_REPRESENTATION:
            raise ValueError("filopodium must remain explicit bundled F-actin")
        if self.n_crosslink_capacity <= 0:
            raise ValueError("filopodium crosslink capacity must be positive")
        _validate_common_actor_arrays(
            component=self.component,
            n_filament_capacity=self.n_filament_capacity,
            position_d=self.position_d,
            force_d=self.force_d,
            fiber_offset_d=self.fiber_offset_d,
            active_filament_d=self.active_filament_d,
            active_node_d=self.active_node_d,
            persistent_filament_id_d=self.persistent_filament_id_d,
            refinement_parent_id_d=self.refinement_parent_id_d,
            refinement_level_d=self.refinement_level_d,
            actor_generation_d=self.actor_generation_d,
            topology_epoch_d=self.topology_epoch_d,
        )
        for label, array in (
            ("bundle_id_d", self.bundle_id_d),
            ("polarity_d", self.polarity_d),
            ("barbed_state_d", self.barbed_state_d),
        ):
            _validate_device_array(
                array,
                label=f"filopodium.{label}",
                dtype=wp.int32,
                shape0=self.n_filament_capacity,
            )
        _validate_device_array(
            self.crosslink_pairs_d,
            label="filopodium.crosslink_pairs_d",
            dtype=wp.int32,
            rank=2,
            shape0=self.n_crosslink_capacity,
            shape1=2,
        )
        _validate_device_array(
            self.crosslink_active_d,
            label="filopodium.crosslink_active_d",
            dtype=wp.int32,
            shape0=self.n_crosslink_capacity,
        )
        mutable = (
            self.position_d,
            self.fiber_offset_d,
            self.active_filament_d,
            self.active_node_d,
            self.persistent_filament_id_d,
            self.refinement_parent_id_d,
            self.refinement_level_d,
            self.bundle_id_d,
            self.polarity_d,
            self.crosslink_pairs_d,
            self.crosslink_active_d,
            self.barbed_state_d,
            self.actor_generation_d,
            self.topology_epoch_d,
        )
        if len({_storage_key(array) for array in mutable}) != len(mutable):
            raise ValueError("filopodium authoritative arrays must not alias")
        _validate_transaction_coverage(
            self.transaction,
            component=self.component,
            required_channels=FILOPODIUM_EVENT_CHANNELS,
            mutable_arrays=mutable,
        )
        if isinstance(self.mechanics, FilopodiumFascinBundleMechanics):
            self.mechanics.assert_bound_to(
                crosslink_pairs_d=self.crosslink_pairs_d,
                crosslink_active_d=self.crosslink_active_d,
                n_crosslink=self.n_crosslink_capacity,
            )

    def geometry(self) -> ProtrusionActorView:
        """Expose non-owning bundle geometry and persistent refinement identity."""
        return ProtrusionActorView(
            component=self.component,
            representation=self.representation,
            position_d=self.position_d,
            force_d=self.force_d,
            fiber_offset_d=self.fiber_offset_d,
            active_filament_d=self.active_filament_d,
            active_node_d=self.active_node_d,
            persistent_filament_id_d=self.persistent_filament_id_d,
            refinement_parent_id_d=self.refinement_parent_id_d,
            refinement_level_d=self.refinement_level_d,
            actor_generation_d=self.actor_generation_d,
            topology_epoch_d=self.topology_epoch_d,
            n_filament_capacity=self.n_filament_capacity,
        )

    def accumulate(self) -> None:
        """Launch explicit bundled-filament mechanics against actor-owned arrays."""
        self.mechanics.accumulate(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger.accumulate_ledger(ledger)


@dataclass(frozen=True, slots=True)
class ProtrusionContractProposal:
    """Local proposal for common-graph additions; does not mutate the reference graph."""

    components: tuple[ComponentContract, ...]
    connectors: tuple[ConnectorContract, ...]


def protrusion_contract_proposal() -> ProtrusionContractProposal:
    """Return the proposed component and edge names required by the protrusion facade."""
    components = (
        ComponentContract(
            LAMELLIPODIUM_COMPONENT,
            ComponentRole.ACTIVE_LOAD_PATH,
            LAMELLIPODIUM_REPRESENTATION,
            "adaptive branch-graph mechanics",
        ),
        ComponentContract(
            FILOPODIUM_COMPONENT,
            ComponentRole.ACTIVE_LOAD_PATH,
            FILOPODIUM_REPRESENTATION,
            "adaptive bundle-graph mechanics",
        ),
    )
    connectors = (
        ConnectorContract(
            LAM_MEMBRANE,
            ConnectorFamily.CONTACT,
            LAMELLIPODIUM_COMPONENT,
            MEMBRANE_COMPONENT,
            True,
            True,
            endpoint_role_a="growing barbed-end material point",
            endpoint_role_b="membrane contact quadrature",
            chemistry_card="actin_membrane_brownian_ratchet_contact",
        ),
        ConnectorContract(
            LAM_CORTEX,
            ConnectorFamily.TRANSIENT_ACTIN,
            LAMELLIPODIUM_COMPONENT,
            CORTEX_COMPONENT,
            True,
            True,
            endpoint_role_a="lamellipodial network rear seam",
            endpoint_role_b="cortex shell material point",
            chemistry_card="transient_actin_crosslink",
        ),
        ConnectorContract(
            LAM_CYTOSOL,
            ConnectorFamily.IMMERSED_TRANSFER,
            LAMELLIPODIUM_COMPONENT,
            CYTOSOL_COMPONENT,
            False,
            False,
            endpoint_role_a="active filament quadrature and barbed-end sink",
            endpoint_role_b="porous drag and G-actin transport stencil",
        ),
        ConnectorContract(
            LAM_NASCENT_FA,
            ConnectorFamily.ACTIN_ANCHOR,
            LAMELLIPODIUM_COMPONENT,
            FA_COMPONENT,
            True,
            True,
            endpoint_role_a="lamellipodial actin material point",
            endpoint_role_b="nascent FA actin-side state",
            chemistry_card="actin_talin_integrin_nascent",
        ),
        ConnectorContract(
            FILO_MEMBRANE,
            ConnectorFamily.CONTACT,
            FILOPODIUM_COMPONENT,
            MEMBRANE_COMPONENT,
            True,
            True,
            endpoint_role_a="bundled barbed-end tip",
            endpoint_role_b="membrane tip contact quadrature",
            chemistry_card="actin_membrane_brownian_ratchet_contact",
        ),
        ConnectorContract(
            FILO_CORTEX,
            ConnectorFamily.TRANSIENT_ACTIN,
            FILOPODIUM_COMPONENT,
            CORTEX_COMPONENT,
            True,
            True,
            endpoint_role_a="filopodium bundle root",
            endpoint_role_b="cortex shell material point",
            chemistry_card="formin_fascin_root_coupling",
        ),
        ConnectorContract(
            FILO_CYTOSOL,
            ConnectorFamily.IMMERSED_TRANSFER,
            FILOPODIUM_COMPONENT,
            CYTOSOL_COMPONENT,
            False,
            False,
            endpoint_role_a="bundle drag quadrature and barbed-end sink",
            endpoint_role_b="porous drag and G-actin transport stencil",
        ),
        ConnectorContract(
            FILO_NASCENT_FA,
            ConnectorFamily.ACTIN_ANCHOR,
            FILOPODIUM_COMPONENT,
            FA_COMPONENT,
            True,
            True,
            endpoint_role_a="filopodium base or shaft actin material point",
            endpoint_role_b="nascent FA actin-side state",
            chemistry_card="actin_talin_integrin_nascent",
        ),
    )
    return ProtrusionContractProposal(components=components, connectors=connectors)


def proposed_protrusion_architecture(base: CellArchitecture) -> CellArchitecture:
    """Build an idempotent validation graph containing the proposed additions.

    Existing names are retained so the helper remains valid after Lead promotion into the
    common graph.  The resulting graph is still checked against the protrusion contract;
    an existing but incompatible declaration therefore fails rather than being shadowed.
    This function always returns a new value and never mutates ``base``.
    """
    proposal = protrusion_contract_proposal()
    component_names = {component.name for component in base.components}
    connector_names = {connector.name for connector in base.connectors}
    architecture = CellArchitecture(
        components=(
            *base.components,
            *(component for component in proposal.components if component.name not in component_names),
        ),
        connectors=(
            *base.connectors,
            *(connector for connector in proposal.connectors if connector.name not in connector_names),
        ),
    )
    _validate_architecture(architecture)
    return architecture


_REQUIRED_CONNECTORS = {
    LAM_MEMBRANE: (LAMELLIPODIUM_COMPONENT, MEMBRANE_COMPONENT, ConnectorFamily.CONTACT),
    LAM_CORTEX: (LAMELLIPODIUM_COMPONENT, CORTEX_COMPONENT, ConnectorFamily.TRANSIENT_ACTIN),
    LAM_CYTOSOL: (LAMELLIPODIUM_COMPONENT, CYTOSOL_COMPONENT, ConnectorFamily.IMMERSED_TRANSFER),
    LAM_NASCENT_FA: (LAMELLIPODIUM_COMPONENT, FA_COMPONENT, ConnectorFamily.ACTIN_ANCHOR),
    FILO_MEMBRANE: (FILOPODIUM_COMPONENT, MEMBRANE_COMPONENT, ConnectorFamily.CONTACT),
    FILO_CORTEX: (FILOPODIUM_COMPONENT, CORTEX_COMPONENT, ConnectorFamily.TRANSIENT_ACTIN),
    FILO_CYTOSOL: (FILOPODIUM_COMPONENT, CYTOSOL_COMPONENT, ConnectorFamily.IMMERSED_TRANSFER),
    FILO_NASCENT_FA: (FILOPODIUM_COMPONENT, FA_COMPONENT, ConnectorFamily.ACTIN_ANCHOR),
}


def _connector_contract(architecture: CellArchitecture, name: str) -> ConnectorContract:
    for connector in architecture.connectors:
        if connector.name == name:
            return connector
    raise ValueError(f"architecture must register proposed connector {name!r}")


def _validate_architecture(architecture: CellArchitecture) -> None:
    """Validate the local proposed graph without modifying the common reference graph."""
    for component_name, representation in (
        (LAMELLIPODIUM_COMPONENT, LAMELLIPODIUM_REPRESENTATION),
        (FILOPODIUM_COMPONENT, FILOPODIUM_REPRESENTATION),
    ):
        component = architecture.component(component_name)
        if component.role is not ComponentRole.ACTIVE_LOAD_PATH:
            raise ValueError(f"{component_name} must have ComponentRole.ACTIVE_LOAD_PATH")
        if component.representation != representation:
            raise ValueError(f"{component_name} cannot use a lumped active-stress representation")
        if not component.owns_geometry or not component.dynamically_evolving:
            raise ValueError(f"{component_name} must own evolving explicit geometry")
    for name, (component_a, component_b, family) in _REQUIRED_CONNECTORS.items():
        connector = _connector_contract(architecture, name)
        if {connector.component_a, connector.component_b} != {component_a, component_b}:
            raise ValueError(f"{name} has incorrect component endpoints")
        if connector.family is not family:
            raise ValueError(f"{name} must use ConnectorFamily.{family.name}")
        if not connector.bidirectional or not connector.adjoint_transfer_required:
            raise ValueError(f"{name} must be bidirectional with adjoint transfer")
        kinetic = family is not ConnectorFamily.IMMERSED_TRANSFER
        if connector.kinetics is not kinetic or connector.commit_on_accept is not kinetic:
            raise ValueError(f"{name} has incorrect accepted-step kinetics semantics")


@dataclass(frozen=True, slots=True)
class ProtrusionStepBindings:
    """Exactly one graph-owned runtime for every required external protrusion edge."""

    connectors: tuple[ProtrusionGraphConnector, ...]

    def __post_init__(self) -> None:
        by_name: dict[str, ProtrusionGraphConnector] = {}
        for connector in self.connectors:
            if connector.name in by_name:
                raise ValueError(f"duplicate protrusion connector runtime {connector.name!r}")
            by_name[connector.name] = connector
        missing = set(_REQUIRED_CONNECTORS) - set(by_name)
        extra = set(by_name) - set(_REQUIRED_CONNECTORS)
        if missing or extra:
            raise ValueError(
                f"protrusion connector bindings mismatch; missing={sorted(missing)}, extra={sorted(extra)}"
            )
        for name, (component_a, component_b, _) in _REQUIRED_CONNECTORS.items():
            connector = by_name[name]
            if {connector.component_a, connector.component_b} != {component_a, component_b}:
                raise ValueError(f"runtime connector {name!r} has incorrect component endpoints")

    def for_actor(self, component: str) -> tuple[ProtrusionGraphConnector, ...]:
        """Return external connectors incident to one protrusion actor."""
        if component not in {LAMELLIPODIUM_COMPONENT, FILOPODIUM_COMPONENT}:
            raise KeyError(component)
        return tuple(
            connector
            for connector in self.connectors
            if component in {connector.component_a, connector.component_b}
        )


@dataclass(frozen=True, slots=True)
class ProtrusionActors:
    """Stateless facade over distinct lamellipodium and filopodium state owners."""

    architecture: CellArchitecture
    lamellipodium: LamellipodiumStateOwner
    filopodium: FilopodiumStateOwner

    def __post_init__(self) -> None:
        _validate_architecture(self.architecture)
        if self.lamellipodium.component != LAMELLIPODIUM_COMPONENT:
            raise ValueError("lamellipodium facade slot has the wrong state owner")
        if self.filopodium.component != FILOPODIUM_COMPONENT:
            raise ValueError("filopodium facade slot has the wrong state owner")
        if self.lamellipodium is self.filopodium:
            raise ValueError("lamellipodium and filopodium must be distinct state owners")
        if _storage_key(self.lamellipodium.position_d) == _storage_key(self.filopodium.position_d):
            raise ValueError("protrusion actors must not share position storage")
        if _storage_key(self.lamellipodium.force_d) == _storage_key(self.filopodium.force_d):
            raise ValueError("protrusion actors must not share force storage")

    def accumulate_mechanics(self, bindings: ProtrusionStepBindings) -> None:
        """Assemble explicit actor mechanics then registered external connector loads."""
        self.lamellipodium.accumulate()
        self.filopodium.accumulate()
        for actor in (self.lamellipodium, self.filopodium):
            view = actor.geometry()
            for connector in bindings.for_actor(actor.component):
                connector.accumulate_actor(view)

    def transaction_participants(
        self,
        bindings: ProtrusionStepBindings,
    ) -> tuple[TransactionParticipant, ...]:
        """Return both actors and all connector states exactly once."""
        return self.lamellipodium, self.filopodium, *bindings.connectors

    def snapshot_candidate(self, bindings: ProtrusionStepBindings) -> None:
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: ProtrusionStepBindings, accepted: wp.array) -> None:
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self,
        bindings: ProtrusionStepBindings,
        accepted: wp.array,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, bindings: ProtrusionStepBindings, ledger: object) -> None:
        """Forward actor and connector ledgers without owning acceptance."""
        self.lamellipodium.accumulate_ledger(ledger)
        self.filopodium.accumulate_ledger(ledger)
        for connector in bindings.connectors:
            connector.accumulate_ledger(ledger)


__all__ = [
    "CORTEX_COMPONENT",
    "CYTOSOL_COMPONENT",
    "FA_COMPONENT",
    "FASCIN_AXIAL_REPEAT_UM",
    "FASCIN_INTERFILAMENT_SPACING_UM",
    "FASCIN_K_CROSSLINK_PN_PER_UM",
    "FILOPODIUM_COMPONENT",
    "FILOPODIUM_EVENT_CHANNELS",
    "FILOPODIUM_REPRESENTATION",
    "FILO_CORTEX",
    "FILO_CYTOSOL",
    "FILO_MEMBRANE",
    "FILO_NASCENT_FA",
    "FilopodiumFascinBundleMechanics",
    "FilopodiumStateOwner",
    "KernelLauncher",
    "LAM_BRANCH_K_THETA_PN_UM",
    "LAM_BRANCH_THETA0_RAD",
    "LAMELLIPODIUM_COMPONENT",
    "LAMELLIPODIUM_EVENT_CHANNELS",
    "LAMELLIPODIUM_REPRESENTATION",
    "LamellipodiumBranchAngleMechanics",
    "LAM_CORTEX",
    "LAM_CYTOSOL",
    "LAM_MEMBRANE",
    "LAM_NASCENT_FA",
    "LamellipodiumStateOwner",
    "MEMBRANE_COMPONENT",
    "ProtrusionActorView",
    "ProtrusionActors",
    "ProtrusionContractProposal",
    "ProtrusionGraphConnector",
    "ProtrusionStepBindings",
    "ProtrusionTransaction",
    "RefinementClosureReport",
    "fascin_bundle_pair_reference",
    "proposed_protrusion_architecture",
    "protrusion_contract_proposal",
]
