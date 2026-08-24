"""First non-invasive runtime seam for the dynamic microtubule structural rig.

The landed :mod:`aleph.components.solid.microtubule` component validates MT bending, but it is a static bead-chain
aster.  This module does not relabel that adapter as dynamic instability.  It fixes the state and connector
ownership needed to add the missing mechanisms without changing the monolithic runtime:

* the MT component owns CUDA rod, MTOC, plus-end and accepted-topology arrays;
* cortex capture, nucleus LINC and motor states remain graph-owned step bindings;
* the transaction delegate must explicitly cover every mutable authoritative MT array; and
* the old bending primitive can be reused only with local ``node_off == 0`` indexing.

All mechanics hooks receive device arrays or existing Warp contributors.  There is no host physics path.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Protocol, runtime_checkable

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
)
from aleph.engine.runtime import (
    CytosolFieldEndpoint,
    ImmersedTransferConnector,
    LedgerContributor,
    MechanicsContributor,
    TransactionParticipant,
)

MICROTUBULE_COMPONENT = "microtubule"
CORTEX_COMPONENT = "cortex"
NUCLEUS_COMPONENT = "nucleus"
SF_COMPONENT = "sf_arc"
CYTOSOL_COMPONENT = "cytosol"
MT_CORTEX_CAPTURE = "mt_cortex_capture"
MT_NUCLEUS_LINC = "mt_nucleus_linc"
MT_SF_SPECTRAPLAKIN = "mt_sf_spectraplakin"
MT_CYTOSOL_TRANSFER = "mt_cytosol_transfer"
STATIC_ADAPTER_STATUS = "STATIC_BENDING_ADAPTER_DYNAMIC_TOPOLOGY_PENDING"


def _storage_key(array: object) -> tuple[object, ...]:
    """Return a storage identity without reading device data."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


def _validate_device_array(
    array: object,
    *,
    label: str,
    dtype: object,
    length: int | None = None,
) -> None:
    """Validate CUDA residence, scalar/vector dtype and one-dimensional shape."""
    device = getattr(array, "device", None)
    if not bool(getattr(device, "is_cuda", False)):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != dtype:
        raise TypeError(f"{label} must have dtype {dtype}")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 1:
        raise ValueError(f"{label} must be one-dimensional")
    if int(shape[0]) <= 0:
        raise ValueError(f"{label} must be non-empty")
    if length is not None and int(shape[0]) != length:
        raise ValueError(f"{label} must have length {length}, got {int(shape[0])}")


@runtime_checkable
class MicrotubuleTransaction(TransactionParticipant, Protocol):
    """Transaction delegate that proves which authoritative MT arrays it snapshots/restores."""

    component_name: str

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return authoritative arrays covered by snapshot/rollback; no device read is permitted."""


@dataclass(frozen=True, slots=True)
class MicrotubuleRigView:
    """Non-owning CUDA view exposed to graph connectors and motor ports."""

    position_d: wp.array
    force_d: wp.array
    mtoc_position_d: wp.array
    mtoc_force_d: wp.array
    fiber_offset_d: wp.array
    active_count_d: wp.array
    plus_end_node_d: wp.array
    phase_d: wp.array
    topology_epoch_d: wp.array
    n_mt: int


@dataclass(frozen=True, slots=True)
class ExternalMechanicalEndpoint:
    """Non-owning graph endpoint view for cortex or nucleus mechanics."""

    component: str
    position_d: wp.array
    force_d: wp.array

    def __post_init__(self) -> None:
        if not self.component.strip():
            raise ValueError("external endpoint component must be non-empty")
        _validate_device_array(
            self.position_d,
            label=f"{self.component}.position_d",
            dtype=wp.vec3d,
        )
        _validate_device_array(
            self.force_d,
            label=f"{self.component}.force_d",
            dtype=wp.vec3d,
            length=int(self.position_d.shape[0]),
        )
        if _storage_key(self.position_d) == _storage_key(self.force_d):
            raise ValueError(f"{self.component} position and force arrays must not alias")
        if str(self.position_d.device) != str(self.force_d.device):
            raise ValueError(f"{self.component} position and force arrays must share a CUDA device")


@runtime_checkable
class MicrotubuleGraphConnector(Protocol):
    """Graph-owned joint connecting an MT material/MTOC endpoint to another component."""

    name: str
    component_a: str
    component_b: str

    def accumulate_rig(
        self,
        rig: MicrotubuleRigView,
        endpoint: ExternalMechanicalEndpoint,
    ) -> None:
        """Scatter equal-and-opposite connector loads through public device views."""

    def snapshot_candidate(self) -> None:
        """Snapshot graph-owned connector state."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore under the scheduler-owned rejected predicate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit capture/LINC kinetics only under the accepted predicate."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Add connector force/work/topology terms without deciding acceptance."""


@runtime_checkable
class MicrotubuleMotorPort(Protocol):
    """Graph-owned motor attached at a live polar MT material coordinate."""

    name: str
    component_a: str
    component_b: str

    def accumulate_motor(self, rig: MicrotubuleRigView) -> None:
        """Apply motor/cargo reactions through the live rod interpolation view."""

    def snapshot_candidate(self) -> None:
        """Snapshot motor binding/abscissa/epoch state."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore motor state under rejection."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit stepping/detachment only under acceptance."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Add motor work/binding terms without deciding acceptance."""


@dataclass(frozen=True, slots=True)
class MicrotubuleStepBindings:
    """External graph endpoints and joints supplied for one coupled candidate step."""

    cortex: ExternalMechanicalEndpoint
    nucleus: ExternalMechanicalEndpoint
    sf_arc: ExternalMechanicalEndpoint
    cytosol: CytosolFieldEndpoint
    cortex_capture: MicrotubuleGraphConnector
    nucleus_linc: MicrotubuleGraphConnector
    sf_spectraplakin: MicrotubuleGraphConnector
    cytosol_transfer: ImmersedTransferConnector
    motor_ports: tuple[MicrotubuleMotorPort, ...] = ()

    def __post_init__(self) -> None:
        if self.cortex.component != CORTEX_COMPONENT:
            raise ValueError("cortex endpoint must be the registered cortex component")
        if self.nucleus.component != NUCLEUS_COMPONENT:
            raise ValueError("nucleus endpoint must be the registered nucleus component")
        if self.sf_arc.component != SF_COMPONENT:
            raise ValueError("SF endpoint must be the registered sf_arc component")
        _validate_runtime_connector(
            self.cortex_capture,
            expected_name=MT_CORTEX_CAPTURE,
            endpoints={MICROTUBULE_COMPONENT, CORTEX_COMPONENT},
        )
        _validate_runtime_connector(
            self.nucleus_linc,
            expected_name=MT_NUCLEUS_LINC,
            endpoints={MICROTUBULE_COMPONENT, NUCLEUS_COMPONENT},
        )
        _validate_runtime_connector(
            self.sf_spectraplakin,
            expected_name=MT_SF_SPECTRAPLAKIN,
            endpoints={MICROTUBULE_COMPONENT, SF_COMPONENT},
        )
        _validate_runtime_connector(
            self.cytosol_transfer,
            expected_name=MT_CYTOSOL_TRANSFER,
            endpoints={MICROTUBULE_COMPONENT, CYTOSOL_COMPONENT},
        )
        for motor in self.motor_ports:
            if MICROTUBULE_COMPONENT not in {motor.component_a, motor.component_b}:
                raise ValueError(f"motor port {motor.name!r} must have a microtubule endpoint")


def _validate_runtime_connector(
    connector: object,
    *,
    expected_name: str,
    endpoints: set[str],
) -> None:
    """Validate a graph-owned runtime view without adopting its state."""
    if getattr(connector, "name", None) != expected_name:
        raise ValueError(f"expected graph connector {expected_name!r}")
    if {
        getattr(connector, "component_a", None),
        getattr(connector, "component_b", None),
    } != endpoints:
        raise ValueError(f"{expected_name} has incorrect component endpoints")


@dataclass(frozen=True, slots=True)
class MicrotubuleRigStateOwner:
    """Authoritative device state for one dynamic polar MT rig.

    The first slice owns the complete topology seam but may delegate mechanics to the static bending adapter.
    ``status`` prevents callers from reporting that adapter as completed dynamic instability.
    """

    n_mt: int
    position_d: wp.array
    force_d: wp.array
    mtoc_position_d: wp.array
    mtoc_force_d: wp.array
    fiber_offset_d: wp.array
    active_count_d: wp.array
    plus_end_node_d: wp.array
    phase_d: wp.array
    topology_epoch_d: wp.array
    mechanics: MechanicsContributor
    transaction: MicrotubuleTransaction
    ledger: LedgerContributor
    status: str = STATIC_ADAPTER_STATUS

    def __post_init__(self) -> None:
        if self.n_mt <= 0:
            raise ValueError("n_mt must be positive")
        if self.status != STATIC_ADAPTER_STATUS:
            raise ValueError("the first MT slice cannot be relabelled as completed dynamic topology")
        _validate_device_array(self.position_d, label="microtubule.position_d", dtype=wp.vec3d)
        n_nodes = int(self.position_d.shape[0])
        _validate_device_array(
            self.force_d,
            label="microtubule.force_d",
            dtype=wp.vec3d,
            length=n_nodes,
        )
        _validate_device_array(
            self.mtoc_position_d,
            label="microtubule.mtoc_position_d",
            dtype=wp.vec3d,
            length=1,
        )
        _validate_device_array(
            self.mtoc_force_d,
            label="microtubule.mtoc_force_d",
            dtype=wp.vec3d,
            length=1,
        )
        _validate_device_array(
            self.fiber_offset_d,
            label="microtubule.fiber_offset_d",
            dtype=wp.int32,
            length=self.n_mt + 1,
        )
        for label, array in (
            ("active_count_d", self.active_count_d),
            ("plus_end_node_d", self.plus_end_node_d),
            ("phase_d", self.phase_d),
        ):
            _validate_device_array(
                array,
                label=f"microtubule.{label}",
                dtype=wp.int32,
                length=self.n_mt,
            )
        _validate_device_array(
            self.topology_epoch_d,
            label="microtubule.topology_epoch_d",
            dtype=wp.int32,
            length=1,
        )
        device_arrays = (
            self.position_d,
            self.force_d,
            self.mtoc_position_d,
            self.mtoc_force_d,
            self.fiber_offset_d,
            self.active_count_d,
            self.plus_end_node_d,
            self.phase_d,
            self.topology_epoch_d,
        )
        if len({str(array.device) for array in device_arrays}) != 1:
            raise ValueError("all authoritative MT arrays must share one CUDA device")
        mutable = (
            self.position_d,
            self.mtoc_position_d,
            self.active_count_d,
            self.plus_end_node_d,
            self.phase_d,
            self.topology_epoch_d,
        )
        if len({_storage_key(array) for array in mutable}) != len(mutable):
            raise ValueError("authoritative MT geometry/topology arrays must not alias")
        if _storage_key(self.position_d) == _storage_key(self.force_d):
            raise ValueError("microtubule position and force arrays must not alias")
        if _storage_key(self.mtoc_position_d) == _storage_key(self.mtoc_force_d):
            raise ValueError("MTOC position and force arrays must not alias")
        n_expected = getattr(self.mechanics, "n_nodes", None)
        if n_expected is not None and int(n_expected) != n_nodes:
            raise ValueError(f"MT mechanics expects {int(n_expected)} nodes but state owns {n_nodes}")
        if self.transaction.component_name != MICROTUBULE_COMPONENT:
            raise ValueError("MT transaction delegate must be owned by the microtubule component")
        covered = {_storage_key(array) for array in self.transaction.owned_arrays()}
        missing = {_storage_key(array) for array in mutable} - covered
        if missing:
            raise ValueError("MT transaction must cover rod, MTOC, endpoint, phase and topology epoch arrays")

    def geometry(self) -> MicrotubuleRigView:
        """Return a non-owning device view for registered connectors."""
        return MicrotubuleRigView(
            position_d=self.position_d,
            force_d=self.force_d,
            mtoc_position_d=self.mtoc_position_d,
            mtoc_force_d=self.mtoc_force_d,
            fiber_offset_d=self.fiber_offset_d,
            active_count_d=self.active_count_d,
            plus_end_node_d=self.plus_end_node_d,
            phase_d=self.phase_d,
            topology_epoch_d=self.topology_epoch_d,
            n_mt=self.n_mt,
        )

    def accumulate(self) -> None:
        """Launch the currently selected Warp rod mechanics against component-owned arrays."""
        self.mechanics.accumulate(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        """Snapshot geometry and accepted topology through the verified delegate."""
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Restore geometry/topology under the scheduler-owned predicate."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit proposed dynamic topology only under final acceptance."""
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Add MT mechanics/topology ledger terms without deciding acceptance."""
        self.ledger.accumulate_ledger(ledger)


def _find_connector(architecture: CellArchitecture, name: str) -> ConnectorContract:
    for connector in architecture.connectors:
        if connector.name == name:
            return connector
    raise ValueError(f"architecture must register {name!r}")


def _validate_build_connector(
    architecture: CellArchitecture,
    *,
    name: str,
    family: ConnectorFamily,
    endpoints: set[str],
    kinetic: bool = True,
) -> None:
    connector = _find_connector(architecture, name)
    if connector.family is not family:
        raise ValueError(f"{name} must use ConnectorFamily.{family.name}")
    if {connector.component_a, connector.component_b} != endpoints:
        raise ValueError(f"{name} has incorrect component endpoints")
    if not connector.bidirectional or not connector.adjoint_transfer_required:
        raise ValueError(f"{name} must be bidirectional with adjoint transfer")
    if connector.kinetics is not kinetic or connector.commit_on_accept is not kinetic:
        qualifier = "kinetic and accepted-step committed" if kinetic else "non-kinetic"
        raise ValueError(f"{name} must be {qualifier}")


def _validate_architecture(architecture: CellArchitecture) -> None:
    component = architecture.component(MICROTUBULE_COMPONENT)
    if component.role is not ComponentRole.STRUCTURAL_RIG:
        raise ValueError("microtubule must have ComponentRole.STRUCTURAL_RIG")
    if not component.owns_geometry or not component.dynamically_evolving:
        raise ValueError("microtubule must own geometry and be dynamically evolving")
    _validate_build_connector(
        architecture,
        name=MT_CORTEX_CAPTURE,
        family=ConnectorFamily.MOTOR,
        endpoints={MICROTUBULE_COMPONENT, CORTEX_COMPONENT},
    )
    _validate_build_connector(
        architecture,
        name=MT_NUCLEUS_LINC,
        family=ConnectorFamily.LINC,
        endpoints={MICROTUBULE_COMPONENT, NUCLEUS_COMPONENT},
    )
    _validate_build_connector(
        architecture,
        name=MT_SF_SPECTRAPLAKIN,
        family=ConnectorFamily.SPECTRAPLAKIN,
        endpoints={MICROTUBULE_COMPONENT, SF_COMPONENT},
    )
    _validate_build_connector(
        architecture,
        name=MT_CYTOSOL_TRANSFER,
        family=ConnectorFamily.IMMERSED_TRANSFER,
        endpoints={MICROTUBULE_COMPONENT, CYTOSOL_COMPONENT},
        kinetic=False,
    )


@dataclass(frozen=True, slots=True)
class MicrotubuleRig:
    """Registered MT component facade; graph connector state is always supplied externally."""

    architecture: CellArchitecture
    state: MicrotubuleRigStateOwner

    def __post_init__(self) -> None:
        _validate_architecture(self.architecture)

    def accumulate_mechanics(self, bindings: MicrotubuleStepBindings) -> None:
        """Launch rod bending and every supplied graph-owned load path."""
        self.state.accumulate()
        view = self.state.geometry()
        bindings.cortex_capture.accumulate_rig(view, bindings.cortex)
        bindings.nucleus_linc.accumulate_rig(view, bindings.nucleus)
        bindings.sf_spectraplakin.accumulate_rig(view, bindings.sf_arc)
        bindings.cytosol_transfer.accumulate_transfer(view, bindings.cytosol)
        for motor in bindings.motor_ports:
            motor.accumulate_motor(view)

    def transaction_participants(
        self,
        bindings: MicrotubuleStepBindings,
    ) -> tuple[TransactionParticipant, ...]:
        """Expose component plus graph-owned topology-changing joints to the global scheduler."""
        return (
            self.state,
            bindings.cortex_capture,
            bindings.nucleus_linc,
            bindings.sf_spectraplakin,
            bindings.cytosol_transfer,
            *bindings.motor_ports,
        )

    def snapshot_candidate(self, bindings: MicrotubuleStepBindings) -> None:
        """Snapshot all participants without adopting connector state."""
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: MicrotubuleStepBindings, accepted: wp.array) -> None:
        """Forward one rejected predicate to state and graph connectors."""
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self,
        bindings: MicrotubuleStepBindings,
        accepted: wp.array,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        """Commit dynamic topology and connector kinetics exactly once on acceptance."""
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, bindings: MicrotubuleStepBindings, ledger: object) -> None:
        """Collect component and connector ledgers without owning the verdict."""
        self.state.accumulate_ledger(ledger)
        bindings.cortex_capture.accumulate_ledger(ledger)
        bindings.nucleus_linc.accumulate_ledger(ledger)
        bindings.sf_spectraplakin.accumulate_ledger(ledger)
        bindings.cytosol_transfer.accumulate_ledger(ledger)
        for motor in bindings.motor_ports:
            motor.accumulate_ledger(ledger)


class _LegacyMicrotubuleCompartment(Protocol):
    node_off: int
    n_nodes: int
    n_mt: int

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the existing Cytosim MT bending kernel."""


@dataclass(frozen=True, slots=True)
class LegacyMicrotubuleBendingAdapter:
    """Reuse the landed bending primitive without claiming dynamic topology or graph capture."""

    compartment: _LegacyMicrotubuleCompartment

    def __post_init__(self) -> None:
        if int(self.compartment.node_off) != 0:
            raise ValueError("legacy MT adapter requires local indexing with node_off == 0")
        if int(self.compartment.n_nodes) <= 0 or int(self.compartment.n_mt) <= 0:
            raise ValueError("legacy MT adapter requires a non-empty aster")

    @property
    def n_nodes(self) -> int:
        return int(self.compartment.n_nodes)

    @property
    def status(self) -> str:
        return STATIC_ADAPTER_STATUS

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Delegate bending to the existing Warp primitive."""
        self.compartment.accumulate(pos, force)


# --------------------------------------------------------------------------------------------------
# KERNEL_BOUND increment: MTOC-as-live-body anchor reaction + accepted-step dynamic-instability topology
# --------------------------------------------------------------------------------------------------
#
# The seam above owns arrays but launches no dynamic MT physics.  This section binds real Warp kernels for
# the two mechanisms the plan (MICROTUBULE_PLAN.md §§2,4) flags as genuinely absent from the static
# bead-chain adapter:
#
#   1. the MTOC minus-end anchor reaction — a Newton-3rd equal-and-opposite spring between each arm's
#      innermost bead and the single MTOC node, so the MTOC is a live mechanical body rather than a host
#      reference position; and
#   2. the accepted-step dynamic-instability transaction — continuous plus-end length accumulation plus
#      discrete terminal-span activation/deactivation, proposed into candidate arrays and committed ONLY
#      under the scheduler-owned acceptance predicate (rejection restores committed topology bit-exactly).
#
# Bending physics is NOT re-derived here: it is reused through an injected ``MechanicsContributor`` delegate
# (the landed ``cytosim_bending_kernel`` via ``LegacyMicrotubuleBendingAdapter``).  Every physical constant
# (anchor stiffness ``k_hub``, rest length, per-MT ``v_grow``/``v_shrink``) is INJECTED by the caller; the
# only literal length is the geometric segment spacing ``seg_um``.  Phase labels are an ``IntEnum``, not a
# physics constant.
#
# STATUS, corrected 2026-07-29 — the half of this note that said the RNG transitions were deferred was stale.
# Catastrophe/rescue is IMPLEMENTED and stochastic: :func:`_mt_di_stochastic_proposal_kernel` draws one
# uniform per MT and switches phase with probability ``1-exp(-f·dt)``, commit/rollback go through
# :func:`_mt_topology_commit_kernel` / :func:`_mt_topology_rollback_kernel`, and
# :meth:`MicrotubuleRodRuntime.propose_dynamic_instability_stochastic` binds it.  What IS still deferred is
# **load dependence** — the kernel takes ``f_catastrophe_per_s`` / ``f_rescue_per_s`` as SCALARS and no force
# array, so the rates cannot yet respond to tip load.  The rate VALUES stay an unset card by design: they are
# priors to be inferred, not constants to be signed (``rusan_2001_llcpk_proxy`` supplies the candidate).


class MicrotubulePhase(IntEnum):
    """Dynamic-instability phase of one MT plus end (Mitchison-Kirschner 1984 state labels, not rates)."""

    GROWING = 0
    SHRINKING = 1
    PAUSED = 2


_PHASE_GROWING = wp.constant(int(MicrotubulePhase.GROWING))
_PHASE_SHRINKING = wp.constant(int(MicrotubulePhase.SHRINKING))
_MIN_ACTIVE_NODES = wp.constant(2)  # a rod needs >=2 nodes for one bending triple; not a tuning constant


@dataclass(frozen=True, slots=True)
class DynamicInstabilityRateCard:
    """Catastrophe/rescue Poisson rates for the plus-end phase switch — an UNSET PI-GAP slot.

    ``f_catastrophe_per_s`` (GROWING→SHRINKING) and ``f_rescue_per_s`` (SHRINKING→GROWING) have NO default:
    the MCF7-specific dynamic-instability rates are a PI GAP (``PI_GAP_LITERATURE_SOURCING_2026-07-23`` §2), so
    a caller must supply a sourced card.  ``provenance`` records where the numbers came from and ``ratified``
    is set True only by PI — the LLCPK-1α proxy factory ships it ``False`` so an unratified proxy cannot
    silently drive a native run.  Pause entry/exit rates are a SEPARATE NOT-FOUND slot and are deliberately
    absent here: this card is the observed 2-state catastrophe/rescue switch only (the ``PAUSED`` phase is
    inert under it), and whether to model a 3rd stochastic state is an open PI decision.
    """

    f_catastrophe_per_s: float
    f_rescue_per_s: float
    provenance: str
    ratified: bool = False

    def __post_init__(self) -> None:
        for label, value in (
            ("f_catastrophe_per_s", self.f_catastrophe_per_s),
            ("f_rescue_per_s", self.f_rescue_per_s),
        ):
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{label} must be finite and nonnegative [s^-1]")
        if not self.provenance.strip():
            raise ValueError("provenance must record the rate source (or an explicit PI-GAP marker)")

    @classmethod
    def rusan_2001_llcpk_proxy(cls) -> DynamicInstabilityRateCard:
        """LLCPK-1α interphase plus-end proxy (Rusan 2001, MBoC 12:971) — CANDIDATE, PI ratification pending.

        f_cat = 0.026 s⁻¹, f_rescue = 0.175 s⁻¹ (live-cell GFP-α-tubulin, interphase, plus-end).  MCF7-specific
        rates are NOT FOUND; this in-vivo mammalian-epithelial set is the recommended proxy but is NOT
        auto-adopted — ``ratified`` stays ``False`` until PI signs off (mirrors the HeLa cortical-tension proxy
        convention).  Growth/shrink SPEEDS remain per-MT injected arrays, not part of this rate card.
        """
        return cls(
            f_catastrophe_per_s=0.026,
            f_rescue_per_s=0.175,
            provenance="Rusan 2001 MBoC 12:971 LLCPK-1a interphase plus-end PROXY (PI ratification pending)",
            ratified=False,
        )


@wp.kernel
def _mtoc_anchor_force_kernel(
    position: wp.array(dtype=wp.vec3d),
    rod_force: wp.array(dtype=wp.vec3d),
    mtoc_position: wp.array(dtype=wp.vec3d),
    mtoc_force: wp.array(dtype=wp.vec3d),
    fiber_offset: wp.array(dtype=wp.int32),
    k_hub: wp.float64,
    rest: wp.float64,
    anchor_load: wp.array(dtype=wp.float64),
    anchor_energy: wp.array(dtype=wp.float64),
) -> None:
    """One thread per MT: Hookean minus-end anchor between the arm base and the single MTOC node.

    Force on the arm base points from base toward the MTOC in extension; the equal-and-opposite reaction is
    scattered onto the shared MTOC node (atomic accumulate over all arms).  Gather and scatter share the same
    endpoints, so nodal virtual work equals -dU.
    """
    m = wp.tid()
    base = fiber_offset[m]
    delta = mtoc_position[0] - position[base]
    length = wp.length(delta)
    pair_force = wp.vec3d(0.0, 0.0, 0.0)
    extension = length - rest
    if length > wp.float64(1.0e-15):
        pair_force = k_hub * extension * delta / length
    wp.atomic_add(rod_force, base, pair_force)
    wp.atomic_add(mtoc_force, 0, -pair_force)
    anchor_load[m] = wp.abs(k_hub * extension)
    anchor_energy[m] = wp.float64(0.5) * k_hub * extension * extension


@wp.kernel
def _mt_di_proposal_kernel(
    phase: wp.array(dtype=wp.int32),
    active_count: wp.array(dtype=wp.int32),
    plus_end_node: wp.array(dtype=wp.int32),
    partial_length: wp.array(dtype=wp.float64),
    fiber_offset: wp.array(dtype=wp.int32),
    v_grow: wp.array(dtype=wp.float64),
    v_shrink: wp.array(dtype=wp.float64),
    seg_um: wp.float64,
    dt_phys: wp.float64,
    candidate_phase: wp.array(dtype=wp.int32),
    candidate_active_count: wp.array(dtype=wp.int32),
    candidate_plus_end: wp.array(dtype=wp.int32),
    candidate_partial_length: wp.array(dtype=wp.float64),
) -> None:
    """Propose one plus-end length increment into candidate arrays; never mutate committed topology.

    Length accumulates continuously from the injected per-MT velocities; when the partial terminal span
    crosses +-``seg_um`` exactly one dormant node is activated or the last node deactivated, so topology is
    not quantised to a full segment per physical step.  Bounds: ``active_count`` stays in
    ``[2, capacity]``.
    """
    m = wp.tid()
    ph = phase[m]
    v = wp.float64(0.0)
    if ph == _PHASE_GROWING:
        v = v_grow[m]
    elif ph == _PHASE_SHRINKING:
        v = -v_shrink[m]
    new_partial = partial_length[m] + v * dt_phys
    ac = active_count[m]
    pe = plus_end_node[m]
    cap_end = fiber_offset[m + 1]  # exclusive dormant-capacity bound for this arm
    if new_partial >= seg_um:
        if pe + 1 < cap_end:
            ac = ac + 1
            pe = pe + 1
            new_partial = new_partial - seg_um
        else:
            new_partial = seg_um  # clamp at allocated capacity; no host reallocation
    elif new_partial <= -seg_um:
        if ac > _MIN_ACTIVE_NODES:
            ac = ac - 1
            pe = pe - 1
            new_partial = new_partial + seg_um
        else:
            new_partial = -seg_um  # clamp at the minimum viable rod
    candidate_phase[m] = ph
    candidate_active_count[m] = ac
    candidate_plus_end[m] = pe
    candidate_partial_length[m] = new_partial


@wp.kernel
def _mt_di_stochastic_proposal_kernel(
    phase: wp.array(dtype=wp.int32),
    active_count: wp.array(dtype=wp.int32),
    plus_end_node: wp.array(dtype=wp.int32),
    partial_length: wp.array(dtype=wp.float64),
    fiber_offset: wp.array(dtype=wp.int32),
    v_grow: wp.array(dtype=wp.float64),
    v_shrink: wp.array(dtype=wp.float64),
    f_catastrophe_per_s: wp.float64,
    f_rescue_per_s: wp.float64,
    seg_um: wp.float64,
    dt_phys: wp.float64,
    rng_seed: wp.int32,
    candidate_phase: wp.array(dtype=wp.int32),
    candidate_active_count: wp.array(dtype=wp.int32),
    candidate_plus_end: wp.array(dtype=wp.int32),
    candidate_partial_length: wp.array(dtype=wp.float64),
) -> None:
    """Stochastic dynamic instability: Poisson catastrophe/rescue phase switch + plus-end length increment.

    One thread per MT draws a single uniform on ``wp.rand_init(rng_seed, m)``.  A GROWING plus end switches to
    SHRINKING with probability ``1-exp(-f_cat·dt)`` (catastrophe); a SHRINKING plus end switches to GROWING
    with ``1-exp(-f_rescue·dt)`` (rescue).  PAUSED is inert (no PI-sourced pause rate).  The length increment
    then uses the POST-transition phase, so a catastrophe shrinks within the same step.  Proposes into
    candidate arrays only; committed topology is untouched (commit/rollback are the existing kernels).  With a
    zero-rate card this reduces exactly to :func:`_mt_di_proposal_kernel` (the deterministic pre-PI path).
    """
    m = wp.tid()
    ph = phase[m]
    rstate = wp.rand_init(rng_seed, m)
    u = wp.float64(wp.randf(rstate))
    if ph == _PHASE_GROWING:
        p_cat = wp.float64(1.0) - wp.exp(-f_catastrophe_per_s * dt_phys)
        if u < p_cat:
            ph = _PHASE_SHRINKING
    elif ph == _PHASE_SHRINKING:
        p_res = wp.float64(1.0) - wp.exp(-f_rescue_per_s * dt_phys)
        if u < p_res:
            ph = _PHASE_GROWING
    v = wp.float64(0.0)
    if ph == _PHASE_GROWING:
        v = v_grow[m]
    elif ph == _PHASE_SHRINKING:
        v = -v_shrink[m]
    new_partial = partial_length[m] + v * dt_phys
    ac = active_count[m]
    pe = plus_end_node[m]
    cap_end = fiber_offset[m + 1]
    if new_partial >= seg_um:
        if pe + 1 < cap_end:
            ac = ac + 1
            pe = pe + 1
            new_partial = new_partial - seg_um
        else:
            new_partial = seg_um
    elif new_partial <= -seg_um:
        if ac > _MIN_ACTIVE_NODES:
            ac = ac - 1
            pe = pe - 1
            new_partial = new_partial + seg_um
        else:
            new_partial = -seg_um
    candidate_phase[m] = ph
    candidate_active_count[m] = ac
    candidate_plus_end[m] = pe
    candidate_partial_length[m] = new_partial


@wp.kernel
def _mt_topology_commit_kernel(
    accepted: wp.array(dtype=wp.int32),
    candidate_phase: wp.array(dtype=wp.int32),
    candidate_active_count: wp.array(dtype=wp.int32),
    candidate_plus_end: wp.array(dtype=wp.int32),
    candidate_partial_length: wp.array(dtype=wp.float64),
    phase: wp.array(dtype=wp.int32),
    active_count: wp.array(dtype=wp.int32),
    plus_end_node: wp.array(dtype=wp.int32),
    partial_length: wp.array(dtype=wp.float64),
    topology_epoch: wp.array(dtype=wp.int32),
) -> None:
    """Commit proposed topology only under the accepted predicate; on rejection reset candidate scratch."""
    m = wp.tid()
    if accepted[0] != 0:
        phase[m] = candidate_phase[m]
        active_count[m] = candidate_active_count[m]
        plus_end_node[m] = candidate_plus_end[m]
        partial_length[m] = candidate_partial_length[m]
        if m == 0:
            topology_epoch[0] = topology_epoch[0] + wp.int32(1)
    else:
        candidate_phase[m] = phase[m]
        candidate_active_count[m] = active_count[m]
        candidate_plus_end[m] = plus_end_node[m]
        candidate_partial_length[m] = partial_length[m]


@wp.kernel
def _mt_topology_rollback_kernel(
    accepted: wp.array(dtype=wp.int32),
    phase: wp.array(dtype=wp.int32),
    active_count: wp.array(dtype=wp.int32),
    plus_end_node: wp.array(dtype=wp.int32),
    partial_length: wp.array(dtype=wp.float64),
    candidate_phase: wp.array(dtype=wp.int32),
    candidate_active_count: wp.array(dtype=wp.int32),
    candidate_plus_end: wp.array(dtype=wp.int32),
    candidate_partial_length: wp.array(dtype=wp.float64),
) -> None:
    """Reset candidate scratch to committed topology after rejection; committed state is never mutated."""
    m = wp.tid()
    if accepted[0] == 0:
        candidate_phase[m] = phase[m]
        candidate_active_count[m] = active_count[m]
        candidate_plus_end[m] = plus_end_node[m]
        candidate_partial_length[m] = partial_length[m]


@dataclass(frozen=True, slots=True)
class MtocAnchorReference:
    """Host virtual-work oracle for one MTOC minus-end anchor spring."""

    force_on_base: npt.NDArray[np.float64]
    force_on_mtoc: npt.NDArray[np.float64]
    load_pn: float
    energy_pn_um: float


def mtoc_anchor_reference(
    base_position: npt.ArrayLike,
    mtoc_position: npt.ArrayLike,
    k_hub_pn_per_um: float,
    rest_um: float,
) -> MtocAnchorReference:
    """Scatter one Hookean MTOC anchor and return the equal-and-opposite base/MTOC reactions.

    ``force_on_base`` points from the arm base toward the MTOC in extension; the MTOC reaction is its
    negative, so the pair is moment-free about the connecting line and ``load = |k*extension|``.
    """
    base = np.asarray(base_position, dtype=np.float64)
    hub = np.asarray(mtoc_position, dtype=np.float64)
    if base.shape != (3,) or hub.shape != (3,):
        raise ValueError("anchor endpoints must be 3-vectors")
    if not np.isfinite(k_hub_pn_per_um) or k_hub_pn_per_um <= 0.0:
        raise ValueError("anchor stiffness must be finite and positive")
    if not np.isfinite(rest_um) or rest_um < 0.0:
        raise ValueError("anchor rest length must be finite and nonnegative")
    delta = hub - base
    length = float(np.linalg.norm(delta))
    extension = length - rest_um
    force = np.zeros(3, dtype=np.float64)
    if length > 1.0e-15:
        force = k_hub_pn_per_um * extension * delta / length
    return MtocAnchorReference(
        force_on_base=force,
        force_on_mtoc=-force,
        load_pn=float(abs(k_hub_pn_per_um * extension)),
        energy_pn_um=float(0.5 * k_hub_pn_per_um * extension**2),
    )


def dynamic_instability_reference(
    *,
    phase: int,
    active_count: int,
    plus_end_node: int,
    partial_length_um: float,
    fiber_offset_start: int,
    fiber_offset_end: int,
    v_grow_um_per_s: float,
    v_shrink_um_per_s: float,
    seg_um: float,
    dt_phys: float,
) -> tuple[int, int, int, float]:
    """Pure-Python mirror of :func:`_mt_di_proposal_kernel` (candidate phase, active, plus_end, partial).

    This is the CPU-green ground truth for the accepted-topology proposal.  ``fiber_offset_start`` is the
    arm base node; ``fiber_offset_end`` is the exclusive dormant-capacity bound.
    """
    if seg_um <= 0.0 or dt_phys <= 0.0:
        raise ValueError("seg_um and dt_phys must be positive")
    if not fiber_offset_start <= plus_end_node < fiber_offset_end:
        raise ValueError("plus-end node must lie within the arm's allocated capacity")
    v = 0.0
    if phase == int(MicrotubulePhase.GROWING):
        v = v_grow_um_per_s
    elif phase == int(MicrotubulePhase.SHRINKING):
        v = -v_shrink_um_per_s
    new_partial = partial_length_um + v * dt_phys
    ac = active_count
    pe = plus_end_node
    if new_partial >= seg_um:
        if pe + 1 < fiber_offset_end:
            ac += 1
            pe += 1
            new_partial -= seg_um
        else:
            new_partial = seg_um
    elif new_partial <= -seg_um:
        if ac > int(_MIN_ACTIVE_NODES):
            ac -= 1
            pe -= 1
            new_partial += seg_um
        else:
            new_partial = -seg_um
    return int(phase), int(ac), int(pe), float(new_partial)


def dynamic_instability_stochastic_reference(
    *,
    phase: int,
    active_count: int,
    plus_end_node: int,
    partial_length_um: float,
    fiber_offset_start: int,
    fiber_offset_end: int,
    v_grow_um_per_s: float,
    v_shrink_um_per_s: float,
    f_catastrophe_per_s: float,
    f_rescue_per_s: float,
    uniform_draw: float,
    seg_um: float,
    dt_phys: float,
) -> tuple[int, int, int, float]:
    """Pure-Python mirror of :func:`_mt_di_stochastic_proposal_kernel` (candidate phase, active, plus_end, partial).

    Deterministic in ``uniform_draw`` (the same U(0,1) the kernel samples), so it is the CPU-green ground truth
    without any CUDA RNG.  The catastrophe/rescue switch runs first, then the length step delegates to
    :func:`dynamic_instability_reference` under the transitioned phase.
    """
    if not (0.0 <= uniform_draw <= 1.0):
        raise ValueError("uniform_draw must lie in [0, 1]")
    for label, value in (
        ("f_catastrophe_per_s", f_catastrophe_per_s),
        ("f_rescue_per_s", f_rescue_per_s),
    ):
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{label} must be finite and nonnegative")
    ph = int(phase)
    if ph == int(MicrotubulePhase.GROWING):
        if uniform_draw < 1.0 - np.exp(-f_catastrophe_per_s * dt_phys):
            ph = int(MicrotubulePhase.SHRINKING)
    elif ph == int(MicrotubulePhase.SHRINKING):
        if uniform_draw < 1.0 - np.exp(-f_rescue_per_s * dt_phys):
            ph = int(MicrotubulePhase.GROWING)
    return dynamic_instability_reference(
        phase=ph,
        active_count=active_count,
        plus_end_node=plus_end_node,
        partial_length_um=partial_length_um,
        fiber_offset_start=fiber_offset_start,
        fiber_offset_end=fiber_offset_end,
        v_grow_um_per_s=v_grow_um_per_s,
        v_shrink_um_per_s=v_shrink_um_per_s,
        seg_um=seg_um,
        dt_phys=dt_phys,
    )


@dataclass(frozen=True, slots=True)
class MicrotubuleRodBuildSpec:
    """Host-validated build inputs for a CUDA MT rod runtime (validated with no device access).

    Arrays are micrometre/piconewton/second.  ``fiber_offset`` has length ``n_mt + 1`` and partitions the
    node array into per-arm capacity blocks; the innermost node of each arm is the MTOC-anchored base.
    Every physical constant is supplied by the caller — none is defaulted to a biological value here.
    """

    position: npt.NDArray[np.float64]
    fiber_offset: npt.NDArray[np.int32]
    mtoc_position: npt.NDArray[np.float64]
    active_count: npt.NDArray[np.int32]
    phase: npt.NDArray[np.int32]
    v_grow_um_per_s: npt.NDArray[np.float64]
    v_shrink_um_per_s: npt.NDArray[np.float64]
    k_hub_pn_per_um: float
    anchor_rest_um: float
    seg_um: float

    @property
    def n_mt(self) -> int:
        return int(self.fiber_offset.shape[0]) - 1

    @property
    def n_nodes(self) -> int:
        return int(self.position.shape[0])

    def plus_end_node(self) -> npt.NDArray[np.int32]:
        """Derive each plus-end node from base offset + active span (never a free global index)."""
        base = self.fiber_offset[:-1].astype(np.int64)
        return (base + self.active_count.astype(np.int64) - 1).astype(np.int32)

    def __post_init__(self) -> None:
        pos = np.asarray(self.position, dtype=np.float64)
        off = np.asarray(self.fiber_offset, dtype=np.int64)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError("position must be (n_nodes, 3)")
        if off.ndim != 1 or off.shape[0] < 2:
            raise ValueError("fiber_offset must be (n_mt + 1,)")
        if off[0] != 0 or off[-1] != pos.shape[0]:
            raise ValueError("fiber_offset must start at 0 and end at n_nodes")
        if np.any(np.diff(off) < int(_MIN_ACTIVE_NODES)):
            raise ValueError("every MT arm needs at least two allocated nodes")
        n_mt = off.shape[0] - 1
        for label, array in (
            ("active_count", self.active_count),
            ("phase", self.phase),
            ("v_grow_um_per_s", self.v_grow_um_per_s),
            ("v_shrink_um_per_s", self.v_shrink_um_per_s),
        ):
            if np.asarray(array).shape != (n_mt,):
                raise ValueError(f"{label} must have shape ({n_mt},)")
        cap = np.diff(off)
        ac = np.asarray(self.active_count, dtype=np.int64)
        if np.any(ac < int(_MIN_ACTIVE_NODES)) or np.any(ac > cap):
            raise ValueError("active_count must lie in [2, per-arm capacity]")
        phase_vals = set(np.asarray(self.phase, dtype=np.int64).tolist())
        if not phase_vals <= {int(p) for p in MicrotubulePhase}:
            raise ValueError("phase entries must be MicrotubulePhase values")
        if np.any(np.asarray(self.v_grow_um_per_s) < 0.0) or np.any(np.asarray(self.v_shrink_um_per_s) < 0.0):
            raise ValueError("dynamic-instability speeds must be nonnegative")
        if np.asarray(self.mtoc_position, dtype=np.float64).shape != (3,):
            raise ValueError("mtoc_position must be a 3-vector")
        if not np.isfinite(self.k_hub_pn_per_um) or self.k_hub_pn_per_um <= 0.0:
            raise ValueError("k_hub_pn_per_um must be finite and positive")
        if not np.isfinite(self.anchor_rest_um) or self.anchor_rest_um < 0.0:
            raise ValueError("anchor_rest_um must be finite and nonnegative")
        if not np.isfinite(self.seg_um) or self.seg_um <= 0.0:
            raise ValueError("seg_um must be finite and positive")


class MicrotubuleRodRuntime:
    """CUDA-resident MT rod owner that launches real MTOC-anchor + dynamic-instability kernels.

    The runtime satisfies :class:`MicrotubuleTransaction`, so a :class:`MicrotubuleRigStateOwner` can adopt it
    as the delegate that snapshots/rolls back/commits every authoritative topology array — binding these real
    kernels through the existing seam with no private host state path.  Bending is reused through an injected
    :class:`MechanicsContributor` delegate; this runtime never re-derives bending force.
    """

    component_name: str = MICROTUBULE_COMPONENT

    def __init__(self, spec: MicrotubuleRodBuildSpec, *, device: str | None = None) -> None:
        dev = wp.get_device(device)
        if not dev.is_cuda:
            raise RuntimeError("MicrotubuleRodRuntime is Warp-CUDA-only (I0-A)")
        self.device = str(dev)
        self.n_mt = spec.n_mt
        self.n_nodes = spec.n_nodes
        self.k_hub = float(spec.k_hub_pn_per_um)
        self.anchor_rest = float(spec.anchor_rest_um)
        self.seg_um = float(spec.seg_um)
        plus_end = spec.plus_end_node()

        def f64_1d(values: npt.ArrayLike) -> wp.array:
            return wp.array(np.ascontiguousarray(values, np.float64), dtype=wp.float64, device=self.device)

        def i32_1d(values: npt.ArrayLike) -> wp.array:
            return wp.array(np.ascontiguousarray(values, np.int32), dtype=wp.int32, device=self.device)

        with wp.ScopedDevice(self.device):
            self.position_d = wp.array(
                np.ascontiguousarray(spec.position, np.float64), dtype=wp.vec3d, device=self.device
            )
            self.force_d = wp.zeros(self.n_nodes, dtype=wp.vec3d, device=self.device)
            self.mtoc_position_d = wp.array(
                np.ascontiguousarray(spec.mtoc_position, np.float64).reshape(1, 3),
                dtype=wp.vec3d,
                device=self.device,
            )
            self.mtoc_force_d = wp.zeros(1, dtype=wp.vec3d, device=self.device)
            self.fiber_offset_d = i32_1d(spec.fiber_offset)
            self.active_count_d = i32_1d(spec.active_count)
            self.plus_end_node_d = i32_1d(plus_end)
            self.phase_d = i32_1d(spec.phase)
            self.partial_length_d = f64_1d(np.zeros(self.n_mt))
            self.topology_epoch_d = i32_1d(np.zeros(1))
            self.candidate_active_count_d = i32_1d(spec.active_count)
            self.candidate_plus_end_d = i32_1d(plus_end)
            self.candidate_phase_d = i32_1d(spec.phase)
            self.candidate_partial_length_d = f64_1d(np.zeros(self.n_mt))
            self.v_grow_d = f64_1d(spec.v_grow_um_per_s)
            self.v_shrink_d = f64_1d(spec.v_shrink_um_per_s)
            self.anchor_load_d = f64_1d(np.zeros(self.n_mt))
            self.anchor_energy_d = f64_1d(np.zeros(self.n_mt))

    def geometry(self) -> MicrotubuleRigView:
        """Return a non-owning device view (identical schema to :meth:`MicrotubuleRigStateOwner.geometry`)."""
        return MicrotubuleRigView(
            position_d=self.position_d,
            force_d=self.force_d,
            mtoc_position_d=self.mtoc_position_d,
            mtoc_force_d=self.mtoc_force_d,
            fiber_offset_d=self.fiber_offset_d,
            active_count_d=self.active_count_d,
            plus_end_node_d=self.plus_end_node_d,
            phase_d=self.phase_d,
            topology_epoch_d=self.topology_epoch_d,
            n_mt=self.n_mt,
        )

    def zero_force(self) -> None:
        """Clear the rod and MTOC force accumulators on device."""
        self.force_d.zero_()
        self.mtoc_force_d.zero_()

    def accumulate_bending(self, bending: MechanicsContributor) -> None:
        """Reuse the landed Cytosim bending kernel through an injected delegate (no re-derivation)."""
        bending.accumulate(self.position_d, self.force_d)

    def accumulate_mtoc_anchor(self) -> None:
        """Launch the Newton-3rd MTOC minus-end anchor reaction against owned rod + MTOC forces."""
        wp.launch(
            _mtoc_anchor_force_kernel,
            dim=self.n_mt,
            inputs=[
                self.position_d,
                self.force_d,
                self.mtoc_position_d,
                self.mtoc_force_d,
                self.fiber_offset_d,
                wp.float64(self.k_hub),
                wp.float64(self.anchor_rest),
            ],
            outputs=[self.anchor_load_d, self.anchor_energy_d],
            device=self.device,
        )

    def propose_dynamic_instability(self, dt_phys: float) -> None:
        """Evaluate growth/shrink proposals into candidate topology arrays (committed state untouched)."""
        if not np.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        wp.launch(
            _mt_di_proposal_kernel,
            dim=self.n_mt,
            inputs=[
                self.phase_d,
                self.active_count_d,
                self.plus_end_node_d,
                self.partial_length_d,
                self.fiber_offset_d,
                self.v_grow_d,
                self.v_shrink_d,
                wp.float64(self.seg_um),
                wp.float64(dt_phys),
            ],
            outputs=[
                self.candidate_phase_d,
                self.candidate_active_count_d,
                self.candidate_plus_end_d,
                self.candidate_partial_length_d,
            ],
            device=self.device,
        )

    def propose_dynamic_instability_stochastic(
        self, dt_phys: float, rate_card: DynamicInstabilityRateCard, rng_seed: int
    ) -> None:
        """Evaluate stochastic catastrophe/rescue + growth/shrink proposals into candidate topology arrays.

        This is the PI-GAP-gated dynamic-instability path: the mechanism is bound, but the catastrophe/rescue
        rates are supplied by ``rate_card`` (an UNSET slot until PI ratifies).  Committed state is untouched —
        acceptance still runs through :meth:`commit_irreversible`.  A zero-rate card reproduces the
        deterministic :meth:`propose_dynamic_instability` bit-for-bit.
        """
        if not np.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if rng_seed < 0:
            raise ValueError("rng_seed must be nonnegative")
        wp.launch(
            _mt_di_stochastic_proposal_kernel,
            dim=self.n_mt,
            inputs=[
                self.phase_d,
                self.active_count_d,
                self.plus_end_node_d,
                self.partial_length_d,
                self.fiber_offset_d,
                self.v_grow_d,
                self.v_shrink_d,
                wp.float64(rate_card.f_catastrophe_per_s),
                wp.float64(rate_card.f_rescue_per_s),
                wp.float64(self.seg_um),
                wp.float64(dt_phys),
                wp.int32(rng_seed),
            ],
            outputs=[
                self.candidate_phase_d,
                self.candidate_active_count_d,
                self.candidate_plus_end_d,
                self.candidate_partial_length_d,
            ],
            device=self.device,
        )

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Authoritative committed topology arrays covered by snapshot/rollback/commit."""
        return (
            self.position_d,
            self.mtoc_position_d,
            self.active_count_d,
            self.plus_end_node_d,
            self.phase_d,
            self.topology_epoch_d,
        )

    def snapshot_candidate(self) -> None:
        """Reset candidate topology to the committed baseline by device-to-device copy."""
        wp.copy(self.candidate_phase_d, self.phase_d)
        wp.copy(self.candidate_active_count_d, self.active_count_d)
        wp.copy(self.candidate_plus_end_d, self.plus_end_node_d)
        wp.copy(self.candidate_partial_length_d, self.partial_length_d)

    def rollback(self, accepted: wp.array) -> None:
        """Discard scratch topology under a rejected outer-step predicate."""
        wp.launch(
            _mt_topology_rollback_kernel,
            dim=self.n_mt,
            inputs=[
                accepted,
                self.phase_d,
                self.active_count_d,
                self.plus_end_node_d,
                self.partial_length_d,
                self.candidate_phase_d,
                self.candidate_active_count_d,
                self.candidate_plus_end_d,
                self.candidate_partial_length_d,
            ],
            device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit proposed topology and advance the epoch only when ``accepted[0] != 0``.

        ``dt_phys``/``rng_seed`` complete the common transaction hook; catastrophe/rescue RNG consumes the
        seed in a later kinetic layer (plan §4.3) — no biological rate is invented here.
        """
        if not np.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if rng_seed < 0:
            raise ValueError("rng_seed must be nonnegative")
        wp.launch(
            _mt_topology_commit_kernel,
            dim=self.n_mt,
            inputs=[
                accepted,
                self.candidate_phase_d,
                self.candidate_active_count_d,
                self.candidate_plus_end_d,
                self.candidate_partial_length_d,
                self.phase_d,
                self.active_count_d,
                self.plus_end_node_d,
                self.partial_length_d,
                self.topology_epoch_d,
            ],
            device=self.device,
        )


__all__ = [
    "CORTEX_COMPONENT",
    "CYTOSOL_COMPONENT",
    "DynamicInstabilityRateCard",
    "ExternalMechanicalEndpoint",
    "LegacyMicrotubuleBendingAdapter",
    "MICROTUBULE_COMPONENT",
    "MT_CORTEX_CAPTURE",
    "MT_CYTOSOL_TRANSFER",
    "MT_NUCLEUS_LINC",
    "MT_SF_SPECTRAPLAKIN",
    "MicrotubuleGraphConnector",
    "MicrotubuleMotorPort",
    "MicrotubulePhase",
    "MicrotubuleRig",
    "MicrotubuleRigStateOwner",
    "MicrotubuleRigView",
    "MicrotubuleRodBuildSpec",
    "MicrotubuleRodRuntime",
    "MicrotubuleStepBindings",
    "MicrotubuleTransaction",
    "MtocAnchorReference",
    "NUCLEUS_COMPONENT",
    "SF_COMPONENT",
    "STATIC_ADAPTER_STATUS",
    "dynamic_instability_reference",
    "dynamic_instability_stochastic_reference",
    "mtoc_anchor_reference",
]
