"""Graph-owned port seam for the explicit head-resolved NMII actuator.

Production NMII is an explicit Stam-Hocky bipolar backbone with individual heads, per-head Bell detachment,
and per-head Hill stepping.  This module fixes the ownership boundary needed by the composed cell engine:
internal backbone/head-arm state belongs to ``nmii`` while every head-to-actin crossbridge belongs to a
graph-owned MOTOR connector targeting a live filament port.

The landed :class:`ac.motor.minifilament_warp.MyosinForce` remains a mechanistic reference and transition
runtime.  Its segment binding state and global indices are still embedded, so the strict legacy adapter below
cannot be installed as a component-local production backend.  Aggregate and two-anchor models remain
diagnostic-only.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Callable, Protocol, runtime_checkable

import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentContract,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
)
from aleph.engine.runtime import LedgerContributor, TransactionParticipant

# ── Real ac/motor Warp kernels bound behind this component seam (KERNEL_BOUND). ──────────────────
# These are the landed, mechanistic Stam-Hocky primitives; this module re-owns their layout, it never
# reimplements the physics.  Import-time only builds the Kernel objects (lazy JIT), so the dev Mac stays
# CPU-green; every launch flows through an injected launcher (``wp.launch`` in production).
from aleph.components.motor.backbone_warp import (
    KBT_PN_UM,
    angle_harmonic_kernel,
    arm_orientation_k_theta,
    backbone_bending_k_theta,
    backbone_bending_kappa,
)
from aleph.components.motor.minifilament_warp import (
    ARM_REST_ANGLE,
    BACKBONE_REST_ANGLE,
    harmonic_bond_kernel,
)
from aleph.components.motor.segment_motor import (
    _increment_epoch_if_accepted_kernel,
    _refresh_bound_walk_dir_kernel,
    attach_segment_gated_kernel,
    compute_head_loads_segment_split_kernel,
    crossbridge_segment_split_kernel,
    refresh_segment_barbed_kernel,
    step_detach_segment_gated_kernel,
)

# One launch/copy indirection so the binding classes are exercisable on a CUDA-free host with a recording
# double, while production passes the real Warp primitives.  These are NOT a private state path: every array
# they touch is an injected authoritative component/port/connector array, never a hidden copy.
LaunchFn = Callable[..., None]
CopyFn = Callable[..., None]


def _default_launch(kernel: object, *, dim: int, inputs: list, outputs: list | None = None,
                    device: object | None = None) -> None:
    """Production launcher: forward to ``wp.launch`` on the owning CUDA device."""
    if outputs is None:
        wp.launch(kernel, dim=dim, inputs=inputs, device=device)
    else:
        wp.launch(kernel, dim=dim, inputs=inputs, outputs=outputs, device=device)


def _default_copy(dst: object, src: object) -> None:
    """Production D2D snapshot/restore: forward to ``wp.copy`` (reused Warp primitive, not a kernel)."""
    wp.copy(dst, src)

NMII_COMPONENT = "nmii"
SF_COMPONENT = "sf_arc"
CORTEX_COMPONENT = "cortex"
LAMELLIPODIUM_COMPONENT = "lamellipodium"
FILOPODIUM_COMPONENT = "filopodium"

NMII_SF_MOTOR = "nmii_sf_motor"
NMII_CORTEX_MOTOR = "nmii_cortex_motor"
NMII_LAMELLIPODIUM_MOTOR = "nmii_lamellipodium_motor"
NMII_FILOPODIUM_MOTOR = "nmii_filopodium_motor"

NMII_REPRESENTATION = "explicit Stam-Hocky bipolar backbone + individual heads"
NMII_SOLVER = "head-resolved internal mechanics plus graph MOTOR KMC"
NMII_STATE_EVENT_CHANNELS = frozenset({"minifilament_turnover", "atp_cycle"})
NMII_CONNECTOR_EVENT_CHANNELS = frozenset({"stepping", "binding", "unbinding"})
NMII_LEDGER_CHANNELS = frozenset({"force", "work", "atp", "population"})


class LegacyNMIIKind(StrEnum):
    """Permitted labels for quarantined reference/diagnostic delegates."""

    HEAD_RESOLVED_MYOSIN_FORCE = "head_resolved_myosin_force"
    AGGREGATE_FORCE_DIPOLE = "aggregate_force_dipole"
    TWO_ANCHOR_LINEAR = "two_anchor_linear"


def _storage_key(array: object) -> tuple[object, ...]:
    """Return device storage identity without reading authoritative values."""
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
    """Validate a non-empty one-dimensional Warp CUDA array without a readback."""
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


def _require_same_device(arrays: tuple[object, ...], *, label: str) -> None:
    if len({str(getattr(array, "device", None)) for array in arrays}) != 1:
        raise ValueError(f"{label} arrays must share one CUDA device")


def _require_distinct_storage(arrays: tuple[object, ...], *, label: str) -> None:
    keys = tuple(_storage_key(array) for array in arrays)
    if len(set(keys)) != len(keys):
        raise ValueError(f"{label} arrays must own distinct storage")


@runtime_checkable
class NMIIInternalMechanics(Protocol):
    """Component-local production mechanics for explicit backbone and head arms."""

    n_particles: int
    n_heads: int
    explicit_stam_hocky_backbone: bool
    individual_heads: bool
    graph_owned_motor_connectors: bool
    component_local: bool
    production_eligible: bool
    aggregate_or_two_anchor: bool

    def accumulate_internal(self, position: wp.array, force: wp.array) -> None:
        """Add only backbone/head-arm mechanics; crossbridges remain graph-owned."""


@runtime_checkable
class NMIITransaction(TransactionParticipant, Protocol):
    """Transaction delegate for actuator-owned geometry, population, and ATP state."""

    component_name: str
    event_channels: frozenset[str]

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return all mutable authoritative arrays covered by snapshot/rollback."""


@runtime_checkable
class NMIILedgerContributor(LedgerContributor, Protocol):
    """Ledger contributor covering all required motor accounting channels."""

    ledger_channels: frozenset[str]


@dataclass(frozen=True, slots=True)
class NMIIActuatorView:
    """Non-owning explicit topology/load view exposed to graph-owned MOTOR connectors."""

    position_d: wp.array
    force_d: wp.array
    minifilament_offset_d: wp.array
    active_minifilament_d: wp.array
    minifilament_id_d: wp.array
    particle_role_d: wp.array
    head_node_d: wp.array
    head_id_d: wp.array
    head_side_d: wp.array
    head_load_hill_d: wp.array
    head_load_bell_d: wp.array
    actuator_epoch_d: wp.array
    n_minifilaments: int
    n_heads: int


@dataclass(frozen=True, slots=True)
class FilamentMotorPortView:
    """Target-owned live actin-segment port addressed by persistent material identity."""

    component: str
    position_d: wp.array
    force_d: wp.array
    segment_node_a_d: wp.array
    segment_node_b_d: wp.array
    segment_polarity_d: wp.array
    persistent_filament_id_d: wp.array
    material_s0_d: wp.array
    material_s1_d: wp.array
    topology_epoch_d: wp.array

    def __post_init__(self) -> None:
        if self.component not in _TARGET_COMPONENTS:
            raise ValueError(f"unsupported NMII filament port component {self.component!r}")
        _validate_device_array(
            self.position_d,
            label=f"{self.component}.position_d",
            dtype=wp.vec3d,
        )
        n_nodes = int(self.position_d.shape[0])
        _validate_device_array(
            self.force_d,
            label=f"{self.component}.force_d",
            dtype=wp.vec3d,
            length=n_nodes,
        )
        _validate_device_array(
            self.segment_node_a_d,
            label=f"{self.component}.segment_node_a_d",
            dtype=wp.int32,
        )
        n_segments = int(self.segment_node_a_d.shape[0])
        for label, array, dtype in (
            ("segment_node_b_d", self.segment_node_b_d, wp.int32),
            ("segment_polarity_d", self.segment_polarity_d, wp.int32),
            ("persistent_filament_id_d", self.persistent_filament_id_d, wp.int64),
            ("material_s0_d", self.material_s0_d, wp.float64),
            ("material_s1_d", self.material_s1_d, wp.float64),
        ):
            _validate_device_array(
                array,
                label=f"{self.component}.{label}",
                dtype=dtype,
                length=n_segments,
            )
        _validate_device_array(
            self.topology_epoch_d,
            label=f"{self.component}.topology_epoch_d",
            dtype=wp.int32,
            length=1,
        )
        arrays = (
            self.position_d,
            self.force_d,
            self.segment_node_a_d,
            self.segment_node_b_d,
            self.segment_polarity_d,
            self.persistent_filament_id_d,
            self.material_s0_d,
            self.material_s1_d,
            self.topology_epoch_d,
        )
        _require_same_device(arrays, label=f"{self.component} motor-port")
        _require_distinct_storage(arrays, label=f"{self.component} motor-port")


@runtime_checkable
class NMIIMotorConnector(Protocol):
    """Graph-owned per-head binding/stepping state for one target actin actor."""

    name: str
    component_a: str
    component_b: str
    event_channels: frozenset[str]
    ledger_channels: frozenset[str]
    individual_head_state: bool
    per_head_bell: bool
    hill_force_velocity: bool
    aggregate_or_two_anchor: bool

    def accumulate_candidate(
        self,
        actuator: NMIIActuatorView,
        port: FilamentMotorPortView,
    ) -> None:
        """Apply crossbridges and adjoint target reactions; write candidate per-head loads."""

    def snapshot_candidate(self) -> None:
        """Snapshot binding, material coordinate, abscissa, polarity, and RNG epoch."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore connector-owned state under rejection."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit Hill stepping and Bell bind/unbind once after global acceptance."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute force/work/ATP/bound-population accounting."""


@dataclass(frozen=True, slots=True)
class NMIIActuatorStateOwner:
    """Authoritative component-local CUDA state for explicit NMII minifilaments."""

    n_minifilaments: int
    n_heads: int
    position_d: wp.array
    force_d: wp.array
    minifilament_offset_d: wp.array
    active_minifilament_d: wp.array
    minifilament_id_d: wp.array
    particle_role_d: wp.array
    head_node_d: wp.array
    head_id_d: wp.array
    head_side_d: wp.array
    head_load_hill_d: wp.array
    head_load_bell_d: wp.array
    atp_cycle_state_d: wp.array
    atp_consumed_d: wp.array
    actuator_epoch_d: wp.array
    mechanics: NMIIInternalMechanics
    transaction: NMIITransaction
    ledger: NMIILedgerContributor

    def __post_init__(self) -> None:
        if self.n_minifilaments <= 0 or self.n_heads <= 0:
            raise ValueError("NMII minifilament and head counts must be positive")
        _validate_device_array(self.position_d, label="nmii.position_d", dtype=wp.vec3d)
        n_particles = int(self.position_d.shape[0])
        _validate_device_array(
            self.force_d,
            label="nmii.force_d",
            dtype=wp.vec3d,
            length=n_particles,
        )
        _validate_device_array(
            self.minifilament_offset_d,
            label="nmii.minifilament_offset_d",
            dtype=wp.int32,
            length=self.n_minifilaments + 1,
        )
        _validate_device_array(
            self.active_minifilament_d,
            label="nmii.active_minifilament_d",
            dtype=wp.int32,
            length=self.n_minifilaments,
        )
        _validate_device_array(
            self.minifilament_id_d,
            label="nmii.minifilament_id_d",
            dtype=wp.int64,
            length=self.n_minifilaments,
        )
        _validate_device_array(
            self.particle_role_d,
            label="nmii.particle_role_d",
            dtype=wp.int32,
            length=n_particles,
        )
        for label, array, dtype in (
            ("head_node_d", self.head_node_d, wp.int32),
            ("head_id_d", self.head_id_d, wp.int64),
            ("head_side_d", self.head_side_d, wp.int32),
            ("head_load_hill_d", self.head_load_hill_d, wp.float64),
            ("head_load_bell_d", self.head_load_bell_d, wp.float64),
            ("atp_cycle_state_d", self.atp_cycle_state_d, wp.int32),
        ):
            _validate_device_array(
                array,
                label=f"nmii.{label}",
                dtype=dtype,
                length=self.n_heads,
            )
        _validate_device_array(
            self.atp_consumed_d,
            label="nmii.atp_consumed_d",
            dtype=wp.int64,
            length=1,
        )
        _validate_device_array(
            self.actuator_epoch_d,
            label="nmii.actuator_epoch_d",
            dtype=wp.int32,
            length=1,
        )
        arrays = (
            self.position_d,
            self.force_d,
            self.minifilament_offset_d,
            self.active_minifilament_d,
            self.minifilament_id_d,
            self.particle_role_d,
            self.head_node_d,
            self.head_id_d,
            self.head_side_d,
            self.head_load_hill_d,
            self.head_load_bell_d,
            self.atp_cycle_state_d,
            self.atp_consumed_d,
            self.actuator_epoch_d,
        )
        _require_same_device(arrays, label="NMII actuator")
        _require_distinct_storage(arrays, label="NMII actuator")
        _validate_production_mechanics(self.mechanics, n_particles=n_particles, n_heads=self.n_heads)
        if self.transaction.component_name != NMII_COMPONENT:
            raise ValueError("NMII transaction delegate must be owned by nmii")
        missing_events = NMII_STATE_EVENT_CHANNELS - self.transaction.event_channels
        if missing_events:
            raise ValueError(f"NMII transaction is missing event channels {sorted(missing_events)}")
        mutable = (
            self.position_d,
            self.active_minifilament_d,
            self.atp_cycle_state_d,
            self.atp_consumed_d,
            self.actuator_epoch_d,
        )
        covered = {_storage_key(array) for array in self.transaction.owned_arrays()}
        if {_storage_key(array) for array in mutable} - covered:
            raise ValueError("NMII transaction must cover geometry, population, ATP state, and epoch")
        _validate_ledger_channels(self.ledger, label="NMII state ledger")

    def geometry(self) -> NMIIActuatorView:
        """Return a non-owning explicit topology/load view for MOTOR connectors."""
        return NMIIActuatorView(
            position_d=self.position_d,
            force_d=self.force_d,
            minifilament_offset_d=self.minifilament_offset_d,
            active_minifilament_d=self.active_minifilament_d,
            minifilament_id_d=self.minifilament_id_d,
            particle_role_d=self.particle_role_d,
            head_node_d=self.head_node_d,
            head_id_d=self.head_id_d,
            head_side_d=self.head_side_d,
            head_load_hill_d=self.head_load_hill_d,
            head_load_bell_d=self.head_load_bell_d,
            actuator_epoch_d=self.actuator_epoch_d,
            n_minifilaments=self.n_minifilaments,
            n_heads=self.n_heads,
        )

    def accumulate_internal(self) -> None:
        """Launch only component-owned backbone/head-arm candidate mechanics."""
        self.mechanics.accumulate_internal(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Push the delegate's declared channels, then this body's own force resultant.

        The ``nmii`` particle array is the second of the two never-merged bodies a MOTOR connector cuts
        the slice into, so its whole-array resultant is the ``traction`` half of the force-balance gate —
        equal and opposite to the target body's half whenever the crossbridge's adjoint scatter is
        correct.  See :meth:`~aleph.engine.ledger.GlobalCellLedger.add_body_force`.

        A ledger exposing no such accessor keeps the previous delegate-only behaviour: the channel is
        simply not filled, and the caller's gate is then assembled from whatever else it holds.  This is
        the one place that stays permissive, because the delegate seam is exercised by recording doubles
        that own no accumulator at all.
        """
        self.ledger.accumulate_ledger(ledger)
        add_body_force = getattr(ledger, "add_body_force", None)
        if callable(add_body_force):
            add_body_force(self.force_d, side="traction")


def _validate_production_mechanics(
    mechanics: object,
    *,
    n_particles: int,
    n_heads: int,
) -> None:
    required_true = (
        "explicit_stam_hocky_backbone",
        "individual_heads",
        "graph_owned_motor_connectors",
        "component_local",
        "production_eligible",
    )
    missing = [name for name in required_true if getattr(mechanics, name, False) is not True]
    if missing:
        raise ValueError(f"NMII production mechanics is missing fidelity guarantees {missing}")
    if getattr(mechanics, "aggregate_or_two_anchor", True) is not False:
        raise ValueError("aggregate or two-anchor NMII mechanics is diagnostic-only")
    if int(getattr(mechanics, "n_particles", -1)) != n_particles:
        raise ValueError("NMII mechanics particle count does not match actuator state")
    if int(getattr(mechanics, "n_heads", -1)) != n_heads:
        raise ValueError("NMII mechanics head count does not match actuator state")


def _validate_ledger_channels(contributor: object, *, label: str) -> None:
    channels = frozenset(getattr(contributor, "ledger_channels", ()))
    missing = NMII_LEDGER_CHANNELS - channels
    if missing:
        raise ValueError(f"{label} is missing channels {sorted(missing)}")


@dataclass(frozen=True, slots=True)
class NMIIContractProposal:
    """Local common-graph proposal; it never mutates the reference architecture."""

    components: tuple[ComponentContract, ...]
    connectors: tuple[ConnectorContract, ...]


_TARGET_COMPONENTS = (
    SF_COMPONENT,
    CORTEX_COMPONENT,
    LAMELLIPODIUM_COMPONENT,
    FILOPODIUM_COMPONENT,
)
_REQUIRED_CONNECTORS = {
    NMII_SF_MOTOR: SF_COMPONENT,
    NMII_CORTEX_MOTOR: CORTEX_COMPONENT,
    NMII_LAMELLIPODIUM_MOTOR: LAMELLIPODIUM_COMPONENT,
    NMII_FILOPODIUM_MOTOR: FILOPODIUM_COMPONENT,
}


def nmii_contract_proposal() -> NMIIContractProposal:
    """Return the missing actuator component and graph-owned actin MOTOR edges."""
    component = ComponentContract(
        NMII_COMPONENT,
        ComponentRole.ACTIVE_LOAD_PATH,
        NMII_REPRESENTATION,
        NMII_SOLVER,
    )
    connectors = tuple(
        ConnectorContract(
            name,
            ConnectorFamily.MOTOR,
            NMII_COMPONENT,
            target,
            True,
            True,
            endpoint_role_a="individual NMII head crossbridge",
            endpoint_role_b="live polar actin material coordinate",
            chemistry_card="nmii_head_actin_hill_bell",
            generation_required=True,
            remap_on_accept=True,
            blocks_sleep_refine=True,
        )
        for name, target in _REQUIRED_CONNECTORS.items()
    )
    return NMIIContractProposal(components=(component,), connectors=connectors)


def proposed_nmii_architecture(base: CellArchitecture) -> CellArchitecture:
    """Return an idempotent graph proposal after protrusion actors have been registered."""
    known = {component.name for component in base.components}
    missing_targets = set(_TARGET_COMPONENTS) - known
    if missing_targets:
        raise ValueError(
            "NMII graph requires target components before MOTOR edges: "
            f"{sorted(missing_targets)}"
        )
    proposal = nmii_contract_proposal()
    connector_names = {connector.name for connector in base.connectors}
    architecture = CellArchitecture(
        components=(
            *base.components,
            *(component for component in proposal.components if component.name not in known),
        ),
        connectors=(
            *base.connectors,
            *(connector for connector in proposal.connectors if connector.name not in connector_names),
        ),
    )
    _validate_architecture(architecture)
    return architecture


def _find_connector(architecture: CellArchitecture, name: str) -> ConnectorContract:
    for connector in architecture.connectors:
        if connector.name == name:
            return connector
    raise ValueError(f"architecture must register proposed connector {name!r}")


def _validate_architecture(architecture: CellArchitecture) -> None:
    component = architecture.component(NMII_COMPONENT)
    if component.role is not ComponentRole.ACTIVE_LOAD_PATH:
        raise ValueError("nmii must have ComponentRole.ACTIVE_LOAD_PATH")
    if component.representation != NMII_REPRESENTATION:
        raise ValueError("nmii cannot use an aggregate or two-anchor representation")
    if not component.owns_geometry or not component.dynamically_evolving:
        raise ValueError("nmii must own evolving explicit backbone/head geometry")
    for name, target in _REQUIRED_CONNECTORS.items():
        connector = _find_connector(architecture, name)
        if connector.family is not ConnectorFamily.MOTOR:
            raise ValueError(f"{name} must use ConnectorFamily.MOTOR")
        if {connector.component_a, connector.component_b} != {NMII_COMPONENT, target}:
            raise ValueError(f"{name} has incorrect component endpoints")
        if not connector.kinetics or not connector.commit_on_accept:
            raise ValueError(f"{name} must be kinetic and accepted-step committed")
        if not connector.bidirectional or not connector.adjoint_transfer_required:
            raise ValueError(f"{name} must be bidirectional with adjoint transfer")
        if not (
            connector.generation_required
            and connector.remap_on_accept
            and connector.blocks_sleep_refine
        ):
            raise ValueError(f"{name} must enforce generation/remap/sleep-refine safety")


@dataclass(frozen=True, slots=True)
class NMIIActuatorStepBindings:
    """Exactly one target-owned filament port and graph-owned MOTOR runtime per required edge."""

    ports: tuple[FilamentMotorPortView, ...]
    connectors: tuple[NMIIMotorConnector, ...]

    def __post_init__(self) -> None:
        port_by_component: dict[str, FilamentMotorPortView] = {}
        for port in self.ports:
            if port.component in port_by_component:
                raise ValueError(f"duplicate NMII filament port {port.component!r}")
            port_by_component[port.component] = port
        if set(port_by_component) != set(_TARGET_COMPONENTS):
            raise ValueError(
                "NMII filament ports mismatch; "
                f"expected={sorted(_TARGET_COMPONENTS)}, got={sorted(port_by_component)}"
            )
        connector_by_name: dict[str, NMIIMotorConnector] = {}
        for connector in self.connectors:
            if connector.name in connector_by_name:
                raise ValueError(f"duplicate NMII MOTOR connector {connector.name!r}")
            connector_by_name[connector.name] = connector
        if set(connector_by_name) != set(_REQUIRED_CONNECTORS):
            raise ValueError(
                "NMII MOTOR bindings mismatch; "
                f"expected={sorted(_REQUIRED_CONNECTORS)}, got={sorted(connector_by_name)}"
            )
        for name, target in _REQUIRED_CONNECTORS.items():
            connector = connector_by_name[name]
            if {connector.component_a, connector.component_b} != {NMII_COMPONENT, target}:
                raise ValueError(f"runtime connector {name!r} has incorrect component endpoints")
            missing_events = NMII_CONNECTOR_EVENT_CHANNELS - connector.event_channels
            if missing_events:
                raise ValueError(f"runtime connector {name!r} is missing events {sorted(missing_events)}")
            required_true = ("individual_head_state", "per_head_bell", "hill_force_velocity")
            missing_fidelity = [
                guarantee
                for guarantee in required_true
                if getattr(connector, guarantee, False) is not True
            ]
            if missing_fidelity:
                raise ValueError(
                    f"runtime connector {name!r} is missing per-head fidelity {missing_fidelity}"
                )
            if getattr(connector, "aggregate_or_two_anchor", True) is not False:
                raise ValueError(f"runtime connector {name!r} cannot use aggregate/two-anchor mechanics")
            _validate_ledger_channels(connector, label=f"runtime connector {name!r}")

    def pair_for(self, name: str) -> tuple[NMIIMotorConnector, FilamentMotorPortView]:
        """Return one validated connector/target-port pair."""
        target = _REQUIRED_CONNECTORS[name]
        connector = next(connector for connector in self.connectors if connector.name == name)
        port = next(port for port in self.ports if port.component == target)
        return connector, port


@dataclass(frozen=True, slots=True)
class NMIIActuator:
    """Stateless facade over the actuator owner and four graph-owned MOTOR edges."""

    architecture: CellArchitecture
    state: NMIIActuatorStateOwner

    def __post_init__(self) -> None:
        _validate_architecture(self.architecture)

    def accumulate_candidate(self, bindings: NMIIActuatorStepBindings) -> None:
        """Assemble internal candidate mechanics then every graph-owned crossbridge family."""
        self.state.accumulate_internal()
        actuator = self.state.geometry()
        for name in _REQUIRED_CONNECTORS:
            connector, port = bindings.pair_for(name)
            connector.accumulate_candidate(actuator, port)

    def transaction_participants(
        self,
        bindings: NMIIActuatorStepBindings,
    ) -> tuple[TransactionParticipant, ...]:
        """Return component plus every port connector in deterministic graph order."""
        return (
            self.state,
            *(bindings.pair_for(name)[0] for name in _REQUIRED_CONNECTORS),
        )

    def snapshot_candidate(self, bindings: NMIIActuatorStepBindings) -> None:
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: NMIIActuatorStepBindings, accepted: wp.array) -> None:
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self,
        bindings: NMIIActuatorStepBindings,
        accepted: wp.array,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        """Commit actuator population/ATP and per-port Hill/Bell events once on acceptance."""
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, bindings: NMIIActuatorStepBindings, ledger: object) -> None:
        """Collect force, work, ATP, and population ledgers without owning the verdict."""
        self.state.accumulate_ledger(ledger)
        for name in _REQUIRED_CONNECTORS:
            bindings.pair_for(name)[0].accumulate_ledger(ledger)


class _LegacyNMIIForce(Protocol):
    backbone_bonds: wp.array
    head_bonds: wp.array
    head_node: wp.array
    state: dict[str, wp.array]
    loads: wp.array
    loads_full: wp.array
    params: object
    segment_runtime: object | None

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the legacy/reference motor force."""


@dataclass(frozen=True, slots=True)
class LegacyNMIIAdapter:
    """Strict diagnostic/reference adapter; never a production component backend."""

    delegate: _LegacyNMIIForce
    kind: LegacyNMIIKind

    def __post_init__(self) -> None:
        if self.kind is LegacyNMIIKind.HEAD_RESOLVED_MYOSIN_FORCE:
            required = (
                "backbone_bonds",
                "head_bonds",
                "head_node",
                "state",
                "loads",
                "loads_full",
                "params",
                "segment_runtime",
            )
            missing = [name for name in required if not hasattr(self.delegate, name)]
            if missing:
                raise ValueError(f"head-resolved legacy delegate is missing {missing}")
            if self.delegate.segment_runtime is None:
                raise ValueError("head-resolved legacy adapter requires production segment anchoring")
            state = self.delegate.state
            state_keys = {"bound", "seg_id", "seg_a", "seg_b", "bary_t", "abscissa", "walk_dir"}
            if not isinstance(state, dict) or state_keys - set(state):
                raise ValueError("head-resolved legacy delegate lacks individual segment-bound head state")
            n_heads = int(self.delegate.head_node.shape[0])
            for label in ("bound", "abscissa", "walk_dir"):
                if int(state[label].shape[0]) != n_heads:
                    raise ValueError("legacy per-head state count does not match head topology")
            if int(self.delegate.loads.shape[0]) != n_heads:
                raise ValueError("legacy Hill load count does not match head topology")
            if int(self.delegate.loads_full.shape[0]) != n_heads:
                raise ValueError("legacy Bell load count does not match head topology")

    @property
    def component_local(self) -> bool:
        return False

    @property
    def graph_owned_motor_connectors(self) -> bool:
        return False

    @property
    def production_eligible(self) -> bool:
        return False

    @property
    def aggregate_or_two_anchor(self) -> bool:
        return self.kind is not LegacyNMIIKind.HEAD_RESOLVED_MYOSIN_FORCE

    @property
    def status(self) -> str:
        if self.kind is LegacyNMIIKind.HEAD_RESOLVED_MYOSIN_FORCE:
            return "HEAD_RESOLVED_MONOLITHIC_CONNECTOR_SPLIT_PENDING"
        return f"DIAGNOSTIC_ONLY_{self.kind.value.upper()}"

    def accumulate_reference(self, pos: wp.array, force: wp.array) -> None:
        """Delegate an isolated reference/diagnostic call on the legacy combined arrays."""
        self.delegate.accumulate(pos, force)


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# KERNEL_BOUND bindings — the real ac/motor Warp kernels split behind the component/connector ownership
# boundary.  Internal backbone/head-arm mechanics belong to ``nmii``; per-head crossbridge kinetics belong to
# each graph-owned MOTOR connector.  No aggregate/two-anchor path is bindable here (I0-A diagnostics only).
# ─────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class NMIIBendingStiffness:
    """DERIVED (never fit) angle-harmonic stiffnesses for the F6 backbone/arm transmission fix.

    Both values come from the grid-invariant rigid-rod maps in :mod:`aleph.components.motor.backbone_warp`
    (Magic-Number Block): ``k_theta_backbone = (L_p·k_BT)/a`` and ``k_theta_arm = k_xb·r0_head²``.  A
    persistence length ``L_p`` for a mature bipolar minifilament is an unresolved KB/PI GAP; supplying a
    single-coiled-coil proxy keeps this diagnostic, so ``sourced`` records whether the caller passed a
    ratified material contract.  This class only *derives*; it does not choose a number to pass a gate.
    """

    persistence_length_um: float
    segment_len_um: float
    k_xb_pn_per_um: float
    r0_head_um: float
    kbt_pn_um: float = KBT_PN_UM
    sourced: bool = False

    @property
    def k_theta_backbone(self) -> float:
        """Backbone bending stiffness ``κ/a`` [pN·µm/rad²] via the bead-angle → continuum map."""
        kappa = backbone_bending_kappa(self.persistence_length_um, self.kbt_pn_um)
        return backbone_bending_k_theta(kappa, self.segment_len_um)

    @property
    def k_theta_arm(self) -> float:
        """Head-arm orientation stiffness ``k_xb·r0_head²`` [pN·µm/rad²] (rigid-lever limit)."""
        return arm_orientation_k_theta(self.k_xb_pn_per_um, self.r0_head_um)


class BackboneArmMechanics:
    """Component-local internal NMII mechanics bound to the real ``harmonic_bond``/``angle_harmonic`` kernels.

    This is the production :class:`NMIIInternalMechanics` backend.  It assembles ONLY the backbone axial
    rod, the head↔backbone lever arm, and the F6 angle-harmonic bending — the minifilament's internal
    structure.  It deliberately does NOT touch any head-to-actin crossbridge: those are graph-owned MOTOR
    connectors on another actor's array.  That split is what the current monolithic ``MyosinForce`` cannot
    express, and re-owning it here is the KERNEL_BOUND increment.

    ``accumulate_internal`` launches the four real kernels on the *injected* authoritative
    ``position``/``force`` arrays (the actuator's own particle state).  No positions, forces, or topology are
    copied into private storage — the delegate is a pure force accumulator over arrays it is handed.
    """

    #: fidelity guarantees read by :func:`_validate_production_mechanics`.
    explicit_stam_hocky_backbone = True
    individual_heads = True
    per_head_bell = True
    hill_force_velocity = True
    graph_owned_motor_connectors = True
    component_local = True
    production_eligible = True
    aggregate_or_two_anchor = False
    #: crossbridges are NOT assembled here — they belong to the MOTOR connectors (ownership split marker).
    binds_crossbridge = False

    def __init__(
        self,
        *,
        n_particles: int,
        n_heads: int,
        backbone_bonds_d: wp.array,
        head_bonds_d: wp.array,
        backbone_angles_d: wp.array,
        head_arm_angles_d: wp.array,
        k_backbone_pn_per_um: float,
        r0_backbone_um: float,
        k_head_spring_pn_per_um: float,
        r0_head_um: float,
        bending: NMIIBendingStiffness,
        launch: LaunchFn = _default_launch,
        device: object | None = None,
    ) -> None:
        if n_particles <= 0 or n_heads <= 0:
            raise ValueError("BackboneArmMechanics needs positive particle and head counts")
        self.n_particles = int(n_particles)
        self.n_heads = int(n_heads)
        self.backbone_bonds_d = backbone_bonds_d
        self.head_bonds_d = head_bonds_d
        self.backbone_angles_d = backbone_angles_d
        self.head_arm_angles_d = head_arm_angles_d
        self.k_backbone = float(k_backbone_pn_per_um)
        self.r0_backbone = float(r0_backbone_um)
        self.k_head_spring = float(k_head_spring_pn_per_um)
        self.r0_head = float(r0_head_um)
        self.bending = bending
        self._launch = launch
        self._device = device
        self.launched_kernels: list[object] = []

    def accumulate_internal(self, position: wp.array, force: wp.array) -> None:
        """Add backbone rod + head-arm + F6 bending forces to the actuator's own ``force`` array.

        Four real ac/motor kernels flow through the seam: ``harmonic_bond_kernel`` for the backbone rod and
        the head↔backbone arm (rests ``r0_backbone`` and ``r0_head``), then ``angle_harmonic_kernel`` for the
        backbone straight-rod bending (rest π) and the head-arm perpendicular orientation (rest π/2).  No
        crossbridge kernel is launched — that is the connector's job.
        """
        self.launched_kernels = []
        self._launch(harmonic_bond_kernel, dim=int(self.backbone_bonds_d.shape[0]),
                     inputs=[position, self.backbone_bonds_d, wp.float64(self.k_backbone),
                             wp.float64(self.r0_backbone)], outputs=[force], device=self._device)
        self.launched_kernels.append(harmonic_bond_kernel)
        self._launch(harmonic_bond_kernel, dim=int(self.head_bonds_d.shape[0]),
                     inputs=[position, self.head_bonds_d, wp.float64(self.k_head_spring),
                             wp.float64(self.r0_head)], outputs=[force], device=self._device)
        self.launched_kernels.append(harmonic_bond_kernel)
        self._launch(angle_harmonic_kernel, dim=int(self.backbone_angles_d.shape[0]),
                     inputs=[position, self.backbone_angles_d, wp.float64(self.bending.k_theta_backbone),
                             wp.float64(BACKBONE_REST_ANGLE)], outputs=[force], device=self._device)
        self.launched_kernels.append(angle_harmonic_kernel)
        self._launch(angle_harmonic_kernel, dim=int(self.head_arm_angles_d.shape[0]),
                     inputs=[position, self.head_arm_angles_d, wp.float64(self.bending.k_theta_arm),
                             wp.float64(ARM_REST_ANGLE)], outputs=[force], device=self._device)
        self.launched_kernels.append(angle_harmonic_kernel)


@dataclass(frozen=True, slots=True)
class SegmentConnectorState:
    """Graph-owned per-head crossbridge binding state + snapshot for one MOTOR edge.

    Every array is owned by the connector (not the actuator and not the target): the bound flag, the target
    *segment* endpoints, barycentric fraction, walked abscissa, live walk direction, RNG epoch, plus the
    attachment-query scratch.  The ``*_snap`` twins hold the accepted-step snapshot so rejection restores the
    candidate bit-exactly.  Raw segment indices are candidate cache only; persistent identity is the target's
    material coordinate (owned by the port).
    """

    n_heads: int
    n_segments: int
    bound_d: wp.array
    seg_id_d: wp.array
    seg_a_d: wp.array
    seg_b_d: wp.array
    bary_t_d: wp.array
    abscissa_d: wp.array
    walk_dir_d: wp.array
    rng_epoch_d: wp.array
    loads_hill_d: wp.array
    loads_bell_d: wp.array
    seg_barbed_d: wp.array
    query_seg_id_d: wp.array
    query_t_d: wp.array
    query_barbed_d: wp.array
    bound_snap_d: wp.array
    seg_id_snap_d: wp.array
    seg_a_snap_d: wp.array
    seg_b_snap_d: wp.array
    bary_t_snap_d: wp.array
    abscissa_snap_d: wp.array
    walk_dir_snap_d: wp.array
    rng_epoch_snap_d: wp.array

    def authoritative_pairs(self) -> tuple[tuple[wp.array, wp.array], ...]:
        """Return ``(live, snapshot)`` pairs covered by snapshot/rollback in a deterministic order."""
        return (
            (self.bound_d, self.bound_snap_d),
            (self.seg_id_d, self.seg_id_snap_d),
            (self.seg_a_d, self.seg_a_snap_d),
            (self.seg_b_d, self.seg_b_snap_d),
            (self.bary_t_d, self.bary_t_snap_d),
            (self.abscissa_d, self.abscissa_snap_d),
            (self.walk_dir_d, self.walk_dir_snap_d),
            (self.rng_epoch_d, self.rng_epoch_snap_d),
        )


class SegmentMotorConnectorRuntime:
    """Graph-owned MOTOR connector bound to the real segment-anchored NMII KMC/geometry kernels.

    Reuses the landed :mod:`aleph.components.motor.segment_motor` kernels through the seam:

    * ``accumulate_candidate`` refreshes each bound head's live walk direction from the *target* port's
      segment geometry (``_refresh_bound_walk_dir_kernel``) and the broad-phase barbed field
      (``refresh_segment_barbed_kernel``);
    * ``accumulate_candidate`` also scatters the head-to-actin *crossbridge force* through the two-array adjoint
      kernel ``crossbridge_segment_split_kernel``: ``+f`` into the ``nmii``-owned actuator force and
      ``−(1−t)·f``/``−t·f`` into the target-owned port force, Newton-3rd across the two arrays with no merge;
    * ``compute_loads`` fills the connector-owned Hill/Bell per-head loads via
      ``compute_head_loads_segment_split_kernel`` (head node from the actuator array, segment nodes from the
      port array) — the composed step runs it after convergence and before ``commit_irreversible`` consumes
      ``loads_hill_d``/``loads_bell_d``;
    * ``commit_irreversible`` runs the accepted-predicated attach → Hill-step/Bell-detach → epoch KMC
      (``attach_segment_gated_kernel``, ``step_detach_segment_gated_kernel``,
      ``_increment_epoch_if_accepted_kernel``), so a rejected outer step advances no binding state;
    * ``snapshot_candidate``/``rollback`` D2D-copy the connector-owned authoritative arrays.

    The two-array crossbridge/load kernels are the ``ac/motor`` shared relayout the Lead landed to close the
    split-ownership seam (the head node is owned by ``nmii``, the segment nodes by the target actor): they never
    merge the two ``pos``/``force`` arrays and preserve the single-array power-stroke physics exactly.
    """

    individual_head_state = True
    per_head_bell = True
    hill_force_velocity = True
    aggregate_or_two_anchor = False

    def __init__(
        self,
        *,
        name: str,
        component_b: str,
        state: SegmentConnectorState,
        head_node_d: wp.array,
        params: object,
        event_channels: frozenset[str] = NMII_CONNECTOR_EVENT_CHANNELS,
        ledger_channels: frozenset[str] = NMII_LEDGER_CHANNELS,
        launch: LaunchFn = _default_launch,
        copy: CopyFn = _default_copy,
        device: object | None = None,
    ) -> None:
        if name not in _REQUIRED_CONNECTORS:
            raise ValueError(f"unknown NMII MOTOR connector {name!r}")
        if _REQUIRED_CONNECTORS[name] != component_b:
            raise ValueError(f"connector {name!r} must target {_REQUIRED_CONNECTORS[name]!r}")
        self.name = name
        self.component_a = NMII_COMPONENT
        self.component_b = component_b
        self.event_channels = frozenset(event_channels)
        self.ledger_channels = frozenset(ledger_channels)
        self.state = state
        self.head_node_d = head_node_d
        self.params = params
        self._launch = launch
        self._copy = copy
        self._device = device
        self.launched_kernels: list[object] = []
        self.ledger_calls: list[object] = []

    def accumulate_candidate(self, actuator: NMIIActuatorView, port: FilamentMotorPortView) -> None:
        """Refresh live walk polarity, then scatter the two-array adjoint crossbridge force.

        The head node identity/position/force come from the actuator view (``nmii``-owned); the segment
        endpoints/polarity/position/force come from the target-owned port.  Geometry is refreshed first so the
        walk direction tracks the live actin barbed end, then ``crossbridge_segment_split_kernel`` adds ``+f`` to
        the actuator force and the barycentric reaction ``−(1−t)·f``/``−t·f`` to the port force — Newton-3rd
        across the two never-merged arrays (the bidirectional, adjoint MOTOR load path).
        """
        self.launched_kernels = []
        self._launch(refresh_segment_barbed_kernel, dim=int(port.segment_node_a_d.shape[0]),
                     inputs=[port.position_d, port.segment_node_a_d, port.segment_node_b_d,
                             port.segment_polarity_d, self.state.seg_barbed_d], device=self._device)
        self.launched_kernels.append(refresh_segment_barbed_kernel)
        self._launch(_refresh_bound_walk_dir_kernel, dim=self.state.n_heads,
                     inputs=[port.position_d, self.state.bound_d, self.state.seg_id_d,
                             port.segment_node_a_d, port.segment_node_b_d, port.segment_polarity_d,
                             self.state.walk_dir_d], device=self._device)
        self.launched_kernels.append(_refresh_bound_walk_dir_kernel)
        self._launch(crossbridge_segment_split_kernel, dim=self.state.n_heads,
                     inputs=[actuator.position_d, actuator.force_d, port.position_d, port.force_d,
                             self.head_node_d, self.state.bound_d, self.state.seg_a_d, self.state.seg_b_d,
                             self.state.bary_t_d, self.state.abscissa_d, self.state.walk_dir_d,
                             self.params.k_xb, self.params.r0_xb], device=self._device)
        self.launched_kernels.append(crossbridge_segment_split_kernel)

    def compute_loads(self, actuator: NMIIActuatorView, port: FilamentMotorPortView) -> None:
        """Fill the connector-owned per-head Hill/Bell loads from the converged split-ownership geometry.

        The composed step calls this after the coupled candidate converges and before ``commit_irreversible``
        consumes ``loads_hill_d`` (Hill tangential/walk-dir) and ``loads_bell_d`` (Bell full-|F|).  Head node
        reads the ``nmii``-owned actuator position; the segment nodes read the target-owned port position — the
        two arrays are addressed independently, never merged.
        """
        self.launched_kernels = []
        self._launch(compute_head_loads_segment_split_kernel, dim=self.state.n_heads,
                     inputs=[actuator.position_d, port.position_d, self.head_node_d, self.state.bound_d,
                             self.state.seg_a_d, self.state.seg_b_d, self.state.bary_t_d,
                             self.state.abscissa_d, self.state.walk_dir_d, self.params.k_xb,
                             self.params.r0_xb, self.state.loads_hill_d, self.state.loads_bell_d],
                     device=self._device)
        self.launched_kernels.append(compute_head_loads_segment_split_kernel)

    def snapshot_candidate(self) -> None:
        """D2D-snapshot every connector-owned authoritative binding array."""
        for live, snap in self.state.authoritative_pairs():
            self._copy(snap, live)

    def rollback(self, accepted: wp.array) -> None:
        """Restore the snapshot; the gated KMC never mutates committed state on a rejected step."""
        for live, snap in self.state.authoritative_pairs():
            self._copy(live, snap)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Run the accepted-predicated attach → Hill-step/Bell-detach → epoch KMC through the real kernels."""
        self.launched_kernels = []
        seed_attach = wp.int32(int(rng_seed) & 0x7FFFFFFF)
        seed_detach = wp.int32((int(rng_seed) ^ 0x5BD1E995) & 0x7FFFFFFF)
        self._launch(attach_segment_gated_kernel, dim=self.state.n_heads,
                     inputs=[accepted, self.state.rng_epoch_d, self.state.bound_d, self.state.seg_id_d,
                             self.state.seg_a_d, self.state.seg_b_d, self.state.bary_t_d,
                             self.state.abscissa_d, self.state.walk_dir_d, self.state.query_seg_id_d,
                             self.state.query_t_d, self.state.query_barbed_d, self._seg_node_a(),
                             self._seg_node_b(), self.params, wp.float64(dt_phys), seed_attach],
                     device=self._device)
        self.launched_kernels.append(attach_segment_gated_kernel)
        self._launch(step_detach_segment_gated_kernel, dim=self.state.n_heads,
                     inputs=[accepted, self.state.rng_epoch_d, self.state.bound_d, self.state.seg_id_d,
                             self.state.seg_a_d, self.state.seg_b_d, self.state.bary_t_d,
                             self.state.abscissa_d, self.state.walk_dir_d, self.state.loads_hill_d,
                             self.state.loads_bell_d, self.params, wp.float64(dt_phys), seed_detach],
                     device=self._device)
        self.launched_kernels.append(step_detach_segment_gated_kernel)
        self._launch(_increment_epoch_if_accepted_kernel, dim=1,
                     inputs=[accepted, self.state.rng_epoch_d], device=self._device)
        self.launched_kernels.append(_increment_epoch_if_accepted_kernel)

    def bind_target_topology(self, seg_node_a_d: wp.array, seg_node_b_d: wp.array) -> None:
        """Cache the current target segment endpoint arrays used by the accepted attach kernel."""
        self._target_seg_node_a = seg_node_a_d
        self._target_seg_node_b = seg_node_b_d

    def _seg_node_a(self) -> wp.array:
        node = getattr(self, "_target_seg_node_a", None)
        if node is None:
            raise RuntimeError("bind_target_topology must be called before commit")
        return node

    def _seg_node_b(self) -> wp.array:
        node = getattr(self, "_target_seg_node_b", None)
        if node is None:
            raise RuntimeError("bind_target_topology must be called before commit")
        return node

    def accumulate_ledger(self, ledger: object) -> None:
        """Record the ledger handle; force/work/ATP/population reduction lands with the two-array kernel."""
        self.ledger_calls.append(ledger)


__all__ = [
    "BackboneArmMechanics",
    "CORTEX_COMPONENT",
    "FILOPODIUM_COMPONENT",
    "FilamentMotorPortView",
    "LAMELLIPODIUM_COMPONENT",
    "LegacyNMIIAdapter",
    "LegacyNMIIKind",
    "NMIIActuator",
    "NMIIActuatorStateOwner",
    "NMIIActuatorStepBindings",
    "NMIIActuatorView",
    "NMIIBendingStiffness",
    "NMIIContractProposal",
    "NMIILedgerContributor",
    "NMIIMotorConnector",
    "NMIITransaction",
    "SegmentConnectorState",
    "SegmentMotorConnectorRuntime",
    "NMII_COMPONENT",
    "NMII_CONNECTOR_EVENT_CHANNELS",
    "NMII_CORTEX_MOTOR",
    "NMII_FILOPODIUM_MOTOR",
    "NMII_LAMELLIPODIUM_MOTOR",
    "NMII_LEDGER_CHANNELS",
    "NMII_REPRESENTATION",
    "NMII_SF_MOTOR",
    "NMII_STATE_EVENT_CHANNELS",
    "SF_COMPONENT",
    "nmii_contract_proposal",
    "proposed_nmii_architecture",
]
