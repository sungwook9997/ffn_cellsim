"""Graph-owned LINC joint population shared by actin-cap, MT, and IF edges.

This module is the first structural connector slice.  Filament and nucleus actors retain
their geometry; a LINC population owns only joint identity, endpoint ``PortRef`` fields,
binding state, kinetic transaction state, and force/work outputs.  Nuclear reactions are
applied through one barycentric socket projector, which selects native surface scatter or
the reduced ``J.T`` path without changing the joint representation.

The legacy nearest-cortex fixed spring is not a runtime backend for this connector.

Sanity Gate:
    * all authoritative SoA state is CUDA-resident, non-aliased, and fixed-capacity;
    * host resolution rejects stale actor/entity generations and wrong edge roles;
    * device generation validation runs before mechanics without a host readback;
    * segment and triangle gather/scatter use the same weights, closing force and work;
    * reduced projection uses ``f_q = J.T f_socket`` and preserves virtual work;
    * bind/unbind, age, RNG epoch, and population epoch commit only under the one
      scheduler-owned acceptance predicate; rejection restores the candidate snapshot.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import IntEnum, StrEnum
from typing import Protocol, runtime_checkable

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.engine.contracts import CellArchitecture, ConnectorFamily
from aleph.engine.load_path import ActorRecord, ActorRegistry, ElementKind, PortRef
from aleph.engine.runtime import TransactionParticipant

ACTIN_CAP_LINC = "actin_cap_linc"
MT_NUCLEUS_LINC = "mt_nucleus_linc"
IF_NUCLEUS_LINC = "if_nucleus_linc"
NUCLEUS_COMPONENT = "nucleus"

LINC_EDGES = frozenset({ACTIN_CAP_LINC, MT_NUCLEUS_LINC, IF_NUCLEUS_LINC})
LINC_EVENT_CHANNELS = frozenset({"bind", "unbind"})

_EDGE_COMPONENTS = {
    ACTIN_CAP_LINC: "sf_arc",
    MT_NUCLEUS_LINC: "microtubule",
    IF_NUCLEUS_LINC: "intermediate_filament",
}
_EDGE_CHEMISTRY = {
    ACTIN_CAP_LINC: "nesprin_actin_linc",
    MT_NUCLEUS_LINC: "microtubule_motor_linc",
    IF_NUCLEUS_LINC: "nesprin3_plectin_if_linc",
}


class LincJointState(IntEnum):
    """Committed state of one preallocated LINC joint record."""

    UNBOUND = 0
    BOUND = 1


class LincFilamentRole(StrEnum):
    """Typed filament-side endpoint roles for the three LINC graph edges."""

    PERINUCLEAR_ACTIN_CAP = "perinuclear_actin_cap"
    MT_PLUS_END_OR_LATTICE = "mt_plus_end_or_lattice"
    IF_JUNCTION_OR_MATERIAL = "if_junction_or_material"


class NuclearProjectionMode(StrEnum):
    """Nuclear reaction backend selected without changing socket identity."""

    NATIVE_SURFACE = "native_surface"
    REDUCED_JT = "reduced_jt"


_EDGE_ROLES = {
    ACTIN_CAP_LINC: LincFilamentRole.PERINUCLEAR_ACTIN_CAP,
    MT_NUCLEUS_LINC: LincFilamentRole.MT_PLUS_END_OR_LATTICE,
    IF_NUCLEUS_LINC: LincFilamentRole.IF_JUNCTION_OR_MATERIAL,
}


def _storage_key(array: object) -> tuple[object, ...]:
    """Return device-storage identity without copying authoritative state."""
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
    """Validate CUDA residence, dtype, rank, and selected extents from metadata."""
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


def _validate_same_device(arrays: tuple[object, ...], *, label: str) -> None:
    if len({str(getattr(array, "device", None)) for array in arrays}) != 1:
        raise ValueError(f"all {label} arrays must share one CUDA device")


@dataclass(frozen=True, slots=True)
class LincFilamentPort:
    """Generation-aware filament material point used by one LINC joint."""

    ref: PortRef
    role: LincFilamentRole

    def __post_init__(self) -> None:
        if self.ref.element_kind is not ElementKind.SEGMENT:
            raise ValueError("the first LINC slice requires a segment filament PortRef")
        _ = self.ref.segment_u


@dataclass(frozen=True, slots=True)
class NuclearSurfaceSocket:
    """Generation-aware barycentric material point on one nuclear surface face."""

    ref: PortRef

    def __post_init__(self) -> None:
        if self.ref.component != NUCLEUS_COMPONENT:
            raise ValueError("LINC nuclear socket must belong to nucleus")
        if self.ref.element_kind is not ElementKind.TRIANGLE:
            raise ValueError("LINC nuclear socket must address a surface triangle")
        bary = self.barycentric
        if any(weight < 0.0 or weight > 1.0 for weight in bary):
            raise ValueError("nuclear socket barycentric weights must be in [0, 1]")
        total = math.fsum(bary)
        roundoff = np.finfo(np.float64).eps * max(1.0, abs(total)) * len(bary)
        if abs(total - 1.0) > roundoff:
            raise ValueError("nuclear socket barycentric weights must sum to one")
        if self.ref.local_coordinates[3] != 0.0:
            raise ValueError("the unused fourth nuclear socket coordinate must be zero")

    @property
    def barycentric(self) -> tuple[float, float, float]:
        """Return the three triangle interpolation weights."""
        return (
            float(self.ref.local_coordinates[0]),
            float(self.ref.local_coordinates[1]),
            float(self.ref.local_coordinates[2]),
        )


@dataclass(frozen=True, slots=True)
class LincJointSpec:
    """One typed graph joint before fixed-capacity SoA packing."""

    joint_id: int
    joint_generation: int
    edge_name: str
    filament: LincFilamentPort
    nucleus: NuclearSurfaceSocket
    initial_state: LincJointState = LincJointState.UNBOUND

    def __post_init__(self) -> None:
        if self.joint_id < 0 or self.joint_generation < 0:
            raise ValueError("LINC joint ID and generation must be nonnegative")
        if self.edge_name not in LINC_EDGES:
            raise ValueError(f"unknown LINC graph edge {self.edge_name!r}")


def validate_linc_architecture(architecture: CellArchitecture) -> None:
    """Require all three declared kinetic, bidirectional LINC edges."""
    by_name = {connector.name: connector for connector in architecture.connectors}
    missing = LINC_EDGES - set(by_name)
    if missing:
        raise ValueError(f"architecture is missing LINC edges {sorted(missing)}")
    for name in sorted(LINC_EDGES):
        connector = by_name[name]
        expected = {_EDGE_COMPONENTS[name], NUCLEUS_COMPONENT}
        if {connector.component_a, connector.component_b} != expected:
            raise ValueError(f"{name} has incorrect component endpoints")
        if connector.family is not ConnectorFamily.LINC:
            raise ValueError(f"{name} must use ConnectorFamily.LINC")
        if not connector.kinetics or not connector.commit_on_accept:
            raise ValueError(f"{name} must have accepted-step kinetics")
        if not connector.bidirectional or not connector.adjoint_transfer_required:
            raise ValueError(f"{name} must be bidirectional and require adjoint transfer")
        if connector.chemistry_card != _EDGE_CHEMISTRY[name]:
            raise ValueError(f"{name} has the wrong chemistry card")


def resolve_linc_joint(
    spec: LincJointSpec,
    registry: ActorRegistry,
    architecture: CellArchitecture,
) -> LincJointSpec:
    """Validate graph edge, filament role, identity, generation, and socket geometry."""
    validate_linc_architecture(architecture)
    registry.validate_port(spec.filament.ref)
    registry.validate_port(spec.nucleus.ref)
    expected_component = _EDGE_COMPONENTS[spec.edge_name]
    if spec.filament.ref.component != expected_component:
        raise ValueError(f"{spec.edge_name} requires a {expected_component} filament port")
    if spec.filament.role is not _EDGE_ROLES[spec.edge_name]:
        raise ValueError(f"{spec.edge_name} has the wrong filament endpoint role")
    if spec.nucleus.ref.component != NUCLEUS_COMPONENT:
        raise ValueError("LINC socket must belong to nucleus")
    return spec


@dataclass(frozen=True, slots=True)
class LincPopulationView:
    """Non-owning CUDA SoA view consumed by mechanics, projection, and ledgers."""

    edge_name: str
    capacity: int
    active_d: wp.array
    joint_id_d: wp.array
    joint_generation_d: wp.array
    state_d: wp.array
    candidate_state_d: wp.array
    age_d: wp.array
    rng_epoch_d: wp.array
    filament_actor_id_d: wp.array
    filament_actor_generation_d: wp.array
    filament_entity_id_d: wp.array
    filament_entity_generation_d: wp.array
    filament_element_id_d: wp.array
    filament_u_d: wp.array
    nucleus_actor_id_d: wp.array
    nucleus_actor_generation_d: wp.array
    nucleus_entity_id_d: wp.array
    nucleus_entity_generation_d: wp.array
    socket_face_id_d: wp.array
    socket_barycentric_d: wp.array
    rest_length_d: wp.array
    stiffness_d: wp.array
    stiffening_d: wp.array
    load_d: wp.array
    energy_d: wp.array
    force_on_filament_d: wp.array
    force_on_nucleus_d: wp.array
    population_epoch_d: wp.array


@dataclass(frozen=True, slots=True)
class FilamentGeometryView:
    """Non-owning live filament geometry and device generations."""

    record: ActorRecord
    position_d: wp.array
    force_d: wp.array
    segments_d: wp.array
    actor_generation_d: wp.array
    entity_generation_d: wp.array

    def __post_init__(self) -> None:
        _validate_device_array(self.position_d, label="filament.position_d", dtype=wp.vec3d)
        _validate_device_array(
            self.force_d,
            label="filament.force_d",
            dtype=wp.vec3d,
            shape0=int(self.position_d.shape[0]),
        )
        _validate_device_array(
            self.segments_d,
            label="filament.segments_d",
            dtype=wp.int32,
            rank=2,
            shape0=self.record.n_elements,
            shape1=2,
        )
        _validate_device_array(
            self.actor_generation_d,
            label="filament.actor_generation_d",
            dtype=wp.int32,
            shape0=1,
        )
        _validate_device_array(
            self.entity_generation_d,
            label="filament.entity_generation_d",
            dtype=wp.int32,
            shape0=1,
        )
        arrays = (
            self.position_d,
            self.force_d,
            self.segments_d,
            self.actor_generation_d,
            self.entity_generation_d,
        )
        _validate_same_device(arrays, label="filament geometry")
        if _storage_key(self.position_d) == _storage_key(self.force_d):
            raise ValueError("filament position and force arrays must not alias")


@dataclass(frozen=True, slots=True)
class NuclearSurfaceView:
    """Native surface geometry plus the optional reduced generalized-force target."""

    record: ActorRecord
    mode: NuclearProjectionMode
    surface_position_d: wp.array
    surface_force_d: wp.array
    faces_d: wp.array
    actor_generation_d: wp.array
    entity_generation_d: wp.array
    generalized_force_d: wp.array | None = None

    def __post_init__(self) -> None:
        if self.record.component != NUCLEUS_COMPONENT:
            raise ValueError("nuclear surface view must belong to nucleus")
        _validate_device_array(
            self.surface_position_d,
            label="nucleus.surface_position_d",
            dtype=wp.vec3d,
        )
        _validate_device_array(
            self.surface_force_d,
            label="nucleus.surface_force_d",
            dtype=wp.vec3d,
            shape0=int(self.surface_position_d.shape[0]),
        )
        _validate_device_array(
            self.faces_d,
            label="nucleus.faces_d",
            dtype=wp.int32,
            rank=2,
            shape0=self.record.n_elements,
            shape1=3,
        )
        for label, array in (
            ("actor_generation_d", self.actor_generation_d),
            ("entity_generation_d", self.entity_generation_d),
        ):
            _validate_device_array(
                array,
                label=f"nucleus.{label}",
                dtype=wp.int32,
                shape0=1,
            )
        arrays: tuple[object, ...] = (
            self.surface_position_d,
            self.surface_force_d,
            self.faces_d,
            self.actor_generation_d,
            self.entity_generation_d,
        )
        if self.mode is NuclearProjectionMode.REDUCED_JT:
            if self.generalized_force_d is None:
                raise ValueError("reduced LINC projection requires generalized_force_d")
            _validate_device_array(
                self.generalized_force_d,
                label="nucleus.generalized_force_d",
                dtype=wp.float64,
            )
            arrays = (*arrays, self.generalized_force_d)
        elif self.generalized_force_d is not None:
            raise ValueError("native LINC projection must not supply generalized_force_d")
        _validate_same_device(arrays, label="nuclear surface")
        if _storage_key(self.surface_position_d) == _storage_key(self.surface_force_d):
            raise ValueError("nuclear surface position and force arrays must not alias")


@runtime_checkable
class LincMechanics(Protocol):
    """Injected Warp mechanics writing pair loads and filament/socket reactions."""

    def accumulate(
        self,
        joints: LincPopulationView,
        filament: FilamentGeometryView,
        nucleus: NuclearSurfaceView,
    ) -> None:
        """Gather both endpoints and scatter the filament reaction adjointly."""


@runtime_checkable
class LincPortGenerationValidator(Protocol):
    """Device validation of packed PortRef generations and element bounds."""

    def validate(
        self,
        joints: LincPopulationView,
        filament: FilamentGeometryView,
        nucleus: NuclearSurfaceView,
    ) -> None:
        """Set device failure state for a stale port; never read it on the host here."""


@runtime_checkable
class NuclearSocketProjector(Protocol):
    """Apply socket reactions through native barycentric scatter or reduced ``J.T``."""

    def accumulate_reaction(
        self,
        joints: LincPopulationView,
        nucleus: NuclearSurfaceView,
    ) -> None:
        """Scatter ``force_on_nucleus_d`` using the mode declared by ``nucleus``."""


@runtime_checkable
class LincKinetics(Protocol):
    """Device candidate generator for force-dependent bind/unbind transitions."""

    event_channels: frozenset[str]

    def propose_candidate(
        self,
        joints: LincPopulationView,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        """Write candidate state only; do not advance accepted RNG or physical time."""


@runtime_checkable
class LincTransaction(TransactionParticipant, Protocol):
    """Accepted-step transaction covering state, topology/remap, and epoch arrays."""

    event_channels: frozenset[str]

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return every authoritative array that creation, removal, or remap may mutate."""


@runtime_checkable
class LincLedgerContributor(Protocol):
    """Force/work and active/bound/unbound population ledger contribution."""

    def accumulate_linc_ledger(self, joints: LincPopulationView, ledger: object) -> None:
        """Accumulate without deciding global acceptance."""


@dataclass(frozen=True, slots=True)
class LincJointPopulation:
    """Fixed-capacity CUDA SoA owned by one of the three LINC graph edges."""

    edge_name: str
    capacity: int
    active_d: wp.array
    joint_id_d: wp.array
    joint_generation_d: wp.array
    state_d: wp.array
    candidate_state_d: wp.array
    snapshot_state_d: wp.array
    age_d: wp.array
    snapshot_age_d: wp.array
    rng_epoch_d: wp.array
    snapshot_rng_epoch_d: wp.array
    filament_actor_id_d: wp.array
    filament_actor_generation_d: wp.array
    filament_entity_id_d: wp.array
    filament_entity_generation_d: wp.array
    filament_element_id_d: wp.array
    filament_u_d: wp.array
    nucleus_actor_id_d: wp.array
    nucleus_actor_generation_d: wp.array
    nucleus_entity_id_d: wp.array
    nucleus_entity_generation_d: wp.array
    socket_face_id_d: wp.array
    socket_barycentric_d: wp.array
    rest_length_d: wp.array
    stiffness_d: wp.array
    stiffening_d: wp.array
    load_d: wp.array
    energy_d: wp.array
    force_on_filament_d: wp.array
    force_on_nucleus_d: wp.array
    population_epoch_d: wp.array
    snapshot_population_epoch_d: wp.array
    mechanics: LincMechanics
    generation_validator: LincPortGenerationValidator
    projector: NuclearSocketProjector
    kinetics: LincKinetics
    transaction: LincTransaction
    ledger: LincLedgerContributor

    def __post_init__(self) -> None:
        if self.edge_name not in LINC_EDGES:
            raise ValueError(f"unknown LINC population edge {self.edge_name!r}")
        if self.capacity <= 0:
            raise ValueError("LINC population capacity must be positive")
        typed_arrays = (
            ("active_d", self.active_d, wp.int32),
            ("joint_id_d", self.joint_id_d, wp.int64),
            ("joint_generation_d", self.joint_generation_d, wp.int32),
            ("state_d", self.state_d, wp.int32),
            ("candidate_state_d", self.candidate_state_d, wp.int32),
            ("snapshot_state_d", self.snapshot_state_d, wp.int32),
            ("age_d", self.age_d, wp.float64),
            ("snapshot_age_d", self.snapshot_age_d, wp.float64),
            ("rng_epoch_d", self.rng_epoch_d, wp.int64),
            ("snapshot_rng_epoch_d", self.snapshot_rng_epoch_d, wp.int64),
            ("filament_actor_id_d", self.filament_actor_id_d, wp.int64),
            ("filament_actor_generation_d", self.filament_actor_generation_d, wp.int32),
            ("filament_entity_id_d", self.filament_entity_id_d, wp.int64),
            ("filament_entity_generation_d", self.filament_entity_generation_d, wp.int32),
            ("filament_element_id_d", self.filament_element_id_d, wp.int32),
            ("filament_u_d", self.filament_u_d, wp.float64),
            ("nucleus_actor_id_d", self.nucleus_actor_id_d, wp.int64),
            ("nucleus_actor_generation_d", self.nucleus_actor_generation_d, wp.int32),
            ("nucleus_entity_id_d", self.nucleus_entity_id_d, wp.int64),
            ("nucleus_entity_generation_d", self.nucleus_entity_generation_d, wp.int32),
            ("socket_face_id_d", self.socket_face_id_d, wp.int32),
            ("socket_barycentric_d", self.socket_barycentric_d, wp.vec3d),
            ("rest_length_d", self.rest_length_d, wp.float64),
            ("stiffness_d", self.stiffness_d, wp.float64),
            ("stiffening_d", self.stiffening_d, wp.float64),
            ("load_d", self.load_d, wp.float64),
            ("energy_d", self.energy_d, wp.float64),
            ("force_on_filament_d", self.force_on_filament_d, wp.vec3d),
            ("force_on_nucleus_d", self.force_on_nucleus_d, wp.vec3d),
        )
        for label, array, dtype in typed_arrays:
            _validate_device_array(
                array,
                label=f"{self.edge_name}.{label}",
                dtype=dtype,
                shape0=self.capacity,
            )
        for label, array in (
            ("population_epoch_d", self.population_epoch_d),
            ("snapshot_population_epoch_d", self.snapshot_population_epoch_d),
        ):
            _validate_device_array(
                array,
                label=f"{self.edge_name}.{label}",
                dtype=wp.int64,
                shape0=1,
            )
        arrays = tuple(array for _, array, _ in typed_arrays) + (
            self.population_epoch_d,
            self.snapshot_population_epoch_d,
        )
        _validate_same_device(arrays, label=f"{self.edge_name} population")
        if len({_storage_key(array) for array in arrays}) != len(arrays):
            raise ValueError("LINC population arrays must own distinct storage")
        missing_kinetics = LINC_EVENT_CHANNELS - self.kinetics.event_channels
        if missing_kinetics:
            raise ValueError(f"LINC kinetics is missing event channels {sorted(missing_kinetics)}")
        missing_transaction = LINC_EVENT_CHANNELS - self.transaction.event_channels
        if missing_transaction:
            raise ValueError(
                f"LINC transaction is missing event channels {sorted(missing_transaction)}"
            )
        mutable = (
            self.active_d,
            self.joint_id_d,
            self.joint_generation_d,
            self.state_d,
            self.candidate_state_d,
            self.snapshot_state_d,
            self.age_d,
            self.snapshot_age_d,
            self.rng_epoch_d,
            self.snapshot_rng_epoch_d,
            self.filament_actor_id_d,
            self.filament_actor_generation_d,
            self.filament_entity_id_d,
            self.filament_entity_generation_d,
            self.filament_element_id_d,
            self.filament_u_d,
            self.nucleus_actor_id_d,
            self.nucleus_actor_generation_d,
            self.nucleus_entity_id_d,
            self.nucleus_entity_generation_d,
            self.socket_face_id_d,
            self.socket_barycentric_d,
            self.rest_length_d,
            self.stiffness_d,
            self.stiffening_d,
            self.population_epoch_d,
            self.snapshot_population_epoch_d,
        )
        covered = {_storage_key(array) for array in self.transaction.owned_arrays()}
        if {_storage_key(array) for array in mutable} - covered:
            raise ValueError("LINC transaction does not cover all mutable/remap population arrays")

    def view(self) -> LincPopulationView:
        """Return the non-owning arrays needed by mechanics, projection, and ledgers."""
        return LincPopulationView(
            edge_name=self.edge_name,
            capacity=self.capacity,
            active_d=self.active_d,
            joint_id_d=self.joint_id_d,
            joint_generation_d=self.joint_generation_d,
            state_d=self.state_d,
            candidate_state_d=self.candidate_state_d,
            age_d=self.age_d,
            rng_epoch_d=self.rng_epoch_d,
            filament_actor_id_d=self.filament_actor_id_d,
            filament_actor_generation_d=self.filament_actor_generation_d,
            filament_entity_id_d=self.filament_entity_id_d,
            filament_entity_generation_d=self.filament_entity_generation_d,
            filament_element_id_d=self.filament_element_id_d,
            filament_u_d=self.filament_u_d,
            nucleus_actor_id_d=self.nucleus_actor_id_d,
            nucleus_actor_generation_d=self.nucleus_actor_generation_d,
            nucleus_entity_id_d=self.nucleus_entity_id_d,
            nucleus_entity_generation_d=self.nucleus_entity_generation_d,
            socket_face_id_d=self.socket_face_id_d,
            socket_barycentric_d=self.socket_barycentric_d,
            rest_length_d=self.rest_length_d,
            stiffness_d=self.stiffness_d,
            stiffening_d=self.stiffening_d,
            load_d=self.load_d,
            energy_d=self.energy_d,
            force_on_filament_d=self.force_on_filament_d,
            force_on_nucleus_d=self.force_on_nucleus_d,
            population_epoch_d=self.population_epoch_d,
        )

    def accumulate(self, filament: FilamentGeometryView, nucleus: NuclearSurfaceView) -> None:
        """Validate generations, assemble pair force, then project the nuclear reaction."""
        expected_component = _EDGE_COMPONENTS[self.edge_name]
        if filament.record.component != expected_component:
            raise ValueError(f"{self.edge_name} requires filament component {expected_component}")
        view = self.view()
        self.generation_validator.validate(view, filament, nucleus)
        self.mechanics.accumulate(view, filament, nucleus)
        self.projector.accumulate_reaction(view, nucleus)

    def snapshot_candidate(self) -> None:
        """Forward the candidate snapshot to the transaction owner."""
        self.transaction.snapshot_candidate()

    def propose_candidate(self, dt_phys: float, rng_seed: int) -> None:
        """Generate bind/unbind candidates without committing physical state."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if rng_seed < 0:
            raise ValueError("rng_seed must be nonnegative")
        self.kinetics.propose_candidate(self.view(), dt_phys, rng_seed)

    def rollback(self, accepted: wp.array) -> None:
        """Restore candidate/snapshot state under the scheduler-owned predicate."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit bind/unbind, age, RNG, and population epoch only on acceptance."""
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute force/work and active/bound/unbound population counts."""
        self.ledger.accumulate_linc_ledger(self.view(), ledger)


@dataclass(frozen=True, slots=True)
class LincRigConnectorAdapter:
    """Expose one LINC population through the MT/IF graph-connector interface.

    The rig facades deliberately pass compact public views that do not contain the segment and
    nuclear-face generation arrays required by LINC.  This adapter owns no state: it binds those
    richer, generation-aware views once at build time, verifies that the facade supplied the same
    position/force storage, and forwards mechanics/transaction/ledger calls to one population.
    """

    population: LincJointPopulation
    filament: FilamentGeometryView
    nucleus: NuclearSurfaceView

    def __post_init__(self) -> None:
        expected = _EDGE_COMPONENTS[self.population.edge_name]
        if self.filament.record.component != expected:
            raise ValueError(
                f"{self.population.edge_name} adapter requires filament component {expected}"
            )
        if self.nucleus.record.component != NUCLEUS_COMPONENT:
            raise ValueError("LINC rig adapter requires a nucleus surface view")

    @property
    def name(self) -> str:
        """Graph connector name consumed by MT/IF binding validation."""
        return self.population.edge_name

    @property
    def component_a(self) -> str:
        """Filament-side registered component."""
        return _EDGE_COMPONENTS[self.population.edge_name]

    @property
    def component_b(self) -> str:
        """Nuclear-side registered component."""
        return NUCLEUS_COMPONENT

    def accumulate_rig(self, rig: object, endpoint: object) -> None:
        """Verify facade storage identity and accumulate the generation-aware LINC pair load."""
        rig_position = getattr(rig, "position_d", None)
        rig_force = getattr(rig, "force_d", None)
        endpoint_position = getattr(endpoint, "position_d", None)
        endpoint_force = getattr(endpoint, "force_d", None)
        if (
            _storage_key(rig_position) != _storage_key(self.filament.position_d)
            or _storage_key(rig_force) != _storage_key(self.filament.force_d)
        ):
            raise ValueError(f"{self.name} rig view does not match its bound filament geometry")
        if (
            _storage_key(endpoint_position) != _storage_key(self.nucleus.surface_position_d)
            or _storage_key(endpoint_force) != _storage_key(self.nucleus.surface_force_d)
        ):
            raise ValueError(f"{self.name} endpoint does not match its bound nuclear surface")
        self.population.accumulate(self.filament, self.nucleus)

    def snapshot_candidate(self) -> None:
        """Forward candidate snapshot to the single population owner."""
        self.population.snapshot_candidate()

    def propose_candidate(self, dt_phys: float, rng_seed: int) -> None:
        """Forward reversible bind/unbind proposal generation."""
        self.population.propose_candidate(dt_phys, rng_seed)

    def rollback(self, accepted: wp.array) -> None:
        """Forward rejected-candidate restoration."""
        self.population.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Forward accepted bind/unbind commit."""
        self.population.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward force/work/population accounting exactly once."""
        self.population.accumulate_ledger(ledger)


@dataclass(frozen=True, slots=True)
class LincPairReference:
    """Native segment↔triangle equal-opposite and virtual-work oracle."""

    filament_force: npt.NDArray[np.float64]
    nucleus_force: npt.NDArray[np.float64]
    force_on_filament: npt.NDArray[np.float64]
    force_residual_pn: npt.NDArray[np.float64]
    energy_pn_um: float
    virtual_work_pn_um: float | None


def segment_triangle_linc_reference(
    filament_position: npt.ArrayLike,
    segment: tuple[int, int],
    u: float,
    nucleus_position: npt.ArrayLike,
    face: tuple[int, int, int],
    barycentric: tuple[float, float, float],
    *,
    rest_length_um: float,
    stiffness_pn_per_um: float,
    stiffening_per_um2: float,
    filament_displacement: npt.ArrayLike | None = None,
    nucleus_displacement: npt.ArrayLike | None = None,
) -> LincPairReference:
    """Evaluate the legacy tension-only law with adjoint segment/triangle scatter.

    The law is a reference oracle only. Its parameters have no defaults because current
    LINC stiffness/kinetic magnitudes remain source/PI gaps.
    """
    fp = np.asarray(filament_position, dtype=np.float64)
    np_ = np.asarray(nucleus_position, dtype=np.float64)
    if fp.ndim != 2 or np_.ndim != 2 or fp.shape[1:] != (3,) or np_.shape[1:] != (3,):
        raise ValueError("positions must be (N, 3)")
    if not 0.0 <= u <= 1.0:
        raise ValueError("segment coordinate must be in [0, 1]")
    if rest_length_um < 0.0 or not math.isfinite(rest_length_um):
        raise ValueError("rest length must be finite and nonnegative")
    if stiffness_pn_per_um <= 0.0 or not math.isfinite(stiffness_pn_per_um):
        raise ValueError("stiffness must be finite and positive")
    if stiffening_per_um2 < 0.0 or not math.isfinite(stiffening_per_um2):
        raise ValueError("stiffening must be finite and nonnegative")
    weights = np.asarray(barycentric, dtype=np.float64)
    if weights.shape != (3,) or np.any(weights < 0.0) or not np.isclose(weights.sum(), 1.0):
        raise ValueError("triangle barycentric weights must be nonnegative and sum to one")
    i0, i1 = segment
    j0, j1, j2 = face
    xf = (1.0 - u) * fp[i0] + u * fp[i1]
    xn = weights[0] * np_[j0] + weights[1] * np_[j1] + weights[2] * np_[j2]
    delta = xn - xf
    length = float(np.linalg.norm(delta))
    extension = max(length - rest_length_um, 0.0)
    pair = np.zeros(3, dtype=np.float64)
    if length > 0.0 and extension > 0.0:
        magnitude = stiffness_pn_per_um * extension * (
            1.0 + stiffening_per_um2 * extension**2
        )
        pair = magnitude * delta / length
    ff = np.zeros_like(fp)
    nf = np.zeros_like(np_)
    ff[i0] += (1.0 - u) * pair
    ff[i1] += u * pair
    nf[j0] -= weights[0] * pair
    nf[j1] -= weights[1] * pair
    nf[j2] -= weights[2] * pair
    work = None
    if filament_displacement is not None or nucleus_displacement is not None:
        if filament_displacement is None or nucleus_displacement is None:
            raise ValueError("both displacement arrays are required for virtual work")
        fd = np.asarray(filament_displacement, dtype=np.float64)
        nd = np.asarray(nucleus_displacement, dtype=np.float64)
        if fd.shape != fp.shape or nd.shape != np_.shape:
            raise ValueError("displacements must match position arrays")
        work = float(np.sum(ff * fd) + np.sum(nf * nd))
    energy = stiffness_pn_per_um * (
        0.5 * extension**2 + 0.25 * stiffening_per_um2 * extension**4
    )
    return LincPairReference(
        filament_force=ff,
        nucleus_force=nf,
        force_on_filament=pair,
        force_residual_pn=ff.sum(axis=0) + nf.sum(axis=0),
        energy_pn_um=float(energy),
        virtual_work_pn_um=work,
    )


@dataclass(frozen=True, slots=True)
class ReducedProjectionReference:
    """Reference result for one reduced nuclear socket ``J.T`` projection."""

    generalized_force: npt.NDArray[np.float64]
    socket_displacement: npt.NDArray[np.float64]
    work_residual_pn_um: float


def reduced_jt_projection_reference(
    socket_jacobian: npt.ArrayLike,
    socket_force: npt.ArrayLike,
    generalized_displacement: npt.ArrayLike,
) -> ReducedProjectionReference:
    """Apply ``J.T`` and report native-socket versus reduced virtual-work closure."""
    jacobian = np.asarray(socket_jacobian, dtype=np.float64)
    force = np.asarray(socket_force, dtype=np.float64)
    dq = np.asarray(generalized_displacement, dtype=np.float64)
    if jacobian.ndim != 2 or jacobian.shape[0] != 3:
        raise ValueError("socket Jacobian must have shape (3, n_generalized)")
    if force.shape != (3,) or dq.shape != (jacobian.shape[1],):
        raise ValueError("force/displacement dimensions do not match the socket Jacobian")
    generalized_force = jacobian.T @ force
    socket_displacement = jacobian @ dq
    residual = float(generalized_force @ dq - force @ socket_displacement)
    return ReducedProjectionReference(generalized_force, socket_displacement, residual)


@dataclass(frozen=True, slots=True)
class LincLedgerReport:
    """Host closeout record for force/work and population conservation gates."""

    capacity: int
    active: int
    bound: int
    unbound: int
    force_residual_pn: float
    work_residual_pn_um: float

    def assert_closed(
        self,
        *,
        force_tolerance_pn: float,
        work_tolerance_pn_um: float,
    ) -> None:
        """Reject invalid counts or an equal-opposite/adjoint closure failure."""
        if self.capacity <= 0:
            raise ValueError("LINC capacity must be positive")
        if min(self.active, self.bound, self.unbound) < 0:
            raise ValueError("LINC population counts must be nonnegative")
        if self.active > self.capacity or self.bound + self.unbound != self.active:
            raise ValueError("LINC bound/unbound counts do not close to active population")
        if force_tolerance_pn < 0.0 or work_tolerance_pn_um < 0.0:
            raise ValueError("ledger tolerances must be nonnegative")
        if abs(self.force_residual_pn) > force_tolerance_pn:
            raise ValueError("LINC equal-opposite force ledger does not close")
        if abs(self.work_residual_pn_um) > work_tolerance_pn_um:
            raise ValueError("LINC adjoint work ledger does not close")


# ---------------------------------------------------------------------------
# KERNEL_BOUND slice: real Warp kernels + concrete delegates + CUDA SoA factory
# ---------------------------------------------------------------------------
#
# The SEAMED slice above proves ownership/transaction coverage through injected
# protocol delegates.  The block below advances the *force-transferring* path
# (generation gate -> segment/triangle pair mechanics -> native socket scatter ->
# accepted-step state machine) from SEAMED to KERNEL_BOUND by binding real
# ``@wp.kernel`` functions through the exact same protocol seams, mirroring the
# already-KERNEL_BOUND ``load_path.py`` runtime.  It ports the in-module
# tension-only oracle (:func:`segment_triangle_linc_reference`) into one Warp SoA
# kernel; it introduces no new constitutive law and no sourced constant.  Every
# stiffness/rest/stiffening magnitude is a caller-supplied explicit input, exactly
# as the oracle requires (no defaults; LINC magnitudes remain source/PI gaps).
#
# Deliberately still SEAMED (not bound here) and reported as gaps:
#   * force-dependent bind/unbind *candidate generation* — no sourced on/off-rate
#     cards exist, so :class:`SourceGatedLincKinetics` refuses to fabricate a
#     hazard law rather than invent rate constants (LC-3);
#   * reduced ``J.T`` socket projection — the Core exposes no per-socket Jacobian
#     object yet (plan gap 2), so only the native barycentric backend is bound;
#   * global-ledger merge — the world has no shared LINC reduction target yet
#     (plan gap 5); :class:`WarpLincLedgerReducer` reduces into its own device
#     accumulators and reports at closeout, never inside the physical loop.

_LINC_STATE_BOUND = wp.constant(int(LincJointState.BOUND))
_LINC_LENGTH_FLOOR = wp.constant(wp.float64(1.0e-15))


def _linc_launch(kernel: object, dim: int, inputs: list, outputs: list | None, device: str) -> None:
    """Single indirection for every LINC device launch.

    Concrete delegates route their ``wp.launch`` through this helper so the seam
    binding (kernel identity + authoritative SoA argument order + device) is
    observable to a CPU-green structural test without a CUDA device, while a real
    CUDA run dispatches the identical kernel.
    """
    if outputs is None:
        wp.launch(kernel, dim=dim, inputs=inputs, device=device)
    else:
        wp.launch(kernel, dim=dim, inputs=inputs, outputs=outputs, device=device)


@wp.kernel
def _linc_pair_force_kernel(
    filament_position: wp.array(dtype=wp.vec3d),
    filament_force: wp.array(dtype=wp.vec3d),
    filament_segments: wp.array(dtype=wp.int32, ndim=2),
    nucleus_position: wp.array(dtype=wp.vec3d),
    nucleus_faces: wp.array(dtype=wp.int32, ndim=2),
    active: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    filament_element: wp.array(dtype=wp.int32),
    filament_u: wp.array(dtype=wp.float64),
    socket_face: wp.array(dtype=wp.int32),
    socket_barycentric: wp.array(dtype=wp.vec3d),
    rest_length: wp.array(dtype=wp.float64),
    stiffness: wp.array(dtype=wp.float64),
    stiffening: wp.array(dtype=wp.float64),
    load: wp.array(dtype=wp.float64),
    energy: wp.array(dtype=wp.float64),
    force_on_filament: wp.array(dtype=wp.vec3d),
    force_on_nucleus: wp.array(dtype=wp.vec3d),
) -> None:
    """One thread per joint: gather segment ``(1-u,u)`` and triangle barycentrics,
    apply the tension-only nonlinear law, scatter the filament reaction adjointly,
    and store the equal-opposite socket reaction for the projector."""
    t = wp.tid()
    if active[t] == 0 or state[t] != _LINC_STATE_BOUND:
        load[t] = wp.float64(0.0)
        energy[t] = wp.float64(0.0)
        force_on_filament[t] = wp.vec3d(0.0, 0.0, 0.0)
        force_on_nucleus[t] = wp.vec3d(0.0, 0.0, 0.0)
        return
    e = filament_element[t]
    i0 = filament_segments[e, 0]
    i1 = filament_segments[e, 1]
    w1 = filament_u[t]
    w0 = wp.float64(1.0) - w1
    xf = w0 * filament_position[i0] + w1 * filament_position[i1]
    f = socket_face[t]
    j0 = nucleus_faces[f, 0]
    j1 = nucleus_faces[f, 1]
    j2 = nucleus_faces[f, 2]
    b = socket_barycentric[t]
    xn = b[0] * nucleus_position[j0] + b[1] * nucleus_position[j1] + b[2] * nucleus_position[j2]
    delta = xn - xf
    length = wp.length(delta)
    pair = wp.vec3d(0.0, 0.0, 0.0)
    magnitude = wp.float64(0.0)
    ext = length - rest_length[t]
    if length > _LINC_LENGTH_FLOOR and ext > wp.float64(0.0):
        magnitude = stiffness[t] * ext * (wp.float64(1.0) + stiffening[t] * ext * ext)
        pair = (magnitude / length) * delta
    # Filament point is pulled toward the nuclear socket; scatter is the transpose
    # of the (1-u,u) gather so nodal virtual work equals pair virtual work.
    wp.atomic_add(filament_force, i0, w0 * pair)
    wp.atomic_add(filament_force, i1, w1 * pair)
    force_on_filament[t] = pair
    force_on_nucleus[t] = -pair
    load[t] = magnitude
    energy[t] = stiffness[t] * (
        wp.float64(0.5) * ext * ext + wp.float64(0.25) * stiffening[t] * ext * ext * ext * ext
    )


@wp.kernel
def _linc_native_socket_scatter_kernel(
    active: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    socket_face: wp.array(dtype=wp.int32),
    socket_barycentric: wp.array(dtype=wp.vec3d),
    force_on_nucleus: wp.array(dtype=wp.vec3d),
    nucleus_faces: wp.array(dtype=wp.int32, ndim=2),
    nucleus_force: wp.array(dtype=wp.vec3d),
) -> None:
    """Scatter the per-joint socket reaction to the three triangle vertices using
    the same barycentric weights that located the socket (native backend)."""
    t = wp.tid()
    if active[t] == 0 or state[t] != _LINC_STATE_BOUND:
        return
    f = socket_face[t]
    j0 = nucleus_faces[f, 0]
    j1 = nucleus_faces[f, 1]
    j2 = nucleus_faces[f, 2]
    b = socket_barycentric[t]
    fs = force_on_nucleus[t]
    wp.atomic_add(nucleus_force, j0, b[0] * fs)
    wp.atomic_add(nucleus_force, j1, b[1] * fs)
    wp.atomic_add(nucleus_force, j2, b[2] * fs)


@wp.kernel
def _linc_generation_validate_kernel(
    active: wp.array(dtype=wp.int32),
    filament_actor_generation: wp.array(dtype=wp.int32),
    filament_entity_generation: wp.array(dtype=wp.int32),
    filament_element: wp.array(dtype=wp.int32),
    nucleus_actor_generation: wp.array(dtype=wp.int32),
    nucleus_entity_generation: wp.array(dtype=wp.int32),
    socket_face: wp.array(dtype=wp.int32),
    live_filament_actor_generation: wp.array(dtype=wp.int32),
    live_filament_entity_generation: wp.array(dtype=wp.int32),
    live_nucleus_actor_generation: wp.array(dtype=wp.int32),
    live_nucleus_entity_generation: wp.array(dtype=wp.int32),
    n_filament_segments: wp.int32,
    n_nucleus_faces: wp.int32,
    stale_count: wp.array(dtype=wp.int32),
) -> None:
    """Compare each packed joint's generations/element bounds against the actors'
    live one-element generation arrays; latch stale joints without any host read."""
    t = wp.tid()
    if active[t] == 0:
        return
    stale = int(0)
    if filament_actor_generation[t] != live_filament_actor_generation[0]:
        stale = 1
    if filament_entity_generation[t] != live_filament_entity_generation[0]:
        stale = 1
    if nucleus_actor_generation[t] != live_nucleus_actor_generation[0]:
        stale = 1
    if nucleus_entity_generation[t] != live_nucleus_entity_generation[0]:
        stale = 1
    if filament_element[t] < 0 or filament_element[t] >= n_filament_segments:
        stale = 1
    if socket_face[t] < 0 or socket_face[t] >= n_nucleus_faces:
        stale = 1
    if stale == 1:
        wp.atomic_add(stale_count, 0, 1)


@wp.kernel
def _linc_rollback_kernel(
    accepted: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    candidate_state: wp.array(dtype=wp.int32),
) -> None:
    """Reset scratch candidates after rejection; committed state is never mutated."""
    t = wp.tid()
    if accepted[0] == 0:
        candidate_state[t] = state[t]


@wp.kernel
def _linc_commit_kernel(
    active: wp.array(dtype=wp.int32),
    accepted: wp.array(dtype=wp.int32),
    candidate_state: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    age: wp.array(dtype=wp.float64),
    rng_epoch: wp.array(dtype=wp.int64),
    population_epoch: wp.array(dtype=wp.int64),
    dt_phys: wp.float64,
) -> None:
    """Commit bind/unbind, age, RNG epoch, and population epoch only when the
    scheduler-owned acceptance scalar is non-zero; restore scratch on rejection."""
    t = wp.tid()
    if active[t] != 0 and accepted[0] != 0:
        previous = state[t]
        state[t] = candidate_state[t]
        if candidate_state[t] == previous:
            age[t] = age[t] + dt_phys
        else:
            age[t] = wp.float64(0.0)
        rng_epoch[t] = rng_epoch[t] + wp.int64(1)
    elif accepted[0] == 0:
        candidate_state[t] = state[t]
    if t == 0 and accepted[0] != 0:
        population_epoch[0] = population_epoch[0] + wp.int64(1)


@wp.kernel
def _linc_ledger_reduce_kernel(
    active: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    energy: wp.array(dtype=wp.float64),
    force_on_filament: wp.array(dtype=wp.vec3d),
    force_on_nucleus: wp.array(dtype=wp.vec3d),
    active_count: wp.array(dtype=wp.int32),
    bound_count: wp.array(dtype=wp.int32),
    energy_sum: wp.array(dtype=wp.float64),
    force_residual: wp.array(dtype=wp.vec3d),
) -> None:
    """Per-edge device reduction of population counts, elastic energy, and the
    equal-opposite force residual (filament + socket reaction, exactly zero per
    bound joint by construction)."""
    t = wp.tid()
    if active[t] == 0:
        return
    wp.atomic_add(active_count, 0, 1)
    if state[t] == _LINC_STATE_BOUND:
        wp.atomic_add(bound_count, 0, 1)
        wp.atomic_add(energy_sum, 0, energy[t])
        wp.atomic_add(force_residual, 0, force_on_filament[t] + force_on_nucleus[t])


@dataclass(frozen=True, slots=True)
class WarpLincMechanics:
    """Concrete :class:`LincMechanics` binding the segment/triangle pair kernel."""

    def accumulate(
        self,
        joints: LincPopulationView,
        filament: FilamentGeometryView,
        nucleus: NuclearSurfaceView,
    ) -> None:
        _linc_launch(
            _linc_pair_force_kernel,
            joints.capacity,
            inputs=[
                filament.position_d,
                filament.force_d,
                filament.segments_d,
                nucleus.surface_position_d,
                nucleus.faces_d,
                joints.active_d,
                joints.state_d,
                joints.filament_element_id_d,
                joints.filament_u_d,
                joints.socket_face_id_d,
                joints.socket_barycentric_d,
                joints.rest_length_d,
                joints.stiffness_d,
                joints.stiffening_d,
            ],
            outputs=[
                joints.load_d,
                joints.energy_d,
                joints.force_on_filament_d,
                joints.force_on_nucleus_d,
            ],
            device=str(joints.active_d.device),
        )


@dataclass(frozen=True, slots=True)
class WarpNativeSocketProjector:
    """Concrete :class:`NuclearSocketProjector` for the native barycentric backend."""

    def accumulate_reaction(
        self,
        joints: LincPopulationView,
        nucleus: NuclearSurfaceView,
    ) -> None:
        if nucleus.mode is not NuclearProjectionMode.NATIVE_SURFACE:
            raise ValueError(
                "native socket projector requires NuclearProjectionMode.NATIVE_SURFACE; "
                "reduced J.T projection is a Core-Jacobian seam gap"
            )
        _linc_launch(
            _linc_native_socket_scatter_kernel,
            joints.capacity,
            inputs=[
                joints.active_d,
                joints.state_d,
                joints.socket_face_id_d,
                joints.socket_barycentric_d,
                joints.force_on_nucleus_d,
                nucleus.faces_d,
                nucleus.surface_force_d,
            ],
            outputs=None,
            device=str(joints.active_d.device),
        )


@dataclass(frozen=True, slots=True)
class WarpLincGenerationValidator:
    """Concrete :class:`LincPortGenerationValidator`; latches stale ports on device.

    The failure latch is a device counter set inside the physical loop and read
    only at closeout, never inside the loop.
    """

    stale_count_d: wp.array

    def __post_init__(self) -> None:
        _validate_device_array(
            self.stale_count_d, label="generation.stale_count_d", dtype=wp.int32, shape0=1
        )

    def validate(
        self,
        joints: LincPopulationView,
        filament: FilamentGeometryView,
        nucleus: NuclearSurfaceView,
    ) -> None:
        _linc_launch(
            _linc_generation_validate_kernel,
            joints.capacity,
            inputs=[
                joints.active_d,
                joints.filament_actor_generation_d,
                joints.filament_entity_generation_d,
                joints.filament_element_id_d,
                joints.nucleus_actor_generation_d,
                joints.nucleus_entity_generation_d,
                joints.socket_face_id_d,
                filament.actor_generation_d,
                filament.entity_generation_d,
                nucleus.actor_generation_d,
                nucleus.entity_generation_d,
                wp.int32(filament.record.n_elements),
                wp.int32(nucleus.record.n_elements),
                self.stale_count_d,
            ],
            outputs=None,
            device=str(joints.active_d.device),
        )


@dataclass(frozen=True, slots=True)
class SourceGatedLincKinetics:
    """Blocked :class:`LincKinetics`: declares channels, refuses to fabricate rates.

    The force-dependent bind/unbind on/off-rate cards for the three LINC edges are
    unresolved source/PI gaps.  This delegate satisfies the population's channel
    contract but raises rather than invent a hazard law, keeping the LC-3 kinetics
    milestone honestly blocked instead of filled with convenient constants.
    """

    event_channels: frozenset[str] = LINC_EVENT_CHANNELS

    def propose_candidate(self, joints: LincPopulationView, dt_phys: float, rng_seed: int) -> None:
        raise NotImplementedError(
            "LINC bind/unbind kinetics is source-gated: no validated on/off-rate card "
            "exists for nesprin/plectin LINC edges (LC-3). Surface to PI before binding."
        )


@dataclass(frozen=True, slots=True)
class WarpLincTransaction:
    """Concrete :class:`LincTransaction` binding the commit/rollback state machine.

    Snapshot is a device-to-device copy of committed state/age/RNG/population epoch
    (and a reset of the scratch candidate).  Commit and rollback select the accepted
    and rejected branches on device without an authoritative host read.  All mutable
    and remap-capable population arrays are reported through :meth:`owned_arrays` so
    the population's coverage guard is satisfied even though topology remap is a
    later (LC-3) milestone.
    """

    active_d: wp.array
    state_d: wp.array
    candidate_state_d: wp.array
    snapshot_state_d: wp.array
    age_d: wp.array
    snapshot_age_d: wp.array
    rng_epoch_d: wp.array
    snapshot_rng_epoch_d: wp.array
    population_epoch_d: wp.array
    snapshot_population_epoch_d: wp.array
    remap_arrays: tuple[wp.array, ...]
    event_channels: frozenset[str] = LINC_EVENT_CHANNELS

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return every authoritative array a commit/rollback/remap may mutate."""
        return (
            self.active_d,
            self.state_d,
            self.candidate_state_d,
            self.snapshot_state_d,
            self.age_d,
            self.snapshot_age_d,
            self.rng_epoch_d,
            self.snapshot_rng_epoch_d,
            self.population_epoch_d,
            self.snapshot_population_epoch_d,
            *self.remap_arrays,
        )

    def snapshot_candidate(self) -> None:
        """Capture committed state on device and reset the scratch candidate."""
        wp.copy(self.snapshot_state_d, self.state_d)
        wp.copy(self.snapshot_age_d, self.age_d)
        wp.copy(self.snapshot_rng_epoch_d, self.rng_epoch_d)
        wp.copy(self.snapshot_population_epoch_d, self.population_epoch_d)
        wp.copy(self.candidate_state_d, self.state_d)

    def rollback(self, accepted: wp.array) -> None:
        _linc_launch(
            _linc_rollback_kernel,
            int(self.active_d.shape[0]),
            inputs=[accepted, self.state_d, self.candidate_state_d],
            outputs=None,
            device=str(self.active_d.device),
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if rng_seed < 0:
            raise ValueError("rng_seed must be nonnegative")
        _linc_launch(
            _linc_commit_kernel,
            int(self.active_d.shape[0]),
            inputs=[
                self.active_d,
                accepted,
                self.candidate_state_d,
                self.state_d,
                self.age_d,
                self.rng_epoch_d,
                self.population_epoch_d,
                wp.float64(dt_phys),
            ],
            outputs=None,
            device=str(self.active_d.device),
        )


@dataclass(frozen=True, slots=True)
class WarpLincLedgerReducer:
    """Concrete :class:`LincLedgerContributor` reducing per-edge closeout terms.

    Reduces into its own device accumulators (population counts, elastic energy,
    equal-opposite force residual).  The global-ledger merge target does not exist
    yet (plan gap 5); the passed ``ledger`` is accepted for the future merge but the
    per-edge device reduction is authoritative and read only at closeout.
    """

    active_count_d: wp.array
    bound_count_d: wp.array
    energy_sum_d: wp.array
    force_residual_d: wp.array

    def __post_init__(self) -> None:
        _validate_device_array(self.active_count_d, label="ledger.active_count_d", dtype=wp.int32, shape0=1)
        _validate_device_array(self.bound_count_d, label="ledger.bound_count_d", dtype=wp.int32, shape0=1)
        _validate_device_array(self.energy_sum_d, label="ledger.energy_sum_d", dtype=wp.float64, shape0=1)
        _validate_device_array(self.force_residual_d, label="ledger.force_residual_d", dtype=wp.vec3d, shape0=1)

    def accumulate_linc_ledger(self, joints: LincPopulationView, ledger: object) -> None:
        _linc_launch(
            _linc_ledger_reduce_kernel,
            joints.capacity,
            inputs=[
                joints.active_d,
                joints.state_d,
                joints.energy_d,
                joints.force_on_filament_d,
                joints.force_on_nucleus_d,
                self.active_count_d,
                self.bound_count_d,
                self.energy_sum_d,
                self.force_residual_d,
            ],
            outputs=None,
            device=str(joints.active_d.device),
        )


def _linc_i32(values: list[int], device: str) -> wp.array:
    return wp.array(np.asarray(values, dtype=np.int32), dtype=wp.int32, device=device)


def _linc_i64(values: list[int], device: str) -> wp.array:
    return wp.array(np.asarray(values, dtype=np.int64), dtype=wp.int64, device=device)


def _linc_f64(values: list[float], device: str) -> wp.array:
    return wp.array(np.asarray(values, dtype=np.float64), dtype=wp.float64, device=device)


def _linc_scalar_field(value: float | npt.ArrayLike, capacity: int, name: str) -> list[float]:
    """Broadcast one required constitutive parameter to a per-joint list.

    The parameter has no default; a scalar broadcasts, a sequence must match the
    packed joint count.  This keeps LINC stiffness/rest/stiffening magnitudes an
    explicit caller responsibility (they remain source/PI gaps).
    """
    array = np.atleast_1d(np.asarray(value, dtype=np.float64))
    if array.ndim != 1:
        raise ValueError(f"LINC {name} must be a scalar or one-dimensional sequence")
    if array.size == 1:
        array = np.repeat(array, capacity)
    if array.size != capacity:
        raise ValueError(f"LINC {name} must broadcast to {capacity} joints")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"LINC {name} must be finite")
    return [float(v) for v in array]


def build_linc_joint_population(
    edge_name: str,
    specs: tuple[LincJointSpec, ...],
    *,
    rest_length_um: float | npt.ArrayLike,
    stiffness_pn_per_um: float | npt.ArrayLike,
    stiffening_per_um2: float | npt.ArrayLike,
    device: str,
) -> LincJointPopulation:
    """Allocate one edge's fixed-capacity CUDA SoA and wire the Warp delegates.

    The constitutive parameters are required (no defaults) and are the caller's
    responsibility until sourced LINC cards land.  This factory requires a CUDA
    device; the population's own validators reject any non-CUDA storage.
    """
    if edge_name not in LINC_EDGES:
        raise ValueError(f"unknown LINC population edge {edge_name!r}")
    if not specs:
        raise ValueError("LINC population needs at least one joint spec")
    if any(spec.edge_name != edge_name for spec in specs):
        raise ValueError("every joint spec must belong to the population edge")
    ids = tuple(spec.joint_id for spec in specs)
    if len(ids) != len(set(ids)):
        raise ValueError("LINC joint IDs must be unique")
    capacity = len(specs)
    rest = _linc_scalar_field(rest_length_um, capacity, "rest_length_um")
    stiff = _linc_scalar_field(stiffness_pn_per_um, capacity, "stiffness_pn_per_um")
    if any(k <= 0.0 for k in stiff):
        raise ValueError("LINC stiffness must be positive")
    stiffen = _linc_scalar_field(stiffening_per_um2, capacity, "stiffening_per_um2")
    if any(r < 0.0 for r in rest) or any(s < 0.0 for s in stiffen):
        raise ValueError("LINC rest length and stiffening must be nonnegative")

    with wp.ScopedDevice(device):
        bary = np.asarray(
            [list(spec.nucleus.barycentric) for spec in specs], dtype=np.float64
        )
        arrays = {
            "active_d": _linc_i32([1] * capacity, device),
            "joint_id_d": _linc_i64([spec.joint_id for spec in specs], device),
            "joint_generation_d": _linc_i32([spec.joint_generation for spec in specs], device),
            "state_d": _linc_i32([int(spec.initial_state) for spec in specs], device),
            "candidate_state_d": _linc_i32([int(spec.initial_state) for spec in specs], device),
            "snapshot_state_d": _linc_i32([int(spec.initial_state) for spec in specs], device),
            "age_d": _linc_f64([0.0] * capacity, device),
            "snapshot_age_d": _linc_f64([0.0] * capacity, device),
            "rng_epoch_d": _linc_i64([0] * capacity, device),
            "snapshot_rng_epoch_d": _linc_i64([0] * capacity, device),
            "filament_actor_id_d": _linc_i64([spec.filament.ref.actor_id for spec in specs], device),
            "filament_actor_generation_d": _linc_i32(
                [spec.filament.ref.actor_generation for spec in specs], device
            ),
            "filament_entity_id_d": _linc_i64([spec.filament.ref.entity_id for spec in specs], device),
            "filament_entity_generation_d": _linc_i32(
                [spec.filament.ref.entity_generation for spec in specs], device
            ),
            "filament_element_id_d": _linc_i32(
                [spec.filament.ref.element_id for spec in specs], device
            ),
            "filament_u_d": _linc_f64([spec.filament.ref.segment_u for spec in specs], device),
            "nucleus_actor_id_d": _linc_i64([spec.nucleus.ref.actor_id for spec in specs], device),
            "nucleus_actor_generation_d": _linc_i32(
                [spec.nucleus.ref.actor_generation for spec in specs], device
            ),
            "nucleus_entity_id_d": _linc_i64([spec.nucleus.ref.entity_id for spec in specs], device),
            "nucleus_entity_generation_d": _linc_i32(
                [spec.nucleus.ref.entity_generation for spec in specs], device
            ),
            "socket_face_id_d": _linc_i32([spec.nucleus.ref.element_id for spec in specs], device),
            "socket_barycentric_d": wp.array(bary, dtype=wp.vec3d, device=device),
            "rest_length_d": _linc_f64(rest, device),
            "stiffness_d": _linc_f64(stiff, device),
            "stiffening_d": _linc_f64(stiffen, device),
            "load_d": _linc_f64([0.0] * capacity, device),
            "energy_d": _linc_f64([0.0] * capacity, device),
            "force_on_filament_d": wp.zeros(capacity, dtype=wp.vec3d, device=device),
            "force_on_nucleus_d": wp.zeros(capacity, dtype=wp.vec3d, device=device),
            "population_epoch_d": _linc_i64([0], device),
            "snapshot_population_epoch_d": _linc_i64([0], device),
        }
        remap_arrays = tuple(
            arrays[name]
            for name in (
                "joint_id_d",
                "joint_generation_d",
                "filament_actor_id_d",
                "filament_actor_generation_d",
                "filament_entity_id_d",
                "filament_entity_generation_d",
                "filament_element_id_d",
                "filament_u_d",
                "nucleus_actor_id_d",
                "nucleus_actor_generation_d",
                "nucleus_entity_id_d",
                "nucleus_entity_generation_d",
                "socket_face_id_d",
                "socket_barycentric_d",
                "rest_length_d",
                "stiffness_d",
                "stiffening_d",
            )
        )
        transaction = WarpLincTransaction(
            active_d=arrays["active_d"],
            state_d=arrays["state_d"],
            candidate_state_d=arrays["candidate_state_d"],
            snapshot_state_d=arrays["snapshot_state_d"],
            age_d=arrays["age_d"],
            snapshot_age_d=arrays["snapshot_age_d"],
            rng_epoch_d=arrays["rng_epoch_d"],
            snapshot_rng_epoch_d=arrays["snapshot_rng_epoch_d"],
            population_epoch_d=arrays["population_epoch_d"],
            snapshot_population_epoch_d=arrays["snapshot_population_epoch_d"],
            remap_arrays=remap_arrays,
        )
        generation_validator = WarpLincGenerationValidator(
            stale_count_d=wp.zeros(1, dtype=wp.int32, device=device)
        )
        ledger = WarpLincLedgerReducer(
            active_count_d=wp.zeros(1, dtype=wp.int32, device=device),
            bound_count_d=wp.zeros(1, dtype=wp.int32, device=device),
            energy_sum_d=wp.zeros(1, dtype=wp.float64, device=device),
            force_residual_d=wp.zeros(1, dtype=wp.vec3d, device=device),
        )
    return LincJointPopulation(
        edge_name=edge_name,
        capacity=capacity,
        mechanics=WarpLincMechanics(),
        generation_validator=generation_validator,
        projector=WarpNativeSocketProjector(),
        kinetics=SourceGatedLincKinetics(),
        transaction=transaction,
        ledger=ledger,
        **arrays,
    )


__all__ = [
    "ACTIN_CAP_LINC",
    "IF_NUCLEUS_LINC",
    "LINC_EDGES",
    "LINC_EVENT_CHANNELS",
    "MT_NUCLEUS_LINC",
    "FilamentGeometryView",
    "LincFilamentPort",
    "LincFilamentRole",
    "LincJointPopulation",
    "LincRigConnectorAdapter",
    "LincJointSpec",
    "LincJointState",
    "LincLedgerReport",
    "LincPairReference",
    "LincPopulationView",
    "NuclearProjectionMode",
    "NuclearSurfaceSocket",
    "NuclearSurfaceView",
    "ReducedProjectionReference",
    "SourceGatedLincKinetics",
    "WarpLincGenerationValidator",
    "WarpLincLedgerReducer",
    "WarpLincMechanics",
    "WarpLincTransaction",
    "WarpNativeSocketProjector",
    "build_linc_joint_population",
    "reduced_jt_projection_reference",
    "resolve_linc_joint",
    "segment_triangle_linc_reference",
    "validate_linc_architecture",
]
