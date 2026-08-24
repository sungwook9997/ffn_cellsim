"""Canonical component facade for stress fibers and transverse/dorsal arcs.

The load-path module supplies typed segment joints, but without a state-owning ``sf_arc`` caller several
common graph edges were unreachable from a whole-cell candidate.  This facade owns only SF/arc geometry and
internal state.  Cortex coupling, actin-cap LINC, porous-fluid transfer, and the dynamic dorsal--arc joint are
graph-owned bindings and are dispatched exactly once here.

FA--ECM series mechanics and NMII crossbridges are intentionally absent: their canonical dispatch owners are
the ECM composite-clutch and NMII actuator facades.  This prevents the SF actor from applying either load a
second time.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
    ConnectorScope,
)
from aleph.engine.load_path import ActorRecord, BorrowedSegmentActorView
from aleph.engine.runtime import (
    CytosolFieldEndpoint,
    ImmersedTransferConnector,
    LedgerContributor,
    MechanicsContributor,
    TransactionParticipant,
)

SF_COMPONENT = "sf_arc"
CORTEX_COMPONENT = "cortex"
NUCLEUS_COMPONENT = "nucleus"
CYTOSOL_COMPONENT = "cytosol"
SF_CORTEX_TRANSIENT = "sf_cortex_transient"
ACTIN_CAP_LINC = "actin_cap_linc"
SF_CYTOSOL_TRANSFER = "sf_cytosol_transfer"
DORSAL_ARC_CROSSLINK = "dorsal_arc_crosslink"


def _storage_key(array: object) -> tuple[object, ...]:
    """Return device storage identity without reading device state."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


def _validate_i32_vector(array: object, *, label: str, length: int) -> None:
    device = getattr(array, "device", None)
    if not bool(getattr(device, "is_cuda", False)):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != wp.int32:
        raise TypeError(f"{label} must have dtype wp.int32")
    if getattr(array, "shape", None) != (length,):
        raise ValueError(f"{label} must have shape {(length,)}")


def _validate_vec3_pair(position: object, force: object, *, label: str) -> None:
    for suffix, array in (("position_d", position), ("force_d", force)):
        device = getattr(array, "device", None)
        if not bool(getattr(device, "is_cuda", False)):
            raise ValueError(f"{label}.{suffix} must be a Warp CUDA device array")
        if getattr(array, "dtype", None) != wp.vec3d:
            raise TypeError(f"{label}.{suffix} must have dtype wp.vec3d")
        shape = getattr(array, "shape", None)
        if not isinstance(shape, tuple) or len(shape) != 1 or int(shape[0]) <= 0:
            raise ValueError(f"{label}.{suffix} must be a non-empty vector")
    if position.shape != force.shape:
        raise ValueError(f"{label} position and force arrays must have identical shape")
    if str(position.device) != str(force.device):
        raise ValueError(f"{label} position and force arrays must share one CUDA device")
    if _storage_key(position) == _storage_key(force):
        raise ValueError(f"{label} position and force arrays must not alias")


@dataclass(frozen=True, slots=True)
class StressFiberRigView:
    """Non-owning live segment geometry and generation view for registered connectors."""

    record: ActorRecord
    position_d: wp.array
    force_d: wp.array
    segments_d: wp.array
    actor_generation_d: wp.array
    entity_generation_d: wp.array
    topology_epoch_d: wp.array


@dataclass(frozen=True, slots=True)
class StressFiberEndpoint:
    """Non-owning cortex or nuclear mechanical endpoint."""

    component: str
    position_d: wp.array
    force_d: wp.array

    def __post_init__(self) -> None:
        if self.component not in {CORTEX_COMPONENT, NUCLEUS_COMPONENT}:
            raise ValueError("SF endpoint must be cortex or nucleus")
        _validate_vec3_pair(self.position_d, self.force_d, label=self.component)


@runtime_checkable
class StressFiberGraphConnector(TransactionParticipant, LedgerContributor, Protocol):
    """Graph-owned kinetic SF connector to cortex or nuclear surface."""

    name: str
    component_a: str
    component_b: str

    def accumulate_rig(self, rig: StressFiberRigView, endpoint: StressFiberEndpoint) -> None:
        """Scatter equal-and-opposite connector mechanics through public views."""


@runtime_checkable
class InternalStressFiberConnector(TransactionParticipant, LedgerContributor, Protocol):
    """Graph-owned dorsal-free--transverse-arc joint population."""

    name: str
    component_a: str
    component_b: str

    def accumulate_internal(self, rig: StressFiberRigView) -> None:
        """Add internal dynamic arc-joint mechanics to SF-owned force storage."""


@runtime_checkable
class StressFiberTransaction(TransactionParticipant, Protocol):
    """Transaction delegate proving coverage of evolving SF geometry and topology."""

    component_name: str

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return every authoritative SF array covered by snapshot/rollback."""


@dataclass(frozen=True, slots=True)
class StressFiberStateOwner:
    """Authoritative component-local CUDA state for SF and arc segment graphs."""

    record: ActorRecord
    position_d: wp.array
    force_d: wp.array
    segments_d: wp.array
    actor_generation_d: wp.array
    entity_generation_d: wp.array
    topology_epoch_d: wp.array
    mechanics: MechanicsContributor
    transaction: StressFiberTransaction
    ledger: LedgerContributor

    def __post_init__(self) -> None:
        if self.record.component != SF_COMPONENT:
            raise ValueError("stress-fiber state owner record must belong to sf_arc")
        BorrowedSegmentActorView(self.record, self.position_d, self.force_d, self.segments_d)
        for label, array in (
            ("actor_generation_d", self.actor_generation_d),
            ("entity_generation_d", self.entity_generation_d),
            ("topology_epoch_d", self.topology_epoch_d),
        ):
            _validate_i32_vector(array, label=f"sf_arc.{label}", length=1)
        all_arrays = (
            self.position_d,
            self.force_d,
            self.segments_d,
            self.actor_generation_d,
            self.entity_generation_d,
            self.topology_epoch_d,
        )
        if len({str(array.device) for array in all_arrays}) != 1:
            raise ValueError("all SF arrays must share one CUDA device")
        if len({_storage_key(array) for array in all_arrays}) != len(all_arrays):
            raise ValueError("SF geometry, force, generation, and topology arrays must not alias")
        if self.transaction.component_name != SF_COMPONENT:
            raise ValueError("SF transaction delegate must be owned by sf_arc")
        authoritative = (
            self.position_d,
            self.segments_d,
            self.actor_generation_d,
            self.entity_generation_d,
            self.topology_epoch_d,
        )
        covered = {_storage_key(array) for array in self.transaction.owned_arrays()}
        if {_storage_key(array) for array in authoritative} - covered:
            raise ValueError("SF transaction must cover geometry, generations, and topology")
        n_nodes = getattr(self.mechanics, "n_nodes", None)
        if n_nodes is not None and int(n_nodes) != int(self.position_d.shape[0]):
            raise ValueError("SF mechanics node count does not match owned geometry")

    def geometry(self) -> StressFiberRigView:
        """Expose live SF material-coordinate geometry without transferring ownership."""
        return StressFiberRigView(
            self.record,
            self.position_d,
            self.force_d,
            self.segments_d,
            self.actor_generation_d,
            self.entity_generation_d,
            self.topology_epoch_d,
        )

    def accumulate(self) -> None:
        """Launch component-local SF/arc mechanics."""
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
class StressFiberStepBindings:
    """Exactly the four common edges canonically dispatched by the SF facade."""

    cortex: StressFiberEndpoint
    nucleus: StressFiberEndpoint
    cytosol: CytosolFieldEndpoint
    cortex_transient: StressFiberGraphConnector
    actin_cap_linc: StressFiberGraphConnector
    cytosol_transfer: ImmersedTransferConnector
    dorsal_arc_crosslink: InternalStressFiberConnector

    def __post_init__(self) -> None:
        if self.cortex.component != CORTEX_COMPONENT:
            raise ValueError("SF cortex endpoint must be the registered cortex")
        if self.nucleus.component != NUCLEUS_COMPONENT:
            raise ValueError("SF nucleus endpoint must be the registered nucleus")
        _validate_runtime_connector(
            self.cortex_transient,
            name=SF_CORTEX_TRANSIENT,
            endpoints={SF_COMPONENT, CORTEX_COMPONENT},
        )
        _validate_runtime_connector(
            self.actin_cap_linc,
            name=ACTIN_CAP_LINC,
            endpoints={SF_COMPONENT, NUCLEUS_COMPONENT},
        )
        _validate_runtime_connector(
            self.cytosol_transfer,
            name=SF_CYTOSOL_TRANSFER,
            endpoints={SF_COMPONENT, CYTOSOL_COMPONENT},
        )
        _validate_runtime_connector(
            self.dorsal_arc_crosslink,
            name=DORSAL_ARC_CROSSLINK,
            endpoints={SF_COMPONENT},
        )


def _validate_runtime_connector(connector: object, *, name: str, endpoints: set[str]) -> None:
    if getattr(connector, "name", None) != name:
        raise ValueError(f"expected graph connector {name!r}")
    if {
        getattr(connector, "component_a", None),
        getattr(connector, "component_b", None),
    } != endpoints:
        raise ValueError(f"{name} has incorrect component endpoints")


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
    scope: ConnectorScope = ConnectorScope.INTER_COMPONENT,
    kinetic: bool,
) -> None:
    connector = _find_connector(architecture, name)
    if connector.family is not family or connector.scope is not scope:
        raise ValueError(f"{name} has incorrect family or scope")
    if {connector.component_a, connector.component_b} != endpoints:
        raise ValueError(f"{name} has incorrect component endpoints")
    if connector.kinetics is not kinetic or connector.commit_on_accept is not kinetic:
        raise ValueError(f"{name} has incorrect accepted-step kinetics")
    if not connector.bidirectional or not connector.adjoint_transfer_required:
        raise ValueError(f"{name} must be bidirectional with adjoint transfer")


def _validate_architecture(architecture: CellArchitecture) -> None:
    component = architecture.component(SF_COMPONENT)
    if component.role is not ComponentRole.ACTIVE_LOAD_PATH:
        raise ValueError("sf_arc must have ComponentRole.ACTIVE_LOAD_PATH")
    if not component.owns_geometry or not component.dynamically_evolving:
        raise ValueError("sf_arc must own evolving geometry")
    _validate_build_connector(
        architecture,
        name=SF_CORTEX_TRANSIENT,
        family=ConnectorFamily.TRANSIENT_ACTIN,
        endpoints={SF_COMPONENT, CORTEX_COMPONENT},
        kinetic=True,
    )
    _validate_build_connector(
        architecture,
        name=ACTIN_CAP_LINC,
        family=ConnectorFamily.LINC,
        endpoints={SF_COMPONENT, NUCLEUS_COMPONENT},
        kinetic=True,
    )
    _validate_build_connector(
        architecture,
        name=SF_CYTOSOL_TRANSFER,
        family=ConnectorFamily.IMMERSED_TRANSFER,
        endpoints={SF_COMPONENT, CYTOSOL_COMPONENT},
        kinetic=False,
    )
    _validate_build_connector(
        architecture,
        name=DORSAL_ARC_CROSSLINK,
        family=ConnectorFamily.TRANSIENT_ACTIN,
        endpoints={SF_COMPONENT},
        scope=ConnectorScope.INTERNAL,
        kinetic=True,
    )


@dataclass(frozen=True, slots=True)
class StressFiberActor:
    """Stateless orchestration facade for SF/arc mechanics and its canonical connector edges."""

    architecture: CellArchitecture
    state: StressFiberStateOwner

    def __post_init__(self) -> None:
        _validate_architecture(self.architecture)

    def accumulate_mechanics(self, bindings: StressFiberStepBindings) -> None:
        """Dispatch component, internal joint, cortex, LINC, and porous-fluid mechanics once."""
        self.state.accumulate()
        rig = self.state.geometry()
        bindings.dorsal_arc_crosslink.accumulate_internal(rig)
        bindings.cortex_transient.accumulate_rig(rig, bindings.cortex)
        bindings.actin_cap_linc.accumulate_rig(rig, bindings.nucleus)
        bindings.cytosol_transfer.accumulate_transfer(rig, bindings.cytosol)

    def transaction_participants(
        self, bindings: StressFiberStepBindings
    ) -> tuple[TransactionParticipant, ...]:
        """Return every SF-owned or incident canonical transaction participant."""
        return (
            self.state,
            bindings.dorsal_arc_crosslink,
            bindings.cortex_transient,
            bindings.actin_cap_linc,
            bindings.cytosol_transfer,
        )

    def snapshot_candidate(self, bindings: StressFiberStepBindings) -> None:
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: StressFiberStepBindings, accepted: wp.array) -> None:
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self,
        bindings: StressFiberStepBindings,
        accepted: wp.array,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, bindings: StressFiberStepBindings, ledger: object) -> None:
        self.state.accumulate_ledger(ledger)
        bindings.dorsal_arc_crosslink.accumulate_ledger(ledger)
        bindings.cortex_transient.accumulate_ledger(ledger)
        bindings.actin_cap_linc.accumulate_ledger(ledger)
        bindings.cytosol_transfer.accumulate_ledger(ledger)


__all__ = [
    "ACTIN_CAP_LINC",
    "CORTEX_COMPONENT",
    "CYTOSOL_COMPONENT",
    "DORSAL_ARC_CROSSLINK",
    "InternalStressFiberConnector",
    "NUCLEUS_COMPONENT",
    "SF_COMPONENT",
    "SF_CORTEX_TRANSIENT",
    "SF_CYTOSOL_TRANSFER",
    "StressFiberActor",
    "StressFiberEndpoint",
    "StressFiberGraphConnector",
    "StressFiberRigView",
    "StressFiberStateOwner",
    "StressFiberStepBindings",
    "StressFiberTransaction",
]
