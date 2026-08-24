"""Structural gates for the canonical stress-fiber/arc facade."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine.load_path import ActorRecord
from aleph.engine.runtime import CytosolFieldEndpoint
from aleph.engine.stress_fiber import (
    ACTIN_CAP_LINC,
    DORSAL_ARC_CROSSLINK,
    SF_COMPONENT,
    SF_CORTEX_TRANSIENT,
    SF_CYTOSOL_TRANSFER,
    StressFiberActor,
    StressFiberEndpoint,
    StressFiberStateOwner,
    StressFiberStepBindings,
)


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:test"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    ptr: int
    shape: tuple[int, ...]
    dtype: object
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _MechanicsSpy:
    calls: list[tuple[object, object]] = field(default_factory=list)
    n_nodes: int = 4

    def accumulate(self, position: object, force: object) -> None:
        self.calls.append((position, force))


@dataclass(slots=True)
class _TransactionSpy:
    arrays: tuple[object, ...]
    component_name: str = SF_COMPONENT
    calls: list[tuple[object, ...]] = field(default_factory=list)

    def owned_arrays(self) -> tuple[object, ...]:
        return self.arrays

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted: object) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))


@dataclass(slots=True)
class _LedgerSpy:
    calls: list[object] = field(default_factory=list)

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append(ledger)


@dataclass(slots=True)
class _GraphConnectorSpy:
    name: str
    component_a: str
    component_b: str
    mechanics_calls: list[tuple[object, object]] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_rig(self, rig: object, endpoint: object) -> None:
        self.mechanics_calls.append((rig, endpoint))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: object) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _InternalConnectorSpy(_GraphConnectorSpy):
    def accumulate_internal(self, rig: object) -> None:
        self.mechanics_calls.append((rig, rig))


@dataclass(slots=True)
class _TransferSpy(_GraphConnectorSpy):
    def accumulate_transfer(self, solid: object, cytosol: object) -> None:
        self.mechanics_calls.append((solid, cytosol))


def _state() -> tuple[StressFiberStateOwner, _MechanicsSpy, _TransactionSpy, _LedgerSpy]:
    position = _FakeArray(100, (4,), wp.vec3d)
    force = _FakeArray(101, (4,), wp.vec3d)
    segments = _FakeArray(102, (2, 2), wp.int32)
    actor_generation = _FakeArray(103, (1,), wp.int32)
    entity_generation = _FakeArray(104, (1,), wp.int32)
    topology_epoch = _FakeArray(105, (1,), wp.int32)
    mechanics = _MechanicsSpy()
    transaction = _TransactionSpy(
        (position, segments, actor_generation, entity_generation, topology_epoch)
    )
    ledger = _LedgerSpy()
    state = StressFiberStateOwner(
        ActorRecord(SF_COMPONENT, 7, 2, 11, 3, 2),
        position,
        force,
        segments,
        actor_generation,
        entity_generation,
        topology_epoch,
        mechanics,
        transaction,
        ledger,
    )
    return state, mechanics, transaction, ledger


def _bindings() -> tuple[
    StressFiberStepBindings,
    _GraphConnectorSpy,
    _GraphConnectorSpy,
    _TransferSpy,
    _InternalConnectorSpy,
]:
    cortex = StressFiberEndpoint(
        "cortex",
        _FakeArray(200, (3,), wp.vec3d),
        _FakeArray(201, (3,), wp.vec3d),
    )
    nucleus = StressFiberEndpoint(
        "nucleus",
        _FakeArray(210, (3,), wp.vec3d),
        _FakeArray(211, (3,), wp.vec3d),
    )
    cytosol = CytosolFieldEndpoint(
        "cytosol",
        _FakeArray(220, (2, 2, 2), wp.float64),
        _FakeArray(221, (2, 2, 2), wp.float64),
    )
    cortex_connector = _GraphConnectorSpy(SF_CORTEX_TRANSIENT, SF_COMPONENT, "cortex")
    linc = _GraphConnectorSpy(ACTIN_CAP_LINC, SF_COMPONENT, "nucleus")
    transfer = _TransferSpy(SF_CYTOSOL_TRANSFER, SF_COMPONENT, "cytosol")
    dorsal = _InternalConnectorSpy(DORSAL_ARC_CROSSLINK, SF_COMPONENT, SF_COMPONENT)
    return (
        StressFiberStepBindings(cortex, nucleus, cytosol, cortex_connector, linc, transfer, dorsal),
        cortex_connector,
        linc,
        transfer,
        dorsal,
    )


def test_reference_architecture_exposes_all_four_canonical_sf_edges() -> None:
    state, *_ = _state()
    actor = StressFiberActor(reference_cell_architecture(), state)
    contracts = {connector.name: connector for connector in actor.architecture.connectors}

    assert contracts[SF_CORTEX_TRANSIENT].family is ConnectorFamily.TRANSIENT_ACTIN
    assert contracts[ACTIN_CAP_LINC].family is ConnectorFamily.LINC
    assert contracts[SF_CYTOSOL_TRANSFER].family is ConnectorFamily.IMMERSED_TRANSFER
    assert contracts[DORSAL_ARC_CROSSLINK].component_a == contracts[DORSAL_ARC_CROSSLINK].component_b


def test_state_requires_distinct_cuda_storage_and_complete_transaction() -> None:
    state, *_ = _state()
    with pytest.raises(ValueError, match="must not alias"):
        replace(state, force_d=state.position_d)
    incomplete = replace(state.transaction, arrays=state.transaction.arrays[:-1])
    with pytest.raises(ValueError, match="must cover"):
        replace(state, transaction=incomplete)


def test_bindings_reject_missing_or_misdirected_graph_connectors() -> None:
    bindings, *_ = _bindings()
    with pytest.raises(ValueError, match=ACTIN_CAP_LINC):
        replace(bindings, actin_cap_linc=replace(bindings.actin_cap_linc, name="private_linc"))
    with pytest.raises(ValueError, match="incorrect component endpoints"):
        replace(
            bindings,
            cytosol_transfer=replace(bindings.cytosol_transfer, component_a="microtubule"),
        )


def test_mechanics_reaches_internal_cortex_linc_and_cytosol_paths_once() -> None:
    state, mechanics, *_ = _state()
    bindings, cortex, linc, transfer, dorsal = _bindings()
    actor = StressFiberActor(reference_cell_architecture(), state)
    actor.accumulate_mechanics(bindings)

    assert mechanics.calls == [(state.position_d, state.force_d)]
    assert len(dorsal.mechanics_calls) == 1
    assert len(cortex.mechanics_calls) == 1
    assert len(linc.mechanics_calls) == 1
    assert len(transfer.mechanics_calls) == 1
    rig = cortex.mechanics_calls[0][0]
    assert rig.position_d is state.position_d and rig.force_d is state.force_d
    assert transfer.mechanics_calls[0][1] is bindings.cytosol


def test_transaction_and_ledgers_cover_state_and_all_four_connectors() -> None:
    state, _, transaction, ledger_spy = _state()
    bindings, cortex, linc, transfer, dorsal = _bindings()
    actor = StressFiberActor(reference_cell_architecture(), state)
    accepted = object()
    ledger = object()

    actor.snapshot_candidate(bindings)
    actor.rollback(bindings, accepted)
    actor.commit_irreversible(bindings, accepted, 0.02, 17)
    actor.accumulate_ledger(bindings, ledger)

    expected = [("snapshot",), ("rollback", accepted), ("commit", accepted, 0.02, 17)]
    assert transaction.calls == expected
    assert all(
        connector.transaction_calls == expected for connector in (dorsal, cortex, linc, transfer)
    )
    assert ledger_spy.calls == [ledger]
    assert all(connector.ledger_calls == [ledger] for connector in (dorsal, cortex, linc, transfer))


def test_sf_facade_does_not_own_fa_or_nmii_dispatch() -> None:
    bindings, *_ = _bindings()
    assert set(type(bindings).__slots__) == {
        "cortex",
        "nucleus",
        "cytosol",
        "cortex_transient",
        "actin_cap_linc",
        "cytosol_transfer",
        "dorsal_arc_crosslink",
    }


def test_wrong_common_family_is_rejected_before_runtime() -> None:
    state, *_ = _state()
    architecture = reference_cell_architecture()
    connectors = tuple(
        replace(connector, family=ConnectorFamily.CONTACT)
        if connector.name == SF_CORTEX_TRANSIENT
        else connector
        for connector in architecture.connectors
    )
    with pytest.raises(ValueError, match="incorrect family or scope"):
        StressFiberActor(CellArchitecture(architecture.components, connectors), state)
