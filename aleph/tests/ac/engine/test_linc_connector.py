"""Structural and analytic gates for the graph-owned LINC connector family."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine import linc_connector as linc_module
from aleph.engine.linc_connector import (
    ACTIN_CAP_LINC,
    IF_NUCLEUS_LINC,
    LINC_EDGES,
    LINC_EVENT_CHANNELS,
    MT_NUCLEUS_LINC,
    FilamentGeometryView,
    LincFilamentPort,
    LincFilamentRole,
    LincJointPopulation,
    LincJointSpec,
    LincJointState,
    LincKinetics,
    LincLedgerContributor,
    LincLedgerReport,
    LincMechanics,
    LincPortGenerationValidator,
    LincRigConnectorAdapter,
    LincTransaction,
    NuclearProjectionMode,
    NuclearSocketProjector,
    NuclearSurfaceSocket,
    NuclearSurfaceView,
    SourceGatedLincKinetics,
    WarpLincGenerationValidator,
    WarpLincLedgerReducer,
    WarpLincMechanics,
    WarpLincTransaction,
    WarpNativeSocketProjector,
    build_linc_joint_population,
    reduced_jt_projection_reference,
    resolve_linc_joint,
    segment_triangle_linc_reference,
    validate_linc_architecture,
)
from aleph.engine.load_path import (
    ActorRecord,
    ActorRegistry,
    ElementKind,
    EndpointRole,
    PortRef,
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
    calls: list[tuple[object, ...]]

    def accumulate(self, joints: object, filament: object, nucleus: object) -> None:
        self.calls.append((joints, filament, nucleus))


@dataclass(slots=True)
class _GenerationSpy:
    calls: list[tuple[object, ...]]

    def validate(self, joints: object, filament: object, nucleus: object) -> None:
        self.calls.append((joints, filament, nucleus))


@dataclass(slots=True)
class _ProjectorSpy:
    calls: list[tuple[object, ...]]

    def accumulate_reaction(self, joints: object, nucleus: object) -> None:
        self.calls.append((joints, nucleus))


@dataclass(slots=True)
class _KineticsSpy:
    calls: list[tuple[object, ...]]
    event_channels: frozenset[str] = LINC_EVENT_CHANNELS

    def propose_candidate(self, joints: object, dt_phys: float, rng_seed: int) -> None:
        self.calls.append((joints, dt_phys, rng_seed))


@dataclass(slots=True)
class _TransactionSpy:
    arrays: tuple[object, ...]
    calls: list[tuple[object, ...]]
    event_channels: frozenset[str] = LINC_EVENT_CHANNELS

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
    calls: list[tuple[object, object]]

    def accumulate_linc_ledger(self, joints: object, ledger: object) -> None:
        self.calls.append((joints, ledger))


def _port(
    record: ActorRecord,
    *,
    kind: ElementKind,
    coordinates: tuple[float, float, float, float],
) -> PortRef:
    return PortRef(
        component=record.component,
        actor_id=record.actor_id,
        actor_generation=record.actor_generation,
        entity_id=record.entity_id,
        entity_generation=record.entity_generation,
        element_kind=kind,
        element_id=0,
        local_coordinates=coordinates,
        role=EndpointRole.TRANSVERSE_ARC_MATERIAL,
    )


def _spec(edge_name: str) -> tuple[LincJointSpec, ActorRegistry]:
    component = {
        ACTIN_CAP_LINC: "sf_arc",
        MT_NUCLEUS_LINC: "microtubule",
        IF_NUCLEUS_LINC: "intermediate_filament",
    }[edge_name]
    role = {
        ACTIN_CAP_LINC: LincFilamentRole.PERINUCLEAR_ACTIN_CAP,
        MT_NUCLEUS_LINC: LincFilamentRole.MT_PLUS_END_OR_LATTICE,
        IF_NUCLEUS_LINC: LincFilamentRole.IF_JUNCTION_OR_MATERIAL,
    }[edge_name]
    filament_record = ActorRecord(component, 10, 3, 100, 4, 1)
    nucleus_record = ActorRecord("nucleus", 20, 5, 200, 6, 1)
    registry = ActorRegistry((filament_record, nucleus_record))
    filament = LincFilamentPort(
        _port(
            filament_record,
            kind=ElementKind.SEGMENT,
            coordinates=(0.25, 0.0, 0.0, 0.0),
        ),
        role,
    )
    nucleus = NuclearSurfaceSocket(
        _port(
            nucleus_record,
            kind=ElementKind.TRIANGLE,
            coordinates=(0.2, 0.3, 0.5, 0.0),
        )
    )
    return LincJointSpec(7, 2, edge_name, filament, nucleus), registry


def _population(
    edge_name: str = ACTIN_CAP_LINC,
) -> tuple[
    LincJointPopulation,
    _MechanicsSpy,
    _GenerationSpy,
    _ProjectorSpy,
    _KineticsSpy,
    _TransactionSpy,
    _LedgerSpy,
]:
    capacity = 2
    ptr = iter(range(100, 200))

    def array(dtype: object, shape: tuple[int, ...] = (capacity,)) -> _FakeArray:
        return _FakeArray(next(ptr), shape, dtype)

    arrays = {
        "active_d": array(wp.int32),
        "joint_id_d": array(wp.int64),
        "joint_generation_d": array(wp.int32),
        "state_d": array(wp.int32),
        "candidate_state_d": array(wp.int32),
        "snapshot_state_d": array(wp.int32),
        "age_d": array(wp.float64),
        "snapshot_age_d": array(wp.float64),
        "rng_epoch_d": array(wp.int64),
        "snapshot_rng_epoch_d": array(wp.int64),
        "filament_actor_id_d": array(wp.int64),
        "filament_actor_generation_d": array(wp.int32),
        "filament_entity_id_d": array(wp.int64),
        "filament_entity_generation_d": array(wp.int32),
        "filament_element_id_d": array(wp.int32),
        "filament_u_d": array(wp.float64),
        "nucleus_actor_id_d": array(wp.int64),
        "nucleus_actor_generation_d": array(wp.int32),
        "nucleus_entity_id_d": array(wp.int64),
        "nucleus_entity_generation_d": array(wp.int32),
        "socket_face_id_d": array(wp.int32),
        "socket_barycentric_d": array(wp.vec3d),
        "rest_length_d": array(wp.float64),
        "stiffness_d": array(wp.float64),
        "stiffening_d": array(wp.float64),
        "load_d": array(wp.float64),
        "energy_d": array(wp.float64),
        "force_on_filament_d": array(wp.vec3d),
        "force_on_nucleus_d": array(wp.vec3d),
        "population_epoch_d": array(wp.int64, (1,)),
        "snapshot_population_epoch_d": array(wp.int64, (1,)),
    }
    mutable_names = (
        "active_d",
        "joint_id_d",
        "joint_generation_d",
        "state_d",
        "candidate_state_d",
        "snapshot_state_d",
        "age_d",
        "snapshot_age_d",
        "rng_epoch_d",
        "snapshot_rng_epoch_d",
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
        "population_epoch_d",
        "snapshot_population_epoch_d",
    )
    mechanics = _MechanicsSpy([])
    generation = _GenerationSpy([])
    projector = _ProjectorSpy([])
    kinetics = _KineticsSpy([])
    transaction = _TransactionSpy(tuple(arrays[name] for name in mutable_names), [])
    ledger = _LedgerSpy([])
    population = LincJointPopulation(
        edge_name=edge_name,
        capacity=capacity,
        mechanics=mechanics,
        generation_validator=generation,
        projector=projector,
        kinetics=kinetics,
        transaction=transaction,
        ledger=ledger,
        **arrays,
    )
    return population, mechanics, generation, projector, kinetics, transaction, ledger


def _geometry(
    component: str,
    *,
    mode: NuclearProjectionMode = NuclearProjectionMode.NATIVE_SURFACE,
) -> tuple[FilamentGeometryView, NuclearSurfaceView]:
    filament_record = ActorRecord(component, 10, 3, 100, 4, 1)
    nucleus_record = ActorRecord("nucleus", 20, 5, 200, 6, 1)
    filament = FilamentGeometryView(
        filament_record,
        _FakeArray(300, (4,), wp.vec3d),
        _FakeArray(301, (4,), wp.vec3d),
        _FakeArray(302, (1, 2), wp.int32),
        _FakeArray(303, (1,), wp.int32),
        _FakeArray(304, (1,), wp.int32),
    )
    generalized = (
        _FakeArray(405, (8,), wp.float64)
        if mode is NuclearProjectionMode.REDUCED_JT
        else None
    )
    nucleus = NuclearSurfaceView(
        nucleus_record,
        mode,
        _FakeArray(400, (3,), wp.vec3d),
        _FakeArray(401, (3,), wp.vec3d),
        _FakeArray(402, (1, 3), wp.int32),
        _FakeArray(403, (1,), wp.int32),
        _FakeArray(404, (1,), wp.int32),
        generalized,
    )
    return filament, nucleus


@dataclass(frozen=True, slots=True)
class _CompactEndpoint:
    position_d: object
    force_d: object


def test_reference_graph_declares_exactly_three_kinetic_adjoint_linc_edges() -> None:
    architecture = reference_cell_architecture()
    validate_linc_architecture(architecture)
    selected = {edge.name: edge for edge in architecture.connectors if edge.name in LINC_EDGES}
    assert set(selected) == LINC_EDGES
    assert all(edge.family is ConnectorFamily.LINC for edge in selected.values())
    assert all(edge.kinetics and edge.commit_on_accept for edge in selected.values())
    assert all(edge.bidirectional and edge.adjoint_transfer_required for edge in selected.values())


@pytest.mark.parametrize("edge_name", sorted(LINC_EDGES))
def test_joint_resolution_validates_typed_filament_port_and_nuclear_socket(edge_name: str) -> None:
    spec, registry = _spec(edge_name)
    assert resolve_linc_joint(spec, registry, reference_cell_architecture()) is spec
    stale_ref = replace(spec.filament.ref, actor_generation=spec.filament.ref.actor_generation + 1)
    stale_spec = replace(spec, filament=replace(spec.filament, ref=stale_ref))
    with pytest.raises(ValueError, match="stale actor generation"):
        resolve_linc_joint(stale_spec, registry, reference_cell_architecture())
    wrong_role = next(role for role in LincFilamentRole if role is not spec.filament.role)
    with pytest.raises(ValueError, match="wrong filament endpoint role"):
        resolve_linc_joint(
            replace(spec, filament=replace(spec.filament, role=wrong_role)),
            registry,
            reference_cell_architecture(),
        )


def test_nuclear_socket_requires_triangle_barycentrics() -> None:
    spec, _ = _spec(ACTIN_CAP_LINC)
    bad_sum = replace(spec.nucleus.ref, local_coordinates=(0.2, 0.3, 0.4, 0.0))
    with pytest.raises(ValueError, match="sum to one"):
        NuclearSurfaceSocket(bad_sum)
    with pytest.raises(ValueError, match="surface triangle"):
        NuclearSurfaceSocket(replace(spec.nucleus.ref, element_kind=ElementKind.NODE))


def test_cuda_soa_rejects_aliasing_and_incomplete_transaction_coverage() -> None:
    population, *_, transaction, _ = _population()
    with pytest.raises(ValueError, match="distinct storage"):
        replace(population, energy_d=population.load_d)
    incomplete = replace(transaction, arrays=transaction.arrays[:-1])
    with pytest.raises(ValueError, match="all mutable/remap population arrays"):
        replace(population, transaction=incomplete)


@pytest.mark.parametrize("edge_name", [MT_NUCLEUS_LINC, IF_NUCLEUS_LINC])
def test_rig_adapter_closes_mt_if_connector_interface(edge_name: str) -> None:
    population, mechanics, generation, projector, _, transaction, ledger_spy = _population(
        edge_name
    )
    component = "microtubule" if edge_name == MT_NUCLEUS_LINC else "intermediate_filament"
    filament, nucleus = _geometry(component)
    adapter = LincRigConnectorAdapter(population, filament, nucleus)
    rig = _CompactEndpoint(filament.position_d, filament.force_d)
    endpoint = _CompactEndpoint(nucleus.surface_position_d, nucleus.surface_force_d)

    adapter.accumulate_rig(rig, endpoint)
    accepted = object()
    ledger = object()
    adapter.snapshot_candidate()
    adapter.rollback(accepted)
    adapter.commit_irreversible(accepted, 0.05, 7)
    adapter.accumulate_ledger(ledger)

    assert adapter.name == edge_name
    assert {adapter.component_a, adapter.component_b} == {component, "nucleus"}
    assert len(generation.calls) == len(mechanics.calls) == len(projector.calls) == 1
    assert transaction.calls == [
        ("snapshot",),
        ("rollback", accepted),
        ("commit", accepted, 0.05, 7),
    ]
    assert ledger_spy.calls[-1][1] is ledger


def test_rig_adapter_rejects_unbound_geometry_storage() -> None:
    population, *_ = _population(MT_NUCLEUS_LINC)
    filament, nucleus = _geometry("microtubule")
    adapter = LincRigConnectorAdapter(population, filament, nucleus)
    wrong_rig = _CompactEndpoint(_FakeArray(999, (4,), wp.vec3d), filament.force_d)
    endpoint = _CompactEndpoint(nucleus.surface_position_d, nucleus.surface_force_d)
    with pytest.raises(ValueError, match="bound filament geometry"):
        adapter.accumulate_rig(wrong_rig, endpoint)


@pytest.mark.parametrize(
    ("edge_name", "component"),
    [
        (ACTIN_CAP_LINC, "sf_arc"),
        (MT_NUCLEUS_LINC, "microtubule"),
        (IF_NUCLEUS_LINC, "intermediate_filament"),
    ],
)
def test_common_population_facade_forwards_generation_mechanics_and_projection(
    edge_name: str,
    component: str,
) -> None:
    population, mechanics, generation, projector, *_ = _population(edge_name)
    filament, nucleus = _geometry(component)
    population.accumulate(filament, nucleus)
    assert len(generation.calls) == len(mechanics.calls) == len(projector.calls) == 1
    assert generation.calls[0][0].edge_name == edge_name
    assert projector.calls[0][1].mode is NuclearProjectionMode.NATIVE_SURFACE

    _, reduced = _geometry(component, mode=NuclearProjectionMode.REDUCED_JT)
    population.accumulate(filament, reduced)
    assert projector.calls[-1][1].generalized_force_d is not None


def test_bind_unbind_candidate_snapshot_rollback_commit_and_ledger_forwarding() -> None:
    population, _, _, _, kinetics, transaction, ledger_spy = _population()
    accepted = object()
    ledger = object()
    population.snapshot_candidate()
    population.propose_candidate(0.025, 71)
    population.rollback(accepted)
    population.commit_irreversible(accepted, 0.025, 71)
    population.accumulate_ledger(ledger)
    assert transaction.calls == [
        ("snapshot",),
        ("rollback", accepted),
        ("commit", accepted, 0.025, 71),
    ]
    assert kinetics.calls[0][1:] == (0.025, 71)
    assert ledger_spy.calls[0][1] is ledger


def test_native_segment_triangle_pair_closes_force_and_adjoint_work() -> None:
    filament = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    nucleus = np.array([[2.0, -1.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0]])
    df = np.array([[0.1, 0.2, 0.0], [-0.2, 0.1, 0.3]])
    dn = np.array([[0.0, -0.1, 0.2], [0.2, 0.0, -0.1], [0.1, 0.2, 0.0]])
    kwargs = {
        "rest_length_um": 0.5,
        "stiffness_pn_per_um": 3.0,
        "stiffening_per_um2": 0.4,
    }
    result = segment_triangle_linc_reference(
        filament,
        (0, 1),
        0.25,
        nucleus,
        (0, 1, 2),
        (0.2, 0.3, 0.5),
        filament_displacement=df,
        nucleus_displacement=dn,
        **kwargs,
    )
    assert result.force_residual_pn == pytest.approx(np.zeros(3), abs=1.0e-14)

    eps = 1.0e-7
    plus = segment_triangle_linc_reference(
        filament + eps * df,
        (0, 1),
        0.25,
        nucleus + eps * dn,
        (0, 1, 2),
        (0.2, 0.3, 0.5),
        **kwargs,
    )
    minus = segment_triangle_linc_reference(
        filament - eps * df,
        (0, 1),
        0.25,
        nucleus - eps * dn,
        (0, 1, 2),
        (0.2, 0.3, 0.5),
        **kwargs,
    )
    energy_derivative = (plus.energy_pn_um - minus.energy_pn_um) / (2.0 * eps)
    assert result.virtual_work_pn_um == pytest.approx(-energy_derivative, rel=1.0e-7)


def test_reduced_jt_projection_preserves_socket_virtual_work() -> None:
    jacobian = np.array(
        [
            [1.0, 0.0, 0.2, -0.1],
            [0.0, 1.0, 0.3, 0.4],
            [0.0, 0.0, 1.0, 0.5],
        ]
    )
    result = reduced_jt_projection_reference(
        jacobian,
        np.array([2.0, -3.0, 1.5]),
        np.array([0.2, -0.1, 0.4, 0.3]),
    )
    assert result.generalized_force == pytest.approx(jacobian.T @ np.array([2.0, -3.0, 1.5]))
    assert result.work_residual_pn_um == pytest.approx(0.0, abs=1.0e-15)


def test_force_work_and_population_ledger_gate() -> None:
    LincLedgerReport(8, 6, 4, 2, 1.0e-12, 2.0e-13).assert_closed(
        force_tolerance_pn=1.0e-10,
        work_tolerance_pn_um=1.0e-11,
    )
    with pytest.raises(ValueError, match="counts do not close"):
        LincLedgerReport(8, 6, 4, 1, 0.0, 0.0).assert_closed(
            force_tolerance_pn=0.0,
            work_tolerance_pn_um=0.0,
        )
    with pytest.raises(ValueError, match="force ledger"):
        LincLedgerReport(8, 6, 4, 2, 1.0, 0.0).assert_closed(
            force_tolerance_pn=0.5,
            work_tolerance_pn_um=0.0,
        )


def test_architecture_rejects_wrong_linc_family() -> None:
    architecture = reference_cell_architecture()
    connectors = tuple(
        replace(connector, family=ConnectorFamily.MOTOR)
        if connector.name == ACTIN_CAP_LINC
        else connector
        for connector in architecture.connectors
    )
    with pytest.raises(ValueError, match="ConnectorFamily.LINC"):
        validate_linc_architecture(CellArchitecture(architecture.components, connectors))


def test_bound_state_is_distinct_from_preallocated_unbound_record() -> None:
    assert int(LincJointState.UNBOUND) != int(LincJointState.BOUND)


# ---------------------------------------------------------------------------
# KERNEL_BOUND slice: real Warp delegates bound through the SEAMED protocols.
#
# These CPU-green tests capture every device launch through the module's single
# ``_linc_launch`` indirection, so they prove the seam binding (kernel identity +
# authoritative SoA argument order + device) without a CUDA device.  The final
# CUDA-gated test is the CUDA_UNIT numeric gate; it is skipped on CPU-only hosts.
# ---------------------------------------------------------------------------

_MUTABLE_NAMES = (
    "active_d",
    "joint_id_d",
    "joint_generation_d",
    "state_d",
    "candidate_state_d",
    "snapshot_state_d",
    "age_d",
    "snapshot_age_d",
    "rng_epoch_d",
    "snapshot_rng_epoch_d",
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
    "population_epoch_d",
    "snapshot_population_epoch_d",
)
_REMAP_NAMES = (
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


def _warp_delegate_population(
    edge_name: str = ACTIN_CAP_LINC,
) -> tuple[LincJointPopulation, dict[str, _FakeArray]]:
    """Build a population wired to the real Warp delegates over fake CUDA arrays."""
    capacity = 2
    ptr = iter(range(500, 700))

    def array(dtype: object, shape: tuple[int, ...] = (capacity,)) -> _FakeArray:
        return _FakeArray(next(ptr), shape, dtype)

    arrays = {
        "active_d": array(wp.int32),
        "joint_id_d": array(wp.int64),
        "joint_generation_d": array(wp.int32),
        "state_d": array(wp.int32),
        "candidate_state_d": array(wp.int32),
        "snapshot_state_d": array(wp.int32),
        "age_d": array(wp.float64),
        "snapshot_age_d": array(wp.float64),
        "rng_epoch_d": array(wp.int64),
        "snapshot_rng_epoch_d": array(wp.int64),
        "filament_actor_id_d": array(wp.int64),
        "filament_actor_generation_d": array(wp.int32),
        "filament_entity_id_d": array(wp.int64),
        "filament_entity_generation_d": array(wp.int32),
        "filament_element_id_d": array(wp.int32),
        "filament_u_d": array(wp.float64),
        "nucleus_actor_id_d": array(wp.int64),
        "nucleus_actor_generation_d": array(wp.int32),
        "nucleus_entity_id_d": array(wp.int64),
        "nucleus_entity_generation_d": array(wp.int32),
        "socket_face_id_d": array(wp.int32),
        "socket_barycentric_d": array(wp.vec3d),
        "rest_length_d": array(wp.float64),
        "stiffness_d": array(wp.float64),
        "stiffening_d": array(wp.float64),
        "load_d": array(wp.float64),
        "energy_d": array(wp.float64),
        "force_on_filament_d": array(wp.vec3d),
        "force_on_nucleus_d": array(wp.vec3d),
        "population_epoch_d": array(wp.int64, (1,)),
        "snapshot_population_epoch_d": array(wp.int64, (1,)),
    }
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
        remap_arrays=tuple(arrays[name] for name in _REMAP_NAMES),
    )
    population = LincJointPopulation(
        edge_name=edge_name,
        capacity=capacity,
        mechanics=WarpLincMechanics(),
        generation_validator=WarpLincGenerationValidator(
            stale_count_d=_FakeArray(701, (1,), wp.int32)
        ),
        projector=WarpNativeSocketProjector(),
        kinetics=SourceGatedLincKinetics(),
        transaction=transaction,
        ledger=WarpLincLedgerReducer(
            active_count_d=_FakeArray(702, (1,), wp.int32),
            bound_count_d=_FakeArray(703, (1,), wp.int32),
            energy_sum_d=_FakeArray(704, (1,), wp.float64),
            force_residual_d=_FakeArray(705, (1,), wp.vec3d),
        ),
        **arrays,
    )
    return population, arrays


def _record_launches(monkeypatch) -> list[dict[str, object]]:
    calls: list[dict[str, object]] = []

    def recorder(kernel: object, dim: int, inputs: list, outputs: list | None, device: str) -> None:
        calls.append(
            {
                "kernel": kernel,
                "dim": dim,
                "inputs": list(inputs),
                "outputs": None if outputs is None else list(outputs),
                "device": device,
            }
        )

    monkeypatch.setattr(linc_module, "_linc_launch", recorder)
    return calls


def test_warp_delegates_conform_to_seamed_protocols() -> None:
    assert isinstance(WarpLincMechanics(), LincMechanics)
    assert isinstance(WarpNativeSocketProjector(), NuclearSocketProjector)
    assert isinstance(SourceGatedLincKinetics(), LincKinetics)
    validator = WarpLincGenerationValidator(stale_count_d=_FakeArray(1, (1,), wp.int32))
    assert isinstance(validator, LincPortGenerationValidator)
    ledger = WarpLincLedgerReducer(
        active_count_d=_FakeArray(2, (1,), wp.int32),
        bound_count_d=_FakeArray(3, (1,), wp.int32),
        energy_sum_d=_FakeArray(4, (1,), wp.float64),
        force_residual_d=_FakeArray(5, (1,), wp.vec3d),
    )
    assert isinstance(ledger, LincLedgerContributor)


def test_accumulate_binds_generation_mechanics_and_native_scatter_kernels(monkeypatch) -> None:
    population, arrays = _warp_delegate_population()
    filament, nucleus = _geometry("sf_arc")
    calls = _record_launches(monkeypatch)
    population.accumulate(filament, nucleus)

    assert [call["kernel"] for call in calls] == [
        linc_module._linc_generation_validate_kernel,
        linc_module._linc_pair_force_kernel,
        linc_module._linc_native_socket_scatter_kernel,
    ]
    assert all(call["dim"] == population.capacity for call in calls)
    assert {call["device"] for call in calls} == {str(arrays["active_d"].device)}

    mechanics = calls[1]
    # Authoritative gather order: filament geometry, nucleus geometry, then the SoA.
    assert mechanics["inputs"][:5] == [
        filament.position_d,
        filament.force_d,
        filament.segments_d,
        nucleus.surface_position_d,
        nucleus.faces_d,
    ]
    assert mechanics["inputs"][9:11] == [arrays["socket_face_id_d"], arrays["socket_barycentric_d"]]
    assert mechanics["outputs"] == [
        arrays["load_d"],
        arrays["energy_d"],
        arrays["force_on_filament_d"],
        arrays["force_on_nucleus_d"],
    ]
    scatter = calls[2]
    assert scatter["inputs"][4] is arrays["force_on_nucleus_d"]
    assert scatter["inputs"][6] is nucleus.surface_force_d


def test_native_projector_rejects_reduced_mode_without_socket_jacobian(monkeypatch) -> None:
    population, _ = _warp_delegate_population()
    filament, reduced = _geometry("sf_arc", mode=NuclearProjectionMode.REDUCED_JT)
    _record_launches(monkeypatch)
    with pytest.raises(ValueError, match="reduced J.T projection is a Core-Jacobian seam gap"):
        population.accumulate(filament, reduced)


def test_transaction_commit_and_rollback_bind_state_machine_kernels(monkeypatch) -> None:
    population, arrays = _warp_delegate_population()
    calls = _record_launches(monkeypatch)
    accepted = object()
    population.rollback(accepted)
    population.commit_irreversible(accepted, 0.05, 11)

    assert calls[0]["kernel"] is linc_module._linc_rollback_kernel
    assert calls[0]["inputs"] == [accepted, arrays["state_d"], arrays["candidate_state_d"]]
    assert calls[1]["kernel"] is linc_module._linc_commit_kernel
    assert calls[1]["inputs"][0] is arrays["active_d"]
    assert calls[1]["inputs"][1] is accepted
    assert calls[1]["inputs"][2] is arrays["candidate_state_d"]
    assert calls[1]["inputs"][6] is arrays["population_epoch_d"]


def test_transaction_owned_arrays_cover_every_mutable_population_array() -> None:
    # Constructing the population already asserts coverage; verify the identity set too.
    population, arrays = _warp_delegate_population()
    covered = {linc_module._storage_key(a) for a in population.transaction.owned_arrays()}
    mutable = {linc_module._storage_key(arrays[name]) for name in _MUTABLE_NAMES}
    assert mutable <= covered


def test_commit_rejects_nonpositive_dt_and_negative_seed(monkeypatch) -> None:
    population, _ = _warp_delegate_population()
    _record_launches(monkeypatch)
    accepted = object()
    with pytest.raises(ValueError, match="dt_phys must be finite and positive"):
        population.commit_irreversible(accepted, 0.0, 1)
    with pytest.raises(ValueError, match="rng_seed must be nonnegative"):
        population.commit_irreversible(accepted, 0.05, -1)


def test_ledger_reducer_binds_per_edge_reduction_kernel(monkeypatch) -> None:
    population, arrays = _warp_delegate_population()
    calls = _record_launches(monkeypatch)
    population.accumulate_ledger(object())
    assert len(calls) == 1
    assert calls[0]["kernel"] is linc_module._linc_ledger_reduce_kernel
    assert calls[0]["inputs"][2] is arrays["energy_d"]
    assert calls[0]["inputs"][3] is arrays["force_on_filament_d"]
    assert calls[0]["inputs"][4] is arrays["force_on_nucleus_d"]


def test_source_gated_kinetics_declares_channels_but_refuses_to_fabricate_rates() -> None:
    kinetics = SourceGatedLincKinetics()
    assert LINC_EVENT_CHANNELS <= kinetics.event_channels
    population, _ = _warp_delegate_population()
    with pytest.raises(NotImplementedError, match="source-gated"):
        population.propose_candidate(0.025, 7)


def test_build_linc_joint_population_validates_specs_and_parameters() -> None:
    spec, _ = _spec(ACTIN_CAP_LINC)
    with pytest.raises(ValueError, match="every joint spec must belong"):
        build_linc_joint_population(
            MT_NUCLEUS_LINC,
            (spec,),
            rest_length_um=0.5,
            stiffness_pn_per_um=3.0,
            stiffening_per_um2=0.4,
            device="cuda:0",
        )
    with pytest.raises(ValueError, match="at least one joint spec"):
        build_linc_joint_population(
            ACTIN_CAP_LINC,
            (),
            rest_length_um=0.5,
            stiffness_pn_per_um=3.0,
            stiffening_per_um2=0.4,
            device="cuda:0",
        )
    with pytest.raises(ValueError, match="joint IDs must be unique"):
        build_linc_joint_population(
            ACTIN_CAP_LINC,
            (spec, spec),
            rest_length_um=0.5,
            stiffness_pn_per_um=3.0,
            stiffening_per_um2=0.4,
            device="cuda:0",
        )
    with pytest.raises(ValueError, match="stiffness must be positive"):
        build_linc_joint_population(
            ACTIN_CAP_LINC,
            (spec,),
            rest_length_um=0.5,
            stiffness_pn_per_um=0.0,
            stiffening_per_um2=0.4,
            device="cuda:0",
        )
    with pytest.raises(ValueError, match="must broadcast"):
        build_linc_joint_population(
            ACTIN_CAP_LINC,
            (spec,),
            rest_length_um=(0.5, 0.6),
            stiffness_pn_per_um=3.0,
            stiffening_per_um2=0.4,
            device="cuda:0",
        )


_CUDA_DEVICE = next((device for device in wp.get_devices() if device.is_cuda), None)


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="Warp LINC runtime requires CUDA")
def test_cuda_unit_linc_pair_force_projection_and_accepted_kinetics() -> None:
    """CUDA_UNIT gate: real kernels close force/energy vs the oracle and honour the
    accepted-step transaction on device.  Skipped on CPU-only hosts."""
    device = str(_CUDA_DEVICE)
    spec, _ = _spec(ACTIN_CAP_LINC)
    spec = replace(spec, initial_state=LincJointState.BOUND)
    rest, stiffness, stiffening = 0.5, 3.0, 0.4
    population = build_linc_joint_population(
        ACTIN_CAP_LINC,
        (spec,),
        rest_length_um=rest,
        stiffness_pn_per_um=stiffness,
        stiffening_per_um2=stiffening,
        device=device,
    )

    fil_pos = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    nuc_pos = np.array([[2.0, -1.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0]])
    fil_seg = np.array([[0, 1]], np.int32)
    nuc_face = np.array([[0, 1, 2]], np.int32)
    fil_record = ActorRecord("sf_arc", 10, 3, 100, 4, 1)
    nuc_record = ActorRecord("nucleus", 20, 5, 200, 6, 1)

    with wp.ScopedDevice(device):
        filament = FilamentGeometryView(
            fil_record,
            wp.array(fil_pos, dtype=wp.vec3d),
            wp.zeros(2, dtype=wp.vec3d),
            wp.array(fil_seg, dtype=wp.int32, ndim=2),
            wp.array(np.array([3], np.int32), dtype=wp.int32),
            wp.array(np.array([4], np.int32), dtype=wp.int32),
        )
        nucleus = NuclearSurfaceView(
            nuc_record,
            NuclearProjectionMode.NATIVE_SURFACE,
            wp.array(nuc_pos, dtype=wp.vec3d),
            wp.zeros(3, dtype=wp.vec3d),
            wp.array(nuc_face, dtype=wp.int32, ndim=2),
            wp.array(np.array([5], np.int32), dtype=wp.int32),
            wp.array(np.array([6], np.int32), dtype=wp.int32),
        )

    ref = segment_triangle_linc_reference(
        fil_pos,
        (0, 1),
        spec.filament.ref.segment_u,
        nuc_pos,
        (0, 1, 2),
        spec.nucleus.barycentric,
        rest_length_um=rest,
        stiffness_pn_per_um=stiffness,
        stiffening_per_um2=stiffening,
    )

    population.accumulate(filament, nucleus)
    wp.synchronize_device(device)
    np.testing.assert_allclose(filament.force_d.numpy(), ref.filament_force, atol=1.0e-12)
    np.testing.assert_allclose(nucleus.surface_force_d.numpy(), ref.nucleus_force, atol=1.0e-12)
    total = filament.force_d.numpy().sum(axis=0) + nucleus.surface_force_d.numpy().sum(axis=0)
    np.testing.assert_allclose(total, 0.0, atol=1.0e-12)
    assert population.energy_d.numpy()[0] == pytest.approx(ref.energy_pn_um)

    population.accumulate_ledger(object())
    wp.synchronize_device(device)
    assert population.ledger.active_count_d.numpy()[0] == 1
    assert population.ledger.bound_count_d.numpy()[0] == 1
    np.testing.assert_allclose(population.ledger.force_residual_d.numpy()[0], 0.0, atol=1.0e-12)

    # Reject: committed BOUND state, RNG epoch, and population epoch do not advance.
    population.candidate_state_d.assign(np.array([int(LincJointState.UNBOUND)], np.int32))
    rejected = wp.array(np.array([0], np.int32), dtype=wp.int32, device=device)
    population.commit_irreversible(rejected, 0.05, 3)
    wp.synchronize_device(device)
    assert population.state_d.numpy()[0] == int(LincJointState.BOUND)
    assert population.candidate_state_d.numpy()[0] == int(LincJointState.BOUND)
    assert population.rng_epoch_d.numpy()[0] == 0
    assert population.population_epoch_d.numpy()[0] == 0

    # Accept an unbind: state flips, epochs advance, age resets on transition.
    population.candidate_state_d.assign(np.array([int(LincJointState.UNBOUND)], np.int32))
    accepted = wp.array(np.array([1], np.int32), dtype=wp.int32, device=device)
    population.commit_irreversible(accepted, 0.05, 3)
    wp.synchronize_device(device)
    assert population.state_d.numpy()[0] == int(LincJointState.UNBOUND)
    assert population.rng_epoch_d.numpy()[0] == 1
    assert population.population_epoch_d.numpy()[0] == 1

    # An unbound tension-only tether carries no mechanical load.
    filament.force_d.zero_()
    nucleus.surface_force_d.zero_()
    population.accumulate(filament, nucleus)
    wp.synchronize_device(device)
    assert np.max(np.abs(filament.force_d.numpy())) == 0.0
    assert np.max(np.abs(nucleus.surface_force_d.numpy())) == 0.0
