"""Structural gates for the explicit head-resolved NMII actuator seam."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine.nmii_actuator import (
    CORTEX_COMPONENT,
    FILOPODIUM_COMPONENT,
    LAMELLIPODIUM_COMPONENT,
    NMII_COMPONENT,
    NMII_CONNECTOR_EVENT_CHANNELS,
    NMII_CORTEX_MOTOR,
    NMII_FILOPODIUM_MOTOR,
    NMII_LAMELLIPODIUM_MOTOR,
    NMII_LEDGER_CHANNELS,
    NMII_REPRESENTATION,
    NMII_SF_MOTOR,
    NMII_STATE_EVENT_CHANNELS,
    SF_COMPONENT,
    BackboneArmMechanics,
    FilamentMotorPortView,
    LegacyNMIIAdapter,
    LegacyNMIIKind,
    NMIIActuator,
    NMIIActuatorStateOwner,
    NMIIActuatorStepBindings,
    NMIIActuatorView,
    NMIIBendingStiffness,
    SegmentConnectorState,
    SegmentMotorConnectorRuntime,
    nmii_contract_proposal,
    proposed_nmii_architecture,
)
from aleph.engine.protrusion import proposed_protrusion_architecture
from aleph.components.motor.backbone_warp import (
    angle_harmonic_kernel,
    arm_orientation_k_theta,
    backbone_bending_k_theta,
    backbone_bending_kappa,
)
from aleph.components.motor.minifilament_warp import harmonic_bond_kernel
from aleph.components.motor.segment_motor import (
    _increment_epoch_if_accepted_kernel,
    _refresh_bound_walk_dir_kernel,
    attach_segment_gated_kernel,
    compute_head_loads_segment_split_kernel,
    crossbridge_segment_split_kernel,
    refresh_segment_barbed_kernel,
    step_detach_segment_gated_kernel,
)


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """CUDA metadata double; no host physics is evaluated by structural tests."""

    ptr: int
    shape: tuple[int, ...]
    dtype: object
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _ProductionMechanicsSpy:
    n_particles: int = 10
    n_heads: int = 4
    explicit_stam_hocky_backbone: bool = True
    individual_heads: bool = True
    per_head_bell: bool = True
    hill_force_velocity: bool = True
    graph_owned_motor_connectors: bool = True
    component_local: bool = True
    production_eligible: bool = True
    aggregate_or_two_anchor: bool = False
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate_internal(self, position: wp.array, force: wp.array) -> None:
        self.calls.append((position, force))


@dataclass(slots=True)
class _TransactionSpy:
    covered: tuple[object, ...]
    component_name: str = NMII_COMPONENT
    event_channels: frozenset[str] = NMII_STATE_EVENT_CHANNELS
    calls: list[tuple[object, ...]] = field(default_factory=list)

    def owned_arrays(self) -> tuple[wp.array, ...]:
        return self.covered  # type: ignore[return-value]

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))


@dataclass(slots=True)
class _LedgerSpy:
    ledger_channels: frozenset[str] = NMII_LEDGER_CHANNELS
    calls: list[object] = field(default_factory=list)

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append(ledger)


@dataclass(slots=True)
class _ConnectorSpy:
    name: str
    component_a: str
    component_b: str
    event_channels: frozenset[str] = NMII_CONNECTOR_EVENT_CHANNELS
    ledger_channels: frozenset[str] = NMII_LEDGER_CHANNELS
    individual_head_state: bool = True
    per_head_bell: bool = True
    hill_force_velocity: bool = True
    aggregate_or_two_anchor: bool = False
    mechanics_calls: list[tuple[NMIIActuatorView, FilamentMotorPortView]] = field(
        default_factory=list
    )
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_candidate(
        self,
        actuator: NMIIActuatorView,
        port: FilamentMotorPortView,
    ) -> None:
        self.mechanics_calls.append((actuator, port))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


def _state_arrays() -> dict[str, _FakeArray]:
    return {
        "position_d": _FakeArray(100, (10,), wp.vec3d),
        "force_d": _FakeArray(101, (10,), wp.vec3d),
        "minifilament_offset_d": _FakeArray(102, (3,), wp.int32),
        "active_minifilament_d": _FakeArray(103, (2,), wp.int32),
        "minifilament_id_d": _FakeArray(104, (2,), wp.int64),
        "particle_role_d": _FakeArray(105, (10,), wp.int32),
        "head_node_d": _FakeArray(106, (4,), wp.int32),
        "head_id_d": _FakeArray(107, (4,), wp.int64),
        "head_side_d": _FakeArray(108, (4,), wp.int32),
        "head_load_hill_d": _FakeArray(109, (4,), wp.float64),
        "head_load_bell_d": _FakeArray(110, (4,), wp.float64),
        "atp_cycle_state_d": _FakeArray(111, (4,), wp.int32),
        "atp_consumed_d": _FakeArray(112, (1,), wp.int64),
        "actuator_epoch_d": _FakeArray(113, (1,), wp.int32),
    }


def _state() -> tuple[
    NMIIActuatorStateOwner,
    _ProductionMechanicsSpy,
    _TransactionSpy,
    _LedgerSpy,
]:
    arrays = _state_arrays()
    mechanics = _ProductionMechanicsSpy()
    transaction = _TransactionSpy(
        covered=(
            arrays["position_d"],
            arrays["active_minifilament_d"],
            arrays["atp_cycle_state_d"],
            arrays["atp_consumed_d"],
            arrays["actuator_epoch_d"],
        )
    )
    ledger = _LedgerSpy()
    state = NMIIActuatorStateOwner(
        n_minifilaments=2,
        n_heads=4,
        mechanics=mechanics,
        transaction=transaction,
        ledger=ledger,
        **arrays,  # type: ignore[arg-type]
    )
    return state, mechanics, transaction, ledger


def _port(component: str, ptr: int) -> FilamentMotorPortView:
    return FilamentMotorPortView(
        component=component,
        position_d=_FakeArray(ptr, (8,), wp.vec3d),
        force_d=_FakeArray(ptr + 1, (8,), wp.vec3d),
        segment_node_a_d=_FakeArray(ptr + 2, (3,), wp.int32),
        segment_node_b_d=_FakeArray(ptr + 3, (3,), wp.int32),
        segment_polarity_d=_FakeArray(ptr + 4, (3,), wp.int32),
        persistent_filament_id_d=_FakeArray(ptr + 5, (3,), wp.int64),
        material_s0_d=_FakeArray(ptr + 6, (3,), wp.float64),
        material_s1_d=_FakeArray(ptr + 7, (3,), wp.float64),
        topology_epoch_d=_FakeArray(ptr + 8, (1,), wp.int32),
    )


def _bindings() -> tuple[NMIIActuatorStepBindings, dict[str, _ConnectorSpy]]:
    ports = (
        _port(SF_COMPONENT, 200),
        _port(CORTEX_COMPONENT, 300),
        _port(LAMELLIPODIUM_COMPONENT, 400),
        _port(FILOPODIUM_COMPONENT, 500),
    )
    connectors = {
        NMII_SF_MOTOR: _ConnectorSpy(NMII_SF_MOTOR, NMII_COMPONENT, SF_COMPONENT),
        NMII_CORTEX_MOTOR: _ConnectorSpy(NMII_CORTEX_MOTOR, NMII_COMPONENT, CORTEX_COMPONENT),
        NMII_LAMELLIPODIUM_MOTOR: _ConnectorSpy(
            NMII_LAMELLIPODIUM_MOTOR,
            NMII_COMPONENT,
            LAMELLIPODIUM_COMPONENT,
        ),
        NMII_FILOPODIUM_MOTOR: _ConnectorSpy(
            NMII_FILOPODIUM_MOTOR,
            NMII_COMPONENT,
            FILOPODIUM_COMPONENT,
        ),
    }
    return NMIIActuatorStepBindings(ports, tuple(connectors.values())), connectors


def _architecture() -> CellArchitecture:
    protrusions = proposed_protrusion_architecture(reference_cell_architecture())
    return proposed_nmii_architecture(protrusions)


def _architecture_replacing_connector(replacement: ConnectorContract) -> CellArchitecture:
    architecture = _architecture()
    return CellArchitecture(
        architecture.components,
        tuple(
            replacement if connector.name == replacement.name else connector
            for connector in architecture.connectors
        ),
    )


def test_common_graph_promotes_nmii_proposal_idempotently() -> None:
    base = reference_cell_architecture()
    proposal = nmii_contract_proposal()
    common_components = {component.name for component in base.components}
    common_connectors = {connector.name for connector in base.connectors}

    assert NMII_COMPONENT in common_components
    assert {connector.name for connector in proposal.connectors} <= common_connectors
    assert {LAMELLIPODIUM_COMPONENT, FILOPODIUM_COMPONENT} <= common_components
    promoted = proposed_nmii_architecture(base)
    assert promoted.component(NMII_COMPONENT).representation == NMII_REPRESENTATION
    assert len(promoted.components) == len(base.components)
    assert len(promoted.connectors) == len(base.connectors)


def test_proposed_graph_registers_explicit_nmii_and_four_motor_edges() -> None:
    architecture = _architecture()
    component = architecture.component(NMII_COMPONENT)
    proposal = nmii_contract_proposal()
    contracts = {connector.name: connector for connector in architecture.connectors}

    assert component.role is ComponentRole.ACTIVE_LOAD_PATH
    assert component.representation == NMII_REPRESENTATION
    assert component.owns_geometry and component.dynamically_evolving
    for proposed in proposal.connectors:
        actual = contracts[proposed.name]
        assert actual.family is ConnectorFamily.MOTOR
        assert actual.kinetics and actual.commit_on_accept
        assert actual.bidirectional and actual.adjoint_transfer_required
        assert actual.generation_required and actual.remap_on_accept
        assert actual.blocks_sleep_refine


def test_architecture_rejects_aggregate_representation_and_wrong_motor_family() -> None:
    architecture = _architecture()
    wrong_components = tuple(
        replace(component, representation="two-anchor aggregate force dipole")
        if component.name == NMII_COMPONENT
        else component
        for component in architecture.components
    )
    state, *_ = _state()
    with pytest.raises(ValueError, match="aggregate or two-anchor"):
        NMIIActuator(CellArchitecture(wrong_components, architecture.connectors), state)

    wrong_motor = ConnectorContract(
        NMII_SF_MOTOR,
        ConnectorFamily.TRANSIENT_ACTIN,
        NMII_COMPONENT,
        SF_COMPONENT,
        True,
        True,
    )
    with pytest.raises(ValueError, match="ConnectorFamily.MOTOR"):
        NMIIActuator(_architecture_replacing_connector(wrong_motor), state)


def test_state_exposes_explicit_minifilament_and_individual_head_identity() -> None:
    state, *_ = _state()
    view = state.geometry()

    assert view.position_d is state.position_d
    assert view.minifilament_id_d is state.minifilament_id_d
    assert view.head_node_d is state.head_node_d
    assert view.head_id_d is state.head_id_d
    assert view.head_side_d is state.head_side_d
    assert view.head_load_hill_d is state.head_load_hill_d
    assert view.head_load_bell_d is state.head_load_bell_d
    assert view.n_minifilaments == 2 and view.n_heads == 4
    assert not hasattr(view, "bound_d")
    assert not hasattr(view, "abscissa_d")


def test_state_rejects_cpu_cross_device_and_storage_aliasing() -> None:
    state, *_ = _state()
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        replace(state, position_d=_FakeArray(900, (10,), wp.vec3d, cpu))

    other_gpu = _FakeDevice(alias="cuda:1")
    with pytest.raises(ValueError, match="share one CUDA device"):
        replace(state, force_d=_FakeArray(901, (10,), wp.vec3d, other_gpu))

    with pytest.raises(ValueError, match="distinct storage"):
        replace(
            state,
            head_load_bell_d=_FakeArray(state.head_load_hill_d.ptr, (4,), wp.float64),
        )


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"individual_heads": False}, "fidelity guarantees"),
        ({"graph_owned_motor_connectors": False}, "fidelity guarantees"),
        ({"aggregate_or_two_anchor": True}, "diagnostic-only"),
        ({"n_heads": 3}, "head count"),
    ],
)
def test_production_mechanics_cannot_drop_head_resolved_hill_bell_graph_contract(
    change: dict[str, object],
    message: str,
) -> None:
    state, mechanics, *_ = _state()
    with pytest.raises(ValueError, match=message):
        replace(state, mechanics=replace(mechanics, **change))


def test_state_transaction_and_ledger_require_all_motor_channels() -> None:
    state, _, transaction, ledger = _state()
    with pytest.raises(ValueError, match="missing event channels"):
        replace(
            state,
            transaction=replace(
                transaction,
                event_channels=frozenset({"atp_cycle"}),
            ),
        )
    with pytest.raises(ValueError, match="transaction must cover"):
        replace(state, transaction=replace(transaction, covered=transaction.covered[:-1]))
    with pytest.raises(ValueError, match="missing channels"):
        replace(
            state,
            ledger=replace(ledger, ledger_channels=frozenset({"force", "work", "population"})),
        )


def test_filament_port_validates_persistent_segment_material_mapping() -> None:
    port = _port(SF_COMPONENT, 200)
    assert port.persistent_filament_id_d.shape == port.segment_node_a_d.shape
    assert port.material_s0_d.shape == port.material_s1_d.shape

    with pytest.raises(ValueError, match="unsupported"):
        replace(port, component="ecm")
    with pytest.raises(ValueError, match="length 3"):
        replace(port, material_s1_d=_FakeArray(999, (2,), wp.float64))
    with pytest.raises(ValueError, match="distinct storage"):
        replace(
            port,
            material_s1_d=_FakeArray(port.material_s0_d.ptr, (3,), wp.float64),
        )


def test_bindings_require_exact_ports_connector_events_and_ledger_channels() -> None:
    bindings, connectors = _bindings()
    with pytest.raises(ValueError, match="filament ports mismatch"):
        replace(bindings, ports=bindings.ports[:-1])
    with pytest.raises(ValueError, match="MOTOR bindings mismatch"):
        replace(bindings, connectors=bindings.connectors[:-1])
    with pytest.raises(ValueError, match="incorrect component endpoints"):
        replace(
            bindings,
            connectors=tuple(
                replace(connector, component_b=CORTEX_COMPONENT)
                if connector.name == NMII_SF_MOTOR
                else connector
                for connector in bindings.connectors
            ),
        )
    with pytest.raises(ValueError, match="missing events"):
        replace(
            bindings,
            connectors=tuple(
                replace(connector, event_channels=frozenset({"stepping"}))
                if connector.name == NMII_CORTEX_MOTOR
                else connector
                for connector in bindings.connectors
            ),
        )
    with pytest.raises(ValueError, match="missing per-head fidelity"):
        replace(
            bindings,
            connectors=tuple(
                replace(connector, per_head_bell=False)
                if connector.name == NMII_SF_MOTOR
                else connector
                for connector in bindings.connectors
            ),
        )
    with pytest.raises(ValueError, match="aggregate/two-anchor"):
        replace(
            bindings,
            connectors=tuple(
                replace(connector, aggregate_or_two_anchor=True)
                if connector.name == NMII_CORTEX_MOTOR
                else connector
                for connector in bindings.connectors
            ),
        )
    with pytest.raises(ValueError, match="missing channels"):
        replace(
            bindings,
            connectors=tuple(
                replace(connector, ledger_channels=frozenset({"force", "work"}))
                if connector is connectors[NMII_FILOPODIUM_MOTOR]
                else connector
                for connector in bindings.connectors
            ),
        )


def test_candidate_mechanics_wires_internal_state_and_each_external_port() -> None:
    state, mechanics, *_ = _state()
    actuator = NMIIActuator(_architecture(), state)
    bindings, connectors = _bindings()
    actuator.accumulate_candidate(bindings)

    assert mechanics.calls == [(state.position_d, state.force_d)]
    for name, connector in connectors.items():
        assert len(connector.mechanics_calls) == 1
        view, port = connector.mechanics_calls[0]
        assert view.head_id_d is state.head_id_d
        assert view.head_load_hill_d is state.head_load_hill_d
        expected_target = {
            NMII_SF_MOTOR: SF_COMPONENT,
            NMII_CORTEX_MOTOR: CORTEX_COMPONENT,
            NMII_LAMELLIPODIUM_MOTOR: LAMELLIPODIUM_COMPONENT,
            NMII_FILOPODIUM_MOTOR: FILOPODIUM_COMPONENT,
        }[name]
        assert port.component == expected_target


def test_transaction_and_four_channel_ledger_cover_state_and_all_ports() -> None:
    state, _, transaction, state_ledger = _state()
    actuator = NMIIActuator(_architecture(), state)
    bindings, connectors = _bindings()
    accepted_d = object()
    ledger = object()

    actuator.snapshot_candidate(bindings)
    actuator.rollback(bindings, accepted_d)
    actuator.commit_irreversible(bindings, accepted_d, dt_phys=0.01, rng_seed=131)
    actuator.accumulate_ledger(bindings, ledger)

    expected = [
        ("snapshot",),
        ("rollback", accepted_d),
        ("commit", accepted_d, 0.01, 131),
    ]
    assert transaction.calls == expected
    assert state_ledger.calls == [ledger]
    for connector in connectors.values():
        assert connector.transaction_calls == expected
        assert connector.ledger_calls == [ledger]


@dataclass(slots=True)
class _LegacyHeadResolvedSpy:
    backbone_bonds: _FakeArray = _FakeArray(600, (4,), wp.vec2i)
    head_bonds: _FakeArray = _FakeArray(601, (4,), wp.vec2i)
    head_node: _FakeArray = _FakeArray(602, (4,), wp.int32)
    loads: _FakeArray = _FakeArray(603, (4,), wp.float64)
    loads_full: _FakeArray = _FakeArray(604, (4,), wp.float64)
    params: object = field(default_factory=object)
    segment_runtime: object | None = field(default_factory=object)
    state: dict[str, _FakeArray] = field(
        default_factory=lambda: {
            "bound": _FakeArray(605, (4,), wp.int32),
            "seg_id": _FakeArray(606, (4,), wp.int32),
            "seg_a": _FakeArray(607, (4,), wp.int32),
            "seg_b": _FakeArray(608, (4,), wp.int32),
            "bary_t": _FakeArray(609, (4,), wp.float64),
            "abscissa": _FakeArray(610, (4,), wp.float64),
            "walk_dir": _FakeArray(611, (4,), wp.vec3d),
        }
    )
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


@dataclass(slots=True)
class _AggregateDiagnosticSpy:
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


def test_head_resolved_legacy_adapter_reuses_current_runtime_but_is_not_production() -> None:
    legacy = _LegacyHeadResolvedSpy()
    adapter = LegacyNMIIAdapter(legacy, LegacyNMIIKind.HEAD_RESOLVED_MYOSIN_FORCE)
    pos = _FakeArray(700, (10,), wp.vec3d)
    force = _FakeArray(701, (10,), wp.vec3d)
    adapter.accumulate_reference(pos, force)

    assert legacy.calls == [(pos, force)]
    assert adapter.status == "HEAD_RESOLVED_MONOLITHIC_CONNECTOR_SPLIT_PENDING"
    assert not adapter.component_local
    assert not adapter.graph_owned_motor_connectors
    assert not adapter.production_eligible
    assert not adapter.aggregate_or_two_anchor

    state, *_ = _state()
    with pytest.raises(ValueError, match="fidelity guarantees"):
        replace(state, mechanics=adapter)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "kind",
    [LegacyNMIIKind.AGGREGATE_FORCE_DIPOLE, LegacyNMIIKind.TWO_ANCHOR_LINEAR],
)
def test_aggregate_and_two_anchor_adapters_are_diagnostic_only(kind: LegacyNMIIKind) -> None:
    adapter = LegacyNMIIAdapter(_AggregateDiagnosticSpy(), kind)  # type: ignore[arg-type]
    assert adapter.status.startswith("DIAGNOSTIC_ONLY_")
    assert adapter.aggregate_or_two_anchor
    assert not adapter.production_eligible

    state, *_ = _state()
    with pytest.raises(ValueError, match="fidelity guarantees"):
        replace(state, mechanics=adapter)  # type: ignore[arg-type]


def test_legacy_head_resolved_label_cannot_hide_missing_segment_state() -> None:
    with pytest.raises(ValueError, match="production segment anchoring"):
        LegacyNMIIAdapter(
            _LegacyHeadResolvedSpy(segment_runtime=None),
            LegacyNMIIKind.HEAD_RESOLVED_MYOSIN_FORCE,
        )

    malformed = _LegacyHeadResolvedSpy()
    malformed.state.pop("walk_dir")
    with pytest.raises(ValueError, match="individual segment-bound head state"):
        LegacyNMIIAdapter(malformed, LegacyNMIIKind.HEAD_RESOLVED_MYOSIN_FORCE)


# ── KERNEL_BOUND: real ac/motor kernels bound behind the seam (CPU-green via recording doubles) ──────


@dataclass(slots=True)
class _RecordingLauncher:
    """CUDA-free launcher double: records which real Kernel object flows through the seam with what arrays."""

    calls: list[tuple] = field(default_factory=list)

    def __call__(self, kernel, *, dim, inputs, outputs=None, device=None) -> None:
        self.calls.append((kernel, int(dim), tuple(inputs), tuple(outputs or ()), device))

    @property
    def kernels(self) -> list:
        return [call[0] for call in self.calls]


@dataclass(slots=True)
class _RecordingCopy:
    calls: list[tuple] = field(default_factory=list)

    def __call__(self, dst, src) -> None:
        self.calls.append((dst, src))


def _bending() -> NMIIBendingStiffness:
    # DERIVED inputs (KB/PI GAP-flagged), never fit to a gate; k_xb=1000 pN/µm is the ac/ master knob band.
    return NMIIBendingStiffness(
        persistence_length_um=1.0,
        segment_len_um=0.1,
        k_xb_pn_per_um=1000.0,
        r0_head_um=0.02,
    )


def _topology_arrays() -> dict[str, _FakeArray]:
    return {
        "backbone_bonds_d": _FakeArray(120, (2,), wp.int32),
        "head_bonds_d": _FakeArray(121, (4,), wp.int32),
        "backbone_angles_d": _FakeArray(122, (1,), wp.int32),
        "head_arm_angles_d": _FakeArray(123, (4,), wp.int32),
    }


def _real_mechanics(launch: _RecordingLauncher) -> BackboneArmMechanics:
    return BackboneArmMechanics(
        n_particles=10,
        n_heads=4,
        k_backbone_pn_per_um=1.0e4,
        r0_backbone_um=0.05,
        k_head_spring_pn_per_um=1.0e3,
        r0_head_um=0.02,
        bending=_bending(),
        launch=launch,
        **_topology_arrays(),  # type: ignore[arg-type]
    )


def test_bending_stiffness_is_derived_not_a_magic_number() -> None:
    bending = _bending()
    kappa = backbone_bending_kappa(1.0)
    assert bending.k_theta_backbone == backbone_bending_k_theta(kappa, 0.1)
    assert bending.k_theta_arm == arm_orientation_k_theta(1000.0, 0.02)
    assert not bending.sourced  # coiled-coil L_p proxy stays a GAP until PI/KB ratifies it


def test_backbone_arm_mechanics_binds_real_kernels_without_crossbridge() -> None:
    launch = _RecordingLauncher()
    mechanics = _real_mechanics(launch)
    position = _FakeArray(700, (10,), wp.vec3d)
    force = _FakeArray(701, (10,), wp.vec3d)

    mechanics.accumulate_internal(position, force)

    # exactly the four real internal kernels — backbone rod, head-arm rod, backbone bending, arm orientation.
    assert launch.kernels == [
        harmonic_bond_kernel,
        harmonic_bond_kernel,
        angle_harmonic_kernel,
        angle_harmonic_kernel,
    ]
    # ownership split: NO head-to-actin crossbridge kernel is launched by the internal mechanics.
    assert mechanics.binds_crossbridge is False
    # no private state path: every launch writes the injected authoritative force and reads injected position.
    topo = _topology_arrays()
    for kernel, dim, inputs, outputs, _device in launch.calls:
        assert outputs == (force,)
        assert inputs[0] is position
    # dims track the real component-owned topology array lengths, not a hidden reshape.
    assert launch.calls[0][1] == topo["backbone_bonds_d"].shape[0]
    assert launch.calls[3][1] == topo["head_arm_angles_d"].shape[0]


def test_backbone_arm_mechanics_installs_as_production_state_backend() -> None:
    launch = _RecordingLauncher()
    mechanics = _real_mechanics(launch)
    state = replace(_state()[0], mechanics=mechanics)  # passes _validate_production_mechanics
    actuator = NMIIActuator(_architecture(), state)

    actuator.state.accumulate_internal()
    assert launch.kernels[0] is harmonic_bond_kernel
    assert launch.calls[0][2][0] is state.position_d  # actuator's own particle array, not a copy


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"individual_heads": False}, "fidelity guarantees"),
        ({"n_heads": 7}, "head count"),
    ],
)
def test_real_mechanics_cannot_be_installed_if_fidelity_or_counts_drift(
    change: dict[str, object],
    message: str,
) -> None:
    launch = _RecordingLauncher()
    mechanics = _real_mechanics(launch)
    for name, value in change.items():
        setattr(mechanics, name, value)
    with pytest.raises(ValueError, match=message):
        replace(_state()[0], mechanics=mechanics)


def _connector_state() -> SegmentConnectorState:
    def fa(ptr: int, n: int, dtype: object) -> _FakeArray:
        return _FakeArray(ptr, (n,), dtype)

    return SegmentConnectorState(
        n_heads=4,
        n_segments=3,
        bound_d=fa(800, 4, wp.int32),
        seg_id_d=fa(801, 4, wp.int32),
        seg_a_d=fa(802, 4, wp.int32),
        seg_b_d=fa(803, 4, wp.int32),
        bary_t_d=fa(804, 4, wp.float64),
        abscissa_d=fa(805, 4, wp.float64),
        walk_dir_d=fa(806, 4, wp.vec3d),
        rng_epoch_d=fa(807, 1, wp.int32),
        loads_hill_d=fa(808, 4, wp.float64),
        loads_bell_d=fa(809, 4, wp.float64),
        seg_barbed_d=fa(810, 3, wp.vec3d),
        query_seg_id_d=fa(811, 4, wp.int32),
        query_t_d=fa(812, 4, wp.float64),
        query_barbed_d=fa(813, 4, wp.vec3d),
        bound_snap_d=fa(820, 4, wp.int32),
        seg_id_snap_d=fa(821, 4, wp.int32),
        seg_a_snap_d=fa(822, 4, wp.int32),
        seg_b_snap_d=fa(823, 4, wp.int32),
        bary_t_snap_d=fa(824, 4, wp.float64),
        abscissa_snap_d=fa(825, 4, wp.float64),
        walk_dir_snap_d=fa(826, 4, wp.vec3d),
        rng_epoch_snap_d=fa(827, 1, wp.int32),
    )


class _XbParams:
    """Minimal crossbridge scalar carrier (k_xb/r0_xb) the split kernels read as launch inputs."""

    k_xb = wp.float64(1000.0)   # ac/ master knob band (NOT fit to a gate)
    r0_xb = wp.float64(0.0)


def _real_connector(
    name: str,
    target: str,
    launch: _RecordingLauncher,
    copy: _RecordingCopy,
) -> SegmentMotorConnectorRuntime:
    connector = SegmentMotorConnectorRuntime(
        name=name,
        component_b=target,
        state=_connector_state(),
        head_node_d=_FakeArray(830, (4,), wp.int32),
        params=_XbParams(),
        launch=launch,
        copy=copy,
    )
    connector.bind_target_topology(_FakeArray(840, (3,), wp.int32), _FakeArray(841, (3,), wp.int32))
    return connector


def test_segment_connector_accumulate_binds_live_geometry_kernels() -> None:
    launch = _RecordingLauncher()
    connector = _real_connector(NMII_SF_MOTOR, SF_COMPONENT, launch, _RecordingCopy())
    state, *_ = _state()
    view = state.geometry()
    port = _port(SF_COMPONENT, 200)

    connector.accumulate_candidate(view, port)

    assert launch.kernels == [
        refresh_segment_barbed_kernel,
        _refresh_bound_walk_dir_kernel,
        crossbridge_segment_split_kernel,
    ]
    # geometry is gathered from the TARGET-owned port, walk direction written to CONNECTOR-owned state.
    barbed_inputs = launch.calls[0][2]
    assert barbed_inputs[0] is port.position_d
    assert barbed_inputs[-1] is connector.state.seg_barbed_d
    walk_inputs = launch.calls[1][2]
    assert walk_inputs[0] is port.position_d
    assert walk_inputs[-1] is connector.state.walk_dir_d
    # the two-array adjoint crossbridge: head force into the ACTUATOR array, reactions into the PORT array,
    # with no shared/merged pos or force buffer (the split-ownership Newton-3rd load path).
    xb_inputs = launch.calls[2][2]
    assert xb_inputs[0] is view.position_d      # actuator (nmii-owned) position
    assert xb_inputs[1] is view.force_d         # actuator (nmii-owned) force  ← +f
    assert xb_inputs[2] is port.position_d      # target-owned segment position
    assert xb_inputs[3] is port.force_d         # target-owned segment force   ← −(1−t)f, −t f
    assert xb_inputs[4] is connector.head_node_d
    assert view.force_d is not port.force_d     # the two force arrays are never merged


def test_segment_connector_compute_loads_binds_split_load_kernel() -> None:
    launch = _RecordingLauncher()
    connector = _real_connector(NMII_SF_MOTOR, SF_COMPONENT, launch, _RecordingCopy())
    state, *_ = _state()
    view = state.geometry()
    port = _port(SF_COMPONENT, 200)

    connector.compute_loads(view, port)

    assert launch.kernels == [compute_head_loads_segment_split_kernel]
    load_inputs = launch.calls[0][2]
    assert load_inputs[0] is view.position_d    # head node read from the actuator array
    assert load_inputs[1] is port.position_d    # segment nodes read from the target port array
    # Hill (tangential/walk-dir) and Bell (full-|F|) per-head loads land in CONNECTOR-owned buffers.
    assert load_inputs[-2] is connector.state.loads_hill_d
    assert load_inputs[-1] is connector.state.loads_bell_d


def test_segment_connector_commit_binds_accepted_predicated_kmc_kernels() -> None:
    launch = _RecordingLauncher()
    connector = _real_connector(NMII_CORTEX_MOTOR, CORTEX_COMPONENT, launch, _RecordingCopy())
    accepted = _FakeArray(900, (1,), wp.int32)

    connector.commit_irreversible(accepted, dt_phys=0.01, rng_seed=131)

    assert launch.kernels == [
        attach_segment_gated_kernel,
        step_detach_segment_gated_kernel,
        _increment_epoch_if_accepted_kernel,
    ]
    # every KMC kernel is predicated on the scheduler-owned accepted scalar (first input).
    for _kernel, _dim, inputs, _outputs, _device in launch.calls:
        assert inputs[0] is accepted


def test_segment_connector_snapshot_and_rollback_cover_authoritative_state() -> None:
    copy = _RecordingCopy()
    connector = _real_connector(NMII_SF_MOTOR, SF_COMPONENT, _RecordingLauncher(), copy)
    pairs = connector.state.authoritative_pairs()

    connector.snapshot_candidate()
    assert copy.calls == [(snap, live) for live, snap in pairs]  # snapshot: snap <- live

    copy.calls.clear()
    connector.rollback(_FakeArray(900, (1,), wp.int32))
    assert copy.calls == [(live, snap) for live, snap in pairs]  # rollback: live <- snap


def test_segment_connector_rejects_wrong_target_and_unknown_name() -> None:
    with pytest.raises(ValueError, match="must target"):
        SegmentMotorConnectorRuntime(
            name=NMII_SF_MOTOR,
            component_b=CORTEX_COMPONENT,
            state=_connector_state(),
            head_node_d=_FakeArray(830, (4,), wp.int32),
            params=object(),
        )
    with pytest.raises(ValueError, match="unknown NMII MOTOR connector"):
        SegmentMotorConnectorRuntime(
            name="nmii_membrane_motor",
            component_b="membrane",
            state=_connector_state(),
            head_node_d=_FakeArray(830, (4,), wp.int32),
            params=object(),
        )


def test_real_connectors_drive_the_actuator_through_step_bindings() -> None:
    launch = _RecordingLauncher()
    copy = _RecordingCopy()
    ports = (
        _port(SF_COMPONENT, 200),
        _port(CORTEX_COMPONENT, 300),
        _port(LAMELLIPODIUM_COMPONENT, 400),
        _port(FILOPODIUM_COMPONENT, 500),
    )
    connectors = (
        _real_connector(NMII_SF_MOTOR, SF_COMPONENT, launch, copy),
        _real_connector(NMII_CORTEX_MOTOR, CORTEX_COMPONENT, launch, copy),
        _real_connector(NMII_LAMELLIPODIUM_MOTOR, LAMELLIPODIUM_COMPONENT, launch, copy),
        _real_connector(NMII_FILOPODIUM_MOTOR, FILOPODIUM_COMPONENT, launch, copy),
    )
    bindings = NMIIActuatorStepBindings(ports, connectors)  # passes fidelity/endpoint/ledger validation
    state = replace(_state()[0], mechanics=_real_mechanics(launch))
    actuator = NMIIActuator(_architecture(), state)

    actuator.accumulate_candidate(bindings)

    # four internal kernels + two live-geometry kernels per port × four ports.
    assert launch.kernels.count(harmonic_bond_kernel) == 2
    assert launch.kernels.count(refresh_segment_barbed_kernel) == 4
    assert launch.kernels.count(_refresh_bound_walk_dir_kernel) == 4
