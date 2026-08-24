"""Structural gates for the dynamic microtubule rig ownership seam."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorContract,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine.microtubule_rig import (
    MT_CORTEX_CAPTURE,
    MT_CYTOSOL_TRANSFER,
    MT_NUCLEUS_LINC,
    MT_SF_SPECTRAPLAKIN,
    STATIC_ADAPTER_STATUS,
    DynamicInstabilityRateCard,
    ExternalMechanicalEndpoint,
    LegacyMicrotubuleBendingAdapter,
    MicrotubulePhase,
    MicrotubuleRig,
    MicrotubuleRigStateOwner,
    MicrotubuleRigView,
    MicrotubuleRodBuildSpec,
    MicrotubuleRodRuntime,
    MicrotubuleStepBindings,
    dynamic_instability_reference,
    dynamic_instability_stochastic_reference,
    mtoc_anchor_reference,
)
from aleph.engine.runtime import CytosolFieldEndpoint


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """CUDA metadata double; tests never execute host physics with it."""

    ptr: int
    shape: tuple[int, ...] = (8,)
    dtype: object = wp.vec3d
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _MechanicsSpy:
    n_nodes: int = 8
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


@dataclass(slots=True)
class _TransactionSpy:
    covered: tuple[object, ...]
    component_name: str = "microtubule"
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
    calls: list[object] = field(default_factory=list)

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append(ledger)


@dataclass(slots=True)
class _GraphConnectorSpy:
    name: str
    component_a: str
    component_b: str
    mechanics_calls: list[tuple[MicrotubuleRigView, ExternalMechanicalEndpoint]] = field(
        default_factory=list
    )
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_rig(
        self,
        rig: MicrotubuleRigView,
        endpoint: ExternalMechanicalEndpoint,
    ) -> None:
        self.mechanics_calls.append((rig, endpoint))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _TransferSpy(_GraphConnectorSpy):
    def accumulate_transfer(self, rig: object, cytosol: object) -> None:
        self.mechanics_calls.append((rig, cytosol))  # type: ignore[arg-type]


@dataclass(slots=True)
class _MotorPortSpy:
    name: str = "kinesin_cargo_port"
    component_a: str = "microtubule"
    component_b: str = "nucleus"
    mechanics_calls: list[MicrotubuleRigView] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_motor(self, rig: MicrotubuleRigView) -> None:
        self.mechanics_calls.append(rig)

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
        "position_d": _FakeArray(100),
        "force_d": _FakeArray(101),
        "mtoc_position_d": _FakeArray(102, shape=(1,)),
        "mtoc_force_d": _FakeArray(103, shape=(1,)),
        "fiber_offset_d": _FakeArray(104, shape=(3,), dtype=wp.int32),
        "active_count_d": _FakeArray(105, shape=(2,), dtype=wp.int32),
        "plus_end_node_d": _FakeArray(106, shape=(2,), dtype=wp.int32),
        "phase_d": _FakeArray(107, shape=(2,), dtype=wp.int32),
        "topology_epoch_d": _FakeArray(108, shape=(1,), dtype=wp.int32),
    }


def _state() -> tuple[
    MicrotubuleRigStateOwner,
    _MechanicsSpy,
    _TransactionSpy,
    _LedgerSpy,
]:
    arrays = _state_arrays()
    mechanics = _MechanicsSpy()
    transaction = _TransactionSpy(
        covered=(
            arrays["position_d"],
            arrays["mtoc_position_d"],
            arrays["active_count_d"],
            arrays["plus_end_node_d"],
            arrays["phase_d"],
            arrays["topology_epoch_d"],
        )
    )
    ledger = _LedgerSpy()
    state = MicrotubuleRigStateOwner(
        n_mt=2,
        mechanics=mechanics,
        transaction=transaction,
        ledger=ledger,
        **arrays,  # type: ignore[arg-type]
    )
    return state, mechanics, transaction, ledger


def _bindings() -> tuple[
    MicrotubuleStepBindings,
    _GraphConnectorSpy,
    _GraphConnectorSpy,
    _GraphConnectorSpy,
    _TransferSpy,
    _MotorPortSpy,
]:
    cortex = ExternalMechanicalEndpoint("cortex", _FakeArray(200), _FakeArray(201))
    nucleus = ExternalMechanicalEndpoint("nucleus", _FakeArray(300), _FakeArray(301))
    sf_arc = ExternalMechanicalEndpoint("sf_arc", _FakeArray(310), _FakeArray(311))
    cytosol = CytosolFieldEndpoint(
        "cytosol",
        _FakeArray(320, shape=(2, 2, 2), dtype=wp.float64),
        _FakeArray(321, shape=(2, 2, 2), dtype=wp.float64),
    )
    capture = _GraphConnectorSpy(MT_CORTEX_CAPTURE, "microtubule", "cortex")
    linc = _GraphConnectorSpy(MT_NUCLEUS_LINC, "microtubule", "nucleus")
    spectraplakin = _GraphConnectorSpy(MT_SF_SPECTRAPLAKIN, "microtubule", "sf_arc")
    transfer = _TransferSpy(MT_CYTOSOL_TRANSFER, "microtubule", "cytosol")
    motor = _MotorPortSpy()
    return (
        MicrotubuleStepBindings(
            cortex,
            nucleus,
            sf_arc,
            cytosol,
            capture,
            linc,
            spectraplakin,
            transfer,
            (motor,),
        ),
        capture,
        linc,
        spectraplakin,
        transfer,
        motor,
    )


def _architecture_without(connector_name: str) -> CellArchitecture:
    base = reference_cell_architecture()
    return CellArchitecture(
        base.components,
        tuple(connector for connector in base.connectors if connector.name != connector_name),
    )


def _architecture_replacing(replacement: ConnectorContract) -> CellArchitecture:
    base = reference_cell_architecture()
    return CellArchitecture(
        base.components,
        tuple(
            replacement if connector.name == replacement.name else connector
            for connector in base.connectors
        ),
    )


def test_reference_architecture_registers_kinetic_cortex_capture() -> None:
    architecture = reference_cell_architecture()
    state, *_ = _state()
    rig = MicrotubuleRig(architecture, state)
    capture = next(
        connector for connector in architecture.connectors if connector.name == MT_CORTEX_CAPTURE
    )

    assert rig.architecture is architecture
    assert capture.family is ConnectorFamily.MOTOR
    assert {capture.component_a, capture.component_b} == {"microtubule", "cortex"}
    assert capture.kinetics and capture.commit_on_accept
    assert capture.bidirectional and capture.adjoint_transfer_required


@pytest.mark.parametrize(
    "missing",
    [MT_CORTEX_CAPTURE, MT_NUCLEUS_LINC, MT_SF_SPECTRAPLAKIN, MT_CYTOSOL_TRANSFER],
)
def test_rig_refuses_to_hide_a_required_graph_edge(missing: str) -> None:
    state, *_ = _state()
    with pytest.raises(ValueError, match=missing):
        MicrotubuleRig(_architecture_without(missing), state)


def test_cortex_capture_requires_motor_family_and_accepted_kinetics() -> None:
    state, *_ = _state()
    wrong_family = ConnectorContract(
        MT_CORTEX_CAPTURE,
        ConnectorFamily.LINC,
        "microtubule",
        "cortex",
        True,
        True,
    )
    with pytest.raises(ValueError, match="ConnectorFamily.MOTOR"):
        MicrotubuleRig(_architecture_replacing(wrong_family), state)

    nonkinetic = ConnectorContract(
        MT_CORTEX_CAPTURE,
        ConnectorFamily.MOTOR,
        "microtubule",
        "cortex",
        False,
        False,
    )
    with pytest.raises(ValueError, match="kinetic and accepted-step committed"):
        MicrotubuleRig(_architecture_replacing(nonkinetic), state)


def test_state_owns_independent_cuda_mtoc_rod_and_topology_arrays() -> None:
    state, *_ = _state()
    view = state.geometry()

    assert view.position_d is state.position_d
    assert view.mtoc_position_d is state.mtoc_position_d
    assert view.plus_end_node_d is state.plus_end_node_d
    assert view.topology_epoch_d is state.topology_epoch_d
    assert view.n_mt == 2
    assert state.position_d.ptr != state.mtoc_position_d.ptr

    with pytest.raises(ValueError, match="must not alias"):
        replace(state, mtoc_force_d=state.mtoc_position_d)

    other_gpu = _FakeDevice(alias="cuda:1")
    with pytest.raises(ValueError, match="share one CUDA device"):
        replace(state, force_d=_FakeArray(999, device=other_gpu))


def test_transaction_must_cover_endpoint_phase_and_topology_epoch() -> None:
    state, _, transaction, _ = _state()
    missing_epoch = _TransactionSpy(covered=transaction.covered[:-1])
    with pytest.raises(ValueError, match="transaction must cover"):
        replace(state, transaction=missing_epoch)

    wrong_owner = _TransactionSpy(covered=transaction.covered, component_name="global_monolith")
    with pytest.raises(ValueError, match="owned by the microtubule component"):
        replace(state, transaction=wrong_owner)


def test_rig_has_no_private_capture_linc_or_motor_state() -> None:
    state, *_ = _state()
    rig = MicrotubuleRig(reference_cell_architecture(), state)
    bindings, capture, linc, spectraplakin, transfer, motor = _bindings()

    assert set(type(rig).__slots__) == {"architecture", "state"}
    assert not hasattr(rig, "cortex_capture")
    assert bindings.cortex_capture is capture
    assert bindings.nucleus_linc is linc
    assert bindings.sf_spectraplakin is spectraplakin
    assert bindings.cytosol_transfer is transfer
    assert bindings.motor_ports == (motor,)


def test_mechanics_uses_public_rig_view_and_external_endpoints() -> None:
    state, mechanics, *_ = _state()
    rig = MicrotubuleRig(reference_cell_architecture(), state)
    bindings, capture, linc, spectraplakin, transfer, motor = _bindings()
    rig.accumulate_mechanics(bindings)

    assert mechanics.calls == [(state.position_d, state.force_d)]
    assert len(capture.mechanics_calls) == len(linc.mechanics_calls) == 1
    assert len(spectraplakin.mechanics_calls) == len(transfer.mechanics_calls) == 1
    assert len(motor.mechanics_calls) == 1
    capture_view, cortex_endpoint = capture.mechanics_calls[0]
    linc_view, nucleus_endpoint = linc.mechanics_calls[0]
    assert capture_view.position_d is state.position_d
    assert linc_view.topology_epoch_d is state.topology_epoch_d
    assert motor.mechanics_calls[0].plus_end_node_d is state.plus_end_node_d
    assert cortex_endpoint is bindings.cortex
    assert nucleus_endpoint is bindings.nucleus
    assert spectraplakin.mechanics_calls[0][1] is bindings.sf_arc
    assert transfer.mechanics_calls[0][1] is bindings.cytosol


def test_transaction_and_ledger_cover_component_capture_linc_and_motor() -> None:
    state, _, transaction, state_ledger = _state()
    rig = MicrotubuleRig(reference_cell_architecture(), state)
    bindings, capture, linc, spectraplakin, transfer, motor = _bindings()
    accepted_d = object()
    ledger = object()

    rig.snapshot_candidate(bindings)
    rig.rollback(bindings, accepted_d)
    rig.commit_irreversible(bindings, accepted_d, dt_phys=0.025, rng_seed=73)
    rig.accumulate_ledger(bindings, ledger)

    expected = [
        ("snapshot",),
        ("rollback", accepted_d),
        ("commit", accepted_d, 0.025, 73),
    ]
    assert transaction.calls == expected
    assert capture.transaction_calls == expected
    assert linc.transaction_calls == expected
    assert spectraplakin.transaction_calls == expected
    assert transfer.transaction_calls == expected
    assert motor.transaction_calls == expected
    assert state_ledger.calls == [ledger]
    assert capture.ledger_calls == linc.ledger_calls == [ledger]
    assert spectraplakin.ledger_calls == transfer.ledger_calls == motor.ledger_calls == [ledger]


def test_step_bindings_reject_misdirected_capture_linc_and_motor_ports() -> None:
    bindings, capture, linc, spectraplakin, transfer, _ = _bindings()
    with pytest.raises(ValueError, match=MT_CORTEX_CAPTURE):
        replace(bindings, cortex_capture=replace(capture, name="private_tip_spring"))
    with pytest.raises(ValueError, match="incorrect component endpoints"):
        replace(bindings, nucleus_linc=replace(linc, component_b="cortex"))
    with pytest.raises(ValueError, match=MT_SF_SPECTRAPLAKIN):
        replace(bindings, sf_spectraplakin=replace(spectraplakin, name="private_bridge"))
    with pytest.raises(ValueError, match=MT_CYTOSOL_TRANSFER):
        replace(bindings, cytosol_transfer=replace(transfer, name="private_drag"))
    with pytest.raises(ValueError, match="must have a microtubule endpoint"):
        replace(bindings, motor_ports=(_MotorPortSpy(component_a="nucleus", component_b="cortex"),))


@dataclass(slots=True)
class _LegacyMicrotubuleSpy:
    node_off: int = 0
    n_nodes: int = 8
    n_mt: int = 2
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


def test_legacy_adapter_is_local_and_cannot_claim_dynamic_completion() -> None:
    legacy = _LegacyMicrotubuleSpy()
    adapter = LegacyMicrotubuleBendingAdapter(legacy)
    pos = _FakeArray(1)
    force = _FakeArray(2)
    adapter.accumulate(pos, force)

    assert adapter.status == STATIC_ADAPTER_STATUS
    assert legacy.calls == [(pos, force)]
    with pytest.raises(ValueError, match="node_off == 0"):
        LegacyMicrotubuleBendingAdapter(_LegacyMicrotubuleSpy(node_off=7))
    with pytest.raises(ValueError, match="non-empty aster"):
        LegacyMicrotubuleBendingAdapter(_LegacyMicrotubuleSpy(n_mt=0))

    state, *_ = _state()
    with pytest.raises(ValueError, match="cannot be relabelled"):
        replace(state, status="DYNAMIC_INSTABILITY_COMPLETE")


def test_legacy_adapter_node_count_must_match_component_state() -> None:
    state, *_ = _state()
    adapter = LegacyMicrotubuleBendingAdapter(_LegacyMicrotubuleSpy(n_nodes=7))
    with pytest.raises(ValueError, match="expects 7 nodes"):
        replace(state, mechanics=adapter)


# --------------------------------------------------------------------------------------------------
# KERNEL_BOUND increment: MTOC anchor + dynamic-instability topology (host oracles CPU-green;
# real Warp kernels CUDA-gated, mirroring test_load_path.py).
# --------------------------------------------------------------------------------------------------


def _rod_spec() -> MicrotubuleRodBuildSpec:
    """A 2-arm aster: arm0 = +x nodes 0..3, arm1 = +y nodes 4..7, capacity 4 each, seg = 1 um."""
    position = np.array(
        [
            [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [4.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0], [0.0, 4.0, 0.0],
        ],
        dtype=np.float64,
    )
    return MicrotubuleRodBuildSpec(
        position=position,
        fiber_offset=np.array([0, 4, 8], np.int32),
        mtoc_position=np.array([0.0, 0.0, 0.0], np.float64),
        active_count=np.array([3, 3], np.int32),
        phase=np.array([int(MicrotubulePhase.GROWING), int(MicrotubulePhase.SHRINKING)], np.int32),
        v_grow_um_per_s=np.array([2.0, 0.0], np.float64),
        v_shrink_um_per_s=np.array([0.0, 2.0], np.float64),
        k_hub_pn_per_um=10.0,
        anchor_rest_um=0.0,
        seg_um=1.0,
    )


def test_phase_enum_labels_are_stable() -> None:
    assert (int(MicrotubulePhase.GROWING), int(MicrotubulePhase.SHRINKING), int(MicrotubulePhase.PAUSED)) == (
        0,
        1,
        2,
    )


def test_mtoc_anchor_reference_is_equal_and_opposite_and_signed() -> None:
    ref = mtoc_anchor_reference([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], k_hub_pn_per_um=10.0, rest_um=0.0)
    # base pulled toward the MTOC (−x); MTOC reaction is the exact negative (Newton 3rd).
    np.testing.assert_allclose(ref.force_on_base, [-10.0, 0.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(ref.force_on_mtoc, [10.0, 0.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(ref.force_on_base + ref.force_on_mtoc, [0.0, 0.0, 0.0], atol=1e-12)
    assert ref.load_pn == pytest.approx(10.0)
    assert ref.energy_pn_um == pytest.approx(5.0)
    # at rest length the anchor carries no force
    rest = mtoc_anchor_reference([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], k_hub_pn_per_um=10.0, rest_um=1.0)
    np.testing.assert_allclose(rest.force_on_base, [0.0, 0.0, 0.0], atol=1e-12)
    assert rest.load_pn == pytest.approx(0.0)


def test_mtoc_anchor_reference_rejects_bad_params() -> None:
    with pytest.raises(ValueError, match="stiffness"):
        mtoc_anchor_reference([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], k_hub_pn_per_um=0.0, rest_um=0.0)
    with pytest.raises(ValueError, match="rest length"):
        mtoc_anchor_reference([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], k_hub_pn_per_um=1.0, rest_um=-1.0)
    with pytest.raises(ValueError, match="3-vectors"):
        mtoc_anchor_reference([1.0, 0.0], [0.0, 0.0, 0.0], k_hub_pn_per_um=1.0, rest_um=0.0)


def test_dynamic_instability_reference_growth_activates_one_span() -> None:
    ph, ac, pe, partial = dynamic_instability_reference(
        phase=int(MicrotubulePhase.GROWING),
        active_count=3, plus_end_node=2, partial_length_um=0.0,
        fiber_offset_start=0, fiber_offset_end=4,
        v_grow_um_per_s=2.0, v_shrink_um_per_s=0.0, seg_um=1.0, dt_phys=0.6,
    )
    # 2 um/s * 0.6 s = 1.2 um > 1 um seg -> activate one dormant node, carry 0.2 um remainder.
    assert (ph, ac, pe) == (int(MicrotubulePhase.GROWING), 4, 3)
    assert partial == pytest.approx(0.2)


def test_dynamic_instability_reference_growth_clamps_at_capacity() -> None:
    # plus-end already at the last allocated node (3) of a capacity-4 arm: cannot activate further.
    ph, ac, pe, partial = dynamic_instability_reference(
        phase=int(MicrotubulePhase.GROWING),
        active_count=4, plus_end_node=3, partial_length_um=0.0,
        fiber_offset_start=0, fiber_offset_end=4,
        v_grow_um_per_s=2.0, v_shrink_um_per_s=0.0, seg_um=1.0, dt_phys=0.6,
    )
    assert (ac, pe) == (4, 3)
    assert partial == pytest.approx(1.0)  # clamped at capacity, no host reallocation


def test_dynamic_instability_reference_shrink_and_minimum_rod() -> None:
    ph, ac, pe, partial = dynamic_instability_reference(
        phase=int(MicrotubulePhase.SHRINKING),
        active_count=3, plus_end_node=6, partial_length_um=0.0,
        fiber_offset_start=4, fiber_offset_end=8,
        v_grow_um_per_s=0.0, v_shrink_um_per_s=2.0, seg_um=1.0, dt_phys=0.6,
    )
    assert (ac, pe) == (2, 5)
    assert partial == pytest.approx(-0.2)
    # a 2-node rod cannot shrink below the minimum viable bending rod.
    _, ac2, pe2, partial2 = dynamic_instability_reference(
        phase=int(MicrotubulePhase.SHRINKING),
        active_count=2, plus_end_node=5, partial_length_um=0.0,
        fiber_offset_start=4, fiber_offset_end=8,
        v_grow_um_per_s=0.0, v_shrink_um_per_s=2.0, seg_um=1.0, dt_phys=0.6,
    )
    assert (ac2, pe2) == (2, 5)
    assert partial2 == pytest.approx(-1.0)


def test_dynamic_instability_reference_paused_is_a_noop() -> None:
    ph, ac, pe, partial = dynamic_instability_reference(
        phase=int(MicrotubulePhase.PAUSED),
        active_count=3, plus_end_node=2, partial_length_um=0.3,
        fiber_offset_start=0, fiber_offset_end=4,
        v_grow_um_per_s=5.0, v_shrink_um_per_s=5.0, seg_um=1.0, dt_phys=1.0,
    )
    assert (ph, ac, pe) == (int(MicrotubulePhase.PAUSED), 3, 2)
    assert partial == pytest.approx(0.3)


def test_rod_build_spec_derives_plus_end_and_validates() -> None:
    spec = _rod_spec()
    assert spec.n_mt == 2 and spec.n_nodes == 8
    np.testing.assert_array_equal(spec.plus_end_node(), np.array([2, 6], np.int32))

    with pytest.raises(ValueError, match=r"\[2, per-arm capacity\]"):
        replace(spec, active_count=np.array([5, 3], np.int32))
    with pytest.raises(ValueError, match="start at 0 and end"):
        replace(spec, fiber_offset=np.array([0, 4, 9], np.int32))
    with pytest.raises(ValueError, match="k_hub"):
        replace(spec, k_hub_pn_per_um=0.0)
    with pytest.raises(ValueError, match="speeds must be nonnegative"):
        replace(spec, v_grow_um_per_s=np.array([-1.0, 0.0], np.float64))
    with pytest.raises(ValueError, match="MicrotubulePhase"):
        replace(spec, phase=np.array([7, 1], np.int32))


def test_rod_runtime_requires_cuda_on_the_dev_mac() -> None:
    if any(device.is_cuda for device in wp.get_devices()):
        pytest.skip("CUDA present; requirement gate is only meaningful on a CPU-only host")
    with pytest.raises(RuntimeError, match="Warp-CUDA-only"):
        MicrotubuleRodRuntime(_rod_spec())


_CUDA_DEVICE = next((device for device in wp.get_devices() if device.is_cuda), None)


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="MT rod runtime kernels require CUDA")
def test_cuda_mtoc_anchor_and_accepted_di_topology() -> None:
    device = str(_CUDA_DEVICE)
    spec = _rod_spec()
    runtime = MicrotubuleRodRuntime(spec, device=device)

    # --- MTOC minus-end anchor: Newton-3rd equal-and-opposite reaction ---
    runtime.zero_force()
    runtime.accumulate_mtoc_anchor()
    wp.synchronize_device(device)
    force = runtime.force_d.numpy()
    mtoc_force = runtime.mtoc_force_d.numpy()
    base0 = mtoc_anchor_reference(spec.position[0], spec.mtoc_position, 10.0, 0.0)
    base1 = mtoc_anchor_reference(spec.position[4], spec.mtoc_position, 10.0, 0.0)
    np.testing.assert_allclose(force[0], base0.force_on_base, atol=1e-10)
    np.testing.assert_allclose(force[4], base1.force_on_base, atol=1e-10)
    # every reaction lands on the shared MTOC node; total is the negated arm sum (global force balance).
    np.testing.assert_allclose(
        mtoc_force[0] + force.sum(axis=0), [0.0, 0.0, 0.0], atol=1e-10
    )
    np.testing.assert_allclose(runtime.anchor_load_d.numpy(), [10.0, 10.0], atol=1e-10)

    # --- dynamic instability: propose into candidates, reject leaves committed intact ---
    runtime.snapshot_candidate()
    runtime.propose_dynamic_instability(0.6)
    wp.synchronize_device(device)
    exp0 = dynamic_instability_reference(
        phase=int(MicrotubulePhase.GROWING), active_count=3, plus_end_node=2,
        partial_length_um=0.0, fiber_offset_start=0, fiber_offset_end=4,
        v_grow_um_per_s=2.0, v_shrink_um_per_s=0.0, seg_um=1.0, dt_phys=0.6,
    )
    exp1 = dynamic_instability_reference(
        phase=int(MicrotubulePhase.SHRINKING), active_count=3, plus_end_node=6,
        partial_length_um=0.0, fiber_offset_start=4, fiber_offset_end=8,
        v_grow_um_per_s=0.0, v_shrink_um_per_s=2.0, seg_um=1.0, dt_phys=0.6,
    )
    np.testing.assert_array_equal(
        runtime.candidate_active_count_d.numpy(), [exp0[1], exp1[1]]
    )
    np.testing.assert_array_equal(runtime.candidate_plus_end_d.numpy(), [exp0[2], exp1[2]])

    rejected = wp.array(np.array([0], np.int32), dtype=wp.int32, device=device)
    runtime.commit_irreversible(rejected, 0.6, 11)
    wp.synchronize_device(device)
    np.testing.assert_array_equal(runtime.active_count_d.numpy(), [3, 3])  # committed untouched
    np.testing.assert_array_equal(runtime.plus_end_node_d.numpy(), [2, 6])
    assert runtime.topology_epoch_d.numpy()[0] == 0
    # rejection resets candidate scratch to committed
    np.testing.assert_array_equal(runtime.candidate_active_count_d.numpy(), [3, 3])

    # --- accepted commit advances committed topology + epoch exactly once ---
    runtime.propose_dynamic_instability(0.6)
    accepted = wp.array(np.array([1], np.int32), dtype=wp.int32, device=device)
    runtime.commit_irreversible(accepted, 0.6, 12)
    wp.synchronize_device(device)
    np.testing.assert_array_equal(runtime.active_count_d.numpy(), [exp0[1], exp1[1]])
    np.testing.assert_array_equal(runtime.plus_end_node_d.numpy(), [exp0[2], exp1[2]])
    np.testing.assert_allclose(runtime.partial_length_d.numpy(), [exp0[3], exp1[3]], atol=1e-12)
    assert runtime.topology_epoch_d.numpy()[0] == 1


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="MT rod runtime kernels require CUDA")
def test_cuda_runtime_binds_into_state_owner_seam() -> None:
    """The real-kernel runtime is accepted as the MicrotubuleRigStateOwner transaction delegate."""
    device = str(_CUDA_DEVICE)
    runtime = MicrotubuleRodRuntime(_rod_spec(), device=device)

    @dataclass(slots=True)
    class _BendingSpy:
        n_nodes: int = 8
        calls: list[tuple[object, object]] = field(default_factory=list)

        def accumulate(self, pos: wp.array, force: wp.array) -> None:
            self.calls.append((pos, force))

    @dataclass(slots=True)
    class _Ledger:
        def accumulate_ledger(self, ledger: object) -> None:
            return None

    bending = _BendingSpy()
    owner = MicrotubuleRigStateOwner(
        n_mt=runtime.n_mt,
        position_d=runtime.position_d,
        force_d=runtime.force_d,
        mtoc_position_d=runtime.mtoc_position_d,
        mtoc_force_d=runtime.mtoc_force_d,
        fiber_offset_d=runtime.fiber_offset_d,
        active_count_d=runtime.active_count_d,
        plus_end_node_d=runtime.plus_end_node_d,
        phase_d=runtime.phase_d,
        topology_epoch_d=runtime.topology_epoch_d,
        mechanics=bending,
        transaction=runtime,
        ledger=_Ledger(),
    )
    # reused bending kernel launches through the seam against the runtime-owned arrays
    runtime.accumulate_bending(bending)
    assert bending.calls == [(runtime.position_d, runtime.force_d)]
    # the runtime's real commit kernel drives the owner's transaction path
    owner.snapshot_candidate()
    runtime.propose_dynamic_instability(0.6)
    accepted = wp.array(np.array([1], np.int32), dtype=wp.int32, device=device)
    owner.commit_irreversible(accepted, 0.6, 5)
    wp.synchronize_device(device)
    assert runtime.topology_epoch_d.numpy()[0] == 1


# ── Dynamic-instability catastrophe/rescue: PI-GAP rate slot (SEAMED → KERNEL_BOUND) ─────────────


def test_di_rate_card_is_an_unset_slot_with_no_default() -> None:
    # the card carries no default rates: a caller must supply sourced numbers (MCF7 DI is a PI GAP).
    card = DynamicInstabilityRateCard(
        f_catastrophe_per_s=0.03, f_rescue_per_s=0.2, provenance="PI-ratified card"
    )
    assert card.ratified is False
    with pytest.raises(ValueError, match="finite and nonnegative"):
        DynamicInstabilityRateCard(f_catastrophe_per_s=-1.0, f_rescue_per_s=0.2, provenance="x")
    with pytest.raises(ValueError, match="finite and nonnegative"):
        DynamicInstabilityRateCard(f_catastrophe_per_s=0.0, f_rescue_per_s=float("nan"), provenance="x")
    with pytest.raises(ValueError, match="provenance"):
        DynamicInstabilityRateCard(f_catastrophe_per_s=0.0, f_rescue_per_s=0.0, provenance="   ")


def test_di_rusan_proxy_carries_found_rates_but_is_not_ratified() -> None:
    proxy = DynamicInstabilityRateCard.rusan_2001_llcpk_proxy()
    # FOUND in-vivo mammalian-epithelial plus-end rates (Rusan 2001) — put in, but flagged proxy not ratified.
    assert proxy.f_catastrophe_per_s == pytest.approx(0.026)
    assert proxy.f_rescue_per_s == pytest.approx(0.175)
    assert proxy.ratified is False
    assert "Rusan 2001" in proxy.provenance and "ratification pending" in proxy.provenance.lower()


def test_stochastic_reference_forced_catastrophe_switches_and_shrinks() -> None:
    # a growing plus end with a certain catastrophe (draw below p_cat) flips to shrinking and loses one span.
    ph, ac, pe, partial = dynamic_instability_stochastic_reference(
        phase=int(MicrotubulePhase.GROWING), active_count=3, plus_end_node=2,
        partial_length_um=0.0, fiber_offset_start=0, fiber_offset_end=4,
        v_grow_um_per_s=2.0, v_shrink_um_per_s=2.0,
        f_catastrophe_per_s=100.0, f_rescue_per_s=0.0, uniform_draw=0.0,
        seg_um=1.0, dt_phys=1.0,
    )
    assert ph == int(MicrotubulePhase.SHRINKING)
    assert (ac, pe) == (2, 1)
    assert partial == pytest.approx(-1.0)


def test_stochastic_reference_forced_rescue_switches_and_grows() -> None:
    ph, ac, pe, partial = dynamic_instability_stochastic_reference(
        phase=int(MicrotubulePhase.SHRINKING), active_count=3, plus_end_node=6,
        partial_length_um=0.0, fiber_offset_start=4, fiber_offset_end=8,
        v_grow_um_per_s=2.0, v_shrink_um_per_s=2.0,
        f_catastrophe_per_s=0.0, f_rescue_per_s=100.0, uniform_draw=0.0,
        seg_um=1.0, dt_phys=1.0,
    )
    assert ph == int(MicrotubulePhase.GROWING)
    assert (ac, pe) == (4, 7)
    assert partial == pytest.approx(1.0)


def test_stochastic_reference_zero_rate_matches_deterministic_proposal() -> None:
    # with a zero-rate card the stochastic path must reduce exactly to the deterministic grow/shrink proposal.
    kwargs = dict(
        active_count=3, plus_end_node=2, partial_length_um=0.0,
        fiber_offset_start=0, fiber_offset_end=4,
        v_grow_um_per_s=2.0, v_shrink_um_per_s=2.0, seg_um=1.0, dt_phys=1.0,
    )
    stoch = dynamic_instability_stochastic_reference(
        phase=int(MicrotubulePhase.GROWING),
        f_catastrophe_per_s=0.0, f_rescue_per_s=0.0, uniform_draw=0.5, **kwargs,
    )
    determ = dynamic_instability_reference(phase=int(MicrotubulePhase.GROWING), **kwargs)
    assert stoch == determ


def test_stochastic_reference_high_draw_keeps_phase() -> None:
    # a draw of 1.0 exceeds any p<1, so no transition happens even with a finite catastrophe rate.
    ph, _, _, _ = dynamic_instability_stochastic_reference(
        phase=int(MicrotubulePhase.GROWING), active_count=3, plus_end_node=2,
        partial_length_um=0.0, fiber_offset_start=0, fiber_offset_end=4,
        v_grow_um_per_s=2.0, v_shrink_um_per_s=2.0,
        f_catastrophe_per_s=0.026, f_rescue_per_s=0.175, uniform_draw=1.0,
        seg_um=1.0, dt_phys=1.0,
    )
    assert ph == int(MicrotubulePhase.GROWING)
    with pytest.raises(ValueError, match="uniform_draw"):
        dynamic_instability_stochastic_reference(
            phase=int(MicrotubulePhase.GROWING), active_count=3, plus_end_node=2,
            partial_length_um=0.0, fiber_offset_start=0, fiber_offset_end=4,
            v_grow_um_per_s=2.0, v_shrink_um_per_s=2.0,
            f_catastrophe_per_s=0.0, f_rescue_per_s=0.0, uniform_draw=1.5,
            seg_um=1.0, dt_phys=1.0,
        )
