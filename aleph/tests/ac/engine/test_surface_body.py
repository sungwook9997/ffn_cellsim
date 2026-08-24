"""Structural/runtime-seam gates for the reduced membrane/cortex Surface Body."""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorContract,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine.runtime import CytosolFieldEndpoint
from aleph.engine.surface_body import (
    ERM_CONNECTOR,
    MEMBRANE_FLUID_BOUNDARY,
    SURFACE_POROUS_TRANSFER,
    CortexFilamentMechanics,
    LegacyMembraneCompartmentAdapter,
    MembranePressureFluidBoundary,
    SurfaceBody,
    SurfaceComponentStateOwner,
    SurfaceStepBindings,
)


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """Metadata-only CUDA array double; it never executes authoritative physics."""

    ptr: int
    shape: tuple[int, ...] = (4,)
    dtype: object = wp.vec3d
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _MechanicsSpy:
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


@dataclass(slots=True)
class _TransactionSpy:
    calls: list[tuple[object, ...]] = field(default_factory=list)

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
class _PairConnectorSpy:
    name: str = ERM_CONNECTOR
    component_a: str = "membrane"
    component_b: str = "cortex"
    mechanics_calls: list[tuple[object, object, object, object]] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_pair(
        self,
        pos_a: wp.array,
        force_a: wp.array,
        pos_b: wp.array,
        force_b: wp.array,
    ) -> None:
        self.mechanics_calls.append((pos_a, force_a, pos_b, force_b))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _FluidBoundarySpy:
    name: str = MEMBRANE_FLUID_BOUNDARY
    component_a: str = "membrane"
    component_b: str = "cytosol"
    mechanics_calls: list[tuple[object, object]] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        self.mechanics_calls.append((pos, force))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _PorousTransferSpy:
    name: str = SURFACE_POROUS_TRANSFER
    component_a: str = "cortex"
    component_b: str = "cytosol"
    mechanics_calls: list[tuple[object, object]] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_transfer(self, solid: object, cytosol: object) -> None:
        self.mechanics_calls.append((solid, cytosol))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


def _architecture_with_membrane_fluid_boundary(
    *, endpoint_a: str = "membrane", endpoint_b: str = "cytosol"
) -> CellArchitecture:
    base = reference_cell_architecture()
    boundary = ConnectorContract(
        MEMBRANE_FLUID_BOUNDARY,
        ConnectorFamily.IMMERSED_TRANSFER,
        endpoint_a,
        endpoint_b,
        kinetics=False,
        commit_on_accept=False,
        endpoint_role_a="live membrane quadrature",
        endpoint_role_b="cytosol boundary field",
    )
    connectors = tuple(
        boundary if connector.name == MEMBRANE_FLUID_BOUNDARY else connector
        for connector in base.connectors
    )
    return CellArchitecture(base.components, connectors)


def _owner(name: str, ptr_base: int) -> tuple[
    SurfaceComponentStateOwner, _MechanicsSpy, _TransactionSpy, _LedgerSpy
]:
    mechanics = _MechanicsSpy()
    transaction = _TransactionSpy()
    ledger = _LedgerSpy()
    owner = SurfaceComponentStateOwner(
        name=name,
        position_d=_FakeArray(ptr_base),
        force_d=_FakeArray(ptr_base + 1),
        mechanics=mechanics,
        transaction=transaction,
        ledger=ledger,
    )
    return owner, mechanics, transaction, ledger


def _body_and_bindings() -> tuple[
    SurfaceBody,
    SurfaceStepBindings,
    tuple[_MechanicsSpy, _TransactionSpy, _LedgerSpy],
    tuple[_MechanicsSpy, _TransactionSpy, _LedgerSpy],
    _PairConnectorSpy,
    _FluidBoundarySpy,
    _PorousTransferSpy,
]:
    membrane, membrane_mechanics, membrane_tx, membrane_ledger = _owner("membrane", 100)
    cortex, cortex_mechanics, cortex_tx, cortex_ledger = _owner("cortex", 200)
    body = SurfaceBody(_architecture_with_membrane_fluid_boundary(), membrane, cortex)
    erm = _PairConnectorSpy()
    fluid = _FluidBoundarySpy()
    porous = _PorousTransferSpy()
    cytosol = CytosolFieldEndpoint(
        "cytosol",
        _FakeArray(400, shape=(2, 2, 2), dtype=wp.float64),
        _FakeArray(401, shape=(2, 2, 2), dtype=wp.float64),
    )
    bindings = SurfaceStepBindings(erm, fluid, porous, cytosol)
    return (
        body,
        bindings,
        (membrane_mechanics, membrane_tx, membrane_ledger),
        (cortex_mechanics, cortex_tx, cortex_ledger),
        erm,
        fluid,
        porous,
    )


def test_reference_graph_registers_the_required_membrane_fluid_boundary() -> None:
    membrane, *_ = _owner("membrane", 100)
    cortex, *_ = _owner("cortex", 200)
    architecture = reference_cell_architecture()
    body = SurfaceBody(architecture, membrane, cortex)
    boundary = next(
        connector
        for connector in architecture.connectors
        if connector.name == MEMBRANE_FLUID_BOUNDARY
    )
    assert body.architecture is architecture
    assert boundary.family is ConnectorFamily.FLUID_BOUNDARY
    assert {boundary.component_a, boundary.component_b} == {"membrane", "cytosol"}


def test_membrane_fluid_boundary_must_join_membrane_and_cytosol() -> None:
    membrane, *_ = _owner("membrane", 100)
    cortex, *_ = _owner("cortex", 200)
    architecture = _architecture_with_membrane_fluid_boundary(
        endpoint_a="cortex", endpoint_b="cytosol"
    )
    with pytest.raises(ValueError, match="must join membrane and cytosol"):
        SurfaceBody(architecture, membrane, cortex)


def test_surface_body_enforces_distinct_state_owners_and_kinematics() -> None:
    architecture = _architecture_with_membrane_fluid_boundary()
    membrane, *_ = _owner("membrane", 100)
    cortex, *_ = _owner("cortex", 200)
    body = SurfaceBody(architecture, membrane, cortex)
    assert body.membrane is membrane and body.cortex is cortex
    assert body.membrane.position_d.ptr != body.cortex.position_d.ptr

    cortex_mechanics = _MechanicsSpy()
    aliasing_cortex = SurfaceComponentStateOwner(
        name="cortex",
        position_d=_FakeArray(100),  # distinct object, same device pointer: storage alias
        force_d=_FakeArray(300),
        mechanics=cortex_mechanics,
        transaction=_TransactionSpy(),
        ledger=_LedgerSpy(),
    )
    with pytest.raises(ValueError, match="distinct position kinematics"):
        SurfaceBody(architecture, membrane, aliasing_cortex)


def test_state_owner_rejects_cpu_arrays_and_position_force_aliasing() -> None:
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        SurfaceComponentStateOwner(
            "membrane",
            _FakeArray(1, device=cpu),
            _FakeArray(2, device=cpu),
            _MechanicsSpy(),
            _TransactionSpy(),
            _LedgerSpy(),
        )

    with pytest.raises(ValueError, match="must not alias"):
        SurfaceComponentStateOwner(
            "membrane",
            _FakeArray(1),
            _FakeArray(1),
            _MechanicsSpy(),
            _TransactionSpy(),
            _LedgerSpy(),
        )


def test_facade_does_not_own_erm_state() -> None:
    body, bindings, *_ = _body_and_bindings()
    assert not hasattr(body, "erm_connector")
    assert not hasattr(body, "erm_state")
    assert bindings.erm_connector is body.transaction_participants(bindings)[2]
    assert set(type(body).__slots__) == {"architecture", "membrane", "cortex"}


def test_mechanics_wiring_uses_each_owner_array_and_external_graph_connector() -> None:
    body, bindings, membrane_spies, cortex_spies, erm, fluid, porous = _body_and_bindings()
    body.accumulate_mechanics(bindings)

    assert membrane_spies[0].calls == [(body.membrane.position_d, body.membrane.force_d)]
    assert cortex_spies[0].calls == [(body.cortex.position_d, body.cortex.force_d)]
    assert erm.mechanics_calls == [
        (
            body.membrane.position_d,
            body.membrane.force_d,
            body.cortex.position_d,
            body.cortex.force_d,
        )
    ]
    assert fluid.mechanics_calls == [(body.membrane.position_d, body.membrane.force_d)]
    surface, cytosol = porous.mechanics_calls[0]
    assert surface.position_d is body.cortex.position_d
    assert surface.force_d is body.cortex.force_d
    assert cytosol is bindings.cytosol


def test_transaction_and_ledger_hooks_cover_all_surface_field_paths() -> None:
    body, bindings, membrane_spies, cortex_spies, erm, fluid, porous = _body_and_bindings()
    accepted_d = object()
    ledger = object()

    body.snapshot_candidate(bindings)
    body.rollback(bindings, accepted_d)
    body.commit_irreversible(bindings, accepted_d, dt_phys=0.05, rng_seed=41)
    body.accumulate_ledger(bindings, ledger)

    expected_tx = [
        ("snapshot",),
        ("rollback", accepted_d),
        ("commit", accepted_d, 0.05, 41),
    ]
    assert membrane_spies[1].calls == expected_tx
    assert cortex_spies[1].calls == expected_tx
    assert erm.transaction_calls == expected_tx
    assert fluid.transaction_calls == expected_tx
    assert porous.transaction_calls == expected_tx
    assert membrane_spies[2].calls == [ledger]
    assert cortex_spies[2].calls == [ledger]
    assert erm.ledger_calls == [ledger]
    assert fluid.ledger_calls == [ledger]
    assert porous.ledger_calls == [ledger]


def test_step_bindings_reject_private_or_misdirected_connectors() -> None:
    cytosol = CytosolFieldEndpoint(
        "cytosol",
        _FakeArray(400, shape=(2, 2, 2), dtype=wp.float64),
        _FakeArray(401, shape=(2, 2, 2), dtype=wp.float64),
    )
    with pytest.raises(ValueError, match=ERM_CONNECTOR):
        SurfaceStepBindings(
            _PairConnectorSpy(name="private_erm"),
            _FluidBoundarySpy(),
            _PorousTransferSpy(),
            cytosol,
        )
    with pytest.raises(ValueError, match="membrane and cortex"):
        SurfaceStepBindings(
            _PairConnectorSpy(component_b="cytosol"),
            _FluidBoundarySpy(),
            _PorousTransferSpy(),
            cytosol,
        )
    with pytest.raises(ValueError, match=MEMBRANE_FLUID_BOUNDARY):
        SurfaceStepBindings(
            _PairConnectorSpy(),
            _FluidBoundarySpy(name="hidden_pressure_callback"),
            _PorousTransferSpy(),
            cytosol,
        )
    with pytest.raises(ValueError, match=SURFACE_POROUS_TRANSFER):
        SurfaceStepBindings(
            _PairConnectorSpy(),
            _FluidBoundarySpy(),
            _PorousTransferSpy(name="private_drag"),
            cytosol,
        )


@dataclass(slots=True)
class _LegacyMembraneSpy:
    node_off: int = 0
    n_verts: int = 4
    n_erm: int = 0
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


def test_legacy_membrane_adapter_reuses_only_local_erm_free_mechanics() -> None:
    legacy = _LegacyMembraneSpy()
    adapter = LegacyMembraneCompartmentAdapter(legacy)
    pos = _FakeArray(10)
    force = _FakeArray(11)
    adapter.accumulate(pos, force)
    assert adapter.n_vertices == 4
    assert legacy.calls == [(pos, force)]

    with pytest.raises(ValueError, match="node_off == 0"):
        LegacyMembraneCompartmentAdapter(_LegacyMembraneSpy(node_off=5))
    with pytest.raises(ValueError, match="cannot own ERM state"):
        LegacyMembraneCompartmentAdapter(_LegacyMembraneSpy(n_erm=1))


def test_legacy_membrane_adapter_vertex_count_must_match_owned_state() -> None:
    adapter = LegacyMembraneCompartmentAdapter(_LegacyMembraneSpy(n_verts=5))
    with pytest.raises(ValueError, match="expects 5 vertices"):
        SurfaceComponentStateOwner(
            "membrane",
            _FakeArray(1, shape=(4,)),
            _FakeArray(2, shape=(4,)),
            adapter,
            _TransactionSpy(),
            _LedgerSpy(),
        )


# --- kernel-bound cortex mechanics (link_spring + cytosim_bending) ---------------------------------


@dataclass(slots=True)
class _LaunchRecorder:
    """CUDA-launch double: records ``wp.launch`` calls without executing a kernel on the host."""

    calls: list[dict[str, object]] = field(default_factory=list)

    def __call__(self, kernel, *, dim, inputs, device) -> None:
        self.calls.append({"kernel": kernel, "dim": dim, "inputs": inputs, "device": device})


_LINK_KERNEL = object()
_BENDING_KERNEL = object()


def _cortex_filament_mechanics(recorder: _LaunchRecorder, *, n_vertices: int = 4) -> CortexFilamentMechanics:
    return CortexFilamentMechanics(
        device="cuda:0",
        n_vertices=n_vertices,
        links_d=_FakeArray(500, shape=(3, 2), dtype=wp.int32),
        link_k_d=_FakeArray(501, shape=(3,), dtype=wp.float64),
        link_r0_d=_FakeArray(502, shape=(3,), dtype=wp.float64),
        bend_triples_d=_FakeArray(503, shape=(2, 3), dtype=wp.int32),
        bend_alpha_d=_FakeArray(504, shape=(2,), dtype=wp.float64),
        link_kernel=_LINK_KERNEL,
        bending_kernel=_BENDING_KERNEL,
        launch=recorder,
    )


def test_cortex_mechanics_launches_link_and_bending_kernels_on_cortex_arrays() -> None:
    recorder = _LaunchRecorder()
    mechanics = _cortex_filament_mechanics(recorder)
    pos = _FakeArray(600)
    force = _FakeArray(601)
    mechanics.accumulate(pos, force)

    assert len(recorder.calls) == 2
    link_call, bend_call = recorder.calls
    assert link_call["kernel"] is _LINK_KERNEL
    assert link_call["dim"] == 3  # one thread per link/segment
    assert link_call["inputs"] == [pos, mechanics.links_d, mechanics.link_k_d, mechanics.link_r0_d, force]
    assert link_call["device"] == "cuda:0"
    assert bend_call["kernel"] is _BENDING_KERNEL
    assert bend_call["dim"] == 2  # one thread per bending triple
    assert bend_call["inputs"] == [pos, mechanics.bend_triples_d, mechanics.bend_alpha_d, force]


def test_cortex_mechanics_validates_topology_and_rejects_cpu_or_mismatched_arrays() -> None:
    recorder = _LaunchRecorder()
    with pytest.raises(TypeError, match="dtype wp.int32"):
        CortexFilamentMechanics(
            device="cuda:0",
            n_vertices=4,
            links_d=_FakeArray(500, shape=(3, 2), dtype=wp.float64),
            link_k_d=_FakeArray(501, shape=(3,), dtype=wp.float64),
            link_r0_d=_FakeArray(502, shape=(3,), dtype=wp.float64),
            bend_triples_d=_FakeArray(503, shape=(2, 3), dtype=wp.int32),
            bend_alpha_d=_FakeArray(504, shape=(2,), dtype=wp.float64),
            link_kernel=_LINK_KERNEL,
            bending_kernel=_BENDING_KERNEL,
            launch=recorder,
        )
    with pytest.raises(ValueError, match="match its topology row count"):
        CortexFilamentMechanics(
            device="cuda:0",
            n_vertices=4,
            links_d=_FakeArray(500, shape=(3, 2), dtype=wp.int32),
            link_k_d=_FakeArray(501, shape=(2,), dtype=wp.float64),  # wrong length vs 3 links
            link_r0_d=_FakeArray(502, shape=(3,), dtype=wp.float64),
            bend_triples_d=_FakeArray(503, shape=(2, 3), dtype=wp.int32),
            bend_alpha_d=_FakeArray(504, shape=(2,), dtype=wp.float64),
            link_kernel=_LINK_KERNEL,
            bending_kernel=_BENDING_KERNEL,
            launch=recorder,
        )
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    # The rejection is driven by the cpu ARRAYS (_device_is_cuda checks each array's device), NOT the device
    # string; keep the string on the CUDA lane so the GPU-only static contract does not false-positive on a
    # literal "cpu" while this negative test still proves cpu arrays are refused.
    with pytest.raises(ValueError, match="CUDA device array"):
        CortexFilamentMechanics(
            device="cuda",
            n_vertices=4,
            links_d=_FakeArray(500, shape=(3, 2), dtype=wp.int32, device=cpu),
            link_k_d=_FakeArray(501, shape=(3,), dtype=wp.float64, device=cpu),
            link_r0_d=_FakeArray(502, shape=(3,), dtype=wp.float64, device=cpu),
            bend_triples_d=_FakeArray(503, shape=(2, 3), dtype=wp.int32, device=cpu),
            bend_alpha_d=_FakeArray(504, shape=(2,), dtype=wp.float64, device=cpu),
            link_kernel=_LINK_KERNEL,
            bending_kernel=_BENDING_KERNEL,
            launch=recorder,
        )


def test_cortex_mechanics_wires_into_state_owner_and_surface_body() -> None:
    recorder = _LaunchRecorder()
    membrane, *_ = _owner("membrane", 100)
    cortex = SurfaceComponentStateOwner(
        name="cortex",
        position_d=_FakeArray(200),
        force_d=_FakeArray(201),
        mechanics=_cortex_filament_mechanics(recorder, n_vertices=4),
        transaction=_TransactionSpy(),
        ledger=_LedgerSpy(),
    )
    body = SurfaceBody(_architecture_with_membrane_fluid_boundary(), membrane, cortex)
    body.cortex.accumulate()

    assert [call["dim"] for call in recorder.calls] == [3, 2]
    # every launch targets the cortex-owned force array, never the membrane's
    for call in recorder.calls:
        assert call["inputs"][0] is body.cortex.position_d
        assert call["inputs"][-1] is body.cortex.force_d


def test_cortex_mechanics_bind_native_wires_the_real_ff_kernels() -> None:
    from aleph.laws.forces_warp import cytosim_bending_kernel
    from aleph.laws.network_warp import link_spring_kernel

    recorder = _LaunchRecorder()
    mechanics = CortexFilamentMechanics.bind_native(
        device="cuda:0",
        n_vertices=4,
        links_d=_FakeArray(500, shape=(3, 2), dtype=wp.int32),
        link_k_d=_FakeArray(501, shape=(3,), dtype=wp.float64),
        link_r0_d=_FakeArray(502, shape=(3,), dtype=wp.float64),
        bend_triples_d=_FakeArray(503, shape=(2, 3), dtype=wp.int32),
        bend_alpha_d=_FakeArray(504, shape=(2,), dtype=wp.float64),
        launch=recorder,
    )
    assert mechanics.link_kernel is link_spring_kernel
    assert mechanics.bending_kernel is cytosim_bending_kernel


# --- kernel-bound membrane pressure fluid boundary ------------------------------------------------


@dataclass(slots=True)
class _PressureTractionSpy:
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


def test_membrane_pressure_boundary_forwards_to_traction_kernel() -> None:
    traction = _PressureTractionSpy()
    boundary = MembranePressureFluidBoundary(traction)
    assert boundary.name == MEMBRANE_FLUID_BOUNDARY
    assert {boundary.component_a, boundary.component_b} == {"membrane", "cytosol"}

    pos = _FakeArray(700)
    force = _FakeArray(701)
    boundary.accumulate_boundary(pos, force)
    assert traction.calls == [(pos, force)]

    # non-kinetic field coupler: transaction hooks are no-ops and never touch the traction primitive
    boundary.snapshot_candidate()
    boundary.rollback(object())
    boundary.commit_irreversible(object(), dt_phys=0.05, rng_seed=7)
    assert traction.calls == [(pos, force)]


def test_membrane_pressure_boundary_serves_as_the_facade_fluid_delegate() -> None:
    body, _bindings, *_ = _body_and_bindings()
    cytosol = CytosolFieldEndpoint(
        "cytosol",
        _FakeArray(400, shape=(2, 2, 2), dtype=wp.float64),
        _FakeArray(401, shape=(2, 2, 2), dtype=wp.float64),
    )
    traction = _PressureTractionSpy()
    bindings = SurfaceStepBindings(
        _PairConnectorSpy(),
        MembranePressureFluidBoundary(traction),
        _PorousTransferSpy(),
        cytosol,
    )
    body.accumulate_mechanics(bindings)
    assert traction.calls == [(body.membrane.position_d, body.membrane.force_d)]


def test_membrane_pressure_boundary_rejects_wrong_endpoints_and_bad_traction() -> None:
    with pytest.raises(ValueError, match="must join membrane and cytosol"):
        MembranePressureFluidBoundary(_PressureTractionSpy(), component_b="cortex")
    with pytest.raises(TypeError, match="accumulate"):
        MembranePressureFluidBoundary(object())
