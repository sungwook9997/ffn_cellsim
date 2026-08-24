"""Structural/runtime gates for the FC-1 Fluid Volume + Core Body seam."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorContract,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine.fluid_core import (
    NUCLEUS_FLUID_BOUNDARY,
    BiotSubstrateFluidSolver,
    FluidCandidateSolver,
    FluidCore,
    FluidCoreStepBindings,
    FluidVolumeStateOwner,
    MovingBoundaryDelegate,
    MovingImpermeableFluidBoundaryFacade,
    NativeCoreLaminaConstants,
    NativeCoreSurfaceStep,
    NativeSurfaceCoreForces,
    NucleusPressureAdjointBoundary,
    ReducedCoreMechanics,
    ReducedCoreSettings,
    ReducedCoreStateOwner,
    SurfaceQuadratureState,
    lumped_surface_control_volumes,
)
from aleph.engine.runtime import (
    CytosolFieldEndpoint,
    LedgerContributor,
    MechanicsContributor,
    TransactionParticipant,
)


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """Metadata-only CUDA array double; no authoritative physics executes in these tests."""

    ptr: int
    shape: tuple[int, ...]
    dtype: object
    device: _FakeDevice = _FakeDevice()


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
class _FluidSolverSpy:
    order: list[str]
    calls: list[float] = field(default_factory=list)

    def solve_candidate(self, dt_phys: float) -> None:
        self.order.append("fluid.solve")
        self.calls.append(dt_phys)


@dataclass(slots=True)
class _CoreMechanicsSpy:
    order: list[str]
    n_deformation_modes: int = 4
    reconstruct_calls: list[tuple[object, object, object]] = field(default_factory=list)
    internal_calls: list[tuple[object, object, object]] = field(default_factory=list)
    project_calls: list[tuple[object, object]] = field(default_factory=list)
    solve_calls: list[tuple[object, ...]] = field(default_factory=list)

    def reconstruct_surface(
        self,
        generalized_position: wp.array,
        surface_position: wp.array,
        surface_velocity: wp.array,
    ) -> None:
        self.order.append("core.reconstruct")
        self.reconstruct_calls.append((generalized_position, surface_position, surface_velocity))

    def accumulate_internal(
        self,
        generalized_position: wp.array,
        volume_multiplier: wp.array,
        generalized_force: wp.array,
    ) -> None:
        self.order.append("core.internal")
        self.internal_calls.append((generalized_position, volume_multiplier, generalized_force))

    def project_surface_force(
        self,
        surface_force: wp.array,
        generalized_force: wp.array,
    ) -> None:
        self.order.append("core.project_jt")
        self.project_calls.append((surface_force, generalized_force))

    def solve_candidate(
        self,
        dt_phys: float,
        generalized_position: wp.array,
        generalized_force: wp.array,
        volume_multiplier: wp.array,
    ) -> None:
        self.order.append("core.solve")
        self.solve_calls.append(
            (dt_phys, generalized_position, generalized_force, volume_multiplier)
        )


@dataclass(slots=True)
class _BoundaryDelegateSpy:
    order: list[str]
    geometry_calls: list[tuple[object, object]] = field(default_factory=list)
    traction_calls: list[tuple[object, object]] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def update_geometry(self, surface_position: wp.array, surface_velocity: wp.array) -> None:
        self.order.append("boundary.geometry")
        self.geometry_calls.append((surface_position, surface_velocity))

    def accumulate_boundary(self, surface_position: wp.array, surface_force: wp.array) -> None:
        self.order.append("boundary.traction")
        self.traction_calls.append((surface_position, surface_force))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


def _fluid_owner(
    order: list[str],
    *,
    pressure_ptr: int = 100,
) -> tuple[FluidVolumeStateOwner, _FluidSolverSpy, _TransactionSpy, _LedgerSpy]:
    solver = _FluidSolverSpy(order)
    transaction = _TransactionSpy()
    ledger = _LedgerSpy()
    owner = FluidVolumeStateOwner(
        name="cytosol",
        pressure_d=_FakeArray(pressure_ptr, (5, 6, 7), wp.float64),
        mask_d=_FakeArray(101, (5, 6, 7), wp.int32),
        solver=solver,
        transaction=transaction,
        ledger=ledger,
    )
    return owner, solver, transaction, ledger


def test_fluid_owner_exposes_live_pressure_and_solver_residual_to_transfer() -> None:
    owner, *_ = _fluid_owner([])
    residual_d = _FakeArray(102, owner.pressure_d.shape, wp.float64)
    endpoint = owner.make_immersed_transfer_endpoint(residual_d)

    assert isinstance(endpoint, CytosolFieldEndpoint)
    assert endpoint.pressure_d is owner.pressure_d
    assert endpoint.residual_d is residual_d

    with pytest.raises(ValueError, match="must not alias"):
        owner.make_immersed_transfer_endpoint(owner.pressure_d)


def _core_owner(
    order: list[str],
    *,
    q_ptr: int = 200,
    n_modes: int = 4,
) -> tuple[ReducedCoreStateOwner, _CoreMechanicsSpy, _TransactionSpy, _LedgerSpy]:
    settings = ReducedCoreSettings(
        n_deformation_modes=n_modes,
        reference_volume_um3=500.0,
    )
    mechanics = _CoreMechanicsSpy(order, n_deformation_modes=n_modes)
    transaction = _TransactionSpy()
    ledger = _LedgerSpy()
    owner = ReducedCoreStateOwner(
        name="nucleus",
        settings=settings,
        generalized_position_d=_FakeArray(
            q_ptr,
            (settings.n_generalized_dofs,),
            wp.float64,
        ),
        generalized_force_d=_FakeArray(201, (settings.n_generalized_dofs,), wp.float64),
        volume_multiplier_d=_FakeArray(202, (1,), wp.float64),
        surface_position_d=_FakeArray(203, (42,), wp.vec3d),
        surface_velocity_d=_FakeArray(204, (42,), wp.vec3d),
        surface_force_d=_FakeArray(205, (42,), wp.vec3d),
        mechanics=mechanics,
        transaction=transaction,
        ledger=ledger,
    )
    return owner, mechanics, transaction, ledger


def _fluid_core() -> tuple[
    FluidCore,
    FluidCoreStepBindings,
    list[str],
    tuple[_FluidSolverSpy, _TransactionSpy, _LedgerSpy],
    tuple[_CoreMechanicsSpy, _TransactionSpy, _LedgerSpy],
    _BoundaryDelegateSpy,
]:
    order: list[str] = []
    fluid, fluid_solver, fluid_tx, fluid_ledger = _fluid_owner(order)
    core, core_mechanics, core_tx, core_ledger = _core_owner(order)
    boundary_delegate = _BoundaryDelegateSpy(order)
    boundary = MovingImpermeableFluidBoundaryFacade(boundary_delegate)
    facade = FluidCore(reference_cell_architecture(), fluid, core)
    bindings = FluidCoreStepBindings(boundary)
    return (
        facade,
        bindings,
        order,
        (fluid_solver, fluid_tx, fluid_ledger),
        (core_mechanics, core_tx, core_ledger),
        boundary_delegate,
    )


def _architecture_with_boundary(connector: ConnectorContract) -> CellArchitecture:
    base = reference_cell_architecture()
    connectors = tuple(
        connector if item.name == NUCLEUS_FLUID_BOUNDARY else item for item in base.connectors
    )
    return CellArchitecture(base.components, connectors)


def test_reference_graph_builds_distinct_fluid_and_core_state_owners() -> None:
    facade, bindings, *_ = _fluid_core()
    contract = next(
        item for item in facade.architecture.connectors if item.name == NUCLEUS_FLUID_BOUNDARY
    )
    assert facade.cytosol is not facade.nucleus
    assert facade.cytosol.pressure_d.ptr != facade.nucleus.generalized_position_d.ptr
    assert {contract.component_a, contract.component_b} == {"cytosol", "nucleus"}
    assert contract.family is ConnectorFamily.FLUID_BOUNDARY
    assert bindings.nucleus_fluid_boundary.impermeable
    assert bindings.nucleus_fluid_boundary.adjoint_transfer_required


def test_reduced_core_requires_modes_and_an_exact_volume_constraint() -> None:
    settings = ReducedCoreSettings(24, 510.0)
    assert settings.n_generalized_dofs == 30
    assert settings.exact_volume_constraint

    with pytest.raises(ValueError, match="positive integer"):
        ReducedCoreSettings(0, 510.0)
    with pytest.raises(ValueError, match="finite and positive"):
        ReducedCoreSettings(4, 0.0)
    with pytest.raises(ValueError, match="exact nuclear-volume"):
        ReducedCoreSettings(4, 510.0, exact_volume_constraint=False)
    with pytest.raises(ValueError, match="volume_penalty is forbidden"):
        ReducedCoreSettings(4, 510.0, volume_penalty=1.0e3)


def test_core_state_shape_tracks_six_pose_plus_retained_modes() -> None:
    order: list[str] = []
    settings = ReducedCoreSettings(4, 500.0)
    with pytest.raises(ValueError, match=r"6 \+ n_deformation_modes"):
        ReducedCoreStateOwner(
            name="nucleus",
            settings=settings,
            generalized_position_d=_FakeArray(200, (9,), wp.float64),
            generalized_force_d=_FakeArray(201, (10,), wp.float64),
            volume_multiplier_d=_FakeArray(202, (1,), wp.float64),
            surface_position_d=_FakeArray(203, (42,), wp.vec3d),
            surface_velocity_d=_FakeArray(204, (42,), wp.vec3d),
            surface_force_d=_FakeArray(205, (42,), wp.vec3d),
            mechanics=_CoreMechanicsSpy(order),
            transaction=_TransactionSpy(),
            ledger=_LedgerSpy(),
        )


def test_fluid_and_core_reject_cpu_or_aliased_state() -> None:
    order: list[str] = []
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        FluidVolumeStateOwner(
            "cytosol",
            _FakeArray(100, (5, 6, 7), wp.float64, cpu),
            _FakeArray(101, (5, 6, 7), wp.int32, cpu),
            _FluidSolverSpy(order),
            _TransactionSpy(),
            _LedgerSpy(),
        )

    fluid, *_ = _fluid_owner(order)
    core, *_ = _core_owner(order, q_ptr=100)
    with pytest.raises(ValueError, match="disjoint state storage"):
        FluidCore(reference_cell_architecture(), fluid, core)


def test_boundary_facade_is_bidirectional_and_satisfies_common_protocols() -> None:
    delegate = _BoundaryDelegateSpy([])
    boundary = MovingImpermeableFluidBoundaryFacade(delegate)
    assert isinstance(boundary, MechanicsContributor)
    assert isinstance(boundary, TransactionParticipant)
    assert isinstance(boundary, LedgerContributor)

    with pytest.raises(ValueError, match="relatively impermeable"):
        MovingImpermeableFluidBoundaryFacade(delegate, impermeable=False)
    with pytest.raises(ValueError, match="adjoint transfer"):
        MovingImpermeableFluidBoundaryFacade(
            delegate,
            adjoint_transfer_required=False,
        )


def test_state_owners_satisfy_common_transaction_and_ledger_protocols() -> None:
    facade, *_ = _fluid_core()
    assert isinstance(facade.cytosol, TransactionParticipant)
    assert isinstance(facade.cytosol, LedgerContributor)
    assert isinstance(facade.nucleus, TransactionParticipant)
    assert isinstance(facade.nucleus, LedgerContributor)


def test_candidate_iteration_closes_forward_geometry_and_adjoint_force_hooks() -> None:
    facade, bindings, order, fluid_spies, core_spies, boundary = _fluid_core()
    facade.candidate_iteration(bindings, dt_phys=0.05)

    assert order == [
        "core.reconstruct",
        "boundary.geometry",
        "fluid.solve",
        "core.internal",
        "boundary.traction",
        "core.project_jt",
        "core.solve",
        "core.reconstruct",
        "boundary.geometry",
    ]
    assert fluid_spies[0].calls == [0.05]
    assert boundary.geometry_calls == [
        (facade.nucleus.surface_position_d, facade.nucleus.surface_velocity_d),
        (facade.nucleus.surface_position_d, facade.nucleus.surface_velocity_d),
    ]
    assert boundary.traction_calls == [
        (facade.nucleus.surface_position_d, facade.nucleus.surface_force_d)
    ]
    assert core_spies[0].project_calls == [
        (facade.nucleus.surface_force_d, facade.nucleus.generalized_force_d)
    ]


@pytest.mark.parametrize("dt_phys", [0.0, -0.1, float("nan"), float("inf")])
def test_candidate_iteration_rejects_invalid_physical_time_before_mutation(
    dt_phys: float,
) -> None:
    facade, bindings, order, *_ = _fluid_core()
    with pytest.raises(ValueError, match="finite and positive"):
        facade.candidate_iteration(bindings, dt_phys)
    assert order == []


def test_transaction_covers_fluid_core_and_nonkinetic_boundary_caches() -> None:
    facade, bindings, _, fluid_spies, core_spies, boundary = _fluid_core()
    accepted_d = object()

    facade.snapshot_candidate(bindings)
    facade.rollback(bindings, accepted_d)
    facade.commit_irreversible(bindings, accepted_d, dt_phys=0.05, rng_seed=17)

    expected = [
        ("snapshot",),
        ("rollback", accepted_d),
        ("commit", accepted_d, 0.05, 17),
    ]
    assert fluid_spies[1].calls == expected
    assert core_spies[1].calls == expected
    assert boundary.transaction_calls == expected


def test_component_boundary_ledgers_and_dof_ledger_are_complete() -> None:
    facade, bindings, _, fluid_spies, core_spies, boundary = _fluid_core()
    ledger = object()
    facade.accumulate_ledger(bindings, ledger)

    assert fluid_spies[2].calls == [ledger]
    assert core_spies[2].calls == [ledger]
    assert boundary.ledger_calls == [ledger]
    assert facade.dof_ledger() == replace(
        facade.dof_ledger(),
        allocated_pressure_cells=5 * 6 * 7,
        generalized_core_dofs=10,
        deformation_modes=4,
        surface_vertices=42,
        volume_constraint_dofs=1,
        exact_volume_constraint=True,
    )


def test_architecture_rejects_nonboundary_or_nonadjoint_nucleus_fluid_edge() -> None:
    order: list[str] = []
    fluid, *_ = _fluid_owner(order)
    core, *_ = _core_owner(order)
    wrong_family = ConnectorContract(
        NUCLEUS_FLUID_BOUNDARY,
        ConnectorFamily.LINC,
        "nucleus",
        "cytosol",
        kinetics=False,
        commit_on_accept=False,
    )
    with pytest.raises(ValueError, match="FLUID_BOUNDARY"):
        FluidCore(_architecture_with_boundary(wrong_family), fluid, core)

    legacy_bridge = ConnectorContract(
        NUCLEUS_FLUID_BOUNDARY,
        ConnectorFamily.IMMERSED_TRANSFER,
        "nucleus",
        "cytosol",
        kinetics=False,
        commit_on_accept=False,
    )
    with pytest.raises(ValueError, match="FLUID_BOUNDARY"):
        FluidCore(_architecture_with_boundary(legacy_bridge), fluid, core)

    nonadjoint = ConnectorContract(
        NUCLEUS_FLUID_BOUNDARY,
        ConnectorFamily.FLUID_BOUNDARY,
        "nucleus",
        "cytosol",
        kinetics=False,
        commit_on_accept=False,
        adjoint_transfer_required=False,
    )
    with pytest.raises(ValueError, match="bidirectional with adjoint"):
        FluidCore(_architecture_with_boundary(nonadjoint), fluid, core)


# =====================================================================================================
# KERNEL_BOUND slice — the real ac/fluid + ac/nucleus backends injected behind the delegate Protocols.
# Host-side wiring/validation/quadrature are CPU-green via injected doubles (no real kernel launch on the
# dev Mac's CPU-only Warp); the end-to-end one-candidate closure runs only under CUDA.
# =====================================================================================================

import numpy as np  # noqa: E402


@dataclass(slots=True)
class _StepSpyGrid:
    p: object
    mask: object
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _SubstrateSpy:
    grid: _StepSpyGrid
    steps: list[float] = field(default_factory=list)

    def step(self, dt: float) -> None:
        self.steps.append(dt)


@dataclass(slots=True)
class _DomainSpy:
    content: object
    remaps: int = 0

    def remap(self) -> object:
        self.remaps += 1
        return self.content


@dataclass(slots=True)
class _PressureCouplingSpy:
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, state: object, out_force: object) -> None:
        self.calls.append((state, out_force))


@dataclass(slots=True)
class _NoFluxSpy:
    value: float = 3.5

    def suppressed_flux(self) -> float:
        return self.value


@dataclass(slots=True)
class _NativeForcesSpy:
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, surface_pos: object, surface_force: object) -> None:
        self.calls.append((surface_pos, surface_force))


def _biot_solver_spy() -> tuple[BiotSubstrateFluidSolver, _SubstrateSpy]:
    grid = _StepSpyGrid(p=_FakeArray(400, (5, 6, 7), wp.float64), mask=_FakeArray(401, (5, 6, 7), wp.int32))
    substrate = _SubstrateSpy(grid)
    return BiotSubstrateFluidSolver(substrate), substrate


def _adjoint_boundary_spy() -> tuple[
    NucleusPressureAdjointBoundary, _DomainSpy, _PressureCouplingSpy, _NoFluxSpy, _TransactionSpy, _LedgerSpy
]:
    content = _FakeArray(500, (1,), wp.float64)
    domain = _DomainSpy(content)
    coupling = _PressureCouplingSpy()
    no_flux = _NoFluxSpy()
    transaction = _TransactionSpy()
    ledger = _LedgerSpy()
    node_volume = _FakeArray(501, (42,), wp.float64)
    boundary = NucleusPressureAdjointBoundary(
        domain, coupling, no_flux, node_volume, transaction, ledger
    )
    return boundary, domain, coupling, no_flux, transaction, ledger


def test_lumped_surface_control_volumes_partition_of_unity() -> None:
    from aleph.components.nucleus.geometry import build_oblate_mesh, face_normals_areas

    verts, faces = build_oblate_mesh(2.5, 1.0, 2, (0.0, 0.0, 0.0))
    dx = 0.5
    node_volume = lumped_surface_control_volumes(verts, faces, dx)

    assert node_volume.shape == (verts.shape[0],)
    assert np.all(node_volume > 0.0)
    _, areas = face_normals_areas(verts, faces)
    # Each triangle's area is split equally across its 3 nodes -> sum == total area * shell thickness.
    np.testing.assert_allclose(node_volume.sum(), float(areas.sum()) * dx, rtol=1.0e-12)

    with pytest.raises(ValueError, match="finite and positive"):
        lumped_surface_control_volumes(verts, faces, 0.0)


def test_biot_fluid_solver_binds_step_and_exposes_public_field() -> None:
    solver, substrate = _biot_solver_spy()

    assert isinstance(solver, FluidCandidateSolver)
    assert solver.pressure_d is substrate.grid.p
    assert solver.mask_d is substrate.grid.mask

    solver.solve_candidate(0.05)
    assert substrate.steps == [0.05]

    for bad in (0.0, -0.1, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="finite and positive"):
            solver.solve_candidate(bad)
    assert substrate.steps == [0.05]  # no candidate mutation on invalid dt


def test_biot_fluid_solver_rejects_non_cuda_grid_and_bad_backend() -> None:
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    cpu_grid = _StepSpyGrid(
        p=_FakeArray(400, (4, 4, 4), wp.float64, cpu),
        mask=_FakeArray(401, (4, 4, 4), wp.int32, cpu),
        device=cpu,
    )
    with pytest.raises(ValueError, match="CUDA device"):
        BiotSubstrateFluidSolver(_SubstrateSpy(cpu_grid))
    with pytest.raises(TypeError, match="BiotSubstrate"):
        BiotSubstrateFluidSolver(object())


def test_adjoint_boundary_routes_remap_pressure_traction_and_no_flux() -> None:
    boundary, domain, coupling, no_flux, transaction, ledger = _adjoint_boundary_spy()

    assert isinstance(boundary, MovingBoundaryDelegate)

    pos = _FakeArray(600, (42,), wp.vec3d)
    vel = _FakeArray(601, (42,), wp.vec3d)
    force = _FakeArray(602, (42,), wp.vec3d)

    # Moving-face mass channel is produced by the real Domain.remap and cached.
    assert boundary.moving_face_content() is None
    boundary.update_geometry(pos, vel)
    assert domain.remaps == 1
    assert boundary.moving_face_content() is domain.content

    # Pressure traction is scattered through PressureCoupling onto the nucleus-owned surface force,
    # using the boundary-owned node control volumes (the adjoint of the boundary-velocity map).
    boundary.accumulate_boundary(pos, force)
    assert len(coupling.calls) == 1
    state, out_force = coupling.calls[0]
    assert isinstance(state, SurfaceQuadratureState)
    assert state.node_pos is pos
    assert state.node_volume is boundary.node_volume_d
    assert out_force is force

    # No-flux "teeth" probe forwards to the real suppressed-flux kernel wrapper.
    assert boundary.no_flux_teeth() == pytest.approx(3.5)

    # Transaction + ledger forward to the scheduler-backed participants (cache rollback is not exempt).
    accepted = object()
    boundary.snapshot_candidate()
    boundary.rollback(accepted)
    boundary.commit_irreversible(accepted, 0.05, 11)
    boundary.accumulate_ledger("L")
    assert transaction.calls == [
        ("snapshot",),
        ("rollback", accepted),
        ("commit", accepted, 0.05, 11),
    ]
    assert ledger.calls == ["L"]


def test_adjoint_boundary_validates_node_volume_and_dependencies() -> None:
    domain = _DomainSpy(_FakeArray(500, (1,), wp.float64))
    coupling = _PressureCouplingSpy()
    no_flux = _NoFluxSpy()
    good_nv = _FakeArray(501, (42,), wp.float64)

    with pytest.raises(TypeError, match="node_volume_d must have dtype"):
        NucleusPressureAdjointBoundary(
            domain, coupling, no_flux, _FakeArray(501, (42,), wp.int32), _TransactionSpy(), _LedgerSpy()
        )
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        NucleusPressureAdjointBoundary(
            domain, coupling, no_flux, _FakeArray(501, (42,), wp.float64, cpu),
            _TransactionSpy(), _LedgerSpy(),
        )
    with pytest.raises(TypeError, match="no_flux_bc must provide"):
        NucleusPressureAdjointBoundary(domain, coupling, object(), good_nv, _TransactionSpy(), _LedgerSpy())


def test_surface_quadrature_state_requires_matching_vec3_and_float64_ports() -> None:
    pos = _FakeArray(600, (42,), wp.vec3d)
    vol = _FakeArray(601, (42,), wp.float64)
    state = SurfaceQuadratureState(pos, vol)
    assert state.node_pos is pos and state.node_volume is vol

    with pytest.raises(TypeError, match="node_pos must have dtype"):
        SurfaceQuadratureState(_FakeArray(600, (42,), wp.float64), vol)
    with pytest.raises(ValueError, match="matching shape"):
        SurfaceQuadratureState(pos, _FakeArray(601, (7,), wp.float64))


def test_native_core_lamina_constants_reject_unsourced_or_inverted_regimes() -> None:
    ok = NativeCoreLaminaConstants(
        k_soft_pn_per_um=6.0,
        k_ac_pn_per_um=150.0,
        knee_strain=0.10,
        eps_rupture=0.30,
        k_vol_pn_per_um2=399.0,
        reference_volume_um3=520.0,
        k_linc_pn_per_um=8.0,
    )
    assert ok.linc_stiffening_per_um2 == 0.0

    with pytest.raises(ValueError, match="knee strain must precede"):
        replace(ok, knee_strain=0.30, eps_rupture=0.30)
    with pytest.raises(ValueError, match="k_soft_pn_per_um must be finite and positive"):
        replace(ok, k_soft_pn_per_um=0.0)
    with pytest.raises(ValueError, match="k_linc_pn_per_um must be finite and nonnegative"):
        replace(ok, k_linc_pn_per_um=-1.0)


def test_native_surface_core_and_step_bind_backend_and_reject_bad_arrays() -> None:
    forces = _NativeForcesSpy()
    step = NativeCoreSurfaceStep(forces)
    pos = _FakeArray(700, (42,), wp.vec3d)
    force = _FakeArray(701, (42,), wp.vec3d)
    step.accumulate_internal(pos, force)
    assert forces.calls == [(pos, force)]

    with pytest.raises(ValueError, match="finite and positive"):
        step.candidate(0.0, pos, force)
    with pytest.raises(TypeError, match="NativeSurfaceCoreForces backend"):
        NativeCoreSurfaceStep(object())

    # The real force binder validates its device-array schema without launching kernels.
    with pytest.raises(TypeError, match="faces_d must be a rank-2"):
        NativeSurfaceCoreForces(
            _FakeArray(710, (42,), wp.int32),  # wrong rank
            _FakeArray(711, (80,), wp.float64),
            _FakeArray(712, (80,), wp.int32),
            NativeCoreLaminaConstants(6.0, 150.0, 0.1, 0.3, 399.0, 520.0, 8.0),
            anchor_pos_d=_FakeArray(713, (1,), wp.vec3d),
            anchor_force_d=_FakeArray(714, (1,), wp.vec3d),
            linc_nucleus_idx_d=_FakeArray(715, (0,), wp.int32),
            linc_anchor_idx_d=_FakeArray(716, (0,), wp.int32),
            linc_rest_d=_FakeArray(717, (0,), wp.float64),
        )


def test_native_step_is_not_a_reduced_core_mechanics_delegate() -> None:
    # The native full-surface backend deliberately does NOT satisfy the reduced Core Body delegate:
    # its vec3d surface force cannot flow through the reduced float64 generalized force without the ROM.
    step = NativeCoreSurfaceStep(_NativeForcesSpy())
    assert not isinstance(step, ReducedCoreMechanics)


_CUDA_DEVICE = next((device for device in wp.get_devices() if device.is_cuda), None)


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="fluid_core KERNEL_BOUND gate requires CUDA")
def test_cuda_one_candidate_closes_mass_no_flux_and_adjoint_work() -> None:
    """One native candidate: Biot solve -> pressure adjoint traction -> native Core forces -> overdamped
    step -> remap, asserting (a) moving-face mass channel is finite, (b) the relative-no-flux mask has
    teeth, and (c) the fluid->core pressure traction is a nonzero, finite adjoint back-reaction."""
    from aleph.components.fluid.biot_substrate import BiotSubstrate, PressureCoupling
    from aleph.components.fluid.boundary import NucleusNoFluxBC
    from aleph.components.fluid.domain import (
        Domain,
        StaticSphereMembraneProvider,
        StaticSphereNucleusMaskProvider,
    )
    from aleph.components.fluid.field_grid import FieldGrid
    from aleph.components.nucleus.envelope import build_nucleus
    from aleph.components.nucleus.lamina_analytic import LaminaParams

    device = str(_CUDA_DEVICE)
    with wp.ScopedDevice(device):
        # --- Cytosol field (sourced I0-B1 closure: c_v = mobility / S) --------------------------------
        dx = 0.5
        n = 24
        origin = (-(n * dx) / 2.0,) * 3
        grid = FieldGrid((n, n, n), dx, origin, device=device)
        c_v = 50.0                      # um^2/s (KB-3.B3.2, Moeendarbary 2013)
        storage_S = 1.0                 # normalized storativity; c_v fixes mobility
        mobility = c_v * storage_S
        substrate = BiotSubstrate(grid, mobility=mobility, storage_S=storage_S, alpha=1.0)

        membrane = StaticSphereMembraneProvider((0.0, 0.0, 0.0), 5.0)
        nucleus_mask = StaticSphereNucleusMaskProvider((0.0, 0.0, 0.0), 2.5)
        domain = Domain(grid, membrane, nucleus_mask)
        domain.bind_storage(storage_S)
        domain.classify()

        # Seed a spatial pressure so grad(p) and the cross-envelope discharge are non-uniform.
        p_host = np.zeros((n, n, n), dtype=np.float64)
        xs = origin[0] + dx * np.arange(n)
        p_host += (xs[:, None, None] - origin[0])  # linear ramp in x [Pa]
        grid.set_pressure(p_host)

        fluid_solver = BiotSubstrateFluidSolver(substrate)
        pressure_coupling = PressureCoupling(substrate)
        no_flux = NucleusNoFluxBC(grid, mobility=mobility)

        # --- Native Core Body surface (full lamina/chromatin backend, I0-B2 sourced constants) --------
        params = LaminaParams(
            k_chrom=3.0, k_lamin_b=3.0, k_lamin_ac=150.0,
            knee_strain=0.10, eps_rupture=0.30, kappa_ne=0.0828,
        )
        mesh = build_nucleus(2.5, params, aspect=1.0, subdivisions=2)
        nv = mesh.verts.shape[0]
        surface_pos = wp.array(mesh.verts, dtype=wp.vec3d, device=device)
        surface_force = wp.zeros(nv, dtype=wp.vec3d, device=device)
        faces_d = wp.array(np.ascontiguousarray(mesh.faces, np.int32), dtype=wp.int32, ndim=2, device=device)
        a0_d = wp.array(mesh.A0_face, dtype=wp.float64, device=device)
        ruptured_d = wp.zeros(mesh.faces.shape[0], dtype=wp.int32, device=device)
        node_volume_d = wp.array(
            lumped_surface_control_volumes(mesh.verts, mesh.faces, dx), dtype=wp.float64, device=device
        )
        constants = NativeCoreLaminaConstants(
            k_soft_pn_per_um=params.k_chrom + params.k_lamin_b,
            k_ac_pn_per_um=params.k_lamin_ac,
            knee_strain=params.knee_strain,
            eps_rupture=params.eps_rupture,
            k_vol_pn_per_um2=399.0,           # E_nuc MCF7 in-situ (Fischer 2020, KB-DRAFT-3.B-05)
            reference_volume_um3=mesh.V0,
            k_linc_pn_per_um=8.0,             # LINC tension ~8 pN (KB-DRAFT-3.B-09; magnitude a GAP)
        )
        anchor_pos_d = wp.array(np.array([[0.0, 0.0, 6.0]]), dtype=wp.vec3d, device=device)
        anchor_force_d = wp.zeros(1, dtype=wp.vec3d, device=device)
        forces = NativeSurfaceCoreForces(
            faces_d, a0_d, ruptured_d, constants,
            anchor_pos_d=anchor_pos_d, anchor_force_d=anchor_force_d,
            linc_nucleus_idx_d=wp.zeros(0, dtype=wp.int32, device=device),
            linc_anchor_idx_d=wp.zeros(0, dtype=wp.int32, device=device),
            linc_rest_d=wp.zeros(0, dtype=wp.float64, device=device),
        )
        core_step = NativeCoreSurfaceStep(forces)

        class _NullTx:
            def snapshot_candidate(self): ...
            def rollback(self, accepted): ...
            def commit_irreversible(self, accepted, dt_phys, rng_seed): ...

        class _NullLedger:
            def accumulate_ledger(self, ledger): ...

        boundary = NucleusPressureAdjointBoundary(
            domain, pressure_coupling, no_flux, node_volume_d, _NullTx(), _NullLedger()
        )

        # --- ONE native candidate iteration -----------------------------------------------------------
        dt = 1.0e-4
        teeth_before = boundary.no_flux_teeth()
        fluid_solver.solve_candidate(dt)                       # real biot_pmass_update_kernel
        surface_force.zero_()
        core_step.accumulate_internal(surface_pos, surface_force)   # real lamina/volume/LINC kernels
        f_internal = surface_force.numpy().copy()
        boundary.accumulate_boundary(surface_pos, surface_force)    # real PressureCoupling adjoint traction
        f_total = surface_force.numpy()
        core_step.candidate(dt, surface_pos, surface_force)         # overdamped native surface candidate
        boundary.update_geometry(surface_pos, None)                # real Domain.remap moving-face mass
        content = boundary.moving_face_content().numpy()[0]

    # (a) mass: the moving-face channel is a finite conservative transport scalar.
    assert np.isfinite(content)
    # (b) no-flux teeth: a transmembrane gradient exists, so the mask is actively holding flux at zero.
    assert teeth_before > 0.0
    # (c) adjoint work: the fluid pressure traction is a nonzero, finite back-reaction distinct from the
    #     internal force alone (Newton-3rd Peskin scatter into the nucleus-owned surface force).
    assert np.all(np.isfinite(f_total))
    assert np.linalg.norm(f_total - f_internal) > 0.0
