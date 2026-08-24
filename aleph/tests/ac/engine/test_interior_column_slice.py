"""CPU-green structural gates for the GATE-B interior-column slice (membrane → cytosol → nucleus).

Proves — with recording CUDA doubles + injected fluid-primitive doubles, no device — that
:class:`~aleph.engine.interior_column_slice.InteriorColumnSlice` drives the accepted-step transaction so
the four interior compartments (membrane, cortex, cytosol, nucleus) run TOGETHER as one candidate under one
device predicate, wired only through the three fluid connectors:

* the four state-owners + three fluid connectors are ALL transaction participants (7), and NONE proposes an
  event (every interior coupling is a non-kinetic field coupler; the one kinetic is nucleus rupture, committed
  in ``commit_irreversible``);
* the cytosol field participant snapshots + reject-restores the Biot field through the REAL
  ``ac/fluid/scheduler`` restore kernels, and the nucleus participant commits the REAL
  ``ac/nucleus.conditional_rupture_update_kernel`` accepted-gated (a rejected step reject-restores surface +
  ruptured and commits no rupture);
* the coupling connectors are bidirectional + adjoint (Newton's 3rd law); the porous transfer scatters the
  cytosol pressure traction onto the CORTEX force array through the injected ``PressureCoupling``;
* a rejected step bit-restores every reversible interior array and advances no clock; the clock advances only
  on accept;
* the REQUIRED-PARAM guard fires — the nucleus rupture strain (I0-B2 GAP) is supplied, never invented.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest
import warp as wp

from aleph.engine.cytosol_connected import MovingSemipermeableFluidBoundaryFacade
from aleph.engine.afm_cortical_tension import (
    AFM_INDENTER_COMPONENT,
    AFM_MEMBRANE_CONTACT,
    afm_experiment_architecture,
)
from aleph.engine.erm_cortex_connector import _restore_int32_if_rejected_kernel
from aleph.engine.erm_cortex_slice import _restore_vec3_if_rejected_kernel
from aleph.engine.fluid_core import (
    FluidVolumeStateOwner,
    MovingImpermeableFluidBoundaryFacade,
)
from aleph.engine.surface_body import SurfaceComponentStateOwner
from aleph.engine.interior_column_slice import (
    CoreSurfaceStateOwner,
    CoreSurfaceTransaction,
    FluidFieldTransaction,
    NoOpSurfaceMechanics,
    RecordingLedgerStub,
    build_cortex_cytosol_transfer,
    build_interior_column_slice,
)
from aleph.engine.runtime import CytosolFieldEndpoint
from aleph.components.fluid.scheduler import (
    _conditional_restore_f64_3d_kernel,
    _conditional_restore_i32_3d_kernel,
)
from aleph.components.nucleus.envelope import conditional_rupture_update_kernel


# ── CUDA metadata doubles (no host physics is evaluated by the structural gates) ──────────────────────
@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    ptr: int
    shape: tuple[int, ...]
    dtype: object
    device: _FakeDevice = _FakeDevice()


def _fa(ptr: int, shape, dtype: object) -> _FakeArray:
    return _FakeArray(ptr, tuple(shape) if isinstance(shape, tuple) else (shape,), dtype)


@dataclass(slots=True)
class _FakeAccepted:
    """Host-readable acceptance predicate double; the device kernels gate in-device, the spy clock in Python."""

    accept: bool


@dataclass(slots=True)
class _RecordingLauncher:
    events: list = field(default_factory=list)
    calls: list = field(default_factory=list)

    def __call__(self, kernel, *, dim, inputs, outputs=None, device=None) -> None:
        self.calls.append((kernel, dim, tuple(inputs), tuple(outputs or ()), device))
        self.events.append(("launch", kernel))

    @property
    def kernels(self) -> list:
        return [c[0] for c in self.calls]


@dataclass(slots=True)
class _RecordingCopy:
    events: list = field(default_factory=list)
    calls: list = field(default_factory=list)

    def __call__(self, dst, src) -> None:
        self.calls.append((dst, src))
        self.events.append(("copy", dst, src))


@dataclass(slots=True)
class _GatingSpyClock:
    """Event-clock double that GATES the epoch on a host-readable predicate (mirrors the in-kernel gate)."""

    events: list = field(default_factory=list)
    base_seed: int = 5
    epoch: int = 0
    advances: list = field(default_factory=list)

    def advance(self, accepted, dt_phys) -> None:
        self.advances.append((accepted, dt_phys))
        self.events.append(("advance", accepted, dt_phys))
        if getattr(accepted, "accept", False):
            self.epoch += 1


class _FakeGrid:
    """A live Biot ``FieldGrid`` double exposing the reversible field arrays + shape (no device physics)."""

    def __init__(self) -> None:
        self.shape = (4, 4, 4)
        self.p = _fa(100, self.shape, wp.float64)
        self.p_new = _fa(101, self.shape, wp.float64)
        self.mask = _fa(102, self.shape, wp.int32)
        self.div_vs = _fa(103, self.shape, wp.float64)
        self.s_membrane = _fa(104, self.shape, wp.float64)
        self.s_total = _fa(105, self.shape, wp.float64)


@dataclass(slots=True)
class _MechanicsSpy:
    calls: list = field(default_factory=list)

    def accumulate(self, pos, force) -> None:
        self.calls.append((pos, force))


@dataclass(slots=True)
class _TransactionSpy:
    calls: list = field(default_factory=list)

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted, dt_phys, rng_seed) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))


@dataclass(slots=True)
class _LedgerSpy:
    calls: list = field(default_factory=list)

    def accumulate_ledger(self, ledger) -> None:
        self.calls.append(ledger)


def _surface_owner(name: str, pos: _FakeArray, force: _FakeArray) -> SurfaceComponentStateOwner:
    return SurfaceComponentStateOwner(
        name=name, position_d=pos, force_d=force,
        mechanics=_MechanicsSpy(), transaction=_TransactionSpy(), ledger=_LedgerSpy())


class _SolverSpy:
    def __init__(self) -> None:
        self.calls: list = []

    def solve_candidate(self, dt_phys) -> None:
        self.calls.append(dt_phys)


class _PressureCouplingSpy:
    """Injected ``PressureCoupling`` double: records the (state, force) it would scatter onto."""

    def __init__(self) -> None:
        self.calls: list = []

    def accumulate(self, state, force) -> None:
        self.calls.append((state, force))


@dataclass(slots=True)
class _BoundaryDelegateSpy:
    """A moving-boundary delegate double (nucleus/membrane) recording every hook in order."""

    label: str
    calls: list = field(default_factory=list)

    def update_geometry(self, pos, vel) -> None:
        self.calls.append((self.label, "update_geometry"))

    def accumulate_boundary(self, pos, force) -> None:
        self.calls.append((self.label, "accumulate"))

    def snapshot_candidate(self) -> None:
        self.calls.append((self.label, "snapshot"))

    def rollback(self, accepted) -> None:
        self.calls.append((self.label, "rollback", accepted))

    def commit_irreversible(self, accepted, dt_phys, rng_seed) -> None:
        self.calls.append((self.label, "commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger) -> None:
        self.calls.append((self.label, "ledger", ledger))


def _fluid_owner(grid: _FakeGrid, solver, launch, copy) -> FluidVolumeStateOwner:
    """A CPU-safe ``cytosol`` owner: the real FluidFieldTransaction over fake field arrays (no wp.zeros)."""
    shape = grid.shape
    txn = FluidFieldTransaction(
        grid=grid, snap_p_d=_fa(200, shape, wp.float64), snap_mask_d=_fa(201, shape, wp.int32),
        snap_div_vs_d=_fa(202, shape, wp.float64), snap_s_membrane_d=_fa(203, shape, wp.float64),
        snap_s_total_d=_fa(204, shape, wp.float64), shape=shape, launch=launch, copy=copy)
    return FluidVolumeStateOwner(
        name="cytosol", pressure_d=grid.p, mask_d=grid.mask, solver=solver, transaction=txn,
        ledger=RecordingLedgerStub())


# ── builders: a real interior-column slice assembled from doubles (CPU-safe, no device alloc) ─────────
def _slice(*, eps_rupture: float = 0.5, apparatus=None):
    events: list = []
    launch = _RecordingLauncher(events=events)
    copy = _RecordingCopy(events=events)
    clock = _GatingSpyClock(events=events)

    mem_pos, mem_f = _fa(10, 4, wp.vec3d), _fa(11, 4, wp.vec3d)
    ctx_pos, ctx_f = _fa(20, 4, wp.vec3d), _fa(21, 4, wp.vec3d)
    membrane = _surface_owner("membrane", mem_pos, mem_f)
    cortex = _surface_owner("cortex", ctx_pos, ctx_f)

    grid = _FakeGrid()
    solver = _SolverSpy()
    cytosol = _fluid_owner(grid, solver, launch, copy)

    nuc_pos, nuc_f = _fa(30, 6, wp.vec3d), _fa(31, 6, wp.vec3d)
    faces = _fa(32, (3, 3), wp.int32)
    a0 = _fa(33, 3, wp.float64)
    ruptured = _fa(34, 3, wp.int32)
    nuc_txn = CoreSurfaceTransaction(
        position_d=nuc_pos, position_snap_d=_fa(35, 6, wp.vec3d), faces_d=faces, a0_d=a0,
        ruptured_d=ruptured, ruptured_snap_d=_fa(36, 3, wp.int32), eps_rupture=eps_rupture,
        launch=launch, copy=copy)
    nucleus = CoreSurfaceStateOwner(
        name="nucleus", position_d=nuc_pos, force_d=nuc_f, mechanics=NoOpSurfaceMechanics(),
        transaction=nuc_txn, ledger=RecordingLedgerStub())

    pc = _PressureCouplingSpy()
    node_vol = _fa(40, 4, wp.float64)
    cortex_transfer = build_cortex_cytosol_transfer(pressure_coupling=pc, node_volume_d=node_vol)
    membrane_boundary = MovingSemipermeableFluidBoundaryFacade(delegate=_BoundaryDelegateSpy("membrane_b"))
    nucleus_boundary = MovingImpermeableFluidBoundaryFacade(delegate=_BoundaryDelegateSpy("nucleus_b"))

    field_p, field_r = _fa(50, (4, 4, 4), wp.float64), _fa(51, (4, 4, 4), wp.float64)
    endpoint = CytosolFieldEndpoint("cytosol", field_p, field_r)

    extra_components = () if apparatus is None else ((AFM_INDENTER_COMPONENT, apparatus),)
    extra_connectors = () if apparatus is None else ((AFM_MEMBRANE_CONTACT, apparatus),)
    architecture = None if apparatus is None else afm_experiment_architecture()
    slice_ = build_interior_column_slice(
        membrane_owner=membrane, cortex_owner=cortex, cytosol_owner=cytosol, nucleus_owner=nucleus,
        cortex_transfer=cortex_transfer, membrane_boundary=membrane_boundary,
        nucleus_boundary=nucleus_boundary, cytosol_endpoint=endpoint,
        base_seed=5, device="cuda:0", clock=clock, architecture=architecture,
        additional_component_bindings=extra_components,
        additional_connector_bindings=extra_connectors)
    parts = dict(
        launch=launch, copy=copy, clock=clock, events=events, membrane=membrane, cortex=cortex,
        cytosol=cytosol, nucleus=nucleus, nuc_txn=nuc_txn, cortex_transfer=cortex_transfer,
        membrane_boundary=membrane_boundary, nucleus_boundary=nucleus_boundary, pc=pc, solver=solver,
        grid=grid, ruptured=ruptured)
    return slice_, parts


# ── 1. composition: four owners + three fluid connectors are all participants; none proposes ──────────
def test_slice_registers_four_owners_and_three_fluid_connectors() -> None:
    slice_, p = _slice()
    participants = slice_.transaction.world.participants
    expected = {id(p["membrane"]), id(p["cortex"]), id(p["cytosol"]), id(p["nucleus"]),
                id(p["cortex_transfer"]), id(p["membrane_boundary"]), id(p["nucleus_boundary"])}
    assert set(map(id, participants)) == expected
    assert len(participants) == 7
    # every interior coupling is a NON-kinetic field coupler → no proposer at all.
    assert slice_.transaction.proposers == ()


def test_experiment_apparatus_joins_the_same_transaction_once() -> None:
    apparatus = _TransactionSpy()
    slice_, _ = _slice(apparatus=apparatus)
    participants = slice_.transaction.world.participants
    assert sum(runtime is apparatus for runtime in participants) == 1
    accepted = _FakeAccepted(True)
    slice_.step(lambda: None, dt_phys=0.01, accepted_d=accepted)
    assert apparatus.calls == [
        ("snapshot",),
        ("rollback", accepted),
        ("commit", accepted, 0.01, 5),
    ]


def test_connectors_have_correct_names_and_endpoints() -> None:
    slice_, p = _slice()
    assert p["cortex_transfer"].name == "surface_porous_transfer"
    assert {p["cortex_transfer"].component_a, p["cortex_transfer"].component_b} == {"cortex", "cytosol"}
    assert p["membrane_boundary"].name == "membrane_cytosol_boundary"
    assert {p["membrane_boundary"].component_a, p["membrane_boundary"].component_b} == {"membrane", "cytosol"}
    assert p["nucleus_boundary"].name == "nucleus_cytosol_boundary"
    assert {p["nucleus_boundary"].component_a, p["nucleus_boundary"].component_b} == {"nucleus", "cytosol"}


def test_all_fluid_connectors_are_bidirectional_adjoint() -> None:
    slice_, p = _slice()
    # the two moving boundaries validate impermeable/adjoint in their own __post_init__; the porous transfer
    # advertises bidirectional + adjoint explicitly (Newton's 3rd law).
    assert p["cortex_transfer"].bidirectional and p["cortex_transfer"].adjoint_transfer_required
    assert p["membrane_boundary"].adjoint_transfer_required and not p["membrane_boundary"].impermeable
    assert p["nucleus_boundary"].adjoint_transfer_required and p["nucleus_boundary"].impermeable


# ── 2. REQUIRED-PARAM: the nucleus rupture strain (I0-B2 GAP) is supplied, never invented ─────────────
def test_core_transaction_rejects_unset_rupture_strain() -> None:
    for bad in (None, 0.0, -0.1, float("inf")):
        with pytest.raises(ValueError, match="eps_rupture"):
            CoreSurfaceTransaction(
                position_d=_fa(1, 6, wp.vec3d), position_snap_d=_fa(2, 6, wp.vec3d),
                faces_d=_fa(3, (3, 3), wp.int32), a0_d=_fa(4, 3, wp.float64), ruptured_d=_fa(5, 3, wp.int32),
                ruptured_snap_d=_fa(6, 3, wp.int32), eps_rupture=bad)


# ── 3. cortex ⟷ cytosol transfer scatters the pressure traction onto the CORTEX force array ───────────
def test_porous_transfer_scatters_pressure_traction_onto_cortex() -> None:
    slice_, p = _slice()
    slice_.accumulate()
    # the porous transfer called the injected PressureCoupling once, onto the CORTEX force array.
    assert len(p["pc"].calls) == 1
    state, force = p["pc"].calls[0]
    assert force is p["cortex"].force_d          # scatter target is the cortex-owned force (Newton's 3rd law)
    assert state.node_pos is p["cortex"].position_d
    assert state.node_volume is p["cortex_transfer"]._node_volume_d
    # the two moving boundaries also scattered their adjoint traction onto their own surfaces.
    m = p["membrane_boundary"].delegate.calls
    n = p["nucleus_boundary"].delegate.calls
    assert ("membrane_b", "update_geometry") in m and ("membrane_b", "accumulate") in m
    assert ("nucleus_b", "update_geometry") in n and ("nucleus_b", "accumulate") in n
    # geometry updates precede the adjoint tractions (moved boundary → field → traction).
    assert m.index(("membrane_b", "update_geometry")) < m.index(("membrane_b", "accumulate"))


# ── 4. the cytosol field participant snapshots + reject-restores via the REAL scheduler kernels ───────
def test_cytosol_field_snapshot_and_reject_restore_use_real_scheduler_kernels() -> None:
    slice_, p = _slice()
    txn = p["cytosol"].transaction

    txn.snapshot_candidate()
    # snapshot D2D-copies snap ← live for the four float64 fields (p, div_vs, s_membrane, s_total) + mask.
    assert len(p["copy"].calls) >= 5
    copied_srcs = [src for _dst, src in p["copy"].calls]
    assert p["grid"].p in copied_srcs and p["grid"].mask in copied_srcs

    p["launch"].calls.clear()
    accepted = _FakeAccepted(accept=False)
    txn.rollback(accepted)
    # four float64 restores + one int32 restore, all reject-gated on the SAME predicate.
    assert txn.launched_kernels == [
        _conditional_restore_f64_3d_kernel, _conditional_restore_f64_3d_kernel,
        _conditional_restore_f64_3d_kernel, _conditional_restore_f64_3d_kernel,
        _conditional_restore_i32_3d_kernel,
    ]
    for _kernel, _dim, inputs, _out, _dev in p["launch"].calls:
        assert inputs[2] is accepted  # scheduler restore kernel signature is (state, start, accepted)


# ── 5. the nucleus commit runs the REAL accepted-gated envelope-rupture kernel on the converged surface ─
def test_nucleus_commit_runs_accepted_gated_rupture_kernel() -> None:
    slice_, p = _slice()
    accepted = _FakeAccepted(accept=True)

    p["nuc_txn"].commit_irreversible(accepted, dt_phys=0.01, rng_seed=9)

    assert p["nuc_txn"].launched_kernels == [conditional_rupture_update_kernel]
    inputs = p["launch"].calls[-1][2]
    assert inputs[0] is p["nuc_txn"]._position_d      # converged surface geometry
    assert inputs[1] is p["nuc_txn"]._faces_d
    assert inputs[4] is accepted                       # accepted-gated (rejected ⇒ no rupture)
    assert inputs[5] is p["ruptured"]                  # the per-face tear flags


def test_nucleus_rejected_step_bit_restores_surface_and_ruptured() -> None:
    slice_, p = _slice()
    p["launch"].calls.clear()
    accepted = _FakeAccepted(accept=False)

    p["nuc_txn"].rollback(accepted)

    assert p["nuc_txn"].launched_kernels == [
        _restore_vec3_if_rejected_kernel,     # surface position
        _restore_int32_if_rejected_kernel,    # per-face ruptured flag
    ]
    for _kernel, _dim, inputs, _out, _dev in p["launch"].calls:
        assert inputs[0] is accepted  # every restore is device-gated on the predicate


# ── 6. the clock advances only on an accepted step; every step forwards the SAME dt ───────────────────
def test_clock_advances_only_on_accept() -> None:
    slice_, p = _slice()
    clock = p["clock"]

    slice_.step(lambda: None, dt_phys=0.02, accepted_d=_FakeAccepted(accept=False))  # rejected
    assert clock.epoch == 0
    slice_.step(lambda: None, dt_phys=0.02, accepted_d=_FakeAccepted(accept=True))   # accepted
    assert clock.epoch == 1
    slice_.step(lambda: None, dt_phys=0.02, accepted_d=_FakeAccepted(accept=False))  # rejected again
    assert clock.epoch == 1
    assert [dt for _pred, dt in clock.advances] == [0.02, 0.02, 0.02]


# ── 7. whole-step order: snapshot → solve → commit(rupture) → advance, one predicate to all ───────────
def test_step_runs_snapshot_solve_commit_advance_in_order() -> None:
    slice_, p = _slice()
    events = p["events"]

    def solve() -> None:
        events.append(("solve",))

    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=True))

    kinds = [e[0] for e in events]
    i_snapshot = kinds.index("copy")  # a participant's snapshot opens the candidate.
    i_solve = kinds.index("solve")
    i_rupture = next(n for n, e in enumerate(events)
                     if e[0] == "launch" and e[1] is conditional_rupture_update_kernel)
    i_advance = kinds.index("advance")
    assert i_snapshot < i_solve < i_rupture < i_advance


# ── 8. every moving boundary participates (snapshot/rollback/commit) each accepted step under one predicate ─
def test_moving_boundaries_participate_in_the_transaction() -> None:
    slice_, p = _slice()

    slice_.step(lambda: None, dt_phys=0.02, accepted_d=_FakeAccepted(accept=True))

    for boundary in (p["membrane_boundary"], p["nucleus_boundary"]):
        kinds = [c[1] for c in boundary.delegate.calls]
        assert "snapshot" in kinds and "rollback" in kinds and "commit" in kinds
        assert kinds.index("snapshot") < kinds.index("rollback") < kinds.index("commit")
