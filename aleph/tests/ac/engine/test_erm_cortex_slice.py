"""CPU-green structural + numeric gates for the GATE-B dynamic membrane--cortex ERM slice.

Proves — with recording CUDA doubles + the real NumPy tether oracle, no device — that
:class:`~aleph.engine.erm_cortex_slice.ERMCortexSlice` drives the accepted-step transaction so the
membrane--cortex ERM tether tension EMERGES from the bound-tether population and the Bell-slip shedding is the
mechanistic bleb-onset:

* the two-array Bell force is Newton's 3rd law across the two never-merged owner arrays (membrane ``-f`` /
  cortex ``+f`` sum to zero to ~1e-12), and is force-free at formation (``rest := length``) and under
  compression (a molecular linker is a tether, not a strut);
* the ERM Bell off-rate RISES monotonically with tensile load (pure SLIP bond), so higher tension ⇒ higher
  detach probability ⇒ the bound population sheds (the CUDA-lane stochastic drop is the Lead's native gate);
* a rejected step bit-restores the ``bound``/``rest``/``epoch`` SoA (reject-gated int32 + float64 kernels) and
  advances no clock; the clock advances only on accept;
* the REQUIRED-PARAM guards fire — the four ERM Bell GAP constants are supplied, never invented.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest
import warp as wp

from aleph.components.incumbent.erm_tether import (
    erm_bell_force_pair_kernel,
    erm_bell_force_pair_reference,
    erm_bell_kmc_pair_kernel,
    increment_erm_epoch_if_accepted_kernel,
)
from aleph.engine.erm_cortex_connector import (
    ErmCortexConnector,
    _restore_float64_if_rejected_kernel,
    _restore_int32_if_rejected_kernel,
)
from aleph.engine.erm_cortex_slice import (
    ERMCortexSlice,
    ERMSliceParams,
    build_erm_cortex_slice,
)
from aleph.engine.surface_body import SurfaceComponentStateOwner
from aleph.components.motor.bell_kinetics_analytic import (
    attach_probability,
    bell_off_rate,
    detach_probability,
)


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


def _fa(ptr: int, n: int, dtype: object) -> _FakeArray:
    return _FakeArray(ptr, (n,), dtype)


@dataclass(slots=True)
class _FakeAccepted:
    """Host-readable acceptance predicate double; the device kernel gates in-device, the spy clock in Python."""

    accept: bool


# ── Recording launch / copy + surface-owner + clock spies ─────────────────────────────────────────────
@dataclass(slots=True)
class _RecordingLauncher:
    events: list = field(default_factory=list)
    calls: list = field(default_factory=list)

    def __call__(self, kernel, *, dim, inputs, outputs=None, device=None) -> None:
        self.calls.append((kernel, int(dim), tuple(inputs), tuple(outputs or ()), device))
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


@dataclass(slots=True)
class _GatingSpyClock:
    """Event-clock double that GATES the epoch on a host-readable predicate (mirrors the in-kernel gate)."""

    events: list = field(default_factory=list)
    base_seed: int = 7
    epoch: int = 0
    advances: list = field(default_factory=list)

    def advance(self, accepted, dt_phys) -> None:
        self.advances.append((accepted, dt_phys))
        self.events.append(("advance", accepted, dt_phys))
        if getattr(accepted, "accept", False):
            self.epoch += 1


# ── Builders: a real slice assembled from doubles ─────────────────────────────────────────────────────
def _surface_owner(name: str, pos: _FakeArray, force: _FakeArray) -> SurfaceComponentStateOwner:
    return SurfaceComponentStateOwner(
        name=name, position_d=pos, force_d=force,
        mechanics=_MechanicsSpy(), transaction=_TransactionSpy(), ledger=_LedgerSpy(),
    )


def _connector(
    launch: _RecordingLauncher, copy: _RecordingCopy, mem_pos: _FakeArray, ctx_pos: _FakeArray,
    *, k_on: float = 1.0, k_off0: float = 0.3, bell_force_pn: float = 5.0, capture_radius_um: float = 0.05,
) -> ErmCortexConnector:
    return ErmCortexConnector(
        membrane_idx_d=_fa(800, 4, wp.int32), cortex_idx_d=_fa(801, 4, wp.int32),
        bound_d=_fa(802, 4, wp.int32), rest_d=_fa(803, 4, wp.float64), rng_epoch_d=_fa(804, 1, wp.int32),
        bound_snap_d=_fa(812, 4, wp.int32), rest_snap_d=_fa(813, 4, wp.float64),
        rng_epoch_snap_d=_fa(814, 1, wp.int32),
        detach_events_d=_fa(820, 1, wp.int32), attach_events_d=_fa(821, 1, wp.int32),
        membrane_pos_d=mem_pos, cortex_pos_d=ctx_pos,
        k_erm=4600.0, k_on=k_on, k_off0=k_off0, bell_force_pn=bell_force_pn,
        capture_radius_um=capture_radius_um, launch=launch, copy=copy,
    )


def _slice():
    events: list = []
    launch = _RecordingLauncher(events=events)
    copy = _RecordingCopy(events=events)
    clock = _GatingSpyClock(events=events)
    mem_pos, mem_f = _fa(10, 4, wp.vec3d), _fa(11, 4, wp.vec3d)
    ctx_pos, ctx_f = _fa(20, 4, wp.vec3d), _fa(21, 4, wp.vec3d)
    membrane = _surface_owner("membrane", mem_pos, mem_f)
    cortex = _surface_owner("cortex", ctx_pos, ctx_f)
    connector = _connector(launch, copy, mem_pos, ctx_pos)
    slice_ = build_erm_cortex_slice(
        membrane_owner=membrane, cortex_owner=cortex, connector=connector,
        base_seed=7, device="cuda:0", clock=clock,
    )
    return slice_, launch, copy, clock, events, membrane, cortex, connector


# ── 1. REQUIRED-PARAM discipline: the 4 ERM Bell GAPs are supplied, never invented ─────────────────────
def test_params_reject_unset_gap_rates_and_source() -> None:
    for missing in ("k_on", "k_off0", "bell_force_pn", "capture_radius_um"):
        kw = dict(k_erm=4600.0, k_on=1.0, k_off0=0.3, bell_force_pn=5.0, capture_radius_um=0.05, source="x")
        kw[missing] = None
        with pytest.raises(ValueError, match=missing):
            ERMSliceParams(**kw)
    with pytest.raises(ValueError, match="source"):
        ERMSliceParams(k_erm=4600.0, k_on=1.0, k_off0=0.3, bell_force_pn=5.0, capture_radius_um=0.05, source="  ")
    # a fully-supplied set builds and yields the source-gated ERMBellKinetics contract the KMC consumes.
    p = ERMSliceParams(k_erm=4600.0, k_on=1.0, k_off0=0.3, bell_force_pn=5.0, capture_radius_um=0.05,
                       source="PROVISIONAL PI-GAP (test)")
    bk = p.to_bell_kinetics()
    assert bk.k_on_s == 1.0 and bk.k_off0_s == 0.3 and bk.capture_radius_um == 0.05
    assert bk.rebind_rest_policy == "formation_length"


def test_connector_rejects_nonpositive_or_infinite_kinetics() -> None:
    launch, copy = _RecordingLauncher(), _RecordingCopy()
    mem_pos, ctx_pos = _fa(10, 4, wp.vec3d), _fa(20, 4, wp.vec3d)
    with pytest.raises(ValueError, match="k_on"):
        _connector(launch, copy, mem_pos, ctx_pos, k_on=0.0)
    with pytest.raises(ValueError, match="k_off0"):
        _connector(launch, copy, mem_pos, ctx_pos, k_off0=float("inf"))
    with pytest.raises(ValueError, match="capture_radius_um"):
        _connector(launch, copy, mem_pos, ctx_pos, capture_radius_um=-0.05)


# ── 2. two-array force is Newton's 3rd law (membrane -f + cortex +f = 0) across the split arrays ────────
def test_two_array_force_is_newton_third_law_and_tensile() -> None:
    # two stretched tethers (length 0.2 µm > rest 0.1 µm) between distinct membrane/cortex nodes.
    pos_m = np.array([[0.0, 0.0, 1.2], [0.10, 0.0, 1.25]])
    pos_c = np.array([[0.0, 0.0, 1.0], [0.10, 0.0, 1.05]])
    out = erm_bell_force_pair_reference(
        pos_m, pos_c, membrane_idx=np.array([0, 1]), cortex_idx=np.array([0, 1]),
        bound=np.array([1, 1]), rest=np.array([0.1, 0.1]), k_erm=4600.0)
    # Newton's 3rd law across the two never-merged owner arrays: total force sums to zero.
    total = out["force_m"].sum(axis=0) + out["force_c"].sum(axis=0)
    assert np.linalg.norm(total) < 1e-12
    # tensile: each tether pulls the membrane inward (−z) and the cortex outward (+z).
    assert out["tension"][0] > 0.0 and out["tension"][1] > 0.0
    assert out["force_m"][0][2] < 0.0 and out["force_c"][0][2] > 0.0
    # the per-tether load is exactly k_erm·(length − rest) = 4600·(0.2 − 0.1).
    assert abs(out["tension"][0] - 4600.0 * 0.1) < 1e-9


def test_force_is_zero_at_formation_and_under_compression() -> None:
    pos_c = np.array([[0.0, 0.0, 1.0]])
    pos_m_form = np.array([[0.0, 0.0, 1.1]])
    # formation: rest := current length bit-identically (exactly what the KMC's ``rest[k] = length`` does), so
    # tension is EXACTLY 0 ⇒ force 0 (binding injects no unrecorded prestress — no float-rounding prestress).
    rest = np.linalg.norm(pos_m_form - pos_c, axis=1)
    at_formation = erm_bell_force_pair_reference(
        pos_m_form, pos_c, np.array([0]), np.array([0]), np.array([1]), rest, 4600.0)
    assert at_formation["tension"][0] == 0.0
    assert np.linalg.norm(at_formation["force_m"]) == 0.0
    assert np.linalg.norm(at_formation["force_c"]) == 0.0
    # compression (length < rest) is also force-free — a molecular linker is a tether, not a strut.
    compressed = erm_bell_force_pair_reference(
        np.array([[0.0, 0.0, 1.05]]), pos_c, np.array([0]), np.array([0]), np.array([1]), rest, 4600.0)
    assert compressed["tension"][0] == 0.0
    assert np.linalg.norm(compressed["force_m"]) == 0.0
    # a broken tether (bound == 0) contributes nothing.
    broken = erm_bell_force_pair_reference(
        np.array([[0.0, 0.0, 1.5]]), pos_c, np.array([0]), np.array([0]), np.array([0]), rest, 4600.0)
    assert np.linalg.norm(broken["force_m"]) == 0.0


# ── 3. the ERM KMC sheds under tensile load (pure SLIP: off-rate RISES with |F|) ───────────────────────
def test_kmc_bell_off_rate_rises_with_tensile_load() -> None:
    k_off0, f0, tau, k_on = 0.3, 5.0, 0.01, 1.0
    forces = np.array([0.0, 1.0, 3.0, 6.0])
    p_off = bell_off_rate(forces, k_off0, f0)
    # ERM is a pure SLIP bond: p_off(0) = k_off0, strictly increasing in load |F|.
    assert p_off[0] == k_off0
    assert np.all(np.diff(p_off) > 0.0)
    # higher tension ⇒ higher detach probability ⇒ the bound population sheds (native gate: the stochastic drop).
    p_det = detach_probability(tau, p_off)
    assert np.all(np.diff(p_det) > 0.0)
    # rebind within capture is a Poisson attach in (0, 1) — the geometric-within-capture KMC rebind.
    assert 0.0 < attach_probability(tau, k_on) < 1.0


# ── 4. composition: membrane + cortex are state-owners, ERM is the only connector, only ERM proposes ───
def test_slice_registers_two_surface_owners_and_erm_connector() -> None:
    slice_, _launch, _copy, _clock, _events, membrane, cortex, connector = _slice()
    participants = slice_.transaction.world.participants
    assert set(map(id, participants)) == {id(membrane), id(cortex), id(connector)}
    # only the ERM connector proposes events; the surface owners have no propose_events.
    assert slice_.transaction.proposers == (connector,)
    assert connector.name == "membrane_erm_cortex"
    assert connector.component_a == "membrane" and connector.component_b == "cortex"


def test_propose_events_is_a_pure_no_op() -> None:
    slice_, launch, copy, _clock, _events, *_ = _slice()
    launch.calls.clear()
    copy.calls.clear()
    # ERM pairs are fixed; rebind lives inside the KMC — propose writes no force, state, or candidate buffer.
    result = slice_.connector.propose_events(rates=None, dt_phys=0.01, rng_seed=7, neighbors=None)
    assert result is None
    assert launch.calls == [] and copy.calls == []


# ── 5. accumulate scatters the two-array Bell tether (membrane = a, cortex = b), never merged ──────────
def test_accumulate_launches_two_array_force_membrane_a_cortex_b() -> None:
    slice_, launch, _copy, _clock, _events, membrane, cortex, connector = _slice()

    slice_.accumulate()

    assert launch.kernels == [erm_bell_force_pair_kernel]
    inputs = launch.calls[0][2]
    assert inputs[0] is membrane.position_d   # pos_m  (membrane = component_a)
    assert inputs[1] is membrane.force_d       # force_m ← −f
    assert inputs[2] is cortex.position_d      # pos_c  (cortex = component_b)
    assert inputs[3] is cortex.force_d         # force_c ← +f
    assert inputs[4] is connector.membrane_idx_d
    assert inputs[5] is connector.cortex_idx_d
    assert inputs[6] is connector.bound_d
    assert inputs[7] is connector.rest_d
    assert membrane.force_d is not cortex.force_d  # two force arrays never merged


# ── 6. commit runs the accepted-gated two-array KMC + epoch on the converged geometry ─────────────────
def test_commit_runs_accepted_gated_kmc_and_epoch() -> None:
    slice_, launch, _copy, _clock, _events, _membrane, _cortex, connector = _slice()
    accepted = _FakeAccepted(accept=True)

    connector.commit_irreversible(accepted, dt_phys=0.01, rng_seed=131)

    assert launch.kernels == [erm_bell_kmc_pair_kernel, increment_erm_epoch_if_accepted_kernel]
    kmc_inputs = launch.calls[0][2]
    assert kmc_inputs[0] is connector._membrane_pos_d   # live membrane pos (converged geometry)
    assert kmc_inputs[1] is connector._cortex_pos_d     # live cortex pos
    assert kmc_inputs[2] is connector.membrane_idx_d
    assert kmc_inputs[3] is connector.cortex_idx_d
    assert kmc_inputs[4] is connector.bound_d           # the flipped bound flags
    assert kmc_inputs[6] is connector.rest_d            # rebind sets rest := length (force-free formation)
    assert kmc_inputs[-4] is connector.rng_epoch_d      # RNG folds the epoch
    assert kmc_inputs[-3] is accepted                   # accepted-gated (rejected ⇒ no-op)
    assert kmc_inputs[-2] is connector.detach_events_d
    assert kmc_inputs[-1] is connector.attach_events_d
    # the epoch advance is device-gated on the SAME accepted predicate.
    epoch_inputs = launch.calls[1][2]
    assert epoch_inputs[0] is accepted and epoch_inputs[1] is connector.rng_epoch_d


# ── 7. snapshot copies the SoA; a rejected step bit-restores bound + rest + epoch (reject-gated) ───────
def test_snapshot_and_rejected_step_bit_restore_bound_rest_epoch() -> None:
    slice_, launch, copy, _clock, _events, _membrane, _cortex, connector = _slice()

    connector.snapshot_candidate()
    # snapshot D2D-copies snap ← live for bound (int32), rest (float64), rng_epoch (int32).
    assert [(dst, src) for dst, src in copy.calls] == [
        (connector.bound_snap_d, connector.bound_d),
        (connector.rest_snap_d, connector.rest_d),
        (connector.rng_epoch_snap_d, connector.rng_epoch_d),
    ]

    launch.calls.clear()
    accepted = _FakeAccepted(accept=False)
    connector.rollback(accepted)
    # rollback launches the reject-gated restore kernels: int32 bound, float64 rest, int32 epoch.
    assert launch.kernels == [
        _restore_int32_if_rejected_kernel,      # bound
        _restore_float64_if_rejected_kernel,    # rest
        _restore_int32_if_rejected_kernel,      # rng_epoch
    ]
    for _kernel, _dim, inputs, _outputs, _device in launch.calls:
        assert inputs[0] is accepted  # every restore is device-gated on the predicate
    assert launch.calls[0][2][1] is connector.bound_d and launch.calls[0][2][2] is connector.bound_snap_d
    assert launch.calls[1][2][1] is connector.rest_d and launch.calls[1][2][2] is connector.rest_snap_d
    assert launch.calls[2][2][1] is connector.rng_epoch_d


# ── 8. the clock advances only on an accepted step; every step forwards the SAME predicate + dt ────────
def test_clock_advances_only_on_accept() -> None:
    slice_, _launch, _copy, clock, _events, *_ = _slice()

    def solve() -> None:
        return None

    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=False))  # rejected
    assert clock.epoch == 0
    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=True))   # accepted
    assert clock.epoch == 1
    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=False))  # rejected again
    assert clock.epoch == 1
    assert [dt for _pred, dt in clock.advances] == [0.02, 0.02, 0.02]


# ── 9. whole-step order: snapshot → (propose no-op) → solve → commit(KMC) → advance ────────────────────
def test_step_runs_snapshot_solve_commit_advance_in_order() -> None:
    slice_, _launch, _copy, _clock, events, _membrane, _cortex, _connector = _slice()

    def solve() -> None:
        events.append(("solve",))

    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=True))

    kinds = [e[0] for e in events]
    i_snapshot = kinds.index("copy")  # the connector's SoA snapshot (bound/rest/epoch) opens the candidate.
    i_solve = kinds.index("solve")
    i_kmc = next(n for n, e in enumerate(events)
                 if e[0] == "launch" and e[1] is erm_bell_kmc_pair_kernel)
    i_advance = kinds.index("advance")
    # snapshot before the inner solve; the KMC commit after solve (converged geometry); clock last.
    assert i_snapshot < i_solve < i_kmc < i_advance


# ── 10. the surface state-owners participate (snapshot/rollback/commit each accepted step) ─────────────
def test_surface_owners_participate_in_the_transaction() -> None:
    slice_, _launch, _copy, _clock, _events, membrane, cortex, _connector = _slice()

    slice_.step(lambda: None, dt_phys=0.02, accepted_d=_FakeAccepted(accept=True))

    for owner in (membrane, cortex):
        kinds = [c[0] for c in owner.transaction.calls]
        assert kinds == ["snapshot", "rollback", "commit"]  # forwarded once each per step, in order
