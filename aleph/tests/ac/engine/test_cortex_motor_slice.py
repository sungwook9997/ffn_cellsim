"""CPU-green structural + numeric gates for the GATE-B dynamic cortex-motor vertical slice.

Proves — with recording CUDA doubles + the real NumPy crossbridge oracle, no device — that
:class:`~aleph.engine.cortex_motor_slice.CortexMotorSlice` drives the accepted-step transaction so
cortical NMII tension EMERGES from binding events:

* the transaction PROPOSES attach events (the query accelerator fills the connector's attach scratch and the
  live segment topology is re-bound each step);
* the KMC ATTACH is accepted-predicated (the attach kernel flips ``bound`` under the forced-accept predicate,
  reading ``k_on`` from the hand params and the query scratch);
* the split crossbridge scatters ``+f`` to the head array and ``−(1−t)f``/``−t·f`` to the two cortex port
  nodes (Newton's 3rd law, momentum-conserving to ~1e-12 — checked against the landed host oracle);
* a rejected step bit-restores the binding SoA and advances no clock; the clock advances only on accept.

The CUDA-lane numeric behaviour (bound-fraction growth, γ rise) is the Lead's native gate on the A5000.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest
import warp as wp

from aleph.engine.cortex_motor_slice import (
    CortexMotorConnector,
    CortexMotorParams,
    CortexMotorSlice,
    _ActuatorStateLedger,
    _ActuatorStateTransaction,
    _restore_int32_if_rejected_kernel,
    _restore_int64_if_rejected_kernel,
    _restore_vec3_if_rejected_kernel,
)
from aleph.engine.nmii_actuator import (
    CORTEX_COMPONENT,
    NMII_COMPONENT,
    NMII_CORTEX_MOTOR,
    BackboneArmMechanics,
    FilamentMotorPortView,
    NMIIActuatorStateOwner,
    NMIIBendingStiffness,
    SegmentConnectorState,
    proposed_nmii_architecture,
)
from aleph.engine.protrusion import proposed_protrusion_architecture
from aleph.engine.contracts import reference_cell_architecture
from aleph.engine.actor import CellActor
from aleph.engine.transaction import CellTransaction
from aleph.engine.world import CellWorldTransaction
from aleph.components.motor.backbone_warp import angle_harmonic_kernel
from aleph.components.motor.minifilament_warp import harmonic_bond_kernel
from aleph.components.motor.segment_motor import (
    SegmentDetachKinetics,
    _increment_epoch_if_accepted_kernel,
    _refresh_bound_walk_dir_kernel,
    attach_segment_gated_r0bind_kernel,
    compute_head_loads_segment_split_r0bind_kernel,
    crossbridge_segment_split_r0bind_kernel,
    crossbridge_segment_split_r0bind_reference,
    r0_bind_at_attach_reference,
    refresh_segment_barbed_kernel,
    step_detach_segment_gated_catch_slip_kernel,
    step_detach_segment_gated_kernel,
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


# ── Recording launch / copy + spies ──────────────────────────────────────────────────────────────────
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
class _SpyQuery:
    """Attach-query accelerator double: records the split-ownership view + port it was asked to query."""

    events: list = field(default_factory=list)
    calls: list = field(default_factory=list)

    def fill_attach_scratch(self, actuator, port) -> None:
        self.calls.append((actuator, port))
        self.events.append(("propose",))


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
def _architecture():
    return proposed_nmii_architecture(proposed_protrusion_architecture(reference_cell_architecture()))


def _bending() -> NMIIBendingStiffness:
    # DERIVED (KB/PI GAP-flagged), never fit to a gate.
    return NMIIBendingStiffness(persistence_length_um=1.0, segment_len_um=0.1, k_xb_pn_per_um=1000.0,
                                r0_head_um=0.02)


def _actuator_state(launch: _RecordingLauncher, copy: _RecordingCopy) -> NMIIActuatorStateOwner:
    position_d = _fa(100, 10, wp.vec3d)
    active_d = _fa(103, 2, wp.int32)
    atp_state_d = _fa(111, 4, wp.int32)
    atp_consumed_d = _fa(112, 1, wp.int64)
    epoch_d = _fa(113, 1, wp.int32)
    mechanics = BackboneArmMechanics(
        n_particles=10, n_heads=4,
        backbone_bonds_d=_fa(120, 2, wp.int32), head_bonds_d=_fa(121, 4, wp.int32),
        backbone_angles_d=_fa(122, 1, wp.int32), head_arm_angles_d=_fa(123, 4, wp.int32),
        k_backbone_pn_per_um=1.0e4, r0_backbone_um=0.05, k_head_spring_pn_per_um=1.0e3, r0_head_um=0.02,
        bending=_bending(), launch=launch,
    )
    transaction = _ActuatorStateTransaction(
        position_d=position_d, active_minifilament_d=active_d, atp_cycle_state_d=atp_state_d,
        atp_consumed_d=atp_consumed_d, actuator_epoch_d=epoch_d,
        position_snap_d=_fa(200, 10, wp.vec3d), active_snap_d=_fa(203, 2, wp.int32),
        atp_state_snap_d=_fa(211, 4, wp.int32), atp_consumed_snap_d=_fa(212, 1, wp.int64),
        epoch_snap_d=_fa(213, 1, wp.int32), launch=launch, copy=copy,
    )
    return NMIIActuatorStateOwner(
        n_minifilaments=2, n_heads=4, position_d=position_d, force_d=_fa(101, 10, wp.vec3d),
        minifilament_offset_d=_fa(102, 3, wp.int32), active_minifilament_d=active_d,
        minifilament_id_d=_fa(104, 2, wp.int64), particle_role_d=_fa(105, 10, wp.int32),
        head_node_d=_fa(106, 4, wp.int32), head_id_d=_fa(107, 4, wp.int64), head_side_d=_fa(108, 4, wp.int32),
        head_load_hill_d=_fa(109, 4, wp.float64), head_load_bell_d=_fa(110, 4, wp.float64),
        atp_cycle_state_d=atp_state_d, atp_consumed_d=atp_consumed_d, actuator_epoch_d=epoch_d,
        mechanics=mechanics, transaction=transaction, ledger=_ActuatorStateLedger(),
    )


def _connector_state() -> SegmentConnectorState:
    return SegmentConnectorState(
        n_heads=4, n_segments=3,
        bound_d=_fa(800, 4, wp.int32), seg_id_d=_fa(801, 4, wp.int32), seg_a_d=_fa(802, 4, wp.int32),
        seg_b_d=_fa(803, 4, wp.int32), bary_t_d=_fa(804, 4, wp.float64), abscissa_d=_fa(805, 4, wp.float64),
        walk_dir_d=_fa(806, 4, wp.vec3d), rng_epoch_d=_fa(807, 1, wp.int32),
        loads_hill_d=_fa(808, 4, wp.float64), loads_bell_d=_fa(809, 4, wp.float64),
        seg_barbed_d=_fa(810, 3, wp.vec3d), query_seg_id_d=_fa(811, 4, wp.int32),
        query_t_d=_fa(812, 4, wp.float64), query_barbed_d=_fa(813, 4, wp.vec3d),
        bound_snap_d=_fa(820, 4, wp.int32), seg_id_snap_d=_fa(821, 4, wp.int32),
        seg_a_snap_d=_fa(822, 4, wp.int32), seg_b_snap_d=_fa(823, 4, wp.int32),
        bary_t_snap_d=_fa(824, 4, wp.float64), abscissa_snap_d=_fa(825, 4, wp.float64),
        walk_dir_snap_d=_fa(826, 4, wp.vec3d), rng_epoch_snap_d=_fa(827, 1, wp.int32),
    )


def _cortex_port() -> FilamentMotorPortView:
    return FilamentMotorPortView(
        component=CORTEX_COMPONENT, position_d=_fa(300, 8, wp.vec3d), force_d=_fa(301, 8, wp.vec3d),
        segment_node_a_d=_fa(302, 3, wp.int32), segment_node_b_d=_fa(303, 3, wp.int32),
        segment_polarity_d=_fa(304, 3, wp.int32), persistent_filament_id_d=_fa(305, 3, wp.int64),
        material_s0_d=_fa(306, 3, wp.float64), material_s1_d=_fa(307, 3, wp.float64),
        topology_epoch_d=_fa(308, 1, wp.int32),
    )


class _HandParamsDouble:
    """k_xb/r0_xb + a marker so the KMC-kernel inputs carry the k_on-bearing hand params through the seam."""

    k_xb = wp.float64(1000.0)
    r0_xb = wp.float64(0.0)
    k_on = wp.float64(50.0)


def _connector(
    launch: _RecordingLauncher, copy: _RecordingCopy, query: _SpyQuery,
    actuator: NMIIActuatorStateOwner, port: FilamentMotorPortView,
    detach: SegmentDetachKinetics = SegmentDetachKinetics.CATCH_SLIP,
) -> CortexMotorConnector:
    return CortexMotorConnector(
        state=_connector_state(), head_node_d=actuator.head_node_d, params=_HandParamsDouble(),
        catch_slip=(object() if detach is SegmentDetachKinetics.CATCH_SLIP else None),
        query=query, actuator_view=actuator.geometry, port=port,
        r0_bind_d=_fa(850, 4, wp.float64), r0_bind_snap_d=_fa(851, 4, wp.float64),
        detach_kinetics=detach, launch=launch, copy=copy,
    )


def _slice(detach: SegmentDetachKinetics = SegmentDetachKinetics.CATCH_SLIP):
    events: list = []
    launch = _RecordingLauncher(events=events)
    copy = _RecordingCopy(events=events)
    query = _SpyQuery(events=events)
    clock = _GatingSpyClock(events=events)
    actuator = _actuator_state(launch, copy)
    port = _cortex_port()
    connector = _connector(launch, copy, query, actuator, port, detach)
    actor = CellActor(_architecture())
    actor.bind_component(NMII_COMPONENT, actuator)
    actor.bind_connector(NMII_CORTEX_MOTOR, connector)
    world = CellWorldTransaction(actor, require_complete=False)
    txn = CellTransaction(world, clock)
    slice_ = CortexMotorSlice(transaction=txn, actuator_state=actuator, connector=connector, port=port,
                              clock=clock)
    return slice_, launch, copy, query, clock, events


# ── 1. REQUIRED-PARAM discipline: GAP inputs are supplied, never invented ──────────────────────────────
def _sourced_kwargs() -> dict:
    return dict(v0=0.5, f_stall=2.0, kappa=0.5, k_xb=1000.0, r0_head=0.02, r0_xb=0.0,
               capture_radius=0.05, k_off0=0.35, f0=7.1)


def test_params_reject_unset_k_on_and_catch_slip_gaps() -> None:
    with pytest.raises(ValueError, match="k_on"):
        CortexMotorParams(k_on=None, **_sourced_kwargs())
    # CATCH_SLIP (default) needs the two-pathway constants — they are a GAP, never defaulted.
    with pytest.raises(ValueError, match="k_catch0"):
        CortexMotorParams(k_on=50.0, **_sourced_kwargs())
    # SLIP smoke path builds on the provisional-sourced Bell inputs without the catch-slip GAP.
    slip = CortexMotorParams(k_on=50.0, detach_kinetics=SegmentDetachKinetics.SLIP, **_sourced_kwargs())
    assert slip.detach_kinetics is SegmentDetachKinetics.SLIP
    full = CortexMotorParams(k_on=50.0, k_catch0=0.2, x_catch=5.0e-4, k_slip0=0.35, x_slip=6.0e-4,
                             **_sourced_kwargs())
    assert full.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP


def test_catch_slip_connector_requires_catch_slip_params() -> None:
    events: list = []
    launch, copy, query = _RecordingLauncher(events=events), _RecordingCopy(events=events), _SpyQuery(events=events)
    actuator = _actuator_state(launch, copy)
    with pytest.raises(ValueError, match="REQUIRED-PARAM.*CATCH_SLIP"):
        CortexMotorConnector(
            state=_connector_state(), head_node_d=actuator.head_node_d, params=_HandParamsDouble(),
            catch_slip=None, query=query, actuator_view=actuator.geometry, port=_cortex_port(),
            r0_bind_d=_fa(850, 4, wp.float64), r0_bind_snap_d=_fa(851, 4, wp.float64),
            detach_kinetics=SegmentDetachKinetics.CATCH_SLIP, launch=launch, copy=copy,
        )


# ── 2. Composition: cortex is a PORT; the two participants are nmii + the connector ────────────────────
def test_slice_registers_two_participants_and_one_proposer() -> None:
    slice_, *_ = _slice()
    participants = slice_.transaction.world.participants
    assert set(map(id, participants)) == {id(slice_.actuator_state), id(slice_.connector)}
    # only the connector proposes events; the NMII state-owner has no propose_events (cortex is not bound).
    assert slice_.transaction.proposers == (slice_.connector,)
    assert slice_.port.component == CORTEX_COMPONENT


# ── 3. propose_events fills the attach scratch + re-binds the LIVE segment topology (items 2 & 3) ──────
def test_propose_events_queries_live_port_and_rebinds_segments() -> None:
    slice_, _launch, _copy, query, _clock, _events = _slice()
    connector, port = slice_.connector, slice_.port

    connector.propose_events(rates=None, dt_phys=0.01, rng_seed=7, neighbors=None)

    # item 2: the attach-query accelerator was run against the split-ownership actuator view + cortex port.
    assert len(query.calls) == 1
    view, queried_port = query.calls[0]
    assert view.position_d is slice_.actuator_state.position_d  # nmii-owned head array
    assert queried_port is port
    # item 3: the LIVE cortex segment endpoints are the ones the accepted attach kernel will read.
    assert connector._target_seg_node_a is port.segment_node_a_d
    assert connector._target_seg_node_b is port.segment_node_b_d


# ── 4. accumulate scatters the split crossbridge + the minifilament's own internal mechanics ──────────
def test_accumulate_wires_internal_mechanics_then_split_crossbridge() -> None:
    slice_, launch, _copy, _query, _clock, _events = _slice()

    slice_.accumulate()

    # four internal backbone/head-arm kernels, then the connector's live-geometry + attach-unstrained crossbridge.
    assert launch.kernels == [
        harmonic_bond_kernel,          # backbone axial rod
        harmonic_bond_kernel,          # head ↔ backbone lever arm
        angle_harmonic_kernel,         # backbone bending (F6)
        angle_harmonic_kernel,         # head-arm orientation (F6)
        refresh_segment_barbed_kernel,        # live cortex barbed field
        _refresh_bound_walk_dir_kernel,       # bound-head live walk direction
        crossbridge_segment_split_r0bind_kernel,  # attach-unstrained two-array crossbridge (per-head r0_bind)
    ]
    xb_call = next(c for c in launch.calls if c[0] is crossbridge_segment_split_r0bind_kernel)
    xb_inputs = xb_call[2]
    assert xb_inputs[0] is slice_.actuator_state.position_d   # actuator (nmii) position
    assert xb_inputs[1] is slice_.actuator_state.force_d      # actuator (nmii) force  ← +f
    assert xb_inputs[2] is slice_.port.position_d             # cortex-owned segment position
    assert xb_inputs[3] is slice_.port.force_d                # cortex-owned segment force ← −(1−t)f, −t f
    assert xb_inputs[-1] is slice_.connector._r0_bind_d       # per-head zero-strain reference (attach-unstrained)
    assert slice_.actuator_state.force_d is not slice_.port.force_d  # two force arrays never merged


def test_split_crossbridge_is_newton_third_law_and_momentum_conserving() -> None:
    # The landed host oracle is a bit-for-formula twin of crossbridge_segment_split_r0bind_kernel.
    actuator_pos = np.array([[0.10, 0.02, 0.0], [0.0, 0.0, 0.0]])  # head at index 0
    port_pos = np.array([[0.0, 0.0, 0.0], [0.20, 0.0, 0.0], [0.40, 0.0, 0.0]])  # seg nodes 0,1
    walk_dir = np.array([1.0, 0.0, 0.0])
    out = crossbridge_segment_split_r0bind_reference(
        actuator_pos=actuator_pos, port_pos=port_pos, head_node=0, seg_a=0, seg_b=1,
        bary_t=0.25, abscissa=0.03, walk_dir=walk_dir, k_xb=1000.0, r0_bind=-0.01,
    )
    f_head, f_a, f_b = out["f_head"], out["f_seg_a"], out["f_seg_b"]
    # Newton's 3rd law across the two never-merged arrays: the three scatters sum to zero to ~1e-12.
    assert np.linalg.norm(f_head + f_a + f_b) < 1e-12
    # the reaction splits by the barycentric weights (1−t):t = 0.75:0.25 onto the two real cortex nodes.
    assert np.allclose(f_a, -0.75 * f_head, atol=1e-12)
    assert np.allclose(f_b, -0.25 * f_head, atol=1e-12)
    # +f lands on the head (into the actuator array), the reaction onto the cortex nodes (into the port array).
    assert np.allclose(out["actuator_force"][0], f_head, atol=1e-12)
    assert np.allclose(out["port_force"][0], f_a, atol=1e-12)
    assert np.allclose(out["port_force"][1], f_b, atol=1e-12)


def test_crossbridge_r0bind_is_unstrained_at_attach_then_builds_from_power_stroke() -> None:
    # The attach-unstrained fix: capture r0_bind at bind so the passive placement force is 0, regardless of
    # how far (on the coarse mesh) the head is from the nearest segment; the ACTIVE force then develops only
    # from the subsequent power stroke (abscissa).
    for offset in (0.02, 0.05, 0.1, 0.2):                      # head-to-segment tangential offset [µm]
        actuator_pos = np.array([[offset, 0.03, 0.0], [0.0, 0.0, 0.0]])   # head at index 0
        port_pos = np.array([[0.0, 0.0, 0.0], [0.30, 0.0, 0.0]])          # one segment (nodes 0,1)
        walk_dir = np.array([1.0, 0.0, 0.0])
        k_xb, t = 1000.0, 0.0
        # r0_bind captured at the instant of binding (abscissa = 0)
        r0_bind = r0_bind_at_attach_reference(actuator_pos, port_pos, 0, 0, 1, bary_t=t, walk_dir=walk_dir)
        at_attach = crossbridge_segment_split_r0bind_reference(
            actuator_pos, port_pos, 0, 0, 1, bary_t=t, abscissa=0.0, walk_dir=walk_dir, k_xb=k_xb,
            r0_bind=r0_bind)
        # passive force at the binding instant is ZERO to ~1e-9, no matter the coarse-mesh placement offset.
        assert np.linalg.norm(at_attach["f_head"]) < 1e-9, offset
        # after a prescribed relative slide (the power stroke), |load| == k_xb · slide (self-limited by Hill).
        slide = 0.004
        after = crossbridge_segment_split_r0bind_reference(
            actuator_pos, port_pos, 0, 0, 1, bary_t=t, abscissa=slide, walk_dir=walk_dir, k_xb=k_xb,
            r0_bind=r0_bind)
        assert abs(after["load_tangential"] - k_xb * slide) < 1e-9, offset
        assert np.linalg.norm(after["f_head"] + after["f_seg_a"] + after["f_seg_b"]) < 1e-12


# ── 5. commit: accepted-predicated attach → catch-slip detach → epoch KMC ──────────────────────────────
def test_commit_runs_accepted_predicated_attach_catch_slip_detach_epoch() -> None:
    slice_, launch, _copy, _query, _clock, _events = _slice(SegmentDetachKinetics.CATCH_SLIP)
    accepted = _FakeAccepted(accept=True)

    slice_.connector.commit_irreversible(accepted, dt_phys=0.01, rng_seed=131)

    assert launch.kernels == [
        attach_segment_gated_r0bind_kernel,            # attach-UNSTRAINED (captures the per-head r0_bind)
        step_detach_segment_gated_catch_slip_kernel,   # PHYSIOLOGICAL default (Kovacs 2007), NOT Bell slip
        _increment_epoch_if_accepted_kernel,
    ]
    # every KMC kernel is predicated on the scheduler-owned accepted scalar (first input) — reject = no-op.
    for _kernel, _dim, inputs, _outputs, _device in launch.calls:
        assert inputs[0] is accepted
    # the attach kernel flips `bound`, captures r0_bind, and reads the split-ownership positions + attach scratch.
    attach_inputs = launch.calls[0][2]
    assert attach_inputs[2] is slice_.connector.state.bound_d           # the flipped bound flags
    assert attach_inputs[9] is slice_.connector._r0_bind_d              # per-head zero-strain reference (set at bind)
    assert attach_inputs[10] is slice_.connector.state.query_seg_id_d   # candidate segment from propose_events
    assert attach_inputs[15] is slice_.actuator_state.position_d        # actuator (nmii) position (for r0_bind)
    assert attach_inputs[16] is slice_.port.position_d                  # cortex-owned segment position
    assert attach_inputs[17] is slice_.connector.head_node_d
    assert slice_.connector.params in attach_inputs                     # carries k_on for attach_prob
    # the catch-slip detach carries the two-pathway constants (the physiological off-rate).
    detach_inputs = launch.calls[1][2]
    assert slice_.connector._catch_slip in detach_inputs


def test_commit_slip_mode_uses_bell_slip_detach_with_r0bind_attach() -> None:
    slice_, launch, _copy, _query, _clock, _events = _slice(SegmentDetachKinetics.SLIP)
    slice_.connector.commit_irreversible(_FakeAccepted(accept=True), dt_phys=0.01, rng_seed=131)
    assert launch.kernels == [
        attach_segment_gated_r0bind_kernel,            # attach-UNSTRAINED (both detach modes)
        step_detach_segment_gated_kernel,              # pure Bell slip
        _increment_epoch_if_accepted_kernel,
    ]


# ── 6. rejected step bit-restores the binding SoA (snapshot ← live ; live ← snapshot) ──────────────────
def test_rejected_step_bit_restores_connector_binding_soa_and_r0bind() -> None:
    slice_, _launch, copy, _query, _clock, _events = _slice()
    connector = slice_.connector
    # the binding SoA is the base authoritative pairs PLUS the new per-head r0_bind reference.
    pairs = (*connector.state.authoritative_pairs(),
             (connector._r0_bind_d, connector._r0_bind_snap_d))

    connector.snapshot_candidate()
    assert [(dst, src) for dst, src in copy.calls] == [(snap, live) for live, snap in pairs]  # snap ← live

    copy.calls.clear()
    connector.rollback(_FakeAccepted(accept=False))
    assert copy.calls == [(live, snap) for live, snap in pairs]  # rollback: live ← snap (bit-restore, incl. r0_bind)


# ── 7. the clock advances only on an accepted step; every step forwards the SAME predicate ─────────────
def test_clock_advances_only_on_accept() -> None:
    slice_, _launch, _copy, _query, clock, _events = _slice()

    def solve() -> None:
        return None

    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=False))  # rejected
    assert clock.epoch == 0
    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=True))   # accepted
    assert clock.epoch == 1
    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=False))  # rejected again
    assert clock.epoch == 1
    # the clock was asked to advance once per step, always with the step's predicate + dt.
    assert [dt for _pred, dt in clock.advances] == [0.02, 0.02, 0.02]


# ── 8. whole-step order: snapshot → propose → solve → compute-loads → finalize(commit) → advance ───────
def test_step_runs_propose_solve_computeloads_commit_advance_in_order() -> None:
    slice_, launch, _copy, _query, _clock, events = _slice()
    order: list[str] = []

    def solve() -> None:
        order.append("solve")
        events.append(("solve",))

    slice_.step(solve, dt_phys=0.02, accepted_d=_FakeAccepted(accept=True))

    kinds = [e[0] for e in events]
    i_propose = kinds.index("propose")
    i_solve = kinds.index("solve")
    i_compute = next(n for n, e in enumerate(events)
                     if e[0] == "launch" and e[1] is compute_head_loads_segment_split_r0bind_kernel)
    i_attach = next(n for n, e in enumerate(events)
                    if e[0] == "launch" and e[1] is attach_segment_gated_r0bind_kernel)
    i_advance = kinds.index("advance")
    # propose_events before the inner solve; loads computed after solve; commit KMC after loads; clock last.
    assert i_propose < i_solve < i_compute < i_attach < i_advance
    assert order == ["solve"]


# ── 9. the NMII state-owner transaction: reject-gated restore + accepted-only epoch commit ─────────────
def test_actuator_state_transaction_gated_restore_and_epoch_commit() -> None:
    launch = _RecordingLauncher()
    copy = _RecordingCopy()
    actuator = _actuator_state(launch, copy)
    txn = actuator.transaction

    txn.snapshot_candidate()
    # snapshot D2D-copies every covered array (position + ATP cycle/consumed + active population + epoch).
    assert len(copy.calls) == 5
    covered = {src.ptr for _dst, src in copy.calls}
    assert actuator.position_d.ptr in covered and actuator.actuator_epoch_d.ptr in covered

    accepted = _FakeAccepted(accept=False)
    txn.rollback(accepted)
    # position is restored via the vec3 reject-gated kernel; the counts via the int32/int64 gated kernels.
    restore_kernels = [k for k, *_ in launch.calls]
    assert _restore_vec3_if_rejected_kernel in restore_kernels
    assert _restore_int32_if_rejected_kernel in restore_kernels
    assert _restore_int64_if_rejected_kernel in restore_kernels
    for _kernel, _dim, inputs, _outputs, _device in launch.calls:
        assert inputs[0] is accepted  # every restore is device-gated on the predicate

    launch.calls.clear()
    txn.commit_irreversible(_FakeAccepted(accept=True), dt_phys=0.01, rng_seed=9)
    # commit advances only the actuator epoch, gated on the accepted predicate.
    assert [k for k, *_ in launch.calls] == [_increment_epoch_if_accepted_kernel]
    assert launch.calls[0][2][1] is actuator.actuator_epoch_d
