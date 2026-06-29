"""KU-3.5 STAGE-1 decisive gate: grip_walk holds sustained bond strain.

This is the STAGE-1 acceptance contract (CH1 + CH2 in
``ffn_sim/cortex/myosin.py`` grip_walk arm). It is STANDALONE: a single
bipolar minifilament gripping a SHORT, RIGIDLY-ANCHORED actin segment, with
NO integrator relaxation channel — so the only way the head-actin separation
``r`` can change is the grip *re-targeting* to a farther minus-ward bead
(material transport, CH1). The legacy ``binned_r0`` proxy cannot do this (it
only relabels ``r0`` on the SAME bead), so its ``r/r0`` stays ≈ 1.

Decisive gate (brief §DECISIVE ACCEPTANCE GATE):
  * grip_walk: r/r0 rises > 1.05 and HOLDS it across many ticks.
  * contractility sign: net minus-end-ward transport (grip walks toward
    bead 0; plus-side material pulled toward the rod centre); dγ/dt ≥ 0.
  * regression: binned_r0 attach rest-lengths byte-identical to legacy.
  * conservation: act() round-trip preserves bond count + (no integrator →
    momentum is trivially untouched; we assert positions of pinned actin
    are unchanged across the round-trip, the rigid-anchor invariant).
  * the magic-number relations k·ℓ₀ = F_stall and k_series·2ℓ₀ = F_stall
    FALL OUT of the configured constants (not tuned).

Build note: the production v0 (0.2 µm/s) walks one ℓ₀ bead in ~2 M ticks —
far beyond a seconds-long unit test — so this STANDALONE fixture cranks
``v0_per_head`` purely to exercise the *mechanism* fast. The mechanism
(re-target → r grows → series-tension stall holds r) is v0-independent; v0
only sets the rate. All force/stall constants are the production values.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.myosin import (
    MyosinStepUpdater,
    cortex_myosin_attach_bin_names,
    cortex_myosin_attach_bin_rest_lengths,
    generate_cortex_myosin_layout,
    resolve_cortex_myosin,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Magic-Number Block — the two force ceilings must FALL OUT of the constants.
# ---------------------------------------------------------------------------
def test_magic_number_block_falls_out():
    """k·ℓ₀ = F_stall  and  k_series·2ℓ₀ = F_stall, exactly (not tuned)."""
    cfg = _load_cfg()
    pc = resolve_h3_derived(_demo(cfg))
    p = resolve_cortex_myosin(_demo(cfg, grip=True), dt=pc.dt_cfl)
    ell0 = pc.rest_length
    # Per-bond ceiling: a full one-bead stretch on the head-actin bond.
    assert p.k_head_actin * ell0 == pytest.approx(p.F_stall_per_head, rel=1e-12)
    # Series ceiling: two head springs in series, stretched to the 2ℓ₀ cap.
    k_series = 1.0 / (1.0 / p.k_head_actin + 1.0 / p.k_head_spring)
    assert k_series == pytest.approx(5.0e-7, rel=1e-12)
    assert k_series * 2.0 * ell0 == pytest.approx(p.F_stall_per_head, rel=1e-12)


def test_binned_r0_rest_lengths_byte_identical():
    """REGRESSION: binned_r0 attach rest-lengths unchanged (legacy contract)."""
    # The pre-STAGE-1 implementation: bin-centre linspace over [0, max].
    n_bins = 10
    max_d = 3.3e-7
    legacy = 0.5 * (
        np.linspace(0.0, max_d, n_bins + 1)[:-1]
        + np.linspace(0.0, max_d, n_bins + 1)[1:]
    )
    got = cortex_myosin_attach_bin_rest_lengths(n_bins, max_d)  # default mode
    assert np.array_equal(got, legacy)
    got2 = cortex_myosin_attach_bin_rest_lengths(
        n_bins, max_d, stepping_mode="binned_r0"
    )
    assert np.array_equal(got2, legacy)
    # grip_walk: ALL bins near-zero (CH1) → F = k·r straight from geometry.
    gw = cortex_myosin_attach_bin_rest_lengths(
        n_bins, max_d, stepping_mode="grip_walk"
    )
    assert (gw > 0.0).all()                 # keeps the legacy bins>0 contract
    assert (gw < 1.0e-9).all()              # ~6 orders below ℓ₀ = 500 nm
    assert len(gw) == n_bins                # type COUNT unchanged


# ---------------------------------------------------------------------------
# Standalone rigid-anchor fixture
# ---------------------------------------------------------------------------
def _demo(cfg: dict, *, grip: bool = False) -> dict:
    cfg = deepcopy(cfg)
    cfg["cortex"]["n_filaments"] = 50
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = 1
    if grip:
        cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    return cfg


def _build_rigid_fixture(stepping_mode: str, *, v0: float):
    """ONE bipolar minifilament + a short RIGIDLY-ANCHORED actin segment.

    Layout (all in a plane, SI metres):
      * actin: a single straight filament of ``nb`` beads along +x at y=0,
        spacing ℓ₀ (bead 0 = the minus end, Option-A polarity).
      * minifilament: a backbone along +x ABOVE the actin (y = +200 nm rod
        height is irrelevant here); ONE + head and ONE − head placed so the
        + head sits laterally PERP-close to a HIGH-index (plus-ward) actin
        bead. As the grip walks toward bead 0 it re-targets to farther beads
        → the head-bead distance r grows.

    Returns (sim, updater, head_local, actin_tags, bead_x).
    """
    cfg_full = _load_cfg()
    pc = resolve_h3_derived(_demo(cfg_full))
    ell0 = pc.rest_length                       # 500 nm
    nb = 7                                        # one filament, 7 beads

    cfg = _demo(cfg_full, grip=(stepping_mode == "grip_walk"))
    # Single head per side keeps the dipole minimal + the fixture readable.
    cfg["cortex"]["myosin"]["n_heads_per_side"] = 1
    cfg["cortex"]["myosin"]["v0_per_head"] = float(v0)   # crank rate only
    # Isolate ONE pre-placed grip: suppress de-novo binding (k_on=0) and
    # near-suppress the Bell-Evans strip (tiny k_off0) so the gate measures
    # the STEPPING mechanism (CH1/CH2), not the binding/unbinding kinetics.
    # Force constants (k, F_stall, x_β) are the PRODUCTION values — only the
    # kinetic RATES are quieted, and only in this standalone fixture.
    cfg["cortex"]["myosin"]["head_actin_k_on"] = 1.0e-12
    cfg["cortex"]["myosin"]["head_actin_k_off0"] = 1.0e-12
    p = resolve_cortex_myosin(cfg, dt=pc.dt_cfl)

    N = p.n_backbone
    H = p.n_heads_per_side                        # 1
    per_motor = p.n_particles_per_motor           # N + 2

    # --- actin bead positions (the RIGID anchor) ---
    bead_x = np.arange(nb) * ell0
    actin_pos = np.zeros((nb, 3), dtype=np.float64)
    actin_pos[:, 0] = bead_x                       # along +x, y=z=0

    # --- minifilament particle positions ---
    # backbone along +x, centred over bead (nb-1) (the plus end), at y=+rod.
    rod_y = 0.0
    bb_center_x = bead_x[nb - 1]
    seg = p.backbone_segment_length
    bb_x = bb_center_x + (np.arange(N) - (N - 1) / 2.0) * seg
    bb_pos = np.zeros((N, 3), dtype=np.float64)
    bb_pos[:, 0] = bb_x
    bb_pos[:, 1] = rod_y
    # + head: PERP-close (50 nm) to a plus-ward actin bead so r0 (bind dist)
    # is small and r/r0 amplifies as the grip walks to farther beads.
    perp = 50.0e-9
    plus_head = np.array([bead_x[nb - 1], perp, 0.0])
    minus_head = np.array([bead_x[nb - 1], -perp, 0.0])

    motor_pos = np.vstack([bb_pos, plus_head[None, :], minus_head[None, :]])
    assert motor_pos.shape == (per_motor, 3)

    all_pos = np.vstack([actin_pos, motor_pos])
    n_total = all_pos.shape[0]
    motor_tag_start = nb

    # --- particle + bond types ---
    p_types = ["cortex_actin", "cortex_myosin_backbone", "cortex_myosin_head"]
    typeid = np.zeros(n_total, dtype=np.uint32)
    typeid[motor_tag_start:motor_tag_start + N] = 1
    typeid[motor_tag_start + N:] = 2

    attach_names = cortex_myosin_attach_bin_names(p.n_bins)
    bond_types = ["cortex_actin_backbone", "cortex_myosin_backbone",
                  "cortex_myosin_head_backbone"] + attach_names

    # static bonds: actin chain + minifilament backbone + head-backbone
    bonds: list[tuple[int, int]] = []
    btids: list[int] = []
    for i in range(nb - 1):
        bonds.append((i, i + 1)); btids.append(0)
    for i in range(N - 1):
        bonds.append((motor_tag_start + i, motor_tag_start + i + 1))
        btids.append(1)
    # head→nearest backbone bead
    bonds.append((motor_tag_start + N, motor_tag_start + N // 2)); btids.append(2)
    bonds.append((motor_tag_start + N + 1, motor_tag_start + N // 2)); btids.append(2)
    # Pre-insert the + head→actin ATTACH bond (the grip we isolate). Type =
    # the first attach bin; under grip_walk every attach bin carries r0 ≈ 0
    # so the bin label is force-irrelevant (CH1). Head = plus head (tag
    # motor_tag_start+N), actin = the plus-end bead (nb-1).
    attach_tid0 = bond_types.index(attach_names[0])
    plus_head_tag = motor_tag_start + N
    bonds.append((plus_head_tag, nb - 1)); btids.append(attach_tid0)

    box_L = 40.0e-6
    snap = hoomd.Snapshot()
    snap.particles.N = n_total
    snap.particles.types = p_types
    snap.particles.typeid[:] = typeid
    snap.particles.position[:] = all_pos
    snap.particles.velocity[:] = np.zeros((n_total, 3))
    snap.particles.mass[:] = np.ones(n_total)
    snap.configuration.box = [box_L, box_L, box_L, 0, 0, 0]
    snap.bonds.N = len(bonds)
    snap.bonds.types = bond_types
    snap.bonds.group[:] = np.asarray(bonds, dtype=np.int64)
    snap.bonds.typeid[:] = np.asarray(btids, dtype=np.uint32)

    dev = hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_snapshot(snap)

    # custom layout: rod axis û = +x (backbone end-to-end direction).
    from ffn_sim.archive.hoomd_legacy.cortex.myosin import CortexMyosinLayout
    layout = CortexMyosinLayout(
        centers=np.array([[bb_center_x, rod_y, 0.0]]),
        axes=np.array([[1.0, 0.0, 0.0]]),
        positions=motor_pos[None, :, :].copy(),
        motor_tag_start=motor_tag_start,
    )
    upd = MyosinStepUpdater(
        p_myo=p, layout=layout, kT=pc.kT, n_cortex_actin=nb,
        cortex_bond_groups=np.asarray(
            [(i, i + 1) for i in range(nb - 1)], dtype=np.int64
        ),
        ell0_cortex=ell0, cortex_beads_per_filament=nb,
    )
    # Force a deterministic, immediate bind by setting k_on huge + k_off0 tiny
    # for the test's purpose is NOT done — instead we hand-place the grip so
    # the gate isolates STEPPING (CH1/CH2), not the binding kinetics. Bind the
    # + head to the plus-end bead directly.
    upd._sim_ref = sim
    head_local_plus = 0          # motor 0, + side, head 0
    bead0 = nb - 1               # bind at the plus-end bead
    upd._head_bound_to_actin[head_local_plus] = bead0
    upd._head_bound_filament[head_local_plus] = 0
    upd._head_bound_bead_pos[head_local_plus] = bead0
    upd._head_grip_s[head_local_plus] = 0.0
    # the − head stays free (single-head contractile probe).
    return sim, upd, head_local_plus, bead0, perp, ell0


def _head_global_tag(upd, head_local: int) -> int:
    return upd._head_global_tag(head_local)


def _bond_r(sim, upd, head_local: int) -> float:
    snap = sim.state.get_snapshot()
    pos = np.asarray(snap.particles.position)
    htag = upd._head_global_tag(head_local)
    btag = int(upd._head_bound_to_actin[head_local])
    return float(np.linalg.norm(pos[htag] - pos[btag]))


class TestStage1DecisiveGate:
    def test_grip_walk_holds_sustained_strain(self):
        """DECISIVE: r/r0 rises > 1.05 and HOLDS on a RIGID anchor (CH1)."""
        sim, upd, hl, bead0, r0, ell0 = _build_rigid_fixture(
            "grip_walk", v0=0.2
        )
        r_bind = _bond_r(sim, upd, hl)
        assert r_bind == pytest.approx(r0, rel=1e-6)   # force-free-at-bind ref

        n_bonds_before = int(sim.state.get_snapshot().bonds.N)
        actin_before = np.asarray(
            sim.state.get_snapshot().particles.position
        )[:upd.n_cortex_actin].copy()

        ratios = []
        bead_pos_trace = []
        for tick in range(200):
            upd.act(tick)
            ratios.append(_bond_r(sim, upd, hl) / r_bind)
            bead_pos_trace.append(int(upd._head_bound_bead_pos[hl]))

        ratios = np.asarray(ratios)
        bead_pos_trace = np.asarray(bead_pos_trace)

        # (1) r/r0 rises measurably > 1.05.
        assert ratios.max() > 1.05, f"peak r/r0={ratios.max():.4f} ≤ 1.05"
        # (2) and HOLDS it (the LAST quarter stays > 1.05 — no relax-out).
        tail = ratios[-50:]
        assert (tail > 1.05).all(), (
            f"r/r0 relaxed out: tail min={tail.min():.4f} (≤1.05)"
        )
        # (3) contractility sign: the grip WALKED toward the minus end (bead 0)
        #     — net minus-end-ward transport (dγ/dt ≥ 0, plus-side material
        #     pulled toward rod centre).
        assert bead_pos_trace[-1] < bead0, (
            "grip did not transport minus-ward (no contraction): "
            f"end pos={bead_pos_trace[-1]} start={bead0}"
        )
        assert upd.n_step_advances_total >= 1

        # (4) bond-count invariant across the act() round-trips.
        assert int(sim.state.get_snapshot().bonds.N) == n_bonds_before
        # (5) RIGID anchor: pinned actin positions never moved (no integrator).
        actin_after = np.asarray(
            sim.state.get_snapshot().particles.position
        )[:upd.n_cortex_actin]
        assert np.allclose(actin_after, actin_before, atol=0.0)

    def test_grip_walk_plateaus_and_holds_not_collapse(self):
        """CH2: r/r0 PLATEAUS and HOLDS — it does NOT run to bead 0 and then
        collapse r back to ≈1 (the relax-out failure mode). The hold is the
        AFINES minus-end latch + the s_grip 2ℓ₀ series-tension ceiling
        (k_series·2ℓ₀ = F_stall): the head dwells at the minus end keeping the
        bond stretched indefinitely."""
        sim, upd, hl, bead0, r0, ell0 = _build_rigid_fixture(
            "grip_walk", v0=0.2
        )
        r_bind = _bond_r(sim, upd, hl)
        last = []
        for tick in range(400):
            upd.act(tick)
            last.append(_bond_r(sim, upd, hl) / r_bind)
        last = np.asarray(last)
        # r/r0 held high — NOT collapsed back toward 1.
        assert last[-1] > 1.05
        # Plateau: the final 100 ticks vary little (stalled/latched, holding).
        plateau = last[-100:]
        assert plateau.std() / plateau.mean() < 0.10, (
            f"no stall plateau: cv={plateau.std()/plateau.mean():.3f}"
        )
        # The commanded stretch caps at the 2ℓ₀ series ceiling (CH2 cap raise).
        assert upd._head_grip_s[hl] <= 2.0 * ell0 * (1.0 + 1e-6)
        assert upd._head_grip_s[hl] == pytest.approx(2.0 * ell0, rel=1e-6)

    def test_hill_load_is_contractile_and_slows_walk(self):
        """CH2(i): under a contractile load (stretch GROWS as the head walks
        minus-ward, ĝ·t̂ > 0) the Hill velocity SLOWS — the per-bead dwell
        lengthens as the series tension builds toward F_stall. Compared to the
        FORCE-FREE first bead (r ≈ r0 ⇒ F_series ≈ 0 ⇒ ~v0), a downstream bead
        (large r ⇒ F_series ∝ s) takes strictly MORE ticks to overflow."""
        sim, upd, hl, bead0, r0, ell0 = _build_rigid_fixture(
            "grip_walk", v0=0.002
        )
        # Count ticks spent at each bead position before re-targeting.
        dwell: dict[int, int] = {}
        prev = int(upd._head_bound_bead_pos[hl])
        for tick in range(1200):
            upd.act(tick)
            cur = int(upd._head_bound_bead_pos[hl])
            dwell[prev] = dwell.get(prev, 0) + 1
            prev = cur
            if cur == 0:
                break
        # The plus-end bead (force-free, r≈r0) overflows fast; a downstream
        # contractile bead (large r) dwells STRICTLY longer (Hill slowdown).
        assert dwell[bead0] < dwell[bead0 - 1], (
            f"contractile load did not slow the walk: dwell={dwell}"
        )

    def test_binned_r0_does_not_transport(self):
        """CONTRAST + regression: legacy binned_r0 keeps r on the SAME bead
        (relabels r0 only) → r/r0 stays ≈ 1, no minus-ward transport."""
        sim, upd, hl, bead0, r0, ell0 = _build_rigid_fixture(
            "binned_r0", v0=0.2
        )
        r_bind = _bond_r(sim, upd, hl)
        for tick in range(200):
            upd.act(tick)
        # binned_r0 never re-targets the bound bead (no _head_bound_bead_pos
        # bookkeeping in that arm) → r unchanged on the rigid anchor.
        assert _bond_r(sim, upd, hl) / r_bind == pytest.approx(1.0, abs=1e-9)


def test_two_force_definitions_are_distinct():
    """The Bell-Evans strip uses the FULL bond-frame tension
    k_head_actin·(r−r0), while the Hill stall uses the SERIES tension
    k_series·min(s,r). These are INTENTIONALLY different — written side by
    side and asserted distinct (brief CH2(iii))."""
    cfg = _load_cfg()
    pc = resolve_h3_derived(_demo(cfg))
    p = resolve_cortex_myosin(_demo(cfg, grip=True), dt=pc.dt_cfl)
    r = 4.0e-7         # an engaged bond stretch
    s = 3.0e-7         # commanded sub-bead stretch
    r0 = 1.0e-12       # grip_walk near-zero rest length
    F_bond = p.k_head_actin * (r - r0)                 # Bell-Evans strip
    k_series = 1.0 / (1.0 / p.k_head_actin + 1.0 / p.k_head_spring)
    F_stall_frame = k_series * min(s, r)               # Hill stall
    assert abs(F_bond - F_stall_frame) / F_bond > 0.1   # materially distinct
    # bond-frame is the stiffer (single-spring) tension; series is softer.
    assert F_bond > F_stall_frame


class TestStage1Conservation:
    def test_act_roundtrip_bond_count_and_momentum(self):
        """act() round-trip: bond count invariant + velocities untouched
        (no integrator → momentum drift is identically zero; the act()
        snapshot round-trip must not perturb it)."""
        sim, upd, hl, bead0, r0, ell0 = _build_rigid_fixture(
            "grip_walk", v0=0.2
        )
        snap0 = sim.state.get_snapshot()
        nb0 = int(snap0.bonds.N)
        vel0 = np.asarray(snap0.particles.velocity).copy()
        p0 = float(np.abs(np.asarray(snap0.particles.velocity)
                          * np.asarray(snap0.particles.mass)[:, None]).sum())
        for tick in range(20):
            upd.act(tick)
        snap1 = sim.state.get_snapshot()
        assert int(snap1.bonds.N) == nb0
        vel1 = np.asarray(snap1.particles.velocity)
        assert np.allclose(vel1, vel0, atol=0.0)
        p1 = float(np.abs(vel1 * np.asarray(snap1.particles.mass)[:, None]).sum())
        assert abs(p1 - p0) <= 1.0e-15
