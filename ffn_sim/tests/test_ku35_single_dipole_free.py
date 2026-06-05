"""KU-3.5 assay ladder — GATE 1-FREE: Single Dipole with FREE actin.

GATE 1 (``tests/test_ku35_dipole_gate.py``) is the ISOMETRIC single-stresslet
gate: ONE bipolar minifilament gripping TWO ANTIPARALLEL **rigidly-anchored**
actin filaments. With no integrator the actin cannot move, so it asks only
whether the grip_walk mechanism builds + sustains a contractile *force* dipole.

GATE 1-FREE (this file) removes the anchor: the two antiparallel filaments are
let go and integrated by the **BAOAB (Leimkuhler–Matthews) integrator** at the
**MCF7 physiological cytoplasm viscosity η = 65.9 Pa·s** (Hu 2024; CLAUDE.md
physiological-baseline HARD rule — never water). The question becomes the
minimal *contraction* unit of Miyazaki 2015's minimal reconstitution: when the
bipolar minifilament pulls, do the two free filaments actually MOVE TOWARD each
other? The observable is the inter-filament centre-of-mass SEPARATION; active
(both heads gripping) is compared against a no-grip CONTROL.

Geometry (mirrors GATE 1, û = rod axis = +x, rod centre x = 0):
  * filament 0 (+ side): minus end OUTBOARD at +x, + head bound to its INNERMOST
    bead; grip walks minus-ward so the head pulls the bead INBOARD = contractile.
  * filament 1 (− side): mirror across x = 0.
  ⇒ a complete antiparallel contractile dipole; now the beads are FREE.

HONEST FINDING (what this gate locks in — numbers in the test asserts):
  At η = 65.9 Pa·s the overdamped relaxation time of an actin bead under the
  attach bond is τ = γ/k ≈ 37 s (γ = 6πηR = 3.7·10⁻⁵ N·s/m). A unit-test-
  feasible window (here 400 ticks × batch_steps = 40 000 BAOAB steps ≈ 5·10⁻²
  s of model time) is ~10³× shorter than τ, so the motor moves the actin only
  ~0.1 nm — and the *absolute* COM displacement is BELOW the thermal-noise floor
  (motor signal ≈ −1·10⁻¹¹ m vs control-seed std ≈ 5·10⁻¹¹ m). The minimal free
  dipole is **DRAG-LIMITED**: it holds a clean contractile FORCE but cannot
  produce a displacement-detectable contraction in a test window at physiological
  viscosity. This matches GATE 1 (isometric force, no motion) and the full-cell
  GENERATION-limited verdict — the floor is not the unit's coherence.

  The robust, noise-free way to see the motor IS to use COMMON RANDOM NUMBERS:
  run active + control with the IDENTICAL BAOAB seed so their thermal trajectory
  is bit-identical; the difference (active_Δsep − control_Δsep) isolates the
  motor contribution. That difference is CLEAN, MONOTONE-DECREASING and ≤ 0 at
  EVERY tick (end ≈ −1.1·10⁻¹¹ m) — the motor unambiguously pulls the filaments
  together. The FORCE readout (dipole P) is rock-solid: P_end ≈ −2.1·10⁻¹⁷ N·m,
  contractile, reproducible across seeds to < 0.1 %.

Force constants are the PRODUCTION values; only kinetic RATES are quieted
(k_on, k_off0 → ~0 to isolate the pre-placed grip) and v0 is cranked to exercise
the MECHANISM RATE fast, exactly as the STAGE-1 grip_walk + GATE 1 fixtures do.
No force_scaling / v0-as-force / passive-stiffness tuning; the drag is the real
physiological η, not a convenience.
"""
from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import hoomd
import hoomd.md as md
import yaml

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import (
    MyosinStepUpdater,
    CortexMyosinLayout,
    cortex_myosin_attach_bin_names,
    register_cortex_myosin_bond_params,
    resolve_cortex_myosin,
)
from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.scripts.h3_ku35_stresslet import (
    minifilament_stresslets,
    POL_ANTIPARALLEL,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"

# MCF7 cytoplasm viscosity (Hu 2024; cell/cytoplasm.py table). Physiological
# baseline — NOT water (CLAUDE.md HARD rule). The free actin is integrated in
# this medium.
ETA_MCF7_CYTOPLASM = 65.9          # Pa·s
NB = 7                              # beads per filament
N_TICKS = 400                      # MyosinStepUpdater ticks (each runs batch_steps BAOAB steps)


# --------------------------------------------------------------------------- #
# Config — mirror tests/test_ku35_dipole_gate.py (rates quieted, v0 = rate)
# --------------------------------------------------------------------------- #
def _cfg(grip: bool) -> dict:
    with open(CONFIG_PATH) as f:
        c = yaml.safe_load(f)
    c["cortex"]["n_filaments"] = 50
    c["cortex"]["demo_mode"] = True
    c["cortex"]["myosin"]["n_motors_per_cell"] = 1
    c["cortex"]["myosin"]["n_heads_per_side"] = 1
    if grip:
        c["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    c["cortex"]["myosin"]["head_actin_k_on"] = 1.0e-12     # isolate the pre-placed grip
    c["cortex"]["myosin"]["head_actin_k_off0"] = 1.0e-12
    return c


def _build_free_dipole_fixture(*, v0: float, active: bool, seed: int = 1):
    """ONE bipolar minifilament + TWO antiparallel FREE actin filaments.

    Identical geometry to GATE 1's ``_build_dipole_fixture`` but: (a) the actin
    is NOT pinned — an ``md.Integrator`` (methods=[]) + harmonic bond forces +
    the BAOAB Updater at η = 65.9 Pa·s integrate it; (b) the attach bonds are
    inserted (and the grip hand-bound) ONLY when ``active`` — the ``not active``
    CONTROL is the SAME thermal system with no motor grip.

    Returns ``(sim, upd, baoab_action, p_myo, nb, ell0, gamma_phys)``.
    """
    cfg_full = _cfg(grip=True)
    pc = resolve_h3_derived(cfg_full)
    ell0 = pc.rest_length
    nb = NB
    cfg = deepcopy(cfg_full)
    cfg["cortex"]["myosin"]["v0_per_head"] = float(v0)
    p = resolve_cortex_myosin(cfg, dt=pc.dt_cfl)
    N = p.n_backbone
    gap = ell0                               # filament offset from rod (y)
    perp = 50.0e-9                           # head perp distance to its bound bead

    # filament 0 (+ side): minus end OUTBOARD at +x. bead j at x = (nb - j)*ell0.
    f0 = np.zeros((nb, 3)); f0[:, 0] = (nb - np.arange(nb)) * ell0; f0[:, 1] = +gap
    # filament 1 (− side): mirror across x = 0.
    f1 = np.zeros((nb, 3)); f1[:, 0] = -(nb - np.arange(nb)) * ell0; f1[:, 1] = -gap
    actin_pos = np.vstack([f0, f1])
    n_actin = 2 * nb

    # minifilament backbone along +x centred at x=0; + head INBOARD near fil-0
    # innermost bead, − head mirror.
    seg = p.backbone_segment_length
    bb = np.zeros((N, 3)); bb[:, 0] = (np.arange(N) - (N - 1) / 2.0) * seg
    plus_head = np.array([+ell0, +gap - perp, 0.0])
    minus_head = np.array([-ell0, -gap + perp, 0.0])
    motor_pos = np.vstack([bb, plus_head[None, :], minus_head[None, :]])
    motor_tag_start = n_actin
    all_pos = np.vstack([actin_pos, motor_pos])
    n_total = all_pos.shape[0]

    p_types = ["cortex_actin", "cortex_myosin_backbone", "cortex_myosin_head"]
    typeid = np.zeros(n_total, dtype=np.uint32)
    typeid[motor_tag_start:motor_tag_start + N] = 1
    typeid[motor_tag_start + N:] = 2

    attach_names = cortex_myosin_attach_bin_names(p.n_bins)
    bond_types = ["cortex_actin_backbone", "cortex_myosin_backbone",
                  "cortex_myosin_head_backbone"] + attach_names
    bonds: list[tuple[int, int]] = []
    btids: list[int] = []
    for f in range(2):
        for i in range(nb - 1):
            bonds.append((f * nb + i, f * nb + i + 1)); btids.append(0)
    for i in range(N - 1):
        bonds.append((motor_tag_start + i, motor_tag_start + i + 1)); btids.append(1)
    bonds.append((motor_tag_start + N, motor_tag_start + N // 2)); btids.append(2)
    bonds.append((motor_tag_start + N + 1, motor_tag_start + N // 2)); btids.append(2)
    attach_tid0 = bond_types.index(attach_names[0])
    plus_head_tag = motor_tag_start + N
    minus_head_tag = motor_tag_start + N + 1
    if active:
        # + head -> fil0 innermost bead (nb-1); − head -> fil1 innermost bead.
        bonds.append((plus_head_tag, nb - 1)); btids.append(attach_tid0)
        bonds.append((minus_head_tag, nb + nb - 1)); btids.append(attach_tid0)

    snap = hoomd.Snapshot()
    snap.particles.N = n_total
    snap.particles.types = p_types
    snap.particles.typeid[:] = typeid
    snap.particles.position[:] = all_pos
    snap.particles.velocity[:] = np.zeros((n_total, 3))
    snap.particles.mass[:] = np.ones(n_total)
    box_L = 80.0e-6
    snap.configuration.box = [box_L, box_L, box_L, 0, 0, 0]
    snap.bonds.N = len(bonds)
    snap.bonds.types = bond_types
    snap.bonds.group[:] = np.asarray(bonds, dtype=np.int64)
    snap.bonds.typeid[:] = np.asarray(btids, dtype=np.uint32)

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=seed)
    sim.create_state_from_snapshot(snap)

    # --- forces: harmonic bonds for EVERY bond type present (PRODUCTION k) ----
    bond = md.bond.Harmonic()
    bond.params["cortex_actin_backbone"] = dict(k=pc.bond_k, r0=pc.rest_length)
    register_cortex_myosin_bond_params(bond, p)
    ig = md.Integrator(dt=pc.dt_cfl)
    ig.forces.append(bond)
    sim.operations.integrator = ig                     # methods=[] (BAOAB contract)

    # --- BAOAB at the PHYSIOLOGICAL cytoplasm drag (γ = 6πηR) -----------------
    gamma_phys = 6.0 * math.pi * ETA_MCF7_CYTOPLASM * pc.bead_radius
    gamma_map = {t: gamma_phys for t in p_types}
    baoab_action, baoab_upd = make_baoab_updater(
        kT=pc.kT, gamma=gamma_map, dt=pc.dt_cfl, seed=seed,
    )
    sim.operations.updaters.append(baoab_upd)

    layout = CortexMyosinLayout(
        centers=np.array([[0.0, 0.0, 0.0]]),
        axes=np.array([[1.0, 0.0, 0.0]]),
        positions=motor_pos[None, :, :].copy(),
        motor_tag_start=motor_tag_start,
    )
    cbg = np.asarray(
        [(f * nb + i, f * nb + i + 1) for f in range(2) for i in range(nb - 1)],
        dtype=np.int64,
    )
    upd = MyosinStepUpdater(
        p_myo=p, layout=layout, kT=pc.kT, n_cortex_actin=n_actin,
        cortex_bond_groups=cbg, ell0_cortex=ell0, cortex_beads_per_filament=nb,
    )
    upd._sim_ref = sim
    if active:
        upd._head_bound_to_actin[0] = nb - 1
        upd._head_bound_filament[0] = 0
        upd._head_bound_bead_pos[0] = nb - 1
        upd._head_grip_s[0] = 0.0
        upd._head_bound_to_actin[1] = nb + nb - 1
        upd._head_bound_filament[1] = 1
        upd._head_bound_bead_pos[1] = nb - 1
        upd._head_grip_s[1] = 0.0
    return sim, upd, baoab_action, p, nb, ell0, gamma_phys


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _pos(sim) -> np.ndarray:
    return np.asarray(sim.state.get_snapshot().particles.position)


def _fil_com(pos: np.ndarray, f: int, nb: int) -> np.ndarray:
    return pos[f * nb:(f + 1) * nb].mean(axis=0)


def _separation(sim, nb: int) -> float:
    """Inter-filament centre-of-mass separation (the contraction observable)."""
    pos = _pos(sim)
    return float(np.linalg.norm(_fil_com(pos, 0, nb) - _fil_com(pos, 1, nb)))


def _stresslet(sim, upd, p, nb):
    sl = minifilament_stresslets(
        pos_byTag=_pos(sim), action=upd, p_myo=p, beads_per_filament=nb)
    return sl


def _run_traces(sim, upd, p, nb, *, n_ticks: int):
    """Interleave grip_walk act() with BAOAB sim.run(batch_steps)."""
    sim.run(0)
    sep0 = _separation(sim, nb)
    sep_tr, P_tr, series_tr = [], [], []
    for tick in range(n_ticks):
        upd.act(tick)
        sim.run(p.batch_steps)                 # FREE actin moves under BAOAB
        sep_tr.append(_separation(sim, nb))
        sl = _stresslet(sim, upd, p, nb)
        if sl:
            P_tr.append(sl[0]["dipole_P"])
            series_tr.append(sl[0]["mean_Fseries_over_Fstall"])
    return (sep0, np.asarray(sep_tr),
            np.asarray(P_tr), np.asarray(series_tr))


# --------------------------------------------------------------------------- #
# GATE 1-FREE-a — the fixture builds + runs STABLY under the integrator
# --------------------------------------------------------------------------- #
def test_free_dipole_builds_and_runs_stably():
    """(1) physiological-η BAOAB integrates the free dipole with no blow-up."""
    sim, upd, ba, p, nb, ell0, gphys = _build_free_dipole_fixture(
        v0=0.2, active=True)
    # the drag IS the real physiological η (not water) — guard the setpoint.
    assert gphys == pytest.approx(6.0 * math.pi * 65.9 * 30.0e-9, rel=1e-9)
    sep0, sep_tr, P_tr, series_tr = _run_traces(
        sim, upd, p, nb, n_ticks=N_TICKS)
    pos = _pos(sim)
    assert np.isfinite(pos).all(), "positions went non-finite (integrator blew up)"
    # actin stayed in the box, near its build site (overdamped, drag-limited).
    assert np.linalg.norm(pos[:2 * nb] - _pos_initial(ell0, nb)) < ell0, (
        "free actin drifted more than one ℓ₀ — unstable integration"
    )
    assert len(sep_tr) == N_TICKS
    # the diagnostic still classifies the (now free) dipole as a complete,
    # antiparallel, contractile pair — moving the beads did not break it.
    s = _stresslet(sim, upd, p, nb)[0]
    assert s["complete"] is True
    assert s["polarity"] == POL_ANTIPARALLEL


def _pos_initial(ell0: float, nb: int) -> np.ndarray:
    """The build-site actin positions (for the drift guard)."""
    f0 = np.zeros((nb, 3)); f0[:, 0] = (nb - np.arange(nb)) * ell0; f0[:, 1] = +ell0
    f1 = np.zeros((nb, 3)); f1[:, 0] = -(nb - np.arange(nb)) * ell0; f1[:, 1] = -ell0
    return np.vstack([f0, f1])


# --------------------------------------------------------------------------- #
# GATE 1-FREE-b — the FORCE dipole is contractile + STALL-VALID with free actin
# --------------------------------------------------------------------------- #
def test_free_dipole_force_is_contractile_and_stall_valid():
    """(2)+(3) FREE actin: the dipole P stays contractile (P<0) and the Hill
    SERIES load F_series/F_stall never exceeds 1 (stall-valid) — confirming the
    series cap holds the per-head force physical even when the actin can move."""
    sim, upd, ba, p, nb, ell0, gphys = _build_free_dipole_fixture(
        v0=0.2, active=True)
    _, _, P_tr, series_tr = _run_traces(sim, upd, p, nb, n_ticks=N_TICKS)
    # contractile throughout the back half (no relax-out / sign flip).
    assert (P_tr[-50:] < 0).all(), (
        f"dipole not held contractile with free actin: P_tail max={P_tr[-50:].max():.3e}"
    )
    assert P_tr[-1] < 0
    # the contractile force BUILDS as the grip engages (|P| grows from ~0).
    assert abs(P_tr[-1]) > abs(P_tr[0]) * 1.5
    # STALL-VALID: the Hill series load (the quantity the motor stalls on) is
    # ≤ F_stall at all times (k_series·2ℓ₀ = F_stall, magic-number block).
    assert np.nanmax(series_tr) <= 1.0 + 1e-6, (
        f"series (Hill) load super-stall with free actin: max={np.nanmax(series_tr):.4f}"
    )


def test_free_dipole_P_reproducible_across_seeds():
    """The contractile FORCE readout is rock-solid: P_end is reproducible
    across BAOAB seeds to < 1 % (the force is a clean signal even though the
    DISPLACEMENT is buried in thermal noise — see the next test)."""
    P_ends = []
    for seed in (1, 2, 3, 4):
        sim, upd, ba, p, nb, ell0, gphys = _build_free_dipole_fixture(
            v0=0.2, active=True, seed=seed)
        _, _, P_tr, _ = _run_traces(sim, upd, p, nb, n_ticks=N_TICKS)
        P_ends.append(P_tr[-1])
    P_ends = np.asarray(P_ends)
    assert (P_ends < 0).all(), f"P_end not contractile across seeds: {P_ends}"
    assert abs(P_ends.std() / P_ends.mean()) < 0.01, (
        f"P_end not reproducible across seeds: {P_ends} (cv="
        f"{abs(P_ends.std()/P_ends.mean()):.4f})"
    )


# --------------------------------------------------------------------------- #
# GATE 1-FREE-c — the DISPLACEMENT observable: HONEST drag-limited finding
# --------------------------------------------------------------------------- #
def test_free_dipole_absolute_displacement_is_below_thermal_floor():
    """HONEST NULL (drag-limited): at η = 65.9 Pa·s the motor's ABSOLUTE COM
    displacement over the test window is SMALLER than the thermal-noise floor.
    We lock this in so the limitation is documented, not papered over.

    We compare the active Δsep against the std of the no-grip CONTROL Δsep over
    several thermal seeds. The motor signal (|active mean Δ|) is < the control
    seed-to-seed std → not displacement-distinguishable in this window. (The
    motor IS there — see the common-random-numbers test — it is just drowned by
    Brownian motion at physiological drag in a unit-test-length run.)"""
    ctrl = []
    for seed in (1, 2, 3, 4, 5):
        sim, upd, ba, p, nb, ell0, gphys = _build_free_dipole_fixture(
            v0=0.2, active=False, seed=seed)
        sep0, sep_tr, _, _ = _run_traces(sim, upd, p, nb, n_ticks=N_TICKS)
        ctrl.append(sep_tr[-1] - sep0)
    ctrl = np.asarray(ctrl)
    thermal_floor = ctrl.std()
    assert thermal_floor > 0.0

    act = []
    for seed in (1, 2, 3, 4, 5):
        sim, upd, ba, p, nb, ell0, gphys = _build_free_dipole_fixture(
            v0=0.2, active=True, seed=seed)
        sep0, sep_tr, _, _ = _run_traces(sim, upd, p, nb, n_ticks=N_TICKS)
        act.append(sep_tr[-1] - sep0)
    act = np.asarray(act)

    motor_signal = abs(act.mean() - ctrl.mean())
    # the finding: motor displacement signal is BELOW the thermal floor.
    assert motor_signal < thermal_floor, (
        "unexpected: motor displacement rose ABOVE the thermal floor — the "
        "drag-limited regime changed (signal=%.3e floor=%.3e)"
        % (motor_signal, thermal_floor)
    )
    # and it is sub-nanometre (drag-limited), not a macroscopic contraction.
    assert abs(act.mean()) < 1.0e-9


def test_free_dipole_common_random_numbers_motor_contracts():
    """ROBUST CONTRACTION (the real positive result): with COMMON RANDOM
    NUMBERS — active and control run with the IDENTICAL BAOAB seed so their
    thermal trajectory is bit-identical — the DIFFERENCE isolates the motor.
    That motor-only inter-filament displacement is ≤ 0 at EVERY tick and
    monotone-decreasing: the bipolar minifilament unambiguously pulls the two
    free filaments TOWARD each other (contraction), with thermal noise removed.

    This is the minimal free contraction unit working: separation DECREASES with
    the motor active vs the no-motor control, once the (large, physiological)
    Brownian background is differenced out."""
    seed = 7
    sa, ua, ba_a, p, nb, ell0, _ = _build_free_dipole_fixture(
        v0=0.2, active=True, seed=seed)
    sc, uc, ba_c, _, _, _, _ = _build_free_dipole_fixture(
        v0=0.2, active=False, seed=seed)
    sa.run(0); sc.run(0)
    # common random numbers: same seed ⇒ identical initial BAOAB noise buffers.
    assert np.array_equal(ba_a._prv_rnds, ba_c._prv_rnds)
    sa0 = _separation(sa, nb); sc0 = _separation(sc, nb)
    assert sa0 == pytest.approx(sc0, rel=1e-12)   # identical build geometry

    diff = []
    for tick in range(N_TICKS):
        ua.act(tick); sa.run(p.batch_steps)
        uc.act(tick); sc.run(p.batch_steps)
        diff.append((_separation(sa, nb) - sa0) - (_separation(sc, nb) - sc0))
    diff = np.asarray(diff)

    # (1) motor-only displacement is contractile (≤ 0) at EVERY tick.
    assert (diff <= 1.0e-20).all(), (
        f"motor-only displacement not contractile at some tick: max={diff.max():.3e}"
    )
    # (2) and it BUILDS monotonically (separation keeps shrinking).
    assert (np.diff(diff) <= 1.0e-18).all(), "motor contraction not monotone"
    # (3) the net motor-only approach is a real, negative contraction.
    assert diff[-1] < 0, f"no net contraction: end diff={diff[-1]:.3e}"
    assert diff[-1] == pytest.approx(diff.min(), rel=1e-9)  # end = most contracted
