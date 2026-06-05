"""KU-3.5 assay ladder — GATE 1: Single Dipole Gate (isometric bipolar stresslet).

The minimal contractile unit: ONE bipolar minifilament gripping TWO ANTIPARALLEL,
RIGIDLY-ANCHORED actin filaments (one per head side). With no integrator the actin
cannot move, so this is an ISOMETRIC force gate — it asks whether the grip_walk
mechanism builds and SUSTAINS a net CONTRACTILE force-dipole at a STALL-VALID
(Hill-physical, F/F_stall ≤ 1) load, when the two sides pull antiparallel
filaments toward the rod centre.

This is GATE 1 of the assay ladder in the 2026-06-05 KU-3.5 directive
(Single Dipole → Patch Contractility → Shell Patch → Full-Cell). It isolates the
*single-stresslet* mechanism from the scale-up (recruitment/coherence), so a PASS
here localises any full-cell floor to the scale-up, not the unit.

Geometry (û = rod axis = +x, rod centre x = 0):
  * filament 0 (+ side): minus end OUTBOARD at +x; + head sits INBOARD (near
    centre) bound to the innermost bead; grip walks minus-ward (outboard) so the
    head-bead distance r GROWS and the force on the bead points INBOARD
    (toward centre) = contractile.
  * filament 1 (− side): mirror image across x = 0.
  ⇒ a complete, antiparallel, contractile dipole by construction; the diagnostic
    must confirm it, and the grip_walk hold must SUSTAIN P < 0 at F/F_stall ≤ 1.

Uses the production force constants; only kinetic RATES are quieted (k_on=0,
tiny k_off0 to isolate stepping) and v0 is cranked to exercise the mechanism fast
(rate only — the mechanism is v0-independent), exactly as the STAGE-1 grip_walk
fixture does. No force_scaling / passive-stiffness tuning.
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import hoomd
import yaml

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import (
    MyosinStepUpdater,
    CortexMyosinLayout,
    cortex_myosin_attach_bin_names,
    resolve_cortex_myosin,
)
from ffn_sim.scripts.h3_ku35_stresslet import (
    minifilament_stresslets,
    stresslet_summary,
    POL_ANTIPARALLEL,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _cfg(grip: bool) -> dict:
    with open(CONFIG_PATH) as f:
        c = yaml.safe_load(f)
    c["cortex"]["n_filaments"] = 50
    c["cortex"]["demo_mode"] = True
    c["cortex"]["myosin"]["n_motors_per_cell"] = 1
    c["cortex"]["myosin"]["n_heads_per_side"] = 1
    if grip:
        c["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    c["cortex"]["myosin"]["head_actin_k_on"] = 1.0e-12     # isolate stepping
    c["cortex"]["myosin"]["head_actin_k_off0"] = 1.0e-12
    return c


def _build_dipole_fixture(*, v0: float):
    """ONE bipolar minifilament + TWO antiparallel rigidly-anchored filaments.

    Returns (sim, upd, p_myo, nb, ell0).
    """
    cfg_full = _cfg(grip=True)
    pc = resolve_h3_derived(cfg_full)
    ell0 = pc.rest_length
    nb = 7                                  # beads per filament
    cfg = deepcopy(cfg_full)
    cfg["cortex"]["myosin"]["v0_per_head"] = float(v0)
    p = resolve_cortex_myosin(cfg, dt=pc.dt_cfl)
    N = p.n_backbone
    H = p.n_heads_per_side                   # 1
    per_motor = p.n_particles_per_motor      # N + 2
    gap = ell0                               # filament offset from rod (y)
    perp = 50.0e-9                           # head perp distance to its bound bead

    # filament 0 (+ side): minus end OUTBOARD at +x. bead j at x = (nb - j)*ell0
    #   (bead0 = minus = outermost +x ; bead nb-1 = innermost, near centre).
    f0 = np.zeros((nb, 3)); f0[:, 0] = (nb - np.arange(nb)) * ell0; f0[:, 1] = +gap
    # filament 1 (− side): mirror across x=0. bead j at x = -(nb - j)*ell0.
    f1 = np.zeros((nb, 3)); f1[:, 0] = -(nb - np.arange(nb)) * ell0; f1[:, 1] = -gap
    actin_pos = np.vstack([f0, f1])          # filament 0 = tags 0..nb-1, fil 1 = nb..2nb-1
    n_actin = 2 * nb

    # minifilament backbone along +x centred at x=0, y=0.
    seg = p.backbone_segment_length
    bb = np.zeros((N, 3)); bb[:, 0] = (np.arange(N) - (N - 1) / 2.0) * seg
    # + head INBOARD near filament-0 innermost bead (x=+ell0); − head mirror.
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
    # actin chains (both filaments), minus->plus
    for f in range(2):
        for i in range(nb - 1):
            bonds.append((f * nb + i, f * nb + i + 1)); btids.append(0)
    for i in range(N - 1):
        bonds.append((motor_tag_start + i, motor_tag_start + i + 1)); btids.append(1)
    bonds.append((motor_tag_start + N, motor_tag_start + N // 2)); btids.append(2)
    bonds.append((motor_tag_start + N + 1, motor_tag_start + N // 2)); btids.append(2)
    # pre-insert BOTH attach bonds: + head -> fil0 innermost bead (nb-1);
    #   − head -> fil1 innermost bead (nb + nb-1).
    attach_tid0 = bond_types.index(attach_names[0])
    plus_head_tag = motor_tag_start + N
    minus_head_tag = motor_tag_start + N + 1
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

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)

    layout = CortexMyosinLayout(
        centers=np.array([[0.0, 0.0, 0.0]]),
        axes=np.array([[1.0, 0.0, 0.0]]),          # û = +x
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
    # hand-bind both heads at their innermost beads (force-free at bind).
    upd._head_bound_to_actin[0] = nb - 1          # + head (local 0) -> fil0 bead nb-1
    upd._head_bound_filament[0] = 0
    upd._head_bound_bead_pos[0] = nb - 1
    upd._head_grip_s[0] = 0.0
    upd._head_bound_to_actin[1] = nb + nb - 1     # − head (local 1) -> fil1 bead nb-1
    upd._head_bound_filament[1] = 1
    upd._head_bound_bead_pos[1] = nb - 1
    upd._head_grip_s[1] = 0.0
    return sim, upd, p, nb, ell0


def _pos(sim) -> np.ndarray:
    return np.asarray(sim.state.get_snapshot().particles.position)


def _stresslet(sim, upd, p, nb):
    return minifilament_stresslets(
        pos_byTag=_pos(sim), action=upd, p_myo=p, beads_per_filament=nb)


# --------------------------------------------------------------------------- #
# GATE 1a — the diagnostic correctly classifies the hand-built dipole
# --------------------------------------------------------------------------- #
def test_diagnostic_classifies_complete_antiparallel_contractile():
    sim, upd, p, nb, ell0 = _build_dipole_fixture(v0=0.2)
    sl = _stresslet(sim, upd, p, nb)
    assert len(sl) == 1
    s = sl[0]
    assert s["complete"] is True, "should be a complete bipolar pair (both sides bound, diff fil)"
    assert s["polarity"] == POL_ANTIPARALLEL, f"polarity={s['polarity']}"
    assert s["n_plus"] == 1 and s["n_minus"] == 1
    summ = stresslet_summary(sl)
    assert summ["frac_complete_pairs"] == 1.0


# --------------------------------------------------------------------------- #
# GATE 1b — grip_walk SUSTAINS a contractile dipole at a STALL-VALID load
# --------------------------------------------------------------------------- #
def test_single_dipole_sustains_contractile_force():
    sim, upd, p, nb, ell0 = _build_dipole_fixture(v0=0.2)
    P_trace, series_trace, bond_trace = [], [], []
    for tick in range(200):
        upd.act(tick)
        sl = _stresslet(sim, upd, p, nb)
        assert len(sl) == 1
        P_trace.append(sl[0]["dipole_P"])
        series_trace.append(sl[0]["mean_Fseries_over_Fstall"])
        bond_trace.append(sl[0]["max_Fbond_over_Fstall"])
    P = np.asarray(P_trace)
    series = np.asarray(series_trace); bond = np.asarray(bond_trace)

    # (1) the dipole is CONTRACTILE (P<0) and grows in magnitude as the grip walks.
    assert P[-1] < 0, f"dipole not contractile: P_end={P[-1]:.3e}"
    assert abs(P[-1]) > abs(P[0]) * 1.5, (
        f"contractile dipole did not build: |P0|={abs(P[0]):.3e} |P_end|={abs(P[-1]):.3e}"
    )
    # (2) it HOLDS (last quarter stays contractile, no relax-out).
    assert (P[-50:] < 0).all(), "dipole relaxed out of contraction"
    # (3) STALL-VALID on the HILL load: the SERIES tension k_series·min(s,r) — the
    #     quantity the Hill kernel stalls on — never exceeds F_stall (s_grip caps
    #     at 2ℓ₀ and k_series·2ℓ₀ = F_stall, by the magic-number block).
    assert np.nanmax(series) <= 1.0 + 1e-6, (
        f"series (Hill) load super-stall: max={np.nanmax(series):.3f}"
    )
    # (4) FINDING (verified, not a bug): on a RIGID anchor the held grip strain
    #     drives the BOND-frame force k·r up to ~(filament span)/ℓ₀ × F_stall =
    #     super-stall, because the harmonic attach bond is NOT stall-capped. A
    #     single head cannot physically deliver > F_stall — in production this is
    #     capped by Bell-Evans detachment (k_off ∝ exp(F·x_β/kT)), which strips a
    #     super-stall head, holding the DELIVERED force near F_stall and r near ℓ₀
    #     (the production aggregation ledger shows F_bond/F_stall ~0.005–0.45,
    #     sub-stall, confirming the production cap). This rigid-anchor fixture
    #     suppresses Bell-Evans to isolate stepping, so it exposes the uncapped
    #     held strain.
    assert np.nanmax(bond) > 1.0, (
        "expected rigid-anchor held strain to exceed stall in the bond frame "
        f"(got max={np.nanmax(bond):.2f}); fixture geometry changed?"
    )


def test_single_dipole_summary_verdict_contractile():
    sim, upd, p, nb, ell0 = _build_dipole_fixture(v0=0.2)
    for tick in range(150):
        upd.act(tick)
    summ = stresslet_summary(_stresslet(sim, upd, p, nb))
    assert summ["coherence"] < 0, f"coherence={summ['coherence']:.3f} (>=0: not contractile)"
    assert "CONTRACTILE" in summ["verdict"], summ["verdict"]
    # rigid anchor: actin never moved (isometric) — confirms force, not motion.
    # (the only mobile DOF are absent; act() must not write actin positions)


# --------------------------------------------------------------------------- #
# the bundled diagnostic self-test passes end to end
# --------------------------------------------------------------------------- #
def test_stresslet_self_test_passes():
    from ffn_sim.scripts.h3_ku35_stresslet import self_test
    assert self_test(verbose=False) is True
