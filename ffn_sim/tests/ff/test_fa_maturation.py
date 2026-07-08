r"""Focal-adhesion maturation (ff/fa_maturation) — talin unfolding + vinculin + Hill FA growth.

FF-native validation of the KB-2.x mechanosensor layer that turns a bare integrin clutch into a mature,
reinforced, force-growing focal adhesion. Every gate is an ANALYTIC ground truth (oracle-is-crosscheck):

  1. TALIN Bell parity (KB-2.6) — ``k_unfold(F) = k_u0·exp(F·dx/kBT)`` (force ENHANCES unfolding, del Rio 2009);
     the Warp forward step matches the closed form; refolding is Bell-suppressed ``k_refold·exp(−F·dx/kBT)``.
  2. Hill FIXED POINT (KB-2.17) — net FA area growth ``dA/dt = (k_g(F) − k_d)·A`` vanishes EXACTLY at F = F_th,
     because the KB constraint ``k_g0 = 2·k_d`` puts the Hill half-max at threshold; grows above, disassembles below.
  3. VINCULIN steady state (KB-2.7) — ``dN/dt = k_rec·p·(N_max−N) − k_diss·N`` relaxes to
     ``N* = k_rec·p·N_max/(k_rec·p + k_diss)``, reinforcing ``k_int_eff = k_int·(1 + α·N)``.
  4. clutch load + force-gated disassembly.

⚠️ DEBUG RECORD (PI 2026-07-08): this maturation runs downstream of the ff/fa_clutch catch-slip bond, whose
peak is F*≈7 pN from the recorded Pereverzev params while the KB / Kong-2009 report ~30 pN — surfaced to PI,
kept as-recorded (NOT retuned to a target). Flagged here because the FA-growth threshold F_th=5 pN and the
clutch peak F* sit on the SAME force axis: when the absolute traction magnitude is later reconciled, both move
together. See ff/fa_clutch_warp.py head + tests/ff/test_fa_clutch.py::test_resolve_and_Fstar.
"""

import numpy as np
import pytest
import warp as wp

from ffn_sim.ff.units import KBT
from ffn_sim.ff.fa_maturation import (
    FA_FTH_PN,
    FA_KD,
    FA_KG0,
    FA_N,
    TALIN_DX_UM,
    TALIN_KREFOLD,
    TALIN_KU0,
    VIN_ALPHA,
    VIN_KDISS,
    VIN_KREC,
    VIN_NMAX,
    MaturationParams,
    fa_disassemble_kernel,
    fa_growth_kernel,
    fa_growth_rate,
    clutch_load_kernel,
    k_int_effective,
    talin_unfold_rate,
    talin_vinculin_kernel,
)

D = "cpu"


# ----------------------------------------------------------------------------- talin unfolding (Bell)
def test_talin_bell_closed_form():
    """k_unfold(F) = k_u0·exp(F·dx/kBT); force enhances (del Rio 2009). Monotone increasing."""
    for F in [0.0, 2.5, 5.0, 10.0, 20.0]:
        assert talin_unfold_rate(F) == pytest.approx(TALIN_KU0 * np.exp(F * TALIN_DX_UM / KBT), rel=1e-12)
    assert talin_unfold_rate(0.0) == pytest.approx(TALIN_KU0)
    rates = [talin_unfold_rate(F) for F in (0.0, 5.0, 10.0, 20.0)]
    assert all(b > a for a, b in zip(rates, rates[1:]))


def test_talin_kernel_forward_step_matches_bell():
    """One Warp step from p=0 advances p_unf by dt·k_unfold(F) (Bell forward); the kernel reproduces the ODE."""
    dt, F = 1.0e-3, 8.0                                        # tiny dt so the linearization p≈dt·kf is exact
    load = wp.array([F], dtype=wp.float64, device=D)
    p_unf = wp.zeros(1, dtype=wp.float64, device=D)
    n_vin = wp.zeros(1, dtype=wp.float64, device=D)
    bound = wp.array([1], dtype=wp.int32, device=D)
    p = MaturationParams()
    wp.launch(talin_vinculin_kernel, dim=1,
              inputs=[load, p_unf, n_vin, bound, wp.float64(p.ku0), wp.float64(p.dx_um), wp.float64(p.kT),
                      wp.float64(p.k_refold), wp.float64(p.k_rec), wp.float64(p.k_diss), wp.float64(p.n_max),
                      wp.float64(dt)], device=D)
    kf = TALIN_KU0 * np.exp(F * TALIN_DX_UM / KBT)
    assert p_unf.numpy()[0] == pytest.approx(dt * kf, rel=1e-9)
    # vinculin uses the UPDATED p_unf: one step from N=0 ⇒ dt·k_rec·p·N_max
    assert n_vin.numpy()[0] == pytest.approx(dt * VIN_KREC * (dt * kf) * VIN_NMAX, rel=1e-9)


def test_talin_unbound_refolds():
    """A detached clutch feels no load ⇒ p_unf decays by the zero-force refolding rate k_refold."""
    dt = 1.0e-3
    load = wp.array([50.0], dtype=wp.float64, device=D)        # load present but bound=0 ⇒ F forced to 0 in-kernel
    p_unf = wp.array([1.0], dtype=wp.float64, device=D)
    n_vin = wp.zeros(1, dtype=wp.float64, device=D)
    bound = wp.array([0], dtype=wp.int32, device=D)
    p = MaturationParams()
    wp.launch(talin_vinculin_kernel, dim=1,
              inputs=[load, p_unf, n_vin, bound, wp.float64(p.ku0), wp.float64(p.dx_um), wp.float64(p.kT),
                      wp.float64(p.k_refold), wp.float64(p.k_rec), wp.float64(p.k_diss), wp.float64(p.n_max),
                      wp.float64(dt)], device=D)
    # p=1,F=0 ⇒ dp = dt·(ku0·(1−1) − k_refold·1) = −dt·k_refold
    assert p_unf.numpy()[0] == pytest.approx(1.0 - dt * TALIN_KREFOLD, rel=1e-9)


# ----------------------------------------------------------------------------- Hill FA growth fixed point
def test_hill_self_consistency_and_fixed_point():
    """The KB constraint k_g0 = 2·k_d places the Hill half-max at F_th ⇒ dA/dt = 0 exactly at F = F_th."""
    assert FA_KG0 == pytest.approx(2.0 * FA_KD)                # self-consistent (fixed point AT threshold)
    assert fa_growth_rate(FA_FTH_PN) == pytest.approx(0.0, abs=1e-15)
    assert fa_growth_rate(2.0 * FA_FTH_PN) > 0                 # above threshold → grows
    assert fa_growth_rate(0.5 * FA_FTH_PN) < 0                 # below threshold → disassembles
    # explicit Hill form
    F = 12.0
    kg = FA_KG0 * F**FA_N / (F**FA_N + FA_FTH_PN**FA_N)
    assert fa_growth_rate(F) == pytest.approx(kg - FA_KD, rel=1e-12)


def _step_growth(F, area0, dt, steps):
    area = wp.array([area0], dtype=wp.float64, device=D)
    load = wp.array([F], dtype=wp.float64, device=D)
    for _ in range(steps):
        wp.launch(fa_growth_kernel, dim=1,
                  inputs=[load, area, wp.float64(FA_KG0), wp.float64(FA_KD), wp.float64(FA_N),
                          wp.float64(FA_FTH_PN), wp.float64(dt)], device=D)
    return float(area.numpy()[0])


def test_fa_growth_kernel_parity_and_direction():
    """fa_growth_kernel one step = A·(1+dt·(k_g(F)−k_d)); grows above F_th, shrinks below, stationary AT F_th."""
    dt = 1.0e-2
    F = 15.0
    kg = FA_KG0 * F**FA_N / (F**FA_N + FA_FTH_PN**FA_N)
    assert _step_growth(F, 1.0, dt, 1) == pytest.approx(1.0 * (1.0 + dt * (kg - FA_KD)), rel=1e-12)
    assert _step_growth(3.0 * FA_FTH_PN, 1.0, 1.0, 50) > 1.0            # sustained supra-threshold → grows
    assert _step_growth(0.3 * FA_FTH_PN, 1.0, 1.0, 50) < 1.0            # sustained sub-threshold → shrinks
    assert _step_growth(FA_FTH_PN, 1.0, 1.0, 50) == pytest.approx(1.0, abs=1e-9)   # fixed point holds


def test_fa_disassembles_below_area_floor():
    """A sustained sub-threshold FA shrinks below a_min and the disassembly kernel unbinds it (force-gated)."""
    area = wp.array([0.05], dtype=wp.float64, device=D)        # already below a_min
    bound = wp.array([1, ], dtype=wp.int32, device=D)
    wp.launch(fa_disassemble_kernel, dim=1, inputs=[area, bound, wp.float64(0.1)], device=D)
    assert bound.numpy()[0] == 0 and area.numpy()[0] == pytest.approx(1.0)   # unbound + reset for a future nascent
    # a healthy FA above the floor stays bound
    area2 = wp.array([0.5], dtype=wp.float64, device=D); bound2 = wp.array([1], dtype=wp.int32, device=D)
    wp.launch(fa_disassemble_kernel, dim=1, inputs=[area2, bound2, wp.float64(0.1)], device=D)
    assert bound2.numpy()[0] == 1


# ----------------------------------------------------------------------------- vinculin reinforcement
def test_vinculin_steady_state_and_reinforcement():
    """N_vin relaxes to the analytic steady state N* = k_rec·p·N_max/(k_rec·p + k_diss); reinforcement is
    k_int_eff = k_int·(1 + α·N)."""
    dt, F, steps = 0.5, 25.0, 8000                            # >> both talin (1/kf..kr) and vinculin (1/k_diss) times
    load = wp.array([F], dtype=wp.float64, device=D)
    p_unf = wp.zeros(1, dtype=wp.float64, device=D)
    n_vin = wp.zeros(1, dtype=wp.float64, device=D)
    bound = wp.array([1], dtype=wp.int32, device=D)
    p = MaturationParams()
    for _ in range(steps):
        wp.launch(talin_vinculin_kernel, dim=1,
                  inputs=[load, p_unf, n_vin, bound, wp.float64(p.ku0), wp.float64(p.dx_um), wp.float64(p.kT),
                          wp.float64(p.k_refold), wp.float64(p.k_rec), wp.float64(p.k_diss), wp.float64(p.n_max),
                          wp.float64(dt)], device=D)
    kf = TALIN_KU0 * np.exp(F * TALIN_DX_UM / KBT)
    kr = TALIN_KREFOLD * np.exp(-F * TALIN_DX_UM / KBT)
    p_star = kf / (kf + kr)
    n_star = VIN_KREC * p_star * VIN_NMAX / (VIN_KREC * p_star + VIN_KDISS)
    assert p_unf.numpy()[0] == pytest.approx(p_star, rel=1e-4)
    assert n_vin.numpy()[0] == pytest.approx(n_star, rel=1e-4)
    assert k_int_effective(1000.0, n_star) == pytest.approx(1000.0 * (1.0 + VIN_ALPHA * n_star), rel=1e-12)


# ----------------------------------------------------------------------------- clutch load readout
def test_clutch_load_is_spring_tension():
    """Per-clutch mechanosensor input F = k_int·|L − rest| [pN]; a detached clutch reports 0."""
    L, rest, k_int = 0.18, 0.05, 1000.0
    pos = wp.array(np.array([[L, 0.0, 0.0], [L, 0.0, 0.0]]), dtype=wp.vec3d, device=D)
    ac = wp.array([0, 1], dtype=wp.int32, device=D)
    anchor = wp.array(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]), dtype=wp.vec3d, device=D)
    bound = wp.array([1, 0], dtype=wp.int32, device=D)
    load = wp.zeros(2, dtype=wp.float64, device=D)
    wp.launch(clutch_load_kernel, dim=2,
              inputs=[pos, ac, anchor, bound, wp.float64(k_int), wp.float64(rest), load], device=D)
    assert load.numpy()[0] == pytest.approx(k_int * abs(L - rest))
    assert load.numpy()[1] == 0.0
