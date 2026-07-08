"""Sanity gates for Piece #2 — biphasic poroelastic cytoplasm (drainage + drained solid) and Piece #4
(membrane reservoir K_A upturn) added to ``simulate_whole_cell_compression_on_device``.

Gate protocol: dimensional / boundary-limits / fixed-point / sign / backward-compat bit-identity.
Analytic ground truth (no oracle): the Kedem-Katchalsky drainage eigenvalue τ_osm=(V0−vmin)/(Lp·A·Π_in)
and the van't Hoff fixed point (V0_eff=V_cyto ⇒ dP_osm=dP0). Physical anchors + citations in FF_RESULTS_LOG.
"""
import numpy as np
import pytest

from ffn_sim.ff.gamma_floor import (CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex,
                                     TURGOR_PI_IN0, TURGOR_DP0, VMIN_FRAC)
from ffn_sim.ff.network_warp import simulate_whole_cell_compression_on_device
from ffn_sim.common.compartments import resolve_membrane


def _cortex(nf=800, seed=1):
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=nf, n_xl=nf, n_myo=nf // 10,
                                  rng=np.random.default_rng(seed))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx


# ---------- analytic gates (no sim) ----------

def test_drainage_eigenvalue_tau_osm():
    """τ_osm = (V0−vmin)/(Lp·A·Π_in) — the MEMBRANE drainage clock (NOT poroelastic τ_p). Spec A.5:
    at Lp=1e-7 µm/(s·Pa), R=7.5 µm → τ_osm ≈ 34 s; at the exosmotic default 1.6e-8 → ≈ 212 s."""
    R = 7.5
    V0 = (4.0 / 3.0) * np.pi * R ** 3          # 1767 µm³
    A = 4.0 * np.pi * R ** 2                    # 707 µm²
    vmin = VMIN_FRAC * V0
    for Lp, expect in [(1.0e-7, 34.0), (1.6e-8, 212.0)]:
        tau = (V0 - vmin) / (Lp * A * TURGOR_PI_IN0)
        assert tau == pytest.approx(expect, rel=0.05), f"τ_osm(Lp={Lp}) = {tau:.1f} s"


def test_vanthoff_fixed_point_is_resting_turgor():
    """Fixed point: V0_eff → V_cyto ⇒ dP_osm → dP0 exactly (drainage lands on resting turgor, no over-drain)."""
    V0, V_cyto = 1767.0, 1550.0
    vmin = VMIN_FRAC * V0
    V0_eff = V_cyto                              # drained fixed point
    dP = TURGOR_PI_IN0 * (V0_eff - vmin) / (V_cyto - vmin) - (TURGOR_PI_IN0 - TURGOR_DP0)
    assert dP == pytest.approx(TURGOR_DP0, abs=1e-9)


# ---------- sim boundary-limit gates ----------

def test_backward_compat_all_off_is_pure_turgor():
    """All new params off → dP_solid=0, drained_frac=0, and the run is deterministic (same seed → same F)."""
    _, m = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.1,
                                                     n_steps=300, device="cpu")
    assert m["dP_solid_Pa"] == 0.0
    assert m["drained_frac"] == pytest.approx(0.0, abs=1e-9)      # V0_eff == V0
    assert m["K_drained_Pa"] == 0.0
    _, m2 = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.1,
                                                      n_steps=300, device="cpu")
    assert m2["F_plate_pN"] == pytest.approx(m["F_plate_pN"], rel=1e-12)   # deterministic


def test_drainage_relaxes_turgor_and_softens():
    """Lp→large + long load_time (fully drained) → drained_frac→1, dP→~dP0, and F drops vs undrained."""
    _, m_un = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                        n_steps=400, device="cpu")
    _, m_dr = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                        n_steps=400, Lp_um_s_Pa=1.0e-3, load_time_s=1.0e5,
                                                        device="cpu")
    assert m_dr["drained_frac"] > 0.9                       # essentially fully drained
    assert m_dr["dP_turgor_Pa"] < 0.2 * m_un["dP_turgor_Pa"]  # turgor relaxed toward dP0
    assert m_dr["F_plate_pN"] < m_un["F_plate_pN"]           # softer once drained


def test_undrained_at_fast_rate():
    """Fast ramp (short load_time ≪ τ_osm) → the cell has no time to drain → stays ~undrained."""
    _, m = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                     n_steps=400, Lp_um_s_Pa=1.6e-8, load_time_s=0.1,
                                                     device="cpu")
    assert m["drained_frac"] < 0.05                          # negligible drainage on a 0.1 s ramp


def test_drained_solid_carries_load_in_setpoint_limit():
    """Terzaghi solid: with the osmotic term held at the setpoint, K_drained adds dP_solid=K·(V0−V_cyto)/V0>0
    (referenced to the FIXED rest V0 — survives the drained limit; the bug fix)."""
    _, m0 = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                      n_steps=400, pressure_setpoint=40.0, device="cpu")
    _, mk = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                      n_steps=400, pressure_setpoint=40.0, K_drained_Pa=300.0,
                                                      device="cpu")
    assert m0["dP_solid_Pa"] == 0.0
    assert mk["dP_solid_Pa"] > 0.0                           # solid active even at the drained pore-pressure limit
    assert mk["dP_turgor_Pa"] == pytest.approx(40.0 + mk["dP_solid_Pa"], rel=1e-6)
    assert mk["F_plate_pN"] > m0["F_plate_pN"]               # solid stiffens the drained cell (lifts off γ-floor)


# ---------- membrane reservoir (Piece #4) ----------

def test_membrane_f_excess_band_guard():
    """resolve_membrane accepts f_excess in [0,0.40] and rejects outside (Raucher-Sheetz/Figard band)."""
    assert resolve_membrane(f_excess=0.25).f_excess == 0.25
    assert resolve_membrane().f_excess == 0.0                # default = pure plateau (backward-compat)
    with pytest.raises(ValueError):
        resolve_membrane(f_excess=0.6)


# ---------- cortex crosslink turnover (viscoelastic remodeling) ----------

def test_turnover_off_is_backward_compatible():
    """xl_koff_per_s=None → no turnover → deterministic + identical to the no-turnover run."""
    _, m0 = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                      n_steps=400, device="cpu")
    _, m1 = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                      n_steps=400, device="cpu")
    assert m0["xl_turnover_on"] is False
    assert m1["F_plate_pN"] == pytest.approx(m0["F_plate_pN"], rel=1e-12)


def test_turnover_is_rate_dependent_and_softens():
    """Slow ramp (load_time ≫ 1/k_off) → cortex remodels → softer than a fast ramp (load_time ≪ 1/k_off).
    α-actinin k_off0=0.066/s → 1/k_off≈15 s: load_time=600 s fully relaxes, 0.1 s stays elastic."""
    koff = 0.066
    _, m_fast = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                          n_steps=400, xl_koff_per_s=koff, load_time_s=0.1,
                                                          device="cpu")
    _, m_slow = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                          n_steps=400, xl_koff_per_s=koff, load_time_s=600.0,
                                                          device="cpu")
    assert m_fast["xl_turnover_on"] and m_slow["xl_turnover_on"]
    assert m_slow["F_plate_pN"] < m_fast["F_plate_pN"]        # remodeled cortex is softer at slow rate
    # the fast ramp barely turns over → close to the no-turnover baseline
    _, m_base = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.12,
                                                          n_steps=400, device="cpu")
    assert m_fast["F_plate_pN"] == pytest.approx(m_base["F_plate_pN"], rel=0.10)


# ---------- physical press-speed ramp (η-limited deformation → viscous transient) ----------

def test_press_ramp_off_is_backward_compatible():
    """v_press_um_s=None → instant strain (legacy); n_ramp=0, deterministic."""
    _, m = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.05,
                                                     n_steps=300, device="cpu")
    assert m["v_press_um_s"] is None and m["n_ramp"] == 0 and m["t_ramp_s"] == 0.0


def test_press_speed_viscous_transient_and_dwell_relaxes():
    """Physical drag (η=65.9): a fast press reads a NON-equilibrium viscous transient (F grows with v_press),
    and dwelling to relax drops the force far below the ramp-end reading (the press-speed artifact, PI 2026-07-08)."""
    _, m_fast = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.03,
                                                          v_press_um_s=5.0, dwell_steps=0, device="cpu")
    _, m_slow = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.03,
                                                          v_press_um_s=0.5, dwell_steps=0, device="cpu")
    _, m_eq = simulate_whole_cell_compression_on_device(_cortex(), NMIIA_MINIFIL_STALL_PN, strain=0.03,
                                                        v_press_um_s=5.0, dwell_steps=4000, device="cpu")
    assert m_fast["n_ramp"] > 0 and m_fast["dt_real_s"] > 0.0        # physical clock set from η
    assert m_fast["F_plate_pN"] > m_slow["F_plate_pN"]               # faster press → larger viscous transient
    assert m_eq["F_plate_pN"] < 0.1 * m_fast["F_plate_pN"]           # relaxation collapses the transient
