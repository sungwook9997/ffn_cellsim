"""Virtual parallel-plate (AFM) compression — measurement-protocol-consistent γ (FF_STAGE6U).

Smoke-guards `simulate_compressed_shell_on_device`: the code path runs, returns the expected metrics
with finite non-negative values, and the plate half-gap follows the imposed strain. The scientific
magnitude (pressed γ_apparent reaches the tension band at ~3-5% strain vs resting γ_myo ~2e-4 mN/m,
turgor-borne) is validated on the A5000 (FF_STAGE6U), not here — a shape-holding cortex needs N≳1500,
too heavy for a CPU unit test; this guards the interface + sign conventions at small N.
"""

import numpy as np
import pytest

from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from aleph.laws.network_warp import (
    simulate_compressed_shell_on_device,
    simulate_whole_cell_compression_on_device,
)
from aleph.laws.compartments import resolve_membrane, resolve_nucleus


def test_compressed_shell_runs_and_returns_finite_metrics():
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=800, n_xl=800, n_myo=80,
                                  rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    _, m = simulate_compressed_shell_on_device(cx, NMIIA_MINIFIL_STALL_PN, strain=0.1,
                                               n_steps=300, device="cpu")
    for k in ("strain", "half_gap", "F_plate_pN", "R_eq_um", "contact_radius_um",
              "dP_turgor_Pa", "V_over_V0", "gamma_apparent_pN_um", "gamma_apparent_mN_m"):
        assert k in m and np.isfinite(m[k])
    assert m["F_plate_pN"] >= 0.0                         # plate reaction is a magnitude
    assert m["gamma_apparent_mN_m"] >= 0.0
    assert m["half_gap"] == cx.R0_mean * (1.0 - 0.1)     # plate follows the imposed strain
    assert 0.5 < m["V_over_V0"] < 1.5                     # volume stays physical (turgor holds)


def test_pressure_setpoint_regulation_holds_dP():
    """Volume-regulation mode (pressure_setpoint) holds ΔP at the setpoint (perfect water-flux limit),
    giving a confinement-INDEPENDENT γ_apparent = ΔP_setpoint·R_eq/2 — matching the experimental
    Fischer-Friedrich confinement-independence (FF_STAGE6V). Here just guard that the setpoint is honored."""
    import numpy as np
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=800, n_xl=800, n_myo=80,
                                  rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    _, m = simulate_compressed_shell_on_device(cx, NMIIA_MINIFIL_STALL_PN, strain=0.1,
                                               n_steps=300, pressure_setpoint=400.0, device="cpu")
    assert abs(m["dP_turgor_Pa"] - 400.0) < 1e-6            # ΔP regulated to the setpoint
    assert m["gamma_apparent_mN_m"] == pytest.approx(0.5 * 400.0 * m["R_eq_um"] * 1e-3, rel=1e-6)


def test_whole_cell_nucleus_engagement_geometry():
    """Whole-cell AFM (cortex+turgor+shared nucleus): the plate reaches the nucleus only at strain
    > 1−R_nuc/R_cell (≈67% for R_nuc/R_cell=1/3), so the nucleus is mechanically absent from the
    cortical-tension regime (<50% strain) and engages only under deep compression — the geometric
    finding behind 'the nucleus does not change the cortical-tension AFM' (FF whole-cell wiring)."""
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=900, n_xl=900, n_myo=90,
                                  rng=np.random.default_rng(1))
    R0 = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    cx.R0_mean = R0
    nuc = resolve_nucleus(R_nuc_um=R0 / 3.0, n_beads=60)      # R_nuc/R_cell = 1/3 → contact at ~67% strain
    # low strain (cortical-tension regime): plate does NOT reach the nucleus
    _, m_lo = simulate_whole_cell_compression_on_device(cx, NMIIA_MINIFIL_STALL_PN, strain=0.1,
                                                        nucleus=nuc, n_steps=300, device="cpu")
    assert m_lo["n_nuc"] == 60
    assert m_lo["nucleus_contact"] is False
    # deep strain (past 1−R_nuc/R_cell=0.667): the plate compresses the nucleus directly
    cx2 = build_crosslinked_cortex(CortexParams(), n_filaments=900, n_xl=900, n_myo=90,
                                   rng=np.random.default_rng(1))
    cx2.R0_mean = float(np.linalg.norm(cx2.net.pos - cx2.net.pos.mean(0), axis=1).mean())
    _, m_hi = simulate_whole_cell_compression_on_device(cx2, NMIIA_MINIFIL_STALL_PN, strain=0.7,
                                                        nucleus=nuc, n_steps=300, device="cpu")
    assert m_hi["nucleus_contact"] is True


def test_whole_cell_membrane_is_small_additive_inward_tension():
    """The plasma-membrane channel (buffered-plateau γ_mem) adds a small INWARD tension: it holds γ_mem at the
    operating point, lowers F_plate by a few % (bears a share of the surface tension), and leaves the
    turgor-Laplace γ_apparent unchanged — i.e. the membrane, like the nucleus, does not change the cortical-
    tension measurement. Stable (no divergence)."""
    mem = resolve_membrane()
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=900, n_xl=900, n_myo=90,
                                  rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    _, m0 = simulate_whole_cell_compression_on_device(cx, NMIIA_MINIFIL_STALL_PN, strain=0.1, membrane=None,
                                                      n_steps=400, pressure_setpoint=40.0, device="cpu")
    cx2 = build_crosslinked_cortex(CortexParams(), n_filaments=900, n_xl=900, n_myo=90,
                                   rng=np.random.default_rng(1))
    cx2.R0_mean = float(np.linalg.norm(cx2.net.pos - cx2.net.pos.mean(0), axis=1).mean())
    _, mm = simulate_whole_cell_compression_on_device(cx2, NMIIA_MINIFIL_STALL_PN, strain=0.1, membrane=mem,
                                                      n_steps=400, pressure_setpoint=40.0, device="cpu")
    assert mm["has_membrane"] is True and m0["has_membrane"] is False
    assert mm["gamma_mem_channel_pN_um"] == pytest.approx(10.0)          # buffered baseline present
    assert 0.0 < mm["dP_mem_Pa"] < 6.0                                   # small inward pressure (~2.7 Pa)
    assert np.isfinite(mm["F_plate_pN"]) and mm["F_plate_pN"] > 0.0      # stable
    assert mm["F_plate_pN"] <= m0["F_plate_pN"] * (1.0 + 1e-3)           # inward channel: does NOT raise F_plate
    assert mm["gamma_apparent_mN_m"] == pytest.approx(m0["gamma_apparent_mN_m"], rel=5e-3)  # γ unchanged (robust)


def test_rigid_plate_confines_cortex_exactly():
    """The RIGID plate (hard z-clamp) confines the cortex EXACTLY to the gap — no node floats above the plate
    (the soft penalty under-confines a stiff cortex: PI 2026-07-02 caught the cortex poking through). Reaction
    is finite (the AFM force = removed-displacement/dt)."""
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=1000, n_xl=1000, n_myo=100,
                                  rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    pos, m = simulate_whole_cell_compression_on_device(cx, NMIIA_MINIFIL_STALL_PN, strain=0.35,
                                                       n_steps=1500, pressure_setpoint=40.0,
                                                       rigid_plate=True, device="cpu")
    cor = pos[:cx.net.n_nodes]
    zc = np.abs(cor[:, 2] - cor[:, 2].mean())
    assert zc.max() <= m["half_gap"] + 1e-6            # NO cortex node above the plate (exact confinement)
    assert m["F_plate_pN"] >= 0.0 and np.isfinite(m["F_plate_pN"])


def test_whole_cell_none_nucleus_matches_cortex_only_interface():
    """With nucleus=None the whole-cell path reduces to the cortex-only AFM (same metric keys, n_nuc=0)."""
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=800, n_xl=800, n_myo=80,
                                  rng=np.random.default_rng(2))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    _, m = simulate_whole_cell_compression_on_device(cx, NMIIA_MINIFIL_STALL_PN, strain=0.1,
                                                     nucleus=None, n_steps=200, device="cpu")
    assert m["n_nuc"] == 0 and m["nucleus_contact"] is False
    for k in ("F_plate_pN", "gamma_apparent_mN_m", "dP_turgor_Pa", "V_over_V0"):
        assert np.isfinite(m[k])


def test_compression_half_gap_scales_with_strain():
    """The plate separation shrinks with strain (0% → touching the resting sphere; 20% → 0.8·R0)."""
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=600, n_xl=600, n_myo=60,
                                  rng=np.random.default_rng(2))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    _, m0 = simulate_compressed_shell_on_device(cx, 0.0, strain=0.0, n_steps=150, device="cpu")
    _, m2 = simulate_compressed_shell_on_device(cx, 0.0, strain=0.2, n_steps=150, device="cpu")
    assert m2["half_gap"] < m0["half_gap"]
    assert m0["half_gap"] == cx.R0_mean
