"""Virtual parallel-plate (AFM) compression — measurement-protocol-consistent γ (FF_STAGE6U).

Smoke-guards `simulate_compressed_shell_on_device`: the code path runs, returns the expected metrics
with finite non-negative values, and the plate half-gap follows the imposed strain. The scientific
magnitude (pressed γ_apparent reaches the tension band at ~3-5% strain vs resting γ_myo ~2e-4 mN/m,
turgor-borne) is validated on the A5000 (FF_STAGE6U), not here — a shape-holding cortex needs N≳1500,
too heavy for a CPU unit test; this guards the interface + sign conventions at small N.
"""

import numpy as np
import pytest

from ffn_sim.ff.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from ffn_sim.ff.network_warp import simulate_compressed_shell_on_device


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


def test_compression_half_gap_scales_with_strain():
    """The plate separation shrinks with strain (0% → touching the resting sphere; 20% → 0.8·R0)."""
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=600, n_xl=600, n_myo=60,
                                  rng=np.random.default_rng(2))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    _, m0 = simulate_compressed_shell_on_device(cx, 0.0, strain=0.0, n_steps=150, device="cpu")
    _, m2 = simulate_compressed_shell_on_device(cx, 0.0, strain=0.2, n_steps=150, device="cpu")
    assert m2["half_gap"] < m0["half_gap"]
    assert m0["half_gap"] == cx.R0_mean
