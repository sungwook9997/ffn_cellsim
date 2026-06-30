"""FF → Layer-2 scale bridge (ff/layer2_bridge) — Stage 6i / Decision 3a.

Anchors: the FF single-cell tension splits into a floored actomyosin γ + the state-dependent osmotic
turgor passive γ (~0.2 mN/m at 40 Pa), total positive; and the spheroid A/A0 law a + b/R + c/R²
follows with FF-grounded coefficients (b ∝ γ_cortical, c ∝ β), monotone-decreasing in R — the PI
curvature FORM, grounded in the FF cortex (no fixed/double-booked 0.57 mN/m).
"""

import numpy as np
import pytest

from ffn_sim.ff.gamma_floor import TURGOR_DP0
from ffn_sim.ff.layer2_bridge import aa0_law, ff_single_cell_tension


def test_ff_single_cell_tension_components():
    t = ff_single_cell_tension(n_filaments=60, n_xl=200, n_myo=80, seed=1, n_steps=200)
    assert t.gamma_active_mN_m < 0.05            # actomyosin floored (≪ band)
    assert t.gamma_passive_mN_m == pytest.approx(0.5 * TURGOR_DP0 * 10.0 * 1e-3)  # 40 Pa turgor → 0.2 mN/m
    assert t.gamma_total_mN_m == pytest.approx(t.gamma_active_mN_m + t.gamma_passive_mN_m)
    assert t.gamma_total_mN_m > 0.0


def test_aa0_law_form():
    R = np.linspace(20.0, 120.0, 21)
    aa0 = aa0_law(R, gamma_cortical_mN_m=0.2, beta_mN_m=0.08)
    assert np.all(np.diff(aa0) < 0)              # monotone decreasing in R (a + b/R + c/R²)
    assert aa0[0] > 1.0                          # spreads at small R
    # coefficients scale linearly with the FF inputs (b ∝ γ, c ∝ β)
    aa0_2g = aa0_law(R, gamma_cortical_mN_m=0.4, beta_mN_m=0.08)
    assert (aa0_2g - 1.0)[5] > (aa0 - 1.0)[5]    # doubling γ raises the b/R term
