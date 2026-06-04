"""Tests for the Winklbauer single-cell-γ → aggregate-σ bridge (L2 Track-2).

Verifies the σ = β − β* = β(1−φ) relation (Winklbauer 2015, SE311) and the three
literature anchors used to bracket the predicted aggregate spheroid surface tension
from the in-band passive cortical tension g_rigid = 0.57 mN/m.
"""

from __future__ import annotations

import numpy as np

from ffn_sim.scripts.layer2_winklbauer_bridge import (
    BETA_MN_M,
    DAVID_PHI,
    ROFFAY_PHI_WINDOW,
    winklbauer_sigma,
)
from ffn_sim.validation.oracles.spheroid import surface_tension_bridge as br


def test_beta_is_in_band_g_rigid():
    """β = the in-band structural cortical tension (NOT the floored active channel)."""
    assert BETA_MN_M == br.G_RIGID_NATIVE_MN_M == 0.57
    lo, hi = br.CORTICAL_TENSION_BAND_MN_M
    assert lo <= BETA_MN_M <= hi


def test_sigma_eq_gamma_endpoint():
    """φ=0 (max adhesion / Roffay σ=γ identity endpoint) ⇒ σ = β."""
    assert winklbauer_sigma(BETA_MN_M, 0.0) == BETA_MN_M


def test_david_three_quarter_beta():
    """David 2014 typical β*≈¼β ⇒ σ = ¾β, and that mid estimate lands IN the KU band."""
    sig = winklbauer_sigma(BETA_MN_M, DAVID_PHI)
    assert np.isclose(sig, 0.75 * BETA_MN_M)
    lo, hi = br.CORTICAL_TENSION_BAND_MN_M
    assert lo <= sig <= hi


def test_sigma_decreases_with_adhesion_fraction():
    """σ = β(1−φ) is strictly decreasing in φ (more residual contact tension → lower σ)."""
    phis = np.linspace(0.0, 0.625, 20)
    sig = np.array([winklbauer_sigma(BETA_MN_M, float(p)) for p in phis])
    assert np.all(np.diff(sig) < 0.0)
    assert np.all(sig <= BETA_MN_M + 1e-12)   # σ ≤ β always (Winklbauer)


def test_roffay_window_brackets_below_band():
    """Roffay outer/interior ratio [1.6,2.0] → φ∈[0.5,0.625] → σ∈~[0.21,0.29] (below mid)."""
    lo = winklbauer_sigma(BETA_MN_M, ROFFAY_PHI_WINDOW[1])
    hi = winklbauer_sigma(BETA_MN_M, ROFFAY_PHI_WINDOW[0])
    assert 0.20 < lo < 0.23
    assert 0.28 < hi < 0.30
    assert lo < hi < winklbauer_sigma(BETA_MN_M, DAVID_PHI)
