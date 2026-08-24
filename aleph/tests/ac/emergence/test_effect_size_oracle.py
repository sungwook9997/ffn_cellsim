"""Self-test of the finite-N isotropic null band + the pre-registered effect-size rule (pure NumPy).

Anchors the Monte-Carlo null against the closed form ``E[sum lambda^2] = 3/(2N)``, checks the
``1/sqrt(N)`` scaling of the null bias, and verifies the falsifiability rule has both specificity (no
false emergence on isotropy) and sensitivity (real order clears the bar).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.emergence import synthetic as syn
from aleph.components.emergence.nematic import fiber_axes, nematic_order
from aleph.components.emergence.null_model import (
    EFFECT_SIZE_Z_CRIT,
    effect_size,
    is_ordered,
    nematic_null_band,
    nematic_null_scale_sq,
)


@pytest.mark.parametrize("n", [50, 100, 400, 1000])
def test_mc_sumsq_matches_analytic_3_over_2N(n: int) -> None:
    """The MC mean of sum(lambda^2) matches the closed-form 3/(2N) (the primary analytic cross-check)."""
    band = nematic_null_band(n, n_mc=3000, seed=100 + n)
    assert band.mc_sumsq == pytest.approx(nematic_null_scale_sq(n), rel=0.05)
    assert band.analytic_sumsq == pytest.approx(3.0 / (2.0 * n), rel=1e-12)


def test_null_bias_scales_as_one_over_sqrt_N() -> None:
    """mu_null(N) * sqrt(N) is ~constant: the isotropic S bias shrinks like 1/sqrt(N)."""
    ns = np.array([100, 200, 400, 800])
    scaled = []
    for n in ns:
        band = nematic_null_band(int(n), n_mc=3000, seed=7)
        scaled.append(band.mean * np.sqrt(n))
    scaled = np.array(scaled)
    # constant to ~15% across an 8x range in N
    assert scaled.std() / scaled.mean() < 0.15


def test_null_band_positive_bias_and_std() -> None:
    """The isotropic null has a strictly positive mean (finite-N bias) and a positive spread."""
    band = nematic_null_band(300, n_mc=3000, seed=11)
    assert band.mean > 0.0
    assert band.std > 0.0


def test_specificity_isotropic_not_ordered() -> None:
    """Independent isotropic draws score z < Z_CRIT: no false emergence (specificity)."""
    for seed in range(8):
        pos, off = syn.isotropic_network(400, seed=1000 + seed)
        s = nematic_order(fiber_axes(pos, off))
        band = nematic_null_band(400, n_mc=2000, seed=13)
        assert not is_ordered(s, band)


def test_sensitivity_aligned_far_above_bar() -> None:
    """A globally aligned network clears the bar by a mile (sensitivity): z >> Z_CRIT."""
    pos, off = syn.aligned_network(400, seed=5, jitter_deg=8.0)
    s = nematic_order(fiber_axes(pos, off))
    band = nematic_null_band(400, n_mc=2000, seed=13)
    assert is_ordered(s, band)
    assert effect_size(s, band) > 10.0 * EFFECT_SIZE_Z_CRIT


def test_effect_size_monotone_in_aligned_fraction() -> None:
    """The effect size increases monotonically with the aligned fraction (dose-response)."""
    band = nematic_null_band(2000, n_mc=1500, seed=21)
    z = []
    for p in np.linspace(0.0, 0.6, 7):
        s = nematic_order(fiber_axes(*syn.partially_aligned_network(2000, aligned_frac=p, seed=22)))
        z.append(effect_size(s, band))
    assert np.all(np.diff(z) > 0.0)
    assert z[0] < EFFECT_SIZE_Z_CRIT < z[-1]  # crosses the pre-registered bar as order rises


def test_z_crit_is_five_sigma_convention() -> None:
    """The pre-registered threshold is the 5-sigma discovery convention (documented, not a physics knob)."""
    assert EFFECT_SIZE_Z_CRIT == 5.0


def test_band_is_reproducible_from_seed() -> None:
    """The null band is fully determined by (N, n_mc, seed) — reproducible acceptance."""
    a = nematic_null_band(250, n_mc=1000, seed=99)
    b = nematic_null_band(250, n_mc=1000, seed=99)
    assert a.mean == b.mean and a.std == b.std and a.mc_sumsq == b.mc_sumsq
