"""Tests for the E-cadherin sliding-rebinding catch-bond oracle (L2.5, Rakshit 2012).

Validates the faithful sliding-rebinding implementation against the model's defining
properties and the measured catch peak (~f0 ≈ 29 pN). Pure NumPy, no HOOMD. This is the
closed-form side of the G5 emergent-vs-oracle gate.
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.validation import cadherin_sliding_rebinding as cb


def test_params_are_rakshit_table_s1():
    p = cb.RAKSHIT_W2A
    assert p.k_off0 == pytest.approx(30.4)
    assert p.x == pytest.approx(0.34e-9)
    assert p.k_on1 == pytest.approx(5.3)
    assert p.k_on2 == pytest.approx(1985.9)
    assert p.f0 == pytest.approx(29.2e-12)
    assert p.n == pytest.approx(4.8)


def test_new_interaction_probability_bounds():
    """Pn = 0 below 0, 1 above f0, monotone sin^{2n} ramp in between, =1 at f0."""
    assert cb.new_interaction_probability(-1e-12) == 0.0
    assert cb.new_interaction_probability(0.0) == pytest.approx(0.0)
    assert cb.new_interaction_probability(40e-12) == pytest.approx(1.0)
    assert cb.new_interaction_probability(29.2e-12) == pytest.approx(1.0, abs=1e-9)
    f = np.linspace(0.0, 29.2e-12, 50)
    pn = cb.new_interaction_probability(f)
    assert np.all(np.diff(pn) >= -1e-12)          # monotone non-decreasing
    assert np.all((pn >= 0.0) & (pn <= 1.0))


def test_single_pair_off_rate_is_bell_slip():
    """k-1(f) increases with force (Bell slip) and equals k_off0 at f=0."""
    assert cb.single_pair_off_rate(0.0) == pytest.approx(30.4)
    assert cb.single_pair_off_rate(20e-12) > cb.single_pair_off_rate(0.0)
    with pytest.raises(ValueError):
        cb.single_pair_off_rate(-1e-12)


def test_generator_columns_leak_nonnegative():
    """Probability conservation: every transient-state leak to the absorbing P00 is ≥ 0.

    This is the check that fixes the SI subscript ambiguity (slide term must be k-1, not k+1).
    """
    for f_pN in [0.0, 5.0, 10.0, 20.0, 29.2, 40.0, 60.0]:
        M = cb._generator(f_pN * 1e-12, cb.RAKSHIT_W2A)
        leak = -M.sum(axis=0)                      # rate each transient column sends to P00
        assert np.all(leak >= -1e-9), (f_pN, leak)


def test_mean_lifetime_positive_and_finite():
    for f_pN in [0.0, 10.0, 29.2, 50.0]:
        tau = cb.mean_lifetime(f_pN * 1e-12)
        assert np.isfinite(tau) and tau > 0.0


def test_catch_peak_near_f0():
    """The defining catch signature: lifetime peaks near f0 ≈ 29 pN (measured ~30 pN)."""
    f_star, tau_star = cb.catch_peak_force(cb.RAKSHIT_W2A)
    assert 24e-12 <= f_star <= 32e-12             # peak within a few pN of f0
    # biphasic: lifetime well above f0 is shorter than at the peak (slip)
    assert cb.mean_lifetime(60e-12) < tau_star
    # and the peak is a genuine catch rise above the mid-range minimum
    assert tau_star > cb.mean_lifetime(17e-12)


def test_effective_k_off_is_inverse_lifetime():
    f = 25e-12
    assert cb.effective_k_off(f) == pytest.approx(1.0 / cb.mean_lifetime(f), rel=1e-9)
