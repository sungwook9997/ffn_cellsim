"""Sanity gates for the C-1 scorer — checked against curves whose answer is known in closed form.

Host arithmetic only (no Warp, no device): this is the half of C-1 that CAN be checked without a GPU,
and checking it here is what makes a GPU number afterwards interpretable. If the scorer is wrong, a
correct Biot solve still produces a wrong exponent and nothing downstream can tell.

The controls matter as much as the fits: a curve that never decays must not return its last sample as
a relaxation time, and a two-point sweep must not report r2 = 1.0 as though it had evidence.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aleph.scripts.c1_psr_scoring import (
    ONE_OVER_E,
    Q_L2,
    Q_PEAK,
    gaussian_spreading_tau,
    one_over_e_time,
    power_law_exponent,
)


# ── the multi-method result, in closed form, before any cell is simulated ─────────────────────────────
def test_two_observers_on_one_field_disagree_by_a_factor_of_three() -> None:
    """The C-1 thesis in miniature: same physics, same exponent, magnitudes 2.93x apart.

    An observer integrating the whole field (L2, q=3/4) and one reading a single point (peak, q=3/2)
    watch the SAME diffusing blob. Their 1/e times differ by (e^(4/3)-1)/(e^(2/3)-1), and there is no
    experiment that can call either wrong — they measure different functionals of one state.
    """
    a, d = 1.5, 50.0
    t_l2 = gaussian_spreading_tau(a, d, Q_L2)
    t_peak = gaussian_spreading_tau(a, d, Q_PEAK)
    ratio = (math.exp(4.0 / 3.0) - 1.0) / (math.exp(2.0 / 3.0) - 1.0)
    assert t_l2 / t_peak == pytest.approx(ratio, rel=1e-12)
    assert ratio == pytest.approx(2.947734, rel=1e-6)


def test_the_exponent_transfers_across_observers_even_though_the_magnitude_does_not() -> None:
    """Why C-1 scores an exponent: p = 2 for every observer, while the prefactor is observer-specific."""
    d = 50.0
    a = np.array([0.5, 1.0, 2.0, 4.0])
    prefactors = []
    for q in (Q_L2, Q_PEAK, 1.0, 3.0):
        p, c, r2 = power_law_exponent(a, np.array([gaussian_spreading_tau(v, d, q) for v in a]))
        assert p == pytest.approx(2.0, abs=1e-12)
        assert r2 == pytest.approx(1.0, abs=1e-12)
        prefactors.append(c)
    assert max(prefactors) / min(prefactors) > 3.0, "the prefactors must genuinely disagree"


def test_a_non_decaying_observer_is_refused() -> None:
    with pytest.raises(ValueError, match="q must be positive"):
        gaussian_spreading_tau(1.0, 50.0, 0.0)


# ── the 1/e crossing ─────────────────────────────────────────────────────────────────────────────────
def test_single_exponential_recovers_its_own_tau_exactly() -> None:
    """log-linear interpolation is EXACT for a single exponential — no tolerance is being spent here."""
    tau = 0.37
    t = np.linspace(0.0, 3.0 * tau, 41)
    assert one_over_e_time(t, np.exp(-t / tau)) == pytest.approx(tau, rel=1e-12)


def test_the_gaussian_spreading_law_is_recovered_from_its_own_curve() -> None:
    """The oracle and the crossing agree on the curve the oracle describes.

    ``||p||_2 / ||p_0||_2 = (a^2 / (a^2 + 2Dt))^(3/4)`` is a power law in t, not an exponential, which is
    exactly why the scorer reads a crossing instead of fitting a rate.
    """
    a, d = 1.7, 43.0
    tau = gaussian_spreading_tau(a, d)
    t = np.linspace(0.0, 4.0 * tau, 4001)
    amp = (a * a / (a * a + 2.0 * d * t)) ** 0.75
    assert one_over_e_time(t, amp) == pytest.approx(tau, rel=1e-4)


def test_a_curve_that_never_decays_returns_none_not_its_last_sample() -> None:
    """An un-decayed run is not a fast one. Returning the end time would read as exactly that."""
    t = np.linspace(0.0, 1.0, 20)
    assert one_over_e_time(t, np.full_like(t, 5.0)) is None
    assert one_over_e_time(t, 5.0 * np.exp(-t / 1e6)) is None


def test_a_rising_amplitude_is_refused() -> None:
    """A relaxation that gains amplitude is a defect (CFL violation, or a source left on)."""
    t = np.linspace(0.0, 1.0, 20)
    with pytest.raises(ValueError, match="rises during the decay"):
        one_over_e_time(t, np.exp(+t))


def test_non_monotone_time_is_refused() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        one_over_e_time([0.0, 0.2, 0.1], [1.0, 0.5, 0.2])


def test_the_crossing_is_at_the_declared_level() -> None:
    """Guard against the 1/e vs 1/2 mix-up the constant exists to prevent."""
    t = np.array([0.0, 1.0])
    assert one_over_e_time(t, np.array([1.0, ONE_OVER_E])) == pytest.approx(1.0, rel=1e-12)


# ── the exponent ─────────────────────────────────────────────────────────────────────────────────────
def test_a_clean_power_law_returns_its_exponent_and_prefactor() -> None:
    x = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
    p, c, r2 = power_law_exponent(x, 3.25 * x ** 2.0)
    assert p == pytest.approx(2.0, abs=1e-12)
    assert c == pytest.approx(3.25, rel=1e-12)
    assert r2 == pytest.approx(1.0, abs=1e-12)


def test_the_diffusive_prediction_is_what_the_oracle_produces() -> None:
    """tau ∝ a^2 at fixed D, and tau ∝ D^-1 at fixed a — the two exponents C-1 declares in advance."""
    d = 50.0
    a = np.array([0.5, 1.0, 2.0, 4.0])
    p_a, _, r2_a = power_law_exponent(a, np.array([gaussian_spreading_tau(v, d) for v in a]))
    assert p_a == pytest.approx(2.0, abs=1e-12) and r2_a == pytest.approx(1.0, abs=1e-12)

    ds = np.array([10.0, 25.0, 50.0, 100.0])
    p_d, _, _ = power_law_exponent(ds, np.array([gaussian_spreading_tau(1.0, v) for v in ds]))
    assert p_d == pytest.approx(-1.0, abs=1e-12)


def test_a_two_point_sweep_is_refused_rather_than_reporting_r2_one() -> None:
    """Two points always fit a line exactly. Reporting r2 = 1.0 there is evidence-shaped and empty."""
    with pytest.raises(ValueError, match="at least three points"):
        power_law_exponent([1.0, 2.0], [1.0, 4.0])


def test_r2_falls_when_the_points_are_not_a_power_law() -> None:
    """The fit-quality floor can only work if r2 actually moves; this is the control for that."""
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    _, _, r2_line = power_law_exponent(x, x ** 1.5)
    _, _, r2_bent = power_law_exponent(x, np.array([1.0, 9.0, 2.0, 30.0, 3.0]))
    assert r2_line == pytest.approx(1.0, abs=1e-12)
    assert r2_bent < 0.5


def test_non_positive_values_are_refused() -> None:
    with pytest.raises(ValueError, match="strictly positive"):
        power_law_exponent([1.0, 2.0, 3.0], [1.0, 0.0, 3.0])


def test_the_oracle_refuses_a_zero_length_scale() -> None:
    with pytest.raises(ValueError, match="must be positive"):
        gaussian_spreading_tau(0.0, 50.0)


def test_the_oracle_prefactor_is_the_declared_number() -> None:
    """C = (e^(4/3) - 1)/2 ~= 1.39683 — the prefactor C-1 declares BEFORE the run, not after."""
    assert gaussian_spreading_tau(1.0, 1.0) == pytest.approx((math.exp(4.0 / 3.0) - 1.0) / 2.0, rel=1e-15)
    assert gaussian_spreading_tau(1.0, 1.0) == pytest.approx(1.396834, rel=1e-6)
