"""Self-consistency + known-physics tests for the forward Hertz S1 oracle.

Validates ``aleph.validation.oracles.mechanics.hertz`` against its own analytic
identities (roundtrip, stiffness = dF/ddelta, pressure ratios, between-plates
factor) and closed-form values, so the oracle can be trusted as an acceptance
gate for the S1 force–displacement run. Pure-Python, no GPU.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.validation.oracles.mechanics import hertz

# Representative MCF7-scale inputs (SI): E ~ 1 kPa cell, R = 7.5 um, rigid platen.
E_CELL = 1.0e3          # Pa
NU = 0.45               # ~incompressible cytoplasm
R = 7.5e-6              # m
E_STAR = hertz.rigid_indenter_reduced_modulus(E_CELL, NU)


# --------------------------------------------------------------------------- #
# Reduced modulus
# --------------------------------------------------------------------------- #
def test_rigid_indenter_reduced_modulus_matches_formula():
    assert E_STAR == pytest.approx(E_CELL / (1.0 - NU**2))


def test_youngs_from_reduced_is_inverse():
    assert hertz.youngs_from_reduced(E_STAR, NU) == pytest.approx(E_CELL)


def test_reduced_modulus_two_body_rigid_limit():
    # E2 -> inf must recover the rigid-indenter single-body result.
    e_star_two = hertz.reduced_modulus(E_CELL, NU, float("inf"), 0.5)
    assert e_star_two == pytest.approx(E_STAR)


def test_reduced_modulus_symmetric_and_softer_than_either():
    e_star = hertz.reduced_modulus(1.0e3, 0.45, 2.0e3, 0.3)
    # Compliances add in series -> E* below the stiffer single-body reduced modulus.
    assert e_star < hertz.rigid_indenter_reduced_modulus(2.0e3, 0.3)
    assert e_star == pytest.approx(
        hertz.reduced_modulus(2.0e3, 0.3, 1.0e3, 0.45)
    )


# --------------------------------------------------------------------------- #
# Forward law: boundary, monotonicity, power law
# --------------------------------------------------------------------------- #
def test_zero_indentation_zero_force_and_stiffness():
    assert hertz.force_sphere_plane(E_STAR, R, 0.0) == pytest.approx(0.0)
    assert hertz.contact_stiffness(E_STAR, R, 0.0) == pytest.approx(0.0)


def test_force_is_monotone_increasing():
    d = np.linspace(0.0, 0.02 * R, 50)
    F = hertz.force_sphere_plane(E_STAR, R, d)
    assert np.all(np.diff(F) >= 0.0)


def test_three_halves_power_law():
    # Doubling delta multiplies F by 2^1.5.
    f1 = hertz.force_sphere_plane(E_STAR, R, 0.01 * R)
    f2 = hertz.force_sphere_plane(E_STAR, R, 0.02 * R)
    assert f2 / f1 == pytest.approx(2.0**1.5)


# --------------------------------------------------------------------------- #
# Analytic identities
# --------------------------------------------------------------------------- #
def test_delta_from_force_roundtrip():
    d = np.linspace(1e-9, 0.02 * R, 25)
    F = hertz.force_sphere_plane(E_STAR, R, d)
    d_back = hertz.delta_from_force(E_STAR, R, F)
    np.testing.assert_allclose(d_back, d, rtol=1e-12, atol=0.0)


def test_stiffness_equals_numerical_derivative():
    d0 = 0.01 * R
    h = 1e-4 * R
    num = (
        hertz.force_sphere_plane(E_STAR, R, d0 + h)
        - hertz.force_sphere_plane(E_STAR, R, d0 - h)
    ) / (2.0 * h)
    assert hertz.contact_stiffness(E_STAR, R, d0) == pytest.approx(float(num), rel=1e-5)


def test_stiffness_equals_two_estar_times_contact_radius():
    d = 0.015 * R
    a = hertz.contact_radius(R, d)
    assert hertz.contact_stiffness(E_STAR, R, d) == pytest.approx(2.0 * E_STAR * a)


def test_contact_radius_formula():
    d = 0.01 * R
    assert hertz.contact_radius(R, d) == pytest.approx(np.sqrt(R * d))


def test_max_pressure_equals_3F_over_2pi_a2():
    d = 0.02 * R
    F = float(hertz.force_sphere_plane(E_STAR, R, d))
    a = float(hertz.contact_radius(R, d))
    p0_geom = 3.0 * F / (2.0 * np.pi * a * a)
    assert hertz.max_contact_pressure(E_STAR, R, d) == pytest.approx(p0_geom)


def test_mean_pressure_is_two_thirds_peak():
    d = 0.02 * R
    p0 = float(hertz.max_contact_pressure(E_STAR, R, d))
    assert hertz.mean_contact_pressure(E_STAR, R, d) == pytest.approx((2.0 / 3.0) * p0)


# --------------------------------------------------------------------------- #
# Sphere between plates
# --------------------------------------------------------------------------- #
def test_between_plates_equals_single_contact_at_half_approach():
    D = 0.02 * R
    assert hertz.force_between_plates(E_STAR, R, D) == pytest.approx(
        float(hertz.force_sphere_plane(E_STAR, R, D / 2.0))
    )


def test_between_plates_explicit_sqrt2_over_3_prefactor():
    D = 0.02 * R
    expected = (np.sqrt(2.0) / 3.0) * E_STAR * np.sqrt(R) * D**1.5
    assert hertz.force_between_plates(E_STAR, R, D) == pytest.approx(expected)


def test_between_plates_softer_than_single_contact_at_equal_total_approach():
    D = 0.02 * R
    assert hertz.force_between_plates(E_STAR, R, D) < hertz.force_sphere_plane(E_STAR, R, D)


# --------------------------------------------------------------------------- #
# Inverse / fit
# --------------------------------------------------------------------------- #
def test_fit_recovers_reduced_modulus_from_ideal_curve():
    d = np.linspace(1e-9, 0.03 * R, 40)
    F = hertz.force_sphere_plane(E_STAR, R, d)
    fit = hertz.fit_reduced_modulus(d, F, R)
    assert fit["E_star"] == pytest.approx(E_STAR, rel=1e-9)
    assert fit["r_squared"] == pytest.approx(1.0, abs=1e-12)
    # Ideal Hertzian curve -> pointwise E* is flat (CV ~ 0).
    assert fit["E_star_flatness"] == pytest.approx(0.0, abs=1e-9)


def test_pointwise_inversion_matches_forward():
    d = 0.02 * R
    F = float(hertz.force_sphere_plane(E_STAR, R, d))
    assert hertz.reduced_modulus_from_force(R, d, F) == pytest.approx(E_STAR, rel=1e-12)


def test_fit_flags_nonhertzian_curve_via_flatness():
    # A curve that is NOT delta^1.5 (linear) -> non-flat pointwise E*.
    d = np.linspace(1e-9, 0.03 * R, 40)
    F = 1e-6 * (d / R)            # linear-in-delta fake response
    fit = hertz.fit_reduced_modulus(d, F, R)
    assert fit["E_star_flatness"] > 0.1


# --------------------------------------------------------------------------- #
# Validity window
# --------------------------------------------------------------------------- #
def test_small_strain_window():
    d = np.array([0.01, 0.03, 0.05]) * R
    v = hertz.valid_small_strain(d, R, max_ratio=0.03)
    assert list(v["valid"]) == [True, True, False]
    assert v["all_valid"] is False
    np.testing.assert_allclose(v["ratio"], [0.01, 0.03, 0.05])


# --------------------------------------------------------------------------- #
# Guards
# --------------------------------------------------------------------------- #
def test_negative_delta_raises():
    with pytest.raises(ValueError):
        hertz.force_sphere_plane(E_STAR, R, -1e-9)


def test_nonpositive_modulus_and_radius_raise():
    with pytest.raises(ValueError):
        hertz.force_sphere_plane(0.0, R, 1e-9)
    with pytest.raises(ValueError):
        hertz.force_sphere_plane(E_STAR, 0.0, 1e-9)


def test_reduced_modulus_rejects_incompressible_nu_one():
    with pytest.raises(ValueError):
        hertz.rigid_indenter_reduced_modulus(E_CELL, 1.0)
