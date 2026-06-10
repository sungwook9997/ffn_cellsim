"""Tests for the Slater radial stress(r) estimator (:mod:`ffn_sim.cortex.radial_stress`).

These are pure-numpy unit tests on the geometric primitive — no HOOMD build — so
they run fast and anywhere. The decisive check mirrors the estimator-audit
discipline (``scripts/h3_ku35_estimator_audit.py``): validate on a SYNTHETIC
field whose answer is known analytically.

Synthetic field: ``M`` isotropic radial rays, each a chain of equal-tension
``T₀`` segments pointing radially outward (force-balanced — every segment carries
the same axial force). At any shell radius strictly between two bead radii,
exactly one segment per ray crosses, so the total crossing radial force is
``M·T₀`` independent of ``r``. Hence:

* sphere:   ``σ(r)·4πr²   = M·T₀``  ⇒  ``σ ∝ r⁻²``  (Slater O1, 3D)
* cylinder: ``σ(r)·2πr·h  = M·T₀``  ⇒  ``σ ∝ r⁻¹``  (Slater O1, 2D-like)
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.cortex.radial_stress import (
    fit_power_law,
    radial_crossing_counts,
    radial_stress_profile,
)

T0 = 1.0e-12          # 1 pN per segment (cortex tension scale)
M_RAYS = 100          # isotropic radial rays
R_IN, R_OUT = 1.0e-6, 11.0e-6
N_BEADS = 21          # → 20 segments per ray


def _isotropic_dirs(m: int) -> np.ndarray:
    """``(m, 3)`` near-isotropic unit directions (Fibonacci sphere)."""
    i = np.arange(m, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / m)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.column_stack(
        [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]
    )


def _radial_chain_field(tension: float = T0):
    """Force-balanced radial-ray field. Returns ``(rA, rB, u, T, r_mid)``.

    ``r_mid`` are the segment-midpoint radii where each shell is crossed by
    exactly one segment per ray (so the analytic answer is exact).
    """
    dirs = _isotropic_dirs(M_RAYS)
    rb = np.linspace(R_IN, R_OUT, N_BEADS)
    rA_list, rB_list = [], []
    for d in dirs:
        for j in range(N_BEADS - 1):
            rA_list.append(rb[j] * d)
            rB_list.append(rb[j + 1] * d)
    rA = np.asarray(rA_list)
    rB = np.asarray(rB_list)
    d = rB - rA
    u = d / np.linalg.norm(d, axis=1)[:, None]
    T = np.full(rA.shape[0], tension)
    r_mid = 0.5 * (rb[:-1] + rb[1:])
    return rA, rB, u, T, r_mid


def test_sphere_force_balance_and_inverse_square():
    """σ(r)·4πr² == M·T₀ at every midpoint shell → exact r⁻² scaling."""
    rA, rB, u, T, r_mid = _radial_chain_field()
    sigma = radial_stress_profile(rA, rB, u, T, r_mid, geometry="sphere")
    flux = sigma * 4.0 * np.pi * r_mid ** 2
    assert np.allclose(flux, M_RAYS * T0, rtol=1e-9), flux / (M_RAYS * T0)
    fit = fit_power_law(r_mid, sigma)
    assert fit["exponent"] == pytest.approx(-2.0, abs=1e-6), fit
    assert fit["r2"] == pytest.approx(1.0, abs=1e-9), fit


def test_cylinder_force_balance_and_inverse_first_power():
    """σ(r)·2πr·h == M·T₀ → r⁻¹ scaling for the cylinder geometry."""
    rA, rB, u, T, r_mid = _radial_chain_field()
    h = 5.0e-6
    sigma = radial_stress_profile(
        rA, rB, u, T, r_mid, geometry="cylinder", height=h
    )
    flux = sigma * 2.0 * np.pi * r_mid * h
    assert np.allclose(flux, M_RAYS * T0, rtol=1e-9), flux / (M_RAYS * T0)
    fit = fit_power_law(r_mid, sigma)
    assert fit["exponent"] == pytest.approx(-1.0, abs=1e-6), fit


def test_sign_sense_tension_vs_compression():
    """T>0 → σ>0 (tensile); T<0 → σ<0 (compressive), same magnitude."""
    rA, rB, u, _, r_mid = _radial_chain_field()
    s_ten = radial_stress_profile(rA, rB, u, np.full(rA.shape[0], +T0), r_mid)
    s_com = radial_stress_profile(rA, rB, u, np.full(rA.shape[0], -T0), r_mid)
    assert np.all(s_ten > 0.0)
    assert np.all(s_com < 0.0)
    assert np.allclose(s_ten, -s_com, rtol=1e-12)


def test_ab_labeling_invariance():
    """Swapping a bond's A/B endpoints must not change σ (|û·r̂| removes it)."""
    rA, rB, u, T, r_mid = _radial_chain_field()
    s0 = radial_stress_profile(rA, rB, u, T, r_mid)
    # Flip every other bond's storage order (A↔B, û→−û).
    flip = np.zeros(rA.shape[0], dtype=bool)
    flip[::2] = True
    rA2, rB2, u2 = rA.copy(), rB.copy(), u.copy()
    rA2[flip], rB2[flip] = rB[flip], rA[flip]
    u2[flip] = -u[flip]
    s1 = radial_stress_profile(rA2, rB2, u2, T, r_mid)
    assert np.allclose(s0, s1, rtol=1e-12), np.max(np.abs(s0 - s1))


def test_no_bonds_and_no_crossing_return_zero():
    """Empty input and out-of-range shells give 0.0, never nan."""
    empty = radial_stress_profile(
        np.zeros((0, 3)), np.zeros((0, 3)), np.zeros((0, 3)),
        np.zeros(0), np.array([2.0e-6, 5.0e-6]),
    )
    assert empty.shape == (2,)
    assert np.all(empty == 0.0)

    rA, rB, u, T, _ = _radial_chain_field()
    # Shells well inside / outside the [R_IN, R_OUT] band have no crossing.
    outside = radial_stress_profile(
        rA, rB, u, T, np.array([0.1e-6, 50.0e-6])
    )
    assert np.all(outside == 0.0)


def test_crossing_counts_match_rays():
    """Each midpoint shell is crossed by exactly one segment per ray."""
    rA, rB, _, _, r_mid = _radial_chain_field()
    counts = radial_crossing_counts(rA, rB, r_mid)
    assert np.all(counts == M_RAYS), counts


def test_cylinder_requires_height():
    rA, rB, u, T, r_mid = _radial_chain_field()
    with pytest.raises(ValueError):
        radial_stress_profile(rA, rB, u, T, r_mid, geometry="cylinder")
    with pytest.raises(ValueError):
        radial_stress_profile(rA, rB, u, T, r_mid, geometry="bogus")


def test_fit_power_law_degenerate_returns_nan():
    """Fewer than two usable points → nan exponent, not a crash."""
    fit = fit_power_law(np.array([3.0e-6]), np.array([1.0]))
    assert math.isnan(fit["exponent"])
    assert fit["n_points"] == 1
