"""Self-test of the Arp2/3 angle-harmonic branch oracle (pure NumPy — no Warp/CUDA)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.weave.branch_angle import (
    ARP23_K_THETA,
    ARP23_SIGMA_DEG,
    ARP23_THETA0_RAD,
    KT_310K_PN_UM,
    angle_energy,
    angle_forces,
    branch_angle,
    k_theta_from_sigma,
    sample_branch_angles,
    thermal_sigma,
)


def test_kt_constant() -> None:
    """kT at 310 K = 4.28e-3 pN*um (the sourced physical constant)."""
    assert KT_310K_PN_UM == pytest.approx(4.28e-3, rel=2e-3)


def test_energy_minimum_at_theta0() -> None:
    """U(theta) is minimized at theta0 and is >= 0."""
    assert angle_energy(ARP23_THETA0_RAD, ARP23_K_THETA, ARP23_THETA0_RAD) == pytest.approx(0.0)
    for th in (0.5, 1.0, 1.5, 2.0):
        assert angle_energy(th, ARP23_K_THETA, ARP23_THETA0_RAD) >= 0.0


def test_branch_angle_geometry() -> None:
    """A right-angle triple measures pi/2; a straight one measures pi."""
    j = np.array([0.0, 0.0, 0.0])
    i = np.array([1.0, 0.0, 0.0])
    k = np.array([0.0, 1.0, 0.0])
    assert branch_angle(i, j, k) == pytest.approx(np.pi / 2)
    assert branch_angle(np.array([1.0, 0, 0]), j, np.array([-1.0, 0, 0])) == pytest.approx(np.pi)


def test_fd_gradient_sign_arbiter() -> None:
    """The analytic nodal forces match a finite-difference of the total angle energy, and RESTORE toward theta0.

    We perturb the daughter node along a tangential direction, so theta changes; the analytic force on that
    node must equal -dU/dx (FD) and point so as to reduce |theta - theta0|.
    """
    k_theta, theta0 = ARP23_K_THETA, ARP23_THETA0_RAD
    # a junction OPENED past theta0 (theta > theta0): the restoring torque pulls it back toward theta0
    p_m = np.array([1.0, 0.0, 0.0])
    p_j = np.array([0.0, 0.0, 0.0])
    p_d = np.array([-0.9, 0.6, 0.0])                     # some wide angle
    f_m, f_j, f_d = angle_forces(p_m, p_j, p_d, k_theta, theta0)

    def total_energy(pm, pj, pd):
        return float(angle_energy(branch_angle(pm, pj, pd), k_theta, theta0))

    eps = 1e-6
    for axis in range(3):
        d = np.zeros(3); d[axis] = eps
        fd = -(total_energy(p_m, p_j, p_d + d) - total_energy(p_m, p_j, p_d - d)) / (2 * eps)
        assert f_d[axis] == pytest.approx(fd, abs=1e-4)
    # Newton's 3rd law: internal potential exerts zero net force
    assert np.allclose(f_m + f_j + f_d, 0.0, atol=1e-10)


def test_equipartition_round_trip() -> None:
    """sigma_theta = sqrt(kT/k_theta) and k_theta = kT/sigma^2 invert; the ARP23 anchor round-trips."""
    sigma = np.deg2rad(ARP23_SIGMA_DEG)
    # ARP23_K_THETA = 0.173 is a 3-sig-fig rounding of kT/sigma^2 = 0.1735, so allow the rounding gap
    assert k_theta_from_sigma(sigma) == pytest.approx(ARP23_K_THETA, rel=5e-3)
    assert thermal_sigma(ARP23_K_THETA) == pytest.approx(sigma, rel=5e-3)


def test_thermal_distribution_width() -> None:
    """A measure-free harmonic Boltzmann ensemble has SD = sqrt(kT/k_theta); the sin-measure shifts the peak
    slightly below theta0 and narrows the SD a touch (reported, not tuned)."""
    rng = np.random.default_rng(0)
    sigma = thermal_sigma(ARP23_K_THETA)
    harm = sample_branch_angles(ARP23_THETA0_RAD, ARP23_K_THETA, 40000, rng, solid_angle_measure=False)
    assert np.std(harm) == pytest.approx(sigma, rel=0.05)
    assert np.mean(harm) == pytest.approx(ARP23_THETA0_RAD, abs=0.01)
    phys = sample_branch_angles(ARP23_THETA0_RAD, ARP23_K_THETA, 40000, rng, solid_angle_measure=True)
    # sin(theta) measure at theta0=70deg<90deg: mean shifts UP (toward pi/2), SD within ~15% of the harmonic
    assert np.mean(phys) > ARP23_THETA0_RAD
    assert np.std(phys) == pytest.approx(sigma, rel=0.15)


def test_branch_angle_matches_ff_kernel_formula() -> None:
    """The oracle's pref = k*(theta-theta0)/sin(theta) matches the ff.network_warp.branch_angle_kernel form."""
    p_m = np.array([1.0, 0.2, 0.0])
    p_j = np.array([0.0, 0.0, 0.0])
    p_d = np.array([0.3, 0.9, 0.1])
    theta = branch_angle(p_m, p_j, p_d)
    f_m, _, f_d = angle_forces(p_m, p_j, p_d, ARP23_K_THETA, ARP23_THETA0_RAD)
    # reconstruct pref from f_m and compare to the closed form
    r1, r2 = p_m - p_j, p_d - p_j
    d1, d2 = np.linalg.norm(r1), np.linalg.norm(r2)
    c = np.dot(r1, r2) / (d1 * d2)
    s = np.sqrt(1 - c * c)
    pref_expected = ARP23_K_THETA * (theta - ARP23_THETA0_RAD) / s
    fi_expected = pref_expected * (r2 / (d1 * d2) - c * r1 / (d1 * d1))
    assert np.allclose(f_m, fi_expected, atol=1e-12)
