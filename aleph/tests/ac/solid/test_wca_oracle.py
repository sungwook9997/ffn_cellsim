"""Self-test of the WCA excluded-volume oracle (pure NumPy — runs anywhere, no Warp/CUDA).

Validates the analytic ground truth itself before it is used to gate the Warp-CUDA I2b steric kernel.
This is acceptance-layer math, not simulation, so it is Warp-only-contract exempt by construction.

Gates (per the Session-C brief): pair-potential FD-gradient sign arbiter; zero-overlap == zero-force;
compressed-pair resists collapse (finite steric volume); k_EV <-> epsilon bookkeeping.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.solid.wca_analytic import (
    compressed_pair_separation,
    contact_stiffness,
    epsilon_from_contact_stiffness,
    max_core_stiffness,
    steric_exclusion_volume,
    wca_cutoff,
    wca_energy,
    wca_force_magnitude,
    wca_local_stiffness,
    wca_pair_force,
)

# Representative um-pN-s scales for the tests (values do not gate physics — the laws hold for any >0).
SIGMA = 0.05          # um steric diameter (illustrative coarse-grained shell)
K_EV = 200.0          # pN/um contact stiffness (numerical repulsion scale)
EPS = epsilon_from_contact_stiffness(K_EV, SIGMA)


# --------------------------------------------------------------------------------------------------
# Gate 1 — pair-potential FD-gradient (the SIGN ARBITER)
# --------------------------------------------------------------------------------------------------
def test_force_is_negative_gradient_of_energy() -> None:
    """F(r) = -dU/dr to O(h^2) across the repulsive branch (central finite difference)."""
    r_c = wca_cutoff(SIGMA)
    r = np.linspace(0.75 * SIGMA, 0.999 * r_c, 400)
    h = 1e-8
    fd = -(wca_energy(r + h, SIGMA, EPS) - wca_energy(r - h, SIGMA, EPS)) / (2.0 * h)
    analytic = wca_force_magnitude(r, SIGMA, EPS)
    # relative agreement (forces span orders of magnitude across the branch)
    assert np.allclose(fd, analytic, rtol=1e-5, atol=1e-6 * np.max(np.abs(analytic)))


def test_force_is_repulsive_everywhere_inside_cutoff() -> None:
    """F > 0 (pushes apart) for all 0 < r < r_c; the potential is strictly decreasing there."""
    r_c = wca_cutoff(SIGMA)
    r = np.linspace(0.6 * SIGMA, 0.9999 * r_c, 500)
    f = wca_force_magnitude(r, SIGMA, EPS)
    assert np.all(f > 0.0)
    u = wca_energy(r, SIGMA, EPS)
    assert np.all(np.diff(u) < 0.0)          # monotone decreasing toward the cutoff


def test_pair_force_vector_points_apart() -> None:
    """The vector force on i points from j toward i (separating), with the scalar magnitude."""
    ri = np.array([0.02, 0.0, 0.0])
    rj = np.array([-0.02, 0.0, 0.0])          # separation 0.04 um < r_c
    f = wca_pair_force(ri, rj, SIGMA, EPS)
    assert f[0] > 0.0 and abs(f[1]) < 1e-15 and abs(f[2]) < 1e-15   # pushes i in +x, away from j
    mag = float(wca_force_magnitude(0.04, SIGMA, EPS))
    assert np.linalg.norm(f) == pytest.approx(mag, rel=1e-12)


# --------------------------------------------------------------------------------------------------
# Gate 2 — zero-overlap == zero-force (OFF / regression) + C1 continuity at the cutoff
# --------------------------------------------------------------------------------------------------
def test_zero_beyond_cutoff() -> None:
    """U and F are identically 0 at and beyond r_c (no attractive tail; OFF == bit-identical)."""
    r_c = wca_cutoff(SIGMA)
    r = np.array([r_c, r_c * 1.0001, 1.5 * r_c, 10.0 * r_c])
    assert np.all(wca_energy(r, SIGMA, EPS) == 0.0)
    assert np.all(wca_force_magnitude(r, SIGMA, EPS) == 0.0)
    assert np.all(wca_local_stiffness(r, SIGMA, EPS) == 0.0)


def test_c1_continuity_at_cutoff() -> None:
    """Energy AND force -> 0 continuously as r -> r_c from below (WCA is C1, unlike a hard shell)."""
    r_c = wca_cutoff(SIGMA)
    r = r_c * (1.0 - np.array([1e-3, 1e-4, 1e-5]))
    assert np.all(wca_energy(r, SIGMA, EPS) > 0.0)
    assert wca_energy(r, SIGMA, EPS)[-1] < wca_energy(r, SIGMA, EPS)[0]     # decreasing to 0
    f = wca_force_magnitude(r, SIGMA, EPS)
    assert np.all(f > 0.0) and f[-1] < f[0]                                 # force -> 0 too


def test_lj_minimum_at_cutoff() -> None:
    """The LJ minimum sits at r_min = r_c = 2^(1/6) sigma; U_WCA(r_c) = 0 with U_LJ(r_min) = -eps."""
    r_c = wca_cutoff(SIGMA)
    assert r_c == pytest.approx(2.0 ** (1.0 / 6.0) * SIGMA, rel=1e-14)
    # WCA energy at contact r=sigma equals eps (shifted LJ), and is 0 exactly at the cutoff
    assert float(wca_energy(SIGMA, SIGMA, EPS)) == pytest.approx(EPS, rel=1e-12)
    assert float(wca_energy(r_c * (1.0 - 1e-12), SIGMA, EPS)) == pytest.approx(0.0, abs=1e-6 * EPS)


# --------------------------------------------------------------------------------------------------
# Gate 3 — k_EV <-> epsilon bookkeeping (the CFL / Magic-Number-Block scalar layer)
# --------------------------------------------------------------------------------------------------
def test_stiffness_epsilon_roundtrip() -> None:
    """k_EV -> epsilon -> k_EV is exact; contact_stiffness is the well curvature U''(r_min)."""
    eps = epsilon_from_contact_stiffness(K_EV, SIGMA)
    assert contact_stiffness(eps, SIGMA) == pytest.approx(K_EV, rel=1e-12)
    r_c = wca_cutoff(SIGMA)
    k_num = float(wca_local_stiffness(r_c * (1.0 - 1e-7), SIGMA, eps))     # curvature just inside r_min
    assert k_num == pytest.approx(K_EV, rel=1e-4)


def test_core_stiffness_exceeds_contact_stiffness() -> None:
    """Under a load cap the CFL stiffness is the deeper-overlap U''(r_eq) >= the contact stiffness."""
    f_cap = 50.0     # pN
    k_core = max_core_stiffness(f_cap, SIGMA, EPS)
    assert k_core >= contact_stiffness(EPS, SIGMA)
    # a larger cap => deeper overlap => stiffer
    assert max_core_stiffness(500.0, SIGMA, EPS) > max_core_stiffness(50.0, SIGMA, EPS)


# --------------------------------------------------------------------------------------------------
# Gate 4 — a compressed pair resists collapse (FINITE steric volume)
# --------------------------------------------------------------------------------------------------
def test_compressed_pair_finite_separation() -> None:
    """For any finite load the pair settles at 0 < r_eq < r_c with F_WCA(r_eq) == load (force balance)."""
    r_c = wca_cutoff(SIGMA)
    for f_load in (1.0, 10.0, 100.0, 1.0e3, 1.0e4):
        r_eq = compressed_pair_separation(f_load, SIGMA, EPS)
        assert 0.0 < r_eq < r_c
        assert float(wca_force_magnitude(r_eq, SIGMA, EPS)) == pytest.approx(f_load, rel=1e-6)


def test_compressed_pair_monotone_and_never_collapses() -> None:
    """r_eq decreases with load but stays strictly > 0 (the steric core never fully overlaps)."""
    loads = np.array([1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6])
    r_eq = np.array([compressed_pair_separation(f, SIGMA, EPS) for f in loads])
    assert np.all(np.diff(r_eq) < 0.0)          # stiffer squeeze -> closer, monotonically
    assert np.all(r_eq > 0.0)                    # but never zero -> finite steric volume
    assert np.all(np.array([steric_exclusion_volume(r) for r in r_eq]) > 0.0)


def test_steric_volume_bounded_below_under_extreme_load() -> None:
    """Even a 1e8 pN load leaves a positive separation (the ~r^-13 core diverges faster than any load)."""
    r_eq = compressed_pair_separation(1.0e8, SIGMA, EPS)
    assert r_eq > 0.0
    assert steric_exclusion_volume(r_eq) > 0.0
