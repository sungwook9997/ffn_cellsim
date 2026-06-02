"""Static tests for the Layer-2 spheroid measurement + acceptance-oracle layer.

Synthetic point clouds with KNOWN closed-form answers verify the geometry/measurement
math WITHOUT a simulation (pure numpy/scipy, fast, no HOOMD). Import-isolated — cannot
perturb the live single-cell build. This is the "measurement protocol before the run"
discipline: the A/A0 sweep (G3) and the stable-aggregate gate (G1) both rest on these.

Cross-checks the runtime ``ffn_sim.spheroid.observables`` against the closed-form oracle
``ffn_sim.validation.oracles.spheroid.aa0_law`` (tests may import oracles; runtime may not).
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.spheroid import observables as obs
from ffn_sim.validation.oracles.spheroid import aa0_law as oracle

RNG = np.random.default_rng(20260602)


# --------------------------------------------------------------------------- #
# Synthetic generators
# --------------------------------------------------------------------------- #
def _uniform_disk_xy(n: int, R: float) -> np.ndarray:
    """n points uniform in area over a disk of radius R in the z=0 plane (N x 3)."""
    u = RNG.random(n)
    theta = RNG.uniform(-np.pi, np.pi, n)
    r = R * np.sqrt(u)
    return np.column_stack([r * np.cos(theta), r * np.sin(theta), np.zeros(n)])


def _uniform_ball(n: int, R: float) -> np.ndarray:
    """n points uniform in volume over a ball of radius R (N x 3)."""
    v = RNG.random(n)
    r = R * np.cbrt(v)
    phi = RNG.uniform(-np.pi, np.pi, n)
    cos_t = RNG.uniform(-1.0, 1.0, n)
    sin_t = np.sqrt(1.0 - cos_t**2)
    return np.column_stack([r * sin_t * np.cos(phi), r * sin_t * np.sin(phi), r * cos_t])


def _cubic_lattice(n_side: int, spacing: float) -> np.ndarray:
    """Simple-cubic lattice, n_side^3 points at the given spacing (N x 3)."""
    g = np.arange(n_side) * spacing
    xx, yy, zz = np.meshgrid(g, g, g, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])


# --------------------------------------------------------------------------- #
# Measurement vs analytic references
# --------------------------------------------------------------------------- #
def test_projected_area_matches_disk_pi_r2():
    """Points filling a disk of radius R -> convex-hull area -> pi R^2."""
    R = 7.5e-6  # MCF7 cell radius scale
    xy = _uniform_disk_xy(8000, R)
    area = obs.projected_area(xy)
    assert area == pytest.approx(oracle.uniform_disk_area(R), rel=0.03)


def test_effective_radius_inverts_area():
    R = 1.234e-5
    assert obs.effective_radius(oracle.uniform_disk_area(R)) == pytest.approx(R, rel=1e-12)


def test_radius_of_gyration_solid_ball():
    """Uniform ball radius R -> Rg -> sqrt(3/5) R."""
    R = 5.0e-5
    p = _uniform_ball(20000, R)
    assert obs.radius_of_gyration(p) == pytest.approx(
        oracle.solid_sphere_radius_of_gyration(R), rel=0.03
    )


def test_radius_of_gyration_translation_invariant():
    p = _uniform_ball(2000, 1.0e-5)
    rg0 = obs.radius_of_gyration(p)
    rg1 = obs.radius_of_gyration(p + np.array([3.0e-4, -1.0e-4, 7.0e-4]))
    assert rg1 == pytest.approx(rg0, rel=1e-12)


def test_nearest_neighbor_on_cubic_lattice():
    """Simple-cubic lattice spacing s -> every interior NN distance == s."""
    s = 1.5e-5
    p = _cubic_lattice(6, s)
    stats = obs.nearest_neighbor_stats(p)
    assert stats["median"] == pytest.approx(s, rel=1e-9)
    assert stats["min"] == pytest.approx(s, rel=1e-9)


def test_projected_area_degenerate_returns_zero():
    assert obs.projected_area(np.zeros((2, 3))) == 0.0
    collinear = np.column_stack([np.linspace(0, 1, 5), np.zeros(5), np.zeros(5)])
    assert obs.projected_area(collinear) == 0.0


def test_detached_fraction_blob_plus_outliers():
    """A compact blob + k far isolated points -> detached fraction = k / N."""
    blob = _uniform_ball(200, 1.0e-5)  # radius 10 um
    outliers = np.array([[1.0e-3, 0, 0], [0, 1.0e-3, 0], [0, 0, 1.0e-3]])  # 1 mm away
    p = np.vstack([blob, outliers])
    frac = obs.detached_fraction(p, d_crit=5.0e-5, neighbor_radius=5.0e-6)
    assert frac == pytest.approx(3.0 / len(p), abs=1e-9)


def test_radial_density_profile_runs_and_is_interior_flat():
    """Uniform ball -> interior shells roughly constant density, then drops at edge."""
    R = 5.0e-5
    p = _uniform_ball(30000, R)
    r_c, dens = obs.radial_density_profile(p, n_bins=20)
    assert r_c.shape == dens.shape == (20,)
    interior = dens[2:12]  # skip the noisy innermost + the edge falloff
    assert interior.std() / interior.mean() < 0.25


# --------------------------------------------------------------------------- #
# A/A0 spreading-law oracle
# --------------------------------------------------------------------------- #
def test_aa0_model_baseline_limit():
    """R -> large suppresses b/R, c/R^2 -> A/A0 -> a."""
    a, b, c = 1.0, 2.0e-5, -3.0e-10
    big = oracle.aa0_model(1.0e3, a, b, c)
    assert big == pytest.approx(a, abs=1e-6)


def test_aa0_fit_recovers_coefficients_noise_free():
    """Noise-free synthetic A/A0(R) -> fit recovers (a, b, c) exactly."""
    a, b, c = 0.8, 1.5e-5, -2.0e-10
    R = np.linspace(2.0e-5, 2.0e-4, 12)  # 20-200 um initial radii
    AA0 = oracle.aa0_model(R, a, b, c)
    fit = oracle.fit_aa0(R, AA0)
    assert fit["a"] == pytest.approx(a, rel=1e-6)
    assert fit["b"] == pytest.approx(b, rel=1e-6)
    assert fit["c"] == pytest.approx(c, rel=1e-6)
    assert fit["r_squared"] == pytest.approx(1.0, abs=1e-9)


def test_aa0_fit_requires_three_distinct_radii():
    with pytest.raises(ValueError):
        oracle.fit_aa0(np.array([1e-5, 1e-5, 1e-5]), np.array([1.0, 1.0, 1.0]))


def test_aa0_model_rejects_nonpositive_radius():
    with pytest.raises(ValueError):
        oracle.aa0_model(np.array([1e-5, 0.0]), 1.0, 1.0e-5, 0.0)


def test_term_contributions_dominance_crosses_over():
    """c/R^2 dominates at small R; a dominates at large R (edge-vs-bulk readout)."""
    a, b, c = 1.0, 1.0e-5, 1.0e-9
    small = oracle.term_contributions(1.0e-6, a, b, c)
    large = oracle.term_contributions(1.0e-3, a, b, c)
    assert small["dominant"] == "c"
    assert large["dominant"] == "a"


# --------------------------------------------------------------------------- #
# Connected-component / fragmentation-robust core area (L2.4 spread observable)
# --------------------------------------------------------------------------- #
def test_connected_components_two_separated_clusters():
    """Two clusters far apart → 2 components; label 0 is the larger one."""
    r0 = 1.0
    a = _uniform_disk_xy(40, 2.0 * r0)                       # bigger cluster at origin
    b = _uniform_disk_xy(15, 1.0 * r0) + np.array([50.0, 0.0, 0.0])  # smaller, far away
    pos = np.vstack([a, b])
    labels = obs.connected_components(pos, link_radius=1.6 * r0)
    assert labels.max() == 1                                 # exactly two components
    assert (labels == 0).sum() == 40                         # label 0 = largest


def test_core_area_ignores_drifting_fragment():
    """Core area = the big cluster's hull, NOT the whole-set hull spanning the gap."""
    r0 = 1.0
    big = _uniform_disk_xy(60, 3.0 * r0)
    stray = _uniform_disk_xy(8, 0.5 * r0) + np.array([80.0, 0.0, 0.0])
    pos = np.vstack([big, stray])
    a_core = obs.core_projected_area(pos, link_radius=1.6 * r0)
    a_hull = obs.projected_area(pos)
    assert a_core < 0.3 * a_hull                             # hull is inflated by the gap
    assert a_core == pytest.approx(obs.projected_area(big), rel=1e-9)


def test_connected_components_single_cluster_is_one_label():
    pos = _uniform_disk_xy(50, 2.0)
    labels = obs.connected_components(pos, link_radius=3.0)  # generous link → all one cluster
    assert labels.max() == 0
    assert obs.largest_connected_component(pos, link_radius=3.0).all()


def test_connected_components_rejects_nonpositive_link():
    with pytest.raises(ValueError):
        obs.connected_components(_uniform_disk_xy(5, 1.0), link_radius=0.0)
