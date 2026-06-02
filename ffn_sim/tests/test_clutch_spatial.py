"""Static tests for bridge/clutch_spatial.py (β1-distribution observable).

Synthetic point sets verify the geometry metrics WITHOUT a simulation (pure numpy,
fast, no HOOMD). Import-isolated — cannot perturb the live build. Regression-guards the
observable against docs/BETA1_DISTRIBUTION_OBSERVABLE.md.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.bridge.clutch_spatial import beta1_distribution_metrics

RNG = np.random.default_rng(20260602)


def _uniform_disk(n: int, R: float = 1.0) -> np.ndarray:
    """n points uniform in area over a disk of radius R, centered at origin."""
    u = RNG.random(n)
    theta = RNG.uniform(-np.pi, np.pi, n)
    r = R * np.sqrt(u)
    return np.column_stack([r * np.cos(theta), r * np.sin(theta)])


def _to3d(xy: np.ndarray) -> np.ndarray:
    return np.column_stack([xy, np.zeros(len(xy))])


def test_uniform_disk_f_edge_near_reference():
    """A uniform areal disk → f_edge ≈ 1 − ρ² (0.51 at ρ=0.7) and low angular CV."""
    xy = _uniform_disk(5000)
    m = beta1_distribution_metrics(xy, np.ones(len(xy), bool))
    assert m["f_edge_uniform_ref"] == pytest.approx(0.51, abs=1e-9)
    assert m["f_edge"] == pytest.approx(0.51, abs=0.08)
    assert m["angular_cv"] < 0.15  # uniform angles
    assert m["n_engaged"] == 5000.0


def test_peripheral_ring_high_f_edge():
    """An edge ring (Pre phenotype) → f_edge ≈ 1, much above the uniform reference."""
    theta = RNG.uniform(-np.pi, np.pi, 2000)
    r = 1.0 + RNG.normal(0, 0.01, 2000)  # thin ring at r≈1
    xy = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    m = beta1_distribution_metrics(xy, np.ones(len(xy), bool))
    assert m["f_edge"] > 0.9
    assert m["f_edge"] > m["f_edge_uniform_ref"]
    assert m["angular_cv"] < 0.2  # ring is angularly uniform


def test_central_blob_low_f_edge():
    """A center-concentrated pattern → f_edge below the uniform reference."""
    u = RNG.random(4000)
    theta = RNG.uniform(-np.pi, np.pi, 4000)
    r = u ** 2  # concentrated near 0, sparse outliers set r_max
    xy = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    m = beta1_distribution_metrics(xy, np.ones(len(xy), bool))
    assert m["f_edge"] < m["f_edge_uniform_ref"]


def test_angular_cluster_high_cv():
    """Points confined to one angular wedge (clumped) → high angular CV."""
    theta = RNG.uniform(-0.1, 0.1, 1500)  # one narrow wedge
    r = np.sqrt(RNG.random(1500))
    xy = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    m_clump = beta1_distribution_metrics(xy, np.ones(len(xy), bool))
    m_unif = beta1_distribution_metrics(_uniform_disk(1500), np.ones(1500, bool))
    assert m_clump["angular_cv"] > 1.0
    assert m_clump["angular_cv"] > m_unif["angular_cv"]


def test_ordering_peripheral_uniform_central():
    """The discriminating ordering the observable exists to detect."""
    ring_theta = RNG.uniform(-np.pi, np.pi, 2000)
    ring = np.column_stack([np.cos(ring_theta), np.sin(ring_theta)])
    disk = _uniform_disk(2000)
    u = RNG.random(2000)
    ct = RNG.uniform(-np.pi, np.pi, 2000)
    central = np.column_stack([(u**2) * np.cos(ct), (u**2) * np.sin(ct)])
    f_ring = beta1_distribution_metrics(ring, np.ones(2000, bool))["f_edge"]
    f_disk = beta1_distribution_metrics(disk, np.ones(2000, bool))["f_edge"]
    f_cen = beta1_distribution_metrics(central, np.ones(2000, bool))["f_edge"]
    assert f_ring > f_disk > f_cen


def test_engaged_mask_subselects():
    """Only engaged points count; n_engaged reflects the mask."""
    xy = _uniform_disk(1000)
    mask = np.zeros(1000, bool)
    mask[:300] = True
    m = beta1_distribution_metrics(xy, mask)
    assert m["n_engaged"] == 300.0


def test_3d_positions_use_xy():
    """3D input uses the xy (en-face) projection."""
    xy = _uniform_disk(2000)
    m2 = beta1_distribution_metrics(xy, np.ones(2000, bool))
    m3 = beta1_distribution_metrics(_to3d(xy), np.ones(2000, bool))
    assert m3["f_edge"] == pytest.approx(m2["f_edge"])


def test_too_few_points_returns_nan():
    """Fewer than the minimum engaged points → NaN metrics, count preserved."""
    xy = np.array([[0.0, 0.0], [1.0, 1.0]])
    m = beta1_distribution_metrics(xy, np.ones(2, bool))
    assert m["n_engaged"] == 2.0
    assert math.isnan(m["f_edge"])
    assert math.isnan(m["angular_cv"])


def test_empty_engaged():
    """No engaged points → n_engaged 0, NaN metrics, no crash."""
    xy = _uniform_disk(100)
    m = beta1_distribution_metrics(xy, np.zeros(100, bool))
    assert m["n_engaged"] == 0.0
    assert math.isnan(m["f_edge"])


def test_invalid_args_raise():
    xy = _uniform_disk(10)
    with pytest.raises(ValueError):
        beta1_distribution_metrics(xy, np.ones(9, bool))  # mask misaligned
    with pytest.raises(ValueError):
        beta1_distribution_metrics(xy, np.ones(10, bool), edge_rho=1.5)
    with pytest.raises(ValueError):
        beta1_distribution_metrics(xy, np.ones(10, bool), n_sectors=1)
