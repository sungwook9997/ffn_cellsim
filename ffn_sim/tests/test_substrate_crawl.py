"""Tests for the D-axis substrate-reacted collective in-plane crawl force."""
from __future__ import annotations
import numpy as np
import pytest
from ffn_sim.archive.hoomd_legacy.spheroid.substrate_crawl import substrate_crawl_forces


def test_zero_crawl_is_zero():
    p = np.random.default_rng(0).normal(size=(20, 3))
    f = substrate_crawl_forces(p, f_crawl=0.0, z_substrate=0.0, basal_band=1.0)
    assert np.allclose(f, 0.0)


def test_force_is_in_plane_zero_z():
    p = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    f = substrate_crawl_forces(p, f_crawl=5.0, z_substrate=0.0, basal_band=1.0)
    assert np.allclose(f[:, 2], 0.0)                       # crawl is tangent to the dish


def test_outward_direction_and_magnitude():
    # two basal cells on the x-axis straddling the centroid -> opposite outward x, |f|=f_crawl
    p = np.array([[3.0, 0.0, 0.0], [-3.0, 0.0, 0.0]])
    f = substrate_crawl_forces(p, f_crawl=9.4e-9, z_substrate=0.0, basal_band=1.0)
    assert f[0, 0] > 0 and f[1, 0] < 0                     # outward from centroid
    assert np.linalg.norm(f[0]) == pytest.approx(9.4e-9, rel=1e-9)


def test_apical_cells_dont_crawl():
    # two basal cells (offset from their centroid) crawl; an apical cell (z >> band) feels ~0
    p = np.array([[3.0, 0.0, 0.0], [-3.0, 0.0, 0.0], [3.0, 0.0, 10.0]])
    f = substrate_crawl_forces(p, f_crawl=1.0, z_substrate=0.0, basal_band=1.0)
    assert np.linalg.norm(f[2]) == pytest.approx(0.0, abs=1e-12)   # apical: no crawl
    assert np.linalg.norm(f[0]) > 0.0 and np.linalg.norm(f[1]) > 0.0  # basal: crawl


def test_basal_band_guard():
    with pytest.raises(ValueError):
        substrate_crawl_forces(np.zeros((3, 3)), f_crawl=1.0, z_substrate=0.0, basal_band=0.0)
