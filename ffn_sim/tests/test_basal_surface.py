"""S1 — FLAT basal SURFACE layer (KU-3.5, PI directive (a) 2026-06-09).

The S-layer positions FA particles ON a FLAT 2D triangulated surface mesh (the
adherent ventral surface) and carries patch connectivity — NO force. Checks: FA on
the flat surface (z = z_basal), inside the footprint disk, triangulation connected,
and the resolution-invariance master gate (disk area independent of ring resolution).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.archive.hoomd_legacy.cell.basal_surface import (
    FlatBasalSurface,
    basal_surface_report,
    build_flat_basal_surface,
    disk_area_resolution_invariance,
)

FOOT_R = 5.0e-6
Z_BASAL = -7.0e-6


def _surf(n_fa=60, n_rings=10, seed=1):
    return build_flat_basal_surface(
        footprint_radius=FOOT_R, z_basal=Z_BASAL, n_rings=n_rings, n_fa=n_fa,
        rng=np.random.default_rng(seed),
    )


# --- §5 flat + within disk + §3 connected ---
def test_surface_is_planar_at_z_basal():
    s = _surf()
    assert np.allclose(s.verts[:, 2], Z_BASAL, atol=1e-15)
    assert np.allclose(s.fa_positions[:, 2], Z_BASAL, atol=1e-15)


def test_fa_within_footprint_disk():
    s = _surf()
    radial = np.linalg.norm(s.fa_positions[:, :2], axis=1)
    assert np.all(radial <= FOOT_R + 1e-12)


def test_triangulation_connected():
    s = _surf()
    rep = basal_surface_report(s)
    assert rep["controls"]["connected"]


def test_normal_is_plus_z():
    s = _surf()
    assert np.allclose(s.normal, [0.0, 0.0, 1.0])


def test_report_pass():
    s = _surf()
    rep = basal_surface_report(s)
    c = rep["controls"]
    assert c["planar"] and c["within_disk"] and c["connected"] and c["area_matches"]
    assert c["n_fa_placed"] == 60
    assert rep["verdict"] == "PASS"


def test_area_approaches_disk():
    s = _surf(n_rings=16)
    disk = math.pi * FOOT_R ** 2
    assert 0.9 <= s.area / disk <= 1.0 + 1e-9


# --- §2 boundary ---
@pytest.mark.parametrize("kw", [
    dict(footprint_radius=0.0),
    dict(footprint_radius=-1.0),
    dict(n_fa=0),
])
def test_boundary_raises(kw):
    base = dict(footprint_radius=FOOT_R, z_basal=Z_BASAL, n_fa=10, n_rings=8)
    base.update(kw)
    with pytest.raises(ValueError):
        build_flat_basal_surface(**base)


def test_n_fa_clamped_to_triangle_count():
    s = build_flat_basal_surface(
        footprint_radius=FOOT_R, z_basal=Z_BASAL, n_rings=4, n_fa=100_000,
        rng=np.random.default_rng(3),
    )
    assert s.fa_positions.shape[0] == s.n_tri
    assert s.fa_positions.shape[0] < 100_000


# --- §6 resolution-invariance master gate (VG-1) ---
def test_disk_area_resolution_invariance():
    res = disk_area_resolution_invariance(footprint_radius=FOOT_R, z_basal=Z_BASAL)
    assert res["invariant"], res
    # finest resolution within a few % of the analytic disk area.
    f16 = res["fraction_by_rings"][16]
    assert 0.93 <= f16 <= 1.0 + 1e-9


# --- bright line: the surface owns no force/bond/tension ---
def test_surface_is_geometry_only():
    s = _surf()
    for forbidden in ("k_edge", "edge_springs", "tension", "area_modulus", "params"):
        assert not hasattr(s, forbidden)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
