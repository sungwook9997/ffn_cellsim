"""S1 — basal SURFACE layer (KU-3.5, PI directive (a) 2026-06-09).

The S-layer positions FA particles ON the 2D surface mesh (surface_manifold) and
carries patch connectivity — NO force. Checks: FA on the surface (radial=R), inside
the south-pole basal cap, basal patches connected, and the resolution-invariance
master gate (cap area fraction independent of subdivision level).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.cell.basal_surface import (
    BasalSurface,
    build_basal_surface,
    cap_fraction_resolution_invariance,
    cap_half_angle,
    select_basal_cap,
    basal_surface_report,
)
from ffn_sim.cortex.surface_manifold import SurfaceManifold

R = 7.5e-6
FOOT_R = 5.0e-6


def test_cap_half_angle_formula():
    assert math.isclose(cap_half_angle(FOOT_R, R), math.asin(FOOT_R / R), rel_tol=1e-12)


@pytest.mark.parametrize("bad", [(-1.0, R), (0.0, R), (R + 1e-7, R)])
def test_cap_half_angle_boundary(bad):
    with pytest.raises(ValueError):
        cap_half_angle(*bad)


def test_basal_cap_is_south_and_within_theta():
    m = SurfaceManifold.icosphere(3, R)
    basal, theta = select_basal_cap(m, footprint_radius=FOOT_R)
    assert basal.size > 0
    cent = np.asarray(m.tri_centroids)[basal]
    # all basal centroids are on the −z (south) side and within the cap angle.
    rr = np.linalg.norm(cent, axis=1)
    ang = np.arccos(np.clip(-cent[:, 2] / rr, -1.0, 1.0))
    assert np.all(ang <= theta + 1e-9)
    assert np.all(cent[:, 2] < 0)  # south hemisphere


def test_fa_placed_on_surface_in_cap_connected():
    surf = build_basal_surface(R=R, footprint_radius=FOOT_R, n_fa=60, subdivisions=3,
                               rng=np.random.default_rng(1))
    rep = basal_surface_report(surf, footprint_radius=FOOT_R)
    c = rep["controls"]
    assert c["on_surface"], c
    assert c["in_basal_cap"], c
    assert c["basal_connected"], c
    assert c["n_fa_placed"] == 60
    assert rep["verdict"] == "PASS"


def test_fa_positions_radial_equals_R():
    surf = build_basal_surface(R=R, footprint_radius=FOOT_R, n_fa=40, subdivisions=3,
                               rng=np.random.default_rng(2))
    radial = np.linalg.norm(surf.fa_positions, axis=1)
    assert np.allclose(radial, R, rtol=1e-12)


def test_n_fa_clamped_to_patch_count():
    # n_fa above the basal patch count clamps to the patch count.
    surf = build_basal_surface(R=R, footprint_radius=3.0e-6, n_fa=10_000,
                               subdivisions=3, rng=np.random.default_rng(3))
    assert surf.fa_positions.shape[0] == surf.basal_tris.shape[0]
    assert surf.fa_positions.shape[0] < 10_000


def test_footprint_below_resolution_raises():
    # a footprint smaller than the triangle scale selects 0 patches → raises.
    with pytest.raises(ValueError):
        build_basal_surface(R=R, footprint_radius=0.3e-6, n_fa=10, subdivisions=2)


def test_footprint_too_large_raises():
    with pytest.raises(ValueError):
        build_basal_surface(R=R, footprint_radius=R + 1e-7, n_fa=10)


def test_resolution_invariance_master_gate():
    """VG-1: the basal-cap area fraction is invariant to subdivision level."""
    res = cap_fraction_resolution_invariance(R=R, footprint_radius=FOOT_R,
                                             levels=(2, 3, 4))
    assert res["invariant"], res
    # converges toward the analytic cap fraction (1 − cos θ)/2.
    fr4 = res["fractions_by_level"][4]
    assert abs(fr4 - res["analytic_cap_fraction"]) < 0.05


def test_no_force_no_bond_surface_is_geometry_only():
    """The S-layer is geometry: the manifold owns no force/bond/tension attribute."""
    m = SurfaceManifold.icosphere(2, R)
    for forbidden in ("k_edge", "edge_springs", "tension", "area_modulus", "params"):
        assert not hasattr(m, forbidden)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
