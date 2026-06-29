"""Sanity tests for ``SurfaceManifold.fit_to_cloud`` (slaved-to-beads deformable shell).

The manifold must track a deforming bead cloud (uniform rescale for a sphere,
anisotropic for an ellipsoid) while keeping connectivity + the geometry-gate
invariants (outward normals, orthonormal frames, positive areas), and the region
masks must still resolve on a fitted (deformed) manifold.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.cortex.manifold_regions import basal_ring_region, polarized_patch_region
from ffn_sim.common.surface_manifold import SurfaceManifold

R_CELL = 7.5e-6


def _sphere_cloud(R, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v * R


def _ellipsoid_cloud(a, b, c, n=4000, seed=1):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v * np.array([a, b, c])


def _gates_hold(m: SurfaceManifold) -> None:
    # outward normals (n̂ · ĉ > 0), orthonormal right-handed frames, positive areas.
    cdir = m.tri_centroids / np.linalg.norm(m.tri_centroids, axis=1, keepdims=True)
    assert np.all(np.einsum("ij,ij->i", m.tri_normals, cdir) > 0.0)
    assert np.allclose(np.linalg.norm(m.tri_normals, axis=1), 1.0, atol=1e-9)
    assert np.allclose(np.einsum("ij,ij->i", m.tri_e1, m.tri_e2), 0.0, atol=1e-9)
    cross = np.cross(m.tri_e1, m.tri_e2)
    assert np.allclose(np.einsum("ij,ij->i", cross, m.tri_normals), 1.0, atol=1e-6)
    assert np.all(m.tri_areas > 0.0)


class TestFitToCloud:
    def test_sphere_rescale(self):
        m = SurfaceManifold.icosphere(3, R_CELL)
        tris0 = m.tris.copy()
        m.fit_to_cloud(_sphere_cloud(1.3 * R_CELL))
        # All vertices land at ≈ the new radius (uniform rescale).
        vr = np.linalg.norm(m.verts, axis=1)
        assert np.allclose(vr, 1.3 * R_CELL, rtol=0.02)
        # Connectivity unchanged.
        assert np.array_equal(m.tris, tris0)
        _gates_hold(m)

    def test_ellipsoid_tracking(self):
        m = SurfaceManifold.icosphere(4, R_CELL)
        a, b, c = 1.4 * R_CELL, 1.0 * R_CELL, 0.6 * R_CELL
        m.fit_to_cloud(_ellipsoid_cloud(a, b, c))
        # Vertices nearest each principal axis track that semi-axis.
        v = m.verts
        ax_x = v[np.argmax(np.abs(v[:, 0]))]
        ax_z = v[np.argmax(np.abs(v[:, 2]))]
        assert abs(abs(ax_x[0]) - a) / a < 0.08
        assert abs(abs(ax_z[2]) - c) / c < 0.12
        # The fitted shell is genuinely anisotropic (x-extent > z-extent).
        assert np.ptp(v[:, 0]) > 1.5 * np.ptp(v[:, 2])
        _gates_hold(m)

    def test_total_area_grows_with_radius(self):
        m = SurfaceManifold.icosphere(3, R_CELL)
        a0 = m.total_area()
        m.fit_to_cloud(_sphere_cloud(2.0 * R_CELL))
        # Area scales ≈ 4× for a 2× radius (within the flat-triangle deficit).
        assert m.total_area() / a0 == pytest.approx(4.0, rel=0.05)

    def test_regions_resolve_on_fitted_manifold(self):
        m = SurfaceManifold.icosphere(3, R_CELL)
        m.fit_to_cloud(_ellipsoid_cloud(1.3 * R_CELL, 1.3 * R_CELL, 0.7 * R_CELL))
        # Region masks still produce non-empty, localized selections on the
        # deformed shell (R_cell passed as the nominal scale).
        br = basal_ring_region(m, R_cell=R_CELL, rest_length=0.5e-6,
                               collar_half_angle=0.25)
        pp = polarized_patch_region(m, R_cell=R_CELL, half_angle_azimuth=math.pi / 6)
        assert br.patch_ids.size > 0
        assert pp.patch_ids.size > 0
        # in_plane_dir stays tangent on the deformed patches.
        assert np.all(np.abs(np.einsum("ij,ij->i", br.in_plane_dir, br.normals)) < 1e-9)

    def test_bad_inputs_raise(self):
        m = SurfaceManifold.icosphere(2, R_CELL)
        with pytest.raises(ValueError):
            m.fit_to_cloud(np.empty((0, 3)))
        with pytest.raises(ValueError):
            m.fit_to_cloud(_sphere_cloud(R_CELL), k_dirs=0)
