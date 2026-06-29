"""CI guard for the surface-manifold grid-invariance MASTER GATE.

The registry records that the *mesh-resolution grid-invariance master gate must
pass before any manifold force is enabled*. The full sweep lives in
``scripts/compartment_smoke/surface_manifold.py``; this is a fast pytest subset
that locks the contract in CI:

  * total flat-triangle area converges toward 4πR² as the mesh refines;
  * the DERIVED ``kring_for_reach`` produces a geodesic ring that COVERS every
    triangle within the physical reach (coverage == 1.0) at EVERY resolution —
    the grid-invariance the gate names;
  * the compartment stays GEOMETRY-ONLY (zero force-bearing bonds/particles).

GEOMETRY-ONLY: no force is built; these are pure geometry checks.
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.cell.compartment_registry import REGISTRY
from ffn_sim.common.surface_manifold import SurfaceManifold

_R = 7.5e-6
_REACH = 1.0e-6
_SUBDIVS = [1, 2, 3]   # fast subset (the script sweeps 0-4)


def _centroids(m: SurfaceManifold) -> np.ndarray:
    return m.verts[m.tris].mean(axis=1)


def test_area_converges_to_sphere():
    sphere = 4.0 * np.pi * _R**2
    errs = []
    for s in _SUBDIVS:
        m = SurfaceManifold.icosphere(s, _R)
        errs.append(abs(m.total_area() - sphere) / sphere)
    # monotone non-increasing and strictly better at the finest resolution.
    assert all(errs[i + 1] <= errs[i] + 1e-12 for i in range(len(errs) - 1))
    assert errs[-1] < errs[0]
    assert errs[-1] < 0.02   # <2% by subdiv 3


@pytest.mark.parametrize("subdiv", _SUBDIVS)
def test_kring_reach_coverage_master_gate(subdiv):
    m = SurfaceManifold.icosphere(subdiv, _R)
    k = int(m.kring_for_reach(_REACH))
    cents = _centroids(m)
    rng = np.random.default_rng(0)
    sample = rng.choice(m.n_tri, size=min(25, m.n_tri), replace=False)
    for t in sample:
        within = np.flatnonzero(np.linalg.norm(cents - cents[t], axis=1) <= _REACH)
        ring = set(m.patch_kring(int(t), k).tolist())
        # the geodesic ring must cover EVERY triangle within the reach ball.
        assert all(int(u) in ring for u in within)


def test_broad_phase_frame_converges():
    rng = np.random.default_rng(1)
    dirs = rng.normal(size=(300, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    pts = _R * dirs
    errs = []
    for s in _SUBDIVS:
        m = SurfaceManifold.icosphere(s, _R)
        patch = m.nearest_patch(pts)
        normals = np.array([m.frame(int(t))[0] for t in patch])
        cos = np.abs(np.sum(normals * dirs, axis=1))
        errs.append(float(np.mean(np.arccos(np.clip(cos, -1, 1)))))
    assert errs[-1] < errs[0]   # home-patch normal aligns better as mesh refines


def test_surface_manifold_is_geometry_only():
    spec = REGISTRY.get("surface_manifold")
    assert spec.geometry_only
    assert spec.performance_contract.n_bonds in (0, "0")
    assert not spec.performance_contract.particle_types_added
    REGISTRY.validate_manifold_geometry_only()   # registry-wide invariant
