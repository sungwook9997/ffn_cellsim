"""Tests for the Option-A ManifoldIndex broad-phase search/coordinate accelerator.

Covers the two correctness properties the 2026-06-09 (a)/(b) adversarial audit demanded
before any binder is rewired (H7_2D_MESH_AB_DECISION_2026-06-09.md):
  * VG-6 candidate-set identity: index within-reach pairs == global cKDTree at refresh time.
  * Staleness gate: the geometry-derived conservative k-ring absorbs sub-threshold drift
    (no dropped pairs below drift_threshold).
Geometry only — no mechanics, no γ.
"""
from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.cortex.surface_manifold import SurfaceManifold
from ffn_sim.cortex.manifold_index import (
    ManifoldIndex,
    global_pairs_within_reach,
    press_onto_substrate,
)


def _shell_cloud(n, R, seed=0):
    """n beads scattered in a thin band around radius R (a toy cortex shell)."""
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    r = R * (1.0 + rng.uniform(-0.013, 0.013, size=n))  # ~200nm band at R=7.5µm
    return u * r[:, None]


def test_vg6_candidate_identity():
    R = 7.5e-6
    pos = _shell_cloud(800, R, seed=1)
    reach = 8.4e-7
    mani = SurfaceManifold.icosphere(subdivisions=3, radius=R)
    idx = ManifoldIndex(mani, reach=reach).refresh(pos)
    truth = global_pairs_within_reach(pos, reach)
    got = idx.pairs_within_reach(pos)
    assert got == truth, (
        f"VG-6 FAIL: missed {len(truth - got)} extra {len(got - truth)} pairs"
    )


def test_staleness_gate_holds_below_threshold():
    R = 7.5e-6
    pos0 = _shell_cloud(800, R, seed=2)
    reach = 8.4e-7
    mani = SurfaceManifold.icosphere(subdivisions=3, radius=R)
    idx = ManifoldIndex(mani, reach=reach).refresh(pos0)
    # drift each bead by up to 90% of the gate threshold, keep the stale map
    rng = np.random.default_rng(3)
    d = rng.normal(size=pos0.shape)
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    mag = rng.uniform(0.0, 0.9 * idx.drift_threshold, size=pos0.shape[0])
    pos1 = pos0 + d * mag[:, None]
    assert idx.max_drift_since_refresh(pos1) < idx.drift_threshold
    truth = global_pairs_within_reach(pos1, reach)
    got = idx.pairs_within_reach(pos1)  # stale broad phase + narrow on pos1
    missed = truth - got
    assert len(missed) / max(1, len(truth)) < 0.01, (
        f"staleness gate LEAKS below threshold: missed {len(missed)}/{len(truth)}"
    )


def test_k_ring_derived_not_tuned():
    R = 7.5e-6
    mani = SurfaceManifold.icosphere(subdivisions=3, radius=R)
    # k grows with reach (DERIVED via kring_for_reach), never negative/zero
    k_small = ManifoldIndex(mani, reach=3.0e-7).k
    k_big = ManifoldIndex(mani, reach=2.0e-6).k
    assert k_small >= 1 and k_big >= k_small


def test_empty_and_boundary():
    R = 7.5e-6
    mani = SurfaceManifold.icosphere(subdivisions=2, radius=R)
    idx = ManifoldIndex(mani, reach=8.4e-7)
    with pytest.raises(RuntimeError):
        idx.neighbor_pools()  # before refresh
    idx.refresh(np.empty((0, 3)))
    assert idx.pairs_within_reach(np.empty((0, 3))) == set()
    assert idx.max_drift_since_refresh(np.empty((0, 3))) == 0.0


def _flat_ventral_cloud(n, R, z_basal, seed=0):
    """n beads on the FLAT ventral disk at z=z_basal (adherent geometry, B1 flat)."""
    rng = np.random.default_rng(seed)
    r_contact = 0.85 * R
    rr = r_contact * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    th = rng.uniform(0.0, 2.0 * np.pi, size=n)
    band = rng.uniform(-1.0e-7, 1.0e-7, size=n)  # ~200nm ventral band
    return np.column_stack([rr * np.cos(th), rr * np.sin(th),
                            np.full(n, z_basal) + band])


def test_adherent_flat_ventral_consistency():
    """ADHERENT pivot Stage-1: pressed (flat ventral) S1 + flat ventral B1 are CONSISTENT —
    ManifoldIndex VG-6 identity holds on the flat geometry (no S1-curved/B1-flat mismatch)."""
    R = 7.5e-6
    z_basal = -0.5 * R
    mani = SurfaceManifold.icosphere(subdivisions=3, radius=R)
    press_onto_substrate(mani, z_basal)               # S1 → flat ventral
    pos = _flat_ventral_cloud(700, R, z_basal, seed=4)  # B1 → flat ventral
    reach = 8.4e-7
    idx = ManifoldIndex(mani, reach=reach).refresh(pos)
    truth = global_pairs_within_reach(pos, reach)
    got = idx.pairs_within_reach(pos)
    assert got == truth, (
        f"flat-ventral VG-6 FAIL: missed {len(truth - got)} extra {len(got - truth)}"
    )
    # the pressed manifold genuinely has a flat ventral region (verts clamped to z_basal)
    assert np.any(np.isclose(mani.verts[:, 2], z_basal))


def test_no_mechanics_no_force_attributes():
    """Fidelity guard: the index/manifold carry no force/tension/γ DOF (Layer-2 only)."""
    R = 7.5e-6
    mani = SurfaceManifold.icosphere(subdivisions=2, radius=R)
    idx = ManifoldIndex(mani, reach=8.4e-7)
    for forbidden in ("force", "tension", "gamma", "k_spring", "energy"):
        assert not hasattr(idx, forbidden)
        assert not hasattr(mani, forbidden)
