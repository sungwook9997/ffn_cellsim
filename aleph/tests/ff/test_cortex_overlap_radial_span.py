"""Smoothness-preserving cortex overlap resolution (``overlap_mode='radial_span'``) — CPU reference test.

The default ``'transverse'`` overlap relaxation shoves each sub-cutoff crossing node IN the shell tangent plane,
knocking it off its smooth great-circle arc → a per-node BENDING KINK. Invisible at the 0.5µm mesh (~0.8°) but
5° at the fine 75nm mesh, where cytosim's bending force F=α·d (α=κ/seg³) amplifies it ×288 — the fine-cortex
GATE-A residual plateau. ``'radial_span'`` instead separates crossing fibers OUT-of-plane (sphere-radial) with
the displacement spread over a ±span-node cosine bump, so each fiber stays a low-curvature arc while the two
fibers still separate in 3D past the WCA cutoff.

These tests certify, on a small synthetic FINE-mesh cortex (pure NumPy, no Warp/CUDA):
  (a) the DEFAULT mode is byte-identical to before this change (Gate-1 bit-identical-cortex contract);
  (b) 'radial_span' sharply reduces the max per-node turn-angle at crossing nodes vs 'transverse'; while
  (c) still resolving the interpenetrations (min cross-fiber node separation ≥ most of the WCA contact radius).
"""

from __future__ import annotations

import dataclasses

import numpy as np
from scipy.spatial import cKDTree

from aleph.laws.cortex_assembly import CortexParams, build_cortex_network, count_interpenetrations

_SIGMA_EV = 0.007
_R_C = 2.0 ** (1.0 / 6.0) * _SIGMA_EV        # WCA contact radius [µm]


def _fine_params() -> CortexParams:
    """A FINE 75nm-mesh cortex with physical radial thickness (where the transverse kink actually bites)."""
    return dataclasses.replace(CortexParams(), seg_um=0.075, beads_per_filament=41, cortex_thickness_um=0.2)


def _turn_angles(net) -> np.ndarray:
    """Per interior-node turn angle [deg] (angle between consecutive segments), indexed by GLOBAL node id."""
    off = np.asarray(net.fiber_offsets, np.int64)
    pos = np.asarray(net.pos, np.float64)
    ang = np.zeros(pos.shape[0])
    for f in range(len(off) - 1):
        seg = pos[off[f]:off[f + 1]]
        if seg.shape[0] < 3:
            continue
        v1 = seg[1:-1] - seg[:-2]
        v2 = seg[2:] - seg[1:-1]
        c = np.einsum("ij,ij->i", v1, v2) / (
            np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-30)
        ang[off[f] + 1:off[f + 1] - 1] = np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
    return ang


def _crossing_nodes(net) -> np.ndarray:
    """Global ids of nodes in a cross-fiber pair within r_c (the interpenetrations the relaxation must fix)."""
    off = np.asarray(net.fiber_offsets, np.int64)
    fid = np.repeat(np.arange(len(off) - 1), np.diff(off))
    pr = cKDTree(net.pos).query_pairs(_R_C, output_type="ndarray")
    if pr.shape[0]:
        pr = pr[fid[pr[:, 0]] != fid[pr[:, 1]]]
    return np.unique(pr.ravel()) if pr.shape[0] else np.zeros(0, np.int64)


def _min_cross_separation(net) -> float:
    """Smallest cross-fiber node distance [µm] (≈ r_c when overlaps are resolved)."""
    off = np.asarray(net.fiber_offsets, np.int64)
    fid = np.repeat(np.arange(len(off) - 1), np.diff(off))
    pos = np.asarray(net.pos, np.float64)
    # nearest neighbour on a DIFFERENT fiber, via a small k-NN sweep
    tree = cKDTree(pos)
    dist, idx = tree.query(pos, k=6)
    best = np.inf
    for col in range(1, dist.shape[1]):
        diff = fid[idx[:, col]] != fid
        if diff.any():
            best = min(best, float(dist[diff, col].min()))
    return best


def test_default_transverse_is_bit_identical() -> None:
    """Gate-1: the default overlap_mode reproduces the historical transverse relaxation bit-for-bit."""
    p = _fine_params()
    ref, _ = build_cortex_network(p, rng=np.random.default_rng(3), n_filaments=1500, resolve_overlaps=True)
    dflt, _ = build_cortex_network(p, rng=np.random.default_rng(3), n_filaments=1500, resolve_overlaps=True,
                                   overlap_mode="transverse")
    assert np.array_equal(ref.pos, dflt.pos), "default mode must be byte-identical to the historical relaxation"


def test_radial_span_reduces_crossing_kink_and_resolves_overlaps() -> None:
    """radial_span cuts the crossing-node kink vs transverse while still clearing the interpenetrations."""
    p = _fine_params()
    n = 1500
    raw, _ = build_cortex_network(p, rng=np.random.default_rng(4), n_filaments=n, resolve_overlaps=False)
    trans, mt = build_cortex_network(p, rng=np.random.default_rng(4), n_filaments=n, resolve_overlaps=True,
                                     overlap_mode="transverse")
    rad, mr = build_cortex_network(p, rng=np.random.default_rng(4), n_filaments=n, resolve_overlaps=True,
                                   overlap_mode="radial_span", overlap_span=2)

    xn = _crossing_nodes(raw)
    assert xn.size > 0, "the fine raw build must have cross-fiber interpenetrations to resolve"

    a_trans = _turn_angles(trans)[xn]
    a_rad = _turn_angles(rad)[xn]
    max_trans, max_rad = float(a_trans.max()), float(a_rad.max())
    mean_trans, mean_rad = float(a_trans.mean()), float(a_rad.mean())

    # (b) radial_span keeps crossing arcs far smoother than the in-plane transverse shove
    assert max_rad < 0.5 * max_trans, f"radial_span must cut peak crossing kink (trans {max_trans:.2f}° -> rad {max_rad:.2f}°)"
    assert mean_rad < 0.5 * mean_trans, f"radial_span must cut mean crossing kink (trans {mean_trans:.2f}° -> rad {mean_rad:.2f}°)"

    # (c) overlaps are still resolved: no residual interpenetrations, min separation ≥ most of r_c
    assert mr["n_interpenetrating_nodes"] == 0, f"radial_span must clear interpenetrations (got {mr['n_interpenetrating_nodes']})"
    assert _min_cross_separation(rad) >= 0.9 * _R_C, "radial_span min cross-fiber separation must reach ~r_c"
    # both modes reach the same overlap-free target
    assert count_interpenetrations(trans)["n_interpenetrating_nodes"] == 0

    print(f"crossing kink  max: transverse {max_trans:.2f}° -> radial_span {max_rad:.2f}°  "
          f"(mean {mean_trans:.2f}° -> {mean_rad:.2f}°); min-sep {_min_cross_separation(rad) / _R_C:.3f}·r_c")
