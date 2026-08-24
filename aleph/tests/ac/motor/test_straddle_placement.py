"""Head-on-actin STRADDLE placement (defect#3 fidelity) — pure-NumPy geometry, no Warp/CUDA.

Certifies that the bipolar minifilament placement/orientation degree of freedom
(:meth:`MinifilamentTopology.placed_positions` + :func:`straddle_frame`, and the assemble-side partner
selection) puts the ± heads ON two anti-parallel actin filaments while preserving the rigid minifilament
geometry (so the backbone-bend / head-arm angle springs stay force-free at rest, exactly as the legacy
:meth:`MinifilamentTopology.positions`).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.assemble import (
    NMII_CAPTURE_UM,
    NMII_HEAD_OFFSET_UM,
    NMII_L_BB_UM,
    NMII_N_BB,
    NMII_N_SIDE,
    _actin_node_tangents,
    _select_antiparallel_partner,
)
from aleph.components.motor.minifilament_topology import MinifilamentTopology, straddle_frame


def _topo() -> MinifilamentTopology:
    return MinifilamentTopology(NMII_N_BB, NMII_N_SIDE, NMII_L_BB_UM, NMII_HEAD_OFFSET_UM)


def _pairwise(x: np.ndarray) -> np.ndarray:
    return np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)


# ── placed_positions: a rigid re-orientation of positions() ──────────────────────────────────────────
def test_placed_positions_is_rigid_transform() -> None:
    """placed_positions must preserve every pairwise distance (⇒ all internal geometry / rest springs)."""
    t = _topo()
    local = t.positions()
    placed = t.placed_positions(
        centre=np.array([1.0, -2.0, 0.5]),
        e_x=np.array([0.3, 1.0, -0.4]),           # arbitrary (non-unit) frame
        e_y=np.array([1.0, 0.2, 0.1]),
    )
    assert np.allclose(_pairwise(local), _pairwise(placed), atol=1e-12)


def test_placed_positions_frame_is_orthonormal_and_centred() -> None:
    """+x → e_x (unit), +y → e_y (re-orthogonalised, unit), backbone stays in the e_x–e_y plane."""
    t = _topo()
    centre = np.array([0.0, 0.0, 0.0])
    e_x = np.array([1.0, 0.0, 0.0])
    e_y = np.array([0.3, 1.0, 0.0])               # deliberately not orthogonal to e_x
    p = t.placed_positions(centre, e_x, e_y)
    bb = p[: t.n_bb]
    axis = bb[-1] - bb[0]
    axis /= np.linalg.norm(axis)
    assert np.allclose(axis, [1.0, 0.0, 0.0], atol=1e-12)     # backbone along e_x
    assert np.allclose(bb.mean(axis=0), centre, atol=1e-12)   # centred


def test_placed_positions_backbone_straight_arms_perpendicular() -> None:
    """Rest geometry unchanged vs positions(): backbone collinear, head arms ⟂ backbone (force-free)."""
    t = _topo()
    p = t.placed_positions(np.zeros(3), np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]))
    bb = p[: t.n_bb]
    seg = np.diff(bb, axis=0)
    seg /= np.linalg.norm(seg, axis=1, keepdims=True)
    assert np.allclose(seg, seg[0], atol=1e-12)               # every backbone segment parallel ⇒ straight
    # each head's arm is perpendicular to the local backbone axis (rest θ = π/2)
    for head_idx, bead_idx in t.head_backbone_bonds():
        arm = p[head_idx] - p[bead_idx]
        assert abs(np.dot(arm / np.linalg.norm(arm), seg[0])) < 1e-9


# ── straddle_frame: heads land on the two anti-parallel filament lines ────────────────────────────────
def test_straddle_frame_orthonormal() -> None:
    c, ex, ey = straddle_frame(
        np.array([0.0, 0.2, 0.0]), np.array([1.0, 0.0, 0.0]),
        np.array([0.0, -0.2, 0.0]), np.array([-1.0, 0.0, 0.0]))
    assert np.isclose(np.linalg.norm(ex), 1.0)
    assert np.isclose(np.linalg.norm(ey), 1.0)
    assert abs(np.dot(ex, ey)) < 1e-12


def test_straddle_heads_on_both_filaments_at_2offset() -> None:
    """Two anti-parallel filaments 2·offset apart ⇒ ± heads land ON their lines (perp residual → 0)."""
    t = _topo()
    off = NMII_HEAD_OFFSET_UM
    a = np.array([0.0, off, 0.0]); t_a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, -off, 0.0]); t_b = np.array([-1.0, 0.0, 0.0])
    c, ex, ey = straddle_frame(a, t_a, b, t_b)
    p = t.placed_positions(c, ex, ey)
    hp = p[t.n_bb: t.n_bb + t.n_heads_per_side]        # + heads → line A (y = +off)
    hm = p[t.n_bb + t.n_heads_per_side:]               # − heads → line B (y = −off)
    assert np.max(np.abs(hp[:, 1] - off)) < 1e-9       # + heads exactly on A
    assert np.max(np.abs(hm[:, 1] + off)) < 1e-9       # − heads exactly on B
    assert np.max(np.abs(p[:, 2])) < 1e-9              # whole unit stays in the filament plane


def test_straddle_frame_degenerate_inputs_do_not_raise() -> None:
    """Parallel tangents / partner on A's line ⇒ still a valid orthonormal frame (graceful fallback)."""
    for t_a, t_b, b in [
        (np.array([1.0, 0, 0]), np.array([1.0, 0, 0]), np.array([0.0, 0.4, 0.0])),   # parallel
        (np.zeros(3), np.zeros(3), np.array([0.4, 0.0, 0.0])),                        # zero tangents
        (np.array([1.0, 0, 0]), np.array([-1.0, 0, 0]), np.array([0.4, 0.0, 0.0])),  # b on A's line
    ]:
        c, ex, ey = straddle_frame(np.zeros(3), t_a, b, t_b)
        assert np.isclose(np.linalg.norm(ex), 1.0) and np.isclose(np.linalg.norm(ey), 1.0)
        assert abs(np.dot(ex, ey)) < 1e-9


# ── assemble-side helpers: tangents + anti-parallel partner selection on a synthetic mesh ─────────────
def _two_antiparallel_fibers(sep: float, n: int = 7, seg: float = 0.5):
    """Two straight fibers along ±x, separated `sep` in y, opposite node ordering (anti-parallel)."""
    xs = (np.arange(n) - (n - 1) / 2) * seg
    fib_a = np.column_stack([xs, np.full(n, sep / 2), np.zeros(n)])
    fib_b = np.column_stack([-xs, np.full(n, -sep / 2), np.zeros(n)])   # reversed ⇒ anti-parallel
    pos = np.vstack([fib_a, fib_b])
    offsets = np.array([0, n, 2 * n], np.int64)
    polarity = np.array([1, 1], np.int64)              # barbed at last node of each
    return pos, offsets, polarity, n


def test_actin_node_tangents_barbed_oriented() -> None:
    pos, off, pol, n = _two_antiparallel_fibers(sep=0.4)
    tang = _actin_node_tangents(pos, off, pol)
    assert np.allclose(np.linalg.norm(tang, axis=1), 1.0)
    assert np.allclose(tang[:n], [1.0, 0.0, 0.0], atol=1e-9)     # fiber A points +x
    assert np.allclose(tang[n:], [-1.0, 0.0, 0.0], atol=1e-9)    # fiber B points −x (anti-parallel)


def test_select_partner_picks_antiparallel_at_2offset() -> None:
    from scipy.spatial import cKDTree
    pos, off, pol, n = _two_antiparallel_fibers(sep=2 * NMII_HEAD_OFFSET_UM)   # ideal 2·offset separation
    tang = _actin_node_tangents(pos, off, pol)
    node_fiber = np.repeat(np.arange(2), np.diff(off))
    tree = cKDTree(pos)
    anchor = n // 2                                     # middle node of fiber A
    partner = _select_antiparallel_partner(tree, pos, tang, node_fiber, anchor, NMII_HEAD_OFFSET_UM)
    assert partner >= n                                 # picked a node on fiber B
    assert np.dot(tang[anchor], tang[partner]) < -0.9   # strongly anti-parallel


def test_select_partner_none_when_no_antiparallel_neighbour() -> None:
    """A single fiber (no cross-fiber partner) ⇒ −1 (caller falls back to legacy placement)."""
    from scipy.spatial import cKDTree
    xs = (np.arange(7) - 3) * 0.5
    pos = np.column_stack([xs, np.zeros(7), np.zeros(7)])
    off = np.array([0, 7], np.int64)
    tang = _actin_node_tangents(pos, off, np.array([1], np.int64))
    tree = cKDTree(pos)
    assert _select_antiparallel_partner(tree, pos, tang, np.zeros(7, int), 3, NMII_HEAD_OFFSET_UM) == -1


def _head_segment_distances(heads: np.ndarray, pos: np.ndarray, off: np.ndarray) -> np.ndarray:
    """Min distance of each head to any actin SEGMENT (the runtime capture metric — point-on-segment)."""
    seg_a = np.concatenate([np.arange(off[f], off[f + 1] - 1) for f in range(off.size - 1)])
    P, Q = pos[seg_a], pos[seg_a + 1]
    D = Q - P
    L2 = np.einsum("ij,ij->i", D, D)
    L2[L2 < 1e-18] = 1e-18
    out = np.empty(heads.shape[0])
    for h in range(heads.shape[0]):
        w = np.clip(np.einsum("ij,ij->i", heads[h] - P, D) / L2, 0.0, 1.0)
        proj = P + w[:, None] * D
        out[h] = np.min(np.linalg.norm(heads[h] - proj, axis=1))
    return out


def test_straddle_places_heads_within_capture() -> None:
    """End-to-end on the synthetic mesh: after straddle placement the heads sit within the capture radius.

    Distance is head→nearest-SEGMENT (point-on-segment), the metric the runtime :class:`SegmentQuery` uses —
    NOT head→nearest-node (discrete nodes at 0.5µm spacing would report a head sitting ON the segment as far).
    """
    from scipy.spatial import cKDTree
    t = _topo()
    pos, off, pol, n = _two_antiparallel_fibers(sep=2 * NMII_HEAD_OFFSET_UM)
    tang = _actin_node_tangents(pos, off, pol)
    node_fiber = np.repeat(np.arange(2), np.diff(off))
    tree = cKDTree(pos)
    anchor = n // 2
    partner = _select_antiparallel_partner(tree, pos, tang, node_fiber, anchor, NMII_HEAD_OFFSET_UM)
    c, ex, ey = straddle_frame(pos[anchor], tang[anchor], pos[partner], tang[partner])
    heads = t.placed_positions(c, ex, ey)[t.n_bb:]
    d = _head_segment_distances(heads, pos, off)
    assert (d <= NMII_CAPTURE_UM).mean() > 0.99         # essentially all heads within capture of a segment
