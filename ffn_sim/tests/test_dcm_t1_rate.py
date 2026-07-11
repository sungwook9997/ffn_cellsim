"""Unit tests for the T1 rate/event coarse-graining primitives (``ffn_sim.dcm.dcm_t1_rate``).

Self-contained (synthetic geometry — no npz dependency). Covers the detection (neighbour set, cutoff), the G3
measurement (a trajectory with a KNOWN neighbour exchange), and the grounded rate law + barrier calibration.
"""
import numpy as np
import pytest

from ffn_sim.dcm.dcm_t1_rate import (
    auto_cutoff, cell_neighbor_set, cell_centroids, measure_t1_rate,
    k_t1_arrhenius, fit_barrier_stiffness, t1_swap_partners, t1_geometric_move,
    K0_DEFAULT, S0_STAR_3D,
)


def _square_cells(spacing: float, jitter: np.ndarray | None = None) -> np.ndarray:
    """4 cell centroids on a unit square of side ``spacing`` (+ optional per-cell offset)."""
    c = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], float) * spacing
    return c + (jitter if jitter is not None else 0.0)


def test_auto_cutoff_scales_with_packing():
    c = _square_cells(2.0)
    # nearest-neighbour distance = 2.0 → cutoff = 1.15 × 2.0
    assert auto_cutoff(c) == pytest.approx(1.15 * 2.0, rel=1e-6)


def test_neighbor_set_edges_only():
    c = _square_cells(1.0)
    # cutoff just above the side (1.0) but below the diagonal (√2≈1.414) → 4 edges, no diagonals
    s = cell_neighbor_set(c, cutoff=1.2)
    assert s == {(0, 1), (0, 2), (1, 3), (2, 3)}


def test_measure_t1_rate_counts_a_known_exchange():
    # Build a 2-frame trajectory of 4 single-node "cells". Frame 0: a square (4 edges). Frame 1: pull cells 0,3
    # together so the 0-3 diagonal becomes a contact (gain) — one neighbour change.
    cof = np.arange(4)
    f0 = _square_cells(1.0)
    f1 = f0.copy()
    f1[0] = [0.45, 0.45, 0.0]   # cell 0 moves toward cell 3 → 0-3 now within cutoff
    f1[3] = [0.55, 0.55, 0.0]
    frames = np.stack([f0, f1])
    r = measure_t1_rate(frames, cof, dt_frame=1.0)
    # cutoff from frame 0 (NN=1.0 → 1.15). Frame 1 gains the (0,3) pair; cells 0/3 may lose an edge too.
    assert r["neighbour_changes"] >= 1
    assert r["rate_per_cell_s"] > 0.0
    assert r["n_cells"] == 4


def test_measure_t1_rate_static_is_zero():
    cof = np.arange(4)
    f = _square_cells(1.0)
    frames = np.stack([f, f, f])          # no motion → no exchanges → zero rate (the FROZEN case)
    r = measure_t1_rate(frames, cof, dt_frame=1.0)
    assert r["neighbour_changes"] == 0
    assert r["rate_per_cell_s"] == 0.0


def test_rate_law_saturates_at_threshold_and_suppresses_below():
    # at/above s0* the barrier has vanished → k_T1 = k0 (the junction-turnover gate)
    assert k_t1_arrhenius(S0_STAR_3D) == pytest.approx(K0_DEFAULT)
    assert k_t1_arrhenius(S0_STAR_3D + 0.5) == pytest.approx(K0_DEFAULT)   # clamped (deficit=0)
    # below s0* the rate is strictly suppressed and monotone in s
    assert k_t1_arrhenius(5.0, barrier_stiffness=2.0) < K0_DEFAULT
    assert k_t1_arrhenius(4.5, barrier_stiffness=2.0) < k_t1_arrhenius(5.0, barrier_stiffness=2.0)


def test_t1_swap_partners_picks_common_neighbours():
    # 4 cells in a square: the contacting pair (0,1) [bottom edge] has common neighbours {2,3} (both touch 0 and 1).
    c = _square_cells(1.0)
    nbrs = cell_neighbor_set(c, cutoff=1.2)     # {(0,1),(0,2),(1,3),(2,3)}
    # 0 neighbours {1,2}; 1 neighbours {0,3}; common(0,1) needs cells touching BOTH → none here (square has no
    # cell touching both 0 and 1), so expect None. Add a centre cell 4 touching all → common becomes {4,...}.
    assert t1_swap_partners(nbrs, 0, 1, c) is None
    # now a triangle bipyramid-ish: cells 2 and 3 both neighbour 0 and 1
    nbrs2 = nbrs | {(0, 3), (1, 2)}             # make 2,3 each touch both 0 and 1
    sp = t1_swap_partners(nbrs2, 0, 1, c)
    assert sp == (2, 3)


def test_t1_geometric_move_crosses_threshold_and_conserves_momentum():
    # a,b touching (dist 0.5 < cutoff 1.0) must separate past cutoff; c,d apart (dist 1.5 > cutoff) must approach below.
    cen = np.array([[0, 0, 0], [0.5, 0, 0],       # a,b close (to separate)
                    [0, 3.0, 0], [1.5, 3.0, 0]], float)  # c,d far (to approach)
    cutoff = 1.0
    disp = t1_geometric_move(cen, 0, 1, 2, 3, cutoff, margin=0.05)
    new = cen.copy()
    for k, dv in disp.items():
        new[k] = cen[k] + dv
    dab = np.linalg.norm(new[0] - new[1]); dcd = np.linalg.norm(new[2] - new[3])
    assert dab == pytest.approx(cutoff * 1.05, rel=1e-6)     # a,b just separated
    assert dcd == pytest.approx(cutoff * 0.95, rel=1e-6)     # c,d just joined
    # each pair momentum-neutral → total displacement ≈ 0
    assert np.allclose(sum(disp.values()), 0.0, atol=1e-12)


def test_t1_geometric_move_noop_when_already_correct():
    # a,b already separated (dist 2>cutoff) and c,d already touching (dist 0.5<cutoff) → zero move
    cen = np.array([[0, 0, 0], [2.0, 0, 0], [0, 3, 0], [0.5, 3, 0]], float)
    disp = t1_geometric_move(cen, 0, 1, 2, 3, cutoff=1.0)
    assert all(np.allclose(v, 0.0) for v in disp.values())


def test_fit_barrier_recovers_known_stiffness():
    # generate rate points from a known B, confirm the fit recovers it
    B_true = 2.7
    s = np.array([4.6, 4.9, 5.2, 5.41])
    rate = k_t1_arrhenius(s, k0=K0_DEFAULT, barrier_stiffness=B_true)
    B_fit = fit_barrier_stiffness(s, rate, k0=K0_DEFAULT)
    assert B_fit == pytest.approx(B_true, rel=1e-6)
