"""Unit tests for the T1 rate/event coarse-graining primitives (``ffn_sim.dcm.dcm_t1_rate``).

Self-contained (synthetic geometry — no npz dependency). Covers the detection (neighbour set, cutoff), the G3
measurement (a trajectory with a KNOWN neighbour exchange), and the grounded rate law + barrier calibration.
"""
import numpy as np
import pytest

from ffn_sim.dcm.dcm_t1_rate import (
    auto_cutoff, cell_neighbor_set, cell_centroids, measure_t1_rate,
    k_t1_arrhenius, fit_barrier_stiffness, K0_DEFAULT, S0_STAR_3D,
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


def test_fit_barrier_recovers_known_stiffness():
    # generate rate points from a known B, confirm the fit recovers it
    B_true = 2.7
    s = np.array([4.6, 4.9, 5.2, 5.41])
    rate = k_t1_arrhenius(s, k0=K0_DEFAULT, barrier_stiffness=B_true)
    B_fit = fit_barrier_stiffness(s, rate, k0=K0_DEFAULT)
    assert B_fit == pytest.approx(B_true, rel=1e-6)
