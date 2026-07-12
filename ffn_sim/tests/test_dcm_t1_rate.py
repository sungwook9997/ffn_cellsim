"""Unit tests for the T1 rate/event coarse-graining primitives (``ffn_sim.dcm.dcm_t1_rate``).

Self-contained (synthetic geometry — no npz dependency). Covers the detection (neighbour set, cutoff), the G3
measurement (a trajectory with a KNOWN neighbour exchange), and the grounded rate law + barrier calibration.
"""
import numpy as np
import pytest

from ffn_sim.dcm.dcm_t1_rate import (
    auto_cutoff, cell_neighbor_set, cell_centroids, measure_t1_rate,
    k_t1_arrhenius, fit_barrier_stiffness, t1_swap_partners, t1_geometric_move,
    t1_deformation_move, apply_cell_elongation, t1_kmc_step, K0_DEFAULT, S0_STAR_3D,
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


def test_t1_swap_partners_picks_nontouching_common_neighbours():
    # Proper T1 quartet: firing pair 0,1 touch; c=2 (above) and d=3 (below) each touch BOTH 0 and 1 but NOT each
    # other (they GAIN contact in the T1). cutoff 1.6: c-0,c-1,d-0,d-1 ≈ 1.58 (<1.6, touch); c-d = 3.0 (>1.6, apart).
    c = np.array([[0, 0, 0], [1, 0, 0], [0.5, 1.5, 0], [0.5, -1.5, 0]], float)
    nbrs = cell_neighbor_set(c, cutoff=1.6)          # {(0,1),(0,2),(1,2),(0,3),(1,3)} — NOT (2,3)
    sp = t1_swap_partners(nbrs, 0, 1, c)
    assert sp == (2, 3)                               # the non-touching common-neighbour pair
    # if 2,3 also touched (add (2,3)), they are no longer a swap-in pair → None
    assert t1_swap_partners(nbrs | {(2, 3)}, 0, 1, c) is None


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


def _dense_blob(n_side=4, spacing=1.0):
    """A dense cubic blob of cells (many contacting pairs + common neighbours → real T1 candidates)."""
    g = np.array([[i, j, k] for i in range(n_side) for j in range(n_side) for k in range(n_side)], float)
    return g * spacing


def test_kmc_step_zero_rate_fires_nothing():
    cen = _dense_blob()
    s = np.full(len(cen), 5.41)          # at threshold, but k0-scaled...
    rng = np.random.default_rng(0)
    out = t1_kmc_step(cen, s, cutoff=1.5, dt=1.0, rng=rng, k0=0.0)   # k0=0 → k_T1=0 → p=0
    assert out["n_fired"] == 0
    assert out["disp"] == {}


def test_kmc_step_high_rate_fires_and_moves():
    cen = _dense_blob()
    s = np.full(len(cen), 5.41)          # fluid → barrier 0 → rate=k0
    rng = np.random.default_rng(1)
    out = t1_kmc_step(cen, s, cutoff=1.5, dt=1.0, rng=rng, k0=1e6)   # p≈1 → (almost) every valid pair fires
    assert out["n_fired"] > 0
    # fired T1s produce non-zero displacements on the involved cells
    assert len(out["disp"]) > 0
    assert any(np.linalg.norm(v) > 0 for v in out["disp"].values())


def test_kmc_step_jammed_suppresses_rate():
    cen = _dense_blob()
    s_jammed = np.full(len(cen), 4.5)    # deep jammed → large barrier → far fewer firings than fluid
    rng1 = np.random.default_rng(2); rng2 = np.random.default_rng(2)
    fluid = t1_kmc_step(cen, np.full(len(cen), 5.41), cutoff=1.5, dt=1.0, rng=rng1, k0=0.05, barrier_stiffness=5.0)
    jammed = t1_kmc_step(cen, s_jammed, cutoff=1.5, dt=1.0, rng=rng2, k0=0.05, barrier_stiffness=5.0)
    assert jammed["n_fired"] <= fluid["n_fired"]


def test_apply_cell_elongation_is_volume_preserving_and_raises_aspect():
    rng = np.random.default_rng(0)
    sph = rng.normal(size=(300, 3)); sph /= np.linalg.norm(sph, axis=1, keepdims=True)   # unit sphere shell
    el = apply_cell_elongation(sph, np.array([1.0, 0, 0]), lam=1.4)
    # volume proxy = product of per-axis std (∝ ellipsoid volume) preserved (det of the scaling = 1)
    assert np.prod(el.std(0)) == pytest.approx(np.prod(sph.std(0)), rel=1e-6)
    # elongated along x → x-aspect rises above 1
    assert el[:, 0].std() / el[:, 1].std() > 1.2


def test_deformation_move_grounded_lambda_and_noop_when_touching():
    cen = np.array([[0, 0, 0], [3e-6, 0, 0]], float)     # c,d 3µm apart, cutoff 2µm → gap>0 → elongate
    m = t1_deformation_move(cen, 0, 1, cutoff=2e-6, r_cell=7.5e-6)
    assert set(m.keys()) == {0, 1}
    lam = m[0][1]
    assert lam > 1.0                                     # prolate
    # grounded: λ = 1 + (gap/2)/r_cell, gap = 3µm − 2µm·0.95 = 1.1µm → λ = 1 + 0.55/7.5
    assert lam == pytest.approx(1.0 + (3e-6 - 2e-6 * 0.95) / 2 / 7.5e-6, rel=1e-6)
    # already touching (gap<0) → no-op
    assert t1_deformation_move(np.array([[0, 0, 0], [1e-6, 0, 0]], float), 0, 1, cutoff=2e-6, r_cell=7.5e-6) == {}


def test_kmc_deform_mode_fires_and_returns_elongations():
    cen = _dense_blob() * 1e-6
    out = t1_kmc_step(cen, np.full(len(cen), 5.41), cutoff=1.8e-6, dt=1.0,
                      rng=np.random.default_rng(1), k0=1e6, move_mode="deform", r_cell=0.5e-6)
    assert out["n_fired"] > 0
    assert len(out["elong"]) > 0
    assert len(out["disp"]) == 0                         # deform mode uses elong, not disp
    ax, lam = next(iter(out["elong"].values()))
    assert lam > 1.0


def test_fit_barrier_recovers_known_stiffness():
    # generate rate points from a known B, confirm the fit recovers it
    B_true = 2.7
    s = np.array([4.6, 4.9, 5.2, 5.41])
    rate = k_t1_arrhenius(s, k0=K0_DEFAULT, barrier_stiffness=B_true)
    B_fit = fit_barrier_stiffness(s, rate, k0=K0_DEFAULT)
    assert B_fit == pytest.approx(B_true, rel=1e-6)
