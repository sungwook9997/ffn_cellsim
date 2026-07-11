"""Rate/event-driven T1 coarse-graining for the DCM — the biology-time route (PI 2026-07-12).

Design + Sanity Gates: ``docs/v2_audit/DCM_T1_RATE_COARSEGRAIN_DESIGN_2026-07-12.md``.

The timescale attack (``DCM_TIMESCALE_ATTACK_2026-07-11``) proved the jammed-unjamming biology-time gap is
FUNDAMENTAL: tissue flow/unjamming comes only from fast, sub-timestep T1 rearrangement events, so any large-dt
scheme on the fine-grained mechanics correctly freezes (BDF2 confirmed). This module supplies the slow rearrangement
as a *physical rate process* (KMC of T1 events) layered on the large-dt BDF2 mechanics.

This file holds the data-independent pieces (detection + measurement + the rate law form). The rate law is grounded,
not fitted to an outcome:
  - ``k0`` (attempt/gating frequency) is anchored to the endocytic junction-turnover rate k_endo ≈ 0.01–0.1 s⁻¹
    (KB-4.13) — the physical rate that gates a neighbour swap.
  - the barrier shape (rate → k0 at s0*=5.41, exponentially suppressed when jammed) reproduces the SPV/vertex law
    "T1 frequency high in the fluid phase, zero in the solid phase" (KB-5.12).
  - the barrier stiffness ``B`` is CALIBRATED against the T1 rate MEASURED from small-dt fine-grained references
    (``measure_t1_rate``), i.e. read from the mechanics — see the k_T1(s) calibration sweep.

Units: DCM positions are in METERS; rates in s⁻¹; the 3D shape index s = S/V^(2/3) (unjamming threshold s0*≈5.41).
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

S0_STAR_3D = 5.41          # Merkel–Manning 3D Voronoi unjamming threshold (KB s0* analog of the 2D p0*=3.81)
K_ENDO_LO = 0.01           # KB-4.13 endocytic junction-turnover rate [1/s] — lower bound of the k0 anchor
K_ENDO_HI = 0.10           # KB-4.13 upper bound
K0_DEFAULT = 0.03          # k0 attempt frequency [1/s]; geometric-mean-ish of the k_endo gate (anchored, not fitted)


def cell_centroids(pos: np.ndarray, cof: np.ndarray, n_cells: int) -> np.ndarray:
    """Per-cell centroid (mean node position). pos in metres → returns metres, shape (n_cells, 3)."""
    return np.array([pos[cof == c].mean(0) for c in range(n_cells)])


def auto_cutoff(centroids: np.ndarray, shell_factor: float = 1.15) -> float:
    """First-coordination-shell cutoff for cell–cell contact = shell_factor × median nearest-neighbour
    centroid distance. Robust to the packing scale (no hard-coded length)."""
    tree = cKDTree(centroids)
    dnn, _ = tree.query(centroids, k=2)          # k=2: self + nearest
    return float(shell_factor * np.median(dnn[:, 1]))


def cell_neighbor_set(centroids: np.ndarray, cutoff: float) -> set[tuple[int, int]]:
    """Set of contacting cell pairs (i<j) within ``cutoff`` — the tissue topology (who touches whom)."""
    return cKDTree(centroids).query_pairs(cutoff)


def measure_t1_rate(frames: np.ndarray, cof: np.ndarray, dt_frame: float,
                    shell_factor: float = 1.15) -> dict:
    """G3 grounding: measure the physical T1 (neighbour-exchange) rate from a fine-grained trajectory.

    A T1 = one neighbour lost + one gained. Between consecutive frames the symmetric difference of the contact-neighbour
    set counts (losses + gains) = 2 × T1 events. Returns the per-cell T1 rate [1/cell/s] and the raw counts.

    Args:
        frames: (n_frames, N, 3) node positions in metres.
        cof: (N,) cell-of-node index (cof<0 = dormant, ignored).
        dt_frame: physical time between frames [s].
        shell_factor: contact cutoff = shell_factor × median NN centroid distance.

    Returns:
        dict with ``rate_per_cell_s``, ``t1_events``, ``neighbour_changes``, ``cutoff_m``, ``n_cells``, ``T_total_s``.
    """
    cof = np.asarray(cof)
    n_cells = int(cof.max()) + 1
    cens = [cell_centroids(frames[k], cof, n_cells) for k in range(len(frames))]
    cutoff = auto_cutoff(cens[0], shell_factor)
    sets = [cell_neighbor_set(c, cutoff) for c in cens]
    changes = sum(len(sets[k] ^ sets[k - 1]) for k in range(1, len(sets)))
    t1_events = changes / 2.0
    T_total = dt_frame * (len(frames) - 1)
    rate = t1_events / (n_cells * T_total) if (n_cells * T_total) > 0 else 0.0
    return dict(rate_per_cell_s=rate, t1_events=t1_events, neighbour_changes=changes,
                cutoff_m=cutoff, n_cells=n_cells, T_total_s=T_total)


def k_t1_arrhenius(s: float | np.ndarray, k0: float = K0_DEFAULT, s_star: float = S0_STAR_3D,
                   barrier_stiffness: float = 1.0) -> float | np.ndarray:
    """Physical T1 rate law (NOT Metropolis): k_T1 = k0·exp(−B·(s0*−s)_+).

    - s ≥ s0* (fluid/unjammed): k_T1 = k0 (the junction-turnover gate, KB-4.13) — barrier vanished.
    - s < s0* (jammed): exponentially suppressed by the shape-index barrier (B·(s0*−s)), reproducing KB-5.12
      "T1 frequency high in fluid, zero in solid".

    k0 is anchored to k_endo (KB-4.13); ``barrier_stiffness`` B is calibrated from the measured k_T1(s) curve.
    """
    deficit = np.maximum(np.asarray(s_star, float) - np.asarray(s, float), 0.0)
    return k0 * np.exp(-barrier_stiffness * deficit)


def fit_barrier_stiffness(s_vals: np.ndarray, rate_vals: np.ndarray, k0: float = K0_DEFAULT,
                          s_star: float = S0_STAR_3D) -> float:
    """Calibrate B from measured (shape index, T1 rate) points by least-squares on ln(rate) = ln k0 − B·(s0*−s)_+.
    Uses only points with rate>0 (a jammed 0-rate point carries no barrier info). Grounds B in the mechanics."""
    s_vals = np.asarray(s_vals, float); rate_vals = np.asarray(rate_vals, float)
    m = rate_vals > 0
    deficit = np.maximum(s_star - s_vals[m], 0.0)
    y = np.log(rate_vals[m] / k0)                       # = −B·deficit
    # least-squares slope through the origin: B = −(Σ deficit·y)/(Σ deficit²)
    denom = float((deficit ** 2).sum())
    return float(-(deficit * y).sum() / denom) if denom > 0 else 1.0
