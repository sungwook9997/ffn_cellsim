"""Rate/event-driven T1 coarse-graining for the DCM — the biology-time route (PI 2026-07-12).

Design + Sanity Gates: ``docs/v2_audit/_historical/DCM_T1_RATE_COARSEGRAIN_DESIGN_2026-07-12.md``.

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


def t1_swap_partners(neighbors: set[tuple[int, int]], a: int, b: int,
                     centroids: np.ndarray) -> tuple[int, int] | None:
    """The cell pair (c, d) that swaps IN when the contacting pair (a, b) undergoes a T1.

    A T1 is the 3D analog of the 2D vertex flip: cells a,b share a junction and lose contact, while the two cells that
    neighbour BOTH a and b (the cells straddling that junction) gain contact. Among the common neighbours we pick the
    closest-approaching pair — the two most poised to touch after a,b separate. Returns None if <2 common neighbours
    (no well-defined swap, e.g. a surface junction).

    Args:
        neighbors: current contact-neighbour set (unordered i<j pairs).
        a, b: the contacting cell pair firing the T1.
        centroids: (n_cells, 3) cell centroids [m].
    """
    def nbrs(x):
        return {j for (i, j) in neighbors if i == x} | {i for (i, j) in neighbors if j == x}
    common = sorted(nbrs(a) & nbrs(b) - {a, b})
    if len(common) < 2:
        return None
    # The proper T1 swap-in pair (c,d) straddles the a–b junction but does NOT yet touch each other (they GAIN contact
    # in the T1). So pick the CLOSEST common-neighbour pair that is currently NON-adjacent; fall back to the closest
    # pair if all common neighbours already touch (degenerate — nothing to swap in).
    best, bestd = None, np.inf
    for i in range(len(common)):
        for j in range(i + 1, len(common)):
            ci, cj = common[i], common[j]
            if (min(ci, cj), max(ci, cj)) in neighbors:
                continue                                # already touching → not a swap-in pair
            d = float(np.linalg.norm(centroids[ci] - centroids[cj]))
            if d < bestd:
                bestd, best = d, (ci, cj)
    return best


def t1_geometric_move(centroids: np.ndarray, a: int, b: int, c: int, d: int, cutoff: float,
                      margin: float = 0.05) -> dict[int, np.ndarray]:
    """The discrete T1 geometric move: separate the firing pair (a,b), approach the swap pair (c,d).

    Returns per-cell centroid displacement vectors (applied as a RIGID translation of each cell's nodes — this preserves
    cell volume and shape, Sanity Gate G5). Each pair moves symmetrically (momentum-neutral, no spurious aggregate
    drift). Magnitudes come from the LOCAL geometry (cross the contact threshold by ``margin``), not a free parameter:
      - a,b separate along their axis until |a−b| = cutoff·(1+margin)  (contact just broken);
      - c,d approach along their axis until |c−d| = cutoff·(1−margin)  (contact just made).
    The resulting small overlaps with surrounding cells are re-resolved by the IPC relax in the BDF2 step. The move is a
    NO-OP direction (zero displacement) for a pair already on the correct side of the threshold.

    Args:
        centroids: (n_cells, 3) cell centroids [m].
        a, b: firing (separating) pair. c, d: swap-in (approaching) pair.
        cutoff: contact threshold [m].
        margin: fractional overshoot past the threshold.
    """
    disp = {a: np.zeros(3), b: np.zeros(3), c: np.zeros(3), d: np.zeros(3)}
    # separate a,b
    v = centroids[a] - centroids[b]; dab = float(np.linalg.norm(v))
    if dab > 1e-30:
        n = v / dab
        target = cutoff * (1.0 + margin)
        half = max(target - dab, 0.0) / 2.0        # only push apart (never pull together)
        disp[a] += half * n; disp[b] -= half * n
    # approach c,d
    w = centroids[c] - centroids[d]; dcd = float(np.linalg.norm(w))
    if dcd > 1e-30:
        n = w / dcd
        target = cutoff * (1.0 - margin)
        half = max(dcd - target, 0.0) / 2.0        # only pull together (never push apart)
        disp[c] -= half * n; disp[d] += half * n
    return disp


def t1_deformation_move(centroids: np.ndarray, c: int, d: int, cutoff: float, r_cell: float,
                        margin: float = 0.05) -> dict[int, tuple[np.ndarray, float]]:
    """Deformation-aware (s-RAISING) T1 move — the harden-first refinement of the rigid translate.

    The rigid ``t1_geometric_move`` translates whole cells into overlaps → the IPC relax rounds them → the shape index
    DROPS (the G4 24 s artifact). Real T1s instead ELONGATE the swapping cells as they squeeze past. Here the swap pair
    (c, d) each elongate ALONG the c–d axis toward each other (prolate, VOLUME-PRESERVING: λ along the axis, 1/√λ
    perpendicular), so their tips meet and form the new contact WITHOUT a gross whole-cell overlap — and the elongation
    RAISES their shape index (S/V^{2/3}), reproducing the T1's real cell strain. The separating pair (a, b) is NOT moved:
    inserting c,d between them stretches the a–b junction, which the force-based cadherin de-cohesion then ruptures
    (topology follows geometry). λ is GROUNDED, not free: each cell extends by half the gap so the tips just touch →
    λ = 1 + (gap/2)/r_cell, gap = |c−d| − cutoff·(1−margin). Returns {cell: (axis_unit, λ)}; no-op (empty) if already
    in contact.

    Args:
        centroids: (n_cells, 3) cell centroids [m]. c, d: the swap-in pair. cutoff: contact threshold [m].
        r_cell: cell radius [m] (sets the grounded elongation). margin: fractional contact overshoot.
    """
    w = centroids[c] - centroids[d]
    dcd = float(np.linalg.norm(w))
    gap = dcd - cutoff * (1.0 - margin)
    if dcd < 1e-30 or gap <= 0.0 or r_cell <= 0.0:
        return {}
    axis = w / dcd
    lam = 1.0 + (gap / 2.0) / r_cell          # extend each tip by half the gap → tips meet (grounded, volume-preserving)
    return {c: (axis, lam), d: (axis, lam)}


def apply_cell_elongation(nodes: np.ndarray, axis: np.ndarray, lam: float) -> np.ndarray:
    """Volume-preserving prolate scaling of a cell's nodes about their centroid: ×λ along ``axis``, ×1/√λ
    perpendicular (det = λ·(1/√λ)² = 1). Raises the cell's shape index without changing its volume (Sanity Gate G5)."""
    cen = nodes.mean(0)
    rel = nodes - cen
    a = axis / (np.linalg.norm(axis) + 1e-30)
    along = rel @ a                                    # (n,) component along the axis
    perp = rel - np.outer(along, a)                    # perpendicular part
    scaled = np.outer(along * lam, a) + perp / np.sqrt(lam)
    return cen + scaled


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


def t1_kmc_step(centroids: np.ndarray, shape_index: np.ndarray, cutoff: float, dt: float, rng,
                *, k0: float = K0_DEFAULT, s_star: float = S0_STAR_3D, barrier_stiffness: float = 1.0,
                margin: float = 0.05, move_mode: str = "rigid", r_cell: float = 7.5e-6) -> dict:
    """One KMC pass of T1 events over the current tissue — the coarse-grained rearrangement supplied at large dt.

    For each contacting cell pair (a,b): the local shape index s = min(s_a, s_b) (the pair unjams only if BOTH sides
    can) sets the rate k_T1 (``k_t1_arrhenius``); the pair fires with p = 1−exp(−k_T1·dt); a fired pair executes a T1
    geometric move (separate a,b; approach their swap partners). Displacements are accumulated per cell (a cell touched
    by several T1s in one pass sums them). Returns the per-cell displacement dict + the list of fired (a,b) pairs.

    Args:
        centroids: (n_cells, 3) cell centroids [m].
        shape_index: (n_cells,) per-cell 3D shape index s.
        cutoff: contact threshold [m].
        dt: the KMC cadence interval [s].
        rng: a numpy Generator (seeded by the caller for reproducibility).
        k0, s_star, barrier_stiffness, margin: rate-law + move parameters (k0/B grounded per the design).
    """
    neighbors = cell_neighbor_set(centroids, cutoff)
    disp: dict[int, np.ndarray] = {}
    elong: dict[int, tuple[np.ndarray, float]] = {}      # deform mode: {cell: (axis, λ)}
    fired: list[tuple[int, int]] = []
    for (a, b) in neighbors:
        s_pair = float(min(shape_index[a], shape_index[b]))
        k = float(k_t1_arrhenius(s_pair, k0, s_star, barrier_stiffness))
        p = 1.0 - np.exp(-k * dt)
        if rng.random() < p:
            partners = t1_swap_partners(neighbors, a, b, centroids)
            if partners is None:
                continue
            c, d = partners
            if move_mode == "deform":
                # STABILITY CAP (breaks the elongate→s↑→rate↑→elongate runaway): a cell that has already reached the
                # unjamming threshold s0* is fluid — no further elongation (cells don't get MORE fluid than s0*). This
                # bounds s at ~s0*, the physical transition. Skip the geometric move if either swap cell is unjammed.
                if shape_index[c] >= s_star or shape_index[d] >= s_star:
                    fired.append((a, b)); continue
                move = t1_deformation_move(centroids, c, d, cutoff, r_cell, margin)
                if not move:
                    continue
                # ONE elongation per cell per pass (no compounding — the accumulation was the runaway); keep the largest.
                for cell, (axis, lam) in move.items():
                    if cell not in elong or lam > elong[cell][1]:
                        elong[cell] = (axis, lam)
            else:
                move = t1_geometric_move(centroids, a, b, c, d, cutoff, margin)
                for cell, dv in move.items():
                    disp[cell] = disp.get(cell, np.zeros(3)) + dv
            fired.append((a, b))
    return dict(disp=disp, elong=elong, fired=fired, n_fired=len(fired))


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
