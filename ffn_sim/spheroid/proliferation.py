"""L2.4 contact-inhibited proliferation: the mechanistic, size-dependent spreading driver.

The minimal CBM (cohesion + edge-traction, ``cbm.py`` + ``spreading.py``) is COHESION-LOCKED
— A/A0 ≈ 1 at the measured MCF7 scales and the PI law ``A/A0 = a + b/R + c/R²`` does NOT
emerge (a model-limit finding, L2.3 REPORT). The experiment runs over **days**; spreading is
partly **proliferation-driven**. This module adds proliferation as a fine-grained MECHANISM,
not a fitted term:

* Each cell carries a stochastic **cell-cycle timer** (mean = MCF7 doubling time, per-cell CV).
* A cell **divides** when its timer elapses **and** it sits on a **free surface** — i.e. its
  first-coordination shell is not full (fewer than the sphere kissing number Z=12 neighbours).
  This is the canonical Drasdo-Höhme "free-space" division rule: a buried bulk cell (≈12
  contacts) is **contact-inhibited / quiescent**; a rim cell (a free face) proliferates.
* On division the cell **buds** a daughter into its least-crowded direction at the rest
  separation r0 (no repulsive-core blow-up, no gap).

Because the proliferating cells form a **rim of ~constant thickness**, the proliferating
fraction ∝ surface/volume ∝ **1/R**. That geometric surface-to-volume ratio is the EMERGENT
origin of the law's 1/R (and curvature 1/R²) terms — small spheroids grow more (relatively)
than large ones. Nothing is hard-coded to 1/R; it falls out of where division is allowed.

HOOMD has a fixed particle count per ``State``, so growth runs in **epochs**: relax the
current population for ``epoch_steps``, read positions back, apply the divisions that became
due, then rebuild a fresh HOOMD state (``cbm.build_cbm_simulation(positions=…)``) for the
next epoch. Overdamped BAOAB-limit dynamics have no persistent velocity, so rebuilding loses
no physical state. The frozen integrator + single ``md.pair.Morse`` are reused verbatim.

Sanity Gate
-----------
- Dimensional: positions [m]; ages/targets/times [s]; shell cutoff & split [m]; Z, counts [–].
- Boundary (proliferation OFF): with an infinite cycle time no division fires ⇒ this reduces
  exactly to the G1 stable aggregate (count conserved). Verified in tests.
- Boundary (dilute limit): an isolated cell (0 neighbours < Z) divides every cycle ⇒
  N(t) → 2^(t/τ) exponential. Verified in tests.
- Sign-sense: more crowding (higher first-shell count) ⇒ FEWER cells eligible ⇒ slower growth;
  division is biased to the RIM (above-median radial position), never the buried core.
- Conservation: each division is +1 cell (one parent → two daughters); COM budding shift is
  ≤ r0/2 (one daughter inherits the parent site, the other buds at the rest separation).
- Measurement-protocol: A(t)/A0 uses the same ``projected_area`` convex hull as G1/G3.
"""

from __future__ import annotations

from typing import Any

import hoomd
import numpy as np
import numpy.typing as npt
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from ffn_sim.spheroid.cbm import build_cbm_simulation, get_positions
from ffn_sim.spheroid.observables import (
    projected_area,
    radius_of_gyration,
)
from ffn_sim.spheroid.params import ResolvedL2, ResolvedProliferation

__all__ = [
    "first_shell_counts",
    "candidate_directions",
    "best_bud_direction",
    "sample_cycle_targets",
    "apply_divisions",
    "run_growth",
]

# Number of candidate bud directions probed for free space (a NUMERICAL sampling count, not
# physics). A Fibonacci sphere of 12 points ~ evenly tiles the 4π solid angle at the kissing
# scale — enough to find a free face on a rim cell.
_N_BUD_CANDIDATES = 12


def first_shell_counts(
    positions: npt.NDArray[np.float64], shell_cutoff: float
) -> npt.NDArray[np.int64]:
    """First-coordination-shell neighbour count per cell (excludes self).

    Args:
        positions: (N, 3) cell centers (m).
        shell_cutoff: first-shell radius (m); cells within this distance are "contacts".

    Returns:
        (N,) integer neighbour counts (number of OTHER cells within ``shell_cutoff``).
    """
    if shell_cutoff <= 0.0:
        raise ValueError("shell_cutoff must be strictly positive.")
    p = np.asarray(positions, dtype=np.float64)
    if p.shape[0] < 2:
        return np.zeros(p.shape[0], dtype=np.int64)
    tree = cKDTree(p)
    # query_ball_point with count_only includes self → subtract 1.
    counts = tree.query_ball_point(p, shell_cutoff, return_length=True) - 1
    return counts.astype(np.int64)


def candidate_directions(n: int = _N_BUD_CANDIDATES) -> npt.NDArray[np.float64]:
    """``n`` ~evenly distributed unit vectors on the sphere (deterministic Fibonacci sphere)."""
    if n < 1:
        raise ValueError("n must be >= 1.")
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)              # polar angle, area-uniform
    golden = np.pi * (1.0 + 5.0**0.5)               # golden-angle azimuth
    theta = golden * i
    return np.column_stack(
        [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]
    )


def best_bud_direction(
    positions: npt.NDArray[np.float64],
    i: int,
    neighbor_idx: npt.NDArray[np.int64],
    split: float,
    *,
    candidates: npt.NDArray[np.float64] | None = None,
) -> tuple[npt.NDArray[np.float64], float]:
    """Freest bud direction for cell ``i`` and the free-space gap it achieves.

    Probes a fixed set of candidate directions; for each, the daughter would sit at
    ``positions[i] + split·dir``. Returns the direction MAXIMISING the distance from that
    daughter site to the nearest existing neighbour, together with that distance (the
    "free-space gap"). A buried bulk cell has a small gap (no room); a rim cell a large one.
    With no neighbours the gap is ``+inf`` (anywhere is free) and the first candidate is used.

    Args:
        positions: (N, 3) all cell centers (m).
        i: index of the dividing cell.
        neighbor_idx: indices of ``i``'s nearby cells (the local set to clear; self excluded).
        split: daughter offset from the parent (m), = rest separation r0.
        candidates: optional precomputed candidate unit vectors (defaults to a 12-pt sphere).

    Returns:
        ``(nhat, gap)`` — chosen unit direction and the daughter's nearest-neighbour distance.
    """
    dirs = candidate_directions() if candidates is None else candidates
    if neighbor_idx.size == 0:
        return dirs[0], float("inf")
    cand_sites = positions[i] + split * dirs                 # (K, 3)
    nb = positions[neighbor_idx]                             # (M, 3)
    dmin = cdist(cand_sites, nb).min(axis=1)                 # (K,)
    k = int(np.argmax(dmin))
    return dirs[k], float(dmin[k])


def sample_cycle_targets(
    n: int, prolif: ResolvedProliferation, rng: np.random.Generator
) -> npt.NDArray[np.float64]:
    """Sample ``n`` per-cell cycle-completion times (Gaussian: mean·(1+CV·N(0,1)), floored)."""
    t = prolif.cycle_time_mean * (1.0 + prolif.cycle_time_cv * rng.standard_normal(n))
    # Floor at 10% of the mean so a rare large negative draw cannot give a non-positive target.
    return np.maximum(t, 0.1 * prolif.cycle_time_mean)


def apply_divisions(
    positions: npt.NDArray[np.float64],
    ages: npt.NDArray[np.float64],
    targets: npt.NDArray[np.float64],
    prolif: ResolvedProliferation,
    rng: np.random.Generator,
) -> dict[str, Any]:
    """Apply all due, non-inhibited divisions to the population (one pass).

    A timer-elapsed cell (``age >= target``) divides iff it passes the contact-inhibition
    gate: (PRIMARY) the Drasdo-Höhme FREE-SPACE rule — its freest bud direction leaves the
    daughter at least ``prolif.min_gap`` (the Morse repulsive-core onset) from every existing
    cell, i.e. there is room; AND (SECONDARY) its first shell is below the kissing-number
    ceiling. Each dividing cell becomes two: daughter 1 inherits the parent site, daughter 2
    buds at the rest separation into that freest direction. Both daughters reset age 0 and
    draw fresh cycle targets. Quiescent or not-yet-due cells pass through unchanged.

    The free-space gate (not a neighbour count) is what makes division RIM-localised: the
    settled liquid-like packing has bulk free-gap ≈ 0.7·r0 < min_gap (no room) but rim
    free-gap ≈ 1.3·r0 (room) — see the L2.4 REPORT diagnostic.

    Returns a dict with the new ``positions``/``ages``/``targets`` arrays, the parent indices
    that divided (``divided_idx``), their pre-division radial distances from the population
    COM (``divided_radial``), and the population ``median_radial`` (rim-localisation gate).
    """
    p = np.asarray(positions, dtype=np.float64)
    n = p.shape[0]
    counts = first_shell_counts(p, prolif.shell_cutoff)
    com = p.mean(axis=0)
    radial_all = np.linalg.norm(p - com, axis=1)

    due = np.where((ages >= targets) & (counts < prolif.kissing_number))[0]
    if due.size == 0:
        return {
            "positions": p, "ages": ages, "targets": targets,
            "divided_idx": due, "divided_radial": radial_all[due],
            "median_radial": float(np.median(radial_all)) if n else 0.0,
        }

    cands = candidate_directions()
    tree = cKDTree(p)
    # search radius for the free-space test: a daughter sits at `split`; clearing min_gap
    # around it means any neighbour within split + min_gap could matter.
    search_r = prolif.split_distance + prolif.min_gap
    keep = np.ones(n, dtype=bool)
    divided, new_pos, new_ages, new_targets = [], [], [], []
    for i in due:
        nbr = np.array(
            [j for j in tree.query_ball_point(p[i], search_r) if j != i], dtype=np.int64
        )
        nhat, gap = best_bud_direction(p, int(i), nbr, prolif.split_distance, candidates=cands)
        if gap < prolif.min_gap:
            continue  # no room — contact-inhibited (quiescent); timer keeps running
        new_pos.append(p[i])                                   # daughter 1: parent site
        new_pos.append(p[i] + prolif.split_distance * nhat)    # daughter 2: budded
        new_ages.extend([0.0, 0.0])
        new_targets.extend(sample_cycle_targets(2, prolif, rng).tolist())
        keep[i] = False
        divided.append(int(i))

    divided_idx = np.array(divided, dtype=np.int64)
    if divided_idx.size == 0:
        return {
            "positions": p, "ages": ages, "targets": targets,
            "divided_idx": divided_idx, "divided_radial": radial_all[divided_idx],
            "median_radial": float(np.median(radial_all)),
        }

    pos2 = np.vstack([p[keep], np.array(new_pos, dtype=np.float64)])
    ages2 = np.concatenate([ages[keep], np.array(new_ages, dtype=np.float64)])
    targets2 = np.concatenate([targets[keep], np.array(new_targets, dtype=np.float64)])
    return {
        "positions": pos2, "ages": ages2, "targets": targets2,
        "divided_idx": divided_idx, "divided_radial": radial_all[divided_idx],
        "median_radial": float(np.median(radial_all)),
    }


def run_growth(
    resolved: ResolvedL2,
    prolif: ResolvedProliferation,
    n_cells_init: int = 120,
    *,
    total_time: float,
    epoch_steps: int = 2_000,
    settle_steps: int = 1_500,
    max_cells: int = 4_000,
    device: hoomd.device.Device | None = None,
    seed: int | None = None,
    f_traction: float = 0.0,
    Lp: float = 11.0e-6,
) -> dict[str, Any]:
    """Grow a CBM spheroid by contact-inhibited proliferation; record A(t), N(t).

    Settles a loose blob (G1), then advances the population in epochs: relax → read back →
    divide → rebuild. Optionally superposes edge-directed active-wetting traction
    (``f_traction`` > 0, recomputed per epoch) so growth and motility act together.

    Args:
        resolved: resolved CBM parameters.
        prolif: resolved proliferation parameters.
        n_cells_init: starting cell count.
        total_time: biological time to simulate (s).
        epoch_steps: BAOAB steps between division checks (≈ epoch_steps·dt seconds/epoch).
        settle_steps: pre-growth settle (defines A0).
        max_cells: hard cap (stops the run; logged, not silently truncated).
        device: HOOMD device (CPU default).
        seed: realization key.
        f_traction: optional edge-traction magnitude (N); 0 = proliferation only.
        Lp: edge-traction screening length (m) when ``f_traction`` > 0.

    Returns:
        Dict with time/count/area/Rg series, A0, A/A0 series, rim-localisation diagnostics,
        and init/settled/final positions.
    """
    seed = resolved.seed if seed is None else int(seed)
    rng = np.random.default_rng(seed)
    dt = resolved.dt_cfl

    # --- initial blob + settle (defines A0) ---
    sim, _a, _u, _rc = build_cbm_simulation(resolved, n_cells_init, device=device, seed=seed)
    sim.run(0)
    pos_init = get_positions(sim)
    sim.run(settle_steps)
    pos = get_positions(sim)
    a0 = projected_area(pos)

    ages = np.zeros(pos.shape[0], dtype=np.float64)
    targets = sample_cycle_targets(pos.shape[0], prolif, rng)
    # Desynchronise: start cells at a random phase of their cycle (else all divide at once).
    ages = rng.uniform(0.0, 1.0, pos.shape[0]) * targets

    t = 0.0
    epoch = 0
    ts, ns, areas, rgs = [0.0], [pos.shape[0]], [a0], [radius_of_gyration(pos)]
    rim_frac_mean = []  # fraction of divisions that were above-median-radial (rim) per epoch
    capped = False

    while t < total_time:
        if pos.shape[0] >= max_cells:
            capped = True
            break
        epoch += 1
        # rebuild a fresh state from the current (possibly grown) population
        sim, _a, _u, _rc = build_cbm_simulation(
            resolved, pos.shape[0], device=device, seed=seed + epoch, positions=pos
        )
        if f_traction > 0.0:
            from ffn_sim.spheroid.spreading import SettableForce, edge_outward_forces

            sim.operations.tuners.clear()
            edge = SettableForce(pos.shape[0])
            edge.set_vectors(edge_outward_forces(pos, f_traction=f_traction, Lp=Lp))
            sim.operations.integrator.forces.append(edge)
        sim.run(0)
        sim.run(epoch_steps)
        pos = get_positions(sim)

        dt_epoch = epoch_steps * dt
        ages = ages + dt_epoch
        t += dt_epoch

        div = apply_divisions(pos, ages, targets, prolif, rng)
        pos, ages, targets = div["positions"], div["ages"], div["targets"]
        if div["divided_idx"].size > 0:
            rim = float(np.mean(div["divided_radial"] >= div["median_radial"]))
            rim_frac_mean.append(rim)

        ts.append(t); ns.append(pos.shape[0])
        areas.append(projected_area(pos)); rgs.append(radius_of_gyration(pos))

    areas = np.asarray(areas)
    return {
        "n_cells_init": n_cells_init,
        "a0": a0,
        "t": np.asarray(ts),
        "n_cells": np.asarray(ns),
        "area": areas,
        "area_over_a0": areas / a0 if a0 > 0 else areas,
        "rg": np.asarray(rgs),
        "rim_fraction_mean": float(np.mean(rim_frac_mean)) if rim_frac_mean else float("nan"),
        "n_division_epochs": len(rim_frac_mean),
        "growth_factor": ns[-1] / ns[0] if ns[0] else float("nan"),
        "capped_at_max_cells": capped,
        "f_traction": f_traction,
        "pos_init": pos_init,
        "pos_final": pos,
    }
