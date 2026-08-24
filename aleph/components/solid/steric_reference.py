"""Host-side reference for the filament-filament excluded-volume force (the I2b algorithm oracle).

The Warp-CUDA runtime (``steric_warp.StericForce``) accumulates the WCA steric force over a device
``wp.HashGrid`` neighbour search. That hash grid is a pure ACCELERATOR: it prunes only pairs beyond the
cutoff ``r_c`` that contribute exactly zero, so the accelerated result must equal the brute-force O(N^2)
sum to summation-order tolerance (the same contract ``dcm/dcm_neighbor_warp.py`` documents for its grid
kernels). This module provides BOTH references in pure NumPy — no Warp, no CUDA — so:

  * :func:`steric_force_bruteforce` is the O(N^2) ground truth (every ordered pair evaluated);
  * :func:`steric_force_celllist` is a CPU cell-list — the exact algorithmic analogue of the device
    ``wp.HashGrid`` (bin into cells of size >= r_c, query the 27-cell stencil) — and MUST reproduce the
    brute-force force bit-for-bit-to-tolerance. That equivalence is the Magic-Number-Block grid-invariance
    statement for a neighbour-accelerated pair potential: the neighbour cell size is a numerical accelerator
    knob that does not change the physics. The lead re-certifies the SAME equivalence natively (Warp
    HashGrid vs brute) on the gbook A5000 (see INTEGRATION.md native-gate spec).

Force convention (matches ``dcm_neighbor_warp.cohesion_grid_kernel``): one row per node; node i sums the
repulsion from every non-excluded neighbour j, ``F_i += F_WCA(|r_i-r_j|) (r_i-r_j)/|r_i-r_j|``. Newton's
third law is satisfied automatically because node j's row adds the equal-and-opposite term, so the net
internal force is zero to round-off (momentum-clean). Same-filament pairs are excluded (a filament's own
segment spacing is held by the bending/inextensibility laws, not by steric self-repulsion).

Units: engine runtime um-pN-s (positions um, forces pN, stiffness pN/um). Unit-agnostic in the math.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.components.solid.wca_analytic import (
    contact_stiffness,
    max_core_stiffness,
    wca_cutoff,
    wca_force_magnitude,
)

__all__ = [
    "steric_force_bruteforce",
    "steric_force_celllist",
    "total_steric_energy",
    "cfl_stiffness_term",
]


def _validate(pos: npt.NDArray[np.float64], fiber_id: npt.NDArray[np.int64] | None) -> None:
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError("pos must have shape (N, 3)")
    if fiber_id is not None and fiber_id.shape != (pos.shape[0],):
        raise ValueError("fiber_id must have shape (N,)")


def steric_force_bruteforce(
    pos: npt.ArrayLike,
    sigma: float,
    epsilon: float,
    fiber_id: npt.ArrayLike | None = None,
    *,
    exclude_same_fiber: bool = True,
) -> npt.NDArray[np.float64]:
    """Brute-force O(N^2) WCA steric force on every node — the ground truth for the hash-grid kernel.

    Args:
        pos: Node positions, shape (N, 3) [length].
        sigma: Steric diameter [length]; cutoff ``r_c = 2^(1/6) sigma``.
        epsilon: WCA energy scale [energy].
        fiber_id: Per-node filament id, shape (N,); required if ``exclude_same_fiber``.
        exclude_same_fiber: Skip pairs on the same filament (bonded segments do not sterically self-repel).

    Returns:
        Force on each node, shape (N, 3) [force]. Sums to ~0 (Newton's third law).
    """
    pos = np.asarray(pos, dtype=np.float64)
    fid = None if fiber_id is None else np.asarray(fiber_id, dtype=np.int64)
    _validate(pos, fid)
    if exclude_same_fiber and fid is None:
        raise ValueError("exclude_same_fiber requires fiber_id")
    n = pos.shape[0]
    r_c = wca_cutoff(sigma)
    force = np.zeros((n, 3), dtype=np.float64)
    for i in range(n):
        d = pos[i] - pos                                   # (N,3): r_i - r_j
        dist = np.linalg.norm(d, axis=1)
        mask = (dist > 0.0) & (dist < r_c)
        mask[i] = False
        if exclude_same_fiber:
            mask &= fid != fid[i]
        if not np.any(mask):
            continue
        rr = dist[mask]
        fmag = wca_force_magnitude(rr, sigma, epsilon)
        force[i] = np.sum((fmag / rr)[:, None] * d[mask], axis=0)
    return force


def steric_force_celllist(
    pos: npt.ArrayLike,
    sigma: float,
    epsilon: float,
    fiber_id: npt.ArrayLike | None = None,
    *,
    exclude_same_fiber: bool = True,
    cell_size: float | None = None,
) -> npt.NDArray[np.float64]:
    """Cell-list WCA steric force — the CPU analogue of the device ``wp.HashGrid`` neighbour search.

    Bins nodes into a uniform grid of cells of side ``cell_size`` (>= ``r_c``); each node queries only the
    3x3x3 stencil around its own cell. Because every within-cutoff pair lies in that stencil, the result
    equals :func:`steric_force_bruteforce` to summation-order tolerance for ANY ``cell_size >= r_c`` — the
    grid-invariance property the native Warp hash-grid gate re-certifies.

    Args:
        pos: Node positions, shape (N, 3) [length].
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].
        fiber_id: Per-node filament id, shape (N,); required if ``exclude_same_fiber``.
        exclude_same_fiber: Skip same-filament pairs.
        cell_size: Neighbour cell side [length]; defaults to ``r_c``. Must be ``>= r_c``.

    Returns:
        Force on each node, shape (N, 3) [force] — equal to the brute-force force to tolerance.
    """
    pos = np.asarray(pos, dtype=np.float64)
    fid = None if fiber_id is None else np.asarray(fiber_id, dtype=np.int64)
    _validate(pos, fid)
    if exclude_same_fiber and fid is None:
        raise ValueError("exclude_same_fiber requires fiber_id")
    n = pos.shape[0]
    r_c = wca_cutoff(sigma)
    h = r_c if cell_size is None else float(cell_size)
    if h < r_c * (1.0 - 1e-12):
        raise ValueError("cell_size must be >= r_c (else within-cutoff pairs are missed)")

    origin = pos.min(axis=0)
    cell = np.floor((pos - origin) / h).astype(np.int64)          # (N,3) integer cell coords
    # bucket node indices by cell key
    buckets: dict[tuple[int, int, int], list[int]] = {}
    for idx in range(n):
        buckets.setdefault(tuple(cell[idx]), []).append(idx)

    force = np.zeros((n, 3), dtype=np.float64)
    stencil = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    for i in range(n):
        ci = cell[i]
        neigh: list[int] = []
        for dx, dy, dz in stencil:
            key = (int(ci[0] + dx), int(ci[1] + dy), int(ci[2] + dz))
            neigh.extend(buckets.get(key, ()))
        if not neigh:
            continue
        j = np.array(neigh, dtype=np.int64)
        d = pos[i] - pos[j]
        dist = np.linalg.norm(d, axis=1)
        mask = (dist > 0.0) & (dist < r_c) & (j != i)
        if exclude_same_fiber:
            mask &= fid[j] != fid[i]
        if not np.any(mask):
            continue
        rr = dist[mask]
        fmag = wca_force_magnitude(rr, sigma, epsilon)
        force[i] = np.sum((fmag / rr)[:, None] * d[mask], axis=0)
    return force


def total_steric_energy(
    pos: npt.ArrayLike,
    sigma: float,
    epsilon: float,
    fiber_id: npt.ArrayLike | None = None,
    *,
    exclude_same_fiber: bool = True,
) -> float:
    """Total WCA excluded-volume energy of a node configuration (each unordered pair counted once).

    Args:
        pos: Node positions, shape (N, 3).
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].
        fiber_id: Per-node filament id, shape (N,).
        exclude_same_fiber: Skip same-filament pairs.

    Returns:
        Total steric potential energy [energy] (>= 0; exactly 0 when no pair overlaps).
    """
    from aleph.components.solid.wca_analytic import wca_energy

    pos = np.asarray(pos, dtype=np.float64)
    fid = None if fiber_id is None else np.asarray(fiber_id, dtype=np.int64)
    _validate(pos, fid)
    n = pos.shape[0]
    r_c = wca_cutoff(sigma)
    total = 0.0
    for i in range(n):
        j = np.arange(i + 1, n)
        if j.size == 0:
            continue
        d = pos[i] - pos[j]
        dist = np.linalg.norm(d, axis=1)
        mask = (dist > 0.0) & (dist < r_c)
        if exclude_same_fiber:
            mask &= fid[j] != fid[i]
        if np.any(mask):
            total += float(np.sum(wca_energy(dist[mask], sigma, epsilon)))
    return total


def cfl_stiffness_term(
    sigma: float,
    epsilon: float,
    f_cap: float | None = None,
) -> float:
    """The steric stiffness ``k_EV`` folded into the CFL ``kmax`` (``dt_mu = 0.1 / kmax``).

    Without a load cap this is the contact stiffness ``U''(r_min)`` (the well curvature). With an operating
    force cap ``f_cap`` (the max per-node force the mechanical solve applies), it is the deeper-overlap
    stiffness ``U''(r_eq(f_cap))`` — the stiffness the integrator actually sees, never the divergent core
    at ``r -> 0``. Grid-invariant: a function of (sigma, epsilon, f_cap) only, independent of the field
    grid ``dx`` and the neighbour cell size.

    Args:
        sigma: Steric diameter [length].
        epsilon: WCA energy scale [energy].
        f_cap: Optional operating-load cap [force]; if given, uses the deeper-overlap stiffness.

    Returns:
        Stiffness [force/length] to include as ``kmax = max(kmax, k_EV)``.
    """
    if f_cap is None:
        return contact_stiffness(epsilon, sigma)
    return max_core_stiffness(f_cap, sigma, epsilon)
