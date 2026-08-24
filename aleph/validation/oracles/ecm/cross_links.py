"""Mikado cross-link generation by 2D segment intersection (KU-1.27, KU-1.3, KU-1.28).

This is the *scientific* construction: for every pair of fibers we
solve the line-segment intersection in 2D; for each crossing we anchor
a cross-link to the nearest bead on each fiber. The resulting
coordination ⟨z⟩ is an emergent property of the network — not a
construction target.

Periodic boundaries
-------------------
Phase 1 enumerates intersections of fiber pairs in the **primary copy
only**: a fiber whose segment runs through a box edge has its periodic
images not checked against neighbours' images. With
L_fiber / L_box ≈ 0.05 (Phase 1 default), ~5 % of fibers cross a
boundary and a fraction of their potential crossings are missed; the
emergent ⟨z⟩ is therefore slightly below the infinite-Mikado theory
prediction. The Sanity Gate tracks both the prediction and the
measurement so the gap is explicit. Phase 2+ will sweep periodic
images (9 in 2D) to recover the missing pairs.

Cross-link force / energy
-------------------------
Each link is a harmonic spring (KU-1.28) between two beads:

    E_xl = ½ k_xl (|r_b − r_a| − r_0)²
    F_a  = +k_xl (|d| − r_0) · d̂,    F_b = − F_a,    d = MI(r_b − r_a)

The energy and forces are computed in
:func:`compute_xl_energy_and_forces` and accumulated into the total
mechanics by :func:`aleph.validation.oracles.ecm.fiber_mechanics.compute_forces` when
the optional ``cross_links`` argument is supplied.

Sanity Gate
-----------
1. Dimensional analysis: parametric s,t ∈ [0,1] dimensionless; rest
   length in m; k_xl in N/m; E_xl in J; F in N. No CFL (no time
   integration in this module).
2. Boundary cases: parallel segments (det = 0) → skipped; |d| = 0 in
   force kernel → guarded.
3. Conservation: per-link force is equal and opposite by construction.
4. Numerical sanity: float64; intersection test vectorized; AABB
   pre-filter rejects O(N²) non-overlapping pairs.
5. Sign / sense: stretched (|d| > r_0) link pulls beads together,
   compressed link pushes apart.
6. Measurement protocol: ``measure_coordination`` reports the achieved
   ⟨z⟩ averaged over every bead, matching the KU-1.3 reporting
   convention; the Sanity Gate compares to the KU-1.3 sub-isostatic
   range, not to a constructed target.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aleph.validation.oracles.ecm.fiber_network import FiberNetwork


@dataclass(frozen=True, slots=True)
class CrossLink:
    fiber_a: int
    bead_a: int
    fiber_b: int
    bead_b: int
    rest_length: float   # m, set to the bead-bead MI distance at construction
    stiffness: float     # N/m, KU-1.28


def _backbone_z_mean(beads_per_fiber: int) -> float:
    if beads_per_fiber < 2:
        return 0.0
    return (2.0 * (beads_per_fiber - 2) + 2.0) / beads_per_fiber


def _segment_intersections(
    endpoints: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (i, j, s_a, s_b) for every interior intersection.

    Parameters
    ----------
    endpoints : np.ndarray, shape (N, 2, 2)
        Unwrapped (p1, p2) for each fiber.

    Returns
    -------
    pair_i, pair_j : np.ndarray, shape (M,)
        Indices with i < j.
    s_a, s_b : np.ndarray, shape (M,)
        Parametric positions of the intersection along fibers i and j,
        in [-L_f/2, +L_f/2]. Convert to bead indices upstream.
    """
    N = endpoints.shape[0]
    p1 = endpoints[:, 0, :]
    p2 = endpoints[:, 1, :]
    d = p2 - p1                                          # (N, 2)
    L = np.linalg.norm(d, axis=1)                        # all equal = L_fiber
    half_L = 0.5 * L[0]

    # AABB pre-filter (vectorized).
    lo = np.minimum(p1, p2)                              # (N, 2)
    hi = np.maximum(p1, p2)
    # Indices of upper triangle pairs i<j.
    iu, ju = np.triu_indices(N, k=1)
    overlap = (
        (lo[iu, 0] <= hi[ju, 0]) & (lo[ju, 0] <= hi[iu, 0]) &
        (lo[iu, 1] <= hi[ju, 1]) & (lo[ju, 1] <= hi[iu, 1])
    )
    iu = iu[overlap]
    ju = ju[overlap]
    if iu.size == 0:
        return (np.array([], dtype=np.int64),) * 4

    # Solve t (along i) and u (along j) by Cramer's rule.
    da = d[iu]                                           # (M, 2)
    db = d[ju]
    diff = p1[ju] - p1[iu]                               # (M, 2)
    # det = da_x · (-db_y) − da_y · (-db_x) = -da_x db_y + da_y db_x
    det = da[:, 1] * db[:, 0] - da[:, 0] * db[:, 1]
    nonparallel = np.abs(det) > 1e-30
    iu, ju, da, db, diff, det = (
        iu[nonparallel], ju[nonparallel], da[nonparallel],
        db[nonparallel], diff[nonparallel], det[nonparallel],
    )
    # t = (diff_x · (-db_y) − diff_y · (-db_x)) / det
    #   = (-diff_x db_y + diff_y db_x) / det
    t = (-diff[:, 0] * db[:, 1] + diff[:, 1] * db[:, 0]) / det
    u = (-diff[:, 0] * da[:, 1] + diff[:, 1] * da[:, 0]) / det
    inside = (t > 0.0) & (t < 1.0) & (u > 0.0) & (u < 1.0)
    iu, ju, t, u = iu[inside], ju[inside], t[inside], u[inside]

    # Map t ∈ [0,1] to s_a ∈ [-L/2, +L/2]
    s_a = (t - 0.5) * (2.0 * half_L)
    s_b = (u - 0.5) * (2.0 * half_L)
    return iu.astype(np.int64), ju.astype(np.int64), s_a, s_b


def generate_cross_links(
    net: FiberNetwork,
    stiffness: float = 1.0e-3,
) -> list[CrossLink]:
    """Generate cross-links by 2D segment intersection (Mikado, KU-1.27).

    For each fiber pair that geometrically crosses, link the *nearest*
    bead on each fiber to the intersection point (Phase 1 discrete
    approximation). The cross-link rest length is set to the current
    bead–bead distance under the minimum-image convention so the link
    starts force-free.
    """
    iu, ju, s_a, s_b = _segment_intersections(net.fiber_endpoints)
    if iu.size == 0:
        return []

    beads = net.bead_positions.shape[1]
    L_f = net.fiber_length
    # Bead index nearest to a parametric position s ∈ [-L_f/2, +L_f/2].
    bead_a = np.clip(np.round((s_a + 0.5 * L_f) * (beads - 1) / L_f),
                     0, beads - 1).astype(np.int64)
    bead_b = np.clip(np.round((s_b + 0.5 * L_f) * (beads - 1) / L_f),
                     0, beads - 1).astype(np.int64)

    ra = net.bead_positions[iu, bead_a, :]
    rb = net.bead_positions[ju, bead_b, :]
    delta = rb - ra
    delta -= net.box_size * np.round(delta / net.box_size)
    rest = np.linalg.norm(delta, axis=1)

    return [
        CrossLink(
            fiber_a=int(iu[k]), bead_a=int(bead_a[k]),
            fiber_b=int(ju[k]), bead_b=int(bead_b[k]),
            rest_length=float(rest[k]),
            stiffness=float(stiffness),
        )
        for k in range(iu.size)
    ]


def measure_coordination(
    links: list[CrossLink],
    n_beads_total: int,
    beads_per_fiber: int,
) -> float:
    """Total per-bead coordination ⟨z⟩ = backbone_z + 2 N_links/N_beads."""
    if n_beads_total <= 0:
        return 0.0
    return _backbone_z_mean(beads_per_fiber) + 2.0 * len(links) / n_beads_total


def measure_xl_per_fiber(links: list[CrossLink], n_fibers: int) -> float:
    """Mean number of cross-links attached to a fiber (≈ Mikado crossings)."""
    if n_fibers <= 0:
        return 0.0
    return 2.0 * len(links) / n_fibers


# ------------------------------------------------------------------ #
# Energy / force kernel                                              #
# ------------------------------------------------------------------ #

def compute_xl_energy_and_forces(
    bead_positions: np.ndarray,
    cross_links: list[CrossLink],
    box_size: float,
) -> tuple[float, np.ndarray]:
    """Harmonic cross-link energy and per-bead forces (KU-1.28).

    Returns
    -------
    energy : float
        Σ ½ k_xl (|d| − r_0)²  in joules.
    forces : np.ndarray, shape (n_fibers, n_beads, 2)
        Per-bead force contribution from all cross-links.
    """
    forces = np.zeros_like(bead_positions)
    if not cross_links:
        return 0.0, forces

    fa = np.fromiter((xl.fiber_a for xl in cross_links), dtype=np.int64,
                     count=len(cross_links))
    ba = np.fromiter((xl.bead_a for xl in cross_links), dtype=np.int64,
                     count=len(cross_links))
    fb = np.fromiter((xl.fiber_b for xl in cross_links), dtype=np.int64,
                     count=len(cross_links))
    bb = np.fromiter((xl.bead_b for xl in cross_links), dtype=np.int64,
                     count=len(cross_links))
    r0 = np.fromiter((xl.rest_length for xl in cross_links), dtype=np.float64,
                     count=len(cross_links))
    k = np.fromiter((xl.stiffness for xl in cross_links), dtype=np.float64,
                    count=len(cross_links))

    ra = bead_positions[fa, ba, :]
    rb = bead_positions[fb, bb, :]
    d = rb - ra
    d -= box_size * np.round(d / box_size)
    dist = np.linalg.norm(d, axis=1)
    safe = np.where(dist > 0.0, dist, 1.0)
    extension = dist - r0
    energy = 0.5 * float(np.sum(k * extension * extension))

    fmag = (k * extension) / safe                # (M,)
    # F_a = +k(|d|−r0) d/|d|; F_b = −F_a (Newton 3rd)
    fvec = fmag[:, None] * d
    np.add.at(forces, (fa, ba), fvec)
    np.add.at(forces, (fb, bb), -fvec)
    return energy, forces
