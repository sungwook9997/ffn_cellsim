"""Point-cloud observables for the Layer-2 CBM spheroid (cell-center particle positions).

Measurement protocol: positions (N x 3, metres) -> spread area, effective radius, radius
of gyration, radial density profile, nearest-neighbour statistics, detached ("gas-like")
fraction. Pure numpy/scipy; NO HOOMD, NO physics parameters, NO baked magic numbers — any
threshold (e.g. ``d_crit``) is an explicit caller argument. Written BEFORE the physics run
per the Sanity-Gate measurement-protocol rule (CLAUDE.md). The closed-form ``A/A0`` law and
the analytic cross-checks live in the runtime-import-forbidden oracle
``ffn_sim.validation.oracles.spheroid.aa0_law``.

These observables drive the Layer-2 acceptance gates:
- G1 (stable aggregate): nearest-neighbour median ~ rest separation; bounded gyration
  radius; ~zero detached fraction.
- G3 (A/A0 sweep): ``projected_area`` over an initial-radius sweep -> the spreading curve.

Sanity Gate
-----------
- Dimensional: positions in metres -> areas [m^2], radii/distances [m], density [m^-3].
- Boundary: < 3 (non-degenerate) points -> projected_area returns 0.0 (no hull).
- Reference limits (validated in tests against the oracle): points filling a disk of
  radius R -> projected_area -> pi R^2; a solid ball radius R -> Rg -> sqrt(3/5) R; a
  cubic lattice of spacing s -> nearest-neighbour distance = s.
- COM-invariance: gyration radius and detached fraction are translation-invariant.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.spatial import ConvexHull, QhullError, cKDTree

__all__ = [
    "center_of_mass",
    "projected_area",
    "effective_radius",
    "radius_of_gyration",
    "radial_density_profile",
    "nearest_neighbor_distances",
    "nearest_neighbor_stats",
    "detached_fraction",
    "connected_components",
    "largest_connected_component",
    "core_projected_area",
    "convex_hull_volume",
    "virial_pressure",
]


def _as_positions(positions: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Validate and coerce to an (N, 3) float array."""
    p = np.asarray(positions, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError(f"positions must be (N, 3); got shape {p.shape}.")
    return p


def center_of_mass(positions: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Unweighted centroid of the cell-center positions (metres)."""
    p = _as_positions(positions)
    return p.mean(axis=0)


def projected_area(
    positions: npt.ArrayLike, plane: tuple[int, int] = (0, 1)
) -> float:
    """2D convex-hull area of the cell centers projected onto two axes.

    This is the platform analog of the experiment's segmented spread area A (a 2D
    footprint). Convex hull is the MVP estimator; a concave (alpha-shape) estimator is a
    documented future refinement for strongly non-convex spreading fronts.

    Args:
        positions: (N, 3) cell-center positions (metres).
        plane: the two axis indices to project onto (default x, y).

    Returns:
        Hull area in m^2; ``0.0`` if fewer than 3 points or all points are collinear.
    """
    p = _as_positions(positions)
    pts2d = p[:, list(plane)]
    if pts2d.shape[0] < 3:
        return 0.0
    try:
        hull = ConvexHull(pts2d)
    except QhullError:
        return 0.0  # collinear / degenerate -> no area
    return float(hull.volume)  # scipy: 2D ConvexHull.volume == enclosed area


def effective_radius(area: float) -> float:
    """Effective radius of a given spread area: ``R = sqrt(A / pi)`` (metres)."""
    if area < 0.0:
        raise ValueError("area must be non-negative.")
    return float(np.sqrt(area / np.pi))


def radius_of_gyration(positions: npt.ArrayLike) -> float:
    """Root-mean-square distance of cell centers from their centroid (metres)."""
    p = _as_positions(positions)
    d = p - p.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(d * d, axis=1))))


def radial_density_profile(
    positions: npt.ArrayLike,
    n_bins: int = 20,
    center: npt.ArrayLike | None = None,
    r_max: float | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """3D radial number-density profile about a center (for spheroid-radius inflection).

    Args:
        positions: (N, 3) positions (metres).
        n_bins: number of radial shells.
        center: reference center (metres); defaults to the centroid.
        r_max: outer radius (metres); defaults to the max particle distance.

    Returns:
        ``(r_centers, density)`` where density = counts / spherical-shell-volume [m^-3].
    """
    p = _as_positions(positions)
    c = p.mean(axis=0) if center is None else np.asarray(center, dtype=np.float64)
    r = np.linalg.norm(p - c, axis=1)
    rmax = float(r.max()) if r_max is None else float(r_max)
    if rmax <= 0.0:
        raise ValueError("r_max must be positive (all points coincide with center?).")
    edges = np.linspace(0.0, rmax, n_bins + 1)
    counts, _ = np.histogram(r, bins=edges)
    shell_vol = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    density = counts / shell_vol
    r_centers = 0.5 * (edges[1:] + edges[:-1])
    return r_centers, density


def nearest_neighbor_distances(positions: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Per-particle distance to the nearest other particle (metres), via a KD-tree."""
    p = _as_positions(positions)
    if p.shape[0] < 2:
        raise ValueError("need >= 2 particles for nearest-neighbour distances.")
    tree = cKDTree(p)
    dists, _ = tree.query(p, k=2)  # column 0 is self (distance 0)
    return dists[:, 1]


def nearest_neighbor_stats(positions: npt.ArrayLike) -> dict[str, float]:
    """Summary statistics of the nearest-neighbour distance distribution (metres)."""
    nn = nearest_neighbor_distances(positions)
    return {
        "mean": float(nn.mean()),
        "median": float(np.median(nn)),
        "std": float(nn.std()),
        "min": float(nn.min()),
        "max": float(nn.max()),
    }


def detached_fraction(
    positions: npt.ArrayLike,
    d_crit: float,
    neighbor_radius: float,
    center: npt.ArrayLike | None = None,
) -> float:
    """Fraction of "gas-like" detached cells (Herold 2023 point-cloud criterion).

    A cell counts as detached if it is BOTH farther than ``d_crit`` from the aggregate
    center AND has no neighbour within ``neighbor_radius``. Both thresholds are explicit
    caller arguments (no baked constant); callers typically set them as multiples of the
    rest separation / aggregate radius so the metric is grid- and scale-invariant.

    Args:
        positions: (N, 3) positions (metres).
        d_crit: distance-from-center threshold (metres, > 0).
        neighbor_radius: isolation radius (metres, > 0).
        center: reference center; defaults to the centroid.

    Returns:
        Detached fraction in [0, 1].
    """
    p = _as_positions(positions)
    if d_crit <= 0.0 or neighbor_radius <= 0.0:
        raise ValueError("d_crit and neighbor_radius must be strictly positive.")
    c = p.mean(axis=0) if center is None else np.asarray(center, dtype=np.float64)
    far = np.linalg.norm(p - c, axis=1) > d_crit
    tree = cKDTree(p)
    # neighbours within radius incl. self -> count >= 2 means a real neighbour exists.
    n_within = np.array([len(tree.query_ball_point(pt, neighbor_radius)) for pt in p])
    isolated = n_within < 2
    detached = far & isolated
    return float(detached.sum() / p.shape[0])


def connected_components(
    positions: npt.ArrayLike, link_radius: float
) -> npt.NDArray[np.int64]:
    """Label connected components: cells within ``link_radius`` are in the same cluster.

    Single-linkage clustering via a KD-tree + union-find. Two cells are "linked" (same
    aggregate) if their centre separation is <= ``link_radius`` — a caller argument (no baked
    constant), typically set just past the rest separation (e.g. 1.5-1.6 r0, between the 1st
    and 2nd coordination shell) so touching/cohesive cells link but a detached fragment does
    not. This separates a fragmented population into its drifting pieces — the structure the
    raw convex hull conflates.

    Args:
        positions: (N, 3) cell centers (m).
        link_radius: maximum centre separation for two cells to share a cluster (m, > 0).

    Returns:
        (N,) integer component labels in ``[0, n_components)``, ordered by descending size
        (label 0 is the largest component).
    """
    p = _as_positions(positions)
    n = p.shape[0]
    if link_radius <= 0.0:
        raise ValueError("link_radius must be strictly positive.")
    if n == 0:
        return np.zeros(0, dtype=np.int64)

    parent = np.arange(n)

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]  # path halving
            a = parent[a]
        return int(a)

    tree = cKDTree(p)
    for i, j in tree.query_pairs(link_radius):  # unique unordered pairs within link_radius
        ri, rj = find(int(i)), find(int(j))
        if ri != rj:
            parent[ri] = rj

    roots = np.array([find(i) for i in range(n)])
    # relabel by descending component size so label 0 is the largest cluster
    uniq, counts = np.unique(roots, return_counts=True)
    order = uniq[np.argsort(-counts)]
    remap = {int(r): k for k, r in enumerate(order)}
    return np.array([remap[int(r)] for r in roots], dtype=np.int64)


def largest_connected_component(
    positions: npt.ArrayLike, link_radius: float
) -> npt.NDArray[np.bool_]:
    """Boolean mask of the cells in the largest connected component (the spheroid core)."""
    labels = connected_components(positions, link_radius)
    if labels.size == 0:
        return np.zeros(0, dtype=bool)
    return labels == 0  # label 0 is the largest (connected_components orders by size)


def convex_hull_volume(positions: npt.ArrayLike) -> float:
    """3D convex-hull volume of the cell centers [m^3] (0.0 if < 4 non-coplanar points).

    The aggregate volume used to normalise the virial pressure. For a spheroid this is the
    natural enclosing volume; for a substrate-confined 3D cap it is the cap's hull volume.
    """
    p = _as_positions(positions)
    if p.shape[0] < 4:
        return 0.0
    try:
        hull = ConvexHull(p)
    except QhullError:
        return 0.0  # coplanar / degenerate -> no volume
    return float(hull.volume)  # scipy: 3D ConvexHull.volume == enclosed volume


def virial_pressure(
    positions: npt.ArrayLike,
    pair_indices: npt.ArrayLike,
    pair_forces: npt.ArrayLike,
    volume: float | None = None,
    *,
    kT: float = 0.0,
) -> float:
    """Mechanical (virial) pressure of the cell aggregate [Pa] — the emergent-tension route.

    Configurational virial pressure of a set of cell centers under pairwise interactions::

        P = (N*kT + W/3) / V ,   W = sum over pairs of  r_ij . f_ij

    where ``r_ij = r_i - r_j`` is the centre-to-centre separation of an interacting pair and
    ``f_ij`` is the force on cell ``i`` due to cell ``j`` (Newton's third law: the force on
    ``j`` is ``-f_ij``, so each unordered pair contributes once). Sign convention: a *repulsive*
    pair (``f_ij`` along ``+r_ij``) gives ``W > 0`` (positive pressure); a *cohesive* pair
    (``f_ij`` pulling ``i`` toward ``j``, opposite ``r_ij``) gives ``W < 0`` (tension). The
    athermal overdamped CBM uses ``kT = 0`` for the configurational pressure; pass the
    thermostat ``kT`` (J) to include the ideal-gas kinetic term.

    The aggregate surface tension is then recovered from the interior-vs-exterior pressure
    difference via Young-Laplace in the oracle
    ``validation.oracles.spheroid.surface_tension_bridge.surface_tension_from_pressure``
    (``sigma = dP / (1/R + 1/R')``) — this observable provides the measured ``P``; the oracle
    is the runtime-forbidden acceptance relation. Pure numpy; no HOOMD, no physics constants.

    Args:
        positions: (N, 3) cell-center positions [m].
        pair_indices: (M, 2) integer array of interacting pair indices ``(i, j)`` into
            ``positions`` (each unordered pair listed once).
        pair_forces: (M, 3) force on cell ``i`` due to cell ``j`` for each pair [N].
        volume: aggregate volume [m^3] (> 0). Defaults to ``convex_hull_volume(positions)``.
        kT: thermal energy [J] for the kinetic term (default 0.0 = configurational only).

    Returns:
        Virial pressure P [Pa]. Positive = net repulsive (over-pressured); negative = net
        cohesive (under tension).

    Raises:
        ValueError: on shape mismatch, out-of-range indices, non-positive/zero volume, or
            ``kT < 0``.
    """
    p = _as_positions(positions)
    idx = np.asarray(pair_indices, dtype=np.int64)
    f = np.asarray(pair_forces, dtype=np.float64)
    if kT < 0.0:
        raise ValueError("kT must be non-negative [J].")
    if idx.ndim != 2 or idx.shape[1] != 2:
        raise ValueError(f"pair_indices must be (M, 2); got shape {idx.shape}.")
    if f.shape != idx.shape[:1] + (3,):
        raise ValueError(
            f"pair_forces must be (M, 3) matching pair_indices; got {f.shape}."
        )
    if idx.size and (idx.min() < 0 or idx.max() >= p.shape[0]):
        raise ValueError("pair_indices out of range for positions.")
    V = convex_hull_volume(p) if volume is None else float(volume)
    if V <= 0.0:
        raise ValueError("volume must be strictly positive [m^3] (degenerate hull?).")
    if idx.size:
        r_ij = p[idx[:, 0]] - p[idx[:, 1]]
        W = float(np.sum(np.sum(r_ij * f, axis=1)))
    else:
        W = 0.0
    n = p.shape[0]
    return float((n * kT + W / 3.0) / V)


def core_projected_area(
    positions: npt.ArrayLike, link_radius: float, plane: tuple[int, int] = (0, 1)
) -> float:
    """Projected area of the LARGEST connected component (fragmentation-robust footprint).

    The raw ``projected_area`` convex hull spans the gaps between drifting fragments and so
    explodes when a growing/active aggregate fragments. Restricting the hull to the largest
    connected component measures the spreading footprint of the cohesive core, ignoring
    flung-off pieces — the honest spread observable for a proliferating spheroid.

    Args:
        positions: (N, 3) cell centers (m).
        link_radius: single-linkage cluster threshold (m); see ``connected_components``.
        plane: the two axis indices to project onto (default x, y).

    Returns:
        Core convex-hull area in m^2 (0.0 if the largest component has < 3 non-collinear cells).
    """
    p = _as_positions(positions)
    mask = largest_connected_component(p, link_radius)
    return projected_area(p[mask], plane=plane)
