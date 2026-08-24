"""Local condensation fields + label-blind bundle LOCALIZATION — the real replacement for bundle_count.

``ff/architecture_metrics.bundle_count = n_fibers`` is VACUOUS: it "detects" a bundle by counting the
fibers that were handed to it under the label "bundle". A real emergence detector must find *where* an
ordered, dense patch has CONDENSED out of an isotropic background, reading geometry alone.

A stress fiber (STRESS_FIBER_TARGET) is a patch that is simultaneously ORIENTATIONALLY ORDERED and
SPATIALLY DENSE. So condensation is measured with two per-fiber spatial fields over a scale-adaptive
neighborhood:
  * ``s_local[f]`` — the validated nematic order S of the axes in fiber f's neighborhood (order);
  * ``density[f]`` — the neighbor count in that neighborhood (crowding).
A fiber is flagged CONDENSED when its local order beats the LOCAL finite-N isotropic null at the
pre-registered effect size (a dense-but-isotropic patch has high density yet fails the order test, and is
correctly NOT called a bundle). Flagged fibers are clustered by spatial connectivity into localized
bundles. No step reads a region/type label (the anti-coupling firewall).

Neighborhood radius (Magic-Number Block — ``NEIGHBORHOOD_SCALE``):
  value        : 3.0 (dimensionless multiplier of the median nearest-neighbor centroid spacing).
  derivable    : r = k * s_med with s_med the config's own median NN spacing; a ball of r ~ 3 s_med holds
                 O(10) neighbors in 3D, enough for a stable local Q-tensor without blurring the patch.
  grid-invariant: yes — r scales with the config's own spacing, so a refined/coarsened network gives the
                 same neighbor statistics; the localization is verified stable across a k-sweep (gate).
  NOT tuned to an outcome: chosen from neighbor-count sufficiency, not adjusted so a planted bundle is
                 recovered; the recovery gate must hold across the whole k-plateau.

Sanity Gate (self-tested in tests/ac/emergence/test_condensation_oracle.py):
  * a planted compact aligned bundle is localized (centroid within a spacing; recall/precision high) while
    the isotropic background is not flagged (specificity);
  * localization is rotation-equivariant (rotate the config -> the recovered centroid rotates by the same R);
  * localization is permutation-invariant (relabel fibers -> the same recovered member SET);
  * the recovery holds across the NEIGHBORHOOD_SCALE k-plateau (grid-invariance, not a single tuned k);
  * a dense but ISOTROPIC patch is NOT flagged (order, not mere crowding).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from aleph.components.emergence.nematic import fiber_axes, fiber_centroids, q_tensor
from aleph.components.emergence.null_model import EFFECT_SIZE_Z_CRIT, nematic_null_band

__all__ = [
    "NEIGHBORHOOD_SCALE",
    "MIN_BUNDLE_NEIGHBORS",
    "ALIGN_COS_MIN",
    "LOCAL_NULL_MC",
    "LocalFields",
    "Bundle",
    "neighborhood_radius",
    "compute_local_fields",
    "localize_bundles",
]

# See the Magic-Number Block in the module docstring.
NEIGHBORHOOD_SCALE: float = 3.0
# A bundle needs several co-aligned fibers; a neighborhood with fewer members (incl. self) cannot give a
# meaningful local order and is never flagged. This is a structural sufficiency floor, not a tuned knob.
MIN_BUNDLE_NEIGHBORS: int = 4
# Membership half-cone: a condensed fiber must lie within 45 deg of its neighborhood director, i.e.
# |n . director| > 1/sqrt(2)  (cos^2 > 1/2 -> the fiber sits in the ALIGNED hemisphere-cone, closer to the
# director than a perpendicular one). This separates a true bundle MEMBER from a background fiber that
# merely sits next to the ordered patch (its neighborhood looks ordered, but its OWN axis is random: a
# uniform axis passes with probability only 1 - 1/sqrt(2) ~ 0.29). Geometric, dimensionless, grid-invariant;
# NOT tuned to a recovery score (the recovery gate holds across the k-plateau at this fixed cone).
ALIGN_COS_MIN: float = 1.0 / np.sqrt(2.0)
# MC draws for a LOCAL null band (per distinct neighbor count). Smaller than the global band is fine: the
# local threshold needs sigma to a few %, far below the z-margin. Cached per neighbor count.
LOCAL_NULL_MC: int = 1000
_LOCAL_NULL_SEED: int = 512


@dataclass(frozen=True, slots=True)
class LocalFields:
    """Per-fiber spatial condensation fields.

    Attributes:
        centroids: (F, 3) fiber centroids.
        axes: (F, 3) fiber unit axes.
        radius: the neighborhood radius used.
        density: (F,) neighbor count within ``radius`` (includes the fiber itself).
        s_local: (F,) local nematic order S over the neighborhood axes (nan if under MIN_BUNDLE_NEIGHBORS).
        align: (F,) |n_f . director_local| — how aligned the fiber's OWN axis is with its neighborhood
            director (nan if under MIN_BUNDLE_NEIGHBORS). Separates a member from a neighbor of the patch.
    """

    centroids: npt.NDArray[np.float64]
    axes: npt.NDArray[np.float64]
    radius: float
    density: npt.NDArray[np.int64]
    s_local: npt.NDArray[np.float64]
    align: npt.NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class Bundle:
    """A localized condensed bundle recovered by the detector.

    Attributes:
        members: (k,) fiber indices in this bundle.
        center: (3,) centroid of the member fibers.
        axis: (3,) nematic director of the member fibers.
        s_bundle: nematic order S of the member fibers (whole-bundle, not per-neighborhood).
        size: number of member fibers.
    """

    members: npt.NDArray[np.int64]
    center: npt.NDArray[np.float64]
    axis: npt.NDArray[np.float64]
    s_bundle: float
    size: int


def neighborhood_radius(
    centroids: npt.NDArray[np.float64], scale: float = NEIGHBORHOOD_SCALE
) -> float:
    """Scale-adaptive neighborhood radius ``r = scale * median_nearest_neighbor_spacing``.

    Args:
        centroids: (F, 3) fiber centroids.
        scale: dimensionless multiplier (default :data:`NEIGHBORHOOD_SCALE`).

    Returns:
        The neighborhood radius in the config's length unit.
    """
    if centroids.shape[0] < 2:
        raise ValueError("need >= 2 fibers to define a spacing")
    d, _ = cKDTree(centroids).query(centroids, k=2)
    s_med = float(np.median(d[:, 1]))
    return scale * s_med


def _local_null_bands(neighbor_counts: npt.NDArray[np.int64]) -> dict[int, tuple[float, float]]:
    """Cache the LOCAL isotropic null (mean, std) of S for each distinct neighbor count present."""
    bands: dict[int, tuple[float, float]] = {}
    for m in np.unique(neighbor_counts):
        m_int = int(m)
        if m_int < MIN_BUNDLE_NEIGHBORS:
            continue
        band = nematic_null_band(m_int, n_mc=LOCAL_NULL_MC, seed=_LOCAL_NULL_SEED + m_int)
        bands[m_int] = (band.mean, band.std)
    return bands


def compute_local_fields(
    pos: npt.NDArray[np.float64],
    fiber_offsets: npt.NDArray[np.int64],
    radius: float | None = None,
    scale: float = NEIGHBORHOOD_SCALE,
) -> LocalFields:
    """Per-fiber neighborhood density + local nematic order S (label-blind geometry only).

    Args:
        pos: (N, 3) node positions (FF fiber-array contract).
        fiber_offsets: (F+1,) offsets.
        radius: neighborhood radius; if ``None`` derived from the config spacing via ``scale``.
        scale: neighborhood multiplier used when ``radius`` is None.

    Returns:
        A :class:`LocalFields` with density and local S per fiber.
    """
    centroids = fiber_centroids(pos, fiber_offsets)
    axes = fiber_axes(pos, fiber_offsets)
    r = neighborhood_radius(centroids, scale) if radius is None else float(radius)
    tree = cKDTree(centroids)
    neighbor_lists = tree.query_ball_point(centroids, r)
    n_fibers = centroids.shape[0]
    density = np.empty(n_fibers, dtype=np.int64)
    s_local = np.full(n_fibers, np.nan, dtype=np.float64)
    align = np.full(n_fibers, np.nan, dtype=np.float64)
    for f in range(n_fibers):
        idx = neighbor_lists[f]
        density[f] = len(idx)
        if len(idx) >= MIN_BUNDLE_NEIGHBORS:
            evals, evecs = np.linalg.eigh(q_tensor(axes[idx]))
            s_local[f] = float(evals.max())
            director = evecs[:, int(np.argmax(evals))]
            align[f] = float(abs(axes[f] @ director))
    return LocalFields(
        centroids=centroids, axes=axes, radius=r, density=density, s_local=s_local, align=align
    )


def localize_bundles(
    fields: LocalFields,
    z_crit: float = EFFECT_SIZE_Z_CRIT,
    min_size: int = MIN_BUNDLE_NEIGHBORS,
) -> list[Bundle]:
    """Flag fibers whose LOCAL order beats the local isotropic null at ``z_crit``, cluster them spatially.

    A fiber ``f`` is condensed iff (a) its neighborhood order beats the local isotropic null,
    ``s_local[f] > mu_null(m_f) + z_crit * sigma_null(m_f)`` for its own neighbor count ``m_f``, AND
    (b) the fiber's own axis is aligned with the neighborhood director, ``align[f] > ALIGN_COS_MIN`` (it
    is a MEMBER of the ordered patch, not merely adjacent to it). Both criteria are fully label-blind and
    pre-registered. Flagged fibers within one neighborhood radius of each other are grouped into connected
    clusters; clusters smaller than ``min_size`` are dropped. Bundles are returned largest-first.

    Args:
        fields: the :class:`LocalFields` from :func:`compute_local_fields`.
        z_crit: pre-registered effect size (default :data:`EFFECT_SIZE_Z_CRIT`).
        min_size: minimum member count for a reported bundle.

    Returns:
        Recovered :class:`Bundle` list, sorted by size descending.
    """
    bands = _local_null_bands(fields.density)
    flagged = np.zeros(fields.centroids.shape[0], dtype=bool)
    for f in range(fields.centroids.shape[0]):
        m = int(fields.density[f])
        if m not in bands or not np.isfinite(fields.s_local[f]):
            continue
        mu, sigma = bands[m]
        ordered = sigma > 0.0 and (fields.s_local[f] - mu) / sigma > z_crit
        aligned = fields.align[f] > ALIGN_COS_MIN
        if ordered and aligned:
            flagged[f] = True

    flagged_idx = np.flatnonzero(flagged)
    if flagged_idx.size == 0:
        return []

    # cluster flagged fibers by spatial connectivity (edges = flagged pairs within the neighborhood radius)
    sub = fields.centroids[flagged_idx]
    pairs = cKDTree(sub).query_pairs(fields.radius, output_type="ndarray")
    n_sub = flagged_idx.size
    if pairs.size:
        rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
        cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
        adj = coo_matrix((np.ones(rows.size), (rows, cols)), shape=(n_sub, n_sub))
    else:
        adj = coo_matrix((n_sub, n_sub))
    n_comp, comp = connected_components(adj, directed=False)

    bundles: list[Bundle] = []
    for c in range(n_comp):
        members = flagged_idx[comp == c]
        if members.size < min_size:
            continue
        member_axes = fields.axes[members]
        evals, evecs = np.linalg.eigh(q_tensor(member_axes))
        axis = evecs[:, int(np.argmax(evals))]
        bundles.append(
            Bundle(
                members=members,
                center=fields.centroids[members].mean(axis=0),
                axis=axis,
                s_bundle=float(evals.max()),
                size=int(members.size),
            )
        )
    bundles.sort(key=lambda b: b.size, reverse=True)
    return bundles
