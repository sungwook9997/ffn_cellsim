"""Synthetic fiber-array configs — the oracle fixtures for the emergence detector (pure NumPy).

Every generator emits the FF fiber-array geometry contract ``(pos, fiber_offsets)`` (and nothing else
label-like), so the detector consumes them exactly as it will consume a native network. Where a config
has a KNOWN planted structure (e.g. a bundle at a known place), the generator ALSO returns a ground-truth
member mask / center — but that ground truth is for the *test* to score against; it is never handed to
the detector (the anti-coupling firewall: the detector sees geometry only).

These configs span the acceptance corners:
  * ``isotropic_network``          -> S ~ 0 null band;
  * ``aligned_network``            -> S = 1 (incl. an ``antiparallel_frac`` knob: a graded-polarity SF
                                     bundle still reads S ~ 1 nematically while a polar mean collapses);
  * ``partially_aligned_network``  -> S monotone in the aligned fraction (condensation-fraction knob);
  * ``planted_bundle``             -> a compact aligned patch inside an isotropic background, order
                                     SHUFFLED so the planted set is not a contiguous index block
                                     (localization + specificity);
  * ``random_rotation`` / ``rotate_config`` / ``permute_fibers`` -> the invariance/equivariance fixtures.

All randomness takes an explicit seed/Generator so every fixture is reproducible.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "isotropic_network",
    "aligned_network",
    "partially_aligned_network",
    "planted_bundle",
    "dense_isotropic_patch",
    "random_rotation",
    "rotate_config",
    "permute_fibers",
]

# Minimal straight fiber discretization (>=2 required by the axis contract; 3 exercises the
# multi-node centroid path). These are geometric fixture knobs, not physical parameters.
_DEFAULT_NODES_PER_FIBER: int = 3
_DEFAULT_SEG_LEN: float = 1.0
_DEFAULT_BOX: float = 10.0


def _rand_axes(n: int, rng: np.random.Generator) -> npt.NDArray[np.float64]:
    """``n`` unit vectors uniform on the sphere (normalized Gaussians)."""
    g = rng.standard_normal((n, 3))
    return g / np.linalg.norm(g, axis=1, keepdims=True)


def _jitter_axes(
    axes: npt.NDArray[np.float64], jitter_deg: float, rng: np.random.Generator
) -> npt.NDArray[np.float64]:
    """Perturb each axis by a small Gaussian tangent kick of ~``jitter_deg`` degrees, renormalized."""
    if jitter_deg <= 0.0:
        return axes
    sigma = np.deg2rad(jitter_deg)
    perturbed = axes + sigma * rng.standard_normal(axes.shape)
    return perturbed / np.linalg.norm(perturbed, axis=1, keepdims=True)


def _straight_fibers(
    centroids: npt.NDArray[np.float64],
    axes: npt.NDArray[np.float64],
    seg_len: float,
    n_nodes: int,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """Build straight ``n_nodes``-point fibers from centroids + unit axes + length.

    Fiber ``f`` spans ``centroid_f +/- (seg_len/2) axis_f``; its end-to-end vector is ``seg_len*axis_f``.

    Returns:
        ``(pos, fiber_offsets)`` in the FF fiber-array contract.
    """
    n_fibers = centroids.shape[0]
    ts = np.linspace(-0.5, 0.5, n_nodes)  # (n_nodes,)
    # (F, n_nodes, 3) = centroid + t*len*axis
    nodes = centroids[:, None, :] + ts[None, :, None] * (seg_len * axes)[:, None, :]
    pos = nodes.reshape(n_fibers * n_nodes, 3)
    fiber_offsets = np.arange(n_fibers + 1, dtype=np.int64) * n_nodes
    return pos, fiber_offsets


def _apply_antiparallel(
    axes: npt.NDArray[np.float64], antiparallel_frac: float, rng: np.random.Generator
) -> npt.NDArray[np.float64]:
    """Flip the end-to-end direction of ``antiparallel_frac`` of the axes (nematic-invariant test)."""
    if antiparallel_frac <= 0.0:
        return axes
    n = axes.shape[0]
    flip = rng.random(n) < antiparallel_frac
    out = axes.copy()
    out[flip] *= -1.0
    return out


def isotropic_network(
    n_fibers: int,
    seed: int,
    box: float = _DEFAULT_BOX,
    seg_len: float = _DEFAULT_SEG_LEN,
    n_nodes: int = _DEFAULT_NODES_PER_FIBER,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """Isotropic network: uniform centroids in ``[0, box]^3``, uniform-on-sphere axes (the S~0 null)."""
    rng = np.random.default_rng(seed)
    centroids = rng.uniform(0.0, box, size=(n_fibers, 3))
    axes = _rand_axes(n_fibers, rng)
    return _straight_fibers(centroids, axes, seg_len, n_nodes)


def aligned_network(
    n_fibers: int,
    seed: int,
    director: npt.ArrayLike = (0.0, 0.0, 1.0),
    antiparallel_frac: float = 0.0,
    jitter_deg: float = 0.0,
    box: float = _DEFAULT_BOX,
    seg_len: float = _DEFAULT_SEG_LEN,
    n_nodes: int = _DEFAULT_NODES_PER_FIBER,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """Aligned bundle along ``director`` (-> S=1). ``antiparallel_frac`` flips a fraction of directions
    (graded-polarity SF bundle: still S~1 nematically); ``jitter_deg`` adds small angular spread."""
    rng = np.random.default_rng(seed)
    d = np.asarray(director, dtype=np.float64)
    d = d / (np.linalg.norm(d) + 1e-12)
    axes = np.tile(d, (n_fibers, 1))
    axes = _apply_antiparallel(axes, antiparallel_frac, rng)
    axes = _jitter_axes(axes, jitter_deg, rng)
    centroids = rng.uniform(0.0, box, size=(n_fibers, 3))
    return _straight_fibers(centroids, axes, seg_len, n_nodes)


def partially_aligned_network(
    n_fibers: int,
    aligned_frac: float,
    seed: int,
    director: npt.ArrayLike = (0.0, 0.0, 1.0),
    box: float = _DEFAULT_BOX,
    seg_len: float = _DEFAULT_SEG_LEN,
    n_nodes: int = _DEFAULT_NODES_PER_FIBER,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """Mix of an aligned fraction along ``director`` and an isotropic remainder (S monotone in the
    aligned fraction — the condensation-fraction knob for the monotonicity gate)."""
    if not 0.0 <= aligned_frac <= 1.0:
        raise ValueError("aligned_frac must be in [0, 1]")
    rng = np.random.default_rng(seed)
    d = np.asarray(director, dtype=np.float64)
    d = d / (np.linalg.norm(d) + 1e-12)
    n_aligned = int(round(aligned_frac * n_fibers))
    axes = _rand_axes(n_fibers, rng)
    axes[:n_aligned] = d
    centroids = rng.uniform(0.0, box, size=(n_fibers, 3))
    return _straight_fibers(centroids, axes, seg_len, n_nodes)


def planted_bundle(
    n_background: int,
    n_bundle: int,
    seed: int,
    director: npt.ArrayLike = (0.0, 0.0, 1.0),
    bundle_center: npt.ArrayLike | None = None,
    bundle_radius: float = 1.5,
    antiparallel_frac: float = 0.5,
    jitter_deg: float = 5.0,
    box: float = _DEFAULT_BOX,
    seg_len: float = _DEFAULT_SEG_LEN,
    n_nodes: int = _DEFAULT_NODES_PER_FIBER,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64], npt.NDArray[np.bool_], npt.NDArray[np.float64]]:
    """Isotropic background + one compact aligned bundle patch; fiber order SHUFFLED.

    The bundle is a spatially localized, orientation-ordered (antiparallel-mixed, lightly jittered)
    patch of ``n_bundle`` fibers centred on ``bundle_center`` within ``bundle_radius``; the ``n_background``
    fibers fill ``[0, box]^3`` isotropically. Fiber order is randomly permuted so the planted set is NOT a
    contiguous index block (the detector cannot exploit ordering).

    Args:
        n_background: isotropic background fiber count.
        n_bundle: planted (aligned + localized) fiber count.
        seed: RNG seed.
        director: bundle nematic axis.
        bundle_center: (3,) patch centre; defaults to the box centre.
        bundle_radius: patch radius (compact => localized).
        antiparallel_frac: fraction of bundle fibers with flipped end-to-end direction (SF graded polarity).
        jitter_deg: small angular spread of the bundle axes.
        box, seg_len, n_nodes: geometric fixture knobs.

    Returns:
        ``(pos, fiber_offsets, planted_mask, bundle_center)`` where ``planted_mask`` (F,) is True on the
        planted fibers (ground truth for scoring — NEVER passed to the detector) and ``bundle_center`` is
        the (3,) patch centre.
    """
    rng = np.random.default_rng(seed)
    d = np.asarray(director, dtype=np.float64)
    d = d / (np.linalg.norm(d) + 1e-12)
    center = np.full(3, box / 2.0) if bundle_center is None else np.asarray(bundle_center, dtype=np.float64)

    bg_centroids = rng.uniform(0.0, box, size=(n_background, 3))
    bg_axes = _rand_axes(n_background, rng)

    # bundle centroids: uniform within a ball of bundle_radius around center
    dirs = _rand_axes(n_bundle, rng)
    radii = bundle_radius * rng.random(n_bundle) ** (1.0 / 3.0)  # uniform-in-ball radial law
    bundle_centroids = center[None, :] + radii[:, None] * dirs
    bundle_axes = _apply_antiparallel(np.tile(d, (n_bundle, 1)), antiparallel_frac, rng)
    bundle_axes = _jitter_axes(bundle_axes, jitter_deg, rng)

    centroids = np.concatenate([bg_centroids, bundle_centroids], axis=0)
    axes = np.concatenate([bg_axes, bundle_axes], axis=0)
    mask = np.concatenate([np.zeros(n_background, bool), np.ones(n_bundle, bool)])

    perm = rng.permutation(n_background + n_bundle)
    centroids, axes, mask = centroids[perm], axes[perm], mask[perm]
    pos, fiber_offsets = _straight_fibers(centroids, axes, seg_len, n_nodes)
    return pos, fiber_offsets, mask, center


def dense_isotropic_patch(
    n_background: int,
    n_patch: int,
    seed: int,
    bundle_center: npt.ArrayLike | None = None,
    bundle_radius: float = 1.5,
    box: float = _DEFAULT_BOX,
    seg_len: float = _DEFAULT_SEG_LEN,
    n_nodes: int = _DEFAULT_NODES_PER_FIBER,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64], npt.NDArray[np.bool_], npt.NDArray[np.float64]]:
    """Isotropic background + a spatially DENSE but ORIENTATIONALLY ISOTROPIC patch (the specificity
    control: crowding without order must NOT be called a bundle). Same geometry as ``planted_bundle``
    but the patch axes are random, and order is SHUFFLED.

    Returns:
        ``(pos, fiber_offsets, patch_mask, patch_center)`` (mask is ground truth for the test, never the
        detector).
    """
    rng = np.random.default_rng(seed)
    center = np.full(3, box / 2.0) if bundle_center is None else np.asarray(bundle_center, dtype=np.float64)
    bg_centroids = rng.uniform(0.0, box, size=(n_background, 3))
    dirs = _rand_axes(n_patch, rng)
    radii = bundle_radius * rng.random(n_patch) ** (1.0 / 3.0)
    patch_centroids = center[None, :] + radii[:, None] * dirs
    centroids = np.concatenate([bg_centroids, patch_centroids], axis=0)
    axes = _rand_axes(n_background + n_patch, rng)  # ALL random orientations (no order anywhere)
    mask = np.concatenate([np.zeros(n_background, bool), np.ones(n_patch, bool)])
    perm = rng.permutation(n_background + n_patch)
    centroids, axes, mask = centroids[perm], axes[perm], mask[perm]
    pos, fiber_offsets = _straight_fibers(centroids, axes, seg_len, n_nodes)
    return pos, fiber_offsets, mask, center


def random_rotation(rng: np.random.Generator) -> npt.NDArray[np.float64]:
    """A proper (det=+1) rotation matrix from the QR of a Gaussian matrix (for invariance fixtures)."""
    a = rng.standard_normal((3, 3))
    q, r = np.linalg.qr(a)
    q = q @ np.diag(np.sign(np.diag(r)))  # fix the QR sign ambiguity -> near-Haar
    if np.linalg.det(q) < 0.0:
        q[:, 0] = -q[:, 0]
    return q


def rotate_config(
    pos: npt.NDArray[np.float64], rot: npt.NDArray[np.float64], about: npt.ArrayLike | None = None
) -> npt.NDArray[np.float64]:
    """Rotate all node positions by ``rot`` about the point ``about`` (default: the node centroid)."""
    pos = np.asarray(pos, dtype=np.float64)
    pivot = pos.mean(axis=0) if about is None else np.asarray(about, dtype=np.float64)
    return (pos - pivot) @ rot.T + pivot


def permute_fibers(
    pos: npt.NDArray[np.float64], fiber_offsets: npt.NDArray[np.int64], perm: npt.NDArray[np.int64]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """Reorder whole fibers by ``perm`` (a relabeling), returning a new contiguous ``(pos, offsets)``.

    Args:
        pos: (N, 3) node positions.
        fiber_offsets: (F+1,) offsets.
        perm: (F,) permutation of ``range(F)`` (the new fiber order).

    Returns:
        ``(pos_permuted, fiber_offsets_permuted)`` — the SAME network with fibers relabeled.
    """
    pos = np.asarray(pos, dtype=np.float64)
    fiber_offsets = np.asarray(fiber_offsets, dtype=np.int64)
    perm = np.asarray(perm, dtype=np.int64)
    blocks = [pos[fiber_offsets[f]:fiber_offsets[f + 1]] for f in perm]
    new_pos = np.concatenate(blocks, axis=0)
    new_offsets = np.zeros(len(perm) + 1, dtype=np.int64)
    new_offsets[1:] = np.cumsum([b.shape[0] for b in blocks])
    return new_pos, new_offsets
