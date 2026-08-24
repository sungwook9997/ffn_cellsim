"""Nematic order parameter S — the VALIDATED Q-tensor form, reused against the fiber-array contract.

This is the exact orientation-invariant order parameter already validated in
``ff/architecture_metrics.parallel_order_parameter``:

    Q = (3 <n (x) n> - I) / 2 ,     S = lambda_max(Q)

built from the per-fiber unit end-to-end axes ``n``. ``S = 1`` for a perfectly parallel bundle,
``S = 0`` for an isotropic orientation distribution (in the infinite-fiber limit; finite-N gives a
small positive bias that ``null_model`` quantifies). The measure is NEMATIC, not polar: it is built
from ``n (x) n`` so it is invariant under ``n -> -n`` and therefore reports full order for an
*antiparallel* stress-fiber bundle (graded / mixed polarity) that a polar mean ``|<n>|`` would miss —
exactly the SF regime of STRESS_FIBER_TARGET.

Everything here is a pure function of GEOMETRY (the FF fiber-array contract ``pos`` / ``fiber_offsets``)
— never of region / type / construction labels (the I5 anti-coupling firewall).

Sanity Gate (self-tested in tests/ac/emergence/test_nematic_oracle.py):
  * dimensional/convention: S is dimensionless in [-1/2, 1]; aligned axes -> S = 1 exactly;
    exact-isotropic <n(x)n> = I/3 -> Q = 0 -> S = 0.
  * rotational invariance: S(R.axes) == S(axes) to machine precision (eigenvalues are rotation-invariant).
  * permutation/label invariance: S is invariant to fiber reordering (a symmetric sum).
  * n -> -n (polarity) invariance: flipping any fiber's end-to-end direction does not change S.
  * parity with the validated ``parallel_order_parameter`` on a shared network.

Reference: de Gennes & Prost, The Physics of Liquid Crystals (Q-tensor / scalar order parameter);
matches ``ff/architecture_metrics.parallel_order_parameter`` (the project's validated form).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "fiber_axes",
    "fiber_centroids",
    "fiber_lengths",
    "q_tensor",
    "nematic_order",
    "nematic_order_and_director",
    "AXIS_EPS",
]

# Tiny length used only to avoid 0/0 when a degenerate fiber has coincident end nodes; it never
# changes a well-posed answer (matches the 1e-12 guard in ff/architecture_metrics._fiber_axis).
AXIS_EPS: float = 1.0e-12


def _validate_fiber_arrays(pos: npt.NDArray[np.float64], fiber_offsets: npt.NDArray[np.int64]) -> None:
    """Check the fiber-array contract shapes (label-blind: geometry + topology only)."""
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"pos must be (N, 3); got {pos.shape}")
    if fiber_offsets.ndim != 1 or fiber_offsets.shape[0] < 2:
        raise ValueError(f"fiber_offsets must be (F+1,) with F>=1; got {fiber_offsets.shape}")
    if fiber_offsets[0] != 0 or fiber_offsets[-1] != pos.shape[0]:
        raise ValueError("fiber_offsets must start at 0 and end at N (contiguous node ownership)")
    if np.any(np.diff(fiber_offsets) < 2):
        raise ValueError("every fiber needs >= 2 nodes (an end-to-end axis is undefined otherwise)")


def fiber_axes(
    pos: npt.NDArray[np.float64], fiber_offsets: npt.NDArray[np.int64]
) -> npt.NDArray[np.float64]:
    """Unit end-to-end axis of every fiber (the vectorized ``ff`` ``_fiber_axis``).

    Fiber ``f`` owns nodes ``pos[fiber_offsets[f]:fiber_offsets[f+1]]``; its axis is
    ``pos[last] - pos[first]`` normalized. Note this is a *directed* axis (it carries the head->tail
    sign), but every downstream nematic measure uses ``n (x) n`` and is therefore sign-blind.

    Args:
        pos: (N, 3) node positions [any length unit].
        fiber_offsets: (F+1,) int offsets; fiber ``f`` = ``pos[fiber_offsets[f]:fiber_offsets[f+1]]``.

    Returns:
        (F, 3) unit axes.
    """
    pos = np.asarray(pos, dtype=np.float64)
    fiber_offsets = np.asarray(fiber_offsets, dtype=np.int64)
    _validate_fiber_arrays(pos, fiber_offsets)
    first = fiber_offsets[:-1]
    last = fiber_offsets[1:] - 1
    vecs = pos[last] - pos[first]
    norms = np.linalg.norm(vecs, axis=1)
    return vecs / (norms[:, None] + AXIS_EPS)


def fiber_centroids(
    pos: npt.NDArray[np.float64], fiber_offsets: npt.NDArray[np.int64]
) -> npt.NDArray[np.float64]:
    """Centroid (mean of nodes) of every fiber — the spatial anchor for local condensation.

    A centroid (rather than the middle node) is used so a curved/kinked fiber still gets a stable
    location; for a straight 2-node fiber it coincides with the segment midpoint.

    Args:
        pos: (N, 3) node positions.
        fiber_offsets: (F+1,) int offsets.

    Returns:
        (F, 3) fiber centroids.
    """
    pos = np.asarray(pos, dtype=np.float64)
    fiber_offsets = np.asarray(fiber_offsets, dtype=np.int64)
    _validate_fiber_arrays(pos, fiber_offsets)
    sums = np.add.reduceat(pos, fiber_offsets[:-1], axis=0)
    counts = np.diff(fiber_offsets).astype(np.float64)
    return sums / counts[:, None]


def fiber_lengths(
    pos: npt.NDArray[np.float64], fiber_offsets: npt.NDArray[np.int64]
) -> npt.NDArray[np.float64]:
    """End-to-end length of every fiber (diagnostic; NOT used to weight the validated S)."""
    pos = np.asarray(pos, dtype=np.float64)
    fiber_offsets = np.asarray(fiber_offsets, dtype=np.int64)
    _validate_fiber_arrays(pos, fiber_offsets)
    first = fiber_offsets[:-1]
    last = fiber_offsets[1:] - 1
    return np.linalg.norm(pos[last] - pos[first], axis=1)


def q_tensor(
    axes: npt.NDArray[np.float64], weights: npt.NDArray[np.float64] | None = None
) -> npt.NDArray[np.float64]:
    """Symmetric traceless nematic Q-tensor ``Q = (3 <n (x) n> - I) / 2``.

    With ``weights=None`` every axis contributes equally — the exact validated form. Weights are
    accepted only so a local neighborhood can be down-weighted smoothly; the global gate uses equal
    weights to stay bit-parity with ``parallel_order_parameter``.

    Args:
        axes: (F, 3) axes (need not be exactly unit; the outer product is used as-is).
        weights: optional (F,) non-negative weights; ``None`` => uniform.

    Returns:
        (3, 3) symmetric traceless Q.
    """
    axes = np.asarray(axes, dtype=np.float64)
    if axes.ndim != 2 or axes.shape[1] != 3:
        raise ValueError(f"axes must be (F, 3); got {axes.shape}")
    if axes.shape[0] == 0:
        raise ValueError("need >= 1 axis to form Q")
    if weights is None:
        mean_outer = np.einsum("ni,nj->ij", axes, axes) / axes.shape[0]
    else:
        w = np.asarray(weights, dtype=np.float64)
        if w.shape != (axes.shape[0],):
            raise ValueError("weights must be (F,)")
        wsum = w.sum()
        if wsum <= 0.0:
            raise ValueError("weights must sum to > 0")
        mean_outer = np.einsum("n,ni,nj->ij", w, axes, axes) / wsum
    return (3.0 * mean_outer - np.eye(3)) / 2.0


def nematic_order(
    axes: npt.NDArray[np.float64], weights: npt.NDArray[np.float64] | None = None
) -> float:
    """Scalar nematic order ``S = lambda_max(Q)`` — the validated orientation-invariant order parameter.

    Args:
        axes: (F, 3) fiber axes.
        weights: optional (F,) weights; ``None`` => the validated uniform form.

    Returns:
        ``S`` in [-1/2, 1]; 1 = perfectly parallel, ~0 = isotropic (finite-N bias in ``null_model``).
    """
    return float(np.linalg.eigvalsh(q_tensor(axes, weights)).max())


def nematic_order_and_director(
    axes: npt.NDArray[np.float64], weights: npt.NDArray[np.float64] | None = None
) -> tuple[float, npt.NDArray[np.float64]]:
    """``(S, director)`` where the director is the unit eigenvector of the largest eigenvalue of Q.

    The director sign is arbitrary (nematic n -> -n symmetry); it is fixed only so the first non-zero
    component is non-negative, for reproducible plots.

    Args:
        axes: (F, 3) fiber axes.
        weights: optional (F,) weights.

    Returns:
        Tuple of scalar order ``S`` and the (3,) unit director.
    """
    evals, evecs = np.linalg.eigh(q_tensor(axes, weights))
    director = evecs[:, int(np.argmax(evals))]
    nz = director[np.abs(director) > AXIS_EPS]
    if nz.size and nz[0] < 0:
        director = -director
    return float(evals.max()), director
