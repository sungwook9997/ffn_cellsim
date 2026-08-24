"""Emergence observers: nematic order, a finite-N null band, and bundle condensation.

`ALEPH-PORT-3602`.

A renderer can show *what* the step did. It cannot, by drawing, say whether the structure it draws is
**ordered** or is a picture of noise. This module answers that, and it is built around one fact that
makes the naive version of the answer wrong:

    An isotropic population does not score zero.

For ``N`` directors drawn uniformly on the sphere the nematic order parameter has expectation
``~0.803 / sqrt(N)`` -- at ``N = 16`` that is ``0.21``, at ``N = 64`` it is ``0.10``. So a fixed
threshold ("call it ordered above 0.3") fires on pure noise at small populations and can never fire at
large ones. **The threshold is not a property of the structure; it is a property of the sample size.**
Every number here exists to replace that threshold with a null band that knows ``N``.

What this module is, and is not
-------------------------------
It is an **observer**. Nothing here returns a force, a torque, a stress or an energy; nothing here is
called from an accepted step; and no module under ``aleph/vertical/``, ``aleph/runtime/``,
``aleph/state/`` or ``aleph/scenarios/`` may import it. That last rule is asserted structurally by
``tests/viz/test_emergence.py::test_no_runtime_or_vertical_module_imports_the_observer``, which walks
those trees with ``ast`` rather than trusting anyone to remember. A measurement that can reach the
dynamics is a measurement that can be tuned until the dynamics agree with it.

The three layers
----------------
====================  ===================================================================
:func:`nematic_order`  ``Q = <(3/2) n(x)n - (1/2) I>``, ``S = lambda_max(Q)``. Sign-blind,
                       so an antiparallel bundle reads as ordered -- which it is.
:func:`isotropic_null_band`
                       What ``S`` would have been if the population were disordered. The
                       part that turns a number into a measurement.
:func:`condensed_bundles`
                       Where an ordered *and* dense patch has separated from the
                       background, scored against a **local** null band.
====================  ===================================================================

Two closed forms, and the difference between them matters
---------------------------------------------------------
:func:`isotropic_eigenvalue_scale_sq` is **exact at every N** -- it is algebra, derived in its own
docstring, and checkable by hand at ``N = 1``. It is what keeps the Monte-Carlo sampler honest.

:func:`isotropic_band_asymptote` gives closed forms for the band's own mean and standard deviation.
Those are **asymptotic**: they invoke the central limit theorem, and the measured relative excess of
the sampled mean over the asymptote runs ``+7.6%`` at ``N = 4`` down to ``+0.35%`` at ``N = 4096``.

**So the detector uses the sampled band and never the asymptote.** At small ``N`` the asymptote
understates the mean, which inflates every ``z`` and biases the detector toward *firing* -- the exact
failure this module exists to prevent. The closed form's job is to falsify the sampler, not replace it.

Units: positions in um; every quantity this module returns is dimensionless or a count, except
:func:`neighbourhood_radius` and :func:`segment_lengths`, which are um.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import numpy.typing as npt

__all__ = [
    "ALIGNMENT_COSINE_MIN",
    "MIN_NEIGHBOURS",
    "NEIGHBOURHOOD_SCALE",
    "ORDER_Z_CRITICAL",
    "CondensedBundle",
    "DegenerateSegmentError",
    "EmergenceObservation",
    "LocalFields",
    "NullBand",
    "condensed_bundles",
    "isotropic_band_asymptote",
    "isotropic_eigenvalue_scale_sq",
    "isotropic_null_band",
    "local_fields",
    "nematic_order",
    "nematic_order_and_director",
    "neighbourhood_radius",
    "observe_emergence",
    "order_excess_z",
    "q_tensor",
    "sample_isotropic_directors",
    "segment_directors",
    "segment_lengths",
    "segment_midpoints",
    "smallest_decidable_population",
]


class DegenerateSegmentError(ValueError):
    """A segment whose endpoints coincide, so it has no director.

    Raised rather than regularised. An epsilon in the denominator would turn a defective segment into
    an arbitrary unit vector that then votes in the order parameter, and the vote would be
    indistinguishable from a real one. A degenerate segment in a filament population is a defect in the
    population, and a measurement that hides it is worse than one that stops.
    """


#: Number of null standard deviations a population's order must exceed before it is called ordered.
#:
#: This is a **statistical significance threshold, not a physical constant.** It is declared here, in
#: advance, so that it is a bar the structure must clear rather than a dial turned until a run passes.
#: It carries no units, no length scale and no population size -- the population size is already
#: divided out, which is the whole point of scoring against a band instead of a threshold.
ORDER_Z_CRITICAL: float = 5.0

#: Neighbourhood radius as a multiple of the configuration's own median nearest-neighbour spacing.
#:
#: Expressing it as a multiple rather than in um is what makes the detector grid-invariant: refine or
#: coarsen a filament population and the neighbourhoods hold the same number of neighbours, so the
#: local null band is the same band. A value near 3 puts of order ten neighbours in a ball in 3D --
#: enough for a stable local Q-tensor without blurring one patch into the next.
NEIGHBOURHOOD_SCALE: float = 3.0

#: Half-cone for bundle membership: a member's own director must satisfy ``|n . d| > 1/sqrt(2)``.
#:
#: 45 degrees, and it separates a *member* of an ordered patch from a segment that merely sits beside
#: one. A background segment next to a bundle inherits an ordered neighbourhood but keeps its own random
#: direction, and a uniformly-oriented director clears this cone with probability only ``1 - 1/sqrt(2)``,
#: about 0.29.
ALIGNMENT_COSINE_MIN: float = 1.0 / math.sqrt(2.0)


def smallest_decidable_population(z_critical: float = ORDER_Z_CRITICAL) -> int:
    """Asymptotic **lower bound** on the ``N`` at which "ordered" is a decidable question.

    The decision threshold is ``mu(N) + z * sigma(N)``, and both terms fall as ``1/sqrt(N)`` while ``S``
    is bounded above by 1. So below some population the threshold sits **above the largest value S can
    take**, and the test is not strict -- it is vacuous. No arrangement of 4 directors, however
    perfectly aligned, can clear a 5-sigma bar built from 4 directors.

    Using the asymptotic band of :func:`isotropic_band_asymptote`, the threshold is
    ``(9/(2 sqrt(10 pi)) + z sqrt((29 pi - 81)/(40 pi))) / sqrt(N)``, so ``threshold < 1`` gives
    ``N > (0.802866 + 0.283590 z)^2``, i.e. ``N >= 5`` at ``z = 5``.

    **This is a bound and not the answer, in the unsafe direction.** The asymptote understates the true
    mean at small ``N`` (see :func:`isotropic_band_asymptote`), so it understates the threshold and
    therefore *overstates* how small a decidable neighbourhood can be. Measured against the sampled
    band, ``N = 5`` gives a threshold of ``1.06`` and is still undecidable; the smallest genuinely
    decidable population at ``z = 5`` is **6**. That is the value :data:`MIN_NEIGHBOURS` carries, and
    ``test_the_asymptotic_bound_is_optimistic_and_the_constant_is_the_measured_one`` pins both numbers.

    The discrepancy is worth stating plainly rather than quietly rounding up: it is this module's own
    argument -- that an asymptotic band flatters a small sample -- turning up in its own constant.

    Args:
        z_critical: The significance threshold the population must be able to clear.

    Returns:
        The smallest integer population for which the **asymptotic** threshold is attainable.
    """
    if not math.isfinite(z_critical) or z_critical < 0.0:
        raise ValueError(f"z_critical must be finite and non-negative; got {z_critical!r}")
    mu_coefficient = 9.0 / (2.0 * math.sqrt(10.0 * math.pi))
    sigma_coefficient = math.sqrt((29.0 * math.pi - 81.0) / (40.0 * math.pi))
    return int(math.floor((mu_coefficient + z_critical * sigma_coefficient) ** 2)) + 1


#: Fewest neighbours (including the segment itself) for which a local order test is decidable.
#:
#: **Measured, not assumed.** :func:`smallest_decidable_population` gives the asymptotic bound 5 at
#: ``ORDER_Z_CRITICAL``; the sampled band puts ``N = 5``'s threshold at ``1.06``, above the largest
#: value ``S`` can take, so 5 is still vacuous and the real answer is 6. Below this a neighbourhood's
#: local order is reported as ``nan`` -- undecidable -- and never as an unflagged measurement.
#:
#: :meth:`NullBand.is_decidable` re-checks this at runtime against the sampled band, so a caller who
#: changes ``z_critical`` is protected even though this constant was fixed at ``z = 5``.
MIN_NEIGHBOURS: int = 6


# ==================================================================================================
# Geometry: a filament population, read as segments
# ==================================================================================================
def _as_positions(positions: Any) -> npt.NDArray[np.float64]:
    array = np.asarray(positions, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"positions must be (N, 3) in um; got shape {array.shape}")
    if not np.isfinite(array).all():
        raise ValueError("positions carry a non-finite value")
    return array


def _as_segments(segments: Any, n_nodes: int) -> npt.NDArray[np.int64]:
    array = np.asarray(segments)
    if array.ndim != 2 or array.shape[1] < 2:
        raise ValueError(f"segments must be (M, 2) node-index pairs; got shape {array.shape}")
    pairs = array[:, :2].astype(np.int64)
    if pairs.shape[0] == 0:
        raise ValueError("a filament population with no segments has no orientation to measure")
    if pairs.min() < 0 or pairs.max() >= n_nodes:
        raise ValueError(
            f"segment indices span [{pairs.min()}, {pairs.max()}] but there are {n_nodes} nodes"
        )
    return pairs


def segment_directors(positions: Any, segments: Any) -> npt.NDArray[np.float64]:
    """Unit director of every segment.

    A segment ``(i, j)`` has director ``(p_j - p_i) / |p_j - p_i|``. The sign is arbitrary -- a filament
    has no head -- and every measure downstream is built from ``n (x) n``, which is blind to it.

    Args:
        positions: ``(N, 3)`` node positions [um].
        segments: ``(M, 2)`` node-index pairs.

    Returns:
        ``(M, 3)`` unit directors, dimensionless.

    Raises:
        DegenerateSegmentError: If any segment's endpoints coincide. Not regularised; see the class.
    """
    nodes = _as_positions(positions)
    pairs = _as_segments(segments, nodes.shape[0])
    deltas = nodes[pairs[:, 1]] - nodes[pairs[:, 0]]
    lengths = np.linalg.norm(deltas, axis=1)
    degenerate = np.flatnonzero(lengths == 0.0)
    if degenerate.size:
        first = int(degenerate[0])
        raise DegenerateSegmentError(
            f"segment {first} joins nodes {int(pairs[first, 0])} and {int(pairs[first, 1])}, which "
            f"are at the same point, so it has no direction ({degenerate.size} such segment(s)). "
            "Refusing rather than dividing by an epsilon: a fabricated director would vote in the "
            "order parameter and be indistinguishable from a real one."
        )
    return deltas / lengths[:, None]


def segment_midpoints(positions: Any, segments: Any) -> npt.NDArray[np.float64]:
    """Midpoint of every segment [um] -- the spatial anchor a local neighbourhood is built around."""
    nodes = _as_positions(positions)
    pairs = _as_segments(segments, nodes.shape[0])
    return 0.5 * (nodes[pairs[:, 0]] + nodes[pairs[:, 1]])


def segment_lengths(positions: Any, segments: Any) -> npt.NDArray[np.float64]:
    """End-to-end length of every segment [um]. Diagnostic; the order parameter is unweighted."""
    nodes = _as_positions(positions)
    pairs = _as_segments(segments, nodes.shape[0])
    return np.linalg.norm(nodes[pairs[:, 1]] - nodes[pairs[:, 0]], axis=1)


# ==================================================================================================
# The order parameter
# ==================================================================================================
def _as_directors(directors: Any, *, minimum: int = 2) -> npt.NDArray[np.float64]:
    array = np.asarray(directors, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"directors must be (M, 3); got shape {array.shape}")
    if array.shape[0] < minimum:
        raise ValueError(
            f"need at least {minimum} directors to measure order; got {array.shape[0]}. A single "
            "director scores 1 by construction and is never evidence of anything."
        )
    if not np.isfinite(array).all():
        raise ValueError("directors carry a non-finite value")
    return array


def q_tensor(directors: Any, weights: Any | None = None) -> npt.NDArray[np.float64]:
    """The symmetric traceless nematic tensor ``Q = (3/2) <n (x) n> - (1/2) I``.

    **The normalisation is fixed by two anchors, and nothing else is free.** Write
    ``M = <n (x) n>``, which is symmetric with unit trace.

    * *Perfect alignment*: every ``n = u``, so ``M = u (x) u`` with eigenvalues ``(1, 0, 0)``. Then
      ``Q`` has eigenvalues ``(1, -1/2, -1/2)`` and ``lambda_max(Q) = 1`` **exactly**.
    * *Perfect isotropy*: ``M = I/3`` by symmetry, so ``Q = 0`` and ``lambda_max(Q) = 0``.

    Those two requirements determine the ``3/2`` and the ``1/2``; any other pair of constants breaks
    one of them. ``Q`` is traceless because ``tr M = 1``.

    The construction is built from the outer product rather than from the mean director, so it is blind
    to ``n -> -n``. That is deliberate and is the reason this measure is used: a bundle whose filaments
    alternate in polarity is fully ordered, and a polar measure ``|<n>|`` would report it as disordered.

    Args:
        directors: ``(M, 3)`` directors, taken as given (they need not be exactly unit).
        weights: Optional ``(M,)`` non-negative weights summing to a positive number. Provided so a
            neighbourhood can be down-weighted smoothly; the global measure uses uniform weights.

    Returns:
        ``(3, 3)`` symmetric traceless ``Q``, dimensionless.
    """
    axes = _as_directors(directors)
    if weights is None:
        mean_outer = np.einsum("mi,mj->ij", axes, axes) / axes.shape[0]
    else:
        w = np.asarray(weights, dtype=np.float64)
        if w.shape != (axes.shape[0],):
            raise ValueError(f"weights must be ({axes.shape[0]},); got {w.shape}")
        if np.any(w < 0.0):
            raise ValueError("weights must be non-negative")
        total = float(w.sum())
        if not total > 0.0:
            raise ValueError("weights must sum to a positive number")
        mean_outer = np.einsum("m,mi,mj->ij", w, axes, axes) / total
    q = 1.5 * mean_outer - 0.5 * np.eye(3)
    return 0.5 * (q + q.T)


def nematic_order(directors: Any, weights: Any | None = None) -> float:
    """Scalar nematic order ``S = lambda_max(Q)``, dimensionless and in ``[-1/2, 1]``.

    ``S = 1`` is perfect alignment; ``S = 0`` is exact isotropy; ``S = -1/2`` is a perfectly oblate
    population confined to a plane, which is *ordered* in a different sense and is why the lower bound
    is not zero.

    **``S`` near zero does not mean disordered at finite population.** Score it with
    :func:`order_excess_z` against :func:`isotropic_null_band` before drawing that conclusion.
    """
    return float(np.linalg.eigvalsh(q_tensor(directors, weights)).max())


def nematic_order_and_director(
    directors: Any, weights: Any | None = None
) -> tuple[float, npt.NDArray[np.float64]]:
    """``(S, d)`` where ``d`` is the unit eigenvector of ``Q``'s largest eigenvalue.

    The director's sign is physically meaningless, so it is pinned -- first non-negligible component
    made positive -- purely so a figure redrawn from the same data is byte-identical.

    **The director is worse-conditioned than S, and the tests treat it that way.** When the two largest
    eigenvalues nearly coincide (a nearly planar population) the axis genuinely is not unique, and no
    amount of precision invents one. ``S`` remains well conditioned throughout.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(q_tensor(directors, weights))
    axis = np.asarray(eigenvectors[:, int(np.argmax(eigenvalues))], dtype=np.float64)
    significant = axis[np.abs(axis) > 1.0e-12]
    if significant.size and significant[0] < 0.0:
        axis = -axis
    return float(eigenvalues.max()), axis


# ==================================================================================================
# The finite-N isotropic null band
# ==================================================================================================
def isotropic_eigenvalue_scale_sq(n_directors: int) -> float:
    """``E[sum_i lambda_i^2] = 3 / (2N)`` for ``N`` isotropic directors. **Exact at every N.**

    Derivation, which is short enough to check by hand and is the reason this function exists:

    Let ``t = n (x) n - I/3`` be one sample's traceless part. Its Frobenius norm is **deterministic**,
    not random, because ``(n (x) n)^2 = n (x) n`` for a unit ``n``::

        |t|_F^2 = tr[(n(x)n)^2] - (2/3) tr[n(x)n] + (1/9) tr(I)
                = 1 - 2/3 + 1/3
                = 2/3

    The ``t_k`` are independent and zero-mean, so the sample mean's expected squared norm is
    ``E |M - I/3|_F^2 = (1/N^2) * N * (2/3) = 2/(3N)``. With ``Q = (3/2)(M - I/3)``::

        E[ sum_i lambda_i^2 ] = E |Q|_F^2 = (9/4) * 2/(3N) = 3 / (2N)

    No limit is taken anywhere, so this holds at **every** ``N >= 1``, and it is checkable without a
    computer: at ``N = 1`` a lone director gives ``Q`` eigenvalues ``(1, -1/2, -1/2)``, whose squares
    sum to ``3/2`` -- which is ``3/(2*1)``.

    That is what makes this the cross-check for :func:`isotropic_null_band`. A Monte-Carlo estimate of
    a quantity nobody can compute independently is a number on trust; this one has to match algebra.

    Args:
        n_directors: Population size ``N >= 1``.

    Returns:
        The exact expected sum of squared Q-eigenvalues, dimensionless.
    """
    if int(n_directors) < 1:
        raise ValueError(f"n_directors must be >= 1; got {n_directors!r}")
    return 3.0 / (2.0 * int(n_directors))


def isotropic_band_asymptote(n_directors: int) -> tuple[float, float]:
    """Closed-form ``(mean, std)`` of ``S`` under the isotropic null. **Asymptotic, not exact.**

    The second-moment identity in :func:`isotropic_eigenvalue_scale_sq` bounds the eigenvalue *scale*
    but says nothing about ``lambda_max`` on its own. This goes further, and the result is the part of
    this module that is new rather than re-derived.

    Using ``E[n_i n_j n_k n_l] = (d_ij d_kl + d_ik d_jl + d_il d_jk)/15`` for a uniform director, the
    covariance of ``t = n(x)n - I/3`` is exactly the isotropic traceless-symmetric form::

        Cov(t_ij, t_kl) = (2/15) [ (d_ik d_jl + d_il d_jk)/2 - d_ij d_kl / 3 ]

    so by the central limit theorem ``Q`` converges to the traceless Gaussian orthogonal ensemble with
    variance parameter ``v = (9/4)(1/N)(2/15) = 3/(10N)``, whose eigenvalue density is
    ``p ~ prod_{i<j}|l_i - l_j| exp(-sum l^2 / 2v) delta(sum l)``.

    Parameterise the traceless plane by ``l_k = R sqrt(2/3) cos(psi - 2 pi k/3)``, which makes
    ``sum l_k = 0`` identically and ``sum l_k^2 = R^2``, with ``psi`` the true polar angle. The
    Vandermonde factor collapses -- via the cubic discriminant and
    ``prod_k cos(psi - 2 pi k/3) = cos(3 psi)/4`` -- to ``R^3 |sin 3 psi| / sqrt(2)``, so with the area
    element ``R dR dpsi``::

        p(R, psi) ~ R^4 exp(-R^2 / 2v) * |sin 3 psi|

    **R and psi separate**, which is what makes this tractable. ``R = sqrt(v) chi_5``, and on the
    fundamental domain ``psi in [-pi/3, pi/3]`` the largest eigenvalue is ``R sqrt(2/3) cos psi``. Both
    remaining integrals are elementary -- ``E[chi_5] = 8 sqrt(2)/(3 sqrt(pi))``,
    ``E[cos psi] = 27/32``, ``E[cos^2 psi] = 29/40`` -- and after simplification::

        mean = 9 / (2 sqrt(10 pi N))                 ~= 0.802866 / sqrt(N)
        std  = sqrt( (29 pi - 81) / (40 pi N) )      ~= 0.283590 / sqrt(N)

    **Both invoke the CLT and are therefore asymptotic.** Measured relative excess of the sampled mean
    over this one: ``+7.6%`` at ``N = 4``, ``+1.8%`` at ``N = 64``, ``+0.35%`` at ``N = 4096`` -- an
    ``O(1/sqrt(N))`` approach. So this is the **oracle for the sampler and never a substitute for it**:
    at small ``N`` it understates the mean, which would inflate every ``z`` and bias a detector built on
    it toward firing.

    The substantive claim, independent of either constant: the null bias falls only as ``1/sqrt(N)``.
    That is why one fixed threshold cannot serve two population sizes.

    Args:
        n_directors: Population size ``N >= 1``.

    Returns:
        ``(mean, std)`` of ``S`` under the isotropic null, dimensionless.
    """
    n = int(n_directors)
    if n < 1:
        raise ValueError(f"n_directors must be >= 1; got {n_directors!r}")
    mean = 9.0 / (2.0 * math.sqrt(10.0 * math.pi * n))
    std = math.sqrt((29.0 * math.pi - 81.0) / (40.0 * math.pi * n))
    return mean, std


def sample_isotropic_directors(
    n_directors: int, rng: np.random.Generator, draws: int = 1
) -> npt.NDArray[np.float64]:
    """``draws`` independent populations of ``n_directors`` uniform unit vectors.

    Normalised Gaussians: the standard normal is spherically symmetric, so dividing by its norm gives
    an exactly uniform direction with no rejection step and no polar-angle bias. (Sampling the two
    spherical angles uniformly is the classic wrong answer -- it clusters at the poles.)

    Args:
        n_directors: Population size per draw.
        rng: An explicit generator, so a band is reproducible from its seed.
        draws: Number of independent populations.

    Returns:
        ``(n_directors, 3)`` when ``draws == 1``, else ``(draws, n_directors, 3)``.
    """
    n = int(n_directors)
    if n < 1:
        raise ValueError(f"n_directors must be >= 1; got {n_directors!r}")
    if int(draws) < 1:
        raise ValueError(f"draws must be >= 1; got {draws!r}")
    gaussian = rng.standard_normal((int(draws), n, 3))
    unit = gaussian / np.linalg.norm(gaussian, axis=2, keepdims=True)
    return unit[0] if int(draws) == 1 else unit


@dataclass(frozen=True, slots=True)
class NullBand:
    """What the order parameter would have been if the population were disordered.

    Attributes:
        n_directors: Population size this band applies to. A band is only valid at its own ``N``.
        mean: Sampled mean of ``S`` under the isotropic null -- the finite-N positive bias.
        std: Sampled standard deviation of ``S``.
        draws: Number of independent isotropic populations sampled.
        seed: Generator seed, so the band is reproducible.
        exact_eigenvalue_scale_sq: ``3/(2N)`` from :func:`isotropic_eigenvalue_scale_sq`, exact.
        sampled_eigenvalue_scale_sq: The same quantity, measured. Must match the exact one.
        asymptotic_mean: Closed-form mean from :func:`isotropic_band_asymptote`.
        asymptotic_std: Closed-form standard deviation from the same.
    """

    n_directors: int
    mean: float
    std: float
    draws: int
    seed: int
    exact_eigenvalue_scale_sq: float
    sampled_eigenvalue_scale_sq: float
    asymptotic_mean: float
    asymptotic_std: float

    @property
    def scale_relative_error(self) -> float:
        """``|sampled - exact| / exact`` for the eigenvalue scale. The sampler's own self-check.

        This compares a Monte-Carlo estimate against **algebra that holds at every N**, so a
        disagreement here is a defect in the sampler and never in the theory.
        """
        return abs(self.sampled_eigenvalue_scale_sq - self.exact_eigenvalue_scale_sq) / abs(
            self.exact_eigenvalue_scale_sq
        )

    @property
    def mean_relative_excess(self) -> float:
        """``(mean - asymptotic_mean) / asymptotic_mean``. Positive, and shrinking like ``1/sqrt(N)``.

        Not an error: the asymptote is a limit and the sample is finite, so the gap is the CLT
        correction. It is reported because its *size* is what justifies the detector using the sampled
        band rather than the closed form.
        """
        return (self.mean - self.asymptotic_mean) / self.asymptotic_mean

    def threshold(self, z_critical: float = ORDER_Z_CRITICAL) -> float:
        """The value of ``S`` this band requires before a population is called ordered."""
        return self.mean + float(z_critical) * self.std

    def is_decidable(self, z_critical: float = ORDER_Z_CRITICAL) -> bool:
        """Whether :meth:`threshold` is below 1, i.e. whether any population could clear it.

        See :func:`smallest_decidable_population`: below a certain ``N`` the bar exceeds the largest
        value ``S`` can take, and the test stops being strict and becomes vacuous.
        """
        return self.threshold(z_critical) < 1.0

    def as_json_obj(self) -> dict[str, Any]:
        """The JSON-able band."""
        return {
            "n_directors": self.n_directors,
            "mean": self.mean,
            "std": self.std,
            "draws": self.draws,
            "seed": self.seed,
            "exact_eigenvalue_scale_sq": self.exact_eigenvalue_scale_sq,
            "sampled_eigenvalue_scale_sq": self.sampled_eigenvalue_scale_sq,
            "scale_relative_error": self.scale_relative_error,
            "asymptotic_mean": self.asymptotic_mean,
            "asymptotic_std": self.asymptotic_std,
            "mean_relative_excess": self.mean_relative_excess,
            "threshold_at_z_critical": self.threshold(),
            "decidable": self.is_decidable(),
        }


#: Draws for a band. Enough to pin the standard deviation to a percent or so, which is far finer than
#: the decision margin it is used at; the band is reproducible because the seed is always explicit.
_DEFAULT_DRAWS: int = 4000
_DEFAULT_SEED: int = 20260731
#: Chunk the Monte-Carlo so a large population cannot allocate an unbounded array.
_MAX_CHUNK_ELEMENTS: int = 1 << 22


def isotropic_null_band(
    n_directors: int, draws: int = _DEFAULT_DRAWS, seed: int = _DEFAULT_SEED
) -> NullBand:
    """Sample the isotropic null distribution of ``S`` at ``n_directors``, with its cross-checks.

    Draws ``draws`` independent isotropic populations, measures ``S`` for each, and returns the
    ``(mean, std)`` together with **both** closed forms so the caller can see the sampler agreeing with
    algebra rather than being asked to trust it.

    Args:
        n_directors: Population size ``N``.
        draws: Independent isotropic populations to sample; ``>= 2`` for a standard deviation.
        seed: Generator seed.

    Returns:
        The :class:`NullBand`.
    """
    n = int(n_directors)
    if n < 1:
        raise ValueError(f"n_directors must be >= 1; got {n_directors!r}")
    n_draws = int(draws)
    if n_draws < 2:
        raise ValueError(f"draws must be >= 2 to estimate a standard deviation; got {draws!r}")

    rng = np.random.default_rng(int(seed))
    orders = np.empty(n_draws, dtype=np.float64)
    scales = np.empty(n_draws, dtype=np.float64)
    identity = np.eye(3)
    chunk = max(1, min(n_draws, _MAX_CHUNK_ELEMENTS // max(1, 3 * n)))
    start = 0
    while start < n_draws:
        stop = min(start + chunk, n_draws)
        block = sample_isotropic_directors(n, rng, draws=stop - start)
        block = block.reshape(stop - start, n, 3)
        # One Q per draw, then one batched symmetric eigen-decomposition for the whole block.
        mean_outer = np.einsum("dmi,dmj->dij", block, block) / n
        q = 1.5 * mean_outer - 0.5 * identity
        eigenvalues = np.linalg.eigvalsh(0.5 * (q + np.transpose(q, (0, 2, 1))))
        orders[start:stop] = eigenvalues[:, -1]
        scales[start:stop] = np.sum(eigenvalues**2, axis=1)
        start = stop

    asymptotic_mean, asymptotic_std = isotropic_band_asymptote(n)
    return NullBand(
        n_directors=n,
        mean=float(orders.mean()),
        std=float(orders.std(ddof=1)),
        draws=n_draws,
        seed=int(seed),
        exact_eigenvalue_scale_sq=isotropic_eigenvalue_scale_sq(n),
        sampled_eigenvalue_scale_sq=float(scales.mean()),
        asymptotic_mean=asymptotic_mean,
        asymptotic_std=asymptotic_std,
    )


def order_excess_z(order: float, band: NullBand) -> float:
    """``z = (S - mu_null) / sigma_null`` -- how far above disorder the population actually sits.

    This is the number a threshold should be applied to, and ``S`` is not. ``z`` has the population
    size divided out of it; ``S`` does not.

    Args:
        order: Measured ``S``.
        band: A band computed at the **same** population size.

    Returns:
        The excess in null standard deviations, dimensionless.
    """
    if not band.std > 0.0:
        raise ValueError(
            f"null band at N={band.n_directors} has std {band.std!r}; cannot score against a "
            "degenerate band"
        )
    value = float(order)
    if not math.isfinite(value):
        raise ValueError(f"order must be finite; got {order!r}")
    return (value - band.mean) / band.std


# ==================================================================================================
# Condensation: where an ordered AND dense patch separated from the background
# ==================================================================================================
@dataclass(frozen=True, slots=True)
class LocalFields:
    """Per-segment neighbourhood fields.

    Attributes:
        midpoints: ``(M, 3)`` segment midpoints [um].
        directors: ``(M, 3)`` segment directors.
        radius: Neighbourhood radius [um].
        density: ``(M,)`` neighbours within ``radius``, counting the segment itself.
        local_order: ``(M,)`` local ``S``. **``nan`` where the neighbourhood is too small to decide**,
            never ``0.0`` -- an undecidable measurement and a measurement of isotropy are different
            things, and ``0.0`` would be read as the second.
        alignment: ``(M,)`` ``|n . d_local|``. ``nan`` on the same neighbourhoods.
    """

    midpoints: npt.NDArray[np.float64]
    directors: npt.NDArray[np.float64]
    radius: float
    density: npt.NDArray[np.int64]
    local_order: npt.NDArray[np.float64]
    alignment: npt.NDArray[np.float64]

    @property
    def n_segments(self) -> int:
        """Population size."""
        return int(self.midpoints.shape[0])

    @property
    def n_undecidable(self) -> int:
        """Segments whose neighbourhood was too small to test. Reported, not hidden."""
        return int(np.count_nonzero(~np.isfinite(self.local_order)))


@dataclass(frozen=True, slots=True)
class CondensedBundle:
    """An ordered, dense patch that separated from the background.

    Attributes:
        members: Segment indices in this bundle.
        centre: ``(3,)`` mean midpoint of its members [um].
        axis: ``(3,)`` nematic director of its members.
        order: ``S`` of the members taken as one population.
        size: Member count.
    """

    members: npt.NDArray[np.int64]
    centre: npt.NDArray[np.float64]
    axis: npt.NDArray[np.float64]
    order: float
    size: int

    def as_json_obj(self) -> dict[str, Any]:
        """The JSON-able bundle. Members are listed so a figure can be checked against the data."""
        return {
            "size": self.size,
            "order": self.order,
            "centre_um": [float(v) for v in self.centre],
            "axis": [float(v) for v in self.axis],
            "members": [int(v) for v in self.members],
        }


def neighbourhood_radius(midpoints: Any, scale: float = NEIGHBOURHOOD_SCALE) -> float:
    """``scale`` times the median nearest-neighbour spacing of the midpoints [um].

    Deriving the radius from the configuration's own spacing is what makes the detector
    grid-invariant: a refined population has a smaller spacing and a proportionally smaller
    neighbourhood, so the neighbour counts -- and therefore the local null bands -- are unchanged.
    A radius fixed in um would silently change the measurement when the mesh changed.
    """
    from scipy.spatial import cKDTree  # local import: keeps module import cheap

    points = np.asarray(midpoints, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"midpoints must be (M, 3); got {points.shape}")
    if points.shape[0] < 2:
        raise ValueError("need at least 2 segments to define a spacing")
    if not float(scale) > 0.0:
        raise ValueError(f"scale must be positive; got {scale!r}")
    distances, _ = cKDTree(points).query(points, k=2)
    spacing = float(np.median(distances[:, 1]))
    if not spacing > 0.0:
        raise ValueError(
            "the median nearest-neighbour spacing is zero, so more than half the segments share a "
            "midpoint; no neighbourhood scale can be derived from this configuration"
        )
    return float(scale) * spacing


def local_fields(
    positions: Any,
    segments: Any,
    radius: float | None = None,
    scale: float = NEIGHBOURHOOD_SCALE,
) -> LocalFields:
    """Neighbourhood density and local nematic order for every segment.

    Args:
        positions: ``(N, 3)`` node positions [um].
        segments: ``(M, 2)`` node-index pairs.
        radius: Neighbourhood radius [um]; derived from the configuration's spacing when ``None``.
        scale: Multiplier used when ``radius`` is ``None``.

    Returns:
        The :class:`LocalFields`.
    """
    from scipy.spatial import cKDTree

    directors = segment_directors(positions, segments)
    midpoints = segment_midpoints(positions, segments)
    r = neighbourhood_radius(midpoints, scale) if radius is None else float(radius)
    if not r > 0.0:
        raise ValueError(f"radius must be positive; got {radius!r}")

    tree = cKDTree(midpoints)
    neighbours = tree.query_ball_point(midpoints, r)
    count = midpoints.shape[0]
    density = np.empty(count, dtype=np.int64)
    local_order = np.full(count, np.nan, dtype=np.float64)
    alignment = np.full(count, np.nan, dtype=np.float64)

    # Precompute every outer product once; a neighbourhood's Q is then a mean over a slice of these.
    outers = directors[:, :, None] * directors[:, None, :]
    identity = np.eye(3)
    for index in range(count):
        idx = neighbours[index]
        density[index] = len(idx)
        if len(idx) < MIN_NEIGHBOURS:
            continue
        q = 1.5 * outers[idx].mean(axis=0) - 0.5 * identity
        eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (q + q.T))
        local_order[index] = float(eigenvalues.max())
        axis = eigenvectors[:, int(np.argmax(eigenvalues))]
        alignment[index] = float(abs(directors[index] @ axis))

    return LocalFields(
        midpoints=midpoints,
        directors=directors,
        radius=r,
        density=density,
        local_order=local_order,
        alignment=alignment,
    )


def _local_bands(
    densities: npt.NDArray[np.int64], draws: int, seed: int
) -> dict[int, NullBand]:
    """One band per distinct decidable neighbour count present.

    The band depends on ``N``, so a neighbourhood of 8 and one of 40 are held to **different** bars.
    That is the substance of scoring against a band rather than a threshold: the same local ``S`` means
    something different at two neighbourhood sizes, and a single cut-off would systematically flag the
    small ones.
    """
    bands: dict[int, NullBand] = {}
    for value in np.unique(densities):
        count = int(value)
        if count < MIN_NEIGHBOURS:
            continue
        bands[count] = isotropic_null_band(count, draws=draws, seed=seed + count)
    return bands


def condensed_bundles(
    fields: LocalFields,
    z_critical: float = ORDER_Z_CRITICAL,
    min_size: int = MIN_NEIGHBOURS,
    draws: int = 1000,
    seed: int = 512,
) -> list[CondensedBundle]:
    """Flag segments whose local order beats their **own** null band, then cluster them.

    A segment is condensed when both hold:

    1. ``local_order > mu_null(m) + z_critical * sigma_null(m)`` for its own neighbour count ``m`` --
       it sits in an ordered patch, judged at the right population size; and
    2. ``alignment > ALIGNMENT_COSINE_MIN`` -- its own director lies in the patch's cone, so it is a
       *member* of the bundle rather than a bystander next to one.

    Requiring both is what stops a merely **dense** patch being called a bundle. Density enters only
    through the neighbourhood: crowding raises ``m``, which raises the bar rather than lowering it.

    Flagged segments within one neighbourhood radius of each other are grouped by connectivity, and
    groups below ``min_size`` are dropped. Bundles come back largest first.

    Args:
        fields: From :func:`local_fields`.
        z_critical: Significance threshold.
        min_size: Smallest reportable bundle.
        draws: Draws per local band.
        seed: Base seed for the local bands.

    Returns:
        The recovered bundles, largest first. An empty list is a real answer: no patch condensed.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    bands = _local_bands(fields.density, draws=draws, seed=seed)
    flagged = np.zeros(fields.n_segments, dtype=bool)
    for index in range(fields.n_segments):
        band = bands.get(int(fields.density[index]))
        if band is None or not math.isfinite(fields.local_order[index]):
            continue
        if not band.std > 0.0:
            continue
        if not band.is_decidable(z_critical):
            # The bar for this neighbourhood size exceeds the largest value S can take, so no
            # arrangement of these directors could clear it. Skipping is the honest action: flagging
            # would be impossible and *not* flagging must not be read as evidence of disorder.
            # MIN_NEIGHBOURS already excludes this at the default z; this catches a raised one.
            continue
        ordered = fields.local_order[index] > band.threshold(z_critical)
        aligned = fields.alignment[index] > ALIGNMENT_COSINE_MIN
        if ordered and aligned:
            flagged[index] = True

    members = np.flatnonzero(flagged)
    if members.size == 0:
        return []

    points = fields.midpoints[members]
    pairs = cKDTree(points).query_pairs(fields.radius, output_type="ndarray")
    size = members.size
    if pairs.size:
        rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
        cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
        adjacency = coo_matrix((np.ones(rows.size), (rows, cols)), shape=(size, size))
    else:
        adjacency = coo_matrix((size, size))
    n_groups, labels = connected_components(adjacency, directed=False)

    bundles: list[CondensedBundle] = []
    for group in range(n_groups):
        group_members = members[labels == group]
        if group_members.size < int(min_size):
            continue
        order, axis = nematic_order_and_director(fields.directors[group_members])
        bundles.append(
            CondensedBundle(
                members=group_members,
                centre=fields.midpoints[group_members].mean(axis=0),
                axis=axis,
                order=order,
                size=int(group_members.size),
            )
        )
    bundles.sort(key=lambda bundle: bundle.size, reverse=True)
    return bundles


@dataclass(frozen=True, slots=True)
class EmergenceObservation:
    """The observers' verdict for one filament population.

    **Every field is a geometric or statistical quantity.** There is no force, torque, stress or energy
    here, and there is no path by which one could be added without a test failing -- see
    ``tests/viz/test_emergence.py::test_the_observation_carries_no_force_field``.

    Attributes:
        n_segments: Population size.
        order: Global ``S``, dimensionless.
        director: ``(3,)`` global nematic axis.
        band: The null band ``order`` was scored against.
        z: ``(order - band.mean) / band.std``.
        is_ordered: Whether ``z`` cleared the threshold.
        radius: Neighbourhood radius used for the local fields [um].
        bundles: Recovered bundles, largest first.
        n_undecidable: Segments whose neighbourhood was too small to test.
    """

    n_segments: int
    order: float
    director: npt.NDArray[np.float64]
    band: NullBand
    z: float
    is_ordered: bool
    radius: float
    bundles: list[CondensedBundle] = field(default_factory=list)
    n_undecidable: int = 0

    @property
    def n_bundles(self) -> int:
        """How many patches condensed. Measured, never the segment count handed in."""
        return len(self.bundles)

    @property
    def condensed_fraction(self) -> float:
        """Fraction of segments inside some bundle."""
        if self.n_segments == 0:
            return 0.0
        return sum(bundle.size for bundle in self.bundles) / self.n_segments

    def as_json_obj(self) -> dict[str, Any]:
        """The JSON-able observation."""
        return {
            "n_segments": self.n_segments,
            "order": self.order,
            "director": [float(v) for v in self.director],
            "z": self.z,
            "is_ordered": self.is_ordered,
            "threshold": self.band.threshold(),
            "radius_um": self.radius,
            "n_bundles": self.n_bundles,
            "condensed_fraction": self.condensed_fraction,
            "n_undecidable": self.n_undecidable,
            "band": self.band.as_json_obj(),
            "bundles": [bundle.as_json_obj() for bundle in self.bundles],
        }


def observe_emergence(
    positions: Any,
    segments: Any,
    z_critical: float = ORDER_Z_CRITICAL,
    scale: float = NEIGHBOURHOOD_SCALE,
    draws: int = _DEFAULT_DRAWS,
    seed: int = _DEFAULT_SEED,
    local_draws: int = 1000,
) -> EmergenceObservation:
    """Measure a filament population's order and find where it condensed.

    Reads geometry and nothing else -- no owner name, no material card, no connector state. The
    detector cannot be told what to find.

    Args:
        positions: ``(N, 3)`` node positions [um].
        segments: ``(M, 2)`` node-index pairs.
        z_critical: Significance threshold.
        scale: Neighbourhood multiplier.
        draws: Draws for the global band.
        seed: Seed for the global band.
        local_draws: Draws per local band.

    Returns:
        The :class:`EmergenceObservation`.

    Note:
        A flat result -- ``order`` inside the band, no bundle recovered -- is a **real answer** and is
        reported as one. It is never a reason to lower ``z_critical``.
    """
    directors = segment_directors(positions, segments)
    count = int(directors.shape[0])
    order, director = nematic_order_and_director(directors)
    band = isotropic_null_band(count, draws=draws, seed=seed)
    z = order_excess_z(order, band)
    fields = local_fields(positions, segments, scale=scale)
    bundles = condensed_bundles(fields, z_critical=z_critical, draws=local_draws)
    return EmergenceObservation(
        n_segments=count,
        order=order,
        director=director,
        band=band,
        z=z,
        is_ordered=bool(z > z_critical),
        radius=fields.radius,
        bundles=bundles,
        n_undecidable=fields.n_undecidable,
    )
