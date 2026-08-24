"""Structural identifiability: which parameter combinations an observation can see, and which not.

This is a different question from precision, and confusing the two is how a project ends up quoting
a number for something it never measured.  Precision asks *how well* a parameter is determined by
noisy data.  Structural identifiability asks whether the forward model's output changes **at all**
when the parameter moves.  If it does not, no amount of data helps, the error bar is infinite
rather than large, and the honest output is a refusal naming the degenerate combination --- not a
value.

The test is the null space of the Jacobian.  ``J v = 0`` means moving along ``v`` produces exactly
no change in the prediction, so ``v`` is invisible to this observation.  Working from ``J`` rather
than from the Fisher matrix ``F = J^T C^-1 J`` is deliberate: ``F`` squares the singular values, so
a direction whose singular value sits at ``1e-9`` relative appears at ``1e-18`` in ``F`` --- at the
edge of float64 --- and the rank decision becomes a coin flip.  The SVD of ``J`` keeps the full
dynamic range.

**Log coordinates are usually the right place to ask.**  Multiplicative degeneracies --- the ones
that actually occur in physics --- become *linear* in logs and therefore exactly detectable.  A
model that depends on its parameters only through the product ``a*b`` has, in log coordinates, the
exact null direction ``(1, -1)/sqrt(2)``: raise ``log a``, lower ``log b``, and the prediction does
not move.  In linear coordinates the same degeneracy is a curved manifold and the local null
direction at ``(a, b)`` is the parameter-dependent ``(b, -a)/norm``, which is correct but tells you
less.

That toy is the same shape as the real degeneracy in this project's first target.  A passive
fluctuation spectrum ``P(q) = A / (kappa*q^4 + sigma*q^2)`` with an unknown amplitude calibration
``A`` is unchanged under ``(A, kappa, sigma) -> (f*A, f*kappa, f*sigma)``, an exact one-dimensional
degeneracy whose null direction in log coordinates is ``(1, 1, 1)/sqrt(3)``.  Both are built as
tests in ``tests/infer/test_identifiability.py`` rather than discovered later as a surprise.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
from numpy.typing import NDArray

from .fisher import describe_direction
from .numdiff import jacobian as numeric_jacobian
from .refusal import Unidentifiable

__all__ = [
    "DEFAULT_NULL_RTOL",
    "DegeneracyReport",
    "NullSpace",
    "TradeOff",
    "degeneracy_report",
    "identifiability_of",
    "null_directions",
    "null_directions_of_model",
    "require_identifiable",
]

#: Relative singular-value cutoff for calling a direction null.
#:
#: Chosen against the accuracy of the Jacobian, not against float64 epsilon.  A central-difference
#: Jacobian carries a relative error near ``eps^(2/3) ~ 4e-11``, so an exactly null direction shows
#: up with a residual singular value around ``1e-10`` of the largest.  ``1e-7`` sits safely above
#: that and safely below any real physical stiffness ratio worth reporting.  Supply an analytic
#: Jacobian and this can be tightened by orders of magnitude.
DEFAULT_NULL_RTOL: float = 1e-7

#: Fraction of a queried direction that may lie in the null space before it is refused.
DEFAULT_MAX_NULL_OVERLAP: float = 1e-3


def _default_names(n: int, names: Sequence[str] | None) -> tuple[str, ...]:
    if names is None:
        return tuple(f"theta{i}" for i in range(n))
    if len(names) != n:
        raise ValueError(f"expected {n} parameter names, got {len(names)}")
    return tuple(str(x) for x in names)


def _canonical_sign(vector: NDArray) -> NDArray:
    idx = int(np.argmax(np.abs(vector)))
    return -vector if vector[idx] < 0 else vector


@dataclass(frozen=True)
class TradeOff:
    """A pairwise exchange rate along a degenerate direction.

    ``slope`` answers: to stay invisible to this observation, if ``first`` moves by one unit, how
    far must ``second`` move?  In log coordinates a slope of ``-1`` is the signature of a product
    degeneracy and a slope of ``+1`` of a ratio degeneracy.
    """

    first: str
    second: str
    slope: float

    def __str__(self) -> str:
        return (
            f"{self.first} trades against {self.second} at {self.slope:+.4f} "
            f"(move {self.first} by +1, compensate with {self.slope:+.4f} of {self.second})"
        )


@dataclass(frozen=True)
class NullSpace:
    """Directions in parameter space that the observation cannot see at all.

    Attributes:
        vectors: ``(n_parameters, n_null)``, orthonormal columns, signs canonicalized.
        singular_values: All singular values of the Jacobian, descending.
        null_singular_values: The subset below the tolerance.
        tolerance: Absolute singular-value cutoff actually used.
        rtol: The relative cutoff it came from.
        names: Parameter names.
    """

    vectors: NDArray
    singular_values: NDArray
    null_singular_values: NDArray
    tolerance: float
    rtol: float
    names: tuple[str, ...]

    @property
    def dimension(self) -> int:
        return int(self.vectors.shape[1])

    @property
    def is_empty(self) -> bool:
        return self.dimension == 0

    @property
    def rank(self) -> int:
        """Number of independent combinations the observation *can* see."""
        return len(self.names) - self.dimension

    def projector(self) -> NDArray:
        """``N N^T`` --- projects a direction onto the invisible subspace."""
        if self.is_empty:
            return np.zeros((len(self.names), len(self.names)), dtype=float)
        return self.vectors @ self.vectors.T

    def overlap(self, direction: NDArray) -> float:
        """Fraction of ``direction``'s length that lies in the null space, in ``[0, 1]``."""
        d = np.asarray(direction, dtype=float).ravel()
        if d.size != len(self.names):
            raise ValueError(f"direction must have {len(self.names)} components, got {d.size}")
        norm = float(np.linalg.norm(d))
        if norm == 0.0:
            raise ValueError("the zero vector has no direction")
        return float(np.linalg.norm(self.projector() @ d) / norm)

    def parameter_overlaps(self) -> dict[str, float]:
        """Per-parameter null overlap: how much of each bare parameter is invisible."""
        eye = np.eye(len(self.names))
        return {name: self.overlap(eye[:, i]) for i, name in enumerate(self.names)}

    def unconstrained_monomial(self, index: int) -> str:
        """The monomial that null direction ``index`` leaves free, read in log coordinates.

        If the analysis ran in ``log`` parameters, a null direction ``v`` means the product
        ``prod_i theta_i^{v_i}`` can be changed freely without changing any prediction.  Everything
        orthogonal to ``v`` is what the data actually constrains.  Rendering this makes the result
        quotable: "we did not measure kappa, we measured kappa/sigma" is a scientific statement;
        "kappa is unidentifiable" is only half of one.
        """
        v = np.asarray(self.vectors[:, index], dtype=float)
        terms = [
            f"{name}^{w:+.4f}"
            for name, w in zip(self.names, v, strict=True)
            if abs(w) >= 1e-3
        ]
        return " * ".join(terms) if terms else "(none above threshold)"

    def report(self) -> str:
        if self.is_empty:
            return (
                f"No null directions at rtol={self.rtol:.3g}: every parameter combination changes "
                f"the prediction. Smallest singular value {float(self.singular_values[-1]):.6g} "
                f"against largest {float(self.singular_values[0]):.6g}."
            )
        lines = [
            f"{self.dimension} null direction(s) of {len(self.names)} at rtol={self.rtol:.3g}; "
            f"the observation constrains only {self.rank} combination(s):"
        ]
        for i in range(self.dimension):
            lines.append(
                f"  null[{i}] {describe_direction(self.vectors[:, i], self.names)} "
                f"(singular value {float(self.null_singular_values[i]):.6g}); "
                f"free monomial in log coordinates: {self.unconstrained_monomial(i)}"
            )
        return "\n".join(lines)


def null_directions(
    jacobian: NDArray,
    *,
    names: Sequence[str] | None = None,
    rtol: float = DEFAULT_NULL_RTOL,
    variance: float | NDArray | None = None,
) -> NullSpace:
    """Null space of a forward-model Jacobian, by SVD.

    Args:
        jacobian: ``(n_observations, n_parameters)``.
        names: Parameter names.
        rtol: Singular values below ``rtol * s_max`` are treated as zero.
        variance: Optional per-observation noise variance.  When given, the Jacobian is whitened by
            ``1/sqrt(v)`` first, which makes the singular values dimensionless in units of noise ---
            worth doing when observations have wildly different error bars, since an unwhitened SVD
            lets a loud, uninformative channel dominate the rank decision.
    """
    j = np.asarray(jacobian, dtype=float)
    if j.ndim != 2:
        raise ValueError(f"jacobian must be 2-D, got shape {j.shape}")
    n_par = j.shape[1]
    if variance is not None:
        v = np.atleast_1d(np.asarray(variance, dtype=float))
        if v.size == 1:
            v = np.full(j.shape[0], float(v[0]))
        if v.shape != (j.shape[0],):
            raise ValueError(f"variance must be scalar or length {j.shape[0]}")
        if not np.all(v > 0.0):
            raise ValueError("every variance entry must be strictly positive")
        j = j / np.sqrt(v)[:, None]

    _, s_values, vt = np.linalg.svd(j, full_matrices=True)
    # full_matrices=True so vt has all n_par rows even when the model has fewer observations than
    # parameters. That under-determined case is exactly when a null space exists for a trivial
    # reason, and truncating it away would hide the most obvious failure of all.
    s_full = np.zeros(n_par, dtype=float)
    s_full[: s_values.size] = s_values
    s_max = float(s_full[0]) if s_full.size else 0.0
    tol = rtol * s_max
    null_mask = s_full <= tol
    vectors = np.array(
        [_canonical_sign(np.asarray(vt[i], dtype=float)) for i in np.nonzero(null_mask)[0]]
    ).T
    if vectors.size == 0:
        vectors = np.zeros((n_par, 0), dtype=float)
    return NullSpace(
        vectors=vectors,
        singular_values=s_full,
        null_singular_values=s_full[null_mask],
        tolerance=float(tol),
        rtol=float(rtol),
        names=_default_names(n_par, names),
    )


def null_directions_of_model(
    model: Callable[[NDArray], NDArray],
    theta: NDArray,
    *,
    names: Sequence[str] | None = None,
    rtol: float = DEFAULT_NULL_RTOL,
    variance: float | NDArray | None = None,
    jacobian: Callable[[NDArray], NDArray] | None = None,
) -> NullSpace:
    """Null space of a generic forward model at ``theta``.

    This is a *local* statement: it holds at ``theta``.  A degeneracy that is exact everywhere (a
    model that only ever sees ``a*b``) and one that happens to be degenerate at a single point look
    identical from here.  Evaluate at several points to tell them apart --- if the null direction
    rotates with ``theta``, the degeneracy is local; if it is constant, it is structural.
    """
    theta = np.asarray(theta, dtype=float)
    j = (
        np.asarray(jacobian(theta), dtype=float)
        if jacobian is not None
        else numeric_jacobian(model, theta)
    )
    return null_directions(j, names=names, rtol=rtol, variance=variance)


def _as_direction(query: int | str | Sequence[float] | NDArray, names: Sequence[str]) -> NDArray:
    """Turn a parameter index, a parameter name, or a vector into a unit direction."""
    n = len(names)
    if isinstance(query, str):
        if query not in names:
            raise ValueError(f"unknown parameter {query!r}; known: {list(names)}")
        vec = np.zeros(n)
        vec[list(names).index(query)] = 1.0
        return vec
    if isinstance(query, (int, np.integer)):
        vec = np.zeros(n)
        vec[int(query)] = 1.0
        return vec
    vec = np.asarray(query, dtype=float).ravel()
    if vec.size != n:
        raise ValueError(f"direction must have {n} components, got {vec.size}")
    norm = float(np.linalg.norm(vec))
    if norm == 0.0:
        raise ValueError("the zero vector has no direction")
    return vec / norm


def identifiability_of(
    query: int | str | Sequence[float] | NDArray,
    null_space: NullSpace,
    *,
    max_null_overlap: float = DEFAULT_MAX_NULL_OVERLAP,
) -> Unidentifiable | None:
    """``None`` when the queried quantity is measurable; a typed refusal when it is not.

    The metric is the fraction of the queried direction that lies in the null space.  It is 0 for a
    fully visible direction and 1 for a fully invisible one, and it is a *fraction* rather than a
    distance so the answer does not depend on how the query was scaled.

    The refusal carries the offending null direction and the free monomial, because the useful
    result is almost never "cannot measure it" on its own --- it is "cannot measure it separately
    from that".
    """
    direction = _as_direction(query, null_space.names)
    if null_space.is_empty:
        return None
    overlap = null_space.overlap(direction)
    if overlap <= max_null_overlap:
        return None
    worst = int(np.argmax(np.abs(null_space.vectors.T @ direction)))
    return Unidentifiable(
        reason=(
            f"[{describe_direction(direction, null_space.names)}] lies "
            f"{100.0 * overlap:.4g}% inside the null space of this observation. A number reported "
            "for it would be an artifact of the prior or the optimizer's starting point, not a "
            "measurement."
        ),
        offending_input={
            "query": direction.tolist(),
            "null_overlap": overlap,
            "max_null_overlap": float(max_null_overlap),
        },
        lift_requires=(
            "Add an observation whose Jacobian has support along this direction, fix the "
            "degenerate combination from independent evidence, or report the free monomial "
            f"[{null_space.unconstrained_monomial(worst)}] as unconstrained and quote only the "
            "combinations orthogonal to it."
        ),
        detail={
            "null_dimension": null_space.dimension,
            "rtol": null_space.rtol,
            "free_monomial": null_space.unconstrained_monomial(worst),
        },
        direction=tuple(float(x) for x in null_space.vectors[:, worst]),
        parameter_names=null_space.names,
    )


def require_identifiable(
    query: int | str | Sequence[float] | NDArray,
    null_space: NullSpace,
    *,
    max_null_overlap: float = DEFAULT_MAX_NULL_OVERLAP,
) -> None:
    """Raise :class:`~aleph.infer.refusal.UnidentifiableDirectionError` if ``query`` is invisible.

    The raising form exists for call sites under an optimizer or inside a pipeline, where returning
    a refusal object would just be ignored by the next line.
    """
    refusal = identifiability_of(query, null_space, max_null_overlap=max_null_overlap)
    if refusal is not None:
        refusal.raise_()


@dataclass(frozen=True)
class DegeneracyReport:
    """The geometry of what this observation cannot separate.

    Attributes:
        null_space: The invisible subspace.
        trade_offs: Pairwise exchange rates, one list per null direction.
        unmeasurable: Parameters whose bare direction is substantially inside the null space.
        max_null_overlap: The cutoff used to build ``unmeasurable``.
    """

    null_space: NullSpace
    trade_offs: tuple[tuple[TradeOff, ...], ...]
    unmeasurable: tuple[str, ...]
    max_null_overlap: float

    @property
    def degenerate(self) -> bool:
        return not self.null_space.is_empty

    def render(self) -> str:
        lines = [self.null_space.report()]
        for i, group in enumerate(self.trade_offs):
            if not group:
                continue
            lines.append(f"  trade-offs along null[{i}]:")
            lines.extend(f"    {t}" for t in group)
        if self.unmeasurable:
            lines.append(
                "  individually unmeasurable parameters: " + ", ".join(self.unmeasurable)
            )
        else:
            lines.append("  every individual parameter is measurable at this threshold")
        return "\n".join(lines)


def degeneracy_report(
    null_space: NullSpace,
    *,
    max_null_overlap: float = DEFAULT_MAX_NULL_OVERLAP,
    min_loading: float = 1e-3,
) -> DegeneracyReport:
    """Describe which parameters trade off against which, and at what exchange rate.

    For each null direction the parameters with meaningful loading are paired against the
    largest-loading one, since that is the pivot a reader will hold fixed when reasoning about the
    degeneracy.  The slope is ``v_j / v_pivot``: staying inside the null space means moving along
    ``v``, so a unit move in the pivot must be accompanied by exactly that much of parameter ``j``.
    """
    groups: list[tuple[TradeOff, ...]] = []
    for i in range(null_space.dimension):
        v = np.asarray(null_space.vectors[:, i], dtype=float)
        pivot = int(np.argmax(np.abs(v)))
        pairs: list[TradeOff] = []
        for j in range(v.size):
            if j == pivot or abs(v[j]) < min_loading:
                continue
            pairs.append(
                TradeOff(
                    first=null_space.names[pivot],
                    second=null_space.names[j],
                    slope=float(v[j] / v[pivot]),
                )
            )
        groups.append(tuple(pairs))

    overlaps = null_space.parameter_overlaps()
    unmeasurable = tuple(
        name for name, value in overlaps.items() if value > max_null_overlap
    )
    return DegeneracyReport(
        null_space=null_space,
        trade_offs=tuple(groups),
        unmeasurable=unmeasurable,
        max_null_overlap=float(max_null_overlap),
    )
