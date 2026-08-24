"""Fisher information from the forward-model Jacobian, and what it says cannot be measured.

For a Gaussian observation ``d = m(theta) + n`` with noise covariance ``C`` that does not depend on
``theta``, the Fisher information matrix is exactly

``F = J^T C^-1 J``     with     ``J_ij = d m_i / d theta_j``

That is a closed form, not an approximation, and it is the positive control this module is required
to reproduce to machine precision: for a linear model ``m(theta) = A theta`` the Jacobian *is*
``A``, so ``F = A^T C^-1 A`` with no error term at all.  If that fails, every uncertainty this
project ever quotes is meaningless, so it is tested first and tested exactly.

**Why the eigenbasis is the useful output, not the diagonal.**  Reporting ``sqrt(F^-1_ii)`` per
parameter hides the structure that actually governs the measurement.  A two-parameter fit can have
excellent precision on one combination and none at all on another, and the marginal error bars will
show two mediocre numbers with no hint that the pair is nearly degenerate.  The eigendecomposition
of ``F`` names the combinations directly: large eigenvalue means *stiff*, the data pins that
combination down; small eigenvalue means *sloppy*, the data barely constrains it.  The eigenvalue
spectrum of a sloppy model typically spans many decades, and the honest report of such a measurement
is "we measured these two combinations and could not see those three".

**Cramer-Rao.**  ``Cov(theta_hat) >= F^-1`` for any unbiased estimator, in the positive-semidefinite
sense.  It is a *lower bound*: a real estimator does no better, and usually does worse.  Quoting it
as "the error bar" is optimistic by an unknown factor, so :class:`CramerRaoBound` says so in its own
docstring, and the coverage campaign in :mod:`aleph.infer.posterior` measures the difference.

A singular or near-singular ``F`` does not get a large number returned for it.  It gets
:class:`~aleph.infer.refusal.Unidentifiable`, naming the offending direction.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy import linalg

from .numdiff import jacobian as numeric_jacobian
from .refusal import Unidentifiable

__all__ = [
    "CramerRaoBound",
    "FisherSpectrum",
    "IdentifiabilitySplit",
    "ParameterDirection",
    "cramer_rao_bound",
    "describe_direction",
    "fisher_from_model",
    "fisher_matrix",
    "fisher_spectrum",
    "identifiable_directions",
]

#: Default split between stiff and sloppy: an eigenvalue below this fraction of the largest one is
#: called sloppy.  1e-8 in eigenvalue is 1e-4 in the parameter uncertainty it implies, and it sits
#: comfortably above float64 noise in a Jacobian built from central differences.
DEFAULT_STIFFNESS_THRESHOLD: float = 1e-8


def _default_names(n: int, names: Sequence[str] | None) -> tuple[str, ...]:
    if names is None:
        return tuple(f"theta{i}" for i in range(n))
    if len(names) != n:
        raise ValueError(f"expected {n} parameter names, got {len(names)}")
    return tuple(str(x) for x in names)


def _canonical_sign(vector: NDArray) -> NDArray:
    """Fix the arbitrary sign of an eigenvector so reports and tests are reproducible.

    An eigenvector is defined up to sign, and numpy's choice depends on LAPACK's internals.  The
    convention here is that the component of largest magnitude is positive; ties break toward the
    lowest index, which makes the result independent of the ordering of near-equal components.
    """
    idx = int(np.argmax(np.abs(vector)))
    return -vector if vector[idx] < 0 else vector


def describe_direction(
    vector: NDArray, names: Sequence[str], *, max_terms: int = 4, min_weight: float = 0.05
) -> str:
    """Render a parameter-space direction as a readable combination.

    ``[0.707, -0.707]`` with names ``("log_kappa", "log_sigma")`` becomes
    ``"+0.707 log_kappa -0.707 log_sigma"``, which is a statement someone can act on --- in log
    coordinates that direction is the ratio ``kappa/sigma``.  Components below ``min_weight`` are
    dropped so a direction that is essentially one parameter reads as that parameter.
    """
    vector = np.asarray(vector, dtype=float)
    order = np.argsort(-np.abs(vector))
    parts: list[str] = []
    for j in order[:max_terms]:
        w = float(vector[j])
        if abs(w) < min_weight:
            continue
        parts.append(f"{w:+.3f} {names[j]}")
    if not parts:
        return "(no component above the reporting threshold)"
    return " ".join(parts)


@dataclass(frozen=True)
class ParameterDirection:
    """One eigendirection of the Fisher matrix, named.

    Attributes:
        index: Rank in the spectrum, 0 = stiffest.
        eigenvalue: Fisher eigenvalue along this direction.
        vector: Unit vector in parameter space, sign-canonicalized.
        names: Parameter names, in the order the vector's components appear.
        stiff: True when the eigenvalue is above the stiffness threshold.
        threshold: The threshold used, kept so the label can be re-derived.
    """

    index: int
    eigenvalue: float
    vector: NDArray
    names: tuple[str, ...]
    stiff: bool
    threshold: float

    @property
    def label(self) -> str:
        return "stiff" if self.stiff else "sloppy"

    @property
    def uncertainty(self) -> float:
        """``1/sqrt(eigenvalue)``: the standard deviation along this direction alone.

        This is the *conditional* uncertainty --- what you would get if every other direction were
        known exactly.  It is not the marginal error bar; that is
        :attr:`CramerRaoBound.standard_deviation`, and it is larger whenever directions correlate.
        """
        if self.eigenvalue <= 0.0:
            return float("inf")
        return float(1.0 / np.sqrt(self.eigenvalue))

    @property
    def description(self) -> str:
        return describe_direction(self.vector, self.names)

    def __str__(self) -> str:
        return (
            f"#{self.index} {self.label}: {self.description} "
            f"(eigenvalue {self.eigenvalue:.6g}, sigma along it {self.uncertainty:.6g})"
        )


@dataclass(frozen=True)
class FisherSpectrum:
    """Eigendecomposition of a Fisher matrix, with the directions named.

    Attributes:
        matrix: The Fisher matrix itself.
        eigenvalues: Descending, so index 0 is the stiffest direction.
        eigenvectors: Columns matching ``eigenvalues``, sign-canonicalized.
        names: Parameter names.
        threshold: Relative eigenvalue below which a direction is called sloppy.
    """

    matrix: NDArray
    eigenvalues: NDArray
    eigenvectors: NDArray
    names: tuple[str, ...]
    threshold: float

    @property
    def condition_number(self) -> float:
        """``lambda_max / lambda_min``.  Decades of this are the signature of a sloppy model."""
        lo = float(self.eigenvalues[-1])
        if lo <= 0.0:
            return float("inf")
        return float(self.eigenvalues[0]) / lo

    def direction(self, index: int) -> ParameterDirection:
        lam = float(self.eigenvalues[index])
        cutoff = self.threshold * float(self.eigenvalues[0])
        return ParameterDirection(
            index=int(index),
            eigenvalue=lam,
            vector=_canonical_sign(np.asarray(self.eigenvectors[:, index], dtype=float)),
            names=self.names,
            stiff=bool(lam > cutoff),
            threshold=float(self.threshold),
        )

    @property
    def directions(self) -> tuple[ParameterDirection, ...]:
        return tuple(self.direction(i) for i in range(self.eigenvalues.size))

    @property
    def stiffest(self) -> ParameterDirection:
        return self.direction(0)

    @property
    def sloppiest(self) -> ParameterDirection:
        return self.direction(int(self.eigenvalues.size) - 1)

    def report(self) -> str:
        lines = [
            f"Fisher spectrum over {list(self.names)}: "
            f"condition number {self.condition_number:.6g}"
        ]
        lines.extend(f"  {d}" for d in self.directions)
        return "\n".join(lines)


def _precision_apply(
    residual_operator: NDArray, variance: NDArray | None, covariance: NDArray | None
) -> NDArray:
    """Return ``C^-1 X`` for ``X = residual_operator``, given either a diagonal or a full ``C``."""
    if (variance is None) == (covariance is None):
        raise ValueError("supply exactly one of variance= or covariance=")
    if variance is not None:
        v = np.atleast_1d(np.asarray(variance, dtype=float))
        n = residual_operator.shape[0]
        if v.size == 1:
            v = np.full(n, float(v[0]))
        if v.shape != (n,):
            raise ValueError(f"variance must be scalar or length {n}, got shape {v.shape}")
        if not np.all(v > 0.0):
            raise ValueError("every variance entry must be strictly positive")
        return residual_operator / v[:, None]
    cov = np.asarray(covariance, dtype=float)
    n = residual_operator.shape[0]
    if cov.shape != (n, n):
        raise ValueError(f"covariance must have shape ({n}, {n}), got {cov.shape}")
    # Cholesky both solves and asserts positive definiteness; a covariance that is not positive
    # definite is not a covariance, and silently pseudo-inverting it would hide that.
    try:
        chol = linalg.cho_factor(cov, lower=True, check_finite=True)
    except linalg.LinAlgError as exc:
        raise ValueError("covariance is not positive definite") from exc
    return linalg.cho_solve(chol, residual_operator, check_finite=True)


def fisher_matrix(
    jacobian: NDArray,
    *,
    variance: float | NDArray | None = None,
    covariance: NDArray | None = None,
) -> NDArray:
    """``F = J^T C^-1 J`` --- the closed form, evaluated as written.

    Args:
        jacobian: ``(n_observations, n_parameters)``.
        variance: Diagonal noise, scalar or per-observation.  Exclusive with ``covariance``.
        covariance: Full noise covariance ``(n_obs, n_obs)``.  Exclusive with ``variance``.

    Returns:
        Symmetric ``(n_parameters, n_parameters)`` matrix.  Symmetrized explicitly: ``J^T C^-1 J``
        is symmetric in exact arithmetic, and the asymmetry left by floating point is noise that
        would otherwise leak complex parts into the eigendecomposition.
    """
    j = np.asarray(jacobian, dtype=float)
    if j.ndim != 2:
        raise ValueError(f"jacobian must be 2-D, got shape {j.shape}")
    precision_j = _precision_apply(j, variance, covariance)
    f = j.T @ precision_j
    return 0.5 * (f + f.T)


def fisher_from_model(
    model: Callable[[NDArray], NDArray],
    theta: NDArray,
    *,
    variance: float | NDArray | None = None,
    covariance: NDArray | None = None,
    jacobian: Callable[[NDArray], NDArray] | None = None,
) -> NDArray:
    """Fisher matrix of a generic forward model at ``theta``.

    The Jacobian is analytic when supplied and central-differenced otherwise.  Note the assumption
    this whole construction rests on: the noise covariance does not depend on ``theta``.  When it
    does --- as it does for a Gamma-distributed spectrum, where the variance scales with the
    predicted power --- the ``J^T C^-1 J`` form is missing the term from the parameter dependence of
    ``C``, and the resulting information is an underestimate.
    """
    theta = np.asarray(theta, dtype=float)
    j = (
        np.asarray(jacobian(theta), dtype=float)
        if jacobian is not None
        else numeric_jacobian(model, theta)
    )
    return fisher_matrix(j, variance=variance, covariance=covariance)


def fisher_spectrum(
    matrix: NDArray,
    *,
    names: Sequence[str] | None = None,
    threshold: float = DEFAULT_STIFFNESS_THRESHOLD,
) -> FisherSpectrum:
    """Eigendecompose a Fisher matrix, sorted stiffest first, signs canonicalized."""
    f = np.asarray(matrix, dtype=float)
    if f.ndim != 2 or f.shape[0] != f.shape[1]:
        raise ValueError(f"Fisher matrix must be square, got shape {f.shape}")
    sym = 0.5 * (f + f.T)
    # eigh, not eig: the matrix is symmetric by construction and eigh guarantees real eigenvalues
    # and an orthonormal basis, which the reports below assume.
    values, vectors = np.linalg.eigh(sym)
    order = np.argsort(-values)
    return FisherSpectrum(
        matrix=sym,
        eigenvalues=values[order],
        eigenvectors=vectors[:, order],
        names=_default_names(sym.shape[0], names),
        threshold=float(threshold),
    )


@dataclass(frozen=True)
class CramerRaoBound:
    """A *lower bound* on the covariance of any unbiased estimator.  Not a promise of performance.

    ``Cov(theta_hat) >= F^-1``.  An efficient estimator attains it; a real one usually does not, and
    nothing here checks which case you are in.  That is what
    :func:`aleph.infer.posterior.coverage_campaign` is for.

    Attributes:
        covariance: ``F^-1``.
        variance: Its diagonal --- marginal, so it already includes the cost of not knowing the
            other parameters.
        names: Parameter names.
    """

    covariance: NDArray
    variance: NDArray
    names: tuple[str, ...]

    @property
    def standard_deviation(self) -> NDArray:
        return np.sqrt(self.variance)

    @property
    def correlation(self) -> NDArray:
        sd = self.standard_deviation
        return self.covariance / np.outer(sd, sd)

    def report(self) -> str:
        return "\n".join(
            f"  {n}: sigma >= {s:.6g}"
            for n, s in zip(self.names, self.standard_deviation, strict=True)
        )


def cramer_rao_bound(
    matrix: NDArray,
    *,
    names: Sequence[str] | None = None,
    threshold: float = DEFAULT_STIFFNESS_THRESHOLD,
) -> CramerRaoBound:
    """Invert the Fisher matrix, or refuse and name the direction that made it singular.

    A pseudo-inverse here would be a lie of the most convenient kind: it returns finite error bars
    for a direction the data cannot see at all, and nothing downstream can tell the difference.  So
    a Fisher matrix whose smallest eigenvalue falls below ``threshold * lambda_max`` raises
    :class:`~aleph.infer.refusal.UnidentifiableDirectionError` carrying the offending direction.
    """
    spectrum = fisher_spectrum(matrix, names=names, threshold=threshold)
    worst = spectrum.sloppiest
    cutoff = threshold * float(spectrum.eigenvalues[0])
    if float(spectrum.eigenvalues[0]) <= 0.0 or worst.eigenvalue <= cutoff:
        Unidentifiable(
            reason=(
                "The Fisher matrix is singular to working precision: the combination "
                f"[{worst.description}] has eigenvalue {worst.eigenvalue:.6g} against a largest "
                f"eigenvalue of {float(spectrum.eigenvalues[0]):.6g}. Inverting it would return a "
                "finite error bar for a direction this observation cannot see at all."
            ),
            offending_input={"fisher_eigenvalues": spectrum.eigenvalues.tolist()},
            lift_requires=(
                "Add an observation that responds to this combination, fix one of its parameters "
                "from independent evidence, or report the identifiable combinations instead of the "
                "individual parameters."
            ),
            detail={
                "condition_number": spectrum.condition_number,
                "threshold": float(threshold),
            },
            direction=tuple(float(x) for x in worst.vector),
            parameter_names=spectrum.names,
        ).raise_()

    cov = linalg.inv(spectrum.matrix)
    cov = 0.5 * (cov + cov.T)
    return CramerRaoBound(covariance=cov, variance=np.diag(cov).copy(), names=spectrum.names)


@dataclass(frozen=True)
class IdentifiabilitySplit:
    """What this observation can constrain, and --- equally important --- what it cannot.

    Attributes:
        identified: Stiff directions, stiffest first.
        unidentified: Sloppy directions, i.e. the ones the data barely or never sees.
        threshold: The relative eigenvalue cutoff used.
        spectrum: The full spectrum, so a caller can re-split at a different threshold.
    """

    identified: tuple[ParameterDirection, ...]
    unidentified: tuple[ParameterDirection, ...]
    threshold: float
    spectrum: FisherSpectrum

    @property
    def rank(self) -> int:
        return len(self.identified)

    @property
    def fully_identifiable(self) -> bool:
        return not self.unidentified

    def report(self) -> str:
        lines = [
            f"Identifiable directions at threshold {self.threshold:.3g} "
            f"({self.rank} of {len(self.spectrum.names)}):"
        ]
        lines.extend(f"  CAN measure   {d}" for d in self.identified)
        if self.unidentified:
            lines.extend(f"  CANNOT measure {d}" for d in self.unidentified)
        else:
            lines.append("  (every direction is constrained at this threshold)")
        return "\n".join(lines)


def identifiable_directions(
    matrix: NDArray,
    *,
    names: Sequence[str] | None = None,
    threshold: float = DEFAULT_STIFFNESS_THRESHOLD,
    absolute_threshold: float | None = None,
) -> IdentifiabilitySplit:
    """Split parameter space into what the observation constrains and what it does not.

    Args:
        matrix: Fisher matrix.
        names: Parameter names, used to render the directions.
        threshold: Relative cutoff --- an eigenvalue below ``threshold * lambda_max`` is
            unidentified.  Relative is the right default because it is invariant to rescaling the
            observation's units, which an absolute cutoff is not.
        absolute_threshold: When given, a direction must additionally have an eigenvalue above this
            to count as identified.  Use it to say "I need sigma better than X along every reported
            direction": pass ``1/X**2``.  A relative test alone will happily call a direction stiff
            when the whole measurement is uninformative.
    """
    spectrum = fisher_spectrum(matrix, names=names, threshold=threshold)
    cutoff = threshold * float(spectrum.eigenvalues[0])
    if absolute_threshold is not None:
        cutoff = max(cutoff, float(absolute_threshold))
    identified: list[ParameterDirection] = []
    unidentified: list[ParameterDirection] = []
    for i in range(spectrum.eigenvalues.size):
        base = spectrum.direction(i)
        stiff = bool(base.eigenvalue > cutoff)
        direction = ParameterDirection(
            index=base.index,
            eigenvalue=base.eigenvalue,
            vector=base.vector,
            names=base.names,
            stiff=stiff,
            threshold=float(threshold),
        )
        (identified if stiff else unidentified).append(direction)
    return IdentifiabilitySplit(
        identified=tuple(identified),
        unidentified=tuple(unidentified),
        threshold=float(threshold),
        spectrum=spectrum,
    )
