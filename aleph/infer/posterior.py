"""Posteriors that carry their own provenance, their own validity domain, and their own coverage.

A posterior in this project is not a mean and a covariance.  Four things travel with it, and each
one exists because omitting it is a specific way projects publish wrong numbers.

**Evidence lineage.**  Which artifact did this data come from?  A number whose data cannot be
retrieved is an assertion.  :class:`EvidenceLineage` is a required field, and a
:class:`Posterior` cannot be constructed without one --- synthetic data is allowed but must be
*labelled* synthetic, via :meth:`EvidenceLineage.for_synthetic`, so it can never be mistaken later
for a measurement.

**An out-of-distribution verdict.**  A posterior computed at parameters where the forward model was
never checked has an error bar that describes the noise and says nothing about the extrapolation.
:meth:`Posterior.credible_interval` refuses to hand out an interval for an out-of-distribution fit
unless the caller says, in the call, that they know.

**A representation error.**  The gap between the model's parameterization and the thing being
measured --- discretization, mode truncation, the finite band of ``q`` the observation covers --- is
not sampling noise and does not shrink with more data.  It is a separate field so that it cannot be
quietly folded into the statistical error bar, and :attr:`Posterior.total_standard_deviation` shows
what happens when it is added in quadrature.

**Measured coverage.**  This is the part that is usually skipped.  A 95% credible interval that
contains the truth 60% of the time is not a 95% interval, and nothing about deriving it correctly
prevents that --- a Laplace approximation to a skewed posterior, or a Cramer-Rao bound treated as an
error bar, will both do it.  :func:`coverage_campaign` draws from truth, runs the whole inference,
and measures the frequentist containment rate with a binomial confidence interval on the measurement
itself.  A posterior that has never been through this is decoration.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy import linalg, optimize, stats

from .numdiff import hessian as numeric_hessian
from .refusal import (
    ModelInadequate,
    OutOfDistribution,
    ProvenanceIncomplete,
    Refusal,
    RefusalError,
    Unidentifiable,
)

__all__ = [
    "CoverageReport",
    "CoverageResult",
    "EvidenceLineage",
    "GaussianApproximation",
    "OutOfDistributionVerdict",
    "Posterior",
    "RepresentationError",
    "coverage_campaign",
    "effective_sample_size",
    "gaussian_discrepancy",
    "laplace_approximation",
    "metropolis_sample",
    "wilson_interval",
]

_DIGEST_PATTERN = re.compile(r"^[0-9a-f]{16,128}$")


@dataclass(frozen=True)
class EvidenceLineage:
    """Where the data came from, in a form that can be looked up again.

    Attributes:
        artifact_digest: Content digest of the artifact the data was read from, lowercase hex.
        observation: Name of the observation operator or protocol that produced it.
        synthetic: True when the "data" was generated for a test or a coverage campaign.  Synthetic
            evidence is legitimate and necessary; passing it off as a measurement is not, so the
            flag is structural rather than a naming convention.
        note: Free text for anything else worth carrying.
    """

    artifact_digest: str
    observation: str
    synthetic: bool = False
    note: str = ""

    def __post_init__(self) -> None:
        missing = [
            name
            for name in ("artifact_digest", "observation")
            if not str(getattr(self, name)).strip()
        ]
        if missing:
            ProvenanceIncomplete(
                reason="Evidence lineage is incomplete; no number derived from it may be quoted.",
                offending_input={
                    "artifact_digest": self.artifact_digest,
                    "observation": self.observation,
                },
                lift_requires=(
                    "Record the artifact digest and the observation operator behind the data."
                ),
                missing_fields=tuple(missing),
            ).raise_()
        if not _DIGEST_PATTERN.match(self.artifact_digest):
            ProvenanceIncomplete(
                reason=(
                    f"artifact_digest {self.artifact_digest!r} is not a content digest. A label "
                    "that looks like a digest but is not one is worse than no digest, because it "
                    "survives review."
                ),
                offending_input=self.artifact_digest,
                lift_requires=(
                    "Supply at least 16 lowercase hex characters of a real content digest."
                ),
                missing_fields=("artifact_digest",),
            ).raise_()

    @classmethod
    def for_synthetic(
        cls, description: str, *, observation: str = "synthetic"
    ) -> "EvidenceLineage":
        """Lineage for data this process generated itself, digested over its own description.

        The digest is reproducible from the description, so a synthetic run is as traceable as a
        real one --- and the ``synthetic`` flag makes it impossible to quote as a measurement by
        accident.
        """
        if not description.strip():
            raise ValueError("synthetic evidence still needs a description of how it was generated")
        digest = hashlib.sha256(("aleph-synthetic:" + description).encode("utf-8")).hexdigest()
        return cls(
            artifact_digest=digest, observation=observation, synthetic=True, note=description
        )


@dataclass(frozen=True)
class OutOfDistributionVerdict:
    """Whether this inference sits inside the region where the forward model was ever checked.

    Attributes:
        in_distribution: The verdict.
        statistic: Whatever the check computed --- a Mahalanobis distance, a box violation, a
            chi-squared.  Recorded so the verdict can be re-derived at another threshold.
        threshold: The value it was compared against.
        basis: What was actually checked, in words.
    """

    in_distribution: bool
    statistic: float
    threshold: float
    basis: str

    @classmethod
    def in_box(
        cls,
        point: NDArray,
        lower: NDArray,
        upper: NDArray,
        *,
        basis: str = "parameter box where the forward model was validated",
    ) -> "OutOfDistributionVerdict":
        """Verdict from a validated parameter box.

        The statistic is the largest violation in units of the box's own half-width, so ``0`` means
        dead centre, ``1`` means exactly on a face, and ``>1`` means outside.  Expressing it in
        box-relative units keeps the number comparable across parameters with different scales.
        """
        p = np.asarray(point, dtype=float)
        lo = np.asarray(lower, dtype=float)
        hi = np.asarray(upper, dtype=float)
        if not (p.shape == lo.shape == hi.shape):
            raise ValueError("point, lower and upper must have the same shape")
        if np.any(hi <= lo):
            raise ValueError("every upper bound must exceed its lower bound")
        centre = 0.5 * (lo + hi)
        half = 0.5 * (hi - lo)
        stat = float(np.max(np.abs(p - centre) / half))
        return cls(in_distribution=stat <= 1.0, statistic=stat, threshold=1.0, basis=basis)

    def refusal(self, offending_input: Any = None) -> OutOfDistribution | None:
        if self.in_distribution:
            return None
        return OutOfDistribution(
            reason=(
                f"This inference sits outside {self.basis}: statistic {self.statistic:.4g} against "
                f"threshold {self.threshold:.4g}. Its error bar describes the noise, not the "
                "extrapolation, and the two are unrelated."
            ),
            offending_input=offending_input,
            lift_requires=(
                "Validate the forward model over the region this fit landed in, or restrict the "
                "support and re-fit."
            ),
            statistic=self.statistic,
            threshold=self.threshold,
        )


@dataclass(frozen=True)
class RepresentationError:
    """Error from the parameterization itself, which more data does not reduce.

    Mode truncation, the finite ``q`` band an observation covers, a discretization that is not the
    continuum --- all of these bias the estimate by an amount set by the representation, not by the
    sample size.  Kept as its own field so it cannot be folded into the statistical error bar, where
    it would appear to shrink as ``1/sqrt(N)`` when it does not shrink at all.

    Attributes:
        value: Magnitude in parameter units, per parameter, or ``None`` when it has not been
            quantified.
        basis: What the number is, or --- when unquantified --- why not.
        quantified: Whether ``value`` may be used.
    """

    value: NDArray | None
    basis: str
    quantified: bool

    @classmethod
    def quantified_as(cls, value: NDArray | float, basis: str) -> "RepresentationError":
        arr = np.atleast_1d(np.asarray(value, dtype=float))
        if np.any(arr < 0.0):
            raise ValueError("representation error must be non-negative")
        if not basis.strip():
            raise ValueError("a quantified representation error must say what it is")
        return cls(value=arr, basis=basis, quantified=True)

    @classmethod
    def unquantified(cls, reason: str) -> "RepresentationError":
        if not reason.strip():
            raise ValueError("say why the representation error has not been quantified")
        return cls(value=None, basis=reason, quantified=False)


@dataclass(frozen=True)
class GaussianApproximation:
    """Output of the Laplace step: a mode and a curvature, before any provenance is attached."""

    mean: NDArray
    covariance: NDArray
    log_posterior_max: float
    converged: bool
    n_iterations: int
    message: str

    @property
    def standard_deviation(self) -> NDArray:
        return np.sqrt(np.diag(self.covariance))


@dataclass(frozen=True)
class Posterior:
    """A posterior with everything needed to decide whether its numbers may be quoted.

    Either ``covariance`` (a Gaussian approximation) or ``samples`` must be present; both is fine
    and is how a Laplace fit checked by MCMC is recorded.

    Attributes:
        parameter_names: Names, in the order components appear.
        mean: Posterior mean, or the MAP for a Laplace fit.
        covariance: Gaussian approximation to the covariance, if any.
        samples: ``(n_samples, n_parameters)`` draws, if any.
        support_lower: Lower edge of the region this posterior is claimed to be valid in.
        support_upper: Upper edge of the same.
        lineage: Where the data came from.  Required.
        ood: Out-of-distribution verdict.  Required.
        representation_error: Non-statistical error.  Required, possibly unquantified.
        method: How it was produced, e.g. ``"laplace"`` or ``"metropolis"``.
        diagnostics: Whatever the method wants to record --- acceptance rate, ESS, iterations.
    """

    parameter_names: tuple[str, ...]
    mean: NDArray
    support_lower: NDArray
    support_upper: NDArray
    lineage: EvidenceLineage
    ood: OutOfDistributionVerdict
    representation_error: RepresentationError
    covariance: NDArray | None = None
    samples: NDArray | None = None
    method: str = "unspecified"
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        n = len(self.parameter_names)
        mean = np.asarray(self.mean, dtype=float).ravel()
        if mean.size != n:
            raise ValueError(f"mean has {mean.size} components but {n} parameter names were given")
        object.__setattr__(self, "mean", mean)
        for name in ("support_lower", "support_upper"):
            arr = np.asarray(getattr(self, name), dtype=float).ravel()
            if arr.size != n:
                raise ValueError(f"{name} must have {n} components, got {arr.size}")
            object.__setattr__(self, name, arr)
        if np.any(self.support_upper <= self.support_lower):
            raise ValueError("every support upper bound must exceed its lower bound")
        if self.covariance is None and self.samples is None:
            raise ValueError(
                "a posterior must carry either a covariance or samples; an object with neither "
                "records only a point estimate and must not be called a posterior"
            )
        if self.covariance is not None:
            cov = np.asarray(self.covariance, dtype=float)
            if cov.shape != (n, n):
                raise ValueError(f"covariance must have shape ({n}, {n}), got {cov.shape}")
            object.__setattr__(self, "covariance", 0.5 * (cov + cov.T))
        if self.samples is not None:
            s = np.atleast_2d(np.asarray(self.samples, dtype=float))
            if s.shape[1] != n:
                raise ValueError(f"samples must have {n} columns, got {s.shape[1]}")
            object.__setattr__(self, "samples", s)
        if not isinstance(self.lineage, EvidenceLineage):
            ProvenanceIncomplete(
                reason="A posterior must carry the evidence lineage of the data it was fitted to.",
                offending_input=self.lineage,
                lift_requires=(
                    "Construct an EvidenceLineage naming the artifact digest and observation."
                ),
                missing_fields=("lineage",),
            ).raise_()
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))
        object.__setattr__(self, "parameter_names", tuple(str(x) for x in self.parameter_names))

    @property
    def n_parameters(self) -> int:
        return len(self.parameter_names)

    @property
    def marginal_standard_deviation(self) -> NDArray:
        """Statistical error only.  Excludes representation error, deliberately."""
        if self.samples is not None:
            return np.std(self.samples, axis=0, ddof=1)
        return np.sqrt(np.diag(np.asarray(self.covariance)))

    @property
    def total_standard_deviation(self) -> NDArray:
        """Statistical and representation error in quadrature, when the latter is quantified.

        When it is not quantified this returns the statistical error unchanged --- and that is
        exactly the situation where the reported error bar is an underestimate of unknown size.
        :attr:`RepresentationError.basis` says why it is unquantified; read it before quoting.
        """
        sd = self.marginal_standard_deviation
        rep = self.representation_error
        if not rep.quantified or rep.value is None:
            return sd
        extra = np.broadcast_to(np.asarray(rep.value, dtype=float), sd.shape)
        return np.sqrt(sd**2 + extra**2)

    def credible_interval(
        self, level: float = 0.95, *, allow_out_of_distribution: bool = False
    ) -> tuple[NDArray, NDArray]:
        """Marginal equal-tailed credible interval per parameter.

        Refuses when the fit is out of distribution unless the caller says otherwise in the call.
        ``allow_out_of_distribution=True`` is the correct setting inside a coverage campaign, where
        the point is to measure the machinery, and the wrong setting anywhere a number is reported.

        From samples this is the empirical quantile pair; from a Gaussian it is ``mean +/- z*sigma``
        with ``sigma`` the *total* standard deviation.  **Either way a quantified representation
        error widens the interval** rather than sitting in a field nobody reads.

        On the samples branch the widening is applied in quadrature to each endpoint's distance
        from the median, which keeps whatever asymmetry the samples carry and reduces exactly to
        the Gaussian branch when they are Gaussian.  The two branches agreeing matters more than it
        looks: a posterior may hold both a covariance and samples, and this module recommends MCMC
        as the check on a Laplace fit, so a branch that dropped the term meant running the
        recommended check *narrowed* the published interval.  The careful analyst got the smaller
        error bar, and ``summary()`` still said "representation error: quantified".
        """
        if not 0.0 < level < 1.0:
            raise ValueError(f"level must lie in (0, 1), got {level}")
        if not allow_out_of_distribution:
            refusal = self.ood.refusal(offending_input={"mean": self.mean.tolist()})
            if refusal is not None:
                refusal.raise_()
        tail = 0.5 * (1.0 - level)
        z = float(stats.norm.ppf(1.0 - tail))
        if self.samples is not None:
            lo = np.quantile(self.samples, tail, axis=0)
            hi = np.quantile(self.samples, 1.0 - tail, axis=0)
            rep = self.representation_error
            if not rep.quantified or rep.value is None:
                return lo, hi
            median = np.quantile(self.samples, 0.5, axis=0)
            extra = z * np.broadcast_to(np.asarray(rep.value, dtype=float), median.shape)
            below = np.sqrt((median - lo) ** 2 + extra**2)
            above = np.sqrt((hi - median) ** 2 + extra**2)
            return median - below, median + above
        sd = self.total_standard_deviation
        return self.mean - z * sd, self.mean + z * sd

    def contains(
        self, truth: NDArray, level: float = 0.95, *, allow_out_of_distribution: bool = True
    ) -> NDArray:
        """Boolean per parameter: is ``truth`` inside the marginal credible interval?"""
        lo, hi = self.credible_interval(level, allow_out_of_distribution=allow_out_of_distribution)
        t = np.asarray(truth, dtype=float).ravel()
        return (t >= lo) & (t <= hi)

    def contains_jointly(self, truth: NDArray, level: float = 0.95) -> bool:
        """Is ``truth`` inside the *joint* credible region?

        For a Gaussian posterior the joint region is the ellipsoid
        ``(t - mu)^T Sigma^-1 (t - mu) <= chi2_k(level)``.  This is the honest multi-parameter
        statement and it is not implied by every marginal interval containing the truth --- two
        strongly correlated parameters can each be individually covered while the pair sits far
        outside the joint region.
        """
        t = np.asarray(truth, dtype=float).ravel()
        delta = t - self.mean
        if self.covariance is not None:
            cov = np.asarray(self.covariance, dtype=float)
        else:
            cov = np.cov(np.asarray(self.samples), rowvar=False)
            cov = np.atleast_2d(cov)
        try:
            chol = linalg.cho_factor(cov, lower=True)
        except linalg.LinAlgError:
            Unidentifiable(
                reason="The posterior covariance is singular, so no joint credible region exists.",
                offending_input={"covariance": np.asarray(cov).tolist()},
                lift_requires=(
                    "Constrain or drop the degenerate direction before asking for a region."
                ),
                parameter_names=self.parameter_names,
            ).raise_()
        d2 = float(delta @ linalg.cho_solve(chol, delta))
        return d2 <= float(stats.chi2.ppf(level, self.n_parameters))

    def interval_within_support(self, level: float = 0.95) -> bool:
        """Does the credible interval stay inside the region this posterior is valid in?

        Intervals are never truncated to the support here.  Truncation would hide the situation
        this method exists to report: an interval spilling past the validated region means the
        posterior is being asked about parameters nobody checked.
        """
        lo, hi = self.credible_interval(level, allow_out_of_distribution=True)
        return bool(np.all(lo >= self.support_lower) and np.all(hi <= self.support_upper))

    def summary(self, level: float = 0.95) -> str:
        lo, hi = self.credible_interval(level, allow_out_of_distribution=True)
        head = (
            f"Posterior by {self.method} over {list(self.parameter_names)}\n"
            f"  lineage: {self.lineage.observation} @ {self.lineage.artifact_digest[:12]}"
            f"{' [SYNTHETIC]' if self.lineage.synthetic else ''}\n"
            f"  out-of-distribution: {'NO' if self.ood.in_distribution else 'YES'} "
            f"(statistic {self.ood.statistic:.4g} vs {self.ood.threshold:.4g})\n"
            f"  representation error: "
            f"{'quantified' if self.representation_error.quantified else 'UNQUANTIFIED'} "
            f"-- {self.representation_error.basis}"
        )
        rows = [
            f"  {n}: {m:.6g} [{a:.6g}, {b:.6g}] at {100 * level:g}%"
            for n, m, a, b in zip(self.parameter_names, self.mean, lo, hi, strict=True)
        ]
        return "\n".join([head, *rows])


def _value_and_gradient(target: Any) -> Callable[[NDArray], tuple[float, NDArray]]:
    """Accept either a callable returning ``(value, gradient)`` or an object exposing it."""
    if hasattr(target, "value_and_gradient"):
        return target.value_and_gradient
    return target


def laplace_approximation(
    log_posterior: Any,
    theta_init: NDArray,
    *,
    bounds: Sequence[tuple[float, float]] | None = None,
    tol: float = 1e-10,
    max_iterations: int = 500,
) -> GaussianApproximation:
    """Maximize the log posterior, then take its curvature at the mode.

    The Laplace approximation is the second-order Taylor expansion of ``log p`` at the MAP:
    ``Sigma = (-H)^-1`` with ``H`` the Hessian of the log posterior.  It is exact for a Gaussian
    posterior --- which is the positive control --- and it is an approximation everywhere else,
    with an error that grows with the skew of the true posterior.  That is what
    :func:`metropolis_sample` and :func:`gaussian_discrepancy` are here to measure, and what
    :func:`coverage_campaign` is here to price.

    The Hessian is obtained by central-differencing the *analytic* gradient rather than
    double-differencing the value: one numerical derivative instead of two, so the curvature
    inherits the gradient's accuracy rather than squaring its error.

    A Hessian that is not negative definite at the mode is not a numerical annoyance; it means
    there is no local maximum there --- either the optimizer stopped early or a direction is flat.
    Both get a typed refusal instead of a pseudo-inverse.
    """
    fg = _value_and_gradient(log_posterior)
    theta0 = np.asarray(theta_init, dtype=float).ravel()

    def negative(theta: NDArray) -> tuple[float, NDArray]:
        value, grad = fg(np.asarray(theta, dtype=float))
        return -float(value), -np.asarray(grad, dtype=float)

    method = "L-BFGS-B" if bounds is not None else "BFGS"
    result = optimize.minimize(
        negative,
        theta0,
        jac=True,
        method=method,
        bounds=bounds,
        options={"maxiter": max_iterations, "gtol": tol}
        if method == "BFGS"
        else {"maxiter": max_iterations, "ftol": tol, "gtol": tol},
    )
    mode = np.asarray(result.x, dtype=float)

    def grad_only(theta: NDArray) -> NDArray:
        return np.asarray(fg(theta)[1], dtype=float)

    hess = numeric_hessian(grad_only, mode)
    negative_hessian = -hess
    eigenvalues = np.linalg.eigvalsh(negative_hessian)
    if np.min(eigenvalues) <= 0.0:
        ModelInadequate(
            reason=(
                "The Hessian of the log posterior is not negative definite at the point the "
                f"optimizer stopped: eigenvalues of -H are {eigenvalues.tolist()}. There is no "
                "local maximum here, so there is no Laplace approximation to take."
            ),
            offending_input={"theta": mode.tolist(), "optimizer_message": str(result.message)},
            lift_requires=(
                "Start the optimizer elsewhere, add a proper prior to curve the flat direction, or "
                "drop the direction reported as unidentifiable."
            ),
            detail={"eigenvalues": eigenvalues.tolist(), "converged": bool(result.success)},
        ).raise_()

    covariance = linalg.inv(negative_hessian)
    covariance = 0.5 * (covariance + covariance.T)
    return GaussianApproximation(
        mean=mode,
        covariance=covariance,
        log_posterior_max=float(-result.fun),
        converged=bool(result.success),
        n_iterations=int(getattr(result, "nit", -1)),
        message=str(result.message),
    )


def metropolis_sample(
    log_posterior: Any,
    theta_init: NDArray,
    *,
    proposal_covariance: NDArray,
    n_samples: int,
    rng: np.random.Generator,
    burn_in: int = 1000,
    thin: int = 1,
    scale: float | None = None,
) -> tuple[NDArray, dict[str, Any]]:
    """Random-walk Metropolis, as an independent check on the Laplace approximation.

    Deliberately the simplest correct sampler: a Gaussian proposal, a Metropolis acceptance ratio,
    no adaptation.  Its job is not to be fast, it is to be obviously right, so that when it
    disagrees with the Laplace approximation the disagreement is evidence about the posterior rather
    than about the sampler.

    The proposal is scaled by ``2.38/sqrt(d)`` against the supplied covariance --- the standard
    optimal scaling for a Gaussian target, which puts the acceptance rate near 0.23-0.44 for the
    dimensions used here.  Pass ``proposal_covariance`` from a Laplace fit; a proposal matched to
    the target's own shape is what makes this usable without tuning.

    Returns the post-burn-in, thinned samples and a diagnostics dict carrying the acceptance rate
    and the per-parameter effective sample size.
    """
    fg = _value_and_gradient(log_posterior)
    theta = np.asarray(theta_init, dtype=float).ravel()
    dim = theta.size
    cov = np.asarray(proposal_covariance, dtype=float)
    if cov.shape != (dim, dim):
        raise ValueError(f"proposal_covariance must have shape ({dim}, {dim}), got {cov.shape}")
    if n_samples < 1 or thin < 1 or burn_in < 0:
        raise ValueError("need n_samples >= 1, thin >= 1, burn_in >= 0")
    factor = (2.38 / np.sqrt(dim)) if scale is None else float(scale)
    chol = np.linalg.cholesky(cov * factor**2)

    current = float(fg(theta)[0])
    total = burn_in + n_samples * thin
    kept = np.empty((n_samples, dim), dtype=float)
    accepted = 0
    index = 0
    for step in range(total):
        proposal = theta + chol @ rng.standard_normal(dim)
        try:
            proposed = float(fg(proposal)[0])
        except RefusalError:
            # A proposal the model refuses to evaluate is outside the support. Rejecting it is the
            # correct Metropolis move for a zero-density point, not an error to propagate.
            proposed = -np.inf
        if np.log(rng.random()) < proposed - current:
            theta, current = proposal, proposed
            accepted += 1
        if step >= burn_in and (step - burn_in) % thin == 0 and index < n_samples:
            kept[index] = theta
            index += 1
    ess = np.array([effective_sample_size(kept[:, j]) for j in range(dim)])
    return kept, {
        "acceptance_rate": accepted / total,
        "effective_sample_size": ess,
        "n_steps": total,
        "proposal_scale": factor,
    }


def effective_sample_size(chain: NDArray) -> float:
    """Effective sample size of a 1-D chain, by the initial-positive-sequence rule.

    ``ESS = n / (1 + 2*sum_k rho_k)``, with the autocorrelation sum truncated at the first
    non-positive term.  The truncation matters: summing the full autocorrelation function adds pure
    noise at large lags and routinely produces an ESS larger than the chain, which then gets quoted.
    """
    x = np.asarray(chain, dtype=float).ravel()
    n = x.size
    if n < 2:
        return float(n)
    centred = x - x.mean()
    var = float(np.dot(centred, centred))
    if var <= 0.0:
        return 1.0  # a constant chain carries exactly one independent piece of information
    size = int(2 ** np.ceil(np.log2(2 * n)))
    spectrum = np.fft.rfft(centred, n=size)
    acov = np.fft.irfft(spectrum * np.conjugate(spectrum), n=size)[:n]
    rho = acov / acov[0]
    total = 0.0
    for k in range(1, n):
        if rho[k] <= 0.0:
            break
        total += rho[k]
    return float(n / (1.0 + 2.0 * total))


def gaussian_discrepancy(
    samples: NDArray, mean: NDArray, covariance: NDArray
) -> dict[str, float]:
    """How far a sample cloud departs from a claimed Gaussian --- the non-Gaussianity check.

    Returns, per the worst parameter:

    - ``max_mean_shift_sigma``: distance between sample mean and claimed mean, in units of the
      claimed standard deviation.  Compare it against the Monte Carlo error ``1/sqrt(ESS)``, not
      against zero.
    - ``max_abs_log_variance_ratio``: ``|log(sample variance / claimed variance)|``.
    - ``max_abs_skewness``: sample skewness, which is 0 for any Gaussian and is the cheapest
      indicator that a Laplace approximation is centred on the wrong side of the mass.
    """
    s = np.atleast_2d(np.asarray(samples, dtype=float))
    mu = np.asarray(mean, dtype=float).ravel()
    sd = np.sqrt(np.diag(np.asarray(covariance, dtype=float)))
    sample_mean = s.mean(axis=0)
    sample_var = s.var(axis=0, ddof=1)
    return {
        "max_mean_shift_sigma": float(np.max(np.abs(sample_mean - mu) / sd)),
        "max_abs_log_variance_ratio": float(np.max(np.abs(np.log(sample_var / sd**2)))),
        "max_abs_skewness": float(np.max(np.abs(stats.skew(s, axis=0)))),
    }


def wilson_interval(successes: int, trials: int, confidence: float = 0.95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion.

    The normal approximation ``p +/- z*sqrt(p(1-p)/n)`` is wrong exactly where a coverage campaign
    lives --- near ``p = 0.95``, where it can return an upper limit above 1 and understates the
    interval for modest ``n``.  The Wilson interval is the inversion of the score test and stays
    inside ``[0, 1]`` by construction.
    """
    if trials <= 0:
        raise ValueError("need at least one trial")
    if not 0 <= successes <= trials:
        raise ValueError(f"successes must lie in [0, {trials}], got {successes}")
    z = float(stats.norm.ppf(0.5 * (1.0 + confidence)))
    p = successes / trials
    denom = 1.0 + z**2 / trials
    centre = (p + z**2 / (2 * trials)) / denom
    half = (z / denom) * np.sqrt(p * (1 - p) / trials + z**2 / (4 * trials**2))
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


@dataclass(frozen=True)
class CoverageResult:
    """Measured containment rate for one parameter (or for the joint region).

    Attributes:
        name: Parameter name, or ``"joint"``.
        nominal: The credible level that was claimed.
        n_trials: Trials that produced a posterior.
        n_covered: Trials whose interval contained the truth.
        measured: ``n_covered / n_trials``.
        ci_low, ci_high: Binomial confidence interval on ``measured`` itself.
        ci_confidence: Confidence level of that interval.
    """

    name: str
    nominal: float
    n_trials: int
    n_covered: int
    measured: float
    ci_low: float
    ci_high: float
    ci_confidence: float

    @property
    def consistent(self) -> bool:
        """Is the nominal level inside the measured interval?

        This is the whole test.  ``False`` means the intervals are demonstrably the wrong size ---
        under-covering if ``measured < nominal``, which is the dangerous direction, or
        over-covering if above, which wastes information but does not mislead.
        """
        return self.ci_low <= self.nominal <= self.ci_high

    def __str__(self) -> str:
        verdict = "OK" if self.consistent else "MISMATCH"
        return (
            f"{self.name}: nominal {self.nominal:.3f}, measured {self.measured:.4f} "
            f"[{self.ci_low:.4f}, {self.ci_high:.4f}] at {100 * self.ci_confidence:g}% "
            f"({self.n_covered}/{self.n_trials}) -> {verdict}"
        )


@dataclass(frozen=True)
class CoverageReport:
    """Result of a coverage campaign, including the trials that refused."""

    level: float
    marginal: tuple[CoverageResult, ...]
    joint: CoverageResult | None
    n_requested: int
    n_refused: int
    refusals: tuple[Refusal, ...]

    @property
    def all_consistent(self) -> bool:
        results = list(self.marginal) + ([self.joint] if self.joint is not None else [])
        return all(r.consistent for r in results)

    def render(self) -> str:
        lines = [
            f"Coverage campaign at nominal {self.level:.3f}: "
            f"{self.n_requested - self.n_refused}/{self.n_requested} trials produced a posterior, "
            f"{self.n_refused} refused."
        ]
        lines.extend(f"  {r}" for r in self.marginal)
        if self.joint is not None:
            lines.append(f"  {self.joint}")
        return "\n".join(lines)


def coverage_campaign(
    truth_sampler: Callable[[np.random.Generator], NDArray],
    simulate: Callable[[NDArray, np.random.Generator], NDArray],
    fit: Callable[[NDArray], Posterior],
    *,
    n_trials: int,
    rng: np.random.Generator,
    level: float = 0.95,
    ci_confidence: float = 0.95,
    include_joint: bool = True,
) -> CoverageReport:
    """Draw from truth, run the whole inference, and measure how often the interval contains it.

    This is a frequentist audit of a Bayesian object, and that is the point.  A correctly derived
    posterior under a correctly specified model, averaged over draws from the prior it assumes, has
    calibrated coverage --- so a measured rate that misses nominal is evidence of one of exactly
    three things: the approximation (Laplace on a skewed posterior), the model (misspecified
    forward model or noise), or a bug.  All three are worth knowing about, and none of them is
    visible from the posterior itself.

    Trials whose fit *refuses* are counted separately rather than dropped.  Dropping them would
    quietly condition the coverage on the refusal, which biases the result toward the easy cases ---
    the exact failure mode of reporting only the runs that converged.

    Args:
        truth_sampler: ``rng -> true parameter vector``.  Draw from the prior the posterior assumes;
            drawing from anything else measures a different quantity.
        simulate: ``(truth, rng) -> data``, including observation noise.
        fit: ``data -> Posterior``.  May raise a
            :class:`~aleph.infer.refusal.RefusalError`, which is recorded.
        n_trials: Number of independent trials.
        rng: Generator, for reproducibility.
        level: Nominal credible level being tested.
        ci_confidence: Confidence level of the binomial interval on the measured coverage.
        include_joint: Also measure coverage of the joint credible region.
    """
    if n_trials < 1:
        raise ValueError("need at least one trial")

    names: tuple[str, ...] | None = None
    marginal_hits: NDArray | None = None
    joint_hits = 0
    n_done = 0
    refusals: list[Refusal] = []

    for _ in range(n_trials):
        truth = np.asarray(truth_sampler(rng), dtype=float).ravel()
        data = simulate(truth, rng)
        try:
            posterior = fit(data)
        except RefusalError as exc:
            refusals.append(exc.refusal)
            continue
        if names is None:
            names = posterior.parameter_names
            marginal_hits = np.zeros(len(names), dtype=int)
        elif posterior.parameter_names != names:
            raise ValueError("every trial must fit the same parameters in the same order")
        assert marginal_hits is not None
        covered = posterior.contains(truth, level, allow_out_of_distribution=True)
        marginal_hits += covered.astype(int)
        if include_joint:
            joint_hits += int(posterior.contains_jointly(truth, level))
        n_done += 1

    if n_done == 0 or names is None or marginal_hits is None:
        raise RuntimeError(
            f"every one of {n_trials} trials refused; there is nothing to measure coverage on. "
            f"First refusal: {refusals[0].message() if refusals else 'none recorded'}"
        )

    results = tuple(
        CoverageResult(
            name=name,
            nominal=level,
            n_trials=n_done,
            n_covered=int(hits),
            measured=float(hits) / n_done,
            ci_low=wilson_interval(int(hits), n_done, ci_confidence)[0],
            ci_high=wilson_interval(int(hits), n_done, ci_confidence)[1],
            ci_confidence=ci_confidence,
        )
        for name, hits in zip(names, marginal_hits, strict=True)
    )
    joint = None
    if include_joint:
        lo, hi = wilson_interval(joint_hits, n_done, ci_confidence)
        joint = CoverageResult(
            name="joint",
            nominal=level,
            n_trials=n_done,
            n_covered=joint_hits,
            measured=joint_hits / n_done,
            ci_low=lo,
            ci_high=hi,
            ci_confidence=ci_confidence,
        )
    return CoverageReport(
        level=level,
        marginal=results,
        joint=joint,
        n_requested=n_trials,
        n_refused=len(refusals),
        refusals=tuple(refusals),
    )
