"""Log-likelihoods for the observation noise models this project actually meets, with gradients.

Everything here is written against a generic forward model::

    ForwardModel = Callable[[NDArray], NDArray]   # parameters -> predicted observable

Nothing in this module knows what a membrane is.  That is deliberate: the observation operator and
the analytic oracles are being built in parallel lanes, and a likelihood that imports either of them
becomes untestable on its own and starts validating itself against the thing it is judging.  Every
test in ``tests/infer`` therefore supplies its own closed-form toy model.

**The choice that matters: a-priori variance, never a fitted one.**

For a spectrum whose bins carry ``nu`` chi-squared degrees of freedom, the noise on the log amplitude
is known before any data is collected --- it is ``2/nu`` to leading order, set by the degrees of
freedom and nothing else.  Using that number means a wrong model has nowhere to hide: its residuals
are compared against
an error bar it cannot influence, chi-squared per degree of freedom rises, and
:class:`~aleph.infer.refusal.ModelInadequate` fires.  Fit the variance instead --- profile it out,
or call it a nuisance parameter with a flat prior --- and the maximum-likelihood variance is exactly
the mean squared residual, so reduced chi-squared is 1 *by construction* for every model, correct or
not.  The goodness-of-fit test silently stops existing.  :func:`gaussian_profiled_variance`
implements that anti-pattern explicitly, so the negative control can demonstrate it rather than
describe it.

**Why a Gamma likelihood as well.**  A periodogram bin is Gamma-distributed, not Gaussian:
``P_hat ~ Gamma(shape=nu/2, scale=2P/nu)``, equivalently ``nu*P_hat/P ~ chi2_nu``.  At large ``nu``
the log-amplitude Gaussian is an excellent approximation, but the first target's spectrum has few
averages in its highest-``q`` bins, where the Gamma tail is heavy and the Gaussian approximation
biases the fit low.  This is the one-dimensional case of a Wishart likelihood for spectral density
matrices; the multivariate extension slots in behind the same interface.

**Both spectral factories take ``chi_squared_dof``, and that is a repair.**  They previously took an
argument each spelled ``n_average``, and read it as different degrees of freedom --- ``nu = N`` in
the log-Gaussian, ``nu = 2N`` in the Gamma.  Each is defensible alone (a real field contributes one
dof per snapshot, a complex mode two), and together they made the sentence above false: the two
disagreed on the information in a bin by a factor of two *at every* ``N``, so the approximation
never converged.  Every test graded a convention against itself and nothing failed.  A count of
snapshots does not determine ``nu`` without also stating whether the field is real or complex, so
the count is no longer accepted --- the caller states ``nu``.

Every gradient here is analytic, and every one is checked against central differences in
``tests/infer/test_likelihood.py``.  The finite-difference check is the arbiter --- if the analytic
gradient and the numerical one disagree, the analytic one is wrong, including its sign.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping

import numpy as np
from numpy.typing import NDArray
from scipy import special, stats

from .numdiff import jacobian as numeric_jacobian
from .refusal import ModelInadequate

__all__ = [
    "ChiSquaredReport",
    "ForwardModel",
    "JacobianFn",
    "LogLikelihood",
    "chi_squared_report",
    "gamma_power_spectrum",
    "gaussian_apriori_log_spectrum",
    "gaussian_known_variance",
    "gaussian_profiled_variance",
    "log_amplitude_variance",
]

#: Parameters -> predicted observable.  The only contract this layer has with the physics.
ForwardModel = Callable[[NDArray], NDArray]
#: Parameters -> d(observable_i)/d(parameter_j), shape ``(n_observations, n_parameters)``.
JacobianFn = Callable[[NDArray], NDArray]

_LOG_2PI = float(np.log(2.0 * np.pi))


@dataclass(frozen=True)
class LogLikelihood:
    """A log-likelihood together with its analytic gradient and enough context to audit it.

    Attributes:
        name: Which noise model this is.  Appears in refusals and posterior records.
        n_data: Number of observations entering the sum.
        n_parameters: Dimension of the parameter vector, when the factory could determine it.
        evaluate: ``theta -> (log_likelihood, gradient)``.
        noise_variance_is_fitted: True when the noise scale was estimated from the residuals.  A
            likelihood with this set cannot support a goodness-of-fit test, and
            :func:`chi_squared_report` refuses to pretend otherwise.
        metadata: Free-form context recorded alongside a posterior.
    """

    name: str
    n_data: int
    n_parameters: int | None
    evaluate: Callable[[NDArray], tuple[float, NDArray]]
    noise_variance_is_fitted: bool = False
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", dict(self.metadata))

    def __call__(self, theta: NDArray) -> float:
        return self.value(theta)

    def value(self, theta: NDArray) -> float:
        return self.evaluate(np.asarray(theta, dtype=float))[0]

    def gradient(self, theta: NDArray) -> NDArray:
        return self.evaluate(np.asarray(theta, dtype=float))[1]

    def value_and_gradient(self, theta: NDArray) -> tuple[float, NDArray]:
        return self.evaluate(np.asarray(theta, dtype=float))


def _as_1d(x: Any, what: str) -> NDArray:
    arr = np.atleast_1d(np.asarray(x, dtype=float))
    if arr.ndim != 1:
        raise ValueError(f"{what} must be 1-D, got shape {arr.shape}")
    return arr


def _broadcast_variance(variance: Any, n: int) -> NDArray:
    v = np.atleast_1d(np.asarray(variance, dtype=float))
    if v.size == 1:
        v = np.full(n, float(v[0]))
    if v.shape != (n,):
        raise ValueError(f"variance must be scalar or length {n}, got shape {v.shape}")
    if not np.all(v > 0.0):
        raise ValueError("every variance entry must be strictly positive")
    return v


def _jacobian_for(model: ForwardModel, analytic: JacobianFn | None) -> JacobianFn:
    if analytic is not None:
        return lambda theta: np.asarray(analytic(theta), dtype=float)
    return lambda theta: numeric_jacobian(model, theta)


def _predict(model: ForwardModel, theta: NDArray, n_data: int, name: str) -> NDArray:
    pred = _as_1d(model(theta), f"{name}: forward model output")
    if pred.size != n_data:
        raise ValueError(
            f"{name}: forward model returned {pred.size} values but there are {n_data} data points"
        )
    return pred


def gaussian_known_variance(
    model: ForwardModel,
    data: NDArray,
    variance: float | NDArray,
    *,
    jacobian: JacobianFn | None = None,
    name: str = "gaussian_known_variance",
) -> LogLikelihood:
    """Independent Gaussian noise with a variance that is known and held fixed.

    ``log L = -0.5 * sum_i [ (d_i - m_i(theta))^2 / v_i + log(2*pi*v_i) ]``

    The normalization term is included even though it is constant in ``theta``.  It costs nothing,
    and leaving it out makes the returned value not a log-likelihood --- which matters the moment
    two noise models are compared, or an evidence integral is attempted.

    ``variance`` may be a scalar (homoscedastic) or one value per bin.  The per-bin case is the one
    that matters for a spectrum: it is the a-priori variance, fixed by the number of averages,
    and it is what allows a wrong model to fail rather than to widen its error bars.

    Args:
        model: Parameters -> prediction, same length as ``data``.
        data: Observations.
        variance: Known noise variance, scalar or per-bin.  Strictly positive.
        jacobian: Analytic Jacobian if available; central differences otherwise.
        name: Label carried into refusals and records.

    Returns:
        A :class:`LogLikelihood` whose gradient is ``J^T (r / v)`` with ``r = d - m(theta)``.
    """
    data = _as_1d(data, "data")
    var = _broadcast_variance(variance, data.size)
    jac = _jacobian_for(model, jacobian)
    const = -0.5 * float(np.sum(np.log(2.0 * np.pi * var)))

    def evaluate(theta: NDArray) -> tuple[float, NDArray]:
        pred = _predict(model, theta, data.size, name)
        resid = data - pred
        value = const - 0.5 * float(np.sum(resid * resid / var))
        grad = np.asarray(jac(theta), dtype=float).T @ (resid / var)
        return value, grad

    return LogLikelihood(
        name=name,
        n_data=int(data.size),
        n_parameters=None,
        evaluate=evaluate,
        noise_variance_is_fitted=False,
        metadata={"variance": var.copy()},
    )


def log_amplitude_variance(chi_squared_dof: int, *, exact: bool = False) -> float:
    """A-priori variance of ``log P_hat`` for a spectral bin carrying ``nu`` chi-squared dof.

    ``Var[log P_hat] = polygamma(1, nu/2)``, which is ``2/nu * (1 + 1/(2*nu) + ...)``.

    **The argument is nu, not a snapshot count**, and the difference is a factor of two that no
    amount of averaging removes.  A *real* Gaussian field contributes one dof per snapshot, so
    ``nu = M``; a *complex* mode carries two real quadratic dof, so ``nu = 2M``.  Both spectral
    likelihoods in this module take the same ``chi_squared_dof`` for exactly this reason --- they
    described the same data with a factor of two between them while each argument was spelled
    ``n_average``, and each test graded a convention against itself.  ``validation.analytic.stats``
    states the same warning on its own count argument.

    The asymptotic form is the default because it is the number the protocol quotes and the one a
    reader will check by hand.  ``exact=True`` returns the trigamma value, which differs by about
    ``1/nu`` in relative terms --- 5% at ``nu=20``, 0.5% at ``nu=200``.  At small ``nu`` prefer
    :func:`gamma_power_spectrum`, which needs no Gaussian approximation at all.
    """
    if chi_squared_dof < 1:
        raise ValueError(f"chi_squared_dof must be at least 1, got {chi_squared_dof}")
    if exact:
        return float(special.polygamma(1, chi_squared_dof / 2.0))
    return 2.0 / float(chi_squared_dof)


def gaussian_apriori_log_spectrum(
    power_model: ForwardModel,
    observed_power: NDArray,
    chi_squared_dof: int,
    *,
    jacobian: JacobianFn | None = None,
    debias: bool = True,
    exact_variance: bool = False,
    name: str = "gaussian_apriori_log_spectrum",
) -> LogLikelihood:
    """Gaussian likelihood in log-amplitude with the per-bin variance fixed a priori.

    The forward model predicts *power*; the comparison happens in ``log`` power, where the noise is
    close to Gaussian and, crucially, where its variance is known before the experiment:
    :func:`log_amplitude_variance`.  Nothing about the fit can inflate it.

    Two details that are easy to get wrong and both matter for ``(kappa, sigma)``:

    **The log of an unbiased estimator is biased.**  With ``nu = chi_squared_dof`` degrees of
    freedom, ``E[log P_hat] = log P + digamma(nu/2) - log(nu/2)``, which is negative --- about -0.10
    at ``nu=10``.  Left uncorrected it pulls the whole fitted spectrum down by a constant in log space,
    which is exactly the amplitude rescale that ``(kappa, sigma) -> (kappa/f, sigma/f)`` absorbs.  A
    bias that lands squarely on a degeneracy direction is the worst possible place for one, so
    ``debias=True`` (the default) subtracts it from the data.

    **The chain rule.**  ``d log P / d theta = (1/P) * dP/d theta``.  The Jacobian argument, if
    given, is the Jacobian of *power*, not of log power.

    Args:
        power_model: Parameters -> predicted power per bin, strictly positive.
        observed_power: Measured power per bin, strictly positive.
        chi_squared_dof: ``nu``, the chi-squared degrees of freedom in each bin. One per snapshot
            for a real field, **two** per snapshot for a complex mode. See
            :func:`log_amplitude_variance`.
        jacobian: Analytic Jacobian of the power model.
        debias: Subtract the known ``log``-of-a-chi-squared bias from the data.
        exact_variance: Use the trigamma variance instead of ``2/N``.
        name: Label carried into refusals and records.
    """
    observed = _as_1d(observed_power, "observed_power")
    if not np.all(observed > 0.0):
        raise ValueError("observed_power must be strictly positive to take its logarithm")
    if chi_squared_dof < 1:
        raise ValueError(f"chi_squared_dof must be at least 1, got {chi_squared_dof}")

    var = np.full(observed.size, log_amplitude_variance(chi_squared_dof, exact=exact_variance))
    bias = 0.0
    if debias:
        nu = float(chi_squared_dof)
        bias = float(special.digamma(nu / 2.0) - np.log(nu / 2.0))
    log_data = np.log(observed) - bias
    jac = _jacobian_for(power_model, jacobian)
    const = -0.5 * float(np.sum(np.log(2.0 * np.pi * var)))

    def evaluate(theta: NDArray) -> tuple[float, NDArray]:
        power = _predict(power_model, theta, observed.size, name)
        if not np.all(power > 0.0):
            ModelInadequate(
                reason=f"{name}: the forward model returned a non-positive power.",
                offending_input={"theta": np.asarray(theta, float).tolist()},
                lift_requires=(
                    "Constrain the parameters so the predicted spectrum stays positive, or fit in "
                    "log-parameters where positivity is structural."
                ),
                detail={"min_power": float(np.min(power))},
            ).raise_()
        resid = log_data - np.log(power)
        value = const - 0.5 * float(np.sum(resid * resid / var))
        # d/dtheta of (-0.5 r^2/v) with r = logdata - log P is (r/v) * (1/P) * dP/dtheta.
        grad = np.asarray(jac(theta), dtype=float).T @ (resid / var / power)
        return value, grad

    return LogLikelihood(
        name=name,
        n_data=int(observed.size),
        n_parameters=None,
        evaluate=evaluate,
        noise_variance_is_fitted=False,
        metadata={
            "chi_squared_dof": int(chi_squared_dof),
            "log_amplitude_variance": float(var[0]),
            "log_bias_removed": bias,
        },
    )


def gamma_power_spectrum(
    power_model: ForwardModel,
    observed_power: NDArray,
    chi_squared_dof: int,
    *,
    jacobian: JacobianFn | None = None,
    name: str = "gamma_power_spectrum",
) -> LogLikelihood:
    """The chi-squared / Wishart-flavoured likelihood for an averaged power spectrum.

    A periodogram bin averaged over ``N`` independent realizations is Gamma-distributed with shape
    ``N`` and scale ``P(theta)/N``; equivalently ``2N*P_hat/P ~ chi2_{2N}`` in the complex-mode
    convention.  This is the exact sampling distribution, so unlike the log-Gaussian form it needs
    no large-``N`` argument and stays correct in the sparsely averaged high-``q`` bins.

    ``log L = sum_k [ N log N - lgamma(N) + (N-1) log P_hat_k - N log P_k - N P_hat_k / P_k ]``

    and ``d log L / d P_k = N (P_hat_k - P_k) / P_k^2``, which vanishes at ``P = P_hat`` as it must.

    In one dimension the Gamma density *is* the complex Wishart density; a spectral density matrix
    over several coupled fields would replace ``P_hat/P`` by ``trace(P^-1 P_hat)`` and ``log P`` by
    ``logdet P`` behind this same interface.

    A non-positive predicted power raises :class:`~aleph.infer.refusal.ModelInadequateError` rather
    than returning ``-inf``.  Returning ``-inf`` lets an optimizer wander into a region where the
    model means nothing and then report where it stopped; the refusal says what broke it.
    """
    observed = _as_1d(observed_power, "observed_power")
    if not np.all(observed > 0.0):
        raise ValueError("observed_power must be strictly positive for a Gamma likelihood")
    if chi_squared_dof < 1:
        raise ValueError(f"chi_squared_dof must be at least 1, got {chi_squared_dof}")

    # shape = nu/2, scale = 2P/nu.  Stated in dof rather than in a snapshot count so that this and
    # the log-Gaussian above cannot read one argument as two different distributions.
    nn = 0.5 * float(chi_squared_dof)
    jac = _jacobian_for(power_model, jacobian)
    const = float(
        observed.size * (nn * np.log(nn) - special.gammaln(nn))
        + (nn - 1.0) * np.sum(np.log(observed))
    )

    def evaluate(theta: NDArray) -> tuple[float, NDArray]:
        power = _predict(power_model, theta, observed.size, name)
        if not np.all(power > 0.0):
            ModelInadequate(
                reason=f"{name}: the forward model returned a non-positive power.",
                offending_input={"theta": np.asarray(theta, float).tolist()},
                lift_requires=(
                    "A Gamma likelihood is defined only for positive predicted power. "
                    "Re-parametrize in logs, or bound the parameters away from the sign change."
                ),
                detail={"min_power": float(np.min(power))},
            ).raise_()
        value = const - nn * float(np.sum(np.log(power))) - nn * float(np.sum(observed / power))
        dl_dp = nn * (observed - power) / (power * power)
        grad = np.asarray(jac(theta), dtype=float).T @ dl_dp
        return value, grad

    return LogLikelihood(
        name=name,
        n_data=int(observed.size),
        n_parameters=None,
        evaluate=evaluate,
        noise_variance_is_fitted=False,
        metadata={
            "chi_squared_dof": int(chi_squared_dof),
            "distribution": "gamma(shape=nu/2, scale=2P/nu)",
        },
    )


def gaussian_profiled_variance(
    model: ForwardModel,
    data: NDArray,
    *,
    jacobian: JacobianFn | None = None,
    name: str = "gaussian_profiled_variance",
) -> LogLikelihood:
    """ANTI-PATTERN, kept so it can be demonstrated rather than described.

    This is the Gaussian likelihood with the noise variance profiled out at its maximum-likelihood
    value ``sigma^2 = mean(r^2)``.  The profiled log-likelihood is

    ``log L = -n/2 * [ log(2*pi) + 1 + log(mean(r^2)) ]``

    and its gradient is ``(n / sum(r^2)) * J^T r``, which points in exactly the same direction as
    least squares --- so the best-fit *parameters* are unchanged.  What changes is that the error
    bars and the goodness of fit are now derived from the residuals themselves.  Reduced
    chi-squared becomes 1 for every model, correct or not, and no amount of model error can ever
    produce a bad fit.

    :func:`chi_squared_report` refuses to grade a fit built on this likelihood.  Use it only to show
    what a fitted variance costs; never to report a measurement.
    """
    data = _as_1d(data, "data")
    n = float(data.size)
    jac = _jacobian_for(model, jacobian)

    def evaluate(theta: NDArray) -> tuple[float, NDArray]:
        pred = _predict(model, theta, data.size, name)
        resid = data - pred
        sse = float(np.sum(resid * resid))
        if sse <= 0.0:
            ModelInadequate(
                reason=f"{name}: residuals are exactly zero, so the profiled variance is zero.",
                offending_input={"theta": np.asarray(theta, float).tolist()},
                lift_requires="Use a known a-priori variance; a fitted variance is undefined here.",
            ).raise_()
        value = -0.5 * n * (_LOG_2PI + 1.0 + np.log(sse / n))
        grad = (n / sse) * (np.asarray(jac(theta), dtype=float).T @ resid)
        return value, grad

    return LogLikelihood(
        name=name,
        n_data=int(data.size),
        n_parameters=None,
        evaluate=evaluate,
        noise_variance_is_fitted=True,
        metadata={"warning": "variance fitted from residuals; goodness of fit is not testable"},
    )


@dataclass(frozen=True)
class ChiSquaredReport:
    """Goodness of fit against an a-priori error budget.

    Attributes:
        chi_squared: Sum of squared standardized residuals.
        dof: Degrees of freedom, ``n_data - n_parameters``.
        reduced: ``chi_squared / dof``.
        p_value: Survival function of ``chi2_{dof}`` at ``chi_squared``; small means the model
            cannot produce data like this given the stated noise.
        verdict: ``None`` when the fit is acceptable, otherwise a typed
            :class:`~aleph.infer.refusal.ModelInadequate`.
    """

    chi_squared: float
    dof: int
    reduced: float
    p_value: float
    verdict: ModelInadequate | None

    @property
    def adequate(self) -> bool:
        return self.verdict is None


def chi_squared_report(
    likelihood: LogLikelihood,
    model: ForwardModel,
    data: NDArray,
    variance: float | NDArray,
    theta: NDArray,
    *,
    n_parameters: int | None = None,
    alpha: float = 1e-3,
) -> ChiSquaredReport:
    """Grade a fit against the noise budget that was fixed before the data arrived.

    This is the check that only exists when the variance is a-priori.  If ``likelihood`` fitted its
    own variance, this function raises rather than returning a number, because the number would be
    1 by construction and would be read as agreement.

    Args:
        likelihood: The likelihood used for the fit --- inspected only for
            ``noise_variance_is_fitted``.
        model: Forward model.
        data: Observations, in the same space as ``model``'s output.
        variance: The a-priori per-bin variance.
        theta: Best-fit parameters.
        n_parameters: Fitted parameter count; defaults to ``len(theta)``.
        alpha: Rejection level for the verdict.
    """
    if likelihood.noise_variance_is_fitted:
        ModelInadequate(
            reason=(
                "Refusing to compute a goodness of fit for a likelihood whose noise variance was "
                "fitted to the residuals. Reduced chi-squared would be 1 by construction."
            ),
            offending_input={"likelihood": likelihood.name},
            lift_requires=(
                "Re-fit with an a-priori variance, e.g. 2/N per bin for a log spectrum averaged "
                "over N samples."
            ),
        ).raise_()

    data = _as_1d(data, "data")
    var = _broadcast_variance(variance, data.size)
    theta = np.asarray(theta, dtype=float)
    pred = _predict(model, theta, data.size, "chi_squared_report")
    resid = data - pred
    chi2 = float(np.sum(resid * resid / var))
    k = int(theta.size if n_parameters is None else n_parameters)
    dof = int(data.size - k)
    if dof <= 0:
        raise ValueError(
            f"no degrees of freedom left: {data.size} data points and {k} fitted parameters"
        )
    p_value = float(stats.chi2.sf(chi2, dof))
    verdict = None
    if p_value < alpha:
        verdict = ModelInadequate(
            reason=(
                f"chi-squared {chi2:.4g} on {dof} degrees of freedom (reduced {chi2 / dof:.4g}) "
                f"has p = {p_value:.3g}: the model cannot produce this data given the stated noise."
            ),
            offending_input={"theta": theta.tolist()},
            lift_requires=(
                "Either the forward model is missing physics, or the a-priori variance understates "
                "the real noise. Do not widen the error bars to fix this without evidence that the "
                "noise, not the model, was wrong."
            ),
            statistic=chi2,
            p_value=p_value,
        )
    return ChiSquaredReport(
        chi_squared=chi2, dof=dof, reduced=chi2 / dof, p_value=p_value, verdict=verdict
    )
