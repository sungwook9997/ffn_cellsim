"""Whether a time series supports quoting a mean, and what the error bar on it actually is.

This module is domain-independent. Nothing in it knows what a membrane is, and nothing in it should
learn. It is numerical-experiment hygiene: given evenly spaced samples of *anything*, decide whether
a steady state was reached, and if so report a mean with an error bar that is not a lie.

The one result everything here rests on
---------------------------------------
For a stationary series ``x_1..x_N`` with sample spacing ``dt``, variance ``sigma^2`` and normalised
autocorrelation ``rho_k``, the variance of the sample mean is *exactly*::

    Var(mean) = (sigma^2 / N) * [ 1 + 2 * sum_{k=1..N-1} (1 - k/N) * rho_k ]

The bracket is the variance inflation factor. For summable ``rho`` and large ``N`` it tends to
``1 + 2*sum_k rho_k``, so defining the integrated autocorrelation time **in samples** as::

    tau_samples = 1/2 + sum_{k>=1} rho_k

the bracket is exactly ``2 * tau_samples``, and::

    sem = sigma * sqrt(2 * tau_samples / N) = sigma * sqrt(2 * tau_int / T)

with ``tau_int = tau_samples * dt`` and ``T = N * dt``.

The familiar ``sigma / sqrt(N)`` is the special case ``tau_samples = 1/2``, i.e. an uncorrelated
series — which is why ``1/2``, not ``0``, is the floor for ``tau_samples``.

Why this is not a small correction
----------------------------------
The ratio of the two forms is ``sqrt(2 * tau_samples)``. For an AR(1) series with coefficient
``phi``, ``rho_k = phi^k``, so ``tau_samples = (1+phi) / (2*(1-phi))`` and::

    sem_correct / sem_naive = sqrt( (1+phi) / (1-phi) )

At ``phi = 0.95`` that is **6.24**. Measured over 400 realisations of ``N = 4000``: the empirical
spread of the sample mean is 6.15x the naive error bar, and the corrected one lands within 0.5% of
it. A nominal 95% confidence interval built the naive way covers the truth **22.5%** of the time;
built the correct way, 94.5%. ``tests/observe/test_stationarity.py`` runs exactly that experiment,
and it is the reason this module exists.

An error bar six times too small is not a rounding issue. It is a false claim, and it is false in
the direction that makes a result look significant.

Windowing the sum
-----------------
``rho_hat_k`` carries estimator noise of order ``1/sqrt(N)`` at every lag while the true ``rho_k``
decays. Summing all ``N-1`` lags therefore adds a random walk of ``N`` steps of size ``1/sqrt(N)``:
an ``O(1)`` contribution with no preferred sign, which in practice inflates ``tau`` because the sum
gets truncated wherever the noise happened to be running high. The Sokal automatic window takes the
smallest ``M`` with ``M >= c * tau(M)``, cutting the sum a fixed number of correlation times out, so
the discarded tail is exponentially small while the accumulated noise stays bounded.

``c = 5`` rather than the usual 6, because Aleph's series are short and a larger ``c`` frequently
never closes the window. When it does not close, that is **reported** and ``tau_int`` is labelled a
lower bound, rather than being returned as though it had converged.

The refusal criterion is self-referential
-----------------------------------------
A series is ``TOO_SHORT`` when ``T < min_tau_windows * tau_int``: when it does not contain enough
independent samples *of the thing being measured*. Since ``n_eff = T / (2 * tau_int)``, that rule is
exactly ``n_eff < min_tau_windows / 2`` — at least 25 independent samples, by default.

The alternative, testing against an externally supplied driving time scale, has a failure mode worth
naming: the driving time is a property of the forcing, and the collective relaxation of the
observable can be orders of magnitude slower. A series can clear that bar and still be hopeless for
its own correlation time. A caller who genuinely has an external bound can impose it as well, via
``min_span_s``; both criteria are recorded in the returned contract.

The limitation of that criterion, stated rather than hidden
-----------------------------------------------------------
The gate tests against the *measured* ``tau_int``, and on a short series the measurement is
**biased low**. Measured on AR(1) at ``phi = 0.99`` (true ``tau_samples = 99.5``):

    N = 1,000     tau measured 64      0.64x the truth
    N = 2,000     tau measured 56      0.56x
    N = 5,000     tau measured 98      0.98x
    N = 100,000   tau measured 109     1.10x

The bias is toward *fewer* correlation times, hence *more* apparent independent samples, hence an
error bar that is **still too small** — the same direction as the naive form, just far less
extreme. So the gate is necessary and not sufficient, and ``sem`` on a series that only just clears
it should be read as a lower bound. That is why the default is 50 correlation times rather than a
handful, why ``TauEstimate.window_closed`` is reported, and why every report carries
``tau_reliability`` (``T / tau_int``) so a reader can see how far into the reliable regime the
estimate sits. ``tests/observe/test_stationarity.py::test_tau_is_biased_low_on_a_short_series``
measures the bias directly rather than leaving this paragraph as the only record of it.

Verdict order, and why it is this order
---------------------------------------
1. an explicit ``min_span_s``, if the caller declared one;
2. the transient guard — a window that demonstrably opens on a rise;
3. the drift test;
4. the length gate;
5. otherwise ``STATIONARY``.

The drift test comes **before** the length gate deliberately. A significant slope is a positive
finding at any length, because :func:`drift_significance` is error-barred against the correlated
null and so cannot mistake a short wander for a trend. Running the gate first would hide real
results: a clean ramp has an enormous ``tau_int`` *precisely because* it is a ramp, so the gate
would call it ``TOO_SHORT`` and the caller would go away and run longer to observe the same ramp
again.

What this module does not do
----------------------------
It takes no position on whether a stationary value is *correct*. Whether a mean is meaningful at all
and whether it is the right number are different questions, and this module answers only the first.
It also assumes no fluctuation-dissipation relation: in a driven steady state the fluctuation size is
not set by temperature, so ``sigma`` is measured, never predicted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

import numpy as np

__all__ = [
    "DEFAULT_SOKAL_C",
    "EquilibrationSearch",
    "StationarityReport",
    "StationarityVerdict",
    "TauEstimate",
    "DEFAULT_DRIFT_SIGMA",
    "DEFAULT_MIN_TAU_WINDOWS",
    "DEFAULT_TRANSIENT_SIGMA",
    "assess_stationarity",
    "drift_significance",
    "effective_sample_size",
    "equilibration_start",
    "normalised_autocorrelation",
    "window_opens_on_a_transient",
    "integrated_autocorrelation_time",
    "naive_standard_error",
    "standard_error_of_correlated_mean",
    "variance_inflation_factor",
]

#: Sokal window constant. See the module docstring for why it is 5 and not 6.
DEFAULT_SOKAL_C: float = 5.0

#: Significance threshold, in standard deviations, for both the drift test and the transient guard.
#: Three sigma rather than one, because both statistics are now properly error-barred against the
#: correlated null — a one-sigma threshold on a correctly normalised statistic rejects a third of
#: perfectly stationary series.
DEFAULT_DRIFT_SIGMA: float = 3.0

#: Same threshold, for the transient guard. Named separately because they test different things and
#: a caller may reasonably want them apart.
DEFAULT_TRANSIENT_SIGMA: float = 3.0

#: A window must span this many measured correlation times before a mean may be quoted, i.e. at
#: least 25 independent samples. Set from the measured bias of the tau estimator rather than by
#: taste: below roughly 50 correlation times the estimate is materially low, which makes the error
#: bar too small in the same direction as the naive form. See the module docstring.
DEFAULT_MIN_TAU_WINDOWS: float = 50.0

#: Correlation times the trailing tenth must span before :func:`window_opens_on_a_transient` will
#: return a verdict at all. Below it the block's own ``tau`` estimate is unusable and the test
#: becomes over-sensitive; such series are ``TOO_SHORT`` and the length gate is the right reporter.
_MIN_TRANSIENT_BLOCK_TAUS: float = 10.0

#: An uncorrelated series has ``tau_samples = 1/2``. It is a floor, not a small number.
_TAU_SAMPLES_FLOOR: float = 0.5


class StationarityVerdict(StrEnum):
    """Whether a series supports quoting a mean.

    ``STATIONARY``
        The trend across the analysed window is smaller than the window's own fluctuation, and the
        window is long enough relative to the measured correlation time for that comparison to mean
        something. A mean and an error bar are returned.

    ``DRIFTING``
        Long enough to judge, and still going somewhere. A mean here is the average of a transient.
        No mean is returned.

    ``TOO_SHORT``
        Not enough independent samples to tell the two apart. Kept distinct from ``DRIFTING`` on
        purpose: "we looked and it moved" and "we did not look long enough" are different findings,
        and merging them lets a fraction of a relaxation time be reported as a measurement.
    """

    STATIONARY = "STATIONARY"
    DRIFTING = "DRIFTING"
    TOO_SHORT = "TOO_SHORT"


@dataclass(frozen=True, slots=True)
class TauEstimate:
    """The correlation-time estimate, with the honesty flags that decide how far to trust it.

    Attributes:
        tau_int_s: Integrated autocorrelation time in the time unit ``dt`` was given in.
        tau_samples: The same quantity in samples: ``tau_int_s / dt``. Floored at ``1/2``.
        n_effective: ``T / (2 * tau_int_s) = N / (2 * tau_samples)``. The number of uncorrelated
            samples that would give the same error on the mean.
        window_lag: The lag at which the Sokal sum was cut.
        window_closed: Whether the criterion ``M >= c*tau(M)`` was ever satisfied. When ``False``,
            ``tau_int_s`` is a **lower** bound and ``n_effective`` an upper one — the series is short
            relative to its own correlation time and cannot say how short.
        reliability: ``T / tau_int``, the number of correlation times the series spans. Below
            roughly 50 the ``tau_int`` estimate is materially biased low, so ``n_effective`` is an
            over-estimate and any error bar built on it is a lower bound. Reported so that fact is
            visible in the record rather than only in the module docstring.
        zero_variance: A constant series. ``tau_samples`` is floored and the report says so, rather
            than the estimator dividing by zero.
        note: What happened, in words.
    """

    tau_int_s: float
    tau_samples: float
    n_effective: float
    window_lag: int
    window_closed: bool
    reliability: float
    zero_variance: bool
    note: str


@dataclass(frozen=True, slots=True)
class EquilibrationSearch:
    """Where the transient was judged to end, and how confident that judgement is.

    Attributes:
        start_index: First index of the retained window.
        start_time_s: ``start_index * dt``.
        n_effective: Effective sample size achieved by the retained window — the quantity that was
            maximised.
        opens_on_transient: Whether the retained window *still* begins on a rise. See
            :func:`window_opens_on_a_transient` for why this can happen and why it is dangerous.
        transient_ratio: The statistic that flagged it, in standard deviations.
        n_contaminated_skipped: How many candidate truncations were REJECTED because their window
            opened on a rise. ⚠ **This field exists because filtering turned the warning light off.**
            Before 2026-08-22 ``opens_on_transient`` was computed once, on the winner, so a
            contaminated series announced itself. Once contaminated candidates are SKIPPED the winner
            is almost never contaminated, and the flag then reads clean on a series that is not — the
            contamination did not go away, the search walked around it. Two recorded tests caught
            exactly that. This counts what was walked around, so the signal MOVES instead of vanishing.
            ⚠ **The COUNT is the signal; a non-zero value is not.** A settled series skips candidates
            too — measured on the calm fixtures in ``tests/world/test_observe_gamma.py``, 3-12 of them,
            against 28 on the contaminated pair. Early truncations of any series open near its start,
            where a rise is easiest to find. Read this against a comparable series, never as a flag.
        fell_back_to_contaminated: True when NO candidate was clean and the best contaminated one was
            returned anyway. ``opens_on_transient`` is then True and means what it always meant.
    """

    start_index: int
    start_time_s: float
    n_effective: float
    opens_on_transient: bool
    transient_ratio: float
    n_contaminated_skipped: int = 0
    fell_back_to_contaminated: bool = False


@dataclass(frozen=True, slots=True)
class StationarityReport:
    """One series' verdict, carrying everything needed to re-judge it without the raw data.

    Attributes:
        verdict: See :class:`StationarityVerdict`.
        observable: Name, for the record.
        mean: Steady-state mean over the analysed window. ``None`` unless ``STATIONARY`` — a mean
            withheld is the finding, and returning it "for information" is how it gets quoted.
        sem: ``sigma * sqrt(2*tau_int/T)``, from *independent* samples. ``None`` with the mean.
        sem_naive: ``sigma / sqrt(N)``. Reported **always**, including when it is wrong, so the
            size of the correction is visible in every record rather than only in this docstring.
        sem_inflation: ``sem / sem_naive = sqrt(2 * tau_samples)``. One is the uncorrelated case.
        z_ess: ``(mean - null_value) / sem``. The defensible significance.
        z_naive: ``(mean - null_value) / sem_naive``. The one a naive analysis would have quoted.
            Both are reported so the difference is impossible to miss.
        std: Sample standard deviation of the analysed window.
        n_samples: Samples in the analysed window.
        tau: The correlation-time estimate.
        equilibration: Where the transient was judged to end.
        analysed_time_s: Physical span of the analysed window.
        drift_per_s: **Signed** least-squares slope. Signed on purpose: a quantity falling toward
            its steady state and one rising toward it are different physics, and ``abs`` hides
            which happened from the person who has to decide what to do next.
        drift_sigma_measured: The drift across the window in units of its own null standard
            deviation, ``|slope*T| / (sigma*sqrt(12/n_eff))`` — see :func:`drift_significance`. This
            is the quantity the verdict is thresholded on, and it is a *significance*, not a ratio to
            the scatter. The distinction is the same one this whole module is about: a ratio to the
            scatter has no null distribution, so a threshold on it means nothing on a correlated
            series.
        contract: The caller's pre-declared thresholds, copied here so a later reader judges the
            verdict against the criterion that produced it rather than against today's.
        notes: What happened, in words.
    """

    verdict: StationarityVerdict
    observable: str
    mean: float | None
    sem: float | None
    sem_naive: float
    sem_inflation: float
    z_ess: float | None
    z_naive: float | None
    std: float
    n_samples: int
    tau: TauEstimate | None
    equilibration: EquilibrationSearch | None
    analysed_time_s: float
    drift_per_s: float | None
    drift_sigma_measured: float | None
    contract: dict[str, Any] = field(default_factory=dict)
    notes: dict[str, Any] = field(default_factory=dict)

    @property
    def refused(self) -> bool:
        """``True`` when no mean may be quoted. ``TOO_SHORT`` and ``DRIFTING`` both refuse."""
        return self.verdict is not StationarityVerdict.STATIONARY

    def as_dict(self) -> dict[str, Any]:
        """JSON-able form, for a run record."""
        return {
            "verdict": self.verdict.value,
            "observable": self.observable,
            "mean": self.mean,
            "sem": self.sem,
            "sem_naive": self.sem_naive,
            "sem_inflation": self.sem_inflation,
            "z_ess": self.z_ess,
            "z_naive": self.z_naive,
            "std": self.std,
            "n_samples": self.n_samples,
            "tau_int_s": None if self.tau is None else self.tau.tau_int_s,
            "tau_samples": None if self.tau is None else self.tau.tau_samples,
            "n_effective": None if self.tau is None else self.tau.n_effective,
            "sokal_window_closed": None if self.tau is None else self.tau.window_closed,
            "tau_reliability": None if self.tau is None else self.tau.reliability,
            "equilibration_time_s": (
                None if self.equilibration is None else self.equilibration.start_time_s
            ),
            "analysed_time_s": self.analysed_time_s,
            "drift_per_s": self.drift_per_s,
            "drift_sigma_measured": self.drift_sigma_measured,
            "contract": dict(self.contract),
            "notes": dict(self.notes),
        }


# ---------------------------------------------------------------------------------------------
# Estimators
# ---------------------------------------------------------------------------------------------


def _as_series(series: Any, *, minimum: int = 4) -> np.ndarray:
    x = np.asarray(series, dtype=np.float64).ravel()
    if x.size < minimum:
        raise ValueError(
            f"need at least {minimum} samples to say anything about a correlation time; "
            f"got {x.size}"
        )
    if not np.all(np.isfinite(x)):
        raise ValueError(
            "series contains non-finite samples. They are refused rather than dropped: a NaN in a "
            "time series is a gap, and silently closing a gap changes the sample spacing that every "
            "quantity here is expressed in"
        )
    return x


def _check_dt(dt: float) -> float:
    if not (isinstance(dt, int | float) and math.isfinite(dt) and dt > 0.0):
        raise ValueError(
            f"dt must be a positive finite sample interval in seconds; got {dt!r}. It is required "
            "and has no default on purpose: passing an iteration index in its place returns a "
            "correlation time measured in iterations, which is exactly the confusion this module "
            "exists to prevent"
        )
    return float(dt)


def normalised_autocorrelation(series: np.ndarray) -> np.ndarray:
    """Biased normalised autocorrelation ``rho_k``, ``k = 0..N-1``, by FFT.

    The mean is subtracted first. Without that this measures the mean rather than the fluctuation
    about it, and for a series with a large offset ``rho_k`` would be near 1 at every lag.

    The estimator is the *biased* one (dividing by ``N`` rather than by ``N-k``). That is
    deliberate: the unbiased form has variance that grows without bound as ``k`` approaches ``N``,
    where only a handful of pairs contribute, and feeding those lags into a sum is how a correlation
    time acquires a long random tail. The biased form tapers them toward zero, which is the right
    behaviour for something that is about to be summed.

    Args:
        series: Finite samples, evenly spaced.

    Returns:
        ``rho`` of length ``N`` with ``rho[0] == 1``.
    """
    x = np.asarray(series, dtype=np.float64).ravel()
    n = x.size
    y = x - x.mean()
    # Zero-pad to at least 2N-1 so the circular correlation the FFT computes equals the linear one.
    size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(y, size)
    acf = np.fft.irfft(spectrum * np.conjugate(spectrum), size)[:n].real
    if acf[0] <= 0.0:
        return np.zeros(n, dtype=np.float64)
    return acf / acf[0]


def integrated_autocorrelation_time(
    series: Any, dt: float, *, c: float = DEFAULT_SOKAL_C
) -> TauEstimate:
    """Estimate ``tau_int`` with a Sokal automatic window.

    ``tau_samples = 1/2 + sum_{k=1..M} rho_k``, with ``M`` the smallest lag satisfying
    ``M >= c * tau_samples(M)``. See the module docstring for why the sum must be windowed at all.

    Args:
        series: Finite samples, evenly spaced. At least 4.
        dt: Sample interval [s]. Required, and unitful — see :func:`_check_dt`.
        c: Sokal window constant.

    Returns:
        A :class:`TauEstimate`. When the window never closed, ``window_closed`` is ``False`` and
        the returned ``tau_int_s`` is a lower bound on the truth.

    Raises:
        ValueError: If ``dt`` is not a positive finite number, if there are fewer than 4 samples,
            or if any sample is non-finite.
    """
    step = _check_dt(dt)
    x = _as_series(series)
    n = x.size
    total_time = n * step

    variance = float(np.var(x))
    if variance <= 0.0:
        # A constant series is perfectly "stationary" and carries no information about correlation.
        # Floor tau and say so, rather than dividing by a zero autocorrelation at lag 0.
        return TauEstimate(
            tau_int_s=_TAU_SAMPLES_FLOOR * step,
            tau_samples=_TAU_SAMPLES_FLOOR,
            n_effective=float(n),
            window_lag=0,
            window_closed=True,
            reliability=float(n) / _TAU_SAMPLES_FLOOR,
            zero_variance=True,
            note="series has zero variance; tau floored at the uncorrelated value dt/2",
        )

    rho = normalised_autocorrelation(x)
    running = np.cumsum(rho[1:])  # running[m-1] = sum_{k=1..m} rho_k

    tau_samples = _TAU_SAMPLES_FLOOR
    window = n - 1
    closed = False
    for m in range(1, n):
        tau_samples = _TAU_SAMPLES_FLOOR + float(running[m - 1])
        if m >= c * tau_samples:
            window = m
            closed = True
            break

    tau_samples = max(tau_samples, _TAU_SAMPLES_FLOOR)
    tau_int_s = tau_samples * step
    n_eff = max(1.0, total_time / (2.0 * tau_int_s))

    if closed:
        note = f"Sokal window closed at lag {window} (c={c:g})"
    else:
        note = (
            f"Sokal window never closed within {n} samples: the series is short relative to its own "
            "correlation time, so tau_int is a LOWER bound and n_effective an UPPER one"
        )

    return TauEstimate(
        tau_int_s=tau_int_s,
        tau_samples=tau_samples,
        n_effective=n_eff,
        window_lag=int(window),
        window_closed=closed,
        reliability=float(n) / tau_samples,
        zero_variance=False,
        note=note,
    )


def variance_inflation_factor(tau_samples: float) -> float:
    """``2 * tau_samples`` — the factor by which correlation inflates ``Var(mean)``.

    One for an uncorrelated series. This is the bracket of the module docstring's exact formula, in
    its large-``N`` limit.
    """
    return 2.0 * float(tau_samples)


def effective_sample_size(n_samples: int, tau_samples: float) -> float:
    """``N / (2 * tau_samples)`` — the uncorrelated sample count with the same error on the mean.

    Never above ``n_samples``, with equality only in the uncorrelated case, because ``tau_samples``
    is floored at ``1/2``.
    """
    return max(1.0, float(n_samples) / (2.0 * max(float(tau_samples), _TAU_SAMPLES_FLOOR)))


def drift_significance(
    window: Any, dt: float, *, c: float = DEFAULT_SOKAL_C
) -> tuple[float, float, float]:
    """Test a window for a trend. Returns ``(slope_per_s, significance_sigma, n_eff_residual)``.

    Derivation
    ----------
    For an ordinary least-squares fit over evenly spaced times spanning ``T``, the slope's variance
    is ``Var(b) = s^2 / (N * Var(t))`` with ``Var(t) = T^2/12``, so the *drift across the window*
    ``b*T`` has standard deviation::

        sd(b*T) = s * T / sqrt(N * T^2/12) = s * sqrt(12 / N)

    where ``s`` is the noise about the fitted line and ``N`` must be the number of **independent**
    samples of that noise. Hence::

        drift_sigma = |b * T| / ( s * sqrt(12 / n_eff) )

    Two choices in that formula do the work, and both are easy to get wrong.

    **The scale must be the residual scatter, not the series scatter.** Using the full series'
    ``sigma`` and its ``tau_int`` conflates the trend being tested with the noise it is tested
    against. For a strong trend the trend itself dominates both — it inflates ``sigma`` and it
    inflates ``tau_int`` enormously, since a monotone series is correlated at every lag — and the
    statistic goes *blind to exactly the thing it is looking for*. A clean ramp with white noise
    measured that way scores around 2 sigma, indistinguishable from a stationary wander. Detrending
    first fixes this at the root: the residuals of a ramp with white noise *are* white noise, so
    ``n_eff`` is the full sample count and the ramp scores in the hundreds of sigma.

    **The count must be ``n_eff``, not ``N``.** For a correlated series the residuals are correlated
    too, and using ``N`` would understate ``sd(b*T)`` by ``sqrt(2*tau_samples)`` — the same factor,
    and the same mistake, that this module exists to prevent for the mean.

    What this is *not*: ``|b*T| / sigma``, the ratio of the trend to the scatter. That has no null
    distribution. It says how large a trend is relative to the noise, which is a description rather
    than a test, and a stationary series with ``n_eff = 4`` produces values around 1.7 routinely.
    A threshold on that ratio therefore calls most short correlated series drifting, which converts
    an honest "not enough independent samples to tell" into a confident and wrong claim.

    Args:
        window: The samples to test.
        dt: Sample interval [s].
        c: Sokal window constant for the residual correlation-time estimate.

    Returns:
        The signed slope per second, the significance in sigma, and the residual effective sample
        count. The slope is signed because rising toward a steady state and falling toward one are
        different physics.
    """
    x = np.asarray(window, dtype=np.float64).ravel()
    step = _check_dt(dt)
    n = x.size
    if n < 4:
        return (0.0, 0.0, float(n))
    times = np.arange(n, dtype=np.float64) * step
    slope, intercept = np.polyfit(times, x, 1)
    residual = x - (slope * times + intercept)
    residual_std = float(np.std(residual, ddof=2))
    if not math.isfinite(residual_std) or residual_std <= 0.0:
        return (float(slope), 0.0, float(n))
    residual_tau = integrated_autocorrelation_time(residual, step, c=c)
    n_eff = residual_tau.n_effective
    null_sd = residual_std * math.sqrt(12.0 / max(n_eff, 1.0))
    if null_sd <= 0.0:
        return (float(slope), 0.0, n_eff)
    return (float(slope), abs(float(slope) * n * step) / null_sd, n_eff)


def naive_standard_error(std: float, n_samples: int) -> float:
    """``sigma / sqrt(N)``. Correct **only** for uncorrelated samples.

    Provided so it can be reported alongside the honest one and so the counterexample test can
    call it by name. It is not a fallback and nothing in this module uses it as one.
    """
    return float(std) / math.sqrt(float(n_samples))


def standard_error_of_correlated_mean(std: float, n_samples: int, tau_samples: float) -> float:
    """``sigma * sqrt(2 * tau_samples / N)``. The error on the mean of a correlated series."""
    return float(std) * math.sqrt(
        2.0 * max(float(tau_samples), _TAU_SAMPLES_FLOOR) / float(n_samples)
    )


#: Samples the retained window must keep, so the search cannot truncate to a tail too short for the
#: statistic it is choosing by. ⚠ Declared for a reason independent of any series: below a few hundred
#: samples an integrated autocorrelation time is not worth searching TO, whatever it comes out as. It
#: is deliberately NOT tuned — the three native series whose defect motivated the change retain
#: 6,000-11,000 samples at their optimum, so this floor is nowhere near binding on them, and a floor
#: chosen to admit the cuts one wants is a threshold chosen after seeing the data.
MIN_RETAINED_SAMPLES = 512


def equilibration_start(
    series: Any,
    dt: float,
    *,
    n_candidates: int = 40,
    c: float = DEFAULT_SOKAL_C,
    transient_sigma: float = DEFAULT_TRANSIENT_SIGMA,
) -> EquilibrationSearch:
    """Choose the transient discard point by maximising the retained window's effective sample size.

    Discard too little and the transient inflates both ``sigma`` and ``tau_int``; discard too much
    and ``T`` shrinks. ``n_eff = T / (2 * tau_int)`` trades those against each other, and taking its
    maximum makes the discard point an *output of a rule* rather than something read off a plot.

    ⚠ **THE SEARCH RANGE WAS THE FIRST HALF UNTIL 2026-08-22, AND THAT BOUND WAS BINDING.** The
    original argument for it is quoted here because it is not wrong: *"searching further would allow a
    rule to keep a quarter of the data because that quarter happened to look calm."* What made it a
    defect is that a real series' transient ran PAST the halfway mark — gamma climbing and still rising
    at the midpoint over 30,000 samples — so the search returned the furthest point it could reach, the
    retained window still opened on the rise, and every tau-derived field described the transient
    rather than a correlation time. The chosen cut was the LAST GRID CANDIDATE on two of three seeds
    and index 0 on the third: three endpoints, no interior optimum, which is what a censored search
    looks like.

    The bound is now a floor on the RETAINED TAIL rather than a cap on the DISCARD. Measured before
    changing it: the criterion self-limits at the far end, because ``n_eff`` falls with ``T`` — over
    three native series the argmax sat at 19,000-24,500 of 30,000 and every candidate past 27,000
    scored 6-29 against those maxima of 53-64. So the feared "calm quarter" does not win on these
    series.

    ⚠ **AND THAT IS NOT THE SAME AS THE CONCERN BEING ANSWERED.** A series whose late quarter is
    genuinely calm, with a spuriously small tau, would still win, and the transient filter below does
    NOT catch it — that filter rejects a window opening on a RISE, not one opening on a lull. The
    honest state is: the bound that used to hide this is gone, one failure mode is now filtered, and a
    second remains unguarded. It is guarded downstream only by ``min_tau_windows``, which a spurious
    tau also inflates.

    ⚠ **The transient flag is now a FILTER, not a label.** It used to be computed once, on the winner,
    AFTER the loop — so it never participated in the choice. Measured consequence on a real seed:
    widening the range alone made it score the HIGHEST reliability of three seeds from a window its
    own flag rejected. Filtering costs nothing on a clean series (the same cut wins) and refuses the
    contaminated one. If NO candidate is clean the best contaminated one is returned with the flag
    set, which is the previous behaviour and is reported rather than hidden.

    ⚠ **The second copy at ``aleph/engine/observe/stationarity.py`` is NOT changed.** It moved to the
    ``archive-readonly`` lane on 2026-08-22, so "widen both copies" is no longer executable and no
    longer needed: the archived copy is history. ``test_stationarity_modules_agree.py`` pins
    ``integrated_autocorrelation_time`` and the differing DEFAULTS, neither of which this touches.

    Args:
        series: Finite samples, evenly spaced.
        dt: Sample interval [s].
        n_candidates: How many truncation points to try.
        c: Sokal window constant passed through to the estimator.
        transient_sigma: Significance threshold for :func:`window_opens_on_a_transient`.

    Returns:
        An :class:`EquilibrationSearch`, including the contamination flag described in
        :func:`window_opens_on_a_transient`.
    """
    step = _check_dt(dt)
    x = _as_series(series)

    limit = max(1, x.size - MIN_RETAINED_SAMPLES)
    stride = max(1, limit // max(1, int(n_candidates)))

    # Two passes over one candidate list, so "best clean" and "best at all" are both known and the
    # fallback is a stated branch rather than an accident of iteration order.
    clean: tuple[float, int] = (-1.0, 0)
    n_skipped = 0
    any_: tuple[float, int] = (-1.0, 0)
    for start in range(0, limit, stride):
        tail = x[start:]
        if tail.size < 4:
            break
        estimate = integrated_autocorrelation_time(tail, step, c=c)
        if estimate.n_effective > any_[0]:
            any_ = (estimate.n_effective, start)
        if window_opens_on_a_transient(x, start, step, n_sigma=transient_sigma, c=c)[0]:
            n_skipped += 1
            continue
        if estimate.n_effective > clean[0]:
            clean = (estimate.n_effective, start)

    fell_back = clean[0] < 0.0
    best_n_eff, best_index = any_ if fell_back else clean
    contaminated, ratio = window_opens_on_a_transient(
        x, best_index, step, n_sigma=transient_sigma, c=c
    )
    return EquilibrationSearch(
        start_index=int(best_index),
        start_time_s=float(best_index) * step,
        n_effective=max(best_n_eff, 0.0),
        opens_on_transient=contaminated,
        transient_ratio=ratio,
        n_contaminated_skipped=n_skipped,
        fell_back_to_contaminated=fell_back,
    )


def window_opens_on_a_transient(
    series: Any,
    start_index: int,
    dt: float,
    *,
    n_sigma: float = DEFAULT_TRANSIENT_SIGMA,
    c: float = DEFAULT_SOKAL_C,
) -> tuple[bool, float]:
    """Does the retained window still begin on a rise? ``(flag, significance_in_sigma)``.

    Why this guard is needed
    ------------------------
    :func:`equilibration_start` maximises ``n_eff = T / (2*tau_int)``, on the premise that leaving a
    transient in inflates ``tau_int``. That premise fails on a common shape: a series that rises
    quickly and then sits flat for a long time can score its maximum ``n_eff`` at ``start = 0``,
    because ``T`` is maximal there and the autocorrelation estimator does not penalise a short
    leading trend enough to compensate.

    The failure is worse than it looks, because it is **self-reinforcing**. A drift verdict of the
    form ``|drift| < k * sigma`` is normalised by the window's own scatter, and keeping the transient
    inflates that scatter — so a *larger* transient makes the drift test *easier* to pass. Left
    unguarded, this combination certifies a rise as a steady state, and the bigger the rise the more
    confidently it does so.

    The statistic, and why it is this one
    -------------------------------------
    A transient is a property of a *segment*, not of a point. Testing the first sample against the
    window's scatter fails on any cut placed near the knee, since that cut necessarily leaves its
    first sample sitting in the tail of the rise. So this compares the **mean of the leading tenth**
    against the **mean of the trailing tenth**.

    The denominator has to satisfy two requirements at once, and getting either wrong breaks the
    test in a different direction.

    **It must be correlation-aware.** Normalising by the raw *scatter* of a block is the very error
    this module exists to prevent: for a correlated series the scatter of a block of ``m`` samples
    says almost nothing about the uncertainty in that block's **mean**. Each block mean has standard
    deviation ``sigma * sqrt(2*tau_samples/m)`` (the module's central result), so the difference of
    two independent block means has standard deviation::

        sd(lead - trail) = sigma * sqrt(4 * tau_samples / m)

    Using the raw scatter makes the test wildly over-sensitive on exactly the series this module
    cares about: a strongly correlated random walk wanders by many block-scatters between its ends
    without drifting at all, and would be flagged as a transient every time. Since the correct
    verdict for such a series is ``TOO_SHORT``, an over-sensitive transient test would replace an
    honest "we could not tell" with a confident and wrong "it is drifting".

    **It must be transient-proof.** Both ``sigma`` and ``tau_samples`` are therefore measured on the
    **trailing tenth alone**, which is transient-free by construction. Measuring them on the whole
    window would restore the self-reinforcing loop in a new place: a transient inflates the window's
    ``tau_int`` enormously (a monotone segment is correlated at every lag), the denominator grows
    with it, and a *larger* transient becomes *harder* to detect. Estimating from a segment the
    transient cannot reach closes that loop for good.

    Args:
        series: The full series.
        start_index: The discard point under test.
        dt: Sample interval [s], needed for the trailing-tenth correlation-time estimate.
        n_sigma: Significance threshold.
        c: Sokal window constant for that estimate.

    Returns:
        ``(True, significance)`` when the retained window still opens on a transient, else
        ``(False, significance)``. A window too short to split into tenths returns ``(False, 0.0)``:
        it cannot support the test, and inventing a verdict for it would be worse than declining.
    """
    step = _check_dt(dt)
    x = np.asarray(series, dtype=np.float64).ravel()
    tail = x[int(start_index) :]
    if tail.size < 40:
        # Below this the trailing tenth has too few samples to carry a correlation-time estimate.
        return (False, 0.0)
    tenth = max(4, tail.size // 10)
    lead = float(np.mean(tail[:tenth]))
    trailing = tail[-tenth:]
    trail = float(np.mean(trailing))

    sigma = float(np.std(trailing, ddof=1))
    if not math.isfinite(sigma) or sigma <= 0.0:
        return (False, 0.0)
    tau = integrated_autocorrelation_time(trailing, step, c=c).tau_samples

    # The trailing tenth must itself contain enough independent samples for its tau estimate to
    # mean anything. Below this the estimate collapses toward the uncorrelated floor, the
    # denominator shrinks with it, and the test starts calling ordinary correlated wander a
    # transient. Declining here is correct: such a series is TOO_SHORT, and the length gate will
    # say so. A test that cannot be supported must not return a verdict anyway.
    if tenth / max(tau, _TAU_SAMPLES_FLOOR) < _MIN_TRANSIENT_BLOCK_TAUS:
        return (False, 0.0)

    sd_difference = sigma * math.sqrt(4.0 * tau / tenth)
    if sd_difference <= 0.0:
        return (False, 0.0)
    significance = abs(lead - trail) / sd_difference
    return (significance > float(n_sigma), significance)


# ---------------------------------------------------------------------------------------------
# The verdict
# ---------------------------------------------------------------------------------------------


def _refusal_report(
    *,
    verdict: StationarityVerdict,
    observable: str,
    window: np.ndarray,
    dt: float,
    tau: TauEstimate | None,
    equilibration: EquilibrationSearch | None,
    drift_per_s: float | None,
    drift_sigma_measured: float | None,
    contract: dict[str, Any],
    reason: str,
) -> StationarityReport:
    std = float(np.std(window, ddof=1)) if window.size > 1 else 0.0
    n = int(window.size)
    sem_naive = naive_standard_error(std, n) if n else 0.0
    inflation = math.sqrt(variance_inflation_factor(tau.tau_samples)) if tau is not None else 1.0
    return StationarityReport(
        verdict=verdict,
        observable=observable,
        mean=None,
        sem=None,
        sem_naive=sem_naive,
        sem_inflation=inflation,
        z_ess=None,
        z_naive=None,
        std=std,
        n_samples=n,
        tau=tau,
        equilibration=equilibration,
        analysed_time_s=float(n) * dt,
        drift_per_s=drift_per_s,
        drift_sigma_measured=drift_sigma_measured,
        contract=contract,
        notes={"reason": reason},
    )


def assess_stationarity(
    series: Any,
    dt: float,
    *,
    observable: str,
    min_tau_windows: float = DEFAULT_MIN_TAU_WINDOWS,
    drift_sigma: float = DEFAULT_DRIFT_SIGMA,
    transient_sigma: float = DEFAULT_TRANSIENT_SIGMA,
    null_value: float = 0.0,
    min_span_s: float | None = None,
    c: float = DEFAULT_SOKAL_C,
) -> StationarityReport:
    """Judge one series, returning a mean only if a steady state was actually reached.

    The thresholds are the caller's contract. They must be declared before the run and they are
    copied into the report, because choosing them after seeing the series is how a transient gets
    certified as a steady state.

    Args:
        series: Finite samples, evenly spaced in time. At least 4.
        dt: Sample interval [s]. This is the *sampling* interval, not the integrator step: with a
            recording stride the two differ, and handing the integrator step to this estimator
            understates ``tau_int`` by exactly that stride — which makes ``n_eff`` and the error bar
            too good in the dangerous direction.
        observable: Name, recorded in the report.
        min_tau_windows: A window must span at least this many measured correlation times.
            Equivalent to requiring ``n_eff >= min_tau_windows / 2``.
        drift_sigma: Stationary requires the drift across the window to be below this many standard
            deviations of its own null distribution — see :func:`drift_significance`. This is a
            *significance* threshold, not a ratio to the scatter.
        transient_sigma: Threshold for :func:`window_opens_on_a_transient`.
        null_value: The value ``z_ess`` and ``z_naive`` are measured against.
        min_span_s: Optional *additional* absolute span requirement [s], for a caller who has an
            external time scale (a turnover time, a relaxation time). Recorded in the contract.
        c: Sokal window constant.

    Returns:
        A :class:`StationarityReport`. ``mean`` and ``sem`` are ``None`` unless the verdict is
        ``STATIONARY``.

    Raises:
        ValueError: If ``dt`` is not a positive finite number, if there are fewer than 4 samples, or
            if any sample is non-finite. These are programming errors in the caller, not findings
            about the data, which is why they raise rather than returning a verdict.
    """
    step = _check_dt(dt)
    x = _as_series(series)
    n_total = x.size

    contract: dict[str, Any] = {
        "min_tau_windows": float(min_tau_windows),
        "implied_min_n_effective": float(min_tau_windows) / 2.0,
        "drift_sigma": float(drift_sigma),
        "transient_sigma": float(transient_sigma),
        "drift_basis": "|slope*T| / (sigma*sqrt(12/n_eff)) — a significance, not a ratio to sigma",
        "null_value": float(null_value),
        "min_span_s": None if min_span_s is None else float(min_span_s),
        "sokal_c": float(c),
        "declared": "before the run; recorded so the verdict is judged against its own criterion",
        "sem_basis": "sigma*sqrt(2*tau_int/T), never sigma/sqrt(N)",
    }

    # An absolute span requirement, when the caller supplied one, is checked on the full series
    # before anything else. It is the cheap test and it needs no estimate.
    if min_span_s is not None and n_total * step < float(min_span_s):
        return _refusal_report(
            verdict=StationarityVerdict.TOO_SHORT,
            observable=observable,
            window=x,
            dt=step,
            tau=None,
            equilibration=None,
            drift_per_s=None,
            drift_sigma_measured=None,
            contract=contract,
            reason=(
                f"{n_total * step:.4g} s of series against a caller-declared minimum span of "
                f"{float(min_span_s):.4g} s"
            ),
        )

    search = equilibration_start(x, step, c=c, transient_sigma=transient_sigma)

    if search.opens_on_transient:
        # Checked before the length test on purpose. A window that demonstrably opens on a rise is
        # a positive finding — "we looked and it was still going somewhere" — and the statistic that
        # found it is already error-barred against the correlated null, so it is not an artefact of
        # a short series. Reporting TOO_SHORT here would hide a real result behind a bookkeeping
        # rule.
        kept = x[search.start_index :]
        kept_tau = integrated_autocorrelation_time(kept, step, c=c) if kept.size >= 4 else None
        slope, z, _ = drift_significance(kept, step, c=c)
        return _refusal_report(
            verdict=StationarityVerdict.DRIFTING,
            observable=observable,
            window=kept,
            dt=step,
            tau=kept_tau,
            equilibration=search,
            drift_per_s=slope,
            drift_sigma_measured=z,
            contract=contract,
            reason=(
                f"the equilibration search kept t >= {search.start_time_s:.4g} s, but that window "
                f"still OPENS on a transient: its leading tenth sits "
                f"{search.transient_ratio:.2f} sigma from its trailing tenth, measured against the "
                "standard error of that difference under the correlated null. Averaging across a "
                "rise is what this detector exists to prevent, and keeping one also inflates the "
                "scatter, which would make a naively normalised drift test easier to pass — the "
                "verdict would have been self-reinforcing. Run longer, or sample the transient more "
                "finely so the search has candidates inside it"
            ),
        )

    tail = x[search.start_index :]
    if tail.size < 4:
        return _refusal_report(
            verdict=StationarityVerdict.TOO_SHORT,
            observable=observable,
            window=tail,
            dt=step,
            tau=None,
            equilibration=search,
            drift_per_s=None,
            drift_sigma_measured=None,
            contract=contract,
            reason=(
                f"after discarding {search.start_time_s:.4g} s of transient, {tail.size} samples "
                "remain; the estimator needs at least 4"
            ),
        )

    tau = integrated_autocorrelation_time(tail, step, c=c)
    analysed_s = float(tail.size) * step
    std = float(np.std(tail, ddof=1))
    n_kept = int(tail.size)

    sem_naive = naive_standard_error(std, n_kept)
    sem_ess = standard_error_of_correlated_mean(std, n_kept, tau.tau_samples)
    inflation = math.sqrt(variance_inflation_factor(tau.tau_samples))

    slope, drift_sigma_measured, drift_n_eff = drift_significance(tail, step, c=c)

    # The drift test runs BEFORE the length gate, and the order is a claim. A significant slope is a
    # positive finding at any length, because :func:`drift_significance` is already error-barred
    # against the correlated null — it cannot mistake a short wander for a trend, so it does not
    # need the length gate's protection. Running the gate first would instead hide a real result:
    # a clean ramp has an enormous tau_int precisely *because* it is a ramp, so the length gate
    # would call it TOO_SHORT and a caller would go away and run for longer to observe the same
    # ramp again.
    if drift_sigma_measured >= float(drift_sigma):
        return StationarityReport(
            verdict=StationarityVerdict.DRIFTING,
            observable=observable,
            mean=None,
            sem=None,
            sem_naive=sem_naive,
            sem_inflation=inflation,
            z_ess=None,
            z_naive=None,
            std=std,
            n_samples=n_kept,
            tau=tau,
            equilibration=search,
            analysed_time_s=analysed_s,
            drift_per_s=slope,
            drift_sigma_measured=drift_sigma_measured,
            contract=contract,
            notes={
                "reason": (
                    f"the analysed window drifts at {drift_sigma_measured:.2f} sigma, at or above "
                    f"the declared {float(drift_sigma):g}. The significance is measured against "
                    f"s_residual*sqrt(12/n_eff) with n_eff = {drift_n_eff:.1f} independent samples "
                    "of the residual about the fitted line, not against the raw scatter. The slope "
                    f"is reported signed ({slope:+.4g} per s) because rising toward a steady state "
                    "and falling toward one are different physics"
                ),
                "tau_estimate": tau.note,
            },
        )

    required_s = float(min_tau_windows) * tau.tau_int_s
    if analysed_s < required_s:
        return StationarityReport(
            verdict=StationarityVerdict.TOO_SHORT,
            observable=observable,
            mean=None,
            sem=None,
            sem_naive=sem_naive,
            sem_inflation=inflation,
            z_ess=None,
            z_naive=None,
            std=std,
            n_samples=n_kept,
            tau=tau,
            equilibration=search,
            analysed_time_s=analysed_s,
            drift_per_s=slope,
            drift_sigma_measured=drift_sigma_measured,
            contract=contract,
            notes={
                "reason": (
                    f"{analysed_s:.4g} s of analysable window against a required "
                    f"{required_s:.4g} s = {min_tau_windows:g} x the measured correlation time "
                    f"{tau.tau_int_s:.4g} s. That is n_eff = {tau.n_effective:.2f} independent "
                    f"samples, below the {float(min_tau_windows) / 2.0:.1f} the contract requires. "
                    "A transient cannot be told from a steady state here, so no mean is returned"
                ),
                "tau_estimate": tau.note,
                "naive_sem_would_have_been": sem_naive,
                "correlated_sem_would_have_been": sem_ess,
            },
        )

    mean = float(np.mean(tail))
    z_ess = (mean - float(null_value)) / sem_ess if sem_ess > 0.0 else None
    z_naive = (mean - float(null_value)) / sem_naive if sem_naive > 0.0 else None

    return StationarityReport(
        verdict=StationarityVerdict.STATIONARY,
        observable=observable,
        mean=mean,
        sem=sem_ess,
        sem_naive=sem_naive,
        sem_inflation=inflation,
        z_ess=z_ess,
        z_naive=z_naive,
        std=std,
        n_samples=n_kept,
        tau=tau,
        equilibration=search,
        analysed_time_s=analysed_s,
        drift_per_s=slope,
        drift_sigma_measured=drift_sigma_measured,
        contract=contract,
        notes={
            "tau_estimate": tau.note,
            "sem_basis": (
                f"sigma*sqrt(2*tau_int/T) = {sem_ess:.6g}. The naive sigma/sqrt(N) would be "
                f"{sem_naive:.6g}, i.e. {inflation:.3g}x too small on this series"
            ),
            "z_comparison": (
                f"z_ess = {z_ess:.4g} against z_naive = {z_naive:.4g}"
                if z_ess is not None and z_naive is not None
                else "zero variance; no significance to report"
            ),
            "tau_reliability": (
                f"the window spans {tau.reliability:.1f} correlation times"
                + (
                    ". Below about 50 the tau_int estimator is biased low, so n_eff is an "
                    "over-estimate and this sem should be read as a LOWER bound"
                    if tau.reliability < 50.0
                    else ". Comfortably inside the regime where the tau_int estimate is reliable"
                )
            ),
            "fdt": (
                "not assumed — in a driven steady state the fluctuation size is not set by "
                "temperature, so sigma is measured rather than predicted"
            ),
        },
    )
