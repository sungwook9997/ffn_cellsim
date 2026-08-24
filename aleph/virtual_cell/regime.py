"""What kind of dynamical state is an accepted trajectory actually in — and may we say so at all?

WHY THIS MODULE EXISTS.  Every gate in this repository has scored a trajectory by asking whether a
residual fell toward zero.  That question has a hidden premise: that the object being observed relaxes
to rest.  A cell does not.  Cortical tension is maintained by myosin turnover, blebbing is pulsatile,
and a spreading edge can be periodic or aperiodic — so **"the force converged to zero" is not the
success criterion, and a non-zero residual is not a failure.**  A stochastic non-equilibrium steady
state, a noise-sustained quasicycle and a deterministic limit cycle are all perfectly valid accepted
states; nothing in this module treats any of them as an error, and there is no code path here in which
a non-zero fluctuation, residual or drift is scored as a defect.

What *is* an error is quoting a number that the trajectory's regime does not support: a time-average
over a limit cycle, a finite difference across a bifurcation, an error bar computed from the step count
of a correlated series.  The regime is the precondition for those questions, so it is measured first
and separately.

THE THREE THINGS THIS DOES THAT A CONVERGENCE CHECK DOES NOT.

1. **It separates acceptance from classification.**  Whether a step was admissible is
   :mod:`aleph.virtual_cell.admissibility`'s question; whether the resulting trajectory is a fixed
   point, a NESS, a quasicycle, a limit cycle, quasiperiodic or chaotic is this one's.  Neither answer
   implies the other.
2. **It is allowed to refuse.**  :attr:`Regime.UNDECIDED` is a first-class, frequently-returned answer,
   and the module is deliberately built so that it is returned often: a series too short relative to its
   own measured correlation time, one still drifting, one whose spectral peak falls between the declared
   sharp and broad bands, and one with no temporal structure at all each return ``UNDECIDED`` with a
   named reason and the expansion the refusal opens.  A classifier that always commits would launder a
   guess into a label, and the label would then be used to decide whether a sensitivity may be reported.
3. **It owns no thresholds.**  Every evidence threshold lives in a caller-supplied
   :class:`RegimeEvidenceContract` with **no defaults**, declared before the run and copied into the
   report, so a later reader judges the verdict against the criterion that produced it.  A threshold
   invented inside this module would be a threshold chosen after seeing the data.

AUTHORITY.  This is an **unverified observer**.  It grants no evidence rung, and
:attr:`RegimeReport.evidence_authority` says so in every report it emits.  It imports no Warp, touches
no device, and runs entirely on CPU with numpy alone.

THE DISCRIMINATING STATISTICS, and why each is here rather than a heavier nonlinear-dynamics library:

* **Autocorrelation structure** — the normalised ACF and the Sokal-windowed integrated correlation time
  ``tau_int``.  This is the same estimator as
  :func:`aleph.engine.observe.stationarity.integrated_autocorrelation_time`, re-expressed here
  because (a) that module cannot be imported without initialising Warp, which this CPU-only layer
  forbids, and (b) the classifier needs the ACF *array* — its zero crossings and rebound lobes are what
  separate a decaying relaxation from an oscillation — which the engine helper does not return.  The
  test suite asserts the two ``tau_int`` estimates agree to machine precision, which is a stronger
  statement than sharing the call would have been.
* **Spectral peak sharpness and count** — a Welch periodogram, peak prominence against a rolling-median
  baseline, and peak width **measured in resolution bins** rather than in Hz.  The bin normalisation is
  the point: a deterministic cycle's linewidth is resolution-limited, so ``width_bins ~ 1`` is a
  grid-invariant statement about the signal while a width in Hz is a statement about the window length.
* **Phase-space recurrence** — delay embedding plus the RQA determinism ratio.  The fully developed
  logistic map has a *flat* power spectrum and a delta-function ACF: spectrally it is indistinguishable
  from white noise.  Only recurrence sees that it is deterministic, which is exactly why the module
  cannot be built on the spectrum alone.
* **Surrogate-referenced nonlinear predictability** — a nearest-neighbour analogue predictor scored
  against phase-randomised surrogates of the same series.  This one is here because the two obvious
  determinism statistics were measured on the controls and both failed: the raw Rosenstein Lyapunov
  exponent puts the Ornstein-Uhlenbeck control at 0.68/s against the logistic map's 0.64/s, and raw
  recurrence determinism puts OU at 0.78 against Lorenz's 0.99 while a phase surrogate *of* Lorenz
  scores 0.988.  Referencing against a spectrum-preserving null removes both failures at once, and is
  the statistic the classifier actually gates on.
* **Largest Lyapunov exponent** — Rosenstein's mean-log-divergence slope on the same embedding.
  **Reported, never thresholded**, for the reason just given; it is kept because it is the defining
  quantity of chaos and a reader is entitled to see it, not because it decides anything here.
* **Effective sample size** — ``n_eff = T / (2 tau_int)``.  Not a discriminator but a gate: below the
  declared ``min_n_eff`` no label is issued at all, because every statistic above is then noise.

units: time in SECONDS of physical time, frequency in Hz, never solver iterations.

Sanity Gate:
    * dimensional: ``dt`` is [s] and has no default; ``tau_int`` and every returned time are [s];
      frequencies are [Hz]; ``lyapunov_per_s`` is [1/s]; ``n_eff``, ``width_bins``, ``determinism`` and
      ``spectral_flatness`` are dimensionless.  Passing an iteration index as ``dt`` yields a
      correlation time in iterations and a Lyapunov exponent per iteration — which is why ``dt`` is
      required, and why the value is recorded in the report.
    * boundary: a series shorter than ``min_samples``, one whose analysed span covers fewer than
      ``min_correlation_windows`` of its own measured ``tau_int``, or one with ``n_eff`` below
      ``min_n_eff`` returns ``UNDECIDED`` and never a label.  A constant series has zero variance;
      ``tau_int`` is floored at ``dt/2`` rather than dividing by zero, and the constant series is a
      fixed point by the amplitude test, not by a spectral accident.
    * conservation/invariant: ``n_eff <= n_samples`` always, with equality only for an uncorrelated
      series.  ``0 <= determinism <= 1``, ``0 <= recurrence_rate <= 1``, ``0 <= spectral_flatness <= 1``,
      and every reported peak's power exceeds both its own continuum baseline and the declared share of
      the spectrum's maximum.
    * numerical: the ACF is computed by FFT on the mean-subtracted series (not subtracting the mean
      measures the mean, not the fluctuation), cost ``N log N``.  The Sokal partial sums are cumulative,
      so the window search is ``O(N)`` rather than the ``O(N^2)`` of the naive re-sum.  Recurrence and
      Lyapunov are ``O(M^2)`` in the embedded point count, which is why ``recurrence_max_points`` is a
      declared contract field and the subsample is strided, not random.
    * sign-sense: drift is signed and reported signed.  A positive Lyapunov exponent means divergence
      and is reported as such; the estimator's own noise floor on a non-chaotic control is not zero, so
      ``lyapunov_min_per_s`` is a declared threshold and not a test against zero.
    * measurement-protocol: every threshold is in the caller's :class:`RegimeEvidenceContract`, declared
      before the run and copied verbatim into :attr:`RegimeReport.contract`.  The contract forbids
      ``quasicycle_min_width_bins <= limit_cycle_max_width_bins``, which guarantees a non-empty
      ``UNDECIDED`` band between "sharp" and "broad" rather than a partition that must commit.
    * resolution limit: the classifier cannot resolve a linewidth below the Welch resolution
      ``1 / (segment length)``.  A quasicycle observed for fewer cycles than its own quality factor is
      therefore reported ``LIMIT_CYCLE`` or ``UNDECIDED``; the resolution is recorded in every report so
      that this specific failure is auditable after the fact rather than invisible.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from fractions import Fraction
from typing import Any

import numpy as np

__all__ = [
    "Regime",
    "RegimeEvidenceContract",
    "RegimeReport",
    "RegimeStatistics",
    "SpectralPeak",
    "classify_regime",
    "delay_embedding",
    "embedding_delay",
    "equilibration_point",
    "integrated_autocorrelation_time",
    "largest_lyapunov_exponent",
    "nonlinear_predictability",
    "normalized_autocorrelation",
    "phase_randomized_surrogate",
    "recurrence_determinism",
    "spectral_flatness",
    "spectral_peaks",
    "welch_psd",
]

#: What this layer is allowed to claim.  Mirrors the ``unverified-observer`` receipt vocabulary; a
#: regime label is never native evidence and never promotes an artifact.
EVIDENCE_AUTHORITY = "unverified-observer"


class Regime(StrEnum):
    """The dynamical regime an accepted trajectory is in.

    None of these is a failure.  ``FIXED_POINT`` is not "success" and ``CHAOS`` is not "divergence";
    they are different states a cell may legitimately be accepted in, and the label exists to decide
    which downstream questions are answerable, not whether the run was good.

    * ``FIXED_POINT`` — fluctuation negligible against the caller's declared amplitude scale.
    * ``STOCHASTIC_NESS`` — statistically stationary with a decaying autocorrelation and no significant
      spectral peak.  Sustained fluctuation about a mean; a time-average is meaningful here.
    * ``QUASICYCLE`` — a noise-sustained oscillation: one significant but BROAD spectral peak.  The
      deterministic part is a decaying spiral; the noise keeps re-exciting it.
    * ``LIMIT_CYCLE`` — one significant, resolution-limited-sharp spectral peak.  A time-average over
      this is not the observable; see :func:`aleph.virtual_cell.sweep.phase_decomposition`.
    * ``QUASIPERIODIC`` — two or more incommensurate fundamentals (a ratio not close to a rational with
      a small denominator, after harmonics of the strongest peak are removed).
    * ``CHAOS`` — deterministic by recurrence, aperiodic by spectrum, with a positive largest Lyapunov
      exponent above the declared floor.
    * ``UNDECIDED`` — the evidence does not support any of the above under the declared contract.  This
      is a real answer and is expected to be common.
    """

    FIXED_POINT = "fixed-point"
    STOCHASTIC_NESS = "stochastic-ness"
    QUASICYCLE = "quasicycle"
    LIMIT_CYCLE = "limit-cycle"
    QUASIPERIODIC = "quasiperiodic"
    CHAOS = "chaos"
    UNDECIDED = "undecided"


def _positive(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{what} must be positive and finite; got {value!r}")
    return number


def _nonnegative(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{what} must be nonnegative and finite; got {value!r}")
    return number


def _positive_int(value: object, *, what: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{what} must be an int")
    if value < minimum:
        raise ValueError(f"{what} must be >= {minimum}; got {value!r}")
    return value


def _unit_interval(value: object, *, what: str) -> float:
    number = _nonnegative(value, what=what)
    if number > 1.0:
        raise ValueError(f"{what} must lie in [0, 1]; got {value!r}")
    return number


@dataclass(frozen=True, slots=True)
class RegimeEvidenceContract:
    """The caller's pre-declared evidence thresholds.  **No field has a default, by design.**

    Every number the classifier compares against lives here.  The module invents none of them: a
    threshold chosen inside the classifier is a threshold chosen after seeing the data, and this project
    has retracted enough numbers that were read off a plot.  The contract is copied verbatim into every
    report so the verdict can be re-judged without the raw series.

    Args:
        min_samples: Shortest series that may be labelled at all.
        sokal_window_c: Sokal automatic-window constant; the window closes at the smallest ``M`` with
            ``M >= sokal_window_c * tau_int(M)``.  6 is the usual MD choice; shorter series need less.
        min_correlation_windows: The analysed span must cover at least this many measured ``tau_int``.
        min_n_eff: Minimum effective independent-sample count before any label is issued.
        drift_sigma: A series whose least-squares drift across the analysed window exceeds this many of
            its own standard deviations is a transient, and returns ``UNDECIDED`` rather than a regime.
        fixed_point_relative_amplitude: ``std / amplitude_scale`` below this is a fixed point.
        welch_segment_count: Number of 50 %-overlapping Hann segments in the periodogram.  This sets the
            spectral resolution and therefore what linewidth is resolvable at all.
        peak_prominence_ratio_min: A spectral maximum counts as a peak only if its power exceeds this
            multiple of the local rolling-median baseline.
        baseline_octave_ratio: Half-width of the peak-prominence baseline window as a MULTIPLICATIVE
            factor in frequency: bin ``f`` is compared to the median over ``[f / r, f * r]``.  A
            log-frequency window rather than a bin count, for the reason in :func:`_continuum_baseline`.
        peak_power_fraction_min: Minimum share of the LARGEST PSD value a candidate peak must carry.
        secondary_peak_power_fraction: A second fundamental is considered real only if its power is at
            least this fraction of the strongest peak's.
        harmonic_relative_tolerance: Relative frequency tolerance for calling a peak a harmonic of the
            fundamental, or a frequency ratio commensurate.  Floored at the spectral resolution.
        max_rational_denominator: Largest denominator considered "a small rational".  A frequency ratio
            approximable to within ``harmonic_relative_tolerance`` by ``p/q`` with ``q`` at most this is
            commensurate — one cycle with harmonics, not two independent frequencies.
        limit_cycle_max_width_bins: Peak half-power width at or below this many resolution bins is
            sharp — a limit cycle.
        quasicycle_min_width_bins: Peak width at or above this many bins is broad — a quasicycle.  Must
            be strictly greater than ``limit_cycle_max_width_bins``; the gap between them is the
            ``UNDECIDED`` band, and a contract without such a gap is rejected.
        white_noise_tau_dt_max: ``tau_int / dt`` at or below this is an uncorrelated sequence.
        white_noise_flatness_min: Spectral flatness at or above this is a broadband, featureless
            spectrum.
        determinism_min: Minimum RQA determinism (the fraction of recurrence points lying on diagonal
            lines of at least ``min_diagonal_length``) required before a trajectory may be called
            deterministic.  A necessary condition only — the OU control reaches 0.78 unaided, so this
            threshold is a sanity floor and the surrogate-referenced test below does the separating.
        predictability_z_min: Minimum surrogate-referenced nonlinear-predictability z-score.  This is
            the discriminating statistic for determinism; see :func:`nonlinear_predictability` for the
            measured reason the Lyapunov exponent is reported but not gated on.
        predictability_neighbours: Analogues averaged into each nearest-neighbour prediction.
        n_surrogates: Size of the phase-randomised surrogate ensemble.
        surrogate_seed: RNG seed for that ensemble, so a verdict is reproducible.
        embedding_dimension: Delay-embedding dimension.
        recurrence_rate_target: The recurrence threshold is chosen as the distance quantile that gives
            this recurrence rate, so the statistic does not depend on the observable's units.
        min_diagonal_length: Shortest diagonal line counted as deterministic.
        recurrence_max_points: Cap on embedded points entering the ``O(M^2)`` recurrence and Lyapunov
            computations; the subsample is strided so it spans the trajectory.
        theiler_tau_multiple: Temporal-neighbour exclusion, in multiples of the measured ``tau_int``,
            applied to BOTH the recurrence diagonals and the Lyapunov neighbour search.  Without it a
            point's nearest neighbour is its own predecessor and both statistics measure the sampling
            rate rather than the dynamics.
        lyapunov_fit_tau_multiple: Length of the mean-log-divergence fit region, in multiples of the
            measured ``tau_int``.  Declaring it against the series' own memory time rather than as an
            absolute lag count is what lets one contract cover a map with ``tau_int = dt/2`` and a flow
            with ``tau_int = 15 dt`` without per-control tuning.
        lyapunov_min_fit_lags: Floor on that resolved lag count.
        lyapunov_max_fit_lags: Ceiling on that resolved lag count.
    """

    min_samples: int
    sokal_window_c: float
    min_correlation_windows: float
    min_n_eff: float
    drift_sigma: float
    fixed_point_relative_amplitude: float
    welch_segment_count: int
    peak_prominence_ratio_min: float
    baseline_octave_ratio: float
    peak_power_fraction_min: float
    secondary_peak_power_fraction: float
    harmonic_relative_tolerance: float
    max_rational_denominator: int
    limit_cycle_max_width_bins: float
    quasicycle_min_width_bins: float
    white_noise_tau_dt_max: float
    white_noise_flatness_min: float
    determinism_min: float
    predictability_z_min: float
    predictability_neighbours: int
    n_surrogates: int
    surrogate_seed: int
    embedding_dimension: int
    recurrence_rate_target: float
    min_diagonal_length: int
    recurrence_max_points: int
    theiler_tau_multiple: float
    lyapunov_fit_tau_multiple: float
    lyapunov_min_fit_lags: int
    lyapunov_max_fit_lags: int

    def __post_init__(self) -> None:
        _positive_int(self.min_samples, what="min_samples", minimum=8)
        _positive(self.sokal_window_c, what="sokal_window_c")
        _positive(self.min_correlation_windows, what="min_correlation_windows")
        _positive(self.min_n_eff, what="min_n_eff")
        _positive(self.drift_sigma, what="drift_sigma")
        _positive(self.fixed_point_relative_amplitude, what="fixed_point_relative_amplitude")
        _positive_int(self.welch_segment_count, what="welch_segment_count")
        if self.peak_prominence_ratio_min <= 1.0:
            raise ValueError(
                "peak_prominence_ratio_min must exceed 1 (a peak must beat its baseline)"
            )
        if self.baseline_octave_ratio <= 1.0:
            raise ValueError(
                "baseline_octave_ratio must exceed 1 (a window of zero width is no window)"
            )
        _unit_interval(self.peak_power_fraction_min, what="peak_power_fraction_min")
        _unit_interval(self.secondary_peak_power_fraction, what="secondary_peak_power_fraction")
        _positive(self.harmonic_relative_tolerance, what="harmonic_relative_tolerance")
        _positive_int(self.max_rational_denominator, what="max_rational_denominator", minimum=2)
        _positive(self.limit_cycle_max_width_bins, what="limit_cycle_max_width_bins")
        _positive(self.quasicycle_min_width_bins, what="quasicycle_min_width_bins")
        if self.quasicycle_min_width_bins <= self.limit_cycle_max_width_bins:
            raise ValueError(
                "quasicycle_min_width_bins must exceed limit_cycle_max_width_bins; the gap between "
                "them is the UNDECIDED band, and a contract with no gap forces the classifier to "
                "commit on evidence it does not have"
            )
        _positive(self.white_noise_tau_dt_max, what="white_noise_tau_dt_max")
        _unit_interval(self.white_noise_flatness_min, what="white_noise_flatness_min")
        _unit_interval(self.determinism_min, what="determinism_min")
        _positive(self.predictability_z_min, what="predictability_z_min")
        _positive_int(self.predictability_neighbours, what="predictability_neighbours")
        _positive_int(self.n_surrogates, what="n_surrogates", minimum=2)
        if isinstance(self.surrogate_seed, bool) or not isinstance(self.surrogate_seed, int):
            raise TypeError("surrogate_seed must be an int")
        if self.surrogate_seed < 0:
            raise ValueError("surrogate_seed must be nonnegative")
        _positive_int(self.embedding_dimension, what="embedding_dimension", minimum=2)
        rate = _positive(self.recurrence_rate_target, what="recurrence_rate_target")
        if rate >= 1.0:
            raise ValueError("recurrence_rate_target must lie in (0, 1)")
        _positive_int(self.min_diagonal_length, what="min_diagonal_length", minimum=2)
        _positive_int(self.recurrence_max_points, what="recurrence_max_points", minimum=32)
        _positive(self.theiler_tau_multiple, what="theiler_tau_multiple")
        _positive(self.lyapunov_fit_tau_multiple, what="lyapunov_fit_tau_multiple")
        _positive_int(self.lyapunov_min_fit_lags, what="lyapunov_min_fit_lags", minimum=3)
        _positive_int(self.lyapunov_max_fit_lags, what="lyapunov_max_fit_lags", minimum=3)
        if self.lyapunov_min_fit_lags > self.lyapunov_max_fit_lags:
            raise ValueError("lyapunov_min_fit_lags must not exceed lyapunov_max_fit_lags")

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able contract, recorded next to every verdict it produced."""
        return {
            "min_samples": int(self.min_samples),
            "sokal_window_c": float(self.sokal_window_c),
            "min_correlation_windows": float(self.min_correlation_windows),
            "min_n_eff": float(self.min_n_eff),
            "drift_sigma": float(self.drift_sigma),
            "fixed_point_relative_amplitude": float(self.fixed_point_relative_amplitude),
            "welch_segment_count": int(self.welch_segment_count),
            "peak_prominence_ratio_min": float(self.peak_prominence_ratio_min),
            "baseline_octave_ratio": float(self.baseline_octave_ratio),
            "peak_power_fraction_min": float(self.peak_power_fraction_min),
            "secondary_peak_power_fraction": float(self.secondary_peak_power_fraction),
            "harmonic_relative_tolerance": float(self.harmonic_relative_tolerance),
            "max_rational_denominator": int(self.max_rational_denominator),
            "limit_cycle_max_width_bins": float(self.limit_cycle_max_width_bins),
            "quasicycle_min_width_bins": float(self.quasicycle_min_width_bins),
            "white_noise_tau_dt_max": float(self.white_noise_tau_dt_max),
            "white_noise_flatness_min": float(self.white_noise_flatness_min),
            "determinism_min": float(self.determinism_min),
            "predictability_z_min": float(self.predictability_z_min),
            "predictability_neighbours": int(self.predictability_neighbours),
            "n_surrogates": int(self.n_surrogates),
            "surrogate_seed": int(self.surrogate_seed),
            "embedding_dimension": int(self.embedding_dimension),
            "recurrence_rate_target": float(self.recurrence_rate_target),
            "min_diagonal_length": int(self.min_diagonal_length),
            "recurrence_max_points": int(self.recurrence_max_points),
            "theiler_tau_multiple": float(self.theiler_tau_multiple),
            "lyapunov_fit_tau_multiple": float(self.lyapunov_fit_tau_multiple),
            "lyapunov_min_fit_lags": int(self.lyapunov_min_fit_lags),
            "lyapunov_max_fit_lags": int(self.lyapunov_max_fit_lags),
        }


@dataclass(frozen=True, slots=True)
class SpectralPeak:
    """One resolved spectral maximum, with its width expressed in resolution bins.

    ``width_bins`` rather than ``width_hz`` is what the classifier thresholds on: a deterministic
    cycle's linewidth is limited by the observation window, so a width in Hz says how long the window
    was, while a width in bins says how sharp the signal is.
    """

    frequency_hz: float
    power: float
    prominence_ratio: float
    width_hz: float
    width_bins: float
    quality_factor: float
    width_truncated: bool

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able peak record."""
        return {
            "frequency_hz": self.frequency_hz,
            "power": self.power,
            "prominence_ratio": self.prominence_ratio,
            "width_hz": self.width_hz,
            "width_bins": self.width_bins,
            "quality_factor": self.quality_factor,
            "width_truncated": self.width_truncated,
        }


@dataclass(frozen=True, slots=True)
class RegimeStatistics:
    """Everything measured on the analysed window, independent of what label it produced.

    The statistics are reported whatever the verdict — including for ``UNDECIDED`` — because the point
    of a refusal is to say which statistic was missing, and a refusal with no numbers cannot.
    """

    n_samples: int
    dt_s: float
    analysed_time_s: float
    equilibration_time_s: float
    mean: float
    std: float
    relative_amplitude: float
    drift_per_s: float
    drift_over_std: float
    tau_int_s: float
    n_eff: float
    tau_note: str
    acf_first_zero_crossing_s: float | None
    acf_rebound: float
    spectral_resolution_hz: float
    spectral_flatness: float
    peaks: tuple[SpectralPeak, ...]
    fundamental_frequencies_hz: tuple[float, ...]
    frequency_ratio: float | None
    ratio_is_commensurate: bool | None
    embedding_delay_samples: int
    embedding_dimension: int
    embedded_points: int
    embedded_points_analysed: int
    recurrence_rate: float
    determinism: float
    max_diagonal_length: int
    lyapunov_per_s: float
    lyapunov_fit_r2: float
    lyapunov_fit_lags_used: int
    prediction_horizon_samples: int
    prediction_error: float
    surrogate_prediction_error_mean: float
    surrogate_prediction_error_std: float
    predictability_z: float

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able statistics block."""
        return {
            "n_samples": self.n_samples,
            "dt_s": self.dt_s,
            "analysed_time_s": self.analysed_time_s,
            "equilibration_time_s": self.equilibration_time_s,
            "mean": self.mean,
            "std": self.std,
            "relative_amplitude": self.relative_amplitude,
            "drift_per_s": self.drift_per_s,
            "drift_over_std": self.drift_over_std,
            "tau_int_s": self.tau_int_s,
            "n_eff": self.n_eff,
            "tau_note": self.tau_note,
            "acf_first_zero_crossing_s": self.acf_first_zero_crossing_s,
            "acf_rebound": self.acf_rebound,
            "spectral_resolution_hz": self.spectral_resolution_hz,
            "spectral_flatness": self.spectral_flatness,
            "peaks": [peak.as_dict() for peak in self.peaks],
            "fundamental_frequencies_hz": list(self.fundamental_frequencies_hz),
            "frequency_ratio": self.frequency_ratio,
            "ratio_is_commensurate": self.ratio_is_commensurate,
            "embedding_delay_samples": self.embedding_delay_samples,
            "embedding_dimension": self.embedding_dimension,
            "embedded_points": self.embedded_points,
            "embedded_points_analysed": self.embedded_points_analysed,
            "recurrence_rate": self.recurrence_rate,
            "determinism": self.determinism,
            "max_diagonal_length": self.max_diagonal_length,
            "lyapunov_per_s": self.lyapunov_per_s,
            "lyapunov_fit_r2": self.lyapunov_fit_r2,
            "lyapunov_fit_lags_used": self.lyapunov_fit_lags_used,
            "prediction_horizon_samples": self.prediction_horizon_samples,
            "prediction_error": self.prediction_error,
            "surrogate_prediction_error_mean": self.surrogate_prediction_error_mean,
            "surrogate_prediction_error_std": self.surrogate_prediction_error_std,
            "predictability_z": self.predictability_z,
            "lyapunov_is_not_gated": (
                "reported, never thresholded: on the analytic controls the OU process scores 0.68/s "
                "against the logistic map's 0.64/s, so the raw exponent does not separate them"
            ),
        }


@dataclass(frozen=True, slots=True)
class RegimeReport:
    """One observable's regime verdict, its evidence, and — when refused — where the refusal routes.

    ``evidence_authority`` is ``unverified-observer`` in every report this module emits.  A regime label
    is a statement about a time series, never about physical validity, and it may not be promoted to
    native evidence by any caller.
    """

    regime: Regime
    observable: str
    reason: str
    statistics: RegimeStatistics
    contract: dict[str, Any]
    expansions: tuple[str, ...] = ()
    evidence_authority: str = EVIDENCE_AUTHORITY

    @property
    def is_decided(self) -> bool:
        """Whether a regime label was actually issued."""
        return self.regime is not Regime.UNDECIDED

    @property
    def is_oscillatory(self) -> bool:
        """Whether a time-average is an incomplete description of this trajectory.

        ``True`` for a quasicycle, a limit cycle and a quasiperiodic trajectory.  A caller reporting a
        scalar mean for one of these is reporting the centre of an orbit as if it were the orbit; see
        :func:`aleph.virtual_cell.sweep.phase_decomposition`.
        """
        return self.regime in {Regime.QUASICYCLE, Regime.LIMIT_CYCLE, Regime.QUASIPERIODIC}

    @property
    def dominant_period_s(self) -> float | None:
        """The strongest fundamental's period [s], or ``None`` if there is no resolved fundamental."""
        if not self.statistics.fundamental_frequencies_hz:
            return None
        frequency = self.statistics.fundamental_frequencies_hz[0]
        if frequency <= 0.0:
            return None
        return 1.0 / frequency

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able report."""
        return {
            "regime": self.regime.value,
            "observable": self.observable,
            "reason": self.reason,
            "statistics": self.statistics.as_dict(),
            "contract": dict(self.contract),
            "expansions": list(self.expansions),
            "evidence_authority": self.evidence_authority,
        }


def _as_series(values: Any) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64).ravel()
    if array.size == 0:
        raise ValueError("series must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError("series must contain only finite values")
    return array


def normalized_autocorrelation(series: Any) -> np.ndarray:
    """Return the normalised autocorrelation ``rho[k]`` of a series, with ``rho[0] == 1``.

    The mean is subtracted before correlating; without that the result measures the mean rather than
    the fluctuation.  Computed by FFT, so the cost is ``N log N``.

    Args:
        series: Observable samples, evenly spaced in time.

    Returns:
        A length-``N`` array of normalised autocorrelations.  A zero-variance series returns
        ``[1, 0, 0, ...]``, which is the honest statement that no fluctuation was observed.

    Raises:
        ValueError: If the series is empty, non-finite, or shorter than 2 samples.
    """
    x = _as_series(series)
    if x.size < 2:
        raise ValueError("need >= 2 samples for an autocorrelation")
    y = x - x.mean()
    variance = float(np.dot(y, y))
    if variance <= 0.0:
        acf = np.zeros(x.size, dtype=np.float64)
        acf[0] = 1.0
        return acf
    size = 1 << (2 * x.size - 1).bit_length()
    spectrum = np.fft.rfft(y, size)
    acf = np.fft.irfft(spectrum * np.conjugate(spectrum), size)[: x.size].real
    return acf / acf[0]


def integrated_autocorrelation_time(
    series: Any, dt: float, *, sokal_window_c: float
) -> tuple[float, float, str]:
    """Return ``(tau_int_s, n_eff, note)`` using the Sokal automatic window.

    ``tau_int = dt * (1/2 + sum_{k=1..M} rho_k)`` with ``M`` the smallest lag satisfying
    ``M >= sokal_window_c * tau_int(M)``.  The window matters: each lag of the ACF tail carries noise of
    order ``1/sqrt(N)``, so summing every lag adds a random walk to the answer and typically inflates
    ``tau`` past the signal.

    This is numerically the same estimator as
    :func:`aleph.engine.observe.stationarity.integrated_autocorrelation_time`; the test suite
    asserts the two agree.  It is restated here only because that module cannot be imported without
    initialising Warp, which this CPU-only layer forbids.

    Args:
        series: Observable samples, evenly spaced.
        dt: Physical time between samples [s].  Required and unitful: an iteration index here silently
            returns a correlation time in iterations.
        sokal_window_c: The window constant, from the caller's contract.

    Returns:
        ``tau_int`` [s] floored at ``dt/2``, the independent-sample count ``n_eff = T / (2 tau_int)``,
        and a note naming what happened when the estimate was degenerate.

    Raises:
        ValueError: If ``dt`` is not positive-finite or the series has fewer than 4 samples.
    """
    dt_s = _positive(dt, what="dt")
    _positive(sokal_window_c, what="sokal_window_c")
    x = _as_series(series)
    if x.size < 4:
        raise ValueError(f"need >= 4 samples to estimate a correlation time; got {x.size}")

    n = x.size
    total_time = n * dt_s
    if float(np.var(x)) <= 0.0:
        return dt_s * 0.5, float(n), "series has zero variance; tau_int floored at dt/2"

    acf = normalized_autocorrelation(x)
    # Cumulative partial sums make the window search O(N) instead of the naive O(N^2) re-sum.
    partial = 0.5 + np.cumsum(acf[1:])
    lags = np.arange(1, n, dtype=np.float64)
    closed = np.flatnonzero(lags >= sokal_window_c * partial)
    if closed.size == 0:
        tau = max(float(partial[-1]), 0.5)
        note = (
            f"Sokal window never closed within {n} samples: the series is short relative to its own "
            "correlation time, so tau_int is a LOWER bound and n_eff an upper one"
        )
        return tau * dt_s, max(1.0, total_time / (2.0 * tau * dt_s)), note

    window = int(lags[closed[0]])
    tau = max(float(partial[closed[0]]), 0.5)
    tau_s = tau * dt_s
    return (
        tau_s,
        max(1.0, total_time / (2.0 * tau_s)),
        f"Sokal window closed at lag {window} (c={sokal_window_c:g})",
    )


def equilibration_point(series: Any, dt: float, *, sokal_window_c: float, n_candidates: int) -> int:
    """Return the discard index whose remaining series carries the most independent samples.

    Chodera's rule.  Discarding too little leaves the transient in and inflates the correlation time;
    discarding too much throws away signal.  The maximum trades those off by a rule rather than by eye.

    Args:
        series: Observable samples.
        dt: Physical time between samples [s].
        sokal_window_c: The caller's Sokal window constant.
        n_candidates: How many truncation points to try across the first half of the series.

    Returns:
        The index to start the analysed window from.
    """
    x = _as_series(series)
    best_index, best_n_eff = 0, -1.0
    limit = max(1, x.size // 2)
    step = max(1, limit // max(1, n_candidates))
    for start in range(0, limit, step):
        tail = x[start:]
        if tail.size < 4:
            break
        try:
            _, n_eff, _ = integrated_autocorrelation_time(tail, dt, sokal_window_c=sokal_window_c)
        except ValueError:
            continue
        if n_eff > best_n_eff:
            best_index, best_n_eff = start, n_eff
    return best_index


def welch_psd(series: Any, dt: float, *, n_segments: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Return ``(frequencies_hz, psd, resolution_hz)`` — a Hann-windowed, 50 %-overlap Welch estimate.

    Welch rather than a single periodogram because a raw periodogram's variance does not fall with the
    record length: every bin has a 100 % standard error however long the run, which makes peak
    prominence unusable.  Segment averaging trades resolution for that variance, and the traded-away
    resolution is returned so the caller can see what linewidth is resolvable at all.

    Args:
        series: Observable samples, evenly spaced.
        dt: Physical time between samples [s].
        n_segments: Number of overlapping segments to average.

    Returns:
        One-sided frequencies [Hz] excluding DC, the corresponding power spectral density, and the
        resolution ``1 / (segment length in seconds)`` [Hz].

    Raises:
        ValueError: If the series is too short to form even one usable segment.
    """
    dt_s = _positive(dt, what="dt")
    segments = _positive_int(n_segments, what="n_segments")
    x = _as_series(series)
    n = x.size
    # 50 % overlap: k segments of length L cover L*(k+1)/2 samples.
    length = int(max(8, (2 * n) // (segments + 1)))
    length = min(length, n)
    if length < 8:
        raise ValueError(f"series of {n} samples cannot supply a usable Welch segment")
    step = max(1, length // 2)
    window = np.hanning(length)
    window_power = float(np.dot(window, window))

    starts = list(range(0, n - length + 1, step))
    accumulator = np.zeros(length // 2 + 1, dtype=np.float64)
    for start in starts:
        chunk = x[start : start + length]
        chunk = chunk - chunk.mean()
        spectrum = np.fft.rfft(chunk * window)
        accumulator += (np.abs(spectrum) ** 2) / (window_power / dt_s)
    psd = accumulator / float(len(starts))
    # One-sided: fold the negative half onto the positive, except DC and Nyquist.
    if psd.size > 2:
        psd[1:-1] *= 2.0
    frequencies = np.fft.rfftfreq(length, d=dt_s)
    resolution = 1.0 / (length * dt_s)
    return frequencies[1:], psd[1:], resolution


def spectral_flatness(psd: Any) -> float:
    """Return Wiener entropy: the geometric mean of the PSD over its arithmetic mean.

    1 for a perfectly flat (white) spectrum, near 0 for a spectrum dominated by one line.  Computed in
    log space so a very peaked spectrum does not underflow the geometric mean to zero.

    Args:
        psd: Nonnegative power spectral density values, DC excluded.

    Returns:
        A value in ``[0, 1]``.  A PSD containing a zero returns 0.0 — the geometric mean is genuinely
        zero there, and substituting a floor would be inventing a number.
    """
    power = np.asarray(psd, dtype=np.float64).ravel()
    if power.size == 0:
        raise ValueError("psd must not be empty")
    if np.any(power < 0.0):
        raise ValueError("psd must be nonnegative")
    arithmetic = float(np.mean(power))
    if arithmetic <= 0.0:
        return 0.0
    if np.any(power <= 0.0):
        return 0.0
    geometric = float(np.exp(np.mean(np.log(power))))
    return float(min(1.0, geometric / arithmetic))


def _continuum_baseline(
    frequencies: np.ndarray, psd: np.ndarray, octave_ratio: float
) -> np.ndarray:
    """Return the smooth continuum under a PSD: a rolling median over a LOG-frequency window.

    A fixed-width bin window cannot serve both jobs a baseline has to do.  It must be wide enough not
    to sit inside a broad line — otherwise a noise-sustained peak is compared to itself and reports a
    prominence near 1, which is how the quasicycle control went undetected at 25 bins — and narrow
    enough to track a sloped continuum, or the plateau of a Lorentzian is compared against its own
    high-frequency roll-off and reports a peak that is not there, which is how the Ornstein-Uhlenbeck
    control grew a spurious 0.056 Hz line at 101 bins.

    A window proportional to frequency does both, because the two things a continuum does here — a
    flat plateau and a power-law roll-off — are both slowly varying in log-frequency, while a spectral
    line is narrow on that scale wherever it sits.

    Args:
        frequencies: One-sided frequencies [Hz], DC excluded, ascending.
        psd: The corresponding power spectral density.
        octave_ratio: Half-width of the window as a multiplicative factor: bin ``f`` is compared to the
            median over ``[f / ratio, f * ratio]``.

    Returns:
        The baseline, same shape as ``psd``.
    """
    lower = np.searchsorted(frequencies, frequencies / octave_ratio, side="left")
    upper = np.searchsorted(frequencies, frequencies * octave_ratio, side="right")
    baseline = np.empty(psd.size, dtype=np.float64)
    for i in range(psd.size):
        start, stop = int(lower[i]), int(upper[i])
        if stop - start < 3:
            start, stop = max(0, i - 1), min(psd.size, i + 2)
        baseline[i] = float(np.median(psd[start:stop]))
    return baseline


def spectral_peaks(
    frequencies: np.ndarray,
    psd: np.ndarray,
    resolution_hz: float,
    *,
    prominence_ratio_min: float,
    baseline_octave_ratio: float,
    power_fraction_min: float,
) -> tuple[SpectralPeak, ...]:
    """Return the significant spectral maxima, strongest first.

    A local maximum counts only if its power exceeds ``prominence_ratio_min`` times the
    :func:`_continuum_baseline` beneath it: an absolute power threshold would depend on the
    observable's units, and comparing to the global mean would let one strong line suppress every
    other.

    The half-power width is found by walking down each flank to ``peak / 2`` with linear interpolation
    and is reported in resolution bins.  A flank that reaches the array edge without halving sets
    ``width_truncated``; the width is then a LOWER bound, which matters because a truncated width can
    only make a peak look sharper than it is.

    Args:
        frequencies: One-sided frequencies [Hz], DC excluded.
        psd: The corresponding power spectral density.
        resolution_hz: The spectral resolution [Hz], from :func:`welch_psd`.
        prominence_ratio_min: Minimum power-over-baseline ratio, from the caller's contract.
        baseline_octave_ratio: Log-frequency half-width of the baseline window, from the contract.
        power_fraction_min: A candidate must also carry at least this fraction of the LARGEST PSD
            value.  Prominence alone is a ratio and says nothing about magnitude: in the far tail of a
            smooth flow's spectrum, where the power has fallen nine orders of magnitude and what is
            left is arithmetic noise, ratios above any threshold occur freely.  The Lorenz control
            produced 191 such "peaks", every one of them below 2 % of the spectrum's maximum, and the
            trajectory was reported quasiperiodic on the strength of two of them.

    Returns:
        Peaks ordered by descending power.
    """
    power = np.asarray(psd, dtype=np.float64).ravel()
    freq = np.asarray(frequencies, dtype=np.float64).ravel()
    if power.shape != freq.shape:
        raise ValueError("frequencies and psd must have the same shape")
    if power.size < 3:
        return ()
    baseline = _continuum_baseline(freq, power, baseline_octave_ratio)
    power_floor = float(power_fraction_min) * float(np.max(power))
    interior = np.arange(1, power.size - 1)
    is_max = (power[interior] > power[interior - 1]) & (power[interior] >= power[interior + 1])
    candidates = interior[is_max]

    peaks: list[SpectralPeak] = []
    for index in candidates:
        if float(power[index]) < power_floor:
            continue
        base = float(baseline[index])
        if base <= 0.0:
            continue
        ratio = float(power[index]) / base
        if ratio < prominence_ratio_min:
            continue
        half = float(power[index]) / 2.0
        left, truncated = float(freq[0]), True
        for j in range(index - 1, -1, -1):
            if power[j] <= half:
                span = power[j + 1] - power[j]
                frac = 0.0 if span == 0.0 else (half - power[j]) / span
                left, truncated = float(freq[j] + frac * (freq[j + 1] - freq[j])), False
                break
        right, right_truncated = float(freq[-1]), True
        for j in range(index + 1, power.size):
            if power[j] <= half:
                span = power[j - 1] - power[j]
                frac = 0.0 if span == 0.0 else (half - power[j]) / span
                right, right_truncated = float(freq[j] - frac * (freq[j] - freq[j - 1])), False
                break
        width_hz = max(right - left, float(resolution_hz))
        peaks.append(
            SpectralPeak(
                frequency_hz=float(freq[index]),
                power=float(power[index]),
                prominence_ratio=ratio,
                width_hz=width_hz,
                width_bins=width_hz / float(resolution_hz),
                quality_factor=float(freq[index]) / width_hz,
                width_truncated=bool(truncated or right_truncated),
            )
        )
    peaks.sort(key=lambda peak: peak.power, reverse=True)
    return tuple(peaks)


def embedding_delay(acf: np.ndarray, *, max_delay: int) -> int:
    """Return the delay-embedding lag: the first lag at which the ACF has decayed to ``1/e``.

    The decorrelation lag, not the first zero crossing.  The zero-crossing rule is the more commonly
    quoted one and it fails on exactly the signals this module has to handle: a trajectory that dwells
    in one lobe of an attractor for several revolutions keeps a positive ACF far past its decorrelation
    time.  On the Lorenz control the zero crossing lands at lag 177 against a 16-sample ``1/e`` lag, and
    the resulting embedding is so over-stretched that nonlinear prediction on it performs *worse* than
    on a phase-randomised surrogate (normalised error 1.40 against 1.34) — i.e. the embedding, not the
    dynamics, destroyed the structure.

    Args:
        acf: Normalised autocorrelation with ``acf[0] == 1``.
        max_delay: Largest lag that may be returned.

    Returns:
        A lag of at least 1.  A series whose ACF never falls to ``1/e`` within ``max_delay`` returns
        ``max_delay``, which is reported rather than hidden.
    """
    limit = int(min(max_delay, acf.size - 1))
    if limit < 1:
        return 1
    decayed = np.flatnonzero(acf[1 : limit + 1] <= math.exp(-1.0))
    if decayed.size:
        return int(decayed[0] + 1)
    return max(1, limit)


def phase_randomized_surrogate(series: Any, rng: np.random.Generator) -> np.ndarray:
    """Return a surrogate with the SAME power spectrum and mean but randomised Fourier phases.

    This is the null hypothesis the determinism statistics have to beat: a linear Gaussian process with
    the observed spectrum.  Randomising the phases preserves the autocorrelation exactly — so anything
    a linear predictor could exploit survives — while destroying every nonlinear relation between
    Fourier components.  A statistic that scores the data no better than these surrogates is a statistic
    that has measured the spectrum, not the dynamics.

    Args:
        series: Observable samples.
        rng: The caller's generator, so the surrogate ensemble is reproducible.

    Returns:
        One surrogate realisation of the same length.
    """
    x = _as_series(series)
    spectrum = np.fft.rfft(x - x.mean())
    phases = rng.uniform(0.0, 2.0 * math.pi, spectrum.size)
    phases[0] = 0.0
    if x.size % 2 == 0:
        phases[-1] = 0.0
    return np.fft.irfft(np.abs(spectrum) * np.exp(1j * phases), x.size) + x.mean()


def _prediction_error(
    series: np.ndarray,
    *,
    delay: int,
    dimension: int,
    theiler: int,
    horizon: int,
    neighbours: int,
    max_points: int,
) -> float:
    """Return the nearest-neighbour analogue prediction error, normalised by the observable's std."""
    embedded = delay_embedding(series, delay=delay, dimension=dimension)
    points = _leading_window(embedded, max_points)
    count = points.shape[0]
    usable = count - horizon
    if usable < 32:
        return float("nan")
    diff = points[:, None, :] - points[None, :, :]
    distance = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
    index = np.arange(count)
    distance[np.abs(index[:, None] - index[None, :]) <= max(1, theiler)] = np.inf
    distance[:, usable:] = np.inf
    window = distance[:usable]
    # A row whose every candidate is excluded (a Theiler window wider than the analysed span) has no
    # analogue at all.  ``argsort`` is not stable, so such a row returns an ARBITRARY permutation of
    # column indices — which then indexes past the array.  Drop those rows rather than clipping them,
    # because a clipped index is a fabricated prediction.
    usable_rows = np.flatnonzero(np.count_nonzero(np.isfinite(window), axis=1) >= neighbours)
    if usable_rows.size < 32:
        return float("nan")
    order = np.argsort(window[usable_rows], axis=1)[:, :neighbours]
    target = points[:, -1]
    scale = float(np.std(target))
    if scale <= 0.0:
        return float("nan")
    predicted = target[order + horizon].mean(axis=1)
    truth = target[usable_rows + horizon]
    return float(np.sqrt(np.mean((predicted - truth) ** 2)) / scale)


def nonlinear_predictability(
    series: Any,
    *,
    delay: int,
    dimension: int,
    theiler: int,
    horizon: int,
    neighbours: int,
    max_points: int,
    n_surrogates: int,
    seed: int,
) -> tuple[float, float, float, float]:
    """Return ``(error, surrogate_mean, surrogate_std, z)`` for surrogate-referenced predictability.

    A nearest-neighbour analogue predictor is run on the delay embedding and on ``n_surrogates``
    phase-randomised surrogates of the same series.  ``z = (surrogate_mean - error) / surrogate_std``:
    positive means the data are predictable by a *nonlinear* rule beyond what their own spectrum
    explains, which is the operational content of "deterministic".

    **This, and not the Lyapunov exponent, is the statistic that separates chaos from a stochastic
    steady state on these controls.**  Measured on the analytic set, the raw Rosenstein exponent gives
    the Ornstein-Uhlenbeck control 0.68/s against the logistic map's 0.64/s — the stochastic control
    scores *higher* than the chaotic one, because for a noisy signal the mean log divergence is
    dominated by the noise decorrelating neighbours.  Raw determinism fails the same way once the signal
    is smooth: the OU control reaches 0.78 against Lorenz's 0.99, and a phase surrogate of Lorenz scores
    0.988, so determinism alone certifies the surrogate too.  Referencing a nonlinear statistic against
    a spectrum-preserving null is what removes both failures.

    Args:
        series: Observable samples.
        delay: Embedding lag in samples.
        dimension: Embedding dimension.
        theiler: Temporal-neighbour exclusion in samples.
        horizon: Prediction horizon in samples.  Set from the series' own correlation time by the
            caller: a map with ``tau_int = dt/2`` is predictable one step and a flow with
            ``tau_int = 23 dt`` is predictable for tens, and one fixed horizon cannot serve both.
        neighbours: How many analogues are averaged into each prediction.
        max_points: Cap on embedded points; a contiguous leading window.
        n_surrogates: Size of the surrogate ensemble.
        seed: Surrogate RNG seed, so the ensemble and therefore the verdict are reproducible.

    Returns:
        The data error, the surrogate ensemble mean and standard deviation, and the z-score.  A
        degenerate case (no usable prediction, or a zero-variance surrogate ensemble) returns
        ``z = 0.0``, which is a refusal to claim determinism rather than a claim of its absence.
    """
    x = _as_series(series)
    error = _prediction_error(
        x,
        delay=delay,
        dimension=dimension,
        theiler=theiler,
        horizon=horizon,
        neighbours=neighbours,
        max_points=max_points,
    )
    if not math.isfinite(error):
        return float("nan"), float("nan"), float("nan"), 0.0
    rng = np.random.default_rng(seed)
    scores = []
    for _ in range(max(2, int(n_surrogates))):
        surrogate = phase_randomized_surrogate(x, rng)
        value = _prediction_error(
            surrogate,
            delay=delay,
            dimension=dimension,
            theiler=theiler,
            horizon=horizon,
            neighbours=neighbours,
            max_points=max_points,
        )
        if math.isfinite(value):
            scores.append(value)
    if len(scores) < 2:
        return error, float("nan"), float("nan"), 0.0
    ensemble = np.asarray(scores, dtype=np.float64)
    mean = float(ensemble.mean())
    std = float(ensemble.std(ddof=1))
    z = 0.0 if std <= 0.0 else (mean - error) / std
    return error, mean, std, float(z)


def delay_embedding(series: Any, *, delay: int, dimension: int) -> np.ndarray:
    """Return the ``(M, dimension)`` delay-coordinate embedding of a scalar series.

    Args:
        series: Observable samples.
        delay: Lag between successive coordinates, in samples.
        dimension: Embedding dimension.

    Returns:
        The embedded points, ``M = N - (dimension - 1) * delay``.

    Raises:
        ValueError: If the series is too short for the requested embedding.
    """
    x = _as_series(series)
    lag = _positive_int(delay, what="delay")
    dim = _positive_int(dimension, what="dimension", minimum=2)
    span = (dim - 1) * lag
    points = x.size - span
    if points < 4:
        raise ValueError(
            f"series of {x.size} samples cannot support dimension {dim} at delay {lag}"
        )
    return np.stack([x[i * lag : i * lag + points] for i in range(dim)], axis=1)


def _leading_window(points: np.ndarray, limit: int) -> np.ndarray:
    """Return the first ``limit`` embedded points — a CONTIGUOUS window, never a strided subsample.

    Striding is the obvious way to cap an ``O(M^2)`` cost and it is wrong here.  Both statistics below
    read the trajectory's temporal structure: a stride of ``s`` multiplies the per-sample Lyapunov
    exponent by ``s``, so a well-resolved chaotic flow can be strided past its own predictability
    horizon and come back looking like noise, and the recurrence diagonals it would have formed are
    dismantled the same way.  A contiguous window analyses less of the record but analyses it at the
    resolution it was recorded at; the analysed count is reported so the truncation is visible.
    """
    if points.shape[0] <= limit:
        return points
    return points[:limit]


def recurrence_determinism(
    embedded: np.ndarray,
    *,
    recurrence_rate_target: float,
    min_diagonal_length: int,
    theiler_window: int,
    max_points: int,
) -> tuple[float, float, int]:
    """Return ``(recurrence_rate, determinism, max_diagonal_length)`` from a recurrence plot.

    Determinism is the RQA ratio: of all recurrence points, the fraction lying on diagonal lines of at
    least ``min_diagonal_length``.  A deterministic trajectory revisits a neighbourhood and then *stays*
    parallel to its earlier pass for a while, so its recurrence points line up diagonally; a stochastic
    one recurs in isolated dots.

    This is the statistic the spectrum cannot supply.  The fully developed logistic map has a flat power
    spectrum and a delta-function autocorrelation — spectrally identical to white noise — and only its
    recurrence structure reveals that it is deterministic at all.

    Two normalisations keep the statistic from measuring the wrong thing:

    * **The threshold is a distance QUANTILE**, chosen to achieve ``recurrence_rate_target``, not an
      absolute distance — so the result is invariant to the observable's units and scale.
    * **Diagonals within the Theiler window of the main diagonal are excluded.**  Without this, any
      smoothly varying signal is trivially "deterministic": consecutive samples of an
      Ornstein-Uhlenbeck path are close simply because they are consecutive, and the near-diagonal band
      that produces is pure sampling rate.  The OU control measured determinism 0.78 without the
      exclusion — indistinguishable from the Lorenz control — and 0.09 with it.

    Args:
        embedded: ``(M, d)`` delay-embedded points.
        recurrence_rate_target: Target fraction of recurrent pairs, from the caller's contract.
        min_diagonal_length: Shortest diagonal counted as deterministic.
        theiler_window: Diagonals at offsets at or below this are discarded as trivial recurrence.
        max_points: Cap on ``M`` entering the ``O(M^2)`` computation.

    Returns:
        The achieved recurrence rate over the admitted pairs, the determinism in ``[0, 1]``, and the
        longest admitted diagonal line.
    """
    points = _leading_window(np.asarray(embedded, dtype=np.float64), max_points)
    count = points.shape[0]
    theiler = max(1, int(theiler_window))
    if count < 4 or theiler >= count - 1:
        return 0.0, 0.0, 0
    diff = points[:, None, :] - points[None, :, :]
    distance = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
    index = np.arange(count)
    admitted = np.abs(index[:, None] - index[None, :]) > theiler
    candidate = distance[admitted]
    if candidate.size == 0:
        return 0.0, 0.0, 0
    epsilon = float(np.quantile(candidate, recurrence_rate_target))
    recurrent = (distance <= epsilon) & admitted
    rate = float(np.count_nonzero(recurrent)) / float(candidate.size)

    total_points = 0
    deterministic_points = 0
    longest = 0
    for offset in range(theiler + 1, count):
        line = recurrent.diagonal(offset)
        if not line.any():
            continue
        padded = np.concatenate(([False], line, [False]))
        edges = np.flatnonzero(padded[1:] != padded[:-1])
        lengths = edges[1::2] - edges[0::2]
        total_points += int(lengths.sum())
        deterministic_points += int(lengths[lengths >= min_diagonal_length].sum())
        longest = max(longest, int(lengths.max()))
    determinism = 0.0 if total_points == 0 else float(deterministic_points) / float(total_points)
    return rate, determinism, longest


def largest_lyapunov_exponent(
    embedded: np.ndarray,
    dt: float,
    *,
    theiler_window: int,
    fit_lags: int,
    max_points: int,
) -> tuple[float, float, int]:
    """Return ``(lambda_per_s, fit_r2, fit_lags_used)`` by Rosenstein's mean-log-divergence slope.

    For every embedded point the nearest neighbour at least ``theiler_window`` samples away in time is
    found, and the mean of ``ln||X_{i+k} - X_{j+k}||`` is tracked against ``k``.  For a chaotic
    trajectory that mean rises linearly at the largest Lyapunov exponent and then **saturates** at the
    attractor diameter; for a periodic or relaxing one it is flat from the start.

    Two details are not optional and both were forced by the controls rather than assumed:

    * **The Theiler window.**  Without it the nearest neighbour of every point is its immediate
      predecessor and the slope measures the sampling rate, not the dynamics.
    * **The fit region is DERIVED, not searched and not fixed.**  The linear stretch ends at the
      predictability horizon ``ln(diameter / initial separation) / lambda`` — a handful of samples for a
      strongly chaotic map, hundreds for a slowly diverging flow.  A lag count fixed by the caller
      therefore averages the saturated plateau into the slope: on the logistic control a 20-lag fit
      returned 0.0027 per iteration against a true ``ln(2)``-scale value.  Searching for the most linear
      sub-window is worse, not better: the shortest admissible window is almost always the most linear,
      so the search collapses onto three points and reports the noise-decorrelation jump as the
      exponent, which is how the OU control produced 2.77/s.  The caller instead declares the window in
      units of the series' OWN measured correlation time and this function receives the resolved lag
      count — a rule fixed before the run, evaluated against a measured scale.

    Args:
        embedded: ``(M, d)`` delay-embedded points.
        dt: Physical time between samples [s].
        theiler_window: Minimum temporal separation between a point and its neighbour, in samples.
        fit_lags: Number of lags in the fit, resolved by the caller from the correlation time.
        max_points: Cap on ``M``; a contiguous leading window, never a stride.

    Returns:
        The slope [1/s], the ``R^2`` of the fit, and how many lags it used.  A trajectory with no
        admissible neighbour pairs returns ``(0.0, 0.0, 0)`` — a refusal to estimate, not a claim of
        zero.
    """
    dt_s = _positive(dt, what="dt")
    points = _leading_window(np.asarray(embedded, dtype=np.float64), max_points)
    count = points.shape[0]
    lags = int(min(fit_lags, count // 4))
    if count < 16 or lags < 3:
        return 0.0, 0.0, 0

    diff = points[:, None, :] - points[None, :, :]
    distance = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
    index = np.arange(count)
    excluded = np.abs(index[:, None] - index[None, :]) <= max(1, theiler_window)
    distance[excluded] = np.inf
    neighbour = np.argmin(distance, axis=1)
    valid = np.isfinite(distance[index, neighbour]) & (distance[index, neighbour] > 0.0)
    if not np.any(valid):
        return 0.0, 0.0, 0

    divergence = np.full(lags, np.nan, dtype=np.float64)
    i_all = index[valid]
    j_all = neighbour[valid]
    for k in range(1, lags + 1):
        keep = (i_all + k < count) & (j_all + k < count)
        if not np.any(keep):
            continue
        separation = np.linalg.norm(points[i_all[keep] + k] - points[j_all[keep] + k], axis=1)
        separation = separation[separation > 0.0]
        if separation.size == 0:
            continue
        divergence[k - 1] = float(np.mean(np.log(separation)))

    times = np.arange(1, lags + 1, dtype=np.float64) * dt_s
    finite = np.isfinite(divergence)
    if np.count_nonzero(finite) < 3:
        return 0.0, 0.0, 0
    t = times[finite]
    v = divergence[finite]
    slope, intercept = np.polyfit(t, v, 1)
    residual = float(np.sum((v - (slope * t + intercept)) ** 2))
    total = float(np.sum((v - v.mean()) ** 2))
    r2 = 1.0 if total <= 0.0 else max(0.0, 1.0 - residual / total)
    return float(slope), float(r2), int(lags)


def _commensurate(ratio: float, *, max_denominator: int, tolerance: float) -> bool:
    """Return whether ``ratio`` is within ``tolerance`` (relative) of a rational with a small denominator."""
    if not math.isfinite(ratio) or ratio <= 0.0:
        return False
    approximation = Fraction(ratio).limit_denominator(max_denominator)
    if approximation == 0:
        return False
    return abs(ratio - float(approximation)) / ratio <= tolerance


def _fundamentals(
    peaks: Sequence[SpectralPeak], *, resolution_hz: float, contract: RegimeEvidenceContract
) -> tuple[SpectralPeak, ...]:
    """Return the peaks that are neither harmonics of nor part of the strongest line, strongest first.

    Two reductions, both of which the controls forced:

    * **Harmonics are removed.**  Their tolerance is floored at the spectral resolution: a peak cannot
      be distinguished from ``k * f0`` more finely than the grid it was measured on, and pretending
      otherwise splits one cycle into two "independent" frequencies.
    * **Peaks inside an already-counted line's own half-power width are absorbed into it.**  A broad,
      noise-sustained line is not smooth — the quasicycle control's single hump at 0.30 Hz resolves into
      local maxima at 0.2875 and 0.3063 Hz, whose ratio 1.065 is not a small rational, so without this
      merge a quasicycle is reported as QUASIPERIODIC.  A second frequency is only a second frequency if
      it is resolved *away* from the first.
    """
    if not peaks:
        return ()
    fundamental = peaks[0]
    kept = [fundamental]
    for peak in peaks[1:]:
        if peak.power < contract.secondary_peak_power_fraction * fundamental.power:
            continue
        if any(
            abs(peak.frequency_hz - other.frequency_hz) <= 0.5 * (other.width_hz + peak.width_hz)
            for other in kept
        ):
            continue
        ratio = peak.frequency_hz / fundamental.frequency_hz
        order = round(ratio)
        if order >= 2:
            tolerance = max(
                contract.harmonic_relative_tolerance * order * fundamental.frequency_hz,
                resolution_hz,
            )
            if abs(peak.frequency_hz - order * fundamental.frequency_hz) <= tolerance:
                continue
        kept.append(peak)
    return tuple(kept)


def classify_regime(
    series: Any,
    dt: float,
    *,
    observable: str,
    amplitude_scale: float,
    contract: RegimeEvidenceContract,
    equilibration_candidates: int = 40,
) -> RegimeReport:
    """Classify one accepted trajectory, returning ``UNDECIDED`` whenever the evidence does not decide.

    The decision order is fixed and documented here because a classifier whose order is implicit is a
    classifier nobody can audit:

    1. **Enough samples at all?**  Below ``contract.min_samples``, ``UNDECIDED``.
    2. **Transient discarded, then still drifting?**  A drift exceeding ``contract.drift_sigma`` of the
       window's own standard deviation is a transient, not a regime — ``UNDECIDED``, routing to a longer
       run.  Note this is *not* a convergence test: a drifting trajectory is not a failed one, it is one
       whose regime cannot yet be read.
    3. **Enough independent samples?**  Below ``contract.min_n_eff``, or a window covering fewer than
       ``contract.min_correlation_windows`` of the measured ``tau_int``, ``UNDECIDED``.
    4. **Fluctuation negligible against the declared scale?**  ``FIXED_POINT``.
    5. **Deterministic (by recurrence AND against phase-randomised surrogates) with no resolved
       spectral line?**  ``CHAOS``.  Placed before the featureless-spectrum test on purpose: a fully
       developed chaotic map has a flat spectrum and a delta-function ACF, so testing "no structure"
       first would file it as noise.
    6. **No temporal structure at all?**  Uncorrelated and spectrally flat with no peak — ``UNDECIDED``,
       because an uncorrelated sequence and a NESS sampled far slower than its own correlation time are
       the same series, and choosing between them needs a faster sampling rate, not a bolder classifier.
    7. **A significant spectral peak?**  Two or more incommensurate fundamentals give ``QUASIPERIODIC``;
       one sharp fundamental gives ``LIMIT_CYCLE``; one broad one gives ``QUASICYCLE``; a width between
       the declared bands gives ``UNDECIDED``.
    8. **Otherwise** — correlated, stationary, no peak — ``STOCHASTIC_NESS``.

    Args:
        series: Observable samples from an accepted trajectory, evenly spaced.
        dt: Physical time between samples [s].  This is the SAMPLE spacing, not the solver step: with a
            telemetry stride the two differ, and the wrong one understates every correlation time by
            exactly that stride.
        observable: Name, recorded in the report.
        amplitude_scale: The observable's declared physical scale in its own units — the resting value,
            the intervention amplitude, or whatever the caller pre-registered.  Only the fixed-point
            test uses it, and it must be declared before the run: choosing it afterwards is how a
            fluctuating trajectory becomes a fixed point.
        contract: The caller's pre-declared thresholds.
        equilibration_candidates: Truncation points tried when locating the end of the transient.

    Returns:
        The :class:`RegimeReport`, always carrying full statistics regardless of verdict.

    Raises:
        ValueError: If ``dt`` or ``amplitude_scale`` is not positive-finite, or the series is empty or
            non-finite.
        TypeError: If ``contract`` is not a :class:`RegimeEvidenceContract`.
    """
    if not isinstance(contract, RegimeEvidenceContract):
        raise TypeError("contract must be a RegimeEvidenceContract declared before the run")
    dt_s = _positive(dt, what="dt")
    scale = _positive(amplitude_scale, what="amplitude_scale")
    if not isinstance(observable, str) or not observable.strip():
        raise ValueError("observable must be a non-empty string")
    x = _as_series(series)
    contract_dict = contract.to_dict()

    def _report(
        regime: Regime, reason: str, statistics: RegimeStatistics, expansions: tuple[str, ...] = ()
    ) -> RegimeReport:
        return RegimeReport(
            regime=regime,
            observable=observable.strip(),
            reason=reason,
            statistics=statistics,
            contract=contract_dict,
            expansions=expansions,
        )

    if x.size < contract.min_samples:
        return _report(
            Regime.UNDECIDED,
            f"{x.size} samples against a declared minimum of {contract.min_samples}",
            _degenerate_statistics(x, dt_s, scale, contract),
            ("time-modality: run longer before asking what regime this is",),
        )

    start = equilibration_point(
        x, dt_s, sokal_window_c=contract.sokal_window_c, n_candidates=equilibration_candidates
    )
    tail = x[start:]
    if tail.size < contract.min_samples:
        start, tail = 0, x
    analysed_s = float(tail.size * dt_s)

    std = float(np.std(tail))
    mean = float(np.mean(tail))
    times = np.arange(tail.size, dtype=np.float64) * dt_s
    slope = float(np.polyfit(times, tail, 1)[0]) if std > 0.0 else 0.0
    drift_over_std = abs(slope * analysed_s) / std if std > 0.0 else 0.0
    relative_amplitude = std / scale

    tau_s, n_eff, tau_note = integrated_autocorrelation_time(
        tail, dt_s, sokal_window_c=contract.sokal_window_c
    )

    # `_as_series` guarantees every SAMPLE is finite, which is not the same as every STATISTIC being
    # finite.  A near-overflow explicit scheme produces samples around 1e200 that are all valid
    # doubles, and then the second moment overflows: the variance is `inf`, the normalised
    # autocorrelation is `nan`, and `tau_s` follows.  `nan` does not corrupt the ANSWER here — it
    # reaches `round(nan)` in the Theiler window below and raises `cannot convert float NaN to
    # integer`, sixty lines from the cause and reading like a library bug.  A diverged run is a
    # legitimate thing to hand a classifier and the honest reply is a verdict: the divergence IS the
    # finding.  Found by the S0 invariance battery, which needs exactly this case to flag one.
    if not (math.isfinite(tau_s) and math.isfinite(n_eff)):
        return _report(
            Regime.UNDECIDED,
            f"the correlation time is not finite (tau_int={tau_s}, n_eff={n_eff}) although every "
            f"sample is: the second moment overflowed, so this trajectory has no readable regime",
            _degenerate_statistics(tail, dt_s, scale, contract),
            (
                "solver-design: classify the scheme, not the series — the divergence is the result",
                "numerical-invariance: run the S0 ladder to find the step at which it diverges",
            ),
        )

    acf = normalized_autocorrelation(tail)
    zero = np.flatnonzero(acf <= 0.0)
    first_zero_s = float(zero[0] * dt_s) if zero.size else None
    rebound = float(np.max(np.abs(acf[zero[0] :]))) if zero.size else 0.0

    frequencies, psd, resolution = welch_psd(tail, dt_s, n_segments=contract.welch_segment_count)
    flatness = spectral_flatness(psd)
    peaks = spectral_peaks(
        frequencies,
        psd,
        resolution,
        prominence_ratio_min=contract.peak_prominence_ratio_min,
        baseline_octave_ratio=contract.baseline_octave_ratio,
        power_fraction_min=contract.peak_power_fraction_min,
    )
    fundamentals = _fundamentals(peaks, resolution_hz=resolution, contract=contract)

    ratio: float | None = None
    commensurate: bool | None = None
    if len(fundamentals) >= 2:
        high = max(fundamentals[0].frequency_hz, fundamentals[1].frequency_hz)
        low = min(fundamentals[0].frequency_hz, fundamentals[1].frequency_hz)
        ratio = high / low if low > 0.0 else None
        if ratio is not None:
            commensurate = _commensurate(
                ratio,
                max_denominator=contract.max_rational_denominator,
                tolerance=contract.harmonic_relative_tolerance,
            )

    # Cap the delay so the embedding always fits with room for the estimators.  A repeated-column
    # fallback would be a FABRICATED embedding — every coordinate identical, hence perfect trivial
    # recurrence — so the span is bounded up front instead of being rescued after it fails.
    max_span_delay = max(1, (tail.size - 4) // max(1, contract.embedding_dimension - 1))
    delay = min(embedding_delay(acf, max_delay=max(1, tail.size // 8)), max_span_delay)
    embedded = delay_embedding(tail, delay=delay, dimension=contract.embedding_dimension)
    theiler = int(max(1, round(contract.theiler_tau_multiple * tau_s / dt_s)))
    recurrence_rate, determinism, longest = recurrence_determinism(
        embedded,
        recurrence_rate_target=contract.recurrence_rate_target,
        min_diagonal_length=contract.min_diagonal_length,
        theiler_window=theiler,
        max_points=contract.recurrence_max_points,
    )
    resolved_fit_lags = int(
        min(
            contract.lyapunov_max_fit_lags,
            max(
                contract.lyapunov_min_fit_lags,
                round(contract.lyapunov_fit_tau_multiple * tau_s / dt_s),
            ),
        )
    )
    lyapunov, lyapunov_r2, lyapunov_lags = largest_lyapunov_exponent(
        embedded,
        dt_s,
        theiler_window=theiler,
        fit_lags=resolved_fit_lags,
        max_points=contract.recurrence_max_points,
    )
    horizon = int(max(1, round(tau_s / dt_s)))
    error, surrogate_mean, surrogate_std, predictability_z = nonlinear_predictability(
        tail,
        delay=delay,
        dimension=contract.embedding_dimension,
        theiler=theiler,
        horizon=horizon,
        neighbours=contract.predictability_neighbours,
        max_points=contract.recurrence_max_points,
        n_surrogates=contract.n_surrogates,
        seed=contract.surrogate_seed,
    )

    statistics = RegimeStatistics(
        n_samples=int(x.size),
        dt_s=dt_s,
        analysed_time_s=analysed_s,
        equilibration_time_s=float(start * dt_s),
        mean=mean,
        std=std,
        relative_amplitude=relative_amplitude,
        drift_per_s=slope,
        drift_over_std=drift_over_std,
        tau_int_s=tau_s,
        n_eff=n_eff,
        tau_note=tau_note,
        acf_first_zero_crossing_s=first_zero_s,
        acf_rebound=rebound,
        spectral_resolution_hz=float(resolution),
        spectral_flatness=flatness,
        peaks=peaks,
        fundamental_frequencies_hz=tuple(peak.frequency_hz for peak in fundamentals),
        frequency_ratio=ratio,
        ratio_is_commensurate=commensurate,
        embedding_delay_samples=int(delay),
        embedding_dimension=int(contract.embedding_dimension),
        embedded_points=int(embedded.shape[0]),
        embedded_points_analysed=int(min(embedded.shape[0], contract.recurrence_max_points)),
        recurrence_rate=recurrence_rate,
        determinism=determinism,
        max_diagonal_length=longest,
        lyapunov_per_s=lyapunov,
        lyapunov_fit_r2=lyapunov_r2,
        lyapunov_fit_lags_used=lyapunov_lags,
        prediction_horizon_samples=horizon,
        prediction_error=error,
        surrogate_prediction_error_mean=surrogate_mean,
        surrogate_prediction_error_std=surrogate_std,
        predictability_z=predictability_z,
    )

    # 2. Still a transient.  Not a failure — a trajectory whose regime is not yet readable.
    if drift_over_std >= contract.drift_sigma:
        return _report(
            Regime.UNDECIDED,
            (
                f"drift across the analysed window is {drift_over_std:.3g}x its own fluctuation "
                f"(declared limit {contract.drift_sigma:g}); this is a transient, and a transient has "
                "no regime yet — it is not a failed run"
            ),
            statistics,
            ("time-modality: extend the accepted trajectory past the transient",),
        )

    # 3. Not enough independent evidence for any statistic below to mean anything.
    if n_eff < contract.min_n_eff or analysed_s < contract.min_correlation_windows * tau_s:
        return _report(
            Regime.UNDECIDED,
            (
                f"n_eff={n_eff:.3g} against a declared minimum {contract.min_n_eff:g}, and the window "
                f"covers {analysed_s / tau_s:.3g} correlation times against a declared "
                f"{contract.min_correlation_windows:g}"
            ),
            statistics,
            (
                "time-modality: the run must cover the transient PLUS the window, not the window alone",
            ),
        )

    # 4. Fluctuation negligible against the caller's declared scale.
    if relative_amplitude < contract.fixed_point_relative_amplitude:
        return _report(
            Regime.FIXED_POINT,
            (
                f"std/amplitude_scale = {relative_amplitude:.3g} below the declared "
                f"{contract.fixed_point_relative_amplitude:g}"
            ),
            statistics,
        )

    # 5. Deterministic, aperiodic, diverging.  Before the featureless test: a fully developed chaotic
    #    map is spectrally white, so "no structure" would otherwise capture it.  Spectral FLATNESS is
    #    deliberately not part of this test — the Lorenz control measures 4.1e-9, because a smooth flow
    #    has almost no high-frequency power, so flatness reads a continuous chaotic spectrum as
    #    "extremely peaked".  Aperiodicity is therefore tested as the ABSENCE of a resolved line.
    deterministic = (
        determinism >= contract.determinism_min
        and predictability_z >= contract.predictability_z_min
    )
    if deterministic and not fundamentals:
        return _report(
            Regime.CHAOS,
            (
                f"no resolved spectral line; recurrence determinism {determinism:.3g} >= "
                f"{contract.determinism_min:g} and nonlinear prediction beats its own "
                f"phase-randomised surrogates by z={predictability_z:.3g} >= "
                f"{contract.predictability_z_min:g} (error {error:.3g} against surrogate "
                f"{surrogate_mean:.3g}). Largest Lyapunov {lyapunov:.3g}/s is reported, not gated on"
            ),
            statistics,
        )

    # 6. No temporal structure at all.
    if (
        tau_s <= contract.white_noise_tau_dt_max * dt_s
        and flatness >= contract.white_noise_flatness_min
        and not peaks
    ):
        return _report(
            Regime.UNDECIDED,
            (
                f"tau_int={tau_s / dt_s:.3g} dt and spectral flatness {flatness:.3g} with no resolved "
                "peak: an uncorrelated sequence is indistinguishable from a NESS sampled far slower "
                "than its own correlation time, and no additional statistic on THIS series separates "
                "them"
            ),
            statistics,
            (
                "time-modality: sample faster than the process, not slower",
                "observation-model extension: a second, faster-sampled observable would break the tie",
            ),
        )

    # 7. Oscillation.
    if fundamentals:
        # A second frequency counts toward quasiperiodicity only if BOTH lines are resolution-limited
        # sharp.  A broad line is not a well-defined frequency, and the Welch estimate of a broad
        # resonance breaks into several local maxima whose ratios are arbitrary: on the quasicycle
        # control, sidebands of one 0.29 Hz hump were read as an incommensurate pair at ratio 1.37.
        sharp = [
            peak for peak in fundamentals if peak.width_bins <= contract.limit_cycle_max_width_bins
        ]
        if len(sharp) >= 2:
            high = max(sharp[0].frequency_hz, sharp[1].frequency_hz)
            low = min(sharp[0].frequency_hz, sharp[1].frequency_hz)
            sharp_ratio = high / low if low > 0.0 else None
            if sharp_ratio is not None and not _commensurate(
                sharp_ratio,
                max_denominator=contract.max_rational_denominator,
                tolerance=contract.harmonic_relative_tolerance,
            ):
                return _report(
                    Regime.QUASIPERIODIC,
                    (
                        f"two resolution-limited fundamentals at {sharp[0].frequency_hz:.4g} Hz and "
                        f"{sharp[1].frequency_hz:.4g} Hz with ratio {sharp_ratio:.4g}, not within "
                        f"{contract.harmonic_relative_tolerance:g} of any rational with denominator "
                        f"<= {contract.max_rational_denominator}"
                    ),
                    statistics,
                )
        dominant = fundamentals[0]
        if dominant.width_bins <= contract.limit_cycle_max_width_bins:
            return _report(
                Regime.LIMIT_CYCLE,
                (
                    f"one fundamental at {dominant.frequency_hz:.4g} Hz, half-power width "
                    f"{dominant.width_bins:.3g} resolution bins <= "
                    f"{contract.limit_cycle_max_width_bins:g} (Q={dominant.quality_factor:.3g})"
                ),
                statistics,
            )
        if dominant.width_bins >= contract.quasicycle_min_width_bins:
            if deterministic:
                return _report(
                    Regime.UNDECIDED,
                    (
                        f"a broad line at {dominant.frequency_hz:.4g} Hz "
                        f"({dominant.width_bins:.3g} bins) together with determinism "
                        f"{determinism:.3g} and surrogate-referenced predictability "
                        f"z={predictability_z:.3g}: a noise-sustained quasicycle and phase-coherent "
                        "chaos produce the same broad peak, and no statistic on this passive scalar "
                        "record separates them"
                    ),
                    statistics,
                    (
                        "latent coordinate: a second, mechanistically independent observable would "
                        "resolve whether the phase diffuses or the amplitude folds",
                        "intervention modality: the decay of a phase perturbation separates them "
                        "directly; a passive record does not",
                    ),
                )
            return _report(
                Regime.QUASICYCLE,
                (
                    f"one fundamental at {dominant.frequency_hz:.4g} Hz, half-power width "
                    f"{dominant.width_bins:.3g} resolution bins >= "
                    f"{contract.quasicycle_min_width_bins:g}: noise-sustained, not self-sustained "
                    f"(Q={dominant.quality_factor:.3g})"
                ),
                statistics,
            )
        return _report(
            Regime.UNDECIDED,
            (
                f"peak width {dominant.width_bins:.3g} bins falls in the declared UNDECIDED band "
                f"({contract.limit_cycle_max_width_bins:g}, {contract.quasicycle_min_width_bins:g}): "
                "a noise-sustained and a self-sustained oscillation are not separated at this "
                "resolution"
            ),
            statistics,
            (
                "time-modality: a longer record narrows the resolution and may separate them",
                "intervention modality: a perturbation-decay experiment separates a self-sustained "
                "cycle from a noise-sustained one directly; a passive record does not",
            ),
        )

    # 8. Correlated, stationary, featureless spectrum.
    return _report(
        Regime.STOCHASTIC_NESS,
        (
            f"stationary with a decaying autocorrelation (tau_int={tau_s:.4g} s, n_eff={n_eff:.4g}) "
            f"and no spectral peak above {contract.peak_prominence_ratio_min:g}x its baseline; "
            "sustained fluctuation about a mean is an accepted state, not a convergence failure"
        ),
        statistics,
    )


def _degenerate_statistics(
    x: np.ndarray, dt_s: float, scale: float, contract: RegimeEvidenceContract
) -> RegimeStatistics:
    """Return the statistics block for a series too short to run the estimators on."""
    std = float(np.std(x))
    return RegimeStatistics(
        n_samples=int(x.size),
        dt_s=dt_s,
        analysed_time_s=float(x.size * dt_s),
        equilibration_time_s=0.0,
        mean=float(np.mean(x)),
        std=std,
        relative_amplitude=std / scale,
        drift_per_s=float("nan"),
        drift_over_std=float("nan"),
        tau_int_s=float("nan"),
        n_eff=0.0,
        tau_note="series shorter than the declared minimum; no estimator was run",
        acf_first_zero_crossing_s=None,
        acf_rebound=float("nan"),
        spectral_resolution_hz=float("nan"),
        spectral_flatness=float("nan"),
        peaks=(),
        fundamental_frequencies_hz=(),
        frequency_ratio=None,
        ratio_is_commensurate=None,
        embedding_delay_samples=0,
        embedding_dimension=int(contract.embedding_dimension),
        embedded_points=0,
        embedded_points_analysed=0,
        recurrence_rate=float("nan"),
        determinism=float("nan"),
        max_diagonal_length=0,
        lyapunov_per_s=float("nan"),
        lyapunov_fit_r2=float("nan"),
        lyapunov_fit_lags_used=0,
        prediction_horizon_samples=0,
        prediction_error=float("nan"),
        surrogate_prediction_error_mean=float("nan"),
        surrogate_prediction_error_std=float("nan"),
        predictability_z=float("nan"),
    )
