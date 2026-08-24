"""A parameter sweep that refuses the differences it is not entitled to take.

WHY THIS MODULE EXISTS.  The brute-force grid sweep was abandoned because it is unaffordable, and the
arithmetic in ``ROADMAP.md`` says so by three to five orders of magnitude.  But cost was never the only
thing wrong with it.  A grid sweep computes ``(f(p + h) - f(p)) / h`` between two accepted trajectories
and calls the result a sensitivity, and that number is meaningless in at least four situations the
sweep had no way to detect:

1. **The two points are in different dynamical regimes.**  A difference taken across a Hopf bifurcation
   is not a derivative of anything; it is the gap between the mean of a steady state and the mean of an
   orbit, and it will change sign, magnitude and units with the step size because there is no
   underlying smooth function.  :func:`finite_difference_sensitivity` **REFUSES** and names the
   expansion that refusal opens — partition the parameter space by regime, or find the latent
   coordinate that makes the response smooth again.
2. **The error bar was computed from the step count.**  An accepted trajectory is correlated; its
   ``N`` samples are worth ``n_eff = T / (2 tau_int)`` independent ones, and on this project's own
   time scales the two differ by more than an order of magnitude.  ``sigma / sqrt(N)`` is the standard
   way to publish an interval that is far too narrow, so every interval here is ESS-based and carries
   the widening factor it applied.
3. **The observable was time-averaged over an orbit.**  The mean of a limit cycle is a point the
   trajectory passes through twice per period and never rests at.  :func:`phase_decomposition`
   separates the phase-resolved profile from the scalar mean so a cycle cannot be silently reported as
   its centre.
4. **The ensemble was quietly malformed.**  A member missing an observable, two members sharing a seed,
   a mixture whose weights do not sum to one — each of these produces a number rather than an error,
   and each of the three is a different lie.  :func:`audit_ensemble` detects all three and reports them
   as typed defects rather than averaging around them.

THE JENSEN TERM IS NOT A BUG.  ``E[f(X)] != f(E[X])`` for any nonlinear ``f``.  When this layer replaces
a point state with a distribution, that inequality becomes a real, measurable quantity, and it is the
**price of the closure**, not evidence that the reference dynamics are wrong.  :func:`moment_closure_gap`
measures it, names it, gives it its own declared tolerance, and is verified against a case where the
answer is known in closed form.  Reading a nonzero gap as a defect in the underlying model would be
exactly backwards: the gap is a property of choosing to propagate moments instead of samples.

AUTHORITY.  Like :mod:`aleph.virtual_cell.regime`, this is an **unverified observer**.  It grants no
evidence rung.  Nothing here imports Warp, touches a device, or self-declares native provenance.

units: time in SECONDS of physical time.  A sensitivity carries the observable's units divided by the
parameter's, and both names are recorded because a bare number is not a sensitivity.

Sanity Gate:
    * dimensional: ``dt`` is [s] and required; a sensitivity is [observable unit / parameter unit] and
      both unit strings are carried in the verdict; ``n_eff`` and ``widening_factor`` are dimensionless.
    * boundary: a finite difference with a zero (or non-finite) parameter separation raises rather than
      dividing; an empty ensemble raises rather than returning an empty summary; a phase decomposition
      over less than one period, or with fewer samples than bins, refuses rather than emitting an empty
      bin as a zero; an ESS interval needs at least four samples.
    * conservation/invariant: ``n_eff <= n_samples`` and ``widening_factor >= 1`` always — an ESS
      interval is never narrower than the naive one.  Mixture weights must sum to one within the
      declared tolerance and none may be negative.  The phase-bin sample counts must sum to the number
      of samples used.
    * numerical: intervals are formed as ``mean +/- z * sem_ess``; ``z`` is the caller's declared
      coverage factor, recorded in the interval, so a reader knows whether "the interval" is 1 sigma or
      95 %.
    * sign-sense: a sensitivity is signed and reported signed, with the low and high points named, so
      the sign is interpretable without knowing which point was passed first.
    * measurement-protocol: every threshold — the regime contract, the ESS floor, the coverage factor,
      the closure tolerance, the weight tolerance — is caller-supplied and copied into the verdict.  A
      refusal always names the expansion it opens; a refusal with no route is how a hard problem gets
      quietly dropped instead of scheduled.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

import numpy as np

from aleph.virtual_cell.regime import (
    EVIDENCE_AUTHORITY,
    Regime,
    RegimeEvidenceContract,
    RegimeReport,
    classify_regime,
    integrated_autocorrelation_time,
)

__all__ = [
    "ClosureGapReport",
    "EnsembleDefect",
    "EnsembleDefectRecord",
    "EnsembleMember",
    "EnsembleSummary",
    "ESSInterval",
    "ExpansionRoute",
    "PhaseDecomposition",
    "SensitivityVerdict",
    "SweepPoint",
    "SweepVerdictKind",
    "audit_ensemble",
    "ess_interval",
    "finite_difference_sensitivity",
    "moment_closure_gap",
    "phase_decomposition",
    "summarize_ensemble",
]


class ExpansionRoute(StrEnum):
    """Where a refusal sends the programme.  A failure is never a reason to delete scope.

    The first five mirror the shared failure-routing table; the rest name the routes this lane actually
    reaches.  A :class:`SensitivityVerdict` that refuses without at least one of these is malformed.
    """

    REGIME_PARTITION = "regime-partition"
    LATENT_COORDINATE = "latent-coordinate"
    ENSEMBLE = "ensemble"
    TIME_MODALITY = "time-modality"
    FORCE_MODALITY = "force-modality"
    INTERVENTION_MODALITY = "intervention-modality"
    MISSING_COMPONENT = "missing-component"
    MISSING_LAW = "missing-law"
    SOLVER_DESIGN = "solver-design"
    ADAPTIVE_DESIGN = "adaptive-design"
    OBSERVATION_MODEL_EXTENSION = "observation-model-extension"
    MODEL_INADEQUACY = "model-inadequacy"
    PROVENANCE_REPAIR = "provenance-repair"


class SweepVerdictKind(StrEnum):
    """Whether a sweep quantity was reported, and if not, why not."""

    REPORTED = "reported"
    REFUSED_CROSS_REGIME = "refused-cross-regime"
    REFUSED_UNDECIDED_REGIME = "refused-undecided-regime"
    REFUSED_INSUFFICIENT_ESS = "refused-insufficient-ess"


class EnsembleDefect(StrEnum):
    """The ways an ensemble is malformed in a manner that still produces a number.

    Each of these has a distinct consequence and they must not be collapsed:

    * ``MISSING_OBSERVABLE`` — a member lacks an observable its siblings have.  Averaging over the
      members that DO have it silently changes the estimator's population between observables, so two
      numbers in the same table are then means over different ensembles.
    * ``DUPLICATE_SEED`` — two members ran the same stream.  They are not independent, so the spread is
      too small and the apparent sample size is inflated by exactly the duplication.
    * ``MISSING_SEED`` — a member carries no seed at all.  Its independence is unknowable, which is not
      the same as its being independent.
    * ``NEGATIVE_WEIGHT`` — a mixture weight below zero is not a proportion of anything.
    * ``WEIGHTS_NOT_NORMALIZED`` — weights that do not sum to one silently rescale every mixed mean.
    * ``PARTIAL_WEIGHTS`` — some members weighted and others not; there is no defensible default for
      the rest.
    """

    MISSING_OBSERVABLE = "missing-observable"
    DUPLICATE_SEED = "duplicate-seed"
    MISSING_SEED = "missing-seed"
    NEGATIVE_WEIGHT = "negative-weight"
    WEIGHTS_NOT_NORMALIZED = "weights-not-normalized"
    PARTIAL_WEIGHTS = "partial-weights"


def _positive(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{what} must be positive and finite; got {value!r}")
    return number


def _finite(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite; got {value!r}")
    return number


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _series(values: Any, *, what: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64).ravel()
    if array.size == 0:
        raise ValueError(f"{what} must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{what} must contain only finite values")
    return array


@dataclass(frozen=True, slots=True)
class ESSInterval:
    """A mean and its interval, computed from INDEPENDENT samples rather than from the step count.

    ``sem_ess = std * sqrt(2 * tau_int / T)`` is the standard error of a correlated mean;
    ``sem_naive = std / sqrt(N)`` is what a sweep reports when it forgets that its samples came from a
    trajectory.  ``widening_factor`` is their ratio and is never below 1 — it is the amount by which the
    naive interval was too narrow, and it is carried so that a reader can see it rather than infer it.
    """

    observable: str
    mean: float
    std: float
    n_samples: int
    tau_int_s: float
    n_eff: float
    sem_naive: float
    sem_ess: float
    widening_factor: float
    coverage_z: float
    low: float
    high: float
    tau_note: str

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able interval record."""
        return {
            "observable": self.observable,
            "mean": self.mean,
            "std": self.std,
            "n_samples": self.n_samples,
            "tau_int_s": self.tau_int_s,
            "n_eff": self.n_eff,
            "sem_naive": self.sem_naive,
            "sem_ess": self.sem_ess,
            "widening_factor": self.widening_factor,
            "coverage_z": self.coverage_z,
            "low": self.low,
            "high": self.high,
            "tau_note": self.tau_note,
            "sem_basis": (
                "sigma*sqrt(2*tau_int/T) — the error on a correlated mean. sigma/sqrt(N) assumes "
                "every step is an independent draw, which no accepted trajectory satisfies"
            ),
        }


def ess_interval(
    series: Any,
    dt: float,
    *,
    observable: str,
    sokal_window_c: float,
    coverage_z: float,
) -> ESSInterval:
    """Return the ESS-based mean and interval of one correlated observable series.

    Args:
        series: Observable samples from an accepted trajectory, evenly spaced.
        dt: SAMPLE spacing [s] — not the solver step.  With a telemetry stride the two differ, and
            handing the physical step to this estimator understates ``tau_int`` by exactly that stride,
            which makes the interval too narrow in the dangerous direction.
        observable: Name, recorded in the interval.
        sokal_window_c: The caller's Sokal automatic-window constant.
        coverage_z: Coverage factor: 1.0 for one standard error, 1.96 for a nominal 95 % interval.
            Declared by the caller and recorded, because "the interval" means nothing on its own.

    Returns:
        The :class:`ESSInterval`.

    Raises:
        ValueError: If ``dt`` or ``coverage_z`` is not positive-finite, or the series is too short.
    """
    dt_s = _positive(dt, what="dt")
    z = _positive(coverage_z, what="coverage_z")
    name = _nonempty(observable, what="observable")
    x = _series(series, what="series")
    if x.size < 4:
        raise ValueError(f"need >= 4 samples for an ESS interval; got {x.size}")

    tau_s, n_eff, note = integrated_autocorrelation_time(x, dt_s, sokal_window_c=sokal_window_c)
    total_time = float(x.size) * dt_s
    std = float(np.std(x))
    mean = float(np.mean(x))
    sem_naive = std / math.sqrt(x.size)
    sem_ess = std * math.sqrt(2.0 * tau_s / total_time)
    widening = 1.0 if sem_naive <= 0.0 else max(1.0, sem_ess / sem_naive)
    return ESSInterval(
        observable=name,
        mean=mean,
        std=std,
        n_samples=int(x.size),
        tau_int_s=tau_s,
        n_eff=n_eff,
        sem_naive=sem_naive,
        sem_ess=sem_ess,
        widening_factor=widening,
        coverage_z=z,
        low=mean - z * sem_ess,
        high=mean + z * sem_ess,
        tau_note=note,
    )


@dataclass(frozen=True, slots=True)
class SweepPoint:
    """One parameter value, the accepted trajectory measured there, and its regime.

    The regime is part of the point, not an afterthought applied to the difference: a sweep that stores
    only ``(parameter, mean)`` has already thrown away the information that would let it refuse.
    """

    parameter_name: str
    parameter_value: float
    parameter_unit: str
    observable: str
    observable_unit: str
    series: np.ndarray
    dt_s: float
    regime: RegimeReport
    interval: ESSInterval

    @classmethod
    def measure(
        cls,
        series: Any,
        dt: float,
        *,
        parameter_name: str,
        parameter_value: float,
        parameter_unit: str,
        observable: str,
        observable_unit: str,
        amplitude_scale: float,
        regime_contract: RegimeEvidenceContract,
        coverage_z: float,
    ) -> SweepPoint:
        """Classify the trajectory's regime and compute its ESS interval in one pass.

        Args:
            series: Observable samples from an accepted trajectory.
            dt: Sample spacing [s].
            parameter_name: The swept parameter's name.
            parameter_value: Its value at this point.
            parameter_unit: Its unit, carried into any sensitivity taken from this point.
            observable: The measured observable's name.
            observable_unit: Its unit.
            amplitude_scale: The observable's declared physical scale, for the fixed-point test.
            regime_contract: Pre-declared regime evidence thresholds.
            coverage_z: Interval coverage factor.

        Returns:
            The :class:`SweepPoint`.
        """
        x = _series(series, what="series")
        report = classify_regime(
            x,
            dt,
            observable=observable,
            amplitude_scale=amplitude_scale,
            contract=regime_contract,
        )
        interval = ess_interval(
            x,
            dt,
            observable=observable,
            sokal_window_c=regime_contract.sokal_window_c,
            coverage_z=coverage_z,
        )
        return cls(
            parameter_name=_nonempty(parameter_name, what="parameter_name"),
            parameter_value=_finite(parameter_value, what="parameter_value"),
            parameter_unit=_nonempty(parameter_unit, what="parameter_unit"),
            observable=_nonempty(observable, what="observable"),
            observable_unit=_nonempty(observable_unit, what="observable_unit"),
            series=x,
            dt_s=_positive(dt, what="dt"),
            regime=report,
            interval=interval,
        )


@dataclass(frozen=True, slots=True)
class SensitivityVerdict:
    """The result of asking for a finite-difference sensitivity — a number, or a routed refusal.

    ``value`` is ``None`` unless ``kind`` is ``REPORTED``.  There is no "best effort" number: a
    sensitivity across a regime boundary is not an approximate derivative, it is not a derivative.
    """

    kind: SweepVerdictKind
    parameter_name: str
    observable: str
    unit: str
    low_value: float
    high_value: float
    low_regime: Regime
    high_regime: Regime
    value: float | None = None
    uncertainty: float | None = None
    coverage_z: float | None = None
    reason: str = ""
    expansions: tuple[ExpansionRoute, ...] = ()
    evidence_authority: str = EVIDENCE_AUTHORITY

    def __post_init__(self) -> None:
        if self.kind is SweepVerdictKind.REPORTED:
            if self.value is None or self.uncertainty is None:
                raise ValueError("a REPORTED sensitivity must carry a value and an uncertainty")
        else:
            if self.value is not None:
                raise ValueError("a refused sensitivity must not carry a value")
            if not self.expansions:
                raise ValueError(
                    "a refusal must name at least one expansion route; a refusal with no route is "
                    "how a hard problem gets dropped instead of scheduled"
                )

    @property
    def is_reported(self) -> bool:
        """Whether a number was actually issued."""
        return self.kind is SweepVerdictKind.REPORTED

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able verdict."""
        return {
            "kind": self.kind.value,
            "parameter_name": self.parameter_name,
            "observable": self.observable,
            "unit": self.unit,
            "low_value": self.low_value,
            "high_value": self.high_value,
            "low_regime": self.low_regime.value,
            "high_regime": self.high_regime.value,
            "value": self.value,
            "uncertainty": self.uncertainty,
            "coverage_z": self.coverage_z,
            "reason": self.reason,
            "expansions": [route.value for route in self.expansions],
            "evidence_authority": self.evidence_authority,
        }


def finite_difference_sensitivity(
    low: SweepPoint, high: SweepPoint, *, min_n_eff: float
) -> SensitivityVerdict:
    """Return ``d(observable)/d(parameter)`` between two swept points, or a routed refusal.

    The refusals, in the order they are checked:

    1. **Either point is ``UNDECIDED``.**  The regime is the precondition for the question; an
       undecided point cannot be differenced against anything.  Routes to a longer or better-sampled
       run, not to a bolder classifier.
    2. **The two points are in DIFFERENT regimes.**  This is the central refusal.  Between a stochastic
       steady state and a limit cycle there is no smooth ``f(p)`` to differentiate: the "observable" on
       one side is a mean the trajectory rests near and on the other a mean it passes through, and the
       finite difference is the gap between two different kinds of object.  Routes to
       ``REGIME_PARTITION`` — sweep each regime separately and locate the boundary as its own object —
       and to ``LATENT_COORDINATE``, since a coordinate in which the response is smooth across the
       boundary (an order parameter, an oscillation amplitude) is exactly what the partition is missing.
    3. **Either point has too few independent samples.**  A difference of two means is only as good as
       the means, and an ESS below the declared floor makes both meaningless.

    Uncertainty is propagated from the two ESS intervals in quadrature and divided by the parameter
    separation, so the reported uncertainty inherits the ESS widening rather than the step count.

    Args:
        low: The point at the lower parameter value.
        high: The point at the higher parameter value.
        min_n_eff: Declared minimum effective sample size at each point.

    Returns:
        The :class:`SensitivityVerdict`.

    Raises:
        TypeError: If either argument is not a :class:`SweepPoint`.
        ValueError: If the two points do not share a parameter, observable or units, or if the
            parameter separation is zero — a zero separation is a caller error, not a finding.
    """
    for point in (low, high):
        if not isinstance(point, SweepPoint):
            raise TypeError("finite_difference_sensitivity expects two SweepPoint objects")
    if low.parameter_name != high.parameter_name:
        raise ValueError("both points must sweep the same parameter")
    if low.observable != high.observable:
        raise ValueError("both points must measure the same observable")
    if low.parameter_unit != high.parameter_unit or low.observable_unit != high.observable_unit:
        raise ValueError("both points must carry the same units")
    separation = high.parameter_value - low.parameter_value
    if separation == 0.0 or not math.isfinite(separation):
        raise ValueError("the two points must have a nonzero, finite parameter separation")
    _positive(min_n_eff, what="min_n_eff")

    unit = f"{low.observable_unit}/{low.parameter_unit}"
    common = {
        "parameter_name": low.parameter_name,
        "observable": low.observable,
        "unit": unit,
        "low_value": low.parameter_value,
        "high_value": high.parameter_value,
        "low_regime": low.regime.regime,
        "high_regime": high.regime.regime,
    }

    undecided = [p for p in (low, high) if p.regime.regime is Regime.UNDECIDED]
    if undecided:
        which = " and ".join(f"p={p.parameter_value:g}" for p in undecided)
        return SensitivityVerdict(
            kind=SweepVerdictKind.REFUSED_UNDECIDED_REGIME,
            reason=(
                f"regime UNDECIDED at {which}: {undecided[0].regime.reason}. A difference cannot be "
                "taken against a point whose dynamical state is unknown"
            ),
            expansions=(ExpansionRoute.TIME_MODALITY, ExpansionRoute.OBSERVATION_MODEL_EXTENSION),
            **common,
        )

    if low.regime.regime is not high.regime.regime:
        return SensitivityVerdict(
            kind=SweepVerdictKind.REFUSED_CROSS_REGIME,
            reason=(
                f"{low.parameter_name}={low.parameter_value:g} is {low.regime.regime.value} and "
                f"{low.parameter_name}={high.parameter_value:g} is {high.regime.regime.value}. There "
                "is no smooth response function spanning these two points, so their difference is not "
                "a derivative: it is the gap between two different kinds of object, and it will change "
                "with the step size rather than converge"
            ),
            expansions=(ExpansionRoute.REGIME_PARTITION, ExpansionRoute.LATENT_COORDINATE),
            **common,
        )

    thin = [p for p in (low, high) if p.interval.n_eff < min_n_eff]
    if thin:
        which = ", ".join(f"p={p.parameter_value:g} n_eff={p.interval.n_eff:.3g}" for p in thin)
        return SensitivityVerdict(
            kind=SweepVerdictKind.REFUSED_INSUFFICIENT_ESS,
            reason=(
                f"effective sample size below the declared floor {min_n_eff:g} at {which}. The step "
                "count is not the sample size: these trajectories are correlated"
            ),
            expansions=(ExpansionRoute.TIME_MODALITY, ExpansionRoute.ADAPTIVE_DESIGN),
            **common,
        )

    value = (high.interval.mean - low.interval.mean) / separation
    combined = math.hypot(high.interval.sem_ess, low.interval.sem_ess) / abs(separation)
    coverage = low.interval.coverage_z
    return SensitivityVerdict(
        kind=SweepVerdictKind.REPORTED,
        value=value,
        uncertainty=coverage * combined,
        coverage_z=coverage,
        reason=(
            f"both points are {low.regime.regime.value}; uncertainty propagated from ESS standard "
            f"errors (widening {low.interval.widening_factor:.3g}x and "
            f"{high.interval.widening_factor:.3g}x over the naive form)"
        ),
        **common,
    )


@dataclass(frozen=True, slots=True)
class PhaseDecomposition:
    """A cycle separated into its phase-resolved profile and the scalar mean that hides it.

    ``time_average`` is what a sweep would have reported.  ``phase_bin_means`` is what the trajectory
    actually does.  ``cycle_amplitude`` — the peak-to-peak span of the phase profile — is the size of
    what the scalar threw away, and ``amplitude_over_sem`` says whether that loss is resolvable or
    merely asserted.
    """

    observable: str
    period_s: float
    n_bins: int
    n_cycles: float
    phase_bin_centres: np.ndarray
    phase_bin_means: np.ndarray
    phase_bin_sems: np.ndarray
    phase_bin_counts: np.ndarray
    time_average: float
    phase_average: float
    cycle_amplitude: float
    amplitude_over_sem: float

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able decomposition."""
        return {
            "observable": self.observable,
            "period_s": self.period_s,
            "n_bins": self.n_bins,
            "n_cycles": self.n_cycles,
            "phase_bin_centres": self.phase_bin_centres.tolist(),
            "phase_bin_means": self.phase_bin_means.tolist(),
            "phase_bin_sems": self.phase_bin_sems.tolist(),
            "phase_bin_counts": self.phase_bin_counts.tolist(),
            "time_average": self.time_average,
            "phase_average": self.phase_average,
            "cycle_amplitude": self.cycle_amplitude,
            "amplitude_over_sem": self.amplitude_over_sem,
        }


def phase_decomposition(
    series: Any, dt: float, *, observable: str, period_s: float, n_bins: int
) -> PhaseDecomposition:
    """Fold an oscillatory trajectory onto its own period; the time-average is NOT the observable.

    For a regime that :class:`~aleph.virtual_cell.regime.RegimeReport.is_oscillatory` reports true
    for, a scalar mean is the centre of an orbit — a value the system is at only in passing.  Folding
    every sample onto phase ``(t mod T) / T`` recovers the profile that mean averaged away, and the
    per-bin standard errors say whether the profile is resolved or is noise.

    ``phase_average`` — the unweighted mean of the bin means — is reported next to ``time_average``
    because they differ whenever the bins are unequally populated, which is the numerical fingerprint
    of an incomplete final cycle.

    Args:
        series: Observable samples, evenly spaced.
        dt: Sample spacing [s].
        observable: Name, recorded in the decomposition.
        period_s: The cycle period [s] — normally ``1 / f`` of the dominant fundamental from the regime
            report, so the fold uses a MEASURED period rather than an assumed one.
        n_bins: Number of phase bins.

    Returns:
        The :class:`PhaseDecomposition`.

    Raises:
        ValueError: If the record covers less than one full period, or has fewer samples than bins —
            an empty phase bin is a hole in the profile, not a zero.
    """
    dt_s = _positive(dt, what="dt")
    period = _positive(period_s, what="period_s")
    name = _nonempty(observable, what="observable")
    if isinstance(n_bins, bool) or not isinstance(n_bins, int) or n_bins < 2:
        raise ValueError("n_bins must be an int >= 2")
    x = _series(series, what="series")
    total_time = x.size * dt_s
    if total_time < period:
        raise ValueError(
            f"record covers {total_time:.4g} s of a {period:.4g} s period: a fold over less than one "
            "cycle is not a phase average"
        )
    if x.size < n_bins:
        raise ValueError(f"{x.size} samples cannot fill {n_bins} phase bins")

    times = np.arange(x.size, dtype=np.float64) * dt_s
    phase = np.mod(times, period) / period
    index = np.minimum((phase * n_bins).astype(int), n_bins - 1)
    counts = np.bincount(index, minlength=n_bins)
    if np.any(counts == 0):
        raise ValueError(
            "at least one phase bin received no sample; reduce n_bins or lengthen the record rather "
            "than filling the hole with a zero"
        )
    totals = np.bincount(index, weights=x, minlength=n_bins)
    means = totals / counts
    squares = np.bincount(index, weights=x * x, minlength=n_bins)
    variance = np.maximum(squares / counts - means * means, 0.0)
    sems = np.sqrt(variance / np.maximum(counts - 1, 1))
    amplitude = float(np.max(means) - np.min(means))
    typical_sem = float(np.mean(sems))
    return PhaseDecomposition(
        observable=name,
        period_s=period,
        n_bins=int(n_bins),
        n_cycles=float(total_time / period),
        phase_bin_centres=(np.arange(n_bins) + 0.5) / n_bins,
        phase_bin_means=means,
        phase_bin_sems=sems,
        phase_bin_counts=counts.astype(np.int64),
        time_average=float(np.mean(x)),
        phase_average=float(np.mean(means)),
        cycle_amplitude=amplitude,
        amplitude_over_sem=(amplitude / typical_sem) if typical_sem > 0.0 else float("inf"),
    )


@dataclass(frozen=True, slots=True)
class EnsembleMember:
    """One realisation's scalar observables, with the seed that makes its independence checkable.

    ``seed`` is deliberately ``int | None`` rather than a required field: a member with no recorded
    seed is a real and common defect, and it must be reportable by :func:`audit_ensemble` alongside the
    others rather than raising at construction and hiding the rest of the audit.
    """

    member_id: str
    observables: Mapping[str, float]
    seed: int | None = None
    weight: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "member_id", _nonempty(self.member_id, what="member_id"))
        values: dict[str, float] = {}
        for key, value in dict(self.observables).items():
            values[_nonempty(key, what="observable name")] = _finite(
                value, what=f"observable {key!r}"
            )
        object.__setattr__(self, "observables", dict(values))
        if self.seed is not None:
            if isinstance(self.seed, bool) or not isinstance(self.seed, int):
                raise TypeError("seed must be an int or None")
            if self.seed < 0:
                raise ValueError("seed must be nonnegative")
        if self.weight is not None:
            object.__setattr__(self, "weight", _finite(self.weight, what="weight"))


@dataclass(frozen=True, slots=True)
class EnsembleDefectRecord:
    """One detected defect, naming the members involved and what it would have corrupted."""

    defect: EnsembleDefect
    detail: str
    member_ids: tuple[str, ...]
    expansions: tuple[ExpansionRoute, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able defect record."""
        return {
            "defect": self.defect.value,
            "detail": self.detail,
            "member_ids": list(self.member_ids),
            "expansions": [route.value for route in self.expansions],
        }


@dataclass(frozen=True, slots=True)
class EnsembleSummary:
    """The audit of an ensemble, plus its means only if the audit is clean.

    ``apparent_n`` is how many members were supplied and ``independent_n`` how many distinct seeds they
    carry.  When they differ, the ensemble's effective size is the smaller one, and quoting the larger
    is the seed-provenance failure stated as a number.
    """

    apparent_n: int
    independent_n: int
    observables: tuple[str, ...]
    defects: tuple[EnsembleDefectRecord, ...]
    means: Mapping[str, float] = field(default_factory=dict)
    weighted: bool = False
    evidence_authority: str = EVIDENCE_AUTHORITY

    @property
    def is_clean(self) -> bool:
        """Whether the ensemble may be averaged at all."""
        return not self.defects

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able summary."""
        return {
            "apparent_n": self.apparent_n,
            "independent_n": self.independent_n,
            "observables": list(self.observables),
            "defects": [record.as_dict() for record in self.defects],
            "means": dict(self.means),
            "weighted": self.weighted,
            "is_clean": self.is_clean,
            "evidence_authority": self.evidence_authority,
        }


def audit_ensemble(
    members: Sequence[EnsembleMember],
    *,
    required_observables: Sequence[str] | None = None,
    weight_tolerance: float,
) -> tuple[EnsembleDefectRecord, ...]:
    """Return every detected defect in an ensemble — all of them, not the first.

    The audit is exhaustive on purpose.  A caller that fixes one defect and reruns learns about the
    next one a cycle later; a caller handed the whole list fixes the ensemble once.

    Args:
        members: The realisations.
        required_observables: Observables every member must carry.  When ``None`` (the default) the
            requirement is the UNION of what the members actually have — which is what makes "member 3
            is missing ``tension``" detectable at all.  Passing an explicit list additionally catches
            an observable that no member has.
        weight_tolerance: How far a weight sum may sit from 1 before it is a defect.

    Returns:
        The defect records, possibly empty.

    Raises:
        ValueError: If ``members`` is empty or ``weight_tolerance`` is not finite and nonnegative.
        TypeError: If any entry is not an :class:`EnsembleMember`.
    """
    entries = list(members)
    if not entries:
        raise ValueError("an ensemble audit needs at least one member")
    for entry in entries:
        if not isinstance(entry, EnsembleMember):
            raise TypeError("members must be EnsembleMember objects")
    tolerance = _finite(weight_tolerance, what="weight_tolerance")
    if tolerance < 0.0:
        raise ValueError("weight_tolerance must be nonnegative")

    defects: list[EnsembleDefectRecord] = []

    if required_observables is None:
        required = sorted({name for entry in entries for name in entry.observables})
    else:
        required = sorted(
            {_nonempty(name, what="required observable") for name in required_observables}
        )

    for name in required:
        missing = tuple(entry.member_id for entry in entries if name not in entry.observables)
        if missing:
            defects.append(
                EnsembleDefectRecord(
                    defect=EnsembleDefect.MISSING_OBSERVABLE,
                    detail=(
                        f"{len(missing)} of {len(entries)} members lack {name!r}. Averaging over the "
                        "members that have it would silently change the estimator's population "
                        "between observables, so two numbers in the same table would be means over "
                        "different ensembles"
                    ),
                    member_ids=missing,
                    expansions=(ExpansionRoute.ENSEMBLE, ExpansionRoute.PROVENANCE_REPAIR),
                )
            )

    seedless = tuple(entry.member_id for entry in entries if entry.seed is None)
    if seedless:
        defects.append(
            EnsembleDefectRecord(
                defect=EnsembleDefect.MISSING_SEED,
                detail=(
                    f"{len(seedless)} members carry no seed. Unknown independence is not the same as "
                    "independence, and the ensemble spread assumes the latter"
                ),
                member_ids=seedless,
                expansions=(ExpansionRoute.PROVENANCE_REPAIR,),
            )
        )
    seen: dict[int, list[str]] = {}
    for entry in entries:
        if entry.seed is not None:
            seen.setdefault(entry.seed, []).append(entry.member_id)
    for seed, ids in sorted(seen.items()):
        if len(ids) > 1:
            defects.append(
                EnsembleDefectRecord(
                    defect=EnsembleDefect.DUPLICATE_SEED,
                    detail=(
                        f"seed {seed} appears in {len(ids)} members. They ran the same stream, so they "
                        "are not independent draws: the spread is too small and the apparent sample "
                        f"size is inflated by {len(ids) - 1}"
                    ),
                    member_ids=tuple(ids),
                    expansions=(ExpansionRoute.PROVENANCE_REPAIR, ExpansionRoute.ENSEMBLE),
                )
            )

    weighted = [entry for entry in entries if entry.weight is not None]
    if weighted and len(weighted) != len(entries):
        defects.append(
            EnsembleDefectRecord(
                defect=EnsembleDefect.PARTIAL_WEIGHTS,
                detail=(
                    f"{len(weighted)} of {len(entries)} members carry a mixture weight. There is no "
                    "defensible default for the rest: an implicit 1 changes the mixture and an "
                    "implicit 0 deletes members"
                ),
                member_ids=tuple(entry.member_id for entry in entries if entry.weight is None),
                expansions=(ExpansionRoute.PROVENANCE_REPAIR,),
            )
        )
    elif weighted:
        negative = tuple(entry.member_id for entry in weighted if entry.weight < 0.0)
        if negative:
            defects.append(
                EnsembleDefectRecord(
                    defect=EnsembleDefect.NEGATIVE_WEIGHT,
                    detail="a mixture weight below zero is not a proportion of anything",
                    member_ids=negative,
                    expansions=(ExpansionRoute.PROVENANCE_REPAIR,),
                )
            )
        total = float(sum(entry.weight for entry in weighted))
        if abs(total - 1.0) > tolerance:
            defects.append(
                EnsembleDefectRecord(
                    defect=EnsembleDefect.WEIGHTS_NOT_NORMALIZED,
                    detail=(
                        f"mixture weights sum to {total:.17g}, off by {abs(total - 1.0):.3g} against a "
                        f"declared tolerance {tolerance:g}. Renormalising here would hide whichever "
                        "member is wrong"
                    ),
                    member_ids=tuple(entry.member_id for entry in weighted),
                    expansions=(ExpansionRoute.PROVENANCE_REPAIR,),
                )
            )

    return tuple(defects)


def summarize_ensemble(
    members: Sequence[EnsembleMember],
    *,
    required_observables: Sequence[str] | None = None,
    weight_tolerance: float,
) -> EnsembleSummary:
    """Audit an ensemble and, only if it is clean, return its means.

    A defective ensemble returns a summary with ``means`` empty and ``is_clean`` false.  It does not
    return partial means: a mean over the members that happened to carry an observable is precisely the
    silent averaging this function exists to prevent.

    Args:
        members: The realisations.
        required_observables: See :func:`audit_ensemble`.
        weight_tolerance: See :func:`audit_ensemble`.

    Returns:
        The :class:`EnsembleSummary`.
    """
    entries = list(members)
    defects = audit_ensemble(
        entries,
        required_observables=required_observables,
        weight_tolerance=weight_tolerance,
    )
    names = tuple(sorted({name for entry in entries for name in entry.observables}))
    independent = len({entry.seed for entry in entries if entry.seed is not None})
    weighted = all(entry.weight is not None for entry in entries) and bool(entries)
    if defects:
        return EnsembleSummary(
            apparent_n=len(entries),
            independent_n=independent,
            observables=names,
            defects=defects,
            weighted=weighted,
        )
    means: dict[str, float] = {}
    for name in names:
        values = np.array([entry.observables[name] for entry in entries], dtype=np.float64)
        if weighted:
            weights = np.array([float(entry.weight) for entry in entries], dtype=np.float64)
            means[name] = float(np.dot(weights, values))
        else:
            means[name] = float(np.mean(values))
    return EnsembleSummary(
        apparent_n=len(entries),
        independent_n=independent,
        observables=names,
        defects=(),
        means=means,
        weighted=weighted,
    )


@dataclass(frozen=True, slots=True)
class ClosureGapReport:
    """The measured moment-closure cost of representing a distribution by its mean.

    ``gap = E[f(X)] - f(E[X])``.  **A nonzero gap is expected and is not a defect in the reference
    dynamics.**  It is the price of the closure: the moment you choose to propagate ``E[X]`` instead of
    the samples, every nonlinear observable acquires this error, and its size is a property of the
    representation you adopted, not evidence that the underlying model is wrong.  Reading it the other
    way round — "the point model is buggy" — misattributes a known, quantifiable representation cost to
    the physics.

    What the tolerance is for: it says whether the mean-field closure is good enough **for this
    observable, at this spread**, and a breach routes to a richer representation (a second moment, a
    mixture, the full ensemble), never to a change in the dynamics.
    """

    observable: str
    ensemble_mean_of_f: float
    f_of_ensemble_mean: float
    gap: float
    relative_gap: float
    tolerance: float
    within_tolerance: bool
    n_samples: int
    weighted: bool
    expansions: tuple[ExpansionRoute, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able closure report."""
        return {
            "observable": self.observable,
            "ensemble_mean_of_f": self.ensemble_mean_of_f,
            "f_of_ensemble_mean": self.f_of_ensemble_mean,
            "gap": self.gap,
            "relative_gap": self.relative_gap,
            "tolerance": self.tolerance,
            "within_tolerance": self.within_tolerance,
            "n_samples": self.n_samples,
            "weighted": self.weighted,
            "expansions": [route.value for route in self.expansions],
            "interpretation": (
                "moment-closure cost of the distributional representation; a nonzero gap is expected "
                "and is NOT a defect in the reference dynamics"
            ),
        }


def moment_closure_gap(
    samples: Any,
    transform: Callable[[np.ndarray], np.ndarray],
    *,
    observable: str,
    tolerance: float,
    weights: Any | None = None,
) -> ClosureGapReport:
    """Measure ``E[f(X)] - f(E[X])`` for a nonlinear observable over an ensemble.

    Verified against closed form: for ``f(x) = x**2`` this returns exactly the ensemble's (population)
    variance, ``E[X^2] - (E[X])^2``, which is an identity and not an approximation — so the measurement
    is checked against truth rather than against itself.

    Args:
        samples: The ensemble's values of the underlying quantity ``X``.
        transform: The nonlinear observable ``f``, applied elementwise to an ndarray.  It must be a
            pure elementwise map: applying it to the mean must be meaningful, which is the whole
            comparison being made.
        observable: Name of ``f(X)``, recorded in the report.
        tolerance: Declared limit on ``|gap| / max(|E[f(X)]|, |f(E[X])|)``.  Its own quantity with its
            own budget, separate from any tolerance on the dynamics.
        weights: Optional mixture weights.  When given they must be nonnegative and sum to 1 within
            ``1e-9``; the gap is then taken over the weighted measure.

    Returns:
        The :class:`ClosureGapReport`.

    Raises:
        ValueError: If the samples are empty or non-finite, ``transform`` does not return a matching
            shape, or the weights are malformed.
    """
    name = _nonempty(observable, what="observable")
    limit = _finite(tolerance, what="tolerance")
    if limit < 0.0:
        raise ValueError("tolerance must be nonnegative")
    x = _series(samples, what="samples")
    if weights is None:
        w = np.full(x.size, 1.0 / x.size, dtype=np.float64)
        weighted = False
    else:
        w = _series(weights, what="weights")
        if w.size != x.size:
            raise ValueError("weights must have one entry per sample")
        if np.any(w < 0.0):
            raise ValueError("mixture weights must be nonnegative")
        if abs(float(np.sum(w)) - 1.0) > 1e-9:
            raise ValueError(f"mixture weights must sum to 1; got {float(np.sum(w))!r}")
        weighted = True

    transformed = np.asarray(transform(x), dtype=np.float64).ravel()
    if transformed.shape != x.shape:
        raise ValueError("transform must map the sample array elementwise to the same shape")
    if not np.all(np.isfinite(transformed)):
        raise ValueError("transform produced a non-finite value")

    mean_of_f = float(np.dot(w, transformed))
    mean_x = float(np.dot(w, x))
    f_of_mean = float(np.asarray(transform(np.array([mean_x], dtype=np.float64))).ravel()[0])
    gap = mean_of_f - f_of_mean
    scale = max(abs(mean_of_f), abs(f_of_mean))
    relative = abs(gap) / scale if scale > 0.0 else abs(gap)
    within = relative <= limit
    return ClosureGapReport(
        observable=name,
        ensemble_mean_of_f=mean_of_f,
        f_of_ensemble_mean=f_of_mean,
        gap=gap,
        relative_gap=relative,
        tolerance=limit,
        within_tolerance=within,
        n_samples=int(x.size),
        weighted=weighted,
        expansions=(() if within else (ExpansionRoute.ENSEMBLE, ExpansionRoute.LATENT_COORDINATE)),
    )
