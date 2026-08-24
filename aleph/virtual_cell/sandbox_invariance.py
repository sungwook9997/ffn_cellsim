"""Is the regime a property of the system, or of the integrator that produced it?

WHY THIS MODULE EXISTS.  :mod:`aleph.virtual_cell.regime` will label a trajectory a limit cycle, a
quasicycle or a stochastic steady state.  Nothing in that classifier knows where the trajectory came
from, so it will label a numerical artefact exactly as confidently as it labels physics — and an
explicit integrator run past its stability step manufactures a *beautiful* sustained oscillation.  This
project has repeatedly mistaken a solver artefact for a result, and the master plan's falsifier for the
whole sandbox is one sentence: **"is the oscillation/distribution an artefact of the numerical
method?"**  The failure criterion is that the classified REGIME depends on ``dt``, on the inner-solve
tolerance, on the iteration budget or on the seed.  This module is the battery that asks.

**A battery that only ever passes has not been shown to be able to fail.**  So the module ships the
positive control with the machinery: :func:`damped_oscillator_series` integrated by explicit Euler is a
system whose apparent regime is *genuinely* a discretisation artefact — below the stability step
``dt* = 2 zeta / omega`` the noise-driven oscillator sits at rest, above it the scheme's own
anti-damping sustains and grows an oscillation that is not in the differential equation at all.  The
battery must FLAG that, and the flagged bracket must contain the analytic ``dt*``.  The same system
under :func:`damped_oscillator_series` with the semi-implicit scheme, at the same step, does not flag —
which is what makes the verdict a statement about the method rather than about the system.

FOUR THINGS THIS DOES THAT A CONVERGENCE STUDY DOES NOT.

1. **It compares at matched PHYSICAL time, and refuses when it cannot.**  Halving ``dt`` at a fixed step
   count halves the experiment.  That is not an abstract hazard: `STATE.md` (c) 16 records a
   blebbistatin dose response in which all six points were retracted because each covered 0.60 s of a
   required 7.14 s.  A ladder whose rungs differ in physical duration by more than the declared
   tolerance returns :attr:`InvarianceVerdict.REFUSED_UNMATCHED_DURATION` and no comparison is made; a
   rung covering fewer than the declared number of its own measured correlation times is refused for
   the same reason, because matched wall duration is not matched *statistical* duration.
2. **It compares with ESS-based uncertainty, and reports the naive counterfactual next to it.**  Two
   runs at different ``dt`` have different autocorrelation *in units of steps*, so the finer rung has
   more samples that are individually less informative.  A naive ``sigma/sqrt(N)`` therefore makes the
   finer rung look more precise purely because its samples are more correlated, and a ladder judged
   that way flags invariance failures that are not there.  Every comparison here uses
   :func:`aleph.virtual_cell.sweep.ess_interval`; ``z_naive`` is computed alongside ``z_ess`` and
   reported, never acted on, so the size of the trap is visible in the artifact.
3. **It separates "the regime is seed-invariant" from "the trajectory is seed-invariant".**  The first
   is required.  The second is FALSE for every stochastic system and is asserted nowhere in this
   module: :class:`SeedInvarianceReport` has no trajectory-equality field to set, reports the measured
   across-seed trajectory divergence as evidence that the seeds really were independent, and refuses
   any ensemble below the declared seed count — a single-seed claim is not a claim.
4. **It routes a statistic that tracks the solver rather than the physics.**  The adversarial case this
   project actually hit: a quantity that moves by orders of magnitude with the iteration budget while
   the trajectory underneath is unchanged.  A statistic declared ``PHYSICAL`` whose dynamic range across
   the knob ladder exceeds the declared bound *while another declared-physical statistic on the same
   runs is invariant* is flagged :attr:`InvarianceVerdict.FLAGGED_RESIDUAL_TRACKING` — distinct from an
   ordinary invariance failure, because the repair is different: the statistic is not an observable of
   the physics, and no amount of solver work will make it one.

AUTHORITY.  This is an **unverified observer**.  It grants no evidence rung; a verdict here is a
statement about a time series and a knob, never about physical validity, and may not be promoted to
native evidence.  No Warp, no device, no GPU: numpy and the CPU only.

NOT A CONVERGENCE TEST.  Nothing here scores a residual falling to zero, and "the force converged" is
not a success criterion in any code path.  A ladder of runs that all sustain a large steady fluctuation
is exactly as invariant as a ladder that all sit at rest.

Sanity Gate:
    * dimensional: ``dt_s`` and every duration are [s]; the knob carries the axis' own units (seconds
      for ``TIME_STEP``, the observable's units for ``TOLERANCE``, a dimensionless count for
      ``ITERATION_BUDGET``, a label for ``SEED``).  ``z_ess`` and ``z_naive`` are dimensionless; the
      statistic's units are the caller's and are recorded on the :class:`StatisticSpec`.  ``dt_s`` is
      the SAMPLE spacing, not the solver step — with a telemetry stride the two differ and every
      correlation time is understated by exactly that stride.
    * boundary: a ladder shorter than ``min_ladder_rungs`` is refused (two points cannot show a trend);
      a rung whose series is non-finite is recorded as diverged and never handed to the classifier; a
      rung below ``min_n_eff`` or ``min_correlation_windows`` is refused rather than compared; a seed
      ensemble below ``min_seeds`` or containing a duplicate seed is refused.
    * conservation/invariant: physical duration is ``n_samples * dt_s`` and is checked, not assumed;
      the Brownian path of :func:`shared_brownian_increments` conserves its increments exactly under
      coarsening (a coarse increment is the SUM of the fine increments it spans, so every rung of a
      refinement ladder is driven by the same realisation of the same Wiener process); ``sem_ess >=
      sem_naive`` always, and their ratio is the reported widening factor.
    * numerical: comparisons are made against a caller-DECLARED reference rung, not against the mean of
      the ladder — an average over rungs contaminated by the artefact is contaminated by the artefact.
      Trajectory differences are evaluated on a common physical time grid at the coarser spacing by
      nearest-sample selection; no interpolation, because interpolating one rung onto another's grid
      low-pass filters it and would understate exactly the difference being measured.
    * sign-sense: ``z`` is an absolute discrepancy and is never signed away; a dynamic range is
      ``max|v| / min|v| >= 1`` by construction, and a statistic passing through zero reports an infinite
      range rather than a rescued finite one.
    * measurement-protocol: every threshold lives in a caller-supplied :class:`InvarianceContract` with
      **no defaults**, declared before the run and copied into the report, and the Sokal window constant
      used for the ESS intervals is taken from the caller's :class:`RegimeEvidenceContract` so that the
      correlation time behind the interval and the correlation time behind the regime label can never
      be two different numbers.
    * resolution limit: a ladder resolves the knob value at which a regime changes only to the BRACKET
      between two adjacent rungs.  :attr:`InvarianceReport.regime_change_bracket` returns that bracket
      and never a point estimate.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import numpy as np

from aleph.virtual_cell.contracts import SandboxExperimentCard
from aleph.virtual_cell.regime import (
    Regime,
    RegimeEvidenceContract,
    RegimeReport,
    classify_regime,
)
from aleph.virtual_cell.sweep import ESSInterval, ExpansionRoute, ess_interval

__all__ = [
    "EVIDENCE_AUTHORITY",
    "ImplicitSolveRun",
    "InvarianceAxis",
    "InvarianceContract",
    "InvarianceReport",
    "InvarianceVerdict",
    "LadderRun",
    "RungEvidence",
    "RungStatistic",
    "SeedInvarianceReport",
    "StatisticRole",
    "StatisticSpec",
    "StatisticVerdict",
    "centered_square",
    "damped_oscillator_series",
    "euler_growth_per_step",
    "euler_stability_step_s",
    "identity",
    "implicit_relaxation_run",
    "invariance_experiment_card",
    "knob_invariance",
    "matched_time_difference",
    "ou_refinement_ladder",
    "ou_series",
    "seed_invariance",
    "shared_brownian_increments",
]

#: What this layer is allowed to claim.  An invariance verdict is never native evidence.
EVIDENCE_AUTHORITY = "unverified-observer"


class InvarianceAxis(StrEnum):
    """The numerical knob a ladder varies.  One axis per ladder; a ladder that moves two knobs at once
    cannot attribute what it finds to either."""

    TIME_STEP = "time-step"
    TOLERANCE = "tolerance"
    ITERATION_BUDGET = "iteration-budget"
    SEED = "seed"


class StatisticRole(StrEnum):
    """What the caller claims a statistic IS, declared before the ladder runs.

    * ``PHYSICAL`` — an observable of the modelled system.  It is asserted invariant under every axis in
      this module, and a failure is a finding.
    * ``SOLVER_DIAGNOSTIC`` — a quantity belonging to the solver, not the system: a residual norm, an
      iteration count, a line-search length.  It is EXPECTED to move with the knob, is never flagged,
      and — the point of the distinction — may never be quoted as a physical result.  Declaring a
      residual ``PHYSICAL`` is the adversarial case this battery exists to catch, so the honest
      declaration has to be available for the dishonest one to be a mistake rather than the only option.
    """

    PHYSICAL = "physical"
    SOLVER_DIAGNOSTIC = "solver-diagnostic"


class InvarianceVerdict(StrEnum):
    """The battery's answer.  Refusals are first-class and are expected to be common.

    * ``INVARIANT`` — every rung agreed, under the declared contract, on evidence sufficient to say so.
    * ``NOT_ASSERTED_DIAGNOSTIC`` — a ``SOLVER_DIAGNOSTIC`` statistic: measured and reported, never
      asserted invariant, never quotable as physics.
    * ``REFUSED_UNMATCHED_DURATION`` — the rungs do not cover the same physical time.  No comparison was
      made; a shorter run is a different experiment, not a finer one.
    * ``REFUSED_INSUFFICIENT_EVIDENCE`` — too few rungs, or a rung with too little independent sampling
      to support any statement.
    * ``REFUSED_UNDECIDED_REGIME`` — at least one rung's regime was ``UNDECIDED`` and the decided rungs
      agreed.  Invariance is NOT established: an absent label is not a matching label.
    * ``FLAGGED_STATISTIC_ARTEFACT`` — a declared-physical statistic moved with the knob by more than
      its own ESS interval allows.
    * ``FLAGGED_RESIDUAL_TRACKING`` — a declared-physical statistic moved by orders of magnitude while
      another declared-physical statistic on the same runs did not.  It is tracking the solver.
    * ``FLAGGED_REGIME_ARTEFACT`` — two rungs carry two different DECIDED regime labels.  The falsifier
      fired: the regime is a property of the discretisation.
    * ``FLAGGED_NUMERICAL_DIVERGENCE`` — a rung produced a non-finite trajectory.  The strongest form of
      knob dependence, and the one a classifier can never see because it cannot be handed the series.
    """

    INVARIANT = "invariant"
    NOT_ASSERTED_DIAGNOSTIC = "not-asserted-diagnostic"
    REFUSED_UNMATCHED_DURATION = "refused-unmatched-duration"
    REFUSED_INSUFFICIENT_EVIDENCE = "refused-insufficient-evidence"
    REFUSED_UNDECIDED_REGIME = "refused-undecided-regime"
    FLAGGED_STATISTIC_ARTEFACT = "flagged-statistic-artefact"
    FLAGGED_RESIDUAL_TRACKING = "flagged-residual-tracking"
    FLAGGED_REGIME_ARTEFACT = "flagged-regime-artefact"
    FLAGGED_NUMERICAL_DIVERGENCE = "flagged-numerical-divergence"


#: Severity order used to aggregate per-statistic and regime verdicts into one.  A FLAG outranks a
#: REFUSAL on purpose: a partial refusal beside a positive finding is still a positive finding.
_SEVERITY: tuple[InvarianceVerdict, ...] = (
    InvarianceVerdict.INVARIANT,
    InvarianceVerdict.NOT_ASSERTED_DIAGNOSTIC,
    InvarianceVerdict.REFUSED_UNDECIDED_REGIME,
    InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE,
    InvarianceVerdict.REFUSED_UNMATCHED_DURATION,
    InvarianceVerdict.FLAGGED_STATISTIC_ARTEFACT,
    InvarianceVerdict.FLAGGED_RESIDUAL_TRACKING,
    InvarianceVerdict.FLAGGED_REGIME_ARTEFACT,
    InvarianceVerdict.FLAGGED_NUMERICAL_DIVERGENCE,
)

_EXPANSIONS: dict[InvarianceVerdict, tuple[ExpansionRoute, ...]] = {
    InvarianceVerdict.INVARIANT: (),
    InvarianceVerdict.NOT_ASSERTED_DIAGNOSTIC: (),
    InvarianceVerdict.REFUSED_UNMATCHED_DURATION: (ExpansionRoute.TIME_MODALITY,),
    InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE: (ExpansionRoute.TIME_MODALITY,),
    InvarianceVerdict.REFUSED_UNDECIDED_REGIME: (
        ExpansionRoute.TIME_MODALITY,
        ExpansionRoute.LATENT_COORDINATE,
    ),
    InvarianceVerdict.FLAGGED_STATISTIC_ARTEFACT: (
        ExpansionRoute.SOLVER_DESIGN,
        ExpansionRoute.ADAPTIVE_DESIGN,
    ),
    InvarianceVerdict.FLAGGED_RESIDUAL_TRACKING: (
        ExpansionRoute.OBSERVATION_MODEL_EXTENSION,
        ExpansionRoute.SOLVER_DESIGN,
    ),
    InvarianceVerdict.FLAGGED_REGIME_ARTEFACT: (
        ExpansionRoute.SOLVER_DESIGN,
        ExpansionRoute.ADAPTIVE_DESIGN,
    ),
    InvarianceVerdict.FLAGGED_NUMERICAL_DIVERGENCE: (ExpansionRoute.SOLVER_DESIGN,),
}


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


def _finite(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite; got {value!r}")
    return number


def _positive_int(value: object, *, what: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{what} must be an int")
    if value < minimum:
        raise ValueError(f"{what} must be >= {minimum}; got {value!r}")
    return value


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


@dataclass(frozen=True, slots=True)
class InvarianceContract:
    """The caller's pre-declared invariance criteria.  **No field has a default, by design.**

    A threshold invented inside the battery is a threshold chosen after seeing the ladder, and the
    whole point of the battery is that this project has done exactly that before.  The contract is
    copied verbatim into every report so a verdict can be re-judged without the raw runs.

    Args:
        min_ladder_rungs: Fewest rungs a ladder may have.  Two points cannot show a trend and cannot
            bracket a threshold.
        duration_relative_tolerance: Maximum ``(max - min) / max`` of the rungs' physical durations.
            Above this the comparison is REFUSED: rungs that ran for different physical times are
            different experiments.  The value is a protocol tolerance, not a physics threshold — the
            durations should match by construction, and this catches the case where they do not.
        min_correlation_windows: Each rung's record must cover at least this many of its OWN measured
            ``tau_int``.  Matched wall duration is not matched statistical duration.
        min_n_eff: Each rung's compared statistic must rest on at least this many independent samples.
        statistic_z_max: Maximum ``|v_i - v_ref| / sqrt(sem_ess_i^2 + sem_ess_ref^2)`` before a
            declared-physical statistic is flagged.  A consistency criterion in units of the
            statistic's own stated error, which is what "beyond its own stated error" has to mean.
        residual_tracking_dynamic_range_min: ``max|v| / min|v|`` across the ladder at or above which a
            flagged physical statistic is re-routed as residual-tracking rather than as an ordinary
            invariance failure — provided at least one other declared-physical statistic on the same
            runs is invariant.  Orders of magnitude, not percent: a physical observable does not move
            by this factor when only the solver knob moved.
        min_seeds: Fewest seeds a seed-robustness claim may rest on.  A single-seed claim is not a
            claim, and two seeds cannot measure a spread.
        min_seed_regime_fraction: Fraction of seeds that must carry the modal DECIDED regime before the
            regime is called seed-invariant.
        ess_calibration_ratio_min: Lower bound of the band in which the across-seed spread of the mean
            divided by the mean within-run ESS standard error must fall.  This is the ESS interval
            checking itself against an independent estimate of the same quantity: for an ergodic
            system the two measure the same thing and their ratio is 1.
        ess_calibration_ratio_max: Upper bound of that band.
        coverage_z: Coverage factor recorded on every interval.  1.0 for one standard error, 1.96 for a
            nominal 95 % interval.  Declared, because "the interval" means nothing on its own.
    """

    min_ladder_rungs: int
    duration_relative_tolerance: float
    min_correlation_windows: float
    min_n_eff: float
    statistic_z_max: float
    residual_tracking_dynamic_range_min: float
    min_seeds: int
    min_seed_regime_fraction: float
    ess_calibration_ratio_min: float
    ess_calibration_ratio_max: float
    coverage_z: float

    def __post_init__(self) -> None:
        _positive_int(self.min_ladder_rungs, what="min_ladder_rungs", minimum=3)
        _positive(self.duration_relative_tolerance, what="duration_relative_tolerance")
        if self.duration_relative_tolerance >= 1.0:
            raise ValueError(
                "duration_relative_tolerance must be below 1: a tolerance that admits a rung of "
                "zero duration admits the protocol error it exists to catch"
            )
        _positive(self.min_correlation_windows, what="min_correlation_windows")
        _positive(self.min_n_eff, what="min_n_eff")
        _positive(self.statistic_z_max, what="statistic_z_max")
        if self.residual_tracking_dynamic_range_min <= 1.0:
            raise ValueError(
                "residual_tracking_dynamic_range_min must exceed 1; a range of 1 is no variation"
            )
        _positive_int(self.min_seeds, what="min_seeds", minimum=3)
        fraction = _positive(self.min_seed_regime_fraction, what="min_seed_regime_fraction")
        if fraction > 1.0:
            raise ValueError("min_seed_regime_fraction must lie in (0, 1]")
        _positive(self.ess_calibration_ratio_min, what="ess_calibration_ratio_min")
        _positive(self.ess_calibration_ratio_max, what="ess_calibration_ratio_max")
        if self.ess_calibration_ratio_min >= 1.0 or self.ess_calibration_ratio_max <= 1.0:
            raise ValueError(
                "the ESS calibration band must straddle 1: for an ergodic system the across-seed "
                "spread and the within-run ESS error estimate the same quantity, so a band that "
                "excludes 1 asserts the estimator is wrong before it is measured"
            )
        _positive(self.coverage_z, what="coverage_z")

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able contract, recorded next to every verdict it produced."""
        return {
            "min_ladder_rungs": int(self.min_ladder_rungs),
            "duration_relative_tolerance": float(self.duration_relative_tolerance),
            "min_correlation_windows": float(self.min_correlation_windows),
            "min_n_eff": float(self.min_n_eff),
            "statistic_z_max": float(self.statistic_z_max),
            "residual_tracking_dynamic_range_min": float(self.residual_tracking_dynamic_range_min),
            "min_seeds": int(self.min_seeds),
            "min_seed_regime_fraction": float(self.min_seed_regime_fraction),
            "ess_calibration_ratio_min": float(self.ess_calibration_ratio_min),
            "ess_calibration_ratio_max": float(self.ess_calibration_ratio_max),
            "coverage_z": float(self.coverage_z),
        }


@dataclass(frozen=True, slots=True, eq=False)
class LadderRun:
    """One rung: a recorded observable, its SAMPLE spacing, and the knob value that produced it.

    The series is stored read-only.  A rung is never mutated after construction, because a battery whose
    inputs can change between the regime call and the statistic call is measuring two different things.

    Args:
        knob_value: The varied knob, in its axis' own units.  Must be distinct across the ladder.
        dt_s: Sample spacing [s] — not the solver step, if telemetry was strided.
        series: Observable samples, evenly spaced.  Non-finite values are ALLOWED here and are recorded
            as a divergence; refusing to construct the rung would hide the most important failure a
            step-size ladder can produce.
        seed: The RNG stream, when there was one.  ``None`` means deterministic.
        channels: Further same-length series recorded on the SAME steps — a per-step residual norm, an
            iteration count, a second observable.  They exist so that a physical anchor and a solver
            diagnostic can be compared on one ladder: the residual-tracking route needs to know that the
            trajectory did not change while the suspect quantity did, and it can only know that if both
            came from the same run.  The regime is always classified on ``series`` alone.
    """

    knob_value: float
    dt_s: float
    series: np.ndarray
    seed: int | None = None
    channels: dict[str, np.ndarray] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "knob_value", _finite(self.knob_value, what="knob_value"))
        object.__setattr__(self, "dt_s", _positive(self.dt_s, what="dt_s"))
        array = np.asarray(self.series, dtype=np.float64).ravel()
        if array.size < 4:
            raise ValueError(f"a rung needs >= 4 samples; got {array.size}")
        array = array.copy()
        array.flags.writeable = False
        object.__setattr__(self, "series", array)
        if self.seed is not None and (
            isinstance(self.seed, bool) or not isinstance(self.seed, int)
        ):
            raise TypeError("seed must be an int or None")
        extra: dict[str, np.ndarray] = {}
        for key, values in (self.channels or {}).items():
            channel = np.asarray(values, dtype=np.float64).ravel()
            if channel.size != array.size:
                raise ValueError(
                    f"channel {key!r} has {channel.size} samples against the observable's "
                    f"{array.size}: a channel is recorded on the same steps, and one that is not "
                    "sits on a different time axis than the correlation time used for its interval"
                )
            channel = channel.copy()
            channel.flags.writeable = False
            extra[_nonempty(key, what="channel name")] = channel
        object.__setattr__(self, "channels", extra)

    def channel_series(self, name: str | None) -> np.ndarray:
        """Return the primary observable when ``name`` is ``None``, else the named channel.

        Raises:
            KeyError: If the channel was not recorded on this rung.
        """
        if name is None:
            return self.series
        channels = self.channels or {}
        if name not in channels:
            raise KeyError(f"rung at knob {self.knob_value!r} has no channel {name!r}")
        return channels[name]

    @property
    def n_samples(self) -> int:
        """Number of recorded samples."""
        return int(self.series.size)

    @property
    def physical_duration_s(self) -> float:
        """Physical time the rung covers [s].  ``n_samples * dt_s`` — computed, never assumed."""
        return float(self.series.size) * self.dt_s

    @property
    def is_finite(self) -> bool:
        """Whether every sample is finite.  ``False`` is a divergence, not a malformed input."""
        return bool(np.all(np.isfinite(self.series)))

    @property
    def is_analysable(self) -> bool:
        """Whether the rung can carry a second moment at all.

        Finite samples are not sufficient.  An explicit scheme run past its stability step reaches
        ``1e243`` while every sample is still a valid double — and then the sum of squares overflows,
        so the variance, the autocorrelation and the periodogram all return ``nan``.  Every statistic
        this battery and the classifier compute is a second moment or built on one, so a rung whose
        squares are not representable has diverged in the only sense that matters here, and is
        recorded as diverged rather than fed to an estimator that will answer ``nan``.
        """
        if not self.is_finite:
            return False
        with np.errstate(over="ignore", invalid="ignore"):
            return bool(np.isfinite(np.dot(self.series, self.series)))


def identity(series: np.ndarray) -> np.ndarray:
    """Return the series unchanged: the per-sample transform whose time average is the mean."""
    return series


def centered_square(series: np.ndarray) -> np.ndarray:
    """Return ``(x - mean(x))**2``: the per-sample transform whose time average is the variance.

    The centring uses each rung's OWN mean.  Using a shared mean would fold a difference in the means
    into the variance comparison and report one failure as two.
    """
    return (series - float(np.mean(series))) ** 2


@dataclass(frozen=True, slots=True, eq=False)
class StatisticSpec:
    """A scalar to be compared across the ladder, as the time average of a per-sample transform.

    Expressing the statistic as a transform rather than as a reduction is what makes the ESS interval
    available: :func:`aleph.virtual_cell.sweep.ess_interval` needs the transformed SERIES to measure
    the correlation time of the quantity actually being averaged.  A statistic handed over as a single
    number has already thrown that away, and its error bar can then only be invented.

    Args:
        name: Recorded on every verdict.
        transform: Per-sample map from the observable series to the series whose mean is the statistic.
        role: What the caller claims this quantity is.  See :class:`StatisticRole`.
        units: The statistic's units, recorded so the number is quotable.
        channel: Which recorded series the transform reads.  ``None`` — the default — is the rung's
            primary observable; a name selects a :attr:`LadderRun.channels` entry.
    """

    name: str
    transform: Callable[[np.ndarray], np.ndarray]
    role: StatisticRole
    units: str
    channel: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _nonempty(self.name, what="name"))
        object.__setattr__(self, "units", _nonempty(self.units, what="units"))
        if not callable(self.transform):
            raise TypeError("transform must be callable")
        if not isinstance(self.role, StatisticRole):
            raise TypeError("role must be a StatisticRole")
        if self.channel is not None:
            object.__setattr__(self, "channel", _nonempty(self.channel, what="channel"))


@dataclass(frozen=True, slots=True)
class RungStatistic:
    """One statistic measured on one rung, with the ESS interval and the naive counterfactual.

    ``sem_naive`` and ``widening_factor`` are carried but never acted on.  They are the size of the trap
    this module exists to avoid: the amount by which a ``sigma/sqrt(N)`` interval understates the error
    on a correlated mean, which grows as ``dt`` falls and therefore biases a refinement ladder
    systematically toward flagging its own finest rung.
    """

    name: str
    units: str
    role: StatisticRole
    value: float
    sem_ess: float
    sem_naive: float
    widening_factor: float
    tau_int_s: float
    n_eff: float
    low: float
    high: float
    tau_note: str

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able record."""
        return {
            "name": self.name,
            "units": self.units,
            "role": self.role.value,
            "value": self.value,
            "sem_ess": self.sem_ess,
            "sem_naive": self.sem_naive,
            "widening_factor": self.widening_factor,
            "tau_int_s": self.tau_int_s,
            "n_eff": self.n_eff,
            "low": self.low,
            "high": self.high,
            "tau_note": self.tau_note,
        }


@dataclass(frozen=True, slots=True)
class RungEvidence:
    """Everything measured on one rung, whatever the ladder concluded.

    A refusal with no numbers cannot say which measurement was missing, so the evidence is reported for
    refused and diverged rungs too — with ``regime`` ``None`` where the classifier was never run,
    because a rung that was not classified must not be indistinguishable from one classified
    ``UNDECIDED``.
    """

    knob_value: float
    dt_s: float
    n_samples: int
    physical_duration_s: float
    diverged: bool
    regime: Regime | None
    regime_reason: str
    tau_int_s: float
    correlation_windows: float
    n_eff: float
    trajectory_difference_over_std: float
    statistics: tuple[RungStatistic, ...]

    def statistic(self, name: str) -> RungStatistic:
        """Return the named statistic on this rung.

        Raises:
            KeyError: If the statistic was not measured here.
        """
        for record in self.statistics:
            if record.name == name:
                return record
        raise KeyError(f"rung at knob {self.knob_value!r} has no statistic {name!r}")

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able record."""
        return {
            "knob_value": self.knob_value,
            "dt_s": self.dt_s,
            "n_samples": self.n_samples,
            "physical_duration_s": self.physical_duration_s,
            "diverged": self.diverged,
            "regime": None if self.regime is None else self.regime.value,
            "regime_reason": self.regime_reason,
            "tau_int_s": self.tau_int_s,
            "correlation_windows": self.correlation_windows,
            "n_eff": self.n_eff,
            "trajectory_difference_over_std": self.trajectory_difference_over_std,
            "statistics": [record.as_dict() for record in self.statistics],
        }


@dataclass(frozen=True, slots=True)
class StatisticVerdict:
    """Whether one statistic survived the ladder, with both the ESS and the naive comparison.

    ``max_z_naive`` is the verdict the same ladder would have received from a ``sigma/sqrt(N)`` interval.
    It is reported for one reason: when it exceeds ``max_z_ess`` and the contract's limit, the naive
    comparison would have manufactured an invariance failure out of nothing but correlation.
    """

    name: str
    units: str
    role: StatisticRole
    verdict: InvarianceVerdict
    reason: str
    reference_knob: float
    reference_value: float
    max_z_ess: float
    max_z_naive: float
    max_z_knob: float
    dynamic_range: float
    values_by_knob: tuple[tuple[float, float], ...]
    expansions: tuple[ExpansionRoute, ...]

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able record."""
        return {
            "name": self.name,
            "units": self.units,
            "role": self.role.value,
            "verdict": self.verdict.value,
            "reason": self.reason,
            "reference_knob": self.reference_knob,
            "reference_value": self.reference_value,
            "max_z_ess": self.max_z_ess,
            "max_z_naive": self.max_z_naive,
            "max_z_knob": self.max_z_knob,
            "dynamic_range": self.dynamic_range,
            "values_by_knob": [list(pair) for pair in self.values_by_knob],
            "expansions": [route.value for route in self.expansions],
        }


@dataclass(frozen=True, slots=True)
class InvarianceReport:
    """One ladder's verdict: the regime, every statistic, and where each failure routes."""

    axis: InvarianceAxis
    observable: str
    verdict: InvarianceVerdict
    reason: str
    regime_verdict: InvarianceVerdict
    regime_reason: str
    reference_knob: float
    rungs: tuple[RungEvidence, ...]
    statistic_verdicts: tuple[StatisticVerdict, ...]
    contract: dict[str, Any]
    regime_contract: dict[str, Any]
    expansions: tuple[ExpansionRoute, ...]
    evidence_authority: str = EVIDENCE_AUTHORITY

    @property
    def regime_by_knob(self) -> tuple[tuple[float, str], ...]:
        """``(knob, label)`` for every rung, ``"not-classified"`` where the classifier never ran."""
        return tuple(
            (rung.knob_value, "not-classified" if rung.regime is None else rung.regime.value)
            for rung in self.rungs
        )

    @property
    def decided_regimes(self) -> frozenset[Regime]:
        """The distinct DECIDED regime labels on the ladder.  Two entries is the falsifier firing."""
        return frozenset(
            rung.regime
            for rung in self.rungs
            if rung.regime is not None and rung.regime is not Regime.UNDECIDED
        )

    @property
    def regime_change_bracket(self) -> tuple[float, float] | None:
        """The adjacent knob pair between which the DECIDED label first changes, or ``None``.

        A bracket, never a point estimate: a ladder localises a threshold only to the interval between
        two rungs it actually ran, and quoting the crossing as a single number would claim a resolution
        the experiment does not have.
        """
        ordered = sorted(
            (
                rung
                for rung in self.rungs
                if rung.regime is not None and rung.regime is not Regime.UNDECIDED
            ),
            key=lambda rung: rung.knob_value,
        )
        for low, high in zip(ordered, ordered[1:], strict=False):
            if low.regime is not high.regime:
                return (low.knob_value, high.knob_value)
        return None

    def statistic(self, name: str) -> StatisticVerdict:
        """Return the named statistic's verdict.

        Raises:
            KeyError: If no such statistic was declared.
        """
        for verdict in self.statistic_verdicts:
            if verdict.name == name:
                return verdict
        raise KeyError(f"no statistic {name!r} on this ladder")

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able report."""
        return {
            "axis": self.axis.value,
            "observable": self.observable,
            "verdict": self.verdict.value,
            "reason": self.reason,
            "regime_verdict": self.regime_verdict.value,
            "regime_reason": self.regime_reason,
            "reference_knob": self.reference_knob,
            "regime_by_knob": [list(pair) for pair in self.regime_by_knob],
            "regime_change_bracket": (
                None if self.regime_change_bracket is None else list(self.regime_change_bracket)
            ),
            "rungs": [rung.as_dict() for rung in self.rungs],
            "statistic_verdicts": [verdict.as_dict() for verdict in self.statistic_verdicts],
            "contract": dict(self.contract),
            "regime_contract": dict(self.regime_contract),
            "expansions": [route.value for route in self.expansions],
            "evidence_authority": self.evidence_authority,
        }


def matched_time_difference(reference: LadderRun, run: LadderRun) -> float:
    """Return the RMS difference between two runs at MATCHED PHYSICAL TIMES, over the reference's std.

    Two rungs of a refinement ladder have different sample counts, so a sample-by-sample difference
    compares physical time ``k dt_a`` with physical time ``k dt_b`` — which is the same protocol error
    as comparing a fixed step count across a ``dt`` ladder, one level down.  Both series are therefore
    read on a common grid at the COARSER spacing, over the shorter of the two durations, by
    nearest-sample selection.

    Nearest-sample and not interpolation: interpolating the fine rung onto the coarse grid low-pass
    filters it, and the high-frequency content is precisely where two integrations of the same system
    differ most.  Filtering it away would understate the very quantity being measured.

    Args:
        reference: The rung the comparison is against; its standard deviation sets the scale.
        run: The rung being compared.

    Returns:
        RMS difference divided by the reference's standard deviation, or ``nan`` when either rung is
        not :attr:`LadderRun.is_analysable`, the overlap is shorter than 4 samples, or the reference is
        constant — this is itself a squared quantity, so a rung whose squares overflow has no
        difference to report either.  For two
        independent realisations of the same stationary process this is ``sqrt(2)``; for a trajectory
        compared with itself it is 0.
    """
    if not reference.is_analysable or not run.is_analysable:
        return float("nan")
    duration = min(reference.physical_duration_s, run.physical_duration_s)
    step = max(reference.dt_s, run.dt_s)
    count = int(duration // step)
    if count < 4:
        return float("nan")
    times = np.arange(count, dtype=np.float64) * step
    left = reference.series[
        np.minimum((times / reference.dt_s).astype(int), reference.n_samples - 1)
    ]
    right = run.series[np.minimum((times / run.dt_s).astype(int), run.n_samples - 1)]
    scale = float(np.std(left))
    if scale <= 0.0:
        return float("nan")
    return float(np.sqrt(np.mean((left - right) ** 2)) / scale)


def _rung_statistic(
    run: LadderRun, spec: StatisticSpec, *, sokal_window_c: float, coverage_z: float
) -> RungStatistic:
    """Measure one statistic on one rung, with its ESS interval and the naive counterfactual."""
    transformed = np.asarray(
        spec.transform(run.channel_series(spec.channel)), dtype=np.float64
    ).ravel()
    if transformed.size != run.n_samples:
        raise ValueError(
            f"statistic {spec.name!r} returned {transformed.size} samples from {run.n_samples}: a "
            "per-sample transform must preserve the sampling grid, or the correlation time behind "
            "its interval is measured on a different time axis than the series"
        )
    with np.errstate(over="ignore", invalid="ignore"):
        representable = bool(np.all(np.isfinite(transformed))) and bool(
            np.isfinite(np.dot(transformed, transformed))
        )
    if not representable:
        # Same reason as :attr:`LadderRun.is_analysable`: an interval is a second moment, so a series
        # whose squares overflow has no interval, and computing one anyway would return nan through a
        # trail of overflow warnings rather than saying so.
        return RungStatistic(
            name=spec.name,
            units=spec.units,
            role=spec.role,
            value=float("nan"),
            sem_ess=float("nan"),
            sem_naive=float("nan"),
            widening_factor=float("nan"),
            tau_int_s=float("nan"),
            n_eff=0.0,
            low=float("nan"),
            high=float("nan"),
            tau_note=(
                "the transformed series is non-finite, or its sum of squares is; no interval "
                "was estimated"
            ),
        )
    interval: ESSInterval = ess_interval(
        transformed,
        run.dt_s,
        observable=spec.name,
        sokal_window_c=sokal_window_c,
        coverage_z=coverage_z,
    )
    return RungStatistic(
        name=spec.name,
        units=spec.units,
        role=spec.role,
        value=interval.mean,
        sem_ess=interval.sem_ess,
        sem_naive=interval.sem_naive,
        widening_factor=interval.widening_factor,
        tau_int_s=interval.tau_int_s,
        n_eff=interval.n_eff,
        low=interval.low,
        high=interval.high,
        tau_note=interval.tau_note,
    )


def _discrepancy(value: float, sem: float, reference: float, reference_sem: float) -> float:
    """Return ``|v - v_ref| / sqrt(sem^2 + sem_ref^2)``, ``inf`` when the errors are zero but the
    values differ, and 0 when neither varies."""
    if not (math.isfinite(value) and math.isfinite(reference)):
        return float("inf")
    spread = math.hypot(sem, reference_sem)
    difference = abs(value - reference)
    if spread <= 0.0:
        return 0.0 if difference == 0.0 else float("inf")
    return difference / spread


def _dynamic_range(values: Sequence[float]) -> float:
    """Return ``max|v| / min|v|`` over the ladder; ``inf`` if the statistic passes through zero."""
    magnitudes = [abs(value) for value in values if math.isfinite(value)]
    if not magnitudes:
        return float("nan")
    smallest, largest = min(magnitudes), max(magnitudes)
    if smallest <= 0.0:
        return float("inf") if largest > 0.0 else 1.0
    return largest / smallest


def _validate_ladder(
    runs: Sequence[LadderRun], *, axis: InvarianceAxis, reference_knob: float
) -> tuple[LadderRun, ...]:
    """Return the ladder ordered by knob, after the structural checks that precede any measurement."""
    if not isinstance(axis, InvarianceAxis):
        raise TypeError("axis must be an InvarianceAxis")
    ordered = tuple(sorted(runs, key=lambda run: run.knob_value))
    if not ordered:
        raise ValueError("a ladder needs at least one rung")
    for run in ordered:
        if not isinstance(run, LadderRun):
            raise TypeError("every rung must be a LadderRun")
    knobs = [run.knob_value for run in ordered]
    if len(set(knobs)) != len(knobs):
        raise ValueError(
            "ladder knob values must be distinct; two rungs at the same knob measure repeatability, "
            "not invariance, and averaging them here would hide that"
        )
    if not any(math.isclose(knob, reference_knob, rel_tol=1e-12, abs_tol=0.0) for knob in knobs):
        raise ValueError(
            f"reference_knob {reference_knob!r} is not a rung of this ladder; the reference must be "
            "declared and present, never inferred from the data"
        )
    return ordered


def knob_invariance(
    runs: Sequence[LadderRun],
    *,
    axis: InvarianceAxis,
    observable: str,
    amplitude_scale: float,
    reference_knob: float,
    statistics: Sequence[StatisticSpec],
    regime_contract: RegimeEvidenceContract,
    contract: InvarianceContract,
) -> InvarianceReport:
    """Run the invariance battery over one knob ladder and return the verdict with its evidence.

    The order of the checks is fixed and documented here, because a battery whose order is implicit is
    a battery nobody can audit:

    1. **Structural.**  Enough rungs, distinct knobs, the declared reference present.  A ladder failing
       these raises rather than returning a verdict: it is a malformed experiment, not a negative one.
    2. **Matched physical duration.**  Rungs whose durations differ by more than the declared tolerance
       are REFUSED before anything is measured on them.  This is check 2 and not check 6 on purpose —
       a comparison that is not entitled to be made must not produce a number that could be quoted.
    3. **Per-rung evidence.**  Divergence, regime, correlation time, ESS interval per statistic, and the
       trajectory difference from the reference rung at matched physical times.
    4. **Regime invariance.**  Two different decided labels is the falsifier firing.  A mix of decided
       and undecided is a REFUSAL, not a pass: an absent label is not a matching label.
    5. **Statistic invariance**, in ESS units, against the declared reference rung — and, for a flagged
       declared-physical statistic that moved by orders of magnitude while another declared-physical
       statistic did not, the residual-tracking route instead.

    Args:
        runs: The rungs.  Order is irrelevant; they are sorted by knob.
        axis: Which knob the ladder varies.
        observable: Name of the recorded observable, passed through to the classifier.
        amplitude_scale: The observable's declared physical scale, pre-registered, in its own units.
            Only the classifier's fixed-point test uses it.
        reference_knob: The rung every comparison is made against.  Declared by the caller — for a
            refinement ladder the finest step, for a tolerance ladder the tightest — and never inferred,
            because a reference chosen after seeing the ladder is a reference chosen to pass.
        statistics: The scalars to compare.  At least one; declare the physical anchor as ``PHYSICAL``
            and solver internals as ``SOLVER_DIAGNOSTIC``.
        regime_contract: The classifier's pre-declared thresholds.  Its ``sokal_window_c`` is also used
            for the ESS intervals, so the correlation time behind an interval and the one behind a
            regime label are the same number.
        contract: The invariance criteria, pre-declared.

    Returns:
        The :class:`InvarianceReport`, always carrying full per-rung evidence whatever the verdict.

    Raises:
        TypeError: If ``axis``, a rung, a statistic or either contract has the wrong type.
        ValueError: If the ladder is structurally malformed, or no statistic was declared.
    """
    if not isinstance(contract, InvarianceContract):
        raise TypeError("contract must be an InvarianceContract declared before the run")
    if not isinstance(regime_contract, RegimeEvidenceContract):
        raise TypeError("regime_contract must be a RegimeEvidenceContract declared before the run")
    name = _nonempty(observable, what="observable")
    scale = _positive(amplitude_scale, what="amplitude_scale")
    specs = tuple(statistics)
    if not specs:
        raise ValueError(
            "declare at least one statistic; an invariance claim about nothing is empty"
        )
    for spec in specs:
        if not isinstance(spec, StatisticSpec):
            raise TypeError("every statistic must be a StatisticSpec")
    if len({spec.name for spec in specs}) != len(specs):
        raise ValueError("statistic names must be distinct")
    ordered = _validate_ladder(runs, axis=axis, reference_knob=reference_knob)
    contract_dict = contract.to_dict()
    regime_dict = regime_contract.to_dict()

    def _report(
        verdict: InvarianceVerdict,
        reason: str,
        *,
        regime_verdict: InvarianceVerdict,
        regime_reason: str,
        rungs: tuple[RungEvidence, ...],
        statistic_verdicts: tuple[StatisticVerdict, ...],
    ) -> InvarianceReport:
        return InvarianceReport(
            axis=axis,
            observable=name,
            verdict=verdict,
            reason=reason,
            regime_verdict=regime_verdict,
            regime_reason=regime_reason,
            reference_knob=float(reference_knob),
            rungs=rungs,
            statistic_verdicts=statistic_verdicts,
            contract=contract_dict,
            regime_contract=regime_dict,
            expansions=_EXPANSIONS[verdict],
        )

    # 1. Structural: too short a ladder cannot show a trend or bracket a threshold.
    if len(ordered) < contract.min_ladder_rungs:
        return _report(
            InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE,
            (
                f"{len(ordered)} rungs against a declared minimum of {contract.min_ladder_rungs}: "
                "a two-point ladder cannot distinguish an invariant statistic from a trend"
            ),
            regime_verdict=InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE,
            regime_reason="ladder too short to compare",
            rungs=(),
            statistic_verdicts=(),
        )

    # 2. Matched physical duration, BEFORE any measurement.  A comparison that is not entitled to be
    #    made must not emit a number that could later be quoted without its refusal.
    durations = [run.physical_duration_s for run in ordered]
    longest = max(durations)
    duration_spread = (longest - min(durations)) / longest
    if duration_spread > contract.duration_relative_tolerance:
        return _report(
            InvarianceVerdict.REFUSED_UNMATCHED_DURATION,
            (
                f"physical durations span {min(durations):.6g} s to {longest:.6g} s, a relative "
                f"spread of {duration_spread:.3g} against a declared "
                f"{contract.duration_relative_tolerance:g}. Rungs that cover different physical "
                "times are different experiments: a finer dt at the same STEP COUNT is a shorter "
                "run, which is the protocol error that retracted all six blebbistatin points "
                "(STATE.md (c) 16, 0.60 s against a required 7.14 s)"
            ),
            regime_verdict=InvarianceVerdict.REFUSED_UNMATCHED_DURATION,
            regime_reason="durations not matched; no rung was classified",
            rungs=(),
            statistic_verdicts=(),
        )

    reference_run = next(
        run
        for run in ordered
        if math.isclose(run.knob_value, reference_knob, rel_tol=1e-12, abs_tol=0.0)
    )

    # 3. Per-rung evidence.
    evidence: list[RungEvidence] = []
    for run in ordered:
        records = tuple(
            _rung_statistic(
                run,
                spec,
                sokal_window_c=regime_contract.sokal_window_c,
                coverage_z=contract.coverage_z,
            )
            for spec in specs
        )
        if not run.is_analysable:
            evidence.append(
                RungEvidence(
                    knob_value=run.knob_value,
                    dt_s=run.dt_s,
                    n_samples=run.n_samples,
                    physical_duration_s=run.physical_duration_s,
                    diverged=True,
                    regime=None,
                    regime_reason=(
                        "the trajectory is non-finite, or finite with a sum of squares that is not "
                        "— either way every second moment the classifier computes returns nan. It "
                        "was never handed to the classifier, because a diverged run has no regime "
                        "to read"
                    ),
                    tau_int_s=float("nan"),
                    correlation_windows=float("nan"),
                    n_eff=0.0,
                    trajectory_difference_over_std=float("nan"),
                    statistics=records,
                )
            )
            continue
        report: RegimeReport = classify_regime(
            run.series,
            run.dt_s,
            observable=name,
            amplitude_scale=scale,
            contract=regime_contract,
        )
        tau_s = report.statistics.tau_int_s
        windows = (
            report.statistics.analysed_time_s / tau_s
            if math.isfinite(tau_s) and tau_s > 0.0
            else float("nan")
        )
        evidence.append(
            RungEvidence(
                knob_value=run.knob_value,
                dt_s=run.dt_s,
                n_samples=run.n_samples,
                physical_duration_s=run.physical_duration_s,
                diverged=False,
                regime=report.regime,
                regime_reason=report.reason,
                tau_int_s=tau_s,
                correlation_windows=windows,
                n_eff=report.statistics.n_eff,
                trajectory_difference_over_std=matched_time_difference(reference_run, run),
                statistics=records,
            )
        )
    rungs = tuple(evidence)

    # 4. Regime.
    regime_verdict, regime_reason = _regime_verdict(rungs, contract=contract)

    # 5. Statistics.
    reference_evidence = next(
        rung
        for rung in rungs
        if math.isclose(rung.knob_value, reference_knob, rel_tol=1e-12, abs_tol=0.0)
    )
    statistic_verdicts = _statistic_verdicts(
        rungs, specs, reference=reference_evidence, contract=contract
    )

    considered = [regime_verdict] + [
        verdict.verdict
        for verdict in statistic_verdicts
        if verdict.verdict is not InvarianceVerdict.NOT_ASSERTED_DIAGNOSTIC
    ]
    overall = max(considered, key=_SEVERITY.index)
    if overall is regime_verdict and overall is not InvarianceVerdict.INVARIANT:
        reason = regime_reason
    elif overall is InvarianceVerdict.INVARIANT:
        reason = (
            f"every rung of the {axis.value} ladder agreed on the regime and on every declared "
            f"physical statistic within {contract.statistic_z_max:g} ESS standard errors, at a "
            f"matched physical duration of {longest:.6g} s"
        )
    else:
        driver = next(verdict for verdict in statistic_verdicts if verdict.verdict is overall)
        reason = f"statistic {driver.name!r}: {driver.reason}"

    return _report(
        overall,
        reason,
        regime_verdict=regime_verdict,
        regime_reason=regime_reason,
        rungs=rungs,
        statistic_verdicts=statistic_verdicts,
    )


def _regime_verdict(
    rungs: Sequence[RungEvidence], *, contract: InvarianceContract
) -> tuple[InvarianceVerdict, str]:
    """Return the regime half of the verdict: divergence, then evidence, then label agreement."""
    diverged = [rung for rung in rungs if rung.diverged]
    if diverged:
        knobs = ", ".join(f"{rung.knob_value:g}" for rung in diverged)
        return (
            InvarianceVerdict.FLAGGED_NUMERICAL_DIVERGENCE,
            (
                f"rung(s) at knob {knobs} produced a non-finite trajectory. A regime that exists at "
                "one knob and diverges at another is a property of the discretisation"
            ),
        )
    starved = [
        rung
        for rung in rungs
        if not math.isfinite(rung.n_eff)
        or rung.n_eff < contract.min_n_eff
        or not math.isfinite(rung.correlation_windows)
        or rung.correlation_windows < contract.min_correlation_windows
    ]
    labels = {rung.regime for rung in rungs if rung.regime is not None}
    decided = labels - {Regime.UNDECIDED}
    if len(decided) >= 2:
        listing = ", ".join(
            f"{rung.knob_value:g} -> {rung.regime.value}"
            for rung in sorted(rungs, key=lambda item: item.knob_value)
            if rung.regime is not None
        )
        return (
            InvarianceVerdict.FLAGGED_REGIME_ARTEFACT,
            (
                f"the ladder carries {len(decided)} different DECIDED regime labels at matched "
                f"physical duration: {listing}. The regime is a property of the numerical method, "
                "not of the system"
            ),
        )
    if starved:
        knobs = ", ".join(f"{rung.knob_value:g}" for rung in starved)
        return (
            InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE,
            (
                f"rung(s) at knob {knobs} rest on fewer than {contract.min_n_eff:g} independent "
                f"samples or cover fewer than {contract.min_correlation_windows:g} of their own "
                "correlation times. Matched wall duration is not matched statistical duration"
            ),
        )
    if Regime.UNDECIDED in labels:
        knobs = ", ".join(
            f"{rung.knob_value:g}" for rung in rungs if rung.regime is Regime.UNDECIDED
        )
        return (
            InvarianceVerdict.REFUSED_UNDECIDED_REGIME,
            (
                f"rung(s) at knob {knobs} returned UNDECIDED while the rest agreed. Invariance is "
                "NOT established: an absent label is not a matching label"
            ),
        )
    if not decided:
        return (
            InvarianceVerdict.REFUSED_UNDECIDED_REGIME,
            "no rung produced a decided regime; there is nothing to compare",
        )
    label = next(iter(decided))
    return (
        InvarianceVerdict.INVARIANT,
        f"every rung classified {label.value} at matched physical duration",
    )


def _statistic_verdicts(
    rungs: Sequence[RungEvidence],
    specs: Sequence[StatisticSpec],
    *,
    reference: RungEvidence,
    contract: InvarianceContract,
) -> tuple[StatisticVerdict, ...]:
    """Return one verdict per declared statistic, resolving residual-tracking in a second pass."""
    provisional: list[StatisticVerdict] = []
    for spec in specs:
        reference_record = reference.statistic(spec.name)
        values = [(rung.knob_value, rung.statistic(spec.name).value) for rung in rungs]
        max_z_ess, max_z_naive, max_z_knob = 0.0, 0.0, reference.knob_value
        for rung in rungs:
            record = rung.statistic(spec.name)
            z_ess = _discrepancy(
                record.value, record.sem_ess, reference_record.value, reference_record.sem_ess
            )
            z_naive = _discrepancy(
                record.value, record.sem_naive, reference_record.value, reference_record.sem_naive
            )
            if z_ess > max_z_ess:
                max_z_ess, max_z_knob = z_ess, rung.knob_value
            max_z_naive = max(max_z_naive, z_naive)
        span = _dynamic_range([value for _, value in values])

        if spec.role is StatisticRole.SOLVER_DIAGNOSTIC:
            verdict = InvarianceVerdict.NOT_ASSERTED_DIAGNOSTIC
            reason = (
                f"declared a solver diagnostic: measured across the ladder with a dynamic range of "
                f"{span:.4g} and NOT asserted invariant. A solver diagnostic is not an observable of "
                "the system and may not be quoted as a physical result"
            )
        elif max_z_ess <= contract.statistic_z_max:
            if span >= contract.residual_tracking_dynamic_range_min:  # nan compares False
                # An interval wide enough to reconcile values spanning orders of magnitude has not
                # measured the statistic — it has only failed to resolve it.  Reporting that as
                # invariance would let a rung whose error bar exploded certify the ladder it broke,
                # which is how a diverging run gets quoted as agreeing with a converged one.
                verdict = InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE
                reason = (
                    f"values span a dynamic range of {span:.4g} across the ladder — at or above the "
                    f"declared {contract.residual_tracking_dynamic_range_min:g} — yet the largest "
                    f"discrepancy is only {max_z_ess:.4g} ESS standard errors. The interval is too "
                    "wide to resolve the statistic, which is a refusal and not an invariance"
                )
            else:
                verdict = InvarianceVerdict.INVARIANT
                reason = (
                    f"largest discrepancy {max_z_ess:.4g} ESS standard errors at knob "
                    f"{max_z_knob:g}, within the declared {contract.statistic_z_max:g} "
                    f"(the same comparison on naive sigma/sqrt(N) errors would have read "
                    f"{max_z_naive:.4g})"
                )
        else:
            verdict = InvarianceVerdict.FLAGGED_STATISTIC_ARTEFACT
            reason = (
                f"moved {max_z_ess:.4g} ESS standard errors at knob {max_z_knob:g}, beyond the "
                f"declared {contract.statistic_z_max:g}, with a dynamic range of {span:.4g} across "
                "the ladder"
            )
        provisional.append(
            StatisticVerdict(
                name=spec.name,
                units=spec.units,
                role=spec.role,
                verdict=verdict,
                reason=reason,
                reference_knob=reference.knob_value,
                reference_value=reference_record.value,
                max_z_ess=max_z_ess,
                max_z_naive=max_z_naive,
                max_z_knob=max_z_knob,
                dynamic_range=span,
                values_by_knob=tuple(values),
                expansions=_EXPANSIONS[verdict],
            )
        )

    # Second pass: a declared-physical statistic that moved by ORDERS OF MAGNITUDE while another
    # declared-physical statistic on the SAME runs did not is not an invariance failure of the physics
    # — it is a quantity that was never physical.  Without an invariant anchor the two cases are
    # indistinguishable and the ordinary flag stands, which is the honest answer.
    anchored = any(
        verdict.role is StatisticRole.PHYSICAL and verdict.verdict is InvarianceVerdict.INVARIANT
        for verdict in provisional
    )
    if not anchored:
        return tuple(provisional)
    resolved: list[StatisticVerdict] = []
    for verdict in provisional:
        if (
            verdict.verdict is InvarianceVerdict.FLAGGED_STATISTIC_ARTEFACT
            and math.isfinite(verdict.dynamic_range)
            and verdict.dynamic_range >= contract.residual_tracking_dynamic_range_min
        ):
            anchor = next(
                other.name
                for other in provisional
                if other.role is StatisticRole.PHYSICAL
                and other.verdict is InvarianceVerdict.INVARIANT
            )
            resolved.append(
                StatisticVerdict(
                    name=verdict.name,
                    units=verdict.units,
                    role=verdict.role,
                    verdict=InvarianceVerdict.FLAGGED_RESIDUAL_TRACKING,
                    reason=(
                        f"declared physical, but its dynamic range across the ladder is "
                        f"{verdict.dynamic_range:.4g} (>= the declared "
                        f"{contract.residual_tracking_dynamic_range_min:g}) while {anchor!r}, "
                        "measured on the SAME runs, is invariant. The trajectory did not change; "
                        "this quantity is tracking the solver, and no amount of solver work makes "
                        "it an observable of the system"
                    ),
                    reference_knob=verdict.reference_knob,
                    reference_value=verdict.reference_value,
                    max_z_ess=verdict.max_z_ess,
                    max_z_naive=verdict.max_z_naive,
                    max_z_knob=verdict.max_z_knob,
                    dynamic_range=verdict.dynamic_range,
                    values_by_knob=verdict.values_by_knob,
                    expansions=_EXPANSIONS[InvarianceVerdict.FLAGGED_RESIDUAL_TRACKING],
                )
            )
        else:
            resolved.append(verdict)
    return tuple(resolved)


@dataclass(frozen=True, slots=True)
class SeedInvarianceReport:
    """Whether the REGIME is seed-invariant — and never whether the trajectory is.

    There is deliberately no field on this class asserting trajectory equality, because for a
    stochastic system that assertion is false and a report that can express it will eventually be made
    to.  ``mean_trajectory_difference_over_std`` is reported for the opposite purpose: it is EVIDENCE
    that the seeds were independent, and for a stationary process with independent streams it sits near
    ``sqrt(2)``.  A value near zero would mean the ensemble was one run counted many times.

    ``ess_calibration_ratio`` is the ESS interval auditing itself.  The across-seed standard deviation
    of the per-seed statistic and the mean within-run ESS standard error are two independent estimates
    of the same quantity, so their ratio is 1 for a correctly calibrated interval.
    ``naive_calibration_ratio`` is the same comparison against ``sigma/sqrt(N)`` — the factor by which a
    naive interval understates the truth, measured rather than argued.
    """

    observable: str
    statistic_name: str
    verdict: InvarianceVerdict
    reason: str
    n_seeds: int
    seeds: tuple[int, ...]
    regime_counts: tuple[tuple[str, int], ...]
    modal_regime: str
    modal_fraction: float
    values: tuple[float, ...]
    value_mean: float
    across_seed_std: float
    across_seed_sem: float
    mean_sem_ess: float
    mean_sem_naive: float
    ess_calibration_ratio: float
    naive_calibration_ratio: float
    mean_trajectory_difference_over_std: float
    contract: dict[str, Any]
    expansions: tuple[ExpansionRoute, ...]
    evidence_authority: str = EVIDENCE_AUTHORITY

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able report."""
        return {
            "observable": self.observable,
            "statistic_name": self.statistic_name,
            "verdict": self.verdict.value,
            "reason": self.reason,
            "n_seeds": self.n_seeds,
            "seeds": list(self.seeds),
            "regime_counts": [list(pair) for pair in self.regime_counts],
            "modal_regime": self.modal_regime,
            "modal_fraction": self.modal_fraction,
            "values": list(self.values),
            "value_mean": self.value_mean,
            "across_seed_std": self.across_seed_std,
            "across_seed_sem": self.across_seed_sem,
            "mean_sem_ess": self.mean_sem_ess,
            "mean_sem_naive": self.mean_sem_naive,
            "ess_calibration_ratio": self.ess_calibration_ratio,
            "naive_calibration_ratio": self.naive_calibration_ratio,
            "mean_trajectory_difference_over_std": self.mean_trajectory_difference_over_std,
            "trajectory_invariance": (
                "NOT asserted and not assertable: two seeds of a stochastic system produce different "
                "trajectories by construction. Only the regime and the statistic are asserted"
            ),
            "contract": dict(self.contract),
            "expansions": [route.value for route in self.expansions],
            "evidence_authority": self.evidence_authority,
        }


def seed_invariance(
    runs: Sequence[LadderRun],
    *,
    observable: str,
    amplitude_scale: float,
    statistic: StatisticSpec,
    regime_contract: RegimeEvidenceContract,
    contract: InvarianceContract,
) -> SeedInvarianceReport:
    """Return whether the REGIME is seed-invariant, with the spread the claim rests on.

    Every run must share a ``dt`` and a physical duration — the seed is the only thing varying, and a
    seed ensemble whose members ran for different times measures the duration, not the seed.

    Args:
        runs: One rung per seed.  Every rung must carry a distinct integer ``seed``.
        observable: Name of the recorded observable.
        amplitude_scale: The observable's declared physical scale, pre-registered.
        statistic: The scalar whose across-seed spread is reported.
        regime_contract: The classifier's pre-declared thresholds.
        contract: The invariance criteria, pre-declared.

    Returns:
        The :class:`SeedInvarianceReport`.  Never asserts trajectory invariance.

    Raises:
        TypeError: If a contract or the statistic has the wrong type.
        ValueError: If a run carries no seed, two runs share a seed, or the ensemble is not a
            constant-``dt``, constant-duration ensemble.
    """
    if not isinstance(contract, InvarianceContract):
        raise TypeError("contract must be an InvarianceContract declared before the run")
    if not isinstance(regime_contract, RegimeEvidenceContract):
        raise TypeError("regime_contract must be a RegimeEvidenceContract declared before the run")
    if not isinstance(statistic, StatisticSpec):
        raise TypeError("statistic must be a StatisticSpec")
    name = _nonempty(observable, what="observable")
    scale = _positive(amplitude_scale, what="amplitude_scale")
    members = tuple(runs)
    if any(run.seed is None for run in members):
        raise ValueError(
            "every member of a seed ensemble must carry its seed; a member whose stream is unknown "
            "cannot be asserted independent, which is not the same as being independent"
        )
    seeds = tuple(int(run.seed) for run in members if run.seed is not None)
    if len(set(seeds)) != len(seeds):
        raise ValueError(
            "duplicate seed in the ensemble: two members ran the same stream, so the spread is too "
            "small and the apparent sample size is inflated by exactly the duplication"
        )
    if len({run.dt_s for run in members}) > 1 or len({run.n_samples for run in members}) > 1:
        raise ValueError(
            "a seed ensemble must hold dt and sample count fixed; otherwise the spread measures the "
            "duration as well as the seed"
        )
    contract_dict = contract.to_dict()

    def _report(
        verdict: InvarianceVerdict,
        reason: str,
        **fields: Any,
    ) -> SeedInvarianceReport:
        base: dict[str, Any] = {
            "observable": name,
            "statistic_name": statistic.name,
            "verdict": verdict,
            "reason": reason,
            "n_seeds": len(members),
            "seeds": seeds,
            "regime_counts": (),
            "modal_regime": Regime.UNDECIDED.value,
            "modal_fraction": 0.0,
            "values": (),
            "value_mean": float("nan"),
            "across_seed_std": float("nan"),
            "across_seed_sem": float("nan"),
            "mean_sem_ess": float("nan"),
            "mean_sem_naive": float("nan"),
            "ess_calibration_ratio": float("nan"),
            "naive_calibration_ratio": float("nan"),
            "mean_trajectory_difference_over_std": float("nan"),
            "contract": contract_dict,
            "expansions": _EXPANSIONS[verdict],
        }
        base.update(fields)
        return SeedInvarianceReport(**base)

    if len(members) < contract.min_seeds:
        return _report(
            InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE,
            (
                f"{len(members)} seeds against a declared minimum of {contract.min_seeds}: a "
                "single-seed claim is not a claim, and too few seeds cannot measure a spread"
            ),
        )

    labels: list[str] = []
    records: list[RungStatistic] = []
    differences: list[float] = []
    for run in members:
        record = _rung_statistic(
            run,
            statistic,
            sokal_window_c=regime_contract.sokal_window_c,
            coverage_z=contract.coverage_z,
        )
        records.append(record)
        if not run.is_analysable:
            labels.append("not-classified")
            continue
        report = classify_regime(
            run.series,
            run.dt_s,
            observable=name,
            amplitude_scale=scale,
            contract=regime_contract,
        )
        labels.append(report.regime.value)
    for index, run in enumerate(members):
        for other in members[index + 1 :]:
            differences.append(matched_time_difference(run, other))

    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    decided = {
        label: count
        for label, count in counts.items()
        if label not in {Regime.UNDECIDED.value, "not-classified"}
    }
    modal_label, modal_count = (
        max(decided.items(), key=lambda item: item[1]) if decided else (Regime.UNDECIDED.value, 0)
    )
    modal_fraction = modal_count / len(members)

    values = np.asarray([record.value for record in records], dtype=np.float64)
    finite = values[np.isfinite(values)]
    across_std = float(np.std(finite, ddof=1)) if finite.size >= 2 else float("nan")
    across_sem = across_std / math.sqrt(finite.size) if finite.size >= 2 else float("nan")
    mean_sem_ess = float(np.nanmean([record.sem_ess for record in records]))
    mean_sem_naive = float(np.nanmean([record.sem_naive for record in records]))
    calibration = across_std / mean_sem_ess if mean_sem_ess > 0.0 else float("nan")
    naive_calibration = across_std / mean_sem_naive if mean_sem_naive > 0.0 else float("nan")
    measured = {
        "regime_counts": tuple(sorted(counts.items())),
        "modal_regime": modal_label,
        "modal_fraction": modal_fraction,
        "values": tuple(float(value) for value in values),
        "value_mean": float(np.mean(finite)) if finite.size else float("nan"),
        "across_seed_std": across_std,
        "across_seed_sem": across_sem,
        "mean_sem_ess": mean_sem_ess,
        "mean_sem_naive": mean_sem_naive,
        "ess_calibration_ratio": calibration,
        "naive_calibration_ratio": naive_calibration,
        "mean_trajectory_difference_over_std": float(np.nanmean(differences)),
    }

    if any(label == "not-classified" for label in labels):
        return _report(
            InvarianceVerdict.FLAGGED_NUMERICAL_DIVERGENCE,
            "at least one seed produced a non-finite trajectory at a step every other seed survived",
            **measured,
        )
    if modal_fraction < contract.min_seed_regime_fraction:
        return _report(
            InvarianceVerdict.FLAGGED_REGIME_ARTEFACT,
            (
                f"the modal decided regime {modal_label!r} holds for {modal_fraction:.3g} of "
                f"{len(members)} seeds, below the declared "
                f"{contract.min_seed_regime_fraction:g}: the label depends on the stream, not on "
                f"the system. Across-seed statistic spread {across_std:.4g}"
            ),
            **measured,
        )
    return _report(
        InvarianceVerdict.INVARIANT,
        (
            f"{modal_count} of {len(members)} seeds classified {modal_label} "
            f"({modal_fraction:.3g} >= {contract.min_seed_regime_fraction:g}); the per-seed "
            f"{statistic.name} spread is {across_std:.4g} against a mean within-run ESS error of "
            f"{mean_sem_ess:.4g} (calibration ratio {calibration:.3g}) while the naive "
            f"sigma/sqrt(N) error would be {mean_sem_naive:.4g} (ratio {naive_calibration:.3g}). "
            f"The trajectories are NOT invariant and are not claimed to be: their mean pairwise "
            f"difference is {measured['mean_trajectory_difference_over_std']:.3g} standard "
            "deviations"
        ),
        **measured,
    )


# ------------------------------------------------------------------------------------------------
# Analytic systems.  These are the battery's controls: every one has a closed-form regime and, where
# it matters, a closed-form stability step, so what the battery reports can be checked against an
# oracle rather than against itself.
# ------------------------------------------------------------------------------------------------


def shared_brownian_increments(
    *, duration_s: float, base_dt_s: float, seed: int
) -> tuple[np.ndarray, int]:
    """Return the Wiener increments of ONE realisation at the finest step, and their count.

    A refinement ladder driven by independent noise streams compares two different realisations of the
    process as well as two discretisations, and the sampling difference is typically the larger of the
    two.  Summing consecutive increments of a single fine path to build the coarse ones keeps the
    Wiener path identical across every rung — the standard strong-convergence protocol — so what
    survives the comparison is the discretisation and nothing else.

    Args:
        duration_s: Physical duration [s].
        base_dt_s: The finest step [s]; every coarser rung must be an integer multiple of it.
        seed: RNG seed.

    Returns:
        The increments ``dW`` with variance ``base_dt_s``, and their count.
    """
    duration = _positive(duration_s, what="duration_s")
    base = _positive(base_dt_s, what="base_dt_s")
    steps = int(round(duration / base))
    if steps < 4:
        raise ValueError("need >= 4 base steps to build a refinement ladder")
    rng = np.random.default_rng(seed)
    return rng.normal(0.0, math.sqrt(base), steps), steps


def _coarsen(increments: np.ndarray, factor: int) -> np.ndarray:
    """Return the increments of the SAME Wiener path at ``factor`` times the base step."""
    usable = (increments.size // factor) * factor
    return increments[:usable].reshape(-1, factor).sum(axis=1)


def ou_series(
    *,
    theta_per_s: float,
    sigma: float,
    dt_s: float,
    increments: np.ndarray,
    scheme: str,
) -> np.ndarray:
    """Return an Ornstein-Uhlenbeck path driven by GIVEN Wiener increments.

    ``dx = -theta x dt + sigma dW``.  Two schemes, and the difference between them is the point:

    * ``"exact"`` — the exact discretisation ``x_{n+1} = a x_n + sigma sqrt((1-a^2)/(2 theta)) z`` with
      ``a = exp(-theta dt)``.  Its stationary variance is ``sigma^2 / (2 theta)`` at EVERY step, so any
      step dependence the battery reports on this scheme is a defect of the battery.  (It consumes the
      increments as scaled Gaussian draws ``z = dW / sqrt(dt)``, so the ladder stays coupled to one
      realisation.)
    * ``"euler-maruyama"`` — ``x_{n+1} = x_n - theta dt x_n + sigma dW``.  Stable for ``theta dt < 2``,
      but its stationary variance is ``sigma^2 / (theta (2 - theta dt))``, i.e. inflated by
      ``1 / (1 - theta dt / 2)``.  That is a genuine ``O(dt)`` discretisation bias in a physical
      statistic while the REGIME is unchanged, which is why this control is here: it separates
      "the regime is an artefact" from "this number is biased", and only the first is the falsifier.

    Args:
        theta_per_s: Relaxation rate [1/s].
        sigma: Noise amplitude [x units / sqrt(s)].
        dt_s: Step [s].
        increments: Wiener increments with variance ``dt_s``, from :func:`shared_brownian_increments`
            (coarsened as needed).
        scheme: ``"exact"`` or ``"euler-maruyama"``.

    Returns:
        The path, one sample per increment plus the initial condition dropped, so ``len(increments)``.

    Raises:
        ValueError: On an unknown scheme or a non-positive parameter.
    """
    theta = _positive(theta_per_s, what="theta_per_s")
    noise = _positive(sigma, what="sigma")
    dt = _positive(dt_s, what="dt_s")
    steps = np.asarray(increments, dtype=np.float64).ravel()
    if steps.size < 4:
        raise ValueError("need >= 4 increments")
    path = np.empty(steps.size, dtype=np.float64)
    if scheme == "exact":
        decay = math.exp(-theta * dt)
        amplitude = noise * math.sqrt((1.0 - decay * decay) / (2.0 * theta))
        draws = steps / math.sqrt(dt)
        state = 0.0
        for index in range(steps.size):
            state = decay * state + amplitude * draws[index]
            path[index] = state
    elif scheme == "euler-maruyama":
        state = 0.0
        for index in range(steps.size):
            state = state - theta * dt * state + noise * steps[index]
            path[index] = state
    else:
        raise ValueError(f"unknown scheme {scheme!r}; use 'exact' or 'euler-maruyama'")
    return path


def ou_refinement_ladder(
    *,
    theta_per_s: float,
    sigma: float,
    base_dt_s: float,
    factors: Sequence[int],
    duration_s: float,
    seed: int,
    scheme: str,
) -> tuple[LadderRun, ...]:
    """Return a ``dt`` ladder of OU paths driven by ONE shared Wiener realisation.

    Every rung covers the same physical duration by construction — which is what the battery's
    duration check is there to verify rather than to assume — and is driven by the same Brownian path,
    so the only difference between rungs is the discretisation.

    Args:
        theta_per_s: Relaxation rate [1/s].
        sigma: Noise amplitude.
        base_dt_s: The finest step [s].
        factors: Integer multiples of the base step, one per rung.
        duration_s: Physical duration [s], identical on every rung.
        seed: RNG seed for the shared path.
        scheme: ``"exact"`` or ``"euler-maruyama"``.

    Returns:
        One :class:`LadderRun` per factor, knob = the rung's ``dt`` [s].
    """
    increments, steps = shared_brownian_increments(
        duration_s=duration_s, base_dt_s=base_dt_s, seed=seed
    )
    runs: list[LadderRun] = []
    for factor in factors:
        multiple = _positive_int(factor, what="factor")
        if steps % multiple:
            raise ValueError(
                f"factor {multiple} does not divide the {steps} base steps; a rung that cannot be "
                "built from whole base increments would cover a different physical duration"
            )
        dt = base_dt_s * multiple
        runs.append(
            LadderRun(
                knob_value=dt,
                dt_s=dt,
                series=ou_series(
                    theta_per_s=theta_per_s,
                    sigma=sigma,
                    dt_s=dt,
                    increments=_coarsen(increments, multiple),
                    scheme=scheme,
                ),
                seed=seed,
            )
        )
    return tuple(runs)


def euler_stability_step_s(*, frequency_hz: float, damping_ratio: float) -> float:
    """Return the explicit-Euler stability step ``dt* = 2 zeta / omega`` [s] of a damped oscillator.

    The ORACLE this battery's positive control is checked against.  Forward Euler applied to
    ``x'' + 2 zeta omega x' + omega^2 x = 0`` has per-step amplification
    ``|A|^2 = 1 - 2 zeta omega dt + omega^2 dt^2``.  The scheme therefore injects anti-damping
    ``omega^2 dt / 2`` per unit time, and at ``dt = 2 zeta / omega`` that exactly cancels the physical
    damping: below it the oscillator decays, above it a sustained and growing oscillation appears that
    is nowhere in the differential equation.

    Args:
        frequency_hz: Undamped natural frequency [Hz].
        damping_ratio: ``zeta``, dimensionless.

    Returns:
        The stability step [s].
    """
    omega = 2.0 * math.pi * _positive(frequency_hz, what="frequency_hz")
    return 2.0 * _positive(damping_ratio, what="damping_ratio") / omega


def euler_growth_per_step(*, frequency_hz: float, damping_ratio: float, dt_s: float) -> float:
    """Return the explicit-Euler per-step amplitude factor ``sqrt(1 - 2 zeta omega dt + omega^2 dt^2)``.

    Above 1 the discrete solution grows without any physical mechanism for it.  Reported next to the
    battery's verdict so the artefact has a magnitude and not just a name.

    Args:
        frequency_hz: Undamped natural frequency [Hz].
        damping_ratio: ``zeta``, dimensionless.
        dt_s: Step [s].

    Returns:
        The amplitude factor per step, dimensionless.
    """
    omega = 2.0 * math.pi * _positive(frequency_hz, what="frequency_hz")
    zeta = _positive(damping_ratio, what="damping_ratio")
    dt = _positive(dt_s, what="dt_s")
    return math.sqrt(max(0.0, 1.0 - 2.0 * zeta * omega * dt + omega * omega * dt * dt))


def damped_oscillator_series(
    *,
    frequency_hz: float,
    damping_ratio: float,
    sigma: float,
    dt_s: float,
    increments: np.ndarray,
    scheme: str,
) -> np.ndarray:
    """Return the position of a noise-driven damped oscillator under an explicit or semi-implicit step.

    ``dx = v dt``, ``dv = (-omega^2 x - 2 zeta omega v) dt + sigma dW``.  The continuous system is a
    stable spiral for every ``zeta > 0``: with small noise it sits at rest, and no oscillation is
    sustained.  What the two schemes do with that is the whole positive control.

    * ``"explicit-euler"`` evaluates both derivatives at the old state.  Its amplification factor is
      :func:`euler_growth_per_step`, so above :func:`euler_stability_step_s` it grows a sustained
      oscillation out of the noise floor. **The regime changes, the system does not.**
    * ``"semi-implicit"`` (symplectic Euler) uses the updated velocity in the position update.  It has
      no such anti-damping and stays a decaying spiral for ``omega dt < 2``, at the SAME step where the
      explicit scheme has already invented a cycle — which is what makes the flagged verdict a
      statement about the method rather than about the system.

    Args:
        frequency_hz: Undamped natural frequency [Hz].
        damping_ratio: ``zeta``, dimensionless.  ``zeta < 1`` is underdamped.
        sigma: Noise amplitude on the velocity equation.
        dt_s: Step [s].
        increments: Wiener increments with variance ``dt_s``.
        scheme: ``"explicit-euler"`` or ``"semi-implicit"``.

    Returns:
        The position series, one sample per increment.  May contain non-finite values when the scheme
        diverges; that is a result, not an error, and :class:`LadderRun` accepts it.

    Raises:
        ValueError: On an unknown scheme or a non-positive parameter.
    """
    omega = 2.0 * math.pi * _positive(frequency_hz, what="frequency_hz")
    zeta = _positive(damping_ratio, what="damping_ratio")
    noise = _positive(sigma, what="sigma")
    dt = _positive(dt_s, what="dt_s")
    steps = np.asarray(increments, dtype=np.float64).ravel()
    if steps.size < 4:
        raise ValueError("need >= 4 increments")
    if scheme not in {"explicit-euler", "semi-implicit"}:
        raise ValueError(f"unknown scheme {scheme!r}; use 'explicit-euler' or 'semi-implicit'")
    path = np.empty(steps.size, dtype=np.float64)
    position, velocity = 0.0, 0.0
    with np.errstate(over="ignore", invalid="ignore"):
        for index in range(steps.size):
            acceleration = -omega * omega * position - 2.0 * zeta * omega * velocity
            new_velocity = velocity + dt * acceleration + noise * steps[index]
            if scheme == "explicit-euler":
                position = position + dt * velocity
            else:
                position = position + dt * new_velocity
            velocity = new_velocity
            path[index] = position
    return path


@dataclass(frozen=True, slots=True, eq=False)
class ImplicitSolveRun:
    """A backward-Euler run and the inner solve's own record of what it did.

    Carries both halves on purpose: ``state`` is the physics, ``final_residual`` and ``iterations`` are
    the solver's internal bookkeeping.  The adversarial case this battery exists to catch is a caller
    declaring the second half physical, so both must be available from one run for the mistake to be
    reproducible.
    """

    dt_s: float
    tolerance: float
    max_iterations: int
    state: np.ndarray
    final_residual: np.ndarray
    iterations: np.ndarray

    def ladder_run(self, knob_value: float, *, seed: int | None = None) -> LadderRun:
        """Return the rung: the trajectory as the observable, the solver's record as channels.

        Both on one rung so that a ladder can hold the physical anchor and the suspect quantity
        together; the residual-tracking route is only available when it can see that the trajectory
        did not move while the other quantity did.
        """
        return LadderRun(
            knob_value=knob_value,
            dt_s=self.dt_s,
            series=self.state,
            seed=seed,
            channels={"final_residual": self.final_residual, "iterations": self.iterations},
        )


def implicit_relaxation_run(
    *,
    theta_per_s: float,
    sigma: float,
    dt_s: float,
    increments: np.ndarray,
    tolerance: float,
    max_iterations: int,
) -> ImplicitSolveRun:
    """Integrate ``dx = -theta x dt + sigma dW`` by backward Euler with a fixed-point inner solve.

    The inner iteration is ``y_{k+1} = x_n - theta dt y_k + sigma dW``, which contracts at rate
    ``theta dt`` and therefore converges for ``theta dt < 1``.  It stops when the update falls below
    ``tolerance`` or the budget is spent, and BOTH stopping rules are exposed, because the two knobs
    fail differently: a loose tolerance leaves a residual proportional to itself, while a small budget
    leaves one proportional to ``(theta dt)^budget``.

    Handing the same ``increments`` to several tolerances or budgets couples them to one Wiener path,
    so the physics is identical by construction and any statistic that still moves is moving with the
    solver.

    Args:
        theta_per_s: Relaxation rate [1/s].
        sigma: Noise amplitude.
        dt_s: Step [s].
        increments: Wiener increments with variance ``dt_s``.
        tolerance: Inner-solve stopping tolerance in the state's units.  ``0.0`` means "never
            satisfied", which is how the iteration-budget axis is exercised in isolation.
        max_iterations: Inner-iteration budget per step, at least 1.

    Returns:
        The :class:`ImplicitSolveRun`.

    Raises:
        ValueError: If the fixed-point iteration is not contractive at this step, which would make the
            "converged" branch meaningless.
    """
    theta = _positive(theta_per_s, what="theta_per_s")
    noise = _positive(sigma, what="sigma")
    dt = _positive(dt_s, what="dt_s")
    tol = _nonnegative(tolerance, what="tolerance")
    budget = _positive_int(max_iterations, what="max_iterations")
    steps = np.asarray(increments, dtype=np.float64).ravel()
    if steps.size < 4:
        raise ValueError("need >= 4 increments")
    if theta * dt >= 1.0:
        raise ValueError(
            f"the fixed-point inner solve contracts at rate theta*dt = {theta * dt:g}; at or above 1 "
            "it does not converge at all, and a tolerance ladder on a non-convergent solve measures "
            "the divergence rate rather than the tolerance"
        )
    state = np.empty(steps.size, dtype=np.float64)
    residual = np.empty(steps.size, dtype=np.float64)
    counts = np.empty(steps.size, dtype=np.float64)
    current = 0.0
    for index in range(steps.size):
        drive = current + noise * steps[index]
        guess = current
        update = float("inf")
        used = 0
        for _ in range(budget):
            nxt = drive - theta * dt * guess
            update = abs(nxt - guess)
            guess = nxt
            used += 1
            if update < tol:
                break
        current = guess
        state[index] = current
        residual[index] = update
        counts[index] = float(used)
    return ImplicitSolveRun(
        dt_s=dt,
        tolerance=tol,
        max_iterations=budget,
        state=state,
        final_residual=residual,
        iterations=counts,
    )


def invariance_experiment_card(
    *, manifest_hashes: Sequence[str], budget_class: str
) -> SandboxExperimentCard:
    """Return the pre-registered card for the S0 numerical-invariance battery.

    Pre-registered means the controls and the falsifier are written down before the ladders run, so a
    verdict cannot be reinterpreted after the fact.  The positive control is the case that MUST be
    flagged; a battery with no such case has not been shown to be able to fail.

    Args:
        manifest_hashes: The cell-state manifests this battery is registered against.  Caller-supplied;
            the battery does not invent a provenance for itself.
        budget_class: The declared cost class of the experiment.

    Returns:
        The :class:`aleph.virtual_cell.contracts.SandboxExperimentCard`.
    """
    return SandboxExperimentCard(
        experiment_id="S0-numerical-invariance-battery",
        question=(
            "Is the classified regime a property of the system, or of the numerical method that "
            "produced the trajectory?"
        ),
        manifest_hashes=tuple(manifest_hashes),
        intervention=(
            "Refine the step, tighten the inner-solve tolerance, enlarge the iteration budget and "
            "change the seed, holding the physical duration and the Wiener realisation fixed"
        ),
        observables=(
            "regime label per rung",
            "time-averaged observable with its ESS interval",
            "time-averaged squared deviation with its ESS interval",
            "per-step final inner residual",
            "across-seed spread of the time average",
        ),
        positive_controls=(
            "explicit-Euler damped noisy oscillator stepped past dt* = 2 zeta / omega: the battery "
            "MUST flag a regime change, and the flagged bracket must contain the analytic dt*",
        ),
        negative_controls=(
            "exact-discretisation Ornstein-Uhlenbeck over the same dt ladder: regime invariant while "
            "the trajectories are not",
            "backward-Euler tolerance and budget ladders on one Wiener path: the trajectory statistics "
            "must not move",
        ),
        adversarial_controls=(
            "the per-step final inner residual DECLARED physical: it moves by orders of magnitude "
            "with the tolerance and the budget while the trajectory is unchanged, and must be routed "
            "as residual-tracking rather than as an invariance failure of the physics",
            "the same ladders judged with naive sigma/sqrt(N) errors instead of ESS errors, which "
            "makes the finest rung look most precise purely because its samples are most correlated",
        ),
        held_out_interventions=(
            "the same oscillator at the same step under the semi-implicit scheme, which must NOT flag",
            "a seed ensemble not used to set any threshold in the contract",
        ),
        falsifier=(
            "The classified regime, or a declared-physical statistic beyond its own ESS interval, "
            "changes with dt, inner-solve tolerance, iteration budget or seed at matched physical "
            "duration"
        ),
        failure_expansions=(
            ExpansionRoute.SOLVER_DESIGN.value,
            ExpansionRoute.ADAPTIVE_DESIGN.value,
            ExpansionRoute.TIME_MODALITY.value,
            ExpansionRoute.OBSERVATION_MODEL_EXTENSION.value,
        ),
        budget_class=budget_class,
        wet_lab=False,
    )
