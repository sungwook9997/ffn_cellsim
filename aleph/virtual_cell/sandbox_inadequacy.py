"""S12 model-inadequacy sandbox: when a law is removed on purpose, does the residual find it?

The master plan's S12 row asks exactly one question — *의도적으로 빠진 law를 residual이 찾는가?* — and
names its own falsifier: **a parameter shift hides the missing law.**  That falsifier is not
hypothetical here.  It is the failure mode every calibrated-but-wrong model in this project's history
used to survive: the fit was refreshed, the residual came back down, and the missing mechanism was
never named because nothing in the pipeline distinguished *"the model is right"* from *"the model is
wrong and its free parameter absorbed the difference"*.

This module builds analytic systems where the truth is constructed WITH a law and the model is fitted
WITHOUT it, and measures whether four independent detector channels find the omission.  Everything is
``oracle-only``: no cell observable is predicted, no engine number is read, no gate is closed.  What is
validated is the DETECTOR, on systems whose answer is known by construction.

The central design decision, and the reason the result is a contrast rather than a score
-----------------------------------------------------------------------------------------
Detectability is **not a property of the omitted law alone.**  It is a property of the omitted law
*relative to the intervention that was varied*.  Two constructed omissions make this exact:

``CO_SCALING``
    The omitted process runs at a rate proportional to the intervention, exactly like the process the
    model does represent.  Its whole first-order signature is therefore degenerate with the model's own
    fitted rate.  A refit absorbs it — and, crucially, keeps absorbing it at every amplitude and under
    every additional condition, because the absorption is an algebraic identity, not a coincidence.

``ABSOLUTE_TIME``
    The identical omitted process, at the identical amplitude and shape, except that its rate does NOT
    respond to the intervention — it runs on absolute time.  Now the omission's signature varies from
    condition to condition in a way one shared rate cannot reproduce.

Same law, same size, same shape.  Only the coupling to the knob differs, and that alone decides whether
the omission is findable.  A programme that varies one knob and refits cannot, even in principle, see a
law that co-varies with that knob.

The four channels, and what each one is for
-------------------------------------------
``D1`` **curvature** — a nested F-test adding one quadratic basis function per condition, reported
    alongside ``ΔBIC``.  Detects leftover pointwise shape after per-condition refitting.
``D2`` **serial structure** — the residual's integrated autocorrelation time against a *permutation*
    null.  The permutation preserves the residual's magnitude exactly and destroys only its ordering,
    so this channel cannot fire because a residual is large.  That is the whole point: the plan's
    routing says a *persistent* residual becomes a missing-law candidate, and persistence means
    structured, not big.
``D3`` **shared-parameter inconsistency** — fit all conditions with the parameters they are supposed to
    share, and F-test that against per-condition refitting.  This is the strong test.  Per-condition
    refitting is exactly the operation that hides a missing law; forbidding it and asking whether the
    misfit can still be absorbed *simultaneously* is what refitting cannot survive.
``D4`` **held-out intervention** — fit the shared law on the training interventions, then predict at an
    intervention never fitted, with only the per-experiment amplitude nuisance free.  A model that
    interpolates inside its fitted conditions is not validated.

The family-wise level is split evenly across the channels that actually ran (Bonferroni).  Every
threshold is a field of :class:`InadequacyContract`, which has **no defaults on any field**, matching
``RegimeEvidenceContract`` for the same reason: a number invented inside the detector is a number
chosen after seeing the data.

What this module deliberately does NOT do
-----------------------------------------
It never asserts a missing law.  :class:`MissingLawProposal` is a proposal carrying evidence; its
``evidence_source`` is a property fixed to ``ANALYTIC_ORACLE`` with no field to overwrite, it reports
``decided is False`` unconditionally, and converting it into an archetype-level
:class:`~aleph.virtual_cell.archetypes.MissingLawVerdict` requires an external adjudication reference
and raises :class:`SelfPromotionError` when the detector's own evidence is offered as that reference.
An inadequacy detector proposes; the PI decides.

It also does not import :func:`~aleph.virtual_cell.sweep.finite_difference_sensitivity` or
:func:`~aleph.virtual_cell.regime.classify_regime`.  Both would require this module to instantiate a
``RegimeEvidenceContract``, whose thirty fields have no defaults by design; inventing thirty threshold
values here to reach a descriptive statistic would be exactly the practice the contract exists to
prevent.  The regime statistics that ARE used —
:func:`~aleph.virtual_cell.regime.integrated_autocorrelation_time` and
:func:`~aleph.virtual_cell.sweep.ess_interval` — need only the caller's Sokal constant and coverage
factor, both of which are contract fields.

Sanity Gate (recorded before first execution):

* **Dimensional** — the working observable is ``u = log y``, a pure log-amplitude, so every residual,
  RMS and tolerance in this module is dimensionless.  Time enters only as the dimensionless
  ``x = r_ref * c * t``: each condition is measured over the same number of e-folds of the process the
  model *does* contain.  That protocol is what makes "the shared model forces one slope in ``x``" an
  exact algebraic statement rather than an approximation, and it is declared before any data.
* **Boundary** — a design matrix with fewer rows than columns, a condition list shorter than
  ``min_conditions``, a series shorter than ``min_samples_per_condition``, a zero or non-finite
  intervention value, and a duplicated intervention value all raise rather than returning a fit.  An
  ``F`` statistic at exactly zero returns ``p = 1`` and never a negative variance ratio; the
  permutation p-value is ``(hits + 1) / (trials + 1)``, so it is bounded strictly inside ``(0, 1]`` and
  can never be reported as zero evidence of nothing.
* **Conservation / invariant** — every channel is *scale-invariant by construction*: D1, D3 and D4 are
  ratios of sums of squares and D2 is a permutation test over a fixed multiset of residuals.
  Multiplying every observation's noise by any positive constant therefore leaves all four p-values
  bit-identical.  A magnitude-based inadequacy detector cannot make that claim, and this one is tested
  for it rather than argued for it.  Nested model comparisons satisfy ``SSR_restricted >= SSR_free`` to
  within least-squares round-off, and a violation raises.
* **Measurement protocol** — the false-positive rate is measured over an explicit trial count with the
  seed recorded in the report, against the declared family-wise level and its binomial three-sigma
  band; the band is stated in the contract before the run and never re-thresholded afterwards.  The
  contract also declares the noise model the channels' null distributions assume, and
  :func:`measure_null_false_positive_rate` accepts a *violating* noise model so the cost of that
  assumption is measured rather than asserted.  The binomial band is a statement about INDEPENDENT
  replicates, so the streams are not left to integer arithmetic: :func:`derive_replicate_streams`
  spawns one ``numpy.random.SeedSequence`` per (replicate, role) pair from the recorded root seed,
  and the report carries the derivation and the stream count so the premise is auditable from the
  artifact.  The rule this replaced offset three roles by different integers and, at the committed
  contract where the data seed equalled the permutation seed, handed two roles the same stream —
  the premise failed exactly where the report asserted it held.
* **Numerical** — all fits are ordinary least squares on explicitly built design matrices solved by SVD
  (``rcond=None``); a rank-deficient design raises instead of returning a minimum-norm fit that would
  understate the residual.  The F survival function is the regularized incomplete beta evaluated by the
  Lentz continued fraction with the standard symmetry swap at ``x > (a+1)/(a+b+2)``, which keeps the
  relative error near machine epsilon in both tails; the tests check it against ``scipy.stats.f``,
  which is the plan's "solve the same question in two different mathematical expressions" applied to
  the one closed form this module depends on.  An upper tail that underflows double precision is
  returned as the smallest positive normal double rather than as zero, because zero would assert an
  exact impossibility that the arithmetic cannot support.  Non-finite inputs raise at construction.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

import numpy as np

from aleph.virtual_cell.archetypes import ArchetypeRing, MissingLaw, MissingLawVerdict
from aleph.virtual_cell.contracts import (
    CellStateManifest,
    EvidenceSource,
    SandboxExperimentCard,
)
from aleph.virtual_cell.regime import integrated_autocorrelation_time
from aleph.virtual_cell.source_registry import (
    EvidenceGraph,
    ProvenanceKind,
    SourceOrigin,
    SourceRecord,
    SourceRegistry,
    SourceStatus,
)
from aleph.virtual_cell.sweep import ClosureGapReport, ESSInterval, ExpansionRoute
from aleph.virtual_cell.sweep import ess_interval as _ess_interval
from aleph.virtual_cell.sweep import moment_closure_gap as _moment_closure_gap

__all__ = [
    "CANDIDATE_LAWS",
    "EVIDENCE_AUTHORITY",
    "OMITTED_MODE_WEIGHT",
    "REFERENCE_RATE_HZ",
    "REPLICATE_STREAM_ROLES",
    "S12_MANIFEST",
    "SANDBOX_S12_CARD",
    "STREAM_DERIVATION",
    "WINDOW_EFOLDS",
    "ChannelId",
    "ChannelResult",
    "ConditionSeries",
    "InadequacyContract",
    "MissingLawProposal",
    "NoCandidateReport",
    "NoiseModel",
    "NullControlReport",
    "OmissionKind",
    "ProposalStatus",
    "ReplicateStreams",
    "ResidualStructureReport",
    "SelfPromotionError",
    "analyse_conditions",
    "build_conditions",
    "derive_replicate_streams",
    "detect_missing_law",
    "experiment_fingerprint",
    "f_survival",
    "measure_null_false_positive_rate",
    "parameter_masking_fraction",
    "rate_closure_gap",
    "s12_evidence_graph",
]

# ------------------------------------------------------------------------------------------------
# Protocol constants of the constructed oracle systems.
#
# These define the TRUTH being constructed; they are not thresholds and nothing is gated on them.
# Every threshold lives in InadequacyContract, which has no defaults.
# ------------------------------------------------------------------------------------------------

#: Reference relaxation rate of the process the reduced model DOES contain [1/s].
REFERENCE_RATE_HZ: float = 1.0

#: Weight of the omitted second relaxation mode.  0.3 is a genuine minority mode: large enough that
#: the omission is a real law and not a rounding error, small enough that the observable still looks
#: like a single relaxation to the eye — which is the situation the detector exists for.
OMITTED_MODE_WEIGHT: float = 0.3

#: Observation window in e-folds of the represented process.  Each condition is measured over the same
#: dimensionless span; see the module Sanity Gate for why that is the protocol and not a convenience.
WINDOW_EFOLDS: float = 4.0

#: The evidence authority this module can ever carry.  A constant, not a parameter.
EVIDENCE_AUTHORITY = (
    "oracle-only: the truth is constructed in closed form. Nothing here is native evidence, "
    "and no verdict from this module may close a physics gate."
)


class OmissionKind(StrEnum):
    """Which law was deliberately removed from the model, and how it couples to the intervention."""

    NONE = "none"
    CO_SCALING = "co-scaling"
    ABSOLUTE_TIME = "absolute-time"


class NoiseModel(StrEnum):
    """The observation-noise process.  Only the first is the one the channels' nulls assume.

    The others exist so the cost of that assumption can be MEASURED.  A detector whose stated level
    holds only under its own noise model, and which never checks, is a detector whose level is a claim
    rather than a result.
    """

    HOMOSCEDASTIC_GAUSSIAN = "homoscedastic-gaussian"
    HETEROSCEDASTIC = "heteroscedastic"
    HEAVY_TAILED = "heavy-tailed"
    OUTLIER_CONTAMINATED = "outlier-contaminated"


class ChannelId(StrEnum):
    """The four detector channels."""

    CURVATURE = "D1-curvature"
    SERIAL_STRUCTURE = "D2-serial-structure"
    SHARED_PARAMETER = "D3-shared-parameter"
    HELD_OUT_INTERVENTION = "D4-held-out-intervention"


class ProposalStatus(StrEnum):
    """Whether the detector proposed a candidate.  Neither value is a decision."""

    CANDIDATE_PROPOSED = "candidate-proposed"
    NO_CANDIDATE = "no-candidate"


class SelfPromotionError(RuntimeError):
    """Raised when a proposal is asked to become a verdict on its own evidence."""


# ------------------------------------------------------------------------------------------------
# validation helpers
# ------------------------------------------------------------------------------------------------


def _positive(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{what} must be positive and finite; got {value!r}")
    return number


def _nonneg(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{what} must be nonnegative and finite; got {value!r}")
    return number


def _positive_int(value: object, *, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{what} must be an int")
    if value <= 0:
        raise ValueError(f"{what} must be positive; got {value!r}")
    return value


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _as_seed_sequence(value: object, *, what: str) -> np.random.SeedSequence:
    """Accept a nonnegative int or a ``SeedSequence`` and return a ``SeedSequence``.

    ``np.random.default_rng(n)`` builds ``SeedSequence(n)`` internally, so promoting an integer here
    changes no existing stream; it only makes the stream an addressable object that can be spawned.
    """
    if isinstance(value, np.random.SeedSequence):
        return value
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{what} must be a nonnegative int or a numpy.random.SeedSequence")
    return np.random.SeedSequence(value)


def _stream_fingerprint(stream: np.random.SeedSequence) -> tuple[int, ...]:
    """128 bits of the state a ``SeedSequence`` would hand a generator.

    ``generate_state`` is a pure function of ``(entropy, spawn_key, pool_size)``, so calling it here
    consumes nothing and does not perturb any generator later built from the same object.
    """
    return tuple(int(word) for word in stream.generate_state(4, dtype=np.uint32))


def _array(values: Any, *, what: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64).ravel()
    if array.size == 0:
        raise ValueError(f"{what} must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{what} must contain only finite values")
    return array


# ------------------------------------------------------------------------------------------------
# F distribution, in closed form, so the module keeps numpy as its only dependency
# ------------------------------------------------------------------------------------------------


def _beta_continued_fraction(a: float, b: float, x: float) -> float:
    """Lentz evaluation of the continued fraction for the incomplete beta function."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        num = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + num * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + num / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        num = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + num * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + num / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-16:
            return h
    raise ArithmeticError("incomplete beta continued fraction did not converge")


def _regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    if x < (a + 1.0) / (a + b + 2.0):
        front = math.exp(a * math.log(x) + b * math.log1p(-x) - log_beta)
        return front * _beta_continued_fraction(a, b, x) / a
    front = math.exp(b * math.log1p(-x) + a * math.log(x) - log_beta)
    return 1.0 - front * _beta_continued_fraction(b, a, 1.0 - x) / b


def f_survival(statistic: float, df_numerator: int, df_denominator: int) -> float:
    """Return ``P(F > statistic)`` for the F distribution.

    Written out rather than taken from SciPy so the tests can check the one closed form this module
    depends on against ``scipy.stats.f.sf`` — the same question answered by two different
    implementations, which is the only kind of self-check worth having on a special function.

    The result is floored at :data:`sys.float_info.min`.  A p-value returned as exactly zero would
    assert that the observation is impossible under the null; double precision cannot support that
    claim, and a channel that reported it would be making a stronger statement than its arithmetic.

    Args:
        statistic: The observed variance ratio.  Zero or negative returns 1.0, which is the honest
            reading of "the restricted model fitted at least as well".
        df_numerator: Numerator degrees of freedom.
        df_denominator: Denominator degrees of freedom.

    Returns:
        The upper-tail probability, in ``(0, 1]``.

    Raises:
        ValueError: If either degree-of-freedom count is not positive, or the statistic is not finite.
    """
    d1 = _positive_int(df_numerator, what="df_numerator")
    d2 = _positive_int(df_denominator, what="df_denominator")
    if isinstance(statistic, bool) or not isinstance(statistic, (int, float)):
        raise TypeError("statistic must be a real number")
    value = float(statistic)
    if not math.isfinite(value):
        raise ValueError("statistic must be finite")
    if value <= 0.0:
        return 1.0
    tail = _regularized_incomplete_beta(d2 / 2.0, d1 / 2.0, d2 / (d2 + d1 * value))
    return max(tail, sys.float_info.min)


# ------------------------------------------------------------------------------------------------
# contract
# ------------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class InadequacyContract:
    """Pre-declared detector thresholds.  **No field has a default, by design.**

    Everything the detector compares a statistic against lives here, so the verdict can be re-judged
    from the report alone and so no number can be chosen after seeing a residual.

    Args:
        family_alpha: Family-wise false-positive level across the channels that run.  The per-channel
            level is ``family_alpha / (number of channels that ran)`` — a Bonferroni split, derived
            rather than declared, because the channel count is a property of the run.
        permutation_trials: Permutations in the D2 null.  Must be large enough that the smallest
            attainable p-value ``1 / (trials + 1)`` is below the per-channel level, otherwise D2 can
            never fire and silently contributes nothing; :meth:`per_channel_alpha` checks this.
        permutation_seed: Seed for the D2 permutation ensemble, recorded so a verdict is reproducible.
        sokal_window_c: Sokal automatic-window constant handed to
            :func:`~aleph.virtual_cell.regime.integrated_autocorrelation_time`.
        coverage_z: Coverage factor for the reported ESS interval on the mean residual.
        min_samples_per_condition: Shortest per-condition series that may be analysed at all.
        min_conditions: Fewest training conditions.  D3 needs at least three: with two conditions the
            per-condition and shared fits differ by one degree of freedom and the shared restriction is
            nearly vacuous, so a cross-condition inconsistency test on two conditions is not the test
            it claims to be.
        false_positive_sigma: Multiplier on the binomial standard error used to judge a measured
            false-positive rate against ``family_alpha``.  Declared before any trial is run.
        assumed_noise_model: The observation-noise process under which the channels' null
            distributions are exact.  Recorded in every report so a reader can see what the stated
            level is conditional on.
    """

    family_alpha: float
    permutation_trials: int
    permutation_seed: int
    sokal_window_c: float
    coverage_z: float
    min_samples_per_condition: int
    min_conditions: int
    false_positive_sigma: float
    assumed_noise_model: NoiseModel

    def __post_init__(self) -> None:
        alpha = _positive(self.family_alpha, what="family_alpha")
        if alpha >= 1.0:
            raise ValueError("family_alpha must be below 1")
        object.__setattr__(self, "family_alpha", alpha)
        object.__setattr__(
            self,
            "permutation_trials",
            _positive_int(self.permutation_trials, what="permutation_trials"),
        )
        if isinstance(self.permutation_seed, bool) or not isinstance(self.permutation_seed, int):
            raise TypeError("permutation_seed must be an int")
        if self.permutation_seed < 0:
            raise ValueError("permutation_seed must be nonnegative")
        object.__setattr__(
            self, "sokal_window_c", _positive(self.sokal_window_c, what="sokal_window_c")
        )
        object.__setattr__(self, "coverage_z", _positive(self.coverage_z, what="coverage_z"))
        object.__setattr__(
            self,
            "min_samples_per_condition",
            _positive_int(self.min_samples_per_condition, what="min_samples_per_condition"),
        )
        if self.min_samples_per_condition < 8:
            raise ValueError(
                "min_samples_per_condition below 8 leaves fewer residual degrees of freedom than the "
                "extended model has parameters; the curvature channel would not be a test"
            )
        object.__setattr__(
            self, "min_conditions", _positive_int(self.min_conditions, what="min_conditions")
        )
        if self.min_conditions < 3:
            raise ValueError(
                "min_conditions below 3 makes the shared-parameter channel vacuous; see the field "
                "documentation for why two conditions is not a cross-condition test"
            )
        object.__setattr__(
            self,
            "false_positive_sigma",
            _positive(self.false_positive_sigma, what="false_positive_sigma"),
        )
        if not isinstance(self.assumed_noise_model, NoiseModel):
            raise TypeError("assumed_noise_model must be a NoiseModel")

    def per_channel_alpha(self, n_channels: int) -> float:
        """Return the Bonferroni per-channel level for a run with ``n_channels`` channels.

        Args:
            n_channels: How many channels actually produced a p-value.

        Returns:
            ``family_alpha / n_channels``.

        Raises:
            ValueError: If ``permutation_trials`` is too small for the permutation channel to reach
                that level.  A channel that cannot fire is worse than an absent one, because it is
                counted in the Bonferroni split while contributing no power.
        """
        count = _positive_int(n_channels, what="n_channels")
        alpha = self.family_alpha / count
        smallest = 1.0 / (self.permutation_trials + 1.0)
        if smallest > alpha:
            raise ValueError(
                f"permutation_trials={self.permutation_trials} bounds the D2 p-value below at "
                f"{smallest:.6g}, above the per-channel level {alpha:.6g}: D2 could never fire while "
                "still consuming a share of the family-wise budget"
            )
        return alpha

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able contract, copied verbatim into every report."""
        return {
            "family_alpha": self.family_alpha,
            "permutation_trials": self.permutation_trials,
            "permutation_seed": self.permutation_seed,
            "sokal_window_c": self.sokal_window_c,
            "coverage_z": self.coverage_z,
            "min_samples_per_condition": self.min_samples_per_condition,
            "min_conditions": self.min_conditions,
            "false_positive_sigma": self.false_positive_sigma,
            "assumed_noise_model": self.assumed_noise_model.value,
        }


# ------------------------------------------------------------------------------------------------
# the constructed systems
# ------------------------------------------------------------------------------------------------


def _shape(x: np.ndarray, separation: float) -> np.ndarray:
    """Log-deviation of a two-mode relaxation from the single mode the reduced model contains.

    ``y(x) = (1 - w) exp(-x) + w exp(-separation * x)``, and the returned quantity is
    ``log y(x) + x`` — what is left in the log-signal once the represented mode is divided out.  It is
    zero at ``x = 0``, so the omitted law adds a process and not an offset.
    """
    weight = OMITTED_MODE_WEIGHT
    return np.log((1.0 - weight) * np.exp(-x) + weight * np.exp(-separation * x)) + x


def _deviation(
    x: np.ndarray, *, omission: OmissionKind, separation: float, intervention: float
) -> np.ndarray:
    """Return the unit-RMS-normalised log-deviation the reduced model is missing, at one condition.

    Normalisation is by the REFERENCE condition's RMS — one constant for the whole experiment — so the
    deviation amplitude means the same thing at every condition and the two omission kinds are
    genuinely amplitude-matched.  Normalising per condition would hide the very condition-dependence
    that separates them.
    """
    if omission is OmissionKind.NONE:
        return np.zeros_like(x)
    reference = _shape(x, separation)
    scale = float(np.sqrt(np.mean((reference - reference.mean()) ** 2)))
    if scale <= 0.0:  # pragma: no cover - separation == 1 is refused upstream
        raise ValueError("degenerate shape: the omitted mode is indistinguishable from the model's")
    if omission is OmissionKind.CO_SCALING:
        raw = _shape(x, separation)
    else:
        # The omitted process runs on absolute time: in the condition's own dimensionless clock its
        # rate is divided by the intervention, so its shape differs from condition to condition.
        raw = _shape(x, separation / intervention)
    return (raw - raw.mean()) / scale


@dataclass(frozen=True, slots=True)
class ConditionSeries:
    """One intervention condition: its dimensionless clock and its observed log-signal.

    ``dimensionless_time`` is ``x = r_ref * c * t``.  Sampling every condition on the same ``x`` grid
    is the measurement protocol declared in the module Sanity Gate; it is what makes "the shared model
    forces one slope in ``x``" exact.
    """

    condition_id: str
    intervention: float
    dimensionless_time: np.ndarray
    log_signal: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "condition_id", _nonempty(self.condition_id, what="condition_id"))
        object.__setattr__(self, "intervention", _positive(self.intervention, what="intervention"))
        x = _array(self.dimensionless_time, what="dimensionless_time")
        u = _array(self.log_signal, what="log_signal")
        if x.size != u.size:
            raise ValueError("dimensionless_time and log_signal must have the same length")
        if x.size < 4:
            raise ValueError("a condition needs at least 4 samples")
        spacing = np.diff(x)
        if spacing.size and not np.allclose(spacing, spacing[0], rtol=1e-9, atol=0.0):
            raise ValueError("dimensionless_time must be evenly spaced")
        if spacing.size and spacing[0] <= 0.0:
            raise ValueError("dimensionless_time must be strictly increasing")
        object.__setattr__(self, "dimensionless_time", x)
        object.__setattr__(self, "log_signal", u)

    @property
    def dx(self) -> float:
        """Sample spacing of the dimensionless clock."""
        return float(self.dimensionless_time[1] - self.dimensionless_time[0])


def _noise(
    rng: np.random.Generator, *, model: NoiseModel, sigma: float, x: np.ndarray
) -> np.ndarray:
    n = x.size
    if model is NoiseModel.HOMOSCEDASTIC_GAUSSIAN:
        return rng.normal(0.0, sigma, n)
    if model is NoiseModel.HETEROSCEDASTIC:
        ramp = 0.3 + 1.7 * (x - x[0]) / (x[-1] - x[0])
        return rng.normal(0.0, 1.0, n) * sigma * ramp
    if model is NoiseModel.HEAVY_TAILED:
        # Student t with 3 dof, rescaled to unit variance so the comparison is at matched sigma.
        return sigma * rng.standard_t(3.0, n) / math.sqrt(3.0)
    contaminated = rng.normal(0.0, sigma, n)
    count = max(1, n // 50)
    where = rng.choice(n, count, replace=False)
    contaminated[where] += 8.0 * sigma
    return contaminated


def build_conditions(
    *,
    interventions: Sequence[float],
    omission: OmissionKind,
    separation: float,
    deviation_sigma_multiple: float,
    noise_sigma: float,
    n_samples: int,
    seed: int | np.random.SeedSequence,
    noise_model: NoiseModel = NoiseModel.HOMOSCEDASTIC_GAUSSIAN,
) -> tuple[ConditionSeries, ...]:
    """Construct the truth WITH a law and return what an experiment would observe.

    The reduced model that will be fitted to these series contains only the single relaxation the
    intervention drives.  Whatever ``omission`` adds is, by construction, the law that was removed.

    Args:
        interventions: The intervention levels ``c`` (dimensionless multipliers on the represented
            rate).  Must be distinct and positive.
        omission: Which law is missing from the model, and how it couples to the intervention.
        separation: Rate ratio of the omitted mode to the represented one.  Must exceed 1; a ratio of
            exactly 1 is not a second process.
        deviation_sigma_multiple: Amplitude of the omitted law's log-deviation, in units of
            ``noise_sigma``, measured at the reference condition.  This is what makes two omissions
            amplitude-matched: the same value means the same pre-refit deviation RMS.
        noise_sigma: Observation-noise scale on the log-signal.
        n_samples: Samples per condition.
        seed: RNG seed, or a ``numpy.random.SeedSequence`` naming the stream directly; every series is
            drawn from that one stream so the experiment is reproducible.  Passing a ``SeedSequence``
            is what lets a caller give two roles PROVABLY distinct streams (distinct spawn keys)
            rather than hoping two integers do not collide — see :func:`derive_replicate_streams`.
            ``default_rng(SeedSequence(n))`` and ``default_rng(n)`` are the same generator, so an
            integer caller sees no change.
        noise_model: The noise process.  Anything but ``HOMOSCEDASTIC_GAUSSIAN`` violates what the
            channels assume, which is the point of offering it.

    Returns:
        One :class:`ConditionSeries` per intervention, in the given order.

    Raises:
        ValueError: If the interventions are not distinct positives, ``separation <= 1`` for a real
            omission, or any scale is not positive-finite.
    """
    levels = tuple(_positive(value, what="intervention") for value in interventions)
    if not levels:
        raise ValueError("at least one intervention level is required")
    if len(set(levels)) != len(levels):
        raise ValueError("intervention levels must be distinct")
    if not isinstance(omission, OmissionKind):
        raise TypeError("omission must be an OmissionKind")
    if not isinstance(noise_model, NoiseModel):
        raise TypeError("noise_model must be a NoiseModel")
    ratio = _positive(separation, what="separation")
    if omission is not OmissionKind.NONE and ratio <= 1.0:
        raise ValueError("separation must exceed 1: a ratio of 1 is not a second process")
    amplitude = _nonneg(deviation_sigma_multiple, what="deviation_sigma_multiple")
    sigma = _positive(noise_sigma, what="noise_sigma")
    count = _positive_int(n_samples, what="n_samples")
    if count < 4:
        raise ValueError("n_samples must be at least 4")
    stream = _as_seed_sequence(seed, what="seed")

    x = np.linspace(0.0, WINDOW_EFOLDS, count)
    rng = np.random.default_rng(stream)
    built: list[ConditionSeries] = []
    for index, level in enumerate(levels):
        deviation = _deviation(x, omission=omission, separation=ratio, intervention=level)
        signal = -x + amplitude * sigma * deviation
        observed = signal + _noise(rng, model=noise_model, sigma=sigma, x=x)
        built.append(
            ConditionSeries(
                condition_id=f"c{index}@{level:g}",
                intervention=level,
                dimensionless_time=x,
                log_signal=observed,
            )
        )
    return tuple(built)


def parameter_masking_fraction(
    *, omission: OmissionKind, separation: float, n_samples: int, intervention: float = 1.0
) -> float:
    """Return the fraction of the omitted law's norm that a parameter shift can absorb.

    This is the a-priori answer to the falsifier, computable in closed form before any data exist:
    project the omitted law's log-deviation onto the reduced model's own column space ``span{1, x}``
    and take the norm ratio.  A value near 1 means the omission IS a parameter shift to within its
    remainder, and no amount of data collected at one condition will separate them.

    Args:
        omission: Which law is missing.
        separation: Rate ratio of the omitted mode.
        n_samples: Grid size, since the projection is taken on the sampled grid.
        intervention: Condition at which to evaluate.

    Returns:
        ``||P deviation|| / ||deviation||`` in ``[0, 1]``; ``0.0`` when nothing is omitted.

    Raises:
        ValueError: If the arguments are outside the ranges :func:`build_conditions` accepts.
    """
    if not isinstance(omission, OmissionKind):
        raise TypeError("omission must be an OmissionKind")
    if omission is OmissionKind.NONE:
        return 0.0
    ratio = _positive(separation, what="separation")
    if ratio <= 1.0:
        raise ValueError("separation must exceed 1")
    count = _positive_int(n_samples, what="n_samples")
    level = _positive(intervention, what="intervention")
    x = np.linspace(0.0, WINDOW_EFOLDS, count)
    deviation = _deviation(x, omission=omission, separation=ratio, intervention=level)
    centred = deviation - deviation.mean()
    total = float(centred @ centred)
    if total <= 0.0:
        return 0.0
    design = np.column_stack([np.ones(count), x])
    _, residual_ss, _ = _ols(design, centred)
    return math.sqrt(max(0.0, 1.0 - residual_ss / total))


def rate_closure_gap(
    *, separation: float, tolerance: float, observation_x: float
) -> ClosureGapReport:
    """Measure ``E[exp(-r x)] - exp(-E[r] x)`` for the omitted two-rate mixture.

    The reduced model does not merely mis-fit the two-mode truth; it *closes* the distribution over
    relaxation rates at its mean.  That closure cost is exactly what
    :func:`~aleph.virtual_cell.sweep.moment_closure_gap` measures, and it predicts the size of the
    pre-refit deviation before any series is generated.  Reusing that function rather than restating
    the arithmetic also keeps the interpretation attached: a nonzero gap is the price of the
    representation, not evidence that the reference dynamics are wrong.

    Args:
        separation: Rate ratio of the omitted mode to the represented one.
        tolerance: Declared limit on the relative gap.  Caller-supplied, with no default.
        observation_x: Dimensionless time at which the closure is evaluated.

    Returns:
        The :class:`~aleph.virtual_cell.sweep.ClosureGapReport`.

    Raises:
        ValueError: If ``separation <= 1`` or ``observation_x`` is not positive-finite.
    """
    ratio = _positive(separation, what="separation")
    if ratio <= 1.0:
        raise ValueError("separation must exceed 1")
    time_point = _positive(observation_x, what="observation_x")
    rates = np.array([1.0, ratio], dtype=np.float64)
    weights = np.array([1.0 - OMITTED_MODE_WEIGHT, OMITTED_MODE_WEIGHT], dtype=np.float64)
    return _moment_closure_gap(
        rates,
        lambda r: np.exp(-r * time_point),
        observable=f"exp(-r*x) at x={time_point:g}",
        tolerance=tolerance,
        weights=weights,
    )


# ------------------------------------------------------------------------------------------------
# least squares
# ------------------------------------------------------------------------------------------------


def _ols(design: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    """Solve one least-squares problem, refusing a rank-deficient design."""
    if design.ndim != 2:
        raise ValueError("design must be two-dimensional")
    if design.shape[0] != target.size:
        raise ValueError("design rows must match the target length")
    if design.shape[0] <= design.shape[1]:
        raise ValueError(
            f"design has {design.shape[0]} rows and {design.shape[1]} columns: with no residual "
            "degrees of freedom the fit reports a zero residual that means nothing"
        )
    coefficients, _, rank, _ = np.linalg.lstsq(design, target, rcond=None)
    if rank < design.shape[1]:
        raise ValueError(
            "rank-deficient design: a minimum-norm solution would understate the residual and make "
            "an inadequate model look adequate"
        )
    residual = target - design @ coefficients
    return coefficients, float(residual @ residual), residual


def _per_condition_fit(
    conditions: Sequence[ConditionSeries], *, degree: int
) -> tuple[float, int, dict[str, np.ndarray]]:
    """Fit ``{1, x, ..., x**degree}`` separately in every condition."""
    total = 0.0
    residuals: dict[str, np.ndarray] = {}
    for condition in conditions:
        x = condition.dimensionless_time
        design = np.column_stack([x**power for power in range(degree + 1)])
        _, ssr, residual = _ols(design, condition.log_signal)
        total += ssr
        residuals[condition.condition_id] = residual
    return total, (degree + 1) * len(conditions), residuals


def _shared_fit(
    conditions: Sequence[ConditionSeries],
) -> tuple[float, int, float, dict[str, np.ndarray]]:
    """Fit one shared rate across conditions, with a free amplitude nuisance in each.

    In the dimensionless clock the represented process has slope exactly ``-1`` at every condition, so
    "the conditions share the rate law" is the single-column restriction below.  The per-condition
    amplitude stays free because a normalisation is a property of the measurement, not of the law.
    """
    count = len(conditions)
    blocks: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    for index, condition in enumerate(conditions):
        x = condition.dimensionless_time
        block = np.zeros((x.size, count + 1), dtype=np.float64)
        block[:, index] = 1.0
        block[:, -1] = -x
        blocks.append(block)
        targets.append(condition.log_signal)
    design = np.vstack(blocks)
    target = np.concatenate(targets)
    coefficients, ssr, residual = _ols(design, target)
    residuals: dict[str, np.ndarray] = {}
    offset = 0
    for condition in conditions:
        size = condition.dimensionless_time.size
        residuals[condition.condition_id] = residual[offset : offset + size]
        offset += size
    return ssr, count + 1, float(coefficients[-1]), residuals


def _nested_f(
    *, ssr_restricted: float, ssr_free: float, params_restricted: int, params_free: int, n_rows: int
) -> tuple[float, int, int, float]:
    """Return ``(F, df1, df2, p)`` for a nested linear comparison, refusing an inverted pair."""
    df1 = params_free - params_restricted
    df2 = n_rows - params_free
    if df1 <= 0 or df2 <= 0:
        raise ValueError("a nested comparison needs strictly more free parameters and spare rows")
    scale = max(abs(ssr_restricted), abs(ssr_free), 1e-300)
    if ssr_restricted < ssr_free - 1e-9 * scale:
        raise ValueError(
            "the restricted model fitted strictly better than the free one: the two designs are not "
            "nested, or the solve is broken"
        )
    gain = max(0.0, ssr_restricted - ssr_free)
    if ssr_free <= 0.0:
        raise ValueError("the free model left no residual; the design is saturated")
    statistic = (gain / df1) / (ssr_free / df2)
    return statistic, df1, df2, f_survival(statistic, df1, df2)


# ------------------------------------------------------------------------------------------------
# channels
# ------------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ChannelResult:
    """One detector channel: its statistic, its null, its p-value, and whether it fired."""

    channel: ChannelId
    statistic_name: str
    statistic: float
    p_value: float
    alpha: float
    fired: bool
    null_model: str
    detects: str
    detail: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.channel, ChannelId):
            raise TypeError("channel must be a ChannelId")
        for name in ("statistic_name", "null_model", "detects"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if not math.isfinite(self.statistic):
            raise ValueError("statistic must be finite")
        if not 0.0 < self.p_value <= 1.0:
            raise ValueError("p_value must lie in (0, 1]")
        _positive(self.alpha, what="alpha")
        if not isinstance(self.fired, bool):
            raise TypeError("fired must be bool")
        if self.fired != (self.p_value <= self.alpha):
            raise ValueError("fired must agree with the declared alpha; it is not a free field")
        object.__setattr__(self, "detail", dict(self.detail))

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able channel record."""
        return {
            "channel": self.channel.value,
            "statistic_name": self.statistic_name,
            "statistic": self.statistic,
            "p_value": self.p_value,
            "alpha": self.alpha,
            "fired": self.fired,
            "null_model": self.null_model,
            "detects": self.detects,
            "detail": dict(sorted(self.detail.items())),
        }


def _curvature_channel(
    conditions: Sequence[ConditionSeries], *, alpha: float
) -> tuple[ChannelResult, float, float, dict[str, np.ndarray]]:
    ssr_reduced, params_reduced, residuals = _per_condition_fit(conditions, degree=1)
    ssr_extended, params_extended, _ = _per_condition_fit(conditions, degree=2)
    rows = sum(condition.dimensionless_time.size for condition in conditions)
    statistic, df1, df2, p_value = _nested_f(
        ssr_restricted=ssr_reduced,
        ssr_free=ssr_extended,
        params_restricted=params_reduced,
        params_free=params_extended,
        n_rows=rows,
    )
    delta_bic = (
        rows * math.log(ssr_reduced / rows)
        + params_reduced * math.log(rows)
        - (rows * math.log(ssr_extended / rows) + params_extended * math.log(rows))
    )
    return (
        ChannelResult(
            channel=ChannelId.CURVATURE,
            statistic_name="nested F, one quadratic term per condition",
            statistic=statistic,
            p_value=p_value,
            alpha=alpha,
            fired=p_value <= alpha,
            null_model="the per-condition affine model is correct",
            detects=(
                "leftover pointwise shape after per-condition refitting; blind to any omission whose "
                "signature lies inside the fitted span"
            ),
            detail={
                "df_numerator": float(df1),
                "df_denominator": float(df2),
                "delta_bic_reduced_minus_extended": delta_bic,
                "ssr_reduced": ssr_reduced,
                "ssr_extended": ssr_extended,
            },
        ),
        ssr_reduced,
        float(params_reduced),
        residuals,
    )


def _mean_tau_over_dx(
    residuals: Mapping[str, np.ndarray], *, dx: float, sokal_window_c: float
) -> float:
    values = [
        integrated_autocorrelation_time(residual, dx, sokal_window_c=sokal_window_c)[0] / dx
        for residual in residuals.values()
    ]
    return float(np.mean(values))


def _serial_structure_channel(
    residuals: Mapping[str, np.ndarray],
    *,
    dx: float,
    contract: InadequacyContract,
    alpha: float,
    permutation_stream: np.random.SeedSequence | None,
) -> ChannelResult:
    observed = _mean_tau_over_dx(residuals, dx=dx, sokal_window_c=contract.sokal_window_c)
    stream = (
        np.random.SeedSequence(contract.permutation_seed)
        if permutation_stream is None
        else permutation_stream
    )
    rng = np.random.default_rng(stream)
    hits = 0
    for _ in range(contract.permutation_trials):
        shuffled = {key: rng.permutation(value) for key, value in residuals.items()}
        if _mean_tau_over_dx(shuffled, dx=dx, sokal_window_c=contract.sokal_window_c) >= observed:
            hits += 1
    p_value = (hits + 1.0) / (contract.permutation_trials + 1.0)
    return ChannelResult(
        channel=ChannelId.SERIAL_STRUCTURE,
        statistic_name="mean residual tau_int / dx",
        statistic=observed,
        p_value=p_value,
        alpha=alpha,
        fired=p_value <= alpha,
        null_model=(
            "the residual's ORDER carries no information: the null ensemble is permutations of the "
            "observed residuals, so it has the identical magnitude and differs only in structure"
        ),
        detects=(
            "a residual that is structured rather than large; this is the channel that separates "
            "'persistent residual' from 'noisy measurement'"
        ),
        detail={
            "permutation_trials": float(contract.permutation_trials),
            "permutation_seed": float(contract.permutation_seed),
            # The DECLARED seed above names the contract; this names the stream that actually ran, so
            # a replicate driven by a spawned SeedSequence is not recorded under the root's number.
            "permutation_stream_word0": float(_stream_fingerprint(stream)[0]),
            "hits": float(hits),
        },
    )


def _shared_parameter_channel(
    conditions: Sequence[ConditionSeries],
    *,
    ssr_per_condition: float,
    params_per_condition: float,
    alpha: float,
) -> tuple[ChannelResult, float, int, float]:
    ssr_shared, params_shared, shared_rate, _ = _shared_fit(conditions)
    rows = sum(condition.dimensionless_time.size for condition in conditions)
    statistic, df1, df2, p_value = _nested_f(
        ssr_restricted=ssr_shared,
        ssr_free=ssr_per_condition,
        params_restricted=params_shared,
        params_free=int(params_per_condition),
        n_rows=rows,
    )
    inflation = ssr_shared / ssr_per_condition
    return (
        ChannelResult(
            channel=ChannelId.SHARED_PARAMETER,
            statistic_name="nested F, shared rate vs per-condition rates",
            statistic=statistic,
            p_value=p_value,
            alpha=alpha,
            fired=p_value <= alpha,
            null_model="one rate law describes every condition; only the amplitude is per-condition",
            detects=(
                "an omission whose signature cannot be absorbed by all conditions AT ONCE. "
                "Per-condition refitting is exactly the operation that hides a missing law, so "
                "forbidding it is the strong test"
            ),
            detail={
                "df_numerator": float(df1),
                "df_denominator": float(df2),
                "ssr_shared": ssr_shared,
                "ssr_per_condition": ssr_per_condition,
                "misfit_inflation": inflation,
                "shared_rate": shared_rate,
            },
        ),
        ssr_shared,
        params_shared,
        shared_rate,
    )


def _held_out_channel(
    held_out: ConditionSeries,
    *,
    shared_rate: float,
    ssr_shared: float,
    dof_shared: int,
    alpha: float,
) -> tuple[ChannelResult, float, float]:
    x = held_out.dimensionless_time
    target = held_out.log_signal + shared_rate * x
    _, ssr_out, _ = _ols(np.ones((x.size, 1)), target)
    dof_out = x.size - 1
    in_sample_rms = math.sqrt(ssr_shared / dof_shared)
    held_out_rms = math.sqrt(ssr_out / dof_out)
    statistic = (ssr_out / dof_out) / (ssr_shared / dof_shared)
    p_value = f_survival(statistic, dof_out, dof_shared)
    return (
        ChannelResult(
            channel=ChannelId.HELD_OUT_INTERVENTION,
            statistic_name="held-out / in-sample residual variance ratio",
            statistic=statistic,
            p_value=p_value,
            alpha=alpha,
            fired=p_value <= alpha,
            null_model="the fitted rate law transfers to an intervention it was never fitted at",
            detects=(
                "a model that interpolates inside its fitted conditions and fails outside them; "
                "in-sample fit quality is not validation"
            ),
            detail={
                "in_sample_rms": in_sample_rms,
                "held_out_rms": held_out_rms,
                "error_ratio": held_out_rms / in_sample_rms,
                "held_out_intervention": held_out.intervention,
                "df_numerator": float(dof_out),
                "df_denominator": float(dof_shared),
            },
        ),
        in_sample_rms,
        held_out_rms,
    )


# ------------------------------------------------------------------------------------------------
# report
# ------------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ResidualStructureReport:
    """Everything measured on one experiment, including what did NOT fire."""

    experiment_id: str
    contract: InadequacyContract
    channels: tuple[ChannelResult, ...]
    per_channel_alpha: float
    pre_refit_rms: float
    post_refit_rms: float
    refit_reduction: float
    in_sample_rms: float
    held_out_rms: float
    held_out_error_ratio: float
    residual_interval: ESSInterval
    n_conditions: int
    n_samples_per_condition: int
    evidence_authority: str = EVIDENCE_AUTHORITY

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "experiment_id", _nonempty(self.experiment_id, what="experiment_id")
        )
        if not isinstance(self.contract, InadequacyContract):
            raise TypeError("contract must be an InadequacyContract")
        channels = tuple(self.channels)
        if not channels or any(not isinstance(item, ChannelResult) for item in channels):
            raise TypeError("channels must be a non-empty tuple of ChannelResult")
        if len({item.channel for item in channels}) != len(channels):
            raise ValueError("a channel may not appear twice")
        object.__setattr__(self, "channels", channels)
        if not isinstance(self.residual_interval, ESSInterval):
            raise TypeError("residual_interval must be an ESSInterval")

    @property
    def fired_channels(self) -> tuple[ChannelId, ...]:
        """Channels whose p-value met the per-channel level, in declaration order."""
        return tuple(item.channel for item in self.channels if item.fired)

    @property
    def any_fired(self) -> bool:
        """Whether any channel fired at all."""
        return bool(self.fired_channels)

    def channel(self, which: ChannelId) -> ChannelResult:
        """Return one channel's result.

        Args:
            which: The channel to look up.

        Returns:
            Its :class:`ChannelResult`.

        Raises:
            KeyError: If the channel did not run.
        """
        for item in self.channels:
            if item.channel is which:
                return item
        raise KeyError(f"channel {which.value} did not run")

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able report."""
        return {
            "experiment_id": self.experiment_id,
            "contract": self.contract.to_dict(),
            "per_channel_alpha": self.per_channel_alpha,
            "channels": [item.to_dict() for item in self.channels],
            "fired_channels": [item.value for item in self.fired_channels],
            "pre_refit_rms": self.pre_refit_rms,
            "post_refit_rms": self.post_refit_rms,
            "refit_reduction": self.refit_reduction,
            "in_sample_rms": self.in_sample_rms,
            "held_out_rms": self.held_out_rms,
            "held_out_error_ratio": self.held_out_error_ratio,
            "residual_interval": self.residual_interval.as_dict(),
            "n_conditions": self.n_conditions,
            "n_samples_per_condition": self.n_samples_per_condition,
            "evidence_authority": self.evidence_authority,
        }


# ------------------------------------------------------------------------------------------------
# candidate laws — pre-registered, before any run
# ------------------------------------------------------------------------------------------------

#: Which candidate law each firing pattern proposes.  Written down BEFORE the detector is run, so the
#: mapping from evidence to candidate is a pre-registration and not a post-hoc reading of a residual.
CANDIDATE_LAWS: Mapping[ChannelId, MissingLaw] = {
    ChannelId.SHARED_PARAMETER: MissingLaw(
        law_id="KU-S12.shared-rate-inconsistency",
        statement=(
            "A relaxation process is present whose rate does NOT respond to the intervention the "
            "reduced model attributes all of the response to. One rate law cannot describe every "
            "condition at once, so the omission is in the law, not in the parameter."
        ),
        blocked_by=(
            "the reduced model has one rate channel; representing a second, intervention-independent "
            "process requires a component with its own state, not another fitted constant"
        ),
        nearest_existing="single-rate relaxation with a per-condition amplitude nuisance",
    ),
    ChannelId.HELD_OUT_INTERVENTION: MissingLaw(
        law_id="KU-S12.intervention-transfer-failure",
        statement=(
            "The fitted rate law does not transfer to an intervention level it was not fitted at. "
            "The dependence of the response on the intervention is mis-specified, not mis-calibrated."
        ),
        blocked_by=(
            "no term in the reduced model carries an intervention dependence other than the one "
            "assumed proportional; a new dependence is a new law"
        ),
        nearest_existing="response rate assumed strictly proportional to the intervention",
    ),
    ChannelId.CURVATURE: MissingLaw(
        law_id="KU-S12.unmodelled-curvature",
        statement=(
            "A second relaxation mode is present at a rate separated enough from the represented one "
            "that a single exponential cannot absorb its curvature."
        ),
        blocked_by="the reduced model carries exactly one relaxation mode",
        nearest_existing="single-exponential relaxation",
    ),
    ChannelId.SERIAL_STRUCTURE: MissingLaw(
        law_id="KU-S12.structured-residual",
        statement=(
            "The residual carries serial structure that survives refitting: an unmodelled process "
            "with its own memory time, not observation noise."
        ),
        blocked_by="the reduced model treats every departure from its own trend as independent noise",
        nearest_existing="independent identically distributed observation noise",
    ),
}

_CHANNEL_EXPANSIONS: Mapping[ChannelId, tuple[ExpansionRoute, ...]] = {
    ChannelId.CURVATURE: (ExpansionRoute.MODEL_INADEQUACY, ExpansionRoute.MISSING_LAW),
    ChannelId.SERIAL_STRUCTURE: (
        ExpansionRoute.MODEL_INADEQUACY,
        ExpansionRoute.OBSERVATION_MODEL_EXTENSION,
    ),
    ChannelId.SHARED_PARAMETER: (ExpansionRoute.MISSING_LAW, ExpansionRoute.MISSING_COMPONENT),
    ChannelId.HELD_OUT_INTERVENTION: (
        ExpansionRoute.MISSING_LAW,
        ExpansionRoute.INTERVENTION_MODALITY,
    ),
}


# ------------------------------------------------------------------------------------------------
# verdicts
# ------------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MissingLawProposal:
    """A CANDIDATE missing law with its evidence.  It is never a decision, and it cannot become one.

    The master plan lists the model-inadequacy detector's physical authority as *제안만 가능* — it may
    propose only.  That is encoded rather than documented:

    * ``evidence_source`` is a property returning :data:`EvidenceSource.ANALYTIC_ORACLE`.  There is no
      field to overwrite and no argument through which a caller could spell a native provenance.
    * ``decided`` is a property returning ``False`` unconditionally.
    * :meth:`to_archetype_verdict` is the only route to a
      :class:`~aleph.virtual_cell.archetypes.MissingLawVerdict`, and it refuses unless an external
      adjudication reference is supplied that is not this proposal.
    """

    proposal_id: str
    question: str
    candidate_laws: tuple[MissingLaw, ...]
    evidence: ResidualStructureReport
    expansions: tuple[ExpansionRoute, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "proposal_id", _nonempty(self.proposal_id, what="proposal_id"))
        object.__setattr__(self, "question", _nonempty(self.question, what="question"))
        laws = tuple(self.candidate_laws)
        if not laws or any(not isinstance(law, MissingLaw) for law in laws):
            raise ValueError(
                "a proposal must NAME at least one candidate law; 'something is wrong' is not a "
                "candidate and cannot be scheduled"
            )
        if len({law.law_id for law in laws}) != len(laws):
            raise ValueError("candidate_laws must not repeat a law_id")
        object.__setattr__(self, "candidate_laws", laws)
        if not isinstance(self.evidence, ResidualStructureReport):
            raise TypeError("evidence must be a ResidualStructureReport")
        if not self.evidence.any_fired:
            raise ValueError(
                "a proposal cannot be raised on evidence in which no channel fired; that is a "
                "NoCandidateReport"
            )
        routes = tuple(self.expansions)
        if not routes or any(not isinstance(route, ExpansionRoute) for route in routes):
            raise ValueError("a proposal must name at least one expansion route")
        if len(set(routes)) != len(routes):
            raise ValueError("expansions must not repeat a route")
        object.__setattr__(self, "expansions", routes)

    @property
    def status(self) -> ProposalStatus:
        """Always ``CANDIDATE_PROPOSED``; this type has no other state."""
        return ProposalStatus.CANDIDATE_PROPOSED

    @property
    def evidence_source(self) -> EvidenceSource:
        """Always ``ANALYTIC_ORACLE``.  A property, so there is no field to overwrite."""
        return EvidenceSource.ANALYTIC_ORACLE

    @property
    def decided(self) -> bool:
        """Always ``False``.  An inadequacy detector proposes; it does not decide."""
        return False

    @property
    def law_ids(self) -> tuple[str, ...]:
        """The candidate law identifiers, in declaration order."""
        return tuple(law.law_id for law in self.candidate_laws)

    def to_archetype_verdict(
        self,
        *,
        archetype_id: str,
        ring: ArchetypeRing,
        adjudicated_by: str,
        unsupported_claims: Sequence[str],
    ) -> MissingLawVerdict:
        """Convert this proposal into an archetype-level refusal, after external adjudication.

        Args:
            archetype_id: The archetype the adjudicated verdict attaches to.
            ring: Its ring.  Ring A archetypes do not emit refusals.
            adjudicated_by: An external reference — a PI decision record, a review, a run record.  It
                may not name this proposal or this module.
            unsupported_claims: What may not be claimed while the law is missing.

        Returns:
            A :class:`~aleph.virtual_cell.archetypes.MissingLawVerdict`.

        Raises:
            SelfPromotionError: If ``adjudicated_by`` is empty, names this proposal, or names this
                module.  A detector that can promote its own proposal is a detector that decides.
        """
        reference = adjudicated_by.strip() if isinstance(adjudicated_by, str) else ""
        if not reference:
            raise SelfPromotionError(
                "a candidate missing law needs an external adjudication reference; the detector's own "
                "evidence is what raised the proposal and cannot also settle it"
            )
        lowered = reference.lower()
        forbidden = {self.proposal_id.lower(), __name__.lower(), "sandbox_inadequacy"}
        if any(token in lowered for token in forbidden) or lowered in {
            law.law_id.lower() for law in self.candidate_laws
        }:
            raise SelfPromotionError(
                f"adjudicated_by={reference!r} names this proposal or this module: self-promotion"
            )
        return MissingLawVerdict(
            archetype_id=archetype_id,
            ring=ring,
            missing_components=(),
            missing_laws=self.candidate_laws,
            unsupported_claims=tuple(unsupported_claims),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able proposal."""
        return {
            "proposal_id": self.proposal_id,
            "status": self.status.value,
            "decided": self.decided,
            "evidence_source": self.evidence_source.value,
            "question": self.question,
            "candidate_laws": [law.to_dict() for law in self.candidate_laws],
            "expansions": [route.value for route in self.expansions],
            "evidence": self.evidence.to_dict(),
        }


@dataclass(frozen=True, slots=True)
class NoCandidateReport:
    """No channel fired.  A separate type, so silence can never be read as a proposal.

    ``silence_is_not_absence`` is carried explicitly because this is the report a masked omission
    produces, and the entire lane exists because that report was historically read as "the model is
    fine".
    """

    report_id: str
    evidence: ResidualStructureReport
    silence_is_not_absence: str = (
        "no channel fired. This bounds what these channels can see at this design and sample size; "
        "it is NOT evidence that no law is missing. An omission degenerate with a fitted parameter "
        "produces exactly this report at every amplitude."
    )

    def __post_init__(self) -> None:
        object.__setattr__(self, "report_id", _nonempty(self.report_id, what="report_id"))
        if not isinstance(self.evidence, ResidualStructureReport):
            raise TypeError("evidence must be a ResidualStructureReport")
        if self.evidence.any_fired:
            raise ValueError("channels fired; this evidence belongs to a MissingLawProposal")

    @property
    def status(self) -> ProposalStatus:
        """Always ``NO_CANDIDATE``."""
        return ProposalStatus.NO_CANDIDATE

    @property
    def evidence_source(self) -> EvidenceSource:
        """Always ``ANALYTIC_ORACLE``."""
        return EvidenceSource.ANALYTIC_ORACLE

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able report."""
        return {
            "report_id": self.report_id,
            "status": self.status.value,
            "evidence_source": self.evidence_source.value,
            "silence_is_not_absence": self.silence_is_not_absence,
            "evidence": self.evidence.to_dict(),
        }


# ------------------------------------------------------------------------------------------------
# the detector
# ------------------------------------------------------------------------------------------------


def analyse_conditions(
    conditions: Sequence[ConditionSeries],
    held_out: ConditionSeries,
    *,
    contract: InadequacyContract,
    experiment_id: str,
    permutation_stream: np.random.SeedSequence | None = None,
) -> ResidualStructureReport:
    """Run all four channels and return the measured report, fired or not.

    Args:
        conditions: The training interventions.
        held_out: One intervention never fitted.
        contract: Pre-declared thresholds.
        experiment_id: Recorded in the report.
        permutation_stream: Optional explicit stream for the D2 permutation ensemble.  ``None`` uses
            ``contract.permutation_seed``, which is the declared single-experiment behaviour.  A
            replicated null control passes a spawned ``SeedSequence`` here so the permutation stream
            is provably disjoint from the streams that generated the data, instead of being an
            integer offset that can land on one of them.

    Returns:
        The :class:`ResidualStructureReport`.

    Raises:
        TypeError: If any argument has the wrong type.
        ValueError: If there are too few conditions, the series are too short or unequal in length,
            the held-out intervention duplicates a training one, or the grids differ.
    """
    if not isinstance(contract, InadequacyContract):
        raise TypeError("contract must be an InadequacyContract")
    if permutation_stream is not None and not isinstance(
        permutation_stream, np.random.SeedSequence
    ):
        raise TypeError("permutation_stream must be a numpy.random.SeedSequence or None")
    training = tuple(conditions)
    if any(not isinstance(item, ConditionSeries) for item in training):
        raise TypeError("conditions must be ConditionSeries")
    if not isinstance(held_out, ConditionSeries):
        raise TypeError("held_out must be a ConditionSeries")
    if len(training) < contract.min_conditions:
        raise ValueError(
            f"need at least {contract.min_conditions} training conditions; got {len(training)}"
        )
    sizes = {item.dimensionless_time.size for item in training}
    if len(sizes) != 1:
        raise ValueError("every training condition must have the same sample count")
    size = sizes.pop()
    if size < contract.min_samples_per_condition:
        raise ValueError(
            f"need at least {contract.min_samples_per_condition} samples per condition; got {size}"
        )
    if held_out.dimensionless_time.size != size:
        raise ValueError("the held-out condition must be sampled on the same grid")
    if not np.allclose(held_out.dimensionless_time, training[0].dimensionless_time):
        raise ValueError("the held-out condition must share the training dimensionless clock")
    levels = [item.intervention for item in training]
    if len(set(levels)) != len(levels):
        raise ValueError("training interventions must be distinct")
    if held_out.intervention in set(levels):
        raise ValueError(
            "the held-out intervention is one of the training interventions; that is an in-sample "
            "prediction wearing a held-out label"
        )

    alpha = contract.per_channel_alpha(len(ChannelId))
    dx = training[0].dx

    curvature, ssr_per_condition, params_per_condition, residuals = _curvature_channel(
        training, alpha=alpha
    )
    serial = _serial_structure_channel(
        residuals,
        dx=dx,
        contract=contract,
        alpha=alpha,
        permutation_stream=permutation_stream,
    )
    shared, ssr_shared, params_shared, shared_rate = _shared_parameter_channel(
        training,
        ssr_per_condition=ssr_per_condition,
        params_per_condition=params_per_condition,
        alpha=alpha,
    )
    rows = size * len(training)
    transfer, in_sample_rms, held_out_rms = _held_out_channel(
        held_out,
        shared_rate=shared_rate,
        ssr_shared=ssr_shared,
        dof_shared=rows - params_shared,
        alpha=alpha,
    )

    pooled = np.concatenate([residuals[item.condition_id] for item in training])
    interval = _ess_interval(
        pooled,
        dx,
        observable="pooled per-condition residual",
        sokal_window_c=contract.sokal_window_c,
        coverage_z=contract.coverage_z,
    )
    pre_refit = math.sqrt(
        float(
            np.mean(
                np.concatenate(
                    [(item.log_signal + item.dimensionless_time) ** 2 for item in training]
                )
            )
        )
    )
    post_refit = math.sqrt(ssr_per_condition / (rows - params_per_condition))

    return ResidualStructureReport(
        experiment_id=experiment_id,
        contract=contract,
        channels=(curvature, serial, shared, transfer),
        per_channel_alpha=alpha,
        pre_refit_rms=pre_refit,
        post_refit_rms=post_refit,
        refit_reduction=pre_refit / post_refit,
        in_sample_rms=in_sample_rms,
        held_out_rms=held_out_rms,
        held_out_error_ratio=held_out_rms / in_sample_rms,
        residual_interval=interval,
        n_conditions=len(training),
        n_samples_per_condition=size,
    )


def detect_missing_law(
    conditions: Sequence[ConditionSeries],
    held_out: ConditionSeries,
    *,
    contract: InadequacyContract,
    experiment_id: str,
) -> MissingLawProposal | NoCandidateReport:
    """Run the detector and return a CANDIDATE proposal, or an explicit no-candidate report.

    The composite rule is declared here and nowhere else: a candidate is proposed when at least one
    channel's p-value meets the Bonferroni per-channel level.  **Residual magnitude is not a
    criterion.**  It cannot be: refitting shrinks it, and inflating the observation noise raises it
    without any law being missing.  Every channel is a ratio or a permutation test and is therefore
    exactly invariant under rescaling the noise.

    Args:
        conditions: Training interventions.
        held_out: One intervention never fitted.
        contract: Pre-declared thresholds.
        experiment_id: Recorded in the report and used to build the proposal identifier.

    Returns:
        A :class:`MissingLawProposal` when a channel fired, otherwise a :class:`NoCandidateReport`.
    """
    report = analyse_conditions(
        conditions, held_out, contract=contract, experiment_id=experiment_id
    )
    fired = report.fired_channels
    if not fired:
        return NoCandidateReport(report_id=f"{report.experiment_id}/no-candidate", evidence=report)
    laws = tuple(CANDIDATE_LAWS[channel] for channel in fired)
    routes: list[ExpansionRoute] = []
    for channel in fired:
        for route in _CHANNEL_EXPANSIONS[channel]:
            if route not in routes:
                routes.append(route)
    return MissingLawProposal(
        proposal_id=f"{report.experiment_id}/candidate",
        question=(
            "which constitutive law, absent from the reduced model, would produce this residual "
            "structure?"
        ),
        candidate_laws=laws,
        evidence=report,
        expansions=tuple(routes),
    )


# ------------------------------------------------------------------------------------------------
# null control
# ------------------------------------------------------------------------------------------------

#: The roles a null-control replicate needs, in the order they are spawned.  Two roles that share a
#: stream are not two draws from the null; they are one draw used twice.
REPLICATE_STREAM_ROLES: tuple[str, ...] = ("training", "held-out", "permutation")

#: How :func:`derive_replicate_streams` builds its streams, recorded in every report so the claim of
#: independence is checkable from the artifact rather than taken on the docstring's word.
STREAM_DERIVATION = (
    "numpy.random.SeedSequence(seed).spawn(trials)[i].spawn(3) -> "
    "(training, held-out, permutation); spawn keys (i, 0) / (i, 1) / (i, 2) are pairwise distinct "
    "by construction, so no two (replicate, role) pairs address the same stream"
)


@dataclass(frozen=True, slots=True)
class ReplicateStreams:
    """The three PROVABLY distinct RNG streams one null-control replicate needs.

    WHY this type exists, stated as the counterexample it was written against.  The previous
    derivation was arithmetic on integers: replicate ``i`` drew its training data from ``seed +
    1009 i``, its held-out condition from ``seed + 1009 i + 1``, and its D2 permutation ensemble from
    ``permutation_seed + i``.  With the committed test contract, ``seed == permutation_seed ==
    20260729``, so replicate 0's TRAINING stream and replicate 0's PERMUTATION stream were the same
    stream, and replicate 0's held-out stream ``20260730`` was replicate 1's permutation stream.  Two
    generators built from one seed are one stream, not two; the binomial reading of the measured
    false-positive rate assumes independent replicates, and that premise did not hold in the form the
    code asserted.  ``SeedSequence`` spawning removes the arithmetic: the streams are addressed by
    spawn key, and distinct spawn keys cannot collide however the root seeds are chosen.
    """

    replicate: int
    training: np.random.SeedSequence
    held_out: np.random.SeedSequence
    permutation: np.random.SeedSequence

    def __post_init__(self) -> None:
        if isinstance(self.replicate, bool) or not isinstance(self.replicate, int):
            raise TypeError("replicate must be an int")
        if self.replicate < 0:
            raise ValueError("replicate must be nonnegative")
        for name in REPLICATE_STREAM_ROLES:
            attribute = name.replace("-", "_")
            if not isinstance(getattr(self, attribute), np.random.SeedSequence):
                raise TypeError(f"{attribute} must be a numpy.random.SeedSequence")

    def stream(self, role: str) -> np.random.SeedSequence:
        """Return the stream for ``role``.

        Args:
            role: One of :data:`REPLICATE_STREAM_ROLES`.

        Returns:
            The ``SeedSequence`` addressed by that role.

        Raises:
            KeyError: If ``role`` is not a declared role.
        """
        if role not in REPLICATE_STREAM_ROLES:
            raise KeyError(
                f"unknown stream role {role!r}; expected one of {REPLICATE_STREAM_ROLES}"
            )
        return getattr(self, role.replace("-", "_"))

    def fingerprints(self) -> dict[str, tuple[int, ...]]:
        """Return 128 bits of each role's generator state, keyed by role.

        This is the quantity a caller compares across replicates to CHECK independence rather than
        assume it: equal fingerprints would mean equal streams.
        """
        return {role: _stream_fingerprint(self.stream(role)) for role in REPLICATE_STREAM_ROLES}


def derive_replicate_streams(*, seed: int, trials: int) -> tuple[ReplicateStreams, ...]:
    """Spawn one independent stream per (replicate, role) pair from a single recorded seed.

    Args:
        seed: Root entropy, recorded in the report so the whole ensemble is reproducible.
        trials: Number of replicates.

    Returns:
        ``trials`` :class:`ReplicateStreams`, one per replicate.

    Raises:
        ValueError: If ``seed`` is not a nonnegative int or ``trials`` is not positive.
    """
    count = _positive_int(trials, what="trials")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a nonnegative int")
    root = np.random.SeedSequence(seed)
    streams: list[ReplicateStreams] = []
    for index, child in enumerate(root.spawn(count)):
        training, held_out, permutation = child.spawn(len(REPLICATE_STREAM_ROLES))
        streams.append(
            ReplicateStreams(
                replicate=index,
                training=training,
                held_out=held_out,
                permutation=permutation,
            )
        )
    return tuple(streams)


@dataclass(frozen=True, slots=True)
class NullControlReport:
    """Measured false-positive rate of the detector when NOTHING is missing.

    A detector that always finds a missing law is useless, so this is not an optional extra.  The
    report also records which noise model was simulated: when it differs from the contract's assumed
    model, ``assumption_violated`` is true and the measured rate is the COST of that assumption, not a
    verdict on the detector's design.

    ``binomial_standard_error`` is only meaningful if the replicates are independent, so the report
    also carries the stream derivation and the number of distinct (replicate, role) streams that the
    measurement actually used.  The independence premise is thereby part of the artifact and can be
    audited, which is what it was not when the derivation lived only in a docstring sentence.
    """

    trials: int
    seed: int
    noise_model: NoiseModel
    assumed_noise_model: NoiseModel
    family_alpha: float
    false_positives: int
    false_positive_rate: float
    binomial_standard_error: float
    sigma_multiplier: float
    upper_bound: float
    within_declared_level: bool
    per_channel_false_positives: Mapping[str, int]
    stream_derivation: str
    distinct_streams: int

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "per_channel_false_positives", dict(self.per_channel_false_positives)
        )
        object.__setattr__(
            self, "stream_derivation", _nonempty(self.stream_derivation, what="stream_derivation")
        )
        _positive_int(self.distinct_streams, what="distinct_streams")
        if self.distinct_streams != self.trials * len(REPLICATE_STREAM_ROLES):
            raise ValueError(
                "distinct_streams must equal trials x roles; a smaller count means two roles shared "
                "a stream and the binomial reading of the rate does not hold"
            )

    @property
    def assumption_violated(self) -> bool:
        """Whether the simulated noise model differs from the one the channels assume."""
        return self.noise_model is not self.assumed_noise_model

    @property
    def inflation_over_declared(self) -> float:
        """Measured rate divided by the declared family-wise level."""
        return self.false_positive_rate / self.family_alpha

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able null-control record."""
        return {
            "trials": self.trials,
            "seed": self.seed,
            "noise_model": self.noise_model.value,
            "assumed_noise_model": self.assumed_noise_model.value,
            "assumption_violated": self.assumption_violated,
            "family_alpha": self.family_alpha,
            "false_positives": self.false_positives,
            "false_positive_rate": self.false_positive_rate,
            "binomial_standard_error": self.binomial_standard_error,
            "sigma_multiplier": self.sigma_multiplier,
            "upper_bound": self.upper_bound,
            "within_declared_level": self.within_declared_level,
            "inflation_over_declared": self.inflation_over_declared,
            "per_channel_false_positives": dict(sorted(self.per_channel_false_positives.items())),
            "stream_derivation": self.stream_derivation,
            "distinct_streams": self.distinct_streams,
            "stream_roles": list(REPLICATE_STREAM_ROLES),
        }


def measure_null_false_positive_rate(
    *,
    contract: InadequacyContract,
    trials: int,
    seed: int,
    interventions: Sequence[float],
    held_out_intervention: float,
    noise_sigma: float,
    n_samples: int,
    noise_model: NoiseModel = NoiseModel.HOMOSCEDASTIC_GAUSSIAN,
) -> NullControlReport:
    """Run the detector ``trials`` times on data with NO missing law and count the firings.

    Args:
        contract: Pre-declared thresholds.  Its ``permutation_seed`` still names the declared
            single-experiment stream and is recorded, but the replicates do NOT run on offsets of it.
        trials: Independent replicate experiments.
        seed: Root entropy.  Every stream is spawned from it by :func:`derive_replicate_streams`, so
            the ``3 x trials`` (replicate, role) streams are pairwise distinct by spawn key.  The
            integer-offset rule this replaced is the defect recorded on :class:`ReplicateStreams`:
            with ``seed == contract.permutation_seed == 20260729`` it gave replicate 0 the same
            stream for its training data and its D2 permutation ensemble, and gave replicate 1's
            permutation ensemble replicate 0's held-out stream.
        interventions: Training intervention levels.
        held_out_intervention: The held-out level.
        noise_sigma: Observation-noise scale.
        n_samples: Samples per condition.
        noise_model: Which noise process to simulate.  Passing anything but the contract's assumed
            model measures the cost of that assumption.

    Returns:
        The :class:`NullControlReport`.

    Raises:
        ValueError: If ``trials`` is not positive or the design is malformed.
    """
    count = _positive_int(trials, what="trials")
    streams = derive_replicate_streams(seed=seed, trials=count)
    per_channel = {item.value: 0 for item in ChannelId}
    firings = 0
    for index, replicate in enumerate(streams):
        training = build_conditions(
            interventions=interventions,
            omission=OmissionKind.NONE,
            separation=2.0,
            deviation_sigma_multiple=0.0,
            noise_sigma=noise_sigma,
            n_samples=n_samples,
            seed=replicate.training,
            noise_model=noise_model,
        )
        outside = build_conditions(
            interventions=(held_out_intervention,),
            omission=OmissionKind.NONE,
            separation=2.0,
            deviation_sigma_multiple=0.0,
            noise_sigma=noise_sigma,
            n_samples=n_samples,
            seed=replicate.held_out,
            noise_model=noise_model,
        )[0]
        report = analyse_conditions(
            training,
            outside,
            contract=contract,
            experiment_id=f"null-control/{index}",
            permutation_stream=replicate.permutation,
        )
        for channel in report.fired_channels:
            per_channel[channel.value] += 1
        firings += int(report.any_fired)

    rate = firings / count
    alpha = contract.family_alpha
    standard_error = math.sqrt(alpha * (1.0 - alpha) / count)
    upper = alpha + contract.false_positive_sigma * standard_error
    return NullControlReport(
        trials=count,
        seed=seed,
        noise_model=noise_model,
        assumed_noise_model=contract.assumed_noise_model,
        family_alpha=alpha,
        false_positives=firings,
        false_positive_rate=rate,
        binomial_standard_error=standard_error,
        sigma_multiplier=contract.false_positive_sigma,
        upper_bound=upper,
        within_declared_level=rate <= upper,
        per_channel_false_positives=per_channel,
        stream_derivation=STREAM_DERIVATION,
        distinct_streams=count * len(REPLICATE_STREAM_ROLES),
    )


# ------------------------------------------------------------------------------------------------
# pre-registration
# ------------------------------------------------------------------------------------------------

S12_MANIFEST: CellStateManifest = CellStateManifest(
    manifest_id="s12-inadequacy-oracle",
    species="none",
    lineage="analytic-oracle",
    identity={"system": "two-mode log-relaxation under a rate intervention"},
    biological_state={"note": "no cell is represented; this manifest exists to bind the card"},
    environment={"intervention_axis": "rate multiplier c"},
    components=("represented-relaxation-mode",),
    connectors=(),
    observations=("log-signal vs dimensionless time",),
    missing_components=("omitted-relaxation-mode",),
    missing_laws=(
        "KU-S12.shared-rate-inconsistency",
        "KU-S12.intervention-transfer-failure",
        "KU-S12.unmodelled-curvature",
        "KU-S12.structured-residual",
    ),
    provenance=("constructed-in-closed-form",),
)

SANDBOX_S12_CARD: SandboxExperimentCard = SandboxExperimentCard(
    experiment_id="S12-model-inadequacy-challenge",
    question=(
        "When a constitutive law is deliberately removed from the model, does the residual identify "
        "the omission — or does refitting a parameter hide it?"
    ),
    manifest_hashes=(S12_MANIFEST.manifest_hash,),
    intervention=(
        "vary the rate multiplier c over three fitted levels and one held-out level, while the truth "
        "carries a second relaxation mode the fitted model does not"
    ),
    observables=(
        "pre-refit residual RMS",
        "post-refit residual RMS",
        "nested-F curvature p-value and delta-BIC",
        "residual integrated autocorrelation time against a permutation null",
        "shared-vs-per-condition misfit inflation",
        "held-out / in-sample residual variance ratio",
    ),
    positive_controls=(
        "ABSOLUTE_TIME omission: the omitted mode does not co-scale with the intervention, so the "
        "shared-parameter and held-out channels must fire",
    ),
    negative_controls=(
        "no omission at all: the measured family-wise false-positive rate must sit at or below the "
        "declared level over a recorded seed and trial count",
    ),
    adversarial_controls=(
        "noise inflated with nothing missing: a magnitude-based detector fires, this one must not",
        "heavy-tailed and heteroscedastic noise with nothing missing: the channels' assumed noise "
        "model is violated and the inflation of the false-positive rate is MEASURED, not assumed away",
        "CO_SCALING omission at matched amplitude: the omitted law is degenerate with the fitted "
        "rate, and the shared-parameter and held-out channels must stay silent at every amplitude",
    ),
    held_out_interventions=(
        "one rate multiplier outside the fitted range, fitted only for its amplitude nuisance",
    ),
    falsifier=(
        "a parameter shift hides the missing law: if per-condition refitting drops every channel "
        "below its declared level while a law is in fact missing, the detector FAILED for that case "
        "and the failure is reported as a failure"
    ),
    failure_expansions=(
        "sim-data mismatch -> model inadequacy OR observation-model extension",
        "archetype failure -> missing component / missing law project",
        "image nonidentifiability -> time / force / intervention modality",
    ),
    budget_class="cpu-analytic-oracle",
    wet_lab=False,
)


def s12_evidence_graph() -> EvidenceGraph:
    """Return the provenance DAG for this lane's artifacts.

    Every node is a fixture or a computation.  None has an empirical origin, so
    :meth:`~aleph.virtual_cell.source_registry.EvidenceGraph.citable_as_biological_evidence` is
    empty by construction — which is the correct statement about a lane whose systems were written
    down rather than measured.

    Returns:
        The frozen :class:`~aleph.virtual_cell.source_registry.EvidenceGraph`.
    """
    registry = SourceRegistry()
    registry.register(
        SourceRecord(
            source_id="s12-oracle-system",
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.COMPUTATION,
            status=SourceStatus.FIXTURE,
            description=(
                "two-mode log-relaxation truth constructed in closed form, with the second mode "
                "removed from the fitted model"
            ),
            parents=("s12-protocol-constants",),
        )
    )
    registry.register(
        SourceRecord(
            source_id="s12-protocol-constants",
            provenance=ProvenanceKind.CONVENIENCE,
            origin=SourceOrigin.EXPERT_JUDGEMENT,
            status=SourceStatus.FIXTURE,
            description=(
                "reference rate, omitted-mode weight and observation window of the constructed "
                "system; they define the truth and gate nothing"
            ),
        )
    )
    registry.register(
        SourceRecord(
            source_id="s12-detector-report",
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.COMPUTATION,
            status=SourceStatus.FIXTURE,
            description="channel statistics and p-values measured on the constructed system",
            parents=("s12-oracle-system",),
        )
    )
    registry.register(
        SourceRecord(
            source_id="s12-missing-law-candidate",
            provenance=ProvenanceKind.INFERRED_POSTERIOR,
            origin=SourceOrigin.COMPUTATION,
            status=SourceStatus.FIXTURE,
            description=(
                "candidate missing law proposed from the residual structure; a proposal awaiting "
                "external adjudication, never a verdict"
            ),
            parents=("s12-detector-report",),
        )
    )
    return registry.freeze()


def experiment_fingerprint(report: ResidualStructureReport) -> str:
    """Return a SHA-256 over the report's canonical content, so a claim can be pinned to a run.

    Args:
        report: The measured report.

    Returns:
        A lowercase SHA-256 hex digest.
    """
    payload = json.dumps(
        report.to_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode()
    return hashlib.sha256(payload).hexdigest()
