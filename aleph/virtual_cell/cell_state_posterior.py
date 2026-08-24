"""`CellStatePosterior` — the object master plan section 0 declares as the project's data unit.

The plan's final data unit is not a cell-line name and not a parameter value:

```text
CellStatePosterior
  = component/connector graph
  + population and geometry distribution
  + molecular-parameter posterior
  + accepted dynamical regime
  + intervention-response distribution
  + observation operators
  + evidence provenance
  + known model inadequacy
```

The adversarial audit found mechanically (`audit_checklist` item A13) that this type did not exist,
and diagnosed the consequence exactly: eleven lanes of good parts that are **not yet a circuit**.
This module is the wiring.  Every field below BINDS a type another lane already built rather than
introducing a parallel one — the graph is
:class:`~aleph.virtual_cell.grammar.EngineGrammar`, the priors are
:class:`~aleph.virtual_cell.archetypes.PriorDistribution` (carrying ``PI_GAP`` where nothing is
sourced), the regime is :class:`~aleph.virtual_cell.regime.RegimeReport`, the provenance is
:class:`~aleph.virtual_cell.source_registry.EvidenceGraph`, the model inadequacy is a
:class:`~aleph.virtual_cell.sandbox_inadequacy.MissingLawProposal`, and the observation operators
are :class:`~aleph.virtual_cell.observation_operator.ObservationOperator` implementations.  A field
that invented its own type here would have re-opened the drift the audit exists to catch.

Four design decisions, and the failure each one refuses
------------------------------------------------------
**1. The representation is a declared field, not an assumption.**  A posterior may be a point, an
explicit ensemble, a local Gaussian, a Gaussian mixture, or a tensor-train / latent density, and it
must say WHICH.  ``EXPLICIT_STATE`` — a point — is a legitimate answer and is what the engine gives
today; it is recorded as a point and given no wider interval.  The reverse is also refused: a
posterior declaring ``GAUSSIAN_MIXTURE`` whose budgets carry no mixture, or ``LOCAL_GAUSSIAN`` whose
budgets do, is rejected, because a representation label that does not describe the content is exactly
the over-claim the audit's A2/R2 rows are about.  For the approximate kinds
:class:`~aleph.virtual_cell.contracts.RepresentationDescriptor` already forces a MEASURED
:class:`~aleph.virtual_cell.contracts.ApproximationErrorReport`; this module additionally fails
closed against the caller's declared ``max_representation_error``.

**2. Uncertainty stays decomposed by source.**  Master plan section 2.2 and audit A1 both forbid
collapsing thermal, molecular-event, quenched-structural, cell-to-cell, measurement,
parameter-epistemic and numerical uncertainty into one ``Sigma``.  That is enforced structurally:
:class:`UncertaintyBudget` has **no single-covariance constructor** — it cannot be built from a total
at all — and :class:`~aleph.virtual_cell.contracts.UncertaintySource` deliberately has no ``TOTAL``
member, so a lumped width has no name to be spelled under.  Every total is computed ON DEMAND from
the components, which remain available afterwards.  :meth:`UncertaintyBudget.collapse_to_gaussian`
exists, is explicitly lossy, MEASURES its own error, and REFUSES when that error exceeds the caller's
declared tolerance.  The error it measures is the one A1 already registered: moment-matching a
two-source mixture gets the variance exactly right and misreports the fourth moment by a factor of
1.9608 for the registered ``sigma = 0.1`` / ``sigma = 1.0`` pair.  Note the contrast the same code
measures: for sources that combine by CONVOLUTION rather than as a mixture, the collapse is exact and
the measured error is 0.  Merging is not always wrong — it is wrong in a specific, computable way,
and this refuses it exactly there.

**3. No point estimate along a direction that is not identifiable.**  Audit A5 rejects
image-to-parameter point prediction, and :mod:`aleph.virtual_cell.sandbox_inverse` already
implements the refusal.  This composes with it — :meth:`MolecularParameterPosterior.point_estimate`
raises the SAME :class:`~aleph.virtual_cell.sandbox_inverse.UnidentifiableDirectionError`, and
:meth:`IdentifiabilityDeclaration.from_report` reads the fractions straight out of an
:class:`~aleph.virtual_cell.sandbox_inverse.IdentifiabilityReport`.  It also fails CLOSED: a
parameter whose identifiability was never measured is refused too, because "nobody looked" and "it is
identifiable" are different statements.  Credible intervals stay available throughout — the refusal
is of the point estimate, never of the inference — and they are computed EXACTLY from the mixture
CDF, so an honest interval never has to route through the lossy collapse.

**4. It cannot promote itself.**  ``native_accepted`` and ``evidence_authority`` are properties over
``__slots__`` with no field behind them, following
:class:`~aleph.virtual_cell.surrogate.SurrogatePrediction`: there is nothing for
``dataclasses.replace`` to set.  That alone was NOT enough, and the over-claim is recorded here
rather than quietly repaired: an external review (2026-07-29, finding 4) subclassed the type,
overrode ``native_accepted`` to ``True``, passed ``isinstance``, and serialised
``to_dict()["native_accepted"] == True``.  Both routes are closed now —
:meth:`CellStatePosterior.__init_subclass__` refuses a subclass that redefines either property, and
every serialisation re-checks that the resolved property is still this class's own, which catches an
override installed after the class was created.  The envelope refuses
``NATIVE_ACCEPTED`` independently.  A :class:`RegimeReport` may only be carried when an
:class:`~aleph.virtual_cell.admissibility.AdmissibilityVerdict` says the steps it was measured on
were admissible AND transaction-sound: a regime read off inadmissible steps is a statement about the
solver, and the plan's section 2.1 makes accepted-step wiring the precondition for everything above
it.  A :class:`MissingLawProposal` stays a proposal; converting one into a verdict still requires
external adjudication in its own module.

Every threshold lives in :class:`PosteriorEvidenceContract`, which has **no defaults on any field**,
matching :class:`~aleph.virtual_cell.regime.RegimeEvidenceContract` for the same reason: a threshold
with a default is a threshold nobody declared.

Sanity Gate (recorded before first execution):

* **Dimensional** — every parameter carries its unit; a budget's unit must equal its parameter's, and
  variances are in that unit SQUARED (the module never adds a variance to a location).  Intervals and
  point estimates come back in the parameter's own unit.  Probabilities (levels, weights,
  identifiable fractions) are dimensionless and range-checked.
* **Boundary** — an empty component list, an empty parameter list, mixture weights that do not sum to
  one, a negative variance, a level outside ``(0, 1)``, a fraction outside ``(0, 1]`` and a
  representation kind outside :data:`POSTERIOR_REPRESENTATION_KINDS` all raise.  A zero-variance
  budget yields a DEGENERATE interval ``[x, x]`` rather than an error: for a point posterior that is
  the correct, visibly useless answer.
* **Conservation / invariant** — the decomposition is conserved through every operation:
  ``sum(variance_by_source) == total_variance`` for convolution-only budgets to round-off, and the
  components survive a collapse (:attr:`GaussianCollapse.budget` still holds them).  The exact
  interval satisfies ``cdf(hi) - cdf(lo) == level`` to the bisection tolerance, and the collapse is
  variance-preserving BY CONSTRUCTION — which is precisely why the variance cannot detect it and the
  fourth moment can.
* **Measurement protocol** — the collapse error is ``|E[x^4] / (3 var^2) - 1|``, a dimensionless
  relative excess-kurtosis error, compared against a tolerance the caller declared before the call and
  recorded in the returned :class:`~aleph.virtual_cell.contracts.ApproximationErrorReport` with
  content hashes of both the source budget and the reduced Gaussian.  A posterior's own identity is
  ``record_hash``, the SHA-256 of :meth:`CellStatePosterior.to_dict`, and that record covers EVERY
  constructor field — including the evidence graph's own content address and the accepted-step
  verdict's identity rather than its verdict bit, both of which it omitted until an external review
  (2026-07-29, finding 1) showed two posteriors with different provenance DAGs and different attempted
  step indices hashing identically.
* **Numerical** — the interval is obtained by bisecting a monotone mixture CDF built from
  :func:`math.erf`; 200 halvings of a bracket of 40 total standard deviations reach the float
  resolution of the answer, and the bracket is verified to contain the quantile before bisecting.  No
  matrix is inverted and no covariance is ever materialised, so there is no conditioning question to
  answer.

Nothing here imports Warp, touches a GPU, runs physics, or grants an evidence rung.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from types import MappingProxyType
from typing import Any

from aleph.virtual_cell.admissibility import AdmissibilityVerdict
from aleph.virtual_cell.archetypes import ArchetypeBuild, PriorDistribution, PriorStatus
from aleph.virtual_cell.contracts import (
    ApproximationErrorReport,
    ArtifactEnvelope,
    CellStateManifest,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    UncertaintySource,
)
from aleph.virtual_cell.grammar import ConnectorInstance, EngineGrammar, PopulationWiring
from aleph.virtual_cell.observation_operator import (
    MechanisticState,
    ObservationOperator,
    ObservationResult,
    validate_operator_registration,
)
from aleph.virtual_cell.regime import EVIDENCE_AUTHORITY, RegimeReport
from aleph.virtual_cell.sandbox_inadequacy import MissingLawProposal
from aleph.virtual_cell.sandbox_inverse import (
    PARAMETER_NAMES,
    IdentifiabilityReport,
    UnidentifiableDirectionError,
)
from aleph.virtual_cell.source_registry import EvidenceGraph

__all__ = [
    "POINT_REPRESENTATION",
    "POSTERIOR_REPRESENTATION_KINDS",
    "CellStatePosterior",
    "GaussianCollapse",
    "IdentifiabilityDeclaration",
    "InterventionResponse",
    "LossyMergeRefused",
    "MolecularParameterPosterior",
    "PosteriorEvidenceContract",
    "PosteriorGraph",
    "UncertaintyBudget",
    "UncertaintyCombination",
    "UncertaintyComponent",
]

#: The representation a deterministic accepted-step run gives today.  It is a legitimate posterior
#: representation and is deliberately NOT dressed up as a distribution: a point posterior's credible
#: interval is a point, and that visible uselessness is the honest report.
POINT_REPRESENTATION = RepresentationKind.EXPLICIT_STATE

#: Representation kinds a cell-state posterior may declare.  ``SURROGATE``, ``SYNTHETIC_IMAGE`` and
#: ``SCALAR_PROJECTION`` are excluded on purpose: they are things you can DERIVE from a posterior, not
#: ways of representing one, and allowing them would let a surrogate's output be filed as the state.
POSTERIOR_REPRESENTATION_KINDS = frozenset(
    {
        RepresentationKind.EXPLICIT_STATE,
        RepresentationKind.EMPIRICAL_ENSEMBLE,
        RepresentationKind.LOCAL_GAUSSIAN,
        RepresentationKind.GAUSSIAN_MIXTURE,
        RepresentationKind.TENSOR_TRAIN,
        RepresentationKind.LATENT_DENSITY,
    }
)

_SHA256_CHARS = frozenset("0123456789abcdef")
#: An identifiable fraction is ``sum_k V[j, k]^2`` over an orthonormal eigenvector row, so it equals
#: one exactly in exact arithmetic and misses by a few ULP in floating point.  Values inside this
#: window are clamped; anything outside it is a real defect and raises.
_FRACTION_ROUNDOFF = 1e-9
#: The properties through which a posterior reports its own authority.  Sealed against redefinition by
#: a subclass; see :meth:`CellStatePosterior.__init_subclass__`.
_SEALED_AUTHORITY_PROPERTIES = ("native_accepted", "evidence_authority")
_WEIGHT_SUM_TOLERANCE = 1e-12
_QUANTILE_BISECTIONS = 200
_BRACKET_SIGMAS = 40.0


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _sha256_digest(value: object, *, what: str) -> str:
    text = _nonempty(value, what=what)
    if len(text) != 64 or any(character not in _SHA256_CHARS for character in text):
        raise ValueError(f"{what} must be a lowercase SHA-256 digest")
    return text


def _sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode()).hexdigest()


def _finite(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite")
    return number


def _nonnegative(value: object, *, what: str) -> float:
    number = _finite(value, what=what)
    if number < 0.0:
        raise ValueError(f"{what} must be nonnegative")
    return number


def _unit_interval(value: object, *, what: str, open_lower: bool, open_upper: bool) -> float:
    number = _finite(value, what=what)
    low_ok = number > 0.0 if open_lower else number >= 0.0
    high_ok = number < 1.0 if open_upper else number <= 1.0
    if not (low_ok and high_ok):
        bracket = f"{'(' if open_lower else '['}0, 1{')' if open_upper else ']'}"
        raise ValueError(f"{what} must lie in {bracket}")
    return number


def _fraction(value: object, *, what: str) -> float:
    """Return a fraction in ``[0, 1]``, clamping only the orthonormality round-off.

    The clamp is deliberately tiny and one-sided-in-magnitude: a fraction of 1.0000000000000004 is
    the eigenvector row summing to one in floating point, while 1.02 is a defect that must surface.
    """
    number = _finite(value, what=what)
    if -_FRACTION_ROUNDOFF <= number < 0.0:
        return 0.0
    if 1.0 < number <= 1.0 + _FRACTION_ROUNDOFF:
        return 1.0
    if not 0.0 <= number <= 1.0:
        raise ValueError(f"{what} must lie in [0, 1]; got {number!r}")
    return number


class UncertaintyCombination(StrEnum):
    """How one uncertainty source combines with the others in a budget.

    The distinction is not cosmetic: it decides whether collapsing the budget to one Gaussian is
    exact or lossy, and the two answers differ by a factor of two in the fourth moment for the pair
    the audit register already measured.

    ``CONVOLUTION``
        The source adds independently to the same quantity — thermal jitter on top of whatever else
        is happening.  Variances add, and for Gaussian sources the sum is again exactly Gaussian.
    ``MIXTURE``
        The source selects between sub-populations — cell-to-cell variability, a quenched structural
        draw, a metastable branch.  The result is a mixture, whose fourth moment is NOT that of the
        moment-matched Gaussian, which is the whole content of audit item A1.
    """

    CONVOLUTION = "convolution"
    MIXTURE = "mixture"


@dataclass(frozen=True, slots=True)
class UncertaintyComponent:
    """One mechanistically distinct uncertainty source, with how it was propagated.

    Args:
        source: Which of the eight declared sources this is.
        variance: Variance contributed, in the parameter's unit squared.  ``0.0`` is meaningful and
            common: it says the source EXISTS and is not propagated, which a missing entry would not.
        combination: Whether this source convolves with the others or mixes over them.
        propagation: How the number was obtained or why it is zero.  Required, because "0.0" and
            "0.0 because the explicit runtime carries one trajectory" are different claims.
        evidence_id: The artifact, source-graph node or report this component came from.
        weight: Mixture weight; required for ``MIXTURE`` and forbidden for ``CONVOLUTION``.
        mean_offset: Location of this component relative to the posterior's declared centre, in the
            parameter's unit.
    """

    source: UncertaintySource
    variance: float
    combination: UncertaintyCombination
    propagation: str
    evidence_id: str
    weight: float | None = None
    mean_offset: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.source, UncertaintySource):
            raise TypeError("source must be an UncertaintySource")
        if not isinstance(self.combination, UncertaintyCombination):
            raise TypeError("combination must be an UncertaintyCombination")
        object.__setattr__(self, "variance", _nonnegative(self.variance, what="variance"))
        object.__setattr__(self, "propagation", _nonempty(self.propagation, what="propagation"))
        object.__setattr__(self, "evidence_id", _nonempty(self.evidence_id, what="evidence_id"))
        object.__setattr__(self, "mean_offset", _finite(self.mean_offset, what="mean_offset"))
        if self.combination is UncertaintyCombination.MIXTURE:
            if self.weight is None:
                raise ValueError(
                    f"{self.source.value} mixes over sub-populations and must declare its weight"
                )
            object.__setattr__(
                self,
                "weight",
                _unit_interval(self.weight, what="weight", open_lower=True, open_upper=False),
            )
        elif self.weight is not None:
            raise ValueError(
                f"{self.source.value} combines by convolution and cannot carry a mixture weight"
            )

    @property
    def propagated(self) -> bool:
        """Whether this source actually contributes spread.  A property; no field to overwrite."""
        return self.variance > 0.0

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this component."""
        return {
            "source": self.source.value,
            "variance": self.variance,
            "combination": self.combination.value,
            "propagation": self.propagation,
            "evidence_id": self.evidence_id,
            "weight": self.weight,
            "mean_offset": self.mean_offset,
        }


class LossyMergeRefused(RuntimeError):
    """Raised when collapsing a decomposed budget into one Gaussian loses more than was allowed.

    The measured error is carried on the exception, so the refusal is a measurement rather than an
    opinion, and the message names the expansion it routes to.
    """

    def __init__(
        self,
        *,
        parameter: str,
        measured_error: float,
        fourth_moment_ratio: float,
        total_variance: float,
        tolerance: float,
    ) -> None:
        self.parameter = parameter
        self.measured_error = measured_error
        self.fourth_moment_ratio = fourth_moment_ratio
        self.total_variance = total_variance
        self.tolerance = tolerance
        super().__init__(
            f"refusing to merge the uncertainty budget of {parameter!r} into one Gaussian: the "
            f"moment-matched Gaussian reproduces the variance EXACTLY ({total_variance:.6g}) while "
            f"misreporting E[x^4] by a factor {fourth_moment_ratio:.6g} (relative error "
            f"{measured_error:.6g} > declared tolerance {tolerance:.6g}). Master-plan section 2.2 and "
            f"audit A1 reject the single-Gaussian collapse; route to a GAUSSIAN_MIXTURE "
            f"representation or partition by regime instead of widening one sigma."
        )


@dataclass(frozen=True, slots=True)
class GaussianCollapse:
    """A DECLARED-lossy single-Gaussian summary that keeps the decomposition it summarised.

    The point of returning the budget alongside the collapse is that the decomposition survives: a
    caller can always ask which source dominated, which is the information a single ``Sigma`` destroys.
    """

    parameter: str
    unit: str
    mean_offset: float
    variance: float
    fourth_moment_ratio: float
    error_report: ApproximationErrorReport
    budget: UncertaintyBudget

    def __post_init__(self) -> None:
        object.__setattr__(self, "parameter", _nonempty(self.parameter, what="parameter"))
        object.__setattr__(self, "unit", _nonempty(self.unit, what="unit"))
        object.__setattr__(self, "mean_offset", _finite(self.mean_offset, what="mean_offset"))
        object.__setattr__(self, "variance", _nonnegative(self.variance, what="variance"))
        object.__setattr__(
            self,
            "fourth_moment_ratio",
            _nonnegative(self.fourth_moment_ratio, what="fourth_moment_ratio"),
        )
        if not isinstance(self.error_report, ApproximationErrorReport):
            raise TypeError("error_report must be an ApproximationErrorReport")
        if not isinstance(self.budget, UncertaintyBudget):
            raise TypeError("budget must be an UncertaintyBudget")

    @property
    def lossy(self) -> bool:
        """Always ``True``.  A property with no field behind it: a collapse cannot be marked exact."""
        return True

    @property
    def decomposition(self) -> Mapping[UncertaintySource, float]:
        """The per-source variances this collapse summarised, still available."""
        return self.budget.variance_by_source()


@dataclass(frozen=True, slots=True)
class UncertaintyBudget:
    """The uncertainty of one quantity, factored by mechanistic source.  Never one number.

    There is deliberately NO constructor taking a total variance or a covariance.  The only way to
    build a budget is to name the sources, and :class:`~aleph.virtual_cell.contracts.UncertaintySource`
    has no ``TOTAL`` member, so a lumped width cannot even be spelled.  Totals are computed on demand
    by :meth:`total_variance`, which leaves :meth:`variance_by_source` available afterwards.

    Args:
        parameter: What the uncertainty is about.
        unit: The quantity's unit; variances are in this unit squared.
        components: One entry per contributing source, each source appearing at most once.
    """

    parameter: str
    unit: str
    components: tuple[UncertaintyComponent, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "parameter", _nonempty(self.parameter, what="parameter"))
        object.__setattr__(self, "unit", _nonempty(self.unit, what="unit"))
        components = tuple(self.components)
        if not components:
            raise ValueError(
                f"the uncertainty budget of {self.parameter!r} names no source; an unfactored width "
                "is exactly what master-plan section 2.2 and audit A1 reject"
            )
        if any(not isinstance(item, UncertaintyComponent) for item in components):
            raise TypeError("components must contain UncertaintyComponent values")
        sources = tuple(item.source for item in components)
        if len(set(sources)) != len(sources):
            raise ValueError(
                "each uncertainty source appears at most once; two entries for one source is a "
                "merge wearing a decomposition"
            )
        object.__setattr__(self, "components", components)
        weights = [item.weight for item in components if item.weight is not None]
        if weights and abs(sum(weights) - 1.0) > _WEIGHT_SUM_TOLERANCE:
            raise ValueError(f"mixture weights must sum to 1; got {sum(weights)!r}")

    # -- structure ------------------------------------------------------------------------------

    @property
    def declared_sources(self) -> tuple[UncertaintySource, ...]:
        """The sources this budget names, in declaration order."""
        return tuple(item.source for item in self.components)

    @property
    def unpropagated_sources(self) -> tuple[UncertaintySource, ...]:
        """Sources that are NAMED but carry no spread — the honest form of 'not modelled yet'."""
        return tuple(item.source for item in self.components if not item.propagated)

    @property
    def mixture_components(self) -> tuple[UncertaintyComponent, ...]:
        """Components that select between sub-populations."""
        return tuple(
            item for item in self.components if item.combination is UncertaintyCombination.MIXTURE
        )

    @property
    def convolution_components(self) -> tuple[UncertaintyComponent, ...]:
        """Components that add independently."""
        return tuple(
            item
            for item in self.components
            if item.combination is UncertaintyCombination.CONVOLUTION
        )

    def variance_by_source(self) -> Mapping[UncertaintySource, float]:
        """Return the per-source variances.  This is what a single covariance would have erased."""
        return MappingProxyType({item.source: item.variance for item in self.components})

    # -- moments, computed on demand -------------------------------------------------------------

    def _folded_mixture(self) -> tuple[tuple[float, float, float], ...]:
        """Return ``(weight, mean_offset, variance)`` after folding the convolution part in.

        A mixture convolved with an independent Gaussian is again a mixture: every branch keeps its
        weight, shifts by the convolution mean and inflates by the convolution variance.  That is why
        the credible interval below needs no approximation.
        """
        conv_mean = sum(item.mean_offset for item in self.convolution_components)
        conv_var = sum(item.variance for item in self.convolution_components)
        mixture = self.mixture_components
        if not mixture:
            return ((1.0, conv_mean, conv_var),)
        return tuple(
            (float(item.weight), item.mean_offset + conv_mean, item.variance + conv_var)
            for item in mixture
        )

    def mean_offset(self) -> float:
        """Mean of the deviation this budget describes, relative to the declared centre."""
        return sum(weight * mean for weight, mean, _ in self._folded_mixture())

    def total_variance(self) -> float:
        """Total variance, computed ON DEMAND; the decomposition stays available afterwards.

        Returns:
            The variance in the parameter's unit squared.
        """
        folded = self._folded_mixture()
        mean = sum(weight * component_mean for weight, component_mean, _ in folded)
        return (
            sum(
                weight * (variance + component_mean**2)
                for weight, component_mean, variance in folded
            )
            - mean**2
        )

    def fourth_central_moment(self) -> float:
        """``E[(x - E x)^4]`` of the deviation, in the parameter's unit to the fourth power."""
        folded = self._folded_mixture()
        mean = sum(weight * component_mean for weight, component_mean, _ in folded)
        return sum(
            weight * (delta**4 + 6.0 * delta**2 * variance + 3.0 * variance**2)
            for weight, component_mean, variance in folded
            for delta in (component_mean - mean,)
        )

    def fourth_moment_ratio(self) -> float:
        """``E[x^4]`` divided by the ``3 var^2`` a moment-matched Gaussian would report.

        Returns:
            ``1.0`` exactly when the budget is Gaussian (any number of convolved sources).  The
            registered A1 counterexample returns 1.9607881580237234.  A degenerate zero-variance
            budget returns ``1.0``: a point has no shape to misreport.
        """
        variance = self.total_variance()
        if variance <= 0.0:
            return 1.0
        return self.fourth_central_moment() / (3.0 * variance**2)

    def excess_kurtosis(self) -> float:
        """``E[x^4] / var^2 - 3``; zero for a Gaussian, positive for the A1 mixture."""
        variance = self.total_variance()
        if variance <= 0.0:
            return 0.0
        return self.fourth_central_moment() / variance**2 - 3.0

    # -- the explicitly lossy operation ----------------------------------------------------------

    def collapse_to_gaussian(self, *, max_relative_fourth_moment_error: float) -> GaussianCollapse:
        """Summarise the whole budget as one Gaussian, or REFUSE because it loses too much.

        The collapse is variance-preserving by construction, which is exactly why the variance cannot
        detect it.  The measured error is therefore the relative fourth-moment error,
        ``|E[x^4] / (3 var^2) - 1|``.

        Args:
            max_relative_fourth_moment_error: Declared tolerance in ``(0, 1)``.  Required, with no
                default: a tolerance chosen after seeing the error is not a tolerance.

        Returns:
            A :class:`GaussianCollapse` flagged lossy, carrying a measured
            :class:`~aleph.virtual_cell.contracts.ApproximationErrorReport` and the original budget.

        Raises:
            LossyMergeRefused: If the measured relative fourth-moment error exceeds the tolerance.
        """
        tolerance = _unit_interval(
            max_relative_fourth_moment_error,
            what="max_relative_fourth_moment_error",
            open_lower=True,
            open_upper=True,
        )
        variance = self.total_variance()
        ratio = self.fourth_moment_ratio()
        measured = abs(ratio - 1.0)
        if measured > tolerance:
            raise LossyMergeRefused(
                parameter=self.parameter,
                measured_error=measured,
                fourth_moment_ratio=ratio,
                total_variance=variance,
                tolerance=tolerance,
            )
        mean = self.mean_offset()
        source_hash = _sha256_json(self.to_dict())
        reduced_hash = _sha256_json({"mean_offset": mean, "variance": variance, "family": "normal"})
        report = ApproximationErrorReport(
            report_sha256=_sha256_json(
                {
                    "source": source_hash,
                    "reduced": reduced_hash,
                    "metric_id": "relative-fourth-moment-error",
                    "requested_tolerance": tolerance,
                    "measured_error": measured,
                }
            ),
            source_blob_sha256=source_hash,
            reduced_blob_sha256=reduced_hash,
            metric_id="relative-fourth-moment-error",
            requested_tolerance=tolerance,
            measured_error=measured,
        )
        return GaussianCollapse(
            parameter=self.parameter,
            unit=self.unit,
            mean_offset=mean,
            variance=variance,
            fourth_moment_ratio=ratio,
            error_report=report,
            budget=self,
        )

    # -- exact interval --------------------------------------------------------------------------

    def cdf(self, value: float) -> float:
        """Exact CDF of the deviation distribution at ``value``.

        Zero-variance branches are atoms and contribute their whole weight at and above their
        location; that is what makes a point posterior's interval degenerate rather than undefined.
        """
        point = _finite(value, what="value")
        total = 0.0
        for weight, mean, variance in self._folded_mixture():
            if variance <= 0.0:
                total += weight if point >= mean else 0.0
            else:
                total += weight * 0.5 * (1.0 + math.erf((point - mean) / math.sqrt(2.0 * variance)))
        return total

    def quantile(self, probability: float) -> float:
        """Exact quantile of the deviation distribution, by bisection on the monotone CDF.

        Args:
            probability: A probability in ``(0, 1)``.

        Returns:
            The deviation ``q`` with ``cdf(q) >= probability``, to the float resolution of the answer.
        """
        target = _unit_interval(probability, what="probability", open_lower=True, open_upper=True)
        folded = self._folded_mixture()
        variance = self.total_variance()
        if variance <= 0.0:
            locations = {mean for _, mean, _ in folded}
            if len(locations) == 1:
                return next(iter(locations))
        centre = self.mean_offset()
        spread = math.sqrt(variance) if variance > 0.0 else 0.0
        span = max(spread * _BRACKET_SIGMAS, 1.0)
        low, high = centre - span, centre + span
        while self.cdf(low) >= target:
            span *= 2.0
            low = centre - span
        while self.cdf(high) < target:
            span *= 2.0
            high = centre + span
        for _ in range(_QUANTILE_BISECTIONS):
            middle = 0.5 * (low + high)
            if middle in (low, high):
                break
            if self.cdf(middle) < target:
                low = middle
            else:
                high = middle
        return high

    def central_interval(self, *, level: float) -> tuple[float, float]:
        """Central credible interval of the deviation at ``level``, exact for any mixture.

        Args:
            level: Central probability in ``(0, 1)``.

        Returns:
            ``(low, high)`` deviations in the parameter's unit.
        """
        probability = _unit_interval(level, what="level", open_lower=True, open_upper=True)
        tail = 0.5 * (1.0 - probability)
        return (self.quantile(tail), self.quantile(1.0 - tail))

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able record of this budget."""
        return {
            "parameter": self.parameter,
            "unit": self.unit,
            "components": [item.to_dict() for item in self.components],
        }


@dataclass(frozen=True, slots=True)
class IdentifiabilityDeclaration:
    """How much of a parameter the available observations actually measured.

    Args:
        parameter: The parameter.
        identifiable_fraction: Fraction of the parameter's unit vector inside the measured subspace,
            in ``[0, 1]``.
        eigenvalue_floor: The declared information floor the fraction was computed at.  It belongs on
            the declaration because the fraction is meaningless without it.
        measured_by: What measured it — a report label, a design id, an artifact.
    """

    parameter: str
    identifiable_fraction: float
    eigenvalue_floor: float
    measured_by: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "parameter", _nonempty(self.parameter, what="parameter"))
        object.__setattr__(self, "measured_by", _nonempty(self.measured_by, what="measured_by"))
        object.__setattr__(
            self,
            "identifiable_fraction",
            _fraction(self.identifiable_fraction, what="identifiable_fraction"),
        )
        floor = _finite(self.eigenvalue_floor, what="eigenvalue_floor")
        if floor <= 0.0:
            raise ValueError("eigenvalue_floor must be positive")
        object.__setattr__(self, "eigenvalue_floor", floor)

    @classmethod
    def from_report(
        cls,
        report: IdentifiabilityReport,
        *,
        eigenvalue_floor: float,
        parameters: Sequence[str] = PARAMETER_NAMES,
    ) -> tuple[IdentifiabilityDeclaration, ...]:
        """Read declarations straight out of an S11 image-inverse identifiability report.

        This is the composition point with :mod:`aleph.virtual_cell.sandbox_inverse`: the fractions
        are not recomputed here, so the posterior refuses on exactly the numbers that module measured.

        Args:
            report: The Fisher eigen-decomposition.
            eigenvalue_floor: The caller's declared information floor.
            parameters: Which parameters to read; defaults to the module's own parameter names.

        Returns:
            One declaration per parameter, in the order given.
        """
        if not isinstance(report, IdentifiabilityReport):
            raise TypeError("report must be an IdentifiabilityReport")
        return tuple(
            cls(
                parameter=name,
                identifiable_fraction=report.identifiable_fraction(
                    name, eigenvalue_floor=eigenvalue_floor
                ),
                eigenvalue_floor=eigenvalue_floor,
                measured_by=f"sandbox_inverse.IdentifiabilityReport[{report.label}]",
            )
            for name in parameters
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this declaration."""
        return {
            "parameter": self.parameter,
            "identifiable_fraction": self.identifiable_fraction,
            "eigenvalue_floor": self.eigenvalue_floor,
            "measured_by": self.measured_by,
        }


@dataclass(frozen=True, slots=True)
class PosteriorEvidenceContract:
    """The caller's pre-declared posterior thresholds.  **No field has a default, by design.**

    Same reasoning as :class:`~aleph.virtual_cell.regime.RegimeEvidenceContract`: every number the
    posterior compares against is declared before the run and copied into the record, so a verdict can
    be re-judged and a quietly relaxed threshold cannot masquerade as the original.

    Args:
        min_identifiable_fraction: How much of a parameter must lie in the measured subspace before a
            POINT estimate is allowed at all.  In ``(0, 1]``.
        eigenvalue_floor: Information floor the declarations must have been measured at.  A
            declaration measured at a different floor is refused rather than rescaled.
        max_relative_fourth_moment_error: Tolerance for the explicitly lossy Gaussian collapse.
        credible_level: Central probability for :meth:`CellStatePosterior.credible_interval`.
        max_representation_error: Largest measured approximation error an approximate representation
            may carry and still be admitted.  Fails closed: a reduced representation whose measured
            error exceeds this is refused at construction, not at reporting time.
    """

    min_identifiable_fraction: float
    eigenvalue_floor: float
    max_relative_fourth_moment_error: float
    credible_level: float
    max_representation_error: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "min_identifiable_fraction",
            _unit_interval(
                self.min_identifiable_fraction,
                what="min_identifiable_fraction",
                open_lower=True,
                open_upper=False,
            ),
        )
        floor = _finite(self.eigenvalue_floor, what="eigenvalue_floor")
        if floor <= 0.0:
            raise ValueError("eigenvalue_floor must be positive")
        object.__setattr__(self, "eigenvalue_floor", floor)
        object.__setattr__(
            self,
            "max_relative_fourth_moment_error",
            _unit_interval(
                self.max_relative_fourth_moment_error,
                what="max_relative_fourth_moment_error",
                open_lower=True,
                open_upper=True,
            ),
        )
        object.__setattr__(
            self,
            "credible_level",
            _unit_interval(
                self.credible_level, what="credible_level", open_lower=True, open_upper=True
            ),
        )
        object.__setattr__(
            self,
            "max_representation_error",
            _unit_interval(
                self.max_representation_error,
                what="max_representation_error",
                open_lower=False,
                open_upper=True,
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able contract, copied verbatim into every posterior record."""
        return {
            "min_identifiable_fraction": self.min_identifiable_fraction,
            "eigenvalue_floor": self.eigenvalue_floor,
            "max_relative_fourth_moment_error": self.max_relative_fourth_moment_error,
            "credible_level": self.credible_level,
            "max_representation_error": self.max_representation_error,
        }


@dataclass(frozen=True, slots=True)
class PosteriorGraph:
    """The component/connector graph half of the posterior, bound to the canonical grammar.

    Co-location is never a connection (CLAUDE.md), so the graph is the DECLARED grammar plus the
    instantiated connector edges — never a component list on its own.  A connector instance naming an
    edge the grammar does not declare is refused here rather than discovered later.
    """

    grammar: EngineGrammar
    connectors: tuple[ConnectorInstance, ...] = ()
    population: PopulationWiring | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.grammar, EngineGrammar):
            raise TypeError("grammar must be an EngineGrammar")
        connectors = tuple(self.connectors)
        if any(not isinstance(item, ConnectorInstance) for item in connectors):
            raise TypeError("connectors must contain ConnectorInstance values")
        instance_ids = tuple(item.instance_id for item in connectors)
        if len(set(instance_ids)) != len(instance_ids):
            raise ValueError("connector instance ids must be unique")
        undeclared = sorted(
            {
                item.connector_name
                for item in connectors
                if item.connector_name not in self.grammar.connector_names
            }
        )
        if undeclared:
            raise ValueError(
                "connector instances name edges the grammar does not declare: "
                + ", ".join(undeclared)
            )
        object.__setattr__(self, "connectors", connectors)
        if self.population is not None:
            if not isinstance(self.population, PopulationWiring):
                raise TypeError("population must be a PopulationWiring")
            known = {item.instance_id for item in self.population.connectors}
            orphaned = sorted(set(instance_ids) - known)
            if orphaned:
                raise ValueError(
                    "connector instances are not part of the declared population wiring: "
                    + ", ".join(orphaned)
                )

    @property
    def component_names(self) -> frozenset[str]:
        """Every component the grammar declares."""
        return self.grammar.component_names

    @property
    def connector_instance_ids(self) -> tuple[str, ...]:
        """Instantiated connector edges, in declaration order."""
        return tuple(item.instance_id for item in self.connectors)

    @property
    def cell_cell_junction_count(self) -> int:
        """Number of explicit cell-cell junctions.

        Reported as its own number because Lane E measured that the canonical census contains none:
        a tissue-scale posterior is blocked on a new connector LAW, not on more parameters, and a
        posterior that silently reported zero as 'fine' would hide that.
        """
        if self.population is None:
            return 0
        return len(self.population.junctions)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able summary of the graph."""
        return {
            "grammar_source": self.grammar.source,
            "n_components": len(self.grammar.components),
            "n_declared_connectors": len(self.grammar.connectors),
            "connector_instance_ids": list(self.connector_instance_ids),
            "cell_cell_junction_count": self.cell_cell_junction_count,
        }


@dataclass(frozen=True, slots=True)
class MolecularParameterPosterior:
    """The molecular-parameter posterior: a location, a representation, and a budget per parameter.

    The representation is checked AGAINST the content, in both directions.  A ``GAUSSIAN_MIXTURE``
    with no mixture anywhere, or a ``LOCAL_GAUSSIAN`` whose budget carries several branches, is
    rejected: a label that does not describe the content is the over-claim the audit register exists
    to prevent, and it is exactly as wrong when it under-claims.
    """

    parameters: tuple[str, ...]
    units: Mapping[str, str]
    location: Mapping[str, float]
    budgets: Mapping[str, UncertaintyBudget]
    representation: RepresentationDescriptor
    identifiability: Mapping[str, IdentifiabilityDeclaration] = field(default_factory=dict)

    def __post_init__(self) -> None:
        names = tuple(_nonempty(name, what="parameter name") for name in self.parameters)
        if not names:
            raise ValueError("a parameter posterior must cover at least one parameter")
        if len(set(names)) != len(names):
            raise ValueError("parameter names must be unique")
        object.__setattr__(self, "parameters", names)
        expected = set(names)
        for attribute in ("units", "location", "budgets"):
            mapping = getattr(self, attribute)
            if not isinstance(mapping, Mapping):
                raise TypeError(f"{attribute} must be a mapping")
            if set(mapping) != expected:
                raise ValueError(f"{attribute} must have exactly one entry per parameter")
        object.__setattr__(
            self,
            "units",
            MappingProxyType(
                {name: _nonempty(self.units[name], what=f"unit of {name!r}") for name in names}
            ),
        )
        object.__setattr__(
            self,
            "location",
            MappingProxyType(
                {name: _finite(self.location[name], what=f"location of {name!r}") for name in names}
            ),
        )
        budgets: dict[str, UncertaintyBudget] = {}
        for name in names:
            budget = self.budgets[name]
            if not isinstance(budget, UncertaintyBudget):
                raise TypeError("budgets must contain UncertaintyBudget values")
            if budget.parameter != name:
                raise ValueError(f"budget for {name!r} declares parameter {budget.parameter!r}")
            if budget.unit != self.units[name]:
                raise ValueError(
                    f"budget for {name!r} is in {budget.unit!r} but the parameter is in "
                    f"{self.units[name]!r}"
                )
            budgets[name] = budget
        object.__setattr__(self, "budgets", MappingProxyType(budgets))

        if not isinstance(self.representation, RepresentationDescriptor):
            raise TypeError("representation must be a RepresentationDescriptor")
        if self.representation.kind not in POSTERIOR_REPRESENTATION_KINDS:
            raise ValueError(
                f"{self.representation.kind.value} is not a way of REPRESENTING a posterior; "
                f"allowed: {sorted(kind.value for kind in POSTERIOR_REPRESENTATION_KINDS)}"
            )
        self._check_representation_matches_content()

        if not isinstance(self.identifiability, Mapping):
            raise TypeError("identifiability must be a mapping")
        declarations: dict[str, IdentifiabilityDeclaration] = {}
        for key, declaration in self.identifiability.items():
            name = _nonempty(key, what="identifiability key")
            if name not in expected:
                raise ValueError(f"identifiability names unknown parameter {name!r}")
            if not isinstance(declaration, IdentifiabilityDeclaration):
                raise TypeError("identifiability values must be IdentifiabilityDeclaration")
            if declaration.parameter != name:
                raise ValueError(f"identifiability[{name!r}] declares {declaration.parameter!r}")
            declarations[name] = declaration
        object.__setattr__(self, "identifiability", MappingProxyType(declarations))

    def _check_representation_matches_content(self) -> None:
        kind = self.representation.kind
        branch_counts = {
            name: len(budget.mixture_components) for name, budget in self.budgets.items()
        }
        if kind is POINT_REPRESENTATION:
            spread = sorted(
                name for name in self.parameters if self.budgets[name].total_variance() > 0.0
            )
            if spread:
                raise ValueError(
                    "a POINT (explicit-state) posterior carries no spread, but these parameters do: "
                    + ", ".join(spread)
                    + "; declare EMPIRICAL_ENSEMBLE / LOCAL_GAUSSIAN / GAUSSIAN_MIXTURE instead"
                )
        if kind is RepresentationKind.LOCAL_GAUSSIAN and max(branch_counts.values(), default=0) > 1:
            raise ValueError(
                "a LOCAL_GAUSSIAN cannot represent a multi-branch mixture; the branches would be "
                "merged into one sigma, which audit A1 rejects"
            )
        if (
            kind is RepresentationKind.GAUSSIAN_MIXTURE
            and max(branch_counts.values(), default=0) < 2
        ):
            raise ValueError(
                "a GAUSSIAN_MIXTURE posterior must actually carry a mixture; declaring the richer "
                "representation without the content over-states what was propagated"
            )

    # -- queries ---------------------------------------------------------------------------------

    def budget(self, parameter: str) -> UncertaintyBudget:
        """Return one parameter's decomposed budget.

        Raises:
            KeyError: If the parameter is not covered.
        """
        return self.budgets[parameter]

    def total_variance(self, parameter: str) -> float:
        """Total variance of one parameter, computed on demand from its decomposition."""
        return self.budgets[parameter].total_variance()

    def variance_by_source(self, parameter: str) -> Mapping[UncertaintySource, float]:
        """Per-source variances of one parameter; the decomposition, not a summary of it."""
        return self.budgets[parameter].variance_by_source()

    def unpropagated_sources(self, parameter: str) -> tuple[UncertaintySource, ...]:
        """Sources this parameter NAMES but does not propagate."""
        return self.budgets[parameter].unpropagated_sources

    def credible_interval(self, parameter: str, *, level: float) -> tuple[float, float]:
        """Central credible interval, exact for any declared budget.

        Always available.  For an unidentifiable parameter this is the correct answer and will simply
        be wide; for a point posterior it is degenerate, which is the visible cost of that
        representation rather than an error.

        Args:
            parameter: A covered parameter.
            level: Central probability in ``(0, 1)``.

        Returns:
            ``(low, high)`` in the parameter's own unit.
        """
        centre = self.location[parameter]
        low, high = self.budgets[parameter].central_interval(level=level)
        return (centre + low, centre + high)

    def point_estimate(self, parameter: str, *, min_identifiable_fraction: float) -> float:
        """Return a point estimate, or REFUSE because the observations did not measure it.

        Args:
            parameter: A covered parameter.
            min_identifiable_fraction: How much of the parameter's direction must lie in the measured
                subspace.  Required with no default.

        Returns:
            The declared location.

        Raises:
            KeyError: If the parameter is not covered.
            UnidentifiableDirectionError: If the measured identifiable fraction is below the
                requirement, OR if identifiability was never measured for this parameter.  The second
                case fails CLOSED: "nobody looked" is not "it is identifiable".
        """
        required = _unit_interval(
            min_identifiable_fraction,
            what="min_identifiable_fraction",
            open_lower=True,
            open_upper=False,
        )
        if parameter not in self.budgets:
            raise KeyError(parameter)
        declaration = self.identifiability.get(parameter)
        if declaration is None:
            raise UnidentifiableDirectionError(
                f"{parameter!r} has no identifiability declaration, so no point estimate is "
                f"available. This fails closed on purpose: an unmeasured direction and an "
                f"identifiable one are indistinguishable from the posterior alone. Use "
                f"credible_interval, or measure identifiability (sandbox_inverse.identifiability_"
                f"report) and attach an IdentifiabilityDeclaration."
            )
        if declaration.identifiable_fraction < required:
            raise UnidentifiableDirectionError(
                f"{parameter!r} is not identifiable from the available observations: identifiable "
                f"fraction {declaration.identifiable_fraction:.4g} < required {required:.4g} at "
                f"eigenvalue floor {declaration.eigenvalue_floor:.4g} "
                f"(measured by {declaration.measured_by}). Master-plan A5 rejects a point estimate "
                f"along an unmeasured direction; use credible_interval, or add an observation "
                f"modality (see sandbox_inverse.measure_modality_contrast)."
            )
        return self.location[parameter]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of the parameter posterior."""
        return {
            "parameters": list(self.parameters),
            "units": dict(sorted(self.units.items())),
            "location": dict(sorted(self.location.items())),
            "budgets": {name: self.budgets[name].to_dict() for name in self.parameters},
            "representation_kind": self.representation.kind.value,
            "representation_axes": list(self.representation.axes),
            "approximation_error_bound": self.representation.approximation_error_bound,
            "identifiability": {
                name: declaration.to_dict()
                for name, declaration in sorted(self.identifiability.items())
            },
        }


@dataclass(frozen=True, slots=True)
class InterventionResponse:
    """One intervention's response distribution, bound to the card that pre-registered it.

    ``card_hash`` is required.  A response with no pre-registration is the self-confirming sandbox
    result audit item A7 rejects, and making the binding a required field is cheaper than arguing
    about it afterwards.
    """

    intervention_id: str
    observable: str
    unit: str
    location: float
    budget: UncertaintyBudget
    card_hash: str

    def __post_init__(self) -> None:
        for name in ("intervention_id", "observable", "unit"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        object.__setattr__(self, "location", _finite(self.location, what="location"))
        if not isinstance(self.budget, UncertaintyBudget):
            raise TypeError("budget must be an UncertaintyBudget")
        if self.budget.unit != self.unit:
            raise ValueError(
                f"response budget is in {self.budget.unit!r} but the observable is in {self.unit!r}"
            )
        object.__setattr__(self, "card_hash", _sha256_digest(self.card_hash, what="card_hash"))

    def credible_interval(self, *, level: float) -> tuple[float, float]:
        """Central credible interval of the response, exact for any declared budget."""
        low, high = self.budget.central_interval(level=level)
        return (self.location + low, self.location + high)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this response."""
        return {
            "intervention_id": self.intervention_id,
            "observable": self.observable,
            "unit": self.unit,
            "location": self.location,
            "budget": self.budget.to_dict(),
            "card_hash": self.card_hash,
        }


@dataclass(frozen=True, slots=True)
class CellStatePosterior:
    """The programme's central data unit — master plan section 0, audit item A13.

    Every field binds a type an existing lane owns.  The invariants enforced at construction are the
    ones the plan and the audit already argued for; they are listed in the module docstring and each
    one raises rather than warns.
    """

    posterior_id: str
    manifest: CellStateManifest
    graph: PosteriorGraph
    priors: tuple[PriorDistribution, ...]
    parameters: MolecularParameterPosterior
    contract: PosteriorEvidenceContract
    evidence: EvidenceGraph
    envelope: ArtifactEnvelope
    accepted_step_evidence: AdmissibilityVerdict | None = None
    regime: RegimeReport | None = None
    interventions: tuple[InterventionResponse, ...] = ()
    observation_operators: tuple[ObservationOperator, ...] = ()
    model_inadequacy: tuple[MissingLawProposal, ...] = ()

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Refuse a subclass that redefines an authority property.

        "A property over ``__slots__`` cannot be spelled" was FALSE, and an external review
        (2026-07-29, finding 4) spelled it: ``class ForgedPosterior(CellStatePosterior)`` overriding
        ``native_accepted`` to return ``True`` passed ``isinstance``, constructed from the base
        dataclass fields, and serialised ``to_dict()["native_accepted"] == True`` — an actual
        native self-promotion, in the record.  ``dataclasses.replace`` was never the only route to a
        forged field; inheritance was the other one, and it is closed here.

        Raises:
            TypeError: If the subclass defines ``native_accepted`` or ``evidence_authority``.
        """
        # Two-argument super() on purpose: ``@dataclass(slots=True)`` REPLACES the class object, so
        # the zero-argument form's ``__class__`` cell points at the discarded original and raises.
        super(CellStatePosterior, cls).__init_subclass__(**kwargs)
        forged = sorted(name for name in _SEALED_AUTHORITY_PROPERTIES if name in cls.__dict__)
        if forged:
            raise TypeError(
                f"{cls.__name__} redefines {', '.join(forged)}, which CellStatePosterior seals: the "
                "generic virtual-cell layer is an unverified observer and may not spell a native "
                "rung for itself by any route, inheritance included. Carry the native artifact "
                "through envelope.source_artifact_ids and let the runtime-attested adapter adjudicate."
            )

    def _assert_authority_is_not_forged(self) -> None:
        """Second belt: the resolved properties must still be this class's own.

        ``__init_subclass__`` cannot see an override installed on the class AFTER it was created
        (``ForgedPosterior.native_accepted = property(...)``), so the identity is re-checked at
        construction and again at serialisation, where the forged value would otherwise be written
        into the record.

        Raises:
            TypeError: If either authority property has been replaced on this instance's class.
        """
        for name in _SEALED_AUTHORITY_PROPERTIES:
            if getattr(type(self), name, None) is not getattr(CellStatePosterior, name):
                raise TypeError(
                    f"{type(self).__name__}.{name} is not CellStatePosterior's own property; the "
                    "authority a posterior reports about itself may not be replaced (external "
                    "review 2026-07-29, finding 4)"
                )

    def __post_init__(self) -> None:
        self._assert_authority_is_not_forged()
        object.__setattr__(self, "posterior_id", _nonempty(self.posterior_id, what="posterior_id"))
        for name, kind in (
            ("manifest", CellStateManifest),
            ("graph", PosteriorGraph),
            ("parameters", MolecularParameterPosterior),
            ("contract", PosteriorEvidenceContract),
            ("evidence", EvidenceGraph),
            ("envelope", ArtifactEnvelope),
        ):
            if not isinstance(getattr(self, name), kind):
                raise TypeError(f"{name} must be a {kind.__name__}")

        priors = tuple(self.priors)
        if any(not isinstance(prior, PriorDistribution) for prior in priors):
            raise TypeError("priors must contain PriorDistribution values")
        prior_ids = tuple(prior.prior_id for prior in priors)
        if len(set(prior_ids)) != len(prior_ids):
            raise ValueError("priors must not repeat a prior_id")
        object.__setattr__(self, "priors", priors)

        if self.envelope.cell_state_manifest_hash != self.manifest.manifest_hash:
            raise ValueError(
                "the envelope is bound to a different manifest than the posterior carries; a "
                "provenance wrapper that does not name its own subject is not provenance"
            )
        if self.envelope.representation.kind is not self.parameters.representation.kind:
            raise ValueError(
                "the envelope and the parameter posterior declare different representations"
            )
        bound = self.parameters.representation.approximation_error_bound
        if bound is not None and bound > self.contract.max_representation_error:
            raise ValueError(
                f"measured representation error {bound:.6g} exceeds the declared "
                f"max_representation_error {self.contract.max_representation_error:.6g}; a reduced "
                "representation is used only when its measured error passes the declared tolerance"
            )
        for declaration in self.parameters.identifiability.values():
            if declaration.eigenvalue_floor != self.contract.eigenvalue_floor:
                raise ValueError(
                    f"identifiability of {declaration.parameter!r} was measured at eigenvalue floor "
                    f"{declaration.eigenvalue_floor:.6g} but the contract declares "
                    f"{self.contract.eigenvalue_floor:.6g}; a fraction is meaningless at another "
                    "floor and is not rescaled here"
                )

        undeclared = sorted(set(self.manifest.connectors) - set(self.graph.grammar.connector_names))
        if undeclared:
            raise ValueError(
                "the manifest names connectors the graph's grammar does not declare: "
                + ", ".join(undeclared)
            )

        if self.accepted_step_evidence is not None and not isinstance(
            self.accepted_step_evidence, AdmissibilityVerdict
        ):
            raise TypeError("accepted_step_evidence must be an AdmissibilityVerdict")
        if self.regime is not None:
            if not isinstance(self.regime, RegimeReport):
                raise TypeError("regime must be a RegimeReport")
            if self.regime.evidence_authority != EVIDENCE_AUTHORITY:
                raise ValueError(
                    "a regime report carried here must retain the unverified-observer authority"
                )
            verdict = self.accepted_step_evidence
            if verdict is None:
                raise ValueError(
                    "a regime label may not be carried without accepted-step evidence: a regime read "
                    "off steps that were never adjudicated is a statement about the solver. Attach "
                    "an admissibility.AdmissibilityVerdict, or drop the regime."
                )
            if not (verdict.admissible and verdict.transaction_sound):
                raise ValueError(
                    "the accepted-step evidence says the step was not admissible or not "
                    "transaction-sound, so a regime measured on it is not a regime; rejected by: "
                    + ", ".join(clause.value for clause in verdict.rejected_by)
                )

        interventions = tuple(self.interventions)
        if any(not isinstance(item, InterventionResponse) for item in interventions):
            raise TypeError("interventions must contain InterventionResponse values")
        intervention_ids = tuple(item.intervention_id for item in interventions)
        if len(set(intervention_ids)) != len(intervention_ids):
            raise ValueError("interventions must not repeat an intervention_id")
        object.__setattr__(self, "interventions", interventions)

        # The protocol check this replaces was structural only, so an operator could satisfy it and
        # still lie about what it is.  External review 2026-07-29 finding 3: a subclass of
        # SummaryProjectionOperator declaring modality=TRACTION, under the reserved operator id
        # "unimplemented::traction", registered here and returned 5.0 for a BLOCKED modality.  The
        # refusal belongs at the registry boundary, so this delegates to observation_operator's
        # validator instead of re-deciding it locally.
        operators = tuple(
            validate_operator_registration(item) for item in self.observation_operators
        )
        operator_ids = tuple(operator.operator_id for operator in operators)
        if len(set(operator_ids)) != len(operator_ids):
            raise ValueError("observation operators must not repeat an operator_id")
        object.__setattr__(self, "observation_operators", operators)

        proposals = tuple(self.model_inadequacy)
        if any(not isinstance(item, MissingLawProposal) for item in proposals):
            raise TypeError("model_inadequacy must contain MissingLawProposal values")
        if any(item.decided for item in proposals):
            raise ValueError("a model-inadequacy PROPOSAL cannot arrive already decided")
        object.__setattr__(self, "model_inadequacy", proposals)

    # -- authority: properties with no field behind them ------------------------------------------

    @property
    def native_accepted(self) -> bool:
        """Always ``False``.

        A property over ``__slots__``, so there is no field for ``dataclasses.replace`` to set and no
        keyword through which a caller could spell a native rung.  Same pattern as
        :class:`~aleph.virtual_cell.surrogate.SurrogatePrediction`.
        """
        return False

    @property
    def evidence_authority(self) -> str:
        """Always ``unverified-observer``.  The generic virtual-cell layer observes; it does not rule."""
        return EVIDENCE_AUTHORITY

    @property
    def representation_kind(self) -> RepresentationKind:
        """The declared representation of the parameter posterior."""
        return self.parameters.representation.kind

    @property
    def is_point_representation(self) -> bool:
        """Whether this posterior is a point — the representation the engine gives today."""
        return self.representation_kind is POINT_REPRESENTATION

    @property
    def pi_gap_prior_ids(self) -> tuple[str, ...]:
        """Priors that are open PI gaps: named, unsourced, and carrying no number."""
        return tuple(prior.prior_id for prior in self.priors if prior.status is PriorStatus.PI_GAP)

    @property
    def missing_law_ids(self) -> tuple[str, ...]:
        """Candidate law ids across every carried proposal, in order."""
        return tuple(law_id for proposal in self.model_inadequacy for law_id in proposal.law_ids)

    # -- queries ----------------------------------------------------------------------------------

    def credible_interval(
        self, parameter: str, *, level: float | None = None
    ) -> tuple[float, float]:
        """Central credible interval for a parameter at the contract's declared level.

        Args:
            parameter: A covered parameter.
            level: Override the contract's ``credible_level``.  Rarely wanted; the contract exists so
                the level is declared once, before the run.

        Returns:
            ``(low, high)`` in the parameter's own unit.
        """
        return self.parameters.credible_interval(
            parameter, level=self.contract.credible_level if level is None else level
        )

    def point_estimate(self, parameter: str) -> float:
        """Point estimate at the contract's declared identifiability requirement, or refuse.

        Raises:
            UnidentifiableDirectionError: As :meth:`MolecularParameterPosterior.point_estimate`.
        """
        return self.parameters.point_estimate(
            parameter, min_identifiable_fraction=self.contract.min_identifiable_fraction
        )

    def unidentifiable_parameters(self) -> tuple[str, ...]:
        """Every parameter for which a point estimate is refused, in declaration order."""
        refused: list[str] = []
        for name in self.parameters.parameters:
            try:
                self.point_estimate(name)
            except UnidentifiableDirectionError:
                refused.append(name)
        return tuple(refused)

    def collapse_uncertainty(self, parameter: str) -> GaussianCollapse:
        """Collapse one parameter's budget to a Gaussian at the contract's tolerance, or refuse.

        Raises:
            LossyMergeRefused: If the measured fourth-moment error exceeds the declared tolerance.
        """
        return self.parameters.budget(parameter).collapse_to_gaussian(
            max_relative_fourth_moment_error=self.contract.max_relative_fourth_moment_error
        )

    def operator(self, operator_id: str) -> ObservationOperator:
        """Return one carried observation operator.

        Raises:
            KeyError: If no operator with that id is carried.
        """
        for operator in self.observation_operators:
            if operator.operator_id == operator_id:
                return operator
        raise KeyError(
            f"no observation operator {operator_id!r}; carried: "
            f"{[item.operator_id for item in self.observation_operators]}"
        )

    def observe(self, operator_id: str, state: MechanisticState) -> ObservationResult:
        """Observe a mechanistic state through one of the carried operators.

        Raises:
            KeyError: If the operator is not carried.
            ~aleph.virtual_cell.observation_operator.ObservationRefused: If the operator refuses
                the state, including every :class:`UnimplementedOperator`.
        """
        return self.operator(operator_id).observe(state)

    def assert_observation_excluded(self, source_id: str) -> None:
        """Assert that a synthetic observation source is NOT in the biological-truth corpus.

        Args:
            source_id: A node of the carried evidence graph.

        Raises:
            KeyError: If the source is not registered.
            ValueError: If the graph considers it citable biological evidence, which for an
                observation this layer produced would be the laundering audit A6 forbids.
        """
        if self.evidence.exclusion(source_id) is None:
            raise ValueError(
                f"{source_id!r} is citable as biological evidence; an observation produced by this "
                "layer is synthetic and must stay excluded, transitively"
            )

    def _accepted_step_record(self) -> dict[str, Any] | None:
        """The accepted-step verdict's IDENTITY, not its verdict bit.

        ``admissible`` alone is one bit about a transaction that has an index, a clause vector, a
        rejection list and a rollback audit.  External review 2026-07-29, finding 1: two posteriors
        whose verdicts differed only in ``attempted_step_index`` (7 against 999) produced the same
        ``record_hash``, so two different accepted-step transactions were one content-addressed
        posterior.  What is recorded here is what distinguishes them.
        """
        verdict = self.accepted_step_evidence
        if verdict is None:
            return None
        return {
            "attempted_step_index": verdict.attempted_step_index,
            "admissible": verdict.admissible,
            "transaction_sound": verdict.transaction_sound,
            "authority": verdict.authority,
            "evidence_source": verdict.evidence_source.value,
            "rejected_by": [clause.value for clause in verdict.rejected_by],
            "predicate_results_sha256": verdict.predicate_results_sha256(),
            "rollback_bit_exact": None if verdict.rollback is None else verdict.rollback.bit_exact,
            "expansions": list(verdict.expansions),
        }

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of the whole posterior.

        Every constructor field reaches this record, because :attr:`record_hash` is the SHA-256 of it
        and a field the record omits is a difference the hash cannot see.  ``evidence`` and the
        accepted-step verdict were both such omissions until the external review of 2026-07-29 found
        them (finding 1).
        """
        self._assert_authority_is_not_forged()
        return {
            "posterior_id": self.posterior_id,
            "native_accepted": self.native_accepted,
            "evidence_authority": self.evidence_authority,
            "manifest_hash": self.manifest.manifest_hash,
            "graph": self.graph.to_dict(),
            "priors": [prior.to_dict() for prior in self.priors],
            "pi_gap_prior_ids": list(self.pi_gap_prior_ids),
            "parameters": self.parameters.to_dict(),
            "contract": self.contract.to_dict(),
            # The whole envelope, not just its id: the ancestry and the observation binding are the
            # part a re-used artifact_id would otherwise hide.
            "envelope_artifact_id": self.envelope.artifact_id,
            "evidence_source": self.envelope.evidence_source.value,
            "envelope": {
                "artifact_id": self.envelope.artifact_id,
                "cell_state_manifest_hash": self.envelope.cell_state_manifest_hash,
                "evidence_source": self.envelope.evidence_source.value,
                "representation_kind": self.envelope.representation.kind.value,
                "representation_axes": list(self.envelope.representation.axes),
                "source_artifact_ids": list(self.envelope.source_artifact_ids),
                "projection_hash": self.envelope.projection_hash,
                "observation_provenance_sha256": self.envelope.observation_provenance_sha256,
            },
            # The provenance DAG's own content address.  A posterior that re-used a source id for a
            # different graph was invisible here, which is exactly the 10-lane composition case.
            "evidence_graph_sha256": self.evidence.graph_hash,
            "evidence_source_ids": list(self.evidence.topological_order),
            "accepted_step_admissible": (
                None
                if self.accepted_step_evidence is None
                else self.accepted_step_evidence.admissible
            ),
            "accepted_step_evidence": self._accepted_step_record(),
            "regime": None if self.regime is None else self.regime.as_dict(),
            "interventions": [item.to_dict() for item in self.interventions],
            "observation_operators": [
                {
                    "operator_id": operator.operator_id,
                    "modality": operator.modality.value,
                    "output_units": operator.output_units,
                    "input_requirement": operator.input_requirement.to_dict(),
                }
                for operator in self.observation_operators
            ],
            "model_inadequacy": [proposal.to_dict() for proposal in self.model_inadequacy],
            "missing_law_ids": list(self.missing_law_ids),
        }

    @property
    def record_hash(self) -> str:
        """SHA-256 of the canonical record; two posteriors agree only if their whole record does.

        The claim in that sentence is the one the external review of 2026-07-29 refuted (finding 1):
        with ``evidence`` absent from :meth:`to_dict` and the accepted-step verdict collapsed to one
        boolean, two posteriors carrying different :class:`EvidenceGraph`\\ s and different attempted
        step indices hashed identically.  The test that would have caught it varied only
        ``posterior_id``.  Every constructor field is now covered, and the test is parametrized over
        each of them.
        """
        return _sha256_json(self.to_dict())

    # -- construction from an archetype build -----------------------------------------------------

    @classmethod
    def from_archetype_build(
        cls,
        build: ArchetypeBuild,
        *,
        posterior_id: str,
        graph: PosteriorGraph,
        parameters: MolecularParameterPosterior,
        contract: PosteriorEvidenceContract,
        evidence: EvidenceGraph,
        artifact_id: str,
        evidence_source: EvidenceSource,
        source_artifact_ids: Sequence[str] = (),
        **optional: Any,
    ) -> CellStatePosterior:
        """Build a posterior from a Ring A archetype build, or REFUSE a Ring B/C refusal.

        This is the composition point with :mod:`aleph.virtual_cell.archetypes`: the manifest and
        the population/geometry priors — including every ``PI_GAP`` — come straight from the build,
        so a posterior cannot quietly acquire priors the archetype registry does not have.

        Args:
            build: The archetype build.
            posterior_id: Identifier for the new posterior.
            graph: The component/connector graph.
            parameters: The molecular-parameter posterior.
            contract: The declared thresholds.
            evidence: The provenance graph.
            artifact_id: Identifier for the envelope this constructs.
            evidence_source: Provenance of the posterior.  ``NATIVE_*`` is refused by the envelope.
            source_artifact_ids: Ancestry; required by the envelope for derived representations.
            **optional: Passed through to the constructor (``regime``, ``interventions``,
                ``observation_operators``, ``model_inadequacy``, ``accepted_step_evidence``).

        Returns:
            The assembled :class:`CellStatePosterior`.

        Raises:
            TypeError: If ``build`` is not an :class:`~aleph.virtual_cell.archetypes.ArchetypeBuild`.
            ValueError: If the build REFUSED — a Ring B/C archetype has no manifest, and forcing one
                into the current grammar is precisely what the missing-law verdict exists to stop.
        """
        if not isinstance(build, ArchetypeBuild):
            raise TypeError("build must be an ArchetypeBuild")
        if not build.supported or build.manifest is None:
            verdict = build.verdict
            laws = ", ".join(verdict.law_ids) if verdict is not None else "unnamed"
            raise ValueError(
                f"{build.archetype_id!r} has no manifest: the canonical grammar refused it, naming "
                f"missing laws [{laws}]. A posterior cannot be built by forcing it into the existing "
                f"component set; route to the missing-component/missing-law expansion instead."
            )
        envelope = ArtifactEnvelope(
            artifact_id=artifact_id,
            cell_state_manifest_hash=build.manifest.manifest_hash,
            evidence_source=evidence_source,
            representation=parameters.representation,
            source_artifact_ids=tuple(source_artifact_ids),
        )
        return cls(
            posterior_id=posterior_id,
            manifest=build.manifest,
            graph=graph,
            priors=build.priors,
            parameters=parameters,
            contract=contract,
            evidence=evidence,
            envelope=envelope,
            **optional,
        )
