"""Refusal harness that will gate a learned surrogate, exercised today on analytic oracles only.

Training a neural surrogate on simulation data is LOCKED: no accepted native trace exists, so any
model fitted now would be fitted to numbers that are not yet evidence.  What is *not* locked is the
machinery that must stand between such a model and a conclusion.  This module is that machinery, and
it is deliberately ``oracle-only``:

* the only fitted object here is a least-squares polynomial written with NumPy — it is not a cell
  model, it predicts no cell observable, and it exists so that the gates below can be validated
  against a target function whose truth is available in closed form;
* nothing in this module imports Warp, torch, or jax, and nothing reads a run record;
* every artifact it can emit is pinned to :data:`RepresentationKind.SURROGATE` /
  :data:`EvidenceSource.SURROGATE` with no parameter through which a caller could say otherwise.

The five gates, and the failure each one exists to catch:

1. **Gate authority ban** — a surrogate number closing a physics gate.  ``EvidenceSource`` is not a
   field of :class:`SurrogatePrediction`; it is a constant property that no caller can supply and no
   subclass may redefine, and — the guarantee that does not depend on the property at all — the
   emitted envelope re-derives ``SURROGATE`` rather than reading it.  See
   :class:`SurrogatePrediction` for which of those three is actually load-bearing; until the Codex
   external review of 2026-07-29 this bullet claimed the strongest of the three and held the
   weakest.  :class:`GateEvidenceLedger` closes the indirect route: an envelope that merely *cites*
   a surrogate ancestor is refused too, transitively.
2. **Trust region** — a confident extrapolation.  Outside the training box the harness returns a
   typed refusal carrying the offending axis and the normalized distance, never a number.
3. **Off-manifold rejection** — a query inside the bounding box but in a hole in the training
   support.  A box accepts these silently; a nearest-neighbour covering radius does not.
4. **Calibration** — an uncertainty that is decoration.  Coverage is measured against the nominal
   level with a binomial three-sigma test, and an uncalibrated model refuses every query.
5. **Held-out intervention** — the falsifier for "the surrogate learned the mechanism".  Small
   interpolation error inside the training distribution is not evidence; the number that matters is
   the error under an intervention setting never seen.

A refusal is not a smaller answer.  Every refusal emits a :class:`NativeQueryRequest` naming what to
evaluate natively and which lane the failure opens — native evaluation, solver/preconditioner, or
adaptive design — per the programme's expansionary failure routing.

Sanity Gate (recorded before first execution):

* **Dimensional** — all distances are computed in per-axis box-normalized coordinates, so the
  trust-region excess, the covering radius, and the nearest-neighbour distance are dimensionless and
  comparable across axes with different units.  A training axis of zero extent is refused at fit
  time rather than given an arbitrary scale.
* **Boundary** — a query exactly on a box face is inside (closed box, excess ``0.0``); a query one
  ULP outside is refused.  A covering-radius test is inclusive at the radius.  Fitting with fewer
  training rows than polynomial terms raises instead of returning a zero-residual illusion.
* **Conservation / invariant** — empirical coverage is a frequency and stays in ``[0, 1]``; the
  ledger's surrogate-ancestry closure is transitive and cycle-checked, so provenance cannot be
  laundered by an intermediate artifact.
* **Measurement protocol** — coverage is measured over a caller-supplied draw ensemble whose seed is
  recorded in the report; the null is ``coverage == nominal_level`` with binomial standard error
  ``sqrt(p(1-p)/n)`` and a declared two-sided three-sigma band.  The band is fixed here, before any
  data, and is never re-thresholded after seeing a result.
* **Numerical** — the design matrix is solved by SVD least squares (``rcond=None``); a rank-deficient
  design raises rather than returning a minimum-norm fit that would understate residual sigma.  The
  reported predictive sigma is ``residual_sigma * sqrt(1 + leverage)``, computed from a pseudo-inverse
  Gram matrix, because the residual sigma alone omits coefficient variance and is systematically
  overconfident for a new point.  Non-finite inputs, targets, or predictions raise.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from statistics import NormalDist
from types import MappingProxyType
from typing import Any

import numpy as np

from aleph.virtual_cell.contracts import (
    ApproximationErrorReport,
    ArtifactEnvelope,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    UncertaintySource,
)

__all__ = [
    "CoverageReport",
    "ExpansionRoute",
    "GateAuthorityError",
    "GateEvidenceLedger",
    "GatedSurrogate",
    "InterventionGeneralizationReport",
    "ManifoldSupport",
    "NativeQueryRequest",
    "PolynomialSurrogateModel",
    "SurrogateAuthorityError",
    "SurrogatePrediction",
    "SurrogateResponse",
    "SurrogateVerdict",
    "TrustRegionBox",
    "evaluate_held_out_intervention",
    "measure_interval_coverage",
]

# Declared before any data: a two-sided three-sigma band on the binomial null, and the same
# multiplier used to decide whether a held-out-intervention error refutes the surrogate's own
# stated uncertainty. It is a statistical convention, not a value chosen to make a check pass.
SIGMA_GATE: float = 3.0

_SURROGATE_UNCERTAINTY_SOURCES: tuple[UncertaintySource, ...] = (
    UncertaintySource.STRUCTURAL,
    UncertaintySource.PARAMETER,
    UncertaintySource.REPRESENTATION,
)

# The two names on SurrogatePrediction that state its evidence authority. They are properties over
# __slots__ with no field behind them, which stops a CALLER from supplying one -- and, until the
# Codex external review of 2026-07-29, stopped nothing else. The reviewer's counterexample:
#
#     class Forged(SurrogatePrediction):
#         @property
#         def evidence_source(self): return EvidenceSource.NATIVE_ACCEPTED
#
# built cleanly from the base dataclass fields, satisfied isinstance(), and reported
# `native-accepted`. The final envelope was still SURROGATE because as_artifact_envelope re-derives
# it, so nothing was laundered -- but the docstring claimed the value could not be spelled, and it
# could. _AuthorityLocked closes the subclass route at class creation and at post-hoc assignment.
_AUTHORITY_PROPERTIES: tuple[str, ...] = ("evidence_source", "representation_kind")


def _refuse_authority_override(owner: str, attribute: str) -> None:
    raise TypeError(
        f"{owner} may not redefine {attribute!r}: it states the evidence authority of a surrogate "
        f"value and is a constant property, not an overridable one (Codex external review "
        f"2026-07-29, finding 4)"
    )


class _AuthorityLocked(type):
    """Metaclass refusing post-hoc assignment onto an authority property.

    ``__init_subclass__`` catches a redefinition written in a subclass body; this catches
    ``Subclass.evidence_source = ...`` afterwards, which is the same override one line later.
    Neither is an absolute guarantee — an unbound ``type.__setattr__`` call bypasses any metaclass,
    and a duck-typed class that is not a subclass at all is outside both — which is why
    :meth:`SurrogatePrediction.as_artifact_envelope` re-derives the evidence source rather than
    trusting the property.
    """

    def __setattr__(cls, name: str, value: Any) -> None:
        if name in _AUTHORITY_PROPERTIES:
            _refuse_authority_override(f"{cls.__name__!r}", name)
        super().__setattr__(name, value)


class SurrogateVerdict(StrEnum):
    """Outcome of one gated surrogate query."""

    ANSWERED = "answered"
    REFUSED_OUTSIDE_TRUST_REGION = "refused-outside-trust-region"
    REFUSED_OFF_MANIFOLD = "refused-off-manifold"
    REFUSED_UNCALIBRATED = "refused-uncalibrated"
    REFUSED_INTERVENTION_UNSUPPORTED = "refused-intervention-unsupported"


class ExpansionRoute(StrEnum):
    """Which lane a surrogate failure opens; a failure never shrinks scope."""

    NATIVE_EVALUATION = "native-evaluation"
    SOLVER_PRECONDITIONER = "solver-preconditioner"
    ADAPTIVE_DESIGN = "adaptive-design"
    MODEL_INADEQUACY = "model-inadequacy"


_REFUSAL_EXPANSIONS: Mapping[SurrogateVerdict, tuple[ExpansionRoute, ...]] = MappingProxyType(
    {
        SurrogateVerdict.REFUSED_OUTSIDE_TRUST_REGION: (
            ExpansionRoute.NATIVE_EVALUATION,
            ExpansionRoute.ADAPTIVE_DESIGN,
        ),
        SurrogateVerdict.REFUSED_OFF_MANIFOLD: (
            ExpansionRoute.NATIVE_EVALUATION,
            ExpansionRoute.ADAPTIVE_DESIGN,
        ),
        SurrogateVerdict.REFUSED_UNCALIBRATED: (
            ExpansionRoute.NATIVE_EVALUATION,
            ExpansionRoute.SOLVER_PRECONDITIONER,
            ExpansionRoute.ADAPTIVE_DESIGN,
        ),
        SurrogateVerdict.REFUSED_INTERVENTION_UNSUPPORTED: (
            ExpansionRoute.NATIVE_EVALUATION,
            ExpansionRoute.ADAPTIVE_DESIGN,
            ExpansionRoute.MODEL_INADEQUACY,
        ),
    }
)


class GateAuthorityError(RuntimeError):
    """Raised when an artifact is asked to close a gate it has no authority to close."""


class SurrogateAuthorityError(GateAuthorityError):
    """Raised when surrogate provenance reaches a gate, directly or through ancestry."""


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _axis_names(axes: Sequence[str]) -> tuple[str, ...]:
    if isinstance(axes, str) or not isinstance(axes, Sequence):
        raise TypeError("axes must be a sequence of strings")
    names = tuple(_nonempty(axis, what="axis name") for axis in axes)
    if not names:
        raise ValueError("axes must not be empty")
    if len(set(names)) != len(names):
        raise ValueError("axis names must be unique")
    return names


def _finite_matrix(values: Any, *, n_axes: int, what: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if array.ndim != 2 or array.shape[0] < 1 or array.shape[1] != n_axes:
        raise ValueError(f"{what} must have shape (n_rows, {n_axes})")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{what} must be finite")
    return np.array(array, dtype=np.float64, copy=True)


def _finite_vector(values: Any, *, what: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64).ravel()
    if array.size < 1 or not np.all(np.isfinite(array)):
        raise ValueError(f"{what} must be a non-empty finite vector")
    return np.array(array, dtype=np.float64, copy=True)


def _readonly(array: np.ndarray) -> np.ndarray:
    frozen = np.array(array, dtype=np.float64, copy=True)
    frozen.setflags(write=False)
    return frozen


@dataclass(frozen=True, slots=True)
class TrustRegionBox:
    """Per-axis support box derived from training inputs; it is measured, never declared."""

    axes: tuple[str, ...]
    lower: np.ndarray
    upper: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "axes", _axis_names(self.axes))
        lower = _finite_vector(self.lower, what="lower")
        upper = _finite_vector(self.upper, what="upper")
        if lower.size != len(self.axes) or upper.size != len(self.axes):
            raise ValueError("lower and upper must have one entry per axis")
        if np.any(upper <= lower):
            raise ValueError(
                "every trust-region axis needs a strictly positive extent; a constant training "
                "axis has no measurable scale and must not be silently rescaled"
            )
        object.__setattr__(self, "lower", _readonly(lower))
        object.__setattr__(self, "upper", _readonly(upper))

    @classmethod
    def from_training(cls, axes: Sequence[str], inputs: Any) -> TrustRegionBox:
        """Fit the box to the observed extent of the training inputs.

        Args:
            axes: Ordered input axis names.
            inputs: Training inputs of shape ``(n_rows, len(axes))``.

        Returns:
            The tightest axis-aligned box containing every training row.
        """
        names = _axis_names(axes)
        matrix = _finite_matrix(inputs, n_axes=len(names), what="training inputs")
        return cls(axes=names, lower=matrix.min(axis=0), upper=matrix.max(axis=0))

    @property
    def extent(self) -> np.ndarray:
        """Per-axis width used as the normalization scale for every distance in this module."""
        return np.asarray(self.upper) - np.asarray(self.lower)

    def normalize(self, points: Any) -> np.ndarray:
        """Map points into per-axis box-normalized coordinates where the box is the unit cube."""
        matrix = _finite_matrix(points, n_axes=len(self.axes), what="points")
        return (matrix - np.asarray(self.lower)[None, :]) / self.extent[None, :]

    def normalized_excess(self, point: Any) -> np.ndarray:
        """Return the per-axis distance outside the closed box, normalized by the axis extent."""
        row = self.normalize(point)
        if row.shape[0] != 1:
            raise ValueError("normalized_excess takes a single point")
        values = row[0]
        return np.maximum(np.maximum(-values, values - 1.0), 0.0)

    def contains(self, point: Any) -> bool:
        """Return whether the point lies in the closed box."""
        return bool(np.all(self.normalized_excess(point) == 0.0))

    def distance_outside(self, point: Any) -> float:
        """Return the Euclidean norm of the normalized per-axis excess (0.0 when inside)."""
        return float(np.linalg.norm(self.normalized_excess(point)))

    def worst_axis(self, point: Any) -> tuple[str, float]:
        """Return the axis with the largest normalized excess and that excess."""
        excess = self.normalized_excess(point)
        index = int(np.argmax(excess))
        return self.axes[index], float(excess[index])

    def summary(self) -> Mapping[str, float]:
        """Return a JSON-scalar description suitable for a ``RepresentationDescriptor``."""
        entries: dict[str, float] = {}
        for index, axis in enumerate(self.axes):
            entries[f"{axis}:lower"] = float(self.lower[index])
            entries[f"{axis}:upper"] = float(self.upper[index])
        return MappingProxyType(entries)


@dataclass(frozen=True, slots=True)
class ManifoldSupport:
    """Nearest-neighbour support test for queries a bounding box would accept by mistake.

    The acceptance radius is the sample's own covering radius at the declared neighbour order — the
    largest leave-one-out ``k``-th nearest-neighbour distance among the training rows, in
    box-normalized coordinates.  It is a statistic of the training set, so it introduces no tunable
    threshold: a query is off-manifold exactly when it is further from the training set than the
    training set is from itself.
    """

    box: TrustRegionBox
    normalized_points: np.ndarray
    neighbour_order: int
    covering_radius: float

    def __post_init__(self) -> None:
        if not isinstance(self.box, TrustRegionBox):
            raise TypeError("box must be a TrustRegionBox")
        points = np.asarray(self.normalized_points, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != len(self.box.axes) or points.shape[0] < 2:
            raise ValueError("normalized_points must be (n_rows >= 2, n_axes)")
        if not np.all(np.isfinite(points)):
            raise ValueError("normalized_points must be finite")
        if (
            isinstance(self.neighbour_order, bool)
            or not isinstance(self.neighbour_order, int)
            or self.neighbour_order < 1
            or self.neighbour_order >= points.shape[0]
        ):
            raise ValueError("neighbour_order must be a positive integer below the training count")
        if not math.isfinite(self.covering_radius) or self.covering_radius <= 0.0:
            raise ValueError("covering_radius must be finite and positive")
        object.__setattr__(self, "normalized_points", _readonly(points))

    @classmethod
    def from_training(
        cls,
        box: TrustRegionBox,
        inputs: Any,
        *,
        neighbour_order: int = 1,
    ) -> ManifoldSupport:
        """Measure the training set's covering radius at ``neighbour_order``.

        Args:
            box: Trust region supplying the per-axis normalization.
            inputs: Training inputs of shape ``(n_rows, n_axes)``.
            neighbour_order: ``k`` in the leave-one-out ``k``-th nearest-neighbour distance.

        Returns:
            A fitted support detector.

        Raises:
            ValueError: If the training rows contain a duplicate-collapsed set too small for ``k``.
        """
        normalized = box.normalize(inputs)
        distances = np.linalg.norm(normalized[:, None, :] - normalized[None, :, :], axis=-1)
        np.fill_diagonal(distances, np.inf)
        ordered = np.sort(distances, axis=1)
        if neighbour_order > ordered.shape[1] - 1:
            raise ValueError("neighbour_order exceeds the available leave-one-out neighbours")
        kth = ordered[:, neighbour_order - 1]
        radius = float(np.max(kth))
        if not math.isfinite(radius) or radius <= 0.0:
            raise ValueError("training rows collapse to a single point; support is not measurable")
        return cls(
            box=box,
            normalized_points=normalized,
            neighbour_order=int(neighbour_order),
            covering_radius=radius,
        )

    def distance(self, point: Any) -> float:
        """Return the query's ``k``-th nearest-neighbour distance in normalized coordinates."""
        row = self.box.normalize(point)
        if row.shape[0] != 1:
            raise ValueError("distance takes a single point")
        gaps = np.linalg.norm(np.asarray(self.normalized_points) - row[0][None, :], axis=-1)
        ordered = np.sort(gaps)
        return float(ordered[self.neighbour_order - 1])

    def supports(self, point: Any) -> bool:
        """Return whether the query lies within the measured covering radius (inclusive)."""
        return self.distance(point) <= self.covering_radius


@dataclass(frozen=True, slots=True)
class CoverageReport:
    """Measured interval coverage against a nominal level, with the binomial null test."""

    label: str
    nominal_level: float
    empirical_coverage: float
    n_draws: int
    standard_error: float
    z_score: float
    sigma_gate: float
    seed: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "label", _nonempty(self.label, what="label"))
        if not 0.0 < self.nominal_level < 1.0:
            raise ValueError("nominal_level must lie strictly inside (0, 1)")
        if not 0.0 <= self.empirical_coverage <= 1.0:
            raise ValueError("empirical_coverage must be a frequency in [0, 1]")
        if isinstance(self.n_draws, bool) or not isinstance(self.n_draws, int) or self.n_draws < 2:
            raise ValueError("n_draws must be an integer of at least two")
        if isinstance(self.seed, bool) or not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError("seed must be a nonnegative integer and must be recorded")
        for name in ("standard_error", "z_score", "sigma_gate"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")

    @property
    def calibrated(self) -> bool:
        """Whether the measured coverage is inside the declared binomial sigma band."""
        return abs(self.z_score) <= self.sigma_gate

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able coverage record."""
        return {
            "label": self.label,
            "nominal_level": self.nominal_level,
            "empirical_coverage": self.empirical_coverage,
            "n_draws": self.n_draws,
            "standard_error": self.standard_error,
            "z_score": self.z_score,
            "sigma_gate": self.sigma_gate,
            "seed": self.seed,
            "calibrated": self.calibrated,
        }


def measure_interval_coverage(
    truths: Any,
    predicted_means: Any,
    predicted_sigmas: Any,
    *,
    nominal_level: float,
    seed: int,
    label: str,
) -> CoverageReport:
    """Measure how often a nominal central interval actually contains the truth.

    The intervals are the central Gaussian intervals implied by the reported sigmas.  Nothing is
    fitted here and nothing is repaired: the function reports the frequency and the binomial
    ``z``-score against the null ``coverage == nominal_level``.

    Args:
        truths: Ground-truth draws from the generative process.
        predicted_means: Reported interval centres, same length as ``truths``.
        predicted_sigmas: Reported standard deviations, strictly positive, same length.
        nominal_level: The claimed central probability, e.g. ``0.95``.
        seed: The RNG seed of the draw ensemble; recorded so the number is reproducible.
        label: Human-readable name of the model whose intervals are being measured.

    Returns:
        The coverage report; :attr:`CoverageReport.calibrated` is the fail-closed verdict.

    Raises:
        ValueError: If the arrays disagree in length, are not finite, or any sigma is nonpositive.
    """
    truth = _finite_vector(truths, what="truths")
    mean = _finite_vector(predicted_means, what="predicted_means")
    sigma = _finite_vector(predicted_sigmas, what="predicted_sigmas")
    if truth.size != mean.size or truth.size != sigma.size:
        raise ValueError("truths, predicted_means, and predicted_sigmas must have equal length")
    if truth.size < 2:
        raise ValueError("coverage needs at least two draws")
    if np.any(sigma <= 0.0):
        raise ValueError("predicted_sigmas must be strictly positive")
    if not 0.0 < nominal_level < 1.0:
        raise ValueError("nominal_level must lie strictly inside (0, 1)")

    half_width = NormalDist().inv_cdf(0.5 + nominal_level / 2.0)
    covered = np.abs(truth - mean) <= half_width * sigma
    empirical = float(np.mean(covered))
    n_draws = int(truth.size)
    standard_error = math.sqrt(nominal_level * (1.0 - nominal_level) / n_draws)
    return CoverageReport(
        label=label,
        nominal_level=float(nominal_level),
        empirical_coverage=empirical,
        n_draws=n_draws,
        standard_error=standard_error,
        z_score=(empirical - float(nominal_level)) / standard_error,
        sigma_gate=SIGMA_GATE,
        seed=int(seed),
    )


@dataclass(frozen=True, slots=True)
class PolynomialSurrogateModel:
    """A deliberately imperfect analytic stand-in for a future learned surrogate.

    A total-degree polynomial fitted by SVD least squares.  It is here so the gates can be validated
    against a target whose truth is known in closed form; it is not a model of anything biological
    and must never be described as one.
    """

    axes: tuple[str, ...]
    degree: int
    exponents: tuple[tuple[int, ...], ...]
    coefficients: np.ndarray
    gram_inverse: np.ndarray
    residual_sigma: float
    n_training: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "axes", _axis_names(self.axes))
        if isinstance(self.degree, bool) or not isinstance(self.degree, int) or self.degree < 1:
            raise ValueError("degree must be a positive integer")
        coefficients = _finite_vector(self.coefficients, what="coefficients")
        if coefficients.size != len(self.exponents):
            raise ValueError("coefficients and exponents must agree in length")
        gram_inverse = np.asarray(self.gram_inverse, dtype=np.float64)
        if gram_inverse.shape != (coefficients.size, coefficients.size) or not np.all(
            np.isfinite(gram_inverse)
        ):
            raise ValueError("gram_inverse must be a finite (n_terms, n_terms) matrix")
        if not math.isfinite(self.residual_sigma) or self.residual_sigma <= 0.0:
            raise ValueError("residual_sigma must be finite and positive")
        if (
            isinstance(self.n_training, bool)
            or not isinstance(self.n_training, int)
            or self.n_training < 1
        ):
            raise ValueError("n_training must be a positive integer")
        object.__setattr__(self, "coefficients", _readonly(coefficients))
        object.__setattr__(self, "gram_inverse", _readonly(gram_inverse))

    @staticmethod
    def _exponent_table(n_axes: int, degree: int) -> tuple[tuple[int, ...], ...]:
        table: list[tuple[int, ...]] = []
        for total in range(degree + 1):
            for combination in itertools.combinations_with_replacement(range(n_axes), total):
                exponents = [0] * n_axes
                for axis in combination:
                    exponents[axis] += 1
                table.append(tuple(exponents))
        return tuple(table)

    def _design(self, inputs: np.ndarray) -> np.ndarray:
        columns = [
            np.prod(inputs ** np.asarray(exponent, dtype=np.float64)[None, :], axis=1)
            for exponent in self.exponents
        ]
        return np.stack(columns, axis=1)

    @classmethod
    def fit(
        cls,
        axes: Sequence[str],
        inputs: Any,
        targets: Any,
        *,
        degree: int,
    ) -> PolynomialSurrogateModel:
        """Fit the polynomial and measure its own residual standard error.

        Args:
            axes: Ordered input axis names.
            inputs: Training inputs of shape ``(n_rows, len(axes))``.
            targets: Training targets of length ``n_rows``.
            degree: Total polynomial degree.

        Returns:
            The fitted model.

        Raises:
            ValueError: If the design matrix is rank deficient or has no residual degrees of
                freedom — either case would report an unearned residual sigma.
        """
        names = _axis_names(axes)
        matrix = _finite_matrix(inputs, n_axes=len(names), what="training inputs")
        values = _finite_vector(targets, what="training targets")
        if values.size != matrix.shape[0]:
            raise ValueError("training inputs and targets must have equal row counts")
        if isinstance(degree, bool) or not isinstance(degree, int) or degree < 1:
            raise ValueError("degree must be a positive integer")

        exponents = cls._exponent_table(len(names), degree)
        n_terms = len(exponents)
        if matrix.shape[0] <= n_terms:
            raise ValueError(
                f"a degree-{degree} fit over {len(names)} axes needs more than {n_terms} training "
                "rows to leave residual degrees of freedom"
            )
        design = np.stack(
            [
                np.prod(matrix ** np.asarray(exponent, dtype=np.float64)[None, :], axis=1)
                for exponent in exponents
            ],
            axis=1,
        )
        coefficients, _, rank, _ = np.linalg.lstsq(design, values, rcond=None)
        if int(rank) != n_terms:
            raise ValueError(
                "the polynomial design matrix is rank deficient; a minimum-norm fit would "
                "understate the residual sigma"
            )
        pseudo_inverse = np.linalg.pinv(design)
        gram_inverse = pseudo_inverse @ pseudo_inverse.T
        residuals = values - design @ coefficients
        dof = matrix.shape[0] - n_terms
        sigma = float(math.sqrt(float(residuals @ residuals) / dof))
        if not math.isfinite(sigma) or sigma <= 0.0:
            raise ValueError(
                "the fit reproduces the training targets exactly; an exactly-zero residual sigma "
                "is not a usable uncertainty statement"
            )
        return cls(
            axes=names,
            degree=int(degree),
            exponents=exponents,
            coefficients=coefficients,
            gram_inverse=gram_inverse,
            residual_sigma=sigma,
            n_training=int(matrix.shape[0]),
        )

    def predict(self, inputs: Any) -> np.ndarray:
        """Return raw polynomial values with no gate applied.

        This is the ungated numerical core.  Callers outside this module should use
        :meth:`GatedSurrogate.query`; the only legitimate direct use is diagnosis of what a gate
        prevented, and such a value may never enter an artifact.
        """
        matrix = _finite_matrix(inputs, n_axes=len(self.axes), what="inputs")
        values = self._design(matrix) @ np.asarray(self.coefficients)
        if not np.all(np.isfinite(values)):
            raise ValueError("polynomial evaluation produced a non-finite value")
        return values

    def prediction_sigma(self, inputs: Any) -> np.ndarray:
        """Return the per-point predictive standard deviation ``sigma * sqrt(1 + leverage)``.

        The residual sigma alone describes the spread of a fitted point, not of a NEW observation:
        it omits the variance of the fitted coefficients themselves.  Reporting the residual sigma
        as a predictive sigma is a small, systematic overconfidence — a few percent here — and the
        calibration gate detects it, so the correct statistic is used rather than the convenient one.
        """
        matrix = _finite_matrix(inputs, n_axes=len(self.axes), what="inputs")
        design = self._design(matrix)
        leverage = np.einsum("ij,jk,ik->i", design, np.asarray(self.gram_inverse), design)
        if np.any(leverage < 0.0) or not np.all(np.isfinite(leverage)):
            raise ValueError("leverage evaluation produced a non-finite or negative value")
        return self.residual_sigma * np.sqrt(1.0 + leverage)


@dataclass(frozen=True, slots=True)
class NativeQueryRequest:
    """A typed request for a native evaluation, emitted instead of a guess.

    The request names the point to evaluate, why the surrogate refused, and which lanes the refusal
    opens.  It is the acquisition half of the harness: a refusal that produced no request would be
    indistinguishable from a dropped question.
    """

    request_id: str
    reason: SurrogateVerdict
    rationale: str
    axes: tuple[str, ...]
    query_point: tuple[float, ...]
    diagnostic: Mapping[str, float]
    expansions: tuple[ExpansionRoute, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "request_id", _nonempty(self.request_id, what="request_id"))
        object.__setattr__(self, "rationale", _nonempty(self.rationale, what="rationale"))
        if not isinstance(self.reason, SurrogateVerdict):
            raise TypeError("reason must be a SurrogateVerdict")
        if self.reason is SurrogateVerdict.ANSWERED:
            raise ValueError("an answered query does not request a native evaluation")
        object.__setattr__(self, "axes", _axis_names(self.axes))
        point = tuple(float(value) for value in self.query_point)
        if len(point) != len(self.axes) or any(not math.isfinite(value) for value in point):
            raise ValueError("query_point must hold one finite value per axis")
        object.__setattr__(self, "query_point", point)
        expansions = tuple(self.expansions)
        if not expansions or any(not isinstance(route, ExpansionRoute) for route in expansions):
            raise ValueError("a refusal must open at least one ExpansionRoute")
        if len(set(expansions)) != len(expansions):
            raise ValueError("expansions must not contain duplicates")
        object.__setattr__(self, "expansions", expansions)
        diagnostic: dict[str, float] = {}
        for key, value in dict(self.diagnostic).items():
            name = _nonempty(key, what="diagnostic key")
            number = float(value)
            if not math.isfinite(number):
                raise ValueError(f"diagnostic[{name!r}] must be finite")
            diagnostic[name] = number
        object.__setattr__(self, "diagnostic", MappingProxyType(diagnostic))

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able acquisition record."""
        return {
            "request_id": self.request_id,
            "reason": str(self.reason),
            "rationale": self.rationale,
            "axes": list(self.axes),
            "query_point": list(self.query_point),
            "diagnostic": dict(sorted(self.diagnostic.items())),
            "expansions": [str(route) for route in self.expansions],
        }

    @property
    def request_hash(self) -> str:
        """SHA-256 of the canonical acquisition record."""
        payload = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode()
        return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True, slots=True)
class SurrogatePrediction(metaclass=_AuthorityLocked):
    """One accepted surrogate value, permanently marked as surrogate evidence.

    Three mechanisms, of decreasing generality, and the docstring names them separately because the
    Codex external review of 2026-07-29 (finding 4) found this class claiming the strongest of the
    three when it only had the weakest:

    1. **No field, no keyword.**  There is no ``evidence_source`` field; :attr:`evidence_source` and
       :attr:`representation_kind` are properties over ``__slots__`` returning constants, so a
       CALLER cannot supply one and :func:`dataclasses.replace` raises.
    2. **No override.**  :meth:`__init_subclass__` refuses a subclass that redefines either property
       in its body, and :class:`_AuthorityLocked` refuses the same redefinition assigned onto the
       subclass afterwards.  This was the hole: before the review, ``class Forged(SurrogatePrediction)``
       with an ``evidence_source`` returning ``NATIVE_ACCEPTED`` was created without complaint and
       satisfied ``isinstance``.  Neither guard is absolute — an unbound ``type.__setattr__`` call
       bypasses a metaclass, and a duck-typed class that never subclasses this one is outside both.
    3. **The envelope does not consult either.**  :meth:`as_artifact_envelope` re-derives
       :data:`EvidenceSource.SURROGATE` and :data:`RepresentationKind.SURROGATE` from constants
       rather than reading the properties, so the provenance actually recorded is pinned even if
       (2) is defeated.  This is the guarantee to rely on; (2) makes the forged type harder to
       bring into existence, and (1) makes the honest path the only convenient one.

    What was corrected in the wording: the previous docstring said self-certification as a native
    result was "a value that cannot be spelled".  It could be spelled, by subclassing.  It could
    not be RECORDED, which is a different and weaker sentence, and is the one this class is entitled
    to make.
    """

    surrogate_id: str
    axes: tuple[str, ...]
    query_point: tuple[float, ...]
    mean: float
    predictive_sigma: float
    trust_region: Mapping[str, float]
    support_distance: float
    coverage: CoverageReport

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Refuse a subclass that redefines an authority property in its own body.

        Motivating counterexample (Codex external review 2026-07-29, finding 4): a subclass whose
        ``evidence_source`` property returned :data:`EvidenceSource.NATIVE_ACCEPTED` passed
        ``isinstance`` and reported native provenance.  Subclassing for any other reason stays legal.

        Args:
            **kwargs: Forwarded to :meth:`type.__init_subclass__`.

        Raises:
            TypeError: If the subclass body defines ``evidence_source`` or ``representation_kind``.
        """
        # Explicit two-argument super(): @dataclass(slots=True) REPLACES this class object after the
        # body executes, so the zero-argument form's __class__ cell points at a class that is no
        # longer in any subclass's MRO and raises instead of guarding.
        super(SurrogatePrediction, cls).__init_subclass__(**kwargs)
        for attribute in _AUTHORITY_PROPERTIES:
            if attribute in cls.__dict__:
                _refuse_authority_override(f"subclass {cls.__name__!r}", attribute)

    def __post_init__(self) -> None:
        object.__setattr__(self, "surrogate_id", _nonempty(self.surrogate_id, what="surrogate_id"))
        object.__setattr__(self, "axes", _axis_names(self.axes))
        point = tuple(float(value) for value in self.query_point)
        if len(point) != len(self.axes) or any(not math.isfinite(value) for value in point):
            raise ValueError("query_point must hold one finite value per axis")
        object.__setattr__(self, "query_point", point)
        for name in ("mean", "predictive_sigma", "support_distance"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        if self.predictive_sigma <= 0.0:
            raise ValueError("predictive_sigma must be positive")
        if self.support_distance < 0.0:
            raise ValueError("support_distance must be nonnegative")
        if not isinstance(self.coverage, CoverageReport):
            raise TypeError("coverage must be a CoverageReport")
        if not self.coverage.calibrated:
            raise ValueError("an uncalibrated surrogate must not emit a prediction")
        trust_region = {
            _nonempty(key, what="trust_region key"): float(value)
            for key, value in dict(self.trust_region).items()
        }
        if not trust_region or any(not math.isfinite(value) for value in trust_region.values()):
            raise ValueError("trust_region must be non-empty and finite")
        object.__setattr__(self, "trust_region", MappingProxyType(trust_region))

    @property
    def evidence_source(self) -> EvidenceSource:
        """Always :data:`EvidenceSource.SURROGATE`; there is no field behind this property."""
        return EvidenceSource.SURROGATE

    @property
    def representation_kind(self) -> RepresentationKind:
        """Always :data:`RepresentationKind.SURROGATE`."""
        return RepresentationKind.SURROGATE

    def as_artifact_envelope(
        self,
        *,
        artifact_id: str,
        cell_state_manifest_hash: str,
        source_artifact_ids: Sequence[str],
        error_report: ApproximationErrorReport,
        fallback: str = "native-evaluation",
    ) -> ArtifactEnvelope:
        """Wrap the prediction in a provenance envelope pinned to surrogate evidence.

        Args:
            artifact_id: Identifier for the emitted artifact.
            cell_state_manifest_hash: Manifest the prediction is about.
            source_artifact_ids: Ancestry; required, because a surrogate is always derived.
            error_report: Measured held-out error binding source and reduced blobs.
            fallback: What runs when the surrogate refuses.

        Returns:
            An envelope whose ``evidence_source`` is :data:`EvidenceSource.SURROGATE`.
        """
        if not isinstance(error_report, ApproximationErrorReport):
            raise TypeError("error_report must be an ApproximationErrorReport")
        representation = RepresentationDescriptor(
            kind=RepresentationKind.SURROGATE,
            uncertainty_sources=_SURROGATE_UNCERTAINTY_SOURCES,
            axes=self.axes,
            approximation_error_bound=error_report.measured_error,
            error_metric=error_report.metric_id,
            error_report=error_report,
            fallback=_nonempty(fallback, what="fallback"),
            trust_region=dict(self.trust_region),
        )
        return ArtifactEnvelope(
            artifact_id=artifact_id,
            cell_state_manifest_hash=cell_state_manifest_hash,
            evidence_source=EvidenceSource.SURROGATE,
            representation=representation,
            source_artifact_ids=tuple(source_artifact_ids),
        )


@dataclass(frozen=True, slots=True)
class SurrogateResponse:
    """Exactly one of an accepted prediction or a typed native-evaluation request."""

    verdict: SurrogateVerdict
    prediction: SurrogatePrediction | None = None
    native_request: NativeQueryRequest | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.verdict, SurrogateVerdict):
            raise TypeError("verdict must be a SurrogateVerdict")
        answered = self.verdict is SurrogateVerdict.ANSWERED
        if answered and (self.prediction is None or self.native_request is not None):
            raise ValueError("an answered response carries a prediction and no native request")
        if not answered and (self.prediction is not None or self.native_request is None):
            raise ValueError("a refusal carries a native request and no prediction")
        if self.native_request is not None and self.native_request.reason is not self.verdict:
            raise ValueError("native_request reason must match the response verdict")

    @property
    def answered(self) -> bool:
        """Whether the harness produced a usable number."""
        return self.verdict is SurrogateVerdict.ANSWERED

    @property
    def value(self) -> float:
        """The predicted mean.

        Raises:
            SurrogateAuthorityError: If the query was refused; a refusal has no number, and reading
                one is exactly the degradation this harness exists to prevent.
        """
        if self.prediction is None:
            raise SurrogateAuthorityError(
                f"surrogate refused ({self.verdict}); a refusal has no value — consume "
                "native_request instead"
            )
        return self.prediction.mean


@dataclass(frozen=True, slots=True)
class GatedSurrogate:
    """A fitted model behind the trust-region, off-manifold, and calibration gates."""

    surrogate_id: str
    model: PolynomialSurrogateModel
    box: TrustRegionBox
    support: ManifoldSupport
    coverage: CoverageReport

    def __post_init__(self) -> None:
        object.__setattr__(self, "surrogate_id", _nonempty(self.surrogate_id, what="surrogate_id"))
        if not isinstance(self.model, PolynomialSurrogateModel):
            raise TypeError("model must be a PolynomialSurrogateModel")
        if not isinstance(self.box, TrustRegionBox):
            raise TypeError("box must be a TrustRegionBox")
        if not isinstance(self.support, ManifoldSupport):
            raise TypeError("support must be a ManifoldSupport")
        if not isinstance(self.coverage, CoverageReport):
            raise TypeError("coverage must be a CoverageReport")
        if self.model.axes != self.box.axes or self.support.box.axes != self.box.axes:
            raise ValueError("model, trust region, and support must share one axis order")

    @classmethod
    def fit(
        cls,
        *,
        surrogate_id: str,
        axes: Sequence[str],
        inputs: Any,
        targets: Any,
        degree: int,
        coverage: CoverageReport,
        neighbour_order: int = 1,
    ) -> GatedSurrogate:
        """Fit the model and derive both support tests from the same training inputs.

        Args:
            surrogate_id: Identifier carried into every prediction and request.
            axes: Ordered input axis names.
            inputs: Training inputs of shape ``(n_rows, len(axes))``.
            targets: Training targets of length ``n_rows``.
            degree: Total polynomial degree.
            coverage: A previously MEASURED coverage report; calibration is never assumed.
            neighbour_order: ``k`` for the covering-radius support test.

        Returns:
            The gated surrogate.
        """
        names = _axis_names(axes)
        matrix = _finite_matrix(inputs, n_axes=len(names), what="training inputs")
        box = TrustRegionBox.from_training(names, matrix)
        return cls(
            surrogate_id=surrogate_id,
            model=PolynomialSurrogateModel.fit(names, matrix, targets, degree=degree),
            box=box,
            support=ManifoldSupport.from_training(box, matrix, neighbour_order=neighbour_order),
            coverage=coverage,
        )

    @property
    def axes(self) -> tuple[str, ...]:
        """Ordered input axis names."""
        return self.box.axes

    def ungated_prediction_for_diagnosis(self, point: Any) -> float:
        """Return what the model WOULD have said, ignoring every gate.

        Diagnostic only.  The value is never wrapped in a prediction or an envelope; its sole
        purpose is to let a test or a report quantify the error a refusal avoided.
        """
        return float(self.model.predict(point)[0])

    def _refuse(
        self,
        *,
        verdict: SurrogateVerdict,
        point: tuple[float, ...],
        rationale: str,
        diagnostic: Mapping[str, float],
    ) -> SurrogateResponse:
        request = NativeQueryRequest(
            request_id=f"{self.surrogate_id}:{verdict}",
            reason=verdict,
            rationale=rationale,
            axes=self.axes,
            query_point=point,
            diagnostic=diagnostic,
            expansions=_REFUSAL_EXPANSIONS[verdict],
        )
        return SurrogateResponse(verdict=verdict, native_request=request)

    def query(self, point: Any) -> SurrogateResponse:
        """Answer or refuse one query, fail-closed and in a fixed gate order.

        The order is calibration, then trust region, then manifold support: an uncalibrated model has
        nothing worth localizing, and a point outside the box has no meaningful nearest-neighbour
        statistic relative to the training support.

        Args:
            point: One query of length ``len(self.axes)``.

        Returns:
            A response carrying either a prediction or a typed native-evaluation request.
        """
        row = _finite_matrix(point, n_axes=len(self.axes), what="query point")
        if row.shape[0] != 1:
            raise ValueError("query takes a single point")
        values = tuple(float(value) for value in row[0])

        if not self.coverage.calibrated:
            return self._refuse(
                verdict=SurrogateVerdict.REFUSED_UNCALIBRATED,
                point=values,
                rationale=(
                    f"measured {self.coverage.nominal_level:.0%} interval coverage was "
                    f"{self.coverage.empirical_coverage:.4f} "
                    f"({self.coverage.z_score:+.1f} sigma from nominal); the surrogate's own "
                    "uncertainty statement is refuted, so no query may be answered"
                ),
                diagnostic={
                    "nominal_level": self.coverage.nominal_level,
                    "empirical_coverage": self.coverage.empirical_coverage,
                    "z_score": self.coverage.z_score,
                },
            )

        excess = self.box.distance_outside(values)
        if excess > 0.0:
            axis, axis_excess = self.box.worst_axis(values)
            index = self.axes.index(axis)
            return self._refuse(
                verdict=SurrogateVerdict.REFUSED_OUTSIDE_TRUST_REGION,
                point=values,
                rationale=(
                    f"axis {axis!r} is {axis_excess:.4g} box-widths outside the training extent "
                    f"[{float(self.box.lower[index]):.6g}, {float(self.box.upper[index]):.6g}]; "
                    "extrapolation is refused rather than reported"
                ),
                diagnostic={
                    "normalized_excess": excess,
                    "worst_axis_excess": axis_excess,
                    "axis_lower": float(self.box.lower[index]),
                    "axis_upper": float(self.box.upper[index]),
                },
            )

        distance = self.support.distance(values)
        if distance > self.support.covering_radius:
            return self._refuse(
                verdict=SurrogateVerdict.REFUSED_OFF_MANIFOLD,
                point=values,
                rationale=(
                    f"the query is {distance:.4g} normalized units from its neighbour of order "
                    f"{self.support.neighbour_order}, beyond the training covering radius "
                    f"{self.support.covering_radius:.4g}; it lies inside the bounding box but off "
                    "the training support"
                ),
                diagnostic={
                    "neighbour_distance": distance,
                    "covering_radius": self.support.covering_radius,
                    "normalized_excess": excess,
                },
            )

        return SurrogateResponse(
            verdict=SurrogateVerdict.ANSWERED,
            prediction=SurrogatePrediction(
                surrogate_id=self.surrogate_id,
                axes=self.axes,
                query_point=values,
                mean=self.ungated_prediction_for_diagnosis(values),
                predictive_sigma=float(self.model.prediction_sigma([values])[0]),
                trust_region=self.box.summary(),
                support_distance=distance,
                coverage=self.coverage,
            ),
        )


@dataclass(frozen=True, slots=True)
class InterventionGeneralizationReport:
    """Interpolation error versus held-out-intervention error, and the mechanism falsifier.

    Interpolating inside the training distribution is not evidence of mechanism: a model that has
    only ever seen intervention settings under which a term is dormant will interpolate them
    perfectly and be wrong the moment that term wakes.  ``mechanism_falsified`` compares the
    held-out error with the surrogate's OWN stated sigma, so the criterion is the model's claim, not
    a threshold chosen after the fact.
    """

    surrogate_id: str
    intervention_axis: str
    training_interventions: tuple[float, ...]
    held_out_intervention: float
    interpolation_rmse: float
    held_out_rmse: float
    predictive_sigma: float
    sigma_gate: float
    seed: int
    native_request: NativeQueryRequest | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "surrogate_id", _nonempty(self.surrogate_id, what="surrogate_id"))
        object.__setattr__(
            self, "intervention_axis", _nonempty(self.intervention_axis, what="intervention_axis")
        )
        interventions = tuple(float(value) for value in self.training_interventions)
        if not interventions or any(not math.isfinite(value) for value in interventions):
            raise ValueError("training_interventions must be non-empty and finite")
        object.__setattr__(self, "training_interventions", interventions)
        for name in ("held_out_intervention", "interpolation_rmse", "held_out_rmse"):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"{name} must be finite")
        if self.interpolation_rmse < 0.0 or self.held_out_rmse < 0.0:
            raise ValueError("errors must be nonnegative")
        if self.predictive_sigma <= 0.0:
            raise ValueError("predictive_sigma must be positive")
        if isinstance(self.seed, bool) or not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError("seed must be a nonnegative integer and must be recorded")

    @property
    def mechanism_falsified(self) -> bool:
        """Whether the held-out error exceeds the surrogate's own ``sigma_gate``-sigma claim."""
        return self.held_out_rmse > self.sigma_gate * self.predictive_sigma

    @property
    def degradation_ratio(self) -> float:
        """Held-out RMSE divided by interpolation RMSE; ``inf`` when interpolation is exact."""
        if self.interpolation_rmse == 0.0:
            return math.inf
        return self.held_out_rmse / self.interpolation_rmse

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able generalization record."""
        return {
            "surrogate_id": self.surrogate_id,
            "intervention_axis": self.intervention_axis,
            "training_interventions": list(self.training_interventions),
            "held_out_intervention": self.held_out_intervention,
            "interpolation_rmse": self.interpolation_rmse,
            "held_out_rmse": self.held_out_rmse,
            "predictive_sigma": self.predictive_sigma,
            "sigma_gate": self.sigma_gate,
            "degradation_ratio": self.degradation_ratio,
            "mechanism_falsified": self.mechanism_falsified,
            "seed": self.seed,
        }


def evaluate_held_out_intervention(
    surrogate: GatedSurrogate,
    *,
    intervention_axis: str,
    training_interventions: Sequence[float],
    held_out_intervention: float,
    interpolation_inputs: Any,
    interpolation_truth: Any,
    held_out_inputs: Any,
    held_out_truth: Any,
    seed: int,
) -> InterventionGeneralizationReport:
    """Measure interpolation error and held-out-intervention error on the same surrogate.

    Both errors are measured with the UNGATED model on purpose.  The gate would refuse the held-out
    setting, which is the correct operational behaviour but tells us nothing about how wrong the
    model was; this function exists to put that number on the record.

    Args:
        surrogate: The fitted, gated surrogate.
        intervention_axis: Name of the axis whose setting is held out.
        training_interventions: Intervention settings present in training.
        held_out_intervention: The setting never seen in training.
        interpolation_inputs: Held-out inputs drawn at TRAINING intervention settings.
        interpolation_truth: Closed-form truth for ``interpolation_inputs``.
        held_out_inputs: Inputs at the held-out intervention setting.
        held_out_truth: Closed-form truth for ``held_out_inputs``.
        seed: Seed of the evaluation draws; recorded so both numbers are reproducible.

    Returns:
        The report; :attr:`InterventionGeneralizationReport.mechanism_falsified` is the verdict.

    Raises:
        ValueError: If ``intervention_axis`` is not an axis of the surrogate, or if the held-out
            setting also appears in ``training_interventions``.
    """
    if not isinstance(surrogate, GatedSurrogate):
        raise TypeError("surrogate must be a GatedSurrogate")
    axis = _nonempty(intervention_axis, what="intervention_axis")
    if axis not in surrogate.axes:
        raise ValueError(f"{axis!r} is not an input axis of the surrogate")
    settings = tuple(float(value) for value in training_interventions)
    if not settings:
        raise ValueError("training_interventions must not be empty")
    if any(math.isclose(float(held_out_intervention), value) for value in settings):
        raise ValueError("the held-out intervention must not appear in training")

    def _rmse(inputs: Any, truth: Any) -> float:
        matrix = _finite_matrix(inputs, n_axes=len(surrogate.axes), what="inputs")
        expected = _finite_vector(truth, what="truth")
        if expected.size != matrix.shape[0]:
            raise ValueError("inputs and truth must have equal row counts")
        residual = surrogate.model.predict(matrix) - expected
        return float(math.sqrt(float(residual @ residual) / residual.size))

    interpolation_rmse = _rmse(interpolation_inputs, interpolation_truth)
    held_out_rmse = _rmse(held_out_inputs, held_out_truth)
    report = InterventionGeneralizationReport(
        surrogate_id=surrogate.surrogate_id,
        intervention_axis=axis,
        training_interventions=settings,
        held_out_intervention=float(held_out_intervention),
        interpolation_rmse=interpolation_rmse,
        held_out_rmse=held_out_rmse,
        predictive_sigma=surrogate.model.residual_sigma,
        sigma_gate=SIGMA_GATE,
        seed=int(seed),
    )
    if not report.mechanism_falsified:
        return report

    index = surrogate.axes.index(axis)
    probe = list(np.asarray(surrogate.box.lower) + 0.5 * surrogate.box.extent)
    probe[index] = float(held_out_intervention)
    request = NativeQueryRequest(
        request_id=f"{surrogate.surrogate_id}:held-out-intervention:{axis}",
        reason=SurrogateVerdict.REFUSED_INTERVENTION_UNSUPPORTED,
        rationale=(
            f"held-out {axis}={float(held_out_intervention):g} RMSE {held_out_rmse:.4g} exceeds "
            f"{SIGMA_GATE:g}x the surrogate's own sigma {surrogate.model.residual_sigma:.4g}, "
            f"while interpolation RMSE at trained settings is {interpolation_rmse:.4g}; the "
            "surrogate reproduces the training distribution without the mechanism, so this "
            "setting must be evaluated natively"
        ),
        axes=surrogate.axes,
        query_point=tuple(float(value) for value in probe),
        diagnostic={
            "interpolation_rmse": interpolation_rmse,
            "held_out_rmse": held_out_rmse,
            "predictive_sigma": surrogate.model.residual_sigma,
            "degradation_ratio": report.degradation_ratio,
        },
        expansions=_REFUSAL_EXPANSIONS[SurrogateVerdict.REFUSED_INTERVENTION_UNSUPPORTED],
    )
    return InterventionGeneralizationReport(
        surrogate_id=report.surrogate_id,
        intervention_axis=report.intervention_axis,
        training_interventions=report.training_interventions,
        held_out_intervention=report.held_out_intervention,
        interpolation_rmse=report.interpolation_rmse,
        held_out_rmse=report.held_out_rmse,
        predictive_sigma=report.predictive_sigma,
        sigma_gate=report.sigma_gate,
        seed=report.seed,
        native_request=request,
    )


class GateEvidenceLedger:
    """Provenance closure that decides whether an artifact may close a gate.

    :class:`ArtifactEnvelope` already refuses to spell native provenance directly.  The route it
    cannot see is indirect: an artifact that is not itself a surrogate but whose ancestry contains
    one.  This ledger walks ``source_artifact_ids`` transitively and refuses the whole chain,
    naming the laundering path so the refusal is diagnosable.  An unregistered ancestor is also a
    refusal — unverifiable ancestry fails closed.
    """

    __slots__ = ("_envelopes",)

    def __init__(self) -> None:
        self._envelopes: dict[str, ArtifactEnvelope] = {}

    def register(self, envelope: ArtifactEnvelope) -> None:
        """Record one artifact envelope.

        Args:
            envelope: The envelope to register.

        Raises:
            ValueError: If a different envelope is already registered under the same id.
        """
        if not isinstance(envelope, ArtifactEnvelope):
            raise TypeError("envelope must be an ArtifactEnvelope")
        existing = self._envelopes.get(envelope.artifact_id)
        if existing is not None and existing is not envelope:
            raise ValueError(f"artifact id {envelope.artifact_id!r} is already registered")
        self._envelopes[envelope.artifact_id] = envelope

    def surrogate_ancestry(self, artifact_id: str) -> tuple[str, ...]:
        """Return the first artifact-id path from ``artifact_id`` to surrogate provenance.

        Args:
            artifact_id: Registered artifact to inspect.

        Returns:
            The path, root first, or an empty tuple when no ancestor is surrogate provenance.

        Raises:
            GateAuthorityError: If the artifact or an ancestor is unregistered, or the ancestry
                graph contains a cycle.
        """
        name = _nonempty(artifact_id, what="artifact_id")
        stack: list[tuple[str, tuple[str, ...]]] = [(name, (name,))]
        seen: set[str] = set()
        while stack:
            current, path = stack.pop()
            envelope = self._envelopes.get(current)
            if envelope is None:
                raise GateAuthorityError(
                    f"artifact {current!r} is not registered; unverifiable ancestry fails closed"
                )
            if current in path[:-1]:
                raise GateAuthorityError(f"ancestry cycle through {current!r}")
            if (
                envelope.evidence_source is EvidenceSource.SURROGATE
                or envelope.representation.kind is RepresentationKind.SURROGATE
            ):
                return path
            seen.add(current)
            for parent in envelope.source_artifact_ids:
                if parent not in seen:
                    stack.append((parent, (*path, parent)))
        return ()

    def assert_gate_closable(self, artifact_id: str, *, gate_id: str) -> None:
        """Refuse to let this artifact close a gate.

        Args:
            artifact_id: Registered artifact proposed as gate evidence.
            gate_id: The gate it would close, used in the refusal message.

        Raises:
            SurrogateAuthorityError: If the artifact or any ancestor carries surrogate provenance.
            GateAuthorityError: Otherwise — this generic layer is an unverified observer and cannot
                produce ``NATIVE_ACCEPTED`` evidence at all, so nothing it holds closes a gate.
        """
        gate = _nonempty(gate_id, what="gate_id")
        path = self.surrogate_ancestry(artifact_id)
        if path:
            raise SurrogateAuthorityError(
                f"gate {gate!r} cannot be closed by {artifact_id!r}: surrogate provenance reaches "
                "it through " + " -> ".join(path)
            )
        raise GateAuthorityError(
            f"gate {gate!r} cannot be closed by {artifact_id!r}: the generic virtual-cell layer is "
            "an unverified observer and holds no NATIVE_ACCEPTED evidence"
        )
