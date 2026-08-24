"""Direction 1, measured instead of asserted: what it costs to carry a distribution at a node.

**The engine carries points.**  Every node in ``aleph/ac`` is a position ``x``, every force is
evaluated once at that position, and nothing in this module changes that or claims otherwise.  What
this module does is turn the master plan's §14.5.4 direction-1 assertion —

    force must be computed from the moments of the distribution, not once at the mean

— into a set of numbers on systems whose truth is available in closed form, so that a future engine
change is argued from a measurement rather than from a slogan.  Everything here is
:data:`~aleph.virtual_cell.contracts.EvidenceSource.ANALYTIC_ORACLE`; it is ``oracle-only``, it
grants no evidence rung, it imports no Warp, and it has not measured one quantity produced by the
CUDA runtime.

WHAT IS ACTUALLY MEASURED, AND WHY EACH PIECE IS NEEDED.

1. **The gap itself.**  ``E[F(X)] != F(E[X])`` for nonlinear ``F``.  Stated that way it is a
   textbook inequality and proves nothing about this project.  Here every force law is a polynomial
   plus a sum of exponentials, which is exactly the family whose Gaussian expectation is closed form
   (raw-moment recursion for the polynomial part, the MGF for the exponential part), so the measured
   gap is scored against truth and not against a bigger simulation.  The **affine control must
   return exactly ``0.0``** — not "small", exactly zero — because a machine that invents a gap where
   Jensen forbids one would make every other number here worthless.

2. **How the gap scales, and where the usual excuse for ignoring it stops working.**  The standard
   defence of a point evaluation is that the correction is second order: ``E[F] ~ F(mu) +
   F''(mu) sigma^2 / 2``.  That defence is testable in two directions.  The fitted ``log|gap|`` vs
   ``log sigma`` slope says whether what was measured is the moment correction at all — and the
   predicted slope is itself computable, as the smallest EVEN derivative order that does not vanish
   at ``mu``, because odd Gaussian central moments are zero.  A law whose ``F''(mu) = 0`` therefore
   has slope 4, not 2, and :func:`leading_gap_order` says so before the fit is run.  Then
   :func:`find_closure_breakdown_sigma` measures the width at which the second-order truncation
   itself departs from the exact expectation by more than a declared tolerance.  **That number is
   the deliverable**: it is what a future engine would consult to decide whether a Gaussian closure
   is still allowed at the width it is actually carrying.

3. **Propagating ``(mu, Sigma)`` one step, and the case where it must refuse.**  Three updates are
   compared, and the comparison is arranged so that truth is never a bigger approximation:

   * the **exact** moment update, computed from Gaussian moments of the polynomial/exponential drift;
   * the **Taylor closure** a real implementation would use, ``mu + (F(mu) + F''(mu) Sigma / 2) dt``
     with a linearised covariance;
   * a **deterministic Gauss-Hermite ensemble**, which is a genuinely different computation (a
     weighted sum over nodes) and is exact for polynomials up to degree ``2n - 1``.

   A seeded Monte-Carlo ensemble is included ONLY to show what would have happened had truth been a
   sampled ensemble: its error is many orders of magnitude larger than the effect being measured, so
   validating a closure against a sampled ensemble would have hidden the closure error entirely.

4. **The refusal.**  For a double well the true stationary density is bimodal and a single Gaussian
   reports a mean sitting on the barrier, where the density is a measured fraction of the modal
   density.  The plan prohibits exactly this ("bimodal event state must not be hidden in a single
   Gaussian"), so :func:`propagate_local_gaussian` REFUSES, and the refusal is driven by measured
   statistics — mode count with a declared prominence, the population bimodality coefficient, and
   the density ratio at the closure mean — never by a hand-set flag.  A refused verdict has **no
   mean field to read**: the propagated state lives inside the admitted-only
   :class:`GaussianClosureStep`, so a caller cannot ignore the refusal by reaching past it.
   A nonlinear drift with NO shape evidence supplied is also refused, because applying a Gaussian
   closure while knowing nothing about the shape of the true density is the failure mode itself.

5. **Source-factored, not one Sigma** (§2.2, audit A1).  Two contrasts, both closed form, and they
   fail in different ways on purpose:

   * **placement changes the first moment.**  The same total variance ``v`` placed on the position
     or on the rate of ``A exp(b x)`` gives ``A e^{b mu + b^2 v / 2}`` versus
     ``A e^{b mu + mu^2 v / 2}``.  Neither correction is zero and their ratio is
     ``exp((mu^2 - b^2) v / 2)`` — one number, one total variance, two different mean forces.
   * **placement does NOT change the first moment, and the naive argument fails here.**  Splitting
     the same total variance between a quenched offset shared by ``N`` bonds and thermal noise
     independent per bond leaves ``E[mean force]`` invariant to round-off.  What changes is the
     variance of the aggregate, by exactly a factor ``N`` between the all-quenched and all-thermal
     limits.  So "keep the sources separate because the mean moves" is only half true, and this
     module measures the half that is false as well as the half that is true.

units: unit-agnostic and consistently so.  A law carries its force-unit string; ``dt`` is [s]; a
diffusion matrix is [state-unit^2 / s] so that ``2 D dt`` is a covariance; ``sigma`` and ``mu`` share
the state's unit.  No formula mixes two unit systems, and none of the defaults that would let it do
so exist, because there are no defaults.

Sanity Gate (recorded before first execution):

* **dimensional** — ``2 * diffusion * dt`` is the only place a diffusion enters and it enters as a
  covariance increment, checked square and symmetric against the state dimension; a scalar law's
  ``sigma`` is a standard deviation and ``Sigma`` a variance, and the two are never passed to the
  same argument (the constructors take ``covariance``, the scalar-law entry points take ``sigma``,
  and both names appear in every report so a reader can tell which was declared).
* **boundary** — ``sigma = 0`` is legal and must return exactly ``0.0`` gap for EVERY law, since a
  degenerate distribution is a point; a negative ``sigma`` or a non-PSD covariance raises; a
  zero-width bracket, a non-straddling bracket, a bracket whose low end sits at or below its own
  cancellation floor, or a law whose exact gap is identically zero makes
  the breakdown search REFUSE rather than return an endpoint; a density grid that does not capture
  its own tails is caught by :func:`~aleph.virtual_cell.probability.audit_probability_mass` rather
  than silently renormalised.
* **conservation / invariant** — Jensen: the gap is exactly ``0.0`` for an affine law at any width;
  the Gauss-Hermite weights sum to one and are positive; a propagated covariance is symmetric PSD by
  construction (every update is a congruence plus a PSD term); the reference density's mass is
  audited, never clipped or renormalised into agreement.
* **numerical** — Gaussian raw moments come from the recursion ``m_n = mu m_{n-1} + (n-1) var
  m_{n-2}`` rather than from an alternating explicit sum, which cancels catastrophically at large
  ``mu / sigma``; exponential cross moments use the Gaussian tilt identity ``E[X^n e^{bX}] =
  e^{b mu + b^2 var / 2} m_n(mu + b var, var)`` rather than numerical integration; the matrix square
  root is taken from the symmetric eigendecomposition so a singular covariance is legal.  The gap is
  itself a difference of nearby numbers, so :func:`gap_cancellation_floor` reports the width below
  which its apparent truncation error is round-off rather than truncation — measured on the cubic
  control, whose expansion is an algebraic identity, the spurious "truncation error" is 7.1e-05 at
  ``sigma = 1e-06`` and 4.6e-13 at ``sigma = 1e-02``, growing like ``sigma^-2``.  A breakdown search
  whose lower bracket sits in that band is measuring floating point, and since the Codex external
  review of 2026-07-29 it is refused rather than merely annotated: the floor is a GATE on the
  bracket, not a footnote beside the answer.
* **sign-sense** — the gap is reported SIGNED as ``E[F(X)] - F(E[X])``, with both terms carried, so a
  reader never has to guess the direction; a Morse bond at its own minimum has ``F(mu) = 0`` and a
  strictly POSITIVE (repulsive) mean force, and that sign is the physical content of the measurement.
* **measurement-protocol** — every threshold lives in :class:`MomentClosureContract`, **no field of
  which has a default**, and the contract is copied into every report.  Each refusal names the
  representation it routes to.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType
from typing import Any, Protocol

import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov

from aleph.virtual_cell.contracts import (
    EvidenceSource,
    RepresentationKind,
    SandboxExperimentCard,
    UncertaintySource,
)
from aleph.virtual_cell.probability import ProbabilityMassAudit, audit_probability_mass
from aleph.virtual_cell.regime import EVIDENCE_AUTHORITY
from aleph.virtual_cell.sweep import ClosureGapReport, ExpansionRoute, moment_closure_gap

__all__ = [
    "EVIDENCE_AUTHORITY",
    "EVIDENCE_SOURCE",
    "ENGINE_NODE_REPRESENTATION_TODAY",
    "AnalyticForceLaw",
    "BimodalityReport",
    "ClosureBreakdownRefusal",
    "ClosureBreakdownReport",
    "ClosureComparison",
    "ClosureRefusedError",
    "ClosureVerdict",
    "ClosureVerdictKind",
    "Drift",
    "GapScalingReport",
    "GaussianClosureStep",
    "GaussianState",
    "JensenGapMeasurement",
    "LinearDrift",
    "MomentClosureContract",
    "NodeRepresentation",
    "PlacementContrast",
    "QuenchedThermalSplit",
    "ScalarLawDrift",
    "SourceFactoredCovariance",
    "boltzmann_density",
    "compare_closure_against_ensemble",
    "cubic_anharmonic_spring",
    "exponential_repulsion",
    "find_closure_breakdown_sigma",
    "gap_cancellation_floor",
    "gauss_hermite_nodes",
    "leading_gap_order",
    "linear_spring",
    "make_gaussian_state",
    "measure_bimodality",
    "measure_cubic_placement_contrast",
    "measure_gap_sigma_scaling",
    "measure_jensen_gap",
    "measure_quenched_thermal_split",
    "measure_rate_placement_contrast",
    "monte_carlo_moment_update",
    "morse_bond",
    "propagate_gaussian_ou_exact",
    "propagate_local_gaussian",
    "quadrature_moment_update",
    "quartic_force",
    "stochastic_node_experiment_card",
]

EVIDENCE_SOURCE = EvidenceSource.ANALYTIC_ORACLE
"""Everything this module emits is analytic-oracle evidence.  It is oracle-only."""

ENGINE_NODE_REPRESENTATION_TODAY = (
    "the Warp/CUDA runtime carries node = x, a point, and evaluates every force once at that point; "
    "this module measures what a distributional node would cost and claims nothing about the engine"
)
"""Stated so that no reader of a report can mistake this lane for a change to the runtime."""


# ---------------------------------------------------------------------------------------------
# validation helpers
# ---------------------------------------------------------------------------------------------


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _finite(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite; got {value!r}")
    return number


def _positive(value: object, *, what: str) -> float:
    number = _finite(value, what=what)
    if number <= 0.0:
        raise ValueError(f"{what} must be positive; got {value!r}")
    return number


def _nonnegative(value: object, *, what: str) -> float:
    number = _finite(value, what=what)
    if number < 0.0:
        raise ValueError(f"{what} must be nonnegative; got {value!r}")
    return number


def _positive_int(value: object, *, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{what} must be an int")
    if value <= 0:
        raise ValueError(f"{what} must be positive; got {value!r}")
    return value


def _frozen_matrix(values: Any, *, what: str, dimension: int | None = None) -> np.ndarray:
    matrix = np.array(values, dtype=np.float64, copy=True)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{what} must be a square matrix")
    if dimension is not None and matrix.shape[0] != dimension:
        raise ValueError(f"{what} must be {dimension}x{dimension}")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{what} must be finite")
    matrix.flags.writeable = False
    return matrix


def _frozen_vector(values: Any, *, what: str, dimension: int | None = None) -> np.ndarray:
    vector = np.array(values, dtype=np.float64, copy=True).ravel()
    if vector.size == 0:
        raise ValueError(f"{what} must not be empty")
    if dimension is not None and vector.size != dimension:
        raise ValueError(f"{what} must have {dimension} entries")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{what} must be finite")
    vector.flags.writeable = False
    return vector


def _check_symmetric_psd(
    matrix: np.ndarray,
    *,
    what: str,
    symmetry_tolerance: float,
    minimum_eigenvalue_tolerance: float,
) -> None:
    asymmetry = float(np.max(np.abs(matrix - matrix.T))) if matrix.size else 0.0
    if asymmetry > symmetry_tolerance:
        raise ValueError(f"{what} must be symmetric; max asymmetry {asymmetry:.3g}")
    smallest = float(np.min(np.linalg.eigvalsh(0.5 * (matrix + matrix.T))))
    if smallest < -minimum_eigenvalue_tolerance:
        raise ValueError(f"{what} must be positive semidefinite; min eigenvalue {smallest:.3g}")


def _symmetric_sqrt(matrix: np.ndarray) -> np.ndarray:
    """Symmetric square root via eigendecomposition, so a singular covariance stays legal."""
    eigenvalues, vectors = np.linalg.eigh(0.5 * (matrix + matrix.T))
    clipped = np.clip(eigenvalues, 0.0, None)
    return vectors @ np.diag(np.sqrt(clipped)) @ vectors.T


# ---------------------------------------------------------------------------------------------
# contract
# ---------------------------------------------------------------------------------------------


class NodeRepresentation(StrEnum):
    """What a node is being carried as.  The engine today is :attr:`POINT` and only that."""

    POINT = "point"
    LOCAL_GAUSSIAN = "local-gaussian"
    QUADRATURE_ENSEMBLE = "quadrature-ensemble"
    MONTE_CARLO_ENSEMBLE = "monte-carlo-ensemble"
    GAUSSIAN_MIXTURE = "gaussian-mixture"


class ClosureVerdictKind(StrEnum):
    """Whether a local-Gaussian step was issued, and if not, which measurement stopped it."""

    ADMITTED = "admitted"
    REFUSED_BIMODAL = "refused-bimodal"
    REFUSED_UNVERIFIED_SHAPE = "refused-unverified-shape"
    REFUSED_WIDTH_BEYOND_VALIDATED_RANGE = "refused-width-beyond-validated-range"


class ClosureRefusedError(RuntimeError):
    """Raised when a caller asks a refused :class:`ClosureVerdict` for its propagated state."""


@dataclass(frozen=True, slots=True)
class MomentClosureContract:
    """Caller-declared thresholds.  **No field has a default, by design.**

    Every number this module compares against lives here and is copied verbatim into the reports, so
    a verdict can be re-judged without rerunning it.  A threshold chosen inside the module would be
    a threshold chosen after seeing the data.

    Args:
        quadrature_nodes: Gauss-Hermite nodes per axis.  Exact for polynomials of degree
            ``2 * quadrature_nodes - 1``; a degree-``d`` drift needs ``2d`` for its covariance.
        gap_relative_tolerance: Limit on ``|gap| / max(|E[F]|, |F(mu)|)`` handed to
            :func:`~aleph.virtual_cell.sweep.moment_closure_gap`.  Above it, the mean-field
            closure is not adequate for that observable at that width.
        expansion_relative_error_tolerance: Limit on the relative error of the second-order
            truncation against the exact expectation.  This is what defines "the width at which a
            Gaussian closure stops being safe".
        derivative_zero_tolerance: Magnitude below which a derivative at ``mu`` counts as zero when
            :func:`leading_gap_order` looks for the leading even order.  Scale-dependent, hence
            declared rather than assumed.
        max_admissible_modes: Largest mode count in the reference density a single Gaussian may
            still represent.  Anything above it is a mixture, not a Gaussian.
        mode_prominence_ratio: A local maximum counts as a mode only if its prominence exceeds this
            fraction of the peak density.  Without it, grid noise is a mode.
        bimodality_coefficient_max: Upper limit on ``(skew^2 + 1) / kurtosis``.  The uniform
            reference value is ``5/9``; a Gaussian sits at ``1/3``.
        modal_density_ratio_min: Lower limit on ``rho(closure mean) / max(rho)``.  This is the
            statistic that says the reported mean sits where no probability mass is.
        mass_normalization_tolerance: Passed to
            :func:`~aleph.virtual_cell.probability.audit_probability_mass`; it is what catches a
            density grid that does not contain its own tails.
        mass_negativity_tolerance: Passed to the same audit.
        covariance_symmetry_tolerance: Largest tolerated ``max|Sigma - Sigma^T|``.
        minimum_eigenvalue_tolerance: How negative an eigenvalue of ``Sigma`` may be and still count
            as round-off rather than a non-PSD covariance.
    """

    quadrature_nodes: int
    gap_relative_tolerance: float
    expansion_relative_error_tolerance: float
    derivative_zero_tolerance: float
    max_admissible_modes: int
    mode_prominence_ratio: float
    bimodality_coefficient_max: float
    modal_density_ratio_min: float
    mass_normalization_tolerance: float
    mass_negativity_tolerance: float
    covariance_symmetry_tolerance: float
    minimum_eigenvalue_tolerance: float

    def __post_init__(self) -> None:
        if _positive_int(self.quadrature_nodes, what="quadrature_nodes") < 2:
            raise ValueError("quadrature_nodes must be at least 2")
        _positive_int(self.max_admissible_modes, what="max_admissible_modes")
        for name in (
            "gap_relative_tolerance",
            "expansion_relative_error_tolerance",
            "derivative_zero_tolerance",
            "mode_prominence_ratio",
            "bimodality_coefficient_max",
            "modal_density_ratio_min",
            "mass_normalization_tolerance",
            "mass_negativity_tolerance",
            "covariance_symmetry_tolerance",
            "minimum_eigenvalue_tolerance",
        ):
            _nonnegative(getattr(self, name), what=name)
        for name in (
            "mode_prominence_ratio",
            "bimodality_coefficient_max",
            "modal_density_ratio_min",
        ):
            if getattr(self, name) > 1.0:
                raise ValueError(f"{name} is a ratio and must not exceed one")

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-able contract, so a verdict can be re-judged without the raw series."""
        return {
            "quadrature_nodes": self.quadrature_nodes,
            "gap_relative_tolerance": self.gap_relative_tolerance,
            "expansion_relative_error_tolerance": self.expansion_relative_error_tolerance,
            "derivative_zero_tolerance": self.derivative_zero_tolerance,
            "max_admissible_modes": self.max_admissible_modes,
            "mode_prominence_ratio": self.mode_prominence_ratio,
            "bimodality_coefficient_max": self.bimodality_coefficient_max,
            "modal_density_ratio_min": self.modal_density_ratio_min,
            "mass_normalization_tolerance": self.mass_normalization_tolerance,
            "mass_negativity_tolerance": self.mass_negativity_tolerance,
            "covariance_symmetry_tolerance": self.covariance_symmetry_tolerance,
            "minimum_eigenvalue_tolerance": self.minimum_eigenvalue_tolerance,
        }


# ---------------------------------------------------------------------------------------------
# Gaussian moments and the analytic force-law family
# ---------------------------------------------------------------------------------------------


def _raw_moments(mean: float, variance: float, order: int) -> np.ndarray:
    """Raw Gaussian moments ``E[X^n]`` for ``n = 0..order`` by the stable recursion."""
    moments = np.zeros(order + 1, dtype=np.float64)
    moments[0] = 1.0
    if order >= 1:
        moments[1] = mean
    for n in range(2, order + 1):
        moments[n] = mean * moments[n - 1] + (n - 1) * variance * moments[n - 2]
    return moments


@dataclass(frozen=True, slots=True)
class AnalyticForceLaw:
    """A scalar force law whose Gaussian expectation is available in closed form.

    ``F(x) = sum_n c_n x^n + sum_j A_j exp(b_j x)``.  This family is not a convenience: it is
    precisely the set of laws for which ``E[F(X)]``, ``E[X F(X)]`` and ``E[F(X)^2]`` under
    ``X ~ N(mu, sigma^2)`` are exact — raw-moment recursion for the polynomial part, the moment
    generating function for the exponential part, and the Gaussian tilt identity for the cross
    terms.  It covers a linear spring, an anharmonic (cubic) spring, a quartic force, a screened
    exponential repulsion, a Bell-type exponential load response, and a Morse bond, which is every
    shape this lane needs to make a statement about.

    Attributes:
        name: Identifier carried into every report.
        polynomial_coefficients: ``c_n``, index = power.
        exponential_terms: ``(A_j, b_j)`` pairs.
        unit: Force unit string, carried rather than assumed.
        description: What the law is a model of, in one line.
    """

    name: str
    polynomial_coefficients: tuple[float, ...]
    exponential_terms: tuple[tuple[float, float], ...]
    unit: str
    description: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _nonempty(self.name, what="name"))
        object.__setattr__(self, "unit", _nonempty(self.unit, what="unit"))
        object.__setattr__(self, "description", _nonempty(self.description, what="description"))
        coefficients = tuple(
            _finite(value, what="polynomial coefficient") for value in self.polynomial_coefficients
        )
        object.__setattr__(self, "polynomial_coefficients", coefficients)
        terms = tuple(
            (
                _finite(amplitude, what="exponential amplitude"),
                _finite(rate, what="exponential rate"),
            )
            for amplitude, rate in self.exponential_terms
        )
        object.__setattr__(self, "exponential_terms", terms)
        if not coefficients and not terms:
            raise ValueError("a force law must have at least one term")

    @property
    def is_affine(self) -> bool:
        """Whether ``F`` is affine, in which case Jensen forbids any gap at all."""
        if self.exponential_terms:
            return False
        return all(value == 0.0 for value in self.polynomial_coefficients[2:])

    def evaluate(self, x: Any) -> np.ndarray:
        """Evaluate ``F`` elementwise.

        Args:
            x: Positions, any array-like.

        Returns:
            ``F(x)`` with the same shape as ``x``.
        """
        points = np.asarray(x, dtype=np.float64)
        total = np.zeros_like(points, dtype=np.float64)
        for power, coefficient in enumerate(self.polynomial_coefficients):
            if coefficient != 0.0:
                total = total + coefficient * points**power
        for amplitude, rate in self.exponential_terms:
            total = total + amplitude * np.exp(rate * points)
        return total

    def derivative(self, order: int, x: Any) -> np.ndarray:
        """Evaluate the ``order``-th derivative of ``F`` elementwise.

        Args:
            order: Derivative order; ``0`` returns the law itself.
            x: Positions, any array-like.

        Returns:
            ``F^{(order)}(x)``.

        Raises:
            ValueError: If ``order`` is negative.
        """
        if isinstance(order, bool) or not isinstance(order, int) or order < 0:
            raise ValueError("derivative order must be a nonnegative int")
        if order == 0:
            return self.evaluate(x)
        points = np.asarray(x, dtype=np.float64)
        total = np.zeros_like(points, dtype=np.float64)
        for power, coefficient in enumerate(self.polynomial_coefficients):
            if coefficient != 0.0 and power >= order:
                factor = math.factorial(power) / math.factorial(power - order)
                total = total + coefficient * factor * points ** (power - order)
        for amplitude, rate in self.exponential_terms:
            total = total + amplitude * rate**order * np.exp(rate * points)
        return total

    def force_at_mean(self, mean: float) -> float:
        """Return ``F(E[X])`` — what the engine computes today.

        Args:
            mean: The node position.

        Returns:
            The point-evaluated force.
        """
        return float(self.evaluate(np.asarray(_finite(mean, what="mean"), dtype=np.float64)))

    def exact_expectation(self, mean: float, sigma: float) -> float:
        """Return the closed-form ``E[F(X)]`` for ``X ~ N(mean, sigma^2)``.

        Args:
            mean: Distribution mean.
            sigma: Distribution standard deviation; ``0`` is legal and reduces to a point.

        Returns:
            The exact expectation.
        """
        mu = _finite(mean, what="mean")
        variance = _nonnegative(sigma, what="sigma") ** 2
        total = 0.0
        if self.polynomial_coefficients:
            moments = _raw_moments(mu, variance, len(self.polynomial_coefficients) - 1)
            total += float(np.dot(np.asarray(self.polynomial_coefficients), moments))
        for amplitude, rate in self.exponential_terms:
            total += amplitude * math.exp(rate * mu + 0.5 * rate * rate * variance)
        return total

    def exact_cross_moment(self, mean: float, sigma: float) -> float:
        """Return the closed-form ``E[X F(X)]``, needed for the exact covariance update.

        Args:
            mean: Distribution mean.
            sigma: Distribution standard deviation.

        Returns:
            ``E[X F(X)]``.
        """
        mu = _finite(mean, what="mean")
        variance = _nonnegative(sigma, what="sigma") ** 2
        total = 0.0
        if self.polynomial_coefficients:
            moments = _raw_moments(mu, variance, len(self.polynomial_coefficients))
            for power, coefficient in enumerate(self.polynomial_coefficients):
                total += coefficient * moments[power + 1]
        for amplitude, rate in self.exponential_terms:
            tilt = math.exp(rate * mu + 0.5 * rate * rate * variance)
            total += amplitude * tilt * (mu + rate * variance)
        return total

    def exact_second_moment(self, mean: float, sigma: float) -> float:
        """Return the closed-form ``E[F(X)^2]``.

        Args:
            mean: Distribution mean.
            sigma: Distribution standard deviation.

        Returns:
            ``E[F(X)^2]``.
        """
        mu = _finite(mean, what="mean")
        variance = _nonnegative(sigma, what="sigma") ** 2
        coefficients = np.asarray(self.polynomial_coefficients, dtype=np.float64)
        total = 0.0
        if coefficients.size:
            squared = np.convolve(coefficients, coefficients)
            moments = _raw_moments(mu, variance, squared.size - 1)
            total += float(np.dot(squared, moments))
        for amplitude, rate in self.exponential_terms:
            tilt = math.exp(rate * mu + 0.5 * rate * rate * variance)
            if coefficients.size:
                tilted = _raw_moments(mu + rate * variance, variance, coefficients.size - 1)
                total += 2.0 * amplitude * tilt * float(np.dot(coefficients, tilted))
        for amplitude_j, rate_j in self.exponential_terms:
            for amplitude_k, rate_k in self.exponential_terms:
                rate = rate_j + rate_k
                total += (
                    amplitude_j * amplitude_k * math.exp(rate * mu + 0.5 * rate * rate * variance)
                )
        return total

    def exact_gap(self, mean: float, sigma: float) -> float:
        """Return the signed exact Jensen gap ``E[F(X)] - F(E[X])``.

        Args:
            mean: Distribution mean.
            sigma: Distribution standard deviation.

        Returns:
            The gap.  Exactly ``0.0`` for an affine law and for ``sigma == 0``.
        """
        return self.exact_expectation(mean, sigma) - self.force_at_mean(mean)

    def second_order_gap(self, mean: float, sigma: float) -> float:
        """Return the leading moment correction ``F''(mu) sigma^2 / 2``.

        Args:
            mean: Distribution mean.
            sigma: Distribution standard deviation.

        Returns:
            The truncated second-order correction — the term a closure would carry.
        """
        variance = _nonnegative(sigma, what="sigma") ** 2
        second = float(self.derivative(2, np.asarray(_finite(mean, what="mean"))))
        return 0.5 * second * variance

    def potential(self, x: Any) -> np.ndarray:
        """Return ``U`` with ``F = -dU/dx`` and ``U(0) = 0``.

        Args:
            x: Positions, any array-like.

        Returns:
            The potential, used to build the analytic stationary density of a gradient system.

        Raises:
            ValueError: If an exponential term has zero rate, whose antiderivative is not of this
                form.
        """
        points = np.asarray(x, dtype=np.float64)
        total = np.zeros_like(points, dtype=np.float64)
        for power, coefficient in enumerate(self.polynomial_coefficients):
            if coefficient != 0.0:
                total = total - coefficient * points ** (power + 1) / (power + 1)
        for amplitude, rate in self.exponential_terms:
            if rate == 0.0:
                raise ValueError("an exponential term with zero rate has no exponential potential")
            total = total - amplitude * (np.exp(rate * points) - 1.0) / rate
        return total


def linear_spring(*, stiffness: float, rest_position: float, unit: str) -> AnalyticForceLaw:
    """Build ``F(x) = -k (x - x0)`` — the control whose Jensen gap must be exactly zero.

    Args:
        stiffness: ``k``.
        rest_position: ``x0``.
        unit: Force unit string.

    Returns:
        The affine law.
    """
    k = _finite(stiffness, what="stiffness")
    x0 = _finite(rest_position, what="rest_position")
    return AnalyticForceLaw(
        name="linear-spring",
        polynomial_coefficients=(k * x0, -k),
        exponential_terms=(),
        unit=unit,
        description="Hookean spring; affine, so Jensen forbids a gap at any width",
    )


def cubic_anharmonic_spring(
    *, stiffness: float, cubic_coefficient: float, unit: str
) -> AnalyticForceLaw:
    """Build ``F(x) = -(k x + a x^3)`` — the case where the second-order expansion is EXACT.

    A cubic force has no fourth or higher derivative, and the third-order term multiplies a
    vanishing odd central moment, so ``E[F] = F(mu) + F''(mu) sigma^2 / 2`` holds identically.  That
    makes it the right control for the scaling measurement and the wrong one for a breakdown search.

    Args:
        stiffness: ``k``.
        cubic_coefficient: ``a``; positive stiffens the spring.
        unit: Force unit string.

    Returns:
        The cubic law.
    """
    k = _finite(stiffness, what="stiffness")
    a = _finite(cubic_coefficient, what="cubic_coefficient")
    return AnalyticForceLaw(
        name="cubic-anharmonic-spring",
        polynomial_coefficients=(0.0, -k, 0.0, -a),
        exponential_terms=(),
        unit=unit,
        description="anharmonic spring; the second-order moment correction is exact, not truncated",
    )


def quartic_force(*, coefficient: float, unit: str) -> AnalyticForceLaw:
    """Build ``F(x) = c x^4``, whose gap has a genuine fourth-order term.

    At ``mu = 0`` the second derivative vanishes, so the leading gap is ``3 c sigma^4`` and the
    second-order expansion is wrong by 100 % at EVERY width.  That is the adversarial control for
    the scaling fit.

    Args:
        coefficient: ``c``.
        unit: Force unit string.

    Returns:
        The quartic law.
    """
    c = _finite(coefficient, what="coefficient")
    return AnalyticForceLaw(
        name="quartic-force",
        polynomial_coefficients=(0.0, 0.0, 0.0, 0.0, c),
        exponential_terms=(),
        unit=unit,
        description="pure quartic force; leading gap order is 4 at the origin, not 2",
    )


def exponential_repulsion(*, amplitude: float, decay_length: float, unit: str) -> AnalyticForceLaw:
    """Build ``F(x) = A exp(-x / lambda)`` — screened repulsion, Bell-type load response.

    Args:
        amplitude: ``A``.
        decay_length: ``lambda``; must be positive.
        unit: Force unit string.

    Returns:
        The exponential law, whose exact Gaussian expectation is the MGF.
    """
    a = _finite(amplitude, what="amplitude")
    lam = _positive(decay_length, what="decay_length")
    return AnalyticForceLaw(
        name="exponential-repulsion",
        polynomial_coefficients=(),
        exponential_terms=((a, -1.0 / lam),),
        unit=unit,
        description="screened exponential repulsion; exact Gaussian expectation via the MGF",
    )


def morse_bond(
    *, well_depth: float, width: float, equilibrium_position: float, unit: str
) -> AnalyticForceLaw:
    """Build the Morse bond force ``2 a D (e^{-2a(x-x0)} - e^{-a(x-x0)})``.

    A Lennard-Jones-like bond with a hard core, a minimum and a soft tail, and — unlike LJ itself —
    an exact Gaussian expectation, because it is a sum of two exponentials.  At ``mu = x0`` the point
    force is exactly zero while the distributional force is strictly positive: the hard core is
    sampled more than the soft tail, so a spread node is pushed outward by a force a point node
    never feels.

    Args:
        well_depth: ``D``.
        width: ``a``; must be positive.
        equilibrium_position: ``x0``.
        unit: Force unit string.

    Returns:
        The Morse law.
    """
    depth = _finite(well_depth, what="well_depth")
    a = _positive(width, what="width")
    x0 = _finite(equilibrium_position, what="equilibrium_position")
    scale = 2.0 * a * depth
    return AnalyticForceLaw(
        name="morse-bond",
        polynomial_coefficients=(),
        exponential_terms=(
            (scale * math.exp(2.0 * a * x0), -2.0 * a),
            (-scale * math.exp(a * x0), -a),
        ),
        unit=unit,
        description="Morse bond; LJ-like shape with an exact Gaussian expectation",
    )


# ---------------------------------------------------------------------------------------------
# 1. the gap, in closed form and then measured
# ---------------------------------------------------------------------------------------------


def gauss_hermite_nodes(mean: float, sigma: float, *, nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Return deterministic quadrature nodes and weights for ``N(mean, sigma^2)``.

    The weights are positive and sum to one, so the node set is a legitimate (deterministic)
    ensemble and can be handed straight to
    :func:`~aleph.virtual_cell.sweep.moment_closure_gap` as a weighted ensemble.

    Args:
        mean: Distribution mean.
        sigma: Standard deviation; ``0`` collapses every node onto the mean.
        nodes: Node count; exact for polynomials of degree ``2 * nodes - 1``.

    Returns:
        ``(positions, weights)``.
    """
    mu = _finite(mean, what="mean")
    width = _nonnegative(sigma, what="sigma")
    count = _positive_int(nodes, what="nodes")
    abscissas, raw = np.polynomial.hermite.hermgauss(count)
    weights = raw / math.sqrt(math.pi)
    weights = weights / float(np.sum(weights))
    return mu + math.sqrt(2.0) * width * abscissas, weights


@dataclass(frozen=True, slots=True)
class JensenGapMeasurement:
    """One law at one ``(mu, sigma)``: the closed forms, the estimate, and the estimate's error.

    ``exact_gap`` and ``second_order_gap`` are both truth-side quantities; ``quadrature_*`` is the
    estimator being scored against them.  ``closure_report`` is the reused
    :class:`~aleph.virtual_cell.sweep.ClosureGapReport`, so the number this lane quotes and the
    number the sweep lane quotes are literally the same computation.
    """

    law_name: str
    mean: float
    sigma: float
    force_at_mean: float
    exact_expectation: float
    exact_gap: float
    second_order_gap: float
    second_order_relative_error: float
    quadrature_expectation: float
    quadrature_absolute_error: float
    quadrature_relative_error: float
    quadrature_nodes: int
    representation: NodeRepresentation
    closure_report: ClosureGapReport
    unit: str
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able measurement."""
        return {
            "law_name": self.law_name,
            "mean": self.mean,
            "sigma": self.sigma,
            "force_at_mean": self.force_at_mean,
            "exact_expectation": self.exact_expectation,
            "exact_gap": self.exact_gap,
            "second_order_gap": self.second_order_gap,
            "second_order_relative_error": self.second_order_relative_error,
            "quadrature_expectation": self.quadrature_expectation,
            "quadrature_absolute_error": self.quadrature_absolute_error,
            "quadrature_relative_error": self.quadrature_relative_error,
            "quadrature_nodes": self.quadrature_nodes,
            "representation": self.representation.value,
            "closure_report": self.closure_report.as_dict(),
            "unit": self.unit,
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
            "engine_node_representation_today": ENGINE_NODE_REPRESENTATION_TODAY,
        }


def measure_jensen_gap(
    law: AnalyticForceLaw, *, mean: float, sigma: float, contract: MomentClosureContract
) -> JensenGapMeasurement:
    """Measure ``E[F(X)] - F(E[X])`` against the closed form for one law and one width.

    Args:
        law: The analytic force law.
        mean: Distribution mean.
        sigma: Standard deviation; ``0`` is legal and must give exactly zero gap.
        contract: Declared thresholds; supplies the quadrature node count and gap tolerance.

    Returns:
        The :class:`JensenGapMeasurement`, including the estimator's error against the closed form.

    Raises:
        TypeError: If ``law`` or ``contract`` has the wrong type.
    """
    if not isinstance(law, AnalyticForceLaw):
        raise TypeError("law must be an AnalyticForceLaw")
    if not isinstance(contract, MomentClosureContract):
        raise TypeError("contract must be a MomentClosureContract")
    mu = _finite(mean, what="mean")
    width = _nonnegative(sigma, what="sigma")

    exact_expectation = law.exact_expectation(mu, width)
    at_mean = law.force_at_mean(mu)
    exact_gap = exact_expectation - at_mean
    second_order = law.second_order_gap(mu, width)
    if exact_gap == 0.0:
        second_order_relative = 0.0 if second_order == 0.0 else math.inf
    else:
        second_order_relative = abs(exact_gap - second_order) / abs(exact_gap)

    positions, weights = gauss_hermite_nodes(mu, width, nodes=contract.quadrature_nodes)
    report = moment_closure_gap(
        positions,
        law.evaluate,
        observable=f"{law.name}:F(x)",
        tolerance=contract.gap_relative_tolerance,
        weights=weights,
    )
    quadrature_expectation = report.ensemble_mean_of_f
    absolute_error = abs(quadrature_expectation - exact_expectation)
    scale = abs(exact_expectation)
    relative_error = absolute_error / scale if scale > 0.0 else absolute_error

    return JensenGapMeasurement(
        law_name=law.name,
        mean=mu,
        sigma=width,
        force_at_mean=at_mean,
        exact_expectation=exact_expectation,
        exact_gap=exact_gap,
        second_order_gap=second_order,
        second_order_relative_error=second_order_relative,
        quadrature_expectation=quadrature_expectation,
        quadrature_absolute_error=absolute_error,
        quadrature_relative_error=relative_error,
        quadrature_nodes=contract.quadrature_nodes,
        representation=NodeRepresentation.QUADRATURE_ENSEMBLE,
        closure_report=report,
        unit=law.unit,
    )


# ---------------------------------------------------------------------------------------------
# 2. how the gap scales, and where the expansion breaks down
# ---------------------------------------------------------------------------------------------


def leading_gap_order(
    law: AnalyticForceLaw, *, mean: float, max_order: int, zero_tolerance: float
) -> int | None:
    """Return the predicted ``sigma`` exponent of the gap: the leading nonvanishing EVEN order.

    ``gap = sum_{k>=2} F^{(k)}(mu) E[(X - mu)^k] / k!`` and every odd central moment of a Gaussian is
    zero, so only even ``k`` survive and the smallest such ``k`` with ``F^{(k)}(mu) != 0`` fixes the
    slope.  This is the analytic prediction the fitted exponent is scored against; it is computed
    from the law, never read off the fit.

    Args:
        law: The analytic force law.
        mean: Where the derivatives are taken.
        max_order: Highest even order examined.
        zero_tolerance: Magnitude below which a derivative counts as zero.  Caller-declared because
            it is scale-dependent.

    Returns:
        The leading even order, or ``None`` if every examined even derivative vanishes (an affine
        law, where there is no gap to have an order).

    Raises:
        ValueError: If ``max_order`` is below 2.
    """
    if not isinstance(law, AnalyticForceLaw):
        raise TypeError("law must be an AnalyticForceLaw")
    limit = _positive_int(max_order, what="max_order")
    if limit < 2:
        raise ValueError("max_order must be at least 2")
    tolerance = _nonnegative(zero_tolerance, what="zero_tolerance")
    point = np.asarray(_finite(mean, what="mean"), dtype=np.float64)
    for order in range(2, limit + 1, 2):
        if abs(float(law.derivative(order, point))) > tolerance:
            return order
    return None


@dataclass(frozen=True, slots=True)
class GapScalingReport:
    """The fitted ``log|gap|`` vs ``log sigma`` slope, against the order predicted from the law."""

    law_name: str
    mean: float
    sigmas: tuple[float, ...]
    gaps: tuple[float, ...]
    fitted_exponent: float
    predicted_exponent: int | None
    intercept: float
    max_log_residual: float
    r_squared: float
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able scaling report."""
        return {
            "law_name": self.law_name,
            "mean": self.mean,
            "sigmas": list(self.sigmas),
            "gaps": list(self.gaps),
            "fitted_exponent": self.fitted_exponent,
            "predicted_exponent": self.predicted_exponent,
            "intercept": self.intercept,
            "max_log_residual": self.max_log_residual,
            "r_squared": self.r_squared,
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
        }


def measure_gap_sigma_scaling(
    law: AnalyticForceLaw,
    *,
    mean: float,
    sigmas: Sequence[float],
    contract: MomentClosureContract,
    max_order: int,
) -> GapScalingReport:
    """Fit the exponent of ``|gap| ~ sigma^p`` and compare it with the analytic prediction.

    A slope near 2 is what says the measured quantity is the leading moment correction rather than a
    numerical artefact.  A slope near 4 is what an ``F''(mu) = 0`` law must give, and measuring that
    is how the fit demonstrates it can tell orders apart at all.

    Args:
        law: The analytic force law.
        mean: Distribution mean, held fixed across the sweep.
        sigmas: Strictly increasing positive widths, at least three.
        contract: Declared thresholds; supplies the derivative zero tolerance.
        max_order: Highest even derivative order examined for the prediction.

    Returns:
        The :class:`GapScalingReport`.

    Raises:
        ValueError: If fewer than three widths are given, if they are not strictly increasing and
            positive, or if any gap is exactly zero — an affine law has no exponent to fit and the
            refusal is louder than a fitted ``nan``.
    """
    if not isinstance(law, AnalyticForceLaw):
        raise TypeError("law must be an AnalyticForceLaw")
    if not isinstance(contract, MomentClosureContract):
        raise TypeError("contract must be a MomentClosureContract")
    widths = np.asarray([_positive(value, what="sigma") for value in sigmas], dtype=np.float64)
    if widths.size < 3:
        raise ValueError("a scaling fit needs at least three widths")
    if not np.all(np.diff(widths) > 0.0):
        raise ValueError("sigmas must be strictly increasing")
    mu = _finite(mean, what="mean")

    gaps = np.asarray([law.exact_gap(mu, float(width)) for width in widths], dtype=np.float64)
    if np.any(gaps == 0.0):
        raise ValueError(
            f"law {law.name!r} has an identically zero gap at mu={mu!r}; there is no exponent to "
            "fit, and reporting one would invent a scaling that does not exist"
        )

    log_sigma = np.log(widths)
    log_gap = np.log(np.abs(gaps))
    design = np.vstack([log_sigma, np.ones_like(log_sigma)]).T
    (slope, intercept), *_ = np.linalg.lstsq(design, log_gap, rcond=None)
    fitted = design @ np.array([slope, intercept])
    residuals = log_gap - fitted
    total = log_gap - float(np.mean(log_gap))
    denominator = float(np.dot(total, total))
    r_squared = (
        1.0 - float(np.dot(residuals, residuals)) / denominator if denominator > 0.0 else 1.0
    )

    return GapScalingReport(
        law_name=law.name,
        mean=mu,
        sigmas=tuple(float(value) for value in widths),
        gaps=tuple(float(value) for value in gaps),
        fitted_exponent=float(slope),
        predicted_exponent=leading_gap_order(
            law,
            mean=mu,
            max_order=max_order,
            zero_tolerance=contract.derivative_zero_tolerance,
        ),
        intercept=float(intercept),
        max_log_residual=float(np.max(np.abs(residuals))),
        r_squared=r_squared,
    )


class ClosureBreakdownRefusal(StrEnum):
    """Why a breakdown search issued no width.  ``None`` on the report means a width WAS issued.

    Typed because a caller must be able to tell "there is no breakdown width for this law" from
    "this bracket cannot measure one" without parsing prose:
    :attr:`EXPANSION_IS_EXACT` is a property of the law, :attr:`UNSAFE_ACROSS_RANGE` and
    :attr:`BREAKDOWN_ABOVE_RANGE` say the answer lies outside the declared bracket, and
    :attr:`BELOW_CANCELLATION_FLOOR` says the bracket produced no measurement at all.
    """

    EXPANSION_IS_EXACT = "expansion-is-exact"
    UNSAFE_ACROSS_RANGE = "unsafe-across-range"
    BREAKDOWN_ABOVE_RANGE = "breakdown-above-range"
    BELOW_CANCELLATION_FLOOR = "below-cancellation-floor"


@dataclass(frozen=True, slots=True)
class ClosureBreakdownReport:
    """The width at which the second-order expansion stops describing the exact gap.

    ``breakdown_sigma`` is the quotable number: below it a Gaussian closure reproduces the exact
    expectation to within ``tolerance``; above it, it does not, and a future engine carrying a width
    larger than this has no licence to use one.  ``breakdown_sigma is None`` with
    ``expansion_is_exact`` set is the honest report for a cubic law, where the expansion never
    breaks down because it is not an approximation at all.

    ``refusal_kind`` is the typed form of ``refusal_reason`` and is ``None`` exactly when a width was
    issued.  It exists because the fourth refusal — :attr:`ClosureBreakdownRefusal.BELOW_CANCELLATION_FLOOR`
    — is not a statement about the law or about where the answer lies but about whether the bracket
    measured anything, and a caller that cannot distinguish those three cases will read the third as
    the second.
    """

    law_name: str
    mean: float
    tolerance: float
    breakdown_sigma: float | None
    relative_error_at_breakdown: float | None
    search_low: float
    search_high: float
    relative_error_at_low: float
    relative_error_at_high: float
    cancellation_floor_at_low: float
    cancellation_floor_at_high: float
    bisection_iterations: int
    expansion_is_exact: bool
    refusal_kind: ClosureBreakdownRefusal | None
    refusal_reason: str | None
    expansions: tuple[ExpansionRoute, ...]
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    def __post_init__(self) -> None:
        if (self.refusal_kind is None) != (self.breakdown_sigma is not None):
            raise ValueError(
                "a breakdown report carries a width or a typed refusal kind, never both and never "
                "neither"
            )
        if (self.refusal_kind is None) != (self.refusal_reason is None):
            raise ValueError("refusal_kind and refusal_reason must be present or absent together")

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able breakdown report."""
        return {
            "law_name": self.law_name,
            "mean": self.mean,
            "tolerance": self.tolerance,
            "breakdown_sigma": self.breakdown_sigma,
            "relative_error_at_breakdown": self.relative_error_at_breakdown,
            "search_low": self.search_low,
            "search_high": self.search_high,
            "relative_error_at_low": self.relative_error_at_low,
            "relative_error_at_high": self.relative_error_at_high,
            "cancellation_floor_at_low": self.cancellation_floor_at_low,
            "cancellation_floor_at_high": self.cancellation_floor_at_high,
            "bisection_iterations": self.bisection_iterations,
            "expansion_is_exact": self.expansion_is_exact,
            "refusal_kind": None if self.refusal_kind is None else self.refusal_kind.value,
            "refusal_reason": self.refusal_reason,
            "expansions": [route.value for route in self.expansions],
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
        }


def gap_cancellation_floor(law: AnalyticForceLaw, *, mean: float, sigma: float) -> float:
    """Estimate the relative-error floor imposed by cancellation in ``E[F] - F(mu)``.

    The gap is a difference of two nearby numbers, so its own relative accuracy degrades as the
    width shrinks: roughly ``eps (|E[F]| + |F(mu)|) / |gap|``, which grows like ``sigma^-2`` once the
    gap is second order.  **Measured on the cubic control at ``mu = 0.8``: the apparent truncation
    error is 7.1e-05 at ``sigma = 1e-06`` and 4.6e-13 at ``sigma = 1e-02``, for a law whose
    second-order expansion is an algebraic identity** — every one of those digits is cancellation,
    not truncation.  A breakdown search whose lower bracket sits below this floor is measuring
    floating point, not the closure, which is why :func:`find_closure_breakdown_sigma` REFUSES on
    that bracket rather than issuing a width.  Until the Codex external review of 2026-07-29 the
    floor was only reported alongside the bracket and read by nothing, which let a width be issued
    from a bracket whose measured error was nine times below its own floor.

    Args:
        law: The analytic force law.
        mean: Distribution mean.
        sigma: Standard deviation.

    Returns:
        The estimated relative floor; ``inf`` when the gap is exactly zero.
    """
    exact_expectation = law.exact_expectation(mean, sigma)
    at_mean = law.force_at_mean(mean)
    gap = exact_expectation - at_mean
    if gap == 0.0:
        return math.inf
    eps = float(np.finfo(np.float64).eps)
    return eps * (abs(exact_expectation) + abs(at_mean)) / abs(gap)


def _expansion_relative_error(law: AnalyticForceLaw, mean: float, sigma: float) -> float:
    exact = law.exact_gap(mean, sigma)
    truncated = law.second_order_gap(mean, sigma)
    if exact == 0.0:
        return 0.0 if truncated == 0.0 else math.inf
    return abs(exact - truncated) / abs(exact)


def find_closure_breakdown_sigma(
    law: AnalyticForceLaw,
    *,
    mean: float,
    tolerance: float,
    search_low: float,
    search_high: float,
    bisection_iterations: int,
    exactness_tolerance: float,
) -> ClosureBreakdownReport:
    """Find the smallest width at which the second-order closure error exceeds ``tolerance``.

    The search is a bisection on **exact functions** — the closed-form gap and the closed-form
    truncation — so the answer carries no sampling error.  The bracket must straddle the tolerance;
    a bracket that does not is refused rather than being widened automatically, because an
    automatically widened bracket is a threshold chosen after seeing the data.

    "Exact functions" is not the same as exactly evaluated ones.  The gap is a difference of two
    nearby numbers, so before any bisection the bracket is checked against its own
    :func:`gap_cancellation_floor`: a lower bracket whose floor is not below ``tolerance``, or whose
    measured error does not exceed that floor, is refused with
    :attr:`ClosureBreakdownRefusal.BELOW_CANCELLATION_FLOOR` and issues no width.  This gate was
    added after the Codex external review of 2026-07-29 found the floor being computed, recorded and
    never read — see the comment at the gate for the counterexample and the numbers it produced.

    Args:
        law: The analytic force law.
        mean: Distribution mean.
        tolerance: Declared limit on ``|exact gap - second-order gap| / |exact gap|``.
        search_low: Lower bracket width; its relative error must be at or below ``tolerance``.
        search_high: Upper bracket width; its relative error must exceed ``tolerance``.
        bisection_iterations: Number of halvings; the returned width is accurate to
            ``(high - low) / 2^iterations``.
        exactness_tolerance: Relative error at or below which the expansion counts as EXACT rather
            than merely accurate.  Judged at the UPPER bracket, where a real truncation error is
            largest and the cancellation floor of :func:`gap_cancellation_floor` is smallest; judged
            at the lower bracket it would be contaminated by round-off and an exact expansion would
            look inexact.  Declared rather than assumed because "exact" here means "equal up to the
            round-off of two different evaluation orders".

    Returns:
        The :class:`ClosureBreakdownReport`.  A law whose expansion is exact returns
        ``breakdown_sigma=None`` with ``expansion_is_exact=True``; every other refusal returns
        ``breakdown_sigma=None`` with a typed :class:`ClosureBreakdownRefusal`.

    Raises:
        ValueError: If the bracket is empty or inverted.
    """
    if not isinstance(law, AnalyticForceLaw):
        raise TypeError("law must be an AnalyticForceLaw")
    mu = _finite(mean, what="mean")
    limit = _nonnegative(tolerance, what="tolerance")
    low = _positive(search_low, what="search_low")
    high = _positive(search_high, what="search_high")
    if high <= low:
        raise ValueError("search_high must exceed search_low")
    iterations = _positive_int(bisection_iterations, what="bisection_iterations")
    exactness = _nonnegative(exactness_tolerance, what="exactness_tolerance")

    error_low = _expansion_relative_error(law, mu, low)
    error_high = _expansion_relative_error(law, mu, high)
    floor_low = gap_cancellation_floor(law, mean=mu, sigma=low)
    floor_high = gap_cancellation_floor(law, mean=mu, sigma=high)

    if error_high <= exactness:
        return ClosureBreakdownReport(
            law_name=law.name,
            mean=mu,
            tolerance=limit,
            breakdown_sigma=None,
            relative_error_at_breakdown=None,
            search_low=low,
            search_high=high,
            relative_error_at_low=error_low,
            relative_error_at_high=error_high,
            cancellation_floor_at_low=floor_low,
            cancellation_floor_at_high=floor_high,
            bisection_iterations=0,
            expansion_is_exact=True,
            refusal_kind=ClosureBreakdownRefusal.EXPANSION_IS_EXACT,
            refusal_reason=(
                "the second-order expansion is exact for this law at this mean to within the "
                "declared exactness tolerance at the upper bracket, so there is no width at which "
                "it breaks down"
            ),
            expansions=(),
        )

    # --- the cancellation-floor gate --------------------------------------------------------
    # Codex external review 2026-07-29, finding 5: floor_low and floor_high were computed, recorded
    # in the report, and read by NO branch and no bisection condition. The reviewer's counterexample
    # -- F(x) = 1e6 + exp(-x) at mu = 0, tolerance 0.01, bracket [1e-4, 2.0] -- returned
    # breakdown_sigma = 0.2003352947748983 while cancellation_floor_at_low was 0.08871397861214594,
    # nine times the tolerance being scored, and relative_error_at_low was 0.0011703962790697465,
    # itself below that floor. The constant contributes nothing to the gap, so the truncation error
    # is the pure exponential's; what the constant destroys is the CONDITIONING of E[F] - F(mu), and
    # every digit of that 0.00117 was round-off. Reading it as "within 1 % at the low end" is how a
    # `sigma <= 0.2003 lambda` safety range gets issued with no measurement behind it.
    #
    # Only a bracket whose low end is certified "within tolerance" can produce a width, so both
    # halves of the gate are decidable there, before any bisection:
    #   (a) floor_low >= limit  -- the tolerance itself is finer than what the arithmetic resolves;
    #   (b) error_low <= floor_low -- the measured error is at or under its own round-off floor.
    # Neither is widened automatically: the bracket must be re-declared, which is the same rule the
    # non-straddling refusals below already follow.
    floor_is_coarser_than_tolerance = floor_low >= limit
    measurement_is_round_off = error_low <= floor_low
    if floor_is_coarser_than_tolerance or measurement_is_round_off:
        failed = []
        if floor_is_coarser_than_tolerance:
            failed.append(
                f"the cancellation floor at the lower bracket ({floor_low:.6g}) is not below the "
                f"declared tolerance ({limit:.6g})"
            )
        if measurement_is_round_off:
            failed.append(
                f"the measured relative expansion error at the lower bracket ({error_low:.6g}) "
                f"does not exceed the cancellation floor there ({floor_low:.6g})"
            )
        return ClosureBreakdownReport(
            law_name=law.name,
            mean=mu,
            tolerance=limit,
            breakdown_sigma=None,
            relative_error_at_breakdown=None,
            search_low=low,
            search_high=high,
            relative_error_at_low=error_low,
            relative_error_at_high=error_high,
            cancellation_floor_at_low=floor_low,
            cancellation_floor_at_high=floor_high,
            bisection_iterations=0,
            expansion_is_exact=False,
            refusal_kind=ClosureBreakdownRefusal.BELOW_CANCELLATION_FLOOR,
            refusal_reason=(
                f"the bracket [{low:.6g}, {high:.6g}] cannot measure a breakdown width: "
                + "; and ".join(failed)
                + ". A relative error below the cancellation floor of E[F] - F(mu) is floating "
                "point, not closure accuracy, so no width may be issued from this bracket; "
                "re-declare it with a lower end where the gap is resolvable"
            ),
            expansions=(ExpansionRoute.ENSEMBLE,),
        )

    if error_low > limit:
        return ClosureBreakdownReport(
            law_name=law.name,
            mean=mu,
            tolerance=limit,
            breakdown_sigma=None,
            relative_error_at_breakdown=None,
            search_low=low,
            search_high=high,
            relative_error_at_low=error_low,
            relative_error_at_high=error_high,
            cancellation_floor_at_low=floor_low,
            cancellation_floor_at_high=floor_high,
            bisection_iterations=0,
            expansion_is_exact=False,
            refusal_kind=ClosureBreakdownRefusal.UNSAFE_ACROSS_RANGE,
            refusal_reason=(
                "the expansion already exceeds the tolerance at the lower bracket; the closure is "
                "unsafe across the whole search range and the bracket must be re-declared, not "
                "widened after the fact"
            ),
            expansions=(ExpansionRoute.ENSEMBLE,),
        )
    if error_high <= limit:
        return ClosureBreakdownReport(
            law_name=law.name,
            mean=mu,
            tolerance=limit,
            breakdown_sigma=None,
            relative_error_at_breakdown=None,
            search_low=low,
            search_high=high,
            relative_error_at_low=error_low,
            relative_error_at_high=error_high,
            cancellation_floor_at_low=floor_low,
            cancellation_floor_at_high=floor_high,
            bisection_iterations=0,
            expansion_is_exact=False,
            refusal_kind=ClosureBreakdownRefusal.BREAKDOWN_ABOVE_RANGE,
            refusal_reason=(
                "the expansion is still within tolerance at the upper bracket; the breakdown width "
                "is outside the declared search range"
            ),
            expansions=(),
        )

    lower, upper = low, high
    for _ in range(iterations):
        middle = 0.5 * (lower + upper)
        if _expansion_relative_error(law, mu, middle) > limit:
            upper = middle
        else:
            lower = middle
    breakdown = 0.5 * (lower + upper)
    return ClosureBreakdownReport(
        law_name=law.name,
        mean=mu,
        tolerance=limit,
        breakdown_sigma=breakdown,
        relative_error_at_breakdown=_expansion_relative_error(law, mu, breakdown),
        search_low=low,
        search_high=high,
        relative_error_at_low=error_low,
        relative_error_at_high=error_high,
        cancellation_floor_at_low=floor_low,
        cancellation_floor_at_high=floor_high,
        bisection_iterations=iterations,
        expansion_is_exact=False,
        refusal_kind=None,
        refusal_reason=None,
        expansions=(ExpansionRoute.ENSEMBLE, ExpansionRoute.LATENT_COORDINATE),
    )


# ---------------------------------------------------------------------------------------------
# 3. propagating (mu, Sigma) one step
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class GaussianState:
    """A local Gaussian node: mean vector and covariance matrix, validated symmetric PSD."""

    mean: np.ndarray
    covariance: np.ndarray
    symmetry_tolerance: float
    minimum_eigenvalue_tolerance: float

    def __post_init__(self) -> None:
        mean = _frozen_vector(self.mean, what="mean")
        covariance = _frozen_matrix(self.covariance, what="covariance", dimension=mean.size)
        _check_symmetric_psd(
            covariance,
            what="covariance",
            symmetry_tolerance=_nonnegative(self.symmetry_tolerance, what="symmetry_tolerance"),
            minimum_eigenvalue_tolerance=_nonnegative(
                self.minimum_eigenvalue_tolerance, what="minimum_eigenvalue_tolerance"
            ),
        )
        object.__setattr__(self, "mean", mean)
        object.__setattr__(self, "covariance", covariance)

    @property
    def dimension(self) -> int:
        """State dimension."""
        return int(self.mean.size)

    @property
    def scalar_sigma(self) -> float:
        """Standard deviation of a one-dimensional state.

        Raises:
            ValueError: If the state is not one-dimensional; a multivariate covariance has no single
                width and silently taking ``sqrt(trace)`` would be exactly the kind of pre-merge this
                module exists to refuse.
        """
        if self.dimension != 1:
            raise ValueError("scalar_sigma is only defined for a one-dimensional state")
        return math.sqrt(float(self.covariance[0, 0]))

    def matrix_sqrt(self) -> np.ndarray:
        """Return the symmetric square root of the covariance."""
        return _symmetric_sqrt(self.covariance)

    def with_moments(self, mean: Any, covariance: Any) -> GaussianState:
        """Return a new state with the same tolerances and the given moments.

        Args:
            mean: New mean vector.
            covariance: New covariance matrix.

        Returns:
            The new :class:`GaussianState`.
        """
        return GaussianState(
            mean=mean,
            covariance=covariance,
            symmetry_tolerance=self.symmetry_tolerance,
            minimum_eigenvalue_tolerance=self.minimum_eigenvalue_tolerance,
        )


def make_gaussian_state(
    *, mean: Any, covariance: Any, contract: MomentClosureContract
) -> GaussianState:
    """Build a :class:`GaussianState` using the contract's tolerances.

    Args:
        mean: Mean vector.
        covariance: Covariance matrix.
        contract: Declared thresholds.

    Returns:
        The validated state.
    """
    if not isinstance(contract, MomentClosureContract):
        raise TypeError("contract must be a MomentClosureContract")
    return GaussianState(
        mean=mean,
        covariance=covariance,
        symmetry_tolerance=contract.covariance_symmetry_tolerance,
        minimum_eigenvalue_tolerance=contract.minimum_eigenvalue_tolerance,
    )


class Drift(Protocol):
    """A drift whose one-step Gaussian moment update is available exactly."""

    name: str

    @property
    def dimension(self) -> int:
        """State dimension the drift acts on."""

    @property
    def is_gaussian_closed(self) -> bool:
        """Whether the Gaussian family is closed under this drift (true only for affine drift)."""

    def evaluate(self, points: np.ndarray) -> np.ndarray:
        """Evaluate the drift at ``(m, n)`` points, returning ``(m, n)``."""

    def exact_moment_update(
        self, state: GaussianState, *, dt: float, diffusion: np.ndarray
    ) -> GaussianState:
        """Return the exact one-step Euler-Maruyama moments from a Gaussian input."""

    def taylor_closure_update(
        self, state: GaussianState, *, dt: float, diffusion: np.ndarray
    ) -> GaussianState:
        """Return the second-order Taylor closure a real implementation would carry."""

    def mean_truncation_error(self, state: GaussianState) -> float:
        """Return the relative error of the closure's mean drift against the exact expectation."""


def _validated_dt_and_diffusion(
    dt: float, diffusion: Any, *, dimension: int, contract: MomentClosureContract
) -> tuple[float, np.ndarray]:
    step = _positive(dt, what="dt")
    matrix = _frozen_matrix(diffusion, what="diffusion", dimension=dimension)
    _check_symmetric_psd(
        matrix,
        what="diffusion",
        symmetry_tolerance=contract.covariance_symmetry_tolerance,
        minimum_eigenvalue_tolerance=contract.minimum_eigenvalue_tolerance,
    )
    return step, matrix


@dataclass(frozen=True, slots=True)
class LinearDrift:
    """``F(x) = A x + b``.  The one case where a Gaussian closure is not an approximation."""

    matrix: np.ndarray
    offset: np.ndarray
    name: str = "linear-drift"

    def __post_init__(self) -> None:
        matrix = _frozen_matrix(self.matrix, what="drift matrix")
        offset = _frozen_vector(self.offset, what="drift offset", dimension=matrix.shape[0])
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "offset", offset)
        object.__setattr__(self, "name", _nonempty(self.name, what="name"))

    @property
    def dimension(self) -> int:
        """State dimension."""
        return int(self.matrix.shape[0])

    @property
    def is_gaussian_closed(self) -> bool:
        """Affine drift keeps a Gaussian Gaussian, so the closure is exact."""
        return True

    def evaluate(self, points: np.ndarray) -> np.ndarray:
        """Evaluate the drift at ``(m, n)`` points.

        Args:
            points: Positions.

        Returns:
            ``A x + b`` per row.
        """
        array = np.atleast_2d(np.asarray(points, dtype=np.float64))
        return array @ self.matrix.T + self.offset

    def exact_moment_update(
        self, state: GaussianState, *, dt: float, diffusion: np.ndarray
    ) -> GaussianState:
        """Exact Euler-Maruyama moments, assembled term by term.

        Args:
            state: Input Gaussian.
            dt: Step in seconds.
            diffusion: Diffusion matrix; noise covariance is ``2 D dt``.

        Returns:
            The propagated state.
        """
        mean = state.mean + (self.matrix @ state.mean + self.offset) * dt
        covariance = (
            state.covariance
            + dt * (state.covariance @ self.matrix.T + self.matrix @ state.covariance)
            + dt * dt * (self.matrix @ state.covariance @ self.matrix.T)
            + 2.0 * dt * diffusion
        )
        return state.with_moments(mean, 0.5 * (covariance + covariance.T))

    def taylor_closure_update(
        self, state: GaussianState, *, dt: float, diffusion: np.ndarray
    ) -> GaussianState:
        """The closure form, assembled as a congruence — algebraically the same, computed differently.

        Args:
            state: Input Gaussian.
            dt: Step in seconds.
            diffusion: Diffusion matrix.

        Returns:
            The propagated state.
        """
        jacobian = np.eye(self.dimension) + self.matrix * dt
        mean = state.mean + (self.matrix @ state.mean + self.offset) * dt
        covariance = jacobian @ state.covariance @ jacobian.T + 2.0 * dt * diffusion
        return state.with_moments(mean, 0.5 * (covariance + covariance.T))

    def mean_truncation_error(self, state: GaussianState) -> float:
        """Zero by construction: an affine drift has no second derivative to truncate.

        Args:
            state: Input Gaussian.

        Returns:
            ``0.0``.
        """
        if state.dimension != self.dimension:
            raise ValueError("state dimension does not match the drift")
        return 0.0


@dataclass(frozen=True, slots=True)
class ScalarLawDrift:
    """A one-dimensional drift given by an :class:`AnalyticForceLaw` (mobility folded in)."""

    law: AnalyticForceLaw

    def __post_init__(self) -> None:
        if not isinstance(self.law, AnalyticForceLaw):
            raise TypeError("law must be an AnalyticForceLaw")

    @property
    def name(self) -> str:
        """Drift name, taken from the law."""
        return self.law.name

    @property
    def dimension(self) -> int:
        """Always one."""
        return 1

    @property
    def is_gaussian_closed(self) -> bool:
        """Only an affine law keeps a Gaussian Gaussian."""
        return self.law.is_affine

    def evaluate(self, points: np.ndarray) -> np.ndarray:
        """Evaluate the drift at ``(m, 1)`` points.

        Args:
            points: Positions.

        Returns:
            ``F(x)`` per row.
        """
        array = np.atleast_2d(np.asarray(points, dtype=np.float64))
        if array.shape[1] != 1:
            raise ValueError("a scalar law drift acts on one-dimensional points")
        return self.law.evaluate(array)

    def exact_moment_update(
        self, state: GaussianState, *, dt: float, diffusion: np.ndarray
    ) -> GaussianState:
        """Exact Euler-Maruyama moments from the law's closed-form Gaussian moments.

        ``Var[X + F(X) dt] = Var[X] + 2 dt Cov(X, F) + dt^2 Var(F)`` with every term exact, so the
        result is truth for this step and not a finer approximation.

        Args:
            state: Input Gaussian (one-dimensional).
            dt: Step in seconds.
            diffusion: ``1x1`` diffusion matrix.

        Returns:
            The propagated state.
        """
        mu = float(state.mean[0])
        sigma = state.scalar_sigma
        variance = sigma * sigma
        expectation = self.law.exact_expectation(mu, sigma)
        cross = self.law.exact_cross_moment(mu, sigma) - mu * expectation
        force_variance = self.law.exact_second_moment(mu, sigma) - expectation * expectation
        mean = np.array([mu + expectation * dt], dtype=np.float64)
        covariance = np.array(
            [[variance + 2.0 * dt * cross + dt * dt * force_variance + 2.0 * dt * diffusion[0, 0]]],
            dtype=np.float64,
        )
        return state.with_moments(mean, covariance)

    def taylor_closure_update(
        self, state: GaussianState, *, dt: float, diffusion: np.ndarray
    ) -> GaussianState:
        """The second-order closure: ``mu + (F(mu) + F''(mu) Sigma / 2) dt`` with linearised Sigma.

        Args:
            state: Input Gaussian (one-dimensional).
            dt: Step in seconds.
            diffusion: ``1x1`` diffusion matrix.

        Returns:
            The propagated state.
        """
        mu = float(state.mean[0])
        variance = float(state.covariance[0, 0])
        point = np.asarray(mu, dtype=np.float64)
        drift = (
            float(self.law.evaluate(point)) + 0.5 * float(self.law.derivative(2, point)) * variance
        )
        jacobian = 1.0 + float(self.law.derivative(1, point)) * dt
        mean = np.array([mu + drift * dt], dtype=np.float64)
        covariance = np.array(
            [[jacobian * jacobian * variance + 2.0 * dt * diffusion[0, 0]]], dtype=np.float64
        )
        return state.with_moments(mean, covariance)

    def mean_truncation_error(self, state: GaussianState) -> float:
        """Relative error of the closure's drift against the exact ``E[F(X)]``.

        This is the quantity that ties the breakdown width of section 2 to the propagation gate:
        when it exceeds the contract's tolerance, :func:`propagate_local_gaussian` refuses.

        Args:
            state: Input Gaussian (one-dimensional).

        Returns:
            The relative error, or ``0.0`` when both quantities are exactly zero.
        """
        mu = float(state.mean[0])
        sigma = state.scalar_sigma
        point = np.asarray(mu, dtype=np.float64)
        exact = self.law.exact_expectation(mu, sigma)
        closure = (
            float(self.law.evaluate(point))
            + 0.5 * float(self.law.derivative(2, point)) * sigma * sigma
        )
        if exact == 0.0:
            return 0.0 if closure == 0.0 else math.inf
        return abs(exact - closure) / abs(exact)


def quadrature_moment_update(
    drift: Drift,
    state: GaussianState,
    *,
    dt: float,
    diffusion: Any,
    contract: MomentClosureContract,
) -> GaussianState:
    """Propagate by a deterministic Gauss-Hermite ensemble — a different computation, same truth.

    Exact for a polynomial drift whose transported degree does not exceed ``2 * nodes - 1``, so for
    the polynomial cases this agrees with :meth:`Drift.exact_moment_update` to round-off while being
    computed by an entirely different route (a weighted sum over transformed nodes rather than an
    algebraic moment formula).

    Args:
        drift: The drift.
        state: Input Gaussian.
        dt: Step in seconds.
        diffusion: Diffusion matrix.
        contract: Declared thresholds; supplies the node count.

    Returns:
        The propagated state.
    """
    step, noise = _validated_dt_and_diffusion(
        dt, diffusion, dimension=state.dimension, contract=contract
    )
    abscissas, raw = np.polynomial.hermite.hermgauss(contract.quadrature_nodes)
    weights_1d = raw / math.sqrt(math.pi)
    weights_1d = weights_1d / float(np.sum(weights_1d))
    dimension = state.dimension
    grids = np.meshgrid(*([abscissas] * dimension), indexing="ij")
    weight_grids = np.meshgrid(*([weights_1d] * dimension), indexing="ij")
    unit_points = np.stack([grid.ravel() for grid in grids], axis=1)
    weights = np.prod(np.stack([grid.ravel() for grid in weight_grids], axis=1), axis=1)
    weights = weights / float(np.sum(weights))
    points = state.mean + math.sqrt(2.0) * (unit_points @ state.matrix_sqrt().T)

    transported = points + drift.evaluate(points) * step
    mean = weights @ transported
    centered = transported - mean
    covariance = (centered * weights[:, None]).T @ centered + 2.0 * step * noise
    return state.with_moments(mean, 0.5 * (covariance + covariance.T))


def monte_carlo_moment_update(
    drift: Drift,
    state: GaussianState,
    *,
    dt: float,
    diffusion: Any,
    contract: MomentClosureContract,
    n_samples: int,
    seed: int,
) -> GaussianState:
    """Propagate by a seeded random ensemble.

    Present as a **counter-example, not a comparator**: its error is set by ``1/sqrt(n_samples)`` and
    on any affordable ensemble it is orders of magnitude larger than the closure error being
    measured, so scoring a closure against a sampled ensemble would hide the very effect this module
    exists to quantify.  The exact and quadrature updates are the comparators.

    Args:
        drift: The drift.
        state: Input Gaussian.
        dt: Step in seconds.
        diffusion: Diffusion matrix.
        contract: Declared thresholds.
        n_samples: Ensemble size.
        seed: RNG seed, so the number is reproducible.

    Returns:
        The propagated state estimated from the sample.
    """
    step, noise = _validated_dt_and_diffusion(
        dt, diffusion, dimension=state.dimension, contract=contract
    )
    count = _positive_int(n_samples, what="n_samples")
    rng = np.random.default_rng(_positive_int(seed, what="seed"))
    unit = rng.standard_normal(size=(count, state.dimension))
    points = state.mean + unit @ state.matrix_sqrt().T
    transported = points + drift.evaluate(points) * step
    mean = np.mean(transported, axis=0)
    centered = transported - mean
    covariance = centered.T @ centered / count + 2.0 * step * noise
    return state.with_moments(mean, 0.5 * (covariance + covariance.T))


def propagate_gaussian_ou_exact(
    drift: LinearDrift,
    state: GaussianState,
    *,
    dt: float,
    diffusion: Any,
    contract: MomentClosureContract,
) -> GaussianState:
    """Propagate an Ornstein-Uhlenbeck state by the exact continuous-time solution.

    Uses ``mu(t) = m + e^{At}(mu0 - m)`` and ``Sigma(t) = Sigma_inf + e^{At}(Sigma0 - Sigma_inf)
    e^{A^T t}`` with ``Sigma_inf`` from the Lyapunov equation ``A Sigma + Sigma A^T + 2 D = 0``.
    Comparing the Euler-Maruyama step against this separates DISCRETISATION error from CLOSURE
    error; conflating the two is how a closure gets blamed for a time step, or acquitted by one.

    Args:
        drift: The linear drift; must be stable.
        state: Input Gaussian.
        dt: Elapsed time in seconds.
        diffusion: Diffusion matrix.
        contract: Declared thresholds.

    Returns:
        The exactly propagated state.

    Raises:
        ValueError: If the drift matrix is not stable, in which case the stationary covariance the
            solution is written around does not exist.
    """
    if not isinstance(drift, LinearDrift):
        raise TypeError("the exact OU solution requires a LinearDrift")
    step, noise = _validated_dt_and_diffusion(
        dt, diffusion, dimension=state.dimension, contract=contract
    )
    eigenvalues = np.linalg.eigvals(drift.matrix)
    if float(np.max(eigenvalues.real)) >= 0.0:
        raise ValueError("the exact OU solution requires a stable drift matrix")
    stationary_mean = np.linalg.solve(drift.matrix, -drift.offset)
    stationary_covariance = solve_continuous_lyapunov(drift.matrix, -2.0 * noise)
    propagator = expm(drift.matrix * step)
    mean = stationary_mean + propagator @ (state.mean - stationary_mean)
    covariance = (
        stationary_covariance
        + propagator @ (state.covariance - stationary_covariance) @ propagator.T
    )
    return state.with_moments(mean, 0.5 * (covariance + covariance.T))


@dataclass(frozen=True, slots=True)
class ClosureComparison:
    """Closure vs exact vs quadrature vs Monte Carlo, one step, with every discrepancy measured."""

    drift_name: str
    dt: float
    exact: GaussianState
    closure: GaussianState
    quadrature: GaussianState
    monte_carlo: GaussianState | None
    closure_mean_error: float
    closure_covariance_error: float
    quadrature_mean_error: float
    quadrature_covariance_error: float
    monte_carlo_mean_error: float | None
    monte_carlo_covariance_error: float | None
    relative_closure_covariance_error: float
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able comparison."""
        return {
            "drift_name": self.drift_name,
            "dt": self.dt,
            "exact_mean": self.exact.mean.tolist(),
            "exact_covariance": self.exact.covariance.tolist(),
            "closure_mean": self.closure.mean.tolist(),
            "closure_covariance": self.closure.covariance.tolist(),
            "quadrature_mean": self.quadrature.mean.tolist(),
            "quadrature_covariance": self.quadrature.covariance.tolist(),
            "closure_mean_error": self.closure_mean_error,
            "closure_covariance_error": self.closure_covariance_error,
            "quadrature_mean_error": self.quadrature_mean_error,
            "quadrature_covariance_error": self.quadrature_covariance_error,
            "monte_carlo_mean_error": self.monte_carlo_mean_error,
            "monte_carlo_covariance_error": self.monte_carlo_covariance_error,
            "relative_closure_covariance_error": self.relative_closure_covariance_error,
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
        }


def compare_closure_against_ensemble(
    drift: Drift,
    state: GaussianState,
    *,
    dt: float,
    diffusion: Any,
    contract: MomentClosureContract,
    monte_carlo_samples: int | None,
    monte_carlo_seed: int | None,
) -> ClosureComparison:
    """Score a local-Gaussian step against the exact moments and a deterministic ensemble.

    Args:
        drift: The drift.
        state: Input Gaussian.
        dt: Step in seconds.
        diffusion: Diffusion matrix.
        contract: Declared thresholds.
        monte_carlo_samples: Optional sampled-ensemble size; ``None`` skips it.
        monte_carlo_seed: Seed for that ensemble; required whenever the size is given.

    Returns:
        The :class:`ClosureComparison`.

    Raises:
        ValueError: If a Monte-Carlo size is given without a seed.
    """
    step, noise = _validated_dt_and_diffusion(
        dt, diffusion, dimension=state.dimension, contract=contract
    )
    exact = drift.exact_moment_update(state, dt=step, diffusion=noise)
    closure = drift.taylor_closure_update(state, dt=step, diffusion=noise)
    quadrature = quadrature_moment_update(drift, state, dt=step, diffusion=noise, contract=contract)
    monte_carlo: GaussianState | None = None
    if monte_carlo_samples is not None:
        if monte_carlo_seed is None:
            raise ValueError("a Monte-Carlo ensemble needs an explicit seed to be reproducible")
        monte_carlo = monte_carlo_moment_update(
            drift,
            state,
            dt=step,
            diffusion=noise,
            contract=contract,
            n_samples=monte_carlo_samples,
            seed=monte_carlo_seed,
        )

    def _mean_error(candidate: GaussianState) -> float:
        return float(np.max(np.abs(candidate.mean - exact.mean)))

    def _covariance_error(candidate: GaussianState) -> float:
        return float(np.max(np.abs(candidate.covariance - exact.covariance)))

    scale = float(np.max(np.abs(exact.covariance)))
    closure_covariance_error = _covariance_error(closure)
    return ClosureComparison(
        drift_name=drift.name,
        dt=step,
        exact=exact,
        closure=closure,
        quadrature=quadrature,
        monte_carlo=monte_carlo,
        closure_mean_error=_mean_error(closure),
        closure_covariance_error=closure_covariance_error,
        quadrature_mean_error=_mean_error(quadrature),
        quadrature_covariance_error=_covariance_error(quadrature),
        monte_carlo_mean_error=None if monte_carlo is None else _mean_error(monte_carlo),
        monte_carlo_covariance_error=(
            None if monte_carlo is None else _covariance_error(monte_carlo)
        ),
        relative_closure_covariance_error=(
            closure_covariance_error / scale if scale > 0.0 else closure_covariance_error
        ),
    )


# ---------------------------------------------------------------------------------------------
# the bimodality measurement and the refusal it drives
# ---------------------------------------------------------------------------------------------


def boltzmann_density(law: AnalyticForceLaw, grid: Any, *, noise_intensity: float) -> np.ndarray:
    """Return the analytic stationary density of ``dX = F(X) dt + sqrt(2 D) dW`` on a grid.

    For a one-dimensional gradient system ``F = -U'`` the stationary density is
    ``rho ~ exp(-U(x) / D)``; this is closed form, so the bimodality statistics below are scored
    against truth rather than against a sampled histogram.

    Args:
        law: The force law; its potential is taken in closed form.
        grid: Strictly increasing positions.
        noise_intensity: ``D``.

    Returns:
        The unnormalised density on the grid, shifted so its maximum exponent is zero (which
        removes the overflow that a deep well otherwise guarantees).
    """
    if not isinstance(law, AnalyticForceLaw):
        raise TypeError("law must be an AnalyticForceLaw")
    positions = np.asarray(grid, dtype=np.float64).ravel()
    if positions.size < 3 or not np.all(np.diff(positions) > 0.0):
        raise ValueError("grid must be strictly increasing with at least three points")
    intensity = _positive(noise_intensity, what="noise_intensity")
    exponent = -law.potential(positions) / intensity
    return np.exp(exponent - float(np.max(exponent)))


def _local_maxima_with_prominence(density: np.ndarray) -> tuple[tuple[int, ...], tuple[float, ...]]:
    peak_indices = [
        index
        for index in range(1, density.size - 1)
        if density[index] > density[index - 1] and density[index] > density[index + 1]
    ]
    prominences: list[float] = []
    for index in peak_indices:
        higher_left = [
            other for other in peak_indices if other < index and density[other] > density[index]
        ]
        higher_right = [
            other for other in peak_indices if other > index and density[other] > density[index]
        ]
        left_key = max(higher_left) if higher_left else 0
        right_key = min(higher_right) if higher_right else density.size - 1
        left_min = float(np.min(density[left_key : index + 1]))
        right_min = float(np.min(density[index : right_key + 1]))
        prominences.append(float(density[index]) - max(left_min, right_min))
    return tuple(peak_indices), tuple(prominences)


@dataclass(frozen=True, slots=True)
class BimodalityReport:
    """Measured shape statistics of a reference density, and whether a Gaussian may represent it."""

    mode_count: int
    mode_positions: tuple[float, ...]
    mode_prominences: tuple[float, ...]
    bimodality_coefficient: float
    skewness: float
    kurtosis: float
    density_ratio_at_evaluation_point: float
    evaluation_point: float
    reference_mean: float
    reference_variance: float
    mass_audit: ProbabilityMassAudit
    failed_clauses: tuple[str, ...]
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    @property
    def is_gaussian_representable(self) -> bool:
        """Whether every declared clause passed."""
        return not self.failed_clauses

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able shape report."""
        return {
            "mode_count": self.mode_count,
            "mode_positions": list(self.mode_positions),
            "mode_prominences": list(self.mode_prominences),
            "bimodality_coefficient": self.bimodality_coefficient,
            "skewness": self.skewness,
            "kurtosis": self.kurtosis,
            "density_ratio_at_evaluation_point": self.density_ratio_at_evaluation_point,
            "evaluation_point": self.evaluation_point,
            "reference_mean": self.reference_mean,
            "reference_variance": self.reference_variance,
            "total_mass": self.mass_audit.total_mass,
            "mass_error": self.mass_audit.mass_error,
            "failed_clauses": list(self.failed_clauses),
            "is_gaussian_representable": self.is_gaussian_representable,
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
        }


def measure_bimodality(
    grid: Any,
    density: Any,
    *,
    evaluation_point: float,
    contract: MomentClosureContract,
) -> BimodalityReport:
    """Measure whether a single Gaussian may stand in for a reference density.

    Three statistics, all measured, none of them a flag:

    * the **mode count** with a declared prominence floor;
    * the **bimodality coefficient** ``(skew^2 + 1) / kurtosis``, whose Gaussian value is ``1/3``;
    * the **density ratio at the evaluation point**, ``rho(mu) / max(rho)``, which is the direct
      statement of "the closure's mean sits where the probability mass is not".

    The mass is audited by :func:`~aleph.virtual_cell.probability.audit_probability_mass` on the
    rectangle sum after normalising by the trapezoid rule.  Those two differ by exactly the boundary
    term, so a grid that does not contain its own tails RAISES here instead of being renormalised
    into a wrong answer.

    Args:
        grid: Strictly increasing, uniformly spaced positions.
        density: Nonnegative density values on that grid; need not be normalised.
        evaluation_point: Where a Gaussian closure would put its mean.
        contract: Declared thresholds.

    Returns:
        The :class:`BimodalityReport`.

    Raises:
        ValueError: If the grid is not uniform, the density is negative, or the grid fails the mass
            audit.
    """
    if not isinstance(contract, MomentClosureContract):
        raise TypeError("contract must be a MomentClosureContract")
    positions = np.asarray(grid, dtype=np.float64).ravel()
    values = np.asarray(density, dtype=np.float64).ravel()
    if positions.size < 5 or positions.size != values.size:
        raise ValueError("grid and density must be the same length and have at least five points")
    spacing = np.diff(positions)
    if not np.all(spacing > 0.0):
        raise ValueError("grid must be strictly increasing")
    step = float(spacing[0])
    if float(np.max(np.abs(spacing - step))) > 1e-9 * step:
        raise ValueError("grid must be uniformly spaced for the rectangle/trapezoid mass audit")
    if np.any(values < 0.0):
        raise ValueError("a density must not be negative")
    if not np.all(np.isfinite(values)):
        raise ValueError("density must be finite")

    normalization = float(np.trapezoid(values, positions))
    if normalization <= 0.0:
        raise ValueError("density integrates to zero")
    normalized = values / normalization
    mass = normalized * step
    audit = audit_probability_mass(
        mass,
        normalization_tolerance=contract.mass_normalization_tolerance,
        negativity_tolerance=contract.mass_negativity_tolerance,
    )

    weights = mass / float(np.sum(mass))
    mean = float(np.dot(weights, positions))
    centered = positions - mean
    variance = float(np.dot(weights, centered**2))
    if variance <= 0.0:
        raise ValueError("reference density has zero variance")
    skewness = float(np.dot(weights, centered**3)) / variance**1.5
    kurtosis = float(np.dot(weights, centered**4)) / variance**2
    coefficient = (skewness * skewness + 1.0) / kurtosis

    peak = float(np.max(normalized))
    indices, prominences = _local_maxima_with_prominence(normalized)
    kept = [
        (index, prominence)
        for index, prominence in zip(indices, prominences, strict=True)
        if prominence >= contract.mode_prominence_ratio * peak
    ]
    mode_positions = tuple(float(positions[index]) for index, _ in kept)
    mode_prominences = tuple(float(prominence) for _, prominence in kept)

    point = _finite(evaluation_point, what="evaluation_point")
    if point < positions[0] or point > positions[-1]:
        raise ValueError("evaluation_point lies outside the reference grid")
    ratio = float(np.interp(point, positions, normalized)) / peak

    failed: list[str] = []
    if len(kept) > contract.max_admissible_modes:
        failed.append(
            f"mode-count {len(kept)} > max_admissible_modes {contract.max_admissible_modes}"
        )
    if coefficient > contract.bimodality_coefficient_max:
        failed.append(
            f"bimodality-coefficient {coefficient:.6g} > max {contract.bimodality_coefficient_max:.6g}"
        )
    if ratio < contract.modal_density_ratio_min:
        failed.append(
            f"density at the closure mean is {ratio:.6g} of the modal density, below the declared "
            f"minimum {contract.modal_density_ratio_min:.6g}"
        )

    return BimodalityReport(
        mode_count=len(kept),
        mode_positions=mode_positions,
        mode_prominences=mode_prominences,
        bimodality_coefficient=coefficient,
        skewness=skewness,
        kurtosis=kurtosis,
        density_ratio_at_evaluation_point=ratio,
        evaluation_point=point,
        reference_mean=mean,
        reference_variance=variance,
        mass_audit=audit,
        failed_clauses=tuple(failed),
    )


@dataclass(frozen=True, slots=True)
class GaussianClosureStep:
    """The propagated local Gaussian.  It exists only inside an ADMITTED verdict."""

    state: GaussianState
    dt: float
    drift_name: str
    representation: NodeRepresentation = NodeRepresentation.LOCAL_GAUSSIAN
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY


@dataclass(frozen=True, slots=True)
class ClosureVerdict:
    """Whether the local-Gaussian step was issued, and what stopped it if not.

    A refused verdict carries **no mean**: :attr:`step` is ``None`` and :meth:`require_step` raises.
    That is deliberate — the prohibition this implements is "a bimodal event state must not be
    hidden in a single Gaussian", and a report that refused while still exposing the mean would let
    a caller launder exactly the number the refusal is about.
    """

    kind: ClosureVerdictKind
    drift_name: str
    step: GaussianClosureStep | None
    bimodality: BimodalityReport | None
    mean_truncation_error: float
    reasons: tuple[str, ...]
    expansions: tuple[ExpansionRoute, ...]
    fallback_representation: RepresentationKind | None
    contract: MomentClosureContract
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    @property
    def admitted(self) -> bool:
        """Whether a propagated state was issued."""
        return self.kind is ClosureVerdictKind.ADMITTED

    def require_step(self) -> GaussianClosureStep:
        """Return the propagated step, or raise if the closure refused.

        Returns:
            The admitted :class:`GaussianClosureStep`.

        Raises:
            ClosureRefusedError: If the verdict is any refusal.
        """
        if self.step is None:
            detail = "; ".join(self.reasons) if self.reasons else "no reason recorded"
            raise ClosureRefusedError(f"{self.kind.value}: {detail}")
        return self.step

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able verdict."""
        return {
            "kind": self.kind.value,
            "drift_name": self.drift_name,
            "admitted": self.admitted,
            "mean": None if self.step is None else self.step.state.mean.tolist(),
            "covariance": None if self.step is None else self.step.state.covariance.tolist(),
            "bimodality": None if self.bimodality is None else self.bimodality.as_dict(),
            "mean_truncation_error": self.mean_truncation_error,
            "reasons": list(self.reasons),
            "expansions": [route.value for route in self.expansions],
            "fallback_representation": (
                None if self.fallback_representation is None else self.fallback_representation.value
            ),
            "contract": self.contract.to_dict(),
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
            "engine_node_representation_today": ENGINE_NODE_REPRESENTATION_TODAY,
        }


def propagate_local_gaussian(
    drift: Drift,
    state: GaussianState,
    *,
    dt: float,
    diffusion: Any,
    contract: MomentClosureContract,
    reference_grid: Any | None = None,
    reference_density: Any | None = None,
) -> ClosureVerdict:
    """Propagate ``(mu, Sigma)`` one step, or refuse and say which representation is needed instead.

    Three gates, all measured:

    1. **Shape evidence.**  A drift under which the Gaussian family is not closed requires a
       reference density; without one, the closure is refused.  Applying a Gaussian closure while
       knowing nothing about the true density's shape is not a small approximation, it is the
       untested assumption the whole gate exists to catch.
    2. **Bimodality.**  If the reference density fails any declared clause of
       :func:`measure_bimodality`, the closure is refused and routed to a mixture or an ensemble.
    3. **Width.**  If the closure's own mean drift departs from the exact ``E[F(X)]`` by more than
       ``contract.expansion_relative_error_tolerance``, the state is wider than the range over which
       a second-order closure was validated, and the closure is refused.

    Args:
        drift: The drift.
        state: Input Gaussian.
        dt: Step in seconds.
        diffusion: Diffusion matrix.
        contract: Declared thresholds.
        reference_grid: Optional positions of a reference density.
        reference_density: Optional reference density values on that grid.

    Returns:
        The :class:`ClosureVerdict`.

    Raises:
        ValueError: If exactly one of the reference grid/density is supplied.
    """
    if not isinstance(contract, MomentClosureContract):
        raise TypeError("contract must be a MomentClosureContract")
    if (reference_grid is None) != (reference_density is None):
        raise ValueError("a reference density needs both its grid and its values")
    step, noise = _validated_dt_and_diffusion(
        dt, diffusion, dimension=state.dimension, contract=contract
    )
    truncation = drift.mean_truncation_error(state)

    bimodality: BimodalityReport | None = None
    if reference_grid is not None:
        if state.dimension != 1:
            raise ValueError("the reference-density gate is one-dimensional")
        bimodality = measure_bimodality(
            reference_grid,
            reference_density,
            evaluation_point=float(state.mean[0]),
            contract=contract,
        )

    if not drift.is_gaussian_closed and bimodality is None:
        return ClosureVerdict(
            kind=ClosureVerdictKind.REFUSED_UNVERIFIED_SHAPE,
            drift_name=drift.name,
            step=None,
            bimodality=None,
            mean_truncation_error=truncation,
            reasons=(
                "the Gaussian family is not closed under this drift and no shape evidence was "
                "supplied; a closure applied without any measurement of the true density's shape "
                "is the assumption this gate exists to refuse",
            ),
            expansions=(ExpansionRoute.ENSEMBLE,),
            fallback_representation=RepresentationKind.EMPIRICAL_ENSEMBLE,
            contract=contract,
        )

    if bimodality is not None and not bimodality.is_gaussian_representable:
        return ClosureVerdict(
            kind=ClosureVerdictKind.REFUSED_BIMODAL,
            drift_name=drift.name,
            step=None,
            bimodality=bimodality,
            mean_truncation_error=truncation,
            reasons=bimodality.failed_clauses,
            expansions=(ExpansionRoute.ENSEMBLE, ExpansionRoute.REGIME_PARTITION),
            fallback_representation=RepresentationKind.GAUSSIAN_MIXTURE,
            contract=contract,
        )

    if truncation > contract.expansion_relative_error_tolerance:
        return ClosureVerdict(
            kind=ClosureVerdictKind.REFUSED_WIDTH_BEYOND_VALIDATED_RANGE,
            drift_name=drift.name,
            step=None,
            bimodality=bimodality,
            mean_truncation_error=truncation,
            reasons=(
                f"the second-order closure's mean drift is wrong by {truncation:.6g} relative, "
                f"above the declared {contract.expansion_relative_error_tolerance:.6g}; the state "
                "is wider than the range over which this closure was validated",
            ),
            expansions=(ExpansionRoute.ENSEMBLE, ExpansionRoute.LATENT_COORDINATE),
            fallback_representation=RepresentationKind.EMPIRICAL_ENSEMBLE,
            contract=contract,
        )

    propagated = drift.taylor_closure_update(state, dt=step, diffusion=noise)
    return ClosureVerdict(
        kind=ClosureVerdictKind.ADMITTED,
        drift_name=drift.name,
        step=GaussianClosureStep(state=propagated, dt=step, drift_name=drift.name),
        bimodality=bimodality,
        mean_truncation_error=truncation,
        reasons=(),
        expansions=(),
        fallback_representation=None,
        contract=contract,
    )


# ---------------------------------------------------------------------------------------------
# 4. source-factored, not one Sigma
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SourceFactoredCovariance:
    """Per-source covariance blocks that are never pre-merged.

    The merge is available but it is a named, explicit operation
    (:meth:`merged_total`), so a report that used a single total can be told apart from one that
    propagated the sources separately.  §2.2 and audit A1: merging thermal, quenched-structural,
    cell-to-cell and parameter-epistemic contributions into one ``Sigma`` destroys both the cause and
    the time scale, and the measurements below show it also destroys the answer.
    """

    blocks: Mapping[UncertaintySource, np.ndarray]
    axes: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.blocks, Mapping) or not self.blocks:
            raise ValueError("at least one source block is required")
        axes = tuple(_nonempty(axis, what="axis") for axis in self.axes)
        if len(set(axes)) != len(axes):
            raise ValueError("axes must be unique")
        frozen: dict[UncertaintySource, np.ndarray] = {}
        for source, block in self.blocks.items():
            if not isinstance(source, UncertaintySource):
                raise TypeError("every block key must be an UncertaintySource")
            matrix = _frozen_matrix(block, what=f"block[{source.value}]", dimension=len(axes))
            _check_symmetric_psd(
                matrix,
                what=f"block[{source.value}]",
                symmetry_tolerance=0.0,
                minimum_eigenvalue_tolerance=1e-12,
            )
            frozen[source] = matrix
        object.__setattr__(self, "blocks", MappingProxyType(frozen))
        object.__setattr__(self, "axes", axes)

    @property
    def sources(self) -> tuple[UncertaintySource, ...]:
        """The declared sources, in insertion order."""
        return tuple(self.blocks)

    def merged_total(self) -> np.ndarray:
        """Return the summed covariance — the representation this module argues against.

        Returns:
            The total covariance.  Named so that using it is a visible decision.
        """
        total = np.zeros((len(self.axes), len(self.axes)), dtype=np.float64)
        for block in self.blocks.values():
            total = total + block
        return total


@dataclass(frozen=True, slots=True)
class PlacementContrast:
    """Two uncertainty placements at the SAME total variance, and the force expectations they give."""

    description: str
    total_variance: float
    first_placement: str
    second_placement: str
    first_expectation: float
    second_expectation: float
    point_force: float
    difference: float
    ratio: float
    closed_form_ratio: float
    sources: tuple[UncertaintySource, ...]
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able contrast."""
        return {
            "description": self.description,
            "total_variance": self.total_variance,
            "first_placement": self.first_placement,
            "second_placement": self.second_placement,
            "first_expectation": self.first_expectation,
            "second_expectation": self.second_expectation,
            "point_force": self.point_force,
            "difference": self.difference,
            "ratio": self.ratio,
            "closed_form_ratio": self.closed_form_ratio,
            "sources": [source.value for source in self.sources],
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
        }


def measure_rate_placement_contrast(
    *, amplitude: float, rate: float, position: float, total_variance: float
) -> PlacementContrast:
    """Same total variance on the position or on the rate of ``A exp(b x)``: two different forces.

    Both expectations are exact Gaussian MGFs::

        variance on the position:  A exp(b mu + b^2 v / 2)
        variance on the rate:      A exp(b mu + mu^2 v / 2)

    so the ratio is ``exp((mu^2 - b^2) v / 2)`` and neither correction is zero.  This is the
    concrete argument for source factoring: a single total variance does not determine the mean
    force, because the mean force depends on WHERE the variance enters the nonlinearity.

    Args:
        amplitude: ``A``.
        rate: ``b``.
        position: ``mu``.
        total_variance: ``v``, identical for both placements.

    Returns:
        The :class:`PlacementContrast`.
    """
    a = _finite(amplitude, what="amplitude")
    b = _finite(rate, what="rate")
    mu = _finite(position, what="position")
    v = _nonnegative(total_variance, what="total_variance")

    point = a * math.exp(b * mu)
    on_position = a * math.exp(b * mu + 0.5 * b * b * v)
    on_rate = a * math.exp(b * mu + 0.5 * mu * mu * v)
    return PlacementContrast(
        description=(
            "A*exp(b*x): the same total variance placed on the position (thermal) or on the rate "
            "(parameter-epistemic) gives different mean forces"
        ),
        total_variance=v,
        first_placement="thermal variance on the position x",
        second_placement="epistemic variance on the rate b",
        first_expectation=on_position,
        second_expectation=on_rate,
        point_force=point,
        difference=on_position - on_rate,
        ratio=on_position / on_rate,
        closed_form_ratio=math.exp(0.5 * (b * b - mu * mu) * v),
        sources=(UncertaintySource.THERMAL, UncertaintySource.PARAMETER),
    )


def measure_cubic_placement_contrast(
    *, stiffness: float, cubic_coefficient: float, position: float, total_variance: float
) -> PlacementContrast:
    """Same total variance on the position or on the stiffness of ``-(k x + a x^3)``.

    The force is CUBIC in ``x`` and LINEAR in ``k``, so variance on the position produces the Jensen
    correction ``-3 a mu v`` while variance on the stiffness produces exactly none.  The sharpest
    possible form of the argument: one of the two corrections is identically zero at the same total
    variance, so a pre-merged ``Sigma`` cannot even get the sign of the correction from its
    magnitude.

    Args:
        stiffness: ``k``.
        cubic_coefficient: ``a``.
        position: ``mu``.
        total_variance: ``v``.

    Returns:
        The :class:`PlacementContrast`.
    """
    k = _finite(stiffness, what="stiffness")
    a = _finite(cubic_coefficient, what="cubic_coefficient")
    mu = _finite(position, what="position")
    v = _nonnegative(total_variance, what="total_variance")

    point = -(k * mu + a * mu**3)
    on_position = -(k * mu + a * (mu**3 + 3.0 * mu * v))
    on_stiffness = point
    difference = on_position - on_stiffness
    independent_ratio = (
        (k * mu + a * mu**3 + 3.0 * a * mu * v) / (k * mu + a * mu**3)
        if (k * mu + a * mu**3) != 0.0
        else math.inf
    )
    return PlacementContrast(
        description=(
            "-(k x + a x^3): position variance yields the Jensen correction -3 a mu v; the same "
            "variance on the stiffness yields exactly zero, because the force is linear in k"
        ),
        total_variance=v,
        first_placement="thermal variance on the position x",
        second_placement="epistemic variance on the stiffness k",
        first_expectation=on_position,
        second_expectation=on_stiffness,
        point_force=point,
        difference=difference,
        ratio=on_position / on_stiffness if on_stiffness != 0.0 else math.inf,
        closed_form_ratio=independent_ratio,
        sources=(UncertaintySource.THERMAL, UncertaintySource.PARAMETER),
    )


@dataclass(frozen=True, slots=True)
class QuenchedThermalSplit:
    """A quenched/thermal split of one total variance across ``N`` exponential bonds.

    The quenched offset is shared by every bond and does not average away; the thermal part is
    independent per bond and averages as ``1/N``.  The measured consequence is the point of the
    class: the MEAN aggregate force is independent of the split, so the usual "keep the sources
    separate because the mean moves" argument is FALSE here — while the VARIANCE of the aggregate
    changes by exactly ``N`` between the two limits, so any downstream nonlinear functional of the
    aggregate differs.
    """

    n_units: int
    quenched_variance: float
    thermal_variance: float
    total_variance: float
    mean_aggregate_force: float
    variance_aggregate_force: float
    second_moment_aggregate_force: float
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    def as_dict(self) -> dict[str, Any]:
        """Return the JSON-able split."""
        return {
            "n_units": self.n_units,
            "quenched_variance": self.quenched_variance,
            "thermal_variance": self.thermal_variance,
            "total_variance": self.total_variance,
            "mean_aggregate_force": self.mean_aggregate_force,
            "variance_aggregate_force": self.variance_aggregate_force,
            "second_moment_aggregate_force": self.second_moment_aggregate_force,
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
        }


def measure_quenched_thermal_split(
    *,
    amplitude: float,
    rate: float,
    position: float,
    quenched_variance: float,
    thermal_variance: float,
    n_units: int,
) -> QuenchedThermalSplit:
    """Closed-form mean and variance of the mean force over ``N`` bonds with a split variance.

    With ``x_i = mu + q + t_i``, ``q ~ N(0, vq)`` shared and ``t_i ~ N(0, vt)`` independent, and
    ``F = A exp(b x)``::

        E[Fbar]   = A e^{b mu} e^{b^2 (vq + vt) / 2}
        E[Fbar^2] = A^2 e^{2 b mu} e^{2 b^2 vq} ( e^{2 b^2 vt} / N + (1 - 1/N) e^{b^2 vt} )

    Args:
        amplitude: ``A``.
        rate: ``b``.
        position: ``mu``.
        quenched_variance: ``vq``, shared across bonds.
        thermal_variance: ``vt``, independent per bond.
        n_units: ``N``, the number of bonds averaged.

    Returns:
        The :class:`QuenchedThermalSplit`.
    """
    a = _finite(amplitude, what="amplitude")
    b = _finite(rate, what="rate")
    mu = _finite(position, what="position")
    vq = _nonnegative(quenched_variance, what="quenched_variance")
    vt = _nonnegative(thermal_variance, what="thermal_variance")
    n = _positive_int(n_units, what="n_units")

    mean = a * math.exp(b * mu) * math.exp(0.5 * b * b * (vq + vt))
    second = (
        a
        * a
        * math.exp(2.0 * b * mu)
        * math.exp(2.0 * b * b * vq)
        * (math.exp(2.0 * b * b * vt) / n + (1.0 - 1.0 / n) * math.exp(b * b * vt))
    )
    return QuenchedThermalSplit(
        n_units=n,
        quenched_variance=vq,
        thermal_variance=vt,
        total_variance=vq + vt,
        mean_aggregate_force=mean,
        variance_aggregate_force=second - mean * mean,
        second_moment_aggregate_force=second,
    )


# ---------------------------------------------------------------------------------------------
# pre-registration card
# ---------------------------------------------------------------------------------------------


def stochastic_node_experiment_card(manifest_hashes: Sequence[str]) -> SandboxExperimentCard:
    """Return the pre-registration card for this oracle.

    Args:
        manifest_hashes: Manifest hashes the card is declared against.

    Returns:
        The :class:`~aleph.virtual_cell.contracts.SandboxExperimentCard`.
    """
    return SandboxExperimentCard(
        experiment_id="D1-stochastic-node-moment-oracle",
        question=(
            "if a node carried a distribution instead of a point, how large is the force error made "
            "by evaluating the force once at the mean, how does it scale with the width, at what "
            "width does a Gaussian closure stop being safe, and where must the closure refuse "
            "outright?"
        ),
        manifest_hashes=tuple(manifest_hashes),
        intervention=(
            "replace the point node x by a distribution of declared width and source composition, "
            "and evaluate the force from the distribution's moments instead of at its mean; the "
            "Warp runtime is not touched and carries points throughout"
        ),
        observables=(
            "signed Jensen gap E[F(X)] - F(E[X]) against the closed form",
            "fitted log|gap| vs log sigma exponent against the analytically predicted order",
            "relative error of the second-order truncation against the exact expectation",
            "breakdown width at a declared closure tolerance",
            "closure vs exact vs quadrature vs Monte-Carlo one-step moments",
            "mode count, bimodality coefficient and density ratio at the closure mean",
            "mean and variance of an aggregate force under a quenched/thermal split",
        ),
        positive_controls=(
            "affine force law: the gap must be exactly 0.0, not merely small",
            "Ornstein-Uhlenbeck: the local-Gaussian closure must reproduce the exact moment update "
            "to round-off, because the Gaussian family is closed under affine drift",
            "cubic force law: the second-order expansion is an identity, so its truncation error "
            "must be exactly zero at every width",
        ),
        negative_controls=(
            "quartic force at mu = 0: F''(mu) = 0, so the fitted exponent must be 4 and the "
            "second-order expansion must be wrong by 100 % at every width",
            "Monte-Carlo ensemble: its error must be orders of magnitude larger than the closure "
            "error, demonstrating that a sampled ensemble could not have scored this closure",
            "zero width: every law must return exactly zero gap",
        ),
        adversarial_controls=(
            "double-well potential: the closure mean sits on the barrier and the closure must "
            "REFUSE, driven by measured statistics rather than a hand-set flag",
            "nonlinear drift with no shape evidence supplied: the closure must refuse rather than "
            "assume unimodality",
            "truncated density grid: the probability-mass audit must raise rather than renormalise",
            "two uncertainty sources with identical total variance: the mean forces must differ, and "
            "in the quenched/thermal case they must NOT, which is the harder half of the claim",
        ),
        held_out_interventions=(
            "Morse bond evaluated at its own minimum, where the point force is exactly zero and the "
            "distributional force is strictly positive",
            "quenched/thermal split at fixed total variance across a bond count not used to build "
            "the closed forms",
            "a decay length outside the one used to tabulate the breakdown widths, where the "
            "predicted safe width must still be the same multiple of lambda",
        ),
        falsifier=(
            "the affine control returns a nonzero gap, OR the fitted exponent departs from the "
            "analytically predicted leading even order, OR the local-Gaussian closure fails to "
            "reproduce the exact Ornstein-Uhlenbeck moment update to round-off, OR the double-well "
            "case is ADMITTED"
        ),
        failure_expansions=(
            "closure error above tolerance -> Gaussian mixture, then empirical ensemble",
            "measured bimodality -> mixture or regime partition, never a wider single Gaussian",
            "sources indistinguishable at equal total variance -> a source-factored propagation "
            "law is needed, not a larger Sigma",
        ),
        budget_class="cpu-analytic-oracle",
        wet_lab=False,
    )
