r"""Finite-difference sensitivity, the Fisher information matrix, and the sloppiness spectrum.

WHY THIS IS THE SECOND MEASUREMENT.  Under the 2026-07-25 reframe the project's PURPOSE is to infer
per-cell-type molecular parameters — parameters are OUTPUTS.  Inference is only meaningful for
parameters the observables actually respond to, and ``∂observable/∂parameter`` is that response.  It
needs no posterior, no likelihood and no sampler: a central finite difference on a deterministic
observable is enough, and it answers the questions that decide what to build next.

* **Which parameters are recoverable at all.**  A parameter whose column of ``J`` is at round-off is
  invisible to this observable; no amount of inference machinery will recover it, and no experiment
  using this observable should be designed to try.
* **Which parameters are confounded.**  Two nearly-parallel columns of ``J`` mean only their
  combination is identifiable — the classic sloppy direction.  Reporting the FIM eigenvectors names
  the combination, so the confound is stated rather than discovered later as a fit that will not move.
* **The effective dimension.**  The participation ratio of the FIM eigenvalue spectrum counts how many
  parameter directions the data constrains, which is the honest ceiling on how many parameters may be
  claimed as inferred from it.
* **The reference for any adjoint gradient.**  An adjoint is fast and easy to get wrong; this is the
  slow, obviously-correct thing to check it against.

DIMENSIONLESS BY CONSTRUCTION.  The Jacobian is taken in LOG-LOG form,
``J_ij = ∂ log o_i / ∂ log θ_j``, i.e. the relative change in observable ``i`` per relative change in
parameter ``j``.  Two reasons, both from the plan's §3 lesson 3: it makes the entries comparable across
parameters carrying different units (a stiffness in pN/µm next to a rate in 1/s), and a dimensionless
formulation is what makes the ``10⁶×`` unit-slip class — which has happened twice in this project, once
faking a compaction result — structurally unable to hide.

THE STEP IS DERIVED, NOT TUNED.  For a central difference the total error is
``O(h²·|f'''|) + O(ε·|f|/h)``, minimised at ``h ∝ ε^(1/3)``; we use the relative step ``ε^(1/3)``
(≈ 6.06e-6), the textbook choice (Nocedal & Wright, *Numerical Optimization*, §8.1).  It is not
adjusted per parameter and cannot be chosen to make a sensitivity appear.

Sanity Gate:
    * dimensional: log-log entries are dimensionless by construction; the FIM built from them is
      dimensionless, so its eigenvalues are directly comparable across parameters.
    * boundary: a parameter or observable at exactly zero has no logarithm, so it is rejected with a
      naming error rather than silently regularised — a "small number added for stability" here would
      be an undeclared magic constant in the middle of an inference claim.
    * conservation/invariant: the FIM is symmetric positive semi-definite by construction (``JᵀJ``);
      any negative eigenvalue is numerical and is reported against the same backward-stability floor
      the spectrum module uses.
    * numerical: central differences are second-order; the step is the derived ``ε^(1/3)`` optimum.
      Non-finite observable evaluations abort rather than being dropped, because dropping a failed
      evaluation silently reports a smaller Jacobian.
    * sign-sense: not applicable — no force is computed here.
    * measurement-protocol: the observable callable must be DETERMINISTIC at fixed parameters.  A
      stochastic observable (any KMC path) makes a finite difference measure seed noise, so a
      sensitivity over an event-driven observable requires a frozen event stream, which is the
      caller's responsibility and is recorded in the returned detail.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence

import numpy as np

__all__ = [
    "FisherReport",
    "LogSensitivity",
    "fisher_report",
    "log_log_jacobian",
]

_EPS64: float = float(np.finfo(np.float64).eps)

#: Optimal RELATIVE central-difference step: the ``ε^(1/3)`` minimiser of truncation + round-off error.
CENTRAL_DIFFERENCE_RELATIVE_STEP: float = float(np.cbrt(_EPS64))

#: An observable evaluation: parameter mapping -> a 1-D array of strictly positive observable values.
ObservableFn = Callable[[Mapping[str, float]], np.ndarray]


@dataclass(frozen=True, slots=True)
class LogSensitivity:
    """The log-log Jacobian of one observable vector with respect to one parameter set.

    Attributes:
        jacobian: ``(n_observables, n_parameters)`` dimensionless matrix
            ``∂ log o_i / ∂ log θ_j``.
        parameter_names: Column order.
        baseline: The observable at the unperturbed parameters.
        relative_step: The relative step each parameter was perturbed by.
        column_norms: ``‖J[:, j]‖₂`` per parameter — the single number saying how much this observable
            responds to that parameter at all.
    """

    jacobian: np.ndarray
    parameter_names: tuple[str, ...]
    baseline: np.ndarray
    relative_step: float
    column_norms: np.ndarray

    def as_artifact_fields(self) -> dict[str, Any]:
        """Return the JSON-able record: the response norm per parameter, ranked."""
        ranked = sorted(
            zip(self.parameter_names, (float(v) for v in self.column_norms)),
            key=lambda pair: pair[1],
            reverse=True,
        )
        return {
            "n_observables": int(self.jacobian.shape[0]),
            "n_parameters": int(self.jacobian.shape[1]),
            "relative_step": self.relative_step,
            "log_response_norm_by_parameter": {name: value for name, value in ranked},
        }


def log_log_jacobian(
    observable: ObservableFn,
    parameters: Mapping[str, float],
    *,
    names: Sequence[str] | None = None,
    relative_step: float = CENTRAL_DIFFERENCE_RELATIVE_STEP,
) -> LogSensitivity:
    """Return ``∂ log observable / ∂ log parameter`` by central finite differences.

    Args:
        observable: Deterministic map from a parameter mapping to a 1-D array of STRICTLY POSITIVE
            observable values (a logarithm is taken).  A spectrum is the intended input: it is a field,
            not a summary scalar, and it is positive on a positive-definite operator.
        parameters: The baseline parameter values.  Each perturbed parameter must be nonzero.
        names: Which parameters to differentiate with respect to; defaults to every key of
            ``parameters``, in iteration order.
        relative_step: Relative perturbation.  The default is the derived optimum; overriding it is a
            convergence STUDY (halve it and check the Jacobian does not move), never a tuning knob.

    Returns:
        The :class:`LogSensitivity`.

    Raises:
        ValueError: If a requested name is absent, a parameter is zero, the step is not in ``(0, 1)``,
            or any observable evaluation is non-positive, non-finite, or changes length.
    """
    if not (0.0 < relative_step < 1.0):
        raise ValueError(f"relative_step must lie in (0, 1); got {relative_step!r}")
    keys = tuple(names) if names is not None else tuple(parameters)
    missing = [name for name in keys if name not in parameters]
    if missing:
        raise ValueError(f"parameters missing the requested names {missing}")

    def evaluate(mapping: Mapping[str, float], label: str) -> np.ndarray:
        values = np.asarray(observable(mapping), np.float64).reshape(-1)
        if values.size == 0:
            raise ValueError(f"observable returned an empty vector at {label}")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"observable returned a non-finite value at {label}")
        if np.any(values <= 0.0):
            raise ValueError(
                f"observable returned a non-positive value at {label}; a log-log sensitivity needs a "
                "strictly positive observable, and clamping one here would be an undeclared constant"
            )
        return values

    baseline = evaluate(parameters, "baseline")
    jacobian = np.zeros((baseline.size, len(keys)), np.float64)
    for column, name in enumerate(keys):
        theta = float(parameters[name])
        if theta == 0.0:
            raise ValueError(f"parameter {name!r} is zero, so a relative perturbation is undefined")
        step = abs(theta) * relative_step
        up = dict(parameters)
        up[name] = theta + step
        down = dict(parameters)
        down[name] = theta - step
        high = evaluate(up, f"{name}+")
        low = evaluate(down, f"{name}-")
        if high.size != baseline.size or low.size != baseline.size:
            raise ValueError(f"observable changed length while perturbing {name!r}")
        # d(log o)/d(log θ) = (θ/o)·(do/dθ); in central-difference form the θ factors cancel to
        # (log o_+ − log o_−) / (log θ_+ − log θ_−), which is what is evaluated here.
        jacobian[:, column] = (np.log(high) - np.log(low)) / (
            np.log(theta + step) - np.log(theta - step)
        )

    return LogSensitivity(
        jacobian=jacobian,
        parameter_names=keys,
        baseline=baseline,
        relative_step=float(relative_step),
        column_norms=np.linalg.norm(jacobian, axis=0),
    )


@dataclass(frozen=True, slots=True)
class FisherReport:
    """The Fisher information of a log-log Jacobian and the sloppiness it exposes.

    Attributes:
        parameter_names: Parameter order, matching the FIM rows/columns.
        fisher: ``JᵀJ``, ``(p, p)``, symmetric positive semi-definite, dimensionless.
        eigenvalues: FIM eigenvalues, DESCENDING — the sloppiness spectrum.
        eigenvectors: Column-wise FIM eigenvectors; column ``k`` names the parameter COMBINATION that
            eigenvalue ``k`` constrains.
        resolvable_floor: ``p·ε·λ_max`` — the backward-stability floor below which an eigenvalue is
            not distinguishable from zero.
        n_resolvable_directions: How many eigenvalues sit above that floor: a hard upper bound on the
            number of parameters this observable can constrain.
        effective_dimension: Participation ratio of the normalised eigenvalue spectrum, a continuous
            count of the well-constrained directions (a spectrum with one dominant eigenvalue gives
            ≈1 however many parameters were varied).
        sloppiness_decades: ``log10(λ_max/λ_min_resolvable)``; the classic sloppy-model signature is
            several decades of roughly uniform spacing.
    """

    parameter_names: tuple[str, ...]
    fisher: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    resolvable_floor: float
    n_resolvable_directions: int
    effective_dimension: float
    sloppiness_decades: float | None

    def as_artifact_fields(self, *, n_directions: int = 4) -> dict[str, Any]:
        """Return the JSON-able record, naming the stiffest parameter combinations.

        Args:
            n_directions: How many leading eigen-directions to name by their parameter weights.

        Returns:
            A dict of plain Python types.
        """
        directions = []
        for k in range(min(int(n_directions), self.eigenvalues.size)):
            weights = self.eigenvectors[:, k]
            ranked = sorted(
                zip(self.parameter_names, (float(v) for v in weights)),
                key=lambda pair: abs(pair[1]),
                reverse=True,
            )
            directions.append(
                {
                    "eigenvalue": float(self.eigenvalues[k]),
                    "parameter_weights": {name: value for name, value in ranked},
                }
            )
        return {
            "parameters": list(self.parameter_names),
            "eigenvalues": [float(v) for v in self.eigenvalues],
            "resolvable_floor": self.resolvable_floor,
            "n_resolvable_directions": self.n_resolvable_directions,
            "effective_dimension": self.effective_dimension,
            "sloppiness_decades": self.sloppiness_decades,
            "leading_directions": directions,
        }


def fisher_report(sensitivity: LogSensitivity) -> FisherReport:
    """Build the Fisher information matrix of a log-log Jacobian and diagonalise it.

    With independent observations of equal relative uncertainty the Fisher information in
    log-parameters is ``JᵀJ`` up to that common variance factor.  The factor is deliberately NOT
    supplied: it scales every eigenvalue identically, so it changes no eigenvector, no effective
    dimension and no ratio reported here — and inventing a measurement noise level for observables
    this project has not yet measured experimentally would be exactly the fabricated constant the
    magic-number rule forbids.  What is reported is therefore information UP TO that scale, which is
    what the design questions need.

    Args:
        sensitivity: The :class:`LogSensitivity` to build the information matrix from.

    Returns:
        The :class:`FisherReport`.

    Raises:
        ValueError: If the Jacobian is empty.
    """
    jacobian = np.asarray(sensitivity.jacobian, np.float64)
    if jacobian.size == 0:
        raise ValueError("cannot build a Fisher information matrix from an empty Jacobian")

    fisher = jacobian.T @ jacobian
    fisher = 0.5 * (fisher + fisher.T)
    values, vectors = np.linalg.eigh(fisher)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]

    p = int(values.size)
    lambda_max = float(values[0]) if values.size else 0.0
    floor = float(p) * _EPS64 * abs(lambda_max)
    resolvable = values[values > floor]
    n_resolvable = int(resolvable.size)

    total = float(np.sum(np.abs(values)))
    if total > 0.0:
        share = np.abs(values) / total
        effective = 1.0 / float(np.sum(share * share))
    else:
        effective = 0.0

    decades: float | None = None
    if n_resolvable >= 1 and float(resolvable[-1]) > 0.0 and lambda_max > 0.0:
        decades = float(np.log10(lambda_max / float(resolvable[-1])))

    return FisherReport(
        parameter_names=sensitivity.parameter_names,
        fisher=fisher,
        eigenvalues=values,
        eigenvectors=vectors,
        resolvable_floor=floor,
        n_resolvable_directions=n_resolvable,
        effective_dimension=effective,
        sloppiness_decades=decades,
    )
