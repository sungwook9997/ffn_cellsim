"""Fail-closed probability-mass checks for analytic and archived distribution artifacts."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = ["ProbabilityMassAudit", "audit_probability_mass"]


@dataclass(frozen=True, slots=True)
class ProbabilityMassAudit:
    """Measured mass diagnostics; input is never clipped or renormalized."""

    total_mass: float
    mass_error: float
    minimum_value: float
    maximum_value: float
    normalization_tolerance: float
    negativity_tolerance: float
    measure_kind: str = "discrete-probability-mass"


def audit_probability_mass(
    values: Any,
    *,
    normalization_tolerance: float = 1e-12,
    negativity_tolerance: float = 0.0,
) -> ProbabilityMassAudit:
    """Validate a discrete probability mass without silently repairing it."""
    mass = np.asarray(values, dtype=np.float64)
    if mass.size == 0:
        raise ValueError("probability mass must not be empty")
    if not np.all(np.isfinite(mass)):
        raise ValueError("probability mass must be finite")
    for name, value in (
        ("normalization_tolerance", normalization_tolerance),
        ("negativity_tolerance", negativity_tolerance),
    ):
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0.0
        ):
            raise ValueError(f"{name} must be finite and nonnegative")
    minimum = float(np.min(mass))
    if minimum < -float(negativity_tolerance):
        raise ValueError("probability mass contains a negative value")
    total = float(np.sum(mass, dtype=np.float64))
    error = abs(total - 1.0)
    if error > float(normalization_tolerance):
        raise ValueError(
            f"probability mass is not normalized: total={total:.17g}, error={error:.3g}"
        )
    return ProbabilityMassAudit(
        total_mass=total,
        mass_error=error,
        minimum_value=minimum,
        maximum_value=float(np.max(mass)),
        normalization_tolerance=float(normalization_tolerance),
        negativity_tolerance=float(negativity_tolerance),
    )
