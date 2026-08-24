"""Independent, domain-guarded biophysical theory experts.

This module has no Aleph import and emits no engine parameter.  Each function
uses explicit units in its name/docstring and refuses invalid domains rather
than extrapolating silently.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np


@dataclass(frozen=True, slots=True)
class TheoryRefusal(Exception):
    theory_id: str
    code: str
    reason: str
    detail: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "refused": True,
            "theory_id": self.theory_id,
            "code": self.code,
            "reason": self.reason,
            "detail": self.detail,
        }


def _finite(name: str, value: float, theory_id: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise TheoryRefusal(theory_id, "NONFINITE_INPUT", f"{name} must be finite", {name: value})
    return result


def tether_force(
    bending_rigidity_pn_um: float,
    apparent_tension_pn_per_um: float,
) -> dict[str, Any]:
    """Quasi-static long-tether force, F=2π√(2κσ_app), in pN."""
    theory_id = "membrane_tether_quasistatic"
    kappa = _finite("bending_rigidity_pn_um", bending_rigidity_pn_um, theory_id)
    sigma = _finite("apparent_tension_pn_per_um", apparent_tension_pn_per_um, theory_id)
    if kappa <= 0 or sigma < 0:
        raise TheoryRefusal(
            theory_id, "OUTSIDE_DOMAIN", "requires kappa>0 and sigma_app>=0",
            {"bending_rigidity_pn_um": kappa, "apparent_tension_pn_per_um": sigma},
        )
    force = 2.0 * math.pi * math.sqrt(2.0 * kappa * sigma)
    return {
        "theory_id": theory_id,
        "tether_force_pn": force,
        "identifiable_combination": "bending_rigidity_times_apparent_tension",
        "aleph_dependency": False,
    }


def tether_apparent_tension(
    tether_force_pn: float,
    bending_rigidity_pn_um: float,
) -> dict[str, Any]:
    """Invert long-tether force for σ_app with fixed κ."""
    theory_id = "membrane_tether_quasistatic"
    force = _finite("tether_force_pn", tether_force_pn, theory_id)
    kappa = _finite("bending_rigidity_pn_um", bending_rigidity_pn_um, theory_id)
    if force < 0 or kappa <= 0:
        raise TheoryRefusal(
            theory_id, "OUTSIDE_DOMAIN", "requires force>=0 and kappa>0",
            {"tether_force_pn": force, "bending_rigidity_pn_um": kappa},
        )
    sigma = force * force / (8.0 * math.pi * math.pi * kappa)
    return {
        "theory_id": theory_id,
        "apparent_tension_pn_per_um": sigma,
        "identifiable": ["apparent_tension_given_bending_rigidity"],
        "non_identifiable": [
            "bilayer_tension_vs_membrane_cortex_attachment",
            "bending_rigidity_vs_apparent_tension_if_kappa_is_not_fixed",
        ],
        "aleph_dependency": False,
    }


def tether_ambiguity_set(
    apparent_tension_pn_per_um: float,
    *,
    points: int = 101,
) -> dict[str, Any]:
    """Return exact σ_bilayer+W decompositions, never a fabricated point estimate."""
    theory_id = "membrane_tether_quasistatic"
    sigma_app = _finite("apparent_tension_pn_per_um", apparent_tension_pn_per_um, theory_id)
    if sigma_app < 0 or points < 2:
        raise TheoryRefusal(
            theory_id, "OUTSIDE_DOMAIN", "requires sigma_app>=0 and at least two ambiguity points",
            {"apparent_tension_pn_per_um": sigma_app, "points": points},
        )
    bilayer = np.linspace(0.0, sigma_app, points, dtype=np.float64)
    attachment = sigma_app - bilayer
    return {
        "theory_id": theory_id,
        "relation": "sigma_app = sigma_bilayer + membrane_cortex_attachment_energy_density",
        "status": "exact_non_identifiability",
        "point_count": points,
        "sigma_bilayer_pn_per_um": bilayer,
        "attachment_energy_density_pn_per_um": attachment,
        "may_choose_one_decomposition": False,
        "aleph_dependency": False,
    }


def standard_linear_solid_relaxation(
    time_s: np.ndarray | float,
    initial_response: float,
    long_time_response: float,
    relaxation_time_s: float,
) -> dict[str, Any]:
    theory_id = "standard_linear_solid_single_mode"
    time = np.asarray(time_s, dtype=np.float64)
    y0 = _finite("initial_response", initial_response, theory_id)
    y_inf = _finite("long_time_response", long_time_response, theory_id)
    tau = _finite("relaxation_time_s", relaxation_time_s, theory_id)
    if tau <= 0 or not np.isfinite(time).all() or (time < 0).any():
        raise TheoryRefusal(
            theory_id, "OUTSIDE_DOMAIN", "requires finite time>=0 and tau>0",
            {"relaxation_time_s": tau, "minimum_time_s": float(time.min(initial=0.0))},
        )
    response = y_inf + (y0 - y_inf) * np.exp(-time / tau)
    return {
        "theory_id": theory_id,
        "response": response,
        "response_unit": "same_as_initial_and_long_time_response",
        "aleph_dependency": False,
    }


def persistent_random_walk_msd_2d(
    lag_time_min: np.ndarray | float,
    speed_um_per_min: float,
    persistence_time_min: float,
) -> dict[str, Any]:
    theory_id = "persistent_random_walk_2d"
    lag = np.asarray(lag_time_min, dtype=np.float64)
    speed = _finite("speed_um_per_min", speed_um_per_min, theory_id)
    tau = _finite("persistence_time_min", persistence_time_min, theory_id)
    if speed < 0 or tau <= 0 or not np.isfinite(lag).all() or (lag < 0).any():
        raise TheoryRefusal(
            theory_id, "OUTSIDE_DOMAIN", "requires lag>=0, speed>=0, and persistence_time>0",
            {"speed_um_per_min": speed, "persistence_time_min": tau},
        )
    msd = 2.0 * speed * speed * tau * (lag - tau * (-np.expm1(-lag / tau)))
    return {
        "theory_id": theory_id,
        "mean_squared_displacement_um2": msd,
        "geometry": "isotropic_2D",
        "aleph_dependency": False,
    }


def evaluate_expert(theory_id: str, inputs: dict[str, Any]) -> dict[str, Any]:
    """Typed dispatcher suitable for a future gating network."""
    try:
        if theory_id == "membrane_tether_quasistatic":
            return tether_force(**inputs)
        if theory_id == "standard_linear_solid_single_mode":
            return standard_linear_solid_relaxation(**inputs)
        if theory_id == "persistent_random_walk_2d":
            return persistent_random_walk_msd_2d(**inputs)
        raise TheoryRefusal(
            theory_id, "EXPERT_NOT_IMPLEMENTED", "theory is not an executable expert",
            {"available": [
                "membrane_tether_quasistatic",
                "standard_linear_solid_single_mode",
                "persistent_random_walk_2d",
            ]},
        )
    except TheoryRefusal as refusal:
        return refusal.as_dict()
