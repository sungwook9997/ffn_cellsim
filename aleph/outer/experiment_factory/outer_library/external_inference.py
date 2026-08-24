"""Independent biophysical inference facade over applicability-gated experts."""

from __future__ import annotations

from typing import Any
import math

from theory_runtime import TheoryRefusal, tether_apparent_tension


SCHEMA = "outer.biophysics.inference.v1"


def _base(provenance: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "refused": False,
        "input_provenance": provenance,
        "state_posterior": [],
        "observable_predictions": [],
        "mechanism_support": [],
        "intervention_directions": [],
        "ambiguity_sets": [],
        "uncertainty": {},
        "applicability": {},
        "refusal": None,
        "aleph_dependency": False,
        "may_emit_Aleph_parameter": False,
    }


def _refuse(provenance: dict[str, Any], code: str, reason: str, detail: dict[str, Any]) -> dict[str, Any]:
    result = _base(provenance)
    result["refused"] = True
    result["refusal"] = {"code": code, "reason": reason, "detail": detail}
    result["applicability"] = {"status": "outside_current_expert_scope"}
    return result


def infer_tether_observation(
    *,
    force_pn: float,
    provenance: dict[str, Any],
    bending_rigidity_pn_um: float | None = None,
    force_standard_error_pn: float | None = None,
) -> dict[str, Any]:
    """Infer only identifiable tether combinations; never choose σ versus W."""
    result = _base(provenance)
    force = float(force_pn)
    if not math.isfinite(force) or force < 0:
        return _refuse(provenance, "INVALID_FORCE", "tether force must be finite and non-negative", {"force_pn": force_pn})
    product = force * force / (8.0 * math.pi * math.pi)
    result["state_posterior"].append({
        "latent_id": "bending_rigidity_times_apparent_tension",
        "estimate": product,
        "unit": "pN^2",
        "conditioning": "none_beyond_quasistatic_tether_expert",
        "authority": "theory_conditional",
    })
    if bending_rigidity_pn_um is None:
        result["ambiguity_sets"].append({
            "relation": "kappa * sigma_app = observed_force^2/(8*pi^2)",
            "constant_pn2": product,
            "free_variables": ["bending_rigidity_pn_um", "apparent_tension_pn_per_um"],
            "may_choose_one_point": False,
        })
    else:
        try:
            apparent = tether_apparent_tension(force, bending_rigidity_pn_um)
        except TheoryRefusal as refusal:
            return _refuse(provenance, refusal.code, refusal.reason, refusal.detail)
        sigma_app = apparent["apparent_tension_pn_per_um"]
        result["state_posterior"].append({
            "latent_id": "apparent_membrane_tension",
            "estimate": sigma_app,
            "unit": "pN/micrometre",
            "conditioning": f"bending_rigidity_fixed_at_{bending_rigidity_pn_um}_pN_um",
            "authority": "theory_conditional",
        })
        result["ambiguity_sets"].append({
            "relation": "sigma_app = sigma_bilayer + membrane_cortex_attachment_energy_density",
            "constant_pn_per_um": sigma_app,
            "bounds_under_nonnegative_components": {
                "sigma_bilayer_pn_per_um": [0.0, sigma_app],
                "attachment_energy_density_pn_per_um": [0.0, sigma_app],
            },
            "may_choose_one_point": False,
        })
    if force_standard_error_pn is not None:
        se = float(force_standard_error_pn)
        if not math.isfinite(se) or se < 0:
            return _refuse(provenance, "INVALID_UNCERTAINTY", "force standard error must be finite and non-negative", {"force_standard_error_pn": force_standard_error_pn})
        result["uncertainty"] = {
            "kind": "first_order_measurement_propagation",
            "force_standard_error_pn": se,
            "product_standard_error_pn2": abs(force) * se / (4.0 * math.pi * math.pi),
            "does_not_include_model_discrepancy": True,
        }
    else:
        result["uncertainty"] = {
            "kind": "missing_measurement_uncertainty",
            "calibrated": False,
        }
    result["mechanism_support"] = [{
        "theory_id": "membrane_tether_quasistatic",
        "support": "compatible_not_compared_against_alternative_models",
        "probability": None,
    }]
    result["applicability"] = {
        "status": "conditional",
        "required_assumptions": [
            "long_cylindrical_tether", "quasi_static_plateau",
            "homogeneous_effective_kappa_and_sigma_app",
        ],
        "assumptions_verified_by_input": False,
    }
    return result


def infer(observation: dict[str, Any]) -> dict[str, Any]:
    modality = observation.get("modality")
    provenance = dict(observation.get("provenance", {}))
    if modality == "AFM_membrane_tether_force":
        if observation.get("unit") != "pN":
            return _refuse(provenance, "UNIT_MISMATCH", "tether expert requires pN force", {"unit": observation.get("unit")})
        return infer_tether_observation(
            force_pn=observation.get("value"),
            bending_rigidity_pn_um=observation.get("bending_rigidity_pn_um"),
            force_standard_error_pn=observation.get("standard_error"),
            provenance=provenance,
        )
    return _refuse(
        provenance, "NO_APPLICABLE_EXPERT", "no implemented expert accepts this measurement operator",
        {"modality": modality, "implemented_modalities": ["AFM_membrane_tether_force"]},
    )
