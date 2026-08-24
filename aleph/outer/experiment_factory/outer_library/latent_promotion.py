#!/usr/bin/env python3
"""Audit whether an external latent may enter the Aleph sweep pipeline.

This is a promotion contract, not a parameter mapper.  It deliberately keeps
externally inferred quantities as observable constraints until identifiability,
external validation, units, and an Aleph forward operator are all established.
"""

from __future__ import annotations

from typing import Any


SCHEMA = "aleph.outer_library.latent_promotion.v1"


def evaluate(
    inference: dict[str, Any],
    *,
    external_dataset_holdout_passed: bool,
    external_lab_holdout_passed: bool,
    uncertainty_calibrated: bool,
    operator_mapping: dict[str, Any] | None = None,
    forward_sensitivity_passed: bool = False,
) -> dict[str, Any]:
    """Return an explicit staged decision; never emit a numeric sweep range."""
    mapping = operator_mapping or {}
    ambiguity_sets = inference.get("ambiguity_sets", [])
    ambiguity_resolved = not ambiguity_sets or all(
        bool(item.get("resolved_by_independent_observable", False))
        for item in ambiguity_sets
    )
    mapping_complete = bool(
        mapping.get("operator_status") in {"native_validated", "surrogate_validated"}
        and mapping.get("operator_id")
        and mapping.get("aleph_observable")
        and mapping.get("external_unit")
        and mapping.get("aleph_unit")
        and mapping.get("unit_conversion_verified") is True
        and mapping.get("geometry_matched") is True
        and isinstance(mapping.get("validation_report_sha256"), str)
        and len(mapping["validation_report_sha256"]) == 64
    )
    checks = {
        "inference_not_refused": inference.get("refused") is False,
        "inference_has_no_Aleph_parameter_authority": (
            inference.get("may_emit_Aleph_parameter") is False
        ),
        "external_dataset_holdout_passed": bool(external_dataset_holdout_passed),
        "external_lab_holdout_passed": bool(external_lab_holdout_passed),
        "uncertainty_calibrated": bool(uncertainty_calibrated),
        "latent_identifiability_resolved": ambiguity_resolved,
        "Aleph_observation_operator_mapping_complete": mapping_complete,
        "Aleph_forward_sensitivity_passed": bool(forward_sensitivity_passed),
    }
    evidence_ready = all(checks[name] for name in (
        "inference_not_refused",
        "inference_has_no_Aleph_parameter_authority",
        "external_dataset_holdout_passed",
        "external_lab_holdout_passed",
        "uncertainty_calibrated",
    ))
    operator_ready = evidence_ready and ambiguity_resolved and mapping_complete
    sweep_ready = operator_ready and bool(forward_sensitivity_passed)
    if sweep_ready:
        stage = "eligible_for_bounded_sweep_proposal"
    elif operator_ready:
        stage = "blocked_pending_Aleph_forward_sensitivity"
    elif evidence_ready:
        stage = "external_evidence_only_blocked_by_identifiability_or_operator"
    else:
        stage = "external_training_or_validation_only"
    return {
        "schema": SCHEMA,
        "stage": stage,
        "checks": checks,
        "failed_checks": [name for name, passed in checks.items() if not passed],
        "external_evidence_ready": evidence_ready,
        "operator_constraint_ready": operator_ready,
        "bounded_sweep_proposal_eligible": sweep_ready,
        "numeric_parameter_range": None,
        "may_select_Aleph_parameter": False,
        "may_mutate_physics": False,
        "operator_mapping": mapping,
        "ambiguity_sets": ambiguity_sets,
        "next_stage": (
            "human_and_physics_review_of_bounded_sweep_proposal"
            if sweep_ready else
            "run_unit_matched_Aleph_forward_sensitivity"
            if operator_ready else
            "resolve_identifiability_and_define_geometry_matched_operator"
            if evidence_ready else
            "acquire_independent_dataset_and_lab_validation_with_calibrated_uncertainty"
        ),
    }
