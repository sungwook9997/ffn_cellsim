#!/usr/bin/env python3
"""Exercise the independent inference contract on measured and OOD inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from external_inference import infer


def evaluate() -> dict[str, object]:
    provenance = {
        "dataset_id": "elife-72381-hela-mechano-osmotic",
        "measurement_operator": "AFM_membrane_tether_plateau_force",
        "authority": "author_reported_experiment_mean",
    }
    without_kappa = infer({
        "modality": "AFM_membrane_tether_force",
        "value": 18.6974,
        "unit": "pN",
        "standard_error": 0.30,
        "provenance": provenance,
    })
    conditional_kappa = infer({
        "modality": "AFM_membrane_tether_force",
        "value": 18.6974,
        "unit": "pN",
        "standard_error": 0.30,
        "bending_rigidity_pn_um": 0.082,
        "provenance": provenance,
    })
    wrong_unit = infer({
        "modality": "AFM_membrane_tether_force",
        "value": 18.6974,
        "unit": "nN",
        "provenance": provenance,
    })
    unsupported = infer({
        "modality": "western_blot_band_intensity",
        "value": 1.2,
        "unit": "relative_intensity",
        "provenance": {"dataset_id": "synthetic_test_only"},
    })
    checks = {
        "force_alone_returns_only_kappa_sigma_product": (
            [item["latent_id"] for item in without_kappa["state_posterior"]]
            == ["bending_rigidity_times_apparent_tension"]
        ),
        "fixed_kappa_adds_conditional_apparent_tension": (
            {item["latent_id"] for item in conditional_kappa["state_posterior"]}
            == {"bending_rigidity_times_apparent_tension", "apparent_membrane_tension"}
        ),
        "bilayer_vs_attachment_remains_ambiguous": (
            conditional_kappa["ambiguity_sets"][-1]["may_choose_one_point"] is False
        ),
        "mechanism_probability_not_fabricated": (
            conditional_kappa["mechanism_support"][0]["probability"] is None
        ),
        "wrong_unit_refused": wrong_unit["refused"] and wrong_unit["refusal"]["code"] == "UNIT_MISMATCH",
        "unsupported_operator_refused": unsupported["refused"] and unsupported["refusal"]["code"] == "NO_APPLICABLE_EXPERT",
        "no_Aleph_dependency_or_parameter_output": all(
            item["aleph_dependency"] is False and item["may_emit_Aleph_parameter"] is False
            for item in (without_kappa, conditional_kappa, wrong_unit, unsupported)
        ),
    }
    return {
        "schema": "outer.biophysics.external_inference_evaluation.v1",
        "phase_objective": "train_and_validate_external_observation_inference_while_Aleph_sweeps_are_infeasible",
        "final_program_objective": "promote_validated_observables_into_bounded_Aleph_parameter_sweeps",
        "objective_change": False,
        "execution_order_change": True,
        "checks": checks,
        "all_checks_passed": all(checks.values()),
        "measured_tether_without_kappa": without_kappa,
        "measured_tether_with_assumed_kappa": conditional_kappa,
        "refusal_examples": {"wrong_unit": wrong_unit, "unsupported_operator": unsupported},
        "current_external_training_phase_depends_on_Aleph_sweep": False,
        "final_program_objective_includes_Aleph_sweep": True,
        "required_promotion_contract": "aleph.outer_library.latent_promotion.v1",
        "aleph_dependency": False,
        "may_emit_Aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate()
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["checks"], sort_keys=True))


if __name__ == "__main__":
    main()
