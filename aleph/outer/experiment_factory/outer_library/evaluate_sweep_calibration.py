#!/usr/bin/env python3
"""Validate an Aleph forward-sensitivity bundle before emitting a numeric axis."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


SCHEMA = "aleph.outer_library.forward_sensitivity.v1"
ALLOWED_AXES = {
    "aleph.vertical.sf_arc.MaterialCard.active_tension_pn": "pN",
    "aleph.vertical.nmii.MotorKinetics.crossbridge_stiffness": "pN_per_um",
    "aleph.vertical.nmii.MotorKinetics.attachment_rate": "per_s",
    "aleph.vertical.focal_adhesion.ClutchCard.bond_stiffness_pn_per_um": "pN_per_um",
    "aleph.vertical.focal_adhesion.ClutchCard.anchor_stiffness_pn_per_um": "pN_per_um",
}
OBSERVABLE_UNITS = {
    "tissue_traction_force_per_cell": "nN_per_cell",
    "membrane_tether_plateau_force": "pN",
}
TRACTION_OBSERVABLE = "tissue_traction_force_per_cell"
TETHER_OBSERVABLE = "membrane_tether_plateau_force"


def _typed_refusal(calibration: dict) -> dict | None:
    refusal = calibration.get("refusal")
    if refusal is None:
        return None
    required = ("code", "reason", "operator_id", "detail")
    if refusal.get("refused") is not True or any(
        not refusal.get(field) for field in required
    ):
        return {
            "status": "refused",
            "refusal_is_well_formed": False,
            "failed_checks": ["refusal_is_well_formed"],
            "may_emit_sweep_axis": False,
            "may_mutate_physics": False,
            "may_write_decision_authority": False,
            "authority": "none",
        }
    return {
        "status": "refused",
        "refusal_is_well_formed": True,
        "upstream_refusal": refusal,
        "failed_checks": [f"operator_refused:{refusal['code']}"],
        "may_emit_sweep_axis": False,
        "may_mutate_physics": False,
        "may_write_decision_authority": False,
        "authority": "none",
    }


def _external_target(report: dict, observable: str) -> tuple[list[float], str]:
    if observable == TRACTION_OBSERVABLE:
        return (
            report["matrix_tension_contrasts"][observable]["bootstrap_ci95"],
            "higher matrix-tension traction contrast",
        )
    if observable == TETHER_OBSERVABLE:
        return (
            report["AFM_membrane_tether_force"]["spreading_phase_30_to_90_min"]
            ["paired_experiment_bootstrap_95ci_pN"],
            "Y-27632 minus control early tether-force contrast",
        )
    raise ValueError(f"unsupported observable {observable!r}")


def _interpolate_axis(
    axis_values: np.ndarray,
    effects: np.ndarray,
    target: tuple[float, float],
) -> list[float]:
    # np.interp accepts non-decreasing x. Duplicate response values are kept at
    # the last corresponding axis value, which is conservative at a plateau.
    return [float(np.interp(value, effects, axis_values)) for value in target]


def evaluate(calibration: dict, external_report: dict) -> dict:
    refusal_result = _typed_refusal(calibration)
    if refusal_result is not None:
        return refusal_result
    axis = calibration.get("axis", "")
    observable = calibration.get("observable", "")
    runs = calibration.get("runs", [])
    checks = {
        "schema_is_supported": calibration.get("schema") == SCHEMA,
        "engine_is_project_aleph": calibration.get("engine") == "Project_Aleph",
        "physics_commit_recorded": bool(str(calibration.get("physics_commit", "")).strip()),
        "scenario_recorded": bool(str(calibration.get("scenario", "")).strip()),
        "observable_operator_recorded": bool(
            str(calibration.get("observable_operator", "")).strip()
        ),
        "observable_operator_is_validated": calibration.get(
            "observable_operator_status"
        ) in {"native_validated", "surrogate_validated"},
        "observable_operator_validation_is_pinned": (
            isinstance(calibration.get("observable_operator_validation_sha256"), str)
            and len(calibration["observable_operator_validation_sha256"]) == 64
            and all(character in "0123456789abcdef" for character in calibration[
                "observable_operator_validation_sha256"
            ].lower())
        ),
        "axis_is_declared": axis in ALLOWED_AXES,
        "axis_unit_matches_declaration": (
            axis in ALLOWED_AXES
            and calibration.get("axis_unit") == ALLOWED_AXES[axis]
        ),
        "observable_is_declared": observable in OBSERVABLE_UNITS,
        "observable_unit_matches_experiment": (
            observable in OBSERVABLE_UNITS
            and calibration.get("observable_unit") == OBSERVABLE_UNITS[observable]
        ),
        # No existing Aleph parameter identifies membrane-cortex attachment W.
        # Tether force constrains sigma_bilayer + W, so accepting any current
        # one-axis card here would silently fit a missing term with the wrong
        # parameter.
        "observable_axis_pair_is_identifiable": observable == TRACTION_OBSERVABLE,
        "reference_axis_value_recorded": isinstance(
            calibration.get("reference_axis_value"), (int, float)
        ),
        "run_ids_are_unique": (
            len(runs) > 0
            and len({str(run.get("run_id", "")) for run in runs}) == len(runs)
            and all(str(run.get("run_id", "")).strip() for run in runs)
        ),
        "all_runs_converged": bool(runs) and all(run.get("converged") is True for run in runs),
        "all_runs_finite": bool(runs) and all(
            isinstance(run.get("axis_value"), (int, float))
            and isinstance(run.get("observable_value"), (int, float))
            and np.isfinite(run["axis_value"])
            and np.isfinite(run["observable_value"])
            for run in runs
        ),
    }

    grouped: dict[float, list[float]] = defaultdict(list)
    if checks["all_runs_finite"]:
        for run in runs:
            grouped[float(run["axis_value"])].append(float(run["observable_value"]))
    axis_values = np.asarray(sorted(grouped), dtype=np.float64)
    means = np.asarray([
        np.mean(grouped[value]) for value in axis_values
    ], dtype=np.float64)
    replicate_counts = {
        str(value): len(grouped[value]) for value in axis_values
    }
    checks["at_least_four_axis_levels"] = len(axis_values) >= 4
    checks["at_least_three_repeats_per_level"] = (
        bool(grouped) and min(map(len, grouped.values())) >= 3
    )
    reference = calibration.get("reference_axis_value")
    checks["reference_level_is_present"] = (
        isinstance(reference, (int, float)) and float(reference) in grouped
    )

    effects = np.asarray([], dtype=np.float64)
    monotonic = False
    brackets = False
    numeric_range = None
    target_ci, target_description = _external_target(external_report, observable)
    prerequisites_pass = all(checks.values())
    if prerequisites_pass and len(axis_values):
        reference_mean = float(np.mean(grouped[float(reference)]))
        effects = means - reference_mean
        direction = calibration.get("expected_direction")
        signed = effects if direction == "increases" else -effects
        monotonic = direction in {"increases", "decreases"} and bool(
            np.all(np.diff(signed) >= -1e-12) and np.any(np.diff(signed) > 1e-12)
        )
        if monotonic:
            target = np.asarray(target_ci, dtype=np.float64)
            brackets = bool(signed.min() <= target[0] and signed.max() >= target[1])
            if brackets:
                numeric_range = _interpolate_axis(
                    axis_values, signed, (float(target[0]), float(target[1]))
                )
    checks["response_is_monotonic_in_declared_direction"] = monotonic
    checks["simulated_response_brackets_external_ci"] = brackets
    passed = all(checks.values())
    return {
        "status": (
            "eligible_for_non_authoritative_bounded_axis" if passed else "refused"
        ),
        "axis": axis,
        "axis_unit": calibration.get("axis_unit"),
        "checks": checks,
        "failed_checks": [name for name, value in checks.items() if not value],
        "axis_levels": axis_values.tolist(),
        "replicate_counts": replicate_counts,
        "mean_observable_by_level": means.tolist(),
        "effect_from_reference_by_level": effects.tolist(),
        "external_target_ci95": target_ci,
        "external_target_description": target_description,
        "numeric_axis_range": numeric_range,
        "numeric_axis_range_is_interpolated": numeric_range is not None,
        "within_simulated_hull": brackets,
        "may_emit_sweep_axis": passed,
        "may_mutate_physics": False,
        "may_write_decision_authority": False,
        "authority": "non_authoritative_human_review_required",
        "limitations": [
            "one-axis sensitivity does not identify confounded motor and adhesion parameters",
            "interpolation is valid only inside the simulated response hull",
            "the proposal inherits the scenario and observable operator recorded in the bundle",
            "tether force identifies apparent tension sigma_bilayer + W, not either term alone",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--tendon-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(
        json.loads(args.calibration.read_text(encoding="utf-8")),
        json.loads(args.tendon_report.read_text(encoding="utf-8")),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
