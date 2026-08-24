#!/usr/bin/env python3
"""Evaluate live morphology -> terminal Annexin-V on sealed S-BIAD2515 wells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr


def _fit_ridge(
    X: np.ndarray, y: np.ndarray, alpha: float
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    mean = X.mean(axis=0)
    scale = X.std(axis=0)
    scale[scale < 1e-8] = 1.0
    standardized = (X - mean) / scale
    target_mean = float(y.mean())
    # The dual solution is exact and avoids a 2,336 x 2,336 inverse for ten wells.
    dual = np.linalg.solve(
        standardized @ standardized.T + alpha * np.eye(len(X)),
        y - target_mean,
    )
    coefficients = standardized.T @ dual
    return mean, scale, target_mean, coefficients


def _predict(
    model: tuple[np.ndarray, np.ndarray, float, np.ndarray], X: np.ndarray
) -> np.ndarray:
    mean, scale, target_mean, coefficients = model
    return (X - mean) / scale @ coefficients + target_mean


def _regression_metrics(
    truth: np.ndarray, predicted: np.ndarray
) -> dict[str, float | int | None]:
    residual = predicted - truth
    total = np.sum((truth - truth.mean()) ** 2)
    rho = (
        float(spearmanr(truth, predicted).statistic)
        if np.std(predicted) > 1e-15
        else None
    )
    return {
        "well_count": int(len(truth)),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual ** 2))),
        "r2": float(1.0 - np.sum(residual ** 2) / total),
        "spearman_rho": rho,
    }


def _binary_auc(truth: np.ndarray, score: np.ndarray) -> float:
    positive = truth.astype(bool)
    ranks = rankdata(score, method="average")
    return float(
        (ranks[positive].sum() - positive.sum() * (positive.sum() + 1) / 2)
        / (positive.sum() * (~positive).sum())
    )


def evaluate(tensor_path: Path, checkpoint_path: Path) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        X = np.asarray(data["X"], dtype=np.float64)
        y = np.asarray(data["y"], dtype=np.float64)
        well = np.asarray(data["well"]).astype(str)
        dose = np.asarray(data["dose_nM"], dtype=np.float64)
        apoptosis = np.asarray(data["author_apoptosis_ground_truth"]).astype(str)
        split = np.asarray(data["split"], dtype=np.uint8)
        feature_names = np.asarray(data["feature_name"]).astype(str)
        target_name = str(np.asarray(data["target_name"]))
    if X.shape != (30, 2_336) or len(feature_names) != X.shape[1]:
        raise ValueError("S-BIAD2515 tensor shape changed")
    if any(name.startswith("Terminal_") for name in feature_names):
        raise ValueError("terminal Annexin-V feature leaked into live input")
    masks = {"fit": split == 0, "calibration": split == 2, "test": split == 3}
    if any(mask.sum() != 10 for mask in masks.values()):
        raise ValueError("S-BIAD2515 requires ten dose-balanced wells per split")
    if any(len(set(well[left]) & set(well[right])) for left in masks.values() for right in masks.values() if left is not right):
        raise ValueError("S-BIAD2515 well leakage detected")

    alpha_grid = np.logspace(-4, 6, 41)
    fit_indices = np.flatnonzero(masks["fit"])
    cross_validation_mse = []
    for alpha in alpha_grid:
        predictions = []
        for held_out in range(len(fit_indices)):
            training = np.delete(fit_indices, held_out)
            model = _fit_ridge(X[training], y[training], float(alpha))
            predictions.append(
                float(_predict(model, X[fit_indices[held_out:held_out + 1]])[0])
            )
        cross_validation_mse.append(
            float(np.mean((np.asarray(predictions) - y[fit_indices]) ** 2))
        )
    selected_alpha = float(alpha_grid[int(np.argmin(cross_validation_mse))])
    model = _fit_ridge(X[masks["fit"]], y[masks["fit"]], selected_alpha)
    calibration_prediction = _predict(model, X[masks["calibration"]])
    test_prediction = _predict(model, X[masks["test"]])

    calibration_residual = np.abs(y[masks["calibration"]] - calibration_prediction)
    rank = min(
        len(calibration_residual) - 1,
        int(np.ceil((len(calibration_residual) + 1) * 0.90)) - 1,
    )
    interval_radius = float(np.partition(calibration_residual, rank)[rank])
    test_residual = np.abs(y[masks["test"]] - test_prediction)

    # A dose-only baseline checks whether image morphology adds information
    # beyond the known staurosporine concentration.
    fit_dose_design = np.column_stack(
        (np.ones(masks["fit"].sum()), np.log1p(dose[masks["fit"]]))
    )
    dose_coefficients = np.linalg.lstsq(
        fit_dose_design, y[masks["fit"]], rcond=None
    )[0]
    test_dose_prediction = np.column_stack(
        (np.ones(masks["test"].sum()), np.log1p(dose[masks["test"]]))
    ) @ dose_coefficients
    mean_prediction = np.full(masks["test"].sum(), y[masks["fit"]].mean())

    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        checkpoint_path,
        input_mean=model[0],
        input_scale=model[1],
        target_mean=np.asarray(model[2]),
        coefficients=model[3],
        feature_name=feature_names,
        target_name=np.asarray(target_name),
        alpha=np.asarray(selected_alpha),
        conformal_interval_radius=np.asarray(interval_radius),
    )
    test_truth = y[masks["test"]]
    test_metrics = _regression_metrics(test_truth, test_prediction)
    dose_metrics = _regression_metrics(test_truth, test_dose_prediction)
    mean_metrics = _regression_metrics(test_truth, mean_prediction)
    positive = apoptosis[masks["test"]] == "positive"
    return {
        "schema": "aleph.outer_library.biad2515_apoptosis_evaluation.v1",
        "task": "six_hour_live_multichannel_morphology_to_terminal_AnnexinV_observable",
        "target": target_name,
        "model": {
            "family": "dual_ridge_regression",
            "input_feature_count": X.shape[1],
            "fit_well_count": int(masks["fit"].sum()),
            "alpha_selected_by_fit_only_leave_one_well_out_cv": selected_alpha,
            "minimum_fit_only_cv_mse": float(min(cross_validation_mse)),
        },
        "split": {
            "unit": "dose_matched_independent_well",
            "fit_wells": well[masks["fit"]].tolist(),
            "calibration_wells": well[masks["calibration"]].tolist(),
            "sealed_test_wells": well[masks["test"]].tolist(),
            "well_overlap_count": 0,
            "each_of_ten_doses_present_once_per_split": True,
            "test_used_for_fitting_selection_or_calibration": False,
        },
        "sealed_same_provider_test": {
            **test_metrics,
            "author_positive_vs_other_auroc_from_predicted_AnnexinV": _binary_auc(
                positive, test_prediction
            ),
        },
        "baselines_on_same_sealed_test": {
            "dose_only_log_linear": dose_metrics,
            "fit_mean": mean_metrics,
            "image_model_beats_dose_only_mae": test_metrics["mae"] < dose_metrics["mae"],
        },
        "conformal_90": {
            "unit": "well",
            "calibration_well_count": int(masks["calibration"].sum()),
            "interval_radius_author_normalized_units": interval_radius,
            "sealed_test_coverage": float((test_residual <= interval_radius).mean()),
            "mean_interval_width_author_normalized_units": 2.0 * interval_radius,
            "conditional_coverage_authority": False,
        },
        "raw_image_accession": "S-BIAD2515",
        "author_processed_profile_source": True,
        "terminal_features_excluded_from_inputs": True,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "ood_validated": False,
        "status": "same_provider_well_holdout_training_only",
        "production_eligible": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "limitations": [
            "only three wells per dose and one provider are available",
            "the input is an author-processed 2,336-feature profile rather than raw pixels",
            "the ten sealed test wells establish well transfer, not dataset or laboratory transfer",
            "the 90% conformal interval has only ten calibration wells and no conditional authority",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor, args.checkpoint)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report["sealed_same_provider_test"], sort_keys=True))


if __name__ == "__main__":
    main()
