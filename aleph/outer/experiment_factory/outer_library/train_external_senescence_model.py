#!/usr/bin/env python3
"""Train and externally test a fibroblast senescence marker head."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def _sigmoid(value: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(value, -30.0, 30.0)))


def _representation(
    X: np.ndarray, fit: np.ndarray, mode: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if mode == "logCPM_fit_standardized":
        base = X.astype(np.float64)
    elif mode == "within_sample_rank_fit_standardized":
        order = np.argsort(np.argsort(X, axis=1, kind="stable"), axis=1, kind="stable")
        base = order.astype(np.float64) / max(1, X.shape[1] - 1)
    else:
        raise ValueError(f"unknown representation {mode}")
    mean = base[fit].mean(axis=0)
    scale = base[fit].std(axis=0)
    scale[scale < 1e-6] = 1.0
    return (base - mean) / scale, mean, scale


def _metrics(probability: np.ndarray, labels: np.ndarray) -> dict[str, object]:
    prediction = probability >= 0.5
    tp = int(np.sum(prediction & (labels == 1)))
    tn = int(np.sum(~prediction & (labels == 0)))
    fp = int(np.sum(prediction & (labels == 0)))
    fn = int(np.sum(~prediction & (labels == 1)))
    return {
        "sample_count": int(len(labels)),
        "accuracy": float((prediction == labels).mean()),
        "balanced_accuracy": float((tp / max(1, tp + fn) + tn / max(1, tn + fp)) / 2),
        "f1": float(2 * tp / max(1, 2 * tp + fp + fn)),
        "brier": float(np.mean((probability - labels) ** 2)),
        "confusion_matrix_true_by_predicted": [[tn, fp], [fn, tp]],
    }


def _fit_logistic(
    X: np.ndarray, labels: np.ndarray, fit: np.ndarray, l2: float
) -> tuple[np.ndarray, float]:
    weights = np.zeros(X.shape[1], dtype=np.float64)
    bias = 0.0
    for _ in range(4_000):
        probability = _sigmoid(X[fit] @ weights + bias)
        error = probability - labels[fit]
        weights -= 0.03 * (X[fit].T @ error / fit.sum() + l2 * weights)
        bias -= 0.03 * float(error.mean())
    return weights, bias


def _temperature(logits: np.ndarray, labels: np.ndarray) -> float:
    values = np.linspace(0.35, 4.0, 147)
    losses = []
    for value in values:
        probability = _sigmoid(logits / value)
        losses.append(float(-np.mean(
            labels * np.log(probability + 1e-12)
            + (1 - labels) * np.log(1 - probability + 1e-12)
        )))
    return float(values[int(np.argmin(losses))])


def _quantile_higher(values: np.ndarray, coverage: float) -> float:
    rank = min(len(values) - 1, int(np.ceil((len(values) + 1) * coverage)) - 1)
    return float(np.partition(values, rank)[rank])


def _paired_report(
    probability: np.ndarray, labels: np.ndarray, replicate: np.ndarray
) -> dict[str, object]:
    differences = []
    for value in sorted(set(replicate.tolist())):
        pair = replicate == value
        if pair.sum() != 2 or set(labels[pair].tolist()) != {0, 1}:
            raise ValueError(f"external replicate {value} is not a state pair")
        differences.append(float(probability[pair & (labels == 1)][0] - probability[pair & (labels == 0)][0]))
    positive = sum(value > 0 for value in differences)
    n = len(differences)
    p_value = sum(math.comb(n, k) for k in range(positive, n + 1)) / (2 ** n)
    rng = np.random.default_rng(297406)
    bootstrap = np.asarray([
        np.mean(rng.choice(differences, size=n, replace=True)) for _ in range(10_000)
    ])
    return {
        "pair_count": n,
        "senescent_minus_proliferating_probability": differences,
        "positive_pair_count": positive,
        "mean_difference": float(np.mean(differences)),
        "bootstrap_95_interval": [float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))],
        "exact_one_sided_sign_test_p": float(p_value),
        "confirmatory_p_below_0_05": bool(p_value < 0.05),
    }


def train(tensor_path: Path, checkpoint_path: Path) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        X_raw = np.asarray(data["X"], dtype=np.float64)
        feature_name = np.asarray(data["feature_name"]).astype(str)
        dataset = np.asarray(data["dataset"]).astype(str)
        lab = np.asarray(data["lab"]).astype(str)
        role = np.asarray(data["role"]).astype(str)
        labels = np.asarray(data["state_label"], dtype=np.int8)
        replicate = np.asarray(data["replicate"], dtype=np.int64)
        geometry = np.asarray(data["geometry"]).astype(str)
        sample_id = np.asarray(data["sample_id"]).astype(str)
    if X_raw.shape != (42, 35):
        raise ValueError("external senescence tensor contract changed")
    fit = role == "fit"
    selection = role == "selection"
    calibration = role == "calibration"
    test = role == "sealed_external_lab_test"
    diagnostic = role == "single_donor_context_shift_diagnostic"
    if set(dataset[test]) != {"GSE297406"} or len(set(lab[test])) != 1:
        raise ValueError("sealed external senescence test changed")

    best = None
    for mode in ("logCPM_fit_standardized", "within_sample_rank_fit_standardized"):
        X, mean, scale = _representation(X_raw, fit, mode)
        for l2 in (1e-4, 1e-3, 1e-2, 1e-1, 1.0):
            weights, bias = _fit_logistic(X, labels, fit, l2)
            metrics = _metrics(_sigmoid(X[selection] @ weights + bias), labels[selection])
            score = float(metrics["balanced_accuracy"])
            candidate = (score, -l2, mode, l2, weights, bias, mean, scale, metrics)
            if best is None or candidate[:2] > best[:2]:
                best = candidate
    assert best is not None
    _, _, mode, l2, weights, bias, mean, scale, selection_metrics = best
    X, _, _ = _representation(X_raw, fit, mode)
    logits = X @ weights + bias
    temperature = _temperature(logits[calibration], labels[calibration])
    probability = _sigmoid(logits / temperature)
    calibration_scores = 1.0 - np.where(
        labels[calibration] == 1, probability[calibration], 1.0 - probability[calibration]
    )
    conformal_threshold = _quantile_higher(calibration_scores, 0.90)
    test_sets = np.column_stack([
        1.0 - probability[test] >= 1.0 - conformal_threshold,
        probability[test] >= 1.0 - conformal_threshold,
    ])
    test_labels = labels[test]
    coverage = float(test_sets[np.arange(test.sum()), test_labels].mean())
    test_metrics = _metrics(probability[test], test_labels)
    paired = _paired_report(probability[test], test_labels, replicate[test])

    diagnostic_results = {}
    for target_geometry in ("2D", "3D"):
        mask = diagnostic & (geometry == target_geometry)
        diagnostic_results[target_geometry] = {
            "early_proliferative_probability": float(probability[mask & (labels == 0)][0]),
            "deeply_senescent_probability": float(probability[mask & (labels == 1)][0]),
            "difference": float(probability[mask & (labels == 1)][0] - probability[mask & (labels == 0)][0]),
            "sample_ids": sample_id[mask].tolist(),
        }

    fit_distance = np.sqrt(np.mean(X[fit] ** 2, axis=1))
    calibration_distance = np.sqrt(np.mean(X[calibration] ** 2, axis=1))
    ood_threshold = _quantile_higher(calibration_distance, 0.95)
    external_distance = np.sqrt(np.mean(X[test] ** 2, axis=1))
    diagnostic_distance = np.sqrt(np.mean(X[diagnostic] ** 2, axis=1))

    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        checkpoint_path, feature_name=feature_name, representation=np.asarray(mode),
        input_mean=mean, input_scale=scale, weight=weights, bias=np.asarray(bias),
        temperature=np.asarray(temperature), conformal_threshold=np.asarray(conformal_threshold),
        ood_threshold=np.asarray(ood_threshold),
    )
    return {
        "schema": "aleph.outer_library.external_senescence_model.v1",
        "task": "fibroblast_senescence_marker_state_cross_lab",
        "model": {
            "family": "frozen_35_marker_logistic_regression",
            "representation_selected_on_GSE63577_IMR90_only": mode,
            "selected_l2": float(l2),
            "selection_metrics": selection_metrics,
            "temperature_selected_on_GSE63577_WI38_only": temperature,
        },
        "split": {
            "fit": "GSE63577_BJ_HFF_MRC5",
            "selection": "GSE63577_IMR90",
            "calibration": "GSE63577_WI38",
            "sealed_external_test": "GSE297406_hDF_four_biological_pairs",
            "external_test_labels_used_for_fit_selection_or_calibration": False,
            "dataset_overlap_count": 0,
            "lab_overlap_count": 0,
        },
        "sealed_external_dataset_lab_test": {
            "metrics": test_metrics,
            "paired_effect": paired,
            "conformal_90": {
                "calibration_sample_count": int(calibration.sum()),
                "test_sample_count": int(test.sum()),
                "threshold": conformal_threshold,
                "coverage": coverage,
                "mean_set_size": float(test_sets.sum(axis=1).mean()),
                "small_sample_authority": False,
            },
            "ood_rejection_fraction": float((external_distance > ood_threshold).mean()),
        },
        "GSE282425_single_donor_context_shift_diagnostic": {
            "geometry_results": diagnostic_results,
            "ood_rejection_fraction": float((diagnostic_distance > ood_threshold).mean()),
            "authority": "diagnostic_only_one_donor_no_biological_replicates",
        },
        "ood": {
            "distance": "RMS_fit_standardized_feature_distance",
            "threshold_calibrated_on_GSE63577_WI38": ood_threshold,
            "fit_median": float(np.median(fit_distance)),
        },
        "independent_dataset_holdout": True,
        "independent_lab_holdout": True,
        "independent_lab_direction_reproduced": paired["positive_pair_count"] == paired["pair_count"],
        "external_confirmatory_threshold_passed": paired["confirmatory_p_below_0_05"],
        "production_eligible": False,
        "status": "external_lab_direction_test_small_n_pending_confirmation",
        "observable_evidence_eligible": False,
        "may_select_aleph_parameter": False,
        "aleph_authority": "none",
        "eventual_sweep_role": "transcriptomic_senescence_observable_candidate_pending_confirmatory_external_replication_and_Aleph_operator",
        "limitations": [
            "GSE297406 has four biological pairs, so even four concordant directions give a one-sided sign-test p value of 0.0625",
            "passage-15 aging is not identical to irradiation-induced, deep, or replicative endpoint senescence",
            "the 35-gene panel validates a transcriptomic observable and does not externally validate the SenSCOUT morphology encoder",
            "GSE282425 contains one donor and supplies 2D/3D context-shift diagnostics only",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = train(args.tensor, args.checkpoint)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["sealed_external_dataset_lab_test"], sort_keys=True))


if __name__ == "__main__":
    main()
