#!/usr/bin/env python3
"""Training-only conformal uncertainty and dataset-shift audit for Piezo1 transfer."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from evaluate_piezo1_transfer import _invariant_features, _metrics
from train_calcium_baseline import (
    binary_auroc,
    calibrate_temperature,
    fit_ridge,
    logits,
    probabilities,
)


def _split(labels: np.ndarray, sample_ids: list[str]) -> dict[str, np.ndarray]:
    parts = {"train": [], "validation": [], "calibration": []}
    for label in (0, 1):
        indices = np.flatnonzero(labels == label).tolist()
        indices.sort(key=lambda index: hashlib.sha256(
            sample_ids[index].encode()
        ).hexdigest())
        if len(indices) != 50:
            raise ValueError("expected 50 HeLa cells per class")
        parts["train"].extend(indices[:25])
        parts["validation"].extend(indices[25:35])
        parts["calibration"].extend(indices[35:])
    return {
        name: np.asarray(sorted(indices), dtype=np.int64)
        for name, indices in parts.items()
    }


def _select_penalty(
    features: np.ndarray,
    labels: np.ndarray,
    split: dict[str, np.ndarray],
) -> tuple[float, dict]:
    best = None
    for penalty in np.geomspace(0.001, 1000.0, 25):
        model = fit_ridge(features[split["train"]], labels[split["train"]], 2,
                          float(penalty))
        probs = probabilities(logits(model, features[split["validation"]]), 1.0)
        metrics = _metrics(probs, labels[split["validation"]])
        key = (metrics["balanced_accuracy"], metrics["macro_f1"], -float(penalty))
        if best is None or key > best[0]:
            best = (key, float(penalty), metrics, model)
    assert best is not None
    return best[1], {"metrics": best[2], "model": best[3]}


def _quantile_higher(values: np.ndarray, probability: float) -> float:
    ordered = np.sort(np.asarray(values, dtype=np.float64))
    index = min(max(math.ceil(probability * len(ordered)) - 1, 0), len(ordered) - 1)
    return float(ordered[index])


def evaluate(
    train_manifest: Path,
    external_manifest: Path,
    alpha: float = 0.10,
) -> dict:
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0, 1)")
    x, y, _, sample_ids = _invariant_features(
        train_manifest, "hela", active_only=False
    )
    external_x, external_y, conditions, external_ids = _invariant_features(
        external_manifest, "external", active_only=True
    )
    split = _split(y, sample_ids)
    penalty, selected = _select_penalty(x, y, split)
    model = selected["model"]
    validation_logits = logits(model, x[split["validation"]])
    temperature = calibrate_temperature(
        validation_logits, y[split["validation"]]
    )
    calibration_probs = probabilities(
        logits(model, x[split["calibration"]]), temperature
    )
    calibration_scores = 1.0 - calibration_probs[
        np.arange(len(split["calibration"])), y[split["calibration"]]
    ]
    probability = min(1.0, math.ceil((len(calibration_scores) + 1) * (1 - alpha))
                      / len(calibration_scores))
    pooled_threshold = _quantile_higher(calibration_scores, probability)

    external_probs = probabilities(logits(model, external_x), temperature)
    pooled_membership = (1.0 - external_probs) <= pooled_threshold
    pooled_sizes = pooled_membership.sum(axis=1)
    pooled_covered = pooled_membership[np.arange(len(external_y)), external_y]

    # Mondrian (class-conditional) calibration prevents the majority external
    # class from borrowing the tighter nonconformity distribution of the other
    # class. Thresholds still use HeLa calibration labels only.
    class_thresholds = {}
    for label in (0, 1):
        selected_scores = calibration_scores[y[split["calibration"]] == label]
        class_probability = min(
            1.0,
            math.ceil((len(selected_scores) + 1) * (1 - alpha))
            / len(selected_scores),
        )
        class_thresholds[label] = _quantile_higher(
            selected_scores, class_probability
        )
    membership = np.column_stack([
        (1.0 - external_probs[:, label]) <= class_thresholds[label]
        for label in (0, 1)
    ])
    set_sizes = membership.sum(axis=1)
    covered = membership[np.arange(len(external_y)), external_y]
    singleton = set_sizes == 1
    singleton_predictions = membership.argmax(axis=1)
    singleton_accuracy = (
        float(np.mean(singleton_predictions[singleton] == external_y[singleton]))
        if singleton.any() else None
    )

    train_z = (x[split["train"]] - model["mean"]) / model["scale"]
    calibration_z = (x[split["calibration"]] - model["mean"]) / model["scale"]
    external_z = (external_x - model["mean"]) / model["scale"]
    centroids = np.stack([
        train_z[y[split["train"]] == label].mean(axis=0) for label in (0, 1)
    ])
    calibration_distance = (
        (calibration_z[:, None] - centroids[None]) ** 2
    ).mean(axis=2).min(axis=1)
    external_distance = (
        (external_z[:, None] - centroids[None]) ** 2
    ).mean(axis=2).min(axis=1)
    ood_threshold = _quantile_higher(calibration_distance, 0.95)
    ood_auroc = binary_auroc(calibration_distance, external_distance)
    external_ood_rate = float(np.mean(external_distance > ood_threshold))

    coverage = float(np.mean(covered))
    singleton_fraction = float(np.mean(singleton))
    uncertainty_passed = (
        coverage >= 1 - alpha
        and singleton_fraction >= 0.30
        and singleton_accuracy is not None
        and singleton_accuracy >= 0.70
    )
    return {
        "task": "Piezo1_external_lab_training_only_uncertainty_audit",
        "training_dataset": "zenodo-18495219 HeLa",
        "external_dataset": "zenodo-19890985 EA.hy926",
        "split": {name: len(indices) for name, indices in split.items()},
        "selection_protocol": (
            "penalty and temperature use HeLa train/validation only; the disjoint HeLa "
            "calibration split sets the conformal threshold; no EA.hy926 label tunes anything"
        ),
        "selected_penalty": penalty,
        "temperature": temperature,
        "validation_metrics": selected["metrics"],
        "conformal": {
            "method": "Mondrian_class_conditional_split_conformal_classification_set",
            "alpha": alpha,
            "calibration_count": len(calibration_scores),
            "nonconformity": "1_minus_probability_of_true_class",
            "class_thresholds": {
                str(label): threshold for label, threshold in class_thresholds.items()
            },
            "external_count": len(external_y),
            "external_coverage": coverage,
            "mean_prediction_set_size": float(set_sizes.mean()),
            "singleton_fraction": singleton_fraction,
            "empty_set_fraction": float(np.mean(set_sizes == 0)),
            "singleton_accuracy": singleton_accuracy,
            "target_coverage_met": coverage >= 1 - alpha,
            "class_conditional_external_coverage": {
                str(label): float(np.mean(covered[external_y == label]))
                for label in (0, 1)
            },
            "pooled_diagnostic": {
                "threshold": pooled_threshold,
                "external_coverage": float(np.mean(pooled_covered)),
                "mean_prediction_set_size": float(pooled_sizes.mean()),
                "target_coverage_met": float(np.mean(pooled_covered)) >= 1 - alpha,
            },
        },
        "dataset_shift": {
            "score": "nearest_training_class_centroid_standardized_distance",
            "calibration_95pct_threshold": ood_threshold,
            "external_ood_rate": external_ood_rate,
            "external_vs_calibration_auroc": ood_auroc,
            "shift_detected": ood_auroc >= 0.70 and external_ood_rate >= 0.50,
        },
        "external_conditions": sorted(set(conditions)),
        "sample_overlap_count": len(set(sample_ids) & set(external_ids)),
        "uncertainty_holdout_passed": uncertainty_passed,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "split-conformal coverage is not guaranteed under the detected cross-lab shift",
            "the result applies to batches containing both active nanoswitch treatment arms",
            "biological replicate identities are unavailable for external cells",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train-manifest", type=Path, required=True)
    parser.add_argument("--external-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.10)
    args = parser.parse_args()
    report = evaluate(args.train_manifest, args.external_manifest, args.alpha)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "conformal": report["conformal"],
        "dataset_shift": report["dataset_shift"],
        "uncertainty_holdout_passed": report["uncertainty_holdout_passed"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
