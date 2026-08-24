#!/usr/bin/env python3
"""Train leakage-labelled SenSCOUT morphology heads across biorepeats.

The model intentionally keeps image groups intact.  Bio1 is the only fitting
repeat, Bio2 is a sealed same-provider repeat test, and Bio3 is a frozen shift
diagnostic.  None of these partitions constitutes an external-lab holdout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def _partition_groups(groups: np.ndarray) -> np.ndarray:
    unique = sorted(set(groups.tolist()))
    ranked = sorted(
        unique,
        key=lambda value: hashlib.sha256(
            f"SenSCOUT-Bio1-v1:{value}".encode("utf-8")
        ).digest(),
    )
    count = len(ranked)
    train_end = max(1, int(round(0.60 * count)))
    selection_end = max(train_end + 1, int(round(0.80 * count)))
    mapping = {
        group: (0 if i < train_end else 1 if i < selection_end else 2)
        for i, group in enumerate(ranked)
    }
    return np.asarray([mapping[group] for group in groups], dtype=np.int8)


def _transform(X: np.ndarray, fit: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    signed_log = np.sign(X) * np.log1p(np.abs(X))
    mean = signed_log[fit].mean(axis=0)
    scale = signed_log[fit].std(axis=0)
    scale[scale < 1e-6] = 1.0
    return ((signed_log - mean) / scale).astype(np.float32), mean, scale


def _sigmoid(value: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(value, -30.0, 30.0)))


def _binary_metrics(probability: np.ndarray, labels: np.ndarray) -> dict[str, float | int]:
    prediction = probability >= 0.5
    tp = int(np.sum(prediction & (labels == 1)))
    tn = int(np.sum(~prediction & (labels == 0)))
    fp = int(np.sum(prediction & (labels == 0)))
    fn = int(np.sum(~prediction & (labels == 1)))
    recall_pos = tp / max(1, tp + fn)
    recall_neg = tn / max(1, tn + fp)
    return {
        "row_count": int(len(labels)),
        "accuracy": float((prediction == labels).mean()),
        "balanced_accuracy": float((recall_pos + recall_neg) / 2),
        "f1": float(2 * tp / max(1, 2 * tp + fp + fn)),
        "brier": float(np.mean((probability - labels) ** 2)),
        "true_prevalence": float(labels.mean()),
        "predicted_prevalence": float(probability.mean()),
    }


def _prevalence_by_group(
    probability: np.ndarray, labels: np.ndarray, groups: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    names = np.asarray(sorted(set(groups.tolist())))
    predicted = np.asarray([probability[groups == name].mean() for name in names])
    observed = np.asarray([labels[groups == name].mean() for name in names])
    return names, predicted, observed


def _fit_logistic(
    X: np.ndarray, y: np.ndarray, train: np.ndarray, selection: np.ndarray
) -> tuple[np.ndarray, float, int]:
    candidates = (1e-4, 1e-3, 1e-2, 1e-1)
    best = None
    for l2 in candidates:
        weights = np.zeros(X.shape[1], dtype=np.float64)
        bias = 0.0
        positive_weight = float(np.sum(y[train] == 0) / max(1, np.sum(y[train] == 1)))
        row_weight = np.where(y[train] == 1, positive_weight, 1.0)
        row_weight /= row_weight.mean()
        selected_epoch = 0
        selected_score = -1.0
        selected_params = (weights.copy(), bias)
        for epoch in range(1, 1201):
            probability = _sigmoid(X[train] @ weights + bias)
            error = (probability - y[train]) * row_weight
            weights -= 0.04 * (X[train].T @ error / train.sum() + l2 * weights)
            bias -= 0.04 * float(error.mean())
            if epoch % 20 == 0:
                score = _binary_metrics(
                    _sigmoid(X[selection] @ weights + bias), y[selection]
                )["balanced_accuracy"]
                if float(score) > selected_score:
                    selected_score = float(score)
                    selected_epoch = epoch
                    selected_params = (weights.copy(), bias)
        candidate = (selected_score, -l2, selected_params, selected_epoch, l2)
        if best is None or candidate[:2] > best[:2]:
            best = candidate
    assert best is not None
    weights, bias = best[2]
    return np.append(weights, bias), float(best[4]), int(best[3])


def _temperature(logits: np.ndarray, labels: np.ndarray) -> float:
    choices = np.linspace(0.35, 4.0, 147)
    losses = []
    for temperature in choices:
        probability = _sigmoid(logits / temperature)
        losses.append(float(-np.mean(
            labels * np.log(probability + 1e-12)
            + (1 - labels) * np.log(1 - probability + 1e-12)
        )))
    return float(choices[int(np.argmin(losses))])


def _quantile_higher(values: np.ndarray, coverage: float) -> float:
    rank = min(len(values) - 1, int(np.ceil((len(values) + 1) * coverage)) - 1)
    return float(np.partition(values, rank)[rank])


def _regression_metrics(predicted: np.ndarray, observed: np.ndarray) -> dict[str, float | int]:
    residual = observed - predicted
    denominator = float(np.sum((observed - observed.mean()) ** 2))
    correlation = float(np.corrcoef(predicted, observed)[0, 1]) if len(observed) > 1 else float("nan")
    return {
        "row_count": int(len(observed)),
        "r2": float(1.0 - np.sum(residual ** 2) / denominator) if denominator > 0 else float("nan"),
        "pearson_r": correlation,
        "mae_log1p": float(np.mean(np.abs(residual))),
    }


def _fit_ridge_heads(
    X: np.ndarray,
    repeat: np.ndarray,
    group: np.ndarray,
    biomarker: np.ndarray,
    marker_value: np.ndarray,
    p21_value: np.ndarray,
) -> dict[str, object]:
    targets = {
        "HMGB1": (biomarker == "HMGB1", marker_value),
        "LMNB1": (biomarker == "LMNB1", marker_value),
        "P16": (biomarker == "P16", marker_value),
        "P21": (np.isfinite(p21_value), p21_value),
    }
    reports: dict[str, object] = {}
    for name, (eligible, raw_target) in targets.items():
        fit_domain = eligible & (repeat == "Bio1") & np.isfinite(raw_target)
        partitions = _partition_groups(group[fit_domain])
        indices = np.flatnonzero(fit_domain)
        train_idx, select_idx, calibration_idx = (
            indices[partitions == part] for part in range(3)
        )
        y = np.log1p(raw_target.astype(np.float64))
        candidates = (0.01, 0.1, 1.0, 10.0, 100.0)
        best = None
        for penalty in candidates:
            design = np.column_stack([X[train_idx], np.ones(len(train_idx))])
            regularizer = np.eye(design.shape[1]) * penalty
            regularizer[-1, -1] = 0.0
            coefficient = np.linalg.solve(
                design.T @ design + regularizer, design.T @ y[train_idx]
            )
            predicted = np.column_stack([X[select_idx], np.ones(len(select_idx))]) @ coefficient
            score = float(np.mean(np.abs(predicted - y[select_idx])))
            if best is None or score < best[0]:
                best = (score, penalty, coefficient)
        assert best is not None
        coefficient = best[2]
        predict = lambda idx: np.column_stack([X[idx], np.ones(len(idx))]) @ coefficient
        calibration_residual = np.abs(predict(calibration_idx) - y[calibration_idx])
        interval_radius = _quantile_higher(calibration_residual, 0.90)
        repeat_reports = {}
        for target_repeat in ("Bio2", "Bio3"):
            test_idx = np.flatnonzero(
                eligible & (repeat == target_repeat) & np.isfinite(raw_target)
            )
            predicted = predict(test_idx)
            metrics = _regression_metrics(predicted, y[test_idx])
            metrics["interval_90_coverage"] = float(
                (np.abs(predicted - y[test_idx]) <= interval_radius).mean()
            )
            metrics["same_provider_frozen_repeat"] = True
            repeat_reports[target_repeat] = metrics
        reports[name] = {
            "target_transform": "log1p",
            "selected_ridge_penalty": float(best[1]),
            "bio1_image_group_counts": {
                role: int(len(set(group[index].tolist())))
                for role, index in zip(
                    ("train", "selection", "calibration"),
                    (train_idx, select_idx, calibration_idx), strict=True
                )
            },
            "calibration_interval_90_radius_log1p": interval_radius,
            "frozen_repeat_evaluation": repeat_reports,
        }
    return reports


def train(tensor_path: Path, checkpoint_path: Path) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        X_raw = np.asarray(data["X"], dtype=np.float64)
        repeat = np.asarray(data["biorepeat"]).astype(str)
        group = np.asarray(data["image_group"]).astype(str)
        biomarker = np.asarray(data["biomarker"]).astype(str)
        labels = np.asarray(data["bgal_senescence_label"], dtype=np.int8)
        marker_value = np.asarray(data["measured_biomarker_value"], dtype=np.float64)
        p21_value = np.asarray(data["measured_p21_value"], dtype=np.float64)
        feature_name = np.asarray(data["feature_name"]).astype(str)
    if X_raw.shape != (5000, 87) or len(feature_name) != 87:
        raise ValueError("SenSCOUT morphology tensor contract changed")

    labelled_bio1 = (repeat == "Bio1") & (labels >= 0)
    partition = np.full(len(labels), -1, dtype=np.int8)
    partition[labelled_bio1] = _partition_groups(group[labelled_bio1])
    train_mask = partition == 0
    selection_mask = partition == 1
    calibration_mask = partition == 2
    if any(len(np.unique(labels[mask])) != 2 for mask in (train_mask, selection_mask, calibration_mask)):
        raise ValueError("Bio1 grouped split lost one BGAL class")
    X, mean, scale = _transform(X_raw, train_mask)
    parameter, l2, epoch = _fit_logistic(X, labels, train_mask, selection_mask)
    logits = X @ parameter[:-1] + parameter[-1]
    temperature = _temperature(logits[calibration_mask], labels[calibration_mask])
    probability = _sigmoid(logits / temperature)

    calibration_groups, calibration_predicted, calibration_observed = _prevalence_by_group(
        probability[calibration_mask], labels[calibration_mask], group[calibration_mask]
    )
    prevalence_radius = _quantile_higher(
        np.abs(calibration_predicted - calibration_observed), 0.90
    )
    bio2 = (repeat == "Bio2") & (labels >= 0)
    test_groups, test_predicted, test_observed = _prevalence_by_group(
        probability[bio2], labels[bio2], group[bio2]
    )
    group_correlation = float(np.corrcoef(test_predicted, test_observed)[0, 1])

    train_hidden = X[train_mask]
    calibration_distance = np.sqrt(np.mean(X[calibration_mask] ** 2, axis=1))
    ood_threshold = _quantile_higher(calibration_distance, 0.95)
    shift = {}
    for target_repeat in ("Bio2", "Bio3"):
        distance = np.sqrt(np.mean(X[repeat == target_repeat] ** 2, axis=1))
        shift[target_repeat] = {
            "row_count": int(len(distance)),
            "rejection_fraction": float((distance > ood_threshold).mean()),
            "median_standardized_rms_distance": float(np.median(distance)),
        }

    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        checkpoint_path,
        feature_name=feature_name,
        input_transform=np.asarray("signed_log1p"),
        input_mean=mean,
        input_scale=scale,
        bgal_weight=parameter[:-1],
        bgal_bias=np.asarray(parameter[-1]),
        temperature=np.asarray(temperature),
        prevalence_interval_radius=np.asarray(prevalence_radius),
        ood_threshold=np.asarray(ood_threshold),
    )
    return {
        "schema": "aleph.outer_library.senscout_morphology_model.v1",
        "task": "morphology_to_senescence_biomarker_state",
        "input": {
            "feature_count": 87,
            "feature_source": "author_curated_morphology_panel",
            "transform": "signed_log1p_then_Bio1_train_only_standardization",
            "target_or_author_cluster_features_used": False,
        },
        "split": {
            "fit_repeat": "Bio1",
            "sealed_same_provider_test_repeat": "Bio2",
            "frozen_shift_diagnostic_repeat": "Bio3",
            "unit": "image_group",
            "bio1_group_counts": {
                role: int(len(set(group[mask].tolist())))
                for role, mask in (
                    ("train", train_mask), ("selection", selection_mask),
                    ("calibration", calibration_mask)
                )
            },
            "cross_split_image_group_overlap_count": 0,
            "cells_treated_as_independent_biological_replicates": False,
        },
        "direct_bgal_head": {
            "family": "class_weighted_logistic_regression",
            "selected_l2": l2,
            "selection_epoch": epoch,
            "temperature_selected_on_Bio1_calibration_groups": temperature,
            "Bio2_cell_level_nonindependent_diagnostic": _binary_metrics(
                probability[bio2], labels[bio2]
            ),
            "Bio2_image_group_prevalence": {
                "image_group_count": int(len(test_groups)),
                "pearson_r": group_correlation,
                "mae": float(np.mean(np.abs(test_predicted - test_observed))),
                "conformal_90_radius": prevalence_radius,
                "conformal_90_coverage": float(
                    (np.abs(test_predicted - test_observed) <= prevalence_radius).mean()
                ),
                "calibration_image_group_count": int(len(calibration_groups)),
            },
        },
        "continuous_biomarker_heads": _fit_ridge_heads(
            X, repeat, group, biomarker, marker_value, p21_value
        ),
        "repeat_shift_diagnostic": {
            "threshold": ood_threshold,
            "calibration_source": "Bio1_disjoint_calibration_rows",
            "repeat_results": shift,
            "authority": "diagnostic_only",
        },
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "status": "training_only_same_provider_biorepeat_holdout",
        "production_eligible": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "eventual_sweep_role": "candidate_observable_encoder_pending_external_dataset_lab_validation_and_Aleph_operator_sensitivity",
        "limitations": [
            "Bio1, Bio2, and Bio3 are repeats from one SenSCOUT source rather than independent laboratories",
            "cells within an image group are correlated; cell-level metrics are diagnostic only",
            "BGAL labels are mixed within many image groups, so group prevalence is the inference-unit evaluation",
            "Bio3 substitutes BGAL_FLUOR for categorical BGAL and therefore supplies shift diagnostics, not direct-label confirmation",
            "biomarker heads predict separately sampled marker rows and do not imply same-cell multimodal pairing",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = train(args.tensor, args.checkpoint)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["direct_bgal_head"], sort_keys=True))


if __name__ == "__main__":
    main()
