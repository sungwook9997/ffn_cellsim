#!/usr/bin/env python3
"""Train a compound-unit programmed-cell-death mechanism baseline."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np


ROLES = ("train", "selection", "calibration", "sealed_test")


def _softmax(logits: np.ndarray) -> np.ndarray:
    shifted = logits - logits.max(axis=1, keepdims=True)
    exponential = np.exp(shifted)
    return exponential / exponential.sum(axis=1, keepdims=True)


def _metrics(probability: np.ndarray, labels: np.ndarray, names: np.ndarray) -> dict[str, object]:
    prediction = probability.argmax(axis=1)
    confusion = np.zeros((len(names), len(names)), dtype=np.int64)
    np.add.at(confusion, (labels, prediction), 1)
    present = confusion.sum(axis=1) > 0
    f1, recall = [], []
    for label in range(len(names)):
        tp = int(confusion[label, label]); fp = int(confusion[:, label].sum() - tp); fn = int(confusion[label].sum() - tp)
        f1.append(float(2 * tp / max(1, 2 * tp + fp + fn)) if present[label] else None)
        recall.append(float(tp / max(1, tp + fn)) if present[label] else None)
    return {
        "compound_count": int(len(labels)),
        "accuracy": float((prediction == labels).mean()),
        "balanced_accuracy": float(np.mean([value for value in recall if value is not None])),
        "macro_f1": float(np.mean([value for value in f1 if value is not None])),
        "per_class_f1": {name: value for name, value in zip(names.tolist(), f1, strict=True)},
        "confusion_matrix_true_by_predicted": confusion.tolist(),
    }


def _fit_softmax(
    X: np.ndarray, labels: np.ndarray, class_count: int, l2: float
) -> tuple[np.ndarray, np.ndarray]:
    weights = np.zeros((X.shape[1], class_count), dtype=np.float64)
    bias = np.zeros(class_count, dtype=np.float64)
    counts = np.bincount(labels, minlength=class_count)
    class_weight = np.sqrt(len(labels) / (class_count * np.maximum(counts, 1)))
    row_weight = class_weight[labels]
    row_weight /= row_weight.mean()
    for _ in range(5_000):
        probability = _softmax(X @ weights + bias)
        gradient = probability
        gradient[np.arange(len(labels)), labels] -= 1.0
        gradient *= row_weight[:, None] / len(labels)
        weights -= 0.03 * (X.T @ gradient + l2 * weights)
        bias -= 0.03 * gradient.sum(axis=0)
    return weights, bias


def _temperature(logits: np.ndarray, labels: np.ndarray) -> float:
    choices = np.linspace(0.35, 6.0, 227)
    losses = []
    for value in choices:
        probability = _softmax(logits / value)
        losses.append(float(-np.log(probability[np.arange(len(labels)), labels] + 1e-12).mean()))
    return float(choices[int(np.argmin(losses))])


def _quantile_higher(values: np.ndarray, coverage: float) -> float:
    rank = min(len(values) - 1, int(np.ceil((len(values) + 1) * coverage)) - 1)
    return float(np.partition(values, rank)[rank])


def train(tensor_path: Path, checkpoint_path: Path) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        X_well = np.asarray(data["X"], dtype=np.float64)
        compounds_well = np.asarray(data["compound"]).astype(str)
        labels_well = np.asarray(data["moa"]).astype(str)
        roles_well = np.asarray(data["compound_split"]).astype(str)
        feature_name = np.asarray(data["feature_name"]).astype(str)
    compounds = np.asarray(sorted(set(compounds_well.tolist())))
    X, label_text, roles = [], [], []
    for compound in compounds:
        mask = compounds_well == compound
        if len(set(labels_well[mask])) != 1 or len(set(roles_well[mask])) != 1:
            raise ValueError(f"compound metadata changes across wells: {compound}")
        X.append(X_well[mask].mean(axis=0)); label_text.append(labels_well[mask][0]); roles.append(roles_well[mask][0])
    X = np.asarray(X); label_text = np.asarray(label_text); roles = np.asarray(roles)
    if X.shape != (55, 672) or {role: int(np.sum(roles == role)) for role in ROLES} != {
        "train": 40, "selection": 4, "calibration": 5, "sealed_test": 6,
    }:
        raise ValueError("Figshare compound tensor contract changed")
    class_name = np.asarray(sorted(set(label_text.tolist())))
    class_index = {name: i for i, name in enumerate(class_name)}
    labels = np.asarray([class_index[name] for name in label_text], dtype=np.int64)
    train_mask = roles == "train"; selection = roles == "selection"; calibration = roles == "calibration"; test = roles == "sealed_test"
    mean = X[train_mask].mean(axis=0); scale = X[train_mask].std(axis=0); scale[scale < 1e-6] = 1.0
    standardized = (X - mean) / scale
    _, singular, right = np.linalg.svd(standardized[train_mask], full_matrices=False)

    best = None
    for dimension in (8, 16, 24, 32):
        basis = right[:dimension]
        score = standardized @ basis.T
        for l2 in (1e-3, 1e-2, 1e-1, 1.0):
            weights, bias = _fit_softmax(score[train_mask], labels[train_mask], len(class_name), l2)
            metrics = _metrics(_softmax(score[selection] @ weights + bias), labels[selection], class_name)
            candidate = (float(metrics["macro_f1"]), -dimension, l2, dimension, weights, bias, metrics)
            if best is None or candidate[:3] > best[:3]:
                best = candidate
    assert best is not None
    _, _, l2, dimension, weights, bias, selection_metrics = best
    basis = right[:dimension]
    score = standardized @ basis.T
    logits = score @ weights + bias
    temperature = _temperature(logits[calibration], labels[calibration])
    probability = _softmax(logits / temperature)
    conformal_scores = 1.0 - probability[calibration, labels[calibration]]
    conformal_threshold = _quantile_higher(conformal_scores, 0.90)
    test_sets = probability[test] >= 1.0 - conformal_threshold
    test_labels = labels[test]
    coverage = float(test_sets[np.arange(test.sum()), test_labels].mean())

    train_center = score[train_mask].mean(axis=0); train_scale = score[train_mask].std(axis=0); train_scale[train_scale < 1e-6] = 1.0
    calibration_distance = np.sqrt(np.mean(((score[calibration] - train_center) / train_scale) ** 2, axis=1))
    ood_threshold = _quantile_higher(calibration_distance, 0.95)
    test_distance = np.sqrt(np.mean(((score[test] - train_center) / train_scale) ** 2, axis=1))

    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        checkpoint_path, feature_name=feature_name, input_mean=mean, input_scale=scale,
        pca_basis=basis, pca_singular_value=singular[:dimension], weight=weights, bias=bias,
        class_name=class_name, temperature=np.asarray(temperature),
        conformal_threshold=np.asarray(conformal_threshold), ood_center=train_center,
        ood_scale=train_scale, ood_threshold=np.asarray(ood_threshold),
    )
    return {
        "schema": "aleph.outer_library.figshare_cell_death_model.v1",
        "task": "programmed_cell_death_mechanism_from_DeepProfiler_well_profiles",
        "model": {
            "family": f"train_compound_PCA_{dimension}_softmax",
            "selected_l2": float(l2), "selection_metrics": selection_metrics,
            "temperature_selected_on_calibration_compounds": temperature,
        },
        "split": {
            "unit": "compound_identity", "compound_counts": dict(Counter(roles)),
            "cross_split_compound_overlap_count": 0,
            "test_compounds_used_for_fit_selection_or_calibration": False,
            "well_rows_are_not_independent_training_units": True,
        },
        "sealed_unseen_compound_test": {
            "metrics": _metrics(probability[test], test_labels, class_name),
            "compound_ids": compounds[test].tolist(),
            "conformal_90": {
                "calibration_compound_count": int(calibration.sum()),
                "test_compound_count": int(test.sum()), "threshold": conformal_threshold,
                "coverage": coverage, "mean_set_size": float(test_sets.sum(axis=1).mean()),
                "class_conditional_authority": False,
            },
            "ood_rejection_fraction": float((test_distance > ood_threshold).mean()),
        },
        "independent_dataset_holdout": False, "independent_lab_holdout": False,
        "same_endpoint_as_biad2515_annexin_v": False,
        "status": "training_only_unseen_compound_same_provider",
        "production_eligible": False, "observable_evidence_eligible": False,
        "may_select_aleph_parameter": False, "aleph_authority": "none",
        "limitations": [
            "only 55 labelled compounds exist and the sealed test contains one compound per mechanism class",
            "selection and calibration omit rare mechanism classes, so class-conditional uncertainty is impossible",
            "all profiles share one Cell Painting study and provider pipeline",
            "mechanism-of-action labels are not interchangeable with the S-BIAD2515 Annexin-V endpoint",
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
    print(json.dumps(report["sealed_unseen_compound_test"], sort_keys=True))


if __name__ == "__main__":
    main()
