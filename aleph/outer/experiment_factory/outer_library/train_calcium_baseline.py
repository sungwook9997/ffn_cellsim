#!/usr/bin/env python3
"""Leakage-safe cell-state baseline for the public HeLa calcium time series."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def load_curves(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["observable"] == "intracellular_calcium_F_over_F0":
                grouped[row["sample_id"]].append(row)
    samples, labels, curves = [], [], []
    for sample_id, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: row["coordinate_value"])
        times = np.asarray([row["coordinate_value"] for row in rows], dtype=np.float64)
        values = np.asarray([row["value"] for row in rows], dtype=np.float64)
        if len(values) != 226 or not np.all(np.diff(times) > 0):
            raise ValueError(f"invalid curve for {sample_id}")
        samples.append(sample_id)
        labels.append(rows[0]["condition"])
        curves.append(values)
    classes = sorted(set(labels))
    class_index = {label: index for index, label in enumerate(classes)}
    return (
        np.asarray(samples),
        np.asarray([class_index[label] for label in labels], dtype=np.int64),
        np.stack(curves),
        classes,
    )


def load_unlabelled_curves(path: Path) -> tuple[np.ndarray, list[np.ndarray]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            grouped[row["sample_id"]].append(row)
    samples, curves = [], []
    for sample_id, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: row["coordinate_value"])
        samples.append(sample_id)
        curves.append(np.asarray([row["value"] for row in rows], dtype=np.float64))
    return np.asarray(samples), curves


def curve_features(curves: np.ndarray) -> np.ndarray:
    indices = np.linspace(0, curves.shape[1] - 1, 32).round().astype(int)
    sampled = curves[:, indices]
    summaries = np.column_stack((
        curves.mean(axis=1), curves.std(axis=1), curves.min(axis=1), curves.max(axis=1),
        curves[:, -1], curves[:, -1] - curves[:, 0],
        np.argmin(curves, axis=1) / (curves.shape[1] - 1),
        np.argmax(curves, axis=1) / (curves.shape[1] - 1),
        np.trapezoid(curves, axis=1) / (curves.shape[1] - 1),
    ))
    return np.column_stack((sampled, summaries))


def stratified_hash_split(samples: np.ndarray, labels: np.ndarray) -> dict[str, np.ndarray]:
    parts = {"train": [], "validation": [], "test": []}
    for label in sorted(set(labels.tolist())):
        indices = np.flatnonzero(labels == label).tolist()
        indices.sort(key=lambda i: hashlib.sha256(samples[i].encode()).hexdigest())
        if len(indices) != 50:
            raise ValueError("the public source must contain 50 cells per condition")
        parts["train"].extend(indices[:30])
        parts["validation"].extend(indices[30:40])
        parts["test"].extend(indices[40:])
    return {name: np.asarray(sorted(indices), dtype=np.int64) for name, indices in parts.items()}


def fit_ridge(x: np.ndarray, y: np.ndarray, class_count: int, penalty: float = 1.0) -> dict:
    mean = x.mean(axis=0)
    scale = x.std(axis=0)
    scale[scale < 1e-12] = 1.0
    z = (x - mean) / scale
    design = np.column_stack((np.ones(len(z)), z))
    target = np.eye(class_count)[y]
    regularizer = np.eye(design.shape[1]) * penalty
    regularizer[0, 0] = 0.0
    weights = np.linalg.solve(design.T @ design + regularizer, design.T @ target)
    return {"mean": mean, "scale": scale, "weights": weights}


def logits(model: dict, x: np.ndarray) -> np.ndarray:
    z = (x - model["mean"]) / model["scale"]
    return np.column_stack((np.ones(len(z)), z)) @ model["weights"]


def probabilities(raw_logits: np.ndarray, temperature: float) -> np.ndarray:
    shifted = raw_logits / temperature
    shifted -= shifted.max(axis=1, keepdims=True)
    values = np.exp(shifted)
    return values / values.sum(axis=1, keepdims=True)


def nll(probs: np.ndarray, labels: np.ndarray) -> float:
    return float(-np.log(np.clip(probs[np.arange(len(labels)), labels], 1e-12, 1)).mean())


def calibrate_temperature(raw_logits: np.ndarray, labels: np.ndarray) -> float:
    grid = np.geomspace(0.05, 20.0, 300)
    return float(min(grid, key=lambda value: nll(probabilities(raw_logits, value), labels)))


def macro_f1(predictions: np.ndarray, labels: np.ndarray, class_count: int) -> float:
    scores = []
    for cls in range(class_count):
        tp = np.sum((predictions == cls) & (labels == cls))
        fp = np.sum((predictions == cls) & (labels != cls))
        fn = np.sum((predictions != cls) & (labels == cls))
        scores.append(2 * tp / max(2 * tp + fp + fn, 1))
    return float(np.mean(scores))


def ece(probs: np.ndarray, labels: np.ndarray, bins: int = 10) -> float:
    confidence = probs.max(axis=1)
    correct = probs.argmax(axis=1) == labels
    result = 0.0
    for lower in np.linspace(0, 1, bins, endpoint=False):
        upper = lower + 1 / bins
        selected = (confidence > lower) & (confidence <= upper)
        if selected.any():
            result += selected.mean() * abs(confidence[selected].mean() - correct[selected].mean())
    return float(result)


def binary_auroc(id_scores: np.ndarray, ood_scores: np.ndarray) -> float:
    # Mann-Whitney interpretation, with OOD expected to have the larger score.
    comparisons = (ood_scores[:, None] > id_scores[None, :]).sum()
    ties = (ood_scores[:, None] == id_scores[None, :]).sum()
    return float((comparisons + 0.5 * ties) / (len(id_scores) * len(ood_scores)))


def semantic_ood(curves: np.ndarray, labels: np.ndarray, features: np.ndarray) -> dict:
    distance_scores = []
    confidence_scores = []
    for held_out in sorted(set(labels.tolist())):
        held = labels == held_out
        remaining = sorted(set(labels.tolist()) - {held_out})
        train_indices, id_indices = [], []
        for label in remaining:
            indices = np.flatnonzero(labels == label)
            train_indices.extend(indices[:40])
            id_indices.extend(indices[40:])
        train_indices = np.asarray(train_indices)
        id_indices = np.asarray(id_indices)
        remap = {label: index for index, label in enumerate(remaining)}
        train_labels = np.asarray([remap[label] for label in labels[train_indices]])
        model = fit_ridge(features[train_indices], train_labels, 2)
        id_uncertainty = 1 - probabilities(logits(model, features[id_indices]), 1.0).max(axis=1)
        ood_uncertainty = 1 - probabilities(logits(model, features[held]), 1.0).max(axis=1)
        confidence_scores.append(binary_auroc(id_uncertainty, ood_uncertainty))
        train_z = (features[train_indices] - model["mean"]) / model["scale"]
        id_z = (features[id_indices] - model["mean"]) / model["scale"]
        ood_z = (features[held] - model["mean"]) / model["scale"]
        id_distance = ((id_z[:, None] - train_z[None]) ** 2).mean(axis=2).min(axis=1)
        ood_distance = ((ood_z[:, None] - train_z[None]) ** 2).mean(axis=2).min(axis=1)
        distance_scores.append(binary_auroc(id_distance, ood_distance))
    return {
        "protocol": "leave_one_author_condition_out",
        "primary_score": "nearest_training_curve_distance",
        "distance_fold_aurocs": distance_scores,
        "confidence_fold_aurocs": confidence_scores,
        "mean_auroc": float(np.mean(distance_scores)),
        "external_lab_ood": False,
    }


def train(
    manifest: Path,
    external_ood_manifest: Path | None = None,
    external_task_report: Path | None = None,
) -> dict:
    samples, labels, curves, classes = load_curves(manifest)
    features = curve_features(curves)
    split = stratified_hash_split(samples, labels)
    model = fit_ridge(features[split["train"]], labels[split["train"]], len(classes))
    val_logits = logits(model, features[split["validation"]])
    temperature = calibrate_temperature(val_logits, labels[split["validation"]])
    test_probs = probabilities(logits(model, features[split["test"]]), temperature)
    test_labels = labels[split["test"]]
    predictions = test_probs.argmax(axis=1)
    external_ood = {
        "evaluated": False,
        "reason": "no independent-lab calcium manifest supplied",
        "dataset_count": 0,
        "lab_count": 0,
    }
    if external_ood_manifest is not None:
        external_samples, external_curves = load_unlabelled_curves(external_ood_manifest)
        external_features = np.stack([
            curve_features(curve[None, :])[0] for curve in external_curves
        ])
        id_uncertainty = 1 - probabilities(val_logits, temperature).max(axis=1)
        external_uncertainty = 1 - probabilities(
            logits(model, external_features), temperature
        ).max(axis=1)
        train_features = features[split["train"]]
        train_labels = labels[split["train"]]
        train_z = (train_features - model["mean"]) / model["scale"]
        val_z = (features[split["validation"]] - model["mean"]) / model["scale"]
        external_z = (external_features - model["mean"]) / model["scale"]
        centroids = np.stack([
            train_z[train_labels == label].mean(axis=0) for label in range(len(classes))
        ])
        id_distance = ((val_z[:, None] - centroids[None]) ** 2).mean(axis=2).min(axis=1)
        external_distance = (
            ((external_z[:, None] - centroids[None]) ** 2).mean(axis=2).min(axis=1)
        )
        external_ood = {
            "evaluated": True,
            "protocol": (
                "independent-lab MDCK-II author-aggregate curves; "
                "OOD score=1-max_probability"
            ),
            "sample_count": len(external_samples),
            "primary_score": "nearest_class_centroid_distance",
            "auroc": binary_auroc(id_distance, external_distance),
            "max_probability_auroc": binary_auroc(id_uncertainty, external_uncertainty),
            "dataset_count": 1,
            "lab_count": 1,
            "task_labels_aligned": False,
        }
    task_holdout_count = 0
    task_holdout_passed_count = 0
    task_holdout_summary = None
    if external_task_report is not None:
        task_holdout_summary = json.loads(external_task_report.read_text(encoding="utf-8"))
        task_holdout_count = int(task_holdout_summary["external_lab_holdout"])
        task_holdout_passed_count = int(task_holdout_summary["task_holdout_passed"])
    return {
        "task": "author_condition_classification_from_single_cell_calcium_time_series",
        "classes": classes,
        "sample_count": len(samples),
        "split": {name: len(indices) for name, indices in split.items()},
        "split_unit": "cell_sample_id",
        "sample_overlap_count": 0,
        "temperature": temperature,
        "test": {
            "accuracy": float(np.mean(predictions == test_labels)),
            "macro_f1": macro_f1(predictions, test_labels, len(classes)),
            "nll": nll(test_probs, test_labels),
            "ece_10_bin": ece(test_probs, test_labels),
            "mean_predictive_entropy": float(
                np.mean(-np.sum(test_probs * np.log(np.clip(test_probs, 1e-12, 1)), axis=1))
            ),
        },
        "semantic_ood": semantic_ood(curves, labels, features),
        "external_lab_ood": external_ood,
        "independent_dataset_ood_count": external_ood["dataset_count"],
        "independent_lab_ood_count": external_ood["lab_count"],
        "independent_dataset_task_holdout_count": task_holdout_count,
        "independent_lab_task_holdout_count": task_holdout_count,
        "independent_dataset_task_holdout_passed_count": task_holdout_passed_count,
        "independent_lab_task_holdout_passed_count": task_holdout_passed_count,
        "external_task_holdout": task_holdout_summary,
        "production_eligible": False,
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--external-ood-manifest", type=Path)
    parser.add_argument("--external-task-report", type=Path)
    args = parser.parse_args()
    report = train(args.manifest, args.external_ood_manifest, args.external_task_report)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["test"], sort_keys=True))


if __name__ == "__main__":
    main()
