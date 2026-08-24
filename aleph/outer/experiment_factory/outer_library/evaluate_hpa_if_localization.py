#!/usr/bin/env python3
"""Evaluate an HPA IF localization head under gene/antibody grouped splits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata


SCHEMA = "aleph.outer_library.hpa_if_localization_evaluation.v1"


def _tokens(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def _auc(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.int8)
    positives = y == 1
    negatives = ~positives
    if not positives.any() or not negatives.any():
        return float("nan")
    ranks = rankdata(score, method="average")
    return float(
        (ranks[positives].sum() - positives.sum() * (positives.sum() + 1) / 2)
        / (positives.sum() * negatives.sum())
    )


def _average_precision(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.int8)
    positives = int(y.sum())
    if positives == 0:
        return float("nan")
    order = np.argsort(-score, kind="stable")
    ordered = y[order]
    precision = np.cumsum(ordered) / np.arange(1, len(ordered) + 1)
    return float((precision * ordered).sum() / positives)


def _encode_labels(locations: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    names = np.asarray(sorted({x for value in locations for x in _tokens(str(value))}), dtype="U64")
    index = {name: i for i, name in enumerate(names)}
    y = np.zeros((len(locations), len(names)), dtype=np.int8)
    for row, value in enumerate(locations):
        for label in _tokens(str(value)):
            y[row, index[label]] = 1
    return y, names


def _fit_ridge(X: np.ndarray, y: np.ndarray, train: np.ndarray) -> np.ndarray:
    mean = X[train].mean(axis=0)
    scale = X[train].std(axis=0)
    scale[scale < 1e-8] = 1.0
    train_z = (X[train] - mean) / scale
    design = np.column_stack((np.ones(train.sum()), train_z))
    penalty = np.eye(design.shape[1], dtype=np.float64) * 25.0
    penalty[0, 0] = 0.0
    weights = np.linalg.solve(design.T @ design + penalty, design.T @ y[train])
    prediction = np.column_stack((np.ones(len(X)), (X - mean) / scale)) @ weights
    return np.clip(prediction, 0.0, 1.0)


def _thresholds(y: np.ndarray, probability: np.ndarray, validation: np.ndarray) -> np.ndarray:
    thresholds = np.full(y.shape[1], 0.5, dtype=np.float64)
    grid = np.linspace(0.05, 0.95, 37)
    for label in range(y.shape[1]):
        truth = y[validation, label]
        if truth.min() == truth.max():
            continue
        best = (-1.0, 0.5)
        for threshold in grid:
            predicted = probability[validation, label] >= threshold
            tp = int((predicted & (truth == 1)).sum())
            fp = int((predicted & (truth == 0)).sum())
            fn = int((~predicted & (truth == 1)).sum())
            f1 = 2 * tp / max(2 * tp + fp + fn, 1)
            candidate = (f1, -abs(threshold - 0.5))
            if candidate > (best[0], -abs(best[1] - 0.5)):
                best = (f1, float(threshold))
        thresholds[label] = best[1]
    return thresholds


def _metrics(y: np.ndarray, probability: np.ndarray, predicted: np.ndarray, mask: np.ndarray) -> dict[str, object]:
    truth = y[mask]
    pred = predicted[mask]
    tp = (pred & (truth == 1)).sum(axis=0)
    fp = (pred & (truth == 0)).sum(axis=0)
    fn = ((~pred) & (truth == 1)).sum(axis=0)
    per_label_f1 = 2 * tp / np.maximum(2 * tp + fp + fn, 1)
    eligible = (truth.sum(axis=0) > 0) & (truth.sum(axis=0) < len(truth))
    aucs = np.asarray([_auc(truth[:, i], probability[mask, i]) for i in range(y.shape[1])])
    aps = np.asarray([_average_precision(truth[:, i], probability[mask, i]) for i in range(y.shape[1])])
    micro_tp, micro_fp, micro_fn = int(tp.sum()), int(fp.sum()), int(fn.sum())
    return {
        "image_count": int(mask.sum()),
        "evaluable_label_count": int(eligible.sum()),
        "macro_f1": float(per_label_f1[eligible].mean()),
        "micro_f1": float(2 * micro_tp / max(2 * micro_tp + micro_fp + micro_fn, 1)),
        "macro_auroc": float(np.nanmean(aucs[eligible])),
        "macro_average_precision": float(np.nanmean(aps[eligible])),
        "exact_match_accuracy": float(np.all(pred == truth, axis=1).mean()),
    }


def evaluate(tensor_path: Path) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        required = {
            "X", "file_prefix", "cell_line", "antibody", "genes", "locations",
            "component_id", "split", "source_sha256", "dataset_id", "lab_group",
        }
        if not required.issubset(data.files):
            raise ValueError(f"HPA tensor misses fields: {sorted(required - set(data.files))}")
        X = np.asarray(data["X"], dtype=np.float64)
        locations = np.asarray(data["locations"]).astype(str)
        genes = np.asarray(data["genes"]).astype(str)
        antibody = np.asarray(data["antibody"]).astype(str)
        component = np.asarray(data["component_id"]).astype(str)
        split = np.asarray(data["split"]).astype(str)
        cell_line = np.asarray(data["cell_line"]).astype(str)
        dataset_id = str(np.asarray(data["dataset_id"])[0])
        lab_group = str(np.asarray(data["lab_group"])[0])
        source_sha256 = str(np.asarray(data["source_sha256"])[0])
    if X.shape != (81_007, 128) or not np.isfinite(X).all():
        raise ValueError(f"unexpected HPA embedding tensor shape/content: {X.shape}")
    masks = {name: split == name for name in ("train", "validation", "calibration", "test")}
    if any(mask.sum() == 0 for mask in masks.values()):
        raise ValueError("HPA grouped split is empty")
    entity_splits: dict[str, set[str]] = {}
    for genes_text, antibody_id, split_name in zip(genes, antibody, split):
        for entity in [f"antibody:{antibody_id}", *[f"gene:{x}" for x in _tokens(genes_text)]]:
            entity_splits.setdefault(entity, set()).add(split_name)
    overlap_count = sum(len(values) != 1 for values in entity_splits.values())
    if overlap_count:
        raise ValueError("HPA gene/antibody grouped split leaked")

    y, label_names = _encode_labels(locations)
    probability = _fit_ridge(X, y, masks["train"])
    thresholds = _thresholds(y, probability, masks["validation"])
    predicted = probability >= thresholds
    validation_metrics = _metrics(y, probability, predicted, masks["validation"])
    test_metrics = _metrics(y, probability, predicted, masks["test"])

    # Calibration is a diagnostic only: row exchangeability is weakened by
    # repeated images within a biological reagent component, and there is no
    # independent provider/lab test here.
    calibration_residual = np.abs(y[masks["calibration"]] - probability[masks["calibration"]])
    alpha = 0.10
    quantile_index = min(
        calibration_residual.size - 1,
        int(np.ceil((calibration_residual.size + 1) * (1.0 - alpha))) - 1,
    )
    radius = float(np.partition(calibration_residual.ravel(), quantile_index)[quantile_index])
    test_covered = np.abs(y[masks["test"]] - probability[masks["test"]]) <= radius

    label_metrics: dict[str, object] = {}
    test = masks["test"]
    for index, name in enumerate(label_names):
        positives = int(y[test, index].sum())
        test_auc = _auc(y[test, index], probability[test, index])
        test_ap = _average_precision(y[test, index], probability[test, index])
        label_metrics[str(name)] = {
            "test_positive_images": positives,
            "threshold_selected_on_validation": float(thresholds[index]),
            "test_auroc": float(test_auc) if np.isfinite(test_auc) else None,
            "test_average_precision": float(test_ap) if np.isfinite(test_ap) else None,
        }

    return {
        "schema": SCHEMA,
        "task": "HPA_IF_subcellular_localization_multilabel_representation",
        "dataset_id": dataset_id,
        "lab_group": lab_group,
        "source_sha256": source_sha256,
        "image_count": int(len(X)),
        "cell_line_count": int(len(set(cell_line))),
        "author_label_count": int(len(label_names)),
        "embedding_dimensions": int(X.shape[1]),
        "model": "fixed_HPA_embedding_view_plus_train_only_ridge_multilabel_head",
        "split": {
            "unit": "connected_component_of_gene_antibody_bipartite_graph",
            "image_counts": {name: int(mask.sum()) for name, mask in masks.items()},
            "component_counts": {
                name: len(set(component[mask])) for name, mask in masks.items()
            },
            "gene_or_antibody_cross_split_overlap_count": overlap_count,
            "test_labels_used_for_fitting_or_threshold_selection": False,
        },
        "validation_metrics": validation_metrics,
        "sealed_grouped_test_metrics": test_metrics,
        "per_label_test_metrics": label_metrics,
        "conformal_diagnostic": {
            "alpha": alpha,
            "marginal_probability_interval_radius": radius,
            "test_entrywise_coverage": float(test_covered.mean()),
            "eligible_for_uncertainty_authority": False,
            "reason": "same-provider repeated-image rows are not independent experimental units",
        },
        "dataset_holdout_present": True,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "status": "training_only_grouped_same_provider",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "limitations": [
            "HPA embeddings and localization labels share one provider and annotation pipeline",
            "the test split seals genes and antibodies but is not an independent laboratory",
            "the 128-D view is a compact diagnostic, not a newly trained raw-image encoder",
            "cell-line imbalance prevents interpreting image-level metrics as universal cell-state accuracy",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["sealed_grouped_test_metrics"], sort_keys=True))


if __name__ == "__main__":
    main()
