#!/usr/bin/env python3
"""Evaluate BBBC054 microglia patches under field-grouped splits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr


SCHEMA = "aleph.outer_library.bbbc054_microglia_evaluation.v1"


def _split(field: int) -> str:
    # One field from each repeating ten-field block enters calibration/test,
    # preserving coverage across the full ordered LPS series.
    bucket = field % 10
    if bucket < 6:
        return "train"
    if bucket < 8:
        return "validation"
    if bucket == 8:
        return "calibration"
    return "test"


def _metrics(truth: np.ndarray, predicted: np.ndarray, probability: np.ndarray) -> dict[str, float | int]:
    classes = probability.shape[1]
    recall, f1 = [], []
    for label in range(classes):
        tp = int(((truth == label) & (predicted == label)).sum())
        fp = int(((truth != label) & (predicted == label)).sum())
        fn = int(((truth == label) & (predicted != label)).sum())
        recall.append(tp / max(tp + fn, 1))
        f1.append(2 * tp / max(2 * tp + fp + fn, 1))
    return {
        "cell_count": int(len(truth)),
        "accuracy": float((truth == predicted).mean()),
        "balanced_accuracy": float(np.mean(recall)),
        "macro_f1": float(np.mean(f1)),
        "multiclass_brier": float(np.mean(np.sum((probability - np.eye(classes)[truth]) ** 2, axis=1))),
    }


def _fit_lda(X: np.ndarray, y: np.ndarray, train: np.ndarray, shrinkage: float):
    mean = X[train].mean(axis=0)
    scale = X[train].std(axis=0)
    scale[scale < 1e-8] = 1.0
    z = (X - mean) / scale
    classes = np.unique(y)
    centroids = np.stack([z[train & (y == label)].mean(axis=0) for label in classes])
    residual = np.concatenate([z[train & (y == label)] - centroids[label] for label in classes])
    covariance = np.cov(residual, rowvar=False)
    covariance = (1.0 - shrinkage) * covariance + shrinkage * np.eye(covariance.shape[0])
    inverse = np.linalg.inv(covariance)
    distance = np.stack([
        np.einsum("ij,jk,ik->i", z - centroid, inverse, z - centroid)
        for centroid in centroids
    ], axis=1)
    scores = -0.5 * distance  # equal class prior protects the rare amoeboid class
    scores -= scores.max(axis=1, keepdims=True)
    probability = np.exp(scores)
    probability /= probability.sum(axis=1, keepdims=True)
    return probability


def _field_composition(field: np.ndarray, y: np.ndarray, labels: np.ndarray) -> dict[str, object]:
    fields = np.unique(field)
    fractions = np.zeros((len(fields), len(labels)), dtype=np.float64)
    for row, value in enumerate(fields):
        mask = field == value
        fractions[row] = np.bincount(y[mask], minlength=len(labels)) / mask.sum()
    result: dict[str, object] = {}
    for index, label in enumerate(labels):
        early = fractions[fields < 10, index]
        late = fractions[fields >= 50, index]
        rng = np.random.default_rng(20260806 + index)
        bootstrap = np.asarray([
            rng.choice(late, len(late), replace=True).mean()
            - rng.choice(early, len(early), replace=True).mean()
            for _ in range(5_000)
        ])
        result[str(label)] = {
            "spearman_fraction_vs_author_FOV_index": float(
                spearmanr(fields, fractions[:, index]).statistic
            ),
            "early_fields_0_to_9_mean_fraction": float(early.mean()),
            "late_fields_50_to_59_mean_fraction": float(late.mean()),
            "late_minus_early_fraction": float(late.mean() - early.mean()),
            "field_bootstrap_ci95": [float(x) for x in np.quantile(bootstrap, (0.025, 0.975))],
        }
    return result


def evaluate(tensor_path: Path) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        X = np.asarray(data["X"], dtype=np.float64)
        label_text = data["label"].astype(str)
        field = np.asarray(data["field"], dtype=np.int16)
        sample_id = data["sample_id"].astype(str)
    if X.shape != (58_186, 12) or len(set(sample_id)) != len(sample_id):
        raise ValueError("BBBC054 tensor shape or sample identity contract changed")
    labels = np.asarray(sorted(set(label_text)), dtype="U16")
    y = np.asarray([int(np.flatnonzero(labels == value)[0]) for value in label_text], dtype=np.int8)
    split = np.asarray([_split(int(value)) for value in field])
    masks = {name: split == name for name in ("train", "validation", "calibration", "test")}
    field_sets = {name: set(field[mask]) for name, mask in masks.items()}
    if any(field_sets[left] & field_sets[right] for left in field_sets for right in field_sets if left < right):
        raise ValueError("BBBC054 field leakage")

    candidates: list[dict[str, object]] = []
    for shrinkage in (0.05, 0.20, 0.50, 0.80):
        probability = _fit_lda(X, y, masks["train"], shrinkage)
        metrics = _metrics(y[masks["validation"]], probability[masks["validation"]].argmax(axis=1), probability[masks["validation"]])
        candidates.append({"shrinkage": shrinkage, "validation_metrics": metrics})
    selected = max(candidates, key=lambda item: (item["validation_metrics"]["macro_f1"], -item["shrinkage"]))
    probability = _fit_lda(X, y, masks["train"], float(selected["shrinkage"]))
    predicted = probability.argmax(axis=1)

    calibration_scores = 1.0 - probability[masks["calibration"], y[masks["calibration"]]]
    alpha = 0.10
    order = min(
        len(calibration_scores) - 1,
        int(np.ceil((len(calibration_scores) + 1) * (1.0 - alpha))) - 1,
    )
    threshold = float(np.sort(calibration_scores)[order])
    prediction_set = probability[masks["test"]] >= (1.0 - threshold)
    test_truth = y[masks["test"]]
    covered = prediction_set[np.arange(len(test_truth)), test_truth]

    return {
        "schema": SCHEMA,
        "task": "LPS_microglia_brightfield_morphophenotype_representation",
        "dataset_id": "bbbc054-v1-LPS-microglia-timecourse",
        "image_fields": 60,
        "annotated_cells": int(len(y)),
        "author_labels": labels.tolist(),
        "author_webpage_label_conflict": {
            "webpage_says": ["round", "amoeboid", "stratified"],
            "downloaded_CSV_contains": labels.tolist(),
            "resolution": "downloaded author CSV labels retained verbatim; no stratified-to-ramified rewrite",
        },
        "split": {
            "unit": "whole_image_field",
            "rule": "field_index_modulo_10:0-5 train,6-7 validation,8 calibration,9 test",
            "field_counts": {name: len(values) for name, values in field_sets.items()},
            "cell_counts": {name: int(mask.sum()) for name, mask in masks.items()},
            "field_overlap_count": 0,
            "test_labels_used_for_model_or_shrinkage_selection": False,
            "pristine_confirmation": False,
            "pristine_reason": "an earlier hash split was opened and found to contain only two test fields; this balanced field-count revision is retrospective",
        },
        "model": "12_feature_brightfield_patch_shrinkage_LDA_equal_class_prior",
        "selection_candidates": candidates,
        "selected_shrinkage": selected["shrinkage"],
        "validation_metrics": selected["validation_metrics"],
        "sealed_field_test_metrics": _metrics(
            test_truth, predicted[masks["test"]], probability[masks["test"]]
        ),
        "split_conformal_diagnostic": {
            "alpha": alpha,
            "coverage": float(covered.mean()),
            "mean_prediction_set_size": float(prediction_set.sum(axis=1).mean()),
            "empty_prediction_set_fraction": float((prediction_set.sum(axis=1) == 0).mean()),
            "eligible_for_uncertainty_authority": False,
        },
        "descriptive_LPS_timecourse": _field_composition(field, y, labels),
        "time_coordinate": "ordered 30-minute FOV series; absolute half-hour offset is not explicit in downloaded metadata",
        "independent_biological_replicate_holdout": False,
        "independent_lab_holdout": False,
        "status": "training_only_single_annotated_replicate",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "limitations": [
            "only replicate 1 has author cell phenotype annotations; replicates 2 and 3 cannot be labelled by inheritance",
            "field-grouped test prevents cell leakage but remains one time series from one biological replicate",
            "no matched unstimulated time-course control is provided",
            "the FOV-to-clock offset is not encoded in the downloaded annotation",
            "morphophenotype labels are not direct molecular immune-state measurements",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["sealed_field_test_metrics"], sort_keys=True))


if __name__ == "__main__":
    main()
