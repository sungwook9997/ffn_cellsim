#!/usr/bin/env python3
"""Grouped migration-mode transfer audit for the Dryad 9jh6m row tensor.

The model never sees behavioral features used to define the author motion-mode
labels.  Fibronectin experiment dates are partitioned into train, validation,
and conformal-calibration groups; every talin experiment date and tracked cell
is then evaluated as a same-lab, different-perturbation retrospective holdout.
This is useful representation evidence, but it is explicitly not Lab 4.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from train_calcium_baseline import binary_auroc, macro_f1


def _balanced_metrics(scores: np.ndarray, labels: np.ndarray) -> dict:
    predicted = (scores > 0).astype(np.int64)
    return {
        "count": int(len(labels)),
        "class_counts": {str(label): int(np.sum(labels == label)) for label in (0, 1)},
        "balanced_accuracy": float(np.mean([
            np.mean(predicted[labels == label] == label) for label in (0, 1)
        ])),
        "macro_f1": macro_f1(predicted, labels, 2),
        "auroc": binary_auroc(scores[labels == 0], scores[labels == 1]),
    }


def _fit_ridge(
    train_x: np.ndarray,
    train_y: np.ndarray,
    validation_x: np.ndarray,
    validation_y: np.ndarray,
) -> tuple[np.ndarray, dict, np.ndarray, np.ndarray]:
    median = np.nanmedian(train_x, axis=0)
    median[~np.isfinite(median)] = 0.0

    def impute(values: np.ndarray) -> np.ndarray:
        return np.where(np.isfinite(values), values, median)

    train = impute(train_x)
    mean, scale = train.mean(axis=0), train.std(axis=0)
    scale[scale < 1e-12] = 1.0

    def design(values: np.ndarray) -> np.ndarray:
        standardized = (impute(values) - mean) / scale
        return np.column_stack([np.ones(len(standardized)), standardized])

    train_design = design(train_x)
    validation_design = design(validation_x)
    target = train_y.astype(np.float64) * 2.0 - 1.0
    # Each author mode contributes equal total mass even when the number of
    # repeated cell-time observations differs strongly between modes.
    class_count = np.bincount(train_y, minlength=2).astype(np.float64)
    sample_weight = 0.5 / class_count[train_y]
    weighted_design = train_design * np.sqrt(sample_weight[:, None])
    weighted_target = target * np.sqrt(sample_weight)
    candidates = []
    for ridge in (0.01, 0.1, 1.0, 10.0, 100.0):
        penalty = np.eye(train_design.shape[1]) * ridge
        penalty[0, 0] = 0.0
        weights = np.linalg.solve(
            weighted_design.T @ weighted_design + penalty,
            weighted_design.T @ weighted_target,
        )
        metrics = _balanced_metrics(validation_design @ weights, validation_y)
        candidates.append((metrics["balanced_accuracy"], metrics["macro_f1"], -ridge, weights, metrics))
    _, _, negative_ridge, weights, validation_metrics = max(candidates, key=lambda x: x[:3])
    validation_metrics = {
        **validation_metrics,
        "selected_ridge": float(-negative_ridge),
        "selection_rule": "maximum_validation_balanced_accuracy_then_macro_f1_then_smaller_ridge",
    }
    return weights, validation_metrics, median, np.stack([mean, scale])


def _scores(
    values: np.ndarray,
    weights: np.ndarray,
    median: np.ndarray,
    normalization: np.ndarray,
) -> np.ndarray:
    imputed = np.where(np.isfinite(values), values, median)
    standardized = (imputed - normalization[0]) / normalization[1]
    return np.column_stack([np.ones(len(values)), standardized]) @ weights


def _cell_aggregate(
    scores: np.ndarray,
    labels: np.ndarray,
    cells: np.ndarray,
    dates: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    result_scores, result_labels, result_cells, result_dates = [], [], [], []
    for cell in sorted(set(cells.tolist())):
        selected = cells == cell
        counts = np.bincount(labels[selected], minlength=2)
        if counts[0] == counts[1]:
            continue
        result_scores.append(float(np.median(scores[selected])))
        result_labels.append(int(np.argmax(counts)))
        result_cells.append(cell)
        date_values = set(dates[selected].tolist())
        if len(date_values) != 1:
            raise ValueError(f"tracked cell crosses experimental dates: {cell}")
        result_dates.append(next(iter(date_values)))
    return (
        np.asarray(result_scores), np.asarray(result_labels, dtype=np.int64),
        np.asarray(result_cells), np.asarray(result_dates),
    )


def _mondrian(
    calibration_scores: np.ndarray,
    calibration_labels: np.ndarray,
    test_scores: np.ndarray,
    test_labels: np.ndarray,
    alpha: float = 0.10,
) -> dict:
    sets: list[list[int]] = []
    for score in test_scores:
        prediction_set = []
        for label in (0, 1):
            selected = calibration_scores[calibration_labels == label]
            if not len(selected):
                raise ValueError(f"calibration lacks class {label}")
            calibration_nonconformity = selected if label == 0 else -selected
            test_nonconformity = score if label == 0 else -score
            p_value = (1 + np.sum(calibration_nonconformity >= test_nonconformity)) / (
                len(calibration_nonconformity) + 1
            )
            if p_value > alpha:
                prediction_set.append(label)
        sets.append(prediction_set)
    sizes = np.asarray([len(item) for item in sets])
    singleton = sizes == 1
    coverage_by_class = {
        str(label): float(np.mean([
            label in prediction_set
            for observed, prediction_set in zip(test_labels, sets)
            if observed == label
        ]))
        for label in (0, 1)
    }
    singleton_accuracy = float(np.mean([
        prediction_set[0] == int(observed)
        for observed, prediction_set, is_singleton in zip(test_labels, sets, singleton)
        if is_singleton
    ])) if singleton.any() else 0.0
    coverage = float(np.mean([
        int(observed) in prediction_set
        for observed, prediction_set in zip(test_labels, sets)
    ]))
    return {
        "method": "tracked_cell_level_class_conditional_Mondrian_split_conformal",
        "alpha": alpha,
        "calibration_cells": int(len(calibration_labels)),
        "test_cells": int(len(test_labels)),
        "coverage": coverage,
        "class_coverage": coverage_by_class,
        "mean_set_size": float(np.mean(sizes)),
        "singleton_fraction": float(np.mean(singleton)),
        "singleton_accuracy": singleton_accuracy,
        "empty_set_fraction": float(np.mean(sizes == 0)),
        "same_lab_uncertainty_transfer_passed": bool(
            coverage >= 1 - alpha
            and min(coverage_by_class.values()) >= 1 - alpha
            and float(np.mean(singleton)) >= 0.30
            and singleton_accuracy >= 0.70
            and not np.any(sizes == 0)
        ),
    }


def _date_cluster_interval(
    scores: np.ndarray,
    labels: np.ndarray,
    dates: np.ndarray,
    *,
    draws: int,
    seed: int = 20260805,
) -> list[float]:
    unique_dates = sorted(set(dates.tolist()))
    per_date = []
    for date in unique_dates:
        selected = dates == date
        per_date.append(_balanced_metrics(scores[selected], labels[selected])["balanced_accuracy"])
    values = np.asarray(per_date)
    rng = np.random.default_rng(seed)
    estimates = values[rng.integers(len(values), size=(draws, len(values)))].mean(axis=1)
    return [float(value) for value in np.quantile(estimates, [0.025, 0.975])]


def evaluate(tensor: Path, *, draws: int = 10_000) -> dict:
    with np.load(tensor, allow_pickle=False) as data:
        schema = str(data["schema"])
        if schema != "aleph.outer_library.grouped_cell_rows.v1":
            raise ValueError(f"unexpected tensor schema: {schema}")
        feature_names = data["feature_names"].astype(str)
        behavioral = set(data["behavioral_features"].astype(str).tolist())
        selected_feature = np.asarray([name not in behavioral for name in feature_names])
        values = data["values"][:, selected_feature]
        dataset = data["dataset"].astype(str)
        dates = data["experimental_date"].astype(str)
        cells = data["tracked_cell"].astype(str)
        mode = data["author_motion_mode"].astype(np.int64)

    labelled = np.isin(mode, [1, 2]) & np.isin(dataset, ["fibronectin", "talin"])
    labels = mode[labelled] - 1
    values, dataset, dates, cells = values[labelled], dataset[labelled], dates[labelled], cells[labelled]
    fibronectin_dates = sorted(set(dates[dataset == "fibronectin"].tolist()))
    if len(fibronectin_dates) != 19:
        raise ValueError(f"expected 19 fibronectin dates, found {len(fibronectin_dates)}")
    # These date groups are fixed from the fibronectin source alone to retain
    # both author modes in validation and calibration.  No talin label enters
    # the partition or model selection.
    validation_dates = {"20110617", "20110705", "20110803", "20110810"}
    calibration_dates = {"20110918", "20120110", "20120128", "20120210"}
    train_dates = set(fibronectin_dates) - validation_dates - calibration_dates
    train = (dataset == "fibronectin") & np.isin(dates, list(train_dates))
    validation = (dataset == "fibronectin") & np.isin(dates, list(validation_dates))
    calibration = (dataset == "fibronectin") & np.isin(dates, list(calibration_dates))
    test = dataset == "talin"
    if any(set(cells[left].tolist()) & set(cells[right].tolist()) for left, right in (
        (train, validation), (train, calibration), (train, test),
        (validation, calibration), (validation, test), (calibration, test),
    )):
        raise ValueError("tracked-cell leakage across partitions")

    weights, validation_metrics, median, normalization = _fit_ridge(
        values[train], labels[train], values[validation], labels[validation]
    )
    calibration_scores = _scores(values[calibration], weights, median, normalization)
    test_scores = _scores(values[test], weights, median, normalization)
    test_observation_metrics = _balanced_metrics(test_scores, labels[test])
    calibration_cell = _cell_aggregate(
        calibration_scores, labels[calibration], cells[calibration], dates[calibration]
    )
    test_cell = _cell_aggregate(test_scores, labels[test], cells[test], dates[test])
    test_cell_metrics = _balanced_metrics(test_cell[0], test_cell[1])
    conformal = _mondrian(calibration_cell[0], calibration_cell[1], test_cell[0], test_cell[1])

    # Training-normalized squared distance diagnoses perturbation/domain shift;
    # labels from the talin holdout do not enter this calculation.
    calibration_x = np.where(np.isfinite(values[calibration]), values[calibration], median)
    test_x = np.where(np.isfinite(values[test]), values[test], median)
    calibration_distance = np.mean(((calibration_x - normalization[0]) / normalization[1]) ** 2, axis=1)
    test_distance = np.mean(((test_x - normalization[0]) / normalization[1]) ** 2, axis=1)
    same_lab_transfer_passed = bool(
        validation_metrics["balanced_accuracy"] >= 0.70
        and validation_metrics["macro_f1"] >= 0.70
        and validation_metrics["auroc"] >= 0.70
        and test_cell_metrics["balanced_accuracy"] >= 0.70
        and test_cell_metrics["macro_f1"] >= 0.70
        and test_cell_metrics["auroc"] >= 0.70
        and conformal["same_lab_uncertainty_transfer_passed"]
    )
    return {
        "task": "organizational_features_predict_author_migration_mode_across_perturbations",
        "source": "Dryad 10.5061/dryad.9jh6m / eLife 10.7554/eLife.11384",
        "protocol": {
            "training_dataset": "fibronectin",
            "train_experimental_dates": sorted(train_dates),
            "validation_experimental_dates": sorted(validation_dates),
            "conformal_calibration_experimental_dates": sorted(calibration_dates),
            "retrospective_test_dataset": "talin",
            "retrospective_test_experimental_dates": sorted(set(dates[test].tolist())),
            "split_group": "tracked_cell_nested_in_author_experimental_date",
            "label_defining_behavioral_features_excluded": sorted(behavioral),
            "fibronectin_date_partition_rule": (
                "fixed_source_only_stratification_retaining_both_author_modes;"
                "talin_labels_not_used"
            ),
            "selected_organizational_feature_count": int(selected_feature.sum()),
            "test_labels_used_for_numeric_fitting": False,
        },
        "partition_observation_counts": {
            "train": int(train.sum()), "validation": int(validation.sum()),
            "calibration": int(calibration.sum()), "test": int(test.sum()),
        },
        "partition_tracked_cell_counts": {
            "train": len(set(cells[train].tolist())),
            "validation": len(set(cells[validation].tolist())),
            "calibration": len(set(cells[calibration].tolist())),
            "test": len(set(cells[test].tolist())),
        },
        "validation_observation_metrics": validation_metrics,
        "retrospective_talin_observation_metrics": test_observation_metrics,
        "retrospective_talin_tracked_cell_metrics": test_cell_metrics,
        "retrospective_talin_date_cluster_balanced_accuracy_95ci": _date_cluster_interval(
            test_cell[0], test_cell[1], test_cell[3], draws=draws
        ),
        "conformal_uncertainty": conformal,
        "dataset_shift": {
            "talin_vs_fibronectin_calibration_distance_auroc": binary_auroc(
                calibration_distance, test_distance
            ),
            "interpretation": "domain-shift diagnostic, not a motion-mode performance metric",
        },
        "same_lab_cross_perturbation_transfer_passed": same_lab_transfer_passed,
        "independent_lab_holdout": False,
        "pristine_lab4_confirmation": False,
        "eligible_as_same_lab_representation_training": True,
        "eligible_as_independent_confirmation": False,
        "aleph_parameter_proposal_eligible": False,
        "may_mutate_physics": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=10_000)
    args = parser.parse_args()
    report = evaluate(args.tensor, draws=args.draws)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "same_lab_cross_perturbation_transfer_passed": report[
            "same_lab_cross_perturbation_transfer_passed"
        ],
        "tracked_cell_metrics": report["retrospective_talin_tracked_cell_metrics"],
        "conformal": report["conformal_uncertainty"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
