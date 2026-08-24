#!/usr/bin/env python3
"""Sealed-lab evaluation for paired, multisite live-cell perturbation data.

The experimental unit is an independent person/experiment, never a cell or a
time point.  Lab 1 learns a perturbation direction, lab 2 selects a decision
threshold, and lab 3 is touched once for the reported holdout.  Centering a
batch on its matched controls is part of the declared paired assay design; the
result is not a classifier for samples lacking a matched control arm.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from train_calcium_baseline import binary_auroc, macro_f1


def _load_replicates(path: Path) -> tuple[np.ndarray, list[dict], list[str]]:
    values: dict[tuple[str, str, str, str, str], list[float]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            key = (
                row["lab_group"], row["biological_replicate"],
                row["technical_replicate"], row["condition"], row["observable"],
            )
            values[key].append(float(row["value"]))
    observables = sorted({key[-1] for key in values})
    feature_index = {name: index for index, name in enumerate(observables)}
    grouped: dict[tuple[str, str, str, str], np.ndarray] = {}
    for key, items in values.items():
        replicate = key[:4]
        grouped.setdefault(replicate, np.full(len(observables), np.nan))
        grouped[replicate][feature_index[key[-1]]] = np.mean(items)
    metadata, matrix = [], []
    for (lab, experiment, technical, condition), vector in sorted(grouped.items()):
        if not np.isfinite(vector).all():
            missing = [observables[i] for i in np.flatnonzero(~np.isfinite(vector))]
            raise ValueError(f"incomplete replicate {lab}/{experiment}/{technical}: {missing}")
        metadata.append({
            "lab": lab, "experiment": experiment,
            "technical_replicate": technical, "condition": condition,
        })
        matrix.append(vector)
    if not matrix:
        raise ValueError("manifest contains no observations")
    return np.stack(matrix), metadata, observables


def _paired_design(
    matrix: np.ndarray,
    metadata: list[dict],
    control: str,
    treated: str,
) -> tuple[np.ndarray, np.ndarray, list[dict], dict[str, int]]:
    batches: dict[tuple[str, str], dict[str, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for index, item in enumerate(metadata):
        if item["condition"] in {control, treated}:
            batches[(item["lab"], item["experiment"])][item["condition"]].append(index)
    centered, labels, centered_meta, experiment_counts = [], [], [], defaultdict(int)
    for (lab, experiment), arms in sorted(batches.items()):
        if not arms[control] or not arms[treated]:
            raise ValueError(f"unpaired experiment: {lab}/{experiment}")
        control_center = matrix[arms[control]].mean(axis=0)
        for condition, label in ((control, 0), (treated, 1)):
            for index in arms[condition]:
                centered.append(matrix[index] - control_center)
                labels.append(label)
                centered_meta.append({**metadata[index], "paired_control_count": len(arms[control])})
        experiment_counts[lab] += 1
    return (
        np.stack(centered), np.asarray(labels, dtype=np.int64), centered_meta,
        dict(experiment_counts),
    )


def _effect_vectors(
    centered: np.ndarray, labels: np.ndarray, metadata: list[dict]
) -> tuple[np.ndarray, list[dict]]:
    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, item in enumerate(metadata):
        groups[(item["lab"], item["experiment"])].append(index)
    vectors, result_meta = [], []
    for (lab, experiment), indices in sorted(groups.items()):
        selected = np.asarray(indices)
        if not np.any(labels[selected] == 1):
            raise ValueError(f"experiment lacks treated arm: {lab}/{experiment}")
        vectors.append(centered[selected][labels[selected] == 1].mean(axis=0))
        result_meta.append({"lab": lab, "experiment": experiment})
    return np.stack(vectors), result_meta


def _balanced_metrics(scores: np.ndarray, labels: np.ndarray, threshold: float) -> dict:
    # An exactly zero projection is the declared no-effect origin and belongs
    # to the control class.
    predicted = (scores > threshold).astype(np.int64)
    return {
        "sample_count": int(len(labels)),
        "balanced_accuracy": float(np.mean([
            np.mean(predicted[labels == label] == label) for label in (0, 1)
        ])),
        "macro_f1": macro_f1(predicted, labels, 2),
        "auroc": binary_auroc(scores[labels == 0], scores[labels == 1]),
        "threshold": float(threshold),
        "class_counts": {str(label): int(np.sum(labels == label)) for label in (0, 1)},
    }


def _bootstrap_mean(values: np.ndarray, draws: int, seed: int) -> list[float]:
    rng = np.random.default_rng(seed)
    estimates = values[rng.integers(len(values), size=(draws, len(values)))].mean(axis=1)
    return [float(value) for value in np.quantile(estimates, [0.025, 0.975])]


def _person_cluster_effects(
    scores: np.ndarray, metadata: list[dict]
) -> dict[str, np.ndarray]:
    clusters: dict[tuple[str, str], list[float]] = defaultdict(list)
    for score, item in zip(scores, metadata):
        parts = item["experiment"].split(":")
        persons = [part for part in parts if part.startswith("person_")]
        if len(persons) != 1:
            raise ValueError(
                "biological_replicate must encode exactly one ':person_<id>:' segment"
            )
        clusters[(item["lab"], persons[0])].append(float(score))
    result: dict[str, list[float]] = defaultdict(list)
    for (lab, _person), values in sorted(clusters.items()):
        result[lab].append(float(np.mean(values)))
    return {lab: np.asarray(values) for lab, values in result.items()}


def _mondrian_conformal(
    calibration_scores: np.ndarray,
    calibration_labels: np.ndarray,
    test_scores: np.ndarray,
    test_labels: np.ndarray,
    *,
    alpha: float = 0.10,
) -> dict:
    sets: list[list[int]] = []
    for score in test_scores:
        predicted_set = []
        for label in (0, 1):
            selected = calibration_scores[calibration_labels == label]
            calibration_nonconformity = selected if label == 0 else -selected
            test_nonconformity = score if label == 0 else -score
            p_value = (
                1 + np.sum(calibration_nonconformity >= test_nonconformity)
            ) / (len(calibration_nonconformity) + 1)
            if p_value > alpha:
                predicted_set.append(label)
        sets.append(predicted_set)
    sizes = np.asarray([len(item) for item in sets])
    singleton = sizes == 1
    coverage = float(np.mean([
        int(label) in predicted_set
        for label, predicted_set in zip(test_labels, sets)
    ]))
    class_coverage = {
        str(label): float(np.mean([
            label in predicted_set
            for observed, predicted_set in zip(test_labels, sets)
            if observed == label
        ]))
        for label in (0, 1)
    }
    singleton_accuracy = float(np.mean([
        predicted_set[0] == int(label)
        for label, predicted_set, is_singleton in zip(test_labels, sets, singleton)
        if is_singleton
    ])) if singleton.any() else 0.0
    passed = (
        coverage >= 1 - alpha
        and min(class_coverage.values()) >= 1 - alpha
        and float(np.mean(singleton)) >= 0.30
        and singleton_accuracy >= 0.70
        and not np.any(sizes == 0)
    )
    return {
        "method": "class_conditional_Mondrian_split_conformal",
        "alpha": alpha,
        "calibration_count": int(len(calibration_labels)),
        "test_count": int(len(test_labels)),
        "coverage": coverage,
        "class_coverage": class_coverage,
        "mean_set_size": float(np.mean(sizes)),
        "singleton_fraction": float(np.mean(singleton)),
        "singleton_accuracy": singleton_accuracy,
        "empty_set_fraction": float(np.mean(sizes == 0)),
        "uncertainty_holdout_passed": passed,
    }


def evaluate(
    manifest: Path,
    *,
    control: str = "control",
    treated: str = "ROCK_inhibited",
    draws: int = 10_000,
) -> dict:
    matrix, metadata, observables = _load_replicates(manifest)
    centered, labels, centered_meta, experiment_counts = _paired_design(
        matrix, metadata, control, treated
    )
    labs = sorted(experiment_counts)
    if len(labs) != 3:
        raise ValueError(f"sealed protocol requires exactly three labs, found {labs}")
    train_lab, validation_lab, test_lab = labs

    train = np.asarray([item["lab"] == train_lab for item in centered_meta])
    validation = np.asarray([item["lab"] == validation_lab for item in centered_meta])
    test = np.asarray([item["lab"] == test_lab for item in centered_meta])
    # Scale is learned from lab 1 only. The small floor prevents a nearly
    # constant feature from dominating the direction.
    scale = np.std(centered[train], axis=0)
    floor = max(float(np.median(scale[scale > 0])) * 0.05, 1e-12)
    scale = np.maximum(scale, floor)
    standardized = centered / scale
    direction = standardized[train & (labels == 1)].mean(axis=0)
    norm = float(np.linalg.norm(direction))
    if norm <= 0:
        raise ValueError("training perturbation direction is zero")
    direction /= norm
    scores = standardized @ direction

    # Paired-control centering declares zero as the no-effect origin. Keeping
    # this threshold fixed avoids fitting any model choice on lab 2 and leaves
    # all of lab 2 available for independent validation and calibration.
    threshold = 0.0
    validation_metrics = _balanced_metrics(
        scores[validation], labels[validation], threshold
    )
    test_metrics = _balanced_metrics(scores[test], labels[test], threshold)
    conformal = _mondrian_conformal(
        scores[validation], labels[validation], scores[test], labels[test]
    )

    effects, effect_meta = _effect_vectors(standardized, labels, centered_meta)
    effect_scores = effects @ direction
    lab_effects = {
        lab: effect_scores[[item["lab"] == lab for item in effect_meta]]
        for lab in labs
    }
    person_effects = _person_cluster_effects(effect_scores, effect_meta)
    test_interval = _bootstrap_mean(person_effects[test_lab], draws, seed=20260805)

    raw_train = matrix[[item["lab"] == train_lab for item in metadata]]
    raw_test = matrix[[item["lab"] == test_lab for item in metadata]]
    raw_mean, raw_scale = raw_train.mean(axis=0), raw_train.std(axis=0)
    raw_scale[raw_scale < 1e-12] = 1.0
    raw_distance_train = np.mean(((raw_train - raw_mean) / raw_scale) ** 2, axis=1)
    raw_distance_test = np.mean(((raw_test - raw_mean) / raw_scale) ** 2, axis=1)
    corrected_distance_train = np.mean(standardized[train] ** 2, axis=1)
    corrected_distance_test = np.mean(standardized[test] ** 2, axis=1)
    ood = {
        "raw_lab3_vs_lab1_auroc": binary_auroc(raw_distance_train, raw_distance_test),
        "paired_centered_lab3_vs_lab1_auroc": binary_auroc(
            corrected_distance_train, corrected_distance_test
        ),
    }

    passed = (
        validation_metrics["balanced_accuracy"] >= 0.70
        and validation_metrics["macro_f1"] >= 0.70
        and validation_metrics["auroc"] >= 0.70
        and test_metrics["balanced_accuracy"] >= 0.70
        and test_metrics["macro_f1"] >= 0.70
        and test_metrics["auroc"] >= 0.70
        and test_interval[0] > 0
        and conformal["uncertainty_holdout_passed"]
    )
    return {
        "task": "paired_ROCK_inhibition_live_cell_state_transfer",
        "protocol": {
            "train_lab": train_lab,
            "validation_and_conformal_calibration_lab": validation_lab,
            "external_test_lab": test_lab,
            "selection_unit": "technical_replicate",
            "inference_unit": "independent_experiment_nested_within_person",
            "bootstrap_resampling_unit": "person_cluster_mean_of_independent_experiments",
            "paired_control_centering": True,
            "decision_threshold": 0.0,
            "decision_threshold_rule": "fixed_no_effect_origin_not_label_selected",
            "test_lab_labels_used_for_numeric_fitting": False,
            "test_lab_was_opened_before_current_revision": True,
            "confirmation_status": (
                "retrospective_external_lab_evaluation_requires_new_lab"
            ),
        },
        "observables": observables,
        "observable_count": len(observables),
        "experiment_counts": experiment_counts,
        "person_counts": {lab: int(len(person_effects[lab])) for lab in labs},
        "validation_metrics": validation_metrics,
        "external_lab3_metrics": test_metrics,
        "conformal_uncertainty": conformal,
        "experiment_effect_projection": {
            lab: {
                "count": int(len(values)),
                "mean": float(np.mean(values)),
                "positive_fraction": float(np.mean(values > 0)),
            }
            for lab, values in lab_effects.items()
        },
        "person_cluster_effect_projection": {
            lab: {
                "count": int(len(values)),
                "mean": float(np.mean(values)),
                "positive_fraction": float(np.mean(values > 0)),
            }
            for lab, values in person_effects.items()
        },
        "external_lab3_effect_bootstrap_95ci": test_interval,
        "dataset_shift": ood,
        "external_lab_holdout_passed": passed,
        "pristine_sealed_confirmation": False,
        "publication_confirmatory_holdout_passed": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "the model requires a matched control arm in every experiment",
            "technical replicates are used for prediction metrics but independent experiments are the inference unit",
            "ROCK inhibition is a pathway perturbation and does not identify one Aleph parameter",
            "lab 3 was opened under the earlier threshold-selected revision, so the current fixed-origin conformal revision requires a new confirmation lab",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--control", default="control")
    parser.add_argument("--treated", default="ROCK_inhibited")
    args = parser.parse_args()
    report = evaluate(args.manifest, control=args.control, treated=args.treated)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["external_lab3_metrics"], sort_keys=True))


if __name__ == "__main__":
    main()
