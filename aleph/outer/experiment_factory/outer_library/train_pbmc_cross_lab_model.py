#!/usr/bin/env python3
"""Train and evaluate the preregistered PBMC pre-sweep context models."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import rankdata


CELL_CLASSES = np.asarray([
    "B_cell", "CD4_T_cell", "CD8_T_cell", "NK_cell", "monocyte", "dendritic_cell"
])
STATE_CLASSES = np.asarray(["control", "IFN_stimulated"])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def encode(values: np.ndarray, classes: np.ndarray) -> np.ndarray:
    mapping = {value: index for index, value in enumerate(classes)}
    try:
        return np.asarray([mapping[value] for value in values], dtype=np.int64)
    except KeyError as error:
        raise ValueError(f"unknown class {error.args[0]!r}") from error


def softmax_from_logits(logits: np.ndarray) -> np.ndarray:
    shifted = logits - logits.max(axis=1, keepdims=True)
    probabilities = np.exp(shifted)
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    return probabilities


def logits_from_parameters(x: np.ndarray, parameters: np.ndarray, class_count: int) -> np.ndarray:
    feature_count = x.shape[1]
    matrix = parameters[: feature_count * (class_count - 1)].reshape(feature_count, class_count - 1)
    intercept = parameters[feature_count * (class_count - 1):]
    nonreference = x @ matrix + intercept
    return np.column_stack([nonreference, np.zeros(x.shape[0], dtype=nonreference.dtype)])


def fit_softmax(x: np.ndarray, y: np.ndarray, class_count: int, l2: float) -> tuple[np.ndarray, dict[str, Any]]:
    feature_count = x.shape[1]
    initial = np.zeros(feature_count * (class_count - 1) + class_count - 1, dtype=np.float64)
    n = x.shape[0]

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        logits = logits_from_parameters(x, parameters, class_count)
        probabilities = softmax_from_logits(logits)
        loss = -np.log(np.maximum(probabilities[np.arange(n), y], 1e-300)).mean()
        weights = parameters[: feature_count * (class_count - 1)].reshape(feature_count, class_count - 1)
        loss += 0.5 * l2 * float(np.sum(weights * weights))
        residual = probabilities[:, :-1]
        residual[np.arange(n), np.minimum(y, class_count - 2)] -= (y < class_count - 1)
        gradient_weights = x.T @ residual / n + l2 * weights
        gradient_intercept = residual.mean(axis=0)
        return float(loss), np.concatenate([gradient_weights.ravel(), gradient_intercept])

    result = minimize(
        objective, initial, method="L-BFGS-B", jac=True,
        options={"maxiter": 80, "ftol": 1e-10, "gtol": 1e-6, "maxls": 30},
    )
    return result.x, {
        "success": bool(result.success), "status": int(result.status),
        "iterations": int(result.nit), "final_objective": float(result.fun),
        "message": str(result.message),
    }


def probabilities(x: np.ndarray, parameters: np.ndarray, class_count: int, temperature: float = 1.0) -> np.ndarray:
    return softmax_from_logits(logits_from_parameters(x, parameters, class_count) / temperature)


def confusion_metrics(y: np.ndarray, probability: np.ndarray, classes: np.ndarray) -> dict[str, Any]:
    predicted = probability.argmax(axis=1)
    per_class = {}
    f1s = []
    for index, label in enumerate(classes):
        tp = int(np.sum((predicted == index) & (y == index)))
        fp = int(np.sum((predicted == index) & (y != index)))
        fn = int(np.sum((predicted != index) & (y == index)))
        denominator = 2 * tp + fp + fn
        f1 = 2 * tp / denominator if denominator else 0.0
        f1s.append(f1)
        per_class[str(label)] = {"support": int(np.sum(y == index)), "f1": f1}
    return {
        "cells": int(len(y)), "accuracy": float(np.mean(predicted == y)),
        "macro_f1": float(np.mean(f1s)), "per_class": per_class,
    }


def ece(y: np.ndarray, probability: np.ndarray, bins: int = 15) -> float:
    predicted = probability.argmax(axis=1)
    confidence = probability.max(axis=1)
    result = 0.0
    for lower in np.linspace(0.0, 1.0, bins + 1)[:-1]:
        upper = lower + 1.0 / bins
        mask = (confidence >= lower) & (confidence < upper if upper < 1.0 else confidence <= upper)
        if np.any(mask):
            result += float(mask.mean()) * abs(float(np.mean(predicted[mask] == y[mask])) - float(confidence[mask].mean()))
    return result


def select_temperature(logits: np.ndarray, y: np.ndarray) -> tuple[float, dict[str, Any]]:
    def objective(log_temperature: float) -> float:
        probability = softmax_from_logits(logits / np.exp(log_temperature))
        return float(-np.log(np.maximum(probability[np.arange(len(y)), y], 1e-300)).mean())
    result = minimize_scalar(objective, bounds=(np.log(0.05), np.log(20.0)), method="bounded")
    return float(np.exp(result.x)), {"success": bool(result.success), "negative_log_likelihood": float(result.fun)}


def conformal_thresholds(y: np.ndarray, probability: np.ndarray, class_count: int, alpha: float = 0.1) -> np.ndarray:
    thresholds = np.empty(class_count, dtype=np.float64)
    for class_index in range(class_count):
        scores = 1.0 - probability[y == class_index, class_index]
        if len(scores) == 0:
            raise ValueError(f"calibration donor lacks class {class_index}")
        rank = min(int(np.ceil((len(scores) + 1) * (1.0 - alpha))), len(scores))
        thresholds[class_index] = np.partition(scores, rank - 1)[rank - 1]
    return thresholds


def conformal_metrics(y: np.ndarray, probability: np.ndarray, thresholds: np.ndarray) -> dict[str, float]:
    included = (1.0 - probability) <= thresholds[None, :]
    return {
        "marginal_coverage": float(np.mean(included[np.arange(len(y)), y])),
        "mean_set_size": float(included.sum(axis=1).mean()),
        "empty_set_fraction": float(np.mean(included.sum(axis=1) == 0)),
    }


def auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = labels.astype(bool)
    positive_count = int(positives.sum())
    negative_count = int((~positives).sum())
    if not positive_count or not negative_count:
        raise ValueError("AUROC requires both classes")
    ranks = rankdata(scores, method="average")
    return float((ranks[positives].sum() - positive_count * (positive_count + 1) / 2) / (positive_count * negative_count))


def grouped_select(x: np.ndarray, y: np.ndarray, groups: np.ndarray, classes: np.ndarray, penalties: list[float]) -> tuple[float, dict[str, Any]]:
    class_count = len(classes)
    by_penalty = {}
    for penalty in penalties:
        fold_probabilities = np.zeros((len(y), class_count), dtype=np.float64)
        fold_reports = []
        for group in sorted(set(groups.tolist())):
            validation = groups == group
            parameters, optimization = fit_softmax(x[~validation], y[~validation], class_count, penalty)
            fold_probabilities[validation] = probabilities(x[validation], parameters, class_count)
            fold_reports.append({"held_out_group": str(group), "cells": int(validation.sum()), "optimization": optimization})
        metrics = confusion_metrics(y, fold_probabilities, classes)
        by_penalty[str(penalty)] = {"metrics": metrics, "folds": fold_reports}
    selected = sorted(penalties, key=lambda penalty: (-by_penalty[str(penalty)]["metrics"]["macro_f1"], penalty))[0]
    return selected, {"selection_metric": "grouped_leave_one_donor_out_macro_f1", "selected_l2_coefficient": selected, "candidates": by_penalty}


def evaluate_task(
    x: np.ndarray, y: np.ndarray, roles: np.ndarray, groups: np.ndarray,
    classes: np.ndarray, penalties: list[float],
) -> tuple[np.ndarray, float, np.ndarray, dict[str, Any]]:
    fit = roles == "fit"
    calibration = roles == "calibration"
    sealed = roles == "sealed_test"
    selected, selection_report = grouped_select(x[fit], y[fit], groups[fit], classes, penalties)
    parameters, optimization = fit_softmax(x[fit], y[fit], len(classes), selected)
    calibration_logits = logits_from_parameters(x[calibration], parameters, len(classes))
    temperature, temperature_report = select_temperature(calibration_logits, y[calibration])
    calibration_probability = softmax_from_logits(calibration_logits / temperature)
    thresholds = conformal_thresholds(y[calibration], calibration_probability, len(classes))
    sealed_probability = probabilities(x[sealed], parameters, len(classes), temperature)
    sealed_metrics = confusion_metrics(y[sealed], sealed_probability, classes)
    sealed_metrics["ece_15_bin"] = ece(y[sealed], sealed_probability)
    sealed_metrics["conformal"] = conformal_metrics(y[sealed], sealed_probability, thresholds)
    return parameters, temperature, thresholds, {
        "selection": selection_report, "final_fit_optimization": optimization,
        "calibration": {"cells": int(calibration.sum()), "temperature": temperature, **temperature_report,
                        "class_conditional_nonconformity_thresholds": {str(classes[i]): float(value) for i, value in enumerate(thresholds)}},
        "sealed_internal_donor": sealed_metrics,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    tensor = np.load(args.tensor)
    x = tensor["train_X"]
    roles = tensor["train_role"]
    donors = tensor["train_donor"]
    penalties = [float(value) for value in protocol["model_contract"]["penalties"]]

    cell_y = encode(tensor["train_label"], CELL_CLASSES)
    cell_parameters, cell_temperature, cell_thresholds, cell_report = evaluate_task(
        x, cell_y, roles, donors, CELL_CLASSES, penalties
    )
    state_labels = np.where(tensor["train_condition"] == "ctrl", "control", "IFN_stimulated")
    state_y = encode(state_labels, STATE_CLASSES)
    state_parameters, state_temperature, state_thresholds, state_report = evaluate_task(
        x, state_y, roles, donors, STATE_CLASSES, penalties
    )

    external_ood = tensor["external_is_semantic_ood"].astype(bool)
    external_probability = probabilities(tensor["external_X"], cell_parameters, len(CELL_CLASSES), cell_temperature)
    external_in_domain = ~external_ood
    external_y = encode(tensor["external_label"][external_in_domain], CELL_CLASSES)
    external_in_probability = external_probability[external_in_domain]
    external_metrics = confusion_metrics(external_y, external_in_probability, CELL_CLASSES)
    external_metrics["ece_15_bin"] = ece(external_y, external_in_probability)
    external_metrics["conformal"] = conformal_metrics(external_y, external_in_probability, cell_thresholds)
    external_metrics["by_experiment_method"] = {}
    pair_values = list(zip(tensor["external_experiment"], tensor["external_method"], strict=True))
    for experiment, method in sorted(set(pair_values)):
        mask = external_in_domain & (tensor["external_experiment"] == experiment) & (tensor["external_method"] == method)
        if not np.any(mask):
            continue
        group_y = encode(tensor["external_label"][mask], CELL_CLASSES)
        group_probability = external_probability[mask]
        metrics = confusion_metrics(group_y, group_probability, CELL_CLASSES)
        metrics["ece_15_bin"] = ece(group_y, group_probability)
        metrics["conformal"] = conformal_metrics(group_y, group_probability, cell_thresholds)
        external_metrics["by_experiment_method"][f"{experiment}::{method}"] = metrics

    ood_scores = 1.0 - external_probability.max(axis=1)
    ood_report = {
        "positive_cells": int(external_ood.sum()), "negative_cells": int(external_in_domain.sum()),
        "score": "one_minus_maximum_calibrated_probability",
        "auroc": auroc(external_ood, ood_scores),
        "positive_author_labels": dict(Counter(tensor["external_author_label"][external_ood].tolist())),
    }
    endpoints = protocol["frozen_primary_cell_type_endpoints"]
    gates = {
        "sealed_internal_macro_f1": cell_report["sealed_internal_donor"]["macro_f1"] >= endpoints["sealed_internal_donor_macro_f1_minimum"],
        "external_macro_f1": external_metrics["macro_f1"] >= endpoints["external_lab_macro_f1_minimum"],
        "external_ece": external_metrics["ece_15_bin"] <= endpoints["external_lab_ece_maximum"],
        "external_conformal_coverage": external_metrics["conformal"]["marginal_coverage"] >= endpoints["external_lab_conformal_marginal_coverage_minimum"],
        "semantic_ood_auroc": ood_report["auroc"] >= protocol["OOD_contract"]["AUROC_minimum"],
    }

    args.model.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.model,
        feature_symbols=tensor["feature_symbols"], cell_classes=CELL_CLASSES,
        cell_parameters=cell_parameters, cell_temperature=np.asarray(cell_temperature),
        cell_conformal_thresholds=cell_thresholds, state_classes=STATE_CLASSES,
        state_parameters=state_parameters, state_temperature=np.asarray(state_temperature),
        state_conformal_thresholds=state_thresholds,
    )
    report = {
        "schema": "aleph.outer_library.pbmc_cross_lab_model.v1",
        "aleph_authority": "none",
        "protocol_sha256": sha256(args.protocol), "tensor_sha256": sha256(args.tensor),
        "implementation": {
            "classifier": "reference-class multinomial logistic regression",
            "objective": "mean cross-entropy plus 0.5*l2*sum(nonreference_weight_squared)",
            "calibration": "single scalar temperature on calibration donor",
            "biological_replication_warning": "cell rows are not independent biological replicates",
        },
        "cell_type": {"internal": cell_report, "pristine_external_lab": external_metrics, "semantic_ood": ood_report, "frozen_endpoint_passes": gates},
        "IFN_state": {"internal_donor_holdout": state_report, "cross_lab_direction_confirmation": "reported_separately_from_single_cell_classifier"},
        "overall_external_context_gate_pass": all(gates.values()),
        "authority": {
            "role": "pre_sweep_context_compression_and_candidate_ranking",
            "aleph_parameter_authority": "none", "may_emit_numeric_sweep_range": False,
            "why_no_numeric_range": "requires a validated Aleph observation operator, identifiability, and monotonic forward sensitivity",
        },
        "model_sha256": sha256(args.model),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "sealed_macro_f1": cell_report["sealed_internal_donor"]["macro_f1"],
        "external_macro_f1": external_metrics["macro_f1"], "external_ece": external_metrics["ece_15_bin"],
        "external_coverage": external_metrics["conformal"]["marginal_coverage"],
        "ood_auroc": ood_report["auroc"], "all_gates_pass": all(gates.values()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
