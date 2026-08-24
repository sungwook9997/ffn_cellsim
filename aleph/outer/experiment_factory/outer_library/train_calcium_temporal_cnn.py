#!/usr/bin/env python3
"""Train a small dependency-light temporal CNN and test cross-lab transfer.

The convolution filters and classifier are trained jointly with L-BFGS on HeLa
cells. Architecture and regularization are selected using HeLa validation cells
only. EA.hy926 is evaluated after selection and is explicitly retrospective,
because that dataset had already been opened by earlier baseline revisions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize

from evaluate_piezo1_transfer import ACTIVE_NANOSWITCH_STRATA, _metrics
from train_calcium_baseline import binary_auroc


SCHEMA = "aleph.outer_library.calcium_temporal_cnn.v1"
SEED = 20260805
LENGTH = 128
KERNEL_SIZE = 9


def preprocess_curve(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64).copy()
    values -= values[0]
    amplitude = max(
        float(np.max(np.abs(values))), float(np.ptp(values)), 1e-8
    )
    values /= amplitude
    return np.interp(
        np.linspace(0.0, 1.0, LENGTH),
        np.linspace(0.0, 1.0, len(values)),
        values,
    )


def _curves(
    path: Path,
    source: str,
) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if source == "hela" and row["observable"] != (
                "intracellular_calcium_F_over_F0"
            ):
                continue
            grouped[row["sample_id"]].append(row)
    curves, labels, sample_ids, conditions = [], [], [], []
    for sample_id, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: row["coordinate_value"])
        condition = rows[0]["condition"]
        if source == "hela":
            if condition == "TRPV4_inhibited_HC-067047":
                continue
            label = int(condition == "Piezo1_inhibited_GsMTx4")
        else:
            if not condition.startswith(ACTIVE_NANOSWITCH_STRATA):
                continue
            if "+GsMTx4" in condition:
                label = 1
            elif "+Mag" in condition:
                label = 0
            else:
                continue
        values = np.asarray([row["value"] for row in rows], dtype=np.float64)
        curves.append(preprocess_curve(values))
        labels.append(label)
        sample_ids.append(sample_id)
        conditions.append(condition)
    if not curves:
        raise ValueError(f"no aligned calcium curves in {path}")
    return (
        np.stack(curves), np.asarray(labels, dtype=np.int64),
        sample_ids, conditions,
    )


def _split(labels: np.ndarray, sample_ids: list[str]) -> dict[str, np.ndarray]:
    parts: dict[str, list[int]] = {"train": [], "validation": []}
    for label in (0, 1):
        indices = np.flatnonzero(labels == label).tolist()
        indices.sort(key=lambda index: hashlib.sha256(
            sample_ids[index].encode()
        ).hexdigest())
        if len(indices) != 50:
            raise ValueError("expected exactly 50 aligned HeLa cells per class")
        parts["train"].extend(indices[:30])
        parts["validation"].extend(indices[30:])
    return {
        name: np.asarray(sorted(indices), dtype=np.int64)
        for name, indices in parts.items()
    }


def _unpack(
    vector: np.ndarray,
    filter_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    offset = 0
    count = filter_count * KERNEL_SIZE
    filters = vector[offset:offset + count].reshape(filter_count, KERNEL_SIZE)
    offset += count
    bias = vector[offset:offset + filter_count]
    offset += filter_count
    count = 4 * filter_count
    classifier = vector[offset:offset + count].reshape(2 * filter_count, 2)
    offset += count
    intercept = vector[offset:offset + 2]
    return filters, bias, classifier, intercept


def _loss_gradient(
    vector: np.ndarray,
    windows: np.ndarray,
    labels: np.ndarray,
    filter_count: int,
    regularization: float,
) -> tuple[float, np.ndarray]:
    filters, bias, classifier, intercept = _unpack(vector, filter_count)
    convolution = windows @ filters.T + bias
    activation = np.maximum(convolution, 0.0)
    mean_pool = activation.mean(axis=1)
    maximum_index = activation.argmax(axis=1)
    maximum_pool = activation[
        np.arange(len(activation))[:, None],
        maximum_index,
        np.arange(filter_count)[None, :],
    ]
    embedding = np.column_stack((mean_pool, maximum_pool))
    raw_logits = embedding @ classifier + intercept
    raw_logits -= raw_logits.max(axis=1, keepdims=True)
    exponent = np.exp(raw_logits)
    probability = exponent / exponent.sum(axis=1, keepdims=True)
    loss = float(-np.log(np.clip(
        probability[np.arange(len(labels)), labels], 1e-12, 1.0
    )).mean())
    loss += regularization * float(
        np.sum(filters * filters) + np.sum(classifier * classifier)
    ) / 2.0

    delta_logits = probability
    delta_logits[np.arange(len(labels)), labels] -= 1.0
    delta_logits /= len(labels)
    delta_classifier = embedding.T @ delta_logits + regularization * classifier
    delta_intercept = delta_logits.sum(axis=0)
    delta_embedding = delta_logits @ classifier.T
    delta_activation = np.zeros_like(activation)
    delta_activation += delta_embedding[:, :filter_count, None].transpose(0, 2, 1) / (
        activation.shape[1]
    )
    delta_activation[
        np.arange(len(activation))[:, None],
        maximum_index,
        np.arange(filter_count)[None, :],
    ] += delta_embedding[:, filter_count:]
    delta_convolution = delta_activation * (convolution > 0.0)
    delta_filters = np.einsum(
        "ntf,ntk->fk", delta_convolution, windows
    ) + regularization * filters
    delta_bias = delta_convolution.sum(axis=(0, 1))
    gradient = np.concatenate((
        delta_filters.ravel(), delta_bias,
        delta_classifier.ravel(), delta_intercept,
    ))
    return loss, gradient


def _fit(
    curves: np.ndarray,
    labels: np.ndarray,
    filter_count: int,
    regularization: float,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    windows = np.lib.stride_tricks.sliding_window_view(
        curves, KERNEL_SIZE, axis=1
    )
    parameter_count = filter_count * KERNEL_SIZE + filter_count + 4 * filter_count + 2
    initial = np.random.default_rng(SEED + filter_count).normal(
        0.0, 0.08, parameter_count
    )
    result = minimize(
        _loss_gradient,
        initial,
        args=(windows, labels, filter_count, regularization),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": 300, "ftol": 1e-10},
    )
    filters, bias, classifier, intercept = _unpack(result.x, filter_count)
    return {
        "filters": filters,
        "bias": bias,
        "classifier": classifier,
        "intercept": intercept,
    }, {
        "converged": bool(result.success),
        "iterations": int(result.nit),
        "final_loss": float(result.fun),
        "message": str(result.message),
    }


def _embedding(model: dict[str, np.ndarray], curves: np.ndarray) -> np.ndarray:
    windows = np.lib.stride_tricks.sliding_window_view(
        curves, KERNEL_SIZE, axis=1
    )
    activation = np.maximum(
        windows @ model["filters"].T + model["bias"], 0.0
    )
    return np.column_stack((activation.mean(axis=1), activation.max(axis=1)))


def _probability(model: dict[str, np.ndarray], curves: np.ndarray) -> np.ndarray:
    raw_logits = _embedding(model, curves) @ model["classifier"] + model["intercept"]
    raw_logits -= raw_logits.max(axis=1, keepdims=True)
    exponent = np.exp(raw_logits)
    return exponent / exponent.sum(axis=1, keepdims=True)


def _checkpoint_sha256(path: Path, model: dict[str, np.ndarray]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        schema=np.asarray(SCHEMA),
        filters=model["filters"],
        bias=model["bias"],
        classifier=model["classifier"],
        intercept=model["intercept"],
        embedding_mean=model["embedding_mean"],
        embedding_scale=model["embedding_scale"],
        class_centroids=model["class_centroids"],
        class_names=np.asarray(["non_inhibited_reference", "Piezo1_inhibited"]),
        input_length=np.asarray(LENGTH, dtype=np.int64),
        kernel_size=np.asarray(KERNEL_SIZE, dtype=np.int64),
    )
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _unlabelled_label_shift(
    probabilities: np.ndarray,
    source_prior: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Saerens-style EM using target predictions but never target labels."""
    target_prior = source_prior.astype(np.float64).copy()
    adjusted = probabilities.copy()
    for iteration in range(1_000):
        adjusted = probabilities * (target_prior / source_prior)[None, :]
        adjusted /= adjusted.sum(axis=1, keepdims=True)
        updated = adjusted.mean(axis=0)
        if float(np.max(np.abs(updated - target_prior))) < 1e-12:
            return updated, adjusted, iteration + 1
        target_prior = updated
    return target_prior, adjusted, 1_000


def train(
    train_manifest: Path,
    external_manifest: Path,
    checkpoint: Path,
) -> dict[str, Any]:
    curves, labels, sample_ids, train_conditions = _curves(
        train_manifest, "hela"
    )
    external_curves, external_labels, external_ids, external_conditions = _curves(
        external_manifest, "external"
    )
    if set(sample_ids) & set(external_ids):
        raise ValueError("train and external sample IDs overlap")
    split = _split(labels, sample_ids)
    candidates = []
    best = None
    for filter_count in (4, 8):
        for regularization in (0.01, 0.1):
            model, optimizer = _fit(
                curves[split["train"]], labels[split["train"]],
                filter_count, regularization,
            )
            validation = _metrics(
                _probability(model, curves[split["validation"]]),
                labels[split["validation"]],
            )
            candidate = {
                "filter_count": filter_count,
                "regularization": regularization,
                "validation_metrics": validation,
                "optimizer": optimizer,
            }
            candidates.append(candidate)
            key = (
                validation["balanced_accuracy"], validation["macro_f1"],
                validation["auroc"], -filter_count, -regularization,
            )
            if best is None or key > best[0]:
                best = (key, model, candidate)
    assert best is not None
    model, selected = best[1], best[2]
    external_probability = _probability(model, external_curves)
    external_metrics = _metrics(external_probability, external_labels)
    source_prior = np.bincount(
        labels[split["train"]], minlength=2
    ).astype(np.float64)
    source_prior /= source_prior.sum()
    estimated_prior, shifted_probability, shift_iterations = _unlabelled_label_shift(
        external_probability, source_prior
    )
    shifted_metrics = _metrics(shifted_probability, external_labels)
    actual_external_prior = np.bincount(
        external_labels, minlength=2
    ).astype(np.float64)
    actual_external_prior /= actual_external_prior.sum()

    training_embedding = _embedding(model, curves[split["train"]])
    validation_embedding = _embedding(model, curves[split["validation"]])
    external_embedding = _embedding(model, external_curves)
    mean = training_embedding.mean(axis=0)
    scale = training_embedding.std(axis=0)
    scale[scale < 1e-12] = 1.0
    train_z = (training_embedding - mean) / scale
    validation_z = (validation_embedding - mean) / scale
    external_z = (external_embedding - mean) / scale
    centroids = np.stack([
        train_z[labels[split["train"]] == label].mean(axis=0)
        for label in (0, 1)
    ])
    validation_distance = (
        (validation_z[:, None] - centroids[None]) ** 2
    ).mean(axis=2).min(axis=1)
    external_distance = (
        (external_z[:, None] - centroids[None]) ** 2
    ).mean(axis=2).min(axis=1)
    ood_auroc = binary_auroc(validation_distance, external_distance)
    model["embedding_mean"] = mean
    model["embedding_scale"] = scale
    model["class_centroids"] = centroids
    checkpoint_sha256 = _checkpoint_sha256(checkpoint, model)
    passed = all(external_metrics[name] >= 0.70 for name in (
        "balanced_accuracy", "macro_f1", "auroc",
    ))
    return {
        "schema": SCHEMA,
        "task": "binary_Piezo1_inhibitor_signature_temporal_CNN_transfer",
        "architecture": {
            "input": "baseline_subtracted_amplitude_normalized_curve_resampled_to_128",
            "temporal_convolution": {
                "filter_count": selected["filter_count"],
                "kernel_size": KERNEL_SIZE,
                "activation": "ReLU",
            },
            "pooling": ["global_mean", "global_max"],
            "head": "two_class_linear_softmax",
            "parameter_count": int(
                selected["filter_count"] * KERNEL_SIZE
                + selected["filter_count"]
                + 4 * selected["filter_count"] + 2
            ),
            "filters_are_trainable": True,
            "optimizer": "scipy_L-BFGS-B_exact_gradient",
            "gpu_required": False,
        },
        "selection_protocol": (
            "filter count and L2 regularization selected only on deterministic "
            "HeLa train/validation cells; EA.hy926 labels never select architecture, "
            "weights, threshold, or checkpoint"
        ),
        "split": {
            "train_cells": int(len(split["train"])),
            "validation_cells": int(len(split["validation"])),
            "split_unit": "cell_sample_id",
            "training_dataset": "zenodo-18495219 HeLa",
            "external_dataset": "zenodo-19890985 EA.hy926",
            "external_dataset_and_lab_holdout": True,
            "sample_overlap_count": 0,
            "external_was_opened_before_current_architecture": True,
            "pristine_external_confirmation": False,
        },
        "candidate_selection": candidates,
        "selected_validation_metrics": selected["validation_metrics"],
        "external_retrospective_metrics": external_metrics,
        "external_embedding_ood": {
            "score": "nearest_training_class_centroid_standardized_embedding_distance",
            "external_vs_validation_auroc": ood_auroc,
            "ood_detected": ood_auroc >= 0.70,
        },
        "unlabelled_label_shift_diagnostic": {
            "method": "Saerens_EM_target_prior_from_unlabelled_probabilities",
            "source_prior": source_prior.tolist(),
            "estimated_external_prior_without_labels": estimated_prior.tolist(),
            "iterations": shift_iterations,
            "retrospective_actual_external_prior": actual_external_prior.tolist(),
            "actual_prior_used_for_fitting_or_selection": False,
            "adjusted_external_metrics": shifted_metrics,
            "adaptation_accepted": bool(
                shifted_metrics["macro_f1"] > external_metrics["macro_f1"]
                and shifted_metrics["balanced_accuracy"]
                > external_metrics["balanced_accuracy"]
            ),
            "interpretation": (
                "label-shift assumptions fail: unlabelled EM moves the prior in "
                "the wrong direction and does not repair conditional cross-lab shift"
            ),
        },
        "checkpoint": {
            "filename": checkpoint.name,
            "sha256": checkpoint_sha256,
            "format": "numpy_npz_no_pickle",
            "contains": [
                "trainable_filters_and_classifier",
                "preprocessing_dimensions",
                "training_embedding_standardization",
                "training_class_centroids_for_OOD_distance",
                "class_names",
            ],
        },
        "task_holdout_passed": passed,
        "architecture_accepted": passed,
        "uncertainty_holdout_passed": False,
        "production_eligible": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "EA.hy926 was already inspected by earlier model revisions, so this is retrospective",
            "external biological replicate identities are unavailable per cell",
            "the CNN has no disjoint calibration partition; uncertainty authority remains with the existing conformal audit",
            "a trainable temporal encoder does not by itself remove cross-lab conditional shift",
            "nearest-centroid embedding distance does not detect this failed transfer, so OOD distance cannot substitute for task labels",
            "unlabelled label-shift EM moves the inferred target prior in the wrong direction and is rejected",
            "a new aligned lab is required before this architecture can be confirmatory",
        ],
        "conditions": {
            "training": sorted(set(train_conditions)),
            "external": sorted(set(external_conditions)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train-manifest", type=Path, required=True)
    parser.add_argument("--external-manifest", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = train(args.train_manifest, args.external_manifest, args.checkpoint)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "selected_validation_metrics": report["selected_validation_metrics"],
        "external_retrospective_metrics": report["external_retrospective_metrics"],
        "architecture_accepted": report["architecture_accepted"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
