#!/usr/bin/env python3
"""Train a small raw-image MLP on leakage-safe BBBC048 cell splits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _softmax(logits: np.ndarray) -> np.ndarray:
    shifted = logits - logits.max(axis=1, keepdims=True)
    exp = np.exp(shifted)
    return exp / exp.sum(axis=1, keepdims=True)


def _predict(X: np.ndarray, params: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    hidden = np.maximum(0.0, X @ params["W1"] + params["b1"])
    logits = hidden @ params["W2"] + params["b2"]
    return hidden, logits


def _metrics(probability: np.ndarray, labels: np.ndarray) -> dict[str, object]:
    prediction = probability.argmax(axis=1)
    confusion = np.zeros((7, 7), dtype=np.int64)
    np.add.at(confusion, (labels, prediction), 1)
    f1, recall = [], []
    for label in range(7):
        tp = confusion[label, label]
        fp = confusion[:, label].sum() - tp
        fn = confusion[label].sum() - tp
        f1.append(float(2 * tp / max(1, 2 * tp + fp + fn)))
        recall.append(float(tp / max(1, tp + fn)))
    confidence = probability.max(axis=1)
    correct = prediction == labels
    ece = 0.0
    for lower in np.linspace(0.0, 0.9, 10):
        mask = (confidence >= lower) & (confidence < lower + 0.1 + 1e-12)
        if mask.any():
            ece += mask.mean() * abs(correct[mask].mean() - confidence[mask].mean())
    return {
        "accuracy": float(correct.mean()),
        "balanced_accuracy": float(np.mean(recall)),
        "macro_f1": float(np.mean(f1)),
        "per_class_f1": f1,
        "per_class_recall": recall,
        "ece_10_bin": float(ece),
        "confusion_matrix_true_by_predicted": confusion.tolist(),
    }


def train(tensor_path: Path, checkpoint_path: Path, *, epochs: int = 25) -> dict[str, object]:
    with np.load(tensor_path) as data:
        pixels = data["pixels"]
        labels = data["label"].astype(np.int64)
        split = data["split"]
        sample_ids = data["sample_id"]
    # Average non-overlapping 2x2 blocks: 3x16x16 -> 3x8x8 raw-image input.
    X = pixels.reshape(-1, 3, 8, 2, 8, 2).mean(axis=(3, 5)).reshape(-1, 192)
    train_mask = split == 0
    mean = X[train_mask].mean(axis=0)
    scale = X[train_mask].std(axis=0)
    scale[scale < 1e-6] = 1.0
    X = ((X - mean) / scale).astype(np.float32)
    rng = np.random.default_rng(48017)
    params = {
        "W1": (rng.standard_normal((192, 32)) * np.sqrt(2 / 192)).astype(np.float32),
        "b1": np.zeros(32, dtype=np.float32),
        "W2": (rng.standard_normal((32, 7)) * np.sqrt(2 / 32)).astype(np.float32),
        "b2": np.zeros(7, dtype=np.float32),
    }
    moments = {key: np.zeros_like(value) for key, value in params.items()}
    velocities = {key: np.zeros_like(value) for key, value in params.items()}
    train_indices = np.flatnonzero(train_mask)
    counts = np.bincount(labels[train_mask], minlength=7)
    class_weight = np.minimum(25.0, np.sqrt(train_indices.size / (7 * counts)))
    best_params = {key: value.copy() for key, value in params.items()}
    best_macro_f1 = -1.0
    best_epoch = 0
    step = 0
    for epoch in range(epochs):
        for batch_indices in np.array_split(rng.permutation(train_indices), 84):
            xb, yb = X[batch_indices], labels[batch_indices]
            hidden, logits = _predict(xb, params)
            probability = _softmax(logits)
            row_weight = class_weight[yb]
            dlogits = probability
            dlogits[np.arange(yb.size), yb] -= 1.0
            dlogits *= (row_weight / row_weight.sum())[:, None]
            gradients = {
                "W2": hidden.T @ dlogits + 1e-4 * params["W2"],
                "b2": dlogits.sum(axis=0),
            }
            dhidden = dlogits @ params["W2"].T
            dhidden[hidden <= 0] = 0
            gradients["W1"] = xb.T @ dhidden + 1e-4 * params["W1"]
            gradients["b1"] = dhidden.sum(axis=0)
            step += 1
            for key in params:
                moments[key] = 0.9 * moments[key] + 0.1 * gradients[key]
                velocities[key] = 0.999 * velocities[key] + 0.001 * gradients[key] ** 2
                mhat = moments[key] / (1 - 0.9 ** step)
                vhat = velocities[key] / (1 - 0.999 ** step)
                params[key] -= 1e-3 * mhat / (np.sqrt(vhat) + 1e-8)
        selection = split == 1
        selection_probability = _softmax(_predict(X[selection], params)[1])
        macro_f1 = float(_metrics(selection_probability, labels[selection])["macro_f1"])
        if macro_f1 > best_macro_f1:
            best_macro_f1 = macro_f1
            best_epoch = epoch + 1
            best_params = {key: value.copy() for key, value in params.items()}

    calibration = split == 2
    calibration_logits = _predict(X[calibration], best_params)[1]
    temperatures = np.linspace(0.5, 3.0, 101)
    nll = []
    for temperature in temperatures:
        probability = _softmax(calibration_logits / temperature)
        nll.append(float(-np.log(probability[np.arange(labels[calibration].size), labels[calibration]] + 1e-12).mean()))
    temperature = float(temperatures[int(np.argmin(nll))])
    calibration_probability = _softmax(calibration_logits / temperature)
    scores = 1.0 - calibration_probability[np.arange(labels[calibration].size), labels[calibration]]
    rank = min(scores.size - 1, int(np.ceil((scores.size + 1) * 0.90)) - 1)
    threshold = float(np.partition(scores, rank)[rank])

    test = split == 3
    test_probability = _softmax(_predict(X[test], best_params)[1] / temperature)
    test_metrics = _metrics(test_probability, labels[test])
    prediction_sets = test_probability >= (1.0 - threshold)
    test_coverage = float(prediction_sets[np.arange(labels[test].size), labels[test]].mean())
    mean_set_size = float(prediction_sets.sum(axis=1).mean())
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        checkpoint_path, **best_params, input_mean=mean, input_scale=scale,
        temperature=np.asarray(temperature), conformal_threshold=np.asarray(threshold),
    )
    return {
        "schema": "aleph.outer_library.bbbc048_mlp.v1",
        "model_family": "three_channel_8x8_raw_image_MLP_192_32_7",
        "optimizer": "Adam_class_weighted_cross_entropy",
        "epochs_requested": epochs,
        "best_selection_epoch": best_epoch,
        "best_selection_macro_f1": best_macro_f1,
        "temperature": temperature,
        "test": test_metrics,
        "conformal_90": {
            "calibration_cell_count": int(calibration.sum()),
            "test_marginal_coverage": test_coverage,
            "test_mean_set_size_of_7": mean_set_size,
            "rare_class_conditional_coverage_authority": False,
        },
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "status": "same_provider_training_only",
        "production_eligible": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "split_sample_id_overlap": 0,
        "notes": [
            "test is a held-out cell split from the same provider and acquisition lineage",
            "anaphase, telophase, and metaphase test counts are too small for conditional claims",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=35)
    args = parser.parse_args()
    report = train(args.tensor, args.checkpoint, epochs=args.epochs)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"test": report["test"], "conformal_90": report["conformal_90"]}, sort_keys=True))


if __name__ == "__main__":
    main()
