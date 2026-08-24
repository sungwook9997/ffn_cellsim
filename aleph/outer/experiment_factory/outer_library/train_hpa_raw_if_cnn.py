#!/usr/bin/env python3
"""Train the frozen two-channel HPA CNN before opening OpenCell labels."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image
from scipy.optimize import minimize_scalar
from scipy.stats import rankdata
import torch
from torch import nn


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def map_locations(values: list[str], protocol: dict[str, Any]) -> np.ndarray:
    result = np.zeros(len(protocol["coarse_ontology_in_order"]), dtype=np.float32)
    for value in values:
        lowered = value.casefold()
        for index, item in enumerate(protocol["coarse_ontology_in_order"]):
            if any(token in lowered for token in item["casefolded_substrings"]):
                result[index] = 1.0
    return result


def verify_source(
    manifest_path: Path, receipt_path: Path, image_dir: Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    by_name = {item["local_name"]: item for item in receipt["files"]}
    required = {
        info["local_name"]
        for row in manifest["rows"]
        for info in row["channels"].values()
    }
    if set(by_name) != required:
        raise ValueError("HPA raw-image receipt does not exactly cover the frozen manifest")
    for name, record in by_name.items():
        path = image_dir / name
        if path.stat().st_size != record["size_bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"HPA image differs from receipt: {name}")
    return manifest, receipt


def read_scaled(path: Path) -> np.ndarray:
    with Image.open(path) as image:
        value = np.asarray(image.convert("F"), dtype=np.float32)
    if value.ndim != 2 or not np.isfinite(value).all():
        raise ValueError(f"invalid HPA channel image: {path.name}")
    low, high = np.percentile(value, [1.0, 99.0])
    if not high > low:
        value = np.zeros_like(value)
    else:
        value = np.clip((value - low) / (high - low), 0.0, 1.0)
    side = min(value.shape)
    top = (value.shape[0] - side) // 2
    left = (value.shape[1] - side) // 2
    crop = np.ascontiguousarray(
        value[top:top + side, left:left + side], dtype=np.float32
    )
    resized = Image.fromarray(crop).resize(
        (128, 128), resample=Image.Resampling.BOX
    )
    result = np.array(resized, dtype=np.float32, copy=True)
    if not np.isfinite(result).all() or result.min() < 0.0 or result.max() > 1.0:
        raise ValueError(f"resized HPA channel left the finite [0,1] contract: {path.name}")
    return result


def build_tensor(
    manifest: dict[str, Any], image_dir: Path, protocol: dict[str, Any]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict[str, Any]]]:
    images, labels, splits, rows = [], [], [], []
    for row in manifest["rows"]:
        y = map_locations(row["author_locations"], protocol)
        if not np.any(y):
            continue
        nucleus = read_scaled(image_dir / row["channels"]["blue"]["local_name"])
        target = read_scaled(image_dir / row["channels"]["green"]["local_name"])
        images.append(np.stack((nucleus, target)))
        labels.append(y)
        splits.append(row["split"])
        rows.append(row)
    return (
        np.stack(images).astype(np.float32),
        np.stack(labels).astype(np.float32),
        np.asarray(splits),
        rows,
    )


class RawIFCNN(nn.Module):
    def __init__(self, class_count: int = 9) -> None:
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(2, 16, 3, padding=1), nn.GroupNorm(4, 16), nn.GELU(), nn.MaxPool2d(2),
            nn.Conv2d(16, 32, 3, padding=1), nn.GroupNorm(8, 32), nn.GELU(), nn.MaxPool2d(2),
            nn.Conv2d(32, 64, 3, padding=1), nn.GroupNorm(8, 64), nn.GELU(), nn.MaxPool2d(2),
            nn.Conv2d(64, 96, 3, padding=1), nn.GroupNorm(8, 96), nn.GELU(),
            nn.AdaptiveAvgPool2d(1),
        )
        self.classifier = nn.Linear(96, class_count)

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        return self.classifier(self.features(value).flatten(1))


def average_precision(y: np.ndarray, score: np.ndarray) -> float:
    positives = int(y.sum())
    if positives == 0:
        return float("nan")
    order = np.argsort(-score, kind="stable")
    ordered = y[order]
    precision = np.cumsum(ordered) / np.arange(1, len(ordered) + 1)
    return float((precision * ordered).sum() / positives)


def auroc(y: np.ndarray, score: np.ndarray) -> float:
    positive = y.astype(bool)
    p, n = int(positive.sum()), int((~positive).sum())
    if not p or not n:
        return float("nan")
    ranks = rankdata(score, method="average")
    return float((ranks[positive].sum() - p * (p + 1) / 2) / (p * n))


def thresholds(y: np.ndarray, probability: np.ndarray) -> np.ndarray:
    grid = np.linspace(0.05, 0.95, 37)
    result = np.zeros(y.shape[1], dtype=np.float64)
    for column in range(y.shape[1]):
        best: tuple[float, float, float] | None = None
        for threshold in grid:
            predicted = probability[:, column] >= threshold
            truth = y[:, column].astype(bool)
            tp = int((predicted & truth).sum())
            fp = int((predicted & ~truth).sum())
            fn = int((~predicted & truth).sum())
            f1 = 2 * tp / max(2 * tp + fp + fn, 1)
            candidate = (f1, -abs(float(threshold) - 0.5), float(threshold))
            if best is None or candidate > best:
                best = candidate
        result[column] = best[2]  # type: ignore[index]
    return result


def metrics(
    y: np.ndarray, probability: np.ndarray, cutoffs: np.ndarray
) -> dict[str, Any]:
    predicted = probability >= cutoffs
    tp = (predicted & (y == 1)).sum(axis=0)
    fp = (predicted & (y == 0)).sum(axis=0)
    fn = ((~predicted) & (y == 1)).sum(axis=0)
    f1 = 2 * tp / np.maximum(2 * tp + fp + fn, 1)
    eligible = (y.sum(axis=0) > 0) & (y.sum(axis=0) < len(y))
    aps = np.asarray([average_precision(y[:, i], probability[:, i]) for i in range(y.shape[1])])
    aucs = np.asarray([auroc(y[:, i], probability[:, i]) for i in range(y.shape[1])])
    return {
        "sample_count": int(len(y)),
        "macro_f1": float(np.mean(f1[eligible])),
        "micro_f1": float(2 * tp.sum() / max(2 * tp.sum() + fp.sum() + fn.sum(), 1)),
        "macro_average_precision": float(np.nanmean(aps[eligible])),
        "macro_auroc": float(np.nanmean(aucs[eligible])),
        "exact_match_accuracy": float(np.all(predicted == y, axis=1).mean()),
        "per_class": [
            {
                "support": int(y[:, i].sum()), "f1": float(f1[i]),
                "average_precision": float(aps[i]) if np.isfinite(aps[i]) else None,
                "auroc": float(aucs[i]) if np.isfinite(aucs[i]) else None,
            }
            for i in range(y.shape[1])
        ],
    }


def multilabel_ece(y: np.ndarray, probability: np.ndarray, bins: int = 15) -> float:
    truth = y.ravel()
    score = probability.ravel()
    total = 0.0
    for lower in np.linspace(0.0, 1.0, bins + 1)[:-1]:
        upper = lower + 1.0 / bins
        mask = (score >= lower) & (score < upper if upper < 1 else score <= upper)
        if np.any(mask):
            total += float(mask.mean()) * abs(float(truth[mask].mean()) - float(score[mask].mean()))
    return total


def sigmoid(logits: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(logits, -50, 50)))


def fit_temperature(logits: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    def objective(log_t: float) -> float:
        p = sigmoid(logits / math.exp(log_t))
        return float(-np.mean(y * np.log(np.maximum(p, 1e-12)) + (1 - y) * np.log(np.maximum(1 - p, 1e-12))))
    fit = minimize_scalar(objective, bounds=(math.log(0.05), math.log(20.0)), method="bounded")
    return float(math.exp(fit.x)), float(fit.fun)


def conformal_thresholds(
    y: np.ndarray, probability: np.ndarray, minimum: int, alpha: float = 0.1
) -> tuple[np.ndarray, np.ndarray]:
    result = np.empty(y.shape[1], dtype=np.float64)
    eligible = np.ones(y.shape[1], dtype=bool)
    for column in range(y.shape[1]):
        score = 1.0 - probability[y[:, column] == 1, column]
        if len(score) < minimum:
            result[column] = 1.0
            eligible[column] = False
            continue
        rank = min(int(np.ceil((len(score) + 1) * (1.0 - alpha))), len(score))
        result[column] = np.partition(score, rank - 1)[rank - 1]
    return result, eligible


def predict(model: nn.Module, x: torch.Tensor, batch_size: int = 32) -> np.ndarray:
    model.eval()
    output = []
    with torch.no_grad():
        for start in range(0, len(x), batch_size):
            output.append(model(x[start:start + batch_size]).cpu().numpy())
    return np.concatenate(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol["frozen_before_opencell_annotation_or_image_content_access"]:
        raise ValueError("OpenCell protocol was not frozen before source access")
    manifest, receipt = verify_source(args.manifest, args.receipt, args.image_dir)
    x_array, y_array, split, rows = build_tensor(manifest, args.image_dir, protocol)
    if not np.isfinite(x_array).all() or x_array.min() < 0.0 or x_array.max() > 1.0:
        raise ValueError("HPA image tensor violates the finite [0,1] contract")
    classes = np.asarray([item["canonical"] for item in protocol["coarse_ontology_in_order"]])
    masks = {name: split == name for name in ("train", "validation", "calibration", "test")}
    if any(not np.any(mask) for mask in masks.values()):
        raise ValueError("mapped HPA tensor has an empty frozen split")
    if np.any(y_array[masks["train"]].sum(axis=0) == 0):
        raise ValueError("HPA training split lacks a frozen coarse class")

    seed = int(protocol["model_contract"]["seed"])
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(min(8, max(1, torch.get_num_threads())))
    x = torch.from_numpy(x_array)
    y = torch.from_numpy(y_array)
    model = RawIFCNN(len(classes))
    train_indices = np.flatnonzero(masks["train"])
    positive = y_array[masks["train"]].sum(axis=0)
    negative = masks["train"].sum() - positive
    positive_weight = np.clip(negative / np.maximum(positive, 1), 1.0, 20.0)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.from_numpy(positive_weight.astype(np.float32)))
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-4)

    best_score = -float("inf")
    best_epoch = -1
    best_state: dict[str, torch.Tensor] | None = None
    history = []
    patience = 0
    for epoch in range(int(protocol["model_contract"]["maximum_epochs"])):
        model.train()
        rng = np.random.default_rng(seed + epoch)
        order = rng.permutation(train_indices)
        epoch_loss = 0.0
        seen = 0
        for start in range(0, len(order), 16):
            indices = order[start:start + 16]
            batch = x[indices].clone()
            for position in range(len(batch)):
                rotation = int(rng.integers(0, 4))
                if rotation:
                    batch[position] = torch.rot90(batch[position], rotation, dims=(1, 2))
                if bool(rng.integers(0, 2)):
                    batch[position] = torch.flip(batch[position], dims=(2,))
            target = y[indices]
            optimizer.zero_grad(set_to_none=True)
            logits = model(batch)
            loss = loss_fn(logits, target)
            loss.backward()
            optimizer.step()
            epoch_loss += float(loss.detach()) * len(indices)
            seen += len(indices)
        validation_logits = predict(model, x[masks["validation"]])
        validation_probability = sigmoid(validation_logits)
        aps = [
            average_precision(y_array[masks["validation"], i], validation_probability[:, i])
            for i in range(len(classes))
        ]
        score = float(np.nanmean(aps))
        history.append({"epoch": epoch + 1, "train_loss": epoch_loss / seen, "validation_macro_average_precision": score})
        if score > best_score + 1e-12:
            best_score = score
            best_epoch = epoch + 1
            best_state = {name: value.detach().cpu().clone() for name, value in model.state_dict().items()}
            patience = 0
        else:
            patience += 1
        if patience >= 15:
            break
    if best_state is None:
        raise AssertionError("training produced no checkpoint")
    model.load_state_dict(best_state)

    all_logits = predict(model, x)
    validation_probability = sigmoid(all_logits[masks["validation"]])
    decision_thresholds = thresholds(y_array[masks["validation"]], validation_probability)
    temperature, calibration_nll = fit_temperature(
        all_logits[masks["calibration"]], y_array[masks["calibration"]]
    )
    probability = sigmoid(all_logits / temperature)
    conformal, conformal_eligible = conformal_thresholds(
        y_array[masks["calibration"]], probability[masks["calibration"]],
        int(protocol["uncertainty_contract"]["HPA_calibration_positive_minimum_per_class"]),
    )
    sealed = metrics(y_array[masks["test"]], probability[masks["test"]], decision_thresholds)
    sealed["multilabel_ece_15_bin"] = multilabel_ece(
        y_array[masks["test"]], probability[masks["test"]]
    )
    included = (1.0 - probability[masks["test"]]) <= conformal
    true_positive = y_array[masks["test"]] == 1
    sealed["positive_label_conformal_coverage"] = float(included[true_positive].mean())
    sealed["mean_prediction_set_size"] = float(included.sum(axis=1).mean())

    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save({
        "state_dict": best_state,
        "classes": classes.tolist(),
        "decision_thresholds": decision_thresholds.tolist(),
        "temperature": temperature,
        "conformal_thresholds": conformal.tolist(),
        "conformal_class_eligible": conformal_eligible.tolist(),
        "input_shape": [2, 128, 128],
        "seed": seed,
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "receipt_sha256": sha256(args.receipt),
    }, args.checkpoint)
    report = {
        "schema": "aleph.outer_library.hpa_raw_if_cnn_training.v1",
        "trained_before_opencell_annotation_or_image_content_access": True,
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "receipt_sha256": sha256(args.receipt),
        "checkpoint_sha256": sha256(args.checkpoint),
        "software": {"torch": torch.__version__, "numpy": np.__version__},
        "classes": classes.tolist(),
        "mapped_image_count": len(rows),
        "split_counts": {name: int(mask.sum()) for name, mask in masks.items()},
        "split_positive_counts": {
            name: {classes[i]: int(y_array[mask, i].sum()) for i in range(len(classes))}
            for name, mask in masks.items()
        },
        "forbidden_HPA_channels_used": False,
        "architecture": protocol["model_contract"]["architecture"],
        "training": {
            "selected_epoch": best_epoch,
            "epochs_executed": len(history),
            "selected_validation_macro_average_precision": best_score,
            "positive_weights": {classes[i]: float(positive_weight[i]) for i in range(len(classes))},
            "history": history,
        },
        "calibration": {
            "temperature": temperature,
            "negative_log_likelihood": calibration_nll,
            "decision_thresholds": {classes[i]: float(decision_thresholds[i]) for i in range(len(classes))},
            "conformal_thresholds": {classes[i]: float(conformal[i]) for i in range(len(classes))},
            "class_eligible": {classes[i]: bool(conformal_eligible[i]) for i in range(len(classes))},
        },
        "HPA_sealed_test": sealed,
        "frozen_internal_endpoint_passed": sealed["macro_average_precision"]
        >= protocol["frozen_endpoints"]["HPA_sealed_test_macro_average_precision_minimum"],
        "OpenCell_accessed": False,
        "external_holdout_status": "unopened",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "mapped_images": len(rows), "selected_epoch": best_epoch,
        "validation_macro_ap": best_score,
        "sealed_test_macro_ap": sealed["macro_average_precision"],
        "sealed_test_macro_f1": sealed["macro_f1"],
        "internal_endpoint_passed": report["frozen_internal_endpoint_passed"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
