#!/usr/bin/env python3
"""Train the frozen HPA/OpenCell multiprovider raw-IF model."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import rankdata
import torch
from torch import nn
from torch.nn import functional as F

from train_hpa_raw_if_cnn import sha256


def group_norm(channels: int) -> nn.GroupNorm:
    groups = max(value for value in range(1, 9) if channels % value == 0)
    return nn.GroupNorm(groups, channels)


class ResidualBlock(nn.Module):
    def __init__(self, channels: int) -> None:
        super().__init__()
        self.conv1 = nn.Conv2d(channels, channels, 3, padding=1, bias=False)
        self.norm1 = group_norm(channels)
        self.conv2 = nn.Conv2d(channels, channels, 3, padding=1, bias=False)
        self.norm2 = group_norm(channels)

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        residual = value
        value = F.gelu(self.norm1(self.conv1(value)))
        value = self.norm2(self.conv2(value))
        return F.gelu(value + residual)


class MultiproviderIFCNN(nn.Module):
    def __init__(self, class_count: int = 9) -> None:
        super().__init__()
        self.stem = nn.Sequential(
            nn.Conv2d(2, 16, 3, padding=1, bias=False), group_norm(16), nn.GELU()
        )
        self.stage1 = nn.Sequential(ResidualBlock(16), nn.MaxPool2d(2))
        self.stage2 = nn.Sequential(
            nn.Conv2d(16, 32, 3, padding=1, bias=False), group_norm(32), nn.GELU(),
            ResidualBlock(32), nn.MaxPool2d(2),
        )
        self.stage3 = nn.Sequential(
            nn.Conv2d(32, 64, 3, padding=1, bias=False), group_norm(64), nn.GELU(),
            ResidualBlock(64), nn.MaxPool2d(2),
        )
        self.stage4 = nn.Sequential(
            nn.Conv2d(64, 96, 3, padding=1, bias=False), group_norm(96), nn.GELU(),
            ResidualBlock(96), nn.AdaptiveAvgPool2d(1),
        )
        self.embedding = nn.Sequential(nn.Linear(96, 64), nn.GELU())
        self.classifier = nn.Linear(64, class_count)

    def forward(self, value: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        value = self.stem(value)
        value = self.stage1(value)
        value = self.stage2(value)
        value = self.stage3(value)
        value = self.stage4(value).flatten(1)
        embedding = self.embedding(value)
        return self.classifier(embedding), embedding


def sigmoid(value: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(value, -50.0, 50.0)))


def average_precision(truth: np.ndarray, score: np.ndarray) -> float:
    positives = int(truth.sum())
    if positives == 0:
        return float("nan")
    order = np.argsort(-score, kind="stable")
    ordered = truth[order]
    precision = np.cumsum(ordered) / np.arange(1, len(ordered) + 1)
    return float((precision * ordered).sum() / positives)


def auroc(truth: np.ndarray, score: np.ndarray) -> float:
    positive = truth.astype(bool)
    p, n = int(positive.sum()), int((~positive).sum())
    if not p or not n:
        return float("nan")
    ranks = rankdata(score, method="average")
    return float((ranks[positive].sum() - p * (p + 1) / 2) / (p * n))


def aggregate(
    logits: np.ndarray, embeddings: np.ndarray, labels: np.ndarray, units: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    names = np.asarray(sorted(set(units.astype(str))))
    out_logits, out_embeddings, out_labels = [], [], []
    unit_strings = units.astype(str)
    for name in names:
        mask = unit_strings == name
        group_labels = labels[mask]
        if not np.all(group_labels == group_labels[0]):
            raise ValueError(f"nested projections disagree on labels for {name}")
        out_logits.append(logits[mask].mean(axis=0))
        normalized = embeddings[mask] / np.maximum(
            np.linalg.norm(embeddings[mask], axis=1, keepdims=True), 1e-12
        )
        mean_embedding = normalized.mean(axis=0)
        mean_embedding /= max(np.linalg.norm(mean_embedding), 1e-12)
        out_embeddings.append(mean_embedding)
        out_labels.append(group_labels[0])
    return (
        np.stack(out_logits), np.stack(out_embeddings), np.stack(out_labels), names
    )


def predict(
    model: nn.Module, images: np.ndarray, batch_size: int = 64
) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    logits, embeddings = [], []
    with torch.no_grad():
        for start in range(0, len(images), batch_size):
            batch = torch.from_numpy(images[start : start + batch_size]).float().div_(255.0)
            batch_logits, batch_embeddings = model(batch)
            logits.append(batch_logits.cpu().numpy())
            embeddings.append(batch_embeddings.cpu().numpy())
    return np.concatenate(logits), np.concatenate(embeddings)


def provider_selection_metrics(
    labels: np.ndarray, probability: np.ndarray, minimum_positive: int = 5
) -> dict[str, Any]:
    per_class = []
    supported_values = []
    for index in range(labels.shape[1]):
        truth = labels[:, index]
        positive = int(truth.sum())
        negative = int(len(truth) - positive)
        prevalence = float(truth.mean())
        ap = average_precision(truth, probability[:, index])
        supported = positive >= minimum_positive and negative >= minimum_positive
        normalized = (
            (ap - prevalence) / max(1.0 - prevalence, 1e-12) if supported else None
        )
        if supported:
            supported_values.append((ap, normalized, prevalence))
        per_class.append({
            "positive_units": positive,
            "negative_units": negative,
            "prevalence": prevalence,
            "average_precision": ap,
            "prevalence_normalized_AP": normalized,
            "supported": supported,
        })
    if not supported_values:
        raise ValueError("provider selection has no supported class")
    return {
        "unit_count": int(len(labels)),
        "supported_class_count": len(supported_values),
        "macro_average_precision": float(np.mean([item[0] for item in supported_values])),
        "prevalence_normalized_AP": float(np.mean([item[1] for item in supported_values])),
        "every_supported_class_AP_above_prevalence": bool(all(
            item[0] > item[2] for item in supported_values
        )),
        "per_class": per_class,
    }


def f1_at(truth: np.ndarray, score: np.ndarray, threshold: float) -> float:
    predicted = score >= threshold
    positive = truth.astype(bool)
    tp = int((predicted & positive).sum())
    fp = int((predicted & ~positive).sum())
    fn = int((~predicted & positive).sum())
    return float(2 * tp / max(2 * tp + fp + fn, 1))


def decision_thresholds(
    provider_data: list[tuple[np.ndarray, np.ndarray]], minimum: int = 5
) -> np.ndarray:
    grid = np.linspace(0.05, 0.95, 37)
    result = np.empty(provider_data[0][0].shape[1], dtype=np.float64)
    for column in range(len(result)):
        eligible = [
            (labels[:, column], probability[:, column])
            for labels, probability in provider_data
            if labels[:, column].sum() >= minimum
            and (1 - labels[:, column]).sum() >= minimum
        ]
        if not eligible:
            result[column] = 0.5
            continue
        candidates = []
        for threshold in grid:
            minimum_f1 = min(f1_at(y, p, float(threshold)) for y, p in eligible)
            candidates.append((minimum_f1, -abs(float(threshold) - 0.5), -float(threshold)))
        result[column] = -max(candidates)[2]
    return result


def multilabel_ece(labels: np.ndarray, probability: np.ndarray, bins: int = 15) -> float:
    truth, score = labels.ravel(), probability.ravel()
    total = 0.0
    for index in range(bins):
        lower, upper = index / bins, (index + 1) / bins
        mask = (score >= lower) & (score < upper if index + 1 < bins else score <= upper)
        if np.any(mask):
            total += float(mask.mean()) * abs(float(truth[mask].mean()) - float(score[mask].mean()))
    return total


def fit_temperature(provider_data: list[tuple[np.ndarray, np.ndarray]]) -> tuple[float, float]:
    def objective(log_temperature: float) -> float:
        temperature = math.exp(log_temperature)
        losses = []
        for logits, labels in provider_data:
            probability = sigmoid(logits / temperature)
            loss = -np.mean(
                labels * np.log(np.maximum(probability, 1e-12))
                + (1 - labels) * np.log(np.maximum(1 - probability, 1e-12))
            )
            losses.append(float(loss))
        return float(np.mean(losses))
    fit = minimize_scalar(
        objective, bounds=(math.log(0.05), math.log(20.0)), method="bounded"
    )
    return float(math.exp(fit.x)), float(fit.fun)


def positive_conformal(
    provider_data: dict[str, tuple[np.ndarray, np.ndarray]], minimum: int, alpha: float
) -> tuple[np.ndarray, dict[str, Any]]:
    class_count = next(iter(provider_data.values()))[0].shape[1]
    threshold = np.ones(class_count, dtype=np.float64)
    details: dict[str, Any] = {}
    for column in range(class_count):
        provider_quantiles = {}
        for provider, (labels, probability) in provider_data.items():
            scores = 1.0 - probability[labels[:, column] == 1, column]
            if len(scores) < minimum:
                provider_quantiles[provider] = None
                continue
            rank = min(int(np.ceil((len(scores) + 1) * (1.0 - alpha))), len(scores))
            provider_quantiles[provider] = float(np.partition(scores, rank - 1)[rank - 1])
        eligible = [value for value in provider_quantiles.values() if value is not None]
        if eligible:
            threshold[column] = max(eligible)
        details[str(column)] = {
            "provider_nonconformity_quantiles": provider_quantiles,
            "cross_provider_calibration_eligible": all(
                value is not None for value in provider_quantiles.values()
            ),
            "selected_nonconformity_threshold": float(threshold[column]),
        }
    return threshold, details


def internal_uncertainty(
    labels: np.ndarray, probability: np.ndarray, conformal: np.ndarray
) -> dict[str, float]:
    included = (1.0 - probability) <= conformal
    positive = labels == 1
    return {
        "positive_label_coverage": float(included[positive].mean()),
        "mean_prediction_set_size": float(included.sum(axis=1).mean()),
        "multilabel_ECE_15_bin": multilabel_ece(labels, probability),
    }


def supervised_contrastive_loss(
    embedding: torch.Tensor, labels: torch.Tensor, temperature: float = 0.1
) -> torch.Tensor:
    normalized = F.normalize(embedding, dim=1)
    similarity = normalized @ normalized.T / temperature
    diagonal = torch.eye(len(labels), dtype=torch.bool, device=labels.device)
    intersection = labels @ labels.T
    union = labels.sum(1)[:, None] + labels.sum(1)[None, :] - intersection
    positive = (intersection > 0) & (intersection / union.clamp_min(1.0) >= 0.5) & ~diagonal
    masked_similarity = similarity.masked_fill(diagonal, -torch.inf)
    log_probability = masked_similarity - torch.logsumexp(masked_similarity, dim=1, keepdim=True)
    positive_count = positive.sum(1)
    eligible = positive_count > 0
    if not torch.any(eligible):
        return embedding.sum() * 0.0
    per_anchor = -(log_probability.masked_fill(~positive, 0.0).sum(1) / positive_count.clamp_min(1))
    return per_anchor[eligible].mean()


def augment(batch: torch.Tensor, rng: np.random.Generator) -> torch.Tensor:
    result = batch.clone()
    for index in range(len(result)):
        rotation = int(rng.integers(0, 4))
        if rotation:
            result[index] = torch.rot90(result[index], rotation, dims=(1, 2))
        if bool(rng.integers(0, 2)):
            result[index] = torch.flip(result[index], dims=(2,))
        if bool(rng.integers(0, 2)):
            result[index] = torch.flip(result[index], dims=(1,))
    return result


def unit_view(
    logits: np.ndarray, embeddings: np.ndarray, labels: np.ndarray,
    units: np.ndarray, mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return aggregate(logits[mask], embeddings[mask], labels[mask], units[mask])


def source_hashes(protocol: dict[str, Any], paths: dict[str, Path]) -> None:
    for key, path in paths.items():
        if sha256(path) != protocol["source_contract"][key]:
            raise ValueError(f"multiprovider source hash changed: {key}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--base-protocol", type=Path, required=True)
    parser.add_argument("--hpa-tensor", type=Path, required=True)
    parser.add_argument("--hpa-report", type=Path, required=True)
    parser.add_argument("--opencell-tensor", type=Path, required=True)
    parser.add_argument("--opencell-report", type=Path, required=True)
    parser.add_argument("--allencell-gate", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--smoke-only", action="store_true")
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    source_hashes(protocol, {
        "base_protocol_sha256": args.base_protocol,
        "hpa_tensor_sha256": args.hpa_tensor,
        "hpa_tensor_report_sha256": args.hpa_report,
        "opencell_tensor_sha256": args.opencell_tensor,
        "opencell_tensor_report_sha256": args.opencell_report,
        "allencell_metadata_gate_sha256": args.allencell_gate,
    })
    allen = json.loads(args.allencell_gate.read_text(encoding="utf-8"))
    if allen["Allen_pixels_accessed"] or allen["model_predictions_made"]:
        raise ValueError("Allen confirmation was accessed before model freeze")
    with np.load(args.hpa_tensor, allow_pickle=False) as source:
        hpa = {key: source[key].copy() for key in source.files}
    with np.load(args.opencell_tensor, allow_pickle=False) as source:
        opencell = {key: source[key].copy() for key in source.files}
    classes = hpa["class_names"].astype(str)
    if not np.array_equal(classes, opencell["class_names"].astype(str)):
        raise ValueError("provider class order differs")
    if list(hpa["images"].shape[1:]) != protocol["input_contract"]["shape"]:
        raise ValueError("HPA input shape differs from protocol")
    if list(opencell["images"].shape[1:]) != protocol["input_contract"]["shape"]:
        raise ValueError("OpenCell input shape differs from protocol")

    seed = int(protocol["optimization"]["seed"])
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(int(protocol["optimization"]["torch_threads"]))
    model = MultiproviderIFCNN(len(classes))
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(protocol["optimization"]["learning_rate"]),
        weight_decay=float(protocol["optimization"]["weight_decay"]),
    )

    hpa_split = hpa["split"].astype(str)
    oc_split = opencell["split"].astype(str)
    hpa_train = np.flatnonzero(hpa_split == "train")
    oc_train = np.flatnonzero(oc_split == "train")
    oc_target = opencell["target_name"].astype(str)
    oc_by_target: dict[str, list[int]] = defaultdict(list)
    for index in oc_train:
        oc_by_target[oc_target[index]].append(int(index))
    oc_train_targets = np.asarray(sorted(oc_by_target))

    hpa_unit_labels = hpa["labels"][hpa_train].astype(np.float32)
    oc_unit_labels = np.stack([
        opencell["labels"][oc_by_target[target][0]] for target in oc_train_targets
    ]).astype(np.float32)
    provider_positive_weights = []
    for unit_labels in (hpa_unit_labels, oc_unit_labels):
        count = len(unit_labels)
        positive = unit_labels.sum(axis=0)
        provider_positive_weights.append(
            torch.from_numpy(np.sqrt(count / np.maximum(positive, 1)).astype(np.float32))
        )

    def selection() -> tuple[dict[str, Any], tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]:
        h_logits, h_embedding = predict(model, hpa["images"][hpa_split == "selection"])
        h_labels = hpa["labels"][hpa_split == "selection"].astype(np.float32)
        o_mask = oc_split == "selection"
        o_logits_raw, o_embedding_raw = predict(model, opencell["images"][o_mask])
        o_logits, _, o_labels, _ = aggregate(
            o_logits_raw, o_embedding_raw, opencell["labels"][o_mask].astype(np.float32),
            oc_target[o_mask],
        )
        h_metrics = provider_selection_metrics(h_labels, sigmoid(h_logits))
        o_metrics = provider_selection_metrics(o_labels, sigmoid(o_logits))
        score = min(
            h_metrics["prevalence_normalized_AP"],
            o_metrics["prevalence_normalized_AP"],
        )
        return ({"HPA": h_metrics, "OpenCell": o_metrics, "score": score},
                (h_logits, h_labels), (o_logits, o_labels))

    batch_per_provider = int(protocol["optimization"]["provider_samples_per_batch"])
    steps = int(protocol["optimization"]["steps_per_epoch"])
    if steps * batch_per_provider > len(hpa_train):
        raise ValueError("frozen HPA epoch sampling exceeds train components")

    def train_epoch(epoch: int, step_limit: int | None = None) -> dict[str, float]:
        model.train()
        rng = np.random.default_rng(seed + epoch)
        hpa_order = rng.permutation(hpa_train)[: steps * batch_per_provider]
        target_order = []
        while len(target_order) < steps * batch_per_provider:
            target_order.extend(rng.permutation(oc_train_targets).tolist())
        target_order = target_order[: steps * batch_per_provider]
        occurrence: dict[str, int] = defaultdict(int)
        oc_order = []
        for target in target_order:
            choices = oc_by_target[str(target)]
            offset = int(hashlib.sha256(f"{seed}|{target}".encode()).hexdigest(), 16)
            position = (offset + epoch + occurrence[str(target)]) % len(choices)
            oc_order.append(choices[position])
            occurrence[str(target)] += 1
        totals = defaultdict(float)
        executed = min(steps, step_limit if step_limit is not None else steps)
        for step in range(executed):
            hs = hpa_order[step * batch_per_provider : (step + 1) * batch_per_provider]
            os = np.asarray(oc_order[step * batch_per_provider : (step + 1) * batch_per_provider])
            images = np.concatenate((hpa["images"][hs], opencell["images"][os]))
            labels = np.concatenate((hpa["labels"][hs], opencell["labels"][os])).astype(np.float32)
            provider = np.concatenate((np.zeros(len(hs), dtype=int), np.ones(len(os), dtype=int)))
            shuffle = rng.permutation(len(images))
            batch = torch.from_numpy(images[shuffle]).float().div_(255.0)
            target = torch.from_numpy(labels[shuffle])
            provider_tensor = torch.from_numpy(provider[shuffle])
            batch = augment(batch, rng)
            optimizer.zero_grad(set_to_none=True)
            logits, embedding = model(batch)
            provider_losses = []
            for provider_index in (0, 1):
                mask = provider_tensor == provider_index
                element = F.binary_cross_entropy_with_logits(
                    logits[mask], target[mask], reduction="none"
                )
                weights = torch.where(
                    target[mask] > 0,
                    provider_positive_weights[provider_index][None, :],
                    torch.ones_like(element),
                )
                provider_losses.append((element * weights).mean())
            bce = torch.stack(provider_losses).mean()
            contrastive = supervised_contrastive_loss(embedding, target)
            loss = bce + float(protocol["optimization"]["supervised_contrastive_weight"]) * contrastive
            if not torch.isfinite(loss):
                raise ValueError("nonfinite multiprovider training loss")
            loss.backward()
            gradient_norm = torch.nn.utils.clip_grad_norm_(
                model.parameters(), float(protocol["optimization"]["gradient_clip_global_norm"])
            )
            optimizer.step()
            totals["loss"] += float(loss.detach())
            totals["BCE"] += float(bce.detach())
            totals["contrastive"] += float(contrastive.detach())
            totals["gradient_norm"] += float(gradient_norm)
        return {key: value / executed for key, value in totals.items()}

    if args.smoke_only:
        values = train_epoch(0, step_limit=1)
        print(json.dumps({"status": "smoke_passed_no_artifact_written", **values}, sort_keys=True))
        return

    start_time = time.perf_counter()
    best_score = -float("inf")
    best_epoch = 0
    best_state: dict[str, torch.Tensor] | None = None
    best_selection: dict[str, Any] | None = None
    history = []
    patience = 0
    maximum_epochs = int(protocol["optimization"]["maximum_epochs"])
    minimum_improvement = float(protocol["optimization"]["minimum_improvement"])
    for epoch in range(maximum_epochs):
        train_values = train_epoch(epoch)
        selection_values, _, _ = selection()
        history.append({"epoch": epoch + 1, "training": train_values, "selection": selection_values})
        if selection_values["score"] > best_score + minimum_improvement:
            best_score = float(selection_values["score"])
            best_epoch = epoch + 1
            best_state = {
                name: value.detach().cpu().clone() for name, value in model.state_dict().items()
            }
            best_selection = selection_values
            patience = 0
        else:
            patience += 1
        print(json.dumps({
            "epoch": epoch + 1, "loss": train_values["loss"],
            "HPA_AP": selection_values["HPA"]["macro_average_precision"],
            "OpenCell_AP": selection_values["OpenCell"]["macro_average_precision"],
            "score": selection_values["score"], "best_epoch": best_epoch,
        }, sort_keys=True), flush=True)
        if patience >= int(protocol["optimization"]["early_stopping_patience_epochs"]):
            break
    if best_state is None or best_selection is None:
        raise AssertionError("training produced no checkpoint")
    model.load_state_dict(best_state)

    final_selection, h_selection, o_selection = selection()
    if abs(final_selection["score"] - best_score) > 1e-12:
        raise AssertionError("restored checkpoint differs from selected score")
    h_selection_probability = sigmoid(h_selection[0])
    o_selection_probability = sigmoid(o_selection[0])
    cutoffs = decision_thresholds([
        (h_selection[1], h_selection_probability),
        (o_selection[1], o_selection_probability),
    ])

    h_cal_mask = hpa_split == "calibration"
    h_cal_logits, h_cal_embedding = predict(model, hpa["images"][h_cal_mask])
    h_cal_labels = hpa["labels"][h_cal_mask].astype(np.float32)
    o_cal_mask = oc_split == "calibration"
    o_cal_raw_logits, o_cal_raw_embedding = predict(model, opencell["images"][o_cal_mask])
    o_cal_logits, o_cal_embedding, o_cal_labels, _ = aggregate(
        o_cal_raw_logits, o_cal_raw_embedding,
        opencell["labels"][o_cal_mask].astype(np.float32), oc_target[o_cal_mask],
    )
    temperature, calibration_nll = fit_temperature([
        (h_cal_logits, h_cal_labels), (o_cal_logits, o_cal_labels)
    ])
    h_cal_probability = sigmoid(h_cal_logits / temperature)
    o_cal_probability = sigmoid(o_cal_logits / temperature)
    conformal, conformal_details = positive_conformal(
        {"HPA": (h_cal_labels, h_cal_probability),
         "OpenCell": (o_cal_labels, o_cal_probability)},
        int(protocol["decision_and_uncertainty_contract"]["provider_class_calibration_minimum_positive_units"]),
        float(protocol["decision_and_uncertainty_contract"]["positive_conformal_alpha"]),
    )

    h_train_logits, h_train_embedding = predict(model, hpa["images"][hpa_split == "train"])
    h_train_labels = hpa["labels"][hpa_split == "train"].astype(np.float32)
    h_train_embedding = h_train_embedding / np.maximum(
        np.linalg.norm(h_train_embedding, axis=1, keepdims=True), 1e-12
    )
    o_train_mask = oc_split == "train"
    o_train_logits_raw, o_train_embedding_raw = predict(model, opencell["images"][o_train_mask])
    _, o_train_embedding, o_train_labels, _ = aggregate(
        o_train_logits_raw, o_train_embedding_raw,
        opencell["labels"][o_train_mask].astype(np.float32), oc_target[o_train_mask],
    )
    fit_embedding = np.concatenate((h_train_embedding, o_train_embedding))
    fit_labels = np.concatenate((h_train_labels, o_train_labels))
    centroids, variances = [], []
    for column in range(len(classes)):
        values = fit_embedding[fit_labels[:, column] == 1]
        centroids.append(values.mean(axis=0))
        variances.append(np.maximum(values.var(axis=0), 1e-6))
    centroids_array, variances_array = np.stack(centroids), np.stack(variances)

    def ood_score(embedding: np.ndarray) -> np.ndarray:
        distances = ((
            embedding[:, None, :] - centroids_array[None, :, :]
        ) ** 2 / variances_array[None, :, :]).mean(axis=2)
        return distances.min(axis=1)

    o_selection_mask = oc_split == "selection"
    o_sel_raw_logits, o_sel_raw_embedding = predict(model, opencell["images"][o_selection_mask])
    _, o_sel_embedding, _, _ = aggregate(
        o_sel_raw_logits, o_sel_raw_embedding,
        opencell["labels"][o_selection_mask].astype(np.float32), oc_target[o_selection_mask],
    )
    o_ood_mask = oc_split == "semantic_OOD"
    o_ood_raw_logits, o_ood_raw_embedding = predict(model, opencell["images"][o_ood_mask])
    _, o_ood_embedding, _, _ = aggregate(
        o_ood_raw_logits, o_ood_raw_embedding,
        opencell["labels"][o_ood_mask].astype(np.float32), oc_target[o_ood_mask],
    )
    in_score, out_score = ood_score(o_sel_embedding), ood_score(o_ood_embedding)
    development_ood_auroc = auroc(
        np.concatenate((np.zeros(len(in_score)), np.ones(len(out_score)))),
        np.concatenate((in_score, out_score)),
    )
    calibration_ood_scores = np.concatenate((
        ood_score(h_cal_embedding / np.maximum(np.linalg.norm(h_cal_embedding, axis=1, keepdims=True), 1e-12)),
        ood_score(o_cal_embedding),
    ))
    ood_threshold = float(np.quantile(calibration_ood_scores, 0.95, method="higher"))

    gates = protocol["selection_contract"]["development_open_gates"]
    gate_results = {
        "HPA_selection_macro_AP_minimum": best_selection["HPA"]["macro_average_precision"] >= gates["HPA_selection_macro_AP_minimum"],
        "OpenCell_selection_macro_AP_minimum": best_selection["OpenCell"]["macro_average_precision"] >= gates["OpenCell_selection_macro_AP_minimum"],
        "minimum_provider_prevalence_normalized_AP_minimum": best_score >= gates["minimum_provider_prevalence_normalized_AP_minimum"],
        "every_supported_provider_class_AP_above_prevalence": (
            best_selection["HPA"]["every_supported_class_AP_above_prevalence"]
            and best_selection["OpenCell"]["every_supported_class_AP_above_prevalence"]
        ),
        "OpenCell_semantic_OOD_AUROC_minimum": development_ood_auroc >= gates["OpenCell_semantic_OOD_AUROC_minimum"],
    }
    development_passed = all(gate_results.values())

    h_test_mask = hpa_split == "test"
    h_test_logits, _ = predict(model, hpa["images"][h_test_mask])
    h_test_probability = sigmoid(h_test_logits / temperature)
    h_test_labels = hpa["labels"][h_test_mask].astype(np.float32)
    h_test_metrics = provider_selection_metrics(h_test_labels, h_test_probability, minimum_positive=1)

    checkpoint = {
        "state_dict": best_state,
        "classes": classes.tolist(),
        "decision_thresholds": cutoffs.tolist(),
        "temperature": temperature,
        "positive_conformal_nonconformity_thresholds": conformal.tolist(),
        "OOD_class_centroids": centroids_array.tolist(),
        "OOD_class_diagonal_variances": variances_array.tolist(),
        "OOD_threshold": ood_threshold,
        "input_shape": protocol["input_contract"]["shape"],
        "seed": seed,
        "protocol_sha256": sha256(args.protocol),
        "HPA_tensor_sha256": sha256(args.hpa_tensor),
        "OpenCell_tensor_sha256": sha256(args.opencell_tensor),
        "selected_epoch": best_epoch,
        "development_open_condition_passed": development_passed,
        "Aleph_parameter_authority": "none",
    }
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(checkpoint, args.checkpoint)
    report = {
        "schema": "aleph.outer_library.multiprovider_if_model_report.v1",
        "protocol_sha256": sha256(args.protocol),
        "checkpoint_sha256": sha256(args.checkpoint),
        "source_hashes": protocol["source_contract"],
        "software": {"numpy": np.__version__, "torch": torch.__version__},
        "device": "CPU",
        "deterministic_algorithms": True,
        "architecture": {
            "family": protocol["architecture"]["family"],
            "parameter_count": sum(value.numel() for value in model.parameters()),
            "embedding_dimension": 64,
            "class_count": len(classes),
        },
        "classes": classes.tolist(),
        "training": {
            "epochs_executed": len(history),
            "selected_epoch": best_epoch,
            "selected_score": best_score,
            "history": history,
            "provider_positive_weights": {
                provider: dict(zip(classes.tolist(), weights.tolist()))
                for provider, weights in zip(("HPA", "OpenCell"), provider_positive_weights)
            },
            "elapsed_seconds": time.perf_counter() - start_time,
        },
        "selection": best_selection,
        "decision_thresholds": dict(zip(classes.tolist(), cutoffs.tolist())),
        "calibration": {
            "temperature": temperature,
            "equal_provider_mean_negative_log_likelihood": calibration_nll,
            "positive_conformal": {
                classes[int(index)]: value for index, value in conformal_details.items()
            },
            "internal_provider_checks": {
                "HPA": internal_uncertainty(h_cal_labels, h_cal_probability, conformal),
                "OpenCell": internal_uncertainty(o_cal_labels, o_cal_probability, conformal),
            },
            "external_uncertainty_authority": False,
        },
        "OOD": {
            "OpenCell_selection_in_domain_target_count": len(in_score),
            "OpenCell_semantic_OOD_target_count": len(out_score),
            "development_AUROC": development_ood_auroc,
            "calibration_score_95th_percentile_threshold": ood_threshold,
            "external_OOD_authority": False,
        },
        "HPA_already_opened_internal_test_diagnostic": h_test_metrics,
        "development_open_gates": gate_results,
        "all_development_open_gates_passed": development_passed,
        "next_external_confirmation_pixels_may_be_opened": development_passed,
        "Allen_confirmation_status": "structurally_refused_before_image_or_pixel_access",
        "Allen_pixels_seen": 0,
        "external_provider_validation_complete": False,
        "production_eligible": False,
        "uncertainty_authority": False,
        "observable_evidence_eligible": False,
        "Aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
        "may_remove_sweep_axis": False,
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "selected_epoch": best_epoch, "score": best_score,
        "development_OOD_AUROC": development_ood_auroc,
        "all_development_gates_passed": development_passed,
        "checkpoint_sha256": report["checkpoint_sha256"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
