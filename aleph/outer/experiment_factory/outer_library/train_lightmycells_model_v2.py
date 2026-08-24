#!/usr/bin/env python3
"""Train the endpoint-aligned two-level Light My Cells spatial model."""

from __future__ import annotations

import argparse
from collections import Counter
import copy
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from train_lightmycells_model import (
    batch_tensors,
    group_metrics,
    predict,
    required_target_macro_study_mae,
    sha256,
)


class ResidualBlock(nn.Module):
    def __init__(self, channels: int) -> None:
        super().__init__()
        self.layers = nn.Sequential(
            nn.Conv2d(channels, channels, 3, padding=1), nn.ReLU(),
            nn.Conv2d(channels, channels, 3, padding=1),
        )

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        return F.relu(value + self.layers(value))


class ConditionalResidualUNet(nn.Module):
    def __init__(self, base: int = 16, condition_channels: int = 7) -> None:
        super().__init__()
        self.stem = nn.Sequential(
            nn.Conv2d(1 + condition_channels, base, 3, padding=1), nn.ReLU(),
            ResidualBlock(base),
        )
        self.down1 = nn.Sequential(
            nn.Conv2d(base, 2 * base, 3, stride=2, padding=1), nn.ReLU(),
            ResidualBlock(2 * base),
        )
        self.down2 = nn.Sequential(
            nn.Conv2d(2 * base, 4 * base, 3, stride=2, padding=1), nn.ReLU(),
            ResidualBlock(4 * base),
        )
        self.bottleneck = ResidualBlock(4 * base)
        self.up2_reduce = nn.Conv2d(4 * base, 2 * base, 1)
        self.up2 = nn.Sequential(
            nn.Conv2d(4 * base, 2 * base, 3, padding=1), nn.ReLU(),
            ResidualBlock(2 * base),
        )
        self.up1_reduce = nn.Conv2d(2 * base, base, 1)
        self.up1 = nn.Sequential(
            nn.Conv2d(2 * base, base, 3, padding=1), nn.ReLU(),
            ResidualBlock(base),
        )
        self.output = nn.Conv2d(base, 1, 1)

    def forward(
        self, image: torch.Tensor, target_onehot: torch.Tensor,
        modality_onehot: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        condition = torch.cat((target_onehot, modality_onehot), dim=1)
        condition = condition[:, :, None, None].expand(
            -1, -1, image.shape[-2], image.shape[-1]
        )
        level1 = self.stem(torch.cat((image, condition), dim=1))
        level2 = self.down1(level1)
        deep = self.bottleneck(self.down2(level2))
        embedding = deep.mean(dim=(-2, -1))
        up2 = F.interpolate(deep, size=level2.shape[-2:], mode="bilinear", align_corners=False)
        up2 = self.up2(torch.cat((self.up2_reduce(up2), level2), dim=1))
        up1 = F.interpolate(up2, size=level1.shape[-2:], mode="bilinear", align_corners=False)
        up1 = self.up1(torch.cat((self.up1_reduce(up1), level1), dim=1))
        return torch.sigmoid(self.output(up1)), embedding


def deterministic_flips(
    image: torch.Tensor, truth: torch.Tensor, indices: np.ndarray, epoch: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    horizontal = torch.from_numpy(((indices + epoch) % 2) == 0)
    vertical = torch.from_numpy((((indices // 2) + epoch) % 2) == 0)
    if horizontal.any():
        image[horizontal] = torch.flip(image[horizontal], dims=(-1,))
        truth[horizontal] = torch.flip(truth[horizontal], dims=(-1,))
    if vertical.any():
        image[vertical] = torch.flip(image[vertical], dims=(-2,))
        truth[vertical] = torch.flip(truth[vertical], dims=(-2,))
    return image, truth


def endpoint_loss(
    prediction: torch.Tensor, truth: torch.Tensor,
    weights: dict[str, float],
) -> tuple[torch.Tensor, dict[str, float]]:
    l1 = F.l1_loss(prediction, truth)
    pred_flat = prediction.flatten(1)
    truth_flat = truth.flatten(1)
    pred_center = pred_flat - pred_flat.mean(dim=1, keepdim=True)
    truth_center = truth_flat - truth_flat.mean(dim=1, keepdim=True)
    cosine = (
        (pred_center * truth_center).sum(dim=1)
        / (pred_center.square().sum(dim=1).sqrt()
           * truth_center.square().sum(dim=1).sqrt()).clamp_min(1e-6)
    )
    correlation_loss = (1.0 - cosine).mean()
    mu_x = pred_flat.mean(dim=1); mu_y = truth_flat.mean(dim=1)
    var_x = pred_flat.var(dim=1, unbiased=False); var_y = truth_flat.var(dim=1, unbiased=False)
    covariance = ((pred_flat - mu_x[:, None]) * (truth_flat - mu_y[:, None])).mean(dim=1)
    c1, c2 = 0.01 ** 2, 0.03 ** 2
    ssim = (
        (2 * mu_x * mu_y + c1) * (2 * covariance + c2)
        / ((mu_x.square() + mu_y.square() + c1) * (var_x + var_y + c2))
    )
    ssim_loss = (1.0 - ssim).mean()
    dx = F.l1_loss(
        prediction[:, :, :, 1:] - prediction[:, :, :, :-1],
        truth[:, :, :, 1:] - truth[:, :, :, :-1],
    )
    dy = F.l1_loss(
        prediction[:, :, 1:, :] - prediction[:, :, :-1, :],
        truth[:, :, 1:, :] - truth[:, :, :-1, :],
    )
    gradient = dx + dy
    total = (
        weights["pixel_L1"] * l1
        + weights["one_minus_centered_cosine"] * correlation_loss
        + weights["one_minus_global_SSIM"] * ssim_loss
        + weights["first_difference_L1"] * gradient
    )
    return total, {
        "L1": float(l1.detach()),
        "correlation_loss": float(correlation_loss.detach()),
        "SSIM_loss": float(ssim_loss.detach()),
        "gradient_L1": float(gradient.detach()),
    }


def baseline_comparison(
    predictions: np.ndarray, truth: np.ndarray, target_names: np.ndarray,
    studies: np.ndarray, fit_mean_maps: dict[str, np.ndarray],
) -> dict[str, Any]:
    result = {}
    for target, mean_map in fit_mean_maps.items():
        target_mask = target_names == target
        if not target_mask.any():
            continue
        model_values, baseline_values = {}, {}
        for study in sorted(set(studies[target_mask].tolist())):
            mask = target_mask & (studies == study)
            model_values[study] = float(np.abs(predictions[mask] - truth[mask]).mean())
            baseline_values[study] = float(np.abs(mean_map[None] - truth[mask]).mean())
        model_macro = float(np.mean(list(model_values.values())))
        baseline_macro = float(np.mean(list(baseline_values.values())))
        result[target] = {
            "model_macro_study_MAE": model_macro,
            "fit_mean_map_macro_study_MAE": baseline_macro,
            "relative_MAE_improvement": (
                (baseline_macro - model_macro) / baseline_macro
                if baseline_macro > 0 else None
            ),
            "per_study_model_MAE": model_values,
            "per_study_fit_mean_map_MAE": baseline_values,
        }
    return result


def selection_endpoints(
    metrics: dict[str, Any], baseline: dict[str, Any],
    required_targets: list[str], denominators: dict[str, float],
) -> dict[str, Any]:
    selected = [metrics["per_target"][target] for target in required_targets]
    values = {
        "Pearson": float(np.mean([row["macro_study_Pearson"] for row in selected])),
        "global_SSIM": float(np.mean([
            row["macro_study_global_SSIM"] for row in selected
        ])),
        "relative_MAE_improvement": float(np.mean([
            baseline[target]["relative_MAE_improvement"] for target in required_targets
        ])),
    }
    ratios = {name: values[name] / denominators[name] for name in values}
    return {
        "values": values,
        "endpoint_ratios": ratios,
        "minimum_endpoint_ratio": min(ratios.values()),
        "all_selection_endpoints_passed": all(value >= 1.0 for value in ratios.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if protocol["training_tensor"]["sha256"] != sha256(args.tensor):
        raise ValueError("v2 protocol does not pin supplied tensor")
    v1_report = args.protocol.parent / protocol["v1_diagnostic"]["report_path"]
    if sha256(v1_report) != protocol["v1_diagnostic"]["report_sha256"]:
        raise ValueError("v1 diagnostic changed after v2 protocol freeze")

    seed = int(protocol["optimization"]["seed"])
    np.random.seed(seed); torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(min(8, max(1, torch.get_num_threads())))
    with np.load(args.tensor, allow_pickle=False) as data:
        inputs = data["inputs"]
        truths = data["targets"]
        roles = data["role"]
        studies = data["study_id"]
        target_names = data["target_channel"]
        modality_names = data["input_modality"]
    target_order = ["Nucleus", "Mitochondria", "Tubulin", "Actin"]
    modality_order = ["BF", "PC", "DIC"]
    target_lookup = {name: index for index, name in enumerate(target_order)}
    modality_lookup = {name: index for index, name in enumerate(modality_order)}
    target_index = np.asarray([target_lookup[name] for name in target_names], dtype=np.int64)
    modality_index = np.asarray([modality_lookup[name] for name in modality_names], dtype=np.int64)
    fit_indices = np.flatnonzero(roles == "fit")
    selection_indices = np.flatnonzero(roles == "selection")
    calibration_indices = np.flatnonzero(roles == "calibration")
    fit_target_counts = Counter(target_names[fit_indices].tolist())
    sampling_weights = np.asarray([
        1.0 / math.sqrt(fit_target_counts[target_names[index]]) for index in fit_indices
    ], dtype=np.float64)
    sampling_weights /= sampling_weights.sum()
    fit_mean_maps = {
        target: truths[fit_indices[target_names[fit_indices] == target]].astype(
            np.float32
        ).mean(axis=0)
        for target in target_order if (target_names[fit_indices] == target).any()
    }
    required_targets = protocol["selection"]["required_targets"]
    denominators = protocol["selection"]["endpoint_denominators"]
    batch_size = int(protocol["optimization"]["batch_size"])
    model = ConditionalResidualUNet(int(protocol["model"]["base_channels"]))
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=protocol["optimization"]["learning_rate"],
        weight_decay=protocol["optimization"]["weight_decay"],
    )
    rng = np.random.default_rng(seed)
    history = []
    best: tuple[tuple[float, int], dict[str, Any]] | None = None
    for epoch in range(1, int(protocol["optimization"]["maximum_epochs"]) + 1):
        sampled = rng.choice(
            fit_indices, size=len(fit_indices), replace=True, p=sampling_weights
        )
        model.train(); component_rows = []
        for start in range(0, len(sampled), batch_size):
            batch = sampled[start:start + batch_size]
            image, truth, target_hot, modality_hot = batch_tensors(
                inputs, truths, target_index, modality_index, batch,
                len(target_order), len(modality_order),
            )
            image, truth = deterministic_flips(image, truth, batch, epoch)
            optimizer.zero_grad(set_to_none=True)
            prediction, _ = model(image, target_hot, modality_hot)
            loss, components = endpoint_loss(
                prediction, truth, protocol["optimization"]["loss_components"]
            )
            loss.backward(); optimizer.step()
            component_rows.append({"total": float(loss.detach()), **components})
        selection_predictions, _ = predict(
            model, inputs, truths, target_index, modality_index, selection_indices,
            batch_size, len(target_order), len(modality_order),
        )
        selection_truth = truths[selection_indices].astype(np.float32)
        metrics = group_metrics(
            selection_predictions, selection_truth,
            target_names[selection_indices], studies[selection_indices],
        )
        baseline = baseline_comparison(
            selection_predictions, selection_truth,
            target_names[selection_indices], studies[selection_indices], fit_mean_maps,
        )
        endpoints = selection_endpoints(metrics, baseline, required_targets, denominators)
        row = {
            "epoch": epoch,
            "fit_loss_components": {
                key: float(np.mean([item[key] for item in component_rows]))
                for key in component_rows[0]
            },
            "selection": endpoints,
        }
        history.append(row)
        key = (-endpoints["minimum_endpoint_ratio"], epoch)
        if best is None or key < best[0]:
            best = (key, {
                "epoch": epoch,
                "state_dict": copy.deepcopy(model.state_dict()),
                "selection_metrics": metrics,
                "selection_baseline": baseline,
                "selection_endpoints": endpoints,
            })
        print(json.dumps({
            "epoch": epoch,
            "minimum_endpoint_ratio": endpoints["minimum_endpoint_ratio"],
            **endpoints["values"],
        }, sort_keys=True), flush=True)
    assert best is not None
    selected = best[1]
    model.load_state_dict(selected["state_dict"])

    calibration_predictions, calibration_embeddings = predict(
        model, inputs, truths, target_index, modality_index, calibration_indices,
        batch_size, len(target_order), len(modality_order),
    )
    calibration_truth = truths[calibration_indices].astype(np.float32)
    calibration_targets = target_names[calibration_indices]
    calibration_studies = studies[calibration_indices]
    thresholds, calibration_report = {}, {}
    for target in target_order:
        target_mask = calibration_targets == target
        target_studies = sorted(set(calibration_studies[target_mask].tolist()))
        if len(target_studies) < 2:
            calibration_report[target] = {
                "status": "unsupported", "study_count": len(target_studies)
            }
            continue
        study_quantiles = {}
        for study in target_studies:
            mask = target_mask & (calibration_studies == study)
            residual = np.abs(calibration_predictions[mask] - calibration_truth[mask])
            study_quantiles[study] = float(np.quantile(residual, 0.90))
        threshold = max(study_quantiles.values())
        thresholds[target] = threshold
        pred, truth = calibration_predictions[target_mask], calibration_truth[target_mask]
        lower, upper = np.clip(pred - threshold, 0, 1), np.clip(pred + threshold, 0, 1)
        calibration_report[target] = {
            "status": "calibrated_same_source_studies",
            "study_count": len(target_studies), "pair_count": int(target_mask.sum()),
            "study_90th_residual_quantiles": study_quantiles,
            "threshold": threshold,
            "pixel_coverage": float(((truth >= lower) & (truth <= upper)).mean()),
            "normalized_interval_width": float((upper - lower).mean()),
        }
    _, fit_embeddings = predict(
        model, inputs, truths, target_index, modality_index, fit_indices,
        batch_size, len(target_order), len(modality_order),
    )
    embedding_mean = fit_embeddings.mean(axis=0).astype(np.float32)
    embedding_variance = np.maximum(fit_embeddings.var(axis=0), 1e-6).astype(np.float32)
    checkpoint = {
        "schema": "aleph.outer_library.lightmycells_checkpoint.v2",
        "protocol_sha256": sha256(args.protocol),
        "tensor_sha256": sha256(args.tensor),
        "base_channels": int(protocol["model"]["base_channels"]),
        "selected_epoch": selected["epoch"],
        "target_order": target_order, "modality_order": modality_order,
        "state_dict": model.state_dict(),
        "calibration_thresholds": thresholds,
        "fit_mean_maps": {
            key: torch.from_numpy(value) for key, value in fit_mean_maps.items()
        },
        "embedding_mean": torch.from_numpy(embedding_mean),
        "embedding_variance": torch.from_numpy(embedding_variance),
    }
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(checkpoint, args.checkpoint)
    readiness = selected["selection_endpoints"]["all_selection_endpoints_passed"]
    report = {
        "schema": "aleph.outer_library.lightmycells_model_training.v2",
        "protocol_sha256": sha256(args.protocol),
        "tensor_sha256": sha256(args.tensor),
        "checkpoint_sha256": sha256(args.checkpoint),
        "architecture": {
            "family": protocol["model"]["family"],
            "base_channels": protocol["model"]["base_channels"],
            "selected_epoch": selected["epoch"],
            "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
            "embedding_dimension": int(embedding_mean.size),
        },
        "history": history,
        "selection_metrics": selected["selection_metrics"],
        "selection_baseline_comparison": selected["selection_baseline"],
        "selection_endpoints": selected["selection_endpoints"],
        "uncertainty_calibration": calibration_report,
        "embedding_reference": {
            "fit_pair_count": len(fit_indices),
            "dimension": int(embedding_mean.size),
            "score": "diagonal_Mahalanobis_squared",
            "calibration_score_mean": float(np.mean(np.square(
                calibration_embeddings - embedding_mean
            ) / embedding_variance)),
        },
        "confirmation_open_condition_passed": readiness,
        "confirmation_pairs_seen": 0,
        "status": (
            "checkpoint_frozen_selection_ready_confirmation_pending"
            if readiness else "checkpoint_frozen_selection_endpoints_failed_confirmation_sealed"
        ),
        "production_eligible": False,
        "observable_evidence_eligible": False,
        "uncertainty_authority": False,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "selected_epoch": selected["epoch"],
        "parameters": report["architecture"]["parameter_count"],
        "selection": selected["selection_endpoints"],
        "confirmation_open_condition_passed": readiness,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
