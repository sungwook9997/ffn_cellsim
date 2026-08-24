#!/usr/bin/env python3
"""Train a study-adversarial, style-randomized label-free spatial model.

The consumed Light My Cells confirmation is intentionally inaccessible to this
program.  Only the original fit, selection, and calibration roles are read.
"""

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

from train_lightmycells_model import batch_tensors, group_metrics, predict, sha256
from train_lightmycells_model_v2 import (
    baseline_comparison,
    deterministic_flips,
    endpoint_loss,
    selection_endpoints,
)


class GradientReverse(torch.autograd.Function):
    @staticmethod
    def forward(ctx, value: torch.Tensor, strength: float) -> torch.Tensor:
        ctx.strength = float(strength)
        return value.view_as(value)

    @staticmethod
    def backward(ctx, gradient: torch.Tensor):
        return -ctx.strength * gradient, None


class NormalizedResidualBlock(nn.Module):
    def __init__(self, channels: int) -> None:
        super().__init__()
        self.conv1 = nn.Conv2d(channels, channels, 3, padding=1)
        self.norm1 = nn.InstanceNorm2d(channels, affine=True)
        self.conv2 = nn.Conv2d(channels, channels, 3, padding=1)
        self.norm2 = nn.InstanceNorm2d(channels, affine=True)

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        hidden = F.gelu(self.norm1(self.conv1(value)))
        return F.gelu(value + self.norm2(self.conv2(hidden)))


class DomainInvariantUNet(nn.Module):
    def __init__(
        self, base: int, condition_channels: int, study_count: int,
        reversal_strength: float,
    ) -> None:
        super().__init__()
        self.reversal_strength = float(reversal_strength)
        self.stem_conv = nn.Conv2d(1 + condition_channels, base, 3, padding=1)
        self.stem_norm = nn.InstanceNorm2d(base, affine=True)
        self.stem_block = NormalizedResidualBlock(base)
        self.down1_conv = nn.Conv2d(base, 2 * base, 3, stride=2, padding=1)
        self.down1_norm = nn.InstanceNorm2d(2 * base, affine=True)
        self.down1_block = NormalizedResidualBlock(2 * base)
        self.down2_conv = nn.Conv2d(2 * base, 4 * base, 3, stride=2, padding=1)
        self.down2_norm = nn.InstanceNorm2d(4 * base, affine=True)
        self.bottleneck = NormalizedResidualBlock(4 * base)
        self.up2_reduce = nn.Conv2d(4 * base, 2 * base, 1)
        self.up2_conv = nn.Conv2d(4 * base, 2 * base, 3, padding=1)
        self.up2_norm = nn.InstanceNorm2d(2 * base, affine=True)
        self.up2_block = NormalizedResidualBlock(2 * base)
        self.up1_reduce = nn.Conv2d(2 * base, base, 1)
        self.up1_conv = nn.Conv2d(2 * base, base, 3, padding=1)
        self.up1_norm = nn.InstanceNorm2d(base, affine=True)
        self.up1_block = NormalizedResidualBlock(base)
        self.output = nn.Conv2d(base, 1, 1)
        self.study_head = nn.Sequential(
            nn.Linear(4 * base + condition_channels, 2 * base),
            nn.GELU(),
            nn.Linear(2 * base, study_count),
        )

    def forward(
        self, image: torch.Tensor, target_onehot: torch.Tensor,
        modality_onehot: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        condition_vector = torch.cat((target_onehot, modality_onehot), dim=1)
        condition = condition_vector[:, :, None, None].expand(
            -1, -1, image.shape[-2], image.shape[-1]
        )
        level1 = self.stem_block(F.gelu(self.stem_norm(
            self.stem_conv(torch.cat((image, condition), dim=1))
        )))
        level2 = self.down1_block(F.gelu(self.down1_norm(self.down1_conv(level1))))
        deep = self.bottleneck(F.gelu(self.down2_norm(self.down2_conv(level2))))
        embedding = deep.mean(dim=(-2, -1))
        up2 = F.interpolate(deep, size=level2.shape[-2:], mode="bilinear", align_corners=False)
        up2 = self.up2_block(F.gelu(self.up2_norm(
            self.up2_conv(torch.cat((self.up2_reduce(up2), level2), dim=1))
        )))
        up1 = F.interpolate(up2, size=level1.shape[-2:], mode="bilinear", align_corners=False)
        up1 = self.up1_block(F.gelu(self.up1_norm(
            self.up1_conv(torch.cat((self.up1_reduce(up1), level1), dim=1))
        )))
        return torch.sigmoid(self.output(up1)), embedding

    def domain_logits(
        self, embedding: torch.Tensor, target_onehot: torch.Tensor,
        modality_onehot: torch.Tensor,
    ) -> torch.Tensor:
        reversed_embedding = GradientReverse.apply(
            embedding, self.reversal_strength
        )
        return self.study_head(torch.cat(
            (reversed_embedding, target_onehot, modality_onehot), dim=1
        ))


def deterministic_style_view(
    image: torch.Tensor, indices: np.ndarray, epoch: int, view: int,
    contract: dict[str, Any], seed: int,
) -> torch.Tensor:
    """Apply deterministic, input-only microscope-style perturbations."""
    index = torch.as_tensor(indices, dtype=image.dtype)

    def phase(multiplier: int, offset: int) -> torch.Tensor:
        raw = torch.remainder(
            index * multiplier + epoch * (offset + 37) + view * (offset + 101),
            10007,
        )
        return raw / 10006.0

    def scale(value: torch.Tensor, bounds: list[float]) -> torch.Tensor:
        return bounds[0] + value * (bounds[1] - bounds[0])

    gamma = scale(phase(7919, 17), contract["gamma_range"])
    gain = scale(phase(6151, 29), contract["gain_range"])
    bias = scale(phase(4621, 43), contract["bias_range"])
    noise_sd = scale(
        phase(3571, 61), contract["gaussian_noise_standard_deviation_range"]
    )
    value = image.clamp(0, 1).pow(gamma[:, None, None, None])
    value = value * gain[:, None, None, None] + bias[:, None, None, None]
    generator = torch.Generator(device=image.device)
    generator.manual_seed(seed + epoch * 1009 + view * 100003 + int(indices[0]))
    noise = torch.randn(
        value.shape, dtype=value.dtype, device=value.device, generator=generator
    )
    return (value + noise * noise_sd[:, None, None, None]).clamp(0, 1)


def binary_auroc(negative: np.ndarray, positive: np.ndarray) -> float:
    values = np.concatenate((negative, positive)).astype(np.float64)
    labels = np.concatenate((np.zeros(len(negative)), np.ones(len(positive))))
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=np.float64)
    cursor = 0
    while cursor < len(values):
        end = cursor + 1
        while end < len(values) and values[order[end]] == values[order[cursor]]:
            end += 1
        ranks[order[cursor:end]] = 0.5 * (cursor + 1 + end)
        cursor = end
    rank_sum = ranks[labels == 1].sum()
    return float((rank_sum - len(positive) * (len(positive) + 1) / 2) /
                 (len(positive) * len(negative)))


def leave_study_out_ood(
    fit_embeddings: np.ndarray, fit_studies: np.ndarray,
    calibration_embeddings: np.ndarray, calibration_studies: np.ndarray,
) -> dict[str, Any]:
    study_order = sorted(set(fit_studies.tolist()))
    centroids = {
        study: fit_embeddings[fit_studies == study].mean(axis=0)
        for study in study_order
    }
    residual = np.concatenate([
        fit_embeddings[fit_studies == study] - centroids[study][None]
        for study in study_order
    ])
    variance = residual.var(axis=0)
    pooled = float(np.mean(variance))
    shrinkage_variance = 0.9 * variance + 0.1 * max(pooled, 1e-6)
    shrinkage_variance = np.maximum(shrinkage_variance, 1e-6)

    def nearest(values: np.ndarray, excluded: list[str | None]) -> np.ndarray:
        result = []
        for row, excluded_study in zip(values, excluded, strict=True):
            candidates = [
                np.mean(np.square(row - centroid) / shrinkage_variance)
                for study, centroid in centroids.items() if study != excluded_study
            ]
            result.append(min(candidates))
        return np.asarray(result, dtype=np.float64)

    fit_scores = nearest(fit_embeddings, fit_studies.tolist())
    calibration_scores = nearest(
        calibration_embeddings, [None] * len(calibration_embeddings)
    )
    fit_study_scores = {
        study: float(fit_scores[fit_studies == study].mean()) for study in study_order
    }
    calibration_study_scores = {
        study: float(calibration_scores[calibration_studies == study].mean())
        for study in sorted(set(calibration_studies.tolist()))
    }
    return {
        "score": "shrinkage_whitened_nearest_fit_study_centroid_distance",
        "fit_scoring": "own_study_centroid_excluded",
        "fit_study_scores": fit_study_scores,
        "calibration_study_scores": calibration_study_scores,
        "study_macro_AUROC": binary_auroc(
            np.asarray(list(fit_study_scores.values())),
            np.asarray(list(calibration_study_scores.values())),
        ),
        "authority": "post_selection_internal_shift_diagnostic_only",
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
        raise ValueError("v3 protocol does not pin supplied tensor")
    for name, expected in (
        ("lightmycells_model_training_protocol_v2.json", protocol["v2_artifacts"]["protocol_sha256"]),
        ("results/lightmycells_model_training_report_v2.json", protocol["v2_artifacts"]["report_sha256"]),
        ("results/lightmycells_model_checkpoint_v2.pt", protocol["v2_artifacts"]["checkpoint_sha256"]),
        (protocol["consumed_confirmation"]["report_path"], protocol["consumed_confirmation"]["report_sha256"]),
    ):
        if sha256(args.protocol.parent / name) != expected:
            raise ValueError(f"pinned prior artifact changed: {name}")

    seed = int(protocol["optimization"]["seed"])
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(min(8, max(1, torch.get_num_threads())))
    device = torch.device("cpu")
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
    if np.any(roles == "confirmation"):
        raise ValueError("v3 tensor must contain zero confirmation pairs")
    fit_study_order = sorted(set(studies[fit_indices].tolist()))
    fit_study_lookup = {name: index for index, name in enumerate(fit_study_order)}
    study_index = np.asarray([
        fit_study_lookup.get(name, -1) for name in studies
    ], dtype=np.int64)

    fit_target_counts = Counter(target_names[fit_indices].tolist())
    sampling_weights = np.asarray([
        1.0 / math.sqrt(fit_target_counts[target_names[index]]) for index in fit_indices
    ], dtype=np.float64)
    sampling_weights /= sampling_weights.sum()
    fit_mean_maps = {
        target: truths[fit_indices[target_names[fit_indices] == target]].astype(np.float32).mean(axis=0)
        for target in target_order if (target_names[fit_indices] == target).any()
    }
    generalization = protocol["domain_generalization"]
    model = DomainInvariantUNet(
        int(protocol["model"]["base_channels"]),
        len(target_order) + len(modality_order),
        len(fit_study_order),
        float(generalization["gradient_reversal_strength"]),
    ).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=protocol["optimization"]["learning_rate"],
        weight_decay=protocol["optimization"]["weight_decay"],
    )
    rng = np.random.default_rng(seed)
    batch_size = int(protocol["optimization"]["batch_size"])
    history = []
    best: tuple[tuple[float, int], dict[str, Any]] | None = None
    for epoch in range(1, int(protocol["optimization"]["maximum_epochs"]) + 1):
        sampled = rng.choice(
            fit_indices, size=len(fit_indices), replace=True, p=sampling_weights
        )
        model.train()
        component_rows = []
        for start in range(0, len(sampled), batch_size):
            batch = sampled[start:start + batch_size]
            image, truth, target_hot, modality_hot = batch_tensors(
                inputs, truths, target_index, modality_index, batch,
                len(target_order), len(modality_order),
            )
            image, truth = deterministic_flips(image, truth, batch, epoch)
            view1 = deterministic_style_view(
                image, batch, epoch, 1, generalization["input_only_augmentations"], seed
            )
            view2 = deterministic_style_view(
                image, batch, epoch, 2, generalization["input_only_augmentations"], seed
            )
            domain_y = torch.from_numpy(study_index[batch]).long()
            optimizer.zero_grad(set_to_none=True)
            prediction1, embedding1 = model(view1, target_hot, modality_hot)
            prediction2, embedding2 = model(view2, target_hot, modality_hot)
            supervised1, components1 = endpoint_loss(
                prediction1, truth, protocol["optimization"]["supervised_loss_components"]
            )
            supervised2, _ = endpoint_loss(
                prediction2, truth, protocol["optimization"]["supervised_loss_components"]
            )
            supervised = 0.5 * (supervised1 + supervised2)
            consistency = F.l1_loss(prediction1, prediction2)
            domain_logits1 = model.domain_logits(embedding1, target_hot, modality_hot)
            domain_logits2 = model.domain_logits(embedding2, target_hot, modality_hot)
            domain_loss = 0.5 * (
                F.cross_entropy(domain_logits1, domain_y)
                + F.cross_entropy(domain_logits2, domain_y)
            )
            total = (
                supervised
                + float(generalization["prediction_consistency_loss_weight"]) * consistency
                + float(generalization["study_adversary_loss_weight"]) * domain_loss
            )
            total.backward()
            optimizer.step()
            component_rows.append({
                "total": float(total.detach()),
                "supervised": float(supervised.detach()),
                "style_consistency_L1": float(consistency.detach()),
                "study_adversary_cross_entropy": float(domain_loss.detach()),
                **components1,
            })

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
        endpoints = selection_endpoints(
            metrics, baseline, protocol["selection"]["required_targets"],
            protocol["selection"]["endpoint_denominators"],
        )
        means = {
            key: float(np.mean([row[key] for row in component_rows]))
            for key in component_rows[0]
        }
        history.append({"epoch": epoch, "fit": means, "selection": endpoints})
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
            "style_consistency_L1": means["style_consistency_L1"],
            "study_adversary_cross_entropy": means["study_adversary_cross_entropy"],
        }, sort_keys=True), flush=True)

    assert best is not None
    selected = best[1]
    model.load_state_dict(selected["state_dict"])
    fit_predictions, fit_embeddings = predict(
        model, inputs, truths, target_index, modality_index, fit_indices,
        batch_size, len(target_order), len(modality_order),
    )
    calibration_predictions, calibration_embeddings = predict(
        model, inputs, truths, target_index, modality_index, calibration_indices,
        batch_size, len(target_order), len(modality_order),
    )
    calibration_truth = truths[calibration_indices].astype(np.float32)
    calibration_targets = target_names[calibration_indices]
    calibration_studies = studies[calibration_indices]
    thresholds, uncertainty_report = {}, {}
    for target in target_order:
        target_mask = calibration_targets == target
        target_studies = sorted(set(calibration_studies[target_mask].tolist()))
        if len(target_studies) < 2:
            uncertainty_report[target] = {
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
        pred = calibration_predictions[target_mask]
        truth = calibration_truth[target_mask]
        lower, upper = np.clip(pred - threshold, 0, 1), np.clip(pred + threshold, 0, 1)
        uncertainty_report[target] = {
            "status": "calibrated_internal_studies",
            "study_count": len(target_studies),
            "pair_count": int(target_mask.sum()),
            "study_90th_residual_quantiles": study_quantiles,
            "threshold": threshold,
            "pixel_coverage": float(((truth >= lower) & (truth <= upper)).mean()),
            "normalized_interval_width": float((upper - lower).mean()),
        }
    ood = leave_study_out_ood(
        fit_embeddings, studies[fit_indices],
        calibration_embeddings, calibration_studies,
    )

    checkpoint = {
        "schema": "aleph.outer_library.lightmycells_checkpoint.v3",
        "protocol_sha256": sha256(args.protocol),
        "tensor_sha256": sha256(args.tensor),
        "selected_epoch": selected["epoch"],
        "target_order": target_order,
        "modality_order": modality_order,
        "fit_study_order": fit_study_order,
        "model": {
            "base_channels": int(protocol["model"]["base_channels"]),
            "condition_channels": len(target_order) + len(modality_order),
            "study_count": len(fit_study_order),
            "gradient_reversal_strength": float(generalization["gradient_reversal_strength"]),
        },
        "state_dict": model.state_dict(),
        "calibration_thresholds": thresholds,
        "fit_mean_maps": {key: torch.from_numpy(value) for key, value in fit_mean_maps.items()},
    }
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(checkpoint, args.checkpoint)
    report = {
        "schema": "aleph.outer_library.lightmycells_model_training.v3",
        "protocol_sha256": sha256(args.protocol),
        "tensor_sha256": sha256(args.tensor),
        "checkpoint_sha256": sha256(args.checkpoint),
        "device": str(device),
        "architecture": {
            "family": protocol["model"]["family"],
            "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
            "selected_epoch": selected["epoch"],
            "fit_study_adversary_classes": len(fit_study_order),
        },
        "consumed_confirmation_access": {
            "pairs_read": 0,
            "images_read": 0,
            "metrics_used_for_selection": 0,
            "may_be_re_evaluated": False,
        },
        "history": history,
        "selection_metrics": selected["selection_metrics"],
        "selection_baseline_comparison": selected["selection_baseline"],
        "selection_endpoints": selected["selection_endpoints"],
        "uncertainty_calibration": uncertainty_report,
        "internal_unseen_study_OOD": ood,
        "status": "training_only_new_external_confirmation_required",
        "production_eligible": False,
        "observable_evidence_eligible": False,
        "uncertainty_authority": False,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "selected_epoch": selected["epoch"],
        "selection": selected["selection_endpoints"],
        "internal_unseen_study_OOD_AUROC": ood["study_macro_AUROC"],
        "checkpoint_sha256": report["checkpoint_sha256"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
