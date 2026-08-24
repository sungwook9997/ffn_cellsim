#!/usr/bin/env python3
"""Train the frozen cross-study Light My Cells spatial prediction model."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class ConditionalUNet(nn.Module):
    def __init__(self, base: int, condition_channels: int = 7) -> None:
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Conv2d(1 + condition_channels, base, 3, padding=1), nn.ReLU(),
            nn.Conv2d(base, base, 3, padding=1), nn.ReLU(),
        )
        self.bottleneck = nn.Sequential(
            nn.Conv2d(base, 2 * base, 3, padding=1), nn.ReLU(),
            nn.Conv2d(2 * base, 2 * base, 3, padding=1), nn.ReLU(),
        )
        self.decoder = nn.Sequential(
            nn.Conv2d(3 * base, base, 3, padding=1), nn.ReLU(),
            nn.Conv2d(base, base, 3, padding=1), nn.ReLU(),
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
        encoded = self.encoder(torch.cat((image, condition), dim=1))
        bottleneck = self.bottleneck(F.avg_pool2d(encoded, 2))
        embedding = bottleneck.mean(dim=(-2, -1))
        up = F.interpolate(bottleneck, size=encoded.shape[-2:], mode="bilinear", align_corners=False)
        decoded = self.decoder(torch.cat((encoded, up), dim=1))
        return torch.sigmoid(self.output(decoded)), embedding


def batch_tensors(
    inputs: np.ndarray, targets: np.ndarray, target_index: np.ndarray,
    modality_index: np.ndarray, indices: np.ndarray,
    target_count: int, modality_count: int,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    image = torch.from_numpy(inputs[indices].astype(np.float32))[:, None]
    truth = torch.from_numpy(targets[indices].astype(np.float32))[:, None]
    target_hot = F.one_hot(
        torch.from_numpy(target_index[indices]), target_count
    ).to(torch.float32)
    modality_hot = F.one_hot(
        torch.from_numpy(modality_index[indices]), modality_count
    ).to(torch.float32)
    return image, truth, target_hot, modality_hot


def spatial_loss(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
    pixel = F.mse_loss(prediction, truth)
    dx = F.l1_loss(
        prediction[:, :, :, 1:] - prediction[:, :, :, :-1],
        truth[:, :, :, 1:] - truth[:, :, :, :-1],
    )
    dy = F.l1_loss(
        prediction[:, :, 1:, :] - prediction[:, :, :-1, :],
        truth[:, :, 1:, :] - truth[:, :, :-1, :],
    )
    return pixel + 0.05 * (dx + dy)


def predict(
    model: nn.Module, inputs: np.ndarray, targets: np.ndarray,
    target_index: np.ndarray, modality_index: np.ndarray,
    indices: np.ndarray, batch_size: int, target_count: int, modality_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    predictions, embeddings = [], []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(indices), batch_size):
            batch = indices[start:start + batch_size]
            image, _, target_hot, modality_hot = batch_tensors(
                inputs, targets, target_index, modality_index, batch,
                target_count, modality_count,
            )
            output, embedding = model(image, target_hot, modality_hot)
            predictions.append(output[:, 0].cpu().numpy())
            embeddings.append(embedding.cpu().numpy())
    return np.concatenate(predictions), np.concatenate(embeddings)


def pearson(prediction: np.ndarray, truth: np.ndarray) -> float:
    x = prediction.ravel().astype(np.float64)
    y = truth.ravel().astype(np.float64)
    x -= x.mean(); y -= y.mean()
    denominator = math.sqrt(float(x @ x) * float(y @ y))
    return float((x @ y) / denominator) if denominator > 0 else float("nan")


def global_ssim(prediction: np.ndarray, truth: np.ndarray) -> float:
    x = prediction.astype(np.float64)
    y = truth.astype(np.float64)
    mux, muy = float(x.mean()), float(y.mean())
    vx, vy = float(x.var()), float(y.var())
    covariance = float(((x - mux) * (y - muy)).mean())
    c1, c2 = 0.01 ** 2, 0.03 ** 2
    return float(
        ((2 * mux * muy + c1) * (2 * covariance + c2))
        / ((mux * mux + muy * muy + c1) * (vx + vy + c2))
    )


def group_metrics(
    predictions: np.ndarray, truth: np.ndarray, targets: np.ndarray,
    studies: np.ndarray,
) -> dict[str, Any]:
    rows = []
    for study in sorted(set(studies.tolist())):
        for target in sorted(set(targets[studies == study].tolist())):
            mask = (studies == study) & (targets == target)
            correlations = [
                pearson(predictions[index], truth[index])
                for index in np.flatnonzero(mask)
            ]
            correlations = [value for value in correlations if np.isfinite(value)]
            ssim = [
                global_ssim(predictions[index], truth[index])
                for index in np.flatnonzero(mask)
            ]
            rows.append({
                "study_id": study,
                "target_channel": target,
                "pair_count": int(mask.sum()),
                "MAE": float(np.abs(predictions[mask] - truth[mask]).mean()),
                "MSE": float(np.square(predictions[mask] - truth[mask]).mean()),
                "mean_image_Pearson": float(np.mean(correlations)) if correlations else None,
                "mean_global_SSIM": float(np.mean(ssim)),
            })
    by_target = {}
    for target in sorted(set(targets.tolist())):
        selected = [row for row in rows if row["target_channel"] == target]
        by_target[target] = {
            "study_count": len(selected),
            "pair_count": sum(row["pair_count"] for row in selected),
            "macro_study_MAE": float(np.mean([row["MAE"] for row in selected])),
            "macro_study_Pearson": float(np.mean([
                row["mean_image_Pearson"] for row in selected
                if row["mean_image_Pearson"] is not None
            ])),
            "macro_study_global_SSIM": float(np.mean([
                row["mean_global_SSIM"] for row in selected
            ])),
        }
    return {"per_study_target": rows, "per_target": by_target}


def required_target_macro_study_mae(
    predictions: np.ndarray, truth: np.ndarray, targets: np.ndarray,
    studies: np.ndarray, required_targets: tuple[str, ...],
) -> tuple[float, dict[str, float]]:
    target_values = {}
    for target in required_targets:
        study_values = []
        for study in sorted(set(studies[targets == target].tolist())):
            mask = (targets == target) & (studies == study)
            study_values.append(float(np.abs(predictions[mask] - truth[mask]).mean()))
        if not study_values:
            raise ValueError(f"selection lacks required target {target}")
        target_values[target] = float(np.mean(study_values))
    return float(np.mean(list(target_values.values()))), target_values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if protocol["training_tensor"]["sha256"] != sha256(args.tensor):
        raise ValueError("training protocol does not pin supplied tensor")
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
    target_order = protocol["conditioning"]["target_channel_one_hot_order"]
    modality_order = protocol["conditioning"]["input_modality_one_hot_order"]
    target_lookup = {name: index for index, name in enumerate(target_order)}
    modality_lookup = {name: index for index, name in enumerate(modality_order)}
    target_index = np.asarray([target_lookup[name] for name in target_names], dtype=np.int64)
    modality_index = np.asarray([modality_lookup[name] for name in modality_names], dtype=np.int64)
    fit_indices = np.flatnonzero(roles == "fit")
    selection_indices = np.flatnonzero(roles == "selection")
    calibration_indices = np.flatnonzero(roles == "calibration")
    if (len(fit_indices), len(selection_indices), len(calibration_indices)) != (3744, 61, 569):
        raise ValueError("frozen tensor role counts changed")

    fit_target_counts = Counter(target_names[fit_indices].tolist())
    weights = np.asarray([
        1.0 / math.sqrt(fit_target_counts[target_names[index]])
        for index in fit_indices
    ], dtype=np.float64)
    weights /= weights.sum()
    batch_size = int(protocol["optimization"]["batch_size"])
    candidate_reports = []
    global_best: tuple[tuple[float, int, int], dict[str, Any]] | None = None
    required_selection_targets = ("Nucleus", "Mitochondria")
    for base in protocol["architecture_selection"]["candidate_base_channels"]:
        rng = np.random.default_rng(seed)
        torch.manual_seed(seed + int(base))
        model = ConditionalUNet(int(base))
        optimizer = torch.optim.Adam(
            model.parameters(), lr=protocol["optimization"]["learning_rate"],
            weight_decay=protocol["optimization"]["weight_decay"],
        )
        history = []
        for epoch in range(1, int(protocol["architecture_selection"]["maximum_epochs"]) + 1):
            sampled = rng.choice(fit_indices, size=len(fit_indices), replace=True, p=weights)
            model.train(); losses = []
            for start in range(0, len(sampled), batch_size):
                batch = sampled[start:start + batch_size]
                image, truth, target_hot, modality_hot = batch_tensors(
                    inputs, truths, target_index, modality_index, batch,
                    len(target_order), len(modality_order),
                )
                optimizer.zero_grad(set_to_none=True)
                prediction, _ = model(image, target_hot, modality_hot)
                loss = spatial_loss(prediction, truth)
                loss.backward(); optimizer.step()
                losses.append(float(loss.detach()))
            selection_prediction, _ = predict(
                model, inputs, truths, target_index, modality_index,
                selection_indices, batch_size, len(target_order), len(modality_order),
            )
            selected_truth = truths[selection_indices].astype(np.float32)
            selected_names = target_names[selection_indices]
            score, target_mae = required_target_macro_study_mae(
                selection_prediction, selected_truth, selected_names,
                studies[selection_indices], required_selection_targets,
            )
            history.append({
                "epoch": epoch, "fit_sampled_loss": float(np.mean(losses)),
                "selection_required_target_macro_study_MAE": score,
                "selection_target_macro_study_MAE": target_mae,
            })
            key = (score, int(base), epoch)
            if global_best is None or key < global_best[0]:
                global_best = (key, {
                    "base_channels": int(base), "epoch": epoch,
                    "state_dict": copy.deepcopy(model.state_dict()),
                    "selection_prediction": selection_prediction.copy(),
                })
        candidate_reports.append({"base_channels": int(base), "history": history})
    assert global_best is not None
    best = global_best[1]
    model = ConditionalUNet(best["base_channels"])
    model.load_state_dict(best["state_dict"])

    selection_predictions, selection_embeddings = predict(
        model, inputs, truths, target_index, modality_index, selection_indices,
        batch_size, len(target_order), len(modality_order),
    )
    selection_truth = truths[selection_indices].astype(np.float32)
    selection_metrics = group_metrics(
        selection_predictions, selection_truth,
        target_names[selection_indices], studies[selection_indices],
    )
    fit_mean_maps = {}
    baseline_metrics = {}
    for target in target_order:
        mask = target_names[fit_indices] == target
        if not mask.any():
            continue
        mean_map = truths[fit_indices[mask]].astype(np.float32).mean(axis=0)
        fit_mean_maps[target] = torch.from_numpy(mean_map)
        selection_mask = target_names[selection_indices] == target
        if selection_mask.any():
            truth = selection_truth[selection_mask]
            per_study_model, per_study_baseline = {}, {}
            for study in sorted(set(studies[selection_indices][selection_mask].tolist())):
                study_mask = selection_mask & (studies[selection_indices] == study)
                study_truth = selection_truth[study_mask]
                per_study_model[study] = float(np.abs(
                    selection_predictions[study_mask] - study_truth
                ).mean())
                per_study_baseline[study] = float(np.abs(
                    mean_map[None] - study_truth
                ).mean())
            model_mae = float(np.mean(list(per_study_model.values())))
            baseline_mae = float(np.mean(list(per_study_baseline.values())))
            baseline_metrics[target] = {
                "model_macro_study_MAE": model_mae,
                "fit_mean_map_macro_study_MAE": baseline_mae,
                "per_study_model_MAE": per_study_model,
                "per_study_fit_mean_map_MAE": per_study_baseline,
                "relative_MAE_improvement": (
                    (baseline_mae - model_mae) / baseline_mae if baseline_mae > 0 else None
                ),
            }

    calibration_predictions, calibration_embeddings = predict(
        model, inputs, truths, target_index, modality_index, calibration_indices,
        batch_size, len(target_order), len(modality_order),
    )
    calibration_truth = truths[calibration_indices].astype(np.float32)
    calibration_targets = target_names[calibration_indices]
    calibration_studies = studies[calibration_indices]
    thresholds = {}
    calibration_report = {}
    for target in target_order:
        target_mask = calibration_targets == target
        target_studies = sorted(set(calibration_studies[target_mask].tolist()))
        if len(target_studies) < protocol["uncertainty_calibration"][
            "minimum_calibration_studies_per_target"
        ]:
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
        pred = calibration_predictions[target_mask]
        truth = calibration_truth[target_mask]
        lower = np.clip(pred - threshold, 0.0, 1.0)
        upper = np.clip(pred + threshold, 0.0, 1.0)
        calibration_report[target] = {
            "status": "calibrated_same_source_studies",
            "study_count": len(target_studies),
            "pair_count": int(target_mask.sum()),
            "study_90th_residual_quantiles": study_quantiles,
            "threshold": threshold,
            "pixel_coverage": float(((truth >= lower) & (truth <= upper)).mean()),
            "normalized_interval_width": float((upper - lower).mean()),
        }

    fit_predictions, fit_embeddings = predict(
        model, inputs, truths, target_index, modality_index, fit_indices,
        batch_size, len(target_order), len(modality_order),
    )
    embedding_mean = fit_embeddings.mean(axis=0).astype(np.float32)
    embedding_variance = np.maximum(fit_embeddings.var(axis=0), 1e-6).astype(np.float32)
    checkpoint = {
        "schema": "aleph.outer_library.lightmycells_checkpoint.v1",
        "protocol_sha256": sha256(args.protocol),
        "tensor_sha256": sha256(args.tensor),
        "base_channels": best["base_channels"],
        "selected_epoch": best["epoch"],
        "target_order": target_order,
        "modality_order": modality_order,
        "state_dict": model.state_dict(),
        "calibration_thresholds": thresholds,
        "fit_mean_maps": fit_mean_maps,
        "embedding_mean": torch.from_numpy(embedding_mean),
        "embedding_variance": torch.from_numpy(embedding_variance),
    }
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(checkpoint, args.checkpoint)
    selection_required = [
        selection_metrics["per_target"][target]
        for target in ("Nucleus", "Mitochondria")
    ]
    report = {
        "schema": "aleph.outer_library.lightmycells_model_training.v1",
        "protocol_sha256": sha256(args.protocol),
        "tensor_sha256": sha256(args.tensor),
        "checkpoint_sha256": sha256(args.checkpoint),
        "architecture": {
            "family": "conditional_small_UNet",
            "selected_base_channels": best["base_channels"],
            "selected_epoch": best["epoch"],
            "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
            "embedding_dimension": int(embedding_mean.size),
        },
        "candidate_selection": candidate_reports,
        "selection_metrics": selection_metrics,
        "selection_baseline_comparison": baseline_metrics,
        "selection_required_target_macro": {
            "MAE": float(np.mean([item["macro_study_MAE"] for item in selection_required])),
            "Pearson_correlation": float(np.mean([
                item["macro_study_Pearson"] for item in selection_required
            ])),
            "global_SSIM": float(np.mean([
                item["macro_study_global_SSIM"] for item in selection_required
            ])),
            "relative_MAE_improvement_over_fit_mean_target": float(np.mean([
                baseline_metrics[target]["relative_MAE_improvement"]
                for target in ("Nucleus", "Mitochondria")
            ])),
        },
        "uncertainty_calibration": calibration_report,
        "embedding_reference": {
            "fit_pair_count": len(fit_indices),
            "dimension": int(embedding_mean.size),
            "score": "diagonal_Mahalanobis_squared",
            "selection_score_mean": float(np.mean(np.square(
                selection_embeddings - embedding_mean
            ) / embedding_variance)),
            "calibration_score_mean": float(np.mean(np.square(
                calibration_embeddings - embedding_mean
            ) / embedding_variance)),
        },
        "fit_predictions_used_for_reporting_only": int(len(fit_predictions)),
        "independent_dataset_or_lab_confirmation_performed": False,
        "confirmation_pairs_seen": 0,
        "status": "checkpoint_and_same_source_calibration_frozen_confirmation_pending",
        "production_eligible": False,
        "observable_evidence_eligible": False,
        "uncertainty_authority": False,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
        "limitations": [
            "selection and calibration are study-disjoint but remain inside one multi-center source collection",
            "global SSIM is reported explicitly and is not a local-window SSIM implementation",
            "per-image normalization prohibits absolute fluorescence-intensity claims",
            "single-use confirmation studies remain unopened",
        ],
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "base_channels": best["base_channels"], "epoch": best["epoch"],
        "parameters": report["architecture"]["parameter_count"],
        "selection_required": report["selection_required_target_macro"],
        "calibrated_targets": sorted(thresholds),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
