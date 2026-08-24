#!/usr/bin/env python3
"""Evaluate the frozen Light My Cells checkpoint exactly once on confirmation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from train_lightmycells_model import group_metrics, predict, sha256
from train_lightmycells_model_v2 import ConditionalResidualUNet, baseline_comparison


def binary_auc(negative: np.ndarray, positive: np.ndarray) -> float:
    """Return pairwise AUROC with half credit for ties."""
    comparisons = positive[:, None] - negative[None, :]
    return float((np.count_nonzero(comparisons > 0) + 0.5 * np.count_nonzero(
        comparisons == 0
    )) / comparisons.size)


def study_scores(
    embeddings: np.ndarray, studies: np.ndarray,
    mean: np.ndarray, variance: np.ndarray,
) -> dict[str, float]:
    pair_scores = np.square(embeddings - mean[None]) / variance[None]
    pair_scores = pair_scores.sum(axis=1)
    return {
        study: float(pair_scores[studies == study].mean())
        for study in sorted(set(studies.tolist()))
    }


def coverage_report(
    predictions: np.ndarray, truths: np.ndarray, targets: np.ndarray,
    studies: np.ndarray, thresholds: dict[str, float], required: list[str],
) -> dict[str, Any]:
    rows = []
    for study in sorted(set(studies.tolist())):
        for target in sorted(set(targets[studies == study].tolist())):
            mask = (studies == study) & (targets == target)
            if target not in thresholds:
                rows.append({
                    "study_id": study, "target_channel": target,
                    "pair_count": int(mask.sum()), "status": "uncalibrated",
                })
                continue
            threshold = thresholds[target]
            prediction, truth = predictions[mask], truths[mask]
            lower = np.clip(prediction - threshold, 0.0, 1.0)
            upper = np.clip(prediction + threshold, 0.0, 1.0)
            rows.append({
                "study_id": study, "target_channel": target,
                "pair_count": int(mask.sum()), "status": "evaluated",
                "pixel_interval_coverage": float(
                    ((truth >= lower) & (truth <= upper)).mean()
                ),
                "normalized_interval_width": float((upper - lower).mean()),
            })
    per_target = {}
    for target in required:
        selected = [
            row for row in rows
            if row["target_channel"] == target and row["status"] == "evaluated"
        ]
        per_target[target] = {
            "study_count": len(selected),
            "pair_count": sum(row["pair_count"] for row in selected),
            "macro_study_pixel_interval_coverage": float(np.mean([
                row["pixel_interval_coverage"] for row in selected
            ])),
            "macro_study_normalized_interval_width": float(np.mean([
                row["normalized_interval_width"] for row in selected
            ])),
        }
    return {
        "per_study_target": rows,
        "per_target": per_target,
        "required_target_macro_study_pixel_interval_coverage": float(np.mean([
            per_target[target]["macro_study_pixel_interval_coverage"]
            for target in required
        ])),
    }


def foreground_weighted_mae(
    predictions: np.ndarray, truths: np.ndarray, targets: np.ndarray,
    studies: np.ndarray,
) -> dict[str, dict[str, float]]:
    """Report target-intensity-weighted error as a non-thresholded diagnostic."""
    result = {}
    for target in sorted(set(targets.tolist())):
        values = {}
        for study in sorted(set(studies[targets == target].tolist())):
            mask = (targets == target) & (studies == study)
            error, weight = np.abs(predictions[mask] - truths[mask]), 1.0 + truths[mask]
            values[study] = float((error * weight).sum() / weight.sum())
        result[target] = {
            "macro_study_foreground_weighted_MAE": float(np.mean(list(values.values()))),
            "per_study_foreground_weighted_MAE": values,
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-protocol", type=Path, required=True)
    parser.add_argument("--model-protocol", type=Path, required=True)
    parser.add_argument("--training-tensor", type=Path, required=True)
    parser.add_argument("--confirmation-tensor", type=Path, required=True)
    parser.add_argument("--confirmation-tensor-report", type=Path, required=True)
    parser.add_argument("--acquisition-report", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--model-report", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    base = json.loads(args.base_protocol.read_text(encoding="utf-8"))
    model_protocol = json.loads(args.model_protocol.read_text(encoding="utf-8"))
    model_report = json.loads(args.model_report.read_text(encoding="utf-8"))
    tensor_report = json.loads(args.confirmation_tensor_report.read_text(encoding="utf-8"))
    acquisition = json.loads(args.acquisition_report.read_text(encoding="utf-8"))
    if not model_report["confirmation_open_condition_passed"]:
        raise ValueError("frozen model did not earn confirmation access")
    if model_report["confirmation_pairs_seen"] != 0:
        raise ValueError("model report claims prior confirmation predictions")
    if model_report["protocol_sha256"] != sha256(args.model_protocol):
        raise ValueError("model report does not pin supplied model protocol")
    if model_report["checkpoint_sha256"] != sha256(args.checkpoint):
        raise ValueError("model report does not pin supplied checkpoint")
    if model_protocol["training_tensor"]["sha256"] != sha256(args.training_tensor):
        raise ValueError("model protocol does not pin supplied training tensor")
    if tensor_report["tensor_sha256"] != sha256(args.confirmation_tensor):
        raise ValueError("confirmation tensor report does not pin supplied tensor")
    if tensor_report["checkpoint_sha256"] != sha256(args.checkpoint):
        raise ValueError("confirmation tensor does not pin supplied checkpoint")
    if not tensor_report["confirmation_consumed"]:
        raise ValueError("confirmation tensor is not marked consumed")

    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    model = ConditionalResidualUNet(int(checkpoint["base_channels"]))
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    torch.set_num_threads(min(8, max(1, torch.get_num_threads())))
    target_order = checkpoint["target_order"]
    modality_order = checkpoint["modality_order"]
    target_lookup = {name: index for index, name in enumerate(target_order)}
    modality_lookup = {name: index for index, name in enumerate(modality_order)}

    with np.load(args.confirmation_tensor, allow_pickle=False) as data:
        inputs, truths = data["inputs"], data["targets"]
        studies, targets = data["study_id"], data["target_channel"]
        modalities, groups = data["input_modality"], data["experimenter_group"]
    target_indices = np.asarray([target_lookup[name] for name in targets], dtype=np.int64)
    modality_indices = np.asarray(
        [modality_lookup[name] for name in modalities], dtype=np.int64
    )
    indices = np.arange(len(inputs))
    predictions, confirmation_embeddings = predict(
        model, inputs, truths, target_indices, modality_indices, indices,
        int(model_protocol["optimization"]["batch_size"]),
        len(target_order), len(modality_order),
    )
    truth32 = truths.astype(np.float32)
    metrics = group_metrics(predictions, truth32, targets, studies)
    mean_maps = {
        key: value.numpy() for key, value in checkpoint["fit_mean_maps"].items()
    }
    baselines = baseline_comparison(
        predictions, truth32, targets, studies, mean_maps
    )
    coverage = coverage_report(
        predictions, truth32, targets, studies,
        checkpoint["calibration_thresholds"],
        base["evaluation_contract"]["required_targets"],
    )

    with np.load(args.training_tensor, allow_pickle=False) as data:
        train_inputs, train_truths = data["inputs"], data["targets"]
        train_roles, train_studies = data["role"], data["study_id"]
        train_targets, train_modalities = data["target_channel"], data["input_modality"]
    fit = np.flatnonzero(train_roles == "fit")
    train_target_indices = np.asarray(
        [target_lookup[name] for name in train_targets], dtype=np.int64
    )
    train_modality_indices = np.asarray(
        [modality_lookup[name] for name in train_modalities], dtype=np.int64
    )
    _, fit_embeddings = predict(
        model, train_inputs, train_truths, train_target_indices,
        train_modality_indices, fit,
        int(model_protocol["optimization"]["batch_size"]),
        len(target_order), len(modality_order),
    )
    embedding_mean = checkpoint["embedding_mean"].numpy()
    embedding_variance = checkpoint["embedding_variance"].numpy()
    fit_study_scores = study_scores(
        fit_embeddings, train_studies[fit], embedding_mean, embedding_variance
    )
    confirmation_study_scores = study_scores(
        confirmation_embeddings, studies, embedding_mean, embedding_variance
    )
    ood_auc = binary_auc(
        np.asarray(list(fit_study_scores.values())),
        np.asarray(list(confirmation_study_scores.values())),
    )

    required = base["evaluation_contract"]["required_targets"]
    spatial_required = [metrics["per_target"][target] for target in required]
    endpoint_values = {
        "macro_study_Pearson_correlation": float(np.mean([
            row["macro_study_Pearson"] for row in spatial_required
        ])),
        "macro_study_SSIM": float(np.mean([
            row["macro_study_global_SSIM"] for row in spatial_required
        ])),
        "relative_MAE_improvement_over_fit_mean_target": float(np.mean([
            baselines[target]["relative_MAE_improvement"] for target in required
        ])),
        "pixel_interval_coverage": coverage[
            "required_target_macro_study_pixel_interval_coverage"
        ],
        "embedding_OOD_AUROC": ood_auc,
    }
    frozen = base["evaluation_contract"]["frozen_endpoints"]
    endpoint_passes = {
        "macro_study_Pearson_correlation": endpoint_values[
            "macro_study_Pearson_correlation"
        ] >= frozen["macro_study_Pearson_correlation"]["value"],
        "macro_study_SSIM": endpoint_values["macro_study_SSIM"]
        >= frozen["macro_study_SSIM"]["value"],
        "relative_MAE_improvement_over_fit_mean_target": endpoint_values[
            "relative_MAE_improvement_over_fit_mean_target"
        ] >= frozen["relative_MAE_improvement_over_fit_mean_target"]["value"],
        "pixel_interval_coverage": (
            frozen["pixel_interval_coverage"]["value"][0]
            <= endpoint_values["pixel_interval_coverage"]
            <= frozen["pixel_interval_coverage"]["value"][1]
        ),
        "embedding_OOD_AUROC": endpoint_values["embedding_OOD_AUROC"]
        >= frozen["embedding_OOD_AUROC"]["value"],
    }
    confirmation_passed = all(endpoint_passes.values())
    novel_groups = set(acquisition["novel_confirmation_experimenter_groups"])
    novel_mask = np.asarray([group in novel_groups for group in groups])
    novel_metrics = (
        group_metrics(
            predictions[novel_mask], truth32[novel_mask],
            targets[novel_mask], studies[novel_mask],
        ) if novel_mask.any() else None
    )
    report = {
        "schema": "aleph.outer_library.lightmycells_confirmation_evaluation.v1",
        "base_protocol_sha256": sha256(args.base_protocol),
        "model_protocol_sha256": sha256(args.model_protocol),
        "model_report_sha256": sha256(args.model_report),
        "checkpoint_sha256": sha256(args.checkpoint),
        "training_tensor_sha256": sha256(args.training_tensor),
        "confirmation_tensor_sha256": sha256(args.confirmation_tensor),
        "confirmation_tensor_report_sha256": sha256(args.confirmation_tensor_report),
        "acquisition_report_sha256": sha256(args.acquisition_report),
        "confirmation_pair_count": len(inputs),
        "confirmation_study_count": len(set(studies.tolist())),
        "supported_targets": required,
        "descriptive_only_targets": sorted(set(targets.tolist()) - set(required)),
        "spatial_metrics": metrics,
        "foreground_weighted_MAE_diagnostic": foreground_weighted_mae(
            predictions, truth32, targets, studies
        ),
        "baseline_comparison": baselines,
        "uncertainty": coverage,
        "embedding_OOD": {
            "score": "fit_reference_diagonal_Mahalanobis_squared",
            "aggregation": "mean pair score within study before AUROC",
            "fit_study_count": len(fit_study_scores),
            "confirmation_study_count": len(confirmation_study_scores),
            "fit_study_scores": fit_study_scores,
            "confirmation_study_scores": confirmation_study_scores,
            "AUROC": ood_auc,
        },
        "novel_experimenter_group_partial_holdout": {
            "group_count": len(novel_groups),
            "groups": sorted(novel_groups),
            "pair_count": int(novel_mask.sum()),
            "metrics": novel_metrics,
        },
        "all_confirmation_labs_independent": acquisition[
            "all_confirmation_labs_independent"
        ],
        "frozen_endpoint_values": endpoint_values,
        "frozen_endpoint_passes": endpoint_passes,
        "all_applicable_confirmation_endpoints_passed": confirmation_passed,
        "status": (
            "validated_cross_study_spatial_observation"
            if confirmation_passed else "blocked_external_spatial_transfer"
        ),
        "confirmation_consumed": True,
        "may_be_reused_as_pristine": False,
        "production_eligible_for_this_spatial_task": confirmation_passed,
        "general_cell_state_production_eligible": False,
        "observable_evidence_eligible": confirmation_passed,
        "uncertainty_authority_for_this_spatial_task": confirmation_passed,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
        "purpose_note": (
            "The independent external task exists because Aleph sweeps are not yet "
            "computationally executable; this result does not replace future Aleph "
            "forward-sweep calibration or authorize a parameter value."
        ),
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "endpoint_values": endpoint_values,
        "endpoint_passes": endpoint_passes,
        "confirmation_passed": confirmation_passed,
        "novel_group_count": len(novel_groups),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
