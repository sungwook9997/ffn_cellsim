#!/usr/bin/env python3
"""Train and evaluate the frozen imaging-contract OOD v2 quality head."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize
import torch

from train_hpa_raw_if_cnn import sha256
from train_multiprovider_if_model import MultiproviderIFCNN, auroc


FAMILIES = (
    "zero_nucleus",
    "zero_target",
    "saturated_both",
    "channel_swap",
    "target_pixel_permutation",
    "target_tile_permutation",
    "cross_unit_target_mismatch",
)


def order_key(salt: str, *values: str) -> str:
    return hashlib.sha256("|".join((salt, *values)).encode()).hexdigest()


def seed64(salt: str, *values: str) -> int:
    return int(order_key(salt, *values)[:16], 16)


def role(salt: str, provider: str, unit_id: str) -> str:
    bucket = int(order_key(salt, "role", provider, unit_id), 16) % 20
    if bucket <= 2:
        return "OOD_calibration"
    if bucket <= 5:
        return "OOD_selection"
    return "OOD_fit"


def provider_units(
    hpa: dict[str, np.ndarray], opencell: dict[str, np.ndarray], salt: str
) -> dict[str, list[dict[str, Any]]]:
    output: dict[str, list[dict[str, Any]]] = defaultdict(list)
    hpa_split = hpa["split"].astype(str)
    for index in np.flatnonzero(hpa_split == "train"):
        unit_id = str(hpa["component_id"][index])
        output[role(salt, "HPA", unit_id)].append({
            "provider": "HPA", "unit_id": unit_id,
            "image": hpa["images"][index],
        })
    oc_split = opencell["split"].astype(str)
    targets = opencell["target_name"].astype(str)
    keys = opencell["projection_key"].astype(str)
    for target in sorted(set(targets[oc_split == "train"])):
        indices = np.flatnonzero((oc_split == "train") & (targets == target))
        selected = min(
            indices.tolist(), key=lambda index: order_key(salt, "projection", keys[index])
        )
        output[role(salt, "OpenCell", target)].append({
            "provider": "OpenCell", "unit_id": target,
            "image": opencell["images"][selected],
        })
    return output


def cap_roles(
    units: dict[str, list[dict[str, Any]]], protocol: dict[str, Any]
) -> dict[str, list[dict[str, Any]]]:
    salt = protocol["development_unit_contract"]["selection_salt"]
    limits = protocol["development_unit_contract"]["maximum_units_per_provider_role"]
    output: dict[str, list[dict[str, Any]]] = {}
    for role_name, values in units.items():
        selected = []
        for provider in ("HPA", "OpenCell"):
            provider_values = [row for row in values if row["provider"] == provider]
            provider_values.sort(
                key=lambda row: order_key(salt, "unit", provider, row["unit_id"])
            )
            selected.extend(provider_values[: int(limits[role_name])])
        output[role_name] = selected
    return output


def final_test_units(
    hpa: dict[str, np.ndarray], opencell: dict[str, np.ndarray], salt: str
) -> list[dict[str, Any]]:
    output = []
    hpa_split = hpa["split"].astype(str)
    for index in np.flatnonzero(hpa_split == "selection"):
        output.append({
            "provider": "HPA", "unit_id": str(hpa["component_id"][index]),
            "image": hpa["images"][index],
        })
    oc_split = opencell["split"].astype(str)
    targets = opencell["target_name"].astype(str)
    keys = opencell["projection_key"].astype(str)
    for target in sorted(set(targets[oc_split == "selection"])):
        indices = np.flatnonzero((oc_split == "selection") & (targets == target))
        selected = min(
            indices.tolist(), key=lambda index: order_key(salt, "projection", keys[index])
        )
        output.append({
            "provider": "OpenCell", "unit_id": target,
            "image": opencell["images"][selected],
        })
    return output


def corrupt(
    row: dict[str, Any], family: str, mismatch_target: np.ndarray, salt: str
) -> np.ndarray:
    image = row["image"].copy()
    if family == "zero_nucleus":
        image[0] = 0
    elif family == "zero_target":
        image[1] = 0
    elif family == "saturated_both":
        image[:] = 255
    elif family == "channel_swap":
        image = image[::-1].copy()
    elif family == "target_pixel_permutation":
        rng = np.random.default_rng(seed64(salt, "pixel", row["provider"], row["unit_id"]))
        image[1] = rng.permutation(image[1].reshape(-1)).reshape(128, 128)
    elif family == "target_tile_permutation":
        rng = np.random.default_rng(seed64(salt, "tile", row["provider"], row["unit_id"]))
        tiles = image[1].reshape(8, 16, 8, 16).transpose(0, 2, 1, 3).reshape(64, 16, 16)
        image[1] = tiles[rng.permutation(64)].reshape(8, 8, 16, 16).transpose(
            0, 2, 1, 3
        ).reshape(128, 128)
    elif family == "cross_unit_target_mismatch":
        image[1] = mismatch_target
    else:
        raise ValueError(f"unknown corruption family: {family}")
    return image


def variants(
    units: list[dict[str, Any]], salt: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    by_provider: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in units:
        by_provider[row["provider"]].append(row)
    mismatch: dict[tuple[str, str], np.ndarray] = {}
    for provider, rows in by_provider.items():
        ordered = sorted(
            rows, key=lambda row: order_key(salt, "mismatch", provider, row["unit_id"])
        )
        for index, row in enumerate(ordered):
            mismatch[(provider, row["unit_id"])] = ordered[(index + 1) % len(ordered)][
                "image"
            ][1]
    images, labels, providers, families = [], [], [], []
    for row in units:
        images.append(row["image"])
        labels.append(0)
        providers.append(row["provider"])
        families.append("clean")
        for family in FAMILIES:
            images.append(corrupt(
                row, family, mismatch[(row["provider"], row["unit_id"])], salt
            ))
            labels.append(1)
            providers.append(row["provider"])
            families.append(family)
    return (
        np.stack(images), np.asarray(labels, dtype=np.float64),
        np.asarray(providers), np.asarray(families),
    )


def qc_features(batch: np.ndarray) -> np.ndarray:
    value = batch.astype(np.float32) / 255.0
    rows = []
    for image in value:
        feature = []
        for channel in image:
            feature.extend((
                float(channel.mean()), float(channel.std()),
                float(np.percentile(channel, 1)), float(np.median(channel)),
                float(np.percentile(channel, 99)), float(np.mean(channel == 0)),
                float(np.mean(channel == 1)), float(np.mean(np.abs(np.diff(channel, axis=1)))),
                float(np.mean(np.abs(np.diff(channel, axis=0)))),
            ))
        a, b = image[0].reshape(-1), image[1].reshape(-1)
        a0, b0 = a - a.mean(), b - b.mean()
        denominator = float(np.linalg.norm(a0) * np.linalg.norm(b0))
        feature.append(float(np.dot(a0, b0) / denominator) if denominator > 0 else 0.0)
        rows.append(feature)
    return np.asarray(rows, dtype=np.float32)


def extract_features(
    model: MultiproviderIFCNN, images: np.ndarray, batch_size: int = 64
) -> np.ndarray:
    output = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(images), batch_size):
            raw = images[start : start + batch_size]
            tensor = torch.from_numpy(raw).float().div_(255.0)
            logits, embedding = model(tensor)
            normalized = torch.nn.functional.normalize(embedding, dim=1)
            output.append(np.concatenate((
                qc_features(raw), logits.cpu().numpy(), normalized.cpu().numpy()
            ), axis=1))
    result = np.concatenate(output).astype(np.float64)
    if result.shape[1] != 92 or not np.isfinite(result).all():
        raise ValueError(f"invalid v2 OOD feature tensor: {result.shape}")
    return result


def sample_weights(labels: np.ndarray, providers: np.ndarray, families: np.ndarray) -> np.ndarray:
    weights = np.zeros(len(labels), dtype=np.float64)
    for provider in ("HPA", "OpenCell"):
        provider_mask = providers == provider
        clean = provider_mask & (families == "clean")
        weights[clean] = 0.5 / clean.sum()
        for family in FAMILIES:
            mask = provider_mask & (families == family)
            weights[mask] = 0.5 / len(FAMILIES) / mask.sum()
    weights *= len(weights) / weights.sum()
    return weights


def sigmoid(value: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(value, -50.0, 50.0)))


def fit_logistic(
    features: np.ndarray, labels: np.ndarray, weights: np.ndarray, l2: float
) -> tuple[np.ndarray, dict[str, Any]]:
    design = np.column_stack((np.ones(len(features)), features))
    denominator = weights.sum()

    def objective(parameter: np.ndarray) -> tuple[float, np.ndarray]:
        linear = design @ parameter
        probability = sigmoid(linear)
        loss = np.sum(weights * (
            np.logaddexp(0.0, linear) - labels * linear
        )) / denominator + 0.5 * l2 * np.dot(parameter[1:], parameter[1:])
        gradient = design.T @ (weights * (probability - labels)) / denominator
        gradient[1:] += l2 * parameter[1:]
        return float(loss), gradient

    fit = minimize(
        objective, np.zeros(design.shape[1]), jac=True, method="L-BFGS-B",
        options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not fit.success or not np.isfinite(fit.x).all():
        raise ValueError(f"v2 OOD logistic fit failed: {fit.message}")
    return fit.x, {
        "success": bool(fit.success), "iterations": int(fit.nit),
        "objective": float(fit.fun), "gradient_max_abs": float(np.max(np.abs(fit.jac))),
    }


def standardize_fit(
    features: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.average(features, axis=0, weights=weights)
    variance = np.average((features - mean) ** 2, axis=0, weights=weights)
    scale = np.sqrt(np.maximum(variance, 1e-12))
    scale = np.maximum(scale, 1e-6)
    return (features - mean) / scale, mean, scale


def probabilities(features: np.ndarray, mean: np.ndarray, scale: np.ndarray,
                  parameter: np.ndarray) -> np.ndarray:
    standardized = (features - mean) / scale
    return sigmoid(parameter[0] + standardized @ parameter[1:])


def provider_family_metrics(
    score: np.ndarray, providers: np.ndarray, families: np.ndarray
) -> dict[str, Any]:
    values = {}
    aucs = []
    for provider in ("HPA", "OpenCell"):
        clean = (providers == provider) & (families == "clean")
        for family in FAMILIES:
            corrupt_mask = (providers == provider) & (families == family)
            mask = clean | corrupt_mask
            truth = corrupt_mask[mask].astype(float)
            value = auroc(truth, score[mask])
            values[f"{provider}:{family}"] = {
                "AUROC": value,
                "clean_units": int(clean.sum()),
                "corrupted_units": int(corrupt_mask.sum()),
            }
            aucs.append(value)
    return {
        "minimum_provider_by_family_AUROC": float(min(aucs)),
        "macro_provider_by_family_AUROC": float(np.mean(aucs)),
        "provider_by_family": values,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--model-report", type=Path, required=True)
    parser.add_argument("--v1-diagnostic", type=Path, required=True)
    parser.add_argument("--hpa-tensor", type=Path, required=True)
    parser.add_argument("--opencell-tensor", type=Path, required=True)
    parser.add_argument("--output-head", type=Path, required=True)
    parser.add_argument("--output-report", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    source_paths = {
        "localization_checkpoint_sha256": args.checkpoint,
        "localization_model_report_sha256": args.model_report,
        "v1_failure_diagnostic_sha256": args.v1_diagnostic,
        "hpa_tensor_sha256": args.hpa_tensor,
        "opencell_tensor_sha256": args.opencell_tensor,
    }
    for key, path in source_paths.items():
        if sha256(path) != protocol["source_contract"][key]:
            raise ValueError(f"v2 OOD source hash changed: {key}")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model = MultiproviderIFCNN(len(checkpoint["classes"]))
    model.load_state_dict(checkpoint["state_dict"])
    with np.load(args.hpa_tensor, allow_pickle=False) as source:
        hpa = {key: source[key].copy() for key in source.files}
    with np.load(args.opencell_tensor, allow_pickle=False) as source:
        opencell = {key: source[key].copy() for key in source.files}
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(8)
    salt = protocol["development_unit_contract"]["selection_salt"]
    roles = cap_roles(provider_units(hpa, opencell, salt), protocol)
    final_units = final_test_units(hpa, opencell, salt)
    role_overlaps: dict[str, set[str]] = defaultdict(set)
    for role_name, rows in roles.items():
        for row in rows:
            role_overlaps[f"{row['provider']}:{row['unit_id']}"] .add(role_name)
    if any(len(value) > 1 for value in role_overlaps.values()):
        raise ValueError("v2 OOD unit leaks across fit/selection/calibration")

    arrays = {}
    for role_name in ("OOD_fit", "OOD_selection"):
        images, labels, providers, families = variants(roles[role_name], salt)
        arrays[role_name] = {
            "features": extract_features(model, images), "labels": labels,
            "providers": providers, "families": families,
        }
    calibration_images = np.stack([row["image"] for row in roles["OOD_calibration"]])
    calibration_features = extract_features(model, calibration_images)
    calibration_providers = np.asarray([row["provider"] for row in roles["OOD_calibration"]])
    final_images, final_labels, final_providers, final_families = variants(final_units, salt)
    final_features = extract_features(model, final_images)

    fit_data = arrays["OOD_fit"]
    weights = sample_weights(
        fit_data["labels"], fit_data["providers"], fit_data["families"]
    )
    fit_standardized, feature_mean, feature_scale = standardize_fit(
        fit_data["features"], weights
    )
    candidates = []
    parameters = {}
    for l2 in protocol["quality_head_contract"]["L2_lambda_candidates"]:
        parameter, optimizer_report = fit_logistic(
            fit_standardized, fit_data["labels"], weights, float(l2)
        )
        selection = arrays["OOD_selection"]
        score = probabilities(
            selection["features"], feature_mean, feature_scale, parameter
        )
        metrics = provider_family_metrics(
            score, selection["providers"], selection["families"]
        )
        candidates.append({
            "L2_lambda": float(l2), "optimizer": optimizer_report,
            "selection": metrics,
        })
        parameters[float(l2)] = parameter
    chosen = max(
        candidates,
        key=lambda item: (item["selection"]["minimum_provider_by_family_AUROC"],
                          item["L2_lambda"]),
    )
    selected_l2 = chosen["L2_lambda"]
    parameter = parameters[selected_l2]
    calibration_score = probabilities(
        calibration_features, feature_mean, feature_scale, parameter
    )
    provider_thresholds = {
        provider: float(np.quantile(
            calibration_score[calibration_providers == provider], 0.95, method="higher"
        ))
        for provider in ("HPA", "OpenCell")
    }
    threshold = max(provider_thresholds.values())
    final_score = probabilities(final_features, feature_mean, feature_scale, parameter)
    final_metrics = provider_family_metrics(final_score, final_providers, final_families)
    false_refusal = {
        provider: float(np.mean(
            final_score[(final_providers == provider) & (final_families == "clean")]
            > threshold
        ))
        for provider in ("HPA", "OpenCell")
    }
    gate_contract = protocol["final_v2_development_test_gates"]
    gates = {
        "minimum_provider_by_corruption_family_AUROC": (
            final_metrics["minimum_provider_by_family_AUROC"]
            >= gate_contract["minimum_provider_by_corruption_family_AUROC"]
        ),
        "macro_provider_by_corruption_family_AUROC": (
            final_metrics["macro_provider_by_family_AUROC"]
            >= gate_contract["macro_provider_by_corruption_family_AUROC"]
        ),
        "maximum_clean_false_refusal_rate_each_provider": all(
            value <= gate_contract["maximum_clean_false_refusal_rate_each_provider"]
            for value in false_refusal.values()
        ),
        "all_seven_corruption_families_present_each_provider": (
            set(final_families) == {"clean", *FAMILIES}
            and all(
                sum((final_providers == provider) & (final_families == family)) > 0
                for provider in ("HPA", "OpenCell") for family in FAMILIES
            )
        ),
        "finite_features_probabilities_and_coefficients": bool(
            np.isfinite(final_features).all() and np.isfinite(final_score).all()
            and np.isfinite(parameter).all()
        ),
    }
    passed = all(gates.values())
    args.output_head.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_head,
        feature_mean=feature_mean,
        feature_scale=feature_scale,
        parameter=parameter,
        clean_refusal_threshold=np.asarray(threshold),
        selected_L2_lambda=np.asarray(selected_l2),
        feature_count=np.asarray(92),
        protocol_sha256=np.asarray(sha256(args.protocol)),
        localization_checkpoint_sha256=np.asarray(sha256(args.checkpoint)),
    )
    report = {
        "schema": "aleph.outer_library.multiprovider_if_ood_v2_report.v1",
        "protocol_sha256": sha256(args.protocol),
        "quality_head_sha256": sha256(args.output_head),
        "source_hashes": protocol["source_contract"],
        "v1_semantic_OOD_gate_remains_failed": True,
        "semantic_novelty_identifiable_from_image_alone": False,
        "role_unit_counts": {
            role_name: dict(sorted(Counter(row["provider"] for row in rows).items()))
            for role_name, rows in roles.items()
        },
        "cross_role_unit_overlap_count": sum(
            len(value) > 1 for value in role_overlaps.values()
        ),
        "final_development_test_unit_counts": dict(sorted(Counter(
            row["provider"] for row in final_units
        ).items())),
        "feature_count": 92,
        "corruption_families": list(FAMILIES),
        "candidate_selection": candidates,
        "selected_L2_lambda": selected_l2,
        "selected_optimizer": chosen["optimizer"],
        "calibration_clean_unit_counts": dict(sorted(Counter(calibration_providers).items())),
        "provider_clean_score_95th_percentiles": provider_thresholds,
        "clean_refusal_threshold": threshold,
        "final_v2_development_test": {
            **final_metrics,
            "clean_false_refusal_rate": false_refusal,
        },
        "development_gates": gates,
        "all_v2_development_gates_passed": passed,
        "replacement_external_provider_may_be_preregistered_and_opened_once": passed,
        "external_provider_pixels_accessed": 0,
        "external_provider_validation_complete": False,
        "production_eligible": False,
        "uncertainty_authority": False,
        "observable_evidence_eligible": False,
        "external_OOD_authority": False,
        "Aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
        "may_remove_sweep_axis": False,
    }
    args.output_report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "selected_L2": selected_l2,
        "minimum_AUROC": final_metrics["minimum_provider_by_family_AUROC"],
        "macro_AUROC": final_metrics["macro_provider_by_family_AUROC"],
        "false_refusal": false_refusal,
        "all_gates_passed": passed,
        "external_pixels_accessed": 0,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
