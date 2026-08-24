#!/usr/bin/env python3
"""Evaluate a cross-assay IF translocation representation conservatively.

BBBC013 and BBBC014 differ in protein, perturbation, cell type, microscope,
resolution, and acquisition site.  The sealed BBBC014 dataset therefore tests
whether an image-derived translocation *observable* transfers.  Both datasets
share the same provider lineage, so this is never reported as an independent-
lab or Aleph-parameter validation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr


SCHEMA = "aleph.outer_library.if_translocation_evaluation.v1"
TRAIN_DATASET = "bbbc013-v1-fkhr-translocation"
EXTERNAL_DATASET = "bbbc014-v1-nfkb-translocation"
SCORE_FEATURE = "signal_nuclear_perinuclear_log_ratio"


def _as_text(values: np.ndarray) -> np.ndarray:
    return np.asarray(values).astype(str)


def _auc(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=np.int8)
    scores = np.asarray(scores, dtype=np.float64)
    positives = labels == 1
    negatives = labels == 0
    if positives.sum() == 0 or negatives.sum() == 0:
        raise ValueError("AUROC requires both classes")
    ranks = rankdata(scores, method="average")
    return float(
        (ranks[positives].sum() - positives.sum() * (positives.sum() + 1) / 2)
        / (positives.sum() * negatives.sum())
    )


def _dose_target(
    dataset: np.ndarray,
    treatment: np.ndarray,
    cell_type: np.ndarray,
    dose: np.ndarray,
) -> np.ndarray:
    """Map author doses to [0,1] independently within each assay context."""
    target = np.empty(len(dose), dtype=np.float64)
    contexts = np.asarray(
        [f"{d}|{t}|{c}" for d, t, c in zip(dataset, treatment, cell_type)],
        dtype=str,
    )
    for context in np.unique(contexts):
        mask = contexts == context
        values = dose[mask]
        if (values < 0).any():
            raise ValueError("dose values must be non-negative")
        transformed = np.log10(values + max(float(values.max()) * 1e-6, 1e-30))
        low, high = float(transformed.min()), float(transformed.max())
        if high <= low:
            raise ValueError(f"dose context has no variation: {context}")
        target[mask] = (transformed - low) / (high - low)
    return target


def _extreme_auc(target: np.ndarray, score: np.ndarray) -> dict[str, float | int]:
    selected = (target <= 0.25) | (target >= 0.75)
    labels = (target[selected] >= 0.75).astype(np.int8)
    return {
        "field_count": int(selected.sum()),
        "low_dose_fields": int((labels == 0).sum()),
        "high_dose_fields": int((labels == 1).sum()),
        "auroc": _auc(labels, score[selected]),
    }


def _context_metrics(
    mask: np.ndarray,
    target: np.ndarray,
    score: np.ndarray,
) -> dict[str, object]:
    if mask.sum() < 8:
        raise ValueError("IF evaluation context is underpopulated")
    correlation = float(spearmanr(target[mask], score[mask]).statistic)
    return {
        "field_count": int(mask.sum()),
        "spearman_dose_vs_translocation_score": correlation,
        "extreme_dose_discrimination": _extreme_auc(target[mask], score[mask]),
    }


def _row_bootstrap_ci(
    mask: np.ndarray,
    rows: np.ndarray,
    target: np.ndarray,
    score: np.ndarray,
    *,
    draws: int,
) -> list[float]:
    groups = np.unique(rows[mask])
    if len(groups) < 2:
        return [float("nan"), float("nan")]
    rng = np.random.default_rng(20260805)
    estimates: list[float] = []
    for _ in range(draws):
        sampled = rng.choice(groups, size=len(groups), replace=True)
        indexes = np.concatenate([np.flatnonzero(mask & (rows == group)) for group in sampled])
        value = float(spearmanr(target[indexes], score[indexes]).statistic)
        if np.isfinite(value):
            estimates.append(value)
    if not estimates:
        return [float("nan"), float("nan")]
    return [float(v) for v in np.quantile(estimates, (0.025, 0.975))]


def _fit_fixed_ridge(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = X_train.mean(axis=0)
    scale = X_train.std(axis=0)
    scale[scale < 1e-12] = 1.0
    train_z = (X_train - mean) / scale
    design = np.column_stack((np.ones(len(train_z)), train_z))
    penalty = np.eye(design.shape[1], dtype=np.float64) * 10.0
    penalty[0, 0] = 0.0
    weights = np.linalg.solve(design.T @ design + penalty, design.T @ y_train)
    prediction = np.column_stack((np.ones(len(X)), (X - mean) / scale)) @ weights
    return prediction, mean, scale


def _ood_scores(
    X_train: np.ndarray,
    X: np.ndarray,
) -> np.ndarray:
    mean = X_train.mean(axis=0)
    scale = X_train.std(axis=0)
    scale[scale < 1e-12] = 1.0
    train_z = (X_train - mean) / scale
    z = (X - mean) / scale
    covariance = np.cov(train_z, rowvar=False)
    shrinkage = 0.20
    covariance = (1.0 - shrinkage) * covariance + shrinkage * np.eye(covariance.shape[0])
    inverse = np.linalg.inv(covariance)
    return np.einsum("ij,jk,ik->i", z, inverse, z)


def evaluate(tensor_path: Path, *, bootstrap_draws: int = 2_000) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        required = {
            "X", "feature_names", "dataset_id", "lab_group", "sample_id",
            "plate_row", "column", "treatment", "cell_type", "dose",
            "dose_text", "source_sha256",
            "platemap_sha256",
        }
        if not required.issubset(data.files):
            raise ValueError(f"IF tensor misses fields: {sorted(required - set(data.files))}")
        X = np.asarray(data["X"], dtype=np.float64)
        feature_names = _as_text(data["feature_names"])
        dataset = _as_text(data["dataset_id"])
        labs = _as_text(data["lab_group"])
        sample_id = _as_text(data["sample_id"])
        rows = _as_text(data["plate_row"])
        treatment = _as_text(data["treatment"])
        cell_type = _as_text(data["cell_type"])
        dose = np.asarray(data["dose"], dtype=np.float64)
    if X.shape != (192, 12) or len(np.unique(sample_id)) != 192:
        raise ValueError(f"expected 192 intact paired-channel fields, found {X.shape}")
    if len(np.unique(dataset)) != 2 or len(np.unique(labs)) != 1:
        raise ValueError("BBBC dataset/lab grouping contract changed")
    if not np.isfinite(X).all() or SCORE_FEATURE not in feature_names:
        raise ValueError("invalid IF feature tensor")

    target = _dose_target(dataset, treatment, cell_type, dose)
    raw_score = X[:, int(np.flatnonzero(feature_names == SCORE_FEATURE)[0])]
    training_orientation_mask = (
        (dataset == TRAIN_DATASET) & (treatment == "Wortmannin") & np.isin(rows, ["A", "B", "C"])
    )
    orientation_correlation = float(
        spearmanr(target[training_orientation_mask], raw_score[training_orientation_mask]).statistic
    )
    orientation = 1.0 if orientation_correlation >= 0 else -1.0
    score = raw_score * orientation

    contexts: dict[str, object] = {}
    masks = {
        "bbbc013_wortmannin": (dataset == TRAIN_DATASET) & (treatment == "Wortmannin"),
        "bbbc013_ly294002_compound_transfer": (dataset == TRAIN_DATASET) & (treatment == "LY294002"),
        "bbbc014_mcf7_dataset_holdout": (dataset == EXTERNAL_DATASET) & np.char.startswith(cell_type, "MCF7"),
        "bbbc014_a549_dataset_holdout": (dataset == EXTERNAL_DATASET) & np.char.startswith(cell_type, "A549"),
    }
    for name, mask in masks.items():
        metrics = _context_metrics(mask, target, score)
        metrics["technical_row_bootstrap_spearman_ci95"] = _row_bootstrap_ci(
            mask, rows, target, score, draws=bootstrap_draws
        )
        metrics["bootstrap_unit"] = "author_replica_plate_row_not_biological_replicate"
        contexts[name] = metrics

    train_mask = training_orientation_mask
    calibration_mask = (
        (dataset == TRAIN_DATASET) & (treatment == "Wortmannin") & (rows == "D")
    )
    external_mask = dataset == EXTERNAL_DATASET
    prediction, _, _ = _fit_fixed_ridge(X[train_mask], target[train_mask], X)
    calibration_residuals = np.abs(target[calibration_mask] - prediction[calibration_mask])
    alpha = 0.10
    order = min(
        len(calibration_residuals) - 1,
        int(np.ceil((len(calibration_residuals) + 1) * (1.0 - alpha))) - 1,
    )
    radius = float(np.sort(calibration_residuals)[order])
    covered = np.abs(target[external_mask] - prediction[external_mask]) <= radius

    ood_score = _ood_scores(X[train_mask], X)
    reference_mask = train_mask | calibration_mask
    ood_labels = np.concatenate(
        (np.zeros(reference_mask.sum(), dtype=np.int8), np.ones(external_mask.sum(), dtype=np.int8))
    )
    ood_values = np.concatenate((ood_score[reference_mask], ood_score[external_mask]))
    ood_auroc = _auc(ood_labels, ood_values)

    mcf7 = contexts["bbbc014_mcf7_dataset_holdout"]
    a549 = contexts["bbbc014_a549_dataset_holdout"]
    dataset_holdout_passed = bool(
        mcf7["spearman_dose_vs_translocation_score"] >= 0.70
        and a549["spearman_dose_vs_translocation_score"] >= 0.70
        and mcf7["extreme_dose_discrimination"]["auroc"] >= 0.85
        and a549["extreme_dose_discrimination"]["auroc"] >= 0.85
    )
    uncertainty_passed = bool(
        covered.mean() >= 0.80 and 2.0 * radius <= 0.60 and dataset_holdout_passed
    )
    report: dict[str, object] = {
        "schema": SCHEMA,
        "task": "paired-channel_IF_nuclear_translocation_representation",
        "field_count": 192,
        "observation_unit": "one intact two-channel field per well",
        "feature_count": int(X.shape[1]),
        "score_feature": SCORE_FEATURE,
        "score_orientation_selected_on": (
            "BBBC013 Wortmannin rows A-C only; BBBC014 remained sealed"
        ),
        "training_orientation_spearman": orientation_correlation,
        "contexts": contexts,
        "dataset_holdout": {
            "train_dataset": TRAIN_DATASET,
            "sealed_external_dataset": EXTERNAL_DATASET,
            "different_protein": True,
            "different_perturbation": True,
            "different_cell_types": True,
            "different_acquisition_site_and_instrument": True,
            "shared_provider_lineage": True,
            "dataset_holdout_present": True,
            "independent_lab_holdout": False,
            "translocation_representation_passed": dataset_holdout_passed,
        },
        "uncertainty": {
            "method": "fixed_ridge_split_conformal",
            "train_group": "BBBC013 Wortmannin rows A-C",
            "calibration_group": "BBBC013 Wortmannin row D",
            "external_group": "all BBBC014 fields",
            "nominal_coverage": 0.90,
            "calibration_field_count": int(calibration_mask.sum()),
            "interval_radius_normalized_dose": radius,
            "external_coverage": float(covered.mean()),
            "external_mean_interval_width": 2.0 * radius,
            "uncertainty_holdout_passed": uncertainty_passed,
        },
        "ood": {
            "method": "training-standardized_shrunk_mahalanobis",
            "reference_fields": int(reference_mask.sum()),
            "external_fields": int(external_mask.sum()),
            "dataset_shift_auroc": ood_auroc,
            "dataset_ood_detected": bool(ood_auroc >= 0.80),
        },
        "limitations": [
            "Four author replica rows per context are plate-level technical replicas, not reported biological replicates.",
            "The two datasets share Ilya Ravkin provider lineage and cannot establish an independent-lab holdout.",
            "Dose-transfer validates an IF translocation observable, not a universal cell-state classifier.",
            "No segmentation ground truth is supplied; image-derived field features are not per-cell labels.",
        ],
        "representation_training_eligible": True,
        "observable_evidence_eligible": bool(dataset_holdout_passed and uncertainty_passed),
        "production_eligible": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "may_set_parameter_range": False,
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=2_000)
    args = parser.parse_args()
    report = evaluate(args.tensor, bootstrap_draws=args.bootstrap_draws)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
