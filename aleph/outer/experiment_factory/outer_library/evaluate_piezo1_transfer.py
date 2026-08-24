#!/usr/bin/env python3
"""Independent-lab transfer test for a binary Piezo1-inhibition signature."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from train_calcium_baseline import (
    binary_auroc,
    curve_features,
    fit_ridge,
    logits,
    probabilities,
    macro_f1,
)


ACTIVE_NANOSWITCH_STRATA = ("LooPINS", "CaPINS")


def _curves(
    path: Path,
    source: str,
    *,
    active_only: bool = False,
) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if source == "hela" and row["observable"] != "intracellular_calcium_F_over_F0":
                continue
            grouped[row["sample_id"]].append(row)
    features, labels, conditions, sample_ids = [], [], [], []
    for sample_id, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: row["coordinate_value"])
        condition = rows[0]["condition"]
        if source == "hela":
            if condition == "TRPV4_inhibited_HC-067047":
                continue
            label = int(condition == "Piezo1_inhibited_GsMTx4")
        else:
            if active_only and not condition.startswith(ACTIVE_NANOSWITCH_STRATA):
                continue
            if "+GsMTx4" in condition:
                label = 1
            elif "+Mag" in condition:
                label = 0
            else:
                continue
        values = np.asarray([row["value"] for row in rows], dtype=float)
        values = values - values[0]
        features.append(curve_features(values[None, :])[0])
        labels.append(label)
        conditions.append(condition)
        sample_ids.append(sample_id)
    return np.stack(features), np.asarray(labels, dtype=int), conditions, sample_ids


def _invariant_features(path: Path, source: str, *, active_only: bool) -> tuple[
    np.ndarray, np.ndarray, list[str], list[str]
]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if source == "hela" and row["observable"] != "intracellular_calcium_F_over_F0":
                continue
            grouped[row["sample_id"]].append(row)
    features, labels, conditions, sample_ids = [], [], [], []
    target_axis = np.linspace(0.0, 1.0, 32)
    for sample_id, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: row["coordinate_value"])
        condition = rows[0]["condition"]
        if source == "hela":
            if condition == "TRPV4_inhibited_HC-067047":
                continue
            label = int(condition == "Piezo1_inhibited_GsMTx4")
        else:
            if active_only and not condition.startswith(ACTIVE_NANOSWITCH_STRATA):
                continue
            if "+GsMTx4" in condition:
                label = 1
            elif "+Mag" in condition:
                label = 0
            else:
                continue
        values = np.asarray([row["value"] for row in rows], dtype=np.float64)
        values -= values[0]
        source_axis = np.linspace(0.0, 1.0, len(values))
        sampled = np.interp(target_axis, source_axis, values)
        amplitude = max(float(np.max(np.abs(values))), float(np.ptp(values)), 1e-8)
        shape = sampled / amplitude
        summaries = np.asarray([
            values.mean() / amplitude,
            values.std() / amplitude,
            values.min() / amplitude,
            values.max() / amplitude,
            values[-1] / amplitude,
            np.ptp(values) / amplitude,
            np.argmin(values) / max(len(values) - 1, 1),
            np.argmax(values) / max(len(values) - 1, 1),
            np.trapezoid(values, source_axis) / amplitude,
        ])
        features.append(np.concatenate((shape, summaries)))
        labels.append(label)
        conditions.append(condition)
        sample_ids.append(sample_id)
    return np.stack(features), np.asarray(labels, dtype=int), conditions, sample_ids


def _metrics(probs: np.ndarray, labels: np.ndarray) -> dict:
    pred = probs.argmax(axis=1)
    result = {
        "sample_count": len(labels),
        "accuracy": float(np.mean(pred == labels)),
        "macro_f1": macro_f1(pred, labels, 2),
        "balanced_accuracy": float(np.mean([
            np.mean(pred[labels == label] == label) for label in (0, 1)
        ])),
        "class_counts": {str(label): int(np.sum(labels == label)) for label in (0, 1)},
    }
    result["auroc"] = binary_auroc(probs[labels == 0, 1], probs[labels == 1, 1])
    return result


def _metrics_at_threshold(scores: np.ndarray, labels: np.ndarray, threshold: float) -> dict:
    probs = np.column_stack((1.0 - scores, scores))
    # _metrics uses argmax at 0.5, so shift scores while retaining their ranking.
    shifted = np.clip(0.5 + (scores - threshold) / 2.0, 0.0, 1.0)
    result = _metrics(np.column_stack((1.0 - shifted, shifted)), labels)
    result["auroc"] = binary_auroc(scores[labels == 0], scores[labels == 1])
    result["decision_threshold"] = threshold
    return result


def _select_training_only(
    features: np.ndarray,
    labels: np.ndarray,
    sample_ids: list[str],
) -> tuple[float, float, dict]:
    train_indices: list[int] = []
    validation_indices: list[int] = []
    for label in (0, 1):
        indices = np.flatnonzero(labels == label).tolist()
        indices.sort(key=lambda index: hashlib.sha256(sample_ids[index].encode()).hexdigest())
        if len(indices) != 50:
            raise ValueError("expected 50 HeLa cells per aligned class")
        train_indices.extend(indices[:30])
        validation_indices.extend(indices[30:])
    train_indices = sorted(train_indices)
    validation_indices = sorted(validation_indices)
    best: tuple[tuple[float, float, float, float], float, float, dict] | None = None
    for penalty in np.geomspace(0.001, 1000.0, 25):
        model = fit_ridge(features[train_indices], labels[train_indices], 2, float(penalty))
        scores = probabilities(logits(model, features[validation_indices]), 1.0)[:, 1]
        for threshold in np.linspace(0.1, 0.9, 161):
            metrics = _metrics_at_threshold(scores, labels[validation_indices], float(threshold))
            key = (
                metrics["balanced_accuracy"],
                metrics["macro_f1"],
                -abs(float(threshold) - 0.5),
                -float(penalty),
            )
            if best is None or key > best[0]:
                best = (key, float(penalty), float(threshold), metrics)
    assert best is not None
    return best[1], best[2], best[3]


def evaluate(train_manifest: Path, external_manifest: Path) -> dict:
    x_train, y_train, train_conditions, _ = _curves(train_manifest, "hela")
    x_external, y_external, external_conditions, _ = _curves(external_manifest, "external")
    model = fit_ridge(x_train, y_train, 2, penalty=1.0)
    external_probs = probabilities(logits(model, x_external), temperature=1.0)
    metrics = _metrics(external_probs, y_external)

    aligned_train, aligned_y, _, train_sample_ids = _invariant_features(
        train_manifest, "hela", active_only=False
    )
    aligned_external, aligned_external_y, aligned_conditions, _ = _invariant_features(
        external_manifest, "external", active_only=True
    )
    penalty, threshold, validation_metrics = _select_training_only(
        aligned_train, aligned_y, train_sample_ids
    )
    aligned_model = fit_ridge(aligned_train, aligned_y, 2, penalty)
    aligned_scores = probabilities(logits(aligned_model, aligned_external), 1.0)[:, 1]
    aligned_metrics = _metrics_at_threshold(aligned_scores, aligned_external_y, threshold)
    stratum_metrics = {}
    for stratum in ACTIVE_NANOSWITCH_STRATA:
        selected = np.asarray([condition.startswith(stratum) for condition in aligned_conditions])
        stratum_metrics[stratum] = _metrics_at_threshold(
            aligned_scores[selected], aligned_external_y[selected], threshold
        )
    passed = (
        aligned_metrics["macro_f1"] >= 0.70
        and aligned_metrics["balanced_accuracy"] >= 0.70
        and aligned_metrics["auroc"] >= 0.70
    )
    return {
        "task": "binary_Piezo1_inhibitor_signature_transfer",
        "training_dataset": "zenodo-18495219 HeLa",
        "external_dataset": "zenodo-19890985 EA.hy926",
        "external_lab_holdout": True,
        "split_unit": "entire dataset and lab",
        "sample_overlap_count": 0,
        "training_sample_count": len(y_train),
        "training_conditions": sorted(set(train_conditions)),
        "external_conditions": sorted(set(external_conditions)),
        "aligned_external_conditions": sorted(set(aligned_conditions)),
        "label_mapping": {
            "0": "non-inhibited magnetic/mechanical reference",
            "1": "GsMTx4 Piezo1-inhibited",
        },
        "metrics": metrics,
        "legacy_broad_arm_metrics": metrics,
        "aligned_active_nanoswitch_task": {
            "selection_protocol": (
                "feature, ridge penalty, and decision threshold selected on a deterministic "
                "30/20-per-class HeLa train/validation split; EA.hy926 never used for selection"
            ),
            "representation": "baseline_subtracted_time_and_amplitude_invariant_curve_shape",
            "active_strata": list(ACTIVE_NANOSWITCH_STRATA),
            "selected_penalty": penalty,
            "selected_threshold": threshold,
            "training_validation_metrics": validation_metrics,
            "external_metrics": aligned_metrics,
            "external_stratum_metrics": stratum_metrics,
        },
        "task_holdout_passed": passed,
        "limitations": [
            "cell line and stimulation protocol differ across labs",
            "EA.hy926 biological-replicate IDs are not encoded per cell",
            "mechanism-level labels are aligned; exact treatment labels are not identical",
            "external class prevalence differs sharply and no external labels tune the threshold",
        ],
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train-manifest", type=Path, required=True)
    parser.add_argument("--external-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(args.train_manifest, args.external_manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result["metrics"], sort_keys=True))


if __name__ == "__main__":
    main()
