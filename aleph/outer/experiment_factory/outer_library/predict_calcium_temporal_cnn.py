#!/usr/bin/env python3
"""Read-only inference for the rejected-but-reproducible calcium CNN.

Predictions are diagnostic only.  This module cannot write Aleph parameters,
and a prediction never upgrades the model's failed external validation status.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from train_calcium_temporal_cnn import (
    KERNEL_SIZE,
    LENGTH,
    SCHEMA,
    _embedding,
    _probability,
    preprocess_curve,
)


def load_checkpoint(path: Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with np.load(path, allow_pickle=False) as payload:
        if str(payload["schema"]) != SCHEMA:
            raise ValueError("unsupported calcium temporal CNN schema")
        if int(payload["input_length"]) != LENGTH:
            raise ValueError("checkpoint input length does not match runtime")
        if int(payload["kernel_size"]) != KERNEL_SIZE:
            raise ValueError("checkpoint kernel size does not match runtime")
        model = {
            name: np.asarray(payload[name], dtype=np.float64)
            for name in (
                "filters", "bias", "classifier", "intercept",
                "embedding_mean", "embedding_scale", "class_centroids",
            )
        }
        class_names = [str(value) for value in payload["class_names"].tolist()]
    filter_count = model["filters"].shape[0]
    expected_shapes = {
        "filters": (filter_count, KERNEL_SIZE),
        "bias": (filter_count,),
        "classifier": (2 * filter_count, 2),
        "intercept": (2,),
        "embedding_mean": (2 * filter_count,),
        "embedding_scale": (2 * filter_count,),
        "class_centroids": (2, 2 * filter_count),
    }
    for name, shape in expected_shapes.items():
        if model[name].shape != shape:
            raise ValueError(f"invalid checkpoint shape for {name}: {model[name].shape}")
    if class_names != ["non_inhibited_reference", "Piezo1_inhibited"]:
        raise ValueError("checkpoint class names do not match runtime task")
    return model, {
        "sha256": digest,
        "schema": SCHEMA,
        "class_names": class_names,
        "filter_count": filter_count,
    }


def load_manifest_curves(
    path: Path,
    *,
    allowed_conditions: tuple[str, ...] = (),
) -> tuple[np.ndarray, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["observable"] not in {
                "intracellular_calcium_F_over_F0",
                "intracellular_calcium_delta_F_over_F0",
                "intracellular_calcium_normalized_fluorescence",
            }:
                continue
            grouped[row["sample_id"]].append(row)
    curves, metadata = [], []
    for sample_id, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: row["coordinate_value"])
        condition = rows[0]["condition"]
        if allowed_conditions and condition not in allowed_conditions:
            continue
        coordinates = np.asarray(
            [row["coordinate_value"] for row in rows], dtype=np.float64
        )
        if len(rows) < KERNEL_SIZE or not np.all(np.diff(coordinates) > 0):
            raise ValueError(f"invalid calcium curve for {sample_id}")
        values = np.asarray([row["value"] for row in rows], dtype=np.float64)
        curves.append(preprocess_curve(values))
        metadata.append({
            "sample_id": sample_id,
            "dataset_id": rows[0]["dataset_id"],
            "lab_group": rows[0]["lab_group"],
            "condition": condition,
        })
    if not curves:
        raise ValueError("manifest contains no selected calcium curves")
    return np.stack(curves), metadata


def predict(
    checkpoint: Path,
    manifest: Path,
    *,
    allowed_conditions: tuple[str, ...] = (),
) -> list[dict[str, Any]]:
    model, checkpoint_info = load_checkpoint(checkpoint)
    curves, metadata = load_manifest_curves(
        manifest, allowed_conditions=allowed_conditions
    )
    probabilities = _probability(model, curves)
    embedding = _embedding(model, curves)
    standardized = (
        embedding - model["embedding_mean"]
    ) / model["embedding_scale"]
    distance = (
        (standardized[:, None] - model["class_centroids"][None]) ** 2
    ).mean(axis=2).min(axis=1)
    rows = []
    for item, probability, ood_distance in zip(
        metadata, probabilities, distance
    ):
        predicted = int(np.argmax(probability))
        rows.append({
            **item,
            "model_schema": SCHEMA,
            "model_checkpoint_sha256": checkpoint_info["sha256"],
            "probability_non_inhibited_reference": float(probability[0]),
            "probability_Piezo1_inhibited": float(probability[1]),
            "predicted_class": checkpoint_info["class_names"][predicted],
            "embedding_ood_distance": float(ood_distance),
            "validation_status": "rejected_external_classification",
            "aleph_authority": "none",
            "may_select_aleph_parameter": False,
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--condition", action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = predict(
        args.checkpoint, args.manifest,
        allowed_conditions=tuple(args.condition),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    print(json.dumps({
        "predictions": len(rows),
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
