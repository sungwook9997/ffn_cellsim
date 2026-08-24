#!/usr/bin/env python3
"""Post-freeze diagnosis of the failed multiprovider IF development OOD gate."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from train_hpa_raw_if_cnn import sha256
from train_multiprovider_if_model import (
    MultiproviderIFCNN,
    aggregate,
    auroc,
    predict,
)


def summary(values: np.ndarray) -> dict[str, float]:
    return {
        "minimum": float(values.min()),
        "q25": float(np.quantile(values, 0.25)),
        "median": float(np.median(values)),
        "q75": float(np.quantile(values, 0.75)),
        "maximum": float(values.max()),
        "mean": float(values.mean()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--model-report", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--tensor-report", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    model_report = json.loads(args.model_report.read_text(encoding="utf-8"))
    tensor_report = json.loads(args.tensor_report.read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if model_report["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("model report does not pin training protocol")
    if model_report["checkpoint_sha256"] != sha256(args.checkpoint):
        raise ValueError("model report does not pin checkpoint")
    if tensor_report["tensor_sha256"] != sha256(args.tensor):
        raise ValueError("tensor report does not pin tensor")
    if model_report["all_development_open_gates_passed"]:
        raise ValueError("diagnostic is only valid for the frozen failed gate")
    if model_report["next_external_confirmation_pixels_may_be_opened"]:
        raise ValueError("failed model unexpectedly permits external pixel access")

    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model = MultiproviderIFCNN(len(checkpoint["classes"]))
    model.load_state_dict(checkpoint["state_dict"])
    with np.load(args.tensor, allow_pickle=False) as source:
        tensor = {key: source[key].copy() for key in source.files}
    split = tensor["split"].astype(str)
    targets = tensor["target_name"].astype(str)
    centroids = np.asarray(checkpoint["OOD_class_centroids"], dtype=np.float64)
    variances = np.asarray(
        checkpoint["OOD_class_diagonal_variances"], dtype=np.float64
    )

    def role_scores(role: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        mask = split == role
        logits, embeddings = predict(model, tensor["images"][mask])
        _, unit_embedding, _, unit_names = aggregate(
            logits, embeddings, tensor["labels"][mask].astype(np.float32), targets[mask]
        )
        distances = ((
            unit_embedding[:, None, :] - centroids[None, :, :]
        ) ** 2 / variances[None, :, :]).mean(axis=2)
        return unit_names, distances.min(axis=1), distances.argmin(axis=1)

    in_names, in_scores, _ = role_scores("selection")
    out_names, out_scores, nearest = role_scores("semantic_OOD")
    reproduced_auroc = auroc(
        np.concatenate((np.zeros(len(in_scores)), np.ones(len(out_scores)))),
        np.concatenate((in_scores, out_scores)),
    )
    if abs(reproduced_auroc - model_report["OOD"]["development_AUROC"]) > 1e-7:
        raise ValueError(
            "post-freeze diagnostic did not reproduce OOD AUROC: "
            f"{reproduced_auroc} != {model_report['OOD']['development_AUROC']}"
        )
    by_target = {row["target_name"]: row for row in manifest["rows"]}
    target_rows: list[dict[str, Any]] = []
    for target, score, class_index in sorted(
        zip(out_names.tolist(), out_scores.tolist(), nearest.tolist())
    ):
        row = by_target[target]
        target_rows.append({
            "target_name": target,
            "OOD_score": score,
            "nearest_frozen_coarse_class": checkpoint["classes"][class_index],
            "grade_3_annotations": row["grade_3_annotations"],
            "grade_2_annotations": row["grade_2_annotations"],
            "grade_1_annotations_not_used_for_role_assignment": row["grade_1_annotations"],
        })
    high_grade = Counter(
        annotation
        for row in target_rows
        for annotation in row["grade_3_annotations"] + row["grade_2_annotations"]
    )
    nearest_counts = Counter(row["nearest_frozen_coarse_class"] for row in target_rows)
    grade1_in_domain = sum(bool(row["grade_1_annotations_not_used_for_role_assignment"])
                           for row in target_rows)
    output = {
        "schema": "aleph.outer_library.multiprovider_if_ood_failure_diagnostic.v1",
        "analysis_role": "post-freeze_failure_diagnosis_only",
        "may_change_frozen_gate_or_checkpoint": False,
        "protocol_sha256": sha256(args.protocol),
        "checkpoint_sha256": sha256(args.checkpoint),
        "model_report_sha256": sha256(args.model_report),
        "tensor_sha256": sha256(args.tensor),
        "tensor_report_sha256": sha256(args.tensor_report),
        "manifest_sha256": sha256(args.manifest),
        "reproduced_development_OOD_AUROC": reproduced_auroc,
        "in_domain_selection_target_count": len(in_names),
        "semantic_OOD_target_count": len(out_names),
        "in_domain_score_summary": summary(in_scores),
        "semantic_OOD_score_summary": summary(out_scores),
        "semantic_OOD_high_grade_author_annotation_counts": dict(sorted(high_grade.items())),
        "semantic_OOD_nearest_frozen_coarse_class_counts": dict(sorted(nearest_counts.items())),
        "semantic_OOD_targets_with_nonempty_grade1_annotation_count": grade1_in_domain,
        "semantic_OOD_targets": target_rows,
        "diagnosis": (
            "the frozen semantic-OOD role encoded unmapped high-grade annotation strings, "
            "not necessarily visually out-of-distribution images; nuclear_punctae targets "
            "were closest to nuclear_interior and big_aggregates targets were closest to "
            "vesicular_organelle, so the image-distance score correctly treated them as "
            "coarse-ontology-like rather than rejecting them"
        ),
        "frozen_development_gate_remains_failed": True,
        "next_external_confirmation_pixels_may_be_opened": False,
        "external_provider_validation_complete": False,
        "production_eligible": False,
        "uncertainty_authority": False,
        "observable_evidence_eligible": False,
        "Aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.output.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "AUROC": reproduced_auroc,
        "in_domain_median": output["in_domain_score_summary"]["median"],
        "semantic_OOD_median": output["semantic_OOD_score_summary"]["median"],
        "high_grade_annotations": output["semantic_OOD_high_grade_author_annotation_counts"],
        "nearest_classes": output["semantic_OOD_nearest_frozen_coarse_class_counts"],
        "external_pixels_may_open": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
