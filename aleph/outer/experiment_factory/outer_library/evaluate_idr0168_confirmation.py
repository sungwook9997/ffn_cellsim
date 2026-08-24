#!/usr/bin/env python3
"""Evaluate IDR0168 once against frozen localization and OOD heads."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr
import torch

from train_hpa_raw_if_cnn import sha256
from train_multiprovider_if_model import MultiproviderIFCNN, auroc
from train_multiprovider_if_ood_v2 import (
    FAMILIES,
    extract_features,
    probabilities,
    variants,
)


def family_metrics(
    scores: np.ndarray, families: np.ndarray
) -> dict[str, Any]:
    clean = families == "clean"
    output = {}
    aucs = []
    for family in FAMILIES:
        corrupted = families == family
        selected = clean | corrupted
        truth = corrupted[selected].astype(np.float64)
        value = auroc(truth, scores[selected])
        output[family] = {
            "AUROC": value,
            "clean_fields": int(clean.sum()),
            "corrupted_fields": int(corrupted.sum()),
        }
        aucs.append(value)
    return {
        "minimum_family_AUROC": float(min(aucs)),
        "macro_family_AUROC": float(np.mean(aucs)),
        "family": output,
    }


def correlation(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    result = spearmanr(x, y)
    statistic = float(result.statistic)
    pvalue = float(result.pvalue)
    if not np.isfinite(statistic) or not np.isfinite(pvalue):
        raise ValueError("IDR0168 Spearman endpoint is non-finite")
    return {"Spearman_rho": statistic, "two_sided_pvalue": pvalue}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--tensor-report", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--model-report", type=Path, required=True)
    parser.add_argument("--quality-head", type=Path, required=True)
    parser.add_argument("--quality-report", type=Path, required=True)
    parser.add_argument("--quality-protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    tensor_report = json.loads(args.tensor_report.read_text(encoding="utf-8"))
    model_report = json.loads(args.model_report.read_text(encoding="utf-8"))
    quality_report = json.loads(args.quality_report.read_text(encoding="utf-8"))
    frozen = protocol["frozen_model_contract"]
    paths = {
        "checkpoint_sha256": args.checkpoint,
        "model_report_sha256": args.model_report,
        "quality_head_sha256": args.quality_head,
        "quality_report_sha256": args.quality_report,
        "quality_protocol_sha256": args.quality_protocol,
    }
    for key, path in paths.items():
        if sha256(path) != frozen[key]:
            raise ValueError(f"frozen IDR0168 model artifact changed: {key}")
    if tensor_report["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("IDR0168 tensor report does not pin protocol")
    if tensor_report["tensor_sha256"] != sha256(args.tensor):
        raise ValueError("IDR0168 tensor report does not pin tensor")
    if not quality_report["all_v2_development_gates_passed"]:
        raise ValueError("frozen v2 quality head did not pass development gates")

    with np.load(args.tensor, allow_pickle=False) as source:
        data = {key: source[key].copy() for key in source.files}
    images = data["images"]
    fractions = data["field_level_intranuclear_target_fraction"].astype(np.float64)
    genes = data["gene_symbol"].astype(str)
    cell_lines = data["cell_line"].astype(str)
    unseen = data["exact_gene_unseen_in_development"].astype(bool)
    unit_ids = data["source_name"].astype(str)
    if images.shape != (45, 2, 128, 128) or unseen.sum() != 25:
        raise ValueError("IDR0168 frozen tensor shape or unseen stratum changed")
    if len(set(unit_ids)) != 45:
        raise ValueError("IDR0168 inference units are not unique")

    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    classes = [str(value) for value in checkpoint["classes"]]
    if classes != model_report["classes"]:
        raise ValueError("IDR0168 checkpoint/report class order changed")
    model = MultiproviderIFCNN(len(classes))
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(8)
    with torch.no_grad():
        logits, _ = model(torch.from_numpy(images).float().div_(255.0))
    temperature = float(model_report["calibration"]["temperature"])
    calibrated = torch.sigmoid(logits / temperature).cpu().numpy()
    nuclear_score = calibrated[:, classes.index("nuclear_interior")].astype(np.float64)
    if not np.isfinite(nuclear_score).all():
        raise ValueError("IDR0168 nuclear probabilities are non-finite")

    endpoints = protocol["localization_endpoints"]
    all_fields = correlation(nuclear_score, fractions)
    unseen_fields = correlation(nuclear_score[unseen], fractions[unseen])
    high = fractions > (2.0 / 3.0)
    low = fractions < (1.0 / 3.0)
    supported = high | low
    diagnostic_auroc = (
        auroc(high[supported].astype(np.float64), nuclear_score[supported])
        if high.sum() and low.sum() else float("nan")
    )

    positive_gene = endpoints["author_anchor_positive_gene"]
    negative_gene = endpoints["author_anchor_negative_gene"]
    anchor_model_wins = 0
    anchor_observation_wins = 0
    anchor_pairs = []
    for cell_line in sorted(set(cell_lines)):
        positive = np.flatnonzero((cell_lines == cell_line) & (genes == positive_gene))
        negative = np.flatnonzero((cell_lines == cell_line) & (genes == negative_gene))
        if len(positive) != 1 or len(negative) != 1:
            raise ValueError(f"IDR0168 anchor pair is not unique: {cell_line}")
        p, n = int(positive[0]), int(negative[0])
        model_win = bool(nuclear_score[p] > nuclear_score[n])
        observation_win = bool(fractions[p] > fractions[n])
        anchor_model_wins += int(model_win)
        anchor_observation_wins += int(observation_win)
        anchor_pairs.append({
            "cell_line": cell_line,
            "positive_gene": positive_gene,
            "negative_gene": negative_gene,
            "positive_nuclear_score": float(nuclear_score[p]),
            "negative_nuclear_score": float(nuclear_score[n]),
            "positive_observed_fraction": float(fractions[p]),
            "negative_observed_fraction": float(fractions[n]),
            "model_score_win": model_win,
            "observed_fraction_win": observation_win,
        })

    units = [
        {"provider": "IDR0168", "unit_id": unit, "image": image}
        for unit, image in zip(unit_ids, images)
    ]
    salt = json.loads(args.quality_protocol.read_text(encoding="utf-8"))[
        "development_unit_contract"
    ]["selection_salt"]
    variant_images, _, _, variant_families = variants(units, salt)
    features = extract_features(model, variant_images)
    with np.load(args.quality_head, allow_pickle=False) as source:
        quality = {key: source[key].copy() for key in source.files}
    quality_scores = probabilities(
        features,
        quality["feature_mean"],
        quality["feature_scale"],
        quality["parameter"],
    )
    OOD = family_metrics(quality_scores, variant_families)
    clean = variant_families == "clean"
    threshold = float(quality["clean_refusal_threshold"])
    clean_false_refusal = float(np.mean(quality_scores[clean] > threshold))
    clean_quality_scores = quality_scores[clean]
    if len(clean_quality_scores) != len(images):
        raise ValueError("IDR0168 clean quality scores do not align with fields")

    field_predictions = []
    for index in range(len(images)):
        field_predictions.append({
            "source_name": unit_ids[index],
            "image_id": str(data["image_id"][index]),
            "gene_symbol": genes[index],
            "cell_line": cell_lines[index],
            "exact_gene_unseen_in_development": bool(unseen[index]),
            "observed_intranuclear_target_fraction": float(fractions[index]),
            "frozen_model_nuclear_interior_probability": float(
                nuclear_score[index]
            ),
            "clean_imaging_OOD_probability": float(clean_quality_scores[index]),
            "clean_imaging_refused": bool(
                clean_quality_scores[index] > threshold
            ),
        })

    OOD_endpoints = protocol["imaging_OOD_endpoints"]
    gates = {
        "all_45_field_Spearman": (
            all_fields["Spearman_rho"]
            >= endpoints["primary_all_45_spearman_minimum"]
        ),
        "unseen_gene_25_field_Spearman": (
            unseen_fields["Spearman_rho"]
            >= endpoints["primary_unseen_gene_25_spearman_minimum"]
        ),
        "minimum_supported_nuclear_fields": (
            int(high.sum()) >= endpoints["minimum_supported_nuclear_diagnostic_fields"]
        ),
        "minimum_supported_non_nuclear_fields": (
            int(low.sum()) >= endpoints["minimum_supported_non_nuclear_diagnostic_fields"]
        ),
        "diagnostic_high_vs_low_AUROC": bool(
            np.isfinite(diagnostic_auroc)
            and diagnostic_auroc >= endpoints["diagnostic_high_vs_low_AUROC_minimum"]
        ),
        "author_anchor_model_score_paired_wins": (
            anchor_model_wins
            >= endpoints["author_anchor_model_score_paired_cell_line_wins_minimum"]
        ),
        "author_anchor_observed_fraction_paired_wins": (
            anchor_observation_wins
            >= endpoints["author_anchor_observed_fraction_paired_cell_line_wins_minimum"]
        ),
        "external_minimum_corruption_family_AUROC": (
            OOD["minimum_family_AUROC"]
            >= OOD_endpoints["external_minimum_family_AUROC"]
        ),
        "external_macro_corruption_family_AUROC": (
            OOD["macro_family_AUROC"]
            >= OOD_endpoints["external_macro_family_AUROC"]
        ),
        "maximum_clean_false_refusal_rate": (
            clean_false_refusal
            <= OOD_endpoints["maximum_clean_false_refusal_rate"]
        ),
        "all_fields_predictions_and_quality_scores_finite": bool(
            np.isfinite(calibrated).all()
            and np.isfinite(fractions).all()
            and np.isfinite(quality_scores).all()
        ),
    }
    passed = all(gates.values())
    report = {
        "schema": "aleph.outer_library.idr0168_confirmation_evaluation.v1",
        "protocol_sha256": sha256(args.protocol),
        "tensor_sha256": sha256(args.tensor),
        "tensor_report_sha256": sha256(args.tensor_report),
        "frozen_artifact_hashes": frozen,
        "field_count": len(images),
        "unseen_gene_field_count": int(unseen.sum()),
        "seen_gene_field_count": int((~unseen).sum()),
        "localization": {
            "all_45_fields": all_fields,
            "unseen_gene_25_fields": unseen_fields,
            "nuclear_diagnostic_field_count": int(high.sum()),
            "intermediate_diagnostic_field_count": int((~supported).sum()),
            "non_nuclear_diagnostic_field_count": int(low.sum()),
            "high_vs_low_AUROC": diagnostic_auroc,
            "author_anchor_model_score_paired_cell_line_wins": anchor_model_wins,
            "author_anchor_observed_fraction_paired_cell_line_wins": anchor_observation_wins,
            "author_anchor_pairs": anchor_pairs,
        },
        "imaging_OOD": {
            **OOD,
            "clean_false_refusal_rate": clean_false_refusal,
            "frozen_clean_refusal_threshold": threshold,
        },
        "field_predictions": field_predictions,
        "gates": gates,
        "all_preregistered_gates_passed": passed,
        "status": (
            "validated_external_nuclear_localization_subset"
            if passed else "failed_and_permanently_consumed"
        ),
        "independent_provider_nuclear_localization_transfer": passed,
        "external_imaging_contract_OOD_for_IDR0168_like_inputs": passed,
        "all_nine_class_external_validation": False,
        "external_numeric_uncertainty_authority": False,
        "nuclear_localization_external_evidence_eligible": passed,
        "general_cell_state_model_ready": False,
        "production_eligible": False,
        "Aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
        "may_remove_sweep_axis": False,
        "model_predictions_made": len(images),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "all_gates_passed": passed,
        "all_rho": all_fields["Spearman_rho"],
        "unseen_rho": unseen_fields["Spearman_rho"],
        "diagnostic_AUROC": diagnostic_auroc,
        "OOD_minimum_AUROC": OOD["minimum_family_AUROC"],
        "OOD_macro_AUROC": OOD["macro_family_AUROC"],
        "clean_false_refusal": clean_false_refusal,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
