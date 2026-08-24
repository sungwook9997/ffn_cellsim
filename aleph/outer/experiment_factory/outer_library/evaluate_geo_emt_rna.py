#!/usr/bin/env python3
"""Evaluate replicated, study-separated RNA direction constraints for EMT sweeps."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


CONTRASTS = {
    "GSE17708": ("untreated", "TGFB1_72h"),
    "GSE42373": ("3d_control", "3d_treated"),
    "GSE125369": ("control", "TGFB3_21d"),
    "GSE69667": ("untreated", "TGFB1_96h"),
    "GSE49644": ("A549_parental", "A549_TGFB_EMT_3week"),
}
EXPECTED = {"CDH1": "decrease", "VIM": "increase"}


def _direction(value: float) -> str:
    return "increase" if value > 0 else "decrease" if value < 0 else "tie"


def evaluate(tensor_path: Path, manifest_report_path: Path, protein_report_path: Path) -> dict[str, object]:
    manifest = json.loads(manifest_report_path.read_text(encoding="utf-8"))
    protein = json.loads(protein_report_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "aleph.outer_library.geo_emt_rna_manifest.v1":
        raise ValueError("GEO EMT manifest schema changed")
    if protein.get("schema") != "aleph.outer_library.karacosta_emt_cytof_evaluation.v1":
        raise ValueError("protein EMT report schema changed")
    with np.load(tensor_path, allow_pickle=False) as data:
        study, condition, gene = (data[name].astype(str) for name in ("study", "condition", "gene"))
        value = data["value"].astype(np.float64)
        sample = data["sample_id"].astype(str)
    rows = []
    for study_name, (control_name, treated_name) in CONTRASTS.items():
        for gene_name in EXPECTED:
            control = value[(study == study_name) & (condition == control_name) & (gene == gene_name)]
            treated = value[(study == study_name) & (condition == treated_name) & (gene == gene_name)]
            if len(control) < 2 or len(treated) < 2:
                raise ValueError(f"biological replication disappeared: {study_name}/{gene_name}")
            delta = float(treated.mean() - control.mean())
            pairwise = (treated[:, None] - control[None, :]).ravel()
            rows.append({
                "study": study_name, "gene": gene_name,
                "control_condition": control_name, "treated_condition": treated_name,
                "control_biological_n": len(control), "treated_biological_n": len(treated),
                "control_mean": float(control.mean()), "treated_mean": float(treated.mean()),
                "treated_minus_control": delta, "direction": _direction(delta),
                "all_between_group_pairwise_differences_have_expected_sign": bool(
                    np.all(pairwise < 0) if EXPECTED[gene_name] == "decrease" else np.all(pairwise > 0)
                ),
                "absolute_scale_transfer_permitted": False,
            })
    all_expected = all(row["direction"] == EXPECTED[row["gene"]] for row in rows)
    separated = all(row["all_between_group_pairwise_differences_have_expected_sign"] for row in rows)
    signature_successes = sum(
        all(
            row["direction"] == EXPECTED[row["gene"]]
            for row in rows if row["study"] == study_name
        )
        for study_name in CONTRASTS
    )
    study_sign_p = 0.5 ** len(CONTRASTS) if signature_successes == len(CONTRASTS) else None

    cross_cell_rows = []
    for cell_line in ("HCC827", "NCI-H358"):
        for gene_name in EXPECTED:
            control = value[
                (study == "GSE49644")
                & (condition == f"{cell_line}_parental")
                & (gene == gene_name)
            ]
            treated = value[
                (study == "GSE49644")
                & (condition == f"{cell_line}_TGFB_EMT_3week")
                & (gene == gene_name)
            ]
            delta = float(treated.mean() - control.mean())
            cross_cell_rows.append({
                "cell_line": cell_line, "gene": gene_name,
                "control_biological_n": len(control),
                "treated_biological_n": len(treated),
                "treated_minus_control": delta, "direction": _direction(delta),
                "expected_direction": EXPECTED[gene_name],
                "direction_agrees": _direction(delta) == EXPECTED[gene_name],
                "independent_study_or_lab": False,
            })

    context_rows = {}
    for gene_name in EXPECTED:
        control = value[(study == "GSE42373") & (condition == "2d_control") & (gene == gene_name)]
        treated = value[(study == "GSE42373") & (condition == "2d_treated") & (gene == gene_name)]
        context_rows[gene_name] = {
            "2d_treated_minus_control": float(treated.mean() - control.mean()),
            "direction": _direction(float(treated.mean() - control.mean())),
            "expected_3d_EMT_direction": EXPECTED[gene_name],
        }
    context_shift = any(context_rows[gene]["direction"] != EXPECTED[gene] for gene in EXPECTED)
    endpoint_signatures = {
        study_name: {
            row["gene"]: row["direction"]
            for row in rows if row["study"] == study_name
        }
        for study_name in CONTRASTS
    }
    in_distribution_scores = {
        study_name: sum(
            direction != EXPECTED[gene_name]
            for gene_name, direction in signature.items()
        )
        for study_name, signature in endpoint_signatures.items()
    }
    ood_score = sum(
        context_rows[gene_name]["direction"] != EXPECTED[gene_name]
        for gene_name in EXPECTED
    )
    refusal_threshold = max(in_distribution_scores.values())
    sealed_ood_refused = ood_score > refusal_threshold

    protein_directions = {
        row["protein_observable"]: row["Karacosta_HCC827_CyTOF_0d_to_10d_TGFB_direction"]
        for row in protein["cross_lab_common_marker_direction"]["rows"]
    }
    rna_protein_agree = (
        EXPECTED["CDH1"] == protein_directions["E_cadherin"]
        and EXPECTED["VIM"] == protein_directions["Vimentin"]
    )
    return {
        "schema": "aleph.outer_library.geo_emt_rna_evaluation.v1",
        "task": "replicated_multi_study_RNA_EMT_direction_constraint_before_Aleph_sweep",
        "design_audit": {
            **manifest["counts"], "outer_holdout_unit": "study_and_lab_group",
            "biological_replicates_not_expression_features": True,
            "cross_platform_absolute_scale_pooling_forbidden": True,
        },
        "endpoint_contrasts": rows,
        "replicated_direction": {
            "all_ten_study_gene_directions_expected": all_expected,
            "all_ten_between_group_ranges_completely_separated": separated,
            "study_gene_contrasts": len(rows),
            "independent_studies": 5, "independent_labs": 5,
            "study_level_signature_successes": signature_successes,
            "nominal_exact_one_sided_sign_test_p_for_five_retrospective_study_signatures": study_sign_p,
            "retrospective_direction_replication_passed": bool(
                study_sign_p is not None and study_sign_p < 0.05
            ),
            "study_selection_was_outcome_blind": False,
            "confirmatory_direction_threshold_0_05_passed": False,
            "confirmatory_failure_reason": (
                "GSE69667 and GSE49644 were selected from literature reporting "
                "successful EMT, so the nominal sign test is subject to study-"
                "selection and publication bias"
            ),
        },
        "within_GSE49644_cross_cell_line_transfer": {
            "rows": cross_cell_rows,
            "all_four_marker_directions_agree": all(
                row["direction_agrees"] for row in cross_cell_rows
            ),
            "counts_as_additional_independent_study_or_lab": False,
        },
        "cross_assay_bridge": {
            "RNA_direction_agrees_with_independent_IF_and_CyTOF_protein_direction": rna_protein_agree,
            "RNA_protein_absolute_scale_transfer_permitted": False,
        },
        "semantic_OOD_challenge": {
            "dataset": "GSE42373", "held_context": "2d_TGFB_TNFA",
            "reference_context": "3d_TGFB_TNFA", "rows": context_rows,
            "at_least_one_marker_direction_changes_by_culture_context": context_shift,
            "challenge_defined": context_shift,
            "sealed_direction_hamming_diagnostic": {
                "training_ID_study_scores": in_distribution_scores,
                "threshold_max_training_ID_score": refusal_threshold,
                "held_2d_context_score": ood_score,
                "held_2d_context_refused": sealed_ood_refused,
                "held_OOD_context_count": 1,
            },
            "diagnostic_passed": sealed_ood_refused,
            "OOD_refusal_model_validated": False,
        },
        "representation_training_eligible": True,
        "directional_pre_sweep_filter_eligible": all_expected and separated and rna_protein_agree,
        "independent_study_direction_holdout": True,
        "study_level_direction_inference_passed": False,
        "numeric_uncertainty_promotion_passed": False,
        "OOD_refusal_validated": False,
        "pre_sweep_role": {
            "purpose": "reduce_and_rank_future_Aleph_EMT_sweep_not_replace_it",
            "supported_action": "reject EMT candidates with CDH1-up or VIM-down endpoint response across compatible A549 TGFB contexts",
            "numeric_range": None, "may_emit_sweep_axis": False,
            "remaining_blockers": [
                "preregistered outcome-blind external EMT study evaluated after protocol freeze",
                "trained and sealed validation of culture-context semantic OOD refusal",
                "Aleph CDH1/VIM-compatible observation operators",
                "monotonic Aleph forward sensitivity linking candidate axes to the observables",
            ],
        },
        "aleph_authority": "none", "may_select_aleph_parameter": False,
        "limitations": [
            "five retrospective study-level EMT signatures give nominal one-sided p=0.03125, but literature-based study selection prevents confirmatory inference",
            "RNA expression constrains marker direction but is not interchangeable with protein abundance",
            "GSE42373 combines TGF-beta and TNF-alpha and contributes the only sealed OOD context, so general OOD validation remains unproven",
            "the data can compress a future sweep but cannot issue a numeric Aleph parameter range",
        ],
        "sample_identity_count": len(set(zip(study, sample))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--protein-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor, args.manifest_report, args.protein_report)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"directional_filter": report["directional_pre_sweep_filter_eligible"],
                      "study_direction_inference": report["study_level_direction_inference_passed"],
                      "numeric_uncertainty_promotion": report["numeric_uncertainty_promotion_passed"]}, sort_keys=True))


if __name__ == "__main__":
    main()
