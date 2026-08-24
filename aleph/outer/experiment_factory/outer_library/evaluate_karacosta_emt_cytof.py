#!/usr/bin/env python3
"""Evaluate an independent-lab EMT marker-direction bridge without pseudoreplication."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


TIMEPOINT_LABELS = {
    1: "0d_control", 2: "2d_TGFB", 3: "6d_TGFB", 4: "10d_TGFB",
    5: "2d_withdrawal", 6: "6d_withdrawal", 7: "10d_withdrawal",
}
COMMON_MARKERS = {
    "CD324": "E_cadherin",
    "Vimentin": "Vimentin",
}


def _direction(delta: float) -> str:
    return "increase" if delta > 0 else "decrease" if delta < 0 else "tie"


def evaluate(
    tensor_path: Path, manifest_report_path: Path, exif_report_path: Path
) -> dict[str, object]:
    manifest = json.loads(manifest_report_path.read_text(encoding="utf-8"))
    exif = json.loads(exif_report_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "aleph.outer_library.karacosta_emt_cytof_manifest.v1":
        raise ValueError("Karacosta manifest schema changed")
    if exif.get("schema") != "aleph.outer_library.figshare_exif_emt_evaluation.v1":
        raise ValueError("ExIF EMT report schema changed")
    with np.load(tensor_path, allow_pickle=False) as data:
        raw_X = data["raw_X"].astype(np.float64)
        raw_timepoint = data["raw_timepoint_id"].astype(np.int8)
        raw_cluster = data["raw_cluster_id"].astype(np.int8)
        state_X = data["state_X"].astype(np.float64)
        state_label = data["state_label"].astype(str)
        feature_names = data["feature_names"].astype(str)
        split_group = data["split_group"].astype(str)
    if raw_X.shape != (90_066, 6) or state_X.shape != (29_066, 6):
        raise ValueError("Karacosta tensor shape changed")
    if len(set(split_group)) != 1:
        raise ValueError("Karacosta cells escaped the mandatory experiment group")
    feature_index = {name: index for index, name in enumerate(feature_names)}

    timepoint_stats: dict[str, dict[str, dict[str, float | int]]] = {}
    for timepoint, label in TIMEPOINT_LABELS.items():
        mask = raw_timepoint == timepoint
        marker_stats = {}
        for marker in COMMON_MARKERS:
            values = raw_X[mask, feature_index[marker]]
            marker_stats[marker] = {
                "n_single_cell_events_descriptive_only": int(len(values)),
                "mean_nontransformed_intensity": float(values.mean()),
                "median_nontransformed_intensity": float(np.median(values)),
                "q25": float(np.quantile(values, 0.25)),
                "q75": float(np.quantile(values, 0.75)),
            }
        timepoint_stats[label] = marker_stats

    progressive_rows = []
    withdrawal_rows = []
    for marker, semantic_name in COMMON_MARKERS.items():
        control = timepoint_stats["0d_control"][marker]["median_nontransformed_intensity"]
        day10 = timepoint_stats["10d_TGFB"][marker]["median_nontransformed_intensity"]
        withdrawal = timepoint_stats["10d_withdrawal"][marker]["median_nontransformed_intensity"]
        progressive_rows.append({
            "marker": semantic_name,
            "control_median": control,
            "10d_TGFB_median": day10,
            "10d_TGFB_minus_control": day10 - control,
            "direction": _direction(day10 - control),
        })
        withdrawal_rows.append({
            "marker": semantic_name,
            "10d_TGFB_median": day10,
            "10d_withdrawal_median": withdrawal,
            "withdrawal_minus_10d_TGFB": withdrawal - day10,
            "direction": _direction(withdrawal - day10),
        })

    transformed_state_means = {}
    for state in sorted(set(state_label)):
        mask = state_label == state
        transformed_state_means[state] = {
            marker: float(state_X[mask, feature_index[marker]].mean())
            for marker in COMMON_MARKERS
        }
    state_relation = {
        "epithelial_state_pool": ["E1", "E2", "E3"],
        "mesenchymal_state": "M",
        "E_state_unweighted_mean_of_author_state_means": {},
        "M_state_mean": {},
        "M_minus_E_direction": {},
    }
    for marker in COMMON_MARKERS:
        epithelial = float(np.mean([
            transformed_state_means[state][marker] for state in ("E1", "E2", "E3")
        ]))
        mesenchymal = transformed_state_means["M"][marker]
        state_relation["E_state_unweighted_mean_of_author_state_means"][marker] = epithelial
        state_relation["M_state_mean"][marker] = mesenchymal
        state_relation["M_minus_E_direction"][marker] = _direction(mesenchymal - epithelial)

    exif_rows = {
        row["marker"]: row
        for row in exif["TGF_beta1_marker_direction_check"]["rows"]
    }
    cross_lab_rows = []
    for row in progressive_rows:
        exif_name = row["marker"]
        exif_row = exif_rows[exif_name]
        cross_lab_rows.append({
            "protein_observable": exif_name,
            "Karacosta_HCC827_CyTOF_0d_to_10d_TGFB_direction": row["direction"],
            "ExIF_A549_fixed_IF_control_to_48h_TGFB_direction": exif_row["observed_direction"],
            "directions_agree": row["direction"] == exif_row["observed_direction"],
            "absolute_intensity_transfer_permitted": False,
        })
    cross_lab_passed = all(row["directions_agree"] for row in cross_lab_rows)
    progressive_expected = {
        row["marker"]: row["direction"] for row in progressive_rows
    } == {"E_cadherin": "decrease", "Vimentin": "increase"}
    withdrawal_expected = {
        row["marker"]: row["direction"] for row in withdrawal_rows
    } == {"E_cadherin": "increase", "Vimentin": "decrease"}

    return {
        "schema": "aleph.outer_library.karacosta_emt_cytof_evaluation.v1",
        "task": "independent_lab_common_protein_EMT_direction_and_state_bridge",
        "source": manifest["source"],
        "design_audit": {
            "raw_single_cell_events": manifest["counts"]["raw_single_cell_events"],
            "transformed_state_subset_cells": manifest["counts"]["transformed_downsampled_state_cells"],
            "mandatory_split_groups": len(set(split_group)),
            "single_cells_are_not_independent_biological_replicates": True,
            "timepoints_are_not_independent_biological_replicates": True,
            "released_biological_replicate_identifiers": False,
            "author_states_were_derived_from_the_same_six_marker_features": True,
            "cell_state_classifier_on_these_features_is_label_circular_and_forbidden": True,
        },
        "raw_timepoint_descriptive_statistics": timepoint_stats,
        "progressive_EMT_0d_to_10d_TGFB": {
            "rows": progressive_rows,
            "E_cadherin_down_Vimentin_up_reproduced": progressive_expected,
            "inferential_test_permitted": False,
        },
        "withdrawal_10d_TGFB_to_10dW": {
            "rows": withdrawal_rows,
            "partial_marker_reversal_reproduced": withdrawal_expected,
            "inferential_test_permitted": False,
        },
        "author_transformed_state_relation": state_relation,
        "cross_lab_common_marker_direction": {
            "discovery_dataset": exif["source"]["dataset_id"],
            "discovery_lab": exif["source"]["lab_group"],
            "external_dataset": manifest["source"]["dataset_id"],
            "external_lab": manifest["source"]["lab_group"],
            "dataset_and_lab_disjoint": (
                exif["source"]["dataset_id"] != manifest["source"]["dataset_id"]
                and exif["source"]["lab_group"] != manifest["source"]["lab_group"]
            ),
            "cell_type_shift": "A549_to_HCC827",
            "modality_shift": "fixed_IF_to_CyTOF",
            "rows": cross_lab_rows,
            "independent_lab_common_marker_direction_reproduced": cross_lab_passed,
            "numeric_scale_transfer_permitted": False,
        },
        "state_axis": "EMT_epithelial_partial_mesenchymal_MET_state",
        "representation_training_eligible": True,
        "directional_bridge_eligible": cross_lab_passed,
        "observable_evidence_eligible": False,
        "independent_dataset_direction_holdout": cross_lab_passed,
        "independent_lab_direction_holdout": cross_lab_passed,
        "independent_lab_inferential_holdout": False,
        "uncertainty_holdout_passed": False,
        "OOD_refusal_validated": False,
        "pre_sweep_role": {
            "role": "independent_lab_EMT_direction_filter_before_bounded_Aleph_sweep",
            "newly_supported_action": (
                "reject EMT-state candidates that require E-cadherin increase or "
                "Vimentin decrease during progressive TGF-beta EMT"
            ),
            "candidate_mechanism_families": [
                "cell_cell_adhesion", "intermediate_filament_remodelling",
                "growth_factor_conditioned_state", "EMT_MET_hysteresis",
            ],
            "required_Aleph_operators": [
                "E_cadherin_observation_bridge_or_declared_unavailable",
                "Vimentin_observation_bridge_or_declared_unavailable",
                "cell_cell_contact_and_dispersion_observables",
            ],
            "numeric_range": None,
            "may_emit_sweep_axis": False,
            "required_before_numeric_promotion": [
                "replicate-complete aligned external EMT experiment",
                "cross-modality uncertainty calibration and semantic OOD refusal",
                "Aleph observation operators with exact declared semantics",
                "monotonic forward sensitivity for candidate parameter axes",
            ],
        },
        "limitations": [
            "the released Source Data exposes one primary biological experiment and no replicate identifier",
            "90,066 cells improve representation but do not increase biological n above one",
            "CyTOF and IF share protein identity but not an absolute intensity scale",
            "CCAST state labels were constructed from the same six markers and cannot validate a classifier trained on them",
            "independent-lab direction transfer rejects incompatible states but cannot emit a numeric Aleph range",
        ],
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--exif-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor, args.manifest_report, args.exif_report)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "cross_lab_direction": report["cross_lab_common_marker_direction"][
            "independent_lab_common_marker_direction_reproduced"
        ],
        "inferential_holdout": report["independent_lab_inferential_holdout"],
        "OOD_refusal": report["OOD_refusal_validated"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
