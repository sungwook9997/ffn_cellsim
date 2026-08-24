#!/usr/bin/env python3
"""Descriptive temporal diagnostic for the single released M-TRACK trajectory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _rank(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2.0
        start = end
    return ranks


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx, ry = _rank(x), _rank(y)
    if rx.std() == 0 or ry.std() == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def evaluate(tensor_path: Path, manifest_report_path: Path,
             exif_report_path: Path | None = None) -> dict[str, object]:
    manifest = json.loads(manifest_report_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "aleph.outer_library.mtrack_emt_temporal_raw_vimentin.v1":
        raise ValueError("M-TRACK manifest schema changed")
    with np.load(tensor_path, allow_pickle=False) as data:
        X = data["X"].astype(np.float64)
        names = data["feature_names"].astype(str)
        time = data["time_index"].astype(np.float64)
        split_group = data["split_group"].astype(str)
    if X.shape != (10, 16) or len(set(split_group)) != 1:
        raise ValueError("M-TRACK one-trajectory tensor contract changed")
    indices = {name: int(np.flatnonzero(names == name)[0]) for name in names}
    diagnostic_names = (
        "mask_area_pixels", "mask_boundary_pixels", "mask_bbox_aspect_ratio",
        "VIM_RFP_masked_mean", "VIM_RFP_masked_std",
        "VIM_RFP_foreground_background_difference",
        "VIM_RFP_foreground_horizontal_gradient", "VIM_RFP_foreground_vertical_gradient",
    )
    trends = []
    for name in diagnostic_names:
        values = X[:, indices[name]]
        trends.append({
            "feature": name,
            "spearman_time_correlation": _spearman(time, values),
            "early_three_median": float(np.median(values[:3])),
            "late_three_median": float(np.median(values[-3:])),
            "late_minus_early": float(np.median(values[-3:]) - np.median(values[:3])),
        })
    vim_trend = next(row for row in trends if row["feature"] == "VIM_RFP_masked_mean")
    exif_context = None
    if exif_report_path is not None:
        exif = json.loads(exif_report_path.read_text(encoding="utf-8"))
        rows = exif["TGF_beta1_marker_direction_check"]["rows"]
        exif_vimentin = next(row for row in rows if row["marker"] == "Vimentin")
        exif_context = {
            "source": exif["source"]["dataset_id"],
            "fixed_IF_TGF_beta1_vs_control_vimentin_direction": exif_vimentin["observed_direction"],
            "MTRACK_live_VIM_RFP_late_vs_early_direction": (
                "increase" if vim_trend["late_minus_early"] > 0 else
                "decrease" if vim_trend["late_minus_early"] < 0 else "tie"
            ),
            "directions_agree": (
                exif_vimentin["observed_direction"] == "increase"
                and vim_trend["late_minus_early"] > 0
            ),
            "inferential_cross_assay_equivalence_claimed": False,
        }
    return {
        "schema": "aleph.outer_library.mtrack_emt_temporal_evaluation.v1",
        "task": "live_VIM_RFP_EMT_temporal_representation_diagnostic",
        "source": manifest["source"],
        "design_audit": {
            "timepoints": len(time), "trajectory_split_groups": len(set(split_group)),
            "biological_replicates_released": False, "untreated_control_released": False,
            "frames_are_repeated_technical_measurements": True,
            "frames_permitted_as_train_test_units": False,
            "paper_reports_10_minute_TRITC_interval": True,
        },
        "temporal_descriptive_trends": trends,
        "cross_lab_vimentin_direction_context": exif_context,
        "state_axis": "EMT_epithelial_mesenchymal_state",
        "representation_training_eligible": True,
        "observable_evidence_eligible": False,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "uncertainty_holdout_passed": False,
        "OOD_refusal_validated": False,
        "pre_sweep_role": {
            "role": "future_temporal_state_conditioning_and_sweep_space_reduction",
            "candidate_mechanism_families": [
                "intermediate_filament_remodelling", "cell_shape_remodelling",
                "cell_cell_adhesion", "actomyosin_contractility",
            ],
            "numeric_range": None, "may_emit_sweep_axis": False,
            "required_before_promotion": [
                "untreated and treated biological-replicate-complete trajectories",
                "independent-lab common-task temporal validation",
                "calibrated uncertainty and OOD refusal",
                "Aleph live-image observation operator and forward sensitivity",
            ],
        },
        "limitations": [
            "the public GUI example contains one xy01 trajectory and no untreated control",
            "ten timepoints are correlated technical measurements, not ten independent samples",
            "late-versus-early changes cannot be attributed uniquely to EMT without a control trajectory",
            "live endogenous VIM-RFP intensity is not assumed numerically equivalent to fixed-antibody IF",
            "the released image-data license is not stated",
            "no numeric Aleph parameter range is identified",
        ],
        "aleph_authority": "none", "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--exif-report", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor, args.manifest_report, args.exif_report)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "timepoints": report["design_audit"]["timepoints"],
        "vimentin_direction_context": report["cross_lab_vimentin_direction_context"],
        "may_emit_sweep_axis": report["pre_sweep_role"]["may_emit_sweep_axis"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
