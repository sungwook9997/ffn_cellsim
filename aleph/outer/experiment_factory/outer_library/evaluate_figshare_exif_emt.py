#!/usr/bin/env python3
"""Evaluate the ExIF raw-IF EMT pilot without granting state authority."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


EXPECTED_DIRECTION = {
    "CD44_variant_9": "decrease",
    "E_cadherin": "decrease",
    "N_cadherin": "increase",
    "EpCAM": "decrease",
    "Vimentin": "increase",
}


def _macro_f1(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    scores = []
    for label in (0, 1, 2):
        tp = int(np.sum((y_true == label) & (y_pred == label)))
        fp = int(np.sum((y_true != label) & (y_pred == label)))
        fn = int(np.sum((y_true == label) & (y_pred != label)))
        denominator = 2 * tp + fp + fn
        scores.append(0.0 if denominator == 0 else 2.0 * tp / denominator)
    return float(np.mean(scores))


def _ridge_predict(X_train: np.ndarray, y_train: np.ndarray,
                   X_test: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    mean = X_train.mean(axis=0)
    scale = X_train.std(axis=0)
    scale[scale < 1e-7] = 1.0
    train = (X_train - mean) / scale
    test = (X_test - mean) / scale
    train = np.column_stack((np.ones(len(train)), train))
    test = np.column_stack((np.ones(len(test)), test))
    target = np.eye(3, dtype=np.float64)[y_train]
    penalty = np.eye(train.shape[1], dtype=np.float64) * alpha
    penalty[0, 0] = 0.0
    weights = np.linalg.solve(train.T @ train + penalty, train.T @ target)
    return np.argmax(test @ weights, axis=1)


def _transfer(X: np.ndarray, label: np.ndarray, columns: np.ndarray,
              train_columns: tuple[str, ...], test_columns: tuple[str, ...]) -> dict[str, object]:
    train = np.isin(columns, train_columns)
    test = np.isin(columns, test_columns)
    expected_train = int(train.sum()) // 3
    expected_test = int(test.sum()) // 3
    if np.bincount(label[train], minlength=3).tolist() != [expected_train] * 3:
        raise ValueError("ExIF paired-well train balance changed")
    if np.bincount(label[test], minlength=3).tolist() != [expected_test] * 3:
        raise ValueError("ExIF paired-well test balance changed")
    prediction = _ridge_predict(X[train], label[train], X[test])
    return {
        "train_columns": list(train_columns),
        "test_columns": list(test_columns),
        "train_samples": int(train.sum()),
        "test_samples": int(test.sum()),
        "accuracy": float(np.mean(prediction == label[test])),
        "macro_f1": _macro_f1(label[test], prediction),
        "confusion_rows_true_columns_predicted": [
            [int(np.sum((label[test] == actual) & (prediction == predicted)))
             for predicted in (0, 1, 2)]
            for actual in (0, 1, 2)
        ],
    }


def evaluate(tensor_path: Path, manifest_report_path: Path) -> dict[str, object]:
    manifest_report = json.loads(manifest_report_path.read_text(encoding="utf-8"))
    if manifest_report.get("schema") != "aleph.outer_library.figshare_exif_emt_raw_if.v1":
        raise ValueError("ExIF manifest report schema changed")
    with np.load(tensor_path, allow_pickle=False) as data:
        X = data["X"].astype(np.float64)
        feature_names = data["feature_names"].astype(str)
        labels = data["label"].astype(np.int64)
        markers = data["marker"].astype(str)
        conditions = data["condition"].astype(str)
        sample_ids = data["sample_id"].astype(str)
        split_groups = data["split_group"].astype(str)
        dataset_id = str(data["dataset_id"][0])
        lab_group = str(data["lab_group"][0])
    n = manifest_report["counts"]["field_samples"]
    if X.shape != (n, 60) or len(sample_ids) != n:
        raise ValueError("ExIF feature tensor shape changed")
    if len(set(sample_ids)) != n or len(set(split_groups)) != 48:
        raise ValueError("ExIF sample/split grouping changed")
    if np.bincount(labels, minlength=3).tolist() != [n // 3] * 3:
        raise ValueError("ExIF condition balance changed")

    # General channels only: TD, DAPI, beta-catenin, and F-actin.  The held
    # variable-marker channel is deliberately excluded from this transfer test.
    general = np.asarray([
        not name.startswith("variable_EMT_marker_") for name in feature_names
    ])
    columns = np.asarray([sample.split(":")[1][1:] for sample in sample_ids])
    paired_transfers = [
        _transfer(X[:, general], labels, columns, ("01", "03", "05"), ("02", "04", "06")),
        _transfer(X[:, general], labels, columns, ("02", "04", "06"), ("01", "03", "05")),
    ]

    variable_mean_index = int(np.flatnonzero(
        feature_names == "variable_EMT_marker_raw_mean"
    )[0])
    direction_rows = []
    reproduced = 0
    for marker in sorted(set(markers)):
        mask = markers == marker
        values = X[mask, variable_mean_index]
        marker_conditions = conditions[mask]
        medians = {
            condition: float(np.median(values[marker_conditions == condition]))
            for condition in sorted(set(marker_conditions))
        }
        control = medians["untreated_control"]
        tgf = medians["TGF_beta1_10ng_per_mL_48h"]
        observed = "increase" if tgf > control else "decrease" if tgf < control else "tie"
        expected = EXPECTED_DIRECTION.get(marker)
        agrees = None if expected is None else observed == expected
        reproduced += int(agrees is True)
        direction_rows.append({
            "marker": marker,
            "raw_field_mean_intensity_medians": medians,
            "TGF_beta1_minus_control": tgf - control,
            "observed_direction": observed,
            "paper_expected_direction": expected,
            "direction_reproduced": agrees,
            "n_fields_per_condition": int(np.sum(marker_conditions == "untreated_control")),
        })

    return {
        "schema": "aleph.outer_library.figshare_exif_emt_evaluation.v1",
        "task": "raw_multichannel_IF_EMT_condition_representation",
        "source": {
            "dataset_id": dataset_id,
            "dataset_doi": "10.6084/m9.figshare.26500210.v1",
            "paper_doi": "10.1038/s41467-025-59592-7",
            "lab_group": lab_group,
        },
        "design_audit": {
            "field_samples": n,
            "well_split_groups": len(set(split_groups)),
            "biological_replicates_released": False,
            "plate_dates": 1,
            "fields_are_technical_measurements": True,
            "condition_is_confounded_with_plate_column_pair": True,
            "common_task_with_HPA_localization_labels": False,
        },
        "paired_well_column_transfer_general_channels_only": paired_transfers,
        "mean_paired_transfer_accuracy": float(np.mean([
            row["accuracy"] for row in paired_transfers
        ])),
        "mean_paired_transfer_macro_f1": float(np.mean([
            row["macro_f1"] for row in paired_transfers
        ])),
        "TGF_beta1_marker_direction_check": {
            "paper_expected_marker_count": len(EXPECTED_DIRECTION),
            "directions_reproduced": reproduced,
            "all_expected_directions_reproduced": reproduced == len(EXPECTED_DIRECTION),
            "rows": direction_rows,
            "inferential_test_permitted": False,
            "reason": "two wells and one admitted field per condition/marker are technical measurements from one plate date",
        },
        "state_axis": "EMT_epithelial_mesenchymal_state",
        "state_axis_status": "author_treatment_labelled_raw_IF_training_only",
        "representation_training_eligible": True,
        "observable_evidence_eligible": False,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "uncertainty_holdout_passed": False,
        "OOD_refusal_validated": False,
        "pre_sweep_role": {
            "role": "future_state_conditioning_and_sweep_space_reduction",
            "candidate_mechanism_families": [
                "cell_cell_adhesion", "intermediate_filament_remodelling",
                "actomyosin_contractility", "growth_factor_signalling_context",
            ],
            "numeric_range": None,
            "may_emit_sweep_axis": False,
            "required_before_promotion": [
                "independent-lab raw-IF EMT condition holdout",
                "biological-replicate-complete uncertainty calibration",
                "orthogonal migration or force observables aligned to the same state",
                "Aleph observation operators and forward sensitivity",
            ],
        },
        "limitations": [
            "the author labels are treatment conditions that weakly or strongly drive EMT, not per-cell EMT ground truth",
            "all admitted fields come from one released plate date",
            "paired-column transfer cannot remove plate-position confounding",
            "this source is an independent provider from HPA but does not share HPA localization labels",
            "no numeric Aleph parameter range is identified",
        ],
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor, args.manifest_report)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "mean_paired_transfer_accuracy": report["mean_paired_transfer_accuracy"],
        "directions_reproduced": report["TGF_beta1_marker_direction_check"]["directions_reproduced"],
        "state_axis_status": report["state_axis_status"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
