#!/usr/bin/env python3
"""Build a leakage-safe SenSCOUT morphology/biomarker tensor from Dryad."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import io
import json
from pathlib import Path
import zipfile

import numpy as np


CSV_MEMBER = (
    "SenSCOUT-main/Data/"
    "Sampled_Dataset_Batch_Corrected_UMAPKMEANS_ZSCORE_logBM_transformed.csv"
)
EXPECTED_ARCHIVE_SHA256 = "1a4a70bbab2fc93d444bf062cc8df7619ea6b1c57344b411764be1fae88adb50"
EXPECTED_BIOMARKERS = {"BGAL", "BGAL_FLUOR", "HMGB1", "LMNB1", "P16"}
CONTINUOUS_MARKER_COLUMN = {
    "BGAL_FLUOR": "BGAL_FLUOR",
    "HMGB1": "HMGB1",
    "LMNB1": "LMNB1",
    "P16": "P16",
}


def _float_or_nan(value: str) -> float:
    if value in {"", "nan", "NaN"}:
        return float("nan")
    return float(value)


def build(
    archive: Path,
    archive_receipt: Path,
    feature_panel: Path,
    tensor_path: Path,
    report_path: Path,
) -> dict[str, object]:
    receipt = json.loads(archive_receipt.read_text(encoding="utf-8"))
    if (
        receipt.get("schema") != "aleph.outer_library.senscout_archive_receipt.v1"
        or receipt.get("sha256") != EXPECTED_ARCHIVE_SHA256
        or receipt.get("zip_crc_verified") is not True
    ):
        raise ValueError("SenSCOUT archive receipt is missing or unverified")
    panel = json.loads(feature_panel.read_text(encoding="utf-8"))
    features = panel["features"]
    if panel.get("feature_count") != 87 or len(features) != 87 or len(set(features)) != 87:
        raise ValueError("SenSCOUT curated morphology feature panel changed")

    values: list[list[float]] = []
    row_ids: list[str] = []
    repeats: list[str] = []
    image_groups: list[str] = []
    donor_codes: list[str] = []
    conditions: list[str] = []
    condition_labels: list[str] = []
    biomarkers: list[str] = []
    treated: list[int] = []
    bgal_labels: list[int] = []
    marker_values: list[float] = []
    p21_values: list[float] = []
    clusters: list[int] = []
    with zipfile.ZipFile(archive) as handle, handle.open(CSV_MEMBER) as raw:
        reader = csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8-sig", newline=""))
        if reader.fieldnames is None or len(reader.fieldnames) != 608:
            raise ValueError("SenSCOUT sampled table width changed")
        missing = sorted(set(features) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"SenSCOUT curated morphology features missing: {missing}")
        for row in reader:
            vector = [float(row[name]) for name in features]
            if not np.isfinite(vector).all():
                raise ValueError("SenSCOUT curated morphology contains non-finite values")
            row_id = row["Unnamed: 0"]
            repeat = row_id.split("_", 1)[0]
            condition = row["label"]
            donor, inducer = condition.split("_", 1)
            biomarker = row["Biomarker"]
            bgal = -1
            if biomarker == "BGAL":
                bgal = {"Non Senescent": 0, "Senescent": 1}.get(
                    row["BGAL"].strip(), -1
                )
            marker_column = CONTINUOUS_MARKER_COLUMN.get(biomarker)
            values.append(vector)
            row_ids.append(row_id)
            repeats.append(repeat)
            image_groups.append(row["FileName"])
            donor_codes.append(donor)
            conditions.append(condition)
            condition_labels.append(row["label_actual"])
            biomarkers.append(biomarker)
            treated.append(int(inducer != "CTRL"))
            bgal_labels.append(bgal)
            marker_values.append(
                _float_or_nan(row[marker_column]) if marker_column is not None else float("nan")
            )
            p21_values.append(_float_or_nan(row["P21"]))
            clusters.append(int(row["KMEANS_UPDATED"]))

    X = np.asarray(values, dtype=np.float32)
    if X.shape != (5_000, 87):
        raise ValueError(f"SenSCOUT sampled tensor shape changed: {X.shape}")
    repeat_counts = Counter(repeats)
    biomarker_counts = Counter(biomarkers)
    if repeat_counts != {"Bio1": 1_246, "Bio2": 1_838, "Bio3": 1_916}:
        raise ValueError(f"SenSCOUT biorepeat counts changed: {repeat_counts}")
    if set(biomarker_counts) != EXPECTED_BIOMARKERS or set(biomarker_counts.values()) != {1_000}:
        raise ValueError(f"SenSCOUT biomarker sampling changed: {biomarker_counts}")

    tensor_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        tensor_path,
        X=X,
        feature_name=np.asarray(features),
        row_id=np.asarray(row_ids),
        biorepeat=np.asarray(repeats),
        image_group=np.asarray(image_groups),
        donor_code=np.asarray(donor_codes),
        condition=np.asarray(conditions),
        condition_label=np.asarray(condition_labels),
        biomarker=np.asarray(biomarkers),
        noncontrol_condition=np.asarray(treated, dtype=np.int8),
        bgal_senescence_label=np.asarray(bgal_labels, dtype=np.int8),
        measured_biomarker_value=np.asarray(marker_values, dtype=np.float32),
        measured_p21_value=np.asarray(p21_values, dtype=np.float32),
        author_morphology_cluster=np.asarray(clusters, dtype=np.int8),
    )
    report = {
        "schema": "aleph.outer_library.senscout_morphology_tensor.v1",
        "source_doi": "10.5061/dryad.rbnzs7hp8",
        "archive_sha256": EXPECTED_ARCHIVE_SHA256,
        "rows": int(X.shape[0]),
        "morphology_feature_count": int(X.shape[1]),
        "feature_panel": panel,
        "biorepeat_counts": dict(sorted(repeat_counts.items())),
        "image_group_count": len(set(image_groups)),
        "donor_code_counts": dict(sorted(Counter(donor_codes).items())),
        "condition_counts": dict(sorted(Counter(conditions).items())),
        "biomarker_counts": dict(sorted(biomarker_counts.items())),
        "direct_bgal_labelled_rows": int(sum(value >= 0 for value in bgal_labels)),
        "bgal_unsure_or_missing_rows": int(sum(value < 0 for value in bgal_labels)),
        "continuous_primary_biomarker_rows": int(np.isfinite(marker_values).sum()),
        "paired_p21_rows": int(np.isfinite(p21_values).sum()),
        "input_provenance": "author_processed_sampled_morphology_features",
        "raw_cell_images_in_archive": False,
        "split_hierarchy": ["biorepeat", "image_group", "cell"],
        "cells_are_not_independent_biological_replicates": True,
        "author_cell_random_split_is_not_used": True,
        "target_leakage_exclusions": panel["excluded_target_or_leakage_families"],
        "author_kmeans_role": "stored_for_descriptive_reproduction_only_not_validation_target",
        "cross_biomarker_rows_are_not_claimed_as_same_cell": True,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "eligible_role": "same_lab_biorepeat_grouped_morphology_and_biomarker_training",
        "may_select_aleph_parameter": False,
        "aleph_authority": "none",
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--archive-receipt", type=Path, required=True)
    parser.add_argument("--feature-panel", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(
        args.archive, args.archive_receipt, args.feature_panel, args.tensor, args.output
    )
    print(json.dumps({key: report[key] for key in (
        "rows", "morphology_feature_count", "biorepeat_counts", "image_group_count",
        "direct_bgal_labelled_rows", "continuous_primary_biomarker_rows"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
