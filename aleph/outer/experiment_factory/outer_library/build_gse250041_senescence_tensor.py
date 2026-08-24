#!/usr/bin/env python3
"""Build a leakage-labelled RNA/ADT senescence representation from GSE250041."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np


RNA_MARKERS = (
    "CDKN1A", "CDKN2A", "GLB1", "LMNB1", "HMGB1", "SERPINE1",
    "IL6", "CXCL8", "CXCL1", "CXCL2", "MMP3", "MMP10", "TGFB1",
    "TP53", "RB1", "GADD45A", "DDB2", "BAX", "BBC3",
    "MKI67", "PCNA", "TOP2A", "MCM2", "MCM5", "MCM6", "CCNA2", "CCNB1",
    "COL1A1", "COL1A2", "COL3A1", "VIM", "DCN", "LUM", "FN1", "ACTA2",
)
CONDITIONS = (
    ("Proliferating", "untreated_proliferating", 0),
    ("Senescent", "IR_10_Gy_day_10_senescent", 1),
)


def _features(path: Path) -> list[tuple[str, str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return [tuple(row) for row in csv.reader(handle, delimiter="\t")]


def _read_selected(
    matrix_path: Path,
    features: list[tuple[str, str, str]],
    barcode_count: int,
) -> tuple[np.ndarray, list[str], list[str]]:
    selected_rows: dict[int, int] = {}
    names: list[str] = []
    origins: list[str] = []
    for row_index, (_, name, kind) in enumerate(features, start=1):
        if (kind == "Gene Expression" and name in RNA_MARKERS) or kind == "Antibody Capture":
            selected_rows[row_index] = len(names)
            names.append(name)
            origins.append("RNA" if kind == "Gene Expression" else "ADT")
    missing = sorted(set(RNA_MARKERS) - {name for name, origin in zip(names, origins) if origin == "RNA"})
    if missing:
        raise ValueError(f"GSE250041 marker genes missing from feature axis: {missing}")

    raw = np.zeros((barcode_count, len(names)), dtype=np.float32)
    gene_library = np.zeros(barcode_count, dtype=np.float64)
    with gzip.open(matrix_path, "rt", encoding="ascii") as handle:
        dimensions_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            values = line.split()
            if not dimensions_seen:
                rows, columns, _ = (int(value) for value in values)
                if rows != len(features) or columns != barcode_count:
                    raise ValueError(f"matrix/annotation dimensions disagree: {(rows, columns)}")
                dimensions_seen = True
                continue
            feature_index, barcode_index, count = (int(value) for value in values)
            if features[feature_index - 1][2] == "Gene Expression":
                gene_library[barcode_index - 1] += count
            selected = selected_rows.get(feature_index)
            if selected is not None:
                raw[barcode_index - 1, selected] = count
    if not dimensions_seen or np.any(gene_library <= 0):
        raise ValueError("GSE250041 matrix has missing dimensions or empty RNA libraries")

    values = np.empty_like(raw)
    rna_columns = np.asarray([origin == "RNA" for origin in origins])
    adt_columns = ~rna_columns
    values[:, rna_columns] = np.log1p(
        raw[:, rna_columns] * (10_000.0 / gene_library[:, None])
    )
    adt_log = np.log1p(raw[:, adt_columns])
    values[:, adt_columns] = adt_log - adt_log.mean(axis=1, keepdims=True)
    return values, names, origins


def build(source_dir: Path, output: Path, report_path: Path) -> dict[str, object]:
    features = _features(source_dir / "GSE250041_Proliferating_features.tsv.gz")
    all_values: list[np.ndarray] = []
    labels: list[np.ndarray] = []
    conditions: list[np.ndarray] = []
    feature_names: list[str] | None = None
    feature_origins: list[str] | None = None
    condition_rows: dict[str, int] = {}
    for prefix, condition, label in CONDITIONS:
        with gzip.open(source_dir / f"GSE250041_{prefix}_barcodes.tsv.gz", "rt") as handle:
            barcode_count = sum(1 for _ in handle)
        values, names, origins = _read_selected(
            source_dir / f"GSE250041_{prefix}_matrix.mtx.gz", features, barcode_count
        )
        if feature_names is not None and (names != feature_names or origins != feature_origins):
            raise ValueError("GSE250041 condition feature projections differ")
        feature_names, feature_origins = names, origins
        all_values.append(values)
        labels.append(np.full(barcode_count, label, dtype=np.int8))
        conditions.append(np.full(barcode_count, condition, dtype="U40"))
        condition_rows[condition] = barcode_count

    tensor = np.concatenate(all_values)
    label_array = np.concatenate(labels)
    condition_array = np.concatenate(conditions)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        X=tensor,
        author_condition_label=label_array,
        condition=condition_array,
        capture_id=condition_array,
        feature_name=np.asarray(feature_names),
        feature_origin=np.asarray(feature_origins),
    )
    report = {
        "schema": "aleph.outer_library.gse250041_senescence_tensor.v1",
        "rows": int(tensor.shape[0]),
        "features": int(tensor.shape[1]),
        "condition_rows": condition_rows,
        "rna_marker_count": int(sum(origin == "RNA" for origin in feature_origins or [])),
        "adt_feature_count": int(sum(origin == "ADT" for origin in feature_origins or [])),
        "normalization": {"RNA": "log1p_CPM_10000", "ADT": "cellwise_CLR_log1p"},
        "label_source": "author_condition_IR_10_Gy_day_10_versus_untreated",
        "split_group": "capture_id",
        "valid_dataset_split_count": 0,
        "valid_lab_split_count": 0,
        "classification_metric_authority": "none_capture_is_fully_confounded_with_state",
        "eligible_role": "multimodal_senescence_representation_training_only",
        "external_validation_eligible": False,
        "may_select_aleph_parameter": False,
        "aleph_authority": "none",
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.source_dir, args.output, args.report), sort_keys=True))


if __name__ == "__main__":
    main()
