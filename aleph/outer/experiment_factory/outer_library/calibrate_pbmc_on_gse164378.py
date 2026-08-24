#!/usr/bin/env python3
"""Fit only the preregistered PBMC temperature and conformal thresholds."""

from __future__ import annotations

import argparse
import csv
import gc
import gzip
import hashlib
import json
import re
import tarfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import mmread

from train_pbmc_cross_lab_model import (
    CELL_CLASSES,
    confusion_metrics,
    conformal_metrics,
    conformal_thresholds,
    ece,
    select_temperature,
    softmax_from_logits,
)


MATRIX_MEMBER = "GSM5008737_RNA_3P-matrix.mtx.gz"
FEATURE_MEMBER = "GSM5008737_RNA_3P-features.tsv.gz"
BARCODE_MEMBER = "GSM5008737_RNA_3P-barcodes.tsv.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def verify_receipt(data_dir: Path, receipt_path: Path) -> dict[str, Any]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    for record in receipt["files"]:
        path = data_dir / record["name"]
        if path.stat().st_size != record["size_bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"source differs from acquisition receipt: {path.name}")
    return receipt


def tar_gzip_lines(archive: tarfile.TarFile, member_name: str) -> list[str]:
    member = archive.extractfile(member_name)
    if member is None:
        raise ValueError(f"missing tar member {member_name}")
    with gzip.GzipFile(fileobj=member) as handle:
        return [line.decode("utf-8").rstrip("\r\n") for line in handle]


def map_label(label: str, protocol: dict[str, Any]) -> str | None:
    lowered = label.casefold()
    if "t" in lowered and not any(token in lowered for token in ("cd4", "cd8", "helper", "regulatory", "treg", "cytotoxic")):
        return None
    for item in protocol["coarse_cell_type_mapping_in_order"]:
        if any(re.search(pattern, lowered, flags=re.IGNORECASE) for pattern in item["case_insensitive_patterns"]):
            return item["canonical"]
    return None


def encode(labels: list[str]) -> np.ndarray:
    index = {str(label): position for position, label in enumerate(CELL_CLASSES)}
    return np.asarray([index[label] for label in labels], dtype=np.int64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--source-model", type=Path, required=True)
    parser.add_argument("--output-model", type=Path, required=True)
    parser.add_argument("--output-report", type=Path, required=True)
    args = parser.parse_args()
    receipt = verify_receipt(args.data_dir, args.receipt)
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    source_model = np.load(args.source_model)
    feature_symbols = source_model["feature_symbols"].astype(str)
    parameters = source_model["cell_parameters"]
    original_parameter_hash = array_sha256(parameters)
    feature_count = len(feature_symbols)
    class_count = len(CELL_CLASSES)
    weights = parameters[: feature_count * (class_count - 1)].reshape(feature_count, class_count - 1)
    intercept = parameters[feature_count * (class_count - 1):]

    metadata_path = args.data_dir / "GSE164378_sc.meta.data_3P.csv.gz"
    with gzip.open(metadata_path, "rt", encoding="utf-8-sig", newline="") as handle:
        metadata_rows = list(csv.DictReader(handle))
    metadata_by_cell = {row[""]: row for row in metadata_rows}
    if len(metadata_by_cell) != len(metadata_rows):
        raise ValueError("duplicate GSE164378 metadata cell identifiers")

    tar_path = args.data_dir / "GSE164378_RAW.tar"
    with tarfile.open(tar_path) as archive:
        barcodes = tar_gzip_lines(archive, BARCODE_MEMBER)
        feature_rows = tar_gzip_lines(archive, FEATURE_MEMBER)
        source_symbols = [row.split("\t")[1] for row in feature_rows]
        source_by_symbol: dict[str, int] = {}
        duplicate_symbols = Counter(source_symbols)
        for row_index, symbol in enumerate(source_symbols):
            source_by_symbol.setdefault(symbol, row_index)
        common_model_positions = np.asarray([
            position for position, symbol in enumerate(feature_symbols)
            if symbol in source_by_symbol
        ], dtype=np.int64)
        common_source_rows = np.asarray([
            source_by_symbol[feature_symbols[position]] for position in common_model_positions
        ], dtype=np.int64)
        if len(common_source_rows) < protocol["frozen_input_policy"]["minimum_common_features"]:
            raise ValueError(f"only {len(common_source_rows)} model features occur in GSE164378")
        matrix_member = archive.extractfile(MATRIX_MEMBER)
        if matrix_member is None:
            raise ValueError(f"missing tar member {MATRIX_MEMBER}")
        with gzip.GzipFile(fileobj=matrix_member) as matrix_handle:
            counts = mmread(matrix_handle, spmatrix=True).tocsr()
    if counts.shape != (len(source_symbols), len(barcodes)):
        raise ValueError("GSE164378 count axes do not match features/barcodes")
    if len(metadata_by_cell) != len(barcodes) or any(cell not in metadata_by_cell for cell in barcodes):
        raise ValueError("GSE164378 metadata and RNA barcode axes do not align exactly")
    library_size = np.asarray(counts.sum(axis=0)).ravel().astype(np.float64)
    if np.any(library_size < protocol["frozen_input_policy"]["minimum_library_size"]):
        raise ValueError("GSE164378 contains a forbidden zero-library cell")
    selected = counts[common_source_rows, :].T.tocsr().astype(np.float64)
    del counts
    gc.collect()
    scale = 10000.0 / library_size
    selected.data *= np.repeat(scale, np.diff(selected.indptr))
    np.log1p(selected.data, out=selected.data)
    nonreference_logits = selected @ weights[common_model_positions, :] + intercept
    logits = np.column_stack([
        np.asarray(nonreference_logits), np.zeros(len(barcodes), dtype=np.float64)
    ])
    del selected, nonreference_logits
    gc.collect()

    retained_indices = []
    retained_labels = []
    retained_donors = []
    retained_author_labels = []
    excluded = Counter()
    for cell_index, cell in enumerate(barcodes):
        row = metadata_by_cell[cell]
        author_label = row["celltype.l1"]
        canonical = map_label(author_label, protocol)
        if canonical is None:
            excluded[author_label] += 1
            continue
        retained_indices.append(cell_index)
        retained_labels.append(canonical)
        retained_donors.append(row["donor"])
        retained_author_labels.append(author_label)
    retained_indices_array = np.asarray(retained_indices, dtype=np.int64)
    calibration_logits = logits[retained_indices_array]
    y = encode(retained_labels)
    class_counts = Counter(retained_labels)
    minimum = protocol["calibration_contract"]["minimum_calibration_cells_per_class"]
    if any(class_counts[str(label)] < minimum for label in CELL_CLASSES):
        raise ValueError(f"calibration class below frozen minimum: {class_counts}")
    temperature, temperature_fit = select_temperature(calibration_logits, y)
    calibrated_probability = softmax_from_logits(calibration_logits / temperature)
    thresholds = conformal_thresholds(y, calibrated_probability, class_count)
    metrics = confusion_metrics(y, calibrated_probability, CELL_CLASSES)
    metrics["ece_15_bin"] = ece(y, calibrated_probability)
    metrics["conformal"] = conformal_metrics(y, calibrated_probability, thresholds)

    donor_metrics = {}
    donor_array = np.asarray(retained_donors)
    for donor in sorted(set(retained_donors)):
        mask = donor_array == donor
        item = confusion_metrics(y[mask], calibrated_probability[mask], CELL_CLASSES)
        item["ece_15_bin"] = ece(y[mask], calibrated_probability[mask])
        item["conformal"] = conformal_metrics(y[mask], calibrated_probability[mask], thresholds)
        donor_metrics[donor] = item

    args.output_model.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_model,
        source_model_sha256=np.asarray(sha256(args.source_model)),
        source_parameter_sha256=np.asarray(original_parameter_hash),
        feature_symbols=feature_symbols,
        cell_classes=CELL_CLASSES,
        cell_parameters=parameters,
        cell_temperature=np.asarray(temperature),
        cell_conformal_thresholds=thresholds,
        calibration_accession=np.asarray("GSE164378"),
    )
    if array_sha256(np.load(args.output_model)["cell_parameters"]) != original_parameter_hash:
        raise AssertionError("calibration artifact changed classifier parameters")
    report = {
        "schema": "aleph.outer_library.pbmc_gse164378_calibration.v1",
        "protocol_sha256": sha256(args.protocol),
        "acquisition_receipt_sha256": sha256(args.receipt),
        "source_model_sha256": sha256(args.source_model),
        "source_parameter_sha256": original_parameter_hash,
        "output_parameter_sha256": array_sha256(np.load(args.output_model)["cell_parameters"]),
        "classifier_weights_bit_identical": True,
        "calibration": {
            "source": "GSE164378", "lab_group": "NYGC_Satija_Hao_multimodal_PBMC",
            "total_cells": len(barcodes), "retained_six_class_cells": len(retained_indices),
            "distinct_donors": sorted(set(retained_donors)),
            "class_counts": dict(class_counts), "author_label_counts": dict(Counter(retained_author_labels)),
            "excluded_author_labels": dict(excluded), "common_model_features": len(common_source_rows),
            "duplicate_source_symbols_using_first_row": sum(value > 1 for value in duplicate_symbols.values()),
            "temperature": temperature,
            "class_conditional_thresholds": {str(CELL_CLASSES[i]): float(value) for i, value in enumerate(thresholds)},
            "fit": temperature_fit, "same_calibration_lab_diagnostics": metrics,
            "same_calibration_lab_diagnostics_by_donor": donor_metrics,
        },
        "retired_confirmation_source": {
            "accession": "GSE222647", "status": "retired_before_model_execution",
            "reason": "frozen exact donor/condition column contract did not admit Sample.Name or Study.Classification",
            "expression_used_for_prediction": False, "labels_used_for_metrics_or_calibration": False,
            "may_be_reused_as_pristine": False,
        },
        "next_confirmation_source": "GSE158055",
        "output_model_sha256": sha256(args.output_model),
        "authority": {
            "role": "pre_sweep_context_calibration_only", "aleph_parameter_authority": "none",
            "may_delete_sweep_axis": False, "may_emit_numeric_sweep_range": False,
            "may_select_aleph_parameter": False,
        },
    }
    args.output_report.parent.mkdir(parents=True, exist_ok=True)
    args.output_report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "cells": len(retained_indices), "donors": len(set(retained_donors)),
        "features": len(common_source_rows), "temperature": temperature,
        "same_lab_coverage_diagnostic": metrics["conformal"]["marginal_coverage"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
