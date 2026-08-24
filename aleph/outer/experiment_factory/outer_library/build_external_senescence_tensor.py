#!/usr/bin/env python3
"""Build a cross-dataset fibroblast senescence marker tensor."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np

from build_gse250041_senescence_tensor import RNA_MARKERS


GSE63577_HEADERS = (
    "BJ_Y1", "BJ_Y2", "BJ_Y3", "BJ_OLD_1", "BJ_OLD_2", "BJ_OLD_3",
    "IMR90_Y1", "IMR90_Y2", "IMR90_Y3", "IMR90_O1", "IMR90_O2", "IMR90_O3",
    "WI_38_Y1", "WI_38_Y2", "WI_38_Y3", "WI_38_O1", "WI_38_O2", "WI_38_O3",
    "HFF_PD16_1", "HFF_PD16_2", "HFF_PD16_3", "HFF_PD74_1", "HFF_PD74_2", "HFF_PD74_3",
    "MRC_5_PD32_1", "MRC_5_PD32_2", "MRC_5_PD32_3", "MRC_5_PD72_1", "MRC_5_PD72_2", "MRC_5_PD72_3",
)
GSE282425_CONTEXT = {
    "1": ("DS2D", 1, "2D"),
    "2": ("DS3D", 1, "3D"),
    "3": ("EP2D", 0, "2D"),
    "4": ("EP3D", 0, "3D"),
}
FROZEN_GENE_ALIASES = {"IL8": "CXCL8"}


def _selected_counts_csv(path: Path, sample_headers: tuple[str, ...]) -> np.ndarray:
    selected = np.zeros((len(sample_headers), len(RNA_MARKERS)), dtype=np.float64)
    library = np.zeros(len(sample_headers), dtype=np.float64)
    marker_index = {name: i for i, name in enumerate(RNA_MARKERS)}
    found: set[str] = set()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or tuple(reader.fieldnames[4:]) != sample_headers:
            raise ValueError("GSE63577 converted counts header changed")
        for row in reader:
            values = np.asarray([float(row[name]) for name in sample_headers])
            library += values
            symbol = FROZEN_GENE_ALIASES.get(row["external_gene_id"], row["external_gene_id"])
            target = marker_index.get(symbol)
            if target is not None:
                selected[:, target] += values
                found.add(symbol)
    missing = sorted(set(RNA_MARKERS) - found)
    if missing:
        raise ValueError(f"GSE63577 marker genes missing: {missing}")
    if np.any(library <= 0):
        raise ValueError("GSE63577 contains an empty sample")
    return np.log1p(selected * (1_000_000.0 / library[:, None]))


def _gse63577(source_dir: Path, soffice: Path) -> tuple[np.ndarray, list[dict[str, object]]]:
    archive = source_dir / "GSE63577_counts_rpkm_exvivo_jenage_data.xls.gz"
    with tempfile.TemporaryDirectory(prefix="aleph-gse63577-") as temporary_text:
        temporary = Path(temporary_text)
        xls = temporary / "GSE63577_counts.xls"
        with gzip.open(archive, "rb") as source, xls.open("wb") as output:
            shutil.copyfileobj(source, output)
        subprocess.run(
            [str(soffice), "--headless", "--convert-to", "csv", "--outdir", str(temporary), str(xls)],
            check=True, capture_output=True, text=True,
        )
        csv_path = temporary / "GSE63577_counts.csv"
        if not csv_path.exists():
            raise ValueError("LibreOffice did not produce GSE63577 counts CSV")
        X = _selected_counts_csv(csv_path, GSE63577_HEADERS)
    rows = []
    for name in GSE63577_HEADERS:
        if name.startswith("BJ_"):
            cell_line, young = "BJ", "_Y" in name
        elif name.startswith("IMR90_"):
            cell_line, young = "IMR90", "_Y" in name
        elif name.startswith("WI_38_"):
            cell_line, young = "WI38", "_Y" in name
        elif name.startswith("HFF_"):
            cell_line, young = "HFF", "PD16" in name
        else:
            cell_line, young = "MRC5", "PD32" in name
        split = "fit" if cell_line in {"BJ", "HFF", "MRC5"} else "selection" if cell_line == "IMR90" else "calibration"
        rows.append({
            "dataset": "GSE63577", "lab": "Leibniz_Institute_JenAge",
            "sample_id": name, "biological_group": cell_line,
            "replicate": int(name[-1]), "state_label": int(not young),
            "context": "replicative_senescence_endpoint", "geometry": "2D",
            "role": split, "independent_biological_replicate": True,
        })
    return X, rows


def _gse297406(source_dir: Path) -> tuple[np.ndarray, list[dict[str, object]]]:
    path = source_dir / "GSE297406_hDF_raw_counts_p5_to_p15.csv.gz"
    headers = tuple([f"p5_{i}" for i in range(1, 5)] + [f"p15_{i}" for i in range(1, 5)])
    selected = np.zeros((8, len(RNA_MARKERS)), dtype=np.float64)
    library = np.zeros(8, dtype=np.float64)
    marker_index = {name: i for i, name in enumerate(RNA_MARKERS)}
    found: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or tuple(reader.fieldnames[-8:]) != headers:
            raise ValueError("GSE297406 count header changed")
        for row in reader:
            values = np.asarray([float(row[name]) for name in headers])
            library += values
            target = marker_index.get(row["Gene ID"])
            if target is not None:
                selected[:, target] += values
                found.add(row["Gene ID"])
    missing = sorted(set(RNA_MARKERS) - found)
    if missing:
        raise ValueError(f"GSE297406 marker genes missing: {missing}")
    X = np.log1p(selected * (1_000_000.0 / library[:, None]))
    rows = [{
        "dataset": "GSE297406", "lab": "Sungkyunkwan_University_Cholab",
        "sample_id": name, "biological_group": "human_dermal_fibroblast",
        "replicate": int(name.rsplit("_", 1)[1]), "state_label": int(name.startswith("p15")),
        "context": "passage_15_versus_passage_5", "geometry": "2D",
        "role": "sealed_external_lab_test", "independent_biological_replicate": True,
    } for name in headers]
    return X, rows


def _gse282425(source_dir: Path) -> tuple[np.ndarray, list[dict[str, object]]]:
    with gzip.open(source_dir / "GSE282425_features.tsv.gz", "rt", encoding="utf-8") as handle:
        features = [line.rstrip("\n").split("\t")[1] for line in handle]
    with gzip.open(source_dir / "GSE282425_barcodes.tsv.gz", "rt", encoding="utf-8") as handle:
        barcodes = [line.strip() for line in handle]
    marker_index = {name: i for i, name in enumerate(RNA_MARKERS)}
    selected_rows = {i + 1: marker_index[name] for i, name in enumerate(features) if name in marker_index}
    missing = sorted(set(RNA_MARKERS) - set(features))
    if missing:
        raise ValueError(f"GSE282425 marker genes missing: {missing}")
    selected = np.zeros((len(barcodes), len(RNA_MARKERS)), dtype=np.float64)
    library = np.zeros(len(barcodes), dtype=np.float64)
    with gzip.open(source_dir / "GSE282425_matrix.mtx.gz", "rt", encoding="ascii") as handle:
        dimensions = False
        for line in handle:
            if line.startswith("%"):
                continue
            values = line.split()
            if not dimensions:
                if tuple(map(int, values[:2])) != (len(features), len(barcodes)):
                    raise ValueError("GSE282425 matrix dimensions changed")
                dimensions = True
                continue
            gene, cell, count = map(int, values)
            library[cell - 1] += count
            target = selected_rows.get(gene)
            if target is not None:
                selected[cell - 1, target] += count
    cell_X = np.log1p(selected * (10_000.0 / library[:, None]))
    suffix = np.asarray([barcode.rsplit("-", 1)[1] for barcode in barcodes])
    X, rows = [], []
    for code, (name, label, geometry) in GSE282425_CONTEXT.items():
        mask = suffix == code
        X.append(cell_X[mask].mean(axis=0))
        rows.append({
            "dataset": "GSE282425", "lab": "Queen_Mary_University_London_Bishop_lab",
            "sample_id": name, "biological_group": "single_42_year_old_donor",
            "replicate": 1, "state_label": label,
            "context": "deeply_senescent_versus_early_proliferative",
            "geometry": geometry, "role": "single_donor_context_shift_diagnostic",
            "independent_biological_replicate": False, "cell_count": int(mask.sum()),
        })
    return np.asarray(X), rows


def build(source_dir: Path, soffice: Path, output: Path, report_path: Path) -> dict[str, object]:
    blocks = [_gse63577(source_dir, soffice), _gse297406(source_dir), _gse282425(source_dir)]
    X = np.concatenate([block[0] for block in blocks]).astype(np.float32)
    rows = [row for block in blocks for row in block[1]]
    if X.shape != (42, 35):
        raise ValueError(f"external senescence tensor shape changed: {X.shape}")
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "dataset", "lab", "sample_id", "biological_group", "replicate",
        "state_label", "context", "geometry", "role", "independent_biological_replicate",
    )
    np.savez_compressed(
        output, X=X, feature_name=np.asarray(RNA_MARKERS),
        **{field: np.asarray([row[field] for row in rows]) for field in fields},
    )
    report = {
        "schema": "aleph.outer_library.external_senescence_tensor.v1",
        "shape": list(X.shape), "feature_count": len(RNA_MARKERS),
        "feature_panel_frozen_before_external_test": True,
        "feature_panel_source": "GSE250041_frozen_RNA_MARKERS",
        "gene_aliases": FROZEN_GENE_ALIASES,
        "dataset_rows": {name: sum(row["dataset"] == name for row in rows) for name in ("GSE63577", "GSE297406", "GSE282425")},
        "lab_count": 3,
        "split_contract": {
            "fit_cell_lines": ["BJ", "HFF", "MRC5"],
            "selection_cell_line": "IMR90",
            "calibration_cell_line": "WI38",
            "sealed_external_dataset_lab": "GSE297406/Sungkyunkwan_University_Cholab",
            "external_labels_used_for_model_selection": False,
            "single_donor_context_shift_diagnostic": "GSE282425",
        },
        "gse63577_source_cells": {"sheet": "counts", "header": "A1:AH1", "gene_symbol_column": "B"},
        "gse63577_legacy_xls_visual_verification": "OLE_converted_to_XLSX_and_CSV_then_header_rendered",
        "cells_are_not_counted_as_biological_replicates": True,
        "aleph_authority": "none", "may_select_aleph_parameter": False,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--soffice", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.source_dir, args.soffice, args.output, args.report), sort_keys=True))


if __name__ == "__main__":
    main()
