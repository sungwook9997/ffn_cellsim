#!/usr/bin/env python3
"""Build the frozen donor/external-lab PBMC tensor without label leakage."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import tarfile
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.io import mmread


TRAIN_LABELS = {
    "B cells": "B_cell",
    "CD4 T cells": "CD4_T_cell",
    "CD8 T cells": "CD8_T_cell",
    "NK cells": "NK_cell",
    "CD14+ Monocytes": "monocyte",
    "FCGR3A+ Monocytes": "monocyte",
    "Dendritic cells": "dendritic_cell",
}
EXTERNAL_LABELS = {
    "B cell": "B_cell",
    "CD4+ T cell": "CD4_T_cell",
    "Cytotoxic T cell": "CD8_T_cell",
    "Natural killer cell": "NK_cell",
    "CD14+ monocyte": "monocyte",
    "CD16+ monocyte": "monocyte",
    "Dendritic cell": "dendritic_cell",
    "Plasmacytoid dendritic cell": "dendritic_cell",
}
SEMANTIC_OOD_LABELS = {"Megakaryocyte"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def verify_receipt(data_dir: Path, receipt_path: Path) -> dict[str, Any]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    for record in receipt["official_files"]:
        path = data_dir / record["name"]
        if path.stat().st_size != record["size_bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"acquired source differs from receipt: {path.name}")
    return receipt


def read_lines_gzip(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n\r") for line in handle]


def read_tar_gzip_lines(archive: tarfile.TarFile, name: str) -> list[str]:
    member = archive.extractfile(name)
    if member is None:
        raise ValueError(f"missing tar member {name}")
    with gzip.GzipFile(fileobj=member) as handle:
        return [line.decode("utf-8").rstrip("\n\r") for line in handle]


def read_tar_matrix(archive: tarfile.TarFile, name: str) -> sparse.csc_matrix:
    member = archive.extractfile(name)
    if member is None:
        raise ValueError(f"missing tar member {name}")
    with gzip.GzipFile(fileobj=member) as handle:
        return mmread(handle, spmatrix=True).tocsc().astype(np.float32)


def normalized_training_barcode(barcode: str, condition: str) -> str:
    if condition == "stim" and barcode.endswith("-11"):
        return barcode[:-3] + "-1"
    return barcode


def normalized_geo_cell(cell: str) -> str:
    experiment, method, barcode = cell.split(".", 2)
    method_map = {
        "Smart-seq2": "SM2",
        "CEL-Seq2": "Celseq2",
        "10x-Chromium-v2-A": "10x_v2_A",
        "10x-Chromium-v2-B": "10x_v2_B",
        "10x-Chromium-v2": "10x_v2",
        "10x-Chromium-v3": "10x_v3",
        "Drop-seq": "Drop",
        "inDrops": "inDrops",
        "Seq-Well": "Seqwell",
    }
    if method == "Seq-Well":
        run, barcode = barcode.split("-", 1)
        barcode = f"{barcode}_{run}"
    return f"{experiment.lower()}_{method_map[method]}_{barcode}".replace("-", "_")


def external_group(cell: str) -> tuple[str, str]:
    experiment, method, _ = cell.split(".", 2)
    return experiment, method


def log_normalize(matrix: sparse.csc_matrix) -> sparse.csc_matrix:
    library = np.asarray(matrix.sum(axis=0)).ravel()
    if np.any(library <= 0):
        raise ValueError("zero-library cell encountered")
    result = matrix.copy().astype(np.float32)
    scale = np.float32(10000.0) / library.astype(np.float32)
    result.data *= np.repeat(scale, np.diff(result.indptr))
    np.log1p(result.data, out=result.data)
    return result


def dense_selected(matrix: sparse.csc_matrix, rows: np.ndarray) -> np.ndarray:
    selected = matrix[rows, :].T.toarray().astype(np.float32, copy=False)
    library = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float32)
    selected *= (np.float32(10000.0) / library)[:, None]
    np.log1p(selected, out=selected)
    return selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    receipt = verify_receipt(args.data_dir, args.receipt)
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    genes = [line.split("\t")[1] for line in read_lines_gzip(args.data_dir / "GSE96583_batch2.genes.tsv.gz")]
    with tarfile.open(args.data_dir / "GSE96583_RAW.tar") as archive:
        control = read_tar_matrix(archive, "GSM2560248_2.1.mtx.gz")
        stimulated = read_tar_matrix(archive, "GSM2560249_2.2.mtx.gz")
        control_ids = read_tar_gzip_lines(archive, "GSM2560248_barcodes.tsv.gz")
        stimulated_ids = read_tar_gzip_lines(archive, "GSM2560249_barcodes.tsv.gz")
    if control.shape != (len(genes), len(control_ids)) or stimulated.shape != (len(genes), len(stimulated_ids)):
        raise ValueError("GSE96583 matrix axes do not match released genes/barcodes")

    metadata_path = args.data_dir / "GSE96583_batch2.total.tsne.df.tsv.gz"
    with gzip.open(metadata_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=[
            "cell_id", "tsne1", "tsne2", "ind", "stim", "cluster", "cell", "multiplets"
        ])
        next(reader)  # released header omits the row-name column
        metadata = {(row["stim"], normalized_training_barcode(row["cell_id"], row["stim"])): row for row in reader}

    matrices = []
    train_rows: list[dict[str, str]] = []
    for condition, matrix, barcodes in (("ctrl", control, control_ids), ("stim", stimulated, stimulated_ids)):
        keep = []
        for index, barcode in enumerate(barcodes):
            row = metadata.get((condition, barcode))
            if row is None:
                continue
            if row["multiplets"] != "singlet" or row["cell"] not in TRAIN_LABELS:
                continue
            keep.append(index)
            train_rows.append({
                "cell_id": f"{condition}:{barcode}", "donor": row["ind"],
                "condition": condition, "label": TRAIN_LABELS[row["cell"]],
                "author_label": row["cell"],
            })
        matrices.append(matrix[:, np.asarray(keep, dtype=np.int64)])
    train_counts = sparse.hstack(matrices, format="csc")
    if train_counts.shape[1] != len(train_rows):
        raise AssertionError("training metadata alignment failed")

    donors = sorted({row["donor"] for row in train_rows})
    frozen = protocol["gse96583_split"]
    released_matches_frozen = set(donors) == set(frozen["fit_donors"] + [frozen["calibration_donor"], frozen["sealed_internal_test_donor"]])
    if released_matches_frozen:
        fit_donors = frozen["fit_donors"]
        calibration_donor = frozen["calibration_donor"]
        test_donor = frozen["sealed_internal_test_donor"]
        split_rule = "frozen_identifiers"
    else:
        if len(donors) != 8:
            raise ValueError(f"fallback requires eight released donors, found {donors}")
        fit_donors, calibration_donor, test_donor = donors[:6], donors[6], donors[7]
        split_rule = "frozen_lexical_fallback"
    for row in train_rows:
        row["role"] = "fit" if row["donor"] in fit_donors else ("calibration" if row["donor"] == calibration_donor else "sealed_test")

    fit_mask = np.asarray([row["role"] == "fit" for row in train_rows])
    normalized_fit = log_normalize(train_counts[:, fit_mask])
    means = np.asarray(normalized_fit.mean(axis=1)).ravel()
    second = np.asarray(normalized_fit.power(2).mean(axis=1)).ravel()
    variances = np.maximum(second - means * means, 0.0)
    eligible = np.asarray([not (gene.startswith("MT-") or gene.startswith("RPL") or gene.startswith("RPS")) for gene in genes])
    ranked = sorted(np.flatnonzero(eligible), key=lambda index: (-float(variances[index]), genes[index]))[:512]
    ranked_symbols = [genes[index] for index in ranked]

    external_genes_raw = read_lines_gzip(args.data_dir / "GSE132044_pbmc_hg38_gene.tsv.gz")
    external_symbols = [line.split("_", 1)[1] if "_" in line else line for line in external_genes_raw]
    external_by_symbol: dict[str, int] = {}
    duplicate_external_symbols = Counter(external_symbols)
    for index, symbol in enumerate(external_symbols):
        external_by_symbol.setdefault(symbol, index)
    common_symbols = [symbol for symbol in ranked_symbols if symbol in external_by_symbol]
    if len(common_symbols) < protocol["feature_contract"]["minimum_common_features"]:
        raise ValueError(f"only {len(common_symbols)} frozen features occur in external data")
    train_by_symbol: dict[str, int] = {}
    for index, symbol in enumerate(genes):
        train_by_symbol.setdefault(symbol, index)
    train_feature_rows = np.asarray([train_by_symbol[symbol] for symbol in common_symbols])
    external_feature_rows = np.asarray([external_by_symbol[symbol] for symbol in common_symbols])
    train_x = dense_selected(train_counts, train_feature_rows)

    external_cells = read_lines_gzip(args.data_dir / "GSE132044_pbmc_hg38_cell.tsv.gz")
    geo_key_to_index: dict[str, int] = {}
    for index, cell in enumerate(external_cells):
        key = normalized_geo_cell(cell)
        if key in geo_key_to_index:
            raise ValueError(f"duplicate normalized GSE132044 cell key: {key}")
        geo_key_to_index[key] = index
    scp = json.loads((args.data_dir / "SCP424_Harmony_TSNE_CellType.json").read_text(encoding="utf-8"))["data"]
    selected_external: list[dict[str, Any]] = []
    unmatched_scp = 0
    excluded_external = Counter()
    for cell, author_label in zip(scp["cells"], scp["annotations"], strict=True):
        index = geo_key_to_index.get(cell.replace("-", "_"))
        if index is None:
            unmatched_scp += 1
            continue
        canonical = EXTERNAL_LABELS.get(author_label)
        is_ood = author_label in SEMANTIC_OOD_LABELS
        if canonical is None and not is_ood:
            excluded_external[author_label] += 1
            continue
        experiment, method = external_group(external_cells[index])
        selected_external.append({
            "matrix_index": index, "cell_id": external_cells[index],
            "author_label": author_label, "label": canonical or "",
            "is_semantic_ood": is_ood, "experiment": experiment, "method": method,
        })
    external_indices = np.asarray([row["matrix_index"] for row in selected_external], dtype=np.int64)
    with gzip.open(args.data_dir / "GSE132044_pbmc_hg38_count_matrix.mtx.gz", "rb") as handle:
        external_counts = mmread(handle, spmatrix=True).tocsc().astype(np.float32)
    if external_counts.shape != (len(external_symbols), len(external_cells)):
        raise ValueError("GSE132044 matrix axes do not match released genes/cells")
    external_x = dense_selected(external_counts[:, external_indices], external_feature_rows)

    zenodo_counts = {}
    for path in sorted((args.data_dir / "zenodo_3357167_labels").glob("*.csv")):
        with path.open("rt", encoding="utf-8") as handle:
            zenodo_counts[path.name] = max(sum(1 for _ in handle) - 1, 0)
    geo_method_counts = Counter(external_group(cell)[1] for cell in external_cells)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        train_X=train_x,
        train_cell_id=np.asarray([row["cell_id"] for row in train_rows]),
        train_label=np.asarray([row["label"] for row in train_rows]),
        train_donor=np.asarray([row["donor"] for row in train_rows]),
        train_condition=np.asarray([row["condition"] for row in train_rows]),
        train_role=np.asarray([row["role"] for row in train_rows]),
        external_X=external_x,
        external_cell_id=np.asarray([row["cell_id"] for row in selected_external]),
        external_label=np.asarray([row["label"] for row in selected_external]),
        external_author_label=np.asarray([row["author_label"] for row in selected_external]),
        external_experiment=np.asarray([row["experiment"] for row in selected_external]),
        external_method=np.asarray([row["method"] for row in selected_external]),
        external_is_semantic_ood=np.asarray([row["is_semantic_ood"] for row in selected_external]),
        feature_symbols=np.asarray(common_symbols),
    )
    report = {
        "schema": "aleph.outer_library.pbmc_cross_lab_tensor.v1",
        "protocol_sha256": sha256(args.protocol),
        "acquisition_receipt_sha256": sha256(args.receipt),
        "split": {"rule": split_rule, "fit_donors": fit_donors, "calibration_donor": calibration_donor, "sealed_test_donor": test_donor},
        "training": {
            "cells": len(train_rows), "genes": len(genes), "selected_fit_only_genes": 512,
            "common_frozen_features": len(common_symbols),
            "role_counts": dict(Counter(row["role"] for row in train_rows)),
            "label_counts": dict(Counter(row["label"] for row in train_rows)),
            "condition_counts": dict(Counter(row["condition"] for row in train_rows)),
        },
        "external": {
            "released_expression_cells": len(external_cells), "author_annotation_cells": len(scp["cells"]),
            "matched_and_retained_cells": len(selected_external), "unmatched_author_annotation_cells": unmatched_scp,
            "in_domain_cells": sum(not row["is_semantic_ood"] for row in selected_external),
            "semantic_ood_cells": sum(row["is_semantic_ood"] for row in selected_external),
            "label_counts": dict(Counter(row["author_label"] for row in selected_external)),
            "experiment_method_counts": {f"{key[0]}::{key[1]}": value for key, value in sorted(Counter((row["experiment"], row["method"]) for row in selected_external).items())},
            "excluded_author_labels": dict(excluded_external),
        },
        "secondary_benchmark_label_audit": {
            "decision": "rejected_as_classifier_targets",
            "reason": "selected CSV row counts do not align one-to-one with the official GSE132044 per-method expression cell axes",
            "zenodo_csv_rows": zenodo_counts,
            "official_geo_method_cells": dict(geo_method_counts),
        },
        "duplicate_external_gene_symbols_using_first_released_row": sum(count > 1 for count in duplicate_external_symbols.values()),
        "tensor_sha256": sha256(args.output),
        "authority": {"aleph_parameter_authority": "none", "role": "pre_sweep_context_compression_only"},
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"training_cells": len(train_rows), "external_cells": len(selected_external), "features": len(common_symbols), "tensor_sha256": report["tensor_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
