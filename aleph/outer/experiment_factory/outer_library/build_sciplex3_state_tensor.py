#!/usr/bin/env python3
"""Stream sci-Plex3 sparse counts into context-grouped Reactome features."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import BinaryIO, Iterator

import numpy as np


CONTEXT_FIELDS = ("cell_type", "replicate", "time_point", "treatment", "dose")
CONTEXT_METADATA_FIELDS = (
    "vehicle",
    "catalog_number",
    "pathway_level_1",
    "pathway_level_2",
    "product_name",
    "target",
    "pathway",
)


def _triplet_chunks(handle: BinaryIO, block_bytes: int = 64 * 1024 * 1024) -> Iterator[np.ndarray]:
    remainder = b""
    while block := handle.read(block_bytes):
        payload = remainder + block
        boundary = payload.rfind(b"\n")
        if boundary < 0:
            remainder = payload
            continue
        complete, remainder = payload[:boundary], payload[boundary + 1:]
        values = np.fromstring(complete.decode("ascii"), dtype=np.int64, sep=" ")
        if values.size % 3:
            raise ValueError("sci-Plex3 sparse matrix row width changed")
        yield values.reshape(-1, 3)
    if remainder.strip():
        values = np.fromstring(remainder.decode("ascii"), dtype=np.int64, sep=" ")
        if values.size % 3:
            raise ValueError("sci-Plex3 final sparse row is incomplete")
        yield values.reshape(-1, 3)


def _contexts(
    pdata_path: Path,
) -> tuple[
    np.ndarray,
    list[tuple[str, ...]],
    list[tuple[str, ...]],
    np.ndarray,
    np.ndarray,
]:
    cell_context: list[int] = []
    context_index: dict[tuple[str, ...], int] = {}
    contexts: list[tuple[str, ...]] = []
    context_metadata: list[tuple[str, ...]] = []
    total_umi: list[float] = []
    cell_count: list[int] = []
    with gzip.open(pdata_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=" ", skipinitialspace=True)
        header = ["row_name"] + next(reader)
        for values in reader:
            row = dict(zip(header, values, strict=True))
            context = tuple(row[field] for field in CONTEXT_FIELDS)
            if any(value in {"", "NA"} for value in context):
                cell_context.append(-1)
                continue
            index = context_index.get(context)
            metadata = tuple(row[field] for field in CONTEXT_METADATA_FIELDS)
            if index is None:
                index = len(contexts)
                context_index[context] = index
                contexts.append(context)
                context_metadata.append(metadata)
                total_umi.append(0.0)
                cell_count.append(0)
            elif context_metadata[index] != metadata:
                raise ValueError(f"sci-Plex3 metadata changes within context {context!r}")
            cell_context.append(index)
            total_umi[index] += float(row["n.umi"])
            cell_count[index] += 1
    if len(cell_context) != 799_317 or len(contexts) != 4_896:
        raise ValueError(f"sci-Plex3 context cardinality changed: {len(cell_context)}, {len(contexts)}")
    return (
        np.asarray(cell_context, dtype=np.int32), contexts,
        context_metadata,
        np.asarray(total_umi, dtype=np.float64), np.asarray(cell_count, dtype=np.int32),
    )


def _treatment_splits(
    contexts: list[tuple[str, ...]],
    context_metadata: list[tuple[str, ...]],
) -> tuple[dict[str, int], list[str]]:
    """Create deterministic, pathway-stratified, treatment-disjoint splits.

    A pathway needs at least four distinct compounds to place one compound in
    every supervised split while retaining at least one training compound.
    Rarer pathway labels remain train-only and are not scored as external
    treatment holdouts. Vehicle contexts are a separate OOD reference split.
    """
    treatment_label: dict[str, str] = {}
    treatment_is_vehicle: dict[str, bool] = {}
    for context, metadata in zip(contexts, context_metadata, strict=True):
        treatment = context[3]
        label = metadata[2]
        is_vehicle = metadata[0] == "TRUE"
        if treatment in treatment_label and treatment_label[treatment] != label:
            raise ValueError(f"sci-Plex3 pathway label changes for {treatment}")
        if treatment in treatment_is_vehicle and treatment_is_vehicle[treatment] != is_vehicle:
            raise ValueError(f"sci-Plex3 vehicle flag changes for {treatment}")
        treatment_label[treatment] = label
        treatment_is_vehicle[treatment] = is_vehicle

    by_label: dict[str, list[str]] = defaultdict(list)
    split_by_treatment: dict[str, int] = {}
    for treatment, label in treatment_label.items():
        if treatment_is_vehicle[treatment]:
            split_by_treatment[treatment] = 4
        else:
            by_label[label].append(treatment)

    evaluable_labels: list[str] = []
    for label, treatments in sorted(by_label.items()):
        ordered = sorted(
            treatments,
            key=lambda treatment: hashlib.sha256(
                f"sciplex3-pathway-split-v1:{label}:{treatment}".encode()
            ).digest(),
        )
        if len(ordered) < 4:
            split_by_treatment.update((treatment, 0) for treatment in ordered)
            continue
        evaluable_labels.append(label)
        heldout_per_split = max(1, int(len(ordered) * 0.1 + 0.5))
        heldout_per_split = min(heldout_per_split, (len(ordered) - 1) // 3)
        boundaries = (heldout_per_split, 2 * heldout_per_split, 3 * heldout_per_split)
        for index, treatment in enumerate(ordered):
            if index < boundaries[0]:
                split = 1
            elif index < boundaries[1]:
                split = 2
            elif index < boundaries[2]:
                split = 3
            else:
                split = 0
            split_by_treatment[treatment] = split
    return split_by_treatment, evaluable_labels


def build(
    matrix_path: Path,
    pdata_path: Path,
    panel_path: Path,
    output_path: Path,
) -> dict[str, object]:
    panel = json.loads(panel_path.read_text(encoding="utf-8"))
    symbols = panel["panel_gene_symbols"]
    gene_count = int(panel["sciplex_gene_annotation_count"])
    row_to_panel = np.full(gene_count + 1, -1, dtype=np.int32)
    for panel_index, symbol in enumerate(symbols):
        for row_index in panel["sciplex_matrix_rows_1_based"][symbol]:
            row_to_panel[int(row_index)] = panel_index
    (
        cell_context,
        contexts,
        context_metadata,
        context_total_umi,
        context_cell_count,
    ) = _contexts(pdata_path)
    raw_counts = np.zeros((len(contexts), len(symbols)), dtype=np.uint64)
    triplet_count = 0
    matrix_total_umi = 0
    selected_triplet_count = 0
    maximum_gene_index = 0
    maximum_cell_index = 0
    with gzip.open(matrix_path, "rb") as handle:
        for triplets in _triplet_chunks(handle):
            gene = triplets[:, 0]
            cell = triplets[:, 1]
            value = triplets[:, 2]
            if (gene <= 0).any() or (cell <= 0).any() or (value <= 0).any():
                raise ValueError("sci-Plex3 sparse indices and counts must be positive")
            if gene.max(initial=0) > gene_count or cell.max(initial=0) > cell_context.size:
                raise ValueError("sci-Plex3 sparse index exceeds annotation dimensions")
            triplet_count += triplets.shape[0]
            matrix_total_umi += int(value.sum(dtype=np.int64))
            maximum_gene_index = max(maximum_gene_index, int(gene.max(initial=0)))
            maximum_cell_index = max(maximum_cell_index, int(cell.max(initial=0)))
            panel_index = row_to_panel[gene]
            context_index = cell_context[cell - 1]
            keep = (panel_index >= 0) & (context_index >= 0)
            if keep.any():
                np.add.at(raw_counts, (context_index[keep], panel_index[keep]), value[keep])
                selected_triplet_count += int(keep.sum())
    assigned_pdata_umi = float(context_total_umi.sum())
    selected_cpm = raw_counts / context_total_umi[:, None] * 1_000_000.0
    log_cpm = np.log1p(selected_cpm).astype(np.float32)
    context_array = np.asarray(contexts, dtype="U96")
    context_metadata_array = np.asarray(context_metadata, dtype="U128")
    split_by_treatment, evaluable_labels = _treatment_splits(contexts, context_metadata)
    treatment_split = np.asarray(
        [split_by_treatment[row[3]] for row in contexts], dtype=np.uint8
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        log1p_cpm=log_cpm,
        raw_panel_counts=raw_counts.astype(np.uint32),
        context=context_array,
        context_field=np.asarray(CONTEXT_FIELDS, dtype="U16"),
        context_metadata=context_metadata_array,
        context_metadata_field=np.asarray(CONTEXT_METADATA_FIELDS, dtype="U24"),
        context_total_umi=context_total_umi,
        context_cell_count=context_cell_count,
        gene_symbol=np.asarray(symbols, dtype="U32"),
        treatment_split=treatment_split,
        treatment_split_name=np.asarray(
            ["train_drugs", "selection_drugs", "calibration_drugs", "test_drugs", "vehicle_reference"],
            dtype="U20",
        ),
    )
    unique_treatments = sorted(
        {
            context[3]
            for context, metadata in zip(contexts, context_metadata, strict=True)
            if metadata[0] != "TRUE"
        }
    )
    split_treatments = {
        name: sum(split_by_treatment[treatment] == split_id for treatment in unique_treatments)
        for split_id, name in enumerate(("train", "selection", "calibration", "test"))
    }
    return {
        "schema": "aleph.outer_library.sciplex3_state_tensor.v1",
        "cell_count_in_matrix": int(cell_context.size),
        "assigned_context_cell_count": int((cell_context >= 0).sum()),
        "excluded_unassigned_cell_count": int((cell_context < 0).sum()),
        "context_count": len(contexts),
        "panel_gene_count": len(symbols),
        "sparse_triplet_count": triplet_count,
        "selected_panel_triplet_count": selected_triplet_count,
        "matrix_total_umi_all_cells": matrix_total_umi,
        "pdata_total_umi_assigned_cells": assigned_pdata_umi,
        "maximum_gene_index": maximum_gene_index,
        "maximum_cell_index": maximum_cell_index,
        "annotation_dimensions_bound_all_sparse_indices": True,
        "tensor_shape": list(log_cpm.shape),
        "treatment_disjoint_split_counts": split_treatments,
        "pathway_level_1_evaluable_on_all_holdouts": evaluable_labels,
        "pathway_level_1_evaluable_count": len(evaluable_labels),
        "rare_pathway_labels_are_train_only": True,
        "vehicle_reference_shared_as_supervised_target": False,
        "replicate_holdout_is_same_lab": True,
        "independent_lab_holdout": False,
        "state_label_authority": "drug_pathway_context_not_direct_cell_state_truth",
        "panel_provenance_schema": panel["schema"],
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--pdata", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.matrix, args.pdata, args.panel, args.tensor)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
