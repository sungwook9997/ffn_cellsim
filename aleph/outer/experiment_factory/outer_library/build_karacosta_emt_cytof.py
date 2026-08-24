#!/usr/bin/env python3
"""Build traceable HCC827 EMT CyTOF tensors without cell-level pseudoreplication."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np

from karacosta_xlsx import Workbook
from schema import Observation


DATASET_ID = "nature-2019-karacosta-hcc827-emt-cytof"
LAB_GROUP = "Plevritis_Bendall_Labs_Stanford"
ACCESSION = "10.1038/s41467-019-13441-6"
SOURCE_NAME = "41467_2019_13441_MOESM4_ESM.xlsx"
RAW_SHEET = "Complete CyTOF raw data set"
MARKERS = ("CD324", "Vimentin", "CD44", "CD24", "MUC1", "Twist")
MARKER_NAMES = {
    "CD324": "E_cadherin_CD324",
    "Vimentin": "Vimentin",
    "CD44": "CD44",
    "CD24": "CD24",
    "MUC1": "MUC1",
    "Twist": "Twist",
}
STATE_SHEETS = {
    "E1 state subset": ("E1", 2),
    "E2 state subset": ("E2", 1),
    "E3 state subset": ("E3", 4),
    "pEMT1 state subset": ("pEMT1", 3),
    "pEMT2 state subset raw": ("pEMT2", 6),
    "pEMT3 state subset": ("pEMT3", 7),
    "M state": ("M", 11),
    "MET state": ("MET", 13),
}
CLUSTER_TO_STATE = {cluster: state for state, cluster in STATE_SHEETS.values()}
TIMEPOINT = {
    1: ("0d_control", 0),
    2: ("2d_TGFB", 2),
    3: ("6d_TGFB", 6),
    4: ("10d_TGFB", 10),
    5: ("2d_withdrawal_after_10d_TGFB", 12),
    6: ("6d_withdrawal_after_10d_TGFB", 16),
    7: ("10d_withdrawal_after_10d_TGFB", 20),
}
EXPECTED_DIMENSIONS = {
    RAW_SHEET: (90_071, 54),
    # State sheets retain three formatted blank rows after their last event.
    "E1 state subset": (2_400, 55),
    "E2 state subset": (3_290, 55),
    "E3 state subset": (1_322, 55),
    "pEMT1 state subset": (3_365, 55),
    "pEMT2 state subset raw": (2_949, 55),
    "pEMT3 state subset": (2_669, 55),
    "M state": (9_615, 55),
    "MET state": (3_504, 55),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _float(row: dict[str, str | None], name: str) -> float:
    value = row.get(name)
    if value is None:
        raise ValueError(f"missing CyTOF value: {name}")
    return float(value)


def _figure_3c_means(book: Workbook) -> dict[str, dict[str, float]]:
    state_order = [state for state, _ in STATE_SHEETS.values()]
    label_to_marker = {
        "E-Cadherin": "CD324", "Vimentin": "Vimentin", "CD44": "CD44",
        "CD24": "CD24", "MUC1": "MUC1", "Twist": "Twist",
    }
    result: dict[str, dict[str, float]] = {}
    for _, cells in book.iter_rows("Figure 3c"):
        label = cells.get(0)
        if label not in label_to_marker:
            continue
        marker = label_to_marker[str(label)]
        if marker in result:
            continue
        values = [float(cells[index]) for index in range(1, 9)]
        result[marker] = dict(zip(state_order, values))
    if set(result) != set(MARKERS):
        raise ValueError("Figure 3c marker means changed")
    return result


def build(
    source_dir: Path, receipt_path: Path, manifest_path: Path, tensor_path: Path
) -> dict[str, object]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") != "aleph.outer_library.karacosta_emt_cytof_acquisition.v1":
        raise ValueError("Karacosta acquisition schema changed")
    source = source_dir / SOURCE_NAME
    pinned = receipt["files"][0]
    if (
        not source.is_file()
        or source.stat().st_size != pinned["size"]
        or _sha256(source) != pinned["sha256"]
    ):
        raise ValueError("Karacosta source workbook changed")
    source_sha = pinned["sha256"]
    book = Workbook(source)
    if set(EXPECTED_DIMENSIONS) - set(book.sheets):
        raise ValueError("Karacosta required sheets changed")
    for name, expected in EXPECTED_DIMENSIONS.items():
        observed = (book.sheets[name].max_row, book.sheets[name].max_column)
        if observed != expected:
            raise ValueError(f"Karacosta sheet dimension changed: {name} {observed}")

    raw_n = EXPECTED_DIMENSIONS[RAW_SHEET][0] - 5
    raw_X = np.empty((raw_n, len(MARKERS)), dtype=np.float64)
    raw_cluster = np.empty(raw_n, dtype=np.int8)
    raw_timepoint = np.empty(raw_n, dtype=np.int8)
    raw_elapsed_days = np.empty(raw_n, dtype=np.int8)
    raw_event = np.empty(raw_n, dtype=np.int32)
    raw_state = np.empty(raw_n, dtype="U24")
    raw_counts: Counter[tuple[int, int]] = Counter()
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as manifest:
        index = 0
        for source_row, row in book.records_after_header(
            RAW_SHEET, "Timepoint", unnamed_first_column="Event #"
        ):
            if index >= raw_n:
                raise ValueError("Karacosta raw sheet grew during parse")
            event = int(_float(row, "Event #"))
            timepoint = int(_float(row, "Timepoint"))
            cluster = int(_float(row, "clusterids"))
            if timepoint not in TIMEPOINT or not 1 <= cluster <= 13:
                raise ValueError("Karacosta timepoint/cluster domain changed")
            condition, elapsed = TIMEPOINT[timepoint]
            state = CLUSTER_TO_STATE.get(cluster, f"unretained_cluster_{cluster}")
            raw_X[index] = [_float(row, marker) for marker in MARKERS]
            raw_cluster[index] = cluster
            raw_timepoint[index] = timepoint
            raw_elapsed_days[index] = elapsed
            raw_event[index] = event
            raw_state[index] = state
            raw_counts[(timepoint, cluster)] += 1
            observation = Observation(
                observation_id=f"karacosta-raw-row{source_row}-cluster",
                dataset_id=DATASET_ID,
                lab_group=LAB_GROUP,
                provider="Stanford Plevritis/Bendall labs via Nature Communications Source Data",
                accession=ACCESSION,
                license_id="cc-by-4.0",
                source_sha256=source_sha,
                source_locator=f"{SOURCE_NAME}:{RAW_SHEET}!row{source_row}:clusterids",
                sample_id=f"karacosta_raw_event_{event}",
                biological_replicate=(
                    "Karacosta_HCC827_primary_experiment_biological_replicate_not_released"
                ),
                technical_replicate=f"single_cell_event_{event}",
                condition=condition,
                cell_type="HCC827_human_NSCLC_cell_line",
                cell_state=f"author_CCAST_EMT_state_{state}",
                modality="single_cell_mass_cytometry_CyTOF_nontransformed",
                observable="author_CCAST_cluster_id",
                value=float(cluster),
                unit="categorical_cluster_identifier",
                unit_source_sha256=source_sha,
                unit_source_locator=(
                    f"{SOURCE_NAME}:{RAW_SHEET}!rows1-5_cluster_and_timepoint_notes"
                ),
                coordinate_name="elapsed_experiment_time",
                coordinate_value=float(elapsed),
                coordinate_unit="day",
            )
            observation.validate()
            manifest.write(json.dumps(observation.as_dict(), sort_keys=True) + "\n")
            index += 1
    if index != raw_n or not np.array_equal(
        raw_event, np.arange(1, raw_n + 1, dtype=np.int32)
    ):
        raise ValueError("Karacosta raw event count/identity changed")

    state_rows: list[list[float]] = []
    state_labels: list[str] = []
    state_clusters: list[int] = []
    state_timepoints: list[int] = []
    state_elapsed: list[int] = []
    state_source_rows: list[int] = []
    state_sheet_names: list[str] = []
    state_counts: Counter[str] = Counter()
    for sheet, (state, expected_cluster) in STATE_SHEETS.items():
        for source_row, row in book.records_after_header(sheet, "Event #"):
            cluster = int(_float(row, "clusterids"))
            timepoint = int(_float(row, "Timepoint"))
            elapsed = int(_float(row, "daylabelid"))
            if cluster != expected_cluster or TIMEPOINT[timepoint][1] != elapsed:
                raise ValueError(f"Karacosta state subset mapping changed: {sheet}")
            state_rows.append([_float(row, marker) for marker in MARKERS])
            state_labels.append(state)
            state_clusters.append(cluster)
            state_timepoints.append(timepoint)
            state_elapsed.append(elapsed)
            state_source_rows.append(source_row)
            state_sheet_names.append(sheet)
            state_counts[state] += 1
    state_X = np.asarray(state_rows, dtype=np.float64)
    if state_X.shape != (29_066, 6):
        raise ValueError("Karacosta transformed state subset count changed")

    author_means = _figure_3c_means(book)
    reproduced_means: dict[str, dict[str, float]] = {}
    state_array = np.asarray(state_labels)
    for marker_index, marker in enumerate(MARKERS):
        reproduced_means[marker] = {}
        for state in STATE_SHEETS.values():
            state_name = state[0]
            mean = float(state_X[state_array == state_name, marker_index].mean())
            reproduced_means[marker][state_name] = mean
            if round(mean, 2) != round(author_means[marker][state_name], 2):
                raise ValueError(f"Figure 3c mean not reproduced: {marker}/{state_name}")

    np.savez_compressed(
        tensor_path,
        raw_X=raw_X,
        raw_cluster_id=raw_cluster,
        raw_timepoint_id=raw_timepoint,
        raw_elapsed_days=raw_elapsed_days,
        raw_event_id=raw_event,
        raw_state=raw_state,
        state_X=state_X,
        state_label=state_array,
        state_cluster_id=np.asarray(state_clusters, dtype=np.int8),
        state_timepoint_id=np.asarray(state_timepoints, dtype=np.int8),
        state_elapsed_days=np.asarray(state_elapsed, dtype=np.int8),
        state_source_row=np.asarray(state_source_rows, dtype=np.int32),
        state_source_sheet=np.asarray(state_sheet_names, dtype="U24"),
        feature_names=np.asarray(MARKERS, dtype="U12"),
        split_group=np.asarray([
            "Karacosta_HCC827_primary_experiment_biological_replicate_not_released"
        ] * raw_n, dtype="U80"),
    )
    return {
        "schema": "aleph.outer_library.karacosta_emt_cytof_manifest.v1",
        "source": {
            "dataset_id": DATASET_ID, "lab_group": LAB_GROUP,
            "paper_doi": ACCESSION, "source_sha256": source_sha,
        },
        "counts": {
            "raw_single_cell_events": raw_n,
            "raw_author_retained_state_events": int(np.isin(raw_cluster, list(CLUSTER_TO_STATE)).sum()),
            "raw_author_unretained_cluster_events": int((~np.isin(raw_cluster, list(CLUSTER_TO_STATE))).sum()),
            "transformed_downsampled_state_cells": len(state_X),
            "canonical_observations": raw_n,
            "canonical_samples": raw_n,
            "biological_replicate_groups": 1,
            "timepoints": 7,
        },
        "raw_cluster_by_timepoint_counts": {
            f"timepoint_{timepoint}_cluster_{cluster}": count
            for (timepoint, cluster), count in sorted(raw_counts.items())
        },
        "transformed_state_counts": dict(sorted(state_counts.items())),
        "feature_names": list(MARKERS),
        "author_figure_3c_means": author_means,
        "reproduced_transformed_state_means": reproduced_means,
        "replicate_contract": {
            "split_unit": "entire_released_primary_HCC827_CyTOF_experiment",
            "single_cells_may_be_split_as_biological_replicates": False,
            "timepoints_may_be_split_as_biological_replicates": False,
            "released_source_contains_biological_replicate_identifier": False,
            "paper_mentions_independent_replicate_in_Supplementary_Figure_6": True,
            "released_Source_Data_maps_that_independent_replicate": False,
        },
        "label_contract": {
            "cluster_labels_are_author_CCAST_outputs": True,
            "six_features_were_used_to_construct_the_cluster_labels": True,
            "state_classification_on_these_features_would_be_label_circular": True,
            "per_cell_state_classifier_holdout_permitted": False,
        },
        "manifest_sha256": _sha256(manifest_path),
        "tensor_sha256": _sha256(tensor_path),
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.source_dir, args.receipt, args.manifest, args.tensor)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
