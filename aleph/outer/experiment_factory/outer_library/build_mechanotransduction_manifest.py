#!/usr/bin/env python3
"""Normalize author-labelled HeLa mechanotransduction source workbooks."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

from schema import Observation


SOURCE_SHA256 = "c5fa3d068715e885ade13538f28cd998ca1bc463be4124282e187f6272c712fc"
ACCESSION = "10.5281/zenodo.18495219"
DATASET_ID = "zenodo-18495219-hela-mechanotransduction"
LAB_GROUP = "zenodo-18495219-authors"
LICENSE = "cc-by-4.0"


def _column_name(index: int) -> str:
    result = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(65 + remainder) + result
    return result


def _row(
    *, filename: str, sheet: str, row: int, column: int, sample: str,
    condition: str, cell_type: str, observable: str, value: float, unit: str,
    modality: str, coordinate_name: str = "", coordinate_value: float | None = None,
    coordinate_unit: str = "",
) -> Observation:
    cell = f"{_column_name(column)}{row + 1}"
    locator = f"zip:{filename}::{sheet}!{cell}"
    digest = hashlib.sha256(f"{SOURCE_SHA256}:{locator}".encode()).hexdigest()[:20]
    result = Observation(
        observation_id=f"mechanotransduction:{digest}",
        dataset_id=DATASET_ID,
        lab_group=LAB_GROUP,
        provider="Zenodo",
        accession=ACCESSION,
        license_id=LICENSE,
        source_sha256=SOURCE_SHA256,
        source_locator=locator,
        sample_id=f"{DATASET_ID}:{sample}",
        biological_replicate="unknown",
        technical_replicate=str(column),
        condition=condition,
        cell_type=cell_type,
        cell_state=condition,
        modality=modality,
        observable=observable,
        value=float(value),
        unit=unit,
        coordinate_name=coordinate_name,
        coordinate_value=coordinate_value,
        coordinate_unit=coordinate_unit,
    )
    result.validate()
    return result


def _workbook_by_prefix(payload: dict, prefix: str) -> tuple[str, dict]:
    matches = [w for w in payload["workbooks"] if os.path.basename(w["input"]).startswith(prefix)]
    if len(matches) != 1:
        raise ValueError(f"expected one workbook matching {prefix!r}, found {len(matches)}")
    return os.path.basename(matches[0]["input"]), matches[0]


def _calcium_cells(payload: dict) -> list[Observation]:
    specs = (
        ("[Fig.5d]", "untreated_control"),
        ("[Fig.5e]", "Piezo1_inhibited_GsMTx4"),
        ("[Fig.5f]", "TRPV4_inhibited_HC-067047"),
    )
    rows: list[Observation] = []
    for prefix, condition in specs:
        filename, workbook = _workbook_by_prefix(payload, prefix)
        sheet = workbook["sheets"][0]
        values = sheet["values"]
        for column in range(1, len(values[1])):
            sample = f"calcium:{condition}:cell-{column}"
            for row_index in range(1, len(values)):
                time = values[row_index][0]
                value = values[row_index][column]
                if not isinstance(time, (int, float)) or not isinstance(value, (int, float)):
                    raise ValueError(f"non-numeric calcium series in {filename}")
                rows.append(_row(
                    filename=filename,
                    sheet=sheet["name"],
                    row=row_index,
                    column=column,
                    sample=sample,
                    condition=condition,
                    cell_type="HeLa",
                    observable="intracellular_calcium_F_over_F0",
                    value=value,
                    unit="ratio",
                    modality="live_cell_calcium_imaging_derived",
                    coordinate_name="time",
                    coordinate_value=float(time),
                    coordinate_unit="s",
                ))
    return rows


def _fret_cells(payload: dict) -> list[Observation]:
    filename, workbook = _workbook_by_prefix(payload, "[Fig.4h]")
    sheet = workbook["sheets"][0]
    values = sheet["values"]
    rows: list[Observation] = []
    groups = ((1, 23, "HeLa"), (23, 45, "microvilli-deficient HeLa"))
    for start, stop, cell_type in groups:
        for column in range(start, stop):
            sample = f"fret:{cell_type}:cell-{column - start + 1}"
            for row_index in range(1, len(values)):
                time = values[row_index][0]
                value = values[row_index][column]
                if not isinstance(time, (int, float)) or not isinstance(value, (int, float)):
                    raise ValueError(f"non-numeric FRET series in {filename}")
                rows.append(_row(
                    filename=filename,
                    sheet=sheet["name"],
                    row=row_index,
                    column=column,
                    sample=sample,
                    condition=cell_type,
                    cell_type=cell_type,
                    observable="normalized_FRET_intensity",
                    value=value,
                    unit="ratio",
                    modality="live_cell_FRET_imaging_derived",
                    coordinate_name="time",
                    coordinate_value=float(time),
                    coordinate_unit="min",
                ))
    return rows


def _snapshot_table(
    payload: dict, prefix: str, observable: str, unit: str, modality: str,
) -> list[Observation]:
    filename, workbook = _workbook_by_prefix(payload, prefix)
    sheet = workbook["sheets"][0]
    values = sheet["values"]
    rows: list[Observation] = []
    for column in range(1, len(values[0])):
        condition = str(values[0][column]).strip()
        cell_type = "microvilli-deficient HeLa" if "deficient" in condition else "HeLa"
        for row_index in range(1, len(values)):
            value = values[row_index][column]
            if not isinstance(value, (int, float)):
                continue
            rows.append(_row(
                filename=filename,
                sheet=sheet["name"],
                row=row_index,
                column=column,
                sample=f"{observable}:{condition}:measurement-{row_index}",
                condition=condition,
                cell_type=cell_type,
                observable=observable,
                value=value,
                unit=unit,
                modality=modality,
            ))
    return rows


def build(workbook_export: Path) -> list[Observation]:
    payload = json.loads(workbook_export.read_text(encoding="utf-8"))
    rows = _calcium_cells(payload)
    rows.extend(_fret_cells(payload))
    rows.extend(_snapshot_table(
        payload, "[Fig.4c]", "bead_displacement", "um", "traction_force_microscopy_derived"
    ))
    rows.extend(_snapshot_table(
        payload, "[Fig.4d]", "traction_stress", "Pa", "traction_force_microscopy_derived"
    ))
    rows.extend(_snapshot_table(
        payload, "[Fig.4i]", "normalized_FRET_intensity", "ratio", "live_cell_FRET_imaging_derived"
    ))
    ids = [row.observation_id for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate observation IDs")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook-export", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(args.workbook_export)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({"observations": len(rows), "samples": len({r.sample_id for r in rows})}))


if __name__ == "__main__":
    main()

