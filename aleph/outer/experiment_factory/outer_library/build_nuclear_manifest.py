#!/usr/bin/env python3
"""Normalize two unambiguous author-labelled tables from the nuclear workbook.

Workbook parsing itself is performed by ``export_workbook.mjs`` through the
artifact-tool spreadsheet runtime.  This module consumes that lossless JSON
export.  We intentionally start with tables whose row/column semantics are
explicit, instead of heuristically scraping every numeric cell in the book.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from schema import Observation


SOURCE_SHA256 = "9a5503ba31b22d12417fd70afde7af28c4127993ac6fcfcb2da577decbc13e3e"
ACCESSION = "10.6084/m9.figshare.27241950"
DATASET_ID = "figshare-27241950-nuclear-mechanics"
LAB_GROUP = "figshare-27241950-authors"
LICENSE = "cc-by-4.0"


def _column_name(index: int) -> str:
    result = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(65 + remainder) + result
    return result


def _sheet_map(workbook_export: Path) -> dict[str, list[list[object]]]:
    payload = json.loads(workbook_export.read_text(encoding="utf-8"))
    return {sheet["name"]: sheet["values"] for sheet in payload["sheets"]}


def _observation(
    *, sheet: str, row: int, column: int, sample_id: str, biological_replicate: str,
    technical_replicate: str, condition: str, observable: str, value: float,
    modality: str = "immunofluorescence_derived",
) -> Observation:
    locator = f"xlsx:{sheet}!{_column_name(column)}{row + 1}"
    digest = hashlib.sha256(f"{SOURCE_SHA256}:{locator}".encode()).hexdigest()[:20]
    result = Observation(
        observation_id=f"nuclear:{digest}",
        dataset_id=DATASET_ID,
        lab_group=LAB_GROUP,
        provider="Figshare",
        accession=ACCESSION,
        license_id=LICENSE,
        source_sha256=SOURCE_SHA256,
        source_locator=locator,
        sample_id=f"{DATASET_ID}:{sample_id}",
        biological_replicate=biological_replicate,
        technical_replicate=technical_replicate,
        condition=condition,
        cell_type="source_report_not_encoded",
        cell_state="source_condition_defined",
        modality=modality,
        observable=observable,
        value=float(value),
        unit="source_unit_unknown",
    )
    result.validate()
    return result


def _figure_1_replicates(values: list[list[object]]) -> list[Observation]:
    rows: list[Observation] = []
    blocks = ((0, "H3K9me2"), (6, "H3K9me3"), (12, "H3K27me3"))
    condition_rows = (4, 15, 26, 37)  # zero-based rows containing author labels
    field_specs = (
        (0, "nuclei_count", "counting_derived"),
        (1, "mean_if_intensity", "immunofluorescence_derived"),
        (2, "if_intensity_standard_deviation", "immunofluorescence_derived"),
        (3, "if_intensity_standard_error", "immunofluorescence_derived"),
    )
    for start_col, marker in blocks:
        for condition_row in condition_rows:
            condition = str(values[condition_row][start_col]).strip()
            for replicate_offset in range(6):
                row_index = condition_row + 2 + replicate_offset
                if row_index >= len(values) or not isinstance(values[row_index][start_col], (int, float)):
                    raise ValueError(f"missing Figure 1 replicate at row {row_index + 1}")
                sample = f"figure1:{marker}:{condition}:replicate-{replicate_offset + 1}"
                for field_offset, field, modality in field_specs:
                    value = values[row_index][start_col + field_offset]
                    if not isinstance(value, (int, float)):
                        raise ValueError(f"non-numeric Figure 1 value at row {row_index + 1}")
                    rows.append(_observation(
                        sheet="figure_1_replicates",
                        row=row_index,
                        column=start_col + field_offset,
                        sample_id=sample,
                        biological_replicate=str(replicate_offset + 1),
                        technical_replicate="aggregate_of_reported_nuclei",
                        condition=condition,
                        observable=f"{marker}_{field}",
                        value=value,
                        modality=modality,
                    ))
    return rows


def _figure_5_intensities(values: list[list[object]]) -> list[Observation]:
    rows: list[Observation] = []
    blocks = (
        (1, "chromocenter_size"),
        (7, "chromocenter_intensity"),
        (13, "chromocenter_intensity_per_area"),
        (19, "nuclear_periphery_intensity"),
    )
    for start_col, observable in blocks:
        for condition_offset in range(4):
            column = start_col + condition_offset
            condition = str(values[4][column]).strip()
            replicate = 0
            for row_index in range(5, len(values)):
                value = values[row_index][column] if column < len(values[row_index]) else None
                if value is None:
                    continue
                if not isinstance(value, (int, float)):
                    raise ValueError(f"unexpected text in Figure 5 intensities at row {row_index + 1}")
                replicate += 1
                rows.append(_observation(
                    sheet="figure_5_intensities",
                    row=row_index,
                    column=column,
                    sample_id=f"figure5:{observable}:{condition}:measurement-{replicate}",
                    biological_replicate="unknown",
                    technical_replicate=str(replicate),
                    condition=condition,
                    observable=observable,
                    value=value,
                ))
    return rows


def build(workbook_export: Path) -> list[Observation]:
    sheets = _sheet_map(workbook_export)
    rows = _figure_1_replicates(sheets["figure_1_replicates"])
    rows.extend(_figure_5_intensities(sheets["figure_5_intensities"]))
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
