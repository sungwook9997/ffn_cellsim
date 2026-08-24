#!/usr/bin/env python3
"""Normalize author source-data points for tendon mechanoculture experiments."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from schema import Observation


SOURCE_SHA256 = "ffc93fd3d25f768b1d32e5fbb25ac7d951df175c0b498920f9ba392a9924f36d"
ACCESSION = "10.1038/s41467-026-70395-2"
DATASET_ID = "nature-2026-tendon-mechanoculture"
LAB_GROUP = "hussien-snedeker-2026"


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
    *,
    sheet: str,
    row_index: int,
    column_index: int,
    condition: str,
    cell_type: str,
    cell_state: str,
    modality: str,
    observable: str,
    value: float,
    unit: str,
    point_kind: str,
    biological_replicate: str | None = None,
) -> Observation:
    # Exported sheets begin at Excel row 2 and column A.
    excel_row = row_index + 2
    cell = f"{_column_name(column_index)}{excel_row}"
    locator = f"xlsx:{sheet}!{cell}"
    identifier = hashlib.sha256(f"{SOURCE_SHA256}:{locator}:{observable}".encode()).hexdigest()[:20]
    result = Observation(
        observation_id=f"tendon-multimodal:{identifier}",
        dataset_id=DATASET_ID,
        lab_group=LAB_GROUP,
        provider="Nature Communications Source Data",
        accession=ACCESSION,
        license_id="cc-by-4.0",
        source_sha256=SOURCE_SHA256,
        source_locator=locator,
        sample_id=f"{DATASET_ID}:{sheet}:{cell}",
        biological_replicate=(biological_replicate or (
            f"author_tissue_point:{cell}" if point_kind == "tissue" else "unknown"
        )),
        technical_replicate=cell,
        condition=condition,
        cell_type=cell_type,
        cell_state=cell_state,
        modality=modality,
        observable=observable,
        value=float(value),
        unit=unit,
    )
    result.validate()
    return result


def _paired_columns(
    values: list[list[object]],
    *,
    sheet: str,
    start: int,
    stop: int,
    conditions: tuple[str, str],
    cell_type: str,
    states: tuple[str, str],
    modality: str,
    observable: str,
    unit: str,
    point_kind: str,
) -> list[Observation]:
    rows: list[Observation] = []
    for row_index in range(start, min(stop, len(values))):
        for column_index, condition in enumerate(conditions):
            value = values[row_index][column_index] if column_index < len(values[row_index]) else None
            if isinstance(value, (int, float)):
                rows.append(_observation(
                    sheet=sheet,
                    row_index=row_index,
                    column_index=column_index,
                    condition=condition,
                    cell_type=cell_type,
                    cell_state=states[column_index],
                    modality=modality,
                    observable=observable,
                    value=float(value),
                    unit=unit,
                    point_kind=point_kind,
                ))
    return rows


def _qpcr_gene_pair(
    values: list[list[object]],
    *,
    sheet: str,
    start: int,
    stop: int,
    columns: tuple[int, int],
    gene: str,
) -> list[Observation]:
    rows: list[Observation] = []
    for row_index in range(start, min(stop, len(values))):
        for condition_index, column_index in enumerate(columns):
            value = (
                values[row_index][column_index]
                if column_index < len(values[row_index]) else None
            )
            if not isinstance(value, (int, float)):
                continue
            condition = ("compliant_boundary", "rigid_boundary")[condition_index]
            rows.append(_observation(
                sheet=sheet,
                row_index=row_index,
                column_index=column_index,
                condition=condition,
                cell_type="rat tail tendon stromal cells",
                cell_state=("low_matrix_tension", "high_matrix_tension")[condition_index],
                modality="RT_qPCR_derived",
                observable=f"{gene}_relative_expression",
                value=float(value),
                unit="source_relative_expression_unit_unspecified",
                point_kind="tissue",
            ))
    return rows


def _human_qpcr_block(
    values: list[list[object]],
    *,
    gene: str,
    compliant_row: int,
    rigid_row: int,
) -> list[Observation]:
    rows: list[Observation] = []
    for disease, columns in (("healthy_control", range(1, 7)),
                             ("systemic_sclerosis", range(7, 12))):
        for row_index, boundary in ((compliant_row, "compliant_boundary"),
                                    (rigid_row, "rigid_boundary")):
            for column_index in columns:
                value = (
                    values[row_index][column_index]
                    if column_index < len(values[row_index]) else None
                )
                if not isinstance(value, (int, float)):
                    continue
                rows.append(_observation(
                    sheet="Figure 4G-I",
                    row_index=row_index,
                    column_index=column_index,
                    condition=f"{disease}:{boundary}",
                    cell_type="human tendon fibroblasts",
                    cell_state=("healthy" if disease == "healthy_control"
                                else "fibrotic_systemic_sclerosis"),
                    modality="RT_qPCR_derived",
                    observable=f"{gene}_relative_expression",
                    value=float(value),
                    unit="source_relative_expression_unit_unspecified",
                    point_kind="donor",
                    # The same workbook column is the only author-provided donor
                    # identity shared by compliant/rigid measurements.
                    biological_replicate=f"author_donor_column:{disease}:{_column_name(column_index)}",
                ))
    return rows


def build(workbook_export: Path) -> list[Observation]:
    sheets = _sheet_map(workbook_export)
    rows: list[Observation] = []

    rows.extend(_paired_columns(
        sheets["Figure 3C"], sheet="Figure 3C", start=2, stop=9,
        conditions=("compliant_boundary", "rigid_boundary"),
        cell_type="rat tail tendon stromal cells",
        states=("low_matrix_tension", "high_matrix_tension"),
        modality="cantilever_post_force_derived",
        observable="tissue_traction_force_per_cell", unit="nN_per_cell",
        point_kind="tissue",
    ))

    figure4 = sheets["Figure 4F"]
    rows.extend(_paired_columns(
        figure4, sheet="Figure 4F", start=4, stop=19,
        conditions=("healthy_control", "systemic_sclerosis"),
        cell_type="human tendon fibroblasts",
        states=("healthy", "fibrotic_systemic_sclerosis"),
        modality="cantilever_post_force_derived",
        observable="tissue_traction_force_per_cell", unit="nN_per_cell",
        point_kind="tissue",
    ))
    rows.extend(_paired_columns(
        figure4, sheet="Figure 4F", start=23, stop=47,
        conditions=("healthy_control", "systemic_sclerosis"),
        cell_type="human dermal fibroblasts",
        states=("healthy", "fibrotic_systemic_sclerosis"),
        modality="cantilever_post_force_derived",
        observable="tissue_traction_force_per_cell", unit="nN_per_cell",
        point_kind="tissue",
    ))

    figure2f = sheets["Figure 2F"]
    for start, stop, observable in (
        (2, 10, "alpha_SMA_relative_fluorescence"),
        (14, 22, "F_actin_relative_fluorescence"),
    ):
        rows.extend(_paired_columns(
            figure2f, sheet="Figure 2F", start=start, stop=stop,
            conditions=("vehicle_control", "TGF_beta_1"),
            cell_type="rat tail tendon stromal cells",
            states=("control", "profibrotic_TGF_beta_1_stimulated"),
            modality="immunofluorescence_derived", observable=observable,
            unit="relative_fluorescence_intensity", point_kind="tissue",
        ))

    figure3f = sheets["Figure 3F"]
    for start, stop, observable in (
        (2, 7, "alpha_SMA_relative_fluorescence"),
        (11, 16, "F_actin_relative_fluorescence"),
    ):
        rows.extend(_paired_columns(
            figure3f, sheet="Figure 3F", start=start, stop=stop,
            conditions=("compliant_boundary", "rigid_boundary"),
            cell_type="rat tail tendon stromal cells",
            states=("low_matrix_tension", "high_matrix_tension"),
            modality="immunofluorescence_derived", observable=observable,
            unit="relative_fluorescence_intensity", point_kind="tissue",
        ))

    rows.extend(_paired_columns(
        sheets["Figure 2H"], sheet="Figure 2H", start=2,
        stop=len(sheets["Figure 2H"]),
        conditions=("vehicle_control", "TGF_beta_1"),
        cell_type="rat tail tendon stromal cells",
        states=("control", "profibrotic_TGF_beta_1_stimulated"),
        modality="nuclear_morphology_derived", observable="nuclear_circularity",
        unit="dimensionless", point_kind="nucleus",
    ))

    rows.extend(_qpcr_gene_pair(
        sheets["Figure 3G"], sheet="Figure 3G", start=2, stop=8,
        columns=(0, 1), gene="Acta2",
    ))
    figure3jk = sheets["Figure 3J-K"]
    for start, stop, columns, gene in (
        (3, 9, (0, 1), "Col1a1"),
        (3, 9, (4, 5), "Col3a1"),
        (3, 9, (8, 9), "Dcn"),
        (13, 19, (0, 1), "Scx"),
        (13, 19, (4, 5), "Tnmd"),
        (13, 19, (8, 9), "Mkx"),
    ):
        rows.extend(_qpcr_gene_pair(
            figure3jk, sheet="Figure 3J-K", start=start, stop=stop,
            columns=columns, gene=gene,
        ))

    figure4_qpcr = sheets["Figure 4G-I"]
    for gene, compliant_row, rigid_row in (
        ("COL1A1", 3, 4),
        ("FN_EDA", 9, 10),
        ("IL6", 15, 16),
        ("CCL2", 21, 22),
        ("ACTA2", 27, 28),
        ("ICAM1", 33, 34),
    ):
        rows.extend(_human_qpcr_block(
            figure4_qpcr, gene=gene,
            compliant_row=compliant_row, rigid_row=rigid_row,
        ))
    rows.extend(_paired_columns(
        sheets["Figure 3I"], sheet="Figure 3I", start=2,
        stop=len(sheets["Figure 3I"]),
        conditions=("compliant_boundary", "rigid_boundary"),
        cell_type="rat tail tendon stromal cells",
        states=("low_matrix_tension", "high_matrix_tension"),
        modality="nuclear_morphology_derived", observable="nuclear_circularity",
        unit="dimensionless", point_kind="nucleus",
    ))

    identifiers = [row.observation_id for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("duplicate tendon multimodal observation IDs")
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
    print(json.dumps({
        "observations": len(rows),
        "samples": len({row.sample_id for row in rows}),
        "modalities": sorted({row.modality for row in rows}),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
