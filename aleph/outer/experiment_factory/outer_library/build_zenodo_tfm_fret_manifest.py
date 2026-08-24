#!/usr/bin/env python3
"""Normalize exact source-data cells from the Zenodo 14692589 Prism archive."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

from schema import Observation


DATASET = "zenodo-14692589-tfm-fret"
LAB = "KU_Leuven_Aytekin_Ramon_Van_Oosterwyck"
ACCESSION = "10.5281/zenodo.14692589"
ARTICLE = "10.1038/s42003-026-09514-0"
LICENSE = "cc-by-4.0"
SOURCE_SHA256 = "ecd063c5175356b7c7fb231d94921a139d1725c83ab23af5bcb57e0a99c353d8"
EXPECTED_TABLES = {
    ("SD_Fig2c.prism", "Mean Traction vs Stiffness"):
        "FFDDB966-80F7-405F-A322-53C82EB2163A",
    ("SD_Fig2b_2f.prism", "VinTS vs TSMod"):
        "70780ACB-C630-4080-8FB3-8B4FC80C9DE7",
    ("SD_Fig4d_4e_4f.prism", "NCC 0614 35 Cell 9"):
        "A7DD9749-A8D1-41E4-A000-9A30C6A88745",
    ("SD_Fig4d_4e_4f.prism", "PCC 0719 4.5 Cell 10"):
        "E8A8BD59-7AF7-4BC4-AF9D-5BC0041F4C06",
    ("SD_Fig4d_4e_4f.prism", "UC 0628 35 Cell 7"):
        "2898EC36-CFC9-41C9-A5BC-1700FDD5D9A8",
}


def _title(value: object) -> str | None:
    if isinstance(value, dict):
        return str(value.get("string"))
    return None if value is None else str(value)


def _source(source_zip: Path, receipt_path: Path) -> tuple[zipfile.ZipFile, dict[str, object]]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") != "aleph.outer_library.zenodo_14692589_acquisition.v1":
        raise ValueError("Zenodo TFM--FRET acquisition schema changed")
    if receipt["source_file"]["sha256"] != SOURCE_SHA256:
        raise ValueError("Zenodo TFM--FRET receipt hash changed")
    payload = source_zip.read_bytes()
    if hashlib.sha256(payload).hexdigest() != SOURCE_SHA256:
        raise ValueError("Zenodo TFM--FRET source hash mismatch")
    return zipfile.ZipFile(io.BytesIO(payload)), receipt


def _sheet(
    outer: zipfile.ZipFile, prism_name: str, sheet_title: str
) -> tuple[list[list[str]], list[str], str]:
    member = f"Source-Data/{prism_name}"
    with zipfile.ZipFile(io.BytesIO(outer.read(member))) as prism:
        document = json.loads(prism.read("document.json"))
        for sheet_id in document["sheets"].get("data", []):
            sheet = json.loads(prism.read(f"data/sheets/{sheet_id}/sheet.json"))
            if sheet["title"] != sheet_title:
                continue
            table = sheet["table"]
            expected = EXPECTED_TABLES[(prism_name, sheet_title)]
            if table["uid"] != expected:
                raise ValueError(f"Prism table identity changed: {prism_name}/{sheet_title}")
            dataset_ids = []
            if table.get("xDataSet"):
                dataset_ids.append(table["xDataSet"])
            dataset_ids.extend(table.get("dataSets", []))
            labels = [
                _title(json.loads(prism.read(f"data/sets/{uid}.json")).get("title"))
                for uid in dataset_ids
            ]
            rows = list(csv.reader(io.StringIO(
                prism.read(f"data/tables/{expected}/data.csv").decode("utf-8-sig")
            )))
            return rows, [str(label) for label in labels], expected
    raise ValueError(f"Prism data sheet missing: {prism_name}/{sheet_title}")


def _observation(
    *, observation_id: str, locator: str, sample_id: str, technical: str,
    condition: str, cell_state: str, modality: str, observable: str,
    value: float, unit: str, unit_locator: str = "",
) -> Observation:
    row = Observation(
        observation_id=observation_id,
        dataset_id=DATASET,
        lab_group=LAB,
        provider="Zenodo",
        accession=ACCESSION,
        license_id=LICENSE,
        source_sha256=SOURCE_SHA256,
        source_locator=locator,
        sample_id=sample_id,
        biological_replicate=(
            "author_reports_at_least_three_replicates_but_cell_to_replicate_"
            "mapping_not_released"
        ),
        technical_replicate=technical,
        condition=condition,
        cell_type="vinculin_null_mouse_embryonic_fibroblast",
        cell_state=cell_state,
        modality=modality,
        observable=observable,
        value=float(value),
        unit=unit,
        unit_source_sha256=SOURCE_SHA256 if unit_locator else "",
        unit_source_locator=unit_locator,
    )
    row.validate()
    return row


def _column_values(rows: list[list[str]], column: int) -> list[tuple[int, float]]:
    values = []
    for row_index, row in enumerate(rows, start=1):
        if column >= len(row) or not row[column].strip():
            continue
        values.append((row_index, float(row[column])))
    return values


def _cell_level(outer: zipfile.ZipFile) -> list[Observation]:
    observations: list[Observation] = []
    rows, labels, table = _sheet(
        outer, "SD_Fig2c.prism", "Mean Traction vs Stiffness"
    )
    if labels != ["4.5 kPa", "13 kPa"]:
        raise ValueError(f"TFM stiffness labels changed: {labels}")
    for column, stiffness in enumerate(("4.5", "13")):
        values = _column_values(rows, column)
        if len(values) != (51, 38)[column]:
            raise ValueError("Figure 2c TFM cell count changed")
        for cell_index, (row_index, value) in enumerate(values, start=1):
            sample_id = f"{DATASET}:Fig2c:{stiffness}kPa:cell-{cell_index}"
            observations.append(_observation(
                observation_id=f"{sample_id}:mean_traction",
                locator=(
                    f"zip:Source-Data.zip!Source-Data/SD_Fig2c.prism!"
                    f"data/tables/{table}/data.csv:R{row_index}C{column + 1}"
                ),
                sample_id=sample_id,
                technical=f"author_cell_index_{cell_index}",
                condition=f"PAA_Young_modulus_{stiffness}_kPa",
                cell_state="VinTS_transfected_adherent",
                modality="traction_force_microscopy_derived",
                observable="cell_mask_mean_traction_stress",
                value=value,
                unit="Pa",
                unit_locator=(
                    "zip:Source-Data.zip!Source-Data/SD_Fig2c.prism!"
                    "data/sets:Mean Traction vs Stiffness"
                ),
            ))

    rows, labels, table = _sheet(
        outer, "SD_Fig2b_2f.prism", "VinTS vs TSMod"
    )
    if labels != ["TSMod", "VinTS", "TSMod", "VinTS", "TSMod", "VinTS"]:
        raise ValueError(f"Figure 2f sensor labels changed: {labels}")
    columns = (
        ("4.5", "TSMod", 20), ("4.5", "VinTS", 93),
        ("13", "TSMod", 20), ("13", "VinTS", 85),
        ("glass", "TSMod", 22), ("glass", "VinTS", 26),
    )
    for column, (substrate, sensor, expected) in enumerate(columns):
        values = _column_values(rows, column)
        if len(values) != expected:
            raise ValueError("Figure 2f FLIM-FRET cell count changed")
        condition = (
            "glass_substrate" if substrate == "glass"
            else f"PAA_Young_modulus_{substrate}_kPa"
        )
        for cell_index, (row_index, value) in enumerate(values, start=1):
            sample_id = (
                f"{DATASET}:Fig2f:{substrate}:{sensor}:cell-{cell_index}"
            )
            observations.append(_observation(
                observation_id=f"{sample_id}:FA_FRET",
                locator=(
                    f"zip:Source-Data.zip!Source-Data/SD_Fig2b_2f.prism!"
                    f"data/tables/{table}/data.csv:R{row_index}C{column + 1}"
                ),
                sample_id=sample_id,
                technical=f"author_cell_index_{cell_index}",
                condition=condition,
                cell_state=(
                    "vinculin_tension_sensor" if sensor == "VinTS"
                    else "force_insensitive_TSMod_control"
                ),
                modality="FLIM_FRET_vinculin_tension_sensor_derived",
                observable="cell_averaged_focal_adhesion_FRET_efficiency",
                value=value,
                unit="percent",
                unit_locator=(
                    "zip:Source-Data.zip!Source-Data/SD_Fig2b_2f.prism!"
                    "data/sets:VinTS vs TSMod"
                ),
            ))
    return observations


def _paired_focal_adhesions(outer: zipfile.ZipFile) -> list[Observation]:
    observations: list[Observation] = []
    sheets = (
        ("NCC 0614 35 Cell 9", "directly_correlated_cell", 46),
        ("PCC 0719 4.5 Cell 10", "inversely_correlated_cell", 77),
        ("UC 0628 35 Cell 7", "uncorrelated_cell", 123),
    )
    for title, category, expected in sheets:
        rows, labels, table = _sheet(
            outer, "SD_Fig4d_4e_4f.prism", title
        )
        if labels != ["Tractions (Pa)", "FRET%"]:
            raise ValueError(f"Figure 4 paired labels changed: {title}/{labels}")
        paired = [
            (row_index, float(row[0]), float(row[1]))
            for row_index, row in enumerate(rows, start=1)
            if len(row) > 1 and row[0].strip() and row[1].strip()
        ]
        if len(paired) != expected:
            raise ValueError(f"Figure 4 paired FA count changed: {title}")
        token = category.replace("_cell", "")
        sample_id = f"{DATASET}:Fig4:{token}:representative-cell"
        for fa_index, (row_index, traction, fret) in enumerate(paired, start=1):
            locator = (
                "zip:Source-Data.zip!Source-Data/SD_Fig4d_4e_4f.prism!"
                f"data/tables/{table}/data.csv:R{row_index}"
            )
            common = {
                "locator": locator,
                "sample_id": sample_id,
                "technical": f"author_focal_adhesion_index_{fa_index}",
                "condition": "author_representative_cell_substrate_not_released_in_sheet",
                "cell_state": category,
            }
            observations.append(_observation(
                observation_id=f"{sample_id}:FA-{fa_index}:traction",
                modality="paired_FA_TFM_FLIM_FRET_derived",
                observable="focal_adhesion_mean_traction_normalized_per_cell",
                value=traction,
                unit="author_normalized_fraction",
                **common,
            ))
            observations.append(_observation(
                observation_id=f"{sample_id}:FA-{fa_index}:FRET",
                modality="paired_FA_TFM_FLIM_FRET_derived",
                observable="focal_adhesion_FRET_efficiency_normalized_per_cell",
                value=fret,
                unit="author_normalized_fraction",
                **common,
            ))
    return observations


def build(
    source_zip: Path, receipt_path: Path
) -> tuple[list[Observation], dict[str, object]]:
    outer, receipt = _source(source_zip, receipt_path)
    try:
        rows = _cell_level(outer) + _paired_focal_adhesions(outer)
    finally:
        outer.close()
    if len(rows) != 847 or len({row.sample_id for row in rows}) != 358:
        raise ValueError(
            f"Zenodo 14692589 manifest contract changed: {len(rows)} rows, "
            f"{len({row.sample_id for row in rows})} samples"
        )
    if len({row.observation_id for row in rows}) != len(rows):
        raise ValueError("Zenodo 14692589 observation IDs are not unique")
    report = {
        "schema": "aleph.outer_library.zenodo_14692589_manifest.v1",
        "dataset_id": DATASET,
        "lab_group": LAB,
        "dataset_doi": ACCESSION,
        "article_doi": ARTICLE,
        "license_id": LICENSE,
        "source_sha256": receipt["source_file"]["sha256"],
        "observation_count": len(rows),
        "sample_count": len({row.sample_id for row in rows}),
        "cell_level_TFM_cells": 89,
        "cell_level_FLIM_FRET_cells": 266,
        "representative_paired_cells": 3,
        "paired_focal_adhesions": 246,
        "modalities": {
            name: sum(row.modality == name for row in rows)
            for name in sorted({row.modality for row in rows})
        },
        "author_reports_at_least_three_replicates": True,
        "cell_to_replicate_mapping_released": False,
        "cells_or_focal_adhesions_called_biological_replicates": False,
        "adversarial_source_findings": [
            "the Figure 2b cell-area source column has 90 stiff-substrate values while the paper caption declares n=85; cell area is withheld",
            "Figure 4 Prism columns are titled Tractions (Pa), but the paper and 0-to-1 values establish per-cell normalization; the manifest uses author_normalized_fraction",
            "the paper reports at least three repeated experiments and biologically independent samples but does not map released cells to replicate IDs",
        ],
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }
    return rows, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-zip", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows, report = build(args.source_zip, args.receipt)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"observations": len(rows), "samples": report["sample_count"]}))


if __name__ == "__main__":
    main()
