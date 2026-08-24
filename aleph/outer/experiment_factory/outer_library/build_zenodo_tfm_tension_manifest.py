#!/usr/bin/env python3
"""Normalize paired TFM and optical-tweezer figure tables from Zenodo 7432971."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re

import numpy as np

from schema import Observation


DATASET = "zenodo-7432971-tfm-tether"
LAB = "Cambridge_Franze_Keyser_collaboration"
ACCESSION = "10.5281/zenodo.7432971"
ARTICLE = "10.1093/pnasnexus/pgac299"
LICENSE = "cc-by-4.0"
EXPECTED_SHEETS = {
    "Figure2_Data.xlsx": {
        "OT_Neuron_Data": "A1:F25",
        "OT_Fibroblast_Data": "A1:H49",
        "TFM_Fibroblast_Data": "A1:E301",
        "TFM_Neuron_Pull_100Pa": "A1:D261",
        "TFM_Neuron_Pull_300Pa": "A1:D789",
    },
    "Figure3_Data.xlsx": {"OT_Fibroblast_Blebbistatin": "A1:H49"},
    "FigureS2_Data.xlsx": {"TFM_Neuron_Data": "A1:E1181"},
    "Kreysing_McHugh_key_resources.xlsx": {"Sheet1": "A1:E94"},
}


def _loaded(export_path: Path, receipt_path: Path) -> tuple[dict[str, list[list[object]]], dict[str, str]]:
    export = json.loads(export_path.read_text(encoding="utf-8"))
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    records = receipt.get("extracted_tables", {})
    if receipt.get("schema") != "aleph.outer_library.zenodo_7432971_figure_tables.v1":
        raise ValueError("Zenodo TFM receipt schema changed")
    if receipt.get("source_archive", {}).get("central_directory_range_verified") is not True:
        raise ValueError("Zenodo TFM central directory was not verified")
    sheets: dict[str, list[list[object]]] = {}
    hashes: dict[str, str] = {}
    observed: dict[str, dict[str, str]] = {}
    for workbook in export.get("workbooks", []):
        name = Path(workbook["input"]).name
        observed[name] = {sheet["name"]: sheet["address"] for sheet in workbook["sheets"]}
        hashes[name] = str(records[name]["sha256"])
        for sheet in workbook["sheets"]:
            if sheet["name"] in sheets:
                raise ValueError(f"duplicate sheet name: {sheet['name']}")
            sheets[sheet["name"]] = sheet["values"]
    if observed != EXPECTED_SHEETS:
        raise ValueError(f"Zenodo TFM workbook contract changed: {observed}")
    return sheets, hashes


def _observation(
    *, source_sha: str, locator: str, sample_id: str, biological: str,
    technical: str, condition: str, cell_type: str, state: str,
    modality: str, observable: str, value: float, unit: str,
) -> Observation:
    row = Observation(
        observation_id=f"{sample_id}:{observable}", dataset_id=DATASET,
        lab_group=LAB, provider="Zenodo", accession=ACCESSION,
        license_id=LICENSE, source_sha256=source_sha, source_locator=locator,
        sample_id=sample_id, biological_replicate=biological,
        technical_replicate=technical, condition=condition,
        cell_type=cell_type, cell_state=state, modality=modality,
        observable=observable, value=float(value), unit=unit,
        unit_source_sha256=source_sha, unit_source_locator=f"{locator}:header",
    )
    row.validate()
    return row


def _tfm_fibroblast(values: list[list[object]], source_sha: str) -> list[Observation]:
    if values[0] != ["Stiffness (kPa)", "TFM Mean (Pa)", "TFM Maximum (Pa)", "frame number", "file"]:
        raise ValueError("fibroblast TFM header changed")
    grouped: dict[tuple[int, str], list[tuple[int, float, float]]] = defaultdict(list)
    for stiffness, mean_stress, maximum_stress, frame, path in values[1:]:
        match = re.search(r"([0-9]+kPa)-([12])-Position([0-9]+)", str(path))
        if match is None or int(stiffness) not in {1, 10}:
            raise ValueError(f"unexpected fibroblast TFM path: {path}")
        cell = match.group(0)
        grouped[(int(stiffness), cell)].append((int(frame), float(mean_stress), float(maximum_stress)))
    if len(grouped) != 60 or any(sorted(row[0] for row in rows) != list(range(5)) for rows in grouped.values()):
        raise ValueError("fibroblast TFM cell/frame hierarchy changed")
    observations = []
    for (stiffness, cell), frames in sorted(grouped.items()):
        batch = cell.split("-")[1]
        sample_id = f"{DATASET}:3T3:{cell}"
        for observable, column in (("time_averaged_mean_traction_stress", 1), ("time_averaged_maximum_traction_stress", 2)):
            observations.append(_observation(
                source_sha=source_sha, locator=f"xlsx:TFM_Fibroblast_Data:file={cell}",
                sample_id=sample_id, biological=f"author_file_batch_{batch}_status_unresolved",
                technical="five_time_frames_averaged", condition=f"PAA_shear_modulus_{stiffness}_kPa",
                cell_type="NIH_3T3_fibroblast", state="untreated_adherent",
                modality="traction_force_microscopy_derived", observable=observable,
                value=np.mean([row[column] for row in frames]), unit="Pa",
            ))
    return observations


def _tfm_neuron(
    stress_values: list[list[object]], pull_sheets: dict[int, list[list[object]]], source_hashes: dict[str, str],
) -> list[Observation]:
    if stress_values[0] != ["Stiffness (Pa)", "TFM Mean (Pa)", "TFM Maximum (Pa)", "Frame number", "path"]:
        raise ValueError("neuron TFM stress header changed")
    stress: dict[tuple[int, str], list[tuple[int, float, float]]] = defaultdict(list)
    for stiffness, mean_stress, maximum_stress, frame, path in stress_values[1:]:
        cell = str(path).removesuffix("/results.h5")
        stress[(int(stiffness), cell)].append((int(frame), float(mean_stress), float(maximum_stress)))
    pulls: dict[tuple[int, str], list[tuple[float, float]]] = defaultdict(list)
    for stiffness, values in pull_sheets.items():
        expected_header = ["idx", f"{stiffness}Pa push (N)", f"{stiffness}Pa pull (N)", "path"]
        if values[0] != expected_header:
            raise ValueError(f"neuron pull header changed at {stiffness} Pa")
        for _, push, pull, path in values[1:]:
            cell = str(path).removesuffix("/results.h5")
            pulls[(stiffness, cell)].append((float(push), float(pull)))
    if set(stress) != set(pulls) or len(stress) != 132:
        raise ValueError("neuron TFM stress/net-force cell hierarchy changed")
    observations = []
    for (stiffness, cell), frames in sorted(stress.items()):
        date = re.search(r"_(20[0-9]{6})_", cell)
        if date is None:
            raise ValueError(f"neuron acquisition date missing: {cell}")
        sample_token = hashlib_sha(cell)
        sample_id = f"{DATASET}:RGC:{stiffness}Pa:{sample_token}"
        common = dict(
            sample_id=sample_id,
            biological=f"author_acquisition_date_{date.group(1)}_donor_status_unresolved",
            condition=f"PAA_shear_modulus_{stiffness}_Pa",
            cell_type="Xenopus_laevis_retinal_ganglion_growth_cone",
            state="untreated_adherent",
        )
        for observable, column in (("time_averaged_mean_traction_stress", 1), ("time_averaged_maximum_traction_stress", 2)):
            observations.append(_observation(
                source_sha=source_hashes["FigureS2_Data.xlsx"],
                locator=f"xlsx:TFM_Neuron_Data:path_sha={sample_token}",
                technical=f"{len(frames)}_available_time_frames_averaged", **common,
                modality="traction_force_microscopy_derived", observable=observable,
                value=np.mean([row[column] for row in frames]), unit="Pa",
            ))
        force_rows = pulls[(stiffness, cell)]
        for observable, column in (("time_averaged_growth_cone_push_force", 0), ("time_averaged_growth_cone_pull_force", 1)):
            observations.append(_observation(
                source_sha=source_hashes["Figure2_Data.xlsx"],
                locator=f"xlsx:TFM_Neuron_Pull_{stiffness}Pa:path_sha={sample_token}",
                technical=f"{len(force_rows)}_available_time_frames_averaged", **common,
                modality="traction_force_microscopy_derived", observable=observable,
                value=np.mean([row[column] for row in force_rows]) * 1e12, unit="pN",
            ))
    return observations


def hashlib_sha(value: str) -> str:
    import hashlib
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]


def _ot_sheet(
    values: list[list[object]], source_sha: str, cell_type: str,
    state: str, sheet: str, treated_columns_only: bool = False,
) -> list[Observation]:
    substrates, units, observables = values[:3]
    if any(unit != "pN" for unit in units):
        raise ValueError(f"OT units changed in {sheet}")
    observations = []
    for column, (substrate, observable) in enumerate(zip(substrates, observables, strict=True)):
        treated = "Bleb" in str(substrate)
        if treated_columns_only and not treated:
            continue
        if not treated_columns_only and treated:
            raise ValueError(f"unexpected treatment column in {sheet}")
        substrate_name = str(substrate).replace(" Bleb", "")
        for row_index, row in enumerate(values[3:], start=4):
            if column >= len(row) or not isinstance(row[column], (int, float)):
                continue
            cell_index = row_index - 3
            # Use the same sample index for steady-state/peak columns of one substrate.
            if sheet.endswith("Blebbistatin"):
                substrate_slot = column % 4
            else:
                substrate_slot = column // 2
            sample_id = f"{DATASET}:{sheet}:substrate-{substrate_slot}:cell-{cell_index}"
            observations.append(_observation(
                source_sha=source_sha, locator=f"xlsx:{sheet}:R{row_index}C{column + 1}",
                sample_id=sample_id, biological="unknown_not_reported",
                technical=f"author_cell_index_{cell_index}",
                condition=f"substrate_{substrate_name.replace(' ', '_')}",
                cell_type=cell_type, state=state,
                modality="optical_tweezers_membrane_tether_derived",
                observable=("membrane_tether_steady_state_force" if observable == "Steady State Force" else "membrane_tether_peak_detachment_force"),
                value=float(row[column]), unit="pN",
            ))
    return observations


def build(export_path: Path, receipt_path: Path) -> tuple[list[Observation], dict[str, object]]:
    sheets, hashes = _loaded(export_path, receipt_path)
    rows = _tfm_fibroblast(sheets["TFM_Fibroblast_Data"], hashes["Figure2_Data.xlsx"])
    rows += _tfm_neuron(
        sheets["TFM_Neuron_Data"],
        {100: sheets["TFM_Neuron_Pull_100Pa"], 300: sheets["TFM_Neuron_Pull_300Pa"]},
        hashes,
    )
    rows += _ot_sheet(
        sheets["OT_Fibroblast_Data"], hashes["Figure2_Data.xlsx"],
        "NIH_3T3_fibroblast", "untreated_adherent", "OT_Fibroblast_Data",
    )
    rows += _ot_sheet(
        sheets["OT_Neuron_Data"], hashes["Figure2_Data.xlsx"],
        "Xenopus_laevis_retinal_ganglion_axon", "untreated_adherent", "OT_Neuron_Data",
    )
    rows += _ot_sheet(
        sheets["OT_Fibroblast_Blebbistatin"], hashes["Figure3_Data.xlsx"],
        "NIH_3T3_fibroblast", "blebbistatin_20_uM", "OT_Fibroblast_Blebbistatin",
        treated_columns_only=True,
    )
    if len(rows) != 1_024 or len({row.sample_id for row in rows}) != 380:
        raise ValueError(f"Zenodo TFM/tether manifest contract changed: {len(rows)} rows")
    report = {
        "schema": "aleph.outer_library.zenodo_7432971_manifest.v1",
        "dataset_id": DATASET, "article_doi": ARTICLE,
        "observation_count": len(rows), "sample_count": len({row.sample_id for row in rows}),
        "cell_types": sorted({row.cell_type for row in rows}),
        "modalities": {name: sum(row.modality == name for row in rows) for name in sorted({row.modality for row in rows})},
        "source_workbooks_visual_header_check": True,
        "source_workbooks_formula_error_count": 0,
        "frames_collapsed_to_cell_level": True,
        "cells_are_not_called_biological_replicates": True,
        "independent_dataset_from_figshare_16826740": True,
        "independent_lab_from_figshare_16826740": True,
        "aleph_authority": "none", "may_select_aleph_parameter": False,
    }
    return rows, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook-export", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows, report = build(args.workbook_export, args.receipt)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(json.dumps(row.as_dict(), sort_keys=True) + "\n" for row in rows), encoding="utf-8")
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
