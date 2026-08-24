#!/usr/bin/env python3
"""Normalize replicated C2C12 collective-motion and ECM-alignment evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from schema import Observation


DATASET = "dryad-08kprr59c-c2c12-supracontractility"
LAB = "anseth-white-skillin-lab-cu-boulder"
ACCESSION = "10.5061/dryad.08kprr59c"
ARTICLE = "10.1126/sciadv.adn0235"
EXPECTED_WORKBOOKS = {
    "Data_Fig2.xlsx": ("Data_Fig2", "B1:BO189"),
    "Data_Fig4.xlsx": ("Sheet1", "B1:Y183"),
}


def _column_name(index: int) -> str:
    output = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        output = chr(65 + remainder) + output
    return output


def _cell(row_index: int, column_index: int) -> str:
    # artifact-tool exports the values-only used range, which begins at column B.
    return f"{_column_name(column_index + 1)}{row_index + 1}"


def _finite(value: object) -> float | None:
    if not isinstance(value, (int, float)):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite author value: {value!r}")
    return result


def _load(
    workbook_export: Path, source_audit: Path,
) -> tuple[dict[str, list[list[object]]], dict[str, str], dict[str, set[str]], dict[str, set[str]]]:
    export = json.loads(workbook_export.read_text(encoding="utf-8"))
    audit = json.loads(source_audit.read_text(encoding="utf-8"))
    if audit.get("schema") != "aleph.outer_library.dryad_08kprr59c_source_audit.v1":
        raise ValueError("Dryad source audit schema changed")
    if audit.get("accession") != ACCESSION or audit.get("license") != "cc0-1.0":
        raise ValueError("Dryad source identity or license changed")
    sheets: dict[str, list[list[object]]] = {}
    observed: dict[str, tuple[str, str]] = {}
    for workbook in export.get("workbooks", []):
        name = Path(workbook["input"]).name
        if len(workbook["sheets"]) != 1:
            raise ValueError(f"unexpected sheet count: {name}")
        sheet = workbook["sheets"][0]
        observed[name] = (sheet["name"], sheet["address"])
        sheets[name] = sheet["values"]
    if observed != EXPECTED_WORKBOOKS:
        raise ValueError(f"Dryad workbook contract changed: {observed}")
    hashes = {
        name: audit["sources"][name]["sha256"] for name in EXPECTED_WORKBOOKS
    }
    red = {
        name: {item["cell"] for item in audit["annotations"][name]["author_excluded"]}
        for name in EXPECTED_WORKBOOKS
    }
    green = {
        name: {
            item["cell"]
            for item in audit["annotations"][name]["maximum_density_alignment"]
        }
        for name in EXPECTED_WORKBOOKS
    }
    return sheets, hashes, red, green


def _observation(
    *, source_name: str, source_sha: str, cell: str, sample_id: str,
    biological_replicate: str, technical_replicate: str, condition: str,
    state: str, modality: str, observable: str, value: float, unit: str,
    unit_cell: str, coordinate_name: str, coordinate_value: float,
    coordinate_unit: str, alignment_marker: bool = False,
) -> Observation:
    locator = f"xlsx:{source_name};sheet={EXPECTED_WORKBOOKS[source_name][0]};cell={cell}"
    if alignment_marker:
        locator += ";author_alignment=maximum_cell_density"
    identity = hashlib.sha256(
        f"{source_sha}:{cell}:{sample_id}:{observable}".encode()
    ).hexdigest()[:20]
    row = Observation(
        observation_id=f"{DATASET}:{identity}", dataset_id=DATASET,
        lab_group=LAB, provider="Dryad Digital Repository", accession=ACCESSION,
        license_id="cc0-1.0", source_sha256=source_sha, source_locator=locator,
        sample_id=sample_id, biological_replicate=biological_replicate,
        technical_replicate=technical_replicate, condition=condition,
        cell_type="C2C12_mouse_myoblast_myotube", cell_state=state,
        modality=modality, observable=observable, value=value, unit=unit,
        unit_source_sha256=source_sha,
        unit_source_locator=f"xlsx:{source_name};sheet={EXPECTED_WORKBOOKS[source_name][0]};cell={unit_cell}",
        coordinate_name=coordinate_name, coordinate_value=coordinate_value,
        coordinate_unit=coordinate_unit,
    )
    row.validate()
    return row


def _time_series_block(
    *, values: list[list[object]], source_name: str, source_sha: str,
    red: set[str], green: set[str], time_column: int, first_value_column: int,
    metrics: list[tuple[str, str]], biological: str, sample: str,
    technical: str, condition: str, state: str, modality: str,
    time_name: str = "time_after_confluence", start_row: int = 3,
) -> list[Observation]:
    output: list[Observation] = []
    for row_index in range(start_row, len(values)):
        row = values[row_index]
        if time_column >= len(row):
            break
        time = _finite(row[time_column])
        if time is None:
            break
        for offset, (observable, unit) in enumerate(metrics):
            column = first_value_column + offset
            cell = _cell(row_index, column)
            if cell in red or column >= len(row):
                continue
            value = _finite(row[column])
            if value is None:
                continue
            output.append(_observation(
                source_name=source_name, source_sha=source_sha, cell=cell,
                sample_id=sample, biological_replicate=biological,
                technical_replicate=technical, condition=condition, state=state,
                modality=modality, observable=observable, value=value, unit=unit,
                unit_cell=_cell(2, column), coordinate_name=time_name,
                coordinate_value=time, coordinate_unit="h",
                alignment_marker=cell in green,
            ))
    return output


def _figure2(
    values: list[list[object]], source_sha: str, red: set[str], green: set[str],
) -> list[Observation]:
    if values[0][10] != "Fig. 2C and 2D - green indicates maximum cell density for timeseries alignment; red indicates datapoint not included in analysis due to stage motion artifact or inability to link PIV analysis to subsequent frame":
        raise ValueError("Fig. 2 author annotation changed")
    metrics_m = [
        ("cell_number_density", "cells_per_mm2"),
        ("orientation_order_parameter", "dimensionless"),
        ("spatial_disorder_fraction", "dimensionless"),
        ("velocity_correlation_length", "micrometre"),
        ("cell_speed", "micrometre_per_hour"),
    ]
    metrics_i = metrics_m[:-1]
    output: list[Observation] = []
    for condition, blocks, metrics in (
        ("mLCN_1.0_mechanically_anisotropic", [12, 17, 22, 27, 32], metrics_m),
        ("iLCN_1.0_mechanically_isotropic", [38, 42, 46, 50], metrics_i),
    ):
        for replicate, first_column in enumerate(blocks, start=1):
            biological = f"author_independent_{condition}_replicate_{replicate}"
            sample = f"{DATASET}:fig2CD:{condition}:replicate-{replicate}"
            output += _time_series_block(
                values=values, source_name="Data_Fig2.xlsx", source_sha=source_sha,
                red=red, green=green, time_column=10,
                first_value_column=first_column, metrics=metrics,
                biological=biological, sample=sample,
                technical="author_live_field_PIV_time_series",
                condition=condition,
                state="myoblast_to_myotube_collective_transition",
                modality="particle_image_velocimetry_and_live_actin_derived",
            )
    # Figure 2F is a distinct whole-substrate experiment (N=3 mLCNs).
    metrics_f = [
        ("orientation_order_parameter", "dimensionless"),
        ("nematic_correlation_length", "micrometre"),
        ("velocity_correlation_length", "micrometre"),
    ]
    for replicate, first_column in enumerate((57, 60, 63), start=1):
        condition = "mLCN_1.0_mechanically_anisotropic_whole_substrate"
        biological = f"author_independent_fig2F_mLCN_replicate_{replicate}"
        output += _time_series_block(
            values=values, source_name="Data_Fig2.xlsx", source_sha=source_sha,
            red=red, green=green, time_column=56,
            first_value_column=first_column, metrics=metrics_f,
            biological=biological,
            sample=f"{DATASET}:fig2F:mLCN:replicate-{replicate}",
            technical="author_whole_substrate_live_actin_PIV_time_series",
            condition=condition, state="myotube_collective_flow_and_global_alignment",
            modality="particle_image_velocimetry_and_live_actin_derived",
        )
    return output


def _figure4(
    values: list[list[object]], source_sha: str, red: set[str], green: set[str],
) -> list[Observation]:
    if values[0][0] != "Fig. 4A, B - red indicates datapoint not included in analysis due to stage motion artifact or inability to link PIV analysis to subsequent frame":
        raise ValueError("Fig. 4 author annotation changed")
    metrics = [
        ("orientation_order_parameter", "dimensionless"),
        ("nematic_correlation_length", "micrometre"),
        ("velocity_correlation_length", "micrometre"),
        ("cell_speed", "micrometre_per_hour"),
    ]
    output: list[Observation] = []
    blocks = (
        ("vehicle_control", 1, 1), ("vehicle_control", 2, 5),
        ("blebbistatin_10_micromolar_daily", 1, 9),
        ("blebbistatin_10_micromolar_daily", 2, 13),
    )
    sample_by_block: dict[tuple[str, int], tuple[str, str]] = {}
    for condition, replicate, first_column in blocks:
        biological = f"author_independent_fig4_{condition}_mLCN_replicate_{replicate}"
        sample = f"{DATASET}:fig4:{condition}:replicate-{replicate}"
        sample_by_block[(condition, replicate)] = (sample, biological)
        output += _time_series_block(
            values=values, source_name="Data_Fig4.xlsx", source_sha=source_sha,
            red=red, green=green, time_column=0,
            first_value_column=first_column, metrics=metrics,
            biological=biological, sample=sample,
            technical="author_live_actin_PIV_time_series",
            condition=f"mLCN_1.0_{condition}",
            state=("myoblast_to_myotube_collective_transition" if condition == "vehicle_control" else "nonmuscle_myosin_II_inhibited_collective_transition"),
            modality="particle_image_velocimetry_and_live_actin_derived",
        )

    # Fig. 4C: author frequency values sum to 100 for every replicate.
    for column, (condition, replicate) in enumerate((
        ("vehicle_control", 1), ("vehicle_control", 2),
        ("blebbistatin_10_micromolar_daily", 1),
        ("blebbistatin_10_micromolar_daily", 2),
    ), start=20):
        sample, biological = sample_by_block[(condition, replicate)]
        frequencies = [
            _finite(row[column]) for row in values[3:]
            if column < len(row) and _finite(row[column]) is not None
        ]
        if len(frequencies) != 180 or not math.isclose(sum(frequencies), 100.0, abs_tol=1e-6):
            raise ValueError(f"Fig. 4C frequency contract changed: column {column}")
        for row_index in range(3, 183):
            angle = _finite(values[row_index][19])
            value = _finite(values[row_index][column])
            if angle is None or value is None:
                raise ValueError("missing Fig. 4C angle/frequency")
            output.append(_observation(
                source_name="Data_Fig4.xlsx", source_sha=source_sha,
                cell=_cell(row_index, column), sample_id=sample,
                biological_replicate=biological,
                technical_replicate="author_endpoint_fibronectin_orientation_distribution",
                condition=f"mLCN_1.0_{condition}",
                state=("differentiated_myotube_ECM" if condition == "vehicle_control" else "nonmuscle_myosin_II_inhibited_differentiated_myotube_ECM"),
                modality="immunofluorescence_fibronectin_orientation_derived",
                observable="fibronectin_orientation_frequency", value=value,
                unit="percent_frequency", unit_cell=_cell(2, column),
                coordinate_name="orientation_angle_to_nematic_director",
                coordinate_value=angle, coordinate_unit="degree",
            ))
    return output


def build(workbook_export: Path, source_audit: Path) -> tuple[list[Observation], dict[str, object]]:
    sheets, hashes, red, green = _load(workbook_export, source_audit)
    rows = _figure2(sheets["Data_Fig2.xlsx"], hashes["Data_Fig2.xlsx"], red["Data_Fig2.xlsx"], green["Data_Fig2.xlsx"])
    rows += _figure4(sheets["Data_Fig4.xlsx"], hashes["Data_Fig4.xlsx"], red["Data_Fig4.xlsx"], green["Data_Fig4.xlsx"])
    if len({row.observation_id for row in rows}) != len(rows):
        raise ValueError("duplicate Dryad C2C12 observation identity")
    report = {
        "schema": "aleph.outer_library.dryad_08kprr59c_manifest.v1",
        "dataset_id": DATASET, "accession": ACCESSION, "paper_doi": ARTICLE,
        "observations": len(rows), "samples": len({row.sample_id for row in rows}),
        "biological_replicates": len({row.biological_replicate for row in rows}),
        "modalities": sorted({row.modality for row in rows}),
        "author_excluded_cells_not_emitted": sum(len(cells) for cells in red.values()),
        "maximum_density_alignment_cells_preserved": len(green["Data_Fig2.xlsx"]),
        "replicate_contract": "each author-labelled LCN replicate is one biological sample; time points and observables remain nested technical observations",
        "aleph_authority": "none",
    }
    return rows, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook-export", type=Path, required=True)
    parser.add_argument("--source-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows, report = build(args.workbook_export, args.source_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(json.dumps(row.as_dict(), sort_keys=True) + "\n" for row in rows), encoding="utf-8")
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
