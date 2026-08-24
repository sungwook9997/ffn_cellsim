#!/usr/bin/env python3
"""Normalize experiment-grouped cell spreading and AFM tether-force data.

The source authors report cells (``n``) nested in independent experiments
(``N``).  This adapter preserves both levels and never promotes cells to
biological replicates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from schema import Observation


ACCESSION = "10.7554/eLife.72381"
DATASET_ID = "elife-72381-hela-mechano-osmotic"
LAB_GROUP = "cadart-paluch-mechano-osmotic-study"
FIG2_ARCHIVE_SHA256 = "e3d8035c3c008d7aed5bc88a817b9d94133e96aaf8ca7a0f9ecd574328011fc2"
FIG3_ARCHIVE_SHA256 = "a529b4fb30df09dc4ee44c1e4cb6b9babe301159f051f34d95a53d6a9636ceed"
FIG2_ARCHIVE_SIZE = 694_270
FIG3_ARCHIVE_SIZE = 465_889
CONTROL_INNER_SHA256 = "1010442e0b6eaf9b2f8f75b36ab4bb75c682e0f18450eb2d7ee280c7d68a0e24"
Y27_INNER_SHA256 = "78328f8fa0fe806c1f08d37d3a097e2a15914029ce473302b5f1e2aab4f12b1e"
TETHER_INNER_SHA256 = "90c3f361e274958b7f2bcdc13eb4d53f381b258ba3e2a5d5e7921fa226066a4c"
SPREADING_UNIT_LOCATORS = {
    "initial_spreading_area_rate": (
        "zip:Figure2/Figure2D/data_sum.xlsx!Sheet1!B1",
        "dA/dt i (µm2/min)",
    ),
    "initial_volume_flux": (
        "zip:Figure2/Figure2D/data_sum.xlsx!Sheet1!E1",
        "dV/dt i (µm3/min)",
    ),
}

CONTROL_EXCLUDED_CELL = (
    "20180220 hela spreading 20xPA FN RICM 2",
    1,
    6,
)

TETHER_CONDITIONS = {
    "control": "vehicle_control",
    "Y-27632": "Y27632_100_micromolar",
    "Y-27": "Y27632_100_micromolar",
    "Y-27632+CK-666": "Y27632_100_micromolar_plus_CK666_100_micromolar",
    "Y-27+CK-666": "Y27632_100_micromolar_plus_CK666_100_micromolar",
    "EIPA": "EIPA_50_micromolar",
    "EIPA+Y-27": "EIPA_50_micromolar_plus_Y27632_100_micromolar",
}

TETHER_SHEETS = {
    "30 min": ("spreading_phase_30_to_90_min", 5),
    "4h": ("steady_state_spread_4_to_5_hour", 5),
    "4h_round": ("steady_state_nonspread_micropattern_4_to_5_hour", 2),
}


def _payload(path: Path, expected_sha256: str) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("source_sha256") != expected_sha256:
        raise ValueError(f"unexpected source workbook hash for {path.name}")
    return payload


def _verify_source(path: Path, expected_sha256: str, expected_size: int) -> None:
    if path.stat().st_size != expected_size:
        raise ValueError(f"unexpected source archive size for {path.name}")
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha256:
        raise ValueError(f"unexpected source archive hash for {path.name}")


def _sheet(payload: dict, name: str) -> list[list[object]]:
    for sheet in payload["sheets"]:
        if sheet["name"] == name:
            return sheet["values"]
    raise ValueError(f"missing sheet {name!r}")


def _column_name(index: int) -> str:
    result = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(65 + remainder) + result
    return result


def _identifier(*parts: object) -> str:
    return hashlib.sha256(":".join(map(str, parts)).encode()).hexdigest()[:20]


def _spreading_replicate(experiment: str, condition: str) -> str:
    if condition == "Y27632_100_micromolar":
        match = re.match(r"(20\d{6})", experiment)
        if not match:
            raise ValueError(f"Y-27632 experiment lacks a date: {experiment}")
        # The workbook contains two labels for the 2018-05-17 experiment, while
        # the figure declares N=4.  Date-normalization reproduces the four
        # author experiments without treating the two labels as a fifth N.
        return f"author_experiment_date:{match.group(1)}"
    return f"author_experiment_label:{experiment}"


def _build_spreading(
    export: Path,
    *,
    expected_inner_sha256: str,
    condition: str,
    archive_sha256: str,
    member: str,
) -> list[Observation]:
    rows = _sheet(_payload(export, expected_inner_sha256), "Sheet1")
    if condition == "vehicle_control":
        expected_headers = ("experiment", "position", "cell", "dV/dt i", "dA/dt i")
        volume_column, area_column = 16, 17
    else:
        expected_headers = ("exp", "pos", "cell", "dV/dt i", "dA/dt i")
        volume_column, area_column = 12, 13
    header = rows[0]
    if (header[0], header[1], header[2], header[volume_column], header[area_column]) != expected_headers:
        raise ValueError(f"unexpected spreading schema for {condition}")

    output: list[Observation] = []
    excluded = 0
    for row_index, row in enumerate(rows[1:], start=2):
        if len(row) <= area_column:
            continue
        experiment, position, cell = row[0], row[1], row[2]
        volume, area = row[volume_column], row[area_column]
        if not all(isinstance(value, (int, float)) for value in (position, cell, volume, area)):
            continue
        if condition == "vehicle_control" and (experiment, position, cell) == CONTROL_EXCLUDED_CELL:
            # Figure 2D/data_sum declares control n=194.  This exact negative-
            # spreading-rate row is the sole exclusion that reconciles the raw
            # 195-row table to the published mean for both dA/dt and dV/dt.
            excluded += 1
            continue
        biological = _spreading_replicate(str(experiment), condition)
        sample = (
            f"{DATASET_ID}:spreading:{condition}:"
            f"{_identifier(experiment, position, cell, row_index)}"
        )
        for column, observable, value, unit in (
            (volume_column, "initial_volume_flux", volume, "cubic_micrometre_per_minute"),
            (area_column, "initial_spreading_area_rate", area, "square_micrometre_per_minute"),
        ):
            locator = f"zip:{member}!Sheet1!{_column_name(column)}{row_index}"
            unit_locator, _ = SPREADING_UNIT_LOCATORS[observable]
            observation = Observation(
                observation_id=f"elife-72381:{_identifier(archive_sha256, locator)}",
                dataset_id=DATASET_ID,
                lab_group=LAB_GROUP,
                provider="eLife Figure 2 source data",
                accession=ACCESSION,
                license_id="cc-by-4.0",
                source_sha256=archive_sha256,
                source_locator=locator,
                sample_id=sample,
                biological_replicate=biological,
                technical_replicate=(
                    f"author_cell:position_{position}:cell_{cell}:source_row_{row_index}"
                ),
                condition=condition,
                cell_type="HeLa Kyoto",
                cell_state="initial_spreading_on_fibronectin",
                modality="FXm_RICM_cell_spreading_derived",
                observable=observable,
                value=float(value),
                unit=unit,
                unit_source_sha256=FIG2_ARCHIVE_SHA256,
                unit_source_locator=unit_locator,
            )
            observation.validate()
            output.append(observation)
    expected_cells = 194 if condition == "vehicle_control" else 121
    if len(output) != 2 * expected_cells:
        raise ValueError(f"unexpected {condition} cell count: {len(output) // 2}")
    if condition == "vehicle_control" and excluded != 1:
        raise ValueError("published control exclusion was not found exactly once")
    expected_experiments = 3 if condition == "vehicle_control" else 4
    if len({row.biological_replicate for row in output}) != expected_experiments:
        raise ValueError(f"unexpected {condition} experiment count")
    return output


def _tether_replicate(identifier: object) -> str:
    text = str(int(identifier)) if isinstance(identifier, (int, float)) else str(identifier)
    if ".jpk-force" in text:
        text = text.split("_", 1)[0]
    return f"author_experiment:{text}"


def _build_tether(export: Path) -> list[Observation]:
    payload = _payload(export, TETHER_INNER_SHA256)
    output: list[Observation] = []
    for sheet_name, (state, expected_arms) in TETHER_SHEETS.items():
        rows = _sheet(payload, sheet_name)
        if len(rows) < 6:
            raise ValueError(f"empty tether-force sheet {sheet_name}")
        arm_count = 0
        for identifier_column in range(0, len(rows[4]), 3):
            force_column = identifier_column + 1
            if force_column >= len(rows[4]) or not rows[4][identifier_column]:
                continue
            source_condition = str(rows[4][identifier_column])
            if source_condition not in TETHER_CONDITIONS:
                raise ValueError(f"unknown tether condition {source_condition}")
            condition = TETHER_CONDITIONS[source_condition]
            expected_cells = int(rows[2][force_column])
            expected_experiments = int(rows[3][force_column])
            arm_rows = []
            for row_index, row in enumerate(rows[5:], start=6):
                if force_column >= len(row) or not isinstance(row[force_column], (int, float)):
                    continue
                identifier = row[identifier_column]
                biological = _tether_replicate(identifier)
                locator = (
                    f"zip:Figure3/Figure3E/tether_force_all.xlsx!"
                    f"{sheet_name}!{_column_name(force_column)}{row_index}"
                )
                sample = (
                    f"{DATASET_ID}:AFM:{state}:{condition}:"
                    f"{_identifier(sheet_name, identifier_column, row_index, identifier)}"
                )
                observation = Observation(
                    observation_id=f"elife-72381:{_identifier(FIG3_ARCHIVE_SHA256, locator)}",
                    dataset_id=DATASET_ID,
                    lab_group=LAB_GROUP,
                    provider="eLife Figure 3 source data",
                    accession=ACCESSION,
                    license_id="cc-by-4.0",
                    source_sha256=FIG3_ARCHIVE_SHA256,
                    source_locator=locator,
                    sample_id=sample,
                    biological_replicate=biological,
                    technical_replicate=f"author_cell_row:{row_index}",
                    condition=condition,
                    cell_type="HeLa Kyoto",
                    cell_state=state,
                    modality="AFM_membrane_tether_force_derived",
                    observable=f"membrane_tether_force:{state}",
                    value=float(row[force_column]),
                    unit="pN",
                )
                observation.validate()
                output.append(observation)
                arm_rows.append(observation)
            if len(arm_rows) != expected_cells:
                raise ValueError(
                    f"{sheet_name}/{source_condition}: n={len(arm_rows)}, expected {expected_cells}"
                )
            if len({row.biological_replicate for row in arm_rows}) != expected_experiments:
                raise ValueError(
                    f"{sheet_name}/{source_condition}: replicate grouping does not reproduce N"
                )
            arm_count += 1
        if arm_count != expected_arms:
            raise ValueError(f"unexpected arm count in {sheet_name}: {arm_count}")
    if len(output) != 314:
        raise ValueError(f"unexpected tether-force observation count: {len(output)}")
    return output


def build(
    figure2_source_zip: Path,
    figure3_source_zip: Path,
    control_export: Path,
    y27_export: Path,
    tether_export: Path,
) -> list[Observation]:
    _verify_source(figure2_source_zip, FIG2_ARCHIVE_SHA256, FIG2_ARCHIVE_SIZE)
    _verify_source(figure3_source_zip, FIG3_ARCHIVE_SHA256, FIG3_ARCHIVE_SIZE)
    rows = []
    rows.extend(_build_spreading(
        control_export,
        expected_inner_sha256=CONTROL_INNER_SHA256,
        condition="vehicle_control",
        archive_sha256=FIG2_ARCHIVE_SHA256,
        member="Figure2/Figure2D/control_dVdt_dAdt.xlsx",
    ))
    rows.extend(_build_spreading(
        y27_export,
        expected_inner_sha256=Y27_INNER_SHA256,
        condition="Y27632_100_micromolar",
        archive_sha256=FIG2_ARCHIVE_SHA256,
        member="Figure2/Figure2D/Y27_dVdt_dAdt.xlsx",
    ))
    rows.extend(_build_tether(tether_export))
    if len(rows) != 944 or len({row.sample_id for row in rows}) != 629:
        raise ValueError("unexpected final eLife 72381 manifest dimensions")
    if len({row.observation_id for row in rows}) != len(rows):
        raise ValueError("duplicate observation identifiers")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure2-source-zip", type=Path, required=True)
    parser.add_argument("--figure3-source-zip", type=Path, required=True)
    parser.add_argument("--control-spreading-export", type=Path, required=True)
    parser.add_argument("--y27-spreading-export", type=Path, required=True)
    parser.add_argument("--tether-export", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(
        args.figure2_source_zip,
        args.figure3_source_zip,
        args.control_spreading_export,
        args.y27_spreading_export,
        args.tether_export,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({
        "observations": len(rows),
        "samples": len({row.sample_id for row in rows}),
        "biological_replicates": len({row.biological_replicate for row in rows}),
        "modalities": sorted({row.modality for row in rows}),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
