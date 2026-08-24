#!/usr/bin/env python3
"""Normalize per-cell Figure 4 trajectory observables from eLife 71032."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from schema import Observation


SOURCE_SHA256 = "03f16d1d25bc603eb26c76542dfa2130ce6918fca7eb62a99a9ff5e1f10202f7"
ACCESSION = "10.7554/eLife.71032"
DATASET_ID = "elife-71032-nih3t3-3d-cdm-migration"
LAB_GROUP = "riveline-godeau-strasbourg-curie"
SPEED_SHEET = "Figure 4 Supplementary 1i"
PERIOD_SHEET = "Figure 4 Supplementary 1j"
CONDITIONS = (
    (0, "untreated_control", 24),
    (1, "C8_BPA_100_micromolar", 7),
    (2, "CK666_50_micromolar", 6),
    (3, "ML7_10_micromolar", 7),
    (4, "Y27632_10_micromolar", 9),
)


def _sheet_map(workbook_export: Path) -> dict[str, list[list[object]]]:
    payload = json.loads(workbook_export.read_text(encoding="utf-8"))
    if payload.get("source_sha256") != SOURCE_SHA256:
        raise ValueError(
            "eLife workbook export does not identify the verified Figure 4 source"
        )
    return {sheet["name"]: sheet["values"] for sheet in payload["sheets"]}


def _column_name(index: int) -> str:
    result = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(65 + remainder) + result
    return result


def build(workbook_export: Path) -> list[Observation]:
    sheets = _sheet_map(workbook_export)
    if SPEED_SHEET not in sheets or PERIOD_SHEET not in sheets:
        raise ValueError("eLife Figure 4 speed/period sheets are missing")
    speed, period = sheets[SPEED_SHEET], sheets[PERIOD_SHEET]
    output: list[Observation] = []
    for column, condition, expected in CONDITIONS:
        speed_values = [
            (row_index, row[column])
            for row_index, row in enumerate(speed[1:], start=1)
            if column < len(row) and isinstance(row[column], (int, float))
        ]
        period_values = [
            (row_index, row[column])
            for row_index, row in enumerate(period[1:], start=1)
            if column < len(row) and isinstance(row[column], (int, float))
        ]
        if len(speed_values) != expected or len(period_values) != expected:
            raise ValueError(
                f"unexpected {condition} count: speed={len(speed_values)}, "
                f"period={len(period_values)}, expected={expected}"
            )
        # The source sheets use the same row order for the two observables.
        if [index for index, _ in speed_values] != [index for index, _ in period_values]:
            raise ValueError(f"speed/period row mismatch for {condition}")
        for cell_index, ((speed_row, speed_value), (period_row, period_value)) in enumerate(
            zip(speed_values, period_values), start=1
        ):
            sample = f"{DATASET_ID}:{condition}:author_cell_point_{cell_index}"
            for sheet, row_index, observable, value, unit in (
                (SPEED_SHEET, speed_row, "persistent_speed", speed_value, "micrometre_per_minute"),
                (PERIOD_SHEET, period_row, "speed_oscillation_period", period_value, "minute"),
            ):
                cell = f"{_column_name(column)}{row_index + 1}"
                locator = f"xlsx:{sheet}!{cell}"
                identifier = hashlib.sha256(
                    f"{SOURCE_SHA256}:{locator}:{observable}".encode()
                ).hexdigest()[:20]
                row = Observation(
                    observation_id=f"elife-71032:{identifier}",
                    dataset_id=DATASET_ID,
                    lab_group=LAB_GROUP,
                    provider="eLife Figure 4 source data",
                    accession=ACCESSION,
                    license_id="cc-by-4.0",
                    source_sha256=SOURCE_SHA256,
                    source_locator=locator,
                    sample_id=sample,
                    biological_replicate="not_reported",
                    technical_replicate=f"author_cell_point:{cell_index}",
                    condition=condition,
                    cell_type="NIH3T3 fibroblast",
                    cell_state=(
                        "untreated_3D_CDM_migration" if condition == "untreated_control"
                        else f"pharmacologically_perturbed_3D_CDM_migration:{condition}"
                    ),
                    modality="3D_CDM_cell_migration_trajectory_derived",
                    observable=observable,
                    value=float(value),
                    unit=unit,
                )
                row.validate()
                output.append(row)
    if len(output) != 106 or len({row.sample_id for row in output}) != 53:
        raise ValueError("unexpected final eLife 71032 manifest size")
    return output


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
        "cells": len({row.sample_id for row in rows}),
        "conditions": sorted({row.condition for row in rows}),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
