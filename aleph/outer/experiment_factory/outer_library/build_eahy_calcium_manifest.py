#!/usr/bin/env python3
"""Normalize independent-lab EA.hy926 single-cell calcium time series."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from schema import Observation


SOURCE_SHA256 = "6927001307dfa22e6b8dc25f37cb25528c911d7911e5fcec79fa3cb957d816e4"
ACCESSION = "10.5281/zenodo.19890985"
DATASET_ID = "zenodo-19890985-eahy926-piezo1"


def build(source_csv: Path) -> list[Observation]:
    if hashlib.sha256(source_csv.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("EA.hy926 calcium CSV hash mismatch")
    with source_csv.open(encoding="utf-8-sig", newline="") as handle:
        table = list(csv.reader(handle))
    headers = table[0]
    rows = []
    curve_number = 0
    for column in range(1, len(headers)):
        condition = headers[column].strip()
        if not condition:
            continue
        curve_number += 1
        sample_id = f"{DATASET_ID}:{condition}:cell-{curve_number}"
        point_count = 0
        for row_index, raw in enumerate(table[1:], start=1):
            if column >= len(raw) or not raw[column].strip():
                break
            time = float(raw[0])
            value = float(raw[column])
            locator = f"csv:Figure 4 A ii.csv:R{row_index + 1}C{column + 1}"
            digest = hashlib.sha256(f"{SOURCE_SHA256}:{locator}".encode()).hexdigest()[:20]
            row = Observation(
                observation_id=f"eahy926:{digest}", dataset_id=DATASET_ID,
                lab_group="zenodo-19890985-authors", provider="Zenodo",
                accession=ACCESSION, license_id="cc-by-4.0",
                source_sha256=SOURCE_SHA256, source_locator=locator,
                sample_id=sample_id, biological_replicate="unknown_of_two_reported",
                technical_replicate=str(curve_number), condition=condition,
                cell_type="EA.hy926", cell_state=condition,
                modality="live_cell_calcium_imaging_derived",
                observable="intracellular_calcium_delta_F_over_F0", value=value,
                unit="ratio", coordinate_name="time", coordinate_value=time,
                coordinate_unit="s",
            )
            row.validate()
            rows.append(row)
            point_count += 1
        if point_count < 80:
            raise ValueError(f"curve {curve_number} has only {point_count} points")
    if curve_number != 625:
        raise ValueError(f"expected 625 author-labelled cell curves, found {curve_number}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(args.source_csv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({"observations": len(rows), "samples": len({r.sample_id for r in rows})}))


if __name__ == "__main__":
    main()
