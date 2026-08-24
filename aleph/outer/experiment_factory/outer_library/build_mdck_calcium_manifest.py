#!/usr/bin/env python3
"""Normalize the author-described MDCK-II calcium model-fit source data."""

from __future__ import annotations

import argparse
import hashlib
import json
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.io import loadmat

from schema import Observation


SOURCE_SHA256 = "106a9b5fd0a38f245d989452a39f34d8144216185b1c1687250fd9d0aa3fc775"
ACCESSION = "10.5281/zenodo.8215150"
DATASET_ID = "zenodo-8215150-mdck-mechanocalcium"


def build(source_zip: Path) -> list[Observation]:
    if hashlib.sha256(source_zip.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("MDCK source archive hash mismatch")
    with ZipFile(source_zip) as archive:
        payload = loadmat(BytesIO(archive.read("model_fit.mat")), squeeze_me=True)
    times = np.asarray(payload["Trange"], dtype=float)
    traces = np.asarray(payload["target_trace"], dtype=float)
    if traces.shape != (7, 118) or times.shape != (118,):
        raise ValueError("unexpected MDCK source shape")
    rows = []
    for trace_index, trace in enumerate(traces):
        population = "target_cell" if trace_index < 5 else "neighbor_cell"
        quantile = trace_index + 1 if trace_index < 5 else trace_index - 4
        condition = f"{population}_author_amplitude_quantile_{quantile}"
        sample_id = f"{DATASET_ID}:{condition}"
        for time_index, (time, value) in enumerate(zip(times, trace, strict=True)):
            locator = f"zip:model_fit.mat:target_trace[{trace_index},{time_index}]"
            digest = hashlib.sha256(f"{SOURCE_SHA256}:{locator}".encode()).hexdigest()[:20]
            row = Observation(
                observation_id=f"mdck:{digest}", dataset_id=DATASET_ID,
                lab_group="zenodo-8215150-authors", provider="Zenodo",
                accession=ACCESSION, license_id="cc-by-4.0",
                source_sha256=SOURCE_SHA256, source_locator=locator,
                sample_id=sample_id, biological_replicate="author_aggregate",
                technical_replicate=str(quantile), condition=condition,
                cell_type="MDCK-II", cell_state=population,
                modality="live_cell_calcium_imaging_derived",
                observable="normalized_calcium_signal", value=float(value), unit="ratio",
                coordinate_name="time", coordinate_value=float(time), coordinate_unit="s",
            )
            row.validate()
            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-zip", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(args.source_zip)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({"observations": len(rows), "samples": len({r.sample_id for r in rows})}))


if __name__ == "__main__":
    main()
