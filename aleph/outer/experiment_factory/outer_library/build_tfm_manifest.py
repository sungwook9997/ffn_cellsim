#!/usr/bin/env python3
"""Normalize the author-labelled Figshare TFM MATLAB source data.

The source files are tiny figure-source tables, not raw traction images.  This
adapter therefore publishes only quantities explicitly named in the MAT files.
Units are recovered from the matching publisher Figure 2E axes, not inferred
from numeric ranges.  Biological labels absent from the source remain unknown.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from schema import Observation


ACCESSION = "10.6084/m9.figshare.16826740"
DATASET_ID = "figshare-16826740-dictyostelium-tfm"
LAB_GROUP = "figshare-16826740-authors"
LICENSE = "cc-by-4.0"
ARTICLE_JATS_SHA256 = "361d64a0bad130edc02c7eacbedbea269bc029ccded6a8494348cca0903a16d8"
FIGURE_2_UNIT_LOCATOR = (
    "publisher_figure:https://link.springer.com/article/10.15252/msb.202110505"
    "#Fig2;panel:E;graphic:44320_2021_BFMSB202110505_Fig2_HTML.png"
)

FILE_ROLES = {
    "98ddd918d06a7693d7156e5d206ca6efc7dafd2f951118151cc4ca85416db988": "fan",
    "9bbf3cdc28d6409a38906d2be98f4cff496a9fb3f3fb72035968afab7cec53b8": "oscillatory",
    "4aeb6bd7829b4c7145a3ef0b555ffa0f4cbb17b8d3affa51c990b2359266e473": "amoeboid",
}

ROLE_FIELDS = {
    "fan": {
        "labels": "fanKind",
        "values": {
            "avgAfan": ("projected_area", "square_micrometre"),
            "avgFtotfan": ("mean_total_traction_force", "nN"),
            "avgFxfan": ("mean_abs_traction_force_x", "nN"),
            "avgFyfan": ("mean_abs_traction_force_y", "nN"),
            "avgLmajFan": ("major_axis_length", "micrometre"),
            "avgLminFan": ("minor_axis_length", "micrometre"),
        },
    },
    "oscillatory": {
        "labels": "osciKind",
        "values": {
            "avgAosci": ("projected_area", "square_micrometre"),
            "avgFtotosci": ("mean_total_traction_force", "nN"),
        },
    },
    "amoeboid": {
        "labels": "amoebKind",
        "values": {
            "Aavg": ("projected_area", "square_micrometre"),
            "Ftotavg": ("mean_total_traction_force", "nN"),
            "Fxavg": ("mean_abs_traction_force_x", "nN"),
            "Fyavg": ("mean_abs_traction_force_y", "nN"),
        },
    },
}


def _text_vector(value: object) -> list[str]:
    return [str(item).strip() for item in np.asarray(value, dtype=object).reshape(-1)]


def _object_vector(value: object) -> list[np.ndarray]:
    return [np.asarray(item, dtype=float).reshape(-1) for item in np.asarray(value, dtype=object).reshape(-1)]


def build(input_dir: Path) -> list[Observation]:
    rows: list[Observation] = []
    for sha256, role in FILE_ROLES.items():
        path = input_dir / sha256
        payload = path.read_bytes()
        if hashlib.sha256(payload).hexdigest() != sha256:
            raise ValueError(f"content hash mismatch: {path}")
        source = loadmat(path, squeeze_me=True, struct_as_record=False)
        spec = ROLE_FIELDS[role]
        labels = _text_vector(source[spec["labels"]])
        arrays = {key: _object_vector(source[key]) for key in spec["values"]}
        if any(len(values) != len(labels) for values in arrays.values()):
            raise ValueError(f"condition/value cardinality mismatch in {path}")
        for condition_index, condition in enumerate(labels):
            lengths = {len(values[condition_index]) for values in arrays.values()}
            if len(lengths) != 1:
                raise ValueError(f"unpaired observables for {condition!r} in {path}")
            for replicate_index in range(lengths.pop()):
                sample_id = (
                    f"{DATASET_ID}:{role}:condition-{condition_index}:cell-{replicate_index}"
                )
                for field, (observable, unit) in spec["values"].items():
                    value = float(arrays[field][condition_index][replicate_index])
                    row = Observation(
                        observation_id=f"{sha256[:12]}:{field}:{condition_index}:{replicate_index}",
                        dataset_id=DATASET_ID,
                        lab_group=LAB_GROUP,
                        provider="Figshare",
                        accession=ACCESSION,
                        license_id=LICENSE,
                        source_sha256=sha256,
                        source_locator=f"mat:{field}[{condition_index}][{replicate_index}]",
                        sample_id=sample_id,
                        biological_replicate="unknown",
                        technical_replicate=str(replicate_index),
                        condition=condition,
                        cell_type="Dictyostelium discoideum",
                        cell_state=role,
                        modality="traction_force_microscopy_derived",
                        observable=observable,
                        value=value,
                        unit=unit,
                        unit_source_sha256=ARTICLE_JATS_SHA256,
                        unit_source_locator=FIGURE_2_UNIT_LOCATOR,
                    )
                    row.validate()
                    rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(args.input_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({"observations": len(rows), "samples": len({r.sample_id for r in rows})}))


if __name__ == "__main__":
    main()
