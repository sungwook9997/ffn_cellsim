#!/usr/bin/env python3
"""Build leakage-safe cell summaries and a compact row tensor from Dryad 9jh6m.

The released CSV rows are repeated observations of tracked cells, not
independent samples.  The JSONL manifest therefore emits one median per
date/condition/cell/feature and retains the author motion-mode proportions.
The optional NPZ retains every author row for grouped representation learning;
its cell and experimental-date identifiers are mandatory split groups.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from schema import Observation


ACCESSION = "10.5061/dryad.9jh6m"
SOURCE_SHA256 = "8d8aecf4dec41e058a5bd789fff3fc483ff827ae7111ef856c69d462f9f034d5"
OFFICIAL_SIZE = 5_450_110
LAB_GROUP = "stromblad-lab-karolinska"
PROVIDER = "Dryad Digital Repository"
PREFIX = "Shafqat-Abbasi et al Source Datasets/"
INDEX_COLUMNS = {
    "Perturbation", "Date", "Celltrace", "Motion mode",
    "CellCenterRow", "CellCenterCol",
}
BEHAVIORAL_FEATURES = (
    "Cell Protrusion Area",
    "Cell Retraction Area",
    "Stable Cell Area",
    "Short Lived Cell Area",
    "Cell Speed",
)
DATASETS = {
    "fibronectin": {
        "member": PREFIX + "Fibronectin_Modulation_Dataset.csv",
        "dataset_id": "dryad-9jh6m-fibronectin-modulation",
        "conditions": {
            "1": "fibronectin_2.5_microgram_per_ml",
            "2": "fibronectin_10_microgram_per_ml",
        },
        "expected_rows": 9_579,
        "expected_cells": 162,
        "expected_features": 60,
    },
    "rock": {
        "member": PREFIX + "ROCK_Modulation_Dataset.csv",
        "dataset_id": "dryad-9jh6m-ROCK-modulation",
        "conditions": {"1": "DMSO_control", "2": "Y27632_6_micromolar"},
        "expected_rows": 3_650,
        "expected_cells": 51,
        "expected_features": 54,
    },
    "talin": {
        "member": PREFIX + "Talin_Modulation_Dataset.csv",
        "dataset_id": "dryad-9jh6m-talin-modulation",
        "conditions": {"1": "control_siRNA", "2": "talin_1_siRNA"},
        "expected_rows": 9_418,
        "expected_cells": 154,
        "expected_features": 60,
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _clean_header(value: str) -> str:
    return " ".join(value.strip().split())


def _float(value: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def _unit(observable: str) -> str:
    lowered = observable.lower()
    if any(token in lowered for token in (
        "eccentricity", "solidity", "compactness", "colocalization",
    )):
        return "dimensionless"
    return "author_processed_source_unit_unspecified"


def _state(mode_counts: Counter[str] | None) -> str:
    if mode_counts is None:
        return "author_motion_mode_not_reported"
    observed = {mode for mode, count in mode_counts.items() if count}
    mapping = {
        frozenset({"0"}): "author_motion_mode_undetermined",
        frozenset({"1"}): "author_motion_mode_discontinuous",
        frozenset({"2"}): "author_motion_mode_continuous",
    }
    return mapping.get(frozenset(observed), "author_motion_mode_mixed_over_time")


def _load_archive(source_zip: Path) -> tuple[str, dict[str, dict]]:
    digest = _sha256(source_zip)
    if source_zip.stat().st_size != OFFICIAL_SIZE or digest != SOURCE_SHA256:
        raise ValueError(
            f"official Dryad archive mismatch: size={source_zip.stat().st_size}, sha256={digest}"
        )
    loaded: dict[str, dict] = {}
    with ZipFile(source_zip) as archive:
        names = set(archive.namelist())
        for key, spec in DATASETS.items():
            if spec["member"] not in names:
                raise ValueError(f"archive member missing: {spec['member']}")
            text = io.TextIOWrapper(archive.open(spec["member"]), encoding="utf-8-sig")
            reader = csv.DictReader(text)
            if reader.fieldnames is None:
                raise ValueError(f"CSV header missing: {spec['member']}")
            raw_to_clean = {name: _clean_header(name) for name in reader.fieldnames}
            rows = [
                {raw_to_clean[name]: value for name, value in row.items()}
                for row in reader
            ]
            feature_names = [
                raw_to_clean[name] for name in reader.fieldnames
                if raw_to_clean[name] not in INDEX_COLUMNS
            ]
            cell_keys = {
                (row["Perturbation"], row["Date"], row["Celltrace"])
                for row in rows
            }
            if len(rows) != spec["expected_rows"]:
                raise ValueError(f"unexpected row count for {key}: {len(rows)}")
            if len(cell_keys) != spec["expected_cells"]:
                raise ValueError(f"unexpected tracked-cell count for {key}: {len(cell_keys)}")
            if len(feature_names) != spec["expected_features"]:
                raise ValueError(f"unexpected feature count for {key}: {len(feature_names)}")
            unexpected = sorted({row["Perturbation"] for row in rows} - spec["conditions"].keys())
            if unexpected:
                raise ValueError(f"unexpected conditions for {key}: {unexpected}")
            loaded[key] = {"rows": rows, "features": feature_names, "spec": spec}
    return digest, loaded


def build(source_zip: Path) -> list[Observation]:
    digest, datasets = _load_archive(source_zip)
    output: list[Observation] = []
    for dataset_key, loaded in datasets.items():
        spec = loaded["spec"]
        grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
        for row in loaded["rows"]:
            grouped[(row["Perturbation"], row["Date"], row["Celltrace"])].append(row)
        for (author_condition, date, celltrace), rows in sorted(grouped.items()):
            condition = spec["conditions"][author_condition]
            modes = (
                Counter(row["Motion mode"] for row in rows)
                if "Motion mode" in rows[0] else None
            )
            state = _state(modes)
            biological = f"author_experimental_date:{date}"
            technical = f"tracked_cell:{celltrace}"
            sample = f"{spec['dataset_id']}:{condition}:{biological}:{technical}"
            for feature in loaded["features"]:
                values = np.asarray([_float(row[feature]) for row in rows], dtype=np.float64)
                finite = values[np.isfinite(values)]
                if not len(finite):
                    continue
                locator = (
                    f"zip:{spec['member']};Perturbation={author_condition};Date={date};"
                    f"Celltrace={celltrace};column={feature};statistic=nanmedian;"
                    f"finite_n={len(finite)};row_n={len(rows)}"
                )
                identifier = hashlib.sha256(f"{digest}:{locator}".encode()).hexdigest()[:20]
                observation = Observation(
                    observation_id=f"dryad-9jh6m:{identifier}",
                    dataset_id=spec["dataset_id"],
                    lab_group=LAB_GROUP,
                    provider=PROVIDER,
                    accession=ACCESSION,
                    license_id="cc0-1.0",
                    source_sha256=digest,
                    source_locator=locator,
                    sample_id=sample,
                    biological_replicate=biological,
                    technical_replicate=technical,
                    condition=condition,
                    cell_type="H1299-PL EGFP-Paxillin RubyRed-LifeAct",
                    cell_state=state,
                    modality="multiscale_live_cell_imaging_derived",
                    observable=feature,
                    value=float(np.median(finite)),
                    unit=_unit(feature),
                )
                observation.validate()
                output.append(observation)
            if modes is not None:
                for mode, label in (
                    ("0", "undetermined"), ("1", "discontinuous"), ("2", "continuous")
                ):
                    locator = (
                        f"zip:{spec['member']};Perturbation={author_condition};Date={date};"
                        f"Celltrace={celltrace};column=Motion mode;"
                        f"statistic=fraction_equal_{mode};row_n={len(rows)}"
                    )
                    identifier = hashlib.sha256(f"{digest}:{locator}".encode()).hexdigest()[:20]
                    observation = Observation(
                        observation_id=f"dryad-9jh6m:{identifier}",
                        dataset_id=spec["dataset_id"],
                        lab_group=LAB_GROUP,
                        provider=PROVIDER,
                        accession=ACCESSION,
                        license_id="cc0-1.0",
                        source_sha256=digest,
                        source_locator=locator,
                        sample_id=sample,
                        biological_replicate=biological,
                        technical_replicate=technical,
                        condition=condition,
                        cell_type="H1299-PL EGFP-Paxillin RubyRed-LifeAct",
                        cell_state=state,
                        modality="author_supervised_migration_mode_derived",
                        observable=f"author_motion_mode_fraction_{label}",
                        value=float(modes[mode] / len(rows)),
                        unit="fraction",
                    )
                    observation.validate()
                    output.append(observation)
    return output


def build_tensor(source_zip: Path, output: Path) -> dict:
    digest, datasets = _load_archive(source_zip)
    features = sorted({name for item in datasets.values() for name in item["features"]})
    feature_index = {name: index for index, name in enumerate(features)}
    matrices, dataset_names, conditions, dates, cells, modes, row_indices = [], [], [], [], [], [], []
    for dataset_key, loaded in datasets.items():
        counters: Counter[tuple[str, str, str]] = Counter()
        for row in loaded["rows"]:
            vector = np.full(len(features), np.nan, dtype=np.float64)
            for feature in loaded["features"]:
                vector[feature_index[feature]] = _float(row[feature])
            key = (row["Perturbation"], row["Date"], row["Celltrace"])
            matrices.append(vector)
            dataset_names.append(dataset_key)
            conditions.append(loaded["spec"]["conditions"][row["Perturbation"]])
            dates.append(row["Date"])
            cells.append(f"{dataset_key}:{row['Perturbation']}:{row['Date']}:{row['Celltrace']}")
            modes.append(int(row.get("Motion mode", -1)))
            row_indices.append(counters[key])
            counters[key] += 1
    matrix = np.stack(matrices)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        schema=np.asarray("aleph.outer_library.grouped_cell_rows.v1"),
        source_sha256=np.asarray(digest),
        feature_names=np.asarray(features),
        values=matrix,
        observed=np.isfinite(matrix),
        dataset=np.asarray(dataset_names),
        condition=np.asarray(conditions),
        experimental_date=np.asarray(dates),
        tracked_cell=np.asarray(cells),
        author_motion_mode=np.asarray(modes, dtype=np.int8),
        row_within_tracked_cell=np.asarray(row_indices, dtype=np.int32),
        behavioral_features=np.asarray(BEHAVIORAL_FEATURES),
    )
    return {
        "rows": int(len(matrix)),
        "tracked_cells": int(len(set(cells))),
        "features": int(len(features)),
        "finite_values": int(np.isfinite(matrix).sum()),
        "output_sha256": _sha256(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-zip", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tensor-output", type=Path)
    args = parser.parse_args()
    rows = build(args.source_zip)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    result = {
        "manifest_observations": len(rows),
        "tracked_cells": len({row.sample_id for row in rows}),
        "experimental_dates": len({row.biological_replicate for row in rows}),
        "datasets": len({row.dataset_id for row in rows}),
    }
    if args.tensor_output:
        result["row_tensor"] = build_tensor(args.source_zip, args.tensor_output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
