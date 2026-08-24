#!/usr/bin/env python3
"""Build a provenance-checked hdF collective-velocity representation dataset.

The public release contains spatial velocity vectors and longitudinal author
summaries, but no identifiers for independent cultures.  Vector points,
frames, and time points therefore remain technical observations nested in one
unresolved study group; none is promoted to a biological replicate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from schema import Observation


ACCESSION = "10.25349/D96W58"
DATASET_ID = "dryad-d96w58-hdf-collective-velocity"
LAB_GROUP = "valentine-lab-ucsb"
PROVIDER = "Dryad Digital Repository"
METADATA_SHA256 = "6de2283ac2032abcc10370921b2cff379b5759ca4e7c57b2a4517308791e4acd"
FILES = {
    "Fig5ab.csv": (160_186, "6a371040a1fc7c84f16c1108a4a7fa10", "2a2b3987049bea339c1b0f54b1636beb0abbcc21a9bc73524514ddf165ff695b"),
    "Fig5cd.csv": (257_287, "19e33f9e219e64e32afb503ac3e01967", "fe9f687c55cd76a82980a32f37be20e4f38b74dec9f702794f44bd2c21c35830"),
    "FigS19bc_1uMFAKi.csv": (3_757, "b93e4005ac4e8d0e4573041a6d949ac3", "8ed67281791e2527f150755645630aef0722034da5aea97c6f74765c598cfeb0"),
    "FigS19d_1uMFAKi.csv": (4_746, "97e4d9ad8c55501b13e842a49802d741", "447269a8bbccb86c4844f6ae3b575cd1562886eb44ffe3eedd302a5ba1c4730c"),
    "FigS19fg.csv": (3_790, "edd57c441fe90e47cda26314022aba7c", "b89ecf6ad43499b899a8fbf1b0de316ffc99e5a4718b71678c5b565393ba4664"),
    "FigS7e.csv": (29_040, "300d332dd842c355058104a42dcabdee", "3ab17090b859fa47bb6f5950449605998aa57defd3cd815ee88fa5fbba780eae"),
    "FigS7f_FigS19h.csv": (4_519, "f6ff513d48f07f81e3bb9684e364d656", "2ee2cd1f2528f3119e335f04d50888721e930a175ecc8dfe40174e71ba967ceb"),
    "FigS7g.csv": (3_973, "c1bf781e3e3edddd14142482373a6809", "84a55e2be26fdd1fdf5b9a4db58486f516333b9f7d00dbff84bc3f5ad1ada7a5"),
    "FigS8a.csv": (962, "12c14778c146f9a8e5289e1cd0c6c62c", "93887a0636f903fc1f7bd7a29f13a0086a2edbcc314e460b28c76fcad5709c80"),
    "FigS8b.csv": (224, "65bc3dee007d3edf7decfb667b8f2b15", "823bf67c38b20fcb472b5f89b597a73f0b788892bbc5b0cd542b50645ad0b039"),
    "README.md": (9_538, "20d3f2827a589d4334678d2f4477f013", "b8af694b7a9d88eaca209df844343964ac77bba65fe35857d889e405f7124a94"),
}
EXPECTED_ROWS = {
    "Fig5ab.csv": 2_944,
    "Fig5cd.csv": 4_701,
    "FigS19bc_1uMFAKi.csv": 111,
    "FigS19d_1uMFAKi.csv": 110,
    "FigS19fg.csv": 112,
    "FigS7e.csv": 821,
    "FigS7f_FigS19h.csv": 111,
    "FigS7g.csv": 95,
    "FigS8a.csv": 30,
    "FigS8b.csv": 3,
}
BIOLOGICAL_GROUP = "independent_culture_identifiers_not_reported"
VELOCITY_UNIT = "author_processed_velocity_source_unit_unspecified"


def _digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"missing CSV header: {path.name}")
        rows = list(reader)
    if len(rows) != EXPECTED_ROWS[path.name]:
        raise ValueError(f"unexpected row count for {path.name}: {len(rows)}")
    return rows


def _finite(value: str, locator: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"non-numeric value at {locator}: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"non-finite value at {locator}: {value!r}")
    return result


def verify_sources(source_dir: Path) -> dict[str, str]:
    metadata_path = source_dir / "zenodo_metadata.json"
    metadata_digest = _digest(metadata_path, "sha256")
    if metadata_digest != METADATA_SHA256:
        raise ValueError(f"metadata SHA-256 mismatch: {metadata_digest}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("id") != 8_087_217 or metadata.get("doi") != ACCESSION:
        raise ValueError("unexpected Dryad/Zenodo record identity")
    if metadata["metadata"]["license"]["id"] != "cc-zero":
        raise ValueError("source record is no longer declared CC0")
    official = {item["key"]: item for item in metadata["files"]}
    digests: dict[str, str] = {}
    for name, (size, md5, sha256) in FILES.items():
        path = source_dir / name
        if not path.is_file() or path.stat().st_size != size:
            raise ValueError(f"source size mismatch: {name}")
        if _digest(path, "md5") != md5 or _digest(path, "sha256") != sha256:
            raise ValueError(f"source digest mismatch: {name}")
        record = official.get(name)
        if record is None or record["size"] != size or record["checksum"] != f"md5:{md5}":
            raise ValueError(f"metadata/file mismatch: {name}")
        digests[name] = sha256
    return digests


def _observation(
    *, digests: dict[str, str], filename: str, locator: str, sample_id: str,
    technical_replicate: str, condition: str, cell_state: str, modality: str,
    observable: str, value: float, unit: str, coordinate_name: str = "",
    coordinate_value: float | None = None, coordinate_unit: str = "",
) -> Observation:
    identity = hashlib.sha256(
        f"{digests[filename]}:{locator}:{observable}".encode()
    ).hexdigest()[:20]
    row = Observation(
        observation_id=f"dryad-d96w58:{identity}",
        dataset_id=DATASET_ID,
        lab_group=LAB_GROUP,
        provider=PROVIDER,
        accession=ACCESSION,
        license_id="cc0-1.0",
        source_sha256=digests[filename],
        source_locator=f"csv:{filename};{locator}",
        sample_id=sample_id,
        biological_replicate=BIOLOGICAL_GROUP,
        technical_replicate=technical_replicate,
        condition=condition,
        cell_type="human dermal fibroblasts (hdF)",
        cell_state=cell_state,
        modality=modality,
        observable=observable,
        value=float(value),
        unit=unit,
        coordinate_name=coordinate_name,
        coordinate_value=coordinate_value,
        coordinate_unit=coordinate_unit,
    )
    row.validate()
    return row


def _field_observations(
    source_dir: Path, digests: dict[str, str], filename: str, substrate: str,
) -> list[Observation]:
    rows = _read_csv(source_dir / filename)
    vectors = np.asarray([
        [_finite(row["vx"], f"{filename}:{index}:vx"),
         _finite(row["vy"], f"{filename}:{index}:vy")]
        for index, row in enumerate(rows, start=2)
    ], dtype=np.float64)
    speed = np.hypot(vectors[:, 0], vectors[:, 1])
    moving = speed > 0
    directions = vectors[moving] / speed[moving, None]
    polar_order = np.hypot(*directions.mean(axis=0))
    nematic_order = np.hypot(
        np.mean(directions[:, 0] ** 2 - directions[:, 1] ** 2),
        np.mean(2 * directions[:, 0] * directions[:, 1]),
    )
    metrics = {
        "tracked_vector_count": (len(vectors), "count"),
        "mean_velocity_x": (vectors[:, 0].mean(), VELOCITY_UNIT),
        "mean_velocity_y": (vectors[:, 1].mean(), VELOCITY_UNIT),
        "mean_speed": (speed.mean(), VELOCITY_UNIT),
        "rms_speed": (np.sqrt(np.mean(speed**2)), VELOCITY_UNIT),
        "speed_p90": (np.quantile(speed, 0.9), VELOCITY_UNIT),
        "velocity_polar_order": (polar_order, "dimensionless"),
        "velocity_nematic_order": (nematic_order, "dimensionless"),
    }
    sample = f"{DATASET_ID}:figure5:{substrate}:single_snapshot"
    return [
        _observation(
            digests=digests, filename=filename,
            locator=f"rows=2:{len(rows)+1};statistic={metric}", sample_id=sample,
            technical_replicate=f"single_author_snapshot:{filename}",
            condition=f"{substrate}_substrate", cell_state="high_density_collective_motion",
            modality="single_cell_tracking_velocity_field_derived",
            observable=metric, value=value, unit=unit,
        )
        for metric, (value, unit) in metrics.items()
    ]


def _time_series(
    source_dir: Path, digests: dict[str, str], filename: str, *, condition: str,
    state: str, sample: str, columns: dict[str, tuple[str, str]], time_column: str,
) -> list[Observation]:
    output: list[Observation] = []
    for row_index, source in enumerate(_read_csv(source_dir / filename), start=2):
        time = _finite(source[time_column], f"{filename}:{row_index}:{time_column}")
        for column, (observable, unit) in columns.items():
            output.append(_observation(
                digests=digests, filename=filename,
                locator=f"row={row_index};column={column}", sample_id=sample,
                technical_replicate=f"author_longitudinal_summary:{filename}",
                condition=condition, cell_state=state,
                modality="live_cell_tracking_derived", observable=observable,
                value=_finite(source[column], f"{filename}:{row_index}:{column}"),
                unit=unit, coordinate_name="time", coordinate_value=time,
                coordinate_unit="h",
            ))
    return output


def build(source_dir: Path) -> list[Observation]:
    digests = verify_sources(source_dir)
    output = (
        _field_observations(source_dir, digests, "Fig5ab.csv", "isotropic")
        + _field_observations(source_dir, digests, "Fig5cd.csv", "nematic")
    )

    for row_index, source in enumerate(_read_csv(source_dir / "FigS7e.csv"), start=2):
        substrate = source["substrate"].strip()
        if substrate not in {"isotropic", "nematic"}:
            raise ValueError(f"unexpected substrate at FigS7e row {row_index}")
        day = int(_finite(source["day after"], f"FigS7e:{row_index}:day after"))
        sample = f"{DATASET_ID}:figureS7e:{substrate}:day-{day}:row-{row_index}"
        common = dict(
            digests=digests, filename="FigS7e.csv", sample_id=sample,
            technical_replicate=f"author_frame_summary:FigS7e.csv:row-{row_index}",
            condition=f"{substrate}_substrate", cell_state=f"day_{day}_density_sweep",
            modality="live_cell_tracking_derived", coordinate_name="day_after_seeding",
            coordinate_value=float(day), coordinate_unit="day",
        )
        output.append(_observation(
            **common, locator=f"row={row_index};column=density",
            observable="cell_number_density", value=_finite(
                source["density"], f"FigS7e:{row_index}:density"
            ), unit="cells_per_mm2",
        ))
        output.append(_observation(
            **common, locator=f"row={row_index};column=vx/vy",
            observable="author_velocity_x_to_y_ratio", value=_finite(
                source["vx/vy"], f"FigS7e:{row_index}:vx/vy"
            ), unit="dimensionless",
        ))

    velocity_columns = {
        "vmean": ("ensemble_mean_speed", VELOCITY_UNIT),
        "vxmean": ("ensemble_mean_absolute_velocity_x", VELOCITY_UNIT),
        "vymean": ("ensemble_mean_absolute_velocity_y", VELOCITY_UNIT),
    }
    faki_sample = f"{DATASET_ID}:nematic:FAKi-1uM:single_unresolved_culture"
    control_sample = f"{DATASET_ID}:nematic:no-FAKi:single_unresolved_culture"
    isotropic_sample = f"{DATASET_ID}:isotropic:single_unresolved_culture"
    output += _time_series(
        source_dir, digests, "FigS19bc_1uMFAKi.csv",
        condition="nematic_substrate_FAK_inhibitor_1_micromolar",
        state="FAK_inhibited_collective_motion", sample=faki_sample,
        columns={"density": ("cell_number_density", "cells_per_mm2"),
                 "S_cs": ("cell_substrate_order_parameter", "dimensionless")},
        time_column="t",
    )
    output += _time_series(
        source_dir, digests, "FigS19d_1uMFAKi.csv",
        condition="nematic_substrate_FAK_inhibitor_1_micromolar",
        state="FAK_inhibited_collective_motion", sample=faki_sample,
        columns=velocity_columns, time_column="Time",
    )
    output += _time_series(
        source_dir, digests, "FigS19fg.csv",
        condition="nematic_substrate_no_FAKi_control",
        state="untreated_collective_motion", sample=control_sample,
        columns={"density": ("cell_number_density", "cells_per_mm2"),
                 "S_cs": ("cell_substrate_order_parameter", "dimensionless")},
        time_column="t",
    )
    output += _time_series(
        source_dir, digests, "FigS7f_FigS19h.csv",
        condition="nematic_substrate_no_FAKi_control",
        state="untreated_collective_motion", sample=control_sample,
        columns=velocity_columns, time_column="Time",
    )
    output += _time_series(
        source_dir, digests, "FigS7g.csv", condition="isotropic_substrate",
        state="untreated_collective_motion", sample=isotropic_sample,
        columns=velocity_columns, time_column="Time",
    )

    vaf_specs = {
        "VAF-nematic high density": "nematic_high_density",
        "VAF-nematic low density": "nematic_low_density",
        "VAF-isotropic high density": "isotropic_high_density",
    }
    for row_index, source in enumerate(_read_csv(source_dir / "FigS8a.csv"), start=2):
        lag = _finite(source["lag time"], f"FigS8a:{row_index}:lag time")
        for column, condition in vaf_specs.items():
            output.append(_observation(
                digests=digests, filename="FigS8a.csv",
                locator=f"row={row_index};column={column}",
                sample_id=f"{DATASET_ID}:figureS8a:{condition}",
                technical_replicate=f"author_VAF_summary:{condition}",
                condition=condition, cell_state="velocity_memory_summary",
                modality="live_cell_tracking_derived",
                observable="velocity_autocorrelation_function",
                value=_finite(source[column], f"FigS8a:{row_index}:{column}"),
                unit="dimensionless", coordinate_name="lag_time",
                coordinate_value=lag, coordinate_unit="h",
            ))

    for row_index, source in enumerate(_read_csv(source_dir / "FigS8b.csv"), start=2):
        condition = f"{source['substrate type']}_{source['density']}_density"
        sample = f"{DATASET_ID}:figureS8b:{condition}:author_summary"
        for column, observable, unit in (
            ("avg density", "author_mean_cell_number_density", "cells_per_mm2"),
            ("density deviation", "author_cell_number_density_deviation", "cells_per_mm2"),
            ("mean velocity", "author_mean_speed", VELOCITY_UNIT),
            ("velocity deviation", "author_speed_deviation", VELOCITY_UNIT),
        ):
            output.append(_observation(
                digests=digests, filename="FigS8b.csv",
                locator=f"row={row_index};column={column}", sample_id=sample,
                technical_replicate=f"author_condition_summary:FigS8b:row-{row_index}",
                condition=condition, cell_state="density_stratified_velocity_summary",
                modality="live_cell_tracking_derived", observable=observable,
                value=_finite(source[column], f"FigS8b:{row_index}:{column}"), unit=unit,
            ))

    if len(output) != 3_154 or len({row.sample_id for row in output}) != 832:
        raise ValueError(
            f"unexpected manifest cardinality: {len(output)} observations, "
            f"{len({row.sample_id for row in output})} samples"
        )
    return output


def build_tensor(source_dir: Path, output: Path) -> dict[str, int | str]:
    digests = verify_sources(source_dir)
    arrays, conditions, source_rows = [], [], []
    for filename, condition in (("Fig5ab.csv", "isotropic"), ("Fig5cd.csv", "nematic")):
        for row_index, row in enumerate(_read_csv(source_dir / filename), start=2):
            values = np.asarray([
                _finite(row[column], f"{filename}:{row_index}:{column}")
                for column in ("vx", "vy", "px", "py")
            ] + [
                math.nan if row["velocity angle"] == "#DIV/0!" else _finite(
                    row["velocity angle"], f"{filename}:{row_index}:velocity angle"
                )
            ], dtype=np.float64)
            undefined_stationary_angle = (
                values[0] == 0 and values[1] == 0 and np.isnan(values[4])
            )
            if not undefined_stationary_angle and not np.isclose(
                math.atan2(values[1], values[0]), values[4], atol=1e-7, rtol=0,
            ):
                raise ValueError(f"velocity angle mismatch at {filename}:{row_index}")
            if np.isnan(values[4]) and not undefined_stationary_angle:
                raise ValueError(f"undefined angle for moving vector at {filename}:{row_index}")
            arrays.append(values)
            conditions.append(condition)
            source_rows.append(f"{filename}:{row_index}")
    matrix = np.stack(arrays)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        vectors=matrix,
        columns=np.asarray(["vx", "vy", "px", "py", "velocity_angle"]),
        condition=np.asarray(conditions),
        source_row=np.asarray(source_rows),
        source_sha256=np.asarray([digests["Fig5ab.csv"], digests["Fig5cd.csv"]]),
        split_group=np.asarray([LAB_GROUP] * len(matrix)),
        biological_group=np.asarray([BIOLOGICAL_GROUP] * len(matrix)),
    )
    return {
        "vectors": len(matrix), "features": matrix.shape[1],
        "split_group": LAB_GROUP, "biological_group_count": 1,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tensor-output", type=Path)
    args = parser.parse_args()
    rows = build(args.source_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    result: dict[str, object] = {
        "observations": len(rows), "samples": len({row.sample_id for row in rows})
    }
    if args.tensor_output is not None:
        result["tensor"] = build_tensor(args.source_dir, args.tensor_output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
