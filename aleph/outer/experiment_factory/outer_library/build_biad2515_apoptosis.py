#!/usr/bin/env python3
"""Build a leakage-safe well tensor for S-BIAD2515 apoptosis endpoint prediction."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import duckdb
import numpy as np


TARGET = "Terminal_Cytoplasm_Intensity_UpperQuartileIntensity_AnnexinV"
EXPECTED = {
    "train.parquet": "8a92d3d34980a8f115fac5a1ec48e62a50dc914221b82b06308091373d48eebe",
    "test.parquet": "69fabdce42807fd88ce504f34873d840a65503bc4b7d4cf9a41eb346a1e82777",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_parquet(path: Path) -> tuple[list[str], dict[str, np.ndarray]]:
    connection = duckdb.connect()
    columns = [
        row[0]
        for row in connection.execute(
            "DESCRIBE SELECT * FROM read_parquet(?)", [str(path)]
        ).fetchall()
    ]
    values = connection.execute(
        "SELECT * FROM read_parquet(?)", [str(path)]
    ).fetchnumpy()
    for name, value in values.items():
        if np.ma.isMaskedArray(value) and np.ma.getmaskarray(value).any():
            raise ValueError(f"S-BIAD2515 contains null values in {name}")
    return columns, {name: np.asarray(np.ma.getdata(value)) for name, value in values.items()}


def build(
    source_dir: Path,
    manifest_path: Path,
    tensor_path: Path,
) -> dict[str, object]:
    for name, digest in EXPECTED.items():
        if _sha256(source_dir / name) != digest:
            raise ValueError(f"S-BIAD2515 author profile checksum changed: {name}")
    train_columns, train = _read_parquet(source_dir / "train.parquet")
    test_columns, test = _read_parquet(source_dir / "test.parquet")
    if train_columns != test_columns or len(train_columns) != 2_873:
        raise ValueError("S-BIAD2515 train/test feature schema changed")
    metadata_columns = [name for name in train_columns if name.startswith("Metadata_")]
    terminal_columns = [name for name in train_columns if name.startswith("Terminal_")]
    feature_columns = [
        name for name in train_columns
        if name not in metadata_columns and name not in terminal_columns
    ]
    if (
        len(metadata_columns) != 7
        or len(feature_columns) != 2_336
        or len(terminal_columns) != 530
        or TARGET not in terminal_columns
    ):
        raise ValueError("S-BIAD2515 live/terminal feature partition changed")

    def combined(name: str) -> np.ndarray:
        return np.concatenate((train[name], test[name]))

    X = np.column_stack([combined(name).astype(np.float64) for name in feature_columns])
    y = combined(TARGET).astype(np.float64)
    well = combined("Metadata_Well").astype(str)
    dose = combined("Metadata_dose").astype(np.float64)
    apoptosis = combined("Metadata_apoptosis_ground_truth").astype(str)
    single_cell_count = combined("Metadata_number_of_singlecells").astype(np.int64)
    if X.shape != (30, 2_336) or not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("S-BIAD2515 processed tensor has invalid shape or values")

    # The authors seal one well at every dose as test. Split the remaining two
    # wells deterministically into fit and calibration, again within dose.
    split = np.full(30, 3, dtype=np.uint8)
    for dose_value in sorted(set(dose[:20])):
        indices = np.flatnonzero((np.arange(30) < 20) & (dose == dose_value))
        if len(indices) != 2:
            raise ValueError(f"expected two author-train wells at dose {dose_value}")
        ordered = sorted(
            indices,
            key=lambda index: hashlib.sha256(
                f"biad2515-fit-calibration-v1:{well[index]}".encode()
            ).digest(),
        )
        split[ordered[0]] = 0
        split[ordered[1]] = 2
    if not all(np.sum((split == part) & (dose == value)) == 1 for part in (0, 2, 3) for value in set(dose)):
        raise ValueError("S-BIAD2515 dose-balanced well split failed")

    rows = []
    for index in range(30):
        role = {0: "fit", 2: "calibration", 3: "sealed_author_test"}[int(split[index])]
        sample_id = f"S-BIAD2515:{well[index]}"
        rows.append({
            "schema": "aleph.outer_library.observation.v1",
            "observation_id": f"{sample_id}:terminal_annexinV",
            "dataset_id": "biad2515-hela-staurosporine-apoptosis",
            "lab_group": "Way_Lab_CU_Anschutz_Saguaro_Yokogawa",
            "sample_id": sample_id,
            "split_group": sample_id,
            "modality": "live_multichannel_morphology_to_AnnexinV_endpoint_derived",
            "cell_type": "HeLa_CCL2",
            "cell_state": f"author_{apoptosis[index]}_apoptosis_ground_truth",
            "condition": f"staurosporine_{dose[index]:g}_nM",
            "unit": "author_normalized_AnnexinV_morphology_feature",
            "value": float(y[index]),
            "biological_replicate": str(well[index]),
            "source_locator": f"author_processed_profile:{well[index]}",
            "license_id": "CC0-1.0",
            "training_role": role,
            "external_lab_holdout": False,
            "aleph_authority": "none",
            "may_select_aleph_parameter": False,
        })
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )
    tensor_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        tensor_path,
        X=X.astype(np.float32),
        y=y,
        well=well,
        dose_nM=dose,
        author_apoptosis_ground_truth=apoptosis,
        author_aggregated_single_cell_count=single_cell_count,
        split=split,
        split_name=np.asarray(["fit", "unused", "calibration", "sealed_test"], dtype="U16"),
        feature_name=np.asarray(feature_columns, dtype="U160"),
        target_name=np.asarray(TARGET),
    )
    return {
        "schema": "aleph.outer_library.biad2515_apoptosis_tensor.v1",
        "dataset_id": "biad2515-hela-staurosporine-apoptosis",
        "well_count": 30,
        "dose_count": 10,
        "wells_per_dose": 3,
        "author_aggregated_single_cell_count": int(single_cell_count.sum()),
        "live_feature_count": len(feature_columns),
        "terminal_feature_count_excluded_from_inputs": len(terminal_columns),
        "target": TARGET,
        "split_well_counts": {"fit": 10, "calibration": 10, "sealed_test": 10},
        "dose_present_once_in_each_split": True,
        "cross_split_well_overlap_count": 0,
        "terminal_AnnexinV_features_present_in_X": False,
        "manifest_sha256": _sha256(manifest_path),
        "tensor_sha256": _sha256(tensor_path),
        "raw_image_archive_accession": "S-BIAD2515",
        "raw_images_locally_acquired": False,
        "author_processed_profiles_locally_acquired": True,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "status": "same_provider_well_disjoint_training_only",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.source_dir, args.manifest, args.tensor)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
