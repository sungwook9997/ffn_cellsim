#!/usr/bin/env python3
"""Build a field-grouped raw 4-plex IF EMT pilot from the ExIF release."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import struct

import numpy as np

from schema import Observation


DATASET_ID = "figshare-26500210-exif-a549-emt-4plex"
LAB_GROUP = "Lock_Lab_UNSW_Sydney"
PROVIDER = "Springer Nature Figshare"
ACCESSION = "10.6084/m9.figshare.26500210.v1"
PLATE_LAYOUT_SHA256 = "3413f22510385c2376185cacee3433a83e3859059ac41cb5da1e16efe20473da"
CHANNELS = ("TD", "DAPI", "variable_EMT_marker", "B_catenin", "F_actin")
FILE_CHANNEL = {
    "TD": "TD",
    "ex405nm": "DAPI",
    "ex488nm": "variable_EMT_marker",
    "ex555nm": "B_catenin",
    "ex647nm": "F_actin",
}
ROW_MARKER = {
    "A": "PTEN",
    "B": "CD44_total",
    "C": "CD44_standard",
    "D": "CD44_variant_9",
    "E": "E_cadherin",
    "F": "N_cadherin",
    "G": "EpCAM",
    "H": "Vimentin",
}
COLUMN_CONDITION = {
    "01": ("untreated_control", "author_control", 0.0),
    "02": ("untreated_control", "author_control", 0.0),
    "03": ("EGF_0.1ug_per_mL_48h", "author_weak_EMT_induction", 1.0),
    "04": ("EGF_0.1ug_per_mL_48h", "author_weak_EMT_induction", 1.0),
    "05": ("TGF_beta1_10ng_per_mL_48h", "author_strong_EMT_induction", 2.0),
    "06": ("TGF_beta1_10ng_per_mL_48h", "author_strong_EMT_induction", 2.0),
}
FILE_PATTERN = re.compile(
    r"^Images/raw/ExIF_EMT/220825_([A-H])(0[1-6])_(\d{4})_"
    r"(TD|ex405nm|ex488nm|ex555nm|ex647nm)\.tif$"
)
FEATURE_NAMES = tuple(
    f"{channel}_{feature}"
    for channel in CHANNELS
    for feature in (
        "raw_mean", "raw_std", "raw_q10", "raw_q25", "raw_q50",
        "raw_q75", "raw_q90", "normalized_gradient_x_abs_mean",
        "normalized_gradient_y_abs_mean", "normalized_fraction_gt_025",
        "normalized_fraction_gt_050", "normalized_fraction_gt_075",
    )
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _values(payload: bytes, endian: str, kind: int, count: int,
            value: int) -> tuple[int, ...]:
    sizes = {3: 2, 4: 4}
    formats = {3: "H", 4: "I"}
    if kind not in sizes:
        raise ValueError(f"unsupported TIFF tag type {kind}")
    size = sizes[kind] * count
    raw = struct.pack(endian + "I", value)[:size] if size <= 4 else payload[value:value + size]
    return struct.unpack(endian + formats[kind] * count, raw)


def decode_uncompressed_tiff(payload: bytes) -> np.ndarray:
    if payload[:2] == b"II":
        endian, dtype = "<", "<u2"
    elif payload[:2] == b"MM":
        endian, dtype = ">", ">u2"
    else:
        raise ValueError("invalid TIFF byte order")
    if struct.unpack_from(endian + "H", payload, 2)[0] != 42:
        raise ValueError("invalid TIFF magic")
    offset = struct.unpack_from(endian + "I", payload, 4)[0]
    entry_count = struct.unpack_from(endian + "H", payload, offset)[0]
    tags: dict[int, tuple[int, int, int]] = {}
    for index in range(entry_count):
        tag, kind, count, value = struct.unpack_from(
            endian + "HHII", payload, offset + 2 + 12 * index
        )
        tags[tag] = (kind, count, value)
    width = _values(payload, endian, *tags[256])[0]
    height = _values(payload, endian, *tags[257])[0]
    bits = _values(payload, endian, *tags[258])[0]
    compression = _values(payload, endian, *tags[259])[0]
    samples = _values(payload, endian, *tags[277])[0]
    if (width, height, bits, compression, samples) != (1024, 1024, 16, 1, 1):
        raise ValueError("ExIF TIFF encoding contract changed")
    strip_offsets = _values(payload, endian, *tags[273])
    strip_bytes = _values(payload, endian, *tags[279])
    raw = b"".join(
        payload[start:start + length]
        for start, length in zip(strip_offsets, strip_bytes)
    )
    image = np.frombuffer(raw, dtype=dtype)
    if image.size != width * height:
        raise ValueError("ExIF TIFF strip size mismatch")
    return image.reshape(height, width).astype(np.float32)


def _features(image: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q10, q25, q50, q75, q90 = np.quantile(image, (0.10, 0.25, 0.50, 0.75, 0.90))
    low, high = np.quantile(image, (0.01, 0.995))
    scale = max(1.0, float(high - low))
    normalized = np.clip((image - low) / scale, 0.0, 1.0)
    gx = np.diff(normalized, axis=1)
    gy = np.diff(normalized, axis=0)
    feature = np.asarray([
        image.mean(), image.std(), q10, q25, q50, q75, q90,
        np.abs(gx).mean(), np.abs(gy).mean(),
        (normalized > 0.25).mean(), (normalized > 0.50).mean(),
        (normalized > 0.75).mean(),
    ], dtype=np.float32)
    compact = (
        normalized.reshape(64, 16, 64, 16).mean(axis=(1, 3)) * 255.0
    ).round().astype(np.uint8)
    return feature, compact


def build(raw_dir: Path, receipt_path: Path, manifest_path: Path,
          tensor_path: Path) -> dict[str, object]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") != "aleph.outer_library.figshare_exif_emt_layouts.v1":
        raise ValueError("ExIF acquisition receipt schema changed")
    if receipt["source"]["dataset_doi"] != ACCESSION:
        raise ValueError("ExIF dataset identity changed")
    plate_layouts = [
        item for item in receipt["extracted_layouts"].values()
        if item["archive_member"]
        == "Images/raw/ExIF_EMT/220902_A549_EM_ActExp3_Plate_Layout.xlsx"
    ]
    if len(plate_layouts) != 1 or plate_layouts[0]["sha256"] != PLATE_LAYOUT_SHA256:
        raise ValueError("ExIF plate-layout digest changed")
    members = receipt["raw_image_pilot"]["members"]
    expected_fields = tuple(receipt["raw_image_pilot"]["selection"]["field_ids"])
    grouped: dict[tuple[str, str, str], dict[str, tuple[str, dict[str, object]]]] = {}
    for source_name, source_record in members.items():
        match = FILE_PATTERN.fullmatch(source_name)
        if match is None:
            raise ValueError(f"unexpected ExIF source member {source_name}")
        row, column, field, file_channel = match.groups()
        channel = FILE_CHANNEL[file_channel]
        grouped.setdefault((row, column, field), {})[channel] = (
            source_name, source_record
        )
    expected_group_count = 8 * 6 * len(expected_fields)
    if len(grouped) != expected_group_count:
        raise ValueError("ExIF grouped field count changed")

    observations: list[Observation] = []
    feature_rows, compact_rows = [], []
    sample_ids, split_groups, markers, conditions, labels = [], [], [], [], []
    source_sets = []
    for row, column, field in sorted(grouped):
        channel_records = grouped[(row, column, field)]
        if set(channel_records) != set(CHANNELS):
            raise ValueError(f"ExIF field misses channels: {row}{column}/{field}")
        field_features, field_compact, field_hashes, locators = [], [], [], []
        for channel in CHANNELS:
            source_name, source_record = channel_records[channel]
            source_path = raw_dir / Path(source_name).name
            if not source_path.is_file():
                raise ValueError(f"missing ExIF raw member {source_path.name}")
            observed_sha = _sha256(source_path)
            if observed_sha != source_record["sha256"]:
                raise ValueError(f"ExIF raw member digest changed: {source_path.name}")
            image = decode_uncompressed_tiff(source_path.read_bytes())
            feature, compact = _features(image)
            field_features.append(feature)
            field_compact.append(compact)
            field_hashes.append(observed_sha)
            locators.append(source_name)
        marker = ROW_MARKER[row]
        condition, cell_state, label = COLUMN_CONDITION[column]
        well = row + column
        sample_id = f"ExIF-EMT:{well}:field{field}"
        split_group = f"ExIF-EMT:plate3:well{well}"
        source_set_sha = hashlib.sha256(
            "\n".join(field_hashes).encode("ascii")
        ).hexdigest()
        observation = Observation(
            observation_id=f"exif-emt-{well}-field{field}",
            dataset_id=DATASET_ID,
            lab_group=LAB_GROUP,
            provider=PROVIDER,
            accession=ACCESSION,
            license_id="cc-by-4.0",
            source_sha256=source_set_sha,
            source_locator=";".join(locators),
            sample_id=sample_id,
            biological_replicate="not_released_single_plate_date",
            technical_replicate=f"plate3_well_{well}_field_{field}",
            condition=condition,
            cell_type="A549_human_lung_adenocarcinoma",
            cell_state=cell_state,
            modality="raw_five_channel_4plex_IF_field",
            observable="author_EMT_induction_condition",
            value=label,
            unit="categorical_author_treatment_strength",
            unit_source_sha256=PLATE_LAYOUT_SHA256,
            unit_source_locator=(
                "Sheet1!D107:I114,D127:I134,D147:I154,D187:I194,"
                "D207:I214,D227:I234"
            ),
            coordinate_name="author_field_index",
            coordinate_value=float(field),
            coordinate_unit="index",
        )
        observation.validate()
        observations.append(observation)
        feature_rows.append(np.concatenate(field_features))
        compact_rows.append(np.stack(field_compact))
        sample_ids.append(sample_id)
        split_groups.append(split_group)
        markers.append(marker)
        conditions.append(condition)
        labels.append(int(label))
        source_sets.append(source_set_sha)

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        "".join(json.dumps(row.as_dict(), sort_keys=True) + "\n" for row in observations),
        encoding="utf-8",
    )
    np.savez_compressed(
        tensor_path,
        X=np.stack(feature_rows),
        compact_image=np.stack(compact_rows),
        feature_names=np.asarray(FEATURE_NAMES, dtype="U64"),
        channel_names=np.asarray(CHANNELS, dtype="U32"),
        sample_id=np.asarray(sample_ids, dtype="U40"),
        split_group=np.asarray(split_groups, dtype="U40"),
        marker=np.asarray(markers, dtype="U32"),
        condition=np.asarray(conditions, dtype="U40"),
        label=np.asarray(labels, dtype=np.uint8),
        source_set_sha256=np.asarray(source_sets, dtype="U64"),
        dataset_id=np.asarray([DATASET_ID], dtype="U64"),
        lab_group=np.asarray([LAB_GROUP], dtype="U64"),
    )
    return {
        "schema": "aleph.outer_library.figshare_exif_emt_raw_if.v1",
        "source": {
            "dataset_id": DATASET_ID,
            "dataset_doi": ACCESSION,
            "paper_doi": "10.1038/s41467-025-59592-7",
            "license_id": "cc-by-4.0",
        },
        "counts": {
            "field_samples": len(observations),
            "well_split_groups": len(set(split_groups)),
            "raw_channel_files": len(observations) * len(CHANNELS),
            "marker_panels": len(set(markers)),
            "conditions": len(set(conditions)),
        },
        "field_ids": list(expected_fields),
        "source_image_shape": [1024, 1024],
        "compact_tensor_shape": list(np.stack(compact_rows).shape),
        "feature_tensor_shape": list(np.stack(feature_rows).shape),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "tensor_sha256": hashlib.sha256(tensor_path.read_bytes()).hexdigest(),
        "biological_replicates_reported": False,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "training_status": "raw_IF_EMT_condition_training_only",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(
        args.raw_dir, args.receipt, args.manifest, args.tensor
    )
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
