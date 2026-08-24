#!/usr/bin/env python3
"""Build exact-cell BBBC054 microglia annotations and brightfield features."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import struct
import zipfile

import numpy as np

from schema import Observation


DATASET_ID = "bbbc054-v1-LPS-microglia-timecourse"
LAB_GROUP = "caltech-bbe-bbbc054"
ARCHIVE_SHA256 = "9661739a008117dee26ad1879d68888e4d05a06cc2aa31668869f2c241208085"
ANNOTATION_SHA256 = "1f4acdd04c7b6854e7de76816e0dd4fa1e5ec06f562c6a81aa5c98549022f59c"
LABELS = ("amoeboid", "ramified", "round")
FEATURE_NAMES = (
    "patch_mean", "patch_std", "patch_q10", "patch_q50", "patch_q90",
    "inner_mean", "outer_mean", "inner_minus_outer", "gradient_x_abs_mean",
    "gradient_y_abs_mean", "center_row_std", "center_column_std",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _values(payload: bytes, endian: str, kind: int, count: int, value: int) -> tuple[int, ...]:
    sizes = {3: 2, 4: 4}
    formats = {3: "H", 4: "I"}
    if kind not in sizes:
        raise ValueError(f"unsupported TIFF tag type {kind}")
    size = sizes[kind] * count
    raw = struct.pack(endian + "I", value)[:size] if size <= 4 else payload[value:value + size]
    return struct.unpack(endian + formats[kind] * count, raw)


def decode_uncompressed_tiff(payload: bytes) -> np.ndarray:
    if payload[:2] == b"II":
        endian = "<"
        dtype = "<u2"
    elif payload[:2] == b"MM":
        endian = ">"
        dtype = ">u2"
    else:
        raise ValueError("invalid TIFF byte order")
    if struct.unpack_from(endian + "H", payload, 2)[0] != 42:
        raise ValueError("invalid TIFF magic")
    offset = struct.unpack_from(endian + "I", payload, 4)[0]
    entry_count = struct.unpack_from(endian + "H", payload, offset)[0]
    tags: dict[int, tuple[int, int, int]] = {}
    for index in range(entry_count):
        tag, kind, count, value = struct.unpack_from(endian + "HHII", payload, offset + 2 + 12 * index)
        tags[tag] = (kind, count, value)
    width = _values(payload, endian, *tags[256])[0]
    height = _values(payload, endian, *tags[257])[0]
    bits = _values(payload, endian, *tags[258])[0]
    compression = _values(payload, endian, *tags[259])[0]
    samples = _values(payload, endian, *tags[277])[0]
    if (bits, compression, samples) != (16, 1, 1):
        raise ValueError("BBBC054 TIFF encoding contract changed")
    strip_offsets = _values(payload, endian, *tags[273])
    strip_bytes = _values(payload, endian, *tags[279])
    raw = b"".join(payload[start:start + length] for start, length in zip(strip_offsets, strip_bytes))
    image = np.frombuffer(raw, dtype=dtype)
    if image.size != width * height:
        raise ValueError("BBBC054 TIFF strip size mismatch")
    return image.reshape(height, width).astype(np.float32)


def patch_features(image: np.ndarray, x: float, y: float, radius: int = 16) -> np.ndarray:
    cx = int(np.clip(round(x), 0, image.shape[1] - 1))
    cy = int(np.clip(round(y), 0, image.shape[0] - 1))
    padded = np.pad(image, radius, mode="reflect")
    patch = padded[cy:cy + 2 * radius + 1, cx:cx + 2 * radius + 1]
    if patch.shape != (2 * radius + 1, 2 * radius + 1):
        raise ValueError("invalid BBBC054 coordinate")
    scale = float(np.percentile(patch, 95) - np.percentile(patch, 5))
    if scale < 1e-6:
        scale = 1.0
    normalized = (patch - np.median(patch)) / scale
    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    inner = xx * xx + yy * yy <= (radius // 2) ** 2
    outer = xx * xx + yy * yy >= (3 * radius // 4) ** 2
    gx = np.diff(normalized, axis=1)
    gy = np.diff(normalized, axis=0)
    q10, q50, q90 = np.quantile(normalized, (0.1, 0.5, 0.9))
    inner_mean = float(normalized[inner].mean())
    outer_mean = float(normalized[outer].mean())
    return np.asarray([
        normalized.mean(), normalized.std(), q10, q50, q90,
        inner_mean, outer_mean, inner_mean - outer_mean,
        np.abs(gx).mean(), np.abs(gy).mean(),
        normalized[radius, :].std(), normalized[:, radius].std(),
    ], dtype=np.float32)


def build(archive: Path, annotation: Path) -> tuple[list[Observation], dict[str, np.ndarray]]:
    if _sha256(archive) != ARCHIVE_SHA256 or _sha256(annotation) != ANNOTATION_SHA256:
        raise ValueError("BBBC054 official source hash mismatch")
    with annotation.open(encoding="utf-8", newline="") as handle:
        annotations = list(csv.DictReader(handle))
    if len(annotations) != 58_186 or set(row["label"] for row in annotations) != set(LABELS):
        raise ValueError("BBBC054 author annotation contract changed")
    by_field: dict[int, list[tuple[int, dict[str, str]]]] = {}
    for csv_row, row in enumerate(annotations, start=2):
        field = int(float(row["axis-0"]))
        by_field.setdefault(field, []).append((csv_row, row))
    if set(by_field) != set(range(60)):
        raise ValueError("BBBC054 annotation does not cover fields 0..59")

    observations: list[Observation] = []
    X = np.empty((len(annotations), len(FEATURE_NAMES)), dtype=np.float32)
    labels: list[str] = []
    fields: list[int] = []
    sample_ids: list[str] = []
    cursor = 0
    with zipfile.ZipFile(archive) as source:
        if source.testzip() is not None:
            raise ValueError("BBBC054 replicate 1 ZIP CRC failure")
        for field in range(60):
            member = f"Replicate 1/IMG_20x_{field}.tif"
            image = decode_uncompressed_tiff(source.read(member))
            for csv_row, row in by_field[field]:
                label = row["label"]
                sample_id = f"replicate1-fov{field:02d}-cell{int(row['index']):05d}"
                # Despite the web-page prose calling axis-1 x and axis-2 y,
                # their ranges match image rows (1080) and columns (1280),
                # respectively. Preserve both original values in provenance.
                axis_1, axis_2 = float(row["axis-1"]), float(row["axis-2"])
                X[cursor] = patch_features(image, axis_2, axis_1)
                labels.append(label)
                fields.append(field)
                sample_ids.append(sample_id)
                observation = Observation(
                    observation_id=f"bbbc054-cell-{int(row['index']):05d}",
                    dataset_id=DATASET_ID,
                    lab_group=LAB_GROUP,
                    provider="Broad Bioimage Benchmark Collection",
                    accession="BBBC054v1",
                    license_id="cc-by-3.0",
                    source_sha256=ANNOTATION_SHA256,
                    source_locator=(
                        f"csv:Replicate1annotation.csv!row:{csv_row};image:{member};"
                        f"axis-1:{axis_1};axis-2:{axis_2}"
                    ),
                    sample_id=sample_id,
                    biological_replicate="replicate_1",
                    technical_replicate=f"field_{field}",
                    condition=f"LPS_exposure_author_FOV_{field}",
                    cell_type="immortalized_mouse_microglia_IMG",
                    cell_state=f"author_morphophenotype_{label}",
                    modality="brightfield_single_cell_patch_derived",
                    observable="author_microglia_morphophenotype_annotation",
                    value=1.0,
                    unit="categorical_annotation_presence",
                    coordinate_name="author_FOV_index",
                    coordinate_value=float(field),
                    coordinate_unit="index_30_minute_order_offset_unspecified",
                )
                observation.validate()
                observations.append(observation)
                cursor += 1
    tensor = {
        "X": X,
        "feature_names": np.asarray(FEATURE_NAMES, dtype="U40"),
        "label": np.asarray(labels, dtype="U16"),
        "field": np.asarray(fields, dtype=np.int16),
        "sample_id": np.asarray(sample_ids, dtype="U48"),
        "source_sha256": np.asarray([ARCHIVE_SHA256, ANNOTATION_SHA256], dtype="U64"),
        "dataset_id": np.asarray([DATASET_ID], dtype="U64"),
        "lab_group": np.asarray([LAB_GROUP], dtype="U64"),
    }
    return observations, tensor


def _write_jsonl(rows: list[Observation], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicate1-zip", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tensor-output", type=Path, required=True)
    args = parser.parse_args()
    observations, tensor = build(args.replicate1_zip, args.annotation)
    _write_jsonl(observations, args.output)
    np.savez_compressed(args.tensor_output, **tensor)
    print(json.dumps({"cells": len(observations), "features": len(FEATURE_NAMES)}, sort_keys=True))


if __name__ == "__main__":
    main()
