#!/usr/bin/env python3
"""Decode frozen native-resolution IDR0168 fields into confirmation tensors."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tarfile
from typing import Any
import urllib.request

import numcodecs
import numpy as np
from scipy import ndimage

from multiprovider_if_common import normalize_crop_resize
from train_hpa_raw_if_cnn import sha256


USER_AGENT = "Project-Aleph-research/1.0"
IDR_TABLE_URL = (
    "https://idr.openmicroscopy.org/searchengine/api/v1/resources/"
    "container_data/?container_name=idr0168-zhang-mllocalization/experimentA"
    "&container_type=project&file_type=csv"
)
IDR_IMAGE_DATA = "https://idr.openmicroscopy.org/webgateway/imgData/{image_id}/"


def fetch(url: str, timeout: int = 180) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def byte_sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def array_sha256(value: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(value)
    return hashlib.sha256(memoryview(contiguous)).hexdigest()


def idr_metadata(rows: list[dict[str, Any]]) -> tuple[bytes, list[dict[str, Any]]]:
    raw = fetch(IDR_TABLE_URL)
    table = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    by_name = {item["name"]: item for item in table}
    if len(table) != 277 or len(by_name) != 277:
        raise ValueError("IDR0168 IDR table row or name count changed")

    def image_data(row: dict[str, Any]) -> dict[str, Any]:
        name = row["image_name"]
        if name not in by_name:
            raise ValueError(f"selected image is absent from IDR table: {name}")
        table_row = by_name[name]
        if table_row["Gene Symbol"] != row["gene_symbol"]:
            raise ValueError(f"IDR table gene mismatch: {name}")
        if table_row["Cell Line"] != row["cell_line"]:
            raise ValueError(f"IDR table cell-line mismatch: {name}")
        image_id = table_row["id"]
        data = json.loads(fetch(IDR_IMAGE_DATA.format(image_id=image_id)))
        data["expected_image_id"] = image_id
        return data

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        metadata = list(pool.map(image_data, rows))
    return raw, metadata


def validate_image_metadata(
    row: dict[str, Any], data: dict[str, Any], gene: str
) -> tuple[int, int]:
    if str(data["id"]) != str(data["expected_image_id"]):
        raise ValueError("IDR image-data response ID changed")
    if data["meta"]["imageName"] != row["image_name"]:
        raise ValueError("IDR image-data response name changed")
    if data["meta"]["pixelsType"] != "uint16":
        raise ValueError("IDR0168 source pixels are not uint16")
    expected = [
        "Alexa555-a-Tubulin",
        "Alexa647-Calreticulin",
        "DAPI-Nucleus",
        f"Alexa488-{gene}",
    ]
    labels = [channel["label"] for channel in data["channels"]]
    if labels != expected:
        raise ValueError(f"IDR0168 named channel contract changed: {labels}")
    size = data["size"]
    if int(size["c"]) != 4 or int(size["t"]) != 1:
        raise ValueError("IDR0168 C/T dimensions changed")
    if int(size["width"]) != 2048 or int(size["height"]) != 2048:
        raise ValueError("IDR0168 native XY dimensions changed")
    return 2, 3


def decode_native_projections(
    path: Path, channel_indices: tuple[int, int]
) -> tuple[np.ndarray, dict[str, Any]]:
    array_name = "./0/0/.zarray"
    chunk_pattern = re.compile(r"^\./0/0/0/(\d+)/(\d+)/(\d+)/(\d+)$")
    metadata: dict[str, Any] | None = None
    codec = None
    projections: dict[int, np.ndarray] = {}
    counts = {index: 0 for index in channel_indices}
    selected_payload_hash = hashlib.sha256()
    with tarfile.open(path, mode="r:gz") as archive:
        for member in archive:
            if member.name == array_name:
                handle = archive.extractfile(member)
                if handle is None:
                    raise ValueError(f"missing readable Zarr metadata: {path.name}")
                metadata = json.loads(handle.read())
                if metadata["zarr_format"] != 2:
                    raise ValueError("unsupported IDR0168 Zarr version")
                if metadata["dimension_separator"] != "/":
                    raise ValueError("unsupported IDR0168 chunk separator")
                if metadata["order"] != "C" or metadata["filters"] is not None:
                    raise ValueError("unsupported IDR0168 Zarr array layout")
                shape = tuple(int(value) for value in metadata["shape"])
                chunks = tuple(int(value) for value in metadata["chunks"])
                if shape[:2] != (1, 4) or shape[3:] != (2048, 2048):
                    raise ValueError(f"IDR0168 native array shape changed: {shape}")
                # The native pyramid level is stored in 1024x1024 tiles.  The
                # 256x256 geometry belongs to a downsampled pyramid level and
                # must not be used for the frozen native-resolution decode.
                if chunks[:3] != (1, 1, 1) or chunks[3:] != (1024, 1024):
                    raise ValueError(f"IDR0168 native chunk shape changed: {chunks}")
                codec = numcodecs.get_codec(metadata["compressor"])
                projections = {
                    index: np.zeros(shape[3:], dtype=np.dtype(metadata["dtype"]))
                    for index in channel_indices
                }
                continue
            match = chunk_pattern.match(member.name)
            if match is None or metadata is None or codec is None:
                continue
            channel, z_index, y_index, x_index = map(int, match.groups())
            if channel not in projections:
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ValueError(f"unreadable IDR0168 Zarr chunk: {member.name}")
            payload = handle.read()
            selected_payload_hash.update(member.name.encode("utf-8"))
            selected_payload_hash.update(payload)
            decoded = codec.decode(payload)
            chunk_shape = tuple(int(value) for value in metadata["chunks"])
            array = np.frombuffer(decoded, dtype=np.dtype(metadata["dtype"]))
            if array.size != math.prod(chunk_shape):
                raise ValueError(f"decoded IDR0168 chunk size changed: {member.name}")
            tile = array.reshape(chunk_shape, order="C")[0, 0, 0]
            y0 = y_index * chunk_shape[3]
            x0 = x_index * chunk_shape[4]
            destination = projections[channel][
                y0 : y0 + tile.shape[0], x0 : x0 + tile.shape[1]
            ]
            np.maximum(destination, tile, out=destination)
            counts[channel] += 1
    if metadata is None:
        raise ValueError(f"native Zarr array metadata absent: {path.name}")
    shape = tuple(int(value) for value in metadata["shape"])
    chunks = tuple(int(value) for value in metadata["chunks"])
    expected = shape[2] * math.ceil(shape[3] / chunks[3]) * math.ceil(
        shape[4] / chunks[4]
    )
    if any(value != expected for value in counts.values()):
        raise ValueError(
            f"IDR0168 selected channel chunks incomplete: {counts} != {expected}"
        )
    stacked = np.stack([projections[index] for index in channel_indices])
    # np.stack canonicalizes same-endian inputs to native byte order.  Preserve
    # the uint16 value contract without incorrectly requiring the container's
    # big-endian marker after stacking.
    if (
        stacked.dtype.kind != "u"
        or stacked.dtype.itemsize != 2
        or not np.isfinite(stacked).all()
    ):
        raise ValueError("IDR0168 decoded projections violate uint16 contract")
    return stacked, {
        "native_z_count": shape[2],
        "native_shape": list(shape),
        "native_chunks": list(chunks),
        "selected_chunk_count_per_channel": counts,
        "selected_compressed_payload_sha256": selected_payload_hash.hexdigest(),
        "native_projection_sha256": array_sha256(stacked),
    }


def otsu_threshold(values: np.ndarray) -> int:
    flat = np.asarray(values, dtype=np.uint16).reshape(-1)
    histogram = np.bincount(flat, minlength=65536).astype(np.float64)
    nonzero = np.flatnonzero(histogram)
    if len(nonzero) < 2:
        raise ValueError("IDR0168 DAPI projection has fewer than two intensities")
    indices = np.arange(len(histogram), dtype=np.float64)
    weight_left = np.cumsum(histogram)
    weighted_left = np.cumsum(histogram * indices)
    total_weight = weight_left[-1]
    total_weighted = weighted_left[-1]
    weight_right = total_weight - weight_left
    valid = (weight_left > 0) & (weight_right > 0)
    difference = np.zeros_like(histogram)
    difference[valid] = (
        weighted_left[valid] / weight_left[valid]
        - (total_weighted - weighted_left[valid]) / weight_right[valid]
    )
    between = weight_left * weight_right * difference**2
    return int(np.argmax(between))


def nuclear_fraction(dapi: np.ndarray, target: np.ndarray) -> tuple[float, dict[str, Any]]:
    threshold = otsu_threshold(dapi)
    initial = dapi > threshold
    labelled, component_count = ndimage.label(initial)
    sizes = np.bincount(labelled.reshape(-1))
    keep = sizes >= 64
    keep[0] = False
    mask = keep[labelled]
    background = float(np.percentile(target.astype(np.float64), 1.0))
    signal = np.clip(target.astype(np.float64) - background, 0.0, None)
    denominator = float(signal.sum())
    if denominator <= 0 or not mask.any():
        raise ValueError("IDR0168 field has empty DAPI mask or target denominator")
    value = float(signal[mask].sum() / denominator)
    if not 0.0 <= value <= 1.0 or not np.isfinite(value):
        raise ValueError("IDR0168 field nuclear fraction is invalid")
    return value, {
        "Otsu_threshold": threshold,
        "initial_component_count": int(component_count),
        "retained_component_count": int(np.sum(keep)),
        "retained_mask_fraction": float(mask.mean()),
        "target_first_percentile": background,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    acquisition = json.loads(args.acquisition.read_text(encoding="utf-8"))
    if acquisition["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("IDR0168 acquisition does not pin supplied protocol")
    if acquisition["selected_stack_count"] != 45:
        raise ValueError("IDR0168 acquisition does not contain 45 frozen stacks")
    rows = acquisition["archives"]
    table_raw, metadata = idr_metadata(rows)

    images = np.empty((len(rows), 2, 128, 128), dtype=np.uint8)
    fractions = np.empty(len(rows), dtype=np.float64)
    genes, cell_lines, image_ids, source_names = [], [], [], []
    decode_reports = []
    for index, (row, image_metadata) in enumerate(zip(rows, metadata)):
        path = args.archive_dir / row["local_name"]
        if path.stat().st_size != row["size_bytes"] or sha256(path) != row["sha256"]:
            raise ValueError(f"IDR0168 archive receipt mismatch: {path.name}")
        channel_indices = validate_image_metadata(
            row, image_metadata, row["gene_symbol"]
        )
        projection, decode = decode_native_projections(path, channel_indices)
        fractions[index], observation = nuclear_fraction(projection[0], projection[1])
        for channel in range(2):
            images[index, channel], constant = normalize_crop_resize(
                projection[channel], 128
            )
            if constant:
                raise ValueError(f"IDR0168 constant decoded channel: {path.name}")
        genes.append(row["gene_symbol"])
        cell_lines.append(row["cell_line"])
        image_ids.append(str(image_metadata["id"]))
        source_names.append(row["source_name"])
        decode_reports.append({
            "archive_name": path.name,
            "image_id": str(image_metadata["id"]),
            "gene_symbol": row["gene_symbol"],
            "cell_line": row["cell_line"],
            "named_channel_indices": {
                "DAPI-Nucleus": channel_indices[0],
                f"Alexa488-{row['gene_symbol']}": channel_indices[1],
            },
            **decode,
            "observation": observation,
            "field_level_intranuclear_target_fraction": fractions[index],
        })
        print(json.dumps({
            "decoded": index + 1,
            "total": len(rows),
            "gene": genes[-1],
            "cell_line": cell_lines[-1],
            "nuclear_fraction": fractions[index],
        }, sort_keys=True), flush=True)

    unseen = set(protocol["development_overlap_strata"][
        "exact_gene_unseen_in_frozen_HPA_and_OpenCell_tensors"
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        images=images,
        field_level_intranuclear_target_fraction=fractions,
        gene_symbol=np.asarray(genes),
        cell_line=np.asarray(cell_lines),
        image_id=np.asarray(image_ids),
        source_name=np.asarray(source_names),
        exact_gene_unseen_in_development=np.asarray(
            [gene in unseen for gene in genes], dtype=np.uint8
        ),
    )
    report = {
        "schema": "aleph.outer_library.idr0168_confirmation_tensor.v1",
        "protocol_sha256": sha256(args.protocol),
        "acquisition_sha256": sha256(args.acquisition),
        "IDR_table_csv_url": IDR_TABLE_URL,
        "IDR_table_csv_size_bytes": len(table_raw),
        "IDR_table_csv_sha256": byte_sha256(table_raw),
        "IDR_imgData_metadata_requests": len(rows),
        "tensor_sha256": sha256(args.output),
        "shape": list(images.shape),
        "dtype": str(images.dtype),
        "field_count": len(rows),
        "unseen_gene_field_count": int(sum(gene in unseen for gene in genes)),
        "seen_gene_field_count": int(sum(gene not in unseen for gene in genes)),
        "field_nuclear_fraction_minimum": float(fractions.min()),
        "field_nuclear_fraction_median": float(np.median(fractions)),
        "field_nuclear_fraction_maximum": float(fractions.max()),
        "numcodecs_version": numcodecs.__version__,
        "pixel_archives_opened_for_decode": len(rows),
        "model_predictions_made": 0,
        "decode_reports": decode_reports,
        "external_provider_validation_complete": False,
        "Aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "tensor": str(args.output),
        "shape": report["shape"],
        "fraction_range": [
            report["field_nuclear_fraction_minimum"],
            report["field_nuclear_fraction_maximum"],
        ],
        "predictions": 0,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
