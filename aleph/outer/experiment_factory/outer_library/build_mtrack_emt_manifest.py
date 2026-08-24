#!/usr/bin/env python3
"""Build a single-trajectory M-TRACK EMT temporal representation tensor."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct
import zlib

import numpy as np

from schema import Observation


DATASET_ID = "github-mtrack-a549-tgfb4ng-xy01"
LAB_GROUP = "Xing_Lab_University_of_Pittsburgh_with_ATCC_Cell_Systems"
ACCESSION = "10.1126/sciadv.aba9319"
FEATURE_NAMES = (
    "mask_area_pixels", "mask_area_fraction", "mask_label_count",
    "mask_boundary_pixels", "mask_bbox_aspect_ratio", "mask_centroid_x_fraction",
    "mask_centroid_y_fraction", "VIM_RFP_masked_mean", "VIM_RFP_masked_std",
    "VIM_RFP_masked_q10", "VIM_RFP_masked_q50", "VIM_RFP_masked_q90",
    "VIM_RFP_background_mean", "VIM_RFP_foreground_background_difference",
    "VIM_RFP_foreground_horizontal_gradient", "VIM_RFP_foreground_vertical_gradient",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _tiff_values(payload: bytes, endian: str, kind: int, count: int,
                 value: int) -> tuple[int, ...]:
    sizes, formats = {3: 2, 4: 4}, {3: "H", 4: "I"}
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
    count = struct.unpack_from(endian + "H", payload, offset)[0]
    tags: dict[int, tuple[int, int, int]] = {}
    for index in range(count):
        tag, kind, length, value = struct.unpack_from(
            endian + "HHII", payload, offset + 2 + 12 * index
        )
        tags[tag] = (kind, length, value)
    width = _tiff_values(payload, endian, *tags[256])[0]
    height = _tiff_values(payload, endian, *tags[257])[0]
    bits = _tiff_values(payload, endian, *tags[258])[0]
    compression = _tiff_values(payload, endian, *tags[259])[0]
    samples = _tiff_values(payload, endian, *tags[277])[0]
    if (bits, compression, samples) != (16, 1, 1):
        raise ValueError("M-TRACK TIFF encoding contract changed")
    offsets = _tiff_values(payload, endian, *tags[273])
    byte_counts = _tiff_values(payload, endian, *tags[279])
    raw = b"".join(payload[start:start + size] for start, size in zip(offsets, byte_counts))
    image = np.frombuffer(raw, dtype=dtype)
    if image.size != width * height:
        raise ValueError("M-TRACK TIFF strip size mismatch")
    return image.reshape(height, width).astype(np.float32)


def _paeth(a: int, b: int, c: int) -> int:
    estimate = a + b - c
    pa, pb, pc = abs(estimate - a), abs(estimate - b), abs(estimate - c)
    return a if pa <= pb and pa <= pc else b if pb <= pc else c


def decode_gray16_png(payload: bytes) -> np.ndarray:
    if payload[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("invalid PNG signature")
    cursor, idat, header = 8, [], None
    while cursor < len(payload):
        length = struct.unpack_from(">I", payload, cursor)[0]
        kind = payload[cursor + 4:cursor + 8]
        data = payload[cursor + 8:cursor + 8 + length]
        if kind == b"IHDR":
            header = struct.unpack(">IIBBBBB", data)
        elif kind == b"IDAT":
            idat.append(data)
        elif kind == b"IEND":
            break
        cursor += 12 + length
    if header is None:
        raise ValueError("PNG misses IHDR")
    width, height, bit_depth, color_type, compression, filtering, interlace = header
    if (bit_depth, color_type, compression, filtering, interlace) != (16, 0, 0, 0, 0):
        raise ValueError("M-TRACK mask PNG encoding contract changed")
    packed = zlib.decompress(b"".join(idat))
    stride, bpp = width * 2, 2
    expected = height * (stride + 1)
    if len(packed) != expected:
        raise ValueError("M-TRACK PNG decompressed size changed")
    rows = np.empty((height, stride), dtype=np.uint8)
    cursor = 0
    for y in range(height):
        mode = packed[cursor]
        raw = packed[cursor + 1:cursor + 1 + stride]
        cursor += stride + 1
        row = rows[y]
        prior = rows[y - 1] if y else None
        for x, value in enumerate(raw):
            left = int(row[x - bpp]) if x >= bpp else 0
            up = int(prior[x]) if prior is not None else 0
            up_left = int(prior[x - bpp]) if prior is not None and x >= bpp else 0
            if mode == 0:
                prediction = 0
            elif mode == 1:
                prediction = left
            elif mode == 2:
                prediction = up
            elif mode == 3:
                prediction = (left + up) // 2
            elif mode == 4:
                prediction = _paeth(left, up, up_left)
            else:
                raise ValueError(f"unsupported PNG filter {mode}")
            row[x] = (int(value) + prediction) & 0xFF
    return rows.reshape(height, width, 2).astype(np.uint16).dot(
        np.asarray([256, 1], dtype=np.uint16)
    )


def _features(image: np.ndarray, labels: np.ndarray) -> np.ndarray:
    if image.shape != labels.shape:
        raise ValueError("M-TRACK image/mask shape mismatch")
    mask = labels > 0
    area = int(mask.sum())
    if area == 0:
        raise ValueError("empty M-TRACK segmentation mask")
    y, x = np.nonzero(mask)
    horizontal_edges = mask[:, 1:] != mask[:, :-1]
    vertical_edges = mask[1:, :] != mask[:-1, :]
    boundary = int(horizontal_edges.sum() + vertical_edges.sum())
    width = int(x.max() - x.min() + 1)
    height = int(y.max() - y.min() + 1)
    foreground = image[mask].astype(np.float64)
    background = image[~mask].astype(np.float64)
    q10, q50, q90 = np.quantile(foreground, (0.1, 0.5, 0.9))
    horizontal_pairs = mask[:, 1:] & mask[:, :-1]
    vertical_pairs = mask[1:, :] & mask[:-1, :]
    gx = np.abs(image[:, 1:] - image[:, :-1])[horizontal_pairs]
    gy = np.abs(image[1:, :] - image[:-1, :])[vertical_pairs]
    background_mean = float(background.mean()) if background.size else 0.0
    return np.asarray([
        area, area / mask.size, len(np.unique(labels[mask])), boundary,
        width / height, x.mean() / image.shape[1], y.mean() / image.shape[0],
        foreground.mean(), foreground.std(), q10, q50, q90, background_mean,
        foreground.mean() - background_mean,
        float(gx.mean()) if gx.size else 0.0,
        float(gy.mean()) if gy.size else 0.0,
    ], dtype=np.float64)


def build(source_dir: Path, receipt_path: Path, manifest_path: Path,
          tensor_path: Path) -> dict[str, object]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") != "aleph.outer_library.mtrack_emt_acquisition.v1":
        raise ValueError("M-TRACK acquisition schema changed")
    if not receipt["selection"]["all_frames_form_one_split_group"]:
        raise ValueError("M-TRACK temporal split-group contract changed")
    by_time: dict[int, dict[str, dict[str, object]]] = {}
    for record in receipt["files"]:
        by_time.setdefault(int(record["time_index"]), {})[str(record["kind"])] = record
    observations, feature_rows, sample_ids, source_sets = [], [], [], []
    for time_index in receipt["selection"]["time_indices"]:
        pair = by_time[int(time_index)]
        if set(pair) != {"segmentation_mask", "VIM_RFP"}:
            raise ValueError(f"M-TRACK time {time_index} misses a paired file")
        paths = {kind: source_dir / str(record["local_name"]) for kind, record in pair.items()}
        for kind, path in paths.items():
            if not path.is_file() or path.stat().st_size != pair[kind]["size"]:
                raise ValueError(f"M-TRACK local file changed: {path.name}")
            if _sha256(path) != pair[kind]["sha256"]:
                raise ValueError(f"M-TRACK digest changed: {path.name}")
        image = decode_uncompressed_tiff(paths["VIM_RFP"].read_bytes())
        labels = decode_gray16_png(paths["segmentation_mask"].read_bytes())
        feature_rows.append(_features(image, labels))
        sample_id = f"MTRACK:A549:TGFb4ng:xy01:t{int(time_index):03d}"
        source_set = hashlib.sha256(
            (pair["VIM_RFP"]["sha256"] + "\n" + pair["segmentation_mask"]["sha256"]).encode("ascii")
        ).hexdigest()
        observation = Observation(
            observation_id=f"mtrack-a549-tgfb4ng-xy01-t{int(time_index):03d}",
            dataset_id=DATASET_ID,
            lab_group=LAB_GROUP,
            provider="Science Advances paper-linked GitHub/Google Drive release",
            accession=ACCESSION,
            license_id="not_stated_for_image_data",
            source_sha256=source_set,
            source_locator=(
                f"GoogleDrive:{pair['VIM_RFP']['google_drive_file_id']};"
                f"GoogleDrive:{pair['segmentation_mask']['google_drive_file_id']}"
            ),
            sample_id=sample_id,
            biological_replicate="not_released_single_xy01_trajectory",
            technical_replicate=f"xy01_fluorescence_frame_{int(time_index):03d}",
            condition="TGF_beta_4ng_per_mL",
            cell_type="A549_VIM_RFP_human_lung_adenocarcinoma",
            cell_state="author_TGF_beta_induced_EMT_trajectory",
            modality="live_VIM_RFP_single_trajectory_with_segmentation",
            observable="author_released_EMT_trajectory_membership",
            value=1.0,
            unit="categorical_author_treatment_trajectory",
            coordinate_name="author_fluorescence_frame_index",
            coordinate_value=float(time_index),
            coordinate_unit="index_at_10_minute_imaging_interval",
        )
        observation.validate()
        observations.append(observation)
        sample_ids.append(sample_id)
        source_sets.append(source_set)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        "".join(json.dumps(row.as_dict(), sort_keys=True) + "\n" for row in observations),
        encoding="utf-8",
    )
    X = np.stack(feature_rows)
    np.savez_compressed(
        tensor_path, X=X, feature_names=np.asarray(FEATURE_NAMES, dtype="U64"),
        time_index=np.asarray(receipt["selection"]["time_indices"], dtype=np.int16),
        sample_id=np.asarray(sample_ids, dtype="U64"),
        split_group=np.asarray(["MTRACK:A549:TGFb4ng:xy01"] * len(sample_ids), dtype="U40"),
        source_set_sha256=np.asarray(source_sets, dtype="U64"),
        dataset_id=np.asarray([DATASET_ID], dtype="U64"), lab_group=np.asarray([LAB_GROUP], dtype="U64"),
    )
    return {
        "schema": "aleph.outer_library.mtrack_emt_temporal_raw_vimentin.v1",
        "source": {"dataset_id": DATASET_ID, "paper_doi": ACCESSION, "lab_group": LAB_GROUP},
        "counts": {"timepoints": len(observations), "trajectory_split_groups": 1, "biological_replicates": 0},
        "source_image_shape": list(image.shape), "feature_tensor_shape": list(X.shape),
        "feature_names": list(FEATURE_NAMES),
        "manifest_sha256": _sha256(manifest_path), "tensor_sha256": _sha256(tensor_path),
        "all_frames_form_one_split_group": True,
        "untreated_control_released_in_pilot": False,
        "training_status": "temporal_EMT_representation_diagnostic_only",
        "aleph_authority": "none", "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.source_dir, args.receipt, args.manifest, args.tensor)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
