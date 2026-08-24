#!/usr/bin/env python3
"""Decode BBBC048 channels into a compact cell-grouped image tensor."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np


CHANNEL_ORDER = ("Ch3", "Ch4", "Ch6")
SPLIT_NAMES = ("train", "selection", "calibration", "test")


def _split_indices(labels: np.ndarray, sample_ids: list[str]) -> np.ndarray:
    """Deterministic stratification; all channels of one cell share one row."""
    split = np.empty(labels.size, dtype=np.uint8)
    for label in np.unique(labels):
        indices = np.flatnonzero(labels == label)
        ordered = sorted(indices, key=lambda i: hashlib.sha256(sample_ids[i].encode()).digest())
        n = len(ordered)
        train_end = int(np.floor(0.65 * n))
        selection_end = train_end + int(np.floor(0.10 * n))
        calibration_end = selection_end + int(np.floor(0.10 * n))
        split[ordered[:train_end]] = 0
        split[ordered[train_end:selection_end]] = 1
        split[ordered[selection_end:calibration_end]] = 2
        split[ordered[calibration_end:]] = 3
    return split


def build(manifest_path: Path, archive_path: Path, output_path: Path) -> dict[str, object]:
    try:
        from PIL import Image
    except ImportError as error:
        raise RuntimeError("Pillow is required only for the BBBC048 decode step") from error
    rows = [json.loads(line) for line in manifest_path.read_text(encoding="utf-8").splitlines()]
    sample_ids = [row["sample_id"] for row in rows]
    labels = np.asarray([row["author_label_code"] for row in rows], dtype=np.uint8)
    pixels = np.empty((len(rows), 3, 16, 16), dtype=np.uint8)
    with zipfile.ZipFile(archive_path) as archive:
        for index, row in enumerate(rows):
            members = {item["channel_id"]: item["archive_member"] for item in row["channels"]}
            for channel_index, channel in enumerate(CHANNEL_ORDER):
                payload = archive.read(f"CellCycle/{members[channel]}")
                with Image.open(io.BytesIO(payload)) as image:
                    if image.mode != "L" or image.size != (66, 66):
                        raise ValueError(f"BBBC048 image contract changed: {image.mode} {image.size}")
                    reduced = image.resize((16, 16), resample=Image.Resampling.BOX)
                    pixels[index, channel_index] = np.asarray(reduced, dtype=np.uint8)
    split = _split_indices(labels, sample_ids)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        pixels=pixels,
        label=labels,
        split=split,
        sample_id=np.asarray(sample_ids, dtype="U40"),
        channel_id=np.asarray(CHANNEL_ORDER, dtype="U3"),
        split_name=np.asarray(SPLIT_NAMES, dtype="U10"),
    )
    split_class_counts = {
        SPLIT_NAMES[split_id]: {
            str(label): int(np.sum((split == split_id) & (labels == label)))
            for label in range(7)
        }
        for split_id in range(4)
    }
    return {
        "schema": "aleph.outer_library.bbbc048_image_tensor.v1",
        "cell_count": len(rows),
        "tensor_shape": list(pixels.shape),
        "tensor_dtype": str(pixels.dtype),
        "decoded_source_shape": [66, 66],
        "compact_shape": [16, 16],
        "channel_order": list(CHANNEL_ORDER),
        "split_class_counts": split_class_counts,
        "split_group_overlap": 0,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "training_role": "same_provider_cell_cycle_supervision",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    receipt = build(args.manifest, args.archive, args.output)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
