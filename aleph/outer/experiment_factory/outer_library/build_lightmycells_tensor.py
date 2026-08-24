#!/usr/bin/env python3
"""Build a compact spatial-prediction tensor from frozen Light My Cells pairs."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image
import tifffile


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode(path: Path) -> np.ndarray:
    array = np.asarray(tifffile.imread(path)).squeeze()
    if array.ndim != 2:
        raise ValueError(f"expected one 2-D TIFF plane, got {array.shape}: {path}")
    if not np.issubdtype(array.dtype, np.integer):
        raise ValueError(f"expected integer microscopy pixels, got {array.dtype}: {path}")
    return array


def normalize_resize(array: np.ndarray, size: int) -> tuple[np.ndarray, dict[str, float]]:
    values = array.astype(np.float32, copy=False)
    low, high = np.percentile(values, [1.0, 99.5])
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        normalized = np.zeros(values.shape, dtype=np.float32)
    else:
        normalized = np.clip((values - low) / (high - low), 0.0, 1.0)
    resized = np.asarray(
        Image.fromarray(normalized).resize((size, size), Image.Resampling.BOX),
        dtype=np.float32,
    )
    if not np.isfinite(resized).all() or resized.min() < 0.0 or resized.max() > 1.0:
        raise ValueError("normalized spatial resize left the frozen [0, 1] range")
    return resized, {
        "raw_p01": float(low),
        "raw_p995": float(high),
        "normalized_mean": float(resized.mean()),
        "normalized_std": float(resized.std()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-manifest", type=Path, required=True)
    parser.add_argument("--acquisition-receipt", type=Path, required=True)
    parser.add_argument("--image-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--size", type=int, default=64)
    args = parser.parse_args()
    if args.size < 32 or args.size > 256:
        raise ValueError("tensor side length must be between 32 and 256")

    pairs = [
        json.loads(line) for line in args.pair_manifest.read_text(encoding="utf-8").splitlines()
    ]
    receipt = json.loads(args.acquisition_receipt.read_text(encoding="utf-8"))
    if receipt["pair_manifest_sha256"] != sha256(args.pair_manifest):
        raise ValueError("acquisition receipt does not pin supplied pair manifest")
    if receipt["confirmation_file_count"] or receipt["confirmation_pixels_accessed"]:
        raise ValueError("acquisition receipt contains forbidden confirmation content")
    receipt_files = {item["path"]: item for item in receipt["files"]}
    required_paths = {
        file["path"] for pair in pairs for file in (pair["input"], pair["target"])
    }
    if set(receipt_files) != required_paths:
        raise ValueError("acquisition receipt and pair payload paths differ")

    inputs = np.empty((len(pairs), args.size, args.size), dtype=np.float16)
    targets = np.empty_like(inputs)
    input_stats: list[dict[str, float]] = []
    target_stats: list[dict[str, float]] = []
    decoded_cache: dict[str, np.ndarray] = {}

    def read(file: dict[str, Any]) -> np.ndarray:
        relative = file["path"]
        if relative not in decoded_cache:
            path = args.image_root / relative
            expected = receipt_files[relative]
            if path.stat().st_size != expected["size_bytes"]:
                raise ValueError(f"payload size changed after acquisition: {relative}")
            decoded_cache[relative] = decode(path)
        return decoded_cache[relative]

    for index, pair in enumerate(pairs):
        input_row, input_summary = normalize_resize(read(pair["input"]), args.size)
        target_row, target_summary = normalize_resize(read(pair["target"]), args.size)
        inputs[index] = input_row.astype(np.float16)
        targets[index] = target_row.astype(np.float16)
        input_stats.append(input_summary)
        target_stats.append(target_summary)
        # Keep memory bounded when mostly unique payloads dominate the source.
        if len(decoded_cache) > 64:
            decoded_cache.clear()

    roles = np.asarray([pair["role"] for pair in pairs])
    studies = np.asarray([pair["study_id"] for pair in pairs])
    acquisitions = np.asarray([pair["acquisition_set_id"] for pair in pairs], dtype=np.int32)
    modalities = np.asarray([pair["input_modality"] for pair in pairs])
    target_names = np.asarray([pair["target_channel"] for pair in pairs])
    alignment = np.asarray([pair["alignment_method"] for pair in pairs])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        inputs=inputs,
        targets=targets,
        role=roles,
        study_id=studies,
        acquisition_set_id=acquisitions,
        input_modality=modalities,
        target_channel=target_names,
        alignment_method=alignment,
    )
    report = {
        "schema": "aleph.outer_library.lightmycells_spatial_tensor.v1",
        "pair_manifest_sha256": sha256(args.pair_manifest),
        "acquisition_receipt_sha256": sha256(args.acquisition_receipt),
        "tensor_sha256": sha256(args.output),
        "pair_count": len(pairs),
        "shape": list(inputs.shape),
        "dtype": str(inputs.dtype),
        "role_counts": dict(sorted(Counter(roles.tolist()).items())),
        "study_count": len(set(studies.tolist())),
        "acquisition_set_count": len(set(zip(studies.tolist(), acquisitions.tolist()))),
        "input_modality_counts": dict(sorted(Counter(modalities.tolist()).items())),
        "target_channel_counts": dict(sorted(Counter(target_names.tolist()).items())),
        "alignment_method_counts": dict(sorted(Counter(alignment.tolist()).items())),
        "normalization": {
            "per_image_percentiles": [1.0, 99.5],
            "clip_range": [0.0, 1.0],
            "resize": "Pillow_BOX",
            "spatial_side_length": args.size,
            "absolute_intensity_prediction_permitted": False
        },
        "input_normalized_mean": float(inputs.astype(np.float32).mean()),
        "target_normalized_mean": float(targets.astype(np.float32).mean()),
        "input_constant_image_count": sum(item["raw_p995"] <= item["raw_p01"] for item in input_stats),
        "target_constant_image_count": sum(item["raw_p995"] <= item["raw_p01"] for item in target_stats),
        "confirmation_pair_count": 0,
        "confirmation_pixels_accessed": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "pairs": len(pairs), "shape": list(inputs.shape),
        "studies": report["study_count"],
        "acquisitions": report["acquisition_set_count"],
        "constant_inputs": report["input_constant_image_count"],
        "constant_targets": report["target_constant_image_count"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
