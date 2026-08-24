#!/usr/bin/env python3
"""Build the frozen single-use Light My Cells confirmation tensor."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any

import numpy as np

from build_lightmycells_tensor import decode, normalize_resize, sha256


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-manifest", type=Path, required=True)
    parser.add_argument("--acquisition-report", type=Path, required=True)
    parser.add_argument("--acquisition-receipt", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--image-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--size", type=int, default=64)
    args = parser.parse_args()
    if args.size < 32 or args.size > 256:
        raise ValueError("tensor side length must be between 32 and 256")

    pairs = [
        json.loads(line)
        for line in args.pair_manifest.read_text(encoding="utf-8").splitlines()
    ]
    report = json.loads(args.acquisition_report.read_text(encoding="utf-8"))
    receipt = json.loads(args.acquisition_receipt.read_text(encoding="utf-8"))
    if report["pair_manifest_sha256"] != sha256(args.pair_manifest):
        raise ValueError("acquisition report does not pin supplied pair manifest")
    if receipt["pair_manifest_sha256"] != sha256(args.pair_manifest):
        raise ValueError("acquisition receipt does not pin supplied pair manifest")
    if receipt["report_sha256"] != sha256(args.acquisition_report):
        raise ValueError("acquisition receipt does not pin supplied report")
    if report["checkpoint_sha256"] != sha256(args.checkpoint):
        raise ValueError("acquisition report does not pin frozen checkpoint")
    if receipt["checkpoint_sha256"] != sha256(args.checkpoint):
        raise ValueError("payload receipt does not pin frozen checkpoint")
    if not report["confirmation_consumed"] or not receipt["confirmation_consumed"]:
        raise ValueError("inputs do not identify the consumed confirmation partition")
    if report["model_predictions_made"] != 0:
        raise ValueError("confirmation acquisition must precede every model prediction")
    if any(pair["role"] != "single_use_confirmation" for pair in pairs):
        raise ValueError("pair manifest contains a non-confirmation row")

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
            if sha256(path) != expected["sha256"]:
                raise ValueError(f"payload digest changed after acquisition: {relative}")
            decoded_cache[relative] = decode(path)
        return decoded_cache[relative]

    for index, pair in enumerate(pairs):
        input_row, input_summary = normalize_resize(read(pair["input"]), args.size)
        target_row, target_summary = normalize_resize(read(pair["target"]), args.size)
        inputs[index] = input_row.astype(np.float16)
        targets[index] = target_row.astype(np.float16)
        input_stats.append(input_summary)
        target_stats.append(target_summary)
        if len(decoded_cache) > 64:
            decoded_cache.clear()

    studies = np.asarray([pair["study_id"] for pair in pairs])
    acquisitions = np.asarray(
        [pair["acquisition_set_id"] for pair in pairs], dtype=np.int32
    )
    modalities = np.asarray([pair["input_modality"] for pair in pairs])
    target_names = np.asarray([pair["target_channel"] for pair in pairs])
    alignment = np.asarray([pair["alignment_method"] for pair in pairs])
    experimenter_groups = np.asarray(
        [pair["experimenter_group"] or "unknown" for pair in pairs]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        inputs=inputs,
        targets=targets,
        role=np.asarray([pair["role"] for pair in pairs]),
        study_id=studies,
        acquisition_set_id=acquisitions,
        input_modality=modalities,
        target_channel=target_names,
        alignment_method=alignment,
        experimenter_group=experimenter_groups,
    )
    tensor_report = {
        "schema": "aleph.outer_library.lightmycells_confirmation_tensor.v1",
        "pair_manifest_sha256": sha256(args.pair_manifest),
        "acquisition_report_sha256": sha256(args.acquisition_report),
        "acquisition_receipt_sha256": sha256(args.acquisition_receipt),
        "checkpoint_sha256": sha256(args.checkpoint),
        "tensor_sha256": sha256(args.output),
        "pair_count": len(pairs),
        "shape": list(inputs.shape),
        "dtype": str(inputs.dtype),
        "study_count": len(set(studies.tolist())),
        "acquisition_set_count": len(set(zip(studies.tolist(), acquisitions.tolist()))),
        "input_modality_counts": dict(sorted(Counter(modalities.tolist()).items())),
        "target_channel_counts": dict(sorted(Counter(target_names.tolist()).items())),
        "alignment_method_counts": dict(sorted(Counter(alignment.tolist()).items())),
        "experimenter_group_counts": dict(
            sorted(Counter(experimenter_groups.tolist()).items())
        ),
        "normalization": {
            "per_image_percentiles": [1.0, 99.5],
            "clip_range": [0.0, 1.0],
            "resize": "Pillow_BOX",
            "spatial_side_length": args.size,
            "absolute_intensity_prediction_permitted": False,
        },
        "input_normalized_mean": float(inputs.astype(np.float32).mean()),
        "target_normalized_mean": float(targets.astype(np.float32).mean()),
        "input_constant_image_count": sum(
            item["raw_p995"] <= item["raw_p01"] for item in input_stats
        ),
        "target_constant_image_count": sum(
            item["raw_p995"] <= item["raw_p01"] for item in target_stats
        ),
        "confirmation_pair_count": len(pairs),
        "confirmation_pixels_accessed": len(pairs),
        "confirmation_consumed": True,
        "may_be_reused_as_pristine": False,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(tensor_report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "pairs": len(pairs),
        "shape": list(inputs.shape),
        "studies": tensor_report["study_count"],
        "acquisitions": tensor_report["acquisition_set_count"],
        "constant_inputs": tensor_report["input_constant_image_count"],
        "constant_targets": tensor_report["target_constant_image_count"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
