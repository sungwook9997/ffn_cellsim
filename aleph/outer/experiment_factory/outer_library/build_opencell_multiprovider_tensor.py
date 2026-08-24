#!/usr/bin/env python3
"""Build the frozen compact OpenCell multiprovider development tensor."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import tifffile

from multiprovider_if_common import normalize_crop_resize
from train_hpa_raw_if_cnn import sha256


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--resolution-report", type=Path, required=True)
    parser.add_argument("--acquisition-receipt", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    resolution = json.loads(args.resolution_report.read_text(encoding="utf-8"))
    receipt = json.loads(args.acquisition_receipt.read_text(encoding="utf-8"))
    if manifest["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("OpenCell manifest does not pin supplied protocol")
    if resolution["manifest_sha256"] != sha256(args.manifest):
        raise ValueError("OpenCell resolution does not pin supplied manifest")
    if receipt["manifest_sha256"] != sha256(args.manifest):
        raise ValueError("OpenCell acquisition does not pin supplied manifest")
    if receipt["resolution_report_sha256"] != sha256(args.resolution_report):
        raise ValueError("OpenCell acquisition does not pin supplied resolution")
    if receipt["Allen_image_files_accessed"]:
        raise ValueError("OpenCell acquisition unexpectedly contains Allen pixels")
    side = int(protocol["image_contract"]["input_shape"][1])
    class_names = [item["canonical"] for item in protocol["coarse_ontology_in_order"]]
    class_lookup = {name: index for index, name in enumerate(class_names)}
    receipt_files = {item["key"]: item for item in receipt["files"]}
    entries = [
        (row, item) for row in manifest["rows"]
        for item in row["projection_objects"]
    ]
    if len(entries) != receipt["file_count"]:
        raise ValueError("OpenCell tensor entries differ from acquisition count")
    images = np.empty((len(entries), 2, side, side), dtype=np.uint8)
    labels = np.zeros((len(entries), len(class_names)), dtype=np.uint8)
    target_names, splits, roles, projection_keys, overlap = [], [], [], [], []
    constant_channel_count = 0
    verified_paths: set[str] = set()
    for index, (row, object_record) in enumerate(entries):
        acquired = receipt_files[object_record["key"]]
        path = args.image_dir / acquired["local_name"]
        if path.stat().st_size != acquired["size_bytes"]:
            raise ValueError(f"OpenCell payload size changed: {acquired['key']}")
        if acquired["local_name"] not in verified_paths:
            if sha256(path) != acquired["sha256"]:
                raise ValueError(f"OpenCell payload hash changed: {acquired['key']}")
            verified_paths.add(acquired["local_name"])
        array = np.asarray(tifffile.imread(path))
        if array.ndim != 3 or array.shape[0] != 2:
            raise ValueError(f"expected exact C,Y,X two-channel TIFF: {array.shape}")
        if not np.issubdtype(array.dtype, np.integer):
            raise ValueError(f"expected integer OpenCell pixels: {array.dtype}")
        for channel in range(2):
            images[index, channel], constant = normalize_crop_resize(array[channel], side)
            constant_channel_count += int(constant)
        for label in row["mapped_labels"]:
            labels[index, class_lookup[label]] = 1
        target_names.append(row["target_name"])
        splits.append(row["split"])
        roles.append(row["development_role"])
        projection_keys.append(object_record["key"])
        overlap.append(row["exact_HPA_gene_overlap"])
    if not np.isfinite(images).all():
        raise ValueError("OpenCell tensor contains non-finite values")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        images=images,
        labels=labels,
        class_names=np.asarray(class_names),
        target_name=np.asarray(target_names),
        split=np.asarray(splits),
        development_role=np.asarray(roles),
        projection_key=np.asarray(projection_keys),
        exact_HPA_gene_overlap=np.asarray(overlap, dtype=np.uint8),
    )
    report = {
        "schema": "aleph.outer_library.opencell_multiprovider_tensor.v1",
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "resolution_report_sha256": sha256(args.resolution_report),
        "acquisition_receipt_sha256": sha256(args.acquisition_receipt),
        "tensor_sha256": sha256(args.output),
        "shape": list(images.shape),
        "dtype": str(images.dtype),
        "encoding": "round(normalized_float_0_1 * 255); model decodes by division by 255",
        "class_names": class_names,
        "projection_count": len(entries),
        "target_count": len(set(target_names)),
        "split_projection_counts": dict(sorted(Counter(splits).items())),
        "split_target_counts": {
            split: len({target for target, assigned in zip(target_names, splits) if assigned == split})
            for split in sorted(set(splits))
        },
        "role_projection_counts": dict(sorted(Counter(roles).items())),
        "positive_projection_counts": {
            name: int(labels[:, index].sum()) for index, name in enumerate(class_names)
        },
        "constant_channel_count": constant_channel_count,
        "exact_HPA_overlap_projection_count": int(sum(overlap)),
        "normalization": protocol["image_contract"]["normalization"],
        "channel_order": protocol["image_contract"]["channel_order"],
        "OpenCell_is_pristine_confirmation": False,
        "Allen_pixels_accessed": 0,
        "model_predictions_made": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "shape": report["shape"], "targets": report["target_count"],
        "splits": report["split_target_counts"],
        "constant_channels": constant_channel_count,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
