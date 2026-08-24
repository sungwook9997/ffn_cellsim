#!/usr/bin/env python3
"""Build the compact HPA development tensor under the shared IF contract."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from PIL import Image

from multiprovider_if_common import map_author_labels, normalize_crop_resize
from train_hpa_raw_if_cnn import sha256


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    receipt = json.loads(args.receipt.read_text(encoding="utf-8"))
    if receipt["manifest_sha256"] != sha256(args.manifest):
        raise ValueError("HPA receipt does not pin supplied expanded manifest")
    receipt_files = {item["local_name"]: item for item in receipt["files"]}
    expected = {
        info["local_name"]
        for row in manifest["rows"]
        for info in row["channels"].values()
    }
    if set(receipt_files) != expected:
        raise ValueError("HPA receipt does not exactly cover expanded manifest")

    ontology = protocol["coarse_ontology_in_order"]
    class_names = [item["canonical"] for item in ontology]
    class_lookup = {name: index for index, name in enumerate(class_names)}
    side = int(protocol["image_contract"]["input_shape"][1])
    mapped_rows = [
        (row, map_author_labels(row["author_locations"], ontology))
        for row in manifest["rows"]
    ]
    excluded_rows = [row for row, mapped in mapped_rows if not mapped]
    admitted_rows = [(row, mapped) for row, mapped in mapped_rows if mapped]
    rows = [row for row, _ in admitted_rows]
    images = np.empty((len(rows), 2, side, side), dtype=np.uint8)
    labels = np.zeros((len(rows), len(class_names)), dtype=np.uint8)
    image_ids, components, genes, splits = [], [], [], []
    constant_channel_count = 0
    verified: set[str] = set()

    for index, (row, mapped) in enumerate(admitted_rows):
        for channel_index, channel_name in enumerate(("blue", "green")):
            local_name = row["channels"][channel_name]["local_name"]
            record = receipt_files[local_name]
            path = args.image_dir / local_name
            if path.stat().st_size != record["size_bytes"]:
                raise ValueError(f"HPA payload size changed: {local_name}")
            if local_name not in verified:
                if sha256(path) != record["sha256"]:
                    raise ValueError(f"HPA payload hash changed: {local_name}")
                verified.add(local_name)
            with Image.open(path) as source:
                array = np.asarray(source.convert("F"), dtype=np.float32)
            images[index, channel_index], constant = normalize_crop_resize(array, side)
            constant_channel_count += int(constant)
        for label in mapped:
            labels[index, class_lookup[label]] = 1
        image_ids.append(row["image_id"])
        components.append(row["component_id"])
        genes.append(";".join(row["genes"]))
        split = "selection" if row["split"] == "validation" else row["split"]
        splits.append(split)

    component_splits: dict[str, set[str]] = {}
    for component, split in zip(components, splits):
        component_splits.setdefault(component, set()).add(split)
    overlap = sorted(key for key, value in component_splits.items() if len(value) > 1)
    if overlap:
        raise ValueError("HPA component leaks across multiprovider splits")
    if np.any(labels[np.asarray(splits) == "train"].sum(axis=0) == 0):
        raise ValueError("HPA training split lacks a frozen coarse class")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        images=images,
        labels=labels,
        class_names=np.asarray(class_names),
        image_id=np.asarray(image_ids),
        component_id=np.asarray(components),
        genes=np.asarray(genes),
        split=np.asarray(splits),
    )
    split_array = np.asarray(splits)
    report = {
        "schema": "aleph.outer_library.hpa_multiprovider_tensor.v1",
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "receipt_sha256": sha256(args.receipt),
        "tensor_sha256": sha256(args.output),
        "shape": list(images.shape),
        "dtype": str(images.dtype),
        "encoding": "round(normalized_float_0_1 * 255); model decodes by division by 255",
        "class_names": class_names,
        "image_count": len(rows),
        "source_manifest_image_count": len(manifest["rows"]),
        "excluded_unmapped_image_count": len(excluded_rows),
        "excluded_unmapped_author_location_counts": dict(sorted(Counter(
            location for row in excluded_rows for location in row["author_locations"]
        ).items())),
        "component_count": len(set(components)),
        "split_image_counts": dict(sorted(Counter(splits).items())),
        "split_component_counts": {
            split: len({component for component, assigned in zip(components, splits) if assigned == split})
            for split in sorted(set(splits))
        },
        "positive_image_counts": {
            name: int(labels[:, index].sum()) for index, name in enumerate(class_names)
        },
        "component_overlap_across_splits": len(overlap),
        "constant_channel_count": constant_channel_count,
        "normalization": protocol["image_contract"]["normalization"],
        "geometry": protocol["image_contract"]["geometry"],
        "channel_order": protocol["image_contract"]["channel_order"],
        "HPA_is_external_confirmation": False,
        "Allen_pixels_accessed": 0,
        "model_predictions_made": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "shape": report["shape"], "components": report["component_count"],
        "splits": report["split_component_counts"],
        "constant_channels": constant_channel_count,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
