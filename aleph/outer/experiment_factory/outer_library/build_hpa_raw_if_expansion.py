#!/usr/bin/env python3
"""Build and acquire the frozen HPA raw-IF development expansion."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import urllib.request
from typing import Any

import numpy as np

from train_hpa_raw_if_cnn import map_locations, sha256


def order_key(salt: str, value: str) -> str:
    return hashlib.sha256(f"{salt}|{value}".encode()).hexdigest()


def row_from_arrays(data: dict[str, np.ndarray], index: int) -> dict[str, Any]:
    prefix = str(data["file_prefix"][index])
    image_id = hashlib.sha256(prefix.encode()).hexdigest()[:20]
    return {
        "image_id": image_id,
        "file_prefix": prefix,
        "cell_line": str(data["cell_line"][index]),
        "antibody": str(data["antibody"][index]),
        "genes": [item.strip() for item in str(data["genes"][index]).split(",") if item.strip()],
        "author_locations": [
            item.strip() for item in str(data["locations"][index]).split(",") if item.strip()
        ],
        "component_id": str(data["component_id"][index]),
        "split": str(data["split"][index]),
        "channels": {
            channel: {
                "meaning": meaning,
                "url": f"{prefix}_{channel}.jpg",
                "local_name": f"{image_id}_{channel}.jpg",
            }
            for channel, meaning in (("blue", "nucleus"), ("green", "target_protein"))
        },
    }


def build(protocol: dict[str, Any], tensor_path: Path, base_manifest_path: Path) -> dict[str, Any]:
    expansion = protocol["HPA_raw_development_expansion"]
    quotas = expansion["split_image_quotas"]
    base = json.loads(base_manifest_path.read_text(encoding="utf-8"))
    base_by_prefix = {row["file_prefix"]: row for row in base["rows"]}
    output_rows: list[dict[str, Any]] = []
    with np.load(tensor_path, allow_pickle=False) as data:
        arrays = {
            name: data[name].astype(str)
            for name in (
                "file_prefix", "cell_line", "antibody", "genes", "locations",
                "component_id", "split",
            )
        }
        prefixes = arrays["file_prefix"]
        components = arrays["component_id"]
        splits = arrays["split"]
        locations = arrays["locations"]
        used_components: set[str] = set()
        for split_name in ("train", "validation", "calibration", "test"):
            retained = [
                row for row in base["rows"] if row["split"] == split_name
            ]
            retained.sort(key=lambda row: order_key(expansion["selection_salt"], row["file_prefix"]))
            for row in retained:
                used_components.add(row["component_id"])
                compact = dict(row)
                compact["channels"] = {
                    key: value for key, value in row["channels"].items()
                    if key in {"blue", "green"}
                }
                output_rows.append(compact)
            needed = int(quotas[split_name]) - len(retained)
            if needed < 0:
                raise ValueError(f"base HPA manifest exceeds frozen {split_name} quota")
            if split_name == "test":
                if needed:
                    raise ValueError("frozen expansion forbids adding HPA test rows")
                continue
            candidates: dict[int, list[int]] = {
                i: [] for i in range(len(protocol["coarse_ontology_in_order"]))
            }
            eligible_indices = np.flatnonzero(splits == split_name).tolist()
            eligible_indices.sort(
                key=lambda i: order_key(expansion["selection_salt"], prefixes[i])
            )
            for index in eligible_indices:
                if prefixes[index] in base_by_prefix or components[index] in used_components:
                    continue
                labels = [item.strip() for item in locations[index].split(",") if item.strip()]
                mapped = map_locations(labels, protocol)
                for class_index in np.flatnonzero(mapped):
                    candidates[int(class_index)].append(index)
            offsets = {key: 0 for key in candidates}
            selected: list[int] = []
            while len(selected) < needed:
                progressed = False
                for class_index in range(len(candidates)):
                    values = candidates[class_index]
                    while offsets[class_index] < len(values):
                        index = values[offsets[class_index]]
                        offsets[class_index] += 1
                        if components[index] in used_components:
                            continue
                        selected.append(index)
                        used_components.add(components[index])
                        progressed = True
                        break
                    if len(selected) >= needed:
                        break
                if not progressed:
                    break
            if len(selected) != needed:
                raise ValueError(f"could not fill frozen HPA {split_name} expansion quota")
            output_rows.extend(row_from_arrays(arrays, index) for index in selected)
    split_counts = {
        name: sum(row["split"] == name for row in output_rows) for name in quotas
    }
    if split_counts != quotas:
        raise ValueError(f"expanded split counts differ from protocol: {split_counts}")
    component_splits: dict[str, set[str]] = {}
    for row in output_rows:
        component_splits.setdefault(row["component_id"], set()).add(row["split"])
    overlap = sum(len(value) > 1 for value in component_splits.values())
    if overlap:
        raise ValueError("expanded HPA manifest leaks a component across splits")
    output_rows.sort(key=lambda row: (row["split"], row["image_id"]))
    return {
        "schema": "aleph.outer_library.hpa_raw_if_expanded_manifest.v1",
        "dataset_id": "hpa-v25.1-subcellular-raw-if-development-expansion",
        "protocol_sha256": None,
        "source_tensor_sha256": sha256(tensor_path),
        "base_manifest_sha256": sha256(base_manifest_path),
        "selection": expansion["selection"],
        "split_counts": split_counts,
        "component_overlap_across_splits": overlap,
        "image_count": len(output_rows),
        "file_count": len(output_rows) * 2,
        "rows": output_rows,
        "independent_lab_holdout": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def download_one(url: str, path: Path) -> dict[str, Any]:
    if not path.exists():
        error: Exception | None = None
        for attempt in range(5):
            try:
                request = urllib.request.Request(
                    url, headers={"User-Agent": "Project-Aleph-research/1.0"}
                )
                with urllib.request.urlopen(request, timeout=90) as response:
                    payload = response.read()
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(payload)
                break
            except Exception as exc:
                error = exc
                time.sleep(1.5 * (attempt + 1))
        else:
            raise RuntimeError(f"failed to download {url}: {error}")
    if path.stat().st_size < 1000:
        raise ValueError(f"downloaded HPA image is too small: {path.name}")
    return {
        "url": url, "local_name": path.name,
        "size_bytes": path.stat().st_size, "sha256": sha256(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--base-manifest", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol["HPA_raw_development_expansion"][
        "frozen_before_opencell_annotation_or_image_content_access"
    ]:
        raise ValueError("HPA expansion was not frozen before OpenCell access")
    manifest = build(protocol, args.tensor, args.base_manifest)
    manifest["protocol_sha256"] = sha256(args.protocol)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    jobs = [
        (info["url"], args.image_dir / info["local_name"])
        for row in manifest["rows"] for info in row["channels"].values()
    ]
    files = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(download_one, url, path) for url, path in jobs]
        for future in as_completed(futures):
            files.append(future.result())
    files.sort(key=lambda item: item["local_name"])
    receipt = {
        "schema": "aleph.outer_library.hpa_raw_if_expanded_receipt.v1",
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "file_count": len(files),
        "total_bytes": sum(item["size_bytes"] for item in files),
        "files": files,
    }
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "images": manifest["image_count"], "files": len(files),
        "bytes": receipt["total_bytes"], "splits": manifest["split_counts"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
