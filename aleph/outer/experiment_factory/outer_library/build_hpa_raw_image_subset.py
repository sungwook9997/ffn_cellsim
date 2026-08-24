#!/usr/bin/env python3
"""Select and optionally download a diverse, split-preserving HPA IF subset."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import urllib.request

import numpy as np


SCHEMA = "aleph.outer_library.hpa_raw_if_subset.v1"
CHANNELS = ("blue", "red", "green", "yellow")
CHANNEL_MEANING = {
    "blue": "nucleus",
    "red": "microtubules",
    "green": "target_protein",
    "yellow": "endoplasmic_reticulum",
}
SPLIT_QUOTAS = {"train": 320, "validation": 80, "calibration": 56, "test": 56}


def _tokens(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def _order_key(prefix: str) -> str:
    return hashlib.sha256(("aleph-hpa-raw-v1|" + prefix).encode()).hexdigest()


def select(tensor_path: Path) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        prefix = data["file_prefix"].astype(str)
        cell_line = data["cell_line"].astype(str)
        antibody = data["antibody"].astype(str)
        genes = data["genes"].astype(str)
        locations = data["locations"].astype(str)
        component = data["component_id"].astype(str)
        split = data["split"].astype(str)
        source_sha256 = str(data["source_sha256"][0])
    selected: list[int] = []
    for split_name, quota in SPLIT_QUOTAS.items():
        candidates = np.flatnonzero(split == split_name).tolist()
        candidates.sort(key=lambda i: _order_key(prefix[i]))
        by_label: dict[str, list[int]] = {}
        for index in candidates:
            for label in _tokens(locations[index]):
                by_label.setdefault(label, []).append(index)
        chosen: list[int] = []
        used: set[int] = set()
        used_components: set[str] = set()
        # Round-robin author locations deliberately raises rare-label coverage.
        offsets = {label: 0 for label in by_label}
        labels = sorted(by_label)
        while len(chosen) < quota:
            progressed = False
            for label in labels:
                rows = by_label[label]
                while offsets[label] < len(rows):
                    index = rows[offsets[label]]
                    offsets[label] += 1
                    if index in used or component[index] in used_components:
                        continue
                    chosen.append(index)
                    used.add(index)
                    used_components.add(component[index])
                    progressed = True
                    break
                if len(chosen) >= quota:
                    break
            if not progressed:
                break
        for index in candidates:
            if len(chosen) >= quota:
                break
            if index not in used and component[index] not in used_components:
                chosen.append(index)
                used.add(index)
                used_components.add(component[index])
        if len(chosen) != quota:
            raise ValueError(f"could not fill HPA raw subset quota for {split_name}")
        selected.extend(chosen)

    rows: list[dict[str, object]] = []
    for index in selected:
        image_id = hashlib.sha256(prefix[index].encode()).hexdigest()[:20]
        channels = {
            channel: {
                "meaning": CHANNEL_MEANING[channel],
                "url": f"{prefix[index]}_{channel}.jpg",
                "local_name": f"{image_id}_{channel}.jpg",
            }
            for channel in CHANNELS
        }
        rows.append({
            "image_id": image_id,
            "file_prefix": prefix[index],
            "cell_line": cell_line[index],
            "antibody": antibody[index],
            "genes": _tokens(genes[index]),
            "author_locations": _tokens(locations[index]),
            "component_id": component[index],
            "split": split[index],
            "channels": channels,
        })
    return {
        "schema": SCHEMA,
        "dataset_id": "hpa-v25.1-subcellular-raw-if-subset",
        "source_embedding_archive_sha256": source_sha256,
        "selection": "deterministic location-round-robin, at most one image per gene-antibody component per split",
        "split_quotas": SPLIT_QUOTAS,
        "image_count": len(rows),
        "download_file_count": len(rows) * len(CHANNELS),
        "cell_line_count": len({str(row["cell_line"]) for row in rows}),
        "author_location_label_count": len({label for row in rows for label in row["author_locations"]}),
        "gene_antibody_component_count": len({str(row["component_id"]) for row in rows}),
        "gene_antibody_component_overlap_across_splits": 0,
        "rows": rows,
        "independent_lab_holdout": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def _download_one(item: tuple[str, Path], retries: int = 4) -> dict[str, object]:
    url, path = item
    if path.exists():
        payload = path.read_bytes()
    else:
        error: Exception | None = None
        for attempt in range(retries):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "Project-Aleph-research/1.0"})
                with urllib.request.urlopen(request, timeout=60) as response:
                    payload = response.read()
                break
            except Exception as exc:  # pragma: no cover - network failure path
                error = exc
                time.sleep(1.5 * (attempt + 1))
        else:
            raise RuntimeError(f"failed to download {url}: {error}")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    if len(payload) < 1_000 or not payload.startswith(b"\xff\xd8") or not payload.endswith(b"\xff\xd9"):
        raise ValueError(f"invalid JPEG payload: {url}")
    return {
        "url": url,
        "local_name": path.name,
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def download(manifest: dict[str, object], directory: Path, workers: int) -> list[dict[str, object]]:
    jobs: list[tuple[str, Path]] = []
    for row in manifest["rows"]:  # type: ignore[union-attr]
        for channel in CHANNELS:
            info = row["channels"][channel]
            jobs.append((info["url"], directory / info["local_name"]))
    receipts: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_download_one, job) for job in jobs]
        for future in as_completed(futures):
            receipts.append(future.result())
    receipts.sort(key=lambda item: str(item["local_name"]))
    return receipts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--download-dir", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    manifest = select(args.tensor)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.download_dir:
        if not args.receipt:
            raise ValueError("--receipt is required with --download-dir")
        receipts = download(manifest, args.download_dir, args.workers)
        receipt = {
            "schema": "aleph.outer_library.hpa_raw_if_download_receipt.v1",
            "image_count": manifest["image_count"],
            "file_count": len(receipts),
            "total_bytes": sum(int(row["size_bytes"]) for row in receipts),
            "files": receipts,
        }
        args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({key: receipt[key] for key in ("image_count", "file_count", "total_bytes")}, sort_keys=True))
    else:
        print(json.dumps({key: manifest[key] for key in ("image_count", "download_file_count", "cell_line_count", "author_location_label_count")}, sort_keys=True))


if __name__ == "__main__":
    main()
