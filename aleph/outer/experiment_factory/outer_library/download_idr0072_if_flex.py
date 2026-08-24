#!/usr/bin/env python3
"""Download and checksum the frozen IDR0072 FLEX image manifest."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import threading
import time
import urllib.error
import urllib.request


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def download(row: dict, output_dir: Path, attempts: int = 5) -> dict:
    target = output_dir / row["local_relative_path"]
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        return {"image_id": row["image_id"], "path": row["local_relative_path"], "bytes": target.stat().st_size, "sha256": sha256(target), "resumed": True}
    request = urllib.request.Request(row["source_url"], headers={"User-Agent": "Project-Aleph-research/1.0"})
    part = target.with_suffix(target.suffix + ".part")
    for attempt in range(attempts):
        try:
            digest = hashlib.sha256()
            size = 0
            with urllib.request.urlopen(request, timeout=300) as response, part.open("wb") as handle:
                while block := response.read(4 * 1024 * 1024):
                    handle.write(block)
                    digest.update(block)
                    size += len(block)
                expected = response.headers.get("Content-Length")
            if expected is not None and size != int(expected):
                raise IOError(f"content-length mismatch for image {row['image_id']}: {size} != {expected}")
            part.replace(target)
            return {"image_id": row["image_id"], "path": row["local_relative_path"], "bytes": size, "sha256": digest.hexdigest(), "resumed": False}
        except (urllib.error.URLError, TimeoutError, OSError):
            part.unlink(missing_ok=True)
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    manifest_payload = args.manifest.read_bytes()
    rows = [json.loads(line) for line in manifest_payload.splitlines() if line]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    receipts = []
    completed = 0
    total_bytes = 0
    lock = threading.Lock()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(download, row, args.output_dir): row for row in rows}
        for future in as_completed(futures):
            receipt = future.result()
            receipts.append(receipt)
            with lock:
                completed += 1
                total_bytes += receipt["bytes"]
                if completed % 100 == 0 or completed == len(rows):
                    print(json.dumps({"completed": completed, "files": len(rows), "GB": round(total_bytes / 1e9, 3)}), flush=True)
    receipts.sort(key=lambda item: int(item["image_id"]))
    report = {
        "schema": "aleph.outer_library.idr0072_if_flex_acquisition.v1",
        "status": "complete",
        "source_manifest_sha256": hashlib.sha256(manifest_payload).hexdigest(),
        "file_count": len(receipts),
        "total_bytes": sum(item["bytes"] for item in receipts),
        "resumed_file_count": sum(bool(item["resumed"]) for item in receipts),
        "files": receipts,
        "model_predictions_made_during_acquisition": 0,
        "aleph_parameter_authority": "none"
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
