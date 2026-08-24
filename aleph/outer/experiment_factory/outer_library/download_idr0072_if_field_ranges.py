#!/usr/bin/env python3
"""Range-acquire the exact first-field TIFF bytes from frozen IDR0072 FLEX files."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import logging
from pathlib import Path
import time
import urllib.error
import urllib.request

import tifffile


RANGE_END = 7_999_999
tifffile.logger().setLevel(logging.CRITICAL)


def source_length(url: str, attempts: int = 5) -> int:
    request = urllib.request.Request(
        url, headers={"User-Agent": "Project-Aleph-research/1.0", "Range": "bytes=0-0"}
    )
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                content_range = response.headers.get("Content-Range")
                if response.status == 206 and content_range and "/" in content_range:
                    return int(content_range.rsplit("/", 1)[1])
                return int(response.headers.get("Content-Length", "0"))
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def first_field_end(payload_path: Path) -> int:
    with tifffile.TiffFile(payload_path) as tif:
        pages = [tif.pages[index] for index in (0, 1)]
        if 65200 not in pages[0].tags:
            raise ValueError(f"missing FLEX XML in {payload_path}")
        return max(
            offset + count
            for page in pages
            for offset, count in zip(page.dataoffsets, page.databytecounts, strict=True)
        )


def acquire(row: dict, output_dir: Path, attempts: int = 5) -> dict:
    target = output_dir / row["local_relative_path"]
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        source_size = target.stat().st_size
        end = first_field_end(target)
        if source_size < end:
            target.unlink()
        else:
            payload = target.read_bytes()[:end]
            part = target.with_suffix(target.suffix + ".range.part")
            part.write_bytes(payload)
            part.replace(target)
            original_size = source_length(row["source_url"])
            return {
                "image_id": row["image_id"], "path": row["local_relative_path"],
                "source_total_bytes": original_size, "retained_range": [0, end - 1],
                "retained_bytes": end, "sha256": sha256(payload), "resumed_or_truncated": True,
            }
    request = urllib.request.Request(
        row["source_url"],
        headers={"User-Agent": "Project-Aleph-research/1.0", "Range": f"bytes=0-{RANGE_END}"},
    )
    probe = target.with_suffix(target.suffix + ".range.part")
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                payload = response.read()
                content_range = response.headers.get("Content-Range")
                status = response.status
            if status == 206:
                if not content_range or "/" not in content_range:
                    raise IOError(f"missing Content-Range for {row['image_id']}")
                source_total = int(content_range.rsplit("/", 1)[1])
            elif status == 200:
                source_total = len(payload)
            else:
                raise IOError(f"unexpected HTTP status {status} for {row['image_id']}")
            probe.write_bytes(payload)
            end = first_field_end(probe)
            if end > len(payload):
                raise IOError(f"range too short for first field of {row['image_id']}: {end} > {len(payload)}")
            retained = payload[:end]
            probe.write_bytes(retained)
            probe.replace(target)
            return {
                "image_id": row["image_id"], "path": row["local_relative_path"],
                "source_total_bytes": source_total, "retained_range": [0, end - 1],
                "retained_bytes": end, "sha256": sha256(retained), "resumed_or_truncated": False,
            }
        except (urllib.error.URLError, TimeoutError, OSError, ValueError):
            probe.unlink(missing_ok=True)
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    manifest_payload = args.manifest.read_bytes()
    rows = [json.loads(line) for line in manifest_payload.splitlines() if line]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    receipts = []
    total_retained = 0
    total_source = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(acquire, row, args.output_dir): row for row in rows}
        for completed, future in enumerate(as_completed(futures), start=1):
            receipt = future.result()
            receipts.append(receipt)
            total_retained += receipt["retained_bytes"]
            total_source += receipt["source_total_bytes"]
            if completed % 100 == 0 or completed == len(rows):
                print(json.dumps({"completed": completed, "files": len(rows), "retained_GB": round(total_retained / 1e9, 3), "source_GB_seen": round(total_source / 1e9, 3)}), flush=True)
    receipts.sort(key=lambda item: int(item["image_id"]))
    report = {
        "schema": "aleph.outer_library.idr0072_if_first_field_range_acquisition.v1",
        "status": "complete",
        "source_manifest_sha256": hashlib.sha256(manifest_payload).hexdigest(),
        "file_count": len(receipts),
        "source_total_bytes_sum": sum(item["source_total_bytes"] for item in receipts),
        "retained_bytes_sum": sum(item["retained_bytes"] for item in receipts),
        "bytes_avoided": sum(item["source_total_bytes"] - item["retained_bytes"] for item in receipts),
        "range_probe_end_inclusive": RANGE_END,
        "files": receipts,
        "retained_content": "original TIFF header, FLEX XML, and exact uint16 strips for sublayout Field #1 only",
        "pixel_values_transformed_during_acquisition": False,
        "model_predictions_made_during_acquisition": 0,
        "aleph_parameter_authority": "none"
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
