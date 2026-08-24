#!/usr/bin/env python3
"""Acquire only the frozen OpenCell multiprovider development projections."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request
from typing import Any

from train_hpa_raw_if_cnn import sha256


def valid_tiff(path: Path) -> bool:
    with path.open("rb") as handle:
        prefix = handle.read(4)
    return prefix in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")


def acquire(
    record: dict[str, Any], output_dir: Path, provider: dict[str, Any]
) -> dict[str, Any]:
    local_name = hashlib.sha256(record["key"].encode()).hexdigest() + ".tif"
    path = output_dir / local_name
    expected = int(record["size_bytes"])
    quoted = urllib.parse.quote(record["key"], safe="/")
    url = (
        f"https://{provider['aws_bucket']}.s3.{provider['aws_region']}"
        f".amazonaws.com/{quoted}"
    )
    if not path.exists() or path.stat().st_size != expected:
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_suffix(".tif.part")
        error: Exception | None = None
        for attempt in range(5):
            try:
                request = urllib.request.Request(
                    url, headers={"User-Agent": "Project-Aleph-research/1.0"}
                )
                with urllib.request.urlopen(request, timeout=300) as response, part.open(
                    "wb"
                ) as handle:
                    while chunk := response.read(4 * 1024 * 1024):
                        handle.write(chunk)
                if part.stat().st_size != expected:
                    raise ValueError(
                        f"download size {part.stat().st_size} != {expected}"
                    )
                part.replace(path)
                break
            except Exception as exc:
                error = exc
                part.unlink(missing_ok=True)
                time.sleep(2.0 * (attempt + 1))
        else:
            raise RuntimeError(f"failed OpenCell projection {record['key']}: {error}")
    if path.stat().st_size != expected:
        raise ValueError(f"local OpenCell size mismatch: {record['key']}")
    if not valid_tiff(path):
        raise ValueError(f"OpenCell payload is not TIFF: {record['key']}")
    return {
        **record,
        "url": url,
        "local_name": local_name,
        "sha256": sha256(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--resolution-report", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    report = json.loads(args.resolution_report.read_text(encoding="utf-8"))
    if manifest["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("OpenCell manifest does not pin supplied protocol")
    if report["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("OpenCell resolution does not pin supplied protocol")
    if report["manifest_sha256"] != sha256(args.manifest):
        raise ValueError("OpenCell resolution does not pin supplied manifest")
    if not report["all_development_gates_passed"]:
        raise ValueError("OpenCell development gates did not pass")
    unique = {
        item["key"]: item for row in manifest["rows"]
        for item in row["projection_objects"]
    }
    if len(unique) != report["admitted_projection_count"]:
        raise ValueError("OpenCell projection count changed after resolution")
    if sum(item["size_bytes"] for item in unique.values()) != report[
        "admitted_projection_bytes"
    ]:
        raise ValueError("OpenCell projection bytes changed after resolution")
    provider = protocol["provider_roles"]["OpenCell"]
    files = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(acquire, item, args.output_dir, provider)
            for item in unique.values()
        ]
        for future in as_completed(futures):
            files.append(future.result())
    files.sort(key=lambda item: item["key"])
    parts = list(args.output_dir.rglob("*.part"))
    if parts:
        raise ValueError(f"partial OpenCell files remain: {len(parts)}")
    receipt = {
        "schema": "aleph.outer_library.opencell_multiprovider_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "resolution_report_sha256": sha256(args.resolution_report),
        "file_count": len(files),
        "total_bytes": sum(item["size_bytes"] for item in files),
        "files": files,
        "OpenCell_role": "consumed_development_only",
        "Allen_image_files_accessed": 0,
        "model_predictions_made": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "files": receipt["file_count"],
        "bytes": receipt["total_bytes"],
        "Allen_files": 0,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
