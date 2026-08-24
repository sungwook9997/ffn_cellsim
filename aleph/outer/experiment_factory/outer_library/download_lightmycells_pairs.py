#!/usr/bin/env python3
"""Download only frozen non-confirmation S-BIAD1047 spatial-pair payloads."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import urllib.request
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def valid_tiff_prefix(path: Path) -> bool:
    with path.open("rb") as handle:
        prefix = handle.read(4)
    return prefix in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")


def acquire(file: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    path = output_dir / file["path"]
    expected_size = int(file["size_bytes"])
    if not path.exists() or path.stat().st_size != expected_size:
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_suffix(path.suffix + ".part")
        error: Exception | None = None
        for attempt in range(5):
            try:
                request = urllib.request.Request(
                    file["url"], headers={"User-Agent": "Project-Aleph-research/1.0"}
                )
                with urllib.request.urlopen(request, timeout=300) as response, part.open("wb") as handle:
                    while True:
                        chunk = response.read(4 * 1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                if part.stat().st_size != expected_size:
                    raise ValueError(
                        f"source size mismatch {part.stat().st_size} != {expected_size}"
                    )
                part.replace(path)
                break
            except Exception as exc:
                error = exc
                part.unlink(missing_ok=True)
                time.sleep(2.0 * (attempt + 1))
        else:
            raise RuntimeError(f"failed to download {file['path']}: {error}")
    if path.stat().st_size != expected_size:
        raise ValueError(f"local size mismatch after acquisition: {file['path']}")
    if not valid_tiff_prefix(path):
        raise ValueError(f"selected payload is not classic TIFF or BigTIFF: {file['path']}")
    return {
        "path": file["path"],
        "url": file["url"],
        "size_bytes": expected_size,
        "sha256": sha256(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--resolution-protocol", type=Path, required=True)
    parser.add_argument("--pair-manifest", type=Path, required=True)
    parser.add_argument("--pair-report", type=Path, required=True)
    parser.add_argument("--integrity-protocol", type=Path)
    parser.add_argument("--integrity-report", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()

    pairs = [
        json.loads(line) for line in args.pair_manifest.read_text(encoding="utf-8").splitlines()
    ]
    if any(pair["role"] == "single_use_confirmation" for pair in pairs):
        raise ValueError("pair manifest contains forbidden confirmation payload")
    if any(pair["confirmation_pixels_accessed"] for pair in pairs):
        raise ValueError("pair manifest claims prior confirmation access")
    pair_report = json.loads(args.pair_report.read_text(encoding="utf-8"))
    if pair_report["resolution_protocol_sha256"] != sha256(args.resolution_protocol):
        raise ValueError("pair report does not pin supplied resolution protocol")
    integrity_report = None
    if bool(args.integrity_protocol) != bool(args.integrity_report):
        raise ValueError("integrity protocol and report must be supplied together")
    if args.integrity_report is not None:
        integrity_report = json.loads(args.integrity_report.read_text(encoding="utf-8"))
        if integrity_report["output_manifest_sha256"] != sha256(args.pair_manifest):
            raise ValueError("integrity report does not pin supplied pair manifest")
        if integrity_report["integrity_protocol_sha256"] != sha256(args.integrity_protocol):
            raise ValueError("integrity report does not pin supplied integrity protocol")
    elif pair_report["output_sha256"] != sha256(args.pair_manifest):
        raise ValueError("pair report does not pin supplied pair manifest")
    unique = {
        file["url"]: file for pair in pairs for file in (pair["input"], pair["target"])
    }
    expected_count = (
        integrity_report["unique_admitted_payload_count"]
        if integrity_report else pair_report["unique_selected_payload_count"]
    )
    expected_bytes = (
        integrity_report["admitted_payload_bytes"]
        if integrity_report else pair_report["selected_payload_bytes"]
    )
    if len(unique) != expected_count:
        raise ValueError("unique payload count changed after pair resolution")
    if sum(int(file["size_bytes"]) for file in unique.values()) != expected_bytes:
        raise ValueError("selected payload bytes changed after pair resolution")

    records = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(acquire, file, args.output_dir) for file in unique.values()]
        for future in as_completed(futures):
            records.append(future.result())
    records.sort(key=lambda item: item["path"])
    part_files = list(args.output_dir.rglob("*.part"))
    if part_files:
        raise ValueError(f"partial files remain after acquisition: {len(part_files)}")
    receipt = {
        "schema": "aleph.outer_library.lightmycells_spatial_pair_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "resolution_protocol_sha256": sha256(args.resolution_protocol),
        "pair_manifest_sha256": sha256(args.pair_manifest),
        "pair_report_sha256": sha256(args.pair_report),
        "integrity_protocol_sha256": (
            sha256(args.integrity_protocol) if args.integrity_protocol else None
        ),
        "integrity_report_sha256": (
            sha256(args.integrity_report) if args.integrity_report else None
        ),
        "file_count": len(records),
        "total_bytes": sum(item["size_bytes"] for item in records),
        "files": records,
        "confirmation_file_count": 0,
        "confirmation_bytes": 0,
        "confirmation_pixels_accessed": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "files": receipt["file_count"], "bytes": receipt["total_bytes"],
        "confirmation_files": 0,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
