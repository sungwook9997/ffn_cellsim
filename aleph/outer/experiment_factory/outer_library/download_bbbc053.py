#!/usr/bin/env python3
"""Download and verify the official BBBC053 mitochondrial-stress archive."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import urllib.request
import zipfile


NAME = "FCCP-20220127T153841Z-001.zip"
URL = f"https://data.broadinstitute.org/bbbc/BBBC053/{NAME}"
SIZE = 209_327_519


def _part(start: int, end: int, path: Path) -> None:
    request = urllib.request.Request(
        URL, headers={"Range": f"bytes={start}-{end}", "User-Agent": "Project-Aleph-research/1.0"}
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                payload = response.read()
            break
        except Exception:
            if attempt == 4:
                raise
            time.sleep(0.5 * 2 ** attempt)
    if len(payload) != end - start + 1:
        raise ValueError(f"short BBBC053 range {start}-{end}")
    path.write_bytes(payload)


def acquire(output_dir: Path, *, segments: int, workers: int) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / NAME
    if not destination.exists() or destination.stat().st_size != SIZE:
        part_dir = output_dir / ".parts"
        part_dir.mkdir(exist_ok=True)
        width = (SIZE + segments - 1) // segments
        parts, jobs = [], []
        for index in range(segments):
            start = index * width
            if start >= SIZE:
                break
            end = min(SIZE - 1, start + width - 1)
            part = part_dir / f"part{index:03d}"
            parts.append(part)
            if not part.exists() or part.stat().st_size != end - start + 1:
                jobs.append((start, end, part))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for future in as_completed([pool.submit(_part, *job) for job in jobs]):
                future.result()
        temporary = output_dir / f".{NAME}.assembling"
        with temporary.open("wb") as target:
            for part in parts:
                target.write(part.read_bytes())
        if temporary.stat().st_size != SIZE:
            raise ValueError("BBBC053 assembled size mismatch")
        temporary.replace(destination)
    with zipfile.ZipFile(destination) as archive:
        corrupt = archive.testzip()
        members = archive.infolist()
    if corrupt:
        raise ValueError(f"BBBC053 ZIP CRC failure: {corrupt}")
    hasher = hashlib.sha256()
    with destination.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return {
        "schema": "aleph.outer_library.bbbc053_acquisition_receipt.v1",
        "accession": "BBBC053v1",
        "file": NAME,
        "official_url": URL,
        "size_bytes": SIZE,
        "sha256": hasher.hexdigest(),
        "zip_member_count": len(members),
        "zip_crc_verified": True,
        "license_id": "cc-by-nc-sa-3.0",
        "retrieved_at": "2026-08-06",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--segments", type=int, default=16)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    receipt = acquire(args.output_dir, segments=args.segments, workers=args.workers)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
