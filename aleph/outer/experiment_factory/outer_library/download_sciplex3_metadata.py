#!/usr/bin/env python3
"""Acquire and verify the official GEO sci-Plex3 annotations and pData."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import hashlib
import json
from pathlib import Path
import time
import urllib.request


BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4150nnn/GSM4150378/suppl"
FILES = {
    "GSM4150378_sciPlex3_A549_MCF7_K562_screen_cell.annotations.txt.gz": 2_770_264,
    "GSM4150378_sciPlex3_A549_MCF7_K562_screen_gene.annotations.txt.gz": 806_168,
    "GSM4150378_sciPlex3_A549_MCF7_K562_hashTable_metadata.txt.gz": 64_715,
    "GSM4150378_sciPlex3_hashSampleSheet.txt.gz": 6_441,
    "GSM4150378_sciPlex3_pData.txt.gz": 84_673_357,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _range(url: str, start: int, end: int, path: Path, total: int) -> None:
    expected = end - start + 1
    partial = path.with_name(path.name + ".partial")
    if partial.exists() and partial.stat().st_size > expected:
        raise ValueError(f"oversized partial range {start}-{end}")
    error: Exception | None = None
    for attempt in range(8):
        received = partial.stat().st_size if partial.exists() else 0
        if received == expected:
            partial.replace(path)
            return
        range_start = start + received
        request = urllib.request.Request(
            url,
            headers={
                "Range": f"bytes={range_start}-{end}",
                "User-Agent": "Project-Aleph-research/1.0",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                content_range = response.headers.get("Content-Range", "")
                expected_range = f"bytes {range_start}-{end}/{total}"
                if content_range != expected_range:
                    raise ValueError(f"unexpected Content-Range {content_range!r}")
                with partial.open("ab") as handle:
                    while block := response.read(1024 * 1024):
                        handle.write(block)
            if partial.stat().st_size == expected:
                partial.replace(path)
                return
            error = ValueError(
                f"short resumable range {start}-{end}: {partial.stat().st_size}/{expected}"
            )
        except Exception as caught:  # transient GEO 429/503 and connection resets
            error = caught
        time.sleep(min(8.0, 0.5 * (2 ** attempt)))
    raise RuntimeError(f"failed resumable range {start}-{end}") from error


def _download(url: str, destination: Path, size: int, *, segments: int, workers: int) -> None:
    if destination.exists() and destination.stat().st_size == size:
        return
    parts_dir = destination.parent / f".{destination.name}.parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    width = (size + segments - 1) // segments
    jobs: list[tuple[str, int, int, Path, int]] = []
    parts: list[Path] = []
    for index in range(segments):
        start = index * width
        if start >= size:
            break
        end = min(size - 1, start + width - 1)
        part = parts_dir / f"part{index:03d}"
        parts.append(part)
        if not part.exists() or part.stat().st_size != end - start + 1:
            jobs.append((url, start, end, part, size))
    pending = jobs
    for round_index in range(30):
        if not pending:
            break
        failed: list[tuple[str, int, int, Path, int]] = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_range, *job): job for job in pending}
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception:
                    failed.append(futures[future])
        pending = failed
        if pending:
            time.sleep(min(15.0, 1.0 + round_index))
    if pending:
        raise RuntimeError(f"failed to acquire {len(pending)} ranges for {destination.name}")
    temporary = destination.parent / f".{destination.name}.assembling"
    with temporary.open("wb") as target:
        for part in parts:
            target.write(part.read_bytes())
    if temporary.stat().st_size != size:
        raise ValueError(f"assembled size mismatch for {destination.name}")
    temporary.replace(destination)
    for part in parts:
        part.unlink()
    parts_dir.rmdir()


def acquire(output_dir: Path, *, segments: int, workers: int) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, size in FILES.items():
        _download(f"{BASE}/{name}", output_dir / name, size, segments=segments, workers=workers)
    for name in FILES:
        with gzip.open(output_dir / name, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
    metadata = output_dir / "GSM4150378_sciPlex3_A549_MCF7_K562_hashTable_metadata.txt.gz"
    with gzip.open(metadata, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 4_992:
        raise ValueError(f"sci-Plex3 condition table changed: {len(rows)}")
    dimensions = {
        "condition_rows": len(rows),
        "cell_types": sorted({row["cell_type"] for row in rows}),
        "replicates": sorted({row["replicate"] for row in rows}),
        "time_points_hours": sorted({int(row["time_point"]) for row in rows}),
        "treatment_count_including_vehicle": len({row["treatment"] for row in rows}),
        "dose_levels": sorted({float(row["dose"]) for row in rows}),
    }
    files = [{
        "file": name, "official_url": f"{BASE}/{name}", "size_bytes": size,
        "sha256": _sha256(output_dir / name), "gzip_stream_verified": True,
    } for name, size in FILES.items()]
    return {
        "schema": "aleph.outer_library.sciplex3_metadata_receipt.v1",
        "accession": "GSM4150378",
        "retrieved_at": "2026-08-06",
        "license_status": "no_explicit_reuse_license_on_GEO_record",
        "files": files,
        "total_bytes": sum(FILES.values()),
        "design_dimensions": dimensions,
        "split_warning": "cells sharing replicate-treatment-dose-time context must not cross splits",
        "external_lab_holdout": False,
        "aleph_authority": "none",
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
    print(json.dumps(receipt["design_dimensions"], sort_keys=True))


if __name__ == "__main__":
    main()
