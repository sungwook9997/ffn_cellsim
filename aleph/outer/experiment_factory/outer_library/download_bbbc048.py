#!/usr/bin/env python3
"""Download and verify the official BBBC048 Jurkat cell-cycle dataset."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile


BASE = "https://data.broadinstitute.org/bbbc/BBBC048"
ARCHIVE_NAME = "BBBC048v1.zip"
ARCHIVE_SIZE = 1_611_429_608
GROUND_TRUTH_NAME = "Ground_truth.lst"


def _download_range(start: int, end: int, path: Path) -> Path:
    request = urllib.request.Request(
        f"{BASE}/{ARCHIVE_NAME}",
        headers={"Range": f"bytes={start}-{end}", "User-Agent": "Project-Aleph-research/1.0"},
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        payload = response.read()
        content_range = response.headers.get("Content-Range", "")
    if len(payload) != end - start + 1:
        raise ValueError(f"short BBBC048 range {start}-{end}")
    expected = f"bytes {start}-{end}/{ARCHIVE_SIZE}"
    if content_range and content_range != expected:
        raise ValueError(f"unexpected Content-Range {content_range!r}")
    path.write_bytes(payload)
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(output_dir: Path, *, segments: int, workers: int) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    archive_path = output_dir / ARCHIVE_NAME
    if not archive_path.exists() or archive_path.stat().st_size != ARCHIVE_SIZE:
        parts_dir = output_dir / ".parts"
        parts_dir.mkdir(exist_ok=True)
        width = (ARCHIVE_SIZE + segments - 1) // segments
        jobs, parts = [], []
        for index in range(segments):
            start = index * width
            if start >= ARCHIVE_SIZE:
                break
            end = min(ARCHIVE_SIZE - 1, start + width - 1)
            path = parts_dir / f"{ARCHIVE_NAME}.part{index:03d}"
            parts.append(path)
            if not path.exists() or path.stat().st_size != end - start + 1:
                jobs.append((start, end, path))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_download_range, *job) for job in jobs]
            for future in as_completed(futures):
                future.result()
        temporary = output_dir / f".{ARCHIVE_NAME}.assembling"
        with temporary.open("wb") as target:
            for part in parts:
                with part.open("rb") as source:
                    for block in iter(lambda: source.read(1024 * 1024), b""):
                        target.write(block)
        if temporary.stat().st_size != ARCHIVE_SIZE:
            raise ValueError("BBBC048 assembled archive size mismatch")
        temporary.replace(archive_path)
        for part in parts:
            part.unlink(missing_ok=True)
        parts_dir.rmdir()
    with zipfile.ZipFile(archive_path) as archive:
        corrupt = archive.testzip()
        member_count = len(archive.infolist())
    if corrupt:
        raise ValueError(f"BBBC048 ZIP CRC failure: {corrupt}")

    ground_truth_path = output_dir / GROUND_TRUTH_NAME
    request = urllib.request.Request(
        f"{BASE}/{GROUND_TRUTH_NAME}", headers={"User-Agent": "Project-Aleph-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        ground_truth = response.read()
    if len(ground_truth) < 100_000:
        raise ValueError("BBBC048 ground truth unexpectedly small")
    text = ground_truth.decode("utf-8")
    if not all(label in text for label in (
        "G1", "S", "G2", "Prophase", "Metaphase", "Anaphase", "Telophase"
    )):
        raise ValueError("BBBC048 ground-truth phase contract changed")
    ground_truth_path.write_bytes(ground_truth)

    files = [
        {
            "file": ARCHIVE_NAME,
            "official_url": f"{BASE}/{ARCHIVE_NAME}",
            "size_bytes": ARCHIVE_SIZE,
            "sha256": _sha256(archive_path),
            "zip_member_count": member_count,
            "status": "downloaded_and_verified_zip_crc_clean",
        },
        {
            "file": GROUND_TRUTH_NAME,
            "official_url": f"{BASE}/{GROUND_TRUTH_NAME}",
            "size_bytes": len(ground_truth),
            "sha256": hashlib.sha256(ground_truth).hexdigest(),
            "status": "downloaded_and_verified",
        },
    ]
    return {
        "schema": "aleph.outer_library.bbbc048_acquisition_receipt.v1",
        "accession": "BBBC048v1",
        "license_id": "cc-by-nc-sa-3.0",
        "retrieved_at": "2026-08-06",
        "file_count": len(files),
        "total_bytes": sum(int(row["size_bytes"]) for row in files),
        "files": files,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--segments", type=int, default=64)
    parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()
    receipt = download(args.output_dir, segments=args.segments, workers=args.workers)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: receipt[key] for key in ("file_count", "total_bytes")}, sort_keys=True))


if __name__ == "__main__":
    main()
