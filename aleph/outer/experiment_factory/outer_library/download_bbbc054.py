#!/usr/bin/env python3
"""Download and verify the official BBBC054 LPS microglia archive set."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile


BASE = "https://data.broadinstitute.org/bbbc/BBBC054"
FILES = {
    "Replicate_1.zip": 157_736_365,
    "Replicate_2.zip": 158_003_623,
    "Replicate_3.zip": 173_445_539,
    "Replicate1annotation.csv": 3_195_702,
}


def _range_job(url: str, start: int, end: int, path: Path) -> Path:
    request = urllib.request.Request(
        url,
        headers={
            "Range": f"bytes={start}-{end}",
            "User-Agent": "Project-Aleph-research/1.0",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
        content_range = response.headers.get("Content-Range", "")
    if len(payload) != end - start + 1:
        raise ValueError(f"short range {start}-{end} from {url}")
    if content_range and content_range != f"bytes {start}-{end}/{FILES[Path(url).name]}":
        raise ValueError(f"unexpected Content-Range: {content_range}")
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
    parts_dir = output_dir / ".parts"
    parts_dir.mkdir(exist_ok=True)
    jobs: list[tuple[str, int, int, Path]] = []
    by_file: dict[str, list[Path]] = {}
    for name, size in FILES.items():
        url = f"{BASE}/{name}"
        width = (size + segments - 1) // segments
        paths: list[Path] = []
        for index in range(segments):
            start = index * width
            if start >= size:
                break
            end = min(size - 1, start + width - 1)
            path = parts_dir / f"{name}.part{index:03d}"
            paths.append(path)
            if not path.exists() or path.stat().st_size != end - start + 1:
                jobs.append((url, start, end, path))
        by_file[name] = paths

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_range_job, *job) for job in jobs]
        for future in as_completed(futures):
            future.result()

    receipts: list[dict[str, object]] = []
    for name, size in FILES.items():
        destination = output_dir / name
        temporary = output_dir / f".{name}.assembling"
        with temporary.open("wb") as target:
            for part in by_file[name]:
                with part.open("rb") as source:
                    for block in iter(lambda: source.read(1024 * 1024), b""):
                        target.write(block)
        if temporary.stat().st_size != size:
            raise ValueError(f"assembled size mismatch for {name}")
        temporary.replace(destination)
        if destination.suffix == ".zip":
            with zipfile.ZipFile(destination) as archive:
                corrupt = archive.testzip()
                member_count = len(archive.infolist())
            if corrupt:
                raise ValueError(f"ZIP CRC failure in {name}: {corrupt}")
        else:
            member_count = None
            if destination.read_bytes()[:20].lower().find(b"index") < 0:
                raise ValueError("BBBC054 annotation header contract changed")
        receipts.append({
            "file": name,
            "official_url": f"{BASE}/{name}",
            "size_bytes": size,
            "sha256": _sha256(destination),
            "zip_member_count": member_count,
            "status": "downloaded_and_verified" if member_count is None else "downloaded_and_verified_zip_crc_clean",
        })
    for paths in by_file.values():
        for path in paths:
            path.unlink(missing_ok=True)
    parts_dir.rmdir()
    return {
        "schema": "aleph.outer_library.bbbc054_acquisition_receipt.v1",
        "accession": "BBBC054v1",
        "license_id": "cc-by-3.0",
        "retrieved_at": "2026-08-06",
        "file_count": len(receipts),
        "total_bytes": sum(int(row["size_bytes"]) for row in receipts),
        "files": receipts,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--segments", type=int, default=16)
    parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()
    receipt = download(args.output_dir, segments=args.segments, workers=args.workers)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: receipt[key] for key in ("file_count", "total_bytes")}, sort_keys=True))


if __name__ == "__main__":
    main()
