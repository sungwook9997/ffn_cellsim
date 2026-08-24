#!/usr/bin/env python3
"""Acquire and verify the official sci-Plex3 sparse UMI count matrix."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from download_sciplex3_metadata import BASE, _download


NAME = "GSM4150378_sciPlex3_A549_MCF7_K562_screen_UMI.count.matrix.gz"
SIZE = 3_096_224_606


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_gzip_stream(path: Path) -> tuple[int, int]:
    """Return decompressed bytes/lines without Python per-line iteration."""
    decompressed_bytes = 0
    decompressed_lines = 0
    final_byte = b""
    with gzip.open(path, "rb") as handle:
        while block := handle.read(16 * 1024 * 1024):
            decompressed_bytes += len(block)
            decompressed_lines += block.count(b"\n")
            final_byte = block[-1:]
    if final_byte and final_byte != b"\n":
        decompressed_lines += 1
    return decompressed_bytes, decompressed_lines


def acquire(output_dir: Path, *, segments: int, workers: int) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / NAME
    url = f"{BASE}/{NAME}"
    _download(url, path, SIZE, segments=segments, workers=workers)
    decompressed_bytes, decompressed_lines = _verify_gzip_stream(path)
    if decompressed_lines < 1_000_000:
        raise ValueError("sci-Plex3 matrix stream has unexpectedly few rows")
    return {
        "schema": "aleph.outer_library.sciplex3_matrix_receipt.v1",
        "accession": "GSM4150378",
        "file": NAME,
        "official_url": url,
        "size_bytes": SIZE,
        "sha256": _sha256(path),
        "gzip_stream_verified": True,
        "decompressed_line_count": decompressed_lines,
        "decompressed_bytes": decompressed_bytes,
        "retrieved_at": "2026-08-06",
        "license_status": "no_explicit_reuse_license_on_GEO_record",
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--segments", type=int, default=96)
    parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()
    receipt = acquire(args.output_dir, segments=args.segments, workers=args.workers)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: receipt[key] for key in (
        "size_bytes", "decompressed_line_count", "sha256"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
