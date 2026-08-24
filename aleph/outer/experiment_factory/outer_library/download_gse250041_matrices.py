#!/usr/bin/env python3
"""Acquire and verify the official GSE250041 CITE-seq sparse matrices."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from download_sciplex3_metadata import _download


BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250041/suppl"
FILES = {
    "GSE250041_Proliferating_matrix.mtx.gz": 152_191_341,
    "GSE250041_Senescent_matrix.mtx.gz": 68_344_555,
}
EXPECTED_SHAPES = {
    "GSE250041_Proliferating_matrix.mtx.gz": (36_609, 8_664),
    "GSE250041_Senescent_matrix.mtx.gz": (36_609, 4_949),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix_market_header(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="ascii") as handle:
        if handle.readline().strip() != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"unexpected MatrixMarket banner: {path.name}")
        for line in handle:
            if not line.startswith("%"):
                rows, columns, nonzero = (int(value) for value in line.split())
                return rows, columns, nonzero
    raise ValueError(f"missing MatrixMarket dimensions: {path.name}")


def acquire(output_dir: Path, *, segments: int, workers: int) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    verified: dict[str, object] = {}
    for name, size in FILES.items():
        path = output_dir / name
        _download(f"{BASE}/{name}", path, size, segments=segments, workers=workers)
        shape = _matrix_market_header(path)
        if shape[:2] != EXPECTED_SHAPES[name]:
            raise ValueError(f"GSE250041 matrix shape changed for {name}: {shape}")
        with gzip.open(path, "rb") as handle:
            while handle.read(16 * 1024 * 1024):
                pass
        verified[name] = {
            "official_url": f"{BASE}/{name}",
            "size_bytes": size,
            "sha256": _sha256(path),
            "shape": {"features": shape[0], "barcodes": shape[1], "nonzero": shape[2]},
            "gzip_stream_verified": True,
        }
    return {
        "schema": "aleph.outer_library.gse250041_matrix_receipt.v1",
        "accession": "GSE250041",
        "retrieved_at": "2026-08-06",
        "files": verified,
        "feature_axis_includes": {"gene_expression": 36_601, "antibody_capture": 8},
        "single_capture_per_condition": True,
        "external_validation_eligible": False,
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--segments", type=int, default=64)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    receipt = acquire(args.output_dir, segments=args.segments, workers=args.workers)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt["files"], sort_keys=True))


if __name__ == "__main__":
    main()
