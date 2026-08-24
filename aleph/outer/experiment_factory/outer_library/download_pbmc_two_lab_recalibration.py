#!/usr/bin/env python3
"""Acquire the frozen PBMC calibration and untouched confirmation sources."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import tarfile
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any


GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"
FILES = {
    "GSE164378_RAW.tar": f"{GEO}/GSE164nnn/GSE164378/suppl/GSE164378_RAW.tar",
    "GSE164378_sc.meta.data_3P.csv.gz": f"{GEO}/GSE164nnn/GSE164378/suppl/GSE164378_sc.meta.data_3P.csv.gz",
    "GSE222647_AMD_PBMC_countmatrix.csv.gz": f"{GEO}/GSE222nnn/GSE222647/suppl/GSE222647_AMD_PBMC_countmatrix.csv.gz",
    "GSE222647_AMD_PBMC_metadata.csv.gz": f"{GEO}/GSE222nnn/GSE222647/suppl/GSE222647_AMD_PBMC_metadata.csv.gz"
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size:
        return
    partial = destination.with_suffix(destination.suffix + ".part")
    offset = partial.stat().st_size if partial.exists() else 0
    headers = {"User-Agent": "Mozilla/5.0 (compatible; Project-Aleph/1)"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=180) as response:
        status = getattr(response, "status", None)
        if offset and status != 206:
            partial.unlink()
            return download(url, destination)
        mode = "ab" if offset else "wb"
        with partial.open(mode) as output:
            while block := response.read(4 * 1024 * 1024):
                output.write(block)
            output.flush()
            os.fsync(output.fileno())
    partial.replace(destination)


def download_segmented(
    url: str, destination: Path, workers: int = 4,
    segment_bytes: int = 256 * 1024 * 1024,
) -> None:
    """Resume a sequential prefix and fetch the remainder with exact ranges."""
    if destination.exists() and destination.stat().st_size:
        return
    head = urllib.request.Request(
        url, method="HEAD",
        headers={"User-Agent": "Mozilla/5.0 (compatible; Project-Aleph/1)"},
    )
    with urllib.request.urlopen(head, timeout=180) as response:
        total_size = int(response.headers["Content-Length"])
        if "bytes" not in response.headers.get("Accept-Ranges", ""):
            return download(url, destination)
    prefix = destination.with_suffix(destination.suffix + ".part")
    prefix_size = prefix.stat().st_size if prefix.exists() else 0
    if prefix_size > total_size:
        raise ValueError(f"partial file exceeds remote size: {destination.name}")
    remaining = total_size - prefix_size
    if not remaining:
        prefix.replace(destination)
        return
    ranges = [
        (start, min(start + segment_bytes - 1, total_size - 1))
        for start in range(prefix_size, total_size, segment_bytes)
    ]

    def fetch(item: tuple[int, int]) -> Path:
        start, end = item
        segment = destination.with_name(
            f"{destination.name}.range-{start}-{end}.part"
        )
        expected = end - start + 1
        current = segment.stat().st_size if segment.exists() else 0
        if current > expected:
            segment.unlink()
            current = 0
        attempts = 0
        while current < expected:
            request = urllib.request.Request(
                url,
                headers={
                    "Range": f"bytes={start + current}-{end}",
                    "User-Agent": "Mozilla/5.0 (compatible; Project-Aleph/1)",
                },
            )
            try:
                with urllib.request.urlopen(request, timeout=180) as response:
                    if getattr(response, "status", None) != 206:
                        raise ValueError("server ignored exact byte range")
                    with segment.open("ab") as output:
                        while block := response.read(4 * 1024 * 1024):
                            output.write(block)
                current = segment.stat().st_size
                attempts = 0
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
                attempts += 1
                if attempts >= 8:
                    raise
                time.sleep(min(2 ** attempts, 16))
        if segment.stat().st_size != expected:
            raise ValueError(f"short ranged download {segment.name}")
        return segment

    with ThreadPoolExecutor(max_workers=workers) as executor:
        segments = list(executor.map(fetch, ranges))
    assembling = destination.with_suffix(destination.suffix + ".assembling")
    with assembling.open("wb") as output:
        if prefix_size:
            with prefix.open("rb") as source:
                while block := source.read(4 * 1024 * 1024):
                    output.write(block)
        for segment in segments:
            with segment.open("rb") as source:
                while block := source.read(4 * 1024 * 1024):
                    output.write(block)
        output.flush()
        os.fsync(output.fileno())
    if assembling.stat().st_size != total_size:
        raise ValueError("assembled ranged download has wrong size")
    assembling.replace(destination)
    if prefix.exists():
        prefix.unlink()
    for segment in segments:
        segment.unlink()


def acquire(output_dir: Path, protocol_path: Path) -> dict[str, Any]:
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    if not protocol.get(
        "frozen_before_calibration_or_confirmation_expression_and_cell_annotation_access"
    ):
        raise ValueError("two-lab protocol is not frozen before data access")
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for name, url in FILES.items():
        path = output_dir / name
        download(url, path)
        members = None
        if name.endswith(".gz"):
            with gzip.open(path, "rb") as handle:
                handle.read(1)
        elif name.endswith(".tar"):
            with tarfile.open(path) as archive:
                members = [
                    {"name": member.name, "size_bytes": member.size}
                    for member in archive.getmembers() if member.isfile()
                ]
        records.append({
            "name": name,
            "url": url,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
            "container_members": members,
        })
    return {
        "schema": "aleph.outer_library.pbmc_two_lab_acquisition.v1",
        "protocol_sha256": sha256(protocol_path),
        "protocol_frozen_commit": "994b3eb",
        "files": records,
        "source_roles": {
            "GSE164378": "calibration_only",
            "GSE222647": "single_use_confirmation",
            "GSE132044": "opened_retrospective_diagnostic_only"
        },
        "aleph_authority": "none"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    result = acquire(args.output_dir, args.protocol)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "files": len(result["files"]),
        "bytes": sum(item["size_bytes"] for item in result["files"]),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
