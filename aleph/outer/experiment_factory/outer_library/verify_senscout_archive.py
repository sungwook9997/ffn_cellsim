#!/usr/bin/env python3
"""Verify the CC0 Dryad SenSCOUT archive before any training extraction."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile


DRYAD_DOI = "10.5061/dryad.rbnzs7hp8"
DRYAD_DATASET_ID = 153_325
DRYAD_VERSION_ID = 360_722
DRYAD_FILE_ID = 4_037_184
EXPECTED_SIZE = 103_221_799
EXPECTED_SHA256 = "1a4a70bbab2fc93d444bf062cc8df7619ea6b1c57344b411764be1fae88adb50"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_member(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts and "\\" not in name


def verify(archive: Path) -> dict[str, object]:
    actual = {"size_bytes": archive.stat().st_size, "sha256": _sha256(archive)}
    expected = {"size_bytes": EXPECTED_SIZE, "sha256": EXPECTED_SHA256}
    if actual != expected:
        raise ValueError(f"SenSCOUT archive digest mismatch: {actual}")
    members: list[dict[str, object]] = []
    roots: Counter[str] = Counter()
    suffixes: Counter[str] = Counter()
    with zipfile.ZipFile(archive) as handle:
        corrupt = handle.testzip()
        if corrupt is not None:
            raise ValueError(f"SenSCOUT ZIP CRC failure: {corrupt}")
        for item in handle.infolist():
            if not _safe_member(item.filename):
                raise ValueError(f"unsafe SenSCOUT ZIP member: {item.filename!r}")
            path = PurePosixPath(item.filename)
            if path.parts:
                roots[path.parts[0]] += 1
            if not item.is_dir():
                suffixes[path.suffix.lower() or "<none>"] += 1
            members.append({
                "name": item.filename,
                "uncompressed_bytes": item.file_size,
                "compressed_bytes": item.compress_size,
                "crc32": f"{item.CRC:08x}",
                "is_directory": item.is_dir(),
            })
    if not any("readme" in str(item["name"]).lower() for item in members):
        raise ValueError("SenSCOUT archive has no README member")
    return {
        "schema": "aleph.outer_library.senscout_archive_receipt.v1",
        "source": {
            "doi": DRYAD_DOI,
            "dryad_dataset_id": DRYAD_DATASET_ID,
            "dryad_version_id": DRYAD_VERSION_ID,
            "dryad_file_id": DRYAD_FILE_ID,
            "license": "CC0-1.0",
        },
        **actual,
        "zip_crc_verified": True,
        "safe_member_paths_verified": True,
        "member_count": len(members),
        "root_member_counts": dict(sorted(roots.items())),
        "file_suffix_counts": dict(sorted(suffixes.items())),
        "members": members,
        "retrieved_at": "2026-08-06",
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    receipt = verify(args.archive)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: receipt[key] for key in (
        "size_bytes", "sha256", "member_count", "file_suffix_counts"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
