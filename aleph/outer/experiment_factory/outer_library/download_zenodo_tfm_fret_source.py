#!/usr/bin/env python3
"""Acquire and structurally verify the Zenodo 14692589 TFM--FRET source data."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import urllib.request
import zipfile


RECORD_API = "https://zenodo.org/api/records/14692589"
SOURCE_NAME = "Source-Data.zip"
SOURCE_SIZE = 2_499_140
SOURCE_MD5 = "acf57b9c3e159cc9bb60f4dedfa1f679"
EXPECTED_PRISM = {
    "SD_Fig2b_2f.prism",
    "SD_Fig2c.prism",
    "SD_Fig3d.prism",
    "SD_Fig3e.prism",
    "SD_Fig4d_4e_4f.prism",
    "SD_Fig4g.prism",
    "SD_Fig4i.prism",
    "SD_Fig5e_5f_5g.prism",
    "SD_SI_Fig1.prism",
    "SD_SI_Fig2.prism",
    "SD_SI_Fig3.prism",
    "SD_SI_Fig4a_h.prism",
    "SD_SI_Fig4i_p.prism",
    "SD_SI_Fig7a.prism",
    "SD_SI_Fig7b_7c.prism",
}


def _request(url: str) -> bytes:
    request = urllib.request.Request(
        url, headers={"User-Agent": "Project-Aleph-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read()


def _safe(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts and "\\" not in name


def _verify_archive(payload: bytes) -> dict[str, object]:
    if len(payload) != SOURCE_SIZE:
        raise ValueError(f"Zenodo source size changed: {len(payload)}")
    if hashlib.md5(payload).hexdigest() != SOURCE_MD5:  # noqa: S324 - provider MD5
        raise ValueError("Zenodo provider MD5 mismatch")
    outer_members: list[str] = []
    nested: dict[str, dict[str, object]] = {}
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise ValueError(f"outer ZIP CRC failure: {bad}")
        for info in archive.infolist():
            if not _safe(info.filename):
                raise ValueError(f"unsafe outer ZIP member: {info.filename}")
            if info.is_dir():
                continue
            outer_members.append(info.filename)
            if not info.filename.endswith(".prism"):
                raise ValueError(f"unexpected source-data member: {info.filename}")
            prism_name = PurePosixPath(info.filename).name
            with zipfile.ZipFile(io.BytesIO(archive.read(info))) as prism:
                nested_bad = prism.testzip()
                if nested_bad is not None:
                    raise ValueError(f"nested Prism CRC failure: {nested_bad}")
                members = prism.namelist()
                if any(not _safe(name) for name in members):
                    raise ValueError(f"unsafe nested Prism member: {prism_name}")
                if "document.json" not in members:
                    raise ValueError(f"Prism document metadata missing: {prism_name}")
                data_tables = sum(
                    name.startswith("data/tables/") and name.endswith("/data.csv")
                    for name in members
                )
                nested[prism_name] = {
                    "member_count": len(members),
                    "data_table_count": data_tables,
                    "document_json_present": True,
                    "crc_verified": True,
                }
    observed = {PurePosixPath(name).name for name in outer_members}
    if observed != EXPECTED_PRISM:
        raise ValueError(
            f"Zenodo Prism member contract changed: missing={EXPECTED_PRISM-observed}, "
            f"extra={observed-EXPECTED_PRISM}"
        )
    return {
        "outer_member_count": len(outer_members),
        "outer_crc_verified": True,
        "nested_prism": dict(sorted(nested.items())),
    }


def acquire(output_dir: Path) -> dict[str, object]:
    metadata = json.loads(_request(RECORD_API))
    if metadata.get("metadata", {}).get("license", {}).get("id") != "cc-by-4.0":
        raise ValueError("Zenodo license changed")
    files = {item["key"]: item for item in metadata.get("files", [])}
    source = files.get(SOURCE_NAME)
    if source is None:
        raise ValueError("Zenodo Source-Data.zip is missing")
    if source.get("size") != SOURCE_SIZE or source.get("checksum") != f"md5:{SOURCE_MD5}":
        raise ValueError("Zenodo Source-Data.zip metadata changed")

    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / SOURCE_NAME
    if destination.exists():
        payload = destination.read_bytes()
    else:
        payload = _request(source["links"]["self"])
        destination.write_bytes(payload)
    structure = _verify_archive(payload)
    return {
        "schema": "aleph.outer_library.zenodo_14692589_acquisition.v1",
        "dataset_doi": "10.5281/zenodo.14692589",
        "article_doi": "10.1038/s42003-026-09514-0",
        "title": metadata["metadata"]["title"],
        "license_id": "cc-by-4.0",
        "retrieved_at": "2026-08-06",
        "source_file": {
            "name": SOURCE_NAME,
            "size_bytes": len(payload),
            "provider_md5": SOURCE_MD5,
            "sha256": hashlib.sha256(payload).hexdigest(),
        },
        "structure": structure,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    report = acquire(args.output_dir)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report["source_file"], sort_keys=True))


if __name__ == "__main__":
    main()
