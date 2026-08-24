#!/usr/bin/env python3
"""Range-extract verified figure tables from the 4.9 GB Zenodo TFM archive."""

from __future__ import annotations

import argparse
import binascii
import hashlib
import json
from pathlib import Path, PurePosixPath
import struct
import time
import urllib.error
import urllib.request
import zlib


RECORD_API = "https://zenodo.org/api/records/7432971"
ARCHIVE_URL = "https://zenodo.org/api/records/7432971/files/data_zenodo.zip/content"
ARCHIVE_SIZE = 4_942_831_343
ARCHIVE_MD5 = "45fd97f362c83333452fc78e93592a55"
CENTRAL_DIRECTORY = {"offset": 4_942_163_365, "size": 667_880, "member_count": 3_408}
EXPECTED = {
    "data_zenodo/Kreysing_McHugh_key_resources.xlsx": {
        "crc32": 0x36C09F39, "compressed_size": 8_849,
        "uncompressed_size": 9_720, "local_header_offset": 310,
    },
    "data_zenodo/Figure Data/Figure2_Data.xlsx": {
        "crc32": 0x377026A2, "compressed_size": 70_243,
        "uncompressed_size": 73_775, "local_header_offset": 12_149,
    },
    "data_zenodo/Figure Data/Figure3_Data.xlsx": {
        "crc32": 0xDC2B5ECE, "compressed_size": 7_379,
        "uncompressed_size": 9_841, "local_header_offset": 82_792,
    },
    "data_zenodo/Figure Data/FigureS2_Data.xlsx": {
        "crc32": 0x81B191A1, "compressed_size": 59_446,
        "uncompressed_size": 62_130, "local_header_offset": 1_241_135,
    },
}


def _request(url: str, start: int | None = None, end: int | None = None) -> bytes:
    headers = {"User-Agent": "Project-Aleph-research/1.0"}
    if start is not None and end is not None:
        headers["Range"] = f"bytes={start}-{end}"
    for attempt in range(6):
        try:
            with urllib.request.urlopen(
                urllib.request.Request(url, headers=headers), timeout=120
            ) as response:
                payload = response.read()
                if start is not None and end is not None:
                    expected = end - start + 1
                    content_range = response.headers.get("Content-Range")
                    if len(payload) != expected or content_range != (
                        f"bytes {start}-{end}/{ARCHIVE_SIZE}"
                    ):
                        raise ValueError("Zenodo range response contract changed")
                return payload
        except urllib.error.HTTPError as error:
            if error.code != 429 or attempt == 5:
                raise
            time.sleep(min(60, int(error.headers.get("Retry-After", "5"))))
    raise RuntimeError("unreachable")


def _safe(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts and "\\" not in name


def _central_members(payload: bytes) -> dict[str, dict[str, int]]:
    members: dict[str, dict[str, int]] = {}
    position = 0
    while position < len(payload):
        values = struct.unpack_from("<4s6H3L5H2L", payload, position)
        if values[0] != b"PK\x01\x02":
            raise ValueError("invalid ZIP central-directory signature")
        name_length, extra_length, comment_length = values[10:13]
        name_start = position + 46
        name = payload[name_start:name_start + name_length].decode("utf-8")
        if not _safe(name) or name in members:
            raise ValueError(f"unsafe or duplicate ZIP member: {name}")
        members[name] = {
            "compression_method": values[4], "crc32": values[7],
            "compressed_size": values[8], "uncompressed_size": values[9],
            "local_header_offset": values[16],
        }
        position += 46 + name_length + extra_length + comment_length
    if position != len(payload):
        raise ValueError("ZIP central directory has trailing bytes")
    return members


def _member_bytes(record: dict[str, int]) -> bytes:
    offset = record["local_header_offset"]
    header = _request(ARCHIVE_URL, offset, offset + 29)
    values = struct.unpack("<4s5H3L2H", header)
    if values[0] != b"PK\x03\x04" or values[3] != record["compression_method"]:
        raise ValueError("ZIP local header changed")
    name_length, extra_length = values[-2:]
    data_start = offset + 30 + name_length + extra_length
    compressed = _request(
        ARCHIVE_URL, data_start, data_start + record["compressed_size"] - 1
    )
    if record["compression_method"] == 8:
        payload = zlib.decompress(compressed, -zlib.MAX_WBITS)
    elif record["compression_method"] == 0:
        payload = compressed
    else:
        raise ValueError("unsupported ZIP compression method")
    if len(payload) != record["uncompressed_size"]:
        raise ValueError("ZIP member uncompressed size changed")
    if binascii.crc32(payload) & 0xFFFFFFFF != record["crc32"]:
        raise ValueError("ZIP member CRC32 mismatch")
    return payload


def acquire(output_dir: Path) -> dict[str, object]:
    metadata = json.loads(_request(RECORD_API))
    files = metadata.get("files", [])
    if len(files) != 1 or files[0].get("key") != "data_zenodo.zip":
        raise ValueError("Zenodo record file contract changed")
    if files[0].get("size") != ARCHIVE_SIZE or files[0].get("checksum") != f"md5:{ARCHIVE_MD5}":
        raise ValueError("Zenodo archive metadata changed")
    central = _request(
        ARCHIVE_URL,
        CENTRAL_DIRECTORY["offset"],
        CENTRAL_DIRECTORY["offset"] + CENTRAL_DIRECTORY["size"] - 1,
    )
    members = _central_members(central)
    if len(members) != CENTRAL_DIRECTORY["member_count"]:
        raise ValueError("Zenodo ZIP member count changed")
    for name, expected in EXPECTED.items():
        observed = members.get(name)
        if observed is None or observed["compression_method"] != 8:
            raise ValueError(f"missing expected Zenodo table: {name}")
        for key, value in expected.items():
            if observed[key] != value:
                raise ValueError(f"Zenodo member changed: {name}/{key}")

    output_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    for name in EXPECTED:
        payload = _member_bytes(members[name])
        destination = output_dir / PurePosixPath(name).name
        destination.write_bytes(payload)
        records[destination.name] = {
            "archive_member": name,
            "size_bytes": len(payload),
            "crc32": f"{binascii.crc32(payload) & 0xFFFFFFFF:08x}",
            "sha256": hashlib.sha256(payload).hexdigest(),
        }
    return {
        "schema": "aleph.outer_library.zenodo_7432971_figure_tables.v1",
        "source_doi": "10.5281/zenodo.7432971",
        "license_id": "cc-by-4.0",
        "source_archive": {
            "name": "data_zenodo.zip", "size_bytes": ARCHIVE_SIZE,
            "official_metadata_md5": ARCHIVE_MD5,
            "full_archive_downloaded": False,
            "central_directory_range_verified": True,
            "member_count": CENTRAL_DIRECTORY["member_count"],
        },
        "extracted_tables": records,
        "retrieved_at": "2026-08-06",
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
    args.receipt.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["extracted_tables"], sort_keys=True))


if __name__ == "__main__":
    main()
