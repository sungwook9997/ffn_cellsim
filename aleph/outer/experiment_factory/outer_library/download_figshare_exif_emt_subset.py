#!/usr/bin/env python3
"""Range-extract the ExIF EMT layouts without downloading the 209 GB release.

The Figshare release stores its two payloads as ZIP64 archives.  The central
directories are small enough to verify independently; selected members are
then fetched by their local-header offsets and checked against ZIP CRC32.
"""

from __future__ import annotations

import argparse
import binascii
from concurrent.futures import as_completed, ThreadPoolExecutor
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import struct
import time
import urllib.error
import urllib.request
import zlib


ARTICLE_API = "https://api.figshare.com/v2/articles/26500210"
ARTICLE_ID = 26_500_210
ARTICLE_DOI = "10.6084/m9.figshare.26500210.v1"
PAPER_DOI = "10.1038/s41467-025-59592-7"

ARCHIVES = {
    "images": {
        "file_id": 52_855_610,
        "name": "Images.zip",
        "size": 193_981_784_985,
        "md5": "94b60280816b0e73005050b2973b9851",
        "member_count": 328_681,
        "central_offset": 193_944_261_803,
        "central_size": 37_523_084,
        "central_sha256": "e3e3dfaf71319550a2c53e295414c441f818b83a5233d6f2272833b7631815d5",
    },
    "quantitative": {
        "file_id": 53_517_041,
        "name": "Quantitative_data.zip",
        "size": 15_362_201_302,
        "md5": "dfc80eb6dacdb66844fa6cc31bcdab09",
        "member_count": 6_660,
        "central_offset": 15_361_223_830,
        "central_size": 977_374,
        "central_sha256": "6e2f72c50b074898e8756fd3ee7d68e73af150ee54b44a3b03e7ff22dae7a304",
    },
}

EXPECTED_LAYOUTS = {
    "Images/raw/ExIF_EMT/220902_A549_EM_ActExp3_Plate_Layout.xlsx": {
        "compression_method": 8,
        "crc32": 0x568FFAF9,
        "compressed_size": 15_447,
        "uncompressed_size": 26_756,
        "local_header_offset": 21_549_802_138,
    },
    "Images/raw/Mulltiplexed_EMT/layout.xlsx": {
        "compression_method": 8,
        "crc32": 0x8A1CD030,
        "compressed_size": 6_816,
        "uncompressed_size": 9_281,
        "local_header_offset": 193_944_254_918,
    },
}

RAW_FIELDS = ("0000", "0018", "0036")
RAW_IMAGE_PATTERN = re.compile(
    r"^Images/raw/ExIF_EMT/220825_([A-H])(0[1-6])_"
    r"(0000|0018|0036)_(TD|ex405nm|ex488nm|ex555nm|ex647nm)\.tif$"
)


def _url(file_id: int) -> str:
    return f"https://ndownloader.figshare.com/files/{file_id}"


def _request(url: str, *, start: int | None = None, end: int | None = None,
             total_size: int | None = None) -> bytes:
    headers = {"User-Agent": "Project-Aleph-research/1.0"}
    if start is not None and end is not None:
        headers["Range"] = f"bytes={start}-{end}"
    for attempt in range(6):
        try:
            with urllib.request.urlopen(
                urllib.request.Request(url, headers=headers), timeout=180
            ) as response:
                payload = response.read()
                if start is not None and end is not None:
                    expected = end - start + 1
                    content_range = response.headers.get("Content-Range")
                    if len(payload) != expected or content_range != (
                        f"bytes {start}-{end}/{total_size}"
                    ):
                        raise ValueError("Figshare range-response contract changed")
                return payload
        except urllib.error.HTTPError as error:
            if error.code not in {429, 500, 502, 503, 504} or attempt == 5:
                raise
            time.sleep(min(30, 2 ** attempt))
    raise RuntimeError("unreachable")


def _safe(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts and "\\" not in name


def _zip64_values(extra: bytes) -> list[int]:
    position = 0
    while position < len(extra):
        if position + 4 > len(extra):
            raise ValueError("truncated ZIP extra field")
        field_id, size = struct.unpack_from("<HH", extra, position)
        payload = extra[position + 4:position + 4 + size]
        if len(payload) != size:
            raise ValueError("truncated ZIP extra payload")
        if field_id == 1:
            if size % 8:
                raise ValueError("invalid ZIP64 extra size")
            return [
                struct.unpack_from("<Q", payload, index)[0]
                for index in range(0, size, 8)
            ]
        position += 4 + size
    return []


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
        extra = payload[
            name_start + name_length:name_start + name_length + extra_length
        ]
        if not _safe(name) or name in members:
            raise ValueError(f"unsafe or duplicate ZIP member: {name}")
        record = {
            "compression_method": values[4],
            "crc32": values[7],
            "compressed_size": values[8],
            "uncompressed_size": values[9],
            "local_header_offset": values[16],
        }
        zip64 = iter(_zip64_values(extra))
        for key in ("uncompressed_size", "compressed_size", "local_header_offset"):
            if record[key] == 0xFFFFFFFF:
                try:
                    record[key] = next(zip64)
                except StopIteration as error:
                    raise ValueError(f"missing ZIP64 value for {name}/{key}") from error
        members[name] = record
        position += 46 + name_length + extra_length + comment_length
    if position != len(payload):
        raise ValueError("ZIP central directory has trailing bytes")
    return members


def _central_payload(spec: dict[str, int | str], workers: int) -> bytes:
    offset = int(spec["central_offset"])
    size = int(spec["central_size"])
    workers = max(1, min(workers, size))
    base, remainder = divmod(size, workers)
    ranges = []
    cursor = offset
    for index in range(workers):
        length = base + (1 if index < remainder else 0)
        ranges.append((cursor, cursor + length - 1))
        cursor += length
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                _request, _url(int(spec["file_id"])), start=start, end=end,
                total_size=int(spec["size"]),
            )
            for start, end in ranges
        ]
        payload = b"".join(future.result() for future in futures)
    if len(payload) != size:
        raise ValueError("central-directory size changed")
    if hashlib.sha256(payload).hexdigest() != spec["central_sha256"]:
        raise ValueError("central-directory digest changed")
    return payload


def _member_bytes(spec: dict[str, int | str], record: dict[str, int]) -> bytes:
    url = _url(int(spec["file_id"]))
    total_size = int(spec["size"])
    offset = record["local_header_offset"]
    header = _request(url, start=offset, end=offset + 29, total_size=total_size)
    values = struct.unpack("<4s5H3L2H", header)
    if values[0] != b"PK\x03\x04" or values[3] != record["compression_method"]:
        raise ValueError("ZIP local header changed")
    name_length, extra_length = values[-2:]
    data_start = offset + 30 + name_length + extra_length
    compressed = _request(
        url,
        start=data_start,
        end=data_start + record["compressed_size"] - 1,
        total_size=total_size,
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


def acquire(output_dir: Path, *, workers: int = 8, raw_workers: int = 16,
            field_ids: tuple[str, ...] = RAW_FIELDS) -> dict[str, object]:
    metadata = json.loads(_request(ARTICLE_API))
    if metadata.get("id") != ARTICLE_ID or metadata.get("doi") != ARTICLE_DOI:
        raise ValueError("Figshare article identity changed")
    if metadata.get("license", {}).get("name") != "CC BY 4.0":
        raise ValueError("Figshare license changed")
    files = {entry["name"]: entry for entry in metadata.get("files", [])}
    for spec in ARCHIVES.values():
        observed = files.get(spec["name"])
        if observed is None:
            raise ValueError(f"missing Figshare archive {spec['name']}")
        if (
            observed.get("id") != spec["file_id"]
            or observed.get("size") != spec["size"]
            or observed.get("supplied_md5") != spec["md5"]
        ):
            raise ValueError(f"Figshare archive metadata changed: {spec['name']}")

    central = _central_payload(ARCHIVES["images"], workers)
    members = _central_members(central)
    if len(members) != ARCHIVES["images"]["member_count"]:
        raise ValueError("Images.zip member count changed")
    for name, expected in EXPECTED_LAYOUTS.items():
        if members.get(name) != expected:
            raise ValueError(f"ExIF layout member changed: {name}")
    if not field_ids or not set(field_ids).issubset(RAW_FIELDS):
        raise ValueError(f"field_ids must be a nonempty subset of {RAW_FIELDS}")
    raw_names = sorted(
        name
        for name in members
        if (match := RAW_IMAGE_PATTERN.match(name)) is not None
        and match.group(3) in field_ids
    )
    if len(raw_names) != 8 * 6 * len(field_ids) * 5:
        raise ValueError("ExIF raw pilot cardinality changed")
    well_fields = {
        match.group(1) + match.group(2) + ":" + match.group(3)
        for name in raw_names
        if (match := RAW_IMAGE_PATTERN.match(name)) is not None
    }
    if len(well_fields) != 8 * 6 * len(field_ids):
        raise ValueError("ExIF raw pilot lost a well/field group")

    output_dir.mkdir(parents=True, exist_ok=True)
    extracted = {}
    for name in EXPECTED_LAYOUTS:
        payload = _member_bytes(ARCHIVES["images"], members[name])
        destination = output_dir / PurePosixPath(name).name
        destination.write_bytes(payload)
        extracted[destination.name] = {
            "archive_member": name,
            "size_bytes": len(payload),
            "crc32": f"{binascii.crc32(payload) & 0xFFFFFFFF:08x}",
            "sha256": hashlib.sha256(payload).hexdigest(),
        }

    raw_dir = output_dir / "raw_exif_emt"
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw_records: dict[str, dict[str, object]] = {}
    missing_raw_names = []
    for name in raw_names:
        destination = raw_dir / PurePosixPath(name).name
        if destination.exists():
            payload = destination.read_bytes()
            expected = members[name]
            if (
                len(payload) == expected["uncompressed_size"]
                and binascii.crc32(payload) & 0xFFFFFFFF == expected["crc32"]
            ):
                raw_records[name] = {
                    "size_bytes": len(payload),
                    "crc32": f"{expected['crc32']:08x}",
                    "sha256": hashlib.sha256(payload).hexdigest(),
                }
                continue
        missing_raw_names.append(name)
    raw_workers = max(1, raw_workers)
    with ThreadPoolExecutor(max_workers=raw_workers) as executor:
        future_names = {
            executor.submit(_member_bytes, ARCHIVES["images"], members[name]): name
            for name in missing_raw_names
        }
        for future in as_completed(future_names):
            name = future_names[future]
            payload = future.result()
            destination = raw_dir / PurePosixPath(name).name
            destination.write_bytes(payload)
            raw_records[name] = {
                "size_bytes": len(payload),
                "crc32": f"{binascii.crc32(payload) & 0xFFFFFFFF:08x}",
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
            if len(raw_records) % 50 == 0 or len(raw_records) == len(raw_names):
                print(f"verified raw ExIF members: {len(raw_records)}/{len(raw_names)}", flush=True)

    ordered_raw_records = {name: raw_records[name] for name in raw_names}
    raw_manifest_sha256 = hashlib.sha256(
        (json.dumps(ordered_raw_records, sort_keys=True) + "\n").encode("utf-8")
    ).hexdigest()

    return {
        "schema": "aleph.outer_library.figshare_exif_emt_layouts.v1",
        "source": {
            "article_id": ARTICLE_ID,
            "dataset_doi": ARTICLE_DOI,
            "paper_doi": PAPER_DOI,
            "license_id": "cc-by-4.0",
        },
        "archives": {
            key: {
                "name": spec["name"],
                "size_bytes": spec["size"],
                "official_metadata_md5": spec["md5"],
                "full_archive_downloaded": False,
                "central_directory_range_verified": key == "images",
                "central_directory_sha256": spec["central_sha256"],
                "member_count": spec["member_count"],
            }
            for key, spec in ARCHIVES.items()
        },
        "extracted_layouts": extracted,
        "raw_image_pilot": {
            "selection": {
                "plate_rows": list("ABCDEFGH"),
                "plate_columns": list(range(1, 7)),
                "field_ids": list(field_ids),
                "channels": ["TD", "ex405nm", "ex488nm", "ex555nm", "ex647nm"],
            },
            "field_groups": len(well_fields),
            "member_count": len(ordered_raw_records),
            "uncompressed_size_bytes": sum(
                int(record["size_bytes"]) for record in ordered_raw_records.values()
            ),
            "manifest_sha256": raw_manifest_sha256,
            "members": ordered_raw_records,
        },
        "retrieved_at": "2026-08-06",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--raw-workers", type=int, default=16)
    parser.add_argument("--field-ids", nargs="+", default=list(RAW_FIELDS))
    args = parser.parse_args()
    report = acquire(
        args.output_dir,
        workers=args.workers,
        raw_workers=args.raw_workers,
        field_ids=tuple(args.field_ids),
    )
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report["extracted_layouts"], sort_keys=True))


if __name__ == "__main__":
    main()
