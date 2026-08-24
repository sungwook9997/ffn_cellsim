#!/usr/bin/env python3
"""Resolve IDR0072 originals by the official plate paths and FLEX well/field convention."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.error
import urllib.parse
import urllib.request


REPOSITORY_COMMIT = "85d154650a26f011b57f13f91e4be11d1f9dde12"
RAW_ROOT = f"https://raw.githubusercontent.com/IDR/idr0072-schormann-subcellref/{REPOSITORY_COMMIT}"
FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/IDR/idr0072-schormann-subcellref"


def fetch(url: str, range_header: str | None = None, attempts: int = 5) -> tuple[bytes, dict[str, str], int]:
    headers = {"User-Agent": "Project-Aleph-research/1.0"}
    if range_header:
        headers["Range"] = range_header
    request = urllib.request.Request(url, headers=headers)
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read(), dict(response.headers.items()), response.status
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def load_plate_map(screen: str) -> tuple[dict[str, str], dict]:
    name = f"idr0072-{screen}-plates.tsv"
    url = f"{RAW_ROOT}/{screen}/{name}"
    payload, _, _ = fetch(url)
    mapping = {}
    marker = "/idr0072-schormann-subcellref/"
    for line in payload.decode("utf-8").splitlines():
        if line.strip():
            plate, repository_path = line.split("\t")
            mapping[plate] = repository_path.split(marker, 1)[1].strip("/")
    return mapping, {"url": url, "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(), "plate_count": len(mapping)}


def expected_filename(row: dict) -> str:
    match = re.fullmatch(r"([a-pA-P])(\d{1,2})", str(row["well_name"]).strip())
    if not match:
        raise ValueError(f"invalid 384-well coordinate: {row['well_name']}")
    row_number = ord(match.group(1).casefold()) - ord("a") + 1
    column_number = int(match.group(2))
    if not 1 <= column_number <= 24:
        raise ValueError(f"invalid 384-well column: {row['well_name']}")
    suffix = 1 if int(row["screen_id"]) == 2952 else 0
    return f"{row_number:03d}{column_number:03d}{suffix:03d}.flex"


def probe(row: dict) -> dict:
    _, headers, status = fetch(row["source_url"], range_header="bytes=0-0")
    content_range = headers.get("Content-Range")
    if status != 206 or not content_range or "/" not in content_range:
        raise ValueError(f"source does not support a one-byte range: {row['source_url']}")
    item = dict(row)
    item["source_total_bytes"] = int(content_range.rsplit("/", 1)[1])
    return item


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=32)
    args = parser.parse_args()
    payload = args.manifest.read_bytes()
    rows = [json.loads(line) for line in payload.splitlines() if line]
    maps, receipts = {}, {}
    for screen_id, name in ((2952, "screenA"), (2953, "screenB")):
        maps[screen_id], receipts[name] = load_plate_map(name)
    candidates = []
    for row in rows:
        relative = maps[int(row["screen_id"])][row["plate_name"]]
        filename = expected_filename(row)
        item = dict(row)
        item["flex_filename"] = filename
        item["source_url"] = f"{FTP_ROOT}/{urllib.parse.quote(relative, safe='/')}/{filename}"
        item["local_relative_path"] = f"screen_{row['screen_id']}/plate_{row['plate_id']}/{filename}"
        candidates.append(item)
    if len({row["image_id"] for row in candidates}) != len(candidates):
        raise ValueError("duplicate image IDs")
    if len({row["source_url"] for row in candidates}) != len(candidates):
        raise ValueError("well-derived FLEX URLs are not one-to-one")
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        resolved = list(pool.map(probe, candidates))
    resolved.sort(key=lambda row: row["selection_key"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in resolved), encoding="utf-8")
    report = {
        "schema": "aleph.outer_library.idr0072_if_flex_resolution.v2",
        "status": "one_to_one_flex_sources_resolved_before_prediction",
        "metadata_repository_commit": REPOSITORY_COMMIT,
        "metadata_receipts": receipts,
        "input_manifest_sha256": hashlib.sha256(payload).hexdigest(),
        "input_rows": len(rows),
        "resolved_rows": len(resolved),
        "unique_image_ids": len({row["image_id"] for row in resolved}),
        "unique_source_urls": len({row["source_url"] for row in resolved}),
        "unique_local_paths": len({row["local_relative_path"] for row in resolved}),
        "source_total_bytes_sum": sum(row["source_total_bytes"] for row in resolved),
        "resolved_manifest_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "model_predictions_made": 0,
        "pristine_holdout_authority": False,
        "aleph_parameter_authority": "none"
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "resolved": len(resolved), "unique_urls": report["unique_source_urls"], "source_GB": round(report["source_total_bytes_sum"] / 1e9, 3)}, sort_keys=True))


if __name__ == "__main__":
    main()
