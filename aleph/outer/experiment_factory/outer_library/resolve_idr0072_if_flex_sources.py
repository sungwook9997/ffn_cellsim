#!/usr/bin/env python3
"""Resolve frozen IDR0072 image IDs to original FLEX source URLs without pixels."""

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


def fetch(url: str, attempts: int = 5) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "Project-Aleph-research/1.0"})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def load_plate_map(screen: str) -> tuple[dict[str, str], dict]:
    name = f"idr0072-{screen}-plates.tsv"
    url = f"{RAW_ROOT}/{screen}/{name}"
    payload = fetch(url)
    mapping: dict[str, str] = {}
    for line in payload.decode("utf-8").splitlines():
        if not line.strip():
            continue
        plate, repository_path = line.split("\t")
        marker = "/idr0072-schormann-subcellref/"
        if marker not in repository_path:
            raise ValueError(f"unexpected repository path: {repository_path}")
        relative = repository_path.split(marker, 1)[1].strip("/")
        if plate in mapping:
            raise ValueError(f"duplicate plate name: {plate}")
        mapping[plate] = relative
    return mapping, {"url": url, "bytes": len(payload), "sha256": digest(payload), "plate_count": len(mapping)}


def resolve_filename(image_id: int) -> str:
    payload = json.loads(fetch(f"https://idr.openmicroscopy.org/webclient/imgData/{image_id}/"))
    exception = payload.get("Exception", "")
    match = re.search(r"/([^/]+\.flex)$", exception, re.IGNORECASE)
    if not match:
        raise ValueError(f"image {image_id} did not expose a FLEX repository basename")
    return match.group(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol["frozen_after_one_file_metadata_probe_but_before_pixel_array_access_or_prediction"]:
        raise ValueError("optical channel protocol was not frozen")
    manifest_payload = args.manifest.read_bytes()
    if digest(manifest_payload) != protocol["source_contract"]["manifest_sha256"]:
        raise ValueError("frozen manifest hash mismatch")
    rows = [json.loads(line) for line in manifest_payload.splitlines() if line]
    screen_maps: dict[int, dict[str, str]] = {}
    metadata_receipts = {}
    for screen_id, screen_name in ((2952, "screenA"), (2953, "screenB")):
        screen_maps[screen_id], metadata_receipts[screen_name] = load_plate_map(screen_name)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        filenames = list(pool.map(lambda row: resolve_filename(int(row["image_id"])), rows))
    resolved = []
    missing_plates = []
    for row, filename in zip(rows, filenames, strict=True):
        relative = screen_maps[int(row["screen_id"])].get(row["plate_name"])
        if relative is None:
            missing_plates.append(row["plate_name"])
            continue
        item = dict(row)
        item["flex_filename"] = filename
        item["source_url"] = f"{FTP_ROOT}/{urllib.parse.quote(relative, safe='/')}/{urllib.parse.quote(filename)}"
        item["local_relative_path"] = f"screen_{row['screen_id']}/plate_{row['plate_id']}/{filename}"
        resolved.append(item)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in resolved), encoding="utf-8")
    report = {
        "schema": "aleph.outer_library.idr0072_if_flex_resolution.v1",
        "status": "flex_sources_resolved_before_pixel_access" if len(resolved) == len(rows) else "refused_unresolved_flex_sources",
        "metadata_repository_commit": REPOSITORY_COMMIT,
        "metadata_receipts": metadata_receipts,
        "input_manifest_sha256": digest(manifest_payload),
        "input_rows": len(rows),
        "resolved_rows": len(resolved),
        "missing_plate_names": sorted(set(missing_plates)),
        "resolved_manifest_sha256": digest(args.output.read_bytes()),
        "additional_pixel_files_downloaded": 0,
        "model_predictions_made": 0,
        "pristine_holdout_authority": False,
        "aleph_parameter_authority": "none"
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "resolved": len(resolved), "missing_plates": len(set(missing_plates))}, sort_keys=True))


if __name__ == "__main__":
    main()
