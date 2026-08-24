#!/usr/bin/env python
"""Download a fixed, tiny, openly licensed raw/source-data pilot with receipts only in Git."""

from __future__ import annotations

import hashlib
import json
import os
import urllib.request
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
OFFICIAL = HERE / "official_metadata_receipts.jsonl"
OUTPUT = ROOT / "data" / "external_training" / "experiment_factory" / "datasets" / "pilot"
RECEIPTS = HERE / "pilot_download_receipts.jsonl"
BYTE_CAP = 5_000_000
USER_AGENT = "Project-Aleph-public-dataset-pilot/1.0 (bounded open-data client)"
PLAN = {
    ("Figshare", "10.6084/m9.figshare.16826740"): {"*"},
    ("Figshare", "10.6084/m9.figshare.27241950"): {"*"},
    ("Figshare", "10.6084/m9.figshare.28309337"): {"Figure 4A-bottom panel_2G3-Occludin in ACTG1 KO 2G3 WT Red-ACTG1 Green-occludin snap1.czi"},
    ("Zenodo", "10.5281/zenodo.18779816"): {"*"},
}
OPEN_LICENSES = {"cc-by-4.0", "cc by 4.0", "gpl-3.0-or-later", "mit-license"}


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def load_official() -> dict[tuple[str, str], dict[str, Any]]:
    rows = [json.loads(line) for line in OFFICIAL.read_text(encoding="utf-8").splitlines() if line.strip()]
    return {(row["provider"], row["accession"]): row for row in rows}


def main() -> int:
    official = load_official()
    selected = []
    for key, names in PLAN.items():
        record = official[key]
        metadata = record["metadata"]
        license_id = str(metadata.get("license_identifier") or "")
        if license_id.lower() not in OPEN_LICENSES:
            raise RuntimeError(f"pilot licence is not allowlisted: {key} {license_id!r}")
        for item in metadata.get("files") or []:
            if "*" in names or item.get("name") in names:
                if not item.get("url") or not isinstance(item.get("size"), int):
                    raise RuntimeError(f"missing official URL or size: {key} {item}")
                selected.append((key, license_id, item))
    declared = sum(item[2]["size"] for item in selected)
    if declared > BYTE_CAP:
        raise RuntimeError(f"declared pilot size {declared} exceeds {BYTE_CAP}")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    receipts = []
    consumed = 0
    for (provider, accession), license_id, item in sorted(selected, key=lambda value: (value[0], value[2]["name"])):
        request = urllib.request.Request(item["url"], headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read(BYTE_CAP - consumed + 1)
        consumed += len(payload)
        if consumed > BYTE_CAP:
            raise RuntimeError("actual downloads exceeded byte cap")
        if len(payload) != item["size"]:
            raise RuntimeError(f"size mismatch for {item['name']}: {len(payload)} != {item['size']}")
        digest = hashlib.sha256(payload).hexdigest()
        destination = OUTPUT / digest
        if destination.exists():
            if hashlib.sha256(destination.read_bytes()).hexdigest() != digest:
                raise RuntimeError(f"existing object digest mismatch: {destination}")
        else:
            temporary = destination.with_suffix(".part")
            temporary.write_bytes(payload)
            os.replace(temporary, destination)
        receipts.append({
            "schema": "aleph.external_training.dataset_pilot_receipt.v1",
            "authority_status": "proposed",
            "provider": provider,
            "accession": accession,
            "file_name": item["name"],
            "bytes": len(payload),
            "sha256": digest,
            "license_identifier": license_id,
            "download_url": item["url"],
            "storage": "ignored_content_addressed_local_object",
        })
    RECEIPTS.write_bytes(b"".join(canonical_bytes(row) for row in receipts))
    print(json.dumps({"files": len(receipts), "bytes": consumed, "cap": BYTE_CAP}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
