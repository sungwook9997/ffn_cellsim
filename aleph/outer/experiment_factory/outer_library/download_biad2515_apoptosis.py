#!/usr/bin/env python3
"""Acquire verified S-BIAD2515 metadata and author-processed apoptosis profiles."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
import urllib.request


ANALYSIS_COMMIT = "78d77675abcc056d0d3e7f3dded9ce49609e8f6c"
FILES = {
    "biostudies_metadata.json": {
        "url": "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-BIAD2515",
        "size": 18_380,
        "sha256": "187646cb2e94e3b78bd9f7d1796ab31832a9f9e432139eaa15e95939ca7b1b81",
    },
    "filelist_bioimage.json": {
        "url": "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/515/S-BIAD2515/Files/filelist_bioimage.json",
        "size": 1_042_081,
        "sha256": "0beec09551350ceb4326e029a99a2e57be06f0fd4a4b4dced5b454f9af49d6de",
    },
    "train.parquet": {
        "url": (
            "https://raw.githubusercontent.com/WayScience/"
            f"live_cell_timelapse_apoptosis_analysis/{ANALYSIS_COMMIT}/"
            "5.bulk_timelapse_model/data_splits/train.parquet"
        ),
        "size": 2_869_835,
        "sha256": "8a92d3d34980a8f115fac5a1ec48e62a50dc914221b82b06308091373d48eebe",
    },
    "test.parquet": {
        "url": (
            "https://raw.githubusercontent.com/WayScience/"
            f"live_cell_timelapse_apoptosis_analysis/{ANALYSIS_COMMIT}/"
            "5.bulk_timelapse_model/data_splits/test.parquet"
        ),
        "size": 2_620_033,
        "sha256": "69fabdce42807fd88ce504f34873d840a65503bc4b7d4cf9a41eb346a1e82777",
    },
    "train_test_wells.parquet": {
        "url": (
            "https://raw.githubusercontent.com/WayScience/"
            f"live_cell_timelapse_apoptosis_analysis/{ANALYSIS_COMMIT}/"
            "5.bulk_timelapse_model/data_splits/train_test_wells.parquet"
        ),
        "size": 1_839,
        "sha256": "7255db106bf47358dbe773b08bf6e709889f84eebd2c9f6d8ff96fb19af9024d",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download(url: str, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(
        url, headers={"User-Agent": "Project-Aleph-research/1.0"}
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                payload = response.read()
            temporary.write_bytes(payload)
            temporary.replace(destination)
            return
        except Exception:
            if attempt == 4:
                raise
            time.sleep(0.5 * 2 ** attempt)


def acquire(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    receipts: dict[str, object] = {}
    for name, expected in FILES.items():
        path = output_dir / name
        if (
            not path.exists()
            or path.stat().st_size != expected["size"]
            or _sha256(path) != expected["sha256"]
        ):
            _download(str(expected["url"]), path)
        actual = {"size": path.stat().st_size, "sha256": _sha256(path)}
        if actual != {"size": expected["size"], "sha256": expected["sha256"]}:
            raise ValueError(f"S-BIAD2515 checksum mismatch for {name}: {actual}")
        receipts[name] = {**actual, "url": expected["url"]}

    filelist = json.loads((output_dir / "filelist_bioimage.json").read_text())
    raw_bytes = sum(int(item["size"]) for item in filelist)
    if len(filelist) != 6_240 or raw_bytes != 45_085_609_920:
        raise ValueError("S-BIAD2515 raw-image file list changed")
    return {
        "schema": "aleph.outer_library.biad2515_acquisition_receipt.v1",
        "accession": "S-BIAD2515",
        "doi": "10.6019/S-BIAD2515",
        "release_date": "2026-01-21",
        "license_id": "CC0-1.0",
        "author_analysis_commit": ANALYSIS_COMMIT,
        "files": receipts,
        "official_raw_image_file_count": len(filelist),
        "official_raw_image_total_bytes": raw_bytes,
        "raw_images_locally_acquired": False,
        "author_processed_profiles_locally_acquired": True,
        "retrieved_at": "2026-08-06",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    receipt = acquire(args.output_dir)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
