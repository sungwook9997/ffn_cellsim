#!/usr/bin/env python3
"""Acquire the frozen replacement PBMC confirmation dataset GSE158055."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from download_pbmc_two_lab_recalibration import download, download_segmented, sha256


BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158055/suppl"
FILES = (
    "GSE158055_cell_annotation.csv.gz",
    "GSE158055_covid19_barcodes.tsv.gz",
    "GSE158055_covid19_counts.mtx.gz",
    "GSE158055_covid19_features.tsv.gz",
    "GSE158055_sample_metadata.xlsx",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol.get("frozen_before_GSE158055_expression_or_cell_annotation_access"):
        raise ValueError("replacement confirmation protocol was not frozen")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for name in FILES:
        url = f"{BASE}/{name}"
        path = args.output_dir / name
        if name == "GSE158055_covid19_counts.mtx.gz":
            download_segmented(url, path, workers=args.workers)
        else:
            download(url, path)
        if name.endswith(".gz"):
            with gzip.open(path, "rb") as handle:
                handle.read(1)
        records.append({
            "name": name, "url": url, "size_bytes": path.stat().st_size,
            "sha256": sha256(path), "gzip_header_verified": name.endswith(".gz"),
        })
    receipt = {
        "schema": "aleph.outer_library.gse158055_confirmation_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "protocol_frozen_commit": "469ccb0",
        "files": records,
        "role": "single_use_untouched_confirmation",
        "labels_used_during_acquisition": False,
        "aleph_authority": "none",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "files": len(records), "bytes": sum(row["size_bytes"] for row in records)
    }, sort_keys=True))


if __name__ == "__main__":
    main()
