#!/usr/bin/env python3
"""Acquire the GSE325309 workbook only after the holdout protocol was committed."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request


NAME = "GSE325309_processed_data.xlsx"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325309/suppl/" + NAME
SIZE = 6_266_161
SHA256 = "0d9d636dff0fba25794c41150bbe56be3e57079cdf5bf924548c160faa64d71e"
PROTOCOL_SHA256 = "9eb83ce6f04bf9db7d1b63a0744a9021689de522fe28067778116414d99027e3"
PROTOCOL_COMMIT = "7984e5b"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(output_dir: Path, protocol_path: Path) -> dict[str, object]:
    if _sha256(protocol_path) != PROTOCOL_SHA256:
        raise ValueError("GSE325309 frozen protocol changed after outcome access")
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / NAME
    if not path.is_file() or path.stat().st_size != SIZE:
        temporary = path.with_suffix(".xlsx.part")
        request = urllib.request.Request(URL, headers={"User-Agent": "Project-Aleph-research/1.0"})
        with urllib.request.urlopen(request, timeout=300) as response, temporary.open("wb") as handle:
            while block := response.read(1024 * 1024):
                handle.write(block)
        if temporary.stat().st_size != SIZE:
            raise ValueError("GSE325309 workbook size changed")
        temporary.replace(path)
    digest = _sha256(path)
    if digest != SHA256:
        raise ValueError("GSE325309 workbook digest changed")
    return {
        "schema": "aleph.outer_library.gse325309_pristine_acquisition.v1",
        "accession": "GSE325309",
        "retrieved_at": "2026-08-07",
        "file": {"name": NAME, "url": URL, "size": SIZE, "sha256": digest},
        "freeze_proof": {
            "protocol_sha256": PROTOCOL_SHA256,
            "protocol_commit_before_download": PROTOCOL_COMMIT,
            "expression_downloaded_after_protocol_commit": True,
        },
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    result = acquire(args.output_dir, args.protocol)
    args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"size": SIZE, "sha256": SHA256}, sort_keys=True))


if __name__ == "__main__":
    main()
