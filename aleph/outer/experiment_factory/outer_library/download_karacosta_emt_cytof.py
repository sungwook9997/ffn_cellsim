#!/usr/bin/env python3
"""Download the pinned Karacosta HCC827 EMT CyTOF Source Data workbook."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request


PAPER_DOI = "10.1038/s41467-019-13441-6"
NAME = "41467_2019_13441_MOESM4_ESM.xlsx"
URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41467-019-13441-6/MediaObjects/" + NAME
)
SIZE = 79_066_708
SHA256 = "bd6e239008b40cb91b722f16a43c81f43f043adc34f6807a74d094746bfaadfa"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def acquire(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / NAME
    if not destination.is_file() or destination.stat().st_size != SIZE:
        temporary = destination.with_suffix(".xlsx.part")
        request = urllib.request.Request(
            URL, headers={"User-Agent": "Project-Aleph-research/1.0"}
        )
        with urllib.request.urlopen(request, timeout=300) as response, temporary.open("wb") as handle:
            while block := response.read(1024 * 1024):
                handle.write(block)
        if temporary.stat().st_size != SIZE:
            raise ValueError("Karacosta Source Data size changed")
        temporary.replace(destination)
    observed = _sha256(destination)
    if observed != SHA256:
        raise ValueError("Karacosta Source Data digest changed")
    return {
        "schema": "aleph.outer_library.karacosta_emt_cytof_acquisition.v1",
        "source": {
            "paper_doi": PAPER_DOI,
            "paper_url": "https://www.nature.com/articles/s41467-019-13441-6",
            "pmcid": "PMC6898514",
            "provider": "Springer Nature article Source Data",
            "license_id": "cc-by-4.0",
            "cytobank_experiment": "https://stanford.cytobank.org/cytobank/experiments/26555",
        },
        "files": [{
            "local_name": NAME, "url": URL, "size": SIZE, "sha256": observed,
            "role": "author_Source_Data_raw_and_transformed_single_cell_CyTOF",
        }],
        "access_audit": {
            "article_source_data_anonymous_download_succeeded": True,
            "cytobank_experiment_requires_login_as_of_2026_08_06": True,
            "cytobank_file_inventory_not_used_or_inferred": True,
        },
        "design_contract": {
            "raw_CyTOF_cells": 90_066,
            "author_state_subset_cells": 29_066,
            "timepoints": ["0d", "2d", "6d", "10d", "2dW", "6dW", "10dW"],
            "released_source_contains_biological_replicate_identifier": False,
            "single_cells_are_not_biological_replicates": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    receipt = acquire(args.output_dir)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"files": 1, "size": SIZE}, sort_keys=True))


if __name__ == "__main__":
    main()
