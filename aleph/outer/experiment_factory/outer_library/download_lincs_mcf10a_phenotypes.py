#!/usr/bin/env python3
"""Download pinned author-released LINCS MCF10A phenotype and metadata tables."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request


PAPER_DOI = "10.1038/s42003-022-03975-9"
BASE = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs42003-022-03975-9/MediaObjects"
)
FILES = (
    (
        "42003_2022_3975_MOESM3_ESM.xlsx",
        25788,
        "a1055652224ba1f0e3f1df1e767637ae417ffb5544024131043d0907b133ca26",
        "Supplementary_Data_1_author_well_level_IF_phenotypes",
    ),
    (
        "42003_2022_3975_MOESM25_ESM.xlsx",
        89683,
        "0f74f51269c94dd8bde300ddfbcd20c1195672fc1be380efbd16d93351196814",
        "Supplementary_Data_23_author_sample_metadata",
    ),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for name, expected_size, expected_sha, role in FILES:
        destination = output_dir / name
        if not destination.is_file() or destination.stat().st_size != expected_size:
            temporary = destination.with_suffix(destination.suffix + ".part")
            request = urllib.request.Request(
                f"{BASE}/{name}",
                headers={"User-Agent": "Project-Aleph-research/1.0"},
            )
            with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as handle:
                while block := response.read(1024 * 1024):
                    handle.write(block)
            if temporary.stat().st_size != expected_size:
                raise ValueError(f"LINCS MCF10A size changed for {name}")
            temporary.replace(destination)
        observed_sha = _sha256(destination)
        if observed_sha != expected_sha:
            raise ValueError(f"LINCS MCF10A digest changed for {name}")
        records.append({
            "local_name": name, "url": f"{BASE}/{name}", "role": role,
            "size": expected_size, "sha256": observed_sha,
        })
    return {
        "schema": "aleph.outer_library.lincs_mcf10a_acquisition.v1",
        "source": {
            "paper_doi": PAPER_DOI,
            "paper_url": "https://www.nature.com/articles/s42003-022-03975-9",
            "pmcid": "PMC9546880",
            "provider": "Springer Nature article supplementary data",
            "license_id": "cc-by-4.0",
        },
        "files": records,
        "design_contract": {
            "C1_and_C2_are_separate_OHSU_collection_periods": True,
            "each_collection_has_at_least_three_biological_replicates": True,
            "phenotype_rows_are_author_well_level_summaries": True,
            "paper_reports_eight_biological_replicates": True,
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
    print(json.dumps({"files": len(receipt["files"])}, sort_keys=True))


if __name__ == "__main__":
    main()
