#!/usr/bin/env python3
"""Normalize a predeclared fibrotic marker panel from GSE226374 bulk RNA-seq."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

from schema import Observation


SOURCE_SHA256 = "c5970662289626ab2194c6cdf15f47017a1ecce5031f14d9e58c8d2e3b191cbd"
MARKERS = ("ACTA2", "CCN2", "COL1A1", "FN1", "THBS1")


def _condition(column: str) -> tuple[str, str]:
    tgf = column.startswith("TGFBplus")
    celastrol = "Celastrol_500nM" in column
    mapping = {
        (False, False): ("vehicle_control", "control"),
        (True, False): ("TGF_beta_1", "profibrotic_TGF_beta_1_stimulated"),
        (False, True): ("celastrol_only", "YAP_pathway_antagonized_control"),
        (True, True): (
            "TGF_beta_1_plus_celastrol",
            "YAP_pathway_antagonized_TGF_beta_1",
        ),
    }
    return mapping[(tgf, celastrol)]


def build(path: Path) -> list[Observation]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError(f"official GSE226374 source hash mismatch: {digest}")
    rows: list[Observation] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("GSE226374 normalized-count header is missing")
        sample_columns = reader.fieldnames[9:]
        for source_row, record in enumerate(reader, start=2):
            gene = record["Gene Symbol"]
            if gene not in MARKERS:
                continue
            for column in sample_columns:
                condition, state = _condition(column)
                locator = f"tsv:row={source_row};gene={gene};column={column}"
                identifier = hashlib.sha256(
                    f"{SOURCE_SHA256}:{locator}".encode()
                ).hexdigest()[:20]
                observation = Observation(
                    observation_id=f"gse226374:{identifier}",
                    dataset_id="geo-gse226374-tgfb-celastrol-rnaseq",
                    lab_group="leask-lab-gse226374",
                    provider="NCBI Gene Expression Omnibus",
                    accession="GSE226374",
                    license_id="geo_public_no_explicit_license",
                    source_sha256=SOURCE_SHA256,
                    source_locator=locator,
                    sample_id=f"GSE226374:{column}",
                    biological_replicate=(
                        "author_replicate_2" if "_2_S" in column
                        else "author_replicate_1"
                    ),
                    technical_replicate=column.rsplit("_", 1)[-1],
                    condition=condition,
                    cell_type="CCD-1077Sk human dermal fibroblasts",
                    cell_state=state,
                    modality="bulk_RNA_seq_derived",
                    observable=f"{gene}_normalized_expression",
                    value=float(record[column]),
                    unit="DESeq2_median_ratio_normalized_count",
                )
                observation.validate()
                rows.append(observation)
    if {row.observable.removesuffix("_normalized_expression") for row in rows} != set(MARKERS):
        raise ValueError("GSE226374 marker panel is incomplete")
    if len(rows) != len(MARKERS) * 8:
        raise ValueError(f"expected 40 marker observations, got {len(rows)}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--normalized-counts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(args.normalized_counts)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({
        "observations": len(rows),
        "samples": len({row.sample_id for row in rows}),
        "markers": list(MARKERS),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
