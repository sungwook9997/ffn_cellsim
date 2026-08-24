#!/usr/bin/env python3
"""Evaluate the frozen GSE72502 donor-paired IFN direction endpoint."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def sample_design(path: Path) -> list[dict[str, str]]:
    wanted: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("!Sample_description", "!Sample_geo_accession", "!Sample_title")):
                row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                values = [value.strip('"') for value in row[1:]]
                if row[0] != "!Sample_description" or all(value.startswith("SAM") for value in values):
                    wanted[row[0]] = values
            elif line.startswith("!Sample_characteristics_ch1"):
                row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                values = [value.strip('"') for value in row[1:]]
                if values and values[0].startswith("donor id:"):
                    wanted["donor"] = [value.split(":", 1)[1].strip() for value in values]
                if values and values[0].startswith("treatment:"):
                    wanted["treatment"] = [value.split(":", 1)[1].strip() for value in values]
    required = ["!Sample_description", "!Sample_geo_accession", "donor", "treatment"]
    if any(key not in wanted for key in required):
        raise ValueError(f"incomplete official sample design: {sorted(wanted)}")
    return [
        {"sample": wanted["!Sample_description"][index], "geo_accession": wanted["!Sample_geo_accession"][index],
         "donor": wanted["donor"][index], "treatment": wanted["treatment"][index]}
        for index in range(len(wanted["donor"]))
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--soffice", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    signature = protocol["IFN_state_contract"]["frozen_signature_genes"]
    source = args.data_dir / "GSE72502_PBMC_RPKMs.xls.gz"
    design_source = args.data_dir / "GSE72502_series_matrix.txt.gz"
    design = sample_design(design_source)

    with tempfile.TemporaryDirectory(prefix="aleph-gse72502-") as temporary:
        temporary_path = Path(temporary)
        workbook = temporary_path / "GSE72502_PBMC_RPKMs.xls"
        with gzip.open(source, "rb") as source_handle, workbook.open("wb") as output_handle:
            shutil.copyfileobj(source_handle, output_handle)
        subprocess.run(
            [str(args.soffice), "--headless", "--convert-to", "csv", "--outdir", str(temporary_path), str(workbook)],
            check=True, capture_output=True, text=True,
        )
        csv_path = temporary_path / "GSE72502_PBMC_RPKMs.csv"
        genes = {}
        with csv_path.open("rt", encoding="utf-8-sig", newline="") as handle:
            next(handle)  # released workbook has a time-group row above the field names
            for row in csv.DictReader(handle):
                if row["symbol"] in signature:
                    genes[row["symbol"]] = row
    missing = sorted(set(signature) - set(genes))
    if missing:
        raise ValueError(f"frozen IFN signature genes missing: {missing}")

    by_donor: dict[str, dict[str, dict[str, str]]] = {}
    for sample in design:
        by_donor.setdefault(sample["donor"], {})[sample["treatment"]] = sample
    differences = {}
    gene_differences = {}
    for donor, samples in sorted(by_donor.items()):
        control = samples["none"]["sample"]
        stimulated = samples["interferon-alpha"]["sample"]
        values = {
            gene: math.log1p(float(genes[gene][stimulated])) - math.log1p(float(genes[gene][control]))
            for gene in signature
        }
        gene_differences[donor] = values
        differences[donor] = sum(values.values()) / len(values)
    passed = len(differences) == 3 and all(value > 0 for value in differences.values())
    report = {
        "schema": "aleph.outer_library.pbmc_ifn_direction.v1",
        "protocol_sha256": sha256(args.protocol),
        "source_sha256": sha256(source), "sample_design_sha256": sha256(design_source),
        "transformation": "per-gene log1p(RPKM), then unweighted mean of the 14 frozen gene differences",
        "sample_design": design, "donor_mean_signature_difference_IFN_minus_control": differences,
        "per_gene_log1p_differences": gene_differences,
        "all_three_paired_donor_differences_positive": passed,
        "numeric_scale_transfer_permitted": False,
        "role": "direction_only_cross_lab_confirmation_for_pre_sweep_context; not a single-cell classifier holdout",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"donor_differences": differences, "passed": passed}, sort_keys=True))


if __name__ == "__main__":
    main()
