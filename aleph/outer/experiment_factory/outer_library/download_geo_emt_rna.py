#!/usr/bin/env python3
"""Acquire three biologically replicated GEO A549 EMT expression studies."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request


FILES = {
    "GSE69667_RNASEQ_TPM_rep1.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69667/suppl/GSE69667_RNASEQ_TPM_rep1.txt.gz",
        "size": 1_082_455,
        "sha256": "11fbc9e45f8ea7102c779b711dbf7052c9caa3ed4b572a95fc5f7231b2c33d5a",
    },
    "GSE69667_RNASEQ_TPM_rep2.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69667/suppl/GSE69667_RNASEQ_TPM_rep2.txt.gz",
        "size": 1_058_187,
        "sha256": "9245710a49be62a743ee563add88ce03ad0f734e9c69900a6f9df16be4266c82",
    },
    "GSE69667_family.soft.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69667/soft/GSE69667_family.soft.gz",
        "size": 2_180,
        "sha256": "2256977b39fed4dbb6809a8be62893cb0c80d473b87c599d237209397a2f9d0d",
    },
    "GSE49644_series_matrix.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49644/matrix/GSE49644_series_matrix.txt.gz",
        "size": 2_036_171,
        "sha256": "0e8afb9d4f2aa4d4d39869b08153f5d49211572d556dcd3d0e1f8e5cd9af623c",
    },
    "GSE49644_annotated_results.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49644/suppl/GSE49644_annotated_results.txt.gz",
        "size": 2_204_155,
        "sha256": "3914a982f3428660cb358fac3b135e970aa228af02dbe68c18b4115e48950884",
    },
    "GSE49644_family.soft.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49644/soft/GSE49644_family.soft.gz",
        "size": 15_002_098,
        "sha256": "29afe137525afaaea50a1c211cfa70bc35e3afcb0dd669209687793a43eab4d1",
    },
    "GSE17708_series_matrix.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17708/matrix/GSE17708_series_matrix.txt.gz",
        "size": 4_604_333,
        "sha256": "968e1c6d27551d6713d610780620d0391060a23ffc912285300f1608f0ccb3fa",
    },
    "GSE42373_series_matrix.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42373/matrix/GSE42373_series_matrix.txt.gz",
        "size": 1_615_527,
        "sha256": "5e2ddea1061ddb2b758dcd004099c9f6c00c80e10deec3d4c179dba559864bf2",
    },
    "GSE125369_RAW.tar": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125369/suppl/GSE125369_RAW.tar",
        "size": 3_737_600,
        "sha256": "c0cec52852335486838a25e1caa14c4be174e17789d365299ab62aecfbef07a6",
    },
    "GSE125369_family.soft.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125369/soft/GSE125369_family.soft.gz",
        "size": 2_026,
        "sha256": "e63ce08516044557e2d262317702d3daf1c725958f3027d1f46985aa36d4c6ca",
    },
    "GPL570.annot.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
        "size": 8_471_521,
        "sha256": "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def acquire(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    verified = {}
    for name, contract in FILES.items():
        path = output_dir / name
        if not path.is_file() or path.stat().st_size != contract["size"]:
            temporary = path.with_suffix(path.suffix + ".part")
            request = urllib.request.Request(
                str(contract["url"]), headers={"User-Agent": "Project-Aleph-research/1.0"}
            )
            with urllib.request.urlopen(request, timeout=300) as response, temporary.open("wb") as handle:
                while block := response.read(1024 * 1024):
                    handle.write(block)
            if temporary.stat().st_size != contract["size"]:
                raise ValueError(f"GEO source size changed: {name}")
            temporary.replace(path)
        digest = _sha256(path)
        if digest != contract["sha256"]:
            raise ValueError(f"GEO source digest changed: {name}")
        verified[name] = {**contract, "observed_sha256": digest}
    return {
        "schema": "aleph.outer_library.geo_emt_rna_acquisition.v1",
        "retrieved_at": "2026-08-06",
        "provider": "NCBI Gene Expression Omnibus official FTP",
        "accessions": [
            "GSE17708", "GSE42373", "GSE125369", "GSE69667", "GSE49644", "GPL570"
        ],
        "files": verified,
        "design_contract": {
            "studies": 5,
            "independent_lab_groups": 5,
            "biological_samples": 72,
            "single_samples_or_expression_features_may_not_be_resampled_as_biological_studies": True,
        },
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    receipt = acquire(args.output_dir)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"files": len(receipt["files"]), "studies": 5}, sort_keys=True))


if __name__ == "__main__":
    main()
