#!/usr/bin/env python3
"""Acquire independent GEO senescence panels for cross-lab validation."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

from download_sciplex3_metadata import _download as _range_download


FILES = {
    "GSE297406_hDF_raw_counts_p5_to_p15.csv.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297406/suppl/GSE297406_hDF_raw_counts_p5_to_p15.csv.gz",
        "dataset": "GSE297406",
        "lab": "Sungkyunkwan_University_Cholab",
        "role": "sealed_external_lab_four_biological_replicates_per_state",
        "size": 1_618_701,
    },
    "GSE282425_barcodes.tsv.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282425/suppl/GSE282425_barcodes.tsv.gz",
        "dataset": "GSE282425",
        "lab": "Queen_Mary_University_London_Bishop_lab",
        "role": "single_donor_2D_3D_context_shift_diagnostic",
        "size": 14_862,
    },
    "GSE282425_features.tsv.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282425/suppl/GSE282425_features.tsv.gz",
        "dataset": "GSE282425",
        "lab": "Queen_Mary_University_London_Bishop_lab",
        "role": "single_donor_2D_3D_context_shift_diagnostic",
        "size": 80_033,
    },
    "GSE282425_matrix.mtx.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282425/suppl/GSE282425_matrix.mtx.gz",
        "dataset": "GSE282425",
        "lab": "Queen_Mary_University_London_Bishop_lab",
        "role": "single_donor_2D_3D_context_shift_diagnostic",
        "size": 33_882_522,
    },
    "GSE63577_counts_rpkm_exvivo_jenage_data.xls.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63577/suppl/GSE63577_counts_rpkm_exvivo_jenage_data.xls.gz",
        "dataset": "GSE63577",
        "lab": "Leibniz_Institute_JenAge",
        "role": "fit_selection_calibration_multicellline_biological_replicates",
        "size": 14_273_182,
    },
    "GSE63577_exVivo_counts_RPKMs_II.xls.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63577/suppl/GSE63577_exVivo_counts_RPKMs_II.xls.gz",
        "dataset": "GSE63577",
        "lab": "Leibniz_Institute_JenAge",
        "role": "fit_selection_calibration_multicellline_biological_replicates",
        "size": 9_796_646,
    },
    "GSE301164_BJ_RNA_seq_raw_counts.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301164/suppl/GSE301164_BJ_RNA_seq_raw_counts.txt.gz",
        "dataset": "GSE301164",
        "lab": "Northwest_University_Yan_lab",
        "role": "post_freeze_external_confirmation_replicative_and_doxorubicin",
        "size": 2_261_178,
    },
}

EXPECTED: dict[str, dict[str, object]] = {
    "GSE297406_hDF_raw_counts_p5_to_p15.csv.gz": {
        "size_bytes": 1_618_701,
        "sha256": "807d7a66c9385535352d6231076db56cc5514c914c33a28fc807f1decbe2add4",
    },
    "GSE282425_barcodes.tsv.gz": {
        "size_bytes": 14_862,
        "sha256": "5945b59ae608ebfb2ad55b8527f50a3f254532b6a692c803ce4f49335cb0c681",
    },
    "GSE282425_features.tsv.gz": {
        "size_bytes": 80_033,
        "sha256": "8f978c80872e705d9214ca3705d74e101a0f5b305c4495efc379aeb0c361da7c",
    },
    "GSE282425_matrix.mtx.gz": {
        "size_bytes": 33_882_522,
        "sha256": "7132e6e715a176d05d04d8f749ef77e7717c555ad4a5e541a347c20da7c7e7bc",
    },
    "GSE63577_counts_rpkm_exvivo_jenage_data.xls.gz": {
        "size_bytes": 14_273_182,
        "sha256": "48b525e052a5cb17f157b875cb466cc9603b706ef6e94358e83079ef84ca5375",
    },
    "GSE63577_exVivo_counts_RPKMs_II.xls.gz": {
        "size_bytes": 9_796_646,
        "sha256": "e37fc2da3cce72cabe5c03ef9eab88cfb07945cdc231ffc1e94202a98a96a2a3",
    },
    "GSE301164_BJ_RNA_seq_raw_counts.txt.gz": {
        "size_bytes": 2_261_178,
        "sha256": "1e07430ef9287c3263ecc8e1a5d323cc469a862be8972aaa8c4b62ba8011ca9b",
    },
}


def _download(url: str, path: Path) -> None:
    partial = path.with_suffix(path.suffix + ".partial")
    request = urllib.request.Request(
        url, headers={"User-Agent": "Project-Aleph-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as output:
        while block := response.read(4 * 1024 * 1024):
            output.write(block)
    partial.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acquire(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    for name, declaration in FILES.items():
        path = output_dir / name
        if not path.exists():
            legacy_partial = path.with_suffix(path.suffix + ".partial")
            if legacy_partial.exists():
                legacy_partial.unlink()
            if int(declaration["size"]) > 2_000_000:
                _range_download(
                    str(declaration["url"]), path, int(declaration["size"]),
                    segments=32, workers=16,
                )
            else:
                _download(str(declaration["url"]), path)
        with gzip.open(path, "rb") as handle:
            prefix = handle.read(256)
            while handle.read(16 * 1024 * 1024):
                pass
        if not prefix:
            raise ValueError(f"empty gzip payload for {name}")
        record = {
            **declaration,
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
            "gzip_stream_verified": True,
        }
        if record["size_bytes"] != declaration["size"]:
            raise ValueError(f"GEO size changed for {name}")
        expected = EXPECTED.get(name)
        if expected is not None:
            for key, value in expected.items():
                if record[key] != value:
                    raise ValueError(f"locked GEO object changed for {name}: {key}")
        records[name] = record
    return {
        "schema": "aleph.outer_library.external_senescence_acquisition.v1",
        "retrieved_at": "2026-08-06",
        "files": records,
        "dataset_count": 4,
        "lab_count": 4,
        "gse297406_biological_replicates_per_state": 4,
        "gse63577_declared_samples": 48,
        "gse282425_donors": 1,
        "gse301164_biological_replicates_per_contrast": 3,
        "gse282425_contexts": ["early_proliferative_2D", "early_proliferative_3D", "deeply_senescent_2D", "deeply_senescent_3D"],
        "external_lab_labels_may_be_used_for_model_selection": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    report = acquire(args.output_dir)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({name: {"size": row["size_bytes"], "sha256": row["sha256"]} for name, row in report["files"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
