#!/usr/bin/env python3
"""Audit the official GSE250041 proliferating/senescent CITE-seq design."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path


FILES = {
    "GSE250041_family.soft.gz": {
        "size": 3_920,
        "sha256": "906358053dbc7517fcd95e52b8fe8080318ef197123cc5a53806cd74e66e461b",
    },
    "GSE250041_Proliferating_barcodes.tsv.gz": {
        "size": 44_407,
        "sha256": "d4a9d7945ee894797ed2e79c2cc94cf439fd0b1fa85d0018429a385d7278ce25",
    },
    "GSE250041_Proliferating_features.tsv.gz": {
        "size": 333_557,
        "sha256": "e034835514a64c40721b95aff253ff30d2b7cf4fba73de3327c4813540885e37",
    },
    "GSE250041_Senescent_barcodes.tsv.gz": {
        "size": 26_244,
        "sha256": "d075bad15bef08eb6e2282f68bfd1d5d29ac41e4c643a9054512dded29cbb65e",
    },
    "GSE250041_Senescent_features.tsv.gz": {
        "size": 333_557,
        "sha256": "8ead549950cc2640cdac915ac3de5c57f92a284db0b82bad94bd1e1b6a87ac83",
    },
}
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250041"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _line_count(path: Path) -> int:
    with gzip.open(path, "rb") as handle:
        return sum(block.count(b"\n") for block in iter(lambda: handle.read(1024 * 1024), b""))


def _features(path: Path) -> tuple[list[tuple[str, str, str]], Counter[str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = [tuple(row) for row in csv.reader(handle, delimiter="\t")]
    if any(len(row) != 3 for row in rows):
        raise ValueError(f"GSE250041 feature table width changed: {path.name}")
    return rows, Counter(row[2] for row in rows)


def audit(source_dir: Path, matrix_receipt: Path | None = None) -> dict[str, object]:
    evidence: dict[str, object] = {}
    for name, expected in FILES.items():
        path = source_dir / name
        actual = {"size": path.stat().st_size, "sha256": _sha256(path)}
        if actual != expected:
            raise ValueError(f"GSE250041 source changed for {name}: {actual}")
        evidence[name] = {
            **actual,
            "url": (
                f"{BASE}/soft/{name}" if name.endswith("family.soft.gz")
                else f"{BASE}/suppl/{name}"
            ),
        }

    proliferating_features, proliferating_types = _features(
        source_dir / "GSE250041_Proliferating_features.tsv.gz"
    )
    senescent_features, senescent_types = _features(
        source_dir / "GSE250041_Senescent_features.tsv.gz"
    )
    if proliferating_features != senescent_features:
        raise ValueError("GSE250041 conditions use different feature axes")
    expected_types = Counter({"Gene Expression": 36_601, "Antibody Capture": 8})
    if proliferating_types != expected_types or senescent_types != expected_types:
        raise ValueError("GSE250041 RNA/ADT feature counts changed")

    barcode_counts = {
        "untreated_proliferating": _line_count(
            source_dir / "GSE250041_Proliferating_barcodes.tsv.gz"
        ),
        "IR_10_Gy_day_10_senescent": _line_count(
            source_dir / "GSE250041_Senescent_barcodes.tsv.gz"
        ),
    }
    if barcode_counts != {
        "untreated_proliferating": 8_664,
        "IR_10_Gy_day_10_senescent": 4_949,
    }:
        raise ValueError(f"GSE250041 barcode counts changed: {barcode_counts}")

    with gzip.open(
        source_dir / "GSE250041_family.soft.gz", "rt", encoding="utf-8"
    ) as handle:
        soft = handle.read()
    required_text = (
        "!Series_geo_accession = GSE250041",
        "!Series_status = Public on Mar 23 2026",
        "!Series_pubmed_id = 39143693",
        "cell line: WI-38",
        "cell type: human fetal lung fibroblasts",
        "treatment: Untreated proliferating controls",
        "treatment: IR-induced senescent (10 days)",
        "cultured for an additional 10 days",
        "NIA-IRP, NIH",
    )
    missing = [text for text in required_text if text not in soft]
    if missing:
        raise ValueError(f"GSE250041 SOFT metadata changed: missing {missing}")

    matrix_evidence: dict[str, object] = {}
    if matrix_receipt is not None:
        receipt = json.loads(matrix_receipt.read_text(encoding="utf-8"))
        if receipt.get("schema") != "aleph.outer_library.gse250041_matrix_receipt.v1":
            raise ValueError("wrong GSE250041 matrix receipt schema")
        if receipt.get("accession") != "GSE250041":
            raise ValueError("wrong GSE250041 matrix receipt accession")
        expected_shapes = {
            "GSE250041_Proliferating_matrix.mtx.gz": (36_609, 8_664),
            "GSE250041_Senescent_matrix.mtx.gz": (36_609, 4_949),
        }
        for name, expected_shape in expected_shapes.items():
            item = receipt["files"].get(name)
            if not item or item.get("gzip_stream_verified") is not True:
                raise ValueError(f"unverified GSE250041 matrix receipt item: {name}")
            path = source_dir / name
            actual = {"size_bytes": path.stat().st_size, "sha256": _sha256(path)}
            if actual != {key: item[key] for key in ("size_bytes", "sha256")}:
                raise ValueError(f"GSE250041 matrix receipt does not match local file: {name}")
            shape = item.get("shape", {})
            if (shape.get("features"), shape.get("barcodes")) != expected_shape:
                raise ValueError(f"GSE250041 matrix dimensions changed: {name}")
            matrix_evidence[name] = item

    return {
        "schema": "aleph.outer_library.gse250041_senescence_design.v1",
        "accession": "GSE250041",
        "pubmed_id": "39143693",
        "public_date": "2026-03-23",
        "provider_lab": "NIA_IRP_NIH_Laboratory_of_Genetics_and_Genomics",
        "organism": "Homo sapiens",
        "cell_line": "WI-38",
        "cell_type": "human_fetal_lung_fibroblast",
        "conditions": ["untreated_proliferating", "IR_10_Gy_day_10_senescent"],
        "author_state_label": "ionizing_radiation_induced_senescence",
        "barcode_counts": barcode_counts,
        "total_cell_barcodes": sum(barcode_counts.values()),
        "feature_type_counts": dict(expected_types),
        "feature_axis_identical_across_conditions": True,
        "verified_files": evidence,
        "matrix_files": matrix_evidence or {
            "GSE250041_Proliferating_matrix.mtx.gz": "pending",
            "GSE250041_Senescent_matrix.mtx.gz": "pending",
        },
        "matrix_training_ready": bool(matrix_evidence),
        "biological_replicate_count_per_condition": 1,
        "single_cells_are_not_independent_biological_replicates": True,
        "independent_dataset_from_sciplex3": True,
        "independent_lab_from_sciplex3": True,
        "eligible_role": "direct_senescence_state_representation_training_only",
        "external_validation_eligible": False,
        "reason_not_external_validation": (
            "one pooled CITE-seq capture per condition; cells do not replace biological replicates"
        ),
        "license_status": "no_explicit_reuse_license_on_GEO_record",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--matrix-receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.source_dir, args.matrix_receipt)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        key: report[key] for key in (
            "total_cell_barcodes", "feature_type_counts", "matrix_training_ready"
        )
    }, sort_keys=True))


if __name__ == "__main__":
    main()
