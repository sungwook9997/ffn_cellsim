#!/usr/bin/env python3
"""Audit sci-Plex3 cell metadata and define non-leaking split hierarchy."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path


def audit(
    pdata_path: Path,
    annotations_path: Path,
    matrix_receipt: dict[str, object] | None = None,
) -> dict[str, object]:
    cell_types: Counter[str] = Counter()
    replicates: Counter[str] = Counter()
    time_points: Counter[str] = Counter()
    treatments: Counter[str] = Counter()
    doses: Counter[str] = Counter()
    groups: Counter[tuple[str, str, str, str, str]] = Counter()
    row_count = 0
    mismatches = 0
    unassigned_context = 0
    with (
        gzip.open(pdata_path, "rt", encoding="utf-8", newline="") as pdata,
        gzip.open(annotations_path, "rt", encoding="utf-8", newline="") as annotations,
    ):
        reader = csv.reader(pdata, delimiter=" ", skipinitialspace=True)
        header = next(reader)
        if len(header) != 39 or header[:4] != ["cell", "sample", "Size_Factor", "n.umi"]:
            raise ValueError("sci-Plex3 pData schema changed")
        fields = ["row_name"] + header
        for values, annotation_line in zip(reader, annotations, strict=True):
            if len(values) != len(fields):
                raise ValueError(f"sci-Plex3 row width changed at row {row_count}")
            row = dict(zip(fields, values, strict=True))
            annotation = annotation_line.rstrip("\n").split("\t")
            if len(annotation) != 2 or annotation != [row["row_name"], row["sample"]]:
                mismatches += 1
            cell_types[row["cell_type"]] += 1
            replicates[row["replicate"]] += 1
            time_points[row["time_point"]] += 1
            treatments[row["treatment"]] += 1
            doses[row["dose"]] += 1
            context = (
                row["cell_type"], row["replicate"], row["time_point"],
                row["treatment"], row["dose"],
            )
            if any(value in {"", "NA"} for value in context):
                unassigned_context += 1
            else:
                groups[context] += 1
            row_count += 1
    if row_count != 799_317 or mismatches:
        raise ValueError(f"sci-Plex3 alignment failed: rows={row_count}, mismatches={mismatches}")
    group_sizes = sorted(groups.values())
    matrix_ready = False
    matrix_evidence: dict[str, object] | None = None
    if matrix_receipt is not None:
        required = {
            "schema": "aleph.outer_library.sciplex3_matrix_receipt.v1",
            "accession": "GSM4150378",
            "size_bytes": 3_096_224_606,
            "gzip_stream_verified": True,
            "aleph_authority": "none",
        }
        mismatched = {
            key: (matrix_receipt.get(key), expected)
            for key, expected in required.items()
            if matrix_receipt.get(key) != expected
        }
        if mismatched:
            raise ValueError(f"invalid sci-Plex3 matrix receipt: {mismatched}")
        if int(matrix_receipt.get("decompressed_line_count", 0)) < 1_000_000:
            raise ValueError("sci-Plex3 matrix receipt has too few sparse rows")
        digest = str(matrix_receipt.get("sha256", ""))
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
            raise ValueError("sci-Plex3 matrix receipt lacks a valid SHA-256")
        matrix_ready = True
        matrix_evidence = {
            "size_bytes": matrix_receipt["size_bytes"],
            "sha256": digest,
            "decompressed_line_count": matrix_receipt["decompressed_line_count"],
            "gzip_stream_verified": True,
        }
    report = {
        "schema": "aleph.outer_library.sciplex3_design_audit.v1",
        "accession": "GSM4150378",
        "cell_count": row_count,
        "cell_type_counts": dict(sorted(cell_types.items())),
        "replicate_counts": dict(sorted(replicates.items())),
        "time_point_counts": dict(sorted(time_points.items(), key=lambda item: (item[0] == "NA", item[0]))),
        "treatment_count_including_vehicle": len([key for key in treatments if key]),
        "dose_counts": dict(sorted(doses.items(), key=lambda item: (item[0] == "NA", item[0]))),
        "assigned_context_cell_count": row_count - unassigned_context,
        "unassigned_context_cell_count": unassigned_context,
        "unassigned_context_cells_may_enter_supervised_training": False,
        "mandatory_context_group_count": len(groups),
        "mandatory_context_group_size": {
            "minimum": group_sizes[0],
            "median": group_sizes[len(group_sizes) // 2],
            "maximum": group_sizes[-1],
        },
        "cell_annotation_alignment_mismatches": mismatches,
        "required_split_group": [
            "provider_lab", "cell_type", "replicate", "time_point_hours",
            "treatment", "dose",
        ],
        "replicate_holdout_is_same_lab_not_external_lab": True,
        "expression_matrix_acquired": matrix_ready,
        "training_ready": matrix_ready,
        "reason_not_training_ready": (
            None if matrix_ready else "verified expression matrix has not yet been acquired"
        ),
        "license_status": "no_explicit_reuse_license_on_GEO_record",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }
    if matrix_evidence is not None:
        report["matrix_evidence"] = matrix_evidence
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdata", type=Path, required=True)
    parser.add_argument("--cell-annotations", type=Path, required=True)
    parser.add_argument("--matrix-receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    matrix_receipt = (
        json.loads(args.matrix_receipt.read_text(encoding="utf-8"))
        if args.matrix_receipt else None
    )
    report = audit(args.pdata, args.cell_annotations, matrix_receipt)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "cell_count", "mandatory_context_group_count", "training_ready"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
