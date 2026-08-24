#!/usr/bin/env python3
"""Resolve frozen OpenCell targets and projection objects for development."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any

from train_hpa_raw_if_cnn import sha256


def normalized(value: str) -> str:
    return re.sub(r"[-\s]+", "_", value.casefold().strip())


def map_label(value: str, ontology: list[dict[str, Any]]) -> list[str]:
    exact = normalized(value)
    exact_matches = [
        item["canonical"] for item in ontology
        if any(normalized(token) == exact for token in item["author_tokens"])
    ]
    if exact_matches:
        return exact_matches
    lowered = value.casefold()
    return [
        item["canonical"] for item in ontology
        if any(token.casefold() in lowered for token in item["author_tokens"])
    ]


def split_labels(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def split_role(salt: str, target: str) -> str:
    bucket = int(hashlib.sha256(f"{salt}|{target}".encode()).hexdigest(), 16) % 20
    if bucket <= 2:
        return "calibration"
    if bucket <= 5:
        return "selection"
    return "train"


def order_key(salt: str, value: str) -> str:
    return hashlib.sha256(f"{salt}|{value}".encode()).hexdigest()


def list_projections(
    row: dict[str, Any], provider: dict[str, Any], salt: str, limit: int,
) -> tuple[str, list[dict[str, Any]], str | None]:
    prefix = f"{provider['raw_projection_prefix']}{row['target_name']}_{row['ensg_id']}/"
    query = urllib.parse.urlencode({"list-type": "2", "prefix": prefix})
    url = (
        f"https://{provider['aws_bucket']}.s3.{provider['aws_region']}"
        f".amazonaws.com/?{query}"
    )
    error: Exception | None = None
    for attempt in range(5):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "Project-Aleph-research/1.0"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            root = ET.fromstring(payload)
            objects = []
            for content in root.findall("{*}Contents"):
                key = content.findtext("{*}Key") or ""
                if key.endswith(provider["raw_projection_suffix"]):
                    objects.append({
                        "key": key,
                        "size_bytes": int(content.findtext("{*}Size") or 0),
                        "etag": (content.findtext("{*}ETag") or "").strip('"'),
                    })
            objects.sort(key=lambda item: order_key(salt, item["key"]))
            return row["target_name"], objects[:limit], None
        except Exception as exc:
            error = exc
            time.sleep(1.5 * (attempt + 1))
    return row["target_name"], [], f"{type(error).__name__}: {error}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--annotation-csv", type=Path, required=True)
    parser.add_argument("--annotation-receipt", type=Path, required=True)
    parser.add_argument("--first-refusal", type=Path, required=True)
    parser.add_argument("--hpa-manifest", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    frozen = protocol["frozen_inputs"]
    if sha256(args.annotation_receipt) != frozen[
        "opencell_annotation_receipt_sha256"
    ]:
        raise ValueError("OpenCell annotation receipt changed after protocol freeze")
    if sha256(args.first_refusal) != frozen["opencell_first_attempt_refusal_sha256"]:
        raise ValueError("OpenCell first refusal changed after protocol freeze")
    receipt = json.loads(args.annotation_receipt.read_text(encoding="utf-8"))
    annotation_record = next(
        item for item in receipt["files"]
        if item["local_name"] == "opencell-localization-annotations.csv"
    )
    if (
        args.annotation_csv.stat().st_size != annotation_record["size_bytes"]
        or sha256(args.annotation_csv) != annotation_record["sha256"]
    ):
        raise ValueError("OpenCell annotation CSV differs from frozen receipt")
    hpa = json.loads(args.hpa_manifest.read_text(encoding="utf-8"))
    hpa_genes = {gene for row in hpa["rows"] for gene in row["genes"]}
    with args.annotation_csv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {
            "target_name", "ensg_id", "annotations_grade_3",
            "annotations_grade_2", "annotations_grade_1",
        }
        if set(reader.fieldnames or []) != expected:
            raise ValueError("OpenCell annotation columns changed")
        source_rows = list(reader)
    if len(source_rows) != 1310:
        raise ValueError("OpenCell released target count changed")
    ontology = protocol["coarse_ontology_in_order"]
    salt = protocol["opencell_development_contract"]["split_salt"]
    rows, excluded = [], []
    for source in source_rows:
        high = split_labels(source["annotations_grade_3"]) + split_labels(
            source["annotations_grade_2"]
        )
        mapped_by_label = {label: map_label(label, ontology) for label in high}
        flags = [bool(values) for values in mapped_by_label.values()]
        if not high:
            role = "excluded_no_grade_2_or_3"
        elif all(flags):
            role = "in_domain"
        elif not any(flags):
            role = "semantic_OOD"
        else:
            role = "excluded_mixed_mapped_unmapped"
        row = {
            "target_name": source["target_name"],
            "ensg_id": source["ensg_id"],
            "grade_3_annotations": split_labels(source["annotations_grade_3"]),
            "grade_2_annotations": split_labels(source["annotations_grade_2"]),
            "grade_1_annotations": split_labels(source["annotations_grade_1"]),
            "mapped_by_author_label": mapped_by_label,
            "mapped_labels": sorted({
                canonical for values in mapped_by_label.values()
                for canonical in values
            }),
            "development_role": role,
            "exact_HPA_gene_overlap": source["target_name"] in hpa_genes,
        }
        if role == "in_domain":
            row["split"] = split_role(salt, source["target_name"])
            rows.append(row)
        elif role == "semantic_OOD":
            row["split"] = "semantic_OOD"
            rows.append(row)
        else:
            excluded.append(row)
    rows.sort(key=lambda row: row["target_name"])

    object_map, errors = {}, {}
    provider = protocol["provider_roles"]["OpenCell"]
    limit = int(protocol["opencell_development_contract"][
        "maximum_projections_per_target"
    ])
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(list_projections, row, provider, salt, limit) for row in rows
        ]
        for future in as_completed(futures):
            target, objects, error = future.result()
            object_map[target] = objects
            if error is not None:
                errors[target] = error
    for row in rows:
        row["projection_objects"] = object_map[row["target_name"]]
    admitted = [row for row in rows if row["projection_objects"]]
    missing = [row["target_name"] for row in rows if not row["projection_objects"]]
    split_counts = Counter(row["split"] for row in admitted)
    support = {
        item["canonical"]: {
            split: sum(
                row["split"] == split and item["canonical"] in row["mapped_labels"]
                for row in admitted
            )
            for split in ("train", "selection", "calibration")
        }
        for item in ontology
    }
    contract = protocol["opencell_development_contract"]
    supported = sum(
        values["train"] >= contract["minimum_targets_per_supported_class"]
        for values in support.values()
    )
    gates = {
        "train_target_minimum": split_counts["train"] >= contract[
            "minimum_train_targets"
        ],
        "selection_target_minimum": split_counts["selection"] >= contract[
            "minimum_selection_targets"
        ],
        "calibration_target_minimum": split_counts["calibration"] >= contract[
            "minimum_calibration_targets"
        ],
        "supported_class_minimum": supported >= contract["minimum_supported_classes"],
        "projection_listing_errors_zero": not errors,
    }
    manifest = {
        "schema": "aleph.outer_library.opencell_multiprovider_development_manifest.v1",
        "protocol_sha256": sha256(args.protocol),
        "annotation_receipt_sha256": sha256(args.annotation_receipt),
        "first_refusal_sha256": sha256(args.first_refusal),
        "hpa_manifest_sha256": sha256(args.hpa_manifest),
        "rows": admitted,
        "excluded_rows": excluded,
        "missing_projection_targets": missing,
        "listing_errors": errors,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    report = {
        "schema": "aleph.outer_library.opencell_multiprovider_resolution.v1",
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "released_target_count": len(source_rows),
        "admitted_target_count": len(admitted),
        "admitted_projection_count": sum(
            len(row["projection_objects"]) for row in admitted
        ),
        "admitted_projection_bytes": sum(
            item["size_bytes"] for row in admitted
            for item in row["projection_objects"]
        ),
        "development_role_counts_before_projection": dict(sorted(Counter(
            row["development_role"] for row in rows + excluded
        ).items())),
        "split_target_counts": dict(sorted(split_counts.items())),
        "class_split_support": support,
        "supported_train_class_count": supported,
        "exact_HPA_gene_overlap_count": sum(
            row["exact_HPA_gene_overlap"] for row in admitted
        ),
        "missing_projection_target_count": len(missing),
        "listing_error_count": len(errors),
        "development_gates": gates,
        "all_development_gates_passed": all(gates.values()),
        "OpenCell_is_pristine_confirmation": False,
        "OpenCell_role": "consumed_development_only",
        "Allen_metadata_or_pixels_accessed": 0,
        "model_predictions_made": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "admitted_targets": len(admitted),
        "projections": report["admitted_projection_count"],
        "bytes": report["admitted_projection_bytes"],
        "splits": report["split_target_counts"],
        "supported_classes": supported,
        "gates": gates,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
