#!/usr/bin/env python3
"""Select and acquire the frozen OpenCell raw-IF external holdout."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any

from train_hpa_raw_if_cnn import sha256


def order_key(salt: str, value: str) -> str:
    return hashlib.sha256(f"{salt}|{value}".encode()).hexdigest()


def split_annotations(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def map_annotation(value: str, protocol: dict[str, Any]) -> list[str]:
    lowered = value.casefold()
    return [
        item["canonical"] for item in protocol["coarse_ontology_in_order"]
        if any(token in lowered for token in item["casefolded_substrings"])
    ]


def verify_annotation_receipt(
    annotation_dir: Path, receipt_path: Path, protocol_path: Path
) -> Path:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt["protocol_sha256"] != sha256(protocol_path):
        raise ValueError("OpenCell annotation receipt protocol hash changed")
    for record in receipt["files"]:
        path = annotation_dir / record["local_name"]
        if path.stat().st_size != record["size_bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"OpenCell annotation source differs from receipt: {path.name}")
    return annotation_dir / "opencell-localization-annotations.csv"


def read_targets(
    csv_path: Path, protocol: dict[str, Any], hpa_genes: set[str]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    with csv_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {
            "target_name", "ensg_id", "annotations_grade_3",
            "annotations_grade_2", "annotations_grade_1",
        }
        if set(reader.fieldnames or []) != expected:
            raise ValueError(f"OpenCell wide grade schema does not resolve uniquely: {reader.fieldnames}")
        source_rows = list(reader)
    if len({row["target_name"] for row in source_rows}) != len(source_rows):
        raise ValueError("OpenCell target_name is not unique")
    targets = []
    source_label_counts: dict[str, int] = {}
    for row in source_rows:
        high_labels = split_annotations(row["annotations_grade_3"]) + split_annotations(
            row["annotations_grade_2"]
        )
        grade_1 = split_annotations(row["annotations_grade_1"])
        mapped_by_author_label = {
            label: map_annotation(label, protocol) for label in high_labels
        }
        for label in high_labels:
            source_label_counts[label] = source_label_counts.get(label, 0) + 1
        mapped_labels = sorted({
            canonical for values in mapped_by_author_label.values() for canonical in values
        })
        mapped_flags = [bool(values) for values in mapped_by_author_label.values()]
        if not high_labels:
            role = "excluded_no_grade_2_or_3"
        elif all(mapped_flags):
            role = "in_domain"
        elif not any(mapped_flags):
            role = "semantic_OOD"
        else:
            role = "excluded_mixed_mapped_unmapped"
        if row["target_name"] in hpa_genes:
            role = "excluded_shared_HPA_gene"
        targets.append({
            "target_name": row["target_name"], "ensg_id": row["ensg_id"],
            "grade_3_annotations": split_annotations(row["annotations_grade_3"]),
            "grade_2_annotations": split_annotations(row["annotations_grade_2"]),
            "grade_1_annotations": grade_1,
            "mapped_labels": mapped_labels,
            "mapped_by_author_label": mapped_by_author_label,
            "role_before_sampling": role,
        })
    audit = {
        "released_target_count": len(source_rows),
        "source_grade_2_or_3_label_counts": dict(sorted(source_label_counts.items())),
        "role_counts_before_sampling": {
            role: sum(row["role_before_sampling"] == role for row in targets)
            for role in sorted({row["role_before_sampling"] for row in targets})
        },
        "shared_HPA_gene_count": sum(
            row["role_before_sampling"] == "excluded_shared_HPA_gene" for row in targets
        ),
    }
    return targets, audit


def select_targets(targets: list[dict[str, Any]], protocol: dict[str, Any]) -> list[dict[str, Any]]:
    contract = protocol["external_sampling_contract"]
    salt = contract["selection_salt"]
    in_domain = [row for row in targets if row["role_before_sampling"] == "in_domain"]
    selected_names: set[str] = set()
    selected: list[dict[str, Any]] = []
    for item in protocol["coarse_ontology_in_order"]:
        canonical = item["canonical"]
        candidates = [row for row in in_domain if canonical in row["mapped_labels"]]
        candidates.sort(key=lambda row: order_key(salt, row["target_name"]))
        contributed = 0
        for row in candidates:
            if row["target_name"] in selected_names:
                continue
            selected_names.add(row["target_name"])
            chosen = dict(row)
            chosen["evaluation_role"] = "in_domain"
            selected.append(chosen)
            contributed += 1
            if contributed >= 48:
                break
    ood = [row for row in targets if row["role_before_sampling"] == "semantic_OOD"]
    ood.sort(key=lambda row: order_key(salt, row["target_name"]))
    for row in ood[:128]:
        chosen = dict(row)
        chosen["evaluation_role"] = "semantic_OOD"
        selected.append(chosen)
    return selected


def list_projections(target: dict[str, Any], protocol: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    source = protocol["source_contract"]
    prefix = f"{source['raw_projection_prefix']}{target['target_name']}_{target['ensg_id']}/"
    query = urllib.parse.urlencode({"list-type": "2", "prefix": prefix})
    url = f"https://{source['aws_bucket']}.s3.{source['aws_region']}.amazonaws.com/?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "Project-Aleph-research/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    root = ET.fromstring(payload)
    objects = []
    for content in root.findall("{*}Contents"):
        key = content.findtext("{*}Key") or ""
        if key.endswith(source["raw_projection_suffix"]):
            objects.append({
                "key": key,
                "size_bytes": int(content.findtext("{*}Size") or 0),
                "etag": (content.findtext("{*}ETag") or "").strip('"'),
            })
    objects.sort(key=lambda item: order_key(
        protocol["external_sampling_contract"]["selection_salt"], item["key"]
    ))
    return target, objects[:2]


def download_object(
    record: dict[str, Any], output_dir: Path, protocol: dict[str, Any]
) -> dict[str, Any]:
    source = protocol["source_contract"]
    quoted_key = urllib.parse.quote(record["key"], safe="/")
    url = f"https://{source['aws_bucket']}.s3.{source['aws_region']}.amazonaws.com/{quoted_key}"
    local_name = hashlib.sha256(record["key"].encode()).hexdigest()[:24] + ".tif"
    path = output_dir / local_name
    if not path.exists():
        error: Exception | None = None
        for attempt in range(5):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "Project-Aleph-research/1.0"})
                with urllib.request.urlopen(request, timeout=180) as response:
                    payload = response.read()
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(payload)
                break
            except Exception as exc:
                error = exc
                time.sleep(1.5 * (attempt + 1))
        else:
            raise RuntimeError(f"failed to download {record['key']}: {error}")
    if path.stat().st_size != record["size_bytes"]:
        raise ValueError(f"OpenCell projection size differs from S3 listing: {record['key']}")
    return {
        **record, "url": url, "local_name": local_name,
        "sha256": sha256(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--annotation-dir", type=Path, required=True)
    parser.add_argument("--annotation-receipt", type=Path, required=True)
    parser.add_argument("--hpa-manifest", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    csv_path = verify_annotation_receipt(
        args.annotation_dir, args.annotation_receipt, args.protocol
    )
    hpa_manifest = json.loads(args.hpa_manifest.read_text(encoding="utf-8"))
    hpa_genes = {gene for row in hpa_manifest["rows"] for gene in row["genes"]}
    targets, annotation_audit = read_targets(csv_path, protocol, hpa_genes)
    selected = select_targets(targets, protocol)
    classes = [item["canonical"] for item in protocol["coarse_ontology_in_order"]]
    pre_projection_in_domain = [
        row for row in selected if row["evaluation_role"] == "in_domain"
    ]
    pre_projection_support = {
        label: sum(label in row["mapped_labels"] for row in pre_projection_in_domain)
        for label in classes
    }
    pre_projection_evaluable = sum(
        value >= protocol["external_sampling_contract"]["minimum_targets_per_reported_class"]
        for value in pre_projection_support.values()
    )
    if pre_projection_evaluable < protocol["external_sampling_contract"]["minimum_evaluable_classes"]:
        refusal = {
            "schema": "aleph.outer_library.opencell_if_external_holdout_refusal.v1",
            "status": "refused_before_projection_download",
            "protocol_sha256": sha256(args.protocol),
            "annotation_receipt_sha256": sha256(args.annotation_receipt),
            "frozen_checkpoint_sha256": "a0b9c18516d5997ec5cf64c97a7027ab6911fc21c7b46e5f7a6c0566dc15027e",
            "annotation_audit": annotation_audit,
            "selected_in_domain_target_count_before_projection": len(pre_projection_in_domain),
            "selected_semantic_OOD_target_count_before_projection": sum(
                row["evaluation_role"] == "semantic_OOD" for row in selected
            ),
            "in_domain_class_support_before_projection": pre_projection_support,
            "evaluable_class_count": pre_projection_evaluable,
            "required_evaluable_class_count": protocol["external_sampling_contract"]["minimum_evaluable_classes"],
            "reason": "literal frozen ontology produced fewer sufficiently supported common classes than preregistered",
            "forbidden_posthoc_repairs_not_applied": [
                "er to endoplasmic reticulum",
                "membrane to plasma membrane",
                "cell_contact to cell junction",
                "underscore to space normalization",
                "lowering the six-class requirement"
            ],
            "first_attempt_projection_metadata_listed": True,
            "projection_image_content_downloaded": False,
            "annotation_labels_consumed": True,
            "may_be_reused_as_pristine": False,
            "model_prediction_executed": False,
            "aleph_authority": "none",
            "may_select_aleph_parameter": False,
        }
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(
            json.dumps(refusal, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps({
            "status": refusal["status"], "evaluable_classes": pre_projection_evaluable,
            "required": refusal["required_evaluable_class_count"],
            "projection_images_downloaded": 0,
        }, sort_keys=True))
        return
    listed: list[tuple[dict[str, Any], list[dict[str, Any]]]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(list_projections, target, protocol) for target in selected]
        for future in as_completed(futures):
            listed.append(future.result())
    rows, object_jobs, missing = [], [], []
    for target, objects in listed:
        if not objects:
            missing.append(target["target_name"])
            continue
        row = dict(target)
        row["projection_objects"] = objects
        rows.append(row)
        object_jobs.extend(objects)
    rows.sort(key=lambda row: row["target_name"])
    in_domain_rows = [row for row in rows if row["evaluation_role"] == "in_domain"]
    if len(in_domain_rows) < 100:
        raise ValueError(f"only {len(in_domain_rows)} selected OpenCell in-domain targets have projections")
    support = {
        label: sum(label in row["mapped_labels"] for row in in_domain_rows)
        for label in classes
    }
    evaluable = sum(
        value >= protocol["external_sampling_contract"]["minimum_targets_per_reported_class"]
        for value in support.values()
    )
    if evaluable < protocol["external_sampling_contract"]["minimum_evaluable_classes"]:
        raise ValueError(f"OpenCell selection yields only {evaluable} evaluable frozen classes")
    manifest = {
        "schema": "aleph.outer_library.opencell_if_external_holdout_manifest.v1",
        "protocol_sha256": sha256(args.protocol),
        "annotation_receipt_sha256": sha256(args.annotation_receipt),
        "frozen_checkpoint_sha256": "a0b9c18516d5997ec5cf64c97a7027ab6911fc21c7b46e5f7a6c0566dc15027e",
        "annotation_audit": annotation_audit,
        "selection": protocol["external_sampling_contract"],
        "selected_target_count": len(rows),
        "selected_in_domain_target_count": len(in_domain_rows),
        "selected_semantic_OOD_target_count": sum(
            row["evaluation_role"] == "semantic_OOD" for row in rows
        ),
        "selected_projection_count": len(object_jobs),
        "in_domain_class_support": support,
        "evaluable_class_count": evaluable,
        "missing_projection_targets": sorted(missing),
        "rows": rows,
        "external_labels_used_for_model_selection": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    files = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(download_object, record, args.image_dir, protocol)
            for record in object_jobs
        ]
        for future in as_completed(futures):
            files.append(future.result())
    files.sort(key=lambda item: item["key"])
    by_key = {item["key"]: item for item in files}
    for row in manifest["rows"]:
        row["projection_objects"] = [by_key[item["key"]] for item in row["projection_objects"]]
    args.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema": "aleph.outer_library.opencell_if_projection_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "annotation_receipt_sha256": sha256(args.annotation_receipt),
        "manifest_sha256": sha256(args.manifest),
        "frozen_model_commit": "e755bd9",
        "file_count": len(files),
        "total_bytes": sum(item["size_bytes"] for item in files),
        "files": files,
        "aleph_authority": "none",
    }
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "targets": len(rows), "in_domain": len(in_domain_rows),
        "semantic_OOD": manifest["selected_semantic_OOD_target_count"],
        "classes": evaluable, "projections": len(files),
        "bytes": receipt["total_bytes"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
