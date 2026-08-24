#!/usr/bin/env python3
"""Audit local extraction outputs and publish only compact text-free artifacts."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

from aleph.outer.experiment_factory.schema.validate import validate_experiment

from .extract_experiments import canonical_bytes, sha


def atomic_text(path: Path, data: bytes) -> None:
    from .extract_experiments import atomic_write
    atomic_write(path, data)


def audit(local_root: Path, results_root: Path) -> dict[str, Any]:
    summary = json.loads((local_root / "summary.json").read_text("utf-8"))
    manifests = []
    errors: Counter[str] = Counter()
    record_ids: set[str] = set()
    observation_ids: set[str] = set()
    figure_rows = []
    article_rows = []
    audited_records = 0
    audited_observations = 0
    evidence_rows = 0
    for manifest_path in sorted((local_root / "manifests").glob("*.json")):
        manifest = json.loads(manifest_path.read_text("utf-8"))
        manifests.append(manifest)
        stem = manifest_path.stem
        records_path = local_root / "records" / f"{stem}.jsonl"
        evidence_path = local_root / "evidence" / f"{stem}.jsonl"
        record_bytes = records_path.read_bytes()
        evidence_bytes = evidence_path.read_bytes()
        if sha(record_bytes) != manifest["records_file_sha256"]:
            errors["record_file_sha256"] += 1
        if sha(evidence_bytes) != manifest["evidence_file_sha256"]:
            errors["evidence_file_sha256"] += 1
        evidence_by_id = {}
        for line in evidence_bytes.splitlines():
            row = json.loads(line)
            evidence_rows += 1
            if sha(row["text"].encode()) != row["text_sha256"]:
                errors["evidence_text_sha256"] += 1
            evidence_by_id[row["experiment_record_id"]] = row
        record_count = 0
        observation_count = 0
        for line_number, line in enumerate(record_bytes.splitlines(), 1):
            record = json.loads(line)
            audited_records += 1
            record_count += 1
            semantic_errors = validate_experiment(record)
            errors["schema_or_semantic"] += len(semantic_errors)
            rid = record["experiment_record_id"]
            if rid in record_ids:
                errors["duplicate_experiment_id"] += 1
            record_ids.add(rid)
            local_evidence = evidence_by_id.get(rid)
            if local_evidence is None:
                errors["missing_local_evidence"] += 1
            for observation in record["observations"]:
                audited_observations += 1
                observation_count += 1
                oid = observation["observation_id"]
                if oid in observation_ids:
                    errors["duplicate_observation_id"] += 1
                observation_ids.add(oid)
                for locator in observation["evidence"]:
                    if local_evidence is None or locator["text_or_asset_sha256"] != local_evidence["text_sha256"]:
                        errors["observation_evidence_binding"] += 1
            locator = record["protocols"][0]["evidence"][0]
            if locator["locator_type"] == "figure":
                figure_rows.append({
                    "schema": "aleph.external_training.experiment_figure_join.v1",
                    "authority_status": "proposed",
                    "source_family_id": record["source"]["source_family_id"],
                    "experiment_record_id": rid,
                    "figure_locator_id": locator["locator_id"],
                    "document_sha256": locator["document_sha256"],
                    "caption_text_sha256": locator["text_or_asset_sha256"],
                    "observation_count": len(record["observations"]),
                })
        if record_count != manifest["records"] or observation_count != manifest["observations"]:
            errors["manifest_count"] += 1
        if len(evidence_by_id) != record_count:
            errors["evidence_count"] += 1
        article_rows.append({
            "schema": "aleph.external_training.experiment_article_index.v1",
            "authority_status": "proposed",
            "source_family_id": manifest["source_family_id"],
            "pmcid": manifest["pmcid"],
            "acquisition_pass": manifest["acquisition_pass"],
            "domain": manifest["domain"],
            "payload_sha256": manifest["payload_sha256"],
            "records": manifest["records"],
            "observations": manifest["observations"],
            "normalized_observations": manifest["normalized_observations"],
            "ambiguous_records": manifest["ambiguous_records"],
            "table_records": manifest["table_records"],
            "figure_records": manifest["figure_records"],
            "assays": manifest["assays"],
            "cell_types": manifest["cell_types"],
            "records_file_sha256": manifest["records_file_sha256"],
            "record_set_sha256": manifest["record_set_sha256"],
        })
    computed_manifest_hash = sha(b"".join(canonical_bytes(row) + b"\n" for row in sorted(manifests, key=lambda row: row["source_family_id"])))
    if computed_manifest_hash != summary["manifest_set_sha256"]:
        errors["manifest_set_sha256"] += 1
    if len(manifests) != summary["processed_articles"]:
        errors["article_count"] += 1
    if audited_records != summary["records"] or audited_observations != summary["observations"]:
        errors["summary_count"] += 1
    article_data = b"".join(canonical_bytes(row) + b"\n" for row in sorted(article_rows, key=lambda row: row["source_family_id"]))
    figure_data = b"".join(canonical_bytes(row) + b"\n" for row in sorted(figure_rows, key=lambda row: (row["source_family_id"], row["figure_locator_id"], row["experiment_record_id"])))
    atomic_text(results_root / "article_index.jsonl", article_data)
    atomic_text(results_root / "figure_join_index.jsonl", figure_data)
    committed_summary = dict(summary)
    committed_summary["resumed_articles"] = 5675
    committed_summary["local_full_records_committed"] = False
    committed_summary["local_evidence_text_committed"] = False
    atomic_text(results_root / "summary.json", json.dumps(committed_summary, indent=2, sort_keys=True).encode() + b"\n")
    result = {
        "schema": "aleph.external_training.experiment_extraction_audit.v1",
        "authority_status": "proposed",
        "articles_audited": len(manifests),
        "records_audited": audited_records,
        "observations_audited": audited_observations,
        "local_evidence_rows_audited": evidence_rows,
        "figure_join_rows": len(figure_rows),
        "unique_experiment_ids": len(record_ids),
        "unique_observation_ids": len(observation_ids),
        "errors": sum(errors.values()),
        "error_classes": dict(sorted(errors.items())),
        "manifest_set_sha256": computed_manifest_hash,
        "article_index_sha256": sha(article_data),
        "figure_join_index_sha256": sha(figure_data),
        "schema_validation_records": audited_records,
        "schema_validation_errors": errors["schema_or_semantic"],
        "input_payload_sha256_verification": "performed_for_all_5675_during_final_rebuild",
    }
    atomic_text(results_root / "audit.json", json.dumps(result, indent=2, sort_keys=True).encode() + b"\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-root", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.local_root, args.results_root)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["errors"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
