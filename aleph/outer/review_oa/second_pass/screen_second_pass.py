#!/usr/bin/env python3
"""Run unchanged OA recoverability rules on a stable second-pass extraction."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE.parent / "screen_oa.py"
SPEC = importlib.util.spec_from_file_location("aleph_oa_screen_base", BASE_PATH)
BASE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(BASE)

EXPECTED_ARTICLES = 2974
VERSION = "1"


def validate_extraction(derived: Path) -> tuple[list[tuple[dict[str, Any], Path]], dict[str, Any]]:
    summary_path = derived / "final_summary.json"
    if not summary_path.is_file():
        audit_path = derived / "audit.json"
        repeat_path = derived / "repeat_summary.json"
        if not (audit_path.is_file() and repeat_path.is_file()):
            raise ValueError("second-pass stable audit and repeat summary are not present")
        audit = json.loads(audit_path.read_text("utf-8"))
        if audit.get("audit_errors") != 0 or not audit.get("manifest_set_hash_equal") or not audit.get("summaries_equal_excluding_resume_count"):
            raise ValueError("second-pass extraction audit is not stable")
        summary_path = repeat_path
    summary = json.loads(summary_path.read_text("utf-8"))
    if summary.get("failures") != 0 or summary.get("xml_parsed") != EXPECTED_ARTICLES:
        raise ValueError("second-pass extraction is not final and complete")
    entries = []
    manifests = []
    for path in sorted((derived / "articles").glob("*.manifest.json")):
        manifest = json.loads(path.read_text("utf-8"))
        chunks_path = path.with_name(path.name.replace(".manifest.json", ".chunks.jsonl"))
        if not chunks_path.is_file() or BASE.sha(chunks_path.read_bytes()) != manifest["chunks_file_sha256"]:
            raise ValueError(f"chunk digest mismatch: {chunks_path}")
        entries.append((manifest, chunks_path))
        manifests.append(manifest)
    set_hash = BASE.sha(b"".join(BASE.canonical_bytes(item) + b"\n" for item in sorted(manifests, key=lambda value: value["source_family_id"])))
    if len(entries) != EXPECTED_ARTICLES or set_hash != summary.get("article_manifest_set_sha256"):
        raise ValueError("second-pass manifest set does not reproduce final summary")
    return entries, summary


def load_first_records(first_local_root: Path) -> list[dict[str, Any]]:
    records = [json.loads(path.read_text("utf-8")) for path in sorted((first_local_root / "records").glob("*.json"))]
    if len(records) != 2701 or len({item["source_family_id"] for item in records}) != 2701:
        raise ValueError("first-pass local record set is unavailable or incomplete")
    return records


def write_queue(records: list[dict[str, Any]], path: Path) -> None:
    queued = sorted((item for item in records if item["decision"] != "hold"), key=lambda item: (-item["priority_score"], item["source_family_id"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["rank", "source_family_id", "pmcid", "decision", "exception_reason", "priority_score", "missing_review_fields", "observation_locator_count", "cell_state_tags", "modality_tags", "mechanics_tags"])
        for rank, record in enumerate(queued, 1):
            writer.writerow([
                rank, record["source_family_id"], record["pmcid"], record["decision"], record["exception_reason"], record["priority_score"],
                ";".join(record["missing_review_fields"]), len(record["observation_locators"]),
                ";".join(record["tags"]["cell_type_state"]), ";".join(record["tags"]["measurement_modality"]), ";".join(record["tags"]["mechanics_observable"]),
            ])
        temporary = Path(handle.name)
    os.replace(temporary, path)


def leakage_groups(first: list[dict[str, Any]], second: list[dict[str, Any]]) -> tuple[dict[str, Any], dict[str, int]]:
    members: dict[str, list[dict[str, str]]] = defaultdict(list)
    for pass_name, records in (("first", first), ("second", second)):
        for record in records:
            for group in record["leakage"]["dataset_groups"]:
                members[group["group_hash"]].append({"pass": pass_name, "source_family_id": record["source_family_id"]})
    groups = []
    cross = within_first = within_second = 0
    for group_hash, values in sorted(members.items()):
        unique = sorted({(item["pass"], item["source_family_id"]) for item in values})
        if len(unique) < 2:
            continue
        passes = {item[0] for item in unique}
        cross += passes == {"first", "second"}
        within_first += sum(item[0] == "first" for item in unique) > 1
        within_second += sum(item[0] == "second" for item in unique) > 1
        groups.append({
            "dataset_group_hash": group_hash,
            "source_family_count": len(unique),
            "members": [{"pass": pass_name, "source_family_id": family} for pass_name, family in unique],
            "uncertainty": "public accession co-membership; dataset identity and processing lineage require manual confirmation",
        })
    counts = {"duplicate_groups": len(groups), "cross_pass_groups": cross, "within_first_groups": within_first, "within_second_groups": within_second}
    return {"schema": "aleph.external_training.combined_oa_dataset_leakage_groups.v1", "authority_status": "proposed", "groups": groups}, counts


def process(derived: Path, candidates: Path, prior_review: Path | None, first_local_root: Path, local_root: Path, output_root: Path) -> dict[str, Any]:
    entries, extraction_summary = validate_extraction(derived)
    metadata = {item["source_family_id"]: item for item in BASE.iter_jsonl(candidates)}
    prior = BASE.load_prior_review(prior_review)
    records_dir = local_root / "records"
    records_dir.mkdir(parents=True, exist_ok=True)
    records = []
    compact_records = []
    for manifest, chunks_path in entries:
        family = manifest["source_family_id"]
        record = BASE.screen_article(manifest, chunks_path, metadata.get(family, {}), prior)
        data = BASE.canonical_bytes(record) + b"\n"
        record_sha = BASE.sha(data)
        BASE.atomic_write(records_dir / (BASE.sha(family.encode())[:24] + ".json"), data)
        records.append(record)
        compact_records.append(BASE.compact(record, record_sha))
    compact_records.sort(key=lambda item: item["source_family_id"])
    compact_data = b"".join(BASE.canonical_bytes(item) + b"\n" for item in compact_records)
    BASE.atomic_write(output_root / "article_screen_index.jsonl", compact_data)
    queue_path = output_root / "manual_review_queue.csv"
    write_queue(records, queue_path)

    first_records = load_first_records(first_local_root)
    first_families = {item["source_family_id"] for item in first_records}
    second_families = {item["source_family_id"] for item in records}
    overlap = first_families & second_families
    if overlap:
        raise ValueError(f"source-family overlap across passes: {len(overlap)}")
    combined_leakage, combined_leakage_counts = leakage_groups(first_records, records)
    combined_leakage_data = BASE.canonical_bytes(combined_leakage) + b"\n"
    BASE.atomic_write(output_root / "combined_dataset_leakage_groups.json", combined_leakage_data)
    within_leakage, within_counts = leakage_groups([], records)
    within_leakage_data = BASE.canonical_bytes({
        "schema": "aleph.external_training.second_pass_oa_dataset_leakage_groups.v1",
        "authority_status": "proposed",
        "groups": within_leakage["groups"],
    }) + b"\n"
    BASE.atomic_write(output_root / "dataset_leakage_groups.json", within_leakage_data)

    decisions = Counter(item["decision"] for item in records)
    article_types = Counter(item["article_type"]["status"] for item in records)
    missing = Counter(field for item in records for field in item["missing_review_fields"])
    signals = {name: sum(item["signals"][name]["status"] == "recoverable_signal" for item in records) for name in BASE.SIGNALS}
    corrections = Counter(item["correction_retraction"]["status"] for item in records)
    first_summary = json.loads((HERE.parent / "results" / "summary.json").read_text("utf-8"))
    combined_decisions = Counter(first_summary["decisions"])
    combined_decisions.update(decisions)
    summary = {
        "schema": "aleph.external_training.oa_recoverability_second_pass_summary.v1",
        "authority_status": "proposed",
        "screen_version": VERSION,
        "articles": len(records),
        "decisions": dict(sorted(decisions.items())),
        "automatic_tier_a_promotions": 0,
        "article_types": dict(sorted(article_types.items())),
        "signal_recoverability_articles": dict(sorted(signals.items())),
        "missing_review_fields": dict(sorted(missing.items())),
        "correction_retraction_status": dict(sorted(corrections.items())),
        "articles_with_observation_locator": sum(bool(item["observation_locators"]) for item in records),
        "articles_with_dataset_accession": sum(bool(item["leakage"]["dataset_groups"]) for item in records),
        "within_second_pass_leakage": within_counts,
        "combined": {
            "articles": len(first_records) + len(records),
            "source_family_overlap": 0,
            "decisions": dict(sorted(combined_decisions.items())),
            "automatic_tier_a_promotions": 0,
            "dataset_leakage": combined_leakage_counts,
        },
        "local_record_set_sha256": BASE.sha(b"".join((item["source_family_id"] + ":" + item["record_sha256"] + "\n").encode() for item in compact_records)),
        "compact_index_sha256": BASE.sha(compact_data),
        "manual_queue_sha256": BASE.sha(queue_path.read_bytes()),
        "within_leakage_sha256": BASE.sha(within_leakage_data),
        "combined_leakage_sha256": BASE.sha(combined_leakage_data),
        "input_manifest_set_sha256": extraction_summary["article_manifest_set_sha256"],
        "first_pass_compact_index_sha256": first_summary["compact_index_sha256"],
        "limitations": [
            "regex and vocabulary hits indicate recoverability only, not adequacy or truth",
            "correction/retraction relations are absent for most second-pass candidates",
            "independent replication and affiliation-resolved laboratory groups require manual review",
            "dataset-accession co-membership is a leakage lead, not proof of identical processed data",
            "no decision is Tier A or training authorization",
        ],
    }
    summary_data = json.dumps(summary, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    BASE.atomic_write(output_root / "summary.json", summary_data)
    sums = {path.name: BASE.sha(path.read_bytes()) for path in sorted(output_root.iterdir()) if path.is_file() and path.name != "SHA256SUMS.json"}
    BASE.atomic_write(output_root / "SHA256SUMS.json", json.dumps(sums, indent=2, sort_keys=True).encode("utf-8") + b"\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derived-root", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--prior-review", type=Path)
    parser.add_argument("--first-local-root", type=Path, required=True)
    parser.add_argument("--local-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(process(args.derived_root, args.candidates, args.prior_review, args.first_local_root, args.local_root, args.output_root), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
