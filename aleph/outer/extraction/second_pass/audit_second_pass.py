#!/usr/bin/env python3
"""Audit deterministic second-pass extraction without emitting article text."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def comparable_summary(summary: dict[str, Any]) -> dict[str, Any]:
    result = dict(summary)
    result.pop("resumed_articles", None)
    return result


def audit(first_path: Path, repeat_path: Path, derived_root: Path) -> dict[str, Any]:
    first = json.loads(first_path.read_text("utf-8"))
    repeat = json.loads(repeat_path.read_text("utf-8"))
    summaries_equal = comparable_summary(first) == comparable_summary(repeat)
    manifest_hash_equal = first.get("article_manifest_set_sha256") == repeat.get("article_manifest_set_sha256")

    article_root = derived_root / "articles"
    manifest_paths = sorted(article_root.glob("*.manifest.json"))
    chunk_paths = sorted(article_root.glob("*.chunks.jsonl"))
    chunk_files_by_name = {path.name.removesuffix(".chunks.jsonl"): path for path in chunk_paths}
    counters = Counter()

    for manifest_path in manifest_paths:
        manifest = json.loads(manifest_path.read_text("utf-8"))
        key = manifest_path.name.removesuffix(".manifest.json")
        chunk_path = chunk_files_by_name.get(key)
        if chunk_path is None:
            counters["missing_chunk_files"] += 1
            continue
        chunk_bytes = chunk_path.read_bytes()
        if sha256(chunk_bytes) != manifest.get("chunks_file_sha256"):
            counters["bad_chunk_file_hashes"] += 1
        local_records = 0
        for line in chunk_bytes.splitlines():
            if not line:
                continue
            local_records += 1
            record = json.loads(line)
            if sha256(record["text"].encode("utf-8")) != record.get("text_sha256"):
                counters["bad_text_hashes"] += 1
            basis = dict(record)
            expected_chunk_hash = basis.pop("chunk_sha256", "")
            basis.pop("text", None)
            if sha256(canonical_bytes(basis)) != expected_chunk_hash:
                counters["bad_chunk_hashes"] += 1
            if record.get("authority_status") != "proposed":
                counters["non_proposed_records"] += 1
        if local_records != manifest.get("chunks"):
            counters["manifest_chunk_count_mismatches"] += 1
        counters["records"] += local_records

    counters["orphan_chunk_files"] = len(set(chunk_files_by_name) - {path.name.removesuffix(".manifest.json") for path in manifest_paths})
    errors = sum(value for key, value in counters.items() if key not in {"records"})
    result = {
        "schema": "aleph.external_training.second_pass_extraction_audit.v1",
        "authority_status": "proposed",
        "summaries_equal_excluding_resume_count": summaries_equal,
        "manifest_set_hash_equal": manifest_hash_equal,
        "article_manifest_set_sha256": repeat.get("article_manifest_set_sha256"),
        "article_manifests": len(manifest_paths),
        "article_chunk_files": len(chunk_paths),
        "chunk_records": counters["records"],
        "audit_errors": errors,
        "error_counts": {key: counters[key] for key in sorted(counters) if key != "records"},
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-summary", type=Path, required=True)
    parser.add_argument("--repeat-summary", type=Path, required=True)
    parser.add_argument("--derived-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.first_summary, args.repeat_summary, args.derived_root)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if result["audit_errors"] == 0 and result["summaries_equal_excluding_resume_count"] and result["manifest_set_hash_equal"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

