#!/usr/bin/env python3
"""Audit the external-corpus run without promoting any source or claim."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def load_json(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None
    return json.loads(path.read_text())


def find_first(root: Path, patterns: tuple[str, ...]) -> Path | None:
    for pattern in patterns:
        matches = sorted(root.glob(pattern))
        if matches:
            return matches[0]
    return None


def find_all(root: Path, patterns: tuple[str, ...]) -> list[Path]:
    """Return each matching file once, in a stable order."""
    paths: set[Path] = set()
    for pattern in patterns:
        paths.update(path for path in root.glob(pattern) if path.is_file())
    return sorted(paths)


def normalize_family(value: str) -> str:
    value = value.strip()
    if ":" not in value:
        return value.lower()
    namespace, identifier = value.split(":", 1)
    return f"{namespace.lower()}:{identifier.strip().lower()}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root

    source_files = [root / "sources" / "seed_primary_sources.jsonl"]
    source_files += sorted((root / "snapshots").glob("*.jsonl"))
    source_files += sorted((root / "acquisition" / "second_pass").glob("candidates*.jsonl"))
    rows_by_file = {str(path.relative_to(root)): load_jsonl(path) for path in source_files}
    occurrences: Counter[str] = Counter()
    canonical: dict[str, dict[str, Any]] = {}
    for rows in rows_by_file.values():
        for row in rows:
            family = normalize_family(row.get("source_family_id", ""))
            if not family:
                continue
            occurrences[family] += 1
            merged = canonical.setdefault(family, dict(row))
            merged["is_open_access_provider_flag"] = bool(
                merged.get("is_open_access_provider_flag") or row.get("is_open_access_provider_flag")
            )
            for key in ("doi", "pmcid", "pmid", "year"):
                if not merged.get(key) and row.get(key):
                    merged[key] = row[key]

    candidates = list(canonical.values())
    review_rows = load_jsonl(root / "review" / "article_quality_screen.jsonl")
    review_decisions = Counter(row.get("quality_decision", "UNKNOWN") for row in review_rows)

    acquisition_paths = []
    acquisition_rows: list[dict[str, Any]] = []
    for path in find_all(root, ("acquisition/**/*.jsonl",)):
        rows = load_jsonl(path)
        if rows and all("status" in row and "pmcid" in row for row in rows):
            acquisition_paths.append(path)
            acquisition_rows.extend(rows)
    acquisition_status = Counter(row.get("status", "UNKNOWN") for row in acquisition_rows)
    acquisition_summaries = {
        str(path.relative_to(root)): load_json(path)
        for path in find_all(root, ("acquisition/**/*summary*.json",))
    }
    rag_summaries = {
        str(path.relative_to(root)): load_json(path)
        for path in find_all(
            root,
            ("ragtagcag/results/*.json", "ragtagcag/**/combined_ingest*.json"),
        )
    }
    extraction_summaries = {
        str(path.relative_to(root)): load_json(path)
        for path in find_all(root, ("extraction/**/results/*.json",))
    }
    model_metrics = load_json(root / "model" / "metrics.json")
    oa_content_metrics = load_json(root / "model" / "oa_content" / "metrics.json")
    model_metric_files = {
        str(path.relative_to(root)): load_json(path)
        for path in find_all(root, ("model/metrics.json", "model/**/metrics.json"))
    }
    oa_quality_summaries = {
        str(path.relative_to(root)): load_json(path)
        for path in find_all(root, ("review_oa/results/summary.json", "review_oa/**/results/summary.json"))
    }
    model_compact = None
    if model_metrics:
        model_compact = {
            "authority": model_metrics.get("authority"),
            "source_families": model_metrics.get("input", {}).get("unique_source_families"),
            "splits": model_metrics.get("splits"),
            "neural_macro_f1": model_metrics.get("tagging", {}).get("neural_macro_f1"),
            "nearest_centroid_macro_f1": model_metrics.get("tagging", {}).get("nearest_centroid_macro_f1"),
            "majority_macro_f1": model_metrics.get("tagging", {}).get("majority_macro_f1"),
            "retrieval_neural": model_metrics.get("retrieval", {}).get("neural_hidden"),
            "retrieval_baseline": model_metrics.get("retrieval", {}).get("fixed_hashed_baseline"),
            "runtime": model_metrics.get("runtime"),
        }

    report: dict[str, Any] = {
        "authority_status": "proposed",
        "candidate_files": {name: len(rows) for name, rows in rows_by_file.items()},
        "candidate_unique_source_families": len(candidates),
        "candidate_cross_file_duplicates": sum(count - 1 for count in occurrences.values()),
        "candidate_oa_provider_flag": sum(bool(row.get("is_open_access_provider_flag")) for row in candidates),
        "candidate_with_pmcid": sum(bool(row.get("pmcid")) for row in candidates),
        "candidate_with_doi": sum(bool(row.get("doi")) for row in candidates),
        "candidate_years": dict(sorted(Counter(str(row.get("year", "UNKNOWN")) for row in candidates).items())),
        "review_records": len(review_rows),
        "review_decisions": dict(sorted(review_decisions.items())),
        "acquisition_receipts": [str(path.relative_to(root)) for path in acquisition_paths],
        "acquisition_records": len(acquisition_rows),
        "acquisition_unique_pmcids": len({row.get("pmcid") for row in acquisition_rows if row.get("pmcid")}),
        "acquisition_status": dict(sorted(acquisition_status.items())),
        "acquisition_summaries": acquisition_summaries,
        "ragtagcag_summaries": rag_summaries,
        "extraction_summaries": extraction_summaries,
        "model_metrics": model_compact,
        "oa_content_model_metrics": oa_content_metrics,
        "model_metric_files": model_metric_files,
        "oa_quality_summaries": oa_quality_summaries,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
