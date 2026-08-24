#!/usr/bin/env python3
"""Build a metadata-only, deduplicated Europe PMC candidate snapshot."""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


def fetch(endpoint: str, query: str, page_size: int) -> list[dict[str, Any]]:
    params = urllib.parse.urlencode(
        {"query": query, "format": "json", "resultType": "core", "pageSize": page_size}
    )
    request = urllib.request.Request(
        f"{endpoint}?{params}",
        headers={"User-Agent": "Project-Aleph-external-corpus/1.0"},
    )
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = json.load(response)
            return payload["resultList"]["result"]
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def family_id(record: dict[str, Any]) -> str | None:
    if record.get("doi"):
        return "doi:" + record["doi"].strip().lower()
    if record.get("pmcid"):
        return "pmcid:" + record["pmcid"].strip().upper()
    if record.get("pmid"):
        return "pmid:" + record["pmid"].strip()
    return None


def compact(record: dict[str, Any], domain: str, snapshot_date: str) -> dict[str, Any]:
    journal = record.get("journalInfo") or {}
    pub_types = (record.get("pubTypeList") or {}).get("pubType") or []
    return {
        "source_family_id": family_id(record),
        "title": record.get("title"),
        "year": int(record["pubYear"]),
        "doi": record.get("doi"),
        "pmid": record.get("pmid"),
        "pmcid": record.get("pmcid"),
        "journal": journal.get("journal", {}).get("title"),
        "authors": record.get("authorString"),
        "publication_types": pub_types,
        "cited_by_count_at_snapshot": int(record.get("citedByCount") or 0),
        "is_open_access_provider_flag": record.get("isOpenAccess") == "Y",
        "in_europe_pmc": record.get("inEPMC") == "Y",
        "provider_license": record.get("license"),
        "discovery_domain": domain,
        "discovered_via": "Europe PMC REST API",
        "snapshot_date": snapshot_date,
        "journal_metric": None,
        "journal_metric_status": "pending_licensed_jcr_lookup",
        "quality_status": "pending_article_level_review",
        "authority_status": "proposed",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--journal-metric-output", type=Path)
    parser.add_argument(
        "--exclude",
        type=Path,
        action="append",
        default=[],
        help="JSONL source records whose source_family_id values must be skipped",
    )
    parser.add_argument("--candidate-multiplier", type=int, default=4)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    seen: set[str] = set()
    for excluded_path in args.exclude:
        for line in excluded_path.read_text().splitlines():
            if line:
                excluded = json.loads(line)
                if excluded.get("source_family_id"):
                    seen.add(excluded["source_family_id"])
    selected: list[dict[str, Any]] = []
    domain_counts: dict[str, int] = {}
    for spec in config["domains"]:
        target = int(spec["target"])
        minimum_date = spec.get("minimum_date", f'{config["minimum_year"]}-01-01')
        maximum_date = spec.get("maximum_date", config.get("maximum_date", config["snapshot_date"]))
        minimum_year = int(minimum_date[:4])
        maximum_year = int(maximum_date[:4])
        query = (
            f'SRC:MED AND FIRST_PDATE:[{minimum_date} TO {maximum_date}] AND {spec["query"]}'
        )
        records = fetch(config["provider_url"], query, min(1000, target * args.candidate_multiplier))
        accepted = 0
        for record in records:
            fid = family_id(record)
            try:
                year = int(record.get("pubYear") or 0)
            except ValueError:
                continue
            types = (record.get("pubTypeList") or {}).get("pubType") or []
            if not fid or fid in seen or year < minimum_year or year > maximum_year:
                continue
            if not any(kind.lower() in {"research-article", "journal article"} for kind in types):
                continue
            seen.add(fid)
            selected.append(compact(record, spec["domain"], config["snapshot_date"]))
            accepted += 1
            if accepted == target:
                break
        domain_counts[spec["domain"]] = accepted
        if accepted != target:
            raise RuntimeError(f'{spec["domain"]}: selected {accepted}, expected {target}')
    if len(selected) != config["target_total"]:
        raise RuntimeError(f'selected {len(selected)}, expected {config["target_total"]}')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n" for row in selected))
    if args.journal_metric_output:
        journal_rows: dict[str, list[dict[str, Any]]] = {}
        for row in selected:
            if row["journal"]:
                journal_rows.setdefault(row["journal"], []).append(row)
        args.journal_metric_output.parent.mkdir(parents=True, exist_ok=True)
        with args.journal_metric_output.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "journal", "candidate_count", "earliest_year", "latest_year",
                    "metric_name", "metric_value", "metric_year", "metric_source",
                    "verified_at", "status",
                ],
            )
            writer.writeheader()
            for journal, journal_sources in sorted(journal_rows.items(), key=lambda item: (-len(item[1]), item[0])):
                writer.writerow(
                    {
                        "journal": journal,
                        "candidate_count": len(journal_sources),
                        "earliest_year": min(row["year"] for row in journal_sources),
                        "latest_year": max(row["year"] for row in journal_sources),
                        "metric_name": "JIF",
                        "metric_value": "",
                        "metric_year": "",
                        "metric_source": "",
                        "verified_at": "",
                        "status": "pending_licensed_jcr_lookup",
                    }
                )
    print(json.dumps({"total": len(selected), "domains": domain_counts}, sort_keys=True))


if __name__ == "__main__":
    main()
