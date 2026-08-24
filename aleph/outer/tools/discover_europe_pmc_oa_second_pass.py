#!/usr/bin/env python3
"""Discover a stratified, non-overlapping second pass of Europe PMC OA papers."""

from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator


PMCID_RE = re.compile(r"PMC[1-9][0-9]*\Z")
USER_AGENT = "Project-Aleph-OA-Discovery/2.0 (polite research client; Europe PMC REST)"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_name, path)
    except BaseException:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def family_id(record: dict[str, Any]) -> str | None:
    doi = str(record.get("doi") or "").strip().lower()
    if doi:
        return f"doi:{doi}"
    pmcid = str(record.get("pmcid") or "").strip().upper()
    if pmcid:
        return f"pmcid:{pmcid}"
    pmid = str(record.get("pmid") or "").strip()
    return f"pmid:{pmid}" if pmid else None


def load_exclusions(paths: list[Path]) -> tuple[set[str], set[str]]:
    families: set[str] = set()
    pmcids: set[str] = set()
    for path in paths:
        with path.open("r", encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
                fid = record.get("source_family_id")
                if fid:
                    families.add(str(fid).strip().lower())
                pmcid = str(record.get("pmcid") or "").strip().upper()
                if PMCID_RE.fullmatch(pmcid):
                    pmcids.add(pmcid)
    return families, pmcids


def request_page(
    endpoint: str, query: str, cursor: str, page_size: int, retries: int = 3
) -> dict[str, Any]:
    params = urllib.parse.urlencode(
        {
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": page_size,
            "cursorMark": cursor,
        }
    )
    request = urllib.request.Request(
        f"{endpoint}?{params}",
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    last_error: BaseException | None = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, ConnectionError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(2**attempt)
    assert last_error is not None
    raise last_error


def pages(
    endpoint: str, query: str, page_size: int, max_pages: int, delay: float
) -> Iterator[tuple[int, dict[str, Any]]]:
    cursor = "*"
    for page_number in range(1, max_pages + 1):
        payload = request_page(endpoint, query, cursor, page_size)
        yield page_number, payload
        next_cursor = payload.get("nextCursorMark")
        records = (payload.get("resultList") or {}).get("result") or []
        if not records or not next_cursor or next_cursor == cursor:
            return
        cursor = next_cursor
        if delay:
            time.sleep(delay)


def is_research_article(record: dict[str, Any]) -> bool:
    types = (record.get("pubTypeList") or {}).get("pubType") or []
    return any(str(item).strip().lower() == "research-article" for item in types)


def compact(
    record: dict[str, Any], domain: str, stratum: str, query: str, rank: int, snapshot_date: str
) -> dict[str, Any]:
    journal_info = record.get("journalInfo") or {}
    journal = journal_info.get("journal") or {}
    return {
        "authority_status": "proposed",
        "authors": record.get("authorString"),
        "cited_by_count_at_snapshot": int(record.get("citedByCount") or 0),
        "discovered_at": utc_now(),
        "discovered_via": "Europe PMC REST API",
        "discovery_domain": domain,
        "discovery_query": query,
        "doi": record.get("doi"),
        "in_europe_pmc": record.get("inEPMC") == "Y",
        "is_open_access_provider_flag": record.get("isOpenAccess") == "Y",
        "journal": journal.get("title"),
        "pmcid": str(record.get("pmcid") or "").upper(),
        "pmid": record.get("pmid"),
        "provider_license": record.get("license"),
        "publication_types": (record.get("pubTypeList") or {}).get("pubType") or [],
        "quality_status": "pending_article_level_review",
        "snapshot_date": snapshot_date,
        "source_family_id": family_id(record),
        "stratum_rank": rank,
        "title": record.get("title"),
        "year": int(record["pubYear"]),
        "year_stratum": stratum,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--exclude", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    excluded_families, excluded_pmcids = load_exclusions(args.exclude)
    global_seen_families = set(excluded_families)
    global_seen_pmcids = set(excluded_pmcids)
    strata: dict[str, list[dict[str, Any]]] = {}
    diagnostics: dict[str, dict[str, Any]] = {}

    for domain_spec in config["domains"]:
        for year_spec in config["year_strata"]:
            key = f'{domain_spec["domain"]}::{year_spec["name"]}'
            query = (
                'SRC:MED AND OPEN_ACCESS:Y AND PUB_TYPE:"research-article" '
                f'AND FIRST_PDATE:[{year_spec["minimum_date"]} TO {year_spec["maximum_date"]}] '
                f'AND {domain_spec["query"]}'
            )
            accepted: list[dict[str, Any]] = []
            pages_fetched = 0
            hit_count = 0
            rejected = Counter()
            for page_number, payload in pages(
                config["provider_url"],
                query,
                int(config["page_size"]),
                int(config["max_pages_per_stratum"]),
                float(config["request_delay_seconds"]),
            ):
                pages_fetched = page_number
                hit_count = int(payload.get("hitCount") or 0)
                records = (payload.get("resultList") or {}).get("result") or []
                for record in records:
                    pmcid = str(record.get("pmcid") or "").strip().upper()
                    fid = family_id(record)
                    try:
                        year = int(record.get("pubYear") or 0)
                    except (TypeError, ValueError):
                        rejected["bad_year"] += 1
                        continue
                    if record.get("isOpenAccess") != "Y":
                        rejected["not_provider_oa"] += 1
                    elif not PMCID_RE.fullmatch(pmcid):
                        rejected["no_valid_pmcid"] += 1
                    elif not fid:
                        rejected["no_source_family"] += 1
                    elif not is_research_article(record):
                        rejected["not_research_article"] += 1
                    elif not (int(year_spec["minimum_date"][:4]) <= year <= int(year_spec["maximum_date"][:4])):
                        rejected["outside_year_stratum"] += 1
                    elif pmcid in global_seen_pmcids or fid.lower() in global_seen_families:
                        rejected["excluded_or_cross_query_duplicate"] += 1
                    else:
                        global_seen_pmcids.add(pmcid)
                        global_seen_families.add(fid.lower())
                        accepted.append(
                            compact(
                                record,
                                domain_spec["domain"],
                                year_spec["name"],
                                query,
                                len(accepted) + 1,
                                config["snapshot_date"],
                            )
                        )
                        if len(accepted) >= int(config["target_per_stratum"]):
                            break
                if len(accepted) >= int(config["target_per_stratum"]):
                    break
            strata[key] = accepted
            diagnostics[key] = {
                "accepted": len(accepted),
                "hit_count": hit_count,
                "pages_fetched": pages_fetched,
                "rejected": dict(sorted(rejected.items())),
            }
            print(json.dumps({"stratum": key, **diagnostics[key]}, sort_keys=True), flush=True)
            if float(config["request_delay_seconds"]):
                time.sleep(float(config["request_delay_seconds"]))

    # Balanced round-robin prevents early domains from consuming the global cap.
    selected: list[dict[str, Any]] = []
    ordered_keys = sorted(strata)
    rank = 0
    while len(selected) < int(config["global_cap"]):
        added = 0
        for key in ordered_keys:
            if rank < len(strata[key]):
                selected.append(strata[key][rank])
                added += 1
                if len(selected) >= int(config["global_cap"]):
                    break
        if not added:
            break
        rank += 1

    output = b"".join(
        (json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n").encode("utf-8")
        for record in selected
    )
    summary = {
        "schema": "aleph.external_training.oa_second_pass_discovery_summary.v1",
        "generated_at": utc_now(),
        "snapshot_date": config["snapshot_date"],
        "global_cap": int(config["global_cap"]),
        "selected_count": len(selected),
        "excluded_source_family_count": len(excluded_families),
        "excluded_pmcid_count": len(excluded_pmcids),
        "domain_counts": dict(sorted(Counter(r["discovery_domain"] for r in selected).items())),
        "year_stratum_counts": dict(sorted(Counter(r["year_stratum"] for r in selected).items())),
        "strata": diagnostics,
    }
    atomic_write(args.output, output)
    atomic_write(args.summary, (json.dumps(summary, sort_keys=True) + "\n").encode("utf-8"))
    print(json.dumps(summary, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
