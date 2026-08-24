#!/usr/bin/env python
"""Build the proposed first-batch external-corpus quality review.

Reads only the frozen 2026-08-05 inputs.  Network reads are limited to the
Europe PMC REST service.  Missing or ambiguous evidence is emitted as UNKNOWN.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
CORPUS = ROOT / "corpus" / "external_training"
REVIEW = CORPUS / "review"
SNAPSHOT = CORPUS / "snapshots" / "europe_pmc_2026-08-05.jsonl"
SEEDS = CORPUS / "sources" / "seed_primary_sources.jsonl"
JOURNAL_QUEUE = CORPUS / "queues" / "journal_metrics_2026-08-05.csv"
VERIFIED_AT = "2026-08-05"
JCR_DATASET_NOTE = "Clarivate JCR 2025 metrics; dataset updated 2026-06-17"


JIF = {
    "Scientific reports": ("Scientific Reports", "4.9", "SCI%20REP-UK"),
    "Biophysical journal": ("BIOPHYSICAL JOURNAL", "3.7", "BIOPHYS%20J"),
    "Nature communications": ("Nature Communications", "18.1", "NAT%20COMMUN"),
    "International journal of molecular sciences": ("INTERNATIONAL JOURNAL OF MOLECULAR SCIENCES", "5.6", "INT%20J%20MOL%20SCI"),
    "Acta biomaterialia": ("Acta Biomaterialia", "10.4", "ACTA%20BIOMATER"),
    "Advanced science (Weinheim, Baden-Wurttemberg, Germany)": ("Advanced Science", "14.1", "ADV%20SCI"),
    "Proceedings of the National Academy of Sciences of the United States of America": ("PROCEEDINGS OF THE NATIONAL ACADEMY OF SCIENCES OF THE UNITED STATES OF AMERICA", "9.5", "P%20NATL%20ACAD%20SCI%20USA"),
    "PLoS computational biology": ("PLoS Computational Biology", "3.7", "PLOS%20COMPUT%20BIOL"),
    "Frontiers in cell and developmental biology": ("Frontiers in Cell and Developmental Biology", "5.3", "FRONT%20CELL%20DEV%20BIOL"),
    "Cell reports": ("Cell Reports", "7.7", "CELL%20REP"),
    "Molecular biology of the cell": ("MOLECULAR BIOLOGY OF THE CELL", "3.0", "MOL%20BIOL%20CELL"),
    "iScience": ("iScience", "4.5", "ISCIENCE"),
    "Communications biology": ("Communications Biology", "5.8", "COMMUN%20BIOL"),
    "Journal of cell science": ("JOURNAL OF CELL SCIENCE", "3.9", "J%20CELL%20SCI"),
    "Soft matter": ("Soft Matter", "2.9", "SOFT%20MATTER"),
}


SIGNALS = {
    "biological_system": r"\b(human|mouse|mice|rat|zebrafish|drosophila|yeast|cell line|primary cells?|organoid|tissue)\b",
    "cell_state": r"\b(mitotic|interphase|quiescent|differentiated|stem cell|activated|apoptotic|migrating|confluent)\b",
    "perturbation": r"\b(knockdown|knockout|deplet(?:ed|ion)|inhibit(?:or|ed|ion)|treated with|transfect(?:ed|ion)|overexpress(?:ed|ion)|stimulat(?:ed|ion))\b",
    "geometry": r"\b(2d|3d|two-dimensional|three-dimensional|spherical|cylindrical|micropattern|confined|geometry)\b",
    "substrate": r"\b(substrate|fibronectin|collagen|laminin|hydrogel|polyacrylamide|glass|polystyrene|matrix stiffness)\b",
    "sample_size": r"(?:\bn\s*[=:]\s*\d+\b|\b\d+\s+(?:cells?|animals?|mice|samples?|donors?|experiments?)\b)",
    "biological_independent_replicates": r"\b(?:biological replicates?|independent experiments?|independent biological|separate experiments?)\b",
    "calibration": r"\b(?:calibrat(?:e|ed|ion)|standard curve|pixel size|bead size|instrument response|force constant)\b",
    "units": r"(?:\b(?:nm|µm|um|mm|pa|kpa|mpa|pn|nn|mn|hz|s|min|h)\b|°c)",
    "uncertainty": r"(?:±|\b(?:standard deviation|standard error|s\.d\.|s\.e\.m\.|confidence interval|95% ci|error bars?)\b)",
    "exclusions": r"\b(?:exclud(?:e|ed|ing|sion)|outlier|pre-specified|blinded|randomized)\b",
    "source_data": r"\b(?:data availability|source data|raw data|deposited|accession|available at|available from)\b",
}


def get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "Project-Aleph-proposed-quality-review/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read()


def short_excerpt(text: str, match: re.Match[str]) -> str:
    start = max(0, text.rfind(".", 0, match.start()) + 1)
    end = text.find(".", match.end())
    if end < 0:
        end = min(len(text), match.end() + 120)
    words = re.sub(r"\s+", " ", text[start : end + 1]).strip().split()
    return " ".join(words[:20])


def fulltext_signals(pmcid: str | None) -> tuple[dict, str, str]:
    empty = {k: {"status": "UNKNOWN", "locator": "UNKNOWN", "evidence_excerpt": "UNKNOWN"} for k in SIGNALS}
    if not pmcid:
        return empty, "metadata_only_no_pmcid", "UNKNOWN"
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{urllib.parse.quote(pmcid)}/fullTextXML"
    try:
        root = ET.fromstring(get(url))
    except Exception as exc:  # network/provider failures remain explicit
        return empty, "full_text_access_failed", f"{type(exc).__name__}: {str(exc)[:120]}"
    sections: list[tuple[str, str]] = []
    abstract = " ".join((root.findtext(".//abstract") or "").split())
    if not abstract:
        abstract = " ".join("".join(root.find(".//abstract").itertext()).split()) if root.find(".//abstract") is not None else ""
    if abstract:
        sections.append(("abstract", abstract))
    for index, sec in enumerate(root.findall(".//sec")):
        title_node = sec.find("title")
        title = " ".join("".join(title_node.itertext()).split()) if title_node is not None else f"untitled-{index}"
        text = " ".join("".join(sec.itertext()).split())
        if text:
            sections.append((title[:120], text))
    found = {}
    for field, pattern in SIGNALS.items():
        hit = None
        for title, text in sections:
            match = re.search(pattern, text, flags=re.I)
            if match:
                hit = {"status": "RECOVERED", "locator": f"PMC full text section: {title}", "evidence_excerpt": short_excerpt(text, match)}
                break
        found[field] = hit or {"status": "UNKNOWN", "locator": "UNKNOWN", "evidence_excerpt": "UNKNOWN"}
    return found, "open_full_text_reviewed", "NONE"


def core_metadata(record: dict) -> dict:
    identifier = record.get("pmid") or record.get("doi")
    if not identifier:
        return {"status": "UNKNOWN", "reason": "no PMID or DOI", "relations": []}
    if record.get("pmid"):
        query = f"EXT_ID:{record['pmid']} AND SRC:MED"
    else:
        query = f'DOI:"{record["doi"]}"'
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "resultType": "core", "pageSize": 1}
    )
    try:
        payload = json.loads(get(url))
        rows = payload.get("resultList", {}).get("result", [])
        if not rows:
            return {"status": "UNKNOWN", "reason": "Europe PMC CORE record not found", "relations": []}
        row = rows[0]
        relations = row.get("commentCorrectionList", {}).get("commentCorrection", [])
        adverse = [r for r in relations if re.search(r"retract|withdraw|expression of concern|erratum|correction", r.get("type", ""), re.I)]
        return {
            "status": "FLAGGED" if adverse else "NO_CORRECTION_OR_RETRACTION_LINK_IN_EUROPE_PMC_CORE",
            "reason": "Europe PMC CORE commentCorrectionList checked",
            "relations": adverse,
            "has_data": row.get("hasData", "UNKNOWN"),
            "data_links_tags": row.get("dataLinksTagsList", "UNKNOWN"),
            "checked_url": url,
        }
    except Exception as exc:
        return {"status": "UNKNOWN", "reason": f"{type(exc).__name__}: {str(exc)[:120]}", "relations": [], "checked_url": url}


def exception_reason(seed: dict) -> str | None:
    if int(seed.get("year") or 0) < 2010:
        return "foundational_pre_2010"
    roles = set(seed.get("aleph_roles", []))
    tags = set(seed.get("topic_tags", []))
    if "negative_evidence" in roles:
        return "negative_evidence"
    if "method" in roles:
        return "validated_method"
    if "benchmark" in roles or "dataset" in tags:
        return "open_benchmark"
    return None


def lab_proxy(authors: str | None) -> str:
    if not authors or authors == "UNKNOWN":
        return "UNKNOWN"
    last = authors.rstrip(". ").split(",")[-1].strip().split()[0]
    return re.sub(r"[^a-z0-9]+", "", last.lower()) or "UNKNOWN"


def dataset_ids(signals: dict) -> list[str]:
    text = signals["source_data"].get("evidence_excerpt", "")
    ids = re.findall(r"\b(?:GSE\d+|PXD\d+|PRJNA\d+|E-MTAB-\d+|10\.\d{4,9}/[-._;()/:A-Z0-9]+)\b", text, flags=re.I)
    return sorted(set(ids)) or ["UNKNOWN"]


def write_journals() -> list[dict]:
    with JOURNAL_QUEUE.open(newline="") as handle:
        queue = list(csv.DictReader(handle))[:50]
    fields = [
        "journal", "jcr_title", "candidate_count", "earliest_year", "latest_year", "metric_name",
        "metric_value", "metric_year", "metric_source", "verified_at", "verification_status", "verification_note", "authority_status",
    ]
    out = []
    for row in queue:
        journal = row["journal"]
        if journal in JIF:
            title, value, abbrev = JIF[journal]
            source = f"https://jcr.clarivate.com/jcr-jp/journal-profile?journal={abbrev}&year=2025&fromPage=%2Fjcr%2Fsearch-results"
            status, note = "verified_licensed_jcr", JCR_DATASET_NOTE
        else:
            title = value = source = "UNKNOWN"
            status = "UNKNOWN_ACCESS_FAILED"
            note = "Chrome extension disconnected after 15 profile verifications; no substitute metric used"
        out.append({
            "journal": journal, "jcr_title": title, "candidate_count": row["candidate_count"],
            "earliest_year": row["earliest_year"], "latest_year": row["latest_year"], "metric_name": "JIF",
            "metric_value": value, "metric_year": "2025" if value != "UNKNOWN" else "UNKNOWN",
            "metric_source": source, "verified_at": VERIFIED_AT, "verification_status": status,
            "verification_note": note, "authority_status": "proposed",
        })
    with (REVIEW / "journal_metrics_verified.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(out)
    return out


def choose_articles(journals: list[dict]) -> list[dict]:
    top50 = {r["journal"] for r in journals}
    records = [json.loads(line) for line in SNAPSHOT.read_text().splitlines() if line.strip()]
    primary = [r for r in records if r.get("journal") in top50 and "research-article" in r.get("publication_types", []) and int(r.get("year") or 0) >= 2010]
    def rank(r: dict) -> tuple:
        metric = JIF.get(r.get("journal"), (None, "UNKNOWN", None))[1]
        known_high = metric != "UNKNOWN" and float(metric) >= 5
        return (known_high, bool(r.get("pmcid")), bool(r.get("is_open_access_provider_flag")), int(r.get("cited_by_count_at_snapshot") or 0))
    primary.sort(key=rank, reverse=True)
    chosen = primary[:120]
    seen = {r["source_family_id"] for r in chosen}
    for line in SEEDS.read_text().splitlines():
        seed = json.loads(line)
        if seed["source_family_id"] in seen:
            continue
        seed["journal"] = "UNKNOWN"
        seed["authors"] = "UNKNOWN"
        seed["publication_types"] = ["UNKNOWN"]
        seed["pmid"] = None
        seed["is_open_access_provider_flag"] = seed.get("access_route") == "open_access"
        seed["cited_by_count_at_snapshot"] = 0
        seed["seed_exception_reason"] = exception_reason(seed)
        chosen.append(seed); seen.add(seed["source_family_id"])
    return chosen


def review_one(record: dict) -> dict:
    signals, access, access_error = fulltext_signals(record.get("pmcid"))
    core = core_metadata(record)
    metric = JIF.get(record.get("journal"), ("UNKNOWN", "UNKNOWN", None))[1]
    is_research = "research-article" in record.get("publication_types", [])
    metric_ok = metric == "UNKNOWN" or float(metric) >= 5
    primary_eligible = int(record.get("year") or 0) >= 2010 and is_research and metric_ok
    exc = record.get("seed_exception_reason")
    required = ["sample_size", "biological_independent_replicates", "calibration", "units", "uncertainty", "exclusions", "source_data"]
    gate_failures = [k for k in required if signals[k]["status"] != "RECOVERED"]
    correction_complete = core["status"] != "UNKNOWN"
    if not correction_complete:
        gate_failures.append("correction_retraction_check")
    # Batch 1 deliberately does not pretend that a section-level keyword hit is
    # an exact observation locator or an independent replication adjudication.
    gate_failures.extend(["exact_observation_locator", "follow_up_replication"])
    if exc:
        lane, decision = "exception", "exception"
    elif not primary_eligible:
        lane, decision = "hold", "hold"
    elif access != "open_full_text_reviewed":
        lane, decision = "primary", "hold"
    else:
        lane, decision = "primary", "hold"
    data_ids = dataset_ids(signals)
    return {
        "authority_status": "proposed",
        "source_family_id": record["source_family_id"], "doi": record.get("doi") or "UNKNOWN",
        "pmid": record.get("pmid") or "UNKNOWN", "pmcid": record.get("pmcid") or "UNKNOWN",
        "title": record.get("title") or "UNKNOWN", "journal": record.get("journal") or "UNKNOWN",
        "year": record.get("year") or "UNKNOWN", "publication_types": record.get("publication_types", ["UNKNOWN"]),
        "tier_candidate": "A", "screening_lane": lane, "quality_decision": decision,
        "quality_gate_failures": gate_failures,
        "exception_reason": exc or "NONE", "jif_metric_name": "JIF", "jif_value": metric,
        "jif_year": "2025" if metric != "UNKNOWN" else "UNKNOWN",
        "article_type_status": "research_article" if is_research else "UNKNOWN",
        "full_text_review_status": access, "full_text_access_error": access_error,
        "biological_system": signals["biological_system"], "cell_state": signals["cell_state"],
        "perturbation": signals["perturbation"], "geometry": signals["geometry"], "substrate": signals["substrate"],
        "sample_size": signals["sample_size"], "biological_independent_replicates": signals["biological_independent_replicates"],
        "calibration": signals["calibration"], "units": signals["units"], "uncertainty": signals["uncertainty"],
        "exclusions": signals["exclusions"], "raw_or_source_data": signals["source_data"],
        "correction_retraction_check": core,
        "follow_up_replication": {"status": "UNKNOWN", "reason": "citation-level replication/contradiction not adjudicated in batch 1"},
        "exact_observation_locator": "UNKNOWN_NOT_EXTRACTED_IN_BATCH1",
        "source_family_leakage_group": record["source_family_id"],
        "lab_leakage_group": lab_proxy(record.get("authors")), "lab_group_basis": "last-author surname proxy; requires affiliation review",
        "dataset_leakage_groups": data_ids,
        "evidence_independence": "supplement_repository_deposit_and_replot_share_this_source_family_and_are_not_independent_evidence",
        "screened_at": VERIFIED_AT,
    }


def write_articles(journals: list[dict]) -> list[dict]:
    candidates = choose_articles(journals)
    reviewed = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {pool.submit(review_one, row): row for row in candidates}
        for index, future in enumerate(as_completed(futures), 1):
            try:
                reviewed.append(future.result())
            except Exception as exc:
                row = futures[future]
                reviewed.append({
                    "authority_status": "proposed", "source_family_id": row["source_family_id"],
                    "title": row.get("title", "UNKNOWN"), "quality_decision": "hold", "screening_lane": "hold",
                    "review_failure": f"{type(exc).__name__}: {str(exc)[:180]}", "screened_at": VERIFIED_AT,
                })
            if index % 25 == 0:
                print(f"reviewed {index}/{len(candidates)}", flush=True)
    reviewed.sort(key=lambda r: r["source_family_id"])
    with (REVIEW / "article_quality_screen.jsonl").open("w") as handle:
        for row in reviewed:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
    fields = ["source_family_id", "doi", "title", "year", "screening_lane", "exception_reason", "quality_decision", "rationale", "authority_status"]
    exceptions = []
    for row in reviewed:
        if row.get("screening_lane") == "exception":
            exceptions.append({
                "source_family_id": row["source_family_id"], "doi": row.get("doi", "UNKNOWN"),
                "title": row.get("title", "UNKNOWN"), "year": row.get("year", "UNKNOWN"),
                "screening_lane": "exception", "exception_reason": row["exception_reason"],
                "quality_decision": row["quality_decision"],
                "rationale": "Explicit policy exception; remains proposed and separate from primary lane",
                "authority_status": "proposed",
            })
    with (REVIEW / "exceptions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(exceptions)
    return reviewed


def write_report(journals: list[dict], articles: list[dict]) -> None:
    verified = [r for r in journals if r["verification_status"] == "verified_licensed_jcr"]
    high = [r for r in verified if float(r["metric_value"]) >= 5]
    access_failed = [r for r in journals if r["verification_status"].startswith("UNKNOWN")]
    decisions = {key: sum(r.get("quality_decision") == key for r in articles) for key in ("tier_a", "hold", "reject", "exception")}
    fulltext_failed = sum(r.get("full_text_review_status") != "open_full_text_reviewed" for r in articles)
    primary_candidates = sum(r.get("article_type_status") == "research_article" for r in articles)
    report = f"""# External training corpus quality review — first batch

**Status:** AGENT-PROPOSED. No record in this review is an Aleph claim or training authorization.

## Counts

- Candidate-rich journals attempted: **{len(journals)}**
- JIF profiles actually verified in licensed Clarivate JCR: **{len(verified)}**
- Verified JIF >= 5: **{len(high)}**
- JIF access failures left `UNKNOWN`: **{len(access_failed)}**
- Tier A candidates screened at article level: **{primary_candidates}**
- Article records reviewed: **{len(articles)}**
- Tier A passed: **{decisions['tier_a']}**
- Holds: **{decisions['hold']}**
- Rejects: **{decisions['reject']}**
- Explicit exceptions: **{decisions['exception']}**
- Article full-text access failures or metadata-only reviews: **{fulltext_failed}**

## Method and limits

The JIF batch records the 2025 metric year, value, licensed JCR profile URL, and verification date. It does not substitute CiteScore, SJR, or inferred values. Chrome institutional access disconnected after 15 successful profile checks; the remaining 35 values are `UNKNOWN` rather than guessed.

Article screening is a conservative first pass over at least 100 research-article candidates from the candidate-rich top-50 journals, plus explicit seed exceptions. Recoverability signals come from Europe PMC full-text XML; correction/retraction links come from Europe PMC CORE `commentCorrectionList`. A keyword hit records a short source excerpt and section locator but does not itself prove methodological adequacy. Missing sample size, biological independence, calibration, units, uncertainty, exclusion, or source-data evidence causes a hold. Citation-level independent replication and exact figure/table observation locators remain `UNKNOWN` in batch 1, so JIF alone never produces Tier A.

## Leakage controls

Every paper, supplement, repository deposit, and replot remains in one `source_family_leakage_group` and is not counted as independent evidence. `lab_leakage_group` is only a last-author-surname proxy and therefore cannot authorize a split until affiliations are manually resolved. Extracted accession identifiers form dataset leakage groups; absent identifiers remain `UNKNOWN`.

## Next safe action

Restore Chrome extension connectivity and verify the remaining 35 top-journal JCR profiles. Then manually adjudicate exact figure/table locators, field applicability, affiliation-resolved laboratory groups, dataset identities, and citation-level replication/contradiction before any hold can become Tier A. All records remain `proposed`.
"""
    (REVIEW / "QUALITY_REVIEW_REPORT.md").write_text(report)


def main() -> None:
    REVIEW.mkdir(parents=True, exist_ok=True)
    journals = write_journals()
    articles = write_articles(journals)
    write_report(journals, articles)
    print(f"journals={len(journals)} articles={len(articles)}")


if __name__ == "__main__":
    main()
