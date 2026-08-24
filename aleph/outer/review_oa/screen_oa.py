#!/usr/bin/env python3
"""Conservative deterministic recoverability screening over local OA JATS chunks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

SCHEMA = "aleph.external_training.oa_recoverability_screen.v1"
COMPACT_SCHEMA = "aleph.external_training.oa_recoverability_compact.v1"
VERSION = "1"
MAX_REFS = 8


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        handle.write(data)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def patterns(*pairs: tuple[str, str]) -> tuple[tuple[str, re.Pattern[str]], ...]:
    return tuple((name, re.compile(expr, re.I)) for name, expr in pairs)


SIGNALS = {
    "sample_size": patterns(
        ("n_equals", r"\b[nN]\s*[=:]\s*\d{1,6}\b"),
        ("sample_size_phrase", r"\bsample size\b"),
        ("counted_cells", r"\b\d{2,6}\s+(?:cells|nuclei|animals|mice|patients|donors|subjects|samples)\b"),
    ),
    "biological_replicates": patterns(
        ("biological_replicate", r"\bbiological replicat(?:e|es|ion)\b"),
        ("independent_experiments", r"\b(?:at least\s+)?\d+\s+independent experiments\b"),
        ("independent_donors", r"\bindependent (?:donors|animals|cultures|samples)\b"),
    ),
    "calibration": patterns(
        ("calibrated", r"\bcalibrat(?:e|ed|es|ing|ion)\b"),
        ("calibration_curve", r"\b(?:standard|calibration) curve\b"),
        ("instrument_standard", r"\b(?:calibration bead|force calibration|pixel calibration|reference standard)\b"),
    ),
    "units": patterns(
        ("si_prefixed_unit", r"(?<![A-Za-z])(?:\d+(?:\.\d+)?\s*)?(?:pN|nN|µN|uN|Pa|kPa|MPa|µm|um|nm|mm|mPa[· ]?s|Pa[· ]?s)(?![A-Za-z])"),
        ("rate_unit", r"(?<![A-Za-z])(?:µm|um|nm|mm)\s*(?:/|per)\s*(?:s|min|h)(?![A-Za-z])"),
    ),
    "uncertainty": patterns(
        ("error_statistic", r"\b(?:standard deviation|standard error|SEM|SD|confidence interval|credible interval)\b"),
        ("error_bars", r"\berror bars?\b"),
        ("plus_minus", r"(?:±|\+/-)"),
    ),
    "exclusions": patterns(
        ("exclusion", r"\b(?:exclusion|excluded|excluding|inclusion criteria|outlier removal|removed as outliers?)\b"),
    ),
    "source_data": patterns(
        ("data_availability", r"\b(?:data availability|source data|data are available|data is available|availability of data)\b"),
        ("repository", r"\b(?:Gene Expression Omnibus|Sequence Read Archive|BioProject|ArrayExpress|Zenodo|Figshare|Dryad)\b"),
    ),
}

DATASET_PATTERNS = (
    re.compile(r"\bGSE\d{3,9}\b", re.I),
    re.compile(r"\b(?:SRP|SRA|ERP|DRP)\d{3,9}\b", re.I),
    re.compile(r"\bPRJ(?:NA|EB|DB)\d{3,12}\b", re.I),
    re.compile(r"\bE-MTAB-\d{2,8}\b", re.I),
)
RARE_MODALITIES = {"magnetic_tweezers", "micropipette_aspiration", "optical_tweezers", "piv"}


def iter_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def load_metadata(directory: Path) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*.jsonl")):
        for record in iter_jsonl(path):
            family = record["source_family_id"]
            if family not in result:
                result[family] = record
            else:
                for key, value in record.items():
                    if result[family].get(key) in (None, "", [], {}) and value not in (None, "", [], {}):
                        result[family][key] = value
    return result


def load_prior_review(path: Path | None) -> dict[str, dict[str, Any]]:
    if path is None or not path.is_file():
        return {}
    return {record["source_family_id"]: record for record in iter_jsonl(path)}


def article_type(metadata: dict[str, Any]) -> dict[str, Any]:
    types = [str(value) for value in metadata.get("publication_types", [])]
    lowered = {value.lower() for value in types}
    title = str(metadata.get("title", "")).lower()
    if "research-article" in lowered:
        status = "research_article_metadata"
    elif any("review" in value for value in lowered) or title.startswith(("review", "editorial", "commentary")):
        status = "non_research_metadata"
    else:
        status = "unknown"
    return {"status": status, "publication_types": sorted(types), "source": "corpus_snapshot_metadata"}


def evidence_ref(chunk: dict[str, Any], pattern_id: str) -> dict[str, str]:
    return {
        "locator": chunk["section_locator"],
        "section_class": chunk["section_class"],
        "chunk_sha256": chunk["chunk_sha256"],
        "text_sha256": chunk["text_sha256"],
        "pattern_id": pattern_id,
    }


def correction_status(family: str, metadata: dict[str, Any], prior: dict[str, dict[str, Any]]) -> dict[str, Any]:
    if family in prior:
        raw = prior[family].get("correction_retraction_check", {})
        return {
            "status": str(raw.get("status", "UNKNOWN")).lower(),
            "relation_count": len(raw.get("relations", [])),
            "source": "prior_corpus_review",
            "uncertainty": "relation semantics require manual adjudication" if raw.get("relations") else "check is dated and not a replication review",
        }
    title = str(metadata.get("title", ""))
    if re.search(r"\b(?:retracted|retraction|correction|expression of concern)\b", title, re.I):
        return {"status": "flagged_title_metadata", "relation_count": 0, "source": "corpus_snapshot_title", "uncertainty": "title-only flag requires manual adjudication"}
    return {"status": "metadata_unavailable", "relation_count": 0, "source": "none", "uncertainty": "no corpus correction/retraction relation was available for this source"}


def screen_article(manifest: dict[str, Any], chunks_path: Path, metadata: dict[str, Any], prior: dict[str, dict[str, Any]]) -> dict[str, Any]:
    signal_refs: dict[str, list[dict[str, str]]] = {name: [] for name in SIGNALS}
    tag_refs: dict[str, dict[str, list[str]]] = {name: defaultdict(list) for name in ("cell_type_state", "measurement_modality", "mechanics_observable", "perturbation")}
    observation: dict[tuple[str, str], dict[str, Any]] = {}
    dataset_ids: set[str] = set()

    for chunk in iter_jsonl(chunks_path):
        text = chunk["text"]
        for name, rules in SIGNALS.items():
            if len(signal_refs[name]) >= MAX_REFS:
                continue
            for pattern_id, rule in rules:
                if rule.search(text):
                    signal_refs[name].append(evidence_ref(chunk, pattern_id))
                    break
        for dimension, tags in chunk["tags"].items():
            for tag in tags:
                locators = tag_refs[dimension][tag]
                if chunk["section_locator"] not in locators and len(locators) < MAX_REFS:
                    locators.append(chunk["section_locator"])
        for rule in DATASET_PATTERNS:
            dataset_ids.update(match.upper() for match in rule.findall(text))
        if chunk["section_class"] in {"figure_caption", "table_caption"}:
            relevant = sorted(set(chunk["tags"]["measurement_modality"] + chunk["tags"]["mechanics_observable"]))
            quantitative = any(rule.search(text) for _, rule in SIGNALS["units"] + SIGNALS["uncertainty"] + SIGNALS["sample_size"])
            if relevant or quantitative:
                key = (chunk["section_locator"], chunk["text_sha256"])
                observation[key] = {
                    "locator": chunk["section_locator"],
                    "section_class": chunk["section_class"],
                    "text_sha256": chunk["text_sha256"],
                    "tags": relevant,
                    "quantitative_signal": quantitative,
                }

    signals = {
        name: {"status": "recoverable_signal" if refs else "not_recovered", "count": len(refs), "evidence": refs}
        for name, refs in signal_refs.items()
    }
    tags = {dimension: {tag: locs for tag, locs in sorted(values.items())} for dimension, values in tag_refs.items()}
    observation_locators = [observation[key] for key in sorted(observation)]
    atype = article_type(metadata)
    correction = correction_status(manifest["source_family_id"], metadata, prior)
    source_group = manifest["source_family_id"]
    dataset_groups = [{"accession": value, "group_hash": sha(("dataset:" + value).encode())[:16]} for value in sorted(dataset_ids)]
    authors = str(metadata.get("authors", ""))
    last_author_proxy = authors.rstrip(".").split(",")[-1].strip() if authors else "UNKNOWN"
    lab_proxy_hash = sha(("last-author-proxy:" + last_author_proxy.lower()).encode())[:16] if authors else "UNKNOWN"

    missing = [name for name, value in signals.items() if value["status"] == "not_recovered"]
    if not tags["cell_type_state"]:
        missing.append("cell_type_state")
    if not tags["measurement_modality"]:
        missing.append("measurement_modality")
    if not observation_locators:
        missing.append("exact_observation_locator")
    if correction["status"] in {"metadata_unavailable", "unknown"}:
        missing.append("correction_retraction_check")
    missing.extend(["independent_replication_check", "affiliation_resolved_lab_group"])
    missing = sorted(set(missing))

    score = 0
    score += 3 if atype["status"] == "research_article_metadata" else 0
    weights = {"sample_size": 2, "biological_replicates": 3, "calibration": 2, "units": 1, "uncertainty": 2, "exclusions": 1, "source_data": 2}
    score += sum(weight for name, weight in weights.items() if signals[name]["status"] == "recoverable_signal")
    score += 3 if observation_locators else 0
    score += 2 if tags["mechanics_observable"] else 0
    score += 1 if tags["measurement_modality"] else 0
    score += 1 if tags["cell_type_state"] else 0
    warning = correction["status"] in {"flagged", "flagged_title_metadata"}
    if warning:
        score -= 10

    core = all(signals[name]["status"] == "recoverable_signal" for name in ("sample_size", "biological_replicates", "uncertainty"))
    basis: list[str] = []
    exception_reason = "NONE"
    modality_set = set(tags["measurement_modality"])
    title = str(metadata.get("title", ""))
    if atype["status"] == "research_article_metadata" and core and observation_locators and not warning and score >= 12:
        decision = "manual_priority"
        basis.extend(["research-article metadata", "core recoverability signals", "figure/table observation locator"])
    elif atype["status"] != "research_article_metadata" and observation_locators and modality_set & RARE_MODALITIES and not warning:
        decision = "exception_candidate"
        exception_reason = "rare_modality"
        basis.extend(["non-primary or unknown article type", "rare modality", "figure/table observation locator"])
    elif atype["status"] != "research_article_metadata" and re.search(r"\b(?:method|protocol|tool|pipeline)\b", title, re.I) and observation_locators and not warning:
        decision = "exception_candidate"
        exception_reason = "validated_method"
        basis.extend(["method-title signal", "figure/table observation locator", "manual validation required"])
    else:
        decision = "hold"
        basis.append("article-level quality gates remain incomplete")
        if warning:
            basis.append("correction/retraction metadata warning")

    return {
        "schema": SCHEMA,
        "authority_status": "proposed",
        "source_family_id": manifest["source_family_id"],
        "pmcid": manifest["pmcid"],
        "payload_sha256": manifest["payload_sha256"],
        "input_chunks_sha256": manifest["chunks_file_sha256"],
        "article_type": atype,
        "signals": signals,
        "tags": tags,
        "observation_locators": observation_locators,
        "correction_retraction": correction,
        "leakage": {
            "source_family_group": source_group,
            "dataset_groups": dataset_groups,
            "lab_group_proxy_hash": lab_proxy_hash,
            "lab_group_basis": "last-author string proxy; affiliations not present in derived corpus and manual resolution is required",
        },
        "missing_review_fields": missing,
        "decision": decision,
        "exception_reason": exception_reason,
        "decision_basis": sorted(set(basis)),
        "priority_score": score,
    }


def compact(record: dict[str, Any], record_sha: str) -> dict[str, Any]:
    return {
        "schema": COMPACT_SCHEMA,
        "authority_status": "proposed",
        "source_family_id": record["source_family_id"],
        "pmcid": record["pmcid"],
        "record_sha256": record_sha,
        "decision": record["decision"],
        "exception_reason": record["exception_reason"],
        "priority_score": record["priority_score"],
        "research_article_status": record["article_type"]["status"],
        "signal_counts": {name: value["count"] for name, value in sorted(record["signals"].items())},
        "missing_review_fields": record["missing_review_fields"],
        "tag_counts": {name: len(value) for name, value in sorted(record["tags"].items())},
        "observation_locator_count": len(record["observation_locators"]),
        "dataset_group_hashes": [value["group_hash"] for value in record["leakage"]["dataset_groups"]],
    }


def validate_extraction(derived: Path) -> list[tuple[dict[str, Any], Path, Path]]:
    summary = json.loads((derived / "final_summary.json").read_text("utf-8"))
    if summary["failures"] or summary["xml_parsed"] != 2701:
        raise ValueError("final extraction is incomplete")
    entries = []
    manifests = []
    for path in sorted((derived / "articles").glob("*.manifest.json")):
        manifest = json.loads(path.read_text("utf-8"))
        chunks_path = path.with_name(path.name.replace(".manifest.json", ".chunks.jsonl"))
        if not chunks_path.is_file() or sha(chunks_path.read_bytes()) != manifest["chunks_file_sha256"]:
            raise ValueError(f"chunk digest mismatch: {chunks_path}")
        entries.append((manifest, chunks_path, path))
        manifests.append(manifest)
    manifest_hash = sha(b"".join(canonical_bytes(value) + b"\n" for value in sorted(manifests, key=lambda item: item["source_family_id"])))
    if len(entries) != summary["xml_parsed"] or manifest_hash != summary["article_manifest_set_sha256"]:
        raise ValueError("manifest set does not reproduce final summary")
    return entries


def process(derived: Path, metadata_dir: Path, prior_path: Path | None, local_root: Path, output_root: Path) -> dict[str, Any]:
    entries = validate_extraction(derived)
    metadata = load_metadata(metadata_dir)
    prior = load_prior_review(prior_path)
    records_dir = local_root / "records"
    records_dir.mkdir(parents=True, exist_ok=True)
    compact_records = []
    records = []
    dataset_members: dict[str, list[str]] = defaultdict(list)
    for manifest, chunks_path, _ in entries:
        family = manifest["source_family_id"]
        record = screen_article(manifest, chunks_path, metadata.get(family, {}), prior)
        data = canonical_bytes(record) + b"\n"
        record_sha = sha(data)
        local_name = sha(family.encode())[:24] + ".json"
        atomic_write(records_dir / local_name, data)
        compact_record = compact(record, record_sha)
        compact_records.append(compact_record)
        records.append(record)
        for group in record["leakage"]["dataset_groups"]:
            dataset_members[group["group_hash"]].append(family)

    compact_records.sort(key=lambda item: item["source_family_id"])
    compact_data = b"".join(canonical_bytes(item) + b"\n" for item in compact_records)
    atomic_write(output_root / "article_screen_index.jsonl", compact_data)

    queue_records = sorted(
        (record for record in records if record["decision"] != "hold"),
        key=lambda item: (-item["priority_score"], item["source_family_id"]),
    )
    queue_path = output_root / "manual_review_queue.csv"
    queue_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=queue_path.parent, delete=False) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["rank", "source_family_id", "pmcid", "decision", "exception_reason", "priority_score", "missing_review_fields", "observation_locator_count", "cell_state_tags", "modality_tags", "mechanics_tags"])
        for rank, record in enumerate(queue_records, 1):
            writer.writerow([
                rank, record["source_family_id"], record["pmcid"], record["decision"], record["exception_reason"], record["priority_score"],
                ";".join(record["missing_review_fields"]), len(record["observation_locators"]),
                ";".join(record["tags"]["cell_type_state"]), ";".join(record["tags"]["measurement_modality"]), ";".join(record["tags"]["mechanics_observable"]),
            ])
        temporary = Path(handle.name)
    os.replace(temporary, queue_path)

    leakage = []
    for group_hash, families in sorted(dataset_members.items()):
        if len(families) > 1:
            leakage.append({"dataset_group_hash": group_hash, "source_family_count": len(families), "source_families": sorted(families), "uncertainty": "accession co-membership; dataset reuse relation requires manual confirmation"})
    leakage_data = canonical_bytes({"schema": "aleph.external_training.oa_dataset_leakage_groups.v1", "authority_status": "proposed", "groups": leakage}) + b"\n"
    atomic_write(output_root / "dataset_leakage_groups.json", leakage_data)

    decisions = Counter(record["decision"] for record in records)
    article_types = Counter(record["article_type"]["status"] for record in records)
    missing = Counter(field for record in records for field in record["missing_review_fields"])
    signal_coverage = {name: sum(record["signals"][name]["status"] == "recoverable_signal" for record in records) for name in SIGNALS}
    correction = Counter(record["correction_retraction"]["status"] for record in records)
    summary = {
        "schema": "aleph.external_training.oa_recoverability_summary.v1",
        "authority_status": "proposed",
        "screen_version": VERSION,
        "articles": len(records),
        "decisions": dict(sorted(decisions.items())),
        "article_types": dict(sorted(article_types.items())),
        "signal_recoverability_articles": dict(sorted(signal_coverage.items())),
        "missing_review_fields": dict(sorted(missing.items())),
        "correction_retraction_status": dict(sorted(correction.items())),
        "articles_with_observation_locator": sum(bool(record["observation_locators"]) for record in records),
        "articles_with_dataset_accession": sum(bool(record["leakage"]["dataset_groups"]) for record in records),
        "duplicate_dataset_groups": len(leakage),
        "local_record_set_sha256": sha(b"".join((item["source_family_id"] + ":" + item["record_sha256"] + "\n").encode() for item in compact_records)),
        "compact_index_sha256": sha(compact_data),
        "manual_queue_sha256": sha(queue_path.read_bytes()),
        "leakage_groups_sha256": sha(leakage_data),
        "input_manifest_set_sha256": json.loads((derived / "final_summary.json").read_text("utf-8"))["article_manifest_set_sha256"],
        "limitations": [
            "regex and vocabulary hits indicate recoverability only, not adequacy or truth",
            "correction/retraction metadata is unavailable for most sources in the local corpus",
            "independent replication and affiliation-resolved laboratory groups require manual review",
            "no decision in this lane is Tier A or training authorization"
        ]
    }
    summary_data = json.dumps(summary, indent=2, sort_keys=True) .encode("utf-8") + b"\n"
    atomic_write(output_root / "summary.json", summary_data)
    sums = {path.name: sha(path.read_bytes()) for path in sorted(output_root.iterdir()) if path.is_file() and path.name != "SHA256SUMS.json"}
    atomic_write(output_root / "SHA256SUMS.json", json.dumps(sums, indent=2, sort_keys=True).encode("utf-8") + b"\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derived-root", type=Path, required=True)
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--prior-review", type=Path)
    parser.add_argument("--local-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    summary = process(args.derived_root, args.metadata_dir, args.prior_review, args.local_root, args.output_root)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
