#!/usr/bin/env python
"""Extract and optionally resolve exact public dataset identifiers from the OA JATS corpus."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[4]
CORPUS = ROOT / "corpus" / "external_training"
HERE = Path(__file__).resolve().parent
LOCAL_FIRST = ROOT / "data" / "external_training" / "review_oa" / "records"
LOCAL_SECOND = ROOT / "data" / "external_training" / "review_oa_second_pass" / "records"
GOLD = CORPUS / "experiment_factory" / "goldset" / "candidates_300.jsonl"
FIRST_RECEIPTS = CORPUS / "acquisition" / "europe_pmc_oa_2026-08-05.jsonl"
SECOND_RECEIPTS = CORPUS / "acquisition" / "second_pass" / "receipts_2026-08-05.jsonl"
OFFICIAL_RECEIPTS = HERE / "official_metadata_receipts.jsonl"
PILOT_RECEIPTS = HERE / "pilot_download_receipts.jsonl"
SCAN_TAGS = {"p", "uri", "object-id", "custom-meta-value", "data-title", "related-object"}
CONTAINER_TAGS = {"sec", "fig", "table-wrap", "supplementary-material", "notes", "app"}
USER_AGENT = "Project-Aleph-public-dataset-registry/1.0 (metadata-only research client)"


PATTERNS: list[tuple[str, str, str, re.Pattern[str]]] = [
    ("GEO", "geo_series", "accession", re.compile(r"\bGSE\d{3,}\b", re.I)),
    ("SRA", "sra_study", "accession", re.compile(r"\b(?:SRP|ERP|DRP)\d{3,}\b", re.I)),
    ("SRA", "sra_run_or_sample", "accession", re.compile(r"\b(?:SRR|ERR|DRR|SRS|ERS|DRS|SRX|ERX|DRX)\d{3,}\b", re.I)),
    ("BioProject", "bioproject", "accession", re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.I)),
    ("BioSample", "biosample", "accession", re.compile(r"\bSAM(?:N|EA|D)\d+\b", re.I)),
    ("ArrayExpress", "arrayexpress_study", "accession", re.compile(r"\bE-(?:MTAB|GEOD|MEXP)-\d+\b", re.I)),
    ("PRIDE", "proteomexchange_project", "accession", re.compile(r"\bPXD\d{4,}\b", re.I)),
    ("MassIVE", "massive_dataset", "accession", re.compile(r"\bMSV\d{6,}\b", re.I)),
    ("MetaboLights", "metabolights_study", "accession", re.compile(r"\bMTBLS\d+\b", re.I)),
    ("BioImageArchive", "bioimage_archive_study", "accession", re.compile(r"\bS-BIAD\d+\b", re.I)),
    ("EMPIAR", "empiar_entry", "accession", re.compile(r"\bEMPIAR-\d+\b", re.I)),
    ("IDR", "image_data_resource", "accession", re.compile(r"\bidr\d{4,}\b", re.I)),
    ("FlowRepository", "flow_repository", "accession", re.compile(r"\bFR-FCM-[A-Z0-9-]+\b", re.I)),
    ("dbGaP", "dbgap_study", "accession", re.compile(r"\bphs\d+(?:\.v\d+\.p\d+)?\b", re.I)),
    ("Zenodo", "doi", "doi", re.compile(r"10\.5281/zenodo\.\d+", re.I)),
    ("Zenodo", "record", "url", re.compile(r"zenodo\.org/(?:record|records)/\d+", re.I)),
    ("Figshare", "doi", "doi", re.compile(r"10\.6084/m9\.figshare\.\d+(?:\.v\d+)?", re.I)),
    ("Dryad", "doi", "doi", re.compile(r"10\.5061/dryad\.[A-Za-z0-9._-]+", re.I)),
    ("OSF", "osf_object", "url", re.compile(r"osf\.io/[a-z0-9]{4,6}", re.I)),
    ("MendeleyData", "doi", "doi", re.compile(r"10\.17632/[a-z0-9.]+", re.I)),
    ("Dataverse", "doi", "doi", re.compile(r"10\.7910/DVN/[A-Z0-9]+", re.I)),
    ("GigaDB", "doi", "doi", re.compile(r"10\.5524/\d+", re.I)),
]


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def normalize_text(value: str) -> str:
    return " ".join(value.split())


def canonical_identifier(provider: str, kind: str, value: str) -> str:
    value = value.strip().rstrip(".,;:)]}")
    if kind == "accession":
        return value.upper() if provider != "dbGaP" else value.lower()
    lower = value.lower()
    if provider == "Zenodo":
        number = re.search(r"(?:zenodo\.|records?/)(\d+)", lower).group(1)
        return f"10.5281/zenodo.{number}"
    if provider == "OSF":
        return "osf:" + lower.rsplit("/", 1)[-1]
    return lower


def provider_url(provider: str, accession: str) -> str | None:
    encoded = urllib.parse.quote(accession, safe="-.:/")
    if provider == "GEO": return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={encoded}"
    if provider in {"SRA", "BioProject", "BioSample", "dbGaP"}: return f"https://www.ncbi.nlm.nih.gov/search/all/?term={encoded}"
    if provider == "ArrayExpress": return f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{encoded}"
    if provider == "PRIDE": return f"https://www.ebi.ac.uk/pride/archive/projects/{encoded}"
    if provider == "MassIVE": return f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task={encoded}"
    if provider == "MetaboLights": return f"https://www.ebi.ac.uk/metabolights/{encoded}"
    if provider == "BioImageArchive": return f"https://www.ebi.ac.uk/biostudies/BioImages/studies/{encoded}"
    if provider == "EMPIAR": return f"https://www.ebi.ac.uk/empiar/{encoded.replace('EMPIAR-', '')}/"
    if provider == "IDR": return f"https://idr.openmicroscopy.org/webclient/?show=project-{encoded[3:]}"
    if provider == "FlowRepository": return f"https://flowrepository.org/id/{encoded}"
    if provider == "Zenodo": return f"https://doi.org/{encoded}"
    if provider == "Figshare": return f"https://doi.org/{encoded}"
    if provider == "Dryad": return f"https://doi.org/{encoded}"
    if provider == "OSF": return f"https://osf.io/{accession.split(':', 1)[-1]}/"
    if provider in {"MendeleyData", "Dataverse", "GigaDB"}: return f"https://doi.org/{encoded}"
    return None


def matches(text: str) -> list[tuple[str, str, str]]:
    found: set[tuple[str, str, str]] = set()
    for provider, identifier_type, kind, pattern in PATTERNS:
        for match in pattern.finditer(text):
            found.add((provider, identifier_type, canonical_identifier(provider, kind, match.group(0))))
    return sorted(found)


def walk_blocks(root: ET.Element) -> Iterable[tuple[str, str]]:
    counters: Counter[tuple[str, str]] = Counter()

    def walk(element: ET.Element, containers: list[str]) -> Iterable[tuple[str, str]]:
        tag = local_name(element.tag)
        next_containers = containers
        if tag in CONTAINER_TAGS:
            identifier = element.attrib.get("id") or f"{tag}:{len(containers)}"
            next_containers = containers + [f"{tag}:{identifier}"]
        if tag in SCAN_TAGS:
            text = normalize_text("".join(element.itertext()))
            if text:
                base = "/".join(next_containers) or "article"
                counters[(base, tag)] += 1
                own = element.attrib.get("id")
                locator = f"{base}/{tag}:{own or counters[(base, tag)]}"
                yield locator, text
        for child in element:
            yield from walk(child, next_containers)

    yield from walk(root, [])


def load_quality() -> dict[str, dict[str, Any]]:
    records = {}
    for directory in (LOCAL_FIRST, LOCAL_SECOND):
        for path in sorted(directory.glob("*.json")):
            record = json.loads(path.read_text(encoding="utf-8"))
            records[record["source_family_id"].lower()] = record
    return records


def modality_hints(source_id: str, quality: dict[str, dict[str, Any]], provider: str) -> list[str]:
    tags = set()
    if source_id in quality:
        tags.update(quality[source_id].get("tags", {}).get("measurement_modality", {}).keys())
    provider_hints = {
        "GEO": {"omics", "pcr"}, "SRA": {"sequencing", "omics"},
        "BioProject": {"sequencing", "omics"}, "BioSample": {"omics"},
        "ArrayExpress": {"omics"}, "PRIDE": {"proteomics", "western_blot"},
        "MassIVE": {"proteomics"}, "MetaboLights": {"metabolomics"},
        "BioImageArchive": {"raw_imaging", "immunofluorescence"},
        "EMPIAR": {"raw_imaging"}, "IDR": {"raw_imaging", "immunofluorescence"},
        "FlowRepository": {"flow_cytometry"},
    }
    tags.update(provider_hints.get(provider, set()))
    return sorted(tags)


def load_existing_official_receipts() -> dict[tuple[str, str], dict[str, Any]]:
    return {(row["provider"], row["accession"]): row for row in read_jsonl(OFFICIAL_RECEIPTS)}


def official_endpoint(provider: str, accession: str) -> str | None:
    q = urllib.parse.quote(accession, safe="-.:/")
    if provider == "GEO": return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={q}&targ=self&form=text&view=quick"
    if provider in {"SRA", "BioProject", "BioSample", "dbGaP"}:
        db = {"SRA": "sra", "BioProject": "bioproject", "BioSample": "biosample", "dbGaP": "gap"}[provider]
        return f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={q}[accn]&retmode=json"
    if provider == "ArrayExpress": return f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{q}"
    if provider == "PRIDE": return f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{q}"
    if provider == "MetaboLights": return f"https://www.ebi.ac.uk/metabolights/ws/studies/{q}"
    if provider in {"BioImageArchive", "EMPIAR"}: return f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{q}"
    if provider == "Zenodo": return f"https://zenodo.org/api/records/{accession.rsplit('.', 1)[-1]}"
    if provider == "Figshare":
        match = re.search(r"figshare\.(\d+)", accession)
        return f"https://api.figshare.com/v2/articles/{match.group(1)}" if match else None
    if provider == "Dryad": return f"https://datadryad.org/api/v2/datasets/doi:{q}"
    if provider == "OSF": return f"https://api.osf.io/v2/nodes/{accession.split(':', 1)[-1]}/"
    return None


def parse_official(provider: str, raw: bytes) -> dict[str, Any]:
    title = None
    license_id = None
    file_count = None
    total_size = None
    files: list[dict[str, Any]] = []
    if provider == "GEO":
        text = raw.decode("utf-8", "replace")
        hit = re.search(r"^!Series_title\s*=\s*(.+)$", text, re.M)
        title = hit.group(1).strip() if hit else None
        urls = re.findall(r"^!Series_supplementary_file(?:_\d+)?\s*=\s*(\S+)", text, re.M)
        file_count = len(urls) if urls else None
    else:
        data = json.loads(raw)
        if provider in {"SRA", "BioProject", "BioSample", "dbGaP"}:
            found = int(data.get("esearchresult", {}).get("count", 0))
            if found == 0:
                raise LookupError("official index returned zero records")
        elif provider == "Zenodo":
            title = data.get("metadata", {}).get("title")
            license_id = (data.get("metadata", {}).get("license") or {}).get("id")
            raw_files = data.get("files") or []
            file_count = len(raw_files)
            total_size = sum(int(item.get("size") or 0) for item in raw_files)
            files = [{"name": item.get("key"), "size": item.get("size"), "url": (item.get("links") or {}).get("self")} for item in raw_files]
        elif provider == "Figshare":
            title = data.get("title")
            license_id = (data.get("license") or {}).get("name")
            raw_files = data.get("files") or []
            file_count = len(raw_files)
            total_size = sum(int(item.get("size") or 0) for item in raw_files)
            files = [{"name": item.get("name"), "size": item.get("size"), "url": item.get("download_url")} for item in raw_files]
        elif provider == "Dryad":
            title = data.get("title")
            license_id = data.get("license")
        elif provider == "OSF":
            attributes = (data.get("data") or {}).get("attributes") or {}
            title = attributes.get("title")
        elif provider == "PRIDE":
            title = data.get("title")
        elif provider in {"ArrayExpress", "BioImageArchive", "EMPIAR", "MetaboLights"}:
            title = data.get("title") or data.get("section", {}).get("title")
    return {"title": title, "license_identifier": license_id, "file_count": file_count, "total_size_bytes": total_size, "files": files}


def resolve_official(provider: str, accession: str) -> dict[str, Any]:
    endpoint = official_endpoint(provider, accession)
    if endpoint is None:
        return {"provider": provider, "accession": accession, "status": "not_attempted", "official_url": None, "http_status": None, "retrieved_at": None, "response_sha256": None, "error": "no official metadata adapter", "metadata": {}}
    request = urllib.request.Request(endpoint, headers={"User-Agent": USER_AGENT, "Accept": "application/json,text/plain;q=0.9,*/*;q=0.1"})
    retrieved = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            raw = response.read(8 * 1024 * 1024 + 1)
            if len(raw) > 8 * 1024 * 1024:
                raise ValueError("metadata response exceeds 8 MiB cap")
            metadata = parse_official(provider, raw)
            return {"provider": provider, "accession": accession, "status": "resolved", "official_url": endpoint, "http_status": response.status, "retrieved_at": retrieved, "response_sha256": sha(raw), "error": None, "metadata": metadata}
    except urllib.error.HTTPError as exc:
        return {"provider": provider, "accession": accession, "status": "not_found" if exc.code == 404 else "failed", "official_url": endpoint, "http_status": exc.code, "retrieved_at": retrieved, "response_sha256": None, "error": f"HTTP {exc.code}", "metadata": {}}
    except (urllib.error.URLError, TimeoutError, ValueError, json.JSONDecodeError, LookupError) as exc:
        return {"provider": provider, "accession": accession, "status": "not_found" if isinstance(exc, LookupError) else "failed", "official_url": endpoint, "http_status": None, "retrieved_at": retrieved, "response_sha256": None, "error": type(exc).__name__, "metadata": {}}


def build(resolve_limit: int = 0) -> dict[str, Any]:
    gold = {row["source_family_id"] for row in read_jsonl(GOLD)}
    quality = load_quality()
    receipts = []
    for acquisition_pass, path in (("first", FIRST_RECEIPTS), ("second", SECOND_RECEIPTS)):
        for row in read_jsonl(path):
            if row.get("status") == "success":
                receipts.append((acquisition_pass, row))

    bindings: dict[tuple[str, str, str], dict[str, Any]] = {}
    for acquisition_pass, receipt in receipts:
        payload_path = ROOT / receipt["object_path"]
        root = ET.parse(payload_path).getroot()
        source_id = receipt["source_family_id"].lower()
        for locator, text in walk_blocks(root):
            found = matches(text)
            if not found:
                continue
            text_sha = sha(text.encode("utf-8"))
            for provider, identifier_type, accession in found:
                key = (source_id, provider, accession)
                if key not in bindings:
                    hints = modality_hints(source_id, quality, provider)
                    bindings[key] = {
                        "schema": "aleph.external_training.source_accession.v1",
                        "authority_status": "proposed",
                        "source_family_id": source_id,
                        "pmcid": receipt["pmcid"],
                        "payload_sha256": receipt["sha256"],
                        "gold_candidate": source_id in gold,
                        "provider": provider,
                        "accession": accession,
                        "identifier_type": identifier_type,
                        "canonical_url": provider_url(provider, accession),
                        "leakage_component_id": "ds:" + sha(f"{provider}\x1f{accession}".encode("utf-8"))[:16],
                        "modality_hints": hints,
                        "evidence": [],
                    }
                evidence = {"document_sha256": receipt["sha256"], "locator": locator, "text_sha256": text_sha}
                if evidence not in bindings[key]["evidence"] and len(bindings[key]["evidence"]) < 8:
                    bindings[key]["evidence"].append(evidence)

    binding_rows = sorted(bindings.values(), key=lambda row: (row["provider"], row["accession"], row["source_family_id"]))
    binding_bytes = b"".join(canonical_bytes(row) for row in binding_rows)
    (HERE / "source_accession_registry.jsonl").write_bytes(binding_bytes)

    by_dataset: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in binding_rows:
        by_dataset[(row["provider"], row["accession"])].append(row)
    official = load_existing_official_receipts()
    pilot_rows = read_jsonl(PILOT_RECEIPTS)
    pilot_by_dataset: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in pilot_rows:
        pilot_by_dataset[(row["provider"], row["accession"])].append(row)
    if resolve_limit:
        ranked = sorted(
            by_dataset,
            key=lambda key: (
                -any(row["gold_candidate"] for row in by_dataset[key]),
                key[0] == "GEO",
                key[0], key[1],
            ),
        )
        attempted = 0
        for key in ranked:
            if attempted >= resolve_limit:
                break
            if key in official and official[key].get("status") == "resolved":
                continue
            official[key] = resolve_official(*key)
            attempted += 1
            time.sleep(0.12)
        official_bytes = b"".join(canonical_bytes(official[key]) for key in sorted(official))
        OFFICIAL_RECEIPTS.write_bytes(official_bytes)

    manifests = []
    for (provider, accession), rows in sorted(by_dataset.items()):
        receipt = official.get((provider, accession))
        if receipt is None:
            receipt = {"status": "not_attempted", "official_url": official_endpoint(provider, accession), "http_status": None, "retrieved_at": None, "response_sha256": None, "error": None, "metadata": {}}
        metadata = receipt.get("metadata") or {}
        license_id = metadata.get("license_identifier")
        manifest = {
            "schema": "aleph.external_training.dataset_manifest.v1",
            "authority_status": "proposed",
            "manifest_id": "dataset:" + sha(f"{provider}\x1f{accession}".encode("utf-8"))[:16],
            "leakage_component_id": rows[0]["leakage_component_id"],
            "provider": provider,
            "accession": accession,
            "identifier_type": rows[0]["identifier_type"],
            "canonical_url": rows[0]["canonical_url"],
            "access_state": "public_metadata_resolved" if receipt["status"] == "resolved" else "unavailable" if receipt["status"] == "not_found" else "public_metadata_unresolved",
            "title": metadata.get("title"),
            "license": {"state": "explicit_open" if license_id else "not_exposed", "identifier": license_id, "source": receipt.get("official_url") if license_id else None},
            "file_count": metadata.get("file_count"),
            "total_size_bytes": metadata.get("total_size_bytes"),
            "modality_hints": sorted({hint for row in rows for hint in row["modality_hints"]}),
            "source_count": len({row["source_family_id"] for row in rows}),
            "source_family_ids": sorted({row["source_family_id"] for row in rows}),
            "gold_candidate_source_count": len({row["source_family_id"] for row in rows if row["gold_candidate"]}),
            "metadata_receipt": {key: receipt.get(key) for key in ("status", "official_url", "http_status", "retrieved_at", "response_sha256", "error")},
            "pilot_files": [
                {key: item[key] for key in ("file_name", "bytes", "sha256", "license_identifier", "download_url")}
                for item in sorted(pilot_by_dataset.get((provider, accession), []), key=lambda value: value["file_name"])
            ],
            "limitations": ["accession co-membership does not establish identical processed data", "metadata resolution does not validate scientific content"],
        }
        manifests.append(manifest)
    manifest_bytes = b"".join(canonical_bytes(row) for row in manifests)
    (HERE / "dataset_manifests.jsonl").write_bytes(manifest_bytes)

    target_terms = {"immunofluorescence", "western_blot", "pcr", "piv", "traction_force_microscopy", "afm", "raw_imaging", "flow_cytometry", "proteomics", "omics", "sequencing"}
    queue = []
    for manifest in manifests:
        relevant = sorted(set(manifest["modality_hints"]) & target_terms)
        score = 100 * (manifest["gold_candidate_source_count"] > 0) + 10 * len(relevant) + min(manifest["source_count"], 9)
        queue.append({
            "priority_score": score,
            "manifest_id": manifest["manifest_id"],
            "provider": manifest["provider"],
            "accession": manifest["accession"],
            "access_state": manifest["access_state"],
            "gold_candidate_source_count": manifest["gold_candidate_source_count"],
            "source_count": manifest["source_count"],
            "target_modality_hints": relevant,
            "download_policy": "metadata_then_manual_file_selection",
            "raw_download_status": "not_attempted",
        })
    queue.sort(key=lambda row: (-row["priority_score"], row["provider"], row["accession"]))
    with (HERE / "acquisition_queue.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(queue[0]) if queue else ["priority_score"],
            lineterminator="\n",
        )
        writer.writeheader()
        for row in queue:
            writer.writerow({**row, "target_modality_hints": ";".join(row["target_modality_hints"])})

    source_with = {row["source_family_id"] for row in binding_rows}
    gold_with = source_with & gold
    provider_counts = Counter(row["provider"] for row in manifests)
    modality_counts = Counter(hint for row in manifests for hint in row["modality_hints"])
    resolution_counts = Counter(row["access_state"] for row in manifests)
    summary = {
        "schema": "aleph.external_training.dataset_registry_summary.v1",
        "authority_status": "proposed",
        "corpus_sources": len(receipts),
        "corpus_sources_with_accession": len(source_with),
        "gold_candidate_sources": len(gold),
        "gold_candidate_sources_with_accession": len(gold_with),
        "source_accession_bindings": len(binding_rows),
        "unique_dataset_manifests": len(manifests),
        "provider_counts": dict(sorted(provider_counts.items())),
        "modality_hint_counts": dict(sorted(modality_counts.items())),
        "resolution_counts": dict(sorted(resolution_counts.items())),
        "pilot": {
            "byte_cap": 5_000_000,
            "files_downloaded": len(pilot_rows),
            "bytes_downloaded": sum(row["bytes"] for row in pilot_rows),
            "reason": "fixed four-dataset pilot; every file has an explicit open licence, official API URL, declared size, and local digest",
        },
        "registry_sha256": sha(binding_bytes),
        "manifests_sha256": sha(manifest_bytes),
        "inputs": {
            "first_receipts_sha256": sha(FIRST_RECEIPTS.read_bytes()),
            "second_receipts_sha256": sha(SECOND_RECEIPTS.read_bytes()),
            "gold_candidates_sha256": sha(GOLD.read_bytes()),
        },
        "limitations": [
            "regex matches are exact identifiers but do not prove the current paper generated the dataset",
            "same accession is a leakage component, not proof of identical processing or samples",
            "official metadata APIs expose heterogeneous and often incomplete size, file, and licence fields",
            "no raw dataset payload was committed or used as training evidence",
        ],
    }
    (HERE / "summary.json").write_bytes(canonical_bytes(summary))
    hashes = {name: sha((HERE / name).read_bytes()) for name in ("acquisition_queue.csv", "dataset_manifests.jsonl", "source_accession_registry.jsonl", "summary.json")}
    if OFFICIAL_RECEIPTS.exists(): hashes[OFFICIAL_RECEIPTS.name] = sha(OFFICIAL_RECEIPTS.read_bytes())
    if PILOT_RECEIPTS.exists(): hashes[PILOT_RECEIPTS.name] = sha(PILOT_RECEIPTS.read_bytes())
    (HERE / "SHA256SUMS.json").write_bytes(canonical_bytes(hashes))
    print(json.dumps(summary, sort_keys=True))
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolve-official", type=int, default=0, metavar="LIMIT")
    args = parser.parse_args(argv)
    build(resolve_limit=max(args.resolve_official, 0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
