#!/usr/bin/env python3
"""Independent, text-free verification of the two-pass external OA corpus."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


INITIAL_CANDIDATES = (
    "corpus/external_training/snapshots/europe_pmc_2026-08-05.jsonl",
    "corpus/external_training/snapshots/europe_pmc_cell_states_2010s.jsonl",
    "corpus/external_training/snapshots/europe_pmc_cell_states_2026-08-05.jsonl",
    "corpus/external_training/snapshots/europe_pmc_expansion_2026-08-05.jsonl",
    "corpus/external_training/snapshots/europe_pmc_historical_2010_2021.jsonl",
    "corpus/external_training/sources/seed_primary_sources.jsonl",
)
SECOND_CANDIDATE = "corpus/external_training/acquisition/second_pass/candidates_2026-08-05.jsonl"
RECEIPTS = (
    "corpus/external_training/acquisition/europe_pmc_oa_2026-08-05.jsonl",
    "corpus/external_training/acquisition/second_pass/receipts_2026-08-05.jsonl",
)
DERIVED = (
    "data/external_training/derived",
    "data/external_training/derived_second_pass",
)
COMBINED_INDEX = "corpus/external_training/retrieval/combined/index.jsonl"
COMBINED_RECEIPT = "corpus/external_training/retrieval/combined/receipt.json"
COMBINED_MODEL_ROOT = "corpus/external_training/model/combined_content"
COMBINED_RAG_RECEIPT = "corpus/external_training/ragtagcag/second_pass/combined_ingest_2026-08-05.json"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def records(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def family(value: Any) -> str:
    text = str(value or "").strip()
    prefix, separator, identifier = text.partition(":")
    if separator and prefix.lower() in {"doi", "pmid", "pmcid", "arxiv"}:
        return f"{prefix.lower()}:{identifier.strip().lower()}"
    if not text:
        raise ValueError("missing source_family_id")
    return text


def families(paths: Iterable[Path]) -> tuple[int, set[str]]:
    rows = [row for path in paths for row in records(path)]
    return len(rows), {family(row.get("source_family_id")) for row in rows}


def audit_receipts(root: Path, paths: tuple[Path, ...]) -> dict[str, Any]:
    pass_results = []
    all_success_families: list[set[str]] = []
    all_object_paths: set[Path] = set()
    totals = Counter()
    errors = Counter()
    for path in paths:
        rows = records(path)
        statuses = Counter(str(row.get("status", "unknown")) for row in rows)
        success_rows = [row for row in rows if row.get("status") == "success"]
        success_families = {family(row.get("source_family_id")) for row in success_rows}
        success_payloads = {str(row.get("sha256")) for row in success_rows}
        verified_bytes = 0
        object_paths: set[Path] = set()
        for row in success_rows:
            object_path = Path(str(row.get("object_path", "")))
            resolved = object_path if object_path.is_absolute() else root / object_path
            object_paths.add(resolved.resolve())
            if not resolved.is_file():
                errors["missing_success_object"] += 1
                continue
            digest = sha256_path(resolved)
            if digest != row.get("sha256"):
                errors["success_object_digest_mismatch"] += 1
            size = resolved.stat().st_size
            if size != int(row.get("bytes", -1)):
                errors["success_object_byte_mismatch"] += 1
            verified_bytes += size
        all_success_families.append(success_families)
        all_object_paths |= object_paths
        totals.update(records=len(rows), success=statuses["success"], failed=statuses["failed"], bytes=verified_bytes)
        pass_results.append({
            "path": str(path.relative_to(root)),
            "sha256": sha256_path(path),
            "records": len(rows),
            "status_counts": dict(sorted(statuses.items())),
            "success_families": len(success_families),
            "success_payload_digests": len(success_payloads),
            "verified_success_objects": len(object_paths),
            "verified_bytes": verified_bytes,
        })
    object_root = root / "data/external_training/objects"
    disk_objects = {path.resolve() for path in object_root.iterdir() if path.is_file()}
    return {
        "passes": pass_results,
        "totals": dict(totals),
        "success_family_overlap_between_passes": len(all_success_families[0] & all_success_families[1]),
        "referenced_unique_objects": len(all_object_paths),
        "disk_object_files": len(disk_objects),
        "unreferenced_disk_objects": len(disk_objects - all_object_paths),
        "referenced_objects_missing_from_disk_listing": len(all_object_paths - disk_objects),
        "errors": dict(sorted(errors.items())),
    }


def audit_derived(root: Path, stores: tuple[Path, ...]) -> dict[str, Any]:
    pass_results = []
    totals = Counter()
    errors = Counter()
    family_sets: list[set[str]] = []
    for store in stores:
        article_root = store / "articles"
        manifest_paths = sorted(article_root.glob("*.manifest.json"))
        chunk_paths = sorted(article_root.glob("*.chunks.jsonl"))
        chunk_by_key = {path.name.removesuffix(".chunks.jsonl"): path for path in chunk_paths}
        manifests = []
        families_in_store: set[str] = set()
        chunks = sections = figures = tables = payload_bytes = records_seen = 0
        for manifest_path in manifest_paths:
            manifest = json.loads(manifest_path.read_text("utf-8"))
            manifests.append(manifest)
            families_in_store.add(family(manifest.get("source_family_id")))
            chunks += int(manifest.get("chunks", 0))
            sections += int(manifest.get("sections", 0))
            figures += int(manifest.get("figures", 0))
            tables += int(manifest.get("tables", 0))
            payload_bytes += int(manifest.get("payload_bytes", 0))
            key = manifest_path.name.removesuffix(".manifest.json")
            chunk_path = chunk_by_key.get(key)
            if chunk_path is None:
                errors["missing_chunk_file"] += 1
                continue
            chunk_data = chunk_path.read_bytes()
            if sha256_bytes(chunk_data) != manifest.get("chunks_file_sha256"):
                errors["chunk_file_digest_mismatch"] += 1
            local_count = 0
            for line in chunk_data.splitlines():
                if not line:
                    continue
                local_count += 1
                record = json.loads(line)
                if sha256_bytes(record["text"].encode("utf-8")) != record.get("text_sha256"):
                    errors["text_digest_mismatch"] += 1
                basis = dict(record)
                expected = basis.pop("chunk_sha256", "")
                basis.pop("text", None)
                if sha256_bytes(canonical_bytes(basis)) != expected:
                    errors["chunk_record_digest_mismatch"] += 1
                if record.get("authority_status") != "proposed":
                    errors["non_proposed_chunk"] += 1
            records_seen += local_count
            if local_count != int(manifest.get("chunks", -1)):
                errors["manifest_chunk_count_mismatch"] += 1
        manifest_set_hash = sha256_bytes(b"".join(canonical_bytes(item) + b"\n" for item in sorted(manifests, key=lambda item: item["source_family_id"])))
        orphan_chunks = len(set(chunk_by_key) - {path.name.removesuffix(".manifest.json") for path in manifest_paths})
        errors["orphan_chunk_file"] += orphan_chunks
        family_sets.append(families_in_store)
        result = {
            "path": str(store.relative_to(root)),
            "article_manifests": len(manifest_paths),
            "chunk_files": len(chunk_paths),
            "source_families": len(families_in_store),
            "chunk_records": records_seen,
            "manifest_chunk_sum": chunks,
            "sections": sections,
            "figures": figures,
            "tables": tables,
            "payload_bytes": payload_bytes,
            "article_manifest_set_sha256": manifest_set_hash,
        }
        pass_results.append(result)
        totals.update(families=len(families_in_store), chunks=records_seen, sections=sections, figures=figures, tables=tables, payload_bytes=payload_bytes)
    return {
        "passes": pass_results,
        "totals": dict(totals),
        "source_family_overlap_between_stores": len(family_sets[0] & family_sets[1]),
        "errors": {key: value for key, value in sorted(errors.items()) if value},
    }


def tracked_hygiene(root: Path) -> dict[str, Any]:
    completed = subprocess.run(["git", "ls-files", "-z"], cwd=root, check=True, capture_output=True)
    tracked = [item.decode("utf-8") for item in completed.stdout.split(b"\0") if item]
    forbidden = [
        path for path in tracked
        if path.startswith("data/external_training/")
        or path.endswith(".chunks.jsonl")
        or Path(path).suffix.lower() in {".pdf", ".xml", ".nxml", ".jats"}
    ]
    return {
        "tracked_files": len(tracked),
        "tracked_external_data_root_files": sum(path.startswith("data/external_training/") for path in tracked),
        "tracked_chunk_text_files": sum(path.endswith(".chunks.jsonl") for path in tracked),
        "tracked_pdf_xml_jats_files": sum(Path(path).suffix.lower() in {".pdf", ".xml", ".nxml", ".jats"} for path in tracked),
        "forbidden_tracked_paths": forbidden,
    }


def last_commit(root: Path, path: str) -> str:
    completed = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", path],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def audit_combined_retrieval(root: Path) -> dict[str, Any]:
    index_path = root / COMBINED_INDEX
    receipt_path = root / COMBINED_RECEIPT
    receipt = json.loads(receipt_path.read_text("utf-8"))
    rows = records(index_path)
    row_families = {family(row.get("source_family_id")) for row in rows}
    routable_families = {
        family(row.get("source_family_id")) for row in rows if row.get("neural_routable") is True
    }
    pass_counts = Counter(str(row.get("pass")) for row in rows)
    errors = Counter()
    forbidden_content_fields = {"text", "snippet", "full_text", "abstract", "body"}
    forbidden_fields_seen: set[str] = set()
    derived_roots = json.loads(
        (root / "corpus/external_training/retrieval/combined/config.json").read_text("utf-8")
    )["derived_roots"]
    chunks = 0
    for row in rows:
        forbidden_fields_seen |= forbidden_content_fields & set(row)
        if row.get("authority_status") != "proposed":
            errors["non_proposed_index_row"] += 1
        chunks += int(row.get("chunks", 0))
        pass_name = str(row.get("pass"))
        derived_root = root / derived_roots[pass_name]
        chunk_path = derived_root / str(row.get("chunks_file"))
        if not chunk_path.is_file():
            errors["missing_indexed_chunk_file"] += 1
        elif sha256_path(chunk_path) != row.get("chunks_file_sha256"):
            errors["indexed_chunk_file_digest_mismatch"] += 1

    model_root = root / COMBINED_MODEL_ROOT
    model_hashes = {
        "weights_sha256": sha256_path(model_root / "weights_combined_content.npz"),
        "splits_sha256": sha256_path(model_root / "splits.jsonl"),
        "coarse_taxonomy_sha256": sha256_path(model_root / "taxonomy.json"),
    }
    split_rows = records(model_root / "splits.jsonl")
    split_families = {family(row.get("source_family_id")) for row in split_rows}
    rag_path = root / COMBINED_RAG_RECEIPT
    rag = json.loads(rag_path.read_text("utf-8"))
    rag_hash = sha256_path(rag_path)
    expected_rag = receipt["ragtagcag"]
    binding_checks = {
        "index_sha256": sha256_path(index_path) == receipt.get("index_sha256"),
        "weights_sha256": model_hashes["weights_sha256"] == receipt["model"].get("weights_sha256"),
        "splits_sha256": model_hashes["splits_sha256"] == receipt["model"].get("splits_sha256"),
        "coarse_taxonomy_sha256": model_hashes["coarse_taxonomy_sha256"] == receipt["model"].get("coarse_taxonomy_sha256"),
        "rag_validation_receipt_sha256": rag_hash == expected_rag.get("validation_receipt_sha256"),
        "rag_ledger_head_hash": rag["first"].get("ledger_head_hash") == expected_rag.get("ledger_head_hash"),
        "rag_snapshot_manifest_hash": rag["first"].get("snapshot_manifest_hash") == expected_rag.get("snapshot_manifest_hash"),
        "rag_projection_hash": rag["first"].get("projection_hash") == expected_rag.get("projection_hash"),
        "model_routable_set_matches_index": split_families == routable_families,
    }
    producer_commits = {
        "model_weights": last_commit(root, f"{COMBINED_MODEL_ROOT}/weights_combined_content.npz"),
        "rag_validation_receipt": last_commit(root, COMBINED_RAG_RECEIPT),
        "combined_index": last_commit(root, COMBINED_INDEX),
    }
    claims_verified = (
        len(rows) == len(row_families) == receipt.get("local_source_families") == 5675
        and len(routable_families) == receipt["model"].get("neural_routable_source_families") == 5669
        and chunks == receipt.get("chunk_count") == 284620
        and dict(sorted(pass_counts.items())) == {"first": 2701, "second": 2974}
        and not forbidden_fields_seen
        and receipt.get("text_committed") is False
        and not errors
        and all(binding_checks.values())
        and producer_commits["model_weights"].startswith("5336672")
        and producer_commits["rag_validation_receipt"].startswith("736153b")
        and producer_commits["combined_index"].startswith("3dbb3aa")
    )
    return {
        "index_rows": len(rows),
        "unique_source_families": len(row_families),
        "neural_routable_source_families": len(routable_families),
        "chunk_count": chunks,
        "pass_source_families": dict(sorted(pass_counts.items())),
        "index_sha256": sha256_path(index_path),
        "forbidden_content_fields": sorted(forbidden_fields_seen),
        "text_committed": receipt.get("text_committed"),
        "binding_checks": binding_checks,
        "producer_commits": producer_commits,
        "errors": dict(sorted(errors.items())),
        "claims_verified": claims_verified,
    }


def run(root: Path) -> dict[str, Any]:
    initial_rows, initial_families = families(tuple(root / path for path in INITIAL_CANDIDATES))
    second_rows, second_families = families((root / SECOND_CANDIDATE,))
    candidate = {
        "first_pass_files": len(INITIAL_CANDIDATES),
        "first_pass_rows": initial_rows,
        "first_pass_unique_families": len(initial_families),
        "second_pass_rows": second_rows,
        "second_pass_unique_families": len(second_families),
        "family_overlap": len(initial_families & second_families),
        "second_pass_new_families": len(second_families - initial_families),
        "combined_unique_families": len(initial_families | second_families),
    }
    receipt = audit_receipts(root, tuple(root / path for path in RECEIPTS))
    derived = audit_derived(root, tuple(root / path for path in DERIVED))
    hygiene = tracked_hygiene(root)
    combined_retrieval = audit_combined_retrieval(root)
    claims_verified = (
        candidate == {
            "first_pass_files": 6,
            "first_pass_rows": 7218,
            "first_pass_unique_families": 7213,
            "second_pass_rows": 3000,
            "second_pass_unique_families": 3000,
            "family_overlap": 5,
            "second_pass_new_families": 2995,
            "combined_unique_families": 10208,
        }
        and receipt["totals"] == {"records": 5713, "success": 5675, "failed": 38, "bytes": 850654786}
        and receipt["success_family_overlap_between_passes"] == 0
        and not receipt["errors"]
        and derived["totals"] == {"families": 5675, "chunks": 284620, "sections": 161989, "figures": 38316, "tables": 4450, "payload_bytes": 850654786}
        and derived["source_family_overlap_between_stores"] == 0
        and not derived["errors"]
        and not hygiene["forbidden_tracked_paths"]
        and combined_retrieval["claims_verified"]
    )
    return {
        "schema": "aleph.external_training.independent_verification.v1",
        "authority_status": "proposed",
        "candidate_families": candidate,
        "receipts_and_objects": receipt,
        "derived_chunks": derived,
        "tracked_hygiene": hygiene,
        "combined_retrieval": combined_retrieval,
        "claims_verified": claims_verified,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.repo_root.resolve())
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if result["claims_verified"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
