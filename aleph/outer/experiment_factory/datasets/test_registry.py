import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import build_registry
import download_pilot


HERE = Path(__file__).parent
SHA = re.compile(r"^[0-9a-f]{64}$")


def _jsonl(name):
    return [json.loads(line) for line in (HERE / name).read_text(encoding="utf-8").splitlines() if line.strip()]


def test_identifier_controls_are_normalized_without_substring_false_positives():
    assert ("GEO", "geo_series", "GSE12345") in build_registry.matches("data: GSE12345")
    assert ("BioProject", "bioproject", "PRJNA123456") in build_registry.matches("PRJNA123456")
    assert ("Zenodo", "record", "10.5281/zenodo.42") in build_registry.matches("https://zenodo.org/records/42")
    assert build_registry.matches("AGSE12345 and PRJNAx123 are not accessions") == []


def test_registry_has_exact_digests_and_no_duplicate_binding():
    rows = _jsonl("source_accession_registry.jsonl")
    keys = [(row["source_family_id"], row["provider"], row["accession"]) for row in rows]
    assert len(keys) == len(set(keys))
    assert all(row["authority_status"] == "proposed" for row in rows)
    for row in rows:
        assert SHA.fullmatch(row["payload_sha256"])
        assert row["evidence"]
        assert all(SHA.fullmatch(item["document_sha256"]) and SHA.fullmatch(item["text_sha256"]) for item in row["evidence"])


def test_manifests_are_one_per_provider_accession_and_never_promoted():
    rows = _jsonl("dataset_manifests.jsonl")
    keys = [(row["provider"], row["accession"]) for row in rows]
    assert len(keys) == len(set(keys))
    assert all(row["authority_status"] == "proposed" for row in rows)
    assert all(row["source_count"] == len(row["source_family_ids"]) for row in rows)
    assert {"GEO", "BioProject", "ArrayExpress", "PRIDE", "Zenodo", "Figshare", "Dryad", "OSF"} <= {row["provider"] for row in rows}


def test_summary_accounts_for_full_corpus_and_gold_candidates():
    summary = json.loads((HERE / "summary.json").read_text(encoding="utf-8"))
    rows = _jsonl("source_accession_registry.jsonl")
    manifests = _jsonl("dataset_manifests.jsonl")
    assert summary["corpus_sources"] == 5675
    assert summary["gold_candidate_sources"] == 300
    assert summary["source_accession_bindings"] == len(rows)
    assert summary["unique_dataset_manifests"] == len(manifests)
    assert summary["corpus_sources_with_accession"] == len({row["source_family_id"] for row in rows})
    assert summary["gold_candidate_sources_with_accession"] == len({row["source_family_id"] for row in rows if row["gold_candidate"]})


def test_acquisition_queue_is_sorted_and_preserves_typed_download_state():
    with (HERE / "acquisition_queue.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ordering = [(-int(row["priority_score"]), row["provider"], row["accession"]) for row in rows]
    assert ordering == sorted(ordering)
    assert {row["raw_download_status"] for row in rows} == {"not_attempted"}


def test_pilot_is_bounded_open_and_receipt_only():
    summary = json.loads((HERE / "summary.json").read_text(encoding="utf-8"))
    receipts = _jsonl("pilot_download_receipts.jsonl")
    assert summary["pilot"]["files_downloaded"] == len(receipts)
    assert summary["pilot"]["bytes_downloaded"] == sum(row["bytes"] for row in receipts)
    assert summary["pilot"]["bytes_downloaded"] <= summary["pilot"]["byte_cap"]
    assert all(row["storage"] == "ignored_content_addressed_local_object" for row in receipts)
    assert all(SHA.fullmatch(row["sha256"]) for row in receipts)
    assert all(row["license_identifier"].lower() in download_pilot.OPEN_LICENSES for row in receipts)


def test_checksum_manifest_matches_all_named_artifacts():
    hashes = json.loads((HERE / "SHA256SUMS.json").read_text(encoding="utf-8"))
    for name, expected in hashes.items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == expected
