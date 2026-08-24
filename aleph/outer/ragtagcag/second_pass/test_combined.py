from __future__ import annotations

import hashlib
import json
from pathlib import Path

from aleph.harness import DataRoot
from aleph.outer.ragtagcag.second_pass.combined import run_combined_validation


def _jsonl(path: Path, rows: list[dict]) -> Path:
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
    return path


def _source(family: str, **extra):
    return {
        "source_family_id": family,
        "title": family,
        "authority_status": "proposed",
        "snapshot_date": "2026-08-05",
        "is_open_access_provider_flag": True,
    } | extra


def test_combined_dedupe_failures_payloads_and_idempotence(tmp_path: Path) -> None:
    initial = _jsonl(
        tmp_path / "initial.jsonl",
        [_source("doi:10.1/A"), _source("doi:10.1/B", publication_types=["Retracted Publication"])],
    )
    second = _jsonl(
        tmp_path / "second.jsonl",
        [_source("DOI:10.1/a"), _source("doi:10.1/C")],
    )
    payload = tmp_path / "paper.xml"
    payload.write_bytes(b"<article>combined</article>")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    receipts = _jsonl(
        tmp_path / "receipts.jsonl",
        [
            {
                "source_family_id": "doi:10.1/a",
                "status": "success",
                "sha256": digest,
                "object_path": str(payload),
            },
            {
                "source_family_id": "doi:10.1/c",
                "status": "failed",
                "http_status": 404,
            },
        ],
    )
    result = run_combined_validation(
        initial_manifests=[initial],
        second_pass_manifest=second,
        receipt_paths=[receipts],
        data_root=DataRoot.for_testing(tmp_path / "state"),
    )
    assert result["families"] == {
        "initial_unique": 2,
        "second_pass_candidate_unique": 2,
        "overlap": 1,
        "new_unique": 1,
        "combined_unique": 3,
    }
    assert result["object_audit"]["receipt_success"] == 1
    assert result["object_audit"]["receipt_failed"] == 1
    assert result["object_audit"]["payloads_present_and_verified"] == 1
    assert result["first"]["retracted_entities"] == 1
    assert result["second"]["events_appended"] == 0
    assert result["first"]["snapshot_manifest_hash"] == result["second"]["snapshot_manifest_hash"]
    assert result["first"]["cag_lint_violations"] == 0
