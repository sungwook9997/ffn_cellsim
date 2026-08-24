from __future__ import annotations

import json
from pathlib import Path

import pytest

from aleph.harness import DataRoot, Ledger, ObjectStore
from aleph.harness.snapshot import fold_entities, resolve
from aleph.outer.ragtagcag.ingest import (
    CorpusIngestError,
    ingest_corpus,
    verify_idempotent_ingest,
)


def _write_jsonl(path: Path, rows: list[dict]) -> Path:
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
    return path


def _row(family: str, **extra):
    return {
        "source_family_id": family,
        "title": family,
        "authority_status": "proposed",
        "snapshot_date": "2026-08-05",
        "is_open_access_provider_flag": True,
    } | extra


def test_deduplicates_by_source_family_and_is_idempotent(tmp_path: Path) -> None:
    one = _write_jsonl(tmp_path / "one.jsonl", [_row("doi:10.1/A"), _row("doi:10.1/B")])
    two = _write_jsonl(tmp_path / "two.jsonl", [_row("DOI:10.1/a", journal="J")])
    root = DataRoot.for_testing(tmp_path / "state")

    result = verify_idempotent_ingest([one, two], root)

    assert result["first"]["input_records"] == 3
    assert result["first"]["source_families"] == 2
    assert result["first"]["duplicate_records"] == 1
    assert result["first"]["events_appended"] == 2
    assert result["second"]["events_appended"] == 0
    assert result["first"]["snapshot_manifest_hash"] == result["second"]["snapshot_manifest_hash"]


def test_preserves_proposed_access_and_license_without_promotion(tmp_path: Path) -> None:
    manifest = _write_jsonl(
        tmp_path / "sources.jsonl",
        [
            _row(
                "doi:10.1/licensed",
                access_route="institution_kaist",
                license_status="publisher_subscription_terms_unverified",
                redistribution="local_only",
            )
        ],
    )
    root = DataRoot.for_testing(tmp_path / "state")
    result = ingest_corpus([manifest], root)
    snapshot = resolve(root, result.snapshot_id)
    entity = next(iter(snapshot.entities.values()))
    assert entity.attrs["authority_status"] == "proposed"
    assert entity.attrs["status"] == "proposed"
    assert entity.attrs["access_route"] == "institution_kaist"
    assert entity.attrs["license_status"] == "publisher_subscription_terms_unverified"
    assert entity.attrs["redistribution"] == "local_only"
    assert "accepted" not in entity.event_types


def test_receipt_and_payload_are_content_addressed_and_present(tmp_path: Path) -> None:
    manifest = _write_jsonl(tmp_path / "sources.jsonl", [_row("doi:10.1/payload")])
    payload = tmp_path / "paper.xml"
    payload.write_bytes(b"<article>evidence</article>")
    import hashlib

    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    receipt = _write_jsonl(
        tmp_path / "receipts.jsonl",
        [
            {
                "source_family_id": "doi:10.1/payload",
                "payload_sha256": digest,
                "payload_path": str(payload),
                "access_route": "open_access",
                "license_status": "cc by",
                "redistribution": "permitted_by_recorded_license",
            }
        ],
    )
    root = DataRoot.for_testing(tmp_path / "state")
    result = ingest_corpus([manifest], root, receipt_paths=[receipt])
    assert result.acquisition_receipts == 1
    assert result.acquisition_payloads_present == 1
    assert result.missing_objects == 0
    assert ObjectStore.open(root).verify(digest)


def test_explicit_retraction_is_append_only_state(tmp_path: Path) -> None:
    manifest = _write_jsonl(tmp_path / "sources.jsonl", [_row("doi:10.1/retracted")])
    root = DataRoot.for_testing(tmp_path / "state")
    first = ingest_corpus([manifest], root)
    _write_jsonl(manifest, [_row("doi:10.1/retracted", quality_status="retracted")])
    second = ingest_corpus([manifest], root)
    with Ledger.open(root) as ledger:
        entity = next(iter(fold_entities(ledger).values()))
        assert ledger.count() == 2
        assert entity.retracted is True
        assert entity.event_types[-1].endswith(".retracted")
    assert first.ledger_events == 1
    assert second.ledger_events == 2

    # A later input that merely omits the retraction cannot silently resurrect
    # the source; reinstatement would require an authority workflow the public
    # snapshot fold intentionally does not expose.
    _write_jsonl(manifest, [_row("doi:10.1/retracted")])
    third = ingest_corpus([manifest], root)
    with Ledger.open(root) as ledger:
        entity = next(iter(fold_entities(ledger).values()))
        assert entity.retracted is True
        assert entity.event_types[-1].endswith(".retracted")
    assert third.ledger_events == 3


def test_refuses_authority_promotion(tmp_path: Path) -> None:
    manifest = _write_jsonl(
        tmp_path / "sources.jsonl", [_row("doi:10.1/promoted", authority_status="accepted")]
    )
    with pytest.raises(CorpusIngestError, match="only proposed imports"):
        ingest_corpus([manifest], DataRoot.for_testing(tmp_path / "state"))


def test_unavailable_payload_digest_is_preserved_but_not_referenced(tmp_path: Path) -> None:
    manifest = _write_jsonl(tmp_path / "sources.jsonl", [_row("doi:10.1/missing")])
    receipt = _write_jsonl(
        tmp_path / "receipts.jsonl",
        [{"source_family_id": "doi:10.1/missing", "payload_sha256": "a" * 64}],
    )
    root = DataRoot.for_testing(tmp_path / "state")
    result = ingest_corpus([manifest], root, receipt_paths=[receipt])
    assert result.acquisition_payloads_unavailable == 1
    assert result.missing_objects == 0
