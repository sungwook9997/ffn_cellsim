from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("audit_second_pass.py")
SPEC = importlib.util.spec_from_file_location("audit_second_pass", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def canonical(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()


def test_audit_accepts_replay_with_only_resume_count_changed(tmp_path: Path) -> None:
    derived = tmp_path / "derived" / "articles"
    derived.mkdir(parents=True)
    text = "synthetic retrieval record"
    record = {
        "schema": "aleph.external_training.jats_chunk.v1",
        "authority_status": "proposed",
        "source_family_id": "doi:10.0/example",
        "pmcid": "PMC1",
        "payload_sha256": "a" * 64,
        "article_locator": "pmcid:PMC1",
        "section_locator": "abstract:0",
        "section_title": "Abstract",
        "section_class": "abstract",
        "chunk_ordinal": 0,
        "text_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "text": text,
        "tags": {"cell_type_state": [], "measurement_modality": [], "mechanics_observable": [], "perturbation": []},
    }
    basis = dict(record)
    basis.pop("text")
    record["chunk_sha256"] = hashlib.sha256(canonical(basis)).hexdigest()
    chunk_bytes = canonical(record) + b"\n"
    (derived / "one.chunks.jsonl").write_bytes(chunk_bytes)
    manifest = {"chunks": 1, "chunks_file_sha256": hashlib.sha256(chunk_bytes).hexdigest()}
    (derived / "one.manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    first = {"article_manifest_set_sha256": "b" * 64, "resumed_articles": 0}
    repeat = {"article_manifest_set_sha256": "b" * 64, "resumed_articles": 1}
    first_path = tmp_path / "first.json"
    repeat_path = tmp_path / "repeat.json"
    first_path.write_text(json.dumps(first), encoding="utf-8")
    repeat_path.write_text(json.dumps(repeat), encoding="utf-8")
    result = MODULE.audit(first_path, repeat_path, tmp_path / "derived")
    assert result["summaries_equal_excluding_resume_count"] is True
    assert result["manifest_set_hash_equal"] is True
    assert result["chunk_records"] == 1
    assert result["audit_errors"] == 0


def test_audit_fails_on_text_hash_corruption(tmp_path: Path) -> None:
    derived = tmp_path / "derived" / "articles"
    derived.mkdir(parents=True)
    record = {"text": "changed", "text_sha256": "0" * 64, "chunk_sha256": "0" * 64, "authority_status": "proposed"}
    chunk_bytes = canonical(record) + b"\n"
    (derived / "one.chunks.jsonl").write_bytes(chunk_bytes)
    (derived / "one.manifest.json").write_text(json.dumps({"chunks": 1, "chunks_file_sha256": hashlib.sha256(chunk_bytes).hexdigest()}), encoding="utf-8")
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({"article_manifest_set_sha256": "a" * 64}), encoding="utf-8")
    result = MODULE.audit(summary, summary, tmp_path / "derived")
    assert result["audit_errors"] >= 2
    assert result["error_counts"]["bad_text_hashes"] == 1
