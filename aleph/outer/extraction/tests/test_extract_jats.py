from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

MODULE_PATH = Path(__file__).parents[1] / "extract_jats.py"
SPEC = importlib.util.spec_from_file_location("extract_jats", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


XML = b"""<article><front><article-meta><abstract><p>AFM measured cell stiffness in fibroblasts.</p></abstract></article-meta></front><body>
<sec id="s1" sec-type="methods"><title>Methods</title><p>Cells were treated with an inhibitor.</p><fig id="f1"><caption><p>Atomic force microscopy image.</p></caption></fig></sec>
<sec id="s2"><title>Results</title><p>Cortical tension changed.</p><table-wrap id="t1"><caption><p>Traction force values.</p></caption></table-wrap></sec>
</body></article>"""


def receipt(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(XML).hexdigest()
    return {
        "schema": "aleph.external_training.oa_acquisition_receipt.v1",
        "status": "success",
        "source_family_id": "doi:10.0/example",
        "pmcid": "PMC123",
        "sha256": digest,
        "object_path": str(path),
    }


def test_extracts_locators_counts_tags_and_hashes(tmp_path: Path) -> None:
    object_path = tmp_path / "object"
    object_path.write_bytes(XML)
    vocabulary = MODULE.load_vocabulary(Path(__file__).parents[1] / "tag_vocabulary.json")
    records, counts = MODULE.extract_payload(XML, receipt(object_path), vocabulary)
    assert counts == {"sections": 3, "figures": 1, "tables": 1}
    assert {record["section_class"] for record in records} >= {"abstract", "methods", "results", "figure_caption", "table_caption"}
    assert any("afm" in record["tags"]["measurement_modality"] for record in records)
    assert all(len(record["chunk_sha256"]) == 64 and len(record["text_sha256"]) == 64 for record in records)
    assert all(record["authority_status"] == "proposed" for record in records)


def test_process_is_resumable_and_repeat_hash_is_deterministic(tmp_path: Path) -> None:
    object_path = tmp_path / "object"
    object_path.write_bytes(XML)
    receipts = tmp_path / "receipts.jsonl"
    receipts.write_text(json.dumps(receipt(object_path)) + "\n", encoding="utf-8")
    derived = tmp_path / "derived"
    vocabulary_path = Path(__file__).parents[1] / "tag_vocabulary.json"
    first = MODULE.process(receipts, tmp_path, derived, vocabulary_path)
    second = MODULE.process(receipts, tmp_path, derived, vocabulary_path)
    assert first["failures"] == 0
    assert second["resumed_articles"] == 1
    assert first["article_manifest_set_sha256"] == second["article_manifest_set_sha256"]
    assert not list(Path(__file__).parents[1].glob("*.chunks.jsonl"))


def test_hash_mismatch_fails_closed(tmp_path: Path) -> None:
    object_path = tmp_path / "object"
    object_path.write_bytes(XML + b"corrupt")
    receipts = tmp_path / "receipts.jsonl"
    receipts.write_text(json.dumps(receipt(object_path)) + "\n", encoding="utf-8")
    summary = MODULE.process(receipts, tmp_path, tmp_path / "derived", Path(__file__).parents[1] / "tag_vocabulary.json")
    assert summary["xml_parsed"] == 0
    assert summary["failures"] == 1
    assert summary["failure_classes"] == {"ValueError": 1}
