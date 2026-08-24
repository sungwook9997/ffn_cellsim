from __future__ import annotations

import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("screen_oa", HERE / "screen_oa.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def chunk(text: str, locator: str = "section:m1", section_class: str = "methods") -> dict:
    digest = MODULE.sha(text.encode())
    return {
        "section_locator": locator,
        "section_class": section_class,
        "chunk_sha256": digest,
        "text_sha256": digest,
        "text": text,
        "tags": {"cell_type_state": [], "measurement_modality": [], "mechanics_observable": [], "perturbation": []},
    }


def test_signal_refs_are_digest_and_locator_only():
    item = chunk("Three independent experiments used n = 24 cells; mean ± SEM was reported in pN.")
    ref = MODULE.evidence_ref(item, "test")
    assert ref["locator"] == "section:m1"
    assert "text" not in ref
    assert set(ref) == {"locator", "section_class", "chunk_sha256", "text_sha256", "pattern_id"}


def test_article_type_requires_explicit_metadata():
    assert MODULE.article_type({"publication_types": ["research-article"]})["status"] == "research_article_metadata"
    assert MODULE.article_type({"publication_types": ["Journal Article"]})["status"] == "unknown"
    assert MODULE.article_type({"publication_types": ["Review"]})["status"] == "non_research_metadata"


def test_screen_never_emits_tier_a(tmp_path: Path):
    chunks = tmp_path / "x.chunks.jsonl"
    rows = [
        chunk("We used n = 24 cells in 3 independent experiments. The force instrument was calibrated. Values in pN are mean ± SEM. Data availability is provided; excluded outliers were predeclared."),
        chunk("Traction force microscopy showed 12 pN ± SEM.", "figure:F1", "figure_caption"),
    ]
    rows[1]["tags"]["measurement_modality"] = ["traction_force_microscopy"]
    rows[1]["tags"]["mechanics_observable"] = ["traction_force"]
    chunks.write_text("".join(json.dumps(row) + "\n" for row in rows))
    manifest = {"source_family_id": "doi:10.test/x", "pmcid": "PMC1", "payload_sha256": "a" * 64, "chunks_file_sha256": MODULE.sha(chunks.read_bytes())}
    record = MODULE.screen_article(manifest, chunks, {"publication_types": ["research-article"]}, {})
    assert record["decision"] == "manual_priority"
    assert "tier_a" not in json.dumps(record).lower()
    assert all("text" not in ref for value in record["signals"].values() for ref in value["evidence"])


def test_exception_is_explicit_and_proposed(tmp_path: Path):
    chunks = tmp_path / "x.chunks.jsonl"
    row = chunk("Optical tweezers measured 10 pN.", "figure:F2", "figure_caption")
    row["tags"]["measurement_modality"] = ["optical_tweezers"]
    chunks.write_text(json.dumps(row) + "\n")
    manifest = {"source_family_id": "doi:10.test/y", "pmcid": "PMC2", "payload_sha256": "b" * 64, "chunks_file_sha256": MODULE.sha(chunks.read_bytes())}
    record = MODULE.screen_article(manifest, chunks, {"publication_types": ["Review"]}, {})
    assert record["decision"] == "exception_candidate"
    assert record["exception_reason"] == "rare_modality"
    assert record["authority_status"] == "proposed"
