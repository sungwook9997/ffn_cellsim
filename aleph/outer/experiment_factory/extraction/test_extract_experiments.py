from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .extract_experiments import (
    extract_article,
    normalize_measure,
    process,
    record_for_unit,
    sha,
)


def receipt(payload: bytes) -> dict:
    return {
        "source_family_id": "doi:10.1/example",
        "pmcid": "PMC1",
        "sha256": hashlib.sha256(payload).hexdigest(),
        "acquisition_pass": "test",
        "receipt_sha256": hashlib.sha256(b"receipt").hexdigest(),
    }


def test_exact_decimal_unit_conversion() -> None:
    assert normalize_measure("0.0963", "pN·s/µm")["normalized_value"] == "0.0963"
    assert normalize_measure("2.5", "kPa")["normalized_value"] == "2500"
    assert normalize_measure("37", "°C")["normalized_value"] == "310.15"


def test_missing_context_and_missing_units_are_not_invented() -> None:
    row = {"source_family_id": "doi:x", "pmcid": None, "sha256": "a" * 64, "acquisition_pass": "test"}
    assert record_for_unit(row, 2020, ["Results"], "section", "s1", "The value was 12.") is None
    assert record_for_unit(row, 2020, ["Results"], "section", "s1", "The pressure was 12 kPa.") is None


def test_background_mortality_is_ineligible() -> None:
    row = {"source_family_id": "doi:x", "pmcid": None, "sha256": "a" * 64, "acquisition_pass": "test"}
    value = "Fibroblasts are implicated in IPF, whose 3-year mortality is 50%."
    assert record_for_unit(row, 2020, ["Background"], "section", "background:1", value) is None


def test_uppercase_s_in_rrna_is_not_seconds() -> None:
    row = {"source_family_id": "doi:x", "pmcid": None, "sha256": "a" * 64, "acquisition_pass": "test"}
    value = "In qPCR methods, fibroblast expression was normalized to 18 S ribosomal RNA."
    assert record_for_unit(row, 2020, ["Methods"], "section", "methods:1", value) is None


def test_explicit_table_and_ambiguity_are_preserved() -> None:
    payload = b"""<article><front><article-meta><pub-date><year>2024</year></pub-date></article-meta></front>
    <body><sec id='m'><title>Methods</title><p>Fibroblasts were measured by atomic force microscopy for 5 min at 37 &#176;C (n=8).</p></sec>
    <table-wrap id='t1'><caption><p>AFM stiffness in fibroblasts</p></caption><table><tr><th>control</th><th>treated</th></tr><tr><td>2.5 kPa</td><td>3.0 kPa</td></tr></table></table-wrap></body></article>"""
    records = extract_article(payload, receipt(payload))
    assert any(r["protocols"][0]["evidence"][0]["locator_id"] == "table:t1/row:1" for r in records)
    methods = next(r for r in records if r["protocols"][0]["evidence"][0]["locator_id"].startswith("section:m"))
    assert all(o["measured_entity"].startswith("protocol_parameter:") for o in methods["observations"])
    assert methods["observations"][0]["replicates"]["unit_of_replication"] == "reported_n=8;type_unspecified"
    # No control/treatment arm pairing is invented from the row layout.
    table = next(r for r in records if r["protocols"][0]["evidence"][0]["locator_type"] == "table")
    assert table["comparisons"] == []


def test_process_is_resumable_and_byte_deterministic(tmp_path: Path) -> None:
    payload = b"<article><body><sec id='r'><title>Results</title><p>AFM of fibroblasts reported 2.5 kPa.</p></sec></body></article>"
    object_path = tmp_path / "object"
    object_path.write_bytes(payload)
    row = {
        **receipt(payload), "schema": "aleph.external_training.oa_acquisition_receipt.v1",
        "status": "success", "object_path": str(object_path),
    }
    receipts = tmp_path / "receipts.jsonl"
    receipts.write_text(json.dumps(row) + "\n")
    output = tmp_path / "out"
    first = process([("test", receipts)], tmp_path, output)
    record_bytes = next((output / "records").glob("*.jsonl")).read_bytes()
    second = process([("test", receipts)], tmp_path, output)
    assert first["manifest_set_sha256"] == second["manifest_set_sha256"]
    assert second["resumed_articles"] == 1
    assert record_bytes == next((output / "records").glob("*.jsonl")).read_bytes()
    assert sha(record_bytes) == json.loads(next((output / "manifests").glob("*.json")).read_text())["records_file_sha256"]


def test_case_sensitive_molar_length_and_force_units() -> None:
    row = {"source_family_id": "doi:x", "pmcid": None, "sha256": "a" * 64, "receipt_sha256": "b" * 64, "acquisition_pass": "first"}
    molar = record_for_unit(row, 2020, ["Methods"], "section", "m:1", "Fibroblasts were treated with compound at 100 nM.")
    assert molar is not None
    assert molar["observations"][0]["measured_entity"] == "protocol_parameter:concentration"
    assert molar["observations"][0]["value"]["normalized_value"] == "0.0000001"
    length = record_for_unit(row, 2020, ["Results"], "section", "r:1", "AFM showed fibroblast length of 100 nm.")
    assert length is not None
    assert length["observations"][0]["measured_entity"] == "reported_measurement_candidate:length"
    assert record_for_unit(row, 2020, ["Methods"], "section", "m:2", "Fibroblasts were counted, n = 2 n.") is None
    assert record_for_unit(row, 2020, ["Methods"], "section", "m:3", "Fibroblasts were counted in 2 N samples.") is None
    force = record_for_unit(row, 2020, ["Results"], "section", "r:2", "AFM measured a fibroblast force of 2 N.")
    assert force is not None and force["observations"][0]["measured_entity"].endswith(":force")


def test_lowercase_micrometre_before_drug_is_ambiguous() -> None:
    row = {"source_family_id": "doi:x", "pmcid": None, "sha256": "a" * 64, "receipt_sha256": "b" * 64, "acquisition_pass": "first"}
    record = record_for_unit(row, 2020, ["Methods"], "section", "m:1", "Fibroblasts were treated with 30 μm VER-155008.")
    assert record is not None
    observation = record["observations"][0]
    assert observation["unit"]["normalization_status"] == "unknown"
    assert observation["value"]["normalized_value"] is None
    assert observation["extraction"]["candidates"] == ["length", "concentration"]


def test_repeated_text_at_distinct_locators_has_distinct_observation_ids() -> None:
    payload = b"<article><body><sec id='r1'><title>Results</title><p>AFM of fibroblasts reported 2 kPa.</p></sec><sec id='r2'><title>Results</title><p>AFM of fibroblasts reported 2 kPa.</p></sec></body></article>"
    records = extract_article(payload, receipt(payload))
    ids = [record["observations"][0]["observation_id"] for record in records]
    assert len(ids) == 2
    assert len(set(ids)) == 2
