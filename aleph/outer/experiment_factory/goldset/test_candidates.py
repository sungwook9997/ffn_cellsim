import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import generate_candidates


HERE = Path(__file__).parent


def _records():
    return [json.loads(line) for line in (HERE / "candidates_300.jsonl").read_text(encoding="utf-8").splitlines()]


def test_candidate_set_is_exactly_300_unique_proposed_sources():
    records = _records()
    assert len(records) == 300
    assert len({record["source_family_id"] for record in records}) == 300
    assert {record["authority_status"] for record in records} == {"proposed"}
    assert {record["review_state"] for record in records} == {"gold_candidate"}


def test_all_twelve_domains_have_the_predeclared_broad_quota():
    counts = Counter(record["primary_domain"] for record in _records())
    assert counts == Counter(generate_candidates.DOMAIN_TARGETS)


def test_no_selected_eligible_leakage_component_is_split():
    records = _records()
    selected_counts = Counter(record["leakage_component_id"] for record in records)
    for record in records:
        assert selected_counts[record["leakage_component_id"]] == record["eligible_component_size"]
    summary = json.loads((HERE / "summary.json").read_text(encoding="utf-8"))
    assert summary["leakage"]["split_selected_components"] == []


def test_each_source_has_two_unassigned_annotation_slots():
    with (HERE / "annotation_assignments_600.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 600
    slots = defaultdict(set)
    for row in rows:
        slots[row["source_family_id"]].add(row["annotator_slot"])
        assert row["assignment_status"] == "unassigned"
    assert set(slots) == {record["source_family_id"] for record in _records()}
    assert all(value == {"A", "B"} for value in slots.values())


def test_committed_digests_match_generated_artifacts():
    hashes = json.loads((HERE / "SHA256SUMS.json").read_text(encoding="utf-8"))
    for name, expected in hashes.items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == expected


def test_generator_replays_byte_identically():
    names = ["annotation_assignments_600.csv", "candidates_300.jsonl", "summary.json", "SHA256SUMS.json"]
    before = {name: (HERE / name).read_bytes() for name in names}
    assert generate_candidates.main() == 0
    after = {name: (HERE / name).read_bytes() for name in names}
    assert after == before
