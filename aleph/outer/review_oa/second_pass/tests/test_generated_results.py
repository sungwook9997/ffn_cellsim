from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / "results"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_second_pass_index_is_complete_compact_and_proposed():
    text = (RESULTS / "article_screen_index.jsonl").read_text("utf-8")
    records = [json.loads(line) for line in text.splitlines()]
    assert len(records) == 2974
    assert len({item["source_family_id"] for item in records}) == 2974
    assert records == sorted(records, key=lambda item: item["source_family_id"])
    assert all(item["authority_status"] == "proposed" for item in records)
    assert {item["decision"] for item in records} <= {"hold", "manual_priority", "exception_candidate"}
    assert all("text" not in item and "title" not in item and "authors" not in item for item in records)
    assert "tier_a" not in text.lower()


def test_summary_and_queue_reconcile_exactly():
    summary = json.loads((RESULTS / "summary.json").read_text("utf-8"))
    assert summary["articles"] == 2974
    assert sum(summary["decisions"].values()) == 2974
    assert summary["automatic_tier_a_promotions"] == 0
    assert summary["combined"]["articles"] == 5675
    assert summary["combined"]["source_family_overlap"] == 0
    assert sum(summary["combined"]["decisions"].values()) == 5675
    assert summary["combined"]["automatic_tier_a_promotions"] == 0
    with (RESULTS / "manual_review_queue.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == summary["decisions"].get("manual_priority", 0) + summary["decisions"].get("exception_candidate", 0)
    assert [int(item["rank"]) for item in rows] == list(range(1, len(rows) + 1))
    assert [int(item["priority_score"]) for item in rows] == sorted((int(item["priority_score"]) for item in rows), reverse=True)


def test_hashes_reproduce():
    sums = json.loads((RESULTS / "SHA256SUMS.json").read_text("utf-8"))
    for name, expected in sums.items():
        assert digest(RESULTS / name) == expected
    summary = json.loads((RESULTS / "summary.json").read_text("utf-8"))
    assert summary["compact_index_sha256"] == digest(RESULTS / "article_screen_index.jsonl")
    assert summary["manual_queue_sha256"] == digest(RESULTS / "manual_review_queue.csv")
    assert summary["within_leakage_sha256"] == digest(RESULTS / "dataset_leakage_groups.json")
    assert summary["combined_leakage_sha256"] == digest(RESULTS / "combined_dataset_leakage_groups.json")


def test_combined_leakage_groups_are_cross_source_and_counted():
    raw = json.loads((RESULTS / "combined_dataset_leakage_groups.json").read_text("utf-8"))
    summary = json.loads((RESULTS / "summary.json").read_text("utf-8"))
    groups = raw["groups"]
    assert len(groups) == summary["combined"]["dataset_leakage"]["duplicate_groups"]
    assert all(item["source_family_count"] == len(item["members"]) > 1 for item in groups)
    assert len({item["dataset_group_hash"] for item in groups}) == len(groups)
    cross = sum({member["pass"] for member in item["members"]} == {"first", "second"} for item in groups)
    assert cross == summary["combined"]["dataset_leakage"]["cross_pass_groups"]
