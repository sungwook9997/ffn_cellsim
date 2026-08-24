from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_full_index_is_complete_compact_and_proposed():
    records = [json.loads(line) for line in (RESULTS / "article_screen_index.jsonl").read_text("utf-8").splitlines()]
    assert len(records) == 2701
    assert len({record["source_family_id"] for record in records}) == 2701
    assert records == sorted(records, key=lambda record: record["source_family_id"])
    assert {record["decision"] for record in records} <= {"hold", "manual_priority", "exception_candidate"}
    assert all(record["authority_status"] == "proposed" for record in records)
    assert all("text" not in record and "title" not in record and "authors" not in record for record in records)
    assert "tier_a" not in (RESULTS / "article_screen_index.jsonl").read_text("utf-8").lower()


def test_counts_and_queue_reconcile():
    summary = json.loads((RESULTS / "summary.json").read_text("utf-8"))
    assert summary["articles"] == 2701
    assert sum(summary["decisions"].values()) == 2701
    with (RESULTS / "manual_review_queue.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == summary["decisions"]["manual_priority"] + summary["decisions"]["exception_candidate"]
    assert [int(row["rank"]) for row in rows] == list(range(1, len(rows) + 1))
    scores = [int(row["priority_score"]) for row in rows]
    assert scores == sorted(scores, reverse=True)


def test_recorded_hashes_reproduce():
    sums = json.loads((RESULTS / "SHA256SUMS.json").read_text("utf-8"))
    for name, expected in sums.items():
        assert digest(RESULTS / name) == expected
    summary = json.loads((RESULTS / "summary.json").read_text("utf-8"))
    assert summary["compact_index_sha256"] == digest(RESULTS / "article_screen_index.jsonl")
    assert summary["manual_queue_sha256"] == digest(RESULTS / "manual_review_queue.csv")
    assert summary["leakage_groups_sha256"] == digest(RESULTS / "dataset_leakage_groups.json")


def test_leakage_groups_only_report_cross_source_groups():
    raw = json.loads((RESULTS / "dataset_leakage_groups.json").read_text("utf-8"))
    assert raw["authority_status"] == "proposed"
    assert all(group["source_family_count"] == len(group["source_families"]) > 1 for group in raw["groups"])
    assert len({group["dataset_group_hash"] for group in raw["groups"]}) == len(raw["groups"])
