from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).with_name("audit_run.py")


def _write_jsonl(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))


def test_audit_normalizes_families_and_aggregates_acquisition_passes(tmp_path: Path) -> None:
    _write_jsonl(
        tmp_path / "sources" / "seed_primary_sources.jsonl",
        [{"source_family_id": "DOI:10.1/ABC", "doi": "10.1/ABC", "year": 2020}],
    )
    _write_jsonl(
        tmp_path / "snapshots" / "first.jsonl",
        [{"source_family_id": "doi:10.1/abc", "pmcid": "PMC1", "is_open_access_provider_flag": True}],
    )
    _write_jsonl(
        tmp_path / "acquisition" / "second_pass" / "candidates.jsonl",
        [{"source_family_id": "pmcid:PMC2", "pmcid": "PMC2", "is_open_access_provider_flag": True}],
    )
    _write_jsonl(
        tmp_path / "acquisition" / "first.jsonl",
        [{"source_family_id": "doi:10.1/abc", "pmcid": "PMC1", "status": "success"}],
    )
    _write_jsonl(
        tmp_path / "acquisition" / "second_pass" / "receipts.jsonl",
        [{"source_family_id": "pmcid:PMC2", "pmcid": "PMC2", "status": "failed"}],
    )

    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(tmp_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads(completed.stdout)

    assert report["candidate_unique_source_families"] == 2
    assert report["candidate_cross_file_duplicates"] == 1
    assert report["candidate_with_pmcid"] == 2
    assert report["acquisition_records"] == 2
    assert report["acquisition_unique_pmcids"] == 2
    assert report["acquisition_status"] == {"failed": 1, "success": 1}
    assert report["acquisition_receipts"] == [
        "acquisition/first.jsonl",
        "acquisition/second_pass/receipts.jsonl",
    ]


def test_audit_does_not_treat_candidate_jsonl_as_receipts(tmp_path: Path) -> None:
    _write_jsonl(tmp_path / "sources" / "seed_primary_sources.jsonl", [])
    _write_jsonl(
        tmp_path / "acquisition" / "second_pass" / "candidates.jsonl",
        [{"source_family_id": "pmcid:PMC3", "pmcid": "PMC3"}],
    )

    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(tmp_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads(completed.stdout)

    assert report["candidate_unique_source_families"] == 1
    assert report["acquisition_records"] == 0
    assert report["acquisition_receipts"] == []
