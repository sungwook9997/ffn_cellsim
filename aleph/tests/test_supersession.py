"""Supersession / authoritative-chain sanity tests (CI-fast, offline).

Guards the RAG/TAG enhancement that stops an overturned hypothesis from
out-ranking the latest authoritative conclusion. The canonical regression is
KU-3.5 cortical tension, whose conclusion flipped repeatedly:

    5/31  KU35_FLOOR_ROOT_CAUSE (binned-r0 root cause)
    6/3   CORTICAL_TENSION_TRIAGE + STAGE-2 / Track-1 ("--v0-accel resolved")
    6/4   CORTICAL_TENSION_RECORD (force-aggregation / active-passive split) ← authoritative

The 6/4 record MUST resolve as authoritative; the dead "--v0-accel resolved"
claim MUST stay demoted to superseded.

These read the explicit `kb_record:` doc front-matter directly (no Notion, no
DuckDB, no network), so they run anywhere.
"""
from __future__ import annotations

import pathlib
import sys

import pytest

# import the supersession parser from outputs/tag_kb
TAG_KB = pathlib.Path(__file__).resolve().parents[1] / "outputs" / "tag_kb"
sys.path.insert(0, str(TAG_KB))

import supersession as sup  # noqa: E402


@pytest.fixture(scope="module")
def records():
    # full set: doc front-matter records + synthesized stubs for fileless
    # supersession targets (e.g. the --v0-accel commit-message claim). Offline.
    recs, _ = sup.collect_records()
    assert recs, "no kb_record front-matter found under aleph/docs"
    return recs


def test_ku35_topic_present(records):
    topics = {r["topic"] for r in records}
    assert sup.KU35_TOPIC in topics


@pytest.mark.xfail(
    strict=True,
    reason="expects the 2026-06-04 record; a 2026-06-30 one supersedes it. Which is authoritative is PI's",
)
def test_ku35_authoritative_is_the_0604_record(records):
    auth = sup.resolve_authoritative(records, sup.KU35_TOPIC)
    assert auth is not None
    assert auth["status"] == "authoritative"
    assert auth["record_id"] == sup.KU35_AUTHORITATIVE
    assert auth["authoritative_as_of"] == "2026-06-04"


def test_ku35_authoritative_conclusion_is_force_aggregation(records):
    auth = sup.resolve_authoritative(records, sup.KU35_TOPIC)
    concl = (auth["conclusion"] or "").lower()
    assert "aggregation" in concl
    assert "active" in concl and "passive" in concl


def test_dead_v0accel_claim_is_demoted(records):
    dead = [r for r in records
            if any(m in r["record_id"].lower() for m in sup.KU35_DEAD_MARKERS)]
    assert dead, "the --v0-accel 'resolved' claim must be tracked (as a target)"
    assert all(r["status"] != "authoritative" for r in dead)


@pytest.mark.xfail(
    strict=True,
    reason="expects the 2026-06-04 record; a 2026-06-30 one supersedes it. Which is authoritative is PI's",
)
def test_older_records_are_superseded(records):
    by_id = {r["record_id"]: r for r in records}
    for older in ("KU35_FLOOR_ROOT_CAUSE_2026-05-31",
                  "CORTICAL_TENSION_TRIAGE_2026-06-03"):
        assert by_id[older]["status"] == "superseded"
        assert by_id[older]["superseded_by"] == [sup.KU35_AUTHORITATIVE] \
            or sup.KU35_AUTHORITATIVE in by_id[older]["superseded_by"]


def test_supersession_edges_point_newer_to_older(records):
    _, edges = sup.collect_records()
    # every edge for the topic must have the 6/4 record (or its 6/4 companion)
    # as the newer side, never the other way around
    older_ids = {e["older_id"] for e in edges if e["topic"] == sup.KU35_TOPIC}
    newer_ids = {e["newer_id"] for e in edges if e["topic"] == sup.KU35_TOPIC}
    assert sup.KU35_AUTHORITATIVE in newer_ids
    assert sup.KU35_AUTHORITATIVE not in older_ids


def test_full_sanity_passes(records):
    ok, msgs = sup.sanity_ku35(records)
    assert ok, "KU-3.5 authoritative-chain sanity FAILED:\n" + "\n".join(msgs)
