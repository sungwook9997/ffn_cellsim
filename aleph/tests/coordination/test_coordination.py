"""Gates for the GPU lease, the run index, and the ownership check.

All three are filesystem coordination primitives with no device and no physics, so every gate here
runs against `tmp_path`. That is the point of keeping them free of Warp: the thing that decides who
may touch the GPU must be testable without one.
"""

from __future__ import annotations

import json
import time

import pytest

from aleph.coordination.gpu_lease import (
    Lease,
    LeaseHeld,
    RunIndex,
    RunRecord,
    acquire,
    read_lease,
    release,
)
from aleph.coordination.ownership import (
    OwnershipViolation,
    check_paths,
    load_ownership,
)


# ── lease ────────────────────────────────────────────────────────────────────────────────────────
def test_acquire_then_read_round_trips(tmp_path) -> None:
    lease = acquire("session-a", minutes=5, reason="native gate", path=tmp_path / "gpu.lease")
    stored = read_lease(tmp_path / "gpu.lease")
    assert stored is not None
    assert stored.holder == "session-a"
    assert stored.reason == "native gate"
    assert stored.is_live()
    assert 0 < stored.seconds_remaining() <= 300
    assert lease.holder == stored.holder


def test_a_live_lease_blocks_everyone_including_its_own_holder(tmp_path) -> None:
    p = tmp_path / "gpu.lease"
    acquire("session-a", minutes=5, path=p)
    with pytest.raises(LeaseHeld) as first:
        acquire("session-b", minutes=5, path=p)
    assert "session-a" in str(first.value)
    # Re-acquiring your OWN live lease is two runs in one session, which is the bug this prevents —
    # not an extension.
    with pytest.raises(LeaseHeld):
        acquire("session-a", minutes=5, path=p)


def test_an_expired_lease_is_takeable_so_a_dead_session_cannot_wedge_the_device(tmp_path) -> None:
    p = tmp_path / "gpu.lease"
    dead = Lease(holder="crashed", host="h", pid=1, acquired_at=time.time() - 7200,
                 expires_at=time.time() - 3600, reason="died mid-run")
    p.write_text(json.dumps(dead.__getstate__() if hasattr(dead, "__getstate__") else {
        "holder": dead.holder, "host": dead.host, "pid": dead.pid,
        "acquired_at": dead.acquired_at, "expires_at": dead.expires_at, "reason": dead.reason}))
    taken = acquire("session-b", minutes=5, path=p)
    assert taken.holder == "session-b"
    assert read_lease(p).holder == "session-b"


def test_a_corrupt_lease_reads_as_absent_not_as_held(tmp_path) -> None:
    # A coordination primitive that fails closed on its own corruption turns one bad write into an
    # outage only a human can clear. The other direction is bounded: takeover is recorded.
    p = tmp_path / "gpu.lease"
    p.write_text("{ this is not json")
    assert read_lease(p) is None
    assert acquire("session-b", minutes=1, path=p).holder == "session-b"


def test_release_refuses_someone_elses_live_lease_but_force_clears_it(tmp_path) -> None:
    p = tmp_path / "gpu.lease"
    acquire("session-a", minutes=5, path=p)
    with pytest.raises(LeaseHeld):
        release("session-b", path=p)
    assert release("session-b", path=p, force=True) is True
    assert read_lease(p) is None


def test_release_of_nothing_is_not_an_error(tmp_path) -> None:
    assert release("session-a", path=tmp_path / "gpu.lease") is False


def test_lease_rejects_a_nameless_holder_or_a_nonpositive_duration(tmp_path) -> None:
    p = tmp_path / "gpu.lease"
    with pytest.raises(ValueError, match="holder"):
        acquire("   ", minutes=5, path=p)
    with pytest.raises(ValueError, match="positive-finite"):
        acquire("session-a", minutes=0, path=p)
    with pytest.raises(ValueError, match="positive-finite"):
        acquire("session-a", minutes=-1, path=p)


# ── run index ────────────────────────────────────────────────────────────────────────────────────
def _record(**over) -> RunRecord:
    base = dict(run_label="drv.py", build_commit="abc123", config_sha256="sha256:dead",
                artifact="out/x.json", device="gbook", started_at=time.time(),
                holder="session-a", population="904 nodes")
    base.update(over)
    return RunRecord(**base)


def test_index_matches_on_build_AND_config_together(tmp_path) -> None:
    index = RunIndex(tmp_path / "runs.jsonl")
    index.append(_record())
    assert len(index.find(build_commit="abc123", config_sha256="sha256:dead")) == 1
    # same config, different build -> a different run
    assert index.find(build_commit="other", config_sha256="sha256:dead") == []
    # same build, different config -> obviously a different run
    assert index.find(build_commit="abc123", config_sha256="sha256:beef") == []


def test_index_survives_a_half_written_final_line(tmp_path) -> None:
    # The normal consequence of a process dying mid-append. It must not make the index unreadable.
    p = tmp_path / "runs.jsonl"
    index = RunIndex(p)
    index.append(_record())
    with p.open("a") as h:
        h.write('{"run_label": "truncated"')
    assert len(list(index)) == 1
    assert index.summary()["n_runs"] == 1


def test_missing_index_is_empty_not_an_error(tmp_path) -> None:
    index = RunIndex(tmp_path / "absent.jsonl")
    assert list(index) == []
    assert index.summary()["n_runs"] == 0


def test_index_summary_counts_distinct_configs_and_builds(tmp_path) -> None:
    index = RunIndex(tmp_path / "runs.jsonl")
    index.append(_record())
    index.append(_record(config_sha256="sha256:beef"))
    index.append(_record(build_commit="def456", holder="session-b"))
    s = index.summary()
    assert s["n_runs"] == 3
    assert s["n_distinct_configs"] == 2
    assert s["n_builds"] == 2
    assert s["holders"] == ["session-a", "session-b"]


# ── ownership ────────────────────────────────────────────────────────────────────────────────────
DECL = {
    "sessions": {"cortex": ["aleph/components/incumbent/**"], "sf-motor": ["aleph/components/motor/**"]},
    "shared": {"STATE.md": "aleph/docs/state_rows.yaml"},
    "free": ["aleph/outputs/**"],
}


def test_a_session_may_commit_its_own_paths() -> None:
    own = load_ownership(data=DECL)
    assert check_paths("cortex", ["aleph/components/incumbent/driver.py"], own) == []


def test_a_session_may_not_commit_another_lane() -> None:
    own = load_ownership(data=DECL)
    [complaint] = check_paths("cortex", ["aleph/components/motor/segment_motor.py"], own)
    assert "sf-motor" in complaint


def test_a_shared_file_names_the_append_only_file_to_edit_instead() -> None:
    own = load_ownership(data=DECL)
    [complaint] = check_paths("cortex", ["STATE.md"], own)
    assert "SHARED" in complaint
    assert "state_rows.yaml" in complaint


def test_free_space_is_writable_by_anyone() -> None:
    own = load_ownership(data=DECL)
    assert check_paths("cortex", ["aleph/outputs/ac/observe/x.json"], own) == []


def test_an_unclaimed_path_is_refused_and_says_how_to_claim_it() -> None:
    own = load_ownership(data=DECL)
    [complaint] = check_paths("cortex", ["aleph/laws/relax.py"], own)
    assert "unclaimed" in complaint


def test_an_undeclared_session_owns_nothing() -> None:
    # A typo in the session name must not grant silent write access to the whole tree.
    own = load_ownership(data=DECL)
    [complaint] = check_paths("cortexx", ["aleph/components/incumbent/driver.py"], own)
    assert "not declared" in complaint


def test_every_offending_path_is_reported_in_one_pass() -> None:
    own = load_ownership(data=DECL)
    complaints = check_paths(
        "cortex", ["aleph/components/incumbent/ok.py", "aleph/components/motor/a.py", "STATE.md", "aleph/laws/b.py"],
        own)
    assert len(complaints) == 3


def test_overlapping_claims_are_rejected_at_load_not_at_commit() -> None:
    with pytest.raises(OwnershipViolation, match="claimed by both"):
        load_ownership(data={"sessions": {"a": ["x/**"], "b": ["x/**"]}})


def test_the_checked_in_declaration_is_valid_and_disjoint() -> None:
    # The real file, not a fixture: an invalid claim set would disable enforcement silently.
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[3]
    own = load_ownership(root / "ownership.yaml")
    assert own.sessions, "no lanes declared"
    assert own.is_free("aleph/outputs/ac/observe/x.json")
    # STATE.md is FREE, not shared, and the distinction is deliberate: its tier-(a) block is
    # generated and `render_state_rows.py --check` in the hook protects that mechanically, while the
    # prose cannot be protected without serialising the project on one session. Marking the whole
    # file shared refused every edit including prose-only ones — it fired on the very commit that
    # introduced it.
    assert own.is_free("STATE.md")
    assert own.shared_redirect("STATE.md") is None
    assert own.shared_redirect("aleph/outputs/tag_kb/results_manifest.yaml")
    for lane, globs in own.sessions.items():
        assert globs, f"lane {lane} claims nothing, which would own nothing"
