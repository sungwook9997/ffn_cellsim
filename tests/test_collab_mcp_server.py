from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import pytest


class _FakeFastMCP:
    def __init__(self, _name: str) -> None:
        self.settings = types.SimpleNamespace(host=None, port=None)

    def tool(self):
        def decorator(fn):
            return fn

        return decorator

    def run(self, transport: str) -> None:
        self.transport = transport


class _FakeContext:
    def __init__(self, token: str, author: str, room: str | None = None) -> None:
        headers = {
            "Authorization": f"Bearer {token}",
            "X-Collab-Author": author,
        }
        if room is not None:
            headers["X-Collab-Room"] = room
        self.request_context = types.SimpleNamespace(
            request=types.SimpleNamespace(headers=headers)
        )


def _install_fake_mcp(monkeypatch: pytest.MonkeyPatch) -> None:
    mcp_mod = types.ModuleType("mcp")
    server_mod = types.ModuleType("mcp.server")
    fastmcp_mod = types.ModuleType("mcp.server.fastmcp")
    fastmcp_server_mod = types.ModuleType("mcp.server.fastmcp.server")
    fastmcp_mod.FastMCP = _FakeFastMCP
    fastmcp_server_mod.Context = _FakeContext
    monkeypatch.setitem(sys.modules, "mcp", mcp_mod)
    monkeypatch.setitem(sys.modules, "mcp.server", server_mod)
    monkeypatch.setitem(sys.modules, "mcp.server.fastmcp", fastmcp_mod)
    monkeypatch.setitem(sys.modules, "mcp.server.fastmcp.server", fastmcp_server_mod)


def _load_server(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    _install_fake_mcp(monkeypatch)
    monkeypatch.setenv("COLLAB_MCP_TOKEN", "unit-test-token")
    monkeypatch.setenv("COLLAB_MCP_AUTHORS", "claude,codex")
    monkeypatch.setenv("COLLAB_MCP_DB", str(tmp_path / "inbox.db"))
    monkeypatch.delenv("COLLAB_MCP_LEDGER", raising=False)
    module_path = Path(__file__).parents[1] / "tools" / "collab_mcp" / "server.py"
    spec = importlib.util.spec_from_file_location("collab_mcp_server_under_test", module_path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_send_rejects_unknown_addressee(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    ctx = _FakeContext("unit-test-token", "codex")

    with pytest.raises(ValueError, match="to must be one of"):
        server.send(ctx, to="mallory", topic="handshake", body="hello")


def test_send_allows_pi_addressee(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    ctx = _FakeContext("unit-test-token", "codex", room="implementation-work")

    out = server.send(ctx, to="pi", topic="handshake", body="hello PI", status="ack")

    assert out["room"] == "implementation-work"
    row = server._db().execute(
        "SELECT author, addressee, room, body, status FROM messages WHERE id=?",
        (out["id"],),
    ).fetchone()
    assert row == ("codex", "pi", "implementation-work", "hello PI", "ack")


def test_claim_blocks_conflicting_owner(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude = _FakeContext("unit-test-token", "claude")
    codex = _FakeContext("unit-test-token", "codex")

    first = server.claim(claude, topic="layer-1", intent="draft plan")
    second = server.claim(codex, topic="layer-1", intent="review plan")

    assert first["ok"] is True
    assert second == {
        "ok": False,
        "conflict": True,
        "owner": "claude",
        "claimed_at": first["claimed_at"],
        "room": "design-discussion",
    }


def test_announce_artifact_rejects_empty_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    ctx = _FakeContext("unit-test-token", "codex")

    with pytest.raises(ValueError, match="path must be non-empty"):
        server.announce_artifact(ctx, kind="figure", path=" ")


def test_bootstrap_advances_cursor_by_default(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude = _FakeContext("unit-test-token", "claude")
    codex = _FakeContext("unit-test-token", "codex")

    server.send(claude, to="codex", topic="layer-1", body="hello from claude")
    server.send(codex, to="claude", topic="layer-1", body="ack from codex")
    third = server.send(claude, to="codex", topic="layer-2", body="separate topic")

    boot = server.bootstrap(claude, limit=10)
    assert boot["count"] == 3
    authors = [m["from"] for m in boot["messages"]]
    assert "claude" in authors and "codex" in authors
    assert boot["cursor_advanced_to"] == third["id"]

    follow_up_unread = server.read(claude)
    assert follow_up_unread["count"] == 0


def test_bootstrap_preserves_cursor_when_advance_disabled(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    server = _load_server(monkeypatch, tmp_path)
    claude = _FakeContext("unit-test-token", "claude")
    codex = _FakeContext("unit-test-token", "codex")

    server.send(claude, to="codex", topic="layer-1", body="hello from claude")
    server.send(codex, to="claude", topic="layer-1", body="ack from codex")
    server.send(claude, to="codex", topic="layer-2", body="separate topic")

    boot = server.bootstrap(claude, limit=10, advance_cursor=False)
    assert boot["count"] == 3
    assert boot["cursor_advanced_to"] is None

    follow_up_unread = server.read(claude)
    assert follow_up_unread["count"] == 1
    assert follow_up_unread["messages"][0]["from"] == "codex"


def test_bootstrap_topic_and_since_filters(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude = _FakeContext("unit-test-token", "claude")
    codex = _FakeContext("unit-test-token", "codex")

    server.send(claude, to="codex", topic="layer-1", body="m1")
    second = server.send(codex, to="claude", topic="layer-2", body="m2")
    server.send(claude, to="codex", topic="layer-1", body="m3")

    topic_only = server.bootstrap(claude, topic="layer-1", limit=10)
    assert [m["body"] for m in topic_only["messages"]] == ["m1", "m3"]

    since = server.bootstrap(claude, since_id=second["id"], limit=10)
    assert [m["body"] for m in since["messages"]] == ["m3"]


# ---------------------------------------------------------------------------
# Room (workroom) isolation tests — Task 2
# ---------------------------------------------------------------------------


def test_send_uses_room_from_header(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude_impl = _FakeContext("unit-test-token", "claude", room="implementation-work")

    out = server.send(claude_impl, to="codex", topic="t1", body="hi")
    assert out["room"] == "implementation-work"


def test_send_room_param_overrides_header(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude_design = _FakeContext("unit-test-token", "claude", room="design-discussion")

    out = server.send(claude_design, to="codex", topic="t", body="hi",
                       room="implementation-work")
    assert out["room"] == "implementation-work"


def test_send_rejects_star_room(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude = _FakeContext("unit-test-token", "claude", room="design-discussion")
    with pytest.raises(ValueError, match="must be a single room name"):
        server.send(claude, to="codex", topic="t", body="hi", room="*")


def test_read_filters_by_caller_room(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude_design = _FakeContext("unit-test-token", "claude", room="design-discussion")
    claude_impl = _FakeContext("unit-test-token", "claude", room="implementation-work")
    codex_design = _FakeContext("unit-test-token", "codex", room="design-discussion")
    codex_impl = _FakeContext("unit-test-token", "codex", room="implementation-work")

    server.send(codex_design, to="claude", topic="t", body="design-msg")
    server.send(codex_impl, to="claude", topic="t", body="impl-msg")

    in_design = server.read(claude_design)
    assert [m["body"] for m in in_design["messages"]] == ["design-msg"]
    assert in_design["room"] == "design-discussion"

    in_impl = server.read(claude_impl)
    assert [m["body"] for m in in_impl["messages"]] == ["impl-msg"]


def test_read_room_star_unread_rejected(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude = _FakeContext("unit-test-token", "claude", room="design-discussion")
    with pytest.raises(ValueError, match="cursors are per-room"):
        server.read(claude, room="*", unread_only=True)


def test_per_room_cursors_do_not_starve_other_rooms(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    """The bug codex-work flagged: with single-author cursors, advancing in
    room A would skip lower-id unread messages in room B. With per-(author,
    room) cursors, each room has its own cursor and neither starves."""
    server = _load_server(monkeypatch, tmp_path)
    claude_design = _FakeContext("unit-test-token", "claude", room="design-discussion")
    claude_impl = _FakeContext("unit-test-token", "claude", room="implementation-work")
    codex_design = _FakeContext("unit-test-token", "codex", room="design-discussion")
    codex_impl = _FakeContext("unit-test-token", "codex", room="implementation-work")

    # Interleave so room=impl message has a LOWER id than room=design
    impl_msg = server.send(codex_impl, to="claude", topic="t", body="impl-1")
    design_msg = server.send(codex_design, to="claude", topic="t", body="design-1")
    assert impl_msg["id"] < design_msg["id"]

    # Claude reads design first, advancing the design-room cursor
    design_read = server.read(claude_design)
    assert [m["body"] for m in design_read["messages"]] == ["design-1"]

    # The impl-room message is older (smaller id) but must still be unread
    # in implementation-work, because its cursor is independent
    impl_read = server.read(claude_impl)
    assert [m["body"] for m in impl_read["messages"]] == ["impl-1"]


def test_status_reports_per_room_unread(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    codex_design = _FakeContext("unit-test-token", "codex", room="design-discussion")
    codex_impl = _FakeContext("unit-test-token", "codex", room="implementation-work")

    server.send(codex_design, to="claude", topic="t", body="d1")
    server.send(codex_impl, to="claude", topic="t", body="i1")
    server.send(codex_impl, to="claude", topic="t", body="i2")

    claude_design = _FakeContext("unit-test-token", "claude", room="design-discussion")
    s = server.status(claude_design)
    assert s["room"] == "design-discussion"
    assert "design-discussion" in s["rooms"] and "implementation-work" in s["rooms"]
    # claude's per-room unread: 1 in design, 2 in impl
    assert s["unread_by_room"]["claude"]["design-discussion"] == 1
    assert s["unread_by_room"]["claude"]["implementation-work"] == 2
    # The legacy `unread` field should reflect the caller's room
    assert s["unread"]["claude"] == 1


def test_legacy_db_messages_map_to_default_room(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    """Pre-Task-2 DBs have no room column. Migration must add it with
    DEFAULT_ROOM so existing messages remain visible to design-discussion
    callers and don't silently disappear from the timeline."""
    import sqlite3

    db_path = tmp_path / "inbox.db"
    # Build a pre-Task-2 schema by hand and insert a legacy row
    with sqlite3.connect(db_path) as conn:
        conn.executescript(
            """
            CREATE TABLE messages(
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ts TEXT NOT NULL,
                author TEXT NOT NULL,
                addressee TEXT NOT NULL,
                topic TEXT NOT NULL,
                body TEXT NOT NULL,
                status TEXT NOT NULL,
                refs TEXT NOT NULL DEFAULT '[]'
            );
            CREATE TABLE cursors(
                author TEXT PRIMARY KEY,
                last_seen_id INTEGER NOT NULL DEFAULT 0
            );
            INSERT INTO messages(ts, author, addressee, topic, body, status, refs)
                VALUES ('2026-05-03 10:00 KST', 'codex', 'claude', 'old-topic',
                        'pre-room body', 'FYI', '[]');
            """
        )

    server = _load_server(monkeypatch, tmp_path)
    claude_design = _FakeContext("unit-test-token", "claude", room="design-discussion")
    out = server.read(claude_design, unread_only=False)
    bodies = [m["body"] for m in out["messages"]]
    assert "pre-room body" in bodies
    legacy = next(m for m in out["messages"] if m["body"] == "pre-room body")
    assert legacy["room"] == "design-discussion"


def test_invalid_room_param_rejected(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude = _FakeContext("unit-test-token", "claude", room="design-discussion")
    with pytest.raises(ValueError, match="room must match"):
        server.send(claude, to="codex", topic="t", body="hi", room="BAD ROOM")


def test_claim_is_room_scoped(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Same topic claimed in two rooms must coexist; same room collides."""
    server = _load_server(monkeypatch, tmp_path)
    claude_design = _FakeContext("unit-test-token", "claude", room="design-discussion")
    codex_design = _FakeContext("unit-test-token", "codex", room="design-discussion")
    claude_impl = _FakeContext("unit-test-token", "claude", room="implementation-work")

    # Claude in design claims topic-A
    a = server.claim(claude_design, topic="topic-A")
    assert a["ok"] is True and a["room"] == "design-discussion"

    # Codex in design tries the same topic-A → conflict
    b = server.claim(codex_design, topic="topic-A")
    assert b["ok"] is False and b.get("conflict") is True
    assert b["owner"] == "claude" and b["room"] == "design-discussion"

    # Claude in impl claims topic-A → no conflict, parallel claim
    c = server.claim(claude_impl, topic="topic-A")
    assert c["ok"] is True and c["room"] == "implementation-work"


def test_release_is_room_scoped(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    claude_design = _FakeContext("unit-test-token", "claude", room="design-discussion")
    claude_impl = _FakeContext("unit-test-token", "claude", room="implementation-work")

    server.claim(claude_design, topic="t")
    # release in the wrong room must NOT release the design claim
    miss = server.release(claude_impl, topic="t", summary="oops")
    assert miss["ok"] is False and miss["error"] == "no active claim"

    hit = server.release(claude_design, topic="t", summary="done")
    assert hit["ok"] is True and hit["room"] == "design-discussion"


def test_status_active_claims_filter_to_caller_room(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    server = _load_server(monkeypatch, tmp_path)
    claude_design = _FakeContext("unit-test-token", "claude", room="design-discussion")
    claude_impl = _FakeContext("unit-test-token", "claude", room="implementation-work")

    server.claim(claude_design, topic="design-only-topic")
    server.claim(claude_impl, topic="impl-only-topic")

    s_design = server.status(claude_design)
    topics = [c["topic"] for c in s_design["active_claims"]]
    assert "design-only-topic" in topics
    assert "impl-only-topic" not in topics  # filtered out
    # Cross-room view still available
    impl_in_cross = [c["topic"] for c in s_design["claims_by_room"]["implementation-work"]]
    assert "impl-only-topic" in impl_in_cross


def test_invalid_x_collab_room_header_rejected(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    server = _load_server(monkeypatch, tmp_path)
    # Strict lowercase: an uppercase header is a client bug, not a thing
    # to silently normalize. Same for whitespace.
    for bad in ("UPPER", "bad room", "design.discussion"):
        claude_bad = _FakeContext("unit-test-token", "claude", room=bad)
        with pytest.raises(ValueError, match="X-Collab-Room header must match"):
            server.send(claude_bad, to="codex", topic="t", body="hi")
