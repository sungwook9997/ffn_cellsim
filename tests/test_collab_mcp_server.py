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
    def __init__(self, token: str, author: str) -> None:
        self.request_context = types.SimpleNamespace(
            request=types.SimpleNamespace(
                headers={
                    "Authorization": f"Bearer {token}",
                    "X-Collab-Author": author,
                }
            )
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
