from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_heartbeat(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("COLLAB_MCP_DB", str(tmp_path / "heartbeat.db"))
    module_path = Path(__file__).parents[1] / "tools" / "collab_mcp" / "heartbeat_daemon.py"
    spec = importlib.util.spec_from_file_location("heartbeat_daemon_under_test", module_path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    monkeypatch.setitem(sys.modules, "heartbeat_daemon_under_test", module)
    spec.loader.exec_module(module)
    return module


def test_build_panes_prefixes_workroom(monkeypatch, tmp_path):
    heartbeat = _load_heartbeat(monkeypatch, tmp_path)

    panes = heartbeat.build_panes("design-discussion")

    assert panes[0] == (
        "design-discussion-claude",
        "design-discussion-claude:0.0",
    )
    assert panes[-1] == (
        "design-discussion-codex",
        "design-discussion-codex:0.0",
    )


def test_tick_seeds_missing_status_rows(monkeypatch, tmp_path):
    heartbeat = _load_heartbeat(monkeypatch, tmp_path)
    posted: list[tuple[str, int, str, str]] = []

    monkeypatch.setattr(heartbeat, "read_status", lambda agent: None)
    monkeypatch.setattr(heartbeat, "pane_alive", lambda target: True)
    monkeypatch.setattr(
        heartbeat,
        "post_status",
        lambda agent, percent, activity, topic: posted.append((agent, percent, activity, topic)),
    )

    heartbeat.tick([("implementation-work-codex", "implementation-work-codex:0.0")])

    assert posted == [("implementation-work-codex", 0, "idle", "")]


def test_tick_preserves_existing_status(monkeypatch, tmp_path):
    heartbeat = _load_heartbeat(monkeypatch, tmp_path)
    posted: list[tuple[str, int, str, str]] = []

    monkeypatch.setattr(
        heartbeat,
        "read_status",
        lambda agent: {"percent": 42, "activity": "reviewing", "topic": "v2"},
    )
    monkeypatch.setattr(heartbeat, "pane_alive", lambda target: True)
    monkeypatch.setattr(
        heartbeat,
        "post_status",
        lambda agent, percent, activity, topic: posted.append((agent, percent, activity, topic)),
    )

    heartbeat.tick([("design-discussion-claude", "design-discussion-claude:0.0")])

    assert posted == [("design-discussion-claude", 42, "reviewing", "v2")]
