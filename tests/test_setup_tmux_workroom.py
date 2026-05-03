from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


def _load_setup(monkeypatch):
    module_path = Path(__file__).parents[1] / "tools" / "collab_mcp" / "setup_tmux_workroom.py"
    spec = importlib.util.spec_from_file_location("setup_tmux_workroom_under_test", module_path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    monkeypatch.setitem(sys.modules, "setup_tmux_workroom_under_test", module)
    spec.loader.exec_module(module)
    return module


def test_new_daemon_session_can_start_command_directly(monkeypatch, tmp_path):
    setup = _load_setup(monkeypatch)
    calls = []
    heartbeat_cmd = setup.HEARTBEAT_CMD_TEMPLATE.format(workroom="design-discussion")

    def fake_tmux(*args, check=True):
        calls.append(args)
        if args[:2] == ("has-session", "-t"):
            return subprocess.CompletedProcess(args, 1, "", "")
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(setup, "tmux", fake_tmux)

    created = setup.ensure_session("design-discussion-heartbeat", tmp_path, heartbeat_cmd)

    assert created is True
    assert (
        "new-session", "-d", "-s", "design-discussion-heartbeat", "-c", str(tmp_path), heartbeat_cmd
    ) in calls


def test_ensure_session_tolerates_duplicate_session_race(monkeypatch, tmp_path):
    setup = _load_setup(monkeypatch)
    calls = []
    has_session_calls = 0

    def fake_tmux(*args, check=True):
        nonlocal has_session_calls
        calls.append(args)
        if args[:2] == ("has-session", "-t"):
            has_session_calls += 1
            return subprocess.CompletedProcess(args, 0 if has_session_calls > 1 else 1, "", "")
        if args[:1] == ("new-session",):
            return subprocess.CompletedProcess(args, 1, "", "duplicate session: design-discussion-heartbeat")
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(setup, "tmux", fake_tmux)

    created = setup.ensure_session("design-discussion-heartbeat", tmp_path, "sleep 60")

    assert created is False
    assert ("new-session", "-d", "-s", "design-discussion-heartbeat", "-c", str(tmp_path), "sleep 60") in calls


def test_send_line_clears_partial_prompt_before_enter(monkeypatch):
    setup = _load_setup(monkeypatch)
    calls = []
    codex_cmd = setup.codex_command("design-discussion")

    def fake_tmux(*args, check=True):
        calls.append(args)
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(setup, "tmux", fake_tmux)

    setup.send_line("design-discussion-codex:0.0", codex_cmd)

    assert calls == [
        ("send-keys", "-t", "design-discussion-codex:0.0", "C-c"),
        ("send-keys", "-t", "design-discussion-codex:0.0", codex_cmd, "Enter"),
    ]


def test_setup_starts_new_interactive_sessions_inside_shell(monkeypatch, tmp_path):
    monkeypatch.setenv("COLLAB_MCP_TOKEN", "unit-test-token")
    setup = _load_setup(monkeypatch)
    setup.PER_ROOM_LOG_DIR = tmp_path / "logs"
    calls = []
    heartbeat_cmd = setup.HEARTBEAT_CMD_TEMPLATE.format(workroom="design-discussion")
    codex_cmd = setup.codex_command("design-discussion")

    def fake_tmux(*args, check=True):
        calls.append(args)
        if args[:2] == ("has-session", "-t"):
            return subprocess.CompletedProcess(args, 1, "", "")
        if args[:1] == ("list-panes",):
            return subprocess.CompletedProcess(args, 0, "zsh\n", "")
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(setup, "tmux", fake_tmux)

    actions = setup.setup(
        tmp_path,
        workroom="design-discussion",
        start_chat=True,
        start_work=False,
        start_win_ssh=False,
        start_heartbeat=True,
    )

    assert "started Codex in design-discussion-codex" in actions
    assert "started heartbeat daemon (design-discussion)" in actions
    assert (
        "new-session", "-d", "-s", "design-discussion-codex", "-c", str(tmp_path)
    ) in calls
    assert (
        "new-session", "-d", "-s", "design-discussion-heartbeat", "-c", str(tmp_path), heartbeat_cmd
    ) in calls
    assert ("send-keys", "-t", "design-discussion-codex:0.0", "C-c") in calls
    assert ("send-keys", "-t", "design-discussion-codex:0.0", codex_cmd, "Enter") in calls
