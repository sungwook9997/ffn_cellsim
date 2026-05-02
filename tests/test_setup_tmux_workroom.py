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

    def fake_tmux(*args, check=True):
        calls.append(args)
        if args[:2] == ("has-session", "-t"):
            return subprocess.CompletedProcess(args, 1, "", "")
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(setup, "tmux", fake_tmux)

    created = setup.ensure_session("heartbeat", tmp_path, setup.HEARTBEAT_CMD)

    assert created is True
    assert (
        "new-session", "-d", "-s", "heartbeat", "-c", str(tmp_path), setup.HEARTBEAT_CMD
    ) in calls


def test_send_line_clears_partial_prompt_before_enter(monkeypatch):
    setup = _load_setup(monkeypatch)
    calls = []

    def fake_tmux(*args, check=True):
        calls.append(args)
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(setup, "tmux", fake_tmux)

    setup.send_line("codex-chat:0.0", setup.CODEX_CMD)

    assert calls == [
        ("send-keys", "-t", "codex-chat:0.0", "C-c"),
        ("send-keys", "-t", "codex-chat:0.0", setup.CODEX_CMD, "Enter"),
    ]


def test_setup_starts_new_interactive_sessions_inside_shell(monkeypatch, tmp_path):
    setup = _load_setup(monkeypatch)
    calls = []

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
        start_chat=True,
        start_work=False,
        start_win_ssh=False,
        start_heartbeat=True,
    )

    assert "started Codex in codex-chat" in actions
    assert "started heartbeat daemon" in actions
    assert (
        "new-session", "-d", "-s", "codex-chat", "-c", str(tmp_path)
    ) in calls
    assert (
        "new-session", "-d", "-s", "heartbeat", "-c", str(tmp_path), setup.HEARTBEAT_CMD
    ) in calls
    assert ("send-keys", "-t", "codex-chat:0.0", "C-c") in calls
    assert ("send-keys", "-t", "codex-chat:0.0", setup.CODEX_CMD, "Enter") in calls
