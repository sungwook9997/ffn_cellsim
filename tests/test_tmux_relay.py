from __future__ import annotations

import importlib.util
import sqlite3
import sys
from pathlib import Path


def _load_relay(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("COLLAB_MCP_DB", str(tmp_path / "relay.db"))
    module_path = Path(__file__).parents[1] / "tools" / "collab_mcp" / "tmux_relay.py"
    spec = importlib.util.spec_from_file_location("tmux_relay_under_test", module_path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    monkeypatch.setitem(sys.modules, "tmux_relay_under_test", module)
    spec.loader.exec_module(module)
    return module


def _init_messages(db_path: Path) -> None:
    conn = sqlite3.connect(db_path)
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
        """
    )
    conn.close()


def _insert_message(db_path: Path, author: str, addressee: str, body: str = "secret body") -> int:
    conn = sqlite3.connect(db_path)
    cur = conn.execute(
        "INSERT INTO messages(ts, author, addressee, topic, body, status, refs) "
        "VALUES('now',?,?,?,?, 'FYI', '[]')",
        (author, addressee, "routing", body),
    )
    conn.commit()
    conn.close()
    return int(cur.lastrowid)


def test_pi_message_routes_to_both_by_default(monkeypatch, tmp_path):
    relay = _load_relay(monkeypatch, tmp_path)
    message = relay.WorkroomMessage(1, "pi", "claude,codex", "routing", "FYI")

    assert relay.route_targets(message, include_llm_messages=False) == {"claude", "codex"}


def test_llm_message_is_not_routed_without_opt_in(monkeypatch, tmp_path):
    relay = _load_relay(monkeypatch, tmp_path)
    message = relay.WorkroomMessage(2, "claude", "codex", "routing", "FYI")

    assert relay.route_targets(message, include_llm_messages=False) == set()
    assert relay.route_targets(message, include_llm_messages=True) == {"codex"}


def test_wrapper_prompt_does_not_include_raw_body(monkeypatch, tmp_path):
    relay = _load_relay(monkeypatch, tmp_path)
    message = relay.WorkroomMessage(3, "pi", "claude,codex", "routing", "FYI")

    prompt = relay.wrapper_prompt(message, "codex")

    assert "id=3" in prompt
    assert "read" in prompt.lower()
    assert "secret body" not in prompt


def test_wrapper_prompt_sanitizes_shell_metacharacters(monkeypatch, tmp_path):
    relay = _load_relay(monkeypatch, tmp_path)
    message = relay.WorkroomMessage(
        4,
        "pi",
        "claude,codex",
        "routing; rm -rf / `oops` $(bad)",
        "FYI",
    )

    prompt = relay.wrapper_prompt(message, "codex")

    assert ";" not in prompt
    assert "`" not in prompt
    assert "$" not in prompt
    assert "(" not in prompt
    assert ")" not in prompt


def test_run_once_dry_run_does_not_mark_cursor(monkeypatch, tmp_path, capsys):
    relay = _load_relay(monkeypatch, tmp_path)
    db_path = tmp_path / "relay.db"
    _init_messages(db_path)
    first = _insert_message(db_path, "pi", "claude,codex")
    _insert_message(db_path, "claude", "codex")
    config = relay.RelayConfig(
        relay_name="unit",
        claude_target="claude-pane",
        codex_target="codex-pane",
        include_llm_messages=False,
        dry_run=True,
    )

    routed = relay.run_once(config, limit=20)
    second_routed = relay.run_once(config, limit=20)

    captured = capsys.readouterr().out
    assert routed == 2
    assert second_routed == 2
    assert "claude-pane" in captured
    assert "codex-pane" in captured
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(
        "SELECT last_routed_id FROM relay_cursors WHERE relay_name='unit'"
    ).fetchone()
    conn.close()
    assert cursor is None


def test_init_cursor_now_marks_latest(monkeypatch, tmp_path):
    relay = _load_relay(monkeypatch, tmp_path)
    db_path = tmp_path / "relay.db"
    _init_messages(db_path)
    _insert_message(db_path, "pi", "claude,codex")
    latest = _insert_message(db_path, "pi", "claude")

    initialized = relay.init_cursor_now("unit")

    assert initialized == latest
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(
        "SELECT last_routed_id FROM relay_cursors WHERE relay_name='unit'"
    ).fetchone()[0]
    conn.close()
    assert cursor == latest
