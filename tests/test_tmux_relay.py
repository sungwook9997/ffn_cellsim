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


def _init_messages(db_path: Path, with_room: bool = True) -> None:
    conn = sqlite3.connect(db_path)
    schema = (
        """
        CREATE TABLE messages(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            ts TEXT NOT NULL,
            author TEXT NOT NULL,
            addressee TEXT NOT NULL,
            topic TEXT NOT NULL,
            body TEXT NOT NULL,
            status TEXT NOT NULL,
            refs TEXT NOT NULL DEFAULT '[]',
            room TEXT NOT NULL DEFAULT 'design-discussion'
        );
        """
        if with_room
        else """
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
    conn.executescript(schema)
    conn.close()


def _insert_message(
    db_path: Path,
    author: str,
    addressee: str,
    body: str = "secret body",
    room: str | None = None,
) -> int:
    conn = sqlite3.connect(db_path)
    if room is None:
        cur = conn.execute(
            "INSERT INTO messages(ts, author, addressee, topic, body, status, refs) "
            "VALUES('now',?,?,?,?, 'FYI', '[]')",
            (author, addressee, "routing", body),
        )
    else:
        cur = conn.execute(
            "INSERT INTO messages(ts, author, addressee, topic, body, status, refs, room) "
            "VALUES('now',?,?,?,?, 'FYI', '[]', ?)",
            (author, addressee, "routing", body, room),
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
    assert "priority=immediate" in prompt
    assert "ack_first" in prompt
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


def test_wrapper_prompt_includes_workroom_in_work_pane(monkeypatch, tmp_path):
    """When relay knows its workroom, work_pane must be the prefixed
    physical session (e.g. design-discussion-codex-work) so chat panes
    brief the right sibling and not a non-existent unprefixed session."""
    relay = _load_relay(monkeypatch, tmp_path)
    message = relay.WorkroomMessage(
        7, "pi", "claude,codex", "routing", "FYI", room="design-discussion"
    )
    prompt = relay.wrapper_prompt(message, "codex", workroom="design-discussion")
    assert "work_pane=design-discussion-codex-work" in prompt
    assert "room=design-discussion" in prompt


def test_wrapper_prompt_legacy_unprefixed_when_no_workroom(monkeypatch, tmp_path):
    relay = _load_relay(monkeypatch, tmp_path)
    message = relay.WorkroomMessage(8, "pi", "claude,codex", "routing", "FYI")
    prompt = relay.wrapper_prompt(message, "codex")
    assert "work_pane=codex-work" in prompt
    assert "room=" not in prompt


def test_relay_skips_other_room_dispatch_and_notify(monkeypatch, tmp_path, capsys):
    """A workroom-scoped relay must NOT dispatch OR desktop-notify on a
    message tagged with a different room. The cursor still advances so
    the relay never re-fetches."""
    relay = _load_relay(monkeypatch, tmp_path)
    db_path = tmp_path / "relay.db"
    _init_messages(db_path, with_room=True)
    notify_calls: list[int] = []
    monkeypatch.setattr(
        relay,
        "notify_desktop",
        lambda message, dry_run=False: notify_calls.append(message.id) or False,
    )
    own = _insert_message(db_path, "pi", "claude,codex", room="design-discussion")
    other = _insert_message(db_path, "pi", "claude,codex", room="implementation-work")
    config = relay.RelayConfig(
        relay_name="design-relay",
        claude_target="d-claude",
        codex_target="d-codex",
        include_llm_messages=False,
        dry_run=True,
        workroom="design-discussion",
    )
    routed = relay.run_once(config, limit=20)
    captured = capsys.readouterr().out
    # own room: claude + codex dispatched
    assert routed == 2
    assert "d-claude" in captured and "d-codex" in captured
    # other-room message: only own-room id appears in notify calls
    assert notify_calls == [own]
    assert other not in notify_calls


def test_relay_routes_legacy_room_empty_messages(monkeypatch, tmp_path, capsys):
    """Legacy messages (no room column) must still flow through any
    workroom-scoped relay; otherwise pre-Task-2 corpus becomes invisible."""
    relay = _load_relay(monkeypatch, tmp_path)
    db_path = tmp_path / "relay.db"
    _init_messages(db_path, with_room=False)
    monkeypatch.setattr(relay, "notify_desktop", lambda message, dry_run=False: False)
    _insert_message(db_path, "pi", "claude,codex")
    config = relay.RelayConfig(
        relay_name="impl-relay",
        claude_target="i-claude",
        codex_target="i-codex",
        include_llm_messages=False,
        dry_run=True,
        workroom="implementation-work",
    )
    routed = relay.run_once(config, limit=20)
    assert routed == 2  # legacy room='' is routed even by an impl-room relay


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


def test_init_cursor_if_missing_preserves_existing(monkeypatch, tmp_path):
    relay = _load_relay(monkeypatch, tmp_path)
    db_path = tmp_path / "relay.db"
    _init_messages(db_path)
    first = _insert_message(db_path, "pi", "claude,codex")
    _insert_message(db_path, "pi", "claude")
    relay.mark_routed("unit", first)

    initialized = relay.init_cursor_if_missing("unit")

    assert initialized is None
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(
        "SELECT last_routed_id FROM relay_cursors WHERE relay_name='unit'"
    ).fetchone()[0]
    conn.close()
    assert cursor == first


def test_init_cursor_if_missing_starts_new_relay_at_latest(monkeypatch, tmp_path):
    relay = _load_relay(monkeypatch, tmp_path)
    db_path = tmp_path / "relay.db"
    _init_messages(db_path)
    _insert_message(db_path, "pi", "claude,codex")
    latest = _insert_message(db_path, "pi", "claude")

    initialized = relay.init_cursor_if_missing("new-workroom")

    assert initialized == latest
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(
        "SELECT last_routed_id FROM relay_cursors WHERE relay_name='new-workroom'"
    ).fetchone()[0]
    conn.close()
    assert cursor == latest
