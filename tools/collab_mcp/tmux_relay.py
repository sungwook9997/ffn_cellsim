"""Safe tmux relay for the ACS collab workroom.

The relay watches the shared SQLite message table and injects compact wrapper
prompts into Claude/Codex tmux panes. It never injects raw PI/LLM message
bodies; the target agent must read the exact message through MCP/bootstrap/read
and respond through MCP `send`.
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import pathlib
import re
import sqlite3
import subprocess
import sys
import time
from contextlib import closing
from typing import Iterable

DB_PATH = pathlib.Path(
    os.environ.get(
        "COLLAB_MCP_DB",
        str(pathlib.Path.home() / ".acs-collab" / "inbox.db"),
    )
).expanduser()

DEFAULT_RELAY_NAME = os.environ.get("COLLAB_TMUX_RELAY_NAME", "default")
DEFAULT_POLL_INTERVAL_S = float(os.environ.get("COLLAB_TMUX_POLL_INTERVAL", "0.25"))
DEFAULT_MAX_MESSAGES = int(os.environ.get("COLLAB_TMUX_MAX_MESSAGES", "20"))
SAFE_LABEL_RE = re.compile(r"[^0-9A-Za-z가-힣._/@,+: -]+")

NOTIFY_DESKTOP = os.environ.get("COLLAB_NOTIFY_DESKTOP", "1").lower() not in {
    "0", "false", "no", "off", ""
}
NOTIFY_AUTHORS = {
    a.strip().lower()
    for a in os.environ.get("COLLAB_NOTIFY_AUTHORS", "claude,codex").split(",")
    if a.strip()
}
NOTIFY_BODY_MAX = 200
NOTIFY_TITLE_MAX = 80
NOTIFY_SUBTITLE_MAX = 80


@dataclasses.dataclass(frozen=True)
class WorkroomMessage:
    id: int
    author: str
    addressee: str
    topic: str
    status: str
    body: str = ""
    room: str = ""


@dataclasses.dataclass(frozen=True)
class RelayConfig:
    relay_name: str
    claude_target: str
    codex_target: str
    include_llm_messages: bool
    dry_run: bool
    workroom: str = ""  # empty = legacy mode (no room filter)


def _db() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def _init_schema() -> None:
    with closing(_db()) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS relay_cursors(
                relay_name TEXT PRIMARY KEY,
                last_routed_id INTEGER NOT NULL DEFAULT 0,
                updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
            )
            """
        )


def _load_cursor(conn: sqlite3.Connection, relay_name: str) -> int:
    row = conn.execute(
        "SELECT last_routed_id FROM relay_cursors WHERE relay_name=?",
        (relay_name,),
    ).fetchone()
    return int(row[0]) if row else 0


def _save_cursor(conn: sqlite3.Connection, relay_name: str, message_id: int) -> None:
    conn.execute(
        "INSERT INTO relay_cursors(relay_name, last_routed_id, updated_at) "
        "VALUES(?,?,CURRENT_TIMESTAMP) "
        "ON CONFLICT(relay_name) DO UPDATE SET "
        "last_routed_id=excluded.last_routed_id, updated_at=CURRENT_TIMESTAMP",
        (relay_name, message_id),
    )


def fetch_new_messages(relay_name: str, limit: int) -> list[WorkroomMessage]:
    _init_schema()
    with closing(_db()) as conn:
        last = _load_cursor(conn, relay_name)
        cols = {row[1] for row in conn.execute("PRAGMA table_info(messages)").fetchall()}
        select = "SELECT id, author, addressee, topic, status, body"
        if "room" in cols:
            select += ", room"
        rows = conn.execute(
            f"{select} FROM messages WHERE id>? ORDER BY id ASC LIMIT ?",
            (last, limit),
        ).fetchall()
    out: list[WorkroomMessage] = []
    for row in rows:
        out.append(
            WorkroomMessage(
                id=int(row[0]),
                author=str(row[1]).lower(),
                addressee=str(row[2]).lower(),
                topic=str(row[3]),
                status=str(row[4]),
                body=str(row[5] or ""),
                room=str(row[6]) if len(row) > 6 and row[6] else "",
            )
        )
    return out


def mark_routed(relay_name: str, message_id: int) -> None:
    _init_schema()
    with closing(_db()) as conn:
        conn.execute("BEGIN IMMEDIATE")
        try:
            _save_cursor(conn, relay_name, message_id)
            conn.execute("COMMIT")
        except Exception:
            conn.execute("ROLLBACK")
            raise


def latest_message_id() -> int:
    _init_schema()
    with closing(_db()) as conn:
        row = conn.execute("SELECT COALESCE(MAX(id), 0) FROM messages").fetchone()
    return int(row[0])


def init_cursor_now(relay_name: str) -> int:
    message_id = latest_message_id()
    mark_routed(relay_name, message_id)
    return message_id


def init_cursor_if_missing(relay_name: str) -> int | None:
    """Initialize a relay cursor at the current latest message only once.

    New workroom relays should not replay the entire shared history into
    fresh chat panes. Existing relays must keep their cursor so a restart
    does not skip messages that arrived while they were down.
    """
    _init_schema()
    with closing(_db()) as conn:
        if _load_cursor(conn, relay_name):
            return None
    return init_cursor_now(relay_name)


def wrapper_prompt(message: WorkroomMessage, target_agent: str, workroom: str = "") -> str:
    """Return the only text injected into tmux panes.

    The relay writes into the one Claude/Codex agent session for the
    workroom. The wrapper names the room and asks the target to read the
    exact message through MCP before responding.
    """
    author = safe_label(message.author)
    addressee = safe_label(message.addressee)
    topic = safe_label(message.topic, max_len=80)
    target_agent = safe_label(target_agent, max_len=20)
    workroom_label = safe_label(workroom, max_len=32) if workroom else ""
    room_field = f" room={workroom_label}" if workroom_label and workroom_label != "unknown" else ""
    return (
        f": mcp_msg id={message.id} from={author} to={addressee} "
        f"topic={topic}{room_field} pane=agent "
        f"priority=immediate ack_first "
        f"action=read_acs_collab_mcp_then_send_if_{target_agent}_should_reply"
    )


def safe_label(value: str, max_len: int = 120) -> str:
    cleaned = SAFE_LABEL_RE.sub("x", value.replace("\n", " ").replace("\r", " "))
    cleaned = " ".join(cleaned.split())
    return cleaned[:max_len] or "unknown"


def parse_targets(addressee: str, body_hint: str = "") -> set[str]:
    """Resolve explicit route labels to agent names.

    The relay normally only has DB metadata, not raw bodies. `body_hint` exists
    for tests and future room-level routing labels.
    """
    text = f"{addressee} {body_hint}".lower()
    if "@both" in text or "both" in text or "claude,codex" in text or "codex,claude" in text:
        return {"claude", "codex"}
    targets: set[str] = set()
    if "@claude" in text or "claude" in text:
        targets.add("claude")
    if "@codex" in text or "codex" in text:
        targets.add("codex")
    return targets


def route_targets(message: WorkroomMessage, include_llm_messages: bool) -> set[str]:
    """Return tmux agents that should receive a wrapper for this message."""
    if message.author == "pi":
        return parse_targets(message.addressee) or {"claude", "codex"}
    if not include_llm_messages:
        return set()
    if message.author == "claude":
        return {"codex"} if "codex" in parse_targets(message.addressee) else set()
    if message.author == "codex":
        return {"claude"} if "claude" in parse_targets(message.addressee) else set()
    return set()


def tmux_send_wrapper(target: str, prompt: str, dry_run: bool = False) -> None:
    """Inject a literal one-line prompt into a tmux pane."""
    if dry_run:
        print(f"[dry-run] tmux target={target}: {prompt}")
        return
    subprocess.run(
        ["tmux", "send-keys", "-t", target, "-l", prompt],
        check=True,
    )
    # Codex/Claude TUIs can lag one event loop tick behind long literal input.
    # Give the pane a moment to materialize the text before submitting it.
    time.sleep(0.08)
    subprocess.run(["tmux", "send-keys", "-t", target, "C-m"], check=True)


def message_belongs_to_relay(message: WorkroomMessage, config: RelayConfig) -> bool:
    """True iff this relay should both notify and dispatch the message.

    Workroom isolation: a relay tied to room X only handles messages
    tagged with room X. Legacy messages (``room == ""``) are handled by
    all relays so the pre-room corpus stays reachable. A relay that
    boots without --workroom (``config.workroom == ""``) behaves like
    the pre-multiplex relay and handles everything.
    """
    if not config.workroom:
        return True
    if not message.room:
        return True
    return message.room == config.workroom


def dispatch_message(message: WorkroomMessage, config: RelayConfig) -> set[str]:
    if not message_belongs_to_relay(message, config):
        return set()
    targets = route_targets(message, config.include_llm_messages)
    for agent in sorted(targets):
        tmux_target = config.claude_target if agent == "claude" else config.codex_target
        if not tmux_target:
            raise ValueError(f"missing tmux target for {agent}")
        tmux_send_wrapper(
            tmux_target,
            wrapper_prompt(message, agent, config.workroom),
            config.dry_run,
        )
    return targets


def _applescript_escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


def _truncate(s: str, max_len: int) -> str:
    s = " ".join(s.split())  # collapse newlines/tabs into single spaces
    if len(s) <= max_len:
        return s
    return s[: max_len - 1].rstrip() + "…"


def notify_desktop(message: WorkroomMessage, dry_run: bool = False) -> bool:
    """Raise a macOS desktop notification for a new LLM-authored message.

    Best-effort and non-blocking. Returns True iff a notification was actually
    fired. PI-authored messages are skipped — PI is the human typing, not a
    notification target. Set ``COLLAB_NOTIFY_DESKTOP=0`` (or false/no/off) to
    silence; ``COLLAB_NOTIFY_AUTHORS`` overrides the default
    ``claude,codex`` allowlist.
    """
    if not NOTIFY_DESKTOP:
        return False
    if message.author not in NOTIFY_AUTHORS:
        return False
    if sys.platform != "darwin":
        return False

    title = _applescript_escape(_truncate(f"ACS Collab — {message.author}", NOTIFY_TITLE_MAX))
    subtitle_raw = f"#{message.id} · {message.topic}"
    if message.status and message.status.lower() not in {"fyi", ""}:
        subtitle_raw = f"{subtitle_raw} [{message.status}]"
    subtitle = _applescript_escape(_truncate(subtitle_raw, NOTIFY_SUBTITLE_MAX))
    body = _applescript_escape(_truncate(message.body or "(no body)", NOTIFY_BODY_MAX))
    script = (
        f'display notification "{body}" '
        f'with title "{title}" subtitle "{subtitle}"'
    )
    if dry_run:
        print(f"[dry-run] osascript: {script}")
        return False
    try:
        subprocess.Popen(
            ["osascript", "-e", script],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except OSError as exc:
        print(f"desktop notify failed: {exc}", file=sys.stderr, flush=True)
        return False
    return True


def run_once(config: RelayConfig, limit: int = DEFAULT_MAX_MESSAGES) -> int:
    messages = fetch_new_messages(config.relay_name, limit)
    routed = 0
    for message in messages:
        # Cursor must always advance, even when the message is filtered
        # out by room — otherwise other-room traffic would re-fetch
        # forever. But notification is only fired when the message
        # actually belongs to this relay, so a workroom-specific
        # heartbeat does not flash a desktop banner for a peer room.
        belongs = message_belongs_to_relay(message, config)
        if belongs:
            targets = dispatch_message(message, config)
            routed += len(targets)
            notify_desktop(message, dry_run=config.dry_run)
        if not config.dry_run:
            mark_routed(config.relay_name, message.id)
    return routed


def run_forever(config: RelayConfig, interval_s: float, limit: int) -> None:
    while True:
        try:
            routed = run_once(config, limit=limit)
            if routed:
                print(f"routed {routed} tmux wrapper prompt(s)", flush=True)
        except Exception as exc:
            print(f"relay error: {exc}", file=sys.stderr, flush=True)
        time.sleep(interval_s)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Route ACS workroom messages into tmux panes.")
    parser.add_argument("--relay-name", default=DEFAULT_RELAY_NAME)
    parser.add_argument("--claude-target", default=os.environ.get("COLLAB_TMUX_CLAUDE_TARGET", ""))
    parser.add_argument("--codex-target", default=os.environ.get("COLLAB_TMUX_CODEX_TARGET", ""))
    parser.add_argument(
        "--workroom",
        default=os.environ.get("ACS_WORKROOM", ""),
        help="Workroom name; relay only routes messages tagged with this "
             "room (legacy room=='' messages are still routed). "
             "Empty = legacy mode (route everything).",
    )
    parser.add_argument("--include-llm-messages", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--once", action="store_true")
    parser.add_argument(
        "--init-cursor-now",
        action="store_true",
        help="Set relay cursor to the current latest message id and exit.",
    )
    parser.add_argument(
        "--init-cursor-if-missing",
        action="store_true",
        help="Before polling, set a missing relay cursor to the current "
             "latest message id. Existing cursors are preserved.",
    )
    parser.add_argument("--interval", type=float, default=DEFAULT_POLL_INTERVAL_S)
    parser.add_argument("--limit", type=int, default=DEFAULT_MAX_MESSAGES)
    return parser


def main(argv: Iterable[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    config = RelayConfig(
        relay_name=args.relay_name,
        claude_target=args.claude_target,
        codex_target=args.codex_target,
        include_llm_messages=args.include_llm_messages,
        dry_run=args.dry_run,
        workroom=args.workroom.strip().lower(),
    )
    if args.limit < 1 or args.limit > 200:
        raise ValueError("--limit must be in [1, 200]")
    if args.interval <= 0:
        raise ValueError("--interval must be positive")
    if args.init_cursor_now:
        message_id = init_cursor_now(args.relay_name)
        print(f"initialized relay cursor {args.relay_name!r} at message #{message_id}")
        return
    if args.init_cursor_if_missing:
        message_id = init_cursor_if_missing(args.relay_name)
        if message_id is not None:
            print(f"initialized missing relay cursor {args.relay_name!r} at message #{message_id}")
    if args.once:
        routed = run_once(config, limit=args.limit)
        print(f"routed {routed} tmux wrapper prompt(s)")
        return
    run_forever(config, interval_s=args.interval, limit=args.limit)


if __name__ == "__main__":
    main()
