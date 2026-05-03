"""Browser room for the Claude/Codex/PI collaboration channel.

This is a tiny stdlib-only web UI over the same SQLite database used by
``server.py``. It intentionally exposes no command execution and no file reads:
the browser can inspect recent messages, active claims, artifact pointers, and
post PI messages into the shared queue.
"""

from __future__ import annotations

import datetime as dt
import email.parser
import email.policy
import html
import http.cookies
import hmac
import json
import os
import pathlib
import re
import sqlite3
import sys
import urllib.parse
import uuid
from contextlib import closing
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any

ALLOWED_STATUSES = {
    "FYI",
    "ack",
    "proposal",
    "review",
    "decision-needed",
    "open-question",
    "blocker",
}
ADDRESSEE_PRESETS = ("claude,codex", "claude", "codex")
APPROVAL_STATUSES = {"decision-needed", "blocker"}
APPROVAL_ACTIONS = {
    "approve": ("ack", "Approved"),
    "reject": ("blocker", "Rejected"),
    "changes": ("open-question", "Needs changes"),
}
DB_PATH = pathlib.Path(
    os.environ.get(
        "COLLAB_MCP_DB",
        str(pathlib.Path.home() / ".acs-collab" / "inbox.db"),
    )
).expanduser()
LEDGER_PATH = (
    pathlib.Path(os.environ["COLLAB_MCP_LEDGER"]).expanduser()
    if os.environ.get("COLLAB_MCP_LEDGER")
    else None
)
ROOM_TOKEN = os.environ.get("COLLAB_ROOM_TOKEN", os.environ.get("COLLAB_MCP_TOKEN", ""))
# Programmatic agents (chat panes, work panes, daemons) may already hold
# COLLAB_MCP_TOKEN as a bearer; accepting both removes a footgun where the
# agent has the MCP secret but not the room cookie token.
ROOM_TOKENS_ACCEPTED = tuple(
    t for t in (
        ROOM_TOKEN,
        os.environ.get("COLLAB_MCP_TOKEN", ""),
        os.environ.get("COLLAB_ROOM_TOKEN", ""),
    )
    if t
)
PI_AUTHOR = os.environ.get("COLLAB_ROOM_AUTHOR", "pi").strip().lower() or "pi"
UPLOAD_DIR = pathlib.Path(os.environ.get("COLLAB_ROOM_UPLOAD_DIR", "/tmp/acs-collab/uploads")).expanduser()
MAX_FORM_BYTES = 65536
MAX_UPLOAD_BYTES = int(os.environ.get("COLLAB_ROOM_MAX_UPLOAD_BYTES", str(50 * 1024 * 1024)))
ROOM_NAME_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,31}$")
DEFAULT_ROOM = os.environ.get("COLLAB_DEFAULT_ROOM", "design-discussion").strip().lower() or "design-discussion"
if not ROOM_NAME_RE.match(DEFAULT_ROOM):
    DEFAULT_ROOM = "design-discussion"
ROOM_PRESETS = tuple(
    room
    for room in (
        part.strip()
        for part in os.environ.get("COLLAB_ROOM_PRESETS", "design-discussion,implementation-work").split(",")
    )
    if room and ROOM_NAME_RE.match(room)
)

DB_PATH.parent.mkdir(parents=True, exist_ok=True)
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)


def _db() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def _init_schema() -> None:
    with closing(_db()) as conn:
        conn.executescript(
            f"""
            CREATE TABLE IF NOT EXISTS messages(
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ts TEXT NOT NULL,
                author TEXT NOT NULL,
                addressee TEXT NOT NULL,
                topic TEXT NOT NULL,
                body TEXT NOT NULL,
                status TEXT NOT NULL,
                refs TEXT NOT NULL DEFAULT '[]',
                room TEXT NOT NULL DEFAULT '{DEFAULT_ROOM}'
            );
            CREATE INDEX IF NOT EXISTS idx_messages_id ON messages(id);

            CREATE TABLE IF NOT EXISTS claims(
                room TEXT NOT NULL DEFAULT '{DEFAULT_ROOM}',
                topic TEXT NOT NULL,
                author TEXT NOT NULL,
                claimed_at TEXT NOT NULL,
                released_at TEXT,
                summary TEXT,
                PRIMARY KEY(room, topic)
            );

            CREATE TABLE IF NOT EXISTS cursors(
                author TEXT NOT NULL,
                room TEXT NOT NULL DEFAULT '{DEFAULT_ROOM}',
                last_seen_id INTEGER NOT NULL DEFAULT 0,
                PRIMARY KEY(author, room)
            );

            CREATE TABLE IF NOT EXISTS artifacts(
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ts TEXT NOT NULL,
                author TEXT NOT NULL,
                kind TEXT NOT NULL,
                path TEXT NOT NULL,
                note TEXT NOT NULL DEFAULT ''
            );

            CREATE TABLE IF NOT EXISTS agent_status(
                agent TEXT PRIMARY KEY,
                percent INTEGER NOT NULL DEFAULT 0,
                activity TEXT NOT NULL DEFAULT '',
                topic TEXT NOT NULL DEFAULT '',
                updated_at TEXT NOT NULL,
                last_active_at TEXT
            );

            CREATE TABLE IF NOT EXISTS approval_settings(
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            """
        )
        conn.execute(
            "INSERT OR IGNORE INTO approval_settings(key, value) "
            "SELECT 'enabled_from_id', CAST(COALESCE(MAX(id), 0) AS TEXT) FROM messages"
        )
        # Idempotent migration for databases that pre-date the
        # last_active_at column. Existing rows are seeded from updated_at
        # so heartbeat-only refreshes don't regress live agents into the
        # "asleep" state on first room.py boot after upgrade.
        cols = {row[1] for row in conn.execute("PRAGMA table_info(agent_status)").fetchall()}
        if "last_active_at" not in cols:
            conn.execute("ALTER TABLE agent_status ADD COLUMN last_active_at TEXT")
            conn.execute(
                "UPDATE agent_status SET last_active_at = updated_at "
                "WHERE last_active_at IS NULL"
            )
        message_cols = {row[1] for row in conn.execute("PRAGMA table_info(messages)").fetchall()}
        if "room" not in message_cols:
            conn.execute(
                f"ALTER TABLE messages ADD COLUMN room TEXT NOT NULL DEFAULT '{DEFAULT_ROOM}'"
            )
        conn.execute("CREATE INDEX IF NOT EXISTS idx_messages_room_id ON messages(room, id)")

        claim_cols = {row[1] for row in conn.execute("PRAGMA table_info(claims)").fetchall()}
        if "room" not in claim_cols:
            conn.executescript(
                f"""
                CREATE TABLE claims_new(
                    room TEXT NOT NULL DEFAULT '{DEFAULT_ROOM}',
                    topic TEXT NOT NULL,
                    author TEXT NOT NULL,
                    claimed_at TEXT NOT NULL,
                    released_at TEXT,
                    summary TEXT,
                    PRIMARY KEY(room, topic)
                );
                INSERT OR REPLACE INTO claims_new(room, topic, author, claimed_at, released_at, summary)
                    SELECT '{DEFAULT_ROOM}', topic, author, claimed_at, released_at, summary FROM claims;
                DROP TABLE claims;
                ALTER TABLE claims_new RENAME TO claims;
                """
            )

        cursor_cols = {row[1] for row in conn.execute("PRAGMA table_info(cursors)").fetchall()}
        if "room" not in cursor_cols:
            conn.executescript(
                f"""
                CREATE TABLE cursors_new(
                    author TEXT NOT NULL,
                    room TEXT NOT NULL DEFAULT '{DEFAULT_ROOM}',
                    last_seen_id INTEGER NOT NULL DEFAULT 0,
                    PRIMARY KEY(author, room)
                );
                INSERT OR REPLACE INTO cursors_new(author, room, last_seen_id)
                    SELECT author, '{DEFAULT_ROOM}', last_seen_id FROM cursors;
                DROP TABLE cursors;
                ALTER TABLE cursors_new RENAME TO cursors;
                """
            )


def _now_human() -> str:
    return dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M (%Z)")


def _now_iso_seconds() -> str:
    return dt.datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")


def _require_nonempty(value: str, field_name: str) -> str:
    cleaned = value.strip()
    if not cleaned:
        raise ValueError(f"{field_name} must be non-empty")
    return cleaned


def _validate_room(value: str | None, default: str = DEFAULT_ROOM) -> str:
    raw = (value or "").strip()
    if not raw:
        return default
    if not ROOM_NAME_RE.match(raw):
        raise ValueError(f"room must match {ROOM_NAME_RE.pattern} (got {value!r})")
    return raw


def _ledger_append(
    ts: str, author: str, addressee: str, topic: str, body: str, status: str, refs: list[str]
) -> bool:
    if LEDGER_PATH is None:
        return False
    try:
        LEDGER_PATH.parent.mkdir(parents=True, exist_ok=True)
        block = (
            f"\n## {ts} — {author} — {topic}\n"
            f"**To:** {addressee}\n"
            f"**Refs:** {', '.join(refs) if refs else '(none)'}\n"
            f"**Status:** {status}\n\n"
            f"{body.rstrip()}\n"
        )
        with LEDGER_PATH.open("a", encoding="utf-8") as fh:
            fh.write(block)
        return True
    except Exception as exc:
        sys.stderr.write(f"ledger append failed: {exc}\n")
        return False


def _insert_pi_message(
    topic: str,
    body: str,
    addressee: str,
    status: str,
    refs: list[str],
    room: str = DEFAULT_ROOM,
) -> dict[str, Any]:
    topic = _require_nonempty(topic, "topic")
    body = _require_nonempty(body, "body")
    addressee = _require_nonempty(addressee, "to")
    room = _validate_room(room)
    if status not in ALLOWED_STATUSES:
        raise ValueError(f"status must be one of {sorted(ALLOWED_STATUSES)}")

    ts = _now_human()
    with closing(_db()) as conn:
        cur = conn.execute(
            "INSERT INTO messages(ts, author, addressee, topic, body, status, refs, room) "
            "VALUES (?,?,?,?,?,?,?,?)",
            (ts, PI_AUTHOR, addressee, topic, body, status, json.dumps(refs), room),
        )
        msg_id = cur.lastrowid
    ledger_appended = _ledger_append(ts, PI_AUTHOR, addressee, topic, body, status, refs)
    return {"id": msg_id, "ts": ts, "room": room, "ledger_appended": ledger_appended}


def _addressee_has(addressee: str, target: str) -> bool:
    return target in {part.strip().lower() for part in addressee.split(",")}


def _approval_refs(refs_json: str) -> set[int]:
    try:
        refs = json.loads(refs_json)
    except json.JSONDecodeError:
        return set()
    out: set[int] = set()
    for ref in refs:
        text = str(ref).strip()
        for prefix in ("approval:", "mcp_msg:"):
            if text.startswith(prefix):
                try:
                    out.add(int(text.removeprefix(prefix)))
                except ValueError:
                    pass
    return out


def _approval_enabled_from_id() -> int:
    with closing(_db()) as conn:
        row = conn.execute(
            "SELECT value FROM approval_settings WHERE key='enabled_from_id'"
        ).fetchone()
    if row is None:
        return 0
    try:
        return int(row[0])
    except (TypeError, ValueError):
        return 0


def _fetch_approvals(room: str = DEFAULT_ROOM, limit: int = 50) -> list[dict[str, Any]]:
    room = _validate_room(room)
    enabled_from_id = _approval_enabled_from_id()
    with closing(_db()) as conn:
        candidate_rows = conn.execute(
            "SELECT id, ts, author, addressee, topic, body, status, refs, room "
            "FROM messages WHERE room=? AND id>? AND author<>? AND status IN (?,?) "
            "ORDER BY id DESC LIMIT 200",
            (room, enabled_from_id, PI_AUTHOR, *sorted(APPROVAL_STATUSES)),
        ).fetchall()
        resolution_rows = conn.execute(
            "SELECT refs FROM messages WHERE room=? AND author=?",
            (room, PI_AUTHOR),
        ).fetchall()

    resolved: set[int] = set()
    for (refs_json,) in resolution_rows:
        resolved.update(_approval_refs(refs_json))

    approvals = []
    for row in reversed(candidate_rows):
        msg = _message_to_dict(row)
        if msg["id"] in resolved:
            continue
        if not _addressee_has(msg["addressee"], PI_AUTHOR):
            continue
        approvals.append(msg)
    return approvals[-limit:]


def _approval_request(approval_id: int, room: str = DEFAULT_ROOM) -> dict[str, Any]:
    room = _validate_room(room)
    with closing(_db()) as conn:
        row = conn.execute(
            "SELECT id, ts, author, addressee, topic, body, status, refs, room "
            "FROM messages WHERE id=? AND room=?",
            (int(approval_id), room),
        ).fetchone()
    if row is None:
        raise ValueError(f"approval request #{approval_id} not found in room {room!r}")
    msg = _message_to_dict(row)
    if msg["status"] not in APPROVAL_STATUSES or not _addressee_has(msg["addressee"], PI_AUTHOR):
        raise ValueError(f"message #{approval_id} is not a pending approval request")
    if msg["id"] not in {a["id"] for a in _fetch_approvals(room=room, limit=200)}:
        raise ValueError(f"approval request #{approval_id} is already resolved")
    return msg


def _handle_approval_action(
    approval_id: int,
    action: str,
    note: str = "",
    room: str = DEFAULT_ROOM,
) -> dict[str, Any]:
    room = _validate_room(room)
    action = action.strip().lower()
    if action not in APPROVAL_ACTIONS:
        raise ValueError(f"action must be one of {sorted(APPROVAL_ACTIONS)}")
    request = _approval_request(approval_id, room=room)
    status, label = APPROVAL_ACTIONS[action]
    target = request["author"] if request["author"] in {"claude", "codex"} else "claude,codex"
    note = note.strip()
    body = f"[approval {action} for #{request['id']}] {label}."
    if note:
        body += f"\n\n{note}"
    return _insert_pi_message(
        topic=request["topic"],
        body=body,
        addressee=target,
        status=status,
        refs=[f"mcp_msg:{request['id']}", f"approval:{request['id']}"],
        room=room,
    )


def _insert_artifact(kind: str, path: str, note: str = "", author: str = PI_AUTHOR) -> dict[str, Any]:
    kind = _require_nonempty(kind, "kind")
    path = _require_nonempty(path, "path")
    ts = _now_human()
    with closing(_db()) as conn:
        cur = conn.execute(
            "INSERT INTO artifacts(ts, author, kind, path, note) VALUES (?,?,?,?,?)",
            (ts, author, kind, path, note or ""),
        )
        artifact_id = cur.lastrowid
    return {"id": artifact_id, "ts": ts, "author": author, "kind": kind, "path": path, "note": note or ""}


def _safe_upload_name(filename: str) -> str:
    base = pathlib.PurePath(filename or "upload.bin").name.strip()
    if not base:
        base = "upload.bin"
    safe = re.sub(r"[^A-Za-z0-9._ -]+", "_", base).strip(" .")
    return safe or "upload.bin"


def _store_upload(filename: str, data: bytes, content_type: str) -> dict[str, Any]:
    safe_name = _safe_upload_name(filename)
    stamp = dt.datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    dest = UPLOAD_DIR / f"{stamp}_{uuid.uuid4().hex[:8]}_{safe_name}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    note = f"{filename or safe_name} ({content_type or 'application/octet-stream'}, {len(data)} bytes)"
    artifact = _insert_artifact("pi_upload", str(dest), note=note)
    return {
        "artifact_id": artifact["id"],
        "path": str(dest),
        "name": filename or safe_name,
        "content_type": content_type or "application/octet-stream",
        "size": len(data),
    }


def _format_uploads_for_message(body: str, uploads: list[dict[str, Any]]) -> str:
    body = body.strip() or "첨부 파일을 공유합니다."
    if not uploads:
        return body
    lines = ["", "[첨부 파일]"]
    for upload in uploads:
        lines.append(
            "- "
            f"{upload['name']} -> {upload['path']} "
            f"(artifact #{upload['artifact_id']}, {upload['content_type']}, {upload['size']} bytes)"
        )
    return body.rstrip() + "\n" + "\n".join(lines)


def _insert_pi_message_with_uploads(
    topic: str,
    body: str,
    addressee: str,
    status: str,
    refs: list[str],
    uploads: list[dict[str, Any]],
    room: str = DEFAULT_ROOM,
) -> dict[str, Any]:
    upload_paths = [str(upload["path"]) for upload in uploads]
    combined_refs = refs + [path for path in upload_paths if path not in refs]
    return _insert_pi_message(
        topic=topic,
        body=_format_uploads_for_message(body, uploads),
        addressee=addressee,
        status=status,
        refs=combined_refs,
        room=room,
    )


def _fetch_messages(after_id: int = 0, limit: int = 200, room: str = DEFAULT_ROOM) -> list[tuple[Any, ...]]:
    room = _validate_room(room)
    with closing(_db()) as conn:
        return conn.execute(
            "SELECT id, ts, author, addressee, topic, body, status, refs, room "
            "FROM messages WHERE id > ? AND room=? ORDER BY id ASC LIMIT ?",
            (int(after_id), room, int(limit)),
        ).fetchall()


def _fetch_recent_messages(limit: int = 200, room: str = DEFAULT_ROOM) -> list[tuple[Any, ...]]:
    room = _validate_room(room)
    with closing(_db()) as conn:
        rows = conn.execute(
            "SELECT id, ts, author, addressee, topic, body, status, refs, room "
            "FROM messages WHERE room=? ORDER BY id DESC LIMIT ?",
            (room, int(limit)),
        ).fetchall()
    return list(reversed(rows))


def _fetch_rooms() -> list[str]:
    rooms = {DEFAULT_ROOM, *ROOM_PRESETS}
    with closing(_db()) as conn:
        rooms.update(row[0] for row in conn.execute("SELECT DISTINCT room FROM messages").fetchall())
        rooms.update(row[0] for row in conn.execute("SELECT DISTINCT room FROM claims").fetchall())
        agent_rows = conn.execute("SELECT agent FROM agent_status").fetchall()
    for (agent,) in agent_rows:
        agent_room = _agent_room(agent)
        if agent_room:
            rooms.add(agent_room)
    return sorted(r for r in rooms if r)


def _fetch_sidebar(room: str = DEFAULT_ROOM) -> dict[str, Any]:
    room = _validate_room(room)
    with closing(_db()) as conn:
        claims = conn.execute(
            "SELECT topic, author, claimed_at, summary FROM claims "
            "WHERE released_at IS NULL AND room=? ORDER BY claimed_at",
            (room,),
        ).fetchall()
        artifacts = conn.execute(
            "SELECT id, ts, author, kind, path, note FROM artifacts "
            "ORDER BY id DESC LIMIT 30"
        ).fetchall()
        cursors = conn.execute(
            "SELECT author, last_seen_id FROM cursors WHERE room=? ORDER BY author",
            (room,),
        ).fetchall()
        agents = conn.execute(
            "SELECT agent, percent, activity, topic, updated_at, last_active_at "
            "FROM agent_status ORDER BY agent"
        ).fetchall()
    agents = [
        row for row in agents
        if not _agent_room(row[0]) or _agent_room(row[0]) == room
    ]
    names = {row[0] for row in agents}
    if f"{room}-claude" in names or f"{room}-codex" in names:
        legacy_room_agents = {
            f"{room}-claude-chat",
            f"{room}-claude-work",
            f"{room}-codex-chat",
            f"{room}-codex-work",
            "claude-chat",
            "claude-work",
            "codex-chat",
            "codex-work",
            "claude",
            "codex",
        }
        agents = [row for row in agents if row[0] not in legacy_room_agents]
    return {"claims": claims, "artifacts": artifacts, "cursors": cursors, "agents": agents}


def _parse_iso(ts: str | None, tzinfo) -> dt.datetime | None:
    if not ts:
        return None
    try:
        return dt.datetime.strptime(ts, "%Y-%m-%dT%H:%M:%S%z")
    except (ValueError, TypeError):
        try:
            parsed = dt.datetime.strptime(ts[:16], "%Y-%m-%d %H:%M")
            return parsed.replace(tzinfo=tzinfo)
        except (ValueError, TypeError):
            return None


def _agents_to_dict(agents: list[tuple[Any, ...]]) -> list[dict[str, Any]]:
    out = []
    now = dt.datetime.now().astimezone()
    for row in agents:
        agent, percent, activity, topic, updated_at = row[:5]
        last_active_at = row[5] if len(row) > 5 else None
        parsed = _parse_iso(updated_at, now.tzinfo)
        age_seconds = max(0, int((now - parsed).total_seconds())) if parsed else None
        active_parsed = _parse_iso(last_active_at, now.tzinfo)
        active_age_seconds = (
            max(0, int((now - active_parsed).total_seconds())) if active_parsed else None
        )
        out.append(
            {
                "agent": agent,
                "percent": int(percent),
                "activity": activity or "",
                "topic": topic or "",
                "updated_at": updated_at,
                "age_seconds": age_seconds,
                "last_active_at": last_active_at,
                "active_age_seconds": active_age_seconds,
            }
        )
    return out


_AGENT_VALUES = {
    "claude", "codex", "pi",
    "claude-chat", "claude-work",
    "codex-chat", "codex-work",
}
WORKROOM_RE = ROOM_NAME_RE
_PREFIXED_AGENT_RE = re.compile(
    r"^(?P<room>[a-z0-9][a-z0-9-]{0,31})-"
    r"(?P<logical>claude|codex|claude-chat|claude-work|codex-chat|codex-work)$"
)


def _agent_logical(agent: str) -> str:
    match = _PREFIXED_AGENT_RE.match(agent)
    if match:
        return match.group("logical")
    return agent


def _agent_room(agent: str) -> str:
    match = _PREFIXED_AGENT_RE.match(agent)
    return match.group("room") if match else ""


def _validate_agent(agent: str) -> str:
    raw = _require_nonempty(agent, "agent")
    lowered = raw.lower()
    if lowered in _AGENT_VALUES:
        return lowered
    if raw != lowered:
        raise ValueError(
            "agent must be one of claude/codex/pi or <workroom>-<agent>"
        )
    match = _PREFIXED_AGENT_RE.match(raw)
    if match and WORKROOM_RE.match(match.group("room")):
        return raw
    raise ValueError(
        "agent must be one of claude/codex/pi or <workroom>-<agent>"
    )


def _agent_base(agent: str) -> str:
    """Map a pane-aware agent value to its base agent (claude/codex/pi)."""
    logical = _agent_logical(agent)
    if logical.startswith("claude"):
        return "claude"
    if logical.startswith("codex"):
        return "codex"
    return logical


def _upsert_agent_status(
    agent: str, percent: int, activity: str, topic: str = "", heartbeat: bool = False
) -> dict[str, Any]:
    agent = _validate_agent(agent)
    try:
        percent_i = int(percent)
    except (TypeError, ValueError) as exc:
        raise ValueError("percent must be an integer 0..100") from exc
    if percent_i < 0 or percent_i > 100:
        raise ValueError("percent must be in [0, 100]")
    activity = (activity or "").strip()[:240]
    topic = (topic or "").strip()[:120]
    ts = _now_iso_seconds()
    # last_active_at advances only on real LLM/operator posts. Daemon
    # heartbeat refreshes (heartbeat=true) bump updated_at so the
    # freshness color stays green, but leave last_active_at untouched so
    # the sidebar can flag panes that haven't actually moved.
    with closing(_db()) as conn:
        if heartbeat:
            conn.execute(
                "INSERT INTO agent_status(agent, percent, activity, topic, updated_at, last_active_at) "
                "VALUES(?,?,?,?,?,NULL) "
                "ON CONFLICT(agent) DO UPDATE SET "
                "percent=excluded.percent, activity=excluded.activity, "
                "topic=excluded.topic, updated_at=excluded.updated_at",
                (agent, percent_i, activity, topic, ts),
            )
        else:
            conn.execute(
                "INSERT INTO agent_status(agent, percent, activity, topic, updated_at, last_active_at) "
                "VALUES(?,?,?,?,?,?) "
                "ON CONFLICT(agent) DO UPDATE SET "
                "percent=excluded.percent, activity=excluded.activity, "
                "topic=excluded.topic, updated_at=excluded.updated_at, "
                "last_active_at=excluded.last_active_at",
                (agent, percent_i, activity, topic, ts, ts),
            )
    return {
        "agent": agent,
        "percent": percent_i,
        "activity": activity,
        "topic": topic,
        "updated_at": ts,
        "heartbeat": bool(heartbeat),
    }


def _message_to_dict(row: tuple[Any, ...]) -> dict[str, Any]:
    msg_id, ts, author, addressee, topic, body, status, refs_json = row[:8]
    room = row[8] if len(row) > 8 else DEFAULT_ROOM
    try:
        refs = json.loads(refs_json)
    except json.JSONDecodeError:
        refs = []
    return {
        "id": int(msg_id),
        "ts": ts,
        "author": author,
        "addressee": addressee,
        "topic": topic,
        "body": body,
        "status": status,
        "refs": refs,
        "room": room,
    }


def _sidebar_to_dict(sidebar: dict[str, Any]) -> dict[str, Any]:
    return {
        "claims": [
            {"topic": t, "author": a, "claimed_at": c, "summary": s or ""}
            for (t, a, c, s) in sidebar["claims"]
        ],
        "artifacts": [
            {"id": int(i), "ts": ts, "author": a, "kind": k, "path": p, "note": n or ""}
            for (i, ts, a, k, p, n) in sidebar["artifacts"]
        ],
        "cursors": {a: int(last) for (a, last) in sidebar["cursors"]},
        "agents": _agents_to_dict(sidebar.get("agents", [])),
    }


def _esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def _parse_refs(raw: str) -> list[str]:
    return [part.strip() for part in raw.split(",") if part.strip()]


def _render_login(error: str = "") -> str:
    error_html = f"<p class='error'>{_esc(error)}</p>" if error else ""
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ACS Collab Room Login</title>
  <style>{CSS}</style>
</head>
<body class="login">
  <main class="login-panel">
    <h1>ACS Collab Room</h1>
    <p>Enter the room token for the shared PI / Claude / Codex channel.</p>
    {error_html}
    <form method="post" action="/login">
      <input type="password" name="token" autocomplete="current-password" autofocus>
      <button type="submit">Enter</button>
    </form>
  </main>
</body>
</html>"""


def _render_room(initial_state: dict[str, Any], flash: str = "") -> str:
    """Initial server-side render. JS takes over polling after first paint."""
    current_room = _validate_room(str(initial_state.get("room") or DEFAULT_ROOM))
    rooms = sorted(set(str(r) for r in initial_state.get("rooms", []) if r) | {current_room, DEFAULT_ROOM})
    room_options = "\n".join(
        f"<option value='{_esc(room)}'{' selected' if room == current_room else ''}>{_esc(room)}</option>"
        for room in rooms
    )
    status_options = "\n".join(
        f"<option value='{_esc(status)}'{' selected' if status == 'open-question' else ''}>{_esc(status)}</option>"
        for status in sorted(ALLOWED_STATUSES)
    )
    addressee_options = "\n".join(
        f"<option value='{_esc(preset)}'>{_esc(preset)}</option>"
        for preset in ADDRESSEE_PRESETS
    )
    flash_html = f"<div class='flash'>{_esc(flash)}</div>" if flash else ""
    initial_json = json.dumps(initial_state)
    initial_b64 = _esc(initial_json)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ACS Collab Room</title>
  <style>{CSS}</style>
</head>
<body>
  <header class="topbar">
    <div class="brand">
      <h1>ACS Collab</h1>
      <span class="cursor-strip" id="cursor-strip"></span>
    </div>
    <div class="sync-strip" id="sync-strip" aria-live="polite">Live</div>
    <nav>
      <button type="button" id="approvals-toggle" class="ghost approval-nav">
        Approvals <span id="approvals-count">0</span>
      </button>
      <label class="room-picker">
        <span>Room</span>
        <select id="room-select">{room_options}</select>
      </label>
      <a href="/logout" class="ghost">Logout</a>
    </nav>
  </header>
  <main class="layout">
    <section class="chat" aria-label="Message timeline">
      {flash_html}
      <div class="messages" id="messages" role="log" aria-live="polite"></div>
      <button type="button" id="jump-bottom" class="jump hidden" aria-label="Jump to newest">↓ <span id="jump-count">new</span></button>
      <form class="composer" id="composer" autocomplete="off" enctype="multipart/form-data">
        <div class="composer-meta">
          <label class="meta-field">
            <span>To</span>
            <select name="to" id="composer-to">{addressee_options}</select>
          </label>
          <label class="meta-field">
            <span>Topic</span>
            <input name="topic" id="composer-topic" value="v2-layer-1" required>
          </label>
          <label class="meta-field">
            <span>Status</span>
            <select name="status" id="composer-status">{status_options}</select>
          </label>
          <label class="meta-field grow">
            <span>Refs (optional)</span>
            <input name="refs" id="composer-refs" placeholder="comma-separated paths/topics">
          </label>
          <label class="meta-field file-field">
            <span>Files</span>
            <input type="file" name="attachments" id="composer-files" multiple>
          </label>
        </div>
        <div class="composer-row">
          <textarea name="body" id="composer-body" rows="2"
            placeholder="메시지를 입력하세요. Enter=전송, Shift+Enter=줄바꿈, @claude / @codex / @both 멘션 가능."
            required></textarea>
          <button type="submit" id="composer-send">Send</button>
        </div>
      </form>
    </section>
  </main>
  <section class="approvals-drawer hidden" id="approvals-drawer" aria-label="Approvals">
    <div class="approvals-panel">
      <div class="approvals-head">
        <div>
          <h2>Approvals</h2>
          <p>PI 승인이 필요한 decision-needed / blocker 요청만 모읍니다.</p>
        </div>
        <button type="button" id="approvals-close" class="ghost">Close</button>
      </div>
      <div class="approvals-list" id="approvals-list"></div>
    </div>
  </section>
  <script id="initial-state" type="application/json">{initial_b64}</script>
  <script>{JS}</script>
</body>
</html>"""


class RoomHandler(BaseHTTPRequestHandler):
    server_version = "ACSCollabRoom/0.2"

    def do_GET(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        if path == "/logout":
            self._send_text(
                _render_login("Logged out."),
                headers={"Set-Cookie": "collab_room_token=; Path=/; Max-Age=0; SameSite=Lax"},
            )
            return
        if path == "/messages.json":
            if not self._is_authenticated():
                self._send_json({"error": "unauthorized"}, status=HTTPStatus.UNAUTHORIZED)
                return
            qs = urllib.parse.parse_qs(parsed.query)
            room = self._room_from_query(qs)
            after = 0
            try:
                after = int(qs.get("after", ["0"])[0])
            except (TypeError, ValueError):
                after = 0
            payload = self._build_payload(after_id=after, room=room)
            self._send_json(payload)
            return
        if path == "/approvals.json":
            if not self._is_authenticated():
                self._send_json({"error": "unauthorized"}, status=HTTPStatus.UNAUTHORIZED)
                return
            qs = urllib.parse.parse_qs(parsed.query)
            room = self._room_from_query(qs)
            self._send_json({"room": room, "approvals": _fetch_approvals(room=room)})
            return
        if path != "/":
            self.send_error(HTTPStatus.NOT_FOUND, "not found")
            return
        if not self._is_authenticated():
            self._send_text(_render_login())
            return
        qs = urllib.parse.parse_qs(parsed.query)
        room = self._room_from_query(qs)
        initial = self._build_payload(after_id=0, recent=True, room=room)
        self._send_text(_render_room(initial))

    def do_POST(self) -> None:
        path = urllib.parse.urlparse(self.path).path
        try:
            if path == "/login":
                self._handle_login()
            elif path == "/send":
                self._handle_send()
            elif path == "/agent_status":
                self._handle_agent_status()
            elif path == "/approval_action":
                self._handle_approval_action()
            else:
                self.send_error(HTTPStatus.NOT_FOUND, "not found")
        except ValueError as exc:
            if self._wants_json():
                self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
            else:
                initial = self._build_payload(after_id=0, recent=True)
                self._send_text(_render_room(initial, flash=str(exc)), status=HTTPStatus.BAD_REQUEST)

    def log_message(self, fmt: str, *args: Any) -> None:
        sys.stderr.write(f"[collab-room] {self.address_string()} - {fmt % args}\n")

    def _room_from_query(self, qs: dict[str, list[str]]) -> str:
        return _validate_room(qs.get("room", [DEFAULT_ROOM])[0])

    def _build_payload(self, after_id: int, recent: bool = False, room: str = DEFAULT_ROOM) -> dict[str, Any]:
        room = _validate_room(room)
        if recent:
            rows = _fetch_recent_messages(limit=200, room=room)
        else:
            rows = _fetch_messages(after_id=after_id, limit=200, room=room)
        sidebar = _sidebar_to_dict(_fetch_sidebar(room=room))
        return {
            "room": room,
            "rooms": _fetch_rooms(),
            "messages": [_message_to_dict(r) for r in rows],
            "claims": sidebar["claims"],
            "artifacts": sidebar["artifacts"],
            "cursors": sidebar["cursors"],
            "agents": sidebar["agents"],
            "approvals": _fetch_approvals(room=room),
        }

    def _handle_login(self) -> None:
        fields = self._read_form()
        token = fields.get("token", [""])[0]
        if not self._token_ok(token):
            self._send_text(_render_login("Invalid token."), status=HTTPStatus.UNAUTHORIZED)
            return
        self._redirect("/", cookie=f"collab_room_token={urllib.parse.quote(token)}; Path=/; SameSite=Lax")

    def _handle_send(self) -> None:
        if not self._is_authenticated():
            if self._wants_json():
                self._send_json({"error": "unauthorized"}, status=HTTPStatus.UNAUTHORIZED)
            else:
                self._send_text(_render_login("Please enter the room token."), status=HTTPStatus.UNAUTHORIZED)
            return
        if self._is_json_request():
            payload = self._read_json()
            room = _validate_room(str(payload.get("room", DEFAULT_ROOM)))
            result = _insert_pi_message(
                topic=str(payload.get("topic", "")),
                body=str(payload.get("body", "")),
                addressee=str(payload.get("to", "")),
                status=str(payload.get("status", "open-question")),
                refs=[str(r).strip() for r in payload.get("refs", []) if str(r).strip()],
                room=room,
            )
            self._send_json({"ok": True, **result})
            return
        if self._is_multipart_request():
            fields, uploads = self._read_multipart_form()
            room = _validate_room(fields.get("room", [DEFAULT_ROOM])[0])
            result = _insert_pi_message_with_uploads(
                topic=fields.get("topic", [""])[0],
                body=fields.get("body", [""])[0],
                addressee=fields.get("to", [""])[0],
                status=fields.get("status", ["open-question"])[0],
                refs=_parse_refs(fields.get("refs", [""])[0]),
                uploads=uploads,
                room=room,
            )
            self._send_json({"ok": True, "uploads": uploads, **result})
            return
        fields = self._read_form()
        room = _validate_room(fields.get("room", [DEFAULT_ROOM])[0])
        result = _insert_pi_message(
            topic=fields.get("topic", [""])[0],
            body=fields.get("body", [""])[0],
            addressee=fields.get("to", [""])[0],
            status=fields.get("status", ["open-question"])[0],
            refs=_parse_refs(fields.get("refs", [""])[0]),
            room=room,
        )
        self._redirect(f"/?room={urllib.parse.quote(room)}&sent={result['id']}")

    def _handle_agent_status(self) -> None:
        if not self._is_authenticated():
            self._send_json({"error": "unauthorized"}, status=HTTPStatus.UNAUTHORIZED)
            return
        payload = self._read_json()
        result = _upsert_agent_status(
            agent=str(payload.get("agent", "")),
            percent=payload.get("percent", 0),
            activity=str(payload.get("activity", "")),
            topic=str(payload.get("topic", "")),
            heartbeat=bool(payload.get("heartbeat", False)),
        )
        self._send_json({"ok": True, **result})

    def _handle_approval_action(self) -> None:
        if not self._is_authenticated():
            self._send_json({"error": "unauthorized"}, status=HTTPStatus.UNAUTHORIZED)
            return
        payload = self._read_json()
        result = _handle_approval_action(
            approval_id=int(payload.get("approval_id", 0)),
            action=str(payload.get("action", "")),
            note=str(payload.get("note", "")),
            room=str(payload.get("room", DEFAULT_ROOM)),
        )
        self._send_json({"ok": True, **result, "approvals": _fetch_approvals(room=result["room"])})

    def _read_form(self) -> dict[str, list[str]]:
        length = int(self.headers.get("Content-Length", "0"))
        if length > MAX_FORM_BYTES:
            raise ValueError("request body too large")
        body = self.rfile.read(length).decode("utf-8")
        return urllib.parse.parse_qs(body, keep_blank_values=True)

    def _read_multipart_form(self) -> tuple[dict[str, list[str]], list[dict[str, Any]]]:
        length = int(self.headers.get("Content-Length", "0"))
        if length > MAX_UPLOAD_BYTES:
            raise ValueError(f"upload request too large; limit is {MAX_UPLOAD_BYTES} bytes")
        content_type = self.headers.get("Content-Type", "")
        raw = self.rfile.read(length)
        parser = email.parser.BytesParser(policy=email.policy.default)
        msg = parser.parsebytes(
            b"Content-Type: " + content_type.encode("utf-8") + b"\r\n"
            b"MIME-Version: 1.0\r\n\r\n"
            + raw
        )
        if not msg.is_multipart():
            raise ValueError("multipart body expected")
        fields: dict[str, list[str]] = {}
        uploads: list[dict[str, Any]] = []
        for part in msg.iter_parts():
            name = part.get_param("name", header="content-disposition")
            if not name:
                continue
            filename = part.get_filename()
            payload = part.get_payload(decode=True) or b""
            if filename:
                if not payload:
                    continue
                uploads.append(
                    _store_upload(
                        filename=filename,
                        data=payload,
                        content_type=part.get_content_type(),
                    )
                )
            else:
                charset = part.get_content_charset() or "utf-8"
                fields.setdefault(name, []).append(payload.decode(charset, errors="replace"))
        return fields, uploads

    def _read_json(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        if length > MAX_FORM_BYTES:
            raise ValueError("request body too large")
        raw = self.rfile.read(length).decode("utf-8") or "{}"
        try:
            data = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid json: {exc}") from exc
        if not isinstance(data, dict):
            raise ValueError("json body must be an object")
        return data

    def _is_json_request(self) -> bool:
        return (self.headers.get("Content-Type") or "").lower().startswith("application/json")

    def _is_multipart_request(self) -> bool:
        return (self.headers.get("Content-Type") or "").lower().startswith("multipart/form-data")

    def _wants_json(self) -> bool:
        ctype = (self.headers.get("Content-Type") or "").lower()
        accept = (self.headers.get("Accept") or "").lower()
        return ctype.startswith("application/json") or "application/json" in accept

    def _is_authenticated(self) -> bool:
        if not ROOM_TOKEN:
            return True
        header_token = self.headers.get("X-Collab-Token", "")
        if header_token and self._token_ok(header_token):
            return True
        cookies = http.cookies.SimpleCookie(self.headers.get("Cookie", ""))
        cookie = cookies.get("collab_room_token")
        return bool(cookie and self._token_ok(urllib.parse.unquote(cookie.value)))

    def _token_ok(self, token: str) -> bool:
        if not ROOM_TOKEN:
            return True
        return any(hmac.compare_digest(token, t) for t in ROOM_TOKENS_ACCEPTED)

    def _send_text(
        self,
        body: str,
        status: HTTPStatus = HTTPStatus.OK,
        headers: dict[str, str] | None = None,
    ) -> None:
        payload = body.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        for key, value in (headers or {}).items():
            self.send_header(key, value)
        self.end_headers()
        self.wfile.write(payload)

    def _send_json(self, data: dict[str, Any], status: HTTPStatus = HTTPStatus.OK) -> None:
        payload = json.dumps(data).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(payload)

    def _redirect(self, location: str, cookie: str | None = None) -> None:
        self.send_response(HTTPStatus.SEE_OTHER)
        self.send_header("Location", location)
        if cookie:
            self.send_header("Set-Cookie", cookie)
        self.end_headers()


CSS = """
:root {
  color-scheme: light;
  --bg: #b2c7d9;
  --panel: #ffffff;
  --ink: #1d232f;
  --muted: #6b7787;
  --line: #d9dee8;
  --accent: #ffeb33;
  --accent-ink: #1d232f;
  --pi-bubble: #ffeb33;
  --claude-bubble: #ffffff;
  --codex-bubble: #ffffff;
  --claude-tag: #246bfe;
  --codex-tag: #2a8a4d;
  --pi-tag: #d49a00;
  --shadow: 0 1px 2px rgba(0,0,0,0.06);
}
* { box-sizing: border-box; }
html, body { height: 100%; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font: 14px/1.45 -apple-system, BlinkMacSystemFont, "Segoe UI", "Apple SD Gothic Neo", "Noto Sans KR", sans-serif;
  display: flex;
  flex-direction: column;
}
h1, h2, p { margin: 0; }
h1 { font-size: 16px; font-weight: 700; }
h2 { font-size: 12px; font-weight: 700; color: var(--muted); text-transform: uppercase; letter-spacing: 0.04em; margin-bottom: 8px; }
a, button { font: inherit; }
.topbar {
  align-items: center;
  background: var(--panel);
  border-bottom: 1px solid var(--line);
  display: flex;
  gap: 16px;
  justify-content: space-between;
  padding: 10px 16px;
  z-index: 5;
}
.brand { display: flex; align-items: baseline; gap: 12px; flex-wrap: wrap; }
.cursor-strip {
  align-items: center;
  color: var(--muted);
  display: inline-flex;
  flex-wrap: wrap;
  font-size: 12px;
  gap: 4px;
  min-width: 0;
}
.cursor-strip .pill {
  background: rgba(0,0,0,0.05);
  border-radius: 10px;
  flex: 0 0 auto;
  font-variant-numeric: tabular-nums;
  padding: 2px 8px;
  white-space: nowrap;
}
.sync-strip {
  align-items: center;
  display: flex;
  flex: 1 1 auto;
  justify-content: flex-end;
  min-width: 0;
  color: var(--muted);
  font-size: 12px;
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
}
nav { display: flex; gap: 8px; }
.room-picker {
  align-items: center;
  color: var(--muted);
  display: inline-flex;
  font-size: 12px;
  gap: 6px;
}
.room-picker select {
  background: #f3f5f8;
  border: 1px solid var(--line);
  border-radius: 6px;
  color: var(--ink);
  padding: 6px 8px;
}
.ghost {
  background: transparent;
  border: 1px solid var(--line);
  border-radius: 6px;
  color: var(--ink);
  cursor: pointer;
  padding: 6px 10px;
  text-decoration: none;
}
.ghost:hover { background: rgba(0,0,0,0.04); }
.approval-nav {
  align-items: center;
  display: inline-flex;
  gap: 6px;
}
.approval-nav #approvals-count {
  background: #c92a2a;
  border-radius: 10px;
  color: #fff;
  display: inline-block;
  font-size: 11px;
  font-weight: 700;
  min-width: 18px;
  padding: 1px 6px;
  text-align: center;
}
.approval-nav.clean #approvals-count {
  background: rgba(0,0,0,0.12);
  color: var(--muted);
}
.layout {
  display: grid;
  grid-template-columns: minmax(0, 1fr);
  flex: 1 1 auto;
  min-height: 0;
}
.chat {
  display: flex;
  flex-direction: column;
  min-height: 0;
  position: relative;
}
.messages {
  flex: 1 1 auto;
  overflow-y: auto;
  padding: 16px 18px 8px;
  display: flex;
  flex-direction: column;
  gap: 4px;
  scroll-behavior: smooth;
}
.day-divider {
  align-self: center;
  background: rgba(0,0,0,0.08);
  border-radius: 10px;
  color: #fff;
  font-size: 11px;
  margin: 8px 0;
  padding: 3px 10px;
}
.msg-row {
  display: flex;
  gap: 8px;
  margin-top: 2px;
  max-width: 100%;
}
.msg-row.same-author { margin-top: 1px; }
.msg-row.left { justify-content: flex-start; }
.msg-row.right { justify-content: flex-end; }
.msg-block { display: flex; flex-direction: column; max-width: 78%; min-width: 0; }
.msg-row.right .msg-block { align-items: flex-end; }
.msg-row.left .msg-block { align-items: flex-start; }
.author-line {
  color: var(--ink);
  font-size: 12px;
  font-weight: 600;
  margin: 6px 4px 2px;
}
.author-line.claude { color: var(--claude-tag); }
.author-line.codex { color: var(--codex-tag); }
.author-line.pi { color: var(--pi-tag); }
.bubble-line { display: flex; align-items: flex-end; gap: 6px; max-width: 100%; }
.msg-row.right .bubble-line { flex-direction: row-reverse; }
.bubble {
  background: var(--claude-bubble);
  border-radius: 14px;
  box-shadow: var(--shadow);
  padding: 8px 12px;
  white-space: pre-wrap;
  word-break: break-word;
  max-width: 100%;
}
.bubble.pi { background: var(--pi-bubble); }
.bubble.codex { background: var(--codex-bubble); }
.bubble.claude { background: var(--claude-bubble); }
.bubble .topic-tag {
  color: var(--muted);
  display: block;
  font-size: 11px;
  margin-bottom: 4px;
}
.bubble .topic-tag .status-chip {
  background: rgba(0,0,0,0.06);
  border-radius: 8px;
  margin-left: 4px;
  padding: 1px 6px;
  text-transform: lowercase;
}
.bubble .mention { background: rgba(36,107,254,0.12); border-radius: 4px; color: var(--claude-tag); font-weight: 600; padding: 0 3px; }
.bubble .mention.codex { background: rgba(42,138,77,0.14); color: var(--codex-tag); }
.bubble .mention.both { background: rgba(212,154,0,0.18); color: #8a6300; }
.bubble .refs { color: var(--muted); font-size: 11px; margin-top: 6px; }
.bubble .refs code { background: rgba(0,0,0,0.05); border-radius: 4px; padding: 1px 4px; }
.meta-rail {
  align-self: flex-end;
  color: var(--muted);
  font-size: 10px;
  line-height: 1.1;
  text-align: right;
  white-space: nowrap;
}
.msg-row.left .meta-rail { text-align: left; }
.meta-rail .seen { color: #7aa97a; display: block; }
.meta-rail .seen.pending { color: var(--muted); }
.composer {
  background: var(--panel);
  border-top: 1px solid var(--line);
  display: flex;
  flex-direction: column;
  gap: 6px;
  padding: 10px 12px 12px;
}
.composer-meta {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
}
.meta-field { display: flex; flex-direction: column; gap: 2px; min-width: 140px; }
.meta-field.grow { flex: 1 1 240px; }
.meta-field span { color: var(--muted); font-size: 11px; font-weight: 600; }
.meta-field input, .meta-field select {
  background: #f3f5f8;
  border: 1px solid var(--line);
  border-radius: 6px;
  color: var(--ink);
  padding: 6px 8px;
}
.meta-field.file-field { flex: 1 1 260px; }
.meta-field.file-field input {
  cursor: pointer;
  max-width: 100%;
}
.composer-row { display: flex; gap: 8px; align-items: stretch; }
.composer-row textarea {
  background: #f3f5f8;
  border: 1px solid var(--line);
  border-radius: 10px;
  flex: 1 1 auto;
  font: inherit;
  max-height: 240px;
  min-height: 38px;
  padding: 9px 12px;
  resize: none;
  width: 100%;
}
.composer-row textarea:focus { outline: 2px solid var(--accent); outline-offset: -1px; }
#composer-send {
  align-self: flex-end;
  background: var(--accent);
  border: 0;
  border-radius: 10px;
  color: var(--accent-ink);
  cursor: pointer;
  font-weight: 700;
  padding: 9px 18px;
}
#composer-send:disabled { opacity: 0.5; cursor: progress; }
.jump {
  background: var(--accent);
  border: 0;
  border-radius: 18px;
  bottom: 130px;
  box-shadow: 0 2px 6px rgba(0,0,0,0.2);
  color: var(--accent-ink);
  cursor: pointer;
  font-weight: 700;
  padding: 7px 14px;
  position: absolute;
  right: 18px;
  z-index: 4;
}
.jump.hidden { display: none; }
.flash {
  background: #fff0f0;
  border: 1px solid #ffc9c9;
  border-radius: 6px;
  color: #8b1e1e;
  margin: 12px 18px 0;
  padding: 8px 10px;
}
.approvals-drawer {
  background: rgba(29,35,47,0.22);
  bottom: 0;
  display: flex;
  justify-content: flex-end;
  left: 0;
  position: fixed;
  right: 0;
  top: 0;
  z-index: 20;
}
.approvals-drawer.hidden { display: none; }
.approvals-panel {
  background: var(--panel);
  border-left: 1px solid var(--line);
  box-shadow: -8px 0 24px rgba(0,0,0,0.18);
  display: flex;
  flex-direction: column;
  min-width: 340px;
  width: min(520px, 92vw);
}
.approvals-head {
  align-items: flex-start;
  border-bottom: 1px solid var(--line);
  display: flex;
  gap: 12px;
  justify-content: space-between;
  padding: 14px 16px;
}
.approvals-head h2 {
  color: var(--ink);
  font-size: 15px;
  letter-spacing: 0;
  margin-bottom: 2px;
  text-transform: none;
}
.approvals-head p {
  color: var(--muted);
  font-size: 12px;
}
.approvals-list {
  display: flex;
  flex-direction: column;
  gap: 10px;
  overflow-y: auto;
  padding: 14px 16px 18px;
}
.approval-card {
  background: #f6f8fb;
  border: 1px solid var(--line);
  border-radius: 8px;
  display: flex;
  flex-direction: column;
  gap: 8px;
  padding: 10px 12px;
}
.approval-card .approval-meta {
  color: var(--muted);
  font-size: 11px;
}
.approval-card .approval-title {
  font-weight: 700;
}
.approval-card .approval-body {
  max-height: 180px;
  overflow: auto;
  white-space: pre-wrap;
  word-break: break-word;
}
.approval-card .approval-actions {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
}
.approval-actions button {
  border: 1px solid var(--line);
  border-radius: 6px;
  cursor: pointer;
  padding: 6px 9px;
}
.approval-actions .approve { background: #2a8a4d; color: #fff; }
.approval-actions .changes { background: var(--accent); color: var(--accent-ink); }
.approval-actions .reject { background: #fff0f0; color: #8b1e1e; }
.approval-empty {
  border: 1px dashed var(--line);
  border-radius: 8px;
  color: var(--muted);
  padding: 16px;
  text-align: center;
}
.error {
  background: #fff0f0;
  border: 1px solid #ffc9c9;
  border-radius: 6px;
  color: #8b1e1e;
  margin-bottom: 12px;
  padding: 8px 10px;
}
.login {
  display: grid;
  min-height: 100vh;
  place-items: center;
}
.login-panel {
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 10px;
  display: grid;
  gap: 12px;
  max-width: 420px;
  padding: 24px;
  width: min(92vw, 420px);
}
.login-panel input {
  border: 1px solid var(--line);
  border-radius: 6px;
  padding: 8px 10px;
}
.login-panel button {
  background: var(--accent);
  border: 0;
  border-radius: 6px;
  color: var(--accent-ink);
  cursor: pointer;
  font-weight: 700;
  padding: 9px 14px;
}
@media (max-width: 900px) {
  .layout { grid-template-columns: 1fr; }
}
"""


JS = r"""
(function () {
  var initialEl = document.getElementById("initial-state");
  var initial = {messages: [], claims: [], artifacts: [], approvals: [], cursors: {}};
  if (initialEl) {
    try { initial = JSON.parse(initialEl.textContent || "{}"); } catch (e) {}
  }
  var ROOM_STORAGE_KEY = "acs-collab-current-room";
  var ROOM_RE = /^[a-z0-9][a-z0-9-]{0,31}$/;
  function validRoom(value) {
    value = String(value || "");
    return ROOM_RE.test(value) ? value : "";
  }
  function roomFromUrl() {
    try {
      return validRoom(new URLSearchParams(window.location.search).get("room"));
    } catch (e) {
      return "";
    }
  }
  function storedRoom() {
    try {
      return validRoom(window.localStorage.getItem(ROOM_STORAGE_KEY));
    } catch (e) {
      return "";
    }
  }
  function persistRoom(room) {
    room = validRoom(room);
    if (!room) return;
    try { window.localStorage.setItem(ROOM_STORAGE_KEY, room); } catch (e) {}
  }
  var selectedRoom = roomFromUrl() || storedRoom() || initial.room || "design-discussion";
  selectedRoom = validRoom(selectedRoom) || "design-discussion";
  persistRoom(selectedRoom);
  var state = {
    lastId: 0,
    messages: [],
    cursors: initial.cursors || {},
    agents: initial.agents || [],
    approvals: initial.approvals || [],
    agentsFetchedAt: Date.now(),
    currentRoom: selectedRoom,
    rooms: initial.rooms || [],
    pending: 0,
    pollFailures: 0,
  };
  var messagesEl = document.getElementById("messages");
  var cursorStrip = document.getElementById("cursor-strip");
  var syncStrip = document.getElementById("sync-strip");
  var approvalsToggle = document.getElementById("approvals-toggle");
  var approvalsCount = document.getElementById("approvals-count");
  var approvalsDrawer = document.getElementById("approvals-drawer");
  var approvalsClose = document.getElementById("approvals-close");
  var approvalsList = document.getElementById("approvals-list");
  var claimsList = null;
  var artifactsList = null;
  var jumpBtn = document.getElementById("jump-bottom");
  var jumpCount = document.getElementById("jump-count");
  var composer = document.getElementById("composer");
  var bodyEl = document.getElementById("composer-body");
  var fileEl = document.getElementById("composer-files");
  var sendBtn = document.getElementById("composer-send");
  var roomSelect = document.getElementById("room-select");

  var SELF_AUTHOR = "pi";
  var POLL_FAST_MS = 600;
  var POLL_IDLE_MS = 1500;
  var POLL_HIDDEN_MS = 2500;
  var POLL_BACKOFF_MAX_MS = 8000;
  var pollTimer = null;
  var pollInFlight = false;

  function isNearBottom() {
    return messagesEl.scrollHeight - messagesEl.scrollTop - messagesEl.clientHeight < 80;
  }

  function scrollToBottom() {
    messagesEl.scrollTop = messagesEl.scrollHeight;
    state.pending = 0;
    jumpBtn.classList.add("hidden");
  }

  function dayLabel(ts) {
    if (!ts) return "";
    var m = ts.match(/^(\d{4}-\d{2}-\d{2})/);
    return m ? m[1] : "";
  }

  function timeShort(ts) {
    var m = ts && ts.match(/(\d{2}:\d{2})/);
    return m ? m[1] : "";
  }

  function authorSide(author) {
    return author === SELF_AUTHOR ? "right" : "left";
  }

  function authorClass(author) {
    if (author === "pi") return "pi";
    if (author === "claude") return "claude";
    if (author === "codex") return "codex";
    return "other";
  }

  var MENTION_RE = /(@(?:claude|codex|both|pi))\b/gi;

  function appendBodyWithMentions(target, body) {
    var lines = String(body).split(/\n/);
    lines.forEach(function (line, idx) {
      var lastIndex = 0;
      MENTION_RE.lastIndex = 0;
      var match;
      while ((match = MENTION_RE.exec(line)) !== null) {
        if (match.index > lastIndex) {
          target.appendChild(document.createTextNode(line.slice(lastIndex, match.index)));
        }
        var span = document.createElement("span");
        var who = match[1].slice(1).toLowerCase();
        span.className = "mention " + who;
        span.textContent = match[1];
        target.appendChild(span);
        lastIndex = match.index + match[1].length;
      }
      if (lastIndex < line.length) {
        target.appendChild(document.createTextNode(line.slice(lastIndex)));
      }
      if (idx < lines.length - 1) {
        target.appendChild(document.createElement("br"));
      }
    });
  }

  function buildSeenLabel(msg) {
    var others = ["claude", "codex"].filter(function (n) { return n !== msg.author; });
    var seenBy = others.filter(function (n) { return (state.cursors[n] || 0) >= msg.id; });
    var span = document.createElement("span");
    if (seenBy.length === others.length) {
      span.className = "seen";
      span.textContent = "✓ seen";
    } else if (seenBy.length > 0) {
      span.className = "seen";
      span.textContent = "✓ " + seenBy.join(",");
    } else {
      span.className = "seen pending";
      span.textContent = "· not seen";
    }
    return span;
  }

  function buildMessageNode(msg, prevAuthor) {
    var row = document.createElement("div");
    var side = authorSide(msg.author);
    row.className = "msg-row " + side + (msg.author === prevAuthor ? " same-author" : "");
    row.dataset.id = msg.id;
    row.dataset.author = msg.author;

    var block = document.createElement("div");
    block.className = "msg-block";

    if (msg.author !== prevAuthor) {
      var nameLine = document.createElement("div");
      nameLine.className = "author-line " + authorClass(msg.author);
      nameLine.textContent = msg.author + " → " + msg.addressee;
      block.appendChild(nameLine);
    }

    var line = document.createElement("div");
    line.className = "bubble-line";

    var bubble = document.createElement("div");
    bubble.className = "bubble " + authorClass(msg.author);

    var topicTag = document.createElement("div");
    topicTag.className = "topic-tag";
    topicTag.textContent = "#" + msg.id + " · " + msg.topic;
    var statusChip = document.createElement("span");
    statusChip.className = "status-chip";
    statusChip.textContent = msg.status || "FYI";
    topicTag.appendChild(statusChip);
    bubble.appendChild(topicTag);

    var bodyWrap = document.createElement("div");
    appendBodyWithMentions(bodyWrap, msg.body || "");
    bubble.appendChild(bodyWrap);

    if (msg.refs && msg.refs.length) {
      var refsP = document.createElement("div");
      refsP.className = "refs";
      refsP.appendChild(document.createTextNode("refs: "));
      msg.refs.forEach(function (r, i) {
        var c = document.createElement("code");
        c.textContent = r;
        refsP.appendChild(c);
        if (i < msg.refs.length - 1) refsP.appendChild(document.createTextNode(" "));
      });
      bubble.appendChild(refsP);
    }

    var meta = document.createElement("div");
    meta.className = "meta-rail";
    var t = document.createElement("span");
    t.textContent = timeShort(msg.ts);
    meta.appendChild(t);
    meta.appendChild(buildSeenLabel(msg));

    line.appendChild(bubble);
    line.appendChild(meta);
    block.appendChild(line);
    row.appendChild(block);
    return row;
  }

  function renderMessages(allMessages) {
    messagesEl.innerHTML = "";
    var prevDay = "";
    var prevAuthor = "";
    allMessages.forEach(function (msg) {
      var d = dayLabel(msg.ts);
      if (d && d !== prevDay) {
        var div = document.createElement("div");
        div.className = "day-divider";
        div.textContent = d;
        messagesEl.appendChild(div);
        prevDay = d;
        prevAuthor = "";
      }
      messagesEl.appendChild(buildMessageNode(msg, prevAuthor));
      prevAuthor = msg.author;
    });
    state.messages = allMessages.slice();
    state.lastId = allMessages.length ? allMessages[allMessages.length - 1].id : state.lastId;
  }

  function appendNewMessages(newMessages) {
    if (!newMessages.length) return;
    var nearBottom = isNearBottom();
    var prevAuthor = state.messages.length ? state.messages[state.messages.length - 1].author : "";
    var prevDay = state.messages.length ? dayLabel(state.messages[state.messages.length - 1].ts) : "";
    newMessages.forEach(function (msg) {
      var d = dayLabel(msg.ts);
      if (d && d !== prevDay) {
        var div = document.createElement("div");
        div.className = "day-divider";
        div.textContent = d;
        messagesEl.appendChild(div);
        prevDay = d;
        prevAuthor = "";
      }
      messagesEl.appendChild(buildMessageNode(msg, prevAuthor));
      prevAuthor = msg.author;
      state.messages.push(msg);
      if (msg.id > state.lastId) state.lastId = msg.id;
    });
    if (nearBottom) {
      scrollToBottom();
    } else {
      state.pending += newMessages.length;
      jumpCount.textContent = state.pending + " new";
      jumpBtn.classList.remove("hidden");
    }
  }

  function rerenderSeen() {
    var rows = messagesEl.querySelectorAll(".msg-row");
    rows.forEach(function (row, idx) {
      var msgId = parseInt(row.dataset.id || "0", 10);
      if (!msgId) return;
      var msg = null;
      for (var i = 0; i < state.messages.length; i++) {
        if (state.messages[i].id === msgId) { msg = state.messages[i]; break; }
      }
      if (!msg) return;
      var rail = row.querySelector(".meta-rail");
      if (!rail) return;
      var oldSeen = rail.querySelector(".seen");
      if (oldSeen) oldSeen.remove();
      rail.appendChild(buildSeenLabel(msg));
    });
  }

  function renderCursors() {
    cursorStrip.innerHTML = "";
    var keys = Object.keys(state.cursors).sort();
    if (!keys.length) {
      cursorStrip.textContent = "no cursors yet";
      return;
    }
    keys.forEach(function (k) {
      var pill = document.createElement("span");
      pill.className = "pill";
      pill.textContent = k + ": #" + state.cursors[k];
      cursorStrip.appendChild(pill);
    });
  }

  function renderSync(agents) {
    if (agents) state.agents = agents;
    if (!syncStrip) return;
    var agentCount = state.agents ? state.agents.length : 0;
    var label = "Live · " + state.currentRoom + " · #" + (state.lastId || 0);
    if (agentCount) label += " · " + agentCount + " agents linked";
    if (state.pollFailures) label = "Reconnecting · " + state.currentRoom;
    syncStrip.textContent = label;
  }

  function renderApprovals(approvals) {
    if (approvals) state.approvals = approvals;
    var count = state.approvals ? state.approvals.length : 0;
    if (approvalsCount) approvalsCount.textContent = String(count);
    if (approvalsToggle) approvalsToggle.classList.toggle("clean", count === 0);
    if (!approvalsList) return;
    approvalsList.innerHTML = "";
    if (!count) {
      var empty = document.createElement("div");
      empty.className = "approval-empty";
      empty.textContent = "No pending approvals";
      approvalsList.appendChild(empty);
      return;
    }
    state.approvals.forEach(function (a) {
      var card = document.createElement("div");
      card.className = "approval-card";

      var title = document.createElement("div");
      title.className = "approval-title";
      title.textContent = "#" + a.id + " · " + a.topic + " · " + a.status;
      card.appendChild(title);

      var meta = document.createElement("div");
      meta.className = "approval-meta";
      meta.textContent = a.author + " → " + a.addressee + " · " + a.ts;
      card.appendChild(meta);

      var body = document.createElement("div");
      body.className = "approval-body";
      body.textContent = a.body || "";
      card.appendChild(body);

      if (a.refs && a.refs.length) {
        var refs = document.createElement("div");
        refs.className = "approval-meta";
        refs.textContent = "refs: " + a.refs.join(" ");
        card.appendChild(refs);
      }

      var actions = document.createElement("div");
      actions.className = "approval-actions";
      [
        ["approve", "Approve"],
        ["changes", "Needs changes"],
        ["reject", "Reject"],
      ].forEach(function (item) {
        var btn = document.createElement("button");
        btn.type = "button";
        btn.className = item[0];
        btn.textContent = item[1];
        btn.addEventListener("click", function () { submitApproval(a.id, item[0]); });
        actions.appendChild(btn);
      });
      card.appendChild(actions);
      approvalsList.appendChild(card);
    });
  }

  function submitApproval(approvalId, action) {
    var note = window.prompt("Optional note for " + action + " #" + approvalId + ":", "") || "";
    fetch("/approval_action", {
      method: "POST",
      credentials: "same-origin",
      headers: {"Content-Type": "application/json", "Accept": "application/json"},
      body: JSON.stringify({
        approval_id: approvalId,
        action: action,
        note: note,
        room: state.currentRoom,
      }),
    })
      .then(function (r) { return r.json().then(function (j) { return {ok: r.ok, body: j}; }); })
      .then(function (res) {
        if (!res.ok) {
          alert((res.body && res.body.error) || "approval action failed");
          return;
        }
        renderApprovals(res.body.approvals || []);
        poll(true);
      })
      .catch(function (err) { alert("approval action failed: " + err); });
  }

  function renderClaims(claims) {
    if (!claimsList) return;
    claimsList.innerHTML = "";
    if (!claims.length) {
      var li = document.createElement("li");
      li.className = "empty";
      li.textContent = "None";
      claimsList.appendChild(li);
      return;
    }
    claims.forEach(function (c) {
      var li = document.createElement("li");
      var topic = document.createElement("strong");
      topic.textContent = c.topic;
      li.appendChild(topic);
      var meta = document.createElement("span");
      meta.className = "meta";
      meta.textContent = c.author + " · " + c.claimed_at;
      li.appendChild(meta);
      if (c.summary) {
        var sm = document.createElement("span");
        sm.className = "meta";
        sm.textContent = c.summary;
        li.appendChild(sm);
      }
      claimsList.appendChild(li);
    });
  }

  function renderArtifacts(artifacts) {
    if (!artifactsList) return;
    artifactsList.innerHTML = "";
    if (!artifacts.length) {
      var li = document.createElement("li");
      li.className = "empty";
      li.textContent = "None";
      artifactsList.appendChild(li);
      return;
    }
    artifacts.forEach(function (a) {
      var li = document.createElement("li");
      var head = document.createElement("strong");
      head.textContent = a.kind + " #" + a.id;
      li.appendChild(head);
      var path = document.createElement("code");
      path.textContent = a.path;
      li.appendChild(document.createElement("br"));
      li.appendChild(path);
      var meta = document.createElement("span");
      meta.className = "meta";
      meta.textContent = a.author + " · " + a.ts;
      li.appendChild(meta);
      if (a.note) {
        var nm = document.createElement("span");
        nm.className = "meta";
        nm.textContent = a.note;
        li.appendChild(nm);
      }
      artifactsList.appendChild(li);
    });
  }

  function pollDelay(hasNewMessages) {
    if (state.pollFailures) {
      return Math.min(POLL_BACKOFF_MAX_MS, POLL_IDLE_MS * Math.pow(2, state.pollFailures - 1));
    }
    if (document.hidden) return POLL_HIDDEN_MS;
    if (hasNewMessages || document.hasFocus()) return POLL_FAST_MS;
    return POLL_IDLE_MS;
  }

  function schedulePoll(delay) {
    if (pollTimer) window.clearTimeout(pollTimer);
    pollTimer = window.setTimeout(function () { poll(false); }, delay);
  }

  function poll(force) {
    if (pollInFlight) {
      if (force) schedulePoll(100);
      return;
    }
    if (pollTimer) {
      window.clearTimeout(pollTimer);
      pollTimer = null;
    }
    pollInFlight = true;
    var url = "/messages.json?after=" + encodeURIComponent(state.lastId) +
      "&room=" + encodeURIComponent(state.currentRoom);
    fetch(url, {credentials: "same-origin", headers: {"Accept": "application/json"}})
      .then(function (r) {
        if (!r.ok) throw new Error("status " + r.status);
        return r.json();
      })
      .then(function (data) {
        state.pollFailures = 0;
        state.currentRoom = data.room || state.currentRoom;
        state.rooms = data.rooms || state.rooms;
        var prevCursors = JSON.stringify(state.cursors);
        state.cursors = data.cursors || {};
        var hasNewMessages = !!(data.messages && data.messages.length);
        if (data.messages && data.messages.length) {
          appendNewMessages(data.messages);
        }
        if (JSON.stringify(state.cursors) !== prevCursors) {
          rerenderSeen();
          renderCursors();
        }
        renderClaims(data.claims || []);
        renderArtifacts(data.artifacts || []);
        state.agentsFetchedAt = Date.now();
        renderSync(data.agents || []);
        renderApprovals(data.approvals || []);
        schedulePoll(pollDelay(hasNewMessages));
      })
      .catch(function (err) {
        // network blip — try again next tick
        state.pollFailures = (state.pollFailures || 0) + 1;
        if (window.console && console.warn) console.warn("poll failed", err);
        schedulePoll(pollDelay(false));
      })
      .finally(function () {
        pollInFlight = false;
      });
  }

  function autoresize(el) {
    el.style.height = "auto";
    el.style.height = Math.min(240, el.scrollHeight) + "px";
  }

  function sendMessage() {
    var body = bodyEl.value.trim();
    var files = fileEl && fileEl.files ? Array.prototype.slice.call(fileEl.files) : [];
    if (!body && !files.length) return;
    sendBtn.disabled = true;
    var refsText = document.getElementById("composer-refs").value || "";
    var payload;
    var fetchOptions;
    if (files.length) {
      payload = new FormData();
      payload.append("to", document.getElementById("composer-to").value);
      payload.append("topic", document.getElementById("composer-topic").value);
      payload.append("status", document.getElementById("composer-status").value);
      payload.append("body", bodyEl.value);
      payload.append("refs", refsText);
      payload.append("room", state.currentRoom);
      files.forEach(function (file) { payload.append("attachments", file, file.name); });
      fetchOptions = {
        method: "POST",
        credentials: "same-origin",
        headers: {"Accept": "application/json"},
        body: payload,
      };
    } else {
      payload = {
        to: document.getElementById("composer-to").value,
        topic: document.getElementById("composer-topic").value,
        status: document.getElementById("composer-status").value,
        body: bodyEl.value,
        refs: refsText.split(",").map(function (s) { return s.trim(); }).filter(Boolean),
        room: state.currentRoom,
      };
      fetchOptions = {
        method: "POST",
        credentials: "same-origin",
        headers: {"Content-Type": "application/json", "Accept": "application/json"},
        body: JSON.stringify(payload),
      };
    }
    fetch("/send", {
      method: fetchOptions.method,
      credentials: fetchOptions.credentials,
      headers: fetchOptions.headers,
      body: fetchOptions.body,
    })
      .then(function (r) { return r.json().then(function (j) { return {ok: r.ok, body: j}; }); })
      .then(function (res) {
        if (!res.ok) {
          alert((res.body && res.body.error) || "send failed");
          return;
        }
        bodyEl.value = "";
        if (fileEl) fileEl.value = "";
        autoresize(bodyEl);
        poll(true);
      })
      .catch(function (err) { alert("send failed: " + err); })
      .finally(function () { sendBtn.disabled = false; bodyEl.focus(); });
  }

  composer.addEventListener("submit", function (e) {
    e.preventDefault();
    sendMessage();
  });

  if (roomSelect) {
    roomSelect.value = state.currentRoom;
    roomSelect.addEventListener("change", function () {
      var nextRoom = roomSelect.value || "design-discussion";
      nextRoom = validRoom(nextRoom) || "design-discussion";
      persistRoom(nextRoom);
      state.currentRoom = nextRoom;
      window.location.assign("/?room=" + encodeURIComponent(nextRoom));
    });
  }

  bodyEl.addEventListener("keydown", function (e) {
    if (e.key === "Enter" && !e.shiftKey && !e.isComposing) {
      e.preventDefault();
      sendMessage();
    }
  });
  bodyEl.addEventListener("input", function () { autoresize(bodyEl); });

  jumpBtn.addEventListener("click", scrollToBottom);
  messagesEl.addEventListener("scroll", function () {
    if (isNearBottom()) {
      state.pending = 0;
      jumpBtn.classList.add("hidden");
    }
  });
  if (approvalsToggle && approvalsDrawer) {
    approvalsToggle.addEventListener("click", function () {
      approvalsDrawer.classList.remove("hidden");
      renderApprovals();
    });
  }
  if (approvalsClose && approvalsDrawer) {
    approvalsClose.addEventListener("click", function () {
      approvalsDrawer.classList.add("hidden");
    });
  }
  if (approvalsDrawer) {
    approvalsDrawer.addEventListener("click", function (e) {
      if (e.target === approvalsDrawer) approvalsDrawer.classList.add("hidden");
    });
  }

  document.addEventListener("visibilitychange", function () {
    if (!document.hidden) poll(true);
  });
  window.addEventListener("focus", function () { poll(true); });

  // Initial paint from server-supplied state.
  state.cursors = initial.cursors || {};
  renderCursors();
  renderMessages(initial.messages || []);
  renderClaims(initial.claims || []);
  renderArtifacts(initial.artifacts || []);
  renderSync(initial.agents || []);
  renderApprovals(initial.approvals || []);
  scrollToBottom();
  autoresize(bodyEl);

  poll(true);
})();
"""


def main() -> None:
    _init_schema()
    bind = os.environ.get("COLLAB_ROOM_BIND", "127.0.0.1:7879")
    if ":" in bind:
        host, port_str = bind.rsplit(":", 1)
        port = int(port_str)
    else:
        host, port = bind, 7879
    if not ROOM_TOKEN:
        sys.stderr.write("warning: COLLAB_ROOM_TOKEN/COLLAB_MCP_TOKEN unset; room has no login gate\n")
    httpd = ThreadingHTTPServer((host, port), RoomHandler)
    sys.stderr.write(
        f"acs-collab-room listening on http://{host}:{port} "
        f"(db={DB_PATH}, ledger={LEDGER_PATH or 'disabled'}, author={PI_AUTHOR})\n"
    )
    httpd.serve_forever()


if __name__ == "__main__":
    main()
