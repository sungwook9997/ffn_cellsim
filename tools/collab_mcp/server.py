"""Claude ↔ Codex collaboration MCP server.

Single source of truth for the Claude/Codex async channel. Two-LLM message
queue + topic locking + automatic ledger append. Hosts on the Windows
workstation; clients are Mac Claude Code and Win Codex CLI, both reaching
the server over Tailscale.

Configuration (all via env vars, never committed):
  COLLAB_MCP_BIND      bind address (default ``127.0.0.1:7878``)
                       localhost-only is the safe default. Set
                       ``0.0.0.0:<port>`` only when you intentionally expose
                       the server over Tailscale.
  COLLAB_MCP_TOKEN     shared bearer token; required, no default
  COLLAB_MCP_DB        SQLite path (default ``~/.acs-collab/inbox.db``)
  COLLAB_MCP_LEDGER    optional path to ``docs/claude_codex_log.md`` for
                       automatic append on every send. If unset, the ledger
                       mirror is disabled (DB remains the source of truth).
  COLLAB_MCP_AUTHORS   comma-separated allowed authors
                       (default ``claude,codex``)

Per-request authentication / identity:
  Authorization: Bearer <COLLAB_MCP_TOKEN>
  X-Collab-Author: claude | codex

The server intentionally exposes no command-execution surface. Tools are
limited to message append, message read, topic claim/release, status query,
and artifact-path announcement. GPU job launches must go through a
separately-approved wrapper added in a later iteration.
"""

from __future__ import annotations

import datetime as dt
import hmac
import json
import os
import pathlib
import re
import sqlite3
import sys
from contextlib import closing
from typing import Any

try:
    from mcp.server.fastmcp import FastMCP
    from mcp.server.fastmcp.server import Context
except ImportError as exc:  # pragma: no cover - import-time only
    sys.stderr.write(
        "fatal: 'mcp' package not installed. Run "
        "`pip install 'mcp[cli]>=1.2'` in the env that hosts the server.\n"
    )
    raise

ALLOWED_AUTHORS = {
    a.strip().lower()
    for a in os.environ.get("COLLAB_MCP_AUTHORS", "claude,codex").split(",")
    if a.strip()
}
ALLOWED_STATUSES = {
    "FYI",
    "proposal",
    "review",
    "decision-needed",
    "open-question",
    "blocker",
}
TOKEN = os.environ.get("COLLAB_MCP_TOKEN", "")
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
ROOM_NAME_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,31}$")
ALL_ROOMS_SENTINEL = "*"
_DEFAULT_ROOM_RAW = os.environ.get("COLLAB_DEFAULT_ROOM", "design-discussion").strip().lower() or "design-discussion"
if not ROOM_NAME_RE.match(_DEFAULT_ROOM_RAW):
    sys.stderr.write(
        f"fatal: COLLAB_DEFAULT_ROOM={_DEFAULT_ROOM_RAW!r} does not match "
        f"{ROOM_NAME_RE.pattern}. The value is interpolated into schema "
        "DDL defaults and must be a safe room identifier.\n"
    )
    sys.exit(2)
DEFAULT_ROOM = _DEFAULT_ROOM_RAW

if not TOKEN:
    sys.stderr.write(
        "fatal: COLLAB_MCP_TOKEN is unset. Set a shared secret on both server "
        "and clients before starting the server.\n"
    )
    sys.exit(2)

DB_PATH.parent.mkdir(parents=True, exist_ok=True)


def _db() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def _init_schema() -> None:
    with closing(_db()) as conn:
        # Create-from-scratch path. Existing databases skip the CREATE
        # IF NOT EXISTS lines and fall through to the migration block
        # below.
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
            """
        )

        # Idempotent migration for pre-room databases. SQLite ADD COLUMN
        # backfills existing rows with the DEFAULT, so legacy messages
        # are bucketed into design-discussion (the room where every
        # pre-Task-2 conversation lived).
        existing_cols = {row[1] for row in conn.execute("PRAGMA table_info(messages)").fetchall()}
        if "room" not in existing_cols:
            conn.execute(
                f"ALTER TABLE messages ADD COLUMN room TEXT NOT NULL DEFAULT '{DEFAULT_ROOM}'"
            )

        # Migrate single-author cursor table to per-(author, room) cursors.
        # Pre-Task-2 schema was PRIMARY KEY(author) — that table can't be
        # altered into a composite PK in place, so we rebuild it. The
        # legacy row maps onto (author, DEFAULT_ROOM) since every
        # pre-existing message lived in that room post-messages-migration.
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

        # Now safe in either path: room column is guaranteed.
        conn.execute("CREATE INDEX IF NOT EXISTS idx_messages_room_id ON messages(room, id)")

        # Migrate the claims table to a composite (room, topic) PK so two
        # workrooms can claim the same logical topic name independently.
        # Pre-Task-2 schema was PRIMARY KEY(topic); rebuild and bucket
        # legacy claims into DEFAULT_ROOM.
        claims_cols = {row[1] for row in conn.execute("PRAGMA table_info(claims)").fetchall()}
        if "room" not in claims_cols:
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


_init_schema()


def _request_headers(ctx: Context) -> dict[str, str]:
    """Best-effort header extraction for streamable-http transport.

    FastMCP exposes the underlying ASGI request via the request context. The
    attribute path differs slightly between mcp versions, so we walk a few
    candidates and fall back to an empty mapping for stdio (dev) transport.
    """
    rc = getattr(ctx, "request_context", None)
    if rc is None:
        return {}
    request = getattr(rc, "request", None) or getattr(rc, "http_request", None)
    headers = getattr(request, "headers", None) if request is not None else None
    if headers is None:
        return {}
    try:
        return {k.lower(): v for k, v in headers.items()}
    except Exception:
        return {}


def _check_token(headers: dict[str, str]) -> None:
    raw = headers.get("authorization", "")
    presented = raw[7:].strip() if raw.lower().startswith("bearer ") else raw.strip()
    if not hmac.compare_digest(presented, TOKEN):
        raise PermissionError("invalid or missing bearer token")


def _author(headers: dict[str, str]) -> str:
    a = headers.get("x-collab-author", "").strip().lower()
    if a not in ALLOWED_AUTHORS:
        raise ValueError(
            f"x-collab-author header must be one of {sorted(ALLOWED_AUTHORS)} (got {a!r})"
        )
    return a


def _auth(ctx: Context) -> str:
    headers = _request_headers(ctx)
    _check_token(headers)
    return _author(headers)


def _caller_room(ctx: Context) -> str:
    """Return the workroom name implied by the caller's headers / env.

    Resolution order:
      1. ``X-Collab-Room`` header on the request (per-CLI MCP client config).
      2. Server env ``COLLAB_DEFAULT_ROOM``.
      3. Hard fallback ``design-discussion``.

    The raw header value is validated against ``ROOM_NAME_RE`` (strict
    lowercase) before being accepted — silent case-folding would mask a
    misconfigured CLI shipping an uppercase room name, which is exactly
    the bug the regex is meant to surface.
    """
    headers = _request_headers(ctx)
    raw = headers.get("x-collab-room", "").strip()
    if raw:
        if not ROOM_NAME_RE.match(raw):
            raise ValueError(
                f"X-Collab-Room header must match {ROOM_NAME_RE.pattern} (got {raw!r})"
            )
        return raw
    return DEFAULT_ROOM


def _resolve_room(value: str | None, default: str) -> str:
    """Validate an explicit ``room`` arg or fall back to the caller's room.

    Like ``_caller_room``, validation runs against the raw stripped value
    so a typo with the wrong case (``IMPLEMENTATION-WORK``) is rejected
    rather than silently normalized.
    """
    if value is None:
        return default
    cleaned = value.strip()
    if not cleaned:
        return default
    if cleaned == ALL_ROOMS_SENTINEL:
        return ALL_ROOMS_SENTINEL
    if not ROOM_NAME_RE.match(cleaned):
        raise ValueError(
            f"room must match {ROOM_NAME_RE.pattern} or be '{ALL_ROOMS_SENTINEL}' (got {value!r})"
        )
    return cleaned


def _require_nonempty(value: str, field_name: str) -> str:
    cleaned = value.strip()
    if not cleaned:
        raise ValueError(f"{field_name} must be non-empty")
    return cleaned


def _require_author(value: str, field_name: str) -> str:
    cleaned = _require_nonempty(value, field_name).lower()
    if cleaned not in ALLOWED_AUTHORS:
        raise ValueError(
            f"{field_name} must be one of {sorted(ALLOWED_AUTHORS)} (got {value!r})"
        )
    return cleaned


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


def _now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def _now_human() -> str:
    return dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M (%Z)")


mcp = FastMCP("acs-collab")


@mcp.tool()
def send(
    ctx: Context,
    to: str,
    topic: str,
    body: str,
    status: str = "FYI",
    refs: list[str] | None = None,
    room: str | None = None,
) -> dict[str, Any]:
    """Append a message to the shared queue and (if configured) the ledger file.

    The message is tagged with a ``room`` (workroom) — when omitted, the
    caller's default room is used (see ``_caller_room``). The rate-limit
    guard is scoped per-room so two workrooms running in parallel don't
    starve each other.
    """
    author = _auth(ctx)
    caller_room = _caller_room(ctx)
    to = _require_author(to, "to")
    topic = _require_nonempty(topic, "topic")
    body = _require_nonempty(body, "body")
    if status not in ALLOWED_STATUSES:
        raise ValueError(
            f"status must be one of {sorted(ALLOWED_STATUSES)} (got {status!r})"
        )
    target_room = _resolve_room(room, caller_room)
    if target_room == ALL_ROOMS_SENTINEL:
        raise ValueError("send: room must be a single room name, not '*'")
    refs = list(refs or [])
    ts_human = _now_human()
    with closing(_db()) as conn:
        recent = conn.execute(
            "SELECT author FROM messages WHERE room=? ORDER BY id DESC LIMIT 3",
            (target_room,),
        ).fetchall()
        if len(recent) == 3 and all(r[0] == author for r in recent):
            raise RuntimeError(
                f"rate limit: {author} sent 3 in a row in room {target_room!r}; "
                f"wait for {to} to reply"
            )
        cur = conn.execute(
            "INSERT INTO messages(ts, author, addressee, topic, body, status, refs, room) "
            "VALUES (?,?,?,?,?,?,?,?)",
            (ts_human, author, to, topic, body, status, json.dumps(refs), target_room),
        )
        msg_id = cur.lastrowid
    appended = _ledger_append(ts_human, author, to, topic, body, status, refs)
    return {
        "id": msg_id,
        "ts": ts_human,
        "author": author,
        "room": target_room,
        "ledger_appended": appended,
    }


@mcp.tool()
def read(
    ctx: Context,
    unread_only: bool = True,
    limit: int = 20,
    room: str | None = None,
) -> dict[str, Any]:
    """Read messages. Default returns unread (advances cursor).

    Filters to the caller's current room by default. Pass an explicit
    ``room`` name to read another room, or ``room='*'`` for a cross-room
    view (admin / migration use).
    """
    author = _auth(ctx)
    caller_room = _caller_room(ctx)
    target_room = _resolve_room(room, caller_room)
    if limit < 1 or limit > 200:
        raise ValueError("limit must be in [1, 200]")
    if target_room == ALL_ROOMS_SENTINEL and unread_only:
        # Cursors are per-(author, room). "Unread" across rooms is
        # ambiguous: each room has its own cursor and skipping it would
        # leak messages. Force the caller to either name a room or read
        # the recent timeline (unread_only=False).
        raise ValueError(
            "read: room='*' requires unread_only=False (cursors are per-room)"
        )
    where_room = "" if target_room == ALL_ROOMS_SENTINEL else " AND room=?"
    room_param: tuple[Any, ...] = () if target_room == ALL_ROOMS_SENTINEL else (target_room,)
    with closing(_db()) as conn:
        if unread_only:
            cur = conn.execute(
                "SELECT last_seen_id FROM cursors WHERE author=? AND room=?",
                (author, target_room),
            ).fetchone()
            last = cur[0] if cur else 0
            sql = (
                "SELECT id, ts, author, addressee, topic, body, status, refs, room "
                "FROM messages WHERE id>? AND author<>?" + where_room +
                " ORDER BY id ASC LIMIT ?"
            )
            rows = conn.execute(sql, (last, author, *room_param, limit)).fetchall()
            if rows:
                conn.execute(
                    "INSERT INTO cursors(author, room, last_seen_id) VALUES(?,?,?) "
                    "ON CONFLICT(author, room) DO UPDATE SET last_seen_id=excluded.last_seen_id",
                    (author, target_room, rows[-1][0]),
                )
        else:
            sql = (
                "SELECT id, ts, author, addressee, topic, body, status, refs, room "
                "FROM messages WHERE 1=1" + where_room +
                " ORDER BY id DESC LIMIT ?"
            )
            rows = conn.execute(sql, (*room_param, limit)).fetchall()
            rows = list(reversed(rows))
    return {
        "count": len(rows),
        "room": target_room,
        "messages": [
            {
                "id": r[0],
                "ts": r[1],
                "from": r[2],
                "to": r[3],
                "topic": r[4],
                "body": r[5],
                "status": r[6],
                "refs": json.loads(r[7]),
                "room": r[8],
            }
            for r in rows
        ],
    }


@mcp.tool()
def claim(
    ctx: Context,
    topic: str,
    intent: str = "",
    room: str | None = None,
) -> dict[str, Any]:
    """Atomically claim a work topic in a workroom.

    Claims are scoped per ``(room, topic)``: two workrooms can hold the
    same logical topic name without conflict. ``room`` defaults to the
    caller's workroom (header / env / fallback).
    """
    author = _auth(ctx)
    target_room = _resolve_room(room, _caller_room(ctx))
    if target_room == ALL_ROOMS_SENTINEL:
        raise ValueError("claim: room must be a single room name, not '*'")
    topic = _require_nonempty(topic, "topic")
    now = _now_iso()
    with closing(_db()) as conn:
        conn.execute("BEGIN IMMEDIATE")
        try:
            existing = conn.execute(
                "SELECT author, claimed_at FROM claims "
                "WHERE room=? AND topic=? AND released_at IS NULL",
                (target_room, topic),
            ).fetchone()
            if existing and existing[0] != author:
                conn.execute("ROLLBACK")
                return {
                    "ok": False,
                    "conflict": True,
                    "owner": existing[0],
                    "claimed_at": existing[1],
                    "room": target_room,
                }
            conn.execute(
                "INSERT INTO claims(room, topic, author, claimed_at, summary) "
                "VALUES(?,?,?,?,?) "
                "ON CONFLICT(room, topic) DO UPDATE SET "
                "author=excluded.author, claimed_at=excluded.claimed_at, "
                "released_at=NULL, summary=excluded.summary",
                (target_room, topic, author, now, intent),
            )
            conn.execute("COMMIT")
        except Exception:
            conn.execute("ROLLBACK")
            raise
    return {"ok": True, "topic": topic, "room": target_room, "owner": author, "claimed_at": now}


@mcp.tool()
def release(
    ctx: Context,
    topic: str,
    summary: str,
    room: str | None = None,
) -> dict[str, Any]:
    """Release a topic claim with a summary of what was done.

    ``room`` defaults to the caller's workroom; release only matches a
    claim in the same room.
    """
    author = _auth(ctx)
    target_room = _resolve_room(room, _caller_room(ctx))
    if target_room == ALL_ROOMS_SENTINEL:
        raise ValueError("release: room must be a single room name, not '*'")
    topic = _require_nonempty(topic, "topic")
    summary = _require_nonempty(summary, "summary")
    now = _now_iso()
    with closing(_db()) as conn:
        conn.execute("BEGIN IMMEDIATE")
        try:
            existing = conn.execute(
                "SELECT author FROM claims "
                "WHERE room=? AND topic=? AND released_at IS NULL",
                (target_room, topic),
            ).fetchone()
            if not existing:
                conn.execute("ROLLBACK")
                return {"ok": False, "error": "no active claim", "room": target_room}
            if existing[0] != author:
                conn.execute("ROLLBACK")
                return {"ok": False, "error": f"claim owned by {existing[0]}", "room": target_room}
            conn.execute(
                "UPDATE claims SET released_at=?, summary=? WHERE room=? AND topic=?",
                (now, summary, target_room, topic),
            )
            conn.execute("COMMIT")
        except Exception:
            conn.execute("ROLLBACK")
            raise
    return {"ok": True, "topic": topic, "room": target_room, "released_at": now, "summary": summary}


@mcp.tool()
def status(ctx: Context) -> dict[str, Any]:
    """Active claims plus per-author unread counts.

    Unread counts are reported both for the caller's current room and for
    every distinct room in the DB, so a workroom-aware client can show its
    own number prominently and still see whether another room has traffic.
    """
    _auth(ctx)
    caller_room = _caller_room(ctx)
    with closing(_db()) as conn:
        # Active claims for the caller's room go in `active_claims`;
        # claims_by_room shows the cross-room view for callers wiring up
        # the sidebar / multi-room dashboards.
        active = conn.execute(
            "SELECT topic, author, claimed_at, summary FROM claims "
            "WHERE released_at IS NULL AND room=? ORDER BY claimed_at",
            (caller_room,),
        ).fetchall()
        cross_active = conn.execute(
            "SELECT room, topic, author, claimed_at, summary FROM claims "
            "WHERE released_at IS NULL ORDER BY room, claimed_at"
        ).fetchall()
        all_rooms = sorted(
            {row[0] for row in conn.execute("SELECT DISTINCT room FROM messages").fetchall()}
            | {row[0] for row in cross_active}
            | {caller_room}
        )
        unread: dict[str, int] = {}
        unread_by_room: dict[str, dict[str, int]] = {}
        for a in sorted(ALLOWED_AUTHORS):
            unread_by_room[a] = {}
            for r in all_rooms:
                cur = conn.execute(
                    "SELECT last_seen_id FROM cursors WHERE author=? AND room=?",
                    (a, r),
                ).fetchone()
                last = cur[0] if cur else 0
                rcnt = conn.execute(
                    "SELECT COUNT(*) FROM messages WHERE id>? AND author<>? AND room=?",
                    (last, a, r),
                ).fetchone()[0]
                unread_by_room[a][r] = rcnt
            unread[a] = unread_by_room[a].get(caller_room, 0)
        last_msg = conn.execute(
            "SELECT id, ts, author, topic, room FROM messages "
            "WHERE room=? ORDER BY id DESC LIMIT 1",
            (caller_room,),
        ).fetchone()
    claims_by_room: dict[str, list[dict[str, Any]]] = {r: [] for r in all_rooms}
    for row in cross_active:
        claims_by_room.setdefault(row[0], []).append(
            {"topic": row[1], "owner": row[2], "claimed_at": row[3], "intent": row[4]}
        )
    return {
        "room": caller_room,
        "rooms": all_rooms,
        "active_claims": [
            {"topic": r[0], "owner": r[1], "claimed_at": r[2], "intent": r[3]}
            for r in active
        ],
        "claims_by_room": claims_by_room,
        "unread": unread,
        "unread_by_room": unread_by_room,
        "last_message": (
            {"id": last_msg[0], "ts": last_msg[1], "from": last_msg[2],
             "topic": last_msg[3], "room": last_msg[4]}
            if last_msg
            else None
        ),
    }


@mcp.tool()
def announce_artifact(
    ctx: Context,
    kind: str,
    path: str,
    note: str = "",
) -> dict[str, Any]:
    """Record that an artifact (figure, HDF5, manifest, dashboard) is available.

    Stores a pointer only; never reads the path. The receiver is expected to
    open the path locally (it's the same Syncthing-shared tree on both sides).
    """
    author = _auth(ctx)
    kind = _require_nonempty(kind, "kind")
    path = _require_nonempty(path, "path")
    now = _now_human()
    with closing(_db()) as conn:
        cur = conn.execute(
            "INSERT INTO artifacts(ts, author, kind, path, note) VALUES (?,?,?,?,?)",
            (now, author, kind, path, note),
        )
        artifact_id = cur.lastrowid
    return {"id": artifact_id, "ts": now, "author": author, "kind": kind, "path": path}


@mcp.tool()
def list_artifacts(ctx: Context, limit: int = 20) -> dict[str, Any]:
    """List the most recent artifact pointers."""
    _auth(ctx)
    if limit < 1 or limit > 200:
        raise ValueError("limit must be in [1, 200]")
    with closing(_db()) as conn:
        rows = conn.execute(
            "SELECT id, ts, author, kind, path, note FROM artifacts "
            "ORDER BY id DESC LIMIT ?",
            (limit,),
        ).fetchall()
    return {
        "count": len(rows),
        "artifacts": [
            {
                "id": r[0],
                "ts": r[1],
                "author": r[2],
                "kind": r[3],
                "path": r[4],
                "note": r[5],
            }
            for r in rows
        ],
    }


@mcp.tool()
def bootstrap(
    ctx: Context,
    limit: int = 50,
    topic: str = "",
    since_id: int = 0,
    advance_cursor: bool = True,
    room: str | None = None,
) -> dict[str, Any]:
    """Read recent messages and mark the caller current.

    Used to warm up a fresh CLI/desktop session with prior collab context.
    Unlike `read`, includes the caller's own messages so a returning session
    sees the full timeline. By default advances the caller's read-cursor to
    the latest returned id so the UI's "seen by X" indicator catches up
    immediately — set ``advance_cursor=False`` for a pure read-only peek.

    Filters: optional ``topic`` exact match, optional ``since_id`` to fetch
    only messages newer than a known id, optional ``room`` to switch
    workrooms (default = caller's current room; ``'*'`` = all rooms).
    """
    author = _auth(ctx)
    caller_room = _caller_room(ctx)
    target_room = _resolve_room(room, caller_room)
    if limit < 1 or limit > 200:
        raise ValueError("limit must be in [1, 200]")
    if since_id < 0:
        raise ValueError("since_id must be >= 0")
    sql = (
        "SELECT id, ts, author, addressee, topic, body, status, refs, room "
        "FROM messages WHERE id>?"
    )
    params: list[Any] = [since_id]
    if topic:
        sql += " AND topic=?"
        params.append(topic)
    if target_room != ALL_ROOMS_SENTINEL:
        sql += " AND room=?"
        params.append(target_room)
    sql += " ORDER BY id DESC LIMIT ?"
    params.append(limit)
    cursor_advanced_to: int | dict[str, int] | None = None
    with closing(_db()) as conn:
        rows = conn.execute(sql, params).fetchall()
        rows = list(reversed(rows))
        if advance_cursor and rows:
            # Cursors are per-(author, room). For a single-room bootstrap
            # the answer is the max id seen. For room='*' we walk the
            # rooms touched by the result and advance each cursor to the
            # highest id from that room — preserving the per-room
            # invariant that future unread reads do not skip messages
            # from a different room that happen to have higher ids.
            per_room_max: dict[str, int] = {}
            for row in rows:
                rm = row[8] or DEFAULT_ROOM
                if row[0] > per_room_max.get(rm, 0):
                    per_room_max[rm] = row[0]
            for rm, max_id in per_room_max.items():
                conn.execute(
                    "INSERT INTO cursors(author, room, last_seen_id) VALUES(?,?,?) "
                    "ON CONFLICT(author, room) DO UPDATE SET "
                    "last_seen_id=MAX(cursors.last_seen_id, excluded.last_seen_id)",
                    (author, rm, max_id),
                )
            if target_room == ALL_ROOMS_SENTINEL:
                cursor_advanced_to = per_room_max
            else:
                cursor_advanced_to = per_room_max.get(target_room)
    return {
        "count": len(rows),
        "room": target_room,
        "cursor_advanced_to": cursor_advanced_to,
        "messages": [
            {
                "id": r[0],
                "ts": r[1],
                "from": r[2],
                "to": r[3],
                "topic": r[4],
                "body": r[5],
                "status": r[6],
                "refs": json.loads(r[7]),
                "room": r[8],
            }
            for r in rows
        ],
    }


def main() -> None:
    bind = os.environ.get("COLLAB_MCP_BIND", "127.0.0.1:7878")
    if ":" in bind:
        host, port_str = bind.rsplit(":", 1)
        port = int(port_str)
    else:
        host, port = bind, 7878
    mcp.settings.host = host
    mcp.settings.port = port
    sys.stderr.write(
        f"acs-collab-mcp listening on {host}:{port} "
        f"(db={DB_PATH}, ledger={LEDGER_PATH or 'disabled'}, "
        f"authors={sorted(ALLOWED_AUTHORS)})\n"
    )
    mcp.run(transport="streamable-http")


if __name__ == "__main__":
    main()
