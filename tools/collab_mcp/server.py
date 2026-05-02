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
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS messages(
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                ts TEXT NOT NULL,
                author TEXT NOT NULL,
                addressee TEXT NOT NULL,
                topic TEXT NOT NULL,
                body TEXT NOT NULL,
                status TEXT NOT NULL,
                refs TEXT NOT NULL DEFAULT '[]'
            );
            CREATE INDEX IF NOT EXISTS idx_messages_id ON messages(id);

            CREATE TABLE IF NOT EXISTS claims(
                topic TEXT PRIMARY KEY,
                author TEXT NOT NULL,
                claimed_at TEXT NOT NULL,
                released_at TEXT,
                summary TEXT
            );

            CREATE TABLE IF NOT EXISTS cursors(
                author TEXT PRIMARY KEY,
                last_seen_id INTEGER NOT NULL DEFAULT 0
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
) -> dict[str, Any]:
    """Append a message to the shared queue and (if configured) the ledger file.

    Rate-limit guard: a single author cannot send four consecutive messages
    without the other side replying — protects against runaway loops.
    """
    author = _auth(ctx)
    to = _require_author(to, "to")
    topic = _require_nonempty(topic, "topic")
    body = _require_nonempty(body, "body")
    if status not in ALLOWED_STATUSES:
        raise ValueError(
            f"status must be one of {sorted(ALLOWED_STATUSES)} (got {status!r})"
        )
    refs = list(refs or [])
    ts_human = _now_human()
    with closing(_db()) as conn:
        recent = conn.execute(
            "SELECT author FROM messages ORDER BY id DESC LIMIT 3"
        ).fetchall()
        if len(recent) == 3 and all(r[0] == author for r in recent):
            raise RuntimeError(
                f"rate limit: {author} sent 3 in a row; wait for {to} to reply"
            )
        cur = conn.execute(
            "INSERT INTO messages(ts, author, addressee, topic, body, status, refs) "
            "VALUES (?,?,?,?,?,?,?)",
            (ts_human, author, to, topic, body, status, json.dumps(refs)),
        )
        msg_id = cur.lastrowid
    appended = _ledger_append(ts_human, author, to, topic, body, status, refs)
    return {
        "id": msg_id,
        "ts": ts_human,
        "author": author,
        "ledger_appended": appended,
    }


@mcp.tool()
def read(
    ctx: Context,
    unread_only: bool = True,
    limit: int = 20,
) -> dict[str, Any]:
    """Read messages. Default returns unread (advances cursor)."""
    author = _auth(ctx)
    if limit < 1 or limit > 200:
        raise ValueError("limit must be in [1, 200]")
    with closing(_db()) as conn:
        cur = conn.execute(
            "SELECT last_seen_id FROM cursors WHERE author=?", (author,)
        ).fetchone()
        last = cur[0] if cur else 0
        if unread_only:
            rows = conn.execute(
                "SELECT id, ts, author, addressee, topic, body, status, refs "
                "FROM messages WHERE id>? AND author<>? ORDER BY id ASC LIMIT ?",
                (last, author, limit),
            ).fetchall()
            if rows:
                conn.execute(
                    "INSERT INTO cursors(author, last_seen_id) VALUES(?,?) "
                    "ON CONFLICT(author) DO UPDATE SET last_seen_id=excluded.last_seen_id",
                    (author, rows[-1][0]),
                )
        else:
            rows = conn.execute(
                "SELECT id, ts, author, addressee, topic, body, status, refs "
                "FROM messages ORDER BY id DESC LIMIT ?",
                (limit,),
            ).fetchall()
            rows = list(reversed(rows))
    return {
        "count": len(rows),
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
            }
            for r in rows
        ],
    }


@mcp.tool()
def claim(ctx: Context, topic: str, intent: str = "") -> dict[str, Any]:
    """Atomically claim a work topic. Fails on conflict with another active claim."""
    author = _auth(ctx)
    topic = _require_nonempty(topic, "topic")
    now = _now_iso()
    with closing(_db()) as conn:
        conn.execute("BEGIN IMMEDIATE")
        try:
            existing = conn.execute(
                "SELECT author, claimed_at FROM claims "
                "WHERE topic=? AND released_at IS NULL",
                (topic,),
            ).fetchone()
            if existing and existing[0] != author:
                conn.execute("ROLLBACK")
                return {
                    "ok": False,
                    "conflict": True,
                    "owner": existing[0],
                    "claimed_at": existing[1],
                }
            conn.execute(
                "INSERT INTO claims(topic, author, claimed_at, summary) "
                "VALUES(?,?,?,?) "
                "ON CONFLICT(topic) DO UPDATE SET "
                "author=excluded.author, claimed_at=excluded.claimed_at, "
                "released_at=NULL, summary=excluded.summary",
                (topic, author, now, intent),
            )
            conn.execute("COMMIT")
        except Exception:
            conn.execute("ROLLBACK")
            raise
    return {"ok": True, "topic": topic, "owner": author, "claimed_at": now}


@mcp.tool()
def release(ctx: Context, topic: str, summary: str) -> dict[str, Any]:
    """Release a topic claim with a summary of what was done."""
    author = _auth(ctx)
    topic = _require_nonempty(topic, "topic")
    summary = _require_nonempty(summary, "summary")
    now = _now_iso()
    with closing(_db()) as conn:
        conn.execute("BEGIN IMMEDIATE")
        try:
            existing = conn.execute(
                "SELECT author FROM claims WHERE topic=? AND released_at IS NULL",
                (topic,),
            ).fetchone()
            if not existing:
                conn.execute("ROLLBACK")
                return {"ok": False, "error": "no active claim"}
            if existing[0] != author:
                conn.execute("ROLLBACK")
                return {"ok": False, "error": f"claim owned by {existing[0]}"}
            conn.execute(
                "UPDATE claims SET released_at=?, summary=? WHERE topic=?",
                (now, summary, topic),
            )
            conn.execute("COMMIT")
        except Exception:
            conn.execute("ROLLBACK")
            raise
    return {"ok": True, "topic": topic, "released_at": now, "summary": summary}


@mcp.tool()
def status(ctx: Context) -> dict[str, Any]:
    """Active claims plus per-author unread counts."""
    _auth(ctx)
    with closing(_db()) as conn:
        active = conn.execute(
            "SELECT topic, author, claimed_at, summary FROM claims "
            "WHERE released_at IS NULL ORDER BY claimed_at"
        ).fetchall()
        unread: dict[str, int] = {}
        for a in sorted(ALLOWED_AUTHORS):
            cur = conn.execute(
                "SELECT last_seen_id FROM cursors WHERE author=?", (a,)
            ).fetchone()
            last = cur[0] if cur else 0
            cnt = conn.execute(
                "SELECT COUNT(*) FROM messages WHERE id>? AND author<>?",
                (last, a),
            ).fetchone()[0]
            unread[a] = cnt
        last_msg = conn.execute(
            "SELECT id, ts, author, topic FROM messages ORDER BY id DESC LIMIT 1"
        ).fetchone()
    return {
        "active_claims": [
            {"topic": r[0], "owner": r[1], "claimed_at": r[2], "intent": r[3]}
            for r in active
        ],
        "unread": unread,
        "last_message": (
            {"id": last_msg[0], "ts": last_msg[1], "from": last_msg[2], "topic": last_msg[3]}
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
) -> dict[str, Any]:
    """Read recent messages without advancing the caller's cursor.

    Used to warm up a fresh CLI/desktop session with prior collab context.
    Unlike `read`, this includes the caller's own messages so a returning
    session sees the full timeline. Filters: optional `topic` exact match,
    optional `since_id` to fetch only messages newer than a known id.
    """
    _auth(ctx)
    if limit < 1 or limit > 200:
        raise ValueError("limit must be in [1, 200]")
    if since_id < 0:
        raise ValueError("since_id must be >= 0")
    sql = (
        "SELECT id, ts, author, addressee, topic, body, status, refs "
        "FROM messages WHERE id>?"
    )
    params: list[Any] = [since_id]
    if topic:
        sql += " AND topic=?"
        params.append(topic)
    sql += " ORDER BY id DESC LIMIT ?"
    params.append(limit)
    with closing(_db()) as conn:
        rows = conn.execute(sql, params).fetchall()
    rows = list(reversed(rows))
    return {
        "count": len(rows),
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
