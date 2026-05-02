"""Browser room for the Claude/Codex/PI collaboration channel.

This is a tiny stdlib-only web UI over the same SQLite database used by
``server.py``. It intentionally exposes no command execution and no file reads:
the browser can inspect recent messages, active claims, artifact pointers, and
post PI messages into the shared queue.
"""

from __future__ import annotations

import datetime as dt
import html
import http.cookies
import hmac
import json
import os
import pathlib
import sqlite3
import sys
import urllib.parse
from contextlib import closing
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any

ALLOWED_STATUSES = {
    "FYI",
    "proposal",
    "review",
    "decision-needed",
    "open-question",
    "blocker",
}
ADDRESSEE_PRESETS = ("claude,codex", "claude", "codex")
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
PI_AUTHOR = os.environ.get("COLLAB_ROOM_AUTHOR", "pi").strip().lower() or "pi"

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

            CREATE TABLE IF NOT EXISTS agent_status(
                agent TEXT PRIMARY KEY,
                percent INTEGER NOT NULL DEFAULT 0,
                activity TEXT NOT NULL DEFAULT '',
                topic TEXT NOT NULL DEFAULT '',
                updated_at TEXT NOT NULL
            );
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
) -> dict[str, Any]:
    topic = _require_nonempty(topic, "topic")
    body = _require_nonempty(body, "body")
    addressee = _require_nonempty(addressee, "to")
    if status not in ALLOWED_STATUSES:
        raise ValueError(f"status must be one of {sorted(ALLOWED_STATUSES)}")

    ts = _now_human()
    with closing(_db()) as conn:
        cur = conn.execute(
            "INSERT INTO messages(ts, author, addressee, topic, body, status, refs) "
            "VALUES (?,?,?,?,?,?,?)",
            (ts, PI_AUTHOR, addressee, topic, body, status, json.dumps(refs)),
        )
        msg_id = cur.lastrowid
    ledger_appended = _ledger_append(ts, PI_AUTHOR, addressee, topic, body, status, refs)
    return {"id": msg_id, "ts": ts, "ledger_appended": ledger_appended}


def _fetch_messages(after_id: int = 0, limit: int = 200) -> list[tuple[Any, ...]]:
    with closing(_db()) as conn:
        return conn.execute(
            "SELECT id, ts, author, addressee, topic, body, status, refs "
            "FROM messages WHERE id > ? ORDER BY id ASC LIMIT ?",
            (int(after_id), int(limit)),
        ).fetchall()


def _fetch_recent_messages(limit: int = 200) -> list[tuple[Any, ...]]:
    with closing(_db()) as conn:
        rows = conn.execute(
            "SELECT id, ts, author, addressee, topic, body, status, refs "
            "FROM messages ORDER BY id DESC LIMIT ?",
            (int(limit),),
        ).fetchall()
    return list(reversed(rows))


def _fetch_sidebar() -> dict[str, Any]:
    with closing(_db()) as conn:
        claims = conn.execute(
            "SELECT topic, author, claimed_at, summary FROM claims "
            "WHERE released_at IS NULL ORDER BY claimed_at"
        ).fetchall()
        artifacts = conn.execute(
            "SELECT id, ts, author, kind, path, note FROM artifacts "
            "ORDER BY id DESC LIMIT 30"
        ).fetchall()
        cursors = conn.execute(
            "SELECT author, last_seen_id FROM cursors ORDER BY author"
        ).fetchall()
        agents = conn.execute(
            "SELECT agent, percent, activity, topic, updated_at FROM agent_status "
            "ORDER BY agent"
        ).fetchall()
    return {"claims": claims, "artifacts": artifacts, "cursors": cursors, "agents": agents}


def _agents_to_dict(agents: list[tuple[Any, ...]]) -> list[dict[str, Any]]:
    out = []
    now = dt.datetime.now().astimezone()
    for agent, percent, activity, topic, updated_at in agents:
        age_seconds = None
        parsed = None
        try:
            parsed = dt.datetime.strptime(updated_at, "%Y-%m-%dT%H:%M:%S%z")
        except (ValueError, TypeError):
            try:
                parsed = dt.datetime.strptime(updated_at[:16], "%Y-%m-%d %H:%M")
                parsed = parsed.replace(tzinfo=now.tzinfo)
            except (ValueError, TypeError):
                parsed = None
        if parsed is not None:
            age_seconds = max(0, int((now - parsed).total_seconds()))
        out.append(
            {
                "agent": agent,
                "percent": int(percent),
                "activity": activity or "",
                "topic": topic or "",
                "updated_at": updated_at,
                "age_seconds": age_seconds,
            }
        )
    return out


_AGENT_VALUES = {
    "claude", "codex", "pi",
    "claude-chat", "claude-work",
    "codex-chat", "codex-work",
}


def _agent_base(agent: str) -> str:
    """Map a pane-aware agent value to its base agent (claude/codex/pi)."""
    if agent.startswith("claude"):
        return "claude"
    if agent.startswith("codex"):
        return "codex"
    return agent


def _upsert_agent_status(
    agent: str, percent: int, activity: str, topic: str = ""
) -> dict[str, Any]:
    agent = _require_nonempty(agent, "agent").lower()
    if agent not in _AGENT_VALUES:
        raise ValueError(
            "agent must be one of "
            "claude/codex/pi/claude-chat/claude-work/codex-chat/codex-work"
        )
    try:
        percent_i = int(percent)
    except (TypeError, ValueError) as exc:
        raise ValueError("percent must be an integer 0..100") from exc
    if percent_i < 0 or percent_i > 100:
        raise ValueError("percent must be in [0, 100]")
    activity = (activity or "").strip()[:240]
    topic = (topic or "").strip()[:120]
    ts = _now_iso_seconds()
    with closing(_db()) as conn:
        conn.execute(
            "INSERT INTO agent_status(agent, percent, activity, topic, updated_at) "
            "VALUES(?,?,?,?,?) "
            "ON CONFLICT(agent) DO UPDATE SET "
            "percent=excluded.percent, activity=excluded.activity, "
            "topic=excluded.topic, updated_at=excluded.updated_at",
            (agent, percent_i, activity, topic, ts),
        )
    return {"agent": agent, "percent": percent_i, "activity": activity, "topic": topic, "updated_at": ts}


def _message_to_dict(row: tuple[Any, ...]) -> dict[str, Any]:
    msg_id, ts, author, addressee, topic, body, status, refs_json = row
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
    <div class="agent-strip" id="agent-strip" aria-live="polite"></div>
    <nav>
      <button type="button" id="toggle-sidebar" class="ghost">Sidebar</button>
      <a href="/logout" class="ghost">Logout</a>
    </nav>
  </header>
  <main class="layout">
    <section class="chat" aria-label="Message timeline">
      {flash_html}
      <div class="messages" id="messages" role="log" aria-live="polite"></div>
      <button type="button" id="jump-bottom" class="jump hidden" aria-label="Jump to newest">↓ <span id="jump-count">new</span></button>
      <form class="composer" id="composer" autocomplete="off">
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
        </div>
        <div class="composer-row">
          <textarea name="body" id="composer-body" rows="2"
            placeholder="메시지를 입력하세요. Enter=전송, Shift+Enter=줄바꿈, @claude / @codex / @both 멘션 가능."
            required></textarea>
          <button type="submit" id="composer-send">Send</button>
        </div>
      </form>
    </section>
    <aside class="sidebar" id="sidebar">
      <section class="panel">
        <h2>Agents · 4 panes</h2>
        <div class="agent-cards" id="agent-cards"></div>
      </section>
      <section class="panel">
        <h2>Active Claims</h2>
        <ul id="claims-list"><li class="empty">None</li></ul>
      </section>
      <section class="panel">
        <h2>Artifacts</h2>
        <ul id="artifacts-list"><li class="empty">None</li></ul>
      </section>
    </aside>
  </main>
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
            after = 0
            try:
                after = int(qs.get("after", ["0"])[0])
            except (TypeError, ValueError):
                after = 0
            payload = self._build_payload(after_id=after)
            self._send_json(payload)
            return
        if path != "/":
            self.send_error(HTTPStatus.NOT_FOUND, "not found")
            return
        if not self._is_authenticated():
            self._send_text(_render_login())
            return
        initial = self._build_payload(after_id=0, recent=True)
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

    def _build_payload(self, after_id: int, recent: bool = False) -> dict[str, Any]:
        if recent:
            rows = _fetch_recent_messages(limit=200)
        else:
            rows = _fetch_messages(after_id=after_id, limit=200)
        sidebar = _sidebar_to_dict(_fetch_sidebar())
        return {
            "messages": [_message_to_dict(r) for r in rows],
            "claims": sidebar["claims"],
            "artifacts": sidebar["artifacts"],
            "cursors": sidebar["cursors"],
            "agents": sidebar["agents"],
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
        if self._wants_json():
            payload = self._read_json()
            result = _insert_pi_message(
                topic=str(payload.get("topic", "")),
                body=str(payload.get("body", "")),
                addressee=str(payload.get("to", "")),
                status=str(payload.get("status", "open-question")),
                refs=[str(r).strip() for r in payload.get("refs", []) if str(r).strip()],
            )
            self._send_json({"ok": True, **result})
            return
        fields = self._read_form()
        result = _insert_pi_message(
            topic=fields.get("topic", [""])[0],
            body=fields.get("body", [""])[0],
            addressee=fields.get("to", [""])[0],
            status=fields.get("status", ["open-question"])[0],
            refs=_parse_refs(fields.get("refs", [""])[0]),
        )
        self._redirect(f"/?sent={result['id']}")

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
        )
        self._send_json({"ok": True, **result})

    def _read_form(self) -> dict[str, list[str]]:
        length = int(self.headers.get("Content-Length", "0"))
        if length > 65536:
            raise ValueError("request body too large")
        body = self.rfile.read(length).decode("utf-8")
        return urllib.parse.parse_qs(body, keep_blank_values=True)

    def _read_json(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        if length > 65536:
            raise ValueError("request body too large")
        raw = self.rfile.read(length).decode("utf-8") or "{}"
        try:
            data = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid json: {exc}") from exc
        if not isinstance(data, dict):
            raise ValueError("json body must be an object")
        return data

    def _wants_json(self) -> bool:
        ctype = (self.headers.get("Content-Type") or "").lower()
        accept = (self.headers.get("Accept") or "").lower()
        return ctype.startswith("application/json") or "application/json" in accept

    def _is_authenticated(self) -> bool:
        if not ROOM_TOKEN:
            return True
        cookies = http.cookies.SimpleCookie(self.headers.get("Cookie", ""))
        cookie = cookies.get("collab_room_token")
        return bool(cookie and self._token_ok(urllib.parse.unquote(cookie.value)))

    def _token_ok(self, token: str) -> bool:
        if not ROOM_TOKEN:
            return True
        return hmac.compare_digest(token, ROOM_TOKEN)

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
.cursor-strip { color: var(--muted); font-size: 12px; }
.cursor-strip .pill {
  background: rgba(0,0,0,0.05);
  border-radius: 10px;
  margin-right: 6px;
  padding: 2px 8px;
}
.agent-strip {
  align-items: center;
  display: flex;
  flex: 1 1 auto;
  flex-wrap: wrap;
  gap: 8px;
  justify-content: flex-end;
  min-width: 0;
}
.agent-pill {
  align-items: center;
  background: #f3f5f8;
  border: 1px solid var(--line);
  border-radius: 14px;
  display: inline-flex;
  font-size: 12px;
  gap: 8px;
  max-width: 360px;
  min-width: 0;
  padding: 4px 10px;
}
.agent-pill .name { font-weight: 700; }
.agent-pill .name.claude { color: var(--claude-tag); }
.agent-pill .name.codex { color: var(--codex-tag); }
.agent-pill .name.pi { color: var(--pi-tag); }
.agent-pill .pct {
  background: rgba(0,0,0,0.06);
  border-radius: 8px;
  font-variant-numeric: tabular-nums;
  padding: 1px 6px;
}
.agent-pill .activity {
  color: var(--muted);
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  min-width: 0;
}
.agent-pill .age { color: var(--muted); font-size: 10px; }
.agent-pill.stale { opacity: 0.55; }
.agent-pill.fresh { box-shadow: 0 0 0 2px rgba(36,107,254,0.18); }
nav { display: flex; gap: 8px; }
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
.layout {
  display: grid;
  grid-template-columns: minmax(0, 1fr) 280px;
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
.sidebar {
  background: var(--panel);
  border-left: 1px solid var(--line);
  display: flex;
  flex-direction: column;
  gap: 12px;
  overflow-y: auto;
  padding: 14px 14px 24px;
}
.sidebar.hidden { display: none; }
.agent-cards { display: flex; flex-direction: column; gap: 8px; }
.agent-card {
  background: #f6f8fb;
  border: 1px solid var(--line);
  border-radius: 8px;
  display: flex;
  flex-direction: column;
  font-size: 12px;
  gap: 4px;
  padding: 8px 10px;
  position: relative;
}
.agent-card.placeholder { background: transparent; border-style: dashed; }
.agent-card.fresh { box-shadow: 0 0 0 2px rgba(36,107,254,0.18); }
.agent-card.stale { opacity: 0.55; }
.agent-card .head {
  align-items: baseline;
  display: flex;
  gap: 6px;
  justify-content: space-between;
}
.agent-card .head .name { font-weight: 700; }
.agent-card .head .name.claude { color: var(--claude-tag); }
.agent-card .head .name.codex { color: var(--codex-tag); }
.agent-card .head .name.pi { color: var(--pi-tag); }
.agent-card .head .role {
  background: rgba(0,0,0,0.06);
  border-radius: 6px;
  color: var(--muted);
  font-size: 10px;
  padding: 1px 5px;
  text-transform: uppercase;
}
.agent-card .head .pct {
  font-variant-numeric: tabular-nums;
  font-weight: 700;
}
.agent-card .progress {
  background: rgba(0,0,0,0.07);
  border-radius: 4px;
  height: 6px;
  overflow: hidden;
}
.agent-card .progress > .bar {
  background: var(--claude-tag);
  height: 100%;
  transition: width 0.3s ease;
}
.agent-card.codex .progress > .bar { background: var(--codex-tag); }
.agent-card.pi .progress > .bar { background: var(--pi-tag); }
.agent-card .activity {
  color: var(--ink);
  white-space: pre-wrap;
  word-break: break-word;
}
.agent-card .meta {
  color: var(--muted);
  display: flex;
  font-size: 11px;
  gap: 8px;
  justify-content: space-between;
}
.panel { display: flex; flex-direction: column; }
.panel ul { list-style: none; margin: 0; padding: 0; display: flex; flex-direction: column; gap: 8px; }
.panel li {
  background: #f6f8fb;
  border: 1px solid var(--line);
  border-radius: 8px;
  font-size: 12px;
  padding: 8px 10px;
}
.panel li.empty { background: transparent; border-style: dashed; color: var(--muted); text-align: center; }
.panel li code { background: rgba(0,0,0,0.05); border-radius: 4px; font-size: 11px; padding: 1px 4px; }
.panel li .meta { color: var(--muted); display: block; font-size: 11px; margin-top: 4px; }
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
  .sidebar { border-left: 0; border-top: 1px solid var(--line); max-height: 240px; }
}
"""


JS = r"""
(function () {
  var initialEl = document.getElementById("initial-state");
  var initial = {messages: [], claims: [], artifacts: [], cursors: {}};
  if (initialEl) {
    try { initial = JSON.parse(initialEl.textContent || "{}"); } catch (e) {}
  }
  var state = {
    lastId: 0,
    messages: [],
    cursors: initial.cursors || {},
    agents: initial.agents || [],
    agentsFetchedAt: Date.now(),
    pending: 0,
  };
  var messagesEl = document.getElementById("messages");
  var cursorStrip = document.getElementById("cursor-strip");
  var agentStrip = document.getElementById("agent-strip");
  var claimsList = document.getElementById("claims-list");
  var artifactsList = document.getElementById("artifacts-list");
  var jumpBtn = document.getElementById("jump-bottom");
  var jumpCount = document.getElementById("jump-count");
  var sidebar = document.getElementById("sidebar");
  var toggleSidebar = document.getElementById("toggle-sidebar");
  var agentCardsEl = document.getElementById("agent-cards");
  var composer = document.getElementById("composer");
  var bodyEl = document.getElementById("composer-body");
  var sendBtn = document.getElementById("composer-send");

  var SELF_AUTHOR = "pi";

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

  function ageLabel(seconds) {
    if (seconds == null) return "";
    if (seconds < 60) return seconds + "s ago";
    if (seconds < 3600) return Math.floor(seconds / 60) + "m ago";
    if (seconds < 86400) return Math.floor(seconds / 3600) + "h ago";
    return Math.floor(seconds / 86400) + "d ago";
  }

  function renderAgents(agents) {
    if (agents) state.agents = agents;
    agents = state.agents;
    agentStrip.innerHTML = "";
    if (!agents || !agents.length) {
      var empty = document.createElement("span");
      empty.className = "agent-pill stale";
      empty.textContent = "no agent status yet";
      agentStrip.appendChild(empty);
      return;
    }
    var elapsedSinceFetch = Math.floor((Date.now() - state.agentsFetchedAt) / 1000);
    agents.forEach(function (a) {
      var age = (a.age_seconds == null) ? null : (a.age_seconds + elapsedSinceFetch);
      var pill = document.createElement("div");
      var staleness = "fresh";
      if (age == null || age > 300) staleness = "stale";
      else if (age > 60) staleness = "";
      pill.className = "agent-pill " + staleness;

      var name = document.createElement("span");
      name.className = "name " + (a.agent || "");
      name.textContent = a.agent;
      pill.appendChild(name);

      var pct = document.createElement("span");
      pct.className = "pct";
      pct.textContent = (a.percent || 0) + "%";
      pill.appendChild(pct);

      var act = document.createElement("span");
      act.className = "activity";
      var label = a.activity || "(idle)";
      if (a.topic) label = "[" + a.topic + "] " + label;
      act.textContent = label;
      act.title = label;
      pill.appendChild(act);

      var ageEl = document.createElement("span");
      ageEl.className = "age";
      ageEl.textContent = ageLabel(age);
      pill.appendChild(ageEl);

      agentStrip.appendChild(pill);
    });
  }

  var PANE_ORDER = ["claude-chat", "claude-work", "codex-chat", "codex-work"];

  function paneBase(name) {
    if (name && name.indexOf("claude") === 0) return "claude";
    if (name && name.indexOf("codex") === 0) return "codex";
    if (name === "pi") return "pi";
    return "other";
  }

  function paneRole(name) {
    if (!name) return "";
    var i = name.indexOf("-");
    return i > 0 ? name.slice(i + 1) : "";
  }

  function buildAgentCard(name, agent, elapsedSinceFetch) {
    var card = document.createElement("div");
    var base = paneBase(name);
    var role = paneRole(name);

    var age = null;
    if (agent && agent.age_seconds != null) {
      age = agent.age_seconds + elapsedSinceFetch;
    }
    var staleness = "placeholder";
    if (agent) {
      staleness = "fresh";
      if (age == null || age > 300) staleness = "stale";
      else if (age > 60) staleness = "";
    }
    card.className = "agent-card " + base + (staleness ? " " + staleness : "");

    var head = document.createElement("div");
    head.className = "head";

    var nameEl = document.createElement("span");
    nameEl.className = "name " + base;
    nameEl.textContent = name;
    head.appendChild(nameEl);

    if (role) {
      var roleEl = document.createElement("span");
      roleEl.className = "role";
      roleEl.textContent = role;
      head.appendChild(roleEl);
    }

    var pct = document.createElement("span");
    pct.className = "pct";
    pct.textContent = agent ? ((agent.percent || 0) + "%") : "—";
    head.appendChild(pct);
    card.appendChild(head);

    var progress = document.createElement("div");
    progress.className = "progress";
    var bar = document.createElement("div");
    bar.className = "bar";
    var pctValue = agent ? Math.max(0, Math.min(100, agent.percent || 0)) : 0;
    bar.style.width = pctValue + "%";
    progress.appendChild(bar);
    card.appendChild(progress);

    var activity = document.createElement("div");
    activity.className = "activity";
    activity.textContent = agent ? (agent.activity || "(idle)") : "(no status yet)";
    card.appendChild(activity);

    var meta = document.createElement("div");
    meta.className = "meta";
    var topicEl = document.createElement("span");
    topicEl.textContent = agent && agent.topic ? "[" + agent.topic + "]" : "";
    var ageEl = document.createElement("span");
    ageEl.textContent = ageLabel(age);
    meta.appendChild(topicEl);
    meta.appendChild(ageEl);
    card.appendChild(meta);

    return card;
  }

  function renderAgentCards(agents) {
    if (agents) state.agents = agents;
    if (!agentCardsEl) return;
    agentCardsEl.innerHTML = "";
    var elapsedSinceFetch = Math.floor((Date.now() - state.agentsFetchedAt) / 1000);
    var byName = {};
    (state.agents || []).forEach(function (a) { byName[a.agent] = a; });
    PANE_ORDER.forEach(function (paneName) {
      var direct = byName[paneName];
      // Backward compat: when only the legacy "claude"/"codex" row exists,
      // mirror it onto both chat and work cards so the sidebar shows
      // something useful until commit#5 splits the heartbeats.
      if (!direct) {
        var base = paneBase(paneName);
        if (byName[base]) direct = byName[base];
      }
      agentCardsEl.appendChild(buildAgentCard(paneName, direct, elapsedSinceFetch));
    });
  }

  function renderClaims(claims) {
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

  function poll() {
    var url = "/messages.json?after=" + encodeURIComponent(state.lastId);
    fetch(url, {credentials: "same-origin", headers: {"Accept": "application/json"}})
      .then(function (r) {
        if (!r.ok) throw new Error("status " + r.status);
        return r.json();
      })
      .then(function (data) {
        var prevCursors = JSON.stringify(state.cursors);
        state.cursors = data.cursors || {};
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
        renderAgents(data.agents || []);
        renderAgentCards();
      })
      .catch(function (err) {
        // network blip — try again next tick
        if (window.console && console.warn) console.warn("poll failed", err);
      });
  }

  function autoresize(el) {
    el.style.height = "auto";
    el.style.height = Math.min(240, el.scrollHeight) + "px";
  }

  function sendMessage() {
    var body = bodyEl.value.trim();
    if (!body) return;
    sendBtn.disabled = true;
    var payload = {
      to: document.getElementById("composer-to").value,
      topic: document.getElementById("composer-topic").value,
      status: document.getElementById("composer-status").value,
      body: bodyEl.value,
      refs: (document.getElementById("composer-refs").value || "")
        .split(",").map(function (s) { return s.trim(); }).filter(Boolean),
    };
    fetch("/send", {
      method: "POST",
      credentials: "same-origin",
      headers: {"Content-Type": "application/json", "Accept": "application/json"},
      body: JSON.stringify(payload),
    })
      .then(function (r) { return r.json().then(function (j) { return {ok: r.ok, body: j}; }); })
      .then(function (res) {
        if (!res.ok) {
          alert((res.body && res.body.error) || "send failed");
          return;
        }
        bodyEl.value = "";
        autoresize(bodyEl);
        poll();
      })
      .catch(function (err) { alert("send failed: " + err); })
      .finally(function () { sendBtn.disabled = false; bodyEl.focus(); });
  }

  composer.addEventListener("submit", function (e) {
    e.preventDefault();
    sendMessage();
  });

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

  toggleSidebar.addEventListener("click", function () {
    sidebar.classList.toggle("hidden");
  });

  // Initial paint from server-supplied state.
  state.cursors = initial.cursors || {};
  renderCursors();
  renderMessages(initial.messages || []);
  renderClaims(initial.claims || []);
  renderArtifacts(initial.artifacts || []);
  renderAgents(initial.agents || []);
  renderAgentCards();
  scrollToBottom();
  autoresize(bodyEl);

  setInterval(poll, 3000);
  // Re-render agent strip + cards every 10s so age labels tick without a server roundtrip.
  setInterval(function () { renderAgents(); renderAgentCards(); }, 10000);
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
