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
            """
        )


def _now_human() -> str:
    return dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M (%Z)")


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


def _load_state() -> dict[str, Any]:
    with closing(_db()) as conn:
        messages = conn.execute(
            "SELECT id, ts, author, addressee, topic, body, status, refs "
            "FROM messages ORDER BY id DESC LIMIT 100"
        ).fetchall()
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
    return {
        "messages": list(reversed(messages)),
        "claims": claims,
        "artifacts": artifacts,
        "cursors": cursors,
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


def _render_room(state: dict[str, Any], flash: str = "") -> str:
    cursor_map = {author: int(last_seen) for author, last_seen in state["cursors"]}
    messages_html = "\n".join(_render_message(row, cursor_map) for row in state["messages"])
    claims_html = "\n".join(_render_claim(row) for row in state["claims"]) or "<li>None</li>"
    artifacts_html = "\n".join(_render_artifact(row) for row in state["artifacts"]) or "<li>None</li>"
    cursors_html = "\n".join(
        f"<li><span>{_esc(author)}</span><code>{last_seen}</code></li>"
        for author, last_seen in state["cursors"]
    ) or "<li>No cursors yet</li>"
    flash_html = f"<p class='flash'>{_esc(flash)}</p>" if flash else ""
    status_options = "\n".join(
        f"<option value='{_esc(status)}'{' selected' if status == 'open-question' else ''}>{_esc(status)}</option>"
        for status in sorted(ALLOWED_STATUSES)
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ACS Collab Room</title>
  <style>{CSS}</style>
</head>
<body>
  <header>
    <div>
      <h1>ACS Collab Room</h1>
      <p>Shared PI / Claude / Codex channel over the MCP SQLite queue.</p>
    </div>
    <nav>
      <a href="/">Refresh</a>
      <a href="/logout">Logout</a>
    </nav>
  </header>
  <main class="layout">
    <section class="timeline" aria-label="Message timeline">
      {flash_html}
      <form class="composer" method="post" action="/send">
        <div class="row">
          <label>Topic
            <input name="topic" value="v2-layer-1" required>
          </label>
          <label>To
            <input name="to" value="claude,codex" required>
          </label>
          <label>Status
            <select name="status">{status_options}</select>
          </label>
        </div>
        <label>Message
          <textarea name="body" rows="5" required></textarea>
        </label>
        <label>Refs
          <input name="refs" placeholder="optional comma-separated paths or topics">
        </label>
        <button type="submit">Send as PI</button>
      </form>
      <div class="messages">{messages_html}</div>
    </section>
    <aside>
      <section>
        <h2>Active Claims</h2>
        <ul>{claims_html}</ul>
      </section>
      <section>
        <h2>Artifacts</h2>
        <ul>{artifacts_html}</ul>
      </section>
      <section>
        <h2>Cursors</h2>
        <ul>{cursors_html}</ul>
      </section>
    </aside>
  </main>
  <script>
    window.setTimeout(function () {{
      var active = document.activeElement;
      var composing = active && active.closest && active.closest(".composer");
      if (!composing) window.location.reload();
    }}, 5000);
  </script>
</body>
</html>"""


def _render_message(row: sqlite3.Row | tuple[Any, ...], cursor_map: dict[str, int]) -> str:
    msg_id, ts, author, addressee, topic, body, status, refs_json = row
    try:
        refs = json.loads(refs_json)
    except json.JSONDecodeError:
        refs = []
    refs_html = ""
    if refs:
        refs_html = "<p class='refs'>Refs: " + ", ".join(_esc(ref) for ref in refs) + "</p>"
    seen_by = [
        name
        for name in ("claude", "codex")
        if name != author and cursor_map.get(name, 0) >= int(msg_id)
    ]
    seen_html = (
        f"<p class='seen'>Seen by {', '.join(_esc(name) for name in seen_by)}</p>"
        if seen_by
        else "<p class='seen pending'>Not seen by Claude/Codex yet</p>"
    )
    body_html = "<br>".join(_esc(body).splitlines())
    return f"""<article class="message author-{_esc(author)}">
  <div class="meta">
    <strong>{_esc(author)}</strong>
    <span>to {_esc(addressee)}</span>
    <span>{_esc(status)}</span>
    <span>{_esc(topic)}</span>
    <time>{_esc(ts)}</time>
    <code>#{msg_id}</code>
  </div>
  <p>{body_html}</p>
  {refs_html}
  {seen_html}
</article>"""


def _render_claim(row: sqlite3.Row | tuple[Any, ...]) -> str:
    topic, author, claimed_at, summary = row
    return (
        f"<li><strong>{_esc(topic)}</strong><br>"
        f"<span>{_esc(author)} since {_esc(claimed_at)}</span><br>"
        f"<small>{_esc(summary or '')}</small></li>"
    )


def _render_artifact(row: sqlite3.Row | tuple[Any, ...]) -> str:
    artifact_id, ts, author, kind, path, note = row
    return (
        f"<li><strong>{_esc(kind)}</strong> <code>#{artifact_id}</code><br>"
        f"<code>{_esc(path)}</code><br>"
        f"<span>{_esc(author)} · {_esc(ts)}</span><br>"
        f"<small>{_esc(note or '')}</small></li>"
    )


class RoomHandler(BaseHTTPRequestHandler):
    server_version = "ACSCollabRoom/0.1"

    def do_GET(self) -> None:
        path = urllib.parse.urlparse(self.path).path
        if path == "/logout":
            self._send_text(
                _render_login("Logged out."),
                headers={"Set-Cookie": "collab_room_token=; Path=/; Max-Age=0; SameSite=Lax"},
            )
            return
        if path != "/":
            self.send_error(HTTPStatus.NOT_FOUND, "not found")
            return
        if not self._is_authenticated():
            self._send_text(_render_login())
            return
        self._send_text(_render_room(_load_state()))

    def do_POST(self) -> None:
        path = urllib.parse.urlparse(self.path).path
        try:
            if path == "/login":
                self._handle_login()
            elif path == "/send":
                self._handle_send()
            else:
                self.send_error(HTTPStatus.NOT_FOUND, "not found")
        except ValueError as exc:
            self._send_text(_render_room(_load_state(), flash=str(exc)), status=HTTPStatus.BAD_REQUEST)

    def log_message(self, fmt: str, *args: Any) -> None:
        sys.stderr.write(f"[collab-room] {self.address_string()} - {fmt % args}\n")

    def _handle_login(self) -> None:
        fields = self._read_form()
        token = fields.get("token", [""])[0]
        if not self._token_ok(token):
            self._send_text(_render_login("Invalid token."), status=HTTPStatus.UNAUTHORIZED)
            return
        self._redirect("/", cookie=f"collab_room_token={urllib.parse.quote(token)}; Path=/; SameSite=Lax")

    def _handle_send(self) -> None:
        if not self._is_authenticated():
            self._send_text(_render_login("Please enter the room token."), status=HTTPStatus.UNAUTHORIZED)
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

    def _read_form(self) -> dict[str, list[str]]:
        length = int(self.headers.get("Content-Length", "0"))
        if length > 65536:
            raise ValueError("request body too large")
        body = self.rfile.read(length).decode("utf-8")
        return urllib.parse.parse_qs(body, keep_blank_values=True)

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

    def _redirect(self, location: str, cookie: str | None = None) -> None:
        self.send_response(HTTPStatus.SEE_OTHER)
        self.send_header("Location", location)
        if cookie:
            self.send_header("Set-Cookie", cookie)
        self.end_headers()


CSS = """
:root {
  color-scheme: light;
  --bg: #f5f6f8;
  --panel: #ffffff;
  --ink: #1d232f;
  --muted: #697386;
  --line: #d9dee8;
  --accent: #246bfe;
  --accent-ink: #ffffff;
  --pi: #fff6d9;
  --claude: #eef6ff;
  --codex: #edf9f1;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font: 14px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}
header {
  align-items: center;
  background: var(--panel);
  border-bottom: 1px solid var(--line);
  display: flex;
  gap: 24px;
  justify-content: space-between;
  padding: 16px 24px;
}
h1, h2, p { margin: 0; }
h1 { font-size: 22px; }
h2 { font-size: 15px; margin-bottom: 8px; }
header p, small, .meta, .refs, aside span { color: var(--muted); }
nav { display: flex; gap: 12px; }
a { color: var(--accent); text-decoration: none; }
.layout {
  display: grid;
  gap: 16px;
  grid-template-columns: minmax(0, 1fr) 320px;
  padding: 16px;
}
.timeline, aside section, .composer, .message, .login-panel {
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 8px;
}
.timeline { padding: 16px; }
.composer { margin-bottom: 16px; padding: 14px; }
.row {
  display: grid;
  gap: 10px;
  grid-template-columns: minmax(160px, 1fr) minmax(140px, 220px) minmax(140px, 200px);
}
label { display: grid; gap: 5px; font-weight: 600; }
input, select, textarea {
  border: 1px solid var(--line);
  border-radius: 6px;
  color: var(--ink);
  font: inherit;
  padding: 8px 10px;
  width: 100%;
}
textarea { min-height: 116px; resize: vertical; }
button {
  background: var(--accent);
  border: 0;
  border-radius: 6px;
  color: var(--accent-ink);
  cursor: pointer;
  font: inherit;
  font-weight: 700;
  margin-top: 10px;
  padding: 9px 14px;
}
.messages { display: grid; gap: 10px; }
.message { padding: 12px; }
.author-pi { background: var(--pi); }
.author-claude { background: var(--claude); }
.author-codex { background: var(--codex); }
.meta {
  align-items: center;
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  margin-bottom: 8px;
}
.seen {
  color: #2f7a45;
  font-size: 12px;
  margin-top: 8px;
}
.seen.pending { color: var(--muted); }
code {
  background: rgba(29, 35, 47, 0.06);
  border-radius: 5px;
  padding: 1px 5px;
}
aside { display: grid; gap: 16px; align-content: start; }
aside section { padding: 14px; }
ul { list-style: none; margin: 0; padding: 0; display: grid; gap: 10px; }
.flash, .error {
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
  display: grid;
  gap: 12px;
  max-width: 420px;
  padding: 24px;
  width: min(92vw, 420px);
}
@media (max-width: 900px) {
  header, .layout { padding: 12px; }
  .layout { grid-template-columns: 1fr; }
  .row { grid-template-columns: 1fr; }
}
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
