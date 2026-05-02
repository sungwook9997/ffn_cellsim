# `acs-collab-mcp` — Claude ↔ Codex collaboration channel

Streamable-HTTP MCP server backing the Claude (Mac) ↔ Codex (Windows)
async channel for the ActiveCellSim project. Hosts on the Windows
workstation (24/7, RTX A5000 box); both LLMs connect over Tailscale.

The server replaces the Syncthing-only `docs/claude_codex_log.md`
turn-based protocol with a real-time message queue + topic locking +
artifact-pointer registry. The ledger file is kept as the durable
human-readable mirror — every `send` automatically appends a block to
it when `COLLAB_MCP_LEDGER` is configured.

The first version exposes **no command-execution surface**. Tools are
limited to message append/read, topic claim/release, status query, and
artifact-pointer announcement. GPU job launches must go through a
separately-approved wrapper added in a later iteration.

`room.py` adds a browser-facing PI / Claude / Codex room over the same
SQLite DB. It is also message-only: the browser can view the timeline,
active claims, artifact pointers, and post PI messages, but it cannot
run commands or read artifact contents. The UI is a KakaoTalk-style chat
(right-aligned PI bubbles, left-aligned LLM bubbles, sticky composer,
Enter to send / Shift+Enter for newline, `@claude` / `@codex` / `@both`
mentions, day dividers, seen indicators) and a top-bar agent status
strip showing each agent's current % progress, activity, and freshness.

`launch_workroom_mac.command` opens the same room as a local app-style
Chrome window on macOS. It starts the room server if needed, then opens
`http://127.0.0.1:7879/` in Chrome app mode. `desktop_room.py` is an
experimental Tkinter client kept as a fallback, but macOS system Tk can
render unreliably on some machines.

## Tools

| Tool | Purpose |
| ---- | ------- |
| `send(to, topic, body, status, refs)` | Post a message; optional ledger append |
| `read(unread_only=True, limit=20)` | Fetch messages (advances cursor when unread) |
| `claim(topic, intent)` | Atomically take ownership of a work topic |
| `release(topic, summary)` | Hand back a topic with a one-line summary |
| `status()` | Active claims + per-author unread counts + latest message |
| `announce_artifact(kind, path, note)` | Register a figure / HDF5 / dashboard path |
| `list_artifacts(limit=20)` | Recent artifact pointers |
| `bootstrap(limit=50, topic, since_id)` | Read recent messages without advancing cursor (warm up a fresh session) |

Allowed `status` values: `FYI`, `proposal`, `review`, `decision-needed`,
`open-question`, `blocker`. Allowed authors: `claude`, `codex`
(override via `COLLAB_MCP_AUTHORS`).

A built-in rate guard rejects a 4th consecutive message from the same
author with no reply in between.

## Configuration

All settings come from environment variables; nothing is committed.

| Variable | Required? | Default | Notes |
| -------- | --------- | ------- | ----- |
| `COLLAB_MCP_TOKEN` | yes | — | Shared bearer token. Set identically on server + both clients. |
| `COLLAB_MCP_BIND`  | no  | `127.0.0.1:7878` | Server listen address. Use `0.0.0.0:<port>` only when intentionally exposing over Tailscale. |
| `COLLAB_MCP_DB`    | no  | `~/.acs-collab/inbox.db` | SQLite file (server-local). |
| `COLLAB_MCP_LEDGER`| no  | unset | If set, every `send` appends to this file. Recommended: `<repo>/docs/claude_codex_log.md` on the host machine. |
| `COLLAB_MCP_AUTHORS`| no | `claude,codex` | Comma-separated allowlist. |
| `COLLAB_ROOM_BIND` | no | `127.0.0.1:7879` | Browser room listen address. Use `0.0.0.0:<port>` only behind Tailscale/firewall. |
| `COLLAB_ROOM_TOKEN` | no | `COLLAB_MCP_TOKEN` | Browser login token. Set separately if PI should not use the MCP bearer token. |
| `COLLAB_ROOM_AUTHOR` | no | `pi` | Author label for browser-submitted messages. |

Per-request authentication / identity (clients send these headers):

```
Authorization: Bearer <COLLAB_MCP_TOKEN>
X-Collab-Author: claude   # or 'codex'
```

## Running the server (Windows host)

1. Install. From a Win terminal in this repo:
   ```powershell
   python -m venv .venv-collab
   .venv-collab\Scripts\Activate.ps1
   pip install -r tools\collab_mcp\requirements.txt
   ```
2. Set env vars (PowerShell session, or persist via `setx`):
   ```powershell
   $env:COLLAB_MCP_TOKEN = "<long-random-string>"
   # Use 0.0.0.0 only after restricting access to the Tailscale interface.
   $env:COLLAB_MCP_BIND  = "0.0.0.0:7878"
   $env:COLLAB_MCP_LEDGER = "C:\Users\sw1\ActiveCellSim\docs\claude_codex_log.md"
   ```
3. Run in foreground for first-test:
   ```powershell
   python tools\collab_mcp\server.py
   ```
4. Background service (optional, after foreground works): register via
   NSSM or Task Scheduler. **Do this only after PI confirmation** — see
   `windows_service.md` (Codex to write).

5. Allow inbound TCP 7878 on the Tailscale interface only (not the
   public network). Codex confirms the firewall rule with PI before
   committing it.

## Running the browser room

Run this on the same host as the MCP server so both processes share
`COLLAB_MCP_DB` and `COLLAB_MCP_LEDGER`.

```powershell
$env:COLLAB_MCP_DB = "$HOME\.acs-collab\inbox.db"
$env:COLLAB_MCP_LEDGER = "C:\Users\sw1\ActiveCellSim\docs\claude_codex_log.md"
$env:COLLAB_ROOM_TOKEN = "<browser-login-token>"
$env:COLLAB_ROOM_BIND = "127.0.0.1:7879"
python tools\collab_mcp\room.py
```

For Tailscale browser access from PI devices:

```powershell
# Use only after restricting inbound access to Tailscale/private profile.
$env:COLLAB_ROOM_BIND = "0.0.0.0:7879"
python tools\collab_mcp\room.py
```

Open `http://<windows-tailscale-name>:7879/` in a browser and enter
`COLLAB_ROOM_TOKEN`. PI messages are written as author `pi`, addressed
to `claude,codex`, and are visible to both LLM clients through normal
MCP `read()` calls.

### Live agent progress strip

Claude and Codex publish their own progress to the room via a small JSON
endpoint. PI sees the result as pills in the top bar — name, %, current
activity, and freshness (fresh / aging / stale).

```bash
curl -s -b cookies.txt -H 'Content-Type: application/json' -X POST \
  -d '{"agent":"claude","percent":60,"activity":"writing room.py JS","topic":"ui-redesign-claude"}' \
  http://127.0.0.1:7879/agent_status
```

Validation: `agent` must be `claude`/`codex`/`pi`; `percent` is 0..100;
`activity` capped at 240 chars; `topic` capped at 120 chars. Auth uses
`COLLAB_ROOM_TOKEN` via the same browser cookie or by adding the cookie
header to the curl call. The pill turns gray after 5 minutes without an
update.

## Running the app-style room on macOS

From Finder, double-click:

```text
tools/collab_mcp/launch_workroom_mac.command
```

Or from Terminal:

```bash
tools/collab_mcp/launch_workroom_mac.command
```

Defaults:

```text
URL:   http://127.0.0.1:7879/
Token: acs-room
DB:    ~/.acs-collab/inbox.db
Log:   docs/claude_codex_log.md
```

This is the recommended local surface on the current Mac. It is still
message-only and writes to the same DB/ledger as MCP.

## Running the tmux relay

`tmux_relay.py` is the first group-chat bridge. It polls the shared
SQLite message table and injects only a compact wrapper prompt into
Claude/Codex tmux panes. It never injects raw message bodies, so quotes,
newlines, code blocks, and Korean text stay in SQLite/MCP instead of
being escaped through `tmux send-keys`.

Install tmux once on macOS:

```bash
brew install tmux
```

Recommended layout separates fast PI chat from long-running work:

```text
claude-chat  receives PI workroom messages
codex-chat   receives PI workroom messages
claude-work  long implementation tasks
codex-work   long review/implementation tasks
```

Create the sessions and launch the chat CLIs:

```bash
python3 tools/collab_mcp/setup_tmux_workroom.py --start-chat
```

Start work CLIs too only when you want separate implementation panes:

```bash
python3 tools/collab_mcp/setup_tmux_workroom.py --start-work
```

Inspect targets:

```bash
tmux list-panes -a -F '#{session_name}:#{window_index}.#{pane_index} #{pane_current_command}'
```

Run a dry-run first:

```bash
COLLAB_TMUX_CLAUDE_TARGET=claude-chat:0.0 \
COLLAB_TMUX_CODEX_TARGET=codex-chat:0.0 \
python3 tools/collab_mcp/tmux_relay.py --once --dry-run
```

Before the first live run on an existing conversation, initialize the
relay cursor so old PI messages are not replayed into both panes:

```bash
python3 tools/collab_mcp/tmux_relay.py --init-cursor-now
```

Run the live relay:

```bash
COLLAB_TMUX_CLAUDE_TARGET=claude-chat:0.0 \
COLLAB_TMUX_CODEX_TARGET=codex-chat:0.0 \
python3 tools/collab_mcp/tmux_relay.py
```

Default policy:

- PI-authored messages route to `@claude`, `@codex`, or both.
- Claude/Codex-authored messages are not auto-routed by default, which
  prevents accidental LLM ping-pong.
- Add `--include-llm-messages` only for a deliberate Claude↔Codex
  routing test.
- The relay stores its own cursor in `relay_cursors`, independent of
  Claude/Codex read cursors, so messages are not injected twice.

## Building the macOS .app bundle

For a Dock-pinnable native app surface that wraps the same browser room,
build the `.app` bundle once:

```bash
python3 tools/collab_mcp/build_workroom_app.py
```

Defaults:

```text
output:   ~/Applications/ACSWorkroom.app
name:     ACSWorkroom
bundleId: com.activecellsim.workroom
```

The build is stdlib-only (no Xcode, no extra packages). The launcher
inside the bundle is the same zsh script as `launch_workroom_mac.command`,
with the repo root baked in at build time. Re-run the build script if
the repo is moved.

Override knobs:

```bash
python3 tools/collab_mcp/build_workroom_app.py \
    --output /Applications \
    --name ACSWorkroom \
    --bundle-id com.activecellsim.workroom
```

Once built, double-click the .app from Finder, drag it to the Dock for
a one-click launch, or run `open ~/Applications/ACSWorkroom.app` from
the terminal.

## Running the experimental Tk desktop room

On macOS or Windows, from the repo root:

```powershell
$env:COLLAB_MCP_LEDGER = "C:\Users\sw1\ActiveCellSim\docs\claude_codex_log.md"
python tools\collab_mcp\desktop_room.py
```

On macOS zsh:

```bash
COLLAB_MCP_LEDGER=/Users/sw1/ActiveCellSim/docs/claude_codex_log.md \
python3 tools/collab_mcp/desktop_room.py
```

The desktop app reads `COLLAB_MCP_DB` if set; otherwise it uses
`~/.acs-collab/inbox.db`, matching the MCP server default. It does not
start Claude or Codex by itself. It is the local workroom surface; the
LLM clients still need to read the MCP queue or run a later relay.
On this Mac, the system Tk window can appear blank; use the app-style
Chrome launcher above when that happens.

## Connecting Claude Code (Mac client)

Edit `~/.claude/settings.json` (user-level, applies across sessions):

```json
{
  "mcpServers": {
    "acs-collab": {
      "type": "http",
      "url": "http://<windows-tailscale-name>:7878/mcp/",
      "headers": {
        "Authorization": "Bearer ${COLLAB_MCP_TOKEN}",
        "X-Collab-Author": "claude"
      }
    }
  }
}
```

`COLLAB_MCP_TOKEN` is read from the shell that launches Claude Code, so
export it in `~/.zshenv` (or equivalent). Restart Claude Code after
editing the config so the new server is loaded.

## Connecting Codex CLI (Windows client)

Edit `~/.codex/config.toml` (or wherever the local Codex CLI keeps
MCP servers):

```toml
[mcp_servers.acs-collab]
type = "http"
url  = "http://localhost:7878/mcp/"
headers = { Authorization = "Bearer ${COLLAB_MCP_TOKEN}", "X-Collab-Author" = "codex" }
```

(Exact key names depend on the Codex CLI version Codex is running on
the Win box. Codex confirms the right schema and updates this README
when the connection is verified.)

## Smoke test

After both clients are configured, run an echo round-trip:

1. Mac (Claude): call `send(to="codex", topic="mcp-handshake", body="echo from claude", status="FYI")`.
2. Win (Codex): call `read()`. Should return one message with
   `from=claude, topic=mcp-handshake`.
3. Win (Codex): call `send(to="claude", topic="mcp-handshake", body="ack", status="FYI")`.
4. Mac (Claude): call `read()`. Should return Codex's ack.
5. Either side: call `status()` — should show 0 unread, no active claims.

If the ledger path is configured, `docs/claude_codex_log.md` will have
both messages auto-appended in the standard entry block format.
