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
