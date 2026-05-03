# ACS Collab — 2-Workroom Layout SOP

**Status:** active operating model after the Task 1–5 split
(commits `6a1a134`, `ee9d3a4`, `35c5112`, `244fec0`, `083ae55`,
`1ea6695`, `04df930`, 2026-05-03).
**Replaces:** the single 4-pane workroom that ran Claude+Codex chat/work
inside one tmux scope.

## Why two rooms

The single-workroom layout mixed pure design discussion (PI ↔ Claude ↔
Codex on Layer 1/2 modeling) with implementation work (Claude/Codex
editing files, running tests, talking to the Win SSH pane). That wedged
two distinct cadences into one MCP timeline:

- Discussion is high-frequency, short-ack, model-quality-sensitive.
- Implementation is low-frequency, long-output, claim/lock-sensitive.

When PI sent a discussion message, the work panes kept reciprocal-ack-ing
or polluting `/agent_status` with "idle" rows; when PI dispatched code,
the chat panes kept hot-trying to keep the 5-second cadence and ate
context. Splitting the rooms turns the same MCP DB into two parallel
timelines so each cadence stays clean.

## Rooms

| Room | Purpose | tmux session prefix | When PI talks here |
|---|---|---|---|
| `design-discussion` | Layer 1/2 modeling, sanity-gate review, paper-anchored debate | `design-discussion-` | Any pure discussion / decision / spec work |
| `implementation-work` | Code edits, test runs, lit/data pipelines, Win-A5000 jobs | `implementation-work-` | Any task that produces a diff, a commit, or a measurement |

A workroom name is `^[a-z0-9][a-z0-9-]{0,31}$` (lowercase alphanum +
dash, ≤32 chars, must start with alphanum). The launcher, MCP server,
relay, and heartbeat all enforce the same regex.

Adding a third room only when **both** of the following hold:
1. The work it carries is genuinely independent of either of the two
   above (separate decision context, separate cadence) — usually a
   multi-day independent stream.
2. PI explicitly says so. Don't fork rooms unilaterally.

## What is per-room vs shared

**Per-room** — one instance per workroom, prefixed with `<wr>-`:
- tmux sessions: `<wr>-claude-chat`, `<wr>-codex-chat`,
  `<wr>-claude-work`, `<wr>-codex-work`, `<wr>-win-ssh`,
  `<wr>-relay`, `<wr>-heartbeat`.
- Log dir: `/tmp/acs-collab/<wr>/`. Sentinel files
  (`.mcp_python_unhealthy`, `.mcp_server_unhealthy`,
  `.mcp_token_unset`, `.selfcheck_failed`) live inside.
- Per-pane Claude MCP config: `/tmp/acs-collab/<wr>/claude_mcp.json`,
  chmod 0600, with `X-Collab-Room: <wr>` baked in. Claude CLI is
  launched with `--mcp-config <path> --strict-mcp-config` so the
  global `~/.claude.json` does not override the room header.
- Codex CLI is launched with `-c
  'mcp_servers.acs-collab.http_headers={ "X-Collab-Author" = "codex",
  "X-Collab-Room" = "<wr>" }'` for the same purpose.
- MCP `relay_cursors` row (one per `<wr>-relay`).
- MCP `cursors` rows (per-`(author, room)` key after Task 2 migration).
- MCP `claims` rows (per-`(room, topic)` key after Task 2 migration).

**Shared singletons** — exactly one process across all rooms:
- `mcp` tmux session (server.py, port 7878).
- `room` tmux session (room.py UI, port 7879).
- SQLite DB at `~/.acs-collab/inbox.db` (single file, room column
  isolates the data).
- Shared lock dir `/tmp/acs-collab/locks/` for the `ensure_tmux_session_locked`
  helper that serializes singleton/per-room session creation across
  concurrent launchers.
- Chrome `--app` lock at `/tmp/acs-collab/chrome_app.lock`.
- `docs/claude_codex_log.md` ledger file (cross-room append).

**Hard rule**: do not split the MCP DB or the ledger per room. The
isolation contract lives at the row level (the `room` column on
`messages`, `cursors`, `claims`); creating separate stores would lose
cross-room context that PI sometimes needs (e.g., when an
implementation-work commit references a design-discussion decision).

## Launcher

```sh
tools/collab_mcp/launch_workroom_mac.command [workroom-name]
```

Default: `design-discussion`. The launcher is idempotent at every
tier — re-running with the same room reuses every existing session;
running with a new room name boots a parallel set without touching the
first.

What the launcher creates, in order:

1. Python health gate (refuses system `python3` fallback; needs
   `.venv-collab/bin/python` with `mcp[cli]>=1.2`).
2. Shared `mcp` singleton if missing (lock-protected).
3. Per-room tmux pane sessions via `setup_tmux_workroom.py
   --workroom <wr> --start-chat --start-work --start-win-ssh
   --start-heartbeat`.
4. Per-room `<wr>-relay` session running `tmux_relay.py --relay-name
   <wr> --workroom <wr> --claude-target <wr>-claude-chat:0.0
   --codex-target <wr>-codex-chat:0.0 --include-llm-messages`.
5. Shared `room` UI singleton if `http://127.0.0.1:7879/` is not yet
   answering (lock-protected).
6. Chrome `--app` to the room URL (idempotent across concurrent
   launches via the Chrome lock).
7. Self-check that re-verifies every per-room tmux session and the
   two singletons; failure surfaces via stderr, sentinel file, and
   desktop notification.

Two launchers running side-by-side cannot collide on session create —
the `ensure_tmux_session_locked` helper takes a `shlock`/`mkdir` lock
per session name and re-probes after the lock so the second launcher
either reuses the first's session or fails loud with a non-duplicate
error.

## MCP routing

Every MCP request now carries an implicit room:

- **Server source of truth**: the `X-Collab-Room` HTTP header. The
  per-pane Claude / Codex configs above bake the workroom name in.
- **Fallback (no header)**: server env `COLLAB_DEFAULT_ROOM`, then the
  hard default `design-discussion`. Pre-Task-2 messages migrate into
  `design-discussion`.
- **Per-call override**: `send/read/bootstrap/claim/release/status` all
  accept `room=<wr>`. `room='*'` is the cross-room sentinel; on
  `read` it requires `unread_only=False` because per-`(author, room)`
  cursors cannot honor a global unread query without skipping.

`route_targets`/`dispatch_message` in `tmux_relay.py` filter on the
relay's workroom — a `<wr>-relay` only injects wrappers for messages
tagged with `<wr>` (and the legacy `room=''` corpus). Cross-room
messages still advance the cursor so other rooms' traffic doesn't
re-fetch forever, and `notify_desktop` fires only for messages that
actually belong to the relay.

Wrapper prompt fields (what chat panes see):

```
: mcp_msg id=<n> from=<author> to=<addr> topic=<t> room=<wr>
  pane=chat work_pane=<wr>-<agent>-work priority=immediate ack_first
  action=read_acs_collab_mcp_then_send_if_<agent>_should_reply
```

The `work_pane=` name is the physical prefixed session, so chat panes
brief the right sibling via `/tmp/acs-collab/work_briefing.md`.

## Claim / cursor / status semantics

| API | Default room | `room='*'` |
|---|---|---|
| `send`     | caller's room | rejected (room must be a single name) |
| `read`     | caller's room | OK only with `unread_only=False` |
| `bootstrap`| caller's room | per-room cursor advance (each touched room → its own MAX) |
| `claim`    | caller's room | rejected |
| `release`  | caller's room | rejected |
| `status`   | caller's room (`active_claims`); cross-room `claims_by_room` and `unread_by_room` always present |

This means two rooms can hold a claim on the same logical topic
(`v2-layer-1`) without conflict — Codex reviewing Layer 1 in
`design-discussion` and Claude implementing in `implementation-work`
won't collide.

## Sidebar UI

Topbar `Room` selector lets PI pick the active workroom for compose +
timeline view. JS poll, send, and agent-card lookup all carry the
selected room. Agent cards prefer `<currentRoom>-<logical>` ids and
fall back to bare logical / suffix-match for legacy unprefixed
sessions, so the legacy 4-pane mode and the new 2-room mode both
render usefully.

Empty-DB defaults: the topbar `Room` select shows
`design-discussion` and `implementation-work` even before any traffic
exists, so PI can post into either room from a cold-start.

## Cross-room context — preventing silent split

The single MCP DB + the cross-room sidebar + the shared
`docs/claude_codex_log.md` ledger are the three places that prevent
"context split is silent" failure mode. Operating SOP:

- **Every milestone or decision** that affects another room: send a
  one-line MCP message in your own room mentioning the cross-room
  reference (`refs=[mcp_msg:<id>]`).
- **At the end of a substantive turn** in either room: append a 3–5
  line summary to `docs/claude_codex_log.md`. The MCP server appends
  every send automatically when `COLLAB_MCP_LEDGER` is set; when room
  splits matter, write the cross-room implication explicitly.
- **Status cards are room-scoped**, but `claims_by_room` /
  `unread_by_room` in `status()` are not. Use them to spot a
  starvation pattern (e.g., `implementation-work` has 12 unread for
  Claude while `design-discussion` is at 0).

## Adding a new room

1. Pick a name that matches `^[a-z0-9][a-z0-9-]{0,31}$`.
2. `tools/collab_mcp/launch_workroom_mac.command <name>` from a fresh
   shell with `COLLAB_MCP_TOKEN` exported (the launchd plist already
   exports it for the default boot).
3. Verify the post-launch self-check line: `[SELF-CHECK OK] workroom=<name>
   all <n> tiers up`. Failure leaves a sentinel under
   `/tmp/acs-collab/<name>/` and a desktop notification.
4. Add a one-line entry to this doc's room table and to
   `~/.claude/projects/-Users-sw1-ActiveCellSim/memory/MEMORY.md`
   under the existing `워크룸 운영` entry.

## Removing a room

1. `tmux kill-session -t <wr>-<logical>` for each per-room session
   (`-claude-chat`, `-claude-work`, `-codex-chat`, `-codex-work`,
   `-win-ssh`, `-relay`, `-heartbeat`).
2. `rm -rf /tmp/acs-collab/<wr>/`.
3. Optionally wipe DB rows: `DELETE FROM messages WHERE room='<wr>';`
   plus the matching `cursors` and `claims` rows. Don't do this if the
   ledger references those messages.
4. Remove the row from this doc + memory entry.

## Migration history

| Date | Commit | Note |
|---|---|---|
| 2026-05-03 | `6a1a134` | Task 1+2 + Task 3 server-side prefixed-agent acceptance + composer/sidebar grouping |
| 2026-05-03 | `ee9d3a4` | Task 3 default room presets |
| 2026-05-03 | `35c5112` | Task 4 collision-safe parallel boot |
| 2026-05-03 | `244fec0` | Task 6 2-room SOP + Claude memory entry |
| 2026-05-03 | `083ae55` | Room selector persistence, PI reply path, and new-relay cursor hotfixes |
| 2026-05-03 | runtime | Task 5 live migration: legacy unprefixed sessions retired after both prefixed rooms passed self-check |
| 2026-05-03 | `04df930` | Prefixed heartbeat rows seeded so UI agent cards show per-room panes instead of stale legacy rows |

Pre-Task-2 messages were bucketed into `design-discussion` by the
schema migration. Pre-Task-1 unprefixed tmux sessions were retired
during Task 5 after `design-discussion` and `implementation-work`
both cold-started successfully. The active runtime should now contain
only prefixed per-room sessions plus the shared `mcp` and `room`
singletons.
