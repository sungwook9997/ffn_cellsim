#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="${0:A:h}"
REPO_ROOT="${SCRIPT_DIR:h:h}"
cd "$REPO_ROOT"

# ---------------------------------------------------------------------------
# Workroom selection
# ---------------------------------------------------------------------------
# Usage: launch_workroom_mac.command [workroom-name]
#   default = $ACS_WORKROOM, else "design-discussion".
# Workroom name selects per-room tmux session prefixes and per-room log
# directories so multiple workrooms (e.g. design-discussion +
# implementation-work) can boot side by side without colliding. The
# shared MCP server (`mcp`) and room.py UI (`room`) remain singletons.
WORKROOM_DEFAULT="${ACS_WORKROOM:-design-discussion}"
WORKROOM="${1:-$WORKROOM_DEFAULT}"
if ! [[ "$WORKROOM" =~ ^[a-z0-9][a-z0-9-]{0,31}$ ]]; then
  echo "[FAIL-LOUD] invalid workroom name '$WORKROOM' — must match ^[a-z0-9][a-z0-9-]{0,31}$ (lowercase alphanum + dash, 1-32 chars)." >&2
  exit 2
fi
export ACS_WORKROOM="$WORKROOM"

export COLLAB_MCP_LEDGER="${COLLAB_MCP_LEDGER:-$REPO_ROOT/docs/claude_codex_log.md}"
export COLLAB_ROOM_BIND="${COLLAB_ROOM_BIND:-127.0.0.1:7879}"
export COLLAB_ROOM_TOKEN="${COLLAB_ROOM_TOKEN:-acs-room}"
# Per-room relay targets so the relay knows which chat panes belong to
# this workroom. Override via COLLAB_TMUX_*_TARGET if you need to point
# the relay at a different physical pane.
export COLLAB_TMUX_CLAUDE_TARGET="${COLLAB_TMUX_CLAUDE_TARGET:-${WORKROOM}-claude-chat:0.0}"
export COLLAB_TMUX_CODEX_TARGET="${COLLAB_TMUX_CODEX_TARGET:-${WORKROOM}-codex-chat:0.0}"

LOG_DIR="/tmp/acs-collab/$WORKROOM"
SHARED_LOG_DIR="/tmp/acs-collab"
mkdir -p "$LOG_DIR" "$SHARED_LOG_DIR"

# Sentinels surface boot-time problems to the post-launch self-check
# (and to anyone tailing /tmp/acs-collab/<workroom>/) so a silent failure
# can never masquerade as a healthy boot the way it did on 2026-05-02
# (id=199). Sentinels live under the per-room dir so two workrooms do
# not overwrite each other's diagnostics.
SENTINEL_PYTHON="$LOG_DIR/.mcp_python_unhealthy"
SENTINEL_MCP="$LOG_DIR/.mcp_server_unhealthy"
SENTINEL_TOKEN="$LOG_DIR/.mcp_token_unset"
SENTINEL_SELFCHECK="$LOG_DIR/.selfcheck_failed"
rm -f "$SENTINEL_PYTHON" "$SENTINEL_MCP" "$SENTINEL_TOKEN" "$SENTINEL_SELFCHECK"

# Tmux session creation race-guard. Two launchers can both pass
# `tmux has-session` for a singleton (mcp / room) before either creates
# it; the loser then crashes `tmux new-session "duplicate session: X"`
# under `set -e` and kills the launcher. We serialize the create through
# a per-name lock file (shlock), and additionally swallow the
# duplicate-session error so a concurrent winner does not poison the
# loser. Idempotency invariant: after this function returns, the
# session exists or stderr explains why it could not be created.
SHARED_LOCK_DIR="$SHARED_LOG_DIR/locks"
mkdir -p "$SHARED_LOCK_DIR"
ensure_tmux_session_locked() {
  # $1=session_name $2=cwd $3=command (single string, will be passed to
  # `tmux new-session` as the command argument).
  local name="$1" cwd="$2" cmd="$3"
  if tmux has-session -t "$name" 2>/dev/null; then
    return 0
  fi
  local lock="$SHARED_LOCK_DIR/${name}.lock"
  # Use shlock when available (BSD/macOS); fall back to a best-effort
  # mkdir-based lock everywhere else. Both are advisory — the
  # has-session re-check below is the actual correctness guarantee.
  local locked=0
  if command -v shlock >/dev/null 2>&1; then
    if shlock -p $$ -f "$lock" >/dev/null 2>&1; then
      locked=1
    fi
  else
    if mkdir "${lock}.d" 2>/dev/null; then
      locked=1
    fi
  fi
  if [[ "$locked" -ne 1 ]]; then
    # Another launcher holds the lock. Spin up to ~3 s for them to
    # finish creating the session.
    local _i
    for _i in 1 2 3 4 5 6; do
      if tmux has-session -t "$name" 2>/dev/null; then
        return 0
      fi
      sleep 0.5
    done
    echo "[FAIL-LOUD] could not acquire tmux session lock for $name and session did not appear" >&2
    return 1
  fi
  # We hold the lock. Re-check, then create.
  if ! tmux has-session -t "$name" 2>/dev/null; then
    if ! tmux new-session -d -s "$name" -c "$cwd" "$cmd" 2>/tmp/.acs_tmux_err.$$; then
      local err
      err=$(cat /tmp/.acs_tmux_err.$$ 2>/dev/null || true)
      rm -f /tmp/.acs_tmux_err.$$
      if echo "$err" | grep -q -i "duplicate session"; then
        : # raced; the winner's session is fine
      else
        echo "[FAIL-LOUD] tmux new-session $name failed: $err" >&2
        if command -v shlock >/dev/null 2>&1; then
          rm -f "$lock"
        else
          rmdir "${lock}.d" 2>/dev/null || true
        fi
        return 1
      fi
    fi
    rm -f /tmp/.acs_tmux_err.$$
  fi
  if command -v shlock >/dev/null 2>&1; then
    rm -f "$lock"
  else
    rmdir "${lock}.d" 2>/dev/null || true
  fi
  if ! tmux has-session -t "$name" 2>/dev/null; then
    echo "[FAIL-LOUD] tmux session $name was not present after create attempt" >&2
    return 1
  fi
}

# Best-effort desktop notification. Always returns 0 so a missing
# osascript (e.g., headless ssh) never aborts the launcher.
notify() {
  local title="$1" subtitle="$2" message="$3"
  if command -v osascript >/dev/null 2>&1; then
    osascript -e "display notification \"$message\" with title \"$title\" subtitle \"$subtitle\"" \
      >/dev/null 2>&1 || true
  fi
}

# Python resolution is fail-loud: the boot path that silently fell back
# to `/usr/bin/python3` on 2026-05-02 left `mcp.log` reporting
# `ModuleNotFoundError: No module named 'mcp'` while launchd.err stayed
# empty. We refuse the system-python fallback and signal the failure on
# stderr (captured by launchd.err), via a desktop notification, and via
# a sentinel file so downstream steps and the self-check can react.
COLLAB_PYTHON="${COLLAB_MCP_PYTHON:-}"
if [[ -z "$COLLAB_PYTHON" && -x "$REPO_ROOT/.venv-collab/bin/python" ]]; then
  COLLAB_PYTHON="$REPO_ROOT/.venv-collab/bin/python"
fi

run_with_timeout() {
  local timeout_s="$1"
  shift
  "$@" &
  local pid=$!
  local waited=0
  while kill -0 "$pid" 2>/dev/null; do
    if (( waited >= timeout_s )); then
      kill "$pid" 2>/dev/null || true
      sleep 0.2
      kill -9 "$pid" 2>/dev/null || true
      wait "$pid" 2>/dev/null || true
      return 124
    fi
    sleep 1
    waited=$(( waited + 1 ))
  done
  wait "$pid"
}

PYTHON_HEALTHY=0
if [[ -n "$COLLAB_PYTHON" && -x "$COLLAB_PYTHON" ]]; then
  if run_with_timeout 5 "$COLLAB_PYTHON" -c 'import importlib.util, sys; sys.exit(0 if importlib.util.find_spec("mcp") else 1)' >/dev/null 2>&1; then
    PYTHON_HEALTHY=1
  fi
fi

if [[ "$PYTHON_HEALTHY" -ne 1 ]]; then
  echo "[FAIL-LOUD] acs-collab python unhealthy: COLLAB_MCP_PYTHON='${COLLAB_MCP_PYTHON:-}' resolved='${COLLAB_PYTHON:-<empty>}' — venv missing, 'mcp' module not discoverable, or health check timed out. Refusing system python3 fallback." >&2
  echo "[FAIL-LOUD] Fix: ensure $REPO_ROOT/.venv-collab/bin/python exists with 'pip install mcp[cli]>=1.2' OR set COLLAB_MCP_PYTHON to a python that has mcp." >&2
  printf 'unhealthy at %s; resolved=%s\n' "$(date -Iseconds 2>/dev/null || date)" "${COLLAB_PYTHON:-<empty>}" \
    >"$SENTINEL_PYTHON"
  notify "ACS Collab Workroom [$WORKROOM]" "MCP startup BLOCKED" "venv missing or mcp module not installed — see launchd.err"
fi

# 1) Syncthing — best-effort, never block.
if ! pgrep -x syncthing >/dev/null 2>&1; then
  if command -v brew >/dev/null 2>&1 && brew services list 2>/dev/null | grep -q '^syncthing'; then
    brew services start syncthing >/dev/null 2>&1 || true
  elif command -v syncthing >/dev/null 2>&1; then
    nohup syncthing --no-browser --no-restart >/dev/null 2>&1 &
  fi
fi

# 2) MCP server — must be alive before Claude/Codex panes start so their
#    MCP clients can complete the startup handshake. Skipped when the
#    python health gate above failed, so we never spawn a server.py that
#    is guaranteed to crash with ModuleNotFoundError. The server is a
#    shared singleton across workrooms (one DB, one process).
if [[ "$PYTHON_HEALTHY" -ne 1 ]]; then
  echo "[FAIL-LOUD] skipping MCP server step — python unhealthy (see $SENTINEL_PYTHON)" >&2
  printf 'skipped: python unhealthy\n' >"$SENTINEL_MCP"
elif [[ -z "${COLLAB_MCP_TOKEN:-}" ]]; then
  echo "[FAIL-LOUD] COLLAB_MCP_TOKEN unset in launcher environment — MCP server not started. Set it in the launchd plist (see com.activecellsim.workroom.plist EnvironmentVariables)." >&2
  printf 'token unset at boot\n' >"$SENTINEL_TOKEN"
  notify "ACS Collab Workroom [$WORKROOM]" "MCP startup BLOCKED" "COLLAB_MCP_TOKEN missing in launchd env — see launchd.err"
else
  ensure_tmux_session_locked "mcp" "$REPO_ROOT" \
    "'$COLLAB_PYTHON' tools/collab_mcp/server.py >>'$SHARED_LOG_DIR/mcp.log' 2>&1"
  MCP_UP=0
  for _i in 1 2 3 4 5 6 7 8 9 10; do
    if curl -sS --max-time 1 -o /dev/null "http://127.0.0.1:7878/mcp/" >/dev/null 2>&1; then
      MCP_UP=1
      break
    fi
    sleep 0.5
  done
  if [[ "$MCP_UP" -ne 1 ]]; then
    echo "[FAIL-LOUD] MCP server did not answer on http://127.0.0.1:7878/mcp/ within ~5 s — see $SHARED_LOG_DIR/mcp.log for the crash trace." >&2
    printf 'no response on 7878 within 5s\n' >"$SENTINEL_MCP"
    notify "ACS Collab Workroom [$WORKROOM]" "MCP startup BLOCKED" "server.py started but 7878 never came up — see mcp.log"
  fi
fi

# 3) tmux work sessions — workroom-prefixed:
#      <wr>-claude-chat, <wr>-codex-chat, <wr>-claude-work, <wr>-codex-work,
#      <wr>-win-ssh, <wr>-heartbeat.
#    setup_tmux_workroom is idempotent: never kills or renames existing
#    sessions, and only starts a CLI/SSH/daemon in panes that are
#    currently a plain shell.
python3 -m tools.collab_mcp.setup_tmux_workroom \
  --workroom "$WORKROOM" \
  --start-chat --start-work --start-win-ssh --start-heartbeat \
  >>"$LOG_DIR/setup_tmux_workroom.log" 2>&1 || true

# 4) Relay session — keeps PI/LLM messages flowing into chat panes for
#    THIS workroom. Each workroom gets its own relay so its DB cursor
#    (relay_cursors row) does not collide with another workroom's relay.
RELAY_SESSION="${WORKROOM}-relay"
ensure_tmux_session_locked "$RELAY_SESSION" "$REPO_ROOT" \
  "python3 tools/collab_mcp/tmux_relay.py \
     --relay-name '$WORKROOM' \
     --workroom '$WORKROOM' \
     --claude-target '${WORKROOM}-claude-chat:0.0' \
     --codex-target '${WORKROOM}-codex-chat:0.0' \
     --include-llm-messages \
     --init-cursor-if-missing \
     >>'$LOG_DIR/tmux_relay.log' 2>&1"

# 5) room.py — browser UI server (shared singleton, port 7879).
if ! curl -fsS "http://127.0.0.1:7879/" >/dev/null 2>&1; then
  ensure_tmux_session_locked "room" "$REPO_ROOT" \
    "python3 tools/collab_mcp/room.py >>'$SHARED_LOG_DIR/room.log' 2>&1"
  for _i in 1 2 3 4 5 6 7 8 9 10; do
    if curl -fsS "http://127.0.0.1:7879/" >/dev/null 2>&1; then
      break
    fi
    sleep 0.5
  done
fi

# 6) Chrome --app mode (fallback to default browser).
#    Idempotent across two failure modes that the original probe missed:
#    (a) Chrome not yet running so the AppleScript URL list returns 0
#        and a concurrent kickstart races into a second `open -na`.
#    (b) The first `open -na` returns before the new tab is queryable,
#        so the second invocation also sees 0 matching URLs.
#    A short-lived lock file (CHROME_LOCK_TTL_S) bridges both windows.
#    The lock is shared across workrooms because the UI is one tab.
CHROME_LOCK="$SHARED_LOG_DIR/chrome_app.lock"
CHROME_LOCK_TTL_S=30
chrome_lock_fresh() {
  [[ -f "$CHROME_LOCK" ]] || return 1
  local now lock_mtime age
  now=$(date +%s)
  lock_mtime=$(stat -f '%m' "$CHROME_LOCK" 2>/dev/null || echo 0)
  age=$((now - lock_mtime))
  (( age >= 0 && age < CHROME_LOCK_TTL_S ))
}
if chrome_lock_fresh; then
  : # another launcher run opened Chrome within the last $CHROME_LOCK_TTL_S — skip
elif open -Ra "Google Chrome" >/dev/null 2>&1; then
  EXISTING_TABS=$(osascript \
    -e 'tell application "Google Chrome" to get URL of every tab of every window' \
    2>/dev/null | tr ',' '\n' | grep -c '127.0.0.1:7879' || true)
  if [[ "${EXISTING_TABS:-0}" -eq 0 ]]; then
    : >"$CHROME_LOCK"
    open -na "Google Chrome" --args --app="http://127.0.0.1:7879/"
  fi
else
  : >"$CHROME_LOCK"
  open "http://127.0.0.1:7879/"
fi

# 7) Post-launch self-check. Re-verifies every tier the launcher just
#    touched and surfaces failures via stderr (captured to launchd.err),
#    a sentinel file, and a desktop notification. Without this step a
#    crashed MCP server step (cf. 2026-05-02 boot) leaves launchd.err
#    empty and the only signal lives buried inside mcp.log.
SELFCHECK_FAILS=()
SELFCHECK_OK=()

selfcheck_session() {
  local name="$1"
  if tmux has-session -t "$name" 2>/dev/null; then
    SELFCHECK_OK+=("tmux:$name")
  else
    SELFCHECK_FAILS+=("tmux:$name MISSING")
  fi
}

selfcheck_http() {
  local label="$1" url="$2"
  if curl -sS --max-time 2 -o /dev/null "$url" >/dev/null 2>&1; then
    SELFCHECK_OK+=("http:$label")
  else
    SELFCHECK_FAILS+=("http:$label DOWN ($url)")
  fi
}

# Per-room sessions
for s in claude-chat claude-work codex-chat codex-work win-ssh relay heartbeat; do
  selfcheck_session "${WORKROOM}-${s}"
done
# Shared singletons
for s in room mcp; do
  selfcheck_session "$s"
done
selfcheck_http "room.py" "http://127.0.0.1:7879/"
selfcheck_http "mcp"     "http://127.0.0.1:7878/mcp/"

if (( ${#SELFCHECK_FAILS[@]} > 0 )); then
  echo "[SELF-CHECK FAIL] ${#SELFCHECK_FAILS[@]} item(s): ${SELFCHECK_FAILS[*]}" >&2
  echo "[SELF-CHECK OK ] ${#SELFCHECK_OK[@]} item(s): ${SELFCHECK_OK[*]}" >&2
  printf 'failed at %s\nworkroom: %s\nFAIL: %s\nOK: %s\n' \
    "$(date -Iseconds 2>/dev/null || date)" \
    "$WORKROOM" \
    "${SELFCHECK_FAILS[*]}" \
    "${SELFCHECK_OK[*]}" \
    >"$SENTINEL_SELFCHECK"
  notify "ACS Collab Workroom [$WORKROOM]" "Self-check FAILED (${#SELFCHECK_FAILS[@]})" "$(printf '%s' "${SELFCHECK_FAILS[*]}" | head -c 200)"
else
  echo "[SELF-CHECK OK ] workroom=$WORKROOM all ${#SELFCHECK_OK[@]} tiers up" >&2
fi
