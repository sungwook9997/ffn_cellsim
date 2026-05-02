#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="${0:A:h}"
REPO_ROOT="${SCRIPT_DIR:h:h}"
cd "$REPO_ROOT"

export COLLAB_MCP_LEDGER="${COLLAB_MCP_LEDGER:-$REPO_ROOT/docs/claude_codex_log.md}"
export COLLAB_ROOM_BIND="${COLLAB_ROOM_BIND:-127.0.0.1:7879}"
export COLLAB_ROOM_TOKEN="${COLLAB_ROOM_TOKEN:-acs-room}"
export COLLAB_TMUX_CLAUDE_TARGET="${COLLAB_TMUX_CLAUDE_TARGET:-claude-chat:0.0}"
export COLLAB_TMUX_CODEX_TARGET="${COLLAB_TMUX_CODEX_TARGET:-codex-chat:0.0}"

LOG_DIR="/tmp/acs-collab"
mkdir -p "$LOG_DIR"

# Sentinels surface boot-time problems to the post-launch self-check
# (and to anyone tailing /tmp/acs-collab/) so a silent failure can never
# masquerade as a healthy boot the way it did on 2026-05-02 (id=199).
SENTINEL_PYTHON="$LOG_DIR/.mcp_python_unhealthy"
SENTINEL_MCP="$LOG_DIR/.mcp_server_unhealthy"
SENTINEL_TOKEN="$LOG_DIR/.mcp_token_unset"
SENTINEL_SELFCHECK="$LOG_DIR/.selfcheck_failed"
rm -f "$SENTINEL_PYTHON" "$SENTINEL_MCP" "$SENTINEL_TOKEN" "$SENTINEL_SELFCHECK"

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
  notify "ACS Collab Workroom" "MCP startup BLOCKED" "venv missing or mcp module not installed — see launchd.err"
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
#    is guaranteed to crash with ModuleNotFoundError.
if [[ "$PYTHON_HEALTHY" -ne 1 ]]; then
  echo "[FAIL-LOUD] skipping MCP server step — python unhealthy (see $SENTINEL_PYTHON)" >&2
  printf 'skipped: python unhealthy\n' >"$SENTINEL_MCP"
elif [[ -z "${COLLAB_MCP_TOKEN:-}" ]]; then
  echo "[FAIL-LOUD] COLLAB_MCP_TOKEN unset in launcher environment — MCP server not started. Set it in the launchd plist (see com.activecellsim.workroom.plist EnvironmentVariables)." >&2
  printf 'token unset at boot\n' >"$SENTINEL_TOKEN"
  notify "ACS Collab Workroom" "MCP startup BLOCKED" "COLLAB_MCP_TOKEN missing in launchd env — see launchd.err"
else
  if ! tmux has-session -t mcp 2>/dev/null; then
    tmux new-session -d -s mcp -c "$REPO_ROOT" \
      "'$COLLAB_PYTHON' tools/collab_mcp/server.py >>'$LOG_DIR/mcp.log' 2>&1"
  fi
  MCP_UP=0
  for _i in 1 2 3 4 5 6 7 8 9 10; do
    if curl -sS --max-time 1 -o /dev/null "http://127.0.0.1:7878/mcp/" >/dev/null 2>&1; then
      MCP_UP=1
      break
    fi
    sleep 0.5
  done
  if [[ "$MCP_UP" -ne 1 ]]; then
    echo "[FAIL-LOUD] MCP server did not answer on http://127.0.0.1:7878/mcp/ within ~5 s — see $LOG_DIR/mcp.log for the crash trace." >&2
    printf 'no response on 7878 within 5s\n' >"$SENTINEL_MCP"
    notify "ACS Collab Workroom" "MCP startup BLOCKED" "server.py started but 7878 never came up — see mcp.log"
  fi
fi

# 3) tmux work sessions (claude-chat, codex-chat, claude-work, codex-work,
#    win-ssh, heartbeat). setup_tmux_workroom is idempotent: never kills or
#    renames existing sessions, and only starts a CLI/SSH/daemon in panes
#    that are currently a plain shell.
python3 -m tools.collab_mcp.setup_tmux_workroom \
  --start-chat --start-work --start-win-ssh --start-heartbeat \
  >>"$LOG_DIR/setup_tmux_workroom.log" 2>&1 || true

# 4) Relay session — keeps PI/LLM messages flowing into chat panes.
if ! tmux has-session -t relay 2>/dev/null; then
  tmux new-session -d -s relay -c "$REPO_ROOT" \
    "python3 tools/collab_mcp/tmux_relay.py >>'$LOG_DIR/tmux_relay.log' 2>&1"
fi

# 5) room.py — browser UI server.
if ! curl -fsS "http://127.0.0.1:7879/" >/dev/null 2>&1; then
  if ! tmux has-session -t room 2>/dev/null; then
    tmux new-session -d -s room -c "$REPO_ROOT" \
      "python3 tools/collab_mcp/room.py >>'$LOG_DIR/room.log' 2>&1"
  fi
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
CHROME_LOCK="$LOG_DIR/chrome_app.lock"
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

for s in claude-chat claude-work codex-chat codex-work win-ssh relay heartbeat room mcp; do
  selfcheck_session "$s"
done
selfcheck_http "room.py" "http://127.0.0.1:7879/"
selfcheck_http "mcp"     "http://127.0.0.1:7878/mcp/"

if (( ${#SELFCHECK_FAILS[@]} > 0 )); then
  echo "[SELF-CHECK FAIL] ${#SELFCHECK_FAILS[@]} item(s): ${SELFCHECK_FAILS[*]}" >&2
  echo "[SELF-CHECK OK ] ${#SELFCHECK_OK[@]} item(s): ${SELFCHECK_OK[*]}" >&2
  printf 'failed at %s\nFAIL: %s\nOK: %s\n' \
    "$(date -Iseconds 2>/dev/null || date)" \
    "${SELFCHECK_FAILS[*]}" \
    "${SELFCHECK_OK[*]}" \
    >"$SENTINEL_SELFCHECK"
  notify "ACS Collab Workroom" "Self-check FAILED (${#SELFCHECK_FAILS[@]})" "$(printf '%s' "${SELFCHECK_FAILS[*]}" | head -c 200)"
else
  echo "[SELF-CHECK OK ] all ${#SELFCHECK_OK[@]} tiers up" >&2
fi
