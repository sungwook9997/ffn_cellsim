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

# 1) Syncthing — best-effort, never block.
if ! pgrep -x syncthing >/dev/null 2>&1; then
  if command -v brew >/dev/null 2>&1 && brew services list 2>/dev/null | grep -q '^syncthing'; then
    brew services start syncthing >/dev/null 2>&1 || true
  elif command -v syncthing >/dev/null 2>&1; then
    nohup syncthing --no-browser --no-restart >/dev/null 2>&1 &
  fi
fi

# 2) tmux work sessions (claude-chat, codex-chat, claude-work, codex-work, win-ssh).
#    setup_tmux_workroom is idempotent: never kills/renames existing sessions
#    and only starts a CLI/SSH command in panes that are currently a plain shell.
python3 -m tools.collab_mcp.setup_tmux_workroom --start-chat --start-work --start-win-ssh \
  >>"$LOG_DIR/setup_tmux_workroom.log" 2>&1 || true

# 3) Relay session — keeps PI/LLM messages flowing into chat panes.
if ! tmux has-session -t relay 2>/dev/null; then
  tmux new-session -d -s relay -c "$REPO_ROOT" \
    "python3 tools/collab_mcp/tmux_relay.py >>'$LOG_DIR/tmux_relay.log' 2>&1"
fi

# 4) room.py — browser UI server.
if ! curl -fsS "http://127.0.0.1:7879/" >/dev/null 2>&1; then
  nohup python3 tools/collab_mcp/room.py >>"$LOG_DIR/room.log" 2>&1 &
  for _i in 1 2 3 4 5 6 7 8 9 10; do
    if curl -fsS "http://127.0.0.1:7879/" >/dev/null 2>&1; then
      break
    fi
    sleep 0.5
  done
fi

# 5) Chrome --app mode (fallback to default browser).
if open -Ra "Google Chrome" >/dev/null 2>&1; then
  open -na "Google Chrome" --args --app="http://127.0.0.1:7879/"
else
  open "http://127.0.0.1:7879/"
fi
