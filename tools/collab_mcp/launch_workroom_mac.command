#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="${0:A:h}"
REPO_ROOT="${SCRIPT_DIR:h:h}"
cd "$REPO_ROOT"

export COLLAB_MCP_LEDGER="${COLLAB_MCP_LEDGER:-$REPO_ROOT/docs/claude_codex_log.md}"
export COLLAB_ROOM_BIND="${COLLAB_ROOM_BIND:-127.0.0.1:7879}"
export COLLAB_ROOM_TOKEN="${COLLAB_ROOM_TOKEN:-acs-room}"

if ! curl -fsS "http://127.0.0.1:7879/" >/dev/null 2>&1; then
  mkdir -p /tmp/acs-collab
  nohup python3 tools/collab_mcp/room.py >/tmp/acs-collab/room.log 2>&1 &
  sleep 1
fi

if open -Ra "Google Chrome" >/dev/null 2>&1; then
  open -na "Google Chrome" --args --app="http://127.0.0.1:7879/"
else
  open "http://127.0.0.1:7879/"
fi
