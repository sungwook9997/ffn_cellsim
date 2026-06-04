#!/usr/bin/env bash
# Refresh the Obsidian RAG mirror from Notion (one command).
# Notion = source of truth; this regenerates the vault from scratch.
#   usage:  bash refresh.sh [--open]
# Requires: conda env ffn_sim, and .notion_token next to this script
# (Internal Integration secret connected to the Contract Graph page).
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -f .notion_token ]; then
  echo "ERROR: .notion_token missing. Create it:"
  echo "  echo 'NOTION_TOKEN=ntn_...' > $(pwd)/.notion_token"
  exit 1
fi

# activate env (works whether or not already active)
source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate ffn_sim 2>/dev/null || true

echo "[1/4] Notion -> vault (claims, papers, params, gates, contracts, relations)"
python notion_to_obsidian.py

echo "[2/4] repo -> code/test/doc nodes linked to KB claims"
python add_code_nodes.py

echo "[3/4] Notion Dev Logs board -> day-log nodes cross-linked into the graph"
python devlogs_to_obsidian.py

echo "[4/4] supersession / authoritative-chain edges (docs front-matter)"
python supersession_to_obsidian.py

echo
echo "[sanity] Obsidian mirror summary"
python - <<'PY'
from pathlib import Path

vault = Path("vault")
notes = list(vault.glob("*.md"))
devlogs = list(vault.glob("DL_*.md"))
index = vault / "00_INDEX.md"

print(f"  vault notes: {len(notes)}")
print(f"  dev-log notes: {len(devlogs)}")
print(f"  index present: {'yes' if index.exists() else 'no'}")
PY

echo "done. vault = $(pwd)/vault"
if [ "${1:-}" = "--open" ]; then
  open -a Obsidian "$(pwd)/vault"
fi
