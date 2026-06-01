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

echo "[1/2] Notion -> vault (claims, papers, params, gates, contracts, relations)"
python notion_to_obsidian.py

echo "[2/2] repo -> code/test/doc nodes linked to KB claims"
python add_code_nodes.py

echo "done. vault = $(pwd)/vault"
if [ "${1:-}" = "--open" ]; then
  open -a Obsidian "$(pwd)/vault"
fi
