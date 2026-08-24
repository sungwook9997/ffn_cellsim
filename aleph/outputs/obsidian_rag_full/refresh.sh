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

# `python` is NOT on PATH in this project's environments and `conda activate` is a no-op in a
# non-interactive shell, so every `python` line below died with "command not found" while CLAUDE.md
# tells sessions to run this script. Same defect, same fix, as the Makefile (2026-07-29) and
# tag_kb/refresh.sh — this copy was simply missed. The predicate is "can this interpreter RUN the
# pipeline", not "does a python exist": resolving /usr/bin/python3 would trade a PATH error for a
# ModuleNotFoundError. Override with `PY=/path/to/python bash refresh.sh`.
# The two `conda` lines that stood here are DELETED, not repaired. Under `set -e` a `source` of a
# missing file exits a non-interactive shell IMMEDIATELY, and `|| true` does not catch it — so with
# conda off PATH this script died before its first echo, silently, exit 1 and no output. They are
# also redundant: resolving the env interpreter below is what the Makefile and tag_kb/refresh.sh do,
# and it needs no activation.
if [ -z "${PY:-}" ]; then
  for c in "${CONDA_PREFIX:-}/bin/python" "$HOME/miniconda3/envs/ffn_sim/bin/python" \
           "$(command -v python 2>/dev/null || true)" "$(command -v python3 2>/dev/null || true)"; do
    [ -n "$c" ] && [ -x "$c" ] && "$c" -c 'import requests' >/dev/null 2>&1 && { PY="$c"; break; }
  done
fi
if [ -z "${PY:-}" ]; then
  echo "no interpreter on this machine can import requests — activate the env (conda activate ffn_sim)" >&2
  echo "or pass one explicitly: PY=/path/to/python bash refresh.sh" >&2
  exit 3
fi
echo "[env] interpreter: $PY"

echo "[1/5] Notion -> vault (claims, papers, params, gates, contracts, relations)"
"$PY" notion_to_obsidian.py

echo "[2/5] repo -> code/test/doc nodes linked to KB claims"
"$PY" add_code_nodes.py

echo "[3/5] Notion Dev Logs board -> day-log nodes cross-linked into the graph"
"$PY" devlogs_to_obsidian.py

echo "[4/5] supersession / authoritative-chain edges (docs front-matter)"
"$PY" supersession_to_obsidian.py

echo
echo "[5/5] vault entry point -> 00_INDEX.md"
"$PY" build_index.py

echo "[sanity] Obsidian mirror summary"
"$PY" - <<'PY'
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
