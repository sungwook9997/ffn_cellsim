#!/usr/bin/env bash
#
# sync_to_gbook.sh — one-command code-sync of ffn_sim/ to the gbook GPU box + import smoke.
#
# CLAUDE.md note: gbook code is NOT auto-synced (Syncthing carries outputs, not the working
# tree's code in a guaranteed-current state for a run). A stale-code gbook run produces
# results indistinguishable from current code (which is exactly why layer2_gpu_scaleup stamps
# the git commit into every JSONL). This script makes "push my code to gbook, then prove it
# imports there" a single repeatable step before launching a GPU production sweep.
#
# What it does:
#   1. rsync the ffn_sim/ package to gbook:~/ffn_cellsim/ffn_sim/ (code only — outputs,
#      __pycache__, *.pyc, .git, .claude are excluded so we never ship bulk data / caches /
#      VCS state that Syncthing or git own).
#   2. ssh gbook and run a one-line import smoke against the gbook conda env, confirming the
#      synced code actually imports there (catches a missing dep / syntax error before a run).
#
# Run from the repo root:  bash ffn_sim/scripts/sync_to_gbook.sh
#
set -euo pipefail

# --- config (no hard-coded device IDs; the host + env path are gbook infra, per memory
#     reference-gbook-workstation-specs: Ubuntu, miniconda3 env ffn_sim). -------------------
GBOOK_HOST="gbook"
GBOOK_REPO="~/ffn_cellsim"
GBOOK_PY="~/miniconda3/envs/ffn_sim/bin/python"

# Resolve the LOCAL ffn_sim/ source dir relative to this script (repo-root-runnable: the
# script lives in ffn_sim/scripts/, so its package root is one level up). Trailing slash on
# the source = "copy the contents of ffn_sim/ into the remote ffn_sim/" (rsync semantics).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_PKG="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "[sync_to_gbook] rsync ${LOCAL_PKG}/ -> ${GBOOK_HOST}:${GBOOK_REPO}/ffn_sim/"
if rsync -a \
    --exclude outputs \
    --exclude __pycache__ \
    --exclude '*.pyc' \
    --exclude .git \
    --exclude .claude \
    "${LOCAL_PKG}/" "${GBOOK_HOST}:${GBOOK_REPO}/ffn_sim/"; then
  echo "[sync_to_gbook] rsync OK"
else
  echo "[sync_to_gbook] FAIL — rsync to ${GBOOK_HOST} failed" >&2
  exit 1
fi

# Import smoke on gbook: PYTHONPATH=repo-root so `import ffn_sim.*` resolves; import the
# native-N GPU driver (the thing a sweep launches) to prove the synced tree is importable.
echo "[sync_to_gbook] import smoke on ${GBOOK_HOST} (layer2_gpu_scaleup) ..."
if ssh "${GBOOK_HOST}" \
    "PYTHONPATH=${GBOOK_REPO} ${GBOOK_PY} -c 'import ffn_sim.scripts.layer2_gpu_scaleup; print(\"import OK\")'"; then
  echo "[sync_to_gbook] OK — code synced and imports on ${GBOOK_HOST}"
else
  echo "[sync_to_gbook] FAIL — synced code does not import on ${GBOOK_HOST}" >&2
  exit 1
fi
