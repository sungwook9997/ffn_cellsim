#!/usr/bin/env bash
# gbook DCM production launcher — fire the N≥400 GPU runs the moment gbook is online.
#
# DCM autonomous loop (2026-06-29). gbook (RTX A5000) is the only production GPU; the local Mac
# Warp build has no CUDA. Syncthing does NOT sync source code (outputs/ only), so code must be
# rsync'd before launching. This script: (1) checks gbook is reachable, (2) rsyncs the ffn_sim
# package (code only), (3) launches the two production runs under setsid + logs, (4) prints
# monitor commands. Runs share the GPU with the FF session (PI-sanctioned: "gpu 같이 써줘").
#
# The runs (validated locally at small N; here scaled to N≥400 / cuda:0):
#   1. CONFLUENT FACETING SWEEP — tests the iter-4 conclusion that faceting needs confinement:
#      at N≥400 interior cells are confined by neighbours, so cortical tension should flatten
#      junctions (faceting) instead of freely deflating the cell. K=2500 (SimuCell3D faceting),
#      cohesion 5e7 (lit MCF7), γ ladder crossing γ̃∈[0.02,0.10]. γ = controlled variable.
#   2. DIVISION PRODUCTION — proliferation separation + mesh integrity at production resolution
#      (subdiv 2 resolves the contact the subdiv-1 CPU smoke could not). Mitotic division at the
#      physiological rate; driver-default K (turgor holds volume so daughters inflate). NO remesh
#      (division+remesh co-run is deferred — see DCM_DIVISION_REMESH_CORUN_DESIGN_2026-06-29.md).
#
# Usage:
#   bash aleph/dcm/gbook_production.sh --dry-run     # print the exact commands, run nothing
#   bash aleph/dcm/gbook_production.sh               # rsync + launch on gbook
#
# CONFIRM these against the live gbook before the first real launch (gbook was OFFLINE at authoring):
set -euo pipefail

GBOOK_HOST="${GBOOK_HOST:-gbook}"
REMOTE_DIR="${REMOTE_DIR:-/home/sungwook/ffn_dcm_run}"   # fresh dir (gbook's repos are pre-restructure); rsync target
CONDA_ENV="${CONDA_ENV:-ffn_sim}"
CONDA_SH="${CONDA_SH:-/home/sungwook/miniconda3/etc/profile.d/conda.sh}"
LOCAL_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"        # the worktree root (…/dcm-aggregation)
OUTBASE="aleph/outputs/h_dcm_two_stage"
DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

# --- production commands (relative to REMOTE_DIR; PYTHONPATH=. ) -----------------------------
GAMMAS="0,5e-4,1e-3,2e-3,3e-3,5e-3,8e-3"
FACET_CMD="PYTHONPATH=. python -m aleph.dcm.gamma_sweep --device cuda:0 \
  --n-cells 400 --subdiv 2 --steps 40000 --frames 12 --k-vol 2500 --adh-strength 5e7 \
  --gap 2.1 --gammas ${GAMMAS} --jobs 2 --out ${OUTBASE}/gbook_confluent_facet_N400"
DIV_CMD="PYTHONPATH=. python -m aleph.dcm.dcm_warp_decohesion --device cuda:0 \
  --n-cells 400 --subdiv 2 --steps 60000 --frames 20 --division --div-rate 0.04 \
  --builder fcc --gap 2.1 --save-frames ${OUTBASE}/gbook_division_N400.npz"

run_remote () {  # $1 = description, $2 = command, $3 = remote log path
  echo "### $1"
  echo "$2"
  echo "   log: ${REMOTE_DIR}/$3"
  if [[ $DRY -eq 0 ]]; then
    ssh "$GBOOK_HOST" "cd ${REMOTE_DIR} && source ~/.bashrc && conda activate ${CONDA_ENV} && \
      setsid bash -c '$2 > $3 2>&1' &" || true
  fi
}

echo "== gbook DCM production launcher (DRY=${DRY}) =="
echo "host=${GBOOK_HOST} remote=${REMOTE_DIR} env=${CONDA_ENV} local=${LOCAL_ROOT}"

if [[ $DRY -eq 0 ]]; then
  echo "-- reachability --"
  ssh -o ConnectTimeout=6 -o BatchMode=yes "$GBOOK_HOST" 'echo ALIVE; nvidia-smi \
    --query-gpu=name,memory.used,memory.total --format=csv,noheader' \
    || { echo "gbook UNREACHABLE — aborting"; exit 1; }
  echo "-- rsync code (ffn_sim package; exclude outputs data / caches / git) --"
  rsync -az --delete \
    --exclude 'outputs/' --exclude '__pycache__/' --exclude '.git/' --exclude '*.npz' \
    "${LOCAL_ROOT}/aleph/" "${GBOOK_HOST}:${REMOTE_DIR}/aleph/"
fi

echo
echo "== production runs =="
run_remote "1. CONFLUENT FACETING SWEEP (N=400, K=2500, faceting=confinement test)" "$FACET_CMD" \
  "${OUTBASE}/gbook_confluent_facet_N400.prod.log"
echo
run_remote "2. DIVISION PRODUCTION (N=400, subdiv2, separation + mesh integrity)" "$DIV_CMD" \
  "${OUTBASE}/gbook_division_N400.prod.log"

echo
echo "== monitor (after launch) =="
echo "  ssh ${GBOOK_HOST} 'tail -f ${REMOTE_DIR}/${OUTBASE}/gbook_confluent_facet_N400/*.log'"
echo "  ssh ${GBOOK_HOST} 'nvidia-smi'"
echo "  # results return via Syncthing (outputs/ is synced); render montages locally."
[[ $DRY -eq 1 ]] && echo "(dry-run: nothing executed)"
