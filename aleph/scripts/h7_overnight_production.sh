#!/usr/bin/env bash
# H.7 overnight production launcher (gbook GPU).
#
# Runs, sequentially (single GPU), logging each + writing status markers a remote
# poller can read:
#   1. GATE-B  — authoritative emergent cortical tension (3 channels separate).
#   2. A/A0    — both single-cell lamellipodium spreading geometries (rim+patch).
#
# Usage (on gbook):  setsid bash aleph/scripts/h7_overnight_production.sh &
set -u
cd "$(dirname "$0")/../.." || exit 1            # repo root (~/ffn_cellsim)
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
conda activate ffn_sim

OUT=aleph/outputs/h7/production
mkdir -p "$OUT"
ST="$OUT/overnight.status"
echo "START $(date -Is)" > "$ST"

echo "GATE_B_START $(date -Is)" >> "$ST"
python -m aleph.scripts.h7_gate_b --device gpu --warmup 2000 --sample 1500 \
  > "$OUT/gate_b_full.log" 2>&1
echo "GATE_B_EXIT $? $(date -Is)" >> "$ST"

echo "SPREADING_START $(date -Is)" >> "$ST"
python -m aleph.scripts.h7_spreading_compare --device gpu \
  --warmup 1500 --chunks 6 --chunk-steps 1500 \
  --out "$OUT/../figs/h7_spreading_compare_full.png" \
  > "$OUT/spreading_full.log" 2>&1
echo "SPREADING_EXIT $? $(date -Is)" >> "$ST"

echo "ALL_DONE $(date -Is)" >> "$ST"
