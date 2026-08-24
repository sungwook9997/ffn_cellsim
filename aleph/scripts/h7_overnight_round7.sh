#!/usr/bin/env bash
# H.7 round 7 (gbook GPU): connected-mesh CONTROL at full ×40 scale — does a
# PERCOLATED cortex (z~3.3) change gamma vs the fragmented random-anchor mesh?
# (Local smoke said no: floor holds. This confirms at production scale.)
set -u
cd "$(dirname "$0")/../.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
conda activate ffn_sim
OUT=aleph/outputs/h7/production/round7; mkdir -p "$OUT"; ST="$OUT/round7.status"
echo "START $(date -Is)" > "$ST"
for s in 1 2 3; do
  echo "CM_RIGID_SEED_${s}_START $(date -Is)" >> "$ST"
  python -m aleph.scripts.h7_gate_b --device gpu --no-fa --connected-mesh --warmup 2000 --sample 1500 --seed "$s" > "$OUT/cm_rigid_seed${s}.log" 2>&1
  echo "CM_RIGID_SEED_${s}_EXIT $? $(date -Is)" >> "$ST"
done
for s in 1 2 3; do
  echo "CM_SOFT_SEED_${s}_START $(date -Is)" >> "$ST"
  python -m aleph.scripts.h7_gate_b --device gpu --no-constrained --connected-mesh --warmup 2000 --sample 1500 --seed "$s" > "$OUT/cm_soft_seed${s}.log" 2>&1
  echo "CM_SOFT_SEED_${s}_EXIT $? $(date -Is)" >> "$ST"
done
echo "ALL_DONE $(date -Is)" >> "$ST"
