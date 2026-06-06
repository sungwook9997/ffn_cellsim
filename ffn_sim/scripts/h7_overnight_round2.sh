#!/usr/bin/env bash
# H.7 overnight production ROUND 2 (gbook GPU). Addresses the round-1 +
# direction/novelty-review findings:
#   - GATE-B 5-seed ensemble (round-1 budget): gamma_rigid reproducibility +
#     is the gamma_active floor robust across seeds?
#   - GATE-B 5x-longer single run: does gamma_active rise with run length, or is
#     it truly the generation floor (us-warmup << seconds myosin timescale)?
#   - Spreading with per-chunk frame dump (both geometries, 12 chunks): data for
#     the cell-movement animation (h7_spreading_anim.py).
#
# Usage (gbook):  setsid bash ffn_sim/scripts/h7_overnight_round2.sh &
set -u
cd "$(dirname "$0")/../.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
conda activate ffn_sim

OUT=ffn_sim/outputs/h7/production/round2
mkdir -p "$OUT" ffn_sim/outputs/h7/figs
ST="$OUT/round2.status"
echo "START $(date -Is)" > "$ST"

for s in 1 2 3 4 5; do
  echo "ENS_SEED_${s}_START $(date -Is)" >> "$ST"
  python -m ffn_sim.scripts.h7_gate_b --device gpu --warmup 2000 --sample 1500 --seed "$s" \
    > "$OUT/gate_b_seed${s}.log" 2>&1
  echo "ENS_SEED_${s}_EXIT $? $(date -Is)" >> "$ST"
done

echo "LONG_START $(date -Is)" >> "$ST"
python -m ffn_sim.scripts.h7_gate_b --device gpu --warmup 10000 --sample 5000 --seed 1 \
  > "$OUT/gate_b_long.log" 2>&1
echo "LONG_EXIT $? $(date -Is)" >> "$ST"

echo "SPREAD_START $(date -Is)" >> "$ST"
python -m ffn_sim.scripts.h7_spreading_compare --device gpu \
  --warmup 1500 --chunks 12 --chunk-steps 2000 \
  --out ffn_sim/outputs/h7/figs/h7_spreading_compare_round2.png \
  --frames-out "$OUT/spreading_frames.npz" \
  > "$OUT/spreading_round2.log" 2>&1
echo "SPREAD_EXIT $? $(date -Is)" >> "$ST"

echo "ALL_DONE $(date -Is)" >> "$ST"
