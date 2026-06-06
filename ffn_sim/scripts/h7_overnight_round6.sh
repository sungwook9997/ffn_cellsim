#!/usr/bin/env bash
# H.7 overnight production ROUND 6 (gbook GPU). The findings are consolidated;
# round 6 is the ×40 coarse-graining CONVERGENCE study on the STABLE paths (the
# novelty-defending grid-convergence the constrained+FA path couldn't run):
#   - constrained --no-fa (free pressurised cortex) gamma_rigid at n_filaments
#     {500, 1000, 2000}: is the structural channel grid-converged?
#   - unconstrained-soft gamma_soft at {500, 2000} (1000 = round-4): is the active
#     floor scale-independent?
#
# Usage (gbook):  setsid bash ffn_sim/scripts/h7_overnight_round6.sh &
set -u
cd "$(dirname "$0")/../.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
conda activate ffn_sim

OUT=ffn_sim/outputs/h7/production/round6
mkdir -p "$OUT"
ST="$OUT/round6.status"
echo "START $(date -Is)" > "$ST"

for nf in 500 1000 2000; do
  echo "RIGID_NF_${nf}_START $(date -Is)" >> "$ST"
  python -m ffn_sim.scripts.h7_gate_b --device gpu --no-fa --n-filaments "$nf" \
    --warmup 2000 --sample 1500 --seed 1 \
    > "$OUT/gate_b_rigid_nf${nf}.log" 2>&1
  echo "RIGID_NF_${nf}_EXIT $? $(date -Is)" >> "$ST"
done

for nf in 500 2000; do
  echo "SOFT_NF_${nf}_START $(date -Is)" >> "$ST"
  python -m ffn_sim.scripts.h7_gate_b --device gpu --no-constrained --n-filaments "$nf" \
    --warmup 2000 --sample 1500 --seed 1 \
    > "$OUT/gate_b_soft_nf${nf}.log" 2>&1
  echo "SOFT_NF_${nf}_EXIT $? $(date -Is)" >> "$ST"
done

echo "ALL_DONE $(date -Is)" >> "$ST"
