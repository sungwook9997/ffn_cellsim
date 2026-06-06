#!/usr/bin/env bash
# H.7 overnight production ROUND 5 (gbook GPU). Round-4 gave a clean active-floor
# ensemble (gamma_soft = 0.0001 mN/m, n=5). Round 5 tests the two remaining
# hypotheses about the active floor + gets a clean structural channel:
#   - TURNOVER-ON unconstrained-soft GATE-B (3 seeds, longer warmup): does
#     Chugh-2017 actin turnover un-floor the active channel? (Expected NO: the
#     turnover timescale tau_half=10 s is also >> the feasible us-dt run, so the
#     active floor is a TIMESCALE GAP, not a missing mechanism — this confirms it.)
#   - clean gamma_rigid: constrained, NO FA (free pressurized cortex; no integrin
#     overload) -> the structural backbone tension at the turgor operating point.
#
# Usage (gbook):  setsid bash ffn_sim/scripts/h7_overnight_round5.sh &
set -u
cd "$(dirname "$0")/../.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
conda activate ffn_sim

OUT=ffn_sim/outputs/h7/production/round5
mkdir -p "$OUT"
ST="$OUT/round5.status"
echo "START $(date -Is)" > "$ST"

for s in 1 2 3; do
  echo "TURNOVER_SEED_${s}_START $(date -Is)" >> "$ST"
  python -m ffn_sim.scripts.h7_gate_b --device gpu --no-constrained --turnover \
    --warmup 5000 --sample 3000 --seed "$s" \
    > "$OUT/gate_b_turnover_seed${s}.log" 2>&1
  echo "TURNOVER_SEED_${s}_EXIT $? $(date -Is)" >> "$ST"
done

for s in 1 2 3; do
  echo "RIGID_NOFA_SEED_${s}_START $(date -Is)" >> "$ST"
  python -m ffn_sim.scripts.h7_gate_b --device gpu --no-fa \
    --warmup 2000 --sample 1500 --seed "$s" \
    > "$OUT/gate_b_rigid_nofa_seed${s}.log" 2>&1
  echo "RIGID_NOFA_SEED_${s}_EXIT $? $(date -Is)" >> "$ST"
done

echo "ALL_DONE $(date -Is)" >> "$ST"
