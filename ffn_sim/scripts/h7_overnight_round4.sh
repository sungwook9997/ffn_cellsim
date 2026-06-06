#!/usr/bin/env bash
# H.7 overnight production ROUND 4 (gbook GPU). Round-3 confirmed softstart length
# is NOT the lever; the real failure is the FA integrin catch-bond force spiking
# past the Pereverzev pN range under the constrained path's large dt
# (|F|/F_s = 954 > 700). Round 4:
#   - constrained GATE-B 5-seed ensemble at dt-safety 0.3 (smaller constrained
#     step -> smaller integrin force spikes): does it stabilize -> gamma_rigid mean±sd?
#   - UNCONSTRAINED soft-backbone GATE-B 5-seed ensemble (inherently stable; the
#     M-SHAKE-shunt cross-check): gamma_soft distribution vs the constrained gamma_rigid.
#
# Usage (gbook):  setsid bash ffn_sim/scripts/h7_overnight_round4.sh &
set -u
cd "$(dirname "$0")/../.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
conda activate ffn_sim

OUT=ffn_sim/outputs/h7/production/round4
mkdir -p "$OUT"
ST="$OUT/round4.status"
echo "START $(date -Is)" > "$ST"

for s in 1 2 3 4 5; do
  echo "DTSAFE_SEED_${s}_START $(date -Is)" >> "$ST"
  python -m ffn_sim.scripts.h7_gate_b --device gpu --warmup 2000 --sample 1500 \
    --dt-safety 0.3 --seed "$s" \
    > "$OUT/gate_b_dt03_seed${s}.log" 2>&1
  echo "DTSAFE_SEED_${s}_EXIT $? $(date -Is)" >> "$ST"
done

for s in 1 2 3 4 5; do
  echo "SOFT_SEED_${s}_START $(date -Is)" >> "$ST"
  python -m ffn_sim.scripts.h7_gate_b --device gpu --warmup 2000 --sample 1500 \
    --no-constrained --seed "$s" \
    > "$OUT/gate_b_soft_seed${s}.log" 2>&1
  echo "SOFT_SEED_${s}_EXIT $? $(date -Is)" >> "$ST"
done

echo "ALL_DONE $(date -Is)" >> "$ST"
