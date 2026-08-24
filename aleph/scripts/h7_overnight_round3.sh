#!/usr/bin/env bash
# H.7 overnight production ROUND 3 (gbook GPU). Round-2 found the constrained
# FA-adhered GATE-B build is seed-unstable (4/5 seeds blew the BAOAB guard at
# softstart=warmup//4). Round 3:
#   - GATE-B 5-seed ensemble with a LONG softstart (=warmup) to stabilize the
#     constrained settle -> a real gamma_rigid mean±sd if it holds.
#   - GATE-B ×40 coarse-graining convergence: n_filaments {500, 2000} (1000 is
#     the ensemble) -> is gamma_rigid grid-converged (novelty-defending)?
#
# Usage (gbook):  setsid bash aleph/scripts/h7_overnight_round3.sh &
set -u
cd "$(dirname "$0")/../.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null \
  || source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
conda activate ffn_sim

OUT=aleph/outputs/h7/production/round3
mkdir -p "$OUT"
ST="$OUT/round3.status"
echo "START $(date -Is)" > "$ST"

for s in 1 2 3 4 5; do
  echo "ENS_SEED_${s}_START $(date -Is)" >> "$ST"
  python -m aleph.scripts.h7_gate_b --device gpu --warmup 2000 --softstart 2000 --sample 1500 --seed "$s" \
    > "$OUT/gate_b_ss_seed${s}.log" 2>&1
  echo "ENS_SEED_${s}_EXIT $? $(date -Is)" >> "$ST"
done

for nf in 500 2000; do
  echo "CONV_NF_${nf}_START $(date -Is)" >> "$ST"
  python -m aleph.scripts.h7_gate_b --device gpu --n-filaments "$nf" \
    --warmup 2000 --softstart 2000 --sample 1500 --seed 1 \
    > "$OUT/gate_b_conv_nf${nf}.log" 2>&1
  echo "CONV_NF_${nf}_EXIT $? $(date -Is)" >> "$ST"
done

echo "ALL_DONE $(date -Is)" >> "$ST"
