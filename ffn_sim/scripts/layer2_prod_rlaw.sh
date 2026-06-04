#!/usr/bin/env bash
# Layer-2 B2 PRODUCTION — the R0-law (Bare, catch cohesion, no substrate/traction) swept across
# the full PI R0 range on the gbook RTX A5000. PI-ratified DOE (2026-06-04):
#   R0 100-400 um x7  (N0 = 1250 4000 8500 16000 26000 40000 60000)
#   5 seeds/R0  x  60 h (2 doublings)  =  35 runs, sequential on the one GPU.
# Each run stamps git commit + cfg hash + actual device + the raw-area footprint A/A0.
# RESUMABLE: a (N0,seed) already present in the output JSONL is skipped — safe to re-run after
# an interruption. Run from the repo root on gbook:
#   PYTHONPATH=~/ffn_cellsim bash ffn_sim/scripts/layer2_prod_rlaw.sh
set -uo pipefail

PY="${PY:-$HOME/miniconda3/envs/ffn_sim/bin/python}"
ROOT="${ROOT:-$HOME/ffn_cellsim}"
OUT="$ROOT/ffn_sim/outputs/layer2/prod_rlaw"
JSONL="$OUT/rlaw_sweep.jsonl"
PROG="$OUT/rlaw_sweep.progress"
mkdir -p "$OUT"; touch "$JSONL"

N0S=(1250 4000 8500 16000 26000 40000 60000)
SEEDS=(0 1 2 3 4)

export PYTHONPATH="$ROOT"
echo "[$(date +%F\ %T)] PRODUCTION R0-law start: ${#N0S[@]} sizes x ${#SEEDS[@]} seeds" >> "$PROG"
for N in "${N0S[@]}"; do
  for S in "${SEEDS[@]}"; do
    # resume: skip if a line with this exact "n" and "seed" already landed
    if grep -q "\"n\": $N,[^}]*\"seed\": $S," "$JSONL" 2>/dev/null; then
      echo "[$(date +%T)] skip N0=$N seed=$S (already done)" >> "$PROG"; continue
    fi
    echo "[$(date +%T)] run N0=$N seed=$S" >> "$PROG"
    "$PY" -m ffn_sim.scripts.layer2_gpu_scaleup run "$N" "$S" gpu >> "$JSONL" 2>> "$OUT/rlaw_sweep.err"
  done
done
echo "[$(date +%F\ %T)] PRODUCTION R0-law DONE ($(wc -l < "$JSONL") lines)" >> "$PROG"
# auto-fit at the end (production-driver convention)
"$PY" -m ffn_sim.scripts.layer2_gpu_scaleup --fit "$JSONL" >> "$OUT/rlaw_fit.log" 2>&1 || true
