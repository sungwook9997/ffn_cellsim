#!/usr/bin/env bash
# M-band targeting A/B (loop24d verification) — CPU preliminary.
# For each seed, runs the same-seed paired coherent-traction differential with
# M-band targeting ON vs OFF (random greedy placement). Same geometry as the
# reference h7_sf_2c_sarcomeric_gpu run (n_fil=36/n_beads=144/n_motors=16).
# Per seed the two modes run in parallel (2 procs); seeds run sequentially.
set -u
cd "$HOME/ffn_cellsim"
export PYTHONPATH=.
PY="$HOME/miniconda3/envs/ffn_sim/bin/python"

OUT=ffn_sim/outputs/h7/production/mband_ab
SCRIPT=ffn_sim/scripts/h7_ventral_sf_traction.py
COMMON="--sarcomeric --n-filaments 6 --n-sarcomeres 3 --n-motors 16 \
        --equilibrate 120000 --contract 120000 --n-samples 20"

for SEED in 1 2 3; do
  echo "[$(date +%H:%M:%S)] seed=$SEED  launching on+off ..."
  "$PY" $SCRIPT $COMMON --seed "$SEED" --mband-mode on  \
      --out-json "$OUT/ab_s${SEED}_on.json"  > "$OUT/ab_s${SEED}_on.log"  2>&1 &
  PID_ON=$!
  "$PY" $SCRIPT $COMMON --seed "$SEED" --mband-mode off \
      --out-json "$OUT/ab_s${SEED}_off.json" > "$OUT/ab_s${SEED}_off.log" 2>&1 &
  PID_OFF=$!
  wait $PID_ON $PID_OFF
  echo "[$(date +%H:%M:%S)] seed=$SEED  done"
done
echo "[$(date +%H:%M:%S)] ALL DONE"
