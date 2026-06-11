#!/usr/bin/env bash
# TWO-STAGE DCM spheroid production sweep — runs on the gbook A5000 in a screen.
# Sequential cell-count sweep (aggregation → spreading), scaling toward production N.
# Each N writes outputs/h_dcm_two_stage/two_stage_n{N}.pkl (the MacBook pulls + renders).
set -u
cd ~/ffn_cellsim-platform || exit 1
PY=~/miniconda3/envs/ffn_sim/bin/python
OUTDIR=ffn_sim/outputs/h_dcm_two_stage
mkdir -p "$OUTDIR"
SWEEP="${SWEEP_N:-100 300 600 1000 1500}"
AGG="${AGG_STEPS:-25000}"
SPREAD="${SPREAD_STEPS:-25000}"
FRAMES="${FRAMES:-14}"
echo "[$(date)] SWEEP START: N={$SWEEP} agg=$AGG spread=$SPREAD" >> "$OUTDIR/SWEEP.log"
for N in $SWEEP; do
  echo "[$(date)] N=$N START" >> "$OUTDIR/SWEEP.log"
  PYTHONPATH=. "$PY" ffn_sim/scripts/dcm_two_stage_production.py \
      --n "$N" --agg-steps "$AGG" --spread-steps "$SPREAD" --frames "$FRAMES" \
      >> "$OUTDIR/run_n$N.log" 2>&1
  rc=$?
  echo "[$(date)] N=$N DONE rc=$rc" >> "$OUTDIR/SWEEP.log"
done
echo "[$(date)] SWEEP COMPLETE" >> "$OUTDIR/SWEEP.log"
