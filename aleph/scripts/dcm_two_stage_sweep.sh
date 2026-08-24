#!/usr/bin/env bash
# TWO-STAGE DCM spheroid production sweep v2 (BIOLOGICAL aggregation) — gbook A5000.
# Stage 1 = emergent active-matter coalescence (SPP motility + tent adhesion + soft
# drop), seeded JUST OUT of adhesion range (2.8·R) so cells must motility-coalesce
# (negative control: motility OFF => no aggregation). Stage 2 = spreading from the
# converged aggregate. Each N writes outputs/h_dcm_two_stage/two_stage_n{N}.pkl.
set -u
cd ~/ffn_cellsim-platform || exit 1
PY=~/miniconda3/envs/ffn_sim/bin/python
OUTDIR=aleph/outputs/h_dcm_two_stage
mkdir -p "$OUTDIR"
SWEEP="${SWEEP_N:-100 200 400}"
SPACING="${AGG_SPACING:-2.8}"
AGG_DT="${AGG_DT:-1e-8}"
AGG_MAX="${AGG_MAX_STEPS:-400000}"
SPREAD="${SPREAD_STEPS:-30000}"
FRAMES="${FRAMES:-16}"
FACT="${F_ACTIVE:-1.6e-10}"
echo "[$(date)] SWEEP v2 START: N={$SWEEP} spacing=$SPACING agg_dt=$AGG_DT agg_max=$AGG_MAX f_active=$FACT" >> "$OUTDIR/SWEEP.log"
for N in $SWEEP; do
  echo "[$(date)] N=$N START" >> "$OUTDIR/SWEEP.log"
  PYTHONPATH=. "$PY" aleph/scripts/dcm_two_stage_production.py \
      --n "$N" --agg-spacing "$SPACING" --agg-dt "$AGG_DT" \
      --agg-max-steps "$AGG_MAX" --spread-steps "$SPREAD" --frames "$FRAMES" \
      --f-active "$FACT" \
      >> "$OUTDIR/run_n$N.log" 2>&1
  echo "[$(date)] N=$N DONE rc=$?" >> "$OUTDIR/SWEEP.log"
done
echo "[$(date)] SWEEP COMPLETE" >> "$OUTDIR/SWEEP.log"
