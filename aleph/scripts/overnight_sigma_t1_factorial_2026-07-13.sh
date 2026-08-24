#!/usr/bin/env bash
# ---------------------------------------------------------------------------------------------------
# sigma  X  T1-rate  2x2 FACTORIAL  (native N=400 full-compartment, PI-directed 2026-07-13)
#
# Question: the aggregate-sigma Foty-Steinberg driver gives END-STATE mechanical compaction that
# SATURATES ~100s; the rate/event-driven T1-KMC supplies the slow tissue rearrangement that large-dt
# mechanics correctly FREEZE (timescale gap = FUNDAMENTAL, T1-event-limited). Does layering T1 on top
# of sigma let the aggregate keep packing PAST the mechanical plateau toward biology-time?
#
# Design: ALL FOUR conditions share integrator=BDF2 + accel_dt=8e-3 + the full physiological stack
# (nucleus E_nuc=399, cortex gamma, cadherin bundle-10, IPC-Newton). The ONLY differences are the two
# independent levers:  SIGMA in {0, 5e-3 N/m}  x  T1_RATE in {0, 1}.  This decomposes them cleanly:
#   c00_base : sigma OFF, T1 OFF  -> loose aggregate, large-dt mechanics only (negative control)
#   c10_sig  : sigma ON,  T1 OFF  -> mechanical compaction (expected: drops then plateaus ~100s)
#   c01_t1   : sigma OFF, T1 ON   -> fluidisation/rearrangement WITHOUT mechanical drive
#   c11_both : sigma ON,  T1 ON   -> the combined hypothesis
#
# Single A5000 (cuda:0) => the 4 runs are SEQUENTIAL. Monitor by ~/ff_scratch/_prod_out/factorial_master.log.
# ---------------------------------------------------------------------------------------------------
set -u
cd ~/ff_scratch || exit 1
export PYTHONPATH=~/ff_scratch
PY=~/miniconda3/envs/ffn_sim/bin/python
OUT=~/ff_scratch/_prod_out
mkdir -p "$OUT"
MASTER=$OUT/factorial_master.log

STEPS=${STEPS:-37500}          # biology-time window; 37500 x 8e-3s = 300s physical (override via env)
PHYS=$(awk "BEGIN{printf \"%.0f\", $STEPS*0.008}")
COMMON="NCELLS=400 WARMUP=500 ACCEL_DT=8e-3 GAP=2.4 BUNDLE=10 ENUC=399 RNUC=0.7 GAMMA=5e-4 NUCLEUS=1 FRAMES=40 BDF2=1 T1_MODE=rigid T1_CADENCE=200 T1_SEED=13 DEVICE=cuda:0 STEPS=$STEPS"

echo "[FACTORIAL START] $(date)  steps=$STEPS (phys ${PHYS}s)  integrator=BDF2 accel_dt=8e-3 N=400" | tee -a "$MASTER"

run() {  # $1=tag  $2=SIGMA  $3=T1_RATE
  local tag=$1 sig=$2 t1=$3
  echo "[RUN $tag] SIGMA=$sig T1_RATE=$t1 START $(date)" | tee -a "$MASTER"
  env $COMMON SIGMA="$sig" T1_RATE="$t1" TAG="_$tag" \
    "$PY" aleph/scripts/_gbook_aggregate_compaction.py > "$OUT/fac_$tag.log" 2>&1
  local rc=$?
  grep -E "\[SUMMARY" "$OUT/fac_$tag.log" | tail -1 | sed "s|^|[$tag] |" | tee -a "$MASTER"
  echo "[RUN $tag] END rc=$rc $(date)" | tee -a "$MASTER"
}

run c00_base 0      0
run c10_sig  5.0e-3 0
run c01_t1   0      1
run c11_both 5.0e-3 1

echo "[FACTORIAL DONE] $(date)" | tee -a "$MASTER"
