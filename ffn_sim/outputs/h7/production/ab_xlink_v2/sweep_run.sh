#!/usr/bin/env bash
# ab_xlink_v2 transmission redesign (Slater 2021) — reduced CPU/GPU preliminary.
# Does active γ_soft respond to crosslink DENSITY (n_xl) + Bell off-rate (k_off0)?
# (v1 tested STIFFNESS → γ ~5%; Slater says density+kinetics are the transmission knobs.)
# Fixed reduced scale; varies one knob at a time vs the phase1_h3 anchors
# (n_xl=1000, alpha_k_off0=0.066 Ferrer2008). Sensitivity-only, no production change.
set -u
cd "$HOME/ffn_cellsim"; export PYTHONPATH=.
PY="$HOME/miniconda3/envs/ffn_sim/bin/python"
OUT=ffn_sim/outputs/h7/production/ab_xlink_v2
COMMON="--n-filaments 60 --warmup 800 --contract-steps 8000 --sample-every 4000 --seed 1 --device gpu"
run () {  # tag  xl_n  koff0
  echo "[$(date +%H:%M:%S)] $1  n_xl=$2 koff0=$3"
  "$PY" ffn_sim/scripts/h7_active_force_budget.py $COMMON \
    --xl-n "$2" --xl-koff0 "$3" \
    --out-json "$OUT/$1.json" --out-png "$OUT/$1.png" > "$OUT/$1.log" 2>&1
}
# density axis (16× range) at anchor k_off0
run xln250_k066   250  0.066
run xln1000_k066  1000 0.066
run xln4000_k066  4000 0.066
# Bell off-rate axis (10× slower = more connected) at anchor density
run xln1000_k0066 1000 0.0066
echo "[$(date +%H:%M:%S)] ALL DONE"
