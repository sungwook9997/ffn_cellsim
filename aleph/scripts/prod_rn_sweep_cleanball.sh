#!/usr/bin/env bash
# PRODUCTION R/N spreading sweep — CLEAN FCC ball + WETTING ONLY (proxy OFF).
# The 2026-06-14 breakthrough: the spheroid DOES spread once you (1) start from a
# clean round ball (the driver's own aggregation over-compacts into an irregular
# clump — PI hyp 1) and (2) turn OFF the contracting active rim-traction proxy
# (--spread-f-act 0). N=12 then spreads A/A0 1.0->1.87 peak (vs 0.74 from the broken
# aggregate). This sweep fits A/A0 = a + b/R + c/R^2 across spheroid sizes.
# All foundation fixes ON: K_vol=7.73e5, eversion guard, regime (a) rep=2e8/dt=8e-6,
# substrate-wetting at lit MCF7 W_cs=2.85e-3, node-face contact, subdiv-2.
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim
export CUDA_HOME=$CONDA_PREFIX CUDA_PATH=$CONDA_PREFIX
cd ~/ffn_cellsim-platform
OUT=aleph/outputs/h_dcm_two_stage
ST=/tmp/cleanball_sweep_status.txt; : > "$ST"
for N in 12 30 60 100 200; do
  BALL=/tmp/cleanball_n${N}.npy
  python -m aleph.scripts.make_clean_ball --n "$N" --subdiv 2 \
    --spacing-factor 2.1 --mode 3d --out "$BALL" >> "$ST" 2>&1
  python -u -m aleph.scripts.dcm_two_stage_production \
    --n "$N" --subdiv 2 --spread-from-npy "$BALL" --node-face-contact \
    --substrate-wetting --spread-f-act 0 \
    --k-vol 7.73e5 --rep-strength 2e8 --adh-strength 1e7 --w-cs-jm2 2.85e-3 \
    --spread-dt 8e-6 --spread-steps 60000 --frames 31 \
    --settle-force 4e-10 --seed 7 > /tmp/cleanball_spread_n${N}.log 2>&1
  cp "$OUT/two_stage_n${N}.pkl" "$OUT/cleanball_n${N}.pkl" 2>/dev/null
  AA=$(grep -oE "A/A0=[0-9.]+" /tmp/cleanball_spread_n${N}.log | tail -1)
  PK=$(grep -oE "A/A0=[0-9.]+" /tmp/cleanball_spread_n${N}.log | grep -oE "[0-9.]+" | sort -rn | head -1)
  echo "N=$N final ${AA} peak=${PK} $(date +%H:%M)" >> "$ST"
done
echo "SWEEP DONE $(date +%H:%M)" >> "$ST"
