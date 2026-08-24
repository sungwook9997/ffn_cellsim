#!/usr/bin/env bash
# PRODUCTION R/N spreading sweep (subdiv-2, all 2026-06-14 fixes) — fit A/A0=a+b/R+c/R^2.
# Per N: aggregate a stable cohesive spheroid, then spread it with the mechanistic
# substrate WETTING (CODE-1). All fixes ON:
#   K_vol=7.73e5 (osmotic incompressibility) + V0=mesh-volume + eversion GUARD (turgor
#   signed-volume clamp) + substrate-wetting (in-plane area-gradient) + stiff floor +
#   tight wetting gate + dt=2e-5 (turgor CFL) + rep=2e9 (non-penetration) + subdiv-2.
# NOTE: no `set -u` (conda cuda-nvcc activate.d trips nounset).
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim
export CUDA_HOME=$CONDA_PREFIX CUDA_PATH=$CONDA_PREFIX
cd ~/ffn_cellsim-platform
ST=/tmp/prod_rn_status.txt; : > "$ST"
for N in 12 30 60 100 200; do
  SPH=/tmp/sph_prod_n${N}.npy
  # 1) AGGREGATE (cohesive, stable subdiv-2)
  python -u -m aleph.scripts.dcm_two_stage_production \
    --n "$N" --agg-only --node-face-contact --subdiv 2 \
    --k-vol 7.73e5 --rep-strength 2e9 --adh-strength 5e7 --f-active 0 \
    --agg-dt 2e-5 --agg-spacing 2.3 --agg-max-steps 50000 --frames 8 \
    --save-spheroid "$SPH" --seed 7 > /tmp/aggprod_n${N}.log 2>&1
  echo "AGG  N=$N done $(date +%H:%M) $(grep -oE 'V/V0=[0-9.]+' /tmp/aggprod_n${N}.log | tail -1)" >> "$ST"
  # 2) SPREAD (mechanistic substrate wetting = CODE-1)
  python -u -m aleph.scripts.dcm_two_stage_production \
    --n "$N" --subdiv 2 --spread-from-npy "$SPH" --node-face-contact \
    --substrate-wetting --k-vol 7.73e5 --rep-strength 2e9 --adh-strength 1e7 \
    --w-cs-jm2 2.85e-3 --spread-dt 2e-5 --spread-steps 120000 --frames 40 \
    --settle-force 4e-10 --seed 7 > /tmp/spreadprod_n${N}.log 2>&1
  AA=$(grep -oE "A/A0=[0-9.]+" /tmp/spreadprod_n${N}.log | tail -1)
  VV=$(grep -oE "V/V0=[0-9.]+" /tmp/spreadprod_n${N}.log | tail -1)
  echo "SPREAD N=$N final ${AA} ${VV} $(date +%H:%M)" >> "$ST"
done
echo "SWEEP DONE $(date +%H:%M)" >> "$ST"
