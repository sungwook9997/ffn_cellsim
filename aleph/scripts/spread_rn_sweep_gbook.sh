#!/usr/bin/env bash
# R/N spread sweep — test the PI law A/A0 = a + b/R + c/R² across spheroid sizes.
# For each N: aggregate a cohesive spheroid (node-face, lit 4e7/5e7), then spread it
# with the mechanistic lamellipodium + pressure-triggered de-cohesion + wetting.
# Small spheroids (few layers) should unfold to a monolayer → high A/A0; large ones
# approach the asymptote a. Runs on gbook GPU, sequential (one GPU).
# NOTE: no `set -u` — conda's cuda-nvcc activate.d references unbound
# NVCC_PREPEND_FLAGS and aborts under nounset (known gbook gotcha).
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim
export CUDA_HOME=$CONDA_PREFIX CUDA_PATH=$CONDA_PREFIX
cd ~/ffn_cellsim-platform
ST=/tmp/spread_rn_status.txt
: > "$ST"
for N in 12 30 60 120 250; do
  SPH=/tmp/sph_rn_n${N}.npy
  # 1) aggregate (cohesive) — small N coheres fast
  python -u -m aleph.scripts.dcm_two_stage_production \
    --n "$N" --agg-only --node-face-contact \
    --rep-strength 4e7 --adh-strength 5e7 --f-active 0 \
    --agg-dt 1e-5 --agg-spacing 2.3 --agg-max-steps 300000 --frames 15 \
    --save-spheroid "$SPH" --seed 7 > /tmp/aggrn_n${N}.log 2>&1
  echo "AGG N=$N done $(date +%H:%M)" >> "$ST"
  # 2) spread (lamellipodium + de-cohesion + wetting); weak spread cohesion 1e7
  python -u -m aleph.scripts.dcm_two_stage_production \
    --n "$N" --spread-from-npy "$SPH" --node-face-contact \
    --lamellipodium --junction-switch \
    --rep-strength 4e7 --adh-strength 1e7 --w-cs-jm2 2.85e-3 \
    --spread-dt 1e-4 --spread-steps 60000 --frames 30 \
    --lamel-pool-per-cell 200 --lamel-vfront-umin 12 --lamel-S 10 \
    --lamel-contact-band 3.0 --settle-force 4e-10 --seed 7 \
    > /tmp/spreadrn_n${N}.log 2>&1
  AA=$(grep -oE "A/A0=[0-9.]+" /tmp/spreadrn_n${N}.log | tail -1)
  echo "SPREAD N=$N final ${AA} $(date +%H:%M)" >> "$ST"
done
echo "SWEEP DONE $(date +%H:%M)" >> "$ST"
