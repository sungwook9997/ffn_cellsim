#!/bin/bash
# PI fallback (scheduled ~2 PM): today's node-face two-stage with spread-dt=2e-3.
# Runs aggregate (node-face genuine cohesion, f_active=0) -> spread (dt=2e-3) for N=400.
# Self-contained on gbook (setsid sleeper launches it); survives ssh/wifi.
set +u
export NVCC_PREPEND_FLAGS="${NVCC_PREPEND_FLAGS:-}"
export NVCC_APPEND_FLAGS="${NVCC_APPEND_FLAGS:-}"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim
export CUDA_HOME="$CONDA_PREFIX"   # node-face RawKernel NVRTC needs this (detached)
export CUDA_PATH="$CONDA_PREFIX"
cd ~/ffn_cellsim-platform
ST=/tmp/fallback_2pm_status.txt
echo "FALLBACK START $(date)" > "$ST"
python -m aleph.scripts.dcm_two_stage_production \
  --n 400 --node-face-contact --rep-strength 4e7 --adh-strength 5e7 \
  --f-active 0 --agg-dt 1e-5 --agg-spacing 2.3 --agg-max-steps 30000 \
  --spread-steps 15000 --spread-dt 2e-3 --frames 30 --seed 7 \
  > /tmp/fallback_2pm_run.log 2>&1
echo "FALLBACK DONE rc=$? $(date)" >> "$ST"
