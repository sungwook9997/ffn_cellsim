#!/bin/bash
# Master aggregation sweep for gbook (run under setsid; survives ssh/wifi drops).
# node-FACE cohesion (genuine surface contact) + pure cohesion (f_active=0), per N.
# Writes /tmp/spheroid_nf_n{N}.npy + /tmp/agg_nf_n{N}.log + appends /tmp/agg_master_status.txt.
set +u                                  # conda activate.d (cuda-nvcc) refs unset vars
export NVCC_PREPEND_FLAGS="${NVCC_PREPEND_FLAGS:-}"
export NVCC_APPEND_FLAGS="${NVCC_APPEND_FLAGS:-}"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim
# detached conda leaves CUDA_HOME empty (cuda-nvcc activate.d aborts on the unbound
# NVCC_PREPEND_FLAGS) → the node-face RawKernel NVRTC compile hard-crashes. Set it.
export CUDA_HOME="$CONDA_PREFIX"
export CUDA_PATH="$CONDA_PREFIX"
cd ~/ffn_cellsim-platform
ST=/tmp/agg_master_status.txt
echo "SWEEP START $(date)" > "$ST"
for N in 400 600 800 1000; do
  # more steps for bigger N (convergence early-stops at the plateau regardless)
  STEPS=$(( 30000 + (N - 400) * 50 ))
  echo "START N=$N steps=$STEPS $(date)" >> "$ST"
  python -m ffn_sim.scripts.dcm_two_stage_production \
    --n "$N" --agg-only --node-face-contact \
    --rep-strength 4e7 --adh-strength 5e7 --f-active 0 \
    --agg-dt 1e-5 --agg-spacing 2.3 --agg-max-steps "$STEPS" --frames 30 \
    --save-spheroid /tmp/spheroid_nf_n${N}.npy --seed 7 \
    > /tmp/agg_nf_n${N}.log 2>&1
  RC=$?
  echo "DONE N=$N rc=$RC $(date) $(grep -c AGG-ONLY /tmp/agg_nf_n${N}.log)" >> "$ST"
done
echo "SWEEP ALL DONE $(date)" >> "$ST"
