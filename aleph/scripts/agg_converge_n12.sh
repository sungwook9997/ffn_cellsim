#!/usr/bin/env bash
# PROPER aggregation to TRUE convergence (PI: the prior aggregate was unconverged —
# contact still rising 0.52->0.91, asph WORSENING 0.22->0.42, f_active=0 gave an
# irregular collapsing clump not a round spheroid). Now: motility ON (f_active>0 so
# cells round up), regime (a) rep=2e8/dt=8e-6, long cap + built-in convergence
# detection (contact plateau + Rg plateau + asph<0.10). Save the converged spheroid.
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim
export CUDA_HOME=$CONDA_PREFIX CUDA_PATH=$CONDA_PREFIX
cd ~/ffn_cellsim-platform
echo "=== AGG-CONVERGE N=12 subdiv-2 rep=2e8 dt=8e-6 f_active=1.6e-10 (motility ON) ==="
python -u -m aleph.scripts.dcm_two_stage_production \
  --n 12 --agg-only --node-face-contact --subdiv 2 \
  --k-vol 7.73e5 --rep-strength 2e8 --adh-strength 5e7 --f-active 1.6e-10 \
  --agg-dt 8e-6 --agg-spacing 2.3 --agg-max-steps 160000 --frames 40 \
  --save-spheroid /tmp/agg_converged_n12.npy --seed 7 2>&1
echo "=== AGG-CONVERGE DONE ==="
