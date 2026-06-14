#!/usr/bin/env bash
# N=12 subdiv-2 stability verification (regime (a): rep=2e8, dt=8e-6, CONSISTENT
# agg+spread) — confirm V/V0 stays bounded through BOTH aggregation and wetting
# spread, past the rep=2e9 divergence window (~step 8k). Mirrors the two-stage
# production path exactly (agg-only save npy → spread-from-npy). gbook GPU.
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim
export CUDA_HOME=$CONDA_PREFIX CUDA_PATH=$CONDA_PREFIX
cd ~/ffn_cellsim-platform
REP=2e8; DT=8e-6
SPH=/tmp/verify_sph_n12.npy
echo "=== VERIFY (a) rep=$REP dt=$DT subdiv-2 N=12 ==="
echo "--- STAGE 1 AGG (rep=$REP dt=$DT adh=5e7) ---"
python -u -m ffn_sim.scripts.dcm_two_stage_production \
  --n 12 --agg-only --node-face-contact --subdiv 2 \
  --k-vol 7.73e5 --rep-strength $REP --adh-strength 5e7 --f-active 0 \
  --agg-dt $DT --agg-spacing 2.3 --agg-max-steps 25000 --frames 10 \
  --save-spheroid "$SPH" --seed 7 2>&1
echo "--- STAGE 2 SPREAD (rep=$REP dt=$DT adh=1e7 wetting) ---"
python -u -m ffn_sim.scripts.dcm_two_stage_production \
  --n 12 --subdiv 2 --spread-from-npy "$SPH" --node-face-contact \
  --substrate-wetting --k-vol 7.73e5 --rep-strength $REP --adh-strength 1e7 \
  --w-cs-jm2 2.85e-3 --spread-dt $DT --spread-steps 40000 --frames 20 \
  --settle-force 4e-10 --seed 7 2>&1
echo "=== VERIFY DONE ==="
