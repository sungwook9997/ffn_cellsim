#!/usr/bin/env bash
# N=12 subdiv-2 MECHANISTIC-STACK spread diagnostic: spread the SAME aggregate that
# compacted under wetting+proxy, but with the full-fidelity stack (mechanistic
# lamellipodium REPLACING the contracting body-force proxy + junction-switch
# de-cohesion + substrate wetting). Same stability regime (a) rep=2e8 dt=8e-6 +
# K_vol fix + eversion guard. Answers: does the mechanistic stack spread the
# spheroid now that K_vol is fixed (memory's compaction was at the old soft K_vol)?
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim
export CUDA_HOME=$CONDA_PREFIX CUDA_PATH=$CONDA_PREFIX
cd ~/ffn_cellsim-platform
SPH=/tmp/verify_sph_n12.npy   # the SAME aggregate as the wetting+proxy run
echo "=== MECH-STACK spread (lamellipodium + junction-switch + wetting) rep=2e8 dt=8e-6 ==="
python -u -m aleph.scripts.dcm_two_stage_production \
  --n 12 --subdiv 2 --spread-from-npy "$SPH" --node-face-contact \
  --lamellipodium --junction-switch --substrate-wetting \
  --k-vol 7.73e5 --rep-strength 2e8 --adh-strength 1e7 --w-cs-jm2 2.85e-3 \
  --spread-dt 8e-6 --spread-steps 40000 --frames 20 \
  --lamel-pool-per-cell 200 --lamel-vfront-umin 12 --lamel-S 10 \
  --lamel-contact-band 3.0 --settle-force 4e-10 --seed 7 2>&1
echo "=== MECH-STACK DONE ==="
cp aleph/outputs/h_dcm_two_stage/two_stage_n12.pkl /tmp/mech_stack_n12.pkl
echo "saved /tmp/mech_stack_n12.pkl"
