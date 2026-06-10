#!/usr/bin/env bash
# Reduced-scale H.2 single-filament L_p on GPU vs CPU (same config seed) for the
# P2c sign-off. Full-production H.2 on GPU is impractical at N=21 (GPU per-step
# overhead dominates tiny-N; ~660 steps/s vs CPU) — this reduced run validates
# L_p parity, not speed. Sequential (shared h2_trajectory.npz output → copy each).
set -u
cd "$HOME/ffn_cellsim"
export PYTHONPATH=.
PY="$HOME/miniconda3/envs/ffn_sim/bin/python"
OUT=ffn_sim/outputs/gpu_signoff
ARGS="--n-equilibrate 50000 --n-sample 200000 --sample-interval 1000"
echo "[h2v2 GPU start $(date +%H:%M:%S)]"
"$PY" ffn_sim/scripts/h2_single_filament.py --device gpu $ARGS > "$OUT/h2v2_gpu.log" 2>&1
cp ffn_sim/outputs/h2/h2_trajectory.npz "$OUT/h2_traj_gpu.npz"
echo "[h2v2 GPU done $(date +%H:%M:%S)]"
"$PY" ffn_sim/scripts/h2_single_filament.py --device cpu --allow-cpu-dev $ARGS > "$OUT/h2v2_cpu.log" 2>&1
cp ffn_sim/outputs/h2/h2_trajectory.npz "$OUT/h2_traj_cpu.npz"
echo "[h2v2 ALL DONE $(date +%H:%M:%S)]"
