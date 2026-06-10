#!/usr/bin/env bash
set -u
cd "$HOME/ffn_cellsim"
export PYTHONPATH=.
PY="$HOME/miniconda3/envs/ffn_sim/bin/python"
OUT=ffn_sim/outputs/gpu_signoff
echo "[tier1 GPU start $(date +%H:%M:%S)]"
"$PY" ffn_sim/scripts/h3_lp_gpu_production.py --scale medium --device gpu > "$OUT/tier1_h3lp_gpu.log" 2>&1
echo "[tier1 GPU done $(date +%H:%M:%S)] rc=$?"
echo "[tier1 CPU start $(date +%H:%M:%S)]"
"$PY" ffn_sim/scripts/h3_lp_gpu_production.py --scale medium --device cpu --allow-cpu-dev > "$OUT/tier1_h3lp_cpu.log" 2>&1
echo "[tier1 CPU done $(date +%H:%M:%S)] rc=$?"
echo "[tier1 ALL DONE $(date +%H:%M:%S)]"
