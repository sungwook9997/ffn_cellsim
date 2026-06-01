#!/usr/bin/env bash
# ffn_cellsim — turnkey bootstrap for a fresh CUDA GPU VM (cloud or HPC).
#
# Brings a bare Ubuntu+NVIDIA-driver box (RunPod / Lambda / Vast / Paperspace /
# KAIST-KISTI) to "ready to run production" in ~10-15 min: installs miniconda
# (if absent), creates the `ffn_sim` conda env with the GPU HOOMD build + cupy,
# and runs a GPU smoke so you KNOW it works before launching a multi-hour run.
#
# Prereqs on the box: an NVIDIA driver (`nvidia-smi` works) and the repo present
# at $REPO (clone it or rsync it first — Syncthing does not reach cloud):
#   git clone <remote> ~/ffn_cellsim        # or: rsync -a ./ffn_sim vm:~/ffn_cellsim/ffn_sim
#
# Usage:
#   bash ffn_sim/scripts/cloud_bootstrap.sh
#   # then, e.g.:
#   conda activate ffn_sim
#   PYTHONPATH=$REPO python -m ffn_sim.scripts.h3_ku35_gripwalk_tier1 \
#       --device gpu --arm grip_walk --n-fil 3000 --force-scaling \
#       --dt-factor 2e-4 --v0-accel 100 --n-warmup 12000 --n-sample 20 \
#       --interval 60000 --gsd-period 60000 --out out/native.json
set -euo pipefail

REPO="${REPO:-$HOME/ffn_cellsim}"
ENV_NAME="ffn_sim"
CUPY_PKG="${CUPY_PKG:-cupy-cuda12x}"   # match host driver; 12x is fwd-compat to CUDA 13

echo "== ffn_cellsim cloud bootstrap =="
echo "REPO=$REPO  ENV=$ENV_NAME  CUPY=$CUPY_PKG"

command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader \
  || { echo "!! no nvidia-smi — this box has no GPU/driver. Aborting."; exit 1; }

[ -d "$REPO/ffn_sim" ] || { echo "!! repo not found at $REPO (clone/rsync it first)."; exit 1; }

# --- miniconda ---
if ! command -v conda >/dev/null 2>&1; then
  echo "== installing miniconda =="
  curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/mc.sh
  bash /tmp/mc.sh -b -p "$HOME/miniconda3"
fi
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# --- env (conda-forge auto-selects the GPU hoomd build via __cuda) ---
if ! conda env list | grep -q "^${ENV_NAME} "; then
  echo "== creating env $ENV_NAME from environment-gpu.yml =="
  conda env create -f "$REPO/environment-gpu.yml" \
    || CONDA_OVERRIDE_CUDA=12.4 conda env create -f "$REPO/environment-gpu.yml"
fi
conda activate "$ENV_NAME"

# --- cupy (pip, matched to host CUDA driver) ---
python -c "import cupy" 2>/dev/null || { echo "== pip install $CUPY_PKG =="; pip install "$CUPY_PKG"; }

# --- GPU smoke: assert the GPU build + cupy + a constrained-Action match ---
echo "== GPU smoke =="
python - <<'PY'
import hoomd, cupy
assert hoomd.version.gpu_enabled, "HOOMD is the CPU-only build! re-create env with CONDA_OVERRIDE_CUDA=12.4"
print("HOOMD", hoomd.version.version, "gpu_enabled", hoomd.version.gpu_enabled,
      "| cupy", cupy.__version__, "| CUDA devs", cupy.cuda.runtime.getDeviceCount())
PY
echo "== per-function GPU-vs-CPU match test =="
( cd "$REPO" && PYTHONPATH="$REPO" python -m pytest ffn_sim/tests/test_constrained_baoab_gpu.py -q ) \
  || echo "!! GPU match test failed — investigate before a production run."

echo "== READY. Activate with: conda activate $ENV_NAME ; export PYTHONPATH=$REPO =="
echo "   Launch long runs under tmux/nohup so they survive disconnects."
