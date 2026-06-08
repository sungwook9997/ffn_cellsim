#!/usr/bin/env bash
# H.7 Gate-A device-BAOAB launcher for gbook GPU.
#
# Use this file launcher instead of a long inline ssh/bash -c command.  It
# creates the log before Python starts, activates conda, sets the opt-in
# device BAOAB env var, warms CuPy, and then execs the Gate-A driver with
# unbuffered output.
#
# Typical launch on gbook:
#   cd ~/ffn_cellsim
#   nohup bash ffn_sim/scripts/h7_gate_a_device_launcher.sh >/tmp/h7_gate_a_launch.out 2>&1 &
# or:
#   setsid bash ffn_sim/scripts/h7_gate_a_device_launcher.sh &
#
# Tunable environment variables:
#   N_FIL=300 BATCH_STEPS=200000 TICKS=1000 MEASURE_EVERY=10 MEASURE_MODE=final FREEZE_XLINKS=1 SNAPSHOT_FREE_MYOSIN=1 MYOSIN_WARM_BIND_TICKS=1 SEED=1
#   DEVICE=gpu LOG_BASENAME=my_run

set -Eeuo pipefail

cd "$(dirname "$0")/../.." || exit 1

OUT_DIR="${OUT_DIR:-ffn_sim/outputs/h7/production/device_gate_a}"
FIG_DIR="${FIG_DIR:-ffn_sim/outputs/h7/figs}"
mkdir -p "$OUT_DIR" "$FIG_DIR"

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_BASENAME="${LOG_BASENAME:-gate_a_device_${STAMP}}"
LOG="$OUT_DIR/${LOG_BASENAME}.log"
STATUS="$OUT_DIR/${LOG_BASENAME}.status"
DMON="$OUT_DIR/${LOG_BASENAME}.nvidia_dmon.log"

: > "$LOG"
exec >> "$LOG" 2>&1

echo "BOOT $(date -Is)"
echo "PWD=$(pwd)"
echo "HOST=$(hostname)"
echo "PID=$$"

cleanup() {
  code=$?
  if [[ -n "${DMON_PID:-}" ]]; then
    kill "$DMON_PID" 2>/dev/null || true
  fi
  echo "EXIT ${code} $(date -Is)" | tee -a "$STATUS"
}
trap cleanup EXIT

echo "START $(date -Is)" > "$STATUS"

source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null \
  || source "$HOME/anaconda3/etc/profile.d/conda.sh" 2>/dev/null \
  || source "$HOME/miniforge3/etc/profile.d/conda.sh" 2>/dev/null
conda activate ffn_sim

export PYTHONUNBUFFERED=1
export FFN_GPU_DEVICE_BAOAB="${FFN_GPU_DEVICE_BAOAB:-1}"
export CUPY_CACHE_DIR="${CUPY_CACHE_DIR:-$(pwd)/ffn_sim/outputs/cupy_cache}"
mkdir -p "$CUPY_CACHE_DIR"

echo "ENV $(date -Is)"
echo "python=$(command -v python)"
echo "FFN_GPU_DEVICE_BAOAB=$FFN_GPU_DEVICE_BAOAB"
echo "CUPY_CACHE_DIR=$CUPY_CACHE_DIR"
nvidia-smi || true

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi dmon -s pucm -d 1 -o DT > "$DMON" 2>&1 &
  DMON_PID=$!
  echo "nvidia-smi dmon pid=$DMON_PID log=$DMON"
fi

echo "PYTHON_IMPORT_WARMUP_START $(date -Is)"
python - <<'PY'
import hoomd
print("hoomd", hoomd.version.version, flush=True)
try:
    import cupy as cp
    print("cupy", cp.__version__, flush=True)
    x = cp.random.default_rng(1).standard_normal((100000, 3))
    y = x * 2.0 + 1.0
    cp.cuda.Stream.null.synchronize()
    print("cupy_warmup_ok", y.shape, flush=True)
except Exception as exc:
    print("cupy_warmup_failed", repr(exc), flush=True)
    raise
PY
echo "PYTHON_IMPORT_WARMUP_DONE $(date -Is)"

N_FIL="${N_FIL:-300}"
BATCH_STEPS="${BATCH_STEPS:-200000}"
TICKS="${TICKS:-1000}"
MEASURE_EVERY="${MEASURE_EVERY:-10}"
MEASURE_MODE="${MEASURE_MODE:-final}"
FREEZE_XLINKS="${FREEZE_XLINKS:-1}"
SNAPSHOT_FREE_MYOSIN="${SNAPSHOT_FREE_MYOSIN:-1}"
MYOSIN_WARM_BIND_TICKS="${MYOSIN_WARM_BIND_TICKS:-1}"
SNAPSHOT_FREE_GRIP_CLOCK="${SNAPSHOT_FREE_GRIP_CLOCK:-1}"
SEED="${SEED:-1}"
DEVICE="${DEVICE:-gpu}"

JSON_OUT="$OUT_DIR/${LOG_BASENAME}.json"
FIG_OUT="$FIG_DIR/${LOG_BASENAME}.png"

echo "RUN_START $(date -Is)" | tee -a "$STATUS"
echo "args: n_fil=$N_FIL batch_steps=$BATCH_STEPS ticks=$TICKS measure_every=$MEASURE_EVERY measure_mode=$MEASURE_MODE freeze_xlinks=$FREEZE_XLINKS snapshot_free_myosin=$SNAPSHOT_FREE_MYOSIN myosin_warm_bind_ticks=$MYOSIN_WARM_BIND_TICKS snapshot_free_grip_clock=$SNAPSHOT_FREE_GRIP_CLOCK seed=$SEED device=$DEVICE"

FREEZE_ARGS=()
if [[ "$FREEZE_XLINKS" == "1" ]]; then
  FREEZE_ARGS+=(--freeze-xlinks)
fi

MYOSIN_FAST_ARGS=()
if [[ "$SNAPSHOT_FREE_MYOSIN" == "1" ]]; then
  MYOSIN_FAST_ARGS+=(--snapshot-free-myosin --myosin-warm-bind-ticks "$MYOSIN_WARM_BIND_TICKS")
  if [[ "$SNAPSHOT_FREE_GRIP_CLOCK" != "1" ]]; then
    MYOSIN_FAST_ARGS+=(--no-snapshot-free-grip-clock)
  fi
fi

python -u -m ffn_sim.scripts.h7_gate_a \
  --n-fil "$N_FIL" \
  --batch-steps "$BATCH_STEPS" \
  --ticks "$TICKS" \
  --measure-every "$MEASURE_EVERY" \
  --measure-mode "$MEASURE_MODE" \
  "${FREEZE_ARGS[@]}" \
  "${MYOSIN_FAST_ARGS[@]}" \
  --seed "$SEED" \
  --device "$DEVICE" \
  --json "$JSON_OUT" \
  --out "$FIG_OUT"

echo "RUN_DONE $(date -Is)" | tee -a "$STATUS"
