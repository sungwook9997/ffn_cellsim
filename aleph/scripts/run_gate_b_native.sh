#!/usr/bin/env bash
# H.7 Gate-B full ×40 native production — CONNECTED-MESH (the 2026-06-04 rebuild,
# REQUIRED for a transmission measurement; Gate-A/old runs wrongly used the
# fragmented mesh). 4 conditions sequential, checkpointed, then aggregate.
# NOTE: no `set -u` — conda's activate functions reference unset vars and would abort.
cd ~/ffn_cellsim
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate ffn_sim 2>/dev/null
export PYTHONPATH="$PWD/native/ffn_hoomd_plugin/build:$PWD/native/ffn_hoomd_plugin/python:$PWD"
export FFN_GPU_DEVICE_COMPARTMENTS=1

TAG="${1:-gateB_cm}"
DTS="${2:-1}"          # dt safety — keep 1× (s_grip development is step-bound; dt>1 under-develops it)
ON_TICKS="${3:-400}"   # myoON: s_grip develops to its plateau (~tick 350 in Gate-A)
OFF_TICKS="${4:-80}"   # myoOFF: passive baseline settle
NUC="${5:-no}"         # 'no' → --no-nucleus (validated decoupled from cortical γ); 'yes' → keep
NUCFLAG=""; [ "$NUC" = "no" ] && NUCFLAG="--no-nucleus"
# connected_mesh is ON by default in the driver (do NOT pass --no-connected-mesh).

P="python -u -m aleph.scripts.h7_gate_b_native --tag $TAG --n-fil 1000 --interval 5000 \
   --seed 1 --measure-every 5 --ckpt-every 10 --dt-safety $DTS $NUCFLAG"

echo "=== Gate-B native (CONNECTED MESH) START $(date) tag=$TAG dt=$DTS nuc=$NUC ON=$ON_TICKS OFF=$OFF_TICKS ==="
$P --mode rigid   --myo off --ticks $OFF_TICKS || { echo "FAILED rigid_myoOFF"; exit 1; }
$P --mode relaxed --myo off --ticks $OFF_TICKS || { echo "FAILED relaxed_myoOFF"; exit 1; }
$P --mode rigid   --myo on  --ticks $ON_TICKS  || { echo "FAILED rigid_myoON"; exit 1; }
$P --mode relaxed --myo on  --ticks $ON_TICKS  || { echo "FAILED relaxed_myoON"; exit 1; }

echo "=== aggregating $(date) ==="
python -m aleph.scripts.h7_gate_b_native --aggregate --tag $TAG
echo "=== Gate-B native DONE $(date) ==="
