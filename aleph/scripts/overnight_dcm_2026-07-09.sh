#!/usr/bin/env bash
# DCM 밤샘 큐 2026-07-09 — gbook A5000 순차 실행. halt-free(각 run 실패해도 다음 진행), 각 자체 로그.
# 계획: aleph/docs/v2_audit/DCM_OVERNIGHT_PLAN_2026-07-09.md
# 실행:  ssh gbook 'nohup bash ~/ff_scratch/aleph/scripts/overnight_dcm_2026-07-09.sh >/dev/null 2>&1 &'
# 산출:  npz(각 run) + 단계 로그; 측정/시각화는 아침에 Lead가.
set -uo pipefail
OUT=~/ff_scratch/_prod_out
PY=~/miniconda3/envs/ffn_sim/bin/python
cd ~/ff_scratch
export PYTHONPATH=~/ff_scratch
LOG=$OUT/overnight_master.log
echo "OVERNIGHT START $(date)" | tee "$LOG"
DRV=aleph/scripts/_gbook_aggregate_compaction.py

run() { echo "[$(date +%H:%M)] $1" | tee -a "$LOG"; }

# ── Phase A: N=2000 대형-dt 장시간 (timescale #2 + interpenetration #4) ──────────────
run "A  N=2000 confluent large-dt STEPS=5000 (40s phys)"
BUILDER=confluent NCELLS=2000 INSET=0.18 WARMUP=300 STEPS=5000 ACCEL_DT=8e-3 SIGMA=0 \
  FRAMES=40 TAG=_ovnA $PY $DRV > "$OUT/ovn_A_n2000.log" 2>&1 || run "A FAILED"

# ── Phase B: N=400 압축 σ-sweep (Foty-Steinberg 물리; σ=0 baseline 포함) ─────────────
for SIG in 0 1e-3 5e-3 1e-2 2e-2; do
  run "B  N=400 loose sigma=$SIG STEPS=4000"
  NCELLS=400 WARMUP=500 STEPS=4000 ACCEL_DT=8e-3 SIGMA=$SIG GAP=2.4 FRAMES=30 TAG=_ovnB \
    $PY $DRV > "$OUT/ovn_B_sig${SIG}.log" 2>&1 || run "B $SIG FAILED"
done

# ── Phase C: N=400 장시간 수렴 (convergence #5) ──────────────────────────────────────
run "C  N=400 loose sigma=5e-3 STEPS=12000 (96s phys)"
NCELLS=400 WARMUP=500 STEPS=12000 ACCEL_DT=8e-3 SIGMA=5e-3 GAP=2.4 FRAMES=40 TAG=_ovnC \
  $PY $DRV > "$OUT/ovn_C_n400long.log" 2>&1 || run "C FAILED"

# ── Phase D: confluent 초기화 outlier 최소화 (geometry-only build + 관통 측정) ─────────
for SUB in 2 3; do for EPS in 0.10 0.18 0.25; do
  run "D  confluent build subdiv=$SUB eps=$EPS + pen measure"
  $PY aleph/dcm/confluent_init_prototype.py --n 2000 --subdiv "$SUB" --eps "$EPS" \
     --out "$OUT/ovn_D_confl_sub${SUB}_eps${EPS}.npz" > "$OUT/ovn_D_sub${SUB}_eps${EPS}.log" 2>&1 \
     || { run "D $SUB $EPS BUILD FAILED"; continue; }
  NPZ="$OUT/ovn_D_confl_sub${SUB}_eps${EPS}.npz" RADIUS_R=1.5 \
     $PY aleph/scripts/_gbook_measure_penetration.py >> "$OUT/ovn_D_sub${SUB}_eps${EPS}.log" 2>&1 || true
done; done

echo "OVERNIGHT DONE $(date)" | tee -a "$LOG"
