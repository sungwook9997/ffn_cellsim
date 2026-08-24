---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# DCM 밤샘 실행 계획서 (2026-07-09 야간, gbook A5000)

**작성:** 2026-07-09  **엔진:** DCM (Warp, A5000, gbook)  **브랜치:** dcm/main
**근거:** 이번 세션 성과 — (1) aggregate-σ 압축 N=400 성공(61978fc), (2) N=2000은 이미 사실상
watertight·게이트 pen_frac이 버그(DCM_N2000_WATERTIGHT_2026-07-09), (3) #1 대형-dt IPC.

## 목표 (fidelity red-flag 정면 공략)

밤새 A5000 한 대로 순차 큐를 돌려, 낮에 손대던 세 red-flag를 **긴 물리시간·통제 스윕**으로 밀어붙인다:
- **#2 timescale** — 대형-dt로 N=2000·N=400을 초→분 규모 물리시간까지.
- **#4 interpenetration** — N=2000 관통이 장시간 0.1%로 안정한지(폭발 안 하는지) 시계열.
- **#5 convergence** — 장시간 run에서 residual·drift가 0으로 수렴하는지.
- **물리 검증** — N=400 압축의 Foty-Steinberg σ-의존성(1–20 mN/m)이 단조·물리적인지.

## 큐 (순차, 각 nohup, 산출물 = npz; 측정·시각화는 아침에)

per-step 실측 근거: N=400 σ런 ≈0.45 s/step(1824 s/4000), N=2000 ≈1.5 s/step(추정, cg 350–404).

| # | Run | 설정 | 물리시간 | ETA(추정) | 검증 대상 |
|---|---|---|---|---|---|
| A | **N=2000 대형-dt 장시간** | confluent, INSET=0.18, accel_dt=8e-3, SIGMA=0, STEPS=5000, FRAMES=40 | 40 s | ~2.2 h | #2·#4 (관통 시계열, per-step) |
| B1–B5 | **N=400 압축 σ-sweep** | loose gap-2.4, accel_dt=8e-3, STEPS=4000, σ∈{0,1,5,10,20} mN/m | 32 s×5 | ~2.5 h | 물리 (Foty-Steinberg 단조성) |
| C | **N=400 장시간 수렴** | loose gap-2.4, σ=5, accel_dt=8e-3, STEPS=12000, FRAMES=40 | 96 s | ~1.5 h | #5 (residual·drift→0) |
| D | **confluent outlier 최소화** (geometry-only, physics 없음) | build_confluent 스윕 lloyd_iters∈{6,12,20}×subdiv∈{2,3}×eps∈{0.10,0.18,0.25} | — | ~20 min | #4 (INIT 관통 최소 조합) |

**총 ETA ≈ 6.5 h.** A→B→C→D 순차. 각 단계 자체 로그 파일 + 산출 npz.

## 자율 실행 원칙 (memory: 8h-autonomous-mandate, no-param-tuning, physiological-baseline)

- **halt 금지, document+continue.** 한 run이 발산/실패하면 로그에 기록하고 다음 run으로. 큐 전체가 멈추지 않게.
- **magic-number 금지.** σ는 lit Foty-Steinberg [1,20] mN/m만; accel_dt=8e-3은 pen/G2 게이트에서 유도된
  수치 상한(튜닝 아님); INSET=0.18은 feasible-start용 기하값; nucleus/turgor/η는 생리값 고정.
- **gate-loosening 금지.** 게이트 pen_frac 수정은 **PI 승인 대기** — 밤새 자율로 코드 수정 안 함. D는 게이트를
  건드리지 않고 초기 기하만 스윕(측정).
- **생리 baseline.** 모든 run이 full 컴파트먼트(nucleus E_nuc=399, turgor+K_vol, η=65.9, cortex γ) 켠 채.
- 아침에 Lead(나)가 각 npz에 `_gbook_measure_penetration.py`(관통 시계열) + `dcm_mesh_viewer_html.py`
  (HTML 형태 뷰어, 브라우저 검증) + porosity/Rg 분석을 돌려 종합.

## 아침 체크리스트 (수집 항목)

1. **A**: N=2000 관통 시계열(프레임별 깊은관통%) — 0.1%로 안정 vs 증가? per-step wall 확정. HTML 뷰어.
2. **B**: σ vs (porosity↓, Rg↓) 곡선 — 단조? lit 범위서 물리적 크기? (σ=0 baseline 대비 각 σ의 압축 delta).
3. **C**: residual·drift 시계열 — 0 수렴? (수렴 = thesis-grade 근거; 미수렴이면 원인 문서화).
4. **D**: 초기 깊은관통% 최소 (lloyd/subdiv/eps) 조합 — 0.1% outlier를 더 낮추는 기하 세팅.
5. 실패한 run 로그 + 원인.

## PI 결정 대기 (밤새 자율로 안 하는 것)

- 게이트 pen_frac 메트릭 수정(DCM_N2000_WATERTIGHT §수정제안 1–3) — gate-contract 변경이라 PI 승인 필수.
- subdiv=3 native(1.3M 노드)로의 확장 — 메모리/시간 큼, PI 판단.

## 실행

`aleph/scripts/overnight_dcm_2026-07-09.sh` (gbook, 순차 nohup). 승인 시 시작.

Related: DCM_N2000_WATERTIGHT_2026-07-09, DCM_AGGREGATE_COMPACTION_LARGEDT_2026-07-08,
DCM_IPC_LARGE_DT_PLAN_2026-07-07, memory [[project-dcm-fidelity-remediation]].
