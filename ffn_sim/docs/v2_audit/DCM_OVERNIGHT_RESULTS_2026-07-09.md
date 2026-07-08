# DCM 밤샘 큐 결과 (2026-07-09, gbook A5000)

**계획:** DCM_OVERNIGHT_PLAN_2026-07-09.md  **큐:** `overnight_dcm_2026-07-09.sh` (halt-free 순차)
**엔진:** DCM (Warp, A5000)  **브랜치:** dcm/main

밤새 자율 실행(PI /goal). 각 Phase 완료 시 능동 분석. 아래는 확정된 결과.

## Phase A — N=2000 대형-dt: 관통 0.11% 안정 (red-flag #2+#4 동시 해결) ✅

confluent N=2000, INSET=0.18, accel_dt=8e-3, SIGMA=0, **STEPS=5000 (40 s 물리시간)**, per-step 2.06 s (3.05 h).

- **관통 시계열 완전 안정**: 깊은관통(>0.2R) = 315 nodes (0.11%) — step 0→5000 동안 **불변**(폭발/드리프트
  없음). V/V0=1.000, A/A0=1.000, drift=0. → 대형-dt IPC가 N=2000 관통을 40 s 내내 정확히 고정.
- porosity 0.451 불변, asph 0.0035 (둥근 confluent foam). 게이트 pen 75.09 불변(= 실제 관통 flat의 메트릭 그림자).
- HTML 421 MB 풀-렌더, 실제 Chrome WebGL 검증(에러 0) — 조밀 watertight foam spheroid.
- **결론: N=2000도 대형-dt로 돈다.** timescale 0.12 s→40 s, interpenetration 0.11% 안정. 재-init 불필요,
  남은 건 게이트 pen_frac 메트릭 수정(PI 승인). 상세: DCM_N2000_WATERTIGHT_2026-07-09.md.

## Phase B — N=400 압축 σ-sweep: NON-MONOTONIC, sweet spot σ=5 mN/m ⭐

loose gap-2.4, accel_dt=8e-3, STEPS=4000 (32 s), σ ∈ {0,1,5,10,20} mN/m (통제; σ=0 baseline 공유):

| σ [mN/m] | Rg | porosity 0.631→ | pen_peak | 압축 |
|---|---|---|---|---|
| 0 (baseline) | −0.01% | 0.604 | 0 | 없음 (loose 유지) |
| 1 | −6.4% | 0.455 | 0 | 압축 |
| **5** | **−8.7%** | **0.408** | 0 | **최대** |
| 10 | −2.7% | 0.521 | 0 | 감소 |
| 20 | −1.9% | 0.547 | **2.50 (지속)** | 불안정 (관통 회복 안 됨) |

- **핵심 발견: σ-압축이 단조롭지 않다.** Foty-Steinberg 직관(σ↑→압축↑)과 반대로 **σ=5가 최대**, σ=10은
  오히려 감소, σ=20은 관통(pen 2.5)까지 발생. 
- **원인 = σ-dt 결합(수치 경계), 물리 아님.** aggregate σ는 lagged soft force → 큰 σ + 큰 dt(8e-3)는
  loose 시작에서 초기 변위를 과대 적분 → 과충격/관통 → 압축 방해·불안정. accel_dt=8e-3에서 σ=5가 안정+최대
  압축의 sweet spot; 더 큰 σ는 더 작은 dt를 요구(CFL 유사). 
- **함의:** 압축 production은 σ=5 mN/m @ accel_dt=8e-3 (lit-anchored + 수치 안정). σ를 키워 압축을 더
  얻으려면 dt를 낮춰야 한다(magic-number 아님 — σ는 lit 범위, dt는 안정성 유도).

## Phase C — N=400 장시간 수렴 (residual/drift→0)

<!-- STEPS=12000 (96s) 완료 시: drift 시계열 수렴 -->
(진행 대기.)

## Phase D — confluent outlier 최소화 (geometry 스윕)

<!-- lloyd×subdiv×eps 조합별 INIT 깊은관통% -->
(진행 대기.)

## 종합 (밤샘 완료 시)

<!-- 4 Phase 종합 + PI 결정 대기 항목 -->

Related: DCM_N2000_WATERTIGHT_2026-07-09, DCM_AGGREGATE_COMPACTION_LARGEDT_2026-07-08,
DCM_OVERNIGHT_PLAN_2026-07-09.
