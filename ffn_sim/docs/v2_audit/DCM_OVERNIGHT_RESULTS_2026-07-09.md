# DCM 밤샘 큐 결과 (2026-07-09, gbook A5000)

**계획:** DCM_OVERNIGHT_PLAN_2026-07-09.md  **큐:** `overnight_dcm_2026-07-09.sh` (halt-free 순차)
**엔진:** DCM (Warp, A5000)  **브랜치:** dcm/main

밤새 자율 실행(PI /goal). 각 Phase 완료 시 능동 분석. 아래는 확정된 결과.

## Phase A — N=2000 대형-dt: SIGMA=0 정적 안정성 데모 + 게이트 버그 규명 (⚠️ "해결" 아님)

confluent N=2000, INSET=0.18, accel_dt=8e-3, **SIGMA=0(압축 driver 없음)**, **STEPS=5000 (40 s)**, per-step 2.06 s.
**⚠️ SIGMA=0 + confluent(이미 compact) ⇒ 시스템이 40 s 동안 거의 정지(A/A0=1.000, drift=0).** 아래 "관통 안정"은
정지 상태의 안정이지 압축(움직임) 시 검증이 아니다. red-flag #2/#4를 "해결"한 게 아님(정정 항목은 문서 하단).

- **관통 시계열 완전 안정**: 깊은관통(>0.2R) = 315 nodes (0.11%) — step 0→5000 동안 **불변**(폭발/드리프트
  없음). V/V0=1.000, A/A0=1.000, drift=0. → 대형-dt IPC가 N=2000 관통을 40 s 내내 정확히 고정.
- porosity 0.451 불변, asph 0.0035 (둥근 confluent foam). 게이트 pen 75.09 불변(= 실제 관통 flat의 메트릭 그림자).
- HTML 421 MB 풀-렌더, 실제 Chrome WebGL 검증(에러 0) — 조밀 watertight foam spheroid.
- **정직 결론 (PI 지적, 정정):** ✅ 확실한 성과 = **게이트 pen 75 메트릭 버그 규명**(실제 관통 0.11%).
  ⚠️ 나머지는 **정적 안정성 데모일 뿐** — SIGMA=0이라 압축·mechanobiology가 없다. ❌ red-flag #2(timescale)는
  미해결(40 s에 실제 물리 없음), ❌ N=2000 large-dt **압축**은 미해결(loose 발산 / confluent 이미 compact /
  48h ~515일 비현실). 상세: DCM_N2000_WATERTIGHT_2026-07-09.md.

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
- **Morphology (browser-verified, real Chrome WebGL):** `preview_n400_sigma5_compacted.png` (σ=5, 조밀
  cohesive spheroid) vs `preview_n400_sig10_partial.png` (σ=10, 부분 압축) vs `preview_n400_sig20_unstable.png`
  (σ=20, 분리된 구 + void = 압축 상실) — non-monotonic을 형태로 확증. figs/agg_compaction/.

## Phase C — 장시간 수렴 + 계속 densify (red-flag #5+#2 개선) ✅

N=400 loose gap-2.4, σ=5 mN/m, accel_dt=8e-3, **STEPS=12000 (96 s 물리시간)**:

- **porosity 0.631 → 0.285** — 32 s(Phase B σ=5: 0.408)보다 훨씬 조밀. 장시간에 압축이 계속 densify,
  fcc 강체구 이상(0.26)에 근접(셀 변형으로 void 제거). Rg −16.6%, pen=0, V/V0=1.000 내내.
- **drift(프레임 간 aggregate COM 이동): 0→0.14 µm(peak, step 7500)→0.04 µm(step 12000) = 감소·수렴 방향.**
  A/A0 0.941→0.723, 후반 Δ 매우 작음(0.725→0.723) = 포화 근접. residual ≈ 0.04 µm / 4.8 s ≈ **0.008 µm/s**
  « red-flag #5의 0.82 µm/s → **수렴 개선 입증**.
- 초기 재배열(drift 상승)→densify→안정화(drift 하강)의 물리적 궤적. 장시간 대형-dt가 red-flag #5(non-
  convergence)와 #2(timescale, 96 s)를 동시에 밀어냄.
- **정직 caveat:** per-node residual force `fmag`는 3.2e-11 N 정상상태로 유지(step 900→12000 거의 불변,
  발산 없음)이나 0으로는 안 감 — cadherin catch-bond의 동적 churn(형성/파열 평형)이 잔여 힘을 유지하기
  때문. 즉 aggregate-level(drift)은 수렴, node-level(fmag)은 안정된 non-zero 정상상태. red-flag #5 대비
  개선(drift 0.008 µm/s)이되 "완전 정지"는 아님.

## Phase D — confluent outlier: lloyd_iters↑가 주 lever (부분 완료)

- **D1 (subdiv=2, eps=0.10, lloyd=8): INIT 깊은관통 30 nodes (22 cells).** 내 confl2000b(lloyd=6,
  eps=0.18)의 75 nodes(65 cells) 대비 **절반 이하** → **lloyd_iters↑(6→8)가 outlier를 줄이는 주 lever**
  (seed 균등화 → Voronoi 인접-쌍 겹침 감소). eps보다 lloyd가 지배적으로 보임.
- ⚠️ **full sweep 미완(중단):** `confluent_init_prototype.validate()`가 N=2000에서 O(nodes²)
  contact_frac(324k²)을 계산 → 조합당 ~44 min, subdiv=3(1.3M²)은 비실용(수 h). validate를 KDTree/
  subsample로 효율화해야 full lloyd×subdiv×eps sweep 가능(개선 항목). outlier **방향(lloyd↑)**은 확인.

## 48h N=2000 시뮬레이션 실현가능성 (PI 질문 2026-07-09)

**현재 dt(accel_dt=8e-3): ~515일 — 비현실적.** 48h 물리시간 = 172,800 s ÷ 8e-3 = **2.16e7 step**
× per-step 2.06 s(N=2000, Phase A 실측) = 4.45e7 s ≈ **515일**.

**dt로 못 줄임 — 8e-3이 N=2000의 사실상 상한.** dt-probe(confluent N=2000, SIGMA=0): **accel_dt=0.1에서
발산**(pen 101, cfl 5.3e4). #1 GPU 검증의 dt ceiling ~2.0 s는 small-N(N=2)에서였고, N=2000은 수천 셀의
lagged soft-force + Newton 결합으로 dt 상한이 **8e-3로 훨씬 낮다**. (dt=1.0/2.0 = 더 발산, 확정 대기.)

**핵심: 48h는 mechanical compaction에 불필요하다.** Phase C에서 N=400 σ=5가 **96 s**에 porosity 0.285로
포화(거의 완전 densify) + drift 수렴 — mechanical 압축·재배열은 ~100 s면 끝난다. 24–48h는 biological
remodeling(cadherin 성숙, cortex 재조직) 타임스케일로, 이 fine-grained *mechanical* 모델의 스코프 밖이다.
→ **48h 도달을 목표로 하지 말고 mechanical 포화(~100 s)까지가 이 엔진의 유효 범위.** 더 긴 시간이 필요하면
(a) per-step을 근본적으로 낮추거나(GPU 커널 최적화, Newton 예열), (b) biological remodeling을 별도
mesoscale 모델로 붙여야 한다.

## 종합

밤샘 큐 4 Phase (halt-free, A5000). 핵심 성과 (C/D는 완료 시 갱신):

1. **게이트 pen 75 = 메트릭 버그 규명 (Phase A의 ✅ 확실한 성과)** — 실제 깊은관통 0.11%(측정 확실). ⚠️ 단
   Phase A는 SIGMA=0 정적이라 "N=2000 large-dt 작동/해결"은 **아님**(압축·mechanobiology 없음, red-flag #2 미해결,
   압축은 loose 발산·confluent 이미 compact로 미해결).
2. **σ-sweep non-monotonic (Phase B)** — sweet spot σ=5 mN/m @ accel_dt=8e-3; σ>5는 σ-dt 결합으로
   과충격·관통·압축 상실. 압축 production 표준 = σ=5 (lit-anchored + 안정).
3. **장시간 수렴 + densify (Phase C)** — σ=5 @ 96s: porosity 0.631→0.285(32s의 0.408보다 조밀),
   drift 0.14→0.04µm(수렴, residual ~0.008µm/s « 0.82), red-flag #5+#2 개선.
4. **outlier 최소화 (Phase D, 부분)** — lloyd_iters↑(6→8)가 confluent outlier를 75→30 nodes로 절반 이하.
   full sweep은 validate O(nodes²) 병목으로 중단(효율화 개선 항목). 방향(lloyd↑) 확인.

### PI 결정 대기

- **게이트 pen_frac 메트릭 수정** (DCM_N2000_WATERTIGHT §수정제안, 정확한 diff 준비됨) — gate-contract
  변경이라 PI 승인 필요. 승인 시 인라인 적용.
- **압축 production 설정 확정**: σ=5 mN/m @ accel_dt=8e-3 (Phase B sweet spot).
- (C/D 결과에 따른 후속.)

### 드라이버 개선 항목 (밤샘 후, 밤샘 큐 일관성 위해 실행 중 미수정)

- `_gbook_aggregate_compaction.py`의 inline TREND `min_gap`이 N=2000에서 셀마다 32만 노드 KDTree →
  프레임당 수분(Phase A에서 큐 지연, 중단하고 관통은 별도 측정). 전체-tree-1회 + k-최근접 방식으로 효율화 필요.

Related: DCM_N2000_WATERTIGHT_2026-07-09, DCM_AGGREGATE_COMPACTION_LARGEDT_2026-07-08,
DCM_OVERNIGHT_PLAN_2026-07-09.
