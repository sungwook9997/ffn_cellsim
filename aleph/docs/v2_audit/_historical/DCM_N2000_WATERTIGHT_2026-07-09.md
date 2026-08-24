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

# N=2000 대형-dt: 이미 사실상 watertight — "pen 75"는 게이트 메트릭 버그였다 (2026-07-09)

**Engine:** DCM (Warp, A5000)  **Branch:** dcm/main  **Follows:** #1 large-dt IPC + N=64/400 압축 성공(61978fc)
**PI 선택(2026-07-08):** N=2000 watertight 재-init을 판다.

## TL;DR (반전)

"N=2000을 watertight하게 재-init해야 한다"의 전제였던 **pen_frac≈75 (G2 FAIL)는 실제 관통이 아니라 게이트
메트릭의 버그**였다. 실측(엔진 자체의 `penetration_depth_kernel`, Ericson closest-point) 결과:

| N=2000 (실측, radius 1.5R) | 깊은관통(>0.2R) 노드 | 비율 | 현재-mean_edge 정규화 |
|---|---|---|---|
| **대형-dt IPC-newton** (INSET 0.18) | 300 / 324,000 | **0.09%** | 5.7 |
| 기존 baoab (bundle 100) | 6,800 / 324,000 | 2.1% | 4.6 |

**대형-dt IPC-newton N=2000은 노드의 99.9%가 관통 0인 사실상 watertight 상태**이고, 기존 baoab 대비
깊은관통을 **20× 줄였다**(fidelity red-flag #4를 사실상 해결). 게이트가 뱉는 75는 실제 관통이 아니라
(a) **초기 build의 mean_edge로 정규화**(settled/변형 mesh보다 ~30× 작음, INSET build는 특히) + (b)
**최악 단일 노드(max)만** 보기 때문에 생기는 아티팩트다. → **실제 필요한 작업은 재-init이 아니라 게이트
pen_frac 메트릭 수정**(PI-gated gate-contract 변경)이다.

## 조사 (workflow understand + 코드 확인)

1. **inset은 작동하고 build 자체는 (거의) watertight.** `warp_icosphere_to_voronoi`
   (`confluent_init_prototype.py:57-78`)는 각 셀 icosphere 노드를 seed에서 Voronoi 경계까지 `t`,
   ×`(1-eps)`로 수축한다. `validate()`(`:132-141`)는 `d_other/d_own>1` (노드가 자기 Voronoi 셀 안 =
   멤버십)로 비겹침을 확인 → eps>0이면 통과. INSET=0.05↔0.18로 게이트 pen이 안 변한 건 inset이 죽어서가
   아니라, pen이 (아래) 메트릭 아티팩트 + cadherin 평형이라 초기 gap에 둔감하기 때문.
2. **그러나 validate의 "비겹침"은 노드 멤버십이지 face-face 겹침이 아니다.** 이산 icosphere를 볼록 Voronoi
   다면체로 warp하면 삼각형 face가 이웃 셀 경계를 국소적으로(꼭짓점/모서리 근처, 이웃이 많은 곳) 넘어갈 수
   있다. 실측: INIT에서 깊은관통 75 노드 / 65 cells (324,000 중 0.023%) — build는 99.98% watertight,
   극소수 국소 겹침만.
3. **게이트 pen_frac = max_pen / mean_edge (초기 1회 계산).** `_penetration_frac`
   (`dcm_warp_decohesion.py:1363-1377`)는 `penetration_depth_kernel`의 per-node 깊이 중 **max**를
   `mean_edge`(`:391`, 초기 build edge)로 나눈다. INSET로 셀이 수축된 build는 초기 edge가 작아 분모가
   작음 → 비율 뻥튀기. N2000_CONFLUENT_FINDINGS(2026-07-07)가 이미 "~30× 과대정규화"로 지적한 값.

## 실측 (엔진 커널 재사용, `aleph/scripts/_gbook_measure_penetration.py`)

`penetration_depth_kernel`을 confl2000b(내 IPC-newton) vs fullcomp_assembly_n2000(기존 baoab) npz에
직접 적용, mean_edge 정규화 없이 절대 관통 깊이(µm, R 단위) 분포를 계산:

- **IPC-newton**: INIT 깊은관통 75 nodes(65 cells) → FINAL 300 nodes(196 cells) = 0.1%. 나머지 99.9%
  pen=0. cadherin(bundle 10)이 이웃 막 너머로 mutual-nearest 노드를 약간 당겨 75→300 증가(SIGMA=0이라
  압축 driver 없음에도).
- **baoab**: FINAL 깊은관통 6,800 nodes = 2.1%, 얕은관통(inside) 26.8%. penalty contact가 관통을 못 막음.

## 게이트 pen_frac 수정 제안 (PI-gated — no-gate-loosening 규칙상 인라인 수정 금지, surface to PI)

현재 메트릭은 **물리적으로 틀렸다**(초기 mean_edge는 settled mesh와 무관, max-only는 분포 무시). 이건
게이트를 통과시키려 느슨하게 하는 게 아니라 **거짓 FAIL을 내는 메트릭을 고치는 것**. 세 가지 중 택(또는 조합):

1. **정규화를 현재-프레임 mean_edge로** — `_penetration_frac`이 측정 시점의 mean_edge를 쓰게. (최소 변경)
2. **max → 분포 지표** — `frac(pen > 0.2R)` (깊은관통 노드 비율) 또는 95퍼센타일 깊이. 단일 outlier가
   게이트를 결정하지 않게.
3. **물리 임계(µm/R)** — `max_pen`을 R 대비 절대값으로(예: median_inside < 0.05R AND frac>0.2R < 0.5%).

권장: (2)+(3) 조합 — `frac(pen>0.2R)`를 임계 ~0.5%로. **PI 승인 후** gate-contract 갱신 + `_penetration_frac` 수정.

**⚠️ 적용 후 런타임 검증 (2026-07-09 gate-verify, 새 코드 confluent N=2000 SIGMA=0, dt 8e-3/0.1/1.0):**
- ✅ 게이트 수정 **작동**: 로그 pen이 이제 frac_deep(0.002~0.017), 구 75(메트릭 버그) 사라짐.
- ⚠️ **그러나 confluent N=2000은 여전히 G2 FAIL**: 게이트가 `pen_frac_peak`(궤적 max) 기준인데 confluent build의
  **warmup 과도기 관통 peak = 1.7% > 0.5% 임계**. `pen_final`(0.3%)만 낮다. → 내 "IPC-newton(0.1%) PASS"
  예상은 **final 기준 착오**였고, PEAK가 걸려 FAIL. (dt 무관 — build/warmup 시점 관통이라 8e-3/0.1/1.0 모두 peak 1.7%.)
- **남은 이슈 (PI 결정):** (a) 게이트가 warmup 과도기까지 pen_peak로 잡는 게 맞나(build 초기 관통 vs 운영 중),
  (b) confluent build의 초기 관통(outlier)을 Phase D의 lloyd_iters↑로 peak도 낮출 수 있는가. 게이트 메트릭 자체는
  이제 실제 관통을 잰다(수정 목표 달성), 하지만 confluent N=2000이 게이트를 통과하려면 초기 관통을 더 줄여야 함.

**PI 승인 시 즉시 적용할 정확한 diff** (`dcm_warp_decohesion.py:1377`, 스코프에 `R=p.R_cell`·`edges_a`·
`pos_d` 존재 확인됨):

```python
# 현재 (버그): 초기-build mean_edge(INSET-수축 → ~30× 작음)로 max-only 정규화
    return float(pen_d.numpy().max()) / mean_edge
# 제안 (권장 옵션 2+3): 현재-프레임 mean_edge + 깊은관통 분포 지표
    pen = pen_d.numpy()                                        # per-node 관통 깊이 [m]
    P = pos_d.numpy()
    me_now = float(np.linalg.norm(P[edges_a[:, 0]] - P[edges_a[:, 1]], axis=1).mean())
    frac_deep = float((pen > 0.2 * R).mean())                 # 깊은관통(>0.2R) 노드 비율
    return frac_deep     # 게이트 임계를 0.3 → 0.005 로 (gate-contract 변경, PI 승인)
```

임계를 바꾸므로 `trajectory`의 `pen_frac` 소비처(게이트 판정 0.3)도 0.005로 동반 갱신해야 한다 —
그래서 **인라인 수정이 아니라 gate-contract 변경(PI 승인)**이다. 옵션 1만(최소변경)이면 `return
float(pen.max())/me_now` 한 줄이지만 max-only 문제는 남는다.

## 대형-dt 관통 안정성 (긴 run, 밤샘 Phase A)

N=2000 confluent INSET=0.18, accel_dt=8e-3, SIGMA=0, **STEPS=5000 (40 s 물리시간)**, 40 프레임.
관통 시계열(`_gbook_measure_penetration.py ALL_FRAMES=1`, radius 1.5R):

| 프레임 | step | 깊은관통(>0.2R) | 비율 | max | median |
|---|---|---|---|---|---|
| f00 (warmup 전) | 0 | 75 nodes (65 cells) | 0.023% | 1.26R | 0.28R |
| f01 (warmup 후) | 0 | 315 nodes (207 cells) | 0.110% | 1.32R | 0.40R |
| **f02 … f41** | 125 … 5000 | **315 nodes (207 cells) — 완전 불변** | **0.110%** | 1.32R | 0.40R |

- **대형-dt IPC는 N=2000 관통을 40 s 물리시간 내내 정확히 안정 유지** — step 0→5000 동안 깊은관통이
  315 nodes(0.11%)로 **완전 불변**(증가도 감소도 없음). 폭발/터널링 없음. V/V0=1.000, A/A0=1.000, drift=0.
- warmup(작은 dt)에서 75→315로 늘어난 건 cadherin(bundle 10)이 이웃 막 너머 mutual-nearest 노드를 당기는
  평형(inset-path 조사와 일치); 메인 루프에서 IPC가 그 평형을 그대로 고정.
- **per-step ≈ 2.06 s** (5000 step + 300 warmup, 3.05 h wall). 게이트 pen은 이 내내 75.09 불변(= 실제
  관통 flat의 메트릭 그림자).

**⚠️ 정직 정정 (2026-07-09, PI 지적):** 이 Phase A는 **SIGMA=0(압축 driver 없음) + confluent(이미 compact)**
이라 40 s 동안 시스템이 거의 **정지**(A/A0=1.000, drift=0)했다. 따라서 이건 "N=2000 large-dt가 해결됐다"가
**아니다**:
- ✅ **확실**: 게이트 pen 75 = 메트릭 버그 규명(실제 깊은관통 0.11%, 측정 확실).
- ⚠️ **정적 안정성 데모일 뿐**: 정지 confluent이 대형-dt에서 안 터짐. **압축(움직임) 시 관통이 안 생기는지는 미검증.**
- ❌ **red-flag #2(timescale) 미해결**: 40 s를 "도달"했으나 그 시간에 mechanobiology(압축·재배열)가 없음 →
  timescale의 핵심(min-hr에 걸친 실제 물리)을 충족 못 함. #4도 정적에서만 안정.
- ❌ **N=2000 large-dt 압축(loose→compact) 여전히 미해결**: loose N=2000은 대형-dt 발산, confluent은 이미
  compact(압축할 void 없음), 48h는 ~515일 비현실(dt=0.1 발산→8e-3 상한). N=400 압축(Phase B/C)만 성립.

## 남은 것

- 0.1% outlier(65→196 cells, 이산 icosphere warp의 국소 face 겹침)를 더 줄이려면 `lloyd_iters`↑ 또는
  `subdiv`↑(기하 정확도, outcome-튜닝 아님) — 다만 0.1%는 이미 thesis-grade(baoab 2.1-37.6% 대비). 필수 아님.
- 게이트 수정은 PI 승인 사항.

## Figures (밤샘 Phase A)

- `outputs/h_dcm_two_stage/figs/agg_compaction/n2000_ovnA_largedt.html` — N=2000 대형-dt IPC 40s 형태
  뷰어(2000 cells, 324k nodes, 42 frames, 풀 렌더 421MB; virial σ_vm peak 340 Pa). 실제 Chrome ANGLE
  Metal WebGL 검증 통과(에러 0). 조밀한 watertight confluent foam spheroid — 대형-dt 내내 안정.
  프리뷰: `preview_n2000_ovnA_largedt.png`. (421MB HTML은 로컬만, npz에서 재생성.)

## Files
- `aleph/scripts/_gbook_measure_penetration.py` — 엔진 커널 재사용 절대-관통 측정.
- 실측 npz: `~/ff_scratch/_prod_out/agg_compaction_n2000_*_confl2000b.npz`, `fullcomp_assembly_n2000_n2000.npz`.

Related: `N2000_CONFLUENT_FINDINGS_2026-07-07.md`(pen 30× 과대정규화 최초 지적),
`DCM_IPC_LARGE_DT_PLAN_2026-07-07.md`(#1), memory [[project-dcm-fidelity-remediation]] red-flag #4.
