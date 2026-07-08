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

## 실측 (엔진 커널 재사용, `ffn_sim/scripts/_gbook_measure_penetration.py`)

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

권장: (2)+(3) 조합 — `frac(pen>0.2R)`를 임계 ~0.5%로. 이러면 IPC-newton(0.1%) PASS, baoab(2.1%) FAIL,
물리적으로 옳음. **PI 승인 후** gate-contract 갱신 + `_penetration_frac` 수정.

## 대형-dt 이득 (긴 run)

<!-- 긴 run(600 step, 물리시간 4.8s) 결과: 관통 시계열 안정성 + per-step wall 채움 -->
(진행 중.)

## 남은 것

- 0.1% outlier(65→196 cells, 이산 icosphere warp의 국소 face 겹침)를 더 줄이려면 `lloyd_iters`↑ 또는
  `subdiv`↑(기하 정확도, outcome-튜닝 아님) — 다만 0.1%는 이미 thesis-grade(baoab 2.1-37.6% 대비). 필수 아님.
- 게이트 수정은 PI 승인 사항.

## Files
- `ffn_sim/scripts/_gbook_measure_penetration.py` — 엔진 커널 재사용 절대-관통 측정.
- 실측 npz: `~/ff_scratch/_prod_out/agg_compaction_n2000_*_confl2000b.npz`, `fullcomp_assembly_n2000_n2000.npz`.

Related: `N2000_CONFLUENT_FINDINGS_2026-07-07.md`(pen 30× 과대정규화 최초 지적),
`DCM_IPC_LARGE_DT_PLAN_2026-07-07.md`(#1), memory [[project-dcm-fidelity-remediation]] red-flag #4.
