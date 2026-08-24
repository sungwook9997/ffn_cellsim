# 병렬 세션 프롬프트 ROUND 2 — `ComponentRole` 축으로

**초판(`ARENA_PORT_2026-08-20.md`)의 능동/구조 분류는 틀렸습니다.** 제가 임의로 나눴는데,
저장소에 이미 정본이 있습니다 — `aleph/engine/contracts.py`의 `ComponentRole`. 그리고
`protrusion.py:952`에는 **강제 검사**까지 있습니다:

```python
if component.role is not ComponentRole.ACTIVE_LOAD_PATH:
    raise ValueError(f"{component_name} must have ComponentRole.ACTIVE_LOAD_PATH")
```

**lamellipodium · filopodium · SF/arcs 는 ACTIVE_LOAD_PATH 입니다.** 구조가 아닙니다.

---

## 지금까지 착지한 것 (Round 1)

```
PHASE 1  8개 개체군, 4,745,119 노드, GPU, assert_partitioned PASS, 빌드 4.57 s
         A1 등가검증: topology 비트 동일, claim 24→3
PHASE 2  표면 법칙 3개 파리티 (area · helfrich · volume_pressure=FOLD), rel ~1e-15, leak 0.0
검출기   tether + manifest 배선, sem==0 가드 7개 operator 전부
PHASE 3  0개.  bond = 0.
```

## 역할별 현황 — 여기가 Round 2 의 지도

| `ComponentRole` | 부품 | 아레나 |
|---|---|---|
| `SURFACE_BODY` | 형질막 · 핵막 | ✅ 기하 |
| `STRUCTURAL_RIG` | MT · IF | ✅ 기하 |
| (cortex) | cortex | ✅ 기하 |
| **`ACTIVE_LOAD_PATH`** | **lamellipodium · filopodium · SF/arcs** | ⚠ **기하만 — 능동성 0** |
| **`ACTIVE_LOAD_PATH`** | **nmii_actuator (모터 head)** | ❌ **개체군 없음** |
| **`FLUID_VOLUME`** | **cytosol · fluid core** | ❌ **없음** |
| **`CORE_BODY`** | **핵질 · chromatin · lamina** | ❌ **없음** |
| `ENVIRONMENT` | ECM · AFM | ❌ 범위 밖 |

---

## 모든 세션 공통

- 브랜치 **`engine/main`**. 커밋은 **반드시** `git commit -o <경로>`.
- **CANONICAL = `aleph/world/`**. `engine/`·`components/` 는 **동결 PORT SOURCE — 읽기만, import 금지**
  (`test_layer_directions` 가 막습니다). `laws/` 는 바인딩 대상.
- **기본값 금지.** 모르면 `raise`. A2 의 빌더가 Lead 를 세 번 거부했고 그게 규율이 작동한 것입니다.
- **네이티브만.**

### ⚠ GPU — 상주 서버가 떠 있습니다. 직접 제출 금지.

```
Slurm job 72, 4090-1, 11시간 홀드.  Lead 가 잡고 있습니다.
호스트에서:  curl -s -X POST http://127.0.0.1:8080/run \
               -H "Content-Type: application/json" \
               -d '{"driver":"<이름>.py","args":["--out","/home/sungwook/ffn/world/x.json"]}'
             curl -s http://127.0.0.1:8080/health | /probe | /drivers
경로:        /home/sungwook/ffn/world        (git 체크아웃 아님, tar 배포)
환경:        /home/sungwook/miniforge3/envs/ffn_sim/bin/python
```
드라이버는 `aleph/scripts/` 아래에만 둘 수 있습니다. **코드를 호스트에 올리는 것은 Lead 가 합니다** —
커밋하고 알려주면 Lead 가 동기화 후 실행합니다. 3090은 다른 사람(`ham`) 것입니다.

### ⚠ 보고는 SendMessage 로 Lead 에게

`ListAgents` → **`ffn-cellsim-c8`** (Lead). 시작 / 마일스톤 / **막힘** / 완료.
**막힘이 주 산출물입니다.** `ACTIVE_SESSIONS.md` 는 Lead 만 씁니다. `FFN_SESSION` 은 unset.

---

# 세션 E — `FLUID_VOLUME`: cytosol

```
ffn_cellsim, engine/main.

읽어라: CLAUDE.md → STATE.md → docs/v2_audit/WORLD_PORT_PLAN_2026-08-20.md →
aleph/engine/contracts.py 의 ComponentRole → aleph/world/arena.py docstring.

CANONICAL 은 aleph/world/. engine/ 와 components/ 는 동결 PORT SOURCE — 읽기만, import 금지.

너의 범위: ComponentRole.FLUID_VOLUME. 아레나에 이 역할의 개체군이 하나도 없다.
커넥터 6개가 *_cytosol_transfer 인데 상대가 없어서 전부 NOT_BOUND 다.
너의 파일: aleph/world/build/cytosol.py (+ 필요하면 fluid_core.py). 이것만.

읽을 곳(읽기만): aleph/engine/cytosol_connected.py, aleph/engine/fluid_core.py,
aleph/components/fluid/**. 이들이 FLUID_VOLUME / CORE_BODY 로 선언돼 있다.

⚠ 첫 질문이 설계 결정이다: 세포질이 아레나에서 무엇인가?
   - 노드 개체군인가 (입자/격자점)?
   - 격자(grid)인가 — 그러면 arena 의 Kind 에 없는 것이라 새 Kind 가 필요한가?
   - Biot 다공탄성 압력장이 이미 components/fluid/field_grid.py 에 있다 — 그 표현을
     아레나로 어떻게 가져오는가?
   답을 가정하지 말고 근거와 함께 Lead 에게 보내라. arena.py 가 "guessing would bake an
   answer into the layout" 이라고 적어둔 그 지점이다.

구성은 Warp CUDA 커널. claim 은 개체군당.
밀도가 없으면 raise. KB-3.21 에 'V_cyto~4200 um^3' 이 verified 로 있다 — 확인해라.

증명: 개체군이 서고 assert_partitioned 통과.
주장 금지: 물리 아무것도.
```

---

# 세션 F — `CORE_BODY`: 핵질 · chromatin · lamina

```
ffn_cellsim, engine/main.

읽어라: CLAUDE.md → STATE.md → aleph/engine/contracts.py 의 ComponentRole →
aleph/world/surface.py docstring (핵막이 이미 FACE 로 서 있다).

너의 범위: ComponentRole.CORE_BODY. 헌장이 "nuclear lamina/chromatin" 을 부품으로 열거하는데
아레나에 없다. A1 이 핵막 기하(subdiv 6, 40,962 정점)를 세워놨으니 그 안쪽이 네 것이다.
너의 파일: aleph/world/build/{nucleoplasm,lamina,chromatin}.py 중 필요한 것. 이것만.

읽을 곳(읽기만): aleph/components/nucleus/**, aleph/components/incumbent/compartments.py
(LaminaParams 가 거기 있고 전부 "I0-B2 PI GAP, TEST values" 로 표시돼 있다).

⚠ 알려진 PI-GAP (발명 금지, 라벨 달고 raise):
   r_eq_um 5.0 (provisional) — 단 laws/cell_geometry.py 는 0.68 x R_cell = 5.1 을 sourced 로 들고 있고
   A1 이 laws/ 쪽을 썼다. 같은 값을 써라.
   eps_rupture 0.50 · k_linc 1.0e2 — 소스가 "8 pN is a TENSION not a stiffness" 라고 경고한다.
   lamina meshwork 크기 — A1 이 핵막 subdiv 6 을 이 값 없이 "유추" 로 세웠다고 보고했다. 네가
   그 값을 찾으면 A1 의 유추가 유도가 된다. 찾지 못하면 그렇게 보고해라.

⚠ chromatin 을 어떻게 표현할지가 설계 결정이다 — 노드 개체군? 밀도장? 가정하지 말고 근거와 함께
   Lead 에게 보내라.

구성은 Warp CUDA 커널. claim 은 개체군당. 기본값 금지.
```

---

# 세션 G — `ACTIVE_LOAD_PATH`: 모터 head 개체군

```
ffn_cellsim, engine/main.

읽어라: CLAUDE.md → STATE.md → aleph/engine/contracts.py 의 ComponentRole →
aleph/engine/nmii_actuator.py (ACTIVE_LOAD_PATH 로 선언돼 있다, 읽기만) →
aleph/world/bond.py docstring 전체.

너의 범위: NMII 모터 head 개체군. incumbent 는 442 seed 를 갖고 있고 아레나는 0 이다.
lamellipodium/filopodium/SF 가 ACTIVE_LOAD_PATH 인데 능동성이 전혀 없는 이유가 이것이다.
너의 파일: aleph/world/build/nmii.py. 이것만.

읽을 곳(읽기만): aleph/components/motor/** — 특히 resting_setpoint.py, segment_motor.py,
minifilament_topology.py, hand.py. 여기에 head-resolved Stam-Hocky bipolar minifilament 가
이미 구현돼 있다.

⚠ 이건 개체군이면서 결합이다. 순서에 주의해라:
   1) minifilament 백본 + head 를 노드 개체군으로 세운다  ← 이번 범위
   2) head 가 actin 에 붙는 것은 BOND 이고, bond.py 가 kinetics 를 "the first family that has
      kinetics" 까지 일부러 미뤄놨다. 네가 그 첫 번째다 → PI 에게 올려라, 짓지 마라.

⚠ 알려진 PI-GAP: NMII_N_SIDE=10 (AFINES) vs 28-30 (Billington) 미결 —
   2026-08-20 감사가 Melli 2018 eLife 와 Tripathi 2021 JoVE 에서 "30 motors per filament end" 를
   찾아 claim_b 쪽을 지지한다. NMII_L_BB_UM=0.301 · NMII_HEAD_OFFSET_UM=0.200 도 PI GAP.
   F_stall_head 0.5 vs 2.0 pN 은 2026-08-20 source audit 결과 "어느 쪽도 측정값이 아님".
   전부 라벨 달고 raise. 발명 금지.

구성은 Warp CUDA 커널. claim 은 개체군당.
```

---

# 세션 H — 능동성: lamellipodium · filopodium · SF 를 실제로 능동으로

```
ffn_cellsim, engine/main.

읽어라: CLAUDE.md → aleph/engine/contracts.py 의 ComponentRole →
aleph/engine/protrusion.py 와 stress_fiber.py (둘 다 ACTIVE_LOAD_PATH, 읽기만) →
docs/v2_audit/PER_CELL_STRUCTURE_COUNTS_2026-08-20.md.

너의 범위: 이미 서 있는 세 ACTIVE_LOAD_PATH 개체군에 능동성을 준다. 지금은 기하만 있다.
너의 파일: aleph/world/active/{protrusion,contraction}.py (신규 디렉터리). 이것만 —
aleph/world/build/ 는 A1/A2 것이다.

능동성이란 무엇인가, 저장소가 이미 말하는 것:
   lamellipodium  KB-3.6 Brownian ratchet: v_p=delta[k_on0 c_G exp(-F delta/kBT)-k_off],
                  delta=2.7nm, F_stall/filament ~2-5 pN, v0~200 nm/s
                  ⚠ 이 줄의 원본에 있던 "Hill w~2" 는 취소한다. laws/motility_warp.py 의 ratchet 은
                    v = v0*exp(-f*delta/kT) 이고 Hill 지수 w 를 읽는 자리가 없다. 그리고 NMII
                    force-velocity 는 PI 2026-07-07 확정으로 LINEAR 이다 (laws/myosin_linear.py:7 —
                    "a non-muscle-myosin-IIA Hill hyperbola does NOT exist in the literature").
                    w 를 살리려면 새 법칙이고 laws/ 소관이다 — 바인딩 세션이 쓸 것이 아니다.
   filopodium     KB-3.8 dL/dt = v_p - v_retro - k_cap L
   SF             KB-DRAFT-7-08 prestressed actomyosin cable, 모터가 붙어야 수축
                  ⚠ A2 가 SF 에 모터 station 을 하나도 안 놓았다 — STATE.md (e) 5 의 N_stations: 0
                    "curved sarcomere ... a rigid bipolar minifilament cannot straddle them" 때문.
                    이건 MODEL 사유이고 "do not force it" 이다. 세션 G 와 조율해라.

⚠ 이 세션은 laws/ 에 없는 것을 쓰게 된다. 그러면 그건 새 법칙이고, laws/ 에 두어야 한다
   (world/ 는 법칙을 소유하지 않는다). Lead 와 상의해라.

⚠ lamellipodium 개수가 2.5-10배 과소다 — KB-3.7 이 verified 로 '~100 fil per µm leading edge'
   를 갖고 있고 폭 5-20 µm 이면 500-2,000 인데 현재 200 이다. leading edge 폭은 PI 결정 대기.

주장 금지: 어떤 속도도 힘도. 이번 범위는 능동 채널이 존재하고 부호가 맞는다까지다.
```

---

# 세션 D1~D28 — 커넥터 (초판과 동일, 주소만 정정)

⚠ 초판 프롬프트의 `d38004` 는 틀렸습니다. Lead 는 **`ffn-cellsim-c8`** 입니다.
⚠ **스캐폴드는 이미 있습니다**: `aleph/world/families/__init__.py` — `FamilySpec` ·
`ConnectorGapError` · `DUPLICATE_OF`. `b9de49cd` 에 착지.

한 세션이 파일 하나(`aleph/world/families/<커넥터>.py`)를 소유합니다. 커넥터 28개 목록은
초판 파일 말미에 있습니다. 프롬프트 본문은 초판 그대로 쓰되 Lead 주소만 `ffn-cellsim-c8` 로,
그리고 GPU 는 위 §"상주 서버" 규칙을 따르십시오.

---

# 열려 있는 PI 결정 (세션들이 여기서 막힙니다)

| # | 결정 | 막는 것 |
|---|---|---|
| 1 | **lamellipodium leading edge 폭** (5–20 µm) | 세션 H. `KB-3.7` verified 밀도로 개수 유도됨 |
| 2 | **Mogilner 2005 수집** 또는 filopodia 간격 축 선언 | `KB-3.8` 의 `r_filo` 가 검증 불가 |
| 3 | **MT 개수 축** 250–600, 그리고 상피 프록시를 선암에 써도 되는가 | 세션 A2 산출물 |
| 4 | **IF · SF 개수 축** — 밴드 근거 없음 | 세션 A2 · H |
| 5 | `k_xl` 축 | 처리량 |
| 6 | `nvidia-smi --accounting-mode=1` (root) | `exact_peak_gpu_bytes` |
| 7 | NMII `N_side` 10 vs 28–30 | 세션 G |

---

# 세션 D — 커넥터. ⚠ **지금은 14개만** 열 수 있습니다

28개 중 절반은 붙일 상대 개체군이 아직 없습니다. 실측:

```
양쪽 개체군이 지금 다 있는 것   14   ← 지금 열 수 있음
빠진 개체군을 기다리는 것       14   ← cytosol 9 · nmii 4 · FA 2
```

세션 **E(cytosol)** 가 착지하면 9개가, **G(모터)** 가 착지하면 4개가 더 풀립니다. FA 2개는 아직
아무도 맡지 않았습니다(`ENVIRONMENT` 계열, 범위 밖).

## 지금 여는 14개 — 각각 세션 하나, 파일 하나

| # | 커넥터 | 미리 아는 것 |
|---|---|---|
| D1 | `membrane_erm_cortex` | ⚠ **재작성 아니라 바인딩.** `compartments.py:573` 의 `erm_bell_force_kernel` 이 이미 단일 배열 · 전역 인덱스로 런치한다 (세션 B 확인). pair-배열 변종은 incumbent 의 사적 할당 때문에만 존재. ⚠ 개수는 `erm_density_per_um2` = `ABSENT_FROM_CONTRACT_GRAPH` → 막힘 예상 |
| D2 | `membrane_cortex_contact` | ⚠ D1 과 같은 물리인지 먼저 판정하라. 같으면 `DUPLICATE_OF` 에 넣고 가족을 만들지 마라 |
| D3 | `mt_nucleus_linc` | ⚠ D4·D5 와 함께 **LINC 3형제**. 하나의 물리일 가능성. `k_linc=1.0e2` 는 PI GAP 이고 소스가 *"8 pN is a TENSION not a stiffness"* 라고 경고한다 |
| D4 | `if_nucleus_linc` | 위와 동일 |
| D5 | `actin_cap_linc` | 위와 동일 |
| D6 | `nucleus_cortex_contact` | 핵–cortex 접촉. steric 인지 결합인지 판정 필요 |
| D7 | `mt_cortex_capture` | +TIP capture. `engine/microtubule_rig.py` 는 `STRUCTURAL_RIG` |
| D8 | `mt_sf_spectraplakin` | cytolinker. D9 와 같은 부류인지 확인 |
| D9 | `if_sf_plectin` | cytolinker |
| D10 | `sf_cortex_transient` | `engine/` 에 `sf_cortex_transient` 빌더가 이미 있다(읽기만). ⚠ SF 는 A2 가 모터 station 을 하나도 놓지 않았다 — `STATE.md` (e) 5 |
| D11 | `filopodium_cortex_root` | filopodium 뿌리. A2 의 filopodium 은 polar cap 위에 서 있다 |
| D12 | `filopodium_membrane_tip` | 팁–막 결합 |
| D13 | `lamellipodium_cortex_seam` | seam |
| D14 | `lamellipodium_membrane_contact` | ⚠ lamellipodium 개수가 2.5–10× 과소 (PI 결정 1번 대기) |

## D 세션 프롬프트 — `<<<커넥터>>>` 와 `<<<미리 아는 것>>>` 만 바꿔 쓰십시오

```
ffn_cellsim, 브랜치 engine/main 에서 작업한다.

읽어라: CLAUDE.md → STATE.md → aleph/world/bond.py 의 docstring 전체 →
aleph/world/families/__init__.py → docs/v2_audit/WORLD_PORT_PLAN_2026-08-20.md (§4 PHASE 3).

CANONICAL 은 aleph/world/. aleph/engine/ 와 aleph/components/ 는 동결 PORT SOURCE —
읽기만 하고 import 하지 마라 (test_layer_directions 가 막는다).

너의 커넥터: <<<커넥터>>>
너의 파일:   aleph/world/families/<<<커넥터>>>.py   ← 이 파일 하나만. 14개 세션이 동시에 돈다.

미리 아는 것: <<<위 표의 해당 줄>>>

읽을 곳(읽기만): aleph/engine/ 에서 이 커넥터가 어떻게 선언됐는지,
aleph/outputs/ac/connector_devicerun/native_record.json 의 해당 항목.

⚠ 30개 중 28개가 NOT_BOUND 다 — 런타임 객체가 아예 없다. stub 보다 약하다.
   즉 이건 "이식" 이 아니라 "최초 구현" 이다. 그렇게 대해라.

BondFamily 를 만들려면 BondCount 가 필수이고, 그건 네 질문에 답해야 한다:
   몇 개인가 · 무엇에 대해서인가 (면적? 부피? 필라멘트당?) · 어느 세포인가 · 누구의 권위로인가

답이 없으면 발명하지 마라. families/__init__.py 의 ConnectorGapError 를 raise 하고,
모듈에 SPEC: FamilySpec 을 노출해서 blocked_by 에 답 못한 질문을 이름으로 담아라.
Lead 가 그것들을 결정 큐로 모은다. **그 거부가 네 주 산출물이다.**

⚠ 실제 사례: incumbent 의 ERM 은 커넥터 1개가 개체군 N개를 대신했고, 그래서 아무도 개수를
   답할 필요가 없었고, 답이 알고 보니 icosphere 정점 수(subdiv 3 에서 642)였다. bond.py 가
   이걸 "a mesh number wearing a physiological label" 로 기록했다. 네 커넥터도 그런지 확인해라.

⚠ 중복 판정: 네 커넥터가 다른 커넥터와 **같은 물리**면 DUPLICATE_OF 에 넣고 가족을 만들지 마라.
   이름이 아니라 물리로 판단하고, 논거를 모듈 docstring 에 써라. bond.py: 커넥터는 전역 노드
   인덱스 2개만 저장하고 어느 개체군인지는 arena ID 범위에서 조회 시 유도된다 —
   그래서 "membrane 쪽 NMII" 와 "SF 쪽 NMII" 가 두 가족일 이유가 없다.

만들 것: BondFamily 빌더 + 테스트 + docstring 에 Sanity Gate + SPEC.
만들면 안 되는 것: kinetic law. bond.py 가 "the first family that has kinetics" 까지 일부러
   미뤄뒀다. 네 커넥터가 그 첫 번째라면 PI 에게 올려라, 짓지 마라.

⚠ GPU 직접 제출 금지. Lead 가 4090-1 홀드(job 72)와 상주 서버를 잡고 있다.
⚠ 보고: ListAgents → ffn-cellsim-c8 (Lead) 에게 SendMessage. 시작 / 막힘 / 완료.
⚠ ACTIVE_SESSIONS.md 금지 (14개가 동시에 쓰면 충돌). FFN_SESSION unset.

커밋: git commit -o aleph/world/families/<<<커넥터>>>.py
```

## 나중에 열 14개 — 무엇을 기다리는지

```
cytosol 대기 (세션 E) 9개:
  filopodium_cytosol_transfer · if_cytosol_transfer · lamellipodium_cytosol_transfer
  membrane_cytosol_boundary · mt_cytosol_transfer · nucleus_cytosol_boundary
  sf_cytosol_transfer · surface_porous_transfer · nmii_cytosol_transfer(+nmii)

nmii 대기 (세션 G) 4개:
  nmii_filopodium_motor · nmii_lamellipodium_motor · nmii_sf_motor · nmii_cytosol_transfer

FA 대기 (미배정) 2개:
  filopodium_nascent_fa · lamellipodium_nascent_fa
```

⚠ **cytosol 9개 중 6개가 `*_cytosol_transfer` 로 같은 이름꼴입니다** — 세션 B 가 이미 지적했듯
하나의 물리일 가능성이 높습니다. E 착지 후 **한 세션이 6개를 함께 판정**하는 편이 낫습니다.
같은 이유로 nmii 모터 3개도 함께.
