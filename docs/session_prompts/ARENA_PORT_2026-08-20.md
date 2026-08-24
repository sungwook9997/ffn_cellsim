# 병렬 세션 프롬프트 — 아레나 포팅 2026-08-20

각 블록을 그대로 복사해 새 Claude Code 세션에 붙이면 됩니다.
**세션 0(플랫폼)은 Lead 세션 `d380041d`이 짓고 있습니다.** 나머지는 바로 열 수 있습니다.

---

## 모든 세션 공통 — 각 프롬프트에 이미 들어 있음 (참고용 요약)

- 브랜치 **`engine/main`**. `git commit`은 **반드시 pathspec** (`git commit -o <경로>`).
- **CANONICAL은 `aleph/world/`**. `aleph/engine/`·`aleph/components/`는 **PORT SOURCE, 동결** — 읽기만.
- `aleph/laws/`는 아레나가 **바인딩**하는 커널 라이브러리.
- `world/`는 `engine`/`components`/`dcm`/`scripts`/`validation`/`archive` **import 금지**.
- **기본값 금지 / 네이티브만.** PI-GAP은 표면화하지 발명하지 않습니다.

### ⚠ GPU는 Lead 가 독점합니다 — 직접 제출 금지

PI 승인은 **4090-1 한 장, 12시간 홀드 하나**입니다. 세션이 각자 `gpu-submit` 하면 서로
홀드를 뺏거나 예약 안 한 카드로 흘러갑니다. **Lead 세션이 홀드를 잡고 상주 서버를 띄웁니다.**
세션은 그 서버를 호출하지 자기 job 을 제출하지 않습니다. GPU 가 필요하면 Lead 에게 요청.
3090은 다른 사람(ham) 것입니다.

### ⚠ 보고는 SendMessage 로 Lead 에게 직접

`ListAgents` 로 Lead 세션을 찾으십시오 — **`ffn-cellsim-c8`** 입니다 (ref `bc0d8f`). ⚠ 이 파일의 초판은 `d38004` 라고 적었는데, 그건 Lead 의 세션 UUID 앞자리이고 ListAgents 의 ref 와는 다른 값입니다 — A1 이 Lead 를 못 찾아 두 세션에 중복 발신하고서야 드러났습니다.
`SendMessage({to: "<그 이름>", message: "..."})` 로 아래 시점마다 보고:

1. **시작** — 무슨 파일을 소유할지
2. **마일스톤** — 무엇이 섰고 무엇이 측정됐는지 (숫자 포함)
3. **막힘** — PI-GAP 이나 결정 대기. **이게 제일 중요합니다.** Lead 가 결정 큐로 모읍니다
4. **완료** — 커밋 해시

### ⚠ `docs/ACTIVE_SESSIONS.md` 는 Lead 만 씁니다

여러 세션이 동시에 append 하면 충돌합니다. 세션 행은 **Lead 가 대신 적습니다** — 시작할 때
SendMessage 로 알리기만 하십시오.

### ⚠ `FFN_SESSION` 은 unset 으로 두십시오

설정하면 `check_ownership.py` 가 켜지는데, 지금 레인 배치가 병렬 세션용이 아니라 커밋이
거부됩니다. 보호는 **파일 분리**로 합니다 — 한 세션이 자기 파일만 만집니다.

---

# 세션 A1 — PHASE 1a: cortex · 막 · 핵막 (빌더가 있는 것)

```
ffn_cellsim, 브랜치 engine/main 에서 작업한다.

먼저 읽어라: CLAUDE.md → STATE.md → docs/v2_audit/WORLD_PORT_PLAN_2026-08-20.md
(§4 PHASE 1) → aleph/world/arena.py 와 strand.py 와 surface.py 의 docstring 전체.

CANONICAL 은 aleph/world/ (아레나) 다. aleph/engine/ 와 aleph/components/ 는 PORT SOURCE
로 동결됐으니 읽기만 하고 절대 import 하지 마라 — test_layer_directions.py 가 막는다.

너의 범위: cortex actin · 형질막 · 핵막. 이미 incumbent 에 빌더가 있는 것들이다.
너의 파일: aleph/world/build/cortex.py, membrane.py, envelope.py — 이 셋만 만져라.
세션 A2 가 MT/IF/filopodium/lamellipodium/SF 를 동시에 짓고 있다.

구성은 반드시 Warp CUDA 커널로 한다 (PI 결정 2026-08-20). 호스트 NumPy 로 만들어 업로드
하는 방식은 쓰지 마라. claim 은 개체군당 한 벌이다 — 지금 build_strand 가 가닥당 3개를
청구해서 네이티브 cortex 가 212,058 claim 이 되는데, 개체군당으로 바꾸는 것이 이 페이즈의
설계 결정이다. arena.py 원문: "a population is a NAMED, CONTIGUOUS, HALF-OPEN ID RANGE".

목표 (밀도 x 기하에서 유도할 것 — 이 숫자를 타이핑하지 마라):
  cortex actin  seg_um 을 sourced 값 0.05 µm 로. 현재 0.5 는 감사가 "5-10x TOO COARSE,
                provenance CONVENIENCE" 로 기록했다 → 약 4.9M 노드
  형질막         icosphere subdiv 7 → 163,842 정점
  핵막           
밀도가 다른 값을 주면 그 값이 답이다. 4.9M 은 추정치다.

증명: 세 개체군이 서고 arena.assert_partitioned() 통과.
측정: exact_peak_gpu_bytes (이 프로젝트에서 아무도 잰 적 없다), 빌드 벽시계, 노드당 바이트.
주장 금지: 물리 아무것도. 스텝 처리량 (아직 아무것도 안 구른다).

⚠ GPU 직접 제출 금지. Lead 가 홀드를 잡고 상주 서버를 띄운다. 그 엔드포인트를 써라.
⚠ 보고: ListAgents 에서 ffn-cellsim-c8 (Lead) 을 찾아 SendMessage 로 시작/마일스톤/
   막힘/완료를 보내라. 막힘(PI-GAP)이 제일 중요하다.
⚠ ACTIVE_SESSIONS.md 는 건드리지 마라. Lead 가 적는다.
⚠ FFN_SESSION 은 unset 으로 둬라.

커밋: git commit -o aleph/world/build/<너의 파일>
```

---

# 세션 A2 — PHASE 1b: MT · IF · filopodium · lamellipodium · SF (노드 0개인 것)

```
ffn_cellsim, 브랜치 engine/main 에서 작업한다.

먼저 읽어라: CLAUDE.md → STATE.md → docs/v2_audit/WORLD_PORT_PLAN_2026-08-20.md
(§1, §4 PHASE 1) → aleph/world/arena.py 와 strand.py 의 docstring 전체.

CANONICAL 은 aleph/world/ 다. aleph/engine/ 와 aleph/components/ 는 동결 — 읽기만,
import 금지 (test_layer_directions.py 가 막는다).

너의 범위: 오늘 노드가 0개인 컴파트먼트 — microtubule, intermediate filament,
filopodium, lamellipodium, stress fiber. 목표 10M 의 절반이 여기서 나온다.
너의 파일: aleph/world/build/{microtubule,intermediate_filament,filopodium,
lamellipodium,stress_fiber}.py — 이것들만. 세션 A1 이 cortex/막/핵막을 동시에 짓는다.

참고할 기존 빌더 (읽기만, import 금지):
  aleph/components/solid/microtubule.py           build_microtubule_compartment
  aleph/components/solid/intermediate_filament.py build_if_compartment
  filopodium / lamellipodium / SF 는 흩어져 있으니 grep 으로 찾아라.

구성은 반드시 Warp CUDA 커널. claim 은 개체군당 한 벌.

⚠ 여기가 PI-GAP 이 제일 많이 터질 곳이다. 밀도가 없으면 빌드가 raise 하게 만들고
   무엇이 막는지 정확히 적어 Lead 에게 SendMessage 로 보내라. 발명하지 마라.
   이미 알려진 것: NMII N_side 10 vs 28-30 미결, k_linc 는 "8 pN is a TENSION not a
   stiffness", 자세한 목록은
   docs/v2_audit/PARALLEL_AUDIT_PARAMS_DETECTORS_2026-08-20.md §3.

증명: 각 개체군이 밀도 x 기하에서 유도된 개수로 서고 assert_partitioned 통과.
주장 금지: 물리 아무것도.

⚠ GPU 직접 제출 금지 — Lead 의 상주 서버를 써라.
⚠ 보고: ListAgents → ffn-cellsim-c8 (Lead) 에게 SendMessage.
⚠ ACTIVE_SESSIONS.md 금지. FFN_SESSION unset.

커밋: git commit -o aleph/world/build/<너의 파일>
```

---

# 세션 B — PHASE 2: 법칙 바인딩 + 파리티

```
ffn_cellsim, 브랜치 engine/main 에서 작업한다.

먼저 읽어라: CLAUDE.md → STATE.md → docs/v2_audit/WORLD_PORT_PLAN_2026-08-20.md
(§4 PHASE 2) → aleph/scripts/arena_force_parity.py 전체.

너의 범위: PHASE 2. 세션 A 가 세운 개체군에 laws/ 커널을 붙이고, 각각을 incumbent 와
force-parity 로 검증한다. 세션 A 가 개체군을 하나 끝낼 때마다 그 개체군부터 시작해라 —
PHASE 1 전체를 기다리지 마라.

패턴은 이미 증명돼 있다: 78942fc4 (arena_force_parity.py) 가 laws.network_warp.
link_spring_kernel 과 laws.forces_warp.cytosim_bending_kernel 을 아레나 주소로 그대로
돌려 상대잔차 1.108e-16, 청구 범위 밖 정확히 0.0 pN 을 냈다. 같은 A/B 를 법칙마다 반복해라.
기준은 Higham 합산 한계 (n_terms x eps64) 이고 실행 전에 선언해야 한다.

⚠ 전부 바인딩이 아니다:
  laws/ 에 있음 (바인딩)  : link_spring, cytosim_bending, membrane area, helfrich_bending
  components/ 에 있음 (재작성): turgor traction (incumbent/membrane_pressure.py:48 — 유일한
                              디바이스 구현), ERM 힘 (incumbent/erm_tether.py)
  ∂V/∂x 는 세 벌로 흩어져 있다: components/nucleus/geometry.py:115 (호스트 NumPy),
  incumbent/membrane_pressure.py:48 (디바이스), resting_balance_oracle.py:94 (오라클 사본).
  이건 옮기는 게 아니라 laws/ 로 접는 작업이다.

측정해야 할 것: 법칙마다의 한계 스텝 비용. 지금 처리량 추정치는 외삽 하나뿐인데
(WORLD_SCALE_AND_THROUGHPUT_2026-08-20.md), 이 측정이 그걸 법칙별 예산으로 바꾼다.

주장하면 안 되는 것: 잔차가 내려간다는 것. 그건 PHASE 4 다.

⚠ GPU 직접 제출 금지 — Lead 가 4090-1 홀드를 독점하고 상주 서버를 띄운다.
⚠ 보고: ListAgents 에서 ffn-cellsim-c8 (Lead) 을 찾아 SendMessage 로
   시작/마일스톤/막힘/완료를 보내라. 막힘(PI-GAP)이 제일 중요하다.
⚠ ACTIVE_SESSIONS.md 는 Lead 만 쓴다. FFN_SESSION 은 unset.

커밋은 pathspec: git commit -o <너의 경로>
```

---

# 세션 C — 검출기: 배선과 상세화

```
ffn_cellsim, 브랜치 engine/main 에서 작업한다.

먼저 읽어라: CLAUDE.md → STATE.md →
docs/v2_audit/PARALLEL_AUDIT_PARAMS_DETECTORS_2026-08-20.md 전체.

측정된 사실: aleph/observe/ 는 4,888줄인데 7개 모듈 중 5개가 호출처 0곳이다.
  stationarity 1,030줄 → 3곳    operator 493줄 → 3곳
  fluctuation 904줄 → 0곳       tether 682줄 → 0곳
  kinematics 728줄 → 0곳        traction 619줄 → 0곳
  manifest 432줄 → 0곳

너의 범위: 우선순위대로 배선한다.
  1) tether — AFM 막-테더 forward operator. PI 의 처리량 목표 전체(6,514x)가 3분
     인덴테이션 프로토콜에서 나왔는데 그 계측기를 읽을 operator 에 호출처가 0개다.
     이게 가장 급하다.
  2) manifest — "A number without its protocol is not a measurement, it is a rumour."
     34개 필드, 기본값 없음, NOT_APPLICABLE(주장) 대 NOT_RECORDED(자백) 구분.
     아무것도 이걸 강제하지 않고 있다.
  3) 나머지 셋 (fluctuation, kinematics, traction) 은 각각 무엇을 읽어야 하는지
     분석만 하고 배선은 PI 판단으로 남겨라.

각 검출기마다 답해야 할 것: 무엇을 입력으로 받는가 (아레나의 어느 개체군/배열),
무엇을 내놓는가, 그 출력이 실험과 비교 가능한가, 그리고 지금 세팅이 맞는가.

⚠ stationarity 는 두 벌이 있고 24.8% 불일치한다
(tests/architecture/test_stationarity_modules_agree.py 가 고정해뒀다). 통합은
구현 하나로 접어야 하는 별도 작업이니 여기서 손대지 말고 상태만 확인해라.

⚠ GPU 직접 제출 금지 — Lead 가 4090-1 홀드를 독점하고 상주 서버를 띄운다.
⚠ 보고: ListAgents 에서 ffn-cellsim-c8 (Lead) 을 찾아 SendMessage 로
   시작/마일스톤/막힘/완료를 보내라. 막힘(PI-GAP)이 제일 중요하다.
⚠ ACTIVE_SESSIONS.md 는 Lead 만 쓴다. FFN_SESSION 은 unset.

커밋은 pathspec: git commit -o <너의 경로>
```

---

# 세션 D1~D28 — PHASE 3: 커넥터 하나씩

⚠ **Lead 세션이 스캐폴드(`aleph/world/families/`)를 커밋한 뒤에 여십시오.** 그 전에 열면
붙일 곳이 없습니다.

**한 세션이 파일 하나만 소유합니다** — `aleph/world/families/<커넥터이름>.py`. 그래서 28개가
동시에 돌아도 충돌하지 않습니다.

```
ffn_cellsim, 브랜치 engine/main 에서 작업한다.

먼저 읽어라: CLAUDE.md → STATE.md → docs/v2_audit/WORLD_PORT_PLAN_2026-08-20.md
(§4 PHASE 3) → aleph/world/bond.py 의 docstring 전체 → aleph/world/families/README.md.

너의 커넥터: <<<여기에 커넥터 이름 하나>>>
너의 파일:   aleph/world/families/<<<커넥터이름>>>.py   ← 이 파일 하나만 만지고,
             다른 파일은 읽기만 해라. 28개 세션이 동시에 돌고 있다.

너의 일: incumbent 의 이 커넥터 선언을 아레나의 BondFamily 로 만든다.
읽을 곳: aleph/engine/ 에서 그 커넥터가 어떻게 선언됐는지 (읽기만, import 금지),
그리고 aleph/outputs/ac/connector_devicerun/native_record.json 의 해당 항목.

⚠ 30개 중 28개는 NOT_BOUND 다 — 런타임 객체가 아예 없다. 즉 대부분은 "이식" 이 아니라
"최초 구현" 이다. 그렇게 대해라.

BondFamily 를 만들려면 BondCount 가 필수이고, 그건 네 질문에 답해야 한다:
   몇 개인가 · 무엇에 대해서인가 (면적? 부피? 필라멘트당?) · 어느 세포에서인가 ·
   누구의 권위로인가 (출처)

답이 없으면 발명하지 마라. 대신 그 슬롯이 PI-GAP 을 명시적으로 들고 raise 하게 만들고,
파일 상단 docstring 에 무엇이 막고 있는지 정확히 적어라. Lead 세션이 그것들을 결정
대기 목록으로 모은다.

⚠ 실제 사례: incumbent 의 ERM 은 커넥터 1개가 개체군 N개를 대신했고, 그래서 아무도 개수를
답할 필요가 없었고, 답이 알고 보니 icosphere 정점 수였다 (bond.py 가 이걸 "a mesh number
wearing a physiological label" 로 기록). 네 커넥터도 같은지 확인해라.

⚠ 중복 확인: NMII 모터 4개(cortex/sf/filopodium/lamellipodium)는 같은 물리일 가능성이
높고, cytosol_transfer 6개도 그렇다. 네 커넥터가 다른 커넥터와 같은 물리면 그렇게 적고
가족을 새로 만들지 마라. 이름이 아니라 물리로 판단해라.

만들 것: BondFamily 빌더 + 테스트 + docstring 에 Sanity Gate.
만들면 안 되는 것: kinetic law (bond.py 가 "the first family that has kinetics" 까지
일부러 미뤄뒀다 — 네 커넥터가 그 첫 번째라면 PI 에게 올려라).

⚠ GPU 직접 제출 금지 — Lead 가 4090-1 홀드를 독점한다.
⚠ 보고: ListAgents 에서 ffn-cellsim-c8 (Lead) 을 찾아 SendMessage 로 시작/막힘/완료를
   보내라. 네가 답할 수 없는 BondCount 질문이 Lead 의 결정 큐로 간다 — 그게 네 주 산출물이다.
⚠ ACTIVE_SESSIONS.md 는 건드리지 마라 (28개가 동시에 쓰면 충돌). FFN_SESSION unset.

커밋은 반드시: git commit -o aleph/world/families/<<<커넥터이름>>>.py
```

## D 세션에 배정할 커넥터 28개 (NOT_BOUND, 세포 내부)

```
actin_cap_linc                 filopodium_cortex_root
filopodium_cytosol_transfer    filopodium_membrane_tip
filopodium_nascent_fa          if_cytosol_transfer
if_nucleus_linc                if_sf_plectin
lamellipodium_cortex_seam      lamellipodium_cytosol_transfer
lamellipodium_membrane_contact lamellipodium_nascent_fa
membrane_cortex_contact        membrane_cytosol_boundary
membrane_erm_cortex            mt_cortex_capture
mt_cytosol_transfer            mt_nucleus_linc
mt_sf_spectraplakin            nmii_cytosol_transfer
nmii_filopodium_motor          nmii_lamellipodium_motor
nmii_sf_motor                  nucleus_cortex_contact
nucleus_cytosol_boundary       sf_cortex_transient
sf_cytosol_transfer            surface_porous_transfer
```

**이미 DEVICE_RUN 인 2개** (`nmii_cortex_motor`, `dorsal_arc_crosslink`)는 세션 배정에서
제외 — 이식할 실물이 있으므로 Lead 가 직접 처리합니다.

**제외한 외부 8개** (ECM/기질/매질) — 이번 범위 아님:
`ecm_crosslink` `ecm_far_field_anchor` `integrin_collagen_clutch` `membrane_ecm_contact`
`membrane_medium_traction` `fa_actin_anchor` `filopodium_nascent_clutch`
`lamellipodium_nascent_clutch`
