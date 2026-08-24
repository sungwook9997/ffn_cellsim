# Project Aleph 파일 구조 지도

작성: 2026-08-09 KST, 세션 `5ccfd28d`. `ffn_cellsim` → **Project Aleph** 병합·개명 시점의 목표 구조.

이 문서는 **파일 → 역할** 매핑 인덱스이고, 동시에 **마이그레이션 진행표**다. `상태` 열이 비어 있지
않은 한 이 문서는 목표를 서술하는 것이지 디스크를 서술하는 것이 아니다.

> **이 문서가 통째로 새로 쓰인 이유.** 이전 판은 2026-05-20 본문 위에 날짜 배너 네 개가 쌓인 형태였고,
> `STATE.md` (e) 4가 "존재하지 않는 디렉토리를 가리킨다 — 기록할 것인가 수리할 것인가"로 열어둔
> 항목이었다. 병합이 그 수리다. 이전 판은 git 이력에 있다.
>
> 결정의 근거와 인용은 [`docs/decisions/PROPOSAL-the-merge-and-the-rename.md`](docs/decisions/PROPOSAL-the-merge-and-the-rename.md).

---

## 0. 한 장 요약

`ffn_cellsim`의 **엔진**과 Project Aleph의 **규율**을 합친다. 엔진이 이쪽인 이유는 측정된 것이고
(`ff → ac` 순환 참조 **0건**, 네이티브 551,434 노드 / 1.65M DOF, hot-loop D2H 복사 0), 규율이 저쪽인
이유도 측정된 것이다(타입으로 강제되는 소유권, 코드로 된 증거 사다리, 원장 135건).

```
목표 트리 ≈ 96k 줄        아카이브 ≈ 134k 줄
```

**병합이 아니라 정리에 가깝다.** 버리는 쪽이 더 크다.

## 1. 트리 — 2026-08-09 착지 상태

**✅ = 디스크에 있음** (개명 + 층 재구조화 완료, 전체 스위트 green).
나머지는 목표이고, 그 행은 아직 서술이 아니다.

⚠ **2026-08-22 — `engine/`·`components/`는 ARCHIVED다.** PI 재정 *"아레나가 엔진이고, 다른 거는 모두 이제
아카이브화된 상황"* — 두 트리는 **읽어서 포팅해 나오는 대상이지 작업 대상이 아니고**, `ownership.yaml`의
`archive-readonly` 레인이 같은 커밋에서 그것을 강제한다. **줄 수와 파일 수는 그대로 둔다** — 아카이브는
삭제가 아니고, 이 표는 디스크를 서술한다. 근거와 대가(tier-(a) 행 대부분이 history가 된다)는
[`STATE.md`](STATE.md) (a)와 [`aleph/docs/STATE_LOG.md`](aleph/docs/STATE_LOG.md) 2026-08-22 항목.

| 경로 | 파일 | 줄 수 | 역할 | 상태 |
|---|---:|---:|---|---|
| **`aleph/world/`** | **47** | **18,779** | **아레나 — 고정 용량 노드 할당, 연속 ID 범위 claim. primitive 8종(STRAND/BOND/FACE/GRID_CELL 포함)** | ✅ **CANONICAL 2026-08-20** |
| `aleph/units/` | 3 | 244 | 인용 필수 acceptance band (`ALEPH-PORT-4001`) | ✅ |
| `aleph/laws/` | 49 | 12,510 | 공유 물리 법칙 + Warp 커널. **`world/`가 바인딩하는 라이브러리** | ✅ |
| `aleph/engine/` | 64 | 32,858 | component/connector 조립, 트랜잭션, dispatch | ✅ **ARCHIVED 2026-08-22**(archive-readonly) |
| `aleph/components/` | 119 | 36,112 | 부품·커넥터 (`incumbent`·`motor`·`ecm`·`fluid`·`solid`·`nucleus`·`weave`·`emergence`) | ✅ **ARCHIVED 2026-08-22**(archive-readonly) |
| `aleph/viz/` | 29 | 11,506 | 렌더·figure. **`cell_app.py` = 풀 개체군 네이티브 뷰어** | ✅ |
| `aleph/validation/` | 30 | 5,618 | 해석적 오라클 + `cytosim_parity`(외부) | ✅ |
| `aleph/dcm/` | 53 | 16,480 | DCM 솔버. **아카이브 불가** — `laws/relax.py`가 씀 | ✅ |
| `aleph/virtual_cell/` | 27 | 30,024 | 확률층 — `CellStatePosterior`·`ObservationOperator` | ✅ **동결(archive-readonly), 배선 0** |
| `aleph/archive/` | 206 | 33,882 | `lane_d` · 미사용 `common` 4개 · 선별에서 빠진 스크립트 | ✅ |
| ~~`aleph/common/`~~ | — | — | **해체 완료.** 3개→`laws/`, `sim_realtime`→`dcm/`, 나머지 4개→`archive/` | ✅ |
| `aleph/scripts/` | 96 | 29,586 | 드라이버 + 게이트 | ✅ |
| `aleph/tests/` | 293 | 60,711 | | ✅ |
| `docs/decisions/` | | | `PROPOSAL-*` · `INHERITED-*` | ✅ |
| `docs/ACTIVE_SESSIONS.md` | | | 동시 세션 채널 | ✅ |
| `ports/` | | | 원장 + `ALEPH-PORT-4001` | ✅ |
| `aleph/observe/` | 8 | 4,986 | 관측 연산자 + `stationarity` | ✅ **⬜→✅ 2026-08-21 확인** |
| `aleph/infer/` | 7 | 2,762 | fisher · identifiability · posterior · refusal | ✅ **⬜→✅** |
| `aleph/campaign/` | 4 | 497 | 사전선언 스윕, 좌표 주소 RNG | ✅ **⬜→✅** |
| `aleph/evidence/` | 4 | 1,233 | claim · retraction · digest | ✅ **⬜→✅** |
| `aleph/artifacts/` | 4 | 1,599 | | ✅ **⬜→✅** |
| `aleph/outer/` | 250 | 57,544 | 외부 연구에서 학습하는 층 | ✅ **⬜→✅** |
| `aleph/inner/` | 1 | 37 | 자기 accepted run 에서 학습하는 층 | ⚠ **디렉터리만 — 37줄** |
| `aleph/state/` | — | — | owner-local state, topology sector | ⬜ **미생성.** Warp 재작성 필요 |
| `aleph/kb/` | — | — | 현재 `outputs/tag_kb/`. 살아있는 설정이 런 기록 사이에 있다 | ⬜ **미생성.** 승격 필요 |

⚠ **2026-08-21 재측정.** 이 표는 2026-08-09 에 작성됐고 그때는 **목표 진행표**였다. 12일 뒤 다시
디스크에 대고 세어보니 `⬜` 여섯 줄이 이미 `✅` 였고 — `observe`·`infer`·`campaign`·`evidence`·
`artifacts`·`outer` — **정본 트리 `aleph/world/` 는 이 지도에 아예 없었다.** 진행표는 갱신되지 않으면
진행을 서술하지 않는다.

남은 `⬜` 둘은 stale 이 아니라 **미생성**이다. `aleph/inner/` 는 37줄짜리 디렉터리만 있어 세 번째
상태이고, 그렇게 표시한다 — 있음도 없음도 아닌 것을 ✅ 로 적으면 그게 오늘 감사 메모에서 잡아낸
"쓰인 적 없는 계획 파일"과 같은 종류가 된다.

## 2. 층 규칙 — 테스트로 강제하고 관례로 두지 않는다

```
units → state → laws → engine → components → scripts
                 ↑
                 laws 는 engine 과 components 를 절대 import 하지 않는다
```

- **`laws/ ↛ engine/`, `laws/ ↛ components/`.** ffn이 `ff → ac` 0건으로 지켜온 규칙이고, Aleph가
  `runtime → vertical` **53건**으로 잃은 규칙이다. Aleph의 순환은 파일 네 개가 만들었다 —
  `law_cases.py`(40건), `residency.py`, `resident_world.py`, `resident_cortex.py`. 전부 "법칙"이
  아니라 검증 하네스와 디바이스 측 조립이었고, 좋은 이유로 하나씩 들어왔다. **그래서 안 보인다.**
  CI가 세는 편이 싸다.
- **`aleph/** ↛ validation/**`.** 오라클은 런타임을 판정하고 그 역은 없다. 자기 오라클을 볼 수 있는
  런타임은 오라클에 맞춰지도록 만들 수 있는 런타임이다.
  ⚠ ffn 레이아웃에서는 이게 Aleph보다 어렵다 — Aleph의 `validation/`은 패키지 **바깥 형제**인데
  ffn의 것은 **패키지 안**이었다. 그래서 `validation/`을 최상위로 끌어올린다.
- **`observe/` · `infer/` · `inner/` · `outer/` · `campaign/` · `evidence/` · `artifacts/`는 accept
  경계 바깥**이므로 host 연산이 허용된다. 헌장이 금지하는 것은 경계 **안쪽**의 권위 있는 host
  상태뿐이고, 헌장 자신이 *"Python may configure, orchestrate and postprocess"* 라고 쓴다.
- **`state/`만 저촉된다.** Aleph의 `StateBlock`은 `np.ndarray`를 `ascontiguousarray` 후 read-only로
  붙들고 있고, 이건 권위 있는 host 상태다. 값어치는 numpy가 아니라 타입 규율 — 이름 필수, 단위 필수,
  `ANONYMOUS_BLOCK_NAMES` 거부, bare array 생성자 부재, 마스크를 다이제스트 대상으로, 스칼라 거부 —
  이고 전부 `wp.array` 위에서 성립한다. 다이제스트는 **accept 경계에서만** D2H.

## 3. 아카이브 — 삭제가 아니라 이동

> **이 표는 2026-08-09에 세 번 줄었다.** 초안은 훨씬 많이 버리려 했고, 측정할 때마다 트리가
> 감당 못 하는 양이라는 게 드러났다. 근거는 결정 기록 §1a·§1b·§1c.

**실제로 아카이브된 것은 `aleph/archive/` 2,975줄이 전부다.**

| 대상 | 줄 수 | 왜 |
|---|---:|---|
| `lane_d/` | 2,001 | 15회 import가 전부 자기 테스트, `STATE.md` (c) 12 |
| `checkpoint`·`filament_math`·`integrity`·`production_policy` | 974 | `common/` 잔여. **import하는 곳이 하나도 없다** |
| `scripts/` 194개 | 31,220 | 게이트도, tier-(a) 생산자도, 테스트가 부르는 것도, 그 1-hop도 아님 |
| Aleph `vertical/`·`scenarios/`·`runtime/` | ~78,000 | 안 가져옴 (ff/ac가 같은 대상을 더 잘 가짐) |
| Aleph `represent/` + `learn/` | 12,400 | 안 가져옴. 연산 봉인 선언층 → **문서로** |

**아카이브에서 빠진 것 — 전부 측정이 뒤집었다:**

| 되돌린 것 | 줄 수 | 반박한 측정 |
|---|---:|---|
| `ac/cell/` → `components/incumbent/` | 17,061 | **243회** import, 그중 `ac/engine`. strangle이 안 끝났다 |
| `virtual_cell/` → 잔류 | 29,097 | **26/26 모듈에 테스트**. 미배선이지 미검증이 아니다 |
| `ff/` 23개 모듈 → `laws/` | — | `implicit_ff`(708× 솔버) · `relax` · `cytosim_parity`(외부 오라클)가 포함돼 있었다 |
| `dcm/` → 잔류 | 16,480 | `laws/relax.py`가 implicit step을 호출 |
| `gamma_floor{,_dynamic,_sweep}` → `laws/` | — | **결함이 살아있음을 기록하는 가드 테스트**가 import한다. (c) 2가 막는 건 *숫자*지 모듈이 아니다 |
| `architecture_metrics` → `laws/` | — | `laws/ecm_library`가 호출하고 오라클 테스트가 그 호출을 검사한다 |

## 4. 아직 서술이 아니라 경고인 것

1. **`ff/ENGINE.md`가 자기 헤드라인을 전부 철회한 상태다.** 배너가 구조 결함 5개 중 **3개가 γ-floor
   프로덕션 경로에 살아 있다**고 적고, `n_xl`이 caller 인자인데 빌더가 `ff/architecture_spec.CORTEX`를
   읽지 않아 모든 호출부가 crosslink 밀도 **1.0**으로 돈다고 적는다. `STATE.md` (c) 2가 이미 막고
   있지만, `laws/`로 옮길 14개 모듈은 이 경로와 겹치는지 확인 후에 인용 가능하다.
2. **`crosslinker.py:139`가 이미 `aleph.laws`를 import한다.** 가드
   (`test_exactly_one_module_imports_the_provider`)가 세고 있고, 그 가드는 아직 **열린** `46143f30`
   행의 소유다. 병합이 import 표면을 바꾸면 **보고하고, 고치지 않는다.**
3. **`Project_Aleph` 트리의 강등은 보류다.** `b4fad06b`가 그 안에서 열려 있다. 살아 있는 세션이 든
   트리를 개명하는 것은 `ACTIVE_SESSIONS.md`가 2026-07-30 16:14에 기록한 사고 그 자체다.
4. **`ROADMAP.md`의 예산 숫자 중 외삽에 기댄 것은 재유도 대상이다.** 어느 기획 문서의 `steps^2.08`
   외삽이 GPU 실측에서 **0.97**이었고, `ALEPH-PORT-3665`의 IMEX가 이미 10× `dt`에서 전세포 **28.3×**를
   낸다.
