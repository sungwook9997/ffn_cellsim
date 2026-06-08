# KB Ops-Linkage Plan — 운영 실체(B) → Contract-Graph(A) 연결

**Date:** 2026-06-08
**Author:** Lead session
**Status:** DESIGN (PI 승인 대기)
**Tracking convention:** `docs/v2_audit/` (cf. `GPU_MAIN_PORT_*.md`)

---

## 0. 배경 — 왜 이 일을 하는가

2026-06-08 KB 시스템 평가에서 나온 1순위 결함:

> "Contract-Graph"인데 contract 쪽이 비어 있다. edges 986개 중 **876개(89%)가
> SourceEvidence↔KnowledgeClaim** 순수 문헌-인용. 정작 이 프로젝트의 정체성
> (paper-model wrapper 금지, KB→코드→런 기계적 연결)을 담아야 할 쪽은
> `model_contract=6, code_mapping=2, run_result=1`.

결과적으로 지식 시스템이 **둘로 쪼개져 평행선**을 달림:

- **(A) 형식 그래프** — Notion 8-DB → duckdb/Obsidian. 문헌 중심. TAG/RAG가 다루는 것.
- **(B) 운영 실체** — `outputs/h*/production/*.{json,md}` (86개), 유닛별 `REPORT.md` (8개),
  Dev-Logs day-logs, git 커밋 전체. 실제 실험 이력·결정·현재 상태는 전부 여기 있음.

온톨로지/그래프는 (A)만 덮는다. 본 계획은 **(B)를 (A)에 적재**해서
"어느 파라미터가 실패한 게이트를 먹이고, 그 게이트를 마지막으로 측정한 런/커밋이 무엇이며,
그 게이트의 contract를 구현한 모듈이 무엇인가"를 *문헌이 아니라 실제 런까지* 한 번에
질의 가능하게 만든다.

---

## 1. 설계 원칙 (철칙 준수)

1. **Notion = 유일 편집점 (SoT).** duckdb/vault는 절대 손대지 않는다. 도구는
   **Notion의 RunResult/CodeMapping DB에 행을 upsert**하고, 기존
   `tag_kb/refresh.sh` + `obsidian_rag_full/refresh.sh`가 read-layer를 재생성한다.
   edges는 `notion_to_duckdb`가 관계 프로퍼티에서 자동 생성하므로,
   Notion 행에 `gate`/`implements_contract` 관계만 채우면 edge는 공짜로 생긴다.
2. **No hallucination.** 결정적(deterministic) 필드는 LLM을 거치지 않는다.
   게이트/컨트랙트 링크는 **아티팩트 본문에서 regex로 추출한 id만** 허용 —
   LLM은 순위만 매기고 원문에 없는 id를 새로 만들 수 없다. 미매칭은 빈 관계로
   두고 manifest에 플래그.
3. **Human-gated (no auto-register).** 기존 `SE_REGISTRATION_CANDIDATES_*` 패턴을
   그대로 따른다. 기본은 dry-run → `OPS_HARVEST_CANDIDATES_<date>.md` manifest 산출
   → PI 검토 → `--apply`. (cf. memory: references corpus "no auto-register".)
4. **Idempotent upsert.** 안정 키로 중복 방지: RunResult=`run_id`, CodeMapping=`path`.
   적용 전 Notion DB에서 해당 키 조회 → 있으면 update, 없으면 create.

---

## 2. 데이터 소스 ID (확인됨)

`obsidian_rag_full/notion_to_obsidian.py::DATA_SOURCES`:

| DB | Notion data source id |
|---|---|
| RunResult | `aa1a80911f6649edae956dab2bf96921` |
| CodeMapping | `e6df88ae27334a2280596941803eab2c` |
| ValidationGate | `277e1d5165e5403d959d5bf15ca4a317` |
| ModelContract | `fc987bb7e4d64212928608b889177929` |

**RunResult 필드:** `outcome, run_id, result_snapshot, notes, commit, artifact_path, date,
gate`(→ValidationGate 관계), `title`
**CodeMapping 필드:** `cm_id, path, implements_contract`(→ModelContract 관계), `status,
gate, notes, title`

---

## 3. 파이프라인 — `tag_kb/harvest_ops.py` (신규)

세 소스, 결정적 추출 + 경계된 LLM, Notion upsert.

### 3.1 Source 1 — production 아티팩트 (`outputs/**/production/*.{json,md}`, layer2, h1 root)

→ **RunResult** 한 행씩.

- **결정적 필드 (LLM 無):**
  - `artifact_path` ← 파일 경로 (repo-relative)
  - `commit` ← 그 파일을 마지막으로 건드린 커밋 (`git log -1 --format=%h -- <path>`)
  - `date` ← 해당 커밋 author-date (mtime fallback)
  - `run_id` ← 경로에서 안정 슬러그 (예: `RUN-h7-gate-a-native-prod`)
  - JSON이면 키-값에서 핵심 메트릭 후보 추출 (gamma_total, steps/s, outcome 등)
- **outcome 분류** ∈ {`smoke`, `production`, `feasibility`, `diagnosis`, `refute`,
  `confirm`}: 파일명/JSON 키 휴리스틱 우선, 모호하면 LLM 확인.
- **gate 링크:** 본문/파일명/JSON에서 `KU-x.y`, `Gate-A/B`, `VG-*` regex 추출 →
  ValidationGate 행과 id 매칭. **regex 매칭분만** 관계로 설정.
- **result_snapshot:** 유일한 자유텍스트 LLM 산출 (아티팩트 1개만 컨텍스트에). 1–2문장.
  `notes`에 `machine-generated 2026-06-08` 태그.

### 3.2 Source 2 — 유닛별 `REPORT.md` (h1,h2,h3,h4,h5,layer2,h7,h1_baoab_freeze,h_0_2)

- 구조화된 헤딩/표 파싱 (Sanity-Gate PASS/FAIL 표, acceptance-band 표).
- 유닛 closeout RunResult 1행 (outcome=closeout summary, commit=REPORT 최종 커밋).
- REPORT가 명시하는 스크립트/모듈 → **CodeMapping 후보**로 수집 (3.3로 넘김).
- ⚠️ ValidationGate.status는 건드리지 않는다 (게이트는 contract, read-only). REPORT의
  PASS/FAIL는 cross-check용 manifest 노트로만.

### 3.3 Source 3 — 모듈/스크립트 → **CodeMapping**

- `ffn_sim/{cortex,cell,ecm,cortex,bridge,junction,integrator,native}/*.py` +
  `scripts/*.py` 중 REPORT/configs/RunResult가 참조하는 것.
- `cm_id` ← 안정 슬러그, `path` ← repo-relative, `status` ∈ {draft, implemented}
  (테스트 존재/REPORT 언급으로 판정).
- `implements_contract` ← 모듈 docstring의 KU/contract id를 ModelContract 행과 매칭
  (regex). 미매칭은 빈 관계.

### 3.4 git + Dev-Logs

- Dev-Log day-logs는 이미 Notion SSOT 내러티브 → **중복 생성 안 함**.
  RunResult.commit / artifact_path 정합성 backfill 용으로만 사용.

---

## 4. 연결(edges)은 공짜

`RunResult.gate` (→ValidationGate), `CodeMapping.implements_contract` (→ModelContract)
관계를 채우면 `notion_to_duckdb`가 `run_result↔validation_gate`,
`code_mapping↔model_contract` edge를 자동 생성. 신규 다중-홉 질의가 열림:

> "실패한 Gate-A에 대해 — 마지막으로 측정한 run/commit은? 그 contract를 구현한 모듈은?
> 그 모듈이 쓰는 파라미터의 문헌 출처는?"

run → gate → contract → parameter → source_evidence 까지 한 경로로 traverse.

---

## 5. 유지보수 통합 (권고 #2·#3 일부 해소)

- `harvest_ops.py`는 클린 파이프라인에 둔다. 날짜박이 일회용 스크립트
  (`link_se_claims_20260604.py` 등)는 별도 `tag_kb/_migrations/`로 격리.
- `refresh.sh`에 `--check` 모드: "지난 harvest 이후 신규 아티팩트 N개" 드리프트만 보고
  (자동 적용 X, PI-gate 유지).
- harvest 자체는 dry-run + PI-gate 기본 (citation-integrity 규율과 동일 결).

---

## 6. 단계 (phasing)

- **P0** — `harvest_ops.py` 스켈레톤 + 결정적 추출기 + dry-run manifest. 대상:
  h7 + layer2 (최근·고가치 런)부터. Notion 쓰기는 아직 안 함.
- **P1** — LLM snapshot + gate-regex linker. 86 아티팩트 + 8 REPORT 백필 manifest.
  PI 검토.
- **P2** — `--apply`로 Notion upsert. CodeMapping harvest. refresh 2종 실행.
  edges 검증.
- **P3** — `refresh.sh --check` 드리프트 감지. CLAUDE.md KB 섹션에 절차 문서화.

---

## 7. Acceptance (검증 기준)

- `run_result` 행: 1 → ~30+ (구분되는 production 런마다, dedup).
- `code_mapping` 행: 2 → ~15 (핵심 모듈).
- edges: `run_result→validation_gate` ≥ 15, `code_mapping→model_contract` ≥ 10.
- `tag_query.py`가 기존에 못 풀던 다중-홉을 디스크 아티팩트 인용과 함께 답함
  (§4 예시 질의).
- **할루시네이션 게이트 링크 0건** (manifest 감사로 확인).

---

## 9. P0 RESULT — 2026-06-08 (DONE)

`tag_kb/harvest_ops.py` 구현·실행 완료 (결정적+regex, LLM/Notion 미사용).
산출: `tag_kb/OPS_HARVEST_CANDIDATES_2026-06-08.md` (dry-run manifest).

- **RunResult 후보 53** (production 아티팩트) + **closeout 후보 8** (REPORT.md).
  현 그래프 run_result=1 → 적재 시 ~50×.
- gate-linked 6/53 (KU-3.5 게이트 보유 아티팩트). closeout 중 h3/h5/layer2 →
  `VG-H3-KU35-cortex-tension` + `MC-H3-composite-tension` 매핑.
- **⭐ 핵심 발견 (PI 결정 필요): H.7 Gate-A(6건)/Gate-B(5건)가 ValidationGate에
  미등록.** 최근 작업 대부분이 *contract-graph에 존재하지 않는 게이트*를 참조 중.
  harvester는 자동생성하지 않고 manifest에 플래그(철칙 §2 준수). → PI가
  H.7 Gate-A/B를 정식 ValidationGate 행으로 등록할지 결정해야 적재 시 edge가 연결됨.
- 버그 2건 잡음(개발 중): (a) bare KU-x.y(=KnowledgeClaim)를 "누락 게이트"로
  오플래그 → KU는 positive 매칭에만 사용. (b) `gate_a_native`의 `_`가 word-char라
  `\b` 미작동 → lookaround `(?<![a-z])gate[-_ ]?([ab])(?![a-z])`로 교체.

**다음(P1, PI 체크포인트 후):** LLM result_snapshot + Notion upsert(`--apply`,
idempotent run_id/path 키) + H.7 게이트 등록 결정 반영.

## 10. P1 RESULT — 2026-06-08 (DONE, PI-approved)

PI 결정(AskUserQuestion): **신규 VG-H7-gate-a/b 2행 등록** (KU-3.5에 묶지 않음).

- **게이트 등록** (PI 승인 gate-contract 변경): `VG-H7-gate-a` (Validation,
  status=failing/REFUTE) + `VG-H7-gate-b` (Validation, status=blocked), 둘 다
  Band=[0.35,0.65] mN/m, Tests Contract→`MC-H3-composite-tension`. H.7 결과에 충실히
  Notes 작성 (Gate-A: s_grip 22×↑/γ_soft 16.2×→3.06e-3=1/114 band, transmission-limited;
  Gate-B: relaxed-M-SHAKE 레버, 3채널 분리, γ_active=0). Notion pages
  379120daec5d814a94a7c80f05482933 / 379120daec5d81bb91a1df65d3b3e12a.
- **P1 코드** (`harvest_ops.py` 확장): `--llm`(opt-in result_snapshot, 공유 LLM 백엔드)
  + `--apply`(멱등 RunResult upsert, RUN ID 키). Outcome은 PASS/FAIL/smoke로 매핑,
  verdict 불명이면 생략(가짜 PASS/FAIL 금지). Gate-A/B는 VG-H7-* 등록 후 자동 링크.
- **적재** (bulk Notion write): h7 16건 검증 후 full. **45 created + 16 updated = 61**
  (h7 재실행분 16건 중복 없이 UPDATE → 멱등성 검증). 스냅샷은 결정적(원시 메트릭/verdict,
  할루시네이션 0); LLM은 비용·anti-hallucination 이유로 생략.
- **refresh + 검증:** duckdb + Obsidian(766→840 노드) 재생성.
  - `run_result` 행 **1 → 62** (목표 ~30+ 초과달성)
  - `run→gate` edges **21** (목표 ≥15 초과). VG-H7-gate-a←6런, VG-H7-gate-b←5런.
  - `tag_query`: 예전 답 불가하던 다중-홉 ("Gate-A에 연결된 RunResult+commit+outcome")을
    디스크 아티팩트 인용과 함께 정확 반환 — **PI 약속한 가치 실증.**

## 11. 권고 #2/#3 — DONE (2026-06-08)

- **#2:** 일회용/날짜박이 스크립트 6개 → `tag_kb/_migrations/` (refresh.sh 무영향, README 추가).
  live 파이프라인 = notion_to_duckdb·references_ingest·supersession·tag_query·harvest_ops.
- **#3:** `verify_sources.py --check` (비파괴 citation-integrity 읽기) → `refresh.sh` 배선.
  즉시 발견: source_audit **12/329 부분/stale** + suspect 2건(Funk2021_eLife DOI_DEAD,
  Lindstrom2010 CHECK) → full `verify_sources.py` 재감사 필요(별도).

## 13. P2 RESULT — CodeMapping harvest (DONE, 2026-06-08)

`harvest_ops.py`에 `harvest_code`/`upsert_code` 추가: 메커니즘 모듈 docstring →
`code_mapping`. 결정적 추출 + docstring KU/MC regex → ModelContract 매칭.
**비파괴**: 기존 큐레이션 CM 행은 Path로 skip(클로버 방지).

- 대상 = 패키지 메커니즘 모듈 (cortex/cell/ecm/bridge/junction/integrator/native/common,
  `__init__`/test 제외) **46개**.
- status 휴리스틱: 테스트 스위트가 모듈을 참조하면 implemented, 아니면 draft.
- 적재: **44 created + 2 skipped** (myosin.py + cortex combo 큐레이션 행 보호).
- refresh 후: `code_mapping` **2 → 46**, `code→contract` edges **8** (총 edges 986→1042).
  MC-H3-composite-tension 구현 모듈 6개 정확 연결 (membrane_surface, nucleus, checkpoint,
  active_gel_seam, enclosed_volume, turnover).
- 매처 개선(정당한 correctness fix): 복합 KU 표기 `KU-3.5/3.1`의 "/3.1"을 파싱하도록
  `_ku_tokens` 보강 (nucleus/enclosed_volume의 KU-3.1 누락 해결, 4→6 link).

### ⭐ 발견 (Gate-A/B와 동형, PI 항목): ModelContract 레이어가 희소
`code→contract`가 8에 그치는 이유 = **ModelContract 행이 6개뿐이고 KU 토큰을 가진 건
MC-H3-composite-tension 하나**. FA(bridge/)·ECM·integrator·lamellipodium 등 대부분
메커니즘에 **contract 행 자체가 없음**. ValidationGate에 H.7 게이트가 없던 것과 같은 구조적
공백 — ModelContract도 PI-authored 행 보강이 필요(코드/런은 이제 적재됐으므로 contract만
채우면 전체 traverse가 열림).

## 12. 남은 일 (추가)

- ModelContract 레이어 보강 (PI-authored; FA/ECM/integrator/lamellipodium contract 행).
- syn(NL→SQL) + BM25 retrieval 안정화 (TAG "불안정"의 본질 — gate coverage와 별개).
- source_audit full 재감사 (#3가 stale 12/329로 surfacing).
- P3: `refresh.sh --check` 드리프트 감지 (신규 아티팩트 N개 since last harvest).

## 8. 비-목표 (scope out)

- ValidationGate/ModelContract/Parameter 행 신규 생성·수정 (read-only; contract는 사람이).
- Dev-Log 내러티브 중복 적재.
- 정식 OWL/RDF 온톨로지 엔진 — 현 규모(8 타입/~16 관계)엔 오버엔지니어링.
  DuckDB edges + supersession 으로 충분 (별도 평가 결론).
