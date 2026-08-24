# TAG — Table-Augmented Generation layer for the ffn_cellsim KB

작성: 2026-06-02 KST. PI 비준 결정 2건 반영 (백엔드=로컬 DuckDB, 콘텐츠 깊이=테이블+PDF 풀텍스트 완전체).

> ⚠️ 용어 주의: 여기서 **TAG = Table-Augmented Generation** (Biswal et al. 2024). "태그(label)"가 **아님**.

---

## 0. Context — 왜 이걸 만드나

현재 RAG는 **Notion semantic search + Obsidian 그래프(semantic 링크)** 2종이다. "관련 노트 찾기"엔 좋지만, Contract-Graph 위의 **관계형·계산형 질문**(조인·집계·다중 홉 필터)은 못 푼다. 예:

- *"Direct-measurement 논문이 근거인데 현재 `failing` ValidationGate에 물린 Parameter는?"* — SE→KB→MC→Param→VG 다중 조인 + 상태 필터.
- *"`draft` ModelContract 중 evidence가 Model-derived만 있고 Direct measurement가 없는 것"* — 조인 + 집계 + 의미 필터.
- *"KB-3.x를 지지하는 SE 중 DOI 없는 것 몇 개"* — 집계.

결정적 관찰: **Notion Contract-Graph는 이미 정규화된 관계형 스키마**(데이터 타입별 9개 atomic DB, 양방향 relation)다 — 즉 TAG가 필요로 하는 *테이블 레이어*가 이미 존재한다. 우리는 그 위에 query 엔진만 얹으면 된다.

**TAG (Biswal, Patel, Jha, Kamsetty, Liu, Gonzalez, Guestrin, Zaharia. "Text2SQL is Not Enough: Unifying AI and Databases with TAG." arXiv:2408.14717, 2024)** — 3단계:
1. **syn (합성)**: 자연어 질문 R → 실행 가능한 DB 쿼리 Q
2. **exec (실행)**: Q를 DB에 실행 → 관련 데이터 T
3. **gen (생성)**: (R, T) → 자연어 답변 A (필요 시 본문 semantic 추론 결합)

TAG는 Text2SQL(gen=항등)과 RAG(syn=임베딩 lookup)를 일반화한다. 진짜 질문은 DB만으로 안 되는 *의미 추론*과 RAG만으로 안 되는 *정확한 조인·집계*를 둘 다 요구한다 → 둘을 합친다. (TAG-Bench: 순수 Text2SQL/RAG ~20%, 손수 짠 TAG ~55–65%.) 실행 substrate로 LOTUS식 semantic operators(`sem_filter`/`sem_agg`) 차용.

### 비준된 결정 (PI 2026-06-02)
- **백엔드 = 로컬 DuckDB**, Notion에서 materialize. Notion = SoT, refresh형(옵시디언 미러와 동일 운용). exec가 로컬이라 rate-limit 없음, 분석 SQL 강력.
- **콘텐츠 깊이 = 테이블 + PDF 풀텍스트 (완전체)**. gen이 표 필드뿐 아니라 references PDF 본문까지 semantic 추론. → references 코퍼스 수집·청킹이 v1 critical path.

---

## 1. 아키텍처

```
 Notion 9-DB (Contract-Graph, SoT)                references/*.pdf  (점진 수집)
        │  query_all() (notion API pull, 재사용)         │  fitz 텍스트 추출 + 청킹
        ▼                                                ▼
┌──────────────────────────────────────────────────────────────────┐
│  kb.duckdb  (로컬, 재생성형)                                        │
│   • 9 node tables: source_evidence, knowledge_claim, model_contract,│
│       parameter, validation_gate, code_mapping, run_result,         │
│       decision_ledger, cell_state (행=Notion 행, 컬럼=평탄화 프로퍼티)│
│   • edges(src_id, src_type, rel, dst_id, dst_type)  ← relation 평탄화 │
│   • references(citation_key, title, doi, authors, year, journal,    │
│       path, se_uid, n_pages, ingest_status, sha)                    │
│   • paper_chunks(chunk_id, citation_key, se_uid, page, text[, emb]) │
│       + DuckDB FTS 인덱스 (키워드 검색; 임베딩은 후속 업그레이드)        │
└──────────────────────────────────────────────────────────────────┘
        ▲                                                ▲
   NL 질문 R ── syn(Claude: DDL+few-shot → SQL/plan) ── exec(DuckDB[+semantic op]) ── gen(Claude → 답변 + 인용)
```

### 1.1 테이블 (DuckDB `kb.duckdb`)
- **9 node tables** — Notion 9-DB 1:1. 컬럼 = `notion_to_obsidian.plain()`로 평탄화한 프로퍼티 (title/select/multi_select/url/number/date/unique_id/formula). `notion_id`(PK), `uid`(SE33 등) 포함.
- **edges** (long form) — relation 프로퍼티를 (src, rel, dst) 트리플로. 양방향이라 한쪽만 적재(dedupe). 모든 조인의 기반.
- **references** — PDF 1행. `se_uid`로 source_evidence와 FK. dedup 키 = DOI 우선, 없으면 정규화 title.
- **paper_chunks** — PDF 본문 청크. gen의 semantic 레이어. DuckDB `fts` 확장으로 BM25 키워드 검색; 임베딩(voyage/local)은 후속.

### 1.2 엔진 (`tag_query.py`)
- **syn** — Claude에 DuckDB DDL + few-shot 주고 → (a) 순수 SQL(관계·집계형) 또는 (b) TAG plan = SQL(narrow) + semantic op(선택 행/논문 청크에 LLM map/filter/agg).
- **exec** — DuckDB에 SQL 실행; semantic op는 해당 행/청크 텍스트를 끌어와 LLM 호출.
- **gen** — Claude가 표 결과 + semantic 결과로 최종 NL 답변 작성, **SE/KB/MC/VG ID + DOI + Notion 링크 인용** 부착.
- **LLM 백엔드 (이중)** — `ANTHROPIC_API_KEY` 있으면 anthropic SDK, 없으면 `claude -p`(Claude Code 헤드리스) subprocess fallback → **새 시크릿 불필요**. 키는 선택적 권장.

### 1.3 Refresh (Notion = SoT)
`refresh.sh` — Notion 재-pull → kb.duckdb 재빌드(옵시디언 refresh와 동일 철학). references 테이블 + paper_chunks는 references/manifest + PDF에서 재생성. PI가 "TAG 갱신해줘" 하면 실행.

---

## 2. references 코퍼스 트랙 ("다른 세션처럼" autonomous ingest)

TAG gen이 읽는 **증거 콘텐츠 레이어 + `references` 테이블**을 채운다.

- **신규 ~53편** (이미 `references/`에 있음, spheroid/cancer/PIV 계열) + **기존 243편 SourceEvidence**(점진 다운로드; 페이월은 gbook KAIST 풀텍스트 [[reference-gbook-kaist-fulltext]]).
- per-PDF: 메타데이터 추출(title/DOI/authors/year — fitz + regex/Claude) → SourceEvidence와 dedup(DOI/title) → 본문 청킹 → `paper_chunks` 적재 → manifest 행 작성 → SE 행 링크(없으면 신규 SE 행 생성).
- **자기참조 앵커**: founding TAG 논문(Biswal 2024 + LOTUS Patel 2024 + 후속)을 첫 references로 받아 등록 → TAG 시스템이 자기 근거 논문을 인용.
- 폴더 정리: 중복 3건(`...S0928098725003677 (1)`, `Directed_invasion (1)`, `journal.pone.0264571 (1)`) 제거, `cellpress_bundle/` 혼합 19편 triage(유관 ~6편만 적재, 오프토픽 ~13편 격리).

---

## 3. 빌드 순서

1. `pip install duckdb anthropic`; references/ 정리(중복 + cellpress triage).
2. **`notion_to_duckdb.py`** → kb.duckdb (9 node tables + edges). *표 텍스트만으로 즉시 동작.*
3. **`references_ingest.py`** → references 테이블 + paper_chunks (fitz 추출 + 청킹 + FTS 인덱스).
4. **`tag_query.py`** syn/exec/gen (이중 LLM 백엔드) + few-shot 예시.
5. **eval**: gold NL 질문 세트(TAG-bench식) — Contract-Graph 위 syn/exec/gen 검증.
6. `refresh.sh` + README + 이 문서 유지.

---

## 4. Verification (E2E)

- **materializer**: 행 수가 Notion과 일치(SE 243 / KB 158 / MC 6 / Param 21 / VG 25 / CM 2 / RUN 1 / DEC 3); 조인 spot-check (Param→Source→VG 한 체인).
- **engine**: gold 질문 5–8개 — 생성 SQL 정확성 + 인용 ID가 실제 행으로 resolve 되는지.
- **content**: FTS 쿼리가 알려진 사실의 정답 논문 반환(예: "Bell-Evans off-rate" → Bell1978_Science).

## 5. Open deps

- `ANTHROPIC_API_KEY`: 선택(없으면 `claude -p` fallback).
- 기존 243편 PDF 다운로드: 점진, 페이월은 gbook.

관련 메모리: [[project-rag-v2-contract-graph]](Notion 9-DB IDs), [[reference-gbook-kaist-fulltext]], [[feedback-production-driver-auto-viz]].
