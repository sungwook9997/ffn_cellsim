# 다운로드 필요 논문 — 2026-06-02 (compartment-audit 워크플로 산출)

> **PI 워크플로**: PI가 PDF를 받아 이 폴더(`aleph/references/`)에 그대로 드롭 →
> Lead가 검사(저자/연도/내용 일치) 후 references에 편입(`analysis/` 추출 +
> `tag_corpus.json` 갱신) → 진행. 아래는 **critical-path 우선순위** 정렬.
>
> 출처: `wf_e87b5c75` (compartment-stack-audit) synthesis. 12편 전부 현재 references에
> **없음**(저자-키워드 grep 확인). 접근성 표기: OPEN = 무료, paywall = 구독/회수 필요
> (gbook KAIST 풀텍스트 또는 ResearchGate 저자본). DOI로 검색.

## TIER 1 — γ critical path (지금 필요)

이 4편이 estimator-fix(step 1)와 KU-3.5 band/turnover 근거(step 3)를 *문헌 계약*으로
못박는다. 현재 코드는 method-of-planes 부호 규약을 출처 없는 "Codex recommendation"으로만
인용 중 — 반드시 1차 문헌으로 고정해야 함.

| # | 논문 | DOI / 위치 | 접근 | 왜 필요 |
|---|---|---|---|---|
| 1 | **Irving & Kirkwood 1950**, *J Chem Phys* 18:817 | 10.1063/1.1747782 | paywall (널리 구함) | method-of-planes/virial 응력의 정본 정의 → `tension.py:109/257` 부호 규약 수정의 문헌 계약. **최우선.** |
| 2 | **Todd, Evans & Daivis 1995**, *Phys Rev E* 52:1627 | 10.1103/PhysRevE.52.1627 | paywall | 비균질 유체 method-of-planes (IK와 페어). Custom-force가 virial에 들어가는 방식 포함. |
| 3 | **Salbreux, Charras & Paluch 2012**, *Trends Cell Biol* 22:536 | 10.1016/j.tcb.2012.07.001 | paywall | KU-3.5 [0.35,0.65] mN/m cortex band 1차 앵커. |
| 4 | **Chugh & Paluch 2017**, *Nat Cell Biol* 19:689 | 10.1038/ncb3525 | paywall (RG 저자본) | cortex surface tension 결정인자 + **actin turnover = steady-state tension의 단일 지배인자**(step 3 load-bearing). |

## TIER 2 — composite 재귀속 (step 4, cortex가 band 든 *후*)

| # | 논문 | DOI / 위치 | 접근 | 왜 필요 |
|---|---|---|---|---|
| 5 | **Fischer-Friedrich et al. 2014**, *Sci Rep* 4:6213 | 10.1038/srep06213 | **OPEN** | 단일 유사분열 세포의 surface tension + 내부 압력 정량 → composite(γ_membrane+γ_cortex) + Laplace 경로 근거. |
| 6 | **Sens & Plastino 2015**, *J Phys Condens Matter* 27:273103 | 10.1088/0953-8984/27/27/273103 | paywall | 막 장력 vs cytoskeleton; 막은 소수 기여(~0.003–0.04 mN/m) 확인 → "cortex-only가 맞는 타깃" 뒷받침. |

## TIER 3 — 미래 compartment gate (nucleus/cytoplasm, γ 경로 밖)

| # | 논문 | DOI / 위치 | 접근 | 왜 필요 |
|---|---|---|---|---|
| 7 | **Stephens, Banigan, Marko et al. 2017**, *Mol Biol Cell* 28:1984 | 10.1091/mbc.E16-09-0653 | **OPEN** (PMC) | `nucleus.py`의 two-regime(chromatin/lamin) bilinear 법칙 1차 출처. H.9 magic-number 근거. |
| 8 | **Stephens et al. 2018**, *Nucleus* 9:119 | 10.1080/19491034.2017.1414118 | **OPEN** (PMC) | 위 후속(chromatin vs lamin 역할 분리). |
| 9 | **Hu et al. 2024**, *Nanoscale Adv* | PMC10929591 | **OPEN** (PMC) | cytoplasm 점도 대비(MCF7 ≈ 5× MDA) = H.10 Tier-1/2 검증 타깃. |
| 10 | **Moeendarbary et al. 2013**, *Nat Mater* 12:253 | 10.1038/nmat3517 | paywall (RG) | poroelastic 밴드(Dp 40–60 µm²/s, pore ~14 nm) = H.10 Tier-3. |
| 11 | **Davidson, Denais, Lammerding 2014**, *Cell Mol Bioeng* 7:293 | 10.1007/s12195-014-0342-y | paywall (RG) | 핵 변형성 = 3D 이동 율속 → KU-1.V.3.3 confined-migration arrest 앵커. |
| 12 | **Wolf et al. 2013**, *J Cell Biol* 201:1069 | 10.1083/jcb.201210152 | **OPEN** (PMC) | pore-size 한계(~핵 단면 10%) — 위와 페어. |

## 권장 회수 순서

1. **지금**: #1 Irving-Kirkwood 1950 (estimator-fix 직행). #3/#4 (band/turnover).
2. OPEN(#5,7,8,9,12)은 PI가 빠르게 받거나, 원하면 Lead가 gbook KAIST로 회수 시도.
3. paywall(#1,2,4,6,10,11)은 gbook KAIST 풀텍스트 경로(ssh gbook + 브라우저 UA curl) 가능 —
   PI가 직접 받는 게 빠르면 그쪽, 막히면 Lead가 시도.

> 드롭 후 한 줄만 알려주면(또는 그냥 폴더에 넣으면) Lead가 파일명↔논문 매칭 검사하고
> `analysis/` 추출 + tag_corpus 편입까지 처리한다.
