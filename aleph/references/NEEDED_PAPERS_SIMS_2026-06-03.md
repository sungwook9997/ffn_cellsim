# 다운로드 필요 — actomyosin network-remodeling SIMULATION 논문 (2026-06-03)

> STAGE-2 아키텍처 결정의 근거. 자동 다운로드 실패(PMC Cloudflare 봇차단 + biorxiv 차단,
> 로컬·gbook 둘 다 403/에러페이지). **PI가 수동 다운로드 → `aleph/references/`에 드롭**
> → Lead가 검사 후 편입. 파일명 prefix는 `sim_` 권장.
>
> 맥락: 우리 플랫폼이 막힌 곳(active cortical tension, transmission/ratchet)을 분야가
> 어떻게 시뮬하는지. **핵심 결론: turnover + filament BUCKLING이 disordered-network
> contractility의 두 축 → 우리 rigid(constrained) backbone이 둘 다 차단하고 있었음
> → STAGE-2는 semi-flexible 모드로.**

| # | 논문 | DOI / ID | 접근 | 왜 (STAGE-2 근거) |
|---|---|---|---|---|
| 1 | **Hiraiwa & Salbreux** — Filament turnover tunes force generation & dissipation, long-range flows in a model actomyosin cortex | PMC5757993 / PMID 29253848 | OA(PMC) | **turnover가 압축 filament 교체 → sustained stress; 최적 율 존재** → 우리 τ_half/f_ss sweet-spot의 직접 근거 |
| 2 | **Belmonte, Leptin, Nedelec 2017** — A theory that predicts behaviors of disordered cytoskeletal networks | 10.15252/msb.20177796 / PMID 28954810 | OA(EMBO) | active-gel은 contractility 설명 못 함, motor+filament 개별로 봐야 → **fine-grained 접근 검증**; Cytosim |
| 3 | **Linsmeier et al. 2016** — Disordered actomyosin networks sufficient for cooperative & telescopic contractility | 10.1038/ncomms12615 / PMC5007339 | OA(NatComm) | disordered network → net contraction (buckling 비대칭) |
| 4 | **Popov, Komianos, Papoian 2016** — MEDYAN: mechanochemical sim of contraction & polarity in actomyosin networks | PMC4847874 | OA(PLoS) | 분야 표준 mechanochemical 시뮬레이터 (remodeling 포함) |
| 5 | **Ennomani et al. 2016** — Architecture and Connectivity Govern Actin Network Contractility | PMC4959279 | OA(CurrBiol) | **connectivity(percolation)가 contractility 지배** → 우리 percolation 우려 직결 |
| 6 | **(Floyd/Banerjee 계열) 2021** — Protein Friction & Filament Bending Facilitate Contraction of Disordered Actomyosin Networks | biorxiv 10.1101/2021.02.23.432588 | OA(biorxiv) | **filament BENDING/buckling이 contraction 촉진** → rigid backbone이 막던 그 메커니즘 |
| 7 | **Hiraiwa & Salbreux 2018(?)** — Actomyosin pulsation & flows in an active elastomer with turnover & network remodeling | PMC5783953 | OA(PMC) | turnover + remodeling → pulsation/flow |
| 8 | **(2025)** — Non-uniform filament turnover & mechanically-driven contractility/bundle formation in disordered actomyosin networks | biorxiv 10.1101/2025.11.05.686892 | OA(biorxiv) | 최신 turnover-contractility |

**우선순위**: #1(turnover-tension 직접 근거), #6(buckling), #5(connectivity/percolation), #2(fine-grained 정당화) — 이 4편이 STAGE-2 설계의 핵심.

> 드롭 후 Lead가 파일명↔논문 검사 + `analysis/` 추출 + tag_corpus 편입.
