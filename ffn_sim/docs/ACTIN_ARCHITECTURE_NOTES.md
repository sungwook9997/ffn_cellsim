# Actin cytoskeleton architecture — primary-literature study notes (PI reference papers)

> PI 2026-06-04: "references/에 넣은 논문들을 하나하나 모든 부분 읽고 우리 구조에 넣을 것 +
> 기록 다 해두고 계속 진행." Sequential deep reads (NOT parallel skims) of the PI-added
> reference papers. Per paper: bibliographic · key architecture findings · quantitative
> parameters · **what goes into OUR model structure** · Notion SE/KC status. Feeds the unified
> actin-architecture framework ([[project-unified-actin-architecture]]) + the cortical-mesh
> construction fix. Every source here must be registered into the Notion Contract-Graph
> (SourceEvidence + KnowledgeClaim) per the PI's "all sources → Notion" directive.

**Reading queue (PI reference set, architecture-priority first):**
1. ⏳ Flormann 2024 PNAS — cortex structure/mechanics vs location & adhesion state  ← reading
2. Fritzsche 2016 Sci Adv — actin kinetics shapes cortical network structure & mechanics
3. Taeyoon Kim 2009 MIT thesis — simulation of actin cytoskeleton structure & rheology (network construction methodology)
4. Banerjee/… J Indian Inst Sci 2021 — actomyosin cortex as a thin film of active matter (review)
5. Nat Phys 2024 (s41567-024-02626-6) — energy partitioning in the cell cortex
6. Sakamoto & Murrell 2024 (Cell Rep Phys Sci) — substrate geometry → F-actin reorganization in adherent model cortex
7. Fritzsche 2017 Nat Commun (ncomms14347) — self-organizing actin patterns shape membrane architecture
8. Li, Gao & Xu 2022 Biophys J — nonlinear power-law relaxation of cell cortex (network dynamics)
9. Bächer 2021 Front Phys — 3D numerical model of active cell cortex (viscous limit)
10. Garlick 2022 Sci Rep — quantifying super-resolved cortical actin
11. Ray 2024 — actin capping protein regulates actomyosin contractility (germline architecture)
12. Nat Commun 2024 (46726) — kinetic trapping organizes actin filaments in droplets
13. Merino-Casallo 2022 — cell migration from the cell surface (review)
14. + cortical-tension batch (Chugh emss-72183/ncb3525, Warmt, Dmitrieff, Winklbauer jcs174623, Bohec, Murrell-Gardel) — partly extracted; deepen.

---

## 1. Flormann et al. 2024, PNAS 121(31):e2320372121 — "The structure and mechanics of the cell cortex depend on the location and adhesion state" (`flormann-et-al-2024-...pdf`)

**Cells:** hTERT-RPE1 (interphase) + HeLa confirmation. **Methods:** SEM + FiNTA mesh-tracing
(mesh hole area MHA), expansion microscopy (thickness, side-view), AFM creep-compliance
(stiffness, ~400 nm indentation = cortex not stress fibers), fluorescence (actin/myosin amount).
**Comparison:** suspended vs adhered; within adhered: nuclear vs perinuclear region.

### Key architecture findings
- **Cortex architecture is NOT one fixed thing — it depends on ADHESION STATE + location.**
  - **Suspended**: THICKER cortex, SMALLER mesh, denser actin, FEW bundles, LOWER stiffness.
  - **Adhered**: THINNER cortex, LARGER mesh, MORE bundles (most in perinuclear), HIGHER stiffness.
  - nuclear ≈ perinuclear thickness; perinuclear has larger mesh + more bundles + stiffer.
- **Positive mesh-size ↔ stiffness correlation in LIVING cells (Pearson R=0.96)** — OPPOSITE to the
  naive in-vitro "smaller mesh = stiffer." Cause: in cells a larger mesh comes WITH more
  cross-linking/**bundling** → thicker bundles (higher bending rigidity) → stiffer. **Stiffness is
  governed by BUNDLING (crosslink-induced), not actin density alone.**
- **Myosin is NOT the bundler** (blebbistatin didn't change bundling; suspended cells have MORE
  myosin yet FEWER bundles). Cross-linkers/bundling proteins do it. (Consistent w/ Chugh.)
- Latrunculin A (depolymerization): adhered → larger mesh + softer (breaks thin filaments →
  lowers connectivity); suspended → mesh robust (thick cortex), just thinner.

### Quantitative anchors
- **hξ⁻² = (cortex thickness)/(mesh hole area) = total F-actin length per unit area** =
  **0.07 nm⁻¹ nuclear, 0.06 nm⁻¹ perinuclear, 0.10 nm⁻¹ suspended.** (constructable invariant.)
- Thickness ~ Clark 2013 range (~200 nm; suspended thicker). ⚠️ exact MHA(nm²)/thickness(nm) live
  in Fig 2E/3A panels (SI) — pull from SI if a precise mesh-size band is needed.
- **Semiflexible cross-linked biopolymer theory (MacKintosh 1995; Gardel/Shin in vitro):**
  - `G ~ K_B² / (kB T · ξ⁵)`, bundle bending rigidity `K_B ~ D_B⁴` (D_B = bundle thickness),
    `D_B ~ ([crosslink]/[actin])^0.3`, mesh `ξ ~ D_B / [actin]^{1/2}`.
  - Constant actin ⇒ **G ~ D_B³ ~ ξ³** (the positive mesh–stiffness correlation).

### → WHAT GOES INTO OUR MODEL
1. **Physiological-baseline (adhesion state matters):** the cortex must be built for the cell's
   actual state. MCF7 in a spheroid = cell–cell ADHERED/cohesive (not suspended) → use the
   ADHERED architecture (thinner, larger-mesh, MORE-bundled, stiffer) — NOT a generic shell.
   This is a concrete instance of [[feedback-physiological-baseline]].
2. **Cross-linkers must BUNDLE**, and stiffness must emerge from bundling, not density alone. Our
   crosslinker model should permit multi-crosslink bundle formation (D_B grows with crosslink:actin).
3. **Acceptance ORACLE (not runtime):** `G ~ K_B²/(kBT ξ⁵)` and the **positive G–ξ correlation at
   constant actin** = a validation gate for the rebuilt connected mesh (the fragmented mesh would
   fail this — no bundling, wrong G–ξ sense). Cross-read with Chugh (tension) + Head 2003 (L/lc).
4. **hξ⁻² ≈ 0.06–0.10 nm⁻¹** = a constructable invariant tying thickness · mesh · filament length —
   use to check the rebuilt mesh's areal contour-length density is physiological.
5. Myosin ≠ bundler → keep bundling in the crosslinker layer, not myosin (matches current split).

**Notion:** SourceEvidence row NEEDED (Flormann2024_PNAS; DOI 10.1073/pnas.2320372121; gbook/KAIST
full-text present). Link to KnowledgeClaim: cortex-architecture (mesh/thickness/bundling↔stiffness)
+ the adhesion-state dependence. ⚠️ not yet in corpus — web-verify + register.

**(감상평 below)**

## 2. Fritzsche, Erlenkämper, Moeendarbary, Charras, Kruse 2016, Sci Adv 2:e1501337 — "Actin kinetics shapes cortical network structure and mechanics" (`sciadv.1501337.pdf`)

**Cells:** HeLa (cervical cancer) + M2 melanoma. **Methods:** FSM single-molecule (formin Diaph1,
Arp2/3) + FRAP + stochastic simulation + AFM (Hertz). The bottom-up "how nucleation+turnover set
the cortex architecture" paper.

### Key architecture findings (the filament-length answer)
- Cortex = roughly **ISOTROPIC** semiflexible F-actin network, crosslinked, with motors. Nucleated
  by **TWO pathways: Arp2/3 (branched) + formin Diaph1 (linear)** — these two alone make most
  cortical actin (Bovellan 2014).
- **Filament length is BIMODAL, EXPONENTIALLY distributed** (shorter filaments more abundant, from
  turnover): **Arp2/3 ≈ 120 nm (HeLa) / 60 nm (M2); formin ≈ 1200 nm (HeLa) / 600 nm (M2)** — formin
  filaments are **~10× longer** than Arp2/3. (Single formins add avg 3900 nm HeLa / 2500 nm M2.)
- **<10 % of filaments are formin-nucleated BY COUNT, but they are the MAIN MECHANICAL contributors**
  (because ~10× longer). The long sparse formin filaments dominate cortex mechanics.
- Two subpopulations come from **two nucleation pathways** (turnover rates differ ~20×); severing
  ALONE cannot generate the two distinct timescales.
- Cortex elastic modulus **E = 4.0 ± 1.5 kPa (HeLa control)**, 3.0 ± 1.7 kPa perturbed (AFM/Hertz).
- Kinetic anchors: r_on = 9 s⁻¹mM⁻¹, r_on,F = 45 s⁻¹mM⁻¹, formin detach w_off,F = 0.12 s⁻¹.

### → WHAT GOES INTO OUR MODEL
1. **Filament length = bimodal exponential, NOT uniform.** Seed Arp2/3-short (~120 nm, ~90% count)
   + formin-long (~1200 nm, ~10% count). At ×40 mesoscale, conserve the SHAPE (bimodal exp) +
   total contour (Chugh invariant, cortex-workflow FIX 3).
2. **The LONG formin filaments are the PERCOLATION + mechanical backbone.** A ~1.2 µm filament
   spans many mesh holes → reaches many distinct crosslink partners → ties the network together.
   Our fragmented mesh (z=1.3) likely lacks this long connecting subpopulation (uniform L=3 µm
   mis-represents it). **Including the long-formin subpopulation is part of the percolation fix**,
   complementary to the bridge-different-filaments crosslink rule.
3. Construction is STATIC (physiological-baseline) → seed the exponential bimodal length DIRECTLY
   (the distribution turnover would produce), don't rely on dynamic turnover to grow it.
4. **Cortex E ≈ 3–4 kPa (HeLa)** = stiffness anchor (consistent with Flormann AFM); an acceptance
   cross-check (overlay; HeLa not MCF7).
5. Kinetic rates (r_on, w_off,F) → feed the existing turnover module if dynamic turnover is on.

**감상평 (Lead):** 필라멘트-길이 아키텍처의 정량적 정답. 가장 큰 통찰: **연결망을 잇는 건 "크로스링크를
더"가 아니라 "긴 formin 필라멘트(개수 10%, 길이 ~1.2µm)가 역학·연결의 backbone"**이라는 것. 짧은
Arp2/3(120nm)만으로는 percolate가 약하고, 긴 formin 필라멘트가 여러 mesh hole을 가로질러 서로 다른
필라멘트들을 묶음. 우리 fragmented mesh는 이 long-subpopulation이 없/부족하거나 uniform-length로
뭉개버린 게 원인일 수 있음 → cortex-workflow의 "contour 보존 + bimodal" FIX와 정확히 합치. 즉
**percolation 수정 = (a) bridge-different-filament 규칙 + (b) 긴-formin 포함한 bimodal 길이분포** 둘 다.
또 "formin 10%가 역학 지배"는 Chugh의 "intermediate length가 tension 최적"과도 연결(긴 필라멘트가
force 전달). 한계: HeLa/M2 (MCF7 아님), turnover-kinetics 논문이라 우리가 원하는 건 그 함의(정적
길이분포). E~4kPa는 Flormann과 일치 = 좋은 cross-check.

## 3. Taeyoon Kim 2007, MIT MS thesis — "Simulation of Actin Cytoskeleton Structure and Rheology" (`181655768-MIT.pdf`, 153 pp)

**What it is:** 3D **Brownian-dynamics** model — actin monomers polymerize into filaments, cross-
linked by **two ACP types: PERPENDICULAR (large/long, e.g. filamin → isotropic NETWORKS) vs
PARALLEL (small/short → BUNDLES)**. Evaluates how parameters set network morphology. The
methodological backbone for crosslinked-actin-network simulation (Kim's foundational work; cf.
later Kim 2009 Biophys J, Kim 2014 — the mature "Kim model" / AFINES lineage).

### Key findings (methodology + connectivity logic)
- **ACP binding-site geometry decides bundle vs network**: parallel-binding ACPs → bundles
  (fascin-like); long ACPs forming ~perpendicular cross-links (filamin) → isotropic networks.
- Network morphology (pore size, isotropy, extent of cross-linking) is set by **crosslinker
  concentration + type + actin concentration** (R = ACP:actin ratio, CA = actin conc, Da).
- **CONNECTIVITY + PERCOLATION are explicit, tunable outputs**: distribution of per-filament
  connectivity (counts of crosslinks/filament), "connectivity 2", and a network that "nearly
  percolates the simulation box" at given (Da, CA, R). Pore size Lpore ∝ crosslink spacing;
  Lm (mean segment between crosslinks) = alt pore-size measure.
- F-actin 7–9 nm diameter (excluded volume).

### → WHAT GOES INTO OUR MODEL
1. **Construction logic = ours**: a crosslinked network's connectivity/percolation/pore-size are
   controlled OUTPUTS of crosslinker concentration + type. Validates measuring z + giant-component
   (our viz_cortex_network) and tuning crosslink density to the percolation set-point.
2. **Perpendicular(network, filamin) vs parallel(bundle, fascin) ACP distinction** = the unified
   crosslinker/bundler axis: cortex = perpendicular/isotropic network (filamin/α-actinin);
   filopodia/microvilli = parallel bundles (fascin/espin). One model, ACP-type parameterized.
3. Self-assembly produces percolation at the right (R, CA); our STATIC construction should SEED
   directly at that percolated morphology (physiological baseline), then optionally turn dynamics on.

**감상평 (Lead):** 방법론적으로 가장 가까운 동족. Taeyoon Kim은 actin-network BD 시뮬레이션의 표준을
세운 사람이라(이후 Kim 2009/2014 = 성숙한 모델, AFINES 계보), 우리 접근이 정통임을 확인해 줌. 가장
유용: **perpendicular(network) vs parallel(bundle) ACP 이분법** — 이게 통합 프레임워크의 crosslinker
축 그 자체(cortex=filamin 수직망, filopodia=fascin 평행다발). 그리고 **connectivity/percolation이
crosslinker 농도·종류의 tunable 출력**임을 직접 보여줘 우리 viz의 z·giant-component 측정·튜닝을
정당화. 한계: 2007 MS 논문(방법론·in-silico, 세포-특이 수치 아님) → HOW-TO + acceptance-logic 소스.
**후속 Kim 2009 Biophys J / Kim 2014를 추가로 끌어오면** 성숙한 파라미터(crosslink stiffness, prestrain,
network self-assembly protocol)를 얻을 수 있음 — references에 없으면 gbook/web으로.

## 4. (Banerjee et al.) 2021, J. Indian Inst. Sci. — "The Actomyosin Cortex of Cells: A Thin Film of Active Matter" (`s41745-020-00220-2.pdf`) — REVIEW (continuum active-gel theory)

**What it is:** review of the **hydrodynamic active-gel theory** of the cortex (Kruse-Jülicher-
Joanny-Prost lineage): ATP-driven myosin generates active stress → large-scale mechanical FLOWS +
orientation/mechanochemical PATTERNS. Continuum, not molecular-architecture numbers.

### → WHAT GOES INTO OUR MODEL (mostly framing / acceptance at the continuum limit)
- The cortex is an **active contractile gel**: coarse-grained, our fine-grained network should
  reproduce its active-gel behavior (active stress, cortical flows). This is the CONTINUUM target
  our mechanistic model maps onto — an acceptance-level cross-check, not a parameter source.
- Confirms force-generation = ATP-myosin coordinated → net active stress (ties to the contractility
  question), and that crosslinkers + motors + filaments together set the active-gel parameters.

**감상평 (Lead):** 연속체 active-gel 이론 리뷰 — 우리 fine-grained 모델이 coarse-grain하면 닿아야 할
"정답 거시 거동"(active stress, cortical flow)을 줌. 단 **분자 아키텍처 수치(길이/메시/밀도)는 없음** →
construction 파라미터 소스가 아니라 framing/acceptance 소스. 우선순위는 낮되, "cortex=active gel"
프레임은 통합 framework의 motor/contractility 축 근거로 등록 가치 있음. (Marchetti 2013, Prost-Jülicher-
Joanny 2015 active-gel 원전이 더 1차적.)

## 5. Chen, Seara, … Bement, Murrell 2024, Nat Phys 20:1824 — "Energy partitioning in the cell cortex" (`s41567-024-02626-6.pdf`)

**What it is:** non-equilibrium thermodynamics of the cortex — entropy-production rate of the
CHEMICAL (Rho-GTPase/actin/myosin) vs MECHANICAL subsystems across pattern regimes (pulses →
choppy waves → labyrinthine/spiral), tuned via Rho-GAP. Onsager reciprocity holds at low drive,
breaks at high drive. **Key principle: energy partitioning + chemical↔mechanical coupling are set
by the COMPETING TIMESCALES of chemical reaction vs mechanical relaxation.** (Xenopus/starfish-type
cortex with Rho waves; not molecular architecture.)

### → WHAT GOES INTO OUR MODEL
- Framing/validation, not construction params. The **chemical-reaction-vs-mechanical-relaxation
  timescale competition** is exactly the regime our binding/turnover-vs-BAOAB-relaxation
  accelerated-dynamics probes live in — a principled caution that the chem/mech timescale ratio
  governs the emergent behavior (don't distort it when accelerating). Active/non-equilibrium frame.

**감상평 (Lead):** 멋진 물리지만 우리 construction엔 직접 파라미터 없음(패턴 열역학). 단 한 줄이 값짐:
**"energy partitioning은 화학반응 vs 역학완화 timescale 경쟁이 결정"** — 우리가 binding/turnover를
가속할 때 chem/mech timescale 비를 왜곡하면 emergent 거동이 바뀐다는 경고와 정확히 같은 물리. Murrell
그룹(§6 Sakamoto-Murrell와 같은 lab) active-cortex 라인. 우선순위 낮음, framing/acceptance로 등록.

---

**감상평 (Flormann, §1):** 우리에게 결정적으로 유용. 두 가지가 큼. (1) **"cortex 아키텍처는 단일 상수가
아니라 adhesion-state의 함수"** — 우리가 MCF7 spheroid를 모델링하면서 generic shell을 쓴 게
바로 physiological-baseline 위반의 교과서 사례. spheroid 내 세포는 cell-cell 접착 상태 → 그 상태의
아키텍처(두께/메시/번들)를 써야 함. (2) **living-cell에선 mesh ↑ ↔ stiffness ↑ (양의 상관, R=0.96)**
— in-vitro 단순 직관(메시 작을수록 뻣뻣)과 반대고, 그 이유가 **번들링(crosslink-induced)**이라는 게
우리 fragmented-mesh 진단과 정확히 맞물림: 우리 cortex는 번들도 못 만들고(crosslink가 same-filament
staple) percolate도 안 됨 → 이 논문의 G~ξ³(번들 기반)을 재현할 수 없음. 즉 **이 논문의 G–ξ 양상
재현이 "연결된 그물 + 번들링"의 강력한 acceptance gate**가 됨. 한계: hTERT-RPE1/HeLa (MCF7 아님 —
overlay-only), 절대 mesh/두께 nm는 figure에 있어 SI 필요. semiflexible-network scaling은 in-vitro
유래라 living-cell엔 정성적 일치까지가 정직한 수준 (저자도 "aspects fit").

---
