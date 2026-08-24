# MDA-MB-231 Baseline Contract — 2026-07-13

**Purpose.** Formalize MDA-MB-231 as the **first held-out validation cell** after MCF7. MCF7 stays the
development/calibration anchor; MDA-MB-231 is **never tuned to its own targets**. This document is an
**audit + decision surface only** — no runtime physics (γ / ΔP / adhesion / NMII) is changed here
(PI-gated, on hold). The R / η / E_nuc wiring is already PI-ratified (2026-07-13) and kept.

**Provenance.** The MDA-MB-231 task plan originated from **Codex** (idea-generator, not a verification
authority — [[feedback-codex-debug-not-judge]]). Every number Codex supplied was **re-verified from the
PRIMARY source** in workflow `wf_35f1b9e4-d78` (10 agents, PubMed/PMC full-text opened). Corrections are
listed in §1. Nothing below is asserted without a confirmed DOI.

---

## 1. Codex-claim verification (primary-source)

| Codex claim | Verdict | Correction |
|---|---|---|
| Moazzeni **"Kaufman"** 2021, PRE 103 032409, cortical prestress 10⁻²–10⁻¹ N/m @0.01s | **number ✅ / cite ✗** | No author "Kaufman" — authors are **Moazzeni, Demiryurek, …, Lin**. Number confirmed. |
| Tsujita 2021 MDA tether tension 45–50 pN/µm | **✅ (inferred)** | Paper prints "epithelial ~100 pN/µm, metastatic ~2× lower" → MDA ~50 pN/µm is an **inference from Fig 1**, not a verbatim MDA figure. Raw observable = tether **force ~20–25 pN**. |
| Zhovmer/Tabdanov **2022** ACS Nano, rounded-MDA AFM cortical tension | **✅ / year ✗** | It is **2021**, ACS Nano 15(11):17528, DOI 10.1021/acsnano.1c04435. Control (+DMSO) = **~290–300 pN/µm ≈ 0.30 mN/m** (Fig 4c, n=33). |
| Fischer-Friedrich 2014 HeLa: interphase ΔP≈40 γ≈0.2; metaphase ΔP≈400 γ≈1.6 | **✅** | Verbatim: interphase γ=0.17±0.13 mN/m, ΔP=40±30 Pa; metaphase γ=1.6±0.5, ΔP=400±120. HeLa = **proxy**, not 231. |
| MDA suspended R = **8.8 ± 1.3 µm** | **✗ (unverified aggregate)** | No traceable primary ensemble. Direct: **Cognart 2020 R≈7.5 µm** (dia ~15µm). 8.8±1.3 ≈ arithmetic mean of bracketing values, not a measured ensemble. Preset R=8.0 is fine (band 7.5–10). |
| cortex thickness **~500 nm** | **✗ (wrong + non-MDA)** | No MDA-specific cortex thickness exists; real cortical actin **~190 nm** (Clark/Paluch 2013, range 100–200). 500 nm is ~2.5× too thick. |
| η 12.0±5.7 (pooled) vs 10.7±5.4 (L=3µm) | **✅ both** | Dessard 2024. Paper's **primary quotable = 10.7±5.4** (Table 1); 12.0 = all-wires pooled. Preset uses 12.0 → decision-surface item 3. |
| G 38.6±5.8 Pa; τ 0.57±0.20 s | **✅ both** | Dessard 2024 (τ=η/G, cytoplasm). |
| E_nuc 157.7±78.6; E_cyto 103.4±89.5 Pa | **✅ both** | Fischer, Hayn & Mierke 2020 (over-nucleus / over-cytoskeleton regions, n=68). |
| Warmt "0.24–0.39 Pa" MDA cortical | **✅ / relabel** | Exactly Warmt 2021 Table 1 MDA row (σ_int 0.39 median / 0.24 relaxation). It is an **active contractile STRESS (Pa)**, optical-stretcher, ~5s — **NOT a surface tension (mN/m)**. |

**Headline correction to the prior draft:** the draft's *"MDA cortical tension has no direct value"* is
**partly refuted** — three DIRECT MDA-MB-231 tension measurements exist, but they are **different physical
quantities** (§2). What genuinely does NOT exist is a **minutes-scale slow-equilibrium** MDA γ or a direct
MDA turgor ΔP.

---

## 2. Tension comparison table — DO NOT CONFLATE (deliverable B)

Every "tension" for MDA-MB-231 is a different physical quantity set by its protocol, force-application
time, and cell state. Only the **slow-equilibrium pressure-bearing total surface tension** corresponds to
FF/DCM baseline **γ = ΔP·R/2**.

| Quantity class | Cell | Value | Protocol | Force timescale | State | → FF γ? | Source (DOI, verified) |
|---|---|---|---|---|---|---|---|
| **tether / membrane tension** | MDA | **~0.05 mN/m** (~50 pN/µm; tether force ~20–25 pN) | optical-tweezer tether pull | transient (~0.5s contact, ~5s pull) | adherent | **NO** | Tsujita 2021, 10.1038/s41467-021-26156-4 |
| **subsecond cortical prestress** | MDA (+3 lines, no per-line breakout) | **~10 mN/m** (10⁻² N/m) @0.01s; 10⁻²–10⁻¹ over 0.01–10s | electrodeformation relaxation | **0.01 s** (subsecond) | suspended | **NO** | Moazzeni 2021, 10.1103/PhysRevE.103.032409 |
| **AFM cortical-shell surface tension** | MDA (rounded) | **~0.30 mN/m** (290–300 pN/µm, n=33) | AFM tipless liquid-drop force balance | **seconds** (dynamic force curve) | rounded / weakly-adhered | **~YES** (seconds-scale, apparent) | Zhovmer/Tabdanov 2021, 10.1021/acsnano.1c04435 |
| **slow-equilibrium total surface tension + ΔP** | **HeLa (PROXY)** | interphase **γ 0.17±0.13 mN/m, ΔP 40±30 Pa**; metaphase 1.6 / 400 | AFM parallel-plate confinement (held minutes) | **minutes** (quasi-static) | rounded / confined | **YES** (the FF baseline proxy) | Fischer-Friedrich 2014, 10.1038/srep06213 |
| **active cortical contractile stress** | MDA | σ_int **0.39 Pa** median (0.24 relax); "predominantly passive" | optical stretcher | ~5 s | suspended | **NO** (stress≠tension) | Warmt 2021, 10.1088/1367-2630/ac254e |
| **whole-cell apparent modulus** | MDA | 224 Pa (AFM) / 206 Pa (aspiration) | AFM / micropipette | — | adherent/susp. | **NO** (modulus) | Zbiral 2023 / Lee&Liu 2015 (validation only) |

**Cross-check (sanity anchor).** FF baseline γ = ΔP·R/2 = 40 Pa · 8 µm / 2 = **0.16 mN/m**. The direct
MDA AFM cortical tension (Zhovmer) is **0.30 mN/m ≈ 1.9× higher** — same order, good anchor, and the gap is
a real modelling question (higher ΔP, or active myosin on top of passive γ). **Never** inject Moazzeni's
~10 mN/m or Tsujita's ~0.05 mN/m as runtime γ — they are the wrong physical quantities.

---

## 3. MDA-MB-231 INPUT / VALIDATION / GAP matrix (deliverable C)

**Validation-only targets are NEVER used for runtime parameter fitting.**

### 3A. INPUTS (physiological setpoints — sourced)
| Input | Value | Status | Source |
|---|---|---|---|
| suspended R | **8.0 µm** (band 7.5–10) | ✅ WIRED | Cognart 2020 (dia ~15µm → R7.5); 8.0 within band |
| cytoplasm η | **12.0 Pa·s** (pooled) / 10.7 (L=3µm primary) | ✅ WIRED (12.0) | Dessard 2024 — **decision item 3** |
| cytoplasm G | 38.6±5.8 Pa | informs η/G | Dessard 2024 |
| cytoplasm relaxation τ | 0.57±0.20 s | consistency | Dessard 2024 |
| E_nuc (in-situ over-nucleus) | **157.7±78.6 Pa** | ✅ WIRED (band 200→150 ratified) | Fischer 2020 |
| cortex thickness | **no MDA value**; ~190 nm generic | not wired / GAP | Clark 2013 (proxy) |
| equilibrium γ_cortex | **emergent** = ΔP·R/2 = 0.16 mN/m | GAP (runtime, PI-gated) | direct MDA target 0.30 (Zhovmer) |
| turgor ΔP | **40 Pa (HeLa proxy)**, no MDA value | GAP (runtime, PI-gated) | Fischer-Friedrich 2014 (HeLa) |
| integrin/FA count, adhesion kinetics | unresolved | GAP (runtime, PI-gated) | — |
| NMII / cortex areal density | unresolved | GAP (runtime, PI-gated) | — |
| cadherin identity (N-cad) | unresolved | GAP (runtime, PI-gated) | — |

### 3B. VALIDATION-ONLY OUTPUTS (compare, never fit)
| Target | Value | Source |
|---|---|---|
| whole-cell AFM E | 224 Pa; E_inst/E_inf 262/105 Pa; τ1/τ2 0.18/4.21 s | Zbiral 2023 |
| E_cytoskeleton (over-cytoskeleton) | 103.4±89.5 Pa | Fischer 2020 |
| 2D migration speed | 0.57±0.03 (plastic) / 0.48 µm/min | Holme 2023 / Biswenger 2018 |
| 3D collagen speed | 0.39±0.05 µm/min | Holme 2023 / Kurniawan 2015 |
| persistence | ~0.4 | Kurniawan 2015 |
| 2D traction \|F\| | 90 / 305 / 375 nN @ 1/5/10 kPa | Kraning-Rush 2012 |
| 3D contractility | 47.6±3.4 nN (Steinwachs) / 70.2 (Cóndor) | Steinwachs 2016 / Cóndor 2019 |
| AFM cortical tension (rounded) | 0.30 mN/m | Zhovmer 2021 — **validate emergent γ against this** |
| subsecond cortical prestress | ~10 mN/m @0.01s | Moazzeni 2021 — different-timescale check |
| morphology | AR 4.42, spread 598 µm², 3D vol 4120–5823 | Liew 2024 / RodríguezCruz 2024 |
| MMP / pore | MT1-MMP absolute (cross-linked collagen); arrest <~7 µm² | Sabeh 2009 / Wolf 2013 |

### 3C. GAPS (runtime, ALL PI-gated → on hold per 2026-07-13 directive)
γ_cortex handling · turgor ΔP · contractility (NMII) · adhesion (FA/integrin/clutch) · cadherin · cortex
thickness/areal-density MDA-specific.

---

## 4. Validation ladder G0–G6 (deliverable D — refined from Codex)

Gates are **prediction vs measurement**; migration shortfalls are reported honestly, never tuned away.

| Gate | Scope | Validate against | Pass = | Note |
|---|---|---|---|---|
| **G0** | HeLa rounded-cell absolute mechanics **oracle** | force–shape–Laplace: ΔP & γ jointly | reproduce Fischer-Friedrich 0.17 mN/m / 40 Pa | proxy oracle; unblocks the γ=ΔP·R/2 self-consistency |
| **G1** | MDA **resting mechanics** (R, η, E_nuc only) | equilibrium γ (vs Zhovmer 0.30), relaxation τ 0.57s, whole-cell AFM E 224 Pa | in-band, no migration yet | R/η/E_nuc already wired |
| **G2** | MDA **2D traction** | traction 90/305/375 nN vs stiffness; per-clutch distribution | prediction in band | **no n_fa/myosin/drag sweep to fit** |
| **G3** | MDA **2D migration** | speed 0.48–0.57, persistence ~0.4, morphology AR 4.4 | report gap + decompose | current ~0.22 nm/s = **~30–43× low**; native η-run in progress updates this |
| **G4** | MDA **3D collagen** | density/pore-size, nuclear deformation, MMP ON/OFF, remodeling | prediction | E_nuc 157.7 matters here |
| **G5** | **lineage generalization** | NIH/3T3 migration–traction (independent SF/FA check) | prediction | non-epithelial lineage |
| **G6** | **final challenge** | HT1080: MT1-MMP tunnel initiation + MMP-independent migration in existing tunnel | prediction | 3D proteolytic invasion |

DCM is **not** used for MDA migration until active rest-area generation, substrate adhesion, clutch
traction, and cortical-tension anisotropy are closed.

---

## 5. PI decision surface — **RESOLVED 2026-07-13**

1. **Equilibrium MDA γ → (a) DECIDED.** γ stays **emergent** (=ΔP·R/2=0.16 mN/m); the direct MDA
   **Zhovmer 0.30 mN/m** is a **validation target** (VG-231-G1, KB-6.1.6). No γ runtime change.
2. **MDA ΔP → KEEP.** 40 Pa HeLa-proxy retained (labelled). No ΔP runtime change.
3. **η → 10.7 Pa·s** (Dessard Table-1 primary). **WIRED** (was 12.0).
4. **R → 7.5 µm** (Cognart direct suspended). **WIRED** (was 8.0). 231 suspended size ≈ MCF7.
5. **Ratified into Notion (DONE):** 4 new SourceEvidence (Zhovmer2021, FischerFriedrich2014, Cognart2020,
   Clark2013) + **6 KnowledgeClaims PI-ratified** (KB-6.2.1/6.1.4/6.1.5/6.2.2/6.3.1 + **6.1.6 tension
   taxonomy**) with Evidence links + **7 ValidationGates** (VG-231-G0…G6, not-started).

> **Runtime physics unchanged** for γ/ΔP/adhesion/NMII (decisions 1,2 = no runtime edit). Only the
> already-ratified R/η/E_nuc *value* wiring was updated to the PI picks (7.5 / 10.7 / 157.7).

---

## 6. Stale / misclassified statements (deliverable E) — corrections

**Applied (files this session authored):**
- `cell_type.py` mda_mb_231 header — **overclaim** "builds a real 231 cell rather than MCF7 mechanics"
  qualified: 231-specific for R/η/E_nuc (3 axes); ΔP/γ/contractility/adhesion/cadherin remain MCF7/generic.
- `ff_crawl_on_substrate.py` merge-comment — "R 10.5µm PROVISIONAL vol-derived" → **R 8.0µm** (Cognart).
- `_historical/MDA_MB_231_KB_REGISTRATION_DRAFT_2026-07-13.md` — §1 "spatial-org only / MCF7 mechanics" stale;
  GAP #3 (E_nuc) and GAP #4 (radius) **RESOLVED**; §7 "5 GAP" → **3 remaining** (γ, ΔP, contractility/
  adhesion/cadherin); η note "crawl driver still uses default" → now consumed from profile.

**Recommended (files owned elsewhere — flag, do not edit mid-audit):**
- `CANCER_CELLTYPE_PARAM_MAP.md` — split nucleus row: in-situ over-nucleus 157.7(231)/399(MCF7) Pa vs
  isolated 1–10 kPa; anchor lamin/E_nuc.
- `FF_S6_MIGRATION_AUTONOMOUS_PLAN` — "5 missing MDA params" → R/η/E_nuc wired; open = γ, contractility/
  adhesion, speed.
