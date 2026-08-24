---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# MDA-MB-231 KB Registration Draft — 2026-07-13

**Status: DRAFT — pending PI sign-off before any Notion entry.**
Notion 8-DB Contract-Graph is the SoT; nothing here is written to duckdb/vault directly.
Precedent for this two-step (SE → Notion now, KnowledgeClaim/Contract/Gate → PI-authored):
`project-necrosis-kb-registration` (2026-07-06).

## 0. Provenance

Generated from workflow `wf_9804cba0-fd2` (13 agents): FF-engine parameter/harness audit +
KB-coverage scan + current-migration-state audit + 5 literature dimensions
(migration-kinematics · whole-cell-mechanics · traction-contractility · ecm-proteolysis-invasion ·
morphology-protrusion), each `find → adversarial citation-verify` pipelined.

- **37 verified SourceEvidence candidates** (deduped by DOI), all `verdict=keep`.
- **Citation integrity:** the verify pass **rejected/corrected 6 items** — see §5. This is exactly
  the class of error the 2026-06-02/06-08 audits target; caught pre-registration here.
- These are **literature validation targets / overlay**, NOT MCF7 runtime inputs
  (literature-first rule; `feedback-pi-exp-validation-target`).

## 1. Why this matters (engine-readiness one-liner)

The FF engine **already implements** mesenchymal invasion physics — `--mmp` (MMP secretion →
diffusive degradation halo → segment severing → invasion channel, KB-1.20), two-way plastic ECM
remodel, viscoelastic crosslink turnover, retrograde-flow clutch crawl.

> ⚠️ **SUPERSEDED (2026-07-13, later same day):** the paragraph below described the state BEFORE the
> R/η/E_nuc wiring. `--cell-type mda_mb_231` now sets **R=8.0µm, η=12.0 Pa·s, E_nuc=157.7 Pa** (3 of ~7
> mechanical axes) — a PARTIAL 231 build. Still MCF7/generic + PI-gated: turgor ΔP, γ (=ΔP·R/2), cortical
> density, contractility, adhesion, cadherin. See `MDA_MB_231_BASELINE_CONTRACT_2026-07-13.md`.

**Historical (pre-wiring):** `--cell-type mda`
was a *spatial-organization preset only* ([cell_type.py](../../ff/cell_type.py)); the
mechanical setpoints (R=7.5µm, η=65.9 Pa·s, turgor 40 Pa, cortical density 100/µm², E_nuc 399 Pa)
were hardcoded **MCF7** defaults in [cortex_assembly.py](../../ff/cortex_assembly.py) /
[units.py](../../ff/units.py) / [gamma_floor.py](../../ff/gamma_floor.py) /
[compartments.py](../../common/compartments.py), and the crawl driver read no YAML. So to test
"does 231 run correctly" we need (a) the **input parameters** below and (b) the **validation
targets** below. KB currently holds 231 only as an MCF7 *contrast* (η, whole-cell E, G); the
Parameter and ValidationGate tables have **zero** 231 rows.

---

## 2. KnowledgeClaim proposals (PI assigns final KB-x.y)

Proposed IDs are placeholders in the cancer namespace (KB-6.x). Each claim lists its datums,
units, cell-line context, and anchoring SourceEvidence (DOI). All numbers are **overlay /
validation targets**, not runtime constants.

### KB-6.2.1 (PROPOSED) — MDA-MB-231 migration kinematics

| Datum | Value | Context | Source (DOI) |
|---|---|---|---|
| 2D random speed | **0.48–0.57 µm/min** | plastic / no gradient | Holme 2023 `10.1038/s41598-023-50227-9`; Biswenger 2018 `10.1371/journal.pone.0203040` |
| 2D basal speed (regression) | 0.10–0.11 µm/min | collagen-coated µchannel | Islam & Resat 2017 `10.1039/c7mb00390k` |
| **3D collagen speed** | **0.39–0.52 µm/min** (range 7–114 µm/hr) | sandwiched / 3D collagen | Holme 2023; Kurniawan 2015 `10.3791/52735` |
| 3D collagen speed vs density | 0.4 (1 mg/mL) → 0.3 (3) → 0.07 (6 mg/mL) | mode switch mesench→amoeboid at 6 | Rodríguez-Cruz 2024 `10.3389/fcell.2024.1435708` |
| Directional persistence | **~0.4** (net/total path) | 3D collagen | Kurniawan 2015 |
| MSD diffusion coeff (in-spheroid) | 0.031→0.047 µm²/min (control→compressed) | 3D spheroid | Pandey 2024 `10.1088/1478-3975/ad3ac5` |
| Anomalous exponent α | 2D superdiffusive ~1.3–1.4; 3D subdiffusive <1 | — | Luzhansky 2018 `10.1063/1.5019196` |
| Durotaxis | **+ (toward stiff)**, reversible by talin-1/2 KD | PAA 0.5–22 kPa | Isomursu/Odde 2022 `10.1038/s41563-022-01294-2` |
| Stiffness→speed | monotonic ↑ over 1–20 kPa | PAA 2D | Lin 2018 `10.1021/acsbiomaterials.7b00835` |
| Pore-size optimum | max at pore≈cell size **with stiff fibers**; soft fibers abolish | 3D collagen | Geiger 2019 `10.1371/journal.pone.0225215` |
| Confined fast subpop | G1 > 0.434 µm/min | µchannel | Zhang & Mak 2021 `10.1038/s41598-021-85640-5` |
| Stromal modulation | ASC/adipose spheroids ↑ motile fraction, persistence, speed | 3D | Horder 2024 `10.1088/1758-5090/ad57f7` |
| EGF mesench shift | ↑persistence, motile threshold >0.2 µm/min | 3D | Geum 2016 `10.1140/epjp/i2016-16008-8` |
| Confinement+shear mode | moderate confinement/shear ↑migration; tight→amoeboid | µchannel | Chen 2025 `10.1039/d5lc00563a` |

### KB-6.1.4 (PROPOSED) — MDA-MB-231 traction & contractility (extends KB-6.1.3)

| Datum | Value | Context | Source (DOI) |
|---|---|---|---|
| **3D contractility (net inward)** | **47.6 ± 3.4 nN** | 3D collagen, TFM (Steinwachs) | Steinwachs 2016 `10.1038/nmeth.3685` |
| 3D contractility (control) | **70.2 ± 3.5 nN** | 3D collagen | Cóndor 2019 `10.1016/j.bpj.2019.02.029` |
| **2D total traction \|F\|** | **~300 nN** (90/305/375 @ 1/5/10 kPa) | 2D PAA, stiffness/ligand sweep | Kraning-Rush 2012 `10.1371/journal.pone.0032572` |
| 2D avg traction stress | a few hundred Pa (peak ~455–1280 Pa) | soft 2D | Lin 2018; Kraning-Rush 2012 |
| Strain energy (3D) | 0.08 (median) / 0.51 ± 0.06 (avg) / up to 3.7 pJ | 3D collagen | Cóndor 2019; Wiener 2021 `10.1016/j.jbiomech.2021.110759`; Koch 2012 `10.1371/journal.pone.0033476` |
| Force polarity (dipole frac) | 0.47–0.53 | 3D collagen | Steinwachs 2016; Cóndor 2019 |

### KB-6.1.5 (PROPOSED) — MDA-MB-231 whole-cell & cytoplasm mechanics (extends KB-6.1.1 / KB-PIV-7)

| Datum | Value | Context | Source (DOI) |
|---|---|---|---|
| Whole-cell E (AFM) | **224 Pa** (colloidal); 262/105/83 Pa (E_inst/E_inf/E_1) | AFM force-spec + stress-relax | Zbiral 2023 `10.3390/ijms241512208` |
| Whole-cell E (micropipette) | 206.2 ± 23.1 Pa (MCF10A ref 441 Pa) | micropipette aspiration | Lee & Liu 2015 `10.1039/c4lc01218f` |
| Whole-cell E (AFM, over nucleus) | 9.0 ± 1.53 kPa (MCF7 ref 5.0 kPa) | Sneddon over nucleus (method-dependent) | Kwon 2020 `10.7150/jca.45897` |
| Cytoplasm η | **10.7 ± 5.4 Pa·s** (~5× < MCF7 65.9) | active microrheology | Dessard 2024 `10.1039/d4na00003j` (KB-PIV-7) |
| Cytoplasm shear G | 38.6 ± 5.8 Pa | active microrheology | Dessard 2024 |
| Stress-relax times / η₁,η₂ | τ₁=0.18, τ₂=4.21 s; η₁=13, η₂=340 Pa·s | AFM two-timescale | Zbiral 2023 |
| Cortical stiffness trend | significantly SOFTER than MCF7/HBL-100 | peak-force AFM + OT | Coceano 2016 `10.1088/0957-4484/27/6/065102` |

> Note: whole-cell modulus **ordering is method-dependent** (colloidal/micropipette ~200 Pa vs
> Sneddon-over-nucleus ~9 kPa). Do NOT encode a single global 231 modulus — per-compartment,
> per-method (`project-mcf7-parameter-collection` governing rule).

### KB-6.2.2 (PROPOSED) — MDA-MB-231 morphology & nuclear deformation

| Datum | Value | Context | Source (DOI) |
|---|---|---|---|
| **Aspect ratio (2D)** | **4.42 ± 0.44** (circularity 0.37) | 2D spread | Liew 2024 `10.1007/s12195-024-00811-4` |
| Spread area (2D) | 598.3 ± 36.7 µm² | 2D | Liew 2024 |
| 3D axes / aspect | major 49–66 µm, minor 14–22 µm, AR ~2.2–4.7 | 3D collagen (density-dep) | Rodríguez-Cruz 2024 |
| Cell volume (3D) | 4120–5823 µm³ | 3D collagen | Rodríguez-Cruz 2024 |
| Nuclear width (unconfined) | 11.62 µm | — | Stöberl 2024 `10.1126/sciadv.adm9195` |
| **Confined speed peak** | biphasic, **peaks @ channel ~11–12 µm ≈ nuclear diameter** | µchannel | Stöberl 2024 |
| Nuclear ε / volume in confinement | ε→1.4 @12µm, →0.5 @4µm; ΔV up to −11% | µchannel | Stöberl 2024; Wang 2026 `10.7150/thno.119211` |
| Nuclear non-convexity (spreading) | INCREASES (vs MCF10A ↓), plastic (no relax after dissection) | 2D spreading | Tocco 2017 `10.1002/jcp.26031` |
| Leading-edge protrusions | 2-protrusion phenotype after constriction, persists 60–240 min | µchannel | Zhang & Mak 2021 |

### KB-6.3.1 (PROPOSED) — MDA-MB-231 proteolytic invasion (231-specific; extends generic KB-1.20)

| Datum | Value | Context | Source (DOI) |
|---|---|---|---|
| **MT1-MMP absolute requirement** | invasion of cross-linked 3D collagen REQUIRES MT1-MMP | native collagen | Sabeh 2009 `10.1083/jcb.200807195` |
| **Physical pore threshold** | migration ~linear in pore size; arrest at **~7 µm²** (≈10% nuclear XS), absolute ~5 µm² | 3D, MMP-blocked | Wolf 2013 `10.1083/jcb.201210152` (SE21, already in KB) |
| MAT switch (rate-preserving) | near-total protease block → spindle→round, β1/MT1-MMP decluster, no tracks, squeezes gaps | 3D collagen | Wolf 2003 `10.1083/jcb.200209006` |
| Invadopodial protease | MT1-MMP is key; degradation separable from invadopodium assembly | gelatin | Artym 2006 `10.1158/0008-5472.CAN-05-2177` |
| Collagen-I cleavage (reC1M) | MMP14/collagen-I proteolysis drives 231 invasion | invadopodia | Meng 2024 `10.1242/jcs.261608` |
| Secreted gelatinases | pro-MMP-2 + MMP-9 present (zymography); high vs MCF7 | conditioned medium | Bachmeier 2001 (PMID 11911253, no DOI) |
| Matrigel amoeboid mode | rounded, rear-contractility, MMP-independent, ~3.2 µm/h | 3D Matrigel | Poincloux 2011 `10.1073/pnas.1010396108` |
| β1-integrin clustering | co-clusters with MT1-MMP at fibers; lost in amoeboid switch | 3D collagen | Wolf 2003 |

---

## 3. Parameter proposals

**INPUT (runtime knobs) — one already wired, rest are the GAPs in §4:**

| Param | 231 value | Status | Engine location (MCF7 now) |
|---|---|---|---|
| `eta_cytoplasm` | 12.0 Pa·s | ✅ registered (KB-PIV-7), wired in phase1_h10.yaml | units.py:56 (65.9) — crawl driver still uses default |
| whole-cell E (Hertz band floor) | 224 Pa | ✅ used as band floor | ff_hertz_validation.py:25 |

**VALIDATION-TARGET bands (overlay, feed proposed ValidationGates — PI-authored):**

| Proposed gate | Band | Anchor |
|---|---|---|
| VG-231-mig-3D (PROPOSED) | 3D collagen speed 0.39–0.52 µm/min | KB-6.2.1 |
| VG-231-persistence (PROPOSED) | ~0.4 | KB-6.2.1 |
| VG-231-traction-3D (PROPOSED) | 47–70 nN | KB-6.1.4 |
| VG-231-aspect (PROPOSED) | AR ~4.4 (2D) / 2.2–4.7 (3D) | KB-6.2.2 |
| VG-231-confined-peak (PROPOSED) | speed peak @ channel ≈ nuclear dia | KB-6.2.2 |
| VG-231-MMP-pore (PROPOSED) | arrest <~7 µm² w/o MT1-MMP | KB-6.3.1 |

> ValidationGate / ModelContract rows are **PI-authored, never auto-created** (CLAUDE.md). Listed
> as PROPOSED for PI to ratify.

---

## 4. PI-gated GAPs — 2 of 5 now RESOLVED; **3 remaining**

> **UPDATE 2026-07-13:** GAP #3 (E_nuc) and GAP #4 (radius) are **RESOLVED** by the sourcing dig +
> primary-source re-audit (wf_35f1b9e4-d78). Also GAP #1 was reframed: a direct MDA cortical tension
> DOES exist (Zhovmer 2021 ~0.30 mN/m, seconds-AFM) — see `MDA_MB_231_BASELINE_CONTRACT_2026-07-13.md` §2.

1. **Cortical tension γ** *(remaining, runtime — PI-gated)* — γ is **emergent** (=ΔP·R/2 = 0.16 mN/m).
   A DIRECT MDA target now exists: **Zhovmer 2021 ~0.30 mN/m** (seconds AFM, rounded) — γ sits ~1.9× below
   it. Options: (a) keep emergent, validate against 0.30 *(recommended)*; (b) raise ΔP / add myosin. **PI.**
2. **Turgor / osmotic ΔP** *(remaining, runtime — PI-gated)* — no 231-specific value; **40 Pa HeLa proxy**
   (Fischer-Friedrich 2014, labelled). **PI decision.**
3. ~~**Nucleus E_nuc**~~ **RESOLVED** — **157.7±78.6 Pa** (Fischer, Hayn & Mierke 2020, in-situ over-nucleus;
   231 ~2.5× softer than MCF7 399; band 200→150 PI-ratified). **WIRED.**
4. ~~**Suspended radius**~~ **RESOLVED** — **R≈8.0 µm** (Cognart 2020 direct suspended dia ~15µm, band 7.5–10;
   the earlier vol-derived 10.5 was protrusion-inflated). **WIRED.** (231 is a small cell, ≈ MCF7.)
5. **Contractility multiplier · adhesion (n_fa / integrin identity) · cadherin identity (E→N,
   vimentin)** — traction magnitudes exist (§KB-6.1.4) but their mapping to `contractility_mult` /
   `n_fa` / catch-bond params is a modeling choice. **PI decision** (mechanistic rule: E→N-cadherin
   catch-bond, not a lumped switch — `feedback-junction-switch-fine-grained`).

---

## 5. Exclusions & citation corrections (integrity — do NOT register as-is)

| Item | Problem | Resolution |
|---|---|---|
| "Sozzi et al." (3D volume) | **fabricated author** | corrected → Rodríguez-Cruz 2024 `10.3389/fcell.2024.1435708` |
| "Elson et al." kPa traction (99.4 kPa) | **wrong units (kPa→Pa) + misattributed** | REJECT the kPa figure; real paper Liew 2024 `10.1007/s12195-024-00811-4` |
| "De Belly, Tsujita" membrane tension | wrong author + article no. | corrected → Tsujita 2021 `10.1038/s41467-021-26156-4` (SE213, already in KB); numeric value not extractable this pass |
| Jerrell & Parekh degradation-vs-rigidity | **cell line = SCC-61, NOT 231** | REJECT for 231 (keep as SCC reference only) |
| Mendoza 2013 FA turnover DOI | wrong DOI (`10.1242/jcs.103788`) | correct DOI `10.1242/jcs.119727`; FA-lifetime still UNDER-SOURCED — dedicated primary needed |
| SE231 (Yang/Anseth), VincentEngler (MSC) | matched search terms, **not 231 cancer** | exclude (false positives) |
| KB-PIV-7 in-row cite "Hu 2024" | wrong attribution | reconcile → Dessard/Manneville/Berret 2024 `10.1039/d4na00003j` |
| Lamellipodium width/thickness | general lit, not 231 | keep as generic reference only |
| Per-clutch / FA traction (231-specific) | **not found** | GAP — do not fabricate |

## 6. Existing-KB reconciliation

Already registered (do not duplicate; extend): **KB-PIV-7** (η 12.0, G 38.6), **KB-6.1.1**
(whole-cell E 224), **KB-6.1.3** (traction, contested), **KB-3.B1.5** (membrane ~2× lower).
Already-registered SE: Wolf2013 (SE21), Dessard2024, Zbiral2023, Kraning-Rush2012, Tsujita2021
(SE213), Liu2015 (SE242), Friedl2011 (SE108). Of the 37, these ~7 are in KB; ~30 are new SE.

## 7. Next step (PI)

- [x] **DONE 2026-07-13** — 29 new SourceEvidence rows pushed to Notion (UID **378–406**, all
  `Anchor Status = verified 2026-07-13`). 6 already-present DOIs skipped; 6 fabrication/misattribution
  items (§5) excluded. Jacquemet FiloQuant dropped (method-only, no 231 datum). Claims relation left
  empty pending KnowledgeClaim authoring.
- [ ] PI ratify the **5 proposed KnowledgeClaims** (assign KB-6.2.1/6.1.4/6.1.5/6.2.2/6.3.1) → then
  link each SE's `Claims` relation.
- [x] **2 of 5 GAPs resolved** (E_nuc 157.7 Pa, R 8.0 µm — WIRED). PI supply/approve the **3 remaining**
  (§4): cortical tension γ handling, turgor ΔP, contractility/adhesion/cadherin. **Runtime γ/ΔP/adhesion/
  NMII HELD for PI (2026-07-13 directive).** Full decision surface: `MDA_MB_231_BASELINE_CONTRACT_2026-07-13.md` §5.
- [ ] PI ratify the **6 proposed ValidationGates** (§3).
- After KnowledgeClaim/Gate entry: `python outputs/tag_kb/harvest_ops.py` + `refresh.sh`; `make kb-check`.
