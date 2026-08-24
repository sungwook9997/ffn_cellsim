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

# P2 ratification package — collagen / ERM / Hosseini γ (READY-TO-APPLY, awaiting PI final sign-off)

**Lane:** P2 evidence/SoT · **Date:** 2026-07-22 · **Builds on:** the 3 proposal cards in commit `e1e6e613`
(`P2_CARD{1,2,3}_*_2026-07-22.md`). **Status:** finalized against PRIMARY SOURCES this session; packaged as
apply-ready (a) gate-contract change / (b) Notion SoT registration draft / (c) code+KB fix patch spec.

**HARD constraints honoured.** No gate value or code was edited/committed — this doc is the sign-off package
only (no-gate-loosening). KB read layers (`kb.duckdb`, Obsidian vault) NOT touched — every registration is a
Notion Contract-Graph draft (SoT = Notion). `make kb-check` re-run this session: **[runs] gate OK · [params]
gate OK · no drift** (this package changes no committed KB artifact). Apply order after sign-off: register
Notion → refresh read layers → apply code patches → re-run `make kb-check`.

Primary sources pulled this session (not memory):
- **Collagen:** Yang & Kaufman 2009, *Biophys J* **96(4):1566–1585**, doi:10.1016/j.bpj.2008.10.063 (PMID
  19217873) — Table 1/2 G′N read directly from `paper_chunks` (the KB already holds this paper's text).
- **Hosseini γ:** Hosseini et al. 2020, *Adv Sci* 7(19):2001276, doi:10.1002/advs.202001276 (PMC7539203) —
  full text via PubMed/PMC (observable definition quoted verbatim below).
- **ERM:** Braunger et al. 2014, *J Biol Chem* 289(14):9833–9843, doi:10.1074/jbc.M113.530659 (PMC3975028) —
  abstract/full text via PubMed/PMC (single-molecule unbinding force quoted).

> Attribution: collagen + ERM + Hosseini full texts retrieved via PubMed/PMC. DOIs:
> [collagen](https://doi.org/10.1016/j.bpj.2008.10.063) · [Hosseini 2020](https://doi.org/10.1002/advs.202001276)
> · [Braunger 2014](https://doi.org/10.1074/jbc.M113.530659).

---

# CARD 3 — Hosseini γ (highest-confidence; do this first)

## Primary-source confirmations (this session)

1. **Observable definition (settles the double-count).** Hosseini 2020 Experimental Section, verbatim:
   *"we calculated a normalized change of cortical surface area (area strain) and an **effective cortical
   tension** from the output height and force signal … Cortical tension was estimated as the **time-averaged
   value … during an oscillation period** … assuming a shape of **minimal surface area at constant cell
   volume** … effective cortical tension is interpreted as **stress readout**."* → the observable is the
   **effective (total) surface tension in the Young-Laplace force balance**, not a cortex-minus-membrane
   quantity. **VERDICT — FIRST READING CONFIRMED:** `ΔP = 2·γ_eff / R`; **γ_mem is already inside γ_eff — do
   NOT re-add it** (the cortex-only reading `ΔP = 2(γ_c+γ_m)/R` does not apply to this datum).
2. **Value is figure-locked (confirmed).** The full text reports γ only graphically (Fig 2g) — there is **no
   inline "γ = 0.27 mN/m"** sentence. The figure-read stands: **median ≈ 0.27 mN/m, IQR 0.18–0.40, whiskers
   ≈ 0.05–0.62, n = 27**, MCF-7 suspended/rounded interphase, ±0.03 read precision (per
   `H7_HOSSEINI_GAMMA_CROSSCHECK_2026-06-07.md`, protocol + identity PubMed-verified).
3. **Protocol:** dynamic AFM parallel-plate confinement, wedged tipless cantilever, piezo oscillation
   amplitude 0.25 µm at **1 Hz**, minimal-surface Laplace tension; senior author Fischer-Friedrich (method's
   originator). Hosseini 2021 (BJ 120:3516, ~0.41 mN/m) is an **upper cross-check pending protocol
   reconciliation** (binning / cell-state differences vs the 2020 suspended-interphase read).

### Pressure mapping (R = 7.5 µm, first reading)
| γ_eff (mN/m) | ΔP = 2γ/R (Pa) |
|---|---|
| 0.18 (IQR low) | 48 |
| **0.27 (median)** | **72** |
| 0.40 (IQR high) | 107 |

## (a) Gate-contract change
- **NEW gate `VG-γ-MCF7-suspended`** (Slice-1 physiological initial-state γ): **band 0.18–0.40 mN/m, central
  0.27**, state = suspended/rounded interphase, T=37 °C, dynamic AFM confinement @1 Hz, observable = effective
  Laplace surface tension (`ΔP = 2γ_eff/R`, γ_mem folded in). Grades the **emergent** cortex tension at t0;
  never a runtime lumped substitute.
- **Reconcile with `VG-H3-KU35-cortex-tension` (0.35–0.65 mN/m):** re-scope 0.35–0.65 as a **modeling-target /
  different-state** band (`CORTICAL_TENSION_RECORD_2026-06-30`), NOT the Slice-1 suspended-interphase
  requirement. Overlap is only 0.35–0.40. **γ-floor conclusion unchanged** (0.27 vs 0.35 is ~1.3×; the
  measured-density floor is ~530× — re-scoping does not touch the density-limited verdict).
- **Pressure link (NOT a hand-swap):** Hosseini γ_eff feeds MCF7 `delta_p_hydrostatic_rest` via `ΔP=2γ_eff/R`
  (48–107 Pa, central 72). **40 Pa stays the HeLa diagnostic proxy** until a PI-authored pressure-consistency
  gate over {R, γ_cortex, γ_mem, ΔP_hyd} lands the MCF7 value.

## (b) Notion SoT registration draft — new KnowledgeClaim (paste into Notion Contract-Graph)
```
kb_id:            KB-6.1.7
title:            MCF-7 cortical tension (suspended interphase) — direct AFM value
subtopic:         Unit 6 Cell-type
status:           authoritative-candidate (PI sign-off pending)
value_range_si:   γ_eff = 0.27 mN/m (median; IQR 0.18–0.40; whiskers 0.05–0.62; n=27), MCF-7
                  suspended/rounded interphase, 37 °C, dynamic AFM parallel-plate confinement @1 Hz,
                  effective Laplace surface tension (time-averaged over an oscillation period, minimal-surface
                  at constant volume). Figure-read (Fig 2g), ±0.03. → ΔP = 2γ_eff/R (γ_mem INCLUDED; do NOT
                  re-add). Mitotic MCF-7 ≈ 0.75 (context, not the interphase gate).
unit:             mN/m
evidence:         SE181 Hosseini2020_AdvSci (Direct measurement)
citations_text:   Hosseini K, Taubenberger A, Werner C, Fischer-Friedrich E 2020 Adv Sci 7(19):2001276
                  (10.1002/advs.202001276, PMID 33042748, PMC7539203). Upper cross-check: Hosseini 2021
                  Biophys J 120:3516 (~0.41 mN/m, protocol-reconciliation pending). Def-tag: effective/total
                  Laplace tension → ΔP=2γ/R.
```
- Link edges: `KB-6.1.7 → SE181` (evidence); `KB-6.1.7 → VG-γ-MCF7-suspended` (gate); reconcile-link
  `VG-γ-MCF7-suspended ↔ VG-H3-KU35-cortex-tension`. Update **KB-6.1.3** note: "γ now numerically extracted →
  KB-6.1.7 (was figure-locked/scarce)."

## (c) Code/KB fix patch spec
- **No code value change** (no running config swaps to 72 Pa). `ac/fluid/params_i0b1.yaml`
  `delta_p_hydrostatic_rest` stays **40 Pa (HeLa proxy)**; add a one-line `expiry_trigger` note pointing at the
  new pressure-consistency gate (documentation only — no value edit pre-sign-off).
- KB: register KB-6.1.7 in Notion → `refresh.sh`. No `.duckdb`/vault hand-edit.

---

# CARD 1 — Collagen concentration-resolved gate

## Primary-source confirmation (this session) — Yang & Kaufman 2009 BJ 96:1566, Table 1 (37 °C) + Table 2 (32 °C)
Measured equilibrium storage modulus **G′N** of **pure PureCol collagen**, cone-plate (60 mm, 1°), oscillatory,
**37 °C, 1 Hz, strain amplitude 0.8 %** (linear viscoelastic):

| c (mg/mL) | **G′N @ 37 °C (Pa)** [Table 1] | G′N @ 32 °C (Pa) [Table 2] |
|---|---|---|
| 0.5 | 1.39 ± 0.01 | 1.77 ± 0.06 |
| 1.0 | 7.60 ± 0.07 | 11.5 ± 0.4 |
| **1.5** | **13.14 ± 0.04** | 29.5 ± 0.4 |

**→ 1.5 mg/mL production-gate value = G′N = 13.14 ± 0.04 Pa @ 37 °C (a DIRECT table value — NOT a figure read,
NOT the interpolated ~11 Pa).** c-scaling from Table 1: `n = log(13.14/1.39)/log(3) ≈ 2.05` (in the KB
n≈2.0–2.1 band). **Temperature is the dropped condition behind the code's 30–100 band:** at 32 °C the 1.5
mg/mL value is 29.5 Pa (≈ the code band's 30 floor) — the code silently mixed 32 °C-scale / 3 mg/mL numbers
onto a 1.5 mg/mL 37 °C label.

## (a) Gate-contract change — split into 4 conditional gates
- **`VG-ECM-G0(c)`** — small-strain G′ per concentration; **conditions T=37 °C, ω≈1 Hz (linear plateau),
  strain γ≤0.008 (LVE), 3D-bulk cone-plate rheology** (NOT AFM indentation, KB-1.14):
  | c (mg/mL) | G′ band (Pa) | source | role |
  |---|---|---|---|
  | 0.5 | 1.39 ± 0.01 | Yang-Kaufman 2009 BJ96 Table 1 | reference |
  | 1.0 | 7.60 ± 0.07 | " | reference |
  | **1.5** | **13.14 ± 0.04** (accept [10, 16] incl. source 2–10× spread) | " | **PRODUCTION GATE** |
  | 3.0 | 30–100 (anchor 55) | KB-1.32 / KB-1.V.2.1 / Licup 2015 (separate source, higher c) | 3 mg/mL cross-check |
  | 7.0 | ~342 | KB-1.V.2.1 | high-c cross-check |
- **`VG-ECM-c-scaling`** — `G′ ∝ c^n`, **n ∈ [1.9, 2.3]** (Table 1 37 °C n≈2.05; BJ97:2051 c^2.1; Licup 2015).
  T=37 °C, fit over 0.5–4 mg/mL. **Distinct from the code `conc_exponent` (mesh) — see (c).**
- **`VG-ECM-stiffening`** — differential modulus K rises **>10× above γ_c ~0.1–0.5**, K~σ, **γ_c
  concentration-independent** (KB-1.V.2.5); emergent-only, never tuned.
- **`VG-ECM-normal-stress`** — **negative N1** (semiflexible network pulls inward); sign is the acceptance
  criterion. **Blocking:** register a SourceEvidence row first (candidate Janmey et al. 2007 Nat Mater 6:48).
- **Dispositions:** `5–100 Pa` (`ff_ecm_validate.py:37`) → **exploration envelope** (1–3 mg/mL span, demoted);
  `30–100 Pa` (`ecm_library.py:111`) → **relabel 3 mg/mL** (NOT 1.5). Historical 11–15 Pa results →
  **superseded by the measured 13.14 Pa @ 37 °C** (they were consistent in magnitude; now grounded).

## (b) Notion SoT registration draft
- **Update KB-1.32** value_range_si to carry the Table-1 37 °C anchors explicitly:
  `G′N(37 °C,1 Hz,0.8%): 0.5→1.39, 1.0→7.60, 1.5→13.14 Pa (Yang-Kaufman 2009 BJ96:1566 Table 1); 3 mg/mL
  30–100 (Licup/KB-1.V.2.1); ξ~c^−0.5; G′∝c^2.1 (BJ97:2051).` Keep the c^2.1 cite on **BJ97:2051**.
- **Split the duplicated SourceEvidence** `YangKaufman2009_BiophysJ` (currently TWO DOI rows under one key)
  into two keyed rows:
  ```
  YangKaufman2009a_BJ96   doi 10.1016/j.bpj.2008.10.063  → G′N(c) direct rheology 0.5/1.0/1.5 mg/mL, 37/32 °C
  YangKaufmanLeone2009_BJ97 doi 10.1016/j.bpj.2009.07.035 → G′∝c^2.1 predicted from 2D confocal
  ```
  Re-anchor: 1.5 mg/mL G′N gate → **BJ96**; c^2.1 scaling → **BJ97**.
- New gates VG-ECM-{G0(c), c-scaling, stiffening, normal-stress} as ValidationGate rows per (a).

## (c) Code fix patch spec — `aleph/laws/ecm_library.py` `collagen_I` (lines 109–123) + `ff_ecm_validate.py`
```diff
# ecm_library.py — collagen_I ECMSpec
-        modulus_band_Pa=(30.0, 100.0), modulus_kind="G", poisson=0.2,
+        modulus_band_Pa=(10.0, 16.0), modulus_kind="G", poisson=0.2,   # 1.5 mg/mL 37°C: G'N=13.14 Pa
+                                                                        # (Yang-Kaufman 2009 BJ96 Table 1)
...
-        ref_conc=1.5, ref_conc_unit="mg/mL", ref_mesh_um=2.0, conc_exponent=0.5, strain_stiffens=True,
+        ref_conc=1.5, ref_conc_unit="mg/mL", ref_mesh_um=2.0, conc_exponent=0.5, strain_stiffens=True,
+        # NB conc_exponent is the MESH exponent ξ∝c^-0.5 (KB-1.7/1.V.2.2), NOT the modulus exponent
+        # n≈2.0-2.1 which EMERGES (do not conflate — the card's "second exponent conflict" is a FALSE ALARM).
...
-        sources=("KB-1.1", "KB-1.2", "KB-1.3", "KB-1.7", "KB-1.30", "KB-1.V.2.1",
-                 "Yang-Kaufman 2009 BiophysJ 10.1016/j.bpj.2008.10.063",
+        sources=("KB-1.1", "KB-1.2", "KB-1.3", "KB-1.7", "KB-1.32", "KB-1.V.2.1",
+                 "Yang-Kaufman 2009 BJ96:1566 10.1016/j.bpj.2008.10.063 (G'N(c) direct, Table 1)",
+                 "Yang-Leone-Kaufman 2009 BJ97:2051 10.1016/j.bpj.2009.07.035 (G'∝c^2.1 scaling)",
                  "Licup 2015 PNAS 10.1073/pnas.1504258112"),
-        notes="Emergent G'(c): 1→5, 3→55, 7→342 Pa (n~2.0-2.1). ⟨z⟩~3.2 sub-isostatic bending-dominated."),
+        notes="G'N(37°C,1Hz,0.8%): 0.5→1.39, 1.0→7.60, 1.5→13.14 Pa (Yang-Kaufman BJ96 Table1); emergent "
+              "n~2.05; 3mg/mL 30-100 (Licup). ⟨z⟩~3.2 sub-isostatic bending-dominated."),
```
```diff
# ff_ecm_validate.py — LIT["collagen_I"] (lines 36-38)
-    "collagen_I": {"kind": "G", "band": (5.0, 100.0), "ref_conc": 1.5,
-                   "note": "G'(1.5mg/mL)≈11 Pa (KB-1.V.2.1 c² from 5Pa@1mg/mL); scaling-anchor 30-100 (KB-1.30)"},
+    "collagen_I": {"kind": "G", "band": (10.0, 16.0), "ref_conc": 1.5,
+                   "note": "G'N(1.5mg/mL,37°C,1Hz,0.8%)=13.14±0.04 Pa (Yang-Kaufman 2009 BJ96 Table1); "
+                           "5-100=exploration envelope; 30-100=3mg/mL. n≈2.05 (VG-ECM-c-scaling)"},
```
- **`conc_exponent=0.5` — NO CHANGE** (it is the mesh exponent; the "second conflict" is a false alarm). The
  added comment prevents PI from "fixing" a non-bug.
- **Do NOT commit these diffs until PI signs off** the VG-ECM-G0(c) band values (no-gate-loosening).

---

# CARD 2 — MCF7 ERM density = production HOLD (capacity-gate-only)

## Primary-source confirmation (this session) — Braunger et al. 2014 (PMC3975028)
Abstract, verbatim: *"the measured **individual unbinding forces between ezrin and F-actin** are independent of
the activating parameters, in the range of **approximately 50 piconewtons**. However, the cumulative adhesion
energy greatly increases in the presence of PIP2, demonstrating that a **larger number of bonds** … has
formed."* → confirms:
- **Single ezrin–F-actin unbinding ≈ 50 pN** (activation-independent). Stiffness `k_erm = 4.6 pN/nm = 4600
  pN/µm`, `k_off ≈ 1.3 s⁻¹`, barrier width ~0.7 nm (KB-3.B1.6). Title: *"…Alters the **Number** of Attachment
  Sites…"* → the **density (bond number)** is the PIP2-modulated variable, NOT a fixed constant.
- **`11.4 pN` ≠ single-ERM.** `erm_rupture_force = 2π√(2κ_m(γ_mem+γ_MCA))` is the **continuum membrane
  tube/tether-extraction force scale** (KB-3.B1.4; f_t ~5–40 pN), already diagnostic in code. The resting
  diagnosis mislabel is corrected in `RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md`.
- Braunger being *about the number of attachment sites* is itself evidence that ρ_ERM is a physiological
  variable requiring an MCF7 datum + bound/active fraction — **reinforces the HOLD**, does not open it.

## (a) Gate-contract change
- **`VG-ERM-capacity` — NECESSARY-ONLY:** `ρ_ERM · f_active_per_ERM ≥ transmitted membrane→cortex traction`,
  with `f_active_per_ERM` bounded by the ~50 pN unbinding force × bound fraction. **The back-computed minimum
  ρ_ERM is a lower bound, NOT the physiological density** — do not adopt it.
- **Production ρ_ERM latch stays FALSE (HOLD).** `1/node` + zebrafish `600 µm⁻²` remain **diagnostic-only**.
- **Bundle contract (approve as ONE unit):** areal density ρ_ERM **+** bound/active fraction **+** Bell
  {k_on, k_off=1.3 s⁻¹, F0, capture radius}. Proxy priority: (1) MCF7 quantitative proteomics × cortical/
  membrane localization fraction; (2) other human epithelial direct areal density; (3) other-species embryo
  (diagnostic-only).

## (b) Notion SoT registration draft
- **No new ρ_ERM value** (that is the whole point — HOLD). Register only:
  - `VG-ERM-capacity` ValidationGate (necessary-only; formula + the 50 pN / 4600 pN·µm⁻¹ / k_off 1.3 s⁻¹
    Braunger inputs).
  - Update **KB-3.B1.6** note: add *"single-molecule unbinding ≈ 50 pN activation-independent (Braunger 2014
    abstract, PMC3975028); PIP2 sets bond NUMBER (density), not force — density is an MCF7 SoT GAP."*
  - Cross-link `KB-3.B1.6 (k_erm) ↔ KB-3.B1.4 (f_t tube-extraction)` with the "11.4 pN = tube-extraction, NOT
    single-ERM" disambiguation tag.

## (c) Code/KB fix patch spec
- **No code value change.** `k_erm=4600 pN/µm` (already correct, Braunger). `erm_rupture_force` stays a
  **diagnostic** in `assemble.py`/`preload_contract.py` — add a one-line comment: *"= continuum tube-extraction
  force scale (KB-3.B1.4), NOT the per-ERM unbinding force (~50 pN, KB-3.B1.6)."*
- KB: register `VG-ERM-capacity` + KB-3.B1.6 note update in Notion → refresh.

---

## Sign-off checklist for PI
- [ ] Card 3: register KB-6.1.7 (Hosseini γ 0.27, IQR 0.18–0.40) + `VG-γ-MCF7-suspended` (0.18–0.40); re-scope
      0.35–0.65 as modeling-target; confirm first-reading (γ_mem folded in); keep 40 Pa HeLa proxy.
- [ ] Card 1: adopt VG-ECM-G0(c) with **1.5 mg/mL = 13.14 Pa @ 37 °C**; split the Yang-Kaufman SE into
      BJ96/BJ97; approve the ecm_library.py / ff_ecm_validate.py diffs; register N1 SourceEvidence (Janmey 2007).
- [ ] Card 2: register `VG-ERM-capacity` (necessary-only); keep ρ_ERM latch FALSE; commission the ρ_ERM bundle.
- [ ] After sign-off, apply in order: Notion register → `refresh.sh` → code diffs → `make kb-check` green.

**Nothing above is committed as a gate value or code change. This package is the surface-to-PI artifact.**
