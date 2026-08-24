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

# W_cs cell–substrate adhesion energy — KB registration draft (2026-07-13)

**Status: DRAFT — SourceEvidence registered to Notion; KnowledgeClaim PROPOSED, pending PI sign-off.**
Notion 8-DB Contract-Graph is the SoT; nothing here is written to duckdb/vault directly.
Two-step precedent (SE → Notion now, KnowledgeClaim → PI-authored):
`project-necrosis-kb-registration`, `MDA_MB_231_KB_REGISTRATION_DRAFT_2026-07-13`.

## 0. Why this registration (the trigger)

The DCM spreading session (2026-07-13) needs the **cell–substrate adhesion energy density
`W_cs`** grounded, because it sets the Young–Dupré / Douezan wetting balance
**S = W_cs − 2γ** that decides whether the cell spreads at all (G2 gate). The 2026-07-13
KB-grounding pass found `W_cs` is **NOT registered in the KB** (0 rows in knowledge_claim /
parameter / source_evidence for "adhesion energy" / "W_cs" / "cell-substrate") — it is used
throughout the code as a "measured MCF7" number with **broken provenance** (see §1). Grounding
it is a hard-rule requirement (no magic numbers; citation integrity).

## 1. The provenance problem (what triggered the audit)

`W_cs = 2.85e-3 J/m²` is used across the DCM code, attributed on-disk to **Gil-Redondo et al.
2023** (Microsc Res Tech 86(9):1069, doi 10.1002/jemt.24368). Two problems:

1. **Misattribution — CONFIRMED (PubMed, 2026-07-13).** Gil-Redondo et al. 2023
   (PMID 37345422, PMC10952526, doi 10.1002/jemt.24368) is a **methodology REVIEW/primer**,
   "Measuring (biological) materials mechanics with atomic force microscopy. 5. Traction force
   microscopy," using MCF10-A / MCF-7 / MDA-MB-231 as practical examples for *traction* forces
   (its headline finding: higher-metastatic MDA-MB-231 has LOWER traction). It reports traction,
   not an adhesion-energy density, and is a Review (not a primary adhesion-energy measurement).
   So `W_cs = 2.85e-3 J/m²` has **no valid source** from this citation.
2. **Possible 1000× unit slip.** Gil-Redondo's own strain-energy ÷ spread-area =
   0.0052 pJ / 1822 µm² = **2.85e-6 J/m²** — the same leading digits "2.85" but 1000× smaller.
   Either (a) a mJ↔µJ slip, or (b) the digit match is coincidence and 2.85e-3 came from elsewhere.
   NB: strain-energy density ≠ adhesion-energy density (different physical quantities) — must
   not be conflated. Given the project's 10⁶× DCM unit-slip history (`reference-dcm-meters-ff-microns`),
   this warranted a check before the number drives the spreading verdict.

## 2. Literature resolution — the correct value + primary source (2026-07-13 verify pass)

**Correct value: W_cs ≈ 1×10⁻⁴ – 1×10⁻³ J/m² (0.1–1 mJ/m²)** — a derived RANGE (effective
interfacial adhesion entering S = W_cs − 2γ), NOT a single measured "MCF7-on-fibronectin" datum
(no paper reports that number directly). Best anchor = the integrin first-principles derivation
(grid-invariant, satisfies the Magic-Number Block), corroborated by the tissue-wetting literature.

**First-principles (integrin) cross-check** — W_cs ≈ ρ_integrin · E_bond · f_engaged
(E_bond = 25 kT @ 310 K = 1.07×10⁻¹⁹ J):

| ρ (/µm²) | f_engaged = 0.1 | f_engaged = 1.0 |
|---|---|---|
| 100 | 1.1×10⁻⁶ | 1.1×10⁻⁵ |
| 1000 | 1.1×10⁻⁵ | 1.1×10⁻⁴ |
| 3000 (mature FA) | 3.2×10⁻⁵ | 3.2×10⁻⁴ |

Bare integrin bonds → 10⁻⁶–3×10⁻⁴ J/m²; reaching ~10⁻³ (mJ/m²) needs FA enrichment +
cortical-tension amplification. **Supports the mJ/m² (10⁻⁴–10⁻³) scale; excludes the 10⁻⁶ scale.**

**Primary sources (all DOIs PubMed-verified 2026-07-13):**

| Role | Source | DOI | SE row |
|---|---|---|---|
| Framework S = W_cs − W_cc (W_cs normalized 0–1 only) | Douezan et al. 2011, PNAS 108:7315 | 10.1073/pnas.1018057108 | **SE157** ✓ |
| Magnitude anchor (tissue surface tension ≈ effective adhesion, 1–45 mN/m; breast lines) | Gonzalez-Rodriguez et al. 2012, Science 338:910 | 10.1126/science.1226418 | **SE159** ✓ |
| Raw vs effective adhesion (raw ~10⁻⁷ "too small"; effective = tissue surface tension ~10⁻³ N/m) | Winklbauer 2015, J Cell Sci 128:3687 | 10.1242/jcs.174623 | **NEW → register** |

**Unit verdict:** the "2.85" is a SLIP, not coincidence: 2.85×10⁻³ = 998.6× the Gil-Redondo
strain-energy density (2.85×10⁻⁶ J/m² = 0.0052 pJ ÷ 1822 µm²) — a clean 1000× (mJ↔µJ) magnitude
error STACKED on a quantity error (strain-energy density ≠ adhesion-energy density). The code
comment confirms the derivation: `dcm_two_stage_production.py:427` "lit ~2.85e-3 from Gil-Redondo
strain energy / spread area." Ironically the ×1000 error lands 2.85 mJ/m² *inside* the plausible
band, so the magnitude is defensible (slightly high vs the 0.1–1 mJ/m² centre) but the derivation
is doubly invalid — **do not keep 2.85e-3 with the Gil-Redondo citation.**

**Do NOT bump W_cs to force spreading (outcome-tuning ban).** At the grounded W_cs (~0.5e-3) and
even the slow-band γ (5e-4), the Douezan S = W_cs − 2γ is marginal-to-negative — passive wetting
does NOT by itself make a single cell pancake. That is expected and consistent with this session's
finding: **T-2 active area generation is the real spreading driver**; W_cs sets the passive floor
and matters in the aggregate regime (γ_tissue ~ mN/m, Gonzalez-Rodriguez). Confidence: **Medium**
(range + effective/model-dependent quantity; no direct MCF7-on-FN datum).

## 3. Code state audit (the inconsistency to reconcile)

`W_cs` / `w_cs` is used with **two different values** and **conflates cell–cell with
cell–substrate** adhesion:

| Location | Symbol | Value | Meaning as-labelled | Correct meaning |
|---|---|---|---|---|
| `dcm/geometry.py:113` `ResolvedDCM.W_cs_Jm2` | W_cs | **0.5e-3** | cell–**substrate** (integrin) | cell–substrate ✓ |
| `dcm/geometry.py:114` `ResolvedDCM.W_cc_Jm2` | W_cc | **0.2e-3** | cell–**cell** (cadherin) | cell–cell ✓ |
| `dcm/dcm_warp_decohesion.py:266` `w_cs_jm2=` | w_cs | **2.85e-3** | cell–substrate (wetting/polarized) | cell–substrate |
| `dcm/dcm_interfacial_tension_warp.py:22` | w_cs | **2.85e-3** | "MCF7 **cell–cell** adhesion energy" (DAH foam) | cell–**cell** |
| `dcm/nondim.py:48` | w_cs | **2.85e-3** | "lit MCF7, Gil-Redondo 2023" | (ambiguous) |
| `dcm/dcm_contact_conservative_warp.py:96` | w_cs | (arg) | adh = 4·w_cs/c_adh² (contact) | cell–cell contact |

**Finding:** `2.85e-3` is used as BOTH the cell–cell (interfacial DAH, contact) AND the
cell–substrate (wetting, polarized) adhesion energy — these are DIFFERENT physical energies
(cadherin trans-dimer vs integrin–ECM). Meanwhile `ResolvedDCM` carries a *separate*,
smaller pair (0.5e-3 / 0.2e-3). The registration must (a) pin `W_cs` (cell–substrate) and
`W_cc` (cell–cell) as distinct grounded claims, and (b) flag the code to converge on the
registered values (a follow-up code-reconciliation, PI-gated — do NOT retune to pass a gate).

## 4. KnowledgeClaim proposals (PI assigns final KB-x.y; Unit 2 ECM-Cell)

### KB-2.x (PROPOSED) — cell–substrate adhesion energy density W_cs
- **Value Range SI**: _[from §2]_ J/m²
- **Unit**: Unit 2 ECM-Cell
- **Subtopic**: cell–substrate wetting adhesion (Young–Dupré / Douezan S = W_cs − 2γ)
- **Confidence**: _[from §2]_
- **Assumptions**: integrin–ECM areal adhesion energy = integrin density × bond energy ×
  engaged fraction; sets the wetting number W_cs/γ and the spreading coefficient S = W_cs − 2γ.
- **Evidence**: _[SE rows from §2]_

### KB-2.y (PROPOSED, if distinct) — cell–cell adhesion energy density W_cc
- Only if the literature pass separates cadherin cell–cell from integrin cell–substrate; the
  code currently reuses 2.85e-3 for cell–cell (interfacial DAH), which may need its own anchor
  (cadherin trans-dimer energy, KU-4.x namespace). Flagged for PI.

## 5. Parameter proposals

| Param | value | status | engine location |
|---|---|---|---|
| `W_cs` (cell–substrate) | _[§2]_ | PROPOSED | `dcm_substrate_warp` well/wetting/adhesion-energy; decohesion `w_cs_jm2` |
| `W_cc` (cell–cell) | _[§2 / existing]_ | PROPOSED | `dcm_interfacial_tension_warp`; `dcm_contact_conservative_warp` |

## 6. SourceEvidence — Notion registration

**Already registered (link as Evidence, do NOT duplicate):**
- **SE157** Douezan2011_PNAS (10.1073/pnas.1018057108) — framework S = W_cs − W_cc.
- **SE159** GonzalezRodriguez2012_Science (10.1126/science.1226418) — tissue-tension magnitude anchor.

**Winklbauer 2015 — ALREADY registered as SE311** (dedup correction, 2026-07-13):
- **SE311** Winklbauer2015_JournalOfCellScience (10.1242/jcs.174623, PMID 26471994) — was
  auto-ingested 2026-06-04 from the references TAG corpus (paper "Cell adhesion strength from
  cortical tension"). Its `Short Source` was empty; **enriched** this pass with the W_cs
  magnitude-anchor note, and linked to KB-2.19 (Evidence↔Claims, bidirectional).
- ⚠️ **Integrity note:** I first created a *duplicate* SE (Winklbauer2015_JCellSci, same DOI)
  before finding SE311 — my earlier dedup regex omitted "winklbauer" (the source wasn't known
  until the verify pass returned). The duplicate was **removed from the SourceEvidence DB**
  (moved to workspace-orphan) and KB-2.19's Evidence re-pointed to the canonical SE311. Net: no
  duplicate remains in the KB; dedup must key on DOI, not just name.

**EXCLUSION (citation integrity — never register as W_cs source):**
- **Gil-Redondo 2023** (10.1002/jemt.24368, PMID 37345422) — a TFM methodology REVIEW; reports
  MCF7 traction 102 nN + strain-energy density 2.85×10⁻⁶ J/m², **no adhesion-energy density**.
  The W_cs = 2.85e-3 attribution is a ×1000 unit slip on the strain-energy density + a
  quantity error. Documented here; the on-disk citation must be removed in the code reconciliation.

## 7. Next steps

- [x] **DONE 2026-07-13** — §2 filled from the literature-verify pass (value 0.1–1 mJ/m²,
      unit verdict = ×1000 slip + quantity error, provenance = Gil-Redondo misattribution).
- [x] **DONE 2026-07-13** — Registered to Notion:
  - **SourceEvidence** Winklbauer2015_JCellSci (10.1242/jcs.174623) — page
    `39c120da-ec5d-81ec-9ad8-c272ad275431`, Anchor Status verified 2026-07-13.
  - **KnowledgeClaim KB-2.19** (PROPOSED, **Status = draft**) "Cell-substrate adhesion energy
    density W_cs", Unit 2 ECM-Cell, Value 1e-4–1e-3 J/m², Evidence → SE157 (Douezan) + SE159
    (Gonzalez-Rodriguez) + Winklbauer — page `39c120da-ec5d-8149-ac2d-edbda9462c5c`.
  - Gil-Redondo 2023 **NOT registered** as W_cs source (documented exclusion, §6).
- [ ] **PI ratify KB-2.19** — confirm/assign final KB ID, flip Status draft → PI-ratified.
- [ ] Code reconciliation (PI-gated): converge the code's W_cs / W_cc on the registered
      value + remove the Gil-Redondo citation at the 5 on-disk sites (`nondim.py:48`,
      `dcm_substrate_warp.py:763`, `dcm_interfacial_tension_warp.py:22`,
      `dcm_warp_decohesion.py:268/303/2062`, `dcm_two_stage_production.py:427`). Do NOT
      retune to outcomes (`feedback-no-param-tuning-to-outcome`); W_cs is the passive floor,
      T-2 is the spreading driver.
- [ ] `bash outputs/tag_kb/refresh.sh` (notion→duckdb pull) so KB-2.19 is queryable + Obsidian
      mirrors it; `make kb-check`.
