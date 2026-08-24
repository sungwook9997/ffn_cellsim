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

# P2 · Card 3 — MCF7 γ_cortex = re-register Hosseini 2020 direct MCF-7 measurement (PROPOSAL, PI sign-off)

**Lane:** P2 evidence/SoT · **Date:** 2026-07-22 · **Status:** PROPOSED KnowledgeClaim + gate reconciliation.
No value swapped into a running config. `make kb-check` unaffected. **SoT rule:** Notion Contract-Graph only.

Executes `AC_DECISION_CARDS_2026-07-22.md` Card 3. Hosseini et al. 2020 **is** a direct MCF-7 cortical-tension
measurement, but the KB (KB-6.1.3) only **cites** it as "figure-locked/scarce" with **no numeric γ extracted**.
This is a SoT gap ("we have the datum but never registered its value"), not "no data".

---

## 1. The datum (verified this session; identity + protocol confirmed 2026-06-07)

> **Hosseini K, Taubenberger A, Werner C, Fischer-Friedrich E (2020).** "EMT-Induced Cell-Mechanical Changes
> Enhance Mitotic Rounding Strength." *Adv Sci* 7(19):2001276. DOI 10.1002/advs.202001276 · PMID 33042748 ·
> PMC7539203. SourceEvidence **SE181** (`Hosseini2020_AdvSci`, source_type = Direct measurement).

- **MCF-7, suspended/rounded interphase**, dynamic AFM parallel-plate confinement (wedged tipless cantilever,
  0.2–0.8 N/m), **n = 27**. Senior author Fischer-Friedrich originated the method.
- **Cortical tension = time-averaged effective cortical surface tension over an oscillation period**, measured by
  Laplace's law **γ = F / [A(1/R₁+1/R₂)]** (confirmed from PMC full text, `H7_HOSSEINI_GAMMA_CROSSCHECK_2026-06-07.md`).
- **Value (Fig 2g, figure-read, ±0.03):** median ≈ **0.27 mN/m**, **IQR 0.18–0.40 mN/m**, whiskers ≈ 0.05–0.62.
- Mitotic MCF-7 (same assay) ≈ 0.75 mN/m (rounding-strength context; not the interphase gate).
- **Intrinsically figure-bound:** the paper reports γ only graphically — no inline "γ = 0.27" sentence. The read is
  the best available value; the "figure-read ±0.03" caveat stays attached (a precision note, **not** a verification gap).

**Upper cross-check — Hosseini 2021** (Biophys J 120:3516, PMID 34022239, SE180): whole-cell stiffness grades
invasiveness; the ~**0.41 mN/m** MCF-7 figure is an **upper cross-check pending protocol reconciliation** (binning /
cell-state differences vs the 2020 suspended-interphase read).

---

## 2. Definition / double-count resolution — **FIRST reading confirmed**

Card 3's warning: if Hosseini's observable is *effective* cortical/surface tension, re-adding γ_mem double-counts.

**Resolved to the first reading.** The 2026-06-07 cross-check established Hosseini's observable is the **effective
cortical surface tension in the Young-Laplace balance** — γ is defined *from* the pressure/curvature relation
`γ = F/[A(1/R₁+1/R₂)]`, so it is the **total effective tension** resisting ΔP over the whole rounded surface (cortex
+ in-plane membrane, lumped into the Laplace γ). Therefore:

> **Hosseini γ = total effective tension → use ΔP = 2·γ_eff / R. Do NOT add γ_mem separately (it is already inside γ_eff).**

The cortex-only reading `ΔP = 2(γ_cortex + γ_mem)/R` does **not** apply to this datum.

### Pressure mapping at R = 7.5 µm
| γ_eff (mN/m) | ΔP = 2γ/R (Pa) |
|---|---|
| 0.18 (IQR low) | **48** |
| **0.27 (median)** | **72** |
| 0.40 (IQR high) | **107** |

So the Hosseini-consistent resting **ΔP ≈ 48–107 Pa, central ~72 Pa**.

### Consequence for the current 40 Pa
Current `delta_p_hydrostatic_rest = 40 Pa` ⇒ γ_total = 40·7.5e-6/2 = **0.15 mN/m**, **below** the Hosseini IQR
(0.18–0.40). Per Card 3, **40 Pa stays a HeLa diagnostic proxy only**
([params_i0b1.yaml](../../ac/fluid/params_i0b1.yaml)) until a PI-authored pressure-consistency gate lands the MCF7
value. **Do NOT hand-swap 40 → ~72 Pa;** the pressure-consistency gate over {R, γ_cortex, γ_mem, ΔP_hyd} decides it
(per the yaml `mcf7_consistency_gate`).

---

## 3. Reconciliation with the existing KU-3.5 target (0.35–0.65 mN/m)

- The **0.35–0.65 mN/m** band is a **modeling-contract target** (`CORTICAL_TENSION_RECORD_2026-06-30`,
  VG-H3-KU35-cortex-tension), **not** a Hosseini measurement. It is the band against which the γ-active magnitude
  floor (~530× under at native, density-limited) was diagnosed.
- **Slice-1 is suspended/rounded interphase — exactly Hosseini 2020's cell state.** So for the Slice-1 initial-state
  γ gate, **Hosseini 2020 (0.18–0.40, central 0.27) is the primary anchor**; the 0.35–0.65 target overlaps it only at
  0.35–0.40. The 0.35–0.65 band should be **re-scoped** (state-labeled / modeling-target) and NOT treated as the
  Slice-1 suspended-interphase requirement (the cross-check doc says exactly this: "do not treat [0.35,0.65] as a
  spread-adherent MCF7 requirement").
- **Does this change the γ-floor conclusion?** No. 0.27 vs 0.35 is ~1.3×; the measured-density floor is ~530× under.
  Re-scoping the target does not touch the "force-magnitude / load-engaged motor density limited" verdict — it only
  fixes *which* band the initial-state gate cites for the suspended cell.

---

## 4. Proposed SoT changes (Notion Contract-Graph; PI-authored)

1. **New KnowledgeClaim (dedicated γ datum), proposed KB-6.1.3-derived (e.g. KB-6.1.7):**
   > "MCF-7 cortical tension (suspended/rounded interphase) γ ≈ **0.27 mN/m** (median; **IQR 0.18–0.40**, whiskers
   > 0.05–0.62), dynamic AFM parallel-plate confinement, effective Laplace surface tension, n=27, figure-read ±0.03."
   > Evidence: SE181 Hosseini2020_AdvSci. Cross-check: Hosseini 2021 (~0.41, SE180, protocol-reconciliation pending).
   Attach the definition tag **"effective/total Laplace tension → ΔP = 2γ/R; γ_mem already included."**
2. **New/updated ValidationGate `VG-γ-MCF7-suspended`:** Slice-1 initial-state γ gate = **0.18–0.40 mN/m** (central
   0.27), state = suspended/rounded interphase. Reconcile-link to VG-H3-KU35 (0.35–0.65, re-scoped as
   modeling-target / different-state).
3. **Pressure-consistency link:** register that Hosseini γ_eff feeds the MCF7 `delta_p_hydrostatic_rest` via
   `ΔP = 2γ_eff/R` (48–107 Pa, central 72) — the PI-authored pressure-consistency gate, NOT a hand-swap. 40 Pa
   remains the HeLa diagnostic proxy until that gate lands.
4. **γ is a physiological initial-state gate the explicit cortex must satisfy — NOT a runtime lumped surface-force
   substitute** (Card 3 + the fine-grained hard rule). The gate grades the emergent cortex tension at t0; it never
   replaces the explicit actomyosin with a prescribed γ.

---

## 5. What PI is asked to ratify

1. Register the dedicated Hosseini-2020 MCF-7 γ KnowledgeClaim (§4.1) with the effective-tension definition tag.
2. Adopt `VG-γ-MCF7-suspended` = 0.18–0.40 mN/m (central 0.27) for Slice-1; re-scope 0.35–0.65 as a state-labeled
   modeling target (not the suspended-interphase requirement).
3. Confirm the **first reading** (γ_eff total → ΔP = 2γ/R, no γ_mem re-add) as the double-count resolution.
4. Keep 40 Pa as the HeLa diagnostic proxy; land the MCF7 ΔP (48–107 Pa, central ~72) only via the PI-authored
   pressure-consistency gate over {R, γ_cortex, γ_mem, ΔP_hyd}, never a hand-swap.
5. Hosseini 2021 (~0.41) registered as upper cross-check, protocol-reconciliation pending.

---

## APPLICATION LOG — 2026-07-23 (PI pre-authorized package application)

Applied per PI pre-authorization of `P2_RATIFICATION_PACKAGE_2026-07-22.md`.

**Notion Contract-Graph (SoT):**
- Created KnowledgeClaim **KB-6.1.7** (KnowledgeClaim DB): "MCF-7 cortical tension (suspended interphase)
  — direct AFM value", Value Range SI = γ_eff 0.27 mN/m (IQR 0.18–0.40, whiskers 0.05–0.62, n=27),
  effective/total Laplace tension def-tag (ΔP=2γ_eff/R, γ_mem folded in). Status = verified, Confidence =
  High. Evidence relation → SE181 Hosseini2020_AdvSci. URL `.../3a5120daec5d81469092ea247351d4f3`.
- Created ValidationGate **VG-γ-MCF7-suspended** (band 0.18–0.40, central 0.27; Type Validation; Status
  draft). URL `.../3a5120daec5d81c5b701f3cd50e2200f`.
- Updated **KB-6.1.3** Assumptions: "γ now numerically extracted → KB-6.1.7 (was figure-locked/scarce)".
- Updated **VG-H3-KU35-cortex-tension** Notes: re-scoped 0.35–0.65 as a modeling-target/different-state
  band (reconcile-link to VG-γ-MCF7-suspended); band value NOT changed (no gate-loosening).

**Code/KB (this repo):**
- `aleph/components/fluid/params_i0b1.yaml` — added `delta_p_hydrostatic_rest.hosseini2020_gamma_link` (documentation
  only) + expiry_trigger now references KB-6.1.7 / VG-γ-MCF7-suspended. **Value unchanged: 40 Pa stays the HeLa
  diagnostic proxy** (no hand-swap to 72 Pa; MCF7 dP lands via the PI-authored pressure-consistency gate only).
- `make kb-check`: [runs] gate OK · [params] gate OK · no drift (params_i0b1.yaml is not a manifest-tracked
  constant; no value edit).
