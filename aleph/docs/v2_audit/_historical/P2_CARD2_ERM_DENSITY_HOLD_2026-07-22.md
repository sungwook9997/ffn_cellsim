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

# P2 · Card 2 — MCF7 ERM linker density = production HOLD (+ scientific correction) (PROPOSAL, PI sign-off)

**Lane:** P2 evidence/SoT · **Date:** 2026-07-22 · **Status:** HOLD affirmed + correction surfaced. No production
density adopted. `make kb-check` unaffected. **SoT rule:** Notion Contract-Graph only; never edit `kb.duckdb`/vault.

Executes `AC_DECISION_CARDS_2026-07-22.md` Card 2. The ERM (ezrin/radixin/moesin) membrane–cortex linker areal
**density** for MCF7 is a genuine SoT gap; the resting-baseline diagnosis conflated two different force scales.
Production density stays **HELD**; a necessary-capacity gate is the only thing that may key off the back-computed
minimum, and even that is not the physiological density.

---

## 1. Scientific correction — `11.4 pN` is NOT the single-ERM rupture force (KB-confirmed)

The resting-baseline diagnosis mis-labeled `erm_rupture_force = 2π√(2κ_m(γ_mem+γ_MCA))` as the per-ERM cap. The KB
already isolates the two distinct scales:

- **`11.4 pN` = continuum membrane tube/tether-extraction force scale** (bilayer bending + membrane–cortex
  adhesion). This is the **KB-3.B1.4** form: `f_t = 2π√(2κ_m(T_m+γ_MCA))`, tether force f_t ~5–40 pN, γ_MCA=W
  ~1e-6–1e-4 J/m². It is a *whole-surface* pull-a-tube force, **not** a single-molecule bond force. Already tagged
  diagnostic in code ([assemble.py](../../ac/cell/assemble.py), [preload_contract.py](../../ac/cell/preload_contract.py)).
- **Single ezrin–F-actin bond (KB-3.B1.6, Braunger et al. 2014, J Biol Chem 289(14):9833–9843,
  DOI 10.1074/jbc.M113.530659, PMID 24500715, PMC3975028 — purified reconstitution, NOT MCF7 density):**
  - **stiffness k_erm ≈ 4.6 pN/nm = 4600 pN/µm** (band 2.2–4.6 pN/nm) — this IS the engine's `k_erm=4600 pN/µm`.
  - **unbinding force ≈ 50 pN**, **k_off ≈ 1.3 s⁻¹**, **barrier width ≈ 0.7 nm**.
  - Corroboration (KB-3.B1.6): Korkmazhan & Dunn 2022 Sci Adv 8(31):eabo2779 — perpendicular detach 1–4 pN (a
    *geometry-dependent* detach, not the in-line unbinding force).

**Net:** `k_erm` is correctly sourced (Braunger). What is missing is the **areal density ρ_ERM** and the **bound/
active fraction** — a single-bond stiffness/force cannot be turned into a shell capacity without them. The
correction is folded into `RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md`.

---

## 2. Evidence status of ρ_ERM (verified this session)

**There is NO MCF7 ERM areal-density claim in the KB** (`knowledge_claim` search for ERM/ezrin density,
`600 µm⁻²`, zebrafish, linker-density → empty). The values circulating in the build are:
- `1 tether / node` — a **discretization convenience**, not a measurement. Diagnostic-only.
- `600 µm⁻²` (zebrafish embryonic ectoderm, other-species embryo) — diagnostic-only proxy at best.

So the production density is genuinely un-sourced. This is the HOLD.

---

## 3. Policy (Card 2, affirmed)

1. **`ρ_ERM · F_single ≥ residual pressure` is a NECESSARY-CAPACITY gate ONLY.** It says the linker population *could*
   bear the resting membrane pressure — it does **not** license adopting the back-computed **minimum** density as the
   physiological density. Register as `VG-ERM-capacity` (necessary, not sufficient).
2. **`1/node` and zebrafish `600 µm⁻²` stay diagnostic-only.** Production latch **false** until an MCF7 direct value
   or a PI-approved epithelial proxy is registered.
3. **Approve as ONE bundle, never piecemeal:** areal density **ρ_ERM** + **bound/active fraction** + Bell kinetics
   **{k_on, k_off, F0, capture radius}**. A density without a bound-fraction over-counts load-bearing linkers; Bell
   parameters without density cannot close the transaction. The bundle is the contract unit.
4. **Proxy priority (if an MCF7 direct areal density is unavailable):**
   1. **MCF7 quantitative proteomics × cortical/membrane localization fraction** (copies/cell → areal density via
      surface area, gated by the ERM-at-cortex fraction). Highest-fidelity proxy.
   2. **Another human epithelial direct areal-density measurement** (e.g. a fixed-cell super-res ezrin count).
   3. **Other-species embryo** (zebrafish 600 µm⁻²) — **diagnostic-only**, never production.

---

## 4. The capacity-gate arithmetic (illustrative, NON-production)

At R=7.5 µm the cortex bears a residual net hydrostatic ΔP_hyd (resting proxy 40 Pa, HeLa — see Card 3 for the MCF7
pressure). Membrane-area pressure load ≈ ΔP over the shell; per-node it is ~40 pN/node at the current mesh. The
capacity gate asks: does `ρ_ERM · f_active_per_ERM ≥` the transmitted membrane→cortex traction, with `f_active_per
_ERM` bounded by the ~50 pN Braunger unbinding force and the bound fraction. This bounds ρ_ERM **from below**; the
physiological ρ_ERM is set by biology (proxy hierarchy §3.4), not by saturating the inequality. **Do not read the
minimum as the value.** (This is the same trap as the γ-floor "read the band-implied density as the datum" error.)

---

## 5. Realization-gap link (open resting blocker, not this card's to close)

`params_i0b1.yaml` `delta_p_hydrostatic_rest.realization_gap` records the live blocker: the discrete cortex is
force-free (~0.03 pN/node) and the membrane carries the full ΔP unbalanced because the ERM linkage does not
transmit force (off-radial pairing, one tether/node). That is a **mechanics/connectivity** defect (P0 lane:
FD-vs-operator + ERM tangent), **distinct** from this card's **density** gap. This card does not unblock resting;
it prevents a wrong density being adopted while P0 fixes the coupling.

---

## 6. What PI is asked to ratify

1. Accept the scientific correction: `11.4 pN` = continuum tether-extraction (KB-3.B1.4), not per-ERM; `k_erm=4600
   pN/µm`, unbinding ~50 pN, k_off 1.3/s (KB-3.B1.6 / Braunger 2014) are the single-ERM numbers.
2. Keep the **production ρ_ERM latch FALSE** (HOLD); `1/node` + `600 µm⁻²` remain diagnostic-only.
3. Register `VG-ERM-capacity` as a **necessary-only** gate (not a density source).
4. Commission the ρ_ERM **bundle** (density + bound/active fraction + Bell {k_on,k_off,F0,capture}) via the proxy
   hierarchy §3.4, MCF7-proteomics-first.

---

## APPLICATION LOG — 2026-07-23 (PI pre-authorized package application)

Applied per PI pre-authorization of `P2_RATIFICATION_PACKAGE_2026-07-22.md`. Card 2 is **HOLD-affirming** — no
production density adopted; no code value changes.

**Notion Contract-Graph (SoT):**
- Created ValidationGate **VG-ERM-capacity** (Type=Sanity check, Status=draft): necessary-only capacity condition
  `ρ_ERM·f_active_per_ERM ≥ transmitted traction`, `f_active_per_ERM ≤ ~50 pN` (Braunger 2014) × bound fraction;
  explicitly flags the back-computed minimum as a lower bound, NOT the physiological density. URL
  `.../3a5120daec5d81638c79dadd5f3889f1`.
- Updated **KB-3.B1.6** (k_erm) note: single ezrin–F-actin unbinding ≈ 50 pN activation-independent (Braunger
  2014, PMC3975028); PIP2 sets bond NUMBER (density), not force; ρ_ERM is an MCF7 SoT gap; 11.4 pN = continuum
  tube-extraction (KB-3.B1.4), not per-ERM.
- Updated **KB-3.B1.4** (tether force) note: cross-link ↔ KB-3.B1.6 disambiguation (f_t = continuum
  tube/tether-extraction scale, NOT the single-ERM unbinding force).

**Code (this repo):**
- **No code value change.** `k_erm=4600 pN/µm` already correct (Braunger).
- **HELD for PI/P0 (out of this lane's scope):** the package (c) proposes a one-line diagnostic comment on
  `erm_rupture_force` in `ac/cell/assemble.py` and `ac/cell/preload_contract.py` ("= continuum tube-extraction
  scale (KB-3.B1.4), NOT the per-ERM unbinding force ~50 pN, KB-3.B1.6"). Those are `ac/cell` files owned by the
  P0 lane; per the application-scope constraint (KB/params/ecm-card + Notion only, do NOT edit ac/cell), this
  comment is **NOT applied here** — left for the P0 lead to add. The scientific correction itself is fully
  recorded in the KB (KB-3.B1.6 / KB-3.B1.4 notes above) and in `RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md`.
- `make kb-check`: [runs] gate OK · [params] gate OK · no drift.
