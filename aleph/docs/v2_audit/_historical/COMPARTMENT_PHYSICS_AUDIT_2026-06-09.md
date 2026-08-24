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

# Compartment Physics Audit — 8 missing compartments (2026-06-09)

Adversarial physics re-audit of the 8 default-OFF compartment modules authored
2026-06-08 (`stress_fibers, linc, intermediate_filaments, microtubules,
osmotic_regulation, membrane_reservoir, junction/cadherin, junction/junctional_actin`).

## Method

Two passes, both **ground-truth-pinned by the Lead** (not agent self-report):
1. **Lead single-pass read** of each module's physics core against pinned constants
   (kT=4.28e-21 J, MT EI=2.2e-23→Lp=5.14 mm, F-actin EA≈4.3e-8 N, Young-Laplace,
   van't Hoff, Kedem-Katchalsky, discrete-WLC bridges). Found the stress-fibers
   `k_actin` order-of-magnitude error and the `cortex_myosin_` γ-blocker; verified
   MT/cadherin/osmotic correct (including disproving a suspected osmotic runaway).
2. **Exhaustive multi-agent audit** (workflow `w8w1jqamm`): 8 modules × 5 independent
   physics lenses (dimensional · literature-fidelity · sign/conservation ·
   numerics/magic-number/grid-invariance · HOOMD-API/OFF-identity) → adversarial
   refutation of *every* finding (skeptics prompted to refute, kill false positives)
   → owner-agent applies CONFIRMED findings + re-tests. 111 agents, ~6.7M tokens.
   The Lead then **re-verified every applied diff** against ground truth and
   reconciled 4 fix-agents that died on a transient socket error after editing
   (their edits were on disk, one left a broken test).

The multi-agent pass caught **three crash-on-enable bugs the single read missed**
(below) — the value of the second pass.

## Confirmed findings + resolutions

| # | Module | Sev | Finding | Resolution |
|---|---|---|---|---|
| 1 | stress_fibers | HIGH | `k_actin=1e-2 N/m` placeholder ~2-3 OOM too soft: a *validation target* (Kumar 2006 SF tension) used as a model *input*, applied as a per-bond stiffness ignoring the series-spring factor (n_beads−1). 24-bead 4µm fiber → needs 575% strain to reach 10 nN. | **FIXED**: retired the placeholder. `μ_SF = N_filaments·EA_single` (EA_single=4.3e-8 N, Kojima/Gittes 1993; N_filaments 10-30, Cramer 1997) → k_bond=μ_SF/ℓ₀≈0.8-2.5 N/m, reached at physical few-% strain. `μ_SF` defaults **None** (no single mesoscale N_filaments) → enabled build halts until PI ratifies. Kumar band relegated to validation; grid-invariance tested across n_beads∈{12,24,48,96}; CFL gate now enforced. |
| 2 | stress_fibers | MED | SF NMII reuses `cortex_myosin_*` bond types → would contaminate cortical γ if wired LIVE (those bonds ARE the active-γ signal, not `sf_`-denylisted). | **RECORDED** as an activation-blocker (registry pi_decision + module PI_DECISIONS): needs a distinct `sf_myosin_*` prefix before wiring. Harmless now (OFF). |
| 3 | linc | HIGH | `extend_snapshot_with_linc` read `snap.particles.tag` — a build-time `gsd.hoomd.Frame` has no `.tag` (only a runtime `cpu_local_snapshot` does) → **AttributeError on every enabled build**. | **FIXED**: tag-ordered convention (tag==row, matching FA/nucleus extenders) + bounds check; regression tests on a `.tag`-less Frame. |
| 4 | linc | MED | `capture_radius` default 150 nm (3×gap) gives **0 LINC bonds** at real geometry (nucleus-cloud↔cortex-shell gap = R_cell−R_nuc ≈ 5.6 µm). | **FIXED**: loud `UserWarning` when `capture_radius < R_cell−R_nuc`; regression tests at production geometry. |
| 5 | linc | LOW | `f_rest=2 pN` attributed to Arsenovic 2016, which reported *relative* FRET, not an absolute pN. | **PI-PENDING** (provenance-only — `f_rest` feeds NO force; k_linc still None→raises). Agent re-anchored to Déjardin 2020 (~8 pN) + flagged for PI/KB refresh. **Lead relaxed the test to a physiological range [1,10] pN** rather than endorse a contested exact value before PI sign-off. |
| 6 | intermediate_filaments | HIGH | A standalone `md.bond.Harmonic` — HOOMD 7.0.1 requires params for **every** bond type in the state, so it **crashes any integrated cell** carrying other bond types. | **FIXED**: `register_if_bond_params` writes if_ params onto the cell's single shared `md.bond.Harmonic` (mirrors stress_fibers/linc); standalone attach now refuses foreign bond types. Un-skipped + rewrote the heavy test as a real coexistence test (`cortex_bond + if_backbone` on one Harmonic, `sim.run(0)`). |
| 7 | intermediate_filaments | HIGH | `if_crosslink` used `r0=0` on already-separated cross-filament beads → a contractile spring injecting ~1.95e6 kT spurious pre-tension at t=0 (physiological-baseline violation). | **FIXED**: per-r0-bin force-free `if_crosslink_b{i}` (platform convention). Construction energy 227,499 kT → 126 kT (~O(kT)·n_xl thermal); born force 34.5 pN → 0.7 pN. |
| 8 | intermediate_filaments | LOW | Block 2018 cited as BOTH "Sci. Adv." and "PRL" (self-contradictory). | **FIXED**: single verified citation Sci. Adv. 4(6):eaat1161, DOI 10.1126/sciadv.aat1161. |
| 9 | microtubules | MED | snapshot-extension None-typeid handling + the (correct) stiff-bending CFL was an aspirational docstring claim, not enforced. | **FIXED**: hardened snapshot extension; `attach_microtubule_forces` now RAISES if host dt exceeds the bend CFL (`cfl_strict`). Core physics (EI/Lp/k_angle/τ_bend) was already correct — unchanged. |
| 10 | osmotic_regulation | MED | (defensive) a hand-built dataclass bypassing the resolver gate must still be rejected at Updater construction. | **FIXED** (guard added). Lead fixed the dead-agent's broken test precondition (batch_dt now genuinely > τ_Kvol). Sign/stability re-confirmed correct (stable negative feedback to ΔP_target). |
| 11 | membrane_reservoir | MED | `W_MCA` membrane-cortex adhesion energy unanchored. | **FIXED**: `W_MCA=1e-5 J/m²` in the cited KU-3.B1.4 / Hochmuth 1996 / Derényi 2002 band [1e-6,1e-4]. Bleb `σ_crit` + reservoir `f_excess` remain **None→raise** (PI-pending). γ-denylist `mem_` recorded in the registry. |
| 12 | cadherin_junction | HIGH | Same standalone-`md.bond.Harmonic` crash as IF (#6). | **FIXED**: foreign bond types registered with `dict(k=0,r0=0)` on the shared force. Catch-bond kinetics still correctly delegated to the validated `validation/cadherin_sliding_rebinding.py` oracle. |
| 13 | junctional_actin | DOC | catch form is a parallel-Pereverzev surrogate, not Buckley 2014's sequential two-state. | **CLARIFIED** in PI_DECISIONS (form/source distinction on record). STUB constants remain None→raise. |

## False positives refuted (did NOT change code)
- **osmotic V0-update "runaway"** (Lead's own initial suspicion): linearization gives
  `∂V̇0/∂V0 < 0` → stable. Confirmed correct, no change.
- Several lens findings on MT/cadherin downgraded to NONE by the adversarial verify
  stage (the core derivations were already right).

## Net
- **3 crash-on-enable bugs** (LINC `.particles.tag`; IF + cadherin standalone Harmonic)
  and **1 physiological-baseline violation** (IF r0=0 pre-tension) fixed — none would
  have been caught without actually trying to build the enabled path.
- **1 quantitative error** (stress-fibers k_actin) fixed by re-derivation to None-gated `μ_SF`.
- Citation integrity tightened (Block 2018; LINC Arsenovic→PI-pending; cadherin Iturri SE flag).
- **All physics fixes are literature-cited or grid-invariant-derived; genuinely-unknown
  constants (μ_SF/N_filaments, k_linc, σ_crit, f_excess, junctional catch set) default
  None → enabled path raises.** No invented numbers.
- Verification (Lead, ground-truth): all 8 module suites + registry + manifest +
  cortical_tension + production_policy **green** (290 passed, 1 heavy-skipped); IF heavy
  test now runs (coexistence). No existing runtime file (`cell.py`, `manifest.py`,
  `cortical_tension.py`, `integrator/`) touched.

## Still PI-pending (unchanged activation blockers)
stress-fibers `N_filaments`/μ_SF + `sf_myosin_*` γ-prefix + bundle bending EI · LINC
`k_linc` + `f_rest` attribution · IF nonlinear strain-stiffening · membrane_reservoir
`σ_crit`/`f_excess` · junctional_actin full catch set. Each keeps its compartment
disabled (build raises) until ratified.
