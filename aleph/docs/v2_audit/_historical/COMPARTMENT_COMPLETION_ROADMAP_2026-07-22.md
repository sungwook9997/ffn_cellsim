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

# Compartment-completion roadmap — `ac/` engine (Lead plan, 2026-07-22)

> Lead-session plan for finishing the 10-compartment cell. Grounds every step in the current
> `build_cell` composition, the [`CELL_MECHANICS_FRAMEWORK_2026-07-15.md`](CELL_MECHANICS_FRAMEWORK_2026-07-15.md)
> target, and the hard rules (no magic numbers → source or GAP-to-PI; validate at full native; physiological
> baseline; no gate-loosening). Companion to [`BLEBBING_VALIDATION_PLAN_2026-07-20.md`](BLEBBING_VALIDATION_PLAN_2026-07-20.md).

## 0. Two hard constraints that shape the sequencing

1. **Resting convergence is currently BLOCKED** (verified 2026-07-21): `--from-resting` does not accept —
   projected-force stalls (no-steric 0.34 pN @6k, full-native 2634 pN from 64,517 steric interpenetrations).
   Root = (a) build interpenetration (overlap-free packing needed) + (b) coupled membrane→ERM→fiber→shell
   ill-conditioning (cond ~3e4, needs exact surface/pressure tangent + multilevel coupled preconditioner).
   Codex owns this in `implicit_mechanics.py`. **Composing MORE compartments into the production cell is
   premature until this lands** — new load paths would add to a non-converging solve.
2. **Most missing compartments need SOURCED parameters** (ECM stiffness, α2β1 clutch Bell params, MT κ, IF
   moduli, chromatin WLC, fascin). The magic-number rule forbids inventing them → they must be
   literature-sourced (the k_erm/Braunger pattern) or surfaced to PI, before implementation.

**Consequence — the safe overnight order is:** ① source params (parallel research → KB) → ② build the
net-new compartments as **isolated, unit-tested modules** (not composed) → ③ enhance visualization
(no convergence needed) → ④ compose + validate at full native **after** the convergence fix lands.

## 1. Current production composition (`assemble.py::build_cell`)

`weave_cell([_cortex_region])` + optional `with_{myosin,steric,nucleus,membrane,pressure}`. **Only the CORTEX
region is woven**; every other region builder (SF/lamellipodium/filopodium/arc/cap) exists in
`ac/weave/regions.py` but is never passed in. No `with_stress_fibers/with_ecm/with_adhesion` flag exists.

## 2. Per-compartment status + completion plan

| # | Compartment | Status | Completion action | Blockers |
|---|---|---|---|---|
| 1 | **Cortex** | ✅ DONE | — (short-filament census is a later refinement) | none |
| 2 | **Membrane** | ✅ DONE | k_erm now sourced (KB-3.B1.6). ERM density/kinetics → NG-3 (PI) | ERM density GAP |
| 3 | **Nucleus** | 🟡 PARTIAL | wire **chromatin WLC net** (`chromatin_analytic.py` exists, no kernel launch) + nucleoplasm viscosity; I7 cap→LINC | chromatin/viscosity params (sourcing ▶) |
| 4 | **Fluid** | 🟡 PARTIAL | compose **RAD G-actin monomer transport** (`fluid/transport.py` unimported) | `c_0`, `φ` GAP (PI) |
| 5 | **Myosin/NMII** | 🟡 PARTIAL | drive heads into the **bound contractile state** (`step_myosin_kinetics`/`commit_segment_kinetics` exist) — the central active term | k_xb/f_stall/v0/κ provisional; needs convergence |
| 6 | **Stress fibers** | 🔴 code-not-composed | give `stress_fiber.py` a **Warp tension kernel** (returns `tension=None` now) + weave SF/arc/cap regions + SF→FA→ECM + cap→LINC | SF force magnitude GAP; ECM dep |
| 7 | **Lamellipodium** | 🔴 builder-not-composed | wire the dendritic-Arp2/3 builder + membrane-ratchet force | branch stiffness params (sourcing ▶) |
| 8 | **Filopodia** | 🔴 seed only | net-new **fascin-bundle** module | fascin params (sourcing ▶) |
| 9 | **Microtubules** | 🔴 MISSING | net-new **MT strut** module (reuse `ff/microtubule.py` physics) — tensegrity compression | MT κ/buckling/count (sourcing ▶) |
| 10 | **Intermediate filaments** | 🔴 MISSING | net-new **keratin/vimentin cage** module (EMT contrast) | IF moduli (sourcing ▶) |
| 11 | **ECM + adhesion clutch + FA** | 🔴 MISSING | net-new **collagen network + α2β1-collagen Bell clutch + FA maturation** — the outer BC + traction | ECM/clutch/FA params (sourcing ▶); highest impact |
| 12 | **Crosslinkers** | 🟡 PARTIAL | wire dynamic **catch/slip Bell KMC** (`crosslink_kmc.py` host oracle) into the device driver | needs convergence (turnover in the loop) |

## 3. Prioritized work streams (ranked by impact × tractability)

**Stream A — ECM + adhesion (highest impact, fully net-new).** The framework's outer boundary + traction +
mechanosensing all require it. Sourcing ▶ (collagen-I fibril modulus, α2β1 catch-bond Bell, talin clutch
stiffness, FA traction). Build as: collagen fiber network (reuse `ff/` ECM library physics) + α2β1-collagen
Bell-Evans clutch kernel + FA maturation. Isolated + unit-tested first; compose as `with_ecm`/`with_adhesion`.

**Stream B — Microtubules (tensegrity strut).** Net-new `ac/solid/microtubule.py` reusing `ff/microtubule.py`
κ-physics: MT as compression strut (Lp~1-6 mm, reinforced buckling), dynein point forces, actin↔MT coupling.
Isolated module + Warp bending kernel + CFL + unit tests. Sourcing ▶.

**Stream C — Intermediate filaments (EMT contrast).** Net-new `ac/solid/intermediate_filament.py`:
strain-stiffening keratin(MCF7)/vimentin(MDA-231) cage, huge extensibility, perinuclear cage + IF-nucleus.
Isolated module + nonlinear force-extension kernel + tests. Sourcing ▶.

**Stream D — Complete PARTIAL compartments.** Chromatin WLC kernel + nucleoplasm viscosity (nucleus);
RAD monomer transport composition (fluid); dynamic crosslinker Bell KMC. These need sourced params + (for
the loop-composed ones) the convergence fix.

**Stream E — Myosin activation.** Drive the assembled cell into the bound contractile state (the central
active tension). Gated on convergence (contraction needs accepted outer steps).

**Stream F — Visualization (continuous, no convergence needed).** Extend `ac_cell_assembled_viz.py`: bleb
patch geometry, ERM lines, merged pressure-field scenes, per-compartment colour, and each new compartment as
it lands. Browser-verified. The user-facing "see the cell" deliverable, refreshed at every milestone.

## 4. Sequencing + gating

- **Now → overnight**: Streams A/B/C sourcing (▶) → isolated module builds + unit tests (Mac source-check +
  gbook A5000 kernel gates); Stream F viz enhancements; roadmap + design specs.
- **After convergence fix (codex)**: compose the isolated modules + Streams D/E into the full native cell;
  validate at 70,686 + all compartments; re-run bleb growth; measure emergent tension/traction.
- **PI sign-off gates**: every new sourced param → KB KnowledgeClaim + verdict OK before use; new
  ValidationGates PI-authored; `ffn/foundation` push PI-gated; force-magnitude GAPs surfaced not tuned.

## 5. Live status (2026-07-22, updated)

**Built + A5000-verified tonight (isolated modules, reuse ff/ physics, §1.4 accumulate contract):**
- ✅ **Microtubule strut** `ac/solid/microtubule.py` — 5/5 (Euler buckling, L_p, straight→0, bent→restoring).
- ✅ **Intermediate-filament cage** `ac/solid/intermediate_filament.py` — 3/3 (force-free, anchor, nucleus LINC; keratin/vimentin EMT split).
- ✅ **ECM adhesion-clutch** `ac/solid/adhesion_clutch.py` — 4/4 (α2β1 slip off-rate, engaged traction −z, unbound=0).
- Consolidated: **12/12 on the gbook A5000**, no regressions.

**KB registered (Notion SoT):** KB-3.24 (MT EI/Gittes), KB-3.25 (MT buckling/Brangwynne), KB-3.26 (IF l_p keratin/vimentin), KB-3.27 (IF extensibility), KB-3.28 (α2β1 slip/Attwood), KB-3.29 (clutch κ_c/Chan-Odde) + 6 DOI-verified SourceEvidence.

**Also:** cell viz refreshed + browser-verified; 7-agent DOI dossier; MT+IF+ECM structure figure.

**Convergence recheck (blocks composition):** WCA-harness merge added new inner solvers (analytic_implicit/
block_descent/anderson/rkc1/contact_schwarz/tournament/contact_tournament). Full-native contact_tournament
reduces residual 2634→1483 pN but still does NOT accept (steric build-interpenetration dominates); no-steric
tournament also does not accept (coupled-preload stall). **Resting-convergence fix not yet landed →
composition + bleb growth + full-native validation remain gated (codex active).**

**⏭ Next:** filopodium (fascin bundle) + PARTIAL wiring (chromatin WLC / RAD monomer / dynamic crosslinker KMC,
params sourced) as isolated modules; on convergence landing → compose MT/IF/ECM (with_microtubules/with_if/
with_ecm) + re-run bleb growth + full-native validate.

## 6. Convergence #2 — fiber conditioning: global coarse-mode deflation (Lead, 2026-07-22)

**Landed (isolated-verified, opt-in, gate-neutral):** the first piece of directive ① (coupled multilevel
preconditioner). `ProjectedAnalyticCG` now carries a **global rigid-body + constant-strain (l ≤ 2) coarse
space** — 3 translations + 3 rotations + 6 constant symmetric strains = 12 modes, exactly the
`l = 0` breathing + `l = 1` translation/rotation + `l = 2` ellipsoidal space. This is the classical elastic
near-kernel: geometry-derived, **no magic number** (orthonormalizing only conditions the coarse matrix), and
it couples membrane⊕ERM⊕fiber⊕nucleus through one global smooth displacement — the residual family the
per-fiber block and node-Jacobi cannot see (the "global radial 17%" of the stall decomposition).

Mechanism: a two-level additive **Galerkin deflation** `z += B (Bᵀ A B)⁻¹ Bᵀ r`. `A_c = Bᵀ A B` is rebuilt each
solve from the live analytic operator (`coarse_modes` operator applies), factored by a device dense Cholesky,
and solved on-device (no host round-trip). **Preconditioner only** → it cannot move the fixed point or relax
any residual gate; it only accelerates the coupled global mode. Wiring: `ProjectedAnalyticCG(coarse_modes=)`
→ `make_inner_solve(implicit_coarse_modes=)` → `run_from_resting(...)` / driver `--implicit-coarse-modes`
/ preload probe `--implicit-coarse-modes`. **Default 0 (OFF)** — a numerical accelerator lands opt-in so it
carries zero regression risk against Codex's tournament (mirrors the `overlap_free_cortex` OFF-by-default
landing); the physiological-baseline "turn it ON" rule is about physics setpoints, not numerics.

Isolation verification (like the membrane tangent): device coarse term matches the NumPy Galerkin oracle
`two_level_coarse_correction` to ~1e-11, rigid modes are exact zero-energy of a central-force network (~1e-16),
and `coarse_modes=12` returns the SAME converged displacement as the dense bending oracle (solution unchanged,
confirming gate-neutrality). Tests: `test_{rigid_strain_coarse_basis_is_orthonormal_and_rigid_null,
two_level_coarse_correction_matches_exact_galerkin, warp_global_coarse_matches_two_level_reference_and_preserves_solution}`.

**⏳ Efficacy PENDING synced gbook:** whether deflation breaks the coupled-preload stall must be measured on a
single-commit-aligned gbook (Codex `SYNC NOW`) — re-run the `bqdo6iq4k`-style probe with
`--implicit-coarse-modes 12` and compare projected-force gate + accept vs `0`. gbook is partial-synced now, so
acceptance is not yet trustworthy. **Follow-ups in ①:** overlapping Schwarz (membrane-node⊕ERM⊕fiber) +
Schur/coarse block (②→ builds on the stored `A B = coarse_action_d`); if additive deflation double-counts on
easy steps at native scale, upgrade to the balanced/BNN form (also reuses stored `A B`, no extra apply).
