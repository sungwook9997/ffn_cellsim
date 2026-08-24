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

# Compartment-Platform Activation — Session Closeout (2026-06-09)

**Branch:** `h7/compartment-platform`  ·  **Working model:** autonomous platform
session (sibling Gate-B session on `h7/full-cell-integration` untouched).
**Mandate:** PI 소유권 허용 — graduate the remaining default-OFF compartments
EXPERIMENTAL/STUB → LIVE via the validated 5-step activation procedure.

## TL;DR
**The activation backlog is COMPLETE.** All **8** default-OFF compartments are now
**LIVE**, each with a passing build-time activation gate, no cortical-γ
contamination, force-free construction, and a bit-identical disabled path. No
EXPERIMENTAL/STUB compartment remains. Plus: a new two-cell doublet assembler, a
SimuCell3D-style versioned morphology visualizer (PI request), and several
crash/overflow fixes.

## Compartments activated this session

| # | Compartment | Commit | Gate | One-line |
|---|---|---|---|---|
| ① | `linc` | `d8d6770` | 5/5 PASS | nucleus↔IF-cage nesprin bridges (Option A), per-bond EXACT-r0 force-free |
| ② | `membrane_reservoir` | `c548ada` | 6/6 PASS | own `mem_node` offset layer + static `mem_tether` mesh (BLOCKER-1 fix) |
| ③ | `ventral_stress_fibers` | `990e51c` | 5/5 PASS | PASSIVE backbone, **long-axis-aligned** basal FA pairs (PI-decided) |
| ④ | `cadherin_junction` | `ab70525` | 5/5 PASS | FIRST multicell — new `build_cell_doublet` two-cell assembler |
| ⑤ | `junctional_actin` | `306bb95` | 5/5 PASS | reserved STUB build IMPLEMENTED — α-catenin cadherin↔cortex catch clutch |

(osmotic_regulation · microtubules · intermediate_filaments were activated in the
prior session; with these 5 the full set of 8 is LIVE.)

**Visualization (PI requests):** `d570f66` (first morphology viz) → `ad895b6`
(SimuCell3D-style surfaces + cutaway + **versioned** output) → `214a993` (v02 with
the adherent stress-fiber cell).

## How each was done (the common pattern)
Every activation followed the validated 5-step procedure: (1) None-gated constants
set to DETAILED_PLANS candidates with provenance/confidence labels (module default
stays `None` → un-anchored build still raises); (2) snapshot-extension wiring into
`cell.py`/`manifest.py` (+ LJ `r_cut=0` + cytoplasm `gamma_map` for new particle
types); (3) registry-driven γ-denylist auto-excludes the new bonds; (4) a build-time
activation gate (full physiological baseline OFF vs ON) with PASS/PARTIAL/REFUTE
stated before the run; (5) tests (off-identity + on-assembled + no-contamination) +
`scripts/h7_<name>_activation_gate.py` (auto-viz) + an AUTORUN_LOG line.

**Universal controls that PASSED for all activations:**
- **No cortical-γ contamination** — `γ_soft` OFF == ON identical (registry-driven
  denylist excludes every new bond prefix; cortex/myosin signal preserved).
- **Force-free construction** — per-bond/per-bin EXACT-r0 (max strain ~1e-14).
- **Disabled-path bit-identity** — every compartment OFF reproduces the baseline.

## New infrastructure
- **`cell/doublet.py` — `build_cell_doublet`** (the biggest new piece): two cortex
  shells in one box (offset `2R+gap`, merged frame with bond/angle re-indexing),
  cadherins seeded on matched facing caps, trans-dimers **seeded pre-bound** (the
  engaged-junction baseline; a bare junction needs ~1e3 binder batches to engage),
  cadherin/junctional anchors hold particles on the surface, catch-slip binder
  attached. `with_junctional_actin=True` adds the α-catenin belt per cell.
- **`scripts/compartment_vis.py`** — single entry point; builds the FULLEST cell and
  renders **SimuCell3D-style closed surfaces + a cutaway** exposing the internal
  architecture (nucleus, MT aster, IF cage, LINC, stress fibers). **Versioned**
  output under `outputs/h7/figs/morphology/cell_vNN_<label>.png` (accumulates, never
  overwrites) + `cell_latest.png`.

## Fixes
- **2 nlist exclusion-cap overflows** (LINC + membrane): many acceptors shared one
  cortex/if bead → degree spike. Fixed with degree-aware / unique-acceptor 1:1
  matching (bounds per-acceptor bond degree to +1).
- **cadherin binder `n_sub` overflow guard** (extreme-slip dynamic run RNG crash).
- **SF registry denylist** corrected to `('sf_',)` (dropped `cortex_myosin_`, which
  the estimator never denylists).

## Test status
Platform suite **green**: LINC 45 · membrane 30 · SF 33 · cadherin 34 · junctional
21 · registry · canary · IF · MT, plus the full compartment/cell/cortex suites.
The **only** failing tests in the whole repo are the **8 pre-existing
`test_ku35_patch_contractility.py` failures** — confirmed failing on a clean
(stashed) tree, unrelated to this platform work (a separate KU-3.5
production-policy-baseline issue).

## Figures (per the visualize-at-closeout rule)
- `outputs/h7/figs/morphology/cell_v01_membrane.png` — suspended cell, 5 internal
  LIVE compartments (3D surface + cutaway + radial profile).
- `outputs/h7/figs/morphology/cell_v02_stressfibers_adherent.png` — adherent cell,
  all 6 single-cell LIVE + FA + ventral stress fibers (24,938 particles).
- `outputs/h7/figs/morphology/cell_latest.png` — pointer to the newest.
- `outputs/h7/figs/h7_{linc,membrane_reservoir,stress_fibers,cadherin_junction,
  junctional_actin}_activation_gate.png` — per-compartment gate figures.

## Open items → PI (DEFERRED; these are NOT activations — see `PLATFORM_PI_QUEUE.md`)
All remaining work is active/dynamic-phase physics or PI parameter ratification:
- **Parameter ratifications:** `k_linc` (1e-2 route-A), `N_filaments` (20), the
  junctional catch-set, `σ_crit_bleb` / `f_excess` (membrane bleb/reservoir), and
  **registering the Iturri-2020 SourceEvidence** row (cited but unregistered).
- **Active/dynamic phases (need an equilibrated run — the raw full-cell run trips the
  BAOAB guard, so an equilibration prelude is required):**
  - cortex + SF **NMII** active contraction via the `sf_myosin_*` prefix split in
    `cortex/myosin.py` (then the Kumar 10-30 nN single-SF tension gate);
  - cadherin / junctional **catch-slip maintenance** (the binders running to the
    engaged-fraction equilibrium; the Iturri ~6.5 nN ensemble de-adhesion observable);
  - IF **nonlinear strain-stiffening** Table law (Kreplak/Block);
  - membrane **bleb** rupture + reservoir-release buffering;
  - **dense junctional belt** — the single-particle cadherin is the ectodomain tip at
    the interface (~0.5 µm from cortex), so only tips within the α-catenin reach
    couple → a sparse belt; a faithful dense belt needs a cadherin-tail particle.
- **MTOC multi-bead core** to unlock a dense aster (`n_mt` > 7; nlist exclusion cap).
- **Per-cell full internal stacks on a doublet** (each doublet cell with
  nucleus/MT/IF/membrane/SF) — run the single-cell extenders per cell.

## State stores
- Running log: `PLATFORM_AUTORUN_LOG_2026-06-09.md` (per-step table + Phase-2 summary).
- PI decisions / deferred: `PLATFORM_PI_QUEUE.md`.
- Smoke + figures index: `COMPARTMENT_SMOKE_REPORT_2026-06-09.md`.
- This closeout: `PLATFORM_ACTIVATION_COMPLETE_2026-06-09.md`.
