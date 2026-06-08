# Platform PI Queue — 2026-06-09

Items the autonomous compartment-platform session (branch `h7/compartment-platform`)
could NOT do because they require **PI parameter ratification**, **shared-file edits**
(`cell/cell.py`, `cell/manifest.py`, `cortex/cortical_tension.py`, existing `cortex/*`
runtime, `integrator/`, `native/`), or a **gate-contract change** — all out of this
session's ownership. Each was SKIPPED here and recorded for PI. None block the
smoke-harness plumbing work (which used standalone harnesses + SMOKE-ONLY values).

The ENABLED-PATH smoke harnesses (`scripts/compartment_smoke/`) confirmed every
compartment's enabled build/step path runs (or is honestly blocked); the items below
are what stands between that plumbing and a real, PI-signed ACTIVATION.

## A. Shared-file edits required for activation (PI-gated; not owned here)

1. **cortical_tension γ-denylist extension (ALL γ-contaminating compartments).**
   `cortex/cortical_tension.py` `ADHESION_BOND_TYPE_PREFIXES` is hardcoded to FA
   prefixes only. Activating any of `ventral_stress_fibers (sf_, sf_myosin_)`, `linc
   (linc_)`, `intermediate_filaments (if_)`, `microtubules (mt_)`, `membrane_reservoir
   (mem_)`, `cadherin_junction (cadherin_)`, `junctional_actin (junc_actin)` needs its
   prefix added (cleanest: have the estimator read `REGISTRY.gamma_denylist()`), with a
   contamination regression test. **Shared file → PI sign-off.**

2. **stress_fibers BLOCKER #2 — `sf_myosin_*` prefix split.** SF NMII currently reuses
   the cortex `cortex_myosin_*` bond types (the active-γ signal). A LIVE SF build would
   contaminate cortical γ. Needs `cortex/myosin.py` parameterized to mint a distinct
   `sf_myosin_*` prefix. **Shared file (`cortex/myosin.py`) → PI sign-off.** Until then
   only the PASSIVE SF backbone is runnable (smoke#6 ran passive, NMII off).

3. **5-point wiring into the loader (every compartment).** `manifest.resolve_baseline`
   + `cell.py` (`_extend_snapshot_*`/`attach_*` + `CellBuildOptions` flag + `gamma_map`
   entry) + `configs/mcf7_baseline.yaml` `optional_subsystems` block + registry
   `status→LIVE`/`manifest_path`. `cell/cell.py` + `cell/manifest.py` are not owned.
   **PI-gated** (gate-contract: status flip).

## B. PI parameter ratifications (None-gated; smoke used SMOKE-ONLY placeholders)

These keep each enabled build raising until ratified (DETAILED_PLANS has candidates):
- **ventral_stress_fibers:** `N_filaments` (→ μ_SF); bundle bending EI (deferred).
- **linc:** `k_linc` (route A/B/C; smoke used 1e-2 N/m route-A candidate).
- **intermediate_filaments:** the nonlinear strain-stiffening Table curve (smoke ran the
  linear path only; nonlinear stays NotImplementedError).
- **microtubules:** `Y_stretch` (the dt lever) + `L_mt` host geometry; DI rates (default-OFF).
- **osmotic_regulation:** `Lp` (smoke used 1e-12 MCF7/AQP5 candidate; default 1e-13).
- **membrane_reservoir:** `σ_crit_bleb` + `f_excess` (bleb + reservoir-release paths blocked).
- **cadherin_junction:** k_trans/r_bind derived (OK to run); Iturri-2020 SourceEvidence
  registration before any deliverable cite.
- **junctional_actin:** the α-catenin catch set `k_catch0`/`k_slip0` (ORDER_ESTIMATE).

## C. STUB build-path implementations (not owned / PI-sanctioned release)

1. **junctional_actin HOOMD build path.** `extend_snapshot_with_junctional_actin` raises
   `NotImplementedError` EVEN when anchored (the explicit head-placement + per-r0 coupling
   topology is reserved for the PI-sanctioned release). Smoke#8 ran the SCALAR laws only
   (catch shape F*≈6 pN, engaged f≈0.91) + asserted the build barrier. Implementing the
   build path is a STUB→EXPERIMENTAL graduation; depends on `cadherin_junction` GATE-J
   passing first (junctional_actin requires cadherin_junction).

2. **membrane_reservoir bleb rupture updater.** `MembraneTetherUpdater` raises
   `NotImplementedError` unconditionally (the Bell-Evans `act()` loop is unwritten AND
   σ_crit_bleb is None). The static tether mesh is fully runnable (smoke#7); the rupture
   path needs the act() implementation + σ_crit_bleb ratification.

3. **microtubules dynamic instability updater** (`MTDynamicInstability`) raises
   `NotImplementedError` (plus-end snapshot-rebuild unwired); default-OFF, separate PI task.

## D. Two-cell / doublet build (cadherin activation)

`cadherin_junction` real activation (GATE-J) needs a `build_cell_doublet()` two-cell
assembler — `cell.py` is single-shell. Smoke#5 used a minimal facing-patch standalone
build to validate the binder; the production doublet is a `cell.py` change → PI-gated.

---
*Logged by the autonomous platform session; see `PLATFORM_AUTORUN_LOG_2026-06-09.md`.*
