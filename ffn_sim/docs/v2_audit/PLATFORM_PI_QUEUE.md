# Platform PI Queue — 2026-06-09

## ✅ ACTIVATED under the PI ownership grant (2026-06-09 "소유권 허용 / 모든 것 달려")
- **osmotic_regulation → LIVE** (gate PASS: τ_RVD=3.2s in band, RVD sign, no-contamination,
  off-identity). Post-build attach in manifest.py (no cell.py surgery).
- **microtubules → LIVE** (gate PASS: aster assembled, CFL passed at cytoplasm drag,
  γ_soft IDENTICAL OFF/ON = no contamination, L_p in band). Full snapshot-extension
  wiring in cell.py + manifest.py. **n_mt CAPPED at 7** (see new blocker below).
- **intermediate_filaments → LIVE** (gate PASS: cage assembled +480 part, γ_soft
  IDENTICAL OFF/ON, if_ denylisted). LINEAR small-strain path; the nonlinear
  strain-stiffening Table law stays PI-pending (NotImplementedError) → SECONDARY
  stiffening gate deferred. Fixed 2 crash-on-enable bugs (gsd-None velocity/angles
  guards in the IF extender — masked by the smoke, caught by the real build).
- **linc → LIVE** (gate PASS: bridges formed n=40, FORCE-FREE per-bond EXACT-r0,
  γ_soft IDENTICAL OFF/ON = no contamination, off-identity +0 particles, CFL ~4-order
  headroom). Option A: bonds nucleus_bead → nearest perinuclear if_bead (the LIVE IF
  cage); requires=('nucleus','intermediate_filaments'). Per-bond EXACT-r0 (FA-clutch
  convention) so each bridge is force-free at its as-built separation. Added
  `unique_acceptor` 1:1 matching (one nesprin per IF anchor) to bound per-if_bead bond
  degree — WITHOUT it the nlist per-particle exclusion cap overflows (same class as MTOC).
  **PI items to ratify (LINC is LIVE but these are pending):**
  - **k_linc = 1e-2 N/m** (route-A folded-rod pre-unfolding secant, Rief 1999; ORDER, conf
    MEDIUM). Module default stays None (un-configured enabled build raises); the config
    carries the candidate. PI to ratify route A vs B (WLC entropic ~6e-6) vs C (md.bond.Table
    nonlinear). Also the 2→8 pN f_rest re-anchor + Déjardin-vs-Arsenovic attribution + the 6
    SourceEvidence rows (Rief99, Déjardin20, Arsenovic16, Autore13, Crisp06, Lombardi11) before
    any deliverable cite.
  - ⚠ **Mesoscale-r0 fidelity caveat (FOLLOW-UP).** At the ×40 mesoscale the nucleus↔if_bead
    bead separation is ~0.5-1.5 µm (NOT the 50 nm real nesprin span), so r0 is an EFFECTIVE
    mesoscale coupling length and the 8 pN f_rest oracle is approximate. **Option C** (seed a
    dedicated R_nuc+50 nm perinuclear acceptor cap — a new particle layer) is the faithful-span
    alternative; surfaced per the "fidelity-impact → PI queue" rule. The gate's no-contamination
    + force-free + bridges-formed controls hold regardless; the [2,10] pN <T_linc> resting-tension
    oracle is the deferred physics gate (needs an equilibrated actomyosin run; raw full-cell run
    trips the BAOAB guard — equilibration prelude required, same as the other compartments).
  - ⚠ **GPU many-bond-types perf (FOLLOW-UP).** Per-bond EXACT-r0 mints one bond type per bridge
    (~n_IF_filaments); combined with IF's 16 crosslink bins HOOMD warns "many bond types perform
    poorly / shared-memory errors on the GPU." Fine for the CPU build-time gate; flag for the GPU
    production path (a coarser shared per-r0 bin needs ~hundreds of bins at k_linc=1e-2 to stay
    thermal, so exact-r0 is actually the leaner choice — a native bond plugin or a stiffness-aware
    binning is the real fix).
- **cortical_tension γ-denylist → registry-driven** (`NONCORTICAL_COMPARTMENT_PREFIXES`
  = `REGISTRY.gamma_denylist()` minus `cortex_*`). Unblocks the no-contamination control
  for every γ-contaminating compartment. ⚠️ sibling Gate-B owns this file → merge-coordinate.

## ⚠️ NEW BLOCKER found during activation
- **MTOC single-hub degree vs HOOMD nlist exclusion cap (7).** The aster's single MTOC
  carries n_mt backbone bonds; with the full-cell LJ nlist on (`exclusions=("bond","1-3")`),
  n_mt>7 overflows the compile-time per-particle exclusion cap → `Too many bonds to process
  exclusions`. So a single-hub MTOC supports only n_mt≤7 (a sparse aster, now LIVE). A DENSER
  interphase aster (PI candidate 20, up to ~250) needs a **multi-bead MTOC core** in
  `cell/microtubules.py build_mt_topology` (distribute arms across ceil(n_mt/5) core beads) —
  a topology change that churns the 35 MT tests + 3 smokes. Filed for a focused follow-up.

---
# Original queue — 2026-06-09

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
