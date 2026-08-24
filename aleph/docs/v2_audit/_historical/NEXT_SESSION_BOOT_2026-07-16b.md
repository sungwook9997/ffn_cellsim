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

# Next-session boot — 2026-07-16 (b): internal-coupling build state

Branch `ff/mech-hierarchy`. This session (10 commits) built the nucleus↔interior coupling axis from the
internal-displacement diagnosis. Read this + CLAUDE.md + the linked docs, then `git log --oneline -11`.

## What landed (all committed + native-validated, honest)
1. **Geometry fix** (4c0f1d9): MCF7 R_cell/R_nuc sourced distribution, retired phantom "Moore 2016" cite.
2. **FSI wiring** (063fc4e): spatial Biot pore-pressure two-way FSI into the cortex (6πηR→γ_solid retired).
3. **Internal-displacement DIAGNOSIS** (340e80b, native): ⭐ the interior (nucleus + MT) is mechanically
   DECOUPLED — RAW |Δ| = 0.0 nm under AFM while the cortex deforms. Harness `ff_internal_displacement_diag.py`.
4. **MT-strut fix** (765655a): MT thickness (25nm radius) + reach-cortex → MT UNFROZEN (native RAW 0→8 nm).
5. **FSI → ALL compartments** (5b7e2c9): the Biot fluid couples to nucleus (INWARD shell force) + MT (−∇p,
   negligible). Nucleus responds ∝ M_biot but SUB-nm at physical M — fluid is SECONDARY for the nucleus.
6. **M_biot physical value** (e143329): 300→1000 Pa (undrained-derived, KB-DRAFT). Honest: fluid couples to
   everything but doesn't dominate nucleus deformation.
7. **IF cage** (fef3c78): radial-spoke cage (LINC + cortex-anchor, linear-first). ⚠ NATIVE finding: does NOT
   flatten the nucleus (RAW 0.1 nm) + over-stiffens cortex (F_plate 6×) — HALT/surface (IF_CAGE_2026-07-16.md).

## The honest through-line (report this to PI on boot)
**No coupling has yet FLATTENED the nucleus.** Fluid = secondary (deviatoric ~Pa); MT = unrelated; IF cage =
over-stiff + buckles (polar radial spokes buckle under axial compression; IF has no bending). ⭐ **Key
physics insight: strain 0.10 is a GENTLE regime where the plate (6.93 µm) does NOT reach the nucleus (5.1 µm)
— real AFM flattens the nucleus only at strain > 0.34 (plate reaches it). A small nucleus response at gentle
strain may be PHYSICALLY CORRECT; I may have been trying to force nucleus motion.** Re-frame with the PI.

## Ready to implement (designed, verified, NOT built)
- **Membrane Helfrich bending** (MEMBRANE_HELFRICH_DESIGN_2026-07-16.md): the membrane is NOT OK (production =
  lumped 2γ/R; real sheet un-wired; Helfrich bending absent everywhere — PI hunch confirmed). Exact
  Seung–Nelson dihedral kernel + 4-vertex force formula (Σf=0, flat→0) + `build_membrane_hinges`. Validate by
  FD-gradient (sign arbiter) + sphere→8πκ + bleb finite-width. Wire membrane_surface=None (replace lumped).

## Open PI decisions (surfaced, do NOT tune to pass)
- IF cage: re-test at strain > 0.34 (plate reaches nucleus)? IF/anchor DENSITY = PI magic number (need the
  cortex-100/µm² analog); over-stiff F_plate 6× = design HALT; buckling (add IF bending?).
- M_biot Biot modulus (1 kPa draft), IF n_fil/ratio_xl/k_LINC, membrane RESERVOIR_STRAIN/k_erm — KB registration.
- KB DRAFTs pending PI: SE_..._biot-coupling-modulus, SE_..._mcf7-cell-nucleus-geometry.

## Constraints (HARD): NATIVE+FULL only (gbook A5000, `~/ff_scratch`, rsync + full python path, log-file
monitor); viz = interactive 3D cell HTML; no-param-tuning; pre-stressed baseline; 한글 답변/영어 코드.
KB gates GREEN (verify_runs 31 / verify_params 42, no drift).
