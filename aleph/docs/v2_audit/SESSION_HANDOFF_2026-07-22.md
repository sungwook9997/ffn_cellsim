# Session handoff — ac/ engine (2026-07-22) → next Claude session

Branch `codex/ff-ac-codex` (shared with a Codex/parallel session). Everything below is COMMITTED. Read this
+ the two roadmap docs, then continue. **You are the Lead** (PI directive 2026-07-22): own the work end-to-end;
Codex aligns to you and the PI decides when Codex moves — do NOT defer solver work to Codex, but be aware Codex
actively edits `implicit_mechanics.py` (see the gbook-sync note).

## What this session delivered (all committed + A5000-verified)

1. **3 new isolated compartments** (reuse ff/ physics, §1.4 `accumulate` contract, NOT yet composed into
   build_cell — gated on convergence): `aleph/components/solid/microtubule.py` (MT strut, 5/5), `intermediate_filament.py`
   (keratin/vimentin cage, nucleus-coupling, 3/3), `adhesion_clutch.py` (α2β1 slip + Chan-Odde κ_c, 4/4).
   Full ac/ suite **499 passed** on the A5000.
2. **KB (Notion SoT)**: KB-3.24–3.33 (MT/IF/nucleus) + KB-1.31–1.32 (collagen ECM) + 13 SourceEvidence,
   DOI-verified. Param dossier: `docs/v2_audit/_historical/COMPARTMENT_PARAM_SOURCING_2026-07-22.md`.
3. **⭐ Convergence failure-mode #1 SOLVED — overlap-free cortex.** The zero-thickness shell packed 70,686
   filaments at exactly R → ~64k WCA interpenetrations → steric-only max **2634 pN** (the artifact that
   dominated the resting-convergence residual). Fix: `CortexParams.cortex_thickness_um` (radial dispersion) +
   `resolve_overlaps` (soft-sphere relaxation) → **build starts overlap-free, steric 2634→0 pN** (verified
   full-native on the A5000), crosslinks force-free (rest = resolved positions). Deployed opt-in through
   `ff.weave.weave(overlap_free=)` → `CellConfig.overlap_free_cortex` → driver **`--overlap-free-cortex`**.
   Default OFF (γ-parity untouched). With it, full-native `--from-resting` residual dropped 2634 → the
   no-steric coupled-preload regime.
4. **Membrane/nucleus resolution ratified (PI 2026-07-22)**: `CellConfig.membrane_subdivisions /
   nucleus_subdivisions` + driver `--membrane-subdivisions / --nucleus-subdivisions`. **subdiv 6 = 40,962
   verts (~131nm) is the dynamic BASELINE**, **subdiv 7 = 163,842 (~66nm) the VALIDATION** resolution; actual
   dynamic runs may later use adaptive bleb-site remeshing. Default 3 (resting). Higher also densifies the
   per-node ERM coupling (642 → 40,962). Full-native: subdiv 6 = 591,754 nodes / 592 MB / 13 s; subdiv 7 =
   837,514 / 752 MB / 34 s (A5000 16 GB → memory is a non-issue).
5. **⭐ Convergence failure-mode #2 STARTED — membrane end.** The analytic operator lumped membrane/nucleus/
   pressure stiffness into ONE scalar regularizer, so it couldn't see the coupled membrane→ERM→fiber mode.
   Added the membrane area-tension (γ_mem) as explicit **edge springs** in `ProjectedAnalyticCG`
   (`implicit_mechanics.py`), removed γ_mem from `omitted_regularization_base`. Verified on A5000: a 1-node
   membrane perturbation now couples to its ~5 neighbours through the operator (was 0). The membrane is now a
   connected taut sheet in the operator.

## THE remaining blocker — convergence #2 (the ONLY gate to dynamic sims)

`--from-resting` does not yet reach `outer_accepted=True`. #1 is solved; #2 (coupled-preload conditioning,
~0.34 pN no-steric projected-force stall) still blocks acceptance. Residual decomposition (Codex-verified):
per-fiber translation 35% + fiber internal 44% + global radial 17% + membrane-ERM star 22%.
- **Membrane end (~39%)**: DONE this session (area-tension edge springs). Follow-up: the membrane **Helfrich
  bending** tangent (hinge 4-node Hessian — still in the scalar regularizer) + nucleus tangents.
- **Fiber conditioning (~79%, DOMINANT)**: the real remaining work — a stronger coupled multilevel/coarse
  preconditioner (overlapping Schwarz over membrane-node⊕ERM⊕fibers + Schur/coarse solve + local radial/shell
  coarse modes l=0/1/2). This is genuine multigrid/domain-decomposition engineering. Priority per the
  diagnosis: ① exact membrane/pressure tangent (started) → ② ERM-fiber-crosslink multilevel coupled block →
  ③ radial/shell coarse modes → ④ semismooth active-set. NO physical pins / regularizer inflation / tolerance
  loosening (the residual is real).
- **4-condition GO signal** (all in ONE clean-gbook run): interpenetration 0 (✅ via --overlap-free-cortex) +
  strict projected-force gate (0.21 pN full-native) + inner_converged=True + outer_accepted=True/no rollback.

## Critical operational notes
- **⚠️ gbook is PARTIAL-SYNCED**: `ffn_ac_native` has builder + `implicit_mechanics.py` = Lead branch, but
  `driver.py` + `inner_mechanics.py` = the `ac/wca-fix1-linear-core` 30 pN side-branch. So `--from-resting`
  ACCEPTANCE on gbook is unreliable (wrong gate); operator/steric/build MEASUREMENTS are reliable. Codex will
  send **`SYNC NOW @ <commit>`** with file hashes when production changes merge + verify — only run the
  authoritative GO probe on a single-commit-aligned gbook. The 30 pN branch is EXCLUDED from GO (142× looser).
- **Codex protocol**: PI controls when Codex moves; Codex aligns to Lead. Codex pings when packing or #2
  pieces land. Don't edit shared files assuming Codex won't — but as Lead your changes are authoritative.
- **gbook**: ssh `gbook`, python `~/miniconda3/envs/ffn_sim/bin/python`, repo `~/ffn_ac_native`, syncs via
  Syncthing (rsync ff/+ac/ before a run). Dev Mac = CPU-only Warp (kernels source-only; NumPy refs gate).
- **Background run `bqdo6iq4k` — DONE, INCONCLUSIVE (partial-sync artifact):** it reported
  `residual_start == candidate == end == 47.405` **Pa** unchanged, `inner_converged=False`. That 47 Pa is the
  ERM force-family metric emitted by the side-branch driver/inner_mechanics still on gbook — NOT the projected
  force gate my membrane tangent targets. Lesson: **the membrane-tangent effect cannot be measured until gbook
  is single-commit-aligned** (Codex `SYNC NOW`). Re-run this exact probe post-sync to get the real number.

## Next steps (prioritized)
1. **#2 fiber conditioning** (dominant): implement the coupled multilevel/coarse preconditioner (see priority
   ①②③④ above). This is the gate. Verify operator-coupling in isolation on gbook (like the membrane tangent),
   then the full GO probe on a synced gbook.
2. On the **4-condition GO** (clean gbook): compose MT/IF/ECM into build_cell (add `with_microtubules/with_if/
   with_ecm` flags, thread like overlap_free_cortex) + re-run **bleb growth** (`aleph/components/incumbent/bleb_growth.py`,
   already hardened) at `--overlap-free-cortex --membrane-subdivisions 6`.
3. Fidelity endgame (parallel, ungated): adaptive membrane remeshing at bleb sites; sourced ERM density
   (currently one-per-node); compose the PARTIAL wirings (chromatin WLC, RAD monomer, dynamic crosslinker KMC).
4. Filopodium/SF/lamellipodium are EMERGENT regions (weave more regions, not scripted) — composition-gated.

## Key files
- Convergence: `aleph/components/incumbent/implicit_mechanics.py` (operator + preconditioner), `aleph/components/incumbent/inner_mechanics.py`
  (convergence gate), `aleph/laws/cortex_assembly.py` (overlap-free builder).
- Compartments: `ac/solid/{microtubule,intermediate_filament,adhesion_clutch}.py`.
- Bleb validation: `ac/cell/{bleb_perturbation,bleb_experiment,bleb_growth}.py`,
  `docs/v2_audit/{BLEBBING_VALIDATION_PLAN,COMPARTMENT_COMPLETION_ROADMAP,COMPARTMENT_PARAM_SOURCING}_2026-07-2x.md`.
- Cell viz: `scripts/ac_cell_assembled_viz.py` + `dump_state.py` + `browser_check.py`.
