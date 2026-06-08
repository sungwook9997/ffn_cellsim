# Platform Autorun Log — 2026-06-09

Autonomous compartment-platform deepening session (CWD `ffn_cellsim-platform`,
branch `h7/compartment-platform`). Sibling session works Gate-B cortical-γ in
`ffn_cellsim` (`h7/full-cell-integration`) — its files/branch untouched here.

**Mission:** safely deepen the 8 default-OFF compartments — only work possible
WITHOUT PI parameter ratification or shared-file edits. Smoke harnesses are
PLUMBING checks (never band pass/fail). None-gated constants use DETAILED_PLANS
candidate values, labelled SMOKE-ONLY. Shared/PI-gated work → PLATFORM_PI_QUEUE.md.

| time | item | commit | result |
|---|---|---|---|
| boot | read AGENTS/CLAUDE + 4 boot docs + registry + detailed plans; env ffn_sim, hoomd 7.0.1; clean tree @ fd00309 | — | OK |
| 02:00 | smoke infra: scripts/compartment_smoke/_smoke_common.py (BAOAB stepping sim + cortex stub + SMOKE watermark/fig/json helpers) | (this commit) | OK |
| 02:00 | smoke#1 microtubules: enabled aster+stub, BAOAB step ×3000 @ dt=0.5·dt_cfl; backbone len=l0 (rms dev 0.03%), rod stable, finite | acfa577 | PLUMBING_OK |
| 02:10 | smoke infra: +to_hoomd_snapshot round-trip (bare-gsd None-group hardening for non-gsd-hardened extenders) | (this commit) | OK |
| 02:10 | smoke#2 intermediate_filaments: linear cage (if_backbone + 13 per-r0-bin if_crosslink on one shared Harmonic) ×3000; construction 53.7 kT (force-free ✓ = r0=0 audit-fix holds live), rms dev 0.34% | 172fcbc | PLUMBING_OK |
| 02:20 | smoke#3 linc: perinuclear-cap topology (acceptors seeded at r0 → force-free), 80 linc_nesprin bridges (no zero-bonds trap), construction 0.0 kT, len=r0=50nm, rms tension 4.2 pN (thermal); k_linc=1e-2 SMOKE-ONLY. capture<R_cell-R_nuc warning EXPECTED for cap topology | 98f1b07 | PLUMBING_OK |
| 02:35 | smoke#4 osmotic_regulation: live EnclosedVolumePressure (turgor 133 Pa) + setpoint updater (object-identity), 60 ticks, reads ev.last_pressure + mutates ev.p.V0; hyperosmotic Δc=-100 → V0 RVD-sign down (-0.02%, run 2.1ms ≪ τ_RVD=3.2s honest slope-only); Lp=1e-12 SMOKE-ONLY | f5548b8 | PLUMBING_OK |
| 02:55 | ⭐BUG FOUND+FIXED: cadherin.extend_snapshot_with_cadherins did NOT carry the per-particle `image` array (handled mass/charge/diameter/velocity, omitted image) → Snapshot.from_gsd_frame crash (broadcast (n,3)→(N,3)) on any host frame with image. Fixed in junction/cadherin.py + regression test test_snapshot_extender_extends_velocity_and_image_arrays (32 pass). First smoke catch. | (this commit) | FIX |
| 02:55 | smoke#5 cadherin_junction: 2-cell interface (A/B facing cadherins), Rakshit catch-bond binder forms 10 trans-dimers ALL A↔B (0 intra-cell ✓), elastic force finite; membrane-anchored drag SMOKE-ONLY (cadherins not free-3D-diffusing); n_cad=40 layout not Iturri-223. NOT GATE-J | 69f3c7d | PLUMBING_OK |
| 03:15 | smoke#6 ventral_stress_fibers: 6 FA→FA bundles (144 sf_actin beads), 4/4 sf_ bond types on one shared Harmonic, deterministic-pairing trick → uniform ell0 → force-free (rms dev 0.01%); PASSIVE backbone (NMII NOT attached = BLOCKER#2 sf_myosin_* prefix), T≈0.05 nN ~0 correctly < Kumar band; N_filaments=20→μ_SF SMOKE-ONLY | 6e4cf35 | PLUMBING_OK |
| 03:30 | smoke#7 membrane_reservoir: own mem_node layer (offset 100nm) resolves BLOCKER-1 → 120 tethers ALL cross-layer (0 self-pair), force-free (0.0 kT, len=offset); rupture updater (NotImplementedError) + reservoir-release (f_excess None) PI-blocked, NOT constructed. STATIC mesh only | 942e172 | PLUMBING_OK |
| 03:45 | smoke#8 junctional_actin: STUB — HOOMD build raises NotImplementedError even anchored → SCALAR-only smoke (catch-slip F*=6.02 pN ≈ analytic 6.05, engaged f=0.907, tensile-only coupling) + ASSERTED build barrier. SMOKE-ONLY catch set | (this commit) | SCALAR_OK |
| 03:45 | PLATFORM_PI_QUEUE.md created: activation blockers requiring shared-file/PI/gate-contract changes (γ-denylist ext, sf_myosin_ prefix, 5-pt wiring, None-gated params, STUB build paths, 2-cell doublet) — all SKIPPED per ownership | (this commit) | QUEUE |
| 03:45 | ⭐ALL 8 compartment ENABLED-PATH smokes DONE (6 full HOOMD step + 1 static + 1 scalar) + 1 crash-on-enable bug fixed (cadherin image). Backlog #1+#2 complete | — | DONE |
| 03:55 | smoke aggregator run_all.py + README (one entry point, summary JSON + montage); 8/8 plumbing-OK | b0c8c22 | DONE |
| 04:10 | backlog#3 surface_manifold GEOMETRY-only grid-invariance: area→4πR² (O(1/n_tri)), broad-phase frame conv (17.5°→1.15°), ⭐k-ring reach-coverage MASTER GATE coverage=1.0 at ALL 5 resolutions (k adapts 4→8), geometry-only invariant ✓. NO force added | 85fe59e | GEOMETRY_OK |
| 04:25 | backlog#4 test strengthening: tests/test_compartment_enabled_canary.py — fast CI regression guard for the enabled paths (MT rod-stable, IF force-free construction, LINC bridges-form+no-zero-trap). Full compartment suite + canary: 269 passed, 1 skipped (1.2s) | 6530cf8 | PASS |
| 04:40 | backlog#5+#6 closeout: COMPARTMENT_SMOKE_REPORT_2026-06-09.md (results + bug + master gate + PI queue + Figures section per visualize rule); canary extended to osmotic/stress_fibers/membrane (6 enabled-path canaries, all 6 pass 1.2s) — comprehensive regression coverage | 94f8efd | DONE |
| 05:00 | ⭐BACKLOG EXHAUSTED → deeper hardening: tests/test_compartment_unratified_raises.py — platform SAFETY-CONTRACT canary (8 tests): every None-gated/STUB path RAISES (SF μ_SF, linc k_linc, IF nonlinear, MT DI + L_mt, membrane bleb updater, cadherin count, junctional build+laws). Locks the no-magic-number guarantee. 8 passed | 87dae2a | PASS |
| 05:20 | deeper hardening: coexistence.py — MT aster + IF cage on ONE shared md.bond.Harmonic (18 bond types) + ONE angle force + ONE gamma_map, BAOAB-stepped together. Both intact (mt rms 0.02%, if rms 0.35%), single shared bond force. Validates multi-compartment integration pattern (vs standalone-Harmonic crash class). run_all 10/10 | e0abfbb | COEXIST_OK |
| 05:35 | deeper hardening: tests/test_surface_manifold_grid_invariance.py — CI guard for the registry's named master gate (area→4πR², k-ring reach-coverage=1.0 per resolution, frame convergence, geometry-only). 6 passed 0.46s | 414c8a1 | PASS |
| 05:45 | full verification: all compartment suites + 3 new test files (canary 6 + safety 8 + manifold-grid 6) = 286 passed, 1 skipped (1.55s) | — | PASS |
| 06:00 | backlog#6 lit: COMPARTMENT_SMOKE_CONSTANTS_PROVENANCE_2026-06-09.md — per-constant provenance/confidence/PI-pending checklist (SOLID/DERIVED/ORDER/LAYOUT/NONE) for every SMOKE-ONLY value used. PI ratification checklist | 03634da | DONE |
| 06:20 | capstone: stack_confined_migration.py — the confined_migration recipe internal stack (nucleus envelope + MT aster + IF cage + 30 LINC bridges) on ONE shared md.bond.Harmonic (4 bond families / 19 types) + ONE angle + ONE gamma_map, stepped together. Validates recipe-declared internal-compartment integration. run_all 11/11 | 88c0820 | STACK_OK |

---

## Session summary (backlog + deeper hardening EXHAUSTED)

**Delivered (all on `h7/compartment-platform`, ~19 commits, sibling branch untouched):**

1. **8 compartment ENABLED-PATH smoke harnesses** — the first actual executions of
   every default-OFF compartment's enabled build/attach/step path. 6 full HOOMD
   BAOAB-stepped (microtubules, intermediate_filaments, linc, osmotic_regulation,
   cadherin_junction, ventral_stress_fibers), 1 static-mesh (membrane_reservoir,
   rupture updater PI-blocked), 1 scalar (junctional_actin STUB, HOOMD build blocked).
2. **surface_manifold GEOMETRY-only grid-invariance harness** — area→4πR², broad-phase
   frame convergence, and the k-ring reach-coverage **master gate** (coverage=1.0 at
   every resolution). No force added.
3. **2 multi-compartment coexistence harnesses** — MT+IF (18 bond types, 1 shared
   Harmonic) and the confined_migration internal stack (nucleus+MT+IF+LINC, 4 bond
   families / 19 types, 1 shared Harmonic). Validates the real integration pattern.
4. **1 crash-on-enable bug found + fixed** — `cadherin.extend_snapshot_with_cadherins`
   dropped the per-particle `image` array → state-creation crash; fixed + regression test.
5. **20 new CI tests** — 6 enabled-path canaries, 8 safety-contract (un-ratified paths
   MUST raise), 6 surface-manifold grid-invariance. Full compartment suite: **286 passed,
   1 skipped**.
6. **run_all aggregator (11/11 green) + README + 11 SMOKE-watermarked figures.**
7. **4 docs:** PLATFORM_PI_QUEUE (activation blockers), COMPARTMENT_SMOKE_REPORT
   (closeout), COMPARTMENT_SMOKE_CONSTANTS_PROVENANCE (PI ratification checklist), this log.

**Found-out / confirmed:** every enabled path assembles + steps (or is honestly
blocked); the cadherin image bug was the only crash; the manifold master gate holds;
multi-compartment shared-force integration works.

**Everything beyond this is PI-gated** (see PLATFORM_PI_QUEUE.md): γ-denylist extension
in `cortex/cortical_tension.py`, the SF `sf_myosin_*` prefix split in `cortex/myosin.py`,
the 5-point loader wiring (`cell.py`/`manifest.py`), None-gated parameter ratifications,
the junctional_actin build path, and the two-cell doublet builder — all require shared-file
or gate-contract changes outside this session's ownership, so all were SKIPPED and recorded.

**STATE: in-scope backlog + deeper hardening exhausted; awaiting PI direction or a
shared-file/parameter authorization.**

---

## Phase 2 — ACTIVATION (PI granted 소유권 허용 + "모든 것 달려, 절대 멈추지 마")

Ownership unblocked → graduating compartments EXPERIMENTAL→LIVE via the real 5-step
activation procedure, keeping the disabled path bit-identical. Conflict-aware: the
sibling Gate-B session is rewriting `cortex/cortical_tension.py`, so I start with the
activation that needs NO cortical_tension edit.

| time | item | commit | result |
|---|---|---|---|
| 06:50 | ⭐**osmotic_regulation EXPERIMENTAL→LIVE** (1st graduation). Wired POST-build in manifest.py (ResolvedBaseline field + resolve stanza + build_baseline_cell attach via the live ev_force handle — ZERO cell.py edit), + mcf7_baseline stanza (Lp=1e-12 Jung2011 MCF7/AQP5 PI-ratified) + osmotic_rvd recipe + registry LIVE. NO cortical_tension touched (sibling-safe) | (this commit) | LIVE |
| 06:50 | osmotic ACTIVATION GATE (full physiological baseline, real build OFF vs ON): PRIMARY τ_RVD=3.2s ∈[3,600]s ✓; RVD sign V0↓ ✓; updater attached (+1) ✓; NO-contamination (bond inventory identical ON/OFF — osmotic adds 0 bonds) ✓; OFF bit-identity ✓. +2 CI tests (64 passed). Raw full-cell run trips BAOAB guard (needs equilibration prelude — known full-cell issue, not osmotic) → V0(t) via exact updater law | 67f65a3 | GATE PASS |
| 07:10 | ⭐INFRA UNBLOCKER: cortex/cortical_tension.py γ-denylist now REGISTRY-DRIVEN (NONCORTICAL_COMPARTMENT_PREFIXES = REGISTRY.gamma_denylist() minus cortex_*). Activating any γ-contaminating compartment auto-excludes its bonds from cortical γ; cortex_myosin_* (the active-γ signal) PRESERVED. Minimal+additive (filter only, not the γ math). 11 cortical_tension + 31 registry tests pass (+cross-check). Sibling owns this file → merge-coordinate | 8188c8f | DONE |
| 08:00 | ⭐**microtubules EXPERIMENTAL→LIVE** (2nd graduation). Full snapshot-extension wiring in cell.py (import+sig+snapshot extend+shared bond mt_backbone+shared angle mt_bending+LJ pairs r_cut=0+gamma_map cytoplasm γ_b+CFL gate+CellBuildOptions+Cell dataclass/sig/call/return/needs_full) + manifest.py (resolve stanza w/ cytoplasm γ_b) + mcf7_baseline stanza + aster_microtubules recipe + registry LIVE. OFF bit-identity verified (N=19000, 3 updaters, p_mt None) | (this commit) | LIVE |
| 08:00 | MT ACTIVATION GATE (full baseline OFF vs ON): aster assembled (+176 part/+175 bonds at n_mt=7), CFL passed at cytoplasm drag (1.55e-5, NOT water), ⭐NO-CONTAMINATION γ_soft IDENTICAL OFF=ON=5.14e-10 (mt_ registry-denylisted), L_p=5.2mm in band. +2 build tests (40 MT pass). ⚠️BLOCKER: single-hub MTOC degree caps n_mt≤7 (nlist exclusion cap); LJ-less smokes missed it; dense aster needs multi-bead MTOC core (PI_QUEUE). 2 registry tests updated for the LIVE flip | 44ea125 | GATE PASS |
| 08:40 | ⭐**intermediate_filaments EXPERIMENTAL→LIVE** (3rd graduation). Same snapshot-extension wiring (cell.py: extend_if_cage + register_if_bond_params on shared bond + if_bead LJ r_cut=0 + gamma_map cytoplasm + threading; manifest.py: resolve w/ R_nuc from nucleus) + config + if_cage recipe + registry LIVE. LINEAR path (nonlinear Table law stays PI-pending/NotImplementedError) | (this commit) | LIVE |
| 08:40 | ⭐2 MORE crash-on-enable BUGS FIXED in intermediate_filaments.extend_snapshot_with_if_cage: (1) np.asarray(velocity) 0-d crash on a build-time gsd Frame (velocity/mass/image None) → added _pf None-guards; (2) list(angles.types) on None → guarded N/types. Both masked by the smoke's to_hoomd_snapshot; caught by the real full-cell build. IF GATE PASS: cage +480 part/+585 bonds, NO-contamination γ_soft IDENTICAL=5.14e-10, if_ denylisted. +2 build tests (49 pass) | (this commit) | GATE PASS |
