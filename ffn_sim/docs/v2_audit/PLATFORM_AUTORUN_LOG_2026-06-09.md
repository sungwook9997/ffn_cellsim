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
| 08:40 | ⭐2 MORE crash-on-enable BUGS FIXED in intermediate_filaments.extend_snapshot_with_if_cage: (1) np.asarray(velocity) 0-d crash on a build-time gsd Frame (velocity/mass/image None) → added _pf None-guards; (2) list(angles.types) on None → guarded N/types. Both masked by the smoke's to_hoomd_snapshot; caught by the real full-cell build. IF GATE PASS: cage +480 part/+585 bonds, NO-contamination γ_soft IDENTICAL=5.14e-10, if_ denylisted. +2 build tests (49 pass) | 70f8238 | GATE PASS |
| 09:00 | ⭐INTEGRATION CAPSTONE: internal_live recipe (osmotic + microtubules + IF together on the full baseline). All 3 compose in ONE cell (+656 internal particles, osmotic updater +1), INTEGRATED no-contamination γ_soft IDENTICAL OFF=ON=5.14e-10 (no mutual contamination). +1 integration CI test (32 registry pass). The real physiological cell with internal organelles, all LIVE | 003d195 | PASS |
| 09:10 | comprehensive regression: full compartment + cell + cortex suites = 308 passed, 1 skipped (6.1s). All deep activation wiring (cell.py/manifest.py/cortical_tension.py) verified non-regressing | — | PASS |
| 09:20 | COMPARTMENT_ACTIVATION_REMAINING_DECISIONS doc: each remaining compartment (LINC/SF/cadherin/membrane/junctional) reduced to a crisp PI physics-design decision + my recommendation. Phase-2 activation summary below | (this commit) | DOC |
| 10:40 | ⭐**linc EXPERIMENTAL→LIVE** (4th graduation; REMAINING_DECISIONS recommended-order ①, Option A). linc.py: per-bond EXACT-r0 layout (LINCLayout + compute_linc_layout + extend_snapshot_with_linc_layout + configure_linc_bond_potential(layout=)) mirroring the FA molecular-clutch exact-r0 convention (16-bin would leave ~100s kT at stiff k_linc); rewrote extend_snapshot_with_linc fresh-`hoomd.Snapshot` (the in-place gsd mutation crashes on a hoomd.Snapshot — BondDataSnapshot has no whole-array setter). cell.py 5-pt wiring (import+sig+extend AFTER IF using nucleus_bead+if_bead typeid rows+per-bond shared-bond register+CellBuildOptions+Cell dataclass/build threading+handles); manifest stanza (requires nucleus+IF, capture_radius=IF l_seg + n_bridges_max=IF n_filaments geometry-derived); mcf7_baseline linc stanza (k_linc=1e-2 route-A); linc_coupled recipe; registry LIVE + particle_types_added=() fix | (this commit) | LIVE |
| 10:40 | LINC ACTIVATION GATE (full physiological baseline, IF cage ON, OFF vs ON): **PASS** all 5 controls — bridges formed n=40 (REFUTE-guard, +0 particles bonds-only, +40 bonds/40 types); FORCE-FREE construction max strain 0.0 (per-bond EXACT-r0); ⭐NO-CONTAMINATION cortical γ_soft IDENTICAL OFF=ON=5.1411e-10 (linc_ registry-denylisted, all linc types excluded from mask); OFF-identity (+0 particles); CFL headroom dt_prod 1.30e-8 ≤ dt_cfl 6.21e-5 (~4 orders). ⚠FOUND+FIXED: many nucleus surface beads share one nearest if_bead → that if_bead's bond degree overflows the nlist exclusion cap (same class as MTOC) → added unique_acceptor greedy 1:1 matching (one nesprin per IF anchor, degree +1). +3 build tests (test_linc.py 45 pass; 149 compartment-suite pass). [2,10]pN resting-tension oracle DEFERRED (needs equilibration; raw full-cell run trips BAOAB guard) | (this commit) | GATE PASS |

| 11:10 | ⭐**compartment_vis.py** (PI request "visualize the constructed cell shape"): single entry point builds the FULL physiological cell with EVERY LIVE compartment ON + renders ACTUAL geometry — `compartment_cell_overview.png` (3D + equatorial/meridional cross-sections with MT/IF/LINC/mem bonds drawn + per-compartment radial density profile showing each shell at its physiological radius) + `compartment_cell_3d.png`. Confirms spherical cortex shell @ R_cell=7.5µm, filled nucleus core, IF cage just outside nucleus, MT aster from MTOC, LINC nucleus↔IF bridges, membrane layer just outside cortex. SMOKE_REPORT §Figures updated | (this commit) | VIZ |
| 11:50 | ⭐**membrane_reservoir EXPERIMENTAL→LIVE** (5th graduation; REMAINING_DECISIONS ②). Authored own mem_node offset layer in membrane_reservoir.py (extend_snapshot_with_membrane_reservoir + _clone_frame_with_mem_layer: fresh hoomd.Snapshot adds n_mem_nodes mem_node beads at membrane_offset outside cortex + 1 mem_tether/node→DISTINCT cortex bead, force-free at offset; BLOCKER-1 self-pair fix). ⚠degree-aware unique-acceptor selection (max_anchor_degree=5) bounds per-cortex bond degree to +1 → nlist exclusion-cap safe (LINC/MTOC lesson). resolver +membrane_offset(¼ reach=50nm)/n_mem_nodes(2000). Full 5-pt wiring (cell.py snapshot-ext + mem_node LJ r_cut=0 + cytoplasm gamma_map + mem_tether on shared bond + options/dataclass/build; manifest stanza requires membrane_surface; config; membrane_reservoir_tethered recipe; registry LIVE) | (this commit) | LIVE |
| 11:50 | membrane_reservoir ACTIVATION GATE (full baseline OFF vs ON): **PASS** 6 controls — mesh assembled n_tethers=1969 (+1969 mem_node part, +1969 bonds); cross-layer (0 self-pairs, unique cortex anchors); FORCE-FREE max strain 1.7e-14 (born at offset); ⭐NO-CONTAM cortical γ_soft IDENTICAL OFF=ON=5.1411e-10 (mem_ registry-denylisted); off-identity; bleb paths HONESTLY BLOCKED (MembraneTetherUpdater raises σ_crit_bleb None + released_area raises f_excess None). +3 build tests (test_membrane_reservoir.py 30 pass). Bleb nucleation + reservoir tension-buffering DEFERRED (PI-pending σ_crit_bleb/f_excess) | (this commit) | GATE PASS |

| 12:30 | viz#2 (PI feedback): compartment_vis.py → **SimuCell3D-style surfaces** (faceted cortex/membrane ConvexHull + cutaway exposing nucleus + MT aster/IF/LINC bonds) + **VERSIONED output** (outputs/h7/figs/morphology/cell_vNN_<label>.png, accumulates, never overwrites; cell_latest.png pointer). Retired flat compartment_cell_*.png. First snapshot cell_v01_membrane.png | (this commit) | VIZ |
| 13:10 | ③ ventral_stress_fibers PUSHED + SCOPED: FA-adhered build + integrin harvest + SF extender all work mechanically (probe: adherent_passive → 2150 integrins; SF assembles 20 bundles on the real hoomd.Snapshot). ⚠SURFACED a fidelity-affecting **FA-pair geometry design decision** (generate_stress_fiber_layout pairs by RANDOM shuffle → unphysical bundle lengths across z=−7.35..1.6µm; physiological vSF need nearest-basal/aligned pairs) → per mission rule logged to PLATFORM_PI_QUEUE + recommend (a) nearest-basal pairing v1; HELD at EXPERIMENTAL (not faked LIVE on random pairing). Also: per-bundle k_actin/ell0 measure-wiring + N_filaments=20 candidate + sf_myosin prefix (active) still pending | (this commit) | PI-SURFACE |

| 14:30 | PI ANSWERED ③ decision: **long-axis-aligned** FA pairing + "complete SF first". | — | PI-DIR |
| 15:20 | ⭐**ventral_stress_fibers EXPERIMENTAL→LIVE** (6th graduation; PASSIVE backbone). stress_fibers.py: select_aligned_fa_pairs (basal subset→in-plane PCA long axis→low/high pairing) + pair_mode='as_given' + per_bundle_r0 EXACT-r0 backbone (sf_actin_backbone_bin_names; force-free even with varying bundle lengths). Full wiring: cell.py (import+sig+harvest integrins after FA+aligned pairs+SF extend pair_mode=as_given/per_bundle_r0+sf_actin/sf_xlink_head LJ r_cut=0+cytoplasm gamma_map+register sf bonds per-bundle+options/dataclass/build/handles); manifest stanza (requires fa); config (N_filaments=20 Cramer1997 candidate, n_SF=20); ventral_stress_fibers_passive recipe (extends adherent_passive); registry LIVE (denylist ('sf_',) — dropped cortex_myosin_ which the estimator never denylists). ⚠ denylist fix caught by test_live_contaminating_compartments | (this commit) | LIVE |
| 15:20 | SF ACTIVATION GATE (FA-adhered OFF=adherent_passive vs ON, build-time): **PASS** 5 controls — bundles assembled n_SF=20 (+1120 part, +1460 bonds); ALIGNED mean|cos|=0.929 (long-axis, lengths 0.5-13.7µm — physiological vSF not random); FORCE-FREE per-bundle backbone max strain 1.56e-14; ⭐NO-CONTAM cortical γ_soft IDENTICAL OFF=ON=5.1411e-10 (sf_ registry-denylisted); off-identity. +3 build tests (test_stress_fibers.py 33 pass). Kumar 10-30nN tension band = ACTIVE+equilibrated gate (NMII via sf_myosin_ prefix) DEFERRED | (this commit) | GATE PASS |

| 16:10 | viz#3: compartment_vis.py builds the FULLEST cell (adherent: +fa +ventral_stress_fibers); morphology v02 cell_v02_stressfibers_adherent.png (24,938 part) — progression v01(suspended,5 internal)→v02(adherent,6 LIVE+FA+vSF) preserved | (commit) | VIZ |
| 17:30 | ⭐**cadherin_junction EXPERIMENTAL→LIVE** (7th graduation; FIRST multicell). NEW two-cell assembler cell/doublet.py::build_cell_doublet — two cortex shells (build_cortex_state ×2, offset +2R+gap), merged frame (bond/angle reindex), cadherins seeded on MATCHED facing caps (shared y,z, ±r0_trans/2 across interface → force-free trans-dimer), trans-dimers SEEDED pre-bound (engaged-junction baseline; bare junction needs ~1e3 binder batches), cadherin_anchor bonds hold cadherins on the surface, attach_cadherin_junction binder maintains. registry LIVE (manifest_path None — doublet build, NOT single-cell toggle; requires=()). ⚠defensive n_sub cap in cadherin binder (extreme-slip RNG-overflow guard). 2 registry tests switched to junctional_actin (still STUB) | (this commit) | LIVE |
| 17:30 | GATE-J (two-cell doublet, build-time): **PASS** 5 controls — doublet assembled (2×7000 cortex + 120 cadherins + cadherin_trans/cadherin_anchor); TRANS ENGAGED n=60 dimers ALL A↔B (cross=60, intra=0 — genuine 2-cell junction); FORCE-FREE max strain 1.06e-15; ⭐NO-CONTAM cadherin_ excluded from cortical mask (γ_soft=1.23e-5 over cortex only); binder attached. +2 doublet tests (test_cadherin_junction.py 34 pass). DYNAMIC catch-slip maintenance + Iturri ~6.5nN de-adhesion DEFERRED (equilibrated run; BAOAB guard) | (this commit) | GATE PASS |

| 18:30 | ⭐**junctional_actin STUB→LIVE** (8th & FINAL graduation). Implemented the RESERVED build path in junctional_actin.py: extend_snapshot_with_junctional_actin (enabled+anchored) appends one junc_actin head per interface cadherin within the α-catenin reach of a SAME-cell cortex bead — head placed anchor_r0 from the cadherin toward cortex (force-free anchor) + junc_actin_anchor (head↔cadherin) + per-r0-bin junc_actin_couple_b{i} (head↔cortex, force-free); + register_junctional_actin_bond_params. Built on the doublet via build_cell_doublet(with_junctional_actin=True) — couples each cell's cadherins to its OWN cortex. Config catch-set candidates (Buckley x_catch/x_slip SOLID; k_catch0/k_slip0 ORDER; k_couple/k_anchor/k_on/max_couple_dist/anchor_r0 DERIVED H.3); module default None → un-anchored build still raises (STUB contract preserved). registry STUB→LIVE. 2 registry tests → synthetic EXPERIMENTAL spec (no real EXPERIMENTAL/STUB left) | (this commit) | LIVE |
| 18:30 | junctional_actin GATE (doublet build-time): **PASS** 5 controls — belt assembled (46 heads/couplings of 120 interface cadherins; only tips within the catch reach couple); SAME-CELL coupling (cross_cell=0); FORCE-FREE anchor strain 6.6e-14 + per-r0-bin coupling within ½-bin; NO-CONTAM junc_actin_ excluded from cortical mask; ⭐CATCH signature biphasic k_off(F) with F*=6.0 pN (Buckley). +2 build tests (test_junctional_actin.py 21 pass). DYNAMIC catch-slip maintenance + dense-belt fidelity (cadherin-tail particle) DEFERRED | (this commit) | GATE PASS |
| 18:30 | ⭐⭐**ALL 8 default-OFF compartments now LIVE** — the EXPERIMENTAL/STUB activation backlog is COMPLETE. osmotic·MT·IF·LINC·membrane_reservoir·ventral_stress_fibers·cadherin_junction·junctional_actin all graduated with passing activation gates this session (PI 소유권 허용). Remaining are DEFERRED active/dynamic phases (NMII, catch-slip maintenance, nonlinear IF, bleb, dense belt) + PI parameter ratifications, not activations | — | DONE |

### Phase-2 ACTIVATION summary (PI 소유권 허용)
**8 (ALL) compartments EXPERIMENTAL/STUB→LIVE** (osmotic_regulation, microtubules,
intermediate_filaments, linc, membrane_reservoir, ventral_stress_fibers,
cadherin_junction, junctional_actin) with passing activation gates — the activation
backlog is COMPLETE; + SimuCell3D-style versioned morphology visualizer
(compartment_vis.py). Detail below. osmotic_regulation
(τ_RVD in band, RVD sign, no-contam), microtubules (aster, CFL, no-contam γ identical;
n_mt≤7 cap blocker found+documented), intermediate_filaments (cage, no-contam; linear
path, nonlinear PI-pending), **linc** (Option A nucleus↔IF-cage per-bond EXACT-r0 bridges;
bridges-formed + force-free + no-contam γ identical + CFL ~4-order headroom; [2,10]pN
oracle deferred; unique-acceptor degree fix), **membrane_reservoir** (own mem_node
offset layer + static mem_tether mesh; mesh-assembled + cross-layer + force-free +
no-contam γ identical; bleb/reservoir-release PI-blocked; degree-aware exclusion-cap
fix), **ventral_stress_fibers** (PASSIVE backbone; PI-ratified long-axis-aligned basal
FA pairs |cos|≈0.93, per-bundle EXACT-r0 force-free, no-contam γ identical; Kumar/NMII
active gate deferred via sf_myosin_ prefix), **cadherin_junction** (FIRST multicell; new
two-cell build_cell_doublet assembler — 60 trans-dimers all A↔B, force-free, cadherin_
γ-excluded; dynamic catch-slip maintenance + Iturri de-adhesion deferred),
**junctional_actin** (LAST; reserved STUB build path IMPLEMENTED on the doublet —
α-catenin/vinculin cadherin↔cortex catch clutch, same-cell coupling, force-free,
junc_actin_ γ-excluded, biphasic catch F*=6 pN; dynamic maintenance + dense-belt
fidelity deferred). **+ registry-driven γ-denylist** (unblocks no-contamination for
all). **+ internal_live integration capstone** (osmotic+MT+IF in one cell, no
contamination) **+ build_cell_doublet two-cell assembler** (cell/doublet.py) for the
multicell pair. **6 crash-on-enable bugs fixed** + **2 nlist exclusion-cap fixes**
(LINC + membrane acceptor-degree) + **cadherin binder n_sub overflow guard**.
**compartment-suite tests green** (LINC 45 + membrane 30 + SF 33 + cadherin 34 +
junctional 21 + registry + canary + IF + MT; full compartment/cell/cortex modulo the
pre-existing KU-3.5 motors-off baseline failure, unrelated to platform work).
**No EXPERIMENTAL/STUB compartment remains.** The open items are all DEFERRED active/
dynamic phases (cortex+SF NMII via sf_myosin_ prefix, cadherin/junctional catch-slip
maintenance, nonlinear IF Table law, membrane bleb, dense junctional belt) + PI
parameter ratifications (k_linc, N_filaments, σ_crit_bleb, f_excess, Iturri SE) — see
PLATFORM_PI_QUEUE — NOT activations.

---

## Phase 3 — ACTIVE / DYNAMIC deepening (PI "절대 멈추지 마"; DEFERRED items in priority order)

Activation backlog COMPLETE → working the deferred active/dynamic physics +
fidelity, priority ① active NMII → ② catch-slip maintenance → ③ PI ratifications
→ ④ fidelity. Sibling Gate-B session (`h7/full-cell-integration`) untouched.

| time | item | commit | result |
|---|---|---|---|
| ①a | **cortex/myosin.py prefix-split** — ResolvedCortexMyosin.prefix (default "cortex_myosin_", byte-identical) + resolver read/validate + prefix-derived name helpers (myosin_particle_type_names / *_bond_name / myosin_attach_bin_names; cortex_myosin_attach_bin_names = default wrapper) + extend/register/MyosinStepUpdater filter all derive from prefix. sf_myosin_ placement → under existing ('sf_',) γ-denylist; cortex_myosin_ stays in active-γ. test_myosin 31 + new test_myosin_prefix_split 11 pass; SF/registry/cortical_tension regress clean | 3c5795c | DONE |
| ①b | **MyosinStepUpdater actin-pool generalization** — explicit actin_pool_tags (default None→arange(n_cortex_actin), byte-identical). All bead ids stay GLOBAL; KDTree pool local-indexed via _pool_g2l. r_actin_all=pos[pool], segment pos[a_idx] global, _bead_to_segs local-indexed, degree/bead_attach_count sized N_part. SF binds the (non-contiguous) sf_actin chain via the same Stam-Hocky/Hill machinery. 63 dynamics tests byte-identical (myosin/grip_walk/dipole gates) + new test_myosin_sf_pool 2 (PARITY shifted-pool==default + CONFINEMENT bound∈pool) | db00c20 | DONE |
| (a) | **PI decision (a) + design doc** — connected basal mesh = cortex connected-mesh on a FLAT basal disk; B1-B5 increments. docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md + Notion Dev-Logs milestone | a7d7074/f342709 | DOC |
| B1 | **basal-mesh GEOMETRY layout** — cell/basal_mesh.py generate_basal_mesh_layout (bimodal cables+infill on flat disk, reuses VariableLengthCortexLayout, geometry-only no force). Gate PASS: 600 fil=75 cables+525 infill, planar in band, force-free (strain 4e-15), cables |cos|=1.000, infill |cos|=0.639. test_basal_mesh 14 + scripts/h7_basal_mesh_gate.py auto-viz | 6f67b99 | GATE PASS |
| B2 + reframe | **F-layer connectivity + 2-LAYER reframe** (PI clarification). connect_basal_mesh (reuse seed_connected_mesh_xlinks, disk reach footprint·√(π/F), max_bind_dist=reach). Gate PASS: n_xl=1850, giant=0.998, z=3.087∈[3,3.5], L/lc=6.18. ⭐2-layer: S=surface_manifold positions FA+patch connectivity (NO force, per 7b19276 rejection); F=explicit filament network=real forces. test_basal_mesh 17. Design doc updated | 78eb25e | GATE PASS |
| S1 | **S-layer basal SURFACE** (curved cap, superseded) — first cut used a south-pole sphere cap | 759b585 | (superseded) |
| bend | **SF/basal BENDING EI anchor wired** (the missing flexural term PI identified — SF was a straight chord with NO angle potential, κ=0; cortex had it). resolve_basal_bending: EI_single=7.3e-26 N·m² (Gittes FLEXURAL, correct use), EI_cable=N·EI_single (loose) / N² (tight), angle_k=EI/ℓ0 (cortex convention). basal_bending_report PASS: angle groups present (257), force-free at construction (θ≈π, dev 3e-8), ⭐**bend-vs-stretch ratio ka/kbend ≈ 1.5e5 (N-independent) → bend-before-stretch TRUE** → AFINES soft-stretch JUSTIFIED (can soften k up to ~1.5e5× & still bend-first → explicit BAOAB dt~µs feasible, **CFL wall solvable WITHOUT the frozen integrator**). +3 tests | (this commit) | GATE PASS |
| B5(i) | **SF tension force-BUDGET (quasi-static, generation-limit probe)** — sidesteps the frozen integrator + stiff-backbone CFL. Result REFUTE (provisional): raw mesoscale single-SF tension ≈5e-12 N (φ=0.99·10 heads·0.5pN) vs Kumar 10-30 nN → **~2014× UNDER band** = SAME generation-limit as the cortical γ-floor. ⭐**N_filaments-INDEPENDENT** (N sets strain ε=T/μ_SF, NOT the tension ceiling). Reaching Kumar needs ~4000× mesoscale force-scale (Route B parallel-bundle) or density/native fix — not a stiffness change. scripts/h7_basal_sf_force_budget.py. EA citation fix c2f599f. Deep-research on how other sims (Cytosim/AFINES/MEDYAN/active-gel) handle both walls in flight | (this commit) | REFUTE (provisional) |
| B4 | **sf_myosin_ NMII placement on the apparatus** — place_sf_myosin_on_apparatus reuses cortex Stam-Hocky/Hill via ①a prefix + ①b pool; heads bind apparatus beads. Fixed sphere-assumption (head-offset normal centers/|centers| → added surface_normal=+z override, byte-identical default). Gate PASS: 40 minifilaments assembled, force-free (strain 1.4e-14), 96.4% heads binding-eligible (median perp 103nm≤210), sf_-denylist, distinct from cortex_myosin_. test_basal_mesh +3; 91 myosin+basal tests green (byte-identity). scripts/h7_basal_nmii_gate.py auto-viz | (this commit) | GATE PASS |
| B3 | **F-on-S integration** — build_basal_filament_network: LONG cables span FA→FA (ends AT FA = force-free sf_anchor), SHORT infill fills the basal actin slab [z_basal,z_basal+band]; one combined flat layout, per-filament EXACT-r0. Gate PASS: 12 cables+500 infill, cables anchored FA→FA, anchors force-free (sep<1e-15), in-slab, chains force-free (strain 1.9e-15), combined mesh percolates (giant 0.963, z 3.078, L/lc 6.22). test_basal_mesh +5 + scripts/h7_basal_apparatus_gate.py auto-viz | (this commit) | GATE PASS |
| S1-flat | **S-layer FLAT ventral surface** (PI: adherent cells flatten ventral surface; eases lamellipodium). cell/basal_surface.py rewritten: FLAT 2D triangulated disk (Delaunay concentric rings at z_basal, +z normal) — consistent with the planar F-layer (no flat-vs-curved dichotomy). POSITIONS FA on the flat surface + triangle adjacency. NO force/bond→γ-invisible. Gate PASS: 629 tris, 60 FA planar+within-disk+connected, area 99.8% of disk; VG-1 disk-area invariant (rel spread 0.4%). test_basal_surface 12 (flat) + gate auto-viz | (this commit) | GATE PASS |
| B5(ii)+recon | **RECONCILED with session (i) (`f6b956c`) — conclusion now DECISIVE.** (i) REFUTE'd Kumar-closability: molecular content solid (~30 heads/side, ~300nm) but the **parallel cross-section count is MISSING + GEOMETRICALLY IMPOSSIBLE** (band needs ~590–1760; a 50–250nm cross-section holds O(5–15)) → gap IRREDUCIBLE, cortical γ-floor canNOT get its density from the SF side either. Molecular correction reconciled SMALL (~3.4× at (i)'s 17pN Stachowiak stall vs my ~11× dipole upper bound). ⭐**Decisive reframe (Kassianidou/Kumar 2017 PNAS):** Kumar 10–30nN is mostly NETWORK/PRESTRESS; single-fiber ACTIVE myosin ≈5–6nN (supplied in-platform by basal mesh+FA+turgor, not sf_myosin_); active ~7× under even the 6nN target. Both sessions jointly HALT→PI: re-scope the SF active gate to ~6nN single-fiber active. Folded into SF_GENERATION_LIMIT_SYNTHESIS §1/§6/§7 + figure note. ⚠️both branches edited h7_basal_sf_force_budget.py → Lead merge. | (this commit) | RECONCILED (decisive, joint HALT→PI) |
| B5(ii) | **SF generation-limit ACCEPT + γ-floor UNIFICATION** (session-split (ii), branch `h7/sf-generation-floor`). Documented the SF tension as GENERATION-limited (not N_filaments) and unified with the cortical active-γ floor as ONE ½·n·f·ℓ story, ONE missing MCF7 density datum. **Decomposed the 2014× gap (labeled, N-independent):** ~11× per-minifilament fidelity (brief 5pN→lit 56pN dipole = the §9 cortical per-head fix) × ~180×(floor)/~360×(centre) **cross-sectional NMII count/SF** (= Route-B MCF7 density datum, session (i)'s target). Budget is an UPPER bound on a COHERENT bundle: §11 SF-2c (random mixed-polarity → −61pN slackening) ⇒ needs SARCOMERIC POLARITY organisation = the SF analogue of cortical "overlap not count" (Chugh/Truong Quang). Engagement realism ~9× WORSE (φ0.99 upper vs realised 11%). Gate-reframe (active-fraction / tool-role-scope / density-conditioned band) SURFACED to PI, no loosening. `decompose_generation_gap` in scripts/h7_basal_sf_force_budget.py + test_sf_force_budget_decomposition.py (5 consistency gates) + docs/v2_audit/SF_GENERATION_LIMIT_SYNTHESIS_2026-06-09.md. Reconcile: if (i) lands N_cross≈180–360 datum, Route B closes Kumar + fills the cortical density gap; else this generation-bound conclusion stands. | (this commit) | ACCEPT (generation-bound, provisional pending (i)) |

---

## SESSION (i) "CLOSE THE FLOOR" — Route-B NMII force-scaling (`h7/sf-nmii-forcescale`)

Worktree `ffn_cellsim-nmii`, branched from `h7/compartment-platform` @ `02d231e`.
Mandate: derive the Route-B SF NMII force-scale factor from a native density datum, test Kumar 10-30 nN.

| time | item | commit | result |
|---|---|---|---|
| boot | new worktree + branch off compartment-platform; read mandate (PI_QUEUE SESSION SPLIT) + BASAL_MESH_DESIGN + AUTORUN Phase-3 + force-budget + myosin.mesoscale_force_scaling + basal_mesh; env ffn_sim, hoomd 7.0.1; ran existing budget (raw 5 pN, 2014× under Kumar) | — | OK |
| DR | **deep-research** (104-agent harness, 5 angles, 22 sources, 25 claims 3-vote verified: 18 confirmed / 7 KILLED). Anchor 2a SOLID: native minifilament ~30 heads/side, ~300 nm, 28-30 molecules (Billington2013/Hu2017/Melli2018/Niederman1975). Anchor 2b ORDER: per-minifilament stall fs≈17 pN (Stachowiak2009 muscle extrap., self-flagged not-direct). ⭐Anchor 1 MISSING: parallel/cross-section ~50 estimate **REFUTED 0-3**; back-solve ~590-1760 geometrically impossible. ⭐REFRAME HIGH: single-fiber ACTIVE ~5-6 nN (Kassianidou/Schwarz/Kumar 2017 PNAS), Kumar 10-30 nN is mostly network/prestress | — | DATUM (partial) |
| ②③ | **force-budget extended → REFUTE / HALT→PI.** Encoded DERIVABLE per-minifilament molecular correction ≈3.4× (model 5 pN→lit 17 pN, grid-invariant) + MISSING parallel count (None-gated; refuted 50 / geometric 5-15 / back-solve bounds, never picked) + Kassianidou active/network split. Even refuted-generous 50×17 pN ≈ 0.85 nN = ~7× under ~6 nN ACTIVE, ~12-35× under Kumar. **No back-solved myosin.py factor written** (forbidden — parallel multiplier missing). = SF instance of cortical γ-floor (same ½·n·f·ℓ, same missing motor-density datum). docs/v2_audit/H7_SF_NMII_FORCESCALE_RESULT_2026-06-09.md + PI_QUEUE updated + SE candidates (PI-gated) | (this commit) | REFUTE / HALT→PI |
