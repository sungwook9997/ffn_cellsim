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

# Compartment Activation — Detailed Plans (2026-06-09)

> **PARAMETER RESEARCH ONLY (PI 2026-07-16).** Counts and literature candidates below may inform I0-Bn, but
> every mechanism must be implemented in the Warp-CUDA-only `aleph/ac/` engine. Any HOOMD/CPU execution,
> fallback, built-in, or wiring instruction in this historical document is superseded.

Per-compartment activation plans for the 8 default-OFF compartments. Researched in parallel (workflow `w0nc67eil`, 8 agents) and Lead-curated. Companion to `COMPARTMENT_ACTIVATION_MASTER_PLAN_2026-06-09.md` (the phased roadmap + the universal 5-step activation procedure these plans instantiate).

> **Every proposed value below is a CANDIDATE for PI ratification** — tagged `SOLID_LITERATURE` / `ORDER_ESTIMATE` / `NEEDS_MEASUREMENT` / `DERIVED` + a confidence. Nothing is wired or auto-adopted; the enabled build still raises until PI ratifies the blocking constants. Order-estimates are labelled as such — no invented numbers.

## Index (by activation phase)

| Phase | Compartment | Blocking parameters | Effort |
|---|---|---|---|
| B | `ventral_stress_fibers` | N_filaments; k_anchor; Bundle bending rigidity EI | ~2-3 focused days. BLOCKER #2… |
| D | `linc` | k_linc | MEDIUM. Mechanism is already fully authored… |
| D | `intermediate_filaments` | E_if; Nonlinear constitutive curve —…; ratio_xl | Medium, ~2-4 focused days IF the PI… |
| D | `microtubules` | Y_stretch; L_mt; n_mt | ~2-3 days. Day 1: denylist fix +… |
| C | `osmotic_regulation` | Lp; delta_c | "Low-to-moderate, ~0.5-1 day. The module, its… |
| C | `membrane_reservoir` | sigma_crit_bleb; f_excess | "Denylist wiring + regression test… |
| E | `cadherin_junction` | r_bind; r0_trans | MEDIUM-HIGH (~3-5 focused sessions). The… |
| E | `junctional_actin` | k_catch0; k_slip0 | Medium, ~3-5 focused days once… |


---

## Ventral stress fibers — `cell/stress_fibers.py`  *(Phase B)*

**Goal.** Activate explicit FA-to-FA contractile actomyosin bundles whose emergent single-fiber axial tension lands in the Kumar 2006 10-30 nN band at a physical few-percent strain, WITHOUT contaminating the cortical-gamma estimator, by (1) ratifying N_filaments so mu_SF=N*EA_single is set, (2) minting a distinct sf_myosin_* bond prefix, and (3) honouring the SF-stiffest-element CFL against the physiological-cytoplasm per-bead drag.

- **Dependencies:** fa compartment (bridge/fa.py) — HARD: SF anchors onto integrin/fa_actin_clutch particles; extend_snapshot_with_stress_fibers RAISES without fa_endpoint_positions+tags. FA must be…; cortex/myosin.py — the NMII D5/D6 builder must be prefix-parameterized to mint sf_myosin_* (BLOCKER #2; cross-module edit, PI sign-off).; cortex/cortical_tension.py — must add sf_/sf_myosin to its denylist so the contamination guard holds.; Physiological baseline ON: cytoplasm 65.9 Pa.s, turgor 133 Pa, membrane_surface, nucleus (the cell SF tension is measured FROM this baseline, per HARD rule).; cell/manifest.py + cell/cell.py + configs/mcf7_baseline.yaml — the standard 5-point plug-in (manifest block, resolve_*, _extend_snapshot_*, CellBuildOptions flag, CompartmentSpec…; EA_single=4.3e-8 N (Gittes1993_JCB, audit-OK) — solid, already wired as default.
- **Ready when:** PI ratifies N_filaments (proposed 20, range 10-30) so mu_SF is set; the sf_myosin_* prefix is minted and added to the cortical-tension denylist (BLOCKER #2 cleared); the SF CFL is gated against the physiological-cytoplasm per-bead drag (not water); and on a small FA-adhered, equilibrated, physiological-baseline cell the activation gate PASSES: single-SF T_SF in [10,30] nN with NMII ON / ~0 with NMII OFF, grid-invariant across n_beads, AND the cortical-gamma OFF/ON delta is provably UNCHANGED by SF presence (contamination guard). The band is fixed before the run and not loosened.
- **Effort:** ~2-3 focused days. BLOCKER #2 (prefix-parameterizing cortex/myosin.py + extending the cortical_tension denylist, with the contamination integration test) is the bulk — ~1 day, and it is a cross-module edit needing PI sign-off. The 5-point manifest/cell wiring (mirroring the fa/erm/lamellipodium pattern) + SF CFL gate-with-physiological-drag is ~0.5 day. The smoke driver + equilibrated FA-adhered build + gate run (depends on a settled FA operating point) is ~0.5-1 day, partly gated on FA equilibration wall-time. N_filaments ratification and the k_anchor policy call are PI decisions (minutes once surfaced), not implementation time. EI/buckling is explicitly DEFERRED (not in this activation).

### PI parameters to ratify

**`N_filaments (vSF cross-section filament count -> mu_SF = N_filaments * EA_single)`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *dimensionless (filaments per bundle…*
- Proposed: N_filaments = 20 (range 10-30; mu_SF = 8.6e-7 N at EA_single=4.3e-8 N). For the x40 mesoscale a value scaled to effective filaments may be argued, but the PHYSICAL bundle EA (=N*EA_single) is the invariant the validation band probes, so propose the native cross-section count 20 and let measure_sf_tension PASS/FAIL on merit.
- Citation: Cramer, Siebert & Mitchison 1997, J Cell Biol 136:1287-1305 (graded-polarity vSF actin bundling, ~10-30 filaments per cross-section); cross-check Deguchi, Ohashi & Sato 2006, J Biomech 39:2603-2610 (isolated single SF tensile test) and Lu, Oswald, Ngu & Yin 2008, Biophys J 95:6060 (SF mechanics).
- Notes: BLOCKER #1: with N=None the enabled build HALTS (NotImplementedError in resolve_backbone_k_bond). N=20 gives mu_SF=8.6e-7 N -> k_bond=1.65 N/m on a 12um/24-bead bundle, reaching the Kumar band at 1.2-3.5% strain (physical). EA_single (4.3e-8 N) is already audit-OK (Gittes1993_JCB verdict=OK). The Cramer 1997 count is a native-architecture number, NOT a mesoscale default — that is exactly why it is PI-pending. Surface both options (native-count-20 vs x40-scaled) to PI; do NOT auto-pick.

**`sf_myosin bond-type prefix (currently reuses cortex_myosin_*)`** — `DERIVED` / conf **HIGH** — *bond-type name (string prefix)*
- Proposed: Mint 'sf_myosin_backbone' / 'sf_myosin_head_backbone' / 'sf_myosin_attach_b{i}' by parameterizing the cortex myosin builder's hardcoded prefix (cortex/myosin.py lines 354-355,358,694-699). Add 'sf_myosin' to GAMMA_DENYLIST_PREFIX coverage so the SF motor load path is denylisted DISTINCTLY from cortical motors.
- Citation: Stam et al. 2017 PNAS (D5 bipolar minifilament); Hill 1938 (D6 F-V); reused verbatim, no new physics. Denylist convention: cortex/cortical_tension.py ADHESION_BOND_TYPE_PREFIXES pattern.
- Notes: BLOCKER #2 (the registry's named ACTIVATION BLOCKER): cortex_myosin_* IS the active-gamma signal (gamma(myo ON)-gamma(myo OFF)). If SF reuses it, a LIVE SF build sums SF motor tension into cortical gamma -> contamination. Interim mitigation already available: denylist 'cortex_myosin_' on SF builds (registry denylist_bond_types already lists both 'sf_' and 'cortex_myosin_'), but that ALSO blinds the cortical estimator to real cortical myosin -> unacceptable for a co-build. The clean fix is the distinct prefix. This…

**`k_anchor (SF end-bead <-> FA clutch anchor bond stiffness)`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *N/m*
- Proposed: Default k_anchor = k_bond = mu_SF/ell0 (~1.65 N/m at N=20) for now (the anchor as stiff as the bundle it terminates). PROPOSE instead anchoring to the FA molecular-clutch stiffness already in the build: k_int_bare = 1.0e-3 N/m (1 pN/nm, KU-2.7) from phase1_h4.yaml. The SF-FA series stiffness is then dominated by the SOFTER clutch (~1e-3 N/m), so physically k_anchor ~ k_int_bare, not k_bond.
- Citation: phase1_h4.yaml integrin.k_int_bare=1.0e-3 N/m (KU-2.7, talin/vinculin clutch); k_int^eff=k_int_bare*(1+alpha_vin*N_vin). Talin/vinculin clutch stiffness ~1 pN/nm: Chan & Odde 2008 Science 322:1687; Elosegui-Artola 2016 Nat Cell Biol.
- Notes: The module currently defaults k_anchor to k_bond (~1650x stiffer than the real clutch). Surface to PI: should the SF terminate on a STIFF anchor (idealized rigid FA, k_anchor=k_bond) or on the REAL compliant clutch (k_anchor~k_int_bare=1e-3 N/m)? The latter is more physiological (the FA clutch is the real load-transfer element) but softens the measured SF tension transmission. Either is defensible; PI must pick. Note: a soft anchor relaxes the CFL (k_stiffest=max(k_bond,k_anchor)=k_bond either way, so CFL is…

**`Bundle bending rigidity EI (angle harmonic / buckling)`** — `NEEDS_MEASUREMENT` / conf **LOW** — *N*m^2 (flexural rigidity)*
- Proposed: DEFER for first activation (module builds a straight FA->FA chord, no angle potential). When added: EI_bundle ~ N_filaments * EI_single with EI_single = 7.3e-26 N*m^2 (F-actin, Lp~17um, kappa=Lp*kT). For a TIGHTLY-crosslinked bundle EI can scale up to N^2*EI_single (Bathe 2008); a loosely-coupled bundle scales ~N. Proposed first cut: EI = N_filaments * EI_single = 1.46e-24 N*m^2 (N=20), the lower (decoupled) bound.
- Citation: Gittes, Mickey, Nettleton & Howard 1993, J Cell Biol 120:923-934 (F-actin Lp~17um -> EI_single~7.3e-26 N*m^2); Bathe, Heussinger, Claessens, Bausch & Frey 2008, Biophys J 94:2955 (crosslinked-bundle EI scaling N..N^2).
- Notes: Documented TODO in PI_DECISIONS; NOT an activation blocker for the tension observable (a contractile bundle is in TENSION, so buckling/EI does not gate the Kumar-band test). Flag it so the PI knows the first activation has zero bending rigidity (the bundle cannot buckle under transient compression). Add only if a compression/buckling observable is later required.

**`n_SF / n_beads_per_SF / sarcomere_spacing (geometry)`** — `SOLID_LITERATURE` / conf **HIGH** — *count / count / m*
- Proposed: n_SF = 20 (sparse basal layout, ~1440 particles); n_beads_per_SF = 24 (grid-invariant per the test across {12,24,48,96}); sarcomere_spacing = 0.5e-6 m (lower band of non-muscle vSF alpha-actinin periodicity).
- Citation: Tojkander, Gateva & Lappalainen 2012, J Cell Sci 125:1855-1864; Hotulainen & Lappalainen 2006, J Cell Biol 173:383-394 (sarcomeric alpha-actinin periodicity 0.5-1 um in vSFs).
- Notes: These already have defaults in resolve_stress_fibers and are NOT the blockers. n_SF is a layout choice (PI may scale); sarcomere_spacing default 0.5um is the tighter banding seen in non-muscle vSFs. Grid-invariance of the tension w.r.t. n_beads is asserted in the test, so n_beads is free.

**`EA_single (F-actin single-filament axial rigidity; FIXED input, not PI-pending)`** — `SOLID_LITERATURE` / conf **HIGH** — *N*
- Proposed: 4.3e-8 N (E~1.3-2.6 GPa x A~3.2e-17 m^2)
- Citation: Gittes, Mickey, Nettleton & Howard 1993, J Cell Biol 120:923-934 (audit verdict=OK in source_audit); Kojima, Ishijima & Yanagida 1994, PNAS 91:12962-12966 (F-actin tensile rigidity).
- Notes: Already hardcoded as the default and audit-clean. Listed for completeness so the PI sees the full mu_SF = N_filaments * EA_single provenance: only N_filaments is pending; EA_single is solid.

### Activation gate

- **Observable:** Single ventral-stress-fiber mean AXIAL TENSION T_SF = mean over one bundle's sf_actin_bond family of k_bond*max(0,|dr|-ell0), reported by measure_sf_tension(...). Measured on a fully-built, FA-adhered, EQUILIBRATED cell at the physiological baseline (cytoplasm 65.9 Pa.s, turgor 133 Pa, membrane_surface + nucleus ON), with SF NMII attached and stepped to steady state. SECONDARY observable: cortical-gamma delta must be UNCHANGED by SF presence (contamination guard).
- **Literature band:** Kumar, Maxwell, Heisterkamp, Polte, Lele, Salanga, Mazur & Ingber 2006, Biophys J 90:3762-3773 (single-living-SF laser ablation, Kelvin-Voigt retraction): single-SF contractile tension ~10-30 nN. PASS band = [10 nN, 30 nN]. (This is a VALIDATION target only; mu_SF is an independent EA-derived input, so the band can genuinely fail.)
- **Controls:**
  - OFF/ON differencing on the MEASURE: SF NMII OFF -> T_SF should be ~0 (passive force-free bundle, anchors at r0=0, backbone at rest); SF NMII ON -> T_SF rises into band. The ON-OFF difference is the active SF tension.
  - CONTAMINATION GUARD (mandatory): build the SAME cell with SF OFF and with SF ON; gamma(myosin_ON)-gamma(myosin_OFF) measured over CORTEX bonds must be IDENTICAL (within MC noise) in both cases. If SF presence shifts cortical gamma, the sf_ / sf_myosin_…
  - RIGID/LIMIT PARITY: grid-invariance check (already in test_stress_fibers) — T_SF independent of n_beads in {12,24,48,96} to <1% (confirms k_bond=mu_SF/ell0 convention). Plus a stiff-anchor vs compliant-anchor (k_anchor=k_bond vs k_anchor=k_int_bare) parity…
  - r0 provenance check: measure_sf_tension MUST be called with the registered construction r0 = mean(layout.ell0_actin); assert r0_is_live==0.0 (a live-mean r0 under-reports a uniformly contracted bundle's tension to ~0).
  - BASAL-PLANE-only check: confirm T_SF is summed strictly over sf_actin_bond (chain) and that integrin_ligand / fa_actin_clutch* are NOT in the SF sum (SF tension != FA traction).
- **PASS / PARTIAL / REFUTE:** PASS: mean T_SF in [10,30] nN with NMII ON, ~0 with NMII OFF, AND cortical-gamma delta unchanged by SF presence (<5% / within MC noise), AND grid-invariant (<1% across n_beads). PARTIAL: T_SF within a factor ~2 of the band (5-60 nN) OR in-band only at the edge of the N_filaments=10-30 range -> report as order-consistent, surface N_filaments sensitivity to PI, do NOT loosen the band. REFUTE: T_SF off by >1 order at the ratified N_filaments (mu_SF…
- **No-γ-contamination check:** Two-layer denylist: (1) every bond this module mints is sf_-prefixed (GAMMA_DENYLIST_PREFIX='sf_'); (2) SF NMII must use a NEW sf_myosin_* prefix (BLOCKER #2) so it is denylisted SEPARATELY from cortex_myosin_*. The registry already lists denylist_bond_types=('sf_','cortex_myosin_'); cortical_tension.cortical_bond_typeid_mask must be extended to add 'sf_' and 'sf_myosin' to…

### Wiring steps

1. BLOCKER #1 (PI ratify): set stress_fibers.N_filaments=20 (or an explicit mu_SF) in the SF config block so resolve_stress_fibers fills mu_SF=N_filaments*EA_single=8.6e-7 N; otherwise resolve_backbone_k_bond raises NotImplementedError and the enabled build HALTS. Also set k_actin=mu_SF/ell0 OR leave k_actin=None and supply it for measure_sf_tension (note: measure_sf_tension currently REQUIRES…
2. BLOCKER #2 (cross-module, PI sign-off): parameterize the cortex myosin builder (cortex/myosin.py: BOND_TYPE_MYOSIN_BACKBONE/HEAD_BACKBONE at lines 354-355, cortex_myosin_attach_bin_names at 358, and the hardcoded type strings at 655-656,694-699) to accept a prefix arg, then build SF NMII with prefix 'sf_myosin'. Update sf_bond_type_names() / GAMMA_DENYLIST coverage to include the sf_myosin_*…
3. Extend cortex/cortical_tension.py: add 'sf_' and 'sf_myosin' to ADHESION_BOND_TYPE_PREFIXES (or a new SF denylist tuple) so cortical_bond_typeid_mask / _gamma_soft exclude the SF load path. (Registry already declares denylist_bond_types=('sf_','cortex_myosin_') — make the estimator honour 'sf_'.)
4. Add a manifest optional_subsystems.ventral_stress_fibers block in configs/mcf7_baseline.yaml (manifest_path is currently None in the registry; set it to ('optional_subsystems','ventral_stress_fibers')). Block carries: enabled:false (default-OFF), N_filaments:20, n_SF:20, n_beads_per_SF:24, sarcomere_spacing:0.5e-6, base_config pointer, and k_anchor policy.
5. Add resolve_stress_fibers wiring into cell/manifest.py: import resolve_stress_fibers, add p_stress_fibers to ResolvedBaseline + the compartments() dict, resolve it in resolve_baseline under the opt block (mirror the fa/erm/lamellipodium pattern at manifest.py:199-261), gated on the fa block being enabled (requires=('fa',)).
6. Add to cell/cell.py: a CellBuildOptions.with_stress_fibers flag, an _extend_snapshot_with_stress_fibers call AFTER _extend_snapshot_with_fa (so FA integrin/clutch tags exist), passing fa_endpoint_positions + fa_endpoint_tags harvested from the FA integration result (the FA clutch anchor particles). Then register_stress_fiber_bond_params on the md.bond.Harmonic, and attach the (prefixed) SF NMII…
7. Add the SF CFL gate enforcement: call stress_fiber_dt_cfl(p, gamma_b=<per-bead drag>, ell0=mean(layout.ell0_actin)) and gate dt against it with a cfl_strict flag (mirror erm.py / enclosed_volume.py). CRITICAL: pass the PHYSIOLOGICAL cytoplasm per-bead drag, NOT the cortex water-viscosity gamma_b (see GPU_and_CFL).
8. Wire the basal observable: at measurement time call measure_sf_tension(..., r0=mean(layout.ell0_actin)) and assert r0_is_live==0; emit T_SF + in_band into the production JSON alongside (separately from) the cortical-gamma readout.
9. Update compartment_registry.ventral_stress_fibers: set status EXPERIMENTAL->LIVE only AFTER the gate passes; set manifest_path; keep gamma_contaminating=True until the sf_myosin_ prefix + estimator denylist land, then it can be re-described.

### Test plan

- Run the existing aleph/tests/test_stress_fibers.py (static + analytical): OFF-identity no-op, n_SF=0 no-op, n_beads<2 ValueError, degenerate-bundle ValueError, exact particle/bond counts, sign-sense (stretched chain bond pulls beads together; stretched anchor pulls end bead toward FA), grid-invariance of T=mu_SF*eps across n_beads in…
- New unit test: resolve with N_filaments=20 -> mu_SF=8.6e-7 N and resolve_backbone_k_bond returns ~1.65 N/m on ell0=0.522um; resolve with N_filaments=None and k_actin=None -> resolve_backbone_k_bond raises NotImplementedError (HALT contract preserved).
- New CFL test: stress_fiber_dt_cfl with gamma_b=cytoplasm drag (6*pi*65.9*R_bead) returns dt_CFL in the microsecond range (>> production dt 6.98e-7 s); with gamma_b=water drag it returns sub-ns (a runaway) -> assert the attach helper RAISES under cfl_strict when fed water drag, documenting that the physiological drag MUST be used.
- Contamination integration test (the gate's secondary observable): build a small FA-adhered cell with SF OFF and with SF ON; assert cortical_bond_typeid_mask excludes every sf_*/sf_myosin_* type and includes every cortex_* type; assert gamma_soft identical (within float tol) between the two builds.
- Measurement-protocol test: measure_sf_tension called WITHOUT r0 sets r0_is_live=1.0 and under-reports a uniformly contracted bundle; called WITH the registered r0 reports the correct mu_SF*eps. Assert measure_sf_tension raises NotImplementedError when k_actin is None (current contract).
- Smoke driver (mirror h7_gate_b_native): build FA-adhered + SF cell, equilibrate, step SF NMII to steady state, dump T_SF and the cortical-gamma OFF/ON delta to outputs/h7/production/h7_sf_activation_smoke.json; auto-generate the figure (per the production-driver auto-viz rule).

### GPU / CFL

CFL is the sharp edge here. stress_fiber_dt_cfl returns cfl_safety*gamma_b/k_stiffest where k_stiffest=max(k_bond,k_anchor). With N=20, k_bond=1.65 N/m. The dt is acutely sensitive to WHICH per-bead drag the caller passes: (a) the cortex's gamma_b uses WATER viscosity (6*pi*1e-3*R_bead ~ 9.4e-10 N.s/m at R=50nm) -> dt_CFL ~ 0.06 ns, which is ~6000x BELOW the production dt of 6.98e-7 s (Gate-B 30x lever) and would make the build unrunnable. (b) the PHYSIOLOGICAL cytoplasm drag (HARD rule: 65.9 Pa.s) -> gamma_b ~ 6.2e-5 N.s/m -> dt_CFL ~ 3.8 us, comfortably ABOVE 698 ns. RESOLUTION: the SF backbone lives in the cytoplasm, so the attach helper…


---

## LINC complex — `cell/linc.py`  *(Phase D)*

**Goal.** Activate the explicit nesprin-SUN LINC compartment (linc_nesprin harmonic bonds nucleus_bead<->cortex/SF/MT) so the nucleus is mechanically coupled to the cytoskeleton at the physiological operating point, with a PI-ratified single-nesprin effective stiffness k_linc and the ~8 pN resting-tension oracle, WITHOUT contaminating cortical gamma.

- **Dependencies:** nucleus compartment ON (registry requires=('nucleus',)): LINC bonds attach to nucleus_bead tags; with no nucleus the build finds zero acceptors and no-ops. This conflicts with the…; Full physiological baseline ON (cytoplasm 65.9 Pa.s, turgor 133 Pa, membrane_surface, nucleus at KU-3.B2.1) -- LINC is measured FROM this baseline, never from a floppy/null shell.; cortex/cortical_tension.py denylist must learn the 'linc_' prefix BEFORE any gamma-measuring run with LINC ON (currently only knows 'fa_actin_clutch'). Hard dependency for the…; An acceptor cytoskeletal population near the nuclear envelope: perinuclear actin cap (preferred) OR cortex/SF/MT beads within a geometry-derived capture_radius. Without it the…; PI sign-off on k_linc (no-magic-number gate) AND on the 2->8 pN f_rest re-anchor + Dejardin-vs-Arsenovic re-attribution (flagged ORACLE-TARGET in the module) AND KB SourceEvidence…; Native constrained-BAOAB plugin already evaluates md.bond.Harmonic (validated) -- no plugin change needed for route A/B; route C (md.bond.Table) needs a check that the…
- **Ready when:** PI has (1) ratified a literature-anchored k_linc (route A/B/C) and the 8 pN f_rest re-anchor/re-attribution; (2) the six SourceEvidence rows (Rief99, Dejardin20, Arsenovic16, Autore13, Crisp06, Lombardi11) are registered in Notion and verify_sources verdict=OK; (3) the manifest/cell.py/CellBuildOptions/registry wiring + the 'linc_' cortical-gamma denylist are in place with their tests; (4) a full-baseline build forms n_bridges>0 at a geometry-derived capture_radius (no zero-bonds trap) and the OFF/ON active-gamma differencing shows LINC does not move cortical gamma; at which point the activation-gate run can measure <T_linc> against the [2,10] pN oracle.
- **Effort:** MEDIUM. Mechanism is already fully authored (resolver + pair_linc_bridges + extend_snapshot_with_linc + configure_linc_bond_potential + force/potential/CFL helpers + sanity-gate tests all exist and are correct). Remaining work is integration + ratification, not new physics: ~0.5 day SE/KB registration + verify; ~0.5 day PI ratification of k_linc (route A is a same-session decision; route C md.bond.Table is +1-2 days of nonlinear-potential build if chosen); ~1 day manifest/cell.py/registry wiring + the 'linc_' denylist + tests; ~1 day acceptor-topology decision + perinuclear-cap seeding (the real geometric risk) + the gate run. Total ~3 days for route A (effective harmonic), ~4-5 days if route C (full-fidelity tabulated nesprin potential) is taken. Main risks: the zero-bonds-at-real-geometry acceptor trap (FA precedent) and the missing cortical-gamma denylist wiring -- both must be closed before the gate can pass.

### PI parameters to ratify

**`k_linc (single-nesprin effective linear stiffness)`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *N/m*
- Proposed: 1.0e-2 N/m (folded-rod / pre-unfolding secant; PRIMARY candidate). Bracket: low 2.5e-4 (unfolded-WLC entropic) -- high 1.0e-2 (folded triple-helix pre-unfold slope). Operating-band cross-check 0.8-1.6e-3 from f_rest/extension at delta=5-10nm.
- Citation: Rief, Pascual, Saraste & Gaub (1999) J Mol Biol 286:553 -- single spectrin repeat unfolds at 25-35 pN, unfolded-strand WLC persistence length Lp~=0.5 nm, shallow unfolding potential width ~1.5 nm. Folded-rod pre-unfolding loading slope ~30 pN over ~3 nm => ~1e-2 N/m. Entropic WLC k=3kT/(2 Lp L0):…
- Notes: GENUINELY NO single Hookean N/m exists for nesprin-2G (nonlinear sequential repeat unfolding). This is the load-bearing PI decision. THREE routes, pick one: (A) effective-harmonic at 1.0e-2 N/m = folded-rod pre-unfolding secant (PI-ratifiable now, CFL-safe, simplest, default-OFF until set); (B) WLC entropic k=3kT/(2 Lp L0) with caller-supplied operating contour L0~200nm tail span => ~6e-6 N/m (very soft, lets nucleus drift); (C) md.bond.Table tabulated nonlinear potential reproducing the Rief…

**`f_rest (resting per-nesprin tension; validation ORACLE, not a stiffness)`** — `SOLID_LITERATURE` / conf **HIGH** — *N*
- Proposed: 8.0e-12 N (8 pN)
- Citation: Dejardin et al. (2020) J Cell Biol 219(10):e201908036 (DOI 10.1083/jcb.201908036) -- mini-nesprin-2G CB TSMod FRET biosensor: nesprin-2G CB held under ~8 pN resting tension generated by actomyosin + microtubules, balanced by cell-cell adhesion. ATTRIBUTION FIX: NOT Arsenovic 2016 -- Arsenovic et…
- Notes: Already hard-coded as _NESPRIN_REST_TENSION_N=8.0e-12 with the correct Dejardin-vs-Arsenovic attribution in the docstring. This is the GATE OBSERVABLE oracle (mean equilibrium linc_nesprin bond tension should sit near 8 pN once k_linc and pre-strain are set). The 2 pN->8 pN re-anchor + re-attribution flagged in the module as ORACLE-TARGET still needs PI sign-off + a fresh KB SourceEvidence row before the gate is treated as live (see contamination_check / dependencies). Caveat: Dejardin's 8 pN is a…

**`r0 (linc_nesprin bond rest length)`** — `SOLID_LITERATURE` / conf **HIGH** — *m*
- Proposed: 50e-9 m (perinuclear-gap geometric default) -- OR a caller-supplied (R_cell - R_nuc)-scale rest length when bonding cortex-shell acceptors directly
- Citation: Crisp et al. (2006) J Cell Biol 172:41 -- SUN-KASH bridge spans the ~30-50 nm perinuclear luminal space. GEOMETRIC length (derivable, not a tuned force knob).
- Notes: Already defaulted (_PERINUCLEAR_GAP_M). Physical caveat: 50 nm is the real nesprin span, but at the x40 mesoscale the nearest nucleus_bead<->cortex_bead separation is R_cell - R_nuc ~= 5.6 um (>100x r0). If LINC bonds the cortex shell directly, r0 should be the perinuclear-acceptor span (keep 50 nm) and acceptors must be SEEDED near the envelope; if bonding far cortex beads, r0 stretches and f_rest is meaningless. See wiring_steps re: perinuclear acceptor seeding -- this is the same zero-bonds-at-real-geometry…

**`capture_radius (geometric nucleus<->cytoskeleton pairing radius)`** — `DERIVED` / conf **HIGH** — *m*
- Proposed: Production value MUST be passed from cell geometry, NOT the 150 nm (3x perinuclear-gap) default. Either (a) seed perinuclear/cap acceptors near the envelope and keep ~150-300 nm, or (b) bond cortex shell with capture_radius on the (R_cell - R_nuc) ~ 5.6 um scale.
- Citation: Geometric, grid-invariant (multiple of a length). resolve_linc already warns when capture_radius < (R_cell - R_nuc).
- Notes: Not a magic number; a build-time geometric pairing scale (mirrors fa_clutch_capture_radius). The 150 nm default forms ZERO bonds at real cell geometry (module docstring documents this). Production capture_radius is set by the chosen acceptor topology (perinuclear cap vs cortex shell), per the physiological-baseline rule.

**`n_bridges (LINC bond count at mesoscale)`** — `DERIVED` / conf **HIGH** — *count*
- Proposed: = n_nuc envelope beads (one LINC per envelope bead, k=1 nearest-acceptor); n_bridges_max=1e6 sentinel never binds
- Citation: Lombardi et al. (2011) J Biol Chem 286:26743 (thousands of LINC per nucleus); Crisp 2006. Mesoscale count set by nucleus bead resolution, not the sentinel.
- Notes: Topology-only, no new particle types (see wiring_steps note re: registry particle_types_added mismatch).

### Activation gate

- **Observable:** Mean equilibrium per-bond tension of the linc_nesprin bond family <T_linc> = <k_linc * (|dr| - r0)> measured over all formed LINC bonds on the FULL physiological-baseline cell (cytoplasm 65.9 Pa.s, turgor 133 Pa, membrane_surface ON, nucleus ON at KU-3.B2.1 setpoint), after equilibration. SECONDARY observable: nuclear-cytoskeletal mechanical coupling ratio (Delta of nucleus centroid displacement / cortex displacement under a myosin-ON vs myosin-OFF differencing), a Lombardi-2011-style decoupling check.
- **Literature band:** PRIMARY: <T_linc> in [2, 10] pN, centered on the Dejardin 2020 ~8 pN resting-tension oracle (band widened from the point value because Dejardin is a cell-cell-adhesion-balanced fibroblast/epithelial datum, not a suspended/adherent MCF7 single cell; Arsenovic 2016 TSMod calibrated range ~1-6 pN bounds the sensor sensitivity). SECONDARY: nucleus-cortex coupling ratio > 0 and DECREASES toward 0 when…
- **Controls:**
  - OFF/ON differencing: build the identical cell with linc.enabled=False vs True; the active-gamma number gamma(myosin ON) - gamma(myosin OFF) MUST be UNCHANGED by LINC within noise (LINC is denylisted, so it cannot move cortical gamma -- this is the…
  - Rigid/limit parity: in the k_linc -> infinity (very stiff) limit the nucleus centroid must rigidly co-move with its bonded cortex anchors (coupling ratio -> 1); in the k_linc -> 0 limit the nucleus must decouple (coupling ratio -> 0, recovering the no-LINC…
  - No-contamination check: assert linc_ is on cortex.cortical_tension denylist (GAMMA_DENYLIST_PREFIX='linc_'); assert measure_cortical_tension excludes every linc_nesprin bond from the method-of-planes sum (no linc_ type appears in the cortical bond mask).
  - Geometry/no-op control: confirm n_bridges > 0 at the production capture_radius (guard against the zero-bonds-at-real-geometry trap); a build that silently forms 0 LINC bonds is a REFUTE of the activation, not a pass.
  - Newton-3rd-law / momentum control: total LINC force on the system = 0 (builtin md.bond.Harmonic is equal-and-opposite); assert no net momentum injected (contrast with the lab-frame ERM field).
- **PASS / PARTIAL / REFUTE:** PASS: <T_linc> in [2,10] pN, n_bridges>0, active-gamma unchanged by LINC (|Delta gamma| < gamma noise), rigid/limit parity monotonic, momentum conserved, AND k_linc is a PI-ratified literature-anchored value (not None, not back-solved from the 8 pN target). PARTIAL: <T_linc> order-correct (within [0.5,30] pN) but outside [2,10] -- physics is alive and uncontaminating but the chosen k_linc/pre-strain/acceptor-topology needs revisiting (surface to…
- **No-γ-contamination check:** linc_nesprin bonds carry a nucleus<->cytoskeleton load path that is NOT cortical hoop tension. The cortical-gamma estimator MUST denylist them by the 'linc_' prefix. This is ALREADY declared: module GAMMA_DENYLIST_PREFIX='linc_' and registry CompartmentSpec(name='linc', gamma_contaminating=True, denylist_bond_types=('linc_',)). REQUIRED WIRING (currently MISSING and the gate cannot pass without…

### Wiring steps

1. 1. SE/KB FIRST (hard prerequisite -- no-magic-number + citation-integrity): register SourceEvidence rows in Notion (SoT) for Rief 1999 (J Mol Biol 286:553), Dejardin 2020 (JCB 219:e201908036), Arsenovic 2016 (Biophys J 110:34), Autore 2013 (PLoS ONE 8:e63633), Crisp 2006 (JCB 172:41), Lombardi 2011 (JBC 286:26743). NONE currently exist in the KB (tag_query: 0 rows for…
2. 2. PI ratifies k_linc (the load-bearing decision): pick route A (effective-harmonic 1.0e-2 N/m) / B (WLC entropic w/ operating contour) / C (md.bond.Table nonlinear). Recommend A as v1 activation + C filed as the full-fidelity follow-on. Set k_linc in configs/mcf7_baseline.yaml under optional_subsystems.linc (NOT a code default -- keep module default None so an un-ratified build still raises).
3. 3. manifest.py optional_subsystems wiring: add `linc_b = opt.get('linc'); if _enabled(linc_b): p_linc = resolve_linc(_opt_cfg(linc_b), R_cell=p_cortex.R_cell, R_nuc=<nucleus R_nuc>)` in resolve_baseline (mirror the fa/erm/lamellipodium blocks); add p_linc to ResolvedBaseline dataclass + its dict (mirror p_nucleus). Guard: raise if linc enabled but nucleus disabled (LINC requires nucleus --…
4. 4. cell.py build: add a CellBuildOptions flag (e.g. enable_linc), an _extend_snapshot_with_linc call wiring nucleus_bead tags as nucleus_tags and the chosen acceptor class tags (cortex actin / SF / MT) as cytoskeleton_tags, and a configure_linc_bond_potential call on the single shared md.bond.Harmonic (writes k=k_linc, r0). The extender + potential helpers already EXIST in linc.py -- this is pure…
5. 5. ACCEPTOR TOPOLOGY (the make-or-break geometric step): decide acceptor class + capture_radius to avoid the zero-bonds trap. Option (a) PREFERRED: seed a thin perinuclear actin-cap / acceptor bead shell near the envelope (r ~ R_nuc + ~150 nm) and keep r0=50 nm, capture_radius ~150-300 nm -- physiologically faithful (perinuclear actin cap is real). Option (b): bond the existing cortex shell…
6. 6. CORTICAL-GAMMA DENYLIST (mandatory, currently missing): add 'linc_' to cortex/cortical_tension.py ADHESION_BOND_TYPE_PREFIXES (or, cleaner, make the estimator read registry.gamma_denylist() as the single source). Add a test asserting registry.gamma_denylist() contains 'linc_' AND measure_cortical_tension excludes a built linc_nesprin bond. Without this, LINC on a gamma run is a contamination…
7. 7. Registry hygiene: fix the CompartmentSpec(name='linc') performance contract -- it claims particle_types_added=('linc_anchor',) and n_bonds with a 'linc_anchor' particle, but the module adds ZERO particle types (topology-only, reuses nucleus_bead + cytoskeleton beads). Correct to particle_types_added=() and flip status EXPERIMENTAL->LIVE only after the gate passes + PI sign-off. Add…
8. 8. configs: add optional_subsystems.linc block to mcf7_baseline.yaml (enabled default false; k_linc, capture_radius from geometry, acceptor_class). Production config turns it ON at the physiological value when running the confined_migration / whole-cell recipe (registry expected_on_recipes=('confined_migration',)).

### Test plan

- Static/analytic (extend tests/test_linc.py): assert linc_bond_force_along_axis sign (stretched>r0 -> +, pulls nucleus toward cortex anchor; compressed -> -); assert linc_bond_potential >=0 and =0 at sep=r0; assert linc_cfl_dt_max = 0.1*gamma/k_linc.
- Dimensional/unit: k_linc [N/m] * (sep-r0)[m] = F[N]; 0.5*k*d^2 = [J]. f_rest [N] is NEVER divided by a length to fabricate k_linc (assert no such code path exists).
- Resolver boundary: enabled=False -> all-zero no-op; k_linc=None & enabled -> resolve OK but extend_snapshot_with_linc/configure_linc_bond_potential RAISE NotImplementedError; k_linc<=0 / r0<=0 / capture_radius<=0 -> ValueError; capture_radius>R_cell -> ValueError; capture_radius<(R_cell-R_nuc) -> warns (zero-bonds trap).
- Pairing (pair_linc_bridges): no acceptor in range -> empty (0,2); k=1 nearest-acceptor; shortest-bridge-first cap; empty nucleus or cytoskeleton -> empty.
- Build-integration on full physiological-baseline cell: assert n_bridges>0 at production capture_radius (REFUTE guard); assert one linc_nesprin bond per envelope bead within range; assert linc_nesprin is in snap.bonds.types and existing bond types/groups are carried across unchanged.
- Conservation: sum of LINC forces = 0 (Newton 3rd law, builtin Harmonic); no net momentum injected (contrast ERM field).
- GAMMA CONTAMINATION (critical): build cell with LINC ON, assert measure_cortical_tension's cortical bond mask EXCLUDES every linc_nesprin bond; assert active-gamma = gamma(myo ON)-gamma(myo OFF) is identical (within noise) with LINC ON vs OFF; cross-check registry.gamma_denylist() contains 'linc_'.
- CFL: assert the production dt (~1e-7 s, cortex/nucleus-limited) satisfies dt <= linc_cfl_dt_max(k_linc, gamma_bead) -- verified head-room is ~4 orders even at k_linc=1e-2 N/m (dt_max~1e-3 s), so LINC never becomes the integration bottleneck.
- Rigid/limit parity sweep: k_linc -> large => nucleus rigidly co-moves with anchors (coupling ratio->1); k_linc -> small => decouples (->0); monotonic.
- Oracle gate run: equilibrate, measure <T_linc>, assert in [2,10] pN band; assert in-band WITHOUT having tuned k_linc to the 8 pN target (the chosen k_linc must be the PI-ratified literature value, recorded in the run JSON provenance).

### GPU / CFL

GPU: fully GPU-resident NOW. LINC adds zero particle types and only static linc_nesprin bonds evaluated by the builtin md.bond.Harmonic, which both the stock integrator and the native constrained-BAOAB CUDA plugin already handle on-device. NO md.force.Custom, NO per-step Python, NO cpu_local_snapshot per step. The one CPU cost is a single build-time scipy.cKDTree nucleus<->cytoskeleton pairing (O(n_nuc log n_cyto)), not in the hot loop. native_forcecompute_candidate=NO (already builtin). CFL: dt <= 0.1*gamma_bead/k_linc. With gamma_bead = 6*pi*65.9*a (a~100-250 nm) ~ 1.2-3.1e-4 N.s/m and the PRIMARY k_linc=1e-2 N/m, dt_max ~ 1.2-3.1e-3 s --…


---

## Intermediate filaments / keratin — `cell/intermediate_filaments.py`  *(Phase D)*

**Goal.** Activate the explicit keratin/vimentin IF perinuclear cage with a faithful NONLINEAR strain-stiffening md.bond.Table constitutive law (replacing the linear-only backbone), wired through the manifest -> resolve_baseline -> Cell.build assembler with the if_ bond prefix denylisted from cortical-gamma, validated against single-IF force-extension data (Kreplak 2005 / Block 2018) and the Guo-2013 "VIFs barely touch cortical gamma, dominate intracellular mechanics" signature.

- **Dependencies:** nucleus ON (R_nuc seeds the perinuclear cage shell [R_nuc, R_nuc+thickness]; resolve passes R_nuc). The cage degrades to an inner-band default if nucleus is off, but the…; Physiological baseline ON: cytoplasm 65.9 Pa.s (apply_cytoplasm_drag rescales gamma_if), turgor 133 Pa, membrane_surface — the IF delta-gamma is measured FROM this baseline…; cortex/crosslinkers.py xlink_attach_bin_rest_lengths (already imported) for force-free per-r0 crosslink binning.; cortical_tension.py denylist must be extended to the if_ prefix (wiring step 7) BEFORE any IF-ON gamma run — this is a hard pre-req for the contamination control.; PI ratification of the nonlinear Table curve (the one LOW-confidence pi_parameter) — the nonlinear gate cannot run until build_intermediate_filament_bond_potential(nonlinear=True)…; OPTIONAL future: LINC complex (linc.py) for nucleus<->IF coupling (desmosome/plectin-LINC load path) — not required for this activation, the cage is internally backbone+crosslink…
- **Ready when:** PI has ratified (a) the keratin-vs-vimentin E_if/Lp cell-type pick and (b) the nonlinear Table force-extension curve (r,U,F grid from Kreplak 2005 / Block 2018 / Qin-Buehler), so build_intermediate_filament_bond_potential(nonlinear=True) no longer raises; the if_ prefix is added to the cortical_tension denylist with a passing contamination regression (gamma IF-ON == IF-OFF to machine precision); the manifest/resolve_baseline/Cell.build/gamma_map wiring (steps 3-6) is in place and a short baseline run is numerically stable; and the activation gate PRIMARY (|delta-gamma|<=5% gamma_active) + SECONDARY (stiffening>=3x, extensibility 2-3.5x) both PASS with all five controls green — at which point the CompartmentSpec flips EXPERIMENTAL->LIVE.
- **Effort:** Medium, ~2-4 focused days IF the PI supplies/ratifies the nonlinear curve quickly. Breakdown: ~0.5 day manifest/resolve/Cell.build/gamma_map wiring (well-templated by nucleus/linc, low risk); ~0.5 day cortical_tension denylist extension + contamination regression; ~1-1.5 days the nonlinear md.bond.Table implementation + toe-parity/stiffening unit tests (the real work — digitizing/validating the Kreplak/Block force-extension curve and resolving the Table-vs-Harmonic coexistence on the shared bonded force); ~0.5-1 day the OFF/ON delta-gamma activation-gate runs at the equilibrated FA-adhered baseline + viz. The dominant schedule risk is PI turnaround on the LOW-confidence nonlinear constitutive curve, not engineering; if PI only approves the linear-toe build and defers the nonlinear law, the wiring+denylist+linear-gate subset is ~1 day but the SECONDARY strain-stiffening gate cannot be claimed.

### PI parameters to ratify

**`E_if (small-strain axial Young's modulus, EPITHELIAL KERATIN K8/K18 — MCF7)`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *Pa*
- Proposed: 6.0e6 (keep current default; band [1e6,1e7])
- Citation: Block J. et al. 2018 Sci.Adv. 4:eaat1161 (single-vimentin tensile, low-strain modulus few-MPa regime); Kreplak L. et al. 2005 J.Mol.Biol. 354:569 (single-IF AFM bending/stretch); Guo M. et al. 2013 Biophys.J. 105:1562 (cell-scale VIF regime). Note: the single-IF small-strain *tangent* modulus from…
- Notes: MCF7 is EPITHELIAL -> the physiological network is KERATIN (K8/K18), not vimentin (vimentin is an EMT marker, KB-4.12). Single-keratin small-strain modulus is reported lower/softer than vimentin (Lichtenstern 2012 assembly flexibility). The 6 MPa default sits mid-band and is defensible as an order estimate; PI should ratify a keratin-specific value (likely 2-6 MPa) vs a vimentin variant. Grid-invariant: k_bb = E_if*A_if/l_seg, the sign/dimensional tests do not use E_if.

**`Nonlinear constitutive curve — TABULATED force-extension U(r), F(r) (the PI-pending core…`** — `NEEDS_MEASUREMENT` / conf **LOW** — *U in J vs r in m (md.bond.Table over…*
- Proposed: 3-regime piecewise tangent-modulus tabulation: (I) linear toe E1=E_if for strain eps<eps_unfold~0.1; (II) softening 'unfolding plateau' alpha-helix->beta with reduced tangent modulus E2~0.1-0.3*E_if over eps in [0.1, 0.8] (Kreplak 2005 / Block 2018 plateau); (III) strain-stiffening backbone E3~3-10*E_if for eps>~1.0 up to rupture at eps_max=max_stretch_ratio-1 (i.e. 2-3.5x rest length). Tabulate U by integrating…
- Citation: Kreplak L. et al. 2005 J.Mol.Biol. 354:569-577 (single IF extends 2-3.5x before rupture, force-extension shape); Block J. et al. 2018 Sci.Adv. 4:eaat1161 (nonequilibrium alpha-helix unfolding -> the plateau + tensile memory); Lichtenstern T. et al. 2012 J.Struct.Biol. 177:54 (K8/K18). Qin…
- Notes: THIS IS THE BLOCKER the PI must ratify. md.bond.Table needs (r_grid, U_grid, F_grid). The 3-regime shape (toe -> unfolding plateau -> stiffening) is solid qualitative literature; the exact tangent moduli E2/E3 and plateau force are NOT pinned for the mesoscale bundle and need either digitized Kreplak/Block force-extension curves or a Qin-Buehler MD curve, re-scaled to the x40 bundle cross-section. Until ratified, build_intermediate_filament_bond_potential(nonlinear=True) MUST keep raising NotImplementedError…

**`max_stretch_ratio (finite-extensibility / rupture stretch)`** — `SOLID_LITERATURE` / conf **HIGH** — *dimensionless (r_rupture / l_seg)*
- Proposed: 3.0 (keep default; literature range 2.0-3.5)
- Citation: Kreplak L. et al. 2005 J.Mol.Biol. 354:569-577 (single IF extends ~2-3.5x rest length before rupture); Block 2018 Sci.Adv.
- Notes: This is the ONE solidly-anchored nonlinear number and is already the field default in the module (max_stretch_ratio=3.0). It sets the upper r-bound of the Table grid. Currently provenance-only (linear bond ignores it) -> becomes load-bearing once the Table law is wired as the table's r_max / rupture cutoff.

**`ratio_xl (IF crosslink stiffness / backbone — plectin/filaggrin cross-bridge compliance)`** — `ORDER_ESTIMATE` / conf **LOW** — *dimensionless*
- Proposed: 0.1 (keep default; soft cross-bridge, same convention as actin xlink<<backbone)
- Citation: No single clean single-molecule number at the mesoscale; convention mirrors cortex/crosslinkers.py (xlink k << backbone k). Plectin mechanics: Bonakdar N. et al. 2015 (plectin-deficient cells soften); Na S. et al. 2009 (plectin-vimentin coupling). Wiche G. reviews for plectin as the IF cross-bridge.
- Notes: Genuinely under-determined at this scale. 0.1 is a defensible order estimate (cross-bridges softer than the backbone). PI to anchor or accept as a sensitivity-swept order estimate. Does NOT gate the activation gate observable (cage stiffness is dominated by backbone at large strain).

**`d_if (assembled IF diameter -> cross-section A_if) and Lp (persistence length)`** — `SOLID_LITERATURE` / conf **HIGH** — *m / m*
- Proposed: d_if=10e-9 (band 8-12 nm); Lp=0.5e-6 keratin (band 0.3-0.5 um), 0.4-1.0 um vimentin
- Citation: Mucke N. et al. 2004 J.Mol.Biol. 335:1241 (vimentin Lp 0.4-1 um by AFM, ~10 nm TEM diameter); Lichtenstern T. et al. 2012 J.Struct.Biol. 177:54 (K8/K18 lower Lp ~0.3-0.5 um); Herrmann/Aebi reviews (canonical 10 nm assembled IF).
- Notes: These two are already well-anchored in the module's Magic-Number Block with band guards. No PI action needed beyond the keratin-vs-vimentin Lp pick (couple it to the E_if cell-type pick: keratin -> Lp~0.4 um, E_if~3 MPa).

### Activation gate

- **Observable:** Cortical-gamma NON-contamination + intracellular load-bearing signature. PRIMARY: delta-gamma_cortical from adding the IF cage at the physiological MCF7 baseline (cytoplasm 65.9 Pa.s, turgor 133 Pa, membrane_surface + nucleus ON), measured over EXPLICIT cortex bonds only with the if_ prefix denylisted. SECONDARY (the load-bearing proof): single-IF / cage force-extension stiffening ratio E_tangent(large strain)/E_tangent(small strain) measured by stretching a built if_backbone chain across the Table law.
- **Literature band:** PRIMARY: |delta-gamma_cortical(IF ON vs OFF)| <= 0.05 * gamma_active (i.e. IF must contribute <5% of cortical tension) — Guo 2013 Biophys.J. 105:1562 ('VIFs contribute LITTLE to cortical stiffness'). SECONDARY: tangent-modulus stiffening ratio >= 3x between the toe regime and the >100%-strain regime, with finite extensibility 2.0-3.5x rest length before the Table rupture cutoff (Kreplak 2005;…
- **Controls:**
  - OFF/ON differencing: gamma_cortical(IF OFF) vs gamma_cortical(IF ON) at the IDENTICAL equilibrated FA-adhered baseline operating point — the IF number is the difference, not the absolute.
  - No-contamination check: assert every if_ bond type (if_backbone, if_crosslink_b{0..15}) is excluded by cortical_tension's denylist; cross-check by re-running gamma with the if_ prefix forcibly ALLOWED and confirming a large spurious gamma jump (proves the…
  - Linear/nonlinear parity: at small strain (eps<eps_unfold~0.1) the Table-law tangent force MUST match the existing linear harmonic k_bb*(r-l_seg) to <1% — proves the nonlinear law degrades to the validated linear backbone in the toe regime (no hidden…
  - Rigid/limit parity: in the max_stretch_ratio->inf, E2=E1=E3 limit the Table law must reduce bit-close to the current linear-only build (off-identity of the nonlinearity).
  - Force-free construction control: assert mean |r-r0| ~ 0 over if_crosslink bonds at t=0 (born-separation binning), construction energy thermal-scale.
- **PASS / PARTIAL / REFUTE:** PASS: |delta-gamma| <= 5% gamma_active AND stiffening ratio >= 3x AND finite-extensibility in [2,3.5]x AND all five controls green. PARTIAL: |delta-gamma| in (5%,15%] (IF leaking modest cortical load — investigate cage radius / crosslink reach pushing IF bonds near the cortex) OR stiffening ratio in [2x,3x) (plateau/stiffening moduli need PI re-anchor) — report, do NOT loosen, surface to PI. REFUTE: |delta-gamma| > 15% gamma_active (IF cage is…
- **No-γ-contamination check:** Triple-guarded: (1) GAMMA_DENYLIST_PREFIX='if_' asserted at module import on every IF bond-type name (already in code, lines 521-523); (2) cortical_tension._is_adhesion_bond_type must be EXTENDED to also denylist the 'if_' prefix family (it currently only denylists fa_/integrin_ adhesion prefixes — the if_ family is NOT yet in ADHESION_BOND_TYPE_PREFIXES, this is a REQUIRED wiring edit, see…

### Wiring steps

1. 1. PI ratifies the nonlinear Table constitutive curve (the LOW-confidence pi_parameter): supply/approve the (r_grid, U_grid, F_grid) for build_intermediate_filament_bond_potential(nonlinear=True) — 3-regime toe/unfolding-plateau/stiffening from digitized Kreplak 2005 / Block 2018 (or a Qin-Buehler MD curve), keratin-vs-vimentin E_if pick. Until then the nonlinear path keeps raising…
2. 2. intermediate_filaments.py: implement build_intermediate_filament_bond_potential(nonlinear=True) returning an md.bond.Table-ready (U,F,r_min,r_max) from the ratified curve; add register_if_bond_params(...,nonlinear=True) path that writes the Table potential onto a SHARED md.bond.Table (NOT the shared Harmonic — Table and Harmonic are distinct ForceComputes; the cell would then carry TWO bonded…
3. 3. configs/mcf7_baseline.yaml: add an optional_subsystems.intermediate_filaments block (enabled: false default-off) with n_filaments, beads_per_fil, E_if, Lp, d_if, ratio_xl, max_stretch_ratio, cage_thickness, and a base_config: phase1_if.yaml pointer; mirror the nucleus/linc block shape.
4. 4. cell/manifest.py: in resolve_baseline, after the nucleus/membrane_surface optional block, add `if_b = opt.get('intermediate_filaments'); p_if = resolve_intermediate_filaments(_opt_cfg(if_b), kT=kT, R_cell=R_cell, R_nuc=R_nuc_frac*R_cell) if _enabled(if_b) else None`; add `p_if` to ResolvedBaseline dataclass + its to-dict; thread R_nuc so the cage seeds just outside the nucleus.
5. 5. cell/cell.py: add p_if param to Cell.build and build_baseline_cell signatures (mirror p_nucleus); set `enable_if = p_if is not None`; in the snapshot-extension phase call `snap = extend_snapshot_with_if_cage(snap, p_if, centroid=nuc_centroid, gamma_if=gamma_if, seed=...)`; after the shared bond force is built, call `register_if_bond_params(shared_bond, p_if, nonlinear=<ratified>)` (and the…
6. 6. cell/cell.py gamma_map: add `if enable_if: gamma_map['if_bead'] = gamma_if` BEFORE apply_cytoplasm_drag (BAOAB HALTS if any present particle type lacks a gamma entry). Derive gamma_if = 6*pi*eta_water*R_if_bead Stokes drag, then let apply_cytoplasm_drag rescale it to the 65.9 Pa.s cytoplasm (same path as nucleus_bead).
7. 7. cortex/cortical_tension.py: REQUIRED contamination fix — extend the denylist so if_ bonds are excluded. Either add 'if_' to a new INTRACELLULAR_BOND_TYPE_PREFIXES family consulted by _is_adhesion_bond_type (rename to _is_noncortical_bond_type), or have measure_cortical_tension also reject any name.startswith('if_'). Add a regression test that a built IF-ON cell's gamma equals the IF-OFF cell's…
8. 8. compartment_registry.py: on PASS, PI flips the intermediate_filaments CompartmentSpec status EXPERIMENTAL -> LIVE, sets manifest_path=('optional_subsystems','intermediate_filaments'), and clears the nonlinear PI-decision blocker. (Spec already exists at line 783 with denylist_bond_types=('if_',).)

### Test plan

- Unit (extend aleph/tests/test_intermediate_filaments.py): keep all 18 existing tests green (off-identity, derived SI, force-free crosslink, denylist prefix, nonlinear-not-faked-when-unratified). ADD: Table-law toe-regime parity vs linear k_bb to <1% (control 3); Table-law stiffening ratio >=3x between eps=0.05 and eps=1.5; Table…
- Integration: build_baseline_cell with intermediate_filaments ON at the physiological baseline (cytoplasm/turgor/membrane/nucleus ON) — assert particle count = N_old + n_filaments*beads_per_fil, bond count grows by backbone + crosslink, gamma_map covers if_bead, BAOAB does not HALT, one short run is numerically stable (no BAOAB guard…
- Contamination regression (the gate's core): gamma_cortical(IF ON) == gamma_cortical(IF OFF) to machine precision at the identical operating point; AND a forced-allow run where if_ is wrongly counted shows a large gamma jump (proves the denylist is load-bearing).
- Activation gate run (PRIMARY): OFF/ON delta-gamma on the equilibrated FA-adhered MCF7 build, myosin ON-vs-OFF active-gamma subtraction unchanged by IF presence within 5%.
- Activation gate run (SECONDARY): stretch a single built if_backbone chain across the Table law, record F(r), confirm finite extensibility in [2,3.5]x and stiffening ratio >=3x.
- CFL gate: assert dt <= cfl_safety_factor * gamma_if/k_bb at the cited soft modulus (k_bb~1e-3 N/m -> IF is NOT the binding CFL step; nucleus lamin still dominates); confirm IF does not lower the cell-wide dt headroom.
- Viz at closeout: extend the h{X}_vis.py entry-point with an IF cage render (perinuclear shell, backbone+crosslink bonds colored, the F(r) Table curve overlaid on the Kreplak/Block reference band) into outputs/h{X}/figs/, per the visualize-at-closeout rule.

### GPU / CFL

GPU: the LINEAR backbone+crosslink path is GPU-builtin (md.bond.Harmonic runs natively on the integrator device, CPU or CUDA — gpu_path_now=HOOMD_BUILTIN, device_resident_now=True in the registry). The build-time cKDTree crosslink search is one-shot CPU, NOT in the hot loop. The NONLINEAR path: md.bond.Table is ALSO a GPU-builtin bonded ForceCompute in HOOMD 7.0.1 (tabulated potentials run on-device), so the nonlinear law stays GPU-resident IF expressible as a Table — preferred over a per-step custom ForceCompute (which would force a GPU<->CPU sync and become a native-ForceCompute candidate; registry bottleneck_risk flags exactly this). No…


---

## Microtubules — `cell/microtubules.py`  *(Phase D)*

**Goal.** Activate the explicit MT aster as a default-OFF, physiologically-baselined compartment so the cell carries an explicit compression-load-bearing / centrosomal radial element, with its (over-stated) stiff-CFL risk resolved by the correct cytoplasm-viscosity drag and a constraint-route fallback, and a Dmitrieff-2017-style buckling/marginal-band observable gate.

- **Dependencies:** cell/cytoplasm.py (CytoplasmDrag) — must scale mt_bead/mtoc Stokes drag to η_eff=65.9 Pa·s; this self-consistency is what defuses the CFL (γ_b cytoplasm not water).; cortex/cortical_tension.py — the denylist must learn the 'mt_' prefix BEFORE any γ-with-MTs measurement, else cortical γ is contaminated (the registry already declares it; the…; cell/manifest.py + cell/cell.py — the resolve+extend+attach wiring (steps 4-6); the physiological baseline (cytoplasm/turgor/membrane/nucleus ON) is the state MTs are added FROM.; cell/dt_reconcile.py — per-compartment dt min must include the MT CFL so the global dt is auditable and provably non-binding at the chosen ℓ_0.; Native constrained-BAOAB M-SHAKE plugin (native/, project memory 'native-constrained-integrator') — ONLY if the PI chooses the inextensible-segment route C; reuses existing…; scripts/h2_single_filament.py — supplies the angle_k=EI/ℓ_0 bridge and the tangent-correlation L_p oracle the PRIMARY gate re-runs.; PI sign-off (gate-contract change) to flip the registry status EXPERIMENTAL→LIVE and set manifest_path — without it UnratifiedCompartmentError blocks loader activation (by design).
- **Ready when:** PI has (1) ratified Y_stretch=2.3e-7 N as an E·A KU OR approved the inextensible-segment constraint route (removing Y_stretch from the CFL); (2) confirmed L_mt (≈5µm) and n_mt for MCF7 or accepted the defaults; (3) signed off the registry EXPERIMENTAL→LIVE + manifest_path change AND the cortical_tension denylist 'mt_' addition; (4) the OFF/ON γ-contamination control passes <2% AND the PRIMARY L_p-recovery gate passes in-band on the full physiological baseline; with dynamic instability left DEFAULT-OFF (its in-vitro Walker-1988 rates flagged, its updater still NotImplementedError until a separate PI-gated rebuild task).
- **Effort:** ~2-3 days. Day 1: denylist fix + γ_b-cytoplasm wiring + manifest/cell.py/cytoplasm/dt_reconcile plumbing + unit tests (the module physics + sanity gate are already complete, so this is integration, not new physics). Day 2: OFF/ON γ-contamination integration control + PRIMARY L_p-recovery gate run (~501 particles, cheap) + optional buckling confirmatory run. Day 3 buffer: PI ratification round-trip on Y_stretch/L_mt/constraint-route + figures (aster render, CFL-vs-ℓ_0 plot showing the cytoplasm-viscosity defusing, L_p-recovery, buckling-vs-Euler) per the visualize-at-closeout rule. DI updater wiring is explicitly OUT of scope (separate PI-gated task).

### PI parameters to ratify

**`EI (flexural rigidity)`** — `SOLID_LITERATURE` / conf **HIGH** — *N·m²*
- Proposed: 2.2e-23
- Citation: Gittes, Mickey, Nettleton, Howard 1993, J. Cell Biol. 120(4):923-934
- Notes: κ=22 pN·µm². Independently re-used by Dmitrieff et al. 2017 PNAS 114(17):4418 (Cytosim MT model, κ=22 pN·µm²) — directly analogous fine-grained model. Module default _EI_MT_DEFAULT already set; satisfies EI=kT·L_p to <2%. NOT pending; listed for completeness.

**`L_p (persistence length)`** — `SOLID_LITERATURE` / conf **HIGH** — *m*
- Proposed: 5.2e-3
- Citation: Gittes 1993 JCB 120:923 (5.2 mm at room T)
- Notes: Band guard already in module [1.0e-3, 8.0e-3] m (Pampaloni 2006 PNAS 103:10248 reports length-dependent up to ~6-8 mm). Consistent. NOT pending.

**`Y_stretch (axial 1-D stretch modulus E·A) — THE dt LEVER`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *N*
- Proposed: 2.3e-7 (keep as harmonic-bond value) OR set inextensible via constraint route (preferred)
- Citation: Kis et al. 2002 PRL 89:248101 (E≈1.2 GPa); Pampaloni et al. 2006 PNAS 103:10248
- Notes: PI-PENDING. Verify Kis2002/Pampaloni2006 source_audit verdict=OK before any deliverable. Geometry cross-check confirms 2.3e-7 N is the right order: hollow tube r_out=12.5nm/r_in=9nm, A≈236 nm² ⇒ E·A≈2.8e-7 N (within 20% of the 2.3e-7 placeholder). KEY FINDING: at the physiological cytoplasm viscosity (65.9 Pa·s, NOT water) the stretch CFL is dt_str=cfl·γ_b·ℓ_0/Y ≈ 0.84 µs (ℓ_0=125nm) to 6.8 µs (ℓ_0=1µm) — COMPARABLE TO OR LARGER THAN the current production dt of 0.70 µs. So the explicit harmonic stretch bond…

**`γ_b (MT-bead Stokes drag) — physiological-baseline input`** — `DERIVED` / conf **HIGH** — *N·s/m*
- Proposed: 1.55e-5 (= 6π·65.9·12.5e-9)
- Citation: η_eff MCF7=65.9 Pa·s, Hu et al. 2024 Nanoscale Adv (PMC10929591), KU-3.B3.1; R_mtbead≈12.5 nm (MT tube radius, Gittes/Howard geometry)
- Notes: DERIVED, but the SINGLE most important correction: the host MUST pass γ_b computed at CYTOPLASM viscosity (65.9 Pa·s), not the cortex's water-based γ_b (3.9e-10). Per the physiological-baseline HARD rule MT beads are immersed in cytoplasm. This is what converts the CFL from catastrophic (water) to benign (cytoplasm). The cytoplasm.py CytoplasmDrag must scale mt_bead like the other immersed types.

**`ℓ_0 / beads_per_mt (segment length — discretisation knob, grid-invariant)`** — `DERIVED` / conf **HIGH** — *m*
- Proposed: ℓ_0 = 125e-9 (Dmitrieff/Cytosim choice) ⇒ beads_per_mt = round(L_mt/125nm)+1
- Citation: Dmitrieff et al. 2017 PNAS 114(17):4418 — MTs of L=9-16 µm 'finely represented with a segmentation of 125 nm' (Cytosim)
- Notes: Not a literature constant per se; it is grid-invariant (k_angle=EI/ℓ_0 and k_backbone=Y/ℓ_0 keep the continuum EI/EA fixed). 125 nm matches the canonical Cytosim MT discretisation and gives ℓ_0/r_g≈12.8, well above the bend↔stretch crossover r_g=√(EI/Y)≈9.8 nm. At ℓ_0=125nm and cytoplasm viscosity the stretch dt (0.84µs) just matches production dt — a comfortable choice. Larger ℓ_0 only relaxes the CFL further.

**`L_mt (MT contour length) — REQUIRED when enabled, host geometry`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *m*
- Proposed: 3.0e-6 to 5.0e-6 (a few µm; near R_cell=7.5µm for an aster that reaches the cortex)
- Citation: MCF7 R_cell=7.5µm (Wagner 2011); interphase MT lengths a few µm typical (Howard, Mechanics of Motor Proteins 2001)
- Notes: Resolver REQUIRES L_mt explicitly (raises if None — good). Set from host geometry, do NOT invent. For a cortex-reaching aster L_mt ≈ R_cell - R_nuc ≈ 5.6 µm. Euler buckling F_crit = π²EI/L² = 8.7 pN (5µm) to 24 pN (3µm) — confirms genuine compression load-bearing at cell scale (well above thermal). PI to confirm L_mt for MCF7 or surface a cell-type KU.

**`n_mt (aster arm count)`** — `ORDER_ESTIMATE` / conf **LOW** — *—*
- Proposed: 20 (module default); scale to ~50-250 for a dense interphase array if PI wants
- Citation: Order estimate; interphase centrosomal MT number is cell-type-variable (no MCF7 KU)
- Notes: Biology choice, not a discretisation knob (count does NOT grow with the ×40 actin scale — aster is centrosomal). 20×25=501 particles is cheap (P2). Flag: a true interphase array is denser; PI to set or accept the default-20 sparse aster as a first build.

**`Dynamic-instability rates v_g, v_s, f_cat, f_res`** — `NEEDS_MEASUREMENT` / conf **LOW** — *m/s, m/s, 1/s, 1/s*
- Proposed: v_g=3.3e-8 (2µm/min), v_s=2.8e-7 (17µm/min), f_cat=0.005, f_res=0.044
- Citation: Walker et al. 1988 J. Cell Biol. 107:1437-1448 (in-vitro purified tubulin)
- Notes: PI-PENDING + DEFAULT-OFF + updater NotImplementedError (honest). In-vitro rates; in-vivo/cancer rates differ greatly (MAPs, +TIPs, taxol). DI is NOT needed for the first activation gate (a stable aster is the baseline). Keep OFF; wiring the plus-end snapshot-rebuild is a separate PI-gated step. Listed so the PI knows the anchor provenance. A cell-type-specific MCF7 DI KU is genuinely a new measurement.

### Activation gate

- **Observable:** MT aster bending-rigidity reproduction (PRIMARY, intrinsic) AND optional marginal-band/buckling load-bearing (Dmitrieff-style, CONFIRMATORY). PRIMARY: thermal tangent-tangent correlation of a free MT chain recovers the input L_p=5.2 mm (i.e. the chain is rod-like over its few-µm contour, <r_ee²>/L_c² ≈ 1) — the H.2 persistence-length oracle re-run on mt_ chains. CONFIRMATORY: under an imposed end compression, the chain resists until the Euler force F_crit = π²EI/L² (8.7 pN at L=5µm, 24 pN at 3µm) then buckles, matching Euler within the discretisation (Dmitrieff 2017 used exactly this π²κ/L² check to validate the Cytosim MT). γ-CONTAMINATION: cortical γ measured with MTs ON must equal γ with…
- **Literature band:** PRIMARY: recovered L_p in the Gittes/Howard band [1.0, 8.0] mm (module _LP_MT_BAND_M), i.e. a few-µm MT is essentially straight (end-to-end/contour ratio > 0.999 over 5 µm since L_c/L_p ~ 1e-3). CONFIRMATORY: buckling onset within [0.7×, 1.3×] of Euler π²EI/L² (discretisation tolerance, Dmitrieff recovered Euler in Cytosim). γ-contamination band: |γ_ON − γ_OFF| / γ_OFF < 0.02.
- **Controls:**
  - OFF/ON differencing on cortical γ: build the FULL physiological baseline (cytoplasm 65.9 Pa·s, turgor 133 Pa, membrane_surface, nucleus ON) once with microtubules.enabled=false and once true; the cortical_tension 3-channel γ must be identical to <2% (proves…
  - OFF-identity: extend_snapshot_with_microtubules with enabled=False returns the input snap byte-identical (0 particles/bonds/angles); attach_microtubule_forces returns (None,None) and a pre-MT run is bit-for-bit reproducible.
  - Rigid/limit parity: in the limit Y_stretch→∞ (constraint route) OR k_backbone huge, the explicit-bond stretch CFL result must converge to the constraint-route bending-only CFL (cross-check dt_str→dt_bend as ℓ_0 fixed); and the bending response must be…
  - Force-free construction: built straight (every bond at ℓ_0, every angle at π) ⇒ total internal force = 0 at t=0 (Sanity Gate §3, already proven static).
  - No-contamination of the OTHER passive channels: turgor γ_passive and membrane tension unchanged ON vs OFF (the aster is internal, adds no surface).
- **PASS / PARTIAL / REFUTE:** PASS: L_p recovered in-band AND |γ_ON−γ_OFF|<2% AND (if compression run) buckling within [0.7,1.3]× Euler AND no global dt reduction forced below the run's chosen dt at the ratified ℓ_0. PARTIAL: L_p in-band and γ clean, but the explicit stretch bond forces dt below production (only if PI insisted on water-viscosity drag or ℓ_0<125nm — then activate via the constraint route instead, re-run). REFUTE: γ_ON−γ_OFF > 2% (denylist leak — the mt_…
- **No-γ-contamination check:** The cortical-γ estimator (cortex/cortical_tension.py) is denylist-mode: it sums EVERY bond type EXCEPT the adhesion prefixes. The mt_backbone / mt_bending types are NOT currently denylisted there — they would be COUNTED as cortical tension and contaminate γ. MANDATORY wiring step: add 'mt_' (and the explicit names mt_backbone, mt_bending — note bending is an ANGLE not a bond, but the rigid-λ…

### Wiring steps

1. 1. compartment_registry.py: the CompartmentSpec already exists (line 814, status=EXPERIMENTAL, manifest_path=None, denylist_bond_types=('mt_',)). To activate via the loader the PI must (a) flip status EXPERIMENTAL→LIVE and (b) set manifest_path=('compartments','microtubules') — this is a gate-contract change requiring PI sign-off (UnratifiedCompartmentError guards it today). Until then it stays…
2. 2. cortex/cortical_tension.py: add the MT load path to the denylist. Introduce 'mt_' to ADHESION_BOND_TYPE_PREFIXES (or, cleaner, rename the concept to NON_CORTICAL_BOND_TYPE_PREFIXES and add 'mt_'), sourced from microtubules.GAMMA_DENYLIST_PREFIX. THIS IS THE CONTAMINATION FIX and is mandatory before any γ measurement with MTs ON. Add a unit test asserting _is_adhesion_bond_type('mt_backbone')…
3. 3. configs/mcf7_baseline.yaml: add a compartments.microtubules block, enabled:false (default-OFF per the physiological-baseline code-hygiene convention — but note the rule that production must turn ON physiological compartments AT setpoint; MTs are NOT part of the resting cortical-γ baseline, so OFF is correct for the Gate-A/B suspended observable. Document that.). Keys: enabled, n_mt:20,…
4. 4. cell/manifest.py: import resolve_microtubules; add p_microtubules:Any|None=None to ResolvedBaseline; in resolve_baseline read comp.get('microtubules'), and if enabled call resolve_microtubules(cfg, kT=..., gamma_b=GAMMA_B_CYTOPLASM, cfl_safety_factor=..., mtoc_center=cell_centroid) — CRITICAL: pass γ_b computed at the CYTOPLASM viscosity 65.9 Pa·s and R_mtbead≈12.5 nm (=1.55e-5 N·s/m), NOT the…
5. 5. cell/cell.py: add enable_microtubules flag to CellBuildOptions; in Cell.build, after the lamellipodium/nucleus extension blocks, call snap = extend_snapshot_with_microtubules(snap, p_microtubules) (no-op when disabled, mirrors lamellipodium), and after the integrator is set call attach_microtubule_forces(sim, p_microtubules, cfl_strict=True). The MTOC center = cell centroid; for the…
6. 6. cell/cytoplasm.py (CytoplasmDrag): add 'mt_bead' and 'mtoc' to the set of immersed particle types whose Stokes drag scales to η_eff. This is what makes the γ_b passed in step 4 self-consistent with the BAOAB gamma_map and defuses the CFL.
7. 7. cell/dt_reconcile.py: register the MT CFL (p_microtubules.dt_cfl) in the per-compartment dt reconciliation so the global dt = min over all compartments is computed correctly and the MT contribution is visible/auditable (it should NOT bind at ℓ_0≥250nm / cytoplasm viscosity).
8. 8. DI updater (DEFER): MTDynamicInstability stays NotImplementedError / default-OFF. Wiring the plus-end snapshot-rebuild is a separate PI-gated integration step (coordinate with cell.py single-writer + BAOAB gamma_map). Not needed for the activation gate.

### Test plan

- Unit (extend module tests in tests/test_microtubules.py): OFF-identity — extend_snapshot_with_microtubules(snap, disabled) IS snap; attach returns (None,None). Already covered by Sanity Gate §2; add an explicit byte-identity assert on a built baseline snapshot.
- Unit: denylist coverage — assert cortical_tension._is_adhesion_bond_type('mt_backbone') and ('mt_bending') return True after step 2; assert a synthetic shell γ is unchanged when mt_ bonds are added then denylisted.
- Unit: γ_b cytoplasm wiring — assert resolve_microtubules receives γ_b≈1.55e-5 (cytoplasm) not 3.9e-10 (water) from the manifest path; assert dt_cfl_stretch≈8.4e-7 s at ℓ_0=125nm (the benign value), proving the physiological-baseline drag.
- Unit: CFL gate — build with dt slightly above p.dt_cfl and assert attach_microtubule_forces raises RuntimeError with the AXIAL-STRETCH binding message; with dt below, no raise.
- Unit: grid-invariance — refine ℓ_0 125→62.5 nm, assert k_angle=EI/ℓ_0 doubles but the recovered continuum EI (from a short static bend test) is invariant to <5%.
- Integration (the OFF/ON γ control): build the FULL physiological baseline (allow_no_nucleus=False, cytoplasm/turgor/membrane/nucleus ON) twice — MTs off vs on — and assert the 3-channel cortical γ matches to <2% (the contamination gate). Small N, short equilibration; CPU dev-fallback OK for the build, GPU for any stepping.
- Physics gate run 1 (PRIMARY): free MT chain, equilibrate, run the H.2 tangent-correlation oracle on mt_ chains; assert recovered L_p in [1,8] mm band. ~501 particles, short run.
- Physics gate run 2 (CONFIRMATORY, optional): impose end compression on a single 5µm MT chain, ramp force, assert buckling onset within [0.7,1.3]× Euler 8.7 pN (Dmitrieff-style check).
- Regression: a pre-MT production run with enable_microtubules=False is bit-for-bit identical (OFF-identity at the Cell.build level).

### GPU / CFL

CFL IMPACT (the critical deliverable), computed at the physiological baseline: the elastic forces are HOOMD built-ins (md.bond.Harmonic + md.angle.Harmonic) so GPU-native on CUDA — no port work for the mechanics (registry gpu_path_now=CPU is conservative; the built-ins are device-resident). The DI updater is CPU snapshot-rebuild, default-OFF, out of loop. The dt story, QUANTIFIED: r_g=√(EI/Y)=9.8 nm is the bend↔stretch crossover; at production ℓ_0 (≥125 nm) the AXIAL STRETCH term binds, as the module docstring states. BUT the absolute value depends on drag, and per the physiological-baseline rule the MT beads sit in CYTOPLASM (η=65.9 Pa·s ⇒…


---

## Dynamic osmotic / volume regulation — `cortex/osmotic_regulation.py`  *(Phase C)*

**Goal.** Turn on time-dependent osmotic volume regulation on the already-pressurised MCF7 baseline cell and show that a hyper/hypo-osmotic shock relaxes the reference volume V0(t) back toward osmotic equilibrium on the Hoffmann-2009 seconds-to-minutes RVD/RVI timescale, with NO contamination of the active cortical-gamma channel.

- **Dependencies:** enclosed_volume (HARD requires — the module mutates EnclosedVolumePressure.p.V0/turgor_dP0 in place; CompartmentSpec.requires=('enclosed_volume',)). enclosed_volume is already ON…; cytoplasm at 65.9 Pa·s, membrane_surface ON, turgor at 133 Pa — the physiological-baseline HARD rule: the RVD/RVI must be measured FROM the pressurised baseline, not from a…; cortex/cortical_tension.py — needed for the no-contamination control (active-γ ON vs OFF differencing over explicit cortex bonds).; The post-build attach must receive the SAME live EnclosedVolumePressure instance Cell.build created (object identity, not a copy) — a wiring dependency on cell/cell.py exposing…; PI sign-off on the Lp value (the named PI-pending decision) BEFORE production — recommend Lp=1.0e-12 m/(s·Pa) for AQP5-expressing MCF7 (Jung 2011), with 1e-13 retained only as the…; PI awareness of the Option-A trajectory-vs-reference gap: the implemented constant-dP_target law makes V0(t)'s OWN half-time = τ_Kvol·ln2 (minutes-to-hours), while the…
- **Ready when:** "PASS conditions met: (1) PI ratifies Lp (recommend 1.0e-12 m/(s·Pa) for MCF7 per Jung-2011 AQP5; current 1e-13 default flagged as likely-wrong-for-MCF7); (2) the manifest block + resolve_/attach_ wiring lands and the OFF-identity + resolver + analytical tests pass; (3) the short integrated RVD/RVI run shows τ_RVD ∈ [3 s, 600 s] across the Lp sweep (analytically guaranteed in-band, 3.1-32.3 s) with correct shock sign; (4) the no-contamination control shows active γ=γ(myoON)−γ(myoOFF) UNCHANGED ON vs OFF and identical cortex bond inventory; (5) figure emitted to outputs/h7/figs/ and REPORT.md §Figures updated. The CompartmentSpec status is then promoted EXPERIMENTAL→LIVE under PI sign-off. NOT-ready / surface-to-PI: if the trajectory half-time (τ_Kvol·ln2) is required to be in-band (→ Option A code change), or if a measured MCF7 Lp is demanded over the Olbrich/AQP5 estimate."
- **Effort:** "Low-to-moderate, ~0.5-1 day. The module, its sanity tests, and its CompartmentSpec already exist and are mechanism-complete; this is wiring + a gate run, not new physics. Breakdown: ~1-2 h wiring (manifest block + resolve_baseline stanza + ResolvedBaseline field + post-build attach_ + CompartmentSpec.manifest_path) — the turnover/SubstrateLigandPin patterns are the templates; ~1 h extending tests (resolver wiring test + τ-band assertions); ~1-2 h wall-clock for the gate runs (5 Lp points × 2 shock signs run concurrently on gbook, ≤16 parallel) + the no-contamination control; ~0.5 h auto-viz + REPORT. Blocking item is the single PI decision on Lp (one number), which I have pre-bracketed with a concrete recommended value and a hard MCF7-specific anchor (AQP5)."

### PI parameters to ratify

**`Lp (membrane hydraulic permeability)`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *m/(s·Pa)*
- Proposed: 1.0e-12 (MCF7 production, AQP-expressing high end); band [1.8e-13, 1.1e-12]; current code default 1.0e-13 (aquaporin-poor LOW end) is most likely WRONG for MCF7
- Citation: Olbrich, Rawicz, Needham & Evans 2000, Biophys J 79:321 (Pf 25-150 µm/s → Lp via Lp=Pf·Vw/(R_gas·T), Vw=18e-6 m³/mol → [1.8e-13,1.1e-12]); MCF7-specificity from Jung, Park, Jeon & Kwon 2011, PLoS ONE 6(12):e28492, doi 10.1371/journal.pone.0028492 (PubMed — MCF7 expresses FUNCTIONAL aquaporin AQP5…
- Notes: DECISION-CRITICAL and the named PI-pending item. No direct MCF7 Lp datum exists, but the AQP5 finding (Jung 2011) is a hard qualitative anchor that the module's stated default rationale ('aquaporin-poor membrane', 1e-13 low end) contradicts. Recommend PI ratify Lp=1.0e-12 m/(s·Pa) for MCF7 production (τ_RVD≈3.1 s, fast end of Hoffmann band) and keep 1e-13 only as the documented AQP-poor lower bound for a sensitivity sweep. Lp enters ONLY the timescale (τ_RVD=V0/(Lp·A·Π_osm)); it is grid-invariant. A 6× Lp range…

**`c_phys (intracellular osmolarity, van't Hoff modulus anchor)`** — `SOLID_LITERATURE` / conf **HIGH** — *mol/m³ (= mOsm/L)*
- Proposed: 300 (code default; band 290-300)
- Citation: Lodish, Molecular Cell Biology; Alberts et al., MBoC — standard mammalian intracellular/extracellular osmolarity ~290-300 mOsm/L. Corroborated qualitatively by Stewart 2011 Nature (rounding-pressure osmotic engine) and Stroka 2014 Cell (osmotic-engine model, J=Lp(ΔP−ΔΠ))
- Notes: Sets Π_osm = R_gas·T·c_phys ≈ 7.74e5 Pa at 310 K (confirmed by direct compute). Already defaulted and grid-invariant; no PI action needed unless a tighter MCF7 osmometer datum is wanted. Π_osm is the restoring stiffness of the REFERENCE osmometer mode (τ_RVD), not of the implemented K_vol-driven trajectory.

**`delta_c (imposed osmolar shock — the EXPERIMENT knob, not a baseline constant)`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *mol/m³*
- Proposed: ±100 to ±150 mOsm step (e.g. +100 mM sorbitol hyperosmotic → RVD; -100 mOsm hypotonic → RVI). Use exactly the Jung-2011 MCF7 protocol magnitude (100 mM sorbitol) for the hyperosmotic arm so the perturbation is MCF7-anchored
- Citation: Jung 2011 PLoS ONE (100 mM sorbitol on MCF7, doi 10.1371/journal.pone.0028492); Hoffmann, Lambert & Pedersen 2009 Physiol Rev 89:193 (RVD/RVI evoked by 30-50% osmotic steps); Stewart 2011 Nature (±100-200 mOsm steps, ~3 min recovery)
- Notes: This is the gate's stimulus, set per-run in cfg['delta_c']; it is NOT a magic number and NOT a baseline setpoint. Keep it inside the physiological shock range the literature uses (≤ ~50% osmotic change) so V/V0 never approaches the V0_min singular guard (0.01·V0_ref, R-floor 0.215·R, unreachable).

**`temperature_K`** — `SOLID_LITERATURE` / conf **HIGH** — *K*
- Proposed: 310.15 (37 °C, code default = physiological)
- Citation: Physiological body temperature; van't Hoff factor. Consistent with the resolver default. Olbrich Pf→Lp band reported at 298 and 310 K — use the 310 K conversion [1.75e-13, 1.05e-12] for consistency
- Notes: No action. Affects Π_osm linearly (and thus τ_RVD); keep at 310.15 for production, NOT the 298 K used in some bilayer Lp conversions.

**`A_mem (water-exchange membrane area)`** — `DERIVED` / conf **HIGH** — *m²*
- Proposed: 4π·R_cell² = 7.069e-10 m² at R_cell=7.5e-6 (sphere default; computed)
- Citation: Geometric; defaults to the cell sphere area in the resolver. Cross-check vs membrane_surface compartment area if that is enabled
- Notes: Grid-invariant geometric default. If the membrane_surface shell is ON in the baseline (it is, per physio-baseline), consider passing its measured area for self-consistency, but the sphere area is within a few % and acceptable. A_mem enters τ_RVD and τ_Kvol identically.

**`batch_steps (slow-water-mode integration stride)`** — `DERIVED` / conf **HIGH** — *steps*
- Proposed: 1000 (code default) → batch_dt = 6.98e-4 s at dtc=6.98e-7 s; passes the slow-mode CFL by ~7 orders (batch_dt ≪ τ_Kvol≈1.9e4 s and ≪ τ_RVD≈3-32 s)
- Citation: Derived from the resolver slow-mode CFL gate (batch_dt < τ_Kvol); not a literature constant
- Notes: Numerical resolution knob, grid-invariant CFL-gated. For a clean V0(t) trajectory sampled against a 3-32 s observable, 1000-5000 steps/tick gives O(1e3) ticks over the relaxation — ample. The resolver AND the Updater __init__ both re-assert the CFL, so an under-resolved stride raises ValueError.

### Activation gate

- **Observable:** Volume-regulation relaxation half-time t½ read from the V0(t) trajectory after an imposed osmotic shock (delta_c). TWO numbers are reported and they are DIFFERENT modes, both must be checked: (1) τ_RVD = V0/(Lp·A·Π_osm) — the REFERENCE free-osmometer timescale, the one that must land in the Hoffmann band; (2) τ_Kvol = V0/(Lp·A·K_vol) — the eigen-timescale of the mode the constant-target law ACTUALLY integrates (the code's own V0(t) half-time = τ_Kvol·ln2). At MCF7 baseline (computed): τ_RVD = 3.1 s (Lp=1.1e-12) to 32.3 s (Lp=1e-13); τ_Kvol = 1.7e3 s to 1.9e4 s. Secondary observable: the active cortical-gamma number γ(myo ON)−γ(myo OFF) must be UNCHANGED to within noise vs the…
- **Literature band:** τ_RVD ∈ [seconds, minutes] = [~3 s, ~600 s] per Hoffmann, Lambert & Pedersen 2009 Physiol Rev 89:193; corroborated by Stewart 2011 Nature (~3 min ≈ 180 s RVD recovery in HeLa metaphase) and the Jung-2011 MCF7 sorbitol protocol. Target sub-band for AQP-rich MCF7 (Lp≈1e-12): t½ ≈ τ_RVD·ln2 ≈ 2-25 s. Π_osm must be 7.5e5-8.0e5 Pa (computed 7.74e5 at 300 mOsm/310 K).
- **Controls:**
  - OFF/ON differencing: build the SAME baseline cell with osmotic_regulation OFF (V0 static) and ON; the ON V0(t) must relax while OFF V0 is flat to machine precision (OFF-identity already unit-tested).
  - Rigid/limit parity: in the Lp→0 limit (or batch_steps→∞ relative to run) dV/dt→0 — V0 frozen, recovers the static enclosed_volume result (OFF-identity).
  - Fixed-point control: with delta_c=None (no shock, dP_target=resting turgor 133 Pa) and the cell already at its operating volume, dV/dt≈0 — NO spurious setpoint drift (already unit-tested as the ΔP_mech==ΔP_target fixed point).
  - No-contamination check: GAMMA_DENYLIST_PREFIX='' is correct BECAUSE this module adds zero bonds — assert the cortical_tension estimator's bond inventory is byte-identical ON vs OFF, and that active γ = γ(myoON)−γ(myoOFF) is unchanged within stochastic noise.…
  - Sign/sense control: hyperosmotic shock must drive V0 DOWN (RVD, water leaves); hypotonic up (RVI). Wrong sign = REFUTE.
  - Lp sensitivity sweep: run Lp ∈ {1e-13, 1.8e-13, 5e-13, 1.05e-12, 1.1e-12} and confirm t½ tracks the predicted τ_RVD (32.3, 18.0, 6.5, 3.1, 2.9 s) — establishes that the SINGLE free parameter Lp moves the observable monotonically and stays in-band across the…
- **PASS / PARTIAL / REFUTE:** PASS: τ_RVD lands in [3 s, 600 s] (it does, 3.1-32.3 s across the full Olbrich Lp band — analytically guaranteed at MCF7 baseline), the measured V0(t) half-time equals τ_Kvol·ln2 to <1% (closed-form Euler reproduction, already an analytical test), shock sign is correct (RVD down / RVI up), and active γ is UNCHANGED ON vs OFF. PARTIAL: τ_RVD in-band but the IMPLEMENTED trajectory half-time (τ_Kvol·ln2 ≈ 20-220 min) is what an observer of V0(t)…
- **No-γ-contamination check:** Module adds 0 particles / 0 bonds / 0 forces (GAMMA_DENYLIST_PREFIX='' is correct by construction, not by denylisting). Verify: (a) cortical_tension estimator's explicit-cortex bond set is identical ON vs OFF; (b) the only mutated state is EnclosedVolumePressure.p.V0 (and optionally turgor_dP0) — a passive channel, denylisted from the active-γ subtraction already; (c) the Updater's act() leaves…

### Wiring steps

1. 1. Manifest block: add an `optional_subsystems.osmotic_regulation` block to a NEW recipe config (e.g. aleph/configs/recipes/osmotic_rvd.yaml that extends mcf7_baseline.yaml) — NOT to mcf7_baseline.yaml itself (it must stay default-OFF for all other production). Block carries: enabled:true, base_config:mcf7_baseline.yaml, Lp:1.0e-12 (PI-ratified value), batch_steps:1000, delta_c:<shock>…
2. 2. resolve_*: in aleph/cell/manifest.py resolve_baseline(), after the existing `tov_b = opt.get('turnover')` stanza, add an `osmo_b = opt.get('osmotic_regulation')` stanza that, if _enabled, calls resolve_osmotic_regulation(_opt_cfg(osmo_b), R_cell=p_cortex.R_cell, dt=dtc, p_enclosed_volume=p_enclosed_volume, temperature_K=...). CRITICAL ordering: it MUST run AFTER p_enclosed_volume is resolved…
3. 3. CellBuildOptions flag: add a build flag (e.g. `osmotic_regulation: bool=False` or thread it through the existing optional-subsystem mechanism) so Cell.build / build_baseline_cell knows to attach the updater. Since this module attaches POST-build (it needs the live EnclosedVolumePressure force object, not the snapshot), it does NOT need an _extend_snapshot_* hook — it is an attach_* only.
4. 4. attach_*: in cell/cell.py (or build_baseline_cell), AFTER the simulation is built and the enclosed-volume force is in sim.operations, call attach_osmotic_regulation_to_simulation(sim, p_osmotic, ev_force) — passing the SAME live EnclosedVolumePressure instance Cell.build created (the module mutates its .p.V0 in place). Mirror the existing post-build attach pattern (SubstrateLigandPin /…
5. 5. CompartmentSpec: the spec ALREADY EXISTS in cell/compartment_registry.py (lines 843-872, status=EXPERIMENTAL, requires=('enclosed_volume',), gamma_contaminating=False). Only edit: set manifest_path from None to ('optional_subsystems','osmotic_regulation') once step 1's block is authored, and promote status EXPERIMENTAL→LIVE only after the gate PASSes (PI sign-off — status change is…
6. 6. No GAMMA_DENYLIST edit needed: GAMMA_DENYLIST_PREFIX='' is correct (zero bonds added). Do NOT add a prefix.
7. 7. Auto-viz: per the production-driver-auto-viz rule, the driver must subprocess to a *_vis.py to emit the V0(t)/τ figure into outputs/h7/figs/.

### Test plan

- A. STATIC/analytical (already largely present in aleph/tests/test_osmotic_regulation.py): confirm OFF-identity, dimensional analysis, sign-sense (RVD down / RVI up), one-Euler-step closed-form reproduction (water_flux_volume_step matches the analytic ΔV to float precision), and the slow-mode CFL gate raises on under-resolution. Extend…
- B. Resolver wiring test: load the new osmotic_rvd.yaml recipe through resolve_baseline and assert p_osmotic.enabled, p_osmotic.K_vol==p_enclosed_volume.K_vol (copied), p_osmotic.V0_ref==p_enclosed_volume.V0, dP_target default==turgor_dP0 (133 Pa) when delta_c unset, and the CFL passes at batch_steps=1000.
- C. SHORT integrated run (the activation gate itself), suspended/rounded MCF7 at the FULL physiological baseline (cytoplasm 65.9 Pa·s, turgor 133 Pa, membrane_surface ON; nucleus may use the sanctioned allow_no_nucleus for the cortical-tension observable). Two arms: hyperosmotic (delta_c≈-100 mOsm sorbitol, Jung-2011 magnitude) and…
- D. NO-CONTAMINATION run (the decisive control): identical baseline, osmotic_regulation OFF vs ON, both with myosin ON and OFF; measure active γ=γ(myoON)−γ(myoOFF) via cortex/cortical_tension.py over EXPLICIT cortex bonds. Assert the active-γ number is unchanged ON vs OFF within ensemble noise and the cortex bond inventory is…
- E. Auto-viz: extend aleph/scripts/h1_h2_vis.py (or a h7_osmotic_vis.py) to plot per-realisation V0(t)/V0_ref thin lines + ensemble mean, with the Hoffmann seconds-to-minutes band shaded, τ_RVD and τ_Kvol·ln2 annotated, SI units, no axis truncation. Emit to outputs/h7/figs/h7_osmotic_rvd.png and reference it in REPORT.md §Figures.

### GPU / CFL

"GPU: builtin / nothing to port. The Updater is a P2 batch hoomd.custom.Action that reads the enclosed-volume force's cached scalar diagnostic (last_pressure, already float64) and writes two scalars per tick — NO per-particle loop, NO cpu_local_snapshot, NO GPU↔CPU sync in the hot path (uses_cpu_local_snapshot=False in the contract). It is agnostic to whether EnclosedVolumePressure runs on the CPU Custom force or its enclosed_volume_gpu.py/CuPy backend. Native ForceCompute candidate = NO (it is a boundary-condition scheduler, not a force). CFL: TWO timescales, both far above batch_dt. At the MCF7 baseline (computed): dtc=6.98e-7 s,…


---

## Membrane reservoir / bleb — `cell/membrane_reservoir.py`  *(Phase C)*

**Goal.** Activate the explicit membrane-cortex breakable-tether mesh + area reservoir so a bleb nucleates as a real topology change (Bell-Evans mem_tether rupture) when local hydrostatic detachment load exceeds membrane-cortex adhesion, and the reservoir buffers in-plane tension before the area-elastic term engages — measured FROM the physiological baseline (turgor 133 Pa, cytoplasm 65.9 Pa.s, nucleus + membrane_surface ON), with mem_ bonds denylisted out of gamma_soft.

- **Dependencies:** membrane_surface (HARD requires=('membrane_surface',)): the tethered membrane beads and the reference area A0 the reservoir grows ARE the membrane_surface shell + A0. Enforced at…; cortex/cortical_tension.py (BLOCKER-0, non-owned): must be wired to honour the registry 'mem_' denylist before enablement; currently hardcoded to FA-clutch prefixes only.; Own mem_node bead layer (BLOCKER-1): radially-offset membrane bead layer is a PREREQUISITE for a non-degenerate tether mesh because membrane_surface currently rides the cortex…; cortex/erm.py (conceptual coupling): the ERM radial tether is the lab-frame analog of the per-bead mem_tether; k_tether is the per-bead form of k_ERM. Avoid double-counting the…; Physiological baseline compartments (cytoplasm 65.9 Pa.s, enclosed_volume/turgor 133 Pa, nucleus) must be ON — the bleb/buffering observables are only valid measured FROM that…; Shared md.bond.Harmonic in cell.py (the cortex/xlink/myosin/FA force): mem_tether MUST be registered on it, not a second appended harmonic force (else IncompleteSpecificationError…
- **Ready when:** "READY for the STATIC-mechanism milestone NOW (resolver + builder + builtin-force attach + 21 static tests exist; CompartmentSpec registered). NOT READY for ENABLEMENT until: (1) BLOCKER-0 — cortical_tension honours the registry 'mem_' denylist and the dynamic no-contamination control (active-gamma invariant under OFF/ON toggle) PASSES; (2) BLOCKER-1 — the membrane has its own radially-offset mem_node layer so the tether mesh is non-degenerate; (3) PI ratifies sigma_crit_bleb (the cross-line ~1e-4 N/m anchor + band) AND the cortical-tension -> radial-detachment-force bridge, before the Bell-Evans act() rupture loop is written; (4) PI ratifies f_excess (MCF7 value or the [0.02,0.40] sweep band). The static tether-mesh + reservoir-buffering path can ship and pass its gate ahead of (3)/(4); the bleb-nucleation path stays declare_pending behind (3)."
- **Effort:** "Denylist wiring + regression test (BLOCKER-0): ~0.5 day (one non-owned-file change + PI sign-off). Own mem_node radial bead layer + non-degenerate tether build + _extend_snapshot hook (BLOCKER-1): ~1.5-2 days (geometry + tag-range plumbing + tests). Manifest/config/CellBuildOptions plumbing: ~0.5 day. Static-mechanism + reservoir-buffering gate run on gbook GPU (OFF/ON differenced) + auto-viz: ~1 day wall (run) + 0.5 day analysis. Bell-Evans act() rupture loop + bleb-nucleation gate (gated behind PI sigma_crit + force-bridge): ~2 days once unblocked. TOTAL to static+buffering activation ~3-4 days; +2 days for the bleb path after PI anchors. Pacing note (PI rule): report a wall-clock ETA at each gbook launch."

### PI parameters to ratify

**`sigma_crit_bleb`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *N/m*
- Proposed: 1.0e-4 (0.1 mN/m = 100 pN/um); propose band [4e-5, 4e-4] N/m for the gate
- Citation: Tinevez, Schulze, Salbreux, Roensch, Joanny, Paluch (2009) PNAS 106(44):18581-18586, doi:10.1073/pnas.0903353106 (PMID 19846787, via PubMed). Cross-confirmed by Dmitrieff et al. (2017) PNAS, aleph/references/dmitrieff-...pdf: 'actomyosin cortex of blebbing cell sigma ~ 100 pN/um (25=Tinevez2009)'.
- Notes: PI-PENDING. Tinevez's critical tension is the threshold BELOW which a laser-nucleated bleb cannot expand; their cells are MOUSE melanoma-type lines (MeSH: Cell Line/Mice), NOT MCF7 and NOT human — so this is a cross-line ORDER anchor, not an MCF7 datum. The canonical actomyosin-cortex blebbing tension ~100 pN/um (1e-4 N/m) is the mid-anchor; Tinevez modulated tension over ~3-4x by myosin/actin-turnover drugs, motivating the [4e-5,4e-4] band. CRITICAL CONSISTENCY CHECK the PI must ratify: this is a CORTICAL TENSION…

**`f_excess`** — `ORDER_ESTIMATE` / conf **LOW** — *dimensionless (fraction of A0)*
- Proposed: 0.1 (10%); propose band [0.02, 0.40] for the gate, MCF7 point value PI-pending
- Citation: Raucher & Sheetz (1999) Biophys J 77(4):1992-2002, doi:10.1016/S0006-3495(99)77040-2 (PMID 10512819) — membrane reservoir buffers tether force at constant length over microns of pulled tube. Figard & Sokac (2014) Bioarchitecture 4:39-46, doi:10.4161/bioa.29069 + Figard et al. (2013) Dev Cell…
- Notes: PI-PENDING. The literature gives a RANGE 'a few % to tens of %' (folds+microvilli+caveolae together), not a single MCF7 number; caveola+microvillar reservoirs can buffer up to ~30-40% area before depletion (Sinha 2011 / Figard 2014). No MCF7-specific measurement exists in the corpus or KB (confirmed by tag_query: no f_excess record). 10% is a conservative mid-of-low-band default for the gate; the [0.02,0.40] band brackets the physiological span. PI must either anchor an MCF7 value (ideally from an MCF7…

**`W_MCA`** — `SOLID_LITERATURE` / conf **HIGH** — *J/m^2*
- Proposed: 1.0e-5 (DEFAULT_W_MCA, mid-band); band [1e-6, 1e-4]
- Citation: KU-3.B1.4 (already in module, W_MCA_BAND): Hochmuth et al. (1996) Biophys J 70:358; Derenyi/Julicher/Prost (2002) PRL 88:238101; Brochard-Wyart et al. (2006) PNAS 103:7660; Diz-Munoz/Fletcher/Weiner (2013) Trends Cell Biol 23:47. KB KB-3.B1.1 (DizMunoz2013, PontesGauthier2017): apparent membrane…
- Notes: NOT PI-pending — already anchored and band-guarded in resolve_membrane_reservoir. Sets the tether adhesion energy E=W_MCA*A_bead and rupture force F_c. Carried here so the PI sees the full constant set; no action needed beyond keeping it in-band.

**`k_tether`** — `DERIVED` / conf **HIGH** — *N/m*
- Proposed: 0.1 (default; = ERM radial pinning order)
- Citation: Charras et al. (2008) Biophys J 94:1836 doi:10.1529/biophysj.107.113605 (KB SE116 Charras2008_BJ, k_ERM ~ 0.1 N/m, anchors KB-3.1/KB-3.18). Same order as cortex/erm.py k_ERM.
- Notes: NOT PI-pending. The per-bead bond form of the membrane-cortex ERM link. Already defaulted in the resolver; feeds the CFL gate dt <= alpha*gamma_b/k_tether and the rupture force eq. R. Confidence HIGH on order; exact per-bead stiffness is a discretization of the measured k_ERM.

**`max_tether_dist`** — `SOLID_LITERATURE` / conf **HIGH** — *m*
- Proposed: 200e-9 (200 nm cortex thickness)
- Citation: KU-3.17 cortex thickness ~200 nm (Clark et al. 2013; Salbreux/Charras/Paluch 2012 Trends Cell Biol 22:536). Already in resolver default.
- Notes: NOT PI-pending; physical acceptor reach (membrane sits within a tether length of the ~200 nm cortex). NB: in the CURRENT default build membrane and cortex share the SAME shell tags, so the meaningful acceptor distance is ZERO until the membrane gets its own radially-offset bead layer (see wiring_steps blocker).

### Activation gate

- **Observable:** TWO coupled observables, mirroring the H.7 Gate-A/B contract style. (1) BLEB-NUCLEATION THRESHOLD: with the rupture updater LIVE, the fraction of mem_tether bonds that rupture per unit time as a function of imposed cytoplasmic Laplace pressure P = 2*sigma/R_cell, sweeping sigma across the sigma_crit_bleb band. The nucleation onset (first sustained rupture cluster) must occur near sigma_crit_bleb, i.e. when the per-tether detachment load P*A_bead crosses F_c = sqrt(2 k_tether W_MCA A_bead). (2) RESERVOIR TENSION-BUFFERING: with reservoir release LIVE, in-plane membrane tension gamma_membrane (from membrane_surface) measured vs imposed area strain dA/A0; tension must stay near-constant…
- **Literature band:** (1) Bleb threshold: nucleation onset at sigma in [4e-5, 4e-4] N/m, anchor ~1e-4 N/m (Tinevez 2009 PNAS, cross-line order anchor; Dmitrieff 2017 PNAS confirms ~100 pN/um for blebbing actomyosin cortex). (2) Reservoir plateau: tension-buffered area strain in [2%, 40%] before depletion (Raucher-Sheetz 1999 BJ; Sinha/Nassoy 2011 Cell caveola reservoir; Figard 2014 Bioarchitecture); plateau…
- **Controls:**
  - OFF/ON differencing (MANDATORY): run membrane_reservoir OFF vs ON at otherwise-identical physiological baseline; the bleb/buffering signal is the ON-minus-OFF difference, never the absolute ON number.
  - Rigid/limit parity: in the W_MCA -> high (1e-4) limit no tether should rupture below the highest band sigma (blebbing is suppressed by strong adhesion) — F_c monotonic-increasing in W_MCA is already a STATIC sanity test; the LIVE run must reproduce the trend…
  - NO-CONTAMINATION check (HARD, the load-bearing control): measure gamma_soft (cortical_tension estimator) with membrane_reservoir OFF and ON; the active-gamma number gamma(myosin ON)-gamma(myosin OFF) MUST be invariant to within seed noise across the OFF/ON…
  - Conservation control: ON build must show ΔN_particles = 0 vs OFF in the default (no-own-layer) build, or = +n_mem if the radially-offset own-layer variant is used; net momentum from the symmetric harmonic tether = 0 (Newton 3rd).
  - Reservoir mass-balance control: total released area Sigma A_release must equal f_excess*A0 exactly (no area created/destroyed beyond the reservoir budget).
- **PASS / PARTIAL / REFUTE:** PASS: (1) bleb-nucleation onset falls inside [4e-5,4e-4] N/m AND scales correctly with W_MCA (rigid/limit parity holds), AND (2) reservoir produces a tension plateau over an area strain inside [2%,40%] then the expected post-depletion rise, AND (3) the no-contamination control shows active-gamma invariant under the membrane_reservoir OFF/ON toggle within seed noise. PARTIAL: the tether/topology + reservoir mechanism runs and the no-contamination…
- **No-γ-contamination check:** The cortical-gamma estimator MUST drop every mem_ bond. Precondition (STATIC, already in test_membrane_reservoir.py): every created bond type starts with GAMMA_DENYLIST_PREFIX='mem_'. Guarantee (DYNAMIC, the actual gate control): after wiring cortical_tension to honour the registry denylist, assert gamma_soft is bit-identical with mem_tether bonds present vs absent on a built cell, and assert…

### Wiring steps

1. BLOCKER 0 (PI-gated, non-owned file): wire aleph/cortex/cortical_tension.py to honour the registry denylist. Today ADHESION_BOND_TYPE_PREFIXES (L127) is hardcoded to ('fa_actin_clutch',) and _is_adhesion_bond_type (L147) only checks that. Extend the denylist to read CompartmentSpec.denylist_bond_types from compartment_registry (which already carries ('mem_',) for membrane_reservoir), OR append…
2. BLOCKER 1 (geometry, the real physics blocker): give the membrane its OWN radially-offset bead layer. membrane_surface currently rides the SAME cortex shell tags [0, n_cortex_actin) (membrane_surface.py L18, cell.py L1378 shell_tag_range=(0,n_cortex_actin)). So in the default build membrane_tag_range == cortex_tag_range and build_membrane_tethers' nearest-cortex query returns each bead's nearest…
3. Add an optional_subsystems.membrane_reservoir block to aleph/configs/mcf7_baseline.yaml (enabled:false by default; W_MCA, k_tether, max_tether_dist, sigma_crit_bleb:null, f_excess:null, batch_steps). Set CompartmentSpec.manifest_path from None to ('optional_subsystems','membrane_reservoir').
4. Add resolve_membrane_reservoir(...) call into aleph/cell/manifest.resolve_baseline (mirror resolve_erm at L246: pass R_cell=p_cortex.R_cell, dt). Gate on the membrane_surface-enabled assertion already present at manifest.py L175-181 (requires=('membrane_surface',) is enforced there).
5. Add a _extend_snapshot_membrane_reservoir hook in cell.py that calls build_membrane_tethers(snap, p, membrane_tag_range=<mem_node range>, cortex_tag_range=(0,n_cortex_actin)) AFTER the mem_node layer exists; thread the returned MembraneTetherLayout to attach.
6. Add attach_membrane_tether_force(sim, p, layout, gamma_b=p_cortex.gamma_b, ...) in the cell.py attach phase. CRITICAL: it must register mem_tether on the EXISTING shared md.bond.Harmonic (the function already does this, L791-796) to avoid IncompleteSpecificationError. Add a CellBuildOptions flag (e.g. enable_membrane_reservoir).
7. DEFER the rupture updater: MembraneTetherUpdater.__init__ raises NotImplementedError UNCONDITIONALLY (L865) until the act() Bell-Evans rupture loop is written AND sigma_crit_bleb is PI-anchored. Wiring the static tether mesh + reservoir is fully separable from enabling bleb nucleation; ship the static/buffering path first, gate the rupture path behind PI sign-off on sigma_crit_bleb + the…
8. FLAG registry/module drift to PI (non-blocking): the CompartmentSpec performance_contract says per_step_force=True and uses_cpu_local_snapshot=True, but the actual default build adds NO per-step Python force (mem_tether is a builtin md.bond.Harmonic) and the builder uses no cpu_local_snapshot. Reconcile the contract to the module (per_step_force=False for the tether; updater is per-batch only).

### Test plan

- STATIC (already present: aleph/tests/test_membrane_reservoir.py, 21 tests) — re-run and confirm: dimensional checks (E=W_MCA*A_bead [J]; F_c=sqrt(2 k W A) [N]; Delta_c=F_c/k [m]); enabled=False no-op (snap unchanged, ΔN=0); n_anchored=0 -> zero tethers; W_MCA out-of-band -> ValueError; sigma_crit_bleb None -> updater…
- NEW STATIC: degenerate-self-tether guard — assert that when membrane_tag_range == cortex_tag_range the builder either (a) excludes self-pairs (tag_i != tag_j) or (b) raises, so the BLOCKER-1 geometry bug cannot silently produce zero-length tethers.
- NEW DYNAMIC (denylist guarantee): on a built baseline cell, assert measure_cortical_tension gamma_soft is bit-identical with vs without live mem_tether bonds; and assert active-gamma = gamma(myo ON)-gamma(myo OFF) invariant under membrane_reservoir OFF/ON toggle within seed noise. This is the no-contamination gate control; it can only…
- NEW DYNAMIC (CFL): assert attach_membrane_tether_force raises on dt > 0.1*gamma_b/k_tether and passes at the production dt (mirror erm.py/crosslinkers.py CFL gate).
- NEW DYNAMIC (force-free construction): build tethers at the construction separation, run 0 steps, assert mem_tether bond energy U=0 and net force on the membrane+cortex subsystem = 0 (Newton 3rd / no spurious momentum).
- NEW DYNAMIC (reservoir mass-balance): with f_excess PI-anchored, assert Sigma A_release == f_excess*A0 and that growing A0 by A_release reduces gamma_membrane along the membrane_surface tension-area curve (buffering sign).
- GATE RUN (production, gbook GPU): the two-observable activation gate above, OFF/ON differenced, on the physiological baseline (turgor 133 Pa, cytoplasm 65.9 Pa.s, nucleus + membrane_surface ON). Auto-generate figures per the auto-viz rule into outputs/h8/figs/ (bleb-onset vs sigma; tension-plateau vs area strain; gamma_soft OFF/ON…

### GPU / CFL

"GPU: the mem_tether static force is a HOOMD builtin md.bond.Harmonic (GPU-resident, native ForceCompute) registered on the SHARED bond.Harmonic — no Python per-step custom force, fully GPU-main compliant. One-shot cKDTree acceptor query at construction only (O(n_mem log n_cortex)); no broad-phase in the hot loop (tether pairs fixed once built). The (deferred) MembraneTetherUpdater is a per-BATCH Action that would host-sync via sim.state.get_snapshot() like the D2 crosslinker/integrin updaters — CPU until the GPU-main binder port lands; it is rare (every batch_steps=100) and P2 priority, so low cost. CFL: (1) tether harmonic relax time tau =…


---

## Cadherin cell–cell junction — `junction/cadherin.py`  *(Phase E)*

**Goal.** Stand up the first TWO-CELL build (two cortex-baseline cells in one box with cross-interface cadherin partner search), add a junction_tension observable that reuses the cortical method-of-planes machinery over cadherin_trans bonds, register the Iturri-2020 SourceEvidence anchor, and validate the doublet against the Sim-2015 1-5 pN homeostatic junction-tension band — all FROM the physiological baseline, with cadherin_* bonds denylisted from cortical gamma.

- **Dependencies:** validation/cadherin_sliding_rebinding.py (Rakshit 2012 sliding-rebinding oracle -- effective_k_off / mean_lifetime; already runtime-imported, no change); spheroid/cadherin_bonds.py (the N_cad=223 scale bridge + contact_zone=10nm engagement-length convention this module reuses); cell/manifest.py + cell/cell.py build path (NEEDS a new two-cell assembler build_cell_doublet() -- the single biggest new piece; Cell.build is single-shell); cortex/cortical_tension.py (NEEDS the 'cadherin_'/'junc_actin' denylist-prefix edit -- PI-gated; the junction_tension observable also REUSES its _method_of_planes_gamma +…; cell/compartment_registry.py (cadherin_junction CompartmentSpec exists at EXPERIMENTAL; status flip to LIVE is PI-gated post-GATE-J); configs/recipes/multicell_junction.yaml (already declares cadherin_junction as declare_pending; un-defer post-gate); Notion SoT + outputs/tag_kb (SourceEvidence registration for Iturri-2020 + refresh); junctional_actin (junction/junctional_actin.py) -- a DOWNSTREAM dependent (requires=cadherin_junction) but stays STUB; NOT needed for GATE-J
- **Ready when:** GATE-J passes: an equilibrated suspended MCF7 doublet (built from the physiological baseline -- cytoplasm 65.9 Pa.s, turgor 133 Pa, nucleus, membrane_surface ON in BOTH cells) reports a mean engaged-cadherin_trans tension in the Sim-2015 1-5 pN band, with all four controls clean (OFF/ON differencing OFF==0; rigid-limit and fully-engaged-limit parities within 5%; cortex gamma bit-unchanged by the junction = no contamination; per-dimer lifetimes match the Rakshit oracle; zero COM drift). PLUS the Iturri-2020 SourceEvidence row exists in Notion with verdict OK, the 'cadherin_' denylist prefix is in cortical_tension + cross-checked by the registry test, and PI has signed off the cadherin_junction EXPERIMENTAL->LIVE status flip. THEN junctional_actin (the cortex-coupling belt) becomes the next multicell unit.
- **Effort:** MEDIUM-HIGH (~3-5 focused sessions). The mechanism (resolver, snapshot extender, catch-bond binder, attach helper) is ALREADY BUILT and statically tested -- the module is 'mostly runnable'. The real work is the NEW two-cell assembler (build_cell_doublet: two cortex shells + tag-offsetting + interface-cap cadherin seeding + box sizing + equilibration prelude) which Cell.build does not provide (~1.5-2 sessions); the junction_tension observable reusing the cortical method-of-planes machinery (~0.5 session); the denylist edit + registry cross-check test (~0.25 session, PI-gated); the GATE-J integration run + equilibration tuning + auto-viz (~1 session); and the Iturri SourceEvidence registration + refresh (~0.25 session, Lead/PI registry action). The k_trans/r_bind/k_on PI anchors are flagged for ratify-or-reject but do NOT block the gate (derived defaults are physically self-consistent: predicted 3.54 pN sits inside the band).

### PI parameters to ratify

**`k_trans`** — `DERIVED` / conf **MEDIUM** — *N/m*
- Proposed: 2.92e-3 (2.92 pN/nm)
- Citation: f0/contact_zone_width = 29.2e-12 N / 1.0e-8 m; f0 from Rakshit 2012 PNAS 109(46):18815 SI Table S1; contact_zone = 10 nm engagement length (spheroid scale bridge / test _CONTACT_ZONE).
- Notes: DERIVED, grid-invariant: per-dimer elastic force reaches the catch transition f0 at one contact-zone extension (same force-reaches-f0-over-one-zone construction as spheroid/cadherin_bonds.py). NOT tuned to a gate. PI choice is derivation-vs-direct-measurement: NO clean MEASURED single-molecule axial ectodomain (EC1-5) stiffness exists in the repo or in the searched literature -- the Manibog 2016 AFM '332 pN/nm' is the CANTILEVER/instrument spring, NOT the cadherin's own elasticity, so it does NOT anchor this. If…

**`r_bind`** — `ORDER_ESTIMATE` / conf **LOW** — *m*
- Proposed: 1.0e-8 (= contact_zone_width, 10 nm)
- Citation: EC1 trans-interface strand-swap reach ~ a few nm; defaulted to the contact_zone engagement length reused from the spheroid scale bridge.
- Notes: REQUIRED-with-documented-default capture radius for the cKDTree partner search. No tight single literature value at the x40 meso scale; the EC1 strand-swap dimer interface separation is ~1-4 nm (crystal structures, e.g. Boggon 2002 Science) so 10 nm is a generous-but-order-correct capture window. Governs binding CAPTURE only, never the bound-state stretch. Override via cfg['r_bind'] if a measured EC1 capture radius lands. Flag to PI as an order-estimate, not an invented number.

**`r0_trans`** — `ORDER_ESTIMATE` / conf **LOW** — *m*
- Proposed: 1.0e-8 (= contact_zone_width, 10 nm)
- Citation: Force-free ectodomain bridge span across the interface; defaulted to contact_zone_width.
- Notes: Trans-dimer rest length. Two opposed EC1-5 ectodomains span ~2x ~20 nm = ~40 nm crystal contour, but the relevant force-free SPAN across a settled interface is the engagement length, so 10 nm is the consistent default. Override via cfg['r0_trans']. Same provenance class as r_bind.

**`n_cad_per_cell`** — `SOLID_LITERATURE` / conf **HIGH** — *count*
- Proposed: 223
- Citation: N_cad = F_detach/f0 = 6.5e-9 N / 29.2e-12 N; Iturri et al. 2020 Cells 9(4):935 DOI 10.3390/cells9040935 (de-adhesion) over Rakshit 2012 f0.
- Notes: Particle COUNT only (grid-invariant; does not scale with n_fil). Fully-engaged ensemble n_cad*f0 = 6.51 nN reproduces the Iturri 6.5 nN de-adhesion datum EXACTLY (scale bridge closes). The per-dimer mechanics run off the fully-registered Rakshit catch-bond, so this anchor injects no mechanism, only a count + the ensemble-force report. BLOCKER: the Iturri anchor itself is UNREGISTERED in the Notion SoT -- see activation_gate.contamination_check and wiring_steps; must get a verdict-OK SourceEvidence row before it is…

**`k_on`** — `DERIVED` / conf **MEDIUM** — *1/s*
- Proposed: ~ effective_k_off(0) (rest-symmetric reformation)
- Citation: = 1/mean_lifetime(0, RAKSHIT_W2A); same rest-symmetric reformation convention as spheroid/cadherin_bonds.py.
- Notes: Trans-dimer on-rate defaulted to the rest off-rate so the rest bound fraction phi = k_on/(k_on+k_off(0)) = 0.50 -> ~112 engaged dimers at rest. Override via cfg['k_on'] if a measured cadherin on-rate (e.g. Rakshit k+1/k+2 are the INTRA-dimer rebinding rates, not the trans-capture on-rate) is preferred. Flag: this is a modelling convention, not a measured trans on-rate -- candidate for a PI literature anchor (cadherin association rate ~1e4-1e5 /M/s -> concentration-dependent; the rest-symmetric default sidesteps…

**`junction_tension_band`** — `SOLID_LITERATURE` / conf **HIGH** — *pN (per-dimer, homeostatic)*
- Proposed: 1-5 pN
- Citation: Sim S, et al. 2015 Mol Biol Cell e14-12-1618 (MDCK, E-cadherin FRET tension sensor), PMC4571300; band recorded in docs/PI_EXP_VALIDATION_MAP.md, confidence HIGH.
- Notes: The activation-gate OBSERVABLE band for per-dimer homeostatic junction tension. Model prediction at the operating point: thermal per-dimer force = sqrt(kT*k_trans) = 3.54 pN -> inside band. This is the falsifiable target; written BEFORE the run, not loosened. (Ensemble-level cross-check: Iturri 6.5 nN de-adhesion via n_cad*f0; cluster apparent surface tension 0.15-0.75 mN/m, Pallares 2023 PMC7617391, is a SECONDARY aggregate-scale check, not the single-doublet gate.)

### Activation gate

- **Observable:** GATE-J (junction tension). PRIMARY: per-engaged-dimer homeostatic trans-tension of a settled, equilibrated suspended cell DOUBLET, measured as <F_dimer> = mean over engaged cadherin_trans bonds of k_trans*max(0, L-r0_trans). SECONDARY (ensemble cross-check): the junction holding force = sum of engaged-dimer tensions, and a junction_tension surface-density gamma_junc computed by the SAME method-of-planes / Irving-Kirkwood estimator the cortex uses but restricted to cadherin_trans bonds across the interface plane.
- **Literature band:** PRIMARY band: per-dimer homeostatic junction tension 1-5 pN (Sim 2015 MBoC e14-12-1618, MDCK E-cadherin FRET, HIGH confidence). Model self-consistent prediction at the derived operating point: 3.54 pN (= sqrt(kT*k_trans)), inside band. SECONDARY band: fully-engaged ensemble holding force ~6.5 nN (Iturri 2020), with rest bound-fraction phi~0.5 giving ~3.3 nN at rest -- both order-consistent.
- **Controls:**
  - OFF/ON differencing: run the IDENTICAL two-cell box with junction enabled=false (no cadherin_trans bonds) vs enabled=true; the junction_tension observable must be exactly 0 in OFF and >0 in ON (the differencing isolates the junction load path, mirrors the…
  - Rigid/limit parity: (a) FROZEN-positions limit -- hold the two cells rigid at a fixed interface gap and confirm <F_dimer> = k_trans*(gap-r0) analytically (static bond-force parity, no kinetics). (b) Fully-engaged limit -- force phi->1 (k_on>>k_off) and…
  - No-contamination check: assert the cortical-gamma estimator over the SAME doublet is UNCHANGED by the presence of cadherin_trans bonds, i.e. every cadherin_* bond is filtered out by cortex/cortical_tension.cortical_bond_typeid_mask (denylist). cortex…
  - Detailed-balance / single-pair parity: per-dimer sampled lifetime distribution at a fixed clamped force F must reproduce validation.cadherin_sliding_rebinding.mean_lifetime(F) within Monte-Carlo error (the binder samples the SAME oracle the spheroid module…
  - Momentum control: total linear momentum injected by the junction = 0 (Newton-3rd-law harmonic pair); the doublet center-of-mass must not drift under the junction force alone.
- **PASS / PARTIAL / REFUTE:** PASS: equilibrated doublet's mean engaged-dimer tension in [1, 5] pN AND OFF/ON differencing clean (OFF==0) AND no-contamination (cortex gamma unchanged) AND rigid-limit & fully-engaged-limit parities within 5%. PARTIAL: per-dimer tension within a factor of 2 of the band (0.5-10 pN) with all controls passing -- indicates the operating point (k_on/k_off bound fraction or k_trans derivation) needs a PI parameter anchor but the mechanism is sound.…
- **No-γ-contamination check:** TWO contamination axes. (1) GAMMA contamination (physics): GAMMA_DENYLIST_PREFIX='cadherin_' must be a SUBSET of the cortical-tension denylist. ACTION REQUIRED -- cortex/cortical_tension.py currently denylists only ADHESION_BOND_TYPES + 'fa_actin_clutch' prefix; 'cadherin_' (and 'junc_actin') are NOT yet in ADHESION_BOND_TYPE_PREFIXES. Add them (PI-gated denylist edit) and extend the registry…

### Wiring steps

1. 1. TWO-CELL BUILD (the genuinely new construction). Cell.build is single-cortex-shell. Author a new assembler build_cell_doublet() (or extend manifest with a 'scope: multicell' / 'cells: [A, B]' block) that: (a) builds two cortex-baseline shells (each via the existing resolve_baseline -> build_cortex_full_simulation machinery, at the physiological setpoint: cytoplasm 65.9 Pa.s, turgor 133 Pa,…
2. 2. SURFACE SEEDING. Seed n_cad_per_cell (=223) cadherin tips on EACH cell's interface-facing cap (the spherical-cap region pointing at the partner cell), NOT the whole sphere -- only the contact-facing cadherins can form trans-dimers. Reuse the south-cap seeding idiom from _extend_snapshot_with_fa in cell.py (anchor + capture_radius), reflected to the interface direction. Anchor each cadherin to…
3. 3. CONFIG WIRING (5-point plug-in, mirror an existing optional like fa/lamellipodium). (a) manifest: a 'cadherin_junction' block under a new multicell scope (the multicell_junction.yaml recipe already DECLARES it as declare_pending). (b) resolve_cadherin_junction() already exists -- call it from a new branch in manifest.resolve_baseline (or the doublet resolver), passing kT, dt,…
4. 4. REGISTER cadherin_trans in the bond types at build time. attach_cadherin_junction requires 'cadherin_trans' already in snap.bonds.types (it sets params but the type must exist). The cleanest wiring (per the module's own note) is for the doublet build to register cadherin_trans on the SINGLE host Harmonic; the zero-fill path in attach_cadherin_junction is the safe fallback.
5. 5. FLIP STATUS EXPERIMENTAL->LIVE in compartment_registry.py for cadherin_junction ONLY AFTER GATE-J passes and PI signs off (status flip + recipe wiring is a gate-contract change, PI-authored). junctional_actin stays STUB (its alpha-catenin/vinculin coupling constants are unanchored -- separate PI decision; do NOT wire it for the bare junction-tension gate).
6. 6. DENYLIST EDIT (PI-gated): add 'cadherin_' and 'junc_actin' to cortex/cortical_tension.ADHESION_BOND_TYPE_PREFIXES and extend the registry gamma-denylist cross-check test (see contamination_check axis 1).
7. 7. SourceEvidence registration for Iturri (Lead/PI registry action; see contamination_check axis 2).

### Test plan

- STATIC (extend tests/test_cadherin_junction.py, already covers OFF-identity/dimensional/sign): add a two-cell snapshot fixture (two small cortex shells) and assert extend_snapshot_with_cadherins adds exactly 2*n_cad_per_cell cadherin particles with disjoint cell_a_tags/cell_b_tags, and the cadherin particle type registered once.
- PARTNER-SEARCH UNIT: drive CadherinTransJunctionUpdater on a hand-placed 2-cell snapshot where exactly K A-cadherins sit within r_bind of K B-cadherins; assert <=K trans bonds form, all are A<->B (no intra-cell), and one-dimer-per-cadherin holds (bound count <= min free cadherins).
- RIGID-LIMIT PARITY: freeze a doublet at a fixed interface gap, run the binder, and assert mean cadherin_trans force = k_trans*max(0,gap-r0_trans) analytically (no kinetics) -- the GATE-J rigid control.
- DENYLIST CROSS-CHECK: assert cortex/cortical_tension cortical_bond_typeid_mask returns False for 'cadherin_trans' (after the prefix edit), and that compute_cortical_gamma on a doublet WITH junction bonds == the value WITHOUT them (no-contamination control).
- KINETICS PARITY: clamp a single trans-dimer at force F over many batch ticks; assert the empirical mean lifetime reproduces validation.cadherin_sliding_rebinding.mean_lifetime(F) within MC error (the catch peak near f0=29.2 pN must appear; reuses the existing catch-bond oracle).
- BATCH-CFL: assert resolve_cadherin_junction caps batch_steps to satisfy batch_steps*dt*k_off_max <= 1e-3 (at dt=1.3e-8 the cap is well below the requested 100), and that the adaptive sub-stepper in act() keeps any stretched-dimer break Bernoulli within budget.
- INTEGRATION / GATE-J RUN: build a suspended doublet at the physiological baseline, equilibrate (soft-start, as FA needs -- the cell isn't born relaxed-on-interface), then measure the junction_tension observable. PASS iff mean engaged-dimer tension in [1,5] pN AND all controls clean. Auto-generate the figure (per the production-driver…

### GPU / CFL

Per-step cost is just the native HOOMD md.bond.Harmonic kernel (GPU-resident, no Python per-step force) -- hot-path priority P2. The CadherinTransJunctionUpdater is a per-BATCH Action using sim.state.get_snapshot() (rank-0 gather, like the D2 xlink binder), NOT a per-step cpu_local_snapshot; it runs on the host (CPU) each batch_steps with a scipy cKDTree partner search O(n_cad log n_cad) ~ negligible vs cortex. GPU-main debt: cross-cell binder host-sync + a future cupy/native partner search -- the optimisation lever, NOT required for the single-interface H.7 doublet (LOW bottleneck risk; risk rises only when MANY interfaces are…


---

## Junctional actin coupling — `junction/junctional_actin.py`  *(Phase E)*

**Goal.** Anchor the four PI-pending constant groups (k_couple + k_anchor; α-catenin catch set k_catch0/x_catch/k_slip0/x_slip; k_on; max_couple_dist + anchor_r0) with literature/derivation, wire junctional_actin into the manifest→Cell.build assembler, and gate it on a "junction transmits cadherin load into cortical tension without contaminating γ" activation contract.

- **Dependencies:** cadherin_junction (aleph/junction/cadherin.py) — HARD prerequisite: owns the cadherin particles and supplies n_cad_per_cell (Iturri 2020 N_cad≈223 / KU-4.17 ~100 mature) that…; Cortex baseline (cortex/cortex.py, crosslinkers.py, myosin.py) — supplies the same-cell cortex bead population the belt couples to, the kT (4.28e-21 J), the dt, and the γ_cortex…; cortex/cortical_tension.py — the γ estimator that MUST denylist junc_actin; the no-contamination control depends on its bond-prefix denylist being live.; Physiological baseline stack (cytoplasm 65.9 Pa·s, enclosed_volume turgor 133 Pa, membrane_surface, nucleus) — all ON at setpoint for both ON and OFF runs (HARD rule).; PI sign-off — k_catch0/k_slip0 are ORDER_ESTIMATE not in the KB; the parallel-vs-sequential form caveat must be ratified before the build path is un-gated. The…
- **Ready when:** PI ratifies (a) the two ORDER_ESTIMATE catch rates k_catch0=1.0/k_slip0=0.02 (or sanctions the H.3 filamin transfer) WITH the parallel-surrogate-vs-Buckley-sequential depth caveat explicit, and (b) the DERIVED transfers for k_couple/k_anchor/k_on/max_couple_dist/anchor_r0; AND the upstream cadherin_junction compartment has cleared its own junction-tension gate so n_cad is real; AND the module's enabled+anchored build path (currently the reserved NotImplementedError at junctional_actin.py:693) is implemented; AND the two-cell activation gate (Δγ_junction > 0, catch signature F*∈5-10 pN, engaged fraction 0.91, zero γ contamination, rigid/decoupled limit parity) returns PASS at the physiological baseline. Until all four hold, the compartment stays default-OFF and the enabled path correctly raises NotImplementedError.
- **Effort:** Medium, ~3-5 focused days once cadherin_junction is gated. Breakdown: ~0.5 day PI-anchor ratification + Magic-Number Block authoring (constants are mostly transfers/KB-anchored, so the research is largely done here); ~1.5 days implementing the reserved enabled+anchored build path (extend_snapshot head placement + bond registration + per-r0-binned same-cell cortex cKDTree seeding) and wiring the JunctionalActinCouplingUpdater onto hoomd.custom.Action; ~0.5 day manifest/CellBuildOptions/registry wiring + dependency guard; ~1 day the two-cell activation-gate smoke (build doublet, apply junction load, Δγ measurement + no-contamination + rigid/decoupled parity) + figures; ~0.5 day test expansion. CRITICAL PATH BLOCKER: cadherin_junction must be wired+gated first — junctional_actin is the third leg and is meaningless without the cadherin leg carrying load. The constant-anchoring itself is low-effort because x_catch/x_slip are already KB-anchored (Buckley 2014/KB-4.17) and the rest are sanctioned H.3 transfers; the real work is the build-path implementation + the gate.

### PI parameters to ratify

**`x_catch (catch-pathway force distance)`** — `SOLID_LITERATURE` / conf **HIGH** — *m*
- Proposed: 4.0e-9
- Citation: Buckley et al. 2014 Science 346:1254211, Fig.4 catch arm; recorded in KB-4.17 as dx_catch=4 nm (KB query verbatim: 'catch arm distance 4 nm, Buckley 2014 Fig4')
- Notes: KB-4.17 holds this as the Buckley α-catenin/F-actin catch-arm Bell distance, status=verified/High. This is the single most defensible anchor in the whole set — it is THE molecular bond this compartment represents (α-catenin→F-actin), measured single-molecule, and already in the KB. KB flag 'Phase 2 only, NOT slip distance' means do not reuse it as a slip distance; here it correctly seeds the CATCH pathway.

**`x_slip (slip-pathway force distance)`** — `SOLID_LITERATURE` / conf **HIGH** — *m*
- Proposed: 0.4e-9
- Citation: Buckley et al. 2014 Science 346:1254211, Fig.4 slip arm; KB-4.17 dx_slip=0.4 nm ('Buckley 2014 Fig4, slip arm')
- Notes: KB-4.17 anchored slip-arm distance for the same bond. Pairs with x_catch=4 nm to give x_catch > x_slip (the catch-bond ordering the catch_off_rate law and Sanity-Gate sign test require). Both distances are from the SAME paper/figure for the SAME bond — strongest possible internal consistency.

**`k_catch0 (catch-pathway zero-force off-rate)`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *1/s*
- Proposed: 1.0
- Citation: Derived to match Buckley 2014 Science Fig.4 zero-force lifetime τ(0)~1 s for the cadherin-catenin/F-actin bond (order-of-magnitude); functional form = Pereverzev 2005 PNAS 102:11281 two-pathway parallel sum (same family as H.3 filamin)
- Notes: NOT in KB (KB query explicitly: 'peak force, τ(F), k_off at zero force are absent from KB-4.17 — surface to PI'). Proposed so k_off(0)=k_catch0+k_slip0=1.02/s gives τ(0)~0.98 s, consistent with Buckley's ~1 s zero-force lifetime. With x_catch=4nm/x_slip=0.4nm this yields catch peak F*=6.0 pN. HONEST FORM CAVEAT (already in PI_DECISIONS): the parallel-sum surrogate exaggerates catch DEPTH — it predicts τ_max~26 s at F*, whereas Buckley's sequential two-state fit peaks at only ~1-2 s near 8-10 pN. The SHAPE…

**`k_slip0 (slip-pathway zero-force prefactor)`** — `ORDER_ESTIMATE` / conf **MEDIUM** — *1/s*
- Proposed: 0.02
- Citation: Pereverzev 2005 PNAS 102:11281 two-pathway structure (k_slip0 < k_catch0 so catch dominates at low F); magnitude set to place catch peak F* in the physiological 5-10 pN band per Buckley 2014 Fig.4
- Notes: Chosen so k_slip0 << k_catch0 (catch branch dominates below F*) — the same authoring logic as the H.3 filamin set (filamin_k_slip0=0.02, filamin_k_catch0=0.1 in phase1_h3.yaml). With k_catch0=1.0 gives F*=6.0 pN and a slip envelope k_off~0.33/s at 30 pN (well inside the batch CFL). This is the constant most exposed to the form/source caveat; PI to ratify.

**`k_couple (head↔cortex coupling-clutch stiffness)`** — `DERIVED` / conf **MEDIUM** — *N/m*
- Proposed: 1.0e-6
- Citation: Sanctioned transfer from H.3 dynamic crosslinker attach bond k_attach=1.0e-7 N/m (phase1_h3.yaml, KU-3.19 0.1 pN/μm), scaled up ~10× to reflect the reinforced vinculin-doubled α-catenin arm (le Duc 2010 JCB 189:1107; Yonemura 2010 NCB 12:533)
- Notes: No direct single-molecule α-catenin–actin stiffness in the corpus. Two defensible anchors: (a) the H.3 attach-bond stiffness 1.0e-7 N/m (a like-for-like dynamic actin-attach clutch — most conservative, choose this if PI wants pure transfer); (b) ~1e-6 N/m to encode the vinculin-reinforced second arm. Grid-invariant (per-bond, no n_fil). Bond-relaxation CFL with γ_cortex=6π·65.9·30nm=3.73e-5 N·s/m: at k_couple=1e-6, τ=γ/k=37 s, dt_max=0.1·τ=3.7 s — trivially satisfied at host dt~1e-8 s. Even k_couple=1e-5 leaves…

**`k_anchor (head↔cadherin permanent anchor stiffness)`** — `DERIVED` / conf **MEDIUM** — *N/m*
- Proposed: 1.0e-5
- Citation: Transfer from H.3 intra-crosslinker stiffness k_intra=1.0e-7 N/m (constitutive head-to-head link, KU-3.19), stiffened ~100× because the α-catenin N-terminus↔β-catenin↔cadherin-tail link is a folded constitutive (force-independent) bond, not a dynamic clutch
- Notes: This is the constitutive α-catenin/β-catenin/cadherin-tail link — it should be STIFFER than the dynamic clutch so the head tracks its cadherin parent. Make it ≥10× k_couple. Conservative alternative = reuse k_couple value (a single coupling stiffness throughout). Grid-invariant. PI to pick the ratio; the physics is insensitive as long as k_anchor ≥ k_couple (head stays pinned to its cadherin).

**`anchor_r0 (head↔cadherin rest length)`** — `DERIVED` / conf **LOW** — *m*
- Proposed: 5.0e-9
- Citation: α-catenin/β-catenin complex spans ~5-10 nm (Yonemura 2010 NCB; Pokutta-Weis structural; order estimate); set short and force-free at construction (head placed at the cadherin position along the cadherin→cortex line)
- Notes: Geometric, not a force constant — the head is BORN at the cadherin tail so the anchor is force-free at t=0 regardless of exact r0. Any value in 0-20 nm is acceptable; 5 nm reflects the catenin complex footprint. Low impact; the build places the head force-free.

**`max_couple_dist (cortex-acceptor search radius)`** — `DERIVED` / conf **MEDIUM** — *m*
- Proposed: 6.0e-8
- Citation: Sanctioned transfer from H.3 crosslinker acceptor search radius max_bind_dist=60.0e-9 m (phase1_h3.yaml, brief §Crosslinkers) — the same broad-phase cKDTree reach used for the cortex-internal attach clutch
- Notes: The junctional belt binds to nearby same-cell cortex beads exactly as the crosslinker binds to nearby actin — reuse the proven 60 nm reach. Grid note: if the mesoscale cortex inter-bead spacing exceeds 60 nm the radius must rise to ≥ one bead spacing (use √(A/n) mesoscale_reach from cortex/connected_mesh.py, ~the same family as the cadherin contact_zone_width). Resolver already ValueErrors on ≤0. Seeds the per-r0 bins over (0, max_couple_dist].

**`k_on (single-head binding rate)`** — `DERIVED` / conf **MEDIUM** — *1/s*
- Proposed: 10.0
- Citation: Sanctioned transfer from H.3 crosslinker k_on=10.0 1/s (phase1_h3.yaml); vinculin-recruitment ON step is qualitative (Yonemura 2010 NCB; le Duc 2010 JCB) so the binding rate is transferred, not measured
- Notes: Gives equilibrium engaged fraction k_on/(k_on+k_catch0+k_slip0)=10/(10+1.02)=0.91 at F≈0 — the Sanity-Gate measurement-protocol target (line 139-143). Matches the H.3 near-fully-bound-at-zero-load convention (filamin gives 0.988). Grid-invariant per-head rate. CFL: batch CFL is the binding gate, see gpu_and_cfl.

### Activation gate

- **Observable:** Δγ_junction = cortical tension γ measured on a cortex bead population WITH a tensioned cadherin junction coupled through the junc_actin belt, MINUS γ on the identical cell with the belt OFF (cadherin junction present but uncoupled). Equivalently: the fraction of a known applied junction load F_junc that arrives in the cortex as measurable hoop tension via the belt. Secondary observable: engaged-coupling fraction at F≈0.
- **Literature band:** (1) Load-transmission: a mature adherens junction transmits a substantial fraction of cadherin-borne tension into the cortex — α-catenin is a documented tension transducer (Yonemura 2010 NCB 12:533; le Duc 2010 JCB 189:1107). With N_cad≈100-223 engaged heads (KU-4.17 / Iturri 2020 Cells 9:935) each carrying up to the ~6 pN catch-peak, the belt should transmit O(0.6-1.3 nN) into the cortex — a…
- **Controls:**
  - OFF/ON differencing: belt OFF (enabled=False, bit-for-bit no-op) vs belt ON — Δγ_junction is the ON−OFF difference, never the absolute ON number (mirrors active-γ = γ(myosin ON) − γ(myosin OFF)).
  - Rigid/limit parity: in the k_couple→∞ rigid limit the cortex bead and junction head move as one — Δγ must saturate at the full cadherin load (no load lost in a compliant clutch); in the k_couple→0 limit Δγ→0 (decoupled). Both limits asserted.
  - No-γ-contamination check (HARD Rule 5): assert every junc_actin_* bond is on the cortical-tension estimator denylist (GAMMA_DENYLIST_PREFIX='junc_actin'); run the cortical-γ estimator with the belt ON and confirm the belt bonds contribute ZERO to the…
  - Physiological baseline (HARD): the ON and OFF runs both keep cytoplasm η=65.9 Pa·s, turgor=133 Pa, membrane_surface, nucleus at setpoint; the belt is measured FROM that baseline, not from a floppy shell.
  - Catch-vs-slip control: re-run with the catch pathway disabled (k_catch0=0, pure slip) — the junction lifetime-under-load and Δγ must DROP, isolating the catch reinforcement as the cause.
- **PASS / PARTIAL / REFUTE:** PASS: Δγ_junction > 0 and rises monotonically with applied junction load; engaged fraction at F≈0 within 0.91±0.05; catch signature F* in 5-10 pN; rigid limit saturates at full cadherin load; zero γ contamination. PARTIAL: Δγ_junction > 0 with correct catch SIGN but magnitude below the O(nN)/O(0.6-1.3 nN) band by <10× (consistent with the documented parallel-surrogate depth caveat) — report as shape-faithful-magnitude-bounded, NOT a fail;…
- **No-γ-contamination check:** Static test: junc_actin_bond_type_names() ⊆ cortical-tension denylist (prefix 'junc_actin'); CompartmentRegistry.gamma_denylist() includes 'junc_actin'. Runtime test: with belt ON + junction loaded, the cortex/cortical_tension.py estimator's bond iterator skips every junc_actin_* bond — assert the cortex-internal γ is IDENTICAL (bitwise within fp tol) to the belt-OFF γ when myosin is held fixed,…

### Wiring steps

1. 1. configs/mcf7_baseline.yaml (or the multicell_junction recipe override): add an optional_subsystems.junctional_actin block with enabled:true and the seven anchored constants (k_couple, k_anchor, anchor_r0, max_couple_dist, k_catch0, x_catch, k_slip0, x_slip, k_on) under a Magic-Number Block comment citing each source — mirroring the phase1_h3.yaml filamin block layout. Keep n_bins:10,…
2. 2. cell/manifest.py resolve_baseline(): import resolve_junctional_actin; in the optional_subsystems section add `ja_b = opt.get('junctional_actin'); if _enabled(ja_b): ...`. CRITICAL ordering — junctional_actin REQUIRES cadherin_junction, so resolve the cadherin block FIRST and pass its resolved n_cad as n_cadherin into resolve_junctional_actin(cfg, kT=p_cortex.kT, dt=dtc,…
3. 3. cell/manifest.py ResolvedBaseline dataclass + the return block (manifest.py:265-281): add p_junctional_actin field and p_cadherin field; thread both through.
4. 4. cell/cell.py CellBuildOptions: add a `junctional_actin: bool = False` flag (mirroring the fa/substrate/lamellipodium flags). Add the resolved-param plumbing to Cell.build.
5. 5. cell/cell.py: add _extend_snapshot_with_junctional_actin(snapshot, p_ja, cadherin_tags=..., cortex_positions=..., cortex_tag_start=...) call that delegates to junctional_actin.extend_snapshot_with_junctional_actin — but FIRST the module's enabled+anchored build path (currently the reserved NotImplementedError at line 693) must be implemented: append one junc_actin head per cadherin tag at the…
6. 6. cell/cell.py attach_*: register the JunctionalActinCouplingUpdater (once it subclasses hoomd.custom.Action) on the simulation operations at the D2 batch cadence, after the cadherin CadherinBondUpdater and the crosslinker XlinkBondUpdater (same updater family).
7. 7. cell/compartment_registry.py junctional_actin CompartmentSpec: flip status STUB→EXPERIMENTAL once anchored+gate-passed; fill n_particles_mesoscale='~100-223 belt heads per interface (Iturri N_cad, n_fil-independent)' and n_bonds='≤2·n_head (1 anchor + ≤1 couple per head)'. Keep denylist_bond_types=('junc_actin',), requires=('cadherin_junction',). Do NOT edit the registry until the gate passes…
8. 8. Promote in configs/recipes/multicell_junction.yaml: move cadherin_junction + junctional_actin from declare_pending → enable once both pairwise gates clear.

### Test plan

- Extend aleph/tests/test_junctional_actin.py (already exists): keep the OFF-identity bit-for-bit no-op test; add an ANCHORED-resolve test (all seven constants present → is_anchored True → build no longer raises).
- Dimensional/scalar law tests (analytic, no HOOMD): coupling_force_magnitude(p, Δr) ≥ 0 and = k_couple·Δr for Δr>0, = 0 for Δr≤0 (tensile-only sign-sense). catch_off_rate(p,F): assert biphasic — d k_off/dF < 0 at F=0 (catch), d k_off/dF > 0 above F*, minimum at the analytic F*=(kT/(x_catch+x_slip))·ln((k_catch0·x_catch)/(k_slip0·x_slip))…
- Denylist static test (HARD Rule 5): junc_actin_bond_type_names(n_bins) — every name startswith 'junc_actin'; CompartmentRegistry.gamma_denylist() ⊇ {'junc_actin'}.
- Engaged-fraction test: at F≈0 the two-state steady state k_on/(k_on+k_catch0+k_slip0)=0.91 ± fp.
- CFL gate tests: batch_steps·dt·k_off_max ≤ 1e-3 with k_off_max = catch envelope (k_catch0 + k_slip0·exp(F_credible·x_slip/kT)); bond-relaxation dt ≤ 0.1·γ_cortex/k_couple (non-binding, assert it passes at host dt).
- Newton-3rd-law conservation: after a (built) extend, net momentum injected by junc_actin_anchor + junc_actin_couple bonds = 0 (symmetric Harmonic pairs).
- Two-cell smoke (the activation gate): build a cortex-baseline doublet at the physiological baseline, cadherin_junction ON, junctional_actin OFF then ON; apply a known junction load; measure Δγ_junction via cortex/cortical_tension.py with the junc_actin denylist active; assert PASS/PARTIAL/REFUTE thresholds + the no-contamination runtime…
- Figure: extend aleph/scripts/h1_h2_vis.py (or the H.7 vis entry-point) to emit catch_off_rate(F) with F* marked, the engaged-fraction-vs-load curve, and the Δγ ON−OFF bar with the literature band overlaid (per the visualize-at-closeout rule).

### GPU / CFL

Two CFLs, both stated in the module Sanity Gate. (1) Batch off-rate CFL (shared D2 contract): batch_steps·dt·k_off_max ≤ 1e-3. k_off_max is the catch envelope ≈ k_catch0 + k_slip0·exp(F·x_slip/kT); at the proposed set and F up to 30 pN, k_off_max ≈ 0.33/s, so batch_dt ≤ 1e-3/0.33 = 3.0 ms — with host dt~1e-8 s this allows batch_steps up to ~3e5; the default batch_steps=100 is far inside. (2) Bond-relaxation CFL: dt ≤ cfl_safety·γ_cortex/k_couple. γ_cortex = 6π·η·R_bead = 6π·65.9 Pa·s·30e-9 m = 3.73e-5 N·s/m (physiological cytoplasm, NOT water). At k_couple=1e-6 N/m: τ=γ/k=37 s, 0.1·τ=3.7 s ≫ host dt — non-binding; even k_couple=1e-5 leaves…
