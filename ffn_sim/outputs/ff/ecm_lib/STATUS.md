# FF ECM Library — session status (2026-07-10)

**Goal (PI):** Build a comprehensive, high-precision ECM library on the FF (Cytosim-physics,
Warp) engine: multiple ECM materials (individual + mixed), 2D and 3D, several alignment degrees,
at 5R×5R lateral extent (R=7.5µm → 75×75µm). Verify by *pressing* (indentation) that each gives
the **real literature Pa**. Produce validation figures for every case. Register missing materials
in the KB (create where absent).

## Unit foundation (locked)
FF engine = µm · pN · s. **1 pN/µm² = 1 Pa exactly** → all moduli come out directly in Pa.
κ[pN·µm²]=k_BT·L_p ; EA[pN]=E_fibril[Pa]·π·r_f²[µm²] ; k_seg[pN/µm]=EA/seg.

## Design

### Two constitutive classes (this distinction IS the physics)
- **Fibrillar ECM** (collagen-I, fibrin, agarose): modulus **EMERGES** from a Mikado semiflexible
  fiber network (density → connectivity ⟨z⟩, bending κ, stretch EA, crosslinks). We set the
  microstructure from literature and *validate* the emergent Pa — NO tuning to the modulus.
- **Continuum gel** (polyacrylamide, hyaluronic acid, Matrigel): a flexible chemical/BM gel with
  molecular (~nm) mesh unresolvable at cell scale. E is a **material INPUT** (literature, set by
  chemistry). Represented as a coarse elastic spring lattice whose spring stiffness is the analytic
  inverse of the target E (constitutive discretization, not a fit). Indentation then *validates* that
  pressing returns the input E (harness + discretization consistency), cross-checked against the
  independent affine shear/uniaxial modulus.

### Modules
- `ffn_sim/ff/ecm_library.py` — `ECMSpec` grounded registry (6 materials), nematic-S alignment
  sampler (2D von-Mises / 3D Watson), 2D & 3D Mikado builder, concentration→density, composite builder.
- `ffn_sim/ff/ecm_mechanics.py` — modulus harness, all → Pa: (1) affine **shear** G, (2) affine
  **uniaxial** E (+ axial/transverse anisotropy for aligned nets), (3) spherical **indentation**
  (Hertz inversion) E_eff — the "pressing" test.
- `ffn_sim/scripts/ff_ecm_build_all.py` — build the full type×alignment×dim matrix at 5R×5R.
- `ffn_sim/scripts/ff_ecm_validate.py` — measure moduli, compare to literature bands, tables+figs.
- `ffn_sim/scripts/ff_ecm_viewer.py` — standalone ECM HTML viewer (fibers+crosslinks+indenter).

### Literature validation anchors (from KB)
- collagen-I G'(c): **1→5, 3→55, 7→342 Pa** (KB-1.V.2.1, Yang-Kaufman); G0≈36 Pa @1.5 mg/mL,
  band 30–100 Pa (KB-1.30, Licup2015). Lp=17µm, r_f=50nm, E_fibril=1.1MPa, ⟨z⟩~3.2, ξ~1/√c (KB-1.1/1.2/1.3/1.7).
- PA/PAA gel: **E 0.1–50 kPa** tunable, default 5 kPa, ν≈0.45 (KB-1.21/1.5, Discher2005/Palchesko2012).
- alignment S: healthy stroma S<0.1, tumor 0.3–0.7, TACS-1/2/3 (KB-1.9/1.V.2.4); tendon/muscle S>0.8.
- fibrin / Matrigel / agarose / HA: **not in KB** — being researched (DOI-verified) for registration.

## Progress log
- 2026-07-10: explored FF ECM code. KB queried (collagen+PAA present). Research workflow DONE
  (PA/fibrin/Matrigel/agarose/HA cards, DOI-verified). Wrote `ecm_library.py` + `ecm_mechanics.py`.

### Validation results (2026-07-10, CPU REV)
- **Harness VALIDATED on known-modulus continuum (PA gel E=5000 Pa):** uniaxial E=5000 (calib),
  shear G=1562 (=0.91× E/[2(1+ν)]=1724), **indentation E_eff=3970 Pa (0.79×), flat E-vs-δ (Hertzian).**
  → the 3 measurement routes are mutually consistent and the indentation inversion is correct.
- **Collagen-I bulk modulus at reference conc MATCHES literature:** G(1.0mg/mL)=9 Pa, G(1.5)=10 Pa
  vs KB 5 / 11 Pa; ⟨z⟩=3.09 (KB 3.2) ✓. Energy route (G=2U/Vγ²), finite-EA spring, converges by 5k steps.
- **Collagen c-scaling: G~c^1.0 (measured) vs c^2.0 (KB/Yang-Kaufman).** Documented finding, NOT a bug:
  linear athermal Mikado is density-limited; the steep experimental exponent needs the thermal/nonlinear
  semiflexible physics (Licup 2015). Reference-conc absolute Pa is correct — the primary "real Pa". No tuning.
- **Fibrillar indentation reads the LOCAL modulus (softer than bulk shear)** — real non-affine local
  softening (KB-1.14), not a bug (harness validated on continuum). Needs bead ≫ mesh ξ to approach bulk.

### Key method decisions
- Energy route primary; boundary-reaction kept as cross-check (unreliable on sparse nets).
- axial_mode='spring' (finite EA) default; 'reshape' diverges with affine pre-strain → deprecated here.
- Local Warp is CPU-only (CUDA not in this build); native/large runs → gbook A5000.

## Checkpoint 2026-07-10 (30-min loop)
Core goal COMPLETE + committed (fcd208f): 6/6 materials IN band, alignment/2D/3D/composite, native 5R×5R
GPU, novelty (strain-stiffening + decoupling), 6 figs + 2 viewers (browser-verified), Notion Dev-Log posted.
**Next (this cycle): register the 5 new material SourceEvidence + KnowledgeClaim rows to the Notion
Contract-Graph** (user authorized "없으면 신설하여 진행") — the one stated goal item still only staged.

## Checkpoint 2 (2026-07-10, later loop cycles)
- KB registration DONE (commit 25639e8): KB-1.V.4.1-4.5 in Notion SoT + duckdb (167 claims), status=verified.
- **Extension: durotaxis stiffness-GRADIENT substrate** (new capability). `build_gradient_ecm` in
  ecm_library.py (spatially-graded continuum bond stiffness) + `indentation_modulus(center_xy=...)` probe +
  `ff_ecm_gradient.py` (physiological 1 / pathological 10 / sharp ~100 Pa/µm regimes vs KB-1.V.1.3, probe
  E(x) by indenting along the gradient axis). Running on gbook A5000.
- Note: cron fired ~8× stacked during long turns — handled as ONE check-in (no repeated work).

## Checkpoint 3 (2026-07-10)
- **Extension: emergent viscoelastic stress relaxation** (KB-1.6). `stress_relaxation` in ecm_mechanics.py
  (step shear + crosslink Bell-slip turnover → G(t); bond-virial σ_xz) + `ff_ecm_viscoelastic.py`. Result:
  **τ ∝ 1/k_off EMERGES** (τ·k_off≈0.82 const) — α-actinin ~12s → weak-physical ~200s (KB-1.6 30-1000s ✓)
  → covalent-LOX ~1170s; G∞/G₀≈0.14 (86% relaxation). figs/stress_relaxation.png (eye-verified).
  Extensions this session: native GPU · strain-stiffening · mode-decoupling · durotaxis gradient · viscoelastic.

## Checkpoint 4 (2026-07-10) — QA pass
Core goal + 5 extensions COMPLETE. This cycle: adversarial correctness review of ecm_library.py +
ecm_mechanics.py (units/formulas/signs/numerics — the "high precision" mandate), focused on what the
empirical validations don't cover (energy-route factor-of-2, bond virial, Hertz inversion, continuum
calibration, alignment sampler, re-pin logic, NaN/zero paths). Will fix any confirmed correctness bugs.

## Checkpoint 5 (2026-07-10) — c-scaling root-caused
Tested the c^2 hypothesis after the QA stretch-dominated insight. Result: exponent is a mechanical-regime +
rigidity-percolation effect, BRACKETED: spring(stretch-dominated) n=1.07 ↔ reshape(bending-dominated,
crosses rigidity threshold) n=7.24; real c²≈2 needs a PI-decided ⟨z⟩(c) crosslink law kept above the
sub-isostatic threshold (my rule gives ⟨z⟩=2.18=floppy at 1mg/mL). No tuning. → CSCALING_REGIME_FINDING.md.
This is a modeling decision (collagen crosslink biology) — PI-scoped, not auto-tuned.

## Checkpoint 6 (2026-07-10) — KB registration FULLY completed
Created 21 SourceEvidence rows (DOI-verified) + wired Claim→Evidence relations for all 5 material claims
(KB-1.V.4.1-4.5: 4/4/5/4/4 evidence links). source_evidence 355→376; DOIs now citation-auditable (were
text-only). SoT snapshot refreshed; kb-check gates green. "TAG KB 신설" is now complete end-to-end
(claim + evidence + relations + duckdb mirror).

## Checkpoint 7 (2026-07-10) — physical connectivity control
Added target_z crosslink control (build_fibrillar_ecm(target_z=...)): subsamples fiber-pair crosslinks to
the KB-1.3 literature ⟨z⟩~3.2 instead of letting ⟨z⟩ blow up with density (was 2.18→14.2). Anchored to an
existing datum, NOT tuned to n. At fixed physical ⟨z⟩=3.2 → clean c^1.04 (density-linear, the honest fixed-
connectivity prediction); reference G(1.5mg/mL)=13.4 Pa still matches. Default None (backward-compat).
The c-scaling investigation is now cleanly closed: physical ⟨z⟩ controllable; c²≈2 needs a PI-decided ⟨z⟩(c)
GROWTH law (collagen LOX biology), not an n-fit.

## Checkpoint 8 (2026-07-10) — PI feedback: PA by concentration + indentation visualized
PI asked: (1) is PA gel done per-concentration? (2) how is the pressing measured — visualize it.
- **PA by recipe/concentration**: PA_FORMULATIONS table (8 recipes, acrylamide%/bis% → E, KB-1.V.4.1 Subramani
  + Engler tissue ladder) + pa_formulation_E(); ff_ecm_pa_ladder.py builds each and PRESSES it — indentation
  E_eff reproduces the recipe Pa across **150 Pa–40 kPa** (uniaxial exact on 1:1; indent ~0.6-0.7× continuum
  factor). figs/pa_stiffness_ladder.png.
- **Indentation visualized**: ff_ecm_indent_viz.py — progressive spherical indentation, animated 3D HTML
  (cross-section + full-3D scenes, bead descending into gel) + figs/indent_curve.png (the DIMPLE surface
  profile deepening + the measured F(δ) with Hertz fit F=(4/3)E*√R·δ^1.5). All eye-verified in browser.

## Checkpoint 9 (2026-07-10) — fibrin per-concentration validation
Applied the PI's per-concentration principle to the other fibrillar ECM. ff_ecm_fibrin_conc.py: fibrin
G'(c) sweep 0.5-8 mg/mL — G 17.7→247.6 Pa, ALL in the Piechocka band (0.1-2000 Pa) ✓; exponent n=0.96
(vs lit 2.3) = the same athermal stretch-dominated regime as collagen (CSCALING_REGIME_FINDING.md applies to
both). figs/fibrin_concentration.png. Fibrillar per-concentration validation now symmetric (collagen + fibrin).

## Checkpoint 10 (2026-07-10) — ROADMAP + executing NEAR #1 (stiffness-sensing)
PI asked for a development plan. Wrote ROADMAP.md (5-lens adversarial workflow → 3-horizon synthesis;
ROADMAP_workflow.json). Executing the NEAR top pick: ff_stiffness_sensing.py — a resting full-compartment
cell on a Winkler compliant substrate (k_sub=2Ea/(1−ν²)), stiffness swept across the library's material
moduli (150 Pa brain-PA → 40 kPa muscle-PA). CPU smoke (COARSE nf=800, non-authoritative): physiological
baseline (ΔP=40.5 Pa, γ=0.153, V/V0=0.998); traction rises 0→0.078 nN with a sharp engagement threshold
(bound 0.01@150Pa → 0.80@2kPa) = durotaxis basis. E=input, k_sub=derived, optimum=measured (no tuning).
Biphasic PEAK needs native (denser myosin → higher per-clutch load shifts the optimum into range) — native
--nf 38000 A5000 run LAUNCHED. Path (a) (grip the actual library ECM network) = documented next increment.

## Checkpoint 11 (2026-07-10) — ROADMAP NEAR #3: full virial Cauchy stress tensor
ecm_material_stress(ecm,pos)→σ[3,3]: crosslink+segment central-bond virial + bending atomic virial (full
tensor generalizing the xz functions). shear_modulus now returns G_virial. Method-independence PROVEN:
PA continuum G_energy=G_virial=G_react=1522 Pa (all agree); collagen sparse G_energy=16.1 ≈ G_virial=16.4
(0.98×, agree) while G_reaction=1.4 (11× low = the unreliable route). So the virial (grid-invariant truth)
CONFIRMS the energy route is accurate and identifies reaction as the outlier — resolving the two-readout
inconsistency. Unblocks the nonlinear master curve + N1 + r^-n (all are components of this one tensor).

## Checkpoint 13 (2026-07-11) — README consolidated + NEAR #4 launched
- README brought to full current state (cell-ECM stiffness-sensing, virial tensor, N1, all scripts, KB
  REGISTERED). Commit 3dfb195.
- **NEAR #4 (6-material native modulus atlas)** ff_ecm_native_atlas.py: builds each of the 6 materials at a
  native-density 40µm REV on the A5000, measures the modulus (fibrillar→shear G via virial tensor;
  continuum→calibrate+indentation), checks vs literature band + REV↔native consistency. ⚠️first launch crashed
  KeyError G_virial_Pa = STALE gbook ecm_mechanics.py (virial tensor added locally in NEAR #3/#5 wasn't re-
  synced); re-synced + relaunched (virial-tensor presence verified on gbook). Monitor armed for completion.
- Session so far executed ROADMAP NEAR #1(native-confirmed)/#3/#4(running)/#5. Lesson: re-sync ff/ after every
  local kernel/API change before a gbook run.

## Checkpoint 14 (2026-07-11) — NEAR #4 native atlas: 6/6 IN BAND at native scale
ff_ecm_native_atlas.py on A5000: all 6 materials at a 40µm native-density REV — collagen 15.2 (33.6k nodes),
fibrin 78.7 (74.7k), PA 3654, HA 331, Matrigel 439, agarose 16117 Pa — ALL in the literature band, and each
MATCHES its small-REV value (modulus is intensive → REV↔native consistency confirmed). Fibrillar measured via
the NEAR-#3 virial tensor (G_virial), integrating #3+#4. figs/native_atlas.png. NEAR #4 modulus production
done; per-material native full-extent VIEWERS = remaining half. Roadmap NEAR #1/#3/#4/#5 done this session.

## Checkpoint 15 (2026-07-11) — NEAR #2 Path (a): cell grips the LIVE library ECM network
ff_stiffness_sensing.py `--mode ecm-network`: the resting cell adheres to an actual ecm_library network
(build_ecm PA-gel at E, two-sided clutch_ecm_spring_kernel + ECM explicit substeps), so stiffness is the
ECM's OWN calibrated modulus (deepest library integration; Path b Winkler kept intact). ECM_CFL_SAFETY=0.01
(diagnosed stability margin, physics-invariant). CPU coarse smoke (nf=800, NON-AUTHORITATIVE):
- ⭐CLEAN signal — ECM DEFORMATION tracks stiffness: soft PA(150Pa) 552nm → stiff(40kPa) 6nm (the cell feels
  ECM stiffness through how much it can deform it = the mechanistic durotaxis basis; Path (a) coupling VALIDATED).
- traction curve is ERRATIC/non-monotonic at coarse (2.57→2.54→0.09→0.01→0.87 nN; the 40kPa jump = a small-N
  stochastic artifact, only 6/200 clutches bound loaded to 144pN) — NOT a clean biphasic; needs native
  (more clutches → statistical averaging) + ensemble. Path (a) [clutches loaded ABOVE F*≈7pN, slip side] and
  Path (b) [Winkler a=0.05µm keeps forces BELOW F*, catch side] bracket the Bangasser-Odde biphasic.
- Path (b) confirmed byte-intact (150→0, 8000→0.162). Native Path (a) launched for the authoritative curve.

## Checkpoint 16 (2026-07-11) — NEAR #6: contact-guidance traction anisotropy on ALIGNED collagen (Path a)
ff_contact_guidance_anisotropy.py (new driver) + fa_ecm.traction_director_anisotropy/nematic_order_2d +
sense_ecm_network(alignment_S=, director=) threading. The resting cell grips aligned collagen-I (nematic S from
the Watson sampler, director fixed in-plane) via the SAME isotropic two-sided clutch; anisotropy must EMERGE
from the aligned microstructure (no directional force prefactor — unfitted, S is the swept IV). Design +
KB-grounding via a 4-lens workflow (validation bands: Ray 2017 >3× force anisotropy, Szulczewski 2021 up-to-35×
matrix directional stiffness, Niraula 2025 traction-saturation, Riching 2014 persistence).
Two readouts — a RATIO cancels the few-clutch stochasticity that made the absolute Path-a traction noisy:
- PRIMARY A_F=F∥/F⊥ (engaged clutches; the Ray force anisotropy). SECONDARY R_σ=σ∥/σ⊥ (ECM virial Cauchy
  tensor, clutch-count-INDEPENDENT → authoritative if native engaged counts stay low).
Readout unit-tested (pure-∥→F⊥=0,S_bound=1; isotropic-4dir→A_F=1.0 exact; nematic S(∥)=1/S(⟂)=−1).
COARSE CPU smoke (nf=500, 1 seed, S={0,0.83}, NON-AUTHORITATIVE) — plumbing VALIDATED + a clear diagnosis:
- ✅ aligned collagen builds (S_measured 0.09→0.82 tracks target), sense loop + both readouts run.
- ✅ R_σ (grid-invariant) is the ROBUST directional signal — monotone 0.78→4.26 with alignment even at coarse
  (independent of the 71→56 engaged-clutch count) — exactly the spec's predicted authoritative signal. NOTE
  R_σ is the cell-scale directional STRESS response (Szulczewski regime), smaller in magnitude than the pure-
  shear E∥/E⊥=63 (a different quantity, validated separately by the library).
- ⚠️ A_F (traction) is NOISY at coarse (1.97→0.93, non-monotone; A_F=1.97 at isotropic S=0 IS the single-seed
  few-clutch geometric noise) — needs NATIVE + ≥3-seed ensemble to resolve above the ~60-clutch layout noise
  (exactly as the spec predicted). Whether the cell TRACTION anisotropy A_F follows the matrix STRESS anisotropy
  R_σ is the open question the native run answers.
Native gbook run queued (gbook SSH was transiently down). MT-OFF (same as NEAR #1/#2 sensing family) — open PI item.

## Checkpoint 17 (2026-07-11 night) — NEAR #6 KB-grounding (the KB already covers contact guidance)
KB query (duckdb) for the #6 validation bands: the primary source + the claim to validate against ALREADY exist —
- SourceEvidence `Ray2017_NatCommun` (DOI 10.1038/ncomms14923, force anisotropy >3×) — **source_audit verdict OK**
  (citation-verified → citable in the deliverable). Also present + OK: `Provenzano2006_BMCMed`, `Conklin2011_AmJPathol`,
  `Bredfeldt2014_JBiomedOpt`, `Ray2021_COCellBiol`.
- KnowledgeClaims to cite: **KB-2.14 "Contact guidance via anisotropic FA elongation"** (the existing claim my
  EMERGENT A_F(S) validates), KB-1.V.2.4 "Fiber alignment TACS-1/2/3", KB-1.9 / KB-3.20 (nematic order parameter).
So #6 CITES existing kb_ids (no new SourceEvidence for Ray/Provenzano). NEW (not in KB, PI-gated DRAFT only):
Szulczewski 2021 (35× matrix directional stiffness), Riching 2014 (persistence), Niraula 2025 (traction saturation).
kb-check gate at this commit: verify_runs OK (31, no drift) · verify_params OK (42, no drift). No DRIFT.

## Checkpoint 18 (2026-07-12 morning) — BOTH native jobs LANDED (overnight pipeline completed)
The gbook A5000 chain ran #2 (3.5h) → #6 (7.5h) to completion; GPU now free. Both figures visually verified.
### Native #2 — Path (a) pa_gel stiffness sensing (Nc=266000, 3 E, --tag native_ecmnet)
traction 0.177nN@150Pa → 0.007@2kPa → 0.009@40kPa (bound 0.33→0.01→0.02); ECM-disp 275→4→0 nm.
→ Path (a) traction DESCENDS with stiffness (soft ECM engages+deforms+traction more) = the slip/DESCENDING arm;
with Path (b) Winkler ASCENDING (0→0.090nN, bound 0.02→0.80) they BRACKET the Bangasser-Odde biphasic at
f/clutch≈F*≈7pN. ECM-deformation is the clean monotone signal. `figs/stiffness_sensing_ecmnet_native.png`,
`ecm_stiffness_ecmnet_native.json`. NEAR #2 ✅ DONE native-confirmed.
### Native #6 — contact guidance on aligned collagen (Nc=266000, 4 S × 3 seeds)
| S | A_F=F∥/F⊥ | R_σ=σ∥/σ⊥ | n_eng |
|---|---|---|---|
| 0.00 | 0.91±0.13 | 1.6±0.8 | 97 |
| 0.30 | 0.80±0.06 | 4.5±4.3 | 90 |
| 0.59 | 0.83±0.15 | 6.4±4.3 | 89 |
| 0.83 | 0.67±0.04 | 30.8±23.0 | 75 |
→ ⭐ The grid-invariant **matrix stress anisotropy R_σ EMERGES** with alignment (1.6→30.8×, tracks the library
E∥/E⊥ ladder 1/3.1/10.3/63, ≈ Szulczewski ≤35×) — the matrix directionally feels the alignment. The resting
cell's **clutch TRACTION A_F stays ~isotropic** (0.91→0.67, no alignment trend; per-clutch ≪ F*, weak engagement).
HONEST split: passive quasi-static sensing reads matrix directional STIFFNESS but NOT active traction guidance
(Ray-2017 >3× needs polarized protrusion/contraction — the FF motility layer). Figure KEY corrected to match the
data (was over-claiming "traction 2-4×"; re-plotted from JSON). KB: Ray2017 (verdict OK) + KB-2.14 contact-guidance
+ KB-1.V.2.4 TACS cited (existing). NEAR #6 ✅ DONE native; active-motility follow-up → MID.
### ALL 6 NEAR items now native-confirmed. Next: MID + the #6 active-motility follow-up. kb-check: runs/params OK.

## Checkpoint 19 (2026-07-12) — NEAR #6 interactive HTML viewer (per PI viz directive)
ff_contact_guidance_viewer.py → contact_guidance_viewer.html (4.4 MB, full-res, no downsample): the aligned
collagen-I microstructure across the S ladder (0→0.83) as a multi-scene interactive viewer, each scene labelled
with measured S + library E∥/E⊥ + the native R_σ (1.6/4.5/6.4/30.8×); a green line marks the fixed director. The
morphology counterpart to the comparison figure. Real-browser verified (browser_check.py: rendered + exited, no
JS errors; screenshot eyeballed — isotropic S=0 fibers all-directions + director axis correct). Fibrin
generalization run still in flight on the A5000 (~4/12 at check; A_F~1, R_σ~1 at S=0 so far). Output/viz locations
kept as-is per PI (ffn_sim/outputs/ff/ecm_lib/ + figs/).

## Checkpoint 20 (2026-07-12) — Job D native (3-seed): directional stress propagation (MID, KB-1.10)
ff_ecm_stress_propagation.py native (box=90µm, 382725 nodes, 100k relax steps, 3 seeds; GPU-resident relax, no
host-CG → 786s for 12 runs). Eshelby-type contractile inclusion (ε=20%, frozen — no force tuning), spherical
far-field BC, σ_rr(r) via ecm_mechanics.stress_field_radial (shell virial).
| S | n_exp (|σ_rr|~r⁻ⁿ) |
|---|---|
| 0.00 | 8.12±2.57 |
| 0.30 | 7.96±0.96 |
| 0.59 | 5.69±0.79 |
| 0.83 | 4.72±0.23 |
→ ⭐ n DROPS MONOTONICALLY with alignment (8.1→4.7): aligned collagen CHANNELS contractile stress ~10× farther
(the figure shows purple S=0.83 sustaining ~10× more σ at r~28µm than blue S=0). 3-seed ensemble resolved the
single-seed S=0.30 outlier; error tightens as S rises (±0.23 at S=0.83). Absolute n≈5-8 STEEPER than elastic
(r⁻²⁻³) = sub-isostatic athermal Mikado localizes stress (short-range) — same stretch-dominated limit as
c-scaling/N1; true fibrous long-range (Notbohm r⁻¹) needs PI-gated nonlinear/bending physics. Equilibration
lesson: the big native network needed 100k steps (4k gave n=15 under-relaxed). Figure verified. TO DO: stress-
field interactive HTML (displacement-colored, iso vs aligned). Compares to KB-1.10 as oracle-overlay (not registered).

## Checkpoint 21 (2026-07-12) — Job D stress-field interactive HTML viewer (per PI viz directive)
ff_ecm_stress_propagation_viewer.py → stress_propagation_viewer.html: the displacement field a contractile
inclusion drives through collagen-I (isotropic vs aligned scenes), still matrix faint, the pulled top-20% coloured
blue→red, green director. The morphology behind the r⁻ⁿ channeling. Real-browser verified (renders, no JS errors,
eyeballed: pulled region concentrated near the inclusion along the director — displacement is LOCALIZED, the
short-range n≈5-8 result made visible). Representative 40µm REV (mechanism view; the quantitative n(S) is the
native box=90 figure). Both viz types now exist for Job D (comparison figure + interactive HTML), per PI.

## Checkpoint 22 (2026-07-12) — Job D generalizes to FIBRIN (aligned-channeling is a general fibrous property)
ff_ecm_stress_propagation.py --material fibrin, native (box=90µm, 637875 nodes — fibrin denser than collagen,
100k relax, 3 seeds). n(S): S=0→16.83±3.18, S=0.30→17.94±5.04, S=0.59→8.61±0.96, S=0.83→8.29±0.57.
→ ⭐ Fibrin shows the SAME trend as collagen — n DROPS with alignment (16.8→8.3, a bimodal split at S~0.5):
aligned fibrin channels contractile stress farther. So the aligned-stress-channeling is a GENERAL fibrous-network
property, not collagen-specific. HONEST caveat: fibrin's ABSOLUTE n is higher than collagen's (4.7-8.1) because
fibrin is denser (637k vs 382k nodes) AND the bigger network is less equilibrated at 100k steps (the collagen
lesson: n falls as steps rise) — so fibrin's absolute n is an UPPER bound (range a lower bound); the TREND is the
robust, converged result. Fixed a plot-title bug (was hardcoded "collagen-I"; now uses meta['material']).
`figs/stress_propagation_sp_fibrin_native.png`. Fibrin contact-guidance (R_σ) run still in flight (~half done).

## Checkpoint 23 (2026-07-12) — consolidated native validation dashboard (per PI "완성·비교 피규어")
ff_ecm_summary_dashboard.py → figs/ecm_validation_dashboard.png: one 4-panel figure tying the whole ECM-library
program to its literature/KB anchors (reads committed native JSONs, no re-sim). (1) 6/6 materials vs literature
bands (real Pa). (2) alignment→R_σ(S) tracks the library E∥/E⊥ ladder, approaches Szulczewski ≤35×. (3) stress-
propagation n(S) collagen+fibrin both drop with alignment, both above the KB-1.10 continuum(3)/fibrous(1) refs
(FF steeper = honest sub-isostatic limit). (4) stiffness sensing Path b(Winkler ascending)+Path a(live-ECM
descending) bracket the motor-clutch biphasic. Verified. The "완성" review artifact for the PI's morning.

## ═══ OVERNIGHT SESSION SUMMARY (2026-07-11 night → 07-12) — for PI morning review ═══
The ECM-library program is COMPLETE + native-A5000-confirmed + fully visualized. What landed overnight (all
committed, all figures/viewers real-browser-verified, kb-check gates green, NO tuning):
- **NEAR (6/6 DONE native):** stiffness sensing Path a (live library ECM, descending/slip) + Path b (Winkler,
  ascending/catch) BRACKET the motor-clutch biphasic; contact guidance on aligned collagen — R_σ=σ∥/σ⊥ EMERGES
  (1.6→30.8× tracking library E∥/E⊥) while resting-cell traction A_F stays ~isotropic (passive sensing reads
  matrix STIFFNESS, active traction guidance needs motility). #6 interactive HTML viewer.
- **MID Job D (DONE native):** directional stress propagation from a contractile inclusion — n drops with
  alignment 8.1→4.7 (collagen), aligned channels stress ~10× farther; GENERALIZES to fibrin (16.8→8.3). KB-1.10
  anchored (FF n>3 steeper than continuum n=3 — sub-isostatic limit, honest). Stress-field HTML viewer.
- **c-scaling gap RESOLVED (diagnosis):** KB-anchored as THERMAL (KB-1.30 G0~kBT·ℓp²/ξ⁵ + KB-1.7 ξ~c^-1/2 →
  c^2.5), NOT ⟨z⟩(c) (KB-1.3 ⟨z⟩~3.2 c-independent — ruled out as over-fit). Analytic overlay figure.
- **Consolidated validation dashboard** (`figs/ecm_validation_dashboard.png`) — all results vs literature/KB.
- **In flight:** fibrin contact-guidance R_σ (native ~8/12) — early finding: fibrin R_σ WEAKER than collagen
  (S=0.59: fibrin ~1.7 vs collagen 6.4) → contact-guidance R_σ is a PARTIAL generalization (unlike Job D's full
  one), because the resting cell loads the denser fibrin network less directionally. Land when done (~2.7h).
- **BLOCKED/PI-GATED (deferred, documented):** #6 active-motility follow-up (native-crawl grid-drag blocker);
  thermal-WLC runtime physics (c² resolution); ν-faithful continuum; ECMRemodeler Bell-rate datum. Awaiting PI.
- **Note:** concurrent sessions (DCM T1-KMC, M6 contractility) share dcm/main; my commits are clean specific-adds,
  no cross-contamination. Some ecm_lib figs show as modified by those sessions — left untouched.

## Checkpoint 24 (2026-07-12) — fibrin contact-guidance LANDED: NEAR #6 GENERALIZES (both readouts)
Native fibrin (nf=38000, Nc=266000, 4 S × 3 seeds, 28197s). Final aggregates:
| S | R_σ=σ∥/σ⊥ | A_F=F∥/F⊥ |
|---|---|---|
| 0.00 | 0.9±0.4 | 1.14±0.17 |
| 0.30 | 1.2±0.6 | 0.98±0.12 |
| 0.59 | 4.9±4.6 | 0.89±0.13 |
| 0.83 | 22.1±22.9 | 0.80±0.06 |
→ ⭐ BOTH NEAR #6 findings GENERALIZE from collagen to fibrin: (a) the matrix stress anisotropy R_σ EMERGES with
alignment (fibrin 0.9→22.1 vs collagen 1.6→30.8 — same order, comparable magnitude); (b) the resting cell's clutch
traction A_F stays ~isotropic (fibrin 1.14→0.80, like collagen 0.91→0.67). So contact guidance (matrix feels
alignment, passive traction doesn't) is a GENERAL fibrous-network property, not collagen-specific.
⚠️ HIGH seed variance (R_σ ±22.9 at S=0.83, driven by seed3=54.5) — the matrix stress anisotropy under cell load
is stochastic (fiber realization near the cell); the MEAN trend is robust. LESSON (re-confirmed): the 2-seed
intermediate reading (S=0.59 → "fibrin weaker") was PREMATURE — the 3-seed mean shows full generalization. Always
wait for the full ensemble before concluding. Fixed the plot material-label bug (was hardcoded "collagen-I").
`figs/contact_guidance_anisotropy_cg_fibrin_native.png`. This completes the fibrin generalization program (both
Job D stress-propagation AND #6 contact-guidance now confirmed general across collagen + fibrin).

## Checkpoint 25 (2026-07-12) — mammary-stroma tissue-mimetic: the athermal ceiling at TISSUE scale (scope/limit)
ff_ecm_tissue_mimetic.py: interpenetrating collagen-I+Matrigel composites, literature-anchored (NOT tuned) vs
KB-1.V.1.2. Result (virial G, box 30µm REV): normal mammary stroma G=89 Pa (band 140-400, ~1.6× UNDER); tumor
stroma TACS-3 (dense 5 mg/mL aligned collagen) G=129 Pa (band 5000-10000, ~50× UNDER). HONEST scope finding (not
a success story): the athermal sub-isostatic fibrillar library reproduces individual DILUTE-material moduli (6/6
in band) but UNDER-shoots realistic DENSE/tumor TISSUE moduli — because desmoplastic stiffening relies on the
thermal/nonlinear collagen stiffening the athermal Mikado lacks (the c-scaling/N1/stress-propagation ceiling, now
at tissue scale). CRUCIALLY continuum gels DO reach kPa (agarose 14k, PA 40k in band), so the ceiling is SPECIFIC
to the fibrillar collagen route → the runtime fix is the PI-gated thermal-WLC. This BOUNDS the library scope:
trustworthy for individual materials + soft tissues + the mechanistic cell-ECM coupling; needs thermal-WLC for
dense tumor stroma. `figs/tissue_mimetic_mammary.png`. Consistent with the c-scaling diagnosis (checkpoint 22).

## Checkpoint 26 (2026-07-12) — DIRECTIONAL stress channeling: aligned collagen is a fiber WAVEGUIDE (novel)
ecm_mechanics.stress_field_directional + ff_ecm_stress_propagation directional readout (n∥ = decay ALONG the
director cone, n⊥ = decay ACROSS). Native collagen (box=90µm, 100k relax, 3 seeds):
| S | n∥ (along fibers) | n⊥ (across fibers) |
|---|---|---|
| 0.00 | 12.8 | 12.0 |
| 0.30 | 7.1 | 14.8 |
| 0.59 | 5.2 | 21.1 |
| 0.83 | 4.5 | 21.6 |
→ ⭐ the DIRECT channeling signature: at S=0 isotropic (n∥≈n⊥≈12); as alignment rises n∥ DROPS (12.8→4.5, stress
propagates FAR along fibers, toward the elastic n=3) while n⊥ RISES (12.0→21.6, stress dies fast across fibers).
Aligned collagen is a stress WAVEGUIDE — it funnels contractile force along the director and blocks it
perpendicular. This EXPLAINS the isotropic-average n(S) drop (8.1→4.7, checkpoint 20): the average falls because
the ∥ channel opens. The mechanistic heart of TACS-3 directional long-range force transmission (aligned tumor
collagen transmits contraction along invasion highways = mechanical cell-cell communication, KB-1.10). Unfitted,
no directional force prefactor — emerges from the aligned microstructure + isotropic bond virial. Coarse-smoke +
native-confirmed. `figs/stress_propagation_directional_sp_dir_native.png`. Novel Job D refinement (checkpoint 20).

## Checkpoint 27 (2026-07-12) — 2D stress-propagation: added --dim, but 2D directional is statistics-limited (honest)
Threaded `--dim {2,3}` through ff_ecm_stress_propagation (2D planar vs 3D bulk; build_fibrillar_ecm already
supports dim=2). Coarse 2D smoke (box=44) HONEST finding: a 2D fibrillar sheet is TOO SPARSE (~2000 nodes vs 44730
in 3D at the same box — nodes scale as area + fewer fiber crossings), so the directional sector-shell bins have too
few bonds → noisy/unphysical exponents (e.g. n⊥=−7.17, negative = insufficient statistics). A clean 2D directional
result would need a much larger box (≥200µm) — deferred as marginal vs the clean, comprehensive 3D result (the 3D
directional waveguide, checkpoint 26, is the authoritative channeling result). The --dim feature is committed for
future use. No 2D conclusion drawn (correctly — the smoke showed it's not ready).

## Checkpoint 28 (2026-07-12) — 2D directional channeling: the waveguide is DIMENSION-GENERAL (box=200, honest)
Following up checkpoint 27: box=200µm 2D gives enough nodes (42000, vs the too-sparse box=44 smoke) to test the
directional channeling in 2D. Native 2D collagen (dim=2, box=200, 100k relax, 3 seeds):
| S | iso-n | n∥ | n⊥ |
|---|---|---|---|
| 0.00 | 20.0 | 17.9 | 30.3 |
| 0.30 | 12.4 | 13.2 | 39.5 |
| 0.59 | 11.1 | 10.9 | 33.7 |
| 0.83 | 7.2 | 6.5 | 10.6 |
→ the channeling GENERALIZES to 2D: n∥ DROPS cleanly with alignment (17.9→6.5) and n∥<n⊥ at EVERY S (the fiber
waveguide holds in 2D). BUT honest caveat: n⊥ is NON-MONOTONE with huge error bars (30→39→34→10.6) — a 2D aligned
network has most fibers ∥ so the ⊥ sector is inherently bond-sparse → n⊥ unreliable; and box=200 is under-
equilibrated at 100k steps (absolute n high, like the 3D 4k-step lesson) — 2D's larger extent needs more relax
steps. So 2D CONFIRMS the waveguide qualitatively (via the clean n∥) but the 3D directional result (checkpoint 26,
n∥ 12.8→4.5 / n⊥ 12→21.6) remains the authoritative quantitative one. `figs/stress_propagation_directional_sp_2d_native.png`.
This completes "2D 3D 모두" for the flagship channeling (dimension-general, honestly bounded).

## ═══════ FINAL COMPLETE SUMMARY (2026-07-12, supersedes the partial one above) — PI REVIEW ═══════
Overnight autonomous ECM-library program, FULLY EXHAUSTED. All native A5000-confirmed (Nc=266000 for cell runs),
all figures + interactive HTML viewers real-browser-verified, kb-check gates green, ZERO tuning, memory updated.
Concurrent DCM(T1-KMC)+M6(contractility) sessions share dcm/main — my commits are clean specific-adds only.

**1. MATERIALS (from 07-10, standing):** 6/6 real-Pa in band (collagen 12, fibrin 63, PA 3.4k, HA 309, Matrigel
411, agarose 14k), 2D+3D, alignment S→E∥/E⊥ 1/3/10/63, composites. 5 novel (strain-stiffening, mode-decoupling,
viscoelastic τ∝1/koff, durotaxis gradient, virial tensor+N1). KB-1.V.4.1-4.5 registered.

**2. CELL-IN-ECM (07-12 native):**
- Stiffness sensing — Path b Winkler (ASCENDING/catch) + Path a live-ECM (DESCENDING/slip) BRACKET the Bangasser-
  Odde biphasic at f/clutch≈F*≈7pN. ECM-deform 275→0nm.
- Contact guidance — R_σ=σ∥/σ⊥ EMERGES with alignment (collagen 1.6→30.8, fibrin 0.9→22.1: GENERALIZES); resting-
  cell traction A_F stays isotropic (needs active motility, blocked). KB-2.14 + Ray2017(OK).

**3. Job D STRESS PROPAGATION (07-12 native, MID/KB-1.10):**
- |σ|~r⁻ⁿ, n DROPS with alignment (collagen 8.1→4.7, fibrin 16.8→8.3): aligned channels stress ~10× farther,
  GENERALIZES. FF n>3 steeper than continuum (sub-isostatic ceiling).
- ⭐DIRECTIONAL WAVEGUIDE — n∥ (along fibers) 12.8→4.5 vs n⊥ (across) 12→21.6: aligned collagen funnels stress
  along the director. DIMENSION-GENERAL (2D box=200 confirms via n∥; 2D n⊥ ⊥-sparsity-noisy). The mechanistic
  heart of TACS-3 directional force transmission.

**4. LIMITS honestly bounded:**
- c-scaling gap = THERMAL (KB-1.30 G0~kBT·ℓp²/ξ⁵ + KB-1.7 ξ~c⁻¹ᐟ² → c^2.5; KB-1.3 rules out ⟨z⟩(c)); runtime fix
  = PI-gated thermal-WLC.
- tissue-mimetic: athermal fibrillar UNDER-shoots dense/tumor tissue (normal 89 vs 140-400, tumor 129 vs 5-10k Pa);
  continuum gels DO reach kPa → the ceiling is the collagen route → thermal-WLC.

**5. DELIVERABLES:** ff_{stiffness_sensing, contact_guidance_anisotropy, contact_guidance_viewer, ecm_stress_
propagation(+_viewer), ecm_summary_dashboard, ecm_cscaling_thermal, ecm_tissue_mimetic}.py + ecm_mechanics.{stress_
field_radial, stress_field_directional} + fa_ecm.{traction_director_anisotropy, nematic_order_2d}. Dashboard +
3 interactive viewers. README/ROADMAP/STATUS/OVERNIGHT_PLAN + memory current.

**6. AWAITING PI (gated, all documented — pick a direction & I resume immediately):**
(a) thermal-WLC runtime physics — closes c² + tumor stiffness (biggest lever); (b) #6 active-motility follow-up —
needs the native-crawl grid-drag/Stokes unblocker; (c) ν-faithful continuum; (d) ⟨z⟩(c) has no literature datum.

## Checkpoint 29 (2026-07-12) — CAPSTONE: a BM-gel continuum SHORT-CIRCUITS the collagen stress waveguide (novel)
ff_ecm_stress_propagation --composite-with matrigel (tumor-stroma composite = aligned collagen + Matrigel BM gel,
interpenetrating). Native (box=90µm, 480061 nodes, 100k relax, 3 seeds):
| S | iso-n | n∥ | n⊥ |
|---|---|---|---|
| 0.00 | 2.76 | 2.7 | 2.8 |
| 0.30 | 2.75 | 2.7 | 2.7 |
| 0.59 | 2.75 | 2.7 | 2.7 |
| 0.83 | 2.74 | 2.7 | 2.7 |
→ ⭐⭐ DEFINITIVE + ultra-clean (±0.01): the composite propagates stress like a LINEAR-ELASTIC CONTINUUM (n≈2.75,
alignment-INDEPENDENT) and PERFECTLY ISOTROPIC (n∥≈n⊥ at every S) — the Matrigel continuum COMPLETELY SHORT-
CIRCUITS the collagen fiber waveguide. Compare pure collagen (checkpoint 26): iso-n 8.1→4.7, n∥ 12.8→4.5 / n⊥
12→21.6 (strong waveguide). BIOLOGICAL PREDICTION (novel, testable): contact-guided DIRECTIONAL long-range force
transmission (TACS-3 invasion highways, mechanical cell-cell communication) is a property of SPARSE FIBER-DOMINATED
collagen — it is LOST when the interstitium is filled with basement-membrane gel (or dense continuum ECM), because
the isotropic continuum provides a dominant non-directional stress path. This UNIFIES the tissue-mimetic finding
(checkpoint 25: composite modulus continuum-dominated) with the channeling finding — the BM gel dominates the
composite in BOTH modulus and stress-transmission directionality. `figs/stress_propagation_directional_sp_comp_native.png`.
Unfitted, physical. Ties tissue-mimetic + Job D + directional into one biological prediction. --composite-with feature added.

## Checkpoint 30 (2026-07-12) — the fiber-waveguide CROSSOVER: gel stiffness vs the collagen channel (novel, quantitative)
ff_ecm_stress_propagation --composite-E sweep (aligned collagen S=0.83 + gel of stiffness E), native (box=80µm,
337721 nodes, 2 seeds each). n∥/n⊥ vs gel stiffness E:
| E_gel (Pa) | n∥ | n⊥ | waveguide |
|---|---|---|---|
| 5 | 4.7 | 11.5 | STRONG (n⊥≫n∥) |
| 15 | 4.3 | 6.7 | weakening |
| 50 | 3.4 | 3.8 | nearly lost |
| 150 | 2.9 | 2.9 | LOST (isotropic) |
| 500 | 2.7 | 2.7 | LOST (continuum) |
→ ⭐ quantitative CROSSOVER completing the capstone (checkpoint 29): the fiber waveguide (n⊥≫n∥) SURVIVES only when
the interstitial gel is SOFTER than the collagen network itself (E_gel ≲ collagen bulk G≈12 Pa); any stiffer gel
progressively SHORT-CIRCUITS it, fully lost by E≈150 Pa (both → isotropic continuum n≈3). Since real basement-
membrane gels (Matrigel 411 Pa) and tissue gels are MUCH stiffer than dilute collagen, the waveguide is generally
destroyed in composite tissue — directional long-range force transmission requires SPARSE, gel-poor fibrous regions
(the biological threshold, quantified). `figs/waveguide_crossover.png` + waveguide_crossover.json. Unfitted (E is
the swept input, the crossover EMERGES at the collagen-modulus scale).

## Checkpoint 31 (2026-07-12) — THERMAL-WLC IMPLEMENTED + VALIDATED: the c-exponent is fixed (PI-authorized, default-off)
The PI-authorized big lever (from the "너가 알아서 결정" mandate). Full Marko-Siggia extensible-WLC option
`axial_mode='wlc'` for fibrillar segments (crosslinks stay linear), DEFAULT-OFF, PROVISIONAL. Built to the design
spec: ff/wlc.py (law, 6 unit tests), wlc_spring_kernel (Warp, bit-parity with host to 1e-9), seg_Lc from the
EMERGENT geometric mesh ξ~c^(-0.500 verified) [the crux — seg_rest would give c¹], and axial_mode threaded through
ecm_mechanics (_links_for/_wlc_block/_to_device/_force_pass/_cfl_dt/_elastic_energy with BASELINE SUBTRACTION/
ecm_material_stress WLC tension). All constants derived (kBT=U.KBT 310K, Lp/EA from spec, ξ emergent, x_max=0.99756
from f'_WLC=EA) — NO tuning.
VALIDATION (collagen G'(c), native REV, spring vs wlc):
| c (mg/mL) | G_spring | G_wlc (Pa) |
|---|---|---|
| 1 | 9.3 | 0.17 |
| 2 | 18.4 | 0.72 |
| 4 | 44.9 | 5.9 |
| 7 | 76.9 | 25.3 |
→ ⭐⭐⭐ the c-EXPONENT RISES from spring n=1.11 to **WLC n=2.63** — CONFIRMS the KB-1.30 thermal diagnosis
(G0~kBT·L_p²/ξ⁵~c^2.5; checkpoint 22). The athermal→thermal physics is the c-scaling fix, validated. Anti-pre-
stress gate PASSES (RMS node drift 0.0001µm ≪ seg_rest — the thermal slack sits the network at its physiological
baseline, not floppy/pre-stressed). ⚠️ HONEST CAVEAT (the spec predicted this): the WLC absolute modulus drops
BELOW band at low c (0.17-0.72 Pa vs 11-13) — replacing the stiff EA rod with the ~33 pN/µm entropic spring makes
the network bending+entropic-co-dominated, softer than the EA-dominated athermal one. So WLC fixes the SCALING
(physics) but the reference-concentration ABSOLUTE anchoring (bending κ / target_z / ξ prefactor) is a PI-scoped
re-calibration — NOT auto-tuned (KB-1.3 forbids ⟨z⟩(c); the spec forbids fitting). Spring mode unchanged (n=1.11).
Remaining validation (next): strain-stiffening K~σ^{3/2} + N1 sign flip. `ff/wlc.py`, `tests/ff/test_wlc.py`.

## Checkpoint 32 (2026-07-12) — thermal-WLC validation complete: strain-stiffening ✓, N1 partial (honest)
Remaining WLC validation (collagen ref, spring vs wlc): STRAIN-STIFFENING — differential modulus K/K0 rises to
1.17× (spring, weak) vs **4.46× (WLC)** at γ=0.2 → the WLC produces the strong semiflexible nonlinear stiffening
(K~σ^{3/2} regime), a second independent confirmation of the thermal physics. N1 (first normal-stress): spring
N1/σxz@20%=0.08 (stays +); WLC ratio rises to 0.20 but N1 STAYS POSITIVE (no flip to the Janmey negative). Honest:
the WLC softens the axial stretch but the entropic tension still carries stretch, so it does NOT fully reach the
bending-dominated negative-N1 regime — consistent with the stretch-dominated diagnosis (CSCALING). So WLC fixes
2 of 3 thermal signatures (c-exponent 1.11→2.63 ✓, strain-stiffening 1.17→4.46× ✓) but not the N1 sign.
`figs/wlc_validation.png` (spring vs WLC c-scaling + strain-stiffening). THE ONE PI-SCOPED ITEM: the WLC absolute
modulus at reference c is below band (entropic softer than the EA rod) — re-anchoring (bending κ / confinement-tube
prefactor / target_z) is a modeling decision, NOT auto-tuned (hard rule). Next: WLC tumor tissue-mimetic (does the
nonlinear wall reach kPa at high strain?).

## Checkpoint 33 (2026-07-13) — WLC E.5 tumor check: does the nonlinear wall reach kPa? (honest: no, reinforces re-anchoring)
Dense aligned collagen (c=7, S=0.59) tangent modulus K at high strain, spring vs wlc: spring K=84→97 Pa (γ 0.1→0.5),
WLC K=71→94 Pa — BOTH ~50-100× below the tumor band (5000-10000 Pa). The WLC's nonlinear enthalpic wall does NOT
reach tumor kPa at these densities/strains — the low-c absolute-modulus softness carries through to dense collagen.
So thermal-WLC fixes the SCALING (c-exponent + strain-stiffening) but not the ABSOLUTE tissue stiffness; reaching
the kPa tumor range needs the reference re-anchoring (the one PI-scoped item) OR higher crosslink density than the
KB-1.3 ⟨z⟩=3.2 physiological value (which must not be tuned). thermal-WLC is COMPLETE as a validated default-off
option: it confirms the c-scaling physics is thermal (the documented gap's root cause), with the absolute anchoring
left as a clean PI modeling decision. This closes the PI-authorized thermal-WLC thrust.
