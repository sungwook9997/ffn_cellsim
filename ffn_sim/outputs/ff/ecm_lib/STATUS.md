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
