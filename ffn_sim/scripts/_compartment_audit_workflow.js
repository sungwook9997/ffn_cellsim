export const meta = {
  name: 'compartment-physics-audit',
  description: 'Exhaustive 5-lens adversarial physics audit of 8 compartments + apply confirmed fixes',
  phases: [
    { title: 'Audit', detail: '5 independent physics lenses per module (dimensional, literature, sign/conservation, numerics/magic-number, API/OFF-identity)' },
    { title: 'Verify', detail: 'adversarial refutation of every finding (kill false positives)' },
    { title: 'Fix', detail: 'one owner-agent per module applies CONFIRMED findings + re-runs tests' },
  ],
}

const REPO = '/Users/sw1/ffn_cellsim'
const ENV = 'source ~/.zshrc; conda activate ffn_sim; cd ' + REPO

const GROUND_TRUTH = `
PINNED GROUND TRUTH (use these; flag any module value that disagrees):
- kT(310 K) = 4.28e-21 J. k_B = 1.38e-23 J/K.
- F-actin: persistence length 17 µm; Young's modulus ~1.3-2.6 GPa; cross-section ~3.2e-17 m²;
  single-filament EA ~ 4.3e-8 N; rupture strain a few %.
- Microtubule: EI = 2.2e-23 N·m² (Gittes 1993); Lp = EI/kT = 5.14 mm; E ~ 1.2 GPa; A ~ 190 nm².
- Intermediate filaments (vimentin/keratin): Lp ~ 0.3-1 µm; diameter 10 nm; small-strain E ~ MPa;
  STRONGLY nonlinear strain-stiffening, extends 2-3.5x before rupture (Kreplak 2005).
- E-cadherin: Rakshit 2012 X-dimer sliding-rebinding catch-bond, f0 ~ 29 pN catch peak.
- Discrete WLC bending (HOOMD md.angle.Harmonic, U=½k(θ-θ0)²): k_angle = EI/ℓ0 [J/rad²], θ0=π straight.
- Discrete axial spring (HOOMD md.bond.Harmonic, U=½k(r-r0)²): k = E·A/ℓ0 [N/m].
- Overdamped CFL: bond τ = γ/k, dt ≤ safety·τ; bending τ_bend = γ·ℓ0³/EI (project convention, manifest.py:134).
- Young-Laplace: ΔP = 2γ/R. van't Hoff: Π = R·T·c. Kedem-Katchalsky water flux: Jv = Lp·(ΔP_hyd - σ·ΔΠ).
- Series springs: N identical bonds of stiffness k in a chain → effective fiber stiffness k/N.
  A per-bond stiffness is NOT the whole-fiber stiffness (off by the bead count).
- Bell slip: k_off(F) = k0·exp(F·x_β/kT). Catch needs a force-WEAKENING term (decreasing then increasing).
- HOOMD 7.0.1: md.bond.Harmonic / md.angle.Harmonic are builtin GPU ForceComputes; md.force.Custom
  needs set_forces(timestep) writing cpu_local_force_arrays; updaters are hoomd.custom.Action with act().
`

const AUDIT_RULES = `
You are auditing ONE physics module of ffn_cellsim, a fine-grained mechanistic HOOMD-blue
single-cell simulator. ${ENV}. You are an ADVERSARIAL physics reviewer: your job is to find
REAL errors, not to bless the code. Read the module file IN FULL. Where a literature value is
cited, VERIFY it (grep the repo references/: ls ${REPO}/ffn_sim/references; the TAG KB:
python ${REPO}/ffn_sim/outputs/tag_kb/tag_query.py "<question>" ; and your own knowledge).

${GROUND_TRUTH}

For EVERY issue, return a finding with: lens, severity (HIGH = wrong physics / wrong sign /
order-of-magnitude-wrong constant / fabricated citation / would crash when enabled; MEDIUM =
conceptually unsound derivation, missing factor, inconsistent with platform convention, untested
enabled path; LOW = docstring/unit-comment error, style), a precise location (file:line), the
CLAIM (what is wrong), the EVIDENCE (the numbers/derivation/quote proving it), and a concrete
recommended_fix. If a constant is genuinely unknown and the module sets it None + raises, that is
CORRECT (not a finding). Do NOT invent literature values. If you find NOTHING for your lens, return
an empty findings array — do not manufacture issues.

Quantify wherever possible. "k seems low" is useless; "k=1e-2 N/m as a per-bond stiffness on a
24-bead 4µm fiber gives k_fiber=4.3e-4 N/m, needing 23µm (575% strain) to reach the cited 10 nN —
~500x too soft" is a finding. Plug in real numbers.
`

const MODULES = [
  {
    name: 'ventral_stress_fibers', module: 'ffn_sim/cell/stress_fibers.py', test: 'ffn_sim/tests/test_stress_fibers.py',
    cites: 'Kumar 2006 Biophys J 90:3762 (SF tension 10-30 nN); Tojkander 2012 JCS / Hotulainen-Lappalainen 2006 JCB (sarcomeric α-actinin 0.5-1 µm); Furuike 2001 / Ferrer 2008 (α-actinin); KU-3.19.',
    seed: 'LEAD PRELIMINARY (confirm/quantify/refute): default k_actin=1e-2 N/m is derived from the MEASURED SF tension (a VALIDATION TARGET used as a model INPUT) and used as a PER-BOND stiffness ignoring the series-spring factor (n_beads-1). Worked: 24-bead 4µm fiber → k_fiber=1e-2/23=4.3e-4 N/m → 10 nN needs 23µm (575% strain) → ~500x too soft vs EA_bundle/ℓ0~5 N/m. Also SF NMII reuses cortex_myosin_* bond types → would contaminate cortical γ if wired LIVE (needs sf_myosin_* prefix).',
  },
  {
    name: 'linc', module: 'ffn_sim/cell/linc.py', test: 'ffn_sim/tests/test_linc.py',
    cites: 'Crisp 2006 (50 nm NE gap); Autore 2013 (nesprin spectrin repeats); Arsenovic 2016 (nesprin resting tension ~2 pN FRET).',
    seed: 'LEAD PRELIMINARY: k_linc=None → enabled build raises (correct). capture_radius default = 3.0·gap (50nm); the 3x multiplier is geometric/overridable. Verify the harmonic bond is builtin and the diagnostic force sign is a consistent restoring convention. Check the nucleus↔cytoskeleton pairing is geometric and does not create spurious long bonds.',
  },
  {
    name: 'intermediate_filaments', module: 'ffn_sim/cell/intermediate_filaments.py', test: 'ffn_sim/tests/test_intermediate_filaments.py',
    cites: 'Mücke 2004 J Mol Biol (vimentin Lp~1µm); Kreplak 2005 J Mol Biol (extend 2-3.5x); Block 2018 PRL (strain-stiffening); Herrmann (10 nm diameter); Guo 2013.',
    seed: 'LEAD PRELIMINARY: k_bb = E_if·A_if/ℓ_seg = 6e6·π(5nm)²/0.5µm = 9.4e-4 N/m (units OK, magnitude OK). Nonlinear strain-stiffening (the DEFINING IF property) raises NotImplementedError (honest). CHECK: does the model build a bending angle at all, or only backbone+crosslink? With Lp~1µm and ℓ_seg~0.5µm (ℓ_seg<Lp) bending matters — omitting it under-stiffens. Is E_if=6 MPa the right small-strain modulus or is it a large-strain value?',
  },
  {
    name: 'microtubules', module: 'ffn_sim/cell/microtubules.py', test: 'ffn_sim/tests/test_microtubules.py',
    cites: 'Gittes 1993 JCB 120:923 (EI=2.2e-23 N·m², Lp=5.2 mm); Kis 2002 PRL / Pampaloni 2006 PNAS (E~1.2 GPa, A~190 nm²); Walker 1988 JCB 107:1437 (DI rates); Mitchison-Kirschner 1984.',
    seed: 'LEAD PRELIMINARY: appears solid — k_angle=EI/ℓ0 and k_backbone=Y_stretch/ℓ0 wired to the right types (verified not swapped at lines 894-899), WLC identity EI=kT·Lp guarded to 2x, θ0=π, τ_bend=γℓ0³/EI matches cortex. Y_stretch=2.3e-7 N is a flagged placeholder. CHECK the DI updater (Walker rates), the MTOC anchor bond, and whether the stiff-bending CFL is correctly computed and whether ANY test exercises the enabled bonded-force build.',
  },
  {
    name: 'osmotic_regulation', module: 'ffn_sim/cortex/osmotic_regulation.py', test: 'ffn_sim/tests/test_osmotic_regulation.py',
    cites: 'Olbrich 2000 (membrane Lp~1e-12..1e-13 m/(s·Pa)); Hoffmann 2009 Physiol Rev (RVD/RVI s-min); vant Hoff; Kedem-Katchalsky.',
    seed: 'LEAD PRELIMINARY: V0-setpoint update dV0/dt = -Lp·A·(ΔP_mech-ΔP_target). I checked stability: ∂(ΔP_mech)/∂V0 = +K_vol·V/V0² > 0 (from ΔP = Π0 - K_vol(V-V0)/V0), so ∂V̇0/∂V0 < 0 → STABLE negative feedback to ΔP_target. I believe the SIGN IS CORRECT — try to REFUTE that. Also verify: τ_RVD = V0/(Lp·A·Π_osm) lands in seconds-minutes for Olbrich Lp + c_phys~300 mOsm; the updater reads ev_force.last_pressure correctly; OFF-identity (updater not attached when disabled).',
  },
  {
    name: 'membrane_reservoir', module: 'ffn_sim/cell/membrane_reservoir.py', test: 'ffn_sim/tests/test_membrane_reservoir.py',
    cites: 'Charras 2008 Nat Rev MCB; Tinevez 2009 PNAS (critical cortical tension for blebs); Dai & Sheetz 1999 (MCA); Raucher-Sheetz 1999; Figard 2014 (excess area).',
    seed: 'LEAD PRELIMINARY: rupture-force energy bridge ½k_tether·Δc²=W_MCA·A_bead → F_c=sqrt(2·k_tether·W_MCA·A_bead), grid-invariant (A_bead∝1/N, ΣE=W_MCA·4πR² intensive) — looks sound. σ_crit_bleb and f_excess are None → bleb updater + released_area() raise (honest). VERIFY: W_MCA value vs Dai&Sheetz; is the Bell-Evans rupture (dynamic) actually implemented or only static energy criterion; does the tether topology pick the nearest cortex bead correctly.',
  },
  {
    name: 'cadherin_junction', module: 'ffn_sim/junction/cadherin.py', test: 'ffn_sim/tests/test_cadherin_junction.py',
    cites: 'Rakshit 2012 (X-dimer sliding-rebinding, f0~29 pN); Iturri 2020 (MCF7 de-adhesion 6.5 nN). Reuses validation/cadherin_sliding_rebinding.py + spheroid/cadherin_bonds.py.',
    seed: 'LEAD PRELIMINARY: delegates catch-bond k_off to the validated validation/cadherin_sliding_rebinding.py oracle (no re-implementation — good). Scale bridges: N_cad=F_detach/f0, k_trans=f0/contact_zone_width, k_on=1/mean_lifetime(0). VERIFY: is k_on=1/mean_lifetime(0) actually = effective_k_off(0)? (the comment asserts it). Is the 2-cell partner search (cKDTree) correct and symmetric (no double-binding)? Is k_trans (per-dimer) consistent with the spheroid ensemble k=N_cad·f0/zone?',
  },
  {
    name: 'junctional_actin', module: 'ffn_sim/junction/junctional_actin.py', test: 'ffn_sim/tests/test_junctional_actin.py',
    cites: 'Yonemura 2010 Nat Cell Biol; Buckley 2014 Science (α-catenin/F-actin catch bond); Yao 2014; le Duc 2010.',
    seed: 'LEAD PRELIMINARY: STUB — k_couple + α-catenin catch constants (k_catch0,x_catch,k_slip0,x_slip) + k_on/max_couple_dist all None → builder raises NotImplementedError (correct integrity). VERIFY: the resolver/topology is coherent (planned bonds make sense), OFF-identity holds, and the raise actually fires on EVERY enabled build path (no path that silently runs with None).',
  },
]

const LENSES = [
  { key: 'dimensional', focus: 'DIMENSIONAL ANALYSIS: check the SI units of EVERY formula, derived quantity, and dataclass field. Verify k [N/m], k_angle [J/rad²], energies [J], forces [N], rates [1/s], CFL [s]. Flag any unit mismatch or a field whose comment unit disagrees with how it is computed/used.' },
  { key: 'literature', focus: 'LITERATURE FIDELITY: for each cited constant, verify the VALUE and the FORMULA against the cited paper (use references/ + tag_query.py + knowledge). Flag fabricated/misattributed citations (this project has a history of 3 hallucinated sources), values off from the paper, or a formula that is not the cited papers. Confirm cited values match PINNED GROUND TRUTH.' },
  { key: 'sign_conservation', focus: 'SIGN-SENSE & CONSERVATION: verify every force restores correctly (stretched bond pulls inward), no wrong-sign feedback, Newton-3rd-law respected, no spurious net force/momentum injection (except intentional external fields, which must be documented). Verify catch vs slip direction. Trace at least one concrete stretch/compress case numerically.' },
  { key: 'numerics_magic', focus: 'NUMERICS, MAGIC NUMBERS, GRID-INVARIANCE: hunt for (a) empirical magic numbers (any constant not cited/derived and not None), (b) coarse-graining/series-spring factor errors (per-bond vs whole-fiber, per-bead vs intensive), (c) grid/resolution dependence (does an observable scale with bead count?), (d) CFL correctness and whether stiff elements (MT bending, stiff bonds) are gated. The stress_fibers k_actin is the archetype — find others like it.' },
  { key: 'api_offidentity', focus: 'HOOMD API & OFF-IDENTITY: verify the disabled path is a true no-op (resolve(enabled=False) zeros; builder returns input unchanged; 0 particles/0 bonds). Verify the ENABLED path would actually run on HOOMD 7.0.1 (correct md.bond/angle/force API, correct snapshot extension, bond-type registration) OR correctly raises NotImplementedError when a constant is None. Run a quick `python -c` import + disabled-resolve check. Flag any enabled path that would silently crash or run with a None constant.' },
]

const FINDING_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['module', 'lens', 'findings'],
  properties: {
    module: { type: 'string' }, lens: { type: 'string' },
    findings: { type: 'array', items: {
      type: 'object', additionalProperties: false,
      required: ['severity', 'title', 'location', 'claim', 'evidence', 'recommended_fix'],
      properties: {
        severity: { type: 'string', enum: ['HIGH', 'MEDIUM', 'LOW'] },
        title: { type: 'string' }, location: { type: 'string' },
        claim: { type: 'string' }, evidence: { type: 'string' },
        recommended_fix: { type: 'string' },
      },
    } },
  },
}

const VERDICT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['title', 'verdict', 'final_severity', 'reasoning'],
  properties: {
    title: { type: 'string' },
    verdict: { type: 'string', enum: ['CONFIRMED', 'REFUTED', 'UNCERTAIN'] },
    final_severity: { type: 'string', enum: ['HIGH', 'MEDIUM', 'LOW', 'NONE'] },
    reasoning: { type: 'string' },
    refined_fix: { type: 'string' },
  },
}

const FIX_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['module', 'tests_passed', 'tests_failed', 'fixes_applied', 'deferred_to_pi', 'any_existing_file_touched'],
  properties: {
    module: { type: 'string' },
    fixes_applied: { type: 'array', items: { type: 'object', additionalProperties: false,
      required: ['title', 'what_changed'], properties: { title: { type: 'string' }, what_changed: { type: 'string' } } } },
    deferred_to_pi: { type: 'array', items: { type: 'string' } },
    tests_passed: { type: 'integer' }, tests_failed: { type: 'integer' },
    test_tail: { type: 'string' },
    any_existing_file_touched: { type: 'boolean' },
  },
}

const results = await pipeline(
  MODULES,
  // STAGE 1 — AUDIT: 5 independent lenses per module (barrier within module to collect all lenses)
  async (m) => {
    const lensResults = await parallel(LENSES.map((L) => () => agent(
      `${AUDIT_RULES}\n\n=== MODULE: ${m.name} (${m.module}) ===\nCited literature: ${m.cites}\n${m.seed}\n\n` +
      `=== YOUR LENS: ${L.key} ===\n${L.focus}\n\nRead ${m.module} in full now and return your lens findings.`,
      { label: `audit:${m.name}:${L.key}`, phase: 'Audit', schema: FINDING_SCHEMA }
    )))
    const all = lensResults.filter(Boolean).flatMap((r) => (r.findings || []).map((f) => ({ ...f, lens: r.lens })))
    return { module: m.name, findings: all }
  },
  // STAGE 2 — VERIFY: adversarially refute each finding (kill false positives)
  async (audited, m) => {
    const findings = (audited && audited.findings) || []
    // de-dup by title; verify all HIGH/MEDIUM + a sample of LOW
    const seen = new Set()
    const toVerify = findings.filter((f) => {
      const k = f.title.toLowerCase().slice(0, 60)
      if (seen.has(k)) return false
      seen.add(k); return f.severity !== 'LOW' || true
    })
    const verdicts = await parallel(toVerify.map((f) => () => agent(
      `You are a SKEPTIC verifying a claimed physics bug in ffn_cellsim. ${ENV}.\n` +
      `Module: ${m.module}\nClaimed finding [${f.severity}]: ${f.title}\nLocation: ${f.location}\n` +
      `Claim: ${f.claim}\nEvidence offered: ${f.evidence}\n\n` +
      `${GROUND_TRUTH}\n\nYour job is to REFUTE this finding. Read the actual code at the location. ` +
      `Re-derive the numbers yourself. A finding is CONFIRMED only if you cannot refute it after a ` +
      `genuine attempt. Default to REFUTED if the claim is a misunderstanding of the code, a ` +
      `convention the code documents, or quantitatively wrong. Default to UNCERTAIN if you cannot ` +
      `tell without running. Give the FINAL severity (downgrade if overstated) and a refined_fix if CONFIRMED.`,
      { label: `verify:${m.name}:${f.title.slice(0, 24)}`, phase: 'Verify', schema: VERDICT_SCHEMA }
    )))
    const confirmed = verdicts.filter(Boolean).filter((v) => v.verdict === 'CONFIRMED' && v.final_severity !== 'NONE')
    return { module: m.name, modulePath: m.module, testPath: m.test, confirmed, allVerdicts: verdicts.filter(Boolean) }
  },
  // STAGE 3 — FIX: one owner-agent per module applies CONFIRMED findings + re-runs tests
  async (verified, m) => {
    const confirmed = (verified && verified.confirmed) || []
    if (confirmed.length === 0) {
      return { module: m.name, fixes_applied: [], deferred_to_pi: [], tests_passed: -1, tests_failed: 0,
               test_tail: 'no confirmed findings — module unchanged', any_existing_file_touched: false }
    }
    const list = confirmed.map((c, i) => `${i + 1}. [${c.final_severity}] ${c.title}\n   fix: ${c.refined_fix || '(see reasoning)'}\n   why: ${c.reasoning}`).join('\n')
    return agent(
      `You OWN ${m.module} and ${m.test} for ffn_cellsim. ${ENV}. Apply the CONFIRMED physics fixes ` +
      `below. Edit ONLY these two files (no other existing file — the Lead handles the registry).\n\n` +
      `${GROUND_TRUTH}\n\nCONFIRMED FINDINGS TO FIX:\n${list}\n\n` +
      `RULES: keep the compartment DEFAULT-OFF and OFF-identity intact. NO empirical magic numbers — ` +
      `every new constant MUST be literature-cited (author+year) or grid-invariant-derived; if a needed ` +
      `value is genuinely unknown, set it None + add to PI_DECISIONS + raise NotImplementedError on the ` +
      `enabled path (do NOT invent it) and list it in deferred_to_pi. For the stress_fibers k_actin class ` +
      `of bug: re-derive the stiffness from an elastic-modulus bridge (E·A/ℓ0) with cited material ` +
      `properties (actin bundle: cite the filament modulus + filaments-per-bundle), or set None if the ` +
      `bundle cross-section is unknown. Update the docstring + any affected test. Then run ` +
      `\`python -m pytest ${m.test} -q\` and iterate until green. Confirm via \`git status --porcelain\` ` +
      `that ONLY your two files changed. Return the structured summary.`,
      { label: `fix:${m.name}`, phase: 'Fix', schema: FIX_SCHEMA }
    )
  },
)

return results.map((fix, i) => ({
  module: MODULES[i].name,
  confirmed_count: undefined,
  fix,
}))
