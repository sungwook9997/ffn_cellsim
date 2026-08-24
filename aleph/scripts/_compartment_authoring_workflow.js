export const meta = {
  name: 'compartment-authoring',
  description: 'Author 8 missing default-OFF mechanobiology compartments + adversarial verify',
  phases: [
    { title: 'Author', detail: 'one agent per compartment: module + resolver + sanity gate + perf contract + OFF-identity test' },
    { title: 'Verify', detail: 'adversarial audit: magic-number scan, OFF-identity, dimensional sanity, run the test' },
  ],
}

const REPO = '/Users/sw1/ffn_cellsim'

const SHARED_RULES = `
You are authoring ONE new mechanobiology compartment module for ffn_cellsim, a
fine-grained mechanistic HOOMD-blue 7.0.1 single-cell simulator. conda env: ffn_sim
(activate with: source ~/.zshrc; conda activate ffn_sim). Repo root: ${REPO}.

ABSOLUTE HARD RULES (violating any = task failure):
1. FULL-FIDELITY EXPLICIT MECHANISM ONLY. Every element is an explicit HOOMD
   particle/bond/angle. NO lumped proxy force that replaces explicit particles.
   NO mesh-as-physics (no edge/area springs standing in for filaments).
2. DEFAULT-OFF. The resolver MUST accept enabled=False (or a config with
   enabled: false) and return a dataclass with .enabled=False. The topology/force
   builder MUST early-return a no-op (snapshot/sim unchanged) when .enabled is
   False. The module must NOT be wired into cell/cell.py or cell/manifest.py —
   the Lead wires the registry separately. DO NOT EDIT cell/cell.py,
   cell/manifest.py, cell/compartment_registry.py, cortex/cortical_tension.py, or
   ANY existing file. You create ONLY your new module file and your new test file.
3. NO EMPIRICAL MAGIC NUMBERS. Every physical constant MUST have a real
   literature citation (author+year, ideally with the measured value) in a comment,
   OR a grid-invariant derivation. If a constant is genuinely unknown, set the
   field to None, add a one-line entry to a module-level PI_DECISIONS list, and make
   the code path that needs it raise a clear NotImplementedError (keep it disabled).
   DO NOT invent a number to make something work. Honesty over completeness.
4. SI UNITS everywhere, stated in the dataclass field comments (N/m, Pa, m, s, N).
5. NO gamma contamination: your bonds carry a non-cortical load path. Name every
   bond type with your declared prefix so the cortical-tension estimator can
   denylist them. State the prefix in a module-level GAMMA_DENYLIST_PREFIX constant.

REQUIRED MODULE CONTENTS (Google-style docstrings, type hints everywhere):
- Module docstring with: one-line purpose; the explicit mechanism (what particles/
  bonds); a "Sanity Gate" section (dimensional analysis, boundary cases,
  conservation/no-net-force where relevant, sign-sense, CFL note if stiff); a
  "Compartment Performance Contract" section (particle types added; particle count
  at n_fil=1000 and native ~38000 scale; bond/angle count; per-step force yes/no;
  per-batch updater yes/no; uses cpu_local_snapshot yes/no; uses cKDTree/broad-phase
  yes/no; hot-path priority P0/P1/P2/P3; GPU path now CPU/CuPy/native/builtin; native
  ForceCompute candidate yes/no; bottleneck risk); a "References" section.
- A frozen @dataclass(slots=True) ResolvedX with .enabled: bool plus SI-unit fields.
- def resolve_x(cfg: dict, *, <context kwargs like kT, R_cell, dt as needed>) -> ResolvedX
  with validation (raise ValueError on bad/negative params). When cfg has
  enabled: false (or missing), return ResolvedX(enabled=False, ...) with zeros.
- The explicit builder(s): a snapshot-extension function and/or a force/updater
  attach function, mirroring the existing patterns. STUDY these templates first:
  cortex/erm.py (force-only md.force.Custom + resolve + attach + Sanity Gate — THE
  pattern to imitate), cell/nucleus.py (_extend_snapshot_with_nucleus particle
  append + NucleusConfinement), cortex/crosslinkers.py (dynamic bond family +
  Bell-Evans updater), cell/cell.py SubstrateLigandPin (custom Action updater).
  Your builder must be callable in isolation (takes a gsd/HOOMD snapshot or a built
  sim) and be a no-op when disabled.
- A module-level PI_DECISIONS: list[str] (open decisions; empty if none) and
  GAMMA_DENYLIST_PREFIX: str.

REQUIRED TEST FILE (pytest, in ffn_sim/tests/):
- test OFF-identity: resolve with enabled=False returns .enabled False; the builder
  on a disabled config returns its input unchanged / adds zero particles+bonds.
- test dimensional/parameter sanity: resolver computes derived SI quantities
  correctly; negative/zero params raise ValueError.
- test sign-sense where a force exists (e.g. stretched bond pulls inward).
- a no-gamma-contamination assertion: GAMMA_DENYLIST_PREFIX is non-empty and every
  bond type your builder would create starts with it (if your compartment adds
  bonds).
Keep tests import-light and fast (resolver-level + tiny snapshot; do NOT build a
full Cell). Mark any test needing a full HOOMD sim with @pytest.mark.skip if it is
slow/heavy — the resolver + OFF-identity + sign tests must run cheaply.

WORKFLOW:
1. Activate env. Read cortex/erm.py fully and skim the other templates named above.
2. Write the module file and the test file (ONLY these two new files).
3. Run: source ~/.zshrc; conda activate ffn_sim; cd ${REPO}; python -m pytest <your test file> -q
   Iterate until your own test passes (or until a heavy test is correctly skipped).
4. Return the structured summary. Your final text IS the return value (data, not prose).

Return STRICTLY this info: module_path, test_path, n_tests_passed, n_tests_failed,
test_output_tail (last ~15 lines), params (list of {name, value_or_None, unit,
citation}), pi_decisions (list), gamma_denylist_prefix, perf_contract (the contract
fields as an object), optimization_debt (string), and any_existing_file_touched
(MUST be false — if you had to touch one, explain why).
`

const COMPARTMENTS = [
  {
    name: 'ventral_stress_fibers',
    module: 'ffn_sim/cell/stress_fibers.py',
    test: 'ffn_sim/tests/test_stress_fibers.py',
    resolve_fn: 'resolve_stress_fibers',
    prefix: 'sf_',
    brief: `Ventral stress fibers: contractile actomyosin BUNDLES spanning focal
adhesion to focal adhesion along the basal plane. Explicit mechanism: each SF is a
bundle of explicit actin filament beads (sf_actin) cross-linked by alpha-actinin in
a sarcomeric periodicity and loaded by explicit NMII bipolar minifilaments
(sf_myosin) — REUSE the existing minifilament + alpha-actinin mechanism vocabulary
(cortex/myosin.py, cortex/crosslinkers.py) rather than inventing a new lumped
spring. Endpoints anchor at FA clutch particles (requires the fa compartment). The
SF tension is a BASAL-PLANE observable, NOT cortical hoop gamma. Literature anchors
to verify+cite: sarcomeric alpha-actinin periodicity ~0.5-1 um (Tojkander 2012 JCS,
Hotulainen-Lappalainen 2006 JCB); single SF contractile force ~ tens of nN
(Kumar 2006 Biophys J ~10-30 nN; Balaban 2001 traction ~5.5 nN/um^2); NMII per SF.
Provide a builder that lays out n_SF bundles between given FA endpoint positions.
hot-path P2. Flag a PI decision for any anchor you cannot cite.`,
  },
  {
    name: 'linc',
    module: 'ffn_sim/cell/linc.py',
    test: 'ffn_sim/tests/test_linc.py',
    resolve_fn: 'resolve_linc',
    prefix: 'linc_',
    brief: `LINC complex: nesprin-SUN molecular bridges across the nuclear envelope
that mechanically COUPLE the nuclear lamina (nucleus_bead particles) to the cytoplasmic
cytoskeleton (cortex actin / stress fibers / microtubules). Explicit mechanism:
harmonic linker bonds (linc_nesprin) between selected nucleus_bead tags and selected
cytoskeletal particle tags. Requires the nucleus compartment. Literature: nesprin-2
giant is a spectrin-repeat molecular spring (Autore 2013, Arsenovic 2016 FRET tension
~2 pN per nesprin); LINC density thousands per nucleus. The per-bond stiffness of a
single nesprin spectrin-repeat spring is genuinely UNCERTAIN — set k_linc=None,
add a PI decision, and make the enabled path that needs it raise. Provide the
topology (which nucleus beads bond to which cytoskeleton beads, geometry-based pairing)
and a HOOMD-builtin harmonic bond family. hot-path P2.`,
  },
  {
    name: 'intermediate_filaments',
    module: 'ffn_sim/cell/intermediate_filaments.py',
    test: 'ffn_sim/tests/test_intermediate_filaments.py',
    resolve_fn: 'resolve_intermediate_filaments',
    prefix: 'if_',
    brief: `Intermediate filaments (keratin in epithelial MCF7; vimentin): an
extensible, strain-STIFFENING perinuclear cage. Explicit mechanism: if_bead chains
with backbone bonds + if_crosslink bonds, very LOW persistence length (flexible) and
NONLINEAR strain-stiffening (IFs extend >100% before rupture, unlike actin). Literature
to verify+cite: vimentin Lp ~ 1 um, keratin Lp ~ 0.3-0.5 um (Mucke 2004 J Mol Biol,
Lichtenstern 2012); strain-stiffening + extensibility (Block 2018 PRL, Kreplak 2005
J Mol Biol — single IF extends 2-3.5x); small-strain modulus ~ few MPa. The faithful
nonlinear bond needs a tabulated/custom potential; provide the small-strain LINEAR
harmonic bond now (cited modulus) and FLAG the nonlinear strain-stiffening extension
as a documented TODO/PI item (do not fake the nonlinearity). hot-path P2.`,
  },
  {
    name: 'microtubules',
    module: 'ffn_sim/cell/microtubules.py',
    test: 'ffn_sim/tests/test_microtubules.py',
    resolve_fn: 'resolve_microtubules',
    prefix: 'mt_',
    brief: `Microtubules: stiff hollow tubes radiating from the MTOC (centrosome),
load-bearing in COMPRESSION, with optional dynamic instability. Explicit mechanism:
mt_bead chains anchored at a single mtoc particle, with a STIFF backbone bond + a
STIFF bending angle (high flexural rigidity). Literature to verify+cite: persistence
length Lp ~ 1-6 mm (Gittes 1993 JCB 5.2 mm), flexural rigidity EI ~ 2.2e-23 N*m^2
(Gittes 1993); dynamic instability growth/catastrophe/rescue (Mitchison-Kirschner
1984) as an OPTIONAL batch updater. CRITICAL: the very high bending stiffness sets a
TIGHT CFL — compute and REPORT the bending CFL dt and warn it may dominate the cell
dt. Derive k_angle from EI and segment length (grid-invariant, not a magic number).
hot-path P2; flag the CFL/dt risk prominently.`,
  },
  {
    name: 'osmotic_regulation',
    module: 'ffn_sim/cortex/osmotic_regulation.py',
    test: 'ffn_sim/tests/test_osmotic_regulation.py',
    resolve_fn: 'resolve_osmotic_regulation',
    prefix: '',
    brief: `Dynamic osmotic / volume regulation: a time-dependent setpoint on the
EXISTING enclosed_volume turgor force (cortex/enclosed_volume.py) — NOT a new force.
Explicit mechanism: a batch Updater (hoomd.custom.Action) that updates V0(t) and/or
Pi_0(t) of the EnclosedVolumePressure compartment over time, modeling ion-pump +
aquaporin water flux (regulatory volume decrease/increase, RVD/RVI). Read
cortex/enclosed_volume.py to see what setpoint to modulate. Literature to verify+cite:
van't Hoff Pi = R*T*delta_c; membrane water permeability Lp ~ 1e-12 to 1e-13
m/(s*Pa) (Olbrich 2000, dvorak); RVD/RVI timescales seconds-minutes (Hoffmann 2009
Physiol Rev). Requires enclosed_volume. NO new particles, NO new bonds (gamma prefix
is empty — it modulates an existing force setpoint). Provide the updater that steps
the setpoint along a cited flux law; OFF-identity = updater not attached / a no-op
when disabled. hot-path P2.`,
  },
  {
    name: 'membrane_reservoir',
    module: 'ffn_sim/cell/membrane_reservoir.py',
    test: 'ffn_sim/tests/test_membrane_reservoir.py',
    resolve_fn: 'resolve_membrane_reservoir',
    prefix: 'mem_',
    brief: `Membrane reservoir + bleb machinery: the plasma membrane has EXCESS area
(folds, microvilli, caveolae) buffering tension, and BLEBS nucleate when the
membrane detaches from the cortex under pressure. Explicit mechanism: breakable
cortex-membrane tether bonds (mem_tether) whose rupture under load nucleates a bleb,
plus a membrane area reservoir that releases area to buffer tension. Couples to
membrane_surface (requires it) and conceptually to ERM. Literature to verify+cite:
cortex-membrane adhesion energy / bleb nucleation when hydrostatic pressure exceeds
membrane-cortex adhesion (Charras 2008 Nat Rev Mol Cell Biol, Tinevez 2009 PNAS
critical sigma); membrane reservoir excess-area fraction ~ a few % to tens of %
(Figard 2014, Raucher-Sheetz 1999 membrane tension). The bleb nucleation threshold
and reservoir excess-area fraction are uncertain for MCF7 — set them None + PI
decision and keep the bleb path disabled; provide the breakable-tether topology with
a cited adhesion energy where available. hot-path P2.`,
  },
  {
    name: 'cadherin_junction',
    module: 'ffn_sim/junction/cadherin.py',
    test: 'ffn_sim/tests/test_cadherin_junction.py',
    resolve_fn: 'resolve_cadherin_junction',
    prefix: 'cadherin_',
    brief: `Explicit E-cadherin trans-dimer CATCH-BOND cell-cell junction between two
cells. Explicit mechanism: cadherin particles on each cell's surface; dynamic
cadherin_trans bonds across the interface whose lifetime follows the Rakshit 2012
sliding-rebinding catch-slip model. CRITICAL: REUSE the existing faithful
implementations — read validation/cadherin_sliding_rebinding.py (the closed-form
catch-slip k_off, runtime-importable) and spheroid/cadherin_bonds.py
(CadherinBondUpdater, the existing single-cell-aggregate runtime catch-bond updater)
and tests/test_cadherin_catch_bond.py. Adapt them to an explicit TWO-CELL interface
(partner search across the two cells' cadherin particles). NO slip-only shortcut, NO
lumped line tension (full KU-4.2 catch-bond per the architectural rule). Literature
already in repo: Rakshit 2012 SI (k-1^0=30.4/s, x=0.34nm, k+1=5.3, k+2=1985.9,
f0=29.2pN, n=4.8); N_cad ~ Iturri 6.5nN/f0. Put the module in ffn_sim/junction/
(currently empty). hot-path P2. Cite the existing modules you reuse.`,
  },
  {
    name: 'junctional_actin',
    module: 'ffn_sim/junction/junctional_actin.py',
    test: 'ffn_sim/tests/test_junctional_actin.py',
    resolve_fn: 'resolve_junctional_actin',
    prefix: 'junc_actin',
    brief: `Junctional actin coupling: the alpha-catenin / vinculin belt that
mechanically couples the cadherin junction to each cell's actin cortex. Explicit
mechanism: a coupling clutch (junc_actin bonds) between cadherin particles and nearby
cortex actin, with alpha-catenin's FORCE-DEPENDENT vinculin recruitment (a CATCH
behavior: alpha-catenin unfolds under tension and recruits vinculin, reinforcing the
link). Requires cadherin_junction. Literature to verify+cite: alpha-catenin
force-dependent vinculin binding / mechanosensing (Yonemura 2010 Nat Cell Biol,
Buckley 2014 Science — the alpha-catenin/F-actin catch bond, Yao 2014). The catch-bond
coupling constants between cadherin and cortex are NOT well anchored for a
quantitative model — this is a STUB: provide the explicit topology + the resolver +
the Sanity Gate + Performance Contract, set the coupling/catch constants to None with
PI decisions, and make the enabled path raise NotImplementedError. Cite Buckley 2014
for the catch mechanism. hot-path P2.`,
  },
]

const AUTHOR_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['name', 'module_path', 'test_path', 'n_tests_passed', 'n_tests_failed',
             'any_existing_file_touched', 'gamma_denylist_prefix', 'pi_decisions',
             'params', 'optimization_debt', 'self_assessment'],
  properties: {
    name: { type: 'string' },
    module_path: { type: 'string' },
    test_path: { type: 'string' },
    n_tests_passed: { type: 'integer' },
    n_tests_failed: { type: 'integer' },
    test_output_tail: { type: 'string' },
    any_existing_file_touched: { type: 'boolean' },
    gamma_denylist_prefix: { type: 'string' },
    pi_decisions: { type: 'array', items: { type: 'string' } },
    params: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['name', 'unit', 'citation'],
        properties: {
          name: { type: 'string' },
          value: { type: ['number', 'string', 'null'] },
          unit: { type: 'string' },
          citation: { type: 'string' },
        },
      },
    },
    perf_contract: { type: 'object', additionalProperties: true },
    optimization_debt: { type: 'string' },
    self_assessment: { type: 'string', description: 'honest note on what is real vs stubbed/flagged' },
  },
}

const VERIFY_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['name', 'verdict', 'magic_number_violations', 'off_identity_ok',
             'touched_existing_files', 'test_rerun_passed', 'issues'],
  properties: {
    name: { type: 'string' },
    verdict: { type: 'string', enum: ['PASS', 'PASS_WITH_NOTES', 'FAIL'] },
    magic_number_violations: { type: 'array', items: { type: 'string' } },
    off_identity_ok: { type: 'boolean' },
    touched_existing_files: { type: 'array', items: { type: 'string' } },
    test_rerun_passed: { type: 'boolean' },
    issues: { type: 'array', items: { type: 'string' } },
    fixes_applied: { type: 'array', items: { type: 'string' } },
  },
}

const results = await pipeline(
  COMPARTMENTS,
  // Stage 1: author
  (c) => agent(
    `${SHARED_RULES}\n\n=== YOUR COMPARTMENT: ${c.name} ===\n` +
    `Module file to create: ${c.module}\nTest file to create: ${c.test}\n` +
    `Resolver function name: ${c.resolve_fn}\nBond-type GAMMA_DENYLIST_PREFIX: '${c.prefix}'\n\n` +
    `${c.brief}\n\n` +
    `Remember: default-OFF, explicit mechanism, cite or flag every constant, create ONLY ` +
    `your two new files, run your test, return the structured summary.`,
    { label: `author:${c.name}`, phase: 'Author', schema: AUTHOR_SCHEMA }
  ),
  // Stage 2: adversarial verify
  (authored, c) => agent(
    `You are an ADVERSARIAL reviewer for a new default-OFF compartment module in ` +
    `ffn_cellsim (${REPO}). Activate env: source ~/.zshrc; conda activate ffn_sim.\n\n` +
    `Module: ${c.module}\nTest: ${c.test}\nDeclared gamma prefix: '${c.prefix}'\n` +
    `Author self-assessment: ${authored ? authored.self_assessment : '(author failed)'}\n\n` +
    `Audit HARD against these rules and REPORT honestly (and apply small fixes if safe):\n` +
    `1. MAGIC NUMBERS: grep the module for every numeric literal. Each physical constant ` +
    `MUST have a literature citation (author+year) in a nearby comment OR be a documented ` +
    `derivation. List any bare/uncited tuning constant as a violation. (Numerical/array ` +
    `indices, 0.0, 1.0, dimensionless tolerances, dt-safety factors are fine.)\n` +
    `2. OFF-IDENTITY: confirm resolve_*(enabled=False) returns .enabled False and the ` +
    `builder is a verified no-op when disabled (read the code; ideally run a quick python ` +
    `-c check).\n` +
    `3. NO EDITS TO EXISTING FILES: run \`cd ${REPO}; git status --porcelain\` and confirm ` +
    `the ONLY changed/added files are ${c.module} and ${c.test} (plus pycache). List any ` +
    `other touched tracked file as a violation.\n` +
    `4. EXPLICIT MECHANISM: confirm the compartment uses explicit HOOMD particles/bonds, ` +
    `not a lumped proxy. Confirm any bond type starts with '${c.prefix}' (gamma denylist).\n` +
    `5. Re-run the test: \`cd ${REPO}; python -m pytest ${c.test} -q\` and report pass/fail.\n` +
    `If you find a SMALL, safe, clearly-correct fix (typo, missing citation comment you can ` +
    `source, wrong sign, OFF-identity gap), APPLY it to the module/test ONLY (never to other ` +
    `files) and list it in fixes_applied. Do NOT invent physical constants. Return the verdict.`,
    { label: `verify:${c.name}`, phase: 'Verify', schema: VERIFY_SCHEMA }
  ),
)

return results.map((v, i) => ({ compartment: COMPARTMENTS[i].name, verify: v }))
