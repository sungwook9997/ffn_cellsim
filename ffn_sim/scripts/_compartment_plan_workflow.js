export const meta = {
  name: 'compartment-activation-plans',
  description: 'Research + draft a detailed activation plan for each of the 8 default-OFF compartments',
  phases: [
    { title: 'Plan', detail: 'one agent per compartment: research the PI-pending constants + design the activation gate + wiring + tests' },
  ],
}

const REPO = '/Users/sw1/ffn_cellsim'
const ENV = 'source ~/.zshrc; conda activate ffn_sim; cd ' + REPO

const SHARED = `
You are drafting a DETAILED ACTIVATION PLAN for ONE default-OFF compartment of ffn_cellsim
(fine-grained mechanistic HOOMD-blue single-cell simulator). ${ENV}. This is PLANNING, not
implementation — do NOT edit any code. Read the module IN FULL first (esp. its module-level
PI_DECISIONS list and the Compartment Performance Contract docstring), then research and design.

The platform context you must honour (from AGENTS.md / the H.7 contracts):
- Every compartment is wired through configs/mcf7_baseline.yaml -> cell/manifest.resolve_baseline
  -> Cell.build (the hand-wired assembler). A new compartment plugs in via: a manifest
  optional_subsystems block, a resolve_*(), an _extend_snapshot_* and/or attach_* in cell.py, a
  CellBuildOptions flag, and a CompartmentSpec in cell/compartment_registry.py.
- The cortical-gamma estimator (cortex/cortical_tension.py) measures over EXPLICIT cortex bonds.
  Any non-cortical load path (adhesion/junction/SF/LINC/IF/MT/bleb) MUST be denylisted by bond
  prefix. The active-gamma number is gamma(myosin ON) - gamma(myosin OFF); turgor + membrane
  tension are SEPARATE passive channels.
- Physiological-baseline HARD rule: production keeps cytoplasm (65.9 Pa.s), turgor (133 Pa),
  membrane_surface, nucleus ON at setpoint. A new compartment is measured FROM that baseline.
- No empirical magic numbers: a constant must be literature-cited or grid-invariant-derived; if
  genuinely unknown it stays None and the enabled path raises (the current state of several
  constants). Your job is to PROPOSE candidate literature anchors for those None constants so the
  PI has something concrete to ratify or reject — clearly marking confidence and whether the value
  is solid, an order-estimate, or genuinely needs a new measurement.
- Activation gate style (mirror the H.7 Gate-A/B contract): a gate names an OBSERVABLE, a
  LITERATURE BAND, mandatory CONTROLS (OFF/ON differencing, rigid/limit parity, no-contamination
  check), and PASS / PARTIAL / REFUTE thresholds. It is written BEFORE the run and not loosened.

Research tools: ls ${REPO}/ffn_sim/references ; grep the references PDFs;
python ${REPO}/ffn_sim/outputs/tag_kb/tag_query.py "<question>" (the project KB);
and your own domain knowledge. Cite real papers (author+year+journal); never fabricate.

Return a structured plan. Be concrete and quantitative.
`

const COMPARTMENTS = [
  { name: 'ventral_stress_fibers', module: 'ffn_sim/cell/stress_fibers.py',
    note: 'PI-pending: N_filaments (sets mu_SF=N*EA_single, currently None->halt); a distinct sf_myosin_* bond prefix (currently reuses cortex_myosin_* => gamma contamination); bundle bending EI. Requires the fa compartment (FA anchors). Observable = single-SF tension band (Kumar 2006 ~10-30 nN), a BASAL-plane diagnostic, NOT cortical gamma.' },
  { name: 'linc', module: 'ffn_sim/cell/linc.py',
    note: 'PI-pending: k_linc (single-nesprin spectrin-repeat spring stiffness — nesprin force-extension is nonlinear repeat-unfolding, so propose how to anchor a usable effective stiffness or a tabulated potential); f_rest attribution (Arsenovic 2016 relative FRET vs Dejardin 2020 ~8 pN). Requires nucleus. Couples nuclear lamina to cortex/SF/MT.' },
  { name: 'intermediate_filaments', module: 'ffn_sim/cell/intermediate_filaments.py',
    note: 'PI-pending: the NONLINEAR strain-stiffening constitutive law (currently linear-only; nonlinear raises) — propose the tabulated md.bond.Table curve from Kreplak 2005 / Block 2018 force-extension; E_if small-strain modulus (keratin vs vimentin); ratio_xl. Perinuclear cage.' },
  { name: 'microtubules', module: 'ffn_sim/cell/microtubules.py',
    note: 'PI-pending: Y_stretch (axial E*A placeholder); dynamic-instability rates for the cell type (Walker 1988 in-vitro vs in-vivo). CRITICAL: the very stiff bending sets a tight CFL — your plan MUST quantify the dt impact on the full cell and propose how to handle it (segment length choice, sub-stepping, or accept a smaller dt).' },
  { name: 'osmotic_regulation', module: 'ffn_sim/cortex/osmotic_regulation.py',
    note: 'PI-pending: Lp (membrane water permeability — no MCF7 datum; propose a value/source). It modulates the EXISTING enclosed_volume turgor setpoint over time (RVD/RVI). Requires enclosed_volume. Observable = volume-regulation timescale tau_RVD vs Hoffmann 2009 (seconds-minutes).' },
  { name: 'membrane_reservoir', module: 'ffn_sim/cell/membrane_reservoir.py',
    note: 'PI-pending: sigma_crit_bleb (Tinevez 2009 critical cortical tension, MCF7-specific unknown) and f_excess (membrane excess-area fraction). Requires membrane_surface. Observable = bleb nucleation threshold + reservoir tension-buffering. The mem_ gamma-denylist needs wiring into cortical_tension when LIVE.' },
  { name: 'cadherin_junction', module: 'ffn_sim/junction/cadherin.py',
    note: 'Mostly runnable (catch-bond delegates the validated Rakshit oracle; k_trans/r_bind derived). Plan the 2-CELL build (two cells in one box, partner search across the interface), the junction_tension observable, and the SourceEvidence registration for the Iturri anchor. This is the multicell entry point.' },
  { name: 'junctional_actin', module: 'ffn_sim/junction/junctional_actin.py',
    note: 'STUB: k_couple + alpha-catenin catch set (k_catch0,x_catch,k_slip0,x_slip) + k_on/max_couple_dist all None. Propose anchoring from Buckley 2014 (alpha-catenin/F-actin catch) or a sanctioned transfer from the H.3 filamin Pereverzev set. Requires cadherin_junction. Couples the junction to each cell cortex.' },
]

const SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['compartment', 'one_line_goal', 'pi_parameters', 'activation_gate', 'wiring_steps', 'test_plan', 'gpu_and_cfl', 'dependencies', 'ready_when', 'effort_estimate'],
  properties: {
    compartment: { type: 'string' },
    one_line_goal: { type: 'string' },
    pi_parameters: { type: 'array', items: { type: 'object', additionalProperties: false,
      required: ['name', 'status', 'unit', 'proposed_value', 'citation', 'confidence', 'notes'],
      properties: {
        name: { type: 'string' },
        status: { type: 'string', enum: ['SOLID_LITERATURE', 'ORDER_ESTIMATE', 'NEEDS_MEASUREMENT', 'DERIVED'] },
        unit: { type: 'string' },
        proposed_value: { type: 'string' },
        citation: { type: 'string' },
        confidence: { type: 'string', enum: ['HIGH', 'MEDIUM', 'LOW'] },
        notes: { type: 'string' },
      } } },
    activation_gate: { type: 'object', additionalProperties: false,
      required: ['observable', 'literature_band', 'controls', 'pass_partial_refute', 'contamination_check'],
      properties: {
        observable: { type: 'string' }, literature_band: { type: 'string' },
        controls: { type: 'array', items: { type: 'string' } },
        pass_partial_refute: { type: 'string' }, contamination_check: { type: 'string' },
      } },
    wiring_steps: { type: 'array', items: { type: 'string' } },
    test_plan: { type: 'array', items: { type: 'string' } },
    gpu_and_cfl: { type: 'string' },
    dependencies: { type: 'array', items: { type: 'string' } },
    ready_when: { type: 'string' },
    effort_estimate: { type: 'string' },
  },
}

const plans = await parallel(COMPARTMENTS.map((c) => () => agent(
  `${SHARED}\n\n=== COMPARTMENT: ${c.name} (${c.module}) ===\n${c.note}\n\n` +
  `Read ${c.module} now (esp. PI_DECISIONS + the Performance Contract), research the PI-pending ` +
  `constants, and return the detailed activation plan.`,
  { label: `plan:${c.name}`, phase: 'Plan', schema: SCHEMA }
)))

return COMPARTMENTS.map((c, i) => ({ compartment: c.name, plan: plans[i] }))
