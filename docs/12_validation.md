# 12 — Validation Principles & Assumption Master List

## Validation Hierarchy

The simulation's validity rests on **three pillars**, in order of priority:

1. **Internal consistency** — conservation laws, numerical convergence, dimensional analysis
2. **Literature anchoring** — every physical parameter and force law traces to published literature (preferably IF ≥ 15)
3. **Cross-comparison with PI's experimental data** — overlay only, no fitting

The simulation is publishable based on Pillar 1 + 2 alone. Pillar 3 is supplementary.

## Pillar 1: Internal Consistency Checks

These run automatically and pass/fail in CI:

### Per-Stage Validation (see `10_dev_roadmap.md`)
- **Stage 1a**: Mass conservation (exact in MLS-MPM); momentum conservation (no external force); energy monotonic dissipation under viscoelastic relaxation; spheroid relaxes to round shape.
- **Stage 1a+**: Hertz contact prediction matches initial contact area scaling; spreading proceeds smoothly without instability.
- **Stage 1a++**: Leader cell fraction in 1–5% range; lamellipodia event rate physically reasonable.
- **Stage 1b**: ULA-like vs pV4D4-like phenotypes show qualitative differences in φ trajectory; Δφ between conditions > 0.5 by 24 sim hr.
- **Stage 1c**: Density gradient emerges (edge ρ > core ρ); cortex stiffness gradient consistent.
- **Stage 1d**: Internal vorticity nonzero; Marangoni driving force nonzero; Q-tensor order at edges.
- **Stage 1e**: Sim A reduces to Sim B in highly-symmetric, slow-spreading regimes (validation of derivation).

### Convergence Tests
- Δt → Δt/2: result changes < 5%
- Grid → 2× grid: result changes < 5%
- N points → 2× N points: aggregate metrics (A/A₀, etc.) change < 10%

### Dimensional Analysis Sanity Checks
- All output values within physically reasonable ranges (not 10⁻²⁰ or 10²⁰)
- Spreading rate within range observed in literature for tissue spreading

### Sanity Gate Protocol — mandatory pre-execution review for every new physics/numerics module

This protocol is part of academic honesty in computational work: a numerical scheme that violates dimensional analysis or stability conditions can produce plausible-looking output that is, in fact, meaningless. The protocol enforces that we *prove* — to ourselves and to the reader — that the discretisation is sound *before* a single simulation step runs. The duplicate copy of this protocol in `CLAUDE.md` is the executable directive; this copy in `12_validation.md` is the academic justification.

**When it runs**: after writing or modifying any `acs/physics/*.py`, `acs/boundary/*.py`, `acs/adhesion/*.py`, or numerical-integrator code, **before** the first execution. Modifications to a single line still require re-running the affected sub-checks.

**The five checks**:
1. **Dimensional analysis** — units annotated on every input/output; characteristic scales L, T, S, V, M computed; the dominant non-dimensional numbers (Re, Ca, De, Pe, Ma) reported; for explicit time-stepping schemes the CFL bound must be re-derived against the configured `dt` and a PASS recorded as `dt/dt_CFL = W < 1`.
2. **Boundary cases** — the code must remain well-defined at extreme limits (`N→0`, `N→∞`, `Δt→0`, `Δt→large`, `γ→0`, `K→∞`).
3. **Conservation invariants** — explicit declaration of what is conserved, where dissipation is intentional, and where leaks are possible.
4. **Numerical sanity** — `dx`/`dt`/precision choices justified against the slowest and fastest physically meaningful scales.
5. **Sign / sense check** — every force term annotated with the intuitive direction it pushes.

**Failure handling**: a FAIL halts further code work on that module. The author surfaces the issue to the PI with at least three concrete remediation options (e.g., reduce Δt, switch to implicit, switch to overdamped) and waits for direction before proceeding. The first `mlsmpm.py` draft (April 2026) caught a 16,000× CFL violation via this protocol; that experience is the reason the protocol is now permanent in `CLAUDE.md` Hard Rules.

### Magic-Number Block

The protocol forbids introducing empirical scaling factors / ad-hoc tuning constants ("magic numbers") into physics or numerics code. Before adding any such factor, all three tests below must pass; failing any blocks the change and forces the author to surface the underlying scheme issue to the PI.

1. **Derivable** — can the value be derived from peer-reviewed literature (preferably IF ≥ 15) or from first principles (dimensional analysis, conservation law, asymptotic expansion)?
2. **Grid-invariant** — does the value remain valid as `dx`, `dt`, `grid_n`, `n_particles`, or other discretisation parameters change? Or does each parameter sweep require re-tuning?
3. **Fitting** — was the value chosen to make a specific simulated number match a target (a gate threshold, an experimental datapoint, a prior result)? "Yes" here is automatic disqualification.

The April 2026 `csf_kappa_scale` episode is the documented anti-pattern. After CSF-driven radius drift exceeded the 5% gate, three successive empirical reductions (κ-scale 0.25 → 0.05, then capillary number 0.01 → 0.001) each made the gate "pass" by weakening the surface-tension force, but every one of them failed the three tests above. The principled fix — a finite-difference curvature operator `κ = -∇·n̂` per Brackbill (1992) — has a literature anchor, is grid-invariant up to standard discretisation error, and was not chosen to hit any target. The static-sphere validation (`κ_measured` vs `2/R`, ±10%) became a permanent gate to lock the scheme rather than the parameters.

**Cousin rule — gate semantics are immutable per run**: a gate's tolerance, normalisation, and check window are part of the validation contract authored before the run starts. They are not tuneable in response to a failing result. If a gate appears to be genuinely incorrect (wrong physical interpretation, not just a tight tolerance), the author halts and surfaces the contract change to the PI rather than editing the gate inline.

**Recording the result**: the five-check report lives in the module's docstring under a `Sanity Gate` heading (or as a sibling `<module>_sanity.md` file when too long). It is part of the deliverable, not optional documentation.

## Pillar 2: Literature Anchoring

Every parameter in `02_force_models.md` has at least one IF ≥ 15 reference. Where IF ≥ 15 is not available, lower-IF references are explicitly flagged and noted as a limitation.

### Master Parameter Audit
A spreadsheet (or markdown table) listing:
| Parameter | Value | Reference | IF | Confidence |
|---|---|---|---|---|
| ... |

This is the single source of truth. **No parameter without reference.** If a parameter is needed but no published value exists, document it in this section as a "to-be-determined-by-sweep" parameter.

### Mid-Range vs Extreme Values
Use mid-range literature values as defaults. Use extremes only for sensitivity analysis or specific phenotype matching (e.g., "stiff cortex" upper bound for pV4D4-like).

## Pillar 3: PI's Experimental Data

### Allowed Uses
- **Visualization overlay**: simulation A/A₀ curve plotted alongside PI's experimental Bare/Pre/Lam4 curves
- **Qualitative comparison**: "the simulation produces stronger spreading for high-φ conditions, qualitatively matching PI's Lam4 vs Bare"
- **Order-of-magnitude check**: e.g., "simulation predicts A/A₀ ≈ 5–15 at 80 hr, experimental range 4–33"

### Forbidden Uses
- **Parameter optimization**: gradient descent on simulation parameters to minimize error against experimental data
- **Model selection by fit quality**: choosing between physics implementations based on which matches experimental data better
- **Initial-condition tuning**: adjusting R₀ or other inputs to match specific experimental spheroids

### Why This Matters
The PI's framing: *"the simulation must be defensible independently of my experimental results."* Fitting to data would defeat this — any model can fit any data with enough free parameters. Parameter-free prediction is the gold standard.

### Cross-Check (Distinct from Fitting)
After Stage 1 completes, optionally:
- Run simulation with default literature parameters
- Compare predicted A/A₀ trajectory to experimental envelope
- If within range → publication-strong "independent prediction"
- If outside range → identify missing physics (academic insight)

This is comparison, not fitting.

## Master Assumption List

(Cross-reference: `00_project_vision.md` "Explicit Assumptions & Limitations". This is the consolidated master list.)

### Stage 1 Assumptions (Active)
1. Newtonian-Cauchy continuum mechanics ✓
2. Constant temperature 37°C, constant pH and osmolality ✓
3. No external mechanical perturbations ✓
4. Cell-equivalent material points (~1 point ≈ 1 cell) ✓
5. Constant cell mass (no proliferation/death) ✓
6. All cells viable (no necrotic core) ✓
7. Uniform cell properties (no heterogeneity) ✓
8. Rigid substrate (no compliance) ✓
9. Static ECM (no remodeling) ✓
10. No chemical signaling ✓
11. No nucleus mechanics ✓
12. Mechano-osmotic Tier 2 (phenomenological) ✓
13. External medium: Stokes drag only ✓
14. Initial spheroid radius < 250 μm (otherwise necrosis caveats) ✓
15. Slow spreading vs cortical actin turnover ✓
16. Substrate fully wettable (S > 0) ✓

### Reserved Stages (Inactive Now, Hooks Present)
- Stage 2: Layer 6 chemistry, necrosis, cell heterogeneity, Lam-supplemented substrate
- Stage 3: Other epithelial cell lines, ECM remodeling, substrate stiffness sweep
- Stage 4: Mesenchymal cells (model boundary test), chemical signaling, optional nucleus mechanics

### Permanent Exclusions (Out of Project Scope)
- Quantum / relativistic mechanics
- Sub-molecular detail
- 3D in-vivo invasion (different substrate, different study)
- Drug response / pharmacokinetics
- Single-cell genome-level dynamics

## Citation Standard
- Every claim in code comments → markdown reference like `[Marchetti 2013]`
- Every default parameter → reference in `02_force_models.md`
- Any deviation from defaults → documented in run config + run log
- New phenomenology added during development → added to assumption list with reasoning

## Validation Reports
After each major stage completion, generate `validation_reports/stage_{X}_report.md` summarizing:
- Conservation law check results
- Convergence test results
- Comparison with stage-specific predictions
- New limitations discovered
- Updated assumption list

This builds an audit trail for the project's academic integrity.
