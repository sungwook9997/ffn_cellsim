# ActiveCellSim v2 Phase 1 Plan — Codex Independent Draft

Date: 2026-05-03 KST  
Author: Codex design-discussion session  
Status: Alpha draft for Claude cross-review and PI synthesis  

## 1. Executive Position

Phase 1 should build a cell/ECM mechanobiology reference engine through
Tier 2, not a full spheroid production engine and not a reduced Marangoni
solver. The target deliverable is a mechanistically explicit, measurable,
and visualizable small-scale simulator that can later provide closure
variables for a reduced droplet/Marangoni/active-gel surrogate.

The first implementation target is conservative but useful: a runnable v2
skeleton with measurement-safe geometry, artifact-aware data contracts,
single-cell active-contour numerics, ECM and cluster state schemas, one
validated cortex-tension force term, frame dumps, and a reviewable 3D
visualization stub. Biology modules beyond that are staged, not ignored.

PI experimental data are not used for fitting, anchoring, or plan selection
in this Phase 1 plan. Literature and simulation references define the model
choices. PI data can later be used only as an explicitly separated overlay or
held-out validation target when Rule 11 measurement matching is satisfied.

## 2. Scientific Framing

The central Phase 1 question is:

> Can explicit single-cell and small-cluster mechanobiology on a remodeling
> Col1-like substrate produce measurable outputs that explain when droplet
> microfluidics / Marangoni-like continuum analogies are legitimate, and when
> they fail because of filopodia, focal adhesions, cell-cell contacts, and ECM
> memory?

This framing avoids two failure modes:

- It does not reduce filopodia, focal adhesions, or ECM remodeling to a
  surface-tension gradient before proving such a reduction exists.
- It does not make brute-force cell-resolved simulation the final sweep
  engine. The explicit model is a reference engine for extracting future
  closure variables.

The Phase 1 scope is therefore:

1. Single-cell mechanics on a local ECM patch.
2. Lamellipodia and filopodia event state, with filopodia treated as a core
   observable.
3. Focal adhesion and traction state.
4. Coarse evolving ECM environment: stiffness, density, orientation, and
   remodeling memory.
5. Small-cluster cell-cell contacts and junction state.
6. Minimal cell-cycle and mechanosignaling hooks, default off unless a staged
   activation gate is passed.
7. Output and visualization sufficient for scientific inspection and future
   coarse-graining.

Out of scope for Phase 1:

- Full 1,000+ cell spheroid production.
- Full reduced Marangoni / active-gel / MPM solver.
- Full individual collagen-fiber graph.
- Full oxygen/nutrient reaction-diffusion and necrosis dynamics.
- Detailed Notch/Wnt/TGF-beta pathway networks.
- Any parameter fit to PI data.

## 3. Architecture Overview

```text
acs/v2/
  measurement_boundary.py       # locked design; canonical measurement polygon
  data_contract.py              # artifact-aware MetricSpec / dataset contract
  metrics.py                    # pure metric registry
  single_cell.py                # per-cell state; existing file extended
  ecm_substrate.py              # ECM field / tensor / memory state
  cell_cluster.py               # Tier 2 cell collection and contact graph
  junction.py                   # pairwise junction state
  closure_dump.py               # HDF5 frame / closure output schema
  dynamics/
    active_contour.py           # numerical boundary evolution
    cortex.py                   # cortex / line tension force
    protrusion.py               # lamellipodia / filopodia events
    focal_adhesion.py           # FA state and traction
    ecm_remodeling.py           # traction-to-ECM update
    contact_inhibition.py       # junction/contact to protrusion/polarity bias
    cell_cycle.py               # default off, scalar hook
    mechanosignal.py            # default off, YAP/TAZ proxy hook
  viz/
    plotly_stub.py              # overnight / smoke visualization
    blender_pipeline.md         # later confocal-quality pipeline spec
```

The first coding wave should not try to activate every dynamics module. It
should establish the state, validation, output, and one numerically sane
mechanics path.

## 4. Solver / Representation Choice

### Recommended Phase 1 Primary Solver

Use a 2D/2.5D active-contour / vertex-boundary model for single cells and
small clusters, with a separate ECM field underneath.

Rationale:

- The measurement contract is polygon-based. Active contour naturally emits
  `MeasurementBoundary` without a lossy projection step.
- Filopodia can be represented by local vertex refinement, narrow protrusion
  events, and event metadata rather than full actin-filament simulation.
- Small clusters can share the same boundary representation while adding a
  contact graph.
- It is simpler to Sanity Gate than a phase-field model because area,
  perimeter, and forces are explicit geometric quantities.
- It avoids the lattice anisotropy of Cellular Potts for the first custom
  mechanobiology implementation.

### Solver Alternatives

| Option | Strength | Weakness | Phase 1 decision |
|---|---|---|---|
| Active contour / vertex boundary | Direct measurement boundary, controllable forces, easy HDF5 geometry | Requires careful topology and vertex spacing; filopodia need dynamic refinement | Primary implementation |
| Phase field | Handles topology and smooth interfaces; natural PDE coupling | Grid cost; harder measurement projection; diffuse-interface parameter can become magic | Design reference, not first implementation |
| Cellular Potts / CompuCell3D | Mature cell-cell / growth / chem fields; many biological behaviors | Lattice artifacts; Monte Carlo time mapping; less direct mechanical force accounting | Benchmark/reference, not custom core |
| PhysiCell | Efficient off-lattice multicell, cell cycle/death/microenvironment | Cell shape is too coarse for filopodia and boundary deformation | Future Tier 3 or hybrid bulk engine |
| Chaste / vertex frameworks | Strong tissue mechanics and testing culture | Integration burden and C++ stack; custom ECM/filopodia still needed | Reference for testing and architecture |
| Full discrete fiber / cytoskeletal simulation | Mechanistically rich | Too expensive and too many parameters for Phase 1 | Future focused module only |

Important: biology papers justify signs, timescales, and parameter ranges. They
do not justify numerical stability. The active-contour scheme must pass its own
convergence, symmetry, and area-drift gates.

## 5. State Schema

### 5.1 SingleCellState

Extend the existing `acs/v2/single_cell.py` state with minimal day-one fields.

Core geometry and mechanics:

- `cell_id: str`
- `time_s: float`
- `measurement_boundary: MeasurementBoundary`
- `height_um: Optional[float]`
- `polarity_xy: Optional[tuple[float, float]]`
- `protrusions: list[ProtrusionEvent]`
- `adhesions: list[FocalAdhesionState]`

Day-one slow-biology hooks:

- `cell_state: Literal["alive", "dead"] = "alive"`
- `cell_age_s: float = 0.0`
- `cell_cycle_phase: Optional[Literal["G0", "G1", "S", "G2", "M"]] = None`
- `division_count: int = 0`
- `parent_cell_id: Optional[str] = None`
- `mechanosignal_yap_taz: Optional[float] = None`
- `neighbor_cell_ids: tuple[str, ...] = ()`

Validation:

- age non-negative
- division count non-negative
- YAP/TAZ proxy in `[0, 1]` if present
- neighbor ids non-empty, unique, and not equal to `cell_id`
- dead cells keep geometry for visualization; dynamics modules skip protrusion
  and active force updates unless a later debris model explicitly says
  otherwise

Do not add apoptotic/necrotic states yet. `alive/dead` is enough for Phase 1
schema; apoptosis and necrosis need criteria and are deferred.

### 5.2 ECMSubstrateState

Initial ECM is a coarse local environment, not a full fiber graph.

Fields:

- `origin_um_xy: tuple[float, float]`
- `spacing_um: float`
- `stiffness_kpa: np.ndarray` on an `(nx, ny)` grid
- `ligand_density: np.ndarray` normalized `[0, 1]`
- `fiber_density: np.ndarray` normalized `[0, 1]`
- `orientation_tensor: np.ndarray` shape `(nx, ny, 2, 2)`
- `remodeling_memory: np.ndarray` normalized scalar or tensor field
- `source: Literal["synthetic", "image_derived", "literature_default"]`

Validation:

- positive grid spacing
- all arrays finite and same grid shape
- stiffness non-negative
- densities in `[0, 1]`
- orientation tensor symmetric and bounded

Phase 1 minimal dynamics:

- cell traction aligns `orientation_tensor`
- sustained traction increases remodeling memory and optionally stiffness /
  fiber density within bounded ranges
- ECM state biases protrusion and FA formation

### 5.3 FocalAdhesionState

Existing schema is a useful seed but needs traction and state identity.

Fields:

- `adhesion_id`
- `cell_id`
- `position_um_xyz`
- `age_s`
- `state: Literal["unbound", "nascent", "mature", "slipping", "released"]`
- `maturity: float`
- `bound_fraction: float`
- `traction_force_nN_xy: tuple[float, float]`
- `source: Literal["simulated", "marker_detected", "manual"]`

Phase 1 should express transition rates as ranges or exploratory parameters.
No single universal FA rate constant should be hard-coded.

### 5.4 ProtrusionEvent

Existing schema is close. Add explicit subtype and parent boundary relation
when implemented.

Fields:

- `event_type: Literal["lamellipodium", "filopodium", "retraction"]`
- `boundary_angle_rad`
- `length_um`
- `lifetime_s`
- `width_rad` or `angular_width_rad`
- `tip_position_um_xy` when available
- `associated_adhesion_ids: tuple[str, ...]`

Filopodia are most important for PI. They should not be averaged away in the
reference engine. Their event distribution must be saved.

### 5.5 JunctionState

Pairwise cell-cell contact state:

- `junction_id`
- `cell_a_id`, `cell_b_id`
- `contact_length_um`
- `contact_age_s`
- `junction_maturity: float`
- `cadherin_proxy: Optional[float]`
- `tension_proxy_nN: Optional[float]`
- `contact_inhibition_signal: Optional[float]`

Validation:

- distinct cell ids
- non-negative contact length and age
- scalar proxies in `[0, 1]` when normalized

### 5.6 CellClusterState

Tier 2 state:

- `cells: dict[str, SingleCellState]`
- `junctions: dict[str, JunctionState]`
- `ecm: ECMSubstrateState`
- `time_s`
- optional role labels: leader / follower / internal

Validation:

- all junction ids reference existing cells
- neighbor symmetry
- no self-contact
- cell ids unique
- cluster-level projected area and packing metrics are measurement-derived,
  not substrate-contact proxies

## 6. Dynamics Modules and Gates

### 6.1 Active Contour Core

Equation form: overdamped boundary motion,

```text
xi_i dx_i/dt = F_cortex + F_area + F_protrusion + F_FA + F_contact + F_noise
```

where each term must be separately togglable.

Initial active terms:

1. cortex / line tension
2. area or volume proxy restoration
3. substrate drag / mobility

Later Phase 1 terms:

4. protrusion force
5. FA traction
6. ECM guidance feedback
7. cell-cell contact force

Gates:

- regular polygon relaxation remains centered and symmetric
- area drift under pure cortex is bounded and explained
- zero-force state stays stationary
- force signs documented
- time step stability bound derived from largest stiffness and mobility

### 6.2 Protrusion and Filopodia

Use an event-based stochastic hazard model in Phase 1, not full actin
biochemistry. Hazards may depend on boundary angle, polarity, ECM orientation,
local FA state, and contact state.

Outputs:

- event count per cell/time
- inter-arrival distribution
- angular distribution
- lifetime distribution
- spatial clustering / first peak of angular correlation

Gates:

- no negative event length or lifetime
- refractory period, if introduced, must be derived or exploratory-labeled
- contact edge protrusion hazard decreases when contact inhibition is active
- filopodia event statistics remain available for future stochastic surrogate
  reduction

### 6.3 Focal Adhesion / Molecular Clutch

Start with a state machine plus traction vector, not a detailed molecular
network. FA formation probability depends on ECM ligand density, stiffness, and
local protrusion/contact state. Mature adhesions transmit traction to ECM.

Gates:

- traction is zero for unbound/released states
- traction opposes or anchors cell motion consistently with sign conventions
- FA lifetime distribution is reported, not tuned to projected area
- transition rates are ranges or exploratory-labeled if condition-specific

### 6.4 ECM Remodeling

Phase 1 ECM is a tensor/scalar field with memory.

Candidate update:

```text
dQ_ecm/dt = alignment_rate * f(traction, current_Q) - relaxation_rate * Q_ecm
dE/dt     = stiffening_rate * g(|traction|, memory) - relaxation_or_saturation
```

No exact formula should be implemented before a short ECM Sanity Gate note.

Gates:

- stiffness and density remain non-negative
- normalized fields remain bounded
- zero traction produces no artificial remodeling
- sustained unidirectional traction aligns ECM in the traction direction
- remodeling memory is explicitly stored so the reduced model does not lose
  history dependence

### 6.5 Cell-Cell Contact and Junctions

Small clusters use a contact graph. Contacts alter protrusion hazards and
polarity. Junctions carry maturity/contact-age state.

Gates:

- contact graph symmetric
- contact inhibition only suppresses contact-edge protrusions
- junction maturity does not change without contact
- cell-cell adhesion force is attractive at contact and does not create
  unbounded overlap

### 6.6 Slow Biology Hooks

Cell cycle and YAP/TAZ are schema-first.

Phase 1 default:

- cell age increments
- cycle phase may be absent
- YAP/TAZ proxy may be absent
- no division unless a later cytokinesis geometry gate is approved

Do not activate necrosis, oxygen, nutrient, or full pathway dynamics in the
overnight alpha implementation.

## 7. Reference Map

### 7.1 Biology Anchors

| Module | Primary references | What they justify |
|---|---|---|
| Active matter / future reduction | Marchetti et al. 2013; Prost et al. 2015 | continuum language and future active-gel closure, not solver CFL |
| Cell shape / edge dynamics | Keren et al. 2008; Raynaud et al. 2016 | shape/time coupling and stochastic edge observables |
| Cortex tension | Salbreux et al. 2012; Chugh et al. 2017; Stewart et al. 2011 | cortex is contractile and shape-restoring; pressure/cortex relation |
| Actin protrusion | Pollard & Borisy 2003; Machacek et al. 2009; Tkachenko et al. 2011 | lamellipodial protrusion biology and signaling coupling |
| Filopodia | Mattila & Lappalainen 2008; Mogilner & Rubinstein 2005; Chan & Odde 2008 | architecture, protrusion mechanics, filopodial traction |
| FA / clutch | Plotnikov et al. 2012; Choi et al. 2008; Kanchanawong et al. 2010; Case & Waterman 2015; Elosegui-Artola et al. 2016 | FA state identity, rigidity sensing, clutch thresholds |
| ECM mechanics / remodeling | Levental et al. 2009; Sopher et al. 2014; Hall et al. 2016; Notbohm et al. 2015 | collagen stiffness, remodeling, long-range force transmission |
| Cell-cell mechanics | Maître et al. 2012; Cai et al. 2014; Bazellières et al. 2015 | junction/tension observables, not universal thresholds |
| Contact inhibition | Carmona-Fontaine et al. 2008; Mayor & Carmona-Fontaine reviews | contact-edge protrusion suppression and polarity repulsion |
| Collective traction / jamming | Trepat et al. 2009; Park et al. 2015; Bi et al. 2015 | cluster-level traction, packing, shape-index observables |
| Mechanosignaling | Dupont et al. 2011; Aragona et al. 2013 | YAP/TAZ proxy depends on stiffness, spread area, and cell density |

### 7.2 Simulation / Numerical Anchors

| Candidate / method | Reference examples | Phase 1 use |
|---|---|---|
| Active contour / spreading mechanics | Kass et al. 1988; Stolarska & Rammohan 2017 | geometric/numerical inspiration only |
| Phase-field cell migration | Marth & Voigt 2014; Moure & Gomez 2017; phase-field spreading models | alternative solver and future topology handling |
| Cellular Potts / CPM | Scianna et al. 2013; Niculescu et al. 2015 / CompuCell3D; Wortel et al. 2021 | benchmark, not custom core |
| CPM + ECM fibers/FA | hybrid CPM/MD ECM models; recent collagen-fibril alignment models | evidence that cell-ECM reciprocity matters |
| PhysiCell | Ghaffarizadeh et al. 2018 | future off-lattice multicell / microenvironment engine |
| Chaste | Mirams et al. 2013 | test architecture and cell-based modeling reference |
| Morpheus / CompuCell3D | Starruß et al. 2014; Swat et al. 2012 | rapid virtual tissue modeling reference |
| SimuCell3D | Runser et al. 2024 | future 3D cell-surface mechanics reference |
| Taichi / MPM | Hu et al. 2018; Hu et al. 2019 | v1 reduced-model and performance infrastructure |

References must be checked and expanded before implementing parameter values.
This draft uses references to choose architecture and observables, not to lock
numerical constants.

## 8. Data Output and Closure Dump

Each simulation frame should be serializable to HDF5:

```text
frame_t/
  cells/
    id
    position_xy_um
    velocity_xy_um_s
    boundary_vertices_xy_um
    signed_area_um2
    projected_area_um2
    perimeter_um
    height_um
    polarity_xy
    cell_state
    cell_age_s
    cell_cycle_phase
    mechanosignal_yap_taz
    role_label
  protrusions/
    event_id
    cell_id
    event_type
    angle_rad
    length_um
    width_rad
    lifetime_s
    tip_position_xy_um
  focal_adhesions/
    adhesion_id
    cell_id
    position_xy_um
    state
    age_s
    maturity
    traction_force_nN_xy
  ecm/
    stiffness_kpa
    ligand_density
    fiber_density
    orientation_tensor
    remodeling_memory
  junctions/
    junction_id
    cell_a_id
    cell_b_id
    contact_length_um
    contact_age_s
    maturity
    contact_inhibition_signal
  statistics/
    filopodia_rate
    filopodia_lifetime_hist
    fa_lifetime_hist
    traction_radial_tangential
    velocity_correlation_length
  meta/
    config_hash
    git_commit
    dt_s
    frame_interval_s
    random_seed
```

These outputs are not optional. They make the future active-gel/Marangoni
surrogate derivable rather than post-hoc fitted.

## 9. Visualization Plan

Alpha overnight visualization is a review stub, not final confocal rendering.

Immediate:

- Plotly or Vispy single-cell boundary mesh
- ECM stiffness/density/orientation overlay
- protrusion event markers
- FA traction vectors
- time slider or exported frames
- metric side panel: area, perimeter, filopodia count, FA count

Later:

- Blender headless render pipeline for confocal-like 3D visual quality
- separate layers for cell body, nucleus/height proxy, ECM fibers, FAs,
  lamellipodia, filopodia, and junctions
- side-by-side mechanistic reference and reduced surrogate comparison

Visualization must not drive physics. It consumes HDF5 frame dumps.

## 10. Implementation Milestones

### Overnight Alpha

1. Finish or consume `MeasurementBoundary`.
2. Add artifact-aware data-contract fields if not already complete.
3. Add `metrics.py` skeleton.
4. Add state schemas: `ECMSubstrateState`, `JunctionState`,
   `CellClusterState`, and Phase 1 `SingleCellState` hooks.
5. Add `active_contour.py` with Sanity Gate docstring and one cortex force.
6. Add frame dump scaffold.
7. Add visualization stub.
8. Run available tests; report blockers via MCP.

### Week 1

- Lock state schemas and data contracts.
- Add example synthetic single-cell dataset contract.
- Add active-contour unit tests.
- Add ECM validation tests.
- Add HDF5 round-trip tests.

### Week 2

- Add protrusion/filopodia event model with statistics.
- Add FA state machine and traction vectors.
- Add single-cell ECM feedback loop.

### Week 3

- Add 2-cell contact and junction model.
- Add contact inhibition.
- Add small-cluster state and tests.

### Week 4+

- Add proliferation/cell-cycle activation only after cytokinesis geometry is
  designed.
- Add higher-quality visualization.
- Add closure extraction analysis.

## 11. Sanity Gate Matrix

| Module | Section 1 dimensional | Section 2 boundary | Section 3 invariants | Section 4 numerical | Section 5 sign | Section 6 measurement protocol |
|---|---|---|---|---|---|---|
| MeasurementBoundary | um, um2, um | N>=3, simple polygon | geometry round-trip | float64, min edge | CCW positive area | top-down polygon |
| Active contour | force/mobility/dt | zero force, high curvature | area/perimeter accounting | dt stability | cortex contracts | emits boundary |
| Protrusion | um, s, event rate | zero events, high rate | no negative lifetimes | hazard bounded | protrusion pushes outward | event stats saved |
| FA traction | nN, um, s | unbound/mature/released | no traction when unbound | rate stability | traction anchors/pulls | FA/TFM compatible |
| ECM remodeling | kPa, density, tensor | zero traction, saturated fields | positivity | bounded update | traction aligns/stiffens | ECM fields saved |
| Junction/contact | um, s, proxy | no contact, full contact | graph symmetry | overlap control | contact suppresses edge protrusion | junction metrics |
| Cell cycle hook | s, phase | no phase, dead cell | lineage count | discrete events | no force by itself | diagnostic only |

Any FAIL stops that module and posts `status=blocker` to PI per the approval
queue rule.

## 12. Magic-Number Risk List

High-risk parameters:

- baseline protrusion hazard `lambda0`
- filopodia angular width and lifetime distribution
- refractory time
- cortex tension / line tension
- area-restoring stiffness
- FA on/off/maturation/slip rates
- traction coupling gain
- ECM remodeling alignment/stiffening rates
- contact-inhibition threshold
- YAP/TAZ proxy mapping
- cell-cycle duration and division threshold

Handling:

- Literature-fixed if a defensible range exists.
- Exploratory if not literature-fixed, with explicit label.
- Never fitted to PI experimental area curves.
- No gate tolerance change without PI decision.

## 13. Alpha Coding-Team Handoff Criteria

Do not hand off a design debate. Hand off a compact implementation brief with:

- module name and file path
- exact state fields
- validation rules
- toggles/default activation
- P0 tests
- Sanity Gate text required before first execution
- blocked questions
- what not to implement yet

For overnight alpha, the first brief should be limited to:

1. state schema objects
2. data-contract/metric skeleton
3. active-contour core with cortex force
4. frame dump and stub visualization

Do not hand off filopodia, FA, ECM remodeling dynamics, or contact inhibition
until the plan comparison with Claude converges.

## 14. Open Decisions for PI / Next Design Units

1. ECM representation: coarse tensor field first, full fiber graph later. This
   plan recommends coarse tensor field for Phase 1.
2. Active contour vs phase-field fallback: active contour first, phase-field
   reserved for topology problems.
3. Cytokinesis geometry: separate design unit before activating proliferation.
4. Confocal-quality visualization: separate design unit after the stub.
5. Future reduced surrogate naming: avoid claiming true Marangoni unless the
   stress-gradient mapping is demonstrated; use "heterogeneous active-gel /
   Marangoni-like surrogate" until proven.

## 15. Bottom Line

Phase 1 should be ambitious in state and outputs, conservative in activated
dynamics, and strict in validation. The explicit reference engine is valuable
only if it records enough cell/ECM information to derive a reduced model later.
The overnight alpha can establish the skeleton and visual review path without
pretending to complete the full biology.

## Appendix A. Reference Candidates To Verify Before Parameterization

This list is intentionally split into biological anchors and simulation /
numerical anchors. Biological papers can justify signs, observables, and
parameter ranges. Numerical papers and our own convergence tests justify
discretization and stability.

### A.1 Biology / Mechanobiology

- Marchetti et al. 2013, "Hydrodynamics of soft active matter", Reviews of
  Modern Physics.
- Prost, Julicher, and Joanny 2015, "Active gel physics", Nature Physics.
- Keren et al. 2008, "Mechanism of shape determination in motile cells",
  Nature.
- Raynaud et al. 2016, single-cell edge / shape dynamics, Nature Physics.
- Salbreux, Charras, and Paluch 2012, "Actin cortex mechanics and cellular
  morphogenesis", Trends in Cell Biology.
- Chugh et al. 2017, cortex architecture and tension, Nature Cell Biology.
- Stewart et al. 2011, pressure and cortex in mitotic rounding, Nature.
- Pollard and Borisy 2003, actin-driven motility, Cell.
- Machacek et al. 2009, Rho/Rac/Cdc42 protrusion cycles, Nature.
- Tkachenko et al. 2011, PKA-wave coupling to protrusion, Nature Cell Biology.
- Mattila and Lappalainen 2008, filopodia architecture, Nature Reviews
  Molecular Cell Biology.
- Mogilner and Rubinstein 2005, physics of filopodial protrusion, Biophysical
  Journal.
- Chan and Odde 2008, filopodial traction/probing, Science.
- Plotnikov et al. 2012, FA force fluctuations and rigidity sensing, Cell.
- Choi et al. 2008, nascent FA assembly / maturation, Nature Cell Biology.
- Kanchanawong et al. 2010, FA molecular architecture, Nature.
- Case and Waterman 2015, molecular clutch organization, Nature Cell Biology.
- Elosegui-Artola et al. 2016, talin / clutch rigidity threshold, Nature Cell
  Biology.
- Levental et al. 2009, ECM crosslinking and stiffness in cancer, Cell.
- Modeling / simulation papers on fibrous ECM remodeling by contractile cells
  and long-range force transmission in collagen networks should be added as
  ECM-specific anchors before any remodeling rate is fixed.
- Maître et al. 2012, cell-cell contact mechanics, Science.
- Cai et al. 2014 and Bazellières et al. 2015, junction force / direction
  sensing, Cell / Nature Cell Biology.
- Carmona-Fontaine et al. 2008, contact inhibition of locomotion, Nature.
- Trepat et al. 2009, collective traction, Nature Physics.
- Park et al. 2015 and Bi et al. 2015, jamming / shape index, Nature Materials
  / physical modeling literature.
- Dupont et al. 2011 and Aragona et al. 2013, YAP/TAZ mechanotransduction,
  Nature / Cell.

### A.2 Simulation / Numerical / Software

- Kass, Witkin, and Terzopoulos 1988, active contours ("snakes"), IJCV:
  numerical ancestor only, not a biological cell model.
- Stolarska and Rammohan 2017, cell spreading with focal adhesion placement,
  PLOS ONE.
- Marth and Voigt 2014, phase-field cell motility with signaling networks,
  Journal of Mathematical Biology.
- Moure and Gomez 2017, 3D phase-field cell migration in fibrous networks,
  Computer Methods in Applied Mechanics and Engineering.
- Scianna, Preziosi, and Wolf 2013, Cellular Potts migration in matrix
  environments.
- Wortel et al. 2021, Cellular Potts migration with actin dynamics, PLOS
  Computational Biology / open-access literature.
- Hybrid Cellular Potts + molecular-dynamics ECM models with mechanosensitive
  focal adhesions are important comparison targets for our ECM coupling.
- Ghaffarizadeh et al. 2018, PhysiCell, PLOS Computational Biology.
- Swat et al. 2012, CompuCell3D / multiscale virtual tissue modeling.
- Mirams et al. 2013, Chaste, PLOS Computational Biology.
- Starruß et al. 2014, Morpheus, Bioinformatics.
- Runser, Vetter, and Iber 2024, SimuCell3D, Nature Computational Science.
- Hu et al. 2018 MLS-MPM and Hu et al. 2019 Taichi remain v1 / performance
  infrastructure anchors, not the Phase 1 cell model.

Before any coding brief uses a parameter value from this appendix, add exact
BibTeX entries to `docs/references.bib` or cite the source in the module's
Sanity Gate note.
