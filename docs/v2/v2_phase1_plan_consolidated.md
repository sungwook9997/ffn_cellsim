# ActiveCellSim v2 Phase 1 Consolidated Plan

Date: 2026-05-03 KST
Authors: Codex + Claude design-discussion consolidation
Inputs:
- `docs/v2/v2_phase1_plan_codex.md`
- `docs/v2/v2_phase1_plan_claude.md`
- `docs/v2/10_dev_roadmap_v2.md`
- `AGENTS.md` project rules

Status: alpha consolidation for implementation-work handoff. This document is
the current design source of truth for Project alpha unless superseded by a
later PI decision or blocker.

## 0. Consolidation Result

Codex and Claude converged on the core architecture.

Phase 1 builds a cell/ECM mechanobiology reference engine through small
clusters. It does not build a full spheroid production solver, does not build
the reduced Marangoni/active-gel surrogate yet, and does not use PI
experimental data for fitting, anchoring, or parameter selection.

The explicit reference engine is still important to the droplet/Marangoni
vision. Its job is to generate closure variables that a future reduced model
can try to coarse-grain: traction radial/tangential components, ECM memory,
protrusion statistics, FA lifetimes, contact inhibition strength, and
cell-cell junction observables.

## 1. Certainty Tiers

This plan separates decisions by certainty so implementation does not
accidentally lock weak constants under overnight pressure.

### Tier A: Locked for alpha implementation

- v2 code lives under `acs/v2/`.
- `MeasurementBoundary` remains the canonical measurement polygon:
  `world_um_y_up`, CCW, float64, read-only.
- Simulation uses a mutable `ActiveContourState` that can export a
  `MeasurementBoundary`; measurement state and simulation state stay separate.
- Phase 1 primary representation is active contour / adaptive vertex boundary.
- ECM is present from day one as a coarse continuous field: stiffness, ligand
  density, fiber density, alignment tensor, and accumulated traction history.
  Full explicit collagen-fiber graph is deferred.
- Filopodia remain explicit events in state/output. They are not hidden inside
  a mean protrusion coefficient.
- Focal adhesions and per-FA traction are central Phase 1 biology, but their
  dynamics are staged after schema and core numerics.
- Tier 2 small clusters use a contact graph plus junction state; cells do not
  merge polygons when they contact.
- Cell-cycle and YAP/TAZ hooks exist as schema fields. They are default OFF
  until their own gates exist.
- O2/nutrient, necrosis/apoptosis-specific dynamics, full signaling pathways,
  and the reduced Tier 3 surrogate are deferred.
- HDF5 frame/closure dump is a first-class output.
- Visualization consumes frame dumps. It must not drive physics.

### Tier B: Design inventory / literature-range candidates

These are scientifically plausible terms and parameters, but not accepted
module constants.

- Cortex/line tension
- Area/volume restoring pressure
- Substrate drag / mobility
- Lamellipodia event hazard and force
- Filopodia event hazard, length, lifetime, and tip force
- FA clutch transition rates and traction magnitudes
- ECM stiffening, density change, and alignment rates
- Junction adhesion/tension and contact-inhibition strength
- Polarity persistence / reorientation timescale
- YAP/TAZ mechanosignal coupling

Each entry must be represented as a named candidate range or exploratory
parameter until module implementation writes the required Sanity Gate.

### Tier C: Deferred until module Sanity Gate

Do not implement executable dynamics for these until the module-specific
Sanity Gate is written and passes:

- Numeric timestep values such as `dt_cell = 10 ms`.
- ECM unit conversions between kPa, pN/um2, traction density, and grid-cell
  storage.
- ECM remodeling rate constants such as `k_stiffen`, `k_align`, and isotropic
  relaxation.
- Contact inhibition thresholds.
- Stochastic protrusion baseline rates and refractory times.
- Any `epsilon_*`, `lambda_*`, or gain-like coefficient.
- Any gate tolerance or acceptance window.

If a value fails the Magic-Number Block, implementation pauses that module and
sends `status=blocker` to PI. It must not be tuned to pass a test or match
PI data.

## 2. Scientific Scope

Phase 1 asks:

1. Can a single active-contour cell on a remodeling Col1-like substrate produce
   measurable shape, protrusion, traction, and ECM-memory observables?
2. Can explicit filopodia statistics be preserved without making the solver a
   full actin-filament simulation?
3. Can FA traction and substrate remodeling be recorded in a form usable by a
   future reduced active-gel/Marangoni surrogate?
4. Can small-cluster contact graphs capture cell-cell contact, junction age,
   contact inhibition, and leader/follower structure?
5. Which variables are actually coarse-grainable, and which remain
   cell-resolved failure modes for a radial or droplet-like approximation?

Phase 1 non-goals:

- No full 1,000+ cell spheroid production.
- No MPM or full hydrodynamic Tier 3 solver.
- No full collagen fiber graph.
- No necrotic core, nutrient diffusion, or apoptosis-specific biology.
- No experimental-data fitting.

## 3. Recommended Package Structure

```text
acs/v2/
  measurement_boundary.py
  data_contract.py
  metrics.py
  single_cell.py
  active_contour.py
  ecm_substrate.py
  focal_adhesion.py
  protrusion.py
  junction.py
  cell_cluster.py
  output/
    frame_dump.py
  dynamics/
    active_contour.py
    cortex.py
    protrusion.py
    focal_adhesion.py
    ecm_remodeling.py
    contact_inhibition.py
    cell_cycle.py
    mechanosignal.py
  viz/
    stub3d.py
    blender_pipeline.md
```

The first alpha coding wave should prioritize schemas, validation, frame dump,
and a reviewable visualization stub. Only active-contour cortex/area numerics
are eligible for executable physics in alpha, and only if the Sanity Gate can
be written before first execution.

## 4. Solver Decision

Use active contour / adaptive vertex boundary as the Phase 1 primary solver.

Reasons:

- It directly emits measurement polygons.
- It keeps area, perimeter, curvature, and boundary forces inspectable.
- It supports local refinement for filopodia without making the whole domain
  high resolution.
- It makes small-cluster contact graph construction explicit.
- It avoids early lattice/time-calibration issues from Cellular Potts methods.

Alternatives remain comparison references:

| Method | Use in Phase 1 |
|---|---|
| Phase field | Design reference for smooth topology and diffuse-interface ECM coupling. Not first core because interface width can become a magic parameter. |
| Cellular Potts / CompuCell3D | Biological benchmark/reference. Not first core because Monte Carlo time and lattice anisotropy complicate force accounting. |
| PhysiCell | Future bulk/Tier 3 comparison, not filopodia-resolving. |
| Chaste / vertex frameworks | Testing and architecture reference. Integration burden too high for alpha. |
| Full fiber/cytoskeleton simulation | Deferred focused module only. |

Numerical stability is not inherited from the biology papers. The active
contour scheme must pass convergence, symmetry, area-drift, and timestep gates
using its own implemented force law and mobility.

### 4.1 Adaptive vertex density

Adaptive vertices are part of the representation, not a later optimization.

Alpha schema should allow:

- baseline boundary spacing for ordinary cell edge segments
- local refinement tags around active filopodia and sharp protrusions
- vertex insertion when local edge length exceeds the registered maximum
- vertex removal when refinement is no longer required and geometry remains
  simple
- a per-cell vertex cap to prevent runaway refinement

Exact spacing values are not locked. Spacing is derived from the smallest
physical feature the module claims to resolve and from the active-contour
stability gate.

## 5. State Schema

### 5.1 `SingleCellState`

Core fields:

- `cell_id`
- `time_s`
- `measurement_boundary: MeasurementBoundary`
- `height_um: float | None`
- `polarity_xy: tuple[float, float] | None`
- `protrusions: tuple[ProtrusionEvent, ...]`
- `adhesions: tuple[FocalAdhesionState, ...]`

Slow-biology hooks:

- `cell_state: Literal["alive", "dead"]`
- `cell_age_s: float`
- `cell_cycle_phase: Literal["G0", "G1", "S", "G2", "M"] | None`
- `division_count: int`
- `parent_cell_id: str | None`
- `mechanosignal_yap_taz: float | None`
- `neighbor_cell_ids: tuple[str, ...]`

Validation:

- age and division count non-negative
- YAP/TAZ in `[0, 1]` if present
- neighbor ids unique and not self
- dead cells retain geometry for visualization but active dynamics skip them

### 5.2 `ActiveContourState`

Mutable simulation-only state:

- `cell_id`
- `vertices_um_xy: np.ndarray`
- `vertex_rest_spacing_um`
- `target_area_um2`
- `height_um`
- optional per-vertex mobility, force accumulator, and refinement tags

Required methods:

- `validate()`
- `to_measurement_boundary()`
- `from_measurement_boundary()`
- adaptive insert/remove helpers

The state may be mutable for time stepping. Its exported
`MeasurementBoundary` is immutable and measurement-safe.

### 5.3 `ECMSubstrateState`

Initial continuous field representation:

- `origin_um_xy`
- `spacing_um`
- `stiffness_kpa[nx, ny]`
- `ligand_density[nx, ny]` in `[0, 1]`
- `fiber_density[nx, ny]` in `[0, 1]`
- `orientation_tensor[nx, ny, 2, 2]`
- `accumulated_traction_nNs_per_um2[nx, ny]` or tensor form
- `source: Literal["synthetic", "image_derived", "literature_default"]`

Validation:

- finite arrays with matching grid shape
- positive grid spacing
- non-negative stiffness
- bounded densities
- orientation tensor symmetric, finite, and bounded
- accumulated traction finite and non-negative for scalar storage

No exact ECM update coefficient is locked here. ECM dynamics are staged after
the ECM Sanity Gate documents units and Rule 10 comparisons cleanly.
`accumulated_traction_nNs_per_um2` is deliberately more concrete than
`remodeling_memory`: it records the history variable's unit and prevents a
dimensionless catch-all state from hiding a scaling choice.

### 5.4 `ProtrusionEvent`

Fields:

- `event_id`
- `cell_id`
- `event_type: Literal["lamellipodium", "filopodium"]`
- `start_time_s`
- `end_time_s | None`
- `boundary_angle_rad`
- `root_position_um_xy`
- `tip_position_um_xy | None`
- `direction_um_xy`
- `length_um`
- `width_um | None`
- `width_rad | None`
- `state: Literal["growing", "stalled", "retracting", "ended"]`
- `force_candidate_nN | None`
- `associated_adhesion_ids: tuple[str, ...]`
- `source: Literal["simulated", "detected", "manual"]`

Filopodia need explicit lifetime, length, direction, and attachment/traction
metadata. The event-rate model is a simulation ansatz and must be marked as
such in code.

### 5.5 `FocalAdhesionState`

Fields:

- `adhesion_id`
- `cell_id`
- `position_um_xy`
- `age_s`
- `state: Literal["unbound", "nascent", "mature", "slipping", "released"]`
- `maturity`
- `bound_fraction`
- `traction_force_nN_xy`
- `linked_protrusion_id | None`
- `source: Literal["simulated", "marker_detected", "manual"]`

Transition rates are condition-dependent ranges, not single universal
constants. Plotnikov/Choi/Kanchanawong/Case-Waterman/Elosegui-Artola anchor
state identity and force-sensing logic, not a complete parameter table.

### 5.6 `JunctionState` and `CellClusterState`

`JunctionState`:

- `junction_id`
- pair of cell ids
- contact edge indices or arc spans
- contact length
- age
- maturity
- `cadherin_proxy: float | None`
- `tension_proxy_nN: float | None`
- `contact_inhibition_signal: float | None`
- state: free/contacting/nascent/cortex-coupled/remodeling/separating

`CellClusterState`:

- `cells: dict[str, SingleCellState]`
- `ecm: ECMSubstrateState`
- junction list
- symmetric neighbor graph
- optional role labels such as leader/follower/internal
- lineage metadata if division is active

Cells remain separate polygons. Contact creates a graph edge and junction
state, not a geometry merge.

## 6. Dynamics Inventory

The following dynamics are design inventory. Implementation order is staged.

### 6.1 Active contour cortex and area pressure

Eligible for alpha executable implementation after Sanity Gate:

- line-tension/cortex term: shape-restoring, contractile along boundary
- area/2.5D volume restoring term: outward when area below target, inward when
  area above target
- overdamped update using explicit mobility/friction

The timestep must be derived from implemented stiffness and mobility. Any
fixed value such as 10 ms is a test candidate only.

### 6.2 Protrusion and filopodia

Staged after alpha:

- lamellipodia: broad edge event
- filopodia: narrow probing event, adaptive vertex refinement, explicit
  event statistics
- event hazards biased by polarity, free/contacted edge state, ECM alignment,
  and mechanosignal if active

The inhomogeneous Poisson process is a modeling ansatz. It is acceptable only
with pre-registered observables: inter-arrival distribution, lifetime
distribution, angular distribution, and spatial autocorrelation.

### 6.3 Focal adhesion / traction

Staged after active contour gates:

- FA state machine
- per-FA traction vector
- substrate reaction-force deposition
- traction radial/tangential decomposition for closure dumps

No traction without attachment. Action-reaction with the substrate must be
explicitly accounted for in diagnostics.

### 6.4 ECM remodeling

Staged after ECM units pass Rule 10:

- sustained traction increases remodeling memory
- traction can align ECM tensor
- optional bounded stiffness/density change
- remodeled ECM biases future protrusion and FA formation

All kPa-to-force-density conversions must be written as a clean unit chain in
the ECM Sanity Gate. Do not include a self-correcting or partially checked
conversion in executable code.

### 6.5 Cell-cell contact and junctions

Staged after single-cell plus FA/ECM skeleton:

- contact detection
- junction creation and aging
- contact-edge protrusion suppression
- polarity reorientation / contact inhibition
- neighbor graph metrics

Junction lifetimes and contact-inhibition thresholds are measurement- and
cell-type-dependent. Use ranges and observable distributions, not universal
hard thresholds.

### 6.6 Slow biology hooks

Default OFF:

- cell age
- phase label
- division count / parent id
- YAP/TAZ scalar in `[0, 1]`

Do not activate proliferation or YAP feedback until their Sanity Gates define
timescales, bounds, measurement protocol, and lineage accounting.

## 7. Output and Closure Dump

Use HDF5 frames plus JSON metadata. Every run must store full config and git
commit hash when executable simulations begin.

Proposed structure:

```text
runs/<run_id>/
  config.yaml
  metadata.json
  frames/
    frame_000000.h5
    frame_000001.h5
  metrics.csv
  stats/
    protrusion_events.csv
    fa_lifetimes.csv
    junction_events.csv
```

Frame HDF5:

```text
/meta/time_s
/cells/id
/cells/state
/cells/centroid_um_xy
/cells/polarity_xy
/cells/height_um
/cells/boundary_vertices_flat_um_xy
/cells/boundary_offsets
/cells/area_um2
/cells/perimeter_um
/protrusions/*
/focal_adhesions/*
/ecm/stiffness_kpa
/ecm/ligand_density
/ecm/fiber_density
/ecm/orientation_tensor
/ecm/accumulated_traction_nNs_per_um2
/junctions/*
/closure/traction_radial_nN
/closure/traction_tangential_nN
/closure/ecm_alignment_order
/closure/protrusion_angle_histogram
```

CSR-style flat boundary storage is preferred for variable vertex counts.

## 8. Data Contract and Artifact Kinds

The previously locked artifact enum had 8 members:

- `CSV_TABLE`
- `RAW_FRAMES`
- `SEGMENTATION_MASK`
- `BOUNDARY_CONTOURS`
- `TRACKING_TABLE`
- `MARKER_CHANNEL`
- `TFM_FIELD`
- `EVENT_ANNOTATIONS`

Phase 1 expands this enum to 9 members:

- `ECM_FIELD`

`ECM_FIELD` means an ECM substrate state field artifact, either simulation
output or a controlled in vitro / in silico ECM input. This is a deliberate
change from the previous 8-member contract, so implementation must update
tests and docs explicitly instead of treating it as an incidental enum edit.

## 9. Visualization

Alpha visualization is a reviewable stub:

- input: HDF5 frame dump
- render: single-cell or small-cluster polygons extruded by height
- substrate: floor plane with ECM scalar/tensor overlay if available
- output: PNG sequence and/or lightweight HTML

Confocal-quality Blender visualization is a separate design unit. It should
consume the same HDF5 frame dumps and must not add physics assumptions.

## 10. Sanity Gate Matrix

Every new or modified physics/numerics module must record the six Sanity Gate
items before first execution.

| Module | Gate focus |
|---|---|
| `measurement_boundary` | coordinate convention, orientation, area/perimeter, measurement protocol |
| `data_contract` | artifact requirements, duplicate metric keys, dataset/metric compatibility |
| `metrics` | pure functions, source modality, boundary validation |
| `ecm_substrate` | units, positivity, tensor bounds, Rule 10 unit chain, ECM artifact protocol |
| `active_contour` | force units, timestep bound, symmetry relaxation, area drift, sign of each force |
| `protrusion` | event distributions, RNG reproducibility, length/lifetime bounds, outward/retraction signs |
| `focal_adhesion` | state transitions, no traction unbound, action-reaction, dt vs transition times |
| `junction` / `cell_cluster` | symmetric neighbor graph, no self-neighbor, contact force sign, graph consistency |
| `cell_cycle` | age monotonicity, phase transitions, lineage accounting, area split if division active |
| `mechanosignal` | `[0,1]` bounds, timescale separation, no PI-data tuning |
| `frame_dump` | schema fidelity, round-trip, frame atomicity, complete artifact metadata |

Any Sanity Gate failure pauses that module and is reported as `status=blocker`.

## 11. Magic-Number Policy

Consolidated correction to both drafts:

- Do not mark planned constants as `PASS` until the actual module Sanity Gate
  derives or cites the value and verifies grid/timestep independence.
- Use labels:
  - `literature-direct`
  - `literature-range candidate`
  - `exploratory default`
  - `requires Sanity Gate derivation before execution`
  - `blocked: fails Magic-Number Block`
- Never tune any value to PI area trajectories or other PI experimental data.
- Never change a gate tolerance to make a run pass.

Known high-risk constants:

- baseline protrusion event rate
- protrusion angular width
- protrusion refractory time
- filopodia force/lifetime/length distributions
- FA transition rates
- FA slip thresholds
- cortex and area stiffnesses
- substrate drag/mobility coefficient
- ECM stiffening/alignment/density rates
- contact inhibition strength
- polarity persistence timescale
- YAP/TAZ coupling strength
- every timestep and every gate tolerance

Candidate numeric inventory, preserved from the independent drafts but not
locked as code defaults:

| Quantity | Candidate label | Notes |
|---|---|---|
| filopodia tip force | `literature-direct` | If using direct measured pN-scale force papers, cite exact source in module gate. |
| cortex tension | `literature-range candidate` | Chugh/Salbreux-style cortex anchor; exact value is cell-state and assay dependent. |
| FA traction / stiffness sensing | `literature-range candidate` | FA 5-stack anchors state identity and ranges; transition rates stay ranges. |
| area stiffness `K_A` | `requires Sanity Gate derivation before execution` | Must be derived from chosen area constraint and mobility, not copied from tissue vertex models. |
| substrate drag coefficient | `requires Sanity Gate derivation before execution` | Must be dimensionally checked with implemented overdamped update. |
| ECM stiffening / alignment / density rates | `requires Sanity Gate derivation before execution` | Rule 10 unit chain required before any executable update. |
| protrusion angular width / refractory time | `exploratory default` | Modeling choices; pre-register observables before use. |
| contact inhibition strength | `exploratory default` | Use signal scalar and hazard ratio; no universal threshold. |
| polarity persistence | `requires Sanity Gate derivation before execution` | Gate should use autocorrelation / displacement-polarity cosine window, not fitted `tau_p`. |
| `dt_cell` | `requires Sanity Gate derivation before execution` | Candidate test points are allowed only after deriving `dt < mobility/stiffness` bound for the implemented force law. |
| `dt_ecm` | `requires Sanity Gate derivation before execution` | Operator-splitting timescale must be checked against ECM update rates. |

## 12. References

Use these as anchors for signs, observables, and candidate ranges. Exact
parameter extraction must verify the paper and be cited in the module Sanity
Gate or `docs/references.bib`.

Layer 1:

- Marchetti et al. 2013, Reviews of Modern Physics: active-matter scaling and
  future continuum language.
- Prost, Julicher, and Joanny 2015, Nature Physics: active-gel physics and
  future surrogate language.
- Keren et al. 2008, Nature; Raynaud et al. 2016, Nature Physics: single-cell
  shape and stochastic edge dynamics.
- Chugh et al. 2017, Nature Cell Biology: single-cell cortex architecture.
- Salbreux, Charras, and Paluch 2012, Trends in Cell Biology: cortex mechanics
  review support.
- Machacek et al. 2009, Nature; Tkachenko et al. 2011, Nature Cell Biology:
  protrusion-cycle signaling anchors.
- Mattila and Lappalainen 2008, Nature Reviews Molecular Cell Biology;
  Chan and Odde 2008, Science: filopodia architecture and traction.
- Plotnikov et al. 2012, Cell; Choi et al. 2008, Nature Cell Biology;
  Kanchanawong et al. 2010, Nature; Case and Waterman 2015, Nature Cell
  Biology; Elosegui-Artola et al. 2016, Nature Cell Biology: FA/clutch stack.
- Levental et al. 2009, Cell: collagen/ECM stiffness and cancer context.

Layer 2:

- Cai et al. 2014, Cell; Bazellieres et al. 2015, Nature Cell Biology:
  junction force/direction-sensing observables.
- Maitre et al. 2012, Science: cell-cell contact mechanics, not single-cell
  cortex default.
- Carmona-Fontaine et al. 2008, Nature: contact inhibition of locomotion.
- Trepat et al. 2009, Nature Physics: collective traction.
- Park et al. 2015, Nature Materials; Bi et al. 2014/2015, physical modeling
  literature: jamming/shape-index context.
- Dupont et al. 2011, Nature; Aragona et al. 2013, Cell: YAP/TAZ hooks.

Numerical/simulation references:

- Kass, Witkin, and Terzopoulos 1988, IJCV: active contours as numerical
  ancestor only.
- Stolarska and Rammohan 2017: active-contour-like cell spreading with focal
  adhesion placement, comparison reference only.
- Marth and Voigt 2014; Moure and Gomez 2017: phase-field migration/ECM
  references, comparison only.
- Scianna/Preziosi/Wolf 2013 and Wortel et al. 2021: Cellular Potts migration
  and actin-coupled comparison references.
- Ghaffarizadeh et al. 2018: PhysiCell, future bulk/Tier 3 comparison.
- Swat et al. 2012: CompuCell3D, virtual-tissue comparison reference.
- Mirams et al. 2013: Chaste, architecture/testing reference.
- Starruss et al. 2014: Morpheus, simulation-platform comparison reference.
- Runser/Vetter/Iber 2024: SimuCell3D candidate comparison reference; verify
  exact citation before formal use.
- Phase-field, Cellular Potts, PhysiCell, CompuCell3D, Chaste, Morpheus, and
  SimuCell3D papers do not authorize our timestep, vertex spacing, or force
  constants.

References still requiring verification before parameterization:

- ECM plastic remodeling rate papers beyond Levental/Hall/Trichet candidates.
- Any paper used to claim a numeric ECM stiffening/alignment rate.
- Any exact cell-cycle duration review for MCF7-like cells.
- Any active-contour-with-internal-pressure cellular paper not already in
  `docs/references.bib`.
- Eichinger 2020 ECM remodeling candidate, Pavlov 2024 cell-cycle candidate,
  and Wang 2019 active-contour/internal-pressure candidate until exact venue
  and title are verified.

## 13. Alpha Implementation Handoff

Alpha handoff is split into P0 and P1. P0 can proceed without executable
physics. P1 proceeds only if its Sanity Gate is written cleanly before first
execution.

### P0: schema, contract, output, and visualization

1. Update `acs/v2/data_contract.py` with artifact-aware `MetricSpec`,
   `ImagingDatasetSpec`, dataset-level validation, and duplicate key
   `(name, measurement_modality)`.
2. Add/verify `acs/v2/metrics.py` with pure registry functions for boundary
   projected area and perimeter.
3. Add schema-only `ECMSubstrateState`, `FocalAdhesionState`,
   `ProtrusionEvent`, `JunctionState`, and `CellClusterState` validation.
4. Extend `SingleCellState` with slow-biology hooks and neighbor ids.
5. Add minimal HDF5 frame dump writer for one cell plus ECM metadata.
6. Add stub visualization that reads the frame dump and renders synthetic or
   schema-only frames.
7. Add focused tests for validation and frame round-trip.

P0 must not implement filopodia, FA, ECM remodeling, junction dynamics, or
contact inhibition dynamics.

### P1: active-contour core, only if Sanity Gate passes

1. Add `ActiveContourState` with export/import to `MeasurementBoundary`.
2. Add cortex/line-tension and area-restoring terms only.
3. Add Sanity Gate docstring before execution:
   - units for vertices, force, mobility, time
   - timestep bound from actual stiffness and mobility
   - boundary cases N<3, NaN, self-intersection, tiny edges
   - area drift and symmetry preservation
   - force sign/sense check
   - measurement protocol for exported boundary
4. Add tests for symmetric polygon relaxation, area drift, and timestep
   stability.

If the timestep bound or unit derivation fails, P1 stops and P0 remains the
alpha deliverable.

### Explicit alpha non-scope

- filopodia event scheduler
- FA molecular clutch dynamics
- ECM remodeling dynamics
- cell-cell junction dynamics
- contact inhibition dynamics
- cell-cycle or YAP feedback
- Blender confocal-quality render
- PI experimental data comparison

## 14. Coding-Team Escalation Rules

- Sanity Gate FAIL: send `status=blocker` to PI and pause that module.
- Magic-Number Block FAIL: send `status=blocker` to PI and pause that module.
- The `ArtifactKind.ECM_FIELD` expansion must be implemented as a visible
  design change with tests and mentioned in handoff/status.
- Need commit or destructive git operation: do not do it without PI approval.
- Test failure that cannot be diagnosed within the allotted implementation
  window: report blocker with exact failing test and current hypothesis.
- Ambiguous physics: ask design-discussion, do not invent a fitted constant.

## 15. Acceptance Criteria for Morning Alpha

Minimum acceptable alpha:

- Consolidated plan present.
- P0 schemas/contracts/frame-dump/viz path drafted or briefed to
  implementation-work.
- No PI experimental data used.
- No new hidden numeric constants.
- No gate tolerance edits.
- Any blocked P1 physics is reported clearly rather than forced through.

Strong alpha:

- P0 implemented with tests.
- P1 active-contour state and cortex/area tests implemented with Sanity Gate.
- One synthetic frame dump renders through the stub visualization.
- Implementation handoff records changed files and skipped/deferred modules.

## 16. Bottom Line

Use active contour plus coarse ECM field as the Phase 1 reference engine.
Implement schemas, validation, output, and visualization first. Add executable
physics only when its Sanity Gate has already made the units, signs, timestep,
and measurement protocol explicit. Keep every numeric value in candidate
status until that module earns it.
