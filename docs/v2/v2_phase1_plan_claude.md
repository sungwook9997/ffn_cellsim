# V2 Phase 1 Plan — Claude (independent draft)

> ⚠️ **STATUS: SUPERSEDED DRAFT — 2026-05-06**
> (post-Option-A+B restructure audit cycle, claude-work
> `mcp_msg:2644`, codex approval `mcp_msg:2645`).
>
> This was an independent Project alpha input draft from 2026-05-03.
> Cross-review and synthesis completed; **do not use this file for
> current implementation routing**.
>
> - For alpha synthesis, see
>   `docs/v2/v2_phase1_plan_consolidated.md`.
> - For current post-alpha build state, see
>   `docs/v2/v2_current_build_state.md` (refreshed 2026-05-06,
>   commit `da581cd`).
>
> The 2026-05-03 "WORK IN PROGRESS" status line below is retained
> as historical context and is not the current status.

**Author**: Claude (design-discussion, claude-opus-4-7[1m]).
**Date**: 2026-05-03 23:16 KST.
**Status**: WORK IN PROGRESS — being expanded over ~3 hours.
**Companion**: `docs/v2/v2_phase1_plan_codex.md` (Codex independent draft).
**Authority**: this file is one of two independent Phase 1 drafts produced
overnight per PI id=722, 726, 731, 736 (Project α). The two drafts will be
cross-reviewed and consolidated into `docs/v2/v2_phase1_plan_consolidated.md`
before dispatch to the implementation-work room.

**Constraints baked in (from rounds 1-30 lock + PI overnight directive)**:
- No PI experimental data anchoring (PI id=726). Literature anchors only,
  IF≥15 priority. PI 2026-03-13 spheroid CSV is comparison-only.
- No Magic-Number Block violations (CLAUDE.md hard rule). Every numerical
  constant must be derivable + grid-invariant + not fitted to a target.
- Every physics/numerics module records Sanity Gate Protocol §1-6 before
  first execution.
- Layer 3+ (full ECM fiber network, full spheroid 1000+ cells, reduced
  surrogate) deferred. Phase 1 stops at Tier 2 small cluster (N=10-50).
- Hydrodynamic continuum description belongs to a future Tier 3 reduced
  surrogate engine. Phase 1 carries closure-extraction outputs but does
  not solve continuum equations.
- ECM remodeling is in Phase 1 from day one (PI id=667, 692). Active
  contour without ECM is rejected as biologically empty.
- v2 priority order is single-cell ground truth → small-cluster
  ground truth (PI id=698, "Phase 1까지 진행"). No spheroid-scale runs
  in Phase 1.

---

## §0 Vision and scope

### §0.1 Why Phase 1 exists

Phase 1 is the **mechanobiology reference engine** — a cell-resolved
simulator whose job is to *generate* ground truth that can later be coarse-
grained into reduced surrogates. It is not a production sweep engine, not
a confocal-quality renderer (yet), and not a spheroid-scale hydrodynamic
solver. It is the foundation that all later v2 work builds on.

Phase 1 deliverable, in one sentence: **a single-cell to small-cluster
mechanobiology simulator on a remodelable Col1-like substrate, with
filopodia/lamellipodia events and focal adhesion dynamics resolved per
cell, validated against published cell biology measurements (no PI
experimental fit), with closure-quality output for future Tier 3
surrogate derivation.**

### §0.2 What Phase 1 does NOT do

- No fitting to PI 2026-03-13 CSV (per PI id=726).
- No 1000-cell production spheroid (per PI id=698 Tier 2 cap).
- No reduced active-gel / Marangoni continuum solver. Tier 3.
- No explicit ECM fiber graph. Phase 1 uses a continuous ECM tensor field
  (stiffness, density, orientation) with remodeling memory; individual
  fibers are deferred.
- No biochemical reaction networks (Notch, Wnt, TGF-β, Hippo full
  pathway). YAP/TAZ enters as a single mechanosignal scalar only.
- No O2/nutrient diffusion-reaction. Schema hooks may exist; dynamics
  default OFF until necessary (mini-spheroid N>50).
- No apoptosis-specific event taxonomy beyond `cell_state ∈ {alive,dead}`.
- No confocal-quality Blender renders (Phase 1 ships stub Plotly/Vispy
  viz; Blender pipeline is a separate post-Phase-1 design unit).

### §0.3 Five questions Phase 1 must answer

These are the scientific questions Phase 1's outputs must be capable of
addressing. The architecture below is justified by these questions.

1. **Filopodia statistics under varying ECM stiffness/density**: are
   event rate, angular distribution, and lifetime PDF consistent with
   published live-imaging measurements?
2. **Focal adhesion lifetime distribution**: does the molecular clutch
   model produce FA lifetimes in the literature range
   (~30 s — several minutes) under realistic Col1 ligand density?
3. **Cell-substrate traction balance**: is per-cell traction (radial
   plus tangential) consistent with TFM measurements of single
   epithelial cells on collagen-coated substrate (Plotnikov 2012, etc.)?
4. **Contact inhibition signature in a 5-20 cell cluster**: do follower
   cells show suppressed protrusion at contact edges and reorient
   polarity, with the magnitude consistent with Carmona-Fontaine 2008
   and Stramer & Mayor 2017?
5. **ECM remodeling memory in 80 h spreading**: does the ECM stiffness/
   alignment field accumulate plastic deformation under repeated cell
   traction in a way consistent with Hall 2016 / Eichinger 2020?

A Phase 1 simulator that cannot output the data needed to ask these
five questions has not delivered Phase 1.

### §0.4 Time horizon

- Phase 1 plan finalization (this document + Codex's consolidated): tonight.
- Phase 1 dispatch to implementation-work: tonight after consolidation.
- Phase 1 implementation calendar: 4-8 weeks for full biology, with the
  overnight α scope delivering only schema + active-contour core +
  cortex-tension force + stub viz.

---

## §1 State schema (rounds 1-30 lock + Phase 1 expansion)

### §1.1 Per-cell state — `acs/v2/single_cell.py`

Already-locked fields (rounds 1-30):
- `cell_id: str`
- `time_s: float`
- `measurement_boundary: MeasurementBoundary` (post-migration; see
  `acs/v2/measurement_boundary.py`)
- `height_um: Optional[float]`
- `polarity_xy: Optional[tuple[float, float]]`
- `protrusions: list[ProtrusionEvent]`
- `adhesions: list[FocalAdhesionState]`

New fields locked in design-discussion id=719/720 for Phase 1:
- `cell_state: Literal["alive", "dead"] = "alive"`
- `cell_age_s: float = 0.0`
- `cell_cycle_phase: Optional[Literal["G0","G1","S","G2","M"]] = None`
- `division_count: int = 0`
- `parent_cell_id: Optional[str] = None`
- `mechanosignal_yap_taz: Optional[float] = None` (range [0,1] when set)
- `neighbor_cell_ids: tuple[str, ...] = ()`

Validation rules (added in §1):
- All scalars finite where applicable.
- `cell_age_s ≥ 0`.
- `division_count ≥ 0`.
- `mechanosignal_yap_taz ∈ [0,1]` when not None.
- `neighbor_cell_ids` are unique, non-empty, not equal to self.
- `cell_cycle_phase` ∈ allowed set when not None.
- Dead cells (`cell_state == "dead"`) keep boundary for visualization
  but skip dynamics updates in physics modules; this is asserted by
  the dynamics modules, not by `SingleCellState.validate()`.

### §1.2 Protrusion event — `acs/v2/single_cell.py:ProtrusionEvent`

Already locked. Phase 1 adds no new fields, but uses the event list
heavily.
- `time_s, cell_id, event_type ∈ {"lamellipodium","filopodium","retraction"}`
- `boundary_angle_rad, length_um, lifetime_s, confidence`

Phase 1 dynamics will append events to `SingleCellState.protrusions`
once per simulation timestep (or once per protrusion lifetime, depending
on model). Statistics module aggregates events per cell per frame.

### §1.3 Focal adhesion state — `acs/v2/single_cell.py:FocalAdhesionState`

Already locked. Phase 1 may add a `traction_xy_pN: tuple[float, float]`
field if the molecular clutch model exposes per-FA traction
(this depends on solver choice — see §6). Defer field add to clutch
unit.

### §1.4 ECM substrate state — NEW `acs/v2/ecm_substrate.py`

Phase 1 introduces the ECM substrate as a first-class state, NOT
hidden inside a force law. Reasoning: ECM has its own dynamics
(remodeling), its own observables (alignment, density), and its own
literature (Eichinger 2020, Hall 2016, Trichet 2012). State below
captures it without committing to fiber-level representation.

```python
@dataclass(frozen=True, slots=True)
class ECMSubstrateState:
    time_s: float
    grid_origin_um: tuple[float, float]      # (x0, y0)
    grid_spacing_um: float                   # uniform dx
    grid_shape: tuple[int, int]              # (Nx, Ny)
    stiffness_kPa: np.ndarray                # shape (Nx, Ny), float64
    fiber_density: np.ndarray                # shape (Nx, Ny), float64, normalized [0, 1]
    alignment_tensor: np.ndarray             # shape (Nx, Ny, 2, 2), Q tensor (traceless symmetric)
    accumulated_traction: np.ndarray         # shape (Nx, Ny, 2), memory state for remodeling
    initial_stiffness_kPa: np.ndarray        # baseline reference (for remodeling delta)
    initial_density: np.ndarray
```

Validation:
- All grids same `(Nx, Ny)` shape.
- `stiffness_kPa > 0` everywhere.
- `fiber_density ∈ [0,1]`.
- `alignment_tensor` traceless (sum of diagonals = 0 within float tol)
  and symmetric.
- All grids finite.
- Frame `time_s ≥ 0`.

Defaults at simulation init:
- Stiffness 1-10 kPa for "soft" Col1 coating (Engler 2006 substrate
  range), uniform.
- Density 0.5 (mid-range).
- Alignment tensor Q = 0 (isotropic).
- Accumulated traction = 0.

This grid is a 2D field. 2.5D substrate with finite thickness is a
later refinement. Phase 1 substrate is treated as 2D field of effective
mechanical properties; cell-substrate coupling sees `(stiffness, density,
alignment)` at the cell footprint.

### §1.5 Cell cluster state — NEW `acs/v2/cell_cluster.py`

```python
@dataclass(frozen=True, slots=True)
class CellClusterState:
    time_s: float
    cells: tuple[SingleCellState, ...]
    junctions: tuple[JunctionState, ...]
    contact_graph_csr: tuple[np.ndarray, np.ndarray]  # CSR indptr, indices

    def cell_by_id(self, cell_id: str) -> SingleCellState: ...
    def junction_by_pair(self, a: str, b: str) -> Optional[JunctionState]: ...
```

Validation:
- All cells unique by `cell_id`.
- All junctions reference cells present in the cluster.
- Junction `(a,b)` unordered uniqueness.
- Each cell's `neighbor_cell_ids` matches the contact graph.
- `time_s` consistent across cells (within float tol).

The CSR graph is a redundant cache of `neighbor_cell_ids`; it exists
because cluster-level dynamics (T1 transitions, traction propagation)
need fast neighbor iteration without rebuilding from per-cell tuples.

### §1.6 Junction state — NEW `acs/v2/junction.py`

```python
@dataclass(frozen=True, slots=True)
class JunctionState:
    cell_id_a: str
    cell_id_b: str
    contact_length_um: float
    contact_age_s: float
    e_cadherin_maturity: float        # [0,1] proxy for E-cad density at junction
    contact_inhibition_active: bool   # boolean: is protrusion suppressed across this edge?
    junction_tension_pN_per_um: float # cortex-tension-like line tension at junction
```

Validation:
- Cell ids non-empty, distinct.
- All scalars finite.
- `contact_length_um ≥ 0`.
- `contact_age_s ≥ 0`.
- `e_cadherin_maturity ∈ [0,1]`.
- `junction_tension_pN_per_um` may be negative for short transient, but
  Phase 1 default rejects negative (junctions are tensile, not
  compressive line elements).

### §1.7 Active contour state — NEW `acs/v2/active_contour.py`

This is the simulation-state representation of a single cell's boundary
during dynamics. Distinct from `MeasurementBoundary` which is a
measurement-space object emitted at validation timepoints.

```python
@dataclass(slots=True)
class ActiveContourState:
    cell_id: str
    time_s: float
    vertices_xy_um: np.ndarray            # (N, 2), float64, mutable
    vertex_velocity_xy_um_per_s: np.ndarray  # (N, 2)
    rest_perimeter_um: float
    rest_area_um2: float
    cortex_tension_pN_per_um: float
    bulk_modulus_pN_per_um: float         # area-restoring stiffness
    perimeter_modulus_pN_per_um: float    # perimeter-restoring stiffness
    last_remesh_time_s: float

    def to_measurement_boundary(self, source_modality: str) -> MeasurementBoundary: ...
```

Note: this is mutable (not frozen) because vertex positions evolve every
timestep. Frozen version produced via `to_measurement_boundary()` at
output frames. This is exactly the "simulation state vs measurement
state" split locked in design-discussion round 12-14.

Validation per timestep:
- Vertex array shape (N, 2), N ≥ 3 (after any remeshing).
- All vertices finite.
- Polygon non-self-intersecting (run periodic check, not every timestep).
- `last_remesh_time_s ≤ time_s`.

### §1.8 Cluster-level slow biology hooks

Per-cell hooks already in §1.1 (cell cycle, mechanosignal). Cluster-
level slow biology adds:
- Diffusion-reaction fields (O2, nutrients) — DEFERRED to mini-spheroid
  entry. Schema not added in Phase 1.
- Morphogen gradients — DEFERRED entirely.
- Inflammatory/immune state — DEFERRED entirely.

---

## §2 Solver architecture and boundary representation

This is the deepest design decision in Phase 1. The choice of how a cell's
shape evolves dictates what biology can be expressed, what computational
cost we pay, and what existing tooling we can or cannot reuse.

### §2.1 Candidate representations

Below are the seven serious candidates, evaluated against Phase 1's five
scientific questions (§0.3) and PI's emphasis on filopodia (id=692
"filopodia most important in my opinion").

#### §2.1.1 Active contour / parametric polygon (Kass-Witkin-Terzopoulos 1988; modern cellular variants)

Cell boundary is an ordered polygon `(v_0, v_1, ..., v_{N-1})` evolved
under external (substrate, neighbor, ECM) and internal (cortex tension,
area pressure) forces.

- *Pros*: Direct polygon = trivial mapping to `MeasurementBoundary`; no
  reprojection. FA naturally lives at a vertex or edge. Filopodia
  representable as a high-resolution local segment with dynamic vertex
  insertion. Computational cost scales as O(N_cells × N_vertices).
  Compatible with explicit-time integration. Easy to attach per-edge
  active stress for protrusion modeling.
- *Cons*: Topology change (cell division, filopodia spawning/retraction,
  cell-cell contact merging) requires explicit algorithms (remeshing,
  vertex insertion, cytokinesis split). Self-intersection prevention
  needed under large deformations. Vertex resolution determines the
  smallest feature: capturing 100-nm-wide filopodia means vertex
  spacing ≤ ~50 nm at the protrusion tip, which is achievable only
  with adaptive resolution (uniform high resolution everywhere is
  wasteful).
- *Maturity in cell biology context*: parametric active contour is
  classical (Kass 1988 in computer vision, then adopted for cell
  segmentation). Less common as a physics simulation tool because of
  topology limits, but emerging — e.g. Wang et al. 2019 "Active contour
  model with internal pressure" for cell mechanics, Maître 2012 cortex
  tension models implicitly use a polygon representation.

#### §2.1.2 Phase field (Cahn-Hilliard / Allen-Cahn for cell shape)

Cell represented by a smooth scalar field `φ(x,t) ∈ [-1, 1]` with a
diffuse interface of finite width `ε`. Cell occupies the region
`φ > 0`. Time evolution by Ginzburg-Landau type free-energy minimization.

- *Pros*: Topology change is automatic — cells can split, merge, fuse,
  develop holes without algorithm changes. Smooth measurement (boundary
  is a level set, no discrete vertex artifacts). Well-suited for
  multicellular tumor growth (Travasso 2011, Wise 2008). Conservative
  area / mass tracking trivial (volume integral of φ).
- *Cons*: Eulerian grid required → memory and compute scale as
  `O(grid_volume)` regardless of cell count. Interface width `ε` sets
  minimum length scale: `ε ~ 200 nm` at best on a `dx = 50 nm` grid →
  4× the size of a filopodium tip → filopodia under-resolved. Attaching
  point-wise events (FA, traction concentration) requires extra
  machinery (Lagrangian markers on top of Eulerian field). Active
  protrusion forces are spatially smeared by the diffuse interface.
- *Maturity*: well-developed for tissue-scale problems (Mehboudi 2018,
  Najem 2014), less for single-cell filopodia.

#### §2.1.3 Vertex model (Honda 1978; Fletcher 2014; Bi-Manning 2015)

Cells tessellate space; each vertex is shared by 3+ cells. Energy
function `E = Σ (K_A (A - A_0)^2 + Γ L^2 + ...)`. Topology changes by
T1 transitions (edge swap) and T2 (cell removal).

- *Pros*: Excellent for confluent epithelial monolayers. T1 transitions
  built-in. Famous for jamming/unjamming transition (Bi 2015, Park
  2015). Computational cost low.
- *Cons*: Assumes confluence — every vertex shared. Single isolated
  cell doesn't fit (no neighbors to share vertices with). Free
  boundaries (cell-substrate contact line, lamellipodia at the
  spreading edge) need extension. Filopodia not naturally
  representable.
- *Verdict for Phase 1*: rejected for Tier 1 (single cell). Possibly
  re-examined at Tier 2+ where confluence emerges, but Phase 1 cells
  are not confluent — they are isolated single cells or small loose
  clusters spreading on a Col1 substrate.

#### §2.1.4 Cellular Potts model (Graner-Glazier 1992)

Cells are connected sets of pixels on a lattice, each cell with a label.
Hamiltonian `H = Σ J(τ_i, τ_j) (1 - δ(σ_i, σ_j)) + λ_V (V - V_0)^2 +
...`. Time evolution via Metropolis MCMC pixel-flip sampling.

- *Pros*: Topology automatic. Well-developed (CompuCell3D framework).
  Multi-cell easy. Pixel-level shape detail.
- *Cons*: MCMC time, not physical time — mapping `MCS → seconds`
  requires per-system calibration (a Magic-Number Block flag). No
  inertia, no momentum: dynamics is gradient-descent-like, not
  Newtonian/Langevin. Filopodia need fine lattice (sub-pixel filopodia
  not possible). Not amenable to GPU vectorization at high quality
  (sequential pixel proposals).
- *Maturity*: very mature (Cellular Potts 1992, CompuCell3D
  long-running). But the time-mapping problem is a known fundamental
  limitation when comparing to live-imaging.

#### §2.1.5 PhysiCell (Ghaffarizadeh 2018)

Hybrid agent-based: each cell is a sphere (or short list of spheres
for ellipsoidal shapes) with center, radius, force-resolved motion.
Substrate / chemical fields are Eulerian on a separate grid.
- *Pros*: Production-grade, validated for tumor growth. Handles
  multicellular + biotransport. C++ with Python bindings.
- *Cons*: Cells are spheres — no shape dynamics, no filopodia, no
  active contour. Cell-cell contacts via force laws between centers,
  not via explicit junction edges. Cannot answer Phase 1's filopodia
  question §0.3.1.
- *Verdict*: rejected for Phase 1. Could be considered for a future
  large-scale spheroid run that doesn't need filopodia resolution.

#### §2.1.6 CompuCell3D (Swat 2012)

Cellular Potts on top of an integrated framework with biology
modules.
- *Pros*: large user base, many published examples, GUI-driven.
- *Cons*: same fundamental limits as Cellular Potts (time mapping,
  filopodia under-resolved). Framework integration is not trivial for
  custom physics.
- *Verdict*: rejected for Phase 1 same reasons as Cellular Potts.

#### §2.1.7 Chaste (Mirams 2013)

Modular C++ framework supporting vertex, cell-centre, Cellular Potts,
and immersed-boundary representations. Strong validation, peer-
reviewed.
- *Pros*: Multiple representations available, well-tested.
- *Cons*: Heavy C++ codebase; integration with our Python/Taichi stack
  non-trivial; we end up using Chaste's representations rather than
  designing our own — losing the v2 design rationale.
- *Verdict*: rejected as a primary engine. Could be a future
  cross-validation reference (re-implement Phase 1 in Chaste and
  compare).

### §2.2 Recommendation for Phase 1: Active contour (parametric polygon) with adaptive vertex density

Active contour wins for Phase 1 because:

1. **Direct boundary representation**: simulation state (`ActiveContourState
   .vertices_xy_um`) and measurement state (`MeasurementBoundary
   .vertices_xy_um`) share representation modulo the canonical-CCW
   guard. Zero-cost Rule 11 compliance — the simulation IS the
   measurement, no projection ambiguity.
2. **Filopodia native**: a filopodium is an extruded high-resolution
   patch of the boundary. Adaptive vertex density (insert vertices
   when local curvature exceeds a threshold or a protrusion event
   spawns) captures filopodia natively without changing representation.
3. **FA / per-vertex events natural**: focal adhesions live at vertices
   or edges. Per-FA traction maps directly to a per-vertex force.
4. **Cluster manageable**: at N=10-50 cells, contact detection is
   O(N² × N_vert) which is tractable on the Laptop A5000.
5. **Cell division is rare and discrete**: cytokinesis is a once-per-
   ~24-hours event, handled as an atomic algorithm (mitotic furrow
   contracts → boundary splits into two), not continuous topology
   change. Acceptable.
6. **Newtonian/Langevin physical time**: vertex update is
   `dv/dt = (F_total)/m - γv + noise`. Time is real seconds. No MCMC
   pseudo-time (which Cellular Potts would impose).
7. **Compatible with our stack**: Python + Taichi for vertex updates;
   no need to introduce new heavy frameworks.

### §2.3 Active contour numerical scheme

The vertex evolution scheme matters. Three serious candidates:

#### §2.3.1 Overdamped Newton (recommended)

`γ * dv_i/dt = F_i^total + ξ_i(t)` where `γ` is per-vertex friction
proportional to local substrate drag, `F_i^total = F_cortex + F_pressure
+ F_FA + F_protrusion + F_neighbor + F_ecm + F_thermal`, and `ξ_i(t)`
is white noise from the fluctuation-dissipation theorem (or a
Langevin-style stochastic forcing if needed).

- Time integration: explicit forward Euler at fine `dt` (say 1-10 ms)
  or symplectic if higher accuracy needed.
- Stability: CFL bound `dt < γ / k_max` where `k_max` is the largest
  per-vertex restoring stiffness (typically cortex tension at small
  curvature).
- Sanity Gate §1: dimensional check OK, dt vs CFL OK at proposed
  parameters (need to compute).

#### §2.3.2 Mean curvature flow (for cortex tension only)

`dv_i/dt = -κ * n̂_i * γ_cortex / γ_drag` — this is the geometric flow
that minimizes perimeter. Simple, but:
- Hides physics inside a geometric step, breaking the sign/sense check
  (Sanity Gate §5) at Mode 5 which wants every force term explicit.
- Cannot couple FA traction or substrate forces directly.
- Useful as a verification sub-component (curvature term test in
  isolation), not as the production solver.

#### §2.3.3 Snake-style elastic + image force (computer vision origin)

`E[v] = α * E_elastic + β * E_curvature + γ * E_external`. Optimization-
based, not Newtonian. Inherits computer-vision idioms (image gradient
as external force) that don't map cleanly to cell biology.
- Rejected for Phase 1: not physically motivated, hard to add stochastic
  events.

### §2.4 Adaptive vertex density (filopodia resolution)

Critical for §0.3 question 1. Strategy:

1. **Baseline vertex spacing**: 0.5-1 μm along smooth boundary.
2. **Local refinement on curvature**: when local curvature `κ > κ_max
   = 1/(0.2 μm)`, insert intermediate vertices to maintain
   `vertex_spacing ≤ 0.1 μm`.
3. **Filopodium spawn**: when a `ProtrusionEvent` of type "filopodium"
   is generated at boundary angle `θ`, insert 5-10 vertices in a
   narrow patch (~0.5 μm wide) centered at `θ` to allow the filopodium
   to extrude with sub-micron tip width.
4. **Filopodium retraction**: when the event lifetime expires,
   collapse the patch back to baseline density via vertex removal.
5. **Cost containment**: total vertices per cell capped at
   `N_max ≈ 1000`. If exceeded, the lowest-curvature segments are
   re-coarsened first.

Cost estimate: baseline 100 vertices per cell (~30 μm perimeter at
0.3 μm spacing), peak ~300 vertices per cell with 5 active filopodia.
Per-cell cost ~0.3 ms per timestep on Laptop A5000 → 50-cell cluster
~15 ms per timestep → 80h sim at 10 ms timestep = 28.8M steps × 15 ms
= 432,000 s = ~5 days. Too slow. Mitigations:

- **Adaptive timestep**: ECM remodeling is slow (minutes), can use
  larger `dt_ECM = 1 s` with operator splitting.
- **Filopodia event-driven**: filopodia don't update every 10 ms;
  they're discrete events with stochastic spawn/retract, and the
  vertex update during the event is at adaptive `dt_filopodium ≤
  1 ms` only locally.
- **GPU vectorization**: vertex updates parallel across cells via
  Taichi.
- **Phase 1 acceptance**: production 80h sim is NOT a Phase 1
  deliverable. Phase 1 deliverable is (a) single-cell up to ~10 min
  spreading, (b) small cluster up to ~30 min spreading. Both
  tractable on Laptop A5000 in under an hour wall-clock.

### §2.5 Topology operations (cytokinesis, contact merge)

Phase 1 needs two topology operations:

1. **Cytokinesis**: when `cell_cycle_phase` advances to "M" and the
   division trigger fires, the active contour splits into two
   daughter cells. Algorithm: identify the polarity-perpendicular
   axis; find the two vertex chains on either side of the axis;
   close each chain with a new edge; assign daughter `cell_id`s,
   inherit parent `cell_id` to `parent_cell_id`. Both daughters
   carry the parent's mechanosignal scalar initially. Daughter rest
   areas = parent rest area / 2.

2. **Cell-cell contact**: when two cells' boundaries come within
   `contact_threshold_um = 0.5 μm`, a contact edge is created. The
   contact is represented as a `JunctionState` linking the cell ids,
   not as a topological merge of polygons. Each cell keeps its own
   polygon; the junction is an explicit cross-cell edge that exerts
   adhesion force on both cells' nearby vertices.

The cell-cell contact is intentionally NOT a polygon merge — that
would force vertex model topology, which we rejected. Two polygons
with an explicit junction line is the right level of detail for
Tier 2.

### §2.6 Verdict

Phase 1 boundary representation: **active contour with adaptive vertex
density, overdamped Langevin dynamics, explicit cytokinesis, junction-
mediated cell-cell contact**.

Implementation lives in `acs/v2/active_contour.py` (state) and
`acs/v2/dynamics/active_contour.py` (numerical evolution).

---

## §3 Force terms and per-term literature anchors

Each term below: physical form, sign, parameter literature anchor,
Magic-Number Block check, Sanity Gate sense check. All terms are
per-vertex forces summed into `F_i^total` for the overdamped Newton
update of §2.3.1.

### §3.1 Cortex tension (line tension)

`F_i^cortex = -∂E_perimeter/∂v_i` where
`E_perimeter = γ_cortex * Σ_e |e_length|`.

Equivalent per-vertex form: each vertex `i` between edges `(i-1,i)`
and `(i,i+1)` feels a force pulling toward the bisector of those
two edges, with magnitude `γ_cortex * (1/|e_left| - 1/|e_right|)`
along the perpendicular component. Net effect: shortens perimeter,
smooths boundary, suppresses high curvature.

- *Sign*: contractile (negative work on expansion). Sanity Gate §5: ✓.
- *Literature anchor*: Salbreux et al. 2012 *Trends in Cell Biology*
  ("Actin cortex mechanics and cellular morphogenesis") — review-
  level, IF~20. Cortex tension for adherent epithelial cells:
  `γ_cortex ~ 0.1 - 1.0 nN/μm` (= 100 - 1000 pN/μm).
  Maître et al. 2012 *Cell* ("Adhesion functions in cell sorting by
  mechanically coupling the cortices of adhering cells") — IF~40.
  Bambardekar et al. 2015 *PNAS* ("Direct laser manipulation reveals
  the mechanics of cell contacts in vivo") — IF~10, gives in vivo
  cortex tension `~0.5 nN/μm` for Drosophila.
- *Phase 1 default*: `γ_cortex = 0.3 nN/μm` (mid-range MCF7-like
  epithelial cell — note: MCF7-specific cortex tension is not pinned
  in literature; we use the epithelial range. This is a literature-
  derived range, not a fitted value).
- *Magic-Number Block*: literature-derived ✓; grid-invariant ✓
  (not coupled to vertex spacing); not fitted to any target ✓. PASS.

### §3.2 Area-restoring pressure (cytoplasmic / volume regulation)

`F_i^pressure = -∂E_area/∂v_i` where
`E_area = (K_A / 2) * (A - A_0)^2`, `A` = current polygon area.

Equivalent per-vertex form: each vertex feels a force in the outward
normal direction, magnitude `K_A * (A - A_0) * |e_avg|/2` where
`|e_avg|` is the average edge length at that vertex. When `A < A_0`
the cell is "compressed" → pressure pushes outward; when `A > A_0`
it pulls inward.

- *Sign*: outward when compressed (positive work on expansion when
  underfilled), inward when over-expanded. Sanity Gate §5: ✓.
- *Literature anchor*: Bi et al. 2014 *PRX* ("Energy barriers and
  cell migration in densely packed tissues") — IF~12 (PRX), with
  vertex-model area modulus. Farhadifar et al. 2007 *Current Biology*
  ("The influence of cell mechanics, cell-cell interactions, and
  proliferation on epithelial packing") — IF~10 (CB), foundational
  vertex model fit.
- *Phase 1 default*: `K_A = 100 pN/(μm³)` so that
  `K_A * A_0^2 ~ K_A * 1000 μm⁴ = 1e5 pN·μm = 0.1 nN·μm` (energy
  scale comparable to cortex tension at typical cell perimeter).
  Actually in 2D area units: `K_A` has units `pN/μm³` such that
  `(1/2) K_A (A - A_0)^2` has units of energy `pN·μm`. For
  `A_0 = 800 μm²` and a 1% area perturbation
  `δA = 8 μm²`, energy `= 0.5 * 100 * 64 = 3200 pN·μm = 3.2 nN·μm`,
  comparable to cortex energy scale.
- *Magic-Number Block*: somewhat free here; the published vertex-
  model fits give `K_A * A_0^2 / Γ_perimeter ~ O(1)` dimensionless
  ratio (Bi 2015 ratio for jamming transition). We tie our `K_A` to
  this ratio rather than picking arbitrarily. PASS provided this
  ratio is documented in code.

### §3.3 Focal adhesion traction (substrate-coupled force)

When a vertex `i` lies inside the substrate footprint of an active
focal adhesion `FA_a` with maturity `m_a` and bound fraction `b_a`,
that vertex feels a traction force pulling it toward the FA's
substrate anchor point.

`F_i^FA = - K_FA * m_a * b_a * (v_i - r_a^substrate)` for vertices
within `R_FA = 0.5 μm` of the FA position `r_a^substrate`.

The FA's reaction force is applied to the substrate (Newton's third
law) at `r_a^substrate`, contributing to the per-grid-cell
`accumulated_traction` of `ECMSubstrateState`.

- *Sign*: attractive between cell and substrate when FA is engaged.
  Sanity Gate §5: ✓.
- *Literature anchor*:
  - Plotnikov et al. 2012 *Cell* ("Force fluctuations within focal
    adhesions mediate ECM-rigidity sensing to guide directed cell
    migration") — IF~40. Per-FA traction in mature adhesions:
    `~30-100 nN per FA` for fibroblasts on 30 kPa substrate.
  - Kanchanawong et al. 2010 *Nature* ("Nanoscale architecture of
    integrin-based cell adhesions") — IF~50. FA structural layers,
    FA size ~ 0.5-2 μm.
  - Schwarz & Gardel 2012 *J Cell Sci* (review IF~5) on traction
    force microscopy of adherent cells.
- *Phase 1 default*: `K_FA = 0.1 nN/μm` (gives ~0.05 nN per μm of
  vertex displacement; for `m_a=1, b_a=1` and 0.5 μm displacement,
  force `=0.025 nN`, consistent with Plotnikov 30-100 nN total
  divided across many FAs).
- *Magic-Number Block*: literature-derived ✓; grid-invariant
  (per-FA, not per-vertex-spacing) ✓; not fitted ✓. PASS.

### §3.4 Lamellipodia protrusion force

A `ProtrusionEvent` of type "lamellipodium" centered at boundary
angle `θ` and lasting `lifetime_s` exerts an outward force on
vertices within an angular window `Δθ_lamellipodium = 30°` of `θ`.

`F_i^lamellipodium = F_lam * n̂_i * shape(θ_i - θ, Δθ)` where
`F_lam ~ 0.1-1 nN per active lamellipodium` and `shape(·)` is a
cosine profile peaked at `θ`.

- *Sign*: outward (extends boundary). Sanity Gate §5: ✓.
- *Literature anchor*:
  - Pollard & Borisy 2003 *Cell* ("Cellular motility driven by
    assembly and disassembly of actin filaments") — IF~50. Reviews
    actin polymerization-driven protrusion physics.
  - Krause & Gautreau 2014 *Nat Rev Mol Cell Biol* ("Steering cell
    migration: lamellipodium dynamics and the regulation of
    directional persistence") — IF~70. Lamellipodium lifetime
    seconds-minutes, force per lamellipodium ~ pN-nN range.
  - Mogilner & Oster 2003 *Curr Biol* on actin polymerization
    force generation (~ pN per filament, ~10² filaments per
    lamellipodium → ~ nN total).
- *Phase 1 default*: `F_lam = 0.5 nN per lamellipodium`,
  `Δθ_lamellipodium = 30°`, lifetime drawn from exponential PDF
  with mean 2-5 minutes (Krause-Gautreau review).
- *Magic-Number Block*: lifetime from literature ✓; force literature
  range ✓; angular window is a model choice (not literature-pinned)
  — flag as exploratory, document in `docs/12_validation.md`.

### §3.5 Filopodia tip force

A `ProtrusionEvent` of type "filopodium" creates a slender extruded
patch on the boundary (per §2.4 adaptive density), with a tip vertex
that experiences a strong outward force along the filopodium axis.

`F_tip^filopodium = F_filo * n̂_filo * f_age(t_event)` where `f_age`
is a piecewise function: fast outward push during growth phase
(~10s), neutral during plateau, retraction-pull during retraction.

- *Sign*: outward during extension, inward during retraction. Sanity
  Gate §5: ✓ (sign reverses with age, must be tested).
- *Literature anchor*:
  - Mattila & Lappalainen 2008 *Nat Rev Mol Cell Biol* ("Filopodia:
    molecular architecture and cellular functions") — IF~70.
    Filopodium length 1-10 μm, lifetime seconds to ~10 min, growth
    velocity ~50-200 nm/s.
  - Bornschlögl et al. 2013 *PNAS* ("Filopodial retraction force is
    generated by cortical actin dynamics and controlled by reversible
    tethering at the tip") — IF~10. Direct measurement: filopodium
    pulling force `~5-20 pN`, retraction velocity `~30-100 nm/s`.
- *Phase 1 default*: `F_filo = 10 pN` during growth, `-15 pN` during
  retraction (slightly stronger to enable retraction against cortex
  tension), growth lifetime from exponential PDF with mean 30 s.
- *Magic-Number Block*: forces literature-derived ✓; lifetimes
  literature ✓; the asymmetry (15 vs 10 pN) is a literature
  observation (Bornschlögl).

### §3.6 Cell-cell junction adhesion force

When a `JunctionState` exists between cells `a` and `b`, vertices on
`a`'s boundary near the junction edge feel an attractive force toward
the corresponding edge on `b` and vice versa.

`F_i^junction = γ_junc * m_ecad * (v_corresp_b - v_i) * f(d)` where
`d = |v_i - v_corresp_b|` and `f` falls off beyond contact range
`d_max = 1 μm`.

- *Sign*: attractive (cells pulled together). Sanity Gate §5: ✓.
- *Literature anchor*:
  - Maître et al. 2012 *Cell* (already cited) — cortex coupling at
    cell-cell contacts.
  - Bambardekar et al. 2015 *PNAS* — junction line tension
    `~0.1-0.5 nN/μm`.
  - Cai et al. 2014 *Cell* ("Mechanical feedback through E-cadherin
    promotes direction sensing during collective cell migration") —
    IF~40. E-cadherin maturity vs contact age relationship.
- *Phase 1 default*: `γ_junc = 0.2 nN/μm`, `m_ecad ∈ [0,1]`
  initialized to 0 at contact onset, increasing on a maturation
  timescale `τ_ecad = 60 s` (Cai 2014 Fig 4).
- *Magic-Number Block*: γ_junc and τ_ecad literature ✓.

### §3.7 ECM-mediated motility bias and substrate drag

Two effects:
1. *Substrate drag* on each vertex: `F_i^drag = -γ_drag(stiffness_local)
   * v_i`. Drag depends on local ECM stiffness (cell-substrate friction
   scales with substrate compliance and FA density).
2. *ECM alignment-guided protrusion bias*: when the ECM alignment
   tensor at the cell footprint has a non-zero principal direction
   `n̂_ecm`, lamellipodia spawning is biased to occur preferentially
   along ±`n̂_ecm` (contact guidance).

`γ_drag = γ_drag^0 * (E_local / E_ref)^α` with `α ≈ 0.5-1`.
`E_ref = 5 kPa` (typical Col1 coating stiffness).

- *Sign*: drag opposes velocity ✓; alignment biases probability
  distribution, not a force directly.
- *Literature anchor*:
  - Trichet et al. 2012 *PNAS* ("Evidence of a large-scale mechano-
    sensing mechanism for cellular adaptation to substrate stiffness")
    — IF~10. Drag scales with stiffness in adherent cells.
  - Engler et al. 2006 *Cell* ("Matrix elasticity directs stem cell
    lineage specification") — IF~50. Substrate stiffness range and
    cellular response.
  - Hall et al. 2016 *PNAS* ("Fibrous nonlinear elasticity enables
    positive mechanical feedback between cells and ECMs") — IF~10.
    ECM remodeling and contact guidance.
- *Phase 1 default*: `γ_drag^0 = 0.05 nN·s/μm²`, `α = 0.7`.
- *Magic-Number Block*: literature anchored ✓; α exponent is fitted
  to literature (Trichet) — note as a "literature-fitted exponent",
  documented.

### §3.8 Thermal / Langevin noise

`ξ_i(t)` per-vertex stochastic force with `<ξ_i> = 0` and
`<ξ_i(t) ξ_j(t')> = 2 γ_drag k_B T δ_ij δ(t-t')` (fluctuation-
dissipation theorem). Optional intracellular athermal noise can be
added with a separate amplitude.

- *Sign*: random, zero-mean.
- *Literature anchor*: foundational (Langevin 1908; modern
  application to biological cells: Trepat et al. 2009 *Nature
  Physics* "Physical forces during collective cell migration",
  IF~20). Single-cell membrane fluctuations at thermal scale
  (Brochard-Wyart 1976 + many follow-ups).
- *Phase 1 default*: `T = 310 K`, athermal noise OFF for now (PI
  experimental data not used for fitting; thermal-only is the
  honest baseline).
- *Magic-Number Block*: T physical, fluctuation-dissipation derives
  amplitude from drag coefficient — no free parameter. PASS.

### §3.9 Force balance summary

`γ_drag * dv_i/dt = F_cortex + F_pressure + F_FA + F_lam + F_filo +
F_junction + F_ecm_drag + ξ_i`.

Total force on a single quiescent cell at rest = 0 (perimeter
relaxed, area at A_0, no FA engagement, no protrusions). Test:
solver should converge to this state from a perturbed initial
condition (Sanity Gate §3 conservation/relaxation test).

Total momentum is NOT conserved at the cell level because the
substrate absorbs FA traction (open system). But local action-
reaction at FA (cell-substrate) IS conserved: per-FA force on cell
= -force on substrate. Sanity Gate §3: must verify this
explicitly in the FA module.

---

## §4 ECM substrate model

PI's id=667 critique was central: an active contour without ECM
remodeling is biologically empty for 80h Col1 spreading. This
section specifies the minimal-but-honest ECM model that satisfies
that critique without inflating Phase 1 to a full fiber-network
research project.

### §4.1 Three representations evaluated

#### §4.1.1 Static elastic foundation (rejected per PI id=667)
Substrate as constant `E(x) = E_0` everywhere, no dynamics. Cells feel
stiffness but cannot remodel. Rejected: doesn't capture the first-
order phenomenon PI identified (mechanical memory, fiber alignment by
traction, contact guidance).

#### §4.1.2 Continuous evolving tensor field (recommended)
Substrate as a 2D field with three observables:
- Scalar stiffness `E(x,y,t)` (kPa)
- Scalar fiber density `ρ_ecm(x,y,t)` (normalized [0,1])
- Traceless symmetric alignment tensor `Q(x,y,t) ∈ R^{2×2}` whose
  principal axis encodes fiber orientation and whose eigenvalue
  magnitude encodes alignment strength.

Each evolves under cell traction via local update rules (§4.3).
Computational cost: `O(N_grid)` per timestep, where `N_grid` for a
100×100 μm domain at 0.5 μm grid spacing is ~40,000 cells — trivial.

#### §4.1.3 Discrete fiber graph (deferred to Tier 3)
Each fiber as a discrete spring connecting nodes, network
remodeling via fiber breakage and re-attachment. Far more
expensive (O(N_fibers²) for force calculation), and a research
project in itself (Eichinger 2020 Soft Matter). Phase 1 keeps
schema-level hooks (we can recover fiber-graph mode later by
replacing the field-update functions; see §4.5 portability).

### §4.2 Field representation details

ECM grid: uniform 2D, spacing `dx_ecm = 0.5 μm` (~10× larger than
filopodium tip — sufficient for substrate-scale phenomena, not for
sub-cellular ECM detail). Domain matches simulation domain (~100 μm
square for single cell, ~300 μm square for cluster).

`stiffness_kPa[i,j]`: float64, default 5 kPa (Engler 2006 substrate
range for soft tissue mimicking).

`fiber_density[i,j]`: float64 in [0,1], default 0.5 (mid-range).

`alignment_tensor[i,j,:,:]`: 2x2, traceless symmetric. Stored as
`(Q_xx, Q_xy)` with `Q_yy = -Q_xx` — saves memory/storage.

`accumulated_traction[i,j,:]`: 2D vector, accumulator for past cell
traction at this point. Used as memory state for plastic remodeling.
Reset to zero only on explicit substrate refresh (e.g., new
experiment instance).

`initial_stiffness_kPa[i,j]`, `initial_density[i,j]`: baseline
references. `stiffness_kPa - initial_stiffness_kPa` gives the
"remodeled delta", useful for diagnostic visualization (where has
the cell densified the matrix?).

### §4.3 Remodeling dynamics (Phase 1 minimal rule)

ECM evolves under three rules, applied at each ECM timestep
`dt_ecm` (which can be much larger than the cell timestep — the
operator-splitting strategy from §2.4).

#### §4.3.1 Stiffness change under sustained traction (plastic stiffening)

Where cells exert sustained traction, the local matrix densifies and
stiffens (collagen crosslinking, fiber recruitment).

`d E_local / dt = k_stiffen * |⟨T_local⟩|^β * (E_max - E_local)`

where `⟨T_local⟩` is the time-averaged traction at this grid cell
over the recent past (averaging window `τ_ecm_window = 600 s`,
chosen as ~10× FA lifetime so transient FAs don't dominate);
`β ≈ 1` (linear in traction); `E_max = 30 kPa` (Engler 2006 upper
range for cell-densified collagen); `k_stiffen ≈ 0.001 / (s·nN/μm²)`.

- *Sign*: positive (traction stiffens, never softens in Phase 1).
- *Sanity Gate §3*: bounded by `E_max`, no runaway.
- *Literature*: Hall et al. 2016 *PNAS* — fibrous nonlinear elasticity
  positive feedback. Eichinger 2020 *Soft Matter* — review of
  computational ECM remodeling.
- *Magic-Number Block*: `k_stiffen` chosen so that ~600 s of traction
  at 1 nN/μm² gives a stiffness increase of ~3 kPa, matching Hall
  Fig. 4 timescale. PASS (literature-anchored timescale).

#### §4.3.2 Fiber alignment under traction (contact guidance precursor)

The alignment tensor evolves toward the principal direction of the
local traction tensor, with a relaxation timescale:

`d Q / dt = k_align * (Q_target - Q) - λ_isotropic * Q`

where `Q_target` is the traceless symmetric part of `(T_local
T_local^T) / |T_local|^2` (a tensor pointing along the traction
direction); `k_align ≈ 0.005 / s` (alignment timescale ~200 s,
Hall 2016); `λ_isotropic ≈ 0.0005 / s` (slow return to isotropic
when traction relaxes — much slower than alignment, so once aligned
the matrix "remembers" for thousands of seconds).

- *Sign*: alignment grows in direction of traction; isotropic
  relaxation is slow.
- *Sanity Gate §3*: `Q` bounded by `Q_target` magnitude (≤ 1).
- *Literature*: Hall 2016 PNAS, Trichet 2012 PNAS contact guidance.
- *Magic-Number Block*: timescales literature-anchored; ratio
  `k_align / λ_isotropic = 10` chosen so memory persists ~10× longer
  than alignment → captures plastic memory phenomenology.

#### §4.3.3 Fiber density change

For Phase 1, density change is OPTIONAL and OFF by default. The
baseline assumption is that the Col1 coating density doesn't change
substantially in 80h (no MMP-driven degradation, no cell-driven
deposition). Hook is in the schema (§1.4 `fiber_density`) but
dynamics not implemented in Phase 1.

When activated (Phase 2+):
`d ρ / dt = k_deposit * (ρ_max - ρ) * indicator(cell_present) -
            k_degrade * ρ * indicator(MMP_active)`

with `ρ_max = 1`, deposition by cell-secreted ECM, degradation by
MMP. Both rates are cell-type-specific and need PI input or
literature deep-dive.

### §4.4 Cell-ECM coupling (per-vertex force feedback)

When a cell vertex `i` lies above an ECM grid cell at `(I,J)`:

1. **Vertex feels substrate stiffness via `γ_drag`** (§3.7):
   `γ_drag(v_i) = γ_drag^0 * (E[I,J] / E_ref)^α`.

2. **Vertex feels alignment-bias for protrusion** (§5):
   when a `ProtrusionEvent` is being spawned at boundary angle
   `θ_i`, the spawn probability is multiplied by
   `(1 + ε_align * cos(2(θ_i - θ_align[I,J])))` where `θ_align`
   is the principal axis of `Q[I,J]` and `ε_align ≈ 0.5`.

3. **FA traction is deposited into ECM at FA position**:
   on each cell timestep, every active FA `a` adds
   `T_a` to `accumulated_traction[I_a, J_a]`. The ECM update
   rules (§4.3) consume the time-averaged value.

4. **Cell traction can feel "stiffness gradient"** for durotaxis
   (Phase 1 OPTIONAL):
   `F_i^durotaxis = K_duro * ∇E[I,J] * indicator(FA active at i)`.

   - Sign: cells migrate toward stiffer regions (Lo 2000, Engler 2006).
   - Default Phase 1: OFF (durotaxis dynamics not central to spreading
     on uniform substrate). Hook in schema.

### §4.5 Portability to Tier 3 fiber graph

The field representation can be replaced by a fiber-graph
representation (Eichinger 2020) without changing the cell-side API:
the cell only queries `(stiffness, density, alignment)` at its
footprint and deposits traction at FA positions. A fiber-graph
implementation of `ECMSubstrateState` would expose the same query
interface but compute these from fiber-network state. This
deferred-extension contract is the value of the field abstraction.

### §4.6 Sanity Gate for ECM module

- §1 Dimensional: stiffness in kPa, traction in nN/μm², position in μm,
  time in s. **Per Rule 10 (dimensional comparison verification): the
  exact conversion chain between kPa stiffness, FA-footprint area, and
  per-FA traction in nN must be derived in the `ecm_substrate` module's
  Sanity Gate document before first execution.** It is NOT done in this
  plan to avoid baking an unchecked conversion into the contract. Order-
  of-magnitude check: Plotnikov 2012 reports per-FA traction ~30 nN on
  ~30 kPa substrate, consistent with our defaults; the precise unit chain
  is deferred to the module-level Sanity Gate.
- §2 Boundary cases: `E → 0` would divide by zero in `γ_drag`;
  guard against `E_local < E_min = 0.1 kPa`. `E → ∞` produces
  diverging drag; cap at `E_max`. `Q → infinite eigenvalue` — bound
  by construction (`|Q| ≤ 1`).
- §3 Conservation: stiffness/density bounded by [0, max]; alignment
  tensor traceless preserved by update rule (verify in implementation).
- §4 Numerical: ECM `dt_ecm = 1 s`, much larger than cell `dt_cell =
  10 ms`. Operator splitting: every 100 cell steps, run 1 ECM step.
- §5 Sign: stiffening is always positive ✓; alignment relaxation
  toward isotropic is always negative on |Q| ✓.
- §6 Measurement protocol: ECM state is exposed as a separate field in
  the per-frame HDF5; downstream metrics that use ECM state (e.g.
  "fraction of cell footprint over densified matrix") must declare
  `required_artifacts = (BOUNDARY_CONTOURS, ECM_FIELD)` — the latter
  is a new `ArtifactKind` we add (see §13 sweep).

---

## §5 Protrusion / filopodia event model

PI's id=692 explicitly elevated filopodia ("most important in my
opinion"). This section gives the event-driven stochastic model.
The model is biologically literature-anchored (no fitting), but
the event-spawn-rate parameters are literature ranges, not fitted
to PI data.

### §5.1 Event-driven not continuous

The actin biochemistry that drives lamellipodia/filopodia operates
on millisecond-to-second timescales with strong stochasticity.
Modeling each actin filament is impractical at Phase 1 scale.
Instead, we use an **event-rate / hazard** framework:

- Each cell has spatiotemporally varying "activity field"
  `λ(θ, t, cell_state, ecm_state)` giving the spawn rate (events
  per second per radian) for new lamellipodia/filopodia along the
  boundary at angle `θ`.
- A Poisson process draws the next event time and angle.
- Events have stochastic `length`, `lifetime`, and `direction`
  drawn from per-event-type PDFs.
- Events apply forces (§3.4, §3.5) for their lifetime, then
  retract.

### §5.2 Lamellipodium spawn rate

Base spawn rate `λ_lam^0 = 0.05 events / (cell · s · radian)`. This
gives roughly 1 lamellipodium per minute for a fully exposed cell
with no contact inhibition (consistent with Krause-Gautreau 2014
lamellipodium lifetimes minutes).

Modulated by:
- **Polarity**: `λ_lam(θ) = λ_lam^0 * (1 + ε_pol * cos(θ - θ_polarity))`
  with `ε_pol = 1.5` (strong polarity bias when cell has nonzero
  polarity).
- **ECM alignment**: `× (1 + ε_align * cos(2(θ - θ_align)))` with
  `ε_align = 0.5` (modest contact guidance).
- **Contact inhibition**: `× (1 - ε_inh * indicator(θ ∈ contact_arc))`
  with `ε_inh = 0.8` (strong suppression on contact edges; not
  perfect, leaves residual ~20% protrusion rate per Stramer & Mayor
  2017 noted leak).
- **Mechanosignal (YAP/TAZ)**: optionally
  `× (1 + ε_mech * mechanosignal_yap_taz)` with `ε_mech = 0.3`
  (mild positive correlation with YAP nuclear localization,
  Aragona 2013).

### §5.3 Filopodium spawn rate

Base spawn rate `λ_filo^0 = 0.2 events / (cell · s · radian)` (4×
higher than lamellipodia — filopodia are more frequent and
shorter-lived).

Modulated by same factors as lamellipodia, plus:
- **Filopodia preferentially spawn at lamellipodium tips**
  (literature: filopodia emerge from existing actin networks).
  Implementation: when a lamellipodium event is active at angle
  `θ_lam`, filopodium spawn rate within ±15° of `θ_lam` is doubled.

### §5.4 Lifetime distributions (literature)

- Lamellipodium lifetime: exponential PDF, mean = 120 s, range
  30 s - 10 min (Krause-Gautreau 2014, Krause 2014).
- Filopodium lifetime: exponential PDF, mean = 30 s, range 5 s - 5
  min (Mattila & Lappalainen 2008).

### §5.5 Length distributions

- Lamellipodium "length" (extension distance from boundary):
  Gamma PDF, mean = 1.5 μm, shape parameter 2 (variance ~ mean²).
- Filopodium length: Gamma PDF, mean = 3 μm, shape parameter 2,
  truncated at max = 15 μm (Mattila 2008).

### §5.6 Implementation interface

```python
# acs/v2/dynamics/protrusion.py

class ProtrusionScheduler:
    def __init__(self, rng_seed: int): ...

    def step(
        self,
        cell: SingleCellState,
        ecm: ECMSubstrateState,
        cluster: Optional[CellClusterState],
        dt_s: float,
    ) -> list[ProtrusionEvent]:
        """Sample new events at this timestep.

        Returns the list of new events spawned. Caller appends to
        cell.protrusions.
        """
        ...

    def evolve_active_events(
        self,
        active_events: list[ProtrusionEvent],
        dt_s: float,
    ) -> tuple[list[ProtrusionEvent], list[ProtrusionEvent]]:
        """Age active events; return (still_active, retracted)."""
        ...
```

The scheduler is a stateful module. The active-contour dynamics
loop calls `step()` to spawn new events and `evolve_active_events()`
to age them, then queries each active event for the per-vertex force
contribution (via §3.4 / §3.5).

### §5.7 Validation against literature

For Phase 1 deliverable:
- Spawn rate per cell = `Σ_θ λ * Δθ` integrated over boundary ≈
  literature-consistent (lamellipodia 1-2/min, filopodia 4-8/min).
- Lifetime histograms match exponential distribution with the
  documented means.
- Angular distribution at quiescence (no polarity, no ECM align,
  no contact) is uniform.
- Angular distribution under polarity matches polarity vector
  direction (positive correlation `r > 0.5`).

These are gates the protrusion module must pass at first execution.

---

## §6 Focal adhesion molecular clutch model

### §6.1 Clutch concept (Chan & Odde 2008)

A focal adhesion is modeled as a "molecular clutch" connecting the
actin cytoskeleton (which retrogrades inward at velocity `v_actin`)
to the substrate (stationary). When the clutch engages, FA traction
is generated proportional to the relative velocity. Bond strengthens
with traction up to a slip threshold; beyond threshold, slip-bond
kinetics increase off-rate, the clutch breaks, and the FA either
ends (release) or matures further.

### §6.2 FA state machine

States and transitions:
- `nascent` (m=0.1, b=0.2): just nucleated, weak. Lifetime 30 s.
  Transitions: → `mature` if traction sustained > threshold for 60s,
  → `released` if traction lost.
- `maturing` (m linearly 0.1→1.0 over τ_mature=120s): under sustained
  traction, FA strengthens. Bond fraction grows.
- `mature` (m=1.0, b=0.8-0.95): fully engaged. Lifetime distribution:
  exponential with mean ~5 min (Plotnikov).
- `slipping` (high b but slip-rate elevated): under high traction,
  individual bonds break stochastically. Effective traction drops
  as b drops. Recovers if traction relaxes.
- `released` (m=0, b=0): FA dissolved. Vertex no longer attached.

### §6.3 FA dynamics module

```python
# acs/v2/dynamics/focal_adhesion.py

class FocalAdhesionScheduler:
    def step(
        self,
        cell: SingleCellState,
        ecm: ECMSubstrateState,
        active_contour: ActiveContourState,
        dt_s: float,
    ) -> tuple[list[FocalAdhesionState], list[str]]:
        """One cell timestep.

        Returns:
          updated_adhesions: new state of all FAs (existing matured
            +newly nucleated +released ones removed).
          released_fa_ids: ids of FAs that ended this step.
        """
```

Per-step operations:
1. **Nucleation**: at each active lamellipodium, with probability
   `p_nucleate = 0.05 / s` per μm of active boundary, spawn a new
   nascent FA at a random vertex inside the active angular window.
2. **Maturation**: existing FAs evolve `maturity` and `bound_fraction`
   per state machine.
3. **Traction calculation**: each active FA generates traction
   `T_a = K_FA * m_a * b_a * (v_actin - v_substrate)`. v_substrate=0.
   v_actin estimated from cortex retrograde flow (default 50 nm/s,
   from Choi 2008).
4. **Force application**: per §3.3, traction force is applied to
   the cell vertex.
5. **Substrate deposit**: `accumulated_traction[I_a, J_a] += T_a *
   dt`. Eaten by ECM remodeling rule.
6. **Slip-bond break**: high-traction FAs lose bonds stochastically.
7. **Release**: when `m_a → 0` or vertex pulls away beyond
   `R_FA = 0.5 μm`, FA released.

### §6.4 Literature anchors

- Chan & Odde 2008 *Science* "Traction dynamics of filopodia on
  compliant substrates" — IF~50, foundational clutch model.
- Plotnikov 2012 *Cell* (cited in §3.3) — FA force fluctuations.
- Wolfenson et al. 2015 *Annu Rev Cell Biol* "Steps in mechanotransduction"
  — IF~25, FA review.
- Elosegui-Artola et al. 2016 *Nat Cell Biol* "Mechanical regulation
  of a molecular clutch defines force transmission and transduction
  in response to matrix rigidity" — IF~25.
- Case & Waterman 2015 *Nat Cell Biol* "Integration of actin dynamics
  and cell adhesion by a three-dimensional, mechanosensitive molecular
  clutch" — IF~25.

### §6.5 Parameter defaults (literature-anchored)

- `K_FA = 0.1 nN/μm` (per §3.3)
- `v_actin = 50 nm/s` (Choi 2008)
- `p_nucleate = 0.05 events / (s · μm)` (active lam area)
- `τ_nascent = 30 s`, `τ_mature = 120 s` (Plotnikov, Wolfenson)
- `mean_FA_lifetime = 300 s` (Plotnikov)
- `slip_threshold = 5 nN per FA` (literature pinpointed
  ambiguous, exploratory flag).

### §6.6 Sanity Gate

- §1 Dimensional: traction in pN, position μm, time s. ✓
- §3 Conservation: per-FA force on cell = -force on substrate;
  verified in test (sum of cell-side FA force + substrate
  accumulated_traction added per step = 0 within numerical precision).
- §5 Sign: traction always pulls cell vertex toward FA position
  (substrate anchor); substrate always feels traction in the cell
  retrograde direction. ✓

---

## §7 Cell-cell junction and contact inhibition

### §7.1 Junction creation

When boundaries of two cells `a` and `b` come within
`d_contact = 0.5 μm`, a JunctionState is created. Junction nucleates
with `m_ecad = 0`, matures over `τ_ecad = 60 s` (Cai 2014).

Multiple junctions per cell pair: junction can have multiple "edges"
along the contact arc, but Phase 1 simplifies to one junction per
pair, with `contact_length_um` aggregating the total contact arc.

### §7.2 Junction force on cells

Per §3.6, junction exerts attractive force on the two cells'
boundary vertices near the contact arc. Force magnitude scales with
`m_ecad * γ_junc`. Direction: each cell's contact vertices pulled
toward the corresponding vertex on the other cell.

### §7.3 Contact inhibition rule

When `m_ecad > 0.5` (i.e., junction is mature), the
`contact_inhibition_active` flag is set True for that junction.

For each cell that has any active inhibition junction:
- Lamellipodia/filopodia spawn rate within ±15° of the contact arc
  midpoint angle is reduced by `(1 - ε_inh) = 0.2` (per §5.2).
- The cell's polarity vector is reoriented toward the OPPOSITE of the
  contact arc midpoint (cells migrate AWAY from contact), with a
  reorientation timescale `τ_polar = 300 s` (Stramer-Mayor 2017).

### §7.4 Junction dissolution

Junction dissolves when:
- `contact_length_um < d_contact * 0.5` (cells separate).
- Or: one cell dies (`cell_state == "dead"`).

Mature junction dissolution is irreversible within the same
contact event (no instant re-nucleation; needs `> 60 s` clear
separation before re-nucleating).

### §7.5 Literature anchors

- Cai et al. 2014 *Cell* (cited in §3.6).
- Carmona-Fontaine et al. 2008 *Nature* "Contact inhibition of
  locomotion in vivo" — IF~50, foundational.
- Stramer & Mayor 2017 *Nat Rev Mol Cell Biol* "Mechanisms and in vivo
  functions of contact inhibition of locomotion" — IF~70.
- Mayor & Carmona-Fontaine 2010 *Trends Cell Biol* — IF~20 review.

### §7.6 Sanity Gate

- §1 Dimensional: tension in pN/μm, length in μm, time in s. ✓
- §2 Boundary cases: zero-contact-length junctions auto-dissolve.
- §3 Conservation: junction force action-reaction (cell a pulled
  toward b → cell b pulled toward a, equal magnitude). Verified in
  test.
- §5 Sign: attractive when m_ecad > 0; never negative.

---

## §8 Slow biology hooks (cell cycle + mechanosignal)

### §8.1 Cell cycle scalar (default OFF)

Cell cycle dynamics module: state advances `cell_cycle_phase` over
time per a stochastic timer (mean G1=12h, S=8h, G2=4h, M=1h —
literature for MCF7-like, but flagged as exploratory; Pavlov 2024
review).

When `cell_cycle_phase` enters "M" (mitosis), the cytokinesis
geometry algorithm fires (§2.5), splitting the cell into two
daughters.

Phase 1 default: OFF. Activated by config flag
`enable_cell_cycle: bool`. When activated, division events appear
in the per-frame HDF5 frame dump.

### §8.2 Mechanosignal (YAP/TAZ)

`mechanosignal_yap_taz` (per cell, [0,1]) updates per timestep based
on local ECM stiffness at the cell footprint:
`d (yap) / dt = k_yap_on * (E_local / E_ref - 1) * (1 - yap) -
                k_yap_off * yap`.

`k_yap_on = 0.001 / s`, `k_yap_off = 0.0005 / s` (slow scale, hours;
Aragona 2013).

Phase 1 default: hook is active but the value only feeds back into
protrusion spawn rate (§5.2). YAP-driven cycle modulation, transcription
changes, etc. are NOT in Phase 1.

### §8.3 Sanity Gate

- §1: yap dimensionless, time in s, k in 1/s. ✓
- §3: yap bounded in [0,1] by construction.
- §5: positive E_local pushes yap up, drift back to 0 on relaxation. ✓

---

## §9 3D visualization architecture

### §9.1 Two pipelines (PI id=636 confocal-quality requirement)

#### Stub pipeline (Phase 1 overnight α deliverable)

- Tool: Plotly or Vispy 3D
- Output: PNG sequence, one per simulation frame (stride 30s sim time)
- Content: cells as colored polygon meshes (active contour vertices
  + extruded by `height_um`); ECM as 2D field underneath; FAs as
  small spheres at FA positions; protrusions as line segments.
- Color coding: cell color = cell_id index (categorical),
  optionally colored by `mechanosignal_yap_taz`.
- Frame export: matplotlib backend or Plotly+kaleido for headless.

#### Confocal-quality pipeline (Phase 2 separate design unit)

- Tool: Blender Cycles via `blender --background --python` (already
  in CLAUDE.md tech stack)
- Volume rendering with fluorescent shaders for each "channel":
  cell cytoplasm (faint), cell membrane (bright), focal adhesions
  (point-like bright spots), ECM (dim background)
- Per-cell mesh from active contour, extruded as 2.5D height
- Camera path scripted (rotate, zoom, time-lapse)
- ffmpeg pipeline for video output

### §9.2 Data-render separation (architectural rule)

Both pipelines read from the same HDF5 frame dump (§10). The
simulation never invokes a renderer; the renderer never modifies
simulation state. This decouples visualization development from
simulation development and supports post-hoc rendering of past runs.

### §9.3 Phase 1 stub viz scope (overnight α)

The overnight deliverable includes:
- `acs/v2/viz/stub3d.py` — Plotly-based viewer
- One demo: single cell on uniform substrate with one lamellipodium
  event, rendered as 60-frame sequence over 30s sim time.
- Documented expectation: this is for development checkup ONLY.
  Confocal quality requires Blender pipeline (post-Phase-1 unit).

---

## §10 Output HDF5 frame dump

### §10.1 File structure

```
runs/<UTC>_<sim_id>/
├── config.yaml          # full input config + git commit hash
├── frames/
│   ├── 0000.h5
│   ├── 0001.h5
│   └── ...              # one per output stride (e.g., 30s sim time)
├── statistics.h5        # accumulated PDFs, correlation lengths
├── events.h5            # divisions, FA nucleation, etc.
└── log/
    ├── stdout.txt
    └── sanity_gate_log.json
```

### §10.2 Frame schema (per `frames/NNNN.h5`)

```
/meta
  frame_index: int
  time_s: float
  config_hash: str
  sim_version: str
/cells
  cell_id[]: str
  state[]: enum(alive, dead)
  age_s[]: float
  cycle_phase[]: enum(None, G0, G1, S, G2, M)
  division_count[]: int
  parent_id[]: str
  mechanosignal[]: float
  position_xy_um[]: (N_cells, 2)
  polarity_xy[]: (N_cells, 2)
  height_um[]: float
  /boundary
    polygon_offsets[]: int (N_cells+1)  # CSR-like
    vertices_xy_um[]: (sum_N_vert, 2)
/protrusions
  cell_id[]: str
  type[]: enum(lam, filo, retract)
  angle_rad[]: float
  length_um[]: float
  lifetime_s[]: float
  age_s[]: float
/focal_adhesions
  fa_id[]: str
  cell_id[]: str
  position_xy_um[]: (N_FA, 2)
  maturity[]: float
  bound_fraction[]: float
  age_s[]: float
  traction_xy_pN[]: (N_FA, 2)
/junctions
  cell_a_id[]: str
  cell_b_id[]: str
  contact_length_um[]: float
  contact_age_s[]: float
  e_cad_maturity[]: float
  contact_inhibition_active[]: bool
/ecm
  grid_origin_um: (2,)
  grid_spacing_um: float
  grid_shape: (2,)
  stiffness_kPa: (Nx, Ny)
  fiber_density: (Nx, Ny)
  alignment_xx: (Nx, Ny)
  alignment_xy: (Nx, Ny)
  accumulated_traction_xy: (Nx, Ny, 2)
```

### §10.3 Statistics dump

`statistics.h5` accumulates over the run:
- Filopodia event rate per cell, per time bin
- Lamellipodia event rate per cell
- FA lifetime histogram
- Filopodia angular distribution
- Cell-cell contact age distribution
- Velocity correlation length per time bin
- ECM stiffness change per region

These are exactly the closure-extraction outputs needed for
future Tier 3 reduced surrogate (per design-discussion id=691).

### §10.4 Compression and I/O

HDF5 chunked storage with gzip compression level 4. Estimated frame
size: 50-cell cluster ≈ 10 MB per frame, 30s stride over 80h sim =
9600 frames = ~96 GB. For Phase 1 deliverable (single cell, 5-10
frames), well under 1 GB.

---

## §11 Implementation milestones

### §11.1 Overnight α — split into P0 (safe) and P1 (conditional)

**Updated per cross-review id=756 RF3**: original 9-patch list split
to prevent overnight time pressure from weakening Sanity Gates.

**P0 (overnight, must-deliver, schema/IO scope)**:
1. data_contract patch (artifact-aware, including new `ECM_FIELD`
   ArtifactKind member — see RF4 below)
2. metrics skeleton patch
3. ECMSubstrateState + CellClusterState + JunctionState + FocalAdhesionState
   extension dataclasses with validation + tests
4. SingleCellState slow-biology hook fields (cell_state, cell_age_s,
   cell_cycle_phase, division_count, parent_cell_id,
   mechanosignal_yap_taz, neighbor_cell_ids) + validation
5. ActiveContourState dataclass + `to_measurement_boundary()` (state
   only, no dynamics)
6. HDF5 frame dump writer (single-cell scope, accepts synthetic
   frame data — no live sim)
7. Plotly stub 3D viewer (reads synthetic frame, renders single
   cell as polygon mesh)
8. Integration smoke test: synthetic frame round-trip (write →
   read → render → schema valid)

**P1 (overnight, conditional on Sanity Gate clarity)**:
9. Active contour numerics core: vertex update with cortex tension
   + area pressure ONLY. Includes Sanity Gate document:
   §1 dimensional check, §2 boundary cases, §3 area conservation,
   §4 dt CFL derivation per-config, §5 cortex contracts / pressure
   restores sign tests, §6 emit valid MeasurementBoundary at output.
   IF Sanity Gate document cannot be cleanly derived overnight →
   defer to next-day work, do not force.

**NOT in overnight α** (any tier): filopodia event scheduler, FA
molecular clutch dynamics, ECM remodeling dynamics, cell-cell
junction dynamics (force coupling), contact inhibition dynamics,
cell cycle dynamics, YAP signaling dynamics, Blender confocal viz.

**ECM_FIELD ArtifactKind extension (RF4)**: prior round 22-26 lock
of 8-member `ArtifactKind` is **explicitly extended** by Phase 1 to
9 members with `ECM_FIELD` added. Semantic: "ECM substrate state
field artifact (simulation output, or in vitro/silico-controlled
ECM input)". Required by ECM-state-derived metrics. Patch 1
(`data_contract`) implements this extension.

### §11.2 Week-by-week post-α

| Week | Module | Sanity Gate | Lit anchor |
|---|---|---|---|
| 2 | Substrate drag + Langevin noise integration | dt stability w/ noise | Trepat 2009 |
| 3 | Lamellipodia event scheduler + hazard model | event rate uniform-prior baseline | Krause 2014 |
| 4 | Filopodia event + adaptive vertex density | 100 nm tip resolution test | Mattila 2008 |
| 5 | FA molecular clutch + traction | Newton 3 conservation | Plotnikov 2012, Chan & Odde 2008 |
| 6 | ECM remodeling (stiffening, alignment) | bounded growth, monotonicity | Hall 2016 |
| 7 | Cell-cell contact + junction | action-reaction, contact age | Cai 2014, Stramer 2017 |
| 8 | Contact inhibition + polarity reorientation | sign of CIL response | Carmona-Fontaine 2008 |
| 9 | Cell cycle scalar + cytokinesis | parent area = sum daughter area | Pavlov 2024 |
| 10 | YAP mechanosignal | bounded [0,1], slow timescale | Aragona 2013 |
| 11 | Cluster validation suite (2, 5, 10, 20 cell tests) | qualitative match to literature | composite |
| 12 | Blender confocal viz pipeline | visual QA against confocal images | (separate review unit) |

### §11.3 Phase 2 entry criteria (deferred)

Move to Phase 2 only after:
- All Phase 1 module Sanity Gates passed
- 5-cell cluster simulation produces reasonable contact-inhibition
  signature
- Filopodia statistics match literature ranges
- ECM remodeling shows expected memory persistence (Hall 2016 fig)
- Documented closure-dump output is consumable by a notional Tier 3
  surrogate consumer

---

## §12 Sanity Gate matrix per module

| Module | §1 Dim | §2 BC | §3 Cons | §4 Num | §5 Sign | §6 MeasProto |
|---|---|---|---|---|---|---|
| `measurement_boundary` | μm/μm² | N≥3, no NaN | rule 11 | float64, ε scale-aware | shoelace sign canonical | top-down xy match exp |
| `data_contract` | non-empty unit | required_artifacts non-empty | N/A | N/A | N/A | artifact subselectors |
| `metrics` | unit per metric | empty registry → error | N/A | numpy/scipy float | metric monotonicity tests | input-source declared |
| `ecm_substrate` | kPa, dimensionless ρ, tensor Q | E≥E_min, ρ∈[0,1] | bounded, traceless Q | dt_ecm vs k_align | stiffening positive, alignment relaxation negative | ECM_FIELD artifact |
| `active_contour` (numerics) | μm vertex, pN force, s time | N≥3, no NaN, no self-intersect | area drift bounded | dt < γ/k_max | each force term sign documented | TBD |
| `protrusion` (events) | s lifetime, μm length, rad angle | length≥0, lifetime>0 | event count balanced | RNG seed reproducible | force outward growth, inward retraction | EVENT_ANNOTATIONS |
| `focal_adhesion` (clutch) | pN traction, μm displacement | maturity∈[0,1] | Newton 3 verified | dt_fa < τ_nascent | attractive cell-substrate | TFM_FIELD |
| `cell_cluster` + `junction` | μm contact length, pN/μm tension | non-self-junction | action-reaction | dt < τ_ecad | attractive at m_ecad>0 | TBD |
| `cell_cycle` | s phase duration | phase ∈ allowed | division: 1 parent → 2 daughter; areas conserved | dt < min phase | monotone phase transition | N/A |
| `mechanosignal` | dimensionless yap | yap ∈ [0,1] | bounded | dt < 1/k_yap_on | positive E → positive yap | N/A |
| `frame_dump` | unit per field | empty cluster OK | no double-write | HDF5 atomicity | N/A | per-frame artifact list |

Each module's `validate()` or first-execution Sanity Gate
docstring must record explicitly all six rows for that module.

---

## §13 Magic-Number Block sweep

Sweep of numerical constants in Phase 1 plan, applying
CLAUDE.md Magic-Number Block tests (derivable / grid-invariant /
not-fitted).

**Updated per cross-review id=756 RF1**: All "PASS" labels replaced
with one of four certainty tiers. No value is "locked" by this
document — these are candidate values to be carried into module
Sanity Gate documents and verified there before first execution.

Tier definitions:
- **L (literature-direct)**: a single named paper reports this exact
  value or a tight range that contains this value.
- **R (literature-range candidate)**: literature reports a range; this
  value is the chosen mid-point or representative.
- **D (derived, requires Sanity Gate)**: must be computed per-config
  from physical/numerical constraints before run.
- **X (exploratory default)**: not directly literature-anchored, model
  choice. Must be flagged in code and in `docs/12_validation.md`.

| Constant | Candidate value | Primary lit ref | Tier |
|---|---|---|---|
| γ_cortex | 0.3 nN/μm | Salbreux 2012 (range 0.1–1.0) | R |
| K_A | 100 pN/μm³ | Bi 2014 / Farhadifar 2007 (dimensionless ratio constraint) | D (ratio derivation) |
| K_FA | 0.1 nN/μm | Plotnikov 2012 (per-FA total / displacement) | R |
| F_lam | 0.5 nN | Pollard 2003 / Mogilner 2003 (range pN–nN) | R |
| F_filo growth | 10 pN | Bornschlögl 2013 (range 5–20 pN) | L |
| F_filo retract | 15 pN | Bornschlögl 2013 (asymmetry observed) | L |
| Δθ_lamellipodium | 30° | model choice | X |
| τ_lam mean | 120 s | Krause 2014 (range 30 s–10 min) | R |
| τ_filo mean | 30 s | Mattila 2008 (range 5 s–5 min) | R |
| γ_junc | 0.2 nN/μm | Bambardekar 2015 (range 0.1–0.5) | R |
| τ_ecad | 60 s | Cai 2014 Fig 4 | L |
| ε_inh contact suppression | 0.8 | Stramer 2017 (qualitative — strong but not perfect) | X |
| τ_polar (CIL reorient) | 300 s | Stramer 2017 | R |
| γ_drag^0 | 0.05 nN·s/μm² | Trichet 2012 | R |
| α drag exponent | 0.7 | Trichet 2012 (literature-fitted) | R |
| E_ref | 5 kPa | Engler 2006 substrate range | R |
| E_max | 30 kPa | Engler 2006 (densified) | R |
| k_stiffen | 0.001 /(s·nN/μm²) | Hall 2016 timescale (derivation) | D (timescale derivation) |
| k_align | 0.005 /s | Hall 2016 / Trichet 2012 timescale | R |
| λ_isotropic | 0.0005 /s | model choice (k_align/λ=10 ratio) | X |
| dt_cell | initial 10 ms test point only | per-config CFL `dt < γ_drag/k_max` | D |
| dt_ecm | initial 1 s test point only | k_align timescale separation, verify per-config | D |
| min_edge_um | 1e-3 μm | numerical floor (already in MeasurementBoundary) | L |
| E_min (substrate) | 0.1 kPa | numerical floor | X |
| ε_pol | 1.5 | model choice | X |
| ε_align (protrusion) | 0.5 | Hall 2016 contact guidance | X |
| ε_mech (yap → protrusion) | 0.3 | Aragona 2013 qualitative | X |

**Total candidate count by tier**: L=4, R=12, D=4, X=8.

**Eight X-tier (exploratory) values** require explicit flagging in
code AND in `docs/12_validation.md` AND in the module Sanity Gate
note. These are the highest-risk Magic-Number Block surfaces in
Phase 1 and must NOT be silently changed during runs.

**No PASS labels**: this document does not lock any numerical value.
All values become candidate inputs to module Sanity Gate documents,
where the Sanity Gate Protocol §1-§6 is run before first execution
and the value is either confirmed, adjusted (with documented reason),
or escalated to PI as a `decision-needed` approval queue item.

---

## §14 Risks and mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Filopodia at sub-μm resolution slows sim | High | Medium | Adaptive vertex density only at active filopodium; cap N_vertices per cell |
| ECM-cell timescale mismatch (s vs ms) | High | Medium | Operator splitting, dt_ecm ≫ dt_cell |
| Cell-cell contact detection O(N²) at large N | Medium | Low (Phase 1 N≤50) | Spatial hashing if needed |
| Cytokinesis splits self-intersecting polygons | Medium | High | Post-split simple-polygon check + vertex repair, or reject division attempt |
| Magic-Number Block exploratory params accumulate | Medium | High | Document all 5 in code + `docs/12_validation.md`; periodic review |
| Sanity Gate skipped under time pressure (overnight α) | Medium | High | Hard rule: blocker status if any module fails; PI approval queue escalation |
| Confocal-quality viz expectation gap | Low (post-warning) | Medium | Stub viz documented as checkup-only; Blender pipeline as separate post-Phase-1 unit |
| Implementation-work agents diverge from this plan | Low | High | Self-contained handoff brief (§15) with explicit boundaries |
| MeasurementBoundary commit not approved before dispatch | Medium | Low (working tree workable) | Coding team uses uncommitted tree; PI commit batch in morning |
| Codex plan and Claude plan diverge on key decision | Low (mostly converged in MCP) | High | Cross-review + consolidation round before dispatch |

---

## §15 Coding-team handoff brief (self-contained for implementation-work)

This section is the brief that gets dispatched to the
implementation-work room as the source of truth for overnight α
coding. It is INTENTIONALLY redundant with the rest of this
document — implementation-work agents may not have time to read the
full plan; this section is what they execute against.

### §15.1 Patch list (overnight α)

Patch 1: **MeasurementBoundary commit** (already drafted, impl id=611
brief, awaiting PI approval). PRIOR ART, no new code.

Patch 2: **data_contract artifact-aware**
- File: `acs/v2/data_contract.py` (modify)
- Add: `ArtifactKind` StrEnum (8 members), `ImagingDatasetSpec`
  fields `available_artifacts`, `available_csv_columns`,
  `MetricSpec` fields `required_artifacts` (non-empty),
  `csv_columns_required`, `channels_required`, removal of
  `source_channel`. Update `MetricSpec.validate(dataset)`.
  Update `default_single_cell_contract` to auto-filter.
  Update `V2DataContract.validate` for new dup key `(name, modality)`.
- Tests: `tests/v2/test_v2_contract.py` migrate; new test cases for
  cross-consistency, dup, missing artifacts, empty contract.
- Sanity Gate: docstring records §1-§6 explicitly.

Patch 3: **metrics skeleton**
- File: `acs/v2/metrics.py` (new)
- Content: registry pattern keyed by `(name, ArtifactKind)`,
  `projected_area` and `perimeter` boundary metrics, defensive
  `boundary.validate()` calls, `V2ContractError(ValueError)` with
  `failure_kind`.
- Tests: `tests/v2/test_v2_metrics.py` new; positive/negative
  registry lookups, metric value correctness.

Patch 4: **ECMSubstrateState + CellClusterState + JunctionState
schema**
- Files: `acs/v2/ecm_substrate.py` (new), `acs/v2/cell_cluster.py`
  (new), `acs/v2/junction.py` (new)
- Content: per §1.4, §1.5, §1.6 of this plan. Frozen dataclasses
  with validate(), no dynamics.
- Tests: schema validation positive/negative cases.

Patch 5: **SingleCellState slow biology hooks**
- File: `acs/v2/single_cell.py` (modify)
- Add: `cell_state, cell_age_s, cell_cycle_phase, division_count,
  parent_cell_id, mechanosignal_yap_taz, neighbor_cell_ids` per §1.1.
- Tests: validation for each new field.

Patch 6: **ActiveContourState + numerics core (cortex tension only)**
- File: `acs/v2/active_contour.py` (new), `acs/v2/dynamics/active_contour.py`
  (new)
- Content: `ActiveContourState` dataclass with mutable vertices and
  `to_measurement_boundary()`. Numerics: overdamped Newton update
  with cortex tension (§3.1) and area pressure (§3.2) ONLY. No FA,
  no protrusion, no ECM coupling, no junction.
- Sanity Gate: explicit §1-§6 docstring; tests for area
  conservation under perturbation, symmetric polygon relaxation
  to circle, dt stability.
- Tests: `tests/test_active_contour.py` new.

Patch 7: **HDF5 frame dump (single-cell scope)**
- File: `acs/v2/output/frame_dump.py` (new)
- Content: minimal frame writer for single-cell sim — `/meta`,
  `/cells/boundary`, no protrusion/FA/junction yet.
- Tests: write-read roundtrip, schema fidelity.

Patch 8: **Plotly stub 3D viewer**
- File: `acs/v2/viz/stub3d.py` (new)
- Content: read frame HDF5, render single-cell polygon as 3D mesh
  (extruded by height), substrate as floor plane, output PNG sequence
  for a 60-frame demo.
- Test: smoke test via subprocess on demo data.

Patch 9: **Integration test**
- File: `tests/v2/test_v2_integration_overnight_alpha.py` (new)
- Content: 1-cell sim with cortex tension only, 100 timesteps,
  area conservation within 1%, frame dump valid, viz produces
  PNGs.

### §15.2 Hard prohibitions

- NO commit unless PI explicitly approves (CLAUDE.md +
  implementation_workflow.md SOP).
- NO Magic-Number Block violations: any new constant must pass
  the 3 tests of §13 above.
- NO Sanity Gate skip: each new physics module records its §1-§6
  notes.
- NO PI experimental data fitting (per PI id=726).
- NO scope expansion beyond patches 1-9 listed above. Filopodia,
  FA, ECM, junction dynamics are NOT in overnight α.
- NO destructive operations (force push, branch delete, etc.).
- NO new external dependencies beyond numpy, scipy, h5py, plotly,
  pytest (all already present).

### §15.3 Escalation rules

- Sanity Gate FAIL → status=blocker → PI approval queue. Module
  pauses; other modules continue.
- Magic-Number Block violation → status=blocker → PI queue.
- Test failure that cannot be diagnosed in 30 min → status=blocker.
- Architectural ambiguity → status=open-question to design-discussion
  room (claude or codex), DO NOT improvise.

### §15.4 Acceptance criteria for overnight α

By morning:
1. All 9 patches drafted as uncommitted changes in working tree.
2. All tests pass (or specifically blocked tests reported).
3. Integration test produces a valid frame dump and PNG render.
4. Sanity Gate notes for every new physics module (active contour
   and ECM substrate validation).
5. No Magic-Number Block violation in any patch.
6. PI approval queue has at most one item: "ready to commit
   patches 1-9 as a single commit batch" (or per-patch if PI
   prefers).

---

## §16 References (consolidated, IF≥15 prioritized)

### §16.1 Cell biology / mechanobiology
- Salbreux G, Charras G, Paluch E. "Actin cortex mechanics and
  cellular morphogenesis" *Trends in Cell Biology* 2012; 22(10):
  536-545.
- Maître JL et al. "Adhesion functions in cell sorting by
  mechanically coupling the cortices of adhering cells" *Cell*
  2012; 338(6104): 253-256.
- Plotnikov SV et al. "Force fluctuations within focal adhesions
  mediate ECM-rigidity sensing to guide directed cell migration"
  *Cell* 2012; 151(7): 1513-1527.
- Wolfenson H, Yang B, Sheetz MP. "Steps in mechanotransduction
  pathways" *Annu Rev Cell Biol* 2015; 31: 463-489.
- Schwarz US, Gardel ML. "United we stand: integrating the actin
  cytoskeleton and cell-matrix adhesions in cellular
  mechanotransduction" *J Cell Sci* 2012; 125: 3051-3060.
- Pollard TD, Borisy GG. "Cellular motility driven by assembly and
  disassembly of actin filaments" *Cell* 2003; 112(4): 453-465.
- Krause M, Gautreau A. "Steering cell migration: lamellipodium
  dynamics and the regulation of directional persistence" *Nat Rev
  Mol Cell Biol* 2014; 15(9): 577-590.
- Mattila PK, Lappalainen P. "Filopodia: molecular architecture and
  cellular functions" *Nat Rev Mol Cell Biol* 2008; 9(6): 446-454.
- Bornschlögl T et al. "Filopodial retraction force is generated by
  cortical actin dynamics" *PNAS* 2013; 110(45): 18928-18933.
- Cai D et al. "Mechanical feedback through E-cadherin promotes
  direction sensing during collective cell migration" *Cell* 2014;
  157(5): 1146-1159.
- Carmona-Fontaine C et al. "Contact inhibition of locomotion in vivo
  controls neural crest directional migration" *Nature* 2008;
  456(7224): 957-961.
- Stramer B, Mayor R. "Mechanisms and in vivo functions of contact
  inhibition of locomotion" *Nat Rev Mol Cell Biol* 2017; 18(1):
  43-55.
- Engler AJ, Sen S, Sweeney HL, Discher DE. "Matrix elasticity
  directs stem cell lineage specification" *Cell* 2006; 126(4):
  677-689.
- Trichet L et al. "Evidence of a large-scale mechano-sensing
  mechanism for cellular adaptation to substrate stiffness" *PNAS*
  2012; 109(18): 6933-6938.
- Hall MS et al. "Fibrous nonlinear elasticity enables positive
  mechanical feedback between cells and ECMs" *PNAS* 2016; 113(49):
  14043-14048.
- Eichinger JF et al. "A computational framework for modeling cell-
  matrix interactions in soft biological tissues" *Biomechanics
  and Modeling in Mechanobiology* 2020. (Soft Matter / BMM venue)
- Kanchanawong P et al. "Nanoscale architecture of integrin-based
  cell adhesions" *Nature* 2010; 468(7323): 580-584.
- Chan CE, Odde DJ. "Traction dynamics of filopodia on compliant
  substrates" *Science* 2008; 322(5908): 1687-1691.
- Choi CK et al. "Actin and α-actinin orchestrate the assembly and
  maturation of nascent adhesions in a myosin II motor-independent
  manner" *Nat Cell Biol* 2008; 10(9): 1039-1050.
- Elosegui-Artola A et al. "Mechanical regulation of a molecular
  clutch defines force transmission and transduction in response
  to matrix rigidity" *Nat Cell Biol* 2016; 18(5): 540-548.
- Case LB, Waterman CM. "Integration of actin dynamics and cell
  adhesion by a three-dimensional, mechanosensitive molecular
  clutch" *Nat Cell Biol* 2015; 17(8): 955-963.
- Trepat X et al. "Physical forces during collective cell
  migration" *Nature Physics* 2009; 5(6): 426-430.
- Aragona M et al. "A mechanical checkpoint controls multicellular
  growth through YAP/TAZ regulation by actin-processing factors"
  *Cell* 2013; 154(5): 1047-1059.
- Dupont S et al. "Role of YAP/TAZ in mechanotransduction" *Nature*
  2011; 474(7350): 179-183.
- Bambardekar K et al. "Direct laser manipulation reveals the
  mechanics of cell contacts in vivo" *PNAS* 2015; 112(5):
  1416-1421.
- Mayor R, Carmona-Fontaine C. "Keeping in touch with contact
  inhibition of locomotion" *Trends in Cell Biology* 2010; 20(6):
  319-328.
- Friedl P et al. "Collective cell migration in morphogenesis,
  regeneration and cancer" *Cell* 2009; (review piece reference).
- Cheung KJ et al. "Collective invasion in breast cancer requires
  a conserved basal epithelial program" *PNAS* 2013; 110(48):
  19490-19495.
- Pavlov M et al. "Mechanical regulation of the cell cycle" *Annual
  Review* 2024 (review type).

### §16.2 Active matter / continuum theory
- Marchetti MC et al. "Hydrodynamics of soft active matter" *Reviews
  of Modern Physics* 2013; 85(3): 1143-1189.
- Prost J, Jülicher F, Joanny JF. "Active gel physics" *Nature
  Physics* 2015; 11(2): 111-117.
- Bi D et al. "Energy barriers and cell migration in densely packed
  tissues" *Physical Review X* 2014.
- Bi D, Yang X, Marchetti MC, Manning ML. "Motility-driven glass and
  jamming transitions in biological tissues" *Physical Review X*
  2016; 6: 021011.
- Park JA et al. "Unjamming and cell shape in the asthmatic airway
  epithelium" *Nature Materials* 2015; 14(10): 1040-1048.

### §16.3 Vertex / cellular models
- Honda H. "Description of cellular patterns by Dirichlet domains"
  *J Theoretical Biology* 1978; 72(3): 523-543.
- Fletcher AG et al. "Vertex models of epithelial morphogenesis"
  *Biophysical Journal* 2014; 106(11): 2291-2304.
- Farhadifar R et al. "The influence of cell mechanics, cell-cell
  interactions, and proliferation on epithelial packing" *Current
  Biology* 2007; 17(24): 2095-2104.
- Graner F, Glazier JA. "Simulation of biological cell sorting using
  a two-dimensional extended Potts model" *Phys Rev Lett* 1992;
  69(13): 2013-2016.

### §16.4 Simulation / numerical
- Kass M, Witkin A, Terzopoulos D. "Snakes: active contour models"
  *International Journal of Computer Vision* 1988; 1(4): 321-331.
- Ghaffarizadeh A et al. "PhysiCell: an open source physics-based
  cell simulator for 3-D multicellular systems" *PLoS Comput Biol*
  2018; 14(2): e1005991.
- Swat MH et al. "Multi-scale modeling of tissues using CompuCell3D"
  *Methods Cell Biol* 2012; 110: 325-366.
- Mirams GR et al. "Chaste: an open source C++ library for
  computational physiology and biology" *PLoS Comput Biol* 2013;
  9(3): e1002970.

### §16.5 Closure-extraction (for future Tier 3)
- Travasso RDM et al. "Tumor angiogenesis and vascular patterning:
  a mathematical model" *PLoS One* 2011 / *Phys Biol*.
- Wise SM et al. "Three-dimensional multispecies nonlinear tumor
  growth—I model and numerical method" *J Theoretical Biology* 2008.

---

[STATUS 23:48 KST: §0-§16 drafted. Plan complete. ~1900 lines.
Proceeding to cross-review with Codex.]

---

## Footer — Corrections applied from cross-review id=756

After cross-review with Codex's `docs/v2/v2_phase1_plan_codex.md`,
the following sections were updated in this file:

- **§4.6 ECM Sanity Gate**: removed self-correcting unit-conversion
  prose (Rule 10 unit-chain risk). Conversion deferred to
  `ecm_substrate` module Sanity Gate document.
- **§11.1 Milestones**: split overnight α into P0 (must-deliver
  schema/IO) and P1 (conditional active-contour numerics).
  Original 9-patch unsplit list replaced.
- **§11.1 ECM_FIELD enum extension**: explicit note that prior
  round 22-26 8-member `ArtifactKind` lock is extended by Phase 1
  to 9 members. Patch 1 implements.
- **§13 Magic-Number Block sweep**: replaced "PASS"/"PASS w/ note"
  labels with four-tier certainty (L/R/D/X). No value is locked
  by this document; all become candidate inputs to module Sanity
  Gate.
- **References (§16)**: most references confirmed; a candidate
  subset (Eichinger 2020 venue, Pavlov 2024 review type, Wang 2019
  active contour internal pressure, Stewart 2011 venue, some IF
  claims) is to be verified before final BibTeX entry. Marked in
  consolidated plan, not here.

The consolidated `docs/v2/v2_phase1_plan_consolidated.md` (Codex
in-progress) is the source of truth for implementation-work
dispatch. This file is preserved as Claude's independent draft for
audit.

[FINAL STATUS 23:53 KST: Claude plan complete with cross-review
corrections applied. Awaiting Codex consolidated plan for final
review and dispatch.]




