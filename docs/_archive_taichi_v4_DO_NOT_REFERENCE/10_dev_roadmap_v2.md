# 10 v2 — Development Roadmap

This roadmap supersedes the spheroid-first roadmap for new development. The
original roadmap remains as the v1 continuum-prototype history in
`docs/v1/10_dev_roadmap.md`.

## Phase V2-0 — Freeze V1 And Open V2

Goal: preserve the existing continuum model and start a separate v2 code path.

Deliverables:

- [x] `docs/v2/00_project_vision_v2.md`
- [x] `docs/v1/v1_continuum_backup.md`
- [x] `docs/claude_codex_log.md` moved to root `docs/`
- [x] `acs/v2/` package opened
- [ ] decide whether to create a dedicated git branch for v2 work
- [ ] update AGENTS.md / CLAUDE.md to point new work to the v2 roadmap

Validation:

- v1 files remain present and unmodified unless intentionally patched
- v2 package imports without requiring Taichi/GPU

## Phase V2-1 — Imaging Data Contract

Goal: define exactly what visual/physical data enters the model.

Deliverables:

- [x] initial Python data-contract objects in `acs/v2/data_contract.py`
- [ ] canonical file layout for confocal/live-imaging datasets
- [ ] required metadata fields: voxel size, frame interval, channels,
  segmentation provenance, calibration/validation split
- [ ] metric registry for single-cell and spheroid morphology
- [ ] example YAML/JSON data contract

Validation:

- contract validation catches missing voxel size, frame interval, channels,
  empty metric names, and calibration/validation overlap
- no model calibration starts before each dataset has an explicit contract

## Phase V2-2 — Single-Cell Minimal Mechanobiology

Goal: model single-cell spreading before spheroids.

Deliverables:

- [x] initial single-cell state/event objects in `acs/v2/single_cell.py`
- [ ] 2D/2.5D boundary representation
- [ ] lamellipodia/filopodia event extraction from imaging
- [ ] focal adhesion state machine
- [ ] polarity and traction representation
- [ ] single-cell validation dashboard

Validation:

- reproduces event timescale distributions
- reproduces projected area and boundary roughness trajectories
- separates fitted, literature-fixed, and exploratory parameters

### Layer 1 Sanity-Gate Reference Anchors

These are the literature anchors that V2-2 dynamics must satisfy *before*
first execution, per the Sanity Gate Protocol in `CLAUDE.md`. Codex review
`mcp_msg:289` (2026-05-02) tightened the original seed list — the
corrections below are the contract.

| Sub-channel | Primary anchor (IF≥15 unless noted) | Notes / what it actually justifies |
|---|---|---|
| Active-matter / time-scale separation | Marchetti et al. 2013 *Rev Mod Phys* (~50); Prost/Jülicher/Joanny 2015 *Nat Phys* (~22) | Sanity Gate 1 (dimensional) anchors only. Do **not** import CFL/stability bounds — `dt_crit` must be derived from our own force law + mobility after implementation. |
| Stochastic edge / shape coupling | Keren et al. 2008 *Nature* (~50); Raynaud et al. 2016 *Nat Phys* (~22) | Single-cell shape/time-scale coupling and stochastic edge switching. Polarity persistence belongs here (Goal: Sanity Gate 6 — polarity-vector autocorrelation `C(Δt)` and displacement-polarity cosine over a *pre-registered* window, not a fitted `tau_p`). |
| Cortex / surface tension | **Chugh et al. 2017 *Nat Cell Biol*** (~21) for single-cell cortex architecture → surface tension; Maître et al. 2012 *Science* (~56) is **cell-cell contact only**, moves to Layer 2 / V2-3. | Sanity Gate 5 (sign/sense): cortex tension is contractile / shape-restoring at the boundary. Col1-spreading-induced cortex changes have no clean IF≥15 anchor; flag exploratory if introduced. |
| Protrusion / lamellipodia events | Machacek et al. 2009 *Nature* (~50) for protrusion-cycle Rho/Rac/Cdc42 signaling; Tkachenko et al. 2011 *Nat Cell Biol* (~21) for PKA wave coupling. | The inhomogeneous Poisson hazard is a **simulation ansatz**, not a measured law. Validate via inter-arrival distribution / CV and spatial `g(r)` first peak — *never* fit `λ0` to PI A(t). |
| Focal adhesion state machine | Plotnikov et al. 2012 *Cell* (~64) — force fluctuations / rigidity sensing; Choi et al. 2008 *Nat Cell Biol* (~21) — nascent assembly + maturation lifetimes; Kanchanawong et al. 2010 *Nature* (~50) — FA molecular architecture (state identity); Case & Waterman 2015 *Nat Cell Biol* (~21) — molecular clutch organization; Elosegui-Artola et al. 2016 *Nat Cell Biol* (~21) — talin/clutch rigidity threshold. | 5-stack. Plotnikov **alone** does not give all rates. Transition rates `unbound→nascent→mature→slipping→released` are condition/cell-type dependent — express as ranges, never single constants. |
| Filopodial traction (if introduced) | Chan & Odde 2008 *Science* (~47) | Justifies filopodial vs lamellipodial traction dynamics if a separate channel is needed. |
| Vertex / active-contour numerics | (none — biology refs do **not** validate discretization) | Run our own convergence / symmetry / area-drift tests. Bi/Park 2014/2015 *Phys Rev X* / *Nat Mater* are valid Layer 2 *tissue* references, **insufficient for Layer 1 single-cell active-contour numerics**. |
| Input segmentation | Cellpose 2.0 / SAM (method papers) | **Input segmentation only.** Does **not** justify solver vertex spacing, `dx`, or any numerical resolution choice. |

Magic-number risk list (declare ranges + literature anchors before
the run; never tune to PI A(t) — see CLAUDE.md Magic-Number Block):
`λ0` (baseline protrusion event rate), `τ_ref` (refractory), `ξ_p`
(protrusion angular width / correlation length), `γ_cortex`, area /
perimeter stiffnesses, FA `k_on` / `k_off`, FA slip threshold,
polarity memory `τ_p`, traction-coupling gain.

Measurement-protocol note (Sanity Gate 6): top-down polygon area
(`SingleCellState.projected_area_um2()` shoelace) is the correct
analog for single-cell top-down segmentation. Protrusion/FA gates
require live boundary tracking / marker / TFM channels — if PI data
do not measure raw `λ`, compare inter-arrival distribution and
spatial correlation only.

## Phase V2-3 — Cell-Resolved Multi-Cell Assembly

Goal: build spheroids from individually simulated cells.

Deliverables:

- [ ] cell-cell junction model
- [ ] cell-ECM adhesion heterogeneity
- [ ] neighbor exchange / rearrangement metrics
- [ ] small-cluster tests before full spheroid simulation

Validation:

- simulated clusters match confocal morphology beyond scalar area
- condition-specific adhesion changes produce interpretable differences

### Layer 2 Sanity-Gate Reference Anchors

Per codex review `mcp_msg:289` (2026-05-02), Layer 2 references
should *define measurable junction force / direction-sensing metrics
first*, before pre-stating any single junction-lifetime threshold —
those are imaging-modality-dependent and become magic numbers
otherwise.

| Sub-channel | Primary anchor (IF≥15) | Notes |
|---|---|---|
| Cell-cell junction force / direction sensing | Cai et al. 2014 *Cell* (~64); Bazellières et al. 2015 *Nat Cell Biol* (~21) | Define `junction length`, `junctional tension proxy`, `E-cad/catenin intensity`, `T1-like neighbor exchange events`, `contact duration before protrusion reorientation` as the measurement contract. Lifetime thresholds (e.g., "≥ X min = mature") are condition-dependent. |
| Cell-cell contact mechanics (interfacial tension) | Maître et al. 2012 *Science* (~56) | Moves here from Layer 1 — anchors the cortex-vs-adhesion balance at cell-cell contacts, not single-cell cortex. |
| Collective traction transmission | Trepat et al. 2009 *Nat Phys* (~22) | Anchors collective traction sign + spatial pattern. **Not** a Layer 1 FA-kinetics reference. |
| Tissue jamming / packing | Park et al. 2015 *Nat Mater* (~47); Bi et al. 2014/2015 *Phys Rev X* (~20) | Density / crowding / shape-index gate. Used as Layer 2 tissue mechanics anchor, not Layer 1 single-cell. |
| Contact inhibition of locomotion | Carmona-Fontaine et al. 2008 *Nature* (~50) | Anchors neighbor-contact protrusion suppression / polarity repulsion bias. |
| Cell-ECM substrate (Col1) | Levental et al. 2009 *Cell* (~64) for ECM stiffness/crosslinking; Kanchanawong 2010 *Nature* for FA architecture (cross-listed with Layer 1) | Substrate rigidity range for Col1-coated: ~0.5–5 kPa. |

Layer 2 gate observables (preferred over absolute `λ`): junction
lifetime distribution, T1-event rate, CIL protrusion suppression
hazard ratio (free-edge λ vs contact-edge λ), collective traction
correlation length, jamming/shape-index. A/A0 is a **secondary
projection metric only**, not a primary Layer 2 gate.

## Phase V2-4 — ECM Fiber Network

Goal: model ECM as explicit fibers when visual data requires it.

Deliverables:

- [ ] fiber graph/filament data structure
- [ ] traction-induced fiber alignment
- [ ] degradation/remodeling events
- [ ] ECM-guided migration metrics

Validation:

- simulated fiber orientation and remodeling tracks match imaging metrics

## Phase V2-5 — Brute-Force Spheroid Reference

Goal: simulate spheroid spreading with cell-resolved mechanics.

Deliverables:

- [ ] full cell-resolved spheroid simulator
- [ ] compute-cost benchmark on 16 GB A5000 baseline
- [ ] realism/error report against held-out imaging

Validation:

- cell-resolved simulation is the reference model for reduction tests

## Phase V2-6 — Reduced Surrogate Comparison

Goal: test when v1-style continuum/CFD/Marangoni models replace brute force.

Deliverables:

- [ ] coarse-grained observables from cell-resolved runs
- [ ] MPM/active-gel/CFD candidate reductions
- [ ] error-vs-speed comparison

Validation:

- surrogate models are accepted only where their error is quantified against
  the cell-resolved reference

