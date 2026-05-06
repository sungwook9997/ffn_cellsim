# 00 v2 — Image-Constrained Cell-Resolved Mechanobiology Framework

## Core Thesis

Build a **3D mechanobiology toolchain that starts from single-cell spreading physics and scales upward to spheroid spreading**, using confocal/live-imaging visual data and physical measurements as first-class constraints. The existing continuum/MPM/Marangoni spheroid model is preserved as a **reduced surrogate branch**, not as the primary ground-truth model.

The new center of gravity is:

> Establish the cell-level mechanisms first, then test whether spheroid-scale continuum approximations can faithfully replace brute-force cell-resolved simulation.

This is a major framing shift from the original Stage 1 plan. The original plan treated the spheroid as the starting object and asked whether radial/continuum hydrodynamic approximations explain the experiment. The v2 plan treats **single-cell protrusion, adhesion, traction, cell-cell interaction, and ECM remodeling** as the primary mechanistic substrate; spheroid spreading becomes an emergent multi-cell outcome.

## Why The Shift Matters

The current model produced useful lessons: top-down measurement matching matters, gate contracts matter, and continuum active-fluid terms can create plausible spreading trajectories. But it also exposed a structural ceiling: spheroid-scale MPM can saturate without actually resolving the cellular events that dominate real epithelial spreading.

The missing mechanobiology is not a small correction term. It includes:

- lamellipodia extension/retraction cycles
- filopodia probing and ECM sensing
- focal adhesion nucleation, maturation, catch/slip lifetime, and release
- actin/cortex-driven polarity persistence
- cell-cell junction remodeling
- cell rearrangement, neighbor exchange, and leader/follower behavior
- ECM fiber deformation, alignment, degradation, and remodeling

These should be modeled directly before asking whether a coarse-grained Marangoni or active-fluid description can replace them.

## New Goal Hierarchy

### Goal 1 — Single-Cell Spreading Ground Truth

Build a single-cell mechanobiology calibration platform.

Inputs:

- confocal and live-imaging time series
- cell outline, height, projected area, and 3D shape
- lamellipodia and filopodia event statistics
- actin/cortex intensity fields where available
- focal adhesion location, lifetime, and maturation markers where available
- migration velocity and polarity persistence
- traction force microscopy / substrate deformation if available

Outputs:

- protrusion/retraction event model
- focal adhesion turnover model
- active boundary traction model
- cortex/volume/shape mechanics model
- image-derived validation metrics, not just scalar A/A0

Success criterion:

The model should reproduce visual and physical single-cell spreading modes across time, not merely match final area.

### Goal 2 — Cell-Resolved Spheroid Buildup

Use calibrated single-cell modules to build spheroids under different cell-cell and cell-ECM adhesion conditions.

Key state variables:

- cell-cell adhesion strength and junction maturity
- cell-ECM adhesion strength and focal adhesion state
- per-cell polarity, protrusion state, traction, and cortical stiffness
- local crowding, compression, and neighbor exchange
- phenotype/condition-dependent adhesion composition

Data role:

Confocal morphology is used to compare whether the simulated spheroid looks structurally similar to the real spheroid: surface roughness, sheet thickness, leading-edge fingering, cell packing, height profile, and projected shape.

Important policy change:

The original project forbade fitting PI experimental data. In the v2 framework, this must be revised. The new aim is not independent-only prediction; it is **image-constrained digital-twin calibration**. Parameter adjustment against imaging data is allowed only when:

- the adjusted parameter has a biological/mechanical interpretation
- the target measurement modality is explicitly matched
- fitted parameters are separated from literature-fixed parameters
- held-out conditions or time windows are used for validation
- no scalar fudge factor is introduced without passing the Magic-Number Block

### Goal 3 — Explicit ECM Fiber Remodeling

Replace scalar substrate adhesion with an ECM fiber-network view where needed.

Model objects:

- collagen fibers as graph/beam/filament elements
- fiber orientation, density, crosslinking, and stiffness
- cell traction transmitted to local fibers
- fiber alignment, bundling, degradation, and deposition
- ECM-guided migration and durotaxis/contact guidance

Success criterion:

The simulation should reproduce not only cell trajectories but also the visual evolution of ECM architecture: aligned fibers, remodeled tracks, local densification, and invasion paths when present.

### Goal 4A — Brute-Force Cell-Resolved Simulation

Treat each cell as an individually simulated object with its own state.

Purpose:

- establish the most mechanistically faithful reference model
- quantify computational cost
- determine whether the cell-level rules reproduce spheroid-scale behavior without spheroid-level fitting
- provide ground truth for coarse-graining

This branch is allowed to be expensive. It is the standard against which reduced models are judged.

### Goal 4B — Continuum / CFD / Marangoni Surrogate

Use the brute-force simulator to derive and validate coarse-grained models.

Candidate reduced models:

- MPM active viscoelastic continuum
- poroelastic or active gel PDEs
- cellular Marangoni surface-tension-gradient models
- phase-field / active-fluid approximations
- CFD-style solvers for continuum reductions

Open-source candidates to evaluate:

- PhysiCell for 3D multicellular agent-based structure
- CompuCell3D for Cellular Potts / adhesion-driven multicellular morphology
- Chaste for computational biology and tissue modeling
- Cytosim for cytoskeletal/filament mechanics inspiration
- FEniCS for PDE/FEM continuum mechanics
- OpenFOAM or Basilisk for CFD-style reduced solvers, if the final coarse-grained equations justify them

Important distinction:

OpenFOAM/CFD is not the core cell model. It is a possible acceleration layer after the cell-resolved model tells us what continuum equations are legitimate.

## Proposed Architecture

```text
Layer 0: Imaging and physical-data ingestion
  - confocal stacks, live imaging, segmentation, tracking, morphology metrics

Layer 1: Single-cell mechanobiology
  - protrusion, filopodia, focal adhesions, polarity, cortex, volume, traction

Layer 2: Cell-cell and cell-ECM spheroid assembly
  - junctions, adhesion composition, neighbor exchange, collective spreading

Layer 3: ECM fiber network
  - fiber mechanics, remodeling, degradation, alignment, contact guidance

Layer 4A: Brute-force cell-resolved simulator
  - expensive mechanistic reference model

Layer 4B: Reduced surrogate simulator
  - MPM / active gel / CFD / Marangoni / radial approximations

Layer 5: Calibration, validation, and comparison dashboard
  - image metrics, physical metrics, held-out validation, uncertainty
```

## Relationship To The Existing Codebase

The existing codebase should not be discarded. It should be relabeled:

- `v1`: continuum spheroid-first prototype
- `v2`: image-constrained cell-resolved framework

Existing useful assets:

- HDF5/CSV/metadata output discipline
- gate report concept
- top-down projection measurement correction
- GPU portability layer
- visualization/dashboard structure
- current MPM solver as a reduced-model baseline
- existing docs on validation failures and measurement-protocol mismatch

Likely refactor direction:

- keep `acs/physics/mlsmpm.py` as reduced-model code
- introduce new modules for image ingestion, single-cell mechanics, cell-resolved simulation, ECM fibers, and calibration
- avoid adding more biology into the monolithic MPM file

## Revised Validation Philosophy

The v2 framework needs two validation modes.

### Mechanistic Validation

Used for single-cell and cell-resolved modules.

Questions:

- Does the lamellipodium extend in the right direction and timescale?
- Do focal adhesions form and release with realistic lifetimes?
- Does traction align with ECM deformation?
- Does the model reproduce visual morphology in held-out imaging?
- Are cell-cell and cell-ECM adhesion effects separable?

### Reduction Validation

Used for continuum/CFD/Marangoni approximations.

Questions:

- Does the reduced model reproduce brute-force cell-resolved behavior?
- Which observables survive coarse-graining?
- Where does radial symmetry break?
- When does Marangoni-like flow replace explicit protrusion dynamics?
- How much compute is saved per unit error?

## Immediate Roadmap Proposal

### Phase V2-0 — Freeze And Label The Current Branch

- mark the current MPM/Marangoni implementation as `v1 continuum prototype`
- document what it explains and where it saturates
- stop adding new biological detail into the monolithic MPM solver until the v2 architecture is agreed

### Phase V2-1 — Data Contract

- define required imaging inputs
- define segmentation/tracking outputs
- define morphology metrics for single cells and spheroids
- define train/calibration vs held-out validation split

### Phase V2-2 — Single-Cell Minimal Model

- 2D/2.5D cell boundary representation
- lamellipodia/filopodia event model
- focal adhesion state machine
- polarity and traction model
- compare against single-cell live imaging

### Phase V2-3 — Multi-Cell Adhesion Model

- add cell-cell junctions
- add cell-ECM adhesion heterogeneity
- simulate small clusters before full spheroids
- validate against confocal morphology

### Phase V2-4 — ECM Fiber Network

- explicit fiber graph/filament substrate
- traction-induced alignment
- remodeling/degradation
- migration response to remodeled ECM

### Phase V2-5 — Brute-Force Spheroid

- full cell-resolved spheroid spreading
- quantify realism and compute cost
- identify emergent coarse-grained observables

### Phase V2-6 — Reduced Surrogate

- derive continuum/CFD/MPM/Marangoni approximations from brute-force outputs
- validate surrogate error
- build acceleration strategy for parameter sweeps

## Non-Negotiable Rules Carried Forward

- Measurement modality must match experimental modality.
- Any dimensional comparison must use consistent units.
- Any new physics/numerics module needs a Sanity Gate before execution.
- Any empirical multiplier must pass the Magic-Number Block or be explicitly classified as a fitted calibration parameter.
- Fitted parameters, literature-fixed parameters, and PI-chosen exploratory parameters must be labeled separately.
- Visual similarity is not enough; physical interpretation must survive inspection.

## Bottom Line

The v2 framework is stronger if the project accepts a two-engine structure:

1. a mechanistic, image-constrained, cell-resolved engine that tries to be biologically faithful
2. a reduced continuum/CFD/Marangoni engine that tries to be fast

The scientific contribution is then not merely "we simulated spheroid spreading." It becomes:

> We built a cell-resolved mechanobiology reference model, then derived and tested when continuum hydrodynamic approximations can replace explicit cellular computation.

