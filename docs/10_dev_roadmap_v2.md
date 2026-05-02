# 10 v2 — Development Roadmap

This roadmap supersedes the spheroid-first roadmap for new development. The
original roadmap remains as the v1 continuum-prototype history in
`docs/10_dev_roadmap.md`.

## Phase V2-0 — Freeze V1 And Open V2

Goal: preserve the existing continuum model and start a separate v2 code path.

Deliverables:

- [x] `docs/00_project_vision_v2.md`
- [x] `docs/v1_continuum_backup.md`
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

