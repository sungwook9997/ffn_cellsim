# V2 Phase 1 Alpha Implementation Brief

Date: 2026-05-03 KST
Source of truth: `docs/v2/v2_phase1_plan_consolidated.md`
Scope: overnight Project alpha handoff to implementation-work

## Objective

Implement a reviewable v2 skeleton for the Phase 1 cell/ECM reference engine.
The alpha target is schema, contract, output, and visualization scaffolding,
with active-contour numerics only if its Sanity Gate can be written before
first execution.

Do not use PI experimental data. Do not fit any parameter. Do not add hidden
numeric defaults.

## P0: Required Alpha Work

P0 can proceed without executable physics.

1. `acs/v2/data_contract.py`
   - Add `ArtifactKind` as a 9-member enum:
     `CSV_TABLE`, `RAW_FRAMES`, `SEGMENTATION_MASK`, `BOUNDARY_CONTOURS`,
     `TRACKING_TABLE`, `MARKER_CHANNEL`, `TFM_FIELD`, `EVENT_ANNOTATIONS`,
     `ECM_FIELD`.
   - Add `ImagingDatasetSpec.available_artifacts`.
   - Add `ImagingDatasetSpec.available_csv_columns`.
   - Replace metric `source_channel` semantics with:
     `required_artifacts`, `csv_columns_required`, `channels_required`.
   - Validate a metric against the full dataset.
   - Duplicate key is `(name, measurement_modality)`.
   - `default_single_cell_contract(dataset)` must auto-filter to available
     artifacts and must not claim unavailable protrusion/roughness metrics.

2. `acs/v2/metrics.py`
   - Add pure metric registry keyed by `(metric_name, ArtifactKind)`.
   - Implement boundary `projected_area` and `perimeter` using
     `MeasurementBoundary`.
   - Call `boundary.validate()` defensively.

3. Schema-only modules
   - Add or update:
     - `acs/v2/ecm_substrate.py`
     - `acs/v2/focal_adhesion.py` or existing `single_cell.py` class
     - `acs/v2/protrusion.py` or existing `single_cell.py` class
     - `acs/v2/junction.py`
     - `acs/v2/cell_cluster.py`
   - Include validation only, no dynamics.

4. `SingleCellState` slow-biology hooks
   - `cell_state`, `cell_age_s`, `cell_cycle_phase`, `division_count`,
     `parent_cell_id`, `mechanosignal_yap_taz`, `neighbor_cell_ids`.

5. `acs/v2/output/frame_dump.py`
   - Write a minimal HDF5 frame for one or more cells plus ECM metadata.
   - Use flat boundary vertices plus offsets for variable vertex counts.
   - Include enough metadata for visualization round-trip.

6. `acs/v2/viz/stub3d.py`
   - Read the frame dump.
   - Render a synthetic/schema-valid single-cell or small-cluster scene.
   - Output PNG sequence and/or lightweight HTML.

7. Tests
   - Data contract validation.
   - Metric registry and boundary metrics.
   - Schema positive/negative validation.
   - HDF5 write/read round-trip.
   - Stub visualization smoke test on synthetic data.

## P1: Conditional Alpha Work

P1 proceeds only if the Sanity Gate can be written before any execution.

1. Add `acs/v2/active_contour.py`
   - `ActiveContourState`
   - mutable `vertices_um_xy`
   - adaptive vertex-density metadata
   - `to_measurement_boundary()`
   - `from_measurement_boundary()`

2. Add `acs/v2/dynamics/active_contour.py` and/or `acs/v2/dynamics/cortex.py`
   - cortex / line-tension term
   - area-restoring term
   - no protrusion, FA, ECM remodeling, junction, contact inhibition, cell
     cycle, or YAP dynamics

3. Sanity Gate required before first execution
   - unit chain for vertices, force, mobility, and time
   - timestep bound from implemented stiffness and mobility
   - boundary cases: N<3, NaN, self-intersection, tiny edges
   - conservation: area drift and symmetry preservation
   - sign/sense check for each force
   - measurement protocol for exported `MeasurementBoundary`

If the unit derivation or timestep bound fails, stop P1 and report blocker.
P0 remains valid alpha progress.

## Explicit Non-Scope

- Filopodia event scheduler
- FA molecular clutch dynamics
- ECM remodeling dynamics
- Cell-cell junction dynamics
- Contact inhibition dynamics
- Cell-cycle or YAP feedback
- Blender confocal-quality rendering
- Any PI experimental data comparison
- Any production/sweep run
- Any commit without PI approval

## Escalation

- Sanity Gate FAIL: send `status=blocker` to PI and pause that module.
- Magic-Number Block FAIL: send `status=blocker` to PI and pause that module.
- Need a new numeric default: label it as `literature-direct`,
  `literature-range candidate`, `exploratory default`, or
  `requires Sanity Gate derivation before execution`.
- Ambiguous physics: ask design-discussion. Do not invent a fitted constant.

## Morning Acceptance

Minimum:

- P0 drafted or implemented.
- Tests added for any implemented code.
- No hidden numeric defaults.
- No PI experimental data use.
- Blocked P1 physics reported rather than forced.

Strong:

- P0 tests pass.
- P1 active-contour state and cortex/area tests pass with Sanity Gate.
- One synthetic frame dump renders through the stub visualization.
