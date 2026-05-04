# V2 Phase 1 — Cytokinesis Design Note (schema-only, NOT executable)

**Status**: design note only, **NOT a Sanity Gate, NOT a code
commitment**. Per
`docs/v2_phase1_forward_roadmap.md` line 51 ("Cytokinesis
design-note only, no executable implementation requirement") and
the "Phase 1 Out / Schema-Only" classification (line 53–61).

**Author**: implementation-work Claude, drafted as a design-lock-
independent idle unit during the 2026-05-04 KST session
(post-watchdog id=1272). PI ratify status: full delegation per PI
id=939/1008. Codex parallel work: V2-1 imaging data contract
fixture (id=1274), no collision.

**Hard contract**: nothing in this note authorizes a code commit
that executes cytokinesis. Phase 1 explicitly defers executable
cytokinesis to Phase 1.5 / Phase 2 (forward roadmap line 56–57).
Any future code that *implements* cytokinesis must write its own
Sanity Gate document and design-discussion lock.

---

## 0. Why this design note exists

The forward roadmap commits Phase 1 to producing a design note for
cytokinesis without an executable implementation. Without this
note, the schema fields touched by cytokinesis (cell_id, FA links,
protrusion links, ECM trace at division site) carry no documented
contract, and a future Phase 1.5 cytokinesis-execution PR would
re-litigate decisions that are already constrained by Phase 1
schema choices.

This note documents:
- which v2 schema fields cytokinesis interacts with;
- what shape a divide-event records in those fields;
- which decisions are deferred to Phase 1.5 / Phase 2 and require
  their own design-discussion round before executable code lands;
- what activation criteria future executable cytokinesis must
  satisfy.

It does **not** propose a constitutive law for cytokinesis timing,
cleavage-furrow mechanics, or daughter-cell volume partitioning.
Those are explicitly Phase 1.5+ scope.

---

## 1. Schema interaction surface (v2 Phase 1 state)

Cytokinesis affects the following canonical v2 schemas. Each row
records the field and the *kind of change* a divide event would
introduce; the note does NOT prescribe how that change is computed.

| Schema | Field(s) | Cytokinesis impact (deferred) |
|---|---|---|
| `acs/v2/single_cell.py::SingleCellState` | `cell_id` | divide-event introduces a *child* `cell_id`; the *parent* `cell_id` either retires or is reassigned to one of the daughters. Relationship recorded out-of-schema (event log) — see §3 below. |
| `acs/v2/single_cell.py::SingleCellState` | `measurement_boundary` | parent boundary partitions into two daughter `MeasurementBoundary` polygons via a cleavage curve; cleavage curve geometry deferred. |
| `acs/v2/single_cell.py::SingleCellState` | `cell_state` | divide event transitions parent through a brief `dividing` literal (NOT YET in the `CellStateLabel` enum — see §4 schema-only TBD). |
| `acs/v2/focal_adhesion.py::FocalAdhesionState` | `cell_id` | each FA's `cell_id` reassigned to whichever daughter the FA's `position_um_xy` falls within. FAs straddling the cleavage curve: deferred (lifetime-end vs reassign-to-nearest is a Phase 1.5 decision). |
| `acs/v2/focal_adhesion.py::FocalAdhesionState` | `linked_protrusion_id` | preserved verbatim if the linked protrusion is still active and on the same side; otherwise transitively follow the protrusion's reassignment. |
| `acs/v2/protrusion.py::ProtrusionEvent` | `cell_id` | each protrusion's `cell_id` reassigned by the same partition rule; protrusions straddling the cleavage curve are *ended* by Phase 1.5 default (cleavage destroys local cytoskeletal continuity). |
| `acs/v2/protrusion.py::ProtrusionEvent` | `associated_adhesion_ids` | follows FA reassignment transitively. |
| `acs/v2/ecm_substrate.py::ECMSubstrateState` | all fields | unchanged by cytokinesis — the ECM is the substrate, not the cell. The cumulative `accumulated_traction_nNs_per_um2` at the cleavage site is *preserved*; if the cleavage produces an immediate displacement of FA traction sources, that is deferred to whatever Phase D scaffolding lands. |
| `acs/v2/cell_cluster.py::CellClusterState` | `cells` mapping | parent `cell_id` removed, two daughter `cell_id`s inserted; total `cells` count increases by 1. |

This table is descriptive. The exact partition rule (cleavage
curve geometry, FA reassignment policy, protrusion termination
policy) is Phase 1.5 scope.

---

## 2. Event-log representation (out-of-schema)

A cytokinesis event is recorded out-of-schema — neither
`SingleCellState` nor `CellClusterState` has a built-in division
history field in Phase 1. A Phase 1.5 commit that wants to
preserve the parent→daughters lineage can either:

- (a) extend `SingleCellState` with optional `parent_cell_id:
  Optional[str]` and `division_time_s: Optional[float]` fields
  (schema-additive, validate per existing
  `validate()` patterns), or
- (b) write the lineage to an external event-log artifact (e.g.,
  a new `CytokinesisEvent` schema in `acs/v2/cytokinesis.py`,
  parallel to `ProtrusionEvent`).

Phase 1 does not commit to either; both are valid future paths.

The note's recommendation (NOT a lock): option **(b)** — a separate
`CytokinesisEvent` schema mirrors the `ProtrusionEvent` precedent
and keeps `SingleCellState` lean. Option (a) creates a perpetual
parent-pointer tax on every cell record even if division never
fires.

---

## 3. What is *deferred* to Phase 1.5 / Phase 2

Each item below requires its own design-discussion round before
any code commit. Listed roughly in order of dependency.

### 3.1 Division trigger criterion (Phase 1.5)

What triggers `cell_state` transition to `dividing`? Candidates
(none of which are Phase 1 commitments):
- **Volume / area threshold**: `projected_area_um2() ≥ A_divide`.
- **Internal accumulated mechanical stress**: `Σ_FA |traction|`
  exceeds threshold.
- **Time-since-birth**: `time_s - parent.division_time_s ≥ τ_cycle`.
- **External signal** (Phase 2+, requires signaling layer).

Phase 1.5 lock decides one or a combination, with literature
reference (per Magic-Number Block: any threshold parameter must
be derivable, grid-invariant, and not chosen to fit a specific
result).

### 3.2 Cleavage curve geometry (Phase 1.5)

Given a parent `MeasurementBoundary` polygon, what curve
partitions it into two daughter polygons? Candidates:
- **Equal-area straight cut**: a chord that bisects the polygon
  area, perpendicular to the longest principal axis (PCA on
  vertices).
- **Centroid-aligned cut**: a chord through the centroid,
  perpendicular to a chosen direction (e.g., spindle direction
  in 3D, projected to 2D).
- **Constriction-mechanics cut**: a curve that follows the
  current cleavage-furrow position, derived from a constriction
  ring model (heaviest path, requires a separate force law).

Phase 1 schema accommodates any 2D cleavage curve as long as the
two daughter polygons each pass `MeasurementBoundary.validate()`.

### 3.3 FA / protrusion reassignment policy (Phase 1.5)

Given a cleavage curve, how are FAs and protrusions reassigned?

- **Position-based** (default candidate): FA whose
  `position_um_xy` falls strictly inside daughter A's boundary →
  daughter A's `cell_id`; same for daughter B; FAs *on* the
  cleavage curve: ended.
- **Fractional reassignment**: an FA on the curve is split into
  two with `bound_fraction` halved each — schema-additive change
  (might require relaxing the `unbound`-and-zero-traction
  invariant, careful).

Phase 1 schema fields support the position-based policy without
modification (FA `cell_id` is mutable through the schema's normal
`replace`-style construction).

### 3.4 Daughter-cell field initialization (Phase 1.5)

What does each daughter inherit from the parent?
- `maturity` / `bound_fraction` per FA: preserved.
- `linked_protrusion_id` per FA: follows protrusion reassignment.
- `time_s` per `SingleCellState`: parent's `time_s` (no clock
  reset).
- `cell_state`: typically `alive` (or `nascent` if Phase 1.5 adds
  that literal — see §4).
- ECM cumulative traction at the cleavage site: preserved
  (substrate state is independent of cell identity).

### 3.5 Cytokinesis force on the active contour (Phase 2)

The constriction-ring force that *drives* cytokinesis on the
active contour mechanics is explicitly Phase 2. Phase 1's active
contour (cortex + area + protrusion force + FA traction) does not
include a constriction term. Phase 1.5 may add a *kinematic*
cytokinesis (instantaneous boundary partition without a force
term); a *force-driven* cytokinesis is Phase 2.

---

## 4. Schema-only TBD: `dividing` `CellStateLabel` literal

The current `CellStateLabel` enum in `acs/v2/single_cell.py`
(per the existing schema) does NOT include a `dividing` literal.
Adding it is a schema-additive change that future executable
cytokinesis will need; it is also a candidate for inclusion *now*
as a schema-only field even before executable code lands, mirroring
the pattern already used for cell-cycle-dynamics fields ("Phase 1
Out / Schema-Only").

This note does NOT add the literal; it flags it as a candidate.
Decision deferred to the cytokinesis Phase 1.5 design round, or to
a separate schema-only PR if the parallel cell-cycle-dynamics
schema work needs the literal first.

---

## 5. Activation criteria (must be satisfied before executable code)

A future PR that *implements* cytokinesis must:

1. **Pass its own Sanity Gate** with the 6-item template applied
   to the chosen division-trigger + cleavage-curve + reassignment
   law, including:
   - §1 dimensional analysis on whatever quantity drives division
     (area μm², stress nN/μm², time s);
   - §2 boundary cases for empty cells, single-vertex cells,
     overlapping daughters;
   - §3 conservation: total cell count increases by exactly 1 per
     event; total FA count is conserved (or explicitly accounted
     for if FAs on the curve are ended); no cumulative ECM
     traction drift.
   - §6 measurement-protocol consistency: the division event must
     be detectable in the simulation output by the same modality
     PI's experimental data uses (typically, a single connected
     cell becomes two connected cells on top-down projection).
2. **Be design-discussion-locked** for the trigger criterion (§3.1)
   and cleavage-curve geometry (§3.2), each with literature
   reference per Magic-Number Block.
3. **Not introduce closed-loop ECM coupling** unless the
   closed-loop ECM gate (Phases A-E per
   `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`) is
   satisfied. A division event that triggers a remodeling response
   on the ECM bypasses Phase D/E entirely; that is not allowed.
4. **Preserve the v2 schema invariants** documented in §1: every
   schema's `validate()` must still pass on each daughter, every
   FA's traction contract holds, every protrusion's lifecycle
   field set is consistent.

If any of (1)–(4) is unmet, the executable cytokinesis PR halts
with a `status=blocker` to PI and three concrete options
(reduce trigger / use kinematic-only / split into smaller PRs).

---

## 6. What this note is *not*

- Not a Sanity Gate. The Sanity Gate is written by the future
  Phase 1.5 cytokinesis-executable PR.
- Not a design lock. Each deferred item (§3.1–§3.5) gets its own
  design-discussion round.
- Not a schema change. No field is added to any v2 schema by this
  note. The `dividing` `CellStateLabel` literal flagged in §4 is
  TBD, not added.
- Not a code commitment. No `cytokinesis_step` function is
  proposed; no module path is committed.
- Not authoritative for biological parameters. Reference values
  (cell-cycle duration, division area threshold, etc.) require
  literature-first extraction per `docs/v2_phase1_forward_roadmap.md`
  Magic-Number Block.

---

## 7. References

- Forward roadmap: `docs/v2_phase1_forward_roadmap.md` lines 51,
  53–61 (Phase 1 Out / Schema-Only).
- v2 schemas (current state):
  - `acs/v2/single_cell.py::SingleCellState`
  - `acs/v2/focal_adhesion.py::FocalAdhesionState`
  - `acs/v2/protrusion.py::ProtrusionEvent`
  - `acs/v2/ecm_substrate.py::ECMSubstrateState`
  - `acs/v2/cell_cluster.py::CellClusterState`
- 6.3a static FA preflight (FA `cell_id` mutability precedent):
  `docs/v2_focal_adhesion_dynamics_sanity_gate.md`.
- 6.3b protrusion-coupled FA dynamics (link semantics precedent):
  `docs/v2_63b_protrusion_coupling_locked.md`.
- Closed-loop ECM gate phased plan (deferred-blocker pattern
  precedent): `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`.
- Memory rules:
  - `rule10_unit_derivation_in_docs.md` — applies to any cell
    volume / cytokinesis duration / division-area unit comparison
    in this note's deferred sections.
  - `hard_rule_11_wording_boundary_meta_test.md` — applies if any
    future commit treats this design note as gate-item satisfaction.
