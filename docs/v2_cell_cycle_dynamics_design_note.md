# V2 Phase 1 — Cell-Cycle Dynamics Design Note (schema-only, NOT executable)

**Status**: design note only, **NOT a Sanity Gate, NOT a code
commitment**. Per
`docs/v2_phase1_forward_roadmap.md` ("Cell-cycle dynamics:
schema-only" entry in the "Phase 1 Out / Schema-Only" section).

**Author**: implementation-work Claude, drafted as a parallel idle
unit alongside the cytokinesis design note (`ef2806c` /
`c8ec9cc` / `3be6fa1`). PI ratify status: full delegation per PI
id=939/1008.

**Hard contract**: nothing in this note authorizes a cell-cycle
dynamics commit. Phase 1 explicitly defers cell-cycle *dynamics*
to Phase 1.5 / Phase 2 while keeping the *schema* fields
(`cell_age_s`, `cell_cycle_phase`, `mechanosignal_yap_taz`)
already present in `acs/v2/single_cell.py` for downstream
validation and persistence. Any future code that *implements*
cell-cycle dynamics must write its own Sanity Gate document and
design-discussion lock.

---

## 0. Why this design note exists

The forward roadmap classifies cell-cycle dynamics as
**schema-only** for Phase 1: the schema carries the fields, but
no dynamics step writes to them. Without this note, the existing
schema fields (`cell_age_s`, `cell_cycle_phase`,
`mechanosignal_yap_taz`) carry no documented contract for what a
future cell-cycle dynamics layer must do, and a future Phase 1.5
PR would re-litigate questions already constrained by the schema
choices.

This note documents:
- the existing schema fields and their validators;
- which dynamics are deferred and require their own
  design-discussion round before executable code lands;
- the relationship to cytokinesis (cell-cycle phase M is the
  natural cytokinesis trigger *candidate*; exact M↔cytokinesis
  coupling deferred to Phase 1.5 design-lock — see §2.5);
- activation criteria for future cell-cycle dynamics.

It does **not** propose constitutive equations for phase
transitions, age-increment timing, or YAP/TAZ feedback. Those are
explicitly Phase 1.5+ scope.

---

## 1. Existing schema (Phase 1, schema-only)

Verified directly via `grep` on `acs/v2/single_cell.py`:

| Field | Type | Default | Validator |
|---|---|---|---|
| `cell_age_s` | `float` | `0.0` | finite and non-negative |
| `cell_cycle_phase` | `Optional[CellCyclePhase]` | `None` | when provided, must be one of `{"G0", "G1", "S", "G2", "M"}` |
| `mechanosignal_yap_taz` | `Optional[float]` | `None` | when provided, must be in `[0, 1]` |
| `division_count` | `int` | `0` | non-negative integer (also used by cytokinesis design note §2) |
| `parent_cell_id` | `Optional[str]` | `None` | non-empty when provided, must differ from `cell_id` |

The `CellCyclePhase` literal is the canonical Phase 1 phase
enumeration:

```python
CellCyclePhase = Literal["G0", "G1", "S", "G2", "M"]
```

Together these five fields form the Phase 1 cell-cycle schema.
**No Phase 1 code mutates them as part of a dynamics step**; they
are caller-supplied (e.g., from imaging metadata) or set at cell
construction and remain constant.

---

## 2. What is *deferred* to Phase 1.5 / Phase 2

Each item below requires its own design-discussion round before
any code commit. Listed roughly in order of dependency.

### 2.1 Phase-transition dynamics (Phase 1.5)

What governs the `cell_cycle_phase` transitions
`G0 → G1 → S → G2 → M`? Candidates (none of which are Phase 1
commitments):

- **Time-based**: each phase has a duration `τ_<phase>`
  (literature-derived, e.g., MCF7 has total cycle ~24 h with G1
  ≈ 11 h, S ≈ 8 h, G2 ≈ 4 h, M ≈ 1 h per typical references).
  Transition fires when `cell_age_s_in_phase` crosses `τ_<phase>`.
- **Mechanosensing-modulated**: phase durations scale with local
  ECM stiffness or `mechanosignal_yap_taz` value (literature-
  documented: stiff substrate → faster cycling).
- **Stochastic / hazard-based**: each phase has a hazard rate;
  transition is sampled each step. Requires RNG + reproducibility
  contract.
- **External signal driven** (Phase 2+): full signaling layer
  (Notch / Wnt / TGF-β feedback) — out of scope for Phase 1.5.

Phase 1.5 lock decides one or a combination, with literature
reference (per Magic-Number Block: any threshold parameter must
be derivable, grid-invariant, and not chosen to fit a specific
result).

### 2.2 Age increment timing (Phase 1.5)

Does `cell_age_s` advance every simulation step (`+= dt_cell_s`)
or only on dynamics-relevant ticks (e.g., per minute of
simulation time)? Candidates:

- **Per-step increment**: `cell_age_s += dt_cell_s` each call to
  the cell step. Couples cell-cycle clock to integration step
  size; care needed if `dt_cell_s` is variable.
- **Wall-clock-style increment**: `cell_age_s = time_s -
  birth_time_s` derived on demand from a separate `birth_time_s`
  field (would require a schema-additive change adding
  `birth_time_s: Optional[float]` next to `cell_age_s`).

Phase 1.5 lock decides; either path keeps the existing
`cell_age_s` field as the surface contract.

### 2.3 YAP/TAZ feedback dynamics (Phase 2)

The `mechanosignal_yap_taz` field is schema-only in Phase 1. A
Phase 2 dynamics layer could:

- compute `mechanosignal_yap_taz` from local ECM stiffness or
  cumulative traction stimulus (open-loop), or
- close the loop: `mechanosignal_yap_taz` → cell-cycle phase
  duration scaling → division rate → traction generation
  feedback.

Both are explicitly out of Phase 1. The closed-loop form (second
bullet) requires the closed-loop ECM gate
(`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`) to be
satisfied first, plus its own per-function Sanity Gate.

### 2.4 G0 / quiescence semantics (Phase 1.5)

Does a cell that enters `G0` exit it on its own, or only under
external signal? Phase 1 schema simply allows `G0` as a valid
phase value; Phase 1.5 dynamics must decide whether the
phase-transition rules above apply to G0 (likely no, by
biological convention) and what triggers re-entry to G1.

### 2.5 M-phase to cytokinesis link (Phase 1.5)

The cell-cycle phase `M` is the mitotic phase that includes
cytokinesis. The exact relationship between `cell_cycle_phase`
transitioning out of `M` and the cytokinesis event firing
(per `docs/v2_cytokinesis_design_note.md`) needs to be locked:

- **Tight coupling**: cytokinesis event fires on `M → G1`
  transition. Two daughters share birth time.
- **Loose coupling**: cytokinesis event fires when its own
  trigger criterion is met (`docs/v2_cytokinesis_design_note.md`
  §3.1), possibly during M-phase, possibly later. The cell-cycle
  layer marks `M`-end independently of cytokinesis firing.

Phase 1 schema allows either; the Phase 1.5 lock (which encompasses
both cytokinesis and cell-cycle dynamics) must pick one.

---

## 3. Relationship to cytokinesis design note

`docs/v2_cytokinesis_design_note.md` is the paired Phase 1
design-only deliverable for cell division. The two notes are
linked:

- Cell-cycle phase `M` is the natural trigger candidate for the
  cytokinesis event (§2.5 above + cytokinesis note §3.1).
- The `division_count` and `parent_cell_id` schema fields are
  shared lineage breadcrumbs (cytokinesis note §2 documents how
  they record post-event lineage identity; cell-cycle note §1
  lists them as part of the schema surface).
- A Phase 1.5 PR that implements cell-cycle dynamics likely lands
  alongside the cytokinesis-execution PR if Phase 1.5 picks tight
  M↔cytokinesis coupling (one of two candidates per §2.5).
  Loose coupling is also a valid Phase 1.5 candidate; either path
  must respect both notes' activation criteria.

**No Phase 1 code couples the two layers**. Either future
implementation must respect both notes' activation criteria.

---

## 4. Activation criteria (must be satisfied before executable code)

A future PR that *implements* cell-cycle dynamics must:

1. **Pass its own Sanity Gate** with the 6-item template applied
   to the chosen phase-transition + age-increment + YAP/TAZ
   dynamics, including:
   - §1 dimensional analysis on `cell_age_s` (s) and any phase
     duration parameters (s) and YAP/TAZ feedback (dimensionless);
   - §2 boundary cases: zero-age cell, `cell_age_s == τ_phase`
     boundary, missing `cell_cycle_phase` value (None
     transitions), `G0` re-entry semantics;
   - §3 conservation: total cell-cycle time accumulated across
     phases must match `cell_age_s` (or the documented
     birth-time-derived alternative);
   - §4 numerical: `dt_cell_s · (1/τ_phase) ≤ safety margin` for
     time-based transitions; the existing
     `_DT_RATE_SAFETY_MARGIN = 0.5` is the obvious reuse target;
   - §6 measurement-protocol consistency: phase transitions
     should be detectable in simulation output via the same
     modality experimental imaging uses (e.g., M-phase via mitotic
     marker fluorescence, cell area change at division).
2. **Be design-discussion-locked** for the phase-transition
   criterion (§2.1), age-increment timing (§2.2), and M→cytokinesis
   coupling (§2.5), each with literature reference per
   Magic-Number Block.
3. **Not introduce closed-loop ECM coupling** unless the
   closed-loop ECM gate (Phases A-E per
   `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`) is
   satisfied. A YAP/TAZ feedback that maps cumulative ECM traction
   to phase duration scaling is a closed-loop coupling and is
   blocked behind Phase D / E.
4. **Preserve the existing schema invariants**:
   - `cell_age_s` finite and non-negative
   - `cell_cycle_phase` in `{G0, G1, S, G2, M}` when provided
   - `mechanosignal_yap_taz` in `[0, 1]` when provided
   - `division_count` non-negative integer
   - `parent_cell_id` non-empty and ≠ `cell_id` when provided
5. **Not modify the schema** unless the modification is
   schema-additive (adding a field, never removing or changing
   the type of an existing field) — Phase 1.5 may add
   `birth_time_s` per §2.2 but cannot remove `cell_age_s`.

If any of (1)–(5) is unmet, the executable cell-cycle dynamics PR
halts with a `status=blocker` to PI and three concrete options
(reduce dynamics scope / add a missing Sanity Gate item / split
into smaller PRs).

---

## 5. What this note is *not*

- Not a Sanity Gate. The Sanity Gate is written by the future
  Phase 1.5 cell-cycle-dynamics-executable PR.
- Not a design lock. Each deferred item (§2.1–§2.5) gets its own
  design-discussion round.
- Not a schema change. No field is added to any v2 schema by this
  note. The existing fields (`cell_age_s`, `cell_cycle_phase`,
  `mechanosignal_yap_taz`, `division_count`, `parent_cell_id`)
  are documented but unchanged.
- Not a code commitment. No `cell_cycle_step` function is
  proposed; no module path is committed.
- Not authoritative for biological parameters. Reference values
  (cell-cycle phase durations, YAP/TAZ-stiffness coupling
  coefficients) require literature-first extraction per
  `docs/v2_phase1_forward_roadmap.md` Magic-Number Block.

---

## 6. References

- Forward roadmap: `docs/v2_phase1_forward_roadmap.md`
  ("Cell-cycle dynamics: schema-only" entry in the "Phase 1 Out
  / Schema-Only" section).
- Existing schema: `acs/v2/single_cell.py` (`cell_age_s`,
  `cell_cycle_phase`, `CellCyclePhase` literal,
  `mechanosignal_yap_taz`, `division_count`, `parent_cell_id`).
- Paired design note: `docs/v2_cytokinesis_design_note.md`
  (M-phase is the natural cytokinesis trigger candidate; exact
  M↔cytokinesis coupling deferred to Phase 1.5 design-lock).
- Closed-loop ECM gate phased plan (gates YAP/TAZ feedback when
  it goes closed-loop):
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`.
- Magic-Number Block discipline: `CLAUDE.md` Magic-Number Block
  + `docs/v2_phase1_forward_roadmap.md`.
- Memory rules:
  - `rule10_unit_derivation_in_docs.md` — applies to any
    `cell_age_s` / phase duration / YAP/TAZ unit comparison in
    this note's deferred sections.
  - `hard_rule_11_wording_boundary_meta_test.md` — applies if any
    future commit treats this design note as gate-item or
    Sanity-Gate satisfaction.
