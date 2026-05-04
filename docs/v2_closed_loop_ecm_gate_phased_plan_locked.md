# V2 Phase 1 Closed-Loop ECM Gate — Phased Plan (LOCKED)

**Date**: 2026-05-04 KST
**Authors**: Claude + Codex design-discussion (3-round adversarial lock,
PI id=809 aggressive debate posture, PI id=1008/1057 autonomy)
**Source unit**: design-discussion `topic=v2-layer-2-closed-loop-ecm-gate`,
MCP id 1223-1230
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for Phase A immediate evidence + Phase D/E blocker resolution
sequencing.

---

## 0. Scope

This document phases the closed-loop ECM gate (forward roadmap §6.4
"Closed-Loop ECM Gate") into 5 sub-phases (A-E). It identifies which
gate items are provable now (open-loop side) vs which require
constitutive-law / interface design lock first (closed-loop side).

**Hard discipline**: open-loop evidence does NOT satisfy closed-loop
gate items. Hard Rule 11 (measured value matches claim) requires we
NOT advertise open-loop accumulator monotonicity as Item 1 response
monotonicity.

---

## 1. Phase A-E Final Lock

### Phase A — Immediate evidence note (0 design blocker)

**Scope**: Item 3 (no response on zero traction) open-loop evidence +
Item 6 (Rule 10 unit chain) open-loop side evidence.

- Item 3: trivial PASS (zero traction → zero accumulator delta), tests
  already in `tests/test_v2_ecm_open_loop.py`. Evidence note cites
  existing tests.
- Item 6: open-loop unit chain `[nN/μm²]·[s] = [nN·s/μm²]` already in
  `acs/v2/dynamics/ecm_open_loop.py` module docstring. Evidence note
  cites docstring.

**Wording lock (Hard Rule 11 protection)**:
- "Item 3 trivial PASS, open-loop side"
- "Item 6 open-loop unit chain documented"
- **NOT** "Items 3, 6 of closed-loop gate satisfied"

**Deliverable**: small docs commit, e.g.
`docs/v2_closed_loop_ecm_gate_phase_a_evidence.md` cross-referencing
existing tests/docstrings. Code commit 0건.

### Phase B — Stimulus monotonicity precursor tests (0 design blocker)

**Scope**: prescribed-traction stimulus accumulator monotonicity
regression tests. Open-loop side only.

**Wording lock (Hard Rule 11 protection)**:
- "Item 1 precursor evidence: prescribed-traction stimulus monotonicity"
- **NOT** "Item 1 satisfied", **NOT** "ECM response monotone"
- Item 1 remains UNSATISFIED until response field + constitutive
  direction locked (Phase E)

**Reasoning**: `accumulated_traction_nNs_per_um2` is a stimulus/memory
field, not an ECM response field. Measuring its monotonicity does not
measure remodeling response monotonicity (Hard Rule 11 modality
matching).

**Deliverable**: tests-only commit in `tests/test_v2_ecm_open_loop.py`
or new `tests/test_v2_closed_loop_ecm_phase_b_precursor.py`.

### Phase C — Open-loop sweep harness baseline (0 design blocker)

**Scope**: open-loop sensitivity sweep over `(grid_spacing, dt_ecm)`
across the 4 ECM-OL preflight functions (traction, stiffness, density,
orientation). New script + Sanity Gate.

**Wording lock**:
- "Open-loop sweep harness and baseline sensitivity evidence"
- **NOT** "Item 5 satisfied"
- Item 5 remains OPEN for closed-loop side (closed-loop grid/dt
  sensitivity once scattering + response law exist)

**Deliverable**: new `scripts/run_ecm_ol_sensitivity_sweep.py` +
`runs/<UTC>_ecm_ol_sensitivity/` directory (PI visible deliverable
pattern, P1 alpha / ECM-OL precedents).

### Phase D — No-op closed-loop scaffolding (BLOCKED on #3 + #4)

**Scope**: Schema/API hooks for FA→ECM and ECM→FA paths, default OFF
or identity, no behavior change.

**Entry blocker**: Phase D code (`step_fa_to_ecm_response`,
`step_ecm_to_fa_bias` function names are API commitments) cannot land
until blockers #3 (FA→ECM scattering geometry interface) and #4
(ECM→FA bias target interface) are design-locked.

Blockers #1/#2/#5 can remain unresolved for active constitutive
behavior, but Phase D no-op interface MUST have #3/#4 resolved or it
bakes the wrong abstraction.

**Deliverable**: code commit only after design lock for #3 + #4.
3 new Sanity Gates per function.

### Phase E — Active closed-loop response (BLOCKED on all 5 + effective_stiffness law)

**Scope**: Items 1, 2, 5, 6 closed-loop side + Lyapunov-like metric
(Item 4 candidate).

**Entry blocker**: All 5 hard blockers + effective_stiffness law
decision (if used).

**Deliverable**: full closed-loop ECM remodeling implementation.
Constitutive law directly determines update equations.

---

## 2. Hard Blockers (5)

Each requires its own design round before the corresponding gate item
can be claimed satisfied.

### Blocker #1 — Constitutive-law direction (Item 1 closed-loop)

Question: which schema field(s) respond to traction stimulus, in what
direction, with what saturation form?

Default direction proposed (deferred to #1 round, see §4):
- FA traction → accumulator/memory only (Item 3 trivial preserved)
- Remodeling response (in Phase E):
  - `fiber_density`: bounded asymptotic toward 1 under traction exposure
    (literature-backed)
  - `orientation_tensor`: bounded asymptotic alignment toward scattering
    direction (after #3 lock)
- raw `stiffness_kpa` unchanged in first closed-loop law

**Resolution**: separate adversarial design round with literature-first
table (Hall 2016 / Trichet 2012 / Stylianopoulos 2018 seed, see §7).

### Blocker #2 — Saturation form (Item 2 closed-loop)

Question: how does response saturate?

Options:
- (a) Hard cap: `field = min(field + delta, max_value)` — looks like hidden clamping
- (b) Asymptotic-to-bound: `field += alpha · (max_value - field) · stimulus` — derived form

**Lock direction**: asymptotic-to-bound for any bounded schema field
(`density ∈ [0,1]`, orientation bound). Hard caps look like hidden
clamping unless derived. Final form per-field decision in #2 round.

**Resolution**: separate design round, may merge with #1 round if
constitutive direction commits to specific bounded fields.

### Blocker #3 — FA→ECM scattering geometry (Phase D entry)

Question: 1 FA at position `(x, y)` → which ECM grid cells receive its
traction signal? With what kernel?

Options:
- (a) Nearest-neighbor: 1 FA → 1 ECM cell
- (b) Bilinear: 1 FA → 4 surrounding ECM cells, weighted by distance
- (c) Gaussian kernel: 1 FA → N ECM cells within radius, weighted by
  Gaussian
- (d) Other (literature-derived)

**Resolution**: separate design round. Interface lock (input/output
shapes, units, sign conventions) required for Phase D no-op
scaffolding.

**Phase D interface contract** (must be defined even if default
implementation is identity):
- Inputs: `(adhesions, ecm_state, scattering_params)`
- Outputs: `traction_density_nN_per_um2[grid_shape]` (units locked)
- Default no-op: returns zero field (identity)

### Blocker #4 — ECM→FA bias target (Phase D entry + side-door risk)

Question: ECM state biases FA dynamics how?

Options:
- (a) `traction_scale_nN` scale (FA traction magnitude depends on local
  stiffness)
- (b) Rate multiplier (FA `k_maturity_per_s` boosted by local stiffness)
- (c) FA nucleation rate (more FAs in stiff regions)
- (d) Combination

**Resolution**: separate design round. Interface lock required for
Phase D no-op scaffolding.

**Phase D interface contract** (must be defined even if default
implementation is neutral):
- Inputs: `(adhesions, ecm_state, bias_params)`
- Outputs: per-FA bias multipliers (typed dict, similar to 6.3b
  `ProtrusionStateMultipliers` pattern)
- Default no-op: returns 1.0 multiplier for all FAs (neutral)

**Side-door bias risk** (see §3 effective_stiffness guard): even
read-only helpers can implement ECM→FA bias if they consume geometry
fields. Must be gated by this blocker.

### Blocker #5 — Lyapunov-like metric (Item 4 closed-loop)

Question: bounded-feedback metric for single-cell loop stability?

Candidates (deferred to #5 round, see §5):
1. Total accumulated traction energy
2. Total fiber alignment magnitude
3. Cell-substrate work rate

**Resolution**: separate design round, MUST be after #1/#2/#3/#4
locked (measurement-protocol-dependent on what loop actually updates).

---

## 3. Effective Stiffness Helper Guard

`ECMSubstrateState.stiffness_kpa` is **raw/background substrate
stiffness**. Closed-loop ECM remodeling (first pass) does NOT mutate
raw `stiffness_kpa`.

**Effective stiffness helper rules**:
- Phase A-C/D: NO `effective_stiffness()` helper unless it returns raw
  `stiffness_kpa` exactly (identity).
- Once helper consumes `fiber_density` or `orientation_tensor`, it
  becomes ECM→FA mechanosensing law and MUST satisfy #4 interface
  lock + Rule 10/sign/saturation gate.
- Schema split 0건 (`stiffness_kpa` retains raw meaning, no
  `effective_stiffness_kpa` stored field).

**Reason**: even read-only helpers, plugged into FA bias path, become
ECM→FA bias by another name. Side-door risk.

---

## 4. Constitutive Direction Seed (deferred to #1 round)

Opening default position for the upcoming #1 design round (literature
table required before code):

1. **Stimulus memory ↔ remodeling response separation**:
   `accumulated_traction_nNs_per_um2` = stimulus memory (current).
   Remodeling response = different field(s).
2. **First closed-loop response targets collagen geometry**:
   `fiber_density` and/or `orientation_tensor` before stiffness.
3. **Raw `stiffness_kpa` unchanged in first closed-loop law**.
4. **Asymptotic-to-bound saturation** for bounded fields (NOT hard
   clipping). Derived form required.

This is the SEED for the #1 design round, NOT a lock. Literature
extraction first.

---

## 5. Lyapunov Metric Seed (deferred to #5 round)

Candidate metrics for the upcoming #5 design round (selection
measurement-protocol-dependent on #1-#4 locks):

1. **Total accumulated traction energy**:
   `0.5 · Σ_grid |T|² · grid_volume` — positive definite, decreasing
   under relaxation if response = relaxation.
2. **Total fiber alignment magnitude**:
   `Σ_grid Tr(Q²)^0.5 · grid_volume` — bounded by `N_grid · max(|Q|)`.
3. **Cell-substrate work rate**:
   `Σ_FA F_i · (v_substrate_at_FA - v_cell_at_FA)` — sign decision +
   boundedness depends on FA contract.

This is the SEED for the #5 round, NOT a lock.

---

## 6. Sanity Gate Matrix per Phase

| Phase | §1 Dim | §2 Bound | §3 Cons | §4 Num | §5 Sign | §6 MeasProto |
|---|---|---|---|---|---|---|
| A (Item 3/6 evidence) | inherit ECM-OL | inherit | inherit | N/A (no new code) | inherit | wording: "open-loop side" |
| B (stimulus monotonicity) | inherit | inherit | tests preserve order | inherit | inherit | "precursor not satisfaction" |
| C (open-loop sweep) | inherit | grid/dt boundary | sweep range bounded | sweep dt validation | inherit | "open-loop sweep" labeling |
| D (no-op scaffolding) | per-function units | per-function | per-function | per-function | per-function | per-function (#3/#4 interfaces) |
| E (active closed-loop) | full constitutive law | full saturation | full Lyapunov | full feedback dt | full sign/sense | full Item 1/2/4/5/6 closed-loop |

---

## 7. Reference Biology Table (literature seed for #1 round, NOT lock)

| Reference | Year | Venue | Finding | Relevance |
|---|---|---|---|---|
| Hall MS et al. | 2016 | PNAS | Cell traction induces collagen fiber alignment; stiffness perception is secondary | Geometry-first response ✓ |
| Trichet L et al. | 2012 | PNAS | Substrate stiffness as cell INPUT (sensing), not output of cell traction | Raw stiffness unchanged ✓ |
| Stylianopoulos T et al. | 2018 | Nat Rev Cancer | Tumor ECM remodeling: collagen alignment time-scale faster than stiffness change | Time-scale separation ✓ |
| Eichinger J et al. | 2020 | Soft Matter | Computational fiber-network ECM remodeling | Reference simulation lit |
| Notbohm J et al. | 2015 | TBD | Long-range force transmission in collagen | Scattering geometry hint |

Each reference must be verified before parameterization (per
`docs/v2_phase1_plan_consolidated.md` Magic-Number Block discipline).

---

## 8. Test Catalog per Phase

### Phase A
- No new tests; cite existing
  `tests/test_v2_ecm_open_loop.py::test_zero_traction_no_op`.

### Phase B
- `test_stimulus_accumulator_monotone_under_positive_traction`
- `test_stimulus_accumulator_no_change_when_traction_zero`
- `test_stimulus_accumulator_signed_traction_can_decrease`
- `test_stimulus_monotonicity_does_not_satisfy_response_item_1` (Hard Rule 11 wording check)

### Phase C
- `test_open_loop_sweep_harness_runs_4_channels`
- `test_open_loop_sweep_grid_spacing_variation`
- `test_open_loop_sweep_dt_variation`
- `test_open_loop_sweep_summary_reports_baseline_only` (Hard Rule 11
  wording check: "baseline sensitivity, NOT closed-loop sensitivity")

### Phase D (TBD per #3/#4 lock)
- Per-function no-op tests (interface validation, identity behavior)

### Phase E (TBD per all blockers lock)
- Full closed-loop ECM remodeling tests

---

## 9. Cross-Room Dispatch & Phase A Immediate Go-List

This file is the design-team input to implementation-work for:

### Phase A immediate go-list (no design blocker)
1. Write `docs/v2_closed_loop_ecm_gate_phase_a_evidence.md` citing:
   - `tests/test_v2_ecm_open_loop.py::test_zero_traction_no_op` (Item 3 trivial PASS)
   - `acs/v2/dynamics/ecm_open_loop.py` module docstring `[nN/μm²]·[s] = [nN·s/μm²]` (Item 6 open-loop side)
2. Wording strictly "open-loop side" / "trivial PASS, NOT gate item satisfaction"
3. Small docs commit, no code changes

### Phase B follow-up (no design blocker)
1. Add tests for stimulus monotonicity in `tests/test_v2_ecm_open_loop.py` or new
   `tests/test_v2_closed_loop_ecm_phase_b_precursor.py`
2. Hard Rule 11 wording: "precursor evidence, NOT Item 1 satisfaction"

### Phase C follow-up (no design blocker)
1. New `scripts/run_ecm_ol_sensitivity_sweep.py`
2. Sweep `(grid_spacing, dt_ecm)` over 4 ECM-OL preflight functions
3. Persistent artifact: `runs/<UTC>_ecm_ol_sensitivity/` (PI visible deliverable
   pattern, P1 alpha / ECM-OL precedents)

### Phase D blocker resolution (design-discussion next rounds)
1. #3 FA→ECM scattering geometry interface (separate design round)
2. #4 ECM→FA bias target interface (separate design round)
3. Phase D code waits for both #3 + #4 interface lock

### Phase E blocker resolution (design-discussion later rounds)
1. #1 constitutive-law direction (literature-first table, separate round)
2. #2 saturation form (may merge with #1)
3. #5 Lyapunov metric (after #1-#4 lock)

### Cadence rule
impl-work in-progress reports per `docs/implementation_workflow.md` §8
(promised cadence even when no change).

### Sequencing rule
design-first lock → cross-room dispatch → impl entry (per PI id=1088/1089
sequencing fix).

Rounds 1-3 of the closed-loop ECM gate phased plan lock are MCP id
1223-1230.
