# 6.3a/6.3b FocalAdhesionDynamicsResult Typed Schema — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-05 KST
**Status**: design-discussion brief, **opening position only —
NOT a lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger)
precedent established by Hard Blocker #3 / #4 / Phase D no-op
scaffolding (all locked + Sanity Gated + implemented + Codex
PASS in this session).
**Source**: sister-gate-mirror audit Finding B1 from
`docs/v2_sister_gate_mirror_audit_2026_05_05.md` (commit
`d19670d`); Codex routing `id=1432` directing the B1 finding
into a separate cycle (NOT an audit commit, NOT deferred to
Phase E).
**Rule activation**: first design unit invoked under the
freshly-promoted Step 6 sister-gate-mirror memory rule
(`design_note_pre_commit_batch.md`, 6-step, promoted per Codex
`id=1428`).

**Hard contract**: this brief is the **opening position** for
the adversarial design round on the typed schema discipline +
immutability of `FocalAdhesionDynamicsResult` + its `diagnostics`
field for 6.3a (`step_focal_adhesions_static`) and 6.3b
(`step_protrusion_coupled_focal_adhesions`). It does NOT lock
the typed shape, does NOT commit the implementation, and does
NOT authorize any consumer migration plan. The locked phased
plan §1 Phase D is unchanged; Phase E remains BLOCKED on all 5
Hard Blockers + the effective_stiffness law decision.

---

## 0. Why this brief exists + Codex `id=1432` nuance

Sister-gate-mirror audit Finding B1 surfaced two closely-related
defects in `acs/v2/dynamics/focal_adhesion.py:153`:

```python
@dataclass
class FocalAdhesionDynamicsResult:
    """Outputs of one 6.3a static FA dynamics step."""

    updated_adhesions: tuple[FocalAdhesionState, ...]
    cell_force_nN_xy: np.ndarray
    substrate_reaction_nN_xy: np.ndarray
    radial_components: np.ndarray
    tangential_components: np.ndarray
    diagnostics: dict = field(default_factory=dict)
```

Two issues:
1. `@dataclass` plain — NOT `frozen=True`, NOT `slots=True`. A
   downstream caller can mutate `result.cell_force_nN_xy = ...`
   etc. Sister-pattern divergence from Phase D / HB#4 dataclass
   discipline (all 4 newer dataclasses use `frozen=True,
   slots=True`).
2. `diagnostics: dict = field(default_factory=dict)` is a loose
   untyped dict — exact regression Codex caught in Phase D round
   1 C1 (`id=1395`). Test consumers depend on specific dict keys
   (e.g., `result.diagnostics["aggregate_cell_force_nN_xy"]` in
   `tests/test_v2_focal_adhesion_dynamics.py:323-324`,
   `result.diagnostics["reciprocal_missing"]` in
   `tests/test_v2_protrusion_coupled_focal_adhesion.py:287-288`),
   defeating typed-key precedent.

**Codex `id=1432` nuance preservation** (critical to brief
framing): HB#4 still has its own `diagnostics_dict` surface
(`ECMToFABiasResult.diagnostics_dict: dict[str, float | int |
str]`). The B1 design question is **NOT** "ban all dict
diagnostics" but rather:

> Whether older 6.3a/6.3b result diagnostics need typed schema
> discipline and frozen/slots immutability before Phase E
> composition.

Phase D `FAToECMResponseResult` rejected a NEW loose diagnostics
field at the wrapper layer because the underlying primitives
already have typed channels. HB#4's existing `diagnostics_dict`
predates that rejection and is grandfathered. The 6.3a/6.3b
question is whether to grandfather similarly or bring forward the
typed discipline.

---

## 1. What is already locked (do not redebate here)

- 6.3a `step_focal_adhesions_static` consumer surface
  (`tests/test_v2_focal_adhesion_dynamics.py:323-324`) reads:
  - `result.diagnostics["aggregate_cell_force_nN_xy"]` (tuple of
    2 floats, summed across FAs)
  - `result.diagnostics["aggregate_substrate_reaction_nN_xy"]`
    (tuple of 2 floats, summed across FAs)
  - The keys are constructed in
    `acs/v2/dynamics/focal_adhesion.py` (6.3a).
- 6.3b `step_protrusion_coupled_focal_adhesions` consumer surface
  (`tests/test_v2_protrusion_coupled_focal_adhesion.py:287-288`)
  reads:
  - `result.diagnostics["reciprocal_missing"]` (int)
  - `result.diagnostics["linked_missing"]` (int)
- 6.3b producer
  (`acs/v2/dynamics/protrusion_coupled_focal_adhesion.py:425-434`)
  writes the 6.3a-shared keys plus four 6.3b-extension keys:
  `linked_missing`, `reciprocal_missing`,
  `multiplier_histogram`, `max_effective_rate_per_name`.
- Both 6.3a and 6.3b return the same `FocalAdhesionDynamicsResult`
  type today.
- HB#4 `ECMToFABiasResult.diagnostics_dict` is grandfathered (not
  in scope of this brief).
- Phase D `FAToECMResponseResult` does NOT have a diagnostics
  field (locked round 1 C1).

---

## 2. The four design questions

### Q1 — Typed schema for the 6.3a-shared diagnostics

Current shape: `diagnostics: dict` with **four** 6.3a-shared keys
(verified against both producers
`acs/v2/dynamics/focal_adhesion.py:381-390` and
`acs/v2/dynamics/protrusion_coupled_focal_adhesion.py:422-430`):

- `n_adhesions: int` (count of FAs in the step)
- `aggregate_cell_force_nN_xy: tuple[float, float]` (component
  sum across FAs, in `nN`)
- `aggregate_substrate_reaction_nN_xy: tuple[float, float]`
  (component sum across FAs, in `nN`)
- `max_traction_magnitude_nN: float` (max of L2 norm of per-FA
  cell force, in `nN`; `0.0` if `n_adhesions == 0`)

Two candidate shapes:

- **(a) Typed `FocalAdhesionDynamicsDiagnostics` dataclass**:
  ```python
  @dataclass(frozen=True, slots=True)
  class FocalAdhesionDynamicsDiagnostics:
      aggregate_cell_force_nN_xy: tuple[float, float]
      aggregate_substrate_reaction_nN_xy: tuple[float, float]
      max_traction_magnitude_nN: float
  ```
  Mirrors HB#4 typed-dataclass precedent. Each field has
  documented unit + shape. Type-check-friendly.
- **(b) Stay with `dict` but require `dict[str, float | int |
  tuple]`** (HB#4 grandfathered pattern). Adds typing on the
  values but keeps the per-call-site key schema flexibility.

### Q2 — 6.3b extension shape

Current shape: `diagnostics: dict` with the 4 6.3a-shared keys
above plus 4 6.3b-extension keys (verified against
`acs/v2/dynamics/protrusion_coupled_focal_adhesion.py:431-434`):

- `linked_missing: int` (always 0 in current implementation;
  reserved for future linkage-resolution paths)
- `reciprocal_missing: int` (count of FAs with resolved
  `linked_protrusion_id` whose linked protrusion did not list
  the FA in `associated_adhesion_ids`)
- `multiplier_histogram: dict[str, collections.Counter[float]]`
  (per-rate-name `Counter` of effective multiplier values rounded
  to 12 decimals for float64-round-off-stable comparison)
- `max_effective_rate_per_name: dict[str, float]` (per-rate-name
  maximum effective rate across all FAs, in `[1/s]`)

Three candidate shapes (only relevant if Q1 picks (a)):

- **(a) `ProtrusionCoupledDynamicsDiagnostics(FocalAdhesionDynamicsDiagnostics)`
  subclass**: 6.3b inherits from 6.3a base + adds 4 extension
  fields. Shared field set (3) lives in the base. Inheritance
  pattern.
- **(b) Separate `ProtrusionCoupledDynamicsResult` dataclass**
  distinct from `FocalAdhesionDynamicsResult`: 6.3a returns
  `FocalAdhesionDynamicsResult`; 6.3b returns a different result
  type whose `diagnostics` field is
  `ProtrusionCoupledDynamicsDiagnostics` with all 7 keys.
  Composition + isolation pattern. Mirrors HB#3 / HB#4 / Phase D
  approach where each phase has its own typed dataclass.
- **(c) Flat `FocalAdhesionDynamicsDiagnostics` base + optional
  protrusion fields** (`Optional[int]` for `linked_missing` /
  `reciprocal_missing`, `Optional[dict]` for
  `multiplier_histogram` / `max_effective_rate_per_name`).
  Single type, optional fields default to `None` for 6.3a.
  Less typed discipline; couples 6.3a's API to 6.3b's
  extensions.

### Q3 — Adopt `frozen=True, slots=True` on
`FocalAdhesionDynamicsResult`?

Two candidate positions:

- **(a) Yes, mirror HB#4 / Phase D dataclass discipline**.
  `FocalAdhesionDynamicsResult` becomes
  `@dataclass(frozen=True, slots=True)`. Tests + consumers must
  not mutate `result.cell_force_nN_xy = ...` etc. Verified
  read-only via grep — no consumer mutates fields today, so the
  flip is safe.
- **(b) No, leave plain `@dataclass`**. Keeps 6.3a/6.3b mutable
  for legacy consumer flexibility. Sister-pattern divergence
  remains.

Recommended: (a) per the Step 6 sister-gate-mirror rule and per
no-defect-on-grep verification.

### Q4 — Backward-compat shim for dict-key access

Current consumers use `result.diagnostics["key"]` (dict-key
access). If Q1 picks (a) typed dataclass, the access pattern
changes to `result.diagnostics.key` (attribute access).

Three candidates:

- **(a) Hard break: rewrite test consumer access to attribute
  form** (cleanest; sister-pattern with HB#4 typed `ECMSampledAtFAs`
  attribute-access). Test diff is 4-6 line edits across 2 test
  files.
- **(b) `__getitem__` shim on the typed dataclass**: support both
  `result.diagnostics.key` AND `result.diagnostics["key"]` via
  `__getitem__` proxy. Backward-compat but defeats the typed-key
  purpose and adds maintenance burden (every new field needs the
  shim).
- **(c) Hybrid: support attribute access by default; deprecate
  `__getitem__` with a `DeprecationWarning` for one release**.
  Over-engineering for a research codebase with no external
  consumers.

Recommended: (a) hard break. The codebase has no external
consumers; tests are internal; rewriting is cheap and clean.

---

## 3. Recommended opening positions (NOT a lock)

| Q | Recommendation | Rationale |
|---|---|---|
| Q1 typed schema | (a) typed `FocalAdhesionDynamicsDiagnostics` dataclass | Sister-pattern with HB#4 typed dataclass precedent; Step 6 rule active |
| Q2 6.3b extension | (b) separate `ProtrusionCoupledDynamicsResult` distinct from `FocalAdhesionDynamicsResult` | Mirrors HB#3 / HB#4 / Phase D phase-isolation pattern (each gets its own typed result); avoids inheritance subtleties + optional-field laxness |
| Q3 frozen+slots | (a) adopt | Sister-pattern; no-mutation verified by grep — safe flip |
| Q4 backward-compat | (a) hard break — rewrite test access | Internal codebase, no external consumers, test diff small |

**Caveats / unresolved-by-this-brief**:

- Q2 (b) "separate result types" implies `step_focal_adhesions_static`
  signature unchanged but `step_protrusion_coupled_focal_adhesions`
  signature returns a NEW type. This is a public-API break for any
  caller currently typing the return as `FocalAdhesionDynamicsResult`
  — but verified via grep that the codebase doesn't have such
  callers outside the test files themselves. Phase E composition
  cycle would import the new type cleanly.
- The audit found `FocalAdhesionDynamicsParameters` (`focal_adhesion.py:66`)
  and `ProtrusionStateMultipliers` (`protrusion_coupled_focal_adhesion.py:82`)
  use `frozen=True` without `slots=True`. **Out of scope for this
  brief** — that is Finding B2, "no action" per Codex `id=1432`.
- HB#4 `ECMToFABiasResult.diagnostics_dict` grandfathered surface
  is **explicitly out of scope** per Codex `id=1432` nuance.

---

## 4. What is *out of scope* for this brief

- HB#4 `ECMToFABiasResult.diagnostics_dict` migration (grandfathered;
  Codex `id=1432`).
- A1 `_BOUNDARY_TOL_RELATIVE` `Final[float]` cosmetic fix
  (next-touch list per Codex `id=1432`).
- B2 `frozen-no-slots` on `FocalAdhesionDynamicsParameters` /
  `ProtrusionStateMultipliers` (no action per Codex `id=1432`).
- Any change to `step_focal_adhesions_static` /
  `step_protrusion_coupled_focal_adhesions` *behavior* — this
  brief is shape-only.
- Any change to `FocalAdhesionDynamicsError` or its `failure_kind`
  set.
- Phase D / Phase E composition.

---

## 5. Sanity Gate items the future B1 lock+code will need

Per CLAUDE.md "Sanity Gate Protocol" applied to a typed-shape
refactor:

- **§1 Dimensional**: no new physics; units inherited from
  existing `cell_force_nN_xy` / `substrate_reaction_nN_xy` etc.
- **§2 Boundary**: empty FA list (already covered by 6.3a/6.3b
  tests); diagnostics fields populate to 0 / empty as
  appropriate.
- **§3 Conservation**: 6.3a sum invariants (`cell_force_nN_xy`
  pairs with `substrate_reaction_nN_xy` etc.) must continue to
  hold; the typed-shape refactor cannot alter conservation
  semantics.
- **§4 Numerical**: float64 / int types per field; no new
  tolerance.
- **§5 Sign**: no new sign convention; preserve existing.
- **§6 Measurement-protocol**: typed dataclass means each field
  has documented unit / shape. No scalarization. Runtime
  meta-test:
  `test_focal_adhesion_dynamics_result_does_not_admit_loose_dict_keys`
  to enforce the typed boundary.

**Magic-Number Block**: zero new tunables (typing-only refactor).

---

## 6. What this brief is *not*

- Not a lock. The typed shape is decided by the design-discussion
  adversarial round.
- Not a Sanity Gate doc. The Sanity Gate is written by the future
  B1 code commit, after the lock artifact lands.
- Not a Phase E activation. B1 is a 6.3a/6.3b dataclass-shape
  refactor; nothing in this brief authorizes Phase E behavior.
- Not authoritative for biological parameters (none introduced).

---

## 7. Process expectation (mirrors HB#3 / HB#4 / Phase D
precedents)

1. **Cross-room dispatch** to design-discussion with this brief
   + the 4 design questions + recommended opening positions.
   Status: `decision-needed`.
2. **Adversarial round** in design-discussion (Codex/Claude per
   memory `feedback_aggressive_design_debate.md`). HB#3 took 4
   rounds; HB#4 took 4 rounds; Phase D took 2 rounds. B1 is a
   typing-only refactor with smaller surface — 1-2 rounds
   plausible.
3. **Lock artifact** delivered as
   `docs/v2_focal_adhesion_dynamics_result_typed_locked.md`,
   parallel to HB#3 / HB#4 / Phase D lock artifacts.
4. **Cross-room dispatch back to implementation-work** with the
   locked typed shape.
5. **Sanity Gate doc** drafted by impl-work, Codex review.
6. **Code commit** — modifies
   `acs/v2/dynamics/focal_adhesion.py` (and possibly
   `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py` if Q2
   picks (b)) + test consumer rewrites + new typed dataclass
   exports through `__init__.py` (Step 6 sister-gate-mirror —
   API surface layer).
7. **Codex review** of code commit.
8. **B1 unit complete** = the audit Finding B1 closed.

implementation-work stays idle/review-capable while this design
round runs. Per `id=1432` routing, B1 unit precedes any Phase E
HB#1/#2/#5 design-discussion entries.

---

## 8. References

- Audit note (this unit's parent):
  `docs/v2_sister_gate_mirror_audit_2026_05_05.md` (commit
  `d19670d`)
- Audit MCP thread: `id=1430` (scope), `id=1431` (findings),
  `id=1432` (routing)
- Phase D lock artifacts (sister-pattern source):
  - `docs/v2_phase_d_no_op_scaffolding_locked.md` (commit
    `18b5430` + `7786b20`)
  - `acs/v2/dynamics/closed_loop_phase_d.py` (commit `60f726f` +
    `b1b2b92`)
- HB#4 lock artifacts (typed dataclass + grandfathered
  diagnostics_dict precedent):
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `acs/v2/dynamics/ecm_to_fa_bias.py`
- HB#3 lock artifacts (typed scatter primitive precedent):
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `acs/v2/dynamics/fa_to_ecm_scattering.py`
- 6.3a/6.3b modules under refactor:
  - `acs/v2/dynamics/focal_adhesion.py`
  - `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`
- 6.3a/6.3b test consumers:
  - `tests/test_v2_focal_adhesion_dynamics.py`
  - `tests/test_v2_protrusion_coupled_focal_adhesion.py`
- Memory rules:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, promoted per Codex `id=1428`)
  - `feedback_aggressive_design_debate.md`
  - `rule10_unit_derivation_in_docs.md`
  - `hard_rule_11_wording_boundary_meta_test.md`

---

## 9. Cross-room dispatch instruction (impl-work → design-discussion)

The MCP message that accompanies this brief:
- `to`: `codex` (design-discussion-pane Codex receives + Claude
  pane reads)
- `room`: `design-discussion`
- `topic`: `v2-6-3a-6-3b-focal-adhesion-dynamics-result-typed`
- `status`: `decision-needed`
- `body`: brief 4-question summary + opening recommendations +
  reference to this file
- `refs`: implementation-work `id=1432`, audit note `d19670d`,
  HB#4 / Phase D precedents

design-discussion round produces a lock artifact that supersedes
this brief; this file remains as historical opening-position
context.
