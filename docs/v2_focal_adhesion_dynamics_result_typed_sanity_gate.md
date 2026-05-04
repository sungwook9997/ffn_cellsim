# V2 6.3a / 6.3b FocalAdhesionDynamicsResult Typed Schema — Sanity Gate

**Date**: 2026-05-05 KST
**Status**: pre-execution Sanity Gate for the B1 typed-schema
migration (commit `ed5c0ca` lock + `7924730` citation fix).
Required by CLAUDE.md "Sanity Gate Protocol" before the first
execution of any new physics / numerics module.
**Source lock**: `docs/v2_focal_adhesion_dynamics_result_typed_locked.md`
(commit `ed5c0ca` initial + `7924730` Phase-D-round-1 citation
correction).
**Source brief**: `docs/v2_focal_adhesion_dynamics_result_typed_brief.md`
(commit `64db861` corrected key set; supersedes `06c3a03` initial).
**Parent audit**: `docs/v2_sister_gate_mirror_audit_2026_05_05.md`
(commit `d19670d`).
**Target modules**: `acs/v2/dynamics/focal_adhesion.py` +
`acs/v2/dynamics/protrusion_coupled_focal_adhesion.py` +
`acs/v2/dynamics/__init__.py` + `acs/v2/__init__.py` +
`tests/test_v2_focal_adhesion_dynamics.py` +
`tests/test_v2_protrusion_coupled_focal_adhesion.py`.

This Sanity Gate is the impl-work side's gate before code lands.
B1 is a **typing migration**; the burden is to prove **no
behavior change**, NOT to re-validate the underlying 6.3a / 6.3b
physics (already gated under their respective locks).

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

Per locked §1:

- New typed dataclass `FocalAdhesionDynamicsDiagnostics`
  (4 fields: `n_adhesions`, `aggregate_cell_force_nN_xy`,
  `aggregate_substrate_reaction_nN_xy`,
  `max_traction_magnitude_nN`) with `frozen=True, slots=True`.
- `FocalAdhesionDynamicsResult` flips to
  `@dataclass(frozen=True, slots=True)`; `diagnostics` field
  type `dict` → `FocalAdhesionDynamicsDiagnostics`.
- New typed dataclass
  `ProtrusionCoupledDynamicsDiagnostics(FocalAdhesionDynamicsDiagnostics)`
  (4 base fields + 4 6.3b-specific: `linked_missing`,
  `reciprocal_missing`, `multiplier_histogram: Mapping[str, Counter]`,
  `max_effective_rate_per_name: Mapping[str, float]`) with
  `frozen=True, slots=True`.
- New typed result subclass
  `ProtrusionCoupledDynamicsResult(FocalAdhesionDynamicsResult)`
  with `frozen=True, slots=True`; `diagnostics` field type
  override to `ProtrusionCoupledDynamicsDiagnostics` (plain
  override, no `# type: ignore[assignment]` per Y4).
- Producer updates:
  - `step_focal_adhesions_static` (`focal_adhesion.py:381-399`)
    constructs `FocalAdhesionDynamicsDiagnostics(...)` instead of
    a dict literal.
  - `step_protrusion_coupled_focal_adhesions`
    (`protrusion_coupled_focal_adhesion.py:422-444`) constructs
    `ProtrusionCoupledDynamicsDiagnostics(...)` and returns
    `ProtrusionCoupledDynamicsResult(...)`. Return annotation
    `FocalAdhesionDynamicsResult` →
    `ProtrusionCoupledDynamicsResult`.
- Exports: 3 new names through both
  `acs/v2/dynamics/__init__.py` and `acs/v2/__init__.py` (with
  `__all__` inclusion at both surfaces).
- Test consumer migration:
  `result.diagnostics["key"]` →
  `result.diagnostics.key` across both test files. New B1
  meta-tests + isinstance + frozen + return-type tests per
  locked §4.

### Explicitly out of scope (will FAIL gate if introduced)

Per locked §1 Forbidden:

- 6.3a / 6.3b force / rate algebra change.
- HB#4 `diagnostics_dict` surface change (grandfathered;
  separate unit).
- Nested mapping deep freeze
  (`multiplier_histogram` /
  `max_effective_rate_per_name` stay
  `Mapping[...]` — separate unit).
- `__getitem__` shim, deprecation hybrid, or any
  backward-compat dict access on the result.
- Preemptive `# type: ignore[assignment]` on the diagnostics
  override (Y4 deferred to a future strict-type-check unit).
- `Generic[T_Diag]` parameterization (overengineered for B1).
- Any other audit finding (A1 cosmetic next-touch list, B2
  frozen-no-slots on Parameters/Multipliers).
- Any change to `FocalAdhesionDynamicsError` or its
  `failure_kind` set.
- Phase D / Phase E composition.

If any of these appear, the gate FAILs and the unit halts to
surface the scope creep to PI per CLAUDE.md "Sanity Gate Failure
handling".

---

## 1. Dimensional analysis

The typed migration introduces **no new unit chain**. Every
diagnostics field carries the same unit as the corresponding
dict value carried before B1 (locked §3 item 1):

- `n_adhesions: int` — dimensionless count.
- `aggregate_cell_force_nN_xy: tuple[float, float]` — `nN × nN`
  (2-component force-per-cell sum).
- `aggregate_substrate_reaction_nN_xy: tuple[float, float]` —
  `nN × nN` (Newton-3 partner of cell force, per 6.3a contract).
- `max_traction_magnitude_nN: float` — `nN` (max of L2 norm of
  per-FA cell force).
- `linked_missing: int` — count.
- `reciprocal_missing: int` — count.
- `multiplier_histogram: Mapping[str, Counter[float]]` —
  dimensionless multiplier × count per rate name.
- `max_effective_rate_per_name: Mapping[str, float]` — `1/s` per
  rate name (locked §1 docstring source `1/s`).

The migration also preserves the 5 physical fields on
`FocalAdhesionDynamicsResult` unchanged:
- `updated_adhesions: tuple[FocalAdhesionState, ...]`
- `cell_force_nN_xy: np.ndarray` (`nN`, shape `(N_FA, 2)` float64)
- `substrate_reaction_nN_xy: np.ndarray` (`nN`, shape `(N_FA, 2)` float64)
- `radial_components: np.ndarray` (`nN`, shape `(N_FA,)` float64)
- `tangential_components: np.ndarray` (`nN`, shape `(N_FA,)` float64)

### Status

**PASS**. No new dimensional chain; every field 1:1 with the
prior dict-value unit. No CFL / stability bound applies (no time
integration at this layer; producer-side aggregation only).

---

## 2. Boundary cases

### 2.1 Empty FA list (`adhesions = ()`)

Per locked §3 item 2: `n_adhesions = 0` produces
`aggregate_cell_force_nN_xy = (0.0, 0.0)`,
`aggregate_substrate_reaction_nN_xy = (0.0, 0.0)`,
`max_traction_magnitude_nN = 0.0`. The producer at
`focal_adhesion.py:388` already handles this:
`np.linalg.norm(cell_force, axis=1).max() if n else 0.0` — the
typed dataclass migration preserves this via the
`FocalAdhesionDynamicsDiagnostics(...)` constructor call.

For 6.3b: `linked_missing = 0`, `reciprocal_missing = 0`,
`multiplier_histogram` is an empty dict (or a dict of empty
Counters per the producer's
`{name: Counter(round(v, 12) for v in values) for ...}` with
empty `values`), `max_effective_rate_per_name` returns the
producer's empty-dict default per `_max_effective_rate_per_name`.

Test 1 (`test_focal_adhesion_diagnostics_is_typed_dataclass`)
covers the typed-dataclass return on empty + non-empty.

### 2.2 Single-FA / multi-FA interior

Producer aggregates the same way; the typed dataclass migration
is pass-through. The 5 existing physical fields and the 4 (6.3a)
or 8 (6.3b) diagnostic fields all populate per the prior dict
producer.

### 2.3 Frozen-instance assignment

Test 3 (`test_focal_adhesion_result_frozen`) asserts
`with pytest.raises(dataclasses.FrozenInstanceError):
result.cell_force_nN_xy = something` — the new `frozen=True`
flag must reject post-construction assignment.

### Status

**PASS**. All three boundary cases are explicitly locked with
runtime invariants and/or tests. The producer's
`if n else 0.0` empty branch is unchanged — typed migration
preserves the value, only its container changes.

---

## 3. Conservation invariants

The migration preserves all 6.3a / 6.3b conservation laws by
construction (zero arithmetic change; dataclass packaging only):

### 3.1 Pure-function constraint

The locked §1 forbidden list explicitly forbids force/rate
algebra change. The producer's algebra at
`focal_adhesion.py:381-399` and
`protrusion_coupled_focal_adhesion.py:422-444` is preserved
verbatim; only the final `return FocalAdhesionDynamicsResult(...)`
changes from a dict-arg to a typed-dataclass-arg constructor
call. **No mutation** is introduced — the typed dataclasses are
all `frozen=True, slots=True`.

Test 7 (`test_b1_no_behavior_change_baseline`) is the
regression test: pre-B1 dict values for a fixed FA configuration
and post-B1 attribute values must be numerically identical.

### 3.2 6.3a Newton-3 contract (preserved)

`cell_force_nN_xy[i] = -substrate_reaction_nN_xy[i]` per FA, and
their component-sums populate
`aggregate_cell_force_nN_xy` /
`aggregate_substrate_reaction_nN_xy`. The migration packages
these into the typed diagnostics field unchanged. Existing 6.3a
tests verifying Newton-3 still pass after the
`result.diagnostics.aggregate_*` attribute-access rewrite.

### 3.3 6.3b delegate-to-6.3a structural invariant

6.3b iterates per-FA and delegates to 6.3a via
`step_focal_adhesions_static` for each single-FA call (per
producer pattern). The typed migration preserves this — 6.3b's
own typed `ProtrusionCoupledDynamicsDiagnostics` is built from
the 6.3a-shared field values plus the 4 6.3b extension values.

The IS-A relation (`ProtrusionCoupledDynamicsResult is a
FocalAdhesionDynamicsResult`) is locked and tested by test 4
(`test_protrusion_coupled_result_is_subclass`).

### 3.4 Determinism

The typed migration is pass-through; producer determinism
preserved. No RNG, no global state, no order-dependent
operations introduced.

### Status

**PASS**. Wrapper introduces zero new conservation laws. All
existing 6.3a / 6.3b invariants preserved by construction
(producer algebra unchanged). Test 7 (regression) is the central
no-behavior-change conservation test.

---

## 4. Numerical sanity

### 4.1 Float precision

`np.float64` throughout for the 5 physical fields (unchanged).
`tuple[float, float]` for aggregated 2-component fields (unchanged
— Python `float` is double-precision, matches `np.float64`
mathematically). `int` for counts (unchanged).
`Counter[float]` keys are float values rounded to 12 decimals
per producer
`Counter(round(v, 12) for v in values)` (unchanged).

### 4.2 No new tolerance constants

The migration introduces no new numerical tolerance. The
existing 12-decimal rounding tolerance for
`multiplier_histogram` is unchanged. **Magic-Number Block: zero
new tunables** — see §7.

### 4.3 Per-call work

Per-call work = unchanged producer work + 1 typed-dataclass
constructor call (O(1) field assignment for `slots=True` types).
No new asymptotic complexity.

### 4.4 Stability

No time integration at this layer. CFL / explicit-stepping
bounds do not apply. Sister gates (6.3a `step_focal_adhesions_static`,
6.3b `step_protrusion_coupled_focal_adhesions`) own those checks
at their respective primitive layers.

### Status

**PASS**. Float precision unchanged; no new tolerance; per-call
work unchanged; no time integration at this layer.

---

## 5. Sign / sense check

The migration performs no arithmetic. Sign / sense is inherited
1:1 from the producer:

- `cell_force_nN_xy` / `substrate_reaction_nN_xy`: 6.3a
  Newton-3 sign relation preserved (locked §3 item 5).
- `aggregate_*_nN_xy`: pass-through `tuple(arr.sum(axis=0).tolist())`
  preserves sign per axis.
- `max_traction_magnitude_nN`: pass-through
  `np.linalg.norm(...).max()` is non-negative.
- `linked_missing` / `reciprocal_missing`: counts ≥ 0 per
  producer.
- `multiplier_histogram` keys: rounded multiplier values; sign
  inherited from producer (multipliers may be 0 or positive per
  6.3b lock; never negative).
- `max_effective_rate_per_name`: rates ≥ 0 per
  `_max_effective_rate_per_name` (`Absent rate (None) is treated
  as 0.0`).

### Status

**PASS**. No migration-introduced sign. Producer sign-preservation
inherited.

---

## 6. Measurement-protocol consistency (Hard Rule 11)

The most critical Sanity Gate item for a typed-migration unit,
since the risk is silent value-renaming or implicit-coercion at
the pre-B1 / post-B1 boundary.

### 6.1 1:1 dict-key → attribute mapping

Every renamed/typed field maps 1:1 to the prior dict key:

| Pre-B1 dict key | Post-B1 attribute | Type |
|---|---|---|
| `result.diagnostics["n_adhesions"]` | `result.diagnostics.n_adhesions` | `int` |
| `result.diagnostics["aggregate_cell_force_nN_xy"]` | `result.diagnostics.aggregate_cell_force_nN_xy` | `tuple[float, float]` |
| `result.diagnostics["aggregate_substrate_reaction_nN_xy"]` | `result.diagnostics.aggregate_substrate_reaction_nN_xy` | `tuple[float, float]` |
| `result.diagnostics["max_traction_magnitude_nN"]` | `result.diagnostics.max_traction_magnitude_nN` | `float` |
| `result.diagnostics["linked_missing"]` (6.3b) | `result.diagnostics.linked_missing` | `int` |
| `result.diagnostics["reciprocal_missing"]` (6.3b) | `result.diagnostics.reciprocal_missing` | `int` |
| `result.diagnostics["multiplier_histogram"]` (6.3b) | `result.diagnostics.multiplier_histogram` | `Mapping[str, Counter]` |
| `result.diagnostics["max_effective_rate_per_name"]` (6.3b) | `result.diagnostics.max_effective_rate_per_name` | `Mapping[str, float]` |

No keys dropped, no keys added, no implicit type coercion. Test
7 (`test_b1_no_behavior_change_baseline`) enforces numerical
equality of every value across the migration.

### 6.2 No `__getitem__` shim

Test 2 (`test_focal_adhesion_diagnostics_no_dict_access`) asserts
`with pytest.raises(TypeError): result.diagnostics["n_adhesions"]`.
This is the runtime measurement-protocol guard against a future
refactor that adds a dict-access shim — the locked §1 Forbidden
list explicitly rejects this.

### 6.3 IS-A subclass relation tested

Test 4 (`test_protrusion_coupled_result_is_subclass`) asserts:
- `isinstance(result, FocalAdhesionDynamicsResult)` — Phase E
  composition consumers needing only the base force arrays
  accept the subclass.
- `type(result) is ProtrusionCoupledDynamicsResult` — 6.3b
  returns the subclass, not the base, ensuring 6.3b-specific
  diagnostics are accessible.

### 6.4 Frozen contract tested

Test 3 (`test_focal_adhesion_result_frozen`) asserts a
`FrozenInstanceError` on field assignment. This is the runtime
measurement-protocol guard against future code that introduces
post-construction mutation.

### 6.5 Sister-gate-mirror application (Step 6 of 6-step batch)

This Sanity Gate doc itself is being written under the
freshly-promoted Step 6 sister-gate-mirror rule. The sister
patterns this gate mirrors:

- **Code (validation order)**: B1 introduces no FA position
  validation (delegates fully to 6.3a / 6.3b primitives, which
  themselves pre-date the HB#3+#4 structural-len-check fix).
  Vacuous Step 6 layer for code-level validation order.
- **Design (typed dataclass)**: HB#4 `ECMSampledAtFAs` /
  `ECMToFABiasResult` and Phase D `FAToECMResponseResult` /
  `PhaseDNoOpStepResult` precedent uses `frozen=True, slots=True`.
  B1 mirrors this for the new types ✓.
- **API surface**: 3 new exports through both
  `acs/v2/__init__.py` and `acs/v2/dynamics/__init__.py` per
  locked §1 (Phase D `b1b2b92` precedent for the export
  pattern).
- **Failure-kind discipline**: B1 introduces no new error class
  or `failure_kind` (typed migration only). Existing
  `FocalAdhesionDynamicsError` continues unchanged.

### Status

**PASS**. Five layers of measurement-protocol guard:
- Locked §1 forbidden list (text-level).
- §6.1 1:1 dict-key → attribute mapping (table above).
- §6.2 runtime meta-test (`__getitem__` raises TypeError).
- §6.3 IS-A isinstance + type-equality test.
- §6.4 frozen-instance assignment test.
- §6.5 Step 6 sister-gate-mirror application across all four
  layers (code / design / API surface / failure-kind).

---

## 7. Magic-Number Block check

**Zero new tunables introduced at any layer.**

The migration preserves:
- Producer algebra (no new constants).
- 12-decimal rounding for `multiplier_histogram` (6.3b lock,
  unchanged).
- Existing `_DT_RATE_SAFETY_MARGIN` and other 6.3a / 6.3b
  numerics (out of scope per locked §1).

Magic-Number Block 3-test verdict:
1. **Derivable**: PASS by inheritance — B1 has no constants to
   derive.
2. **Grid-invariant**: PASS by inheritance — typed migration
   doesn't depend on grid resolution.
3. **Fitting**: PASS — no numbers chosen to fit any target.

### Status

**PASS**. Zero migration-introduced magic numbers.

---

## 8. Test catalog (B1-specific additions per locked §4 + dict-attribute rewrites)

Owned by `tests/test_v2_focal_adhesion_dynamics.py` and
`tests/test_v2_protrusion_coupled_focal_adhesion.py` (post-B1
modifications). Each new test maps to a locked invariant in
`docs/v2_focal_adhesion_dynamics_result_typed_locked.md` §4.

### Existing test rewrites (4-6 line edits in each file)

- `test_v2_focal_adhesion_dynamics.py:323-324` (and any other
  dict-key access): `result.diagnostics["aggregate_cell_force_nN_xy"]`
  → `result.diagnostics.aggregate_cell_force_nN_xy`,
  `result.diagnostics["aggregate_substrate_reaction_nN_xy"]` →
  `result.diagnostics.aggregate_substrate_reaction_nN_xy`.
- `test_v2_protrusion_coupled_focal_adhesion.py:287-288` (and
  any other dict-key access):
  `result.diagnostics["reciprocal_missing"]` →
  `result.diagnostics.reciprocal_missing`,
  `result.diagnostics["linked_missing"]` →
  `result.diagnostics.linked_missing`.

### New B1 tests (7 new, per locked §4)

1. `test_focal_adhesion_diagnostics_is_typed_dataclass` —
   `isinstance(result.diagnostics, FocalAdhesionDynamicsDiagnostics)`
   in 6.3a; `isinstance(result.diagnostics,
   ProtrusionCoupledDynamicsDiagnostics)` in 6.3b.
2. `test_focal_adhesion_diagnostics_no_dict_access` —
   `with pytest.raises(TypeError):
   result.diagnostics["n_adhesions"]`.
3. `test_focal_adhesion_result_frozen` —
   `with pytest.raises(dataclasses.FrozenInstanceError):
   result.cell_force_nN_xy = something`.
4. `test_protrusion_coupled_result_is_subclass` —
   `isinstance(protrusion_coupled_result,
   FocalAdhesionDynamicsResult)` AND
   `type(protrusion_coupled_result) is
   ProtrusionCoupledDynamicsResult`.
5. `test_protrusion_coupled_diagnostics_is_subclass` —
   `isinstance(diag, FocalAdhesionDynamicsDiagnostics)` AND
   `type(diag) is ProtrusionCoupledDynamicsDiagnostics`.
6. `test_protrusion_coupled_diagnostics_n_adhesions_inherited`
   — assert the 6.3b diagnostics still surfaces `n_adhesions`
   correctly via the inherited attribute.
7. `test_b1_no_behavior_change_baseline` (regression) — pre-B1
   dict values for a fixed FA configuration and post-B1
   attribute values must be numerically identical (the central
   no-behavior-change gate).

If a regression case surfaces during code commit (e.g., Codex
review catches a subclass-override edge case), an 8th test is
added with explicit lock reference; the count is not capped at
7.

---

## 9. Gate verdict

| Item | Status | Notes |
|---|---|---|
| §1 Dimensional | **PASS** | Every typed field 1:1 with prior dict-value unit |
| §2 Boundary | **PASS** | Empty-FA / single-FA / frozen-instance assignment all locked |
| §3 Conservation | **PASS** | Producer algebra unchanged; 6.3a Newton-3 + 6.3b delegate-to-6.3a preserved |
| §4 Numerical | **PASS** | float64 unchanged; no new tolerance; O(1) typed-dataclass overhead |
| §5 Sign | **PASS** | No arithmetic; producer sign-preservation inherited |
| §6 Measurement-protocol | **PASS** | Five-layer guard (locked text + 1:1 mapping + meta-test + IS-A test + frozen test + Step 6 sister-gate-mirror) |
| §7 Magic-Number Block | **PASS** | Zero new tunables |

**Overall**: gate PASS. impl-work is clear to commit the typing
migration after Codex review of this Sanity Gate doc, mirroring
HB#3 / HB#4 / Phase D Sanity Gate review precedent.

### Outstanding before code lands

- Codex review of this Sanity Gate doc.
- Code commit must reproduce the locked §1 forbidden list at the
  module docstring level (text-level guard layered on top of
  the runtime tests).
- All 7 new B1 tests + dict→attribute rewrites must pass at
  first commit; no "TODO test_X" placeholders.
- `pytest tests/test_v2_focal_adhesion_dynamics.py
  tests/test_v2_protrusion_coupled_focal_adhesion.py
  tests/test_v2_fa_to_ecm_scattering.py
  tests/test_v2_ecm_to_fa_bias.py
  tests/test_v2_closed_loop_phase_d.py` must PASS at commit time
  (no sister-gate regression).

---

## 10. References

- B1 lock: `docs/v2_focal_adhesion_dynamics_result_typed_locked.md`
  (commit `ed5c0ca` initial + `7924730` Phase-D-round-1 citation
  fix).
- B1 brief (superseded by lock):
  `docs/v2_focal_adhesion_dynamics_result_typed_brief.md`
  (commit `06c3a03` + `64db861` correction).
- Parent audit:
  `docs/v2_sister_gate_mirror_audit_2026_05_05.md` (commit
  `d19670d`).
- Sister-gate Sanity Gate precedents:
  - `docs/v2_fa_to_ecm_scattering_sanity_gate.md` (HB#3)
  - `docs/v2_ecm_to_fa_bias_sanity_gate.md` (HB#4)
  - `docs/v2_phase_d_no_op_scaffolding_sanity_gate.md` (Phase D)
- 6.3a / 6.3b modules under refactor:
  - `acs/v2/dynamics/focal_adhesion.py`
  - `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`
- Test consumers under migration:
  - `tests/test_v2_focal_adhesion_dynamics.py`
  - `tests/test_v2_protrusion_coupled_focal_adhesion.py`
- Memory rules informing this gate:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, per Codex `id=1428`)
  - `feedback_aggressive_design_debate.md`
  - `cadence_promise_must_send_even_when_idle.md`
- Design-discussion thread: `id=1434` → `1441` (4-round lock + ack).
- impl-work review thread: `id=1447` (citation fix request) →
  `7924730` (citation fix commit) → `id=1449` (Sanity Gate
  entry ack).
