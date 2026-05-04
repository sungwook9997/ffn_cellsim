# V2 6.3a / 6.3b FocalAdhesionDynamicsResult Typed Schema — Locked Design

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (4-round adversarial lock,
PI id=809 aggressive debate posture applied, PI id=1008/1057/1161
autonomy + visible deliverable focus)
**Source unit**: design-discussion `topic=v2-6-3a-6-3b-focal-adhesion-dynamics-result-typed`,
MCP id 1434–1441 (rounds 1–4 + ack)
**Parent audit**: `docs/v2_sister_gate_mirror_audit_2026_05_05.md` (commit `d19670d`),
Finding B1 routed to separate cycle by Codex `id=1432`
**Brief**: `docs/v2_focal_adhesion_dynamics_result_typed_brief.md`
(opening commit `06c3a03`, Q1/Q2 key-set/type completeness correction
`64db861` — does not contradict this lock; Y2 already accepts
`n_adhesions` in base diagnostics, §1 already pins the 6.3b extension
diagnostics shape)
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for the B1 Sanity Gate doc + code entry.

---

## 0. Scope

B1 is a **typed-schema migration** of the `FocalAdhesionDynamicsResult`
return type for both 6.3a (`step_focal_adhesions_static`) and 6.3b
(`step_protrusion_coupled_focal_adhesions`). It is a **shape-only
refactor** — NO behavior change to rate algebra, traction algebra,
state transitions, or any 6.3a / 6.3b numerical contract.

**Why now**: sister-pattern divergence audit (Finding B1) caught two
defects relative to HB#3 / HB#4 / Phase D locked patterns:
1. `@dataclass` plain (NOT `frozen=True, slots=True`) — diverges from
   HB#3 / HB#4 / Phase D result types
2. `diagnostics: dict = field(default_factory=dict)` — loose untyped
   dict, the exact regression caught at Phase D round 1 C1 (`id=1395`)

HB#4 `diagnostics_dict` keeps its grandfathered surface (separate unit,
out of B1 scope). 6.3a / 6.3b is older code that needs typed-schema
discipline + immutability before Phase E composition.

---

## 1. Final Lock Summary

### Module changes

#### `acs/v2/dynamics/focal_adhesion.py`

**New typed diagnostics dataclass** (4 fields):

```python
@dataclass(frozen=True, slots=True)
class FocalAdhesionDynamicsDiagnostics:
    """Typed diagnostics for one 6.3a static FA dynamics step."""
    n_adhesions: int
    aggregate_cell_force_nN_xy: tuple[float, float]
    aggregate_substrate_reaction_nN_xy: tuple[float, float]
    max_traction_magnitude_nN: float
```

**`FocalAdhesionDynamicsResult` change** — `@dataclass` → `@dataclass(frozen=True, slots=True)`,
diagnostics field type `dict` → `FocalAdhesionDynamicsDiagnostics`:

```python
@dataclass(frozen=True, slots=True)
class FocalAdhesionDynamicsResult:
    """Outputs of one 6.3a static FA dynamics step."""
    updated_adhesions: tuple[FocalAdhesionState, ...]
    cell_force_nN_xy: np.ndarray
    substrate_reaction_nN_xy: np.ndarray
    radial_components: np.ndarray
    tangential_components: np.ndarray
    diagnostics: FocalAdhesionDynamicsDiagnostics
```

**6.3a producer** (`step_focal_adhesions_static`, lines 381–399): replace
the dict literal with a `FocalAdhesionDynamicsDiagnostics(...)` constructor
call (same field values, same units).

#### `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`

**New typed diagnostics subclass** (base 4 fields + 4 6.3b-specific fields):

```python
@dataclass(frozen=True, slots=True)
class ProtrusionCoupledDynamicsDiagnostics(FocalAdhesionDynamicsDiagnostics):
    """Typed diagnostics for one 6.3b protrusion-coupled FA step."""
    linked_missing: int
    reciprocal_missing: int
    multiplier_histogram: Mapping[str, Counter]
    max_effective_rate_per_name: Mapping[str, float]
```

**New typed result subclass** — plain subclass override on the
`diagnostics` field (no `# type: ignore`):

```python
@dataclass(frozen=True, slots=True)
class ProtrusionCoupledDynamicsResult(FocalAdhesionDynamicsResult):
    """Outputs of one 6.3b protrusion-coupled FA step. IS-A 6.3a result."""
    diagnostics: ProtrusionCoupledDynamicsDiagnostics
```

**`step_protrusion_coupled_focal_adhesions` change**: return annotation
`FocalAdhesionDynamicsResult` → `ProtrusionCoupledDynamicsResult`. Producer
(lines 422–444) replaces the dict literal with a
`ProtrusionCoupledDynamicsDiagnostics(...)` constructor call and the
result constructor with `ProtrusionCoupledDynamicsResult(...)`.

#### Exports

`acs/v2/dynamics/__init__.py` adds 3 new names to imports + `__all__`:
- `FocalAdhesionDynamicsDiagnostics`
- `ProtrusionCoupledDynamicsDiagnostics`
- `ProtrusionCoupledDynamicsResult`

`acs/v2/__init__.py` mirrors the 3 additions.

#### Test consumer migration

- `tests/test_v2_focal_adhesion_dynamics.py` — rewrite
  `result.diagnostics["key"]` → `result.diagnostics.key` (4–6 line edits
  near lines 323–324 + any other dict access).
- `tests/test_v2_protrusion_coupled_focal_adhesion.py` — same pattern
  near lines 287–288 + any other dict access.
- **New meta-test (each file)**: assert
  `with pytest.raises(TypeError): result.diagnostics["n_adhesions"]` —
  enforces no `__getitem__` shim.
- **New isinstance test (6.3b only)**:
  `assert isinstance(result, FocalAdhesionDynamicsResult)` — guarantees
  IS-A relation for Phase E composition.
- **New return-type test (6.3b only)**:
  `assert type(result) is ProtrusionCoupledDynamicsResult` — guarantees
  6.3b returns the subclass, not the base.

### Forbidden in B1

- 6.3a / 6.3b force / rate algebra change
- HB#4 `diagnostics_dict` surface change (grandfathered, separate unit)
- Nested mapping deep freeze (`multiplier_histogram` /
  `max_effective_rate_per_name` stay `Mapping[...]` — separate unit)
- `__getitem__` shim, deprecation hybrid, or any backward-compat dict
  access on the result
- Preemptive `# type: ignore[assignment]` on the diagnostics override
  (handle in a future tool-specific unit if a strict type checker is
  introduced)
- Generic[T_Diag] parameterization (overengineered for B1)
- Any other audit finding (A1 cosmetic next-touch list, B2 frozen-no-slots
  on Parameters/Multipliers — separate units)

---

## 2. Reasoned-acceptance trace (Y1–Y5)

The lock converged after 4 rounds. Each Y is a Claude concession with
the round in which it was accepted, preserving the adversarial-debate
audit trail (PI id=809 posture).

**Y1 (round 2, accepted Codex Q2 = subclass approach)**:
- Claude opening lean: separate `ProtrusionCoupledDynamicsResult`
  duplicating the 5 physical fields.
- Codex challenge: subclass `ProtrusionCoupledDynamicsResult(FocalAdhesionDynamicsResult)`
  + `ProtrusionCoupledDynamicsDiagnostics(FocalAdhesionDynamicsDiagnostics)`
  preserves shared physics + isolates diagnostics, prevents drift in the
  5 shared physical fields, and naturally supports Phase E composition
  (consumers needing only base force arrays accept the base type).
- Resolution: subclass for both result and diagnostics.

**Y2 (round 2, accepted Codex Challenge 2 = `n_adhesions` in base)**:
- Brief Q1 listing omitted `n_adhesions` from the proposed base
  diagnostics fields.
- Codex catch: producers in both 6.3a and 6.3b currently emit
  `n_adhesions` — dropping it would be a behavior regression masquerading
  as a typing migration.
- Resolution: `n_adhesions: int` is field 1 of the base diagnostics
  dataclass.

**Y3 (round 2, accepted Codex Challenge 3 = loose typing for nested mappings)**:
- Open question whether to also freeze the inner `Counter` /
  `dict[str, float]` structures.
- Codex argument: B1 is shape-only. Deep-freezing nested mappings is a
  separate concern that may need a wrapper type or normalization step
  — out of scope here.
- Resolution: `multiplier_histogram: Mapping[str, Counter]`,
  `max_effective_rate_per_name: Mapping[str, float]`. Future unit can
  freeze if it becomes a real issue.

**Y4 (round 4, accepted Codex covariant-override = (B) plain override)**:
- Claude lean (round 2): option (A) `# type: ignore[assignment]` on the
  subclass diagnostics override, single line.
- Codex pushback (round 3): option (B) plain override, no preemptive
  ignore. Repo doesn't currently run strict mypy/pyright; preemptive
  ignore is noise and can hide real issues. Defer the ignore to whatever
  future unit introduces a strict type-check job.
- Resolution: plain subclass override. Tests assert runtime contract
  via `isinstance` + `type(result) is …`.

**Y5 (round 4, accepted Codex lock-artifact = short doc not commit-only)**:
- Claude floated (round 2): commit-message-only might suffice for a
  small typing migration.
- Codex pushback (round 3): API-shape change crossing 6.3a, 6.3b,
  exports, and tests; HB#3 / HB#4 / Phase D precedent all produced
  short lock docs; future sister-pattern audit re-runs need a doc to
  reference.
- Resolution: this document.

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

B1 is a typing migration. The Sanity Gate is correspondingly narrow —
the burden is to prove **no behavior change**, not to re-validate the
underlying physics (already gated under 6.3a / 6.3b locks).

1. **Units**: every diagnostics field carries the same unit as the
   corresponding dict value previously did:
   - `n_adhesions`: count (dimensionless int)
   - `aggregate_cell_force_nN_xy`: 2-vector, each component in nN
   - `aggregate_substrate_reaction_nN_xy`: 2-vector, each component in nN
   - `max_traction_magnitude_nN`: nN
   - `linked_missing`: count
   - `reciprocal_missing`: count
   - `multiplier_histogram`: dimensionless multiplier counts per rate name
   - `max_effective_rate_per_name`: 1/s per rate name
2. **Boundary**: `n_adhesions=0` produces
   `aggregate_cell_force_nN_xy=(0.0, 0.0)`,
   `aggregate_substrate_reaction_nN_xy=(0.0, 0.0)`,
   `max_traction_magnitude_nN=0.0` exactly as the dict producer did
   (line 388: `np.linalg.norm(cell_force, axis=1).max() if n else 0.0`).
3. **Conservation**: 6.3a Newton-3 contract unchanged (no producer
   algebra change), 6.3b delegates to 6.3a per existing per-FA loop.
4. **Numerical**: float widths unchanged (tuple of float, ndarray
   dtypes unchanged).
5. **Sign**: no force-direction sign change; diagnostics fields are
   pass-through aggregations of the same arrays.
6. **Measurement-protocol**: every renamed/typed field maps 1:1 to the
   prior dict key; existing tests (after attribute-access rewrite)
   still measure the same quantity. Add the meta-test
   (`__getitem__` raises) and the `isinstance` test (6.3b IS-A 6.3a)
   to lock the contract.

---

## 4. Test Catalog (B1-specific additions)

Beyond the dict-access → attribute-access rewrites, add the following
new tests:

- `test_focal_adhesion_diagnostics_is_typed_dataclass` — assert
  `isinstance(result.diagnostics, FocalAdhesionDynamicsDiagnostics)` in
  6.3a and `isinstance(result.diagnostics, ProtrusionCoupledDynamicsDiagnostics)`
  in 6.3b.
- `test_focal_adhesion_diagnostics_no_dict_access` — assert
  `with pytest.raises(TypeError): result.diagnostics["n_adhesions"]`.
- `test_focal_adhesion_result_frozen` — assert
  `with pytest.raises(dataclasses.FrozenInstanceError): result.cell_force_nN_xy = something`.
- `test_protrusion_coupled_result_is_subclass` — assert
  `isinstance(protrusion_coupled_result, FocalAdhesionDynamicsResult)`
  AND `type(protrusion_coupled_result) is ProtrusionCoupledDynamicsResult`.
- `test_protrusion_coupled_diagnostics_is_subclass` — assert
  `isinstance(diag, FocalAdhesionDynamicsDiagnostics)` AND
  `type(diag) is ProtrusionCoupledDynamicsDiagnostics`.
- `test_protrusion_coupled_diagnostics_n_adhesions_inherited` —
  assert the 6.3b diagnostics still surfaces `n_adhesions` correctly.
- `test_b1_no_behavior_change_baseline` (regression) — pre-B1 dict
  values for a fixed FA configuration and post-B1 attribute values must
  be numerically identical.

---

## 5. Files

- `docs/v2_focal_adhesion_dynamics_result_typed_brief.md` (existing,
  opening brief, commit `06c3a03`)
- `docs/v2_focal_adhesion_dynamics_result_typed_locked.md` (this file,
  source of truth)
- `acs/v2/dynamics/focal_adhesion.py` (modify lines 153–162 +
  381–399; add `FocalAdhesionDynamicsDiagnostics` near line 153)
- `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py` (modify
  imports + return annotation at 218 + producer at 422–444; add
  `ProtrusionCoupledDynamicsDiagnostics` and
  `ProtrusionCoupledDynamicsResult` near top)
- `acs/v2/dynamics/__init__.py` (3 new exports + alphabetical reorder
  of `__all__`)
- `acs/v2/__init__.py` (3 new mirror exports)
- `tests/test_v2_focal_adhesion_dynamics.py` (attribute migration +
  new tests)
- `tests/test_v2_protrusion_coupled_focal_adhesion.py` (attribute
  migration + new tests)

---

## 6. References

- Brief (opening position): `docs/v2_focal_adhesion_dynamics_result_typed_brief.md`
  (opening `06c3a03`, Q1/Q2 completeness correction `64db861`)
- Sister-gate-mirror audit: `docs/v2_sister_gate_mirror_audit_2026_05_05.md`
- HB#3 lock (sister precedent): `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
- HB#4 lock (sister precedent): `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
- Phase D no-op scaffolding lock: `docs/v2_phase_d_no_op_scaffolding_locked.md`
- 6.3a static FA traction preflight: design-discussion
  `topic=v2-layer-2-separated-dynamics`, MCP id 1070–1081
- 6.3b protrusion-coupled lock: `docs/v2_63b_protrusion_coupling_locked.md`
- Closed-loop ECM gate phased plan: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Forward roadmap: `docs/v2_phase1_forward_roadmap.md`
- Adversarial debate posture: memory `feedback_aggressive_design_debate.md`,
  PI id=809
- Cadence rule: memory `cadence_promise_must_send_even_when_idle.md`,
  PI id=1025/1031
- Step 6 pre-commit batch rule: memory `design_note_pre_commit_batch.md`,
  PI/Codex id=1428

---

## 7. Cross-room dispatch

This file is the design-team input to implementation-work for:
1. impl Claude writes B1 Sanity Gate doc
   (`docs/v2_focal_adhesion_dynamics_result_typed_sanity_gate.md`) from
   this lock — narrow scope per §3.
2. impl Codex review (5 focus per Codex `id=1441`):
   - no behavior changes to 6.3a/6.3b force/rate algebra
   - no dict diagnostics on 6.3a/6.3b result
   - tests migrated to attribute access AND dict access fails
   - 6.3b returns `ProtrusionCoupledDynamicsResult` AND remains
     `isinstance(..., FocalAdhesionDynamicsResult)`
   - exports include all 3 new typed classes through both
     `acs/v2/dynamics/__init__.py` and `acs/v2/__init__.py`
3. On Sanity Gate PASS: B1 code commit (typing migration, no behavior
   change, all existing tests rewrite to pass + new B1-specific tests
   per §4).
4. After commit: design-discussion idle, awaiting either impl backflow
   or Phase E hard-blocker entry (HB#1 constitutive-law direction, HB#2
   saturation form, HB#5 Lyapunov metric, effective_stiffness law
   decision).

Rounds 1–4 of the design lock (+ ack) are MCP id 1434–1441.
