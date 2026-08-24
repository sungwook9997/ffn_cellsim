# Port ledger — ALEPH-PORT-303 · phase-ordered dispatch with evaluation-proof coverage

| Field | Value |
|---|---|
| Aleph target | `aleph/runtime/pipeline.py` |
| Lane | L3 runtime |
| Written | 2026-07-30 (before the code landed, per PLAN §0.2.5) |
| Status | LANDED |
| Evidence class | STRUCTURAL. This proves dispatch, not physics. |
| Relation to PLAN | §1.1 row 5, the single most important structural finding of the audit |

## What was read in `ffn_cellsim` (read-only)

- `ffn_sim/ac/engine/dispatch.py` — `CandidatePhase` (7 phases), `CandidateTask`
  (`runtime` + `method_name` + `component_claims` + `connector_claims`), `CanonicalFacadeClaim`, and
  the exact-once coverage gate at `dispatch.py:296-302`.
- `ffn_sim/ac/engine/composed_native.py:296` — `_excluded_facade_noop`, whose own docstring says it
  "launches nothing" and exists "solely so ... the exact-once dispatch-coverage gate still fires".
- `ffn_sim/ac/engine/dispatch.py:105-128` — `membrane_cortex_contact` and four other connector claims
  assigned to `SurfaceBody.accumulate_mechanics`. One of the comments there is explicit that a claim
  is "DECLARED ONLY — no exterior solve exists behind it".

## The fault, stated precisely

The reference gate computes, at build time, the union of `connector_claims` over all scheduled tasks
and requires it to equal the architecture's connector set exactly once. Every input to that
computation is a *string written next to a task*. Therefore the gate is sound with respect to the
proposition "every connector has an owner" and says nothing whatsoever about the proposition "every
connector is evaluated". Two independent mechanisms make the gap reachable:

1. A claim can name a method whose body has no slot for it (`membrane_cortex_contact` →
   `SurfaceBody.accumulate_mechanics`, whose slots are at `surface_body.py:194,364-389`).
2. A claim can name a method that is deliberately empty (`_excluded_facade_noop`).

Both pass. Neither computes a force. Downstream artifacts carry the missing term silently.

## What Aleph does instead — three independent conditions

Re-derived from scratch; no reference code, no reference identifiers, no reference phase names.

| Condition | Mechanism | Defeats |
|---|---|---|
| It was **scheduled** | Registration and binding are separate acts. `assert_schedule_complete()` compares them. | A connector declared but never put in a phase. |
| It was **reached** | The dispatcher writes an `InvocationReceipt` for every bound slot it runs. | A slot that exists but is never executed. |
| It **evaluated something** | The receipt brackets the dispatch with `StepContext.witness_count` (backend op launches + RNG draws). A slot whose witness delta is zero is a `NO_WITNESS` fault. | `_excluded_facade_noop`. A function that launches nothing cannot move a counter it does not touch. |

`assert_schedule_complete()` is documented in the source as **necessary but not sufficient**, with an
explicit note that it is the check the reference had and the check that passed for a connector that is
never evaluated. Keeping it, and labelling it, is more useful than deleting it.

The escape hatch is `Inclusion.EXCLUDED`: it requires a written reason, it is not counted as covered,
and it **cannot be bound into the schedule at all**. The reference project's third option — a stub in
the schedule that satisfies the gate — is not expressible.

## Controls

- Positive control: `test_fully_wired_step_passes_coverage`.
- Negative control (must fail): `test_registered_connector_never_invoked_raises_coverage_violation`.
  In full, the deliberately failing set — this is the lane's reason to exist:
  - `test_registered_connector_never_invoked_raises_coverage_violation` — the direct analogue of
    `membrane_cortex_contact`.
  - `test_noop_placeholder_cannot_satisfy_coverage` — the direct analogue of `_excluded_facade_noop`:
    a stub that returns a fully-populated, entirely truthful-looking `WorkReceipt` and launches
    nothing, rejected by the witness check.
  - `test_static_schedule_gate_alone_is_insufficient` — asserts, in code, that the build-time gate
    passes for the very stub the runtime gate rejects. This test is the audit finding, executable.
  - `test_excluded_participant_cannot_be_scheduled`,
    `test_excluded_registration_requires_a_written_reason`,
    `test_double_binding_in_one_phase_is_refused`.

## Provenance and envelope

- **Source repository / commit read:** `ffn_cellsim` at commit `be0e5876`, read-only.
- **Discarded prose:** all original comments and docstrings discarded. The reference's own docstring
  for `_excluded_facade_noop` is quoted once, above, as *evidence of the defect* — it is cited in the
  ledger, never carried into Aleph source.
- **Units, domains, invariants:** dimensionless. Invariants: phases dispatch in `PHASE_ORDER`
  regardless of binding order; one `(name, phase)` slot may be bound at most once; an object
  registered under two names is snapshotted, committed and rolled back exactly once; an `EXCLUDED`
  registration has an empty required-phase set and cannot be bound.
- **Numerical and precision envelope:** integer counters only. The witness is an exact integer
  difference, so there is no tolerance and no rounding to reason about — deliberately, because a
  coverage gate that needed a tolerance would be a coverage gate that could be tuned into passing.
- **Residency:** host. The witness reads a backend counter, so it works unchanged against a device
  backend that counts its own launches.

## Independence

No `ffn_cellsim` import. Phase names, fault codes, receipt types and the witness mechanism have no
counterpart in the reference. The one thing carried across is the *finding*, which is recorded in
`PLAN.md` §1.1 and restated in this file in Aleph's own words.
