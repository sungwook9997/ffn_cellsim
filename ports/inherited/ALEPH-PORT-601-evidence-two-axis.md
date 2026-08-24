# Port ledger entry ALEPH-PORT-601 — two-axis evidence labelling

| Field | Value |
|---|---|
| Lane | L6 (provenance, artifacts, evidence) |
| Aleph target | `aleph/evidence/ladder.py` |
| Reference read | `/Users/sw1/ffn_cellsim/ffn_sim/ac/engine/contracts.py` lines 50–272 |
| Classification | **CONCEPT_ONLY — not a code port** |
| Written | 2026-07-30, before any code landed |

## What was taken

One idea, and nothing else: *structural progress and quotability are two axes, not one ladder.*
The reference arrived at this after a driver stamped a rung as a string constant and after an
artifact was hand-written as `MECHANISM PASS` + `quantitative_claim_status: BLOCKED` because no
vocabulary existed for the honest case.

## What was NOT taken

- The reference's ten-member rung vocabulary. Aleph's scope (ALEPH-DQ-101, membrane–cortex–pressure)
  has no `SEAMED` / `KERNEL_BOUND` / `OPTIMISED` distinction to make, and inheriting rungs Aleph
  cannot earn would recreate the registry-ahead-of-physics defect the audit found (PLAN §1.2).
- The `StrEnum` + side-table rank encoding. Aleph uses an `IntEnum` with sparse ranks so the
  ordering *is* the type and cannot drift from a parallel dict.
- Every docstring, comment, and rationale string. All prose here is Aleph's.
- The reference's `LADDER_FLOOR` concept of "below-ladder wiring markers". Aleph's `ASSUMED` is on
  the ladder as its floor; it is a rung, not a dump marker.

## Re-derivation

Aleph's rungs were chosen from the *execution facts Aleph can actually establish* for the ratified
initial scope, weakest to strongest:

| Rung | The fact it asserts |
|---|---|
| `ASSUMED` | A human wrote it down. No execution. |
| `ANALYTIC_ORACLE` | A closed form was evaluated. No simulation ran. |
| `CPU_UNIT` | Aleph code executed on CPU at unit scale. |
| `GPU_UNIT` | The same executed on the ratified backend at unit scale. |
| `GPU_NATIVE_SLICE` | Ran at native population, but the step was not accepted. |
| `GPU_NATIVE_ACCEPTED` | Ran at native population and the transaction accepted the step. |

The orthogonal axis is `QuantitativeStatus` — `BLOCKED` / `PROVISIONAL` / `QUOTABLE`. The invariant
enforced in code and proved by test is that **no rung, however high, moves the quantitative axis**.

## Controls shipped

- Positive: `test_classify_claim_quotable_with_artifact`.
- Negative (must refuse): `test_quotable_without_artifact_is_refused`,
  `test_no_rung_implies_quotable`, `test_blocking_reason_overrides_requested_status`.

## Independence

Imports `enum` and `dataclasses` only. Passes with `ffn_cellsim` absent from `sys.path`.
