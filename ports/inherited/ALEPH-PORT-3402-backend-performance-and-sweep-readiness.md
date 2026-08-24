# ALEPH-PORT-3402 — backend performance measurement and sweep readiness

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3402` |
| Lane | `GPU/performance — backend performance and sweep budgeting` |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Zero lines, zero identifiers, zero constants taken from the reference. |
| Authorises | `aleph/runtime/perf.py` |
| Controls | `tests/runtime/test_perf.py` |

---

## 1. Aleph API

```python
from aleph.runtime.perf import (
    cross3, scatter_add_rows,
    assert_bit_identical_cross, assert_bit_identical_scatter,
    OrderSensitivity, measure_order_sensitivity,
    HotSpot, RelaxationProfile, profile_relaxation,
    SweepBudget, TessellationCost, budget_sweep,
)
```

This module is **measurement and bit-exact substitution only**. It owns no physical state, applies
no force, and is imported by nothing in `aleph/vertical/**`. It is a tool the sweep scripts and the
GPU lane use to decide *what to optimise* and to *prove an optimisation changed nothing*.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY, never modified) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Nearest source paths | none — the reference has no parity, ULP-envelope, or sweep-budget module |
| Lines taken | **0** |

The reference was consulted for a performance-measurement harness and has none. Its profiling is ad
hoc timing inside run scripts, not a module with controls. There was nothing to port.

## 3. Why not clean-room — and why the answer is "clean-room won"

**RE-DERIVED. Nothing was ported.** Every function here is derived from two facts that are in this
repository rather than in any reference:

- `aleph/runtime/warp_kernels.py`'s docstring states that the atomic scatter path's non-determinism
  "should be a measured quantity rather than a rumour" and nobody had measured it.
- `PLAN.md` §2.5 records that step count depends on the *tessellation*, not the resolution — an
  icosphere at V=162 converges in 3,426 steps where a Fibonacci hull at the same vertex count needs
  189,419.

## 4. The law represented, and the derivation

Three separate identities, each independently checkable.

**(a) The cross product of two 3-vectors is a closed-form expression.**
`c₀ = a₁b₂ − a₂b₁`, `c₁ = a₂b₀ − a₀b₂`, `c₂ = a₀b₁ − a₁b₀`. `numpy.cross` evaluates exactly this
expression, in exactly this association, for the `(N,3) × (N,3)` case. Writing the three lines out
therefore produces the *same sequence of floating-point operations* and is bit-identical by
derivation, not by luck. What it removes is `numpy.cross`'s axis-normalisation and `moveaxis`
bookkeeping, which is per-call overhead independent of `N`.

**(b) A row-grouped sum is order-determined, and two implementations that visit contributions in
the same order agree bitwise.** `numpy.add.at(dest, idx, vals)` adds contributions to each
destination row in ascending contribution index. `numpy.bincount(idx, weights=w)` accumulates by
sweeping `w` once in index order. Same order, same float64 accumulator, therefore same bits. The
speedup comes entirely from `bincount` being a compiled single pass where `ufunc.at` is the
unbuffered generic path.

**(c) A varying summation order perturbs a float sum by a bounded amount.** For `m` contributions
to one row, reassociating the sum moves the result by at most about `(m−1)·u·Σ|contributions|`.
This is what an atomic scatter costs on a device where thread order varies. It is *measurable on a
CPU* by permuting the contribution order explicitly — which is the substance of §7.

## 5. Units, domains, singular cases, invariants

| Quantity | Unit | Domain |
|---|---|---|
| `cross3` operands | any | `(..., 3)`, last axis exactly 3 |
| `scatter_add_rows` indices | dimensionless | `[0, rows)`, integer |
| `OrderSensitivity.ulp_spread` | float32 ULP | `≥ 0` |
| `RelaxationProfile` timings | seconds | `≥ 0` |
| `SweepBudget` totals | seconds | `≥ 0` |

**Singular cases handled explicitly rather than by accident**: zero-length arrays (`n = 0`), a
single destination row, an index array that names one row only, and an empty contribution set. Each
has a control.

**The invariant this module exists to protect**: a substitution offered here is *bit-identical* to
the operation it replaces, or it is not offered. `assert_bit_identical_*` are the checks, and they
are the point of the module rather than decoration — `CLAUDE.md` §3 forbids a speedup that changes
a residual, and the only way to know is to compare bits.

## 6. Evidence class and known retractions

| Claim | Evidence class |
|---|---|
| `cross3` is bit-identical to `numpy.cross` | **MEASURED** — 400 randomised cases, float32 and float64 |
| `scatter_add_rows` is bit-identical to `numpy.add.at` | **MEASURED** — 800 randomised cases across 4 adversarial index patterns |
| Ordered warp scatter is deterministic on CPU | **MEASURED** — 16 repeats, zero spread |
| Atomic warp scatter is deterministic **on CUDA** | **UNVERIFIED** — cannot be measured without a GPU; see §7 |
| Order-sensitivity envelope 2–5 float32 ULP | **MEASURED on CPU**, and it is a *proxy*, not a CUDA observation |
| The 55× tessellation factor | **QUOTED** from `PLAN.md` §2.5; the icosphere half (3,426) is re-measured here |

**No retraction.** One correction this entry makes to the record: `docs/design/GPU_POLICY.md` and
`HANDOFF.md` both state that warp is not installed and that `aleph/` contains no `@wp.kernel`
definitions. Both were true when written and are now false — warp 1.15.0 is installed and
`warp_kernels.py` defines 18 kernels.

## 7. Independent oracle

Each function is checked against an oracle that does not share its implementation:

| Function | Oracle |
|---|---|
| `cross3` | `numpy.cross` — bitwise equality, not a tolerance |
| `scatter_add_rows` | `numpy.add.at` — bitwise equality, not a tolerance |
| `measure_order_sensitivity` | an explicitly permuted `numpy.add.at`, which *constructs* the reassociation a GPU atomic would produce by accident |
| `profile_relaxation` | wall-clock totals must sum to within the measured whole; a profile that does not account for its own runtime is rejected |

**The oracle for the un-measurable claim.** Atomic determinism on CUDA cannot be observed here.
Rather than leave it unknown, `measure_order_sensitivity` bounds it: it takes the *same*
contributions, permutes their order explicitly, and reports the spread. A CUDA atomic scatter cannot
produce a spread larger than the worst reassociation, so the CPU-measured envelope is an upper bound
on the CUDA quantity. This converts "unknown until a GPU is authorised" into "at most N ULP, and
here is N".

## 8. Positive control

`tests/runtime/test_perf.py::test_cross3_is_bit_identical_to_numpy_cross` and
`::test_scatter_add_rows_is_bit_identical_to_add_at` — randomised, multi-dtype, multi-pattern.
`::test_order_sensitivity_detects_a_real_reassociation` asserts the envelope is strictly positive on
input constructed to cancel, proving the measurement responds to the thing it measures.

## 9. Negative control — the deliberately failing kind

A measurement that cannot fail is not a measurement. Three counterexamples:

- `::test_bit_identity_checker_rejects_a_deliberately_wrong_cross` — a cross product with two
  components transposed must be *caught*, not tolerated.
- `::test_bit_identity_checker_rejects_a_reordered_scatter` — a scatter that sums in a permuted
  order is bit-different and must be reported as such, even though it is physically equivalent and
  agrees to 5 ULP. **This is the control that matters**: it proves the checker tests bits and not a
  tolerance, so it would catch exactly the class of silent change `CLAUDE.md` §3 forbids.
- `::test_order_sensitivity_is_zero_for_a_single_contribution_per_row` — with nothing to reassociate
  the envelope must be exactly zero, so a nonzero reading is a defect in the harness.

## 10. Numerical and precision envelope

Round-off, tolerance and ULP are the subject of this module rather than a caveat on it.

- Bit-identity claims are `array_equal`, **tolerance zero**. There is no tolerance to loosen.
- The order-sensitivity envelope is reported in **float32 ULP** relative to a constituent scale
  (`Σ|contributions|` per row), never relative to the resultant. A residual reported against the
  resultant looks catastrophic wherever the resultant cancels — the `subtract[cancel]` row in
  `aleph/runtime/parity.py`'s report reads `5.4e-02` against the resultant and `0.45 u32` against
  the constituents, and only the second number means anything.
- float64 accumulation is preserved everywhere the replaced operation used it. `bincount`'s
  `weights` accumulator is float64 regardless of the input dtype, which is what makes (b) hold.

## 11. Production-backend residency and transfer

**Host-side, CPU, numpy only.** This module imports no accelerator, allocates on no device, and
performs no host-to-device transfer. It is subject to
`tests/runtime/test_import_boundaries.py::test_no_module_level_accelerator_import` like every other
`aleph/runtime/*.py`, and it satisfies it trivially by never importing warp at all.

The GPU relevance is indirect and is the point: `profile_relaxation` identifies which operations a
`WarpBackend` substitution would actually accelerate, and `measure_order_sensitivity` bounds what
the atomic path would cost in reproducibility *before* a GPU is authorised, so that the first
authorised GPU hour is spent confirming a prediction rather than discovering a question.

## 12. Comments and docstrings to discard

None — no source prose was read, so none was carried. Every comment and docstring in the module is
new prose written for this entry. Where a number appears in the module it is either a measured value
with its measurement named, or a documented default with the reason for the default.

## 13. Acceptance result, reviewer, rollback

| Field | Value |
|---|---|
| Acceptance result | Status is `PROPOSED`. The controls run and pass; the entry is not `ACCEPTED` because no reviewer other than its author has read it. |
| Reviewer | none yet — **this entry has not been reviewed** |
| Rollback | Delete `aleph/runtime/perf.py` and `tests/runtime/test_perf.py`. Nothing imports the module, so removal is a two-file delete with no callers to repair. |

**Rollback is genuinely free and that is deliberate.** The module is additive and unreferenced by
the physics; it can be removed without touching a single line of `aleph/vertical/**`.
