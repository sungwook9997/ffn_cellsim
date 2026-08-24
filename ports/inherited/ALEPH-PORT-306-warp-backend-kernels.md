# ALEPH-PORT-306 — Warp kernels, ordered scatter-add, and the backend parity harness

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-306` |
| Lane | `L22 accelerator backend` |
| Status | `ACCEPTED` |
| Acceptance scope | **CPU-Warp only.** See §10 — the CUDA path is written, guarded, and unexercised. No kernel has ever run on a GPU. |
| Written | 2026-07-30 (entry written before the kernels landed, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` (plus one `DEFECT_STUDY` finding, §3) |
| Evidence class | `STRUCTURAL` — a substitution measured against the CPU reference. Not physics, and not a GPU result. |

---

## 1. Aleph API

```python
from aleph.runtime.backend import ScatterMode, WarpBackend      # extended; guard untouched
from aleph.runtime.warp_kernels import CHUNK, load, loaded_warp
from aleph.runtime.parity import (
    EPS32, U32, OpCase, OpParity, ParityReport, ScatterDeterminism,
    default_cases, measure_scatter_determinism, parity_report,
    perturbed_backend, report_as_dict,
)
```

Aleph targets: `aleph/runtime/warp_kernels.py`, `aleph/runtime/parity.py`, `aleph/runtime/backend.py`,
`scripts/gpu_smoke.py`.

`require_gpu_authorization`, `check_gpu_authorization` and `write_authorization_record` are **not**
touched by this entry. The one change to the authorization surface is the addition of
`WarpBackend(device="cpu")`, which is described in §7 and is flagged there as needing the PI's
ratification rather than presented as settled.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`), 2026-07-29 23:57:26 +0900 |
| Source path | `ffn_sim/dcm/radial_shell_warp.py`, `ffn_sim/dcm/baoab_warp.py`, `ffn_sim/dcm/dcm_cadherin_gpu.py` |
| Source symbol(s) | the module docstrings' parity contracts; `wp.atomic_add` reduction and stream-compaction sites |
| Read from | working tree |
| Working tree == commit? | **no** — the working tree carries 31 uncommitted changes at `be0e5876`. What was read is the working tree, and that is what this row records. |

Nothing was read for the purpose of copying. The kernels in Aleph are original; what was taken from
the reference is a **negative** finding about the shape of its parity discipline (§3), and that
finding shaped Aleph's design decision in the opposite direction.

## 3. Why source-derived porting beats clean-room — it does not, and the reason is the finding

Clean-room wins outright. Every kernel here is a few lines of arithmetic that any competent author
derives directly from the operation's definition; there is no empirically discovered failure mode and
no delicate case analysis to inherit. So no code was ported.

The reference *was* worth reading, for one thing, and it is a `DEFECT_STUDY`:

`ffn_sim/dcm/radial_shell_warp.py` declares its own parity contract in prose, and the contract says
that the path the gate runs is not the path the full port runs:

- `reduce="host"` is described as **"the GATED path"**: the reduction happens in numpy, only the
  per-bead force law runs in the kernel, and it achieves ~1e-13 relative agreement.
- `reduce="warp"` is described as **"the FULL-port diagnostic"**, reduced with `wp.atomic_add`, and is
  "graded only at the .cu's own reduction-order tolerance and reported as a diagnostic, **not the
  committed bit-parity claim**".

The reference is honest about this in its prose, and the structural consequence still stands: the
configuration that carries the committed claim is not the configuration a full GPU port would run. It
is the same shape as the WCA finding recorded elsewhere in this ledger — the parity gate runs uncapped
while production runs capped — and in both cases the gate is green while the thing it certifies is not
the thing that executes.

**Aleph's response is not a tighter tolerance. It is to remove the gap.** Rather than gate a
host-reduced path and ship an atomic one, `ScatterMode.ORDERED` makes the *default and production*
path deterministic by construction, so the configuration the parity table measures is the
configuration that runs. The atomic path still exists, because a fast path that is never measured is
worse than one that is — but it is opt-in and it is measured
(`measure_scatter_determinism`), and its result carries a caveat field that refuses to let a CPU
measurement be read as evidence about CUDA.

## 4. Physical or mathematical law represented

None. This is a compute substrate, and saying so plainly matters: nothing in this entry is a physical
claim, and no number produced by these kernels is evidence about a membrane.

What is re-derived here is **numerical**, and there are two derivations worth stating because both
are load-bearing.

**(a) Why float64 accumulation changes the accounting, and where the error actually comes from.**
Summing `n` values in floating point with a running accumulator has a worst-case error bound
`≈ n·u·Σ|xᵢ|`, where `u` is the accumulator's unit round-off. In float32 (`u ≈ 6e-8`) at `n = 10⁶`
that bound exceeds the sum itself for a cancelling sum, which is exactly the case a force-closure
residual is. Widening only the *accumulator* to float64 (`u ≈ 1.1e-16`) drops that term by eight orders
of magnitude and leaves a different error dominant: the one-time rounding of each *input* into
float32, which is `u₃₂·Σ|xᵢ|` and cannot be removed without storing float64. So the residual error is
one rounding per input rather than one per addition, and the conditioning factor `Σ|x|/|Σx|` is what
converts that into a relative bound. `aleph/runtime/parity.py` computes that factor from the data and
reports it as the `amp` column, which is why the bound widens on an ill-conditioned sum instead of
failing on correct arithmetic.

**(b) Why the ordered scatter is bit-reproducible, and why that is not merely tidiness.**
Float addition is not associative, so a sum's value depends on its order. An atomic scatter fixes the
*multiset* of contributions but not their order, so two identical launches can differ in the last
bits. The fix is to fix the order. Grouping contributions by destination row (a CSR built by a
**stable** sort on the host) and giving one thread sole ownership of one destination row means no two
threads write the same element, no atomic is involved, and the addition sequence is a function of the
index array alone. Because the sort is stable, that sequence is ascending contribution index — which
is also the order `numpy.add.at` uses, so the ordered mode agrees with the CPU reference **bit for
bit** in float64 rather than merely to round-off. `test_ordered_scatter_reproduces_the_references_summation_order_exactly`
pins that, and it is deliberately an equality rather than an `approx`: if it ever has to be loosened,
the ordering guarantee has been lost and the determinism argument has gone with it.

## 5. Units, domains, singular cases, invariants

The backend is dimensionless; units belong to the caller. Domain facts and invariants this entry is
responsible for:

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| every array element | dimensionless (caller's) | dimensionless | float32 or float64 only; any other dtype raises `TypeError` |
| `indices` in `scatter_add` | index | — | `0 ≤ i < dest.shape[0]`; outside raises `IndexError` |
| `CHUNK` | elements/thread | — | constant `256`; the block count is a pure function of the element count |

Invariants, each with a test:

- `op_count` advances identically to `NumpyBackend` for the same call sequence, and is **not** advanced
  by `to_host` or `synchronize`. The dispatch coverage gate reads this counter, so two backends that
  disagreed about what a step cost would make the same schedule pass on one and fail on the other.
- Reductions return float64 **outputs**, not just float64 accumulators. Narrowing the output would
  discard the widening on the way out.
- `scatter_add` accumulates into `dest`; it never assumes a zeroed buffer.
- An empty contribution set is a no-op, not an error.
- A CUDA request never degrades to the CPU device, and a CPU request never reaches CUDA. Both
  directions refuse.

Singular cases handled explicitly: zero-length arrays (`sum`/`norm` → `0.0`, `max_abs` → `0.0`,
matching the reference backend); a reduction axis the kernels do not implement **refuses** rather than
guessing, because a silently transposed reduction is a plausible-looking wrong answer.

## 6. Numerical and precision envelope

`ALEPH-DQ-107` made concrete: float32 compute, float64 accumulation where the accounting needs it.

| Channel | dtype | Members |
|---|---|---|
| compute | float32 (default) | `zeros`, `array`, `add`, `subtract`, `multiply`, `scale` |
| accounting | float64 always | every reduction accumulator and output; the row accumulator in the ordered scatter |

No tolerance is chosen at this layer. The parity harness reports a **measured** max-absolute and
max-relative difference per operation, together with a round-off bound `u₃₂·(amp+1)` derived from the
data (§4a), where `u₃₂ = 5.96e-08`. Relative difference is defined as
`max|candidate − reference| / max|reference|`; the denominator is the array's scale rather than each
element's own value, because elementwise relative error is unusable on a force array whose components
legitimately pass through zero.

Measured on CPU-Warp (warp-lang 1.15.0, macOS arm64, float32 compute vs the float64 reference), all 29
default cases sit inside their bound. Worst ratios of measured to bound: `multiply` 0.73, `subtract`
0.64, `scale` 0.61. The `subtract[cancel]` case carries `amp ≈ 2e6` and a correspondingly wide bound —
that case exists so the `amp` column is demonstrably doing work rather than reading 1.00 everywhere.
Observed reduction errors are far *inside* the worst-case bound (e.g. `sum` at n=250,000: measured
1.15e-08 against a bound of 7.70e-05), because the input roundings are independent and accumulate like
`√n` rather than `n`. The bound is worst-case and is not tightened to match the observation, because a
bound fitted to the runs that motivated it is not a bound.

Scatter-add reproducibility, measured over 8 repeats of 24,576 contributions into 4,096 rows:
`ORDERED` bit-identical, spread exactly `0.0`. `ATOMIC` also bit-identical **on the CPU device**, and
that number is reported with a caveat stating that a warp CPU launch is a serial loop, so it is a
property of that schedule and is **not** evidence that atomics are deterministic on CUDA.

## 7. Residency and transfer

Device-resident throughout: allocation, elementwise arithmetic, reductions and both scatter paths run
on the warp device. Host-side, deliberately, and each for a stated reason:

- **The CSR build** (`np.argsort(..., kind="stable")`, `np.bincount`, `np.cumsum`) for the ordered
  scatter. This is the cost of the determinism decision. It is cached by the content of the index
  array (blake2b-16, bounded to 8 entries so it cannot grow into device memory), and for Aleph it
  amortizes to nothing because the index set is the mesh connectivity and is constant across the steps
  of a run.
- **The final combine of reduction partials.** `⌈n/256⌉` float64 values are read back and summed in a
  fixed order. A device-side tree reduction would have been faster and would have made the result
  depend on the device's block count, which is the determinism property traded away.

`to_host` is the only transfer in the Protocol and is uncounted, because a host read is not evaluation
and a participant must not be able to manufacture an evaluation witness by reading.

**The `device="cpu"` carve-out, stated as a decision rather than buried.** `WarpBackend(device=None)`
— the default — means a CUDA device and is fully guarded: `require_gpu_authorization` runs before
`warp` is imported. `WarpBackend(device="cpu")` selects Warp's CPU device, touches no GPU, and does
not require a GPU grant. Three properties keep that from being a bypass: it is never a default and
never a fallback; the CUDA path raises `BackendUnavailable` rather than degrading to it; and a
CPU-mode backend asserts after device resolution that the device is not CUDA and refuses if it is.
`require_gpu_authorization` itself is unchanged. **This is an agent proposal and wants the PI's
ratification** — it is what makes CPU-Warp parity measurable on a laptop, and it is the only
judgement call in this entry that the written policy does not already settle.

## 8. Comments and docstrings discarded

No comment, docstring, or identifier from the reference survives, because no code was taken. The one
piece of reference *prose* that mattered — `radial_shell_warp.py`'s parity contract — is quoted only
here, in the ledger, where naming the source is required; it appears nowhere in `aleph/**`, where
naming it is forbidden. All prose in the landed modules is new and written for Aleph's situation.
`tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package` enforces
this mechanically and passes.

## 9. Independent oracle and controls

The oracle is `NumpyBackend`, which is Aleph's own and is the *definition* of every operation in the
Protocol — not a second opinion about it. Correctness is additionally separated from precision:
correctness tests run **float64 on both backends**, so a disagreement beyond float64 round-off is a
wrong kernel rather than a narrow mantissa. A float32-vs-float64 comparison cannot tell those apart,
and an indexing bug hiding under "well, it's float32" is exactly what survives to production.

| Control | Test | Result |
|---|---|---|
| Positive | `tests/runtime/test_warp_backend.py::test_elementwise_matches_the_reference` | exact equality, float64 both sides |
| Positive | `tests/runtime/test_warp_backend.py::test_ordered_scatter_reproduces_the_references_summation_order_exactly` | bit-identical to `numpy.add.at` |
| Positive | `tests/runtime/test_parity.py::test_cpu_warp_is_within_float32_roundoff_on_every_op` | all 29 cases inside bound |
| Positive | `tests/runtime/test_warp_backend.py::test_the_counter_advances_exactly_as_the_reference_does` | both backends +10 |
| Negative (must fail) | `tests/runtime/test_parity.py::test_a_biased_op_is_caught` | a 1e-5 relative bias on any of 7 ops is flagged |
| Negative (must fail) | `tests/runtime/test_parity.py::test_a_dropped_scatter_contribution_is_caught` | one dropped contribution is flagged |
| Negative (must fail) | `tests/runtime/test_parity.py::test_a_perturbed_warp_backend_is_still_caught` | the control runs against the **shipping** backend, not only against numpy |
| Negative (must fail) | `tests/runtime/test_gpu_smoke.py::test_a_parity_failure_exits_four` | the script's failure exit path is reachable |
| Negative (must fail) | `tests/runtime/test_warp_backend.py::test_a_cuda_request_never_degrades_to_the_cpu` | a GPU request with no GPU refuses |
| Specificity | `tests/runtime/test_parity.py::test_only_the_perturbed_op_is_flagged` | exactly one op flagged, not all |
| Anti-over-fitting | `tests/runtime/test_parity.py::test_a_perturbation_below_the_bound_is_not_flagged` | 1e-9 relative is correctly *not* flagged |

The last two rows are the ones that keep the gate honest in both directions. A harness that flagged
everything would "catch" the bug and locate nothing, and the repair applied under time pressure would
be to loosen the bound until nothing failed — which is how a gate stops working.

`test_a_perturbed_warp_backend_is_still_caught` is the direct answer to §3: the negative control is run
on the configuration that actually ships.

## 10. Acceptance result, reviewer, rollback

**Acceptance is scoped to the CPU-Warp substitution and to the structural properties above. It is not
a GPU result.** No CUDA device was reached and no kernel has ever executed on one, because no GPU
authorization record exists for any host and an agent may not create one. Every GPU-requiring test
skips with a reason naming the hardware (`RTX A5000 Laptop` on `gbook`, or an `rtx4090_0`/`rtx3090`
Slurm license on `gpu-sungwook`) and the grant it lacks — never silently.

Unverified, and listed so it is not mistaken for tested:

- kernel compilation for `sm_86` / `sm_89`, and CUDA kernel launch of any kind;
- atomic-mode non-determinism magnitude under real thread contention (the CPU measurement says
  nothing about it, and its caveat field says so);
- device memory reporting, `compute_capability`, and free-memory numbers on a real device;
- whether the ordered scatter's determinism survives to CUDA. It is deterministic *by construction*
  and there is no mechanism by which it should not, but it has not been observed there.

| Field | Value |
|---|---|
| Reviewer | none — agent-authored, awaiting PI review |
| Verified by | `tests/runtime` 189 passed / 3 skipped (GPU); full suite green |
| Rollback | delete `aleph/runtime/warp_kernels.py`, `aleph/runtime/parity.py`, `scripts/gpu_smoke.py` and revert the `WarpBackend` body to its refusal stub. `NumpyBackend`, the Protocol and the guard are untouched, so nothing else depends on this. |
| Open question for the PI | ratify or overturn the `device="cpu"` carve-out (§7) |
