# ALEPH-PORT-3634 — the backend seam, the silent write it hid, and what a crossing costs

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3634` |
| Lane | `46143f30` Lane W2 (lead) |
| Status | `PROPOSED` |
| Written | `2026-08-05` — **AFTER the code, and that is a lapse.** `CLAUDE.md` §3 requires the ledger first. It is recorded as a repair rather than backdated, because a ledger that claims to have preceded code it followed is worth less than no ledger. |
| Port class | `NOT A PORT` — nothing is taken from `/Users/sw1/ffn_cellsim`; no law changes and no physics is derived. This is the substrate seam, not a law. |
| Aleph target | `aleph/runtime/backend.py` (`WarpBackend` only), `aleph/scenarios/world.py` and `aleph/scenarios/whole_cell.py` (one defaulted parameter each), `tests/scenarios/test_whole_cell_backend_parity.py` (new) |
| Exists because | The whole cell could not be handed a `WarpBackend` at all. It now can, it produces **the same trajectory to the bit**, and the price of doing it that way is a counted number rather than an impression. |

---

## 0. Scope claim, because this file is not in this lane's registry row

`aleph/runtime/backend.py` was owned by the `5ad3a7ba` (S2) row, which is **CLOSED** with its work
uncommitted in the working tree. No live row claims it: the `e874a7fe` (Lane G10) row names it under
*will not touch* and defers it to S2. So this is a **takeover of a closed row's path**, not a
collision, and `docs/ACTIVE_SESSIONS.md` carries the corresponding note. Recording it here as well
because a scope decision that lives in only one file is a scope decision nobody can check.

## 1. What did not work, in three layers

`build_whole_cell(backend=WarpBackend(...))` failed at step time, and fixing each failure exposed
the next. The order matters because the third was invisible and the first two were not.

| # | symptom | cause |
|---|---|---|
| 1 | `TypeError: Unable to infer the type of argument 'a', got numpy.ndarray` | `_flat` handed host arrays straight to `wp.launch`. The message named the argument and not the boundary it had crossed. |
| 2 | `ValueError: shape mismatch: (30, 3) vs (30, 1)` | callers rely on **NumPy broadcasting**, which `NumpyBackend` inherits from `np.multiply` for free and the `Backend` Protocol never promised. |
| 3 | *nothing* — the run completed and the answer was wrong | see §2. |

## 2. The defect that raised nothing

Fixing (1) by uploading host arrays inside `_flat` was right for **inputs**, which are only read.
It was wrong for `out`, which is written:

    scale(a, factor, out=<host array>)
        -> _flat(out) uploads a COPY of `out` to the device
        -> the kernel writes the result into the copy
        -> the copy is dropped; the caller's array still holds its old values
        -> the method returns the right shape and raises nothing

**A leniency that removes an exception can remove what the exception was protecting.** Before the
change, a host `out` was fatal. After it, a host `out` was a no-op.

`_elementwise` had a guard for this. `scale` did not, because `scale` builds its own `wp.launch`
for the scalar. One omission in one method.

### What it cost, measured

| after | relative energy gap, numpy vs warp-cpu |
|---|---|
| 1 step | `0.000e+00` — bit-identical |
| 2 steps | `8.036e-04` |
| 3 steps | `2.408e-03` |
| 5 steps | `8.014e-03` |

Per-owner after 2 steps: `microtubule` 1.53e-3 µm, `lamellipodium` 1.33e-3, `filopodium` 9.4e-4,
`sf_arc` 6.4e-4, `ecm` 6.3e-5, `membrane` 1.3e-5, `nucleus` 4.8e-22. **`cortex` and `cytosol` were
bit-identical**, which is why a whole-cell scalar would have been a poor detector: the largest
owner in the cell agreed perfectly while four others drifted.

### The two hypotheses that were wrong, and the controls that killed them

Recorded because they are the reason the search took as long as it did, and because each was
plausible enough to have been written up as the answer.

1. **"An iterative solver's stopping rule."** The cytosol's Biot step is MINRES with a relative
   tolerance, and a tolerance-level disagreement between substrates is ordinary. *Refuted:*
   tightening `SOLVE_RELATIVE_TOLERANCE` from `1e-14` through `1e-16` to `1e-18` moved the gap by
   **zero digits** — `8.036e-04` all three times. A gap that does not respond to the tolerance is
   not made of tolerance.
2. **"Round-off amplified by a stiff, non-smooth model."** The full reductions genuinely do differ
   between substrates (`sum` 2.3e-14, `dot` 7.1e-15 absolute — chunked parallel reduction has a
   different summation order and cannot not), so amplification was the comfortable story.
   *Refuted:* perturbing one owner's positions by `1e-12` **relative** and running the same two
   steps on NumPy alone moved the final energy by `8.8e-16` relative. The model does not amplify.
   A `8e-4` gap therefore could not be made of `1e-14` seeds, and had to be a wrong answer.

The second control is the one that mattered. Without it, "the system is sensitive" would have been
recorded as an explanation, the defect would have stayed, and every future substrate disagreement
would have had a ready excuse waiting for it.

## 3. The change

Four things in `WarpBackend`, and nothing outside it that changes a number:

1. **`_on_device(a)`** — a host array becomes a device array here, via `_upload`, which counts it.
   `_flat` is now `_on_device` plus a flatten, so the two axis-wise reduction kernels (`sum_axis0`,
   `sum_axis1`, `norm_axis1`), which took their operand raw, go through it too.
2. **`_resolve_out(out, shape, dtype)` / `_deliver(target, host_out)`** — the single path an `out`
   may take. A host `out` gets a device temporary and an explicit copy back, and the copy is
   counted. Every launch that accepts an `out` calls both; that is what stops the next one
   repeating `scale`'s omission.
3. **Host-side broadcasting** — when shapes differ, the host operand is expanded with
   `np.broadcast_to` *before* its upload, where the expansion costs nothing extra: it changes what
   is copied, not how many times. A **device** operand needing expansion is refused by name,
   because materialising one needs a kernel this backend does not have and doing it quietly would
   make the seam's cost depend on a shape nobody wrote down.
4. **`_upload` / `_download`** — the only two places this class copies across the bus, and both
   count the crossing and its bytes. `host_uploads`, `host_downloads`, `host_upload_bytes` and
   `host_download_bytes` are the public readings. §5 records that the first version of this counted
   in two *other* places and saw a third of the traffic.

`aleph/scenarios/world.py` and `whole_cell.py` each take `backend: Any = None`, defaulting to
`NumpyBackend()`. **Passing one does not make the owners device-resident.** They still hold host
NumPy; what changes is which substrate the connectors' scatters and the pipeline's witness run on.
That is the seam, and it is not the port.

## 4. Acceptance: `==`, not `approx`

`tests/scenarios/test_whole_cell_backend_parity.py`, 14 controls, all passing:

- 5 steps of the whole cell, `numpy` vs `warp-cpu` float64, **equal energies as `==`**;
- and **every owner's positions equal as `np.array_equal`**, because a scalar can agree by
  cancellation and the defect's signature was per-owner;
- `scale` and the three elementwise operations each write a **host** `out`, over three shapes,
  with an explicit assertion that the array is not still zeros — a dropped device write leaves
  exactly that;
- broadcast `multiply` over three shape pairs against NumPy's own answer, so an expansion along
  the wrong axis fails rather than looking plausible;
- the crossing counters are asserted **non-zero**, which is deliberate and reads oddly on purpose:
  it is the honest statement of where the port is, and the day an owner becomes device-resident
  that control fails and must be rewritten to say so;
- and `test_no_crossing_escapes_the_two_counters` spies on `warp.array.__init__` and
  `warp.array.numpy` themselves, then demands the backend's counters equal what the substrate was
  actually asked to do. **It asserts no fixed number**, so it does not merely pin today's traffic:
  a new uncounted crossing fails it the day it is written. It found two sites beyond the ones the
  audit had listed — `scatter_add`'s index array and `gather`'s — on its first run.

A tolerance anywhere in the first two would have passed the defect: `8e-3` relative after five
steps is inside any `rtol` a reasonable person would have typed.

**Why per-kernel ULP parity could not have caught this.** The `k1` sweep grades each kernel against
its NumPy definition inside a ULP budget and passed 80 of 85 mutants on the RTX 4090 (job 24). It
grades kernels that have *already been handed device arrays*. The defect was in how host arrays
reach them. Every kernel involved was correct; the caller's array was never written. **A gate is
blind to the layer beneath the one it grades**, and that is a general property of this harness
worth stating once here rather than rediscovering.

## 5. The cost, which is the PI's question

Same process, warm Warp kernel cache, 5 steps of the whole cell at level 2:

> **CORRECTED 2026-08-05 ~03:3x. The first version of this table was wrong by 3.3× on uploads,
> and an adversarial audit of this lane's own work is what caught it.** The counters were
> incremented in `_on_device` and `_deliver` only, while `scatter_add`'s values AND its index
> array, `_csr`'s order and offsets, `gather`'s indices, `_partials`' download of its reduction
> partials, and `to_host` all crossed the bus uncounted. Every crossing now goes through one of
> two funnels, `_upload` and `_download`, and a new control (§4) spies on `warp.array.__init__`
> and `warp.array.numpy` and fails if the backend's counters disagree with what the substrate was
> actually asked to do. That control immediately found two more sites than the audit had listed.

| backend | wall | ops | host uploads | host downloads |
|---|---:|---:|---:|---:|
| `numpy` | 0.52 s | 860 | 0 | 0 |
| `warp-cpu` float64 | 0.56 s | 860 | **962 (3.873 MB)** | **275 (1.797 MB)** |
| *first version, wrong* | — | — | *290* | *170* |

**1.08×, and it is slower rather than faster.** That is the expected and correct result for this
state of the port, and it is the measurement rather than an argument: with the owners on the host,
every operation that touches the device pays to get there and pays again to come back. A download
is the worse half — the host cannot read a buffer the device has not finished writing, so each one
**synchronises**.

**Why the undercount mattered more than the number.** The uncounted sites were *structural*, not
residency-dependent: `scatter_add` builds its index array from the host unconditionally, `_csr`
always builds from host NumPy, `_partials` always downloads. So making an owner device-resident
would have driven the published counters to 0/0 while every scatter still round-tripped its
indices and every reduction still downloaded its partials — and the §4 control written to *fail
loudly on exactly that day* asserts only `> 0`, so it would have failed for the wrong reason and
reported an all-clear. A counter that measures a third of its subject is worse than no counter,
because a number gets quoted and an absence gets investigated.

**A number that must not be quoted: an earlier reading of `5.78 s` against `1.45 s` (4.0×).** That
run compiled the Warp kernels inside the timed region. It is one-time cost, it is not the steady
state, and it is recorded here so that if it appears in any earlier note it is known to be wrong.

The PI's stated experience — *"중간의 검측기랑 cpu-gpu host time이 제일 많이 깎아먹었었어"* — is
what the counters were built to make checkable, and on this measurement the host time is real and
the seam is not where the win is.

**Coordinates for that quote, and they are weaker than they should be** (added 2026-08-05 after an
audit found it uncited, per `CLAUDE.md` §2 rule 2). The sentence appears in this lane's transcript
`~/.claude/projects/-Users-sw1-Project-Aleph/46143f30-d1ab-4f8a-9c66-b63204a3b7a4.jsonl` **only inside a compaction summary** — record uuid
`72f47109-8979-451b-889a-5d81e90f1137`, timestamp `2026-08-04T16:50:03.819Z`. The original PI
message predates a context compaction and is **not** in the file as a user record; searching it for
the distinctive fragment `검측기` returns that one summary record and nothing else.

So this is a quote of a *summary of* a PI message, not of the message. It is labelled that way
rather than dressed as a direct citation, because the point of the coordinates rule is that a
session which did not receive the message can check what is there — and what is there is a summary.
Nothing in this entry rests on it: every number in §5 is a measurement, and the quote is context for
why the counters exist rather than authority for anything.

## 6. What this does NOT do

- ~~**It does not port anything to the GPU.**~~ **SUPERSEDED 2026-08-05 ~04:0x — see §13.** The
  whole cell now steps on an RTX 4090 and is bit-identical to NumPy there. Getting it to do so
  required two more fixes, both of which CPU-Warp's leniency had concealed.
- **It does not make any owner device-resident**, which is the actual port and the actual win. The
  counters name the distance to it: 460 crossings per 5 steps is what "the owners are on the host"
  costs, and the target for that number is zero.
- **It does not change any law.** `NumpyBackend` is untouched; the CPU form is still the definition.

## 7. Units, domains, singular cases, invariants

No units enter: this is a substrate, and it moves whatever it is handed.

| singular case | treatment |
|---|---|
| zero-size array | no launch; the allocated target is returned unwritten, matching NumPy |
| `out` is also an input | correct — inputs are uploaded (read) before the target is written |
| `out` shape differs from the broadcast result | `reshape` raises rather than silently ravelling |
| device operand needing broadcast | refused by name (§3.3) |
| an integer index array | goes through `scatter_add`'s own path, which was measured exact and is not touched here |

| invariant | how it holds |
|---|---|
| the substrate does not change the answer | asserted `==` on energies and `array_equal` on every owner's positions |
| an `out` is always the array the caller passed | one `_resolve_out`/`_deliver` path; `scale`'s omission is the reason it is one path and not a convention |
| a crossing is always counted | every host↔device copy goes through `_upload` or `_download`, and `test_no_crossing_escapes_the_two_counters` spies on the substrate itself and fails if the counters disagree. **The first version asserted this invariant naming `_on_device` and `_deliver` as the only two sites, and it was false** — see §5. |

## 8. Known remaining seam gaps, named rather than left

- **`gather` raises on a host `src`** (`Unable to infer the type of argument 'src'`). Its *index*
  array now uploads through `_upload`, but the source array does not, and the whole cell does not
  reach that path. A gap, not a defect today. **Note this entry originally said the whole cell does
  not reach `gather` at all; it does — `cytosol.internal_face_flux_um_per_s` calls it, and that is
  how the new counter control found `gather`'s uncounted index upload.**
- **Full reductions differ by ~1e-14** between substrates (`sum` 2.3e-14, `dot` 7.1e-15 absolute)
  and this is **not** fixable without giving up the chunked parallel reduction. It is inherent to
  a different summation order. It did not move the whole-cell trajectory at all here — the five-step
  energies are equal to the bit — and §2's second control is why that is a measurement rather than
  luck.
- **`float32`.** Every figure here is float64. The seam at float32 is a precision question and has
  not been measured; the class default remains float32 and a caller taking it gets an unmeasured
  configuration.

## 9. Numerical envelope

Host `float64` throughout. The parity claim is exact (`==` / `np.array_equal`), which is available
precisely because the elementwise kernels and the axis-wise reductions were each measured
bit-identical to NumPy before the whole-cell claim was made. Full reductions are the one place
exactness is *not* claimed, and §8 says so rather than leaving it inside a passing test.

## 10. Source identity, why not clean-room, discarded prose

No source repository, path, symbol or commit. `/Users/sw1/ffn_cellsim` was not opened for this
entry. There is no foreign prose to discard.

## 11. Acceptance, reviewer, rollback

Status `PROPOSED`. Reviewer: the PI. Rollback is the default value in both scenario builders: pass
no `backend` and the cell is built on `NumpyBackend` exactly as before. The `WarpBackend` changes
are additive and no `NumpyBackend` path is touched, so the reference implementation cannot have
moved.

## 12. Open

**Which owner goes device-resident first, and what is the acceptance for it?** The earlier survey
named `cortex` as the best candidate on size. This entry deliberately does not decide it: making an
owner device-resident changes where its state *lives*, which is a different and much larger change
than which substrate its arithmetic runs on, and the parity control in §4 is the thing that will
have to be rewritten when it happens — its non-zero-crossings assertion is designed to fail loudly
on that day rather than quietly keep passing.


## 13. The seam on an actual RTX 4090 — the measurement §6 called not established

Job 31, `sm_89`, float64, `subdivision_level=2`, 5 steps, inside a Slurm allocation.

| run | wall | energy final [pN·µm] | uploads | downloads |
|---|---:|---:|---:|---:|
| `numpy` | 0.574 s | 34714.790796179 | 0 | 0 |
| `numpy` again | 0.556 s | 34714.790796179 | 0 | 0 |
| `warp` CUDA, first | **6.448 s** | 34714.790796179 | 1272 (4.983 MB) | 560 (2.491 MB) |
| `warp` CUDA, warm | **0.860 s** | 34714.790796179 | 1272 (4.983 MB) | 560 (2.491 MB) |

**All four energies agree to the last digit.** `bit_identical_to_numpy: true`,
`relative_energy_gap: 0.0`. That is the claim worth having; the wall clock is the price of it.

**1.55× slower than NumPy, warm.** At `subdivision_level=2` the cell is 794 nodes, the device is
starved, and 1,832 crossings totalling 7.47 MB per 5 steps dominate. **This is the number that
decides whether device residency is worth doing**, and it now exists for the target hardware
rather than being extrapolated from a Mac.

**The first CUDA run took 6.448 s and the warm one 0.860 s — a 5.6 s one-time compilation.** That
independently confirms why §5's original 4.0× reading was withdrawn: it was measuring exactly this.

### Two defects the card found and the Mac could not

Both were live at the commit where §4's fourteen controls passed and this entry claimed the whole
cell runs on the Warp seam. **That claim was true of CPU-Warp and false of CUDA**, and no amount of
CPU testing could have separated them, because Warp's CPU device accepts host arrays where a
`wp.array` is expected — same address space, so leniency costs it nothing.

| # | where | what CUDA said |
|---|---|---|
| 1 | `WarpBackend.scatter_add` | `argument 'dest' expects an array of type array(ndim=2, dtype=float64), but passed value has type ndarray` — `dest` is written and was never routed through `_resolve_out`/`_deliver`. |
| 2 | seven owners' `SOLVE` phase | `RuntimeError: Invalid indexing in slice` — `np.asarray(step)` on a value returned by `backend.scale`. A view on host memory; on device memory it raises. |

(2) was enumerated rather than chased: `cortex_filaments`, `filopodium`, `intermediate_filament`,
`lamellipodium`, `microtubule`, `sf_arc`, `ecm_network`. `cytosol.py` already used
`backend.to_host` correctly at two sites, so the right pattern was in the tree and seven owners did
not use it. The NumPy path is unchanged, measured by stashing the six files around the same run:
`34239.72643872188` both ways, because `NumpyBackend.to_host` *is* `np.asarray`.

**The pattern is worth naming once.** On this project, *"runs on Warp"* measured only on CPU-Warp is
not evidence that it runs on CUDA. The permissive substrate was mistaken for a correct program —
the same shape as §2's `scale` defect, where a leniency that removed an exception also removed what
the exception protected.

### Still not established

- **Whether the device ever wins.** Every figure here is `subdivision_level=2`. The crossings scale
  with the number of operations and the arithmetic scales with the number of nodes, so a larger
  cell should close the gap — but that is an argument, not a measurement, and no run has been made.
- **float32 on CUDA.** Untouched. §8 still applies.

---

## 14. The same ending through the dtype channel — session `e33e377d` (S-BACKEND)

| Field | Value |
|---|---|
| Lane | `e33e377d` S-BACKEND, appended to this entry on PI instruction |
| Scope | `aleph/runtime/backend.py` (`WarpBackend` only), append-only controls in the existing `tests/runtime/test_warp_backend.py`, and this section |
| PI citation | transcript `~/.claude/projects/-Users-sw1-Project-Aleph/e33e377d-49a6-44ed-8cc0-2b16ce5ebe6a.jsonl`, record uuid `5d68087c-76ae-4ac0-9e16-a6e39f8d8033`, `2026-08-05T09:17:04.752Z` (18:17 KST). A citation, never a grant — `ALEPH-PORT-3618`. |
| Substrate | `WarpBackend(device="cpu")` throughout. **No GPU, no Slurm allocation, no card.** §14.8 says what that leaves unestablished. |

§2 of this entry recorded a write that was uploaded, computed, and dropped: the caller's array kept
its old values and nothing raised. §14 is the same ending reached through a different channel. Not
residency this time — **dtype**.

### 14.1 Two defects, and the split between the silent one and the loud one

**(a) `_upload` took its element type from `self.dtype`, not from the array being uploaded.**
`self.dtype` is the backend's *construction* default and is float32. So every host float64 operand
was narrowed on its way to the device, by an argument nobody passed.

| witness | float64 | as uploaded | relative loss |
|---|---|---|---|
| one element | `1.0000000009313226` | `1.0` | `9.3132e-10` |
| `dot(y, y)`, `y` four copies of it | `4.0000000074505806` (NumPy) | `4.0` (Warp) | `1.8626e-09` |

**(b) `scale` chose its kernel with `is32 = a.dtype is wp.float32`.** That identity can only be true
of an array that is **already on the device**. A host array's `.dtype` is `np.dtype('float32')` —
equal in meaning to `wp.float32`, not the same object — so the branch was unreachable from the host
side and `scale_f64` was selected for every host operand, on every backend, always.

**The two defects met in `_resolve_out`, and that is what made the worst symptom.** `_resolve_out`
took the element type from the operand and was *right*; `_upload` took it from `self.dtype` and was
*wrong*. One operation consulted both. So a float64 target was allocated and float32 inputs were
handed to the kernel:

```
add(a, b, out=<host float64>)   on WarpBackend(device="cpu")   # the class default, float32
  -> RuntimeError: array_store() value argument type (float32) must be of the same type
                   as the array (float64)
  -> the caller's `out` still holds array([0., 0., 0., 0.])
```

**Where the two call sites meet you get a crash; where only `_upload` is involved you get a wrong
number and silence.** The reduction channel (`sum`, `dot`, `norm`, `max_abs`) reaches the device
through `_flat` alone and its accumulator is float64 whatever the input, so there is no second call
site to disagree with the narrowing. That is the `dot` row of the table above: no exception, right
shape, right sign, wrong in the ninth digit.

**On the cache-state dependence.** The brief this lane was given reported that with a warm warp
kernel cache the same call returns the caller's array unchanged and raises nothing. **This lane
could not reproduce the silent variant.** Measured on warp 1.15.0, macOS, `device="cpu"`, with a
`WARP_CACHE_PATH` pointed at an empty directory and then at the same directory a second time, the
call raised in both, and a third arrangement — a successful `k_add` instantiation compiled first,
then the bad call — raised as well. What *did* vary with process state is **which** error: a first
bad call raises the codegen error above, a later one in the same interpreter raises
`Failed to find forward kernel '_define__locals__k_add'` after warp has printed
`Failed to lookup symbol`. The caller's `out` is untouched in every arrangement. So the recorded
claim is corrected to: **the message depends on process and cache state, the dropped write does
not.** `test_the_failure_mode_does_not_depend_on_the_warp_kernel_cache` pins the cold-and-warm pair
out of process so the question stays answered rather than remembered.

### 14.2 Why the fourteen parity controls passed straight through it

`tests/scenarios/test_whole_cell_backend_parity.py:39` constructs
`WarpBackend(device="cpu", dtype=np.float64)`, and the file says why: *"float64 rather than the
class default float32 because the claim below is bit-identical."* That reasoning is right and the
coverage it produces is one dtype. On a float64 backend:

- coercing an operand to `self.dtype` is the identity, so (a) cannot fire; and
- `scale_f64` — the branch the broken identity always picks — is the *correct* branch, so (b) is
  right by being wrong.

The `k1` ULP sweep cannot see it either, for the reason §13 already gives about the seam: it grades
kernels that have already been handed device arrays, and every kernel here was correct.

**Two controls, both green, both blind to the same defect, for two different reasons.** The new
controls in `tests/runtime/test_warp_backend.py` use the class's **default** dtype deliberately.

### 14.3 The change

One question — *what element type does this operation compute in?* — answered once, in one place,
and fed to both the allocation and the uploads.

- **`_declared_dtype(operand)`** — the float32/float64 element type an operand declares for itself,
  or `None`. A list, a scalar, or an integer array declares nothing and still falls back to
  `self.dtype`, so nothing that worked before behaves differently.
- **`_operand_dtype(*operands)`** — NumPy's promotion, restricted to the two element types that have
  kernels; `self.dtype` only when no operand declares anything. The rule it encodes is the one this
  class's own docstring already implied: **a construction asks the backend, an operation asks its
  operands.** `NumpyBackend.dtype` is consulted by `zeros` and `array` and by nothing else.
- **`_as_device_array(data, dtype=None)`** — `None` now means *the dtype this array already has*.
  An explicit dtype still wins, which is what `scatter_add` passes when it uploads `values` at
  `dest.dtype`.
- **`_elementwise` / `scale` / `dot`** — resolve once, allocate the target with it, upload both
  operands with it. `scale`'s kernel is selected by `dt == np.dtype(np.float32)`: a value
  comparison, on an object that answers it.

`array()` and `zeros()` are **not** changed. They are constructions, they have no operand to ask,
and `NumpyBackend` coerces in exactly the same place — which is what makes the two backends agree.
`test_an_explicit_construction_still_takes_the_backends_dtype` pins that, because the tempting
"simplification" of this change is to make `array()` preserve too, and that would break parity with
the reference in the other direction.

### 14.4 A third defect, of §13's family, found while reading

`gather` passed its **source** array to `wp.launch` raw. Warp's CPU device accepts a host `ndarray`
where a `wp.array` is expected; CUDA refuses it. That is §13's row 1 — `scatter_add`'s `dest`,
fixed one commit earlier — recurring in the dual operation, and **no CPU test could have seen it**.
It is now routed through `_on_device`, and the observable is the upload counter:
`test_gather_uploads_a_host_source_instead_of_handing_it_to_the_launch` fails on the pre-fix code
because the counter does not move.

Found by reading rather than by running, because the CUDA acceptance run this lane is asking for
(§14.8) would otherwise have spent itself finding it.

### 14.5 A rule conflict, reported and NOT settled — this is a PI question

After (a) and (b), `_elementwise` can be handed a **host float64** operand and a **device float32**
one. Warp's generic kernels require both `Any` parameters to resolve to one element type
(*"Input types must be the same, got ['float64', 'float32']"*). A host operand can be converted for
free — it is being copied across the bus anyway. A device operand cannot, without either a cast
kernel this backend does not have or a round trip through the host, which would charge the operation
the exact crossing this class exists to count.

**Two rules in this file disagree about what to do, and `CLAUDE.md` §2.5 says that is a finding.**

| rule | source | what it says here |
|---|---|---|
| `NumpyBackend` is the definition of every operation | `WarpBackend` docstring | promote `float64 + float32` to float64, as `np.add` does |
| float32 is the *compute* channel | `ALEPH-DQ-107`, the precision table in the same docstring | promoting every mixed pair to float64 makes that channel vacuous, because every owner's host state is float64 |

The two candidate resolutions, with their costs:

- **(B) promote like NumPy.** Needs a device-side cast. `wp.utils.array_cast` exists in warp 1.15.0
  and was measured working for 1D in both directions and **failing for 2D**
  (`Couldn't find function overload for 'float64' that matched inputs with types:
  [array(ndim=1, dtype=float32)]`), so it would have to be applied to the flattened view or replaced
  by two explicit kernels in `warp_kernels.py` in the style of `k_scale_f32`/`k_scale_f64`. Cost:
  one appended kernel pair and one extra launch in the mixed case. Consequence: on a float32
  backend almost all arithmetic is actually performed in float64, because the owners hold float64.
- **(C) let the device operand's dtype decide; the host side converts on upload.** No new kernel, no
  new launch, and `ALEPH-DQ-107` means what it says. Consequence: the precision of `add(x, y)`
  depends on which operand happened to be uploaded earlier — an answer that depends on residency,
  which is the class of thing this project's controls exist to prevent.

**Neither is taken here.** The seam now **refuses by name**, in the same shape as the existing
device-broadcast refusal, with a message that states the likely cause and points at this section.
It is inert in every configuration the tree currently runs: at `dtype=np.float64` every device array
a whole-cell step constructs is float64 and no pair reaches it — `tests/runtime` and
`tests/scenarios` are green with it in place (§14.9).

### 14.6 The float32 path: what the fix bought, and where it still stops

`build_whole_cell(subdivision_level=1, backend=WarpBackend(device="cpu", dtype=np.float32))`

| | before | after |
|---|---|---|
| construction | builds | builds |
| first step | `RuntimeError: Error launching kernel 'k_scale_f64', argument 'a' expects an array with dtype=float64 but passed array has dtype=float32` at `aleph/vertical/cortex_filaments.py:1751` | that line passes |
| next | — | stops at **`cortex_filaments.py:1752`** with §14.5's named refusal |

So defect (b) was the blocker the brief named, and it is gone. **The remaining blocker is not in the
seam.** Line 1752 is
`backend.add(self.forces, backend.array(self.internal_forces()), out=self.forces)`: the owner asks
the backend to *construct* an array — which correctly takes the backend's dtype, float32 — and adds
it to `self.forces`, which the owner keeps as host **float64** regardless of the backend. The owner
is mixing "let the backend pick the element type" with "my state is float64", and the seam is where
the two meet.

**This is a pattern, not a line.** `backend.array(` appears at **61 call sites across 15 modules**
under `aleph/vertical/**` — `cortex_filaments`, `cortex`, `membrane`, `nucleus`, `cytosol`,
`microtubule`, `intermediate_filament`, `lamellipodium`, `filopodium`, `focal_adhesion`, `sf_arc`,
`ecm_network`, `nmii`, `medium`, `pressure`. All of `aleph/vertical/**` belongs to other lanes, so
this is **reported and not patched**, per `CLAUDE.md` §1. Whichever way §14.5 is decided, those
owners are where the float32 path is finished.

**Unchanged by all of this, and measured both ways:** `subdivision_level=1`, seed 7, 5 steps —
`numpy` and `WarpBackend(device="cpu", dtype=np.float64)` both give `energy_final_pn_um =
29488.36676227611`, before the change and after it. The float64 seam did not move.

### 14.7 Gate G11 — scoped only, deliberately not executed

G11 asks that the float32 path be *graded*. Scoping it, and nothing more:

1. **A ULP parity harness over the float32 path**, scoring the nineteen ledgers `-3601`…`-3619`.
   Re-baselining those is a **separate decision** and is not proposed here.
2. **It is a physics decision before it is a backend one.** `ALEPH-PORT-3603` §13 records that
   float32 descent converges in **0 of 12 settings** and worsens recovery sigma by **331×**. A
   graded-but-divergent path is not a path.
3. **§14.5 has to be decided first**, because it changes what "the float32 path" computes.

**Nothing in G11 was run.** This section is a scope, and it is the PI's to schedule.

### 14.8 UNVERIFIED — every number in §14 is CPU-Warp

`§13` is this file's own record that CPU-Warp passing is not evidence about CUDA: `scatter_add`'s
`dest` and seven owners' `np.asarray(step)` were accepted by CPU-Warp and refused by the card. §14.4
is a third instance of exactly that, found by reading rather than by a device.

So: **the fix in §14.3 is UNVERIFIED on CUDA.** What a device run would settle, and nothing else
can:

- that the resolved dtype reaches a CUDA kernel as the same element type it resolved to;
- that `gather`'s now-uploaded source is accepted where the raw host array was refused;
- that `WarpBackend(device="cpu", dtype=np.float64)` and the same on CUDA still agree to the bit on
  the five-step whole cell.

The request, with the numbers a PI decision needs, is in §14.10. **No job has been submitted and
none will be before the PI states a card and a duration in session** — `ALEPH-PORT-3618`: the record
is a citation, never a grant.

### 14.9 Acceptance

Nine new controls in `tests/runtime/test_warp_backend.py`, appended to the existing file because
`e874a7fe` (Lane G10) claims **new** files under `tests/runtime/**`. Each was run against the
pre-fix `backend.py` in an isolated copy of the tree before being run against the fixed one:

| control | pre-fix |
|---|---|
| `test_an_uploaded_operand_keeps_its_own_dtype_rather_than_the_backends_default` | FAIL — `4.0` vs `4.0000000074505806` |
| `test_a_host_out_is_written_even_when_the_backend_defaults_to_a_narrower_dtype` | FAIL — codegen error, `out` left at zeros |
| `test_scale_selects_its_kernel_by_value_and_not_by_object_identity` ×4 | FAIL ×3; **the `float64`/`float64` cell PASSES pre-fix** — that cell is the whole-cell parity suite's configuration, and it passing for the wrong reason is §14.2 in one line |
| `test_the_float32_backend_survives_the_call_that_kills_the_whole_cell` | FAIL |
| `test_a_mixed_dtype_device_operand_is_refused_by_name_rather_than_by_warp` | FAIL |
| `test_gather_uploads_a_host_source_instead_of_handing_it_to_the_launch` | FAIL — the upload counter does not move |
| `test_the_failure_mode_does_not_depend_on_the_warp_kernel_cache` | FAIL, out of process, cold cache and warm |
| `test_an_explicit_construction_still_takes_the_backends_dtype` | **PASS** — it pins behaviour the change preserves |

Suites, after: see §14.11.

### 14.10 What a CUDA acceptance run needs, stated so the PI can decide it

Requested, **not** claimed:

| field | value |
|---|---|
| card | one of `4090-0`, `4090-1`, `3090` — any; nothing here is memory-bound |
| wall clock | `00:20:00` |
| memory | `MEM_GB=16`, `CPUS_PER_TASK=4` |
| command | `tests/runtime/test_warp_backend.py` and `tests/scenarios/test_whole_cell_backend_parity.py`, plus the five-step whole cell at `subdivision_level=1` on CUDA float64 compared to `29488.36676227611` |
| submitted through | `MEM_GB=16 CPUS_PER_TASK=4 gpu-submit <card> 00:20:00 <command>` |

**One run should be enough**, which is why §14.4 was fixed by reading first rather than left for the
device to find.

### 14.11 Suites

| suite | result | tree state |
|---|---|---|
| `tests/runtime/test_warp_backend.py` + `test_backend.py` + `test_backend_krylov.py` | **72 passed, 3 skipped.** 11 of the 72 are new here; the 3 skips are the GPU tests, with their absent hardware named in the reasons. | **Isolated.** These three import `aleph.runtime` only, so no other lane's in-flight file can reach them. |
| `tests/scenarios/test_whole_cell_backend_parity.py` | **14 passed** — the controls this change had to not break. | Imports `aleph.vertical.**`; re-run at 18:5x KST against the tree as it then stood. |
| `tests/runtime` + `tests/scenarios` | **2 failed**, both foreign and both proven so — §14.11a. **No pass count is quoted from that run**, deliberately: it finished at 18:47:17 KST and the live `926ddf90` (S-CORTEX) row wrote `aleph/vertical/assembly.py` at 18:41:06 and `aleph/vertical/cortex_surface_coupling.py` at 18:42:40, so it overlapped another lane's edits. `CLAUDE.md` §1 says a whole-repo run can fail on a file that is mid-write; a count taken across one is not this lane's to report. The two failures stand anyway, because each was reproduced separately against an isolated copy at `HEAD`. | mixed |

### 14.11a Two red tests, both foreign, both proven so rather than asserted

`CLAUDE.md` §1 says to report a foreign breakage rather than fix it. Both were re-run against an
**isolated copy of the tree carrying `backend.py` at `HEAD`** — every other file identical — and
both fail there too, so neither is this lane's:

| test | owner | proof |
|---|---|---|
| `tests/runtime/test_reference_input_guard.py::test_every_registered_case_reference_is_sensitive_to_position_precision` | live `e874a7fe` (Lane G10) — its subject is `law_cases.py`/`law_kernels.py`, which that row owns and is editing | fails identically with pre-fix `backend.py`. Already recorded as foreign by Lane W2 at the bottom of `docs/ACTIVE_SESSIONS.md`, and by `HANDOFF.md` §C-0 before that. |
| `tests/scenarios/test_whole_cell.py::TestTheExteriorMediumResistsTheRigidBodyModes::test_the_exterior_medium_actually_dissipates_work_during_a_run` | live `46143f30` (Lane W2) — `aleph/scenarios/**`, and `aleph/vertical/medium.py` | fails identically with pre-fix `backend.py`. It cannot be this lane's by construction either: it calls `build_whole_cell(seed=7)` with the default backend, which is `NumpyBackend`, and **every hunk of this change is inside `class WarpBackend`** — `NumpyBackend` is byte-identical to `HEAD`. |

### 14.12 The parity harness has the same blind spot, and this is G11's problem

`aleph/runtime/parity.py` belongs to this lane and **was not changed**. It did not need to be for
the fix — and reading it produced a finding that §14.7 needs.

Every case in `default_cases()` wraps its inputs in `be.array(...)`:

```python
cases.append(OpCase("add", (n,), (a, b),
                    lambda be, xs: be.add(be.array(xs[0]), be.array(xs[1])), _amp_elementwise))
```

`be.array(...)` is the **explicit construction** path, so every case hands the operation a device
array already at the backend's own dtype. **The parity harness therefore cannot see the seam** — for
exactly the reason §13 gives about the `k1` ULP sweep, which was thought to be a property of that
sweep and is in fact a property of the whole grading apparatus. Neither instrument has a case in
which a host array reaches an operation, which is the only way the whole cell actually calls this
class.

**Not fixed here, deliberately.** `tests/runtime/test_warp_backend.py::test_gpu_parity_is_within_float32_roundoff`
asserts `parity_report(...).all_within_roundoff` **on CUDA**, and adding cases to the set it grades
would put untested-on-device cases into the one acceptance run §14.10 asks for. Host-operand cases
belong in G11's scope, alongside the ULP harness, and they are named here so that scope is written
from a measurement rather than from an assumption.

### 14.13 Who this section hands work to

Three lanes, and none of the work is done here:

1. **`926ddf90` (S-CORTEX)** owns `aleph/vertical/cortex_filaments.py` and registered while this
   lane was running. §14.6 is theirs: line 1752's `backend.add(self.forces,
   backend.array(self.internal_forces()), out=self.forces)` mixes a backend-typed construction with
   host float64 state, and the seam now refuses it by name instead of failing in warp's codegen.
   The same shape appears at **61 call sites across 15 modules**, so it is a pattern to decide once
   and not a line to patch. `cortex_filaments.py` was last written at 14:11 KST, before that row
   opened, so the line numbers here are against the committed file and will move.
2. **The PI** owns §14.5 — promote like NumPy, or let the device operand's dtype decide — and
   §14.7's G11 schedule, and the §14.10 card and duration. None of the three is taken here.
3. **Whoever runs the acceptance job** gets §14.8's list. §14.4 was fixed by reading rather than by
   the device precisely so that one run is enough.

### 14.14 `tests/ports` — five red, none of them this append's

Checked rather than assumed, because this section is an append to a ledger entry and the port guards
read ledger entries.

| guard | why it is red | this lane's part |
|---|---|---|
| `test_ledger_entries_are_complete` | `-3632` and **`-3634`** are missing mandatory header fields (`source_path`, `independent_oracle`, `negative_control`, `positive_control`) | **none.** The `\| Field \| Value \|` header table this guard parses is **byte-identical between `HEAD` and this working copy** — verified by comparing the two directly. §14 appends sections below it and does not touch it. Completing another lane's header would be writing their record. |
| `test_index_is_not_stale`, `test_the_index_is_not_missing_a_tracked_entry` | `ports/ledger/INDEX.md` does not list `-3634`, which is **committed at `HEAD`** | **none, and deliberately not repaired.** `regenerate_index.py` reads the working tree, which right now carries the live `926ddf90` row's uncommitted `aleph/vertical/**` edits. Every row in `docs/ACTIVE_SESSIONS.md` since 2026-07-31 records this file as red and not regenerated, for the same reason. |
| `test_named_controls_resolve_to_real_tests` | fails on a pristine `HEAD` checkout | pre-existing |
| `test_no_provider_vocabulary_leaks_into_the_package` | `aleph/harness/ingest/adapters.py`, the open `7bc39358` ingest row | foreign, already reported twice in `docs/ACTIVE_SESSIONS.md` |

**No guard was edited.** `CLAUDE.md` §2.4: a guard that blocks you is a guard working.

---

## 15. Three PI answers, transcribed — session `e33e377d`

> `d/넣자/cuda 수용 실행`
>
> — PI, `2026-08-05T10:06:05.492Z` (19:06 KST)
> `~/.claude/projects/-Users-sw1-Project-Aleph/e33e377d-49a6-44ed-8cc0-2b16ce5ebe6a.jsonl`
> record uuid `433e123e-beca-4814-ab41-7d8fe81a6f83`

Answering, in order, the three questions §14.5, §14.7 and §14.10 put. **A citation, never a grant**
(`ALEPH-PORT-3618`) — it records what was said and is refuted by opening the transcript. The
preceding message in that transcript is the one carrying the option tables, so what "d" and "넣자"
select is checkable rather than asserted. **No `decided_by` field appears anywhere in this entry.**

### 15.1 `d` — §14.5's mixed-dtype pair is fixed in the owners, not in the seam

Option (D): the owners hold their state at a dtype consistent with the backend they were handed, and
then the mixed pair never forms. **The seam's refusal stays exactly as §14.5 left it** — it is not a
stopgap, it is the seam correctly reporting that a caller asked for two incompatible things in one
line.

Neither (B) nor (C) is taken, and the reason each was rejected is worth keeping:

- **(B) promote like NumPy** would have made `dtype=np.float32` compute in float64 for essentially
  every operation, because 25 of 25 owner state arrays are float64 — measured, §14.6. The mode would
  have been a label with nothing behind it.
- **(C) let the device operand decide** would have made the precision of `add(x, y)` depend on which
  operand happened to be uploaded first. An answer that depends on residency history rather than on
  the arithmetic is the class of thing this project's controls exist to prevent.

**Nothing further is done here.** The work is 61 call sites across 15 modules in `aleph/vertical/**`
and belongs to the lanes that own them — for `cortex_filaments.py`, the live `926ddf90` (S-CORTEX)
row.

### 15.2 `넣자` — §14.7's host-operand cases go into the harness, and they needed a second fix

Landed as **`aleph/runtime/parity.seam_cases()`**. Nine rows, every one handing an operation a
**host** array: `add[host]`, `multiply[host]`, `scale[host]`, `dot[host]`, `sum[host]`,
`norm_axis1[host]`, `add[host,out]`, `scatter_add[host]`, `gather[host]`.

**Adding the cases was not enough, and that is the finding.** `compare_case` computed
`bound = U32 * (amp + 1)` for every case unconditionally. The narrowing defect's relative error is
`9.3132e-10`, which is **below `U32 = 5.960e-08`** — so a host case added under the old bound would
have been graded *within round-off* and the report would have said all-clear over a broken seam.
Measured on the pre-fix backend:

| row | measured | at `U32` | at `U64` |
|---|---|---|---|
| `dot[host]` | `8.1441e-09` | bound `3.9378e-06` → **ok** | bound `7.3347e-15` → **OVER BOUND** |
| `sum[host]` | `4.1614e-08` | bound `5.4104e-06` → **ok** | bound `1.0078e-14` → **OVER BOUND** |

So `OpCase` gained `unit_roundoff`, defaulting to `U32` — every existing case is graded exactly as
before — and the seam rows declare `U64`.

**A guard stopped the first attempt, and it was right.** The host cases were first written straight
into `default_cases()`, and `tests/runtime/test_parity.py::test_a_perturbation_below_the_bound_is_not_flagged`
went red: it asserts a `1e-9` relative perturbation is *inside* the bound, which is true of a float32
set and false of a float64 one. `test_only_the_perturbed_op_is_flagged` went red too, on the
one-row-per-op assumption. **Neither test was edited** (`CLAUDE.md` §2.4). The set was split instead,
which is the better design anyway: `default_cases()` asks *how far does float32 drift*, `seam_cases()`
asks *does the operand arrive unchanged*, and one set cannot answer two questions.

That split is itself pinned by a paired control,
`test_the_seam_bound_is_the_thing_that_makes_the_set_able_to_fail`: one `1e-9` perturbation on `add`,
**not** flagged by `default_cases()` and **flagged** by `seam_cases()`. Same number, two correct
verdicts.

`tests/runtime/test_parity.py` + `test_warp_backend.py`: **69 passed, 4 skipped** (the skips are the
GPU tests, including the new `test_gpu_seam_parity_is_within_float64_roundoff`).

### 15.3 `cuda 수용 실행` — §14.10's acceptance run

The PI directed the run. **The card and the wall clock were not restated independently**, so they are
taken as §14.10 proposed them and **not extended**: `00:20:00`, `MEM_GB=16 CPUS_PER_TASK=4`, and one
card from the three §14.10 offered. That reading is recorded here rather than assumed, so a session
that disagrees can say so against the same transcript record. Outcome in §16.

---

## 16. The acceptance run, on an RTX 4090 — §14.8's `UNVERIFIED` closed, except for one half

| Field | Value |
|---|---|
| Job | Slurm `3`, `MEM_GB=16 CPUS_PER_TASK=4 gpu-submit 4090-0 00:20:00` |
| Device | **NVIDIA GeForce RTX 4090**, `cuda:0`, **sm_89**, 128 SMs, 24.0 GiB |
| Host | `DESKTOP-EEUU4RQ`, WSL2 `5.15.167.4-microsoft-standard`, python 3.12.3 |
| Substrate | **warp 1.14.0**, numpy 2.5.1 — note the warp version differs from the 1.15.0 the CPU figures were taken on, and the fix holds on both |
| Wall clock | `2026-08-05T10:23:59Z` → `10:24:23Z`, **24 s** of the 20 min asked for |
| Citation | `86f712d6-5a83-47bc-b8c7-612e15a1157e`, well inside the 2 h stated. Not extended. |
| Artefact | `runs/cuda_accept_3634_seam.json` on the GPU host; `scripts/cuda_accept_3634_seam.py` produced it |

### 16.1 The six defect checks, on the card

| check | ledger | result |
|---|---|---|
| `upload_preserves_host_float64` | §14.1(a) | **PASS** — `dot(y,y)` = `4.000000007450581` on the device, equal to NumPy, `rel = 0.0` |
| `host_out_is_written` | §14.1, §2 | **PASS** |
| `scale_selects_by_value[float32]` | §14.1(b) | **PASS** |
| `scale_selects_by_value[float64]` | §14.1(b) | **PASS** |
| `gather_uploads_a_host_source` | §14.4 | **PASS** — *and this is the one the Mac could not have answered* |
| `scatter_add_writes_a_host_dest` | §13 row 1 | **PASS** |

The last two are the point of having run this at all. §13 records that CPU-Warp accepted a host
`ndarray` where CUDA refused one, so those two rows were `UNVERIFIED` in the strict sense — not
"probably fine", but *untested on the only substrate that can fail them*. They now pass on sm_89.

### 16.2 Both parity sets, on the card

`default_cases()` — **within bound**, so the float32 kernel grading did not regress.

`seam_cases()` — **within bound**, nine rows, on `cuda:0` at `compute_dtype=float32`:

```
op                        shape         n     max_abs     max_rel       amp       bound
add[host]                 10000     10000   0.000e+00   0.000e+00  1.04e+00   2.270e-16  ok
multiply[host]            10000     10000   0.000e+00   0.000e+00  1.00e+00   2.220e-16  ok
scale[host]               10000     10000   0.000e+00   0.000e+00  1.00e+00   2.220e-16  ok
dot[host]                 10000     10000   2.842e-14   1.635e-15  3.68e+02   4.094e-14  ok
sum[host]                 10000     10000   5.684e-14   7.287e-16  1.02e+02   1.144e-14  ok
norm_axis1[host]         4096x3     12288   8.882e-16   1.938e-16  1.00e+00   2.220e-16  ok
add[host,out]             10000     10000   0.000e+00   0.000e+00  1.04e+00   2.270e-16  ok
scatter_add[host]        4096x3     12288   0.000e+00   0.000e+00  1.64e+00   2.934e-16  ok
gather[host]             8192x3     24576   0.000e+00   0.000e+00  1.00e+00   2.220e-16  ok
```

**Six of the nine are exactly `0.000e+00`** — a float32 backend handed host float64 operands returns
bit-identical answers to NumPy on a GPU, which is the whole claim of §14.3 stated as a measurement.
The three that are not zero are the reductions and `norm_axis1`, where the chunked decomposition
sums in a different order than NumPy's pairwise reduction; each is inside its own float64 bound.

### 16.3 The half that did NOT run, and it is not a defect in this change

`whole_cell` is **absent from the artefact**. The job errored building it:

```
File "aleph/scenarios/whole_cell.py", line 1065, in _nmii
File "aleph/vertical/nmii.py", line 211, in <module>
    from scipy.integrate import quad
ModuleNotFoundError: No module named 'scipy'
```

The GPU host's warp environment (`~/wp`) carries numpy and warp-lang and **not** scipy. So:

- **`§14.6`'s figure `29488.36676227611` remains UNVERIFIED on CUDA.** The claim that the composed
  cell produces a bit-identical trajectory on the device is not established by this run and is not
  claimed by it.
- **It is a finding for other lanes, not only for this one.** `aleph/vertical/nmii.py` imports
  `scipy.integrate.quad` at **module scope**, so *any* lane that wants to build the whole cell on
  the GPU host hits this before it reaches a kernel. That file belongs to the closed `8448bfd4`
  (S-MOTOR) row and is **reported, not patched**.
- **No package was installed to route around it.** Installing scipy into a shared machine's
  environment is a change nobody asked for, and this lane's own commit message for
  `scripts/cuda_accept_3634_seam.py` says so about pytest. Put to the PI instead.

### 16.4 What is established now, and what is not

**Established on hardware:** the dtype channel — `_operand_dtype`, the value-comparison kernel
selection, `gather`'s uploaded source, `scatter_add`'s host destination, and all nine `seam_cases()`
rows — on an RTX 4090, sm_89, warp 1.14.0.

**Still `UNVERIFIED`:**

1. The whole cell on CUDA (§16.3), blocked on scipy in `~/wp`.
2. **Whether the device ever wins.** §13 already said this and the run does not change it: 24 s of
   wall clock says nothing about throughput, no timing comparison was made, and the seam still
   uploads host arrays per operation. Nothing here is a performance claim.
3. **float32 as a physics configuration.** §14.7 stands untouched: 0 of 12 settings converge and
   recovery sigma worsens 331× (`ALEPH-PORT-3603` §13). This run grades the *seam* at float32; it
   says nothing about whether float32 should be switched on.

---

## 17. §16.3's gap closed, and the acceptance run passes whole — job 5

The PI directed installing scipy into the GPU host's warp environment. Done and recorded:
`scipy 1.18.0` into `~/wp`, **`numpy` unchanged at `2.5.1`** (pip reported the requirement already
satisfied), so the §16 figures and these are comparable. Reversible with
`~/wp/bin/python3 -m pip uninstall scipy`. Before: `numpy==2.5.1`, `pip==26.1.2`, `warp-lang==1.14.0`
— three packages.

Verified first that scipy was the *only* gap: a five-step whole cell loads exactly one third-party
module beyond numpy and warp. `duckdb` and `ffn_sim` are lazy imports and this path does not reach
them; the module-scope scipy imports are `aleph/vertical/nmii.py:237-238` and
`aleph/state/locality.py:87`.

### 17.1 Job 4 failed, and the defect was in this lane's own harness

`whole_cell:not_bit_identical` — CUDA float64 `29522.828995645967` against the hardcoded
`29488.36676227611`. **The device was right.** Between the constant being measured at `698bab5` and
the job running at `0713140`, `42a360e` gave the nmii motor a clock (`bound_fraction`
`0.000000 → 0.375000`) and the energy moved. Re-measured in an isolated checkout at the run's own
commit:

| substrate | five-step energy at `0713140` |
|---|---|
| `NumpyBackend` | `29522.828995645967` |
| `WarpBackend(device="cpu", dtype=np.float64)` | `29522.828995645967` |
| **`WarpBackend(cuda:0, dtype=np.float64)`** | **`29522.828995645967`** |

Three substrates, bit-identical, and one stale constant.

**This is worth more than the pass it delayed.** A hardcoded expectation is a gate that certifies a
*different configuration* from the one that runs, which is precisely what `aleph/runtime/parity.py`'s
module docstring says the predecessor project shipped — and this lane wrote one anyway, in a file
whose whole purpose is acceptance. The failure is also **asymmetric in the wrong direction**: the
device gets blamed for the constant's age, on a substrate this project already has documented reason
to distrust. A reader would have started looking at CUDA.

Fixed by removing the constant, not by updating it: `_whole_cell` now builds both worlds from the
same builder with the same seed **in the same process**, one on `NumpyBackend` and one on the device.
The only difference between them is the substrate, which is the entire claim. The old figure survives
as `HISTORICAL_ENERGY_PN_UM_AT_698BAB5` with the drift reported in the artefact, because how fast
this number moves is itself worth seeing: **59.7 % in about two hours**, across four other lanes'
commits.

### 17.2 Job 5 — the run, whole

| Field | Value |
|---|---|
| Job | Slurm `5`, `4090-0`, `00:20:00` asked, **8 s** used |
| Commit | `9d20f35` |
| Device | RTX 4090, `cuda:0`, sm_89, warp 1.14.0, numpy 2.5.1, scipy 1.18.0 |
| Citation | `86f712d6-…`, inside the 2 h stated, not extended |

**`passed: true`, `failures: []`.**

| | result |
|---|---|
| defect checks | **6 / 6** — including `gather`'s uploaded host source and `scatter_add`'s host dest, the two CPU-Warp cannot answer |
| `default_cases()` | within bound |
| `seam_cases()` | within bound, nine rows |
| **whole cell** | **bit-identical**: `11895.02466168799` on CUDA float64 against the in-process `NumpyBackend` reference, `relative = 0.0` |

`host_uploads = 3552`, `host_downloads = 1154` for five steps — **the same counts CPU-Warp reports**,
which is the expected result while the owners hold host NumPy and is a cost figure, not a speed one.

### 17.3 What §14.8's list looks like now

1. ~~The whole cell on CUDA~~ — **closed here.** Bit-identical at `9d20f35`.
2. **Whether the device ever wins** — still open, and untouched. 8 s of wall clock is not a
   throughput measurement, no timing comparison was made, and the seam still uploads per operation.
   Nothing in §16 or §17 is a performance claim.
3. **float32 as a physics configuration** — still open. §14.7 stands: 0 of 12 settings converge and
   recovery sigma worsens 331× (`ALEPH-PORT-3603` §13). This run grades the *seam* at float32 and
   says nothing about whether float32 should be switched on.

---

## 18. §14.5's named refusal, inventoried — session `26c58bd8` (S-RESIDENCY)

**Append-only. Nothing in `aleph/vertical/**` was edited by this session; every live row that owns
one of those files still owns it.** This section is a work order and a finding. It grades no gate.

### 18.1 What was counted, and why the criterion is narrower than the last one

§14.5 records option **(D)** — the mixed-dtype pair is resolved in the owners, not in the seam — so
`_agreed_dtype` keeps raising by name. The question this session was given is *which file, which
line*.

The refusal fires when an operation receives one operand already on the device at element type `D1`
and resolves `D2 ≠ D1` from the pair. Therefore **the work is the set of allocator calls that mint a
device array without a dtype**, not the set of backend calls. Measured over `aleph/vertical/**`:

| | count |
|---|---:|
| `backend.<method>(...)` calls, all methods | **148**, in 19 modules |
| of those, calls that mint a device array (`array`, `zeros`) | **26** |
| of those, calls that pass **no** dtype and inherit `WarpBackend.dtype` | **23**, in **16 modules** |
| calls that already decide their dtype | 3 (`cytosol.py:1997`, `:2012`, `:2030`) |

The other 122 calls consume arrays minted by those 23. Fix the mint and the consumer follows, which
is why the list is 23 lines rather than 148.

**S-BACKEND's closing figure of "61 call sites in 15 `aleph/vertical/**` modules" is not reproduced
here and this session cannot reconcile it**, because no script was left for it. The criterion above
is stated exactly so that the next reader is not left to guess at this one too. The **16 modules**
agrees; the site count does not, and the difference is the criterion.

### 18.2 The list

Full table, per module and per line, with which sites an L1 whole cell actually reaches:
`docs/results/2026-08-05-float32-owner-inventory/README.md` §2. Regenerate it with

    OPENBLAS_NUM_THREADS=1 <aleph python> docs/results/2026-08-05-float32-owner-inventory/verify.py

Head of the distribution: `cytosol.py` **6**, `focal_adhesion.py` **2**, `sf_arc.py` **2**, and one
site each in `cortex.py`, `cortex_filaments.py`, `ecm_gel.py`, `ecm_network.py`, `filopodium.py`,
`intermediate_filament.py`, `lamellipodium.py`, `medium.py`, `membrane.py`, `microtubule.py`,
`nmii.py`, `nucleus.py`, `pressure.py`. Fifteen of the sixteen modules are literally the same line:

```python
backend.add(self.forces, backend.array(self.internal_forces()), out=self.forces)
```

**Confirmed in situ rather than inferred.** A whole cell on `WarpBackend(device="cpu",
dtype=float32)` raises on its first step; the deepest frame inside `aleph/vertical/**` is
**`cortex_filaments.py:1877`**, the first entry in the list, under `backend.py:773 → :759 → :551`
(`add → _elementwise → _agreed_dtype`). This is the same file S-BACKEND named as where the float32
path dies, now with the line and the stack. `verify.py` asserts it.

**Two of the sixteen modules are being edited by other live sessions while this was measured**;
`microtubule.py`'s site moved from line 1895 to 2031 between two runs twenty minutes apart. The line
numbers are a snapshot and the script is the authority.

### 18.3 The decision this session did not take

Two repairs both remove the refusal and they are not equivalent:

- **(a)** pin each mint to float64 — 23 one-line edits, removes every refusal, and makes the float32
  *compute* channel vacuous at those sites, which is the exact tension `_agreed_dtype`'s docstring
  and §14.5 already name;
- **(b)** give the owners device-resident state at the backend's dtype — the real repair, priced as
  *not worth it now* by `docs/results/2026-08-05-engine-cost-audit/README.md` §5.7 below ~120,000
  cortex nodes.

(a) makes the float32 path **run** without making it **compute in float32**. Reporting that as "the
float32 path works" would be a claim above its evidence, and §14.7 stands — 0 of 12 settings converge
and recovery sigma worsens 331× (`ALEPH-PORT-3603` §13). **Which of (a) and (b) a sweep lane should
do is left to the PI.**

### 18.4 One site option (D) does not cover — reported, and deliberately not fixed

`aleph/runtime/backend.py:1189`, inside `jacobi_preconditioner`:

```python
inverse = backend.array(1.0 / d)              # no dtype → the backend default
def apply(residual):
    return backend.multiply(residual, inverse)  # meets whatever `residual` is
```

Reached from `cytosol.py:2073` during a real step, so it is live. **Same defect shape as the 23, but
it is in the seam** — and option (D) says the seam does not resolve mixed dtypes, so the rule as
recorded does not name who fixes it. `backend.py` is this session's path and this session could have
changed the line; it did not, because choosing the dtype there is the same precision decision as
§18.3 and taking it quietly in a file nobody is reading is the drift option (D) exists to prevent.

---

## 19. The residency invariant — a wrong reading, and which side was fixable

`docs/results/2026-08-05-engine-cost-audit/README.md` §5.5 reports *"the residency invariant is
already violated"* and cites the L3 whole cell at 82 `to_host` calls/step, 110,190 elements. Measured
(`docs/results/2026-08-05-residency-invariant-scope/verify.py`, no arguments), the sentence had no
scope clause and permitted three readings:

| Reading | Measured | Verdict |
|---|---|---|
| `ResidentMembraneState`, which owns `push_positions`/`pull_positions` | L3, 3 steps: 9 crossings, sizes `5,3,3` repeating, **0 mesh-sized** | **holds** |
| the whole module, i.e. including `resident_membrane_forces_and_energy` | L2, 3 steps: **4 mesh-sized crossings, 486 elements each** | **violated, and its own `TransferCounts` already reported it exactly** |
| the L3 whole cell — the audit's subject | **`ResidentMembraneState` constructed 0 times**; `push_positions` never called | **out of scope** |

The `5, 3, 3` are the two-stage reductions' float64 partials (`ceil(1280/CHUNK)`, `ceil(642/CHUNK)`),
which is what `TransferCounts` declares it models — not arrays.

**The declaration was the fixable side and it was fixed.** Both violating readings need edits inside
`aleph/vertical/**` / `aleph/scenarios/**`, which live rows own and which `residency.py`'s own *"What
is NOT here"* section already declares out of scope; the whole-cell one is additionally priced as not
worth it now by the cost audit's own §5.7. The invariant block in `aleph/runtime/residency.py` now
carries its scope, both carve-outs and the measured numbers.

**Where this session and the audit disagree, and it is not resolved.** Re-measuring the audit's own
subject: **14.0 `Backend.to_host` calls/step pulling 47,741 elements**, and 180.7 `_download`
calls/step pulling 120,264 elements one level deeper. The audit reports 82 / 110,190. **The largest
single readback agrees exactly at 8,652 elements**, so both instruments watch the same phenomenon and
the audit's underlying finding — the whole cell pulls mesh-sized arrays back — is confirmed. The call
counts differ because the instruments attach at different depths. The audit left no reproduction
script, so **the gap cannot be closed from this side**; reported, not reconciled, and not counted
against either lane.

**One control was not filed.** The reading-2 carve-out wants a control in
`tests/runtime/test_residency.py` beside the existing C5. That file is not among this session's
paths, so the control is specified in
`docs/results/2026-08-05-residency-invariant-scope/README.md` and left for the owning lane;
`verify.py` asserts the same two facts meanwhile, and it fails in **both** directions so the
docstring's carve-out is flagged if it ever goes stale.
