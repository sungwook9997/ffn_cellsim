# ALEPH-PORT-3628 — the per-site Python call in the step

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3628` |
| Lane | `1afc4d25` — S5, the per-site Python call in the step |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **before the code**, per `PLAN.md` §0.2.5 and `CLAUDE.md` §3 |
| Port class | `NOT A PORT` — nothing is taken from `/Users/sw1/ffn_cellsim`. This entry changes no law and derives no physics. It removes a Python loop from in front of arithmetic that already exists and must come out the other side **bit-identical**. It is filed as a ledger entry because the repository requires one before code touches `aleph/**`, and because the bit-identity argument in §4 is the only thing that makes the change safe. |
| Aleph target | `aleph/vertical/cortex_filaments.py`, `aleph/vertical/ecm_network.py`, `aleph/vertical/sf_arc.py` — **append-only for the new method, plus the in-module callers of the old one rewired**. `tests/vertical/test_material_point_batch.py` (new). |
| Supersedes | Nothing. `resolve_material_point` survives **unchanged, in all three owners**, and becomes the oracle this entry is graded against. |
| PI decision cited | None. The brief is an engineering instruction, not a model choice; no `decided_by` field appears anywhere in this entry. |

---

## 0. The brief re-measured, and one of its two diagnoses withdrawn

The brief supplied two groups of measured facts and required both to be re-measured before any code
was written. **The call counts reproduce exactly. The timing table does not reproduce at all, and
the diagnosis built on it is wrong** — though the conclusion drawn from it happens to be right, for
a different and stronger reason.

### 0.1 The call counts — exact, all three

Reproduced with the brief's own command, `subdivision_level=1`, `max_steps=6`, `sort='tottime'`:

| Brief | Re-measured (this session) | Verdict |
|---|---|---|
| `9600 calls  0.052s  cortex_filaments.py:1066  resolve_material_point` | **9600 calls**, `tottime` 0.041 s, `cumtime` 0.078 s | ✅ count exact |
| `920 calls  0.012s  ecm_network.py:910  resolve_material_point` | **920 calls** — 720 from `crosslink_arrays`, 200 from `endpoint_sink` | ✅ count exact |
| `40 calls  0.119s(cum)  cortex_filaments.py:1226  crosslink_arrays` | **40 calls**, `cumtime` 0.096 s | ✅ count exact |
| `resolve_material_point` is **#1** by `tottime` | **#2.** #1 is `numpy.linalg.inv`, 6 calls, `tottime` **0.214–0.317 s** of a 0.707 s profile | ❌ see §0.3 |

All 9600 cortex calls come from **one** caller, `crosslink_arrays` (`print_callers` attributes
9600/9600 to it). Nothing else in the cortex reaches `resolve_material_point` during a step. That
narrows the whole task to a single loop and is why §5 is as small as it is.

### 0.2 The scaling table — node counts exact, milliseconds not reproducible

| Level | Brief nodes | Re-measured nodes | Brief ms/step | Re-measured ms/step (3 repeats, fresh world each) |
|---|---|---|---|---|
| 1 | 310 | **310** ✅ | 158 | 236.0, 163.4, **47.4** |
| 2 | 790 | **790** ✅ | 170 | 175.2, 431.1, **149.3** |
| 3 | 2710 | **2710** ✅ | 881 | **1356.4**, 4345.4, 1708.0 |

The node counts match to the integer, so this is the same world measured twice and not a different
scenario. **The milliseconds spread 5× at fixed configuration** — level 1 ranges 47→236 ms/step and
level 3 ranges 1356→4345 ms/step, on repeats of the identical recipe. Every one of the brief's three
figures falls inside the corresponding spread.

The cause is not mysterious and is recorded in `CLAUDE.md` §1: **this repository is worked by several
concurrent sessions**, and `docs/ACTIVE_SESSIONS.md` carried six open rows while these numbers were
taken. A wall clock on a contended machine measures the machine. §7 states what this entry does
about it, because the acceptance criterion is a before/after millisecond table and a 5× noise floor
would make that table meaningless.

**A second measurement defect, found by making the first mitigation and watching it fail.** §7's
first attempt reported `time.process_time` instead of wall clock, on the reasoning that CPU time
cannot be inflated by a neighbouring process. That reasoning is wrong here: `process_time` sums
**all threads**, and this environment's BLAS busy-waits its worker threads after a call finishes, so
`process_time` charged the spin-wait to the measurement. It reported ~1420 ms/step at level 1 where
wall clock reported ~47 ms — inflated 30×, and the before/after deltas it produced (1.00×, 1.10×,
1.03×) were noise. The mitigation that actually works is **pinning the thread count to one**
(`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`), after which CPU
time and wall clock agree and the measurement is repeatable. Every number in §7 is taken that way.
This is recorded rather than quietly fixed because §0.4 below was *also* measured before the pin,
and it was wrong by 4× in the other direction.

### 0.3 The diagnosis is withdrawn: the cost is quadratic in the mesh, not fixed

> "노드가 2.5배인데 시간이 1.08배라는 것은 시간이 문제 크기가 아니라 고정 파이썬 오버헤드에 지배된다는 뜻이다."

**This is refused.** The 1.08× it rests on is the noise in §0.2. And the code says the opposite:

```python
def _filament(self, filament_id: int) -> Filament:
    for filament in self.filaments:            # ← linear scan, every call
        if filament.filament_id == filament_id:
            return filament
    raise KeyError(f"the cortex has no filament with id {filament_id}")
```

`resolve_material_point` opens with `self._filament(filament_id)`, so **one resolve is `O(F)` in the
filament count**, and `crosslink_arrays` — which resolves two endpoints per crosslink — is
`O(2·C·F)`. Both `C` and `F` grow linearly with the tessellation. Measured directly:

| Level | Cortex nodes | Filaments `F` | Crosslinks `C` | Scan iterations per `crosslink_arrays` call (`2·C·F/2`) |
|---|---|---|---|---|
| 1 | 126 | 42 | 120 | ~5,040 |
| 2 | 486 | 162 | 480 | ~77,760 |
| 3 | 1926 | 642 | 1920 | **~1,232,640** |

That is a **245× growth** in Python-interpreter work between level 1 and level 3 for a 15× growth in
nodes — quadratic, not flat. The brief's own level-3 observation ("레벨 3에서야 산술이 보이기
시작한다") is better explained by this than by arithmetic appearing: what appears at level 3 is a
quadratic term that was invisible at level 1.

**Why this correction matters rather than being pedantry.** Under the brief's diagnosis the fix
yields a constant saving and the before/after table should improve by a similar factor at every
level. Under the measured diagnosis the saving *grows with the level*, and a table that improved
level 1 only would mean the fix had missed. §7 predicts the shape before measuring it, so the
prediction is falsifiable.

### 0.4 A term that looked like the ceiling, and was an artefact of the machine

**This section is a correction of itself, kept rather than rewritten, because the wrong version is
the more useful warning.**

The first profile taken here — multi-threaded, on the contended machine of §0.2 — put
`aleph/vertical/medium.py::_rebuild_mobility` at **74% of a level-3 step** (`numpy.linalg.inv`
6.632 s + `cholesky` 3.952 s of a 15.604 s profile) and this entry's target subtree,
cortex `crosslink_arrays`, at **17.3%**. On those numbers this task's ceiling was ~17% and the real
bottleneck was a dense `O(n³)` solve in another lane's file.

**Both figures were wrong, in opposite directions, and the cause is the threads.** LAPACK spreads
`inv` and `cholesky` across cores; when five other sessions hold those cores, the LAPACK call spends
most of its wall clock *waiting*, and `cProfile` charges the wait to the call. Re-profiled with the
thread count pinned to one on an otherwise idle machine, the same level-3 step reads:

| Level 3, single-threaded, **before** this entry's change | cumtime | of the 4.497 s profiled step |
|---|---|---|
| cortex `crosslink_arrays` | **3.121 s** | **69.4%** |
| ↳ `resolve_material_point`, 153,600 calls | 2.815 s | 62.6% |
| ↳ ↳ `_filament` alone, 153,600 calls | 1.069 s | 23.8% |
| `medium.py::_rebuild_mobility` (with `inv` at 0.172 s) | 0.286 s | **6.4%** |

So the ceiling is not 17%: **`crosslink_arrays` was more than two-thirds of a level-3 step**, and the
dense mobility rebuild is a twentieth of it. The `O(N·F)` scan of §0.3 is not a detail of the
profile, it *is* the profile — `_filament`'s linear scan alone is a quarter of the whole step.

**What generalises past this entry.** A profile taken on this repository's shared machine, with
BLAS threads unpinned, will systematically over-weight every LAPACK call and under-weight every
Python loop — by 11× and 4× respectively in this one measurement. Both of the brief's diagnoses and
the first version of this section were produced that way. `medium.py` is still not edited here (it
belongs to the PI's own row) and is no longer claimed to be the bottleneck.

### 0.5 A fourth fact the brief did not have: `sf_arc` has nothing to batch

The whole cell builds `SFArcGraph` with **3 bundles and 0 internal joints**. `joint_arrays` therefore
loops zero times, and all 100 of `sf_arc`'s `resolve_material_point` calls come from `endpoint_sink`,
which resolves **one** site by construction and cannot be batched from inside this owner. The brief
names `sf_arc` in the owned set, so the batch method is written and tested there for API symmetry and
for the joint population that a non-whole-cell scenario supplies — but **its contribution to the
before/after table is expected to be exactly zero**, and claiming otherwise later would be
unfalsifiable. The same is true of the 200 `ecm_network` and 100 `sf_arc` `endpoint_sink` calls:
batching an endpoint means batching at the *connector*, which is another lane's file.

Net: of the 10,620 resolve calls in a level-1 step, **9600 + 720 = 10,320 are reachable** from
inside the three owned files, and 300 are not.

---

## 1. What is being changed, in one sentence

Each of the three owners gains a `resolve_material_points` (plural) that takes an array of ids and an
array of material coordinates and returns arrays of `node_a`, `node_b`, `weight` and `position`; the
in-module loops that call the singular form once per site are rewired onto it; **and the singular
form is not touched.**

## 2. Why the singular form stays

Three independent reasons, any one sufficient:

1. **It is the oracle.** §4's acceptance is `np.array_equal` between the batch output and a Python
   loop over the singular form. If the singular form changed in the same commit, the test would be
   comparing two versions of the same new code and would prove nothing.
2. **It has callers this lane does not own** — `aleph/vertical/cortex_surface_coupling.py`,
   `aleph/vertical/connectors_protrusion.py` (the `83c83640` row), `aleph/scenarios/spread.py`
   (Lane W2, live). Changing its signature or its behaviour would reach into three other lanes.
3. **`aleph/vertical/**` is recorded by four closed lanes as the frozen parity reference for the CUDA
   port** (`ALEPH-PORT-3601` §14.7, `-3603` §14.6, and the Lane Track C and `e874a7fe` rows). Every
   ULP figure in `-3601`…`-3619` was measured against these NumPy laws. A lane that edits one is
   contaminating an oracle it cannot see the users of.

## 3. Why the singular form is slow, precisely

Per call, in order:

| Step | Cost |
|---|---|
| `self._filament(id)` | **`O(F)` Python loop.** At level 3, 642 iterations of attribute access and integer compare. |
| `filament.node_indices` | a fresh `np.arange` allocation, per call |
| `self.material_coordinate_um[index]` | a fancy-index gather into a fresh array, per call |
| `np.searchsorted(table, s, 'right')` | one NumPy dispatch (~1 µs of interpreter and dispatch for a 3-element array) |
| `MaterialPoint(...)` | one frozen-dataclass construction |

The arithmetic — one subtract, one divide, two multiplies and an add on three components — is
perhaps 20 ns. **Everything above it is overhead**, and the first row is the one that grows.

## 4. Why the batch result must be bit-identical, operation by operation

This is the load-bearing section. The claim is not "the answers agree to floating-point tolerance".
It is that **the batch performs the same IEEE-754 double operations on the same operands**, so the
results are equal as bit patterns, and `np.array_equal` rather than `pytest.approx` is the correct
assertion.

| Singular | Batch | Why identical |
|---|---|---|
| `k = int(np.searchsorted(table, s, side='right')) - 1` | `k = (table_row <= s).sum(axis=1) - 1` | For a sorted array, `searchsorted(..., side='right')` **is by definition** the count of elements `<= s`. The batch performs the identical `float64 <= float64` comparisons on the identical values. **No arithmetic is involved in choosing `k`**, so there is no rounding that could differ. Padding uses `+inf`, and `inf <= s` is `False` for every finite `s`, contributing 0. |
| `k = min(max(k, 0), node_count - 2)` | `np.clip(k, 0, node_count - 2)` | integer operations |
| `span = float(table[k+1] - table[k])` | `span = t1 - t0` | one `float64` subtraction, same two operands |
| `weight = (s - float(table[k])) / span` | `weight = (s - t0) / span` | one subtraction and one division, same operands. The singular's `float()` calls are `np.float64 → Python float`, which is **exact** — both are IEEE binary64. |
| `position = (1.0-w)*p[a] + w*p[b]` | `position = (1.0-w)[:,None]*p[a] + w[:,None]*p[b]` | NumPy evaluates this as four separate ufunc passes with materialised intermediates in **both** cases — `1.0-w`, then two multiplies, then one add. There is no fused-multiply-add contraction across separate ufunc calls, so each component is the same four binary64 operations on the same operands. |

**The one thing that would break this, and is therefore not done.** If the batch introduced a
*reduction* over sites — a sum, a dot, an `np.add.at` — the summation order would differ from the
loop's and equality would legitimately drop to `approx`. The batch introduces none: every output
element is a function of exactly one input site. That is what makes `array_equal` claimable, and if a
later change adds a reduction here, this row is where the claim has to be revisited.

**Where the batch is *not* free to differ, and a deviation is recorded rather than hidden.** The
singular form raises in a fixed order — `KeyError` from `_filament` first, then the range
`ValueError`. A loop over N sites therefore fails on the *first* offending site in index order, with
that site's error. The batch reproduces exactly this: it computes both failure masks, takes the
lowest failing index, and raises the error that index would have raised. §6 is the control.

## 5. The algorithm

Given `ids` (N,) and `coords` (N,):

1. **Build the per-population lookup once**, not once per site: `starts`, `counts` and an
   `id → row` map, from one `O(F)` pass over `self.filaments`. This replaces `N` scans of length `F`
   with **one** pass of length `F`. It is rebuilt on every call and **not cached**, because the
   population is a plain mutable list with no epoch that is guaranteed to move when it changes, and a
   stale cache here is a wrong force. §8 records the cache as measured-if-needed, not assumed.
2. `rows = map[ids]`; unknown id → the `KeyError` of §4.
3. `total = material_coordinate_um[starts + counts - 1]`; out-of-range → the `ValueError` of §4.
4. **Segmented `searchsorted` without a Python loop.** Because every population here stores its nodes
   as a *contiguous* range (`node_start`, `node_count`; verified for `Filament`, `Fibre` and
   `Bundle`), site `i`'s table is `M[start_i : start_i+count_i]`, and
   `k_i = (that slice <= s_i).count() - 1`. With `L = max(counts)` this is one `(N, L)` boolean
   reduction. **`L` is 3 for the cortex, 5 for the ECM and 4 for `sf_arc`** (§0.3's table), so the
   `(N, L)` intermediate at level 3 is 3840×3 — trivial, and the reason this simple formulation is
   chosen over anything cleverer.
5. Gather, subtract, divide, interpolate — all elementwise, per §4.

**The rejected alternative, and why.** A global monotone table built by adding a per-filament offset
would allow a single `np.searchsorted` over all nodes, `O(N log N)` instead of `O(N·L)`. It is
rejected because the offset addition is *floating-point*: two coordinates that compare distinctly at
their own scale can round to the same value once a large offset is added, silently moving a site to
the neighbouring segment. That is a wrong weight, not a slow one, and it cannot be ruled out by
testing. `O(N·L)` with `L ≤ 5` costs nothing worth that risk.

## 6. Acceptance — every criterion is a number or a raised exception

| # | Control | Criterion |
|---|---|---|
| A1 | batch vs. singular loop, cortex | `np.array_equal` on `node_a`, `node_b`, `weight`, `position` — **equal, not approx** |
| A2 | as A1, `ecm_network` | `np.array_equal` |
| A3 | as A1, `sf_arc` | `np.array_equal` |
| A4 | the whole cell's own crosslink arrays, before vs. after | `np.array_equal` on all six returned arrays |
| A5 | **negative control** — unknown id | batch raises `KeyError` with the **same message** as the singular form |
| A6 | **negative control** — coordinate above `total`, and below `0` | batch raises `ValueError` with the **same message** |
| A7 | **negative control** — ordering | with an unknown id at index `j` and an out-of-range coordinate at index `i < j`, the batch raises what a **loop** would raise, i.e. index `i`'s `ValueError` |
| A8 | empty input | returns empty arrays of the same dtypes; no exception |
| A9 | existing `tests/vertical/**` | all green, unchanged |
| A10 | ms/step, levels 1/2/3, before and after | reported as a table. **If it is not faster, this entry says it is not faster.** |

A4 is modelled on
`tests/scenarios/test_whole_cell.py::test_the_backend_scatter_matches_the_owners_own_material_point_sink`,
which the brief names as the example of the method.

**Result: A1–A9 pass. A10 is §7.2 and §7.3.** `tests/vertical/test_material_point_batch.py` is **24
controls green**; the whole of `tests/vertical` is **1,479 collected, all passed**, with no failure
attributable to this entry. The end-to-end bit-identity check is §7.5. One foreign red, in
`tests/scenarios`, is reported in §9.

## 7. The measurement, and the result

### 7.1 Method, fixed before the numbers were taken

1. **Threads pinned to one** (§0.2). Without this the measurement is dominated by BLAS spin-wait and
   by neighbouring sessions, and both of §0.2's and §0.4's wrong numbers came from not doing it.
2. **The "before" arm is a monkeypatch, not a `git stash`.** Five other sessions had uncommitted work
   in this tree and two commits landed on `main` mid-session, so stashing would compare two trees
   differing in more than this entry's change — and un-stashing writes files this lane does not own.
   Rebinding the three rewired methods to their pre-change loop bodies isolates *exactly* this edit:
   every other byte in the process, including every other lane's working-tree state, is identical
   across the two arms.
3. **The arms alternate** rather than running in blocks, so a machine that gets busier partway
   through penalises both equally.
4. **Minimum over repeats is the headline**, not the mean: contention can only add time.
5. **A correctness gate runs before any timing.** The harness asserts the two arms produce a
   bit-identical `RelaxationReport` and refuses to print a table if they do not — a speed comparison
   between two different answers is not a speed comparison.

**The prediction, recorded before measurement** (§0.3): the saving is a *quadratic* term, so it must
**grow with level**. A result that improved level 1 and not level 3 would mean the fix had missed.

### 7.2 ms/step, whole cell, 9 alternated repeats × 6 steps

| Level | Nodes | Before | After | **Speedup** | Saved |
|---|---|---|---|---|---|
| 1 | 310 | 32.9 ms/step | 23.8 ms/step | **1.38×** | 9.1 ms |
| 2 | 790 | 95.6 ms/step | 50.8 ms/step | **1.88×** | 44.8 ms |
| 3 | 2710 | 374.6 ms/step | 155.7 ms/step | **2.41×** | 218.9 ms |

The speedup grows monotonically with the tessellation — 1.38× → 1.88× → 2.41× — which is the
prediction §0.3 made, and the *opposite* of what the brief's "fixed Python overhead" diagnosis
predicts. An independent earlier run of the same harness gave 1.34× / 1.70× / 2.34×, so the third
digit is not meaningful but the trend and the magnitude are.

### 7.3 `crosslink_arrays` alone, which is what actually changed

The step table is diluted by everything else in the step. This is the call itself, 7 alternated
repeats, with the two arms gated on `np.array_equal` over all six returned arrays first:

| Level | Filaments `F` | Crosslinks `C` | Before | After | **Speedup** |
|---|---|---|---|---|---|
| 1 | 42 | 120 | 1.249 ms/call | 0.099 ms/call | **12.6×** |
| 2 | 162 | 480 | 5.901 ms/call | 0.293 ms/call | **20.1×** |
| 3 | 642 | 1920 | 35.788 ms/call | 1.152 ms/call | **31.1×** |

**The two tables agree, which is the check that they mean something.** A 6-step relaxation makes 40
`crosslink_arrays` calls (§0.1), i.e. 6.67 per step. At level 3 that predicts a saving of
`6.67 × (35.788 − 1.152) = 231 ms/step`; the step table independently measured **218.9 ms/step**. The
rest of the step is untouched, as it should be.

### 7.4 What did not get faster, stated as plainly as what did

- **`sf_arc` contributed nothing to these tables**, exactly as §0.5 said in advance: the whole cell
  builds it with zero internal joints, so the rewired `joint_arrays` loops zero times in this world.
  Its batch method is exercised by its own control and by nothing in the step table.
- **`endpoint_sink` is unchanged.** 300 of a level-1 step's 10,620 resolves go through it one site at
  a time; batching those means batching at the *connector*, which is another lane's file.
- **Level 1 gains 9 ms.** The brief's framing implies fixed overhead is the whole problem at small
  sizes. It is not: at level 1 the loop is about a quarter of the step, and 1.38× is what removing it
  is worth.

### 7.5 Bit-identity, end to end

Beyond the ten controls in §6, the whole cell was relaxed 6 steps before and after the change and the
`RelaxationReport` compared field by field. **Every float is identical to the last bit** —
`energy_initial_pn_um`, `energy_final_pn_um`, `max_residual_force_pn`, `final_step` and the whole
`residual_history` tuple, under `==` and not `approx`. The A/B harness re-checks this on every run
and aborts before printing a timing table if it ever fails.

## 8. Open, and deliberately not taken here

1. **The per-call lookup rebuild is itself `O(F)`.** It is `1/N`-th of what it replaces and is not
   the bottleneck, but caching it against a topology epoch would remove it. Not done: see §5 step 1.
2. **`crosslink_arrays` returns five arrays that do not depend on `positions` at all** — `node_a`,
   `node_b`, `weight`, stiffness and rest length are functions of topology and the material-coordinate
   table only. In principle the whole call could be computed once per topology change instead of once
   per force evaluation, which would beat any batching. It is **not** done here because it is a
   caching question with an invalidation contract, the brief asked for a batch API, and getting
   invalidation wrong produces a wrong force rather than a slow one. Recorded as the larger remaining
   win.
3. **Four other owners carry the identical `resolve_material_point` and the identical `_filament`
   linear scan** — `lamellipodium.py`, `filopodium.py`, `microtubule.py`,
   `intermediate_filament.py`. They are outside this lane's grant and are untouched. Their measured
   cost in the whole cell is small today (260 calls total at level 1) because their populations are
   small, but the `O(N·F)` shape is the same and will surface the same way if those populations grow.
4. **`refresh_steric_pairs` is now the largest single term in a level-3 step** — 0.922 s of `tottime`
   in the 3.741 s post-change profile, where `medium.py`'s mobility rebuild is 0.703 s and
   `crosslink_arrays` no longer appears in the top eight. It is in this lane's own
   `cortex_filaments.py` but it is **not** a `resolve_material_point` caller and so is outside the
   brief's scope; it is named here as the next thing to measure rather than quietly left.

## 9. One foreign red, reported and not fixed

`tests/scenarios/test_whole_cell.py::TestTheRefusalsAreMeasurements::test_nothing_is_refused_any_more_and_the_machinery_that_drove_a_refusal_still_works` **Renamed 2026-08-06:** the control was called the name it carried while `REFUSED` still had an entry in it until `REFUSED` became empty on 2026-08-04; the machinery it drives is unchanged and the test now says so in its name. The entry's citation had drifted, which is the failure `test_named_controls_resolve_to_real_tests` exists to catch.
fails with `DID NOT RAISE CoverageViolation`. It asserts that scheduling `cytosol` trips the coverage
gate — which is exactly what session `5ad3a7ba` (S2) is removing, per `ALEPH-PORT-3625`.
`aleph/vertical/cytosol.py` was clean at this lane's start and is modified in the working tree now.

**Measured rather than assumed**, per `CLAUDE.md` §1: this lane's three files were reverted to `HEAD`
and the test **fails identically**, so the failure is not this entry's. The files were then restored
and the batch controls re-run green. Nothing in `aleph/scenarios/**` or `aleph/vertical/cytosol.py`
was edited. The branch also moved twice mid-session (`7375e9e`, `cf0a72e`); neither commit touches
this lane's three files, verified with `git log -- <path>` before the restore.

---

## Sections completed by Lane W2 (`46143f30`), 2026-08-04 22:xx

**Written by a different lane from the one that did the work**, and that is stated rather than
hidden. `1afc4d25` closed with these mandatory fields unfilled; the entry was otherwise complete and
its code and controls are green. Every claim below is read off the code and the tests in the same
commit, not supplied from outside them — a field that could only be filled by the authoring lane
would have been left open with a note instead.

### Source repository identity, source path and symbol

**None, and the header already says so.** Port class is `NOT A PORT`: `/Users/sw1/ffn_cellsim` was
not opened for this entry. There is no source path, no symbol and no commit, because nothing was
taken. The field is answered rather than skipped so that a reader cannot mistake an absent citation
for an unrecorded one.

### Why source-derived porting beats clean-room

**It does not, and here it is not even the question.** This entry adds a batched resolver beside an
existing one in Aleph's own modules. There is no external implementation to derive from and none was
consulted. The registry's own accepted answer to this field — "RE-DERIVED / nothing was ported" —
applies with the stronger form: nothing was written from any reference at all.

### Positive control

`test_the_batch_resolver_is_bit_identical_to_a_loop_over_the_single_point_one`, with
`test_the_cortex_crosslink_arrays_match_a_loop_over_the_single_point_resolver` and its `ecm` and
`sf_arc` siblings.

The batched resolver and a loop over the surviving single-point resolver are run on the same input
and the results must be **equal**, not close. That is the right bar and not a strict-for-its-own-sake
one: the two paths compute the same expression, so any difference at all means the batch is a second
rule rather than the same rule read faster. `test_the_batch_resolver_holds_bitwise_after_the_geometry_moves`
carries the same claim past a configuration change, which is where a cached index would break first.

### Deliberately failing negative control

Four, and they are failure-*parity* controls rather than mutants — which is the correct shape here,
because this entry's whole risk is that the fast path disagrees with the slow one somewhere:

    test_an_unknown_element_id_raises_the_same_error_as_the_single_point_resolver
    test_an_out_of_range_coordinate_raises_the_same_error_as_the_single_point_resolver
    test_the_batch_reports_the_failure_a_loop_would_have_reached_first
    test_an_empty_batch_returns_empty_arrays_of_the_right_shape_and_dtype

The third is the sharpest. A batch that raises on *some* bad input is not equivalent to a loop; it
must raise on the **first** one a loop would have reached, or an error message names the wrong
element and sends a reader to the wrong filament. The fourth pins the degenerate case, where a
vectorised path most often returns the wrong dtype or a `(0,)` where `(0, 3)` was promised.

### Comments and docstrings to discard

**None.** Nothing was copied, so there is no foreign prose to strip. `resolve_material_point` and its
docstring survive unchanged in all three owners — deliberately, since it is the oracle this entry is
graded against, and rewriting the thing you are grading against is how a parity claim becomes
circular.
