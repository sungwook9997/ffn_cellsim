# ALEPH-PORT-3659 — the cell list was a Python loop over cells, and it is half the step

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3659` |
| Lane | `b4fad06b` |
| Status | `AUDITED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Follows | `ALEPH-PORT-3632`, which replaced the O(N²) all-pairs scan with this cell list |
| Found by | `cProfile` over 10 steps of the sourced whole cell, for `docs/design/ENGINE_PERFORMANCE_TARGET.md` |

---

## 1. Aleph API

**No signature changes and no new argument.** `CortexFilamentNetwork.refresh_steric_pairs` keeps its
name, its return type and — the whole point — **its exact output**. Only the neighbour sweep's
implementation moves, from a Python loop over occupied cells to array operations.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** — no provider commit is behind this |
| Source path and symbol | **not applicable** — no `ffn_sim` path and no provider symbol was read |
| Read from | this tree's own profile, and `ALEPH-PORT-3632` which wrote the loop |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. **RE-DERIVED** — a rewrite of this repository's own code
against its own measurement.

## 4. Physical or mathematical law represented

**None, and that is the entire claim.** A uniform cell list is a *search*: the pairs it returns are
filtered by a measured distance, and both branches of the steric law have **finite support**
(`inside = distance < cutoff`, and outside it the energy and force are `np.where`-ed to exactly
`0.0`). So the pair list is an index set, not a physical object, and any construction that produces
the same set in the same order produces bit-identical forces.

## 5. Units, domains, singular cases, invariants

- **I1 — the output is BIT-IDENTICAL, as an array, to the loop it replaces.** Not "the same set" and
  not "close": `np.array_equal` on the returned `(P, 2)` int64 array. The final `np.lexsort` fixes a
  total order that does not depend on how the candidates were generated, and the `i < j` filter plus
  one-cell-per-node means no pair can be produced twice — so the order is well defined by the set
  alone.
- **I2 — the small-node path is untouched.** Below `_CELL_LIST_MIN_NODES = 512` the all-pairs scan
  still runs, and it is the oracle.
- **I3 — the empty case still returns `(0, 2)` int64.**
- **I4 — no allocation exceeds the pair count.** The loop it replaces built one `np.repeat` and one
  `np.tile` **per matched cell pair**; the replacement builds each concatenated block once.

## 6. Source evidence class and known retractions

Not a literature claim. Evidence class: **profile and equality test on this tree.**

## 7. Independent oracle or derivation

**`refresh_steric_pairs_all_pairs` already exists for exactly this purpose** and `ALEPH-PORT-3632`
kept it unchanged as the grader for the same claim. This entry is graded against it a second time,
and additionally against the loop it replaces.

**Measured before the code, `cProfile`, 10 steps of the sourced whole cell (93,679 nodes,
ρ = 100 /µm², sweeps = 2), host:**

| | cumulative over 10 steps | per step | share |
|---|---:|---:|---|
| whole relaxation | 22.92 s | 2.29 s | (profiler-inflated from 1.35 s) |
| **`refresh_steric_pairs`** | **11.62 s** | **1.16 s** | **51 %** |
| … of which its own Python (`tottime`) | 5.99 s | 0.60 s | |
| `dict.get` | 2.89 s | 0.29 s | **25,245,758 calls** |
| `np.tile` | 1.42 s | 0.14 s | 946,924 calls |

**Where the 25 million dictionary lookups come from, exactly.** The sweep is

```text
for dx, dy, dz in 27 offsets:
    for key, (begin, end) in starts.items():      # ~93,000 occupied cells
        neighbour = starts.get(key + shift)       # <- one Python dict lookup, each time
```

`27 × 93,000 ≈ 2.5 × 10⁶` lookups **per step**, and `10 × 2.5e6 = 2.5e7` over the profiled run —
which is the measured count to two significant figures. The same loop calls `np.repeat` and
`np.tile` once per *matched* cell pair, which is where the 946,924 `tile` calls come from.

**So the prediction is that this is arithmetic on Python objects, not on arrays**, and vectorising
the sweep removes it. **Predicted: `refresh_steric_pairs` falls to ≲ 0.2 s/step and the whole step
falls by ≈ 1.7–1.9×.** A measured gain below 1.3× refutes the diagnosis, not the code — it would
mean the cost is in the pair *volume* rather than in the Python.

## 8. Positive control

`tests/vertical/test_steric_cell_list.py`:

- `test_the_vectorised_sweep_matches_the_all_pairs_oracle` — I1 against
  `refresh_steric_pairs_all_pairs`, `np.array_equal`, at a node count above the cell-list threshold.
- `test_the_pair_list_is_bit_identical_across_densities` — the same, at three densities, because a
  cell list's failure mode is a boundary case at one occupancy and not a uniform error.
- `test_the_small_node_path_is_unchanged` — I2.
- `test_the_forces_are_bit_identical` — the claim that matters downstream: the accumulated steric
  force on every node, `np.array_equal`, oracle versus cell list. A pair list that is right as a set
  and wrong as a sequence would pass the two controls above and fail this one.

## 9. Deliberately failing negative control

- `test_the_oracle_and_the_cell_list_can_disagree` — a **vacuity control**. It shrinks the steric
  cutoff so the two constructions are asked for genuinely different sets and asserts they then
  differ. Without it, every control above would pass against a `refresh_steric_pairs` that simply
  called the oracle.
- `test_an_empty_pair_list_is_still_shaped` — I3, on a population spread far enough apart that
  nothing is in range.

## 10. Numerical and precision envelope

**No tolerance anywhere in this entry.** Every comparison is `np.array_equal`. Distances are still
computed in float64 on the same `positions` array with the same expression, so the `keep` mask is
the same mask, bit for bit — the change is only in which candidate pairs are handed to it, and the
candidate set is a superset produced by the same integer cell arithmetic.

## 11. Production-backend residency and transfer

**None.** This is host-side index construction. The pair array's dtype, shape convention and
ordering are unchanged, so every device kernel that consumes it is unaffected and the existing
host/device parity controls still grade it.

## 12. Comments and docstrings to discard

Nothing copied. `refresh_steric_pairs`'s docstring keeps `ALEPH-PORT-3632`'s reasoning about *why* a
cell list is exact — that reasoning is unchanged and is what licenses this entry — and gains the
measured cost of the loop being removed.

## 13. Acceptance

`AUDITED`. §8 and §9 pass — `tests/vertical/test_steric_cell_list.py`, **8 passed**, every assertion
`np.array_equal`.

**Re-measured, and §7's magnitude was optimistic. Recorded as measured.**
`scripts/engine_clock_run.py --steps 100 --construction sourced`, host, 93,679 nodes:

| | s/step at 10 | s/step at 25 | speedup |
|---|---:|---:|---:|
| before | 1.3291 | 1.3205 | — |
| after | **0.9280** | **0.9351** | **1.41×** |

**Predicted 1.7–1.9×, measured 1.41×.** §7 set the refutation line at 1.3× and this clears it, so
the diagnosis — that the cost was Python-object arithmetic in the sweep — stands; the *share* the
sweep held was over-attributed, because `cProfile` inflates Python-call-heavy code relative to array
code and the loop was exactly the kind of code it inflates. **A profiler's percentages are a ranking,
not a budget**, and this entry is the measurement of that.

**And the physics did not move, which is the claim that matters.** The engine clock after the change
is `1.045848e-07`, `1.950294e-06`, `9.823787e-06` s at 1, 10 and 25 steps — **bit-identical to the
run before it**, on a trajectory that had already been recorded independently for
`ALEPH-PORT-3658`.

## 14. Honest limits

- **This is lever L1 and L1 is bounded.** `ENGINE_PERFORMANCE_TARGET.md` §3 measures that the sweep
  needs ≈20× and that the method — not the wall time per step — is where most of it has to come
  from. **A 1.8× is a real 1.8× and it is not the target.**
- **The Verlet skin is still off.** `steric_skin_um = 0.0` rebuilds every step. The mechanism to
  reuse the list already exists (`steric_pairs_are_stale`) and this entry does **not** turn it on,
  because the step cap is 2 nm against a ~10 nm steric cutoff — a node crosses a fifth of the cutoff
  in one step, so a skin large enough to buy several steps also enlarges the search radius enough to
  cost what it saves. **That trade-off is arithmetic, not measurement, and it deserves its own
  measured entry rather than a guess here.**
- **The candidate volume is untouched.** A cell of edge `radius` yields 27 cells of candidates where
  a sphere of that radius needs `4π/3 ≈ 4.2` cells' worth; most candidates are rejected by the
  distance test. Reducing that is a different optimisation.
