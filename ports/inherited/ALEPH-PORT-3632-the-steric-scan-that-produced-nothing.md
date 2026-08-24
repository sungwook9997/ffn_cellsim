# ALEPH-PORT-3632 — the O(N²) steric scan, and the fixture where it produced nothing

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3632` |
| Lane | `46143f30` Lane W2 (lead) |
| Status | `PROPOSED` |
| Written | `2026-08-05` |
| Port class | `NOT A PORT` — a uniform cell list is textbook; nothing was read from `/Users/sw1/ffn_cellsim` |
| Aleph target | `aleph/vertical/cortex_filaments.py`, `tests/vertical/test_cortex_steric_cell_list.py` (new) |
| Supersedes | Nothing. The all-pairs scan survives **unchanged** as `refresh_steric_pairs_all_pairs` and is the oracle this entry is graded against. |

---

## 1. What it cost

`refresh_steric_pairs` built every node pair with `np.triu_indices(n, k=1)` and then filtered. Its
own docstring said so: *"Exact and O(N²): every node pair is tested. There is no neighbour list and
no hash grid."*

| cortex nodes | pairs allocated | int64 index memory | time |
|---:|---:|---:|---:|
| 486 | 1.2e5 | 1.9 MB | 5.2 ms |
| 1,926 | 1.9e6 | 30 MB | 90.6 ms |
| 7,686 | 3.0e7 | 472 MB | **1,275.7 ms** |
| 30,726 | 4.7e8 | **7.5 GB** | — |
| 491,526 | 1.2e11 | **1.9 TB** | `MemoryError` in `__post_init__` |

Called from three places — the constructor, `commit`, and `accumulate(DISCOVER)` — so it is on the
step path, not only the build path.

## 2. The measurement that makes it absurd rather than merely expensive

**On the whole cell's own fixture the scan returns ZERO pairs at every level.**

The steric cutoff is `0.0112 µm`. Nodes closer than that are all on the *same filament* and excluded
by design. Different filaments sit ~0.5 µm apart, forty times the cutoff.

So at level 4 the engine spent **1.28 seconds per call allocating 472 MB to produce an empty array.**

That is worth naming as its own finding: the cost was not buying accuracy. It was buying nothing.

## 3. Why a cell list is exact here and not an approximation

**The steric law is finite support in BOTH branches**, and that was checked by reading the law
rather than inferred from its name (`cortex_filaments.py:699-712`):

    cutoff = 2.5σ if allow_steric_attraction else wca_cutoff_um(σ)
    inside = distance < cutoff
    energy = np.where(inside, ..., 0.0)
    magnitude = np.where(inside, ..., 0.0)

Outside the cutoff the energy and the force are `0.0` exactly — `np.where`, not a small number, in
the WCA branch and in the truncated-attraction branch alike. A cell list built at the same radius
therefore drops only terms that were already exactly zero.

`test_the_steric_law_is_finite_support_which_is_why_this_is_exact` drives a pair at
`cutoff × 1.000001` and requires energy and force to be exactly zero, so the premise is measured and
a future law with a tail turns it red.

## 4. Bit-identity, and why the set is not enough

The cell list returns the **same pairs in the same order**. Matching the set would not be enough:
`np.add.at` accumulates in list order and float addition is not associative, so a different sequence
gives forces that differ in the last bits — close, and not the same.

Row-major `(i, j)` order is what `triu_indices` produces, and `np.lexsort((j, i))` reproduces it.

Below `_CELL_LIST_MIN_NODES = 512` the small path **calls the oracle itself**, so the two agree by
construction rather than by argument, and a control asserts that the small path takes that route.

## 5. Measured

Real cortex, all bit-identical:

| level | cortex nodes | pairs | cell list | all pairs | speedup |
|---:|---:|---:|---:|---:|---:|
| 2 | 486 | 0 | 4.98 ms | 5.24 ms | 1.1× |
| 3 | 1,926 | 0 | 17.79 ms | 90.64 ms | **5.1×** |
| 4 | 7,686 | 0 | 57.14 ms | 1,275.72 ms | **22.3×** |

A synthetic **dense** configuration, where 2.16M pairs actually survive, gives only 3.7× at 15,000
nodes — and that is the honest shape of the result: the cell list wins by not building what is
thrown away, so it wins most where most is thrown away. Both regimes are reported rather than only
the flattering one.

## 6. Units, domains, singular cases, invariants

Units: length µm. The cell edge is `max(cutoff + skin, 1e-3 µm)`; the floor stops a degenerate cutoff
asking for an unbounded number of cells.

| invariant | how it holds |
|---|---|
| same pairs as the scan | 27-cell neighbour sweep covers every pair within one cell edge, and the edge **is** the search radius |
| same order | `np.lexsort((j, i))` is row-major `(i, j)`, which is `triu_indices` order |
| each pair once | `i < j` filter after the sweep, which halves rather than filters |
| same-filament pairs excluded | unchanged, `filament_of_node` |
| small case identical by construction | the small path calls the oracle |

**Singular case:** an empty cell list — every node isolated — returns `(0, 2)` rather than a
zero-length 1-D array, so downstream indexing is unchanged.

## 7. Positive control, negative control

**Positive** — `test_the_pair_list_is_bit_identical_to_the_all_pairs_scan` at 200/400/800 filaments,
and `test_the_forces_are_bit_identical_too`, which drives the actual force law through both lists
and requires `np.array_equal`. Both refuse to run on an empty comparison: the fixture is asserted to
produce pairs first, because the whole cell's own cortex produces none and a control written against
it would compare two empty arrays and prove nothing.

**Negative** — `test_same_filament_pairs_are_still_excluded`, and the finite-support control in §3,
which is the premise rather than the conclusion.

## 8. Numerical envelope

Host `float64`. The identity claims are exact (`np.array_equal`, `==` on the energy), because the
arithmetic is the same expression on the same pairs in the same order. A tolerance would hide
exactly the reordering this entry has to rule out.

## 9. Residency

Host NumPy. The neighbour sweep is a gather over sorted cell keys and is the shape a device port
wants, but nothing here is on the device and none is claimed.

## 10. Open

**The skin is still unused.** `steric_skin_um` is added to the radius and the list is rebuilt on
every accepted step regardless. Verlet reuse — rebuild only when a node has moved more than
`skin/2` — is the next win, and it interacts with the transaction: a rejected candidate must not
update the reference. The medium's `committed_surface_reference_um` is the rule to follow.
