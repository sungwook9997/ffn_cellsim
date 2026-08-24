# ALEPH-PORT-202 — locality queries and resolution-invariant neighbourhoods

| Field | Value |
|---|---|
| Lane | L2 (geometry and topology primitives) |
| Target | `aleph/state/locality.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/common/surface_manifold.py` (`nearest_patch`, `patch_kring`, `kring_for_reach`) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` |
| Source file digest | `sha256:1c070b60f8e36789e42dbb7d915e380ce0a0602386a6b509268fd2d5df34f730` |
| Source read | **The commit, not the working tree.** `git status --porcelain` reports the file clean, and the worktree digest equals `git show be0e5876:...` byte for byte. The repository as a whole carries 31 uncommitted changes, so this was checked per file rather than assumed. |
| Source test status | **NO LIVE TESTS.** |
| Port class | **RE-DERIVED**, with the reach rule re-stated and the coverage claim converted from an assertion into a test |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. What the reference does

Three things: a `cKDTree` nearest-centroid lookup; a BFS k-ring over face adjacency; and

    k = ceil((reach + 2 * max_circumradius) / min_centroid_step) + safety_rings

justified in prose as "the smallest centroid-to-neighbour-centroid step is a lower bound on the
geodesic advance per face-hop, so k hops guarantee at least `k * min_centroid_step` of geodesic
radius".

## 2. Independent re-derivation

### 2.1 The reference's stated justification does not hold

The quoted argument is backwards. Hop distance `h` from face `f` to face `g` means there *exists* a
path of `h` hops, which bounds the centroid separation *above* by `h * step_max`, not below by
`h * step_min`. What coverage actually requires is the converse statement — every face whose
centroid lies within `reach` of the seed centroid is reachable in at most `k` hops — and that is a
quasi-geodesic property of the dual graph, not something the minimum step implies.

The two padding terms are still the right ones and the ceiling-plus-safety-ring structure is right.
What is wrong is the claim that the result is *proved* — and, as §2.2 shows, the per-hop advance it
divides by is also too large. Aleph therefore ships the rule as a **calibrated bound with an
executable postcondition**, not as a theorem.

### 2.2 The rule that was written first, and why it was wrong

The obvious per-hop advance is the *typical* centroid step. Two triangles sharing an edge have
centroids separated by `2 * (1/3) * h_tri`; for an equilateral triangle of side `l` that is
`l / sqrt(3) ~ 0.5774 l`, and the measured mean centroid step on an icosphere converges to exactly
that (ratio 0.9928, 0.9973, 0.9988, 0.9994 at levels 2 to 5). This was implemented, with the padding
below, and its justification was written into this ledger.

**The coverage test then failed it.** At level 5 with a 1 um reach, a ring sized from the mean step
left four faces at 0.9989 um outside a neighbourhood claiming to cover 1.179 um. The failure is
recorded here rather than quietly corrected, because the whole reason the postcondition is
executable is that a docstring cannot catch this.

The error: hop distance is not proportional to Euclidean distance along the *worst* path. A shortest
dual-graph path can be forced to rotate about a vertex instead of running across a triangle strip,
and rotating around a valence-6 vertex costs three hops to net roughly one edge length. Measuring the
true worst case — the minimum of (Euclidean distance)/(hop count) over every face pair within range,
over many seeds — gives **0.3647, 0.3421, 0.3339, 0.3326** times the mean edge length at levels 2 to
5, converging to `1/3`.

### 2.3 Aleph's rule as shipped

    hop_step              = mean_edge_length / 3
    ring_for_reach(reach) = ceil((reach + 2 * max_face_radius) / hop_step) + safety_rings  [default 1]
    covered_reach(k)      = k * hop_step - 2 * max_face_radius

`max_face_radius` is the exact bound on how far a point in a face can be from that face's centroid
(a triangle is the convex hull of its vertices, so the maximum is attained at a vertex). Both the
query point and the target can sit that far off their own centroids, hence the factor 2. The
reference called this quantity a circumradius; it equals the circumradius only for an equilateral
triangle, so Aleph does not use that name.

`hop_step` and `max_face_radius` are measured from the mesh; `reach` is supplied by the caller
(physics). Nothing is tunable.

**The price is real and is not hidden.** Because the bound is worst-case while a typical direction
advances `sqrt(3)` times faster, the ring covers about 3x the area of the reach ball in the limit
(measured 24.2, 10.1, 5.1, 3.2 times the exact spherical-cap area at levels 2 to 5, converging to
`(sqrt 3)^2 = 3`). For a broad-phase candidate search that is the correct way to be wrong: a missed
neighbour is a silent physics error, an extra candidate is arithmetic.

### 2.4 Why this is resolution-invariant

`hop_step` and `max_face_radius` both scale as `l ~ 2^-n`, so `k` grows as `2^n` and
`covered_reach(ring_for_reach(r)) - r` — the overshoot — shrinks as `O(l)` (measured 0.616, 0.310,
0.159, 0.083 um at levels 2 to 5 for a 1 um reach, always inside the guaranteed bound of
`(1 + safety_rings) * hop_step`). A coarse mesh covers the same physical ball in fewer, larger hops;
a fine mesh in more, smaller hops; the covered *physical* radius converges to `reach` from above and
the ring's physical area converges to a fixed multiple of the reach ball's. That is what the tests
measure, at fixed `reach`, across subdivision levels.

### 2.5 The coverage claim is a test, not a comment

`verify_reach_coverage()` takes seed faces and checks directly that every face whose centroid lies
within `reach` of the seed centroid is in `k_ring(seed, ring_for_reach(reach))`. The test suite runs
it over sampled seeds at several levels. If the bound is ever violated on a mesh, the test says so
rather than a docstring asserting it cannot happen.

### 2.6 `nearest_triangle` is nearest-*centroid*, and says so

The KD-tree returns the face with the closest centroid, which is not always the face containing the
closest surface point — the discrepancy is bounded by `max_face_radius` and matters only for points
near a face boundary. The reference buried this in an application-specific aside about cortex shell
thickness. Aleph states the bound in the docstring and exposes the quantity, so a caller can decide
whether it matters. No cell-biology assumption is baked in.

## 3. Controls shipped

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_locality.py::test_nearest_vertex_recovers_the_vertex_itself` | querying a vertex returns its own index |
| Positive | `...::test_nearest_triangle_recovers_the_centroid_itself` | querying a centroid returns its own face |
| Positive | `...::test_nearest_triangle_of_a_slightly_displaced_centroid_is_stable` | robustness to a normal offset |
| Positive | `...::test_k_ring_zero_is_the_seed_alone` | `k_ring(f, 0) == {f}` |
| Positive | `...::test_k_ring_one_is_the_seed_plus_three_neighbours` | 3-regular dual graph |
| Positive | `...::test_k_ring_is_monotone_and_saturates` | nested, and reaches the whole mesh |
| Positive | `...::test_vertex_k_ring_one_is_the_incident_vertices` | vertex-graph companion |
| Positive | `...::test_typical_hop_step_matches_the_measured_centroid_step` | `l / sqrt(3)` is the typical step — the trap, verified as such |
| Positive | `...::test_worst_case_hop_advance_is_one_third_of_an_edge` | the §2.2 measurement, re-derived inside the test |
| Positive | `...::test_reach_coverage_holds_on_sampled_seeds` | the executable postcondition, at four levels |
| Positive | `...::test_covered_reach_is_at_least_the_requested_reach` | `covered_reach(ring_for_reach(r)) >= r` |
| Positive | `...::test_ring_overshoot_is_bounded_by_two_hop_steps` | overshoot is `O(mean edge length)` and shrinks 7.4x from level 2 to 5 |
| Positive | `...::test_ring_physical_area_converges_across_levels` | ring area converges to ~3x the exact spherical-cap area, with shrinking increments |
| **Negative** | `...::test_fixed_ring_count_is_not_resolution_invariant` | the naive fixed-`k` alternative loses >3x of covered physical area per level; the derived one converges |
| **Negative** | `...::test_fixed_ring_count_undercovers_the_reach_ball` | fixed `k` misses 1052 faces inside the reach ball at level 5, where the derived count misses none |
| **Negative** | `...::test_negative_k_is_rejected` | typed refusal |
| **Negative** | `...::test_out_of_range_seed_is_rejected` | typed refusal |
| **Negative** | `...::test_non_positive_reach_is_rejected` | typed refusal, including `0.0`, `nan` and `inf` |
| **Negative** | `...::test_bad_point_shape_is_rejected` | typed refusal |

## 4. Compliance with PLAN §0.2

1. No source text carried; vocabulary check in `tests/state/test_provenance_hygiene.py`. ✔
2. Law re-derived; the reference's stated justification found unsound and replaced (§2.1). ✔
3. New prose. ✔
4. Positive and failing negative controls (§3). ✔
5. Entry written before the code. ✔
6. No `ffn_cellsim` import. ✔

## 5. Verdict on the reference

**Shape of the formula usable, per-hop advance too optimistic, justification unsound.** The two
padding terms are the right ones and the ceiling-plus-safety-ring structure is right. Two problems:

1. The justification proves something other than what coverage requires (§2.1).
2. The per-hop advance is too large. The reference divides by the *minimum centroid step*, measured
   at `0.494` of the mean edge length on an icosphere, where the true worst-case advance is `0.333`
   of it — optimistic by a factor of 1.49. Choosing the minimum *step* is the right instinct applied
   to the wrong quantity: the binding constraint is not the shortest single hop but the least
   productive sequence of hops, and a path rotating about a vertex makes full-length hops that go
   almost nowhere.

Whether the reference's rule undercovers in any given case depends on how much its
`2 * max_circumradius` padding happens to absorb, which is not a property anyone was tracking.
Nothing at source tested it. It was found in Aleph only because the coverage claim was written as an
enumeration rather than as a sentence.
