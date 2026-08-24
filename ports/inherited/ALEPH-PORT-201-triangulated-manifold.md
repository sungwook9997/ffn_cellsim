# ALEPH-PORT-201 — closed triangulated surface manifold (icosphere)

| Field | Value |
|---|---|
| Lane | L2 (geometry and topology primitives) |
| Target | `aleph/state/manifold.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/common/surface_manifold.py` |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (worktree HEAD, 2026-07-30) |
| Source file digest | `sha256:1c070b60f8e36789e42dbb7d915e380ce0a0602386a6b509268fd2d5df34f730` |
| Source read | **The commit, not the working tree.** `git status --porcelain` reports the file clean, and the worktree digest equals `git show be0e5876:...` byte for byte. The repository as a whole carries 31 uncommitted changes, so this was checked per file rather than assumed. |
| Source test status | **NO LIVE TESTS.** `test_surface_manifold.py` was deleted with the retired subsystem. Nothing in the reference is validated. |
| Port class | **RE-DERIVED** (no source text carried; algorithm independently reconstructed and, in two places, replaced) |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. What the reference does

Builds an icosphere by recursive 4-into-1 subdivision from a hard-coded 12-vertex / 20-face
icosahedron table; derives per-face centroids, normals, a tangent frame, and flat areas; builds a
face-adjacency array by hashing undirected edges; exposes a nearest-centroid KD-tree lookup and a
BFS k-ring; and derives a k-ring count from a physical reach. Normals are made outward by flipping
each face normal individually whenever `n . c_hat < 0`.

It carries no mechanics, which is the one structural decision Aleph keeps.

## 2. Independent re-derivation, and where Aleph departs

### 2.1 Base icosahedron — constructed, not tabulated

The reference hard-codes 12 coordinate triples and a 20-row face table. A table is exactly the kind
of artifact that cannot be checked by reading it, and it is untested at source. Aleph instead
*derives* the icosahedron:

- The 12 vertices are the cyclic permutations of `(0, +/-1, +/-phi)` with `phi = (1+sqrt 5)/2`. These
  are generated programmatically from that rule and normalised to the unit sphere.
- Adjacency is recovered geometrically: on the regular icosahedron every vertex has exactly 5
  nearest neighbours at a common distance, so the 5 smallest pairwise distances per vertex give the
  edge set. The construction asserts 30 edges.
- The 20 faces are the mutually-adjacent vertex triples, enumerated combinatorially. The
  construction asserts exactly 20.

Nothing is trusted to a literal; every count is checked.

### 2.2 Orientation — global and consistent, not per-face

**This is a correctness upgrade, not a restyle.** The reference orients each face independently by
the sign of `n . c_hat` (normal against the ray from the origin to the face centroid). That test:

- assumes the surface is star-shaped about the origin, which the icosphere is but a *deformed*
  membrane (the actual use case for this lane) need not be;
- flips the computed *normal vector* while leaving the *vertex winding* untouched, so the stored
  triangle winding and the stored normal can disagree, and any downstream consumer that recomputes a
  normal from `tris` gets the opposite sign;
- never establishes that neighbouring faces agree, so it cannot detect a locally inverted triangle.

Aleph does it in two stages, both of which apply to any closed orientable surface:

1. **Consistency.** Breadth-first traversal of the face-adjacency graph, flipping a face's winding
   whenever it traverses a shared edge in the same direction as its already-visited neighbour. On a
   closed orientable surface the result is a globally consistent winding: every interior edge is
   traversed exactly once in each direction. This is checked, not assumed.
2. **Sense.** The signed volume `V = (1/6) sum_f v0 . (v1 x v2)` is then computed once. If `V < 0`
   the whole mesh is reversed; if `V == 0` the mesh is rejected. Outwardness is a single global
   property of a consistently-oriented closed surface, so one sign decides it.

Normals are then taken *from the winding only*. There is no per-face flip anywhere in Aleph's
manifold, so winding and normal can never disagree.

**Independent confirmation, and the size of the defect.** The ports lane audited the reference by
running it and reached the same verdict from the other direction: translating the shell by `3R` — a
pure isometry, which cannot change which way is out — flips 106 of its 320 normals, and its
`fit_to_cloud` centres on an arbitrary cloud centroid, so any cell not sitting at the origin gets a
silently corrupted normal field. Their verdict is **RE-DERIVE, do not port**.

Aleph reproduces that number exactly, on its own mesh, as a test. Applying the rejected
`sign(n . c_hat)` rule to a level-2 icosphere translated by `3R` disagrees with the true outward
normal on **106 of 320 faces (33.1%)**, rising to 122 (38.1%) at `5R`; at level 3 it is 428 of 1280
and 516 of 1280. Aleph's own normals flip on **none** of them and the enclosed volume changes by
`0` to `3.4e-16` relative. Roughly a hemisphere's worth of faces is what the origin-ray rule always
gets wrong once the surface is displaced — it is asking which way is away from the *origin*, not
which way is out of the *surface*.

The reason this shipped is that the rejected rule *is* correct on a sphere centred at the origin,
which is the only configuration an unmoved test would ever build. Hence
`test_orientation_survives_a_large_translation`, and its companion
`test_the_origin_ray_orientation_test_would_have_failed_that_translation`, which measures the
disagreement so the first test cannot quietly become vacuous.

### 2.3 The divergence-theorem identity is exact, not approximate

For a closed surface, `V = (1/3) * integral over S of (x . n) dA`. On a flat triangle `n` is
constant and `(x - c) . n = 0` for every `x` in the triangle, so `integral (x . n) dA = (c . n) A`
*exactly*. Therefore

    (1/3) sum_i (c_i . n_i) A_i  ==  (1/6) sum_f v0 . (v1 x v2)

holds to floating-point round-off, not to discretisation error. Aleph asserts it at `1e-12`
relative and uses it as the outward-orientation invariant. The reference states no such identity.

### 2.4 Per-vertex area — barycentric, and why

Aleph's primary `vertex_areas` is barycentric: `A_v = (1/3) sum_{f containing v} A_f`, because as a
*partition of area* it is the better object:

1. It partitions the total area exactly: `sum_v A_v == sum_f A_f` to round-off (measured: 0.0 at
   every level tested). The Voronoi construction does not, once obtuse triangles force the Meyer
   "mixed" fallback.
2. It is unconditionally positive. The unmodified Voronoi area is negative for obtuse triangles.

`vertex_areas_mixed_voronoi()` is shipped alongside, and **an initial claim written into this ledger
before the code was measured turned out to be wrong.** The claim was that barycentric is the mass
matrix that pairs with the cotangent Laplacian, and that the two schemes agree to within about 1% on
an icosphere. Both parts are false and the measurement is recorded here rather than quietly fixed:

- The two schemes agree at valence-6 vertices at `O(h^2)` (relative difference 9.9e-3, 2.3e-3,
  5.9e-4, 1.5e-4 at levels 2 to 5) but differ by ~14.6% at the twelve valence-5 vertices, and that
  gap does **not** shrink with refinement.
- Pairing the cotangent sum with mixed-Voronoi areas rather than barycentric areas improves the
  worst-vertex mean-curvature error from a stalled 14.5% to 1.8e-6 at level 5. Meyer's mixed area is
  *derived* by integrating the mean-curvature normal over the Voronoi cell; the cotangent sum and
  that area are two halves of one derivation, whereas barycentric lumping agrees with it only where
  the 1-ring is locally symmetric.

So `aleph/state/curvature.py` defaults to `mixed_voronoi` and the manifold's own `vertex_areas`
stays barycentric. The tests assert the O(h^2) valence-6 agreement and the non-convergence at
valence-5 explicitly, so neither fact can quietly regress. This is ALEPH-PORT-204's core measurement
and is cross-referenced there.

### 2.5 Dropped from the reference

- `fit_to_cloud` — an inverse-angular-distance shell fit to a bead cloud. It is a *modelling*
  operation, not geometry, it silently assumes a star-shaped cloud about its own centroid, and it
  belongs to whichever lane owns membrane/cortex coupling. Out of scope for L2 and not ported.
- `min_centroid_step` / `max_circumradius` as manifold attributes — moved to
  `aleph/state/locality.py`, which is the only consumer (see ALEPH-PORT-202).
- The `1.0e-300` denominators used to make degenerate faces silently survive. Aleph rejects a
  degenerate face at construction instead of producing a meaningless frame for it.

## 3. Units

Lengths in micrometres, areas in square micrometres. This is the convention the
membrane-cortex-pressure vertical needs so that wavevectors come out in 1/um, matching the interface
agreed with lane L4 (`validation/analytic/helfrich.py`). The module states this and does not import
anything from `validation/`.

## 4. Controls shipped

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_manifold.py::test_face_and_vertex_counts_match_closed_form` | `20*4^n` faces, `10*4^n + 2` vertices, levels 0-5 |
| Positive | `...::test_euler_characteristic_is_two` | `V - E + F == 2` at every level |
| Positive | `...::test_every_edge_has_exactly_two_faces` | manifoldness |
| Positive | `...::test_face_adjacency_is_three_regular_symmetric_connected` | dual graph is 3-regular, symmetric, one component |
| Positive | `...::test_divergence_theorem_recovers_enclosed_volume` | `(1/3) sum (c.n)A == signed volume > 0` to 1e-12 rel |
| Positive | `...::test_enclosed_volume_converges_to_the_sphere_volume` | converges to `4/3 pi R^3` from below, order ~2 |
| Positive | `...::test_total_area_converges_to_four_pi_r_squared` | converges to `4 pi R^2` from below, order ~2 |
| Positive | `...::test_orientation_survives_a_large_translation` | translating by 3R, 5R and -8R flips no normal and changes no volume or area — the test the reference lacked |
| Positive | `...::test_orientation_survives_a_rotation` | the other isometry: volume and area invariant, normals rotate with the surface |
| Positive | `...::test_frames_are_right_handed_orthonormal` | `e1 x e2 == n`, all unit, all mutually orthogonal |
| Positive | `...::test_barycentric_vertex_areas_partition_the_total_area` | `sum_v A_v == sum_f A_f` |
| **Negative** | `...::test_single_inverted_triangle_is_rejected` | one flipped winding raises `OrientationError` |
| **Negative** | `...::test_the_origin_ray_orientation_test_would_have_failed_that_translation` | the rejected `sign(n . c_hat)` rule disagrees on 106 of 320 faces at 3R — the reference's shipped defect, reproduced |
| **Negative** | `...::test_a_reflection_is_rejected` | an improper isometry inverts the sense and must not validate |
| **Negative** | `...::test_globally_reversed_winding_is_rejected` | inside-out sphere raises `OrientationError` |
| **Negative** | `...::test_non_manifold_edge_is_rejected` | edge with 3 faces raises `NonManifoldError` |
| **Negative** | `...::test_genus_one_surface_is_rejected` | a torus (`chi = 0`) raises `TopologyError` — the Euler check has teeth on a real non-genus-0 surface |
| **Negative** | `...::test_open_surface_is_rejected` | boundary edge raises `NonManifoldError` |
| **Negative** | `...::test_degenerate_triangle_is_rejected` | zero-area face raises `DegenerateFaceError` |
| **Negative** | `...::test_orientation_repair_recovers_the_inverted_case` | `orient_faces_outward` fixes what the check rejected |

## 5. Compliance with PLAN §0.2

1. All original comments and docstrings stripped; nothing carried. Verified by
   `tests/state/test_provenance_hygiene.py`, which greps the shipped modules for reference-only
   vocabulary and identifiers (`H.7`, `Sanity Gate`, `HOOMD`, `bead`, `SurfaceManifold`,
   `patch_kring`, `nearest_patch`, `kring_for_reach`, `fit_to_cloud`, `tri_adj`,
   `min_centroid_step`, `max_circumradius`, `ffn_sim`, `ffn_cellsim`). The word "cortex" is
   deliberately *not* on that list: it names Aleph's own first vertical
   (membrane-cortex-pressure), so banning it would be banning Aleph's vocabulary rather than the
   reference's. ✔
2. Law re-derived independently (§2). Two defects in the reference identified and fixed (§2.2). ✔
3. New prose written for Aleph's situation. ✔
4. Positive and deliberately failing negative controls shipped (§4). ✔
5. This entry written before the code landed. ✔
6. No import of `ffn_cellsim` anywhere; passes with it absent from `sys.path`. Asserted by
   `tests/state/test_provenance_hygiene.py::test_no_forbidden_imports`. ✔

## 6. Verdict on the reference

**Partially correct, and unsafe to trust as-is.** The subdivision scheme, the edge-hash adjacency,
and the mechanics-free scope are sound. The orientation handling is wrong for the deformed-membrane
case this lane exists to serve, and it cannot detect an inverted triangle at all — the exact failure
Aleph's negative control is built around. Given the source has no live tests, this was not
detectable at source.
