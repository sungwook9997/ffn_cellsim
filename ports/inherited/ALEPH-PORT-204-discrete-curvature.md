# ALEPH-PORT-204 — discrete mean and Gaussian curvature

| Field | Value |
|---|---|
| Lane | L2 (geometry and topology primitives) |
| Target | `aleph/state/curvature.py` |
| Reference consulted | **None.** `ffn_cellsim/ffn_sim/common/surface_manifold.py` carries no curvature; `filament_math.py` is 1-D filament statistics and is unrelated. |
| Port class | **CLEAN ROOM** — nothing ported |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. The two operators have different error characters, and that is the point

### 1.1 Gaussian curvature by angle defect is exact in the integrated sense

`K_v = (2 pi - sum of incident angles) / A_v`. Summing the *defect* (not the density) over a closed
triangulation:

    sum_v (2 pi - sum angles at v) = 2 pi V - pi F

because the interior angles of every triangle sum to `pi`. A closed triangulation has `3F = 2E`, so
`E = 3F/2`, and `V - E + F = 2` gives `V = 2 + F/2`. Substituting:

    2 pi (2 + F/2) - pi F = 4 pi + pi F - pi F = 4 pi

**exactly**, for any closed genus-0 triangulation, at any resolution, however irregular. The result
does not depend on the vertex areas, the vertex positions beyond the angles, or the mesh quality.
Gauss-Bonnet is therefore a machine-precision test here, not a convergence test — measured residual
`7.1e-15` at level 1 rising only to `-2.5e-12` at level 5, which is accumulated round-off over 10242
vertices and nothing else. What converges with refinement is the pointwise density, because that is
the only place `A_v` enters.

This is worth having as a control precisely because it cannot be passed by accident: a topology bug
that changed `V`, `E` or `F` would move the sum off `4 pi` by a multiple of `pi`, not by a small
amount.

### 1.2 Mean curvature by the cotangent Laplacian has no such statement

    L x_i = (1 / (2 A_i)) sum_{j in N(i)} (cot alpha_ij + cot beta_ij) (x_j - x_i)  ->  Delta_S x = -2 H n

so `H = -(1/2) (L x) . n`, positive for a convex body with outward normals. The manifold guarantees
outward normals, which is what makes the sign a testable claim rather than a convention.

## 2. The measurement that changed the design

The intended pairing was cotangent stiffness with **barycentric** vertex areas, on the usual grounds
that barycentric is the lumped mass matrix. Measurement on an icosphere of radius 5 um against the
analytic `H = 1/R` and `K = 1/R^2`:

| level | h [um] | barycentric H median | order | barycentric H max | mixed-Voronoi H max | barycentric K max | mixed-Voronoi K max | order (mV) |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.9114 | 3.214e-02 | — | 1.024e-01 | 2.78e-16 | 1.795e-01 | 7.98e-02 | — |
| 2 | 1.4967 | 9.905e-03 | 1.77 | 1.345e-01 | 1.15e-04 | 1.538e-01 | 2.03e-02 | 2.13 |
| 3 | 0.7537 | 2.337e-03 | 2.10 | 1.430e-01 | 2.88e-05 | 1.479e-01 | 5.50e-03 | 2.09 |
| 4 | 0.3775 | 5.928e-04 | 1.98 | 1.452e-01 | 7.21e-06 | 1.464e-01 | 1.41e-03 | 2.01 |
| 5 | 0.1888 | 1.488e-04 | 2.00 | 1.457e-01 | 1.80e-06 | 1.460e-01 | 3.55e-04 | 2.00 |

With barycentric areas the **median** converges at order 2.00 but the **maximum stalls at 14.5% and
never improves**. The stall is localised: at level 4 the valence-6 vertices have max error 3.06e-03
while the twelve valence-5 vertices have max error 1.452e-01. Those twelve are the icosahedron's
original vertices, cone points of the combinatorics, and refinement does not dilute them because
there are always exactly twelve.

With mixed-Voronoi areas there is no stall at all: the maximum error is nearly five orders of magnitude
smaller at level 5 (0.1457 against 1.80e-6) and still falling, and for Gaussian curvature the maximum converges cleanly at order 2.00
with maximum and median essentially equal.

**Why.** Meyer's mixed-Voronoi area is *derived* by integrating the mean-curvature normal over the
Voronoi cell of the vertex; the cotangent sum and that area are two halves of one derivation.
Barycentric lumping agrees with it to `O(h^2)` where the 1-ring is locally symmetric and disagrees at
`O(1)` where it is not — which is exactly the pentagonal vertices. This corrects a claim written into
ALEPH-PORT-201 §2.4 *before* the measurement was taken; that ledger has been amended in place rather
than silently fixed.

`aleph/state/curvature.py` therefore defaults to `vertex_area_scheme="mixed_voronoi"` while
`TriangulatedSurface.vertex_areas` stays barycentric, because a partition of area (exact sum,
unconditional positivity) and a curvature mass matrix are different jobs.

## 3. Controls shipped

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_curvature.py::test_mean_curvature_converges_to_one_over_r` | `H -> 1/R`, order ~2, measured and asserted `> 1.7` |
| Positive | `...::test_gaussian_curvature_converges_to_one_over_r_squared` | `K -> 1/R^2`, order ~2 |
| Positive | `...::test_gauss_bonnet_is_exact_on_the_icosphere` | `integral K dA == 4 pi` to 1e-9 absolute at every level |
| Positive | `...::test_gauss_bonnet_is_exact_on_a_triangulated_cube` | `4 pi` on a mesh that is not an icosphere and not smooth |
| Positive | `...::test_mean_curvature_is_positive_with_outward_normals` | sign convention |
| Positive | `...::test_cotangent_weights_are_positive_on_an_acute_mesh` | assembly sanity |
| Positive | `...::test_mixed_voronoi_beats_barycentric_at_valence_five` | the §2 measurement, locked in as a test |
| **Negative** | `...::test_barycentric_max_error_does_not_converge` | the failure that motivated the default — asserted to persist, so a future "improvement" that hides it fails |
| **Negative** | `...::test_inverted_sphere_never_reaches_curvature` | globally reversed winding raises at construction, so `H = -1/R` is unreachable |
| **Negative** | `...::test_torus_is_rejected_before_gauss_bonnet_can_lie` | a genus-1 mesh (`chi = 0`) raises `TopologyError`; the Euler check has teeth on a real non-genus-0 surface |
| **Negative** | `...::test_unknown_area_scheme_is_refused` | typed refusal |

## 4. Compliance with PLAN §0.2

Items 1-3 vacuous (clean room). Items 4-6 hold: positive and failing negative controls, entry
written before the code, no `ffn_cellsim` import.
