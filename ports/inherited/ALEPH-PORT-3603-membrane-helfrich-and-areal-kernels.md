# ALEPH-PORT-3603 — the membrane's bending and areal energies as Warp kernels

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3603` |
| Lane | `46143f30` Lane G2 (CUDA track), Track G task **G2** |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Depends on | `ALEPH-PORT-3601` (the law-level parity harness). G2 **uses** that harness and does not build a second one. |
| Exists because | G1 put one law on a kernel and stated its own limit (`ALEPH-PORT-3601` §14.2): the tether is *"a central, pairwise, unilateral spring with no connectivity … it does not exercise a mesh operator, a curvature stencil, or anything with a face loop — so nothing here is evidence that the harness's shape survives the membrane's Helfrich term."* This entry is that evidence, or its absence. It is also the first kernel in this project with a **stencil** rather than a per-site loop. |

---

## 1. Aleph API

```python
from aleph.runtime.law_kernels import (
    ArealVariant,               # TRUE | WRONG_AREA_HALF | WRONG_CORNER_GRADIENT
    BendingVariant,             # TRUE | WRONG_ABS_COTANGENT | WRONG_COTANGENT_HALF
                                #      | WRONG_VERTEX_AREA_LUMPING | WRONG_FROZEN_COTANGENT
    PositionPrecision,          # GLOBAL_F32 | LOCAL_F32 | POSITIONS_F64   -- see section 10
    load, loaded_warp,          # unchanged, extended to compile the membrane kernels
)

from aleph.runtime.law_cases import (
    MEMBRANE_BENDING_RIGIDITY_PN_UM,   # 0.086 pN.um, the vertical's own value
    DEFAULT_AREAL_ULP_BUDGET,
    DEFAULT_BENDING_ULP_BUDGET,
    areal_law_case,             # sigma * A_total, per face, as a LawCase
    bending_law_case,           # (kappa/2) sum_i |K_i|^2 / A_i, as a LawCase
    membrane_state,             # reproducible closed triangulated surface + material
    membrane_kernel_forces,     # both kernels' assembled vertex force [pN], for the oracles
    recovered_tension_from_kernel_pn_per_um,   # the dilational virial, read through the kernel
    SWEEP_GRID,                 # the twelve (level, sigma) configurations of the tension sweep
)
```

Nothing outside this list is covered. In particular this entry does **not** authorise any change to
`aleph/vertical/**` — **`aleph/vertical/membrane.py` is the frozen parity reference and is read,
never edited** — nor to `aleph/runtime/parity.py`, `aleph/runtime/backend.py`,
`aleph/runtime/warp_kernels.py`, `aleph/scenarios/**` or `aleph/viz/**`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read.** No file under that tree was opened by this lane. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

The law being ported is **Aleph's own**: `aleph/vertical/membrane.py`, whose own provenance is
`ALEPH-PORT-1101-helfrich-surface-mechanics.md`. This entry moves that already-ported law onto a
kernel; it does not re-port it from anywhere. The reference is the module, bit for bit, called
rather than re-implemented.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported from another project.** The discretisation is fixed by
`aleph/vertical/membrane.py` — the barycentric lumping, the cotangent identity, the association
order of every product — because the whole purpose is to agree with it. A clean-room second
implementation would be a different discretisation and would disagree by discretisation error, which
is not what a parity gate measures.

## 4. Physical or mathematical law represented

Two energies, deliberately split into **two separate `LawCase`s** rather than one, because they have
different exact identities and a single total would hide both.

### (a) The areal response, per face

```
A_f = (1/2) |e1_f x e2_f|                          face area [um^2]
E_area = sigma * sum_f A_f                          [pN.um]
f_v = -sigma * sum_{f,k : tri(f,k)=v} grad_{p_k} A_f
grad_{p_k} A_f = (1/2) (p_{k+1} - p_{k+2}) x n_f    [um]
```

Degree-**two** homogeneous in the vertex positions. That is the property Laplace's law is read from.

### (b) The Helfrich bending term, cotangent-Laplacian form

```
K_i = (1/2) sum_{j in N(i)} (cot alpha_ij + cot beta_ij) (x_i - x_j)     ~ 2 H_i A_i n_i   [um]
A_i = (1/3) sum_{f > i} A_f                                              barycentric lumping [um^2]
E_bend = (kappa/2) sum_i |K_i|^2 / A_i                                   [pN.um]
```

with the exact gradient the module derives — the chain rule through `grad A_f`, through
`grad_{p_m} cot_k = (1/|N|) d(u_k.v_k)/dp_m - (cot_k/|N|) 2 grad_{p_m} A_f`, and the direct term.

**The three assumptions carried over verbatim from the reference**, because a kernel that quietly
dropped one would be a different law: the vertex is the quadrature point; the surface is embedded and
well-shaped, with degenerate triangles refused rather than regularised; topology is fixed, which is
what licenses dropping the Gaussian term.

**The lumped area is barycentric and NOT mixed-Voronoi, and that is not an accident to be improved
here.** The reference states why: the mixed construction branches on whether a triangle is obtuse,
so the energy is only C^0 across that branch and the force is discontinuous. This kernel reproduces
the barycentric lumping exactly. A "more accurate" kernel would be a different law and would fail its
own gradient check.

### What is deliberately **not** ported: the spontaneous-curvature sector

`c0 != 0` adds `E_cross = -kappa c0 sum_i K_i . n_i` and `E_offset = (kappa c0^2/2) A_total`, and a
third gradient stage through the vertex normal. **It is refused, not silently ignored:**
`membrane_state` raises `ValueError` for `c0 != 0`. §14.3 records this as a named gap. The vertical's
own controls exercise only `c0 = 0` (`membrane.py`: *"Zero is the default and the only value this
vertical's controls exercise"*), and `assembly.recovered_tension_pn_per_um` **refuses** `c0 != 0`
outright, so the sector that G2's oracle can judge is exactly the sector ported.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| vertex positions `p` | µm | m | finite |
| face area `A_f`, vertex area `A_i` | µm² | m² | `> 0` strictly |
| curvature vector `K_i` | µm | m | finite |
| bending rigidity `kappa` | pN·µm | J | `>= 0` finite |
| tension `sigma` | pN/µm | N/m | `>= 0` finite |
| energy | pN·µm | J | `E_area >= 0`, `E_bend >= 0` |
| force | pN | N | finite |
| ULP figure | float32 ULP of the field's scale | dimensionless | `>= 0`; `inf` if the reference field is identically zero and the candidate is not |

Singular and boundary cases, each with the behaviour Aleph requires:

- **A degenerate triangle (`|N_f| = 0`).** The reference raises `DegenerateGeometryError`. The case
  builder rejects such a mesh before either side runs; the kernel does not regularise, and no
  epsilon is added to any denominator. A geometric failure that becomes a large finite force
  pointing in an arbitrary direction is strictly worse than stopping.
- **A vertex with zero lumped area.** Same: refused by the builder.
- **`c0 != 0`.** Refused by the builder (§4). Not approximated, not ignored.
- **An open surface.** Refused: the state builder requires every edge to be shared by exactly two
  faces. The Gaussian-term argument and the volume identity both depend on closure.
- **Empty mesh.** Refused. A gate that passes vacuously on no faces is the failure this whole track
  exists to prevent.
- **A negative cotangent (an obtuse corner).** **Not** singular and **not** refused — it is a
  legitimate branch of the law that a wrong kernel can clamp away. See I7.

Invariants that must hold, each with the test that asserts it:

- **I1.** A membrane law compared with itself is `0.00` ULP on every field, exactly —
  `test_a_membrane_law_compared_with_itself_is_exactly_zero_ulp`.
- **I2.** `E_bend` computed **by the kernel** is invariant under uniform scaling of the mesh, to
  float32 — `test_the_kernel_bending_energy_is_invariant_under_uniform_scaling`.
- **I3.** `E_bend` computed **by the kernel** approaches `8 pi kappa` as the icosphere refines, and
  does so monotonically — `test_the_kernel_bending_energy_approaches_the_sphere_limit`.
- **I4.** The kernel's assembled membrane force recovers the declared tension through the dilational
  virial, `sigma = -(sum_i f_i . x_i) / (2 A)`, at every one of the tension sweep's twelve
  configurations — `test_the_kernel_forces_recover_the_declared_tension_across_the_sweep_grid`.
- **I5.** `F = -grad E` for the kernel's own energy and its own force, by central difference —
  `test_the_kernel_force_is_minus_the_gradient_of_the_kernel_energy`.
- **I6.** Every scatter destination is allocated through the backend, never as a host ndarray —
  `test_every_membrane_scatter_destination_is_allocated_through_the_backend`.
- **I7.** The default state contains obtuse corners (negative cotangents) and both icosphere
  valences — `test_the_default_membrane_state_visits_both_cotangent_signs_and_both_valences`.
- **I8.** The parity report names the position-precision mode alongside the scatter mode, the device
  and the compute dtype — `test_the_report_names_the_position_precision_and_the_scatter_mode`.
- **I9.** The reference law is driven on **float64** positions and is measurably not a float32 law —
  `test_the_membrane_reference_is_driven_in_float64_positions`.

## 6. Source evidence class and known retractions

Nothing is inherited from another project, so there is no inherited evidence class. What is carried:

- **`ALEPH-PORT-1101`** — the reference law's own ledger entry, in this repository. Its controls
  (the central-difference force check in `tests/vertical/test_membrane.py`, the exact scale
  invariance I5) are the reason the reference can be treated as frozen.
- **`docs/results/2026-07-31-tension-sweep/README.md`** — Laplace recovery to a worst relative error
  of **1.566e-05** across σ = 10…60 pN/µm and two resolutions, twelve configurations, all converged.
  Read in this repository. **Two properties of that oracle are load-bearing and are used here:**
  the relative error *is* the dilational virial of the residual force, matched to 2.2e-16; and
  `max_residual_force_pn` **cannot** gate it, because the virial is a signed sum and twenty stopping
  points on one fixed mesh span 603× in |error| with five sign changes.
- **`docs/design/GPU_STATE_2026-07-31.md` and commit `636b0c8`** — `ORDERED` = 0.00 ULP /
  `ATOMIC` = 10.00 ULP on the A5000. `636b0c8` is itself a retraction of `278e6e5`, which reported a
  GPU-only backend defect that was a misuse of `WarpBackend.scatter_add` (a raw host ndarray as
  `dest`). Every scatter destination in this entry comes from `backend.zeros(...)`, and I6 is the
  control.
- **`ALEPH-PORT-3601` §14.7** — the float32-position finding, **raised to the PI in
  `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` and not taken here.**
- **No retraction was searched for in `/Users/sw1/ffn_cellsim`**, because nothing from that tree is
  used. No magnitude from that project appears in this entry or in the code it authorises.

## 7. Independent oracle or derivation

Four oracles, of decreasing strength, and the reason there are four is the lesson Lane A's nucleus
work produced tonight: **a quantity produced by heavy cancellation is not a parity target.** Lane A's
net internal force is a 15-order cancellation (376 pN → ~4e-13) where a single-ULP change moves the
result by 40%. For the membrane the exact identities are available, so they are used.

**O1 — the frozen NumPy law**, `aleph.vertical.membrane.helfrich_energy_and_forces`, in float64, on
the host, with no backend. Called, never re-implemented. It checks the kernel; it cannot check
itself, which is what O2–O4 are for.

**O2 — exact scale invariance of the discrete bending energy.** `E_bend(lambda x) = E_bend(x)`
*identically*, verified here in float64 to the last bit (measured: bit-identical at
lambda ∈ {0.5, 1, 2, 7.3}). `K_i` is degree one and `A_i` degree two in the positions, so
`|K|^2/A_i` is degree zero. **This is a strong anchor and a documented blind spot at the same time:
it is invariant to a wrong constant multiplying the whole term**, so a kernel that lost the
barycentric `1/3` passes it. That is asserted as a control (§9), not left implicit.

**O3 — the sphere limit, which is what closes O2's blind spot.** On a sphere of any radius,
`E_bend = (kappa/2)(2/R)^2 (4 pi R^2) = 8 pi kappa`, **radius-independent**, so a wrong
normalisation cannot hide inside a radius. Measured in float64 on the reference at radius 5 µm:

| icosphere level | vertices | `E_bend / (8 pi kappa)` |
|---|---|---|
| 1 | 42 | 0.931402 |
| 2 | 162 | 0.982400 |
| 3 | 642 | 0.995560 |

Converging monotonically to 1 with refinement, which is the discretisation error of the one-point
quadrature and is the reason this is an oracle with a stated tolerance rather than an identity.
A factor-of-3 or factor-of-4 normalisation error is 3.0 or 4.0 on this scale, four orders of
magnitude outside the discretisation error, so this oracle separates them without ambiguity.

**O4 — the dilational virial, which is the tension sweep's own derivation.** For the membrane's
total energy at `c0 = 0`,

```
E(lambda x) = E_bend(x) + lambda^2 sigma A(x)
=>  sum_i (grad E)_i . x_i = dE/dlambda |_{lambda=1} = 2 sigma A
=>  sigma = -(sum_i f_i . x_i) / (2 A)          f = -grad E
```

**exactly, at any configuration** — not only at equilibrium, and with no continuum limit anywhere.
Verified in float64 on the reference: `|sigma_recovered - sigma| / sigma <= 1.8e-16` at
σ ∈ {10, 30, 60} on levels 1, 2 and 3. This is the identity
`assembly.recovered_tension_pn_per_um` rests on, with the transmitted load replaced by the membrane's
own internal force — the two differ, at equilibrium, exactly by the residual virial, which is what
the sweep's 1.566e-05 *is*.

So the sweep reproduces through the kernels as I4: the same twelve `(level, sigma)` configurations,
σ read back from the kernel's own assembled force. §14.5 states precisely what that does and does
not establish relative to the published table.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_true_areal_kernel_agrees_with_the_frozen_numpy_law` | The float32 Warp areal kernel, on warp's CPU device with `ScatterMode.ORDERED`, agrees with `helfrich_energy_and_forces` on energy and on the assembled vertex force, within each field's declared ULP budget. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_true_bending_kernel_agrees_with_the_frozen_numpy_law` | The same for the cotangent-Laplacian bending term — **the first stencil kernel in this project**, with a cotangent weight per interior edge and a lumped vertex area. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_a_membrane_law_compared_with_itself_is_exactly_zero_ulp` | I1, the harness's zero point on these cases. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_kernel_bending_energy_is_invariant_under_uniform_scaling` | I2 / O2. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_kernel_bending_energy_approaches_the_sphere_limit` | I3 / O3, including monotone convergence across three levels. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_kernel_forces_recover_the_declared_tension_across_the_sweep_grid` | I4 / O4 — **the oracle that decides whether this port is correct**, at all twelve sweep configurations. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_kernel_force_is_minus_the_gradient_of_the_kernel_energy` | I5, by central difference on the kernel's own energy, with the step chosen from the measured float32 energy noise rather than assumed. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_every_membrane_scatter_destination_is_allocated_through_the_backend` | I6, by type. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_default_membrane_state_visits_both_cotangent_signs_and_both_valences` | I7 — the vacuity guard on the *state*, which is the limit G1 named and could not close. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_report_names_the_position_precision_and_the_scatter_mode` | I8. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_membrane_reference_is_driven_in_float64_positions` | I9 — the mechanical form of the hazard the proposal raises: a reference secretly given float32 input would agree with a float32 kernel and both would be wrong together. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_local_coordinates_recover_the_accuracy_a_world_offset_destroys` | `PositionPrecision.LOCAL_F32` recovers, on a mesh whose centre sits 50 µm from the world origin, the agreement that `GLOBAL_F32` loses there — the measured basis for the proposal's option B. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_float64_positions_beat_float32_positions_on_a_refined_mesh` | `PositionPrecision.POSITIONS_F64` measurably improves on `GLOBAL_F32` at level 3, where local coordinates cannot help — the measured basis for the proposal's option D. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_position_precision_is_a_parameter_and_not_a_constant` | All three modes are reachable through the public case builders and produce **three different** measurements. A parameter that could not change the answer would be decoration, and a PI decision under the proposal would change nothing. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_reference_splits_the_two_energies_exactly` | The split into two cases (`kappa = 0` for one, `sigma = 0` for the other) reproduces the combined law to 1e-14. If it did not, every figure below would be measuring the split rather than the kernel. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_the_declared_membrane_budgets_sit_inside_what_roundoff_explains` | No declared budget exceeds its own round-off explanation, so a red gate cannot be turned green by widening a number past what the mesh's conditioning could excuse. |
| Positive | `tests/runtime/test_membrane_law_parity.py::test_every_bending_mutant_is_caught_on_the_default_state` | All four bending mutants — one per stage of the four-stage chain — fail on one state. |

## 9. Deliberately failing negative control

Six wrong kernels ship in `aleph/runtime/law_kernels.py`, beside the right ones, compiled by the same
`load()` and launched down the same driver. All six are *plausible hand-translation errors*.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `::test_an_areal_kernel_that_dropped_the_half_in_its_area_is_caught` | `ArealVariant.WRONG_AREA_HALF` writes `A_f = |N_f|` where the law says `|N_f|/2`. **Its forces are bit-identical to the true kernel.** Caught on the energy field alone. |
| Negative (must fail) | `::test_an_areal_kernel_that_reused_one_corner_gradient_is_caught` | `ArealVariant.WRONG_CORNER_GRADIENT` computes `grad A_f` for corner 0 and uses it for all three — a plausible loop bug. **Its energy is bit-identical to the true kernel.** Caught on the force field alone. Together with the row above this is the argument for grading fields separately, in both directions. |
| Negative (must fail) | `::test_a_bending_kernel_that_dropped_the_half_in_the_cotangent_laplacian_is_caught` | `BendingVariant.WRONG_COTANGENT_HALF` drops the `1/2` in `K_i`, so `E_bend` is 4× and the force is 4×. The stencil's own constant. |
| Negative (must fail) | `::test_a_bending_kernel_that_forgot_the_barycentric_third_is_caught` | `BendingVariant.WRONG_VERTEX_AREA_LUMPING` uses `A_i = sum_{f>i} A_f` without the `1/3`. `E_bend` is 3× too small. |
| Negative (must fail) | `::test_a_bending_kernel_that_froze_the_cotangent_weights_is_caught` | `BendingVariant.WRONG_FROZEN_COTANGENT` omits the `grad cot` term from the force — treating the cotangent weights as constants under vertex motion. **The single most common real error in a cotangent-Laplacian force code**, and its energy is bit-identical to the true kernel. |
| Negative (must fail) | `::test_a_bending_kernel_that_lost_the_cotangent_sign_is_caught_on_an_obtuse_mesh` | `BendingVariant.WRONG_ABS_COTANGENT` takes `|cot|`, which is the mesh analogue of G1's missing unilateral clamp. |
| Negative (must fail) | `::test_a_zero_ulp_budget_fails_the_true_membrane_kernels` | The vacuity guard on the verdict: at `ulp_budget = 0.0` the *correct* kernels fail, because float32 storage is not bitwise float64. |
| Negative (must fail) | `::test_no_membrane_mutant_escapes_through_the_atomic_floor` | Neither the widened `ATOMIC` budget nor the scatter floor is an amnesty for any of the six. |
| Negative (must fail) | `::test_a_wrong_bending_kernel_breaks_the_tension_recovery` | **The vacuity guard on the oracle in §7.** A tension-recovery gate that passed whatever the bending kernel did would be grading the areal kernel alone and reporting it as a statement about both. `WRONG_COTANGENT_HALF` must move the recovered sigma. |

And the three blind spots the negative controls expose, recorded as controls because they are the
honest limits of the method:

| Control | Location | Asserts |
|---|---|---|
| Blind spot (must hold) | `::test_scale_invariance_is_blind_to_a_wrong_bending_normalisation` | `WRONG_VERTEX_AREA_LUMPING` **passes** O2's exact scale-invariance anchor, because a constant factor is invisible to a homogeneity argument. It is O3, the sphere limit, that kills it — measured `E_bend/(8 pi kappa) = 3` instead of `~0.98`. **An exact identity is not automatically a strong test; it is only as strong as what it is not invariant to.** |
| Blind spot (must hold) | `::test_the_frozen_cotangent_mutant_is_invisible_in_the_energy` | `WRONG_FROZEN_COTANGENT` is bit-identical to the true kernel on the energy field and is caught only by the force field and by I5. An energy-only gate certifies it. |
| Blind spot (must hold) | `::test_the_lost_cotangent_sign_is_invisible_on_an_all_acute_mesh` | On an icosphere, whose corners are all acute, `WRONG_ABS_COTANGENT` is bit-identical to the true kernel and the gate passes it at 0.00 ULP. **This is G1's all-taut blind spot in mesh form**, and it is why the default state is a *perturbed* mesh with a stated obtuse-corner count and I7 guards that count. |
| Blind spot (must hold) | `::test_a_pristine_sphere_is_not_a_valid_parity_state_for_the_bending_force` | On a pristine icosphere the **correct** bending kernel measures 77.75 ULP against the same budget it meets at 8.64 on the jittered mesh. A sphere has uniform curvature, so it is nearly a critical point of the bending energy, the per-vertex force is a near total cancellation, and the ULP denominator collapses with it. **That is Lane A's lesson — a quantity produced by heavy cancellation is not a parity target — arriving in the membrane**, and it says that the obvious state, the clean icosphere everyone reaches for, is the wrong one to grade a force on. |

## 10. Numerical and precision envelope

**Working precision.** Positions enter in the precision the case declares (below). Everything derived
from them — edge vectors, normals, areas, cotangents, the area gradient, the cotangent gradient, the
curvature vector, the coefficient fields and the assembled force — is **float32**, `ALEPH-DQ-107`'s
compute channel. Per-face and per-vertex **energies are widened to float64 inside the kernel** before
the product is formed, and the total is reduced by `Backend.sum`'s two-stage float64 accumulation, so
the energy carries one rounding per input rather than one per addition. The reference is float64 end
to end and I9 asserts it.

**Position precision is a declared parameter, not a buried constant.** This is the lane's response to
`ALEPH-PORT-3601` §14.7 and to the build plan's instruction that no lane may pick silently:

| `PositionPrecision` | What it does | Device cost |
|---|---|---|
| `GLOBAL_F32` | positions uploaded float32 in world coordinates, differenced in float32 | the status quo; what G1 used |
| `LOCAL_F32` | the mesh centroid subtracted **in float64 on the host**, positions uploaded float32 as offsets, differenced in float32 | none — one host subtraction per upload |
| `POSITIONS_F64` | positions uploaded float64; the **only** float64 arithmetic on the device is the subtraction that forms the edge vectors, whose result is rounded to float32 once | 2× position storage, one FP64 subtract per edge component |

**The design point that makes the parameter cheap.** Exactly **one** kernel in this entry reads
positions: `k_face_edges_*`, which forms the three opposite-edge vectors of each face. Every other
quantity in both laws — areas, normals, cotangents, all three gradient chains, the curvature vector —
is a function of those edge vectors and of nothing else. So the precision decision is localised to a
single ~20-line kernel that exists in two dtype variants, and a PI decision under
`docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` flips a keyword argument.

**Where the disagreement comes from, derived rather than assumed.** Writing `u = eps32/2` for the
unit round-off, each position is rounded once before anything happens, so an edge inherits both
roundings:

```
|d e_{f,k}| <~ u (|p_a| + |p_b|)
cancellation = max_{f,k} (|p_a| + |p_b|) / |e_{f,k}|      dimensionless, a property of the mesh
```

`cancellation` grows as the mesh refines — the edge halves per level while `|p|` does not — and grows
linearly with the distance of the body from the world origin. Both effects are measured in
`docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` §3 and reproduced in §13 here. The
per-field amplification is `cancellation` times a factor counting the differencings in the chain that
feeds that field, computed **from the host state by the case** so that it is a property of the data
and not a number chosen after seeing the answer; `roundoff_ulp_bound = (amp + 1)/2` converts it to
ULP exactly as `parity.py` already does. The declared budget is the measurement plus a margin, and
`LawFieldParity.budget_within_roundoff_explanation` reports whether it sits inside its own
justification — a budget wider than its round-off explanation can hide a defect and is flagged.

**One association-order difference from the reference, stated rather than discovered.** The reference
accumulates the gradient with several separate `np.add.at` calls, so each vertex sums its
contributions across faces in a fixed but term-by-term order. The kernel fuses a face's three
gradient sub-terms **in registers** before scattering one vector per corner, which is the point of a
kernel. In float32 this is a different summation order and it contributes to the measured ULP. It is
not a different law: both compute the same sum of the same terms.

**The atomic floor.** `ATOMIC_SCATTER_ULP_FLOOR = 10.0` (`parity.py`, measured on the A5000). The
assembled vertex forces and the lumped vertex area are scatter-assembled and take the floor in
`ATOMIC` mode; the energies are reduced by `Backend.sum`, which is ordered by construction, and do
not. §13's measurement shows the floor is **not** the binding term for this law either.

**Outside the envelope.** One device — and it is **warp's CPU device**, on a Mac whose warp build
reports *"CUDA not enabled in this build"* and `wp.get_devices() == ['cpu']`. `ATOMIC` has never been
measured here, and `measure_scatter_determinism` already refuses to let a CPU zero be read as
evidence about CUDA. This entry inherits that caveat verbatim rather than restating it weaker.

## 11. Production-backend residency and transfer

The kernel driver is production code and runs on the device. Per case:

| Array | Where | Precision | Transfer |
|---|---|---|---|
| vertex positions | device | float32 or float64, per `PositionPrecision` | uploaded once per case by `backend.array` |
| triangle indices | device | int32 | uploaded once per case |
| edge vectors, normals, areas, cotangents, both gradient tables | device | float32 | never leave the device |
| per-face / per-vertex energy | device | float64 | never leaves the device before `backend.sum` reduces it |
| lumped vertex area, curvature vector, assembled vertex force | device, **allocated by `backend.zeros`** | float32 | read to host once, by the harness, for the comparison |

**The rule inherited from `ALEPH-PORT-3601` §11 and kept by I6.** `WarpBackend.scatter_add` takes a
`dest` allocated *through the backend*; a raw host ndarray works on warp's CPU device and fails on
CUDA. `278e6e5` reported that as a GPU-only backend defect and `636b0c8` retracted it — the API was
right and the caller was wrong. This entry has **three** scatter destinations (lumped vertex area,
curvature vector, assembled force), which is three chances to make that mistake, and the control
asserts by type rather than by comment.

The harness itself is host-side and is not production physics. It reads each graded field to the host
exactly once per case.

## 12. Comments and docstrings to discard

Nothing was read from another project, so there is no source prose to discard. The list is kept as a
positive statement of what must **not** appear in the modules this entry authorises, with
`tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package` as the
mechanical half:

- no provider repository name, module path, kernel name, gate name or branch name;
- no absolute path into another project;
- no provider datum as a literal;
- no claim that a number was measured on hardware this lane did not run on. Every ULP figure in the
  package's docstrings is either measured here on warp's CPU device and says so, or is cited to
  `docs/design/GPU_STATE_2026-07-31.md` / `636b0c8` and says that.

What replaces them: the derivations in §4 and §10 restated in the modules' own docstrings, and the
measured numbers from this repository's own run records.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Measured 2026-07-31 on **warp's CPU device** (this machine's warp build reports *"CUDA not enabled in this build"* and `wp.get_devices() == ['cpu']`), warp 1.15.0, float32 compute, `ScatterMode.ORDERED`, `PositionPrecision.GLOBAL_F32` unless stated. `tests/runtime/test_membrane_law_parity.py`: **33 passed**. `tests/runtime` as a whole: **302 passed, 3 skipped** (was 269 + 3 before this entry). `tests/vertical`: **1,358 passed**, unchanged — the frozen reference was not edited. Mutation study: **7 of 7 killed** (§13a). **The tension sweep reproduces through the kernels at a worst relative error of 1.5610e-05 against the published 1.5663e-05, twelve of twelve converged.** The status stays `PROPOSED` because no CUDA device was touched and no human has reviewed it. |
| Reviewer | **Agent-proposed. Unratified.** Written by a subagent of session `46143f30`; no PI review, no `decided_by` field, and none may be added by an agent. |
| Rollback | Revert the appended blocks of `aleph/runtime/law_kernels.py` and `aleph/runtime/law_cases.py` and delete `tests/runtime/test_membrane_law_parity.py`. `ALEPH-PORT-3601`'s tether path is untouched by construction and keeps working. Nothing under `aleph/vertical/` changes in either direction. |

### 13.1 The parity measurement

Default state: a level-2 icosphere at radius 5 µm, jittered, **162 vertices / 320 faces, 6 obtuse
corners out of 960, valences {5, 6}**.

**The areal response.**

| kernel | `energy_pn_um` | `forces_pn` | verdict |
|---|---|---|---|
| `TRUE` | **0.01** | **5.03** | pass |
| `WRONG_AREA_HALF` | **8,388,608.03** | 5.03 *bit-identical* | **fail — energy only** |
| `WRONG_CORNER_GRADIENT` | 0.01 *bit-identical* | **11,448,867.64** | **fail — force only** |
| declared budget | 4 | 12 | — |
| round-off explains | 11.09 | 16.38 | — |

**The Helfrich bending term.**

| kernel | `energy_pn_um` | `forces_pn` | verdict |
|---|---|---|---|
| `TRUE` | **0.09** | **8.64** | pass |
| `WRONG_ABS_COTANGENT` | 22,350.64 | 333,525.89 | fail — both |
| `WRONG_COTANGENT_HALF` | 25,165,824.36 | 8,050,802.50 | fail — both |
| `WRONG_VERTEX_AREA_LUMPING` | 5,592,405.31 | 5,527,183.77 | fail — both |
| `WRONG_FROZEN_COTANGENT` | 0.09 *bit-identical* | **1,377,764.29** | **fail — force only** |
| declared budget | 8 | 24 | — |
| round-off explains | 26.97 | 64.04 | — |

Four things to read off these two tables.

1. **Three of the six mutants are invisible in one of the two fields, and they are invisible in
   different ones.** `WRONG_AREA_HALF` is energy-only; `WRONG_CORNER_GRADIENT` and
   `WRONG_FROZEN_COTANGENT` are force-only. A gate reporting one aggregated number certifies at
   least one of them, and a gate grading only the forces — which is what most parity gates grade,
   because forces are what the integrator consumes — certifies the first.
2. **Every mutant clears its budget by three to six orders of magnitude**, so the 2.4–2.8× margin
   between measurement and declared budget is nowhere near the detection threshold.
3. **Every declared budget sits below what round-off explains**, asserted by
   `test_the_declared_membrane_budgets_sit_inside_what_roundoff_explains`.
4. On the **all-acute** pristine sphere the true kernel measures 0.59 / 77.75 and
   `WRONG_ABS_COTANGENT` measures **exactly the same two numbers** — bit-identical, gate passes.
   That is §9's blind spot and §14.4's limit, observed rather than argued.

### 13.2 Position precision — the measured basis for the proposal this lane routed upward

Assembled **bending force** ULP, true kernel, by mode and by mesh:

| mesh | body at the origin | body at (50, 0, 0) µm |
|---|---|---|
| L1 (42v) | `GLOBAL_F32` 5.94 / `LOCAL_F32` 5.43 / `POSITIONS_F64` **2.58** | 35.32 / **5.43** / **2.58** |
| L2 (162v) | 8.64 / 8.13 / **1.61** | 62.52 / **8.13** / **1.61** |
| L3 (642v) | 21.86 / 23.56 / **1.68** | 139.05 / **23.56** / **1.68** |

Assembled **areal force** ULP, same layout: L1 1.77 / 1.58 / 0.99 and 11.13 / 1.58 / 0.99;
L2 5.03 / 4.28 / 0.95 and 30.76 / 4.28 / 0.95; L3 13.06 / 12.75 / 1.90 and 85.90 / 12.75 / 1.90.

**This confirms, in a kernel, the two-term structure the proposal argues from.** `LOCAL_F32`
recovers the at-origin figure from the displaced one *exactly* — 139.05 → 23.56 at level 3 — and
buys essentially nothing when the body is already centred. `POSITIONS_F64` is the only mode whose
figure does **not** grow with refinement: 2.58 → 1.61 → 1.68 against 5.94 → 8.64 → 21.86.

The membrane's float32 penalty is far milder than the tether's ~100 ULP, and the reason is
structural rather than lucky: a mesh stencil differences *neighbouring* vertices, whose separation
is a fraction of the radius, while a tether differences two points 0.03 µm apart on surfaces 5 µm
from the origin. **The cancellation ratio is the whole story in both cases**, which is the
proposal's point.

### 13.3 The exact anchors

| oracle | measured through the kernel | reference / expected |
|---|---|---|
| scale invariance of `E_bend` | **bit-identical** at λ ∈ {0.25, 0.5, 1, 2, 4} | exact, degree zero |
| sphere limit `E_bend / 8 pi kappa`, L1 / L2 / L3 | **0.931402 / 0.982400 / 0.995560** | same to six digits in float64; → 1 |
| the same, `WRONG_VERTEX_AREA_LUMPING` | **0.327467** | exactly `0.982400 / 3` |
| `F = -grad E`, central difference | error ratios **4.06 / 4.05 / 4.00** for h = 0.8→0.1 | order 2 |
| `F = -grad E`, `WRONG_FROZEN_COTANGENT` | **38% off** the numeric derivative | caught |

The scale-invariance row and the sphere-limit row are the pair that matters:
`WRONG_VERTEX_AREA_LUMPING` **passes** the first bit-identically and lands on exactly a third of the
second. **An exact identity is only as strong as what it is not invariant to**, and that is asserted
as a control rather than left as a remark.

The `F = -grad E` order-2 window is measured, not assumed: the truncation term `~h^2` and the
float32 energy's own round-off `~eps/h` cross at `h ≈ 0.1` on this mesh, so the ratios are clean
above it and stall below it (1.84e-04 → 4.53e-05 → 1.12e-05 → 2.79e-06, then 1.15e-06 at h = 0.05).
**A demonstration of order 2 below the noise floor would have been a demonstration of nothing**, and
the step is chosen from the measurement rather than from habit.

### 13.4 The oracle that decides the port: the tension sweep, three ways

`docs/results/2026-07-31-tension-sweep/README.md`: σ recovered to a worst relative error of
**1.566e-05** across σ = 10…60 pN/µm and two resolutions, twelve configurations, all converged.

**(a) The identity, with no relaxation at all.** `sigma = -(sum_i f_i . x_i)/(2A)` holds at any
configuration, so the twelve `(level, sigma)` configurations were run directly through the kernels on
two mesh families — the pristine icosphere and a jittered non-sphere. **Worst relative error over
all 24: 3.5253e-08.** This is `test_the_kernel_forces_recover_the_declared_tension_across_the_sweep_grid`.

**(b) The identity on the sweep's own relaxed configurations.** The twelve configurations were
relaxed by the unmodified NumPy engine — reproducing the published table row for row, worst
**1.5663e-05**, which is the check that this is the same measurement and not a different one that
happens to agree — and σ was then read back through the kernels' force:

| position mode | worst relative error over the twelve relaxed configurations |
|---|---|
| `GLOBAL_F32` | **1.1962e-07** |
| `LOCAL_F32` | 1.1962e-07 |
| `POSITIONS_F64` | **6.0022e-08** |

These are **two orders of magnitude tighter than the sweep's own number, and that is not a claim
that the kernel beats the relaxation.** The sweep's 1.566e-05 is the dilational virial of the
*residual* force — a statement about where the run stopped — while (b) is the identity on the
membrane's internal force, which carries no residual. The two measure different things and the
difference between them *is* the stopping error.

**(c) The literal reproduction: the relaxation driven by the kernels.** `internal_forces` supplied
entirely by the G2 kernels (`POSITIONS_F64`), `relax_to_equilibrium` otherwise unmodified,
`max_steps = 40,000`:

| descent energy | converged | worst relative error |
|---|---|---|
| float64 (the frozen NumPy energy) | **12 / 12** | **1.5610e-05** |
| the kernels' own float32-term energy | **0 / 12** | 5.1755e-03 |

**The first row is the port's acceptance: 1.5610e-05 against a published 1.5663e-05, twelve of
twelve converged, with every membrane force on the mesh produced by a Warp kernel.**

**The second row is a finding and §14.6 states it as one.** Isolated by holding one channel fixed at
a time on `L1_sigma10`:

| membrane force | descent energy | accepted | converged | max\|F\| [pN] | recovered σ | rel |
|---|---|---|---|---|---|---|
| kernel | float64 | 2,725 | **True** | 9.27e-04 | 9.99994105 | 5.90e-06 |
| NumPy | float64 | 2,689 | True | 9.87e-04 | 9.99995176 | 4.82e-06 |
| kernel | float32 terms | 573 | **False** | 2.78e-01 | 10.00056190 | 5.62e-05 |
| NumPy | float32 terms | 527 | **False** | 3.30e-01 | 10.00728332 | 7.28e-04 |

**The force channel is not what stalls it; the energy channel is.** A NumPy force with a float32
energy stalls just as hard. `MonotonePotentialDescent` compares two energies, and at
`E ≈ 2,929 pN·µm` one float32 ULP is `≈ 3.5e-04 pN·µm`; once the per-step decrease falls below that
the controller can no longer tell a descent from noise, the step collapses to `~5e-13`, and the run
stops two orders of magnitude short of the residual the same mesh reaches in float64.

### 13a. Mutation study — 7 planted defects, 7 killed

Run with `PYTHONDONTWRITEBYTECODE=1` (a same-length mutant with a sub-second edit/revert leaves a
`.pyc` CPython considers valid, which has already contaminated one study in this repository).
Baseline green, each mutant applied alone, reverted before the next.

| # | Planted defect | Killed by |
|---|---|---|
| M1 | the bending driver ignores the requested variant and always launches the true stages | 12 tests |
| M2 | the areal driver launches `areal_true` whatever was asked | 3 tests |
| M3 | a bending scatter destination allocated as a host ndarray instead of `backend.zeros` | 25 tests |
| M4 | `LOCAL_F32` stops subtracting the reference point (silently becomes `GLOBAL_F32`) | 2 tests |
| M5 | `POSITIONS_F64` silently uploads float32 and launches the float32 edge kernel | 2 tests |
| M6 | the bending gradient is stored without the negation, so the "force" is `+grad E` | 5 tests |
| M7 | **the frozen reference is handed float32-rounded positions** | 5 tests |

**M7 is the one worth a sentence.** It is the exact hazard
`docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` §4 names — a reference quietly
given the same float32 input as the kernel, so both lose the same digits and the gate reports a
small healthy number. It is killed by
`test_the_membrane_reference_is_driven_in_float64_positions`, which is the mechanical form of that
hazard at case level. M4 and M5 exist because a declared parameter that no test can see the effect
of is decoration.

## 14. Honest limits

What this entry does **not** establish:

1. **No CUDA device was touched.** This machine's warp build reports *"CUDA not enabled in this
   build"* and `wp.get_devices()` returns `['cpu']`; `~/.aleph_data/gpu_authorization.json` is
   expired and names another host, and an agent may never write it. Every number here is from
   **warp's CPU device**, where a launch is a serial loop. It is evidence that the kernels' algebra
   matches the reference and that the identities hold; it is `UNVERIFIED` as evidence about CUDA. In
   particular the `ATOMIC` rows are reproducible here by accident of the schedule, exactly as
   `measure_scatter_determinism`'s caveat already says.
2. **The `c0 != 0` sector is not ported.** The spontaneous-curvature cross and offset terms, and the
   third gradient stage through the vertex normal, are absent. `membrane_state` refuses `c0 != 0`
   rather than ignoring it. Nothing here is evidence about a membrane with spontaneous curvature,
   and the dilational-virial oracle would not be exact for one anyway
   (`assembly.recovered_tension_pn_per_um` refuses it for the same reason).
3. **One mesh family.** Every state is an icosphere, jittered or not. `PLAN.md` §2.5 measured that a
   Fibonacci hull at the same vertex count needed **189,419** relaxation steps against 3,426 — a
   factor of 55 from tessellation alone — so a tessellation this entry has not run is a
   tessellation this entry says nothing about. The kernels contain no assumption about valence
   beyond what the scatter maps carry, and the default state does visit both icosphere valences, but
   that is an argument and not a measurement.
4. **The negative controls are as strong as the state they run on, and two of them are provably
   blind somewhere.** `WRONG_ABS_COTANGENT` is bit-identical to the true kernel on any all-acute
   mesh; `WRONG_VERTEX_AREA_LUMPING` passes the exact scale-invariance identity. Both are asserted
   by controls. **This method cannot tell you that a case's inputs visit every branch of a law** —
   choosing the state remains a human judgement, and the improvement over G1 is that the branch
   counts are now *fields of the state* that a test can assert on, not that the judgement was
   removed.
5. **The declared budgets are engineering numbers and they belong to one mesh.** They are
   measurements at the default state plus a margin. §10's cancellation grows with refinement, so a
   case on another mesh must declare its own budget; `DEFAULT_*_ULP_BUDGET` says so on the constant.
   No energy, force or tension in this entry may be reported as a property of a cell: the evidence
   class is `STRUCTURAL`, and `ANALYTIC_ORACLE` for the Laplace recovery, which is an identity of
   the model and not a measurement of anything.
6. **A float32 energy channel is not usable as a monotone-descent criterion, and this entry does not
   fix it.** §13.4(c) measures it: 0 of 12 configurations converge and the worst recovered-σ error
   degrades 331× when the descent reads the kernels' own energy. The kernels already widen each
   per-face and per-vertex energy to float64 before the product and reduce with `Backend.sum`'s
   float64 accumulator; what remains float32 is the *term* — the face area, the curvature vector —
   and that is `ALEPH-DQ-107`'s ratified compute channel. **Widening it is a decision above this
   lane**, of the same family as the position-precision question, and it is raised in
   `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` §8 rather than taken here. The
   alternative — giving the descent predicate a relative tolerance instead of an exact comparison —
   is a change to `aleph/vertical/relax.py`, which this lane does not own.
7. **The energy is compared as a scalar total.** A kernel wrong on two faces in opposite directions
   would cancel and pass. A per-face and per-vertex energy comparison would catch it and is not
   implemented; inherited from `ALEPH-PORT-3601` §14.6 and named again rather than assumed fixed.
8. **The association order differs from the reference by construction** (§10): the kernel fuses a
   face's three gradient sub-terms in registers before scattering. That is a float32 difference, it
   is inside the measured budget, and it is not a different law — but it means the parity figure is
   not attributable to position storage alone.
9. **`POSITIONS_F64` is option D of the proposal, not option A.** Positions are float64 and the
   differencing is float64; everything after it is float32. A float64-end-to-end kernel is **not**
   implemented and no figure here bounds what it would buy.
10. **Nothing here measures throughput.** The FP64 cost of `POSITIONS_F64`, the cost of three
    scatters versus one, and the cost of the four-stage chain against a fused kernel are all
    unmeasured, on any device. The choice between the position modes is presented on accuracy alone,
    which is half of the decision the PI is being asked to make.
