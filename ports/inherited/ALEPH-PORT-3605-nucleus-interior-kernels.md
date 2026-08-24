# ALEPH-PORT-3605 — the nucleus interior on kernels, and the self-force term that is not optional

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3605` |
| Lane | `46143f30` Lane G4 (CUDA track) |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | `ALEPH-PORT-3501` landed the nucleus's interior in NumPy and the PI redirected the project to CUDA-only at 21:44. That entry is now the **frozen parity oracle** for this one. Two of its properties make the port harder than G1–G3's and are the reason it is a task of its own: the volume law's geometry **belongs to another module** and the package deliberately has exactly one implementation of it, and the owner applies a **live-centroid gradient correction** worth 2.978 pN of spurious self-force if it is dropped. |

---

## 1. Aleph API

```python
from aleph.runtime.law_kernels import (
    ChromatinVariant,           # TRUE | WRONG_ENERGY_HALF | WRONG_UNWEIGHTED
                                #      | WRONG_UNNORMALISED_DIRECTION
    SelfForceVariant,           # TRUE | WRONG_NO_CORRECTION | WRONG_SCALAR_MEAN
    VolumeVariant,              # TRUE | WRONG_AREA_TIMES_RADIUS | WRONG_GRADIENT_CYCLE
                                #      | WRONG_ENERGY_HALF | WRONG_FORCE_SIGN
    VOLUME_GEOMETRY_KERNELS,    # precision -> the ONE device volume-geometry kernel name
    VOLUME_STAGE_KERNELS,       # variant -> which stage kernel it swaps
    VOLUME_TRUE_STAGES,         # the default kernel of each volume stage, in launch order
    CHROMATIN_STAGE_KERNELS,
    CHROMATIN_TRUE_STAGES,
    SELF_FORCE_KERNELS,
)

from aleph.runtime.law_cases import (
    DEFAULT_CHROMATIN_ULP_BUDGET,
    DEFAULT_VOLUME_ULP_BUDGET,
    NUCLEUS_BULK_MODULUS_PN_PER_UM2,
    NUCLEUS_SHEAR_MODULUS_PN_PER_UM2,
    chromatin_law_case,             # the fixed-anchor interior network as a LawCase
    chromatin_owner_law_case,       # ... with the owner's live-centroid correction applied
    nucleus_kernel_interior,        # the kernels' own energies and forces, no reference involved
    nucleus_state,                  # a reproducible nucleus, with its branch counts as fields
    volume_geometry_kernel_name,    # the shared-kernel accessor both volume paths must use
    volume_kernel_volume_um3,       # the kernel's OWN enclosed volume, for the exact oracles
    volume_law_case,                # the nucleoplasm volume constraint as a LawCase
)
```

Nothing outside this list is covered. In particular this entry does **not** authorise any change to
`aleph/runtime/parity.py`, `aleph/runtime/backend.py`, `aleph/runtime/warp_kernels.py`, or to any
module under `aleph/vertical/` — the NumPy laws are the frozen reference and are read, never edited.
It also does not authorise regenerating `ports/ledger/INDEX.md`, which is already red from other
lanes (§14.8).

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read by this lane.** No file under that tree was opened. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

The law being ported is **Aleph's own**, in `aleph/vertical/nucleus_interior.py` and
`aleph/vertical/nucleus.py`, written by Lane A of this same session and recorded in
`ALEPH-PORT-3501`. That entry did read the provider tree and records what it took (§2, §6 there); the
present entry inherits its conclusions by citation and adds nothing from that tree. What this entry
ports is Aleph NumPy → Aleph Warp, which is a translation of code this repository owns.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported from another project.** The kernel is determined by three
Aleph surfaces: `aleph.vertical.nucleus_interior`'s two energy functions fix the law bit for bit,
`aleph.runtime.backend.Backend` fixes what a kernel driver may call, and `ALEPH-PORT-3601`'s
`LawCase` fixes the shape of the gate. A second implementation would need no help.

### 3.1 The one place a second implementation is a live hazard, and it is inside Aleph

`ALEPH-PORT-3501` §3.1 records that the first version of `nucleus_interior.py` carried **its own copy
of the enclosed-volume law**, that the copy disagreed with `aleph.vertical.membrane`'s by up to
`2.0e-13 µm³` on the same mesh, and that the fix was to delete it and import — with
`test_the_volume_law_has_exactly_one_implementation_in_the_package` asserting **function-object
identity**, `nucleus_interior.enclosed_volume_um3 is membrane.enclosed_volume_um3`, rather than
equality of output.

A port to kernels is exactly the occasion on which that fix is undone by accident, because the
device has no volume kernel to import: G2 ported the membrane's *areal* and *bending* terms and left
`enclosed_volume_um3` / `volume_gradient` on the host. So this entry introduces the device volume
geometry, and it introduces **one** of it:

- the kernels `k_volume_face_f32` / `k_volume_face_f64` are declared in the **membrane block** of
  `law_kernels.py` and are labelled as `ALEPH-PORT-1101`'s law, not the nucleus's;
- `VOLUME_GEOMETRY_KERNELS` maps a `PositionPrecision` to the one name, and
  `law_cases.volume_geometry_kernel_name(precision)` is the single accessor;
- the nucleus's volume driver resolves its geometry kernel through that accessor and through nothing
  else, and **a control records the kernel objects a candidate evaluation actually launches** and
  requires the shared object to be among them (§8, `test_the_device_volume_law_has_exactly_one_kernel_and_the_nucleus_launches_it`).

So the device statement is the same statement as the NumPy one, in the form the device admits:
*the object identity of the thing that ran*, not the equality of two numbers that agree today. §9.1
M1 is the mutant that gives the nucleus its own duplicate volume kernel, and it is the reason this
control exists rather than a comment.

**What is NOT claimed:** the kernel is not the same *implementation* as the NumPy law — it is a
port, and the whole point of the parity gate is that this is a second implementation being checked
against the first. The property preserved is narrower and is the one that was actually lost once:
there is one volume law **per layer**, so a future membrane-side volume case and this nucleus case
cannot drift apart on the device the way two NumPy copies drifted apart on the host.

### 3.2 What the constitutive split means on a device

`nucleus_interior.volume_energy_and_forces`'s own docstring states the division: the geometry is
`membrane`'s (`V` and `dV/dx`), and *"what this function adds is the CONSTITUTIVE part — the strain
measure, the energy and the chain rule — and nothing else."* The kernel chain keeps that division as
a **stage boundary**:

| stage | kernel | what it owns |
|---|---|---|
| geometry | `k_volume_face_*` | `p₀·(p₁×p₂)` per face and `(p₁×p₂)/6` etc. per corner — `ALEPH-PORT-1101` |
| reduction | `Backend.scatter_add` into a length-1 float64 slot | `Σ_f` — the global reduction §11 predicted would cost a synchronisation |
| constitutive | `k_volume_constitutive_*` | `ε_V`, `E_vol`, `dE/dV` — `ALEPH-PORT-3501` §4.2, and **nothing geometric** |
| chain rule | `k_volume_force` + `Backend.scatter_add` | `−(dE/dV)·∂V/∂x` assembled per vertex |

That is not decoration. It is what lets a mutant be attributed: `WRONG_ENERGY_HALF` and
`WRONG_FORCE_SIGN` live in the constitutive stage and are exactly `ALEPH-PORT-3501` §9.1's M1 and M5;
`WRONG_AREA_TIMES_RADIUS` and `WRONG_GRADIENT_CYCLE` live in the geometry stage; and the two are
graded in different channels.

## 4. Physical or mathematical law represented

Three laws, and the third is the one this task exists for.

### 4(a) The nucleoplasm volume constraint

For a closed, oriented, outward-wound triangulated surface (`ALEPH-PORT-3501` §4.1, `ALEPH-PORT-1101`):

```
V        = (1/6) Σ_f p₀·(p₁ × p₂)                                   [µm³]
∂V/∂p₀   = (p₁ × p₂)/6 ,  ∂V/∂p₁ = (p₂ × p₀)/6 ,  ∂V/∂p₂ = (p₀ × p₁)/6
ε_V      = V/V_ref − 1                                              [1]
E_vol    = (K_V/2) · V_ref · ε_V²                                   [pN·µm]
dE/dV    = K_V ε_V                                                  [pN/µm²]
f_a      = −(dE/dV) · ∂V/∂x_a                                       [pN]
```

`K_V = 0` gives **exactly** zero energy and exactly zero force — the inert model is a point in the
same family and not a separate code path. In the reference that is an early `return`; in the kernel
it needs no branch at all, because `K_V ε_V` is exactly `0.0` and `0.0 · finite` is exactly `0.0`.
That is a stronger statement than the reference's and it is asserted with `==` (§8).

**Three exact identities this form has, and each becomes an oracle:**

- **`det(M)` scaling.** Under any linear map about the centroid, `(Ma × Mb)·Mc = det(M)(a×b)·c` per
  face, so `V` scales by exactly `det(M)`. This is the oracle **with discriminating power**.
- **`s³` uniform scaling.** A special case, and `ALEPH-PORT-3501` §9.1 M10 is the measurement that it
  is *only* a special case: an area × mean-radius impostor passes it at every scale factor to twelve
  digits, because scaling about the centroid cannot distinguish the enclosed volume from **any**
  translation-invariant degree-three-homogeneous formula.
- **Reflection.** A reflection through a coordinate plane is an isometry of every triangle, so it
  preserves every face area **bit-identically** while negating `V` **bit-identically** — negation and
  the products of negated values are exact in floating point. So `V(reflected) == −V(original)` is
  testable with `==`, and the membrane's areal energy over the same surface is testable with `==`
  too, which is `ALEPH-PORT-3501` §7 oracle 2 arriving on a kernel.

### 4(b) The chromatin network

One radial strand per envelope vertex, anchored at the nuclear centroid, rest length `r₀`
(`ALEPH-PORT-3501` §4.3):

```
d_v      = x_v − c ,   r_v = |d_v|                                  [µm]
ε_v      = (r_v − r₀)/r₀                                            [1]
E_chr    = (G/2) Σ_v w_v ε_v²                                       [pN·µm]
f_v      = −(G w_v ε_v / r₀) · d_v/r_v                              [pN]
```

`w_v` is the vertex's share of the **reference** surface area and is held fixed: recomputing it from
the live positions drops a whole term of the gradient, and `ALEPH-PORT-3501` §10.3 records the
error curve of that mistake as a **plateau at 2.5942e-03** rather than a slope. The weights arrive
as a device array and are never recomputed inside a kernel.

The centroid is a **parameter**, for the same reason: a finite-difference probe must be able to hold
it still. `r_v = 0` makes the strand direction undefined and the reference **raises**. A kernel
cannot raise, so the device path refuses on the host in float64 before the launch **and** flags a
per-site `violation` the driver reads after it, exactly as `ALEPH-PORT-3604` §5 decided for the
contact. It never regularises.

### 4(c) The live-centroid correction — the highest-risk item in this port

`NucleusOwner.chromatin_term` uses the **live** centroid `c = (1/N) Σ_k x_k` and owes the gradient a
correction for it. With `E = Σ_v φ(|x_v − c|)`:

```
dE/dx_k = φ'(r_k) u_k + (dc/dx_k)ᵀ (dE/dc)
        = φ'(r_k) u_k − (1/N) Σ_v φ'(r_v) u_v
```

since `dc/dx_k = (1/N) I` and `dE/dc = −Σ_v φ'(r_v) u_v`. The first term is what the fixed-anchor
function returns, so the whole correction is

```
f = f_fixed − mean(f_fixed)                                          [pN]
```

— subtract the **column** mean, which is the projection that leaves the body with no net internal
force. This is a Newton's-third-law term. `ALEPH-PORT-3501` §4.5 measures what dropping it costs:
the chromatin term alone exerts **2.9777 pN** of net force on itself, against a peak single-vertex
chromatin force of `3.4230` pN — **87% of the peak**. It is not a rounding term and it does not
average away over a trajectory; it is a systematic drift of the whole body.

On the device the mean is a **column reduction of a force array inside a force evaluation**, which
is the second global reduction in this port. It is taken with `Backend.sum(f, axis=0)`, whose
`sum_axis0` kernel accumulates in float64 and is ordered by construction, and the subtraction is one
kernel. `WRONG_NO_CORRECTION` is `ALEPH-PORT-3501` §9.1's M9 as a kernel; `WRONG_SCALAR_MEAN`
subtracts `mean(f)` over *all* entries instead of per column, which is the one-character version of
the same mistake (`f.mean()` for `f.mean(axis=0)`) and is the more plausible of the two.

### 4(d) What is NOT ported here, stated so the gap is named

- **The lamina's two surface terms** (`lamina_areal_energy_and_forces`, `curvature_energy_and_forces`)
  stay on the host. G4 is the *interior*. So no fully-kernel `NucleusOwner` exists after this entry,
  and every assembled-owner figure below is a **hybrid**: interior from kernels, lamina from NumPy.
  Labelled as such everywhere it appears.
- **`reduced_coordinates`** — declared by the contract, unimplemented in the reference, and this
  entry does not invent it.
- **The drag / mobility** — a rate response, not a potential, and no kernel here consumes it.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| vertex positions, centroid | µm | m | finite; `(N, 3)` |
| `enclosed_volume` `V`, `V_ref` | µm³ | m³ | `V_ref > 0`; `V > 0` for an outward-wound closed surface |
| `∂V/∂x` | µm² | m² | finite |
| `bulk_modulus_pn_per_um2` (`K_V`) | pN/µm² | Pa | `≥ 0`; `0` is the inert model |
| `shear_modulus_pn_per_um2` (`G`) | pN/µm² | Pa | `≥ 0`; `0` is the inert model |
| `reference_radius_um` (`r₀`) | µm | m | `> 0` |
| `reference_vertex_area_um2` (`w_v`) | µm² | m² | `≥ 0`, one per vertex |
| energies | pN·µm | J | finite |
| forces | pN | N | finite |
| ULP figure | float32 ULP of the field's own scale | 1 | `≥ 0`; `inf` if the reference field is identically zero and the candidate is not |

Singular and boundary cases, each with the behaviour Aleph requires:

- **A vertex coincident with the centroid (`r_v = 0`).** The reference raises `ValueError`. The
  device path **refuses twice and never regularises**: a host float64 predicate before the launch,
  and the kernel's own `violation` flag afterwards. The second is not redundant — the host predicate
  is float64 and the kernel's is not.
- **`K_V = 0` / `G = 0`.** Exactly zero energy and exactly zero force, on **both** sides, asserted
  with `==` and with a parity budget of `0.0`.
- **A degenerate face.** `V` and `∂V/∂x` are polynomials and are finite on one; the volume law does
  not refuse it and neither does the kernel. The state builder refuses one, because a zero-area face
  makes the chromatin weights and every areal oracle meaningless.
- **A mesh that is not closed.** Refused by the state builder. The divergence-theorem form is
  meaningless on an open surface and a gate on one would be measuring nothing.
- **The rest configuration.** *Not* a singular case of the law and *is* a singular case of the
  **gate**: both interior energies are zero there, so `ulp_at_scale` divides by a reference scale of
  exactly zero and reports `inf` for a correct kernel. §10.5, and it has a control.

Invariants that must hold, each with the test that asserts it:

- **I1.** The device volume law has exactly one kernel per dtype in the package, and the nucleus's
  volume driver launches **that object** —
  `test_the_device_volume_law_has_exactly_one_kernel_and_the_nucleus_launches_it`.
- **I2.** `V_kernel(Mx) = det(M) · V_kernel(x)`, and the impostor kernel is required to **fail** it
  inside the same control — `test_a_linear_map_scales_the_kernel_volume_by_its_determinant`.
- **I3.** `V_kernel(s·x) = s³ · V_kernel(x)`, and the impostor kernel **passes** it — asserted, so
  the narrowness of the cube law is a measurement here too —
  `test_the_kernel_volume_falls_as_the_cube_and_the_impostor_passes_that_too`.
- **I4.** A reflection negates the kernel volume **bit-identically** and leaves the membrane areal
  kernel's energy **bit-identical** — `test_a_reflection_negates_the_kernel_volume_exactly`.
- **I5.** `F = −∇E` at observed order 2 for both interior terms, through the kernels, with the step
  read off a measured curve — `test_the_kernel_volume_force_is_minus_the_gradient_of_the_kernel_energy`.
- **I6.** `K_V = 0` and `G = 0` give exactly zero on both sides, at a declared budget of `0.0` —
  `test_an_inert_interior_is_exactly_zero_through_the_kernels`.
- **I7.** The chromatin response is mesh-resolution independent because the weights are areas —
  `test_the_kernel_chromatin_response_is_mesh_resolution_independent`.
- **I8.** The assembled body exerts no net force on itself once the correction is applied, and the
  claim is the **ratio to the peak** and not the residual —
  `test_the_kernel_self_force_correction_leaves_the_body_with_no_net_internal_force`.
- **I9.** A non-positive strand radius is refused, on the host path and on the device path —
  `test_a_vertex_on_the_centroid_is_refused_before_any_kernel_launches`,
  `test_the_position_precision_changes_whether_the_strand_refusal_fires`.
- **I10.** The frozen reference is driven in **float64 positions**, never float32-rounded input —
  `test_the_nucleus_references_are_driven_in_float64_positions`. This is
  `PROPOSAL-position-precision-for-the-cuda-port.md` §4's hazard, and G2 and G3 both shipped the
  case-level version; it is repeated here because the hazard is per case builder.
- **I11.** The state's branch counts are fields of the state and are asserted —
  `test_the_default_nucleus_state_visits_both_strain_signs_and_both_strand_branches`.

## 6. Source evidence class and known retractions

Nothing is inherited from the provider tree by this lane, so there is no inherited evidence class to
carry. What **is** carried, and where it was checked in this repository:

- **`ALEPH-PORT-3501` §4.5's self-force table**, including its own two corrections: an early draft
  quoted `3.3e-15` against `3.423` (the chromatin term in isolation) as though it were the assembled
  owner's figure, and its replacement `4.26e-14` **reproduced under no code path that entry ever
  shipped**. The reproducible assembled figure is `4.0856e-13` pN against a peak of `376.31` pN.
  Read here, not quoted from memory.
- **`ALEPH-PORT-3501` §9.1 M3 and M10.** M3 (subtract the centroid first) is an **equivalent
  mutant**, proved so in §9.2 by the closed-surface identity `Σ_f n_f dA = 0`. That is the same
  identity that makes `PositionPrecision.LOCAL_F32` legitimate for this law, so this port's
  local-coordinate mode has a proof behind it and not a hope (§10.2). M10 (area × mean radius) is
  the impostor that survived the entire `R³` oracle, and it is ported here as a **shipped wrong
  kernel** rather than as a docstring claim.
- **`docs/design/GPU_STATE_2026-07-31.md` and commit `636b0c8`** — the `ORDERED` 0.00 / `ATOMIC`
  10.00 ULP measurement on the A5000, and the retraction that made "every scatter destination comes
  from `backend.zeros`" a rule. `636b0c8` retracts `278e6e5`, which reported a GPU-only backend
  defect that was a caller error.
- **No retraction was searched for in `/Users/sw1/ffn_cellsim`**, because nothing from that tree is
  used. `STATE.md` §(c)'s not-quotable list is not engaged: no magnitude from that project appears
  in this entry or in the code it authorises.
- **Evidence class of everything this entry produces: `STRUCTURAL`.** Numerical agreement between
  two implementations of Aleph's own law. Nothing here is a measurement of a nucleus, and the cards
  remain `UNSOURCED` exactly as `ALEPH-PORT-3501` §14 says.

## 7. Independent oracle or derivation

Five, and they are deliberately of different kinds, because `ALEPH-PORT-3501` §9.1 M10 is this
project's clearest demonstration that an oracle can be real, passing, and blind.

1. **The frozen NumPy law**, called and never re-implemented: `volume_energy_and_forces`,
   `chromatin_energy_and_forces`, and — for the correction — `NucleusOwner.chromatin_term` itself,
   on an owner constructed from the case's own state. It is the definition of the physics and it
   already carries 47 controls of its own.
2. **`det(M)`, with the impostor evaluated inside the control and required to fail.** This is the
   oracle with discriminating power and §9.1 M10 is the evidence that the cube law is not.
3. **`s³`, with the impostor required to PASS.** Shipped as an oracle *and* as a demonstration of
   its own narrowness — the same discipline the frozen control adopted after M10 survived it.
4. **Exact identities that need no tolerance:** the reflection negation (`==`), the inert model
   (`==`), and the refusal domain.
5. **The two anchors `ALEPH-PORT-3501` measured and this port must reproduce through the kernel:**
   the chromatin radial stiffness `k_chromatin = 13316.236483` pN/µm on a level-2 sphere at
   `R = 3 µm`, `G = 120 pN/µm²`, `r₀ = 3 µm`, and the sphere volume `4πR³/3` recovered **from below**
   at level 4.

**What is deliberately NOT an oracle**, and this is the sharpest instruction this port carries:

> **`sum_to_zero`'s net force is NOT a parity target.** `ALEPH-PORT-3501` §4.5 measured that it is a
> fifteen-order cancellation — peak `376.31` pN, net `~4e-13` pN — in which **a single ULP in
> `dE_dV` moves the net between `2.9e-13` and `4.0e-13`**. A kernel's different summation order will
> legitimately give a different residual, and float32 arithmetic will give one many orders larger.
> The stable claim is **the ratio to the peak**, and the frozen control's own tolerance is
> `3.76e-07` pN (`1e-9 × peak`). §10.4 states what the kernel figure may and may not be compared
> against.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/runtime/test_nucleus_law_parity.py::test_the_true_volume_kernel_agrees_with_the_frozen_numpy_law` | Every field of the volume case within its declared ULP budget, on the deformed default state, with the measured figure reported. |
| Positive | `::test_the_true_chromatin_kernel_agrees_with_the_frozen_numpy_law` | The same for the fixed-anchor chromatin case. |
| Positive | `::test_the_true_corrected_chromatin_kernel_agrees_with_the_owners_own_term` | The same for the **owner's** `chromatin_term`, i.e. with the live-centroid correction, against `NucleusOwner.chromatin_term` itself. |
| Positive | `::test_a_nucleus_law_compared_with_itself_is_exactly_zero_ulp` | The harness's zero point on all three cases. |
| Positive | `::test_the_device_volume_law_has_exactly_one_kernel_and_the_nucleus_launches_it` | **I1** — the shared-kernel identity, asserted by recording the kernel objects a candidate evaluation launches, and the NumPy function-object identity restated beside it. |
| Positive | `::test_a_linear_map_scales_the_kernel_volume_by_its_determinant` | **I2** — three maps, two volume-preserving, with the area-ratio vacuity guard, and the impostor kernel required to fail. |
| Positive | `::test_the_kernel_volume_falls_as_the_cube_and_the_impostor_passes_that_too` | **I3** — and the impostor passing it is asserted, not described. |
| Positive | `::test_a_reflection_negates_the_kernel_volume_exactly` | **I4** — `==`, both halves. |
| Positive | `::test_the_kernel_volume_recovers_four_thirds_pi_r_cubed_from_below` | The level-4 sphere, two-sided so the sign is pinned. |
| Positive | `::test_the_kernel_chromatin_reproduces_the_recorded_radial_stiffness` | §7 anchor: `k_chromatin = 13316.236483` pN/µm through the kernel energy. |
| Positive | `::test_the_kernel_chromatin_response_is_mesh_resolution_independent` | **I7** — refining by one level moves the energy of a fixed uniform radial strain by under 1%. |
| Positive | `::test_the_kernel_volume_force_is_minus_the_gradient_of_the_kernel_energy` | **I5** — order 2 for both terms, step from the curve below. |
| Positive | `::test_the_central_difference_step_is_chosen_from_a_measured_error_curve` | The curve itself, with its order-2 window and its floor, printed rather than assumed. |
| Positive | `::test_the_kernel_self_force_correction_leaves_the_body_with_no_net_internal_force` | **I8** — the corrected net as a **ratio to the peak**, with the uncorrected net reported beside it. |
| Positive | `::test_the_default_nucleus_state_visits_both_strain_signs_and_both_strand_branches` | **I11** — branch coverage as numbers. |
| Positive | `::test_the_nucleus_references_are_driven_in_float64_positions` | **I10** — the proposal §4 hazard, at case level. |
| Positive | `::test_every_nucleus_scatter_destination_is_allocated_through_the_backend` | The `636b0c8` rule, by type. |
| Positive | `::test_the_nucleus_report_names_the_position_precision_and_the_scatter_mode` | A stored ULP figure never travels without its accumulation path or its position mode. |
| Positive | `::test_the_position_precision_is_a_parameter_and_not_a_constant_for_these_laws` | All three modes run and give **different** figures, so the parameter is live. |
| Positive | `::test_the_declared_nucleus_budgets_sit_inside_what_roundoff_explains` | No budget exceeds its own round-off explanation. |

## 9. Deliberately failing negative control

Nine wrong kernels ship in `aleph/runtime/law_kernels.py` beside the true ones, compiled by the same
`load()` and launched down the same drivers. Each is a *plausible hand-translation error*, and they
are chosen so that the three graded channels each catch a different one — which is G2's and G3's
finding carried forward.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `::test_a_volume_kernel_that_used_area_times_mean_radius_is_caught` | `VolumeVariant.WRONG_AREA_TIMES_RADIUS` — **`ALEPH-PORT-3501` §9.1 M10 as a kernel.** Exact for a sphere, translation invariant, degree-three homogeneous, and **not the enclosed volume of anything else**. Caught by the parity gate in both channels *and*, more importantly, by `det(M)` while `s³` passes it. |
| Negative (must fail) | `::test_a_volume_kernel_that_cycled_its_gradient_wrongly_is_caught` | `WRONG_GRADIENT_CYCLE` — corners 1 and 2 receive each other's cross product. **Energy bit-identical**; caught on the force alone. |
| Negative (must fail) | `::test_a_volume_kernel_that_dropped_the_half_in_its_energy_is_caught` | `WRONG_ENERGY_HALF` — M1 as a kernel. **Forces bit-identical**; caught on the energy alone. |
| Negative (must fail) | `::test_a_volume_kernel_that_pushes_where_the_law_pulls_is_caught` | `WRONG_FORCE_SIGN` — M5 as a kernel, `+dE/dV·∂V/∂x`. **Energy bit-identical**; caught on the force alone, and by `F = −∇E`. |
| Negative (must fail) | `::test_a_chromatin_kernel_that_dropped_the_half_in_its_energy_is_caught` | `ChromatinVariant.WRONG_ENERGY_HALF`. **Forces bit-identical**; energy alone. |
| Negative (must fail) | `::test_a_chromatin_kernel_that_weighted_every_strand_equally_is_caught` | `WRONG_UNWEIGHTED` — M4 as a kernel, `w_v := 1`. Caught by the parity gate **and** by the mesh-resolution-independence oracle, which is the control M4 was written for. |
| Negative (must fail) | `::test_a_chromatin_kernel_that_forgot_to_normalise_its_direction_is_caught` | `WRONG_UNNORMALISED_DIRECTION` — the force keeps `d_v` instead of `d_v/r_v`. **Energy bit-identical**; force alone. |
| Negative (must fail) | `::test_a_kernel_that_dropped_the_self_force_correction_is_caught` | `SelfForceVariant.WRONG_NO_CORRECTION` — M9 as a kernel. Caught against the owner's own term, and by the net-force control at **2.978 pN against a 3.42 pN peak**. |
| Negative (must fail) | `::test_a_kernel_that_subtracted_a_scalar_mean_is_caught` | `WRONG_SCALAR_MEAN` — `f.mean()` where the law says `f.mean(axis=0)`. The one-character version, and the more plausible of the two. |
| Negative (must fail) | `::test_a_zero_ulp_budget_fails_the_true_nucleus_kernels` | The vacuity guard on the verdict: with `ulp_budget = 0.0` the **correct** kernels fail, because float32 storage is not bitwise float64. |
| Negative (must fail) | `::test_no_nucleus_mutant_escapes_through_the_atomic_floor` | The 10 ULP atomic floor is a floor, not an amnesty. |
| Negative (must fail) | `::test_a_vertex_on_the_centroid_is_refused_before_any_kernel_launches` | **I9 host half** — the same `ValueError` the reference raises, in float64, before anything launches. |
| Negative (must fail) | `::test_the_position_precision_changes_whether_the_strand_refusal_fires` | **I9 device half** — the kernel's own violation flag, on a configuration float64 calls legal. |

### 9a. Declared blind spots — which states cannot see which mutant

G1's finding, sharpened by G3: *a parity case is only as strong as the branches its state visits*,
and G3 found a mutant invisible on **the state production actually evaluates**. This law has the
same shape of problem and it is worse, because the blind state is the one the constructor returns.

| Blind spot (must hold) | Location | Asserts |
|---|---|---|
| **The rest nucleus is blind to every constitutive mutant, and is not a valid parity state at all.** | `::test_the_rest_nucleus_is_blind_to_every_constitutive_mutant_and_is_not_a_parity_state` | `build_spherical_nucleus` takes `V_ref` and `r₀` from *its own* configuration, so both interior energies are zero at construction — `E_vol` **exactly** `0.0` and `E_chr` ~`1e-28` against a deformed `~66` pN·µm. A factor-of-two error in either energy moves zero to zero. And the gate cannot be run there: `ulp_at_scale` divides by a reference scale of exactly zero and returns `inf` for the **correct** kernel. Asserted, so the state is refused rather than met later on a kernel that mattered. |
| Blind spot (must hold) | `tests/runtime/test_nucleus_interior_law_parity.py::test_the_kernel_volume_falls_as_the_cube_and_the_impostor_passes_that_too` | The `s³` oracle passes `WRONG_AREA_TIMES_RADIUS` at every scale factor. This is M10's finding reproduced on the device, and it is why I2 exists. |
| Blind spot (must hold) | `::test_the_unweighted_chromatin_is_invisible_on_a_uniform_weight_state` | On a state whose reference vertex weights are all equal, `WRONG_UNWEIGHTED` differs from the true kernel by a constant factor that the *shape* of the force cannot see; only the magnitude does, and on a state where the weights are normalised to 1 it is bit-identical. |

## 10. Numerical and precision envelope

**Working precision.** Positions arrive in the declared `PositionPrecision`. Per-face triple
products are widened to float64 before reduction, exactly as every energy in G1–G3 is; the corner
gradients, the chromatin per-site forces and every assembled force array are float32 —
`ALEPH-DQ-107`'s compute channel. The constitutive scalars `ε_V`, `E_vol` are float64 inside
`k_volume_constitutive_*`; `dE/dV` is rounded to float32 once, at the stage boundary, because it
multiplies a float32 gradient. The reference is float64 end to end.

### 10.1 Where the disagreement comes from, derived rather than assumed

Two mechanisms, and they are not the same as G2's:

**The volume geometry cancels against the world origin, cubically.** Each face's triple product is
`O(|p|³)` while their sum is `6V = O(R³)` for a body of radius `R` at distance `d` from the origin.
So the amplification is

```
amp_volume ≈ Σ_f |p₀·(p₁×p₂)| / |Σ_f p₀·(p₁×p₂)|
```

which is `O(1)` for a centred body and grows like `F·d³/(6V)` for a displaced one. This is
`ALEPH-PORT-3501` §10.4 — measured there in **float64**, where the relative error reaches `5.2e-11`
at a `1e3 µm` offset and `4.5e-01` at `1e6 µm`. Float32 has `2^29 ≈ 5.4e8` times the unit round-off,
so the same curve arrives about nine orders of magnitude earlier: **a cell at 50 µm is already deep
into it.** §13.4 measures it.

**The chromatin offset cancels linearly**, in the ordinary way: `|d(x_v − c)| ≲ u(|x_v| + |c|)`, and
the strain divides by `r₀`. This is the same shape as G1's tether and is mild by comparison.

Both amplifications are computed **from the host state by the case**, so they are properties of the
data and not numbers chosen after seeing the answer.

### 10.2 `LOCAL_F32` for this law is an equivalent mutant with a proof

`ALEPH-PORT-3501` §9.2 proves that subtracting a common point before forming the triple products
changes no value in exact arithmetic, from the closed-surface identity `Σ_f n_f dA = 0`:

```
(1/6) Σ_f ((a−c₀) × (b−c₀))·(c−c₀)  =  (1/6) Σ_f (a×b)·c      for every c₀
```

That entry calls the shipped global form **the worse-conditioned of the two** and declines to change
it, because changing shipped numerics for a regime nothing uses is not a porting lane's business.
Here it is not a change to shipped numerics: it is a **declared parameter of a new device path**, and
it is exactly the better-conditioned form that entry identified and could not adopt. So
`PositionPrecision.LOCAL_F32` is not a heuristic for this law, it is M3.

### 10.3 The two global reductions, and the synchronisation `ALEPH-PORT-3501` §11 predicted

That entry says, of a future device backend: *"the volume reduction is the term that needs attention:
it is a global scalar reduction over faces inside a force evaluation, which is the shape that costs a
synchronisation."* It is right, and there are **two** of them here:

| reduction | how | why not a host round trip |
|---|---|---|
| `Σ_f p₀·(p₁×p₂)` | `Backend.scatter_add` into a length-1 float64 slot | the constitutive kernel consumes it **on the device**; reading it to the host would put the strain and the energy on the host and make the constitutive stage untestable by a kernel mutant |
| `Σ_v f_v` per column | `Backend.sum(f, axis=0)` → a float64 device array of 3 | the mean-subtraction kernel consumes it on the device |

Both are ordered by construction (the CSR scatter and `sum_axis0`), so neither is exposed to the
atomic floor, and neither is a host synchronisation in this implementation. The cost on a real device
is a kernel-launch dependency chain of four stages per volume evaluation and three per chromatin
evaluation, and **that is not measured here** — §14.

### 10.4 What the self-force residual may and may not be compared against

`ALEPH-PORT-3501` §4.5, restated because it is the instruction this port most easily gets wrong:

- The **uncorrected** chromatin net is `2.9777` pN against a chromatin peak of `3.4230` pN. That
  ratio, 87%, is stable, is a property of the physics, and **is** a target.
- The **corrected** net is a cancellation residual. In float64 it is `3.33e-15` pN for the chromatin
  term and `4.0856e-13` pN for the assembled owner, and a single ULP of `dE_dV` moves the latter
  between `2.9e-13` and `4.0e-13`. In float32 it will be many orders larger and that is **not a
  defect**. The control asserts the **ratio**, with a budget derived from `N · eps32` and stated.
- The frozen control's `1e-9 × peak` tolerance is a **float64** property of a **float64** code path.
  Quoting it as a target for a float32 kernel would be an over-claim of exactly the kind
  `ALEPH-PORT-3501` corrected twice.

### 10.5 The rest configuration is not a parity state, and it is the shipped one

`ulp_at_scale` returns `inf` when the reference field is identically zero and the candidate is not.
`build_spherical_nucleus` returns a nucleus whose volume energy is **exactly** `0.0` and whose
volume force is exactly zero, because `V_ref` is taken from that very configuration. So the correct
kernel — whose float32 volume differs from `V_ref` in the seventh digit — reports `inf` there.

This is the same class of finding as G2's `test_a_pristine_sphere_is_not_a_valid_parity_state_for_the_bending_force`
and G3's `WRONG_MIDPOINT_ONLY`, and it is the sharpest instance of the three: **the blind state is
what the constructor hands you.** It is asserted rather than left to be discovered.

### 10.6 The atomic floor

`ATOMIC_SCATTER_ULP_FLOOR = 10.0`, measured on the A5000 (`636b0c8`). The volume force, the chromatin
force and the corrected force are scatter-assembled and take the floor in atomic mode; the energies
are reduced by ordered paths and keep their declared budgets. The volume's own `Σ_f` reduction is a
scatter into one slot and is therefore **also** exposed to the floor in atomic mode — recorded here
because it is the one place in this port where an *energy* depends on a scatter, which is a departure
from G1–G3 and is stated rather than silently inherited.

## 11. Production-backend residency and transfer

| Array | Where | Precision | Transfer |
|---|---|---|---|
| vertex positions | device | float32 or float64 per `PositionPrecision` | uploaded once per case by `backend.array` |
| connectivity | device, `wp.array` int32 | int32 | not a scatter destination; `WarpBackend` accepts float dtypes only, and connectivity carries no precision question |
| reference vertex weights `w_v` | device | float32 | uploaded once; **never recomputed in a kernel** |
| per-face triple product | device | float64 | reduced on the device |
| total triple, `E_vol`, `dE/dV` | device, **allocated by `backend.zeros`** | float64 / float64 / float32 | never leave the device between stages |
| per-corner volume gradient | device | float32 | never leaves the device |
| per-site chromatin energy | device | float64 | reduced by `Backend.sum` |
| chromatin `violation` flags | device, `wp.array` int32 | int32 | read to host once per launch, by the driver, to decide the refusal |
| assembled vertex forces | device, **allocated by `backend.zeros`** | float32 | read to host once, by the harness, for the comparison |

**The one rule this entry keeps from `636b0c8`:** every scatter destination comes from
`backend.zeros(...)`. A raw host ndarray works on warp's CPU device and fails on CUDA; `278e6e5`
reported that as a GPU-only backend defect and `636b0c8` retracted it — the API was right and the
caller was wrong. `test_every_nucleus_scatter_destination_is_allocated_through_the_backend` asserts
it by type.

The harness itself is host-side and is not production physics. It reads each graded field to the
host exactly once per case.

## 12. Comments and docstrings to discard

Nothing was read from the provider tree by this lane, so there is no source prose to discard. The
list is kept as a positive statement of what must **not** appear in the code this entry authorises,
and `tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package` is the
mechanical half:

- no provider repository name, module path, kernel name, gate name or branch name;
- no absolute path into another project;
- no provider datum as a literal;
- no claim that a number was measured on hardware this lane did not run on. Every ULP figure in the
  package's docstrings is measured here on warp's CPU device and says so, or is cited to
  `docs/design/GPU_STATE_2026-07-31.md` and `636b0c8` and says that.

`ALEPH-PORT-3501` §12 additionally discarded the provider's device vocabulary from the NumPy port —
the `wp.kernel` / `wp.tid()` / `wp.atomic_add` framing and the "one thread" note — on the grounds
that Aleph's version was a host reduction. **This entry is the point at which that vocabulary becomes
Aleph's own, arrived at independently**: the kernels here are written from §4's derivation, and the
one-thread-per-face decomposition is forced by the law rather than taken from anybody.

What replaces the discarded prose: §4's derivation restated in each kernel's own comment, §10's
conditioning restated on the budget constants, and the three findings this port inherits from G1–G3
restated where a reader will meet them.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Measured 2026-08-01 on **warp's CPU device** (this machine's warp build reports *"CUDA not enabled in this build"*), warp 1.15.0, float32 compute, `ScatterMode.ORDERED`, `PositionPrecision.GLOBAL_F32` unless stated. `tests/runtime/test_nucleus_law_parity.py`: **37 passed**. `tests/runtime` as a whole: **375 passed, 3 skipped** (was 338 + 3 before this entry). `tests/vertical`: **1,358 passed**, unchanged — the frozen references were read and never edited. Mutation study: **10 of 10 killed** (§13a). **Every exact anchor `ALEPH-PORT-3501` recorded reproduced through the kernels**, at float32 precision, and the self-force correction survived. The status stays `PROPOSED` because no CUDA device was touched and no human has reviewed it. |
| Reviewer | **Agent-proposed and unratified.** Written and measured by a subagent of session `46143f30`; no PI review, no `decided_by` field, and none may be added by an agent. |
| Rollback | Delete `tests/runtime/test_nucleus_law_parity.py` and the `ALEPH-PORT-3605` blocks appended to `aleph/runtime/law_kernels.py` and `aleph/runtime/law_cases.py`. Nothing else imports them. What is lost: the nucleus returns to having no device physics, and — because this entry introduced the **only** device implementation of the enclosed volume — so does any future membrane-side volume case. |

### 13.1 The parity measurement

Default state: the level-2, `R = 3 µm` nucleus at `(1.12, 0.93, 1.04)` stretch and `0.03 µm` jitter,
seed `20260731` — **`tests/vertical/test_nucleus_owner.py::build_test_nucleus`'s own configuration**,
reproduced bit-for-bit and asserted to be so. 162 vertices / 320 faces, `ε_V = +8.34%`, 76 compressed
and 86 stretched strands.

| kernel | `energy_pn_um` | `forces_pn` | verdict |
|---|---|---|---|
| volume `TRUE` | **4.76** | **2.90** | pass |
| volume `WRONG_AREA_TIMES_RADIUS` | 7,086,813 | 3,005,129 | **fail, both** |
| volume `WRONG_GRADIENT_CYCLE` | 4.76 *bit-identical* | **11,238,337** | **fail, force only** |
| volume `WRONG_ENERGY_HALF` | **8,388,618** | 2.90 *bit-identical* | **fail, energy only** |
| volume `WRONG_FORCE_SIGN` | 4.76 *bit-identical* | **16,777,218** | **fail, force only** |
| volume declared budget | 16 | 12 | — |
| volume round-off explains | 52.99 | 27.80 | — |
| chromatin `TRUE` | **1.29** | **7.10** | pass |
| chromatin `WRONG_ENERGY_HALF` | **8,388,611** | 7.10 *bit-identical* | **fail, energy only** |
| chromatin `WRONG_UNWEIGHTED` | 3,831,837 | 4,002,208 | **fail, both** |
| chromatin `WRONG_UNNORMALISED_DIRECTION` | 1.29 *bit-identical* | **19,863,250** | **fail, force only** |
| owner term `TRUE` (`SelfForceVariant.TRUE`) | **1.29** | **6.45** | pass |
| owner term `WRONG_NO_CORRECTION` | 1.29 *bit-identical* | **45,049** | **fail, force only** |
| owner term `WRONG_SCALAR_MEAN` | 1.29 *bit-identical* | **28,460** | **fail, force only** |
| chromatin declared budget | 8 | 11 | — |
| chromatin round-off explains | 46.27 | 11.97 | — |

**Six of the nine wrong kernels are bit-identical to the true one in the energy channel, and two are
bit-identical in the force channel.** An energy-only gate certifies six of them; a force-only gate —
which is what most parity gates are, because forces are what an integrator consumes — certifies two.
That is G2's and G3's finding reproduced in a third and a fourth law, and it is the whole argument
for `compare_law_case` grading fields separately.

**One budget is tight and it is tight by necessity.** The chromatin force measures 7.10 against a
round-off explanation of 11.97, i.e. 59% of it, so the declared 11 carries a margin of 1.55× where
`ALEPH-PORT-3601` could afford 2.5–3×. Widening it further would put the budget above its own
justification, which `LawFieldParity.budget_within_roundoff_explanation` exists to prevent. Recorded
rather than smoothed over.

### 13.2 The exact anchors, and what "exact" means through a float32 kernel

| oracle | frozen float64 | through the kernel | verdict |
|---|---|---|---|
| `det(M)`, volume-preserving stretch `diag(2, 1, 0.5)` | `1.1e-16` | **`2.5e-08`** | identity survives |
| `det(M)`, volume-preserving shear | `1.1e-16` | **`8.1e-09`** | identity survives |
| `det(M)`, general contraction `diag(1.5, 0.56, 1)`, `det = 0.84` | `1.1e-16` | **`2.9e-08`** | identity survives |
| …the impostor kernel on the same three maps | — | off by **0.663 / 0.106 / 0.334** | **caught** |
| `s³` uniform scaling, six factors | `1.6e-15` | **`3.4e-08`** worst | identity survives |
| …the impostor kernel on the same six factors | — | **`5.4e-08`** worst — **passes** | **blind, as designed** |
| reflection: `V(reflected) == -V` | exact | **exact, bitwise** | `==` |
| reflection: the areal kernel's energy | exact | **exact, bitwise** | `==` |
| sphere limit, level 4 | under by `0.216%` | **under by `0.216%`** | from below |
| `k_chromatin` radial stiffness | `13316.236483` pN/µm | **`13315.864710`** (`2.8e-05`) | reproduced |
| `K_V = 0`, `G = 0` | exactly `0.0` | **exactly `0.0`, budget `0.0`** | exact |
| `F = −∇E`, volume, `POSITIONS_F64` | order 2 | **1.9984 / 1.9938 / 1.9754** | order 2 |

**Two rows are the point of the table.** The `det(M)` identity is exact in float64 and holds to
`~3e-08` through a float32 kernel — the identity *survives the port at float32 precision*, and its
**discriminating power is undiminished**: the gap between the true kernel and the impostor is seven
orders wide. And the `s³` row is the one this entry is proudest of asserting: the impostor **passes**
it, on the device, exactly as `ALEPH-PORT-3501` §9.1's M10 passed it on the host, so the narrowness
of the cube law is a measurement here and not a sentence inherited from a docstring.

### 13.3 The gradient curves, and the one that has no order-2 window

Measured through the kernels at `POSITIONS_F64` on a strongly deformed state (`scale=(1.45, 0.78,
1.10)`, `jitter=0.06 µm`), because a finite difference against a float32 *energy* needs a signal well
above that energy's own round-off floor.

**Volume, against the kernel's own energy:**

```
h          error         order
2.0e-01    1.065250e-03    --
1.0e-01    2.666078e-04   1.9984
5.0e-02    6.694047e-05   1.9938
2.5e-02    1.702316e-05   1.9754
1.25e-02   4.543770e-06   1.9055
6.25e-03   1.424508e-06   1.6734   <- leaving the clean band
3.0e-03    6.233420e-07   1.1261
1.0e-03    4.099139e-07   0.3815   <- the floor
```

`VOLUME_GRADIENT_STEP_UM = 5e-2` is taken from the middle of that band, not from a formula.

**`POSITIONS_F64` and not the shipped default, and that is a finding rather than a convenience.**
Under `GLOBAL_F32` the volume energy carries ~5 ULP whose absolute size exceeds the finite-difference
signal at every step, and **no order-2 window exists at all**. That is
`docs/design/BUILD_PLAN-2026-07-31-engine-gaps.md`'s second precision channel — "it is the energy
channel, not the force" — arriving in a third law, and it is the reason this entry puts no relaxation
on device and drives no descent from a kernel energy.

**Chromatin, differencing the LIVE-centroid energy:**

```
h          corrected      uncorrected
2.0e-01    1.742382e-04   1.194961e-02
1.0e-01    6.721297e-05   1.205664e-02
5.0e-02    3.307760e-04   1.245463e-02
2.5e-02    1.364568e-04   1.226031e-02
1.25e-02   1.002859e-04   1.222414e-02
6.25e-03   3.500923e-04   1.247394e-02
```

The uncorrected column is **flat to within 5% over a 32× reduction in `h`** — a plateau, which
`ALEPH-PORT-3501` §10.3 identifies as the signature of a gradient wrong by a fixed amount, as against
round-off, which *grows*. The corrected column is two orders below it. **That shape, not a single
point, is what the control asserts**, which is the discipline
`test_the_error_curve_falls_rather_than_plateauing` established on the host.

**A trap this lane walked into and is recording so nobody repeats it.** Differencing the
*live-centroid* energy against the *fixed-anchor* force produces exactly the same plateau, at exactly
the same magnitude, and looks precisely like a wrong correction. It is a wrong *probe*: the two laws
differentiate different energies. `_restate(..., live_anchor=...)` in the test file exists to make
the choice explicit at every call site rather than implicit in a default.

### 13.4 The three position modes — and a counterexample to option D's shape

Force ULP against the frozen float64 law, at the origin and with the same nucleus at `(50, 0, 0) µm`:

| | volume E | volume F | chromatin E | chromatin F |
|---|---|---|---|---|
| origin, `GLOBAL_F32` | 4.76 | 2.90 | 1.29 | 7.10 |
| origin, `LOCAL_F32` | 2.25 | 2.57 | 1.43 | 5.47 |
| origin, `POSITIONS_F64` | **0.00** | **1.01** | 1.43 | 5.47 |
| 50 µm, `GLOBAL_F32` | **40.84** | **53.02** | 4.66 | **65.61** |
| 50 µm, `LOCAL_F32` | 2.25 | **2.57** | 1.43 | 5.47 |
| 50 µm, `POSITIONS_F64` | **0.00** | 11.37 | 1.43 | 5.47 |

**Two rows are findings and both are raised, not settled.**

**(1) `LOCAL_F32` recovers the at-origin figure exactly, in every column — and for this law that has
a proof rather than a measurement behind it.** `ALEPH-PORT-3501` §9.2 shows that subtracting a common
point before forming the triple products changes no value in exact arithmetic, from the closed-surface
identity `Σ_f n_f dA = 0`; it is that entry's **equivalent mutant M3**. That entry called the shipped
global form the worse-conditioned of the two and declined to change it, because changing shipped
numerics for a regime nothing uses is not a lane's business. Here it is not a change to shipped
numerics — it is a declared parameter of a new device path.

**(2) `POSITIONS_F64` is WORSE than `LOCAL_F32` on the volume force: 11.37 against 2.57.** This is a
counterexample to the *shape* of option D in
`docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` §5, which is defined as *"the only
float64 arithmetic on the device is the subtraction that forms edge vectors and site separations"*
and is recommended in §6 on the grounds that *"the cancellation is removed at the point where it
happens, which is the only place it happens."*

> **The volume law has no such subtraction.** Its cancellation is a **summation**: the per-face triple
> products are `O(d³)` and their sum is `O(R³)`, and the per-corner gradient rows are `O(d²)` while the
> assembled vertex gradient is `O(A)`. Both cancellations happen in a **reduction**, not in a
> difference, so widening the differencing does not reach them. Option B does, exactly and for free.

That is one measured law against the reasoning behind one recommendation, not a refutation of option
D — for the membrane's stencil and both connectors, D behaves as §5 says it does. It is added to the
record here because `PROPOSAL-...` §7 already names lane ledgers as where measured mode costs are
reported, and because a PI answering §6 should have it. **This lane takes no position on the answer.**

### 13.5 The self-force correction

| variant | net (component max) | peak | net / peak |
|---|---|---|---|
| frozen NumPy, fixed anchor (uncorrected) | `2.977746` pN | `3.4233` pN | `0.8698` |
| frozen NumPy, `NucleusOwner.chromatin_term` | `3.330669e-15` pN | `3.4233` pN | `9.7e-16` |
| **kernel `SelfForceVariant.TRUE`** | **`2.462868e-06` pN** | `3.4230` pN | **`7.19e-07`** |
| kernel `WRONG_NO_CORRECTION` | `2.977710` pN | `3.4233` pN | `0.8698` |
| kernel `WRONG_SCALAR_MEAN` | `1.881086` pN | `3.4301` pN | `0.5484` |

**The correction survived the port.** Three things are worth reading off it:

1. The uncorrected kernel reproduces the frozen float64 defect to five digits — `2.977710` against
   `2.977746` — so the kernel has ported the *defect* faithfully as well as the fix, which is what
   makes the comparison a measurement rather than a coincidence.
2. The corrected ratio is `7.19e-07` against a budget **derived** as `N · eps32 = 162 × 1.19e-7 =
   1.93e-05`, not borrowed. The frozen control's `1e-9 × peak` is a float64 property of a float64 path
   and is **eight orders below what float32 can reach**; quoting it as a target here would be the
   over-claim `ALEPH-PORT-3501` §4.5 corrected twice.
3. `WRONG_SCALAR_MEAN` is the instructive one. It *looks* like the correction — it removes a net
   quantity and the residual falls from `2.98` to `1.88` — and it removes the wrong projection, so the
   body still drifts and each component is now contaminated by the other two.

The assembled figure, for completeness and labelled as the **hybrid** it is (kernel interior + NumPy
lamina, because `ALEPH-PORT-3605` does not port the lamina): net `4.055910e-06` pN against a peak of
`376.3079` pN, ratio `1.08e-08`. The all-float64 reference on the same configuration is
`4.085621e-13` pN at ratio `1.06e-15`, which **reproduces `ALEPH-PORT-3501` §4.5's `4.0856e-13`
exactly** — confirming that this entry's state really is that entry's configuration and that the
comparison above is apples to apples.

### 13a. Mutation study — 10 planted defects, 10 killed

Run with `PYTHONDONTWRITEBYTECODE=1`. Each mutant is applied to the **shipped** source at
`7775eae`, the suite is run, the mutant is reverted with `git checkout`, and both target files are
confirmed byte-identical with `git diff --exit-code` **before and after every mutant**, so a
concurrent lane's edit can neither be mistaken for a result nor destroyed by the study. **No test
file was modified at any point** — mutating a control measures the mutation and not the suite.

| # | Planted defect | Killed by |
|---|---|---|
| M1 | **the nucleus launches its OWN duplicate volume kernel** — a byte-identical copy added to `_define()` and launched instead of the shared object | **1 test** |
| M2 | `POSITIONS_F64` silently uploads float32 and launches the float32 volume kernel | 3 tests |
| M3 | the volume driver ignores the requested variant and always runs the true stages | 6 tests |
| M4 | the chromatin driver ignores the requested variant | 5 tests |
| M5 | the self-force driver always applies the TRUE correction whatever was asked for | 5 tests |
| M6 | **the kernel's own violation flag is never read, so the device-side refusal is dropped** | **1 test** |
| M7 | **the frozen volume reference is handed float32-rounded positions** | 3 tests |
| M8 | a volume scatter destination allocated as a host ndarray instead of `backend.zeros` | 2 tests |
| M9 | **the chromatin weights come from the DEFORMED positions, not the frozen rest ones** | **1 test** |
| M10 | the self-force correction divides by 1 instead of by the vertex count | 4 tests |

**M1 is this entry's own contribution to the list, and it is the reason §3.1's control is written by
watching rather than by reading.** The duplicate computes *exactly the same numbers* — it is a
copy — so every parity figure, every oracle and every mutant control stays green. The only thing that
changes is **which object ran**, and
`test_the_device_volume_law_has_exactly_one_kernel_and_the_nucleus_launches_it` is the one control
that can see it. That is the device form of the property `ALEPH-PORT-3501` §3.1 lost once on the host
and got back with a function-object identity assertion.

**M9 is the finding of this study, and it is about the method rather than the code.** It corrupts the
*state* — the chromatin weights — and both the reference and the kernel read the same state, so
**the parity gate is structurally incapable of seeing it**: every ULP figure stays exactly where it
was. It is killed by the recorded anchor `2.9777 pN` in the self-force control, and by nothing else.

> **A law-parity gate cannot see a defect in the state it grades on. Only an anchor to an
> independently recorded number can.** That is the argument for §7 oracle 5 and for quoting
> `ALEPH-PORT-3501`'s measured figures rather than re-deriving them here, and it generalises to every
> case builder G5 will write.

**Three mutants are killed by exactly one control each (M1, M6, M9), and that is a finding rather
than a comfort** — `ALEPH-PORT-3604` §13a made the same observation about its own M4. A control that
is the sole killer of a real defect is one deletion away from that defect surviving, and all three
are named here so a future refactor meets the fact before it deletes the test.

**M6 is inherited from `ALEPH-PORT-3604` deliberately.** It is the failure mode the twice-refusing
design exists to prevent: the host refuses in float64, the kernel disagrees, and if nobody reads the
flag the run continues on a configuration the law calls impossible — with exactly zero force at those
sites, which looks like a strand at rest rather than like a fault.

**M7 is `ALEPH-PORT-3603`'s and `-3604`'s M7 repeated, again deliberately**, because it is the hazard
`PROPOSAL-position-precision-for-the-cuda-port.md` §4 names and it is invisible in the artefact: a
reference quietly given the same float32 input as the kernel loses the same digits and the gate
prints a small healthy number.

## 14. Honest limits

What this entry does **not** establish:

1. **No CUDA device was touched.** A grant exists for `gpu-sungwook` device 0 until
   2026-08-01 09:00 and CUDA there is confirmed working inside a Slurm allocation, **but that host
   has no Python environment**, and an agent may never write the authorization record
   (`CLAUDE.md` §3, `GPU_POLICY.md`). Every number here is from **warp's CPU device**, where a kernel
   launch is a serial loop. It is evidence that the kernels' *algebra* matches the reference, and it
   is `UNVERIFIED` as evidence about CUDA. In particular the `ORDERED`/`ATOMIC` distinction is
   reproducible here only by accident of the schedule, as `measure_scatter_determinism` already says.
2. **The lamina is not ported.** `ALEPH-PORT-3605` is the *interior*. There is no fully-kernel
   `NucleusOwner` after this entry, and every assembled-owner figure in §13.5 is a **hybrid** —
   kernel interior, NumPy lamina — and is labelled as one. A future entry that ports the lamina must
   re-measure the assembled net rather than quote §13.5's.
3. **The two global reductions are not benchmarked.** `ALEPH-PORT-3501` §11 predicted that the volume
   reduction *"is the shape that costs a synchronisation"*, and this entry implements it as an
   ordered scatter into one slot and a `sum_axis0` — device-resident and deterministic, so it is not
   a host round trip *in this implementation*. What it costs on a real device is a kernel-launch
   dependency chain of four stages per volume evaluation and three per chromatin evaluation, and
   **that has not been measured on any device**.
4. **The declared budgets are engineering numbers, not physics**, and they belong to **one state**.
   The volume budget in particular depends on `1/|ε_V|`, which diverges as the nucleus approaches its
   rest volume, and on the distance from the world origin **cubically**. A case on another state must
   pass its own budget.
5. **`K_V`, `G`, `r₀` and `V_ref` remain `UNSOURCED`.** No stiffness computed from them — including
   the `13316.236483` pN/µm anchor reproduced in §13.2 — may be reported as a property of a real
   nucleus, in any document, at any rung. That prohibition is `ALEPH-PORT-3501` §14's and it travels
   with the numbers.
6. **The negative controls are as strong as the state they run on**, and two of them are provably
   blind on a state this package can build: `WRONG_UNWEIGHTED` is bit-identical on uniform weights,
   and the impostor volume is invisible to the cube law. Both are asserted rather than left to be
   met later. **And the state the constructor returns is blind to every constitutive mutant in both
   laws** — §10.5 — which is the sharpest instance of G1's limit this project has produced.
7. **A parity gate cannot see a defect in its own state.** §13a M9 is the measurement of that, and it
   is not fixed here: nothing structurally prevents a case builder from corrupting a material
   quantity that both sides then read. The recorded anchors are the only defence and they are not
   mechanically enforced.
8. **`ports/ledger/INDEX.md` was not regenerated.** It is already red from other lanes for
   `ALEPH-PORT-3601`–`-3604` and the scenario lane's uncommitted entries; `-3605` joins that
   already-red assertion as a fifth *name*, not as a new red test. Reported, not fixed, on
   instruction and per `CLAUDE.md` §1.
9. **The float32 `WRONG_GRADIENT_CYCLE` and `WRONG_AREA_TIMES_RADIUS` kernels have no float64
   twin.** They are graded at `GLOBAL_F32` only. The `det(M)` and `s³` oracles they exist for are
   exact in any precision, so a float64 copy would add a second wrong kernel without adding a second
   question — but it does mean the `POSITIONS_F64` path has two fewer negative controls than the
   `GLOBAL_F32` path, and that is stated rather than left to be counted.
10. **§13.4's counterexample is one law.** That `POSITIONS_F64` is worse than `LOCAL_F32` for the
    volume force is measured, and the mechanism — a cancellation in a reduction rather than in a
    difference — is derived. Whether it generalises to other reduction-shaped laws is **not**
    established here, and nothing in this entry answers
    `PROPOSAL-position-precision-for-the-cuda-port.md`.
