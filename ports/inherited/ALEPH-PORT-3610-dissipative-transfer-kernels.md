# ALEPH-PORT-3610 â the dissipative transfer family as Warp kernels

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3610` |
| Lane | `46143f30` Lane G5b (CUDA track), Track G task **G5b** |
| Status | `PROPOSED` |
| Written | `2026-08-01` â **before the code**, per `PLAN.md` Â§0.2.5 |
| Port class | `RE-DERIVED` |
| Depends on | `ALEPH-PORT-3601` (the law-level parity harness and its `LawCase` contract), `ALEPH-PORT-3604` (the per-site launch shape, the branch-count discipline, and the Newton-mirror measurement this entry extends), `ALEPH-PORT-3606` (the rows-vs-kernels census discipline, and the vacuity finding). All three are **reused, never duplicated**. |
| Exists because | Fourteen of the thirty wired rows are on kernels after `-3604` and `-3606`. **Sixteen remain**, and every law ported so far has been a *spring* â a potential in a separation. This entry takes the first family that is not: a **resistance in a velocity**, with no stored energy at all. |

---

## 1. Aleph API

```python
from aleph.runtime.law_kernels import (
    BlockDiagonalDragVariant,   # TRUE | WRONG_ISOTROPIC | WRONG_AXIS_SWAPPED
                                #      | WRONG_UNWEIGHTED | WRONG_SIGN
    DragDissipationVariant,     # TRUE | WRONG_RAYLEIGH_HALF
    DenseResistanceVariant,     # TRUE | WRONG_NODE_DIAGONAL | WRONG_COMPONENT_STRIDE
                                #      | WRONG_SIGN
    load, loaded_warp,          # unchanged, extended to compile these kernels
)

from aleph.runtime.law_cases import (
    DISSIPATIVE_TRANSFER_BATCH,          # the batch census: class -> rows -> kernel family
    DEFAULT_BLOCK_DIAGONAL_DRAG_ULP_BUDGET,
    DEFAULT_DENSE_RESISTANCE_ULP_BUDGET,
    ANISOTROPIC_PARALLEL_DRAG_PN_S_PER_UM2,
    ANISOTROPIC_RATIO,
    ISOTROPIC_DRAG_PN_S_PER_UM,
    MEDIUM_VISCOSITY_PN_S_PER_UM2,
    MEDIUM_RADIUS_UM,
    transfer_batch_connector,            # build one member through ITS OWN frozen builder
    anisotropic_drag_state,              # reproducible tensor-drag geometry + velocity
    isotropic_drag_state,                # the same builder at ratio = 1, L = 1
    anisotropic_drag_law_case,
    isotropic_drag_law_case,
    block_diagonal_drag_kernel_site_forces,
    dense_resistance_state,              # a real ExteriorStokesMedium and its operator
    dense_resistance_law_case,
    dense_resistance_kernel_site_forces,
    refuse_non_unit_tangents,            # the host refusal, on the LAUNCH path â see Â§5
)
```

Nothing outside this list is covered. This entry does **not** authorise any change to
`aleph/vertical/**` â `connectors_fluid.py`, `connectors_transfer.py`, `medium.py` and `wiring.py`
are the frozen parity reference and are **read, never edited** â nor to `aleph/runtime/parity.py`,
`backend.py`, `warp_kernels.py`, `residency*.py`, `aleph/scenarios/**`, `aleph/viz/**`, or
`ports/ledger/INDEX.md`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) â cited for the record only |
| Source path / symbol | **none named, none read by this lane.** No file under that tree was opened. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable â nothing was read |

The laws being ported are Aleph's own:
`aleph.vertical.connectors_fluid.ImmersedDragConnector.evaluate_sites`,
`aleph.vertical.connectors_transfer.AnisotropicImmersedDragConnector.{drag_tensors_pn_s_per_um,
evaluate_sites,dissipation_pn_um_per_s}`,
`aleph.vertical.connectors_transfer.ExteriorStokesTractionConnector.{evaluate_sites,
dissipated_power_pn_um_per_s}` and
`aleph.vertical.medium.ExteriorStokesMedium.{_apply_resistance,resistance_force_pn}`. This entry
moves already-ported laws onto kernels; it re-ports nothing.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported from another project.** As in `-3603`, `-3604` and `-3606`:
the whole purpose of a parity gate is agreement with the frozen module, so the association order of
every product is fixed by that module. One instance here is load-bearing rather than stylistic and
is called out in Â§4.2: `drag_tensors_pn_s_per_um` writes the tensor as an **isotropic base plus an
anisotropic excess** rather than as `Î¶â¥ tt^T + Î¶â¥(I â tt^T)`, and its own docstring says why â at a
ratio of exactly one the excess is exactly `0.0` and the transverse response is an *algebraic* zero
rather than a residual of order 1e-17. A clean-room kernel would almost certainly write the other
form, and it would destroy the exact reduction this entry's leverage claim rests on.

## 4. Physical or mathematical law represented

### 4.0 The scope decision, the criterion, and the two counts

Task G5b is "the remaining wired connectors". The unit is again **not** position in the wiring
table. The criterion is:

> **Every wired connector whose force is a linear resistance in a velocity** â that is, whose
> `evaluate_sites` returns `energy_start_pn_um == energy_end_pn_um == 0.0` *structurally*, and
> whose two force blocks satisfy `f(Î±V) = Î± f(V)` and `f(0) == 0` **exactly**.

That criterion is measurable rather than editorial, and
`test_every_member_of_the_batch_is_a_linear_resistance_with_no_stored_energy` asserts it by
*evaluating* each member through its own frozen builder: the energy fields with `==`, homogeneity
and additivity to round-off, and `f(0) == 0.0` with `==`. It is disjoint by construction from
`-3606`'s batch, whose criterion was `site_pair_forces` with a `UnilateralSpring`, and
`test_the_two_connector_batches_are_disjoint` asserts that against `CENTRAL_PAIR_BATCH`.

It selects **three classes carrying five of the thirty wired registry rows**:

| class | module | operator `Z` | rows | kernel family |
|---|---|---|---|---|
| `AnisotropicImmersedDragConnector` | `connectors_transfer` | block-diagonal, `3Ã3` per site | **2** | **D (new)** |
| `ImmersedDragConnector` | `connectors_fluid` | block-diagonal, `Î¶ I` | **2** | *D, by an **exact** reduction â see Â§4.3* |
| `ExteriorStokesTractionConnector` | `connectors_transfer` | **dense**, `3NÃ3N` | **1** | **E (new)** |

> ### â  ROWS AND KERNELS ARE DIFFERENT NUMBERS AND QUOTING EITHER ALONE MISREPRESENTS THE PORT
>
> `wiring.py` counts **rows**. This entry writes **kernels**. The honest report is
> **5 registry rows, 3 classes, 2 new kernel families, 12 new kernels (3 true, 9 wrong on purpose,
> 1 shared reduction), and 2 rows carried with no kernel of their own** â
> `nmii_cytosol_transfer` and `sf_cytosol_transfer`, which family D covers by the bitwise reduction
> of Â§4.3 rather than by a second implementation.
>
> `test_the_batch_covers_five_wired_registry_rows_with_two_new_kernel_families` asserts both
> numbers against `wiring.WIRED` directly.

> ### â  THE BRIEF'S CENSUS WAS WRONG IN THREE PLACES, AND THE CORRECTIONS ARE MEASURED
>
> The G5b brief states "12 of 30 wired registry rows are ported" and "about 18 rows remain", and
> lists `BidirectionalLinker` and `StrainStiffeningCableLink` as unported members of
> `connectors_frame.py`, and `AnisotropicImmersedDragConnector` as serving "three cytoskeletal
> cytosol transfers". Measured against `wiring.WIRED` and `law_cases.CENTRAL_PAIR_BATCH`:
>
> | claim | brief | measured |
> |---|---|---|
> | rows already ported | 12 | **14** (`-3604` two + `-3606` twelve, no overlap) |
> | rows remaining | ~18 | **16** |
> | `BidirectionalLinker` | unported | **ported** â it is a member of `CENTRAL_PAIR_BATCH` with `kernel_family: "crosslink"`, and `-3606` Â§4.0's table names its two rows |
> | `AnisotropicImmersedDragConnector` rows | 3 | **2** (`if_cytosol_transfer`, `mt_cytosol_transfer`; `nmii_` and `sf_cytosol_transfer` are `ImmersedDragConnector`) |
>
> Recorded rather than silently used. `-3606` found the build plan had undercounted one class by
> two rows and reported it the same way; this is the same defect class arriving in a brief instead
> of a plan, and the answer is the same â **assert the census against `wiring.WIRED`, never against
> a sentence.**

**After this entry: 19 of 30 wired rows are on kernels, and 11 remain** â
`SeriesJointConnector` (4), `NmiiMotorConnector` (4), `BrownianRatchetConnector` (2),
`StrainStiffeningCableLink` (1). `test_the_ported_row_census_matches_the_wiring_registry` asserts
that total across both batches.

### 4.1 The block-diagonal resistance â family D

Per site `i`, with unit tangent `t_i`, quadrature length `L_i` [Âµm], parallel drag `Î¶â¥`
[pNÂ·s/ÂµmÂ²] and anisotropy ratio `Ï = Î¶â¥/Î¶â¥ â¥ 1`:

```
Î    = Î¶â¥ â Î¶â¥ = (Ï â 1) Î¶â¥                                  â¥ 0
Z_i  = L_i [ Î¶â¥ I + Î (I â t_i t_iáµ) ]                       [pNÂ·s/Âµm]
f_i^(a) = â Z_i v_i                                          [pN]
f_i^(b) = â Z_i (â v_i)                                      [pN]
D       = Î£_i v_i Â· Z_i v_i                                  [pNÂ·Âµm/s]
U       â¡ 0     at both ends of the step, exactly
```

`Z_i / L_i` is transversely isotropic about `t_i` by construction: eigenvalue `Î¶â¥` once, along
`t_i`, and `Î¶â¥` twice in the plane normal to it. It is symmetric positive definite for `Î¶â¥ > 0` and
`Ï â¥ 1`, so `D > 0` for every `v â  0` â which is the statement that makes the element **passive**,
and the one the energy ledger accepts it on.

### 4.2 The dense resistance â family E

```
f^(a) = â R V,     R â â^{3NÃ3N} symmetric positive definite
f^(b) = â R (â V)
P     = â Î£_i f_i^(a) Â· V_i = Váµ R V                         [pNÂ·Âµm/s]
U     â¡ 0                                                     exactly
```

**`R` is not ported and this entry says so rather than implying otherwise.** `R = L^{-T} L^{-1}`
from the Cholesky factor of the regularised-Stokeslet mobility, built by
`ExteriorStokesMedium._rebuild_mobility`. That is a **factorisation**, not a force law: it is where
the medium *proves* its mobility is positive definite, and moving a Cholesky onto a kernel would
move a refusal onto a device that cannot raise. The kernel receives `R` as a device array and
performs the matvec, which is the whole of the connector's own arithmetic. `-3606` recorded the
analogous decision for `k_pair_offset_*` reuse; this is the mirror of it â a piece of the pipeline
**declared as staying on the host**.

### 4.3 The reduction that carries two extra rows with no extra kernel, and why it is *bitwise*

`ImmersedDragConnector` computes `f^(a) = âÎ¶ v` with a scalar. Family D at `Ï = 1` and `L_i = 1`
computes `Z_i = 1Â·(Î¶ I + 0Â·(I â t t^T)) = Î¶ I` and then a matvec. The two are **bit-identical**,
not merely equal to round-off, and the reason is the form the frozen reference chose:

* `excess = Î¶â¥ â Î¶â¥` is exactly `0.0` when `Ï` is exactly `1.0` â the reference's own docstring
  says this is why it writes `Î¶â¥ I + Î(I â ttáµ)` and not `Î¶â¥ ttáµ + Î¶â¥(I â ttáµ)`;
* the off-diagonal entries of `Z_i` are then exactly `0.0`, and `0.0 Â· v_k` is exactly `0.0` for
  every finite `v_k`;
* `x + 0.0 + 0.0 == x` exactly in IEEE-754 for every finite `x`.

**Measured before this entry was written**: 2,000 random draws in float64,
`max |(âÎ¶ v) â (âZ v)| = 0.0`. The kernel-side version of the same statement is asserted by
`test_the_isotropic_drag_is_the_bitwise_reduction_of_the_tensor_kernel`, which compares the two
*kernel* force fields with `np.array_equal` â the number, not a tolerance.

**A mutant that broke the exact form would be caught by that control and by nothing else**, because
every parity budget in this entry is wide enough to absorb a 1e-17 residual. That is the whole
reason the control asserts bit-equality.

### 4.4 What is structurally different about this family, and it is not presentational

**1. THERE IS NO ENERGY CHANNEL.** Every law ported in G1âG6 was a potential, and every lane's
sharpest finding was about grading energy separately from force. Here `U â¡ 0` at both ends of the
step, *structurally* â there is no potential, so there is no branch for stored energy to leak into.
The energy field is therefore an **exact zero over the whole law**, and the parity gate grades it
at a declared budget of literally `0.0` ULP. `-3604` could do that on one *branch* (an all-separated
contact); this is the first law in the port where it is true of every configuration.

**2. WHAT REPLACES IT IS A RATE, NOT AN ENERGY, AND THE DIFFERENCE HAS A MUTANT.** The physical
content of a purely dissipative element is `D = Î£ vÂ·ZÂ·v` [pNÂ·Âµm/s]. It is not conserved, it is not
a potential, and it must be **strictly positive**. `medium.py` ships `dissipation_function` (`Â½VáµRV`,
the potential the resistance is the gradient of *in velocity space*) and `dissipated_power`
(`VáµRV`) as two separate methods precisely because the factor of two gets confused. So this family's
analogue of `WRONG_ENERGY_HALF` is **`WRONG_RAYLEIGH_HALF`** â the dissipation *function* returned
where the *power* is wanted â and, exactly like `WRONG_ENERGY_HALF`, **its force fields are
bit-identical to the true kernel** and it is caught on the dissipation channel alone.

**3. THE POSITION-PRECISION DECISION DOES NOT REACH THIS FAMILY AT ALL.** Both connectors are handed
two position blocks and **verify** them; neither differences them for the law. The input is a
*velocity*, which is `O(1)` and carries no origin, so the ~300-fold cancellation that has set the
achievable agreement of every figure since `ALEPH-PORT-3601` Â§10 **does not exist here**. This entry
adds **no** kernel that reads a position â the module's count of three position-reading families
(`k_face_edges_*`, `k_pair_offset_*`, `k_volume_face_*`) is unchanged, and
`test_the_position_precision_is_inert_for_the_dissipative_family` asserts the force fields are
**bitwise identical** under `GLOBAL_F32`, `LOCAL_F32` and `POSITIONS_F64`.

> That is a contribution to
> `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md`, and it is stated here rather
> than written into the proposal, per `-3605` Â§13.4's boundary. G4 measured `POSITIONS_F64` **worse**
> than `LOCAL_F32`; G5 measured it **48Ã/191Ã better**; both agreed the discriminator is whether a
> law's cancellation is a **subtraction** or a **reduction**. This family adds a third category the
> question's option set does not currently have a slot for: **a law with neither**, for which every
> option is bitwise the same and the decision is inert. A per-law rule (option C) must be able to say
> "this law does not care", or it will assign a mode to five registry rows on no evidence.

**4. THE NEWTON-PAIR MUTANT IS PROVABLY UNKILLABLE FOR ALL THREE CLASSES, BY A STRONGER ARGUMENT
THAN `-3604`'s.** See Â§9a. **No mirror kernel is shipped**, on `-3606`'s precedent.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `Î¶â¥`, `Î¶â¥` | pNÂ·s/ÂµmÂ² | 10â»â¶ PaÂ·s | finite, `> 0` |
| `Ï = Î¶â¥/Î¶â¥` | â | â | finite, `â¥ 1` |
| `L_i` (quadrature length) | Âµm | 10â»â¶ m | finite, `> 0` |
| `t_i` (tangent) | â | â | `\|t_i\| = 1` to 1e-12 |
| `Î¶` (isotropic) | pNÂ·s/Âµm | 10â»Â¹Â² NÂ·s/m | finite, `> 0` |
| `v_i`, `V_i` (relative / surface velocity) | Âµm/s | 10â»â¶ m/s | finite |
| `Z_i`, `R` | pNÂ·s/Âµm | 10â»Â¹Â² NÂ·s/m | symmetric positive definite |
| `f^(a)`, `f^(b)` | pN | 10â»Â¹Â² N | finite |
| `D`, `P` | pNÂ·Âµm/s | 10â»Â¹â¸ W | `> 0` unless `V = 0` |
| `U` | pNÂ·Âµm | 10â»Â¹â¸ J | **exactly `0.0`** |

Singular and boundary cases:

- **`Ï = 1` exactly.** Not singular: `Î = 0.0` exactly and the tensor is exactly `Î¶â¥ L I`. This is
  the case that carries two extra registry rows (Â§4.3) and it is asserted with `==`.
- **`v = 0` exactly.** The force is exactly `0.0` and `D` is exactly `0.0` â no branch intervenes,
  so both are testable with `==`. `medium.py`'s `resistance_force_pn` docstring makes the same
  promise for the dense operator and this entry asserts it through the kernel.
- **A non-unit tangent.** *Refused*, never normalised: `|t| â  1` silently rescales `Î¶â¥` by `|t|Â²`
  and leaves `Î¶â¥` alone, which is an anisotropy nobody asked for. **A Warp kernel cannot raise**, so
  the refusal happens on the host in float64 â and it happens on the **launch** path, not only in
  the state builder, because a caller who assembles their own state mapping would otherwise bypass
  it and receive a plausible-looking wrong tensor. `refuse_non_unit_tangents` is that function.
  `-3606` Â§5 records the identical hole being real in its first draft.
- **A zero or negative quadrature length.** Refused on the host, by the same function, for the same
  reason: a site standing for zero length transmits nothing and would pass every closure check
  vacuously.
- **A non-positive-definite mobility.** Refused by `ExteriorStokesMedium` itself, on the host,
  before this entry sees an `R` at all. Not re-implemented here.

Invariants:

- **I1.** `Z_i / L_i` has eigenvalues `{Î¶â¥, Î¶â¥, Î¶â¥}` â the tensor's *structure*, independent of
  either coefficient's value. Asserted on the **kernel's** tensor by
  `test_the_kernel_drag_tensor_is_transversely_isotropic_about_the_tangent`.
- **I2.** `U == 0.0` exactly, at both ends of the step, for every member. Budget `0.0` ULP.
- **I3.** `f(0) == 0.0` exactly, per site, for every member.
- **I4.** `D > 0` for every `v â  0` (passivity), and `D == 0.0` exactly for `v = 0`.
- **I5.** `D == â Î£_i v_i Â· f_i^(a)` to round-off â the *closure* identity between the two channels,
  which is what a drag element has instead of `F = ââU`. It is an oracle, not a parity field, and it
  is what `WRONG_SIGN` violates.
- **I6.** `f^(b) == â f^(a)` **bitwise**, for all three classes. This is the blind spot, not the
  oracle: see Â§9a.
- **I7.** The reduction of Â§4.3, asserted with `np.array_equal` on the kernel force fields.

## 6. Source evidence class and known retractions

No `ffn_cellsim` artefact is cited and none was read; Â§2 records that. The Aleph-side provenance:

* `ImmersedDragConnector` â `ALEPH-PORT-3201`-class work by the connector lane; live controls in
  `tests/vertical/test_connectors_fluid.py`, which exercise the law (dissipation positivity, the
  reaction actually reaching the fluid, the refusal on an unset velocity) and not only the plumbing.
* `AnisotropicImmersedDragConnector` and `ExteriorStokesTractionConnector` â `connectors_transfer.py`
  with live controls in `tests/vertical/test_connectors_transfer.py`, including an
  **eigen-decomposition** of the drag tensor rather than an inspection of the expression that built
  it. That control is the reason I1 above can be asserted on the kernel: the reference's own oracle
  is structural.
* `ExteriorStokesMedium` â `aleph/vertical/medium.py`, with the positive-definiteness refusal and
  the six rigid-body modes.

**Two self-declared limits in the frozen reference, carried forward rather than repaired:**

1. `connectors_transfer.py`'s module docstring records that `mt_cytosol_transfer` and
   `if_cytosol_transfer` each declare **two** couplings on one edge â immersed drag *and* a
   barbed-end monomer sink â and that **only the drag half is wired**. So a kernel for this class
   covers the wired half of a two-part contract. **Stated, not implied**: this entry's coverage
   claim is about `wiring.WIRED` rows, and those rows are wired to the drag.
2. `nmii_cytosol_transfer`'s builder records that the connector **does not read the duty ratio**, so
   nothing in it is evidence that the unbound motor population is being dragged. Unchanged by this
   port and repeated here so a reader of the coverage table does not take it for more.

No retraction was found attached to any of these symbols. Looked in: each module's docstring, the
`ALEPH-PORT-32xx`/`33xx` ledgers, `HANDOFF.md` Â§F, and `docs/design/CONNECTOR_WIRING_MAP.md`.

## 7. Independent oracle or derivation

Five oracles, none of which is "agreement with the reference":

* **O1 â the eigen-structure.** `Z_i/L_i` must have eigenvalues `{Î¶â¥, Î¶â¥, Î¶â¥}` with the parallel
  eigenvector `Â± t_i`. Measured by `numpy.linalg.eigh` on the **kernel's** assembled tensor, so it
  checks the kernel's algebra against the definition of transverse isotropy rather than against the
  reference's expression.
* **O2 â resolution independence.** Splitting one quadrature site of length `L` into two of length
  `L/2` at the same tangent, position and velocity leaves the *total* drag unchanged, exactly to
  round-off. This is what the `L_i` factor is *for*, and it is the oracle that
  `WRONG_UNWEIGHTED` fails.
* **O3 â passivity and the closure identity (I4, I5).** `D > 0`, and `D = âÎ£ vÂ·f^(a)`.
  Two independently computed quantities that must agree; `WRONG_SIGN` breaks the second while
  leaving the first's *magnitude* alone.
* **O4 â the exact zeros (I2, I3).** `U == 0.0` and `f(0) == 0.0`, asserted with `==`, per site.
* **O5 â Stokes' law for a translating sphere.** The exterior operator's total resistance to a
  **uniform translation** of a sphere of effective radius `a` approaches `6ÏÎ¼a` â a closed-form
  analytic result the discretisation does not contain. Measured through the *kernel*, so the dense
  matvec is checked against physics and not against `medium.py`. **Measured, at three resolutions:
  3.760e-02 relative at `N = 24`, 2.567e-02 at 48, 1.759e-02 at 96 â an oracle with a convergence,
  not an identity.**

> ### â  AND O5 SCORES THE **WRONG** KERNEL FOURTEEN ORDERS OF MAGNITUDE BETTER THAN THE RIGHT ONE
>
> `medium.py`'s `isotropic_node_drag_pn_s_per_um` distributes exactly `6ÏÎ¼a` over the nodes **by
> construction** â its own docstring calls that "the whole trick, and the whole trap".
> `DenseResistanceVariant.WRONG_ISOTROPIC_NODE_DRAG` applies exactly that lumped drag in place of
> the operator, and on a uniform translation its resultant reproduces Stokes' law to **4.2e-16 at
> every resolution** while the true operator is 3.8% off and falling.
>
> **So a closed-form oracle, read as "how close is the resultant to `6ÏÎ¼a`", certifies the mutant
> and fails the true kernel â at every resolution, by 14 orders of magnitude in float64 and 6 in
> float32.** What separates them is the **convergence rate**: the true operator's error falls with
> `N`, the impostor's does not move. A single-resolution accuracy comparison cannot see that, and
> the parity gate â which grades the per-node force *field* against the frozen reference â catches
> the impostor at 2.2e6 ULP without needing to.
>
> **An analytic oracle is not automatically a stronger check than parity. It is a different check,
> and a discretisation that satisfies the closed form by construction is unfalsifiable by it.**
> `test_the_analytic_stokes_oracle_prefers_the_impostor_to_the_true_kernel` asserts both halves,
> and `test_the_true_dense_kernel_converges_to_stokes_law_and_the_impostor_does_not` asserts the
> discriminator that does work.
>
> **This corrects a claim in an earlier draft of this entry.** That draft attributed the
> resultant-blindness to `WRONG_NODE_DIAGONAL` â retaining only `R`'s node-diagonal blocks â which
> is a **different object** and is **not** resultant-blind: measured, its resultant is off by a
> factor of **6.4**. The claim was written from the reference's docstring rather than from a
> measurement, and it was wrong. Corrected here and in the kernel's own comment rather than quietly
> repaired, because a blind spot attributed to the wrong mutant is a blind spot nobody re-measures.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| P1 | `tests/runtime/test_transfer_law_parity.py::test_the_true_anisotropic_drag_kernel_agrees_with_the_frozen_numpy_law` | Every field of the anisotropic case inside `DEFAULT_BLOCK_DIAGONAL_DRAG_ULP_BUDGET`, on warp's CPU device, float32, `ORDERED`. |
| P2 | `â¦::test_the_true_isotropic_drag_kernel_agrees_with_the_frozen_numpy_law` | The same for the isotropic case, graded against `ImmersedDragConnector`'s own reference. |
| P3 | `â¦::test_the_true_dense_resistance_kernel_agrees_with_the_frozen_numpy_law` | The same for the exterior traction, against a real `ExteriorStokesMedium`. |
| P4 | `â¦::test_a_transfer_law_compared_with_itself_is_exactly_zero_ulp` | The harness's zero point: `0.00` ULP on every field of all three cases. |
| P5 | `â¦::test_the_isotropic_drag_is_the_bitwise_reduction_of_the_tensor_kernel` | Â§4.3: `np.array_equal` on the two kernels' force fields. The leverage claim, as a number. |
| P6 | `â¦::test_the_kernel_drag_tensor_is_transversely_isotropic_about_the_tangent` | O1, by `eigh` on the kernel's tensor. |
| P7 | `â¦::test_the_kernel_drag_is_resolution_independent_under_site_splitting` | O2. |
| P8 | `â¦::test_the_kernel_dissipation_equals_minus_the_power_of_its_own_force_field` | I5/O3, both families. |
| P9 | `â¦::test_the_kernel_exterior_resistance_recovers_stokes_law_for_a_translating_sphere` | O5, through the kernel. |
| P10 | `tests/runtime/test_connector_law_parity.py::test_a_zero_ulp_budget_fails_the_true_connector_kernels` | I2 at budget `0.0`. |
| P11 | `â¦::test_a_zero_velocity_carries_exactly_zero_force_and_zero_dissipation` | I3/I4, with `==`. |
| P12 | `â¦::test_the_position_precision_is_inert_for_the_dissipative_family` | Â§4.4(3), bitwise across all three modes. |
| P13 | `â¦::test_the_batch_covers_five_wired_registry_rows_with_two_new_kernel_families` | The census, against `wiring.WIRED`. |
| P14 | `â¦::test_the_ported_row_census_matches_the_wiring_registry` | 19 of 30 ported, 11 remaining, both batches, against `wiring.WIRED`. |

## 9. Deliberately failing negative control

Nine wrong kernels ship in the package, compiled by the same `load()` and launched down the same
driver. **Never select a `WRONG_*` variant in a run.**

| Variant | What it is wrong about | Caught on | Control |
|---|---|---|---|
| `BlockDiagonalDragVariant.WRONG_ISOTROPIC` | `Î` deleted: `Z = L Î¶â¥ I`. The defect `connectors_transfer.py`'s docstring names in as many words â "an isotropic Î¶ deletes thatâ¦ no transverse component and no reorienting torque". | `forces_a_pn`, `forces_b_pn`, `dissipation_pn_um_per_s` | `test_a_drag_kernel_that_lost_its_anisotropy_is_caught` |
| `â¦WRONG_AXIS_SWAPPED` | `Z = L(Î¶â¥ I + Î ttáµ)`: the high-drag axis put **along** the filament. The reference's `__post_init__` calls this "not a limiting case of a filament, it is a sign error". | `forces_a_pn`, `forces_b_pn`, `dissipation_pn_um_per_s` | `test_a_drag_kernel_with_its_drag_axes_swapped_is_caught` |
| `â¦WRONG_UNWEIGHTED` | `L_i` dropped, so the drag integral stops being resolution-independent. | `forces_*`, `dissipation` | `test_a_drag_kernel_that_dropped_the_quadrature_length_is_caught` |
| `â¦WRONG_SIGN` | `f = +Z v`: the medium becomes a source of power at every site. | `forces_*` **only** â see Â§9b | `test_a_drag_kernel_whose_sign_makes_the_medium_a_power_source_is_caught` |
| `DragDissipationVariant.WRONG_RAYLEIGH_HALF` | `Â½ vÂ·ZÂ·v`, the Rayleigh *function* where the *power* is wanted. **Force fields bit-identical.** | `dissipation_pn_um_per_s` **alone** | `test_a_dissipation_kernel_that_returned_the_rayleigh_function_is_caught` |
| `DenseResistanceVariant.WRONG_NODE_DIAGONAL` | Only each node's own `3×3` block of `R` retained: the hydrodynamic interaction deleted. **Its resultant is 6.4× wrong too** — an earlier draft of this table claimed it was resultant-blind, and that claim belonged to the row below. | `forces_*`, `dissipation`, and the resultant | `test_a_dense_kernel_that_kept_only_the_node_diagonal_is_caught` |
| `…WRONG_ISOTROPIC_NODE_DRAG` | `f_i = -γ_i V_i` with `γ_i = 6πμa w_i/A` — `medium.py`'s own shipped break, whose coefficients sum to `6πμa` **by construction**. **On a uniform translation it satisfies Stokes' law exactly at every resolution while the true operator is 3.8% off and converging**, so O5 ranks it first. Takes an extra kernel argument the true law provably does not need. | `forces_*`, `dissipation` — **and NOT O5**, §7 | `test_the_analytic_stokes_oracle_prefers_the_impostor_to_the_true_kernel`, `test_the_true_dense_kernel_converges_to_stokes_law_and_the_impostor_does_not` |
| `â¦WRONG_COMPONENT_STRIDE` | `R` indexed component-major (`cÂ·N + i`) where it is node-major (`3i + c`). The likeliest real slip in a `(3N, 3N)` flatten. | `forces_*`, `dissipation` | `test_a_dense_kernel_that_flattened_the_operator_the_wrong_way_is_caught` |
| `â¦WRONG_SIGN` | as above, for the dense operator. | `forces_*` **and** `dissipation` â see Â§9b | `test_a_dense_kernel_whose_sign_makes_the_medium_a_power_source_is_caught` |

Plus the budget-vacuity control `test_a_zero_ulp_budget_fails_the_true_transfer_kernels` (a gate
that cannot fail is not a gate), the atomic-floor control
`test_no_transfer_mutant_escapes_through_the_atomic_floor`, and the round-off-explanation control
`test_the_declared_transfer_budgets_sit_inside_what_roundoff_explains`.

### 9a. The Newton-pair control that is NOT shipped, and why the argument is stronger here

`HANDOFF.md` Â§F-1 and `-3604` Â§9a established that `force_b := âforce_a` is unkillable in an
*exactly linear* law. `-3604`'s case was a central law with a scalar magnitude and the evidence was
2,000 draws. Here the argument closes:

> **The whole family is exactly linear in the velocity, and `f^(b)` is the same operator evaluated
> at `âV`.** In IEEE-754 with any rounding mode that is symmetric about zero â which round-to-nearest
> is â negation is exact and commutes with every operation in the chain: `fl(âx) = âfl(x)`,
> `fl((âa)Â·b) = âfl(aÂ·b)`, `fl((âp) + (âq)) = âfl(p + q)`. So `Z(âv)` is the **exact bitwise
> negation** of `Z v` for *any* evaluation order, any blocking, any FMA contraction, on any device.
> `f^(b) = âZ(âv) = +Z v = âf^(a)` **bitwise**, and no measurement of any kind can separate the
> mutant from the true law.

This is not a claim that the reference's care is pointless: the reference writes the two evaluations
out precisely so that a resistance which *became* state-dependent or asymmetric would separate them,
and this argument is conditional on linearity, which is a property a future edit can remove.

**Measured anyway, because a proof about arithmetic is worth one measurement**: 2,000 draws per
family, `max |f^(b) + f^(a)| == 0.0` exactly, asserted by
`test_the_newton_mirror_mutant_is_provably_unkillable_for_the_whole_transfer_family` for the tensor
contraction and the dense matvec both. **A `WRONG_NEWTON_MIRROR` kernel is not shipped for this
family**, on `-3606`'s precedent, and this section is the record of why.

### 9b. The sign mutant is caught in different channels in the two families, and that is a finding

`WRONG_SIGN` is the same physical defect in both families â the medium becomes a power source. It is
caught in **different** channels, because the two references compute their dissipation differently
and the port matches each:

| family | how the reference computes the rate | `WRONG_SIGN` visible in `dissipation`? |
|---|---|---|
| **D** (block-diagonal) | `Î£ vÂ·ZÂ·v` â **independently** of the force, from the tensor | **No.** The forces move, the rate does not. |
| **E** (dense) | `â Î£ f^(a)Â·V` â **from the force it delivered** | **Yes.** Both channels move together. |

`medium.py`'s own docstring explains E's choice: *"Computed from the force that was actually
delivered rather than from the medium's own dissipation functionâ¦ reporting a number derived from
one while applying the other is how a work ledger stops describing the forces that were used."*
So D's independence buys channel separation and E's dependence buys agreement â **the two are a
genuine trade and neither is a defect.** `test_the_sign_mutant_is_blind_in_the_block_diagonal_rate_and_visible_in_the_dense_one`
asserts both halves, which is what keeps Â§9's "caught on" column from being a guess.

## 10. Numerical and precision envelope

Working precision: **float32** for the tensor, the velocity, the matvec and the force â the compute
channel `ALEPH-DQ-107` ratified. Accumulation: **float64**, per site, before `Backend.sum`'s float64
chunked reduction, for the energy (identically zero) and the dissipation partials.

**This family's conditioning is the best in the port, and the reason is structural.** Every figure
from G1 onward has been dominated by float32 *position storage*: a short element between surfaces
5 Âµm from the origin cancels ~300-fold before any kernel runs. **There is no such cancellation
here.** The law's input is a velocity of order 0.1â1 Âµm/s with no origin, the tensor is a bounded
`O(Î¶ L)` quantity, and the only reductions are a `3`-term matvec (family D) or a `3N`-term row sum
(family E). So:

* **family D's amplification is `O(1)`** â the derived round-off bound is `(amp + 1)/2` with
  `amp â 3` from the three-term contraction, and the budget is set from *that* rather than inherited;
* **family E's amplification is `O(â(3N))` at worst** â a `3N`-term dot product with terms of mixed
  sign. The budget scales with the node count and the case's own state carries `N`, so a larger
  sphere does not silently widen the gate.

**No budget in this entry is inherited from a float64 tolerance.** `-3605` Â§14 recorded the trap
(the frozen `1e-9 Ã peak` is a float64 property of a float64 path, eight orders below float32's
reach). `test_the_declared_transfer_budgets_sit_inside_what_roundoff_explains` refuses any budget
above what `roundoff_ulp_bound` explains, so a number chosen after seeing the answer is a red test.

Outside the envelope: a non-unit tangent, a non-positive quadrature length, a non-finite velocity
and a non-positive-definite mobility are **refused on the host, in float64, before any upload** â
never silently degraded, because a kernel cannot raise.

## 11. Production-backend residency and transfer

Per evaluation, family D uploads the tangents `(N,3)`, the quadrature lengths `(N,)` and the
velocity `(N,3)` as float32, plus two float32 scalars, and reads back two `(N,3)` force fields and
two float64 reductions. Family E additionally uploads `R` as `(3N, 3N)` float32 â **once per
mobility rebuild, not per evaluation**, because `R` changes only when the surface moves, and this
entry's launch takes `R` as an already-resident array so a caller can keep it there.

**This entry does not make anything resident.** `ALEPH-PORT-3609` owns that layer, and the batching
of `R` across a descent is exactly the kind of thing that belongs there rather than here. Stated as
an absence: a caller who rebuilds `R` every step will re-upload `(3N)Â²` float32 every step, and
nothing in this entry stops them.

**No relaxation is put on device by this entry.** `-3603` Â§"the second precision channel" and
`-3609` Â§14.6 are unresolved and above a porting lane.

## 12. Comments and docstrings to discard

Nothing to discard: no source prose crossed a boundary, because no source file was read (Â§2).

What replaces the absent provenance, so the new code is not unattributed: every kernel names the
frozen Aleph symbol it reproduces and the association order it matches; every mutant names the
sentence in the *reference's own* docstring that describes the defect it plants, since in three of
the eight cases the reference explicitly warns a reader against exactly that error
(`WRONG_ISOTROPIC`, `WRONG_AXIS_SWAPPED`, `WRONG_RAYLEIGH_HALF`); and Â§9a records a control that is
deliberately **absent** with the argument for its absence, rather than shipping a control that
cannot fail.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | 2026-08-01, warp CPU device, float32, `ORDERED`. `tests/runtime` **505 passed / 3 skipped / 1 failed** (the failure is `-3608`'s guard, Â§14.11); `tests/runtime` + `tests/vertical` **1,882 / 3 / 1**. `aleph/vertical/**` byte-for-byte untouched. **17 planted mutants, 17 killed**, plus one equivalent mutant recorded. Every field of every true case inside its declared budget, and every budget inside what round-off explains. |
| Reviewer | **Agent-proposed. Unratified.** No `decided_by` field appears anywhere in this entry. |
| Rollback | Delete `BlockDiagonalDragVariant`, `DragDissipationVariant`, `DenseResistanceVariant` and their twelve kernels from `law_kernels.py`; delete the `DISSIPATIVE_TRANSFER_BATCH` block and the three case builders from `law_cases.py`; delete `tests/runtime/test_transfer_law_parity.py`. **Nothing else breaks**: no existing line of either module is changed by this entry, `aleph/vertical/**` is untouched, and the five registry rows return to being NumPy-only, which is what they are today. |

### 13.1 Measured â warp CPU device, float32, `ORDERED`, 2026-08-01

**The census, against `wiring.WIRED`, as two numbers:**

| | count |
|---|---|
| wired registry rows **covered by this entry** | **5** |
| implementing classes | **3** |
| new kernel **families** | **2** |
| new **kernels** | **12** (3 true, 9 wrong on purpose — counted off the `_define()` dispatch table, not off this sentence) |
| rows carried with **no kernel of their own** | **2** (`nmii_`/`sf_cytosol_transfer`, by the bitwise reduction) |
| wired rows on kernels **after** this entry | **19 of 30** |

**The ULP table, with the catching channel per wrong kernel.** Budgets: family D `forces 2.0`,
`dissipation 2.0`, `energy 0.0`; family E `forces 16.0`, `dissipation 16.0`, `energy 0.0`.

#### Family D â the anisotropic state (96 sites: 16 parallel, 6 transverse, 74 oblique)

| kernel | `forces_a` | `forces_b` | `dissipation` | `energy` | caught by |
|---|---|---|---|---|---|
| **`bdiag_drag_true`** | **0.64** | **0.64** | **0.03** | **0.00** | â passes |
| `bdiag_drag_wrong_isotropic` | 2,782,579.71 | 2,782,579.71 | 2,358,261.66 | 0.00 | all three |
| `bdiag_drag_wrong_axis_swapped` | 4,728,935.54 | 4,728,935.54 | 193,763.52 | 0.00 | all three **+ the eigenvector of O1** |
| `bdiag_drag_wrong_unweighted` | 3,138,982.50 | 3,138,982.50 | 1,894,520.55 | 0.00 | all three **+ the splitting oracle O2** |
| `bdiag_drag_wrong_sign` | 16,777,216.40 | 16,777,216.40 | *bit-identical* | 0.00 | **force alone** |
| `drag_dissipation_wrong_rayleigh_half` | *bit-identical* | *bit-identical* | 4,194,303.98 | 0.00 | **dissipation alone** |

> **Two mutants, opposite blindnesses, in the same law.** `WRONG_SIGN` leaves the rate bit-identical
> and `WRONG_RAYLEIGH_HALF` leaves both force fields bit-identical. Neither channel alone grades
> this family, which is G2's "grade energy and force separately" arriving in a law that **has no
> energy** â the channels are force and *rate*, and the argument is identical.

#### Family D â the isotropic state (`ratio = 1`, `L = 1`): three declared blind spots

| kernel | `forces_a` | `forces_b` | `dissipation` | caught by |
|---|---|---|---|---|
| **`bdiag_drag_true`** | **0.57** | **0.57** | **0.11** | â passes |
| `bdiag_drag_wrong_isotropic` | *bit-identical* | *bit-identical* | *bit-identical* | **nothing â declared blind** |
| `bdiag_drag_wrong_axis_swapped` | *bit-identical* | *bit-identical* | *bit-identical* | **nothing â declared blind** |
| `bdiag_drag_wrong_unweighted` | *bit-identical* | *bit-identical* | *bit-identical* | **nothing â declared blind** |
| `bdiag_drag_wrong_sign` | 16,777,216.47 | 16,777,216.47 | *bit-identical* | force alone |

> **This is `ALEPH-PORT-3601` Â§14.7 in its sharpest form yet: the two `ImmersedDragConnector` rows
> sit on a configuration that cannot grade three of the four force mutants**, because at `ratio = 1`
> and `L = 1` those three *are* the true kernel. That is not a weakness of the kernel â it is the
> statement that the anisotropic state is what grades them â and
> `test_three_drag_mutants_are_invisible_on_the_isotropic_unit_length_configuration` asserts it with
> `np.array_equal` so nobody later reads the isotropic row's pass as covering them.

#### Family E â the dense operator, `N = 24`, `72Ã72`

| kernel | `forces_a` | `forces_b` | `dissipation` | caught by |
|---|---|---|---|---|
| **`dense_resistance_true`** (mixed `V`) | **1.60** | **1.60** | **0.13** | â passes |
| **`dense_resistance_true`** (uniform translation) | **4.00** | **4.00** | **0.06** | â passes |
| `dense_resistance_wrong_node_diagonal` | 2,105,191.01 | 2,105,191.01 | 183,152.07 | all three |
| `dense_resistance_wrong_isotropic_node_drag` | 7,392,991.54 | 7,392,991.54 | 7,096,935.19 | all three â **and O5 *prefers* it, Â§7** |
| `dense_resistance_wrong_component_stride` | 3,282,764.05 | 3,282,764.05 | 560,969.59 | all three |
| `dense_resistance_wrong_sign` | 16,777,214.82 | 16,777,214.82 | 16,777,216.13 | **force *and* rate** â contrast family D, Â§9b |

**The budgets, and that round-off explains them:**

| | worst true kernel, over every loaded state | tightest round-off explanation | budget | margin over / under |
|---|---|---|---|---|
| family D forces | **0.806** | **3.00** | **2.0** | 2.5Ã / 1.5Ã |
| family D rate | **0.253** | **3.00** | **2.0** | 7.9Ã / 1.5Ã |
| family E forces | **3.998** | **61.17** | **16.0** | 4.0Ã / 3.8Ã |
| family E rate | **0.132** | **91.86** | **16.0** | 121Ã / 5.7Ã |

> **These are the tightest force budgets in the port, by two to three orders of magnitude** â 2 and
> 16 against 256 for a crosslink, 2048 for a soft-core contact and ~100 ULP for the first tether â
> **and it is structural, not care.** There is no position difference anywhere in this law, so
> `ALEPH-PORT-3601` Â§10's ~300-fold cancellation has nothing to act on.

**The other measured anchors:**

| claim | measured |
|---|---|
| self-comparison, all three cases, every field | **0.00 ULP** |
| `U` at both ends of the step, every kernel including every mutant | **exactly `0.0`** |
| `f(V = 0)` and `D(V = 0)`, per site, all three classes | **exactly `0.0`** (`np.array_equal`) |
| `f^(b) == âf^(a)` bitwise, through both kernel families | **True** |
| `f^(b) + f^(a)`, 2,000 host draws of the tensor contraction, 4,965 sign-mixed sites | **exactly `0.0`** |
| isotropic reduction vs `âÎ¶ v` in float32 | **bit-identical** (`np.array_equal`) |
| off-diagonal entries of `Z` at `ratio = 1` | **exactly `0.0`** |
| eigenvalues of `Z/L` vs `{Î¶â¥, Î¶â¥, Î¶â¥}` | 1.55e-07 relative |
| `â¨parallel eigenvector, tâ©` | 0.9999999999999994 |
| `D` vs `âÎ£ vÂ·f^(a)` (the closure identity I5) | agrees to 1e-5 relative; `D > 0` at every site |
| splitting oracle O2, total drag under 2Ã refinement | unchanged to 1e-5 relative |
| Stokes' law through the kernel, `N = 24 / 48 / 96` | 3.760e-02 / 2.567e-02 / 1.759e-02 â **converging** |
| the impostor's Stokes error, same resolutions | 5.5e-08 (float32); 4.2e-16 in float64 â **flat** |
| `PositionPrecision` across `GLOBAL_F32` / `LOCAL_F32` / `POSITIONS_F64` | **bitwise identical** |
| `WRONG_NODE_DIAGONAL`'s resultant vs the true one | **6.4Ã wrong** â *not* resultant-blind |

**The mutation study: 17 planted, 17 killed, plus one equivalent mutant recorded.** Nine are the
shipped wrong kernels above, each killed by a named control; `test_a_drag_kernel_that_lost_its_
anisotropy_is_caught` was itself re-pointed at the TRUE variant to confirm the control can fail.
Eight more were planted in the shipped modules, confirmed to redden a named test, and reverted â
`PYTHONDONTWRITEBYTECODE=1`, from a bytecode tree verified empty under `aleph/` and `tests/`:

| # | planted defect | reddens |
|---|---|---|
| M9 | `refuse_non_unit_tangents` removed from `_drag_launch` (left in the state builder only) | `test_a_non_unit_tangent_is_refused_before_any_kernel_launches` |
| M10 | `isotropic_drag_state` built at `ratio = 1 + 2e-16` instead of exactly `1.0` | `test_the_isotropic_drag_is_the_bitwise_reduction_of_the_tensor_kernel` |
| M11 | the exact-zero control's state given `at_rest=False`, emptying its mask | `test_a_zero_velocity_carries_exactly_zero_force_and_zero_dissipation` |
| M12 | `np.array(...)` â `np.asarray(...)` in both per-site helpers (11 sites) â the `-3609` Â§14.5 view hazard | `test_a_pulled_transfer_field_owns_its_memory_and_does_not_view_the_device` |
| M13 | `DEFAULT_BLOCK_DIAGONAL_DRAG_ULP_BUDGET` forces/rate widened 2.0 â 4.0 | `test_the_declared_transfer_budgets_sit_inside_what_roundoff_explains` |
| M14 | the TRUE kernel's **published tensor** scaled by 1.02, forces left correct | `test_the_kernel_drag_tensor_is_transversely_isotropic_about_the_tangent` (and the parity gate, via the rate) |
| M15 | the rate kernel's `tensor` argument fed `site_force_a` â warp does **not** shape-check `array3d` against `array2d` | `test_the_sign_mutant_is_blind_in_the_block_diagonal_rate_and_visible_in_the_dense_one` |
| M16 | `TensileLinker` (a spring) added to `DISSIPATIVE_TRANSFER_BATCH` | `test_the_two_connector_batches_are_disjoint`, `test_the_batch_covers_five_wired_registry_rowsâ¦` |

> ### â  TWO OF THESE FOUND REAL DEFECTS IN THIS LANE'S OWN CONTROLS, AND BOTH ARE NEW SHAPES
>
> **M12 survived its first planting, and the reason is a sixth way a control of G6's shape can be
> silent.** The control was written on `force_a_pn` â a **float32** device array pulled to float64 â
> and the *dtype conversion already copies*, so `np.asarray(view, dtype=np.float64)` and
> `np.array(...)` are indistinguishable there. **A narrowing pull hides the aliasing.** The hazard
> is live only on the fields whose device dtype already matches the host dtype â here the three
> float64 reductions â and a view control written on any other field is *structurally* unable to
> see it. The control was moved to `dissipation_pn_um_per_s` and now checks `owndata`/`base`
> directly; replanted, killed.
>
> **M15 survived its first planting because the control was an equality with no anchor.** It
> asserted "the sign mutant's rate equals the true kernel's rate", which is preserved by any defect
> that corrupts *both* arms identically â and M15 is exactly that. **A relative assertion is only
> evidence once one side is independently anchored**, so the control now requires the true rate to
> be inside its own parity budget first. Replanted, killed. This is `-3606` M8's finding in a third
> place: G5 found a gate cannot see a defect in its reference; here a gate could not see a defect in
> *both* of the things it was comparing.

**One equivalent mutant, recorded rather than counted as a kill or as a survivor.** **M17**
re-associates the tensor form from `excessÂ·(base â t_r t_c)` to `excessÂ·base â excessÂ·t_r t_c`.
`ALEPH-PORT-3610` Â§3 calls the reference's association load-bearing, and this measures exactly how
far that goes: the **exactness at `ratio = 1` survives** (both forms give `0.0Â·base â 0.0 = 0.0`),
so the bitwise-reduction control correctly does not fire, and the anisotropic force ULP moves
**0.643 â 1.381** against a budget of 2.0. So Â§3's claim is true only of the *exactness*, and the
association's effect on round-off consumes 69% of the budget rather than 32%. Recorded because a
future widening of that budget would remove the only pressure on it.

## 14. Honest limits

1. **`UNVERIFIED` on a CUDA device.** Every figure is warp's CPU device. `-3609` measured that a
   CUDA device changes the float32 noise floor enough to move a descent window by a binary order;
   nothing here is licensed to be quoted as a device result. A device run is optional for this lane
   and was not taken.
2. **`R` is not on a kernel and this entry does not claim it is.** Â§4.2. The Cholesky, the
   positive-definiteness refusal and the mobility assembly stay on the host.
3. **Family D's dissipation is graded, and family D's *reference* computes it from the tensor.**
   So the dissipation channel and the force channel share the tensor. A defect in the tensor
   assembly moves both and the two are not independent evidence â which is why Â§9's "caught on"
   column lists both for three of the four D mutants rather than claiming channel separation it does
   not have. Genuine separation exists in exactly two places: `WRONG_RAYLEIGH_HALF` (rate only) and
   `WRONG_SIGN` in family D (force only).
4. **The Newton-pair control is absent by argument, not by measurement alone.** Â§9a. If a future
   edit makes any of these three resistances state-dependent, the argument lapses and a mirror
   control becomes both possible and necessary. Nothing mechanically detects that lapse.
5. **A parity gate still cannot see a defect in the state it grades on (`-3605` Â§13a), in the
   reference it grades against (`-3606` M8), or when its asserted-over set is empty (`-3606` O2).**
   This entry inherits all three. The mitigations it ships are the eigen-structure oracle O1 and
   Stokes' law O5, which are checks on *physics* rather than on the reference, and a non-emptiness
   assertion on every `==` control's mask.
6. **`connectors_transfer.py` records that two of the five rows are half-wired** (drag yes, monomer
   sink no) and that `nmii_cytosol_transfer` does not read the duty ratio. Â§6. The coverage claim is
   about wired rows and inherits both limits unchanged.
7. **11 of 30 wired rows remain unported after this entry** â `SeriesJointConnector` (4),
   `NmiiMotorConnector` (4), `BrownianRatchetConnector` (2), `StrainStiffeningCableLink` (1). None
   is started here and none is blocked by anything in this entry.
8. **The exact-zero energy channel is cheap evidence and is labelled as such.** `U â¡ 0` is asserted
   at budget `0.0`, and it is true of *any* kernel that writes zeros into that field â including
   every mutant here. It is a real property of the law and it grades nothing. The force and
   dissipation channels are the ones that carry the port.
9. **`ports/ledger/INDEX.md` is stale and red before this entry existed**, naming `-3606`/`-3607`/
   `-3608`/`-3609`. `-3610` adds a fifth *name*, not a third red test. **Not regenerated**, on
   instruction; reported here and in `docs/ACTIVE_SESSIONS.md`.

10. **Â§7's blind-spot claim was wrong in the first draft of this entry and is corrected in place.**
    It attributed the resultant-blindness to `WRONG_NODE_DIAGONAL`; the resultant-blind object is
    `WRONG_ISOTROPIC_NODE_DRAG`, and `WRONG_NODE_DIAGONAL`'s resultant is 6.4Ã wrong. The claim was
    written from the reference's docstring rather than from a measurement. Both mutants now ship,
    both are asserted, and the correction is recorded in `law_kernels.py` beside the kernel as well
    as here â a blind spot attributed to the wrong object is a blind spot nobody re-measures.

11. ### â â  A `-3608` GUARD IS RED ON THIS ENTRY'S ISOTROPIC CASE, AND THIS LANE DID NOT TOUCH IT

    `tests/runtime/test_reference_input_guard.py::test_every_registered_case_reference_is_sensitive_to_position_precision`
    requires every registered `LawCase` to move by **more than 1.0 ULP** when its reference's input
    is rounded to float32 â the half of the reference-contamination check that stops it being empty.
    `isotropic_drag_law_case` moves by **0.443 ULP** and fails it.

    **The guard is right and the finding is structural, not a defect in either.**
    `ImmersedDragConnector`'s law is a **single multiplication**, `f = âÎ¶ v`. One float32 multiply
    of a once-rounded input cannot move its result by more than about 1 ULP, at any state, ever. So
    if that case's reference were accidentally handed float32, the gate would print 0.44 ULP against
    a budget of 2.0 and nothing would notice â which is exactly what the guard says.

    Two facts bound how much that matters. The **mechanical** half of the mechanism,
    `parity.refuse_narrowed_reference_input`, is a *dtype* check run inside `compare_law_case`
    before either side runs; it is unconditional and this case cannot evade it. The guard is the
    *numerical visibility* check on a failure the dtype check already makes impossible.

    **Three candidate repairs, and this lane takes none of them** (`CLAUDE.md` Â§2 rules 4 and 5):

    | | repair | why not this lane's |
    |---|---|---|
    | (a) | `-3608` widens the guard to accept a case whose law provably has no amplification â it already accepts a case that *refuses* the rounded input, and one that returns exact zeros | edits `-3608`'s guard |
    | (b) | `-3610` withdraws `isotropic_drag_law_case` from the registry and covers those two rows by the bitwise-reduction control alone | removes the only grading of `ImmersedDragConnector`'s **own** reference, to make a guard green â routing around it |
    | (c) | the threshold moves from `> 1.0` ULP to something a single multiplication can reach | weakens the guard for every other case |

    **Reported, not resolved.** The anisotropic case passes the same guard *by refusal* â rounding
    its tangents to float32 breaks the unit-vector constraint and the frozen reference raises, which
    the guard explicitly accepts as the strongest form of sensitivity â and the dense case passes on
    measured sensitivity. It is one case of three.

12. **The tightest budgets in the port are also the ones with the least headroom.** M17 measured that
    a *legitimate* re-association of the tensor expression moves the anisotropic force figure from
    0.643 to 1.381 ULP against a budget of 2.0. The budget is inside what round-off explains (3.00),
    but a second such change would not be. Named so that a future edit to the tensor form is
    understood as consuming margin rather than as free.

13. **`warp` does not shape-check a `wp.array3d` parameter against a `wp.array2d` argument** at
    launch â discovered by M15, which ran clean and produced meaningless numbers. Nothing in this
    entry depends on it, and no guard here would catch the same slip in a *different* launch.
    Recorded as a property of the substrate rather than repaired, because repairing it belongs to
    whoever owns the launch seam.
