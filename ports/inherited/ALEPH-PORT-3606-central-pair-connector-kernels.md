# ALEPH-PORT-3606 — the central-pair connector family as Warp kernels

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3606` |
| Lane | `46143f30` Lane G5 (CUDA track), Track G task **G5** |
| Status | `PROPOSED` |
| Written | `2026-08-01` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Depends on | `ALEPH-PORT-3601` (the law-level parity harness), `ALEPH-PORT-3603` (`PositionPrecision`), and `ALEPH-PORT-3604` (`k_pair_offset_f32` / `_f64`, the state-branch-count discipline, and the Newton-mirror measurement this entry repeats for two more laws). All three are **reused rather than duplicated**. |
| Exists because | Task G5 is "the remaining 28 wired connectors", and 28 lanes is the wrong unit. `aleph/vertical/wiring.py` maps **30 wired rows onto 13 implementing classes**, so the leverage is in the classes. This entry takes the largest coherent class family and covers **12 rows with 2 new kernel families**. |

---

## 1. Aleph API

```python
from aleph.runtime.law_kernels import (
    CrosslinkVariant,           # TRUE | WRONG_UNILATERAL_CLAMP | WRONG_ENERGY_HALF
                                #      | WRONG_ABS_STRETCH
    CentralContactVariant,      # TRUE | WRONG_LINEAR_CORE | WRONG_ADHESIVE_GATE
                                #      | WRONG_DERIVATIVE_DENOMINATOR   (renamed -- see 9c)
    PositionPrecision,          # unchanged, from ALEPH-PORT-3603
    load, loaded_warp,          # unchanged, extended to compile these kernels
)

from aleph.runtime.law_cases import (
    CROSSLINK_STIFFNESS_PN_PER_UM,      # 50 pN/um
    CROSSLINK_REST_LENGTH_UM,           # 0.050 um
    CENTRAL_CONTACT_STANDOFF_UM,        # 0.020 um
    CENTRAL_CONTACT_STIFFNESS_PN_PER_UM,# 2300 pN/um
    DEFAULT_CROSSLINK_ULP_BUDGET,
    DEFAULT_CENTRAL_CONTACT_ULP_BUDGET,
    CENTRAL_PAIR_BATCH,          # the batch's own census: class -> registry rows -> kernel family
    crosslink_state,             # reproducible two-sided geometry, both senses, both step kinds
    crosslink_rest_length_state, # the origin-anchored geometry the exact anchor NEEDS -- see 13.6
    crosslink_law_case,          # the two-sided linear spring as a LawCase
    refuse_coincident_central_pairs,  # the host refusal, on the LAUNCH path -- see 5
    crosslink_kernel_site_forces,# the kernel's own PER-SITE fields, for the exact anchors
    central_contact_state,       # reproducible compressive line-of-centres geometry
    central_contact_law_case,    # the central SOFT_CORE law as a LawCase
    central_contact_kernel_site_forces,
)
```

Nothing outside this list is covered. In particular this entry does **not** authorise any change to
`aleph/vertical/**` — **`connectors.py`, `connectors_crosslink.py`, `connectors_solid.py` and
`connectors_frame.py` are the frozen parity reference and are read, never edited** — nor to
`aleph/runtime/parity.py`, `backend.py`, `warp_kernels.py`, `aleph/scenarios/**`, `aleph/viz/**`,
or `ports/ledger/INDEX.md`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read by this lane.** No file under that tree was opened. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

The laws being ported are Aleph's own: `aleph.vertical.connectors.site_pair_forces` driven by
`aleph.vertical.connectors.UnilateralSpring`, whose provenance is `ALEPH-PORT-1103`, and the four
connector classes that select its parameters — `ALEPH-PORT-3301` (`connectors_frame`) and the
crosslink/solid entries. This entry moves already-ported laws onto kernels; it re-ports nothing.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported from another project.** As in `-3603` and `-3604`: the whole
purpose of a parity gate is agreement with the frozen module, so the association order of every
product and the exact placement of every clamp are fixed by that module. A clean-room second
implementation would put its `max(·, 0)` somewhere else and would disagree at the gate, which is the
one thing this gate exists to measure.

## 4. Physical or mathematical law represented

### 4.0 The scope decision, and the criterion that produced it

The build plan offers "4 to 6 connectors" and asks for a stated selection criterion. The criterion
used here is **not** proximity in the wiring table. It is:

> **Every wired connector whose reference evaluation is `site_pair_forces(a, b, spring, bound=)`
> with a `UnilateralSpring`** — one reference function, one constitutive class, and three declared
> parameters (`engages_in_tension`, `bidirectional`, `law`) distinguishing the members.

That criterion is checkable rather than editorial, and
`test_every_class_in_the_batch_evaluates_through_the_same_frozen_reference_function` asserts it by
resolving `type(connector).evaluate_sites` to `_DistributedUnilateralConnector.evaluate_sites` for
every member — a **function-object identity** check, `ALEPH-PORT-3501` §3.1's discipline applied to
a class family instead of to a geometry helper.

It selects four classes, and they carry **12 of the 30 wired registry rows**:

| class | module | `engages_in_tension` | `bidirectional` | `law` | rows | kernel family |
|---|---|---|---|---|---|---|
| `CrosslinkConnector` | `connectors_crosslink` | `True` | **`True`** | `LINEAR` | **7** | **A (new)** |
| `BidirectionalLinker` | `connectors_frame` | `True` | **`True`** | `LINEAR` | **2** | **A (new)** |
| `NonAdhesiveContact` | `connectors_solid` | `False` | `False` | **`SOFT_CORE`** | **2** | **B (new)** |
| `TensileLinker` | `connectors_solid` | `True` | `False` | `LINEAR` | **1** | *`-3604`'s, reused* |

> ### ⚠ ROWS AND KERNELS ARE DIFFERENT NUMBERS AND CONFLATING THEM WOULD OVERSTATE THIS PORT
>
> **`wiring.py` counts rows. This entry writes kernels. One class serves many rows.**
> `CrosslinkConnector` alone serves **seven** wired contracts — `dorsal_arc_crosslink`,
> `ecm_crosslink`, `filopodium_cortex_root`, `if_sf_plectin`, `lamellipodium_cortex_seam`,
> `mt_sf_spectraplakin`, `sf_cortex_transient` — because a crosslink is one *element form* whose
> contracts differ only in their parameters and in which owners their two ends land on. So the
> honest report of this lane is **12 registry rows, 2 new kernel families, 8 new kernels (2 true,
> 6 wrong on purpose), and 1 row carried by a kernel another lane already landed.**
>
> `test_the_batch_covers_twelve_wired_registry_rows_with_two_new_kernel_families` asserts both
> numbers against `wiring.WIRED` directly, so neither can drift and neither can be quoted alone.
>
> **The build plan's own count was five, and it is seven.** The plan names five contracts for
> `CrosslinkConnector`; `mt_sf_spectraplakin` and `sf_cortex_transient` are also wired to it and
> were not in that list. Recorded here rather than silently used, because a leverage claim that
> nobody can reproduce from the registry is the same defect class as a tolerance chosen after
> seeing the answer.

**`TensileLinker` is claimed as a covered row and NOT as a new kernel, and the distinction is the
point.** Its spring is `(True, False, LINEAR)`, which is `ErmTether`'s spring exactly;
`ALEPH-PORT-3604` already put that law on `k_tether_step_true`. So `actin_cap_linc` is covered by
running the connector's *own* spring through `-3604`'s case and gate, and
`test_the_tensile_linker_row_is_covered_by_the_kernel_the_previous_lane_landed` is that run. Writing
a second kernel for it would have produced a larger kernel count and no more coverage.

### 4.1 Family A — the two-sided linear spring (9 rows)

```
s0 = |a0 - b0|            s1 = |a1 - b1|                              [um]
x  = s - r                       -- SIGNED, and there is NO max(x, 0)
U(s) = (k/2) x^2                                                      [pN.um]
dU/ds = +k x                                                          [pN]

               { -(U(s1) - U(s0)) / (s1 - s0)     |s1 - s0| >  1e-13   DISCRETE GRADIENT
F_s(s0, s1) =  {
               { -k x_mid                         |s1 - s0| <= 1e-13   MIDPOINT
F_s := 0 wherever `bound` is False                                     RUPTURE

m = ((a0 - b0) + (a1 - b1))/2,   f^a = F_s m/|m|,   f^b = F_s (-m)/|-m|
```

**The one line that is deleted relative to `-3604`'s tether is the whole physics of nine registry
rows.** `UnilateralSpring.engaged_stretch` returns the *signed* offset when `bidirectional` is set,
and the reference says in three separate modules why a reader must not "repair" that: a crosslink
that only pulled would let two crosslinked fibres pass through each other, and `mt_nucleus_linc` is
"the only LINC path that can deliver a compressive nuclear load."

Three consequences follow and each one changes what this gate can and cannot do.

1. **There is no inactive branch, so there is no exact zero to grade.** `-3604`'s strongest oracle
   does not exist for family A. `is_engaged` is true at every separation except `s == r` exactly.
   The two exact anchors that *do* survive are named in §7 (O2), and one of them is new:
   the **rest-length zero**, a site placed at exactly `s = r`, where `U` and `dU/ds` are both
   exactly `0.0` for an algebraic reason and not a small-number reason.
2. **The energy channel is even in `x` and the force channel is odd in it.** `U = (k/2)x²` cannot
   see the sign of the stretch; `dU/ds = kx` is nothing but the sign. This is the sharpest instance
   in the project so far of `ALEPH-PORT-3603`'s "grade the channels separately", and §9's
   `WRONG_ABS_STRETCH` is built directly on it.
3. **`load_carried_in_forbidden_sense_pn` is a structural `0.0` here**, and
   `connectors_crosslink.py` says so itself: the accessor that catches a tether-become-a-strut
   catches nothing on a crosslink, because a crosslink has no forbidden sense. **The audit that
   protects the other family is inert on this one**, and its replacement is an *inverted* control —
   `test_the_kernel_crosslink_carries_load_in_both_senses` asserts the kernel really does push on
   an approached site, so a crosslink silently degraded into a tether fails rather than reporting a
   comfortable zero.

### 4.2 Family B — the compressive soft core along the **line of centres** (2 rows)

```
x  = max(r - s, 0)                                                    [um]
           { (k r / 2) x^2 / s      x > 0      engaged
U(s)    =  { 0        EXACTLY       x = 0      separated
           { -(k r / 2)(r^2 - s^2)/s^2   x > 0
dU/ds   =  { 0        EXACTLY            x = 0
```

then the same discrete-gradient step average, the same rupture mask and the same line-of-centres
direction as family A.

**The scalar law `u(s)` is the one `ALEPH-PORT-3604` ported; the combination is new, and the new
part is where the difficulty is.** `-3604` put `SOFT_CORE` on a kernel through
`normal_contact_forces`, which acts along the cortex vertex normal and **takes no end positions at
all** — it has no step, so it never forms a difference quotient. Family B routes the same divergent
`u(s)` through `scalar_force_pn`, so for the first time in this project a **divergent** energy is
differenced across a step. `U ~ 1/s` near the core, so `U1 - U0` on a moving engaged site is a
cancellation between two large numbers, and it is graded by §10's amplification rather than assumed
benign.

`connectors_solid.py` records why these two contacts are central at all while
`membrane_cortex_contact` is not: neither `membrane_ecm_contact` nor `nucleus_cortex_contact` has a
triangulated side whose normal is the right one to push along at *both* ends. That is a modelling
choice of the frozen module, carried across unchanged, and `normal_contact_forces`'s own docstring —
which shows a central contact is unconditionally unstable tangentially at `-f/s` — is a limit this
entry inherits and repeats in §14 rather than repairs.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| site positions | µm | m | finite |
| separation `s` | µm | m | `> 0` strictly |
| rest length / standoff `r` | µm | m | `> 0` strictly |
| signed stretch `x` (family A) | µm | m | **any sign**, finite |
| clamped stretch `x` (family B) | µm | m | `>= 0`, exactly `0.0` when separated |
| stiffness `k` | pN/µm | N/m | `>= 0` finite |
| energy `U` | pN·µm | J | `>= 0`; exactly `0.0` when separated (B) or at `s = r` (A) |
| force | pN | N | finite; exactly `0.0` when ruptured, when separated (B), or at `s = r` (A) |
| `_DISCRETE_GRADIENT_FLOOR_UM` | µm | m | `1.0e-13`, the reference's own constant, passed in |
| ULP figure | float32 ULP of the field's scale | dimensionless | `>= 0`; `inf` if the reference field is identically zero and the candidate is not |

### Singular cases, and what the device path does

- **A coincident pair** at either end of the step or at the midpoint. `site_pair_forces` raises
  `ConnectorGeometryError` through `_unit_directions`; a Warp kernel cannot raise. **The state
  builders refuse such a configuration on the host in float64 before anything is uploaded**, and the
  kernels add no epsilon to `|m|`. This is `-3604` §5's decision reused; unlike the normal contact,
  neither family here has a *second* refusal predicate that float32 can disagree with the host
  about, because the quantity being tested is a norm and not a signed projection — a norm that is
  positive in float64 is positive in float32 unless it underflows, which at µm scale it does not.
  **So this entry ships one refusal, not two, and that is a difference from `-3604` rather than an
  omission.** `test_a_coincident_central_pair_is_refused_before_any_kernel_launches` drives it.

  > **Corrected during implementation, and the correction is the interesting part.** The refusal was
  > first written inside the state builder only. That left every caller who assembles their own state
  > mapping — which is what the hand-computed anchors of §7 O5 do, and what any downstream user of
  > `crosslink_kernel_site_forces` would do — reaching a kernel with a coincident pair and receiving
  > a **NaN force field that looks like a number**. It now lives in
  > `refuse_coincident_central_pairs`, called by the builder *and* by `_central_pair_launch` before
  > anything is uploaded, and the control drives the launch path rather than the builder.
- **`s -> 0` in family B.** Not refused: it is the divergence the soft core exists to have. The
  state keeps `s` bounded away from zero and declares the bound; §14 says what that does not cover.
- **`s == r` exactly, family A.** Not singular: a legitimate configuration on which the law returns
  exactly `0.0`, and the default state contains such sites deliberately.
- **`bound = False`.** Not singular: exactly zero however stretched or compressed.
- **The ambiguous band** of the discrete gradient — `1e-13 µm < |s1 - s0| <` one float32 ULP of `s`
  — is `-3604` §10's, unchanged, refused by the builder and asserted against.

### Invariants, each with the test that asserts it

- **I1.** Each law compared with itself is `0.00` ULP on every field, exactly —
  `test_a_central_pair_law_compared_with_itself_is_exactly_zero_ulp`.
- **I2.** **Family B's separated branch is exactly `0.0` through the kernel**, asserted with `==`,
  per site, on energy and both force fields —
  `test_the_kernel_central_contact_carries_exactly_zero_on_every_separated_site`.
- **I3.** A ruptured site carries exactly `0.0` through the kernel in **both** families, however
  stretched — `test_a_ruptured_kernel_crosslink_site_carries_exactly_zero_however_stretched`.
- **I4.** **Family A's rest-length zero survives the port**, asserted with `==` on a site placed at
  exactly `s = r` — `test_a_crosslink_site_at_exactly_the_rest_length_carries_exactly_zero`.
- **I5.** **Family A carries load in both senses through the kernel**, with the compressive sites'
  force pointing *apart* and the tensile sites' *together* — the inverted control that replaces the
  inert `load_carried_in_forbidden_sense_pn` — `test_the_kernel_crosslink_carries_load_in_both_senses`.
- **I6.** `F = -grad E` for each kernel's own energy and its own force, by central difference with an
  **observed order 2**, the step read off a measured error curve —
  `test_the_kernel_crosslink_force_is_minus_the_gradient_of_the_kernel_energy`,
  `test_the_kernel_central_contact_force_is_minus_the_gradient_of_the_kernel_energy`.
- **I7.** Every scatter destination is allocated through the backend, never as a host ndarray —
  `test_every_central_pair_scatter_destination_is_allocated_through_the_backend`.
- **I8.** **Both states declare which branches AND which senses they visit, as numbers a test can
  assert on** — `test_the_default_crosslink_state_visits_both_senses_and_both_step_kinds`,
  `test_the_default_central_contact_state_visits_both_the_engaged_and_the_separated_branch`.
- **I9.** The parity report names the position-precision mode alongside the scatter mode —
  `test_the_central_pair_report_names_the_position_precision_and_the_scatter_mode`.
- **I10.** The reference laws are driven on **float64** positions and are measurably not float32
  laws — `test_the_central_pair_references_are_driven_in_float64_positions`.
- **I11.** **The batch's coverage claim is a measurement against `wiring.WIRED`, not a sentence** —
  `test_the_batch_covers_twelve_wired_registry_rows_with_two_new_kernel_families`, and the family
  membership itself is `test_every_class_in_the_batch_evaluates_through_the_same_frozen_reference_function`.

## 6. Source evidence class and known retractions

Nothing is inherited from another project, so there is no inherited evidence class. What is carried,
and it is six findings from four lanes rather than a re-derivation of any of them:

- **`ALEPH-PORT-3601` §14.7** — the accuracy limit is **float32 position storage**, not the scatter
  mode. Honoured: `PositionPrecision` is a declared parameter here too and §13.4 measures all three
  modes at two body offsets.
- **`ALEPH-PORT-3601` §14.2 / `-3604` §14.5** — a case is only as strong as the branches its state
  visits. Honoured: every branch count *and every sense count* is a field of the state.
- **`ALEPH-PORT-3603` §13.4(c) / §14.6 and the build plan's second-precision-channel block** — a
  float32 energy read by a descent criterion gives 0/12 convergence. **No relaxation is put on
  device by this entry**, and no oracle here needs one.
- **`ALEPH-PORT-3603` / `-3604` / `-3605`** — grade energy and each force field **separately**.
  Honoured, and §9 records that 4 of this entry's 6 killable mutants are invisible in some channel.
- **`ALEPH-PORT-3604` §9a and `HANDOFF.md` §F-1** — the `force_b := -force_a` mirror mutant is
  unkillable in a central law with a scalar magnitude. **Re-measured for both new laws before any
  kernel was written** (§9a below): 2,000 draws each, exactly `0` bitwise, including 674 gate
  crossings in family B. **So no mirror kernel is shipped by this entry.**
- **`ALEPH-PORT-3605` §13a M9 and §13.4** — a parity gate cannot see a defect in the *state* it
  grades on, and more precision is not monotone. Both are §14 limits here, and §7's O5 exists
  specifically because of the first.
- **`docs/design/GPU_STATE_2026-07-31.md` and commit `636b0c8`** — `ORDERED` = 0.00 ULP /
  `ATOMIC` = 10.00 ULP on the A5000; `636b0c8` retracts a reported GPU-only backend defect that was
  a raw host ndarray passed as a `scatter_add` destination. This entry has **four** scatter
  destinations and I7 is the control.
- **No retraction was searched for in `/Users/sw1/ffn_cellsim`**, because nothing from that tree is
  used. No magnitude from that project appears in this entry or in the code it authorises.

## 7. Independent oracle or derivation

**A quantity produced by heavy cancellation is not a parity target.** Every claim below is an exact
identity, a ratio, or a declared budget with a round-off justification asserted by a test.

**O1 — the frozen NumPy law**, `aleph.vertical.connectors.site_pair_forces`, in float64, on the
host, with no backend, driven by the **connector's own spring object** obtained from the frozen
builders (`build_ecm_crosslink`, `build_mt_nucleus_linc`, `build_actin_cap_linc`,
`build_nucleus_cortex_contact`). Called, never re-implemented, and obtained from the builder rather
than reconstructed so that a change to a builder's parameters is a test failure and not a drift.

**O2 — the exact zeros that survive, which are not the ones `-3604` had.** Family A has no inactive
branch, so:

| anchor | family | why it is exactly zero |
|---|---|---|
| `bound = False` | A and B | the mask multiplies a finite magnitude by `0.0` |
| `s == r` exactly | **A** | `x = 0` algebraically, so `U = (k/2)·0² = 0` and `dU/ds = k·0 = 0` |
| `s >= r` | **B** | `max(r - s, 0)` is the number `0.0`, and `where(x > 0, core, 0.0)` returns the literal |

All three are asserted with `==`, **per site**, because a per-vertex sum would destroy the property
before anything read it.

> ### ⚠ THE REST-LENGTH ANCHOR IS NOT CONSTRUCTIBLE WHERE THIS ENTRY FIRST PUT IT
>
> This entry planned the `s == r` anchor on `crosslink_state`, whose near-rest sites have an
> **intended** separation of exactly `r`. That does not survive the geometry. The state places
> `a = b + s·n̂` and both sides then recover `s` as `|a − b|`; at a 5 µm anchor radius the round trip
> returns **0.05000019 µm against a 0.05 µm target**, in float32 *and* in float64. A mask of
> `separation == rest` therefore selected **no sites at all**, and every assertion under it passed
> on an empty array — `np.all([]) is True`.
>
> **That is `ALEPH-PORT-3605` §13a M9 in a worse form.** M9 made a parity gate silent by corrupting
> the state it graded; this made an *exact-anchor* control vacuous by constructing a state in which
> the anchor did not exist. A silent gate reports a number nobody can interpret; a vacuous one
> reports a pass.
>
> `crosslink_rest_length_state()` is the geometry in which the anchor survives: **origin-anchored
> and axis-aligned**, `b = 0`, `a = (r, 0, 0)`, so the subtraction is exact and the norm of an
> axis-aligned vector is exactly its one non-zero component. `at_rest_sites` is now computed from
> the **reconstructed** separation everywhere, reports `0` on the general state honestly, and
> `test_the_rest_length_anchor_needs_a_geometry_the_general_state_cannot_provide` is the record.
> §13.6 carries what that geometry then revealed about `LOCAL_F32`.

**O3 — `F = -grad E` by central difference, on each kernel's own energy and its own force, with
the step chosen from a measured error curve.** `PLAN.md` §6.1 records this project losing a night to
steps below the round-off floor. The curve is measured over a decade sweep, the order-2 window is
read off it, and the step is taken from inside that window —
`test_the_central_pair_difference_step_is_chosen_from_a_measured_error_curve`. For family A this
oracle is evaluated on the **stationary** state, where the discrete gradient reduces to `-U'(s)`
exactly; that reduction is a stated property of the rule and this is where it is checked.

**O4 — the sense identity, which is family A's own exact anchor and has no analogue in `-3604`.**
`U(r + δ) = U(r − δ)` **exactly** for any δ, because the energy is even in the stretch, while
`F_s(r + δ) = −F_s(r − δ)` exactly, because the derivative is odd in it. So a compressed site and
its mirror-stretched partner must produce **bit-identical energies and exactly negated forces**
through the kernel. This is an identity of the law at an arbitrary configuration, not an equilibrium
property, and it is the oracle that makes `WRONG_ABS_STRETCH` a defect rather than a taste.

**It is exact in the law and only conditionally exact through the kernel.**
`x = fl32(r ± δ) − fl32(r)` does not give `±` the same magnitude for generic `r` and `δ`: at
`r = 0.05`, `δ = 0.018` the two energies come out **0.008100002655 against 0.008099999302**. The
control therefore asserts the bitwise identity on `r = 2⁻⁴`, `δ = 2⁻⁶`, where both subtractions are
exact, **and** asserts that the generic case holds only to round-off. Recorded rather than quietly
avoided: the law's exact symmetry reaches the kernel only where the arithmetic is exact, which is
the same lesson as the anchor above.

**O5 — the hand-computable values, because a parity number alone says only that two implementations
agree, which is also what two copies of one mistake do.** `-3605` M9 made that concrete: a corrupted
state leaves every ULP figure untouched.

| law | configuration | derived by hand |
|---|---|---|
| A | `k = 50` pN/µm, `r = 0.050` µm, `s = 0.070` µm | `x = +0.020`, `U = ½·50·4e-4 = 0.010` pN·µm, `F_s = −k x = −1.0` pN (pulls together) |
| A | the same at `s = 0.030` µm | `x = −0.020`, `U = 0.010` pN·µm **identically**, `F_s = +1.0` pN (pushes apart) |
| B | `k = 2300` pN/µm, `r = 0.020` µm, `s = 0.010` µm | `U = (kr/2)(r−s)²/s = 23·1e-4/0.01 = 0.23` pN·µm, `−U′ = (kr/2)(r²−s²)/s² = 23·3e-4/1e-4 = 69` pN (pushes apart) |

The two family-A rows are one anchor doing two jobs: they are hand-checkable numbers *and* they are
O4's identity at a concrete pair of configurations.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/runtime/test_central_pair_law_parity.py::test_the_true_crosslink_kernel_agrees_with_the_frozen_numpy_law` | The float32 Warp two-sided-spring kernel, warp CPU device, `ScatterMode.ORDERED`, agrees with `site_pair_forces` on both step energies and both assembled force fields within each field's declared ULP budget — on a state visiting stretched, compressed, at-rest, sense-crossing, moving, stationary and ruptured sites. |
| Positive | `::test_the_true_central_contact_kernel_agrees_with_the_frozen_numpy_law` | The same for the compressive line-of-centres `SOFT_CORE` law, including the divergent energy differenced across a moving engaged step. |
| Positive | `::test_a_central_pair_law_compared_with_itself_is_exactly_zero_ulp` | I1, the harness's zero point on these cases. |
| Positive | `::test_the_kernel_central_contact_carries_exactly_zero_on_every_separated_site` | **I2 / O2**, with `==`, on energy *and* both force fields. |
| Positive | `::test_a_ruptured_kernel_crosslink_site_carries_exactly_zero_however_stretched` | **I3 / O2**, with `==`, on sites far from the rest length in both senses. |
| Positive | `::test_a_crosslink_site_at_exactly_the_rest_length_carries_exactly_zero` | **I4 / O2 — family A's own exact anchor**, the one that replaces the inactive branch it does not have. |
| Positive | `::test_the_kernel_crosslink_carries_load_in_both_senses` | **I5 / §4.1(3)** — the inverted control. Compressed sites push apart, stretched sites pull together, both through the kernel, so a crosslink silently degraded into a tether fails instead of reporting a comfortable zero. |
| Positive | `::test_the_kernel_crosslink_energy_is_even_and_its_force_is_odd_in_the_stretch` | **O4**, the sense identity, bitwise on the energy and exactly negated on the force. |
| Positive | `::test_the_kernel_crosslink_matches_the_hand_computed_mirror_pair` | **O5**, family A's two rows, which are also O4 at a concrete configuration. |
| Positive | `::test_the_kernel_central_contact_matches_the_hand_computed_single_site_value` | **O5**, family B's row. |
| Positive | `::test_the_kernel_crosslink_force_is_minus_the_gradient_of_the_kernel_energy` | I6 / O3 on the stationary state, where the discrete gradient reduces to `-U'(s)`. |
| Positive | `::test_the_kernel_central_contact_force_is_minus_the_gradient_of_the_kernel_energy` | I6 / O3 for the divergent law. |
| Positive | `::test_the_central_pair_difference_step_is_chosen_from_a_measured_error_curve` | The curves themselves — a decade sweep showing the order-2 window *and* the round-off floor, so the step is a measurement rather than a habit. |
| Positive | `::test_an_all_separated_central_contact_passes_at_a_zero_ulp_budget` | The exact zero at **case** level: on an all-separated configuration the whole law is identically zero and the gate passes at a budget of literally `0.0` ULP. |
| Positive | `::test_the_default_crosslink_state_visits_both_senses_and_both_step_kinds` | **I8 for family A** — stretched, compressed, at-rest, sense-crossing, moving, stationary and ruptured counts are all `> 0` and are fields of the state. |
| Positive | `::test_the_default_central_contact_state_visits_both_the_engaged_and_the_separated_branch` | **I8 for family B**. |
| Positive | `::test_every_central_pair_scatter_destination_is_allocated_through_the_backend` | I7, by type, by watching every `scatter_add` call. |
| Positive | `::test_the_central_pair_report_names_the_position_precision_and_the_scatter_mode` | I9. |
| Positive | `::test_the_central_pair_references_are_driven_in_float64_positions` | I10. |
| Positive | `::test_the_position_precision_is_a_parameter_and_not_a_constant_for_the_central_pair_laws` | All three modes reachable through the public builders and producing three different measurements. |
| Positive | `::test_the_declared_central_pair_budgets_sit_inside_what_roundoff_explains` | No declared budget exceeds its own round-off explanation, on every state it grades. |
| Positive | `::test_the_crosslink_state_stays_out_of_the_ambiguous_discrete_gradient_band` | `-3604` §10's band, asserted rather than hoped. |
| Positive | `::test_a_coincident_central_pair_is_refused_before_any_kernel_launches` | §5 — the host float64 refusal, raising the same `ConnectorGeometryError` type the reference raises. |
| Positive | `::test_every_class_in_the_batch_evaluates_through_the_same_frozen_reference_function` | **I11 / §4.0** — the selection criterion, by function-object identity rather than by assertion in prose. |
| Positive | `::test_the_batch_covers_twelve_wired_registry_rows_with_two_new_kernel_families` | **I11 / §4.0** — the coverage claim measured against `wiring.WIRED`: 12 rows, 4 classes, 2 new kernel families. Rows and kernels are separate assertions so neither can be quoted alone. |
| Positive | `::test_the_tensile_linker_row_is_covered_by_the_kernel_the_previous_lane_landed` | §4.0's honest half — `actin_cap_linc`'s own spring, from its own frozen builder, driven through `-3604`'s tether case and gate. One row, zero new kernels, and it says so. |
| Positive | `::test_the_crosslink_kernel_serves_all_seven_wired_crosslink_contracts` | The leverage claim itself: all seven contracts wired to `CrosslinkConnector` resolve to one spring parameterisation and one kernel family. |
| Positive | `::test_the_rest_length_anchor_needs_a_geometry_the_general_state_cannot_provide` | **The correction of §7 O2, kept as a control rather than quietly fixed.** The general state's `at_rest_sites` is `0` and must stay `0`; the anchored state's is not; and the anchor holds under `GLOBAL_F32` and `POSITIONS_F64` and **fails under `LOCAL_F32`**. |
| Positive | `::test_the_degenerate_states_really_are_degenerate` | The blind-spot states are what they claim **before** any bit-identity is asserted on them. A blind-spot control running on a state that was not actually all-stretched would assert an identity that held for the wrong reason. |

## 9. Deliberately failing negative control

**Six** wrong kernels ship in `aleph/runtime/law_kernels.py` beside the two right ones, compiled by
the same `load()` and launched down the same drivers. All are *plausible hand-translation errors*.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/runtime/test_central_pair_law_parity.py::test_a_crosslink_kernel_that_reintroduced_the_unilateral_clamp_is_caught` | `CrosslinkVariant.WRONG_UNILATERAL_CLAMP` puts `max(x, 0)` back, turning nine registry rows of two-sided physics into tethers. **This is the headline mutant of this entry** — it is the exact "repair" that `connectors.py`, `connectors_crosslink.py` and `connectors_frame.py` each warn a reader against in as many words, which is how it gets to be a *likely* error rather than sabotage. Caught on energy and on both forces. |
| Negative (must fail) | `::test_a_crosslink_kernel_that_dropped_the_half_in_its_energy_is_caught` | `CrosslinkVariant.WRONG_ENERGY_HALF` drops the `1/2` from `U`. **Forces bit-identical.** Caught on the energy field alone. |
| Negative (must fail) | `::test_a_crosslink_kernel_that_took_the_magnitude_of_the_stretch_is_caught` | `CrosslinkVariant.WRONG_ABS_STRETCH` uses `|x|` where the derivative needs the signed `x` — the hand translation of a tether kernel in which `x` was non-negative by construction. **Energy bit-identical always** (the energy is even in `x`), and **bit-identical on every *moving* site too**, because the discrete gradient recovers the sign from the energies. It is visible **only on stationary sites** — and `_DistributedUnilateralConnector.accumulate` evaluates every production configuration stationary. See §9b: this is `-3604`'s blind-spot finding with its sign reversed. |
| Negative (must fail) | `::test_a_central_contact_kernel_that_used_the_linear_law_is_caught` | `CentralContactVariant.WRONG_LINEAR_CORE` computes `(k/2)x²` instead of the soft core — the substitution the reference calls a correctness regression rather than a refinement, because a linear contact carrying load `f` presents an interpenetration barrier of only `f²/2k`. Caught on energy and on both forces. |
| Negative (must fail) | `::test_a_central_contact_kernel_that_dropped_the_unilateral_gate_is_caught` | `CentralContactVariant.WRONG_ADHESIVE_GATE` omits the `max(r − s, 0)`, so a separated contact **pulls**. The exact-zero breaker for family B. |
| Negative (must fail) | `::test_a_central_contact_kernel_that_divided_by_the_rest_length_is_caught` | `CentralContactVariant.WRONG_DERIVATIVE_DENOMINATOR` drops the square from the denominator of the **derivative only**, so the reported force stops being the gradient of the reported energy. **Energy bit-identical at every configuration; both force fields bit-identical on every *moving* site.** Caught on the force fields, on stationary sites, and by `F = −grad E`. **Its exact zero is intact and both unilateral-gate controls pass it unchanged**, so nothing else in this entry would have run that could see it. Renamed from the entry's originally-planned `WRONG_CORE_DENOMINATOR` — §9c is why. |
| Negative (must fail) | `::test_a_zero_ulp_budget_fails_the_true_central_pair_kernels` | The vacuity guard on the verdict: at `ulp_budget = 0.0` the *correct* kernels fail, because float32 storage is not bitwise float64. |
| Negative (must fail) | `::test_no_central_pair_mutant_escapes_through_the_atomic_floor` | Neither the widened `ATOMIC` budget nor the scatter floor is an amnesty for any of the six mutants. |

And the blind spots, recorded as controls because they are the honest limits of the method:

| Control | Location | Asserts |
|---|---|---|
| Blind spot (must hold) | `::test_the_reintroduced_clamp_is_invisible_on_an_all_stretched_configuration` | `WRONG_UNILATERAL_CLAMP` is bit-identical to the true kernel when every site is stretched past `r`. G1's all-taut blind spot arriving in the two-sided element, and the reason the default state's compressed-site count is a field a test asserts on. |
| Blind spot (must hold) | `::test_the_magnitude_of_the_stretch_is_invisible_in_the_energy_and_on_every_moving_site` | `WRONG_ABS_STRETCH` is bit-identical on the energy channel at every configuration **and** on both force fields at every moving site. Two independent blindnesses in one mutant, and §9b is why that matters more than it looks. |
| Blind spot (must hold) | `::test_the_adhesive_gate_is_invisible_on_an_all_engaged_central_contact` | `WRONG_ADHESIVE_GATE` is bit-identical when every site is inside the standoff. |
| Blind spot (must hold) | `::test_the_newton_mirror_mutant_is_provably_unkillable_for_the_whole_batch` | §9a: `force_b := −force_a` is bit-identical to the true law for **both** new families, measured on the frozen reference over 2,000 draws each, so no mirror kernel is shipped. |

### 9a. The Newton-pair control this entry **does not** ship, and the measurement that decided it

The build plan requires this lane to check the `HANDOFF.md` §F-1 trap for *its own* laws before
writing a Newton-pair control, and warns specifically that crosslinks are linear springs.

**Measured on the frozen reference before any kernel was written**, 2,000 draws per law, with both
configurations moving:

| law | engaged at start | gate / sense crossings | largest per-site `f^a + f^b` | force scale | `f^b` bitwise `== −f^a`? |
|---|---|---|---|---|---|
| family A, two-sided linear | 2,000 / 2,000 | 0 (there is no gate) | **0.0** | 2.682 pN | **yes** |
| family B, central soft core | 727 / 2,000 | **674** | **0.0** | 168.4 pN | **yes** |
| both, fixed configuration | — | — | **0.0** | — | **yes** |

The reason is structural and not luck, and it is `-3604` §9a's: `site_pair_forces` builds
`f^b = F_s (b − a)/|b − a|`; `|b − a|` is bit-identical to `|a − b|` because the components are
squared, and `(−d)/L` is bit-identical to `−(d/L)` because negation is exact. **It holds for every
member of this batch at every precision, and family B's 674 gate crossings show the divergent law
does not escape it either.**

> **So no mirror kernel is shipped, and no Newton-pair control is shipped.** `-3604` shipped a
> mirror *kernel* for the tether so the claim would be checkable by running it; this entry does not
> ship a second copy, because the lines under test — the two `_store_row` calls that form each
> side's own unit direction — are *the same lines* in these kernels as in `k_tether_step_true`, and
> a kernel that cannot fail is a kernel a future reader has to spend time re-deciding about. What
> ships instead is the measurement, as a control, over both new laws.
>
> **A Newton-pair closure check is a real check only where the two sides can differ**, and nowhere
> in this batch can they. `-3604`'s normal contact was the one place in the connector tree where
> they could.

### 9c. A seventh mutant that was planned, measured, and is not shipped

The entry as written planned `WRONG_CORE_DENOMINATOR`: divide the soft-core **energy** by the rest
length instead of by the separation. It was implemented, and it is **not a mutant**:

> `(k r / 2) x² / r` **is** `(k / 2) x²`.

It is `WRONG_LINEAR_CORE`, algebraically, and it measured as `WRONG_LINEAR_CORE` — the same
4,007,671.63 / 5,354,892.31 / 6,894,574.34 / 7,007,153.57 to two decimal places on every field. A
second kernel that is a rewriting of the first is not a second negative control; it is a control that
cannot distinguish anything, which is the defect this whole track exists to refuse.

What ships in its place is `WRONG_DERIVATIVE_DENOMINATOR`, which touches the derivative and not the
energy and is therefore genuinely independent of every other mutant here. The rejected version is
recorded in the kernel's own comment and asserted as arithmetic inside
`test_a_central_contact_kernel_that_divided_by_the_rest_length_is_caught`, so nobody re-invents it.

### 9b. Two mutants, two opposite blindnesses, and the state the pipeline actually runs

`-3604` §13.1 found that `WRONG_MIDPOINT_ONLY` is invisible on a **stationary** configuration, and
that every configuration `accumulate` evaluates is stationary — a wrong kernel of that shape would
ship undetected today.

`WRONG_ABS_STRETCH` is the mirror image. It is invisible on every **moving** site, because the
discrete gradient reconstructs the force from two energies and the energy is even in the stretch.
It is visible **only** where `-3604`'s mutant is invisible.

> **Neither state alone can grade a two-sided connector.** A case built from the pipeline's own use
> of the law (stationary, `end == start`) is blind to the dropped discrete gradient; a case built
> only from moving steps is blind to the sign of the stretch. The default state here carries both
> kinds and asserts both counts, and that is the requirement rather than a nicety.

## 10. Numerical and precision envelope

**Working precision.** Positions enter in the precision the case declares. Everything derived from
them — the site offsets, the separations, the stretches, the energies used to form the difference
quotient, the midpoint direction and both assembled force fields — is **float32**,
`ALEPH-DQ-107`'s compute channel. Per-site **energies are widened to float64 inside the kernel**
before the product is formed, and the total is reduced by `Backend.sum`'s two-stage float64
accumulation. The references are float64 end to end and I10 asserts it.

**The step-average force is formed in the compute channel**, from the kernel's **float32** energies,
exactly as `-3604` §10 chose and for the same reason: the quotient's output is a force.

**Position precision is a declared parameter, not a buried constant** — `-3603`'s enum, reused
unchanged, for `-3601` §14.7's reason. Raised in
`docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` and **not settled here**.

**The design point that keeps the parameter cheap is preserved.** Exactly **one** kernel in this
entry reads a position: `k_pair_offset_f32` / `_f64`, `ALEPH-PORT-3604`'s kernel, **reused and not
re-written**. Both families are functions of its output. This entry therefore adds **zero** new
places the position-precision decision reaches, which is the first time in the CUDA track that a
lane has added a law without adding one.

**One association-order difference from the reference, inherited and stated again.** The reference
forms the midpoint direction as `(a0 + a1)/2 − (b0 + b1)/2`; the kernels form it as
`((a0 − b0) + (a1 − b1))/2`, the same number in exact arithmetic and better conditioned. `-3604`
§10 records it; it is repeated here because these kernels inherit it rather than because it is new.

### Where the disagreement comes from, derived rather than assumed

Each position is rounded once before anything happens, so with `u = eps32/2` and
`m_i = |a_i| + |b_i|`:

```
both:       |d s_i|   <~ u m_i

family A:   |d U_i|   <~ u k |x_i| m_i                       x SIGNED, so |x| and not the clamp
            |d F_i|   <~ u [ k m_i (1 + 2 |x_i| / s_i) ]                        stationary sites
            |d F_i|   <~ u [ (U0_i + U1_i + |F_i| m_i)/|change_i| + 2 |F_i| m_i / s_i ]   moving

family B:   |d U_i|   <~ u |U'(s_i)| m_i
            |d F_i|   <~ u [ |U''(s_i)| m_i + 2 |F_i| m_i / s_i ]               stationary sites
            |d F_i|   <~ u [ (U0_i + U1_i + |F_i| m_i)/|change_i| + 2 |F_i| m_i / s_i ]   moving
```

with `U''(s) = k r³/s³` on the engaged branch and `0` elsewhere. Family A's energy amplification is
`2 Σ |x_i| m_i / Σ x_i²` — `-3601` §10's expression with the clamp removed, which is the *only*
change the two-sidedness makes to the conditioning. Each amplification is computed **from the host
state by the case**, so it is a property of the data and not a number chosen after seeing the
answer; `roundoff_ulp_bound = (amp + 1)/2` converts it to ULP exactly as `parity.py` already does.

**Family B is the worse-conditioned of the two and the reason is the same as `-3604`'s contact**: a
0.02 µm gate coordinate between surfaces 5 µm from the origin, then divided by `s²` and `s³`. What
is new is that its *step* is now differenced too, so a moving engaged site adds the quotient's own
cancellation on top.

**The atomic floor.** `ATOMIC_SCATTER_ULP_FLOOR = 10.0` (`parity.py`, measured on the A5000). Both
force fields of both families are scatter-assembled and take the floor in `ATOMIC` mode; the
energies are reduced by `Backend.sum`, which is ordered by construction, and do not.

**Outside the envelope.** One device — **warp's CPU device**, on a Mac whose warp build reports
*"CUDA not enabled in this build"*. `ATOMIC` has never been measured here, and
`measure_scatter_determinism` already refuses to let a CPU zero be read as evidence about CUDA.

## 11. Production-backend residency and transfer

| Array | Where | Precision | Transfer |
|---|---|---|---|
| site positions, start and end, both bodies | device | float32 or float64, per `PositionPrecision` | uploaded once per case by `backend.array` |
| `bound` rupture flags | device (`wp.array`) | int32 | uploaded once per case |
| site offsets | device | float32 | never leave the device |
| per-site energy, start and end | device | float64 | never leaves the device before `backend.sum` reduces it |
| assembled per-vertex forces, both sides, both families | device, **allocated by `backend.zeros`** | float32 | read to host once, by the harness |

**Four** scatter destinations across the two cases — two force fields each — which is four chances
to make the mistake `636b0c8` retracted, and I7 is the control that watches every `scatter_add`
call by type rather than by reading the source.

The `bound` array is allocated with `wp.array` and **not** through the backend, for the reason
`-3604` §11 already gives: `WarpBackend` accepts float32 and float64 only, and adding an int32 dtype
there means adding a kernel instantiation and a parity row for it. It is not a scatter destination,
which is what §11's rule is about.

## 12. Comments and docstrings to discard

Nothing was read from another project, so there is no source prose to discard. Kept as a positive
statement of what must **not** appear in the modules this entry authorises, with
`tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package` as the
mechanical half: no provider repository name, module path, kernel name, gate name or branch name; no
absolute path into another project; no provider datum as a literal; and **no claim that a number was
measured on hardware this lane did not run on**. Every ULP figure in the package's docstrings is
either measured here on warp's CPU device and says so, or is cited to
`docs/design/GPU_STATE_2026-07-31.md` / `636b0c8` and says that.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Measured 2026-08-01 on **warp's CPU device** (this machine's warp build reports *"CUDA not enabled in this build"*), float32 compute, `ScatterMode.ORDERED`, `PositionPrecision.GLOBAL_F32` unless stated. `tests/runtime/test_central_pair_law_parity.py`: **41 passed**. `tests/runtime` as a whole: **416 passed, 3 skipped** (was 375 + 3 before this entry). `tests/vertical`: **1,358 passed**, unchanged — the frozen references were read and never edited, and `git status aleph/vertical/` is empty. Mutation study: **11 of 11 killed** (§13a). **12 wired registry rows covered by 2 new kernel families / 8 new kernels**, with a thirteenth row's worth of law (`TensileLinker`) carried by `ALEPH-PORT-3604`'s kernel and claimed as **zero** new kernels. The status stays `PROPOSED` because no CUDA device was touched and no human has reviewed it. |
| Reviewer | **Agent-proposed. Unratified.** Written by a subagent of session `46143f30`; no PI review, no `decided_by` field, and none may be added by an agent. |
| Rollback | Revert the appended blocks of `aleph/runtime/law_kernels.py` and `aleph/runtime/law_cases.py` and delete `tests/runtime/test_central_pair_law_parity.py`. `-3601`, `-3603`, `-3604` and `-3605`'s paths are untouched by construction and keep working. Nothing under `aleph/vertical/` changes in either direction. |

### 13.0 Coverage — the two numbers, and why they are two

| | count | what it is |
|---|---|---|
| **wired registry rows covered** | **12** of 30 | `wiring.WIRED` contracts whose implementing class is in this batch |
| implementing classes | 4 of 13 | `CrosslinkConnector`, `BidirectionalLinker`, `NonAdhesiveContact`, `TensileLinker` |
| **new kernel families** | **2** | family A (two-sided linear), family B (central soft core) |
| **new kernels** | **8** | 2 true, 6 wrong on purpose |
| rows carried by an existing kernel | 1 | `actin_cap_linc`, on `-3604`'s `k_tether_step_true` |
| new position-reading kernels | **0** | `k_pair_offset_*` reused |

`CrosslinkConnector` alone carries **seven** rows. **The build plan's own figure was five** — it
lists five contracts for that class, and `mt_sf_spectraplakin` and `sf_cortex_transient` are also
wired to it. `test_the_batch_covers_twelve_wired_registry_rows_with_two_new_kernel_families` and
`test_the_crosslink_kernel_serves_all_seven_wired_crosslink_contracts` assert both figures against
`wiring.WIRED`, as separate assertions, so neither can be quoted alone.

### 13.1 Family A — the two-sided spring, nine registry rows

Default state: 96 sites over 32+32 vertices — **54 stretched / 33 compressed at the start**, 9
near-rest, 15 sense-crossing, 66 moving / 30 stationary, 14 ruptured.

| kernel | `energy_start` | `energy_end` | `forces_a` | `forces_b` | verdict |
|---|---|---|---|---|---|
| **TRUE** | **3.54** | **12.13** | **65.65** | **49.80** | pass |
| `WRONG_UNILATERAL_CLAMP` | 2,330,084.45 | 2,581,193.53 | 5,859,817.54 | 3,702,429.77 | fail — **all four** |
| `WRONG_ENERGY_HALF` | **8,388,600.92** | **8,388,583.73** | *bit-identical* | *bit-identical* | fail — **energy only** |
| `WRONG_ABS_STRETCH` | *bit-identical* | *bit-identical* | **8,313,435.29** | **6,195,648.23** | fail — **forces only** |
| declared budget | 64 | 64 | 256 | 256 | — |
| round-off explains | 410.67 | 321.04 | 1,074.61 | 1,074.61 | — |

On the degenerate states, which is where the blind spots live:

| state | TRUE | `WRONG_UNILATERAL_CLAMP` | `WRONG_ABS_STRETCH` |
|---|---|---|---|
| all-stretched | 1.89 / 0.16 / 24.53 / 24.30 | ***bit-identical*** | ***bit-identical*** |
| stationary (`end == start`) | 3.54 / 3.54 / 62.73 / 46.40 | 2,330,084 / 5,217,337 / 3,803,314 | **10,434,620 / 7,606,643** |

### 13.2 Family B — the central soft core, two registry rows

Default state: 96 sites, **32 engaged / 64 separated at the start**, 16 gate-crossing, 72 moving / 24
stationary, 14 ruptured, `s ∈ [0.0060, 0.0439]` µm against `r = 0.020` µm, `k = 2300` pN/µm.

| kernel | `energy_start` | `energy_end` | `forces_a` | `forces_b` | verdict |
|---|---|---|---|---|---|
| **TRUE** | **92.19** | **75.73** | **326.66** | **335.49** | pass |
| `WRONG_LINEAR_CORE` | 4,007,671.63 | 5,354,892.31 | 6,894,574.34 | 7,007,153.57 | fail — **all four** |
| `WRONG_ADHESIVE_GATE` | 11,116,817.64 | 5,589,275.86 | 1,364,124.39 | 1,575,338.34 | fail — **all four** |
| `WRONG_DERIVATIVE_DENOMINATOR` | *bit-identical* | *bit-identical* | **3,916,957.93** | **4,022,762.71** | fail — **forces only** |
| declared budget | 512 | 512 | 2,048 | 2,048 | — |
| round-off explains | 1,641.71 | 1,560.44 | 3,858.40 | 3,858.40 | — |

All-engaged state (96 engaged / 0 separated): TRUE 42.37 / 126.78 / 271.56 / 395.15, and
**`WRONG_ADHESIVE_GATE` is bit-identical to TRUE there.**

All-separated state: **every field exactly `0.00`, and the gate passes at a declared budget of `0.0`
ULP.** `WRONG_ADHESIVE_GATE` scores `inf` there (the reference is identically zero and it is not);
**`WRONG_LINEAR_CORE` and `WRONG_DERIVATIVE_DENOMINATOR` both pass at `0.0` too**, because they
still clamp. A budget of zero is the strongest tolerance available and it is still only as strong as
the state it grades.

Stationary state: TRUE 92.19 / 92.19 / **973.10 / 981.02** against round-off bounds
1,641.23 / 1,641.23 / **2,717.67**. **That row set the force budget.** A first pass wrote 4,096,
which no round-off argument on any graded state supports; 2,048 is the number that does, at a 2.1×
margin over the worst true measurement against mutants that clear it by three to four orders.

### 13.3 Six things to read off those two tables

1. **Four of the six killable mutants are invisible in some channel, and in three different ones.**
   `WRONG_ENERGY_HALF` is energy-only; `WRONG_ABS_STRETCH` and `WRONG_DERIVATIVE_DENOMINATOR` are
   force-only *and* stationary-only; `WRONG_UNILATERAL_CLAMP` and `WRONG_ABS_STRETCH` are both
   invisible on an all-stretched state. A gate reporting one aggregated number certifies at least
   three of them.
2. **Two of them are invisible on every moving site, by two unrelated defects.** §9b. The cause is
   the discrete gradient and not either law.
3. **Every killable mutant clears its budget by three to six orders of magnitude**, so the 2.1–5×
   margin between measurement and declared budget is nowhere near the detection threshold.
4. **Every declared budget sits below what round-off explains on every state it grades**, asserted
   by `test_the_declared_central_pair_budgets_sit_inside_what_roundoff_explains` over seven cases.
5. **`ORDERED` vs `ATOMIC` is again not the binding constraint.** Under `ATOMIC` on warp's CPU device
   the crosslink measures 3.54 / 12.13 / 65.65 / **49.73** against 49.80, and the contact is
   unchanged; every mutant still fails. (`measure_scatter_determinism` already refuses to let a CPU
   result be read as evidence about CUDA.)
6. **The two-sided law is better conditioned than the compressive one by roughly 5×**, and the reason
   is structural: family A never divides by the separation, and family B divides by `s²` and `s³`.

### 13.4 Position precision — and a **counterexample to G4's counterexample**

Assembled `forces_a` ULP, true kernels, by mode and by where the body sits:

| mode | crosslink @ origin | crosslink @ (50,0,0) µm | contact @ origin | contact @ (50,0,0) µm |
|---|---|---|---|---|
| `GLOBAL_F32` | 65.65 | **562.06** | 326.66 | **2,948.49** |
| `LOCAL_F32` | 67.41 | **67.41** | 257.54 | **257.54** |
| `POSITIONS_F64` | **1.37** | **1.37** | **1.71** | **1.71** |

> **`POSITIONS_F64` is 48× better than `GLOBAL_F32` on family A and 191× better on family B, at the
> origin, and the gap widens with the offset.** That is the **opposite** of what `ALEPH-PORT-3605`
> §13.4 measured on the enclosed volume, where `POSITIONS_F64` was *worse* than `LOCAL_F32`.
>
> **The two results agree about the mechanism and that is the useful part.** Option D of
> `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` is defined as float64 arithmetic
> *in the subtraction that forms edge vectors and site separations*. Both laws here **are** that
> subtraction, and option D buys one to two orders of magnitude. The enclosed volume has no such
> subtraction — its cancellation is a **reduction** — and option D buys nothing there and costs.
> So the finding for a PI answering §6 is not "float64 helps" or "float64 does not help"; it is
> **option D's benefit is predicted by whether the law's cancellation is a subtraction or a
> reduction**, and both signs of that prediction have now been measured. Raised here, not settled.

`LOCAL_F32` is marginally *worse* than `GLOBAL_F32` for family A at the origin (67.41 against 65.65),
which is `-3604` §13.5's observation repeated: the body is already centred, so the centroid
subtraction is a small translation with no cancellation to remove.

### 13.5 `F = −grad E` — the step read off a curve, and off the **right probe**

Both curves measured before either control was written, `POSITIONS_F64`, on the stationary states.

**Family A, on its COMPRESSED sites** (`s ≈ 0.020` µm against `r = 0.050`):

| `h` [µm] | 3.0e-3 | 1.0e-3 | 3.0e-4 | 1.0e-4 | 3.0e-5 | 1.0e-5 |
|---|---|---|---|---|---|---|
| relative | 1.806e-02 | 2.020e-03 | **2.226e-04** | 1.290e-04 | 1.627e-04 | 1.421e-03 |
| ratio | — | **8.94** | **9.08** | 1.73 | 0.79 | 0.11 |

Order 2 for a 3× step reduction is 9.0; measured **8.94 and 9.08**, with the round-off floor visible
from `h = 1e-4` down. Step used: **3e-4**.

> **The first draft of this control read its step off the wrong probe and was 36× outside its own
> window.** The curve was measured on the four most heavily loaded sites, which at `s ≈ 0.088` µm are
> all *stretched* — and `WRONG_ABS_STRETCH` is bit-identical to the true kernel on every stretched
> site, so the mutant half of the control was also passing for the wrong reason. The truncation term
> of a central difference through `s = |a − b|` carries the geometric curvature `~(h/s)²`, so the
> same 3e-3 step that gives 8.5e-4 on the stretched sites gives **1.8e-2** on the compressed ones.
> **A step read off a curve measured on a different probe than the control uses is not a measured
> step.** Both halves now run on the compressed probe.

`WRONG_ABS_STRETCH` measures a relative error of **exactly 2.0000** against the same numeric
derivative — `|−F − F| / |F|`. A flipped sign has a signature, and the control asserts the signature
rather than "large".

**Family B:** ratios **10.73** at `h = 3e-4` and **9.10** at `h = 1e-4` (order 2 is 10.7 and 9.0 for
those reductions), floor at `h = 1e-5`. Step used: **1e-4**, relative error 2.58e-4.

### 13.6 The exact anchors, and the one that moved

| oracle | measured through the kernel | reference / expected |
|---|---|---|
| family A, ruptured | **exactly `0.0`** on all 14 sites, both force fields, and both senses present among them | exact |
| family A, `s == r` exactly | **exactly `0.0`** on all 8 anchored sites, energy and both forces, under `GLOBAL_F32` **and** `POSITIONS_F64` | exact |
| family A, `s == r` under `LOCAL_F32` | **NOT exact — 1.863e-07 pN** of spurious force | see below |
| family B, separated at both ends | **exactly `0.0`** on every such site, energy and both forces | exact |
| family B, all-separated state | **every field exactly `0.0`**; gate passes at budget `0.0` | exact |
| family A, sense identity | `U(r+δ)` **bitwise** `U(r−δ)` and `F` exactly negated, at `r = 2⁻⁴`, `δ = 2⁻⁶` | exact |
| per-site `f^a + f^b`, both laws | **bitwise zero, 2,000 draws each**, incl. 674 gate crossings in family B | §9a: unkillable |
| hand-computed family A, `s = 0.070` / `0.030`, `r = 0.050`, `k = 50` | `U = 0.010` pN·µm **at both**, `F_s = ∓1.0` pN | derived in §7 O5 |
| hand-computed family B, `s = 0.010`, `r = 0.020`, `k = 2300` | `U = 0.23` pN·µm, `−U′ = 69.0` pN | derived in §7 O5 |
| `F = −grad E`, family A, compressed | ratios **8.94 / 9.08** | order 2 |
| `F = −grad E`, family B | ratios **10.73 / 9.10** | order 2 |

> **The `LOCAL_F32` row is a finding and not a defect.** `LOCAL_F32` subtracts a centroid in float64
> and uploads offsets, which destroys the axis alignment the exactness depends on. So a statement
> the *law* makes exactly — a crosslink at its rest length carries no load — **survives the port
> under two position representations and not under the third**.
>
> `ALEPH-PORT-3604` §13.5 found that the position representation moves **which configurations the
> model accepts**. This is the same decision moving **which exact statements the model can make**.
> A third consequence hanging on §6 of the open proposal, raised here and not settled.

### 13a. Mutation study — 11 planted defects, 11 killed

Run with `PYTHONDONTWRITEBYTECODE=1` **and with every `__pycache__` under `aleph/` and `tests/`
deleted first**, then re-run from a verified-empty bytecode tree — a same-length mutant with a
sub-second edit/revert leaves a `.pyc` CPython considers valid, and it has already contaminated one
study in this repository. Baseline green, each mutant applied alone, reverted before the next, green
again after the last.

| # | Planted defect | Killed by |
|---|---|---|
| M1 | the driver ignores the requested variant and always launches the true kernel | 11 tests |
| M2 | a scatter destination allocated as a host ndarray instead of `backend.zeros` | 17 tests |
| M3 | **the coincidence refusal is dropped from the *launch* path** (left in the state builder) | 1 test |
| M4 | the rupture mask is uploaded as all-bound | 3 tests |
| M5 | **`at_rest_sites` reported from the *placed* separation** — the vacuous-pass bug, replanted | 2 tests |
| M6 | the exact-anchor geometry loses its axis alignment, so `s == r` stops being exact | 2 tests |
| M7 | the frozen reference is handed float32-rounded positions | 2 tests |
| M8 | **the *reference* spring loses its two-sidedness**, so both sides are wrong together | 7 tests |
| M9 | the **true** family-A kernel's midpoint force scaled by 1.0001 | 5 tests |
| M10 | the **true** family-B kernel's unilateral gate leaks (`> −1e-4` instead of `> 0`) | 2 tests |
| M11 | the census claims `TensileLinker` needed a new kernel — the leverage claim inflated | 1 test |

**M8 is this entry's own contribution to the list, and it is `ALEPH-PORT-3605` M9's lesson aimed at
the half a parity gate really cannot see.** M9 corrupted the *state*; M8 corrupts the **reference**,
by giving the frozen `UnilateralSpring` `bidirectional=False`. Both sides then compute a tether, the
kernel and the reference agree beautifully, and nine registry rows of physics have quietly gone. It
is killed only because three controls ask questions that are not "do the two sides agree": the
inverted both-senses control, the mutants' expected *channels*, and the round-off budgets. **A pure
agreement gate is silent on it**, which is the sharpest available statement of why §7's exact anchors
and §8's inverted control are not decoration.

**M5 and M6 are the two halves of the §7 O2 correction, replanted as mutants** so the repair cannot
regress silently — M5 restores the accounting that made the anchor vacuous, M6 destroys the geometry
that makes it real.

**M11 is unusual and deliberate: it is a mutant of a *claim* rather than of a computation.** This
lane's headline number is a coverage figure, and a coverage figure is exactly the kind of thing that
drifts into prose and stops being checked. It is killed by the assertion against `wiring.WIRED`.

## 14. Honest limits

What this entry does **not** establish:

1. **No CUDA device was touched.** Every number is from **warp's CPU device**, where a launch is a
   serial loop. It is evidence that the kernels' algebra matches the reference and that the
   identities hold; it is `UNVERIFIED` as evidence about CUDA.
2. **12 rows is not 12 connectors' worth of physics.** The rows share a class because they share an
   element *form*; what differs between `ecm_crosslink` and `mt_sf_spectraplakin` is their
   parameters, their owners and their kinetics, and **none of the three is ported here**. A passing
   gate on this batch is evidence about one force law evaluated at one parameterisation, not about
   seven contracts.
3. **No kinetics.** `bound` is a flag nothing in this vertical flips, and the kernels read it rather
   than evolving it. Nothing here is evidence about crosslink lifetime, LINC lifetime, cortical
   capture duration, plectin or spectraplakin turnover, or force-dependent rupture.
4. **`if_nucleus_linc` is deliberately excluded and it is the interesting exclusion.**
   `StrainStiffeningCableLink` is the one class in `connectors_frame` this batch does *not* cover,
   because its spring is `StrainStiffeningCable` — a duck-typed stand-in that forwards to
   `intermediate_filament.tension_pn`, a three-branch card law with a knee and a plateau. It is a
   different force-law family by the §4.0 criterion and belongs to a lane that can grade its branch
   structure properly. Excluding it is the criterion working, not the criterion being convenient.
5. **The negative controls are as strong as the state they run on, and three of the six are provably
   blind somewhere** — §9's blind-spot table names each one and the configuration that hides it.
   **`WRONG_ABS_STRETCH` is blind on every moving site and in the energy channel at every
   configuration**, which leaves exactly one channel on one kind of site able to see it.
6. **A Newton-pair control is impossible for this whole batch, not omitted** (§9a). That is a limit
   of the law's own structure and it holds at every precision.
7. **A parity gate cannot see a defect in the state it grades on, and this lane met the sharper
   form of it** (`-3605` §13a M9). Both state builders here construct their own geometry, so a
   defect in a builder is invisible to the gate and is caught only by O5's hand-computed values and
   by the branch- and sense-count assertions. **That is why §7 carries three hand-computed rows
   rather than none**, and it is still not a mechanical guarantee.

   Worse than silence is available, and this entry produced it: an exact-anchor control written
   against a state in which the anchor did not exist **passed on an empty mask** (§7 O2). A silent
   gate reports a number nobody can interpret; a vacuous control reports a pass. M5 and M6 replant
   both halves so the repair cannot regress, but nothing mechanically prevents the next empty mask.
8. **A pure agreement gate cannot see a defect in the *reference* either** (§13a M8). Giving the
   frozen spring `bidirectional=False` makes both sides compute a tether, in perfect agreement, with
   nine registry rows of physics gone. Only the non-agreement controls catch it. **This is the limit
   that most directly threatens a lane whose value proposition is leverage**, because one wrong
   parameter would silently retract all nine rows at once.
9. **`LOCAL_F32` does not preserve family A's exact rest-length anchor** (§13.6). That is reported,
   not repaired: which representation the port carries is a PI decision, and this entry adds a
   consequence to it rather than choosing.
10. **The declared budgets are engineering numbers and they belong to one state each.** The
   cancellation grows with the distance of the body from the origin and with site density, so a case
   on another state must declare its own budget.
11. **The energy is compared as a scalar total.** A kernel wrong on two sites in opposite directions
   would cancel and pass — and family A makes this *sharper* than in any previous entry, because its
   energy is even in the stretch, so a kernel that mirrored a site's stretch about `r` would be
   invisible in the energy total by construction. The per-site exact-zero and sense checks are not
   subject to this; the totals are. Inherited from `-3601` §14.6 and named again.
12. **The ambiguous band of the discrete gradient is not ported** (`-3604` §14.2, unchanged). The
    states stay out of it and a control asserts that.
13. **Family B's tangential instability is inherited, not repaired.** `normal_contact_forces`'s
    docstring shows a central contact transmitting an outward load `f` across a separation `s` has a
    transverse stiffness of `−f/s` and is unconditionally unstable. `membrane_ecm_contact` and
    `nucleus_cortex_contact` are central by the frozen module's choice; putting that law on a kernel
    makes it faster and no more stable, and **nothing here is evidence that either contact is
    usable in a descent**.
14. **No relaxation is on device and none is proposed**, per `-3603` §14.6 and the build plan's
    second-precision-channel block. No oracle here needs one: every identity used holds at an
    arbitrary configuration.
15. **Nothing here measures throughput.** The cost of these kernels against the NumPy path, and the
    FP64 cost of `POSITIONS_F64`, are unmeasured on any device.
16. **Evidence class.** `STRUCTURAL` for the parity figures; `ANALYTIC_ORACLE` for the exact zeros,
    the sense identity and `F = −grad E`. No energy or force in this entry may be reported as a
    property of a cell.
