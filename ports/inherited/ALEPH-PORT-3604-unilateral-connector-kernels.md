# ALEPH-PORT-3604 — the two unilateral connectors as Warp kernels

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3604` |
| Lane | `46143f30` Lane G3 (CUDA track), Track G task **G3** |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Depends on | `ALEPH-PORT-3601` (the law-level parity harness) and `ALEPH-PORT-3603` (the membrane kernels, whose `PositionPrecision` parameter and `k_face_edges_*` stage this entry **reuses rather than duplicates**). |
| Exists because | G1 ported a smooth per-site expression and G2 a smooth per-face and per-hinge one. The vertical's load path is carried by two **unilateral** elements, and their defining property is not smoothness but a *bit-level* one: the inactive branch returns **exactly `0.0`**, the number and not a small number. `PLAN.md` §6.1 is explicit that this is testable with `==` because `dU/ds` is identically zero on the whole branch rather than merely small near the gate. That is the property a kernel can break silently, and no gate in this project has ever measured it on a device path. |

---

## 1. Aleph API

```python
from aleph.runtime.law_kernels import (
    TetherStepVariant,          # TRUE | WRONG_MIDPOINT_ONLY | WRONG_IGNORES_RUPTURE
                                #      | WRONG_NEWTON_MIRROR   (a shipped BLIND SPOT -- see 9a)
    ContactVariant,             # TRUE | WRONG_LINEAR_CORE | WRONG_ENERGY_HALF
                                #      | WRONG_ADHESIVE_GATE | WRONG_FROZEN_NORMAL
    CONTACT_TRUE_STAGES,        # the default kernel of each contact stage, in launch order
    CONTACT_STAGE_KERNELS,      # which stage each ContactVariant swaps
    PositionPrecision,          # unchanged, from ALEPH-PORT-3603
    load, loaded_warp,          # unchanged, extended to compile the connector kernels
)

from aleph.runtime.law_cases import (
    CONTACT_STANDOFF_UM,             # 0.020 um, the vertical's own default
    CONTACT_STIFFNESS_PN_PER_UM,     # 2300 pN/um, the scale build_vertical derives
    TETHER_STEP_REST_LENGTH_UM,      # 0.050 um, the vertical's own default
    TETHER_STEP_STIFFNESS_PN_PER_UM, # 50 pN/um, likewise
    DEFAULT_CONTACT_ULP_BUDGET,
    DEFAULT_TETHER_STEP_ULP_BUDGET,
    contact_law_case,           # the CortexMembraneContact law as a LawCase
    contact_state,              # reproducible cortex surface + paired membrane sites
    contact_kernel_forces,      # the kernels' own energy and both force fields, for the oracles
    tether_step_law_case,       # the full ErmTether law -- rupture mask and discrete gradient
    tether_step_state,          # reproducible two-configuration tether geometry
)
```

Nothing outside this list is covered. In particular this entry does **not** authorise any change to
`aleph/vertical/**` — **`aleph/vertical/connectors.py` and `aleph/vertical/membrane.py` are the
frozen parity reference and are read, never edited** — nor to `aleph/runtime/parity.py`,
`aleph/runtime/backend.py`, `aleph/runtime/warp_kernels.py`, `aleph/scenarios/**`, `aleph/viz/**`,
or `ports/ledger/INDEX.md`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read.** No file under that tree was opened by this lane. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

The laws being ported are Aleph's own: `aleph/vertical/connectors.py`, whose provenance is
`ALEPH-PORT-1103-tensile-tether-and-compressive-contact.md`, and the two mesh helpers
`_face_terms` / `_vertex_terms` it imports from `aleph/vertical/membrane.py`
(`ALEPH-PORT-1101`). This entry moves already-ported laws onto kernels; it re-ports nothing.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported from another project.** As in `ALEPH-PORT-3603`: the
discretisation, the association order of every product and the exact placement of every clamp are
fixed by the frozen module, because the whole purpose is to agree with it. A clean-room second
implementation of a unilateral law would put its `max(·, 0)` in a different place and would disagree
on the gate, which is the one thing a parity gate for *this* law exists to measure.

## 4. Physical or mathematical law represented

Two laws, deliberately kept as **two separate `LawCase`s** and **accumulated independently**. The
tether is tensile-only and the contact compressive-only; fusing them would make the force-closure
check a tautology and would hide which of the two a defect lives in.

### (a) `ErmTether` — the full connector law, not only the fixed configuration

G1 ported the fixed-configuration special case (`end == start`, no rupture mask). The law
`site_pair_forces` actually implements has **two more branch families**, and this entry ports them:

```
s0 = |a0 - b0|                  s1 = |a1 - b1|                          [um]
x  = max(s - r, 0)              U(s) = (k/2) x^2                        [pN.um]

                { -(U(s1) - U(s0)) / (s1 - s0)      |s1 - s0| >  1e-13    DISCRETE GRADIENT
F_s(s0, s1) =   {
                { -U'((s0 + s1)/2) = -k x_mid       |s1 - s0| <= 1e-13    MIDPOINT
F_s := 0 wherever `bound` is False                                        RUPTURE

m   = (a0 + a1)/2 - (b0 + b1)/2                                          [um]
f^a = F_s m/|m|                 f^b = F_s (-m)/|-m|
```

`_DISCRETE_GRADIENT_FLOOR_UM = 1.0e-13` is a module constant of the reference and is passed to the
kernel rather than re-chosen. The discrete gradient is the rule that **closes the energy books
across a gate crossing**, which a midpoint force does not, and it is the reason
`WRONG_MIDPOINT_ONLY` (§9) is a plausible error rather than sabotage.

### (b) `CortexMembraneContact` — `ContactLaw.SOFT_CORE` along the cortex vertex normal

This is the connector's **default** law (`CortexMembraneContact.build` selects `SOFT_CORE` unless
the deliberately-adhesive variant is requested), and it is a different shape from anything ported so
far. It is not central, its gate coordinate is a **signed** normal standoff rather than a distance,
and its energy **diverges** as the standoff closes.

```
N_f  = (p1 - p0) x (p2 - p0)                                             [um^2]
m_i  = (1/2) sum_{f > i} N_f          n_i = m_i / |m_i|                  area-weighted vertex normal
d_i  = x^m_i - x^c_i                  s_i = d_i . n_i                    SIGNED standoff [um]

           { (k r / 2) (r - s)^2 / s        s < r        engaged (compressive)
U(s_i) =   {
           { 0            EXACTLY           s >= r       separated
           { -(k r / 2)(r^2 - s^2)/s^2      s < r
U'(s_i) =  {
           { 0            EXACTLY           s >= r

f^m_i = -U'(s_i) n_i
f^c_i = +U'(s_i) n_i  -  g_i          with  g = grad_{x^c} [ sum_i U(s_i) ]  through n_i(x^c):

  t_i  = d_i - s_i n_i                            tangential misregistration [um]
  c_i  = U'(s_i) t_i / |m_i|                                              [pN/um]
  C_f  = sum_{k} c_{tri(f,k)}
  g_{tri(f,k)} += (1/2) C_f x (p_{k+2} - p_{k+1})
```

**Three properties of this law are load-bearing and are ported verbatim rather than improved.**

1. **The standoff is signed and the law refuses a non-positive one.** `normal_contact_forces` raises
   `ConnectorGeometryError`. §5 states exactly what the device path does instead, because a kernel
   cannot raise.
2. **The soft core is not a refinement of the linear contact, it is a correctness fix**, and the
   reference says why: a linear unilateral contact carrying a load `f` presents an interpenetration
   barrier of only `f^2/2k`, which a descent finds and punches through one site at a time. A kernel
   that quietly used the linear law would produce a model whose global minimum has the membrane
   inside the cortex. That is `ContactVariant.WRONG_LINEAR_CORE`.
3. **The force is *not* a Newton pair site by site**, and `central_force` is reported `False`
   honestly. The closures that do hold are the **resultant** ones, by translation and rotation
   invariance of the potential. §7 and §9a state what that can and cannot catch.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| site / vertex positions | µm | m | finite |
| separation `s`, standoff `s` | µm | m | separation `> 0`; **standoff `> 0` strictly** |
| rest length `r`, standoff `r` | µm | m | `> 0` strictly |
| stiffness `k` | pN/µm | N/m | `>= 0` finite |
| energy `U` | pN·µm | J | `>= 0`, and **exactly `0.0`** on the inactive branch |
| force | pN | N | finite; **exactly `0.0`** on the inactive branch |
| `_DISCRETE_GRADIENT_FLOOR_UM` | µm | m | `1.0e-13`, the reference's own constant |
| ULP figure | float32 ULP of the field's scale | dimensionless | `>= 0`; `inf` if the reference field is identically zero and the candidate is not |

### The singular case a kernel cannot express, and the decision this lane took

**A non-positive normal standoff.** `normal_contact_forces` raises `ConnectorGeometryError`; a Warp
kernel has no exception. The build plan requires this lane to decide and state what the device path
does, and to ship a control that drives whichever it chose.

> **Decision: the device path refuses, twice, and never clamps.**
>
> 1. **Before launch.** `contact_state()` and the candidate driver compute the signed standoff on the
>    host in float64 and raise **`ConnectorGeometryError` — the same exception type the frozen
>    reference raises** — if any site is non-positive. So the ordinary bad configuration never
>    reaches a kernel, and the two sides agree on the *refusal* as well as on the numbers.
> 2. **After launch.** The site kernel writes `violation[i] = 1` whenever the standoff **it** formed,
>    in the position precision the case declared, is non-positive; for that site it writes exactly
>    `0.0` energy and `0.0` force rather than dividing by it. The driver reads the flag array and
>    raises. This second refusal is not belt-and-braces: the host predicate is evaluated in float64
>    and the kernel's is not, so **a standoff that is positive in float64 can be non-positive in
>    float32**, and whether it is depends on `PositionPrecision`. §13.5 measures that.
>
> **Not chosen, and why.** A *sentinel* return (NaN, or `inf` energy) can be ignored by accident and
> would travel into an energy ledger as a number; an exception cannot. *Clamping* the standoff at a
> small positive epsilon replaces a geometric failure with a large finite force pointing in an
> arbitrary direction, which is precisely what `ConnectorGeometryError`'s own docstring refuses, and
> it would make the divergent core — the whole reason `SOFT_CORE` exists — finite again.

Other singular and boundary cases, each with the behaviour Aleph requires:

- **Coincident tether sites** (`s = 0` at either end of the step, or at the midpoint). The reference
  raises `ConnectorGeometryError`. The state builder rejects such a state before either side runs;
  the kernel adds no epsilon to `|m|`.
- **A degenerate cortex triangle**, or a vertex whose incident face normals cancel (`|m_i| = 0`).
  The reference raises `DegenerateGeometryError`. Refused by the builder.
- **`bound = False`.** Not singular: a legitimate branch that must return **exactly zero**, and the
  default state contains ruptured sites so `WRONG_IGNORES_RUPTURE` is detectable.
- **A step whose `|s1 - s0|` lands between `1e-13` µm and one float32 ULP of `s`.** *Not* refused,
  but declared out of scope and asserted against — see §10, "the ambiguous band".

Invariants that must hold, each with the test that asserts it:

- **I1.** Each law compared with itself is `0.00` ULP on every field, exactly —
  `test_a_connector_law_compared_with_itself_is_exactly_zero_ulp`.
- **I2.** **The inactive branch is exactly `0.0` through the kernel**, asserted with `==` and not a
  tolerance, on both laws and on every inactive site —
  `test_the_kernel_tether_carries_exactly_zero_on_every_slack_site`,
  `test_the_kernel_contact_carries_exactly_zero_on_every_separated_site`.
- **I3.** A ruptured tether site carries exactly `0.0` through the kernel whatever its separation —
  `test_a_ruptured_kernel_tether_site_carries_exactly_zero_however_taut`.
- **I4.** `F = -grad E` for the kernel's own energy and its own force, by central difference with an
  **observed order 2**, and the step chosen from a measured error curve rather than a formula —
  `test_the_kernel_contact_force_is_minus_the_gradient_of_the_kernel_energy`.
- **I5.** The contact's resultant closure `sum f^m + sum f^c = 0` holds through the kernel, judged
  against `sum |f|` and not against the vanishing resultant —
  `test_the_kernel_contact_resultant_closure_holds_against_the_constituent_scale`.
- **I6.** Every scatter destination is allocated through the backend, never as a host ndarray —
  `test_every_connector_scatter_destination_is_allocated_through_the_backend`.
- **I7.** **Both states declare which branches they visit, as numbers a test can assert on**, and
  both branches of each unilateral element are visited —
  `test_the_default_tether_step_state_visits_both_branches_and_both_step_kinds`,
  `test_the_default_contact_state_visits_both_the_engaged_and_the_separated_branch`.
- **I8.** The parity report names the position-precision mode alongside the scatter mode —
  inherited from `ALEPH-PORT-3603` I8 and re-asserted on these cases by
  `test_the_connector_report_names_the_position_precision_and_the_scatter_mode`.
- **I9.** The reference laws are driven on **float64** positions and are measurably not float32 laws
  — `test_the_connector_references_are_driven_in_float64_positions`.
- **I10.** A non-positive standoff is **refused, not clamped**, on the host path *and* on the device
  path — `test_a_non_positive_standoff_is_refused_before_any_kernel_launches`,
  `test_a_standoff_the_host_calls_positive_and_the_kernel_calls_zero_is_refused_after_launch`.

## 6. Source evidence class and known retractions

Nothing is inherited from another project, so there is no inherited evidence class. What is carried:

- **`ALEPH-PORT-1103`** — the reference laws' own ledger entry, in this repository, and
  `tests/vertical/test_connectors.py`'s exact-zero controls
  (`test_compressed_tether_carries_exactly_zero`, and its partner
  `test_bidirectional_tether_is_caught_carrying_compression` which proves that check has teeth).
  Those are the reason these laws can be treated as frozen.
- **`ALEPH-PORT-3601` §14.7 / §14.2** — the float32-position finding, and G1's stated limit that the
  missing-clamp mutant is provably invisible on an all-taut state. Both are honoured here rather
  than rediscovered: the position precision is a declared parameter, and **every case declares its
  branch counts**.
- **`ALEPH-PORT-3603` §13.4(c) / §14.6** — the second precision channel. **No relaxation is put on
  device by this entry**, and no oracle here needs one.
- **`docs/design/GPU_STATE_2026-07-31.md` and commit `636b0c8`** — `ORDERED` = 0.00 ULP /
  `ATOMIC` = 10.00 ULP on the A5000. `636b0c8` retracts `278e6e5`, which reported a GPU-only backend
  defect that was a raw host ndarray passed as a `scatter_add` destination. This entry has **three**
  scatter destinations and I6 is the control.
- **`HANDOFF.md` §F-1** — the measured statement that in an exactly linear law a
  `force_b := -force_a` mutant is unkillable, because `-(zeta v)` and `zeta (-v)` are bit-identical
  in IEEE-754. §9a re-measures it for *this* law before shipping any Newton-pair control.
- **No retraction was searched for in `/Users/sw1/ffn_cellsim`**, because nothing from that tree is
  used. No magnitude from that project appears in this entry or in the code it authorises.

## 7. Independent oracle or derivation

**A quantity produced by heavy cancellation is not a parity target** (Lane A's nucleus net force is a
15-order cancellation where one ULP moves the answer 40%). For these two connectors the *exact*
anchors are unusually good, so they are used in preference to agreement.

**O1 — the frozen NumPy laws**, `aleph.vertical.connectors.site_pair_forces` and
`normal_contact_forces`, in float64, on the host, with no backend. Called, never re-implemented.

**O2 — the exact zero, which is the property this task exists for.** `dU/ds` is identically zero on
the whole inactive branch, so a correct kernel returns `0.0` bitwise on every slack tether site,
every separated contact site and every ruptured site. This is asserted with `==`. It is the strongest
oracle here because it admits no tolerance to widen: a kernel that computes the active expression and
multiplies by a 0/1 mask passes it (`0.0 * finite == 0.0`); a kernel that evaluates `max(x, 0)` on a
value that is `-1e-8` rather than `0` does not.

**O3 — `F = -grad E` by central difference, on the *kernel's own* energy and the *kernel's own*
force, with the step chosen from a measured error curve.** `PLAN.md` §6.1 records this project
losing a night to steps below the round-off floor, and Lane A found the clean band three orders of
magnitude from the generic estimate. So the curve is measured over a decade sweep, the order-2
window is read off it, and the step is taken from inside that window. **This is the only oracle that
sees the contact's normal-gradient term**, because that term is invisible in the energy and in the
membrane-side force.

**O4 — the resultant closures, which hold by symmetry rather than by construction.** The contact's
potential is translation- and rotation-invariant, so
`sum_i f^m_i + sum_i f^c_i = 0` and the total moment vanishes — exactly, at any configuration, with
no equilibrium anywhere. Measured on the reference in float64: **6.36e-14 pN against a constituent
scale of 1,626 pN**. Judged against `sum_i |f_i|` and never against the resultant itself, for the
reason `connectors.py` gives on `constituent_force_pn`: on a concentric pair the resultant cancels to
round-off by symmetry, and a check judged against it gets stricter the more correct the connector is.

**O5 — the hand-computable value.** One contact site with the membrane exactly along the normal at
`s = 0.010 µm`, `r = 0.020 µm`, `k = 2300 pN/µm`: `U = (k r/2)(r - s)^2/s = 23 * 1e-4 / 0.01`
`= 0.23 pN·µm` and `-U' = (k r/2)(r^2 - s^2)/s^2 = 23 * 3e-4 / 1e-4 = 69 pN`, pushing outward. A
parity number alone says only that two implementations agree, which is also what two copies of one
mistake do.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/runtime/test_connector_law_parity.py::test_the_true_tether_step_kernel_agrees_with_the_frozen_numpy_law` | The float32 Warp tether kernel, warp CPU device, `ScatterMode.ORDERED`, agrees with `site_pair_forces` on both step energies and both assembled force fields, within each field's declared ULP budget — on a state that visits taut, slack, gate-crossing, stationary and ruptured branches. |
| Positive | `::test_the_true_contact_kernel_agrees_with_the_frozen_numpy_law` | The same for the `SOFT_CORE` normal-standoff contact, including the normal-gradient term of the cortex force. |
| Positive | `::test_a_connector_law_compared_with_itself_is_exactly_zero_ulp` | I1, the harness's zero point on these cases. |
| Positive | `::test_the_kernel_tether_carries_exactly_zero_on_every_slack_site` | **I2 / O2 for the tether**, with `==`. |
| Positive | `::test_the_kernel_contact_carries_exactly_zero_on_every_separated_site` | **I2 / O2 for the contact**, with `==`, on energy *and* on both force fields. |
| Positive | `::test_a_ruptured_kernel_tether_site_carries_exactly_zero_however_taut` | I3, with `==`, on a site whose separation is far above the rest length. |
| Positive | `::test_the_kernel_contact_force_is_minus_the_gradient_of_the_kernel_energy` | I4 / O3 — both force fields, order 2 observed, step read off the measured curve. |
| Positive | `::test_the_kernel_tether_force_is_minus_the_gradient_of_the_kernel_energy` | I4 for the tether, on the **stationary** state, where the discrete gradient reduces to `-U'(s)` exactly. That reduction is a stated property of the rule and this is where it is checked. |
| Positive | `::test_the_central_difference_step_is_chosen_from_a_measured_error_curve` | The curve itself: a decade sweep with the order-2 window and the round-off floor both visible, so the step is a measurement rather than a habit. |
| Positive | `::test_an_all_separated_contact_passes_at_a_zero_ulp_budget` | **The exact zero at case level, and the strongest form the claim can take.** On an all-separated configuration the whole law is identically zero, so the parity gate passes at a budget of literally `0.0` ULP. No other law in this project can be graded that way. |
| Positive | `::test_the_position_precision_changes_whether_the_geometric_refusal_fires` | The finding of §13.5: at a `1e-7` µm standoff the *same* configuration is refused under `GLOBAL_F32` and runs clean under `POSITIONS_F64`. The position representation moves not only a ULP figure but **which configurations the model accepts at all**. |
| Positive | `::test_the_kernel_contact_resultant_closure_holds_against_the_constituent_scale` | I5 / O4. |
| Positive | `::test_the_kernel_contact_matches_the_hand_computed_single_site_value` | O5. |
| Positive | `::test_every_connector_scatter_destination_is_allocated_through_the_backend` | I6, by type. |
| Positive | `::test_the_default_tether_step_state_visits_both_branches_and_both_step_kinds` | **I7 for the tether** — taut, slack, both-slack, gate-crossing, stationary, moving and ruptured counts are all `> 0` and are fields of the state. |
| Positive | `::test_the_default_contact_state_visits_both_the_engaged_and_the_separated_branch` | **I7 for the contact** — engaged and separated counts both `> 0`, and the tangential misregistration is non-zero so the normal-gradient term is live. |
| Positive | `::test_the_connector_report_names_the_position_precision_and_the_scatter_mode` | I8. |
| Positive | `::test_the_connector_references_are_driven_in_float64_positions` | I9. |
| Positive | `::test_a_non_positive_standoff_is_refused_before_any_kernel_launches` | **I10, host half** — the same `ConnectorGeometryError` type the reference raises. |
| Positive | `::test_a_standoff_the_host_calls_positive_and_the_kernel_calls_zero_is_refused_after_launch` | **I10, device half** — the kernel's own violation flag, on a configuration float64 calls legal. |
| Positive | `::test_the_position_precision_is_a_parameter_and_not_a_constant_for_these_laws` | All three modes reachable through the public builders and producing three different measurements. |
| Positive | `::test_the_declared_connector_budgets_sit_inside_what_roundoff_explains` | No declared budget exceeds its own round-off explanation. |
| Positive | `::test_the_state_stays_out_of_the_ambiguous_discrete_gradient_band` | §10's stated limit, asserted rather than hoped: every moving site moves by many float32 ULP of its own separation, so the `moved` predicate has the same value in float32 and float64. |

## 9. Deliberately failing negative control

Seven wrong kernels ship in `aleph/runtime/law_kernels.py` beside the right ones, compiled by the
same `load()` and launched down the same drivers. All are *plausible hand-translation errors*.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `::test_a_tether_kernel_that_forgot_the_discrete_gradient_is_caught` | `TetherStepVariant.WRONG_MIDPOINT_ONLY` always takes the midpoint force. **Bit-identical to the true kernel on a stationary configuration** — which is every configuration this vertical's pipeline evaluates — and caught only on a moving step. |
| Negative (must fail) | `::test_a_tether_kernel_that_ignores_rupture_is_caught` | `TetherStepVariant.WRONG_IGNORES_RUPTURE` drops the `bound` mask, so a broken linker keeps pulling. **Bit-identical on a state with no ruptured site.** |
| Negative (must fail) | `::test_a_contact_kernel_that_used_the_linear_law_is_caught` | `ContactVariant.WRONG_LINEAR_CORE` computes `(k/2) x^2` instead of the soft core — the substitution the reference calls a correctness regression, not a refinement. Caught on energy and on both forces. |
| Negative (must fail) | `::test_a_contact_kernel_that_dropped_the_half_in_its_energy_is_caught` | `ContactVariant.WRONG_ENERGY_HALF` drops the `1/2` from `U`. **Forces bit-identical.** Caught on the energy field alone. |
| Negative (must fail) | `::test_a_contact_kernel_that_dropped_the_unilateral_gate_is_caught` | `ContactVariant.WRONG_ADHESIVE_GATE` omits the `max(·, 0)`, so a separated contact **pulls**. This is the exact-zero breaker and the headline mutant of this entry. |
| Negative (must fail) | `::test_a_contact_kernel_that_froze_the_surface_normal_is_caught` | `ContactVariant.WRONG_FROZEN_NORMAL` omits the normal-gradient term from the cortex force, i.e. treats `n_i` as constant under cortex motion. **Energy bit-identical, membrane force bit-identical.** Caught on `forces_cortex_pn` alone and by O3. The direct analogue of G2's `WRONG_FROZEN_COTANGENT`, and the same class of real error. |
| Negative (must fail) | `::test_a_zero_ulp_budget_fails_the_true_connector_kernels` | The vacuity guard on the verdict: at `ulp_budget = 0.0` the *correct* kernels fail, because float32 storage is not bitwise float64. |
| Negative (must fail) | `::test_no_connector_mutant_escapes_through_the_atomic_floor` | Neither the widened `ATOMIC` budget nor the scatter floor is an amnesty for any of the six killable mutants. |

And the blind spots, recorded as controls because they are the honest limits of the method:

| Control | Location | Asserts |
|---|---|---|
| Blind spot (must hold) | `::test_the_dropped_discrete_gradient_is_invisible_on_a_stationary_configuration` | `WRONG_MIDPOINT_ONLY` is bit-identical to the true kernel when `end == start`. |
| Blind spot (must hold) | `::test_the_ignored_rupture_is_invisible_when_nothing_has_ruptured` | `WRONG_IGNORES_RUPTURE` is bit-identical when `bound` is all-True. |
| Blind spot (must hold) | `::test_the_adhesive_gate_is_invisible_on_an_all_engaged_configuration` | `WRONG_ADHESIVE_GATE` is bit-identical to the true kernel when every site is inside the standoff. **G1's all-taut blind spot in the compressive element**, and the reason the default contact state's separated-site count is a field a test asserts on. |
| Blind spot (must hold) | `::test_the_frozen_normal_is_invisible_in_the_energy_and_in_the_membrane_force` | `WRONG_FROZEN_NORMAL` is bit-identical on two of the three graded fields. An energy-only gate certifies it; so does a gate that grades only the side the integrator reads first. |

### 9a. The Newton-pair control this entry **does not** ship, and the measurement that decided it

`HANDOFF.md` §F-1 records that in an exactly linear drag law a `force_b := -force_a` mutant cannot be
killed by any measurement, because `-(zeta v)` and `zeta (-v)` are bit-identical in IEEE-754. The
build plan requires this lane to check whether the same holds for *these* laws before writing a
Newton-pair control.

**Measured on the frozen reference, 2,000 draws, before any kernel was written:**

| law | largest per-site `f^a + f^b` | `f^b` bitwise equal to `-f^a`? |
|---|---|---|
| `ErmTether`, fixed configuration | **0.0** | **yes** |
| `ErmTether`, moving step with 384 gate crossings in 2,000 | **0.0** | **yes** |
| `CortexMembraneContact` (`SOFT_CORE`, normal law) | 6.45 pN (force scale 91 pN) | **no** |

The reason is structural and not luck. The tether is central with a *scalar* magnitude: the reference
forms `f^b = F_s (b - a)/|b - a|`, and `|b - a|` is bit-identical to `|a - b|` because the
components are squared, while `(-d)/L` and `-(d/L)` are bit-identical because negation is exact. So
`f^b == -f^a` **to the last bit for every input**, and a mutant that wrote the mirror instead of
accumulating the second side independently is invisible to *any* measurement, at *any* precision.

> **So `TetherStepVariant.WRONG_NEWTON_MIRROR` ships as a declared blind spot, not as a killable
> mutant.** Its control is `::test_the_newton_mirror_mutant_is_provably_unkillable_for_the_tether`,
> which asserts it is bit-identical to the true kernel on every field of every shipped state — the
> honest statement — rather than a `pytest.raises`-shaped assertion that would have passed for the
> wrong reason. **A Newton-pair closure check is a real check only where the two sides can differ**,
> and for the contact they can: that is I5, and it is a genuine measurement of a genuinely
> non-mirrored pair.

## 10. Numerical and precision envelope

**Working precision.** Positions enter in the precision the case declares. Everything derived from
them — the site offsets, the separations, the engaged stretches, the face normals, the vertex normal
sums and unit normals, the standoff, the traction, the tangential misregistration, the normal
gradient and both assembled force fields — is **float32**, `ALEPH-DQ-107`'s compute channel.
Per-site **energies are widened to float64 inside the kernel** before the product is formed, and the
total is reduced by `Backend.sum`'s two-stage float64 accumulation. The references are float64 end to
end and I9 asserts it.

**The step-average force is formed in the compute channel, and that is a stated choice.** The
reference's `scalar_force_pn` computes `-(U(s1) - U(s0))/(s1 - s0)` in float64 throughout. The kernel
forms the quotient from its **float32** energies, because the quotient's output is a force and forces
are the compute channel; only the *reported* per-site energies are widened. §13.1 measures the cost.

**Position precision is a declared parameter, not a buried constant.** `PositionPrecision` is
`ALEPH-PORT-3603`'s enum, reused unchanged, and the reason is `ALEPH-PORT-3601` §14.7: the correct
tether kernel sits ~100 float32 ULP from its reference because of float32 **position storage**, an
order of magnitude above what the CUDA atomic scatter costs. It is raised in
`docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` and is **not settled here**.

**The design point that keeps the parameter cheap is preserved.** Exactly **two** kernels in this
entry read a position:

| kernel | what it forms | used by |
|---|---|---|
| `k_pair_offset_f32` / `_f64` | `d_i = a_i - b_i`, one row per site | both laws |
| `k_face_edges_f32` / `_f64` | `e_{f,k} = p_{k+1} - p_{k+2}` | the contact's cortex mesh — **`ALEPH-PORT-3603`'s kernel, reused, not re-written** |

Everything else in both laws is a function of those two outputs. `N_f = e_1 x e_2` reproduces the
reference's `(p1 - p0) x (p2 - p0)` bit for bit, because `p1 - p0 = -e_2` and `p2 - p0 = e_1` and
`cross(-e_2, e_1)` and `cross(e_1, e_2)` differ only by exact negations of exact subtractions. The
lever `p_{k+2} - p_{k+1}` is `-e_{f,k}`, likewise exact.

**One association-order difference from the reference, stated rather than discovered.** The
reference forms the tether's midpoint direction as `(a0 + a1)/2 - (b0 + b1)/2`, summing two
*position-scale* quantities before differencing them. The kernel forms it as
`((a0 - b0) + (a1 - b1))/2`, which is the same number in exact arithmetic and is **better
conditioned** — the cancellation happens once per configuration instead of once at the end. It is
therefore closer to the float64 reference than a literal transcription would be, and that is
recorded here so the parity figure is not read as evidence that literal transcription is fine.

### The ambiguous band, and why it is a limit rather than a bug

The `moved` predicate is `|s1 - s0| > 1e-13` µm. In float32 at the vertical's scale one ULP of a
`0.05` µm separation is `~4e-9` µm, so any real change is either **exactly zero** (when the two
configurations are bit-identical, in which case both precisions agree) or at least one float32 ULP,
which is `~4e4` times the floor. The band between `1e-13` µm and one float32 ULP is where the
predicate can take *different values on the two sides*, and where the difference quotient's own
cancellation is catastrophic: at a change of one ULP, `U1 - U0` cancels ~`10^6`-fold.

**This entry does not enter that band and does not pretend to grade it.** The state builder records
`min_moving_change_um` and the case asserts it is many float32 ULP above the separation's
resolution; §14 names the band as an unported region of the law rather than a solved one.

### Where the disagreement comes from, derived rather than assumed

Each position is rounded once before anything happens, so with `u = eps32/2`:

```
tether:   |d s_i|        <~ u (|a_i| + |b_i|)                         =: u m_i
          |d F_i|        <~ u [ k m_i (1 + 2 x_i / s_i) ]              stationary sites
          |d F_i|        <~ u [ (U0_i + U1_i + |F_i| m_i)/|change_i| + 2 |F_i| m_i / s_i ]   moving
contact:  |d s_i|        <~ u [ m_i + |d_i| * meshcancel ]            standoff, normal included
          |d U_i|        <~ u |U'(s_i)| |d s_i| / u
          |d f_i|        <~ u [ |U''(s_i)| |d s_i| / u + |U'(s_i)| * meshcancel ]
```

with `meshcancel = max_{f,k} (|p_a| + |p_b|) / |e_{f,k}|`, `ALEPH-PORT-3603`'s mesh cancellation
ratio, reused. `U'' = k r^3 / s^3` on the engaged branch and `0` elsewhere. Each field's
amplification is computed **from the host state by the case**, so it is a property of the data and
not a number chosen after seeing the answer; `roundoff_ulp_bound = (amp + 1)/2` converts it to ULP
exactly as `parity.py` already does, and `budget_within_roundoff_explanation` reports whether the
declared budget sits inside its own justification.

**The contact is the worst-conditioned law this project has put on a kernel, and the reason is
structural.** Its gate coordinate is a `0.02` µm standoff between surfaces `5` µm from the origin —
a `250`-fold cancellation before the normal is even applied — and it is then divided by `s^2` and
`s^3`. This is the tether's problem (§ `ALEPH-PORT-3601` §10) made sharper, and it is further
evidence for the proposal, not a defect of the kernel.

**The atomic floor.** `ATOMIC_SCATTER_ULP_FLOOR = 10.0` (`parity.py`, measured on the A5000). The
cortex normal sum and the assembled cortex force are scatter-assembled and take the floor in
`ATOMIC` mode; the energies are reduced by `Backend.sum`, which is ordered by construction, and do
not.

**Outside the envelope.** One device — **warp's CPU device**, on a Mac whose warp build reports
*"CUDA not enabled in this build"*. `ATOMIC` has never been measured here, and
`measure_scatter_determinism` already refuses to let a CPU zero be read as evidence about CUDA.

## 11. Production-backend residency and transfer

| Array | Where | Precision | Transfer |
|---|---|---|---|
| site / vertex positions | device | float32 or float64, per `PositionPrecision` | uploaded once per case by `backend.array` |
| cortex triangle indices | device | int32 | uploaded once per case (`wp.array`; not a scatter destination) |
| `bound` rupture flags | device | int32 | uploaded once per case |
| site offsets, face edges, face normals, unit normals, standoff, traction, tangential, coefficients | device | float32 | never leave the device |
| per-site energy | device | float64 | never leaves the device before `backend.sum` reduces it |
| cortex normal sum, assembled cortex force, assembled membrane force, assembled tether forces | device, **allocated by `backend.zeros`** | float32 | read to host once, by the harness |
| the contact's `violation` flags | device, `wp.array` | int32 | read to host once per launch, by the driver, to decide the refusal |

**Four** scatter destinations — the two tether force fields, the cortex normal sum, and the cortex
force — which is four chances to make the mistake `636b0c8` retracted. I6 asserts by type, by
watching every `scatter_add` call rather than by reading the source.

The `violation` array is allocated with `wp.array` and **not** through the backend, for the reason
`_device_triangles` already gives: `WarpBackend` accepts float32 and float64 only, and adding an
int32 dtype there means adding a kernel instantiation and a parity row for it. It is not a scatter
destination, which is the thing §11's rule is about.

Two `scatter_add` calls write into the **same** destination — the site kernel writes the cortex
force's direct term into `forces_cortex` and the gradient stage's corner rows are then accumulated
on top of it. `WarpBackend.scatter_add` adds into `dest` (and its ORDERED path accumulates each
row's contributions in float64), so this is one array assembled by two passes rather than two
arrays summed afterwards.

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
| Acceptance result | Measured 2026-07-31 on **warp's CPU device** (this machine's warp build reports *"CUDA not enabled in this build"*), warp 1.15.0, float32 compute, `ScatterMode.ORDERED`, `PositionPrecision.GLOBAL_F32` unless stated. `tests/runtime/test_connector_law_parity.py`: **36 passed**. `tests/runtime` as a whole: **338 passed, 3 skipped** (was 302 + 3 before this entry). `tests/vertical`: **1,358 passed**, unchanged — the frozen references were read and never edited. Mutation study: **8 of 8 killed** (§13a). **The exact-zero inactive branch survived the port on every branch of both laws, asserted with `==` per site**, and an all-separated contact passes the parity gate at a declared budget of `0.0` ULP. The status stays `PROPOSED` because no CUDA device was touched and no human has reviewed it. |
| Reviewer | **Agent-proposed. Unratified.** Written by a subagent of session `46143f30`; no PI review, no `decided_by` field, and none may be added by an agent. |
| Rollback | Revert the appended blocks of `aleph/runtime/law_kernels.py` and `aleph/runtime/law_cases.py` and delete `tests/runtime/test_connector_law_parity.py`. `ALEPH-PORT-3601`'s and `-3603`'s paths are untouched by construction and keep working. Nothing under `aleph/vertical/` changes in either direction. |

### 13.1 The tether — the full law, five branch families in one state

Default state: 96 sites over 32+32 vertices — **55 taut / 41 slack at the start, 54 / 42 at the end,
30 slack at *both* ends, 23 crossing the gate, 72 moving / 24 stationary, 14 ruptured.**

| kernel | `energy_start` | `energy_end` | `forces_a` | `forces_b` | verdict |
|---|---|---|---|---|---|
| **TRUE** | **0.91** | **4.60** | **55.88** | **39.44** | pass |
| `WRONG_MIDPOINT_ONLY` | *bit-identical* | *bit-identical* | **666,240.99** | **431,405.62** | fail — **forces only** |
| `WRONG_IGNORES_RUPTURE` | *bit-identical* | *bit-identical* | **4,909,079.44** | **3,969,428.14** | fail — **forces only** |
| `WRONG_NEWTON_MIRROR` | *bit-identical* | *bit-identical* | *bit-identical* | *bit-identical* | **passes — see §9a** |
| declared budget | 32 | 32 | 192 | 192 | — |
| round-off explains | 410.81 | 333.78 | 435.07 | 435.07 | — |

On the two degenerate states, which is where the blind spots live:

| state | TRUE | `WRONG_MIDPOINT_ONLY` | `WRONG_IGNORES_RUPTURE` |
|---|---|---|---|
| stationary (`end == start`) | 0.91 / 0.91 / 100.12 / 42.39 | ***bit-identical*** | 7,270,359 / 5,747,059 |
| all-taut + unbroken | 10.29 / 3.03 / 26.68 / 28.63 | 143,833 / 114,343 | ***bit-identical*** |

**The first row of that table is the finding worth carrying.** `_DistributedUnilateralConnector.accumulate`
calls `site_pair_forces` with **no end positions**, so `end == start` at every step of every run this
vertical has ever done. A parity case built from the law's own production use is therefore *provably
blind* to a kernel that dropped the discrete gradient — the rule which exists specifically to close
the energy books across a gate crossing. This is G1's all-taut blind spot, except that the blind
configuration is the shipped one.

### 13.2 The contact — the worst-conditioned law this project has put on a kernel

Default state: a level-2 jittered cortex, **162 vertices / 320 faces, 70 engaged / 92 separated**,
standoff ∈ [0.0080, 0.0398] µm against `r = 0.020` µm, `k = 2300` pN/µm, mean tangential
misregistration 0.221 µm.

| kernel | `energy_pn_um` | `forces_membrane_pn` | `forces_cortex_pn` | verdict |
|---|---|---|---|---|
| **TRUE** | **49.51** | **539.73** | **496.24** | pass |
| `WRONG_LINEAR_CORE` | 4,319,959.44 | 6,382,083.17 | 6,386,326.98 | fail — all three |
| `WRONG_ENERGY_HALF` | **8,388,508.97** | *bit-identical* | *bit-identical* | fail — **energy only** |
| `WRONG_ADHESIVE_GATE` | 5,899,620.95 | 1,316,356.33 | 1,246,799.85 | fail — all three |
| `WRONG_FROZEN_NORMAL` | *bit-identical* | *bit-identical* | **559,000.46** | fail — **cortex force only** |
| declared budget | 256 | 2,048 | 2,048 | — |
| round-off explains | 1,871.21 | 3,924.80 | 9,486.60 | — |

All-engaged state (162 engaged / 0 separated): TRUE 15.18 / 874.98 / 827.93 against bounds
1,870.35 / 4,666.00 / 11,306.01, and **`WRONG_ADHESIVE_GATE` is bit-identical to TRUE there** —
the compressive element's version of the all-taut blind spot.

All-separated state (0 engaged / 162 separated): **every field exactly `0.00`, and the gate passes
at a declared budget of `0.0` ULP.**

Five things to read off these two tables.

1. **Three of the seven wrong kernels are invisible in some channel, and in different ones.**
   `WRONG_ENERGY_HALF` is energy-only; `WRONG_FROZEN_NORMAL` is cortex-force-only; every tether
   mutant is force-only. A gate reporting one aggregated number certifies at least one of them, and
   a gate grading only the membrane side — the side carrying the turgor, the one a reader reaches
   for first — certifies `WRONG_FROZEN_NORMAL`.
2. **Every killable mutant clears its budget by three to six orders of magnitude**, so the 2.3–5×
   margin between measurement and declared budget is nowhere near the detection threshold.
3. **Every declared budget sits below what round-off explains on every state it grades**, asserted
   by `test_the_declared_connector_budgets_sit_inside_what_roundoff_explains`.
4. **The contact's ~500 ULP is not a defect and is not the scatter.** Its gate coordinate is a
   0.02 µm standoff between surfaces 5 µm from the origin — a 250-fold cancellation before the
   normal is applied — and it is then divided by `s^2` and `s^3`. §13.4 separates the terms.
5. **`ORDERED` vs `ATOMIC` is again not the binding constraint.** On warp's CPU device the contact's
   true kernel measures 47.89 / 526.15 / 478.83 under `ATOMIC` against 49.51 / 539.73 / 496.24 under
   `ORDERED`, and the tether's figures are unchanged; every killable mutant still fails with the
   atomic floor in force. (The CPU device's atomic path is reproducible by accident of the schedule;
   `measure_scatter_determinism` already refuses to let that be read as evidence about CUDA.)

### 13.3 `F = -grad E` — the step read off a curve, not off a formula

`PLAN.md` §6.1 records this project losing a night to central-difference steps below the round-off
floor. So both curves were measured before either control was written.

**The contact, both force fields, relative error at four probe vertices (`POSITIONS_F64`):**

| `h` [µm] | membrane | ratio | cortex | ratio |
|---|---|---|---|---|
| 3.0e-3 | 1.457e-01 | — | 1.503e-01 | — |
| 1.0e-3 | 1.461e-02 | **9.98** | 1.507e-02 | **9.98** |
| 3.0e-4 | 1.288e-03 | **11.34** | 1.327e-03 | **11.35** |
| 1.0e-4 | 1.603e-04 | 8.03 | 2.222e-04 | 5.97 |
| 3.0e-5 | **9.68e-05** | 1.66 | **1.894e-04** | 1.17 |
| 1.0e-5 | 5.32e-04 | **0.18** | 1.138e-03 | **0.17** |
| 3.0e-6 | 9.19e-04 | 0.58 | 4.034e-03 | 0.28 |

Order 2 is 9.0 and 11.1 for step reductions of 3 and 3.33; measured **9.98 and 11.34**. The curve
**bottoms out at `h = 3e-5` and then rises**, which is the float32 round-off floor of the energy
channel divided by an ever smaller `h`. The step the control uses is **`h = 3e-4`**, inside the
window the curve establishes. Under `GLOBAL_F32` the same curve floors an order of magnitude higher
(1.06e-3 at `h = 3e-4`, ratios 10.04 then 13.69), which is why the gradient oracle declares
`POSITIONS_F64`: the order-2 demonstration would otherwise have been taken *at* the floor.

`WRONG_FROZEN_NORMAL` measures **3.5%** against the same numeric derivative — 27× the true kernel's
1.3e-3 — so the oracle sees the term the mutant deletes.

**The tether, on the stationary state (`GLOBAL_F32`), where the discrete gradient reduces to
`-U'(s)`:**

| `h` [µm] | 1.0e-2 | 3.0e-3 | 1.0e-3 | 3.0e-4 | 1.0e-4 | 3.0e-5 | 1.0e-5 |
|---|---|---|---|---|---|---|---|
| relative | 4.449e-03 | **4.061e-04** | 8.459e-05 | 2.027e-04 | 2.146e-04 | 2.449e-03 | 9.359e-03 |
| ratio | — | **10.96** | 4.80 | 0.42 | 0.94 | 0.09 | 0.26 |

10.96 against 11.1 for exact second order, and a floor from `h = 3e-4` down. Step used: **3e-3**.

### 13.4 Position precision — a third law, and the same two-term structure

Assembled contact force ULP, true kernel, by mode and by where the body sits:

| mode | body at the origin | body at (50, 0, 0) µm |
|---|---|---|
| `GLOBAL_F32` | 539.73 / 496.24 | **3,041.18 / 2,937.84** |
| `LOCAL_F32` | 445.34 / 430.41 | **445.34 / 430.41** |
| `POSITIONS_F64` | **17.19 / 16.79** | **17.19 / 16.79** |

(membrane / cortex; the energy runs 49.51 → 145.22 for `GLOBAL_F32`, 29.25 for `LOCAL_F32` and
1.34 for `POSITIONS_F64`, unchanged by the offset in the last two.)

**This confirms in a third law the structure `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md`
argues from.** `LOCAL_F32` removes the origin-distance term *exactly* — 3,041 → 445 — and buys about
18% when the body is already centred. `POSITIONS_F64` removes both terms and is 26× better than
`GLOBAL_F32` at the origin, which is a larger gap than either the tether or the membrane showed,
because this law divides the cancelled quantity by its own square and cube.

### 13.5 The refusal a kernel cannot raise, and a finding for the open decision

`ALEPH-PORT-3604` §5's decision, measured. Sites placed at a standoff in the band where float64 and
float32 disagree about the sign, 162 sites, jittered level-2 cortex:

| standoff placed at | `GLOBAL_F32` | `LOCAL_F32` | `POSITIONS_F64` |
|---|---|---|---|
| 1e-8 µm | **72 sites flagged** | 78 flagged | 3 flagged |
| 1e-7 µm | **33 flagged** | 37 flagged | **runs clean** |
| 1e-6 µm | clean | clean | clean |

Every one of those configurations is **positive in float64 on the host**, so the host predicate
passes and the launch happens; the kernel's own arithmetic then disagrees, flags, writes exactly
zero for the offending sites, and the driver raises `ConnectorGeometryError`.

> **The finding.** At a 1e-7 µm standoff the *same* configuration is **refused under `GLOBAL_F32`
> and accepted under `POSITIONS_F64`.** The position representation therefore moves not only a ULP
> figure but **which configurations the model will accept at all** — a different kind of
> consequence from the ones §3 and §3b of the proposal record, and one that is not in it yet.
> `test_the_position_precision_changes_whether_the_geometric_refusal_fires` is the measurement.
> This lane does not take the decision; it adds one more thing that hangs on it.

`LOCAL_F32` being marginally *worse* than `GLOBAL_F32` here (78 against 72) is not a defect: the
body is already centred on the origin in these states, so the centroid subtraction is a small
translation with no cancellation to remove, and the difference is which sites happen to land on
which side of a float32 rounding.

### 13.6 The exact anchors

| oracle | measured through the kernel | reference / expected |
|---|---|---|
| tether, both ends slack | **exactly `0.0`** on all 30 sites, both force fields, both energies | exact |
| tether, ruptured | **exactly `0.0`** on all 14 sites, and 8 of them are taut | exact |
| contact, separated | **exactly `0.0`** on all 92 membrane sites | exact |
| contact, all-separated state | **every field exactly `0.0`**; gate passes at budget `0.0` | exact |
| contact resultant closure | **1.220e-05 pN** against a constituent scale of 3,151.9 pN (3.9e-09 relative) | exact identity |
| contact per-site `f^m + f^c` | 6.45 pN against a 91 pN force scale — **the two sides genuinely differ** | not a Newton pair |
| tether per-site `f^a + f^b` | **bitwise zero, 2,000 draws incl. 384 gate crossings** | §9a: unkillable |
| hand-computed contact site, `s = 0.010`, `r = 0.020`, `k = 2300` | `U = 0.23` pN·µm/site and a membrane force magnitude of `69.0` pN, matched to 1e-4 / 1e-3 relative | derived in §7 (O5) |
| `F = -grad E`, contact, both sides | ratios **9.98 / 11.34** | order 2 |
| `F = -grad E`, tether, stationary | ratio **10.96** | order 2 |

### 13a. Mutation study — 8 planted defects, 8 killed

Run with `PYTHONDONTWRITEBYTECODE=1` (a same-length mutant with a sub-second edit/revert leaves a
`.pyc` CPython considers valid, which has already contaminated one study in this repository).
Baseline green, each mutant applied alone, reverted before the next, green again after.

| # | Planted defect | Killed by |
|---|---|---|
| M1 | the contact driver ignores the requested variant and always launches the true stages | 7 tests |
| M2 | the tether driver launches `tether_step_true` whatever was asked | 3 tests |
| M3 | a contact scatter destination allocated as a host ndarray instead of `backend.zeros` | 21 tests |
| M4 | `LOCAL_F32` stops subtracting the reference point (silently becomes `GLOBAL_F32`) | 1 test |
| M5 | `POSITIONS_F64` silently uploads float32 and launches the float32 offset kernel | 3 tests |
| M6 | **the kernel's own violation flag is never read, so the device-side refusal is dropped** | 2 tests |
| M7 | **the frozen contact reference is handed float32-rounded positions** | 2 tests |
| M8 | the tether's rupture mask is uploaded as all-bound | 7 tests |

**M6 is this entry's own contribution to the list.** It is the failure mode the §5 decision exists to
prevent: the host refuses in float64, the kernel disagrees, and if nobody reads the flag the run
continues on a configuration the law says is impossible — with exactly zero force at those sites,
which looks like a separated contact rather than like a fault. It is killed by
`test_a_standoff_the_host_calls_positive_and_the_kernel_calls_zero_is_refused_after_launch` and by
`test_the_position_precision_changes_whether_the_geometric_refusal_fires`.

**M7 is `ALEPH-PORT-3603`'s M7 repeated deliberately**, because it is the hazard
`docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` §4 names and it is invisible in
the artefact: a reference quietly given the same float32 input as the kernel loses the same digits,
and the gate prints a small healthy number.

**M4 is killed by exactly one test, and that is a finding rather than a comfort.** A declared
parameter that only one control can see the effect of is one deletion away from being decoration.

## 14. Honest limits

What this entry does **not** establish:

1. **No CUDA device was touched.** Every number is from **warp's CPU device**, where a launch is a
   serial loop. It is evidence that the kernels' algebra matches the reference and that the
   identities hold; it is `UNVERIFIED` as evidence about CUDA.
2. **The ambiguous band of the discrete gradient is not ported.** §10 defines it. A step whose
   separation change is between `1e-13` µm and one float32 ULP of the separation is a configuration
   on which the `moved` predicate can differ between the reference and the kernel, and on which the
   quotient's float32 cancellation is catastrophic. The state stays out of it and a control asserts
   that; nothing here is evidence about a run that enters it.
3. **The `ContactLaw.LINEAR` branch of `CortexMembraneContact` is not ported.** That branch is the
   deliberately-adhesive negative-control variant (`adhesive=True`), which routes back through the
   central `site_pair_forces` path the tether kernel already covers. The **default** `SOFT_CORE`
   normal law is what is ported, which is what a run uses.
4. **No unbinding kinetics.** `bound` is a flag nothing in this vertical flips, and the kernel
   reads it rather than evolving it. Nothing here is evidence about tether lifetime or blebbing.
5. **The negative controls are as strong as the state they run on, and four of them are provably
   blind somewhere** — §9's blind-spot table names each one and the configuration that hides it.
   This method cannot tell you that a case's inputs visit every branch of a law; the improvement
   over G1 is that the branch counts are **fields of the state that a test asserts on**, not that
   the judgement was removed.
6. **A Newton-pair control for the tether is impossible, not omitted** (§9a). That is a limit of the
   law's own structure and it holds at every precision.
7. **The declared budgets are engineering numbers and they belong to one state each.** The
   cancellation grows with the distance of the body from the origin and with mesh refinement, so a
   case on another state must declare its own budget.
8. **The energy is compared as a scalar total.** A kernel wrong on two sites in opposite directions
   would cancel and pass. Inherited from `ALEPH-PORT-3601` §14.6 and `-3603` §14.7, named again
   rather than assumed fixed. The exact-zero check (I2) is per-site and is not subject to this.
9. **No relaxation is on device and none is proposed**, per `ALEPH-PORT-3603` §14.6. No oracle here
   needs one: every identity used holds at an arbitrary configuration.
10. **Nothing here measures throughput.** The cost of the four-stage contact chain against a fused
    kernel, and the FP64 cost of `POSITIONS_F64`, are unmeasured on any device.
11. **Evidence class.** `STRUCTURAL` for the parity figures; `ANALYTIC_ORACLE` for the exact zero,
    the resultant closure and `F = -grad E`. No energy or force in this entry may be reported as a
    property of a cell.
