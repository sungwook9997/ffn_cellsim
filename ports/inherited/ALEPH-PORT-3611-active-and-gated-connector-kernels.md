# ALEPH-PORT-3611 — the last eleven wired connector rows as Warp kernels

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3611` |
| Lane | `46143f30` Lane G5c (CUDA track), Track G task **G5c** |
| Status | `PROPOSED` |
| Written | `2026-08-01` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Depends on | `ALEPH-PORT-3601` (the `LawCase` contract and the law-level parity harness), `ALEPH-PORT-3604` (`k_pair_offset_*`, the per-site launch shape, the exact-zero discipline), `ALEPH-PORT-3606` (the discrete-gradient step kernel and the rows-vs-kernels census discipline), `ALEPH-PORT-3610` (the *rate*-channel pattern for a law with no stored energy, and the two-count reporting rule). All four are **reused, never duplicated**. |
| Exists because | `-3610` left **11 of the 30 wired rows** unported and named the four classes that carry them. This entry takes all four. Two of them are springs the harness already knows how to launch; **two of them are not force laws of a configuration at all** — their force is set by internal state and the position blocks they are handed are *verified rather than used*. That is a structurally new thing to grade and §4.0 states what a parity gate can and cannot say about it. |

---

## 1. Aleph API

```python
from aleph.runtime.law_kernels import (
    SeriesJointVariant,       # TRUE | WRONG_ENERGY_HALF | WRONG_UNILATERAL_CLAMP
                              #      | WRONG_PARALLEL_STIFFNESS | WRONG_OPEN_CLUTCH_LEAKS
    CableVariant,             # TRUE | WRONG_OMIT_RESTIFFENING | WRONG_TOE_COMPRESSION
                              #      | WRONG_ENERGY_BRANCH_OFFSET
    RatchetVariant,           # TRUE | WRONG_LINEAR_LOAD | WRONG_LOAD_PER_SITE | WRONG_SIGN
    CrossbridgeVariant,       # TRUE | WRONG_ENERGY_HALF | WRONG_ABS_STRAIN
                              #      | WRONG_STRAIN_FROM_SEPARATION
    load, loaded_warp,        # unchanged, extended to compile these kernels
)

from aleph.runtime.law_cases import (
    ACTIVE_AND_GATED_BATCH,               # the batch census: class -> rows -> kernel family
    DEFAULT_SERIES_JOINT_ULP_BUDGET,
    DEFAULT_CABLE_ULP_BUDGET,
    DEFAULT_RATCHET_ULP_BUDGET,
    DEFAULT_CROSSBRIDGE_ULP_BUDGET,
    active_batch_connector,               # build one member through ITS OWN frozen builder
    series_joint_state,                   # engaged / open / compressed clutches, placed
    series_joint_law_case,
    series_joint_kernel_site_forces,
    strain_stiffening_cable_state,        # every branch of the card visited, counted
    strain_stiffening_cable_law_case,
    strain_stiffening_cable_kernel_site_forces,
    brownian_ratchet_state,               # v = 0, v = v0, and the interior, counted
    brownian_ratchet_law_case,
    brownian_ratchet_kernel_site_forces,
    crossbridge_state,                    # bound heads at mixed-sign strain, one exactly zero
    crossbridge_law_case,
    crossbridge_kernel_site_forces,
    refuse_non_unit_axes,                 # the host refusal, on the LAUNCH path — see §5
    refuse_nonpositive_series_geometry,   # the host refusal, on the LAUNCH path — see §5
)
```

Nothing outside this list is covered. This entry does **not** authorise any change to
`aleph/vertical/**` — `connectors_adhesion.py`, `connectors_frame.py`, `connectors_protrusion.py`,
`nmii.py`, `intermediate_filament.py`, `focal_adhesion.py` and `wiring.py` are the frozen parity
reference and are **read, never edited** — nor to `aleph/runtime/parity.py`, `backend.py`,
`warp_kernels.py`, `residency*.py`, `aleph/scenarios/**`, `aleph/viz/**`, `aleph/vertical/relax.py`,
`tests/runtime/test_reference_input_guard.py`, or `ports/ledger/INDEX.md`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read by this lane.** No file under that tree was opened. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

The laws being ported are Aleph's own:
`aleph.vertical.connectors_adhesion.SeriesJointConnector.evaluate_sites`,
`aleph.vertical.connectors_frame.StrainStiffeningCable.{energy_pn_um,derivative_pn,scalar_force_pn}`
composed through `aleph.vertical.connectors.site_pair_forces`,
`aleph.vertical.intermediate_filament.{_branch_widths,tension_pn,storage_energy_density_pn}`,
`aleph.vertical.connectors_protrusion.{BrownianRatchet.load_pn,BrownianRatchetConnector.evaluate_sites}`,
and `aleph.vertical.nmii.{NmiiPopulation.head_axial_force_pn,NmiiMotorConnector.evaluate_sites}`.
This entry moves already-ported laws onto kernels; it re-ports nothing.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported from another project.** As in `-3603`, `-3604`, `-3606` and
`-3610`: the purpose of a parity gate is agreement with the frozen module, so the association order
of every product is fixed by that module rather than chosen here.

Two instances are load-bearing rather than stylistic and are called out where they arise:

1. **`intermediate_filament._branch_widths` writes the three-branch law as three *clipped extents*
   `(a, b, c)` and not as a chain of `if` branches**, and its own docstring says why: exactly one of
   the three varies at any strain, so `dW/de` telescopes to `T` with no leftover term and the energy
   is the *exact* integral of the tension rather than an approximation of it. A clean-room kernel
   would almost certainly write `if e <= ek: ... elif e <= ep: ...`, which is the same number in
   exact arithmetic and destroys the property the entry's gradient oracle rests on. The kernel is
   written with the clipped extents.
2. **`SeriesJointConnector` reads the joint compliance from `ClutchEndpoint.series_stiffness_pn_per_um`
   and never from the two half stiffnesses.** The halves are published by `focal_adhesion.py` as a
   *diagnostic*, documented in both modules as not being a load path. The true kernel therefore
   takes the gated series stiffness and **is not given the halves at all**; the two impostors of §9
   take them as an extra argument they cannot be written without. §4.1.

## 4. Physical or mathematical law represented

### 4.0 The scope decision, the criterion, and the two counts

Task G5c is "the remaining wired connector rows". The unit is again **not** position in the wiring
table. `-3610` left four classes; this entry takes all four, split by a criterion a test can
**measure** rather than assert:

> **Does `evaluate_sites` compute its force from the two position blocks it is handed?**

* **Batch F — yes.** `SeriesJointConnector` and `StrainStiffeningCableLink`. Perturb either block
  and the force moves. These go on the existing `k_pair_offset_*` stage.
* **Batch G — no.** `BrownianRatchetConnector` and `NmiiMotorConnector`. The force is set by
  internal state — a growth velocity and a head strain respectively — and both classes' own
  docstrings say so in as many words ("the blocks are verified rather than used"; "a motor is not a
  spring between the two arrays it is handed"). Perturb a block within what the connector accepts
  and **the force does not move at all**.

`test_batch_f_reads_its_configuration_and_batch_g_verifies_it` asserts both halves by *evaluating*:
for batch F a position perturbation changes the force; for batch G every position the connector
accepts gives the same force, and a position it does not accept **raises**. It is disjoint by
construction from `-3606`'s and `-3610`'s batches and
`test_the_three_connector_batches_are_disjoint` asserts that against `CENTRAL_PAIR_BATCH` and
`DISSIPATIVE_TRANSFER_BATCH`.

The batch selects **four classes carrying the last eleven of the thirty wired registry rows**:

| class | module | law | rows | kernel family |
|---|---|---|---|---|
| `SeriesJointConnector` | `connectors_adhesion` | gated two-sided spring on a **series** compliance | **4** | **J (new)** |
| `StrainStiffeningCableLink` | `connectors_frame` | tension-only **three-branch** card law | **1** | **K (new)** |
| `BrownianRatchetConnector` | `connectors_protrusion` | logarithmic **load–velocity** relation | **2** | **L (new)** |
| `NmiiMotorConnector` | `nmii` | crossbridge `k_xb x` in the **head strain** | **4** | **M (new)** |

> ### ⚠ ROWS AND KERNELS ARE DIFFERENT NUMBERS AND QUOTING EITHER ALONE MISREPRESENTS THE PORT
>
> `wiring.py` counts **rows**. This entry writes **kernels**. The honest report is
> **11 registry rows, 4 classes, 4 new kernel families, 17 new kernels (4 true, 13 wrong on
> purpose), and 0 rows carried without a kernel.**
>
> `test_the_batch_covers_the_last_eleven_wired_registry_rows` asserts both numbers against
> `wiring.WIRED` directly, and
> `test_the_ported_row_census_is_now_complete_against_the_wiring_registry` asserts the running total
> — **30 of 30** — as a difference against the three batch tables rather than as a sentence.

### 4.1 Family J — the gated series joint. Four rows.

With `x_a` the actin anchor block, `x_l` the `ecm` ligand block, `d0` the joint's unloaded gap and
`k` the joint's **series** stiffness:

```
d   = |x_l - x_a|                    joint separation [um]
e   = d - d0                         extension [um] — SIGNED, there is no clamp
T   = k e                            joint tension [pN]
U   = sum_i (k_i/2) e_i^2            stored energy [pN.um]
f^(a) = T (x_l - x_a)/d              f^(b) = T (x_a - x_l)/d
```

Three facts decide what a gate on this law can do:

1. **The gate is an exact zero that arrives as data.** `ClutchEndpoint.series_stiffness_pn_per_um`
   returns literally `0.0` for an open clutch, so `T = 0.0 * e` is the number `0.0` for an algebraic
   reason and every component of both force blocks is exactly `0.0`. `-3604`'s strongest oracle
   exists here, on a **branch of the population** rather than of the geometry, and it is asserted
   with `==` per site.
2. **The compliance is a series combination and the halves are not a load path.** `k` is
   `k_a k_l/(k_a + k_l)`, always softer than either member. The two impostors of §9 are exactly the
   two ways a reader reaches for the halves — the parallel combination `k_a + k_l`, and the gate
   ignored — and **neither can be written without an argument the true kernel does not have.**
3. **The law is two-sided.** `extension` is signed; a compressed clutch pushes. The unilateral clamp
   is therefore available as a mutant and is **bit-identical to the true kernel on an all-stretched
   state**, which is `-3601` §14.7's rule arriving in a fifth law. The default state's compressed
   count is a field a test asserts on.

### 4.2 Family K — the three-branch strain-stiffening cable. One row.

On the bond's own rest length `s0`, with the card `(E_toe, ek, E_pl, ep, E_st)`:

```
e       = (s - s0)/s0                                dimensionless — the CARD's coordinate
a       = clip(e, 0, ek)                             toe extent
b       = clip(e - ek, 0, ep - ek)                   plateau extent
c       = max(e - ep, 0)                             re-stiffening extent
T(e)    = E_toe a + E_pl b + E_st c                  [pN] — exactly 0.0 for e <= 0
W(e)    = E_toe a^2/2 + T_k b + E_pl b^2/2 + T_p c + E_st c^2/2
U(s)    = s0 W(e)                                    [pN.um]
dU/ds   = T(e)                                       exactly, by the telescoping of §3
```

driven over a step by `UnilateralSpring.scalar_force_pn`'s discrete gradient with the reference's
own `_DISCRETE_GRADIENT_FLOOR_UM`, the midpoint fallback, the rupture mask, and the two directed
normalisations — all four **inherited from `-3606`'s kernel shape**, not re-decided.

`T_k = E_toe ek` and `T_p = T_k + E_pl (ep - ek)` are the branch **offsets**: they are what makes
`W` continuous at the two knees. Dropping them integrates each branch from zero and produces an
energy that is wrong everywhere past the knee while leaving `T` — and therefore every force —
bit-identical. That is §9's energy-only mutant and it is a very ordinary transcription slip.

### 4.3 Family L — the Brownian ratchet. Two rows.

The first **transcendental** in this module, and the first law whose independent variable is a
growth velocity:

```
arg   = eps + (1 - eps) v/v0                  eps in (0,1) strictly, v > v_min
f_bar = -(kT/d) log(arg)                      load per barbed end [pN]
F_i   = n_i f_bar_i                           site load [pN], n_i barbed ends share it
f^(a) = -F axis                               f^(b) = +F axis            U = 0 EXACTLY
```

Three structural facts:

1. **`U` is identically zero over the whole law**, as in `-3610`'s family D, and it is again *cheap
   evidence*: every mutant here preserves it. The rate that plays the accounting role is the
   **active power** `P = sum_i F_i v_i`, which vanishes at *both* ends of the curve and is strictly
   positive between them. That pair of facts is the physical content of a rectifier.
2. **`v = v0` gives `log(1) = 0` and therefore a load of exactly `-0.0` in the reference.** An exact
   anchor that is not a gate and not a mask — it is the transcendental's own fixed point. **The
   float32 kernel loses it by one ULP of association** and §14.10 is that finding, including the
   reassociation that would restore it and the reason this lane does not take it.
3. **`f_bar(0) = (kT/d) ln(1/eps)` is the stall force, and it contains no rate.** Scaling `v0` by
   three decades leaves it unchanged, which is the oracle a leading edge that had degraded into a
   stiffness cannot pass. §7 O3.

### 4.4 Family M — the crossbridge. Four rows.

```
axial_h = k_xb x_h                            over the BOUND heads of this connector [pN]
f^(a)   = +axial t_hat                        f^(b) = -axial t_hat
U       = (k_xb/2) sum_h x_h^2                [pN.um]
```

`x_h` is the head's own strain, advanced by the population's kinetics — **not** a function of the
anchor/site separation. Two heads at identical separations carry different forces if they are at
different points in the stroke, and the reference refuses a configuration that is not the one its
strains were advanced against. §9's headline mutant for this family is exactly the confusion the
reference's docstring names: computing `k_xb(|a - b| - rest)` from the separation. **It reads a
position the true kernel never touches**, so the position-precision decision reaches the impostor
and not the law.

`x_h = 0` gives exactly `0.0` force and contributes exactly `0.0` to the energy. Asserted with `==`.

### 4.5 The kernel is not the whole connector, in two places, and both are stated rather than implied

* **Family J does not port `ClutchEndpoint.series_stiffness_pn_per_um`.** The series reduction, the
  occupancy weighting, the maturation scale and the stale-handle refusal live in `focal_adhesion.py`
  and stay on the host. The kernel receives the gated stiffness as data. Porting it would move a
  second owner's arithmetic into a connector's kernel.
* **Family M does not port `NmiiPopulation.propose_kinetics`.** Attachment, detachment, the stroke
  and the power budget are the motor's kinetics and are not `evaluate_sites`. The kernel receives
  the strains and the tangents. A kernel that advanced a stroke would be a different port.

## 5. Units, domains, singular cases, invariants

Units throughout: length **µm**, force **pN**, energy **pN·µm**, time **s**, stiffness **pN/µm**,
velocity **µm/s**, power **pN·µm/s**. No other unit appears in any argument or return.

| refusal | where the frozen reference does it | where this entry does it |
|---|---|---|
| coincident joint ends (`d <= 0`) | `SeriesJointConnector.evaluate_sites` → `AdhesionCouplingError` | `refuse_nonpositive_series_geometry`, in float64 on the **launch path** |
| non-unit ratchet axis (> 1e-12) | `BrownianRatchetConnector.set_axes` → `ProtrusionConnectorError` | `refuse_non_unit_axes`, in float64 on the **launch path** |
| growth velocity off the branch (`arg <= 0`) | `BrownianRatchet.load_pn` → `RatchetDomainError` | same function, called by the state builder **and** by the launch path |
| non-unit actin tangent (> 1e-12) | `ActinBindingSites.__post_init__` → `NmiiGeometryError` | inherited: the state is built through the frozen class |

**A Warp kernel cannot raise.** So each refusal happens on the host, in float64, *on the path a
launch takes* and not only on the path the state builder takes — `-3606` §5 and `-3610` §5 both
record that hole being real in a first draft, and
`test_the_host_refusals_run_on_the_launch_path` asserts it by handing a hand-assembled state
mapping straight to the per-site helper.

Invariants this entry checks rather than assumes:

* **I1** `U` and both force blocks are exactly `0.0` at every open clutch (family J), on every slack
  site (family K), and at `x_h = 0` (family M). All with `==`, all with a non-emptiness assertion on
  the mask. **Family L's anchor at `v = v0` is exact in the reference and is LOST by the kernel**,
  by exactly one float32 ULP of the argument's association — measured, not assumed, and §14.10 is
  the finding rather than a repair.
* **I2** `f^(b) = -f^(a)` **bitwise** in all four families. §9a: this makes a Newton-pair control
  structurally unable to fail, and **no mirror mutant and no closure control ships for any of the
  four.** The bit-identity is asserted instead.
* **I3** family K's `dU/ds` equals `T(e)` to round-off at every strain, in both directions across
  both knees — the telescoping of §3.1.
* **I4** family L's active power vanishes at `v = 0` and at `v = v0` and is strictly positive
  between.
* **I5** family J's `k` is strictly below both halves at every engaged joint — the series
  combination is softer than either member, and the parallel one is stiffer than both.

## 6. Source evidence class and known retractions

Everything here is `ANALYTIC_ORACLE` rung and `BLOCKED` for quoting. No parameter in any state
builder is a cell measurement: the ratchet card, the IF material card, the clutch card and the motor
kinetics are all fixture values chosen to make the algebra sharp and to place a site on every branch.
They are labelled as such in the module and none acquires an evidence class by being used here.

Three limits are inherited from the frozen references and are **not** repaired by this entry:

* `BrownianRatchetConnector.snapshot`/`commit`/`rollback` are empty: the class ships **no kinetics**,
  so nothing adds a subunit and the growth velocity is set by the caller each step. The registered
  contract is `[K]` with `commit_on_accept=True`. The coverage claim is about the wired rows and
  inherits this unchanged.
* `IFMaterialCard.citation` is `None` in this entry's fixture, so the response may be reported only
  as a **linear tangent limit** and never as "strain-stiffening". Nothing here reports it.
* `SeriesJointConnector` routes exactly two force blocks; `focal_adhesion` receives none, because it
  has no degrees of freedom to receive it with. A three-member series joint is a different object and
  the frozen module refuses one at import.

## 7. Independent oracle or derivation

Five, each of which is a check on **physics** rather than on the reference, and each of which can
disagree with the parity gate:

* **O1 — the series combination is softer than either member.** `1/k = 1/k_a + 1/k_l` gives
  `k < min(k_a, k_l)`, while the parallel combination gives `k > max(k_a, k_l)`. Measured through
  the kernel against the halves the state carries. `WRONG_PARALLEL_STIFFNESS` fails it by a factor
  of `(k_a + k_l)^2/(k_a k_l)`, which for this fixture is 8.17.
* **O2 — the cable's tension is the derivative of its energy.** A centred finite difference of
  `U(s)` against the kernel's own `T`, at strains placed inside each of the three branches, must
  show order 2. `WRONG_ENERGY_BRANCH_OFFSET` fails it; `WRONG_OMIT_RESTIFFENING` **passes it
  exactly**, because dropping a branch leaves the law a perfect gradient of its own energy — which
  is the reference's own recorded warning and the reason O5 exists.
* **O3 — the stall force contains no rate.** Scale `v0` by 1e-3 and by 1e+3 and the load at `v = 0`
  must not move. A ratchet that had degraded into a linear spring or a drag fails this, because a
  spring's blocking force is proportional to its stiffness. `WRONG_LINEAR_LOAD` **passes it**, and
  that is the point of O5 below.
* **O4 — the crossbridge energy is the elastic energy of the strains and nothing else.** `2U/k_xb`
  must equal `sum x_h^2` to round-off, computed from the state without touching the kernel.
* **O5 — the oracle that needs more than one point, and the reason it is here.** `-3610` measured
  that an analytic oracle can score the wrong kernel better than the right one. Both of this entry's
  *conservative* mutants are that shape again, and the ratchet's is sharper:

  > **`WRONG_LINEAR_LOAD` matches the true law EXACTLY at both ends of the load–velocity curve, by
  > construction.** At `v = 0` it *is* the stall force; at `v = v0` it is exactly zero. So it passes
  > O3, it passes the free-velocity anchor, and it passes the vanishing-power check at both ends —
  > and every one of those is an *endpoint* of this curve.

  Measured: the two curves agree at the endpoints to `< 1e-6` of the stall force and disagree by
  **25.0% at `v/v0 = 0.336`**, with 22.7% at the midpoint. So the oracle is evaluated at **three**
  interior points straddling the peak, and
  `test_a_linear_ratchet_matches_the_true_law_at_both_ends_of_the_curve_and_nowhere_between`
  asserts the endpoint agreement as well as the interior gap, so nobody reads the endpoint pass as
  evidence.

  `WRONG_OMIT_RESTIFFENING` is the other shape and it is worse in one respect: it is a **perfect
  gradient**, measured to `< 1e-6` relative by a centred difference of the reference's own
  `omit_restiffening_branch` flag. No gradient probe, no closure check and no energy ledger can
  express a preference at all. Only the parity gate against the true card, or a strain-limit probe
  past `ep`, can — which is why §9's "caught on" column does **not** list O2 for it.

**Control names:** `test_the_series_stiffness_is_softer_than_either_half_through_the_kernel` (O1),
`test_the_cable_tension_is_the_derivative_of_its_energy_on_every_branch` (O2),
`test_the_stall_force_does_not_move_when_the_free_velocity_is_scaled` (O3),
`test_the_crossbridge_energy_is_the_elastic_energy_of_the_strains` (O4),
`test_a_linear_ratchet_matches_the_true_law_at_both_ends_of_the_curve_and_nowhere_between` and
`test_a_cable_kernel_that_dropped_the_restiffening_branch_is_caught` (O5).

## 8. Positive control

`test_the_series_joint_kernel_reproduces_its_frozen_reference`,
`test_the_cable_kernel_reproduces_its_frozen_reference_on_all_three_branches`,
`test_the_ratchet_kernel_reproduces_its_frozen_reference`,
`test_the_crossbridge_kernel_reproduces_its_frozen_reference` — each drives the TRUE kernel and the
frozen NumPy law on identical inputs through `parity.law_parity_report` and requires every field
inside its declared budget.

`test_the_declared_budgets_sit_inside_what_roundoff_explains` requires each declared number to be
below the round-off explanation computed from the state's own conditioning, so a budget cannot be
widened to fit a measurement.

`test_every_batch_member_is_built_through_its_own_frozen_builder` builds each of the four classes
through the vertical's own builder rather than constructing it here, so a change to a builder's
parameterisation is a red test rather than a drift.

## 9. Deliberately failing negative control

Thirteen wrong kernels ship in `law_kernels.py`, compiled by the same `load()` and launched down the
same driver. Each is a plausible hand-translation error and each is killed by a **named** control.

Measured on the warp CPU device, float32, `ORDERED` scatter, `GLOBAL_F32`, on each family's default
state. Every figure is the worst field in that channel, in float32 ULP.

| family | variant | the defect | energy | force | caught on |
|---|---|---|---|---|---|
| **J** | **TRUE** | — | **5.24** / 64 | **68.68 / 59.94** / 128 | pass |
| J | `WRONG_ENERGY_HALF` | `U = k e²` | 8,388,606 | *bit-identical* | **energy alone** |
| J | `WRONG_UNILATERAL_CLAMP` | `e = max(d − d0, 0)`; a compressed clutch goes slack | 2,488,110 | 6,483,010 | both; **bit-identical on an all-stretched state** |
| J | `WRONG_PARALLEL_STIFFNESS` | `k = k_a + k_l`; needs the halves the true kernel is not given | 82,680,000 | 78,590,000 | both, and O1 |
| J | `WRONG_OPEN_CLUTCH_LEAKS` | the gate ignored; an open clutch transmits | 3,131,000 | 5,329,000 | both, and **family J's exact-zero control** |
| **K** | **TRUE** | — | **45.34 / 32.33** / 64 | **8.04 / 4.71** / 128 | pass |
| K | `WRONG_ENERGY_BRANCH_OFFSET` | `T_k b` and `T_p c` dropped from `W` | 3,063,000 | 505,900 | both — **and force is bit-identical on a STATIONARY state** |
| K | `WRONG_OMIT_RESTIFFENING` | the plateau extended forever | 6,655,000 | 7,725,000 | both; **passes O2 exactly** — it is a perfect gradient |
| K | `WRONG_TOE_COMPRESSION` | the toe extended below zero strain; a slack cable pushes | 1,347,000 | 909,800 | both, and the exact-zero control |
| **L** | **TRUE** | — | **0.0** (exact, grades nothing) | **0.40 / 0.40** / 2, power **0.031** / 2 | pass |
| L | `WRONG_LINEAR_LOAD` | `F = F_stall(1 − v/v0)` | 0.0 | 5,245,000, power 5,433,000 | force and power; **exact at both endpoints**, visible only between |
| L | `WRONG_LOAD_PER_SITE` | `n_i` dropped; a bundle tip stalls like one filament | 0.0 | 6,800,000 | force and power; bit-identical where every `n = 1` |
| L | `WRONG_SIGN` | the ratchet pulls the membrane inward | 0.0 | 16,780,000 | **force alone — the power channel is bit-identical** |
| **M** | **TRUE** | — | **0.245** / 1.5 | **0.146 / 0.314** / 1.5 | pass |
| M | `WRONG_ENERGY_HALF` | `U = k_xb Σx²` | 8,388,606 | *bit-identical* | **energy alone** |
| M | `WRONG_ABS_STRAIN` | `\|x\|`; a negatively strained head pulls the wrong way | *bit-identical* | 17,440,000 | **force alone — the energy is even in `x`** |
| M | `WRONG_STRAIN_FROM_SEPARATION` | `k_xb(\|a − b\| − rest)`; the motor made a spring | 398,100,000 | 36,230,000 | both; **it reads a position the true kernel does not** |

**Self-comparison is exactly `0.00` ULP on every field of all four families.**

> ### ⚠ THE §9 TABLE'S FIRST DRAFT WAS WRONG IN ONE ROW AND THE MEASUREMENT CORRECTED IT
>
> `WRONG_ENERGY_BRANCH_OFFSET` was written as "energy alone; every force bit-identical", on the
> argument that the tension is untouched. It is not: **the shipped force on a moving site is a
> discrete gradient of two energies**, so a wrong energy moves it. The bit-identity holds exactly on
> **stationary** sites, where the force comes off the derivative instead, and
> `test_a_cable_kernel_whose_energy_lost_the_branch_offsets_is_caught` asserts it there.
>
> Recorded rather than silently fixed, because a claim written from a derivation and never measured
> is exactly what `ALEPH-PORT-3610` §14.10 says a blind spot attributed to the wrong object costs.
> **Five of thirteen mutants are invisible in some channel, in four different channels**, and the
> two energy-half mutants and `WRONG_ABS_STRAIN` are the clean ones.

### 9a. A Newton-pair control is not shipped for ANY of the four, and this is a measurement

`HANDOFF.md` §F-1 and `ALEPH-PORT-3604` §9a measured `force_b := -force_a` unkillable for a central
linear law. Here it closes as an **argument that covers all four families at once**, and one of them
contradicts a sentence in the frozen source:

* **L and M assemble both sides from one scalar and one unit vector**, `∓F t`. Negation is exact in
  IEEE-754, so the mirror is bit-identical by construction.
* **J and K normalise twice, deliberately, and it buys nothing.** `SeriesJointConnector.evaluate_sites`
  says in as many words that the two blocks are built "from two independent normalisations, not from
  one negation", and that "the consequence is that closure holds to round-off rather than
  bit-exactly, which is the price of the check being a check."

> **Measured: that consequence does not hold.** `np.linalg.norm(x)` and `np.linalg.norm(-x)` are
> bit-identical, because the squares are identical; dividing exact negations by an identical
> positive denominator gives exact negations. So `force_b` **is** the bitwise negation of `force_a`
> at every site, and the closure check the second normalisation was written to preserve is a
> tautology anyway.

`test_the_newton_mirror_is_bitwise_unkillable_for_all_four_families` asserts the bit-identity with
`np.array_equal` over a state that visits every branch, and
`test_the_two_normalisations_of_a_series_joint_are_bitwise_one_negation` asserts the specific
contradiction so it is checkable rather than argued. **No mirror mutant and no force-closure control
ships for any of the four**, and the reason is recorded rather than the control being written to
pass.

### 9b. Two mutants are CONSERVATIVE, and that is a different kind of blindness

`WRONG_OMIT_RESTIFFENING` and (in the reference) `allow_cable_compression` both leave `F = -grad U`
**exactly true**. A finite-difference gradient control, a force-closure check and an energy-ledger
check are all structurally unable to see them: the physics is wrong and the calculus is perfect.
Only a *constitutive* probe — a strain-limit sweep, or the parity gate against the true card — can.
This is the same shape as `-3610`'s Stokes-law trap, arriving through a different door: there the
oracle preferred the wrong kernel, here the oracle cannot express a preference at all.

## 10. Numerical and precision envelope

`ALEPH-DQ-107`: float32 compute, float64 where the accounting needs it. Per site, in every family:
positions, separations, extensions, strains, directions and forces are **float32**; each per-site
energy is **widened to float64 before the product is formed**, so a square is never rounded into
float32 before `Backend.sum`'s float64 accumulator sees it.

**The position-precision decision reaches families J and K and does not reach families L and M.**
J and K consume `k_pair_offset_*` — `-3604`'s kernel, **reused**, so this entry adds **zero** new
places the decision reaches. L and M read no position at all in their true kernels, and
`test_the_position_precision_is_inert_for_the_two_active_families` asserts their force fields are
bitwise identical under all three modes. That is the **second** measured instance of `-3610` §14's
fourth shape for the open question: *a per-law rule must be able to say "this law does not care"*,
and it now covers six more registry rows.

**`WRONG_STRAIN_FROM_SEPARATION` is the exception and it is instructive.** It is the only kernel in
family M that reads a position, so it is the only one `PositionPrecision` moves — a defect that
*acquires* a sensitivity the law does not have. Recorded because a mode sweep that found family M
sensitive would be finding the mutant, not the law.

Measured figures are in §13.

## 11. Production-backend residency and transfer

Every array is allocated through `Backend.zeros` and every reduction through `Backend.sum`; no raw
host `ndarray` is handed to a launch, per commit `636b0c8`'s retraction. Scatter is `ORDERED`.

`-3609` §14.5's hazard applies and is controlled: **`warp.array.numpy()` on warp's CPU device returns
a zero-copy view.** `-3610` M12 then found that a control written on a *narrowing* pull cannot see
the aliasing, because the dtype conversion already copies. So
`test_a_pulled_field_owns_its_memory_and_does_not_view_the_device` is written on a **float64** field
whose device dtype already matches the host dtype — the per-site energy — and checks `owndata` and
`base` directly rather than comparing values.

Nothing in this entry is resident across a step; residency is `-3609`'s layer and is untouched.

## 12. Comments and docstrings to discard

Not applicable: nothing was read from `/Users/sw1/ffn_cellsim`. Every docstring in the new code is
written for the Aleph law it describes.

## 13. Acceptance

`tests/runtime` **550 passed / 3 skipped / 1 failed** (was 505/3/1); `tests/runtime` +
`tests/vertical` **1,927 / 3 / 1** (was 1,882/3/1). The one failure is `-3608`'s guard and is
**foreign to this entry's work** — §14.7. `aleph/vertical/**`, `aleph/viz/**` and `aleph/state/**`
are byte-for-byte untouched, verified by `git status`.

### 13.1 The measured anchors

| | |
|---|---|
| self-comparison, every field of all four families | **0.00 ULP** |
| `f^(b) == −f^(a)` **bitwise**, all four families | **True** — §9a |
| open clutch, per site, `force_a`, `force_b`, `energy` | **exactly `0.0`** on 4 of 12 joints |
| all-slack cable, whole population, budget declared `0.0` | **0.00 ULP**, every field |
| head at exactly zero strain | **exactly `0.0`** force and energy |
| series stiffness vs its halves | `k < min(k_a, k_l)` at every engaged joint; parallel impostor **8.17× / 14.08×** too stiff |
| `dU/ds = T(e)` inside each of the three branches | agrees to `< 1e-8` relative |
| stall force under `v0 × 1e-3` and `× 1e+3` | **bit-identical** — it contains no rate |
| `2U/k_xb` vs `Σ x²` | agrees to `< 1e-6` relative |
| linear-ratchet impostor vs truth, at the two endpoints | `< 1e-6` of the stall force — **it passes every endpoint oracle** |
| the same, at `v/v0 = 0.336` | **25.0%** of the stall force |
| `WRONG_OMIT_RESTIFFENING` as a gradient | perfect to `< 1e-6` — **no gradient probe can see it** |
| `PositionPrecision` across all three modes, families L and M | **bitwise identical** |
| `PositionPrecision` on the family-M *impostor*, 50 µm from the origin | **differs** — it reads a position the law does not |

### 13.2 Position modes on families J and K — and a third category for the open question

`ALEPH-PORT-3605` measured `POSITIONS_F64` **worse** than `LOCAL_F32`; `ALEPH-PORT-3606` measured it
48×/191× **better**; both agreed the mechanism is whether the cancellation is a *subtraction* or a
*reduction*. Family K splits a case they did not distinguish:

| | `GLOBAL_F32` | `LOCAL_F32` | `POSITIONS_F64` |
|---|---|---|---|
| J energy | 5.24 | 4.67 | **0.61** |
| J force (a) | 68.68 | 64.69 | **1.05** |
| K energy (start) | 45.34 | 12.05 | **2.25** |
| K force (a) | 8.04 | 10.91 | **8.42** |

> **Family K's ENERGY improves 20× under float64 positions and its FORCE does not improve at all.**
> The energy is a subtraction; the force is a *difference quotient of two energies*, and its dominant
> error is downstream of the subtraction. So a law can be sensitive to the representation in one
> output field and inert in another — a **third** category the open proposal's per-law option (C)
> does not currently have a shape for, alongside `-3610`'s "this law does not care".
> `test_the_position_precision_is_live_and_not_monotone_for_the_two_configuration_families`.

### 13.3 The mutation study: 25 planted, 25 killed

Thirteen are the shipped wrong kernels of §9, each killed by a named control. Twelve more were
planted in the shipped modules, confirmed to redden a named test, and reverted —
`PYTHONDONTWRITEBYTECODE=1`, from a bytecode tree verified empty under `aleph/` and `tests/`, with
every revert from a file copy rather than from `git checkout --` (`ALEPH-PORT-3610`'s recorded near-
miss).

| # | planted defect | reddens |
|---|---|---|
| M14 | `refuse_nonpositive_series_geometry` removed from `_series_joint_launch` | `test_the_host_refusals_run_on_the_launch_path` |
| M15 | `refuse_non_unit_axes` removed from `_ratchet_launch` | the same |
| M16 | `np.array` → `np.asarray` in `_per_site` — the `-3609` §14.5 view hazard | `test_a_pulled_field_owns_its_memory_and_does_not_view_the_device` |
| M17 | every clutch engaged, emptying the open-joint mask | `test_an_open_clutch_carries_exactly_zero_load_through_the_kernel` |
| M18 | `DEFAULT_CROSSBRIDGE_ULP_BUDGET` force widened 1.5 → 4.0 | `test_the_declared_budgets_sit_inside_what_roundoff_explains` |
| M19 | `_amp_crossbridge_force` returns the constant `1.0` — the *pre-repair* derivation | the same |
| M20 | the re-stiffening band deleted from `strain_stiffening_cable_state`'s placement | `test_the_cable_state_visits_every_branch_of_the_card` |
| M21 | the zero-strain head set to `1e-18` instead of `0.0` | `test_a_head_at_exactly_zero_strain_carries_exactly_zero_force_and_energy` |
| M22 | `TensileLinker` (a spring) added to `ACTIVE_AND_GATED_BATCH` | `test_the_three_connector_batches_are_disjoint` |
| M23 | the true family-J launch fed `half_actin` where the series stiffness goes | `test_the_series_joint_kernel_reproduces_its_frozen_reference` |
| M24 | the family-K launch fed `plateau_end_strain` where `knee_strain` goes | `test_the_cable_kernel_reproduces_its_frozen_reference_on_all_three_branches` |
| M25 | `k_series_joint_true`'s `force_b` written from `force_a`'s own direction | `test_the_newton_mirror_is_bitwise_unkillable_for_all_four_families` |

### 13.4 ⚠ Two of these found real defects in this lane's own work, and one was in a DERIVATION

**M18/M19 did not need to be planted — the control fired on a CORRECT kernel first.** The first
draft's family-L and family-M amplifications divided a per-site error bound by the largest *per-site*
force, following `ALEPH-PORT-3606`'s note that "a vertex receiving several sites has a larger force,
so this can only widen the bound". **That is true only when the sites add constructively.** These two
families scatter along independent unit vectors, so a vertex can receive several sites that largely
cancel while their errors add. Measured before the repair: a 48-site ratchet at **1.259 ULP against a
1.197 explanation**, and a 10-head crossbridge at **1.624 against 1.000** — both the true kernel,
both over-running a bound that under-explained them.

> **A round-off derivation inherited from another lane can be wrong for a law with a different
> assembly, and the symptom is a CORRECT kernel failing.** The repair was to the derivation
> (`_assembled_amplification`), never to the budget, which is the direction
> `test_the_declared_budgets_sit_inside_what_roundoff_explains` exists to force. M18 and M19 replant
> both halves — the widened budget and the pre-repair derivation — and both are killed.

**M20 survived its first planting, and the reason is a new shape.** Deleting the re-stiffening band
from the cable state's *placement* left `restiffening_sites > 0` anyway, because the step pushes some
plateau sites past `ep` on its own. **The coverage control was satisfied by an accident of the step
rather than by the placement it was checking.** Fixed by reporting the branch counts at the **start
configuration alone** as well as over both ends, and asserting on both. Replanted, killed. This is
`ALEPH-PORT-3606` O2's finding in a new place: there a control asserted over an empty set; here it
asserted over a set filled by something other than the thing under test.

### 13.5 A ULP budget is a property of the state's ASSEMBLY, not of the kernel

Family M's law is `f = k_xb x t` — a product with no subtraction, no reduction and no transcendental
— so its per-site conditioning is exactly **1** and its entire budget is the scatter. A 10-head,
5-vertex assembly puts the **true** kernel at **1.62 ULP against the declared 1.5**, while staying
comfortably inside its own round-off explanation of **3.65**. Neither number is wrong and the budget
was not widened to hide it:
`test_the_crossbridge_budget_is_at_its_own_roundoff_floor_and_a_denser_assembly_exceeds_it` asserts
both halves, so the declared 1.5 is understood as a statement about `crossbridge_state()`'s assembly.

### 13.6 The census

**30 of 30 wired registry rows are on kernels. 0 remain.** `81` kernels in `law_kernels.py` (was 64).
Rows and kernels are asserted separately against `wiring.WIRED` by
`test_the_batch_covers_the_last_eleven_wired_registry_rows` and
`test_the_ported_row_census_is_now_complete_against_the_wiring_registry`.

## 14. Honest limits

1. **`UNVERIFIED` on a CUDA device.** Every figure is warp's CPU device. `-3609` measured that a
   CUDA device changes the float32 noise floor enough to move a descent window by a binary order,
   and **family L computes a logarithm**, which is precisely the class of operation CUDA implements
   differently from libm. So family L's budget is the one in this entry with the least claim to
   transfer, and it is named here rather than discovered later. A device run is optional for this
   lane and was not taken.
2. **Two owners' arithmetic stays on the host and this entry does not claim otherwise.** §4.5: the
   clutch's series reduction and the motor's kinetics.
3. **Batch G's forces are functions of state this port does not compute.** A parity gate cannot see
   a defect in the state it grades on (`-3605` §13a) — and for families L and M the state *is* the
   physics. A growth velocity or a head strain that was advanced wrongly produces two agreeing
   evaluations of a wrong number. The mitigations shipped are O3 and O4, which are checks on the
   law's own structure, and the exact anchors of I1.
4. **A parity gate still cannot see a defect in the reference it grades against** (`-3606` M8) **or
   when its asserted-over set is empty** (`-3606` O2). Both inherited; every `==` control here
   carries a non-emptiness assertion on its mask.
5. **The energy channel of family L is an exact zero and grades nothing.** As in `-3610` §14.8: it
   is true of any kernel that writes zeros there, including every mutant. The force channel and the
   active power carry that family's port.
6. **`ports/ledger/INDEX.md` is stale and red before this entry existed**, naming `-3606` through
   `-3610`. `-3611` adds a sixth *name*, not another red test. **Not regenerated**, on instruction;
   reported here and in `docs/ACTIVE_SESSIONS.md`.
7. ### ⚠⚠ THE `-3608` GUARD IS RED ON **TWO** CASES, NOT ONE, AND THE SECOND WAS INVISIBLE

   `tests/runtime/test_reference_input_guard.py::test_every_registered_case_reference_is_sensitive_to_position_precision`
   is red before this entry, on `-3610`'s `isotropic_drag_law_case` at **0.443 ULP**, and that is a
   PI decision (`HANDOFF.md` §C-0, `-3610` §14.11). **This lane did not touch that file and takes
   none of the three repairs `-3610` named.**

   Driving the guard's own procedure over every registered case by hand, which is not what the test
   does, produced something the test cannot report:

   | case | worst field, ULP | verdict |
   |---|---|---|
   | `isotropic_drag` | **0.443** | RED — known, PI decision |
   | **`dense_exterior_resistance`** | **0.262** | **RED — and `-3610` §14.11 states it passes** |
   | `anisotropic_drag`, `brownian_ratchet`, `crossbridge` | — | pass **by refusal** |
   | the other 12, including `series_joint` (68.26) and `strain_stiffening_cable` (45.73) | 2.07 … 519.48 | pass on measured sensitivity |

   > **The assertion is inside the loop, so the test stops at the first failing case and every case
   > after it is never evaluated.** `dense_exterior_resistance` sorts after `isotropic_drag` in the
   > registry, so its redness has never been printed by anything. `-3610` §14.11's sentence *"the
   > dense case passes on measured sensitivity"* is **wrong**, and it was written from a run that
   > could not have measured it.

   **And the two share a mechanism, which makes this a class rather than an anomaly.**
   `f = −ζv` and `f = −RV` are both **linear in the rounded input**, so rounding it moves the output
   by exactly the input's own rounding — under 0.5 ULP, always, at any state. Family M's law
   `f = k_xb x t` is a third of the same kind and would be a fourth red case if it did not refuse
   the narrowed configuration for an unrelated reason (§14.8).

   > The guard's threshold assumes every law amplifies. **A law with unity conditioning cannot clear
   > a bar set above its own conditioning**, and there are now at least three such laws in the
   > registry carrying seven wired rows. That is the shape of the decision, and it is larger than the
   > one `-3610` raised.

   **Reported, not resolved.** Two of `-3610`'s three candidate repairs edit `-3608`, which the PI
   ratified; the third routes around the guard. `CLAUDE.md` §2 rules 4 and 5. **The one thing this
   entry adds is that a repair which special-cases a single law would leave the other two red.**
8. **Two of this entry's four cases pass that guard by REFUSAL rather than by a number, and the
   reason is structural.** Rounding a state's `(N, 3)` float arrays to float32 moves a unit vector
   off unit length by ~3.9e-08, and both batch-G references hold a unit-vector constraint at 1e-12 —
   `BrownianRatchetConnector.set_axes` and `ActinBindingSites.__post_init__` — so the frozen
   reference raises. The guard explicitly accepts a refusal as the strongest form of sensitivity and
   `-3610`'s anisotropic case already passes it the same way. Recorded because a *set* in which too
   many cases refuse leaves that control asserting nothing about any number, and the guard's own
   `measured >= 10` floor is what watches for it. The two batch-F cases are numeric.
9. **The cable is one row and it is the only row in the port whose constitutive card has three
   branches.** Its ULP figures are therefore branch-dependent in a way no other family's are, and a
   single worst-case number for it is less informative than the per-branch table in §13. Both are
   reported.

10. ### ⚠⚠ AN EXACT ANCHOR THIS PORT LOSES, AND TWO OF ITS OWN RULES POINT OPPOSITE WAYS

    At `v = v0` a barbed end carries **exactly** `-0.0`: the argument of the logarithm is exactly
    `1.0` and `log(1) = 0`. The float64 reference produces that number exactly. **The float32 kernel
    does not**, and the cause is one ULP of *association*:

    | | argument | `log` | load |
    |---|---|---|---|
    | float64 reference | `1.0` | `0.0` | **`-0.0` exactly** |
    | float32, the reference's association `((1−ε)·v)/v0` | **`0.99999994`** | `-5.96e-08` | `+9.45e-08 pN` |
    | float32, reassociated `(1−ε)·(v/v0)` | `1.0` | `0.0` | **`-0.0` exactly** |

    against a stall force of 3.361 pN — a relative `2.8e-08`, which is float32's own floor and not a
    physics error. But it is an exact anchor, and this port's exact anchors have been its sharpest
    evidence.

    > **§3 requires the kernel to match the frozen reference's association order, so a reported ULP
    > is the representation and not a gratuitously different arithmetic. §5 I1 requires every exact
    > anchor to be preserved and asserted with `==`. Here the two rules point opposite ways, and the
    > distance between them is exactly one float32 ULP of the argument.**

    This lane keeps §3, reports the anchor as **lost**, measures the residual and measures the
    reassociation that would restore it — **and takes neither repair**, per `CLAUDE.md` §2 rule 5.
    `test_the_free_velocity_exact_anchor_is_LOST_by_one_ulp_of_association` asserts all four facts,
    including that reassociating still restores it, so if warp's `log` or the association ever
    changes the control goes red rather than the finding going stale.

    **The decision is not this lane's**, and it is not only about this law: any kernel whose exact
    anchor sits at a value the reference's association reaches through an intermediate rounding has
    the same choice. It is the first time in G1–G5c that a *porting rule* and an *evidence rule*
    have conflicted rather than a rule and a measurement.

11. **`WRONG_ENERGY_BRANCH_OFFSET`'s channel claim was wrong in the first draft of §9 and is
    corrected in place** — see §9's own block. It was written from a derivation ("the tension is
    untouched, so the forces are untouched") and never measured; the shipped force on a *moving*
    site is a discrete gradient of two energies. The claim is now what the measurement says and the
    bit-identity is asserted where it actually holds.

12. **Family M's parity gate is the tightest in the port and it is graded on state this port does
    not compute.** §14.3 already says the second half; the first half matters too, because a budget
    of 1.5 ULP has no room to absorb a change in the assembly (§13.5) and none at all to absorb a
    device. Of the four families this is the one whose numbers would need re-measuring first on
    CUDA, with family L second because of `wp.log`.

13. **This entry ports no kinetics and no owner arithmetic, and the two absences are asymmetric.**
    Family J's missing half (the clutch's series reduction) is *pure* — it takes no step and could be
    put on a kernel by a later lane with no transaction question. Family M's missing half (attachment,
    detachment, the stroke) is *stateful and transactional*, and putting it on a device is a
    different kind of task from anything G1–G5c has done. A plan that treats "the remaining NMII
    work" as one item would be treating those two as the same size.
