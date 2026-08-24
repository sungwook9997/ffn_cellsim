# ALEPH-PORT-3624 — the adhesion rows: one series joint, dispatched whole

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3624` |
| Lane | `46143f30` Lane W (the whole cell, for ablation) |
| Status | `PROPOSED` |
| Written | `2026-08-04` |
| Port class | `RE-DERIVED` — nothing is copied; this wires Aleph's own `focal_adhesion` into Aleph's own world |
| Aleph target | `aleph/scenarios/whole_cell_couplings.py`, `aleph/scenarios/whole_cell.py` |
| Depends on | `ALEPH-PORT-3621` (the whole cell), `ALEPH-PORT-3623` (the couplings), `ALEPH-PORT-3303` (the focal adhesion owner) |
| Exists because | `focal_adhesion` was the last compartment in the whole cell that no connector reached. It is the only owner with **no geometry at all**, so every routing rule the other twenty-eight couplings use is inapplicable to it, and the four declared adhesion rows sat in the registry with no instance. |
| What this entry does **not** claim | That any traction number here is physics. Every position is where the fixture's constructor put it, and the resulting pre-load is large (§6). Evidence class `UNVERIFIED` throughout. |

---

## 1. The thing that makes this different from every other coupling

`focal_adhesion` owns no node array, no position array and no material point. That is enforced by the
owner rather than being an accident of the fixture — `assert_owns_no_material_points` is a method on
`FocalAdhesionClutchGraph`. Its two ends are **names**: an owner name and a site id on each side.

So the load path is

    actin owner  --(k_actin)--  plaque  --(k_ligand)--  ecm

and the plaque in the middle has no degrees of freedom. It supplies a compliance and an engagement
gate, and it receives nothing, because there is nothing for a force to be scattered into.

`SeriesJointConnector.endpoints` therefore reports `(actin_owner, "ecm")` and **not** the registered
endpoint pair of any row it serves. Every row names `focal_adhesion` as one endpoint.

## 2. The dispatch unit is the series group, not the registry row

`fa_actin_anchor` and `integrin_collagen_clutch` are two registry rows describing **one** series
joint. The owner refuses half a group at connector construction:

> composite group `alpha2beta1_collagen_series` was dispatched as `['fa_actin_anchor']`, missing
> `['integrin_collagen_clutch']`. Its members are one series joint: they share one load and their
> extensions add […] Evaluated apart, each half carries the load its own extension implies, the two
> do not balance, and the plaque — which has no degrees of freedom — becomes a hidden force source.
> The joint also reports the parallel stiffness, which is stiffer than either member, and stiffer
> joints read as higher traction.

That refusal was **received**, not read: the first construction attempt in this lane dispatched
`fa_actin_anchor` alone and was refused. The guard is the reason this entry exists in the shape it
does.

### 2.1 The refusal as a number

Driven through `dispatch_rows_independently`, on this world's clutches
(`k_actin = 1200`, `k_ligand = 400`, series `300` pN/µm):

| variant | max &#124;F_actin&#124; | Newton &#124;ΣF_a + ΣF_b&#124; | direction vs true |
|---|---:|---:|---:|
| **true (whole group)** | 445.87 pN | **0.00e+00** | +1.0000 |
| `dispatch_rows_independently` | 1783.46 pN (**4.00×**) | **1.90e+03 pN** | +1.0000 |
| `transmit_through_an_open_clutch` | 445.87 pN | 0.00e+00 | +1.0000 |
| `flip_joint_sign` | 445.87 pN | 0.00e+00 | **−1.0000** |

Read the second row: splitting the group makes **1.9 nN of force appear from a body with no degrees
of freedom.** That is the docstring's "hidden force source" measured rather than quoted.

Read the fourth row: a pure sign error **closes exactly**. A check that only tested `Fa + Fb == 0`
would pass a connector doing the opposite of physics. Direction is what catches it, and the cosine is
exactly −1.

Controls: `test_splitting_the_series_group_breaks_the_newton_pair_by_1900_pn`,
`test_flipping_the_joint_sign_reverses_every_force_and_still_closes`,
`test_the_series_stiffness_is_below_either_half`.

## 3. A declared blind spot, measured and then removed

`transmit_through_an_open_clutch` is **bit-identical to the true kernel on this world**, and that is
a property of the state rather than a defect in the gate: every clutch in `build_whole_cell` is
`BOUND/BOUND`, so no clutch is open for the gate to close on.

This is the distinction `law_kernels.DECLARED_BLIND_VARIANTS` was built to draw, and the same rule
applies: **measure why it is blind before declaring it.** Measured, the fix turned out to be a better
state rather than a declaration. On a mixed population the mutant is caught cleanly:

| clutch | engagement gate | gate defeated |
|---|---:|---:|
| bound | 165.0 pN | 165.0 pN |
| **released** | **0.0 pN** | **165.0 pN** |

A released clutch is `actin_side=None`, **not** `actin_side_state=UNBOUND` with the reference kept —
the owner refuses that by name: *"an unbound side holding a reference is an attachment that has
quietly not let go."* Which is why a whole-cell clutch cannot simply be flipped open, and why the
control builds its own two-clutch population.

Controls: `test_the_engagement_gate_is_only_visible_on_a_population_that_has_a_released_clutch`, and
`test_every_clutch_in_the_whole_cell_is_engaged_so_the_gate_control_needs_its_own_population`, which
pins the premise so a future fixture change turns it red instead of quietly invalidating the reason.

## 4. Three adapters, each demanded by a guard that refused first

Nothing here was designed ahead of the refusals. Each adapter exists because a guard rejected the
simpler thing, in this order:

| # | Refusal received | Adapter |
|---|---|---|
| 1 | `FocalAdhesionRoleError: clutch 'fa1' does not publish role 'fa_actin_anchor'` | publish both roles on every clutch, through the owner's own gate |
| 2 | `FocalAdhesionSeriesError: composite group … missing ['integrin_collagen_clutch']` | dispatch by group; `_rows_for` reads the group from `composite_groups()` rather than tabulating it |
| 3 | `AdhesionCouplingError: the actin anchor block has 17 rows for 1 clutch handle(s)` | `_ClutchSiteOwner`, so the world's **out-of-step** probes see the clutch block and not the whole node array |
| 4 | `WorldCompositionError: returns a (1, 3) block for 'filopodium', whose force array is (17, 3) … a routing adapter is needed` | `_ClutchRoutedJoint` + `_material_point_scatter` |
| 5 | `CoverageViolation: … was reached but the compute substrate never moved (witness delta 0)` | the scatter moved onto `ctx.backend.scatter_add` |

Refusal 4's message is worth quoting in full because it carries its own history:

> This connector does not deliver to its declared endpoints — a routing adapter is needed, and
> **dropping the block instead would remove the coupling while the world still relaxed**.

and `world.py`'s comment above it names the incident: *a spreading world once lost every adhesion to
exactly that silent drop, and relaxed to a plausible equilibrium with no adhesion in it.* Refusing
rather than skipping is what turned a silent physics loss into a five-minute fix.

## 5. A clutch pulls at a material point, not at a node

Every other coupling in `whole_cell_couplings.py` slices node arrays: `positions[index]` in, scatter
back to `index` out. That is **wrong here**. An adhesion site is at a *material coordinate* — a rest
arclength along a bundle — which in general falls between two nodes. `endpoint_sink(site_id)`
resolves it to a `MaterialPoint` carrying `(node_a, node_b, weight)`.

Snapping to the nearer node would move the load path by up to half a segment and would do it
**silently**: every energy stays finite, closure still holds, and the traction is simply applied
somewhere the clutch is not.

So the block reads `material_point.position` and splits each force `(1−w)·f` at `node_a`, `w·f` at
`node_b`.

### 5.1 One rule read twice, and the reading is checked

Refusal 5 forced the scatter onto `ctx.backend.scatter_add`, which meant expressing
`MaterialPointSink.scatter`'s two lines as an index/value pair. **A reimplementation nobody compares
is a second rule.** `test_the_backend_scatter_matches_the_owners_own_material_point_sink` runs both
paths on the same forces and requires the arrays to be *equal*, not close — the arithmetic is
identical when the rule is. It also asserts that at least two nodes are loaded, which is what would
fail if a site were ever snapped to a node.

## 6. What is now wired, and the one number a reader should be suspicious of

| | before | after |
|---|---:|---:|
| couplings in the world | 14 | **17** |
| adhesion registry rows dispatched | 0 / 4 | **4 / 4** |
| compartments no connector reaches | 2 | **1** (`extracellular_medium`) |
| `tests/scenarios/test_whole_cell.py` | 26 passed | **34 passed** |

Force transfer, measured per connector at the fixture's own configuration:

| connector | clutches | k_series | max &#124;F&#124; | Newton residual |
|---|---:|---:|---:|---:|
| `filopodium_nascent_fa` | 1 | 300 pN/µm | 170.05 pN | 0.00e+00 |
| `lamellipodium_nascent_fa` | 1 | 300 pN/µm | 169.98 pN | 0.00e+00 |
| `alpha2beta1_collagen_series` | 2 | 300 pN/µm | 445.87 pN | 0.00e+00 |

**The suspicious number is the force, and it is the fixture's fault rather than the connector's.**
The clutch rest gap is `0.05 µm` and the fixture leaves the actin and ECM sites `0.62–1.54 µm` apart,
so the joint starts stretched by ~12–30× its rest gap and the world carries a large pre-load at
`t = 0`. Per clutch at occupancy 4 that is ~56 pN per bond, which is at the upper end of a single
integrin bond's rupture force — so a kinetics run would break these clutches rather than hold them.

That is a **fixture-geometry question, not a wiring one**, and it is left open rather than tuned
away: choosing where a stress fibre meets the collagen is a scenario decision, and picking a number
that makes the pre-load look reasonable would be choosing the answer. Flagged for the PI in §9.

## 7. The census had to change, and the change is a real distinction

`WholeCell.uncoupled` was computed from `slot.endpoint_names()` — where the **force blocks** land.
For `focal_adhesion` that is permanently the wrong question: it can never receive a block, so it
would have been reported uncoupled while four connectors ran their whole load through its compliance
and its gate. **An ablation reading that census would conclude removing it costs nothing, which is
the exact opposite of the truth.**

So the census now also reads the registry rows a connector serves (`_rows_endpoints`), and an owner
named by a live row is in the load path whether or not it receives anything. The two questions —
*in the load path* and *receives a force block* — are different, and only the first belongs in an
ablation census.

## 8. Scenario state this entry changed

`_focal_adhesion` carried two clutches, both naming `sf_arc`. It now carries four, adding one on
`lamellipodium` and one on `filopodium`. The registry declares `lamellipodium_nascent_fa` and
`filopodium_nascent_fa` as connectors in their own right — a nascent adhesion under a protrusion is a
different load path from the matured traction spine, and it is the one a spreading cell forms first.
With the population naming only `sf_arc`, those two rows were declared and unpopulated: nothing
failed, and the protrusions simply transmitted no traction to the substrate.

This is a change to the *fixture's* clutch population, not to any compartment's physics.

## 9. Open for the PI

1. **The pre-load in §6.** The fixture's site placement gives ~56 pN per integrin bond at `t = 0`.
   Either the fixture geometry should place the actin and ECM sites near the clutch rest gap, or the
   whole cell should be relaxed before any adhesion measurement is read from it. Both are scenario
   decisions and neither is this lane's to take.
2. **`extracellular_medium` is the last uncoupled compartment.** It negotiates through
   `publish_endpoint` rather than the site protocols every other owner uses; wiring it is the
   remaining item.

---

## 10. Source repository identity, source path and symbol

**None, and this is not an omission.** Port class is `RE-DERIVED`: `/Users/sw1/ffn_cellsim` was not
opened for this entry and no source path or symbol was consulted. The clutch graph, the series joint
and the registry rows are all Aleph's own (`ALEPH-PORT-3303`); this entry instantiates them in
Aleph's own world. There is nothing to compare against a provider file, so a *source commit* would
be a citation to something that played no part.

## 11. Units, domains, singular cases, invariants

Units: length µm, force pN, stiffness pN/µm, energy pN·µm. Evidence class `UNVERIFIED`.

**Domain.** A clutch is defined for separation `d > 0`. `d = 0` is not a small-`d` limit but a point
where the joint axis genuinely does not exist — `(x_l − x_a)/|x_l − x_a|` is `0/0` — and
`SeriesJointConnector` refuses a coincident pair by name rather than regularising the denominator.
That refusal is correct and is why `_ClutchSiteBlock` never invents a position: a regularised
denominator would return a finite force in a direction chosen by round-off.

**Invariants held, and each is measured in §2.1 or by a named control:**

| invariant | how it is held | tolerance |
|---|---|---|
| `ΣF_actin + ΣF_ligand = 0` | two independent normalisations, **not** `F_b := −F_a` | **exactly 0** |
| `k_series ≤ min(k_actin, k_ligand)` | series law, asserted against the reciprocal identity | `1e-9` relative |
| open clutch transmits nothing | engagement gate | **exactly 0** |
| a site's load reaches both bracketing nodes | material-point interpolation, weights `(1−w, w)` | exact equality between the two scatter paths |

`AdjointPair` forbids `force_b := −force_a` precisely so that force closure is a check rather than a
tautology. The consequence is that closure would in general hold only to round-off; on this world it
holds **exactly**, which is stronger than required and is asserted as `== 0.0` rather than
`approx(0)` so a drift into round-off would be visible instead of absorbed.

## 12. Positive control, negative control

**Positive control** — `test_the_backend_scatter_matches_the_owners_own_material_point_sink`. The
backend scatter and the owner's own `MaterialPointSink` are run on the same forces and the resulting
arrays must be **equal**, not close. It also asserts that ≥2 nodes are loaded, which is what fails if
a material point is ever snapped to a node.

**Negative controls** — three deliberate breaks shipped on the connector, each driven and each
required to change the answer:

| control | break | required outcome |
|---|---|---|
| `test_splitting_the_series_group_breaks_the_newton_pair_by_1900_pn` | `dispatch_rows_independently` | Newton residual `> 1.0e3` pN and actin block `> 3×` |
| `test_flipping_the_joint_sign_reverses_every_force_and_still_closes` | `flip_joint_sign` | cosine **exactly** −1, closure still exact |
| `test_the_engagement_gate_is_only_visible_on_a_population_that_has_a_released_clutch` | `transmit_through_an_open_clutch` | released clutch `0.0 → 165.0` pN |

The third carries its own premise-check
(`test_every_clutch_in_the_whole_cell_is_engaged_so_the_gate_control_needs_its_own_population`), so a
fixture that later releases a clutch turns the reason red rather than leaving a hand-built population
in place for no stated cause. **A negative control that cannot fail is not a control**, and §3 is the
record of finding one that could not and fixing the state rather than declaring the blindness.

## 13. Numerical and precision envelope

Host `float64` throughout; no `float32` path and no ULP budget, because nothing here is a law kernel
graded against a reference implementation. The two quantities that are asserted **bit-exactly** are
named above: the Newton residual (`== 0.0`) and the agreement between the two scatter paths
(`np.array_equal`). Both are exact because the arithmetic on each side is the same expression, not
because a tolerance was chosen generously — a tolerance here would hide precisely the reimplementation
drift §5.1 exists to catch.

The reciprocal series identity is compared with `np.allclose` (default relative `1e-5`), since
`1/(1/1200 + 1/400)` is not exactly representable.

## 14. Production-backend residency and transfer

Positions and forces live on the **host** as NumPy arrays; the owners are host objects and their
material-point sinks write host memory. The one operation that crosses to the backend is the scatter,
`ctx.backend.scatter_add(owner.forces, indices, values)`, and it crosses because the coverage gate
requires the compute substrate to move (§4, refusal 5) — a step whose only adhesion work was a host
`+=` advanced `ctx.witness_count` by zero and the pipeline refused the whole step.

The backend is read from the **live** `StepContext` inside `resolve` rather than captured at build
time: a world may step on a different device than the one it was assembled with, and a captured
backend would scatter into host memory while the owner's array lived on the device. No device
residency is claimed for this coupling beyond that single scatter; there is no device-resident clutch
state and no GPU kernel in this entry.

## 15. Comments and docstrings discarded

None to discard — nothing was copied. All prose in the touched modules is new and written for Aleph's
vocabulary. No provider identifier, path or comment appears in `whole_cell_couplings.py` or
`whole_cell.py`.

## 16. Acceptance, reviewer, rollback

**Status `PROPOSED`.** Reviewer: the PI. Acceptance evidence to date is
`tests/scenarios/test_whole_cell.py` — **34 passed**, up from 26, the eight new ones being the
controls in §12 and §5.1.

Rollback is one line: remove the `build_adhesion_couplings(owners)` call from `build_whole_cell`. The
four adhesion rows return to `withdrawn` with their reasons, the census reports `focal_adhesion`
uncoupled again, and every other coupling is untouched. The clutch population added in §8 would also
be reverted to two, since the two protrusion rows would have nothing to dispatch.
