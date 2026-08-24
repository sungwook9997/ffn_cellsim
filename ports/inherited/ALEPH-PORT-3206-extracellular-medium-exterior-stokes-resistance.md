# ALEPH-PORT-3206 — extracellular medium: an exterior Stokes resistance on the live cell surface

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3206` |
| Lane | `L32 owner modules — extracellular_medium` |
| Status | `PROPOSED` |
| Written | `2026-07-30`, **before the code**, per PLAN §0.2.5 |
| Port class | `RE-DERIVED` — nothing was ported, nothing was read |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

Aleph target file: `aleph/vertical/medium.py`.
Controls: `tests/vertical/test_medium_controls.py`.

```python
from aleph.vertical.medium import (
    COMPONENT_NAME,
    MEDIUM_ENDPOINT_ROLES,
    OWNED_STATE_KEYS,
    PUBLISHABLE_ENDPOINT_ROLES,
    REGULARISATION_RATIO,
    ExteriorDragChannel,
    ExteriorStokesMedium,
    MediumCard,
    MediumCouplingError,
    MediumEndpointRole,
    MediumMobilityError,
    MediumOwnershipError,
    MediumQuadratureError,
    MediumRoleError,
    MediumStaleHandleError,
    SurfaceQuadrature,
    SurfaceTractionStencil,
    assert_medium_is_not_darcy,
    assert_no_duplicate_medium_owner,
    assert_state_keys_disjoint_from,
    default_medium_card,
    exterior_mobility_um_per_pn_s,
    fibonacci_sphere_quadrature,
    owned_state_keys,
    regularisation_um,
    stokes_rotation_resistance_pn_um_s,
    stokes_translation_resistance_pn_s_per_um,
)
```

**This is a state owner, and it wires zero connectors.** It publishes the medium-side endpoint of the
one declared connector that names `extracellular_medium` — `membrane_medium_traction` — and
implements the connector itself nowhere. A `Binding` row in `aleph/vertical/wiring.py` belongs to the
commit that wires a connector, never to this one, and this lane added none.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — **recorded from Aleph's own prior ledger entries, not verified by this lane** |
| Source path | none identified — **no exterior Stokes solve is known to this lane to exist in the reference at all** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened for this module. |
| Working tree == commit? | not applicable — nothing was read, so there is no revision to have mis-cited |

The commit hash is carried forward from `ALEPH-PORT-2801` and `-3201..-3205` so that a later auditor
has one identifier to work from. Recording it is **not** a claim that anything was consulted, and
this lane did not verify that it resolves.

Unlike the sibling entries, this one cannot name a plausible source path. The registered connector
contract for `membrane_medium_traction` says in its own words that *"no such solve is declared
anywhere in this lane"*, and this lane found no reason to believe one exists on the other side
either. That is stated as ignorance, not as a finding about the reference.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** The module is written from two of Aleph's own documents: the
registered `extracellular_medium` contract in `aleph/state/census_environment_surface.py:184`, and the
registered `membrane_medium_traction` connector contract in
`aleph/state/connectors_surface_traction.py:419`.

The mathematics is textbook and re-derivable in an afternoon: the Stokes Green's function, a
single-layer collocation on a surface quadrature, a blob regularisation whose width is *fixed by the
quadrature rather than chosen*, and the two closed-form sphere limits. Per `ports/TEMPLATE.md` §3
that is a reason to write clean-room and cite prior art, not a reason to port. The one part that
would have been worth inheriting — an empirically discovered failure mode — is exactly the part no
source is known to have.

## 4. Physical or mathematical law represented

Units µm / pN / s throughout; viscosity in pN·s/µm² (1 Pa·s = 1 pN·s/µm² exactly).

### 4.1 The problem

Zero-Reynolds Stokes flow on the **unbounded** region exterior to a closed surface `S`:

```
grad p = mu * laplacian(u),   div u = 0   outside S
u -> 0 at infinity,           u = V(x) prescribed on S
```

No inertia, no wall, no far-field forcing, one body. `V` is the surface velocity, and it is a
difference of two configurations — the live surface and a *committed* reference — which is the whole
reason this component owns state at all.

### 4.2 Mobility, then resistance

The Stokes Green's function (Oseen tensor) is

```
G(r) = (1 / (8 pi mu)) * ( I / |r|  +  r r^T / |r|^3 )
```

so `u(x) = sum_j G(x - x_j) F_j` solves the Stokes equations for point forces `F_j`. Collocating that
representation at the `N` surface quadrature nodes gives the **mobility**

```
u_i = sum_j M_ij F_j,        M_ij = G(x_i - x_j)      [um / (pN s)]
```

which maps surface point forces to surface velocities — the registered representation's own sentence.
The **resistance** is its inverse, and the force the medium exerts back on the surface is

```
f = - R V,      R = M^-1                                [pN]
```

`R` is a `3N x 3N` operator, not a scalar and not a per-node coefficient. That is the entire content
of the registered `mechanical_role`: the six rigid-body modes carry zero internal stiffness, so
nothing else in the registry can load them, and a per-node coefficient that could would be a number
with no dependence on the cell's size or shape.

### 4.3 Why the kernel must be regularised, and why the width is not a knob

`G` is singular at `r = 0`, so `M_ii` needs a finite value, and the discrete operator must be
**symmetric positive definite** — `V^T R V` is the dissipated power, and a resistance with a negative
eigenvalue is a medium that can push a cell along.

Replace each point force by a smooth radial blob of unit weight and width `eps`. The regularised
Stokeslet for the algebraic blob is

```
S_eps(r) = (1 / (8 pi mu)) * ( (|r|^2 + 2 eps^2) I  +  r r^T ) / (|r|^2 + eps^2)^(3/2)
```

which tends to `G` as `eps -> 0` and is finite at the origin: `S_eps(0) = I / (4 pi mu eps)`.

**`eps` is then fixed by the quadrature, not chosen.** A node stands for a surface patch of area
`w`, and its self-mobility must be that patch's mobility. Integrating the *exact* Oseen tensor over
an equal-area disc of radius `rho = sqrt(w / pi)` in the tangent plane, with outward normal `n`:

```
integral_disc  I / |r|      dS = 2 pi rho * I
integral_disc  r^ r^^T      dS =   pi rho * (I - n n^T)
=> integral_disc G dS = (rho / (8 mu)) * (3 I - n n^T)
=> M_ii^patch = (rho / (8 mu w)) * (3 I - n n^T)      [per unit force at the node]
```

Matching the isotropic part (one third of the trace; `tr(3I - nn^T) = 8`) to the blob's
`S_eps(0)`:

```
(1/3) tr M_ii^patch = rho / (3 mu w)
(1/3) tr S_eps(0)   = 1 / (4 pi mu eps)
=>  eps = 3 sqrt(w) / (4 sqrt(pi)) = 0.42314218766... * sqrt(w)
```

That constant is `REGULARISATION_RATIO`. It is a **derived** number, it carries no fitted content,
and it vanishes with the quadrature under refinement — which is what makes the sphere limits below a
convergence statement rather than a calibration. For unequal patches the pair width is
`eps_ij^2 = (eps_i^2 + eps_j^2) / 2`, which keeps `M` exactly symmetric.

**Symmetry of `M` is Lorentz reciprocity at the discrete level**, and it is exact rather than
approximate: it holds to round-off, in every configuration, and a control asserts it.

### 4.4 The two closed-form sphere limits, derived here

For a rigid sphere of radius `a` translating at `U` in unbounded Stokes flow the surface traction is
uniform, `t = -3 mu U / (2 a)`, so the total force is `-3 mu U / (2a) * 4 pi a^2 = -6 pi mu a U`.

For the same sphere rotating at `Omega` the traction is `t = -3 mu (Omega x x) / a`, and

```
integral_S  x  x  (Omega x x) dS = integral_S ( Omega |x|^2 - x (x . Omega) ) dS
                                 = 4 pi a^4 Omega - (4 pi a^4 / 3) Omega
                                 = (8 pi a^4 / 3) Omega
=> torque = (3 mu / a) * (8 pi a^4 / 3) * Omega = 8 pi mu a^3 Omega
```

So:

```
translational resistance  =  6 pi mu a      [pN s / um]
rotational   resistance  =  8 pi mu a^3     [pN um s]
```

**Both are checked, and checking only the first is the specific mistake this module is built to make
visible.** A per-node isotropic drag `gamma_i = 6 pi mu a w_i / A` reproduces the translational limit
*exactly, to machine precision, by construction*, and gets the rotational limit wrong by a factor of
exactly two:

```
sum_i x_i x gamma_i (Omega x x_i) = gamma * N * (a^2 - a^2/3) * Omega = 4 pi mu a^3 Omega
```

That break is shipped as `use_isotropic_node_drag`, and the pair of assertions on it — translational
exact, rotational off by 2 — is the sharpest control in this module.

### 4.5 Exact geometric similarity

The scheme has no length scale of its own: `eps` scales with `sqrt(w)`, so scaling a quadrature by
`s` scales `M` by `1/s` and hence `R` by `s`. Therefore, for two geometrically similar quadratures,

```
translational resistance ratio = s^1      exactly (to round-off)
rotational    resistance ratio = s^3      exactly (to round-off)
```

This is stronger than the convergence statement and independent of it: it holds at *any* resolution,
including a coarse one where both absolute values are several per cent off. It is the direct
refutation of "a mobility unrelated to the cell's size or shape".

### 4.6 Dissipation, and why the medium banks work

The Rayleigh dissipation function is

```
D(V) = (1/2) V^T R V           [pN um / s]
f    = - grad_V D(V) = - R V   [pN]
P    = V^T R V = 2 D(V)        [pN um / s], the dissipated power
```

`f = -grad_V D` is this module's version of `F = -grad E`: the gradient identity holds in **velocity**
space, because the medium is dissipative and has no potential energy. `D` is exactly quadratic, so a
central difference in `V` is exact to round-off — see §10, which says what that does to the observed
order and what is asserted instead.

Over one accepted step of length `dt`, with `dx = x_live - x_committed_reference` and
`V = dx / dt`,

```
banked work  dW = P dt = dx^T R dx / dt  >= 0
```

and `cumulative_dissipated_work_pn_um` is the running sum over accepted steps. Non-negativity is not
assumed: it follows from `R` positive definite, which is *enforced* by the Cholesky factorisation the
solve uses. A quadrature whose mobility fails to factorise is refused with `MediumMobilityError`
rather than solved, because the alternative is a medium that can inject energy into the cell and
report a plausible traction while doing it.

### 4.7 Darcy is not an option outside the cell

The screened (Brinkman/Darcy) operator adds `-(mu / k) u` to the momentum balance, with screening
length `sqrt(k)` set by the solid volume fraction. Outside the cell the solid fraction is **zero**, so
`k -> infinity` and `sqrt(k) -> infinity`: the screening term vanishes identically and there is no
small parameter left to expand in. Darcy outside the cell is therefore **wrong rather than
approximate** — it is a different equation whose defining scale does not exist here, not a coarse
version of this one. `assert_medium_is_not_darcy` refuses a zero or negative solid fraction for
exactly this reason and returns quietly for a positive one, which is the *intracellular* porous case
and is `cytosol`'s business, not this component's.

### 4.8 One dissipation channel on the outer surface, and only one

The registered `mechanical_role` states that this component is the sole owner of
velocity-proportional dissipation on the outer surface, and that *"a per-node drag applied in addition
to the medium is not a refinement, it is the same dissipation counted twice"*.

`declare_exterior_drag_channel(owner, channel)` is the exterior analogue of
`aleph.vertical.cytosol.CytosolField.declare_drag_channel`, which guards the same class of defect on
the interior. It refuses at **declaration** time, before any number exists, because double-counted
velocity-proportional dissipation shows up as an error that grows with speed — which reads as a
velocity-dependent constitutive law rather than as a defect. The surface owner this medium was
constructed against may only ever declare `EXTERIOR_STOKES`; a `PER_NODE_DRAG` declaration from it is
refused outright.

## 4a. Mutation testing of the controls

Three mutants applied to `aleph/vertical/medium.py`, the controls run against each, the source then
restored and verified byte-identically. Recorded because a control suite no mutant can break is
decoration.

| # | Mutation | Tests killed | Survived? |
|---|---|---|---|
| M1 | `regularisation_um`: drop the derived `3 / (4 sqrt(pi))` factor, returning `sqrt(w)` | **7** | no |
| M2 | `exterior_mobility_um_per_pn_s`: use `eps_i^2` alone instead of the symmetric mean `(eps_i^2 + eps_j^2) / 2` | **3** | no |
| M3 | `commit`: bank the dissipated work but do **not** advance the committed surface reference | **2** | no |

Source restored and verified **byte-identical**: SHA-256
`354ec41bf4dfba191838318c430358c944705a1a5d3ce30f6b7145825ece4897` before and after, taken and
re-taken by the same script that applied the mutants. The hash is the check rather than `git diff`,
because the module's committed revision is not this lane's finished one — see §13b.

M1 kills: both sphere-limit controls, the isotropic-drag negative control, the derived-width
control, the written-out kernel block, the Oseen limit, and the non-uniform-quadrature limit.

M2 kills: the written-out kernel block, the symmetric-pair-width control, and the mobility symmetry
on a non-uniform quadrature.

M3 kills: the reference-advances control and the constant-rate dissipation control.

**The exercise found two real gaps, which is the reason to run it.**

**M2 was written before its controls existed and would have survived them.** Every fixture in the
first draft of the control file was a `fibonacci_sphere_quadrature`, and that helper gives every node
the **same** patch area — on equal weights the symmetric mean `(eps_i^2 + eps_j^2)/2` and the
one-sided `eps_i^2` are the same number, so the mutation is invisible and *every one of the 75
controls passed*. A suite that only ever sees equal weights cannot see a defect in how unequal
patches are combined. Two controls were added to close it: a Gauss-Legendre-by-uniform-azimuth sphere
quadrature whose weights vary by more than a factor of three, and a direct block-by-block comparison
of the kernel against the closed form on two nodes with different weights. The second is the sharper
one, and it is the kind of control the sphere oracles cannot substitute for: an integral quantity can
absorb an error in one factor of the kernel and still converge.

**M3 was killed by exactly one control**, which is not enough for a defect this quiet. A `commit`
that banks the work but never moves the reference produces a measured velocity that grows linearly
with the step index and a dissipation that grows quadratically, and every individual force still
looks like an ordinary drag. `::TestTheMediumIsAStateOwnerAndNotABoundaryCondition::
test_a_surface_moving_at_constant_velocity_dissipates_at_a_constant_rate` was added afterwards and
states the physics directly: a body translating at constant velocity through a Newtonian medium
dissipates at a constant rate, step after step, and that is true only if the reference tracks the
accepted configuration.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| quadrature point | µm | 1e-6 m | finite |
| quadrature weight (patch area) | µm² | 1e-12 m² | `> 0` |
| outward normal | — | — | unit, finite |
| `viscosity_pn_s_per_um2` | pN·s/µm² | 1 Pa·s | `> 0` |
| mobility `M` | µm/(pN·s) | m/(N·s) | symmetric, positive definite |
| resistance `R = M^-1` | pN·s/µm | N·s/m | symmetric, positive definite |
| translational resistance | pN·s/µm | N·s/m | `> 0` |
| rotational resistance | pN·µm·s | N·m·s | `> 0` |
| `surface_traction_pn` | pN | 1e-12 N | finite |
| `cumulative_dissipated_work_pn_um` | pN·µm | 1e-18 J | monotone non-decreasing |

`surface_traction_pn` carries a **force** unit, not a stress unit, and the registered key name says
so. It is the traction integrated over the node's quadrature patch — a nodal force at a quadrature
site. The name is reproduced verbatim from the registry, including that mild infelicity, because a
module that silently improves a registered key name is a module that has stopped agreeing with the
registry.

Singular and boundary cases, each with the behaviour Aleph requires:

- **Two coincident quadrature nodes**: refused with `MediumQuadratureError`. Two nodes at one point
  are one point counted twice; the mobility acquires two identical rows and the resistance does not
  exist.
- **A non-positive quadrature weight**: refused. A node standing for zero area has no patch, and
  `eps = 0` puts the singular kernel back on the diagonal.
- **A surface that is not closed** (`|sum_i w_i n_i| / A` above `closure_tolerance`): refused. An open
  sheet does not separate an interior from an exterior, and the exterior Stokes problem is posed on
  the complement of a body. A hemisphere scores `0.5` on this measure and a closed sphere scores
  `<= 2.3e-2` at every resolution tested, so the discrimination is not marginal.
- **A quadrature enclosing non-positive volume** (`(1/3) sum_i w_i x_i . n_i <= 0`): refused. Inward
  normals would make the medium resist motion in the wrong sense everywhere.
- **A mobility that fails to factorise**: refused with `MediumMobilityError`, never regularised
  further. See §4.6.
- **Zero surface velocity**: returns **exactly** `0.0` force at every node and exactly `0.0`
  dissipated power. Exactly, testable with `==`, because `R` applied to a zero vector is zero and no
  branch intervenes.
- **`dt <= 0`**: refused by `StepContext` upstream and independently here.
- **A traction stencil requested for an owner this medium was not built against**: refused with
  `MediumRoleError`. The resistance is an operator on *that* surface.

Invariants, each with the test that asserts it:

- **I1.** `f = -grad_V D(V)`, exactly (a quadratic form differentiated by an exact central difference)
  — `tests/vertical/test_medium_controls.py::TestTheResistanceIsTheGradientOfItsDissipationFunction::test_the_force_is_minus_the_velocity_gradient_of_the_dissipation_function`
- **I2.** `M` and the `6x6` rigid-body resistance are symmetric — Lorentz reciprocity —
  `::TestTheResistanceOperatorIsWhatItClaimsToBe::test_the_mobility_is_symmetric`,
  `::test_the_rigid_body_resistance_is_symmetric`
- **I3.** `R` positive definite, so the dissipated power is `> 0` for every non-zero velocity —
  `::TestTheResistanceOperatorIsWhatItClaimsToBe::test_the_dissipated_power_is_positive_for_every_velocity`
- **I4.** Both sphere limits converge under surface refinement at a measured order —
  `::TestTheSphereLimitsAreReproduced::test_the_translational_resistance_converges_to_six_pi_mu_a`,
  `::test_the_rotational_resistance_converges_to_eight_pi_mu_a_cubed`
- **I5.** Exact geometric similarity: `s^1` and `s^3` —
  `::TestTheSphereLimitsAreReproduced::test_the_resistance_scales_with_the_first_and_third_power_of_the_radius`
- **I6.** Rollback restores every owned block bit-identically —
  `::TestOwnershipAndTheAcceptedStepContract::test_a_rejected_step_restores_every_owned_block_bit_identically`
- **I7.** Purity and determinism: no evaluation mutates its inputs or the owner's arrays
  (`np.array_equal`, not `allclose`) —
  `::TestTheEvaluationsArePureAndDeterministic::test_no_evaluation_mutates_its_arguments`
- **I8.** Declared state keys equal the registered contract, in registered order —
  `::TestTheDeclarationsAgreeWithTheCode::test_the_owned_state_keys_are_exactly_the_registered_contract`
- **I9.** No state key of this owner is declared by any other registered owner —
  `tests/state/test_census_key_ownership.py::test_no_state_key_has_two_registered_owners`
- **I10.** The mobility is symmetric on a quadrature with **unequal** patch areas, and the kernel
  block equals the closed form of §4.3 written out — `::TestTheRegularisedKernelIsTheOneDerived::
  test_the_mobility_block_matches_the_regularised_stokeslet_written_out`,
  `::TestTheOperatorDoesNotDependOnTheQuadratureRule::
  test_the_mobility_is_symmetric_on_a_non_uniform_quadrature`. Added in response to mutant M2; see
  §4a.
- **I11.** A surface translating at constant velocity banks a constant amount of work per accepted
  step — `::TestTheMediumIsAStateOwnerAndNotABoundaryCondition::
  test_a_surface_moving_at_constant_velocity_dissipates_at_a_constant_rate`. Added in response to
  mutant M3; see §4a.

## 6. Source evidence class and known retractions

**No claim is made about the reference implementation, because none of it was read**, and in this
case the ignorance is unusually complete: this lane could not even name a candidate source file. Its
evidence class is **unknown here**. Where I looked: nowhere inside `/Users/sw1/ffn_cellsim` —
deliberately, so this derivation could not be contaminated.

The registered `extracellular_medium` contract carries `citation_status: UNSOURCED`, and nothing here
changes that. `default_medium_card()` returns the viscosity of water at room temperature expressed in
Aleph units, which is a textbook constant and not a measurement of any culture medium; the card says
so and no run may read a medium property off it.

The registered `membrane_medium_traction` contract carries `implemented=False`, and **this entry does
not change that flag and did not touch that file**. This module supplies the exterior solve the
connector's `mechanical_interpretation` says does not exist; the connector itself remains unwritten
and unwired.

## 7. Independent oracle or derivation

Four, none of which is the reference implementation:

1. **`6 pi mu a`**, the translational Stokes drag, derived in §4.4 from the exact traction field. It
   is a property of the whole operator, so no rescaling of one matrix entry reproduces it.
2. **`8 pi mu a^3`**, the rotational Stokes drag, derived the same way. It is the sector the
   translational check cannot see, and §4.4 exhibits a wrong law that passes the first exactly and
   fails this one by a factor of two.
3. **Exact geometric similarity**, `s^1` and `s^3` (§4.5). Independent of resolution and of both
   absolute values, so it survives at a coarse quadrature where 1 and 2 are only a few per cent.
4. **Central finite differences** of the module's own dissipation function in velocity space — an
   oracle for the gradient claim that depends on no constant in the law.

Lorentz reciprocity (`M = M^T`) and positive definiteness are structural checks rather than oracles:
they can be satisfied by a wrong operator, which is why they are listed under §5 invariants and not
here.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_medium_controls.py::TestTheSphereLimitsAreReproduced::test_the_translational_resistance_converges_to_six_pi_mu_a` | Measured translational resistance approaches `6 pi mu a` monotonically at `N = 50, 200, 800`, with observed order `> 0.8` in `h = sqrt(A/N)` |
| Positive | `tests/vertical/test_medium_controls.py::TestTheSphereLimitsAreReproduced::test_the_rotational_resistance_converges_to_eight_pi_mu_a_cubed` | Same, against `8 pi mu a^3`. Reproducing only the translational limit leaves this whole sector unmeasured |
| Positive | `tests/vertical/test_medium_controls.py::TestTheSphereLimitsAreReproduced::test_the_resistance_scales_with_the_first_and_third_power_of_the_radius` | Doubling the radius multiplies the two resistances by exactly `2` and exactly `8`, to `1e-12` relative |

Measured numbers are recorded in §8a after the run. Nothing in this section is quoted from memory.

## 8a. Measured, this session, on this tree

`/Users/sw1/miniconda3/envs/aleph/bin/python`, numpy float64, 2026-07-30/31. Every number below was
produced by a run in this session; none is quoted from memory. Fixture: `mu = 4e-3 pN s/um^2`,
`a = 3 um`, so `6*pi*mu*a = 0.226194671 pN s/um` and `8*pi*mu*a^3 = 2.714336053 pN um s`.

Convergence under surface refinement, `h = sqrt(A / N)`:

| N | h (µm) | translational (pN·s/µm) | rel. err | rotational (pN·µm·s) | rel. err |
|---|---|---|---|---|---|
| 50 | 1.5040 | 0.231857 | 0.02503 | 2.975368 | 0.09617 |
| 200 | 0.7520 | 0.228752 | 0.01131 | 2.844281 | 0.04787 |
| 800 | 0.3760 | 0.227324 | 0.00499 | 2.778481 | 0.02363 |

**Observed order in `h`: 1.147 and 1.179 for translation (mean 1.163); 1.006 and 1.019 for rotation
(mean 1.012).** Both sequences are monotone and both approach from above. The scheme is first order,
which is what the flat-patch matching argument in §4.3 predicts, and the absolute error at the finest
quadrature is 0.5% in translation and 2.4% in rotation — reported, not tuned away.

Exact geometric similarity, `N = 200`, radius doubled from 3 µm to 6 µm:

| Quantity | measured ratio | exact |
|---|---|---|
| translational | `2.000000000000000` | 2 |
| rotational | `7.999999999999997` | 8 |

The isotropic-node-drag break, same quadrature:

| Quantity | measured / exact |
|---|---|
| translational | `1.000000000000003` |
| rotational | `0.500000000000000` |

That is the whole argument for checking both limits, in two lines: the wrong law reproduces the first
to fifteen digits and is out by a factor of exactly two on the second.

Structural properties, `N = 200`:

| Property | measured |
|---|---|
| mobility asymmetry, relative | `0.000e+00` (exactly symmetric) |
| 6×6 grand-resistance asymmetry, relative | `1.467e-16` |
| 6×6 eigenvalues | `0.228746 0.228752 0.228759 2.843898 2.844250 2.844696` — all six strictly positive |
| translation–rotation coupling / diagonal | `1.812e-04` |
| smallest mobility eigenvalue | `1.981e+01` (positive definite) |

The six positive eigenvalues are the registered contract's headline claim measured: three
translational and three rotational rigid-body modes, all loaded, none driven.

Gradient identity in velocity space, `N = 50`, `D = 0.689624 pN·µm/s`, force scale `0.040733 pN`:

| step (µm/s) | max&#124;f + ∇_V D&#124; (pN) | relative |
|---|---|---|
| 4e-4 | 3.189e-13 | 7.828e-12 |
| 2e-4 | 9.392e-13 | 2.306e-11 |
| 1e-4 | 1.165e-12 | 2.860e-11 |

The residual **grows** as the step shrinks, which is the signature of a difference that is exact and
round-off limited rather than truncation limited — and is exactly why no observed order is asserted
here: computed from these three numbers it would be *negative*. See §10.

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_medium_controls.py::TestTheNegativeControlsFail::test_an_isotropic_node_drag_passes_the_translational_limit_and_fails_the_rotational_one` | With `use_isotropic_node_drag=True` the translational resistance is `6 pi mu a` to `1e-12` relative **and** the rotational resistance is half of `8 pi mu a^3` to `1e-4`; with the flag off the same quadrature gets both within the discretisation error |
| Negative (must fail) | `tests/vertical/test_medium_controls.py::TestTheNegativeControlsFail::test_a_flipped_resistance_sign_makes_the_medium_inject_energy` | With `flip_resistance_sign=True` the banked step work is strictly negative; with the flag off it is strictly positive, and the two differ only in sign |
| Negative (must fail) | `tests/vertical/test_medium_controls.py::TestTheNegativeControlsFail::test_omitting_the_committed_reference_leaves_the_rigid_body_modes_unloaded` | With `omit_committed_reference=True` a translating and rotating surface produces **exactly** `0.0` force at every node and banks **exactly** `0.0` work; with the flag off the same motion produces a finite force and a strictly positive banked work |

All three breaks are flags on the **shipped** owner, so each control drives the real code path. Each
is asserted in **both** directions: a control that only checks the broken branch misbehaves would
pass for an owner that misbehaves unconditionally.

## 10. Numerical and precision envelope

Working and accumulation precision: float64 throughout; no reduced-precision path exists. The solve
is a dense Cholesky factorisation of the `3N x 3N` mobility, `O(N^3)`, and the largest quadrature any
control builds is `N = 800` (a `2400 x 2400` factorisation, 46 MB, well under a second).

**The dissipation function is exactly quadratic in the velocity**, so its central difference is exact
and the finite-difference residual sits at round-off for every step size — measured at `8e-12` to
`3e-11` relative to the force scale (§8a), which is the round-off of `D` itself divided by a force
scale an order of magnitude below `D`, not a truncation error.
Reporting an "observed order" from a sequence of round-off-level errors would be reporting noise, so
the gradient control asserts the *residual* against a round-off bound (`< 1e-10` relative), which is a
strictly stronger statement than order two, and says in its own docstring why the order statistic is
absent. **This is a deliberate departure from the house per-term gradient control** and is recorded
here rather than left for a reader to notice. The house control's other requirement is kept and is
load-bearing: `D > 0` is asserted before the comparison, so a dissipation function that returned zero
cannot agree with its own finite difference and pass.

The measured convergence order that *is* reported comes from the sphere limits under surface
refinement (§8), where the error is a genuine discretisation error and the order is meaningful.

Tolerances and why they are not looser:

- Sphere limits: asserted as monotone decrease plus observed order `> 0.8`, against a scheme whose
  error is `O(h)`. The absolute error at `N = 800` is a few tenths of a per cent to a couple of per
  cent depending on sector, which is why the assertion is on the *trend* and not on a single value —
  a single-value tolerance loose enough to pass at `N = 50` would pass mutant M1 as well.
- Geometric similarity: `1e-12` relative. It is an exact algebraic identity of the scheme, not a
  physical approximation, so anything looser would be hiding a defect behind a physics-shaped
  tolerance.
- Symmetry and reciprocity: `1e-14` of the operator's own scale. Same argument.
- Positive definiteness: the smallest eigenvalue is required strictly positive, and the Cholesky
  factorisation in the shipped code is the enforcement rather than the test.

Outside the envelope: every singular case in §5 refuses with a typed exception. There is no input
range over which this module returns a plausible number it cannot stand behind.

## 11. Production-backend residency and transfer

Host, numpy, float64. The module imports `EntityCensus`, `Phase`, `StepContext` and `WorkReceipt`
from `aleph.runtime.participant` and uses all four. `accumulate` uses `ctx.backend` for the array
operations in `ACCUMULATE_FIELD` so the coverage witness advances; the dense factorisation itself is
numpy and is **not** routed through the `Backend` protocol, which has no linear-solve entry point.
That is a real limitation and it is named here rather than in a docstring: a GPU backend would today
see this owner's factorisation as a host round-trip.

No GPU work of any kind was run by this lane, and no authorization was sought or held.

Cost note for whoever assembles: the mobility is dense and `O(N^2)` in memory. A live membrane
quadrature with tens of thousands of nodes is **not** reachable with this implementation, and a fast
summation method is the obvious next port rather than a tuning exercise.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered this module, because no reference file was read.
`tests/vertical/test_medium_controls.py::TestTheDeclarationsAgreeWithTheCode::test_no_reference_project_vocabulary_survives`
asserts the module text carries none of the provider tokens, and
`::test_the_module_imports_nothing_from_validation` asserts the `aleph/** -> validation/**` firewall.

One piece of Aleph's own prose is corrected by this module landing: `aleph/vertical/ecm.py:77` refers
to `aleph.vertical.medium` as though it existed, and until now that cross-reference dangled. **This
lane did not edit `ecm.py`** — the reference simply stops being broken.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 84 controls pass (`tests/vertical/test_medium_controls.py`, 2026-07-31, 5.15 s). Three mutants applied, killing 7 / 3 / 2 tests; none survived; source restored byte-identically. Status stays `PROPOSED` — see §13a and §14. |
| Reviewer | Agent-proposed (lane L32). Unratified. No PI review. |
| Rollback | Delete `aleph/vertical/medium.py` and `tests/vertical/test_medium_controls.py`. Nothing breaks: no module in the tree imports either, and no connector reaches them. `aleph/vertical/ecm.py:77`'s cross-reference returns to dangling, which is where it was. |

## 13a. Acceptance result, and the runs behind it

```text
tests/vertical/test_medium_controls.py            84 passed in 5.15s
tests/state/ tests/ports/                        865 passed, 1 failed in 7.80s
```

**The one remaining failure is not this lane's.** `test_index_is_not_stale` reports that
`ports/ledger/INDEX.md` no longer matches the entries on disk, because six new entries — `-3201`
through `-3206` — landed from six concurrent lanes without a regeneration. `INDEX.md` is a shared
file that belongs to no session, so this lane deliberately did not regenerate it: whoever commits
this batch should run `python ports/regenerate_index.py` **once** for all six rather than six times
in six commits, each of which would then conflict with the last.

Measured at the start of this session, before any of this entry's files existed, `tests/ports/` was
already at `2 failed, 864 passed` — the second failure being
`test_named_controls_resolve_to_real_tests`, which the sibling owner lanes cleared during this
session by landing the control files their entries had named. **`ALEPH-PORT-3206` was never among
them**: every control this entry names is defined in `tests/vertical/test_medium_controls.py` and
resolves.

### Deviations from the plan this entry set out before the code

Recorded because an entry written before the code that never says where the code diverged is an
entry nobody can check against.

1. **§4a's predicted kill counts were wrong, and the corrected table is above.** M1 killed 7 rather
   than 5, M2 killed 3 rather than 4, M3 killed 2 rather than 4 — and M3 killed only **1** until a
   control was added in response. The numbers in the first draft were estimates written before the
   controls existed, which is the honest cost of writing §4a in advance; they are replaced by
   measurements rather than left standing.
2. **The mutation exercise changed the control suite**, in the way §4a says it is supposed to. Two
   controls (the non-uniform quadrature and the constant-rate dissipation) exist because a mutant
   would otherwise have survived or been under-detected. Both are now in §5's invariant list.
3. **No API listed in §1 was renamed or dropped**, and nothing was added to it.

## 13b. This module was committed while its author was still writing it

Recorded here because the mechanism costs more than the incident, and because the ledger is where a
future reader will look to find out why the committed revision and the finished one differ.

`7564509 feat(vertical): the six missing compartments` and `b299a3a fix(tests): repair
test_medium_controls.py, which I committed mid-write` both landed **during** this lane's run. The
lead session ran the control file, saw 84 pass, and committed — and this lane edited the file
between that run and the `git add`, so what landed was a mid-edit state with an unterminated
f-string. `b299a3a` repaired it, arriving independently at the same two-line fix this lane had
already made in its own working tree.

Nothing was lost and nothing is broken. The working tree at the end of this lane is the finished
version: it is a strict superset of `b299a3a`, it carries the line-length pass and one restored
conjunction in a refusal message that the committed revision does not, and it is the revision the
84 controls and the three mutants above were measured against. The next commit of these paths should
take the working tree as it stands.

The lead's own commit message states the rule this suggests, and it is worth keeping: **verify and
stage in one step, or do not commit a file whose author has not reported finished.** Staging by
explicit path does not help here — the hazard is timing, not staging. `docs/ACTIVE_SESSIONS.md`
records the same shape twice already, both times arriving through `git add -A`; this is the third
occurrence and the first through a different door.

## 14. Honest limits — what this entry does NOT establish

**Status is `PROPOSED`, not `ACCEPTED`.** The controls pass and the oracles are Aleph's own, but no PI
has reviewed this, the connector it unblocks is not wired, and the owner has never been stepped inside
a world.

**Wired connector count contributed by this module: zero.** `membrane_medium_traction` still carries
`implemented=False` and this lane did not change it. What this module supplies is the *medium-side*
endpoint and the exterior solve the connector needs; the membrane-side endpoint, the connector object,
the adjoint pair and the `Binding` row are all still absent. A registry reader must not read this
entry as the connector having landed.

**No number here is a property of any real culture medium.** Evidence rung `ANALYTIC_ORACLE`,
quantitative status `BLOCKED`. The viscosity default is water; the cell radii are chosen to make the
algebra sharp. `citation_status` for the `extracellular_medium` contract stays `UNSOURCED`.

**The convergence is first order and the absolute error at the resolutions tested is not small.** At
`N = 800` on a sphere the two limits are reproduced to a few per cent, not to a few parts in a
thousand. That is the scheme's honest accuracy, it is reported rather than tuned away, and it means
this operator is fit for checking a law and not yet for quoting a drag coefficient.

**Positive definiteness is enforced, not proved.** Every configuration this lane built factorises, and
the Cholesky refusal is the guard. This entry does **not** claim a theorem that the regularised
single-layer collocation is positive definite for every closed surface and every quadrature. A
surface that fails is refused, which is the correct behaviour, but it would also be a finding.

**The quadrature is an input and nothing in Aleph produces one.** `fibonacci_sphere_quadrature` is a
build helper for spheres. A live membrane surface quadrature — the thing the registered
representation actually calls for — comes from the membrane owner, which does not publish one. So
"evaluated on the live cell surface quadrature" is implemented as *"evaluated on whatever quadrature
it is handed"*, and the live half is untested because there is nothing live to hand it.

### 14.1 The registered `unsupported_claims`, reproduced as stated absences

Every one of the six is absent from this implementation, in the registry's own terms:

1. **No wall and no half-space image system.** The Green's function used is the free-space Stokeslet.
   There is no image system, therefore no near-substrate lubrication, no cell-substrate hydrodynamic
   film, and no asymmetric exterior of any kind. A cell approaching a substrate in this model feels
   exactly the same medium it feels in free space, which is wrong near a wall and is not corrected.
2. **No inertia.** The representation is the zero-Reynolds limit. There is no unsteady term, no
   added mass, and no acoustic response. Nothing here has a time constant of its own; the only time
   in the operator is the `dt` used to convert a configuration difference into a velocity.
3. **No non-Newtonian medium.** The viscosity is a single scalar constant. No shear thinning, no
   viscoelastic exterior, no polymer depletion layer, and no dependence of the operator on the
   deformation rate.
4. **No thermal fluctuation, and therefore no fluctuation-dissipation statement.** This is the
   absence most likely to be misread, because the module *does* own a positive-definite resistance and
   a banked dissipated work, which is exactly the pair a fluctuation-dissipation relation needs. It
   is not one. There is no noise term, no temperature, and no stochastic integrator here, and this
   module makes no claim about the amplitude of any fluctuation. A passive fluctuation observation
   must get its noise from a declared stochastic integrator; taking it from this operator would be
   inventing a theorem this entry does not have.
5. **No medium chemistry.** No osmolarity, no solute transport, no ion content. The osmotic condition
   is a declared experimental context elsewhere and is not a property of this component; nothing in
   this module can be read as an osmotic statement.
6. **No free surface, no air interface, and no far-field flow.** The domain is unbounded with
   `u -> 0` at infinity. There is no imposed shear, no Poiseuille background, no meniscus and no
   interface of any kind.

### 14.2 Two more absences worth naming

**Only one body.** The exterior problem is posed around one closed surface. Two cells in this medium
would interact hydrodynamically and nothing here represents that; a second quadrature would need to
enter the same mobility, not a second medium.

**No coupling to the interior fluid.** `cytosol` is a separate owner with a separate representation,
and the two are not connected by anything this module supplies. The surface is a boundary for both,
and whether their tractions are consistent at that surface is a question no test in this lane asks.

**Not run:** any GPU job (zero), and the full suite. `tests/vertical/test_medium_controls.py`,
`tests/state/` and `tests/ports/` were exercised. Five other lanes were writing concurrently, so a
whole-repo run would conflate a foreign breakage with mine — and `tests/ports/` already carries two
failures that are not this lane's, recorded in §13a.
