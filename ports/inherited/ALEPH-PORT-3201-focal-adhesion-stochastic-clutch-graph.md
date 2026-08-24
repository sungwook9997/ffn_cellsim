# ALEPH-PORT-3201 — `focal_adhesion` as a geometry-less stochastic clutch/joint graph

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3201` |
| Lane | `L32 focal_adhesion` (rebuild under a new id; `ALEPH-PORT-2901` is `REJECTED`) |
| Status | `PROPOSED` |
| Written | `2026-07-30`, **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |

---

## 0. Why this entry exists at all, and why the id is new

`ALEPH-PORT-2901` is `REJECTED`. It described 94 KB of `focal_adhesion` physics and named six
control tests **in two test files that never existed**, and the port-discipline gate did not catch
it because the gate then skipped entries below `ACCEPTED`. The code is in
`_archive/focal_adhesion_unverified_2026-07-30/` and its README forbids copying from it. **Nothing in
this entry, and nothing in the module it authorises, was taken from that archive.** The module was
written against the registered `focal_adhesion` contract in `aleph/state/census_actomyosin.py:670`
and the two connector contracts in `aleph/state/connectors_surface_traction.py:491,518`, which are
Aleph's own documents and were not what was wrong.

Every control named below is a test that exists in
`tests/vertical/test_focal_adhesion_controls.py`, and the gate
`tests/ports/test_port_discipline.py::test_named_controls_resolve_to_real_tests` now checks that at
every status, which is the specific hole 2901 fell through.

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.vertical.focal_adhesion import (
    ACTIN_SIDE_OWNERS,
    ALPHA2BETA1_SERIES_GROUP,
    ActinSideReference,
    BindingState,
    COMPONENT_NAME,
    Clutch,
    ClutchCard,
    ClutchEndpoint,
    ClutchEndpointRole,
    ClutchTransition,
    FA_ENDPOINT_ROLES,
    FocalAdhesionClutchGraph,
    FocalAdhesionOwnershipError,
    FocalAdhesionRoleError,
    FocalAdhesionSeriesError,
    FocalAdhesionStaleHandleError,
    FocalAdhesionStateError,
    LIGAND_SIDE_OWNER,
    LigandSideReference,
    MATURATION_ORDER,
    MaturationState,
    OWNED_STATE_KEYS,
    PUBLISHABLE_ENDPOINT_ROLES,
    ROLE_MATURATION_STATES,
    assert_maturation_is_bound,
    assert_series_group_is_dispatched_whole,
    assert_state_keys_disjoint_from,
    default_clutch_card,
    owned_state_keys,
    series_joint_energy_and_forces,
    series_stiffness_pn_per_um,
)
```

This list is `__all__` verbatim, and
`::TestTheDeclarationsAgreeWithTheCode::test_every_exported_name_resolves` plus
`::test_every_public_symbol_is_exported` pin it in both directions — `ecm.py` shipped an `__all__`
naming five symbols that did not exist, and to a reader scanning an export list each one read as
shipped code.

`aleph/vertical/focal_adhesion.py` is the Aleph target. It is a **state owner**: it allocates the six
registered state blocks, implements `snapshot` / `commit` / `rollback`, `owned_entities` and
`accumulate`, and publishes endpoints. It wires **zero** connectors — a `Binding` row in
`aleph/vertical/wiring.py` is added in the same commit that wires one, never before, and this lane
does not touch that file.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — verified to resolve, read-only, by `git rev-parse` only |
| Source path | `ffn_sim/ac/engine/adhesion*.py` — **named, not read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened by this lane. Nor was any file under `_archive/focal_adhesion_unverified_2026-07-30/` other than its `README.md`. |
| Working tree == commit? | not applicable — nothing was read, so there is no revision to have mis-cited. The provider working tree reports 31 dirty paths, which is recorded only so a later auditor knows the tree and the commit are not the same thing. |

The source path is recorded so a later auditor can establish independently whether Aleph's
re-derivation agrees with the reference. Listing it is **not** a claim that it was consulted.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** Two linear springs in series is undergraduate mechanics; a
Markov state on a bond is a two-line definition. Per `ports/TEMPLATE.md` §3 a law anyone competent
re-derives in an afternoon is a reason to write clean-room, not a reason to port. The interesting
content here is not the algebra — it is *which* algebra the registered contract forbids, and that is
an Aleph document.

## 4. Physical or mathematical law represented

Units µm / pN / pN·µm / s throughout.

### 4.1 The representation is geometry-less, and that is load-bearing

The registered `numerical_representation` is emphatic: *"A GEOMETRY-LESS STOCHASTIC CLUTCH/JOINT
GRAPH. There is no adhesion mesh, no adhesion particle, and no adhesion volume."* So this owner has
**no node array, no position array and no material point of its own.** A clutch carries a *binding
state*, not a place. Every position this module ever touches arrives as an argument from another
owner's endpoint and is never stored.

The consequence for the load path is exact and it is the whole design: a body with no degrees of
freedom cannot hold a force. The plaque is a **massless, position-less middle**. Whatever load enters
the actin side must leave the ligand side, because there is nowhere in between for a difference to
be stored.

### 4.2 The series joint, re-derived

One clutch spans an actin-side anchor point `x_a` (a material point owned by `sf_arc`,
`lamellipodium` or `filopodium`) and a ligand-side point `x_l` (a material point owned by `ecm`).
Neither point is this owner's. Write `d = |x_l − x_a|`, `u = (x_l − x_a)/d`, and `e = d − d0` for a
rest gap `d0 > 0`.

The clutch is two Hookean elements in series — the talin/integrin anchoring on the actin side with
stiffness `k_a`, the α2β1/collagen bond on the ligand side with stiffness `k_l`. Series means two
statements, and both are needed:

1. **One load.** `T = k_a e_a = k_l e_l`.
2. **Additive extensions.** `e_a + e_l = e`.

Eliminating `e_a`, `e_l`:

```
1/k_series = 1/k_a + 1/k_l          k_series = k_a k_l / (k_a + k_l)
T = k_series e                      e_a = e k_l/(k_a+k_l)   e_l = e k_a/(k_a+k_l)
E = ½ k_series e²
```

`k_series < min(k_a, k_l)` strictly, for any positive pair: **a series joint is softer than either
member, and it is governed by whichever member is softer.** The forces are the exact gradient:

```
f_anchor = +T u        f_ligand = −T u
```

so the pair sums to exactly zero and carries exactly zero moment about any origin. The moment
argument is worth writing out, because a first draft of this entry got it wrong and §8a records the
correction: the net moment is `(T_a x_a − T_l x_l) × u`, which collapses to
`T (x_a − x_l) × u = −T d (u × u) = 0` **only because the two loads are equal**. Moment closure is
therefore not a free consequence of the two forces being parallel to the joint axis; it is a
consequence of the series law itself, and it fails when the series law does.

### 4.3 What the independent-spring evaluation does, in closed form

This is the defect the registry spends two paragraphs forbidding, and this module's job is to make
it *fail a control* rather than be warned about in a docstring.

If `fa_actin_anchor` and `integrin_collagen_clutch` are dispatched as two independent springs, each
half is evaluated against its own kinematics. The middle has no position, so each half sees the same
pair of endpoints and therefore the same extension `e`:

```
T_a = k_a e        T_l = k_l e        E_indep = ½ (k_a + k_l) e²
```

Three consequences, each a separately checkable number:

- **The two halves carry different loads.** `T_a − T_l = (k_a − k_l) e ≠ 0` whenever `k_a ≠ k_l`.
  The load applied to the actin side and the load applied to the ligand side no longer balance, so
  the net force on the isolated joint is `(k_a − k_l) e`, injected into the world from a body with
  no degrees of freedom. **The plaque has become a hidden force source**, which is the registry's own
  phrase.
- **The joint is stiffer than either member.** The energy curvature is `k_a + k_l`, the *parallel*
  combination, and `k_a + k_l > max(k_a, k_l) > k_series`. The overstatement factor is
  `(k_a + k_l)² / (k_a k_l) ≥ 4`, with equality only at `k_a = k_l`. Stiffer joints read as higher
  traction, so this is not a bookkeeping error — it is a traction number that is wrong upward and
  cannot be wrong downward.
- **The gate stops gating.** With the clutch disengaged the series path is open and must carry
  **exactly** zero. Evaluated independently, the actin half still reports `k_a e` while the ligand
  half reports nothing: a traction reading with nothing on the other end of it.

Closed-form oracle used by the controls, at `k_a = 100 pN/µm`, `k_l = 300 pN/µm`:
`k_series = 75`, `k_a + k_l = 400`, overstatement `= 5.3333…`.

### 4.4 Occupancy, maturation, and why neither is decorative

A clutch is a *bundle* of molecular bonds, not one bond. `focal_adhesion_occupancy` is an integer:
how many integrin/collagen bonds of that clutch are currently engaged, in `[0, bond_capacity]`.
Bonds in a bundle are in **parallel**, so both stiffnesses are proportional to occupancy:

```
k_l = n · bond_stiffness_pn_per_um
k_a = n · anchor_stiffness_pn_per_um · m(maturation)
```

Two properties follow, and both are asserted:

- **`k_series` is linear in `n`.** So splitting one clutch of occupancy `2n` into two clutches of
  occupancy `n` at the same pair of endpoints gives exactly the same total load. This is the
  clutch-graph analogue of resolution invariance: the mechanics is a statement about how many bonds
  are engaged and not about how they were grouped into clutch records.
- **Maturation is bound to the mechanics.** `m(NASCENT) < m(MATURING) < m(MATURE)` strictly, so
  advancing an adhesion's maturation state changes the load it carries at fixed geometry and fixed
  occupancy. `assert_maturation_is_bound` refuses a graph where it does not. Without this, a
  maturation state that nothing reads conserves energy perfectly and is invisible to every gradient,
  closure and conservation check — the same class of defect as the cortex's unbound axial law
  (`ALEPH-PORT-2301`).

**`m` is a placeholder and is labelled `UNSOURCED` in the code.** The registered contract owes
evidence on the maturation transition rates; nothing here supplies any. What is established is that
the state is *wired to something*, not that the scale factors are the ones talin has.

### 4.5 The stochastic transitions, and why they queue

Every transition — binding an actin side, engaging or releasing integrins, advancing maturation —
is **proposed** into a candidate buffer and applied only by `commit()`. `rollback()` discards them.
A bond formed inside a candidate step that is then rejected leaves the graph remembering a bond that
never formed, and the graph is then stiffer than its own accepted history says it is. No force check
can see this, because the force check is computed from the topology it is trying to check.

`focal_adhesion_accepted_binding_topology` is that committed graph, and `topology_epoch` is bumped by
any commit that changes it. Maturation is **irreversible and monotone**: a regression is refused, and
disassembly is expressed as releasing bonds rather than as un-maturing.

### 4a. Mutation testing of the controls

Three mutants applied to `aleph/vertical/focal_adhesion.py`, the suite run against each, the source
restored and verified byte-identically by SHA-256 and `git diff`. Recorded here rather than in prose,
because a control suite no mutant can break is decoration.

| # | Mutation | Tests killed | Survived? |
|---|---|---|---|
| M1 | `series_stiffness_pn_per_um`: `k_a k_l/(k_a+k_l)` → `k_a + k_l` (series becomes parallel) | **8** | no |
| M2 | `series_joint_energy_and_forces`: `engaged = actin_on & ligand_on` → all-true (the gate removed) | **3** | no |
| M3 | `Clutch.actin_side_stiffness_pn_per_um`: maturation scale forced to `1.0` (maturation unbound) | **3** | no |

M1 is the headline defect written as a one-token edit. Its eight kills come entirely from the closed
forms; **not one finite-difference control noticed it**, which is the measured version of §7's
argument — the parallel law is a perfectly conservative energy, so `F = −grad E` holds exactly for
the wrong stiffness and a gradient control and an oracle are not substitutes for each other.

**The exercise found two real gaps, which is the reason to run it.** Both mutants that killed only
two tests killed them one call *below* the API a consumer will use:

- **M2** was caught only by controls that call the free function directly. Nothing exercised the gate
  through `FocalAdhesionClutchGraph.joint_energy_and_forces`, which is the path every consumer takes,
  and a gate tested only one level below the API is a gate a caller can be routed around.
  `::TestTheJointTransmitsAndTheTransmissionCloses::test_an_open_joint_delivers_exactly_zero_through_the_owner`
  was added, with an open clutch beside a closed one so it cannot be satisfied by a graph that
  returns zero for everything.
- **M3** was caught only by controls that read a *stiffness*. Nothing asserted that maturation
  reaches the force actually delivered to the two borrowed endpoints — and a state that reaches a
  diagnostic but not the load is still a state no load reads.
  `::TestOccupancyAndMaturationAreBound::test_maturation_reaches_the_load_delivered_to_the_endpoints`
  was added and checks the delivered load against the closed form at two maturation states.

The maturation control also kills M1, which is how a control written for one reason earns its keep.
The open-joint control does not, and that is correct: M1 changes what a *closed* joint reports and
leaves an open one at zero either way.

Source restored and verified **byte-identical**: SHA-256
`dc940e5ad8076376d4786db1e9b27bc5c582189969559639c55d55a4a90afe11` before and after all three
mutations. The table above is the re-run against the final file, not the first pass — the exercise
was repeated after the two new controls landed so the counts describe the suite that ships.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| endpoint position (borrowed, never stored) | µm | 1e-6 m | finite |
| `bond_stiffness_pn_per_um`, `anchor_stiffness_pn_per_um` | pN/µm | 1e-6 N/m | `> 0` |
| `rest_gap_um` | µm | 1e-6 m | `> 0` |
| occupancy `n` | dimensionless count | — | integer in `[0, bond_capacity]` |
| maturation scale `m` | dimensionless | — | `> 0`, strictly increasing in maturation |
| joint load `T` | pN | 1e-12 N | finite |
| energy | pN·µm | 1e-18 J | `>= 0` |

Singular and boundary cases, each with the behaviour Aleph requires:

- **Coincident endpoints** (`d = 0`): refused with `ValueError`. The joint direction is genuinely
  undefined and a silent epsilon returns a finite force whose direction came from round-off.
- **`rest_gap_um <= 0`**: refused. A joint with no rest gap has no unloaded state.
- **Occupancy outside `[0, bond_capacity]`**: refused with `FocalAdhesionStateError`.
- **Occupancy `0` with ligand side `BOUND`, or occupancy `> 0` with ligand side `UNBOUND`**: refused.
  The two are one physical fact recorded twice and they may not disagree.
- **A disengaged clutch**: energy and both forces are **exactly** `0.0` — not small, identically
  zero on the whole disengaged branch, which is why the control asserts it with `==`.
- **Maturation regression** (`MATURE → NASCENT`): refused with `FocalAdhesionStateError`.
- **A duplicate clutch id**: refused with `FocalAdhesionOwnershipError`.
- **A nascent-only role published on a matured clutch**: refused with `FocalAdhesionRoleError`.
- **A handle used after a committed topology change**: refused with
  `FocalAdhesionStaleHandleError`.
- **A partial dispatch of `alpha2beta1_collagen_series`**: refused with `FocalAdhesionSeriesError`.
- **Empty clutch set**: returns exactly `(0.0, zeros, zeros)`, not a crash.
- **`k_a + k_l == 0`**: refused; the series combination is undefined there rather than zero.

Invariants, each with the test that asserts it, all in
`tests/vertical/test_focal_adhesion_controls.py`:

- **I1.** `F = −grad E` for the series law, order 2, in both endpoints separately —
  `::TestTheSeriesLawIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order`
- **I2.** One load, additive extensions —
  `::TestTheSeriesLawAgainstItsClosedForm::test_one_load_and_additive_extensions`
- **I3.** `k_series = k_a k_l/(k_a+k_l)` exactly, and strictly softer than either member —
  `::TestTheSeriesLawAgainstItsClosedForm::test_the_effective_stiffness_is_the_closed_form`
- **I4.** Force and moment closure through the joint —
  `::TestTheJointTransmitsAndTheTransmissionCloses::test_the_two_loads_are_exactly_equal_and_opposite`,
  `::test_the_joint_carries_no_net_moment`
- **I5.** Occupancy-grouping invariance —
  `::TestOccupancyAndMaturationAreBound::test_splitting_one_clutch_into_two_changes_no_load`
- **I6.** Maturation is bound to the mechanics —
  `::TestOccupancyAndMaturationAreBound::test_advancing_maturation_changes_the_load_it_carries`
- **I7.** Purity: no function mutates its arguments (`np.array_equal`, not `allclose`) —
  `::TestTheLawIsAPureFunction::test_no_call_mutates_its_arguments`
- **I8.** Rollback restores bit-identically —
  `::TestOwnershipAndTheAcceptedStep::test_rollback_restores_every_block_bit_identically`
- **I9.** Declared state keys equal the registered contract —
  `::TestTheDeclarationsAgreeWithTheCode::test_the_owned_state_keys_are_exactly_the_registered_contract`
- **I10.** This owner's endpoint is deliberately **not** a `CoupledSiteEndpoint` —
  `::TestTheEndpointShape::test_the_clutch_endpoint_is_not_a_coupled_site_endpoint`

## 5a. The endpoint shape — the open question, answered here and escalated

`docs/design/ENDPOINT_CONTRACT.md` §5 and `HANDOFF.md` §C-⑤ ask whether `focal_adhesion` publishes a
`CoupledSiteEndpoint` at all.

**It cannot, and this module does not.** `CoupledSiteEndpoint` requires `positions_um` — an `(N, 3)`
block of µm — and `scatter_forces_pn`. A geometry-less clutch has neither. Any position it returned
would be invented, and `as_site_positions` would accept the invention without complaint, producing an
endpoint that imports cleanly and reports a coordinate no registered state block contains. That is
the exact failure mode `PUBLISHABLE_ENDPOINT_ROLES` exists to prevent, one protocol lower down.

This module therefore publishes `ClutchEndpoint`, which carries **binding state and compliance**:
`is_engaged`, `series_stiffness_pn_per_um`, `rest_gap_um`, `issued_at_epoch`, and
`propose_accepted_load_pn` for the load-dependent unbinding law. It has no `positions_um` and no
`scatter_forces_pn`, and a control asserts `isinstance(handle, CoupledSiteEndpoint) is False`.

This is not a refusal to be wired. The series joint needs three things — the actin-side material
point (a `CoupledSiteEndpoint` from `sf_arc`), the ligand-side material point (a
`CoupledSiteEndpoint` from `ecm`), and the compliance and gate in between (this `ClutchEndpoint`) —
and force lands on the two owners that have degrees of freedom. **Nothing is withheld from the
connector that the connector could use**; what is withheld is a position that does not exist.

Written up, with the options and the reversal cost, in
`docs/decisions/PROPOSAL-focal-adhesion-endpoint-shape.md` (`Status: AGENT-PROPOSED`).

## 5b. Two findings about the registry, reported and not fixed

Neither is this lane's to change, and `CLAUDE.md` §2 rule 5 says a contradiction is a finding.

1. **The composite group is one joint with three owners, and the schema has room for two.**
   `fa_actin_anchor` (`sf_arc`↔`focal_adhesion`) and `integrin_collagen_clutch`
   (`focal_adhesion`↔`ecm`) share `composite_group="alpha2beta1_collagen_series"` and must be
   dispatched as one joint. That joint spans **three** owners and delivers force to **two** of them.
   A `ConnectorContract` has exactly two endpoints, so the object that must be dispatched cannot be
   named by any single row of the registry — it exists only as the group. `CONNECTOR_PROTOCOL`'s
   `evaluate_sites → SitePairForces` likewise assumes two force blocks per *connector*, where here
   there are two force blocks per *group*.
2. **`lamellipodium_nascent_fa` and `filopodium_nascent_fa` reach the same actin side with
   `composite_group=None`.** Their prose says force reaches the substrate only through the
   adhesion's own series path, but the field that would enforce it is unset, so the registry permits
   for a nascent adhesion exactly the independent dispatch it forbids for a mature one. This module
   closes the hole on its own side — the series law is the only route from an actin-side reference
   to a ligand-side reference, whichever connector supplied the reference — but the registry
   asymmetry is real and is not mine to edit.

## 6. Source evidence class and known retractions

No claim is made about the reference implementation's adhesion model, because none of it was read.
Its evidence class is **unknown here**. Where I looked: nowhere inside it.

No claim is carried forward from `_archive/focal_adhesion_unverified_2026-07-30/`. In particular the
archived lane's reported traction peak at a substrate stiffness of 1000 pN/µm is **not** reproduced,
not cited, and not assumed; the archive README marks it as a lane's report rather than a measurement,
and this lane did not re-run it. **This module makes no traction claim of any kind** — the registered
contract lists traction magnitude under `unsupported_claims` because the binding rate and the
load-dependent unbinding law are both unsourced, and nothing here changes that.

`citation_status` for the `focal_adhesion` contract remains `UNSOURCED`. `default_clutch_card()`
returns placeholder magnitudes and says so in its own docstring.

## 7. Independent oracle or derivation

Three, none of which is the reference implementation and none of which is the module's own
vocabulary:

1. **Central finite differences** of the module's own energy, in each endpoint separately. An oracle
   for the gradient claim that is independent of every constant in the law.
2. **The series closed form** `k_series = k_a k_l/(k_a+k_l)`, derived in §4.2 from the two series
   statements alone. This catches what the finite difference cannot: the parallel law
   `k_a + k_l` is *also* a perfectly conservative energy, so `F = −grad E` holds exactly for the
   wrong stiffness. Mutant M1 is precisely that, and it is invisible to every gradient check.
3. **The additive-extension identity** `e_a + e_l = e` with `T = k_a e_a = k_l e_l`, checked as two
   independent equalities rather than as the combination they imply. A single check of `T = k e`
   cannot distinguish a right stiffness reached through wrong extensions.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_focal_adhesion_controls.py::TestTheSeriesLawIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order` | `energy > 0` first, then max&#124;F+∇E&#124;/max&#124;F&#124; at steps 4e-4/2e-4/1e-4 µm, observed order > 1.6, per endpoint |
| Positive | `tests/vertical/test_focal_adhesion_controls.py::TestTheSeriesLawAgainstItsClosedForm::test_the_effective_stiffness_is_the_closed_form` | Measured `T/e` equals `k_a k_l/(k_a+k_l)` to 1e-12 relative, and is strictly below `min(k_a, k_l)` |
| Positive | `tests/vertical/test_focal_adhesion_controls.py::TestTheSeriesLawAgainstItsClosedForm::test_one_load_and_additive_extensions` | The two halves carry bit-identical load and their extensions sum to the joint extension |

**Measured numbers are in §8a**, filled from a run in this session. No number in this entry is quoted
from memory.

### 8a. Measured, this session, on this tree

Interpreter `/Users/sw1/miniconda3/envs/aleph/bin/python`. Every number below was produced by a run
in this session; none is quoted from memory. Geometry: three clutches, endpoints of order 1 µm, joint
axis `(2, −1, 2)/3` so that no force component is zero by construction, gap `0.30 µm`, rest gap
`0.05 µm`, extension `0.25 µm`.

Finite-difference agreement of the series law, steps 4e-4 / 2e-4 / 1e-4 µm:

| Perturbed endpoint | Energy (pN·µm) | Force scale (pN) | relative err at 1e-4 | observed order |
|---|---|---|---|---|
| anchor | 7.031250 | 12.5000 | 6.173e-09 | **2.000** |
| ligand | 7.031250 | 12.5000 | 6.173e-09 | **2.000** |

Series closed form, `k_a = 100`, `k_l = 300` pN/µm, one clutch:

| Quantity | Predicted | Measured | relative |
|---|---|---|---|
| `k_series` | 75.0 | 75.000000000000 | 1.895e-16 |
| `T_a − T_l` | 0.0 | `0.0` | **exact** |
| `e_a + e_l − e` | 0.0 | −5.551e-17 | 2.220e-16 |
| independent-branch curvature | 400.0 | 400.000000 | 0.000e+00 |
| overstatement `(k_a+k_l)²/(k_a k_l)` | 5.333333 | 5.333333 | 0.000e+00 |

Closure through the joint, four clutches, residual ÷ the **constituent** scale (`Σ|f_i|` for force,
`Σ|r_i||f_i|` for moment — never the resultant):

| Branch | force residual | moment residual |
|---|---|---|
| series (healthy) | **0.000e+00** (bit-exact) | 7.869e-18 |
| independent (broken) | **5.000e-01** | **9.479e-02** |

The broken branch's force residual of exactly one half is the closed form
`|T_a − T_l| / (|T_a| + |T_l|) = |k_a − k_l|/(k_a + k_l) = 200/400`.

**A correction this lane made to its own entry.** A first draft of §4.2 argued that the moment closes
in both branches, because both loads are parallel to the joint axis, and concluded that a moment
check is blind to this defect. That is wrong, and the control that was written to assert it failed.
The net moment is `(T_a x_a − T_l x_l) × u`, which collapses to `T (x_a − x_l) × u = 0` **only when
the two loads are equal**. So moment closure is a consequence of the series law rather than of the
geometry, and it is a second independent witness of the defect rather than a check that cannot see
it. §4.2 is corrected, and
`::TestTheJointTransmitsAndTheTransmissionCloses::test_the_broken_branch_fails_force_closure_against_the_constituent_scale`
now asserts both residuals in both branches. This is exactly the case PLAN §0.2.5 anticipates —
writing the code made a section of the entry wrong, and the entry moved.

## 9. Deliberately failing negative control

Every break below is a flag on the **shipped** function or owner, default `False`, so the control
drives the real code path rather than a hand-edited copy.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `::TestTheIndependentSpringEvaluationFailsTheSeriesControl::test_the_independent_evaluation_reports_a_joint_stiffer_than_either_member` | With `evaluate_halves_independently=True` the energy curvature is `k_a + k_l = 400` pN/µm, stiffer than either member and `5.3333…×` the series value; with the flag off it is exactly `75` |
| Negative (must fail) | `::TestTheIndependentSpringEvaluationFailsTheSeriesControl::test_the_independent_evaluation_makes_the_plaque_a_hidden_force_source` | Broken branch: net force on the isolated joint is `(k_a − k_l) e ≠ 0`. Healthy branch: the two forces sum to **exactly** `0.0`, asserted with `np.array_equal` |
| Negative (must fail) | `::TestTheIndependentSpringEvaluationFailsTheSeriesControl::test_an_open_clutch_still_reports_load_when_the_halves_are_split` | Broken branch: a disengaged clutch's actin half still reports `k_a e`. Healthy branch: **exactly** `0.0` |
| Negative (must fail) | `::TestTheNegativeControlsFail::test_a_flipped_joint_sign_is_caught_by_the_gradient_test` | `flip_joint_sign=True` drives the FD disagreement above 0.5 relative; measured 2.0 |
| Negative (must fail) | `::TestTheNegativeControlsFail::test_an_open_clutch_that_carries_load_is_visible_in_the_load_it_reports` | `allow_open_clutch_to_carry_load=True` gives a disengaged clutch a finite load; with the gate in place the same configuration is **exactly** `0.0` |
| Negative (must fail) | `::TestOccupancyAndMaturationAreBound::test_an_unbound_maturation_law_is_caught_by_the_guard` | `omit_maturation_law=True` makes advancing maturation change nothing at all, and `assert_maturation_is_bound` raises |

Each is asserted in **both** directions on purpose. Asserting only that the broken branch misbehaves
passes for a function that misbehaves unconditionally; asserting only that the healthy branch is zero
passes for a function that returns zero always. The `PLAN.md` §6.1 defect — a cold cache read making a
deliberately broken strut look innocent — is the reason this is a rule rather than a preference.

## 10. Numerical and precision envelope

float64 throughout; no reduced-precision path exists.

Finite-difference steps 4e-4 / 2e-4 / 1e-4 µm on coordinates of order 1 µm. The floor for a central
difference is `eps^(1/3) ≈ 6e-6` relative; the smallest step is 17× above it. Below the floor,
cancellation in the energy evaluation dominates, the apparent order goes negative, and a correct
gradient is indistinguishable from a wrong one.

Tolerances and why they are not looser:

- The series closed form is asserted at **1e-12 relative** against a measured `0.0`. It can be this
  tight because `k_a k_l/(k_a+k_l)` is three floating-point operations, not a discretisation; mutant
  M1 perturbs it by a factor of 5.3, so there is no tolerance at which M1 survives.
- Force closure in the healthy branch is asserted with `np.array_equal` against exact zeros, not with
  a tolerance. The forces are constructed as `+pull` and `−pull` from one array, so their sum is
  bit-exact zero; accepting a tolerance here would accept an implementation that computed the two
  sides separately, which is the defect.
- The disengaged branch is asserted with `==` for the same reason: the gate returns identical zero on
  the whole branch rather than a small number near it.
- The FD relative agreement is asserted at 1e-7 against a measured 4.0e-09.

Outside the envelope: every singular case in §5 refuses with a typed exception rather than degrading
silently. There is no input range over which this module returns a plausible number it cannot stand
behind.

## 11. Production-backend residency and transfer

Host, numpy, float64. `accumulate` uses `ctx.backend` only for the reductions it reports; the joint
law is `O(C)` in the clutch count with no host round-trip per step and no device allocation. No GPU
work of any kind was run by this lane, on either host, and no authorization was sought or held.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered this module, because no reference file was read, and
no archived-module prose entered it either. Asserted, not asserted-in-prose:

- `::TestTheDeclarationsAgreeWithTheCode::test_no_reference_project_vocabulary_survives` — the module
  text carries none of `ffn_cellsim`, `ffn_sim`, `dcm_contact`, `surface_body`.
- `::TestTheDeclarationsAgreeWithTheCode::test_the_module_imports_nothing_from_validation` — the
  `aleph/** → validation/**` firewall.
- `::TestTheDeclarationsAgreeWithTheCode::test_the_ledger_entry_the_module_cites_exists` — the file
  you are reading. `ALEPH-PORT-2901` cited six control tests in two files that never existed; a
  docstring citing a ledger that does not exist is the same failure with the arrow reversed.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 80 controls pass (`tests/vertical/test_focal_adhesion_controls.py`, 2026-07-30, 0.22 s). Three mutants applied, killing 8 / 3 / 3 tests; none survived; source restored byte-identically by SHA-256. Status stays `PROPOSED` because zero connectors are wired to this owner and the endpoint-shape question in §5a is a PI decision that has not been taken. |
| Reviewer | Agent-proposed (lane L32). Unratified. No PI review. |
| Rollback | Delete `aleph/vertical/focal_adhesion.py`, `tests/vertical/test_focal_adhesion_controls.py`, and `docs/decisions/PROPOSAL-focal-adhesion-endpoint-shape.md`, then set this entry `REJECTED`. Nothing else breaks: no module in the tree imports the owner, `aleph/vertical/__init__.py` does not export it, and no connector reaches it. |

## 14. Honest limits — what this entry does NOT establish

**Wired connector count contributed by this module: zero.** Four connectors name `focal_adhesion` —
`fa_actin_anchor`, `integrin_collagen_clutch`, `lamellipodium_nascent_fa`, `filopodium_nascent_fa` —
and none of them is wired by this lane. This module makes their `focal_adhesion` endpoint *exist*;
whether a connector can be written against a `ClutchEndpoint` under the current
`CONNECTOR_PROTOCOL` is §5b finding 1 and is unresolved.

**No stochastic transition rates exist.** The registered contract owes evidence on the
integrin/collagen binding rate, its load-dependent unbinding law, the actin-side anchoring kinetics
and the maturation transition rates. **This module supplies none of them and deliberately does not
default them.** `propose_*` takes the transition to make as an argument; it does not decide when to
make one. So the graph is stochastic in *structure* — transitions are events committed only on an
accepted step — and there is no kinetic model behind them. A consumer that needs a rate must refuse,
which the registered `approximation` explicitly requires for the ligand side. **`UNVERIFIED`:
nothing here is evidence about adhesion lifetimes, bound fractions, or turnover.**

**No traction number, and no adhesion count.** Both are in the contract's `unsupported_claims`. The
module has no notion of adhesion size, area, shape or spatial growth, because a geometry-less
representation does not have them; a consumer asking for one gets a refusal rather than a value.

**The maturation scale factors are invented.** `m(NASCENT) = 1.0`, `m(MATURING) = 2.0`,
`m(MATURE) = 4.0` are placeholders chosen to make the "maturation is bound" control sharp. They are
labelled `UNSOURCED` in the code. What is established is that maturation is wired to the compliance,
**not** that talin reinforcement has this magnitude or this ordering-by-factor-of-two.

**No transaction has ever run.** `snapshot`/`commit`/`rollback` and `accumulate` are exercised by
this file's controls in isolation. This owner has never been stepped inside a world, is not in
`aleph/vertical/assembly.py`, and no acceptance predicate has ever read its receipt. The
rejected-step integrity control drives `rollback()` directly rather than through
`aleph.runtime.transaction`.

**The independent-spring defect is established as an arithmetic fact, not as a simulation result.**
The controls show that the parallel evaluation gives `k_a + k_l` where the series law gives
`k_a k_l/(k_a+k_l)`, and that the difference injects a net force. They do **not** show what that does
to a cell — no relaxation was run, no traction was measured, and evidence rung is `ANALYTIC_ORACLE`
with quantitative status `BLOCKED`.

**The reference implementation is unaudited by this lane**, and so is the archived Aleph attempt
beyond its README. Whether either agrees with this derivation is unknown, and this entry deliberately
does not guess.

**Not run:** any GPU job (zero); the full suite. Only `tests/vertical/test_focal_adhesion_controls.py`
and `tests/ports/` were exercised — other sessions are writing this tree concurrently, and `CLAUDE.md`
§1 says to report a foreign breakage rather than conflate it with one's own. Two pre-existing failures
in `tests/ports/` are reported in this lane's summary and were **not** touched.
