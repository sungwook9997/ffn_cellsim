# ALEPH-PORT-3303 — the four adhesion connectors, dispatched as series joints and never as halves

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3303` |
| Lane | lead session `f70de564` — adhesion connector lane |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Builds on | `ports/ledger/ALEPH-PORT-3201-focal-adhesion-stochastic-clutch-graph.md` |
| Decision it depends on | `docs/decisions/PROPOSAL-focal-adhesion-endpoint-shape.md` (Option C, `AGENT-PROPOSED`) |
| Decision it raises | `docs/decisions/PROPOSAL-a-composite-group-is-one-connector-object.md` |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
# new module aleph.vertical.connectors_adhesion
AdhesionCouplingError
SeriesJointConnector
LIGAND_SIDE_OWNER
NASCENT_ROLE_FOR_ACTIN_OWNER
SERIES_GROUP_ROWS
build_alpha2beta1_collagen_series      # ONE object, TWO registry rows
build_filopodium_nascent_fa
build_lamellipodium_nascent_fa
```

**Three builders for four connectors, and that is the finding rather than an omission.**
`fa_actin_anchor` and `integrin_collagen_clutch` share `composite_group="alpha2beta1_collagen_series"`
and are the two halves of one series joint. There is no builder for either half alone, because a half
is not a dispatchable object — see §4. `ALEPH-PORT-3101` §1 records the same discipline in the
opposite direction: an absent name is a visible gap and a stub is an invisible one.

Nothing in `aleph/vertical/wiring.py`, `aleph/vertical/focal_adhesion.py` or `aleph/state/**` is
modified by this entry. The `Binding` rows it earns are **reported** to the lead session and added by
the single writer of `wiring.py`, in the same commit that wires them.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | **none named and none read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under the reference repository was opened for this work. |
| Working tree == commit? | not applicable — nothing was read |

The archived `_archive/focal_adhesion_unverified_2026-07-30/` tree in **this** repository was also not
opened; its own `README.md` says not to copy from it, and `ALEPH-PORT-2901` is `REJECTED`.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** Two linear springs in series have compliances that add:

```
1/k = 1/k_a + 1/k_l   =>   k = k_a k_l / (k_a + k_l)
```

That is secondary-school statics, it is already derived and gradient-checked inside
`aleph/vertical/focal_adhesion.py` under `ALEPH-PORT-3201`, and per `ports/TEMPLATE.md` §3 a
convention any competent author reaches independently is a reason to write clean-room. What this
entry adds is **not a law**. It is the interface: the object a schedule can evaluate, snapshot and
roll back, and the routing that decides which owner each of the two force blocks lands on.

## 4. Physical or mathematical law represented

**The claim in one line: the dispatchable object is the series joint, and the series joint is named
by two registry rows, so a connector written against one row cannot exist.**

### 4a. The joint has exactly two force blocks, and neither of them is on `focal_adhesion`

The load path is

```
actin material point  ->  talin/integrin state  ->  collagen ligand material point
   (sf_arc |                (focal_adhesion,           (ecm, has a position)
    lamellipodium |          NO position, NO
    filopodium)              degrees of freedom)
```

`focal_adhesion` is registered as a **geometry-less** clutch graph: no adhesion mesh, no adhesion
particle, no adhesion volume. A body with no degrees of freedom holds no force, so whatever load
enters the actin side leaves the ligand side. The joint therefore delivers exactly two blocks of
site force — one to the actin owner, one to `ecm` — which is precisely the shape of
`aleph.vertical.connectors.SitePairForces`.

With `d = |x_l - x_a|`, `u_al = (x_l - x_a)/d`, `u_la = -u_al` computed independently, `e = d - d0`:

```
k = k_a k_l / (k_a + k_l)      (0.0 exactly when the clutch is open)
T = k e        E = (k/2) e^2
f_actin  = T u_al              f_ligand = T u_la
```

`k` is read from `ClutchEndpoint.series_stiffness_pn_per_um`, the one number that owner publishes as
a load path. `half_stiffnesses_pn_per_um` is documented there as a diagnostic and **is never
multiplied by an extension** in the shipped path; the only code here that touches it is the
deliberate break in §9.

### 4b. Why `CONNECTOR_PROTOCOL` fits the group and not the row

`CONNECTOR_PROTOCOL = ("evaluate_sites", "accumulate")` and `evaluate_sites` returns
`SitePairForces` — two force blocks per **connector**. `PROPOSAL-focal-adhesion-endpoint-shape.md`
recorded the arity problem: the joint spans three owners and is named by two rows, so "the object
that must be dispatched is the *group*, not either row".

That statement is right, and the conclusion drawn from it here is narrower than "the protocol cannot
express a composite group". The protocol is satisfied **by the group object**:

| Question | Answer |
|---|---|
| How many force blocks does the joint produce? | Two. `SitePairForces` fits without being widened. |
| Which owners do they land on? | The actin owner and `ecm` — the two with degrees of freedom. |
| How many registry rows name the joint? | Two (`fa_actin_anchor`, `integrin_collagen_clutch`). |
| Can one row be evaluated alone? | **No.** Its `endpoint_b` is `focal_adhesion`, which has no position to evaluate at and no array to receive a reaction. |

So the one thing that genuinely does not hold is the assumption **one `Binding` row ↔ one
dispatchable object**. Two rows here bind to one object. That is a claim about the wiring table and
not about the physics, and it is raised as
`docs/decisions/PROPOSAL-a-composite-group-is-one-connector-object.md` rather than decided here.

### 4c. What an independent dispatch would do, in arithmetic

Evaluating each row against its own kinematics gives each half the load its own extension implies.
The middle has no position, so both halves see the same `e`:

1. `T_a = k_a e` and `T_l = k_l e` differ, so the joint applies unequal loads to its two ends and
   injects `(k_a - k_l) e` into the world. The plaque becomes a hidden force source.
2. The energy curvature becomes `k_a + k_l`, the **parallel** combination, stiffer than either member
   by `(k_a + k_l)^2 / (k_a k_l) >= 4`. Stiffer joints read as higher traction, so the error is
   one-sided and upward.
3. Each half applies its own gate, so a clutch with its ligand side released still has its actin half
   reporting load — a traction reading with nothing on the other end of it.

All three are checked against closed forms in §8, and all three are driven by a flag on the
**shipped** class rather than by a hand-edited copy.

### 4d. The registry asymmetry, and the position this module takes

`lamellipodium_nascent_fa` and `filopodium_nascent_fa` reach the same actin-side state with
`composite_group=None`. Their prose says force reaches the substrate only through the adhesion's own
series path, but the field that would enforce it is unset: **the registry permits for a nascent
adhesion exactly the independent dispatch it forbids for a mature one.** Reported, not edited —
`PROPOSAL-focal-adhesion-endpoint-shape.md` §"A second finding" raised it and `CLAUDE.md` §2 rule 5
says which rules a module may move.

The position taken here, stated because the design had to take one: **`SeriesJointConnector` has no
branch that evaluates a half.** A nascent connector is built from the same class as the mature group,
so the series discipline the registry asks of the composite group is applied to the nascent rows
whether the registry demands it or not. The only way to obtain a half evaluation is the deliberate
break, which is off by default and exists so the control can catch it. If the PI decides the nascent
rows should be independent springs, that is a one-line change and a deleted control; if the PI
decides the registry should carry `composite_group` on all four rows, nothing here changes at all.

### 4e. What "the connector owns the binding state" means here

Both composite rows carry `commit_on_accept=True` and mechanism prose saying the anchoring state is
"owned and committed by the connector", while the census gives `focal_adhesion` the state keys
`focal_adhesion_actin_side_state` and `focal_adhesion_ligand_side_state`. Two registered artefacts,
one fact, and they do not agree about who holds it. **Reported, not resolved** — `CLAUDE.md` §2
rule 5.

What this module does in the meantime is the reading that contradicts neither: the array stays with
the owner, and the connector reaches it only through the published channel. After delivering load the
connector calls `ClutchEndpoint.propose_accepted_load_pn`, which **queues** the load into the owner's
candidate buffer; the owner applies it on `commit()` and discards it on `rollback()`. No transition is
proposed by this module, because **no kinetic rate law exists anywhere in this project** — the
registered contract owes all four and `ALEPH-PORT-3201` §14 says so. Absent, not implied.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI | Domain |
|---|---|---|---|
| anchor / ligand position | µm | 1e-6 m | finite |
| series stiffness `k` | pN/µm | 1e-6 N/m | `>= 0`, exactly `0.0` when open |
| rest gap `d0` | µm | 1e-6 m | `> 0` |
| joint tension `T` | pN | 1e-12 N | finite |
| joint energy `E` | pN·µm | 1e-18 J | `>= 0` |

Singular and boundary cases, each with the behaviour required:

- **Coincident anchor and ligand.** Refused. The joint axis is genuinely undefined and a regularised
  denominator would return a finite force whose direction came from round-off in the last
  subtraction.
- **Open clutch** (either side unbound). Contributes **exactly** `0.0` to energy and to both force
  blocks — the number, not a small number — because the branch is identically zero, not near zero.
- **Zero-length position block.** Refused: a connector with no sites transmits nothing, and a
  zero-site endpoint is indistinguishable downstream from a broken one.
- **Mismatched block sizes** (anchors, ligands, clutch handles). Refused.
- **Non-finite position.** Refused.
- **A clutch handle held across a committed binding change.** Refused by
  `FocalAdhesionStaleHandleError`, which the owner already raises; this connector holds handles and
  therefore inherits the guard rather than defeating it by caching numbers.
- **A composite group dispatched with one row missing.** Refused by
  `focal_adhesion.assert_series_group_is_dispatched_whole`, called at construction. The guard is
  imported, not reimplemented.

Invariants, each with the control that asserts it:

- **I1.** The connector's site forces equal the owner's own `joint_energy_and_forces` on the same
  configuration, to round-off. Two independently written expressions, one law.
- **I2.** `sum f_actin + sum f_ligand` closes against the **constituent** scale `sum|f_i|`, never
  against the resultant.
- **I3.** `F = -grad E` against a central finite difference, observed order `> 1.6`, with `E > 0`
  asserted first — a term returning zero energy and zero force agrees with its own finite difference
  perfectly.
- **I4.** The reported joint stiffness equals the closed form `k_a k_l/(k_a+k_l)` and is strictly
  **less** than either member.
- **I5.** An open clutch transmits exactly `0.0`, in both force blocks and in the energy.
- **I6.** Every row this object serves is a declared connector whose registered endpoints match what
  the object claims; and the two force owners are **not** the registered endpoint pair, which is
  asserted rather than left as a surprise.
- **I7.** The load delivered is queued, not applied: visible in the owner's pending buffer, applied
  on `commit`, discarded on `rollback`.

## 6. Source evidence class and known retractions

No claim is made about the reference implementation, and nothing there was read. Its own record, as
quoted in `PLAN.md` §7, is that `SOLVE_COUPLED` has a declaration and zero production tasks.

Within this repository: `ALEPH-PORT-2901` (`alpha2beta1-collagen-series-joint`) is **`REJECTED`** and
its code is archived; nothing here derives from it. `ALEPH-PORT-3201` is `PROPOSED` and supplies the
series law this entry consumes. The `focal_adhesion` contract's `citation_status` is `UNSOURCED`,
every magnitude on `ClutchCard` is `UNVERIFIED`, and nothing in this entry changes either. **No
traction, adhesion-lifetime or bound-fraction conclusion may be drawn from anything here.**

## 7. Independent oracle or derivation

1. **The series closed form.** `k = k_a k_l/(k_a+k_l)`, computed in the control from the two half
   stiffnesses read off the endpoint's diagnostic accessor, against the `k` the connector actually
   used. The parallel combination differs from it by a factor `>= 4`, so no tolerance separates them
   and a finite-difference check cannot tell them apart — the parallel law is also perfectly
   conservative and `F = -grad E` holds exactly for the wrong stiffness. Only the closed form catches
   it.
2. **The owner's own gradient**, `FocalAdhesionClutchGraph.joint_energy_and_forces`, evaluated on the
   same configuration by a different code path (population arrays and half stiffnesses) than the
   connector uses (per-clutch published series stiffness).
3. **Central finite differences** of the connector's own energy, with a measured convergence order.
4. **Newton's third law**, from two independently normalised direction vectors. `force_b := -force_a`
   appears nowhere; writing it would turn the closure check into a tautology.

## 8. Positive control

Measured on this tree, this session, with `/Users/sw1/miniconda3/envs/aleph/bin/python`,
macOS/darwin, float64. Every number below came from a run in this session.

Fixture, so every number below can be recomputed: `k_bond = 1.0e2`, `k_anchor = 3.0e2`,
`d0 = 5.0e-2` µm, occupancy `2`, maturation `MATURING` (scale `2.0`). Hence
`k_a = 2·300·2 = 1200`, `k_l = 2·100 = 200`, `k_series = 1200·200/1400 = 171.42857142857142` pN/µm.
The joint offset is `(0.12, 0.09, 0.14)` µm — **deliberately not axis-aligned**, so no coordinate
makes the joint energy exactly quadratic and the finite-difference control measures a truncation
order rather than round-off.

| Control | Location | Asserts | Measured |
|---|---|---|---|
| Positive | `tests/vertical/test_connectors_adhesion.py::TestItAgreesWithTheOwnersOwnLaw::test_the_connector_site_force_equals_the_owner_s_own_joint_gradient` | I1 | max abs difference **0.0** pN on both blocks; energies identical at **6.19244112119575** pN·µm |
| Positive | `tests/vertical/test_connectors_adhesion.py::TestItAgreesWithTheOwnersOwnLaw::test_the_reported_stiffness_is_the_series_closed_form_and_softer_than_either_member` | I4 | relative error **0.0**; `k = 171.42857142857142` pN/µm against `k_a = 1200.0`, `k_l = 200.0` |
| Positive | `tests/vertical/test_connectors_adhesion.py::TestTheGradientIsTheForce::test_the_joint_force_is_minus_the_energy_gradient` | I3 | observed orders **2.0000, 2.0000**; errors 6.268e-06 → 1.567e-06 → 3.917e-07 pN, at `E = 4.128294080797167` pN·µm |
| Positive | `tests/vertical/test_connectors_adhesion.py::TestClosureAndInvariance::test_the_pair_closes_against_the_constituent_scale` | I2 | residual **0.0** pN against a constituent scale of **212.82** pN — see §10, this is exact and the reason matters |
| Positive | `tests/vertical/test_connectors_adhesion.py::TestClosureAndInvariance::test_the_joint_is_translation_invariant_and_deterministic` | I2 | determinism `np.array_equal` **True**; translation drift/scale **4.50e-15** under a `(7, −3, 11)` µm shift |
| Positive | `tests/vertical/test_connectors_adhesion.py::TestClosureAndInvariance::test_the_two_force_blocks_are_accumulated_independently` | the source discipline behind I2 | two normalisations, no negation |
| Positive | `tests/vertical/test_connectors_adhesion.py::TestTheOpenClutchCarriesExactlyZero::test_an_open_clutch_contributes_exactly_zero_to_both_blocks` | I5 | `== 0.0` exactly, both blocks and the energy |
| Positive | `tests/vertical/test_connectors_adhesion.py::TestTheRowsAgreeWithTheRegistry::test_every_row_is_declared_with_the_endpoints_this_object_claims` | I6 | 4 rows across 3 objects |
| Positive | `tests/vertical/test_connectors_adhesion.py::TestTheAcceptedStepContract::test_the_delivered_load_is_queued_and_becomes_real_only_on_commit` | I7 | pending → committed; and `::test_a_rejected_step_discards_the_delivered_load` for the other branch |

**Suite for this entry: 35 passed** (`tests/vertical/test_connectors_adhesion.py`, 0.38 s,
2026-07-31).

`tests/ports/` re-run alongside: **15 passed, 2 failed, and both failures are foreign.**
`test_named_controls_resolve_to_real_tests` fails on `ALEPH-PORT-3302` and `ALEPH-PORT-3304`, two
sibling lanes' entries that are mid-write; `test_index_is_not_stale` fails because
`ports/ledger/INDEX.md` predates `ALEPH-PORT-3301`, `-3302`, `-3303` and `-3304`. **INDEX.md was not
regenerated by this lane**: doing so would write three other sessions' in-progress entries into a
shared generated file, which is the failure `docs/ACTIVE_SESSIONS.md` records twice in one afternoon.
Reported to the lead session rather than repaired — `CLAUDE.md` §1.

The I1 result is the one worth reading twice. The difference is **exactly 0.0** and not merely small,
across two expressions written independently: the owner assembles population arrays of `k_a` and
`k_l` and combines them, while the connector reads one published `k` per clutch through an endpoint
handle. They agree bit-for-bit because the same three floating-point operations end up in the same
order, which is a fact about this fixture and not a guarantee; the control asserts `< 1e-12` so a
reordering does not turn a correct connector red.

## 9. Deliberately failing negative control

Each drives a flag on the **shipped** class, default `False`, so the control exercises the real code
path rather than a copy of it.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_connectors_adhesion.py::TestTheNegativeControlsFail::test_dispatching_the_two_rows_independently_breaks_the_newton_pair` | With `dispatch_rows_independently=True` the joint injects `(k_a − k_l) e` into the world; the healthy branch closes at `<= 1e-14` relative |
| Negative (must fail) | `tests/vertical/test_connectors_adhesion.py::TestTheNegativeControlsFail::test_dispatching_the_two_rows_independently_reports_the_parallel_stiffness` | Broken curvature is `k_a + k_l`, measured `(k_a+k_l)^2/(k_a k_l)` times the series value and `>= 4` |
| Negative (must fail) | `tests/vertical/test_connectors_adhesion.py::TestTheNegativeControlsFail::test_an_open_clutch_that_still_reports_load_is_caught` | With `transmit_through_an_open_clutch=True` a released clutch transmits; the healthy branch is **exactly** `0.0` |
| Negative (must fail) | `tests/vertical/test_connectors_adhesion.py::TestTheNegativeControlsFail::test_a_flipped_joint_sign_breaks_the_gradient_check` | `flip_joint_sign=True` makes a stretched clutch push, and the finite-difference control detects it |
| Negative (must fail) | `tests/vertical/test_connectors_adhesion.py::TestTheNegativeControlsFail::test_building_half_of_the_composite_group_is_refused` | Constructing an object for `fa_actin_anchor` without `integrin_collagen_clutch` raises `FocalAdhesionSeriesError` |

The third is the sharpest, and it is the one the registry singles out. A control that only checked
the closed clutch would pass a connector that returned zero for everything; asserting **both**
directions — the broken branch transmits, the healthy branch is exactly inert — is what makes the
zero mean something. `PLAN.md` §6.1 records a negative control that could not fail because a cold
cache returned `0.0` for the correct element and the broken one alike.

## 10. Numerical and precision envelope

float64 throughout; no reduced-precision path.

The open-clutch branch is asserted with `==` and not a tolerance, because `k` is exactly `0.0` on the
whole disengaged branch rather than small near it — `Clutch.series_stiffness_pn_per_um` returns the
literal `0.0` when `is_engaged` is `False`.

Closure (I2) is judged against `sum_i |f_i|`, never against the resultant. `PLAN.md` §2.5 records a
resultant-based tolerance rejecting 40,000 consecutive steps of a correct relaxation; a two-block
central pair cancels to round-off by construction, so a residual judged against the resultant gets
stricter the more correct the physics is.

The finite-difference control (I3) uses steps `(4e-4, 2e-4, 1e-4)` µm and asserts an observed order
`> 1.6` rather than `≈ 2`. Below about `1e-5` µm the cancellation error in the difference quotient
exceeds the truncation error and the order collapses; above about `1e-2` µm the quadratic energy's
higher derivatives are no longer negligible against the step. `E > 0` is asserted first.

**The measured closure residual is exactly `0.0`, and that is a limitation rather than a result.**
The two blocks are built from two independently normalised direction vectors — `ligand - anchor` and
`anchor - ligand`, each divided by its own length — and in IEEE-754 negation and division are both
exact on the sign, so the two directions come out bit-for-bit opposite. The numerical closure check
therefore **cannot fail for this element**, which means it is not the check that carries the claim.

Written down rather than left as a comfortable zero, because a residual of `0.0` reads like strong
evidence and here it is none. What actually enforces the discipline is
`test_the_two_force_blocks_are_accumulated_independently`, which reads the source: two
`np.linalg.norm` calls, and no `force_b = -force_a`. The mutation table below confirms it — the
mutant that assigns the negation is killed by that control and by nothing else.

The same is true of `aleph/vertical/connectors.py::site_pair_forces`, which normalises
`midpoint_a - midpoint_b` and `midpoint_b - midpoint_a` separately for the same reason and gets the
same exact cancellation. This module follows that precedent rather than `focal_adhesion`'s own
`series_joint_energy_and_forces`, which constructs `+pull` and `-pull` from one array and asserts
bit-exactness deliberately. Recorded so a reader meeting the two disciplines side by side does not
"repair" one into the other.

Translation invariance is asserted at `4.50e-15` relative and **not** with `==`. The invariance is
exact as a property of the law; the arithmetic that evaluates it is not, because forming
`ligand - anchor` after shifting both by `(7, −3, 11)` µm loses low-order bits to cancellation.
Determinism *is* asserted with `np.array_equal`, because the same input evaluated twice has no reason
to differ in any bit.

## 11. Production-backend residency and transfer

Host, numpy, float64. **No GPU work of any kind was run by this lane, and no authorization was sought
or held.** The joint is `C` independent three-vectors and one normalisation per site; it is a natural
device kernel and nothing here is on the critical path for that. No host round-trip is required per
step.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered, because no reference file was read. No provider
vocabulary, module path, gate name or absolute source path appears in
`aleph/vertical/connectors_adhesion.py`; `tests/ports/test_port_discipline.py` group 1 checks it.

One piece of **this project's own** prose is corrected rather than discarded.
`PROPOSAL-focal-adhesion-endpoint-shape.md` §"What this costs the connector lane" says "A connector
written against one row cannot honour the composite rule, because it can only see one half." That
remains true and is the reason this module has no per-row builder. What the same section leaves open
— whether the protocol can express the group at all — is answered in §4b: it can, by making the group
the object. The proposal is not edited; the new one cites it.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 35 controls pass (`tests/vertical/test_connectors_adhesion.py`, 0.38 s, 2026-07-31); three mutants applied and all three killed (§4a-mutation). `tests/ports/` is 15 passed / 2 failed and both failures are foreign — see §8. Status stays `PROPOSED`: two of the four connectors this entry names are wired only under a reading of the wiring table that has not been ratified (§14 and the proposal it cites), and no `Binding` row has been added, because `aleph/vertical/wiring.py` has a single writer and this is not it. |
| Reviewer | Agent-proposed. Unratified. No PI review. |
| Rollback | Delete `aleph/vertical/connectors_adhesion.py` and `tests/vertical/test_connectors_adhesion.py`, and remove the four `Binding` rows if they were added. Nothing else imports either file; `focal_adhesion.py`, `wiring.py` and `aleph/state/**` are untouched by this entry. |

### §4a-mutation — mutation testing

Three mutants applied to the shipped module, controls run against each, module restored and verified
byte-identical with `git diff`.

| # | Mutant | Tests killed | Survived? |
|---|---|---|---|
| M1 | shipped path uses `k_a + k_l` (the parallel combination) instead of `k_a k_l/(k_a+k_l)` | **6** | no |
| M2 | the open-clutch gate removed: the stiffness accessor always ignores `is_engaged` | **3** | no |
| M3 | `direction_la = -direction_al`, i.e. the ligand block built as the negation of the actin block | **1** | no |

Killed by M1: `test_the_connector_site_force_equals_the_owner_s_own_joint_gradient`,
`test_dispatching_the_two_rows_independently_breaks_the_newton_pair`,
`test_dispatching_the_two_rows_independently_reports_the_parallel_stiffness`,
`test_an_open_clutch_contributes_exactly_zero_to_both_blocks`,
`test_an_open_clutch_that_still_reports_load_is_caught`,
`test_the_forbidden_sense_audit_is_a_measurement_and_not_a_structural_zero`.
M1 also kills the two negative controls, which is the right shape: with the shipped path already
parallel, the "broken" branch is no longer distinguishable from it, and a negative control that can
no longer tell the two apart must fail.

M3's kill count of one is honest rather than reassuring. The two constructions agree **bit-for-bit**
(§10), so no numerical control can see the difference; only
`test_the_two_force_blocks_are_accumulated_independently`, which reads the source, catches it. That
control exists because `AdjointPair`'s contract forbids the assignment for a reason no numerical
check can express — it makes the closure check unable to fail.

Restored afterwards and verified byte-identical against a pre-mutation copy (`diff`, empty); the
suite returns to **35 passed**.

## 14. Honest limits — what this entry does NOT establish

**No kinetics.** Nothing here binds, unbinds, matures or turns over. All four connectors are family
`[K]` and every one of them is wired **without** the kinetic law its family names, because no rate law
exists in this project: the registered contract owes the integrin/collagen binding rate, its
load-dependent unbinding law, the actin-side anchoring kinetics and the maturation transition rates.
`propose_accepted_load_pn` is a channel with nothing on the far end of it. Absent, not implied.

**No traction claim.** `focal_adhesion`'s contract lists adhesion size, area, shape, growth, count and
traction under `unsupported_claims`, and delivering a force through a joint does not remove any of
them. Every stiffness in the controls is chosen to make the algebra sharp; none is a cell number.

**The two composite rows are wired only under an unratified reading of the wiring table.** One
object serves both. Whether two `Binding` rows may name one dotted path for one *object* — as opposed
to two instances of one class, which `ecm_crosslink`/`dorsal_arc_crosslink` and
`nmii_sf_motor`/`nmii_cortex_motor` already do — is what
`docs/decisions/PROPOSAL-a-composite-group-is-one-connector-object.md` asks. If the PI decides
otherwise, the honest count from this lane is **two** and not four.

**The `endpoints` a row declares are not the owners its forces land on**, for all four rows. The
registry declares `(sf_arc, focal_adhesion)`, `(focal_adhesion, ecm)`,
`(lamellipodium, focal_adhesion)` and `(filopodium, focal_adhesion)`; the force blocks land on the
actin owner and on `ecm`. That is a consequence of `focal_adhesion` having no degrees of freedom and
is asserted by a control rather than left to be discovered, but it means the registry's `endpoint_b`
is not, for these rows, "the thing that receives the reaction". Reported.

**The nascent rows are given a discipline the registry does not require of them** (§4d). If that is
wrong, it is wrong in the conservative direction — the module refuses an evaluation the registry
would have permitted — but it is a position and not a derivation.

**Nothing here is stepped inside a world.** No `assembly.py` schedule runs any of these connectors,
so nothing establishes that an adhesion carries load in a relaxation, only that the object computes
and delivers what the law says. Evidence rung `ANALYTIC_ORACLE`, quantitative status `BLOCKED`.

**The maturation magnitudes are `UNVERIFIED`.** `ALEPH-PORT-3201` says the maturation scale factors
are placeholders; this connector multiplies through them and inherits that status exactly.
