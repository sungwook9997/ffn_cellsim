# ALEPH-PORT-2901 — the `alpha2beta1_collagen_series` composite joint and the geometry-less clutch graph

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2901` |
| Lane | `L29 focal_adhesion` |
| Status | `REJECTED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Rejected | `2026-07-30 18:0x KST` by session `bb6ec49e`, PI instruction item P3 |
| Port class | `RE-DERIVED` |

---

> ## REJECTED — the controls named below were never written
>
> `REJECTED` is `ports/README.md` §6's status for an entry that will not be ported, where "the
> reason is the deliverable". The reason is this box. Note what it does **not** mean here: the
> focal-adhesion capability is still wanted and will be rebuilt. It is *this entry* that is dead,
> because its evidence never existed. The rebuild gets a new port ID, not a revival of this one.
>
> **Do not read anything past this box as evidence.** Every section is retained exactly as written
> so the rewrite can see what the first attempt intended, but no claim in this entry is supported.
>
> This entry names six controls in two files:
> `tests/vertical/test_focal_adhesion.py` and `tests/vertical/test_focal_adhesion_motor_clutch.py`.
> **Neither file has ever existed in this repository.** The modules the entry authorises,
> `aleph/vertical/focal_adhesion.py` and `aleph/vertical/focal_adhesion_motor_clutch.py`, were
> 94 KB of code that no test anywhere imported. They are archived at
> `_archive/focal_adhesion_unverified_2026-07-30/` with a README explaining the rewrite rules.
>
> The port-discipline gate did not catch this because `test_named_controls_resolve_to_real_tests`
> skipped entries that were not `ACCEPTED`, and this one was `PROPOSED`. A `PROPOSED` entry naming a
> test that does not exist is exactly the case worth catching, since it is the state in which
> unverified work is easiest to mistake for verified work. The gate now checks named controls at
> every status.
>
> **The registered `focal_adhesion` component contract is not withdrawn.** It lives in
> `aleph/state/census_actomyosin.py` and is unaffected; so are the `fa_actin_anchor` and
> `integrin_collagen_clutch` connector contracts. The census was not what was wrong. A rewrite gets
> a **new** port ID with its ledger written before its code — not a revival of this one.
>
> Numbers appearing below, in particular the "clear interior peak at 1000 pN/µm", were carried on a
> lane's report and were never re-measured. Treat them as unverified.

---

## 1. Aleph API

The exact public surface this entry authorises.

```python
from aleph.vertical.focal_adhesion import (
    ActinSideState,
    Alpha2Beta1CollagenSeries,
    ClutchGraph,
    ClutchKinetics,
    ClutchSiteField,
    HalfCompliance,
    HalfKind,
    LigandSideState,
    MaturationState,
    SeriesJoint,
    SeriesLawViolation,
    SiteReference,
    assert_dispatched_as_one_joint,
    assert_geometry_less,
    assert_series_law,
    assert_single_series_path,
    owned_state_keys,
    parallel_over_series_ratio,
    series_compliance_um_per_pn,
    series_stiffness_pn_per_um,
)
from aleph.vertical.focal_adhesion_motor_clutch import (
    MotorClutchParameters,
    MotorClutchResult,
    StiffnessSweep,
    motor_clutch_run,
    motor_clutch_stiffness_sweep,
)
```

Nothing outside this list is covered by this entry.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/load_path.py`; `ffn_sim/ac/motor/two_filament_reference.py` |
| Source symbol(s) | `resolve_fa_series_group`, `resolve_fa_series_joint`, `MechanicalGroupResolution`, `CompositeFAJointSpec`, `ResolvedJoint`; and `series_stiffness` in the motor module |
| Read from | working tree — **stated explicitly** |
| Working tree == commit? | `yes` — verified with `git diff be0e5876 -- <each path>`, both empty. The repository as a whole carries 23 modified and 8 untracked paths at that commit, so this row is a per-file statement and not a repository-wide one. |

## 3. Why source-derived porting beats clean-room

**It does not.** Nothing was ported. The class is `RE-DERIVED` and the reason is specific rather
than dogmatic:

* The **law** this lane needed — two elastic elements in series share one load and add their
  extensions, so `1/k = 1/k1 + 1/k2` — is one line of secondary-school mechanics. Deriving it takes
  less time than reading somebody else's spelling of it. There is no numerically delicate form, no
  case analysis, and no empirically discovered failure mode in it.
* The **prohibition** — the two semantic edges must resolve to one mechanical joint — is already
  registered in Aleph's own contract (`aleph/state/connectors_surface_traction.py`,
  `ALPHA2BETA1_COLLAGEN_SERIES`, `composite_groups()`), which was written by Aleph's own census lane
  from Appendix A. It did not need to come across the boundary either.

What the source reading **did** establish, and what makes the entry worth keeping, is a defect
finding. It is recorded in §6.

## 4. Physical or mathematical law represented

Two compliant elements `1` and `2` connected end to end between a driving point `A` and an anchor
point `B`, with a massless junction `J` between them (the adhesion plaque — which in a geometry-less
representation is exactly a junction and nothing else).

Force balance on the massless junction:

    F_1 = F_2 = F                                          (i)

because the junction stores no momentum and has no other attachment. This is the entire physical
content of "one series joint": there is one load, not two.

Kinematics along the load path:

    x_total = x_1 + x_2                                     (ii)

Constitutive law for each half, written in compliance form:

    x_i = c_i F                                             (iii)

Substituting (iii) into (ii) with (i):

    x_total = (c_1 + c_2) F

so the joint's compliance is the **sum of the half compliances**, and its stiffness is

    k = 1 / (c_1 + c_2) = 1 / (1/k_1 + 1/k_2)                (iv)

Consequences that are the reason the composite rule exists, each derived rather than asserted:

* **(iv) is a harmonic mean, so `k <= min(k_1, k_2)`.** The joint is *softer* than either half. A
  series joint is governed by whichever half is softer.
* **Evaluating the halves independently gives `k_parallel = k_1 + k_2`**, because that is what two
  springs sharing one pair of endpoints do. The ratio is

      k_parallel / k_series = (k_1 + k_2)^2 / (k_1 k_2) >= 4

  with equality only at `k_1 = k_2`. So the *smallest possible* error from the mistake is a factor of
  four in stiffness, and it grows without bound as the two halves separate in scale. This is the
  quantity `parallel_over_series_ratio` reports.
* **Disengagement is the sharp case.** An unbound half has `c = +inf`. Then (iv) gives `k = 0`
  exactly: the load path is open and the joint carries exactly zero. The independent-spring
  evaluation gives `k_parallel = k_engaged + 0 = k_engaged`, i.e. the actin anchor still reports its
  full load while the integrin/collagen clutch is open. A traction reading with nothing on the other
  end of it. The series form gets this right for free; the parallel form gets it wrong by 100% of the
  reported traction.
* **Load sharing is a check, not a definition.** Given the joint's total extension, `F = k x_total`
  and then `x_1 = c_1 F`, `x_2 = c_2 F`. The identity `x_1 + x_2 = x_total` is an independent
  assertion on the implementation because it is (ii) recovered from (iv), not (iv) restated.

The **stochastic clutch** part is a separate law and is not a mechanics identity. Per bond, per
candidate step of length `dt`, with a state-dependent hazard `k`:

    P(transition in dt) = 1 - exp(-k dt)

which is the exact first-event probability of a Poisson process of constant rate over the step, not
the linearisation `k dt`. Unbinding is load-dependent through a Bell hazard,

    k_off(F) = k_off0 * exp(F / F_b)

which is the standard single-barrier form: the applied load tilts the escape barrier by `-F d`, and
`F_b = k_B T / d` is the force scale absorbing the reaction-coordinate length `d`. Aleph carries
`F_b` as the primitive because `d` is not measurable here and `k_B T` may not be defaulted (L1's
standing rule: both temperatures are `ASSUMED`, so no entry point defaults to `k_B T`).

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| half stiffness `k_i` | pN/µm | N/m (×1e-6) | `(0, +inf)` — a zero-stiffness half is `UNBOUND`, expressed as infinite compliance, not as `k=0` |
| half compliance `c_i` | µm/pN | m/N (×1e6) | `(0, +inf]`, `+inf` exactly for a disengaged half |
| joint stiffness `k` | pN/µm | — | `[0, min(k_1,k_2)]`, exactly `0.0` when either half is open |
| separation, extension | µm | m (×1e-6) | `> 0` for separation; extension may be either sign |
| force, rupture force scale | pN | N (×1e-12) | force any sign; `F_b > 0` |
| rate | 1/s | 1/s | `>= 0`; `0` means the transition is declared not to occur |
| retrograde speed | µm/s | m/s (×1e-6) | `>= 0` |

Singular and boundary cases, each with the behaviour Aleph requires:

- **Either half disengaged.** `total_compliance = +inf`, `stiffness == 0.0` exactly,
  `scalar_force_pn == 0.0` exactly. Testable with `==`, because it is zero on the whole branch
  rather than small near a gate.
- **Both halves disengaged.** Same, and the two `inf` values add to `inf` rather than to `nan`.
- **A half constructed with `k <= 0` or non-finite.** Refused at construction. A zero-stiffness bound
  half and an unbound half are different physical statements and must not share a spelling.
- **Coincident endpoints.** Refused: there is no line of centres, so the force has no direction. Not
  regularised.
- **`dt` non-positive or non-finite.** Refused by the kinetics.
- **A rate with no declared value.** Refused: `ClutchKinetics` has **no defaults on any field**,
  because the registered contract says a consumer that needs a rate law the census does not carry
  must refuse rather than default.
- **A site reference given as an integer.** Refused. An integer index into another owner's array is
  ownership by another name; references are opaque strings.

Invariants that must hold, each with the test that asserts it:

- **I1. The series law, exactly.** `1/k == 1/k_1 + 1/k_2` to floating-point round-off, over decades
  of stiffness ratio. `tests/vertical/test_focal_adhesion.py::TestSeriesLaw::test_series_law_is_exact`
- **I2. Load sharing and extension additivity.** `x_1 + x_2 == x_total` and `F_1 == F_2`.
  `...::TestSeriesLaw::test_one_load_and_additive_extensions`
- **I3. Softer than either half.** `k <= min(k_1, k_2)`.
  `...::TestSeriesLaw::test_series_joint_is_softer_than_either_half`
- **I4. An open half opens the joint, exactly.** `...::TestClutchOpensTheLoadPath::test_open_half_gives_exactly_zero`
- **I5. A half cannot be evaluated on its own.** Asking a `HalfCompliance` for a stiffness, a force,
  or an energy raises `SeriesLawViolation`.
  `...::TestHalvesCannotBeEvaluatedIndependently`
- **I6. Geometry-lessness.** No array owned by `ClutchGraph` carries a spatial axis, and the owner
  exposes no positions. `...::TestGeometryLess`
- **I7. Rejected-step integrity.** After a rejected candidate in which bonds demonstrably flipped,
  every clutch array and the RNG stream position are bit-identical to before.
  `...::TestRejectedStepIntegrity`
- **I8. The state block matches the registered contract.** `owned_state_keys()` equals
  `census_actomyosin.FOCAL_ADHESION.owned_state`. `...::TestContract::test_owned_state_matches_the_registered_contract`

## 6. Source evidence class and known retractions

**Where I looked:** `ffn_sim/ac/engine/load_path.py` in full; `grep` for `series`, `compliance`,
`series_stiffness`, and `1.0 / (1.0 /` across every `*.py` in the source tree; the source's own test
file for the motor reference; and `ports/audit/*.md` for a prior Aleph audit of this asset (there is
none — this lane is the first to read `load_path.py`).

What the source claims for this code: its module header declares a five-part sanity gate including a
`double-count` clause — "the ACTIN_ANCHOR and FA_CLUTCH semantic edges resolve to exactly one spring
energy and one force scatter". The symbol is reachable from the source's declared VS-0 seam and has
live tests. It self-labels as "not a complete FA maturation or partner-search implementation".

**The defect finding, which is the reason this entry exists.**

The source enforces the *prohibition* and never states the *law*.

* `resolve_fa_series_group` validates the group's membership, families, scope and chemistry cards,
  then **measures** how many mechanical joints the declared connectors resolve to and raises
  `"FA semantic edges would double-count mechanical stiffness"` when that count is not one. That
  guard is real and it is well built — its own comment records that the field it reads previously
  carried a default which made the guard structurally unable to fire, and that counting by group
  name instead of by component pair would have been a second dead guard. Both observations are
  correct and both are about *reachability of a guard*, which is the same class of defect Aleph's
  coverage lane found elsewhere in the same tree.
* `resolve_fa_series_joint` then constructs the one composite joint from
  `CompositeFAJointSpec.stiffness_pn_per_um` — **a single scalar supplied by the caller.** The two
  halves never appear as two compliances anywhere. There is no actin-side stiffness, no ligand-side
  stiffness, and no composition. `1/k = 1/k1 + 1/k2` is never evaluated for this joint.

So the prohibition is enforced against a quantity the model does not compute. Nothing in that path
can be wrong about the series combination, because nothing in that path performs one — the composite
stiffness is an input, which is exactly the "stiffness fudge factor" outcome the composite rule
exists to prevent. A guard that forbids double-counting a number you were handed is not a series
joint.

**One correction to the brief that dispatched this lane**, because it overstates the finding and the
overstatement is checkable: the brief says the series law "appears nowhere in its tree". It does
appear — `ffn_sim/ac/motor/two_filament_reference.py:188` computes
`1.0 / (1.0 / k_xb + 1.0 / k_head_arm)` for the NMII crossbridge-plus-arm chain, with a docstring
that reasons correctly about compliance not being a force loss, and a source test asserting it to
`rel=1e-12`. The accurate statement is narrower and worse: **the law is present in the source and is
not connected to the focal-adhesion composite joint.** The tree contains an author who understood
series compliance, in a different module, while the FA path takes its stiffness as a given. That is a
wiring defect rather than a knowledge gap, and it is the more instructive version of the finding.

No retraction, supersession or `BLOCKED` label is attached to `load_path.py` itself. The whole-cell
artifact the engine composes into (`outputs/ac/cell_assembled/composed_world_v2.json`) self-labels
`quantitative_claim_status: BLOCKED`, which is the enclosing claim, not this symbol's.

## 7. Independent oracle or derivation

Four independent checks, none of which is agreement with the source.

1. **The closed form (iv), derived in §4 from force balance and kinematics.** `1/k = 1/k_1 + 1/k_2`
   is asserted to floating-point round-off across ten decades of stiffness ratio. The residual is
   measured and reported rather than merely bounded.
2. **An algebraic identity the implementation does not use.** `x_1 + x_2 == x_total` recovers the
   kinematic premise (ii) from the composed stiffness, so it fails if the composition is right by
   coincidence.
3. **The `k <= min(k_1, k_2)` inequality**, a property of the harmonic mean, which the parallel form
   violates in the opposite direction for every input. It is therefore a discriminating check and not
   a sanity check.
4. **A qualitative regime prediction with a known sign structure**, and it is the strongest available
   check on the stochastic part: a load-and-fail motor–clutch has traction that is **non-monotone**
   in substrate stiffness. The prediction is structural — soft substrates cannot build force at all,
   stiff substrates load bonds fast enough that Bell unbinding empties the clutch population — so a
   monotone curve falsifies the coupling between load and unbinding regardless of parameter values.
   No parameter was tuned to produce the peak, and the sweep reports where it landed rather than
   asserting where it should land.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_focal_adhesion.py::TestSeriesLaw::test_series_law_is_exact` | `abs(1/k - (1/k1 + 1/k2)) / (1/k1 + 1/k2) <= 4 * eps` for 121 stiffness pairs spanning `k1, k2` in `[1e-3, 1e7]` pN/µm |
| Positive | `tests/vertical/test_focal_adhesion.py::TestSeriesLaw::test_one_load_and_additive_extensions` | `x_1 + x_2 == x_total` to `<= 8 * eps` relative, and both halves carry the same load to exactly `0.0` difference |
| Positive | `tests/vertical/test_focal_adhesion.py::TestClutchOpensTheLoadPath::test_open_half_gives_exactly_zero` | `stiffness == 0.0` and `scalar_force_pn(...) == 0.0`, exactly, on the whole open branch |
| Positive | `tests/vertical/test_focal_adhesion.py::TestRejectedStepIntegrity::test_rejected_candidate_restores_every_bond_and_the_stream` | ≥1 bond flipped in the candidate; after rejection every one of the six state arrays is bit-identical and the RNG bit-generator state and draw count match exactly |
| Positive | `tests/vertical/test_focal_adhesion_motor_clutch.py::test_traction_is_non_monotone_in_substrate_stiffness` | the sweep's peak is interior to the swept range, and the peak exceeds both endpoints by a stated margin |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_focal_adhesion.py::TestNegativeControlIndependentSprings::test_parallel_stiffness_is_rejected_by_the_series_check` | `assert_series_law(k1 + k2, k1, k2)` raises `SeriesLawViolation`, and the raised message carries the measured ratio. Fails if the check stops discriminating. |
| Negative (must fail) | `...::TestNegativeControlIndependentSprings::test_a_half_refuses_to_report_a_stiffness` | `HalfCompliance.stiffness_pn_per_um`, `.force_pn(...)` and `.energy_pn_um(...)` each raise `SeriesLawViolation`. Fails if a half becomes independently evaluable. |
| Negative (must fail) | `...::TestNegativeControlIndependentSprings::test_open_clutch_with_parallel_evaluation_reports_phantom_traction` | with the ligand half open, the series joint reports exactly `0.0` while the parallel evaluation reports the full actin-side load; the test asserts the guard catches the parallel number. Fails if the open-path case stops being distinguishable. |
| Negative (must fail) | `...::TestDispatchGuards::test_dispatching_the_two_semantic_edges_separately_is_refused` | `assert_dispatched_as_one_joint(("fa_actin_anchor", "integrin_collagen_clutch"))` raises; and `assert_single_series_path` raises when a rival connector spans an already-spanned component pair. |
| Negative (must fail) | `...::TestRejectedStepIntegrity::test_a_kinetics_that_writes_authoritative_state_is_caught` | a deliberately broken kinetics that writes the accepted arrays instead of the candidate arrays survives a rejection, and the integrity assertion fails on it. Proves the integrity test is not vacuous. |
| Negative (must fail) | `...::TestGeometryLess::test_a_planted_position_array_is_caught` | `assert_geometry_less` raises when a `(N, 3)` array is planted on the graph. Fails if the geometry-less check stops looking. |

## 10. Numerical and precision envelope

Working and accumulation precision: IEEE-754 float64 throughout; no float32 anywhere in this port.

* **The series law** is asserted at `4 * eps` relative (`8.9e-16`). That number and not a looser one
  because the expression `1/(1/k1 + 1/k2)` costs three roundings and its inverse a fourth, so four
  ulp is the arithmetic budget with nothing spare. Measured worst case over the 121-point grid is
  reported in the test output and in the lane report. The stated tolerance holds for
  `k1, k2` in `[1e-3, 1e7]` pN/µm, a ten-decade span chosen so the harmonic mean is dominated by
  each half in turn.
* **Conditioning.** The compliance form is the well-conditioned one. `c1 + c2` is a sum of positives,
  so there is no cancellation at any ratio; the stiffness form `k1 k2 / (k1 + k2)` would overflow for
  `k1 k2 > 1e308` while the compliance form is exact in the same regime. This is why
  `HalfCompliance` stores compliance and not stiffness — the choice is numerical as well as
  type-level.
* **Extension additivity** is asserted at `8 * eps` relative: the round trip stiffness → force →
  two extensions → sum costs twice the series law's budget.
* **The open branch** is asserted at exactly `0.0`. Permitted because `+inf` compliance makes the
  reciprocal exactly `0.0` in IEEE-754 and the force is a product with that exact zero — it is zero
  on the whole branch, not small near a gate.
* **The kinetics** uses `1 - exp(-k dt)` rather than `k dt`. For `k dt` near `1` the linearisation
  is wrong by tens of percent and can exceed `1`, which silently becomes "always fires". Outside the
  domain — `dt <= 0`, non-finite `dt`, negative rate — the kinetics **refuses** rather than clamping.
* **The Bell factor** `exp(F / F_b)` overflows for `F / F_b > 709`. Aleph refuses a non-finite hazard
  rather than saturating it, because a saturated hazard is indistinguishable from a very large one in
  the output and a rupture force scale that small is a parameter error.
* **The motor–clutch sweep** is an explicit-Euler-in-time, quasi-static-in-force integration. Its
  time step is validated against two conditions and refused if either fails: the per-step transition
  probability must stay below a declared ceiling, and the per-step load increment must stay below a
  declared fraction of `F_b`. It carries **no** convergence claim; the peak's *location* is reported
  to within the sweep's own decade resolution and is explicitly not a converged number.

## 11. Production-backend residency and transfer

Host, float64 numpy, for all of it. The clutch graph is six flat arrays of length `n_clutches` plus
one integer generation counter — there is no spatial array to place, which is what geometry-less
buys. Per candidate step the graph consumes six vector RNG draws of length `n_clutches` and writes
six candidate arrays; nothing is transferred to a device and no host round-trip exists because there
is no device copy.

When the ratified backend acquires kernels, the residency answer changes for exactly one thing: the
per-clutch force evaluation would move device-side alongside the endpoint owners' position arrays,
and the six state arrays would follow it. The composition `1/(c1 + c2)` is two divisions and an add
per clutch and is not a reason to move anything. The kinetics draw count per step is fixed and
independent of how many bonds flip, which is a property worth keeping on a device: a variable draw
count makes the stream position depend on the trajectory and destroys bit-reproducibility across
backends.

`aleph/vertical/focal_adhesion_motor_clutch.py` is an experiment harness and will never run on the
production backend. It owns the 1-D scalar kinematics the component deliberately does not own, which
is also the reason it is a separate module: putting a position in the component would break the
geometry-less contract in the one place a reader would not look.

## 12. Comments and docstrings to discard

None of the source prose survives, because none of it was carried. Specifically discarded:

* The seam label `VS-0` and every "VS-0 mechanics requires…" phrasing. It is a provider milestone
  name and Aleph's scope discipline is `E`/`I`/`H`/`B`/`X`, decided independently.
* The five-part `Sanity Gate:` block in the module header, and its clause names
  (`topology` / `generation` / `double-count` / `conservation` / `transaction`). Aleph's equivalents
  are the acceptance predicate's named conditions, which already exist and already report structure.
* The provider's family and role vocabulary (`ACTIN_ANCHOR`, `FA_CLUTCH`, `EndpointRole`,
  `ConnectorScope.INTER_COMPONENT`, `JointKind.FA_COMPOSITE_SERIES`, `ElementKind.SEGMENT`,
  `chemistry_card`). Aleph's registry already names the same two edges from Appendix A, in Appendix
  A's words.
* The structure-of-arrays / device-upload commentary and the `_d` device-array suffix convention.
  Aleph has no device arrays here and will not inherit a naming convention for buffers it does not
  allocate.
* `fa_cluster_id`, and the fixed-capacity-SoA framing around it. A capacity is a build-time constant
  and this component's population is a dynamic output; the registered contract says so explicitly.
* Every reference to the source's own module paths.

What replaces them: the derivation in §4 written into the module docstring in Aleph's own terms; the
open-load-path consequence stated where the guard is, because that is where a reader meets it; and
the four-way regime argument for the motor–clutch prediction written where the sweep is defined.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Not accepted. `PROPOSED`. The controls in §8 and §9 pass, and the numbers are in the lane report, but this entry is agent-written and no reviewer has ratified it. Two substantive items are unresolved and are named in §14. |
| Reviewer | none — agent-proposed, unratified |
| Rollback | Delete `aleph/vertical/focal_adhesion.py`, `aleph/vertical/focal_adhesion_motor_clutch.py`, `tests/vertical/test_focal_adhesion.py`, `tests/vertical/test_focal_adhesion_motor_clutch.py`. Nothing imports them: `aleph/vertical/__init__.py` is owned by another lane and was not modified, so the modules are reachable only by explicit import. Removing them loses the only implementation of the series law in the tree and returns the composite rule to the state this entry documents in §6 — a prohibition with nothing behind it. |

## 14. Honest limits

What this entry does **not** establish:

1. **No rate is sourced.** Every kinetic constant is `UNSOURCED`, and the registered contract already
   says evidence is owed on the integrin/collagen binding rate, its load-dependent unbinding law, the
   actin-side anchoring kinetics, and the maturation transition rates. `ClutchKinetics` therefore has
   no defaults at all — a caller must state every rate — but stating a number is not sourcing it.
   **No claim about traction magnitude, adhesion lifetime, or maturation timing is licensed by this
   port.** Only the mechanics identities and the qualitative non-monotonicity are.
2. **The motor–clutch sweep is a regime demonstration, not a converged result.** One seed per
   stiffness point, one time step, one parameter set, no convergence study, and its parameters are
   assumed values chosen to sit in the load-and-fail regime. The peak's location carries a decade of
   resolution at best. Compare the level-3 Laplace episode: a single point is not a numerical result.
3. **Two open questions the docs do not answer, raised rather than settled.**
   * *What `occupancy` is.* Appendix A lists `occupancy` alongside the actin-side state, the
     ligand-side state and the maturation state, and does not define it. This port reads it as the
     per-clutch record of whether the series path is closed, stores it as owned state, and asserts on
     every read that it agrees with the conjunction of the two side states. That reading is a choice.
     It is not the only one: "occupancy" could equally mean an integrin count per plaque, or a
     population-level bound fraction. The invariant makes the choice checkable but does not make it
     right, and the alternative readings differ in what state exists.
   * *Whether the clutch bears compression.* Appendix A states the series path and does not state the
     sense. A molecular linkage argues for tensile-only; the motor–clutch literature this lane's
     regime prediction comes from uses a two-sided linear spring. This port implements both behind an
     explicit `tension_only` flag with **no default**, so the caller must say which, and neither is
     presented as the registry's answer. The sweep states which it used.
4. **The two halves' compliances are effective, not resolved.** The registered contract already
   concedes this: bond mechanics is one effective compliance per side rather than talin domain
   unfolding. So the joint is the right *topology* with an assumed *constitutive law* on each half.
5. **No connection to the rest of the traction spine.** `ecm_crosslink` and `ecm_far_field_anchor`
   are the other two members of connector group B and neither exists. So this joint's ligand side
   attaches to whatever the caller hands it, and there is no matrix behind it. The substrate spring
   in the sweep is a stand-in for that missing matrix and is not the `ecm` owner.
6. **Maturation has mechanical consequences that are declared, not derived.** Reinforcement enters as
   two multipliers — on the actin-side stiffness and on the off-rates — because a maturation state
   with no mechanical consequence is a flag, and this lane's brief was explicit that it must not be
   one. The *existence* of reinforcement is standard; **these multipliers are assumed** and no
   measurement constrains either.
7. **Nothing here has run on the ratified backend**, which ships no kernels. No GPU was used and none
   is authorised.
8. **`aleph/vertical/__init__.py` was deliberately not modified**, so these modules are absent from
   the package's public surface. That is a lane-boundary consequence, not a design decision, and
   whoever owns that file should decide whether they belong in it.
