# ALEPH-PORT-3302 — the six protrusion connectors: two Brownian ratchets, two transient seams, two motors

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3302` |
| Lane | protrusion connector lane, under lead session `f70de564` |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Registry this implements | `aleph/state/connectors_protrusion_motor.py` — `PROTRUSION_CONNECTORS`, `NMII_MOTOR_TARGETS` |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
# new module aleph.vertical.connectors_protrusion
ProtrusionConnectorError                    # base refusal
RatchetDomainError                          # a load-velocity pair off the physical branch

BrownianRatchet                             # the load-velocity element; frozen value type
BrownianRatchetConnector                    # the connector; evaluate_sites + accumulate
ProtrusionActinTarget                       # an owner's NMII binding sites, plus the way back

build_lamellipodium_membrane_contact        # -> BrownianRatchetConnector
build_filopodium_membrane_tip               # -> BrownianRatchetConnector
build_lamellipodium_cortex_seam             # -> CrosslinkConnector  (reused, not re-derived)
build_filopodium_cortex_root                # -> CrosslinkConnector  (reused, not re-derived)
lamellipodium_actin_binding_sites           # -> ProtrusionActinTarget
filopodium_actin_binding_sites              # -> ProtrusionActinTarget
build_nmii_lamellipodium_motor              # -> NmiiMotorConnector  (reused, not re-derived)
build_nmii_filopodium_motor                 # -> NmiiMotorConnector  (reused, not re-derived)
```

**Four of the six connectors get no new class, and that is the finding rather than a shortcut.**

* `lamellipodium_cortex_seam` and `filopodium_cortex_root` are `CrosslinkConnector` from
  `aleph/vertical/connectors_crosslink.py`, already wired for `sf_cortex_transient`. §4 gives the
  argument that it is the right element and the one property that would make it the wrong one.
* `nmii_lamellipodium_motor` and `nmii_filopodium_motor` are `NmiiMotorConnector` from
  `aleph/vertical/nmii.py`, already wired for `nmii_sf_motor` and `nmii_cortex_motor`. That class is
  parameterised by the actin owner it reaches and `MOTOR_CONNECTOR_TARGETS` already lists both new
  targets. What was missing was not a class; it was that no lamellipodium and no filopodium existed
  to be a target. They exist now. Wiring them is therefore a **projection plus a control**, not a
  wrapper, and inventing a subclass to make the diff look like new physics would be indirection for
  bookkeeping.

`aleph/vertical/wiring.py` is **not** written by this lane. The `Binding` rows this entry earns are
reported to the lead session, which is that file's single writer.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | **none named and none read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened for this work. |
| Working tree == commit? | not applicable — nothing was read |

`ALEPH-PORT-2101` §6 audited all six protrusion edges in the reference project and found every one a
**SEAM**: a declared contract, a name constant, a `Protocol` with one method, and no concrete
implementation anywhere outside a test double. There is nothing there to port even if the policy
allowed it, and the audit is quoted rather than re-verified — this lane read no reference file.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** The ratchet's load–velocity relation is the standard
polymerisation-ratchet form: the rate of monomer addition is suppressed by the Boltzmann factor for
the work the load does over one monomer step, and the reverse rate is not. Everything else in this
entry reuses an Aleph class that already exists.

## 4. Physical or mathematical law represented

### 4a. The Brownian ratchet — the load-velocity relation *is* the element

The registered mechanism, verbatim from `connectors_protrusion_motor.py`:

> "The thermal excursion of the membrane opens the gap a monomer needs; adding the monomer rectifies
> the excursion. Growth rate and carried load are therefore one coupled quantity, not a kinematic
> rate plus a separate contact force."

and the interpretation, on why the cheap substitution is not available:

> "replacing it with a generic non-penetration contact removes the coupling that makes protrusion
> generate force at all: a plain contact resists interpenetration but cannot convert polymerisation
> into work against a load, so the leading edge would push with a stiffness and not with a
> polymerisation velocity."

So a ratchet is **not** a `_DistributedUnilateralConnector` and it is not forced into one. A
unilateral spring's force is a function of a *separation*; the ratchet's is a function of a
*velocity*, and it has no rest length, no engagement gate and no stored elastic energy at all. The
nearest shape already in this package is `ImmersedDragConnector`, whose force is also set by a rate
rather than by a configuration; the ratchet is written against that precedent and not against
`ErmTether`.

The element, with `f` the load per **barbed end** [pN]:

```
v(f) = v0 * ( exp(-f d / kT) - e ) / (1 - e)              (R1)
f(v) = -(kT / d) * ln( e + (1 - e) v / v0 )               (R2)  the exact inverse of (R1)
```

| Symbol | Meaning | Unit | Where it comes from |
|---|---|---|---|
| `kT` | `k_B T` for this run | pN·µm | **the caller.** Class `ASSUMED`; no entry point defaults it |
| `d` | monomer step: how far a barbed end advances per subunit | µm | the caller |
| `v0` | load-free elongation velocity | µm/s | the caller — the owner **refuses** it (§6) |
| `e` | `k_off / (k_on C)`, the inverse supersaturation, in `(0, 1)` | — | the caller |

(R1) is a rectification statement and not a fit. The forward step must do work `f d` against the
load, so its rate carries `exp(-f d / kT)`; the reverse step — losing a subunit — does not, because
the load helps it. Normalising by the load-free value gives (R1) with `v(0) = v0` exactly.

**The stall force is the oracle, and it is closed-form:**

```
f_stall = (kT / d) * ln(1/e)                              (R3)
f_stall * d = kT * ln(1/e)                                (R4)
```

(R4) is the thermodynamic statement (R3) is worth having: the work the ratchet can do against a load
in one step equals the free energy released by adding one subunit. **`f_stall` does not contain
`v0`.** A ratchet made a thousand times faster stalls at exactly the same load, because `v0` is a
rate and stall is a thermodynamic limit. That is the independent check in §7, and it is the check a
model that had quietly become a stiffness could not pass.

### 4b. Where the filopodial case differs — the load lands, not the law

From the registry:

> "Concentrating the ratchet at a tip is the mechanical difference from the lamellipodium: the load
> per barbed end is set by how many filaments the bundle carries to the tip, and the membrane sees a
> near-point load with a high local curvature cost rather than a distributed pressure."

So the *element* is identical and the **carrier count** differs. One connector class serves both,
with a per-site barbed-end count `n_i`:

```
site load     F_i = n_i * f( v_i )                        (R5)
site velocity v_i is the velocity of one barbed end, not of the site
front stall   F_stall,total = (sum_i n_i) * f_stall       (R6)
```

The lamellipodial front is many sites of `n_i = 1` — one growing barbed end each, spread across a
branched leading edge. The filopodial tip is **one** site with `n_i` equal to the number of filaments
the bundle carries to the tip, which is exactly `FilopodiumBundle.barbed_end_count`. (R6) is the
sharp consequence and it is testable: a bundle of `n` filaments stalls at `n` times the single-
filament stall, and its *velocity* at a given total tip load is the single-filament velocity at
one `n`-th of it. A model that divided the load anywhere else would give a bundle that stalls like
one filament.

**The membrane's local curvature cost is the membrane's**, not this connector's. `HelfrichMembrane`
owns bending; a near-point load produces a high curvature cost there by construction. Nothing here
computes a curvature and nothing here is evidence about one — see §14.

### 4c. The seam and the root — the same transient element, twice

Both contracts carry the disjoint-population rule verbatim, and both are emphatic about what they are
not: *"NOT a shared filament, NOT a shared node, NOT a permanent weld."*

`CrosslinkConnector` is a two-sided spring between two material points of two owners, carrying a
`bound` flag, snapshot/commit/rollback and `accumulate`. It is already wired across owners for
`sf_cortex_transient`, whose registered mechanism is "transient actin crosslinking" between the
stress-fibre and cortical populations. The seam and the root are the same object:

* two disjoint filament populations that stay disjoint — the connector never touches a node index of
  the far side, only a material-point sink the far owner issues;
* a **finite stiffness**, so the lamellipodium cannot hold the cortex rigidly, which the contract
  says a transient connection cannot do;
* a `bound` flag whose clearing gives **exactly zero** transmitted load, which is what "can release"
  means when written as arithmetic.

**Two-sided rather than tension-only, and the reason is not symmetry with `sf_cortex_transient`.** A
filopodium under tip load transmits compression down the bundle to its root. A tension-only root
would transmit exactly nothing in that state, and the bundle would push through the cortex — which
is the same defect as the shared-node one, arriving from the other side. The rear seam is the same
argument with the sign reversed: the lamellipodial network is pushed rearward by the membrane it
loads, so the seam is loaded in compression precisely when protrusion is happening.

**What would make it the wrong element**, stated so this reuse can be refuted: the registry says both
have *"their own binding and unbinding state"*. `CrosslinkConnector` carries the `bound` flag and
**no kinetics** — nothing in it binds or unbinds. So the element is right and the *population
dynamics are absent*, exactly as they are absent for `sf_cortex_transient`. Nothing in this entry is
evidence about seam lifetime, about a root's detachment rate under load, or about how either
population turns over. Absent, not implied.

### 4d. The two motor connectors

`NmiiMotorConnector.evaluate_sites` computes `k_xb * x` in the **head's own strain**, which is
population state advanced by the kinetics, and *verifies* the configuration it is handed rather than
ignoring it — two heads at identical anchor-to-site separations carry different forces at different
points in the stroke. Nothing about that is specific to which actin owner is on the far side, which
is why one class serves four connectors.

What this lane adds is the **projection**: an owner publishes `NMII_HEAD`-role attachment sites, and
a motor connector needs `ActinBindingSites` — positions, **unit** polar tangents, material
coordinates and filament ids. `ProtrusionActinTarget` builds that block from the owner's own
published sites and keeps the site ids, so `deliver_forces_pn()` can scatter the accumulated actin-
side forces back through the owner's `endpoint_sink` handles, with the owner's own interpolation
weights and its own stale-handle guard. Without that last step the motor's reaction lands in a
detached array and the load path stops at the boundary.

## 5. Units, domains, singular cases, invariants

Length µm, force pN, energy pN·µm, time s. Velocity µm/s, `k_B T` pN·µm, stiffness pN/µm.

| Quantity | Aleph unit | Domain |
|---|---|---|
| `kb_t_pn_um` | pN·µm | finite, `> 0` |
| `monomer_step_um` | µm | finite, `> 0` |
| `free_velocity_um_per_s` | µm/s | finite, `> 0` |
| `reverse_ratio` | — | `0 < e < 1` |
| `barbed_ends_per_site` | — | integer `>= 1`, one per site |
| `growth_velocity_um_per_s` | µm/s | `v_min < v <= v0`, with `v_min = -e v0 / (1 - e)` |

Singular and boundary cases, each with the behaviour required:

- `e = 0` — **refused.** A ratchet with no reverse step never stalls: (R3) diverges and the element
  would push arbitrarily hard. A model that cannot stall is not a ratchet, it is a source.
- `e >= 1` — **refused.** The barbed end shrinks at zero load; there is no protrusion to speak of and
  (R2) is negative at `v = v0`.
- `v > v0` — **refused** unless the deliberate-break flag `allow_ratchet_to_pull` is set. `f(v) < 0`
  there: the ratchet would *pull* the membrane inward. A ratchet is a one-way device, and a negative
  load is the signature of it having become an ordinary spring.
- `v <= v_min` — **refused.** The log's argument is non-positive; there is no finite load that
  produces it. Depolymerisation cannot exceed the bare off-rate however hard you push.
- `kb_t_pn_um` absent — **there is no default and no module-level constant.** Every entry point takes
  it. `PLAN.md` §6 classes both registered temperatures `ASSUMED`, so a defaulted `kT` would be an
  assumption entering through a signature nobody reads.
- `free_velocity_um_per_s` absent — same: no default. `lamellipodium.rate_prior(...)` **raises** for
  `free_barbed_end_elongation_rate_um_per_s`, and this module does not route around its own owner's
  refusal.
- `evaluate_sites` before any velocity is set — **refused.** Returning zeros would make "nothing was
  evaluated" and "stalled" the same answer, and at `v = 0` this element is at its *maximum* load, so
  the confusion is in the dangerous direction. `PLAN.md` §6.1 records the afternoon a cold cache
  returning `0.0` cost.
- A site block of the wrong length — refused, naming the site count the connector holds.
- A non-unit ratchet axis — refused. The axis *is* the filament polarity and it alone sets the
  direction of the force, so an axis of length 0.9 would rescale the load by 10% and read as a soft
  ratchet.
- A bundle with zero barbed ends at the tip — refused by the owner
  (`FilopodiumBundle.load_per_barbed_end_pn`) and refused here: a ratchet with nothing to rectify
  against is not a weak ratchet.
- A motor builder pointed at the wrong owner — refused. Naming a connector `nmii_lamellipodium_motor`
  while its target's `owner` is `cortex` would put load on one population under another's name.

Invariants, each with the control that asserts it:

- **I1.** `v(f(v)) == v` and `f(v(f)) == f` to round-off — the two directions are one curve —
  `test_the_load_and_the_velocity_are_one_curve_not_two`.
- **I2.** `v(0) == v0` exactly — `test_the_zero_load_velocity_is_the_free_velocity`.
- **I3.** `f(0) == f_stall == (kT/d) ln(1/e)`, against the closed form —
  `test_the_stall_force_matches_the_closed_form`.
- **I4.** `f_stall` is invariant under scaling `v0` over three decades —
  `test_the_stall_force_does_not_depend_on_any_rate`.
- **I5.** `f_stall * d == kT ln(1/e)`, the thermodynamic identity — `test_the_work_per_monomer_at_stall_is_the_free_energy_per_monomer`.
- **I6.** A bundle of `n` barbed ends stalls at `n` times a single filament's stall —
  `test_the_bundle_stalls_at_n_times_the_single_filament_stall`.
- **I6b.** The stall the connector *reports* is the load it *delivers* at zero velocity — two code
  paths for one fact — `test_the_reported_stall_is_the_load_actually_delivered_at_zero_velocity`.
  Added because mutation testing showed I6 alone did not cover the multi-site front; see §13a.
- **I7.** The two sides are equal and opposite per site, assembled separately —
  `test_the_two_sides_are_equal_and_opposite_per_site`.
- **I8.** The ratchet's force does not depend on the separation it is handed —
  `test_the_ratchet_force_does_not_depend_on_the_separation`.
- **I9.** The active power vanishes at both ends of the curve and is positive between —
  `test_the_active_power_vanishes_at_both_ends_and_is_positive_between`.
- **I10.** The seam and the root carry load in **both** senses —
  `test_the_seam_and_the_root_carry_load_in_both_senses`.
- **I11.** An unbound seam or root transmits **exactly** `0.0` —
  `test_an_unbound_seam_or_root_transmits_exactly_zero`.
- **I12.** The owner's own root gate and the connector's `bound` flag agree —
  `test_the_owners_detached_root_and_the_connectors_unbound_flag_agree`.
- **I13.** Every one of the six carries both `evaluate_sites` **and** `accumulate` —
  `test_every_protrusion_connector_carries_both_protocol_methods`.
- **I14.** The motor's reaction reaches the owner's nodes with the right resultant **and** the right
  moment — `test_the_motor_reaction_reaches_the_owner_with_force_and_moment`.
- **I15.** Each connector names the declared contract with the declared endpoints —
  `test_each_connector_names_a_declared_contract_with_the_declared_endpoints`.
- **I16.** Four of the six are literally the reused classes, checked as types rather than asserted
  in prose — `test_four_of_the_six_reuse_a_class_that_already_existed`.

## 6. Source evidence class and known retractions

No claim is made about the reference implementation beyond `ALEPH-PORT-2101` §6's audit, which is
quoted, not re-verified.

Evidence rung `ANALYTIC_ORACLE`, quantitative status `BLOCKED`. **No magnitude in this module or its
controls is a measurement of a cell.** `k_B T` is `ASSUMED` in this project and is never defaulted
here; the monomer step, the free elongation velocity and the reverse ratio are all supplied by the
caller. `ALEPH-PORT-3205` §14 records that all five of the lamellipodium's rate priors are evidence
owed and that `rate_prior` refuses every one of them, and **this module ships no rate**: it does not
carry, default, or infer a polymerisation velocity. `ALEPH-PORT-3204` §4b records the same for the
filopodial tip elongation rate.

## 7. Independent oracle or derivation

1. **The stall force (R3), and its independence from every rate.** At stall the velocity is zero, so
   `v0` has cancelled out of the problem entirely: what remains is the thermodynamic maximum work a
   subunit addition can do, `kT ln(1/e)`, divided by the step `d`. A control scales `v0` by `10^-3`
   and `10^3` and asserts the stall force is unchanged to round-off. **A stiffness cannot do this**:
   a spring's blocking force is proportional to its stiffness, so anything that had degraded into a
   contact would move.
2. **The two ends of one curve.** `v(0) = v0` and `f(0) = f_stall` are the same relation evaluated at
   its two limits, and `v(f(v)) = v` is an algebraic identity that holds independently of every
   constant. It cannot be satisfied by a law that is right at one end and fitted at the other.
3. **The thermodynamic identity (R4).** `f_stall * d = kT ln(1/e)` is dimensionally an energy and is
   independent of `d` once expressed per monomer. A model that got the exponent's argument wrong —
   `f d / 2kT`, say, which is a real and common slip — reproduces neither (R3) nor (R4).
4. **The carrier-count scaling (R6)** is an independent structural oracle: it involves no constant of
   the element at all, only how many barbed ends share the load.
5. **The transient element against its owner's own gate.** `FilopodiumBundle.root_transmitted_force_pn`
   returns an exact zero array when the root is detached. The connector's `bound[:] = False` must give
   the same exact zero, and the two implementations are independent.

## 8. Positive control

Named here; **measured numbers are recorded in §13 after the run**, never from memory.

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheRatchetIsALoadVelocityRelation::test_the_load_and_the_velocity_are_one_curve_not_two` | I1 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheRatchetIsALoadVelocityRelation::test_the_zero_load_velocity_is_the_free_velocity` | I2 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheStallForceIsTheOracle::test_the_stall_force_matches_the_closed_form` | I3 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheStallForceIsTheOracle::test_the_stall_force_does_not_depend_on_any_rate` | I4 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheStallForceIsTheOracle::test_the_work_per_monomer_at_stall_is_the_free_energy_per_monomer` | I5 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheFilopodialCaseDiffersInWhereTheLoadLands::test_the_bundle_stalls_at_n_times_the_single_filament_stall` | I6 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheFilopodialCaseDiffersInWhereTheLoadLands::test_the_reported_stall_is_the_load_actually_delivered_at_zero_velocity` | I6b |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheRatchetConnectorDelivers::test_the_two_sides_are_equal_and_opposite_per_site` | I7 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheRatchetConnectorDelivers::test_the_ratchet_force_does_not_depend_on_the_separation` | I8 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheRatchetConnectorDelivers::test_the_active_power_vanishes_at_both_ends_and_is_positive_between` | I9 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheSeamAndTheRootAreTransientNotWelds::test_the_seam_and_the_root_carry_load_in_both_senses` | I10 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheSeamAndTheRootAreTransientNotWelds::test_the_owners_detached_root_and_the_connectors_unbound_flag_agree` | I12 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheMotorConnectorsReachTheNewOwners::test_the_motor_reaction_reaches_the_owner_with_force_and_moment` | I14 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheWiringClaimIsReal::test_every_protrusion_connector_carries_both_protocol_methods` | I13 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheWiringClaimIsReal::test_each_connector_names_a_declared_contract_with_the_declared_endpoints` | I15 |
| Positive | `tests/vertical/test_connectors_protrusion.py::TestTheWiringClaimIsReal::test_four_of_the_six_reuse_a_class_that_already_existed` | I16 |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_connectors_protrusion.py::TestTheNegativeControlsFail::test_a_ratchet_that_pushes_with_a_stiffness_is_caught_by_the_stall_oracle` | A stiffness-based leading edge — the exact substitution the contract forbids — has a blocking force proportional to its stiffness and therefore **fails I4**. Built and shown to fail. |
| Negative (must fail) | `tests/vertical/test_connectors_protrusion.py::TestTheNegativeControlsFail::test_a_ratchet_allowed_to_pull_is_caught` | Driving the shipped `allow_ratchet_to_pull` flag, the connector delivers a negative load and pulls the membrane inward; the healthy branch **refuses** instead. |
| Negative (must fail) | `tests/vertical/test_connectors_protrusion.py::TestTheNegativeControlsFail::test_a_seam_degraded_to_a_tether_stops_carrying_the_compression` | A tension-only seam built from the same machinery carries exactly zero in compression, so I10 distinguishes the two rather than passing on both. |
| Negative (must fail) | `tests/vertical/test_connectors_protrusion.py::TestTheNegativeControlsFail::test_a_weld_with_no_release_is_distinguishable_from_a_transient_root` | Driving the owner's shipped `allow_detached_root_to_transmit_load`, a detached root transmits load — the weld defect — where the healthy branch is exactly inert. |
| Negative (must fail) | `tests/vertical/test_connectors_protrusion.py::TestTheNegativeControlsFail::test_the_carrier_count_dropped_would_change_the_bundle_stall` | The counterexample for (R6): a one-barbed-end tip and a sixteen-barbed-end tip must not report the same stall, so the carrier-count control is known to be able to fail. |

Each negative control asserts in **both** directions: the broken branch misbehaves **and** the healthy
branch is exactly inert (`== 0.0`, the number, not a small number).

## 10. Numerical and precision envelope

float64 throughout; no reduced-precision path. The three refusals above (`e` out of range, `v` out of
range, missing velocity) keep the `exp` and the `ln` in (R1)/(R2) on branches where neither
overflows: `f d / kT` is bounded above by `ln(1/e)` on the physical branch, so the exponential's
argument lies in `[-ln(1/e), 0]`.

`v(f(v)) = v` is asserted at `rel=1e-12`, not `==`: the round trip crosses one `exp` and one `log`
and a bit-identical claim would be a claim about libm rather than about the physics. The **exact**
assertions are reserved for the places where exactness is the physics: an unbound seam transmits
`0.0` and a detached root transmits `0.0`, both `np.array_equal` against zeros, because "cannot
transmit" and "transmits a little" are different models.

Force closure on the ratchet is reported against the **constituent** scale
(`sum_i |f_i|`), never the resultant. `PLAN.md` §2.5 records a resultant-based tolerance rejecting
40,000 consecutive steps of a correct relaxation, and a leading edge whose barbed ends point in many
directions has exactly that shape: large per-site forces whose vector sum is small.

`central_force=False` is reported for the ratchet, honestly. The force acts along the **filament
axis**, not along the line of centres between the barbed end and the membrane quadrature point.
Declaring it central would put the moment residual under a check it must fail for correct physics.
The per-site moment is instead a diagnostic of registration quality: it vanishes exactly when the
membrane contact lies ahead of the barbed end along the axis, and grows with the tangential
misregistration.

## 11. Production-backend residency and transfer

Host, numpy, float64. **No GPU work of any kind was run by this lane, and no authorization was sought
or held.** The ratchet is elementwise over sites and is a natural device kernel; the scatter back into
an owner goes through that owner's material-point sinks one site at a time, which is host-side and is
a property of the existing endpoint contract rather than of this module.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered, because no reference file was read. All prose in the
new module is written for this entry.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** Measured this session (`/Users/sw1/miniconda3/envs/aleph/bin/python`, darwin, float64, 2026-07-31): **65 passed** for `tests/vertical/test_connectors_protrusion.py` in 0.57 s, and **82 passed** for `tests/vertical/test_connectors_protrusion.py` together with `tests/ports/`. Two `tests/ports/` failures seen mid-session belonged to sibling lanes and were cleared by their owners rather than by this one — see §13b. Three mutants were applied to the shipped module and all three were killed; the table is in §13a. Status stays `PROPOSED`: the entry has had no PI review, four of the six connectors reuse classes whose own entries are also `PROPOSED`, and no connector here has ever been stepped inside an assembled world. |
| Reviewer | Agent-proposed. Unratified. No PI review. |
| Rollback | Delete `aleph/vertical/connectors_protrusion.py` and `tests/vertical/test_connectors_protrusion.py`, and remove the six `Binding` rows from `aleph/vertical/wiring.py`. Nothing else imports either file; no existing module was modified by this lane. |

### 13a. Mutation testing

Three mutants applied to `aleph/vertical/connectors_protrusion.py`, controls run against each, then
restored and verified byte-identical with `git diff`.

| # | Mutant | Tests killed | Survived? |
|---|---|---|---|
| 1 | `stall_force_pn`: `(kT/d) ln(1/e)` → `(kT/d) / e` — the same force scale, wrong dimensionless factor | 5 | no |
| 2 | `velocity_um_per_s`: `exp(-f d / kT)` → `exp(-f d / (2 kT))` — the factor-of-two slip in the exponent | 2 | no |
| 3 | the carrier count dropped from the delivered load: `n_i * f(v_i)` → `f(v_i)`, in both `site_load_pn` and `evaluate_sites` | 3 | no |

The module was restored after each mutant and verified **byte-identical** by `sha256`
(`17d98107…8fa527b` before and after) and by `git status`.

Mutant 1 is the one worth reading: it keeps `kT/d` — the whole force scale — and changes only the
dimensionless factor, so every dimensional check and every closure check still passes. It is killed
by the closed-form oracle, by the thermodynamic identity, and by the depolymerisation limit.

**Mutant 3 exposed a real gap on its first run and the gap was closed rather than reported around.**
`stall_force_pn` multiplies one barbed end's stall by the total carrier count, while `site_load_pn`
multiplies each site's per-barbed-end load by that site's count — two code paths for one fact.
Dropping the count from the *delivered* load therefore left the *reported* stall untouched, and on the
first run this mutant killed only two controls, neither of them
`test_the_bundle_stalls_at_n_times_the_single_filament_stall`, which is the control the mutant was
aimed at. `test_the_reported_stall_is_the_load_actually_delivered_at_zero_velocity` was written to pin
the two paths together and is what raised the kill count to three. **A ratchet that reports a stall
force it cannot deliver is exactly the failure this entry's oracle exists to catch, and until that
control existed the oracle only covered the single-site case.**

### 13b. Two transient failures in `tests/ports/`, neither this entry's, both since resolved

Recorded rather than dropped, because the reasoning is the part worth keeping. Three sibling lanes
landed `ALEPH-PORT-3301`, `-3303` and `-3304` into `ports/ledger/` in the same window as this entry,
and for part of this session two shared-file checks were red:

* `test_named_controls_resolve_to_real_tests` failed on **`ALEPH-PORT-3304`**, which cited three
  control tests not yet defined under `tests/`. Not this entry's file and not this lane's to edit.
* `test_index_is_not_stale` failed because `ports/ledger/INDEX.md` is generated from every entry on
  disk and four new entries had appeared. **Not regenerated here.** Doing so would have swept three
  other lanes' in-progress entries into a shared file this lane does not own — the precise move
  `docs/ACTIVE_SESSIONS.md` records as having cost this project two commits' worth of provenance in
  one afternoon.

Both were cleared by their owning lanes while this entry was being finished, and the final measured
run is **82 passed** across `tests/vertical/test_connectors_protrusion.py` and `tests/ports/`, with
`INDEX.md` now listing this entry. The reason for leaving them alone stands independently of that
outcome: waiting cost nothing and a regeneration could not have been undone.

### 13c. One guard is red until this file is staged, and it is right to be

`tests/vertical/test_package_exports.py::TestTheExportListMatchesThePackage::test_every_listed_name_is_a_module_git_actually_tracks`
fails while `aleph/vertical/connectors_protrusion.py` is untracked: commit `ceb267f` already lists
`connectors_protrusion` in the `aleph.vertical` package `__init__`'s `__all__`, and that guard
compares the list against `git ls-files` rather than against the directory, precisely so that a
declaration describing one machine's disk is caught before a fresh clone catches it.

**Not worked around.** The guard is doing exactly its job, `CLAUDE.md` §2 rule 4 forbids editing one
to unblock a lane, and that package `__init__` is a shared file this lane does not own. It clears
with `git add aleph/vertical/connectors_protrusion.py` by the session that commits. The rest of
`tests/vertical/` is **1,191 passed** and `tests/firewall/` is **26 passed**, both measured this
session.

## 14. Honest limits — what this entry does NOT establish

**No kinetics, anywhere in the six.** Nothing here binds, unbinds, caps, nucleates or turns over. The
seam and the root carry a `bound` flag and no law that changes it; the ratchet carries a velocity and
no law that advances the monomer count. So nothing in this entry is evidence about seam lifetime,
root detachment under load, capping, or the barbed-end population. Absent, not implied.

**No polymerisation rate is shipped, so no protrusion velocity is predicted.** The zero-load velocity
`v0` is a required argument with no default and no fallback. Every velocity this module reports is a
function of a number the caller supplied. `ALEPH-PORT-3205` §14 states the same limit from the owner's
side, and this module does not route around it.

**`k_B T` is `ASSUMED` and this module never supplies it.** Every stall force below is therefore an
`ASSUMED`-class quantity in its magnitude, and `ANALYTIC_ORACLE` only in its form.

**The membrane's curvature cost under a near-point load is not computed here and is not checked.**
The registry's distinction between a distributed pressure and a near-point load is implemented as
*where the force lands and how many barbed ends share it* — (R5), (R6) — and the curvature response
is `HelfrichMembrane`'s. Whether a single-site tip load produces the right local curvature is a
membrane control that **is not written**.

**The ratchet is evaluated at a frozen configuration inside a candidate**, like every other connector
here, so its `AdjointPair` reports zero endpoint work and zero active work over the step. The *rate*
is real and `active_power_pn_um_per_s()` reports it; booking a per-second quantity in a pN·µm slot
would be a unit error dressed as accounting. **There is therefore no control that the ratchet's energy
books close over a real step**, because no real step exists to close them over.

**Nothing here has been stepped inside an assembled world.** `assembly.py` contains none of these
owners. Six connectors that evaluate correctly in isolation is what is established; a protruding cell
is not.

**The seam and the root reuse a class whose own entry is `PROPOSED`**, and the two-sidedness argument
in §4c is this lane's reasoning about the registered mechanism, not a decision anybody ratified. If a
transient actin seam should be tension-only, the change is one argument to one builder and the control
`test_the_seam_and_the_root_carry_load_in_both_senses` is the thing that would have to be inverted —
which is where it should be, rather than buried.

**The motor projection assumes the owner's published `NMII_HEAD` sites are the right binding sites.**
`ProtrusionActinTarget` projects every published site of that role, in publication order. Which sites
an owner *should* publish — how many, how spaced, whether occupancy excludes neighbours — is the
owner's question and this module answers none of it. `nmii.py`'s own §14 records that the
`ActinBindingSites` shape is provisional and that reconciling it with the actin owners is future work;
this entry is the first use of it against a real owner and does not close that item.

**`aleph/vertical/wiring.py` was not written by this lane** and the connector count in it is not this
entry's to change.
