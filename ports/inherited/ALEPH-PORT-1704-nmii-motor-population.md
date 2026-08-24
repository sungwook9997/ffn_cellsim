# ALEPH-PORT-1704 — nmii: a head-resolved motor population, never an aggregate tension field

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1704` |
| Lane | `L17 registry — actomyosin load-path compartments` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `CONCEPT_ONLY` — **verdict: RE-DERIVE the code, and this is the one compartment where a future PORT is worth re-examining.** Nothing crossed today. |

---

## 1. Aleph API

```python
from aleph.state.census_actomyosin import NMII, MOTOR_CONNECTOR_RULE
```

Aleph target file: `aleph/state/census_actomyosin.py`. This entry authorises the `nmii` contract entry
and the motor-connector rule constant. No motor kinetics, no crossbridge law and no force scale is
authorised.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/nmii_actuator.py`, `ffn_sim/ac/motor/segment_motor.py`, `ffn_sim/ac/motor/hand.py`, `ffn_sim/ac/motor/powerstroke_analytic.py`, `ffn_sim/ac/motor/hill_fv_analytic.py`, `ffn_sim/ac/motor/minifilament_topology.py`, `ffn_sim/ac/motor/bell_kinetics_analytic.py` |
| Source symbol(s) | `NMIIActuatorStateOwner`, `NMIIActuator`, `BackboneArmMechanics`, `NMIIBendingStiffness`, `SegmentMotorConnectorRuntime`, `bell_off_rate`, `catch_slip_off_rate`, `hill_velocity`, `hill_force`, `crossbridge_force`, `stall_abscissa`, `working_stroke_strain`, `head_newton_residual`, `MinifilamentTopology`, `straddle_frame` |
| Read from | **`git show be0e5876:<path>`** — the immutable commit. |
| Working tree == commit? | `yes` — `git diff be0e5876 -- ffn_sim/ac/engine/nmii_actuator.py` and `-- ffn_sim/ac/motor/segment_motor.py` are both empty. |

| Path | `sha256:` | Lines |
|---|---|---|
| `ffn_sim/ac/engine/nmii_actuator.py` | `sha256:1523ed040a259971` | 1169 |
| `ffn_sim/ac/motor/segment_motor.py` | `sha256:95a004d55448a71c` | 904 |

## 3. Why source-derived porting beats clean-room

This is the one compartment of the five where the honest answer is not a flat RE-DERIVE, so it is worth
being precise.

**What is genuinely there.** Head-resolved attachment and detachment with both Bell and catch–slip
laws, a working stroke as a strain offset with a Newton residual on the head, a Hill force–velocity
pair with its power and step-work companions, a bipolar backbone with a derived bending stiffness rather
than a magic number, and — the part that matters most — a **two-sided crossbridge scatter** that writes
into two separately-owned force arrays. It is the only compartment with native evidence at
non-trivial population, and its own audit calls its SF motor connector the best-evidenced connector in
that engine. Its bound-head fraction is reported as *emergent from attachment events* rather than
imposed, which is the right shape.

**Why nothing crossed anyway, today.** Three reasons, in descending order of force:

1. **The force scale is split-brained, and the lane that produced every banked number is the one
   without guards.** The engine layer refuses to default the crossbridge stiffness, stall force,
   unloaded velocity, bending rigidity and attachment rate. A second assembly path defaults **all of
   them, unguarded**, and its own audit records that every banked headline number came through that
   path. So the code is disciplined and the results are not, and a port that took the code without the
   guards would inherit the appearance of rigour.
2. **The headline γ is not the power stroke.** Its own audit records the crossbridge contribution as
   `0.0191` against a network total of `3.7043` — the active term is **0.51%** of the quoted number.
   Whatever that number measures, it is not predominantly motor activity.
3. **Every step was force-accepted**, its acceptance array set to all-ones, with no reference row ever
   measured. An accepted-step transaction whose predicate always returns true is a transaction in name.
   And Aleph's acceptance predicate is `ALEPH-DQ-104`, reserved for the PI — so porting anything whose
   behaviour is entangled with an always-accept predicate would pre-empt that decision.

There is also a provenance defect worth recording because it is the shape this project exists to
prevent: a run stamped its own provenance as source-grounded while supplying a backbone persistence
length explicitly labelled a non-production diagnostic fixture. Its own audit's summary is the right
one — the guard works and the ledger lies.

**Verdict: RE-DERIVE the code. Revisit for a genuine PORT when a motor vertical is authorised**, at
which point `segment_motor.py`'s two-sided scatter and `hand.py`'s catch–slip branch structure are the
two things worth a fresh, narrowly-scoped entry. A branch structure whose case analysis is easy to get
subtly wrong is exactly what `ports/TEMPLATE.md` §3 admits as legitimate. That entry is not this one and
must not inherit this one's audit.

## 4. Physical or mathematical law represented

No kinetics was ported. Two declarative laws, derived here.

**A motor population is an owner; an active tension field is not.** Consider the two candidate
representations of myosin activity acting on an actin owner `a`.

*As a connector.* NMII owns heads. Head `h` attached to actin material point `x` at strain `ε` exerts
`f` on `a` and `−f` on the minifilament backbone. Two accumulations, two owners, `Σf_a + Σf_nmii = 0` as
a **check**. Momentum has a source and a sink, both named. Removing the connector removes exactly that
active load path.

*As a field.* An aggregate active tension `σ_a(x)` is added into `a`'s force array. Then:

- there is no second accumulation, so closure is not a check — there is nothing for the force to
  balance against;
- momentum enters the system from nothing, so the assembly's net force no longer vanishes and the one
  global control that would catch it has been disabled by construction;
- it cannot be switched off by removing a connector, because it never was one, so the standard negative
  control is unavailable.

So "NMII acts through a separate MOTOR connector" is the condition under which active loading is
falsifiable. It is not a preference about code organisation. `MOTOR_CONNECTOR_RULE` states it and
`validate_group` requires it on the entry.

**One connector per actin owner, and the count is forced by disjointness.** Since cortex, `sf_arc`,
lamellipodium and filopodium own disjoint populations, and a connector's endpoints are two named
owners, NMII needs a distinct MOTOR connector for each — four, in this registry. This is not
duplication: each is separately schedulable, so removing one isolates exactly one active load path,
which is what makes a per-compartment control possible at all. A single "NMII acts on actin" connector
would be ill-typed, because "actin" is not an owner.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| head strain | µm | m | signed |
| crossbridge stiffness | pN/µm | N/m | > 0, **unsourced** |
| stall force | pN | N | > 0, **unsourced** |
| unloaded sliding velocity | µm/s | m/s | > 0, **unsourced** |
| working-stroke displacement | µm | m | > 0, **unsourced** |
| attachment / detachment rate | 1/s | 1/s | ≥ 0, **unsourced** |
| bound-head fraction | 1 | 1 | [0, 1], an **output** |

Singular cases:

- **Zero bound heads.** Legal, and the transmitted active load must be exactly zero. A minifilament
  whose heads are all unbound still experiences cytosolic drag, so it is not inert — which is why drag
  is a separate connector and not a term inside the motor connector.
- **A bound-head fraction supplied as an input.** Refused. It is an output of the rates; reproducing an
  observed fraction by construction and then citing the agreement is circular, and
  `NMII.unsupported_claims` names that specific circularity.
- **An aggregate contractility parameter requested.** Refused. The component does not offer one, and
  `NMII.unsupported_claims` says so explicitly rather than leaving the absence to be discovered.
- **A per-isoform mechanical conclusion.** Refused. Isoform labels exist; per-isoform kinetics do not.

Invariants:

- **I1.** `nmii` records the MOTOR-connector rule verbatim.
  `tests/state/test_census_actomyosin.py::test_nmii_reaches_actin_only_through_a_separate_motor_connector`.
- **I2.** It owns no actin filament population, and says so verbatim.
  `test_nmii_owns_no_actin_filament_population`.
- **I3.** Isoform, activation and spatial-population labels are its own state, so they are
  subpopulations rather than fields. `test_nmii_owns_head_resolved_subpopulations_and_not_a_tension_field`.
- **I4.** Its numerical representation names "aggregate tension field" as the thing it is not.
  Same test.
- **I5.** Head and minifilament counts and the bound fraction are declared outputs, with five rate
  priors named as owed. `test_each_dynamic_population_records_that_its_count_is_an_output`.
- **I6.** No state key it owns is claimed by another owner.
  `test_no_state_key_is_claimed_by_two_owners`.

## 6. Source evidence class and known retractions

I looked in that repository's audit directory, its gap-evidence cards, and its rolling roadmap, and
there is more recorded against this compartment than any other.

- Its component audit rates `nmii` at its **second-highest execution rung, force-real**, with a real
  native population and an emergent bound fraction. Its SF motor connector is called the
  best-evidenced connector in that engine.
- The same audit records that on the lane that ran, **both the actuator and its port were the same
  position and force arrays** — so the "two never-merged arrays" claim in its own docstring is false on
  the lane that produced the numbers. That is the single most important entry in this section, and it is
  its own audit's finding, not mine.
- The crossbridge contribution is `0.0191` of a network total of `3.7043`.
- Every step was force-accepted; no reference row was measured.
- Its gap cards record the SF NMII magnitude as `SEAMED` and the traction magnitude as **`INVALID`**.
- Its roadmap records the motor connector as mechanism-passing and **magnitudes solver-blocked**.
- Its own audit records **no resolvable DOI** for the catch–slip law's default parameterisation, under
  both the motor and adhesion uses.

Reachability: **live**, genuinely. Constructed on non-test paths and exercised by a non-test driver.
This is the only one of the five compartments where that is true.

Tests: 901 lines, 24 tests, 64 doubles — but **8 real device allocations**, and several tests touch
physics: a derived-not-magic bending stiffness, real-kernel binding with fidelity and count-drift
checks, live geometry kernels, the two-sided load scatter, acceptance-predicated kinetics, and a
bit-exact snapshot/rollback. Plus separate analytic oracle suites for the Bell, Hill and power-stroke
laws. So the tests are the best of the five and still leave the magnitudes untested, because the
magnitudes have no sources.

**Retraction in substance, found:** the provenance stamp described in §3 — source-grounded claimed while
a diagnostic fixture was supplied — and the `INVALID` traction magnitude. Both are that project's own
records.

## 7. Independent oracle or derivation

For this entry's declarative content: an exact statement. Under the connector representation the
assembly's net force vanishes identically; under the field representation it does not. So "is NMII a
connector or a field" is decidable from a net-force measurement on an isolated assembly, and Aleph
already owns that control — `test_isolated_closed_surface_has_no_net_force` in the vertical measured
`1.46e-13 pN` against a force scale of `1.95e4 pN`. This entry does not add a new oracle; it records
that the existing one is the discriminator, which is better than inventing one.

For the kinetics a later lane would write, the oracles available to Aleph without the reference open
are all closed forms: detailed balance on a two-state attachment cycle at zero load fixes the ratio of
rates exactly; the Bell law reduces to a bare rate at zero force; the Hill relation has exact endpoints
at zero load and at stall; and the working stroke gives an exact energy-per-cycle bound. Four exact
endpoint checks, and none needs a comparison against that repository — which matters, because
`ports/TEMPLATE.md` §7 is explicit that agreement with it would not be acceptance.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_actomyosin.py::test_nmii_owns_head_resolved_subpopulations_and_not_a_tension_field` | Six head- and population-resolved state blocks are owned, and the representation names the aggregate field as excluded. |
| Positive | `tests/state/test_census_actomyosin.py::test_nmii_reaches_actin_only_through_a_separate_motor_connector` | The MOTOR-connector rule is present verbatim. |
| Positive | `tests/state/test_census_actomyosin.py::test_nmii_owns_no_actin_filament_population` | It owns no `*_filament_population` block and is absent from the filament-owner list. |
| Positive | `tests/state/test_census_actomyosin.py::test_each_dynamic_population_names_the_rate_priors_evidence_is_owed_on` | Five rate priors are named as owed. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_nmii_without_the_motor_connector_rule_is_refused` | Replacing `nmii.mechanical_role` with only the no-filament rule is refused, because the entry would then no longer exclude an aggregate tension field. This is the control that guards §4. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_a_participant_that_drops_the_no_filament_rule_is_refused` | A non-filament participant whose role omits the no-ownership rule is refused. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_an_entry_with_no_stated_limits_is_refused` | Emptying `nmii.unsupported_claims` is refused — the four disclaimers, including the circular bound-fraction one, are mandatory. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_dropping_the_count_is_an_output_rule_is_refused` | Removing the count rule is refused. |

## 10. Numerical and precision envelope

No arithmetic crosses this boundary; controls are exact string and set comparisons and no tolerance is
stated, because the properties are discrete.

Two precision facts are recorded because they bear on a future motor lane and would otherwise have to
be rediscovered. First, a head-level kinetic scheme accumulates a very large number of small force
contributions, so the accumulation precision matters more than the working precision — `float32`
compute with `float64` accumulation is the configuration `ALEPH-DQ-107` proposes, and a motor lane is
the place where the accumulation half earns its keep. Second, an active connector's power is a signed
sum, and a signed sum is exactly what a max-norm residual gate cannot see: PLAN §11 established that
for the dilational virial, and the same argument applies verbatim to active power. A motor lane must
not gate on a max residual.

## 11. Production-backend residency and transfer

Host-side declarations here; no device residency, no transfer, not in the step loop.

The residency note this entry leaves for a motor lane is the substantive one. A two-sided crossbridge
scatter requires **both** force accumulators resident simultaneously, and the head-to-actin binding map
changes on accepted steps only — so the map is device-resident state that may be mutated exactly once
per accepted step, and never during a rejected trial. That is the coupling between residency and the
transaction, and it is where the always-accept defect in §3 did its damage: with every step accepted,
the distinction between "mutate on accept" and "mutate always" never showed up.

## 12. Comments and docstrings to discard

Read, and **not** carried: every class, kernel, runtime and analytic-function name; the module paths;
the evidence-rung vocabulary; the gap-card identifiers; the roadmap commit references; the catch–slip
default parameterisation and its unresolvable citation; the persistence-length fixture label; every
force-scale constant.

Special mention, because it is the most instructive thing I read: a docstring there asserts two
never-merged arrays while the lane that ran used one. **That is a docstring making a claim its own
runtime contradicts**, and PLAN §9 records me being misled by exactly that class of prose inside
Aleph's own tree earlier in this project. It is the reason §12 of the template exists. What replaces it:
§4 above, which derives *why* two arrays are required from the falsifiability of force closure, so the
claim is checkable rather than asserted — and Aleph's existing net-force control, which is the check.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** `PROPOSED`. The contract entry and its controls pass (`tests/state/test_census_actomyosin.py`, 88 passed). The substantive claims are the connector-versus-field derivation in §4 and the RE-DERIVE-for-now verdict in §3, neither reviewed by a human. |
| Reviewer | agent-proposed, **unratified**. |
| Rollback | Delete the `NMII` entry, `MOTOR_CONNECTOR_RULE`, and their tests. `validate_group` then fails on the missing participant and on the missing motor rule. Nothing outside the module depends on them. |

## 14. Honest limits

- **Unratified.**
- **No motor kinetics exists in Aleph.** This is a declaration. Nothing here is evidence that Aleph
  represents NMII.
- **Every magnitude is `UNSOURCED`** and no literature was read. Five priors are named as owed and none
  is supplied.
- **This entry's verdict is explicitly provisional in one direction.** I judged RE-DERIVE partly because
  no motor vertical is authorised, which is a scope judgement rather than a quality judgement. If a
  motor vertical is authorised, the two-sided scatter and the catch–slip branch structure deserve a
  fresh entry, and it must not inherit this verdict.
- **The `0.51%` and the always-accept findings are that project's own audit, read, not measured by me.**
  I did not run its suite, launch its kernels, or reproduce its γ. Evidence class `AUDIT_READ`.
- **I did not verify the claim that its guarded and unguarded paths disagree** by running both. I read
  the audit that says so and read enough of both files to find it plausible. That is weaker than
  measurement and it is one of the two loads this entry's §3 verdict rests on.
- **Nothing here addresses whether four MOTOR connectors is the right number.** It follows from the
  registry's four actin owners. If a fifth actin owner is added — a lamella, say, under the re-entry
  condition recorded on the `LAMELLA` entry in `aleph/state/census_actomyosin.py` — a fifth motor
  connector follows, and nothing currently checks that implication.
