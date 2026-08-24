# ALEPH-PORT-3621 — the whole cell: every constructible compartment in one world

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3621` |
| Lane | `46143f30` Lane W (the whole cell, for ablation) |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` — nothing is ported by this entry; it composes Aleph's own owners |
| Aleph target | `aleph/scenarios/whole_cell.py` (new) |
| Exists because | Thirteen compartments declare state and `python -m aleph.scenarios.audit` reports 13/13 implemented. `build_vertical` assembles **three** of them, and the two scenarios that exist — `indent.py`, `spread.py` — use only those three. So the project has thirteen compartments and no world that contains them. The PI's method is **ablation**: start from a complete cell, remove one compartment at a time, and measure whether the answer moves. That method has no baseline until a complete cell exists, and an ablation measured against an incomplete baseline attributes to the removed compartment whatever the *silently missing* one was already doing. |
| What this entry does **not** claim | That the assembled world is a cell, that it is coupled, or that any number it produces is physics. See §14.1 and §14.2, which are the two sentences a reader in a hurry should take. |

---

## 1. Aleph API

```python
# CREATED — authorised by this entry, and nothing else is written by this lane
aleph/scenarios/whole_cell.py
tests/scenarios/test_whole_cell.py
```

```python
# CONSUMED, READ-ONLY — every one of these is another lane's file and none is edited
from aleph.scenarios.world import (
    CellWorld, ConnectorSlot, OwnerSlot, WorldCompositionError, build_cell_world,
    site_field_resolver,
)
from aleph.scenarios.audit import COMPARTMENT_MODULES      # the registry of the thirteen
from aleph.vertical.assembly import build_vertical         # membrane, cortex, osmotic envelope
from aleph.vertical.connectors import SiteField
from aleph.vertical.relax import relax_to_equilibrium
from aleph.runtime.participant import Phase
# and the ten remaining owner modules under aleph/vertical/, constructors only
```

**The list of thirteen is imported, never retyped.** `COMPARTMENT_MODULES` in
`aleph/scenarios/audit.py` is the registry, and `build_whole_cell` derives its census from
`COMPARTMENT_MODULES.keys()`. A hand-written list in this module would be a second registry, and the
first time somebody added a fourteenth compartment the census would report 13/13 complete while
missing it — the failure mode is *silence*, which is the same one the whole entry exists against.
`test_the_compartment_list_is_the_registrys_and_not_this_modules` grades it.

**No compartment module is modified to make it compose.** That is the load-bearing constraint on
this lane. A compartment that will not assemble, or that assembles and cannot be scheduled, stays
**out of the world with its reason carried as data** — see §4 — rather than being edited, wrapped,
stubbed, or given a shim that makes the count reach thirteen.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` — cited for the record only |
| Source path | none consulted for this entry |
| What was taken | **Nothing.** No line, no identifier, no constant, no parameter. Every object composed here is Aleph's own, and the composition machinery (`build_cell_world`) is Aleph's own from `ALEPH-PORT-3404`. |

The one fact carried forward from the reference project is a **negative** one, and it is already in
`PLAN.md` §7: that project shipped a whole-cell solve with a declaration and **zero production
tasks**. This entry is written so that the same thing cannot happen quietly here — §4's census and
§8's stepping control are both aimed at exactly that failure.

## 3. Why source-derived porting beats clean-room

It does not, here. Nothing was ported; the class is `RE-DERIVED` and the entry is a **composition**
entry. The rule that applies instead is `ALEPH-PORT-3613` §3's: a composition entry earns its keep
only if it could have come back negative. This one did, twice, and both are recorded rather than
smoothed over — one compartment does not compose (§4b) and eight owner/phase pairs do not schedule
(§4c). An entry whose only possible outcome was "13 of 13, complete" would not have been worth
writing.

## 4. Physical or mathematical law represented

**None.** No law is implemented, derived or claimed by this entry. What it establishes is a
**bookkeeping invariant**, and the invariant is the deliverable:

> **I1 — the census is total and disjoint.** For every world this module builds,
> `set(present) | set(absent.keys()) == set(COMPARTMENT_MODULES)` and
> `set(present) & set(absent.keys()) == {}`. Every one of the thirteen is either in the world or
> carries a written reason for not being, and none is both or neither.

That is a weak-looking claim and it is the entire point. An ablation study measures a *difference*
between a baseline and a knockout. If the baseline is silently missing a compartment, the difference
is attributed to the compartment that was removed on purpose, and there is no observable anywhere in
the experiment that would say otherwise. `aleph/state/full_census.py` records the same trap one level
up: an empty manifest finalizes cleanly while the real one cannot, so it looks *healthier*.

### 4a. Why the world is composed and not written

`build_cell_world` (`ALEPH-PORT-3404`) already refuses every structurally wrong configuration —
an empty world, a duplicate owner name, a connector whose endpoint is absent, a connector registered
and never bound, a coupling whose force lands nowhere. Reimplementing any of that here would be a
second place for those refusals to be right. This module supplies **owners, phases, energies, loads
and a census**; it supplies no acceptance rule, no descent, and no refusal that `world.py` already
makes.

### 4b. The compartment that does not compose, and why it is not repaired here

**`cytosol` is refused, and the reason is measured rather than inferred.** `CytosolField` constructs
(an 8³ grid at `dx=0.5` with a nucleus inclusion gives 1,308 internal and 408 boundary faces), and it
computes: `ACCUMULATE_FIELD` reports a divergence-theorem residual, `SOLVE` reports one monolithic
backward-Euler step of the coupled Biot system, `EVALUATE_GATES` reports a constituent Darcy flux.
All three return non-empty `WorkReceipt`s.

All three also advance `ctx.witness_count` by **exactly 0**, because the module computes with NumPy
directly rather than through `ctx.backend`. The witness is `backend.op_count + rng.draw_count`
(`aleph/runtime/participant.py`), so the coverage gate cannot see any of that work. Scheduling it
raises, at `propose` time, before any predicate runs:

```
aleph.runtime.pipeline.CoverageViolation: coverage FAILED: 3 finding(s);
2/3 required participants evaluated
  [no_witness] participant 'cytosol' in phase accumulate_field: ... (witness delta 0)
  [no_witness] participant 'cytosol' in phase evaluate_gates:   ... (witness delta 0)
  [no_witness] participant 'cytosol' in phase solve:            ... (witness delta 0)
```

A world containing `cytosol` therefore takes **no step at all** — not a rejected step, no step. The
pipeline's own message names the remedy — *"register it EXCLUDED with a reason instead of scheduling
a stub"* — and `absent["cytosol"]` is that remedy one level up, where `build_cell_world` has no
`EXCLUDED` channel to offer.

**This is not a defect report against `cytosol` and this lane does not repair it.** Routing its
solve through `ctx.backend` is a change to `aleph/vertical/cytosol.py`, which this lane does not own
(`CLAUDE.md` §1), and the choice between *"route the Biot solve through the backend"* and *"give the
pipeline a way to admit a host-only field owner"* is an architecture decision above this level
(`CLAUDE.md` §2 rule 5). Both are named so the PI chooses from a list; neither is taken.

### 4c. The eight owner/phase pairs that are scheduled nowhere, and are recorded

The same witness rule removes **phases**, not only compartments, and a phase dropped in silence is
the same class of failure as a compartment dropped in silence. Measured, one owner at a time, on a
fresh context:

| Owner | Phase | Witness delta | What the receipt said it did |
|---|---|---|---|
| `cortex` | `DISCOVER` | 0 | steric pair list rebuilt: 0 cross-filament pairs in range |
| `filopodium` | `DISCOVER` | 0 | crosslinker topology resolved: 1 spacing element over 3 filaments |
| `lamellipodium` | `DISCOVER` | 0 | population: 3 filaments, 12 nodes, 1 branch junction |
| `ecm` | `DISCOVER` | 0 | 9 committed crosslinks resolved at topology epoch 0 |
| `microtubule` | `DISCOVER` | 0 | topology rebuilt: 8 segments, 6 bending stencils |
| `extracellular_medium` | `DISCOVER` | 0 | mobility factorisation |
| `ecm` | `PROPOSE_KINETICS` | 0 | damage proposed on 0 of 24 segments |
| `focal_adhesion` | `PROPOSE_KINETICS` | 0 | 0 clutch transitions queued |

Every one is a **pure topology or bookkeeping pass in Python**, which is why it launches nothing.
Each is carried in `WholeCell.unscheduled_phases` with the measurement, and
`test_the_phases_that_witness_nothing_are_recorded_rather_than_silently_dropped` grades that the
record is non-empty and names the gate. Dropping `DISCOVER` costs nothing for a static baseline —
none of the five topologies changes when nothing is inserted or removed — and that is a claim about
*this* fixture, not about those owners.

**`extracellular_medium` is the one worth stating separately**: its `ACCUMULATE_FIELD` was measured
to work with **no** prior `DISCOVER` at all (witness delta 3, 4 declared ops), so dropping `DISCOVER`
does not silently disarm it. That was checked by running rather than assumed, because "the
factorisation is lazy" and "the field is now zero" produce the same absence of an exception.

### 4d. What ablation means in this module, and what it does not

`build_whole_cell(without=("nucleus",))` removes an owner from the world. It does **not** remove a
term from a coupled system, because ten of the twelve present compartments are not coupled to
anything (§14.2). Two consequences are implemented rather than left to the caller, and both are
recorded in `WholeCell.withdrawn` — a **separate** channel from `absent`, so `absent`'s keys stay
compartments and only compartments and I1 can remain an equality against the registry:

* **A connector whose endpoint is ablated is withdrawn**, with the ablation that caused it named in
  the reason. Leaving it would make `build_cell_world` refuse the world outright — which is correct
  behaviour from `world.py` and an unhelpful error from here.
* **A field owner whose load target is ablated is withdrawn**, on the same terms. The osmotic
  envelope loads the membrane and nothing else; without a membrane its load "lands nowhere", which
  `CellWorld.evaluate_forces` refuses by design.

Both are *derived* removals and both say so in their reason text. A derived removal that read like a
caller's choice would be a silent second ablation inside the first.

### 4e. Ablating either coupling endpoint leaves a world that cannot step — MEASURED

Found by running, not predicted. Both connectors in this world name `membrane` and `cortex`, so
ablating **either** leaves a world with no coupling at all. Such a world composes, reports a finite
energy and a finite residual, and **accepts no step**:

```
FAIL balance.connectors_recorded  no connector recorded an adjoint pair, so the balance gate has
                                  nothing to check. An empty ledger is not a balanced one.
```

`BalancePredicate` is right, and `relax_to_equilibrium` composes it with the coverage and descent
gates, so every candidate is rejected. The rejection is step-size independent, so the controller
breaks out immediately and the report reads `accepted=0, rejected=1` rather than looking like a
stiff configuration ground down to `min_step`.

**It is not routed around.** Handing the driver a predicate with the balance gate removed would step
happily and would be measuring a world whose books nobody closed — `CLAUDE.md` §2 rule 4. Instead
`WholeCell.steppable` and `WholeCell.unsteppable_because` carry it as data, and
`test_a_world_with_no_coupling_is_reported_unsteppable_and_the_balance_gate_is_why` drives it in
both directions.

This is the strongest available argument for §14.2: a world with twelve compartments and two
couplings has a **single point of failure for steppability**, and the ablation method is what
surfaced it.

## 5. Units, domains, singular cases, invariants

| Quantity | Unit | Domain |
|---|---|---|
| length | µm | the fixture spans ≈ 0 – 6 µm |
| force | pN | |
| energy | pN·µm | must be finite; the descent requires it non-increasing |
| descent step | dimensionless (`relax_to_equilibrium`'s `dt` is a step size, not a time) | > 0 |
| owner count | owners | ≥ 1, refused at 0 by `build_cell_world` |
| compartment census | compartments | `len(present) + len(absent) == 13`, always |

**Singular cases, each handled explicitly rather than by an exception nobody reads:**

* **`without=` naming something that is not a compartment.** Refused with `WholeCellError`. A typo
  that silently ablated nothing would make an ablation study report "removing `nucelus` changed
  nothing", which is true and worthless.
* **`without=` naming a compartment that is already absent.** Accepted, and the reason records
  *both* facts, because "you removed it" and "it was never going to be there" are different and the
  second is the one that matters to the reader.
* **Ablating every compartment.** `build_cell_world` refuses an empty world; the error propagates
  unchanged.
* **Ablating every owner that reports `internal_forces`.** `CellWorld.evaluate_forces` refuses;
  the error propagates unchanged. Both refusals belong to `world.py` and are not re-implemented.
* **An ablation that leaves the world with no coupling.** Composes, reports finite numbers, and
  accepts no step — §4e. Carried as `steppable=False` with a written reason rather than raised,
  because the world is still a legitimate object to hold and read; what it is not is drivable.
* **The pace owner.** `relax_to_equilibrium` reads `world.membrane.mobility_um_per_pn_s`. If the
  membrane is ablated, the first present owner carrying a mobility is used and **named** in the
  report, rather than the world failing at step time on an attribute error.

## 6. Source evidence class and known retractions

1. **Every parameter in this module's fixtures is `UNVERIFIED` and most are `ASSUMED`.** The
   material cards, geometries, node counts and populations are assembly fixtures chosen so the
   owners construct; not one was read from a primary source by this lane. Where a module ships a
   default (`default_microtubule_card`, `default_material_cards`, `default_lamellipodium_card`) that
   default is used, and it carries whatever evidence its own ledger entry gives it — no more.
2. **The geometry is not anatomy.** The compartments are placed where their constructors put them.
   The nucleus is inside the membrane by arithmetic (radius 3 µm inside 5 µm) and nothing else is
   positioned with respect to anything. There is no collision detection between owners and no claim
   that the arrangement resembles a cell.
3. **`intermediate_filament`'s five parameters per population have `citation=None`** and the module
   deliberately ships no default card. They are placeholders and are marked as such at the call
   site.
4. **`nmii`'s target is a `RigidActinTrack`, which is an instrument and not a model of actin.** It
   exists so the motor is not inert. A real assembly hands `register_target` an `ActinBindingSites`
   published by an actin owner; this one does not, and therefore the motor pulls against a rigid
   fiction.
5. **Retracted before it is made:** any reading of "12 of 13 compartments" as "a whole cell was
   simulated". Twelve owners were composed, scheduled, and stepped together. §14 says what that is
   and is not.

## 7. Independent oracle or derivation

Three, and they are independent of one another:

* **The verified vertical.** Reducing the whole cell to its membrane and cortex — `without=` the
  other eleven — must reproduce `build_vertical(subdivision_level=1)`'s own total potential energy,
  because it then contains exactly the same three owners and the same two connectors. This compares
  the composition against the one physics result the project has verified, and it is the strongest
  control in the file. It is `test_the_vertical_trio_alone_reproduces_the_verified_verticals_energy`.
* **`aleph/scenarios/audit.py`'s registry** is the oracle for the *list*: the census is derived from
  `COMPARTMENT_MODULES`, so this module cannot agree with itself about what thirteen means.
* **The pipeline's coverage verdict** is the oracle for *stepping*: it is computed from the witness
  delta per slot, which a stub cannot fabricate without launching work. A world that assembles,
  steps and reports a plausible energy while a slot witnessed nothing is exactly the reference
  project's whole-cell solve, and the verdict is the only thing that distinguishes them.

## 8. Positive control

All in `tests/scenarios/test_whole_cell.py`.

| Control | Asserts |
|---|---|
| `test_the_census_is_total_and_disjoint` | I1 — present ∪ absent is the thirteen, present ∩ absent is empty, for the baseline and for three ablations. |
| `test_the_compartment_list_is_the_registrys_and_not_this_modules` | §1 — the census keys equal `COMPARTMENT_MODULES`'s keys, and the module's own text does not carry a second hard-coded list of thirteen names. |
| `test_every_absent_compartment_carries_a_non_empty_reason` | §4 — an absence with no reason is indistinguishable from an oversight. |
| `test_the_whole_cell_takes_an_accepted_step_and_the_clock_advances` | §4 — at least one **accepted** step, a finite energy, and a clock and step index that moved. Assembling without stepping proves only that constructors run. |
| `test_the_coverage_verdict_is_complete_after_the_descent` | §7 — every scheduled owner and connector witnessed work. |
| `test_the_coverage_verdict_is_incomplete_before_anything_has_run` | §7 — the other direction, so `ok` cannot be a constant. |
| `test_the_vertical_trio_alone_reproduces_the_verified_verticals_energy` | §7 — the composition against the verified vertical, to within the `fsum`-versus-`plus` few-ULP difference. |
| `test_removing_one_compartment_changes_the_owner_count_and_the_world_still_steps` | §4d — the ablation contract, for `nucleus`. |
| `test_every_present_compartment_can_be_removed_singly_and_the_world_still_builds` | §4d — every one of the twelve, one at a time, each still composing and reporting a finite energy. |
| `test_removing_the_cortex_withdraws_the_couplings_that_touch_it_and_records_them` | §4d — a derived absence is recorded, not silent. |
| `test_removing_the_membrane_withdraws_the_osmotic_envelope_that_loads_it` | §4d — the load-lands-nowhere case, recorded rather than raised. |
| `test_the_phases_that_witness_nothing_are_recorded_rather_than_silently_dropped` | §4c — the eight owner/phase pairs are in the report with the gate that removed them, **and** the world does not schedule what the record says is unscheduled. Prose beside a schedule that disagreed with it would be worse than no record. |
| `test_every_scheduled_phase_is_one_the_owner_declares` | §4c — driven rather than reasoned: an owner scheduled into a phase it declares no work for raises out of `accumulate` during `propose`, so a world that steps has proved its own schedule. |
| `test_the_report_names_the_owner_that_paces_the_descent` | §5 — the pace owner is named rather than assumed to be the membrane, and a membrane-less world's pace owner carries a mobility. |
| `test_the_world_contains_exactly_what_the_census_says` | §4 — `present` is a claim about the world and `owner_names()` is the world; the only permitted difference is the osmotic envelope, named rather than tolerated as slack. |
| `test_the_refusal_table_and_the_builders_do_not_disagree` | §4b — everything refused must still be **buildable**, or the counterexample cannot be driven and the recorded reason becomes unfalsifiable. |

In the checkable form the port-discipline parser reads:

| Positive | `tests/scenarios/test_whole_cell.py::test_the_census_is_total_and_disjoint` | Every one of the thirteen is present or absent-with-a-reason, and never both. |
| Positive | `tests/scenarios/test_whole_cell.py::test_the_whole_cell_takes_an_accepted_step_and_the_clock_advances` | The composed world takes at least one accepted step; energy finite, clock advanced. |
| Positive | `tests/scenarios/test_whole_cell.py::test_the_vertical_trio_alone_reproduces_the_verified_verticals_energy` | The composition reduced to the vertical reproduces the vertical's own energy. |

## 9. Deliberately failing negative control

| Control (must fail) | Asserts |
|---|---|
| `test_nothing_is_refused_any_more_and_the_machinery_that_drove_a_refusal_still_works` | §4b — the counterexample. `schedule_refused=("cytosol",)` puts it in the world; the first `propose` raises `CoverageViolation` naming `no_witness` in all three of its phases. **The reason recorded in `absent` is therefore checkable rather than asserted**, which is the difference between a documented exclusion and a rumour.  **Renamed 2026-08-06:** the control was called the name it carried while `REFUSED` still had an entry in it until `REFUSED` became empty on 2026-08-04; the machinery it drives is unchanged and the test now says so in its name. The entry's citation had drifted, which is the failure `test_named_controls_resolve_to_real_tests` exists to catch. |
| `test_ablating_a_compartment_that_does_not_exist_is_refused` | §5 — `without=("nucelus",)` raises `WholeCellError`, rather than ablating nothing and reporting no change. |
| `test_a_census_that_lost_a_compartment_would_be_caught` | §4 — the invariant is non-vacuous: a census with a compartment dropped from both `present` and `absent` fails I1, and one that appears in both fails it too. |
| `test_the_absent_reason_for_the_refused_compartment_names_the_gate_that_fires` | §4b — the reason text names the witness/coverage gate **and** the argument that reruns it, so a reader is pointed at the mechanism rather than at an opinion. |
| `test_a_world_with_no_coupling_is_reported_unsteppable_and_the_balance_gate_is_why` | §4e — ablating either coupling endpoint leaves a world that accepts no step; the balance gate is driven in both directions, with the baseline asserted steppable so the check cannot pass for any world. |
| `test_ablating_every_compartment_is_refused_by_the_world_and_not_by_a_wrapper` | §5 — the empty-world refusal is `build_cell_world`'s and still reaches the caller; the census machinery does not swallow it on the way out. |
| `test_forcing_a_compartment_that_is_not_refused_is_refused` | §5 — `schedule_refused` is for driving a recorded refusal, not a general override. |

In the checkable form the port-discipline parser reads:

| Negative (must fail) | `tests/scenarios/test_whole_cell.py::test_nothing_is_refused_any_more_and_the_machinery_that_drove_a_refusal_still_works` | Scheduling `cytosol` makes the world raise `CoverageViolation` at propose time, in all three of its declared phases. |
| Negative (must fail) | `tests/scenarios/test_whole_cell.py::test_ablating_a_compartment_that_does_not_exist_is_refused` | A misspelled `without=` name is refused rather than silently ablating nothing. |
| Negative (must fail) | `tests/scenarios/test_whole_cell.py::test_a_census_that_lost_a_compartment_would_be_caught` | The totality invariant is non-vacuous against a planted omission. |
| Negative (must fail) | `tests/scenarios/test_whole_cell.py::test_a_world_with_no_coupling_is_reported_unsteppable_and_the_balance_gate_is_why` | A world with no coupling accepts no step, the balance gate is why, and the baseline is steppable so the check is not vacuous. |

## 10. Numerical and precision envelope

Everything below is **float64 on `NumpyBackend`**, measured on this Mac before the module was
written, from `build_vertical(subdivision_level=1)` and the fixtures §6 describes.

### 10.1 The composed world at rest — MEASURED

| Quantity | Value |
|---|---|
| compartments | **12 of 13 present**, `cytosol` absent (§4b) |
| owners in the world | 13 (12 compartments + the osmotic envelope, which is a field and not one of the thirteen) |
| connectors | 2 — `cortex_membrane_contact`, `erm_tether` |
| **uncoupled compartments** | **10 of 12** — see §14.2 before reading anything else here |
| build wall clock | 0.46 s |
| total potential energy | `9441.890678061878` pN·µm |
| max residual force | `443.60915712737756` pN |
| constituent force scale | `1262.5682555743924` pN (a sum of magnitudes, never a resultant — `PLAN.md` §2.5) |
| largest per-owner residual | `microtubule`, 4.4361e+02 pN — the default rest spacing does not match the fixture geometry, so its axial term is live at construction |
| owners reporting `internal_forces` | 9 of 13 (`nmii`, `focal_adhesion`, `extracellular_medium` and the envelope carry none) |
| owners contributing a potential energy | 10 of 13 |

### 10.2 A bounded descent — MEASURED

40 steps at the ratified float64 slack (`2 · u64 = 2.220446e-16` relative):

| Quantity | Value |
|---|---|
| accepted / rejected | **40 / 0** |
| energy | `9.441891e+03` → `9.365882e+03` pN·µm |
| max residual at the end | `1.548419e+02` pN |
| final step | `4.526e-03` |
| coverage verdict | **complete** |
| clock | `time = 3.4933985203468423e-04`, `step_index = 40` |
| wall clock | 0.82 s |

**`converged` is `False` at 40 steps and that is expected.** This is a bounded descent that
establishes the world *steps*, not that it reaches equilibrium. No equilibrium number is read off
it, and none is quoted anywhere in this entry.

### 10.3 What the descent cannot see

Three owners contribute **no** potential energy and are still scheduled in `SOLVE` or
`ACCUMULATE_FIELD`: `nmii` (motor heads), `extracellular_medium` (exterior Stokes resistance) and
`focal_adhesion` (no material points of its own). The monotone-potential predicate therefore cannot
observe what they do, and a step that increased their internal work would be accepted on the
strength of the other ten. That is a property of composing a dissipative element into an
energy-descent driver, it is stated rather than fixed, and it is why §14.4 exists.

### 10.4 The two ablations that cannot step — MEASURED

Each built at the same fixture, then driven with the same predicate for up to 12 steps:

| World | couplings | energy | accepted / rejected |
|---|---|---|---|
| baseline | 2 | `9.441891e+03` pN·µm | ≥ 1 accepted |
| `without=("cortex",)` | **0** | `9.441891e+03` pN·µm | **0 / 1**, then the controller breaks out |
| `without=("membrane",)` | **0** | `6.654395e+02` pN·µm | **0 / 1**, then the controller breaks out |

The rejected condition is `balance.connectors_recorded` in both, and `coverage.full_coverage` and
both descent conditions pass in both — including `potential_did_not_increase`, so the candidate
these worlds cannot accept is a candidate that was descending correctly. §4e.

## 11. Production-backend residency and transfer

**Host only.** Every owner is built and stepped on `NumpyBackend`; nothing in this entry touches a
device, allocates a warp array, or reads `~/.aleph_data/gpu_authorization.json`. No GPU run is
performed or requested by this lane. The composed **device-resident** world is a different object
with a different builder (`aleph/runtime/resident_world.py`, `ALEPH-PORT-3614`/`-3616`) and covers
two owners, not twelve; nothing here is evidence about it, and nothing here is a proposal to extend
it.

## 12. Comments and docstrings to discard

Nothing is copied from the reference project, so there is nothing to discard.

One habit of the reference project is deliberately **not** carried: its whole-cell solve declared a
compartment set and produced no work for it. The shape of this module — a census that must be total
and disjoint, plus a coverage verdict read after a real step — is the inverse of that, and §8's two
stepping controls are what make it inverse in fact rather than in intent.

## 13. Acceptance

| | |
|---|---|
| `tests/scenarios/test_whole_cell.py` | **26 passed**, 3.6 s |
| `tests/scenarios/test_cell_world.py` | **16 passed**, 122 s — the composition machinery this entry builds on, re-run unchanged |
| `tests/firewall` + `tests/vertical/test_package_exports.py` | **34 passed** |
| `tests/ports` | **20 passed, 1 failed** before regenerating `INDEX.md`, and the one failure is `test_index_is_not_stale` caused by *this* entry. Verified by moving this file out of the ledger directory and re-running: **21 passed**. Regenerated per §14.6. |
| `ruff check` on both new files | clean |
| `tests/indent` / `tests/spread` inside `tests/scenarios` | **not waited on.** Both are the scenario lane's long-running uncommitted work; neither imports anything this entry touches. `CLAUDE.md` §1: a foreign breakage is reported, not fixed, and a foreign *slowness* is not a reason to hold a commit. |
| `aleph/**` outside `aleph/scenarios/whole_cell.py` | **byte-for-byte untouched** — no compartment, connector, runtime or guard file is edited by this lane, and `aleph/scenarios/__init__.py` is deliberately not extended (it is the scenario lane's file, and `whole_cell` imports fine without it) |
| reviewer | not yet reviewed; `PROPOSED` |
| rollback | delete the two new files and this entry, and regenerate `INDEX.md`. Nothing else in the tree depends on them: no existing module imports `whole_cell`, and the module writes no state, no artefact and no record. |

## 14. Honest limits

### 14.1 Twelve compartments were composed. A cell was not simulated.

Twelve owners were constructed, registered, scheduled, stepped forty times together, and every one
of them witnessed work. That is a statement about **assembly and dispatch**. It is not a statement
that the result resembles a cell, that the parameters are right, that the geometry is anatomical, or
that any number here would survive contact with a measurement. Evidence class `UNVERIFIED`
throughout; quantitative status `BLOCKED`.

### 14.2 The world has twelve compartments and **two** couplings

This is the most misreadable fact in the entry, so it is stated flatly: the only connectors in the
composed world are the vertical's own `cortex_membrane_contact` and `erm_tether`. **Ten of the
twelve present compartments exchange force with nothing.** They occupy the same world, are dispatched
in the same phases, and are integrated in the same transaction, but no connector delivers force
between them.

Two consequences follow and neither is small:

1. **An ablation of an uncoupled compartment changes the other owners' trajectories not at all.** It
   changes the owner count, the total energy and the residual — nothing else. `WholeCell.uncoupled`
   carries exactly which compartments those are, so an ablation result cannot be read without
   meeting the fact.
2. **This is a composition baseline, not a coupled-cell baseline.** The couplings exist —
   `aleph/vertical/wiring.py` carries the declared set and `aleph/vertical/connectors*.py` the
   implementations — but wiring twelve compartments' worth of them requires resolvers, endpoint
   pairings and force-routing adapters per pair, several of which are open work (the spreading lane
   needed a routing adapter for one). That is a separate entry and this one does not pre-empt it.
3. **The world has a single point of failure for steppability**, and §4e measured it: both
   couplings name `membrane` and `cortex`, so ablating either leaves a world that composes and
   accepts no step. That is not a defect in any gate — it is what "two couplings" costs, stated as
   a number rather than as a worry, and it is the strongest argument in this entry for doing the
   coupling work next.

### 14.3 One of the thirteen is absent, and the count is 12 rather than 13 on purpose

`cytosol` is out for the measured reason in §4b. Twelve honest owners beat thirteen with a stub in
one of them, because the stub would satisfy every count in the project while contributing nothing,
and the census would then read *complete*. That is the reference project's failure with better
provenance and it is precisely what a baseline must not do.

### 14.4 The descent is blind to three of the twelve

§10.3. The acceptance predicate is `MonotonePotentialDescent`, which reads a total potential; three
owners contribute none. This is a limitation of driving a world that contains dissipative and
kinetic elements with an energy-descent scaffold, and `aleph/vertical/relax.py`'s own module
docstring says at length that it is scaffolding and must not become the production solver. Nothing
here changes that and nothing here should be read as proposing an acceptance rule — `ALEPH-DQ-104`
is reserved for the PI.

### 14.5 The eight dropped phases are a claim about this fixture

§4c. `DISCOVER` is dropped for six owners because it witnesses nothing **in a static baseline where
no topology changes**. A scenario that inserts a filament, branches a lamellipodium, or breaks a
crosslink needs those passes, and this module's schedule would then be wrong rather than merely
economical. The measurement is carried in the report so the next lane meets it instead of
rediscovering it.

### 14.6 `ports/ledger/INDEX.md` is regenerated by this lane, and the reason changed under it

The first draft of this section said the index was already stale from another lane and would not be
touched. **That was measured and found false**: `tests/ports` is `21 passed` with this entry moved
out of the ledger directory and `20 passed, 1 failed` with it in. The staleness is *this* entry's,
so leaving it would be leaving a red test somebody else has to explain.

It is therefore regenerated, under the rule
`tests/ports/test_index_names_only_tracked_files.py` encodes — **a generated artefact that is
committed must be derivable from what is committed** — which means the order matters: this entry is
staged *first*, the ledger directory is checked for any other lane's untracked entry *second*
(there was none at generation time), and only then is the index generated and committed in the same
commit. Regenerating from a working tree holding another lane's in-progress entry is the defect that
file records happening four times in one night, and the sequence above is what avoids repeating it.

### 14.7 One measurement in this entry was made twice, and the first version was wrong

An earlier survey by this same session reported that five compartments had "no owner class" and that
`sf_arc` had "no module". Both were **false**: the survey searched for class names ending in
`Owner`/`Field`/`Network` and therefore missed `NmiiPopulation`, `FilopodiumBundle`,
`FocalAdhesionClutchGraph`, `ExteriorStokesMedium` and `SFArcGraph`. Every compartment in §4c and
§10 was subsequently **constructed in a live process and driven through its participant protocol**
before being written down here. The lesson is kept in the entry rather than in a session's memory:
enumerate the classes, do not pattern-match the names.
