# ALEPH-PORT-2102 — head-resolved NMII motor connector contracts (group F)

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2102` |
| Lane | `L21 connector contracts — protrusion, NMII motors, absences` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

```python
from aleph.state.connectors_protrusion_motor import (
    NMII_MOTOR_CONNECTORS,          # the four group-F contracts
    NMII_MOTOR_TARGETS,             # ((connector name, actin component), ...) x4
    ConnectorContract,              # resolved from aleph.state.schema (lane L15)
    ConnectorFamily,
    build_connector_contract,
    by_name,
    declared_but_unimplemented,
    validate_contract,
)
```

Aleph target file: `aleph/state/connectors_protrusion_motor.py`.
The six protrusion contracts in the same module are authorised by `ALEPH-PORT-2101`, not here.

The four contracts: `nmii_sf_motor` (nmii ↔ sf_arc), `nmii_cortex_motor` (nmii ↔ cortex),
`nmii_lamellipodium_motor` (nmii ↔ lamellipodium), `nmii_filopodium_motor` (nmii ↔ filopodium).
All `K`, all `commit_on_accept=True`, all `implemented=False`, all `composite_group=None`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/nmii_actuator.py`, `ffn_sim/ac/engine/cortex_motor_slice.py`, `ffn_sim/ac/engine/sf_motor_slice.py`, `ffn_sim/ac/engine/dispatch.py`, `ffn_sim/ac/engine/contracts.py`, `ffn_sim/ac/engine/cortex_state.py` |
| Source symbol(s) | `SegmentMotorConnectorRuntime`, `FilamentMotorConnector`, `CortexMotorConnector` (alias), `FilamentMotorPortView`, `NMIIActuatorView`, `_REQUIRED_CONNECTORS`, the four `NMII_*_MOTOR` name constants, `build_cortex_motor_port`, `canonical_facade_claims()` |
| Read from | **working tree** — `sed`/`grep` on the checked-out files, not `git show` |
| Working tree == commit? | **Five of six yes, one no.** `git diff be0e5876 -- <path>`: `nmii_actuator.py` 0 lines, `cortex_motor_slice.py` 0, `sf_motor_slice.py` 0, `dispatch.py` 0, `contracts.py` 0, **`cortex_state.py` 52 diff lines — read from the working tree, which differs from the commit.** |
| Source test status | Live tests exist, and here they do exercise law as well as plumbing — for the two implemented edges. `test_nmii_actuator.py` 24 tests / 901 lines, `test_cortex_motor_slice.py` 13 / 515, `test_sf_motor_slice.py` 9 / 409. **But the two seamed edges are tested only through fabricated ports:** `test_nmii_actuator.py:889-890` constructs the real connector class under both protrusion names against ports that production never builds, so those two names are covered by tests and uncovered by any production path. That is the sharpest form of the declaration-versus-reality gap in this group — a green test for an edge that cannot run. |

The `cortex_state.py` divergence was inspected rather than assumed harmless. The uncommitted change
adds an `INCUMBENT_CHANNEL_NAME` map and a derived channel-omission helper for a double-counting
guard between two force accumulators. It does not touch `build_cortex_motor_port`, which is the only
symbol this entry cites from that file. So the audit conclusion below is unaffected — but the file
was read from the working tree and this row says so, because a ledger entry that cites a commit and
was written from the working tree is citing a revision it never read.

**Per-file digest of the exact bytes read** (`shasum -a 256`, working tree). Recorded so a later
reader can prove what was in front of the author, which a commit id cannot do when the tree is
dirty — and for `cortex_state.py` the two digests differ, which is the whole reason the row exists:

| Source file | Digest |
|---|---|
| `ffn_sim/ac/engine/nmii_actuator.py` | `sha256:1523ed040a259971a55cc34281b6b6c763f7b9a299fa0901e03ebf5d79a2ceae` |
| `ffn_sim/ac/engine/cortex_motor_slice.py` | `sha256:1ab66a85256d4f4db5a39fb866c80a00b52f3a120b2e5b902284280d8c0cc8c4` |
| `ffn_sim/ac/engine/sf_motor_slice.py` | `sha256:0812858af4006d5c92f087e166e9c3ce273d09614f1a1a2661fc489d45fd55d7` |
| `ffn_sim/ac/engine/dispatch.py` | `sha256:325307d38f7df0989c1bf71f421f5fb8ea9cf5341f3dd739e2ad890f70d2311a` |
| `ffn_sim/ac/engine/contracts.py` | `sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72` |
| `ffn_sim/ac/engine/cortex_state.py` — **read (working tree)** | `sha256:1e5075143eac031df29680e5c396ebfcdb023d591253a03af0c00ec20782cef3` |
| `ffn_sim/ac/engine/cortex_state.py` — **at `be0e5876` (not read)** | `sha256:d28f254496f7398169ec0f14cd78970a60f7b0a830772a1ef39c2dbad938eab3` |

## 3. Why source-derived porting beats clean-room

**It does not, and no code was ported.** Two things crossed the boundary, and neither is code:

1. **An audit result** — which of the four motor connectors has a runtime behind it. That is a fact
   about the reference, measurable only by reading its call graph.
2. **One design decision, adopted with its reason re-derived**: NMII acts on each actin component
   through a *separate* connector. That is not a convention a competent author would necessarily
   arrive at — the obvious implementation is one active-stress or active-tension field over all
   actin, and it is cheaper. §4 states why the separate-connector form is the one that can carry
   head-resolved kinetics, and that argument was written here rather than quoted.

The mechanism text and the mechanical interpretations are new prose. What was **not** taken: any
kernel, any crossbridge or attach/detach law, any parameter, and — deliberately — the reference's
`MOTOR` family name. Aleph's contract has two families and the manuscript legend marks all four
entries `[K]`, so `MOTOR` would have been an inherited third family Aleph has not decided to have.

## 4. Physical or mathematical law represented

**A structural commitment about representation, plus the crossbridge cycle it presupposes. No
constitutive law is written down here.**

The commitment: an individual NMII head binds a **live polar actin material coordinate**, executes a
working stroke along that filament's polarity, responds to load through a force-velocity relation,
and detaches at a load-dependent rate. Four connectors, one per actin population.

Why four and not one, derived rather than asserted:

- **An aggregate field has no head.** Binding and detachment are events in a discrete population;
  a continuous tension field has no population to be discrete over, so per-head kinetics cannot be
  expressed in it at all. The load dependence of the duty ratio — bound fraction rising or falling
  under sustained tension — then has to be supplied by hand as a constitutive choice, which means
  the model's most interesting behaviour becomes an input rather than an outcome.
- **An aggregate field loses the accounting.** With one field, "how much active work went into
  cortical actin rather than into stress fibres" has no answer, because the four populations were
  summed before the question was asked. Four connectors keep the ledger separable by construction.
- **The populations are disjoint, so the coordinate must name its owner.** The cortex, sf_arc, the
  lamellipodium and the filopodium own disjoint filament populations. A head bound to "actin" is
  bound to nothing in particular; a head bound to a *material coordinate on a named owner's
  filament* is a well-posed attachment. This is the same ownership discipline as `ALEPH-PORT-2101`
  §4, arriving from the motor side.
- **"Live" is load-bearing.** The coordinate moves and its polarity is re-evaluated as the filament
  deforms. Binding to a fixed lab-frame position would give a stroke with no defined direction and a
  bond that resists filament motion rather than walking along it.

The transfer rule follows from the endpoints being two owners: `+f` on the head, and the equal
opposite reaction distributed onto the actin segment the head is bound to, **accumulated into the
two owners' force arrays independently and never merged**. Two arrays rather than one is what makes
the force-closure check non-vacuous: if both sides were written into one array, closure would hold
by construction and would test nothing.

`composite_group` is `None` for all four, and that is a positive decision. A composite group means
"dispatch these as one series joint". Grouping the four motors would rebuild the aggregate this
section argues against, in the dispatch layer, where it would be harder to see.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| — | — | — | This module declares no numeric quantity. |

No stall force, no zero-load velocity, no attachment rate, no working-stroke length is declared
here. Those are a runtime's parameters and each will need its own sourced constant with provenance
under Aleph's units discipline. Declaring any of them in a contract registry would smuggle in a
parameter with no source.

Singular and boundary cases — as `ALEPH-PORT-2101` §5, plus:

- The four contracts must address four *distinct* `endpoint_b` values. A repeat would mean two
  connectors compete for one actin population with no rule for which one loads it.
- `endpoint_a` is `nmii` for all four. A motor connector whose head side is not the motor owner has
  its adjoint the wrong way round.

Invariants that must hold, each with the test that asserts it:

- **I1.** Four contracts, endpoints `nmii ↔ {sf_arc, cortex, lamellipodium, filopodium}`, four
  distinct targets —
  `tests/state/test_connectors_protrusion_motor.py::test_nmii_acts_through_four_separate_connectors_and_no_aggregate_group`.
- **I2.** No motor contract declares a composite group — same test. This is the aggregate-field
  refusal, pinned.
- **I3.** Endpoint roles are exactly `individual NMII head crossbridge` and `live polar actin
  material coordinate` — `::test_motor_endpoint_roles_are_head_and_live_polar_actin_coordinate`.
- **I4.** Each mechanism records all four elements of the cycle — binding, working stroke,
  force-velocity, load-dependent detachment —
  `::test_motor_mechanism_records_the_full_head_resolved_cycle`. Parametrised over the four, since
  any three of the four can be implemented without the fourth and still look like a motor.
- **I5.** All four are `K` with `commit_on_accept=True` —
  `::test_every_connector_is_kinetic_and_commits_only_on_accept`.
- **I6.** All four bidirectional with the adjoint required —
  `::test_every_connector_is_bidirectional_and_requires_the_adjoint`.
- **I7.** None claims to be implemented in Aleph —
  `::test_no_connector_claims_to_be_implemented_in_aleph`.

## 6. Source evidence class and known retractions

**Audited by reading the call graph. The earlier finding was right in its verdict and wrong in its
reason for two of the four, and the correction matters.**

| Connector | Verdict | Evidence |
|---|---|---|
| `nmii_sf_motor` | **IMPLEMENTED** | `sf_motor_slice.py:492` constructs `FilamentMotorConnector` and `:203` builds a real `FilamentMotorPortView` over `sf_arc`-owned arrays. The runtime scatters a two-array adjoint crossbridge force and runs an accepted-predicated attach/detach cycle. |
| `nmii_cortex_motor` | **IMPLEMENTED** | `cortex_motor_slice.py:414` `FilamentMotorConnector`, whose `name` defaults to this connector; the cortex port is built by `cortex_state.py:454` `build_cortex_motor_port`, and `composed_native.py` binds the pair. |
| `nmii_lamellipodium_motor` | **SEAM** | Name constant `nmii_actuator.py:85`; required-edge entry `:523`; contract `contracts.py:724`. **No production call site constructs a connector for it and no lamellipodium port exists.** |
| `nmii_filopodium_motor` | **SEAM** | Name constant `nmii_actuator.py:86`; required-edge entry `:524`; contract `contracts.py:732`. Same absence. |

**The correction.** The earlier audit recorded `nmii_lamellipodium_motor` and
`nmii_filopodium_motor` as "Protocol seams with no concrete runtime". The verdict SEAM is right; the
stated reason is not, and the difference changes what would have to be built:

- There **is** a concrete runtime class that accepts them. `SegmentMotorConnectorRuntime`
  (`nmii_actuator.py:977`) validates `name` against `_REQUIRED_CONNECTORS`, which contains all four
  names, and `FilamentMotorConnector` was made target-agnostic — every kernel it launches addresses
  the target only through the `FilamentMotorPortView` interface, and `name` / `component_b` are
  constructor parameters. Its own docstring anticipates the protrusion edges.
- The reference's own tests construct that real runtime under both names
  (`tests/ac/engine/test_nmii_actuator.py:889-890`), against fabricated ports.
- So what is missing is **not** a connector implementation. It is a `FilamentMotorPortView` for the
  lamellipodium and the filopodium. Only two ports are built in production — cortex
  (`cortex_state.py`, plus one script) and sf_arc (`sf_motor_slice.py`). The protrusion state
  owners expose a `ProtrusionActorView` with node-level position/force arrays and **no segment
  table** (`segment_node_a_d` / `segment_node_b_d` / `segment_polarity_d`), which is exactly what
  the motor kernels address the target through.

That is a more useful and a more demanding finding than "no runtime". A branched lamellipodial graph
and a fascin-bundled filopodial graph have no segment/polarity table of the kind a segment-anchored
motor binds to, so the missing piece is a representation decision about the protrusion actors, not a
copy of an existing connector. Recorded here because a reader who believed "no concrete runtime"
would plan the wrong work.

**None of this makes any of the four implemented in Aleph.** Aleph has no NMII owner, no actin
segment table, and no motor kernel. All four contracts carry `implemented=False` and the flag
answers one question — does Aleph evaluate it.

**The coverage-gate defect applies to two of the four.** `dispatch.py:180-190` claims all four for
`NMIIActuator.accumulate_candidate` and `validate_canonical_facade_claims` checks only single
non-duplicate ownership. So `nmii_lamellipodium_motor` and `nmii_filopodium_motor` pass the coverage
gate with nothing behind them — the `membrane_cortex_contact` pattern again.

**Retractions: found, and they are substantive.** Unlike the protrusion group, this one has live
claims in the reference's own status documents, and two are withdrawn:

- `STATE.md:35` — "the `nmii_sf_motor` force field is non-conservative", quotable as attribution
  and scaling exponents only, **no magnitude**, at 904 nodes / 320 heads.
- `STATE.md:38` — the acceptance predicate tests **adjoint closure, not convergence**, so a
  non-converged step is not what it rejects. Recorded because "the predicate said yes" and "the
  step converged" are different facts there.
- `STATE.md:87` / `STATE_NONQUOTABLE.md:64` — the "≈20× cost multiplier" for `nmii_sf_motor` is
  **withdrawn**, confounded by solver budget rather than residual. Their own re-reading points the
  other way and is the better result.
- Both `nmii_sf_motor` claims are at **0.18% of native population**. Any cost or conditioning
  number taken from them would be an extrapolation across nearly three decades of population.

Nothing was taken from any of these, so none of the retractions is inherited. They are recorded
because an entry that says "no retraction found" without looking is worthless, and here there was
something to find.

## 7. Independent oracle or derivation

**No numerical oracle. There is no number in this module.** As `ALEPH-PORT-2101` §7: a registry of
declarations cannot be checked against a closed form, and claiming otherwise would be the promotion
PLAN §0.5 forbids.

The independent checks used are the manuscript appendix — endpoints, both endpoint roles, and all
four elements of the mechanism asserted field-by-field against the extracted text, which is
independent of the reference source — and shape agreement with lane L15's contract type verified
through `dataclasses.fields()` rather than assumed. Agreement with `ffn_cellsim` is not used as
evidence anywhere in this entry.

The oracle a *runtime* would need is worth naming now so nobody reaches for a weaker one when it
lands: a single head against a rigid filament at zero load must walk at the declared unloaded
velocity; at the stall force it must not walk; the work delivered by one full working stroke against
a known spring must equal the spring's stored energy to round-off; and the adjoint must close
exactly — `Σ` head forces `+ Σ` actin reactions `= 0` in the two never-merged arrays. Those are
Aleph-owned checks and none of them requires the reference.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_connectors_protrusion_motor.py::test_nmii_acts_through_four_separate_connectors_and_no_aggregate_group` | Four contracts; `endpoint_b` is exactly `("sf_arc", "cortex", "lamellipodium", "filopodium")`; four distinct targets; `endpoint_a == "nmii"` for all four; `composite_group is None` for all four. |
| Positive | `tests/state/test_connectors_protrusion_motor.py::test_motor_mechanism_records_the_full_head_resolved_cycle` | Each of the four records binding, working stroke, force-velocity response, and load-dependent detachment. Parametrised, so a mechanism that loses one element fails on the connector it affects. |
| Positive | `tests/state/test_connectors_protrusion_motor.py::test_motor_endpoint_roles_are_head_and_live_polar_actin_coordinate` | Both roles exactly: `individual NMII head crossbridge`, `live polar actin material coordinate`. Rules out a bundle-level force and a lab-frame binding site. |
| Positive | `tests/state/test_connectors_protrusion_motor.py::test_every_connector_is_kinetic_and_commits_only_on_accept` | All four `K`, all four `commit_on_accept=True`. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_connectors_protrusion_motor.py::test_construction_rejects_bidirectional_false` | A motor contract with `bidirectional=False` cannot be constructed. For a motor this is the sharpest case: a head that pushes actin without the reaction on the head is a self-propelling motor, and it produces smooth plausible motion. |
| Negative (must fail) | `tests/state/test_connectors_protrusion_motor.py::test_validation_rejects_a_kinetic_connector_that_commits_early` | A `K` motor with `commit_on_accept=False` is refused — a rejected step would leave heads bound that the accepted world never bound. |
| Negative (must fail) | `tests/state/test_connectors_protrusion_motor.py::test_validation_rejects_a_connector_that_does_not_require_the_adjoint` | `adjoint_required=False` is refused. |
| Negative (must fail) | `tests/state/test_connectors_protrusion_motor.py::test_no_connector_in_this_lane_declares_a_composite_group` | Fails if any of the ten acquires a composite group. For the four motors this is the aggregate-field refusal: grouping them for dispatch is exactly how the separate-connector decision would be undone.\* |
| Negative (must fail) | `tests/state/test_connectors_protrusion_motor.py::test_no_connector_claims_to_be_implemented_in_aleph` | Fails the moment any contract claims `implemented=True` without an evaluator behind it — the declaration-versus-reality gap, caught in Aleph's own tree. |

\* This one is a guard rather than a rejection test: it asserts a current property and fails on
change. It is listed as a negative control because what it protects against is a future edit, and
because `test_construction_rejects_bidirectional_false` above already covers the
gate-watched-failing requirement with an exception that has been observed.

Positive complement, so the gates are not vacuous:
`test_a_continuous_connector_without_kinetic_commit_is_accepted`.

## 10. Numerical and precision envelope

None. No arithmetic, no accumulation, no tolerance, no conditioning. Every assertion is exact:
string identity, tuple identity, `is True` / `is False`, integer counts, exception types.

Recorded deliberately rather than left blank, because the *runtime* this contract describes will
have a demanding envelope and it should not inherit an empty section from here. A head-resolved
motor accumulates a force from every bound head into two owners' arrays; with `float32` compute the
accumulation order matters at native head counts, and `ALEPH-DQ-107` already says accounting-
critical accumulation is `float64`. That decision belongs to the runtime's own ledger entry. This
entry establishes nothing about it.

## 11. Production-backend residency and transfer

**This code will never run on the production backend.** Host-side declarative registry: four frozen
dataclasses built once at import, resident in CPU memory. No device array, no kernel, no per-step
work, nothing transferred in either direction. `float32`/`float64` do not arise. Pure standard
library — no `numpy`, no `scipy`, no backend import — which is also what keeps it importable under
gate R1 with the reference tree absent.

The residency question for a real head-resolved motor is a genuinely hard one and this entry
answers none of it: split ownership means the head node lives with the motor owner and the segment
endpoints with the actin owner, so a crossbridge kernel reads two owners' arrays and writes two
owners' arrays without merging them. Whoever writes that runtime inherits that problem, not a
solution to it.

## 12. Comments and docstrings to discard

Nothing was carried. Read and deliberately left behind:

- Provider module paths and package namespace in every form.
- Its runtime-status vocabulary — `SEAMED`, `KERNEL_BOUND`, `magnitudes GAP`, `TARGET-AGNOSTIC
  (2026-07-25)`, `DECLARED ONLY` — and its gate labels. Its self-assessment scale is not Aleph's.
- Named-author attributions for specific kinetic laws, and the laws themselves. A physiological
  catch-slip detachment law and a Hill-type force-velocity law are cited by name in the source;
  neither the citation nor the law crosses into Aleph. When Aleph needs a detachment law it will be
  read from the literature and sourced under Aleph's units provenance discipline, not inherited from
  a docstring.
- Dated "this landed on <date>" assertions and PI-ratification notes naming decisions Aleph's PI
  has not taken.
- Internal document references (integration notes, plan files, lane and item numbers).
- Its `MOTOR` connector family — see §3, discarded as an inherited classification decision.
- Its device-residency preconditions, per §11.

**What replaces them:** new prose whose organising principle is that each structural commitment is
stated together with what is lost if it is dropped. The four motor contracts share one mechanism
string and one interpretation suffix, composed once at module level, so the four cannot drift into
four subtly different descriptions of the same cycle — which is the failure mode of writing the
same paragraph four times.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Control file green on 2026-07-30: `tests/state/test_connectors_protrusion_motor.py`, all motor invariants I1–I7 asserted, 0 failed. Status remains `PROPOSED`. |
| Reviewer | **Agent-proposed, unratified.** Written by the L21 lane agent during autonomous operation; no human has reviewed it. The §6 audit — including the correction to the earlier "Protocol seams with no concrete runtime" verdict — is evidence about the reference tree, not about this code. |
| Rollback | Remove `NMII_MOTOR_CONNECTORS`, `NMII_MOTOR_TARGETS` and the two shared mechanism strings from `aleph/state/connectors_protrusion_motor.py`, and the motor tests from its control file. Nothing imports them today. What is lost is the four-connectors-not-one-field decision and the reason for it; the aggregate form is the natural default, so its absence is silent. |

Why `PROPOSED`: there is no oracle to qualify against (§7), and a registry passing tests about its
own declarations is a consistency result rather than a correctness one.

## 14. Honest limits

What this entry does **not** establish:

- **No motor physics whatsoever.** No crossbridge, no working stroke, no force-velocity relation,
  no attachment or detachment rate, no stall force. Four contracts declared is not a motor. Aleph
  evaluates none of them and none can generate a newton of force.
- **The four elements of the cycle are recorded as requirements, not derived.** `test_motor_mechanism_records_the_full_head_resolved_cycle`
  checks that the words are present. It cannot check that any of them is correct, and no test in
  this lane could.
- **`INHERITED_UNVERIFIED`: the argument for four connectors is a design argument, not a
  measurement.** §4 argues that an aggregate field cannot carry per-head kinetics or population-
  resolved accounting. That is sound as reasoning and is not evidence: nobody has run both forms and
  measured a difference, here or in the reference. If a measurement ever shows the aggregate
  reproduces the head-resolved answer within the observable's error, this decision costs
  performance for nothing and should be revisited.
- **The §6 audit rests on a static read, not on running the reference.** The 05:05 ports lane
  audited by *running*, which is stronger. A dynamically constructed connector or a port built in a
  script outside the searched paths would not appear in a grep. What is solid regardless: no
  production `FilamentMotorPortView` exists for the lamellipodium or the filopodium, and the
  protrusion actor views expose no segment/polarity table for one to be built from.
- **One cited file was read at a revision that is not the commit.** `cortex_state.py` differs from
  `be0e5876` by 52 diff lines in the working tree. The diff was inspected and does not touch the
  cited symbol, but the read was from the working tree and §2 says so.
- **`nmii_cytosol_transfer` is the fifth NMII edge and is not in this lane.** The reference claims
  it on the same facade, and it is what stops an all-heads-unbound minifilament being a free rigid
  body. It belongs to the cytosol-transfer group (lane L20). Named here so its absence from this
  entry is not read as an absence from the registry.
- **No cross-lane count.** Whether the four motor contracts plus the other 32 sum to 36 with no
  duplicate name and no missing edge is not checked by any test in the tree, because six lanes are
  landing concurrently. That control has to be written above the lanes.
- **No GPU, and none possible.** §11: nothing here has a device side. Zero GPU jobs were run for
  this entry.
