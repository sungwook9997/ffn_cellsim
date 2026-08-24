# ALEPH-PORT-2101 — protrusion connector contracts (group E) and the absence register

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2101` |
| Lane | `L21 connector contracts — protrusion, NMII motors, absences` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

```python
from aleph.state.connectors_protrusion_motor import (
    PROTRUSION_CONNECTORS,          # the six group-E contracts
    ConnectorContract,              # resolved from aleph.state.schema (lane L15)
    ConnectorContractError,
    ConnectorFamily,
    SCHEMA_SOURCE,
    build_connector_contract,
    by_name,
    connector_names,
    declared_but_unimplemented,
    validate_contract,
)
```

Aleph target file: `aleph/state/connectors_protrusion_motor.py`.
The four NMII motor contracts in the same module are authorised by `ALEPH-PORT-2102`, not here.

### 1b. The absence register — a companion module, and a different provenance

```python
from aleph.state.absent_connectors import (
    ABSENT_CONNECTORS,              # 19 typed records, each with reason + re-entry condition
    AbsentConnector, AbsenceRecordError,
    CHEMICAL_FLUX_LIMIT, chemical_flux_connector_count,
    ScopeDomain, ScopeRefusal, UnsupportedClaim,
    absences_for_domain, absent_connector_names, by_name,
    model_scope_statement, require_in_scope,
    unsupported_claim_domains, unsupported_claims,
)
```

Aleph target file: `aleph/state/absent_connectors.py`.

This module is listed under this entry because it landed in the same change and shares its
reviewer, **but its provenance is different and weaker in the direction that matters**: nothing in
it was derived from the reference tree at all. Its content comes from the manuscript appendix
already extracted into Aleph's own tree at
`docs/manuscripts/extracted/APPENDIX_A_mechanics_registry.txt` — the `INTENTIONALLY ABSENT
CONNECTORS`, `CURRENT CHEMICAL-CONNECTOR LIMIT` and `MODEL-SCOPE STATEMENT` sections. No provider
file was read for it, so there is no source symbol to cite and no comment to strip. It gets an
entry anyway, because "nobody consulted the reference" and "nobody wrote down whether they
consulted the reference" are the same absence in the ledger and very different facts.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/protrusion.py`, `ffn_sim/ac/engine/contracts.py`, `ffn_sim/ac/engine/dispatch.py`, `ffn_sim/ac/engine/cortex_population.py` |
| Source symbol(s) | `ProtrusionGraphConnector` (Protocol), `ProtrusionStepBindings`, `ProtrusionActors.accumulate_mechanics`, `_REQUIRED_CONNECTORS`, the eight `LAM_*` / `FILO_*` name constants, `canonical_facade_claims()` |
| Read from | **working tree** — `sed`/`grep` on the checked-out files, not `git show` |
| Working tree == commit? | **yes for every file this entry cites.** Verified per file with `git diff be0e5876 -- <path>`: `protrusion.py` 0 lines, `contracts.py` 0, `dispatch.py` 0, `cortex_population.py` 0. |
| Source test status | Live tests exist and exercise the **plumbing, not the law**: `ffn_sim/tests/ac/engine/test_protrusion.py`, 17 tests / 726 lines, and it is the *only* file in the tree that implements the connector interface or constructs the bindings — with test doubles. It checks binding completeness, endpoint agreement, ownership disjointness and storage identity. No test asserts a protrusion connector force, because no production connector exists to produce one. |

The provider working tree carries 31 uncommitted changes at `be0e5876`. Of the files read for this
lane only `ffn_sim/ac/engine/cortex_state.py` differs (52 diff lines), and it is cited by
`ALEPH-PORT-2102`, not here. Every file cited above is byte-identical to the commit, so citing the
commit is accurate for this entry rather than merely conventional.

**Per-file digest of the exact bytes read** (`shasum -a 256`, working tree; identical to the commit
blob for all four, per the row above). Recorded so a later reader can prove what was in front of the
author, which a commit id alone cannot do when the tree is dirty:

| Source file | Digest |
|---|---|
| `ffn_sim/ac/engine/protrusion.py` | `sha256:98300c76d61b3a6d572d91fa838b27fe127a50bb43eb19040fd8e81baba8fc63` |
| `ffn_sim/ac/engine/contracts.py` | `sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72` |
| `ffn_sim/ac/engine/dispatch.py` | `sha256:325307d38f7df0989c1bf71f421f5fb8ea9cf5341f3dd739e2ad890f70d2311a` |
| `ffn_sim/ac/engine/cortex_population.py` | `sha256:04d4542cf0c07a1f35bf1600eda4b40828e60e0ecebb60d71388ac363e903b70` |

## 3. Why source-derived porting beats clean-room

**It does not, and no code was ported.** What crossed the boundary is one *audit result*: which of
the eight declared protrusion connectors has a runtime behind it. That is a fact about the
reference, not a design, and it cannot be re-derived from first principles — only measured, by
reading the tree and following the call graph.

Everything else was written from the manuscript appendix, which is Aleph's own extracted document.
Endpoint names, roles and mechanisms are re-stated in Aleph's vocabulary; the mechanical
interpretations are new prose that states the failure mode each ownership rule prevents, which the
appendix does not do.

There is a positive reason not to port here: the reference's family vocabulary is *finer* than
Aleph's and disagrees with the manuscript. Its `ConnectorFamily` classifies
`lamellipodium_membrane_contact` as `CONTACT`, `lamellipodium_cortex_seam` as `TRANSIENT_ACTIN`,
`lamellipodium_nascent_fa` as `ACTIN_ANCHOR`, and the cytosol edges as `IMMERSED_TRANSFER`.
Aleph's contract has two families, `K` and `C`, and the appendix's connector legend marks all six
group-E entries `[K]`. Inheriting the finer taxonomy would import a classification decision Aleph
has not made, so the two-family reading was taken from the legend and the finer names discarded.

## 4. Physical or mathematical law represented

**There is no law here. What crosses the boundary is an ownership discipline plus two mechanisms,
and it is worth saying which is which.**

The *discipline* is not physics; it is a constraint on how physics may be assembled:

1. Every physical filament belongs to exactly one component. The cortex, the stress-fibre/arc
   system, the lamellipodium and the filopodium own disjoint filament populations.
2. No physical connection is created by spatial overlap, array co-location, or shared indexing.
   Components exchange force only through a declared connector.
3. Every mechanical connector is bidirectional and transfers through an adjoint reaction.
4. A kinetic connector commits binding, unbinding, capture, remapping or topology change only
   after an accepted physical step.

Rule 3 is the one with a derivation behind it. A transfer that adds `+f` to one endpoint and
nothing to the other is a momentum source: `Σ F ≠ 0` over the closed assembly by exactly the
missing reaction. The force-closure check cannot catch it, because the missing term was never
accumulated — there is no residual to be out of balance. So bidirectionality is refused at
construction, where it is still visible, rather than checked downstream where it is not.

The two *mechanisms* are physics, and both are recorded rather than generalised:

- **Brownian ratchet** (`lamellipodium_membrane_contact`, `filopodium_membrane_tip`). The
  membrane's thermal excursion opens the gap a monomer needs; inserting the monomer rectifies the
  excursion. Polymerisation velocity and carried load are one coupled quantity — a load-dependent
  growth rate — not a kinematic rate plus a separate contact force. A generic non-penetration
  contact has no such coupling: it resists interpenetration and cannot convert polymerisation into
  work against a load, so substituting it leaves a leading edge that pushes with a stiffness
  instead of with a velocity.
- **Dynamic transient attachment** (`lamellipodium_cortex_seam`, `filopodium_cortex_root`). A bond
  population with its own binding and unbinding state, joining two filament populations that stay
  disjoint. The distinction from a weld is not stylistic: a welded seam can hold the cortex rigidly
  and can never release, and a spliced filopodial root inherits the cortex's crosslink topology and
  material coordinates, so a retracting filopodium would have to be implemented as a topology edit
  of the cortex.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| — | — | — | This module declares no numeric quantity. |

A contract registry carries names, roles, families and flags. There is nothing to convert and no
tolerance to state. Every quantity a runtime attaches to one of these contracts belongs to the
owner on each side, in Aleph's pN / µm / s convention.

Singular and boundary cases, each with the behaviour Aleph requires:

- `bidirectional=False` — **refused at construction.** Not corrected, not warned about.
- A `K` contract with `commit_on_accept=False` — refused by `validate_contract`.
- A `C` contract with `commit_on_accept=True` — refused; it claims kinetic state it does not own.
- `adjoint_required=False` — refused.
- An empty name — refused; an anonymous contract cannot be looked up, so it can only be dead weight.
- A duplicate name in the registry — refused at import, not at first use.
- `endpoint_a == endpoint_b` with `internal_to=None` — refused. None of the six is internal.
- An unknown name passed to `by_name` — raises `KeyError`. A `None` return would be checked by some
  callers and not others.
- An absence record with a blank or token reason or re-entry condition — refused
  (`AbsenceRecordError`). Minimum 24 characters, which rejects `"n/a"`, `"none"` and `"scope"`
  without demanding an essay.
- A `ScopeDomain` member with no absence behind it — refused at import. An empty refusal list reads
  exactly like being in scope.
- An unrecognised domain name passed to a scope query — `ValueError`, deliberately *not*
  `ScopeRefusal`. A typo must not be reported as "nothing is missing there".

Invariants that must hold, each with the test that asserts it:

- **I1.** Six group-E contracts, all `K`, all `commit_on_accept=True` —
  `tests/state/test_connectors_protrusion_motor.py::test_every_connector_is_kinetic_and_commits_only_on_accept`.
- **I2.** All bidirectional with the adjoint required —
  `::test_every_connector_is_bidirectional_and_requires_the_adjoint`.
- **I3.** The rear seam is not a weld — `::test_the_rear_seam_is_not_a_weld`.
- **I4.** The filopodial root goes through a connector, not shared nodes —
  `::test_the_filopodial_root_goes_through_a_connector_not_shared_nodes`.
- **I5.** Both membrane-contact connectors record the ratchet —
  `::test_the_membrane_contact_connectors_are_brownian_ratchets`.
- **I6.** Nothing claims to be implemented in Aleph —
  `::test_no_connector_claims_to_be_implemented_in_aleph`.
- **I7.** Nineteen absences, each with a reason and a re-entry condition —
  `tests/state/test_absent_connectors.py::test_every_absence_carries_a_reason_and_a_reentry_condition`.
- **I8.** No registered connector uses `CHEMICAL_FLUX` —
  `::test_no_registered_connector_uses_the_chemical_flux_family`.
- **I9.** The unsupported-claims query is non-empty, for every domain —
  `::test_every_domain_has_at_least_one_unsupported_claim`.
- **I10.** The contract type has the shape lane L15 declares —
  `::test_the_contract_type_has_the_agreed_shape`.

## 6. Source evidence class and known retractions

**Audited by reading the call graph, not by trusting the register.** Verdict per connector, for all
eight edges the reference declares on the protrusion facade — the six this lane owns plus the two
cytosol-transfer edges owned by lane L20, included because they share the same interface and the
same absence of a runtime:

| Connector | Verdict | Evidence |
|---|---|---|
| `lamellipodium_membrane_contact` | **SEAM** | Contract declared `contracts.py:656`; name constant `protrusion.py:44`; interface is the `ProtrusionGraphConnector` Protocol `protrusion.py:126`. No concrete implementation. |
| `lamellipodium_cortex_seam` | **SEAM** | `contracts.py:663`, `protrusion.py:45`. A cortex-side endpoint domain *is* registered (`cortex_population.py:348`, node domain) — so the target has a declared endpoint and no connector to reach it through. |
| `lamellipodium_nascent_fa` | **SEAM** | `contracts.py:676`, `protrusion.py:47`. |
| `filopodium_membrane_tip` | **SEAM** | `contracts.py:683`, `protrusion.py:48`. |
| `filopodium_cortex_root` | **SEAM** | `contracts.py:690`, `protrusion.py:49`; cortex endpoint domain `cortex_population.py:352`. |
| `filopodium_nascent_fa` | **SEAM** | `contracts.py:703`, `protrusion.py:51`. |
| `lamellipodium_cytosol_transfer` (L20) | **SEAM** | `contracts.py:670`, `protrusion.py:46`. |
| `filopodium_cytosol_transfer` (L20) | **SEAM** | `contracts.py:697`, `protrusion.py:50`. |

**Confirmed: all eight are seams.** `grep -rn accumulate_actor` over the whole tree returns three
hits — the Protocol's own method definition, the one call site inside
`ProtrusionActors.accumulate_mechanics`, and a string in an `__all__`-style hook tuple in
`actor.py`. The only file that constructs a `ProtrusionStepBindings` or implements
`accumulate_actor` is `ffn_sim/tests/ac/engine/test_protrusion.py`, with test doubles. Nothing in
production ever builds the bindings, so the loop that would call the connectors never runs.

**Not seams, and worth recording as the honest other half:** the two protrusion *state owners* are
real. `LamellipodiumBranchAngleMechanics` launches an Arp2/3 angle-harmonic kernel and
`FilopodiumFascinBundleMechanics` a fascin crosslink kernel, both through an injectable launch
seam. So the components have mechanics and the *couplings between components* do not — which is
the shape this project exists to notice, since a component that computes internal forces looks
alive from the outside.

**The coverage-gate defect reproduces here, on all eight.** `dispatch.py:163-176` registers a
`CanonicalFacadeClaim` assigning every one of the eight to
`ProtrusionActors.accumulate_mechanics`, and `validate_canonical_facade_claims` checks only that
each connector has exactly one non-duplicate owner. Ownership is asserted; evaluation is not
checked. This is the same defect the 03:00 audit found for `membrane_cortex_contact` — a connector
that passes the dispatch coverage gate and is never evaluated — and it applies to eight more
connectors than that finding named. Aleph's runtime lane already refuses this: coverage there needs
scheduled **and** reached **and** evaluated-something, with a witness count of backend launches
(`tests/runtime` — `test_static_schedule_gate_alone_is_insufficient`).

Retractions searched for and not found: `STATE.md` and `STATE_NONQUOTABLE.md` were grepped for all
eight names. No row mentions any of them — no claim, and therefore no retraction. The one
protrusion-adjacent artifact reference is `virtual_cell/archetypes.py:1000`, which cites
`lamellipodium_nascent_fa` as a `nearest_existing` neighbour for an archetype lookup; that is in a
tree the 03:00 audit classified as unreachable from production, and it makes no mechanical claim.

The whole-cell artifact `outputs/ac/cell_assembled/composed_world_v2.json` self-labels
`evidence: CENSUS-WIRED` with `quantitative_claim_status: BLOCKED`, which is the reference's own
statement that its census is ahead of its physics. Nothing here contradicts that label.

## 7. Independent oracle or derivation

**There is no numerical oracle, because there is no number.** Saying otherwise would be the
promotion this project's honesty rules forbid — a registry cannot be checked against a closed form.

What is available, and is used, is an independent *documentary* check plus a structural one:

1. **The manuscript appendix**, extracted into Aleph's tree before this lane ran, is independent of
   the reference source: the six connectors' endpoints, endpoint roles, mechanisms and ownership
   rules were taken from it and each is asserted field-by-field in
   `tests/state/test_connectors_protrusion_motor.py`. Where the appendix and the source disagree —
   the family taxonomy, §3 — the appendix's connector legend won and the disagreement is recorded
   rather than reconciled silently.
2. **Shape agreement with lane L15**, checked by `dataclasses.fields()` rather than by assumption.
   The contract type is imported defensively; if the shared module's field set does not match, the
   module falls back to a local definition of identical shape and *records that it did*
   (`SCHEMA_SOURCE`), and a test fails naming the drift. At the time of writing L15's
   `aleph/state/schema.py` did not exist and `SCHEMA_SOURCE` read `local-fallback`; it landed
   during the lane and now reads `aleph.state.schema`, with the field order matching
   `REQUIRED_CONTRACT_FIELDS` exactly.

Neither of these is agreement with `ffn_cellsim`, which would not be acceptable evidence.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_connectors_protrusion_motor.py::test_all_ten_connectors_are_registered_with_the_expected_names` | The registry holds exactly ten contracts, six protrusion and four motor, under the ten expected names in registry order. |
| Positive | `tests/state/test_connectors_protrusion_motor.py::test_endpoints_are_the_declared_component_pair` | Each of the ten joins the declared component pair — parametrised, so a swapped endpoint fails on the one connector it affects. |
| Positive | `tests/state/test_connectors_protrusion_motor.py::test_the_rear_seam_is_not_a_weld` | `lamellipodium_cortex_seam` carries "not a shared filament", "not a shared node", "not a permanent weld" and "disjoint filament populations", and its mechanism records "dynamic transient connection". |
| Positive | `tests/state/test_absent_connectors.py::test_the_scope_statement_says_what_it_is_and_what_it_is_not` | The scope statement reports kind `single-cell mechanochemical mechanics world`, `is_complete_cell_biology_world=False`, 19 absences, 36 registered connectors, 0 chemical-flux connectors, 13 resolved load paths, and names all three domain qualifiers in its closing condition. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_connectors_protrusion_motor.py::test_construction_rejects_bidirectional_false` | Constructing a contract with `bidirectional=False` raises. Both the contract type directly and `build_connector_contract` are exercised, so a gate removed from either is caught. |
| Negative (must fail) | `tests/state/test_connectors_protrusion_motor.py::test_validation_rejects_a_kinetic_connector_that_commits_early` | A `K` contract with `commit_on_accept=False` is refused. Constructed through the dataclass so the type's own gate is bypassed and the validator is the thing under test. |
| Negative (must fail) | `tests/state/test_connectors_protrusion_motor.py::test_validation_rejects_a_continuous_connector_that_claims_kinetic_state` | The mirror case: `C` with `commit_on_accept=True` is refused. |
| Negative (must fail) | `tests/state/test_absent_connectors.py::test_an_absence_with_no_reason_is_rejected` | An absence with an empty or token reason raises `AbsenceRecordError`. This is the control that makes the file a register rather than a comment block: an absence with no reason cannot be distinguished from a connector somebody forgot. |
| Negative (must fail) | `tests/state/test_absent_connectors.py::test_an_absence_with_no_reentry_condition_is_rejected` | Same for the re-entry condition. |
| Negative (must fail) | `tests/state/test_absent_connectors.py::test_absences_for_domain_refuses_an_unknown_domain` | A misspelled domain raises rather than returning an empty tuple, so a typo cannot read as "nothing is missing there". |

Each negative control has a positive complement in the same file —
`test_a_continuous_connector_without_kinetic_commit_is_accepted` and
`test_a_well_formed_absence_is_accepted` — because a validator that raised unconditionally would
satisfy every negative control above and the suite would report a gate that admits nothing.

## 10. Numerical and precision envelope

No floating-point arithmetic, no accumulation, no tolerance. Every assertion in both control files
is exact: string identity, tuple identity, `is True` / `is False` on bools, integer counts, and
exception types. There is no envelope to state and no input range outside which behaviour degrades,
because there is no numeric input.

The one quantity with a chosen threshold is `AbsentConnector.MIN_JUSTIFICATION_CHARS = 24`, and it
is a text-length floor rather than a numerical tolerance. It was chosen to reject the specific
non-answers a hurried author writes — `"n/a"`, `"none"`, `"scope"`, `"excluded"` — while admitting
a one-sentence reason. It is not tuned and it is not evidence of anything; a longer minimum would
buy nothing, since prose length is not a proxy for prose quality. Outside that, the record refuses
rather than degrading: there is no truncation, no default, and no silently accepted blank.

## 11. Production-backend residency and transfer

**This code will never run on the production backend.** It is a host-side declarative registry: two
modules of frozen dataclasses built once at import, resident in CPU memory, with no device array,
no kernel, and no per-step work. Nothing is transferred, in either direction, at any point in a
step. `float32`/`float64` do not arise.

A runtime that eventually evaluates one of these six connectors will hold its own device state and
will read these contracts once at assembly time to learn which endpoints it is allowed to touch.
That read is a build-time lookup, not a per-step host round-trip.

Neither module imports `numpy`, `scipy`, or any backend. Both are pure standard library
(`dataclasses`, `enum`, `typing`), which is also what keeps them importable under gate R1 with the
reference tree absent.

## 12. Comments and docstrings to discard

Nothing was carried, so the list below is what was **read and deliberately left behind** rather
than what was stripped from copied text:

- Provider module paths and package namespace in every form.
- Its runtime-status vocabulary — `SEAMED`, `KERNEL_BOUND`, `PI-GAP-gated`, `DECLARED ONLY`,
  `CENSUS-WIRED`, `magnitudes GAP` — and its gate labels. These are the provider's own
  self-assessment scale, and importing the scale would import the authority to grade with it.
- Dated in-code status assertions of the form "this landed on <date>", and PI-ratification notes
  naming decisions Aleph's PI has not taken.
- Internal document references (plan files, integration notes, lane names).
- Its finer `ConnectorFamily` taxonomy (`CONTACT`, `TRANSIENT_ACTIN`, `ACTIN_ANCHOR`,
  `IMMERSED_TRANSFER`) — see §3. Discarded as an inherited classification decision, not as noise.
- Its backend and device-residency assumptions: the reference's protrusion seam validates that
  every array is a CUDA array at construction. Aleph's registry is host-side and declares no
  residency at all, per §11.

**What replaces them:** new prose written for Aleph's situation, whose organising principle is that
every ownership rule is stated together with the failure it prevents and why that failure is
invisible. The seam entry does not only say "not a weld" — it says a shared array would transmit
load with no connector, produce plausible forces, and leave the binding state with nothing to bind,
so no test about forces could detect it. That sentence is the reason the rule exists, and it is
absent from both the appendix and the source.

The two `mechanical_interpretation` fields for the rear seam and the bundle root are the load-bearing
prose in this entry, and their contents are asserted by I3 and I4 rather than left to review.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Both control files green on 2026-07-30: `tests/state/test_connectors_protrusion_motor.py` and `tests/state/test_absent_connectors.py`, 177 tests passed, 0 failed, all assertions exact. Status remains `PROPOSED`, not `ACCEPTED`. |
| Reviewer | **Agent-proposed, unratified.** Written by the L21 lane agent during autonomous operation; no human has reviewed it. The audit verdicts in §6 were produced by reading the reference tree and following its call graph, which is evidence about *that* tree and not about this code. |
| Rollback | Delete `aleph/state/connectors_protrusion_motor.py` and `aleph/state/absent_connectors.py` with their two test files. Nothing imports either module today, so nothing else breaks. What is lost is the ownership discipline for the rear seam and the bundle root, and the queryable scope limits — neither of which any other module currently carries, and both of which are silent when absent. |

Why `PROPOSED` and not `ACCEPTED`: `ACCEPTED` would assert that the port is qualified, and there is
no oracle here to qualify it against (§7). A registry of declarations passing tests about its own
declarations is a consistency result, not a correctness one.

## 14. Honest limits

What this entry does **not** establish:

- **No physics.** No force is computed, no load path is evaluated, nothing is integrated. Six
  contracts declared is not six connectors working, and Aleph evaluates none of them. That is the
  point of `implemented=False` and of `declared_but_unimplemented()`.
- **The mechanisms are recorded, not derived.** The Brownian-ratchet mechanism is stated as the
  mechanism these two connectors must implement. No ratchet law is written down here, no
  force-velocity relation is derived, and no literature was read for this lane. When a runtime
  lands it needs its own ledger entry with a real oracle — a stall-force limit, a zero-load
  velocity, or a fluctuation-dissipation check.
- **`INHERITED_UNVERIFIED`: the eight-seams verdict rests on a static read.** It was established by
  grep and call-graph following, not by running the reference. A concrete
  `ProtrusionGraphConnector` constructed dynamically — by name, through a factory, or in a script
  outside the paths searched — would not appear in a grep for `accumulate_actor`. The 05:05 ports
  lane audited its ten candidates by *running* the reference, which is stronger; this entry did
  not, and the verdict should be read at that lower class. What is solid regardless is the
  positive half: `ProtrusionStepBindings` is constructed nowhere outside the reference's own test
  file, so the loop that calls the connectors cannot execute in production.
- **The manuscript appendix is not authority.** PLAN §1.3: use the manuscript for science, never
  for authority. It self-describes as the 29 July seed scope requiring reconciliation before
  citation. The 19 absences and the scope statement are that seed scope, faithfully registered.
  Whether Aleph's *ratified* scope has the same 19 absences is a PI decision this entry cannot
  make and does not claim to have made.
- **The 36-connector count is inherited, not counted by Aleph.** `REGISTERED_CONNECTOR_COUNT = 36`
  comes from the appendix. Six lanes are registering the 36 concurrently; until every one has
  landed, no test in the tree adds them up. The `CHEMICAL_FLUX`-is-zero claim is therefore verified
  for the ten contracts this lane owns and *asserted* for the other 26. A cross-lane test that
  sums the registries and re-checks the family count is the missing control, and it cannot be
  written from inside one lane.
- **The absence register cannot enforce itself where it matters most.** `require_in_scope` refuses
  only when a caller chooses to route through it. Nothing stops a caller from computing an
  organelle-specific answer and never asking. Making the refusal unavoidable needs an observation
  layer that consults the register before it emits a claim; that layer does not exist.
- **The 19 absences are not exhaustive.** They are the 19 the appendix lists. A connector absent
  from both the registry and this register is invisible to every check here, and no procedure in
  this lane could have found one.
- **No GPU, and no GPU is possible.** §11: nothing here has a device side.
