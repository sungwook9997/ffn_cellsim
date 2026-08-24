# ALEPH-PORT-2001 — cytoskeletal cross-system connector contracts (group C, 8 kinetic edges)

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2001` |
| Lane | `L20 cross-system and cytosol connector contracts` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

```python
from aleph.state.connectors_crosssystem_cytosol import (
    CYTOSKELETAL_CROSS_SYSTEM,      # the eight group-C contracts, in manifest order
    LINC_CONNECTORS,                # the three nuclear-envelope edge names
    PRIOR_ART_NOTES,                # per-connector basis for the audit verdict
    PRIOR_ART_REALIZATION,          # per-connector audit verdict
    PriorArtRealization,            # name_only / seam_only / kernel_bound_no_caller / …_reachable
    ConnectorContractError,
    build_contract,
    connector,
    connectors_touching,
    validate_contract,
)
```

Also exported by the same module but authorised by `ALEPH-PORT-2002`:
`CYTOSOL_IMMERSED_TRANSFER`. `LANE_CONNECTORS`, `SCHEMA_SOURCE` and `SCHEMA_MISMATCH` span both
entries. Nothing outside these lists is covered.

The eight contracts: `sf_cortex_transient`, `actin_cap_linc`, `mt_nucleus_linc`,
`mt_cortex_capture`, `if_nucleus_linc`, `if_sf_plectin`, `mt_sf_spectraplakin`,
`dorsal_arc_crosslink`. All eight are `family=K`, `commit_on_accept=True`, `bidirectional=True`,
`adjoint_required=True`, `implemented=False`. `dorsal_arc_crosslink` additionally sets
`internal_to="sf_arc"`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/contracts.py`, `ffn_sim/ac/engine/dispatch.py`, `ffn_sim/ac/engine/sf_mechanics.py`, `ffn_sim/ac/engine/stress_fiber.py`, `ffn_sim/ac/engine/sf_population.py`, `ffn_sim/ac/engine/microtubule_rig.py`, `ffn_sim/ac/engine/intermediate_filament_rig.py`, `ffn_sim/ac/engine/linc_connector.py`, `ffn_sim/ac/engine/cortex_population.py`, `ffn_sim/ac/engine/composed_native.py`, `ffn_sim/ac/engine/sf_motor_slice.py`, `ffn_sim/scripts/ac_synth_dump.py` |
| Source symbol(s) | `ConnectorContract`, `canonical_facade_claims`, `SF_CONNECTOR_BINDING_STATUS`, `SFInternalArcJointConnector`, `StressFiberStepBindings`, `MicrotubuleStepBindings`, `IntermediateFilamentStepBindings`, `MicrotubuleGraphConnector`, `IntermediateFilamentGraphConnector`, `LincJointPopulation`, `LincRigConnectorAdapter`, `build_linc_joint_population`, `validate_linc_architecture`, `LincFilamentRole`, `SFArcStateOwner`, `build_sf_arc_state_owner`, `compose_native_cell_world`, `_excluded_facade_noop`, `build_sf_motor_slice`, `build_connectors` |
| Read from | **working tree** — verified equal to the commit for every file cited below except where noted |
| Working tree == commit? | **yes for every file read.** `git status --short` reports 31 changes at this commit; none of them is a file cited in this entry. The modified engine files are `cortex_state.py`, `microtubule_rig.py`, `observe/artifact.py`; of these only `microtubule_rig.py` is cited here, and `git diff be0e5876 -- ffn_sim/ac/engine/microtubule_rig.py` was checked before any claim about it was written. |

> Reading the working tree rather than `git show` was a deliberate choice: the audit question is
> "what is on that machine right now", not "what was committed". The row above is what makes that
> answer reproducible.

**File digests of every source file read for this entry** (`shasum -a 256` of the working-tree copy,
2026-07-30). A commit identifies a revision; a digest identifies the bytes actually read, which is the
only provenance that survives a dirty tree:

```
sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72  ffn_sim/ac/engine/contracts.py
sha256:325307d38f7df0989c1bf71f421f5fb8ea9cf5341f3dd739e2ad890f70d2311a  ffn_sim/ac/engine/dispatch.py
sha256:e1bb8e05b47c3881534aa6c68983a4cb3ea81af6565f3aeaa46cd99a163474a7  ffn_sim/ac/engine/sf_mechanics.py
sha256:209982c9def60d7602d815a402dc74dc2ec60879a4bd6826eda5603fbfd3a990  ffn_sim/ac/engine/stress_fiber.py
sha256:4a6426848493fe0b7b8c52fae88adededede268bd16d1b602e0c0adfb1cd931c  ffn_sim/ac/engine/sf_population.py
sha256:8728a9060ea7a5eadc7803c176cd300c81fde41dc55bcd1e6ea04a13186c0d0e  ffn_sim/ac/engine/microtubule_rig.py
sha256:8560fe7dd8b499821732df15a5a83d9327cc93b3e7813ebc7945e5d4a338362a  ffn_sim/ac/engine/intermediate_filament_rig.py
sha256:5d3a01f80c55491b01e13c69542c6b85f04cc66157e986a31d77d74cbbd38093  ffn_sim/ac/engine/linc_connector.py
sha256:04d4542cf0c07a1f35bf1600eda4b40828e60e0ecebb60d71388ac363e903b70  ffn_sim/ac/engine/cortex_population.py
sha256:d3115b37bdb2d3110841661f557435c5ee8b5a14ae71c5874c4d8b39aeb9fc63  ffn_sim/ac/engine/composed_native.py
sha256:0812858af4006d5c92f087e166e9c3ce273d09614f1a1a2661fc489d45fd55d7  ffn_sim/ac/engine/sf_motor_slice.py
sha256:a9fd194329e64cd42285293d6aadffb8a413b1be5fea955b5fd78a7e487a6ed7  ffn_sim/scripts/ac_synth_dump.py
```

`microtubule_rig.py` is the one file above that differs from `be0e5876`; its digest is therefore the
working-tree digest and not the committed one. Every other digest is both.

**Correction to a standing assumption of this ledger.** The template's warning says a working-tree
read may be citing a revision it never read. In this case `HEAD` *is* `be0e5876`, so the commit and
the working tree share a base and the 31 changes are uncommitted deltas on top of it, not a different
revision. The distinction matters for `microtubule_rig.py` and for nothing else in this entry.

## 3. Why source-derived porting beats clean-room

**It does not, and no code was ported.** The eight declarations were re-derived from Aleph's own
manuscript extract `docs/manuscripts/extracted/APPENDIX_A_mechanics_registry.txt` (connector legend
plus group C), which is Aleph material, not provider material. The endpoint owners, endpoint roles and
mechanisms come from there; the mechanical interpretations, the ownership argument, the LINC
non-collapse argument and every line of validation were written here.

What the archive was read **for** is a different question, and the answer is the audit in §6: which of
these eight connectors has an implementation, which has only a typed seam, and which is a name. That
is a measurement about the predecessor, not a design input, and it is the reason this entry exists at
all rather than being folded into the module docstring.

Two things were *considered* for porting and rejected:

* **The connector-contract dataclass shape.** Any competent author converges on
  (name, family, endpoints, roles, mechanism, commit semantics); it is not a liftable asset. Aleph's
  shape is lane L15's, and it differs in substance anyway — L15 splits construction-time from
  finalisation-time validation, and drops the provider's `chemistry_card`, `generation_required`,
  `remap_on_accept` and `blocks_sleep_refine` fields, none of which this lane needs to declare.
* **The internal-joint kernel binding** (`SFInternalArcJointConnector`). Real, small, and correct as
  far as it goes — a Hookean joint over an index pair table. Rejected because Aleph has no stress-fibre
  component to bind it to, no ratified backend surface for it yet, and because a spring over a pair
  list is precisely the case where clean-room wins outright.

## 4. Physical or mathematical law represented

**There is no law here.** What crosses the boundary is an *ontology*: a statement about which pairs of
state owners may exchange force and through what kind of endpoint. That is a structural claim, and it
should be said plainly rather than dressed as physics.

The one structural theorem the group rests on is worth stating because it is what makes the
declarations falsifiable. Let two owners hold disjoint node sets `A` and `B`, with force accumulators
`f_A` and `f_B`. A connector evaluating a pair `(i in A, j in B)` is *adjoint* when the force it adds
satisfies

    f_A[i] += g,    f_B[j] += -g

for one shared `g`. Then for the isolated pair of owners

    sum(f_A) + sum(f_B) = 0

**identically**, at every configuration, converged or not — it is a statement about the scatter, not
about equilibrium. That is the check the ownership discipline buys: it holds if and only if the
scatter is two-sided, unsigned-consistent and not double counted, and it fails loudly otherwise.

Merge `A` and `B` into one array and the same sum is `0` whatever the connector does, because both
contributions land in one accumulator and cancel by construction. **The check survives but stops
testing anything.** That, and not tidiness, is why `sf_cortex_transient`, `if_sf_plectin` and
`mt_sf_spectraplakin` must join disjoint populations, and why `dorsal_arc_crosslink` — which genuinely
is intra-owner — has to declare `internal_to` so nobody later reads its trivially-zero closure as
evidence.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| (none) | — | — | This entry declares no numeric quantity. No stiffness, rate, length or force is asserted by any of the eight contracts. |

Deliberate: an endpoint role names *where* a connector may attach, never *how stiffly*. A stiffness
declared alongside an unimplemented connector would be an unsourced constant with nothing to check it.

Singular and boundary cases, each with the behaviour Aleph requires:

- **Both endpoints name the same owner.** Legal only with `internal_to` set to that owner. Otherwise
  refused — see `test_a_self_joining_edge_that_omits_internal_to_is_rejected`.
- **`internal_to` set on an inter-owner edge**, or naming an owner that is not the shared endpoint.
  Refused; the second case would send an ownership check to the wrong accumulator.
- **`bidirectional=False`.** Refused at construction by the contract type *and* by this lane's gate.
  There is no legal value other than `True`, so the field exists only to be asserted against.
- **`adjoint_required=False`.** Refused for every family. Group C is where it looks harmless, because a
  kinetic bond feels like it belongs to one side; it does not.
- **`commit_on_accept=False` on a `K` connector.** Refused. A rejected candidate step must leave no
  trace, so a bond committed during one is indistinguishable from a legitimate bond afterwards.
- **An unknown connector name passed to `connector()`.** Raises `KeyError`. Returning `None` would let
  a caller proceed with no load path, which is the failure this layer exists to prevent.

Invariants that must hold, each with the test that asserts it:

- **I1.** Exactly eight contracts, in manifest order —
  `test_group_c_is_the_eight_named_connectors_in_order`.
- **I2.** Every group-C connector is `family=K` — `test_every_group_c_connector_is_kinetic`.
- **I3.** Every group-C connector sets `commit_on_accept=True` —
  `test_every_group_c_connector_commits_only_on_accept`.
- **I4.** `dorsal_arc_crosslink` sets `internal_to="sf_arc"` and is the only internal edge here —
  `test_dorsal_arc_crosslink_is_internal_to_sf_arc`,
  `test_dorsal_arc_crosslink_is_the_only_internal_connector_here`.
- **I5.** Every other contract joins two distinct owners —
  `test_every_other_connector_joins_two_distinct_owners`.
- **I6.** The three LINC edges are distinct in filament owner, filament role, envelope role and
  mechanism — `test_their_filament_side_owners_are_three_different_components`,
  `test_their_filament_side_roles_are_all_different`,
  `test_their_envelope_side_roles_are_all_different`,
  `test_their_mechanisms_name_three_different_molecular_routes`.
- **I7.** `sf_cortex_transient` states the disjoint-population rule in its own text —
  `test_sf_cortex_transient_states_the_disjoint_population_rule`.
- **I8.** Every endpoint owner is one of the nine registered names —
  `test_every_owner_named_is_one_of_the_registered_nine`.
- **I9.** Nothing claims to be implemented —
  `test_nothing_in_this_lane_claims_to_be_implemented`.

## 6. Source evidence class and known retractions

**This section is the deliverable of the entry.** Verdicts were reached by following call paths, not by
reading status strings.

| connector | verdict | what it rests on |
|---|---|---|
| `sf_cortex_transient` | **SEAMED** | Typed slot `StressFiberStepBindings.cortex_transient: StressFiberGraphConnector`. No concrete implementer of `accumulate_rig` exists for it. `SF_CONNECTOR_BINDING_STATUS` self-labels it `"SEAMED (needs live cortex port + two-array adjoint scatter)"`. The endpoint data exists on both sides — `SFArcPopulation.cortex_local` / `cortex_sites`, and a `CortexEndpointDomain` entry — so the endpoints are materialised and the transfer is not. |
| `actin_cap_linc` | **KERNEL-BOUND, NO PRODUCTION CALLER** | `linc_connector.py` holds 7 real `@wp.kernel`s (`_linc_pair_force_kernel`, `_linc_native_socket_scatter_kernel`, `_linc_generation_validate_kernel`, `_linc_rollback_kernel`, `_linc_commit_kernel`, `_linc_ledger_reduce_kernel`, plus the `_linc_launch` indirection) and a concrete `LincRigConnectorAdapter`. **Verified unreachable:** the only importer of the module anywhere in the tree is `tests/ac/engine/test_linc_connector.py`; `LincJointPopulation(` is constructed at `linc_connector.py:1537` (its own builder) and at two places in that test file, nowhere else. Its own comment block records that bind/unbind candidate generation, reduced-Jacobian projection and the global-ledger merge are still open. |
| `mt_nucleus_linc` | **KERNEL-BOUND, NO PRODUCTION CALLER** | Same module, same kernels, same absence of a caller. `MicrotubuleStepBindings.nucleus_linc` is a typed slot with no concrete implementer wired in. |
| `mt_cortex_capture` | **SEAMED** | Typed slot `MicrotubuleStepBindings.cortex_capture`, validated at construction for name/endpoints, no concrete implementer. Cortex-side `CortexEndpointDomain` metadata exists; the capture kinetics do not. |
| `if_nucleus_linc` | **KERNEL-BOUND, NO PRODUCTION CALLER** | Same LINC module and same unreachable path. `IntermediateFilamentStepBindings.nucleus_linc` is an unfilled typed slot. |
| `if_sf_plectin` | **SEAMED** | Typed slot `IntermediateFilamentStepBindings.sf_plectin`. `SF_CONNECTOR_BINDING_STATUS` self-labels it `"SEAMED (needs intermediate_filament component; owned by the IF rig)"`. |
| `mt_sf_spectraplakin` | **SEAMED** | Typed slot `MicrotubuleStepBindings.sf_spectraplakin`. Self-labelled `"SEAMED (needs microtubule component; owned by the MT rig)"`. |
| `dorsal_arc_crosslink` | **KERNEL-BOUND AND REACHABLE** | `SFInternalArcJointConnector.accumulate` launches a real link-spring kernel over an `(J, 2)` joint table. Reachable: `build_sf_arc_state_owner(..., with_internal_arc=True)` binds it whenever the topology has arc joints, `SFArcStateOwner.accumulate` launches it, and that owner is constructed by `build_sf_motor_slice` — which two driver scripts (`scripts/ac_gate_b_sf_motor_native.py`, `scripts/ac_observe_sf_operator_native.py`) call — and by `build_native_composed_cell_world`, called by `scripts/ac_gate_b_composed_native.py`. Self-labelled `"KERNEL_BOUND"` and, unusually, the label is accurate. |

**Summary: 1 of 8 reachable, 3 of 8 kernel-bound but unreachable, 4 of 8 typed seams.** The earlier
audit's expectation of mostly bad news is confirmed and slightly sharpened: `dorsal_arc_crosslink` is
better than expected and is the *internal* edge, which is the opposite of a counterexample. It was
cheap precisely because no cross-owner scatter was needed.

**Retractions and self-labelled gaps found (looked in `STATE.md`, `STATE_NONQUOTABLE.md`,
`outputs/tag_kb/run_audit_report.md`, `outputs/ac/cell_assembled/composed_world_v2.json`, the module
docstrings and the in-code status dicts):**

* `outputs/ac/cell_assembled/composed_world_v2.json` self-labels `evidence: CENSUS-WIRED`,
  `quantitative_claim_status: BLOCKED`. No number from any of these eight is quotable.
* `SF_CONNECTOR_BINDING_STATUS` is an in-code retraction and is honest about three of the four SF
  edges. Its accuracy was checked against the call graph rather than trusted, and it held.
* The LINC module's own comment block names three open gaps. Also accurate.

**Two defects found that are not self-labelled:**

1. **The facades that claim six of these eight edges are constructed nowhere outside tests.**
   `canonical_facade_claims()` assigns the group-C edges to `StressFiberActor`, `MicrotubuleRig` and
   `IntermediateFilamentRig`. None of those three classes is instantiated anywhere outside
   `ffn_sim/tests/`. In the one composed native path, `compose_native_cell_world` substitutes
   `_excluded_facade_noop` for the microtubule, IF, protrusion, surface and fluid facade types and
   drives only `SFArcStateOwner.accumulate()` — which launches component-local mechanics plus the
   internal arc joint, and no group-C inter-owner edge. So the exact-once dispatch-coverage gate is
   satisfied by classes that production never builds. This is the same shape as the previously
   recorded `membrane_cortex_contact` finding, one layer up: there the method existed and had no slot,
   here the class exists and has no caller.
2. **`scripts/ac_synth_dump.py::build_connectors` writes wrong endpoints while claiming they are
   real.** Its docstring says "the TOPOLOGY (component names, roles, connector families, endpoint
   components) is the REAL `reference_cell_architecture()` graph". Three of its fourteen rows
   contradict the architecture it names: `sf_cytosol_transfer` is given `("sf_arc", "nucleus")` where
   the contract says `sf_arc <-> cytosol`; `nmii_sf_motor` is given `("sf_arc", "cortex")` where the
   contract says `nmii <-> sf_arc`; `integrin_collagen_clutch` is given `("cortex", "ecm")` where the
   contract says `focal_adhesion <-> ecm`. Nothing validates the table against the architecture — the
   loop does `if a not in by or b not in by: continue`, so a wrong or absent owner is silently
   skipped. The output is written through the same writer real dumps use and is byte-schema-identical
   to one, so a viewer renders joints between the wrong owners under a correct connector name. The
   file is untouched in the working tree, so this is a defect at the commit.

**Reachability from the source's production path:** three of the eight (the LINC trio) are dead code
by the strict definition — nothing but tests reaches them. Four are interfaces with no implementer.
One is live.

**Live tests in the source, and whether they exercise the law:** `tests/ac/engine/test_linc_connector.py`
is substantial and does exercise force and adjointness against an in-module reference. But an
in-module reference is not an independent oracle, and the code under test is on a path nothing else
takes. The other seven are covered by contract/architecture-validation tests and by binding-status
assertions (`test_sf_mechanics.py` asserts the strings above start with `"SEAMED"`), i.e. tests of the
plumbing and of the honesty labels — not of any law.

## 7. Independent oracle or derivation

For declarations, the oracle is Aleph's own manuscript extract plus internal consistency, and both are
checkable without the archive:

* **Aleph's manuscript extract** (`APPENDIX_A_mechanics_registry.txt`, connector legend and group C)
  fixes the eight names, their endpoint owners, their endpoint roles and their mechanisms. The tests
  pin all four per connector, so a declaration drifting from the manuscript fails.
* **The adjoint-closure identity of §4** is the mathematical statement that gives the ownership rule
  teeth. It is not exercisable here because nothing is implemented; it is recorded as the check any
  future implementation must satisfy, and `ALEPH-PORT-2002` §7 says the same for group D.
* **Internal consistency** is machine-checked: family determines `commit_on_accept` for all fourteen,
  endpoint equality determines `internal_to`, every owner name is drawn from the registered nine, and
  the three LINC edges differ pairwise on four independent attributes.

Agreement with the archive is explicitly **not** the check. Where the archive and this module differ
they differ on purpose — it carries a `chemistry_card` string per edge and Aleph does not, because a
card naming a mechanism that has no rate law attached is a declaration that reads like a citation.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_connectors_crosssystem_cytosol.py::TestAllFourteenAreConstructible::test_group_c_is_the_eight_named_connectors_in_order` | The eight names appear exactly once each, in manifest order. |
| Positive | `…::TestCommitSemantics::test_every_group_c_connector_commits_only_on_accept` | All eight set `commit_on_accept=True`. |
| Positive | `…::TestOwnershipDiscipline::test_dorsal_arc_crosslink_is_internal_to_sf_arc` | `internal_to == "sf_arc"` and both endpoints are `sf_arc`. |
| Positive | `…::TestTheThreeLincEdgesStayDistinct::test_their_envelope_side_roles_are_all_different` | Three distinct LINC-facing envelope roles, each naming LINC. |
| Positive | `…::TestBidirectionalCannotBeFalse::test_all_fourteen_pass_the_lane_gate` | `validate_contract` accepts all fourteen — the control that stops a gate which rejects everything from looking correct. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_connectors_crosssystem_cytosol.py::TestBidirectionalCannotBeFalse::test_the_contract_type_itself_rejects_bidirectional_false` | `bidirectional=False` is refused at construction. Fails if the contract type stops refusing. |
| Negative (must fail) | `…::TestBidirectionalCannotBeFalse::test_the_lane_gate_also_rejects_bidirectional_false` | The lane gate refuses it independently, probed on a duck-typed stub so construction cannot mask it. |
| Negative (must fail) | `…::TestBidirectionalCannotBeFalse::test_the_lane_gate_rejects_adjoint_required_false` | `adjoint_required=False` is refused. |
| Negative (must fail) | `…::TestValidationCanFail::test_a_kinetic_connector_that_commits_before_acceptance_is_rejected` | A `K` connector with `commit_on_accept=False` is refused. |
| Negative (must fail) | `…::TestValidationCanFail::test_a_self_joining_edge_that_omits_internal_to_is_rejected` | The rule keeping `dorsal_arc_crosslink` honest, probed on a synthetic edge so the eight real contracts cannot satisfy it by accident. |
| Negative (must fail) | `…::TestValidationCanFail::test_internal_to_naming_the_wrong_owner_is_rejected` | `internal_to` pointing at a component that is not the shared endpoint is refused. |

The negative controls are probed on a synthetic `probe_edge`, not on a real contract, and
`test_a_valid_probe_declaration_constructs` establishes that the probe is well formed. So each
rejection differs from an accepted declaration by exactly one field, which is what makes it a control
rather than a coincidence.

## 10. Numerical and precision envelope

No arithmetic. Working precision, accumulation precision, conditioning and tolerance are all
inapplicable: every assertion in this entry is an exact comparison of strings, booleans, enum members
and tuple order. There is no tolerance to loosen and therefore no envelope to leave.

The one thing that *is* exact and worth naming: the tests compare tuple order with `==`, not set
membership, so inserting a connector in the wrong position fails. Order is part of the manifest.

Outside the declared vocabulary the behaviour is refusal, not degradation: an unknown name raises
`KeyError`, an out-of-vocabulary family value raises at construction, and an inconsistent
`internal_to` raises `ConnectorContractError`.

## 11. Production-backend residency and transfer

**This code never runs on the production backend.** It is a declaration registry: frozen dataclasses
of strings, booleans and enum members, resident on the host, read at build time to decide what must be
dispatched and to refuse what may not be. No device array is allocated, nothing is transferred, no
host round-trip occurs per step, and there is no per-step cost of any kind.

That is a residency answer, and it is the one that keeps this module importable on a machine with no
GPU and no `warp` installed — which is the condition every test in this lane runs under. When the
eight connectors are eventually implemented, the implementations will have a residency answer of their
own and will need their own ledger entries; this entry authorises none of them.

## 12. Comments and docstrings to discard

Nothing from the provider survives. Specifically discarded:

* **Provider module paths and symbol references.** Every `ffn_sim.ac.engine.*` path, every
  `:class:`/`:func:` cross-reference into that tree, every `ac/engine` and `ac/cell` mention. They live
  in §2 and §6 of this entry and nowhere in `aleph/**`.
* **Provider status vocabulary as code.** `SEAMED`, `KERNEL_BOUND`, `PI-GAP`, `CENSUS-WIRED`,
  `LANDED`, `verify:binding hygiene flag`, and gate labels. Replaced by Aleph's own
  `PriorArtRealization` enum whose four members are defined in terms of *what was found* — a name, an
  interface, kernels with no caller, kernels with a caller — rather than by inheriting a label.
* **Dated authority claims** such as "PI-ratified 2026-07-25", "PI D5-A, 2026-07-28", "native-validated
  2026-07-25". Those are the provider's decision record, not Aleph's; Aleph's decisions live in
  `docs/decisions/`.
* **`chemistry_card` strings** (`"transient_actin_crosslink"`, `"actin_talin_integrin_nascent"`, …).
  Deliberately not carried: a card naming a mechanism with no rate law attached reads like a citation.
  The mechanism is stated in prose in `mechanism` instead.
* **Sourced constants.** No stiffness, persistence length or rate crosses this boundary. The α-actinin
  crosslink stiffness and the actin bending modulus appear in the archived internal-joint code; they
  are not inherited here, and any future use of either needs its own ledger entry with its own
  literature citation.
* **Plan-document and branch references** (`SURFACE_BODY_PLAN.md:257`, `T2`, `T10`, `Card-5`, `LC-3`,
  `N1-N9`, `plan gap 2`/`plan gap 5`). Unreadable without the archive open, which is the second half
  of the porting rule and the half usually skipped.

Replaced by: a module docstring that derives the ownership argument, the LINC non-collapse argument
and the commit-semantics argument from first principles; per-contract `mechanical_interpretation`
prose written here; and `PRIOR_ART_NOTES`, which states the audit basis in Aleph's vocabulary and
points at this entry for the paths.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED`. As of 2026-07-30, `tests/state/test_connectors_crosssystem_cytosol.py` passes 61/61 covering both this entry and `ALEPH-PORT-2002`. That is evidence the declarations are self-consistent and that the validation can fail; it is not evidence that any connector works, because none is implemented. |
| Reviewer | **Agent-proposed, unratified.** Written by lane L20 during autonomous operation. No PI review. |
| Rollback | Delete `aleph/state/connectors_crosssystem_cytosol.py` and `tests/state/test_connectors_crosssystem_cytosol.py`. Nothing breaks: nothing imports either as of this writing. The manifest loses eight of its thirty-six connector declarations, which blocks a state-schema completeness check and therefore blocks the representation lanes' axis-order work — the exact dependency `ALEPH-PD-002` identified. |

## 14. Honest limits

What this entry does **not** establish:

* **No physics.** Not one of the eight transfers a force in Aleph. `implemented=False` is the whole
  status. A green test suite here means the declarations are consistent, and consistency is not
  correctness.
* **No stiffness, rate, or magnitude is declared**, so nothing here constrains a number. Every
  future implementation needs its own sourced constants and its own controls.
* **The adjoint-closure identity of §4 is unexercised.** It is derived and it is the right check, and
  no test in this lane runs it, because there is nothing to run it against. Carried as `UNVERIFIED`.
* **The audit in §6 is a reading of a call graph, not an execution.** It was done by grep and by
  following constructor call sites on a CPU-only Mac. Two consequences: a runtime-registered
  implementation (a plugin, a dynamic import, an entry point) would be invisible to it — none was
  found, but absence by grep is weaker than absence by running; and the reachability claim for
  `dorsal_arc_crosslink` says a driver script constructs the owner, **not** that the kernel was
  observed to launch. That would need a GPU, and this session runs zero GPU jobs.
* **The three LINC verdicts are one verdict, not three.** All three edges share one module, so
  "kernel-bound, no caller" is a property of that module. If the three were implemented separately the
  verdicts could diverge.
* **`microtubule_rig.py` is modified in the working tree.** The diff against `be0e5876` was checked
  before any claim about it was made, but three of the eight verdicts touch that file and a reader
  reproducing this on a differently-dirty tree may see something else.
* **Whether the shared state schema is the live contract type is a runtime fact, not a claim of this
  entry.** `SCHEMA_SOURCE` reports it. At the time of writing lane L15's shared state-schema module had
  landed and is the live type; the local fallback of identical shape remains for the case where it is
  not, and a test refuses a silent fallback. **This entry's Aleph target is one file only** — the
  connector module named in §1 — and it authorises nothing in the shared schema, which is L15's to
  ledger. Two other lanes' census modules failed to import against L15's `citation_status` vocabulary
  at that moment, which is a live cross-lane disagreement this entry does not resolve and does not
  depend on.
* **`composite_group` is `None` for all eight**, which asserts that none of them is a member of a
  series joint. That follows from the manuscript extract, which places the only composite group in
  another connector block. If that changes, this claim is wrong and the tests will not catch it,
  because they assert the current answer.
* **Nothing here is evidence about `ALEPH-DQ-104`.** The commit-on-accept semantics assume an
  acceptance predicate exists; which predicate is the PI's reserved decision.
