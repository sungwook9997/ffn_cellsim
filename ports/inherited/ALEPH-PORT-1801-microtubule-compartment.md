# ALEPH-PORT-1801 — microtubule compartment contract, and the audit of the reference rod rig

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1801` |
| Lane | `L18 census — internal frame / internal fluid / nucleus` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Zero lines taken. Two ideas recorded as prior art, both re-derivable in an afternoon. |

---

## 1. Aleph API

This entry authorises the `microtubule` and `centrosome_mtoc` entries of:

```python
from aleph.state.census_frame_fluid_nucleus import (
    FRAME_FLUID_NUCLEUS_COMPONENTS,
    INTERNAL_FRAME_COMPONENTS,
    by_name,
    validate_registry,
)
```

Nothing outside those two registry entries and the validator that checks them is covered here. In
particular **no microtubule mechanics is authorised by this entry** — there is no rod, no bending
law, no dynamic-instability kinetics, and no MTOC anchor in `aleph/`. What lands is a contract that
says what a microtubule component would own and what it would not be allowed to claim.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/microtubule_rig.py`, `ffn_sim/ac/solid/microtubule.py`, `ffn_sim/ff/microtubule.py` |
| Source symbol(s) | `STATIC_ADAPTER_STATUS`, `LegacyMicrotubuleBendingAdapter`, `_mtoc_anchor_force_kernel`, `mtoc_anchor_reference`, `MtocAnchorReference`, `_mt_di_proposal_kernel`, `_mt_di_stochastic_proposal_kernel`, `_mt_topology_commit_kernel`, `_mt_topology_rollback_kernel`, `MicrotubuleRodBuildSpec`, `MicrotubuleRodRuntime`, `DynamicInstabilityRateCard`, `MicrotubulePhase` |
| Read from | **BOTH.** `git show be0e5876:ffn_sim/ac/engine/microtubule_rig.py` and the working tree, then the diff between them. |
| Working tree == commit? | **NO.** `git diff be0e5876 -- ffn_sim/ac/engine/microtubule_rig.py` reports 10 insertions / 2 deletions, 1311 → 1319 lines. The other two paths are clean. |

### 2a. What the divergence actually is, since it changes an earlier audit finding

The entire diff is one comment block. The committed revision says the dynamic-instability RNG
transitions are a deferred kinetic layer; the working tree replaces that with a correction saying
the transitions are implemented and stochastic, and that what remains deferred is only their **load
dependence** (the transition-rate kernel takes scalar rates and no force array).

Reading the code rather than either comment: the working tree is right. A per-microtubule uniform
draw switches phase with probability `1 - exp(-f dt)`, and the switch is proposed into candidate
arrays with separate commit and rollback kernels. So a prior audit finding reached me as
*"microtubule dynamics there is a name, not a behaviour"*, and that is **too strong** and is
corrected here:

- **The plus-end kinetics are a behaviour.** Stochastic catastrophe/rescue with commit/rollback
  exists, in kernels, on candidate arrays.
- **The bending is a name.** `STATIC_ADAPTER_STATUS = "STATIC_BENDING_ADAPTER_DYNAMIC_TOPOLOGY_PENDING"`
  is a module constant, is the default `status` of the state owner, is returned by the bending
  adapter's `status` property, and is re-exported. Present in **both** revisions at the same four
  line numbers, so it is not a working-tree artefact. The bending adapter wraps a fixed bead-chain
  aster: it is validated at a frozen topology, and its own label says the dynamic-topology case is
  pending. Bending and topology are therefore not consistent with each other in that tree.
- **The rate values are unset by design**, labelled priors to be inferred rather than constants to
  be signed. Aleph agrees with that stance and inherits neither the values nor the label.

An audit finding that arrived at the right verdict through a wrong reason is worth correcting even
when the verdict does not move, which is why this subsection is longer than the verdict.

## 3. Why source-derived porting beats clean-room

**It does not. This is `RE-DERIVED`, and there is nothing here that clean-room does not win.**

Weighed and rejected as porting reasons:

- *The MTOC anchor.* A Hookean spring between a rod base and a hub point, gather and scatter on the
  same two endpoints. Any competent author writes this correctly the first time. Per the template's
  own rule, a convention anyone would arrive at independently is prior art, not a port.
- *The candidate/commit/rollback split for topology changes.* Genuinely good, and Aleph already has
  it: `ALEPH-PORT-305` (accepted-step transaction) and `ALEPH-PORT-303` (evaluation-proof coverage)
  cover the same ground on Aleph's own derivation, and `ALEPH-PORT-303` was hardened *against* the
  defect that tree's coverage gate has. Taking it now would be importing an architecture Aleph has
  already re-derived better.
- *The 1319-line rig itself.* It is an ownership seam over four Protocols plus five kernels, written
  against a connector-contract system Aleph deliberately does not share. Master plan §5.3 rejects
  inheriting that census-as-schema structure, and `ALEPH-PD-002` replaces it with registered data.
  Porting the seam would re-import the thing the decision removed.

**Recorded as prior art, unported:** (i) a `topology_generation` counter, so a length change is an
observable committed event rather than a quiet array resize — Aleph registers the state slot and
derives its semantics from its own transaction contract; (ii) keeping transition rates out of the
component and in an injected card, which is the same conclusion Aleph's `require_sourced` policy
reaches from the other direction.

## 4. Physical or mathematical law represented

**None, and saying so is the point.** What crosses the boundary in this entry is neither a law nor
code. It is a **scope argument**: which microtubule facts are owned state, which are internal to the
owner, and which are claims the representation cannot support. The template asks for this case to be
named explicitly, so it is named: this is a contract.

Two derivations Aleph does own, stated here because the contract's fields are meaningless without
them:

**(a) Why the MTOC is internal and not a component.** A component in this ontology is a legal
connector endpoint. The MTOC is loaded only by the minus ends of rods that anchor to it; no other
compartment reaches it except through those rods. Its degrees of freedom are therefore not
independently addressable from outside the microtubule owner, and a connector to it would have to
route through the owner anyway. `INTERNAL` is not a demotion — it is the statement that the hub has
no external endpoint, which is checked by `validate_registry` rule R6/R8.

**(b) Why `mtoc_torque` is registered *and* refused.** Take a hub at `X` and a rod base at `x`, with
an anchor potential `U(|X - x|)`. The force on the base is `-dU/dx`, the reaction on the hub is its
negative, and both act along `X - x`. The moment of that pair about the hub is `(x - X) x f`, and `f`
is parallel to `(X - x)`, so the cross product vanishes identically — for *any* central `U`, not
just a Hookean one, and at any extension. A single-point hub anchored by central forces can
therefore never acquire torque: the array can exist and can only contain zeros. Registering the
slot without registering that fact is exactly how a zero-filled array becomes a reported result, so
the contract carries both. The same argument applies to the rotational half of `mtoc_pose`.

The reference confirms the consequence independently: `torque` appears **zero times** in
`microtubule_rig.py` at either revision, and the anchor kernel's own note states the pair is
moment-free. So the seed registry declares MTOC torque as owned state and no implementation
anywhere populates it. That is a registry-versus-reality gap, and it is now written down where a
reader of the contract meets it.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `name`, `group`, `parent` | — (identifier) | — | lower snake_case; `parent` resolves to a registered `E` entry |
| `scope` | — (`ScopeTag`) | — | one of `E I H B X` |
| `owned_state` | — (state-slot names) | — | non-empty for `E`, empty for `H`/`B`/`X` |
| `unsupported_claims` | — (prose) | — | non-empty for every entry in this group |

This entry registers no dimensional quantity. That is deliberate: a census entry that carried a
stiffness would be a parameter card wearing a contract's name, and the value would then be inherited
without a source. Dimensional quantities belong to `aleph/units` (`ALEPH-PORT-101`) under
`require_sourced`.

Singular and boundary cases, with the behaviour Aleph requires:

- A microtubule entry with empty `owned_state`: refused (R5). An explicit owner that owns nothing is
  a placeholder a connector would later be attached to.
- `centrosome_mtoc` with no parent, or a parent that is not a registered `E` entry: refused (R6/R8).
  A hub with no owner is a top-level component wearing the wrong tag.
- `centrosome_mtoc` naming itself as parent: refused (R11), the one-cycle no other rule catches.
- An unknown group name or an unknown component name: `KeyError` from `in_group` / `by_name`, never
  `None`. A missing component is a scope question and the caller must not walk past it on a falsy
  value.

Invariants, each with the test that asserts it:

- **I1.** `microtubule` is `E`, owns seven state slots, and has no parent —
  `tests/state/test_census_frame_fluid_nucleus.py::test_every_entry_carries_its_expected_scope_tag`
  and `::test_only_explicit_entries_own_state`.
- **I2.** `centrosome_mtoc` is `I`, parented to `microtubule`, and owns exactly the five declared
  slots — `::test_the_mtoc_is_an_internal_hub_of_the_microtubule_owner`.
- **I3.** The torque and orientation channels are registered as unsupported, with the moment-free
  reason attached — `::test_the_mtoc_declares_its_torque_and_orientation_channels_as_unsupported`.
- **I4.** Every parent in the registry resolves to a registered `E` entry —
  `::test_every_internal_entry_names_a_registered_explicit_parent`.
- **I5.** Contracts are frozen, so a scope cannot be edited at runtime —
  `::test_entries_are_frozen_so_a_registry_cannot_be_edited_in_place`.

## 6. Source evidence class and known retractions

- **Self-declared status.** `STATIC_ADAPTER_STATUS` is a status label the source ships in its own
  public surface. A module that exports the string `"..._PENDING"` as part of its API is telling the
  reader it is incomplete; that is the source's own evidence class for the bending path, and it is
  low.
- **A retraction found, and it is in the working tree.** The comment correction described in §2a is
  a retraction of the committed revision's claim that the RNG transitions were deferred. It is
  uncommitted, so a reader of `be0e5876` alone gets the retracted claim. This is precisely the
  hazard the "read from" row exists for, and it fired on the first candidate this lane examined.
- **Reachability.** `MicrotubuleRodRuntime.__init__` raises unless the device is CUDA, so nothing in
  the rod runtime is reachable on a CPU-only host. The host-side oracles (`mtoc_anchor_reference`,
  `dynamic_instability_reference`, `dynamic_instability_stochastic_reference`) are reachable and are
  what that tree's own tests exercise.
- **Live tests, and what they exercise.** `ffn_sim/tests/ac/engine/test_microtubule_rig.py` exists
  and, on a CPU host, can only be testing the host oracles and the structural seams — the kernels
  cannot launch. So its green status is evidence about the plumbing and the oracles, and not about
  the device path. Aleph draws no conclusion from it either way, having run nothing.
- **Where I looked for retractions:** the committed-vs-working-tree diff, the module's own status
  constants and comment blocks, `ffn_sim/ac/solid/microtubule.py`, `ffn_sim/ff/microtubule.py`, and
  the test file's name and location. I did **not** read that tree's planning documents, and no claim
  here rests on them.

## 7. Independent oracle or derivation

The contract's oracle is the **scope discipline itself**, and it is not the source in any part:

1. The moment-free derivation in §4(b) is a two-line exact argument from the definition of a central
   force. It refutes a registry claim without reference to any implementation, and it is the reason
   `mtoc_torque` carries a refusal.
2. `validate_registry` is checked against a **local stand-in object that performs no validation of
   its own** (`::test_the_validator_itself_refuses_independently_of_the_contract_type`), so the
   refusals are proven to belong to the validator rather than to whichever dataclass supplied the
   fields. Its premise is separately asserted by
   `::test_the_stand_in_object_really_performs_no_validation_of_its_own`.
3. Agreement with `ffn_cellsim` is not used as evidence anywhere in this entry. Where the two
   disagree — MTOC torque — this entry sides with the derivation and records that the registry claim
   is unsupported.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_the_registry_satisfies_the_scope_discipline` | the shipped 16-entry registry passes all eleven rules |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_the_mtoc_is_an_internal_hub_of_the_microtubule_owner` | `I`, parent `microtubule`, exactly the five declared state slots |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_the_mtoc_declares_its_torque_and_orientation_channels_as_unsupported` | the moment-free refusal and the moment-arm re-entry condition are present |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_every_internal_entry_names_a_registered_explicit_parent` | all six `I` entries resolve to one of the four registered `E` owners |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_an_internal_entry_with_no_parent` | R6: an internal substructure with no owner is refused, not defaulted to a top-level component |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_parent_that_is_not_a_registered_explicit_component` | R8, both ways: an unknown parent, and a parent that exists but is itself internal |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_an_explicit_entry_that_names_a_parent` | R7: an owner nested inside another owner is refused |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_self_parented_entry` | R11: a one-cycle, which no other rule catches |

Each builds its malformed input with `dataclasses.replace` on a **known-good** entry, so exactly one
field differs and the test cannot pass by tripping a different rule than the one it names.

## 10. Numerical and precision envelope

No floating-point arithmetic is introduced by this entry, so the usual envelope is empty and the
honest statement is why:

- **Working precision:** not applicable. The module holds strings, tuples of strings, an enum and
  `None`. There is no accumulation, no conditioning question, and no tolerance to choose.
- **Tolerances in the controls:** every assertion is exact — set equality, identity comparison on
  enum members, and substring containment. A tolerance would be a defect here: a contract that
  matched its declared scope "to within 1e-9" would mean nothing.
- **The numerical statement that *is* being made** is about a future implementation and is recorded
  as a refusal rather than as a number: the torque channel is exactly zero for a central-force
  anchor, algebraically and not to round-off, which is why the contract refuses the claim outright
  instead of bounding it.
- **Outside the envelope:** malformed input raises `CensusContractError` or `KeyError`. There is no
  degraded path and no default entry.

## 11. Production-backend residency and transfer

Host-side registry data, resident in the Python process for the life of the interpreter, immutable
(`frozen=True, slots=True`), and never transferred to a device. No per-step cost, because it takes
part in no step.

This is a residency answer and not an exemption: the contract is what a future backend implementation
is *written against*, so it must be readable without a device present. Constructing it must never
touch a GPU. That constraint is sharpened by something found during this audit and recorded in §6 of
`ALEPH-PORT-1803`: 48 modules in the reference call `wp.init()` at module scope, so importing one of
them initialises a GPU runtime as an import side effect. Under PLAN §0.1 an Aleph module that did
that would violate the zero-GPU policy by being imported. Nothing in this lane imports a backend.

## 12. Comments and docstrings to discard

No source text survives, because none was taken. Specifically discarded rather than translated:

- The status constant `STATIC_ADAPTER_STATUS` and its `"..._PENDING"` string. It is that tree's
  status vocabulary; Aleph's status vocabulary is the ports ledger and the evidence ladder.
- Internal plan-section references of the form "§4.3", gate labels, and lane names. Every one is a
  reference into a document Aleph does not have and must not depend on.
- Working-tree assumptions, including the corrected comment block itself. It is quoted **in this
  ledger entry** as audit evidence and appears nowhere in `aleph/`.
- The literature-proxy rate-card name attached to the dynamic-instability rates. Inheriting a
  citation label without reading the citation is the failure `require_sourced` exists to prevent.

**What replaces them:** prose written for Aleph's situation. The module docstring explains why the
MTOC is internal and why a subpopulation may own nothing while another internal entry owns five
slots; the contract fields carry the derivations of §4 in Aleph's vocabulary; and the moment-free
argument is stated from first principles rather than cited.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED` on 2026-07-30. The registry entries construct and 40 controls pass in `tests/state/`, but the shared state schema this data registers into was authored concurrently by another lane and was **not on disk** when this entry was written (see §14). Acceptance waits on that seam closing. |
| Reviewer | **Agent-proposed, unratified.** No PI review. The `RE-DERIVE` verdict and the correction to the earlier "dynamics is a name" finding are both agent judgements. |
| Rollback | Delete `aleph/state/census_frame_fluid_nucleus.py` and `tests/state/test_census_frame_fluid_nucleus.py`. Nothing imports either as of writing, so nothing breaks — which is also the honest measure of how much this entry currently establishes. |

## 14. Honest limits

- **The contract seam is unverified.** `aleph/state/schema.py` did not exist when this ran, so the
  module imports the contract type defensively and falls back to a local dataclass of identical
  shape. `SCHEMA_SOURCE` records which branch is live and reported `local_fallback`. Field names,
  order, types and defaults were matched to the specification given to this lane, **not** to code —
  so "identical shape" is `ASSUMED` until that module lands.
- **Nothing was executed in the reference.** The rod runtime is CUDA-only and this session ran zero
  GPU jobs by policy. Every statement in §2 and §6 comes from reading source and diffs. The
  strongest claims here — that `STATIC_ADAPTER_STATUS` is present at both revisions, that `torque`
  appears zero times — are textual facts and are stated as such; nothing is claimed about what the
  kernels compute.
- **`INHERITED_UNVERIFIED`:** none. No value, constant, or citation crosses from the reference.
- **Believed but not tested:** that a future Aleph microtubule component will find the seven
  registered state slots sufficient. A registry entry is a hypothesis about what an implementation
  will need, and this one has never been implemented against.
- **Needs a PI decision:** whether `mtoc_torque` and the rotational half of `mtoc_pose` should be
  registered at all. Registering a channel that no declared law can load is defensible — it makes
  the gap visible — and it is also exactly the shape of thing that later gets filled with a
  plausible number. The alternative is to drop both slots and let the re-entry condition carry the
  requirement. This lane chose visibility and flags the choice rather than settling it.
- **Not established:** anything about microtubule physics. No bending law, no dynamic-instability
  kinetics, no anchor force, and no oracle for any of the three exists in Aleph. This entry
  authorises two rows of registry data and a validator.
