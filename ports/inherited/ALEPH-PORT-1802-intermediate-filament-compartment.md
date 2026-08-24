# ALEPH-PORT-1802 — intermediate-filament compartment contract, and the audit of the reference cage

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1802` |
| Lane | `L18 census — internal frame / internal fluid / nucleus` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Zero lines taken. The one genuinely valuable thing in the source is a *negative* result, and a negative result is inherited by re-deriving the experiment, not by copying the code. |

---

## 1. Aleph API

This entry authorises the `intermediate_filament`, `keratin_population` and `vimentin_population`
entries of:

```python
from aleph.state.census_frame_fluid_nucleus import (
    INTERNAL_FRAME_COMPONENTS,
    by_name,
    validate_registry,
)
```

No intermediate-filament mechanics is authorised. There is no cable, no strain-stiffening law, no
crosslink kinetics and no material card in `aleph/`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/intermediate_filament_rig.py`, `ffn_sim/ff/intermediate_filaments.py`, `ffn_sim/ac/solid/intermediate_filament.py` |
| Source symbol(s) | `LEGACY_REFERENCE_STATUS`, `LegacyIntermediateFilamentCageAdapter`, `_LegacyIntermediateFilamentCage`, `IntermediateFilamentTransaction`, `IntermediateFilamentGraphConnector`, `IntermediateFilamentFluidTransfer`, `IF_LINEAR_TANGENT_REGIME`, `IF_WLC_STRAIN_STIFFENING_REGIME`, `LinearBackboneCableAdapter`, `NonlinearWlcCableAdapter`, `NonlinearCableCard`, `keratin_vimentin_backbone_stiffness`, `keratin_vimentin_persistence_length`, `IntermediateFilamentTurnoverCard`, `BackboneTurnoverKinetics`, `resolve_intermediate_filaments` |
| Read from | **working tree.** |
| Working tree == commit? | **YES** for all three paths — verified with `git diff be0e5876 -- <path>` on each; all report no change. |

The reference working tree carries 31 uncommitted changes overall, so this row is a verification and
not an assumption. Of the files this lane read across all four entries, exactly one differed from the
commit, and it is not one of these three (see `ALEPH-PORT-1801` §2a).

## 3. Why source-derived porting beats clean-room

**It does not. This is `RE-DERIVED`.** Three specific reasons, in descending order of how close the
call was:

1. **The nonlinear law that the registry declares is not implemented there.** The seed registry's
   numerical representation for this component is an "explicit nonlinear strain-stiffening cable
   graph". In the source, the nonlinear branch is a named-but-unfilled regime — the constant
   `IF_WLC_STRAIN_STIFFENING_REGIME` is literally the string
   `"WLC_STRAIN_STIFFENING_UNFOLDING_CARD_REQUIRED"` — and what runs is the **linear tangent** limit,
   with the deferral stated in the module's own prose and the strain-stiffening parameters withheld
   pending registration. So the thing worth porting is absent, and what is present is a Hookean
   spring. There is nothing to gain from a Hookean spring.
2. **The real cage is a spoke model, and its geometry is a finding rather than an asset.** The load
   path is a set of radial bead chains from just outside the nucleus to just inside the cortex, with
   a LINC bond at the inner bead and a cortex anchor at the outer bead. The source states plainly
   that the earlier *tangential perinuclear shell* did not couple the nucleus to the cortex at all
   and was decoration. That is a real and useful negative result about model topology — a
   perinuclear cage that surrounds the nucleus without radially connecting it transmits nothing —
   and it is inherited correctly by **re-deriving the test that detects it**, not by importing the
   spokes. A geometry is cheap to write; knowing which geometry is inert is what took the effort, and
   that fact now lives in this ledger entry.
3. **Confirmed: 999 lines, zero `@wp.kernel`, four Protocols.** An earlier audit reported this and it
   verifies exactly. `grep -c '@wp.kernel'` returns `0`; the four Protocol classes are
   `IntermediateFilamentTransaction`, `IntermediateFilamentGraphConnector`,
   `IntermediateFilamentFluidTransfer` and `_LegacyIntermediateFilamentCage`. Every device kernel it
   uses is imported from a shared network module — the file is an ownership seam over borrowed
   kernels. A seam written against a connector-contract system Aleph deliberately does not share
   (master plan §5.3, `ALEPH-PD-002`) has negative porting value: it would re-import the structure
   the decision removed.

**Recorded as prior art, unported:** (i) the radial-spoke-versus-tangential-shell negative result
above; (ii) keeping keratin and vimentin as distinct **material cards** under a single owner rather
than as two components — the source does this and Aleph's registry independently reaches the same
conclusion from the ownership rule in §4, so it is corroboration rather than inheritance;
(iii) refusing to fit the strain-stiffening parameters to make a coupling appear. That last stance
matches Aleph's `require_sourced` policy and PLAN §0.5, and it is the most creditable thing in the
file.

## 4. Physical or mathematical law represented

**No law. A scope argument and an ownership rule**, named explicitly as the template requires.

**Why two populations and one owner.** The census's global ownership rule is that every physical
filament belongs to exactly one component, and that co-location is never a mechanical connection.
Keratin and vimentin filaments differ in stiffness and persistence length but are mechanically
interchangeable in *kind*: both are tension-only cables, both crosslink through the same cytolinker
family, both terminate on the same LINC sites. Nothing in the load path distinguishes them
topologically. Two components would therefore duplicate a graph, duplicate three connectors, and
create an ownership question — which component owns a filament that crosslinks to both — with no
mechanical fact to settle it. One owner with per-population material cards has none of those
problems and loses nothing, because a material card is selected per filament by label.

The consequence the contract records is this: **a population is an identity, not a store.** Both
population entries are `INTERNAL` with `owned_state = ()`, while another internal entry in the same
group (`centrosome_mtoc`, `ALEPH-PORT-1801`) owns five slots. Both are internal; only one holds
state. Separating identity from storage is the mechanism that stops a subpopulation from quietly
becoming a second owner, and the registry now exhibits both cases side by side so the distinction is
demonstrated rather than asserted.

**Why the honest default is the linear tangent limit.** A strain-stiffening cable is soft below a
knee strain and stiff above it. The entire mechanical consequence of the nonlinearity is *where the
knee is*, so a nonlinear law with an assumed knee is not a nonlinear law — it is a linear law plus a
chosen crossover, and any stiffness read from it is the assumption restated. The contract therefore
registers the nonlinear representation as declared, and lists a specific refusal: no strain-stiffening
response may be reported at any strain until the material card is sourced **per population**, since
keratin and vimentin do not share one.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `name`, `group`, `parent` | — (identifier) | — | lower snake_case; `parent` resolves to a registered `E` entry |
| `scope` | — (`ScopeTag`) | — | one of `E I H B X` |
| `owned_state` | — (state-slot names) | — | non-empty for `E`; empty here for both `I` populations |
| `unsupported_claims` | — (prose) | — | non-empty for every entry |

No dimensional quantity is registered, deliberately. A persistence length or an areal modulus in a
census entry would be an unsourced parameter card wearing a contract's name; those belong in
`aleph/units` under `require_sourced` (`ALEPH-PORT-101`). This matters more here than elsewhere,
because the source's material cards are the one place a number would have been tempting to carry
across — and carrying it would be an unrecorded inheritance whether or not any code moved.

Singular and boundary cases:

- Either population without a parent, or with a parent that is not a registered `E` entry: refused
  (R6/R8).
- A population parented to the *other* population rather than to the owner: refused by R8, because
  the parent must be `EXPLICIT` and a population is not.
- `intermediate_filament` with empty `owned_state`: refused (R5).
- An `I` entry with empty `owned_state`: **allowed**, and that is intentional. R4 forbids state on
  `H`/`B`/`X`; it does not require state on `I`. An identity-only substructure is legitimate and
  this entry contains the two canonical examples.

Invariants, each with the test that asserts it:

- **I1.** `intermediate_filament` is `E`, owns six slots including `population_identity`, and has no
  parent — `tests/state/test_census_frame_fluid_nucleus.py::test_keratin_and_vimentin_are_distinct_identities_under_one_owner`.
- **I2.** Both populations are `I`, parented to `intermediate_filament`, own nothing, and each names
  a distinct material card — same test.
- **I3.** Neither population is a legal connector endpoint: neither appears in
  `EXPLICIT_COMPONENT_NAMES` — `::test_every_internal_entry_names_a_registered_explicit_parent`.
- **I4.** The owner refuses an unsourced nonlinear response and names the linear tangent limit as the
  honest default — `::test_the_intermediate_filament_owner_refuses_an_unsourced_nonlinear_response`.
- **I5.** Each population's re-entry condition requires a **second explicit owner**, so a distinct
  keratin network cannot be created by relabelling — `::test_keratin_and_vimentin_are_distinct_identities_under_one_owner`.

## 6. Source evidence class and known retractions

- **Self-declared status, twice.** `LEGACY_REFERENCE_STATUS = "MONOLITHIC_HOOKEAN_REFERENCE_NONLINEAR_TURNOVER_PENDING"`
  and `IF_WLC_STRAIN_STIFFENING_REGIME = "WLC_STRAIN_STIFFENING_UNFOLDING_CARD_REQUIRED"`. Both are
  exported constants. The source is telling its reader, in its public surface, that the nonlinear
  law is unfilled and the cage it wraps is Hookean. That is the source's own evidence class for this
  component and it is low, honestly labelled.
- **Provisional parameters, self-flagged.** The cage module marks its filament count, its
  crosslink/backbone stiffness ratio, its LINC stiffness and its anchor stiffness as provisional and
  pending registration, and states that a weak coupling at the anchored backbone stiffness is a
  finding rather than a knob. Aleph inherits none of the four values. Two persistence lengths and two
  moduli in the seam module do carry literature attributions; those citations were **not read** by
  this lane and no value crosses, so nothing here is `INHERITED_UNVERIFIED` — it is simply absent.
- **A retraction found.** The tangential perinuclear shell is described as not coupling the nucleus
  to the cortex and as decoration, superseded by the radial spokes. That is a retraction of an
  earlier model topology in the same tree, stated by the source about itself.
- **Reachability.** The seam's cable adapters bind device kernels from a shared network module, so
  the device path is unreachable on a CPU-only host. The host-side helpers
  (`keratin_vimentin_backbone_stiffness`, `keratin_vimentin_persistence_length`,
  `backbone_turnover_reference`) are reachable.
- **Live tests, and what they exercise.** `ffn_sim/tests/ac/engine/test_intermediate_filament_rig.py`
  and `ffn_sim/tests/ac/solid/test_intermediate_filament.py` exist. Given zero kernels in the seam
  and a CUDA-gated device path, what they can exercise on a CPU host is the Protocol conformance, the
  stiffness algebra, and the turnover reference — plumbing and host oracles, not the cage under load.
  Aleph ran none of them and draws no conclusion.
- **Where I looked for retractions:** the three source files' own prose and status constants, their
  exported names, and `git diff be0e5876` on each. That tree's planning documents were not read and
  nothing here depends on them.

## 7. Independent oracle or derivation

The oracle is Aleph's own scope discipline, and the ownership argument in §4 is derived rather than
inherited: it follows from the census's "one owner per filament" and "co-location is not a
connection" rules, and it reaches the same one-owner/two-cards conclusion the source reached
independently. Agreement with `ffn_cellsim` is **corroboration recorded in this ledger**, not
evidence, and no control asserts it.

`validate_registry` is additionally checked against a local stand-in object that performs no
validation of its own
(`::test_the_validator_itself_refuses_independently_of_the_contract_type`), so the refusals are
proven to be the validator's and not the dataclass's, with the premise separately asserted by
`::test_the_stand_in_object_really_performs_no_validation_of_its_own`.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_keratin_and_vimentin_are_distinct_identities_under_one_owner` | two `I` populations, one `E` owner, distinct names, distinct material cards, neither owning state, `population_identity` on the owner |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_the_intermediate_filament_owner_refuses_an_unsourced_nonlinear_response` | the strain-stiffening refusal names the material card and the linear tangent limit |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_every_internal_entry_names_a_registered_explicit_parent` | both populations resolve to the registered owner, and no population is an endpoint |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_only_explicit_entries_own_state` | the owner owns state; the populations own none |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_an_internal_entry_with_no_parent` | R6: an orphaned population is refused rather than promoted to a component |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_parent_that_is_not_a_registered_explicit_component` | R8: a population parented to a non-owner is refused — the exact mistake that would make one population the parent of the other |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_an_explicit_owner_that_owns_nothing` | R5: an owner reduced to a label is refused |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_duplicate_names_and_bad_name_shapes` | R1: two populations that collide on a name make ownership unprovable |

## 10. Numerical and precision envelope

No arithmetic is introduced. Working precision, accumulation precision and conditioning are all
not applicable; the module holds strings, tuples of strings, an enum and `None`. Every control
asserts exactly — set equality, enum identity, substring containment — because a contract that
matched its declared scope only to within a tolerance would assert nothing.

The one numerical statement made is a **refusal**, and it is stated as a structural fact rather than
a bound: the mechanical content of a strain-stiffening law is the location of its knee, so no
tolerance on a reported stiffness is meaningful while the knee is assumed. Outside the envelope,
malformed input raises `CensusContractError` or `KeyError`; there is no degraded path and no default
material card.

## 11. Production-backend residency and transfer

Host-side registry data, resident in the Python process, immutable, never transferred to a device,
zero per-step cost. Constructing it must not require a backend or a GPU — it is what a future
implementation is written against, so it has to be readable with no device present. Nothing in this
lane imports a backend module.

## 12. Comments and docstrings to discard

No source text survives; none was taken. Discarded rather than translated:

- Both status strings (`MONOLITHIC_HOOKEAN_REFERENCE_NONLINEAR_TURNOVER_PENDING`,
  `WLC_STRAIN_STIFFENING_UNFOLDING_CARD_REQUIRED`) — that tree's status vocabulary, replaced by the
  ports ledger and the evidence ladder.
- Internal audit-document references and dated diagnosis filenames cited in the cage module's
  docstring. They point into a tree Aleph does not have; the *finding* is restated in §3 of this
  entry in Aleph's own words, with no path.
- Gate names, lane names, and plan-section numbers.
- The literature attributions on the persistence lengths and moduli. Carrying a citation label
  without having read the citation is the failure `require_sourced` exists to prevent, so both the
  values and their labels stay out.
- The provisional-parameter warning block. Aleph makes the same point structurally instead: the
  contract lists the unsourced nonlinear response as an unsupported claim, so the constraint is
  enforced by a test rather than by a comment.

**What replaces them:** the module docstring's own explanation of identity-versus-storage, written
for Aleph, and the ownership derivation in §4 restated in the contract's `approximation` and
`reentry_condition` fields.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED` on 2026-07-30. The three registry entries construct and the controls pass, but the shared state schema this data registers into was not on disk when this entry was written (§14). |
| Reviewer | **Agent-proposed, unratified.** No PI review. The `RE-DERIVE` verdict and the one-owner/two-cards argument are agent judgements. |
| Rollback | Delete `aleph/state/census_frame_fluid_nucleus.py` and `tests/state/test_census_frame_fluid_nucleus.py`. Nothing imports either, so nothing breaks — which is the honest measure of what this establishes. |

## 14. Honest limits

- **The contract seam is unverified.** `aleph/state/schema.py` did not exist when this ran; the
  module falls back to a local dataclass of identical shape and `SCHEMA_SOURCE` reported
  `local_fallback`. "Identical shape" is `ASSUMED` against a specification, not checked against code.
- **Nothing was executed in the reference.** All statements in §2, §3 and §6 come from reading source
  and diffs. The textual claims — 999 lines, zero `@wp.kernel`, four Protocols, the two status
  strings — are verifiable by `grep` and are stated as textual facts. Nothing is claimed about what
  the cage computes under load.
- **`INHERITED_UNVERIFIED`:** none. No value, constant or citation crosses.
- **Believed but not tested:** that one owner with per-population cards will still be the right
  factorisation once a real crosslink topology exists between the two populations. If keratin and
  vimentin turn out to require distinct network topologies rather than distinct cards, the re-entry
  condition on both entries is the correct trigger — and it has never been exercised.
- **Not established:** any intermediate-filament physics. No cable law, no strain-stiffening
  response, no turnover kinetics, no material card, and no oracle for any of them exists in Aleph.
  The radial-spoke negative result of §3 is recorded as **prior art that Aleph has not reproduced**;
  it is an interesting claim from an unrun tree, and it should be re-derived before anything is
  built on it.
- **Needs a literature reading that has not happened:** the keratin and vimentin material cards. Both
  populations sit at `citation_status = UNSOURCED` and must stay there until read.
