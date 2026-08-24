# ALEPH-PORT-1804 — nucleus compartment contract, the nucleoplasm and nuclear-pore scope limits, and the audit of the reference nucleus

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1804` |
| Lane | `L18 census — internal frame / internal fluid / nucleus` |
| Status | `PROPOSED` |
| Port class | `RE-DERIVED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Verdict | **RE-DERIVE.** Zero lines taken. This is the closest call in the lane: the reference nucleus is the best-matched asset this lane examined, and it still loses to clean-room on two specific grounds. |

---

## 1. Aleph API

This entry authorises the `nucleus`, `nuclear_envelope`, `nuclear_lamina`, `chromatin`,
`nucleoplasm` and `nuclear_pores` entries of:

```python
from aleph.state.census_frame_fluid_nucleus import (
    DUAL_SCOPE_TAGS,
    NUCLEUS_COMPONENTS,
    by_name,
    validate_registry,
)
```

No nuclear mechanics is authorised. There is no envelope surface, no lamina law, no chromatin
network, no volume constraint and no drag coefficient in `aleph/`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/nucleus/envelope.py`, `ffn_sim/ac/nucleus/lamina_analytic.py`, `ffn_sim/ac/nucleus/chromatin_analytic.py`, `ffn_sim/ac/nucleus/geometry.py`, `ffn_sim/ac/nucleus/mask_provider.py`, `ffn_sim/ff/nucleus_envelope.py`, `ffn_sim/ac/engine/contracts.py` (for the connector-family vocabulary), `ffn_sim/ac/engine/monomer_flux.py` |
| Source symbol(s) | `NucleusMesh`, `build_nucleus`, `face_reference_areas`, `lamina_areal_tension_kernel`, `rupture_update_kernel`, `conditional_rupture_update_kernel`, `nucleoplasm_volume_reduce_kernel`, `nucleoplasm_volume_force_kernel`, `nucleoplasm_viscosity_kernel`, `linc_tether_kernel`, `chromatin_wlc_kernel`, `LaminaParams`, `lamina_tension`, `lamina_tangent_modulus`, `calibrate_kappa_tilde`, `dihedral_bending_energy`, `dihedral_bending_forces`, `sphere_willmore_energy`, `is_ruptured`, `wlc_tension`, `wlc_energy`, `wlc_small_strain_stiffness`, `mesh_volume`, `mesh_volume_gradient`, `oblate_volume`, `ConnectorFamily`, `MonomerFluxConnector` |
| Read from | **working tree.** |
| Working tree == commit? | **YES** for every path above — verified individually with `git diff be0e5876 -- <path>`; all report no change. |

## 3. Why source-derived porting beats clean-room

**It does not, and this one was close.** The reference nucleus matches the registry's declared
structure more exactly than any other asset in this lane: an envelope mesh with per-face reference
areas, a three-regime areal lamina tension with a rupture criterion, a dihedral bending term, a
worm-like-chain chromatin network, and — precisely as the registry declares for the `H/I` nucleoplasm
— an incompressible volume penalty plus an overdamped viscous drag, and **nothing else**. Two
independent grounds still decide against porting:

1. **Its parameters are the whole model, and they are unsourced.** The lamina response is nonlinear
   with a knee, a plateau and a re-stiffening branch; the mechanical content of the model *is* where
   those transitions sit. The bending stiffness is set by a calibration against a reference surface
   integral rather than being an independent measurement. Porting the code means porting the shape of
   the response, and its shape without a sourced card is a fitted assumption. Aleph would inherit an
   unverified constitutive choice while gaining only a discretisation it can write itself.
2. **The pure-NumPy analytic modules are oracles, and an oracle must be independent of what it
   grades.** `lamina_analytic.py` and `chromatin_analytic.py` are the strongest artefacts here — and
   they are strongest *as oracles*, which means Aleph must not take both them and the kernels they
   grade. Taking only the oracle is the standard oracle-porting error: an oracle ported from the same
   author as the implementation it will grade is not independent, and PLAN §0.2 is explicit that
   matching the source numerically is not evidence. Aleph's own route is already built —
   `validation/analytic/` holds closed-form oracles under a firewall that `aleph/runtime/**` may not
   import — and a Willmore energy for a sphere, a worm-like-chain tension, and a divergence-theorem
   volume gradient are all short independent derivations. `ALEPH-PORT-204` already re-derived discrete
   curvature and found the source's own defect on the way, which is direct evidence about how this
   would go.

**Recorded as prior art, unported:**

- **The nucleoplasm reduction is corroborated, not inherited.** That an incompressible volume
  constraint plus a nodal viscous drag is *sufficient* to stand in for a nuclear fluid domain is a
  modelling claim, and the reference implements exactly those two terms and no third. That is
  independent support for the registry's `H/I` treatment, and it is recorded here as corroboration
  because Aleph reached the same decomposition from the two-effects argument in §4(b).
- **Per-face independent rupture.** The reference tears faces individually rather than through a single
  global tension scalar, and states that the single-global-sigma predecessor was a scaffold. That is a
  real observation about discretisation: a global scalar cannot localise a tear, so it can only report
  that a nucleus ruptured somewhere, which is not a mechanical statement.
- **A mask-provider protocol** so the fluid domain and the nucleus never co-edit each other's state.
  Aleph's one-owner rule reaches the same place independently; see `ALEPH-PORT-1803` §3.

## 4. Physical or mathematical law represented

No law is ported. Three derivations Aleph owns, since the six contract entries are unreadable without
them:

**(a) Why the nucleus is one owner with four internal substructures.** A component is a legal
connector endpoint. Count the connectors that reach the nucleus in the seed census: LINC from the
perinuclear actin cap, LINC from microtubules, LINC from intermediate filaments, non-adhesive
compression from the cortex, and the fluid boundary with the cytosol — five, and **every one of them
lands on the envelope surface**. The lamina is a constitutive law evaluated on that surface's faces
and hinges; the chromatin is a network interior to it; the nucleoplasm is a constraint on the volume
it encloses. None of the four has an endpoint an external connector could address that is not already
the envelope's. Three components would therefore overlap geometrically, share the same endpoints, and
create an ownership question — which of them carries a cortical compression load — with no mechanical
fact to settle it. One owner, four internal substructures, one endpoint set. This is the rule that
keeps a nucleus from becoming three components that overlap, and `validate_registry` enforces it as
R6/R8.

**(b) Why the nucleoplasm is `H/I` and what the encoding costs.** A nuclear fluid domain would
contribute two mechanical effects at the scale the census resolves: a **volumetric** one — the
nucleus resists compression at constant surface area because its contents are nearly incompressible —
and a **rate-limiting** one — its relaxation takes time because material must move through a viscous
medium. Replace the field by an incompressible volume constraint on the enclosed volume and a viscous
drag on the nodes that move through it, and both effects are reproduced with two parameters. Nothing
else is: there is no nuclear pressure field, no nuclear grid, no flow. So the honest tag is composite,
and both halves are load-bearing — it is homogenized (an effective law with a declared inadequacy) and
internal (it has no endpoint of its own).

`ScopeTag` is single-valued, so a decision was required — and **the first decision was wrong, which is
recorded here because being caught is the useful part.** The original encoding was `HOMOGENIZED` with
`parent = "nucleus"`, reasoning that the tag would carry the homogenization so the empty-owned-state
rule would bite, while the parent carried the internality.

The shared state schema refused it. Registering these sixteen entries into a shared manifest produced
exactly one problem — `V-PARENT-ON-NON-INTERNAL` — whose argument is short and correct: a parent *means*
the entity is resolved inside another, which is precisely what `I` says, so a parent on any other tag
asserts two different things at once. The final encoding is therefore:

* `scope = INTERNAL` with `parent = "nucleus"`, so the tag and the parent agree and the entry registers
  cleanly;
* the homogenized half held in `DUAL_SCOPE_TAGS` and enforced by **rule R12** of `validate_registry`,
  which requires no owned state and a non-empty declared `approximation` for any entry whose composite
  tag contains an `H`, whichever single tag it actually carries.

That is a better design than the original and the reason is worth stating: picking `HOMOGENIZED` got the
empty-owned-state check *for free*, and "for free" was the whole appeal. A free check that depends on a
tag meaning something it does not is how a rule ends up enforcing the wrong thing — and had the
homogenized half been merely documented, re-tagging to satisfy the schema would have silently dropped
the obligation. R12 costs six lines, names what it enforces, and has its own negative control.

**(c) Why nuclear pores are `X` and not `H`, and why that is the sharper claim.** A homogenized pore
population would be a small effective envelope permeability standing in for many pores. That is not
what happens: the nucleus-cytosol boundary is **mechanically impermeable**, the cross-boundary
volumetric flux is exactly zero rather than small, and there is no coefficient anywhere that could be
read as a transport rate. `EXCLUDED` therefore says something stronger than `HOMOGENIZED` would: not
"this is approximated coarsely" but "this channel is absent, and any nucleocytoplasmic-transport
conclusion is refused rather than approximated". Nothing is silently absorbed, which is exactly what
the `X` tag exists to guarantee.

The re-entry condition names **two** requirements, and neither substitutes for the other: (1) an
explicit nuclear-pore population, giving pores an identity, a count and a location on the envelope;
and (2) `CHEMICAL_FLUX` connectors between the nucleus and the cytosol, because a pore with no
connector transports nothing. The second is the one that is easy to underestimate. The
`CHEMICAL_FLUX` family exists in the architecture vocabulary and **none of the 36 registered
connectors of the seed census uses it** — the family has zero registered members. So this is not an
extension of an existing transport channel; it is the family's first registered use, and the contract
says so. Adding a permeability to the envelope instead would be homogenizing an excluded structure,
which the entry exists to forbid.

**Verified in the source, with one correction to the manuscript — see §6a.** Nuclear pores are absent
from the reference tree: a case-insensitive search for `nuclear[_ ]?pore|nucleoporin|nucleocytoplasmic`
across all Python returns **zero** hits, and the fluid boundary at the envelope is implemented as a
relative no-flux condition, i.e. exactly impermeable. The registry's `X` tag and the impermeability
claim are both confirmed.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `name`, `group`, `parent` | — (identifier) | — | lower snake_case; `parent` resolves to a registered `E` entry |
| `scope` | — (`ScopeTag`) | — | one of `E I H B X` |
| `owned_state` | — (state-slot names) | — | nine slots on `nucleus`; **empty** on all five others |
| `reentry_condition` | — (prose) | — | non-empty and mandatory for `X` (R9) |

No dimensional quantity is registered. This is deliberate for the reason §3 point 1 gives: the lamina
knee, the plateau, the re-stiffening slope, the volume modulus and the drag coefficient are the model,
and putting any of them in a census entry would inherit a fitted constitutive choice without a source.
Dimensional quantities belong in `aleph/units` under `require_sourced` (`ALEPH-PORT-101`).

Singular and boundary cases:

- Any of the five substructures acquiring owned state: refused (R4) for `nucleoplasm` and
  `nuclear_pores`; for the three `I` entries it is *allowed by the rules* and refused by a control
  (`::test_nucleus_is_one_owner_and_its_four_substructures_are_internal`), because "the envelope's
  geometry lives on the nucleus owner's arrays" is a design fact about this group rather than a
  universal rule — the MTOC in `ALEPH-PORT-1801` is an internal entry that legitimately owns five
  slots.
- `nucleoplasm` gaining owned state: refused by **R12**, not by R4 — its tag is `I`, and R4 only fires
  on `H`, `B` and `X`. This is the asymmetry the composite tag creates, and it is why R12 exists rather
  than being left as prose.
- `nucleoplasm` losing its declared `approximation`: refused by R12 for the same reason.
- `nuclear_pores` with an empty `reentry_condition`: refused (R9). An exclusion with no stated way back
  is indistinguishable from an omission.
- `nuclear_pores` given a parent: not refused by any rule, and asserted `None` by a control — an
  excluded structure is not internal to anything, it is absent.

Invariants, each with the test that asserts it:

- **I1.** `nucleus` is the only `E` entry in the group, owns nine slots, and the other five own
  nothing — `tests/state/test_census_frame_fluid_nucleus.py::test_nucleus_is_one_owner_and_its_four_substructures_are_internal`.
- **I2.** Envelope, lamina, chromatin and nucleoplasm all name `nucleus` as parent, and none is `E` —
  same test.
- **I3.** The envelope is declared the moving impermeable boundary of the cytosol, with the flux stated
  as exactly zero — `::test_the_nuclear_envelope_is_the_moving_impermeable_boundary_of_the_cytosol`.
- **I4.** `nucleoplasm` is `I` with parent `nucleus`, owns nothing, records `H/I` in `DUAL_SCOPE_TAGS`,
  carries a declared approximation, and declares no three-dimensional fluid domain —
  `::test_nucleoplasm_is_dual_tagged_and_carries_no_cfd_domain`. The composite tag's own premise — that
  it names a registered entry — is asserted by
  `::test_the_composite_tag_is_enforced_rather_than_documented`.
- **I4b.** All sixteen entries register into the shared state manifest with **zero** validation
  problems — `::test_the_registry_registers_into_the_shared_state_manifest_without_a_problem`. This is
  the control that caught the original encoding.
- **I5.** `nuclear_pores` is `X`, not `H`; the boundary is mechanically impermeable; the entry says
  explicitly that it is excluded and not homogenized —
  `::test_nuclear_pores_are_excluded_and_the_boundary_is_mechanically_impermeable`.
- **I6.** The re-entry condition names `CHEMICAL_FLUX`, an explicit pore population, the 36 registered
  connectors, zero registered members, and first use of the family —
  `::test_the_nuclear_pore_reentry_condition_names_chemical_flux_connectors`.
- **I7.** The nucleus group contains exactly the six expected names —
  `::test_nucleus_is_one_owner_and_its_four_substructures_are_internal`.

## 6. Source evidence class and known retractions

- **Self-declared execution status.** `ac/nucleus/envelope.py` states that its kernels are net-new,
  are not executed on the authoring machine because it has no CUDA device, and are gated by the
  pure-NumPy oracles in the `*_analytic` modules with the native gates run elsewhere. So the device
  path's evidence class in that tree is "written, oracle-gated, run by someone else" — which is
  honestly labelled and is not "validated".
- **A retraction found, and it is about discretisation.** The single-global-tension predecessor is
  described as a scaffold, superseded by per-face independent rupture. That is a retraction of an
  earlier model in the same tree, stated by the source about itself, and it is the useful kind.
- **A calibration that is not a measurement.** The bending stiffness is set by matching a reference
  surface integral rather than measured independently. Not a defect — it is a legitimate
  normalisation — but it means a bending number taken from that tree would carry the calibration
  convention with it, which is the shape of the 1-ULP and 300 K-convention inheritances the ports
  discipline test already guards against.
- **Reachability.** All eight nucleus kernels are CUDA-gated and unreachable on a CPU-only host. The
  four analytic/geometry modules are pure NumPy and reachable.
- **Live tests, and what they exercise.** Seven test modules under `ffn_sim/tests/ac/nucleus/` name
  chromatin, lamina area, lamina bending, LINC capstan, mask provider, rupture and volume. On a CPU
  host these exercise the analytic oracles and the host geometry — that is the law, not merely the
  plumbing — and they are the reason §3 rates this asset the best-matched in the lane. Aleph ran none
  of them and claims nothing from their existence.
- **Where I looked for retractions:** the four nucleus modules' own prose and status notes,
  `git diff be0e5876` on every path in §2, a tree-wide search for nuclear-pore vocabulary, and the
  connector-family enum together with every use of `CHEMICAL_FLUX`. That tree's planning documents
  were not read.

### 6a. A correction to the source material — the manuscript's chemical-connector limit is stale

Appendix A states, under "CURRENT CHEMICAL-CONNECTOR LIMIT", that the `CHEMICAL_FLUX` family exists in
the architecture vocabulary, that **none of the current 36 registered connectors uses it**, and that
G-actin transport is instead carried through the protrusion-cytosol transfer endpoints (connectors 24
and 25).

The first two statements are true **as statements about the 36-connector registry**, and that is how
the contract registers them. The third is stale with respect to the reference tree. There is a
G-actin conserved-pool connector runtime there whose declared family is
`ConnectorFamily.CHEMICAL_FLUX`, with a host-side acceptance oracle and structural tests, and it is a
cytosol-facing connector that is **not one of the 36** in Appendix A. So in the reference the family
does have a member, obtained by adding a connector outside the registry rather than by re-typing one
of the 36.

Three consequences, recorded rather than acted on:

1. The contract's wording is scoped precisely to *the 36 registered connectors*, which remains true,
   and does not claim the family is unused everywhere. The distinction matters: "zero registered
   members" is a fact about the registry, and it is the fact the re-entry condition needs.
2. Appendix A's own header warns that the registry is the 29 July seed scope and must be reconciled
   with the current composition contract before citation. This is a concrete instance of that warning
   firing, found by reading code rather than prose.
3. It slightly *weakens* the "first use of the family" framing — someone has since built one, elsewhere
   — while *strengthening* the underlying point: the family's first use arrived as an off-registry
   connector, which is exactly how a census stops describing the thing it is a census of. Aleph's
   answer to that is `ALEPH-PD-002` — connectors are registered data, so adding one extends a manifest
   and is visible.

## 7. Independent oracle or derivation

Nothing here is checked against `ffn_cellsim`:

1. **The one-owner argument** in §4(a) is derived by counting the connectors that reach the nucleus in
   the seed census and observing that all five land on the same surface. It settles the ownership
   question with no implementation in view.
2. **The two-effects argument** in §4(b) derives the `H/I` reduction from what a nuclear fluid
   contributes at the resolved scale, and predicts that exactly two terms suffice. The reference
   implements exactly two terms and no third — which is a **prediction confirmed by inspection**, and
   is recorded as corroboration in §3 rather than used as a control.
3. **The `X`-versus-`H` argument** in §4(c) turns on exact zero versus small, which is a structural
   distinction rather than a numerical one and needs no reference at all.
4. `validate_registry` is checked against a local stand-in object performing no validation of its own
   (`::test_the_validator_itself_refuses_independently_of_the_contract_type`), premise asserted by
   `::test_the_stand_in_object_really_performs_no_validation_of_its_own`. So R9 — an exclusion must
   state a way back — is proven to be enforced by the validator and not by the dataclass.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_nucleus_is_one_owner_and_its_four_substructures_are_internal` | one `E` owner, four substructures parented to it, none `E`, none owning state, and exactly six names in the group |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_the_nuclear_envelope_is_the_moving_impermeable_boundary_of_the_cytosol` | the moving impermeable boundary role is stated on the envelope, with exactly-zero flux |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_nucleoplasm_is_dual_tagged_and_carries_no_cfd_domain` | `I` + parent `nucleus` + `DUAL_SCOPE_TAGS["nucleoplasm"] == "H/I"` + declared approximation, volume constraint and viscous drag, no fluid domain |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_the_registry_registers_into_the_shared_state_manifest_without_a_problem` | all sixteen entries register into the shared state manifest with zero validation problems |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_nuclear_pores_are_excluded_and_the_boundary_is_mechanically_impermeable` | `X` not `H`; mechanically impermeable; not folded into a permeability |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_the_nuclear_pore_reentry_condition_names_chemical_flux_connectors` | `CHEMICAL_FLUX`, explicit pore population, 36 registered connectors, zero registered members, first use |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_an_excluded_entry_with_no_reentry_condition` | R9: emptying the nuclear-pore re-entry condition is rejected — an exclusion with no way back is an omission |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_an_excluded_entry_that_owns_state` | R4 on `X`: giving nuclear pores a permeability array is rejected, so the exclusion cannot become a homogenization by the back door |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_an_internal_entry_with_no_parent` | R6: an orphaned nuclear lamina is rejected rather than promoted to a component |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_parent_that_is_not_a_registered_explicit_component` | R8: chromatin parented to the envelope — a substructure of a substructure — is rejected |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_self_parented_entry` | R11: a nucleus that is its own parent, the one-cycle no other rule catches |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_composite_h_tagged_entry_that_owns_state` | R12: giving the nucleoplasm a pressure field is rejected even though its tag is `I` — the control that proves the composite tag is a rule and not a comment |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_composite_h_tagged_entry_with_no_approximation` | R12's other half: the declared inadequacy is required despite the `I` tag |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_the_dual_scope_tag_presence_check_refuses_a_tag_with_no_entry` | a composite tag naming a component that does not exist is rejected, so the rule cannot end up guarding nothing |

### Mutation evidence that R12 is not decoration

R12 was deliberately disabled in `validate_registry` and the suite re-run. **Exactly two controls
failed** — `test_validation_refuses_a_composite_h_tagged_entry_that_owns_state` and
`test_validation_refuses_a_composite_h_tagged_entry_with_no_approximation` — and nothing else. The rule
was then restored and the suite verified green.

Two failures and no more is the right sensitivity: the rule is load-bearing for exactly the two claims
that name it, and no other control was silently depending on it. Note what the first failure means in
detail — with R12 off, giving the nucleoplasm a nuclear pressure field passes validation, because R4
does not reach an `I` entry. That is the encoding hazard of §4(b) made concrete: had the homogenized
half of `H/I` been documented rather than enforced, the re-tag from `H` to `I` would have silently
dropped the obligation and no test would have noticed.

## 10. Numerical and precision envelope

No arithmetic is introduced. Working precision, accumulation precision and conditioning are not
applicable; the controls assert exactly — set equality, enum identity, substring containment — because
a contract that matched its declared scope only to within a tolerance would assert nothing.

Two numerically substantive statements are made, both as structural facts rather than bounds:

- **The envelope flux is exactly zero, not small.** This is the reason the pores are `X` rather than
  `H`, and it is a statement no tolerance applies to: an impermeable boundary contributes an identically
  zero cross-boundary term, on the whole branch, not a term below a threshold. Aleph has met this
  distinction before and it earned its keep — `ALEPH-PORT-1103` records a tensile-only tether whose
  inactive branch returns exact `0.0`, testable with `==` precisely because the derivative is
  identically zero on the whole branch rather than merely small near the gate.
- **No error bar is meaningful on the lamina or nucleoplasm parameters while their cards are
  unsourced**, for the reason in §3 point 1: the shape of the nonlinearity *is* the model, so an
  interval around an assumed knee would imply the central value estimates something.

Outside the envelope, malformed input raises `CensusContractError` or `KeyError`. There is no degraded
path, no default parameter, and no fallback entry.

## 11. Production-backend residency and transfer

Host-side registry data, resident in the Python process, immutable, never transferred to a device, zero
per-step cost. A future Aleph nucleus **is** a device-resident object with a real transfer cost, and
this contract is what it will be written against, so the contract must be constructible with no backend
imported and no device present. Nothing in this lane imports a backend; the census module imports
`dataclasses` and, conditionally, the shared state schema, and nothing else.

## 12. Comments and docstrings to discard

No source text survives; none was taken. Discarded rather than translated:

- Machine names, device names, and statements about which host runs which gate. They are working-tree
  and hardware assumptions about a setup Aleph does not share, and PLAN §0.1 makes Aleph's GPU
  situation a policy question rather than a comment.
- Plan-section references, framework-numbered feature labels, gate identifiers and slice names.
- Every literature attribution on the lamina regimes, the rupture strain, the volume modulus, the
  nuclear drag coefficient and the chromatin persistence length. None of these citations was read by
  this lane, so neither the values nor their labels cross. Carrying a citation one has not read is the
  failure `require_sourced` exists to prevent.
- The calibration convention on the bending stiffness. If Aleph ever needs a nuclear bending stiffness
  it must derive its own normalisation and state it, because an inherited convention is exactly how a
  systematic offset becomes invisible.
- Cross-references to the membrane bending kernel that the nucleus module re-exports. Aleph's membrane
  bending is `ALEPH-PORT-1101`, re-derived, and the two must not be linked by inherited prose.

**What replaces them:** the module docstring's own account of the two composite-tag encodings and why
`parent` is not exclusive to `INTERNAL`, plus the three derivations of §4 restated in the
`mechanical_role`, `approximation`, `unsupported_claims` and `reentry_condition` fields in Aleph's
vocabulary.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED` on 2026-07-30. The six registry entries construct and the controls pass, but the shared state schema was not on disk when this entry was written (§14). |
| Reviewer | **Agent-proposed, unratified.** No PI review. The `RE-DERIVE` verdict, the `H/I` encoding decision, and the scoping of the chemical-connector statement in §6a are agent judgements. |
| Rollback | Delete `aleph/state/census_frame_fluid_nucleus.py` and `tests/state/test_census_frame_fluid_nucleus.py`. Nothing imports either, so nothing breaks. |

## 14. Honest limits

- **The contract seam is unverified.** `aleph/state/schema.py` did not exist when this ran; the module
  falls back to a local dataclass of identical shape and `SCHEMA_SOURCE` reported `local_fallback`.
  "Identical shape" is `ASSUMED` against a specification, not checked against code.
- **The `H/I` encoding is a lane decision that may not survive, and its first form did not.** If the
  shared schema later admits a composite tag, or a set of tags, `nucleoplasm` should move to it and
  `DUAL_SCOPE_TAGS` plus rule R12 should be deleted. The current form has a real cost: a reader who
  greps for `ScopeTag.HOMOGENIZED` will not find the nucleoplasm, and only R12 and its two negative
  controls stop that from mattering. Flagged for the PI and the schema lane rather than settled.
- **Nothing was executed in the reference.** Every statement in §2, §3 and §6 comes from reading source,
  diffs and text searches. The nuclear-pore absence is a search result: it establishes that the
  vocabulary is absent, which for a pore population is strong, but a search is not a proof that no
  equivalent exists under another name.
- **`INHERITED_UNVERIFIED`:** none. No lamina regime, rupture strain, volume modulus, drag coefficient,
  persistence length, bending stiffness or citation crosses.
- **Believed but not tested:** that nine state slots are the right decomposition for a nucleus owner,
  and that `reduced_coordinates` belongs on the owner rather than on a separate reduction object. Both
  are hypotheses about an implementation that does not exist.
- **The modal-reduction validation obligation is registered and unenforced.** Both `nucleus` and
  `chromatin` refuse results from an unvalidated reduction, and there is no mechanism in Aleph that
  could detect one — no reduction, no native representation, nothing to compare. It is a contract term
  awaiting a runtime, and it should be read as a requirement on future work rather than as a guarantee.
- **§6a is a finding about the source material, not a fix.** Appendix A's chemical-connector limit is
  stale and this entry scopes around it. Nobody has reconciled the seed registry with the reference
  tree's current composition, and Appendix A's own header asks for exactly that before citation. That
  reconciliation is not this lane's and has not happened.
- **Needs a literature reading that has not happened:** every nuclear material parameter. All six
  entries sit at `citation_status = UNSOURCED` and must stay there until read.
- **Not established:** any nuclear physics. No envelope surface, no lamina law, no chromatin network, no
  volume constraint, no drag, and no oracle for any of them exists in Aleph.
