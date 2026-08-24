# ALEPH-PORT-1803 — cytosol compartment contract, the organelle homogenization, and the audit of the reference fluid

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1803` |
| Lane | `L18 census — internal frame / internal fluid / nucleus` |
| Status | `PROPOSED` |
| Port class | `RE-DERIVED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Verdict | **RE-DERIVE.** Zero lines taken. The physics is textbook Terzaghi consolidation; the non-obvious part is a moving-boundary bookkeeping rule, and Aleph must derive that against its own transaction contract or the guarantee does not hold. |

---

## 1. Aleph API

This entry authorises the `cytosol` entry and the four homogenized organelle entries of:

```python
from aleph.state.census_frame_fluid_nucleus import (
    INTERNAL_FLUID_COMPONENTS,
    ORGANELLE_HOMOGENIZATION_NOTE,
    ORGANELLE_REENTRY_CONDITION,
    ORGANELLE_UNSUPPORTED_CLAIMS,
    UNDOCUMENTED_HOMOGENIZATION_MAP,
    UNDOCUMENTED_HOMOGENIZATION_MAP_ENTRIES,
    CITATION_STATUS_RATIONALE,
    by_name,
    validate_registry,
)
```

No fluid solver is authorised. There is no pressure field, no Darcy flux, no control-volume
discretisation and no permeability value in `aleph/`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/fluid_core.py`, `ffn_sim/ff/biot_fluid_warp.py`, `ffn_sim/ac/fluid/domain.py`, `ffn_sim/ac/fluid/velocity.py`, `ffn_sim/ac/fluid/transport.py`, `ffn_sim/ac/fluid/darcy_analytic.py`, `ffn_sim/ac/fluid/transport_reference.py`, `ffn_sim/ac/cell/ng4_payload.py` |
| Source symbol(s) | `FluidCore`, `FluidVolumeStateOwner`, `ReducedCoreStateOwner`, `MovingImpermeableFluidBoundaryFacade`, `NucleusPressureAdjointBoundary`, `BiotSubstrateFluidSolver`, `lumped_surface_control_volumes`, `biot_diffusion_kernel`, `Domain`, `NucleusMaskProvider`, `StaticSphereNucleusMaskProvider`, `PHI_PROVISIONAL` |
| Read from | **working tree.** |
| Working tree == commit? | **YES** for every path above — verified individually with `git diff be0e5876 -- <path>`; all report no change. |

## 3. Why source-derived porting beats clean-room

**It does not. This is `RE-DERIVED`**, and the reasoning is different from the other three entries in
this lane, because here the source is *good* and the answer is still no.

1. **`fluid_core.py` is a re-composition, and it says so.** An earlier audit rated it substantial;
   verified, and the rating is about craft rather than content. Its own opening states that it does
   not introduce another fluid or nuclear constitutive law and instead injects existing solvers behind
   state-owner, moving-boundary, mechanics, transaction and ledger interfaces. It is 1030 lines, one
   `@wp.kernel`, five Protocols. So what is on offer is an **ownership seam**, not a fluid — and it is
   a seam over the connector-contract system `ALEPH-PD-002` deliberately replaces with registered
   data. Porting it would re-import the structure the decision removed.
2. **The law underneath is textbook.** The field solver integrates Biot/Terzaghi consolidation,
   `dp/dt = c_v grad^2 p`, on a device grid with an explicit stencil, `c_v = k M / mu`. That is a
   linear diffusion equation with a named coefficient, attributed in the source to Darcy 1856, Biot
   1941 and a standard porous-media text. Per the template's own rule, a form any competent author
   arrives at independently is prior art and not a reason to port. Aleph additionally has the better
   route: `validation/analytic/` already holds closed-form oracles and the L4 lane's machinery, and a
   1-D Terzaghi series solution is exactly the kind of thing that belongs there rather than in a
   ported kernel.
3. **The one genuinely non-obvious asset must be re-derived to be worth anything.** The conservative
   remap under a moving domain boundary is the real intellectual content: when a control volume
   changes class between interior fluid and exterior as the membrane or nuclear envelope moves, its
   content must be *transported* across the moving face rather than gained or lost, so the total
   changes only by genuine boundary flux plus that transport. That guarantee is a statement about
   Aleph's accepted-step transaction — which volumes are reclassified when, and what is committed
   together — and a ported implementation would be guaranteeing something about a transaction it was
   not written for. The guarantee does not survive the port even if the arithmetic does. Aleph has to
   derive it against `ALEPH-PORT-305`.
4. **A defect that alone would block a port.** 48 modules in that tree call `wp.init()` at module
   scope, including `ff/biot_fluid_warp.py`, the file that holds the pressure-field kernels. Importing
   it initialises a GPU runtime as an **import side effect**. Under PLAN §0.1 and
   `docs/design/GPU_POLICY.md` an Aleph module that did that would violate the zero-GPU policy simply
   by being imported, and `scripts/gpu_preflight.py` — which refuses by default and exits 2 — would be
   bypassed by an `import`. Any port of that file would have to strip the call, at which point the
   remaining content is the diffusion stencil of point 2.

**Recorded as prior art, unported:** (i) the moving-boundary content-transport identity of point 3,
as a **specification for an oracle Aleph should write**, not as code; (ii) the observation that a
lumped nodal drag and a resolved Darcy drag double-count the same viscosity if both are active — the
source names this as its single biggest correctness hazard and it is worth writing down, because it is
a trap Aleph could walk into independently the moment it has both a filament drag and a fluid field;
(iii) treating the nuclear-envelope mask as a *provider protocol* so the fluid domain and the nucleus
never co-edit each other's state, which is the same conclusion Aleph's one-owner rule reaches from the
other direction.

## 4. Physical or mathematical law represented

No law is ported. Two derivations Aleph owns, because the contract is unreadable without them:

**(a) Why the membrane and nuclear-envelope couplings are boundary connectors and not contacts.**
A contact transmits a normal force at the points where two bodies touch; its endpoint set is
determined by proximity and it carries no information when the bodies are apart. A moving fluid
boundary is different in kind: the surface imposes a **kinematic condition** on the field over its
whole extent — the normal fluid velocity matches the surface's normal velocity, and for an
impermeable surface the relative flux is exactly zero — and returns a **pressure traction**
distributed over every face, whether or not anything is "touching". So the two couplings differ in
their endpoint sets, in what they transfer, and in whether they can be inactive. Registering both as
contacts would give the cytosol a load path that switches off when the geometry separates, which for
an enclosed fluid is not a limiting case but a contradiction. The contract records the distinction on
`cytosol.mechanical_role`, and it is the reason `domain_cell_classification` is **owned state**: which
control volumes are interior fluid is a function of where the moving boundaries are, so it changes
during the run and must be committed with the step.

**(b) Why the organelle homogenization is a documented gap and not a modelling choice.**
Four organelle families — endoplasmic reticulum, mitochondria, Golgi apparatus, and the
endolysosomal/vesicle populations — occupy a large fraction of the cytoplasmic volume and are
absorbed into three effective coefficients of the cytosol field: porosity, permeability, viscosity.
Absorption is defensible. What makes it a *gap* is that the map is missing: there is no written
relation `organelle inventory -> (phi, k, mu)`. Without that map the three coefficients are not
underdetermined-but-meaningful, they are **uninterpretable**: no organelle census constrains them, and
no measured value of them constrains an organelle census. They are calibration parameters of a lumped
field. So a reader who sees a cytosol permeability and no note beside it will read it as a measured
property of cytoplasm, and be wrong in a way nothing in the model contradicts.

That distinction needed a label, and here is where the lane hit a real cross-lane constraint rather
than a preference. The shared state schema restricts `citation_status` to a closed vocabulary —
`UNSOURCED`, `PROPOSED`, `SOURCED`, `INHERITED_UNVERIFIED`, `NOT_APPLICABLE` — and refuses anything
else at construction, deliberately, so that a typo cannot invent a status that looks stronger than it
is. That is a good rule and this lane does not get to edit around it. **The schema lane's module was
not on disk when this registry data was written and appeared partway through, at which point the
attempted `UNDOCUMENTED_HOMOGENIZATION_MAP` status was refused at construction.**

Of the five allowed values, `UNSOURCED` is the only honest one: the coefficients genuinely have no
source. `NOT_APPLICABLE` would be false and `INHERITED_UNVERIFIED` would imply something was carried
in from elsewhere, which nothing here was. So all five entries carry `UNSOURCED`, and the stronger
fact travels beside it in two machine-readable places that a test asserts:
`UNDOCUMENTED_HOMOGENIZATION_MAP_ENTRIES`, the set of entries whose status cannot be promoted by a
literature search alone, and `CITATION_STATUS_RATIONALE`, which names the reason per entry.
**The cytosol is in that set**, because a component that absorbs four undocumented approximations and
then reports a clean citation status has merely relocated the invisible claim.

This is a recorded request to the schema lane, in §14, and not a workaround: the right fix is for the
shared vocabulary to gain a value for "no interpretation exists yet".

**Verified independently in the source, and this is the strongest confirmation in the lane.**
A case-insensitive search of the whole reference Python tree for `organelle|crowding|mitochondri|golgi|endoplasmic|endolysosom`
returns **16 hits, none of which is an organelle representation** — they are a neighbour-count in an
unrelated aggregate model, a bundle-detection density statistic, and comments in steric-parameter
derivations. There is no ER, no mitochondrion, no Golgi and no vesicle population anywhere in that
tree. Meanwhile porosity is repeatedly annotated as an unfilled gap requiring PI authorship, and one
assembly payload carries a literal provisional porosity value flagged in an inline comment as
unsourced. So the manuscript's "missing formal contract" is not a caution about future work: the map
is absent, the coefficient it would define is a placeholder, and the four organelle families that
justify the lumping do not appear in the code at all. Registering that is the single most useful thing
this lane does.

The six claims the approximation does not support are registered verbatim on **each** of the four
entries — organelle-specific deformation, organelle-specific force transmission, spatial mitochondrial
or ER stress, vesicle transport trajectories, vesicle-mediated membrane-area delivery, and spatially
resolved perinuclear crowding. Repeated per entry rather than stated once in a note, because a claim is
made about a specific organelle and the refusal has to be findable from that organelle's own contract.
The four share one text by construction (a single module-level constant), so they cannot drift apart.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `name`, `group`, `parent` | — (identifier) | — | lower snake_case; `parent` resolves to a registered `E` entry |
| `scope` | — (`ScopeTag`) | — | one of `E I H B X` |
| `owned_state` | — (state-slot names) | — | non-empty for `E`; **exactly empty** for all four `H` entries |
| `citation_status` | — (label) | — | `UNSOURCED` on all five, each with an entry in `CITATION_STATUS_RATIONALE` |

No dimensional quantity is registered. This is most important in this entry: the three coefficients
under discussion — porosity, permeability, viscosity — are exactly the numbers whose meaning is
missing, so registering a value for any of them would be the defect. The contract names them and
refuses them. `phi`, `k` and `mu` do not appear in `aleph/` with a value.

Singular and boundary cases:

- Any of the four organelle entries acquiring owned state: refused (R4). This is the specific defect
  where a lumped structure becomes a state owner without being declared one, and it is the reason R4
  exists.
- An organelle entry with an empty `approximation`: refused (R10). A homogenization with no declared
  inadequacy is the invisible claim.
- `cytosol` with empty `owned_state`: refused (R5).
- An organelle entry given a `parent`: refused by the shared schema, whose `V-PARENT-ON-NON-INTERNAL`
  rule holds that a parent means "resolved inside another", which is what `I` says. That is correct
  here: these four have `parent = None` because they are lumped into a *sibling* component and are not
  internal to one. It is also what forced the nucleoplasm encoding in `ALEPH-PORT-1804` §4(b) to
  change.

Invariants, each with the test that asserts it:

- **I1.** All four organelle entries are `H`, own nothing, and have no parent —
  `tests/state/test_census_frame_fluid_nucleus.py::test_all_four_organelle_families_are_homogenized_and_own_no_state`.
- **I2.** All four carry the shared homogenization note verbatim, name all three coefficients, contain
  the phrase that they cannot be interpreted biologically, and carry
  the shared homogenization note verbatim —
  `::test_all_four_organelle_families_carry_the_undocumented_homogenization_map_note`; and each
  appears in `CITATION_STATUS_RATIONALE` at `UNSOURCED` with the undocumented map named as the reason
  — `::test_the_undocumented_map_is_recorded_as_a_citation_rationale_not_swallowed`.
- **I3.** All six unsupported claims appear on all four entries, identical by construction —
  `::test_the_organelle_entries_refuse_all_six_named_claims`.
- **I4.** No organelle family is a legal connector endpoint —
  `::test_no_organelle_family_appears_as_a_connector_endpoint_candidate`.
- **I5.** `cytosol` is Biot/Darcy, names both moving boundaries, states that the couplings are
  boundary connectors and not contacts, and owns the domain classification —
  `::test_cytosol_is_a_porous_fluid_field_with_two_moving_boundaries`.
- **I6.** `cytosol` inherits the undocumented-map citation status and calls its three coefficients
  calibration parameters — `::test_cytosol_inherits_the_undocumented_map_as_its_own_limitation`.

## 6. Source evidence class and known retractions

- **Self-labelled scope limit.** `fluid_core.py` states in its own opening that it is a structural
  optimisation slice and not the authoritative first native-population run, and that production use
  additionally requires response/force/work gates at the physiological operating point. That is the
  source's own evidence class for the seam, and it is explicitly below "validated".
- **Self-labelled gap, repeatedly.** Porosity is annotated as an unfilled gap requiring PI authorship
  in at least five separate fluid modules, and the assembly payload's provisional value carries an
  inline unsourced marker. The source is not hiding this; the manuscript's "missing formal contract"
  is the same fact stated at registry level.
- **A defect not labelled by the source:** the 48 module-scope `wp.init()` calls of §3 point 4. The
  source does not flag these as a problem — under its own policy Warp is the mandated runtime, so an
  init at import is unremarkable there. It is a defect *relative to Aleph's* GPU policy, which is the
  frame that matters for a port, and it is recorded here rather than in `aleph/` because it is a fact
  about the source.
- **Reachability.** The device solvers are CUDA-gated and unreachable on a CPU-only host. The
  host-side analytic references (Darcy, transport, manufactured solutions, the moving-interval
  identity, the finite-volume reference) are reachable and are pure NumPy. Those are the strongest
  artefacts in the fluid area, and notably they are **oracles rather than solvers** — which is
  consistent with Aleph's own conclusion that the oracle belongs in `validation/analytic/`.
- **Live tests, and what they exercise.** Nine test modules under `ffn_sim/tests/ac/fluid/` name
  oracles: consolidation, Green's function, Darcy, transport, manufactured solutions, osmotic
  response, coupled boundary. On a CPU host these exercise the analytic references, which is the law
  and not merely the plumbing — this is the one area in the lane where the source's tests appear to
  test physics. Aleph ran none of them and claims nothing from their existence.
- **Where I looked for retractions:** the seam and solver docstrings, the porosity annotations across
  the fluid package, `git diff be0e5876` on all eight paths, and a tree-wide search for organelle
  vocabulary. That tree's planning documents were not read.

## 7. Independent oracle or derivation

Nothing here is checked against `ffn_cellsim`, and one claim is checked *against* it:

1. **The absent-organelle count is an independent measurement of the source**, not a reading of its
   prose: a tree-wide case-insensitive search returns no organelle representation. It confirms the
   manuscript's "missing formal contract" by a method the manuscript did not use, and it is stronger
   evidence than the manuscript's own sentence.
2. **The boundary-versus-contact derivation** in §4(a) is an argument from the definition of a
   kinematic boundary condition. It settles a registry question — why connectors 04 and 05 of the seed
   census are `FLUID_BOUNDARY` and not `CONTACT` — with no implementation in view.
3. `validate_registry` is checked against a local stand-in object performing no validation of its own
   (`::test_the_validator_itself_refuses_independently_of_the_contract_type`), with its premise
   asserted by `::test_the_stand_in_object_really_performs_no_validation_of_its_own`. So R4 — the rule
   that keeps a lumped organelle from owning state — is proven to be enforced by the validator.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_all_four_organelle_families_are_homogenized_and_own_no_state` | four `H` entries, `owned_state == ()` on each, no parent, all in `internal_fluid` |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_all_four_organelle_families_carry_the_undocumented_homogenization_map_note` | the shared note verbatim, all three coefficients named, and the distinct citation status on each |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_the_organelle_entries_refuse_all_six_named_claims` | all six refusals present on all four entries, identical by construction |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_cytosol_is_a_porous_fluid_field_with_two_moving_boundaries` | Biot/Darcy, both moving boundaries, boundary-connector-not-contact, and the owned domain classification |
| Positive | `tests/state/test_census_frame_fluid_nucleus.py::test_cytosol_inherits_the_undocumented_map_as_its_own_limitation` | the consequence lands on the cytosol too, rather than being relocated |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_homogenized_entry_that_owns_state` | R4, the load-bearing refusal of this lane: giving mitochondria one state slot is rejected, so a lumped organelle cannot become an undeclared owner |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_a_homogenized_entry_with_no_declared_approximation` | R10: emptying the Golgi entry's approximation is rejected — the invisible-claim rule, stated mechanically |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_refuses_an_explicit_owner_that_owns_nothing` | R5: a cytosol that owns no field is rejected |
| Negative (must fail) | `tests/state/test_census_frame_fluid_nucleus.py::test_validation_reports_every_violation_rather_than_only_the_first` | three simultaneous violations are all reported, so one bad edit does not become several fix rounds |

### Mutation evidence that these controls are not decoration

R4 was deliberately disabled in `validate_registry` (the scope check short-circuited to `False`) and
the suite re-run. **Five controls failed and no others**:
`test_validation_refuses_a_homogenized_entry_that_owns_state`,
`test_validation_refuses_an_excluded_entry_that_owns_state`,
`test_validation_reports_every_violation_rather_than_only_the_first`,
`test_the_validator_itself_refuses_independently_of_the_contract_type`, and
`test_the_stand_in_object_really_performs_no_validation_of_its_own`. The rule was then restored and the
file verified byte-for-byte against its pre-mutation copy, with the suite green again.

That is the check worth running on a negative control: a control that keeps passing when the rule it
guards is removed is an assertion about nothing. The last two failures matter most — the
type-independent controls fail too, which shows R4 is enforced by the validator rather than by the
shared schema's own `__post_init__`.

## 10. Numerical and precision envelope

No arithmetic is introduced, so working precision, accumulation precision and conditioning are not
applicable and the controls assert exactly — set equality, enum identity, substring containment. A
contract that matched its declared scope to within a tolerance would assert nothing.

The numerically substantive statement in this entry is a **refusal of a tolerance**, and it belongs
here rather than in §4: while the organelle homogenization map is undocumented, no error bar on a
cytosol porosity, permeability or viscosity is meaningful, because an error bar is a statement about
a quantity whose meaning is fixed. Reporting `phi = 0.5 +/- 0.05` would be more misleading than
reporting nothing — the interval implies the central value is an estimate of something. The contract
therefore records those coefficients as calibration parameters, with no value and no interval, and
`require_sourced` in `aleph/units` is the mechanism that keeps it that way.

Outside the envelope, malformed input raises `CensusContractError` or `KeyError`. There is no
degraded path and no default coefficient.

## 11. Production-backend residency and transfer

Host-side registry data, resident in the Python process, immutable, never transferred to a device,
zero per-step cost.

The residency answer has teeth in this entry specifically. A future Aleph cytosol field **is** a
device-resident object with a real transfer cost, and this contract is what it will be written
against — so the contract must be constructible with no backend imported and no device present.
That requirement is not hypothetical: as recorded in §3 point 4, the reference's pressure-field module
initialises a GPU runtime at import, which would let an `import` bypass a preflight that otherwise
refuses by default and exits 2. Nothing in this lane imports a backend module, and the census module
imports nothing beyond `dataclasses` and, conditionally, the shared state schema.

## 12. Comments and docstrings to discard

No source text survives; none was taken. Discarded rather than translated:

- Every internal gate identifier, plan-section reference, slice label and dated-document citation in
  the fluid modules' docstrings. They point into a tree Aleph does not have.
- The provisional porosity value and its inline unsourced marker. **The value is deliberately not
  recorded anywhere in `aleph/`** — carrying a placeholder number across, even labelled, would be an
  unrecorded inheritance of exactly the kind PLAN §0.2 forbids, and it is precisely the number this
  entry argues has no meaning yet.
- The literature attributions on the consolidation coefficient and the cytoplasmic diffusivity. The
  law is re-derived in §3 point 2 as a linear diffusion equation; the numbers and their citations stay
  out until read.
- The runtime-mandate statements and backend policy claims. `ALEPH-DQ-107` chose Aleph's backend on
  Aleph's own evidence and is a `PROPOSED` decision awaiting the PI.
- The naming of a specific cell line in the viscosity provenance. That datum is already on the ports
  discipline test's provider-datum list and must not appear in `aleph/` without an entry; it does not
  appear at all.

**What replaces them:** the module docstring's own account of why this group carries the census's most
consequential homogenizations, and the derivations of §4 restated in the `approximation`,
`mechanical_role` and `unsupported_claims` fields in Aleph's vocabulary.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED` on 2026-07-30. The five registry entries construct and the controls pass, but the shared state schema was not on disk when this entry was written (§14). |
| Reviewer | **Agent-proposed, unratified.** No PI review. The `RE-DERIVE` verdict, the boundary-versus-contact derivation, and the decision to give the cytosol the same citation status as the organelles are agent judgements. |
| Rollback | Delete `aleph/state/census_frame_fluid_nucleus.py` and `tests/state/test_census_frame_fluid_nucleus.py`. Nothing imports either, so nothing breaks. |

## 14. Honest limits

- **The contract seam is unverified.** `aleph/state/schema.py` did not exist when this ran; the module
  falls back to a local dataclass of identical shape and `SCHEMA_SOURCE` reported `local_fallback`.
  "Identical shape" is `ASSUMED` against a specification, not checked against code.
- **Nothing was executed in the reference.** All statements in §2, §3 and §6 come from reading source,
  diffs, and text searches. The organelle-absence result is a search over the tree and is stated as
  such; it establishes that the vocabulary is absent, which for a mesh, a membrane or a particle
  population is conclusive, but a search is not a proof that no equivalent exists under another name.
- **`INHERITED_UNVERIFIED`:** none. No value, coefficient or citation crosses. In particular no
  porosity, permeability, viscosity or consolidation coefficient appears in `aleph/`.
- **The moving-boundary conservation guarantee is unwritten.** §3 point 3 identifies it as the one
  asset worth having and explicitly does **not** deliver it: there is no Aleph derivation, no oracle,
  and no test. It is recorded as a specification for future work, and until it exists Aleph has no
  moving-boundary conservation property at all.
- **The double-counted-drag trap is recorded, not guarded.** Nothing in Aleph currently prevents a
  future lumped filament drag and a future resolved Darcy drag from both being active. There is no
  test, because there is nothing yet to test. It is on this ledger so that whoever builds the second
  of the two meets it.
- **Believed but not tested:** that seven state slots are the right decomposition for the cytosol
  field, and that `domain_cell_classification` is the correct place to put the moving-boundary
  bookkeeping. Both are hypotheses about an implementation that does not exist.
- **Open cross-lane request, not settled here:** `UNDOCUMENTED_HOMOGENIZATION_MAP` should be a
  first-class member of the shared `citation_status` vocabulary. It is not, so all five entries carry
  `UNSOURCED` — which on its own understates the problem — with the reason held in
  `CITATION_STATUS_RATIONALE` beside it. That is the best available encoding and it is not the right
  one: a reader who checks the status and not the rationale gets the understated version, and nothing
  stops them. It is arguably a *quotability* statement rather than a citation one, the same two-axis
  distinction `ALEPH-PORT-601` draws, which is a further reason it needs a decision rather than a
  lane-local constant. Flagged for the PI and the schema lane; not settled by this entry.
- **A second lane hit the same wall independently.** A sibling census module was refused at
  construction for `citation_status = 'INFERENCE_TARGET'` on a membrane-reservoir entry. Two lanes
  needing a value the closed vocabulary lacks is evidence about the vocabulary rather than about the
  lanes, and it is recorded here because a single instance would have been easy to dismiss.
- **Not established:** any fluid physics. No pressure field, no Darcy law, no consolidation solution,
  no permeability, and no oracle for any of them exists in Aleph.
