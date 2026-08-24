# ALEPH-PORT-1901 — Surface, pressure and fluid-boundary connector contracts (8)

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1901` |
| Lane | `L19 connector contracts — surface / pressure / fluid boundary` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `DEFECT_STUDY` — no code and no prose crossed. One structural defect crossed as an argument, and it was verified rather than inherited. |

---

## 1. Aleph API

```python
from aleph.state.connectors_surface_traction import (
    CONTRACT_FIELD_DEFAULTS,
    CONTRACT_FIELD_ORDER,
    ConnectorContract,
    ConnectorContractError,
    ConnectorFamily,
    MEMBRANE_CORTEX_CONTACT,
    MEMBRANE_CYTOSOL_BOUNDARY,
    MEMBRANE_ECM_CONTACT,
    MEMBRANE_ERM_CORTEX,
    MEMBRANE_MEDIUM_TRACTION,
    NUCLEUS_CORTEX_CONTACT,
    NUCLEUS_CYTOSOL_BOUNDARY,
    REGISTERED_OWNERS,
    SCHEMA_SOURCE,
    SURFACE_POROUS_TRANSFER,
    SURFACE_PRESSURE_FLUID_CONNECTORS,
    build_contract,
    continuous_contracts,
    contracts_by_name,
    declared_but_unimplemented,
    internal_contracts,
    kinetic_contracts,
)
```

Aleph target file: `aleph/state/connectors_surface_traction.py`.
Controls: `tests/state/test_connectors_surface_traction.py`.

This entry covers the eight Group-A contracts and the shared validation machinery.
`ALEPH-PORT-1902` covers the four Group-B contracts in the same module.

`ConnectorContract` and `ConnectorFamily` are declared by lane L15 in `aleph/state/schema.py`.
At the time this entry was written that module was **not on disk**, so the L19 module carries a
fallback definition of identical shape and exports `SCHEMA_SOURCE` to say which one is live. The
field names, order and defaults are pinned in `CONTRACT_FIELD_ORDER` /
`CONTRACT_FIELD_DEFAULTS` and asserted, so the two definitions cannot silently diverge.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — verified `git rev-parse HEAD`, branch `codex/ff-ac-codex` |
| Source path | `ffn_sim/ac/engine/contracts.py`, `dispatch.py`, `surface_body.py`, `erm_cortex_connector.py`, `erm_cortex_slice.py`, `interior_column_slice.py`, `cytosol_connected.py`, `fluid_core.py`, `medium_exterior.py`, `ecm_world.py`, `cortex_population.py`, `composed_native.py`, `ac/cell/erm_tether.py`, `ac/fluid/biot_substrate.py`, plus `scripts/ac_gate_b_erm_cortex_native.py`, `scripts/ac_gate_b_interior_column_native.py`, `scripts/ac_gate_t10_medium_native.py` |
| Source symbol(s) | `ConnectorContract` rows at `contracts.py:564-583, 756-810`; `CanonicalFacadeClaim` at `dispatch.py:105-134`; `SurfaceStepBindings` (`surface_body.py:193-225`); `SurfaceBody.accumulate_mechanics` (`surface_body.py:364-389`); `ErmCortexConnector` (`erm_cortex_connector.py:114`); `CortexCytosolPorousTransfer` (`interior_column_slice.py:283`); `MembranePressureFluxAdjointBoundary` (`cytosol_connected.py:141`); `NucleusPressureAdjointBoundary` (`fluid_core.py:745`); `ExteriorStokesMedium` (`medium_exterior.py:393`); `ECMMembraneContact` Protocol (`ecm_world.py:320-346`) |
| Read from | **working tree**, at a checkout whose `HEAD` is `be0e5876` |
| Working tree == commit? | **yes, for every cited file.** `git diff be0e5876 --stat -- <path>` returned empty output for all of them. The tree carries **31** uncommitted changes overall (`git status --short`); of the engine files it touches — `cortex_state.py`, `microtubule_rig.py`, `observe/artifact.py` — only `microtubule_rig.py` is cited here at all, and only in the `accumulate_*` census in §6, where the two symbols counted were checked present in both revisions. |

> The read-provenance rows are not bureaucracy. Lane L9 found that at least one earlier audit
> candidate differed from the commit it was named against. Here the working tree and the commit
> agree on every cited path, and that was checked per path rather than assumed from a clean
> `git status` — which the tree does not have.

## 3. Why source-derived porting beats clean-room

**For code, it does not, and no code was taken.** The eight contracts are transcribed from Aleph's
own Appendix A mechanics registry (`docs/manuscripts/extracted/APPENDIX_A_mechanics_registry.txt`,
Group A) and restated in Aleph's vocabulary. Endpoint names, roles and mechanisms come from that
document, not from the reference tree. No identifier, no class name, no dispatch-slot name and no
prose crossed.

**One thing crossed as an argument, and only because it was verified by reading the source rather
than believed from a summary:**

> A compressive membrane/cortex contact was declared as a connector, was named in a facade's
> canonical claim list, and has no evaluation site anywhere. The structural coverage gate passes
> for it, because the gate proves the connector was *claimed by something callable* — not that
> anything calls it.

That is a case analysis a competent author would not independently think to guard against, it was
discovered by somebody else's build, and `ports/TEMPLATE.md` §3 names exactly that as a legitimate
reason to let something cross. The Aleph consequence is stated as physics and as a schema
obligation, not as a bug report:

- A tensile molecular tether and a compressive contact are **two connectors**, not two modes of
  one. Collapsing them removes the only path outward turgor has to the membrane.
- A registry must be able to answer "which declared connectors have no evaluator" as a **list**.
  `declared_but_unimplemented()` is that list, and all twelve L19 contracts are on it.

## 4. Physical or mathematical law represented

**There is no law here, and saying so is the substance of this entry.** What crosses is a
*contract*: a typed declaration of what a mechanical connection between two owners must satisfy.
Three of its clauses are physics-derived and are enforced at construction rather than documented:

1. **Bidirectionality is not optional.** A connector that applies a force to one endpoint and no
   reaction to the other is a momentum source. `dp/dt = ΣF` over the pair would then be non-zero
   for an isolated pair, so `bidirectional=False` is refused rather than configured.
2. **The adjoint reaction is not optional.** For a transfer built from an interpolation operator
   `S` (fluid → structure), the reaction must be applied through `Sᵀ`. If spread and interpolate
   use different stencils or different weights the pair does work on nothing: momentum enters or
   leaves the fluid at a rate no per-owner force check can see, because each owner's own books
   still balance. This is why `surface_porous_transfer` states "spread must be the transpose of
   interpolate" in its interpretation and not in a comment somewhere else.
3. **Commit semantics follow from what a connector owns.** A `K` connector owns a bond, a capture
   or a topology. Those changes are irreversible, so they may commit only after a step is accepted;
   committing inside a candidate step that is then rejected leaves accepted state describing a
   configuration that no longer exists, and a force check computed *from* that topology cannot
   detect it. A `C` connector owns no such state, so `commit_on_accept` must be `False` — a
   reserved commit slot that is never used is a slot somebody will eventually put something in.

**Sign restriction, which is the physics this group is about.** For a scalar separation `s` and a
rest separation `s₀`, a tether stores `U = (k/2)·max(0, s − s₀)²` and a non-adhesive contact stores
`U = (k/2)·max(0, s₀ − s)²`. Each is `C¹` at `s = s₀` and each is **identically zero** on its
inactive branch — so a compressed tether carries exactly `0`, not a small number. The two branches
are disjoint, so no single element can supply both, which is the derivation behind
"`membrane_erm_cortex` and `membrane_cortex_contact` are not duplicates". The same argument applies
unchanged at the nuclear envelope: LINC tethers occupy the tensile branch, so a compressive channel
requires `nucleus_cortex_contact` as its own connector.

**Moving boundary versus contact.** `membrane_cytosol_boundary` and `nucleus_cytosol_boundary` have
no gap and therefore no engagement gate: the surface's normal velocity *is* the fluid's boundary
condition, and the fluid's pressure and viscous traction load the surface. There is no `s₀` to
choose and no stiffness with a physical referent. Modelling them as contacts would invent both.

## 5. Units, domains, singular cases, invariants

This module holds **no numeric quantity at all** — it is a registry of declarations. The units
table is therefore about what the *declared* connectors will carry when evaluators exist, and is
recorded so a later lane cannot pick a different convention silently.

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| separation, standoff | µm | m | `> 0` for contact and tether laws |
| force on an endpoint | pN | N | finite |
| traction on a quadrature site | pN/µm² | Pa | finite |
| pressure (osmotic, cytosolic) | pN/µm² | Pa | finite |
| stored energy | pN·µm | J | `≥ 0` for a unilateral element |
| drag coefficient | pN·s/µm | N·s/m | `> 0` |

Singular and boundary cases, each with the behaviour Aleph requires:

- **An endpoint naming an owner that is not registered** — refused with
  `ConnectorContractError`. A contract with one real end looks like coverage of a load path that
  cannot exist.
- **An empty endpoint role** — refused. The role says what quantity crosses the boundary at that
  end; without it there is nothing checkable to declare.
- **`bidirectional=False` / `adjoint_required=False`** — refused, as §4.1 and §4.2.
- **A `K` contract with `commit_on_accept=False`** — refused.
- **A `C` contract with `commit_on_accept=True`** — refused.
- **`internal_to` disagreeing with the endpoints in either direction** — refused. Unset on a
  same-owner connector presents it as a cross-owner load path; set on a cross-owner connector
  hides it from cross-owner dispatch.
- **A duplicate connector name in a registry projection** — refused rather than last-wins, because
  which of the two survives would otherwise be a dictionary-ordering accident.
- **A frozen contract** — mutation raises. A registry that was audited must be the registry that
  ran.

Invariants, each with the test that asserts it (all in
`tests/state/test_connectors_surface_traction.py`):

- **I1.** All eight Group-A contracts are constructible and registered in Appendix-A order.
  `::test_all_twelve_contracts_are_constructible_and_registered`
- **I2.** Every endpoint is a registered owner. `::test_every_endpoint_is_a_registered_owner`
- **I3.** Every contract is bidirectional and requires the adjoint.
  `::test_every_contract_is_bidirectional_and_requires_the_adjoint`
- **I4.** Every `K` contract commits only on an accepted step, and there are exactly four.
  `::test_every_kinetic_connector_commits_only_on_an_accepted_step`
- **I5.** No `C` contract claims committable state.
  `::test_no_continuous_connector_claims_committable_state`
- **I6.** The ERM tether is declared tensile-only and is not described as compressive.
  `::test_the_erm_tether_is_declared_tensile_only`
- **I7.** The membrane/cortex contact is declared compressive-only, carries the osmotic-pressure
  path, and says it is not a duplicate.
  `::test_the_membrane_cortex_contact_is_declared_compressive_only`
- **I8.** The two are distinct contracts with the same endpoints, different families and different
  interpretation texts.
  `::test_the_tensile_and_compressive_paths_are_distinct_contracts`
- **I9.** `nucleus_cortex_contact` records that LINC tethers alone cannot supply the compressive
  channel. `::test_the_nucleus_cortex_contact_records_that_linc_tethers_cannot_supply_it`
- **I10.** The two fluid-boundary contracts declare kinematics-plus-traction and say explicitly
  they are not contacts.
  `::test_the_fluid_boundary_connectors_are_kinematics_and_traction_not_contacts`
- **I11.** `surface_porous_transfer` declares immersed drag with velocity interpolation and an
  adjoint reaction transfer.
  `::test_the_porous_transfer_is_immersed_drag_with_an_adjoint_reaction`
- **I12.** `membrane_ecm_contact` is contact and explicitly not integrin-mediated adhesion.
  `::test_the_membrane_ecm_contact_is_contact_and_not_adhesion`
- **I13.** `membrane_medium_traction` records that it requires an exterior Stokes solve and is
  `implemented=False`. `::test_the_medium_traction_records_that_it_needs_an_exterior_solve`
- **I14.** All twelve are `implemented=False` and `declared_but_unimplemented()` returns all of
  them. `::test_every_contract_is_honestly_unimplemented`
- **I15.** That projection is not a constant: a contract with `implemented=True` drops out of it.
  `::test_the_unimplemented_projection_shrinks_when_something_is_implemented`

## 6. Source evidence class and known retractions

A read-only audit of `be0e5876` classified each of the eight as `IMPLEMENTED` (a real force or
traction routine a production driver calls), `SEAMED` (a Protocol, a claim or a registration with
no production body reached) or `NAME_ONLY`. Evidence class of the audit itself: `AUDIT_READ` — the
source was read, and in the reference's own words where those words are quoted. No numeric value
was inherited, so there is nothing here that could be retracted numerically.

| Connector | Verdict at `be0e5876` | Decisive evidence |
|---|---|---|
| `membrane_erm_cortex` | **IMPLEMENTED** | `erm_cortex_connector.py:205-224` launches a real kernel; `ac/cell/erm_tether.py:276-277` does `atomic_add(force_m, −f)` / `atomic_add(force_c, +f)`; reached by `erm_cortex_slice.py:307-310` and `scripts/ac_gate_b_erm_cortex_native.py:200`. It is also unilateral — it returns early when tension `≤ 0`. |
| `membrane_cortex_contact` | **SEAMED** | No class, no Protocol method, no binding field, no call site. See §6.1. |
| `surface_porous_transfer` | **IMPLEMENTED** | `interior_column_slice.py:327-338` delegates to a real pressure-gradient scatter (`ac/fluid/biot_substrate.py:346`); called at `interior_column_slice.py:563`. |
| `membrane_cytosol_boundary` | **IMPLEMENTED** | `cytosol_connected.py:193-202` — a Kedem–Katchalsky source plus a real adjoint traction; called at `interior_column_slice.py:560`. |
| `nucleus_cytosol_boundary` | **IMPLEMENTED** | `fluid_core.py:791-803`; called at `interior_column_slice.py:561`. Claimed under a *different* facade (`dispatch.py:129-134`) from the other three. |
| `membrane_medium_traction` | **SEAMED (component runtime real)** | See §6.2 — this is where the brief I was given is wrong. |
| `membrane_ecm_contact` | **SEAMED** | `ecm_world.py:320-346` is a `Protocol` whose `accumulate_contact` (`:329-334`) is docstring-only; the only implementation in the tree is a test double at `tests/ac/engine/test_ecm_world.py:192`. |
| `nucleus_cortex_contact` | **SEAMED** | Contract at `contracts.py:779-789`, claim at `dispatch.py:120`, census at `cortex_population.py:92`. No class, no Protocol method, no binding field, no call site. |

Nothing in Group A is `NAME_ONLY`: every one of the eight has at least a contract row, a claim and
a census entry, which is precisely what makes the seams hard to see.

### 6.1 The declared-and-never-evaluated contact, confirmed with corrected line numbers

The claim I was handed was that `dispatch.py:105-128` assigns `membrane_cortex_contact` to a
facade method with no slot to call it, at `surface_body.py:194, 364-389`. **Confirmed, and the
exact lines are:**

- `dispatch.py:105-128` is the surface facade's canonical claim. The claimed method name is at
  `:108`; the connector tuple runs `:110-127`; `"membrane_cortex_contact"` is at **`:119`**.
- The bindings dataclass is decorated at `surface_body.py:193` and declared at `:194`. It has
  **exactly four fields**, `:197-200`. There is no contact slot — and none for
  `nucleus_cortex_contact` or `membrane_medium_traction` either.
- The facade's accumulation method, `surface_body.py:364-389`, makes **five** calls: two component
  accumulations and **three** connector dispatches (`:372-377`, `:378-381`, `:382-389`). Three of
  the six connector edges claimed at `dispatch.py:110-127` are dispatched; the three at `:119`,
  `:120` and `:126` are not.
- Stronger than the original claim: the same three bindings are all that appear in
  `transaction_participants` (`:391-407`), `snapshot_candidate` (`:409`), `rollback` (`:414`),
  `commit_irreversible` (`:419`) and `accumulate_ledger` (`:430-436`). The contact edges are absent
  from **every** lifecycle hook, not only from force assembly. So there is no candidate/commit path
  they could have been reached through either.
- The reference says so itself, at `contracts.py:804-805`: *"Still declared-only — this records
  WHICH connector owns the path, not that it runs."* The declaration is honest **in the contract
  row**; it is the coverage gate that reads it as covered.

### 6.2 A claim in my own brief that the source refutes

I was told `membrane_medium_traction` is *"declared only — no exterior solve exists behind it,
admitted inline in their dispatch"*. **The admission exists verbatim. The claim it makes is stale.**

The comment is at `dispatch.py:124-125`, immediately above the connector name at `:126`, and reads:
*"DECLARED ONLY — no exterior solve exists behind it until T10, so this claim reserves ownership
and does not assert that any medium traction is evaluated."*

But the exterior solve now exists. `medium_exterior.py` is a real regularised-Stokeslet
single-layer operator: a matvec kernel at `:143`, a matrix-free conjugate-gradient resistance solve
at `:534`, a force scatter at `:203`, and an accumulation entry point at `:620-662`. What is
missing is the **connector wiring**, not the physics — and two other admissions in the tree say so
accurately: `medium_exterior.py:49-51` ("No connector wiring") and
`scripts/ac_gate_t10_medium_native.py:517` ("declared and NOT dispatched; rung cannot exceed
CUDA_UNIT").

**The correct verdict is: component runtime IMPLEMENTED, connector edge SEAMED.** Recorded here
because a port ledger that repeated the brief would have carried a false statement about the source
into Aleph's own records, and the porting policy exists to stop exactly that. It also matters for
Aleph: the reason `MEMBRANE_MEDIUM_TRACTION.implemented` is `False` in Aleph is that **Aleph** has
no exterior solve, which is a true statement about Aleph and is the only kind this registry is
allowed to make.

### 6.3 Two cross-cutting findings, both about coverage rather than about any one connector

- **The exact-once dispatch manifest is satisfiable by a no-op.** `composed_native.py:296-303`
  defines a facade stand-in that, in the reference's own words, *"launches nothing"*, and
  `:453-463` maps every facade type outside a set of three onto it. So the surface facade's six
  connector claims pass the coverage check while evaluating nothing.
- **The bindings types are test-only.** Construction sites for the surface bindings exist only in
  `tests/ac/engine/test_surface_body.py` (`:218, 362, 369, 376, 383, 604`). Since the facade's
  accumulation method takes a bindings object as its only argument, **that facade's connector
  dispatch has never run outside the test suite.** The four `IMPLEMENTED` verdicts above therefore
  rest on a different lane entirely — two *slice* drivers (`erm_cortex_slice.py:300`,
  `interior_column_slice.py:546`) invoked from two CUDA gate scripts. Judged against the composed
  whole-cell pipeline instead, all eight Group-A connectors are seams.

  This is the finding with the longest reach, and it is a warning about Aleph's own coverage gate
  rather than about the reference: *a connector can be genuinely implemented and genuinely
  unreachable at the same time, depending on which entry point you call production.* Aleph's
  runtime already requires scheduled + reached + evaluated-something
  (`ALEPH-PORT-303`); this says the definition of "reached" must name the entry point.

**A claim in the brief that the source also refutes, recorded for the same reason as §6.2.** I was
told only three production `accumulate_*` implementations exist in the whole tree. The census at
`be0e5876` finds **179** `def accumulate_<suffix>` definitions and **98** bare `def accumulate(`,
**277** by AST, of which **181** are outside `tests/`; roughly 33 of those are docstring-only
Protocol declarations and about 6 are no-ops. Narrowed to *connector-level force or traction
accumulators with real numeric bodies*, the count is **9**; narrowed further to *those reached by a
live driver*, it is **4** (the four `IMPLEMENTED` rows above). Under no reading is it three. The
"three" figure appears to have been a count of one particular kind of accumulator that was then
repeated as a count of the tree, and it should not be cited again.

## 7. Independent oracle or derivation

The oracles are Aleph's own, and none of them is agreement with the reference:

- **Appendix A itself**, `docs/manuscripts/extracted/APPENDIX_A_mechanics_registry.txt` Group A and
  the connector legend, which is Aleph's manuscript rather than the reference's code. Endpoints,
  roles, mechanisms and the two mechanical interpretations are checked against it line by line.
- **`ALEPH-PD-002`**, which fixes the contract field set and the rule that the census is a
  registered manifest rather than the schema. `CONTRACT_FIELD_ORDER` and
  `CONTRACT_FIELD_DEFAULTS` are asserted against that decision, so a schema authored by another
  lane cannot drift from it without a test failing.
- **The unilateral identity** `dU/ds ≡ 0` on the inactive branch, which is what makes the
  tensile/compressive split a statement about disjoint branches rather than a preference. It is
  derived in §4 and is the reason I6–I9 assert the *presence and distinctness* of the two
  interpretations rather than a numeric tolerance.
- **Mutation.** Three deliberate mutations were applied to the module and each was caught: see
  §9.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_connectors_surface_traction.py::test_all_twelve_contracts_are_constructible_and_registered` | Eight Group-A and four Group-B contracts construct, in Appendix-A order, with the exact names. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_every_contract_carries_its_declared_endpoints_roles_and_family` | Every contract's family letter, both endpoints, both roles and its mechanism match an independently written expectation table. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_contract_shape_matches_the_registered_field_set` | The dataclass field names, order and defaults are the fourteen `ALEPH-PD-002` specifies, whichever module supplied the class. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_erm_tether_is_declared_tensile_only` | The ERM contract states tensile-only, states it is not a compression strut, names turgor, and does not describe itself as compressive. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_membrane_cortex_contact_is_declared_compressive_only` | The contact contract states compressive-only, names the osmotic-pressure path, and states it is not a duplicate of the tether. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_membrane_cortex_contact_records_why_the_honest_flag_matters` | The contract that the reference declared and never evaluated carries, in its own text, the reason it is `implemented=False` and names the projection that surfaces it. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_nucleus_cortex_contact_records_that_linc_tethers_cannot_supply_it` | The nuclear compressive channel is declared, and the tether-cannot-push argument is recorded at the connector. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_fluid_boundary_connectors_are_kinematics_and_traction_not_contacts` | Both fluid boundaries declare moving-boundary kinematics plus traction and state they are not contacts. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_porous_transfer_is_immersed_drag_with_an_adjoint_reaction` | Immersed drag, velocity interpolation, adjoint reaction, and the spread-is-the-transpose-of-interpolate requirement. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_every_contract_is_honestly_unimplemented` | All twelve are `implemented=False` and the unimplemented projection returns all twelve. |

## 9. Deliberately failing negative control

Every registry invariant has a control that constructs a violating contract and asserts the
refusal. `::test_the_validator_is_not_vacuous` accepts the same keyword arguments each of them
mutates, so a validator that raised unconditionally fails rather than looking thorough.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_a_one_way_connector_is_refused` | `bidirectional=False` is refused. A force with no reaction is a momentum source, not a connector. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_dropping_the_adjoint_reaction_is_refused` | `adjoint_required=False` is refused. Without the adjoint, force closure is a tautology. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_a_kinetic_connector_that_does_not_commit_on_accept_is_refused` | A `K` contract with `commit_on_accept=False` is refused, and the corrected version is accepted — so the refusal is about that one field. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_a_continuous_connector_that_claims_committable_state_is_refused` | A `C` contract may not reserve a commit slot it can never fill. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_an_unregistered_endpoint_owner_is_refused` | An endpoint outside the owner registry is refused rather than accepted and ignored. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_an_empty_endpoint_role_is_refused` | A contract with no endpoint role declares nothing checkable. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_an_empty_mechanism_is_refused` | Same, for the mechanism. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_an_unnamed_contract_is_refused` | A nameless contract cannot be looked up, dispatched or reported. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_an_unrecognised_family_is_refused` | Only `K` and `C` exist; a third family would have undefined commit semantics. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_an_empty_composite_group_is_refused` | A whitespace group name would silently create a group nothing else can join. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_a_duplicate_connector_name_is_refused` | Two contracts under one name are refused rather than resolved by dictionary order. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_contracts_are_frozen` | A registered contract cannot be mutated at runtime, so the registry that was audited is the registry that runs. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_the_unimplemented_projection_shrinks_when_something_is_implemented` | The honest-coverage projection is not a constant that happens to equal the registry. |

**Mutation evidence, run on 2026-07-30 rather than asserted.** Three mutations were applied to
`aleph/state/connectors_surface_traction.py`, each run, and the file restored byte-for-byte
afterwards:

| Mutation | Result |
|---|---|
| Neuter the `bidirectional` refusal (`if fields["bidirectional"] is not True:` → `if False:`) | **2 failed, 41 passed** — `test_a_one_way_connector_is_refused`, `test_a_refusal_is_catchable_as_a_value_error` |
| Neuter the kinetic commit refusal **and** drop `commit_on_accept=True` from one `K` contract | **2 failed, 41 passed** — `test_every_kinetic_connector_commits_only_on_an_accepted_step`, `test_a_kinetic_connector_that_does_not_commit_on_accept_is_refused` |
| Rewrite the ERM interpretation as a bidirectional strut ("it can push and pull") | **1 failed, 42 passed** — `test_the_erm_tether_is_declared_tensile_only` |

The third is the one that matters: it is the reference's own defect, planted in Aleph's tree, and it
was caught.

## 10. Numerical and precision envelope

**No floating-point arithmetic occurs in this module.** Every field is a string, a boolean, a
`StrEnum` member or `None`, so there is no working precision, no accumulation precision and no
conditioning to state — and that is the honest answer rather than an omission.

The controls are therefore exact rather than toleranced. Booleans are compared with `is True` /
`is False` rather than truthily, so `1` or `"yes"` in a flag position fails; family membership is
compared against enum members; the interpretation controls are substring containments on lower-cased
text, which is a presence assertion and not a similarity score. `::test_every_contract_states_a_mechanical_interpretation`
uses a length floor of 80 characters, which is a crude bound chosen only to refuse a one-word
placeholder — it is not a quality measure and is not presented as one.

Tolerances become relevant the moment any of these connectors acquires an evaluator. The one that
should be recorded now, because it is a property of the law and not of an implementation: the
inactive branch of a unilateral element must be asserted with `==` against exact zero, not with a
tolerance. A tolerance there admits a small wrong force in exactly the place where a small wrong
force is the whole defect.

## 11. Production-backend residency and transfer

Host-side Python only, and permanently so. This module is a registry of immutable frozen
dataclasses constructed at import time; nothing here is transferred to a device, and there is no
per-step cost because nothing here runs per step. It will be read by whatever assembles the
dispatch schedule — on the host, once — and the contracts themselves never cross the host/device
boundary.

When an evaluator lands for one of these connectors, the residency answer belongs to that
evaluator's ledger entry, not to this one. What this entry fixes is what must be transferred *per
connector* regardless of backend: the adjoint reaction. A transfer whose reaction is applied
through a different stencil than its interpolation cannot be repaired at the accounting layer.

`aleph/state/**` may not import `validation/**`, and this module imports only `dataclasses`, `enum`
and `typing`.

## 12. Comments and docstrings to discard

No source prose was carried, so none survives. The reference's vocabulary is deliberately
**absent** from `aleph/state/connectors_surface_traction.py`: no repository name, no package
namespace, no module paths, no facade or slot names, no gate labels, no class names, no branch
names. The module refers to "the predecessor project this one replaces" and describes the defect in
terms anyone can check against Aleph's own code with that repository absent from disk, which is
gate R1's requirement and master plan §6's second half.

Specifically not carried into `aleph/**`, though quoted here where quoting is the evidence:

- the facade method name that the coverage claim referred to;
- the bindings dataclass name and its field names;
- the gate/rung labels and the milestone tag in the reference's declared-only comments;
- the reference's contract-row wording, including the "declared-only" comment quoted in §6.1;
- the connector-family vocabulary of the reference's own enum, which differs from Appendix A's
  two-family legend.

What replaces them: §3, §4 and §6 of this entry, and the module docstring, which states the defect
as *"a declared load path that is never evaluated, passing a gate that counts registrations"* —
a description with no proper nouns in it.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`. The code named in §1 has landed and its controls pass — 43 tests in `tests/state/test_connectors_surface_traction.py`, plus the three mutations in §9 — but what is being registered here is a set of *contracts*, and a contract is accepted when the physics it declares has an evaluator and a control, not when the declaration parses. Every contract in this entry is `implemented=False`, so there is nothing to accept yet, and moving this row before then would be the declaration-versus-reality gap in Aleph's own ledger. |
| Reviewer | agent-proposed, **unratified**. No human has reviewed the tensile/compressive argument in §4, the audit verdicts in §6, or the two corrections in §6.2 and §6.3 — one of which contradicts the brief this lane was given. |
| Rollback | Delete `aleph/state/connectors_surface_traction.py` and `tests/state/test_connectors_surface_traction.py`. Nothing imports them yet, so nothing outside those two files breaks. What is lost is the record that eight declared load paths have no evaluator — which would revert to being an absence rather than a list. |

## 14. Honest limits

- **Unratified**, and one of its substantive claims contradicts the brief the lane was given
  (§6.2, `membrane_medium_traction`). If the brief is right and my read is wrong, §6.2 is the first
  thing to check.
- **No physics is established by this entry.** Nothing here evaluates a force, and no contract in
  it is evidence that the load path it declares can be assembled, is well posed, or is stable.
  `implemented=False` on all eight is the whole content of that statement.
- **The interpretation controls are text assertions.** They assert that the tensile-only and
  compressive-only statements are *present and distinct*. They cannot assert that an evaluator
  obeys them, because there is no evaluator. A future evaluator must ship its own exact-zero
  control; this entry does not substitute for it.
- **`SCHEMA_SOURCE` is `"local fallback"` as written.** Lane L15's `aleph/state/schema.py` was not
  on disk. The two definitions are asserted to agree in field names, order and defaults, but they
  have never actually been exercised against each other, and L15's `declared_but_unimplemented()`
  has never been called by this module. That reconciliation is untested and is the first integration
  risk in this lane.
- **The owner allowlist is a transcription.** `REGISTERED_OWNERS` is the fourteen top-level
  components of Appendix A, taken from that document with no independent check that fourteen is the
  right number or that these are the right cuts. `ALEPH-PD-002` explicitly refuses to treat the
  census as mandatory ontology, so this is one registered manifest and not a structure — but it is
  an inherited manifest, and it is carried at `ASSUMED`.
- **Six of the eight connectors name owners no Aleph module represents.** `cytosol`, `nucleus`,
  `ecm` and `extracellular_medium` have no state, no geometry and no mechanics in Aleph. Their
  contracts are declarations about owners that do not exist yet, which is exactly what
  `ALEPH-PD-002` asked for and is also a reason not to read this registry as a description of
  anything runnable.
- **The audit is a read, not a run.** The `IMPLEMENTED` / `SEAMED` verdicts in §6 come from reading
  the source and its import and construction graph. Nothing in the reference was executed for this
  entry, so a body that exists and is reached could still be numerically wrong, and a seam could
  conceivably be reached by a path the read missed. §6.3's finding — that the answer depends on
  which entry point counts as production — is a direct consequence of that limit.
- **No claim is made about the reference's present state.** The verdicts are pinned to
  `be0e5876`, and the tree carries 31 uncommitted changes. Every cited file was checked identical
  to that commit, but the tree will move.
- **`membrane_ecm_contact` is a candidate-search problem this entry ignores.** Contact between a
  surface and a collagen segment requires a proximity search, and there is no correspondence
  between membrane sites and segments that could be fixed by construction. The contract declares
  the connector; it says nothing about how the pairs are found, and that is a real gap between this
  declaration and anything implementable.
