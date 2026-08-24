# Port ledger — ALEPH-PORT-1501 · the CellState contract registry

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1501` |
| Lane | `L15 state contract` |
| Status | `ACCEPTED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Evidence class | STRUCTURAL. A schema is not a measurement and this entry claims no physics. |

---

## 1. Aleph API

The exact public surface this entry authorises. Aleph targets:
`aleph/state/schema.py`, `aleph/state/context.py`, `aleph/state/cell_state.py`,
`aleph/state/census.py`.

```python
from aleph.state.schema import (
    ScopeTag, ConnectorFamily,
    ComponentContract, ConnectorContract,
    StateFieldDeclaration, TopologySectorDeclaration, UnitDeclaration, MigrationRule,
    StateSchemaManifest, FrozenStateSchema, ForcePath, UnimplementedInventory,
    SchemaValidationError, DuplicateRegistrationError, ContextBoundaryError,
    stable_id,
)
from aleph.state.context import (
    ExperimentContextManifest, ExternalField, BoundaryCondition, ScheduleEntry,
    CONTEXT_OWNED_QUANTITIES,
)
from aleph.state.cell_state import (
    TopologySector, SectorJumpMap, TopologyAtlas, TopologyTransitionError,
    StateBlock, OwnerLocalState, CellStateSnapshot,
    RefusalKind, ScopeRefusal, ScopeRefusalError,
)
from aleph.state.census import (
    SEED_CENSUS_2026_07_29, CensusRegistration,
    register_census_group, seed_census_groups, register_seed_census, new_seed_manifest,
)
```

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/contracts.py` |
| Source file digest | `sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72` |
| Source symbol(s) | `ComponentContract`, `ConnectorContract`, `ConnectorFamily`, `ConnectorScope`, `CellArchitecture`, `reference_cell_architecture` |
| Read from | working tree, then confirmed identical to the commit |
| Working tree == commit? | `yes` — `git diff be0e5876 -- ffn_sim/ac/engine/contracts.py` is empty; `git status --porcelain` reports the path clean |
| Source test status | The reference's own suite passes on this file, and its tests exercise the dataclass guards only. No test there asserts that a *declared* connector was ever evaluated — which is the defect ALEPH-PORT-303 records and which §6 below expands on. |

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** Not one line, identifier set, field list, or docstring
crossed the boundary. This entry exists because the reference was *read*, and PLAN §0.2 requires a
record of that reading whether or not code moved.

The reference's registry is a fixed tuple of components plus a fixed tuple of connectors, validated
in `CellArchitecture.__post_init__` and produced by one hardcoded factory function. That shape is
precisely what ALEPH-PD-002 and master plan §5.3 forbid: it makes the census the *structure* of the
schema, so adding an owner means editing a schema rather than extending a manifest. Aleph therefore
inverts the ownership — an open registry with a `finalize()` adjudication step — which is a different
design, not a re-expression of the same one. Re-deriving was the cheaper path as well as the required
one.

Two ideas were confirmed by the reading and are carried as *findings restated in Aleph's words*,
never as text:

1. A connector that is not bidirectional is not a load path, so refusing `bidirectional=False` at
   construction is right. Aleph adopts the refusal and states its own reason: a one-way force is a
   momentum source, and a momentum source inside a closed mechanical assembly makes the force-closure
   check unfalsifiable.
2. A composite series joint must be dispatched as a group, because evaluating a series of two
   compliances as two independent springs computes the wrong stiffness — the series compliance adds,
   the parallel stiffness does not.

## 4. Physical or mathematical law represented

**There is no physical law here, and saying so is the point.** What crosses the boundary is an
*argument about ownership*, and the entry would be dishonest if it implied otherwise.

The mathematical content is a typed graph and two hash functions:

* The owner graph is `G = (V, E, tau_V, tau_E)` with `tau_V : V -> {E,I,H,B,X}` and
  `tau_E : E -> {K,C}`. Every declared load path is a walk in `G` whose edges are named connectors,
  which is the formal content of the ontology rule that co-location is never a mechanical
  connection: adjacency in space is not adjacency in `G`.
* `force_path` is breadth-first search on `G` restricted to vertices with `tau_V = E`. BFS returns a
  path of minimal edge count, which for a mechanical graph means the fewest declared connectors a
  load must traverse. It is a reachability statement and it carries no magnitudes; a returned path
  proves the transmission route was declared, never that any force is transmitted.
* The scope partition is an identity that is derived rather than asserted:
  `owned_state(v) != ()` implies `tau_V(v) in {E, I}`. An `H` field is an effective law and holds no
  degrees of freedom; a `B` condition lives in the context manifest; an `X` entity is outside the
  boundary. So a stateless-scope entity that owns state is a contradiction in terms rather than a
  policy violation, and `finalize()` reports it as one.
* The schema hash and the context hash are SHA-256 over Aleph's existing canonical JSON encoding,
  domain-separated by a contract identifier so that a hash computed under a different encoding can
  never share a namespace with one computed under this encoding.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| Declared field unit | any registered symbol (`pN`, `um`, `s`, `pN*um`, `pN/um`, `K`, `1`, `count`, ...) | as declared in `UnitDeclaration.si_unit` | must resolve in the unit registry at `finalize()` |
| Temperature (context) | `K` | `K` | `> 0`; no default is available, because `kB*T` inherits an `ASSUMED` evidence class and no entry point may supply it silently |
| Topology generation | `count` | dimensionless | monotone, `>= 0` |
| Schema / context hash | `sha256:<64 hex>` | — | opaque identity |

Singular and boundary cases, each with the behaviour Aleph requires:

- An empty manifest finalizes successfully and hashes to a stable value. It is a legitimate state:
  six lanes register into one manifest and the empty one is what they start from.
- A component named but never populated in a snapshot is **not** an out-of-scope component. The two
  produce different refusal kinds (`DECLARED_BUT_ABSENT` vs `OUT_OF_DECLARED_SCOPE`) because
  conflating them is how a missing owner becomes a zero.
- A connector whose two endpoints are the same owner is not an error — it is an internal joint, and
  it must say so via `internal_to`. Silence there is the error.
- A `composite_group` with one member is refused: a series of one is not a series, and a group of one
  is the shape a half-finished split leaves behind.
- Registration is deliberately permissive and `finalize()` is strict. A manifest half-built by six
  independent lanes is legitimately inconsistent; adjudicating at every `register_*` call would make
  registration order load-bearing.

Invariants that must hold, each with the test that asserts it:

- **I1.** Every connector endpoint resolves to a registered component —
  `test_unknown_endpoint_is_refused_at_finalize`.
- **I2.** Every connector endpoint has scope `E` — `test_internal_component_cannot_be_a_connector_endpoint`,
  `test_boundary_component_cannot_be_a_connector_endpoint`.
- **I3.** An `I` component has a registered `E` parent — `test_internal_without_parent_is_refused`,
  `test_internal_parent_must_be_explicit`.
- **I4.** `H` / `B` / `X` own no state — `test_stateless_scope_owning_state_is_refused`.
- **I5.** A same-owner connector sets `internal_to` — `test_self_edge_without_internal_to_is_refused`.
- **I6.** Every `K` connector has `commit_on_accept=True` —
  `test_kinetic_connector_must_commit_on_accept`.
- **I7.** `bidirectional=False` is refused at construction —
  `test_one_way_connector_is_refused_at_construction`.
- **I8.** Duplicate names are refused at registration —
  `test_duplicate_component_name_is_refused`, `test_duplicate_connector_name_is_refused`.
- **I9.** A composite group is returned whole — `test_mechanical_group_returns_the_whole_series`.
- **I10.** A context quantity cannot become a component's owned state —
  `test_context_quantity_cannot_be_registered_as_owned_state`.
- **I11.** The accepted-state digest sees RNG position — `test_snapshot_digest_sees_rng_position`.
- **I12.** A sector change requires a declared jump map — `test_undeclared_sector_transition_is_refused`.

## 6. Source evidence class and known retractions

Looked at: the reference file itself at `be0e5876`, the audit records under `ports/audit/`, the
port-ledger entries `ALEPH-PORT-302` and `ALEPH-PORT-303`, and `PLAN.md` §1.1 and §7.

What the reference claims for this code: nothing quantitative. It is a build-time declaration layer
and its own artifacts label the assembled world `evidence: CENSUS-WIRED` with
`quantitative_claim_status: BLOCKED`.

The defect that matters, already recorded by this project: the reference registry declares 36
connectors, of which roughly ten transfer force, nineteen are protocol seams and three are names
only — and at least one declared connector satisfies a dispatch coverage gate while having no
evaluation slot at all. **The census is far ahead of the physics, and the gate could not see it.**
That is why `implemented: bool = False` is a field on both contracts here and why
`declared_but_unimplemented()` exists as a first-class query rather than as a report someone might
run. A contract in Aleph asserts that something was *declared*; it asserts nothing about evaluation.

`bootstrap/registries/legacy_ffn_mechanics_snapshot.json` carries a projection of the reference
registry with `must_not_seed_aleph: true` and an open provenance defect (`ALEPH-PROV-001`: the
projection was taken without capturing an immutable commit). No row of it is used here. The seed
census content is owned by the six registry lanes and is deliberately absent from
`aleph/state/census.py`.

No retraction of the reference's dataclass guards was found, and its `bidirectional=False` refusal
is the one guard this reading endorses on its own merits.

## 7. Independent oracle or derivation

A schema has no analytic oracle, so the check is **mutation**: every validation must be demonstrated
able to fail. `tests/state/test_schema_validation_can_fail.py` builds, for each rule, a manifest that
violates exactly that rule and asserts the specific problem code is reported — and asserts the same
manifest with the violation repaired finalizes clean. A validation suite that only ever sees valid
input is decoration, which is the failure mode this project was restarted over.

Three further independent checks:

* The schema hash is invariant to registration order and changes on any content change
  (`test_schema_hash_is_order_invariant`, `test_schema_hash_changes_when_a_field_changes`). Order
  invariance is what allows six lanes to register concurrently and still agree on an identity.
* `force_path` is checked against a hand-computed path on a graph small enough to enumerate by eye,
  and against a deliberately disconnected graph where it must return `None` rather than a guess.
* The accepted-state digest is **not** re-implemented. `CellStateSnapshot` calls
  `aleph.artifacts.digest.accepted_state_digest`, which already covers schema hash, context hash,
  topology generation, owner state, RNG position and clock, and already carries the counterexample
  proving two states differing only in RNG position hash differently. Writing a second digest would
  have created a second answer to "what is this state", which is the defect the first one exists to
  prevent.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_schema.py::test_a_valid_manifest_finalizes_and_answers_graph_queries` | A 9-component / 8-connector reference manifest — 5 `E`, 1 `I`, 1 `H`, 1 `B`, 1 `X`, one kinetic edge, one internal self edge, one two-member series — finalizes with zero problems; `force_path("membrane","nucleus")` returns the unique 3-owner walk with its 2 connectors named in order; `mechanical_group` returns both members of the series; `owners_of_scope("E")` returns exactly the 5 explicit owners; the internal self edge is absent from `neighbors("cortex")`; and the schema hash is stable across two independent build orders. |
| Positive | `tests/state/test_cell_state.py::test_snapshot_digest_matches_the_artifact_digest` | The snapshot's `accepted_state_digest` equals `aleph.artifacts.digest.accepted_state_digest` called directly with the same components — so the snapshot is a caller of the one digest, not a second implementation. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_schema_validation_can_fail.py::test_every_declared_validation_is_provably_able_to_fail` | Iterates the declared rule table; each entry's mutated manifest must report that rule's code and the repaired manifest must not. If any validation stops being able to fail, this test fails. |
| Negative (must fail) | `tests/state/test_schema_validation_can_fail.py::test_a_declared_but_unimplemented_connector_cannot_pass_a_coverage_assertion` | A manifest whose connector is `implemented=False` finalizes clean — a contract may be unimplemented — but `declared_but_unimplemented()` names it, so an assertion that coverage is complete fails. This is the reference's dispatch-gate defect reproduced in Aleph's own vocabulary and caught. |
| Negative (must fail) | `tests/state/test_cell_state.py::test_out_of_scope_owner_refuses_instead_of_returning_zero` | Asking for an owner outside the declared scope raises a refusal that names the scope. The test also asserts the refusal is not falsy-equal to a zero-length state, so a caller cannot treat it as empty data. |
| Negative (must fail) | `tests/state/test_context.py::test_context_quantity_cannot_be_registered_as_owned_state` | Registering `temperature_K` as a component's owned state is refused, naming the context manifest as its home. |

## 10. Numerical and precision envelope

Working precision: this layer holds no floating-point physics. The only numeric values it touches
are declared shapes and populations (Python `int`, arbitrary precision, refused if negative) and the
context temperature (`float`, refused if non-finite or non-positive).

Hashing: SHA-256 over `aleph.artifacts.stamp.canonical_json`, whose float encoding is `repr()` — the
shortest decimal that round-trips to the identical double. So two contexts one ULP apart in
temperature hash differently, deliberately: they are different experiments. `-0.0` and `0.0` also
hash differently, and NaN and infinity are refused before they can acquire an identity.

The tolerance on every control in §8 and §9 is **exact equality**, and that is the right number
rather than a looser one because every assertion here is over strings, tuples, integers and hashes.
There is no discretisation and therefore nothing for a tolerance to absorb. The one place a
tolerance would have been tempting — comparing two schema hashes — is exactly the place where a
tolerance would be meaningless.

Outside the envelope: an unregistered unit, an unresolvable endpoint, a negative population and a
non-finite temperature all refuse. None of them degrades to a default.

## 11. Production-backend residency and transfer

Host only, and permanently. This is a declaration and validation layer: it runs once at assembly
time, produces two hashes and a frozen graph, and never appears in a step loop. Nothing here is
transferred to a device and no host round-trip is introduced per step.

`OwnerLocalState.blocks` holds the *shape, unit, name and mask* of each state block and is the
host-side description of arrays that will live on the device under the ratified backend. The arrays
themselves are numpy on the host in this layer; the block description is what survives to the device
side, so a future device-resident owner state changes where `StateBlock.array` points without
changing this contract. `CellStateSnapshot` is built from accepted state after a step has committed,
so it is off the hot path by construction.

No GPU is required, none was used, and no code path here can launch one.

## 12. Comments and docstrings to discard

Discarded in full — nothing from the source's prose survives in any form:

- The pressure-ownership decision note attached to the reference's envelope constant, including its
  option labels and its internal decision dates. Aleph's turgor ownership is settled by
  `ALEPH-PORT-1102` on Aleph's own measurement.
- The reference's gate vocabulary, evidence-rung names, and `quantitative_claim` labels. Aleph's
  evidence ladder is `aleph/evidence/` and is a two-axis construction, not this one.
- Every reference module path, symbol citation and working-tree assumption. `aleph/state/**` must be
  readable with no other repository checked out.
- The reference's role vocabulary (`SURFACE_BODY`, `FLUID_VOLUME`, `CORE_BODY`, ...) and its
  `ConnectorScope` enum. Aleph's `mechanical_role` is free text carried for the reader, and internal
  scope is expressed by `internal_to` naming the owner rather than by a second enum that can
  disagree with the endpoints.
- Its factory function's claim that the graph "declares allowed component-level load paths". Aleph's
  `force_path` docstring states the weaker and true thing: a returned path proves a route was
  declared and proves nothing about evaluation.

What replaces them: new prose written for Aleph's situation, stating in each module why the census is
data rather than structure, and carrying `declared_but_unimplemented()` as the standing answer to
"how do we know the registry has not run ahead of the physics again".

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | 2026-07-30. All lane L15 tests pass and the full suite stays green; every validation in §5 is demonstrated able to fail by the mutation table in §7. Exact-equality assertions throughout. |
| Reviewer | **Agent-proposed (lane L15). Not ratified by the PI.** Two design decisions need PI review and are named in §14. |
| Rollback | Delete `aleph/state/{schema,context,cell_state,census}.py` and `tests/state/test_{schema,schema_validation_can_fail,context,cell_state,census}.py`. Nothing else imports them yet, so the existing 168 geometry tests and the vertical are unaffected. What breaks is the six registry lanes and the tensor-network and neural lanes, all of which are unstartable without this contract — which is the whole content of ALEPH-PD-002. |

## 14. Honest limits

What this entry does **not** establish:

- **It establishes no physics.** Every contract here is a declaration. `implemented=False` is the
  default on both contract types precisely so that nothing in this layer can be mistaken for an
  evaluated force path.
- **`E`-only connector endpoints is a design decision, not a derived necessity.** The manuscript
  ontology says an `I` entity is "not an independent connector endpoint" and a `B` condition "lives
  in the context manifest", and this layer enforces both as hard refusals. That is the documented
  discipline, but it constrains the six registry lanes: a lane that wants a connector to a far-field
  clamp must model the clamp as a context boundary condition instead. **This needs PI review**, and
  loosening it later is an additive change to the validation table rather than a schema change.
- **Requiring `reentry_condition` on `H` as well as on `X` is stricter than the manuscript's letter.**
  The manuscript demands "declared inadequacy" for `H`; this layer reads that as also requiring a
  statement of what would force promotion. **This needs PI review.**
- The context-quantity guard is a **name-based** check against `CONTEXT_OWNED_QUANTITIES` plus an
  exact-key check against a bound context. A lane that spells the same physical quantity differently
  (`temp_kelvin` for `temperature_K`) defeats the name list. The bound-context check is the strong
  half and it only fires once a context is bound.
- No migration has ever been executed. `MigrationRule` records intent and is validated for
  completeness; the machinery that would replay a schema change across stored snapshots does not
  exist and is not claimed.
- The seed census content is **absent by design** and this entry therefore says nothing about
  whether the 14/36 registry is right, complete, or physical. It says only that registering it is an
  extension of a manifest rather than a change to a schema.
- `TopologyAtlas` validates that a jump map's population delta is consistent with the two sectors it
  joins. It does not simulate a transition, does not allocate anything, and has never been driven by
  a stepper. Whether the ratified transaction layer commits a sector change the way this contract
  describes is **untested across the seam** and needs the runtime lane.
- `citation_status` defaults to `"UNSOURCED"` on every component, and nothing in this layer promotes
  it. Sourcing the census is the registry lanes' work and no count here should be read as sourced.
