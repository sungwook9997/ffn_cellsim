# The Aleph tensor language and the MTG-PN reader stack — implementation specification

| Field | Value |
|---|---|
| Status | `SPEC` — the shared text the representation and learning lanes implement |
| Written | 2026-07-30 12:40 KST |
| Authority for this document | `ALEPH-PD-002` (PI): design work "lands as design documents"; plus the PI's 12:31 instruction to specify the tensor language and the reader design |
| May this exist as code? | **Yes, ratified 14:24** — `docs/decisions/RATIFICATION-2026-07-30-seal-is-a-computation-boundary.md`. Declaration and planning are open; qualification and training stay sealed until R5. The withdrawn `ALEPH-PD-003` stays withdrawn: its substance was ratified through the front door, its method was not. |
| Implementation state | complete. `aleph/represent/`: `index`, `blocks`, `cuts`, `families`, `plan`, `errors`, `qualification`. `aleph/learn/`: `outputs`, `lanes`, `architecture`, `losses`, `evaluation`, `training`. Full suite 2669, zero failures. |
| Science source | manuscript §3, §3.2, §4, §10, §10.1–10.3 — extracted to `docs/manuscripts/extracted/`. **Reference, never authority.** |
| Domain | the full compartment census: 14 components / 36 connectors, as *one registered manifest* |
| Scope of this document | L0 declaration and L1 planning. L2 qualification is specified so that it can be refused precisely, and is not built. |

---

## 0. What this document is for, in one paragraph

A tensor network is a factorization over a state space, so the language that describes the
factorization has to be written in the vocabulary of the state space rather than in the vocabulary of
arrays. Everything below therefore names owners, connectors, sectors, scopes and units — the objects
the `CellState` contract already defines — and nothing below names a rank, a width, or a tolerance
value. Where a number would go, there is a symbol and a slot. The slots are filled at L2, after an
accepted native trace exists; the language is finished before then, and that is the point: **training,
when it comes, is expressed in this language rather than bolted beside it.**

---

## 1. Three levels, and the seal between them

| Level | Question it answers | Status | Enforcement |
|---|---|---|---|
| **L0 declaration** | what is a legal index, block, bond, cut, lane, loss, refusal? | OPEN | — |
| **L1 plan** | for *this* frozen schema, which family factors which block, in what axis order, across which cuts, at what symbolic cost? | OPEN | plan is hashable and reproducible |
| **L2 qualification** | how large is each error channel against an accepted native trace, and does it pass? | **SEALED** | every entry point raises `RepresentationSealed`, gate `R5` named in the message |

The boundary is enforced statically, and as of the 14:24 ratification it is in force: **no module
under `aleph/represent/` or
`aleph/learn/` may import `numpy`, `torch`, `jax`, `scipy`, or `warp`**, at any scope. A pure-Python
declaration cannot fit a rank or train a weight, so the seal holds by construction rather than by
review. A second rule bans magic numbers: no rank, width, token count, or threshold literal may appear
in either package (§10.3 — provisional hyperparameters are preregistered only after pilot
measurements). Structural bounds that are *definitions* rather than tuning (`ndim >= 1`) are exempt and
must be named as such.

---

## 2. The index algebra — the smallest unit of the language

Everything is built from one object. An index (a tensor leg) is **not** an integer dimension; it is a
declaration of what varies along an axis and who owns that variation.

```python
@dataclass(frozen=True, slots=True)
class Index:
    name: str                     # stable, unique within a TensorDeclaration
    owner: str                    # census component or connector name; "" only when kind is BOND
    kind: IndexKind               # see below
    extent: Extent                # SYMBOLIC. never an int
    unit: str                     # a symbol resolvable in the unit registry; "1" for dimensionless
    basis: BasisKind              # NODAL | QUADRATURE | SPECTRAL | MODAL | CATEGORICAL | NONE
    sector_dependent: bool        # does the extent change across a topology jump?
    context_conditioned: bool     # is this a B-scope conditioning axis rather than state?
```

`IndexKind` — `ENTITY` (particles, filaments, vertices, elements), `SPATIAL` (a component of a vector
or tensor field, extent 2 or 3 by definition), `FIELD` (a named field channel), `TIME`, `SECTOR`
(a topology sector label), `MODE` (a spectral or reduced coordinate), `EVENT` (a position in a marked
stream), `ENDPOINT` (a connector endpoint role, extent 2 for a binary connector), `BOND` (an internal
contraction leg introduced by a factorization, owned by no one).

`Extent` is symbolic: `Extent.of("n_vertices", owner="membrane")`, `Extent.literal_dimension(3)` for a
spatial dimension, `Extent.unbounded()` for an append-only event axis. There is no `Extent.of(4096)`.
An extent that is unknown at schema time is not an error; an extent that is *asserted* as a number is.

### 2.1 Invariants — each is a test, each has a named refusal

| # | Invariant | Refusal when violated |
|---|---|---|
| I1 | Every non-`BOND` index names exactly one owner, and that owner exists in the frozen schema | `UNKNOWN_OWNER` |
| I2 | An index whose owner has scope `H`, `B` or `X` cannot be state | `SCOPE_HAS_NO_STATE` |
| I3 | A `B`-scope quantity may appear only as `context_conditioned=True`, never as a state axis | `CONTEXT_AS_STATE` |
| I4 | Two indices may be contracted only if unit, basis and extent symbol all agree | `UNIT_MISMATCH` / `BASIS_MISMATCH` / `EXTENT_MISMATCH` |
| I5 | A `sector_dependent` index may not be an axis of a factorization declared over more than one sector | `SECTOR_EXTENT_UNKNOWN` |
| I6 | An index for a declared-but-absent owner is a refusal, never an extent of zero | `ABSENT_NOT_ZERO` |
| I7 | Index names are unique within a declaration and stable across replans | `DUPLICATE_INDEX` |

I2 and I3 are where the E/I/H/B/X discipline of §4 becomes machine-checkable in the representation
layer instead of remaining prose. I6 is the manuscript's rule that absence is never encoded as a zero
measurement, hoisted from the observation layer up into the algebra, where it costs nothing to enforce.

---

## 3. State blocks — the census enters as a manifest, not as a structure

§10 partitions the state into six blocks. The census supplies *instances*; the blocks are the
structure. Adding a fifteenth component extends a manifest and touches no code here — that is
`ALEPH-PD-002`'s requirement carried into the representation layer.

| Block | Native form (§10) | Which census entries land here |
|---|---|---|
| `OWNER_PARTICLES` | stable local IDs, coordinates, velocities, forces, binding/topology state | components whose `owned_state` holds per-entity arrays: cortex, sf_arc, lamellipodium, filopodium, MT, IF, NMII, FA |
| `MESH_SURFACE` | vertices, elements, material fields, connectivity, boundary labels | membrane, nucleus envelope, ECM network, substrate surface where explicit |
| `CONTINUUM_FIELD` | grid/FEM coefficients with units, domain, boundary conditions | cytosol pressure and flow, homogenized drag fields where the law owns coefficients |
| `EVENT_STREAM` | append-only marked stream with accepted step and owner | binding/unbinding, severing, nucleation, topology jumps |
| `CONNECTOR` | endpoint references, force/work/energy/dissipation, kinetic state, active-power ledger | all 36 connectors, by family |
| `OBSERVATION` | raw protocol output plus calibration and analysis lineage | observation operators declared in `aleph/observe/` |

**Block assignment is derived, not hand-written.** A component's `owned_state` field names, together
with its `scope`, determine its block(s); an explicit override manifest exists for entries the
derivation cannot classify, and every override records why. A component may span blocks (the membrane
owns a mesh *and* participates in events); a component may own no block at all.

Non-state scopes still register, because declaring what is not modelled is the discipline:

- `H` → `HomogenizedLawSlot`: carries the effective law's name, its declared inadequacy, and its
  re-entry condition. It owns no index and can be a *coefficient* of one.
- `B` → `ContextConditioningAxis`: lives in `ExperimentContextManifest`, enters the network as a
  conditioning input, never as a factorized axis.
- `X` → `ExcludedRegion`: a typed record with its re-entry condition, so "what this model cannot
  answer" stays a query rather than a comment.

---

## 4. Graph cuts and bond structure — derived from the connector graph, never from space

This is the single most important rule in the document, and it comes straight from §4: **co-location
is never a mechanical connection.** Two owners that are adjacent in space but share no declared
connector have no interaction to factorize across, so a cut between them costs nothing; two owners on
opposite sides of the cell that share a connector are strongly coupled. A factorization that used
spatial adjacency for its bond structure would be optimizing the wrong graph.

A cut is a partition of the owner graph `G = (V, E, tau_V, tau_E)`. Its declared properties:

```python
@dataclass(frozen=True, slots=True)
class GraphCut:
    name: str
    left: frozenset[str]          # owner names
    right: frozenset[str]
    crossing: tuple[str, ...]     # connector names crossing the cut, sorted
    bond_budget: BondBudget       # SYMBOLIC upper bound; never a number
```

### 4.1 Legality predicates — checkable now, at L1

| # | Rule | Why | Refusal |
|---|---|---|---|
| C1 | A composite series group may not be split | evaluating a series of two compliances as two independent springs computes the wrong stiffness; the group must dispatch together (`ALEPH-PORT-1501` §3) | `CUT_SPLITS_COMPOSITE` |
| C2 | An `adjoint_required` connector may be split only if the plan declares an exact `+f/-f` architectural layer across that bond | a one-way force is a momentum source, and a momentum source inside a closed assembly makes force closure unfalsifiable | `CUT_SPLITS_ADJOINT_PAIR` |
| C3 | A cut may not separate an `internal_to` connector from its host component | an `I`-scope object is not an independent endpoint | `CUT_SPLITS_INTERNAL` |
| C4 | Every crossing connector must be `implemented`, or the cut is declared over fiction | the ledger must not let a declared-but-unevaluated connector pass | `CUT_CROSSES_UNIMPLEMENTED` |
| C5 | Cut size counts a composite group once, not once per member | otherwise the cost model prefers cuts that split composites | — (definition) |

### 4.2 What is *not* decidable now, and must refuse

Bond dimension is the rank of the interaction across the cut. That is a measured quantity. At L1 a cut
carries a `BondBudget` symbol and an ordering — cut A is not wider than cut B when A's crossing set is
a subset of B's — and any call that asks for a number raises `BondDimensionUnmeasured`. This is the
place where a project quietly invents `rank=16` and never revisits it.

---

## 5. Axis order is derived, so that it is reproducible

§3.2 requires axis order to follow from "the owner graph, connector cut sizes, topology sector, event
locality, and query family." The derivation is a total function, and its determinism is what makes a
plan hashable:

1. Take the **query family** (which owners and observables the plan must serve) as the root set.
2. Breadth-first traverse `G` from the root set over `E`-scope vertices, so distance is *connector
   count*, not spatial distance.
3. Order owners by (BFS depth, **coupling to the already-ordered set, descending**, stable id) —
   `TieBreak.STRONG_COUPLING_ADJACENT`, PI-decided 2026-07-30 ~14:56. The bond dimension across a cut
   is set by the correlation crossing it, so two owners that share several connectors must sit next
   to each other; separated, they force *every* cut between them to carry that coupling. Weakly
   coupled owners therefore fall to the ends of the chain, where a truncation costs least.

   **This reverses what an earlier draft of this step said**, and the reversal is a correction rather
   than a preference: the old rule sorted cut size *ascending*, which put weakly coupled owners next
   to the root — the opposite of the rationale printed beside it. Rule and reason contradicted each
   other and the contradiction survived until the V1 example was worked by hand. The alternative
   remains available as `TieBreak.WEAK_COUPLING_ADJACENT`, and the choice is part of the plan hash,
   so a future measurement can overturn it visibly.
4. Within an owner, order blocks in the fixed sequence `MESH_SURFACE, OWNER_PARTICLES,
   CONTINUUM_FIELD, CONNECTOR, EVENT_STREAM, OBSERVATION`; within a block, order indices by kind in the
   fixed sequence `SECTOR, ENTITY, SPATIAL, FIELD, MODE, EVENT, ENDPOINT, TIME`, then by name.
5. Sector-dependent axes come first within their owner, so a sector jump re-plans a prefix rather than
   the whole chain.

Ties are broken by `stable_id`, never by dictionary order. Two runs of step 1–5 over the same frozen
schema and the same query family produce the same order, and `plan_hash` is SHA-256 over the canonical
JSON of the ordered declaration, domain-separated by a contract id — the same discipline the schema
hash already uses.

---

## 6. Factorization families and where each is admissible

| Family | Admissible blocks | Note |
|---|---|---|
| `DENSE` | all | the always-legal fallback; never needs qualification |
| `SPARSE` | particles, connectors, events | |
| `SVD_POD` | mesh, field | a linear-subspace baseline that TT must beat to be worth its complexity |
| `TT_MPS` | field, mesh, particles (**local only**) | §10: local TT/MPS on particles only after measured qualification |
| `TTN` | mesh, field | hierarchical cuts |
| `GRAPH_TN` | mesh, connectors | when the cut structure is not a chain or a tree |
| `MPO` | operators, not states | an operator representation; declared here so an operator plan is expressible |
| `SPECTRAL_BASIS` | mesh, field | the membrane mode spectrum already used by the observation lane |
| `NEURAL_OPERATOR` | field | learned; L2 only |
| `EQUIVARIANT_PARTICLE` | particles | learned; L2 only |
| `MARKED_PROCESS` | events | Hawkes/marked encoder; L2 only |
| `TYPED_EDGE` | connectors | must carry the exact `+f/-f` or energy-consistent layer |
| `PARTICLE_ENSEMBLE` | all | the fallback form when a channel exceeds tolerance |

Selecting a family is L1 and free. Claiming a family *works* is L2 and sealed. A plan that names a
learned family is legal to write and refuses to instantiate.

---

## 7. The error vector — eight channels, never summed

§3.2 lists the channels and requires them separate. Summing them into one scalar is prohibited,
because the sum cannot say which invariant broke, and a fallback trigger has to know which invariant
broke to choose a fallback.

| Channel | Measures | Fallback if exceeded |
|---|---|---|
| `BASIS` | the basis cannot represent the state at all | preserve native block |
| `RANK_TRUNCATION` | discarded singular/bond weight | raise budget, else native |
| `MASS` | violated conservation of a conserved quantity | native block; never rescale silently |
| `NEGATIVITY` | a non-negative quantity went negative | native block |
| `OBSERVATION` | error in the predicted observation, not in the state | re-plan for the query family |
| `DYNAMICS` | error accumulated over stepping, not at one instant | shorten re-plan interval, else native |
| `TOPOLOGY_TAIL` | error concentrated near a sector jump | sector-local native block |
| `BYTE_DECODE` | serialization round-trip is lossy | hard failure; never a tolerance |

Each channel owns a `ToleranceSlot` that is **empty** at L1. Reading an unfilled slot raises
`ToleranceUnset`. Filling one requires a decision record — a threshold chosen from the runs that
motivated it is a gate fitted to its data, which the project has already refused once
(`ALEPH-DQ-104`).

**Fallback is recorded, never silent.** A plan that falls back carries a `FallbackRecord` with the
channel, the trigger, and the replacement form, and the record travels with any artifact the plan
produced.

---

## 8. MTG-PN — the neural side

### 8.1 The three layers are the three reader lanes — PI-decided 2026-07-30 ~14:56

The instruction's "three layers" means **the three reader lanes**, `MTG_PN_SIM`, `MTG_PN_EXP` and
`MTG_PN_HYBRID`, which is §10.1's own table. The question is settled and this document no longer
carries two readings of it.

**The encoder hierarchy is not the other reading and is not deleted.** H1–H3 is the stack all three
lanes share — one encoder, three provenance-separated ways of feeding and trusting it — so it is a
component of the answer rather than a competitor to it. It stays because §10.1's architecture needs
it and `aleph/learn/architecture.py` implements it. What is gone is the framing that treated the two
as rival interpretations of the same phrase.

**(a) Three reader lanes — the horizontal split, by evidence provenance.**

| Lane | Inputs | Permitted output | Ceiling | Gate |
|---|---|---|---|---|
| `MTG_PN_SIM` | accepted native states, interventions, traces, operators, representation errors | simulation posterior, regime, proposal, acceleration; outside the trust region it **requests a native run** | may not exceed the authority of the traces it trained on | an accepted native trace digest must exist — **unsatisfiable today** |
| `MTG_PN_EXP` | images, force/time series, omics, identity/context, protocol manifests | empirical state prior/likelihood/predictive observation; microscopic mechanism stays broad where unidentifiable | **permanent non-native ceiling** — it can never certify a native claim | dataset identity, protocol, calibration and split-leakage resolved |
| `MTG_PN_HYBRID` | both lanes plus operator-aligned paired or unpaired evidence | joint posterior and counterfactual response | bounded by the weaker of its two parents; must carry native ancestry | both parent gates, plus a domain-shift/OOD refusal path |

**(b) The encoder hierarchy — the vertical stack, shared by all three lanes.**

```
H1  node / head / segment                     leaf entities, equivariant
H2  filament / mesh patch / event neighborhood   local pooling, geometry-aware
H3  component token + local field representation  one token per census owner
--- below the line: not hierarchy levels ---
C   representation-qualified tensor core      factorized per the L1 plan   [PI-added 14:56]
G   typed connector graph + accepted-state identity
L   observation likelihood factors
P   posterior / predictive observation / intervention response
V   uncertainty + OOD verdict + native-query request
```

**C is where the tensor language enters the reader**, and it was missing until the PI restored it
(finding 21). It sits between H3 and G because that is where the state is owner-shaped and not yet
edge-shaped — a factorization is defined over owner state, and the graph then passes messages between
cores rather than between raw tokens. It is not a fourth pooling level: it does not summarise its
input into coarser tokens, it factorizes the state those tokens already describe, and **the structure
comes from the plan rather than from the encoder** — family per block, axis order, cuts and bond
budget all arrive from an L1 `FactorizationPlan`. *Qualified* is the load-bearing half of the name: a
core may be used only at a rank whose eight error channels were measured against an accepted native
trace, so the stage is unreachable until R5. It consumes no `BOND` index — a bond is what it
introduces, and a core reading bonds would be reading its own output.

H1–H3 are declared as `EncoderStage` records: input block, output token space, equivariance group,
pooling rule, and which representation-plan indices they consume. The plan from §5 is what tells the
encoder its axis order — this is the join between the two halves of the document, and it is why the
language had to be finished first.

### 8.2 One output schema for all three lanes

All lanes emit the same six types, so that a consumer never has to know which lane answered:

`PosteriorSamples`, `PredictiveObservation`, `InterventionResponse`, `UncertaintyBreakdown`,
`OODVerdict`, `NativeQueryRequest`.

Two rules with teeth:

- **Missing modality uses an explicit mask** and a calibrated product-of-experts fusion. Absence is
  never encoded as a zero measurement — the same rule as index invariant I6, in a different layer.
- **`NativeQueryRequest` is a first-class return value, not an exception.** A reader that cannot
  answer within its trust region succeeds by asking for a simulation. This is the mechanism that keeps
  a surrogate from silently extrapolating, and it is why the type sits in the shared schema.

Every output carries `schema_hash`, `context_hash`, `accepted_state_digest`, `lane`, `authority_ceiling`
and `representation_plan_hash`. An output whose plan hash does not resolve is not interpretable.

### 8.3 Losses — the exact and the optimized are different kinds of thing

| Architectural / exact — **not** loss terms | Optimized / measured — loss terms |
|---|---|
| one owner per entity | state reconstruction |
| typed endpoints | event likelihood |
| accepted-state digest | observation likelihood |
| one clock | neural posterior estimation |
| missing-modality masks | physics residuals |
| evidence ceiling | TN approximation error |
| fallback / refusal | calibration |
| passive action–reaction (or energy-consistent parameterization) | force/work discrepancy |
| separate active-power ledger | posterior coverage, intervention error |

The left column is enforced by the architecture and must be *impossible* to violate, not penalized
when violated. A soft penalty on action–reaction is a statement that momentum conservation is
negotiable at a price, and the registry already refuses that in the connector layer. The loss registry
therefore types the two columns differently: an `ExactConstraint` has no weight field at all, and
attempting to give one a weight is a refusal rather than a configuration.

### 8.4 Evaluation — the split rule is a hard constraint

Hold out **entire** intervention, regime and archetype for simulation; **entire** donor, lab, batch and
perturbation for experiment. **Random neighbouring-frame splits are prohibited** — neighbouring frames
are the same state observed twice, so a random split measures interpolation and reports it as
generalization. This is a checkable predicate over a split declaration, and it is written now so that
the first training run cannot be the moment it gets negotiated.

Reported for every evaluation: posterior coverage, calibration, support, OOD refusal rate, uncertainty
source decomposition, representation error (by channel, unsummed), held-out intervention response.

Ablations, each a named record: Dense vs TT; SVD/POD vs TT; event stream on/off; field encoder on/off;
connector graph on/off; simulation prior on/off; real-data likelihood on/off; native fallback on/off.

A model failure is classified into exactly one of five kinds, and the classification is part of the
result rather than a discussion afterwards: `WRONG_PRIOR`, `NON_IDENTIFIABILITY`,
`INSUFFICIENT_OBSERVATION`, `REPRESENTATION_FAILURE`, `MISSING_MECHANISM_OR_CONNECTOR`.

---

## 9. Module map and lane ownership

```
aleph/represent/
  __init__.py        seal, DECLARATION_ONLY_MODULES, ScopeFirewallError, RepresentationSealed
  index.py           Index, IndexKind, Extent, BasisKind, contraction legality       [L23]
  blocks.py          StateBlock, block derivation, H/B/X slots, census binding       [L23]
  cuts.py            GraphCut, BondBudget, legality predicates C1–C5                 [L24]
  families.py        FactorizationFamily, admissibility, cost ordering               [L24]
  plan.py            TensorDeclaration, FactorizationPlan, axis order, plan_hash     [L25]
  errors.py          ErrorChannel, ToleranceSlot, FallbackRecord                     [L26]
  qualification.py   the L2 protocol, every entry point refusing                     [L26]

aleph/learn/
  __init__.py        seal, DECLARATION_ONLY_MODULES
  outputs.py         the six shared output types + provenance fields                 [L27]
  lanes.py           the three lanes, admissibility, ceilings, gates                 [L27]
  architecture.py    EncoderStage H1–H3, connector graph layer, exact layers         [L28]
  losses.py          ExactConstraint vs MeasuredLoss registry                        [L28]
  evaluation.py      split rules, ablation matrix, failure taxonomy                  [L29]
```

Every lane writes tests under `tests/represent/` or `tests/learn/` mirroring its module, and every
lane's tests must include at least one **counterexample test**: a declaration that must be refused,
with the refusal type asserted. A rule with no counterexample test is a comment.

## 10. What this specification deliberately does not decide

- No rank, width, token count, or threshold. §10.3 preregisters those after pilot measurements.
- No choice between conditional flow and another calibrated density estimator for the hybrid posterior.
- No adaptive-campaign acquisition function. Expected decision-relevant information gain is named as
  the criterion; the estimator is not chosen.
- No claim that any family will pass qualification for any block. Every admissibility entry above says
  *may be planned*, never *works*.
</content>
</invoke>

---

## 11. Defects in this specification, found by implementing it once

Two lanes implemented §2–§4 against the real census before the session was stopped. Their code is
archived and unratified; what follows is worth more than the code — the places where building it
showed the specification to be wrong, ambiguous, or missing something. **These are unfixed.** Anyone
implementing this must resolve them rather than route around them.

### 11.1 The index algebra (§2)

1. **`Extent.literal_dimension(3)` is a bad name and invites the thing it forbids.** `Extent.spatial(n)`
   says what it is.
2. **There is no constructor for an endpoint arity.** §2 says `ENDPOINT` has "extent 2 for a binary
   connector" but names only `of`, `literal_dimension` and `unbounded`. Reusing a *spatial* extent for
   a count conflates a space with a cardinality. `Extent.arity(n)` is needed.
3. **Two `unbounded` extents are not contractible, and the spec does not say so.** I4 requires the
   extent symbol to agree, but `unbounded()` carries no stream identity, so two append-only event axes
   agree only if they are the same stream — and nothing records which stream.
4. **I2 and I3 overlap on `B` scope and the spec does not order them.** Both refuse a `B`-scope state
   axis; `CONTEXT_AS_STATE` is the more actionable refusal, so I3 must run first, or the message says
   what is wrong without saying what to do.
5. **`H` scope is stronger than I2 states.** §3 says a homogenized law "owns no index", which also
   forbids it a *conditioning* axis — something "cannot be state" does not cover.
6. **I6's "absent" is an instance fact, not a schema fact.** A frozen schema cannot know that a
   particular cell grew no filopodia. Absence must be supplied by the caller alongside the schema; the
   spec reads as though the schema knows.
7. **"A unit symbol resolvable in the unit registry" is a comment, not an invariant.** It is missing
   from the I-table, and `aleph.units` is unimportable from a declaration-only module anyway because
   it pulls in numpy transitively. `FrozenStateSchema.unit_declarations` is the pure path, and the
   spec never connects the two.

### 11.2 Block derivation (§3)

8. **"Block assignment is derived" is under-determined.** Deriving a block from `owned_state` name
   text is keyword matching, and rule *order* carries real weight — `surface_traction_pn` matches both
   a connector rule and a mesh rule. The honest design reports what it cannot classify instead of
   guessing, which makes the override manifest load-bearing rather than an escape hatch.
9. **The three census manifests spell their aggregate tuples differently** — `CENSUS_ENVIRONMENT_SURFACE`,
   `FRAME_FLUID_NUCLEUS_COMPONENTS`, `ACTOMYOSIN_LOAD_PATH` — so binding by module attribute name is
   fragile. Bind through the frozen schema's queries instead.
10. **`H`/`B`/`X` records cannot be manufactured from census text.** A `HomogenizedLawSlot` requires a
    declared inadequacy and a re-entry condition; synthesising those from whichever contract field
    happens to be non-empty fabricates exactly the content the record exists to require. They must be
    authored, or refused.

### 11.3 Axis order and cuts (§4–§5), from the V1 worked example

11. **DECIDED 2026-07-30 ~14:56 — strong coupling stays adjacent.** The PI chose the cortex side, and
    it is also the standard chain-ordering heuristic, so §5 step 3 is reversed and the old rule
    survives as `TieBreak.WEAK_COUPLING_ADJACENT`. The finding turned out to be sharper than first
    written: §5's rule and its own printed rationale contradicted each other. Original text: **Cut
    size and query relevance disagree at equal BFS depth, and §5 chose without noticing.** In
    V1, `cortex` is reached by two connectors and `cytosol` by one, so cut-size ordering puts the
    weakly coupled fluid *closer* to the membrane than the shell that carries its load — while the
    query is a membrane spectrum, for which `cortex` is the relevant neighbour. §5 says cut size wins.
    That may be right for truncation cost and is not obviously right for a query-served plan.
    **Open question for the PI.**
12. **Every V1 cut refuses today**, because rule C4 requires each crossing connector to be
    `implemented` and the registry's default is `False`. This is correct behaviour, and it means the
    representation layer's first green integration test arrives with the runtime, not with the
    language.

### 11.4 The error vector and qualification (§6–§7), found while building them

13. **The five qualification targets are not in this document.** Analytic limits, local projections,
    worst ROIs, events, and conservation/error budgets exist only in
    `docs/manuscripts/extracted/SECTION_10_tensor_and_mtgpn.txt`. §9 gives `qualification.py` one
    line. A manuscript extract is a reference, never authority, so the targets must be absorbed here
    before qualification can cite anything.
14. **§6 and §7 contradict each other on fallback forms.** §6 calls `PARTICLE_ENSEMBLE` "the fallback
    form when a channel exceeds tolerance"; §7's table never mentions it and names five per-channel
    forms instead. Either it is a sixth form or it is what `PRESERVE_NATIVE_BLOCK` means for a
    particle block. **Needs a ruling.**
15. **`NEGATIVITY` has no forbidden repair and should have one.** §7 names "never rescale silently"
    for `MASS` alone. Clipping a negative to zero destroys the same evidence and is at least as
    tempting.
16. **`BYTE_DECODE` has no stated trigger.** A channel with no tolerance is compared against nothing.
    The implied trigger is "non-zero at all" and no document says it.
17. **CLOSED 2026-07-30 15:05, PI-delegated.** A tolerance now has one form per channel, fixed in
    `errors.ToleranceForm` and `TOLERANCE_FORM`:
    `MIXED_ABSOLUTE_RELATIVE` — `|error| <= atol + rtol*|reference|` — for the five approximation
    channels, because a purely relative bound explodes where the reference passes through zero and a
    purely absolute one changes meaning under a unit change or a mesh refinement;
    `ROUNDOFF_BUDGET` for `MASS` and `NEGATIVITY`, which are bug detectors rather than accuracy
    knobs, so their bound is a statement about floating-point accumulation and a percent-level
    tolerance there would hide the defect it exists to catch; and `EXACT` for `BYTE_DECODE`, which
    owns no slot at all. A mixed bound must name a **reference scale symbol** resolved from the plan
    — never the measured value, since a bound taken against what it is bounding always passes. No
    number was set. Original text: **Nothing decides what a tolerance *is*** — absolute, relative,
    per-sector, or a pair of bounds.
    §7 requires a decision record and is silent on the value's shape, so two lanes filling slots under
    different conventions would later compare incomparable numbers. **One decision fixes this.**
18. **`ANALYTIC_LIMIT` does not structurally need an accepted native trace.** Its ground truth is
    closed-form, which is also the only case that can fail because the *reference* is wrong. So one
    fifth of the qualification protocol is in principle runnable before R5. Starting it early is a
    seal change and a decision record, not a code change.
19. **Enum value casing is not settled across the declaration layer.** `IndexKind`, `BasisKind` and
    `StateBlock` use lowercase values; `ErrorChannel`, `LaneId` and `FactorizationFamily` names are
    uppercase. It has already produced two false reconciliation failures, each absorbed by
    case-folding at the boundary. It should be settled once rather than absorbed a third time.

### 11.5 Two defects that are structural, not editorial

20. **§10's native-form column is not disjoint, so no faithful rule table can be single-valued.**
    "Forces" appears under both `OWNER_PARTICLES` and `CONNECTOR`; "binding/topology state" under both
    `OWNER_PARTICLES` and `EVENT_STREAM`. Finding 8 called this an ordering problem; it is not. Order
    was the symptom and the overlapping vocabulary is the cause, so a derivation that returns one
    block for every name is *necessarily* lying about one of these.

    Measured against the real seed census — 36 components, 104 `owned_state` names — the honest
    derivation gives **54 derived, 21 ambiguous, 29 unclassified: 48% refuses**, and 17 of the 21
    ambiguities are the `EVENT_STREAM`/`OWNER_PARTICLES` overlap alone. `surface_traction_pn` comes
    out ambiguous, reproducing this document's own worked example.

    That is the correct behaviour and it is also a workload: roughly fifty owner-state names need an
    authored adjudication. **No lane may write those** — each is a modelling decision about who owns
    a quantity — so the override manifest ships empty and the list is the question.

21. **CLOSED 2026-07-30 ~14:56 by PI decision — the core goes between H3 and the connector graph.**
    Implemented as `StageId.C_REPRESENTATION_QUALIFIED_TENSOR_CORE`; see §8.1(b) above. Original
    text: **§8.1(b) drops the representation-qualified tensor core.** Manuscript §10.1's prose lists nine
    components including "representation-qualified tensor core", but its own arrow chain omits it and
    §8.1(b) copied the chain. So the encoder stack has no stage for the thing the entire tensor
    language exists to make possible, and the implementation followed the spec rather than patching
    it. **This needs a decision:** either the core is a stage between H3 and the connector graph, or
    the language feeds the readers some other way and §8 should say how.

### 11.6 Smaller findings from the same pass

22. **§6's admissibility column mixes domains.** `MPO`'s entry is "operators, not states", which is
    not a set of state blocks. Modelled as an explicit "not over states" rather than an empty set,
    because an empty set reads as "nothing matched".
23. **`OBSERVATION` has no census source.** Its instances are the operators in `aleph/observe/`, so
    deriving it from a component's owned state is refused — which leaves only `DENSE` and
    `PARTICLE_ENSEMBLE` admissible for it. That may not be intended.
24. **§6's "local only" on TT/MPS over particles is a footnote with no mechanism.** Made load-bearing:
    admissibility refuses unless the caller acknowledges the restriction.
25. **§8.3's last measured cell packs two quantities** — posterior coverage and intervention error.
    One weight over a coverage statistic and a response error would be one decision record standing in
    for two measurements, so they are registered separately and the table's row-pairing is cosmetic.
26. **§8.3 drops "predictive checks"** from the manuscript's measured column, unless it is
    `calibration` under another name. Representation error is correctly absent — it is the §7
    eight-channel vector, an L2 quantity, never a loss.
27. ~~**A `TrainingPlan` validates only that its stage and constraint lists are non-empty.**~~
    **CLOSED 14:40.** `TrainingPlan.__post_init__` now delegates to
    `architecture.assert_plan_matches_architecture`, so a plan naming stages that do not exist is
    refused at construction rather than at a check a caller might skip. Original text: The
    registries to check them against now exist (`ARCHITECTURE_STAGE_NAMES`, `EXACT_CONSTRAINT_NAMES`,
    and `assert_plan_matches_architecture`), but nothing forces a caller to run them, so a plan naming
    meaningless stages still hashes cleanly. Closing it is a call from `TrainingPlan.__post_init__`
    behind a local import — deliberately not done at handoff time, because it changes what an existing
    valid plan is.
