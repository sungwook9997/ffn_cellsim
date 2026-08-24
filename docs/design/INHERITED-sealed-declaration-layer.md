# INHERITED — the sealed declaration layer, as a specification

`Project_Aleph/aleph/learn/` (6,163 lines) and `represent/` (6,237) are **computation-sealed by PI
ratification 2026-07-30**: no array or autodiff import, no float literal, no integer outside
{0,1,2,3,64}, and every compute entry point raises a typed refusal naming gate R5. Verified on
2026-08-09: **15 files, 12,400 lines, zero `numpy` / `torch` imports.**

So there is no code to port — there is a specification, and decision **A4** says it travels as one.
This file is that specification, extracted: every module's docstring and every public symbol, with
the arithmetic that does not exist left out because it does not exist.

> ⚠ **`AGENT-PROPOSED` and unratified.** That is a reason it may not be used to **reject** an
> implementation. It is not a reason to leave it in a tree that is becoming an archive. Where this
> tree already has something on the same axis, this tree's is the one that runs — see
> `docs/decisions/PROPOSAL-the-merge-and-the-rename.md` §1a, where `virtual_cell/observation_operator.py`
> turned out to be the stronger of the two.

Source: `Project_Aleph` at `9ef077ba`, `aleph/{learn,represent}/`.


---

# `learn/`

## `learn/architecture.py` — 963 lines

The MTG-PN stack of §8.1(b) as declarations: what each stage reads, emits and preserves.

An architecture module that contains layers contains a model, and a model is sealed. What this
module contains instead is the *specification* each layer would have to satisfy: for every stage,
which state blocks it reads, which representation-plan index kinds it consumes, what token space it
emits, which symmetry it preserves, and how it pools. That is enough to check the stack for the
mistakes that are expensive later and free now — a stage reading an index its block does not own, a
stage in the equivariant part of the hierarchy that quietly fixes a frame, a state block that no
stage consumes and which therefore silently never enters the model.

Seven stages, not three
-----------------------
§8.1(b) draws H1, H2, H3 as a hierarchy and then draws four more arrows. The four are not hierarchy
levels and are typed as a different `StageKind`, because conflating them is how "the encoder" ends up
meaning both "the thing that pools segments into filaments" and "the thing that decides whether to
refuse". H1–H3 pool; G relates; L brings measurements in; P answers; V says how much to believe the
answer and when not to.

The exact layers
----------------
The left column of §8.3 is enforced here, by construction, and `losses.py` is where it is declared.
`EXACT_ARCHITECTURAL_LAYERS` names, for each constraint in that column, the stage that makes it
impossible to violate and the mechanism that does so. An import-time check asserts the cover is
total: an exact constraint with no enforcing layer is a rule that was written down and not built.

These layers carry no weight field, and `ExactArchitecturalLayer.weigh` refuses with the same
`ExactConstraintNotWeightable` the registry uses. There is one refusal type for the category error,
not two, so a caller catching it catches it wherever they reached for it.

Where the representation plan joins
-----------------------------------
Every stage declares `consumes_index_kinds`, and index kinds are `aleph.represent.index.IndexKind` —
the same objects the L1 plan is written in, imported rather than mirrored. The axis order a stage
reads is derived in §5 from the owner graph and the cut sizes; a stage cannot know it without the
plan, which is why `build` and `instantiate` take a `representation_plan_hash` and why an output
whose plan hash does not resolve is uninterpretable (§8.2).

Nothing here instantiates. `build`, `instantiate` and `forward` carry full signatures and raise
`TrainingSealed`. Per `ALEPH-PD-003` the signature is the design and the refusal is the evidence the
design was not quietly executed; an absent function proves nothing, because it looks the same whether
the work was deferred or forgotten.

**Public surface:** `ArchitectureError`, `IndexKindNotInBlock`, `BondIndexNotConsumable`, `EquivarianceUndeclared`, `PoolingRuleMissing`, `StageChainBroken`, `BlockNotConsumed`, `UnknownEncoderStage`, `ConstraintUnenforced`, `StageKind`, `StageId`, `TokenSpace`, `TokenCardinality`, `EquivarianceGroup`, `PoolingRule`, `EncoderStage`, `ExactArchitecturalLayer`, `stage`, `assert_declared_stages`, `assert_plan_matches_architecture`, `build`, `instantiate`, `forward`

## `learn/evaluation.py` — 1915 lines

How an MTG-PN reader is evaluated — the split rule, the reported-metric contract, the ablation
matrix and the failure taxonomy, all written before the first training run.

`ALEPH-TN-LANGUAGE-SPEC` section 8.4 is the specification this module implements; manuscript section
10.3 is its reference source. Nothing here computes, and nothing here can: every entry point that
would *run* an evaluation raises :class:`~aleph.learn.TrainingSealed`, and every number the design
will eventually need is an :class:`EmptySlot` rather than a default.

Three things live here, and each exists because of a specific way an evaluation goes wrong.

**The split rule, as a predicate rather than a habit.** Hold out *entire* intervention, regime and
archetype for the simulation lane; *entire* donor, lab, batch and perturbation for the experimental
lane. Random neighbouring-frame splits are prohibited, and the reason is not statistical fussiness:
neighbouring frames are the same state observed twice, so a random split puts a state in the training
set and its near-copy in the held-out set. The resulting number measures interpolation and gets
reported as generalization. The check is written now, before a training run exists, because the
moment somebody needs a number is exactly the moment this rule becomes negotiable.

**The reported-metric contract, as named empty slots.** Seven metrics are reported for every
evaluation, each with a stated meaning and a stated unit-or-dimensionless status. All of them are
empty, and reading an unfilled one is a refusal rather than a zero. Two of them have no scalar form
at all: representation error is per channel and uncertainty is per source, and asking either for a
total raises :class:`ErrorChannelsNotSummable`. Section 7 of the spec gives the reason — a total
cannot say which invariant broke, and a fallback trigger has to know which invariant broke in order
to pick a fallback.

**The ablation matrix and the failure taxonomy, as records.** Eight ablations, each naming what it
removes and which claim would be unsupported without it. Five failure kinds, exhaustive, and the
classification is a *field of the result* rather than a paragraph written afterwards — a failed
evaluation that has not been classified cannot be constructed.

What is reused from elsewhere in Aleph, and what is deliberately not:

* :data:`aleph.evidence.claim.CLAIM_ID_PATTERN` is imported and used to check the form of the claim
  identifier each ablation protects. The claim itself need not exist yet; the pattern check is the
  same "check form, not existence" discipline that
  :func:`aleph.evidence.ladder.is_artifact_digest` already uses.
* :mod:`aleph.infer.refusal`'s vocabulary is referenced *by wire label* in the failure taxonomy, so
  that a classification of ``NON_IDENTIFIABILITY`` here and an ``Unidentifiable`` refusal there are
  visibly the same finding. :func:`assert_refusal_kinds_resolve` checks the labels against the real
  module on demand; the module is not imported at module scope, so this file stays importable in
  isolation.
* :mod:`aleph.infer.identifiability` is **not** imported, and cannot be: it depends on ``numpy``, and
  `ALEPH-PD-003` forbids an array library anywhere under ``aleph/learn/``. The link is by name only.
  A ``NON_IDENTIFIABILITY`` classification is expected to cite a degeneracy report produced there,
  which is why :class:`FailureClassification` requires stated evidence and will not accept a bare
  enum member.
* :mod:`aleph.learn.outputs` and :mod:`aleph.learn.lanes` are imported lazily inside the
  ``assert_*`` helpers only, so this module does not depend on their landing order.

**Public surface:** `EvaluationDeclarationError`, `IllegalSplit`, `MetricContractError`, `SlotUnfilled`, `ErrorChannelsNotSummable`, `UnclassifiedFailure`, `LaneName`, `HoldoutAxis`, `SplitMethod`, `SplitRefusal`, `SplitViolation`, `EvaluationUnit`, `SplitDeclaration`, `check_split`, `require_legal_split`, `EmptySlot`, `MetricSlot`, `MetricDimensionality`, `MetricShape`, `MetricDeclaration`, `RepresentationErrorChannel`, `RepresentationErrorReport`, `empty_representation_error_report`, `metric_contract`, `FailureKind`, `FailureKindDeclaration`, `FailureClassification`, `classify_failure`, `Ablation`, `AblationDeclaration`, `EvaluationOutcome`, `EvaluationResult`, `assert_lane_names_match_learn_lanes`, `assert_error_channels_match_representation_layer`, `assert_uncertainty_sources_resolve`, `assert_refusal_kinds_resolve`, `run_evaluation`, `run_ablation`, `run_ablation_matrix`, `measure_representation_error`, `measure_metric`

## `learn/lanes.py` — 899 lines

The three MTG-PN reader lanes of §8.1(a): what each may read, emit, claim — and why none may run.

A lane is a *split by evidence provenance*, not by architecture. All three share the encoder stack of
§8.1(b) and all three emit the identical six-type schema of `outputs.py`, so what actually
distinguishes them is the only thing a consumer cannot recover from the numbers: where the evidence
came from, and therefore how strong a claim the numbers are allowed to be.

This module declares that. It trains nothing, instantiates nothing, and has no path to doing either.

The gates are unsatisfiable today, by construction rather than by policy
------------------------------------------------------------------------
`MTG_PN_SIM` and `MTG_PN_HYBRID` require an accepted native trace digest. The registry of accepted
native trace digests is `ACCEPTED_NATIVE_TRACE_DIGESTS`, and it is empty — not because nobody has
filled it in yet, but because no accepted-step transaction exists that could prove a trajectory was
accepted rather than force-accepted. So `LaneGate.evaluate` cannot return a satisfied verdict for
those two lanes for *any* evidence a caller can construct: handing it a well-formed digest still fails
the membership test against an empty set. That is what "unsatisfiable by construction" means here, and
`test_lanes.py` asserts it against fabricated evidence rather than against the absence of evidence.

`MTG_PN_EXP`'s gate is different in kind and worth keeping distinct: dataset identity, protocol,
calibration and split leakage are *unresolved*, not impossible. Its status is `UNRESOLVED`. What is
permanent for `MTG_PN_EXP` is not its gate but its ceiling.

The ceiling algebra has teeth
-----------------------------
`MTG_PN_EXP` sits permanently at `NON_NATIVE`: no quantity of experimental evidence certifies a claim
about an accepted native state, because the experiment never saw one. `MTG_PN_SIM` sits at
`TRACE_BOUNDED` — it may not exceed the authority of the traces it read. §8.1 says the hybrid is
"bounded by the weaker of its two parents", so the hybrid ceiling is *computed* here with
`AuthorityCeiling.meet` and never typed in, and the arithmetic says `NON_NATIVE`.

`compose_authority` and `compose_hybrid_provenance` are where that becomes load-bearing: composing an
`MTG_PN_EXP` result into a joint result cannot produce a native-authority claim, and asking for one
raises `AuthorityEscalationError` instead of being quietly clamped. Clamping would be the friendlier
behaviour and the wrong one — a caller who asked for native authority has a belief about their
evidence that is false, and silently downgrading the label leaves the belief in place.

Native ancestry is separately enforced: §8.1 requires the hybrid to carry it, so composing only
experimental contributions is a `ProvenanceIncomplete` refusal rather than a hybrid result with the
sentinel in its accepted-state field.

Entry points exist and refuse
-----------------------------
`instantiate`, `fit`, `train` and `predict` carry full signatures and raise `TrainingSealed` naming
what they would need. Per `ALEPH-PD-003`: the signature is the design, and the refusal is the evidence
the design was not quietly executed. An absent method proves nothing — it looks the same whether the
work was deferred or forgotten.

**Public surface:** `InadmissibleModalityError`, `LaneGateUnsatisfiedError`, `InadmissibleModality`, `LaneGateUnsatisfied`, `GateRequirement`, `GateStatus`, `GateVerdict`, `GateEvidence`, `LaneGate`, `ClaimKind`, `LaneReader`, `LaneSpec`, `lane_spec`, `compose_authority`, `compose_hybrid_provenance`

## `learn/losses.py` — 663 lines

The §8.3 registry: the exact column and the measured column, typed so they cannot be confused.

The table in §8.3 has two columns and it is tempting to read it as one list with a flag. It is not.
The two columns are different kinds of thing, and the difference is the whole content of the table.

A **measured loss** is a number the optimizer trades against other numbers. It needs a weight, the
weight is a measured quantity, and §10.3 preregisters it only after a pilot. So a measured loss
carries a `WeightSlot` — the same slot `training.py` already defined — and the slot is empty.

An **exact constraint** is not a number at all. "One owner per entity" is not a quantity that can be
partially satisfied; it is a property the architecture either has or does not have. Giving it a
weight is not a stricter or looser setting of the same dial, it is a category error with a physical
consequence: a soft penalty on passive action–reaction is a statement that momentum conservation is
negotiable at a price, and once it is negotiable the force closure of a connector assembly stops
being falsifiable. The connector registry already refuses that (`ALEPH-PORT-1501` §3), and this
module refuses it in the same way — by making the weight *unrepresentable*.

So `ExactConstraint` has no weight field. Not an optional one, not one defaulting to infinity: none.
`dataclasses.replace(constraint, weight=...)` is a `TypeError` before any check of ours runs, and the
paths that could attach a weight by name (`weight_slot_for`, `assert_weight_slots_are_measured`)
raise `ExactConstraintNotWeightable` and say why.

Where the architecture comes in
-------------------------------
This module says *what* is exact. `architecture.py` says *where each exact constraint is enforced by
construction*, and asserts at import that every constraint in this registry has an enforcing layer.
An exact constraint with no enforcing layer would be a rule stated and then not implemented, which is
indistinguishable from a comment.

Nothing here computes a loss. There is no tensor, no reduction and no weight value anywhere in the
module; a loss registry that could evaluate a loss would be the training layer, and that is sealed.

**Public surface:** `LossRegistryError`, `ExactConstraintNotWeightable`, `UnknownLossName`, `ColumnConfusion`, `IncompleteRegistry`, `LossColumn`, `ExactConstraintId`, `MeasuredLossId`, `ExactConstraint`, `MeasuredLoss`, `LossRegistry`, `weight_slots`, `weight_slot_for`, `assert_declared_exact_constraints`, `assert_weight_slots_are_measured`

## `learn/outputs.py` — 1222 lines

The six shared output types every MTG-PN lane emits, and the provenance that makes them readable.

Specification §8.2 fixes one output schema for all three reader lanes — `PosteriorSamples`,
`PredictiveObservation`, `InterventionResponse`, `UncertaintyBreakdown`, `OODVerdict`,
`NativeQueryRequest` — so that a consumer never has to know which lane answered. This module is that
schema and nothing else. It declares *what a result carries and what makes it interpretable*; it
computes no result, and under `ALEPH-PD-003` it may not import an array library to do so.

Why the provenance block is part of the type
--------------------------------------------
An output that cannot be traced back to a schema, a context, an accepted state, a lane, a ceiling and
a representation plan is a number with a story attached, and stories survive review. So the six
fields of §8.2 are a required constructor argument, not a logging convention, and a plan hash that
does not resolve is a refusal at construction time rather than a warning at read time. The refusal
type reused for that is `aleph.infer.refusal.ProvenanceIncomplete`, which already means exactly this
and already requires a remedy.

The two rules with teeth
------------------------
**Missing modality is an explicit mask, never a zero.** A `MaskedChannel` marked `ABSENT` may carry
nothing at all. A caller who substitutes a zero — the ordinary, well-intentioned way an absence turns
into a measurement of zero force, zero intensity, zero expression — gets `MissingModalityFilled`,
whose message names the substitution. This is index invariant I6 in the neural layer: the mask is the
statement, and the fusion rule (§10.1's calibrated product of experts, or a named equivalent) is what
consumes it. `UncertaintyBreakdown` extends the same rule one step: an output whose mask declares an
absent modality must attribute uncertainty to `MISSING_MODALITY`, or the mask has been recorded and
then ignored.

**`NativeQueryRequest` is a return value, not an exception.** A reader that cannot answer inside its
trust region *succeeds* by asking for a native run. So `NativeQueryRequest` is a `LaneOutput` like the
other five: same base class, same provenance block, same place in the union a caller switches on. It
is not an `Exception` subclass and a test asserts that it never becomes one. Its executable payload is
`aleph.infer.refusal.NativeQueryRequest`, imported here as `SimulationRequest` — that type already
states which observable, at which parameters, how many samples and what result would lift the
refusal, and re-deriving it here would have produced two request formats for one scheduler to read.
What this module adds around it is the lane provenance the inference-layer type has no reason to
carry, and the rule that a request certifies nothing: its authority ceiling must be `NO_CLAIM`.

The authority ceiling is an ordered lattice, not a label
--------------------------------------------------------
`AuthorityCeiling` is totally ordered from `NO_CLAIM` up to `NATIVE_CERTIFIED`, and `meet` returns the
weaker of two. That single operation is what makes §8.1's "bounded by the weaker of its two parents"
arithmetic rather than prose: composing an `MTG_PN_EXP` result — permanently `NON_NATIVE` — into any
joint result yields `NON_NATIVE`, and no configuration flag can raise it. `aleph.learn.lanes` builds
the hybrid ceiling with this operation instead of typing a value in.

The accepted-state digest, and the one place §8.2 is underspecified
-------------------------------------------------------------------
§8.2 requires every output to carry `accepted_state_digest`. An `MTG_PN_EXP` output has no accepted
native state to point at — that is the whole content of its permanent non-native ceiling — so the
requirement as written is unsatisfiable for the one lane that could otherwise proceed. Rather than let
that field be filled with a plausible-looking placeholder, this module gives absence a name,
`NO_ACCEPTED_STATE`, and makes it load-bearing: an output that declares it may not carry a ceiling
that can certify a native claim. The sentinel is not hex and cannot be mistaken for a digest.

**Public surface:** `MissingModalityFilledError`, `AuthorityEscalationError`, `UninterpretableOutputError`, `UncertaintySourceOmittedError`, `MissingModalityFilled`, `AuthorityEscalation`, `UninterpretableOutput`, `UncertaintySourceOmitted`, `LaneId`, `AuthorityCeiling`, `Modality`, `ModalityPresence`, `FusionRule`, `UncertaintySource`, `OutputKind`, `OutputProvenance`, `MaskedChannel`, `ModalityMask`, `LaneOutput`, `PosteriorSamples`, `PredictiveObservation`, `InterventionResponse`, `UncertaintyBreakdown`, `OODVerdict`, `NativeQueryRequest`, `is_lane_output`, `is_native_query_request`, `output_kinds`

## `learn/training.py` — 425 lines

A training run is an object written in the tensor language, not a script that consumes it.

The failure this module exists to prevent is not dramatic. A representation layer ships; someone
writes a loader; the loader needs an axis order, so it hardcodes one; the memory budget needs a rank,
so it picks one; the split defaults to random. Nothing in that script is a lie, and afterwards the
model's provenance is a shell history, its axis order is not the one the plan derived, its rank was
never qualified, and its generalization number measures interpolation. Every one of those is
unrecoverable after the run and free to prevent before it.

So a training run is declared, hashed, and validated as data, and three of its properties carry the
whole argument:

* ``representation_plan_hash`` is required. A model whose axis order is not the derived one cannot be
  compared against one whose axis order is.
* Loss weights are slots, never defaults. This module cannot even hold a weight — a weight is a
  measured, preregistered quantity that lives in a decision record, and what a slot holds is the
  requirement that such a record exist.
* The authority ceiling is derived from the admitted evidence rather than asserted by the run.

Nothing here trains. `train()` exists with its full signature and refuses; the signature is the
design and the refusal is the evidence that the design was not quietly executed. See
`docs/design/ALEPH-TRAINING-ON-THE-LANGUAGE.md` and `ALEPH-PD-003`.

**Public surface:** `TrainingPlanError`, `EvidenceInadmissible`, `WeightSlotUnset`, `MalformedDigest`, `PreconditionUnmet`, `TrainingLane`, `WeightSlot`, `TraceDigest`, `AdmittedTrace`, `EvidenceAdmission`, `TrainingPlan`, `PreconditionStatus`, `TrainingPrecondition`, `unmet_preconditions`, `train`

---

# `represent/`

## `represent/blocks.py` — 1375 lines

State blocks — the census enters as a manifest, and what cannot be classified is reported.

Section 10 partitions the state into six blocks; the census supplies instances. This module owns
the partition and owns none of the instances, so adding a fifteenth component extends a manifest in
``aleph/state/`` and touches nothing here.

Why the derivation reports instead of deciding
----------------------------------------------
Specification section 3 says block assignment "is derived, not hand-written", from a component's
``owned_state`` field names together with its ``scope``. Deriving a block from a field *name* is
keyword matching, and section 11.2 item 8 records the defect that follows: rule order carries real
weight, because ``surface_traction_pn`` matches a mesh rule on ``surface`` and a connector rule on
``traction``, and whichever rule is written first silently wins.

Implementing it makes the defect sharper than "order matters". The keyword vocabulary is not free —
it comes from section 10's own *native form* column — and **that column is not disjoint**. It lists
"forces" under both ``OWNER_PARTICLES`` and ``CONNECTOR``, and "binding/topology state" under both
``OWNER_PARTICLES`` and ``EVENT_STREAM``. No rule set faithful to the specification can therefore be
single-valued, and any rule set that *is* single-valued got there by an ordering decision nobody
recorded.

So the derivation here has no order at all. Every rule is evaluated against every name, the matches
are collected as a set, and the outcome is one of three things:

* exactly one block matched — :attr:`DerivationOutcome.DERIVED`;
* nothing matched — :attr:`DerivationOutcome.UNCLASSIFIED`;
* more than one block matched — :attr:`DerivationOutcome.AMBIGUOUS`.

The last two are *reported*, never guessed. Adding a rule can turn ``DERIVED`` into ``AMBIGUOUS``
and can never quietly change one block into another, which is the property that makes the rule table
safe to extend. The override manifest is where an unresolved name is adjudicated, and it is
load-bearing rather than an escape hatch: an override that merely agrees with a confident derivation
is *refused* (:class:`OverrideNotLoadBearing`), so the manifest holds only decisions the derivation
could not make.

Why nothing here names a census module
--------------------------------------
Section 11.2 item 9: the three census manifests spell their aggregate tuples differently —
``CENSUS_ENVIRONMENT_SURFACE``, ``FRAME_FLUID_NUCLEUS_COMPONENTS``, ``ACTOMYOSIN_LOAD_PATH`` — so
binding by module attribute name is fragile in the worst way, because a renamed tuple fails as a
*missing owner* rather than as a missing import. Everything here binds through the frozen schema's
queries (``component_names``, ``scope_of``, ``owners_of_scope``, ``composite_groups``,
``incident_connectors``) and a module handed in where a schema was expected is refused by name
(:class:`CensusBoundByModuleAttribute`). ``tests/represent/test_blocks.py`` asserts that none of the
three spellings appears in this file's source.

Why the H/B/X records must be authored
--------------------------------------
Section 11.2 item 10. A :class:`HomogenizedLawSlot` requires a declared inadequacy and a re-entry
condition. Those cannot be lifted off whichever contract field happens to be non-empty: the record
exists precisely to require that someone wrote down what the lumping forbids, and synthesising it
from the census text fabricates the content instead of demanding it. Each of the three record types
therefore has a ``from_contract`` classmethod that exists only to refuse
(:class:`RecordNotAuthored`), because the constructor a reader reaches for should say why it is the
wrong one rather than not be there.

What this module refuses
------------------------
* a state name it cannot classify, or can classify two ways (:class:`BlockDerivationUnclassified`,
  :class:`BlockDerivationAmbiguous`) — never a guess, and never a default block;
* an override that repeats the derivation (:class:`OverrideNotLoadBearing`), or names a state the
  owner does not own (:class:`OverrideUnknownState`);
* state on a scope that holds none (:class:`StatelessScopeOwnsState`);
* an ``OBSERVATION`` block derived from the census (:class:`ObservationBlockNotDerivable`) — an
  observation operator is declared in ``aleph/observe/`` and is not a component's owned state;
* an ``H``/``B``/``X`` record manufactured from a contract (:class:`RecordNotAuthored`);
* a ``B``-scope conditioning axis asked to be state (:class:`ContextAxisIsNotState`);
* a claim made about an excluded region (:class:`ExcludedRegionClaimed`);
* a census module handed in where a schema was expected (:class:`CensusBoundByModuleAttribute`).

Nothing here imports ``aleph.state``. The schema is read through the structural protocols below, the
same way ``cuts.py`` reads a connector contract, so the two lanes can land in either order and there
is no second copy of a contract to drift.

**Public surface:** `StateBlock`, `BlockRefusal`, `BlockError`, `BlockDerivationUnresolved`, `BlockDerivationUnclassified`, `BlockDerivationAmbiguous`, `MalformedOverride`, `MalformedRule`, `CompositeGroupsDisagree`, `DuplicateOverride`, `OverrideNotLoadBearing`, `OverrideUnknownState`, `StatelessScopeOwnsState`, `UnknownScope`, `ObservationBlockNotDerivable`, `RecordNotAuthored`, `ContextAxisIsNotState`, `ExcludedRegionClaimed`, `CensusBoundByModuleAttribute`, `ComponentLike`, `SchemaQueries`, `BlockRule`, `state_name_tokens`, `DerivationOutcome`, `RuleMatch`, `StateNameAssignment`, `assert_resolved`, `BlockOverride`, `BlockOverrideManifest`, `derive_state_name`, `ComponentBlocks`, `derive_component_blocks`, `ConnectorBlockBinding`, `connector_block_for_owner`, `connector_block_binding`, `observation_block_members`, `SchemaBlockReport`, `derive_schema_blocks`, `assert_fully_classified`, `HomogenizedLawSlot`, `ContextConditioningAxis`, `ExcludedRegion`, `assert_component_scopes_are_recorded`

## `represent/cuts.py` — 813 lines

Graph cuts and symbolic bond structure — read off the connector graph, never off space.

Why this module refuses to know where anything is
-------------------------------------------------
A factorization introduces a bond wherever it splits the state, and the width of that bond is set by
how strongly the two sides interact. The whole design of this module rests on one claim about what
"interact" means here, taken from manuscript section 4 and restated by `ALEPH-PORT-1501` section 4:
**co-location is never a mechanical connection.**

Two owners may share a boundary in space and exchange nothing, because no connector between them was
declared. Cutting between them separates nothing that was ever joined, so the bond is empty and the
cut is free. Two owners may sit at opposite poles of the cell and be one declared connector apart —
a stress-fibre end and the collagen it pulls on, coupled through an adhesion — and cutting between
them severs a load path, so the bond is not empty and the cut is not free. Spatial distance predicts
neither case.

The consequence for this file is structural rather than stylistic: `OwnerGraph` has no coordinate,
no radius, no neighbour list and no way to acquire one. Its only notion of adjacency is
"a declared connector names both of you", and its only notion of distance is connector count — the
same distance the force-path BFS of `ALEPH-PORT-1501` section 4 walks. A module that could see
positions could be tempted to use them, and the temptation is worth designing out rather than
resisting: a factorization whose bonds followed spatial adjacency would be optimizing a graph the
physics does not have.

What a bond budget is, and what it deliberately is not
------------------------------------------------------
Bond dimension is the rank of the interaction across a cut. A rank is a measured quantity: you learn
it by decomposing an accepted native state and looking at where the spectrum falls off. No accepted
native trace exists, so no rank exists, and every function here that would return one raises
`BondDimensionUnmeasured` instead. This is the exact place where a project quietly types `rank=16`
into a constructor default and never revisits it; the refusal is cheaper than the archaeology.

What survives without a number is an *ordering*. `BondBudget` carries the set of connectors that
generate the bond, and one budget is not wider than another when its crossing set is a subset of the
other's. That is a partial order, not a total one — two cuts crossing disjoint connector sets are
simply incomparable, and saying so is more useful than inventing a tie-break. `BondBudget` therefore
implements `<=` and `>=` and deliberately does not implement `<`, `>` or `int()`; asking a budget for
a magnitude is the mistake, so the magnitude accessors refuse by name.

The five legality predicates
----------------------------
Each of C1 through C5 in `docs/design/ALEPH-TN-LANGUAGE-SPEC.md` section 4.1 is a separate
`check_cN_*` function raising a separate exception type, because a rule whose refusal cannot be
distinguished from the neighbouring rule's refusal cannot be counterexample-tested. `check_cut_legal`
runs them in order and raises the first; the individual predicates stay public so that a test can
put a cut in front of exactly one of them.

Nothing here owns a connector. The contracts live in `aleph/state/`, this module reads the fields it
needs off them through a structural protocol, and `owner_graph_from_contracts` indexes them without
copying them — so a change to a contract cannot leave a stale duplicate behind in the cut layer.

**Public surface:** `ConnectorLike`, `CutRefusal`, `CutIllegal`, `CutSplitsComposite`, `CutSplitsAdjointPair`, `CutSplitsInternal`, `CutCrossesUnimplemented`, `CutSizeSplitsComposite`, `CutNotAPartition`, `CutNamesUnknownOwner`, `CutNamesUnknownConnector`, `CutCrossingDisagreesWithGraph`, `BondDimensionUnmeasured`, `BondBudget`, `OwnerGraph`, `owner_graph_from_contracts`, `connector_distance`, `GraphCut`, `crossing_connector_names`, `bond_budget_for`, `build_cut`, `check_cut_is_a_partition`, `check_crossing_is_derived_from_the_graph`, `check_c1_composite_group_not_split`, `check_c2_adjoint_layer_declared`, `check_c3_internal_not_separated_from_host`, `check_c4_crossing_is_implemented`, `cut_size`, `check_c5_size_counts_a_composite_once`, `LegalityPredicate`, `check_cut_legal`

## `represent/errors.py` — 1042 lines

The error vector — eight channels kept apart, and the tolerances nobody has earned yet.

Why the channels are separate, and why that is a type rather than a convention
------------------------------------------------------------------------------
`docs/design/ALEPH-TN-LANGUAGE-SPEC.md` section 7 prohibits summing the eight channels into one
scalar. The reason is operational rather than aesthetic: the sum cannot say *which* invariant broke,
and the fallback a broken invariant licenses is different for each one. A representation that lost
mass must be replaced by the native block and must never be silently rescaled; a representation whose
truncation was too aggressive may simply be given a wider bond budget; a representation whose
predicted observation drifted needs a re-plan for the query family and may be perfectly good for
every other query. One number cannot select between those three, so a caller holding one number
either picks the wrong repair or picks none.

So the prohibition is enforced by the type. :class:`ErrorVector` has no ``total``, no ``float``, no
``__add__`` and no maximum, and each of those refuses by name rather than being absent — an absent
method raises ``TypeError: unsupported operand``, which reads as a Python limitation rather than as a
design decision. Underneath that, the stronger fact: a :class:`ChannelReading` carries no magnitude
field at all. There is nothing to add even before anything refuses to add it, because a magnitude is
an L2 quantity measured against an accepted native trace and no accepted native trace exists.

Why a tolerance is a slot and never a default
---------------------------------------------
A tolerance is the line between "this representation is good enough" and "fall back". Typing one in
while the design is being written sets the project's acceptance criterion by accident, and choosing
one after seeing the runs that motivated it is a gate fitted to its own data — which this project has
refused once already (`ALEPH-DQ-104`). So every channel's :class:`ToleranceSlot` is empty, reading an
empty one raises :class:`ToleranceUnset`, and filling one requires naming the decision record that
set it. `None` is not accepted as "no tolerance": :data:`UNFILLED` is a distinct value, because a
``None`` in a numeric field is the shape that later gets defaulted to zero.

The asymmetry that matters: ``BYTE_DECODE`` has no tolerance at all
-------------------------------------------------------------------
Seven channels are budgets. ``BYTE_DECODE`` is not. A lossy serialization round-trip means the bytes
on disk are not the state that was written, and there is no amount of that which is acceptable —
every downstream number is then computed against a state nobody wrote. Section 7 spells this "hard
failure; never a tolerance", so :class:`ToleranceSlot` refuses to exist for that channel and
:class:`FallbackRecord` refuses to name it. The two refusals are separate types because they are
separate mistakes: one is budgeting the unbudgetable, the other is continuing past it.

Fallback is recorded or it did not happen
-----------------------------------------
A silent fallback is the defect this module exists to make impossible. A plan that quietly reverted
to the native block still produces artifacts, and those artifacts are then read as results of the
plan that was declared rather than of the plan that ran. So :func:`take_fallback` is the only way to
fall back, it requires a :class:`FallbackRecord`, and the record names the channel, the trigger and
the replacement form. :func:`assert_record_travels_with` is the second half of the same rule: a
record that does not travel with the artifact is a record nobody reading the artifact will see.

**Public surface:** `ErrorChannelError`, `ToleranceUnset`, `ToleranceDecisionUnrecorded`, `ToleranceAlreadySet`, `HardFailureHasNoTolerance`, `MalformedToleranceSlot`, `UnrecordedFallback`, `FallbackNotDeclared`, `HardFailureIsNotAFallback`, `VectorMissingChannel`, `ErrorChannelsNotSummable`, `ChannelErrorUnmeasured`, `ErrorChannel`, `ToleranceDiscipline`, `ToleranceForm`, `tolerance_form`, `FallbackForm`, `ChannelDeclaration`, `assert_channel_table_is_exhaustive`, `channel_declaration`, `budgeted_channels`, `ToleranceSlot`, `empty_tolerance_slots`, `ChannelReading`, `ErrorVector`, `empty_error_vector`, `FallbackRecord`, `take_fallback`, `assert_record_travels_with`, `assert_every_artifact_is_covered`

## `represent/families.py` — 861 lines

Factorization families — which family may be *planned* for which block, and in what order of cost.

Section 6 of the language specification is a table of thirteen families and the blocks each may
factor. This module is that table as data, plus the three things the table implies and does not say
out loud.

**Selecting a family is L1 and free; claiming a family works is L2 and sealed.** So every entry point
here answers "may this be planned?" and none answers "does this work?". The distinction is not
bookkeeping: a family becomes admissible for a block by an argument about structure, and it becomes
*correct* for a block by a measured error against an accepted native trace. No such trace exists, so
the second question has no answer and is refused rather than defaulted.

Three things the table implies
------------------------------
**Qualification.** ``DENSE`` is the always-legal fallback and never needs qualification, because a
dense block is the native block re-indexed — it discards nothing, so there is no discarded quantity
to measure. Every other family here approximates: it discards bond weight, projects onto a subspace,
substitutes an ensemble, or learns a map. Each therefore requires qualification before it may be
claimed to work, and :func:`requires_qualification` says so per family rather than leaving it to a
reader's judgement. ``PARTICLE_ENSEMBLE`` is included in that list deliberately: it is the *fallback
form* when a channel exceeds tolerance, which makes it the answer to a failed qualification and not
an exemption from one.

**Learned families are legal to plan and refuse to instantiate.** ``NEURAL_OPERATOR``,
``EQUIVARIANT_PARTICLE`` and ``MARKED_PROCESS`` are marked "learned; L2 only" in the table. A plan
that names one is a legitimate declaration of intent — it says which block the learned map would
serve and lets the plan hash record it — and :func:`instantiate_family` raises
:class:`~aleph.represent.RepresentationSealed` when anyone tries to build it. The refusing signature
is worth more than an absent one, for the reason ``aleph/represent/__init__.py`` gives: it records
what instantiation will need as input and is the evidence that the design was written without being
quietly executed against data that cannot support it.

**Cost is an ordering, never a number.** The same discipline ``cuts.BondBudget`` applies to bond
width applies here. A cost number for a family is a rank times a width times an axis count, and none
of those exist before measurement, so :func:`cost_of` refuses by name. What survives is a partial
order over *declared structure*: one family is not more costly than another when its declaration is
a special case of the other's — a linear subspace is a train with one bond, a train is a tree whose
tree is a chain, a tree is a graph with no cycle. Those are implications, not measurements, and each
is recorded with its ground in :data:`COST_RELATIONS`. Families whose cost is an architecture choice
rather than a structural fact — every learned one — are refused a place in the order
(:class:`CostOrderingUnavailable`) instead of being given a plausible one.

What this module refuses
------------------------
* a family planned for a block section 6 does not admit (:class:`FamilyNotAdmissible`);
* ``MPO`` asked to factor a state block (:class:`FamilyIsNotOverStates`) — the table's admissibility
  entry for it reads "operators, not states", so its admissible *state*-block set is empty and the
  refusal says why rather than reporting an empty set;
* a restricted admissibility used without acknowledging the restriction
  (:class:`RestrictionNotAccepted`) — section 6 admits ``TT_MPS`` on particles "local only", and a
  restriction nobody has to acknowledge is a footnote;
* a learned family instantiated (:class:`LearnedFamilySealed`), and any family instantiated from a
  declaration layer that owns no arrays (:class:`DeclarationLayerCannotInstantiate`);
* a family asked for a cost number (:class:`CostIsNotANumber`) or for its place in an order it has
  no structural claim to (:class:`CostOrderingUnavailable`).

**Public surface:** `FactorizationFamily`, `FamilyRefusal`, `FamilyError`, `FamilyNotAdmissible`, `FamilyIsNotOverStates`, `RestrictionNotAccepted`, `CostIsNotANumber`, `CostOrderingUnavailable`, `LearnedFamilySealed`, `DeclarationLayerCannotInstantiate`, `Admissibility`, `FamilyDeclaration`, `declaration_of`, `is_learned`, `requires_qualification`, `admissibility_for`, `is_admissible`, `assert_admissible`, `admissible_blocks`, `admissible_families`, `CostRelation`, `not_more_costly_than`, `assert_not_more_costly_than`, `are_cost_comparable`, `strictly_cheaper_than`, `minimal_cost_families`, `cost_of`, `FactorizationChoice`, `plan_factorization`, `instantiate_family`

## `represent/index.py` — 828 lines

The index algebra — the layer that refuses to let an axis be a number.

An index is a tensor leg, and the whole reason this module exists is that a tensor leg written as
an integer has already thrown away every fact a factorization plan needs. ``4096`` does not say who
varies along the axis, in what unit, in which basis, or whether the count changes when the topology
jumps; it says only how much memory to allocate. So nothing here can hold a dimension. An extent is
a symbol with an owner, and the two places where a small integer *is* a definition rather than a
measurement — a spatial component and a connector's endpoint arity — are the only places an integer
may appear at all, are restricted to the index kinds that define them, and are named as
definitional in the type itself.

What this module refuses, and why each refusal is a type rather than a comment:

* an extent asserted as a number (:class:`AssertedNumericExtent`) — a rank typed in before it was
  measured is the failure `ALEPH-PD-003` re-cut the seal to make unwritable;
* an index that names no owner, or a bond that names one (:class:`UnknownOwner`,
  :class:`BondHasOwner`) — a leg with no owner cannot be placed in the owner graph, and a bond with
  one pretends a factorization artefact is physical state;
* state on a scope that has none, and context smuggled in as state (:class:`ScopeHasNoState`,
  :class:`ContextAsState`) — this is where the ``E``/``I``/``H``/``B``/``X`` discipline stops being
  prose;
* a contraction between legs whose unit, basis or extent disagree
  (:class:`UnitMismatch`, :class:`BasisMismatch`, :class:`ExtentMismatch`);
* a sector-dependent axis inside a factorization declared over more than one sector
  (:class:`SectorExtentUnknown`) — an extent that changes under the plan is not an extent;
* an index for a declared-but-absent owner (:class:`AbsentNotZero`) — absence is a refusal, never a
  measurement of zero, and hoisting that rule from the observation layer into the algebra costs
  nothing here and is unenforceable later;
* two indices sharing a name (:class:`DuplicateIndex`) — a plan hash over an ambiguous name is a
  hash of whichever one the dictionary happened to keep.

Invariants I1, I2, I3 and I6 need to know what the frozen schema declares, so they are functions
that take a :class:`SchemaView` rather than validation inside ``__post_init__``. The view is a
projection — owner names, their scope letters, which owners this world does not instantiate, and
the declared unit symbols — and it is duck-typed onto the state contract rather than imported from
it, so this module holds no import of ``aleph.state`` at any scope and the two lanes can land in
either order.

**Public surface:** `IndexAlgebraError`, `MalformedIndex`, `AssertedNumericExtent`, `NotADefinitionalSize`, `ExtentKindMismatch`, `BondHasOwner`, `UnknownOwner`, `ScopeHasNoState`, `ContextAsState`, `ContractionRefused`, `UnitMismatch`, `BasisMismatch`, `ExtentMismatch`, `UnboundedNotContractible`, `SectorExtentUnknown`, `AbsentNotZero`, `DuplicateIndex`, `UnitNotDeclared`, `UnitRegistryUnavailable`, `IndexKind`, `BasisKind`, `ExtentKind`, `Extent`, `Index`, `SchemaView`, `check_owner_registered`, `check_scope_permits_state`, `check_context_axis_is_conditioning`, `check_owner_present`, `check_index_against_schema`, `assert_contractible`, `contractible`, `check_sector_axes`, `assert_unique_names`, `check_unit_declared`

## `represent/plan.py` — 559 lines

The L1 plan: which family factors which block, in what axis order, across which cuts.

An axis order that someone *chose* is an axis order nobody can reproduce. Specification section 5
therefore makes the order a total function of the owner graph, the connector cut sizes and the query
family, and this module is that function plus the record it produces. Two runs over the same frozen
schema and the same query yield the same order and the same hash, which is what lets a later result
say *which* plan it was computed under and be checked.

Three things here are deliberate and worth stating, because each is a place where the easy version is
wrong:

* **Distance is connector hops, never micrometres.** Two owners adjacent in space with no declared
  connector have no interaction to factorize across. `aleph.represent.cuts` owns that distance and
  this module never invents another.
* **The tie-break rule is recorded in the plan, not assumed.** At equal BFS depth, cut size and query
  relevance disagree — in the V1 vertical they disagree on the first realistic case — and the
  specification settled it without noticing there was a question. So the plan carries which rule it
  used, the rule is part of the hash, and a future change re-plans visibly instead of silently
  producing different numbers under the same name.
* **Nothing here is a bond dimension.** A rank is measured. `cuts.BondDimensionUnmeasured` is the
  answer to every question of the form "how wide", and this module does not route around it.

Declaration only, per `ALEPH-PD-002` and the unresolved
`docs/decisions/PROPOSAL-resolve-the-seal-contradiction.md`: no array library, no fitted quantity, and
no numeric literal standing in for one.

**Public surface:** `PlanRefusal`, `QueryFamilyEmpty`, `OwnerNotInPlan`, `OwnerUnreachable`, `AxisWithoutFamily`, `BlockUnknown`, `PlanReconciliationFailed`, `TieBreak`, `QueryFamily`, `PlannedAxis`, `BlockAssignment`, `derive_owner_order`, `derive_axis_order`, `FactorizationPlan`, `build_plan`, `assert_block_order_matches_blocks_module`, `assert_families_are_declared`

## `represent/qualification.py` — 668 lines

The L2 qualification protocol — written so it can be refused precisely and executed never.

What this module is for
-----------------------
Qualification is the act of measuring how far a representation is from an accepted native trace and
deciding whether that distance is acceptable. It is the one thing `ALEPH-TN-LANGUAGE-SPEC` section 1
seals: declaration is open, planning is open, *qualification is not*. The seal is not a schedule
slip. An error is a distance from a reference, no accepted native trajectory exists, and a project
that qualifies before it has a reference ends up choosing the reference to suit the representation —
tuning the physics until the factorization looks good, which is the failure that is invisible
afterwards because every number agrees with every other number.

So every entry point here has a full, honest signature and no body. The signature is the design: it
records what qualification will need — which plan, against which trace, on which targets, judged
against which tolerances — and getting that written down now is most of the value, because it is the
list a later reader checks the eventual implementation against. The refusal is the evidence that the
design was written without being quietly run against data that cannot support it. A module with the
computation commented out looks identical to one where the computation was never written; a module
that refuses by name does not.

What is declared as data, and why
---------------------------------
Manuscript section 10.2 names five things representation qualification runs on: analytic limits,
local projections, worst ROIs, events, and conservation/error budgets. They are five *different
questions*, and the list is worth more than its length: each target catches a failure the others
structurally cannot see. A global error budget averages away the error at a topology jump; a worst
ROI is the target that goes looking for it. An analytic limit is the only target that can fail when
the reference itself is wrong, because it is the only one whose ground truth is not the simulation.

:data:`QUALIFICATION_PROTOCOL` therefore records, per target, what it runs on, what it catches that
the others do not, and which section 7 channels it is evidence for. A qualification declared over
fewer than the five is refused by :func:`assert_protocol_is_complete` — not because five is a magic
number, but because dropping a target silently narrows what "qualified" means while leaving the word
unchanged.

What this module does not import
--------------------------------
The analytic-limit target names the closed-form oracle tree as a *string*. `aleph/**` may not import
`validation/**`: an oracle is the ground truth the runtime is judged against, so a runtime that can
import it can validate itself against it. The reference is by name, and resolving that name is the
job of whatever eventually runs the protocol, outside this package.

Nothing here imports a plan type either. `aleph/represent/plan.py` has not landed, and structural
protocols (:class:`PlanLike`, :class:`BlockLike`) name the fields this layer reads rather than the
classes that carry them — the same discipline `cuts.py` uses for connector contracts, and for the
same reason: the lanes then land in either order and there is no second copy of a contract to drift.

**Public surface:** `QualificationDeclarationError`, `QualificationProtocolIncomplete`, `MalformedTraceReference`, `TraceNotAccepted`, `BlockLike`, `PlanLike`, `NativeTraceRef`, `QualificationTarget`, `QualificationTargetDeclaration`, `target_declaration`, `channels_covered_by`, `assert_protocol_is_complete`, `QualificationVerdict`, `qualify`, `measure_channel`, `select_rank`, `compress`, `empty_verdict_shape`
