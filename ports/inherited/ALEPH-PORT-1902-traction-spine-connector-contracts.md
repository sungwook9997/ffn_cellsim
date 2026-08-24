# ALEPH-PORT-1902 — Stress-fiber, focal-adhesion and ECM traction-spine connector contracts (4)

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1902` |
| Lane | `L19 connector contracts — SF / focal-adhesion / ECM traction spine` |
| Status | `PROPOSED` |
| Port class | `CONCEPT_ONLY` — one composite-joint rule crossed as an argument. No code, no prose, no identifier, no constant. |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |

---

## 1. Aleph API

```python
from aleph.state.connectors_surface_traction import (
    ALPHA2BETA1_COLLAGEN_SERIES,
    ECM_CROSSLINK,
    ECM_FAR_FIELD_ANCHOR,
    FA_ACTIN_ANCHOR,
    INTEGRIN_COLLAGEN_CLUTCH,
    TRACTION_SPINE_CONNECTORS,
    composite_groups,
)
```

Aleph target file: `aleph/state/connectors_surface_traction.py`.
Controls: `tests/state/test_connectors_surface_traction.py`.

This entry covers the four Group-B contracts and the `composite_groups()` projection.
`ALEPH-PORT-1901` covers the eight Group-A contracts, the shared validation, and the read
provenance in full — the two entries name the same module because they are two groups of one
registry, and splitting the module would split a registry that has to be queried as one.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — verified `git rev-parse HEAD`, branch `codex/ff-ac-codex` |
| Source path | `ffn_sim/ac/engine/load_path.py`, `ecm_world.py`, `contracts.py`, `dispatch.py`, `sf_population.py`, `ecm_mechanics.py`, `composed_native.py`, `dump_state.py`, `ffn_sim/virtual_cell/grammar.py`, `ffn_sim/scripts/ac_gate_b_ecm_forces_native.py`, and the tests `tests/ac/engine/test_load_path.py`, `test_ecm_world.py`, `test_composition.py`, `test_ledger.py` |
| Source symbol(s) | `ConnectorContract` rows at `contracts.py:584-594, 744-755`; `resolve_fa_series_group` (`load_path.py:278-322`), `CompositeFAJointSpec` (`:325-343`), `resolve_fa_series_joint` (`:381-410`), `JointKind.FA_COMPOSITE_SERIES` (`:108-112`), `LoadPathJointRuntime` (`:651`) and its accumulation at `:724-747`, the segment-pair force kernel at `:523-574`, the host reference at `:459-515`; `CompositeLoadPathECMClutchAdapter` (`ecm_world.py:228-273`); `ECMInternalCrosslink` Protocol (`ecm_world.py:157-179`); `ECMBoundaryAnchorFacade` (`:349-392`) and `FarFieldDirichletBoundaryRuntime` (`:450`) with its reaction-ledger kernel at `:395-412`; `ECMWorld.accumulate_mechanics` dispatch sites at `:992-1002` |
| Read from | **working tree**, at a checkout whose `HEAD` is `be0e5876` |
| Working tree == commit? | **yes, for every cited file.** `git diff be0e5876 --stat -- <path>` returned empty for all of them. The tree carries **31** uncommitted changes overall; none of them touches any path cited in this entry. |

> Same caution as `ALEPH-PORT-1901` §2: this ledger records which revision was *read*, not which
> revision was named, because lane L9 found an earlier audit candidate that differed from the commit
> it was cited against.

## 3. Why source-derived porting beats clean-room

**For code, it does not, and nothing was taken.** The four contracts are transcribed from Aleph's
own Appendix A mechanics registry (Group B) and restated in Aleph's vocabulary. No identifier, no
kernel, no constant, no group-resolution routine and no prose crossed. The one string that appears
in both trees is `alpha2beta1_collagen_series`, and it appears in both because it is the name
Appendix A gives the mechanical group — it is Aleph's manuscript vocabulary, not the reference's
code vocabulary, and it is the *subject* of the contract rather than an inherited symbol.

**One rule crossed as an argument**, and it is the reason the contracts carry a
`composite_group` field at all:

> `fa_actin_anchor` and `integrin_collagen_clutch` form one mechanical group. They must be
> dispatched as one series joint and must not be evaluated as two independent springs.

Appendix A states it as a rule. What made it worth a ledger entry rather than a transcription is
that the reference had already worked out *how a registry lets you get it wrong* — two separately
declared semantic edges that a naive dispatcher will happily evaluate twice — and had built two
independent guards against it. Reading those guards told me what the contract-level expression of
the rule has to be: a shared group name, plus a projection that refuses a group of size one. That
is a case analysis discovered by somebody else's build, which `ports/TEMPLATE.md` §3 names as a
legitimate reason to let something cross.

The derivation of *why* the rule is physics rather than bookkeeping is in §4 and was worked out
here, not read there.

## 4. Physical or mathematical law represented

**The series-joint law.** Two elastic elements in series share one load and add their extensions.
For compliances `c₁ = 1/k₁` and `c₂ = 1/k₂`:

```
F₁ = F₂ = F           δ = δ₁ + δ₂ = (c₁ + c₂)·F           k_series = 1/(c₁ + c₂) = k₁k₂/(k₁+k₂)
```

Three consequences, each of which is a distinct way that evaluating the two halves independently
goes wrong:

1. **The load is shared, not assigned.** Evaluating each half against its own kinematics gives each
   the load *its own* extension implies, so `F₁ ≠ F₂` in general. A joint has no way to store the
   difference `F₁ − F₂`, so the adhesion plaque becomes a hidden force source. Nothing in a
   per-owner force check sees it, because each owner's books balance against the force it was
   handed.
2. **The joint is softer than either half**, `k_series < min(k₁, k₂)`, and it is *governed by the
   softer half*: if `k₁ ≪ k₂` then `k_series → k₁`. Two independent springs evaluated in parallel
   instead give an effective `k₁ + k₂`, which is stiffer than either. Since traction scales with
   joint stiffness at fixed strain, an independent evaluation **over-reports traction**, and it does
   so most severely exactly when one half is compliant — which is the physically interesting case.
3. **The clutch can open.** `integrin_collagen_clutch` is kinetic: its engagement is a discrete
   event, so the joint's load path can close and open. Evaluated separately, the actin anchor can
   report load while the ligand-side clutch is disengaged. That is a traction reading with nothing
   on the other end of it, and it is the failure mode that looks most like working code, because a
   single spring evaluates cleanly and returns a plausible number.

**Why the two are `K` and not `C`.** Both own an accepted kinetic state — anchoring/maturation on
the actin side, bound/unbound on the ligand side. Both therefore carry `commit_on_accept=True`: a
bond formed or broken inside a candidate step that is subsequently rejected leaves the accepted
state describing an adhesion that does not exist.

**`ecm_crosslink`: why an internal connector is still a connector.** Both endpoints are `ecm`, so
`internal_to="ecm"` and it is not a cross-owner load path. It is nonetheless a connector rather
than the matrix's private bookkeeping, because it owns **topology**: a remapped crosslink changes
*which* segments carry load. That change is irreversible and must commit only inside an accepted
step. It is also the one case in this lane where a force check is structurally incapable of catching
the error, because the force check is computed *from* the topology — a matrix connected differently
than the accepted state records will produce a perfectly self-consistent set of forces for the wrong
network.

**`ecm_far_field_anchor`: the only legitimate momentum sink.** A simulated matrix is finite; the
tissue it stands for is not. This connector supplies the reaction the omitted matrix would have
provided at the truncation boundary. Without it the whole assembly is a free body: `Σ F_ext = 0`, so
cell-generated contraction can only translate the assembly and traction is not measurable at all.
With it, `Σ F_ext = −R_far-field ≠ 0` legitimately, which is the one place in this registry where a
global force-closure check must account for an external term explicitly rather than asserting the
assembly is isolated. Getting that wrong in either direction is a live hazard: omit the anchor and
traction reads zero; forget it in the closure check and a correct assembly fails its own audit.

## 5. Units, domains, singular cases, invariants

This module holds no numeric quantity — it is a registry of declarations. The table records what
the declared connectors will carry, so a later lane cannot pick a different convention silently.

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| joint extension | µm | m | finite; sign meaningful |
| joint stiffness `k_series` | pN/µm | N/m | `> 0` |
| joint compliance `1/k` | µm/pN | m/N | `> 0`, additive in series |
| load through the joint | pN | N | finite |
| far-field reaction | pN | N | finite; the assembly's external term |
| clutch bound state | dimensionless | — | `{bound, unbound}`, committed on accept |
| crosslink topology | dimensionless | — | segment-pair identity, committed on accept |

Singular and boundary cases:

- **A composite group with one member** — refused by `composite_groups()`. A lone member is
  dispatched as an independent spring, which is precisely what the rule forbids, and it is the
  refusal most likely to be needed in practice: a member added, renamed or dropped leaves a group of
  one, and a group of one still evaluates.
- **A `K` contract with `commit_on_accept=False`** — refused. All four Group-B connectors except
  the far-field anchor are `K`.
- **`ecm_crosslink` without `internal_to="ecm"`** — refused, and refused in both directions:
  unset presents an internal connector as a cross-owner load path, and pointing at the wrong owner
  is refused too.
- **`internal_to` set on a cross-owner connector** — refused; it would hide a real cross-owner load
  path from cross-owner dispatch.
- **`world_boundary` as an endpoint** — legal for exactly one connector. More than one route out of
  the assembly, and the external term in the force closure is ambiguous.

Invariants, each with the test that asserts it (all in
`tests/state/test_connectors_surface_traction.py`):

- **I1.** The four Group-B contracts are constructible, in Appendix-A order.
  `::test_all_twelve_contracts_are_constructible_and_registered`
- **I2.** `fa_actin_anchor` and `integrin_collagen_clutch` carry the same
  `composite_group`, the group is exactly that pair, and no other contract in the registry sets the
  field. `::test_the_two_clutch_connectors_share_one_composite_group`
- **I3.** Both halves state the series rule in their own interpretation, so a reader who opens
  either one meets the prohibition rather than only the flag.
  `::test_the_series_joint_rule_is_stated_in_both_halves`
- **I4.** A composite group with one member is refused, and the pair is accepted.
  `::test_a_composite_group_with_one_member_is_refused`
- **I5.** `ecm_crosslink` has both endpoints `ecm` and sets `internal_to="ecm"`; it is the only
  internal contract in the registry. `::test_ecm_crosslink_is_internal_to_the_matrix`
- **I6.** No cross-owner contract claims to be internal.
  `::test_no_cross_owner_contract_claims_to_be_internal`
- **I7.** All three Group-B `K` contracts commit only on an accepted step.
  `::test_every_kinetic_connector_commits_only_on_an_accepted_step`
- **I8.** `ecm_far_field_anchor` is the only connector touching `world_boundary`, and its
  interpretation records that it makes the assembly's net force legitimately non-zero.
  `::test_the_far_field_anchor_is_the_only_route_out_of_the_assembly`
- **I9.** Adhesion to the matrix (`integrin_collagen_clutch`, `K`) and contact with the matrix
  (`membrane_ecm_contact`, `C`) are different connectors in different families, and the contact one
  is not in any composite group.
  `::test_the_membrane_ecm_contact_is_contact_and_not_adhesion`

## 6. Source evidence class and known retractions

Evidence class of the audit: `AUDIT_READ` at `be0e5876` — the source and its construction graph
were read; nothing in the reference was executed. No numeric value was inherited.

| Connector | Verdict at `be0e5876` | Decisive evidence |
|---|---|---|
| `fa_actin_anchor` | **SEAMED** | A real adjoint segment-pair kernel exists (`load_path.py:523-574`, scattering `+f`/`−f` by barycentric weight at `:568-571`), but its runtime is constructed **only** in `tests/ac/engine/test_load_path.py:246`. |
| `integrin_collagen_clutch` | **SEAMED** | Same object — one composite joint, so one verdict. Its ECM-side adapter (`ecm_world.py:228-273`) is constructed only at `tests/ac/engine/test_ecm_world.py:368`. |
| `ecm_crosslink` | **SEAMED** | `ecm_world.py:157-179` is a `Protocol` whose `accumulate_internal` (`:166-167`) is docstring-only. The only implementation in the tree is a test spy at `tests/ac/engine/test_ecm_world.py:140`. |
| `ecm_far_field_anchor` | **SEAMED** | A genuinely real kernel-backed Dirichlet closure exists — reaction ledger kernel `ecm_world.py:395-412`, rollback `:415-424`, commit `:427-447`, accumulation `:518-525` — and is instantiated only in tests (`test_ecm_world.py:59, 799`, `test_ledger.py:187`). |

None of the four is `NAME_ONLY`: each has a contract row, a facade claim (`dispatch.py:197-201`)
and, for three of them, real code. That is what makes these seams worth recording — two of them are
seams **with working implementations behind them**, which is a different and more misleading state
than an empty declaration.

### 6.1 The composite series joint: confirmed, and better than the brief said

I was told `load_path.py` implements the composite series joint. **Confirmed, and it is stricter
than "implements" suggests** — there are two independent guards, and the reference's own words:

- The group must be exactly the two edges. `load_path.py:285-288` raises if the resolved group does
  not contain exactly the two semantic edge names.
- The resolution must produce exactly **one** mechanical joint. `load_path.py:398-410` raises
  `"FA semantic edges would double-count mechanical stiffness"` when the mechanical joint count is
  not 1, and the resulting joint carries **both** semantic edges under a single
  `FA_COMPOSITE_SERIES` kind.
- The ECM facade re-asserts the same invariant at build time, independently:
  `ecm_world.py:891-892` requires the resolution to name the group *and* to have a mechanical joint
  count of 1.
- Dispatch happens exactly once. `ecm_world.py:255-257` is a method whose entire body is one call,
  under the docstring *"Dispatch the mechanical group exactly once."*
- The force kernel launches once per joint and scatters one force pair from one spring energy
  (`load_path.py:724-747`, kernel `:523-574`).

So the answer to "one series joint or two independent springs" is: **genuinely one joint**, and the
reference guards it in two places. This is the one place in the whole L19 audit where the source is
stronger than its reputation, and it is recorded here as such.

**Two caveats, both of which change what a reader should conclude:**

- **Caveat 1, and it is a real gap: the joint is composite, not a computed series reduction.** The
  joint's stiffness (`load_path.py:333`) is a single caller-supplied scalar. There is **no**
  `1/k = 1/k₁ + 1/k₂` combination anywhere in that file. The two chemistries are validated against
  each other but never composed numerically. So the "series joint" is a structural guarantee — one
  spring, one dispatch, one force pair — and **not** an implementation of §4's compliance
  addition. The prohibition on double-counting is enforced; the series *law* is delegated to
  whoever supplies the scalar. Aleph must not inherit that: `ALPHA2BETA1_COLLAGEN_SERIES` is a
  contract-level grouping obligation, and any Aleph evaluator that claims it owes a control against
  `k₁k₂/(k₁+k₂)` and against the `k₁ ≪ k₂` limit.
- **Caveat 2: it never runs in production.** The joint runtime is instantiated only in tests. The
  one real composed-world path routes the ECM facade to the ECM *component's* accumulation
  (`composed_native.py:442-443, :453-463`), which is component mechanics and not connector
  dispatch, and the CUDA ECM gate script does the same (`scripts/ac_gate_b_ecm_forces_native.py:169-175`).

### 6.2 The finding with the longest reach

The bindings object that carries all five ECM connector delegates is constructed **only** in
`tests/ac/engine/test_ecm_world.py` (`:314, 537, 545`). Since the ECM facade's accumulation method
takes that bindings object as its only argument, **none of the five ECM connector edges has ever
been dispatched outside the test suite** — including the far-field anchor, whose kernels are real
and complete.

And the coverage gate does not notice, for a reason worth stating precisely: the ECM facade *is* one
of the three facades the composed world treats as real, so its five connector claims are
structurally covered while being numerically absent. The facades outside that set of three are
mapped onto a stand-in that the reference itself documents as *"It launches nothing"*
(`composed_native.py:296-303`, mapping at `:453-463`).

**The lesson for Aleph, which is about Aleph's own gate and not about the reference:** a connector
can be genuinely implemented and genuinely unreachable at the same time, and which it is depends on
*which entry point you call production*. `ALEPH-PORT-303` already requires scheduled + reached +
evaluated-something. This says the definition of "reached" must name the entry point, or a facade
that is real for one driver will certify connectors that no driver reaches.

## 7. Independent oracle or derivation

- **Appendix A Group B and its composite-joint rule**, `docs/manuscripts/extracted/APPENDIX_A_mechanics_registry.txt` — Aleph's own manuscript, checked line by line against the four contracts.
- **The series-compliance law of §4**, derived here. It is what tells you the joint is softer than
  either half and governed by the softer one, which is the quantitative statement any future
  evaluator of `ALPHA2BETA1_COLLAGEN_SERIES` must be controlled against — and which the reference's
  caller-supplied scalar does not provide.
- **`ALEPH-PD-002`**, which fixes the contract field set including `composite_group` and
  `internal_to`.
- **Newton's third law at the truncation boundary**, which is what makes `ecm_far_field_anchor` a
  connector rather than a boundary condition: the reaction is applied to a named second endpoint and
  appears in the assembly's external term.
- **Mutation**, run rather than asserted — see `ALEPH-PORT-1901` §9. The relevant one for this
  entry is the kinetic-commit mutation, which took down both the registry-wide `K` assertion and its
  refusal control.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_two_clutch_connectors_share_one_composite_group` | Both halves carry `alpha2beta1_collagen_series`, the group is exactly that ordered pair, and no other contract sets the field. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_series_joint_rule_is_stated_in_both_halves` | Each half's own interpretation names the series joint, names the group, and states it is not an independent spring. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_ecm_crosslink_is_internal_to_the_matrix` | Both endpoints `ecm`, `internal_to="ecm"`, and it is the only internal contract in the registry. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_every_kinetic_connector_commits_only_on_an_accepted_step` | Exactly four `K` contracts, all with `commit_on_accept=True`; three of the four are in this group. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_far_field_anchor_is_the_only_route_out_of_the_assembly` | The anchor terminates at `world_boundary`, records the non-zero net force, and is the only connector touching that owner. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_the_membrane_ecm_contact_is_contact_and_not_adhesion` | Matrix contact and matrix adhesion are distinct contracts in distinct families, and only the adhesion one is in a composite group. |
| Positive | `tests/state/test_connectors_surface_traction.py::test_all_twelve_contracts_are_constructible_and_registered` | The four Group-B names, in Appendix-A order. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_a_composite_group_with_one_member_is_refused` | A group of one is refused, and the pair is accepted — so the refusal is about the missing half and not about grouping. This is the control for the composite rule: a lone member would be dispatched as an independent spring and would evaluate cleanly. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_a_kinetic_connector_that_does_not_commit_on_accept_is_refused` | A `K` contract that does not commit on accept is refused, and the corrected version is accepted. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_a_same_owner_connector_without_internal_to_is_refused` | An `ecm`↔`ecm` contract with `internal_to` unset is refused, and setting it to `"ecm"` is accepted. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_internal_to_naming_the_wrong_owner_is_refused` | `internal_to="cortex"` on an `ecm`↔`ecm` contract is refused. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_internal_to_on_a_cross_owner_connector_is_refused` | The inverse: a genuine cross-owner load path may not be marked internal and hidden from cross-owner dispatch. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_an_empty_composite_group_is_refused` | A whitespace group name would create a group nothing can join, silently making both halves lone members. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_the_validator_is_not_vacuous` | The base keyword arguments every refusal above mutates are themselves accepted, so a validator that raised unconditionally fails here. |
| Negative (must fail) | `tests/state/test_connectors_surface_traction.py::test_every_contract_is_honestly_unimplemented` | All four are `implemented=False`. This test is designed to break the day one of them acquires an evaluator, forcing the flag, the evaluator and its control to land together. |

## 10. Numerical and precision envelope

**No floating-point arithmetic occurs.** Every field is a string, a boolean, a `StrEnum` member or
`None`, so there is no working precision, no accumulation precision and no conditioning to state.
The controls are exact: booleans compared with `is True` / `is False`, families compared against
enum members, group membership compared as tuple equality against a written-out pair.

The precision statement that belongs here is the one a future evaluator inherits, recorded now
because §6.1 shows it is exactly where the reference stopped: the series reduction
`k = k₁k₂/(k₁+k₂)` is ill-conditioned as written when `k₁ + k₂` is small, and its useful form is
`k = 1/(1/k₁ + 1/k₂)` evaluated in compliance space, which is well conditioned across the whole
stiffness-ratio range and degrades gracefully as either half goes slack. Any evaluator claiming
`ALPHA2BETA1_COLLAGEN_SERIES` owes a control at the `k₁ ≪ k₂` limit — the compliant-half case,
which is the physically interesting one and the one an independent-spring evaluation gets most
wrong.

## 11. Production-backend residency and transfer

Host-side Python only, permanently. Immutable frozen dataclasses constructed at import time; no
device transfer, no per-step cost, nothing that runs inside a step. Whatever assembles the dispatch
schedule reads this registry once, on the host.

Two residency obligations this entry fixes for future evaluators, because they are properties of the
contract rather than of an implementation:

- The composite group must be resolved to one joint **before** the schedule is built, on the host.
  Resolving it per step would put a topology decision inside the step, where a rejected candidate
  could leave it half-applied.
- `ecm_crosslink` owns topology, so a remap changes the connectivity arrays a device kernel reads.
  Its commit must therefore be part of the accepted-step transaction and cannot be a device-side
  side effect — which is what `commit_on_accept=True` declares.

`aleph/state/**` may not import `validation/**`; this module imports only `dataclasses`, `enum` and
`typing`.

## 12. Comments and docstrings to discard

No source prose was carried. Deliberately absent from `aleph/state/connectors_surface_traction.py`:
the repository name, the package namespace, module paths, the joint-resolution and runtime class
names, the joint-kind enum member, the adapter and facade names, the bindings type name, the
Protocol names, the kernel names, the gate and milestone labels, and the reference's error-message
wording — including the double-counting message quoted in §6.1, which is quoted here because the
quote *is* the evidence and must not appear in the package.

Also not carried: the reference's connector-family vocabulary, which differs from Appendix A's
two-family legend; and its caller-supplied joint stiffness, which §6.1 explains is a gap rather
than an asset and which Aleph must not inherit as a convention.

What replaces them: §4's compliance derivation, and the two contracts' own
`mechanical_interpretation` fields, which state the series rule and the reason for it in terms with
no proper nouns — asserted present in both halves by
`::test_the_series_joint_rule_is_stated_in_both_halves`.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`. The code named in §1 has landed and its controls pass — 43 tests in `tests/state/test_connectors_surface_traction.py`, plus three mutations recorded in `ALEPH-PORT-1901` §9 — but all four contracts are `implemented=False`, so there is no physics to accept. A contract is accepted when the load path it declares has an evaluator and a control, not when the declaration parses. |
| Reviewer | agent-proposed, **unratified**. Nobody has reviewed §4's series derivation, the §6.1 finding that the reference's series joint supplies no series reduction, or the §6.2 conclusion that none of the five ECM connector edges has run outside a test. |
| Rollback | Delete `aleph/state/connectors_surface_traction.py` and `tests/state/test_connectors_surface_traction.py` — the same two files as `ALEPH-PORT-1901`, since the two entries cover one module. Nothing imports them yet. What is lost is the record that the traction spine is four declarations with no evaluator. |

## 14. Honest limits

- **Unratified**, and §4's series derivation is the substantive claim nobody has checked.
- **No physics is established.** Nothing here evaluates a force. In particular there is no series
  reduction implemented in Aleph either — this entry declares that the two halves form one joint and
  derives what that joint must satisfy; it does not compute it. `implemented=False` on all four is
  the whole content of that.
- **`composite_group` is a grouping obligation, not a joint.** A shared string and a refusal of
  size-one groups is as far as a contract can go. Nothing in this module can prevent a future
  dispatcher from reading the group and then evaluating the members separately anyway. That guard
  has to live in the dispatcher, and it does not exist yet — this is the single largest gap between
  what this entry declares and what would actually enforce Appendix A's rule.
- **Two of the four owners do not exist in Aleph.** `focal_adhesion`, `sf_arc`, `ecm` and
  `world_boundary` have no state, no geometry and no mechanics anywhere in `aleph/`. All four
  contracts are declarations about owners that are not represented, which is what `ALEPH-PD-002`
  asked for and is also a reason not to read this as describing anything runnable.
- **`ecm_crosslink`'s remodelling semantics are out of scope.** The contract says the connector owns
  binding, unbinding, remapping and remodelling-compatible topology, and commits on accept. It says
  nothing about crosslink lifetime, force-dependent rupture, rebinding rates, or how a remap is
  chosen. Nothing here is evidence about matrix remodelling.
- **`ecm_far_field_anchor` has no truncation policy.** The contract declares that a far-field
  reaction is transferred; it does not say where the boundary is, how far-field it has to be, or how
  the finite-domain error scales. That error is a real systematic on any traction measurement and it
  is untouched by this entry.
- **The audit is a read, not a run.** Nothing in the reference was executed. A body that exists and
  is reached could still be numerically wrong, and the `SEAMED` verdicts rest on a construction-site
  search that a sufficiently indirect factory could evade. §6.2's conclusion is the strongest claim
  here and it is the one most exposed to that limit.
- **The verdicts are pinned to `be0e5876`** and the reference tree carries 31 uncommitted changes.
  Every cited file was verified identical to that commit, but the tree will move and no claim is
  made about its present state.
- **The brief this lane was given was wrong twice, and the corrections are recorded rather than
  quietly applied.** See `ALEPH-PORT-1901` §6.2 (`membrane_medium_traction` — the exterior Stokes
  solve does exist; the reference's own inline admission is stale) and §6.3 (the "only three
  production accumulators" figure is off by roughly two orders of magnitude at the loosest reading
  and is 4 rather than 3 at the strictest). If either correction is wrong, it is mine and it is
  where a reviewer should start.
