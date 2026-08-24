# ALEPH-PORT-1702 — sf_arc: one owner, four subpopulations, and three declared-only connectors

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1702` |
| Lane | `L17 registry — actomyosin load-path compartments` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` — **verdict: RE-DERIVE.** One *argument* crosses; no code, no identifier, no constant. |

---

## 1. Aleph API

```python
from aleph.state.census_actomyosin import (
    SF_ARC,
    VENTRAL_STRESS_FIBER,
    DORSAL_STRESS_FIBER,
    TRANSVERSE_ARC,
    PERINUCLEAR_ACTIN_CAP,
    SF_ARC_SUBPOPULATIONS,
)
```

Aleph target file: `aleph/state/census_actomyosin.py`. Five entries: one `E` owner and its four `I`
subpopulations. No bundle mechanics, no crosslink law and no material card is authorised here.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/sf_mechanics.py`, `ffn_sim/ac/weave/regions.py`, `ffn_sim/ac/weave/woven_cell.py` |
| Source symbol(s) | `SF_CONNECTOR_BINDING_STATUS`, `SFMechanicsTopology`, `build_sf_mechanics_topology`, `SFFilamentMechanics`, `SFInternalArcJointConnector`, `link_spring_force_reference`, `bending_force_reference`, `WovenCell`, `region_slice` |
| Read from | **`git show be0e5876:<path>`** — the immutable commit. Two excerpts were printed and read directly (`sf_mechanics.py:103-124`, `woven_cell.py:213-222`). |
| Working tree == commit? | `yes` — `git diff be0e5876 -- ffn_sim/ac/engine/sf_mechanics.py` and `-- ffn_sim/ac/weave/woven_cell.py` and `-- ffn_sim/ac/weave/regions.py` are all empty. |

| Path | `sha256:` | Lines |
|---|---|---|
| `ffn_sim/ac/engine/sf_mechanics.py` | `sha256:e1bb8e05b47c3881` | 546 |
| `ffn_sim/ac/weave/regions.py` | `sha256:091ac7ea8d1a0645` | 598 |
| `ffn_sim/ac/weave/woven_cell.py` | `sha256:25a61890595df7b0` | 477 |

## 3. Why source-derived porting beats clean-room

**For code, it does not — RE-DERIVE.** `sf_mechanics.py` is the best-written file I read in that
repository: it imports its α-actinin stiffness from a single source-of-truth module instead of
re-typing the literal, it **raises** rather than defaults when the axial stiffness is not supplied, and
it keeps host NumPy reference implementations of its own force laws beside the device path. All of that
is craft Aleph should match. None of it is content Aleph cannot derive: a link spring and a bending
term are the two most standard laws in filament mechanics.

**For one argument, it does.** `SF_CONNECTOR_BINDING_STATUS` at `sf_mechanics.py:105-122` is a
per-connector honesty table, and reading it settled a question I had been given as an unverified claim.
Three of the five `sf_arc` connectors are marked `SEAMED`, and the *reasons* are the transferable part:

| Connector | Recorded status | Why it is stuck |
|---|---|---|
| `dorsal_arc_crosslink` | `KERNEL_BOUND` | Internal to `sf_arc`, so both ends live in one array. |
| `nmii_sf_motor` | `KERNEL_BOUND`, magnitudes gap | Two-array adjoint scatter exists; the force scale does not. |
| `sf_cortex_transient` | **`SEAMED`** | Needs a live cortex endpoint and a **two-array adjoint scatter**. |
| `if_sf_plectin` | **`SEAMED`** | Needs an intermediate-filament owner. |
| `mt_sf_spectraplakin` | **`SEAMED`** | Needs a microtubule owner. |

The pattern is the finding, and it is a *consequence of disjointness*, not an accident of scheduling.
The one connector that got bound is the one whose endpoints share an array. Every connector that must
cross between two separately-owned populations is stuck, because crossing requires scattering a force
into an array the evaluator does not own, with the equal and opposite reaction into its own —
Newton's third law across an ownership boundary. **Disjoint ownership makes cross-component coupling
strictly harder to implement, and that difficulty is the honest cost of the correct ontology.** A
single shared array makes every one of those connectors trivial, which is exactly why a project under
schedule pressure drifts toward one. Recorded here so Aleph pays the cost deliberately.

This is a case analysis discovered by somebody else's implementation attempt, which
`ports/TEMPLATE.md` §3 names as a legitimate thing to let cross. What crosses is the argument. No
code, no identifier, no constant, no status vocabulary.

## 4. Physical or mathematical law represented

No mechanics was ported. Two declarative laws are derived here.

**A subpopulation is a fibre of a labelling, not a component.** Let `S` be `sf_arc`'s filament set and
`ℓ : S → {ventral, dorsal, arc, cap}` a label. The four subpopulations are the fibres `ℓ⁻¹(·)`. They
partition `S`, they are distinguishable in every query, and they are **not** elements of the owner set
`O` from `ALEPH-PORT-1701` §4. Promoting a fibre to an owner does not add information — the fibre was
already resolvable — while it does add four owners that need connectors between them, none of which
the registry declares. So the promotion strictly loses: it costs a connector set and buys a
distinction that already existed.

**An adjoint transfer across an ownership boundary is what a connector is for.** For a connector
between owners `a` and `b`, the transfer must produce `f_a` and `f_b` from two independent
accumulations, with `f_a + f_b = 0` as a *check* rather than as an assignment. If `f_b := −f_a` is
written, closure becomes a tautology and can never fail. This is why the three `SEAMED` connectors are
hard: two arrays means two accumulations, and two accumulations means closure is falsifiable.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| bundle node position | µm | m | finite |
| bundle axial stiffness | pN/µm | N/m | > 0, **unsourced, no default** |
| bundle bending rigidity | pN·µm² | J·m | > 0, **unsourced** |
| crosslinker stiffness | pN/µm | N/m | > 0, **unsourced** |
| subpopulation label | 1 | 1 | one of four |

Singular cases:

- **A subpopulation asked for its own state block.** Refused. It owns nothing; its label is a field of
  the parent's arrays.
- **A subpopulation used as a connector endpoint.** Refused. A connector names `sf_arc` and selects an
  endpoint *role*, which is what the registry's endpoint-role column has always meant.
- **An axial stiffness requested with no source.** Refused, never defaulted. The reference project got
  this right and it is worth matching rather than inheriting.
- **A curved bundle treated as a straight one.** For the transverse arc the bending term is not a small
  correction to the axial term, which is why the material card is per-subpopulation rather than shared.

Invariants:

- **I1.** Each of the four subpopulations is `[I]` and names `sf_arc` as parent.
  `tests/state/test_census_actomyosin.py::test_each_sf_arc_subpopulation_names_sf_arc_as_parent`.
- **I2.** Each owns no independent state.
  `test_each_sf_arc_subpopulation_owns_no_independent_state`.
- **I3.** The parent is an `E` owner in this group and owns the state they are resolved in.
  `test_the_subpopulation_parent_is_an_explicit_owner_in_this_group`.
- **I4.** The distinguishing labels are state the parent owns, so the four are resolvable in a
  snapshot. `test_the_subpopulation_labels_are_state_the_parent_owns`.
- **I5.** The internal set is exactly the four the registry names — no fifth, no omission.
  `test_the_sf_arc_subpopulations_are_exactly_the_four_the_registry_names`.
- **I6.** `sf_arc` records the disjointness rule and does not share a state key.
  `test_each_filament_owner_records_the_disjointness_rule`,
  `test_no_state_key_is_claimed_by_two_owners`.
- **I7.** The perinuclear cap does not become the cortex's on grounds of proximity.
  `test_the_perinuclear_cap_does_not_share_filaments_with_the_cortex`.

## 6. Source evidence class and known retractions

I looked in that repository's audit directory, its roadmap, and its own test assertions.

- `AUDIT_AC_ENGINE_2026-07-25.md:30` rates `sf_arc` **`KERNEL_BOUND`, force-real**, on CUDA at
  **20–40 fibres**, with its own builder self-labelling "NOT the native SF inventory". `:167` records
  that its own "validate at full native" rule is violated in effect for `sf_arc` at 40 fibres.
- `:53` records `dorsal_arc_crosslink` as **`KERNEL_BOUND` but unreachable as a connector** — it
  exposes only an accumulate entry point with no snapshot, rollback, commit or ledger, while its
  contract declares kinetic state that commits on acceptance. Its own audit calls this the cleanest
  contract/runtime divergence in the tree. That is worth naming precisely because it is the *inverse*
  of the usual failure: not a declaration with no implementation, but an implementation that cannot
  satisfy its declaration.
- `PI_GAP_EVIDENCE_CARDS_2026-07-25.md:39` records the SF NMII magnitude as `SEAMED` and the traction
  magnitude as **`INVALID`**.
- Its own test module asserts the three `SEAMED` strings start with `"SEAMED"` — an honesty assertion,
  and a good one.

Reachability: `SFFilamentMechanics` is live at small fibre counts. The three `SEAMED` connectors have
**no runtime class anywhere** — I checked, and each name appears only in a contract declaration, a
name list, a status string, that assertion, and a synthetic visualisation script.

Tests: 11 tests, 266 lines, 35 doubles, no device allocation. Four are real physics (axial tension
under stretch, restoring bending force, the two sourced constants); the rest are plumbing. Mixed, and
the mix is the honest answer.

No retraction of the α-actinin or actin-persistence-length constants was found; I looked in the audit
directory and the gap-evidence cards. Neither constant was inherited anyway.

## 7. Independent oracle or derivation

For this entry's declarative content, the oracle is the same exact set-theoretic identity as
`ALEPH-PORT-1701` §7, applied to the labelling `ℓ`: the four fibres partition `S` iff each `I` entry
owns no state and names the same parent. Exact, so asserted with `==`.

For the mechanics a later lane would write, the oracle available to Aleph and *not* to that repository
is a pair of closed forms Aleph already owns: a taut-string equipartition spectrum
(`ALEPH-PORT-402`) for an axially tensioned bundle, and the discrete-bending scale-invariance identity
established in `ALEPH-PORT-1101` — a bending energy that is exactly scale-invariant contributes exactly
zero to a dilational virial, which gives an exact zero to test against rather than a fitted residual.
Neither needs the reference project open.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_actomyosin.py::test_each_sf_arc_subpopulation_names_sf_arc_as_parent` | All four `[I]` entries name `sf_arc`, and all four are `[I]`. |
| Positive | `tests/state/test_census_actomyosin.py::test_each_sf_arc_subpopulation_owns_no_independent_state` | All four own exactly `()`. |
| Positive | `tests/state/test_census_actomyosin.py::test_the_subpopulation_labels_are_state_the_parent_owns` | The label and material-card blocks are in `sf_arc.owned_state`, so the four are resolvable without being owners. |
| Positive | `tests/state/test_census_actomyosin.py::test_the_group_declares_five_explicit_owners_four_internal_and_one_homogenized` | The group's scope shape, so a silent retag is caught. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_promoting_a_subpopulation_to_a_component_is_refused` | The most tempting wrong move — making `transverse_arc` an `E` owner with its own filament block — is refused, and refused **on the parent rule**, because the contradiction is the tag rather than the arrays. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_an_internal_subpopulation_that_owns_state_is_refused` | An `[I]` entry that acquires a state block is refused; without this, four subpopulations could become four owners of one filament set. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_an_internal_subpopulation_with_no_parent_is_refused` | `[I]` with no parent is refused: internal means internal to something. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_a_subpopulation_whose_parent_is_not_in_this_group_is_refused` | Reparenting the dorsal fibres onto the cortex is refused — that would put `sf_arc`'s state inside another group's owner without that group agreeing. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_an_explicit_owner_that_names_a_parent_is_refused` | Giving `sf_arc` a parent is refused: an `E` owner with a parent is not independent. |

## 10. Numerical and precision envelope

No arithmetic crosses this boundary. Every control is an exact comparison — set equality on the
subpopulation names, identity on the parent string, `()` on the owned state — so no tolerance is
stated and none would be meaningful. The property being checked is discrete; a tolerance on it would
be a defect, not a loosening.

Precision statements are recorded here only for what a later lane inherits: the reference's own device
path runs its link-spring and bending kernels in float64, and its host references match. Aleph's
existing filament-adjacent work is float64 throughout with a mixed absolute/relative residual band.
Nothing in this entry depends on either.

## 11. Production-backend residency and transfer

Host-side declarations. No device array, no transfer, not in the step loop. When the `SEAMED`
connectors' Aleph counterparts are written, the residency question this entry raises for them is
concrete: a cross-ownership adjoint scatter needs **two** force accumulators resident on the backend
simultaneously and two independent scatter passes, which is a real cost and is the reason the
reference stalled. Recording it here means the cost is budgeted rather than discovered.

## 12. Comments and docstrings to discard

Read, and **not** carried: the connector-status vocabulary and every one of its status strings; the
evidence-rung names and their ordering; the topology, mechanics and connector class names; the module
paths; the gap-card identifiers; the roadmap commit references; the two sourced constants and their
citations.

The status *vocabulary* deserves a specific note, because it is the most tempting thing in the file.
It is a genuinely good idea — a per-connector honesty table — and Aleph already has the same idea in a
better form: `implemented: bool` on the contract plus `declared_but_unimplemented()` as a queryable
value. Adopting a second ladder alongside it would let a reader think an Aleph entry had been graded
on that scale. So the idea is kept and the vocabulary is dropped, and §3 above is the replacement
prose.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** `PROPOSED`. The five contract entries and their controls pass (`tests/state/test_census_actomyosin.py`, 88 passed), but the substantive claim of this entry is the argument in §3 about why cross-ownership connectors are hard, and no human has reviewed it. |
| Reviewer | agent-proposed, **unratified**. |
| Rollback | Delete the five entries and their tests. `validate_group` then fails on the missing filament owner and on the empty internal set. Nothing outside the module depends on them. |

## 14. Honest limits

- **Unratified**, and §3 is the part that needs a human.
- **No bundle mechanics exists in Aleph.** Nothing here is evidence that Aleph represents a stress
  fibre, an arc, or a cap.
- **Every stiffness is `UNSOURCED`** and no literature was read. Four material cards are named as owed
  and none is supplied.
- **The `SEAMED` finding is `AUDIT_READ`.** I read the status table and verified that no runtime class
  exists for the three names by searching that tree. I did **not** run its suite or attempt to
  construct the connectors, so "no concrete runtime" rests on a search plus its own audit's
  concurrence, not on an execution.
- **The claim that disjointness makes coupling harder is an argument, not a measurement.** It is
  strongly suggested by which connectors got bound and which did not, and it has an obvious
  alternative explanation — the three `SEAMED` ones also need owners that do not exist yet
  (intermediate filament, microtubule), so scheduling alone could explain them. `sf_cortex_transient`
  is the case that discriminates, because the cortex *does* exist there, and it is still `SEAMED`.
  That is one data point, and I am reporting it as one.
- **Nothing here addresses whether four subpopulations are the right four.** The registry says four; I
  registered four. Whether the biology needs a fifth is not a question this entry asks.
