# ALEPH-PORT-1701 — lamellipodium: branched protrusive actin as a census contract

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1701` |
| Lane | `L17 registry — actomyosin load-path compartments` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` — **verdict: RE-DERIVE.** No code, no identifier, no constant crosses. |

---

## 1. Aleph API

```python
from aleph.state.census_actomyosin import (
    LAMELLIPODIUM,
    DISJOINT_FILAMENT_POPULATION_RULE,
    POPULATION_COUNTS_ARE_OUTPUTS,
    validate_contract,
)
```

Aleph target file: `aleph/state/census_actomyosin.py`. This entry authorises the `lamellipodium`
contract entry and nothing else. No mechanics, no kernel, no rate law is authorised here, because none
was written.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/protrusion.py`, `ffn_sim/ac/weave/lamellipodium.py`, `ffn_sim/ac/weave/membrane_ratchet.py`, `ffn_sim/ac/weave/branch_angle.py`, `ffn_sim/ac/weave/woven_cell.py` |
| Source symbol(s) | `LamellipodiumBranchAngleMechanics`, `LamellipodiumStateOwner`, `build_lamellipodium`, `AdaptiveLamellipodium`, `ratchet_velocity`, `free_polymerization_velocity`, `stall_force`, `angle_energy`, `angle_forces`, `WovenCell` |
| Read from | **`git show be0e5876:<path>`** — the immutable commit, not the working tree. |
| Working tree == commit? | `yes` for every path above — verified with `git diff be0e5876 -- <path>`, empty in all five cases. |

File digests at `be0e5876` (first 16 hex of `sha256` over the blob), recorded so a later reader can
tell whether they are looking at what was read:

| Path | `sha256:` | Lines |
|---|---|---|
| `ffn_sim/ac/engine/protrusion.py` | `sha256:98300c76d61b3a6d` | 1107 |
| `ffn_sim/ac/weave/lamellipodium.py` | `sha256:050d2f910598d436` | 360 |
| `ffn_sim/ac/weave/membrane_ratchet.py` | `sha256:6417399c82bad83f` | 158 |
| `ffn_sim/ac/weave/woven_cell.py` | `sha256:25a61890595df7b0` | 477 |

**Correction to a premise this audit was given.** `be0e5876` **is `HEAD`** in that repository, so the
"31 uncommitted changes" are exactly the working-tree delta against it: 22 modified tracked files plus
9 untracked. None of them is a file cited above. The premise as handed over implied `be0e5876` was an
older commit the tree had moved past; it is not.

## 3. Why source-derived porting beats clean-room

**It does not.** Verdict: **RE-DERIVE**, and the reason is structural rather than a matter of taste.

The candidate assets split into three groups and none of the three is portable:

1. **The engine-side owner and mechanics delegate** (`protrusion.py:303-511`) is well shaped and
   **never runs.** It is constructed only in that repository's own tests; the two references to it on
   the driver path are comments. Its test module carries 45 fake/spy doubles and allocates no device
   array at all, so the delegate has never been code-generated for a backend. There is nothing to
   inherit but a class shape, and a class shape is a convention any competent author reaches
   independently — `ports/TEMPLATE.md` §3 names that as an illegitimate reason.
2. **The branch-angle kernel and its analytic layer** (`branch_angle.py:67-183`) is real numeric code
   and the physics is correct as far as it goes — an angular harmonic with a `k_θ ↔ σ_θ`
   equipartition round-trip. But its own audit records the branch-angle citation as `DOI_DEAD`, and
   the whole point of a rate prior is its source. Porting a formula whose only parameter has a dead
   citation would inherit the unverifiability along with the algebra.
3. **The Brownian-ratchet host code** (`membrane_ratchet.py:68-158`) is correct NumPy and is the one
   asset here I would have expected to port. It is out of Aleph's current scope: no protrusion
   vertical exists, and a ratchet with no membrane to push is a function with no oracle.

Nothing crossed. This entry is a record of **having looked**, which is a different fact from nobody
having looked, and PLAN §0.2.5 wants both facts distinguishable.

## 4. Physical or mathematical law represented

No law was ported. What this entry authorises is an **ontological law** — a statement about ownership
that is checkable — and it is derived here rather than quoted:

**One owner per physical filament.** Let `F` be the set of physical filaments in the model and let
`O` be the set of components. An ownership assignment is a function `o : F → O`. Being a *function*
is the entire content of the rule: a filament has exactly one owner, so `o` is single-valued, and the
populations `o⁻¹(c)` for `c ∈ O` therefore partition `F`. Disjointness is not an extra axiom; it is
what "function" means. What can go wrong is that `o` is never constructed, and a *relation* is used
instead — which is what a shared array with a region label is.

The consequence that makes it a mechanical statement rather than a philosophical one: if two declared
components are slices of one array, then two nodes in different "populations" are already neighbours
in one index space, so any pairwise interaction evaluated over that array transmits load between them
**whether or not a connector was declared**. Force closure cannot detect this, because the merged
network's forces really do balance. Only the ownership map can, and only before allocation.

## 5. Units, domains, singular cases, invariants

This entry authorises declarations, not quantities, so the units table records what the *declared
state keys would* carry when a later lane allocates them. Nothing here holds a number.

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| node position | µm | m | finite |
| filament count | count | 1 | ≥ 0, an output |
| barbed-end elongation velocity | µm/s | m/s | sign set by load |
| branch angle | 1 (radian) | rad | (0, π) |
| branch-nucleation rate | 1/s | 1/s | ≥ 0, **unsourced** |

Singular cases, each with the behaviour Aleph requires:

- **Zero filaments.** A legal state, not an error: a lamellipodium that has not nucleated yet owns an
  empty population. Any consumer that divides by the filament count must refuse rather than return a
  large number.
- **A count requested as an input.** Refused. The count is an output; see §4 of
  `aleph/state/census_actomyosin.py`'s module docstring and `POPULATION_COUNTS_ARE_OUTPUTS`.
- **A rate requested with no source.** Refused, not defaulted. Every rate in the table above is
  `UNSOURCED` and the contract says so.

Invariants that must hold, each with the test that asserts it:

- **I1.** The entry carries the disjointness rule verbatim in `mechanical_role`.
  `tests/state/test_census_actomyosin.py::test_each_filament_owner_records_the_disjointness_rule`.
- **I2.** No state key it owns is claimed by another owner.
  `tests/state/test_census_actomyosin.py::test_no_state_key_is_claimed_by_two_owners`.
- **I3.** It owns a named filament population block, prefixed with its own name, so a collision is
  detectable. `test_each_filament_owner_names_its_own_filament_population_block`.
- **I4.** It records that its count is an output and names the priors evidence is owed on.
  `test_each_dynamic_population_records_that_its_count_is_an_output`.
- **I5.** It declares `implemented=False` and `citation_status="UNSOURCED"`.
  `test_nothing_in_this_group_claims_to_be_implemented`.

## 6. Source evidence class and known retractions

I looked, and there is a great deal to report.

That repository's own audit at `ffn_sim/docs/v2_audit/AUDIT_AC_ENGINE_2026-07-25.md:37` rates
`lamellipodium` **`KERNEL_BOUND`, zero CUDA execution, "declared-only in practice"**, and at `:89`
lists the lamellipodium among components that were **never executed**. Its compartment-validation
track (`COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md:227`) records that the branch-angle harmonic *is*
genuinely native — but on the cortex path, under a cortex owner, never under a lamellipodium owner.
`ROLLING_ROADMAP.md:279-280` still lists lamellipodial branching and the Brownian ratchet as work to
bind. The composed-world artifact self-labels `evidence: CENSUS-WIRED`,
`quantitative_claim_status: BLOCKED`.

Reachability: the engine owner and mechanics delegate are **dead code** outside tests. The branch
kernel is live, but under a different owner. The ratchet runs in no device runtime at all.

Tests in the source: `test_lamellipodium_oracle.py` and `test_adaptive_lamellipodium.py` are real
physics (branch-angle distribution against θ₀/σ, nucleation-capping balance, protrusion slowing under
load). `test_protrusion.py` is 726 lines of plumbing with 45 doubles. So the answer to "does it have
tests" and the answer to "do the tests exercise the law" are different, which is exactly why the
template asks them separately.

Retraction found: the branch-angle citation is recorded **`DOI_DEAD`** by that repository's own audit
(`:160`). That is a retraction in substance — the parameter's provenance does not resolve — and it is
the single strongest reason not to port the angular term with its constant attached.

## 7. Independent oracle or derivation

The ownership law in §4 is checked by an oracle Aleph owns and that needs no reference at all: the
assignment `o : F → O` is single-valued **iff** no state key appears under two owners. That is an
exact set-theoretic identity, so the control is an equality on a name collision and not a tolerance.
`assert_state_keys_are_disjoint` implements it, and its negative control constructs the collision
deliberately.

There is no numerical oracle in this entry because there is no number in it. When a lamellipodial
mechanics lane exists, the oracle available to it is the free-polymerization/stall-force limit pair:
protrusion velocity must approach the load-free elongation rate as load → 0 and reach exactly zero at
the stall force, with the exponential interpolation between them derived from the ratchet's own
Boltzmann factor rather than fitted. That is a two-endpoint exact check and it is Aleph's to write.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_actomyosin.py::test_each_filament_owner_records_the_disjointness_rule` | The disjointness rule is present verbatim in `mechanical_role`. |
| Positive | `tests/state/test_census_actomyosin.py::test_no_state_key_is_claimed_by_two_owners` | All 35 declared state keys across the group have exactly one owner. |
| Positive | `tests/state/test_census_actomyosin.py::test_each_dynamic_population_records_that_its_count_is_an_output` | The count rule and at least three named rate priors are on the entry. |
| Positive | `tests/state/test_census_actomyosin.py::test_registering_the_group_into_a_manifest_yields_ten_unimplemented_components` | Registration works and the manifest reports ten declared, zero implemented. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_a_filament_owner_that_drops_the_disjointness_rule_is_refused` | Replacing `lamellipodium.mechanical_role` with plausible prose that omits the rule is refused. If this stops firing, I1 is unguarded. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_two_owners_sharing_a_state_key_is_refused` | Giving `filopodium` the lamellipodial filament-population key — the declarative form of one shared array — is refused by name. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_dropping_the_count_is_an_output_rule_is_refused` | Removing the count rule from an entry is refused, so the rule is not decoration. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_an_explicit_owner_that_owns_nothing_is_refused` | An `E` entry with an empty `owned_state` is refused: a participant with no state cannot carry a load path. |

## 10. Numerical and precision envelope

No floating-point computation crosses this boundary, so there is no working precision and no
tolerance to justify: every control here is an exact string or set comparison, asserted with `==` and
`in`, and a tolerance would be meaningless. That is a real answer to this field rather than a dodge —
the check is exact because the property is discrete, and making it approximate would be the defect.

The one number that will matter downstream is recorded so a later lane cannot inherit it accidentally:
the branch angle and its spread, the nucleation rate, the capping rate and the elongation rate are all
absent. There is no default, no placeholder value, and no "provisional" constant anywhere in
`aleph/state/census_actomyosin.py`. A consumer that needs one must refuse.

## 11. Production-backend residency and transfer

This module is a registry of declarations. It never runs on a production backend, holds no device
array, and performs no host transfer — `owned_state` is a tuple of `str`. It is host-side Python, CPU
by definition, and that is a residency answer: it belongs on neither side of the accelerator boundary
because it is not in the step loop at all.

`aleph/state/**` must not import `validation/**`; this module imports only `dataclasses`, `enum`, and
`aleph.state.schema`.

## 12. Comments and docstrings to discard

Read, and deliberately **not** carried:

- The provider's class, owner, delegate and kernel names, and its module paths.
- Its evidence-rung vocabulary and its gate labels. Aleph has its own two-axis evidence ladder; a
  second rung vocabulary would let a reader think an Aleph entry had been graded on that scale.
- Its status strings and roadmap identifiers.
- Its branch-angle constant and the citation attached to it, since that citation is recorded dead by
  its own audit.
- Its region-label vocabulary, which is the thing this entry exists to refuse.

What replaces them: §4 above, restated as an ownership law with a set-theoretic derivation, and the
module docstring of `aleph/state/census_actomyosin.py`, which describes the failure mode in terms
anyone can check against Aleph's own code without that repository on disk.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`. The contract entry and its controls exist and pass — `tests/state/test_census_actomyosin.py`, 88 passed — but this entry stays `PROPOSED` because the verdict it records (RE-DERIVE, nothing ported) is an agent judgement about another repository's assets and no human has checked it. |
| Reviewer | agent-proposed, **unratified**. |
| Rollback | Delete the `LAMELLIPODIUM` entry from `aleph/state/census_actomyosin.py` and its tests. `validate_group` then fails on the missing filament owner, which is the correct loud failure. Nothing outside the module depends on it yet. |

## 14. Honest limits

- **Unratified.** The RE-DERIVE verdict is mine, not the PI's.
- **A declaration, not a model.** There is no lamellipodial mechanics in Aleph. This entry cannot be
  cited as evidence that Aleph represents a lamellipodium; it is evidence that Aleph has *declared*
  what a lamellipodium would own.
- **Every rate is `UNSOURCED`** and no literature was read for this entry. Four priors are named as
  owed and none is supplied.
- **The reachability findings are that repository's own self-report plus my reading**, not my
  execution. I did not run its suite, did not construct its owners, and did not launch a kernel. So
  "never executed" is `AUDIT_READ`, corroborated by its audit trail, and not `MEASURED`.
- **The ratchet was judged out of scope rather than judged wrong.** If a protrusion vertical is
  authorised, `membrane_ratchet.py` deserves a fresh audit against a real oracle rather than
  inheriting this entry's verdict.
- **The disjointness law is checked on names, and names can be gamed.** Two owners that both allocate
  from one buffer while declaring distinct keys would pass. Closing that needs an allocation-time
  check, which needs an allocator, which does not exist. Recorded rather than papered over.
