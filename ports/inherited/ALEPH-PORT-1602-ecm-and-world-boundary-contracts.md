# ALEPH-PORT-1602 — The `ecm` and `world_boundary` registry entries, and why one of them owns nothing

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1602` |
| Lane | `L16 registry data — external environment and cell surface` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Nothing was ported. One structural defect in the reference registry is recorded, and the Aleph entry is shaped so the same defect cannot recur. |

---

## 1. Aleph API

The registry entries named `ecm` and `world_boundary`. They are filed together because the second
exists to terminate the first: a finite fibre network with no far-field anchor is not a stiff network,
it is an unconstrained one, and neither entry means anything without the other.

```python
from aleph.state.census_environment_surface import (
    EXTERNAL_ENVIRONMENT,
    ScopeTag,
    contract_for,
    names_at_scope,
)
```

Aleph target file: `aleph/state/census_environment_surface.py`.

Declared data only: scope tag, interpretation, numerical representation, owned-state names, mechanical
role, unsupported claims, re-entry condition. **No constitutive law, no crosslink kinetics, no contact
law and no solver** is authorised by this entry.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/ecm_world.py`, `ffn_sim/ac/engine/contracts.py` |
| Source symbol(s) | `ECMWorldSettings`, `BoundaryAnchorMode`, `ECMInternalCrosslink`, `CompositeECMClutch`, `ECMMembraneContact`, `ECMBoundaryAnchorFacade`, `FarFieldDirichletBoundaryRuntime`, `MembraneContactSurfaceView`, `_far_field_reaction_ledger_kernel`, `_far_field_reaction_commit_kernel`; and the `ComponentContract` declarations of `ecm` and `world_boundary` inside `reference_cell_architecture` |
| Read from | **`git show be0e5876:<path>`.** The working tree was not read for content. |
| Working tree == commit? | **yes** for both files, verified by digest: `ffn_sim/ac/engine/ecm_world.py` = `sha256:e05a3c787bdb33667c8d749cf57cf477fc3a590730fe092e46fc64a49eb24d63`, `ffn_sim/ac/engine/contracts.py` = `sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72`. |
| Working-tree caveat | The tree carries **31 uncommitted changes**, and one file this lane audited (`ffn_sim/ac/engine/cortex_state.py`) does differ from the commit — see `ALEPH-PORT-1604` §2. Per-file digests were therefore computed rather than assumed. |

## 3. Why source-derived porting beats clean-room

**It does not.** RE-DERIVE.

The unit of approval is an enumeration of what a collagen network owns. Any author who sits down to
list the state of an explicit crosslinked fibre network arrives at node positions, material
coordinates, crosslink topology, ligand coordinates, damage and remodelling state — this is a
convention any competent author reaches independently, which `ports/TEMPLATE.md` §3 names explicitly as
**not** a reason to port.

Two things in the source are better than a naive enumeration and are worth recording as *reasons*
rather than as code:

1. **The adhesion endpoint is a barycentric material point on a segment, not a captured node.** That
   matters because it lets an endpoint slide along a fibre and be remapped after a topology change
   without inventing or destroying load. A clean-room author might well have captured a node, which is
   simpler and wrong in a way that only shows up after remeshing.
2. **A topology epoch.** Any cached view of the network taken before the epoch changes is *stale*
   rather than merely old, so the distinction has to be representable. Again an argument, not code.

Both are recorded in the Aleph entry's prose in Aleph's own words, and both are re-derivable from the
requirement that a connector's load survive a topology change — which is where §4 derives them from.

## 4. Physical or mathematical law represented

**No law crosses.** What crosses is an ownership argument, in two parts.

**Why the ECM owns state and the boundary does not.** The ECM's fibres have positions that evolve, a
crosslink topology that changes, and damage that accumulates: it is an integrand of the same
transaction as the cell. The world boundary has a pose that never changes. It has no degrees of
freedom, no constitutive law, nothing to snapshot and nothing to restore. The reaction it receives is
computed by the anchor connector and reported; it is not state belonging to the frame.

This is not bookkeeping pedantry. A registry that lists a reference frame among its state owners has
miscounted its own breadth, and the miscount is invisible because the row looks like every other row.
Aleph's `ScopeTag.BOUNDARY` plus the rule that `H`/`B`/`X` own no state is exactly the machinery that
makes the miscount fail a test.

**Why the anchor is not optional.** Let the modelled matrix volume carry a net traction `T` from the
cell's adhesions. With no anchor the network's centre of mass accelerates freely: the momentum balance
is satisfied by translation rather than by strain, and the compliance the cell's adhesions then measure
is the compliance of a free-floating object — dominated by the domain's size and by nothing physical.
So a finite domain with no far-field constraint does not model a soft substrate; it models no
substrate. The reaction the boundary reports is therefore a **diagnostic of domain truncation**, and
that is what the Aleph entry registers it as, rather than as a tissue-scale force measurement.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `fibre_node_positions_um` | µm | m | finite |
| `fibre_material_coordinates` | µm (arc length) | m | ≥ 0, monotone along a fibre |
| `crosslink_topology` | dimensionless (index sets) | — | valid indices only |
| `topology_epoch` | dimensionless (counter) | — | non-decreasing integer |
| `fibre_damage_state` | dimensionless | — | in `[0, 1]` per element |
| `world_boundary` owned state | — | — | **empty by construction** |

Singular and boundary cases the entries commit to:

- **An unanchored finite matrix volume is not a physiological configuration.** It is refused, not
  approximated by a large stiffness. The registry records this as the boundary's mechanical role rather
  than as a settings default.
- **A topology change invalidates cached views.** After the epoch advances, any view taken earlier is
  stale, and stale is a refusal condition rather than a tolerance.
- **Damage and remodelling commit only after an accepted step.** An irreversible change made during a
  candidate step that is then rejected would leave the network in a state no accepted trajectory
  passed through.

Invariants, each with the control that asserts it:

- **I1.** `ecm` is at `E` and owns non-empty state, including `crosslink_topology`, `topology_epoch` and
  `fibre_damage_state`.
  `tests/state/test_census_environment_surface.py::test_the_ecm_owns_topology_and_an_epoch`.
- **I2.** ~~`world_boundary` is at `B` and owns **exactly nothing**.~~
  **SUPERSEDED 2026-08-04 by PI decision** (`c529fd3`). The frame is at `E` and owns five keys —
  a far-field displacement, a strain schedule, a surround compliance, an applied traction and the
  reaction it already received. The ruling: **the registry lists engine parts, not cell parts**, the
  external environment varies and applies forces, and `ecm` and `extracellular_medium` were already
  `E` while sitting outside the cell. This entry's own `reentry_condition` had named the trigger —
  *"an experiment that loads the far field ... needs a boundary that evolves"* — and it was met.
  `tests/state/test_census_environment_surface.py::test_the_boundary_owns_state_and_evolves`.
- **I3.** The boundary receives the far-field reaction and says so.
  `tests/state/test_census_environment_surface.py::test_the_boundary_receives_the_far_field_reaction`.
- **I4.** The boundary states **why its scope follows its state**, so the rule rather than the
  instance is what a reader meets. ~~why it is *not* counted as a state owner~~ — superseded with
  I2 above.
  `tests/state/test_census_environment_surface.py::test_the_boundary_says_why_its_scope_follows_its_state`.
- **I5.** `ecm` declares that it has no interstitial fluid and therefore no poroelastic response, and
  points at the entry that owns that exclusion.
  `tests/state/test_census_environment_surface.py::test_the_ecm_excludes_the_interstitial_fluid_it_is_paired_with`.

## 6. Source evidence class and known retractions

Read at `be0e5876`; not executed — `ecm_world.py` requires Warp and this session ran **zero GPU jobs**
(PLAN §0.1). Evidence class `AUDIT_READ` at a stated commit. No numerical claim about the source is made.

Retractions looked for in: both module docstrings, `PLAN.md` §1.1 and §7, and `ports/audit/`. The
source's own settings object refuses to build a production ECM without a physiological anchor, without
crosslink kinetics, and without remodelling and damage state — a set of refusals rather than defaults,
which is the right shape and is noted as such.

**Finding — the registry counts a boundary condition among its state owners.** In
`reference_cell_architecture`, `world_boundary` is declared as a `ComponentContract` with
`owns_geometry=False, dynamically_evolving=False`. A `ComponentContract` in that registry *is* the
declaration of a state-owning component; the registry's own header calls its members state owners and
counts fourteen of them. So the headline count includes at least one member that owns no geometry and
never evolves — a reference frame filed as an owner.

Combined with `ALEPH-PORT-1601` §6, the same file also declares `extracellular_medium` with
`owns_geometry=False, dynamically_evolving=False` while its implementation demonstrably owns and
commits state. So the two flags are doing two different jobs in the same table: for `world_boundary`
they are true and mean "this is not an owner", and for `extracellular_medium` they are false in fact.
A reader cannot tell which meaning applies to a given row without opening the implementation.

The root cause is structural and is worth naming precisely, because it is the reason this lane exists:
**that registry has no scope tag.** Its `ComponentContract` fields are `name`, `role`,
`representation`, `solver`, `owns_geometry`, `dynamically_evolving`, `has_events`. There is no `E/I/H/B/X`
field, no owned-state list, no approximation field, no unsupported-claims list and no re-entry
condition. With no way to say "boundary", a boundary has to be filed as a component with its flags
turned off — and with no way to say "homogenized" or "excluded", the six non-explicit entries of this
lane's ten have no representation in that registry at all. They live in prose or nowhere.

Aleph's response: `world_boundary` is registered at `B`, `validate_contract` refuses any `B` entry that
lists owned state, and `test_a_boundary_entry_that_owns_state_is_refused` fails if that refusal stops
working.

## 7. Independent oracle or derivation

For the ECM entry there is no numerical oracle, because the entry computes nothing — the honest answer,
stated rather than dressed up.

What is independently derivable, and is derived in §4 without the source open, is the **necessity of
the anchor**: a momentum-balance argument shows an unanchored finite domain translates rather than
strains, so the boundary is a precondition for traction being a meaningful measurement rather than a
modelling nicety. That derivation is Aleph's own and needs no reference to check.

For a future implementation the available oracles are: linear-elastic effective-modulus limits for an
affine fibre network, exact zero net force on an unloaded network, and the identity that the anchor
reaction must equal minus the total cell-applied traction at equilibrium. **None is exercised here.**

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_environment_surface.py::test_the_ecm_owns_topology_and_an_epoch` | `ecm` owns crosslink topology, a topology epoch and damage state, and its role states that topology changes commit only after an accepted step. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_boundary_owns_state_and_evolves` | `world_boundary` is at `E`, owns non-empty state and declares itself evolving. **Replaced its own negation on 2026-08-04**, per `c529fd3`. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_boundary_carries_every_way_the_surround_can_act` | The five state keys are asserted **by name**, not by count — a count passes when one is swapped for another, which is how `applied_traction` would quietly go missing at five. |
| Positive | `tests/state/test_census_environment_surface.py::test_it_still_refuses_the_claims_a_lumped_compliance_cannot_support` | Giving the frame state did not give it structure: tissue-scale force and spatially resolved surround structure remain unsupported, and the re-entry condition still names a network. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_boundary_receives_the_far_field_reaction` | The far-field reaction is named as what the boundary receives. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_ecm_excludes_the_interstitial_fluid_it_is_paired_with` | `ecm` declares no poroelasticity and names the entry that owns the exclusion. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_interstitial_fluid_records_both_halves_of_its_reassignment` | The interstitial fluid's `H` entry names both destinations of its absorbed contribution and states that no documented homogenization map exists, so the absorption is a declared reassignment rather than a calibration. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_a_boundary_entry_that_owns_state_is_refused` | A `B` entry listing owned state is rejected. This is the reference's exact defect, planted deliberately: the control fails if Aleph ever lets a reference frame be registered as an owner. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_a_homogenized_entry_that_owns_state_is_refused` | An `H` entry listing owned state is rejected, which is the same conflation for the interstitial fluid and the substrate. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_a_homogenized_entry_with_neither_approximation_nor_reentry_is_refused` | An `H` entry that declares neither what is lumped nor what would bring it back is rejected — an omission dressed as a decision. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_a_duplicate_name_is_refused_at_registry_level` | Two entries under one name are rejected, so a manifest cannot describe one compartment twice with different scope. |

## 10. Numerical and precision envelope

No arithmetic; no tolerance introduced; nothing here may be quoted as a numerical result. The entries
sit at the evidence class of a declaration and this field records that plainly rather than borrowing an
envelope from an implementation that does not exist in Aleph.

The one numerical *requirement* the entries record for later work: `topology_epoch` must be an exact
integer comparison, never a tolerance. "The view is one epoch stale" and "the view is current" are
discrete facts, and a float epoch with a comparison tolerance would make staleness a matter of degree,
which it is not.

## 11. Production-backend residency and transfer

Host-side declarative data in `aleph/state/census_environment_surface.py`: frozen dataclasses, no
array, no device allocation, no kernel, no per-step transfer. It never runs on the production backend
and must not import one — the census is what a backend is configured from, so the dependency must not
run the other way.

Imports are standard library and `aleph` only, asserted by
`test_the_module_imports_nothing_but_the_standard_library_and_aleph`. `aleph/state/**` may not import
`validation/**`.

## 12. Comments and docstrings to discard

No source prose survives. Discarded and replaced rather than merely deleted:

- provider module paths, decision identifiers, plan section references, track labels and severity
  labels — replaced by the physical statement of what is owned and why;
- the source's `Sanity Gate` docstring convention — replaced by the invariant list in §5, each item
  naming a test that runs;
- the source's component-role vocabulary, which cannot express `B` and therefore forces a boundary to
  be declared as a component with its flags off — replaced by `ScopeTag.BOUNDARY` plus the empty-state
  rule;
- the source's settings-object switch names. Aleph records the *refusals* those switches encode as
  scope statements in prose, not as a mirrored configuration surface.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`. The data has landed and 83 controls in `tests/state/test_census_environment_surface.py` pass, but the registered content is a modelling declaration no human has reviewed, The contract it registers into has since landed and this data now constructs against it, so the remaining blocker is review rather than a missing dependency. |
| Reviewer | agent-proposed, **unratified**. |
| Rollback | Delete the `ecm` and `world_boundary` entries from `EXTERNAL_ENVIRONMENT`. Breaks the controls in §8 and the `substrate_basement_membrane` entry, which is defined as the combination of these two plus the far-field anchor and would then reference entries that do not exist. |

## 14. Honest limits

- **Unratified**, as above.
- Written against a state-contract type that **was not on disk when this lane ran**; the module imports
  defensively and records which definition is live in `CONTRACT_SOURCE`.
- The reference was **read, never run.** No GPU authorization exists. Everything in §6 is a reading of
  two files at one commit.
- The ECM owned-state list is a declaration. **No Aleph ECM implementation exists**, so the list is
  entirely unexercised and an implementation may find it insufficient — in particular it says nothing
  about how remeshing or refinement state is represented, which the reference does carry and which this
  entry deliberately does not enumerate for lack of a derivation Aleph owns.
- The ECM entry says nothing about the *constitutive law* of a collagen fibre. That is out of V1 scope
  and is not merely unimplemented but undecided.
- "The reaction is a diagnostic of domain truncation" is an argument, not a measurement. How large a
  domain suffices for a given adhesion pattern is unquantified here and would need a convergence study.
- `world_boundary`'s rigidity is a modelling choice with no supporting measurement, carried at
  `UNSOURCED`. Whether a real substrate's surround is well approximated by a rigid frame is an open
  question this entry does not address; it only makes the choice visible.
