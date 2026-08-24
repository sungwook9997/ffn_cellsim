# ALEPH-PORT-1601 — The `extracellular_medium` registry entry: ownership, state, and scope limits

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1601` |
| Lane | `L16 registry data — external environment and cell surface` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** No code and no prose crosses. One structural fact about the reference is recorded as a *contradiction to avoid*, in §6. |

---

## 1. Aleph API

The registry entry named `extracellular_medium`, and nothing else in the module.

```python
from aleph.state.census_environment_surface import (
    EXTERNAL_ENVIRONMENT,
    contract_for,
    scope_of,
)
```

Aleph target file: `aleph/state/census_environment_surface.py`.

This entry covers **declared data**: a scope tag, a biological interpretation, a numerical
representation, the names of the state the component owns, its mechanical role, and the list of claims
it does not support. It authorises **no solver, no operator and no kernel.** An exterior Stokes
implementation is a separate port with its own entry and its own controls, and this entry must not be
cited as evidence for one.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/medium_exterior.py`, `ffn_sim/ac/engine/contracts.py` |
| Source symbol(s) | `ExteriorStokesMedium`, `ExteriorStokesMediumSettings`, `MediumWallMode`, `assert_single_dissipation_owner`, `derive_blob_epsilon`, `build_exterior_stokes_medium`, `_medium_commit_kernel`, `_medium_rollback_kernel`, `_surface_velocity_kernel`; and the `ComponentContract` declaration of `extracellular_medium` inside `reference_cell_architecture` |
| Read from | **`git show be0e5876:<path>`.** The working tree was not read for content; it was compared, below. |
| Working tree == commit? | **yes** for both files. Verified by file digest, not by `git status`: `ffn_sim/ac/engine/medium_exterior.py` = `sha256:fe04594c777bf6c6a83fce62d8681f2b16e027f90ce0333a3fe55ab6110980f1` and `ffn_sim/ac/engine/contracts.py` = `sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72` at the commit, both matching the tree. |
| Working-tree caveat | That tree carries **31 uncommitted changes** overall and at least one file audited by this lane (`ffn_sim/ac/engine/cortex_state.py`, see `ALEPH-PORT-1604`) differs from the commit. So "the tree matches" is a per-file measurement here, not an inherited assumption. |

## 3. Why source-derived porting beats clean-room

**It does not.** RE-DERIVE, for three independent reasons, none of which is a criticism of the source's
craft:

1. **Nothing in the unit of approval is code.** This entry registers a declaration. There is no
   branch structure to get subtly wrong and no conditioned expression to inherit.
2. **The reference's registry cannot express what Aleph's must.** Its `ComponentContract` carries
   `name`, `role`, `representation`, `solver`, `owns_geometry`, `dynamically_evolving`, `has_events` —
   and **no scope tag, no owned-state list, no approximation field, no unsupported-claims list and no
   re-entry condition.** It can register explicit state owners and nothing else. There is therefore no
   source shape to port for six of this lane's ten entries, because the source has no way to say
   "homogenized", "boundary" or "excluded" in the registry at all; those statements live in prose
   elsewhere or nowhere.
3. **The operator itself is out of V1 scope and off-limits to run.** The exterior Stokes resistance is
   Warp/CUDA-resident by construction, and this session holds **no GPU authorization** (PLAN §0.1). An
   asset that cannot be executed cannot be audited by running it, and a read-only audit is not
   evidence about numerical behaviour.

The `assert_single_dissipation_owner` *pattern* — one named owner per dissipation channel, enforced as
a guard rather than as a deletion — is genuinely good, and Aleph already arrived at the same shape
independently for the pressure path (`ALEPH-PORT-1102`, `assert_single_pressure_path`). Convergent
design is prior art, not a reason to port; PLAN §0.2 and `ports/TEMPLATE.md` §3 both say so.

## 4. Physical or mathematical law represented

There is **no law in this entry.** What crosses into Aleph is an *ownership declaration*, and saying
so plainly is the honest answer to this field.

The declaration rests on two derivations Aleph owns, and both are recorded in the entry's prose so a
reader meets them where they matter:

**Why an exterior medium and not a numerical regulariser.** A closed cell has six rigid-body modes
that carry exactly zero internal elastic energy: translating or rotating the whole assembly changes no
bond length, no area and no curvature. So the internal stiffness matrix is singular on a
six-dimensional subspace, and every scheme must do something about it. Adding `a·I` with `a` set from
the step-size bound is a mobility the integrator already computed — it has no dependence on the cell's
size or shape, so no measurement can constrain it and no experiment can be wrong about it. An exterior
fluid resistance is the physically correct occupant of that slot: it is a genuine force law, it scales
with the cell's geometry, and it has closed-form limits for a sphere (translational resistance
`6·π·μ·a`, rotational `8·π·μ·a³`) that an implementation can be checked against.

**Why Stokes and not Darcy.** Outside the cell there is no solid fraction, so the Brinkman screening
length `√k` diverges. Darcy is not an approximation of the exterior problem; it is a description of a
different problem. The registry records that as a scope statement rather than a parameter choice.

**Why the medium owns state.** Its force is `f = −M⁻¹v` with `v = (x_candidate − x_committed)/dt`. A
velocity is a difference of two configurations, and one of them has to be remembered across a step.
That memory *is* state, whatever a declaration says about it.

## 5. Units, domains, singular cases, invariants

The entry is data, so its units live in the names of the owned state. That is deliberate: a state
name carrying its unit cannot be silently reinterpreted.

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `committed_surface_reference_um` | µm | m | finite |
| `surface_traction_pn` | pN | N | finite |
| `cumulative_dissipated_work_pn_um` | pN·µm | J | ≤ 0 for any real motion |
| medium viscosity (a context parameter, not owned state) | pN·s/µm² | Pa·s | > 0 |

`1 pN/µm² = 1 Pa` exactly, which makes the resistance's dimensional check trivial and is worth
remembering rather than re-deriving.

Singular and boundary cases the registry entry commits to:

- **Zero surface velocity gives exactly zero medium force.** Not small — zero. A resistance that
  invents a force at rest is a spurious source, and the exactness is testable because the branch is
  identically zero rather than nearly zero.
- **No wall.** A half-space image system is exact and closed-form and is *absent*, so any
  near-substrate claim is out of scope. Registering "free space" silently while a substrate is present
  would be the worst available outcome, so the entry names the exclusion.
- **Sign.** Hydrodynamic power is strictly negative for any non-zero surface velocity field. A medium
  that adds energy is a sign error, never a physiological result.

Invariants, each with the control that asserts it:

- **I1.** The medium is registered at `E` and owns non-empty state.
  `tests/state/test_census_environment_surface.py::test_every_explicit_entry_owns_state`.
- **I2.** The committed surface reference is among its owned state, so nothing downstream can treat
  the medium as non-evolving.
  `tests/state/test_census_environment_surface.py::test_the_medium_owns_a_committed_surface_reference`.
- **I3.** It is declared the sole owner of velocity-proportional surface dissipation.
  `tests/state/test_census_environment_surface.py::test_the_medium_is_the_sole_owner_of_surface_dissipation`.
- **I4.** The absent wall is declared, with the consequence named.
  `tests/state/test_census_environment_surface.py::test_the_medium_excludes_a_wall_rather_than_defaulting_to_free_space`.

## 6. Source evidence class and known retractions

The two files were read at `be0e5876`; neither was executed, because both require Warp on CUDA and
this session ran **zero GPU jobs** (PLAN §0.1). So the evidence class for everything in this section is
`AUDIT_READ` at a stated commit, and no numerical claim about the source is made or implied.

Looked for retractions in: the module docstrings of both files, `PLAN.md` §1.1 and §7, and
`ports/audit/`. What was found:

- The source's own docstring records that `mu_medium` has no knowledge-base claim and no source
  evidence, and is therefore required from the caller with no default. That is the correct handling and
  Aleph inherits the *policy*, not the parameter: this entry's `citation_status` is `UNSOURCED`.
- The component's declaration carries a note that declaring it does not implement it, and that the
  rigid modes stay numerically regularised until an implementation lands. The declaration is honest
  about its own emptiness, which is unusual and worth recording.

**Finding — the contract and its implementation disagree about whether this component evolves.**
`reference_cell_architecture` declares `extracellular_medium` with `owns_geometry=False,
dynamically_evolving=False`. Its implementation `ExteriorStokesMedium` owns `committed_position_d`,
`committed_traction_d` and a `dissipated_work` accumulator, exposes `snapshot_candidate` and a
reject-gated `_medium_rollback_kernel`, and advances the committed reference and traction inside
`commit_irreversible` via `_medium_commit_kernel` under the accepted predicate. Snapshot, rollback and
commit under an acceptance predicate is the definition of evolving under the transaction.

The direction of the disagreement is the dangerous one. A consumer that trusts
`dynamically_evolving=False` has no reason to snapshot the medium, so after a **rejected** step the
surface rolls back while the medium's committed reference does not — and the retried step then measures
its velocity against the wrong reference. The result is a wrong drag, of plausible magnitude, on a
retried step only. Nothing crashes and no gate on a force norm would see it.

Aleph's response is structural, not a fix to the reference: the committed surface reference is
registered as **owned state**, and `test_the_medium_owns_a_committed_surface_reference` fails if it
ever stops being.

**Finding — the reference registry counts a non-owner as an owner.** See `ALEPH-PORT-1602` §6 for
`world_boundary`, which is the same conflation in its purest form.

## 7. Independent oracle or derivation

For an implementation, the oracles are the closed-form sphere limits `6·π·μ·a` and `8·π·μ·a³`, exact
symmetry of the assembled mobility, positive definiteness on the whole space including the six rigid
modes, and strictly negative hydrodynamic power. **None of them is exercised by this entry**, because
this entry registers no operator. Saying which oracles exist and that this entry does not use them is
the honest form of this field; claiming them would promote a declaration to a measurement.

What *is* independently checkable here is the declaration's internal consistency, and it is checked:
the scope tag against the owned-state rule, the owned state against the velocity derivation in §4, and
the claim list against the two named exclusions (no wall, no inertia).

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_environment_surface.py::test_the_medium_owns_a_committed_surface_reference` | `extracellular_medium` owns `committed_surface_reference_um`, and its role declares it a state owner rather than a boundary condition. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_medium_regularises_the_six_rigid_body_modes_physically` | The role names the rigid-body modes and carries the closed-form sphere limit, so the entry states a law rather than a mechanism. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_medium_is_the_sole_owner_of_surface_dissipation` | Sole-owner declaration present, with the double-count consequence named. |
| Positive | `tests/state/test_census_environment_surface.py::test_every_entry_validates` | The whole manifest passes `validate_registry`, including this entry. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_an_explicit_entry_that_owns_nothing_is_refused` | An `E` entry with an empty `owned_state` is rejected. This is the exact shape the source's declaration has — `owns_geometry=False, dynamically_evolving=False` and no state list — so the control fails if Aleph ever lets that shape be registered as an explicit owner. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_a_boundary_entry_that_owns_state_is_refused` | The inverse: a `B` entry that lists owned state is rejected, so the medium cannot be demoted to a boundary while keeping its committed reference. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_a_contract_type_missing_a_census_field_is_refused` | A contract type without the census fields is rejected rather than silently skipping the invariants that use them — the failure mode a reshaped upstream contract would otherwise produce. |

## 10. Numerical and precision envelope

No arithmetic is performed by this entry, so it introduces no tolerance and no precision requirement.
That is a real answer rather than an evasion, and it has a consequence worth stating: **nothing here
can be validated numerically, so nothing here may be quoted as a numerical result.** The entry sits at
the evidence class of a declaration.

The precision envelope it *records for a future implementation*: float64 throughout, exact symmetry of
the mobility as an IEEE statement rather than a tolerance (the kernel is even in the separation
vector, so the two blocks are the same products in the same order), the regularisation length tied to
the mesh's own spacing so that refining the surface is the lever and lowering the regulariser is not,
and the iterative solve's residual reduction reported rather than asserted. Those are requirements on
a later entry, not claims of this one.

## 11. Production-backend residency and transfer

This entry is **pure host-side declarative data** — module-level frozen dataclasses in
`aleph/state/census_environment_surface.py`, with no array, no device allocation and no kernel. It
never runs on the production backend, and it must not import one: the census is what a backend is
configured *from*, so a dependency in that direction would make the declaration unreadable without the
runtime.

`aleph/state/**` may not import `validation/**`, and this module imports only the standard library and
`aleph`, asserted by
`tests/state/test_census_environment_surface.py::test_the_module_imports_nothing_but_the_standard_library_and_aleph`.

A future exterior-Stokes implementation would be device-resident (surface quadrature, traction,
committed reference), with the viscosity and the time step as the only per-step host scalars. That
belongs to that entry's §11, not this one's.

## 12. Comments and docstrings to discard

No source prose survives, because none was carried: every field of the registered entry was written
here in Aleph's vocabulary from the derivation in §4.

Specifically **not** carried across, and each replaced rather than merely deleted:

- provider package and module paths, decision identifiers, gate labels, checklist item numbers,
  severity labels and internal document citations — replaced by the physical statement of *what*, with
  no reference to *where it was decided*;
- the source's `Sanity Gate` docstring block, which is a convention of that tree — replaced by the
  invariant list in §5, each item naming a test that runs;
- the source's own component-role vocabulary, which has no scope tag and therefore cannot express six
  of this lane's ten entries — replaced by the `E/I/H/B/X` discipline;
- any working-tree assumption. This entry states the digest it read.

The module is readable with the reference archive absent from the machine, which is the half of the
porting policy that is easiest to skip.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`. The code named in §1 has landed and its controls pass (83 tests in `tests/state/test_census_environment_surface.py`), but the entry stays `PROPOSED` because the registered *content* is a modelling declaration that no human has reviewed, The contract it registers into has since landed and this data now constructs against it, so moving to `ACCEPTED` now needs only the declaration to be read by the PI. |
| Reviewer | agent-proposed, **unratified**. |
| Rollback | Delete the `extracellular_medium` entry from `EXTERNAL_ENVIRONMENT`. Breaks `tests/state/test_census_environment_surface.py::test_the_module_registers_exactly_the_expected_entries` and the four controls in §8. Nothing computational depends on it, because it computes nothing. |

## 14. Honest limits

- **Unratified**, as above.
- The contract type **did not exist on disk when this entry was written and landed while this lane
  ran.** The module imports it defensively and falls back to a local dataclass of identical shape;
  which one is live is recorded in `CONTRACT_SOURCE` and asserted by
  `test_the_live_contract_source_is_recorded`. It now reads `aleph.state.schema`, so the real contract
  constructs and validates this data. Field *names* match exactly. If a field's intended semantics
  differ from what was assumed here, this entry's data could be wrong in a way no test here can see.
- The reference was **read, never run.** No GPU authorization exists, so nothing in §6 is a numerical
  claim about the source. The `dynamically_evolving=False` contradiction is a reading of two files at
  one commit; it has not been demonstrated by executing a rejected step.
- The reference's own artifacts for this component self-label as blocked from quoting magnitudes, and
  Aleph inherits nothing numerical from them. `citation_status` is `UNSOURCED` and the medium
  viscosity has no value anywhere in Aleph.
- The owned-state list is a *declaration of what would be owned*, derived from the velocity argument in
  §4. No Aleph implementation owns it yet, so the list is unexercised: an implementation could discover
  it needs more.
- What "sole owner of surface dissipation" means operationally — which components a guard must inspect
  and at what point in assembly — is not specified here and needs its own entry.
