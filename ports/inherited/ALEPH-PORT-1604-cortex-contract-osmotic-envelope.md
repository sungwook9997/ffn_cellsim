# ALEPH-PORT-1604 — The `cortex` registry entry: the osmotic envelope, and a filament population that is its own

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1604` |
| Lane | `L16 registry data — external environment and cell surface` |
| Status | `REJECTED` |
| Superseded by | `ALEPH-PORT-3615`; filament mechanics remain under `ALEPH-PORT-2301` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` — with one **argument** inherited, which makes the inherited part a `DEFECT_STUDY`. No code crosses. |
| Verdict | **RE-DERIVE.** The only thing that crosses the boundary is the description of a mistake, and the entry is shaped so the mistake cannot be made silently. |

---

> **Supersession, 2026-08-03.** The PI-selected explicit cortex is a filament graph and has no
> enclosed volume. The membrane is now the sole osmotic envelope under `ALEPH-PORT-3615`. This entry
> remains as the auditable record of the rejected shell-envelope arrangement; its disjoint filament
> ownership argument survives, but its cortex-volume ownership does not.

## 1. Aleph API

The registry entry named `cortex`, plus the registry-level turgor declaration that entry depends on.

```python
from aleph.state.census_environment_surface import (
    CELL_SURFACE,
    FILAMENT_POPULATIONS_DISJOINT_FROM_CORTEX,
    FORBIDDEN_OSMOTIC_OWNERS,
    OSMOTIC_ENVELOPE,
    OsmoticEnvelopeError,
    assert_single_osmotic_envelope,
    contract_for,
)
```

Aleph target file: `aleph/state/census_environment_surface.py`.

Declared data plus **one guard**. `assert_single_osmotic_envelope` is the only executable thing this
lane contributes, and it exists because the rule it enforces cannot be checked by looking at a result.
No cortical mechanics — no filament bending, no crosslink law, no steric interaction, no motor kinetics
— is authorised here.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/cortex_state.py`, `ffn_sim/ac/engine/cortex_population.py`, `ffn_sim/ac/engine/contracts.py`, `ffn_sim/common/compartments.py` |
| Source symbol(s) | `CORTEX_CHANNELS_NOT_BOUND`, `CortexCompositeMechanics`, `CortexStateLedger`, `CortexStateTransaction`, `assert_component_state_disjoint`, `build_cortex_state_owner`, `cortex_segment_material_coordinates`, `_restore_vec3_if_rejected_kernel`; `CORTEX_COMPONENT`, `CortexPopulation`, `CortexTopology`, `CORTEX_CONNECTORS`; `rsf_apply` (radial-shell law 2, turgor); and the `ComponentContract` declaration of `cortex` |
| Read from | **`git show be0e5876:<path>`** for all four paths. For `cortex_state.py` the working tree was additionally diffed against the commit — see the next two rows, which is the whole reason those rows exist. |
| Working tree == commit? | **NO for `ffn_sim/ac/engine/cortex_state.py`.** Commit digest `sha256:d28f254496f7398169ec0f14cd78970a60f7b0a830772a1ef39c2dbad938eab3`; the tree differs by **+34 lines, 0 deletions**, uncommitted. **Yes** for the other three: `cortex_population.py` = `sha256:04d4542cf0c07a1f35bf1600eda4b40828e60e0ecebb60d71388ac363e903b70`, `contracts.py` = `sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72`, `common/compartments.py` = `sha256:938391b418ca928ca53e465dad00b5f29f8f2a493f43643a6804a9de4383222b`. |
| What the uncommitted delta contains | A double-count guard that **exists only in the working tree**: a map from this owner's force-channel names to the names the incumbent driver knows the same kernels by, and a function that *derives* the set of incumbent channels a caller must suppress from the channels actually bound — refusing on any bound channel with no mapping. Its own comment states the motive: the branch-angle channel appears only for a mixed filament population, so a hand-written suppression set is right for the default cortex and silently doubles that force for the mixed one. An audit that read only `be0e5876` would conclude the guard does not exist. It does; it is just not committed. This is exactly the failure mode `ports/TEMPLATE.md` §2 added those two rows for, and this entry is the first in the ledger to actually trip it. |

## 3. Why source-derived porting beats clean-room

**Not for code.** RE-DERIVE. But for **one argument** it does, and this is the case
`ports/TEMPLATE.md` §3 describes as legitimate: an empirical failure mode discovered by somebody else's
run, which a clean-room author would plausibly not have thought of.

The argument, restated in Aleph's words:

> An endpoint that is a *view* onto another owner's arrays cannot close an adjoint pair. If a
> component's "force array" is the same allocation that other components accumulate into, then the two
> halves of the Newton pair are the same numbers, and a force-closure check on that pair is a tautology
> that passes whatever the physics does.

That is worth having. A clean-room author designing a component registry would very likely let a
component be a view for convenience — it saves an allocation and it looks harmless — and would then
build a force-closure gate that could never fail. Aleph's runtime lane already independently arrived at
the general principle that a coverage gate must witness evaluation and not merely registration
(`ALEPH-PORT-303`), but the *specific* form here, that shared storage makes an adjoint check vacuous, is
the sharper statement and it came from the reference.

So the entry registers `filament_segment_positions_um` and `filament_force_pn` as **privately owned**,
and says why in the entry's own prose. Nothing else crosses: no kernel, no dataclass, no channel name,
no population census, no number.

The second argument that crosses is the turgor-ownership defect, and it is already owned by
`ALEPH-PORT-1102`. This entry records the cortex's half of it and does not restate the reasoning.

## 4. Physical or mathematical law represented

**The cortex is the osmotic envelope.** The osmotic pressure difference is set by the solute imbalance
across the envelope that separates the two compartments. The cortex is the innermost closed surface
bounding the cytoplasm; it is what the interior pushes on; so the osmotic energy is a function of the
cortex's enclosed volume and of nothing else:

```
E_turgor = E_turgor(V_cortex)          and nothing else
```

The membrane is outside it and receives the load **only** through the membrane-cortex connectors — a
compressive non-adhesive contact for the outward push, and a tensile tether that cannot carry
compression, so a slack or broken tether transmits exactly zero rather than a little.

**Why the duplication is dangerous rather than merely wrong.** Applying `ΔP` to both surfaces adds the
same gradient term twice. The total is exactly `2ΔP∇V`. The energy is still a scalar potential, the
force is still its exact gradient, the descent still decreases it monotonically, and the relaxation
still converges to a smooth sphere. Only the radius is wrong — and by an amount indistinguishable from
a different choice of tension. There is no residual to inspect, no non-convergence to notice and no
symmetry that breaks. The failure has no signature at all in the answer, which is why it must be
refused *structurally, at assembly, before any number exists*. That is the entire justification for
`assert_single_osmotic_envelope` being a callable guard rather than a sentence in a docstring.

**Co-location is not ownership.** Every physical filament has exactly one component owner. The cortex's
population is disjoint from those owned by `sf_arc`, `lamellipodium` and `filopodium`. Sharing a region
of space, an array, or an index range creates no mechanical connection — a force exists between two
populations only where a connector declares it, and each such connector joins the populations without
merging their nodes or their filament identities. The consequence is directional and worth stating: a
filament cannot be moved between owners by relabelling it, and an explicit new actin system cannot be
created by reclassifying cortical filaments.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `filament_segment_positions_um` | µm | m | finite |
| `filament_force_pn` | pN | N | finite |
| `filament_material_coordinates` | µm (arc length) | m | ≥ 0, monotone along a filament |
| `filament_polarity` | dimensionless (sign/direction) | — | unit-norm direction |
| `topology_epoch` | dimensionless (counter) | — | non-decreasing integer |
| `myosin_binding_site_occupancy` | dimensionless | — | integer occupancy per site |
| turgor `ΔP` (a context parameter, owned by nobody) | pN/µm² | Pa | any sign; positive inflates |

`1 pN/µm² = 1 Pa` exactly.

Singular and boundary cases:

- **Zero osmotic load owners is refused, not treated as zero pressure.** An unowned turgor is not a
  depressurised cell; it is a world in which the envelope carries no interior push at all, which is a
  different model and must be declared rather than arrived at by omission.
- **Two loads on the cortex are still two loads.** The guard counts applications, not distinct owners,
  because the arithmetic error is the same either way.
- **A load on any owner other than the cortex is refused**, with the reason named for the membrane
  specifically, since the membrane is the load's most likely second home.
- **Shared storage.** A component whose force array is another's allocation is a vacuous adjoint
  endpoint. Registered as private ownership; the check that private really is private belongs to a
  runtime entry, not to this one.

Invariants, each with the control that asserts it:

- **I1.** The declared osmotic envelope is the cortex, and it is registered at `E`.
  `tests/state/test_census_environment_surface.py::test_the_declared_osmotic_envelope_is_the_cortex`.
- **I2.** The cortex's role declares that it owns the turgor load and is the osmotic envelope.
  `tests/state/test_census_environment_surface.py::test_the_cortex_entry_declares_that_it_owns_the_turgor`.
- **I3.** The role forbids applying pressure independently to both, and names the consequence — a
  plausible-looking sphere at the wrong radius, not a crash.
  `tests/state/test_census_environment_surface.py::test_the_cortex_entry_forbids_applying_pressure_to_both`.
- **I4.** The membrane's role states the counterpart: it does not own the osmotic load.
  `tests/state/test_census_environment_surface.py::test_the_membrane_entry_declares_that_it_does_not_own_the_osmotic_load`.
- **I5.** A single load on the cortex is accepted, so the guard is not merely a rejecter.
  `tests/state/test_census_environment_surface.py::test_a_single_load_on_the_cortex_is_accepted`.
- **I6.** The cortex declares a disjoint filament population and names all three other owners
  explicitly, and states that co-location is not ownership.
  `tests/state/test_census_environment_surface.py::test_the_cortex_owns_a_disjoint_filament_population`.
- **I7.** Owning the load is not a permeability claim: the entry states there is no solute transport
  anywhere in the registry and that `ΔP` is a declared parameter.
  `tests/state/test_census_environment_surface.py::test_the_cortex_does_not_claim_to_be_an_osmotic_barrier`.

## 6. Source evidence class and known retractions

Read at `be0e5876`; **not executed** — all four files require Warp on CUDA and this session ran zero GPU
jobs (PLAN §0.1). Evidence class `AUDIT_READ` at a stated commit. `cortex_state.py` was additionally
diffed against the working tree, which is a comparison of two revisions and not an execution.

Retractions and self-labelled status looked for in: the four module docstrings, `PLAN.md` §1.1 and §7,
and `ports/audit/`. Found:

- `cortex_state.py` records that the cortex was previously a *bind-target port* rather than a
  participant: its force array was the same allocation other components accumulated into, so no adjoint
  pair terminating on the cortex could be closed. That is a retraction of a prior arrangement and it is
  the argument §3 inherits.
- The same file records a set of force channels it deliberately does **not** bind, each with a reason,
  and states that a force absent for a stated reason is a different thing from a force forgotten. That
  is the right discipline; Aleph's equivalent is the `unsupported_claims` list, which is a per-entry
  field rather than a module constant.
- PLAN §1.1 reports, from a prior read-only audit, that the compressive membrane-cortex contact is
  unimplemented while passing a dispatch coverage gate. This entry does not re-verify that and does not
  depend on it.

**Finding 1 — the pressure channel is declared not-bound on the cortex.** `cortex_state.py` lists the
pressure coupling among the channels this owner does not bind, attributing it to the cytosol and to the
porous-transfer connector. Taken with `ALEPH-PORT-1603` Finding 2 — where a pore-pressure traction is
applied directly to the membrane faces — the interior load in that tree enters at the **membrane**,
while the declared architecture says the cortex owns it and relays it. Both descriptions are live and
nothing compares them.

This is worth stating carefully, because the reference's choice is defensible: if the cortex is
permeable to water, the fluid pressure genuinely acts on the bilayer, and a resolved pressure field is
more informative than a lumped scalar. The defect is not the physics of either option. It is that a
reader cannot determine which one the tree implements without opening three files, and that the
connector the declared ownership depends on is declared without being evaluated. Aleph's answer is that
the envelope is one named constant and the guard refuses everything else.

**Finding 2 — a lumped turgor is applied by mean radius on tagged particles, with an energy that does
not match its force on one branch.** `common/compartments.py` `rsf_apply` law 2 computes
`ΔP = pa − pb·(V − pc)/pc` from a **mean radius** of tagged particles and applies `+ΔP·A_i·n̂`
outward. Two observations:

1. It is a sphere model, not a surface: the load depends only on a mean radius, so it cannot represent
   any non-spherical envelope, and it is applied to whichever tag range the caller selects — which is a
   configuration choice, not an ownership declaration. Whether it lands on a cortex or on a membrane is
   decided by an argument at the call site.
2. Law 1 of the same kernel has an energy/force inconsistency on its sub-reference branch, established
   by hand differentiation and recorded in full as `ALEPH-PORT-1603` Finding 1.

**Finding 3 — the double-count guard for the cortex is uncommitted.** See §2. The map from this owner's
channels to the incumbent's, and the derived suppression set that refuses on an unmapped channel, exist
only in the working tree at the time of this audit. This matters twice over: an audit citing
`be0e5876` would report the guard as absent, and a *reader* of the tree gets protection that no commit
records. Reported as a provenance finding, not as a criticism of the guard, which is well designed —
deriving the suppression set from what is actually bound is the correct shape, and its own comment names
the exact silent-doubling case a hand-written set would produce.

## 7. Independent oracle or derivation

**The double-application check is an exact arithmetic identity, and that is the oracle.** With one
envelope the pressure force is `F_i = ΔP ∇_i V`. With two, it is `2ΔP ∇_i V`. The difference is exactly
the duplicated term, vertex by vertex, with no tolerance and no continuum limit — so a control can
assert the difference *equals* the duplicate rather than merely differs from zero. That identity is
owned by `ALEPH-PORT-1102` and is exercised there against a real assembly; this entry's contribution is
the structural refusal that makes the assembly unbuildable in the first place, which is checked by
rejection rather than by measurement.

Two further checks that are Aleph's own and require no reference:

- **`∮ n dA = 0` on any closed surface**, so a uniform pressure exerts no resultant force. This holds
  whatever the discretisation does, which makes it a check on the load path rather than on the mesh.
- **The dilational-virial Laplace relation.** From `A(λx) = λ²A`, `V(λx) = λ³V` and the exact scale
  invariance of the discrete bending energy, stationarity of `σA − ΔP V` under uniform dilation gives
  `ΔP = 2σA/(3V)`, which on a sphere is `ΔP = 2σ/R`. A doubled turgor shifts the recovered `σ` by
  exactly a factor of two — so the recovered tension is where the duplication *is* visible, and it is
  visible nowhere else. That is the strongest available statement about why this rule matters, and it
  belongs to `ALEPH-PORT-1102`'s controls, not to this entry's.

For cortical mechanics there is **no oracle here and none claimed.** The cortex's constitutive
representation is not implemented in Aleph.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_environment_surface.py::test_the_declared_osmotic_envelope_is_the_membrane` | The superseding control pins the envelope constant to the closed membrane. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_cortex_entry_declares_that_it_does_not_own_the_osmotic_load` | The superseding role records why the filament graph cannot own pressure. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_cortex_entry_forbids_applying_pressure_to_both` | The role forbids independent application to both and names the consequence — counted twice, wrong radius. |
| Positive | `tests/state/test_census_environment_surface.py::test_a_single_load_on_the_membrane_is_accepted` | The current guard accepts the one correct configuration, so it is not a function that rejects everything. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_cortex_owns_a_disjoint_filament_population` | Disjointness is declared, all three other owners are named, and co-location is denied as a connection. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_registry_does_not_claim_a_resolved_osmotic_barrier` | Pressure ownership remains separate from an unsupported solute-barrier claim. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_pressure_on_both_membrane_and_cortex_is_refused` | The exact defect: both surfaces carrying the load is rejected, with the double-count named. Fails if the guard stops refusing. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_two_loads_on_the_membrane_are_still_two_loads` | Two applications to the current owner are also rejected — the arithmetic error does not care that the owner is correct. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_no_owner_at_all_is_refused_rather_than_treated_as_zero_pressure` | An empty set is refused rather than silently meaning zero pressure, so "nothing was applied" and "no load exists" stay distinguishable. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_an_unrelated_owner_is_refused` | A load on an owner from another group is refused, so the guard is not a membrane special case. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_a_registry_whose_osmotic_envelope_is_not_registered_is_refused` | Removing the membrane envelope from the manifest must fail validation rather than leaving the load unowned. |

## 10. Numerical and precision envelope

The registry data performs no arithmetic. `assert_single_osmotic_envelope` performs **no floating-point
arithmetic at all**: it counts a sequence and compares strings. So it has no tolerance, no precision
requirement and no conditioning, and it cannot degrade — which is the argument for putting the check
here rather than on a computed residual.

That contrast is the substance of this field, and it comes from a measured result of this project. PLAN
§10 and §11 record that `max_residual_force_pn` **cannot** gate the Laplace error: the error equals a
signed dilational virial, the gate is a max norm, and on one fixed mesh twenty stopping points gave five
sign changes and a 603× spread in `|error|`. A doubled turgor is worse than that, because it does not
appear in a residual at any tolerance — it is absorbed into the equilibrium radius. A structural check
with no tolerance is therefore not a weaker instrument than a numerical one here; it is the only
instrument that sees this failure.

The one numerical requirement recorded for later: `topology_epoch` must be an exact integer comparison,
never a tolerance. Stale and current are discrete facts.

## 11. Production-backend residency and transfer

Host-side. The registry entries are frozen dataclasses in
`aleph/state/census_environment_surface.py`; no array, no device allocation, no kernel.

`assert_single_osmotic_envelope` is a **host-side assembly-time structural check** and deliberately
never touches a device. That placement is the point: it runs before the first step, so it cannot be
skipped by a run that fails early, and it cannot be made conditional on a device predicate. A guard
that only fires during stepping would be absent from exactly the runs that never get there.

A future cortical implementation would be device-resident (segment positions, private force array,
crosslink topology), with the turgor as a host scalar. That belongs to that entry.

Imports are standard library and `aleph` only, asserted by
`test_the_module_imports_nothing_but_the_standard_library_and_aleph`. `aleph/state/**` may not import
`validation/**`.

## 12. Comments and docstrings to discard

No source prose survives. Discarded and replaced:

- provider module paths, driver and facade names, dispatch-slot names, execution-plan and track
  references, decision identifiers and gate labels. The inherited argument in §3 is stated in terms of
  **what goes wrong physically** — two halves of a Newton pair being the same numbers — with no
  reference to the other project's identifiers, so it is readable with that archive absent;
- the source's channel-name vocabulary and its not-bound constant. Aleph's equivalent is the per-entry
  `unsupported_claims` list, which attaches the omission to the compartment rather than to a module;
- the source's population census and every filament and node count. **No magnitude crosses**; the entry
  contains no number, asserted by `test_owned_state_is_named_rather_than_valued`;
- the source's `Sanity Gate` docstring convention — replaced by the invariant list in §5, each item
  naming a test that runs;
- the radial-shell turgor law entirely. It is a mean-radius sphere model applied to a caller-selected
  tag range, and nothing about it is carried; it appears here only as Finding 2.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`. The data and the guard have landed and 83 controls in `tests/state/test_census_environment_surface.py` pass, including the seven invariants of §5 and five deliberate rejections. It stays `PROPOSED` because the registered content is a modelling declaration no human has reviewed, The contract it registers into has since landed and this data now constructs against it, so the remaining blocker is review rather than a missing dependency. |
| Reviewer | agent-proposed, **unratified**. The choice of the cortex as the envelope is a modelling decision, not a derivation, and it is the decision most worth the PI's attention in this entry. |
| Rollback | Delete the `cortex` entry, `OSMOTIC_ENVELOPE`, `FORBIDDEN_OSMOTIC_OWNERS` and `assert_single_osmotic_envelope`. Breaks the controls in §8 and §9, and removes the only structural refusal of a doubled turgor that exists at census level — leaving the rule enforced only inside the vertical, where a world assembled by any other path would not meet it. |

## 14. Honest limits

- **Unratified**, as above. In particular, *that* the cortex is the envelope is inherited from the seed
  registry and is a modelling decision. A cortex freely permeable to water would put the load on the
  membrane instead, which is what the reference's field traction does, and that alternative is
  defensible. Aleph registers one answer and makes it refusable; it does not establish that the answer
  is right.
- The state-contract type **was not on disk when this entry was written and landed while this lane
  ran.** The module imports it defensively and records which definition is live in `CONTRACT_SOURCE`;
  that value is now `aleph.state.schema`, so the data is constructed and validated by the real
  contract and not by the fallback. What this does *not* establish is semantic agreement: the field
  *names* match exactly, and if a field's intended meaning differs from what was assumed here, no
  test in this lane can see it.
- The reference was **read, never run.** No GPU authorization exists.
- Finding 3 is a **provenance** finding and depends on the working-tree state at the moment of this
  audit. That state is mutable and will change; the digests recorded in §2 are what make the finding
  reproducible at all.
- `assert_single_osmotic_envelope` takes the names a caller *says* it applied the load to. **It cannot
  verify that claim.** A caller that applies a load and does not report it defeats the guard entirely.
  Closing that gap needs the check to run over a pipeline's actual registrations, which is a runtime
  concern and is where `ALEPH-PORT-1102`'s `assert_single_pressure_path` sits; the two are complementary
  and neither subsumes the other. This is the sharpest limit in this entry.
- The cortex's owned-state list is a declaration. **No Aleph cortical implementation exists** — the
  vertical's cortex is a single areal spring on a closed shell (`ALEPH-PORT-1102`), which owns almost
  none of what is listed here and represents no filament at all.
- Disjointness of filament populations is **declared, not enforced.** Nothing in this lane can check
  that a filament has exactly one owner, because no filament exists yet. When populations exist, the
  check is a storage-disjointness assertion across owners, and it belongs to a runtime entry.
- The `unsupported_claims` list is the set this lane could name. A claim absent from it is not thereby
  permitted.
