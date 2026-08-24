# ALEPH-PORT-1603 — The `membrane` registry entry, and the reservoir that owns no state

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1603` |
| Lane | `L16 registry data — external environment and cell surface` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Nothing was ported. One energy/force inconsistency found in the reference is recorded in §6, and one live contradiction about who owns the interior pressure. |

---

> **Ownership amendment, 2026-08-03.** `ALEPH-PORT-3615` resolves the contradiction named in this
> entry: the closed membrane owns the osmotic volume; the cortex is an explicit filament graph. The
> reservoir and trafficking portions of this entry are unchanged.

## 1. Aleph API

The registry entries named `membrane`, `membrane_reservoir` and `membrane_trafficking`. They are filed
together because the reservoir is registered *under* the membrane and trafficking's only channel into
the model is the reservoir; separating them would put each one's meaning in a different file.

```python
from aleph.state.census_environment_surface import (
    CELL_SURFACE,
    COMPOSITE_TAGS_AS_WRITTEN,
    contract_for,
    scope_of,
)
```

Aleph target file: `aleph/state/census_environment_surface.py`.

Declared data only. **No bending energy, no areal force law, no tension kernel and no contact law** is
authorised here. The membrane's discrete Helfrich mechanics is a separate port with its own controls
(`ALEPH-PORT-1101`); this entry must not be read as evidence for it.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/membrane_area.py`, `ffn_sim/ff/membrane_surface.py`, `ffn_sim/ac/cell/compartments.py`, `ffn_sim/common/compartments.py`, `ffn_sim/ac/cell/membrane_pressure.py`, `ffn_sim/ac/engine/contracts.py` |
| Source symbol(s) | `MembraneAreaCard`, `effective_rest_area`, `area_tension`, `area_stiffness`; `MembraneMesh`, `membrane_area_kernel`, `membrane_area_reduce_kernel`, `reservoir_tension`, `helfrich_bending_kernel`, `erm_tether_kernel`; `MembraneCompartment`, `MEMBRANE_DEFAULTS`; `rsf_apply` (radial-shell law 1); and the `ComponentContract` declaration of `membrane` |
| Read from | **`git show be0e5876:<path>`** for all six paths. The working tree was not read for content. |
| Working tree == commit? | **yes** for all six, verified by digest: `ac/engine/membrane_area.py` = `sha256:a4b0183f7dd1a2c59c1af3ffab84781e5f8770760544b77c639a696acd43396c`; `ff/membrane_surface.py` = `sha256:66210eaeeb4b15d627b8ae739ac18b564e4ad8d21926c995de0fb7d00202457b`; `ac/cell/compartments.py` = `sha256:566905d3c0c8f23218b11739a3c8a0f26e37d09578ad37e6fc8b83aa6532500b`; `common/compartments.py` = `sha256:938391b418ca928ca53e465dad00b5f29f8f2a493f43643a6804a9de4383222b`; `ac/engine/contracts.py` = `sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72`. |
| Working-tree caveat | The tree carries **31 uncommitted changes** and one file this lane audited does differ from the commit (`ALEPH-PORT-1604` §2). Per-file digests were computed, not assumed. |

## 3. Why source-derived porting beats clean-room

**It does not.** RE-DERIVE. This is the entry where the question is closest to interesting, so the
reasoning is given in full rather than asserted.

The clamp `A0_eff = clamp(A, A0, A0·(1 + capacity))` is a real piece of modelling, and its three
properties — an exactly flat plateau inside the reserve, a stiff upturn outside it, and an asymmetric
lower clamp so a sub-rest sheet is slack rather than pre-compressed — are the kind of case analysis
`ports/TEMPLATE.md` §3 says is easy to get subtly wrong. That is a genuine argument for porting.

It fails for a reason specific to Aleph's situation: **the clamp is already in Aleph's own source
material.** `docs/manuscripts/extracted/APPENDIX_A_mechanics_registry.txt` states the formula, the
plateau behaviour and the exhaustion condition directly. So the choice is not "port or re-derive from
nothing"; it is "port code, or transcribe a formula Aleph already holds in a document it owns". The
second is strictly cheaper and carries no provider prose.

Two further reasons:

1. **The unit of approval here owns no state and applies no force.** The reservoir is registered at `H`.
   There is nothing to port but a one-line expression.
2. **The asymmetry of the clamp is derivable in one sentence, and the sentence is the load-bearing
   part.** A reservoir of folds can pay area out and cannot absorb it, so below `A0` the sheet is slack;
   a law that returned a negative tension there would have a bilayer carrying compression, which a
   bilayer cannot do. Once that is said, the `max(A0, min(A, ceiling))` form is forced. An author who
   has the sentence cannot get the clamp wrong, and an author who does not have it would not be saved by
   copying the expression.

So: the formula is transcribed from Aleph's own extracted appendix, the derivation of its shape is
written here, and no code crosses. The reference's *handling* of `capacity` — required from the caller,
no default, provenance mandatory — is independently the right policy and Aleph reaches the same place
from PD-002's rule that anything outside declared scope is a typed refusal rather than a default.

## 4. Physical or mathematical law represented

Three statements, all re-derived here.

**The membrane owns the osmotic load.** The osmotic pressure difference is set across a closed
compartment boundary. The explicit cortical filament graph has no enclosed-volume functional, while
the oriented membrane does, so the osmotic energy is a function of membrane volume alone. Pressure
therefore loads the membrane directly. Contact and ERM tethers exchange mechanics with the cortical
graph but are not a pressure relay. A second field on either body would duplicate or misroute the load
and is refused structurally. `ALEPH-PORT-3615` owns the integration; this entry records the registry
contract.

**The reservoir sets a rest area, not a force.** With `A0_eff = clamp(A, A0, A0·(1+capacity))` the areal
strain is `(A − A0_eff)/A0_eff`. Inside the reserve, `A0_eff = A` and the strain is **identically zero**
— not small, zero — so the tension is the plateau value bit for bit and `dσ/dA = 0` exactly. Once
`A > A0·(1+capacity)`, `A0_eff` sticks at the ceiling and `dσ/dA = K_A/A0_eff > 0`. The law is
continuous at the ceiling, where the tension equals the plateau from both sides; a jump there would be
a spurious force at a configuration nothing physical distinguishes. Both limits remain reachable and
both are meaningful: unbounded capacity is the pure constant-tension plateau, zero capacity is pure
areal elasticity measured from `A0`. So this is a strict generalisation of each rather than a third
model — which is what licenses using it without invalidating a result obtained in either limit.

**Trafficking cannot be recovered from the reservoir.** The reservoir is one dimensionless number with
no kinetics. A trafficking claim is a claim about *events*: rates, load dependence, delivery to a
location. A number with no time dependence cannot be decomposed into a fold contribution and a delivery
contribution, and cannot change during a run in response to load. So `capacity` recovered by inference
is a lumped signature and is **not** a trafficking measurement. This is a statement about what an
inference output may be called, and it is the substitution the `X/H` entry exists to forbid.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `vertex_positions_um` | µm | m | finite |
| `total_area_um2` | µm² | m² | > 0 |
| `membrane_tension_pn_per_um` | pN/µm | N/m | ≥ 0; a bilayer carries no compression |
| `discrete_mean_curvature_per_um` | 1/µm | 1/m | finite |
| `capacity` (a parameter of the reservoir, owned by nobody) | dimensionless | dimensionless | ≥ 0, unbounded above |
| areal strain | dimensionless | dimensionless | ≥ 0 by the clamp |

Singular and boundary cases:

- **`A ≤ 0` or `A0 ≤ 0`** is not a slack membrane, it is a broken surface. Refused, not clamped.
- **`A < A0`**: slack, tension floors at the plateau. Returning `γ + K_A·(A−A0)/A0 < 0` would be a
  bilayer under compression.
- **`capacity = 0`**: legal and meaningful (pure areal elasticity), not a degenerate case.
- **Unbounded capacity**: legal and meaningful (pure plateau).
- **`A` exactly at the ceiling**: the tension equals the plateau from both sides. Tested from both.

Invariants, each with the control that asserts it:

- **I1.** The membrane is at `E` and owns non-empty state including geometry, area, curvature, tension
  and a contact quadrature. `tests/state/test_census_environment_surface.py::test_every_explicit_entry_owns_state`.
- **I2.** The membrane declares that it owns the osmotic load as the sole closed surface, while the
  cortical filament graph receives mechanics through connectors.
  `tests/state/test_census_environment_surface.py::test_the_membrane_entry_declares_that_it_owns_the_turgor`.
- **I3.** Inside the reserve the areal strain is exactly zero, asserted with `== 0.0` rather than a
  tolerance, because the branch is identically zero and not merely small.
  `tests/state/test_census_environment_surface.py::test_inside_the_reservoir_the_rest_area_follows_the_area_so_strain_is_exactly_zero`.
- **I4.** Past the reserve the rest area sticks and the strain becomes strictly positive.
  `tests/state/test_census_environment_surface.py::test_past_the_reservoir_the_rest_area_sticks_and_strain_becomes_positive`.
- **I5.** Below `A0` the sheet is slack, not pre-compressed.
  `tests/state/test_census_environment_surface.py::test_below_rest_area_the_sheet_is_slack_rather_than_pre_compressed`.
- **I6.** The law is continuous at the ceiling.
  `tests/state/test_census_environment_surface.py::test_the_law_is_continuous_at_the_ceiling`.
- **I7.** Both limiting capacities are reachable.
  `tests/state/test_census_environment_surface.py::test_both_limiting_capacities_are_reachable`.
- **I8.** The entry states that `capacity` is an inference target and that no entry point may supply
  a value for it. It is asserted from `approximation`, **not** from `citation_status`: the state
  contract's citation vocabulary is a controlled set (`UNSOURCED`, `PROPOSED`, `SOURCED`,
  `INHERITED_UNVERIFIED`, `NOT_APPLICABLE`) so that a typo cannot invent a status which looks stronger
  than it is, and "this is an inference target" is a rule about the parameter rather than a state of
  its citation. `citation_status` is `UNSOURCED`, which is literally true.
  `tests/state/test_census_environment_surface.py::test_capacity_is_declared_an_inference_target_and_not_a_default`.
- **I9.** The reservoir owns no state and is parented on a registered explicit owner.
  `tests/state/test_census_environment_surface.py::test_the_reservoir_owns_no_state_and_is_parented_on_the_membrane`.

## 6. Source evidence class and known retractions

All six files were read at `be0e5876`; **none was executed** — five require Warp on CUDA, and this
session ran zero GPU jobs (PLAN §0.1). `ac/engine/membrane_area.py` is pure Python and could in
principle have been run; it was not, because running one file of six would produce a claim about that
file that reads as a claim about the group. Evidence class throughout: `AUDIT_READ` at a stated commit.

Retractions looked for in: the six module docstrings, `PLAN.md` §1.1/§7, and `ports/audit/`. Found and
worth recording:

- `ff/membrane_surface.py` carries its own **unsourced-constant warning** against the reservoir strain
  value it defines: it states plainly that the value is not sourced, that no fold-reservoir inventory is
  registered, that the reservoir is fold and microvillus unfolding rather than caveolae, and that runs
  should be treated as provisional. That is a retraction in substance and it is honest. Aleph inherits
  **no number** from it — `capacity` has no value anywhere in Aleph and its `citation_status` is
  `UNSOURCED`, with the stronger statement — an inference target that no entry point may default —
  carried in `approximation` where a test asserts it. See §5 I8 for why it is not a citation status.
- `ac/engine/membrane_area.py` records that the membrane in the production path had `dσ/dA = 0`
  everywhere, i.e. the infinite-reservoir limit, and that its own replacement is a strict generalisation
  of that. Also honest, and the argument is reproduced independently in §4 rather than taken.

**Finding 1 — energy and force disagree on the sub-reference branch of the radial-shell membrane law.**
In `ffn_sim/common/compartments.py`, `rsf_apply` law 1 computes the tension as
`γ_tot = pa + pb·(S − pc)/pc` and the per-bead force from `ΔP = 2γ_tot/R`, but the energy it reports
clamps the tension term: `Uc = pa·(S − pc)`, then `if Uc < 0: Uc = 0`.

Differentiating the unclamped energy `U = pa·(S − pc) + ½·pb·(S − pc)²/pc` with `S = 4πR²` gives
`dU/dR = 8πR·(pa + pb(S − pc)/pc) = 8πR·γ_tot`, and the code's total radial force is
`2γ_tot·S/R = 8πR·γ_tot`. So `F = −dU/dR` **exactly**, as long as the clamp is inactive. When `S < pc`
the clamp fires: the reported energy loses the `pa` term while the force keeps it. On that branch the
reported `U` is not the potential of the reported `F`, and any work-closure or energy-monotonicity check
evaluated there is comparing two different models. A descent that judges progress by this energy while
stepping along this force would see the two disagree only when the shell is below its reference area —
which is exactly the slack regime a reservoir model spends time in.

Established **by reading the kernel at the commit and differentiating it by hand**, not by running it;
no GPU authorization exists. Recorded as a finding, not as a fix: Aleph does not modify the reference.

**Finding 2 — two live descriptions of where the interior pressure enters, and they disagree.**
`ffn_sim/ac/cell/membrane_pressure.py` applies a pore-pressure boundary traction **directly to the live
plasma-membrane faces**, `f_face = (p_inside − p_ext)·A·n_out`, reconstructed from the resolved cytosol
pressure field. `ffn_sim/ac/cell/compartments.py` describes its `MembraneCompartment` as "also the
osmotic envelope". Meanwhile the declared architecture routes cortex-owned osmotic pressure to the
membrane through a compressive membrane-cortex contact, and `ffn_sim/ac/engine/cortex_state.py` lists
the pressure channel among the channels the cortex deliberately does **not** bind, attributing it to the
cytosol.

Two honest qualifications, because this could easily be over-claimed:

1. The membrane-face traction is a **coherent and arguably better** model on its own terms. If the
   cortex is permeable to water — which it is — then the fluid pressure acts on the bilayer, and reading
   it from a resolved field is more informative than a lumped scalar. `membrane_pressure.py` says it
   deliberately replaces a retired lumped turgor. This is not simply a bug.
2. "Osmotic envelope" in `compartments.py` may mean the surface that classifies the fluid domain rather
   than the owner of the pressure load. The phrase is ambiguous and is reported as ambiguous.

What is *not* ambiguous is that a reader of that tree cannot determine which surface owns the interior
load without opening three files, and that PLAN §1.1 separately reports the compressive membrane-cortex
contact as unimplemented while passing a dispatch coverage gate. So the connector that the declared
ownership depends on is declared and not evaluated, while a different path applies the load elsewhere,
and nothing compares the two. That is the same shape as PLAN §1.2's finding that two whole-cell
assemblies are alive at once — located here at named symbols and a commit rather than carried as audit
prose.

Aleph's response is structural: `OSMOTIC_ENVELOPE` is a single named constant, the membrane's
registered role states that it owns the load, the cortex role states that a filament graph cannot,
and `assert_single_osmotic_envelope` refuses a world that applies the load anywhere else or more than
once.

## 7. Independent oracle or derivation

The reservoir clamp is checked against an **independent re-implementation written from the derivation in
§4**, living in the test file as `_a0_eff`, and exercised over all three regimes plus both limits and
the continuity point. That satisfies `ports/TEMPLATE.md` §7's last option — an independent
implementation written from the derivation without the source open — and it is the reason this entry has
a real control rather than a string comparison.

Three properties of the law are asserted exactly rather than to a tolerance, which is stronger than a
numerical agreement: the strain is `== 0.0` inside the reserve, the clamped rest area is `== A0` below
`A0`, and `capacity → 0` gives `A0_eff == A0` at any area above it. Each of those is an identity of the
clamp, not a limit of it.

For the membrane's mechanics the oracles are the Helfrich sphere energy `8πκ`, the exact scale
invariance of the discrete bending energy, and the dilational-virial Laplace identity. **All belong to
`ALEPH-PORT-1101`; none is exercised here**, and this entry claims none of them.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_environment_surface.py::test_the_registered_formula_is_stated_verbatim` | The clamp appears in the registered representation exactly as `A0_eff(t) = clamp(A(t), A0, A0 * (1 + capacity))`, and the entry states that folds are not resolved. |
| Positive | `tests/state/test_census_environment_surface.py::test_inside_the_reservoir_the_rest_area_follows_the_area_so_strain_is_exactly_zero` | Areal strain is exactly `0.0` at four areas inside the reserve — the plateau is an identity, not an approximation. |
| Positive | `tests/state/test_census_environment_surface.py::test_past_the_reservoir_the_rest_area_sticks_and_strain_becomes_positive` | Above the ceiling the rest area equals `A0·(1+capacity)` and the strain is strictly positive at three areas. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_law_is_continuous_at_the_ceiling` | Strain vanishes from both sides at the ceiling, so no spurious force appears there. |
| Positive | `tests/state/test_census_environment_surface.py::test_both_limiting_capacities_are_reachable` | Zero capacity gives pure elasticity from `A0`; a large capacity gives the pure plateau. |
| Positive | `tests/state/test_census_environment_surface.py::test_the_membrane_entry_declares_that_it_owns_the_turgor` | The membrane's role states osmotic ownership and the closed-volume reason. |
| Positive | `tests/state/test_census_environment_surface.py::test_each_named_trafficking_mechanism_is_individually_excluded` | Seven named trafficking mechanisms are each declared absent, not merely covered by a blanket sentence. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_pressure_on_both_membrane_and_cortex_is_refused` | Applying the osmotic load to both is rejected with the double-count named. This is Finding 2's failure mode, planted deliberately. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_pressure_on_the_cortex_alone_is_refused_with_the_reason` | A filament graph with no closed triangulation is rejected as a pressure owner. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_a_homogenized_entry_that_owns_state_is_refused` | The reservoir cannot acquire owned state while staying at `H` — the promotion that would make it look like a resolved fold population without resolving one. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_below_rest_area_the_sheet_is_slack_rather_than_pre_compressed` | The lower clamp holds at two sub-rest areas; the test fails if the clamp is removed and a compressive branch appears. |
| Negative (must fail) | `tests/state/test_census_environment_surface.py::test_trafficking_only_channel_is_the_lumped_reservoir_and_that_is_not_a_measurement` | The entry states that a recovered capacity is not a trafficking measurement; the test fails if that prohibition is dropped. |

## 10. Numerical and precision envelope

float64, and the entry's arithmetic is one `max` of a `min`, so there is no conditioning question: the
clamp is monotone in each argument and introduces no cancellation.

Three assertions are made **exactly** (`==`), and the choice is deliberate rather than lucky. Inside the
reserve `A0_eff` is the identity function on `A`, so `A − A0_eff` is `0.0` bit for bit, not a small
residue; asserting a tolerance there would be weaker than the truth and would hide a clamp that had
started to leak. Below `A0` and at `capacity = 0` the result is the literal `A0` that was passed in, so
again exactness is available.

The continuity check is the one place a tolerance is unavoidable, because it perturbs the ceiling by
`±1e-12` relative and the products `A0·(1+capacity)` are not exactly representable. It is asserted at
`1e-15` on the below-ceiling side, where the strain is an exact zero, and `1e-11` on the above-ceiling
side, where the strain is the perturbation itself and so is `~1e-12` by construction — the tolerance is
one decade above the perturbation, not chosen against an observed value. Outside this envelope the
representation refuses (`A ≤ 0`, `A0 ≤ 0`, negative capacity) rather than degrading silently.

No claim here has a tolerance that was tuned to make it pass; every one is either an identity or one
decade above a perturbation the test itself chose.

## 11. Production-backend residency and transfer

Host-side declarative data. Frozen dataclasses in `aleph/state/census_environment_surface.py`; no array,
no device allocation, no kernel, no per-step transfer. It never runs on the production backend and must
not import one.

`_a0_eff` in the test file is host-side float64 and is a **test oracle**, not production code — which is
also a residency answer, and it is the one that keeps oracles out of the runtime. A future device-side
implementation of the tension law would reduce the total area on device and compute the clamp from a
device scalar; the only host round trip would be the reported tension. That belongs to that entry.

Imports are standard library and `aleph` only, asserted by
`test_the_module_imports_nothing_but_the_standard_library_and_aleph`. `aleph/state/**` may not import
`validation/**`.

## 12. Comments and docstrings to discard

No source prose survives. Discarded and replaced:

- the reference strain constant and every numeric membrane default. **No magnitude crosses**; the
  registry entry contains no number at all, asserted by
  `test_owned_state_is_named_rather_than_valued`;
- provider module paths, decision identifiers, hard-truth numbering, plan and framework references,
  and knowledge-base claim identifiers — replaced by the derivation in §4;
- the source's `Sanity Gate` docstring convention — replaced by the invariant list in §5, each item
  naming a test that runs;
- the source's argument for why lumping is defensible here. Aleph re-derives its own version and keeps
  the scope limit that goes with it: the reservoir is registered at `H` **with** an explicit
  unsupported-claims list and a re-entry condition, so the licence and its price are in the same object;
- the source's radial-shell membrane law entirely. It is a mean-radius sphere model, not a surface, and
  nothing about it is carried; it appears in this entry only as Finding 1.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`. The data has landed and 83 controls in `tests/state/test_census_environment_surface.py` pass, including the nine invariants of §5. It stays `PROPOSED` because the registered content is a modelling declaration no human has reviewed, The contract it registers into has since landed and this data now constructs against it, so the remaining blocker is review rather than a missing dependency. |
| Reviewer | agent-proposed, **unratified**. Finding 1's differentiation has been checked by hand only, and Finding 2 contains an ambiguity that is reported as such rather than resolved. |
| Rollback | Delete the `membrane`, `membrane_reservoir` and `membrane_trafficking` entries from `CELL_SURFACE`. Breaks the controls in §8 and leaves the `cortex` entry's non-ownership counterpart unstated, since the turgor rule needs both halves to be checkable. |

## 14. Honest limits

- **Unratified**, as above.
- The state-contract type **was not on disk when this entry was written and landed while this lane
  ran.** The module imports it defensively and records which definition is live in `CONTRACT_SOURCE`;
  that value is now `aleph.state.schema`, so the data is constructed and validated by the real
  contract and not by the fallback. What this does *not* establish is semantic agreement: the field
  *names* match exactly, and if a field's intended meaning differs from what was assumed here, no
  test in this lane can see it.
- The reference was **read, never run.** Finding 1 is a hand differentiation of a kernel that was not
  executed; it is a strong reading, not a measurement, and it would be settled in minutes by running the
  kernel on a sub-reference configuration — which requires a GPU authorization that does not exist.
- Finding 2 contains an **unresolved ambiguity**: whether "osmotic envelope" in that docstring names the
  pressure owner or the fluid-domain classifier. The finding is reported with the ambiguity intact rather
  than resolved in the direction that makes it sound worse.
- `capacity` has **no value** and no sourced range. Whether it is identifiable from the observation
  Aleph has chosen is an open question this entry does not address, and the entry deliberately makes it
  impossible to answer by supplying a default.
- The membrane's owned-state list is a declaration. `ALEPH-PORT-1101`'s implementation owns a subset of
  it; `contact_quadrature` and `accepted_surface_state` in particular are named here and not yet
  represented by anything on disk.
- **No trafficking claim is supported and none is planned.** The `X/H` entry lists ten prohibitions and
  that list is not exhaustive; it is the set this lane could name, and a claim not on it is not thereby
  permitted.
- Nothing here says how the reservoir's plateau interacts with the membrane's *bending* response. The
  two are separate terms and the entry does not assert they are independent.
