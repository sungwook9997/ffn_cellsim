# ALEPH-PORT-2801 — ECM collagen fibres as a rod/cable element library

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2801` |
| Lane | `L28 ecm / extracellular_medium` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (**after** the code, which is a breach of PLAN §0.2.5 — see §14) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.vertical.ecm import (
    Crosslink,
    ECMOwnershipError,
    ECMTopologyError,
    Fibre,
    FibreCard,
    FibreElement,
    LigandSite,
    MaterialPoint,
    MaterialPointSink,
    WorldBoundaryStateError,
    axial_energy_and_forces,
    bending_energy_and_forces,
    crosslink_energy_and_forces,
    default_collagen_card,
    owned_state_keys,
    straight_fibre_nodes,
)
```

**This is an element library, not a state owner.** The registered `ecm` contract declares seven state
blocks and this module allocates none of them. There is no `CollagenNetwork`, no
`snapshot`/`commit`/`rollback`, no topology epoch, no ligand resolution against a live network, no
remodelling, no far-field anchor and no `world_boundary` object. See §14, which is the section of
this entry that matters most.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — verified to resolve, read-only |
| Source path | `ffn_sim/ac/engine/ecm_mechanics.py`, `ffn_sim/ac/engine/ecm_world.py` — **named, not read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened for this module. |
| Working tree == commit? | not applicable — nothing was read, so there is no revision to have mis-cited |

The two source paths are recorded so a later auditor knows where the reference's ECM lives and can
establish independently whether Aleph's re-derivation agrees or disagrees with it. Listing them is
**not** a claim that they were consulted. The only fact this entry asserts about the reference is
that the commit hash resolves.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** The module was written from the registered `ecm` contract in
`aleph/state/census_environment_surface.py` and from Appendix A's `[E] ecm` entry, both of which are
Aleph's own documents. A Hookean axial spring with a rest-length normalisation, a tangent-difference
bending stencil, and a four-node barycentric spring are laws any competent author re-derives in an
afternoon; per `ports/TEMPLATE.md` §3 that is a reason to write clean-room and not a reason to port.

## 4. Physical or mathematical law represented

Three energy terms, each with its exact gradient. Units µm / pN / pN·µm throughout.

**Axial.** `E = Σ_s ((1 − d_s) k_s / 2) (L_s − L0_s)² / L0_s`, giving segment tension
`T = (1 − d) k (L − L0) / L0`, which pulls the two ends together. Dividing by the rest length makes
the modulus resolution-independent, and that has a checkable consequence: a chain in equilibrium
under an end load `F` carries `T = F` in every segment and extends by `ΔL = F L0 / k` **regardless of
how many segments it was cut into**. A `tension_only` segment below rest length contributes
identically zero to energy and force — not a small number, an identical zero on the whole branch,
which is why the control asserts it with `==`.

**Bending.** On unit tangents `t1`, `t2` of the two segments meeting at interior node `i`:
`E_i = (κ_i / 2h_i) |t2 − t1|²` with `h_i` the **rest** Voronoi length `(L0_prev + L0_next)/2`.
Because `h` is the rest length, a uniformly stretched straight rod has exactly zero bending energy
and exactly zero bending force, so bending and stretching are decoupled and a buckling measurement is
a statement about `κ` alone. Small-angle limit: `|t2 − t1| = 2 sin(θ/2) → θ → C h`, so
`E_i → (κ/2) C² h` and the sum tends to `∫ (κ/2) C² ds` — Euler–Bernoulli with coefficient `κ` and no
leftover factor.

**The consequence that makes the bending law falsifiable.** Linearising a pinned chain about its
straight compressed state with mode `y_i = sin(π x_i / L)`:

- bending: `Σ (κ/2h)(y_{i+1} − 2y_i + y_{i−1})²/h² = (8κ/h³) sin⁴(πh/2L) Σ y²`
- geometric work of the compressive load: `(2F/h) sin²(πh/2L) Σ y²`

Equating gives the **discrete** critical load `F_c = (4κ/h²) sin²(πh/2L)`, which tends to
`π²κ/L²` from below at order `h²` (relative deficit `−π²h²/12L²`). Both numbers are asserted, because
checking only the continuum limit lets a wrong `h` convention hide inside the discretisation error,
and checking only the discrete formula lets a wrong continuum coefficient survive.

**Crosslink.** `E = Σ_c (k_c/2)(d_c − r_c)²` between two interpolated material points
`x = (1−w) p[a] + w p[b]`, one on each of two fibres. The gradient is scattered to all four nodes
with the same weights that built the two points — that is the condition under which the chain rule
closes, so `F = −grad E` holds for an element attached *between* nodes and the four forces sum to
zero. Two-sided by construction: a crosslinker resists separation and resists interpenetration, and
the tension-only discipline belongs to the fibre cables, not here.

## 4a. Mutation testing of the controls

Three mutants were applied to `aleph/vertical/ecm.py`, the suite run against each, and the source
restored and verified with `git diff`. Recorded here rather than in prose because a control suite
that no mutant can break is decoration.

| # | Mutation | Tests killed | Notes |
|---|---|---|---|
| M1 | `bending`: rest Voronoi `h` → **current** segment lengths | **3** | per-term FD (bending), total-force FD, and the scale-invariance control this exercise caused to be written |
| M2 | `axial`: drop the `/ rest` normalisation (`k Δ²` for `k Δ²/L0`) | **8** | all four `test_the_analytic_equilibrium_is_an_equilibrium` cases, resolution-independence, two buckling loads, buckling convergence |
| M3 | `crosslink`: scatter the `b` side with the `a` weights | **4** | per-term FD (crosslink), total-force FD, per-term moment closure, total moment closure |

All three killed, none survived. Source restored and verified **byte-identical** — the SHA-256 of the
restored file equals that of the pre-mutation copy, and `git diff` on the module shows only this
session's `__all__` and docstring repair (§12).

**The exercise found a real gap, which is the reason to run it.** M1 was originally killed by only
**two** tests, both finite-difference checks, and by neither buckling control nor the straight-rod
stretch control. The reason is specific: the buckling base state is straight and compressed by only
`O(F/k) = 2.5e-05`, so rest and current Voronoi lengths differ there by far less than that test's
tolerance; and a straight rod carries zero bending energy under either convention, so that control
was comparing zero with zero.
`::TestTheBendingLawIsTheOneClaimed::test_a_bent_rod_stretched_uniformly_keeps_exactly_its_bending_energy`
was added to close it, and it pins a property that matters well beyond this mutant — scale invariance
of the discrete bending energy is what `PLAN.md` §10 shows the dilational-virial identity rests on.

M2 is the counterpart lesson in the other direction: dropping the rest-length normalisation leaves
`F = −grad E` **exactly true**, so no finite-difference check can see it, and every one of its eight
kills came from the two independent closed forms. A gradient control and an oracle are not
substitutes for each other.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| node position | µm | 1e-6 m | finite |
| `axial_modulus_pn` | pN | 1e-12 N | `> 0` (stiffness × length) |
| `bending_rigidity_pn_um2` | pN·µm² | 1e-24 N·m² | `> 0` for a rod, exactly `0` for a cable |
| `stiffness_pn_per_um` | pN/µm | 1e-6 N/m | `>= 0` |
| damage `d` | dimensionless | — | `[0, 1)` |
| energy | pN·µm | 1e-18 J | `>= 0` |

Singular and boundary cases, each with the behaviour Aleph requires:

- **Collapsed segment** (`L = 0`): refused with `ValueError`. The direction is genuinely undefined.
- **Zero-length crosslink**: refused. A silent epsilon would return a finite force whose direction
  came from the round-off in the last subtraction.
- **`rest_length <= 0`**: refused; the normalisation would divide by zero.
- **`damage >= 1` or `< 0`**: refused; a fully damaged segment has no length scale left.
- **Cable with non-zero bending rigidity, or rod with zero**: refused at card construction. Each is a
  third element type that neither branch of the code implements.
- **Empty element set**: returns exactly `(0.0, zeros_like(positions))`, not a crash.
- **Self-crosslink** (`fibre_a == fibre_b`): refused with `ECMTopologyError`.

Invariants, each with the test that asserts it:

- **I1.** `F = −grad E` per term, order 2 —
  `tests/vertical/test_ecm_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order`
- **I2.** Force and moment closure against the constituent scale —
  `::TestForceAndMomentClosureOnTheIsolatedNetwork::test_each_term_closes_on_its_own`,
  `::test_each_term_carries_no_net_moment`
- **I3.** Translation invariance of the total energy — `::test_translating_the_whole_network_changes_no_energy`
- **I4.** Resolution-independent axial modulus — `::TestTheAxialLawAgainstItsClosedForm::test_the_extension_is_independent_of_the_node_count`
- **I5.** No term mutates its inputs (bit-identical) — `::TestTheTermsArePureFunctions::test_no_term_mutates_its_arguments`
- **I6.** Declared state keys equal the registered contract — `::TestTheDeclarationsAgreeWithTheCode::test_the_owned_state_keys_are_exactly_the_registered_contract`
- **I7.** No state key of this owner is declared by any other registered owner —
  `tests/state/test_census_key_ownership.py::test_no_state_key_has_two_registered_owners`

### Two of this owner's state keys collided with the cortex, and were renamed

Recorded here 2026-07-30 because it was missing: this entry described half of a two-owner defect and
said nothing about it, while `ALEPH-PORT-2301` described the other half and said the opposite of
what was true.

The ECM declared `crosslink_topology` and `topology_epoch`. So did the registered `cortex` contract.
Two owners, one key, twice over — which no gate checked, because `ComponentContract` only refuses a
name repeated inside a *single* contract. Nothing was broken at runtime, since no assembly puts both
owners against one key-addressed store; a store built from the registry later would have handed the
cortex and the ECM one crosslink table and one epoch counter.

Renamed to **`ecm_crosslink_topology`** and **`ecm_topology_epoch`**, with the cortex taking
`cortex_*` on its side. The reasoning, the alternative that was rejected, and the reversal procedure
are in `docs/decisions/PROPOSAL-cortex-ecm-state-key-collision.md`. I7 above is the guard that stops
this class of defect from needing to be discovered by a lane tripping over it.

## 6. Source evidence class and known retractions

No claim is made about the reference implementation's ECM, because none of it was read. Its
`ffn_sim/ac/engine/ecm_mechanics.py` and `ecm_world.py` are unaudited by this lane and their evidence
class is **unknown here**. Where I looked: nowhere inside them — deliberately, so that this module's
derivation could not be contaminated by theirs.

The registered `ecm` contract this module was written against carries `citation_status: UNSOURCED`,
and nothing here changes that. `default_collagen_card()` returns placeholder magnitudes and says so.

## 7. Independent oracle or derivation

Three, none of which is the reference implementation:

1. **Central finite differences** of the module's own energies — an oracle for the gradient claim
   that is independent of every constant in the law.
2. **The cable end-load closed form** `ΔL = F L0 / k`, derived in §4 from the energy alone. This one
   catches a defect the finite difference cannot: a rest-length normalisation error keeps `F = −grad E`
   exactly true while making the modulus resolution-dependent (mutant M2 above).
3. **Euler buckling**, `F_c = π²κ/L²` in the continuum and `(4κ/h²) sin²(πh/2L)` for this
   discretisation. A critical load is an eigenvalue of the whole operator, so an accidental rescaling
   of one term cannot reproduce it the way it can reproduce a single energy value.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_ecm_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order` | Per term, `energy > 0` first, then max &#124;F + ∇E&#124; / max&#124;F&#124; at steps 4e-4/2e-4/1e-4 µm with observed order > 1.6 |
| Positive | `tests/vertical/test_ecm_controls.py::TestTheBendingLawIsTheOneClaimed::test_the_critical_load_matches_the_discrete_prediction` | Measured buckling load equals `(4κ/h²)sin²(πh/2L)` to 2e-4 relative at N = 5, 9, 17 |
| Positive | `tests/vertical/test_ecm_controls.py::TestTheAxialLawAgainstItsClosedForm::test_the_analytic_equilibrium_is_an_equilibrium` | The uniform-tension configuration is an equilibrium and the end node carries exactly the applied load |

**Measured, this session, on this tree** (`/Users/sw1/miniconda3/envs/aleph/bin/python`, numpy 2.5.1,
scipy 1.18.0). Every number below was produced by a run in this session; none is quoted from memory.

Finite-difference agreement, steps 4e-4 / 2e-4 / 1e-4 µm:

| Term | Energy (pN·µm) | Force scale (pN) | max&#124;F+∇E&#124; at 1e-4 (pN) | relative | observed order |
|---|---|---|---|---|---|
| axial | 14.547974 | 63.16 | 1.493e-06 | 2.364e-08 | **2.000** |
| bending | 6.583502 | 20.93 | 3.271e-07 | 1.563e-08 | **2.000** |
| crosslink | 135.128836 | 50.25 | 1.707e-09 | 3.397e-11 | **1.945** |
| total | 156.260313 | — | — | 2.000e-08 | — |

Euler buckling, `L = 4 µm`, `κ = 20 pN·µm²`, axial modulus `5e5 pN`, continuum `π²κ/L² = 12.337006 pN`:

| N | h (µm) | measured F_c (pN) | discrete prediction (pN) | vs discrete | vs continuum |
|---|---|---|---|---|---|
| 5 | 1.0000 | 11.716255 | 11.715729 | 4.489e-05 | 5.0316 % |
| 9 | 0.5000 | 12.179844 | 12.179275 | 4.674e-05 | 1.2739 % |
| 17 | 0.2500 | 12.298001 | 12.297421 | 4.721e-05 | 0.3162 % |

Observed order of convergence to the continuum value: **1.996**. The residual 4.7e-05 disagreement
with the discrete formula is stable across all three resolutions and is the `O(F/k) = 2.5e-05`
compression of the base state that the linearised formula neglects — it is not a discretisation
error, which is why it does not shrink with `h`.

Closure on the isolated network (resultant ÷ constituent scale):

| Term | force | moment |
|---|---|---|
| axial | 4.47e-18 | 1.02e-17 |
| bending | 2.06e-17 | 2.78e-17 |
| crosslink | 2.99e-17 | 9.27e-17 |
| total | **4.98e-18** | **2.56e-17** |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_ecm_controls.py::TestTheNegativeControlsFail::test_a_flipped_axial_sign_is_caught_by_the_gradient_test` | With `flip_axial_sign=True` the FD disagreement rises above **0.5 relative**; measured 2.0 |
| Negative (must fail) | `tests/vertical/test_ecm_controls.py::TestTheNegativeControlsFail::test_a_cable_carrying_compression_is_visible_in_the_tension_it_reports` | With `allow_cable_compression=True` a cable at half its rest length stores energy and carries force; with the gate in place the same nodes carry **exactly** `0.0` |

Both breaks are flags on the shipped function, not hand-edited copies, so the control drives the real
code path. The second is asserted in **both** directions on purpose: asserting only that the healthy
branch is zero would pass for a function that returns zero unconditionally, which is the shape of the
`PLAN.md` §6.1 defect where a cold cache read made a deliberately broken strut look innocent.

## 10. Numerical and precision envelope

Working and accumulation precision: float64 throughout; no reduced-precision path exists.

Finite-difference steps 4e-4 / 2e-4 / 1e-4 µm on coordinates of order 1–6 µm. The floor for a central
difference is `eps^(1/3) ≈ 6e-6` **relative**; the smallest step used is 17× above it. Below that
floor cancellation in the energy evaluation dominates, the apparent order goes negative, and a
correct gradient is indistinguishable from a wrong one — the afternoon `PLAN.md` §6.1 records.

Tolerances and why they are not looser: the per-term relative FD agreement is asserted at `1e-7`
while the measured worst case is `2.4e-08`, a factor of 4 of headroom — tight enough that mutant M3
(a weight swap, which perturbs the gradient at the 1e-1 level) cannot pass. Closure is asserted at
`1e-14` of the constituent scale against a measured `5e-18`; the gap is deliberate, because closure
is an identity of the algebra and its residual is set by summation order, not by the physics. The
buckling load is asserted at `2e-4` relative against the discrete prediction, which is 4× the
systematic `4.7e-05` base-state offset documented above — assert it tighter and the test fails for a
reason that is understood and is not a defect.

Outside the envelope: every singular case in §5 refuses with an exception rather than degrading
silently. There is no input range over which this module returns a plausible number it cannot stand
behind, **with one exception that is not the module's fault**: the damage array is an input, and
nothing in this module or anywhere else produces one, so a caller supplying damage is supplying an
unvalidated model of damage.

## 11. Production-backend residency and transfer

Host, numpy, float64. Nothing here touches a backend, a device, or `aleph/runtime/`'s `Backend`
protocol — the module imports `EntityCensus`, `Phase`, `StepContext` and `WorkReceipt` from
`aleph.runtime.participant` and **uses none of them**, which is a leftover of the owner class that was
never written (see §14). No GPU work of any kind was run by this lane. When an owner does land, these
terms are `O(S)` and `O(T)` gather/scatter kernels with no host round-trip per step.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered this module, because no reference file was read.
`tests/vertical/test_ecm_controls.py::TestTheDeclarationsAgreeWithTheCode::test_no_reference_project_vocabulary_survives`
asserts the module text carries none of `ffn_cellsim`, `ffn_sim`, `dcm_contact`, `surface_body`, and
`::test_the_module_imports_nothing_from_validation` asserts the `aleph/** → validation/**` firewall.

Two pieces of **Aleph's own** prose were discarded, and they are the interesting ones:

- `__all__` exported five names that did not exist — `CollagenNetwork`, `FarFieldAnchor`,
  `WorldBoundary`, `assert_no_fibre_cable_carries_compression`, `assert_world_boundary_owns_no_state`.
  `from aleph.vertical.ecm import *` raised `AttributeError`, and to a reader scanning the export list
  each name read as shipped code. Removed rather than stubbed: an absent name is a visible gap and a
  stub is an invisible one. Replaced by
  `::TestTheDeclarationsAgreeWithTheCode::test_every_exported_name_resolves`, which now pins the list
  in both directions.
- The module docstring described all five in the **present tense**, including the accepted-step commit
  discipline for crosslink topology and the far-field anchor's reaction. Replaced by an
  `IMPLEMENTATION STATUS` block that states what exists, what does not, and that the module wires zero
  connectors. The design sections are kept, relabelled as the design for the follow-on lane.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 60 controls pass (`tests/vertical/test_ecm_controls.py`, 2026-07-30, 0.40 s). Three mutants applied, killing 3 / 8 / 4 tests respectively; none survived; source restored byte-identically. Status stays `PROPOSED` because the registered contract's owner does not exist — see §14. |
| Reviewer | Agent-proposed (lane L28). Unratified. No PI review. |
| Rollback | Delete `aleph/vertical/ecm.py` and `tests/vertical/test_ecm_controls.py`. Nothing breaks: no module in the tree imports either, and no connector reaches them. |

## 14. Honest limits — what this entry does NOT establish

**The ledger entry was written after the code.** PLAN §0.2.5 requires it before. The lane died on an
API error mid-flight and this entry was reconstructed afterwards; that is the reason and it is not an
excuse, since a ledger written after the fact cannot have constrained the design it describes.

**Wired connector count contributed by this module: zero.** Not one endpoint of one connector is
implemented. `membrane_ecm_contact`, `integrin_collagen_clutch`, `ecm_crosslink` and
`ecm_far_field_anchor` all still have an absent `ecm` endpoint. The lane brief's target of six
connectors gated by `ecm` + `extracellular_medium` is **not** advanced by anything here. (Measured
against the registry: five connectors name `ecm` or `extracellular_medium` as an endpoint —
`membrane_medium_traction`, `membrane_ecm_contact`, `integrin_collagen_clutch`, `ecm_crosslink`,
`ecm_far_field_anchor` — not six. Where the sixth comes from is unresolved and is reported rather
than guessed.)

**`extracellular_medium` does not exist at all.** No `aleph/vertical/medium.py` was written. So there
is no exterior Stokes operator, no `6πμRV` drag comparison, no convergence under surface refinement,
no rigid-body-mode regularisation, and no medium snapshot/rollback control — including the negative
control the brief asked for, in which dropping the medium's rollback leaves it advanced after a
rejected step. That was the whole second half of the lane and none of it happened.

**No state, therefore no transaction evidence.** The registered contract declares seven state blocks
(`fibre_node_positions_um`, `fibre_material_coordinates`, `crosslink_topology`, `topology_epoch`,
`ligand_coordinates`, `fibre_damage_state`, `remodelling_state`) and this module allocates **none** of
them. `owned_state_keys()` agrees with the registry, which means the *declaration* is right and the
*implementation* is absent — precisely the configuration that reads as complete from the registry
side. Specifically unestablished:

- **Rollback restoring bit-identically.** There is no state to roll back. The nearest available
  control is that the energy terms are pure functions and do not mutate their inputs
  (`np.array_equal`, not `allclose`), which is necessary for a correct rollback and nowhere near
  sufficient.
- **Crosslink topology commits only on an accepted step.** The commit discipline is described in the
  docstring and implemented nowhere. `Crosslink` is a frozen dataclass with validation; there is no
  propose/commit/rollback and no epoch counter.
- **Ligands survive remodelling.** The design argument is in the docstring and is, I believe, correct
  — material coordinates in a reference frame fixed at construction cannot move when rest lengths are
  rewritten. It is **untested**, because there is no remodelling operation to survive and no resolver
  to resolve a ligand against a network. `LigandSite` is validated and inert.
- **Damage evolution.** The axial law consumes a damage array and is verified to apply `(1 − d)`
  exactly. Nothing produces a damage array, so the irreversibility, the monotonicity and the
  strain envelope in `FibreCard` are declarations only.
- **The far-field anchor and `world_boundary`.** No clamp, no reaction, no frame object. The claim
  that the reaction transferred to the frame equals the negative of what the constraint removed, and
  the global-equilibrium identity that would check it, are both unimplemented and unmeasured.
  `WorldBoundaryStateError` is a defined exception that nothing raises.

**The controls test the elements a network is made of, not a network.** `_MixedNetwork` in the test
file assembles segments, triples and crosslinks by hand. That scaffolding is the test's own, and its
existence is itself a finding: an element library with no owner obliges every consumer to reassemble
the network by hand, and two consumers that do it differently are two different models.

**Physical magnitudes are not established.** Evidence rung `ANALYTIC_ORACLE`, quantitative status
`BLOCKED`. No modulus, rigidity, crosslink stiffness or damage threshold here is a collagen number;
they are chosen to make the algebra sharp. `citation_status` for the `ecm` contract remains
`UNSOURCED` and this entry does not improve it.

**The reference implementation is unaudited by this lane.** Whether Aleph's re-derivation agrees or
disagrees with `ffn_sim/ac/engine/ecm_mechanics.py` is unknown, and this entry deliberately does not
guess.

**Not run:** any GPU job (zero, and no authorization was sought or held), and the full suite. Only
`tests/vertical/` was exercised — five other lanes were writing concurrently, and PLAN §1 says to
report a foreign breakage rather than fix it, which means not conflating one with mine.
