# ALEPH-PORT-3203 — intermediate filament as an explicit strain-stiffening, tension-only cable graph

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3203` |
| Lane | `L32 vertical owners — internal frame` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (**before** the code, per PLAN §0.2.5 — see §13 for the one section that could not be) |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Zero lines taken, no reference file opened. `ALEPH-PORT-1802` already established that the nonlinear law this component's contract declares **does not exist** in the reference, so there is nothing there to port even in principle. |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.vertical.intermediate_filament import (
    IF_ENDPOINT_ROLES,
    OWNED_STATE_KEYS,
    PUBLISHABLE_ENDPOINT_ROLES,
    ROLE_POPULATIONS,
    Cable,
    IFAttachmentRole,
    IFAttachmentSite,
    IFMaterialCard,
    IFOwnershipError,
    IFResponseClaim,
    IFRoleError,
    IFStaleHandleError,
    IFUnsourcedCardError,
    IntermediateFilamentNetwork,
    InternalCrosslink,
    MaterialPoint,
    MaterialPointSink,
    Population,
    assert_populations_own_no_state,
    assert_state_keys_disjoint_from,
    assert_strain_limit_is_bound,
    cable_axial_energy_and_forces,
    crosslink_energy_and_forces,
    owned_state_keys,
    storage_energy_density_pn,
    straight_cable_nodes,
    tangent_modulus_pn,
    tension_pn,
)
```

This **is** a state owner, unlike `ALEPH-PORT-2801`'s element library: it allocates all six registered
state blocks, carries the accepted-step transaction, and publishes endpoints. It wires **zero**
connectors — publishing an endpoint is not wiring one, and `aleph/vertical/wiring.py` is untouched by
this entry.

There is deliberately **no `default_if_card()`**, on the pattern `ALEPH-PORT-2401` established for
`PoroelasticCard`, and a control asserts the name does not exist. §4 says why that is the whole
argument of this module rather than a stylistic choice.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — the commit `ALEPH-PORT-1802` audited and verified to resolve |
| Source path | `ffn_sim/ac/engine/intermediate_filament_rig.py`, `ffn_sim/ff/intermediate_filaments.py`, `ffn_sim/ac/solid/intermediate_filament.py` — **named, not read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened for this module or for this entry. |
| Working tree == commit? | not applicable — nothing was read, so there is no revision to have mis-cited |

The three paths are recorded so a later auditor knows where the reference's cage lives and can
establish independently whether Aleph's re-derivation agrees with it. Listing them is **not** a claim
that they were consulted. Everything this entry says *about* the reference is quoted from
`ALEPH-PORT-1802`, which did read them, and is attributed there rather than re-asserted here.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** `ALEPH-PORT-1802` §3 settled this with a fact rather than a
preference: the reference's strain-stiffening branch is a named-but-unfilled regime whose exported
constant is the string `"WLC_STRAIN_STIFFENING_UNFOLDING_CARD_REQUIRED"`, and what actually runs
there is the **linear tangent** limit — a Hookean spring. The single thing this module exists to
supply is therefore the single thing the reference does not have. A Hookean spring is not worth
porting, and a seam over borrowed device kernels written against a connector-contract system Aleph
deliberately does not share has negative porting value.

The law in §4 was written from the registered `intermediate_filament` contract in
`aleph/state/census_frame_fluid_nucleus.py:407` and from the three connector contracts in
`aleph/state/connectors_crosssystem_cytosol.py`, both of which are Aleph's own documents.

**Recorded as prior art, unported** (all three inherited from `ALEPH-PORT-1802`, not re-established
here): keratin and vimentin as two material cards under one owner; the radial-spoke-versus-tangential-shell
negative result about cage topology; and the reference's refusal to fit strain-stiffening parameters
in order to make a coupling appear. The third is the stance this module implements mechanically.

## 4. Physical or mathematical law represented

Two energy terms, each with its exact analytic gradient. Units µm / pN / pN·µm throughout.

### 4.1 The tension-only strain-stiffening cable

The registered representation is an "explicit nonlinear strain-stiffening cable graph … soft at small
strain and stiffens at large strain, so **the tangent stiffness is state and not a parameter**". The
registered `unsupported_claims` forbids claiming that response *at any specific strain* until the
material card fixing "the knee, the plateau and the re-stiffening slope" is sourced, per population.

So the law is written with exactly those three features and the card carrying them is a **required
constructor argument with no default**. On segment strain `e = (L − L0)/L0`:

```
                 0                                        e ≤ 0        (slack: a cable does not push)
T(e) =   E_toe · e                                        0 < e ≤ ek   (toe)
         T_k + E_pl · (e − ek)                            ek < e ≤ ep  (plateau)
         T_p + E_st · (e − ep)                            e > ep       (re-stiffening)

T_k = E_toe · ek                    T_p = T_k + E_pl · (ep − ek)
```

with `E_toe > 0`, `0 ≤ E_pl < E_toe` (a plateau is *softer* than the toe or it is not a plateau), and
`E_st > E_toe` (the law must actually stiffen or the representation is false). All three are refused
at card construction, not clamped.

The stored energy of a segment is `E_seg = L0 · W(e)` with `W(e) = ∫₀^e T(s) ds`:

```
W(e) = ½E_toe·a² + T_k·b + ½E_pl·b² + T_p·c + ½E_st·c²
a = min(e⁺, ek)      b = clip(e − ek, 0, ep − ek)      c = max(e − ep, 0)      e⁺ = max(e, 0)
```

Then `dE_seg/dL = L0 · W′(e) · (1/L0) = T(e)` exactly, so the axial tension **is** the derivative of
the energy with respect to length and the force on the two ends is `±T · t̂`, positive tension pulling
them together. Two consequences that the controls turn into assertions:

* **Resolution independence.** The energy carries the factor `L0` and the strain is dimensionless, so
  a cable cut into twice as many segments has twice as many springs each twice as stiff. A cable of
  total rest length `L0` under end load `F` extends by the *same* amount however it is discretised —
  a closed form no finite-difference check can see (`ALEPH-PORT-2801` §4a mutant M2 is the same
  lesson).
* **The tangent is state.** `dT/de` is `E_toe`, `E_pl` or `E_st` according to where the segment
  currently sits, so two segments cut from one card report different stiffnesses. That is the
  registered claim, and it is asserted directly rather than inferred from a force.

`W` is C¹ and piecewise quadratic; the tangent is piecewise constant with jumps at `ek` and `ep`.
That is deliberate — it is the minimal law with the three declared features and no fourth invented
parameter — and it has a consequence for the controls that §10 states: a finite-difference probe must
not straddle a knee, so every gradient fixture is built with each segment at least `5e-3` in strain
away from both.

### 4.2 The internal crosslink

`E = Σ_c (k_c/2)(d_c − r_c)²` between two interpolated material points, one on each of two
**different** cables (state key `internal_crosslink_topology`). Two-sided, for the reason
`ALEPH-PORT-2801` gives: the tension-only discipline belongs to the cables, where it is the difference
between a load path and a fiction; applied to a crosslinker it would let two crosslinked cables pass
through each other. The gradient is scattered to all four nodes with the *same* weights that built
the two points, which is the condition under which the chain rule closes.

### 4.3 The law that is deliberately absent

**No bending.** The registered representation is a cable graph and a cable carries tension only, so
there is no bending stencil, no persistence length and no Euler buckling here. `sf_arc` and `ecm`
carry bending because their contracts declare rods; this one does not.

**No active channel.** Nothing in this owner generates tension of its own. `active_power_pn_um_per_s()`
returns exactly `0.0`.

**No turnover kinetics.** `propose_rebinding` is called by something outside this owner; no rate
constant, no RNG and no unbinding law exists here. §14 says what that costs.

### 4.4 Why `strain_history` is a channel that carries physics

The registered contract's own words: "a strain-stiffening cable with turnover is history-dependent —
a rest length that is re-set on rebinding is state, and omitting it would make the law
path-independent by accident."

Implemented as two rest-length frames that are *both read*:

* **The reference frame**, fixed at construction, is what `material_coordinates` is measured in. An
  attachment site is `(cable_id, s)` with `s` a reference rest arclength, so a site resolves to the
  same node pair with the same interpolation weight before and after turnover. An attachment that
  moved when the network remodelled would be an attachment that had silently let go.
* **The current frame** is `segment_rest_um`, and it is what the law in §4.1 measures strain against.
  A committed rebinding re-sets the rest lengths of that cable's segments **to their current
  lengths**, so the accumulated strain is released and the subsequent load path differs.

The consequence is testable and is tested: two networks brought to the *same* geometry by two
different paths — one of which committed a rebinding en route — report different tension. With the
shipped break `freeze_strain_history=True` the two paths agree exactly, which is the contract's
"path-independent by accident" reproduced on purpose so the control drives the shipped code path.

### 4.5 The defect this module exists not to have

The registered `mechanical_role` is that this owner is "the **strain-limiting** one: soft over the
working range and stiff past its knee, so it sets a **ceiling on nuclear deformation** rather than a
resting force."

A cable graph with a linear axial law has **no ceiling**. It is a perfectly conservative model: the
gradient control passes, force and moment closure pass exactly, energy is conserved exactly, and the
nucleus can be deformed arbitrarily at quadratic cost. This is the same shape of defect as
`ALEPH-PORT-2301`'s missing axial law — *a missing energy term conserves energy perfectly* — and no
gradient, conservation or closure check can see it.

`assert_strain_limit_is_bound` is the check that can. Per cable, using that cable's **own** card
(the card is per population, so the probe strains are too), it evaluates the network's axial energy at
three equally spaced strains inside the declared plateau and at three inside the declared
re-stiffening branch, and takes the second difference of each triple:

```
tangent_window = [E(x − d) − 2E(x) + E(x + d)] / d²
```

Within a window that lies entirely inside one linear-tangent branch this is *exactly* that branch's
modulus times the cable's total rest length, because `W` is quadratic there. The ratio
`tangent_high / tangent_low` is therefore `E_st / E_pl` for a bound law and **exactly 1.0** for a law
whose re-stiffening branch is missing. Three details are load-bearing:

* **The second difference, not the energy or its first difference.** A second difference annihilates
  any constant and any linear term, so the other cables' unchanged energy and the network's absolute
  energy offset both drop out, and the number measured is a curvature rather than a level.
* **Per cable, with per-cable probe strains**, because a whole-network probe would let one stiff cable
  hide a network of cables that never stiffen, and because a keratin knee is not a vimentin knee.
* **Through the network's own `axial_term()`**, not through `tension_pn` directly, so the guard fails
  if the law is correct but unbound — which is exactly the reference cortex's failure mode.

The probe positions are laid out by putting each cable straight along its own current end-to-end
direction with segment lengths `(1 + e)·L0`, giving every segment *exactly* the target strain, and the
positions are restored bit-identically afterwards.

## 4a. Mutation testing of the controls

Four mutants applied to `aleph/vertical/intermediate_filament.py`, the control suite (90 tests) run
against each, then the source restored and verified. Recorded here rather than in prose because a
control suite that no mutant can break is decoration.

| # | Mutation | Tests killed | Survived? |
|---|---|---|---|
| M1 | The shipped `omit_restiffening_branch` default flipped to `True`: the re-stiffening branch simply absent from the shipped owner, and the law still a **perfect gradient of its own energy** | **6** | no |
| M2 | `_branch_widths`: the re-stiffening extent `c` zeroed while the plateau extent `b` stays clipped at `ep` | **8** | no |
| M3 | `cable_axial_energy_and_forces`: drop the `L0 ·` factor on the energy (`W` summed directly), leaving `F = −grad E` **false** by one factor of `L0` per segment | **7** | no |
| M4 | `commit`: rebinding writes the *reference* rest lengths back instead of the current lengths, so turnover is a no-op that still bumps the epoch and the counter | **3** | no |

All four killed, none survived. The module was restored from a pre-mutation copy and verified
**byte-identically**: SHA-256 `6ae7027e07b1395200dce441bf2cb04b5e11b7c471d0f5647db3a490cc105c62`
before and after. (`git diff` cannot witness it — the file is untracked, because this lane stages
nothing.) Measured numbers are in §8; the two findings the exercise produced are in §13.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| node position | µm | 1e-6 m | finite |
| segment rest length `L0` | µm | 1e-6 m | `> 0` |
| strain `e` | dimensionless | — | `(−1, ∞)`; response identically zero for `e ≤ 0` |
| `toe_modulus_pn` `E_toe` | pN | 1e-12 N | `> 0` |
| `plateau_modulus_pn` `E_pl` | pN | 1e-12 N | `[0, E_toe)` |
| `restiffening_modulus_pn` `E_st` | pN | 1e-12 N | `> E_toe` |
| `knee_strain` `ek` | dimensionless | — | `> 0` |
| `plateau_end_strain` `ep` | dimensionless | — | `> ek` |
| tension `T` | pN | 1e-12 N | `>= 0` — a cable never reports a negative tension |
| crosslink `stiffness_pn_per_um` | pN/µm | 1e-6 N/m | `>= 0` |
| energy | pN·µm | 1e-18 J | `>= 0` |

Singular and boundary cases, each with the behaviour Aleph requires — refusal, never a clamp and
never a regularisation:

- **Slack cable** (`e ≤ 0`): energy and force are **exactly** `0.0`, testable with `==`, because
  `dW/de` is identically zero on the whole slack branch rather than merely small near the gate.
- **Collapsed segment** (`L = 0`): refused with `ValueError`. The segment direction is undefined.
- **`rest_length <= 0`**: refused; the strain would divide by zero.
- **Zero-length crosslink**: refused. A silent epsilon returns a finite force whose direction came
  from the round-off in the last subtraction.
- **Self-crosslink** (`cable_a == cable_b`): refused with `IFRoleError`. Between two points of one
  cable it is a second, undeclared axial law with its own rest length.
- **Card with `E_st <= E_toe`**: refused at construction. That card describes a law that does not
  stiffen, and calling it strain-stiffening would be the claim the contract forbids.
- **Card with `E_pl >= E_toe`**: refused. A plateau stiffer than the toe is not a plateau.
- **Card with `ep <= ek`**: refused; the plateau would have negative width.
- **A population with no card**: refused at construction. There is no shared card and no default one.
- **Reporting a nonlinear response from an unsourced card**: refused with `IFUnsourcedCardError`,
  naming the linear tangent limit as the honest default and the populations that are missing a
  citation.
- **Material coordinate outside its cable**: refused rather than clamped; clamping moves a load path
  without telling anybody.
- **A handle used after a committed topology or turnover change**: refused with `IFStaleHandleError`.
- **Empty element set**: returns exactly `(0.0, zeros_like(positions))`, not a crash.

Invariants, each with the test that asserts it. All test paths are in
`tests/vertical/test_intermediate_filament_controls.py` unless stated.

- **I1.** `F = −grad E` per term, observed order > 1.6 —
  `::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order`,
  and separately in each branch of the constitutive law —
  `::test_the_gradient_holds_in_every_branch_of_the_constitutive_law`
- **I2.** Force and moment closure against the **constituent** scale, per term and in total —
  `::TestForceAndMomentClosureOnTheIsolatedNetwork::test_every_term_closes_on_its_own`,
  `::test_the_internal_moments_sum_to_zero`
- **I3.** Translation invariance of the total energy —
  `::TestForceAndMomentClosureOnTheIsolatedNetwork::test_translating_the_whole_network_changes_no_energy`
- **I4.** No term mutates its inputs, bit-identically (`np.array_equal`) —
  `::TestForceAndMomentClosureOnTheIsolatedNetwork::test_no_term_mutates_its_arguments`
- **I5.** Resolution-independent extension under an end load —
  `::TestTheConstitutiveLawIsTheOneDeclared::test_the_extension_under_an_end_load_does_not_depend_on_the_node_count`
- **I6.** The tangent stiffness is state, not a parameter —
  `::TestTheConstitutiveLawIsTheOneDeclared::test_the_tangent_modulus_is_state_and_not_a_parameter`
- **I7.** A rebinding changes the subsequent load path —
  `::TestStrainHistoryIsPathDependent::test_a_rebinding_event_changes_the_subsequent_load_path`
- **I8.** Material coordinates survive turnover bit-identically —
  `::TestStrainHistoryIsPathDependent::test_material_coordinates_survive_a_rebinding`
- **I9.** Rollback restores bit-identically and a rejected rebinding leaves no trace —
  `::TestOwnershipAndTheAcceptedStepContract::test_rollback_restores_bit_identically`,
  `::TestStrainHistoryIsPathDependent::test_a_rejected_rebinding_leaves_the_rest_lengths_untouched`
- **I10.** A committed change bumps the epoch and a stale handle is refused —
  `::TestOwnershipAndTheAcceptedStepContract::test_a_committed_topology_change_bumps_the_epoch`,
  `::test_a_handle_held_across_a_commit_is_refused_and_a_fresh_one_works`
- **I11.** The scattered load carries the force **and** the moment of the load applied —
  `::TestTheEndpointSink::test_the_scattered_moment_is_also_conserved`
- **I12.** Declared state keys equal the registered contract and collide with no other owner —
  `::TestTheDeclarationsAgreeWithTheCode::test_the_owned_state_keys_are_exactly_the_registered_contract`,
  `tests/state/test_census_key_ownership.py::test_no_state_key_has_two_registered_owners`
- **I13.** Neither population owns state and neither is a connector endpoint —
  `::TestOwnershipAndTheAcceptedStepContract::test_neither_population_owns_state_or_is_a_connector_endpoint`
- **I14.** The role enum equals the registered connector manifest for this owner —
  `::TestTheDeclarationsAgreeWithTheCode::test_the_role_enum_matches_the_registered_connector_manifest`
- **I15.** A work receipt declares non-zero ops and an unscheduled phase is refused —
  `::TestTheWorkReceiptCannotLieAboutHavingRun`

### A note on the six unprefixed state keys

All six registered keys — `cable_geometry`, `material_coordinates`, `strain_history`,
`population_identity`, `internal_crosslink_topology`, `turnover_state` — are **bare**, where `cortex`,
`ecm` and `sf_arc` all carry an owner prefix. They are reproduced verbatim; renaming module-side would
put the module and the registry into exactly the disagreement `ALEPH-PORT-2801` §5 documents.

**Measured, not assumed: none of the six currently collides.**
`tests/state/test_census_key_ownership.py` passes across the whole census (3 passed, this session), so
no other registered owner declares any of them today. `material_coordinates` is the exposed one — it
is a category name rather than an owner name, and `ecm` and `cortex` both hold the same concept under
`fibre_material_coordinates` and `filament_material_coordinates`. A future owner that declares
`material_coordinates` bare would collide, and that is now a *failing test* rather than a discovery,
which is the whole point of the census-wide guard. No proposal is filed because there is no live
collision to propose a resolution for.

## 6. Source evidence class and known retractions

No claim is made about the reference implementation's intermediate-filament code, because none of it
was read by this lane. Where I looked: nowhere inside it — deliberately, so this module's derivation
could not be contaminated. `ALEPH-PORT-1802` §6 is the audit of record, it found the source's own
status constants declaring the nonlinear branch unfilled and the cage Hookean, and it found one
retraction (the tangential perinuclear shell, superseded by radial spokes as non-coupling
decoration). This entry inherits those findings by citation and re-establishes none of them.

The registered `intermediate_filament` contract carries `citation_status: UNSOURCED`, and **nothing
here changes that.** Every magnitude in the tests is chosen to make the algebra sharp; not one is a
keratin or vimentin number. That is precisely why `IFMaterialCard.citation` is a required field with
no default and why `IFUnsourcedCardError` exists.

## 7. Independent oracle or derivation

Four, none of which is the reference implementation:

1. **Central finite differences** of the module's own energies, per term and per branch — an oracle
   for the gradient claim that is independent of every constant in the law.
2. **The closed-form energy** `W(e)` integrated by hand in §4.1, compared against a numerical
   quadrature of the module's own `tension_pn` over `[0, e]`. This catches an error in the energy
   that leaves `F = −grad E` intact — a mis-integrated branch offset shifts the energy by a constant
   in one branch and the finite difference never sees it.
3. **The end-load closed form** `ΔL = F·L0/E_toe` in the toe branch, which is resolution-independent.
   `ALEPH-PORT-2801` mutant M2 is the recorded proof that this catches what a gradient check cannot.
4. **The exact second difference of a piecewise-quadratic energy**, §4.5: inside one branch the
   window curvature is the branch modulus exactly, so the guard's expected value is an identity
   rather than a fit, and a broken law returns exactly `1.0` rather than something merely small.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_intermediate_filament_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order` | Per term, `energy > 0` first, then max&#124;F + ∇E&#124; / max&#124;F&#124; at 4e-4 / 2e-4 / 1e-4 µm with observed order > 1.6 |
| Positive | `::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_gradient_holds_in_every_branch_of_the_constitutive_law` | The same, separately, with every segment in the toe / plateau / re-stiffening branch |
| Positive | `::TestTheConstitutiveLawIsTheOneDeclared::test_the_law_is_the_exact_integral_of_its_own_tension` | `W(e)` equals a 200 001-point trapezoid of `tension_pn` at 1e-6 relative, at eight strains per card and two cards, spanning all four branches; measured worst case 1.272e-11 |
| Positive | `::TestTheStrainLimitIsActuallyBound::test_a_strain_stiffening_cable_stiffens_past_its_plateau` | Window curvature ratio equals `E_st/E_pl` to 1e-9 relative; measured 66.666667 (keratin) and 62.500000 (vimentin) against the cards' exact 66.666667 and 62.500000 |
| Positive | `::TestStrainHistoryIsPathDependent::test_a_rebinding_event_changes_the_subsequent_load_path` | Same final geometry, two histories, tensions differing by more than 3×; measured 12.653846× |

**Measured, this session, on this tree** (`/Users/sw1/miniconda3/envs/aleph/bin/python`, numpy 2.5.1).
Every number below was produced by a run in this session; none is quoted from memory.

Finite-difference agreement, steps 4e-4 / 2e-4 / 1e-4 µm, on the mixed fixture:

| Term | Energy (pN·µm) | Force scale (pN) | max&#124;F+∇E&#124; at 1e-4 (pN) | relative | observed order |
|---|---|---|---|---|---|
| axial | 3.189128e+02 | 1.0455e+03 | 3.774e-06 | 3.610e-09 | **2.000** |
| crosslink | 3.559358e+00 | 6.8893e+00 | 6.251e-08 | 9.073e-09 | **2.000** |
| total | 3.224721e+02 | 1.0455e+03 | 3.774e-06 | 3.610e-09 | **2.000** |

Per branch of the constitutive law, axial term only:

| Branch | strain | Energy (pN·µm) | relative FD residual at 1e-4 | observed order |
|---|---|---|---|---|
| toe | 0.030 | 7.292170e-01 | 9.581e-08 | **2.000** |
| plateau | 0.200 | 1.554650e+01 | 9.941e-10 | **2.000** |
| re-stiffening | 0.600 | 3.143735e+02 | 6.014e-09 | **2.000** |

Closure on the isolated network (resultant ÷ constituent scale):

| Term | force | moment |
|---|---|---|
| axial | **0.000e+00** | 4.000e-17 |
| crosslink | 8.537e-18 | 1.175e-17 |
| total | 8.800e-18 | 3.996e-17 |

Strain-limit guard, keratin card (`E_toe = 8.0e2`, `ek = 0.05`, `E_pl = 6.0e1`, `ep = 0.35`,
`E_st = 4.0e3`) and vimentin card (`E_toe = 5.0e2`, `ek = 0.09`, `E_pl = 4.0e1`, `ep = 0.50`,
`E_st = 2.5e3`):

| Cable | population | low-window curvature ÷ ΣL0 (pN) | high-window curvature ÷ ΣL0 (pN) | ratio | expected `E_st/E_pl` |
|---|---|---|---|---|---|
| 5 | keratin | 60.000000 | 4000.000000 | **66.666667** | 66.666667 |
| 2 | vimentin | 40.000000 | 2500.000000 | **62.500000** | 62.500000 |
| 9 | keratin | 60.000000 | 4000.000000 | **66.666667** | 66.666667 |

Each window curvature is the branch modulus **to every printed digit**, which is the §7 item 4
identity rather than a fit. `assert_strain_limit_is_bound` returns the minimum over the three cables,
62.500000. With `omit_restiffening_branch=True` the same probe returns **exactly 1.000000** on all
three, and the guard raises.

Path dependence of the load path across a rebinding, one keratin cable, rest length 4.0 µm:

| Path | strain at the shared final geometry | segment tension (pN) |
|---|---|---|
| no rebinding | 0.500000 | 6.580000e+02 |
| rebinding committed at strain 0.20 | 0.250000 | 5.200000e+01 |
| `freeze_strain_history=True`, same rebinding | 0.500000 | 6.580000e+02 |

The third row is exactly the first, bit for bit, which is the contract's "path-independent by
accident" made visible.

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `::TestTheNegativeControlsFail::test_a_flipped_axial_sign_is_caught_by_the_gradient_test` | With `flip_axial_sign=True` the FD disagreement rises above 0.5 relative; measured **2.000** |
| Negative (must fail) | `::TestTheNegativeControlsFail::test_a_cable_carrying_compression_is_visible_and_the_healthy_branch_is_exactly_zero` | With `allow_cable_compression=True` a cable at 0.7× rest length stores `7.292170e+01` pN·µm and every segment reports `-2.400000e+02` pN, i.e. it pushes; with the gate in place the energy is **exactly** `0.0`, every force is **exactly** `0.0` and every tension is **exactly** `0.0` |
| Negative (must fail) | `::TestTheStrainLimitIsActuallyBound::test_the_same_probe_is_free_of_stiffening_when_the_restiffening_branch_is_omitted` | With `omit_restiffening_branch=True` the guard's curvature ratio is exactly `1.000000` and `assert_strain_limit_is_bound` raises; with it bound the ratio is 62.500000 and the guard passes |
| Negative (must fail) | `::TestStrainHistoryIsPathDependent::test_the_frozen_history_break_makes_the_law_path_independent` | With `freeze_strain_history=True` the two histories give bit-identical tension (`6.580000e+02` pN both ways, compared with `==`); with it off they differ by 12.653846× |

All four breaks are flags on the **shipped** owner, default `False`, so each control drives the real
code path rather than a duplicate of it. A negative control against a duplicate proves only that the
duplicate is broken. Each is asserted in **both** directions — the broken branch misbehaves *and* the
healthy branch is exactly inert — because asserting only that the healthy branch is zero would pass
for a function that returns zero unconditionally.

## 10. Numerical and precision envelope

Working and accumulation precision: float64 throughout; no reduced-precision path exists.

Finite-difference steps 4e-4 / 2e-4 / 1e-4 µm on coordinates spanning 0.004–2.44 µm (measured on the
mixed fixture). The floor for a central difference is `eps^(1/3) ≈ 6e-6` **relative**; against the
segment length scale of ≈ 0.5 µm the smallest step is 2e-4 relative, about 30× above it. Below that
floor cancellation in the energy evaluation dominates, the apparent order goes negative, and a correct
gradient is indistinguishable from a wrong one — the afternoon `PLAN.md` §6.1 records. The measured
order is 2.000 at every one of the seven measurements in §8, which is the evidence that the steps are
in fact above the floor rather than the argument that they should be.

**The knee constraint, which is specific to this law.** `W` is C¹ but not C², so a central difference
that straddles `ek` or `ep` measures the tangent jump rather than the truncation error and its
observed order collapses to 1. Every gradient fixture is therefore built with each segment's strain at
least `5e-3` away from both knees, against a largest FD step of `4e-4` µm on segments of rest length
`≈ 0.5 µm`, i.e. a strain perturbation of `≈ 8e-4`. The margin is 6×, chosen rather than fitted, and
it is stated here because a later author who moves a fixture strain onto a knee will see an order of
1 and will not otherwise know why.

Tolerances and why they are not looser: the relative FD agreement is asserted at `1e-6` against a
measured worst case of `9.581e-08` (the toe-branch fixture, whose force scale is the smallest of the
seven measured), a margin of about 10×. That is tight enough that mutant M3 — a whole factor of `L0`
per segment — cannot pass; it did not, and it killed all four gradient controls. Closure is asserted at
`1e-13` of the constituent scale against a measured worst of `4.000e-17`; the gap is deliberate,
because closure is an identity of the algebra and its residual is set by summation order, not by
physics.  The axial term's force closure is measured at **exactly** `0.000e+00`, which is what a
segment law whose two scatters are `+pull` and `-pull` of the same array should give and is asserted
against the same tolerance as the rest rather than being special-cased. The strain-limit curvature ratio is asserted
against the card's own `E_st/E_pl` at `1e-9` relative, which is possible only because the window
curvature is an exact identity inside one branch (§7 item 4).

Outside the envelope: every singular case in §5 refuses with an exception rather than degrading
silently. There is no input range over which this module returns a plausible number it cannot stand
behind — **with one exception that is stated rather than fixed**: a card whose `citation` is `None`
still produces numbers. It must, or the law could not be tested at all. What it may not do is have
those numbers *reported as a strain-stiffening response*, and `assert_nonlinear_response_is_reportable`
is the refusal that enforces the difference.

## 11. Production-backend residency and transfer

Host, numpy, float64. The owner takes a `Backend` through `StepContext` and uses it for the
accumulator clear, the add and the overdamped update, exactly as `cortex_filaments` does; every energy
term is plain numpy on host arrays. No GPU work of any kind was run by this lane and no GPU
authorization was sought or held. When a device path lands, the axial term is an `O(S)` segment
gather/scatter and the crosslink term an `O(C)` four-node gather/scatter, neither needing a host
round-trip per step. The strain-history arrays are `O(S)` and are written only on an accepted step, so
they do not participate in the per-candidate-step traffic.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered this module, because no reference file was read.
`::TestTheDeclarationsAgreeWithTheCode::test_no_reference_project_vocabulary_survives` asserts the
module text carries none of `ffn_cellsim`, `ffn_sim`, `dcm_contact`, `surface_body`, `wlc`, `unfolding`, and
`::test_the_module_imports_nothing_from_validation` asserts the `aleph/** → validation/**` firewall.

Two pieces of vocabulary are deliberately **not** used even though they are tempting and would read as
precise:

- **"worm-like chain" / "WLC".** The law in §4.1 is a three-branch piecewise-linear tangent, not a WLC.
  Calling it a WLC would import a specific molecular picture, a persistence length, and a divergence at
  full extension that this law does not have, and would let a reader infer a parameterisation nobody
  registered.
- **"unfolding".** The plateau branch is a plateau. Attributing it to α-helix unfolding is a mechanism
  claim, and the registered `unsupported_claims` explicitly forbids subunit-level mechanics — "a cable
  graph resolves neither the coiled-coil dimer nor the unit-length-filament assembly step".

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 90 controls pass (`tests/vertical/test_intermediate_filament_controls.py`, re-run 2026-07-31 after a docstring reflow, 0.45 s). Four mutants applied, killing 6 / 8 / 7 / 3 tests; none survived; source restored byte-identically (SHA-256 unchanged). Status stays `PROPOSED`: no connector is wired, no material card is sourced, and no PI review has happened. |
| Reviewer | Agent-proposed (lane L32). Unratified. No PI review. No `decided_by` field appears in this entry by design. |
| Rollback | Delete `aleph/vertical/intermediate_filament.py` and `tests/vertical/test_intermediate_filament_controls.py`, and regenerate `ports/ledger/INDEX.md`. Nothing breaks: no module in the tree imports either, `aleph/vertical/__init__.py` is untouched, and no connector reaches them. |

**The one section written after the code, and why.** §4a and the measured tables in §8 and §9 record
numbers a run produced; they could not exist before the run. Everything that *constrains* the design —
§1, §4, §5, §7, §10, §12 — was written first and is what the module was built against. Two things the
code then made wrong in the first draft, corrected here rather than left standing:

- §4.5 originally proposed a **first**-difference (secant) probe for the strain-limit guard. Working
  through the arithmetic showed that a plateau-only law's secant ratio is `1 + E_pl·d/S₁`, which is
  greater than 1 and grows with the probe width — so the guard would have had a threshold that could
  be tuned to pass. The second difference is exactly `1.0` for the broken law and exactly `E_st/E_pl`
  for the bound one. The entry was rewritten before the guard was implemented.
- §1 originally listed `default_if_card`. Removed: the whole argument of §4 is that this card has no
  default, and a factory named in the ledger would have licensed one.

**What the mutation exercise found**, twice, which is the reason to run it rather than to describe it.

**First: §4.5's argument, reproduced as a measurement.** M1 removes the re-stiffening branch from the
shipped owner and killed **6** of 90 controls. Not one of them was a finite-difference, force-closure,
moment-closure, translation-invariance, purity or determinism control — every one of those still
passes, exactly, on a cage that limits nothing. All six kills came from
`assert_strain_limit_is_bound`, the tangent-is-state control, the kernel-versus-oracle control and the
path-dependence control. **84 of the 90 controls pass on a network whose entire registered mechanical
role is absent.** That is the number this entry exists to publish.

**Second: a mutant I designed to be conservative was not, and the design note was wrong.** M2 zeroes
the re-stiffening extent `c` while leaving the plateau extent `b` clipped at `ep`, which was intended
as a second pure missing-term case. It is not: past `ep` the *energy* goes flat while the *tension*
retains its plateau-end value, so `F = -grad E` fails and the gradient controls kill it. §4a is
recorded as measured rather than as designed, and this paragraph replaces the sentence that would
have claimed M2 was conservative. The genuinely conservative version of that defect is M1, which is
the shipped flag, which is why the flag is shipped.

## 14. Honest limits — what this entry does NOT establish

**No connector is wired. Endpoints are published; that is not the same thing.** This module supplies
the intermediate-filament side of `if_nucleus_linc`, `if_sf_plectin` and `if_cytosol_transfer`. Each
still needs the *other* owner's endpoint and a connector object, and `aleph/vertical/wiring.py` is
untouched by this entry. `wired_connector_count()` returns `0` and a control asserts it.

**`if_cytosol_transfer` is blocked on a decision that is not this lane's.** Its cytosol side is a
`TransferStencil`, which `docs/design/ENDPOINT_CONTRACT.md` §2 classifies as an
`ImmersedTransferEndpoint` scattering a *scalar*, not a force. The solid-to-fluid coupling protocol is
open decision 2 in `HANDOFF.md` §C. The endpoint published here is a force site with a quadrature
weight, which is what the IF *side* needs under any resolution of that decision — but whether it can
be joined to a transfer stencil is unresolved and is **not** asserted by anything here.

**No population role is withheld, and that is a claim about the registry, not about the physics.**
`PUBLISHABLE_ENDPOINT_ROLES` equals `IF_ENDPOINT_ROLES`: all three declared roles resolve to a material
point on a cable, which this representation has. `ROLE_POPULATIONS` is likewise **total** — every role
is available on both populations. That is honest and it is also weak: the registry gives no basis for
narrowing (it never says which population a LINC site or a plectin bond binds), and inventing one
would be a mechanism claim. The mechanism is shipped and exercised so that a *sourced* narrowing later
is a one-line data change rather than a redesign, but today it refuses nothing that occurs in practice.

**Turnover has no kinetics.** `propose_rebinding` is called by a caller this module does not have. No
rate constant, no load-dependent unbinding, no RNG stream, no `PROPOSE_KINETICS` phase. So `turnover_state`
records *that* rebindings happened and *what they did to the rest lengths*; it does not model *when*
they happen. Every path-dependence result in §8 is therefore a statement about the mechanism, not
about a rate.

**Physical magnitudes are not established.** Evidence rung `ANALYTIC_ORACLE`, quantitative status
`BLOCKED`. No modulus, knee strain, plateau width or crosslink stiffness here is a keratin or vimentin
number; they are chosen to make the algebra sharp, and both test cards carry `citation=None` so the
module's own refusal classifies them correctly. `citation_status` for the `intermediate_filament`
contract remains `UNSOURCED` and this entry does not improve it. **No strain-stiffening response at
any specific strain is claimed by anything in this port.**

**No nuclear deformation ceiling is demonstrated.** §4.5 establishes that this owner's law has a
strain limit. It does **not** establish that the limit sets a ceiling on nuclear deformation, because
`if_nucleus_linc` is not wired and there is no nucleus in any test here. The registered
`mechanical_role` is therefore *half* established: the strain-limiting property is measured, the
consequence for the nucleus is not.

**The keratin/vimentin mechanical difference is an input, not a result.** The two test cards differ
because I made them differ. The registered `vimentin_population` contract names this exact trap
("with one card assumed the difference is an input, not a result") and nothing here escapes it. What
*is* established is that the label is read: identical cables under different labels report different
tension, so the card is selected per population rather than shared.

**Not established: rest-length re-setting is the right turnover model.** A rebinding that resets rest
length to current length is *a* history mechanism with the property the contract requires. Whether it
is the mechanism keratin and vimentin turnover actually has is unknown here, and no literature was
read. It is marked `UNVERIFIED` rather than omitted.

**The reference implementation is unaudited by this lane.** Whether Aleph's re-derivation agrees with
the reference's cage is unknown, and this entry deliberately does not guess. `ALEPH-PORT-1802`'s audit
is cited, not extended.

**`ports/ledger/INDEX.md` was never touched by this lane, and that is worth recording because it
failed and then passed without this lane acting.** `INDEX.md` is outside this lane's file allocation
and six `ALEPH-PORT-32xx` entries were being added concurrently by six lanes.
`tests/ports/test_port_discipline.py::test_index_is_not_stale` failed for part of this session and
passes at the end of it, because another session regenerated the index; the file now lists
`ALEPH-PORT-3203`. Whoever takes the commit should regenerate it **last**, after every concurrent
entry has landed, or the same race will reopen. Reported rather than fixed, per `CLAUDE.md` §1.

**Not run:** any GPU job (zero), and the full suite. Only
`tests/vertical/test_intermediate_filament_controls.py`, `tests/state/` and `tests/ports/` were
exercised — other lanes are writing concurrently, and `CLAUDE.md` §1 says to report a foreign breakage
rather than fix it, which means not conflating one with mine.
