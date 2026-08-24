# ALEPH-PORT-3205 — lamellipodium: a locally adaptive branched actin network

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3205` |
| Lane | `L32 owner modules — lamellipodium` |
| Status | `PROPOSED` |
| Written | `2026-07-30`, **before the code**, per PLAN §0.2.5 |
| Port class | `RE-DERIVED` — nothing was ported, nothing was read |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

Aleph target file: `aleph/vertical/lamellipodium.py`.
Controls: `tests/vertical/test_lamellipodium_controls.py`.

```python
from aleph.vertical.lamellipodium import (
    LAMELLIPODIUM_ENDPOINT_ROLES,
    OWNED_STATE_KEYS,
    PUBLISHABLE_ENDPOINT_ROLES,
    RATE_PRIORS_EVIDENCE_OWED,
    ROLE_SITE_KINDS,
    BarbedEndState,
    BranchJunction,
    BranchProposal,
    LamellipodiumAttachmentSite,
    LamellipodiumCard,
    LamellipodiumEndpointRole,
    LamellipodiumEvidenceError,
    LamellipodiumNetwork,
    LamellipodiumOwnershipError,
    LamellipodiumRoleError,
    LamellipodiumStaleHandleError,
    LamellipodiumTopologyError,
    MaterialPoint,
    MaterialPointSink,
    ProtrusionFilament,
    SiteKind,
    assert_branch_angle_is_bound,
    assert_no_shared_node_storage_with,
    assert_population_count_is_an_output,
    assert_state_keys_disjoint_from,
    axial_energy_and_forces,
    bending_energy_and_forces,
    branch_angle_energy_and_forces,
    default_lamellipodium_card,
    junction_link_energy_and_forces,
    owned_state_keys,
    rate_prior,
    straight_filament_nodes,
)
```

**This is a state owner, and it wires zero connectors.** It publishes endpoints for the five declared
connectors that name `lamellipodium`; it implements none of them. A `Binding` row in
`aleph/vertical/wiring.py` belongs to the commit that wires a connector, never to this one, and this
lane added none.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — **recorded, not verified by this lane** |
| Source path | `ffn_sim/ac/engine/protrusion.py`, `ffn_sim/ac/weave/lamellipodium.py`, `ffn_sim/ac/weave/branch_angle.py`, `ffn_sim/ac/weave/membrane_ratchet.py` — **named, not read** |
| Source symbol(s) | **none** |
| Read from | **neither the commit nor the working tree.** No file under `/Users/sw1/ffn_cellsim` was opened by this lane. |
| Working tree == commit? | not applicable — nothing was read, so there is no revision to have mis-cited |

The commit hash and the four paths are transcribed from `ALEPH-PORT-1701` §2, which is an Aleph
document. This entry asserts only that `ALEPH-PORT-1701` records them; it does **not** assert that the
hash resolves today, because checking would have meant touching a tree this lane was told not to open.
That is marked `UNVERIFIED` here rather than omitted.

Everything below was derived from three Aleph documents:
`aleph/state/census_actomyosin.py` (the registered `lamellipodium` contract, line 261),
`aleph/state/connectors_protrusion_motor.py` and `aleph/state/connectors_crosssystem_cytosol.py`
(the five connectors that name this component), and `docs/design/ENDPOINT_CONTRACT.md`.

## 3. Why source-derived porting beats clean-room

**It does not.** `ALEPH-PORT-1701` already audited the reference's lamellipodium and returned the
verdict RE-DERIVE, on three findings this lane inherits rather than re-establishes:

1. The reference's lamellipodium owner and its mechanics delegate are **dead code outside that
   repository's own tests**. There is nothing to inherit but a class shape.
2. Its branch-angle kernel is real numeric code — but its **branch-angle citation is recorded
   `DOI_DEAD`** by that repository's own audit. Porting a formula whose only parameter has a dead
   citation inherits the unverifiability along with the algebra. This is the single strongest reason
   the angular term here ships **with no angle value at all**.
3. Its branch-angle harmonic is genuinely native — **on the cortex path, under a cortex owner, never
   under a lamellipodium owner.**

Finding 3 is not a scheduling detail; it is the defect this module exists not to have, and §4a returns
to it. A branch angle evaluated inside another owner's arrays means the branched population is a slice
of the cortical one. The registered contract forbids exactly that, in its own words: *"spatial overlap,
co-location in a shared array, or a shared index space creates no mechanical connection."* A branched
network that lives in the cortex's array is not a lamellipodium with a wrong number in it — it is a
lamellipodium that does not exist as a separate mechanical body, and every force it reports balances.

An axial spring with a rest-length normalisation, a tangent-difference bending stencil, a barycentric
link and a `(cos θ − cos θ₀)²` angular bond are laws any competent author re-derives in an afternoon.
`ports/TEMPLATE.md` §3 names that as a reason to write clean-room, not a reason to port.

## 4. Physical or mathematical law represented

Four energy terms, each with its exact analytic gradient, plus a population mechanism that carries no
rate. Units µm / pN / pN·µm / s throughout.

**Axial.** `E = Σ_s (k_s / 2) (L_s − L0_s)² / L0_s`, giving segment tension `T = k (L − L0) / L0`,
which pulls the two ends together. Dividing by the rest length makes the modulus
resolution-independent: cut a filament into twice as many segments and you get twice as many springs
each twice as stiff, i.e. the same end-to-end compliance. Two-sided, deliberately — a compressed
lamellipodial segment must push, because that is the only channel through which added monomer becomes
load (see the elongation law below).

**Bending.** On unit tangents `t1`, `t2` of the two segments meeting at interior node `i`:
`E_i = (κ_i / 2h_i) |t2 − t1|²` with `h_i` the **rest** Voronoi length `(L0_prev + L0_next)/2`. Using
the rest length decouples bending from stretching: a uniformly stretched or contracted straight
filament has exactly zero bending energy and exactly zero bending force. That is a property and not an
oversight — it is precisely why the axial law cannot be optional. Small-angle limit
`|t2 − t1| = 2 sin(θ/2) → θ → C h` gives `E_i → (κ/2) C² h`, so the sum tends to `∫ (κ/2) C² ds` with
no coefficient left over.

**Branch junction link.** A daughter filament's pointed-end node is bonded to a material point on its
mother: `E = Σ_j (k_j / 2)(d_j − r_j)²` between the daughter root node and the interpolated mother
point. The gradient is scattered to the three participating nodes with the same weights that built the
mother point, which is the condition under which the chain rule closes.

**Branch angular bond — the term that makes this a branched network rather than a bundle.** With
`u` the mother's unit tangent through the junction segment (oriented pointed → barbed) and `v` the
daughter's first-segment unit tangent (same orientation), and `c = u · v`:

```
E = Σ_j (k_θ,j / 2) (c_j − c0_j)²
```

Written in the **cosine** rather than in the angle. Three consequences, each load-bearing:

1. `d(arccos)/dc = −1/√(1−c²)` diverges at `c = ±1`. An angle-form bond would hand a straightening or
   a folding branch an unbounded force at exactly the configuration the missing-term control drives it
   into, which would make the control measure the singularity instead of the law.
2. The energy depends on `u` and `v` only through their inner product, so it is **degenerate in the
   azimuth** around the mother axis: the daughter is free to sit anywhere on the cone of half-angle
   `θ₀`. That is the right physics for a nucleator that sets a branch angle and not a branch plane,
   and it is asserted rather than assumed.
3. `c0` is a **datum of the junction, supplied by the caller**, not a constant of this module. See
   §4b.

Gradient. With `u_vec = p[m_b] − p[m_a]`, `v_vec = p[d_1] − p[d_0]`, `û = u_vec/|u_vec|`,
`v̂ = v_vec/|v_vec|`:

```
∂c/∂p[m_b] = (I − û ûᵀ) v̂ / |u_vec| = −∂c/∂p[m_a]
∂c/∂p[d_1] = (I − v̂ v̂ᵀ) û / |v_vec| = −∂c/∂p[d_0]
```

and `F = −k_θ (c − c0) ∂c/∂p`. Two exact identities follow and both are asserted:

* **Force closure** is immediate, because each pair is equal and opposite.
* **Moment closure** is not, and it is worth the two lines. The mother pair contributes
  `u_vec × (I − û ûᵀ)v̂ / |u_vec| = û × (v̂ − c û) = û × v̂`, and the daughter pair contributes
  `v̂ × (û − c v̂) = v̂ × û`. They cancel identically, so the angular bond injects no couple. A term
  that closed in force but not in moment would spin the protrusion for free.

**Elongation, and why it is a mechanical statement.** A barbed-end elongation of `Δ` **adds `Δ` to the
terminal segment's rest length and moves no node.** The tip therefore finds itself inside a segment
whose stress-free length just grew, the axial law reads `L − L0 < 0`, and the segment pushes the tip
outward. That is the whole content of "polymerisation does work against a load" expressed in an energy
rather than in a prescribed velocity — and it is the half of the Brownian ratchet this owner is
allowed to have. The other half, *how much* `Δ` arrives per unit time under a given load, is the
`lamellipodium_membrane_contact` connector's, because it needs the free elongation rate (§4b). The
alternative implementation — move the tip node by `Δ` and leave the rest length alone — produces the
same kinematics at exactly zero force, which is a protrusion that grows without pushing.

**Material coordinates are measured from the pointed end.** Rest arclength from the end that
polymerisation does not move. Measuring from the low node index instead would make every published
attachment on a growing filament slide backwards by `Δ` each time a monomer is added — an attachment
that moves when the filament grows is an attachment that has silently detached. The barbed end itself
is therefore addressed **topologically** (by filament id, "the tip, wherever it now is") and never by
a material coordinate, because the tip is the one place on a growing filament where a material address
is wrong.

### 4a. The ownership law, which is the design

Let `F` be the physical filaments and `O` the components. `ALEPH-PORT-1701` §4 derives the rule: the
ownership assignment `o : F → O` is a **function**, so the populations partition `F`, and disjointness
is not an extra axiom but what "function" means. What goes wrong is that `o` is never constructed and
a *relation* is used instead — one array with a region label.

This module makes the failure checkable at two levels rather than asserting it in prose:

* **On names, before allocation.** `assert_state_keys_disjoint_from` refuses any overlap between this
  owner's nine registered keys and another owner's. This is the only point at which the check is
  cheap: once two owners share a buffer, every force-closure test still passes, because the two halves
  of the Newton pair are then the same numbers.
* **On storage, after allocation.** `assert_no_shared_node_storage_with` refuses when another owner's
  position block shares memory with this one (`np.shares_memory`), and **passes** when it merely holds
  identical values at identical coordinates. That asymmetry is the contract's sentence made
  executable: co-location is not a connection, a view is.

### 4b. No rate value is carried, and a consumer that needs one is refused

The registered contract owes evidence on four priors — the branch angle and its spread, the
branch-nucleation rate, the capping rate, and the free barbed-end elongation rate — and says *"No
value for any of the four is carried here."* This module carries none of them either, and the
refusal is executable: `rate_prior(name)` raises `LamellipodiumEvidenceError` for all five names in
`RATE_PRIORS_EVIDENCE_OWED`, naming the prior and what would have to be established. There is no
default, no fallback, and no "reasonable value" anywhere in the module.

What ships instead is the **mechanism**: `propose_branch`, `propose_capping` and `propose_elongation`
take their magnitudes from the caller. A caller that cannot defend a branch angle cannot obtain one
here.

**And there is no target population count anywhere.** The contract's rule is that the count is a
dynamic output — *"a population fixed to the observed population cannot be tested against the observed
population, so the one comparison that could falsify the mechanism is the one the fixing removes."*
So: no constructor argument, no card field, no method and no module constant names a target,
desired, or expected filament count; nothing resamples, rejects or rescales a population toward one;
and `assert_population_count_is_an_output` is shipped so the absence is callable from assembly code
rather than believed. A control scans the module text for the whole family of names and fails on any
of them, because "we did not add one" is not a property a reader can check and a scan is.

## 4c. Mutation testing of the controls

Three mutants applied to `aleph/vertical/lamellipodium.py`, the controls run against each, then the
source restored and verified byte-identical. Recorded here because a control suite no mutant can break
is decoration.

| # | Mutation | Tests killed | Survived? |
|---|---|---|---|
| M1 | `branch_angle_energy_and_forces`: scatter the daughter-side gradient with the **mother's** projector (`I − û ûᵀ` on both sides) | **4** | no |
| M2 | `axial_energy_and_forces`: drop the `/ rest` normalisation (`k Δ²` for `k Δ²/L0`, in the energy and the tension together, so `F = −grad E` stays exactly true) | **3** | no |
| M3 | `_apply_elongation`: move the tip node by `Δ` **as well as** extending the rest length, so growth generates no force | **4** | no |

The controls each mutant killed, and the discussion of what they taught, are in §8b.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| node position | µm | 1e-6 m | finite |
| `axial_modulus_pn` | pN | 1e-12 N | `> 0` (a stiffness × a length) |
| `bending_rigidity_pn_um2` | pN·µm² | 1e-24 N·m² | `>= 0` |
| `junction_link_stiffness_pn_per_um` | pN/µm | 1e-6 N/m | `>= 0` |
| `branch_angular_stiffness_pn_um` | pN·µm | 1e-18 J | `>= 0` |
| `cos_rest_angle` | 1 | 1 | `[-1, 1]`, **supplied by the caller, never defaulted** |
| elongation increment `Δ` | µm | 1e-6 m | `> 0` |
| filament count | count | 1 | `>= 0`, **an output** |
| energy | pN·µm | 1e-18 J | `>= 0` |

Singular and boundary cases, each with the behaviour this module requires:

- **Zero filaments.** A legal state, not an error: a lamellipodium that has not nucleated yet owns an
  empty population and reports energy `0.0` and no forces. Nothing here divides by the count.
- **Collapsed segment** (`L = 0`): refused. The direction is genuinely undefined.
- **Zero-length junction link**: refused. A silent epsilon returns a finite force whose direction came
  from the round-off in the last subtraction.
- **A branch whose mother is itself**: refused. A filament rooted on its own side is a kink
  mis-declared as a topology element, and its angular bond would be a function of one tangent.
- **A material coordinate outside its filament**: refused, never clamped. Clamping an attachment to
  the end of a filament moves a load path without telling anybody.
- **A rate prior requested**: refused, never defaulted (§4b).
- **A target population count requested**: there is nothing to request; see §4b.
- **Capping an already-capped filament**, or **elongating a capped one**: refused. A cap that can be
  applied twice, or that a growth event can ignore, is not a state.
- **A non-positive elongation increment**: refused. Depolymerisation is a different mechanism and this
  owner does not declare one.
- **An endpoint handle used after a committed population change**: refused with
  `LamellipodiumStaleHandleError`.

Invariants, each with the test that asserts it:

- **I1.** `F = −grad E` per term, observed order > 1.6 —
  `tests/vertical/test_lamellipodium_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order`
- **I2.** Force and moment closure per term, against the constituent scale —
  `::TestClosureAndInvariance::test_every_term_closes_on_its_own`
- **I3.** A handle goes stale when the population moves —
  `::TestOwnershipAndTheAcceptedStep::test_a_handle_held_across_a_committed_nucleation_is_refused_and_a_fresh_one_works`
- **I4.** Rollback restores the network bit-identically, including the array **shape** —
  `::TestOwnershipAndTheAcceptedStep::test_rollback_restores_the_network_bit_identically`
- **I5.** An accepted elongation moves no published attachment —
  `::TestOwnershipAndTheAcceptedStep::test_an_elongation_moves_no_existing_attachment`
- **I6.** Declared state keys equal the registered contract, and collide with no other owner —
  `::TestOwnershipAndTheAcceptedStep::test_the_declared_keys_are_exactly_the_registered_ones`,
  `::test_the_declared_keys_collide_with_no_other_registered_owner`
- **I7.** No name in the module carries a target population count —
  `::TestTheDeclarationsAgreeWithTheCode::test_the_module_declares_no_target_population_count`
- **I8.** The role enum is exactly the registered connector set —
  `::TestTheDeclarationsAgreeWithTheCode::test_the_role_enum_matches_the_registered_connectors`

### The endpoint contract, and the one place this owner is stricter than it

`docs/design/ENDPOINT_CONTRACT.md` I3 requires a handle to go stale when the topology moves, and calls
the guard "deliberately conservative" because adding a crosslink does not move a node. **For this
owner it is not conservative, it is load-bearing**, and that is a difference worth recording rather
than inheriting silently: a committed nucleation **reallocates the node arrays**, because the
population size is a state variable. A sink issued before the commit holds a reference to the previous
`forces` buffer, which nothing integrates any more. Its scatter would not be wrong by a weight — it
would be discarded in silence. `::test_a_handle_held_across_a_committed_nucleation_is_refused_and_a_fresh_one_works`
asserts both halves: the stale handle raises, and the array identity really did change, so the hazard
is demonstrated rather than postulated.

## 6. Source evidence class and known retractions

No claim is made about the reference implementation's lamellipodium by this lane, because none of it
was read. `ALEPH-PORT-1701` §6 is the audit of record and it reports, in that repository's own
documents: the compartment rated `KERNEL_BOUND` with zero device execution and "declared-only in
practice"; the owner and its mechanics delegate unreachable outside tests; and **the branch-angle
citation recorded `DOI_DEAD`** — a retraction in substance, since the parameter's provenance does not
resolve.

The registered `lamellipodium` contract carries `citation_status: UNSOURCED`, and nothing here changes
that. `default_lamellipodium_card()` returns four unsourced moduli and says so; the four evidence-owed
rate priors have no value here at all, sourced or otherwise.

Evidence rung `ANALYTIC_ORACLE`; quantitative status `BLOCKED`. Nothing in this entry is evidence that
a protrusion velocity, a branch angle, a network stiffness or a nucleation rate has the value a
lamellipodium has.

## 7. Independent oracle or derivation

Four, none of which is the reference implementation:

1. **Central finite differences** of the module's own energies, per term. Independent of every
   constant in the law, and blind by construction to a missing term.
2. **The rest-length-normalised closed form.** A chain in equilibrium under an end load `F` extends by
   `ΔL = F L0 / k` regardless of how many segments it was cut into. This catches what the finite
   difference cannot: dropping the normalisation leaves `F = −grad E` exactly true while making the
   modulus resolution-dependent (mutant M2).
3. **The angular bond's two closure identities**, derived in §4 from the projector algebra. Moment
   closure in particular is an identity of the gradient's structure, not of its magnitude, so a
   wrongly-scattered angular gradient fails it while still summing to zero in force (mutant M1).
4. **The azimuthal degeneracy.** `E` depends on `u` and `v` only through `u · v`, so rotating the
   daughter around the mother's axis at fixed opening angle must change the energy by exactly zero.
   A bond that secretly picked a branch plane fails this and passes every gradient check.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_lamellipodium_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order` | Per term: `energy > 0` first, then max&#124;F + ∇E&#124; / max&#124;F&#124; at steps 4e-4 / 2e-4 / 1e-4 µm with observed order > 1.6 |
| Positive | `tests/vertical/test_lamellipodium_controls.py::TestTheBranchAngleLawIsActuallyBound::test_folding_a_daughter_onto_its_mother_costs_energy` | Rotating a daughter onto its mother's axis raises the network energy, and the rise comes from the angular term alone |
| Positive | `tests/vertical/test_lamellipodium_controls.py::TestThePopulationCountIsAnOutput::test_nucleation_and_capping_change_the_count_during_a_run` | The filament count and the node count both change across accepted steps, in both directions for the free-barbed-end population, with no target anywhere |

### 8a. Measured, this session, on this tree

Every number below was produced by a run in this session with
`/Users/sw1/miniconda3/envs/aleph/bin/python`; none is quoted from memory.

Finite-difference agreement on the four-term fixture, steps 4e-4 / 2e-4 / 1e-4 µm:

| Term | Energy (pN·µm) | Force scale (pN) | max&#124;F+∇E&#124; at 1e-4 (pN) | relative | observed order |
|---|---|---|---|---|---|
| axial | 4.330838 | 96.13 | 8.192e-06 | 8.522e-08 | **2.000** |
| bending | 1.983229 | 39.78 | 5.145e-06 | 1.293e-07 | **2.000** |
| junction link | 0.011396 | 1.263 | 1.876e-06 | 1.486e-06 | **2.000** |
| branch angle | 0.204459 | 5.088 | 2.051e-06 | 4.031e-07 | **2.000** |
| total | 6.529923 | 95.77 | 4.778e-06 | 4.989e-08 | — |

**The junction link's relative column is the largest and its absolute column is the smallest**, which
is worth naming because it is the reason its tolerance is ten times looser than the others'. All four
terms sit on one truncation floor of a few times `1e-6` pN at `h = 1e-4` µm; what differs by a factor
of 76 is the force each carries in this fixture, and a relative tolerance divides by that. The
discriminating quantity is the order, which is **2.000** for all four — the first draft of this
fixture instead nucleated the branches *before* straining the network, which stretched every junction
link to several times its rest length and pushed the junction term's relative residual to `9.2e-06`.
That was a bad configuration rather than a bad gradient, and it was fixed in the fixture rather than
absorbed into a tolerance.

Closure on the isolated network (resultant ÷ constituent scale, `Σ|f_i|` and `Σ|r_i||f_i|`):

| Term | force | moment |
|---|---|---|
| axial | 8.90e-18 | 1.21e-17 |
| bending | 2.01e-17 | 1.93e-17 |
| junction link | 6.00e-18 | 5.01e-17 |
| branch angle | **0.00e+00** | 9.05e-18 |
| total | 9.04e-18 | 1.05e-17 |

The angular bond's force residual is an exact zero rather than a small number, which is what the
projector form predicts: its four gradients are two equal-and-opposite pairs added with the same
`scale` array, so the cancellation happens in the same floating-point value rather than across a sum.

Suite: **94 passed** in `tests/vertical/test_lamellipodium_controls.py`, 0.96 s.

### 8b. What the mutants taught, with the numbers

Baseline before mutation: 0 failed, 94 passed.

| # | Killed | Which controls |
|---|---|---|
| M1 | **4** | `TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order[branch_angle_term]`, `::test_the_total_internal_force_is_the_gradient_of_the_total_energy`, `TestClosureAndInvariance::test_every_term_closes_on_its_own[branch_angle_term]`, `::test_the_internal_moments_sum_to_zero` |
| M2 | **3** | `TestGrowthGeneratesForce::test_the_energy_stored_by_an_elongation_matches_its_closed_form`, `TestClosureAndInvariance::test_end_to_end_axial_stiffness_does_not_depend_on_node_count`, `::test_translating_the_whole_network_changes_no_energy` |
| M3 | **4** | `TestGrowthGeneratesForce::test_an_accepted_elongation_pushes_the_tip_outward`, `::test_an_elongation_extends_the_rest_length_and_moves_no_node`, `::test_the_energy_stored_by_an_elongation_matches_its_closed_form`, `TestTheNegativeControlsFail::test_growth_generates_no_force_when_the_tip_is_moved_instead` |

**M1 is the case for checking moment closure separately from force closure.** Scattering the daughter
side with the mother's projector leaves the two pairs equal-and-opposite, so force closure still holds
for that term and `test_every_term_closes_on_its_own` catches it only through its moment half — the
`û × v̂ + v̂ × û` cancellation derived in §4 stops cancelling. Force closure alone would have passed
this mutant, and so would every energy check, since the energy is untouched.

**M2 is the case that a gradient control and an oracle are not substitutes.** Dropping the
normalisation from the energy and the tension together leaves `F = −grad E` **exactly true**, so
neither finite-difference control moved. Its kills came from the two independent closed forms. One of
the three is weaker than it looks and is reported as such: `test_translating_the_whole_network_changes_no_energy`
fails under M2 on its **force** comparison at `rtol=1e-12`, because the mutant changes the relative
magnitudes of the terms and so the cancellation round-off after a 4 µm translation, not because the
mutant broke translation invariance. It is a real failure and a weak signal, and calling it the
latter is the difference between three kills and two.

**M3 is the same lesson on the mechanism side.** Moving the tip *and* extending the rest length is
kinematically identical to real growth — same length, same rate, same direction — and produces exactly
zero force, because the configuration it produces is a perfectly valid stress-free one. No gradient,
closure or conservation control can object to it. Only a control that demands growth *push* can, and
three of its four kills are exactly those.

Source restored and verified **byte-identical** after each mutant: the SHA-256 of the restored file
equals that of the pre-mutation copy (`264cca1794a8f5fa…`), and `git diff` on
`aleph/vertical/lamellipodium.py` is empty.

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_lamellipodium_controls.py::TestTheNegativeControlsFail::test_a_flipped_axial_sign_is_caught_by_the_gradient_test` | With `flip_axial_sign=True` the FD disagreement rises above 0.5 relative, while the energy is bit-unchanged |
| Negative (must fail) | `tests/vertical/test_lamellipodium_controls.py::TestTheNegativeControlsFail::test_folding_a_branch_flat_is_free_when_the_angular_bond_is_unbound` | With `omit_branch_angle_law=True` folding a daughter onto its mother costs **exactly** `0.0`, and the shipped guard raises on it; with the law bound the same fold costs energy |
| Negative (must fail) | `tests/vertical/test_lamellipodium_controls.py::TestTheNegativeControlsFail::test_a_capped_barbed_end_still_grows_when_the_cap_is_ignored` | With `allow_capped_barbed_end_growth=True` a capped filament elongates; with the cap honoured the same call raises and the rest length changes by **exactly** `0.0` |

All three breaks are flags on the **shipped** owner, defaulting to `False`, so each control drives the
real code path rather than a hand-edited copy — a negative control against a duplicate proves only
that the duplicate is broken. Each is asserted in **both** directions: asserting only that the broken
branch misbehaves passes for a module that misbehaves always, and asserting only that the healthy
branch is quiet passes for a module that does nothing at all.

## 10. Numerical and precision envelope

Working and accumulation precision: float64 throughout; no reduced-precision path exists and no
backend kernel is invoked by the energy terms.

Finite-difference steps 4e-4 / 2e-4 / 1e-4 µm on coordinates of order 0.1–3 µm. The round-off floor
for a central difference is `eps^(1/3) ≈ 6e-6` **relative**; the smallest step used is more than an
order of magnitude above it. Below that floor, cancellation in the energy evaluation dominates the
truncation error, the apparent convergence order goes negative, and a correct gradient becomes
indistinguishable from a wrong one — the afternoon `PLAN.md` §6.1 records.

Tolerances and why they are not looser: per-term FD agreement is asserted at `1e-6` relative for three
terms against measured worst cases of `8.5e-08`, `1.3e-07` and `4.0e-07`, and at `1e-5` for the
junction link against a measured `1.5e-06`. That one entry is looser for the arithmetic reason given in
§8a — its absolute residual is the *smallest* of the four and its force scale is 76 times smaller than
the axial term's — and not because its gradient is worse. The observed order is required to be above
1.6 against a measured 2.000 for all four; that is what actually separates a wrong gradient from a
stiff one, and it is tight enough that M1, which perturbs the angular gradient at the `1e-1` level,
cannot pass.
Closure is asserted at `1e-12` of the **constituent** scale against a measured `5e-17`. The scale is
`Σ|f_i|` and `Σ|r_i||f_i|`, never the resultant: a resultant-based tolerance gets *stricter the more
correct the physics is*, because a correct Newton pair cancels, and `PLAN.md` records that defect
rejecting forty thousand consecutive steps while reporting a plausible tension.

The angular bond has no small-denominator regime by construction: the projector form divides by
`|u_vec|` and `|v_vec|`, both of which are refused at zero, and the cosine form has no `1/√(1−c²)`.

Outside the envelope: every singular case in §5 refuses with a typed exception rather than degrading
silently. There is no input range over which this module returns a plausible number it cannot stand
behind.

## 11. Production-backend residency and transfer

Host, numpy, float64. The energy terms are `O(S)`, `O(T)` and `O(J)` gather/scatter kernels with no
host round-trip per step, so a device port is mechanical when one is wanted. The owner touches
`ctx.backend` only in `accumulate`, through `scale` / `add` / `max_abs`, exactly as the other owners
do. **No GPU work of any kind was run by this lane**, on either host, and no authorization was sought
or held.

One residency consequence of the adaptive population is worth naming before somebody meets it on a
device: a committed nucleation **reallocates** the node arrays. On host that is a `vstack`; on a device
it is a new allocation and a copy, once per accepted topology change rather than once per step. Any
future device path must not cache a device pointer across a commit, for the same reason no connector
may cache a sink across one.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered this module, because no reference file was read. Two
controls assert this rather than leaving it as a claim —
`::TestTheDeclarationsAgreeWithTheCode::test_no_reference_project_vocabulary_survives` scans the module
text for the provider's identifiers, and `::test_the_module_imports_nothing_from_validation` asserts
the `aleph/** → validation/**` firewall.

The prose that *is* here was written for this module. Where it quotes the registered contract or a
connector contract it says so and quotes exactly, because the paraphrase of an ownership rule is the
place ownership rules go to die.

**One style deviation, reported rather than hidden.** `ruff` reports 79 `E501` lines across the module
and its controls — every one of them a prose line of 101–106 columns against a configured limit of
100, and none of them code. It is not repaired, and the reason is a measurement rather than a
preference: two attempts at an automated re-flow were made and both corrupted text. The first moved a
word out of an f-string fragment into a plain one, so an interpolation stopped interpolating; the
second dropped the trailing space between two concatenated fragments, so a refusal message read
`a shared statekey`. Neither was caught by any control, because no control asserts on the wording of a
refusal. The files were restored to the measured version both times. The module is left with the
long prose lines and this note, which is a smaller defect than a silently mangled message, and it is
the same deviation `cortex_filaments.py`, `ecm.py`, `sf_arc.py` and their controls already carry (78
`E501` across those four). `ruff` is not run by any test in this repository.
One `UP035` is also left standing — `from typing import Iterable, Sequence` rather than
`collections.abc` — because that is verbatim what the three existing owner modules do, and a single
module diverging on an import convention is worse than the convention being old.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 94 controls pass (`tests/vertical/test_lamellipodium_controls.py`, 2026-07-30, 0.96 s). Three mutants applied, killing 4 / 3 / 4 controls; none survived; source restored byte-identically. Status stays `PROPOSED`: the entry is agent-written, no human has reviewed it, and the physical magnitudes it uses are unsourced by construction. |
| Reviewer | Agent-proposed (lane L32). Unratified. No PI review. |
| Rollback | Delete `aleph/vertical/lamellipodium.py` and `tests/vertical/test_lamellipodium_controls.py`. Nothing breaks: no module in the tree imports either, `aleph/vertical/__init__.py` was not touched, and no connector reaches them. |

## 14. Honest limits — what this entry does NOT establish

**Wired connector count contributed by this module: zero.** Five declared connectors name
`lamellipodium` — `lamellipodium_membrane_contact`, `lamellipodium_cortex_seam`,
`lamellipodium_nascent_fa`, `nmii_lamellipodium_motor`, `lamellipodium_cytosol_transfer`. This module
publishes an endpoint for each; it implements none, and `aleph/vertical/wiring.py` was not touched.
The count of wired connectors is unchanged by this work and any claim otherwise should be checked
against that file, which is the only place the number lives.

**The Brownian ratchet is not implemented, and half of it cannot be.** What this owner supplies is the
mechanical half: an accepted elongation increment becomes an outward force through the axial law
(§4). The rate half — how large an increment arrives per unit time under a given load — requires the
free barbed-end elongation rate, which is one of the four priors evidence is owed on. So the
`lamellipodium_membrane_contact` endpoint is publishable and **inert**: a connector can attach to it
today and cannot make it grow. That is the honest state and it is why `rate_prior` refuses rather than
returning a placeholder. **No protrusion velocity of any kind is established here.**

**The barbed-end monomer sink publishes a demand and not a chemistry.** The
`lamellipodium_cytosol_transfer` contract carries two couplings and says why the second is not
optional: without it "the lamellipodium would polymerise from an infinite reservoir". Both roles are
reachable from this side — the drag quadrature points and the barbed-end sinks are published as two
site kinds of one role — but what the sink reports is **the elongation actually accepted, in µm**, not
a monomer count and not a demand rate. Converting µm of F-actin into monomers needs a monomers-per-µm
constant that nobody here has sourced, and inventing one would be the same failure as inventing a rate.
So the depletion applied to the G-actin field is the connector's arithmetic, and **this module does not
establish that the monomer budget closes**, because it never sees the field.

**Nothing here has been stepped inside a world.** The owner implements the participant protocol and
is not in `aleph/vertical/assembly.py`. It has controls; it has never been assembled, never been
relaxed, and never met another owner. Every closure number in §8a is a statement about an isolated
network with no connector attached.

**The population mechanism has no driver.** `propose_branch`, `propose_capping` and
`propose_elongation` are the mechanism, and nothing in Aleph calls them except the controls. The
adaptivity is therefore demonstrated (the count really does change across accepted steps, and the
arrays really are reallocated) and **not** exercised at any scale: the largest population any control
builds is a handful of filaments over a handful of steps. No statement about the statistics of a
branched network — branch density, filament length distribution, the nucleation/capping balance —
follows from anything here, and all three of those would need the priors §4b refuses.

**Debranching, severing and depolymerisation are absent rather than zero.** The registered
representation names nucleation, elongation and capping; it does not name a filament leaving the
mechanical population. So the node count here is monotone non-decreasing within a run, capping removes
a filament from the *growing* population and not from the *mechanical* one, and a run long enough for
turnover to matter would exhaust memory before it exhausted the model. That is a limit of the declared
representation and this module does not paper over it.

**No steric exclusion.** The filaments of this network can pass through each other and through
everything else. The cortex owner carries a WCA core; this one does not, because the registered
`lamellipodium` representation names segments, bending, junctions with an angular bond, and barbed
ends, and does not name excluded volume. Adding one would be a mechanism the census does not declare —
but the consequence is real and is stated rather than left to be discovered: a dense branched network
here has no packing limit.

**Physical magnitudes are not established.** Every number on `default_lamellipodium_card()` is an
unsourced placeholder of roughly plausible order, and every number in the controls was chosen to make
the algebra sharp. Evidence rung `ANALYTIC_ORACLE`, quantitative status `BLOCKED`.

**The reference implementation was not read by this lane**, so whether this re-derivation agrees or
disagrees with it is unknown here and is deliberately not guessed. `ALEPH-PORT-1701` is the audit of
record and this entry adds nothing to it.

**Not run:** any GPU job (zero, and no authorization sought or held), and the full suite.
`tests/vertical/`, `tests/state/` and `tests/ports/` were exercised. Five sibling lanes were writing
the other owner modules concurrently and `CLAUDE.md` §1 says to report a foreign breakage rather than
fix it, which means not conflating one with mine:

* `tests/ports/test_port_discipline.py::test_index_is_not_stale` fails, and **this entry is now one of
  its causes.** `ports/ledger/INDEX.md` is generated, and adding any entry makes it stale. It is a
  shared file this lane was told not to touch, so it was not regenerated; the run that lands these
  entries should regenerate it in the same commit. The test was already failing before this entry
  existed, on the untracked `ALEPH-PORT-3101`.
* `tests/vertical/test_medium_controls.py` fails twice, in another lane's file, mid-write. Reported,
  not touched.

`tests/vertical/` was 999 passed / 2 failed on that basis; this entry's own file was 94 passed, 0
failed.
