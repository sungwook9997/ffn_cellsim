# ALEPH-PORT-3204 — filopodium: a slender bundled strut whose root is a connector and whose effective stiffness refuses

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3204` |
| Lane | owner modules — `filopodium` state owner |
| Status | `PROPOSED` |
| Written | `2026-07-30`, **before the code**, per `PLAN.md` §0.2.5. §4a and §8 measured `2026-07-31` — see §14 |
| Port class | `RE-DERIVED` — nothing was ported; no reference file was opened by this lane |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.vertical.filopodium import (
    FILOPODIUM_ENDPOINT_ROLES,
    OWNED_STATE_KEYS,
    PUBLISHABLE_ENDPOINT_ROLES,
    ROLE_ALLOWED_ZONES,
    BundleCard,
    BundleCrosslink,
    BundleFilament,
    BundleZone,
    FilopodiumAttachmentRole,
    FilopodiumAttachmentSite,
    FilopodiumBundle,
    FilopodiumOwnershipError,
    FilopodiumPopulationError,
    FilopodiumRoleError,
    FilopodiumStaleHandleError,
    FilopodiumUnsourcedParameterError,
    MaterialPoint,
    MaterialPointSink,
    MonomerSinkEndpoint,
    RootAttachmentState,
    TipState,
    TransferCoupling,
    assert_bundle_cannot_contract_for_free,
    assert_monomer_sink_is_declared_with_drag,
    assert_population_count_is_an_output,
    assert_state_keys_disjoint_from,
    axial_energy_and_forces,
    bending_energy_and_forces,
    crosslink_spacing_energy_and_forces,
    owned_state_keys,
    straight_bundle_nodes,
)
```

Aleph target file: `aleph/vertical/filopodium.py`. Controls:
`tests/vertical/test_filopodium_controls.py`.

**Deliberately absent, and the absence is the design:** there is no `default_bundle_card()`, no
`effective_bending_rigidity` value, no buckling load, no critical length, no tip elongation rate, no
root attachment rate, and no parameter anywhere that names a target filopodium count. Every one of
those is a *method that raises* or a *name that does not exist*, and §9 names the controls that hold
each of them to it.

**This module wires zero connectors.** It publishes endpoints. The wiring module is not touched by
this entry and no `Binding` row is added — a `Binding` is added in the same commit that
wires a connector, never before.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — **quoted from `ALEPH-PORT-1703` §2, an Aleph document. Not independently verified by this lane.** |
| Source path | `ffn_sim/ac/engine/protrusion.py`, `ffn_sim/ac/weave/regions.py` — **named, not read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened by this lane, and no `git` command was run against it. |
| Working tree == commit? | not applicable — nothing was read, so there is no revision to have mis-cited |

The paths are recorded so a later auditor can establish independently whether Aleph's re-derivation
agrees or disagrees with the reference. Listing them is **not** a claim that they were consulted.
The only facts this entry asserts about the reference are the ones `ALEPH-PORT-1703` already
recorded, and they are cited to that entry rather than restated as first-hand.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** `ALEPH-PORT-1703` §3 already reached the verdict
`RE-DERIVE` for this compartment and gave the reason: the crosslink there is a Hookean pair force,
which is the least portable object imaginable. This module was written from Aleph's own registered
`filopodium` contract (`aleph/state/census_actomyosin.py:546`), the five registered connector
contracts that name it, and `docs/design/ENDPOINT_CONTRACT.md`.

What `ALEPH-PORT-1703` §3 identified as the one thing worth carrying was a **discipline**, not a
constant: an empty crosslinker-stiffness slot that *refuses at construction* rather than borrowing a
nearby crosslinker's number. §14 of that entry recorded the obligation this way:

> When a mechanics lane lands, it should raise on an unsupplied crosslinker stiffness, and this
> entry is where that obligation is written down.

**This is that lane, and §9 below is the discharge.** The refusal is stronger here than a
constructor raise, because it is placed on the *consumer* rather than on the constructor: a bundle
can be built, stepped, and force-closed without a bundle effective stiffness, and only the questions
that genuinely need one — the effective bending rigidity, the buckling load, the critical length —
refuse. A constructor raise would have made the whole owner unusable and would have been quietly
worked around by supplying a number.

## 4. Physical or mathematical law represented

Three energy terms, each with its exact analytic gradient. Units µm / pN / pN·µm / s throughout.

**Axial, per segment.** `E = Σ_s (k_s / 2) (L_s − L0_s)² / L0_s`, tension
`T = k (L − L0) / L0`, which pulls a segment's two ends together. Dividing by the rest length is what
makes the modulus resolution-independent: cut a filament into twice as many segments and you get
twice as many springs each twice as stiff, so the end-to-end compliance is unchanged and a
convergence study measures the physics rather than the discretisation.

**Two-sided, and for this compartment that is the whole point.** The registered `mechanical_role`
calls a filopodium a *slender compressive strut*. A tension-only axial law would give a strut that
cannot push, which is the one thing this compartment exists to do — it loads the membrane at its tip.
It is also the term whose absence is invisible: a missing energy term conserves energy perfectly, so
no gradient, closure or conservation control can see it. `assert_bundle_cannot_contract_for_free`
exists for exactly that reason and is called from module code rather than living only in a test.

**Bending, per interior node.** On the unit tangents `t1`, `t2` of the two segments meeting at
interior node `i`: `E_i = (κ_i / 2h_i) |t2 − t1|²` with `h_i` the **rest** Voronoi length
`(L0_prev + L0_next)/2`. Using the rest length decouples bending from stretching, so a uniformly
stretched or contracted straight filament carries exactly zero bending energy and exactly zero
bending force. Small-angle limit `|t2 − t1| = 2 sin(θ/2) → θ → C h` gives `E_i → (κ/2) C² h` and the
sum tends to `∫ (κ/2) C² ds` — Euler–Bernoulli with coefficient `κ` and no leftover factor.

Gradient, with `u = p_i − p_{i−1}`, `v = p_{i+1} − p_i`, `w = t2 − t1`:

```
dE/dv = (κ/h) (I − t2 t2ᵀ) w / |v|
dE/du = −(κ/h) (I − t1 t1ᵀ) w / |u|
```

and `grad_{i−1} = −dE/du`, `grad_i = dE/du − dE/dv`, `grad_{i+1} = dE/dv`, which sum to zero
identically. **`κ` here is the rigidity of a single filament, never of the bundle.** That distinction
is the subject of §5's central refusal.

**Crosslink spacing, between two filaments.** `E = Σ_c (k_c / 2) (d_c − a_c)²` between two
interpolated material points, one on each of two *different* filaments, with `a_c` the crosslinker's
rest interfilament spacing. The gradient is scattered to all four nodes with the **same** weights
that built the two points — the condition under which the chain rule closes, so `F = −grad E` holds
for an element attached *between* nodes, and the four forces sum to zero.

**Two-sided, deliberately.** A crosslinker resists the two filaments separating *and* resists their
approaching. A one-sided version would let a bundle collapse laterally through itself, which is a
bundle with no interfilament spacing at all — and "setting the interfilament spacing" is the
registered representation's own description of what this population does.

### 4b. The one law this module does NOT have, stated as a law

A bundle of `n` filaments crosslinked densely enough to shear-lock behaves as one rod with an
effective rigidity `B_eff ~ n² B_1`; crosslinked loosely, `B_eff ~ n B_1`. The crossover between the
two regimes is set by the crosslinker stiffness. `ALEPH-PORT-1703` §4 derived this and drew the
consequence: for a ten-filament bundle the two scalings differ by a factor of ten in the buckling
load.

The crosslinker stiffness is unsourced. So **`B_eff` is not formed anywhere in this module**, and
every quantity downstream of it refuses:

| Asked for | Answer |
|---|---|
| `FilopodiumBundle.effective_bending_rigidity_pn_um2()` | `FilopodiumUnsourcedParameterError` |
| `FilopodiumBundle.buckling_load_pn()` | `FilopodiumUnsourcedParameterError` |
| `FilopodiumBundle.critical_length_um()` | `FilopodiumUnsourcedParameterError` |
| `FilopodiumBundle.tip_elongation_rate_um_per_s()` | `FilopodiumUnsourcedParameterError` |
| `FilopodiumBundle.root_attachment_rate_per_s()` | `FilopodiumUnsourcedParameterError` |

The first three discharge `FILOPODIUM.unsupported_claims`' first item; the last two discharge its
second. Each raise names the specific unsourced quantity and this entry's §5, so the refusal is
actionable rather than decorative.

### 4c. Two couplings on one edge, and why publishing one is a defect

`filopodium_cytosol_transfer` is `[C]` and carries **two** couplings: immersed drag on the bundle,
and a barbed-end G-actin monomer sink. Its registered `mechanical_interpretation` is explicit that
the second is not decoration:

> a filopodium is a narrow bundle far from the cell body, so its barbed ends are at the far end of a
> long diffusive path and monomer transport, not polymerisation kinetics, is what limits its length.
> Declaring drag alone would remove that limit entirely.

A bundle that publishes only the drag coupling therefore polymerises from an infinite reservoir: its
length is limited by nothing. That is not a small error and it is not conservative — it is the most
flattering error a protrusion model can make, because it deletes the constraint that sets the answer.

**This module publishes both**, and ships `assert_monomer_sink_is_declared_with_drag`, which refuses
a bundle carrying an immersed-drag site with no barbed-end monomer-sink site at its tip. The pairing
is a guard rather than a convention, so it cannot be dropped by whoever writes the connector.

## 4a. Mutation testing of the controls

Three mutants applied to `aleph/vertical/filopodium.py`, the controls run against each, then the
source restored and verified. Recorded because a control suite no mutant can break is decoration.

**Measured after the controls ran** — these numbers cannot exist before the code does, and inventing
them in advance is the failure this ledger's own §0.2.5 discipline exists to prevent.

| # | Mutation | Tests killed | Survived? |
|---|---|---|---|
| M1 | `bending`: rest Voronoi `h` → **current** segment lengths | **4** of 79 | no |
| M2 | `axial`: drop the `/ rest` normalisation (`k Δ²` for `k Δ²/L0`) | **5** of 79 | no |
| M3 | `crosslink_spacing`: scatter the `b` endpoint with the `a` weights | **4** of 79 | no |

Detail, because *which* tests died is the informative part, and in two of the three cases it is not
what a reader would guess.

**M1** killed the per-term finite difference on bending, the total-energy finite difference, **and
both free-contraction guard controls** —
`TestTheAxialLawIsActuallyBound::test_the_shipped_guard_catches_the_unbound_bundle_and_passes_the_bound_one`
and `::test_the_guard_looks_per_filament_and_not_at_the_whole_bundle`. The last two were not
anticipated and they are the interesting kills: `assert_bundle_cannot_contract_for_free` works
because the bending energy is exactly scale-invariant, so contracting a filament changes only the
axial term. With `h` taken from the current geometry the bending energy stops being scale-invariant,
the bending-only bundle no longer looks free, and the guard silently stops detecting the defect it
exists for. **The guard's correctness is a consequence of the bending law's convention**, which is
not obvious from either one alone.

**M2** killed the resolution-independence control and all four parametrisations of the closed-form
control, and **killed nothing else** — not one finite-difference control saw it. That is exactly the
expected shape: dropping the normalisation leaves `F = −grad E` exactly true, so every gradient check
passes on it and only an independent closed form can tell. It is the same lesson `ALEPH-PORT-2801`
§4a records for its own M2, reproduced independently here. `test_a_chain_under_an_end_load_carries_the_closed_form_tension`
was added *because of this run*: before it, M2 was killed by a single test, which is a control suite
one rename away from missing a defect that no gradient check can see.

**M3** killed the crosslink finite difference, the total-energy finite difference, the per-term
closure on the crosslink term, and the total **moment** closure. **Force closure survived it
entirely** — `test_the_internal_forces_sum_to_zero` passed under M3 — and that is not a gap, it is
the reason the moment control exists. Scattering both endpoints with the `a` weights still sums to
`(1−wₐ) + wₐ − (1−wₐ) − wₐ = 0`, so the Newton pair closes exactly while the load lands at the wrong
point along the second filament. A suite with force closure and no moment closure would have called
this correct.

Source restored and verified **byte-identical**: `shasum -a 256` of the restored file equals the
pre-mutation value, `4bd56ffcc42867bfcec4330aeb31ca0ed528e7361a856000c2d5d500832f8cdd`, and the
suite returns to `79 passed`. (`git diff` is silent on this file because it is new and untracked, so
the hash is the evidence and not the diff.)

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| node position | µm | 1e-6 m | finite |
| `filament_axial_modulus_pn` | pN | 1e-12 N | `> 0` (stiffness × length) |
| `filament_bending_rigidity_pn_um2` | pN·µm² | 1e-24 N·m² | `> 0`, **single filament only** |
| `interfilament_spacing_um` | µm | 1e-6 m | `> 0` |
| `spacing_stiffness_pn_per_um` | pN/µm | 1e-6 N/m | `>= 0`, caller-supplied, **unsourced** |
| bundle effective rigidity | pN·µm² | — | **absent. No value, no default, no derivation.** |
| tip elongation rate | µm/s | 1e-6 m/s | **absent. Refused.** |
| root attachment rate | 1/s | 1/s | **absent. Refused.** |
| energy | pN·µm | 1e-18 J | `>= 0` |

Singular and boundary cases, each with the behaviour required:

- **Collapsed segment** (`L = 0`): refused with `ValueError`. The direction is genuinely undefined and
  a silent epsilon returns a finite force whose direction came from round-off.
- **Zero-length crosslink**: refused, same reason.
- **`rest_length <= 0`**: refused; the normalisation would divide by zero.
- **Self-crosslink** (`filament_a == filament_b`): refused with `FilopodiumRoleError`. Between two
  points of one filament it is a second, undeclared axial law with its own rest length, and the
  filament's stiffness stops being a function of its card.
- **A detached root.** Legal — a retracting or newly nucleated bundle. It transmits **exactly** `0.0`,
  not a small number, and that is testable with `==` because the gate returns an identical zero on
  the whole branch.
- **A bundle with no filament reaching the tip.** Refused: the load per barbed end is `|F| / n` and
  `n = 0` is not "no load", it is a division by zero dressed as a protrusion.
- **A bundle of one filament.** Legal and degenerate: the crosslink topology is empty and the
  bundle is one filament. A consumer that assumes `n >= 2` must say so.
- **A crosslinker stiffness or effective rigidity requested with no source.** Refused, never
  defaulted, never borrowed from a different crosslinker.
- **A target filopodium count.** Refused by `assert_population_count_is_an_output`, and no such
  parameter exists to be set.
- **Empty element set**: returns exactly `(0.0, zeros_like(positions))`, not a crash.
- **A role published in the wrong zone** (a membrane-tip role at the root, a nascent-adhesion role at
  the tip): refused with `FilopodiumRoleError`.
- **A stale endpoint handle** held across a committed topology change: refused with
  `FilopodiumStaleHandleError`.

Invariants, each with the control that asserts it:

- **I1.** `F = −grad E` per term, observed order > 1.6 —
  `tests/vertical/test_filopodium_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order`
- **I2.** Force and moment closure against the **constituent** scale, per term —
  `::TestClosureAndInvariance::test_every_term_closes_on_its_own`
- **I3.** Translation invariance — `::TestClosureAndInvariance::test_translating_the_whole_bundle_changes_no_energy`
- **I4.** Resolution-independent axial modulus —
  `::TestTheAxialLawIsActuallyBound::test_end_to_end_stiffness_does_not_depend_on_node_count`
- **I5.** No term mutates its inputs, bit-identically (`np.array_equal`) —
  `::TestClosureAndInvariance::test_no_term_mutates_its_arguments`
- **I6.** Declared state keys equal the registered contract —
  `::TestOwnershipAndTheAcceptedStepContract::test_the_declared_keys_are_exactly_the_registered_ones`
- **I7.** No state key of this owner is declared by any other registered owner —
  `tests/state/test_census_key_ownership.py::test_no_state_key_has_two_registered_owners`
- **I8.** A handle issued before a committed topology change refuses to scatter
  (`ENDPOINT_CONTRACT.md` I3) —
  `::TestOwnershipAndTheAcceptedStepContract::test_a_handle_held_across_a_commit_is_refused_and_a_fresh_one_works`
- **I9.** Rollback restores bit-identically —
  `::TestOwnershipAndTheAcceptedStepContract::test_rollback_restores_the_bundle_bit_identically`
- **I10.** A detached root transmits exactly zero —
  `::TestOwnershipAndTheAcceptedStepContract::test_detaching_the_root_is_a_state_change_and_not_a_topology_edit`

## 6. Source evidence class and known retractions

No claim is made about the reference implementation, because none of it was read by this lane. The
facts `ALEPH-PORT-1703` §6 recorded about it — that its filopodium is rated declared-only with zero
device execution, that its motor connector has no bind target, and that its crosslinker stiffness
slot is empty with a construction-time raise — are **cited to that entry**, at that entry's evidence
class, and are not re-asserted here as first-hand observations.

The registered `filopodium` contract this module was written against carries
`citation_status: UNSOURCED`, and nothing here changes that. `BundleCard` has **no defaults at all**,
so there is no unsourced number in this module for a reader to mistake for a measurement; every
magnitude in the controls is supplied by the control and is chosen to make the algebra sharp.

Retractions: none known. Nothing was inherited, so there is nothing to retract.

## 7. Independent oracle or derivation

Three, none of which is the reference implementation:

1. **Central finite differences** of the module's own energies — an oracle for the gradient claim
   that is independent of every constant in the law.
2. **The axial closed form** `ΔL = F L0 / k` for a chain under an end load, derived in §4 from the
   energy alone. This catches the defect the finite difference cannot: a rest-length normalisation
   error keeps `F = −grad E` exactly true while making the modulus resolution-dependent — mutant M2
   in §4a, which no finite-difference control killed.
3. **Force and moment closure on an isolated bundle**, reported against the constituent scale
   `Σ|f_i|` and `Σ|r_i||f_i|`. Never against the resultant: a resultant-based tolerance gets
   *stricter the more correct the physics is*, because a correct Newton pair cancels, and `PLAN.md`
   records that defect rejecting forty thousand consecutive steps of a good relaxation while
   reporting a plausible tension.

**The oracle this module deliberately does NOT use is Euler buckling.** `ALEPH-PORT-1703` §7
proposed measuring the bundle's effective rigidity from a buckling sweep and comparing it against the
`n` versus `n²` scalings, which would turn the missing crosslinker constant into a measurable. That
is a good idea and it is **not implemented here**. It would require the module to form a bundle
effective rigidity, which §4b refuses to do — so the buckling oracle and the refusal are in direct
tension, and resolving that is a decision above this lane's level. §14 records it as an absence.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_filopodium_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order` | Per term, `energy > 0` first, then `max\|F + ∇E\| / max\|F\|` at steps 4e-4 / 2e-4 / 1e-4 µm with observed order > 1.6 |
| Positive | `tests/vertical/test_filopodium_controls.py::TestTheAxialLawIsActuallyBound::test_uniformly_contracting_the_bundle_costs_energy` | Uniform contraction of a straight bundle costs energy, the bending term contributes exactly zero to it, and the axial term contributes all of it |
| Positive | `tests/vertical/test_filopodium_controls.py::TestTheSink::test_the_tip_load_is_divided_among_the_barbed_ends_that_reach_it` | The per-barbed-end load is the tip load divided by the number of filaments reaching the tip, and a bundle with none refuses |

**Measured, this session, on this tree** (`/Users/sw1/miniconda3/envs/aleph/bin/python`). Every
number below was produced by a run in this session; none is quoted from memory. See §13 for the
suite line.

Finite-difference agreement, steps 4e-4 / 2e-4 / 1e-4 µm, on the three-filament oblique fixture:

| Term | Energy (pN·µm) | Force scale (pN) | max\|F+∇E\| at 1e-4 (pN) | relative | observed order |
|---|---|---|---|---|---|
| axial | 1.376267 | 19.744 | 4.718e-07 | 2.389e-08 | **2.000** |
| bending | 0.158354 | 5.099 | 2.598e-07 | 5.095e-08 | **2.000** |
| crosslink spacing | 3.683306 | 26.599 | 5.379e-07 | 2.022e-08 | **2.000** |
| total | 5.217926 | 30.254 | 6.050e-07 | 2.000e-08 | **2.000** |

Closure on the isolated bundle (resultant ÷ constituent scale):

| Term | force | moment |
|---|---|---|
| axial | 3.88e-19 | 3.25e-18 |
| bending | 5.22e-18 | 4.50e-17 |
| crosslink spacing | 2.30e-17 | 1.24e-17 |
| total | **1.02e-17** | **7.35e-18** |

Resolution independence of the axial modulus at identical end-to-end strain, node counts 3 / 5 / 9 /
17: axial energy `12.000000000000021` pN·µm at every resolution, relative spread **exactly
`0.0`** — the four sums are the same floats in a different order and at uniform spacing they
associate identically.

The closed-form control, same four resolutions: a chain of rest length 3.0 µm and modulus 800 pN
stretched to `ΔL = F L0 / k = 0.15` µm carries `40.0` pN in every segment, the two end nodes carry
`+40.0` and `−40.0` pN, and the interior residual is at worst `4.7e-14` of the applied load. Under
mutant M2 the same configuration reports `120/S` pN per segment, where `S` is the segment count.

The free-contraction guard, on the mixed fixture: the smallest self-energy difference between the
two contraction depths is `301.323` pN·µm. Under `omit_axial_law=True` it is exactly `0.0` and the
guard raises.

The flipped-sign negative control: relative finite-difference disagreement **2.000**, against the
0.5 the control demands — and the axial energy is bit-identical between the healthy and flipped
branches (`1.3762669970268657` both), which is the whole point.

## 9. Deliberately failing negative control

Every break below is a **flag on the shipped code**, default `False`, so the control drives the real
path rather than a hand-edited copy. A negative control against a duplicate proves only that the
duplicate is broken.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_filopodium_controls.py::TestTheNegativeControlsFail::test_a_flipped_axial_sign_is_caught_by_the_gradient_test` | With `flip_axial_sign=True` the finite-difference disagreement rises above **0.5 relative** (measured 2.000), while the **energy is bit-unchanged** — the failure an energy-only check cannot see |
| Negative (must fail) | `tests/vertical/test_filopodium_controls.py::TestTheAxialLawIsActuallyBound::test_the_same_contraction_is_free_when_the_axial_law_is_unbound` | With `omit_axial_law=True` a uniform contraction costs **exactly `0.0`** and every internal force is **exactly `0.0`**; with the law bound the same contraction costs a strictly positive energy |
| Negative (must fail) | `tests/vertical/test_filopodium_controls.py::TestTheNegativeControlsFail::test_a_one_sided_spacing_law_lets_the_bundle_collapse_laterally` | With `allow_one_sided_spacing=True` two crosslinked filaments pushed closer than the rest spacing store **exactly `0.0`** energy and feel **exactly `0.0`** force; with the two-sided law they are pushed apart |
| Negative (must fail) | `tests/vertical/test_filopodium_controls.py::TestTheNegativeControlsFail::test_a_detached_root_that_still_transmits_load_is_the_shared_node_defect` | With `allow_detached_root_to_transmit_load=True` a detached root transmits the full load; with the gate in place it transmits **exactly `0.0`** — which is the mechanical content of "the root is a connector, not a weld" |
| Negative (must fail) | `tests/vertical/test_filopodium_controls.py::TestTheRefusalsAreRealRefusals::test_drag_declared_without_the_monomer_sink_is_refused` | A bundle publishing immersed drag with no barbed-end monomer sink is refused by `assert_monomer_sink_is_declared_with_drag`; adding the sink makes it pass |

Each is asserted in **both** directions on purpose. Asserting only that the healthy branch is zero
would pass for a function that returns zero unconditionally, which is the shape of the `PLAN.md`
§6.1 defect where a cold cache read made a deliberately broken strut look innocent.

## 10. Numerical and precision envelope

Working and accumulation precision: float64 throughout. No reduced-precision path exists.

Finite-difference steps 4e-4 / 2e-4 / 1e-4 µm on coordinates of order 0.1–4 µm. The floor for a
central difference is `eps^(1/3) ≈ 6e-6` **relative**; the smallest step is well above it. Below that
floor cancellation in the energy evaluation dominates, the apparent convergence order goes negative,
and a correct gradient is indistinguishable from a wrong one — the afternoon `PLAN.md` §6.1 records.

Tolerances and why they are not looser: per-term relative finite-difference agreement is asserted at
`1e-6` against a measured worst case of `5.1e-08`, ~20× of headroom — tight enough that mutant M3 (a
weight swap, which perturbs the gradient at the 1e-1 level) cannot pass. Closure is asserted at
`1e-12` of the constituent scale against a measured worst case of `4.5e-17`; the gap is deliberate, because closure
is an identity of the algebra and its residual is set by summation order rather than by the physics.
Exact-zero properties — the detached root, the one-sided spacing branch, the unbound axial law — are
asserted with `==` and not with a tolerance, because each is an identical zero on a whole branch
rather than a small number near a gate, and a tolerance there would only hide a sign error.

Outside the envelope: every singular case in §5 refuses with a typed exception rather than degrading
silently. **There is no input range over which this module returns a plausible number it cannot stand
behind**, and the quantities it cannot stand behind at all — the bundle effective rigidity, the
buckling load, the critical length, the tip and root rates — are not clamped, defaulted or
approximated. They raise.

## 11. Production-backend residency and transfer

Host, numpy, float64. The module imports `EntityCensus`, `Phase`, `StepContext` and `WorkReceipt`
from `aleph.runtime.participant` and **uses all four** — `accumulate` reaches the backend through
`StepContext.backend` and returns a `WorkReceipt` with a non-zero declared operation count, so a
participant that returned early is distinguishable from one that computed the forces. No GPU work of
any kind was run by this lane, and none was requested.

The three terms are `O(S)` and `O(T)` gather/scatter kernels with no host round-trip per step. The
crosslink spacing term is a short-range interaction **within one owner's array**, so it needs no
cross-ownership transfer at all. The root is different, and that is the point of the compartment: it
needs a two-array adjoint scatter, which is precisely why it is a connector rather than a weld.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered this module, because no reference file was read by
this lane.
`tests/vertical/test_filopodium_controls.py::TestTheDeclarationsAgreeWithTheCode::test_no_reference_project_vocabulary_survives`
asserts the module text carries none of the provider's identifiers, and
`::test_the_module_imports_nothing_from_validation` asserts the `aleph/** → validation/**` firewall.

Two pieces of prose were **deliberately not written**, and naming them is the useful part:

- **No docstring describes an unimplemented feature in the present tense.** `ALEPH-PORT-2801` §12
  records `ecm.py` shipping a module docstring that described five absent symbols as though they
  existed, and `__all__` exporting all five. Here `__all__` is pinned in *both* directions by
  `::TestTheDeclarationsAgreeWithTheCode::test_every_exported_name_resolves` and
  `::test_every_public_symbol_is_exported`, and everything absent is absent on purpose and named in
  §14.
- **No "provisional" magnitude with a comment promising to replace it.** `BundleCard` has no
  defaults, so there is no number in this module at all that a later reader could mistake for a
  measurement. That is the same discipline `cytosol.py`'s `PoroelasticCard` carries and it is
  enforced the same way — a control asserts the factory name does not exist.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 79 controls pass (`tests/vertical/test_filopodium_controls.py`, `79 passed in 0.24s`). Three mutants applied, killing 4 / 5 / 4 tests respectively; none survived; source restored byte-identically (SHA-256 match). Status stays `PROPOSED`: nothing here has been reviewed by a human, this owner has never been stepped inside a world, and it wires zero connectors. |
| Reviewer | Agent-proposed. **Unratified.** No PI review. |
| Rollback | Delete `aleph/vertical/filopodium.py` and `tests/vertical/test_filopodium_controls.py`. Nothing breaks: no module in the tree imports either, the package `__init__` does not name them, and no connector reaches them. |

## 14. Honest limits — what this entry does NOT establish

**Crosslinker-limited buckling is not resolved, and this module cannot resolve it.** The registered
`approximation` says the bundle's crosslinkers set an *effective* bending rigidity rather than being
resolved as individual bonds under load. This module does neither of those two things exactly: it
resolves the crosslinkers as **spacing** elements — two-sided springs about a rest interfilament
spacing, with a caller-supplied stiffness — and it forms **no bundle effective rigidity at all**. The
two readings agree on the observable that matters, which is why the deviation is reported rather than
papered over: **under either reading, crosslinker-limited buckling is not resolved**, because the
`n` versus `n²` crossover of §4b is never evaluated. What this module can say about a bent bundle is
the sum of its filaments' own bending energies plus the spacing energy of its crosslinks. What it
cannot say is what load buckles it, and it refuses that question rather than answering it from a
rigidity nobody sourced.

**The deviation from the registered `approximation` is a finding, not a licence.** The contract's
wording describes a homogenised implementation; this one is resolved-in-spacing and
refused-in-rigidity. If the PI reads that as a contract disagreement rather than a compatible
narrowing, the contract is the artefact that should move and this entry is not the place to move it.

**The Euler-buckling oracle proposed by `ALEPH-PORT-1703` §7 is not implemented.** It is the natural
independent check on this compartment and it is in direct tension with the refusal above: measuring
`B_eff` from a buckling sweep requires forming `B_eff`. Whether the right resolution is to implement
the sweep as a *measurement* that returns a number without the module ever carrying one as a
*parameter* is a real design question and it is not answered here.

**Wired connector count contributed by this module: zero.** Five registered connectors name
`filopodium` — `filopodium_membrane_tip`, `filopodium_cortex_root`, `filopodium_nascent_fa`,
`nmii_filopodium_motor`, `filopodium_cytosol_transfer` — and this module publishes an endpoint for
every one of them. Publishing an endpoint is not wiring a connector, and `wiring.py` is untouched, so
the wired count is unchanged by this entry.

**`filopodium_cytosol_transfer` publishes both couplings, but only one of them is mechanics.** The
immersed-drag endpoint is a real quadrature point with a position and a force sink. The barbed-end
monomer sink publishes a **count and a position**, and nothing else: the module has no monomer field,
no diffusion, no polymerisation rate and no consumption law, because `tip elongation under load` is
on the evidence-owed list. So the transport limit the connector's contract calls decisive is
**declared and unquantified**, and `tip_elongation_rate_um_per_s()` refuses rather than returning a
number. The absence is visible; the limit is not modelled.

**Physical magnitudes are not established, at all.** Evidence rung `ANALYTIC_ORACLE`, quantitative
status `BLOCKED`. No modulus, rigidity, spacing or stiffness anywhere in the controls is a filopodial
number; they are chosen to make the algebra sharp. `citation_status` for the `filopodium` contract
remains `UNSOURCED` and this entry does not improve it.

**This owner has never been stepped inside a world.** It is not in the assembled vertical, it is not
exported from the package `__init__`, and no relaxation has ever been run with it. The
controls exercise one candidate step's worth of phases directly. `accumulate` is verified to produce
a receipt and to refuse an unscheduled phase; it is **not** verified to converge, to be stable at any
step size, or to compose with any other owner.

**The step-size consequence of the missing stiffness is recorded and not measured.**
`ALEPH-PORT-1703` §10 noted that a bundle crosslink at physiological spacing is a stiff short-range
spring and that the explicit stability limit `dt < 2/(M·k)` scales inversely with it — so the missing
crosslinker stiffness sets the cost as well as the physics. This module takes no position on step
size and runs no stability study.

**The reference implementation is unaudited by this lane**, and deliberately so: no file under
`/Users/sw1/ffn_cellsim` was opened, so whether Aleph's re-derivation agrees or disagrees with it is
unknown here and this entry does not guess.

**Not run:** any GPU job (zero; no authorization was sought or held), and the full suite. Only
`tests/vertical/test_filopodium_controls.py`, `tests/state/` and `tests/ports/` were exercised —
other sessions are writing this tree concurrently, and `CLAUDE.md` §1 says to report a foreign
breakage rather than fix it, which means not conflating one with mine. Two pre-existing failures in
`tests/ports/` are reported in this lane's summary and were **not** touched.

**§4a and §8 were filled in after the controls ran**, on `2026-07-31`; the rest of the entry was
written on `2026-07-30` before the code, per `PLAN.md` §0.2.5, and it constrained the design: the refusal table in §4b, the
both-couplings guard in §4c, the zone allowlist in §5 and the four deliberate-break flags in §9 were
all specified here first and then implemented. A measured number cannot be written before the
measurement exists, and writing one would be the precise failure this discipline prevents.
