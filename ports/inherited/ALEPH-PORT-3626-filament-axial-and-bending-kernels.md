# ALEPH-PORT-3626 — the filament axial and bending laws: six owner modules, one kernel family

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3626` |
| Lane | `G10 owner-interior kernels` (session `e874a7fe`) |
| Status | `PROPOSED` |
| Written | `2026-08-04` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Aleph API | `aleph/runtime/law_kernels.py`, `aleph/runtime/law_cases.py` — both append-only |
| Exists because | 107 kernels pass parity on the RTX 4090 (`ALEPH-PORT-3619`) and **every one of them is a connector or a membrane/nucleus law.** The filament *owners* — the things that actually hold a cell up — had no kernel at all. This entry is the first port of an owner interior. |

---

## 1. Aleph API

Two new law-case builders, `filament_axial_law_case` and `filament_bending_law_case`, registered in
`LAW_CASE_BUILDERS` and in `_true_variants()`; two new variant enums, `FilamentAxialVariant` and
`FilamentBendingVariant`; **fourteen new kernels** in `aleph/runtime/law_kernels.py` — four geometry
stages, five axial and five bending, of which **eight are wrong on purpose**.

**This entry authorises no change to any module under `aleph/vertical/`.** The six NumPy laws are
the frozen parity reference; they were read, driven and never edited. `git status` on
`aleph/vertical/` is empty.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing was taken from `/Users/sw1/ffn_cellsim` (READ ONLY, and not read for this entry) |
| Source commit | not applicable — no provider file was consulted; the reference read is this repository at `9936d20` |
| Source path | internal: `aleph/vertical/{cortex_filaments,ecm,sf_arc,microtubule,filopodium,lamellipodium}.py` |
| Source symbol(s) | `axial_energy_and_forces`, `bending_energy_and_forces` — **six independent definitions of each** |
| Read from | the working tree at `9936d20`, verified clean on those six files before reading |
| Working tree == commit? | yes for all six — `git status aleph/vertical/` is empty |

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported from the provider.** This is `RE-DERIVED` in the strict sense:
the laws are Aleph's own, written by Aleph's own owner lanes, and the kernel is a transcription of
Aleph's float64 NumPy into float32 Warp with the association order matched deliberately so that any
disagreement the harness reports is the float32 *representation* and not a gratuitously different
order of operations.

## 4. Physical or mathematical law represented

### 4.1 Axial — `FilamentAxialVariant`

Per segment `s` joining nodes `a` and `b`, on the rest length `L0` and the axial modulus `k`:

```
d = p_b - p_a      L = |d|      x = L - L0            SIGNED. No clamp on this line.
E = sum_s (k/2) x^2 / L0
T = (k x) / L0     f_a = +T d/L     f_b = -T d/L
```

**The division by `L0` is the design and not a normalisation of convenience.** It makes the modulus
resolution-independent: a filament cut into twice as many segments has twice as many springs each
twice as stiff, and the same end-to-end stiffness. Without it, refining a mesh silently changes the
mechanics and every convergence study measures the discretisation instead of the physics.

**Two-sided, deliberately.** A one-sided version restores exactly the free contraction mode these
six owners exist to remove.

### 4.2 Bending — `FilamentBendingVariant`

Per triple `(i-1, i, i+1)`, with `h` the **rest** Voronoi length `(L0_prev + L0_next)/2`:

```
u = p_i - p_{i-1}   v = p_{i+1} - p_i   t1 = u/|u|   t2 = v/|v|   w = t2 - t1
E = sum (kappa / 2h) |w|^2
dE/dv = (kappa/h) (I - t2 t2^T) w / |v|
dE/du = -(kappa/h) (I - t1 t1^T) w / |u|
f = ( +dE/du , dE/dv - dE/du , -dE/dv )                 sums to zero identically
```

`h` is the **rest** length and not the current one, which is what keeps bending and stretching
decoupled: a uniformly stretched straight filament carries exactly zero bending energy and exactly
zero bending force.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `rest_lengths_um` | µm | m | `> 0`, refused otherwise by the frozen law |
| `axial_modulus_pn` | pN | N | any real; negative is not refused and is not physical |
| `voronoi_length_um` | µm | m | `> 0`, refused otherwise |
| `bending_rigidity_pn_um2` | pN·µm² | N·m² | any real |
| `energy_pn_um` | pN·µm | J | — |
| `forces_pn` | pN | N | — |

**Singular cases.** A zero-length segment is refused by the reference (`ValueError`) and would
divide by zero in the kernel; `filament_state()` cannot produce one. A kernel cannot raise, so the
domain is enforced on the host before launch, exactly as `ALEPH-PORT-3604` §5 records for the
contact.

**Invariants.** Both laws are internal, so `sum(f) = 0` identically — a property the kernel holds on
its own with no reference at all, and two controls check it. The bending law additionally has an
**exact** anchor: a straight filament returns `0.0` on both channels, asserted with `==`.

## 6. Source evidence class and known retractions

No provider source, so no provider retraction applies. The Aleph-internal reference carries one
claim this entry checked rather than inherited: `cortex_filaments.bending_energy_and_forces`'
docstring says it is *"pinned to it numerically by a control"* with respect to `sf_arc`. That claim
is true and is now measured over 400 configurations rather than asserted at one — see §7.

## 7. Independent oracle or derivation

**The oracle is that the law exists six times and the six agree.**

`axial_energy_and_forces` and `bending_energy_and_forces` are defined independently in six owner
modules with six **distinct source bodies** (`inspect.getsource` gives six different digests, and
six different function objects — asserted, because if they were shared the agreement below would be
a tautology).

Measured over 400 random configurations, chain lengths 3–24 nodes, positions ~N(0, 3 µm), rest
lengths 0.05–4 µm, moduli 1–2000 pN, rigidities 0.1–300 pN·µm²:

| Comparison | Result |
|---|---|
| axial, energy identical under `==` and forces under `np.array_equal` | **2000 / 2000 owner-pairs agree** |
| bending, same test | **2000 / 2000 owner-pairs agree** |
| tension-only branch, `ecm` vs `sf_arc` (the two that take a per-element mask) | **400 / 400**, with **237 slack segments** actually exercised |

So **one kernel family serves six owners**, and the number of kernels and the number of owners are
two separate figures, quoted separately, for the reason `ALEPH-PORT-3606` gives: one implementation
serving many rows makes any single count a misrepresentation.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/runtime/test_filament_law_parity.py::test_the_axial_kernel_matches_the_frozen_law_within_budget` | axial parity inside the declared ULP budget |
| Positive | `tests/runtime/test_filament_law_parity.py::test_the_bending_kernel_matches_the_frozen_law_within_budget` | bending parity inside the declared budget |
| Positive | `tests/runtime/test_filament_law_parity.py::test_the_six_owners_axial_law_is_one_law_measured_bitwise` | the §7 measurement, 2000/2000 |
| Positive | `tests/runtime/test_filament_law_parity.py::test_the_six_owners_bending_law_is_one_law_measured_bitwise` | the §7 measurement, 2000/2000 |
| Positive | `tests/runtime/test_filament_law_parity.py::test_a_straight_filament_has_exactly_zero_bending_energy_and_force` | the exact anchor, with `==` |
| Positive | `tests/runtime/test_filament_law_parity.py::test_the_axial_kernel_forces_sum_to_zero` | closure, no reference involved |
| Positive | `tests/runtime/test_filament_law_parity.py::test_the_bending_kernel_forces_sum_to_zero` | closure, no reference involved |
| Positive | `tests/runtime/test_filament_law_parity.py::test_the_budget_covers_a_sweep_of_states_and_not_just_the_default_one` | the budget holds over twelve states, not one — see §13a |

## 9. Deliberately failing negative control

**Eight wrong kernels. Eight caught. 8 / 8.**

| Control (must fail) | Variant | Caught on | Bit-identical to TRUE on |
|---|---|---|---|
| `tests/runtime/test_filament_law_parity.py::test_an_axial_kernel_that_dropped_the_half_in_its_energy_is_caught` | `filament_axial_wrong_energy_half` | `energy_pn_um` | `forces_pn` |
| `tests/runtime/test_filament_law_parity.py::test_an_axial_kernel_that_dropped_the_rest_length_normalisation_is_caught` | `filament_axial_wrong_no_rest_normalisation` | both | none |
| `tests/runtime/test_filament_law_parity.py::test_an_axial_kernel_that_became_tension_only_is_caught` | `filament_axial_wrong_tension_only` | both | none |
| `tests/runtime/test_filament_law_parity.py::test_an_axial_kernel_with_a_flipped_sign_is_caught` | `filament_axial_wrong_sign_flip` | `forces_pn` | `energy_pn_um` |
| `tests/runtime/test_filament_law_parity.py::test_a_bending_kernel_that_dropped_the_half_in_its_energy_is_caught` | `filament_bending_wrong_energy_half` | `energy_pn_um` | `forces_pn` |
| `tests/runtime/test_filament_law_parity.py::test_a_bending_kernel_that_dropped_the_tangent_projection_is_caught` | `filament_bending_wrong_no_projection` | `forces_pn` | `energy_pn_um` |
| `tests/runtime/test_filament_law_parity.py::test_a_bending_kernel_that_used_the_current_voronoi_length_is_caught` | `filament_bending_wrong_current_voronoi` | both | none |
| `tests/runtime/test_filament_law_parity.py::test_a_bending_kernel_that_mis_signed_the_middle_node_is_caught` | `filament_bending_wrong_middle_sign` | `forces_pn` | `energy_pn_um` |

Plus one second, independent channel on the last of those:
`tests/runtime/test_filament_law_parity.py::test_the_mis_signed_middle_node_also_breaks_force_closure`
— closure needs no reference at all, so a port that lost its oracle would still be caught by it.

**`DECLARED_BLIND_VARIANTS` gains no row.** Four of the eight are bit-identical to the true kernel in
**one** channel, which is the point of grading two channels and is not blindness; none is
bit-identical in both. Under `tests/runtime/test_declared_blind_spots.py`'s rule — *measure why a
mutant is invisible before declaring it blind* — there is nothing here to declare, and declaring
anything would be the waiver that guard exists to refuse.

### 9b. The mutant that no gradient, closure or energy check can see

`filament_axial_wrong_no_rest_normalisation` drops `/ L0` from **both** the energy and the tension.
`F = -grad E` stays exactly true, so the mutant is **perfectly conservative**: a finite-difference
gradient probe, a force-closure check and an energy ledger are all structurally unable to see it.
What it destroys is the resolution-independence §4.1 describes — the whole reason the division is
there. Only parity against the frozen law catches it. This is the clearest argument in the port for
why the parity harness exists at all, and it is the same shape as `ALEPH-PORT-3611` §9b's
`WRONG_OMIT_RESTIFFENING`.

### 9c. The mutant the exact anchor cannot catch

`filament_bending_wrong_current_voronoi` replaces the rest Voronoi length with the current one.
It **preserves the straight-filament exact zero** — `w = 0` either way — so the anchor that guards
this law is blind to it and only a bent state grades it. `ALEPH-PORT-3605` §13a M9's lesson,
arriving in a second law: an exact anchor is a strong control and it is not a complete one.

## 10. Numerical and precision envelope

Measured on **warp's CPU device**, float32 compute, float64 energy accumulation, over the twelve
states in `FILAMENT_BUDGET_SWEEP_SEEDS`. Maximum and median ULP against the frozen float64 law:

| Law | Position mode | energy max | energy med | forces max | forces med |
|---|---|---|---|---|---|
| axial | `GLOBAL_F32` | **17.277** | 6.130 | **26.049** | 20.819 |
| axial | `LOCAL_F32` | 3.665 | 0.627 | 5.185 | 3.665 |
| axial | `POSITIONS_F64` | 2.687 | 0.626 | 3.236 | 1.577 |
| bending | `GLOBAL_F32` | 3.243 | 0.923 | 18.140 | 12.186 |
| bending | `LOCAL_F32` | 0.500 | 0.182 | 7.345 | 3.350 |
| bending | `POSITIONS_F64` | 0.334 | 0.125 | 2.920 | 1.597 |

Declared budgets, set at ~1.4× the measured maximum: axial `energy 24.0 / forces 40.0`, bending
`energy 8.0 / forces 32.0`.

**Force closure on the true kernels**, no reference involved:

| Law | `|sum f|` | peak force | ratio |
|---|---|---|---|
| axial | 2.5678e-05 pN | 2.2699e+02 pN | 1.131e-07 |
| bending | 4.7197e-05 pN | 2.5667e+02 pN | 1.839e-07 |

Both are at the float32 accumulation floor for 23 and 22 elements respectively, not at a chosen
tolerance.

**Reference-input sensitivity** (`tests/runtime/test_reference_input_guard.py`'s own measurement,
applied to these two cases): axial moves `energy 7.197 / forces 21.257` ULP and bending
`energy 0.985 / forces 7.752` ULP when the reference's input positions are rounded to float32. Both
clear the guard's `> 1 ULP` requirement, so **neither new case joins the recorded conflict** that
`HANDOFF.md` §C-0 holds open for linear isotropic drag.

### 10a. A measurement for the open position-precision proposal, and a counterpoint to `-3605`

`ALEPH-PORT-3605` §13.4 reported `POSITIONS_F64` **worse** than `LOCAL_F32` on the enclosed-volume
law, because that law's cancellation lives in a *reduction* while option D of
`docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` is defined as float64 arithmetic
only in the subtraction that forms edge vectors and site separations.

The filament laws are the opposite case and they behave the opposite way. Their cancellation lives
**entirely** in the segment subtraction `p_b - p_a`, and there is no global reduction in the force
channel at all — so option D applies exactly as written, and `POSITIONS_F64` is the best of the three
modes on both laws. `LOCAL_F32` recovers most of the same benefit (axial forces 5.185 against
3.236 max ULP) for a host-side centroid subtraction and no float64 storage on the device.

**This lane answers none of that proposal's questions.** It contributes one law family's numbers to
§7 of it, where lane ledgers are asked to report mode costs, and it records that `-3605`'s finding
does not generalise — as `-3605` itself was careful to say, naming one law against one recommendation.

## 11. Production-backend residency and transfer

Every scatter destination and every reduced array comes from `backend.zeros` / `backend.array`, per
`ALEPH-PORT-3603` §11. Connectivity (`segments`, `triples`) is allocated with `wp.array` as int32
and is deliberately **not** backend-allocated: it is not a compute element, carries no precision
question, and is not a scatter destination.

One scatter per law. The axial kernel writes `(2S, 3)` pair forces scattered by
`segments.reshape(-1)`; the bending kernel writes `(3T, 3)` stencil forces scattered by
`triples.reshape(-1)`. The energy is reduced by `Backend.sum` — ordered by construction — so
`scatter_assembled` names `forces_pn` only and the atomic ULP floor does not widen the energy budget.

## 12. Comments and docstrings to discard

None. No provider prose was read, so none survives. Every comment in the fourteen new kernels was
written for Aleph's situation, and the wrong kernels each carry the defect they stand for in their
own comment, so a reader meeting one in a diff learns what it is for rather than mistaking it for a
bug.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | `tests/runtime/test_filament_law_parity.py`: **27 passed**, 2026-08-04, on warp's CPU device. `tests/runtime`: **674 passed, 7 skipped, 1 failed** — the one failure is the pre-existing reference-input guard on `isotropic_drag` (0.443 ULP), verified foreign by re-running with this lane's two files stashed. Mutation study **8 / 8 killed**. ULP table §10. `aleph/vertical/` byte-for-byte untouched. |
| Reviewer | **agent-proposed; unratified.** No `decided_by` field appears in this entry. |
| Rollback | Remove the two `LAW_CASE_BUILDERS` rows, the two `_true_variants()` rows, the two enums, the fourteen kernels and their registry rows, and `tests/runtime/test_filament_law_parity.py`. Nothing else imports them; no `Binding` was added to `aleph/vertical/wiring.py`, because this port adds no connector. |

### 13a. A defect this port found in its own budget, before shipping it

The first `DEFAULT_FILAMENT_AXIAL_ULP_BUDGET` was `energy 8.0 / forces 32.0`, set from the default
state, where the axial energy measures **7.909 ULP**. That reads as a budget with a 1% margin.

It was a coincidence. Over the twelve states in `FILAMENT_BUDGET_SWEEP_SEEDS` the same channel
reaches **17.277 ULP**, and **four of the twelve exceed 8.0** (seeds 2, 101, 5150, 8888 at 11.556,
17.100, 11.101 and 17.277). The budget was already violated by states `filament_state()` itself can
produce, and nothing would have said so until a sweep happened to draw one.

`PLAN.md` §8 records this exact lesson from lane L13 — *"Laplace's law recovers to 0.0004% at one
subdivision level, one tolerance, one seed, one contact stiffness. That is a single point, and a
single point is not a numerical result."* A ULP budget read off the default state is that same
mistake wearing a tolerance. The budget is now set from the distribution, and
`test_the_budget_covers_a_sweep_of_states_and_not_just_the_default_one` is the control; it was
verified able to fail by re-running it against the original budget, where it reports the four
violations above.

## 14. Honest limits

1. **No GPU. Every figure above is warp's CPU device.** No PI citation naming a card and a duration
   was given to this session, and an agent may never write one. What `ALEPH-PORT-3619` established
   for the 107 existing kernels — that CPU-Warp and CUDA-Warp agree on the vast majority and differ
   in ways worth reporting — is **not** established for these fourteen, and this entry claims
   nothing about them on a device. The CUDA row of §10 is `UNVERIFIED`.
2. **Six owners agree bitwise; that is a measurement over 400 configurations of a chain topology.**
   The states swept are open chains with sequential segments and triples. A branched network
   (`lamellipodium`'s junctions, `ecm_network`'s sheet) has the same per-element law but a different
   scatter pattern, and no figure here covers a node shared by more than two segments.
3. **The port covers two of the laws these six owners carry, not all of them.** Untouched and still
   without a kernel: `cortex_filaments.crosslink_energy_and_forces` and `steric_energy_and_forces`,
   `ecm.crosslink_energy_and_forces`, `sf_arc.joint_energy_and_forces`,
   `microtubule.minus_end_anchor_energy_and_forces`, `filopodium.crosslink_spacing_energy_and_forces`,
   `lamellipodium.branch_angle_energy_and_forces` and `junction_link_energy_and_forces`.
4. **Nine further owner modules have no kernel at all** and were named, not ported, by this lane's
   derivation: `cortex`, `cortex_surface_coupling`, `cytosol` (session S2 is rewriting its solver, so
   porting against it would be parity against a moving reference), `ecm_network`, `pressure`, and the
   composition/driver modules `assembly` and `relax`, which carry no force law.
5. **No `Binding` row was added to `aleph/vertical/wiring.py`.** This port adds no connector and
   changes no wiring count. The 30/36 figure is unmoved by it.
6. **`ports/ledger/INDEX.md` not regenerated**, on the standing instruction every lane since G1 has
   followed: the generator reads the working tree, and three other sessions hold uncommitted changes
   right now. `-3626` joins the already-red `test_the_index_is_not_missing_a_tracked_entry` as one
   more *name*, not as a new red test.
