# ALEPH-PORT-3627 — the four connector rows that are neither live nor a cytosol gap

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3627` |
| Lane | `83c83640` (the four leftover connector rows) |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **before the code**, per `CLAUDE.md` §3 and `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` — nothing is copied from `/Users/sw1/ffn_cellsim`; this supplies construction arguments to Aleph's own connector classes |
| Aleph target | `aleph/vertical/connectors_protrusion.py`, `aleph/vertical/connectors_crosslink.py` |
| Depends on | `ALEPH-PORT-3302` (the protrusion connectors), `ALEPH-PORT-3623` (the whole cell's couplings) |
| Exists because | Lane W2's census splits 36 registry rows into 17 live, 9 blocked on the cytosol, 4 that are not gaps, and **4 left over**. These are the four. |
| What this entry does **not** claim | That any coupled number is physics. Every ratchet parameter is `UNVERIFIED` and evidence-owed by its own owner; every fixture geometry is where the constructors put it. |

---

## 1. The four rows, and why "unwired" is the wrong word for three of them

`build_whole_cell()` reports these four in `withdrawn`, quoted verbatim:

| # | Row | Withdrawal reason as recorded |
|---|---|---|
| 1 | `lamellipodium_membrane_contact` | *"BrownianRatchetConnector needs construction arguments this module does not supply (a population, a clutch graph or a solver), so its coupling is left to the lane that owns those objects"* |
| 2 | `filopodium_membrane_tip` | *(identical)* |
| 3 | `if_sf_plectin` | *"only 2 distinct site pairs survived (needed 3): nearest-point put several sites on one partner, or every candidate pair was coincident"* |
| 4 | `ecm_far_field_anchor` | *"endpoint(s) ['world_boundary'] are not owners in this world"* |

**Rows 1–3 are already `Binding` rows in `wiring.py`.** Verified at `9936d20`. They are absent from
the *world*, not from the wiring table, and `whole_cell_couplings.py` states the distinction itself:
*"'Wired' in `wiring.py` means a class resolves and a test exercises it; it has never meant a world
evaluates it."* So this lane adds **no `Binding` row for them** — a row that is already there is not
added twice, and adding one would be the vacuity the wiring guard exists to refuse.

Row 4 carries an `UnwiredReason`, not a `Binding`. This lane's grant is `Binding` rows only.
It is therefore routed to a proposal and **no line of `wiring.py` is written by this lane at all.**

## 2. Rows 1 and 2 — what argument was actually missing

`BrownianRatchetConnector` takes `(name, endpoints, ratchet, barbed_ends_per_site)` and then requires
two pieces of caller state before it will evaluate: `set_axes` and `set_growth_velocity_um_per_s`.
`_require_state` refuses without both, and its refusal message says why the refusal is not optional:

> *"Returning zeros here would be far worse than for a spring: at zero velocity this element is at
> its **maximum** load, so 'nothing was evaluated' and 'stalled' would be the same answer and the
> confusion runs in the dangerous direction."*

`whole_cell_couplings._construct` builds a connector from `(name, endpoints, spring, site_count)` and
`BrownianRatchetConnector` is not in `_SPRING_FAMILIES`, so it returns `None`. **That is correct
behaviour and not a defect**: the generic pairing cannot supply the axes.

### 2.1 Why the axis cannot come from the generic nearest-point pairing

The axis **is the filament polarity**, barbed-end-ward. It is not the line of centres between the
paired sites. `BrownianRatchetConnector.set_axes` says the consequence in one line — *"an axis of
length 0.9 would rescale the delivered load by 10% and read as a soft ratchet"* — and a line-of-centres
direction is not merely mis-scaled, it points somewhere else entirely. Only the protrusion owner
knows the polarity, and it publishes it as `polar_direction(filament_id, s)`.

So the missing object is the same shape as the one `lamellipodium_actin_binding_sites` already
supplies for the motor: **a projection from an owner's published sites into what the connector's
constructor wants.** This lane writes the ratchet's counterpart beside it, in the same module, with
the same delivery discipline — the reaction goes back through `owner.endpoint_sink(site_id)`, with
the owner's interpolation weights, and not through a node index.

### 2.2 The membrane correspondence is *ahead along the axis*, not nearest

Nearest-point can select a membrane vertex **behind** the barbed end. A ratchet that pushes against
a vertex behind its own tip is not a weakly registered ratchet; it is a ratchet delivering its load
into the half-space it cannot reach. `BrownianRatchetConnector.accumulate` already records the
diagnostic that detects it — the per-site moment *"vanish[es] exactly when the membrane contact lies
ahead of the barbed end along the axis"* — so the correspondence this lane writes takes the nearest
vertex with a **positive** axial projection and **refuses** when none exists, rather than falling
back to nearest and leaving the moment to report a registration failure nobody reads.

### 2.3 What is not supplied, and stays not supplied

`kb_t_pn_um`, `monomer_step_um`, `free_velocity_um_per_s`, `reverse_ratio` and the growth velocity
are **required arguments with no defaults**, here as everywhere. `lamellipodium.rate_prior` refuses
`free_barbed_end_elongation_rate_um_per_s` as evidence-owed and this module does not route around its
own owner's refusal. The builder demands them from its caller; whoever passes them owns them, and
they are `UNVERIFIED`.

## 3. Row 3 — the pair collapse is geometry, and the rule change alone makes it worse

Measured before anything was written (§10.1). The `intermediate_filament` fixture has **14** nodes and
`sf_arc` **11**, and the distance matrix between the 8 sampled IF sites and all 11 arc nodes shows the
cause is not the rule:

* the IF network occupies `x ∈ [0, 2.05]`; the arc's bundle 0 places nodes at `x = 0, 2, 4, 6`;
* so **only arc nodes 0 and 1 lie within 2 µm of any IF node**, and columns 2..10 are all ≥ 1.9 µm;
* nearest-point maps all 8 IF sites onto those two, and `np.unique` keeps 2.

**An injective rule alone does not fix this, it launders it.** A greedy minimum-weight matching
recovers 7 pairs — at separations up to 5.6 µm. Plectin is a ~0.2 µm cytolinker. Six of those seven
pairs are crosslinks that cannot exist, and a coupling that carries them reports a stiffness that is
*too high* by a factor nobody can see, which is the same failure as the one being fixed with its sign
reversed.

So this lane writes **both halves**, and neither is silent:

1. **an injective correspondence with a required reach cutoff** (`max_span_um`), in
   `connectors_crosslink.py` — one partner per site by construction, and no pair beyond the linker's
   physical span, with the count rejected on each ground reported separately; and
2. **the measurement of the refinement route** — what the surviving count becomes when the arc is
   discretised finely enough to have nodes where the IF actually is (§10.2).

**The stiffness factor is made visible rather than repaired.** `whole_cell_couplings` uses a per-site
`5.0 pN/µm` fixture, so the population's total is `N × 5.0` and silently depends on how many pairs the
geometry happened to yield. This lane's builder takes the **total** and reports the per-site value it
derived, so the invariant quantity is the declared one and a control can assert that changing `N`
leaves the total where the card put it. §10.3 records the alternative reading — that a population
which genuinely has fewer linkers genuinely is softer — as the thing the PI may overturn. It is not
settled here.

## 4. Row 4 — `world_boundary` has no owner, and that is a decision, not a bug

`ecm_far_field_anchor` names `('ecm', 'world_boundary')`. `world_boundary` is no longer an owner: on
2026-08-04 the PI ruled that the registry lists **engine parts, not cell parts**, and
`aleph/state/full_census.py` records the consequence — the frame moved to `EXPLICIT` and now carries
five state keys (far-field displacement, a strain schedule, a surround compliance, an applied
traction, and the reaction it already received). `V-ENDPOINT-NOT-EXPLICIT` was removed *with* that
decision.

So the endpoint is explicit state with no owner object, and how an external load reaches the ECM is
the real question. It is **not** answered here. See §11 and
`docs/decisions/PROPOSAL-the-far-field-anchor-needs-a-load-carrier.md`, `Status: AGENT-PROPOSED`, with
no `decided_by` field.

## 5. Method — one connector at a time, on PI instruction

The PI specified the order in session: build one connector, call `evaluate_sites`, measure `max|F|`
and the Newton residual `|F_a + F_b|`, then the next. Not build all four and test at the end. The
per-row tables in §10 are in that order and were produced that way.

## 6. Acceptance

| Gate | Where |
|---|---|
| `max|F|` (pN) and Newton residual per row | §10 |
| every deliberate-break flag driven, and the answer measured to move | §10.4 |
| `tests/vertical/test_connectors_protrusion.py`, `test_connectors_crosslink.py` green | §10.5 |
| no `aleph/scenarios/**` byte written | §12 |

## 7. Declared blind spots

Stated here before measuring, so a zero that could not have been anything else is not read as
evidence. `ALEPH-PORT-3604` §9a is the precedent: a Newton-pair control for `ErmTether` *cannot* fail,
because the two sides are the same product with one sign, and it shipped as a declared blind spot with
a control asserting the bit-identity rather than as a passing check.

The ratchet is the same shape. `force_a = -(loads[:, None]) * axes` and `force_b = (loads[:, None]) *
axes`; IEEE-754 gives `(-a) * b == -(a * b)` exactly, so `|F_a + F_b|` is **exactly** `0.0` and no
defect in this connector can make it otherwise. It is reported as `0.0` in §10 **and named here as
uninformative**. The measurements that can fail are the stall-force oracle, the power sign, and the
axial-registration refusal.

## 8. Evidence class

`UNVERIFIED` throughout. No parameter here is sourced; an oracle is not a simulation and a passing
control is not physics.

## 8a. Source identity — nothing was read from the reference tree

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) — **not opened by this lane** |
| Source commit | none |
| Source path | **none.** No `ffn_sim/**` file was read, and no symbol was consulted |
| Source symbol(s) | none |
| Read from | nothing |
| Why not clean-room | it *is* clean-room. `RE-DERIVED`: every object here is constructed from Aleph's own published owner surfaces (`polar_direction`, `free_barbed_ends`, `growing_barbed_end_count`, `endpoint_sink`) and Aleph's own connector classes. The registered mechanism in `aleph/state/**` is the specification; the reference tree contributes nothing to it. `PRIOR_ART_REALIZATION` records `if_sf_plectin` as `SEAM_ONLY` there — a typed slot with no implementation — so there was nothing to port even had this lane wanted to. |

## 8b. Numerical and precision envelope

Everything here is float64 on the CPU backend. **No GPU, no kernel, no device residency** — these are
host-side NumPy constructions feeding connector classes that `ALEPH-PORT-3604` and `-3606` have
already ported to Warp; this entry adds no law and therefore no parity obligation.

| Quantity | Envelope | Basis |
|---|---|---|
| Newton residual `‖F_a + F_b‖` | **exactly `0.0`**, asserted with `==`, not a tolerance | both families form the two sides as `±(scalar × axis)`; `(-a)·b == -(a·b)` bitwise in IEEE-754 |
| `f_stall` under a six-decade `v0` sweep | **0 ULP** — `3.650023777042 pN` in every digit at `v0 = 3e-3 … 3e2` | `v0` cancels out of `(kT/d)·ln(1/e)` analytically |
| bundle law `stall(n) / stall(1)` | exact integers to 12 printed digits (`2.000000000000`, `100.000000000000`) | `n` multiplies a float; no cancellation |
| population total `N·k` | `40.000000000` pN/µm at `N = 1…13`, `rel=1e-14` in the control | one division and one multiplication |
| `two_sided_spring` in compression, break on | **exactly `0.0`** force and energy, asserted with `==` | the unilateral gate returns a structural zero, not a small number |
| axis unit-length | re-normalised in `_ratchet_axes`; `set_axes` refuses above `1e-12` | a non-unit axis rescales every delivered load |
| pair separations | `1e-9` µm coincidence floor (`cytolinker_site_pairs`) | below it `(a−b)/‖a−b‖` is `0/0` and the connectors refuse by name |

**The cancellation this entry does not have.** `ALEPH-PORT-3601` §14.7 records that float32 position
storage costs ~100 ULP on a tether between surfaces 5 µm from the origin. Nothing here stores a
position in float32, and no figure above is a subtraction of two large nearly-equal numbers, so the
position-precision proposal's fork does not reach this entry. Stated so the silence is not read as
having been checked and passed.

## 8c. Comments and docstrings

No prose was carried from anywhere. Both modules are append-only; not one existing docstring was
edited, and the two `__all__` lists plus one new error class are the only pre-existing lines touched.

## 9. Named controls

### 9a. Positive controls — the thing must do what it says

| Control | Location | Asserts |
|---|---|---|
| Positive | `test_connectors_protrusion.py::TestTheRatchetTargetSuppliesWhatTheWorldCouldNot::test_the_lamellipodial_front_is_built_with_the_owner_s_own_polarity` | the delivered axis is `network.polar_direction(...)` **exactly** (`abs=0.0`), not a line of centres |
| Positive | `…::test_the_bundle_tip_carries_the_owner_s_growing_count_and_not_a_caller_s` | the tip's stall is `growing_barbed_end_count ×` a single filament's, at `rel=1e-15` |
| Positive | `…::test_the_two_rows_differ_only_in_where_the_load_lands` | a 4×1 front and a 1×4 tip give the same total load, and the tip's per-site load is 4× the front's |
| Positive | `…::TestTheMembraneContactMustLieAheadOfTheBarbedEnd::test_the_chosen_vertex_has_a_positive_axial_projection` | every coupled membrane vertex is ahead of its own barbed end |
| Positive | `…::TestTheRatchetTargetDeliversItsReactionThroughTheOwner::test_the_reaction_reaches_the_owner_s_nodes` | the scattered reaction reaches the owner's real node array, not a detached buffer |
| Positive | `test_connectors_crosslink.py::TestTheCytolinkerCorrespondenceIsInjectiveAndReachLimited::test_one_partner_serves_at_most_one_site` | the correspondence is injective by construction |
| Positive | `…::test_every_rejected_candidate_is_counted_under_its_own_ground` | `kept + coincident + beyond_reach == requested` — the census is total |
| Positive | `…::TestThePopulationTotalIsTheInvariantQuantity::test_the_total_is_unchanged_by_how_many_pairs_survived` | `N·k == 40.0` pN/µm at `N = 1…13`, `rel=1e-14` |
| Positive | `…::TestTheCrosslinkStillLoadsInCompression::test_the_two_sided_law_carries_load_in_compression` | the cytolinker really does load in compression, so §9b's zero means something |

### 9b. Negative controls — each must fail, and each is driven through the shipped code path

Not hand-edited copies. Every one of these runs the module as it ships, with a declared break flag or
a real geometry, which is the discipline `BrownianRatchetConnector.allow_ratchet_to_pull` exists for
in the first place: *"so the negative control drives the shipped code path rather than a hand-edited
copy."* Measured outcomes are in §10.4.

| Control | Location | Asserts the break changes the answer |
|---|---|---|
| Negative (must fail) | `TestTheDeliberateBreakChangesTheAnswer::test_the_break_flips_the_sign_of_the_load_and_of_the_power` | with `allow_ratchet_to_pull`, the site load, the forbidden-sense total and the active power **all three** change sign or leave zero. If they did not, the healthy `0.0` would measure nothing |
| Negative (must fail) | `TestTheDeliberateBreakChangesTheAnswer::test_the_healthy_branch_refuses_a_super_free_velocity` | with the break **off**, `v > v0` is a `RatchetDomainError` and not a quietly negative load |
| Negative (must fail) | `TestTheCrosslinkStillLoadsInCompression::test_the_break_carries_exactly_nothing_there` | `bidirectional=False` — a cytolinker degraded to a tether — transmits **exactly `0.0`** force and energy in compression, which is two populations passing through each other |
| Negative (must fail) | `TestTheMembraneContactMustLieAheadOfTheBarbedEnd::test_a_vertex_behind_the_tip_is_not_selected_even_when_it_is_nearest` | a counterexample geometry where nearest-point picks the wrong vertex, and the rule does not |
| Negative (must fail) | `…::test_no_vertex_ahead_within_reach_is_a_refusal_and_not_a_fallback` | with only a behind-vertex available, the correspondence raises rather than falling back |
| Negative (must fail) | `TestTheRatchetTargetSuppliesWhatTheWorldCouldNot::test_a_capped_barbed_end_gets_no_ratchet_site` | a fully capped network is refused rather than published as an empty ratchet |
| Negative (must fail) | `TestTheCytolinkerCorrespondenceIsInjectiveAndReachLimited::test_injectivity_alone_manufactures_pairs_the_linker_cannot_span` | **the counterexample this entry turns on**: an injective rule with no reach bound recovers the full count at spans no cytolinker has |

**28 new controls**, 9 positive and 7 negative above plus 12 refusal and determinism checks. Both
files green together at **109 passed** (`0.92 s`).

---

## 10. Measurements

Every figure below is on the whole cell's own fixture at `9936d20`, driven through the public
surface. **Every parameter is `UNVERIFIED`** and none is a measurement of a cell. The ratchet element
is `kT = 4.28e-3 pN·µm`, `d = 2.7e-3 µm`, `v0 = 0.30 µm/s`, `e = 0.10`, giving
`kT/d = 1.585185 pN`, `f_stall = 3.650024 pN`, `w_monomer = 9.855064e-03 pN·µm`.

### 10.0 The table the acceptance criterion asks for

| # | Row | sites × ends | `max‖F‖` (pN) | per-site `‖F_a+F_b‖` | `‖ΣF_a + ΣF_b‖` |
|---|---|---|---|---|---|
| 1 | `lamellipodium_membrane_contact` | 3 × 1 | **9.476824e-01** | 0.0 | 0.0 |
| 2 | `filopodium_membrane_tip` | 1 × 2 | **1.895365e+00** | 0.0 | 0.0 |
| 3 | `if_sf_plectin` | 8 | **6.895451e+00** | 0.0 | 0.0 |
| 4 | `ecm_far_field_anchor` | — | **not built** | — | — |

**Every Newton residual is exactly `0.0` and every one of them is uninformative**, as §7 declared
before measuring. Both connector families assemble their two sides from one scalar block with one
sign — `(-a)·b == -(a·b)` bitwise in IEEE-754 — so no defect in either can make the residual
non-zero. They are reported because the acceptance criterion asks for them, and named as
uninformative because reporting them as passing checks would be a claim above the evidence.
`ALEPH-PORT-3604` §9a is the precedent.

Row 4 has no force column because nothing was built for it; §4 and the proposal say why.

### 10.1 Row 1 — the branched front

3 free barbed ends, 3 sites, front stall **10.950071 pN**. At `v = 0.5 v0`:

| Quantity | Value |
|---|---|
| load per barbed end | `0.94768236` pN at every site |
| site load | `0.94768236` pN |
| active power | `4.264571e-01` pN·µm/s |
| load in the forbidden sense | `0.0` pN |
| stored energy, start / end | `0.0` / `0.0` pN·µm — a ratchet has no potential |
| separations across the contact | `3.005412`, `4.291842`, `3.377035` µm |

**The stall oracle holds bitwise.** `v0` over six decades (`3e-3` … `3e2` µm/s) leaves
`f_stall = 3.650023777042 pN` unchanged in every digit. A spring's blocking force scales with its
stiffness and cannot do this.

**A finding for Lane W2, not a defect here.** At `max_reach_um = 2.0` the correspondence **refuses**:
the whole cell places the lamellipodium's barbed ends **3.0–4.3 µm from the nearest membrane vertex
ahead of them**. A "membrane contact" at that fixture is a contact across 3 µm of cytoplasm. The
figures above are taken at a 20 µm reach so that a number exists to report; the geometry is the
scenario layer's to place and is reported rather than adjusted.

### 10.2 Row 2 — the bundle tip

1 site, 2 growing barbed ends, tip stall **7.300048 pN**, separation `2.003123` µm, active power
`2.843047e-01` pN·µm/s. Load per barbed end `0.94768236` pN — **identical to row 1's**, and the tip
load is `1.895365` pN, exactly twice it. That is the registered distinction between the two rows as
arithmetic: same element, same per-end load, distributed against concentrated.

**The bundle law is exact.** `n` barbed ends stall at `n ×` a single filament's stall:

| `n` | 1 | 2 | 5 | 20 | 100 |
|---|---|---|---|---|---|
| stall (pN) | 3.65002378 | 7.30004755 | 18.25011889 | 73.00047554 | 365.00237770 |
| ratio | 1.000000000000 | 2.000000000000 | 5.000000000000 | 20.000000000000 | 100.000000000000 |

And the other direction of the same curve: a **total** tip load `F` on `n` ends gives exactly the
velocity `F/n` gives one end — `0.144048005`, `0.209827339`, `0.260488850` µm/s at `n = 1, 2, 5`,
agreeing with `v_single(F/n)` to every printed digit.

### 10.3 Row 3 — four routes, and the rule change alone is not the fix

`intermediate_filament` 14 nodes, `sf_arc` 11. Eight pairs wanted.

| Route | Rule | Kept | Separations (µm) |
|---|---|---|---|
| 0 | independent nearest-point + unique — **what runs now** | **2/8** | 0.502, 0.550 |
| A | injective, **no** reach bound | **8/8** | 0.118 … **6.090** |
| B | injective + 0.2 µm reach (plectin) | **1/8** | 0.118 |
| C | refine the arc, then B | see below | |

**Route A is the finding.** It recovers the full count, and **7 of its 8 pairs are longer than a
0.2 µm plectin**, the longest at 6.090 µm. Changing the correspondence alone does not repair the
withdrawal — it inverts it, from a stiffness too low by an invisible factor to one too high by an
invisible factor. The cause is geometric and is stated in the module: the IF network occupies
`x ∈ [0, 2.05]` and the arc's bundle 0 places nodes at `x = 0, 2, 4, 6`, so **only two arc nodes lie
within 2 µm of any IF node at all.**

**Route C — refinement is the route that works, and it is not monotone:**

| bundle-0 nodes | 4 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|
| arc nodes total | 11 | 15 | 23 | 39 | 71 |
| kept at 0.2 µm reach | 1/8 | **0/8** | 3/8 | 4/8 | 4/8 |
| kept at 2.0 µm reach | 3/8 | 5/8 | 7/8 | **8/8** | **8/8** |

The `4 → 8` step *loses* the one surviving pair: at 4 nodes the arc has a node at exactly `x = 2`,
0.118 µm from an IF node, and at 8 nodes the nodes land at `x = 0, 0.857, 1.714, 2.571, …` and that
coincidence is gone. **"Refine until it passes" is not a safe loop** — refinement moves nodes away
from the other population as readily as toward it. Recorded because it is the kind of thing a lane
would otherwise discover by bisecting.

**The coupling, built on route C at 64 arc nodes and a 2.0 µm reach** (8/8 kept, 0.0486…1.5459 µm,
median rest 0.175779 µm), with the population total at `40.0 pN/µm`:

| Quantity | Value |
|---|---|
| per-site `k` derived from the total over 8 sites | `5.000000` pN/µm |
| `max‖F_a‖` = `max‖F_b‖` | `6.895451e+00` pN |
| stored energy | `9.149392e+00` pN·µm |

**The stiffness bookkeeping, and the fork it does not settle.** `distributed_crosslink_spring` holds
the population's total invariant: `N = 2, 3, 5, 8, 13` give `k = 20.0, 13.333333, 8.0, 5.0, 3.076923`
pN/µm and `N·k = 40.000000000` in every row. `whole_cell_couplings` does the opposite — a fixed
per-site `5.0` against a variable `N` — so its `if_sf_plectin` at 2 sites delivers a quarter of what
its 8-site sibling would, and nothing reports the factor.

**Which is correct is a PI decision and is not taken here.** If `N` is a *discretisation* of a
population, holding the total fixed is right. If `N` is a *count of linkers*, a population with fewer
of them genuinely is softer and holding the total fixed invents stiffness. The function makes the
declared quantity the invariant one and says so in its own docstring; what is not defensible is the
current state, where the answer depends on `N` and no caller is told.

### 10.4 Every deliberate-break flag, driven

| # | Break | Healthy | Broken | Did the answer move? |
|---|---|---|---|---|
| 1 | `allow_ratchet_to_pull` | site load `+0.947682` pN, forbidden-sense `0.0`, power `+4.264571e-01` | site load `−0.588997` pN, forbidden-sense `1.766991` pN, power `−7.951460e-01` | **yes, in three channels** |
| 2 | `UnilateralSpring.bidirectional` | `max‖F‖ = 9.250232` pN, `E = 6.845343e+01` pN·µm | **exactly `0.0` / `0.0`** | **yes** |
| 3 | capping, through site publication | 3 sites, front stall `10.950071` pN | 2 sites, front stall `7.300048` pN | **yes, exactly 2/3** |
| 4 | axial registration vs nearest-point | all three projections `+2.9717`, `+2.4440`, `+2.6886` µm | nearest-point gives `+2.9717`, **`−0.2121`**, `+2.6886` µm | **yes** |

Break 1 also confirms the refusal is real: at `v = 1.5 v0` with the flag **off**, the shipped path
raises `RatchetDomainError`; with it **on**, the same velocity returns a negative load and a negative
power — the element absorbing work rather than doing it. The healthy `0.0` in
`load_carried_in_forbidden_sense_pn` therefore measures something, which is the whole point.

**Break 4 is the strongest result in this entry.** It is not a planted mutant: plain nearest-point,
on the whole cell's own fixture, puts **1 of 3 lamellipodial sites behind its own tip**. The
correspondence a generic world would have used mis-registers a third of this front, and every force
it produces is finite, momentum-conserving and closure-passing.

### 10.5 Controls, and four foreign red gates reported rather than fixed

`tests/vertical/test_connectors_protrusion.py` + `tests/vertical/test_connectors_crosslink.py`:
**109 passed** (`0.61 s`), 28 of them new. No pre-existing control was edited. The six vertical files
most exposed to an append — `test_package_exports`, `test_wiring_reasons`, `test_vertical_controls`,
`test_connectors_frame`, `test_lamellipodium_controls`, `test_filopodium_controls` — pass at **247**.

`ruff` at the project's `line-length = 100` is clean on everything this lane wrote. Three findings
remain in the two edited modules and **all three pre-date this lane**, verified against
`git show HEAD:<path>`: `connectors_crosslink.py:6` (E501) and `:33` (`F401 dataclasses.field`), and
`test_connectors_crosslink.py:105` (E501). Not fixed — they are another lane's lines. One baseline
finding *was* removed as a side effect: `connectors_crosslink.py` imported `numpy` without using it,
and this lane's code uses it.

**`tests/ports` has four red gates and every one of them is foreign.** Attributed by name rather
than assumed:

| Gate | Names | Whose |
|---|---|---|
| `test_ledger_entries_are_complete` | `ALEPH-PORT-3625`, `ALEPH-PORT-3628` — both **untracked** in the working tree | two other live lanes |
| `test_the_index_is_not_missing_a_tracked_entry` | `ALEPH-PORT-3626-filament-axial-and-bending-kernels.md`, **committed** and absent from `INDEX.md` | another lane |
| `test_index_is_not_stale` | `INDEX.md` | **verified foreign**: still red with this entry moved out of `ports/ledger/` and the test re-run |
| `test_no_provider_vocabulary_leaks_into_the_package` | `aleph/harness/ingest/adapters.py:147` | session `7bc39358`, already recorded as its known failure in `docs/ACTIVE_SESSIONS.md` |

**`INDEX.md` was not regenerated**, on the precedent this file records five times: the generator reads
the working tree, and `-3625` and `-3628` are uncommitted right now, so regenerating would bake two
other lanes' mid-write ledgers into a committed index. `-3627` joins that already-red assertion as a
name, not as a new red test.

**`-3625`, `-3626` and `-3628` all appeared while this lane was working.** `-3627` was verified free
at 20:1x against both `ports/ledger/` and a tree-wide grep, and it is still free; the neighbours
moving is why the check was done by grep rather than by reading the directory listing alone.

**`tests/vertical` as a whole directory: 11 red, and every one of them is foreign.** A partial run
earlier in this lane saw only two of them and that count was wrong; the full run is the number below.
Recorded rather than quietly replaced, because the under-count is exactly the failure mode
`CLAUDE.md` §1 warns about — a directory being written by four lanes at once does not give a stable
answer, and an early partial read of it invites a confident wrong attribution.

| Failing | Count | Owning lane, and how it was attributed |
|---|---|---|
| `test_cytosol_controls.py::TestTheSolveIsMatrixFreeAndOnTheSubstrate` | 3 | **`5ad3a7ba` (S2)** — `ALEPH-PORT-3625`'s own header names `aleph/vertical/cytosol.py` and this test file as its targets; both are `M` (446 and 240 changed lines) |
| `test_material_point_batch.py::…[sf_arc]` | 7 | a lane mid-write in `aleph/vertical/sf_arc.py` (154 changed lines, `M`). **The test file itself is untracked (`??`) and did not exist in `tests/vertical/` when this lane started** |
| `test_ecm_network_controls.py::TestTheDeclarationsAgreeWithTheCode::test_the_export_list_is_sorted` | 1 | a lane mid-write in `aleph/vertical/ecm_network.py` (161 changed lines, `M`). Verified directly: `ecm_network.__all__` is unsorted at `ECMAttachmentSite, ECMEndpointRole, ECMRoleError, ECMSiteType` right now |

**None of the three failing test files imports `connectors_crosslink` or `connectors_protrusion`**,
verified by grep, and this lane wrote no byte of `cytosol.py`, `sf_arc.py` or `ecm_network.py`. This
is the case `CLAUDE.md` §1 names in advance — *"a whole-repo pytest run can fail on another session's
file that is mid-write … report a foreign breakage rather than fixing it"* — and it is why the
acceptance figure quoted above is the two owned files plus the six most exposed neighbours, and not a
whole-directory number taken while four other lanes are writing into it.

## 11. Raised and not taken

1. **`ecm_far_field_anchor` (row 4).** Three landed artefacts disagree about whether it is a
   connector, and `ecm_network.assert_far_field_frame_owns_no_state` now refuses on a premise the
   2026-08-04 census contradicts. Full statement in
   `docs/decisions/PROPOSAL-the-far-field-anchor-is-a-clamp-and-the-registry-calls-it-a-connector.md`.
   No guard edited, no `aleph/state/**` byte written.
2. **The site-count/stiffness fork** in §10.3. Recorded, not settled.
3. **The whole cell places both protrusions 2–4 µm from the membrane** (§10.1, §10.2). Neither
   `lamellipodium_membrane_contact` nor `filopodium_membrane_tip` can be placed at a physical contact
   distance in that fixture. The geometry belongs to `aleph/scenarios/whole_cell.py`, which is Lane
   W2's; reported, not adjusted.
4. **`whole_cell_couplings._construct` returning `None` is correct and should stay.** The fix is not
   to teach it about ratchets — it cannot know a polarity. It should call the builders in §2 instead.

## 12. What this lane did not touch

`aleph/scenarios/**` (Lane W2 is live on it — `whole_cell.py` and `whole_cell_couplings.py` were
read and driven, never written). `aleph/vertical/wiring.py` — **zero bytes**; rows 1–3 were already
`Binding` rows and row 4 carries an `UnwiredReason`, which is outside this lane's grant.
`aleph/vertical/connectors_frame.py` (`build_if_sf_plectin` lives there and is another lane's),
`aleph/vertical/ecm_network.py`, `aleph/vertical/cytosol.py` (S2), `aleph/runtime/**` (S2/S3),
`aleph/state/**`, `aleph/represent/**`, `aleph/learn/**`, every guard and firewall test,
`ports/ledger/INDEX.md` (already red from other lanes; **not regenerated**). No GPU — every figure
here is force arithmetic on CPU. No `ALEPH-PD-*` file, no `decided_by` field anywhere.

Both edited modules are **append-only**: no existing line of `connectors_protrusion.py` or
`connectors_crosslink.py` was changed except the two `__all__` lists and one added error class.
