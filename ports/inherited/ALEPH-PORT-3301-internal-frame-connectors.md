# ALEPH-PORT-3301 — the five internal-frame connectors: two linkers, one stiffening cable, two cytolinkers

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3301` |
| Lane | connector lane under lead session `f70de564` — internal frame |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Owners this depends on | `ALEPH-PORT-3202` (microtubule), `ALEPH-PORT-3203` (intermediate filament) |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
# new module aleph/vertical/connectors_frame.py
BidirectionalLinker            # a two-sided molecular linker, compression-capable
StrainStiffeningCable          # the constitutive element of if_nucleus_linc — NOT a new law
StrainStiffeningCableLink      # the connector that carries it
build_mt_nucleus_linc
build_mt_cortex_capture
build_if_nucleus_linc
build_if_sf_plectin
build_mt_sf_spectraplakin
```

Test file: `tests/vertical/test_connectors_frame.py`.

**Two classes are reused rather than added.** `if_sf_plectin` and `mt_sf_spectraplakin` are built as
`aleph.vertical.connectors_crosslink.CrosslinkConnector`; §4.4 gives the registry sentences that
make that the right element and not merely a convenient one. Their builders live here because their
*contracts* are internal-frame contracts, not because they need a class of their own.

**No `Binding` row is added by this lane.** `aleph/vertical/wiring.py` has a single writer. The rows
this entry has earned are quoted in §13 for that writer to place.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | **none named and none read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened for this work. |
| Working tree == commit? | not applicable — nothing was read |

The input to this work is `aleph/state/connectors_crosssystem_cytosol.py` lines 447–550 — five
registered `ConnectorContract`s — plus the two owner modules `ALEPH-PORT-3202` and `-3203` already
landed. There is no third artefact and none was wanted.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** Nothing here is a new force law at all. Four of the five
connectors carry a law that already exists in this repository with its own gradient control
(`UnilateralSpring` with `bidirectional=True`, via
`aleph.vertical.connectors_crosslink.two_sided_spring`); the fifth carries the three-branch cable law
that `aleph/vertical/intermediate_filament.py` already implements and that `ALEPH-PORT-3203` already
audited. A port would have had to introduce a *second* copy of one of those, which is the failure
`tests/vertical/test_connectors_crosslink.py::test_it_agrees_with_the_owner_s_own_crosslink_gradient`
was written to make impossible.

## 4. Physical or mathematical law represented

**The claim in one line: a connector adds an interface, not a law — and the element type is read off
the registered mechanism rather than chosen.**

The five contracts differ in three ways that are mechanical and not nominal: whether the coupling can
carry **compression**, whether its stiffness is **constant**, and what its two endpoint **roles** are.
Read in that order the registry decides all five.

### 4.1 `mt_nucleus_linc` — two-sided, and it is the only one the registry says so about

> "A microtubule can also push, which neither of the other two LINC edges does, so this is the only
> LINC path that can deliver a compressive nuclear load."

That sentence has no implementation if the element is tension-only. `actin_cap_linc` is a
`connectors_solid.TensileLinker` precisely because the registry gives *it* no such sentence; using
the same class here would silently delete the one property that makes this edge distinct, and the
deletion would be invisible — a tension-only branch satisfies `F = -grad E` exactly, closes force and
moment exactly, and conserves energy exactly, because a missing branch conserves energy perfectly.
`aleph/vertical/microtubule.py` records the same argument one level down: its axial law is two-sided
and its docstring cites *this* contract as the reason.

So the element is the two-sided spring `U(s) = (k/2)(s - s0)^2`, on every branch:

```
U(s)      = (k/2) (s - s0)^2                 for all s > 0
dU/ds     = k (s - s0)
F_s       = -dU/ds                            positive pushes the two sites apart
```

built through `two_sided_spring`, which is `UnilateralSpring(bidirectional=True, law=LINEAR)` — the
one implementation of that quadratic the wired vertical already uses. Nothing is reimplemented.

**Why a new class and not `CrosslinkConnector`.** The mechanics are identical and the class would be
four lines shorter. The objection is the name: `CrosslinkConnector`'s docstring says "a population of
two-sided **crosslinks**", and a LINC complex is a linker between two components, not a crosslink
within a network. Worse for `mt_cortex_capture`, whose contract distinguishes it *from* a transient
crosslink in so many words (§4.2). A wiring table that binds two contracts to a class whose name
denies them is a table a reader has to un-learn. `BidirectionalLinker` is therefore
`connectors_solid.TensileLinker`'s two-sided counterpart, and the only behaviour it adds over the
base class is the same `load_carried_in_forbidden_sense_pn -> 0.0` override `CrosslinkConnector`
carries, for the same reason: an element with no forbidden sense has no branch for load to leak into,
so the accessor is structurally zero rather than a measurement, and saying so is the honest report.

### 4.2 `mt_cortex_capture` — two-sided, on the registry's own distinguishing sentence

> "The cortex side is a capture site rather than an arbitrary material point, which is what
> distinguishes this edge from a transient crosslink."

The registry locates the difference from `sf_cortex_transient` in the **site identity**, not in the
force law. Had the law been the distinguishing feature — a capture that only pulls — that is what the
sentence would have said, because a one-sided law is a much larger difference than a labelled site.
Read with "a captured plus end is **held** at the cortex" (held is a positional constraint, not a
tether) and with the microtubule being the compression-capable member of the frame, the element is
the same two-sided spring as §4.1.

**This is a reading and it is recorded as one.** The registry does not state the sense of the capture
bond in words of its own. If the reading is wrong the repair is one builder returning a
`TensileLinker`; no law and no test scaffolding moves. §14 carries it as an open question rather than
as a settled fact.

Also, and separately: this connector is the microtubule's **only** declared path to the cortex — the
component has exactly four registered connectors and no steric one. A tension-only capture would
therefore leave a growing rod nothing at all to push the cortex with. That is an argument for
two-sidedness and not a proof, because capture acts only where a bond exists and steric contact acts
everywhere; the missing MT–cortex steric edge is a registry gap either way and §14 names it.

### 4.3 `if_nucleus_linc` — tension-only, and **not** linear

> "The intermediate-filament network is the strain-stiffening cage around the nucleus, so this edge
> dominates the nuclear load at large deformation while contributing little near rest — a load
> character neither the actin cap nor the microtubule edge reproduces."

A linear tether cannot produce that sentence: with a constant stiffness this edge's share of the
nuclear load is the same fraction at every deformation, so "dominates at large deformation while
contributing little near rest" would have no implementation. The load character *is* the three-branch
card law, and that law already exists, audited, in the owner:

```
e(s)  = (s - s0) / s0                                    strain on the bond's rest length
T(e)  = intermediate_filament.tension_pn(e, card)        exactly 0.0 for e <= 0
U(s)  = s0 * intermediate_filament.storage_energy_density_pn(e, card)
dU/ds = s0 * W'(e) * (1/s0) = T(e)                       the tension IS dU/ds
```

`StrainStiffeningCable` is that identity and nothing else: it holds an `IFMaterialCard` and forwards
to the owner's two scalar oracles. It defines **no branch logic, no knee, no plateau and no
integral** — every one of those lives in `intermediate_filament.py` and the positive control pins the
connector's site force against the owner's own `cable_axial_energy_and_forces` kernel at matched
parameters.

It is a *duck-typed stand-in* for `UnilateralSpring`, carrying the five members
`_DistributedUnilateralConnector` and `site_pair_forces` actually read — `energy_pn_um`,
`derivative_pn`, `scalar_force_pn`, `is_engaged`, `rest_separation_um`, `engages_in_tension` — so the
distributed machinery, the `bound` gate, the two independent directed accumulations, the census and
the transaction triple are all inherited unchanged. `scalar_force_pn` is not reimplemented: it calls
`UnilateralSpring.scalar_force_pn(self, ...)`, which touches only `self.energy_pn_um` and
`self.derivative_pn`, so the discrete-gradient step-average rule stays in exactly one place.

Tension-only: `T(e) = 0` identically for `e <= 0`, so the inactive branch is testable with `==`.

### 4.4 `if_sf_plectin` and `mt_sf_spectraplakin` — two-sided cytolinkers, `CrosslinkConnector`

Plectin is a cytolinker, and `intermediate_filament.crosslink_energy_and_forces` — the owner's own
cytolinker — already states the law and the reason: "**Two-sided, and that is not an oversight.** The
tension-only discipline belongs to the cables… applied to a cytolinker it would let two crosslinked
cables pass through each other." `if_sf_plectin` is that element with its far end on `sf_arc` instead
of on a second cable, which is exactly the step `sf_cortex_transient` already took, so
`CrosslinkConnector` is reused rather than duplicated.

Spectraplakin gets it from the registry directly:

> "It is the load path that lets actin contractility bend and buckle microtubules and, **in the other
> direction**, lets microtubule compression resist actin contraction."

Both directions named in one sentence is a two-sided element. Its microtubule side is a **lattice
site** and not a plus end — "its position along the lattice is part of the kinetic state" — which
`microtubule.ROLE_SITE_KINDS` already enforces at the owner; the connector inherits that constraint
rather than restating it, and a control asserts the owner still refuses a plus-end spectraplakin site
so the constraint is known to be live.

### 4.5 What every one of the five inherits, and does not restate

Two independent directed accumulations (`force_b := -force_a` appears nowhere), the discrete-gradient
step average, the `bound` gate that makes an unbound site carry an exact `0.0` at any separation, the
`SiteField` endpoint shape, `snapshot`/`commit`/`rollback` over `bound`, `owned_entities`, and the
`AdjointPair` that `accumulate` returns with its constituent force and moment scales. All of it from
`_DistributedUnilateralConnector`.

## 5. Units, domains, singular cases, invariants

Length µm, force pN, energy pN·µm, stiffness pN/µm. `IFMaterialCard`'s five moduli are pN — tension
at unit strain — so `U = s0 * W(e)` is µm·pN and `dU/ds = T` is pN, with no conversion anywhere.

| Singular case | Behaviour | Why not softened |
|---|---|---|
| Coincident site pair | `ConnectorGeometryError` (inherited) | A regularised denominator returns a large force in a direction taken from round-off. |
| `site_count <= 0` | `ValueError` at construction | A population that transmits nothing should not be registered, not registered empty. |
| `rest_length_um <= 0` | `ValueError` | The strain in §4.3 divides by it. |
| Negative or non-finite stiffness | `ValueError` (inherited from `UnilateralSpring`) | — |
| A card that does not stiffen | `ValueError` from `IFMaterialCard` (the owner's) | `ALEPH-PORT-3203` §4; not restated here. |
| Unbound site, any separation | exactly `0.0`, the number | Rupture is state, not a threshold. |
| `if_nucleus_linc` at or below rest | exactly `0.0` energy and force | A cage cable does not push. |

Invariants asserted by control rather than by comment:

* **I-a** Newton pair closes per site: `force_a + force_b == 0` to 1e-14, from two independently
  normalised direction vectors.
* **I-b** The connector's site force equals the *endpoint owner's own* gradient at matched
  parameters — `microtubule.axial_energy_and_forces` for the three MT edges,
  `intermediate_filament.cable_axial_energy_and_forces` for `if_nucleus_linc`,
  `intermediate_filament.crosslink_energy_and_forces` for `if_sf_plectin`.
* **I-c** A site force pushed through the owner's real `MaterialPointSink` equals minus the finite
  difference of the connector's energy **in node coordinates**. This is the chain-rule closure both
  owners' `MaterialPoint` docstrings claim, checked across the seam rather than inside one module.
* **I-d** `sum f_i` and `sum (x_i - c) × f_i` are judged against the **constituent** scales
  `sum|f_i|` and `sum|r_i||f_i|`, never against the resultant.
* **I-e** Every connector's `name` and `endpoints` equal the registered contract's.

## 6. Source evidence class and known retractions

Evidence rung: `ANALYTIC_ORACLE`. Quantitative status: `BLOCKED`.

**No magnitude here is a measurement.** Every stiffness, rest length and card in the module's
builders is supplied by the caller and every number in the test file was chosen to make the algebra
sharp. Nothing in this entry is evidence that a LINC complex has this stiffness, that a plectin bond
has this rest length, or that the IF cage's knee sits where the fixture card puts it — the fixture
cards are unsourced (`citation=None`) and the owner's `IFResponseClaim` machinery already refuses to
call a response from an unsourced card strain-stiffening.

No retraction is known against the two owner entries this depends on.

## 7. Independent oracle or derivation

Four independent checks, no two of which share an implementation:

1. **The owner's own gradient kernel**, per §5 I-b. `microtubule.axial_energy_and_forces` computes
   `E = (k_ax/2)(L - L0)^2 / L0` with `T = k_ax (L - L0)/L0`; set `k_ax = k s0` and `L0 = s0` and it
   is algebraically the connector's two-sided spring, reached through a completely different code
   path (a segment table and a scatter, not a site-pair normalisation). Agreement to 1e-12 relative
   is therefore a real cross-check and not a tautology.
2. **Central finite difference in node coordinates**, per §5 I-c, at steps `(4e-4, 2e-4, 1e-4)` µm,
   asserting observed order > 1.6 with `energy > 0` asserted first — without that, a term returning
   zero energy and zero force agrees with its own finite difference perfectly.
3. **A closed form for the three-branch law**: at a strain placed in a named branch, `T` and `U` are
   written out by hand from the card's five numbers and compared with the connector's.
4. **The ratio that is the registry's actual claim** for `if_nucleus_linc`: the ratio of this edge's
   tension to a linear tether's, matched at a small strain, must *grow* with strain. That is
   "dominates at large deformation while contributing little near rest" stated as a number.

## 8. Positive control

`tests/vertical/test_connectors_frame.py`. The load-bearing ones, in the order they matter:

* `test_the_mt_linkers_agree_with_the_microtubule_s_own_axial_gradient` — positive control for I-b on
  the three microtubule edges, **including in compression**, which is where the contract's claim lives.
* `test_the_if_linc_agrees_with_the_owner_s_own_cable_gradient` — positive control for I-b on
  `if_nucleus_linc`, against `cable_axial_energy_and_forces` in all three branches.
* `test_the_plectin_coupling_agrees_with_the_owner_s_own_crosslink_gradient` — positive control for
  I-b on `if_sf_plectin`, against the owner's cytolinker kernel summed per side.
* `test_a_scattered_site_force_is_minus_the_gradient_in_node_coordinates` — positive control for I-c,
  through the real `MaterialPointSink` of both owners.
* `test_the_if_linc_share_of_the_load_grows_with_deformation` — positive control for oracle 4.
* `test_every_pair_force_is_equal_and_opposite`, `test_each_names_a_declared_contract_with_the_declared_endpoints`,
  `test_an_unbound_site_carries_exactly_zero_whatever_the_separation`,
  `test_the_transaction_triple_restores_the_bond_state`,
  `test_each_one_carries_accumulate_and_not_only_evaluate_sites`.

## 9. Deliberately failing negative control

Each builds the wrong element on purpose and asserts in **both** directions — the broken variant
misbehaves *and* the shipped one is exactly inert (`== 0.0`, the number).

* `test_a_tension_only_linc_would_miss_the_compressive_nuclear_load` — the counterexample for §4.1.
  Builds a `TensileLinker` on the same fixture and shows it stores exactly zero where the shipped
  connector pushes. Without it, "two-sided" is a word.
* `test_a_linear_spring_cannot_reproduce_the_strain_stiffening_load_character` — the counterexample
  for §4.3 and for oracle 4: a linear tether matched at small strain keeps a *constant* share of the
  load, so a `StrainStiffeningCableLink` silently degraded into a linear one fails.
* `test_a_cable_allowed_to_compress_is_caught_pushing_on_the_nucleus` — drives the owner's shipped
  `allow_cable_compression` break through the connector, so the tension-only claim has something real
  that can fail it.
* `test_a_one_sided_capture_is_distinguishable_from_the_shipped_one` — must fail if §4.2's reading is
  ever changed silently.
* `test_the_owner_still_refuses_a_plus_end_spectraplakin_site` — the lattice-site constraint of §4.4,
  shown to be live at the owner rather than asserted here.

## 4a. Mutation testing

Three mutants applied to `aleph/vertical/connectors_frame.py`, the controls run against each, then
restored and verified byte-identically. Results are in §13.

## 10. Numerical and precision envelope

float64 throughout. The gradient agreements are asserted at `rel=1e-12`; the finite-difference
control at observed order > 1.6 over `(4e-4, 2e-4, 1e-4)` µm, which is the same ladder the two owner
control files use, chosen so truncation still dominates round-off at the fixture's energy scale. The
exact-zero assertions use `== 0.0` and not a tolerance, because the inactive branch is identically
zero over its whole extent rather than small near the gate.

The discrete-gradient floor `_DISCRETE_GRADIENT_FLOOR_UM = 1e-13` is inherited from
`connectors.py` and is not re-tuned here; the strain-stiffening element crosses it on the same rule as
every other, since it borrows `UnilateralSpring.scalar_force_pn` verbatim.

## 11. Production-backend residency and transfer

CPU/NumPy only. `accumulate` scatters through `ctx.backend.scatter_add`, which is the one backend
call in the inherited path, so a future `WarpBackend` moves these five connectors without a change
here. `warp` is not installed and no GPU work was run for this entry.

## 12. Comments and docstrings to discard

**None, because nothing was read.** No comment, identifier, or number in `connectors_frame.py`
originates outside this repository. The reference project is named in this entry, which is required,
and appears nowhere in `aleph/**`, which is enforced by
`tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package`.

## 13. Acceptance

Status stays `PROPOSED` until a `Binding` row lands for each connector.

The rows this entry has earned, for `wiring.py`'s single writer:

```python
Binding(contract="mt_nucleus_linc",
        dotted_path="aleph.vertical.connectors_frame:BidirectionalLinker",
        evidence="tests/vertical/test_connectors_frame.py"),
Binding(contract="mt_cortex_capture",
        dotted_path="aleph.vertical.connectors_frame:BidirectionalLinker",
        evidence="tests/vertical/test_connectors_frame.py"),
Binding(contract="if_nucleus_linc",
        dotted_path="aleph.vertical.connectors_frame:StrainStiffeningCableLink",
        evidence="tests/vertical/test_connectors_frame.py"),
Binding(contract="if_sf_plectin",
        dotted_path="aleph.vertical.connectors_crosslink:CrosslinkConnector",
        evidence="tests/vertical/test_connectors_frame.py"),
Binding(contract="mt_sf_spectraplakin",
        dotted_path="aleph.vertical.connectors_crosslink:CrosslinkConnector",
        evidence="tests/vertical/test_connectors_frame.py"),
```

Measured on completion, `2026-07-31`, with
`/Users/sw1/miniconda3/envs/aleph/bin/python -m pytest … -p no:cacheprovider`:

| | |
|---|---|
| `tests/vertical/test_connectors_frame.py` | **29 passed** |
| `tests/ports/` | 15 passed, 2 failed — **both foreign**, see below |
| `ruff check` on both new files | clean |
| Mutants applied | 3 |
| Mutants that survived | **0** |

| # | Mutant on `connectors_frame.py` | Killed | Tests that caught it |
|---|---|---|---|
| M1 | `build_mt_nucleus_linc` uses a tension-only `UnilateralSpring` instead of `two_sided_spring` | 4 | `..._agree_with_the_microtubule_s_own_axial_gradient`, `..._carry_load_in_both_senses`, `..._report_no_forbidden_sense`, `test_a_tension_only_linc_would_miss_the_compressive_nuclear_load` |
| M2 | `StrainStiffeningCable.energy_pn_um` drops the `rest_separation_um` factor, so `U != s0 W(e)` | 2 | `..._agrees_with_the_owner_s_own_cable_gradient`, `test_a_scattered_site_force_is_minus_the_gradient_in_node_coordinates` |
| M3 | `build_if_sf_plectin` declares its endpoints in the reverse order | 1 | `test_each_names_a_declared_contract_with_the_declared_endpoints` |

M3 is deliberately a *weak* mutant and its single kill is the correct result, not a thin suite: an
endpoint-order swap has no mechanical consequence at all — the pair force is the same — so the only
thing that can catch it is agreement with the registry, and that is exactly the test that fired.
Recorded rather than replaced with a mutant that scores better.

Restored and verified byte-identically: the file's SHA-256 after the third restore equals the
pre-mutation one (`1fccdfc84e0d…`), and the suite is back to 29 passed.

**The two `tests/ports/` failures are foreign to this lane and are reported rather than fixed**, per
`CLAUDE.md` §1:

* `test_named_controls_resolve_to_real_tests` names `ALEPH-PORT-3302` (16 controls) and
  `ALEPH-PORT-3304` (14 controls), both written by sibling lanes at 00:21 and both still
  ledger-before-code. `ALEPH-PORT-3301` is not implicated: its 8 named controls all resolve, checked
  directly against `collect_defined_test_functions()`.
* `test_index_is_not_stale` — `ports/ledger/INDEX.md` is generated and four new entries (3301–3304)
  landed in the same window. Regenerating it from this lane would sweep three other sessions'
  entries into a file none of us owns, so it is left for whoever commits.
  `python ports/regenerate_index.py` is the fix, run once after the lanes close.

## 14. Honest limits — what this entry does NOT establish

1. **No kinetics, on any of the five.** All five contracts are family `[K]`: bind, unbind, capture and
   release are accepted-step events. This lane ships the `bound` flag the element already has and
   **nothing that flips it**. So nothing here is evidence about LINC lifetime, about the dwell time of
   a cortical capture, about plectin or spectraplakin turnover, or about force-dependent rupture on
   any of the five. Absent, not implied.
2. **The kinetic state the registry puts in these contracts is not carried.** `mt_nucleus_linc` and
   `mt_cortex_capture` declare "the identity of the bound lattice site" as accepted-step state, and
   `mt_sf_spectraplakin` declares "its position along the lattice" as kinetic state. A
   `_DistributedUnilateralConnector` holds a fixed site list and a boolean per site; it has no
   representation of a bond that *migrates* along the lattice. That is a real gap between the contract
   and the element, and it is a gap in the connector protocol rather than in either owner — the
   microtubule owner already resolves a lattice site by material coordinate and would answer a moving
   coordinate correctly if anything asked it to.
3. **`mt_cortex_capture`'s two-sidedness is a reading, not a quotation** (§4.2). The registry never
   states the sense of the capture bond. This is the one element choice in this entry that is not
   forced by a sentence, and it is flagged for the PI rather than buried.
4. **There is no declared microtubule–cortex steric contact.** The microtubule has four registered
   connectors and none of them is a contact, so a rod can pass through the cortex anywhere it is not
   captured, whatever `mt_cortex_capture` does. That is a registry-level absence and this lane does
   not invent an edge to cover it.
5. **`mt_cortex_capture` transmits no motor force.** Cortical dynein is a motor and the registry says
   so ("a dynein site", "the motor's bound lattice site"). Family `[K]` is what is implemented here —
   a passive elastic bond whose attachment is discrete — and no active work, no stall force and no
   velocity dependence is present. `AdjointPair.active_work` is `0.0` and that is a true report of
   this element, not of a dynein.
6. **The three-branch card is unsourced.** `IFMaterialCard(citation=None)` is what the fixtures use,
   so by the owner's own `IFResponseClaim` rule the only modulus reportable from this connector is the
   small-strain tangent. `test_the_if_linc_share_of_the_load_grows_with_deformation` establishes the
   *shape* of the load character the registry describes; it establishes no magnitude.
7. **Nothing has been stepped inside a world.** These five connectors have controls; they have never
   been evaluated inside `assembly.py` against a live microtubule network and a live nucleus. A
   connector that agrees with its owner's gradient at a fixed configuration is not a relaxation.
8. **The nucleus, cortex and sf_arc sides are not exercised as owners.** I-c is checked through the
   real `MaterialPointSink` of `microtubule` and `intermediate_filament`, which are this lane's two
   owners. The far endpoint of each connector is exercised as a `SiteField` of bare positions, so
   nothing here tests that `nucleus_lamina`, `cortex_filaments` or `sf_arc` scatter these loads
   correctly.
9. **`central_force=True` is reported by the inherited `accumulate` and is true for all five**, since
   every element here acts along the line of centres. It is inherited, not re-verified in this lane's
   controls beyond the moment closure.
