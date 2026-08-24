# ALEPH-PORT-3623 — the couplings the whole cell was missing

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3623` |
| Lane | `46143f30` Lane W (the whole cell, for ablation) |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **after the code, which is a violation, recorded in §0** |
| Port class | `RE-DERIVED` — nothing is copied; this instantiates Aleph's own connector classes in Aleph's own world |
| Aleph target | `aleph/scenarios/whole_cell_couplings.py`, `aleph/scenarios/whole_cell.py` |
| Depends on | `ALEPH-PORT-3621` (the whole cell) |
| Followed by | `ALEPH-PORT-3624` (the adhesion rows, which needed a different shape entirely) |
| Exists because | `build_whole_cell` assembled thirteen owners carrying **two** couplings. Ten compartments sat in one transaction exchanging force with nothing. |
| What this entry does **not** claim | That any coupled number is physics. Every fixture parameter is `UNVERIFIED` and the geometry is where the constructors put it. |

---

## 0. This entry is late, and that is the first thing in it

`CLAUDE.md` §3 and `PLAN.md` §0.2.5 both say the ledger is written **before** the code. This one was
not: `aleph/scenarios/whole_cell_couplings.py` was written, tested and committed with a docstring
citing `ALEPH-PORT-3623`, and the file that citation named did not exist. It was found by reading the
ledger directory while writing `ALEPH-PORT-3624`, not by a guard.

Recorded rather than quietly backfilled, because the useful part is *how it was found*: a module
citing a ledger is not currently checked against the ledger directory. The port-discipline guard
checks citations **inside** ledgers. A module that cites a ledger nobody wrote is exactly the shape
of a claim above its evidence, and it survived a full-repo test run. **Proposed as a guard**, not
written by this lane: see §9.

## 1. What was actually missing, and it was not physics

Every connector class the registry binds was already written, tested and green. What nobody had
written is the object that turns a class into a coupling in a world:

```python
ConnectorSlot(connector, resolver(owner_a, owner_b), owners=(owner_a, owner_b))
```

`build_vertical` constructs exactly two connector objects — the membrane/cortex contact and the ERM
tether — so those were the two the whole cell could carry. The other twenty-eight existed as classes
with no instance.

**"Wired" in `wiring.py` means a class resolves and a test exercises it. It has never meant a world
evaluates it**, and `scenarios/audit.py` says so in as many words: *schedulable is not correct and it
is not exercised either*. The gap between 30/36 "wired" and 2 couplings in a world is that sentence.

## 2. The one piece that is genuinely per-pair

A connector is handed **paired site arrays** — `evaluate_sites(positions_a, positions_b,
triangles_b)` with `positions_a[i]` coupled to `positions_b[i]`. Producing those pairs is the
correspondence, and it is geometry rather than boilerplate. `build_vertical` gets its pair for free
because `build_radial_cortex_network` builds the cortex *outward from* the membrane and returns the
map as a by-product. **No other pair has that.**

`site_pairs` is the general rule — nearest point — with the two cases it must not get wrong.

### 2.1 A self-coupling must not pair a site with itself

`ecm_crosslink` and `dorsal_arc_crosslink` name one owner on both sides. Nearest-point with the owner
against itself returns the identity, every separation is exactly `0.0`, and a central force has no
direction there: `(a−b)/|a−b|` is `0/0`.

Measured, not reasoned: the connectors refuse it by name —
`ConnectorGeometryError: a coupled site pair is coincident at one end of the step`. The guard was
right and the correspondence was wrong.

Both rows are now reported ABSENT with the reason that they are **internal forces, not couplings**:
`AdjointPair` refuses a connector naming one owner twice. The registry row is real; its home is the
owner's internal law.

### 2.2 A degenerate pair on a small owner must not be smuggled through

`if_sf_plectin` and `mt_sf_spectraplakin` join owners with 11–16 nodes, where nearest-point collapses
several sites onto one partner. Those are dropped **and reported**, not silently deduplicated: *a
coupling that quietly carries three sites where its card says eight has a stiffness that is wrong by
a factor nobody can see.*

## 3. The routing adapter, and the incident behind the guard that demanded it

`CellWorld.evaluate_forces` refuses a force block whose shape does not match the owner's array —
and refuses it **rather than dropping it**:

> This connector does not deliver to its declared endpoints — a routing adapter is needed, and
> **dropping the block instead would remove the coupling while the world still relaxed**.

`world.py`'s comment above that check names the incident it was written for: *a spreading world once
lost every adhesion to exactly that silent drop, and relaxed to a plausible equilibrium with no
adhesion in it.* A coupling acting on 8 of an owner's 126 nodes genuinely does return an 8-row block,
so the adapter belongs here.

`_RoutedConnector` evaluates on the paired sites and places each row at its owner's index. `np.add.at`
rather than `[]=`, because two sites of one coupling may land on one node and an assignment would
keep the last instead of summing.

`PairedOwner` is the same restriction on the **other** read path: `CellWorld` reads a coupling's
geometry twice — through the resolver inside a step, and straight off `slot.owners` for the
out-of-step energy and residual probes. The second cannot go through a resolver, because
`StepContext` refuses to be constructed for anything that is not a real step, and `world.py` records
that refusal as correct.

## 4. Five mistakes this lane made, and the guard that caught each

None of these was found by inspection. Each was a refusal received while trying to do the simpler
thing.

| # | What was attempted | What refused it |
|---|---|---|
| 1 | pair a self-coupled owner by nearest point | `ConnectorGeometryError: a coupled site pair is coincident` |
| 2 | build `membrane_cortex_contact` on top of `build_vertical`'s own | `test_the_vertical_trio_alone_reproduces_the_verified_verticals_energy` — two springs on one pair, every term still looking correct while the total load doubled |
| 3 | hand a connector two arrays of different length | *"all four site arrays must have the same shape"* |
| 4 | return an 8-row block for a 126-node owner | `WorldCompositionError` (§3) |
| 5 | attach motor heads by hand, bypassing `propose_kinetics` | the moment-closure check: residual `4.20e-2`, never closed. Deleting `_attach_heads` closed it. |

Mistake 5 is the one worth keeping. A motor is **not** a spring between two arrays —
`CrossbridgeVariant`'s kernel docstring says so, and it is why `NmiiMotorConnector` takes
`(population, target_index)` rather than `(endpoints, spring, site_count)`. Its force is the head's
own strain through the population's kinetics, so writing `head_attached`/`head_target`/`head_site` by
hand produced a deterministic binding the kinetics never chose, and a net torque with it. **Not
attaching heads at all was the fix**: `propose_kinetics` attaches them inside the step, at its own
rate, into a coordinate it chooses.

## 5. `_positions` was wrong in the function that warned about it

`_positions` tries a list of attribute names. It started at three, and `nmii` publishes
`anchor_positions_um` and `backbone_nodes_um` — neither of which was on the list. `NmiiPopulation`
was therefore reported as *publishing no site positions* while publishing two arrays of them.

That is the same failure the function's own docstring warned about, committed inside the function
carrying the warning. It is the session's recurring error class: **a count derived from a
hand-written pattern cannot see what the pattern does not match.**

The list stays explicit rather than becoming "any `(N, 3)` float array", because the general version
would silently pick up a velocity field or a reference configuration and couple to it. The risk is
therefore stated in the docstring instead of removed: *a new owner publishing under a sixth name is
reported absent, and the report is the thing to read.*

## 6. Result

| | before | after this entry |
|---|---:|---:|
| couplings in the world | 2 | 14 |
| compartments no connector reaches | 10 | 3 |
| compartments singly ablatable | 11 / 13 | 13 / 13 |

The last row is the one that mattered for the PI's method. Both original couplings named the membrane
*and* the cortex, so ablating either left a world with **no coupling at all** — and a world with no
coupling cannot take a step, because `BalancePredicate` refuses a candidate whose ledger recorded no
adjoint pair. Every compartment can now be removed singly and the world still steps.

`ALEPH-PORT-3624` carries this from 14 couplings to 17 and from 3 unreached compartments to 1.

## 7. What is still absent, and why — read this before reading the count

`CouplingReport.absent` carries a written reason for every row not built. They fall into four kinds,
and only the first is a gap in this lane:

1. **needs construction arguments this module does not supply** — `BrownianRatchetConnector` and
   `NmiiMotorConnector` want a population, a clutch graph or a solver. `filopodium_membrane_tip`,
   `lamellipodium_membrane_contact`, `nmii_cortex_motor`, `nmii_sf_motor`.
2. **the endpoint is not in this world** — every `*_cytosol_transfer` row, because `cytosol` is
   refused by the coverage gate (`ALEPH-PORT-3621` §REFUSED).
3. **the row is an internal force, not a coupling** — `ecm_crosslink`, `dorsal_arc_crosslink` (§2.1).
4. **already in the world under another name** — `membrane_cortex_contact`, `membrane_erm_cortex`,
   built by `build_vertical` (§4, mistake 2).

`nmii_cortex_motor` and `nmii_sf_motor` are absent for a fifth reason worth naming separately: those
owners publish `NMII_HEAD` sites and **no module converts them to `ActinBindingSites`**.
`connectors_protrusion` supplies that conversion for the lamellipodium and the filopodium only.
Writing the same helper against `sf_arc`'s and `cortex`'s own site types is a change to a module this
lane does not own.

## 8. Evidence class

`UNVERIFIED` throughout. What is verified is *mechanical*, not biological: force closure on every
built coupling, moment closure on the motors, and that the vertical trio alone still reproduces
`build_vertical`'s energy to the digit. None of that makes a stiffness or a separation a measured
quantity.

## 9. Open, and proposed rather than decided

1. **A module citing a ledger is not checked against the ledger directory.** §0. The existing
   port-discipline guard reads citations inside ledgers; a `ALEPH-PORT-NNNN` in a module docstring is
   unchecked, and one was dangling through a full-repo run. Proposed as a new control in
   `tests/ports/`; not written here, because a guard belongs to the lane that owns the guard file.
2. **The `ActinBindingSites` converters for `sf_arc` and `cortex`** (§7) — two motor rows wait on
   them, and the module is not this lane's.

---

## 10. Source repository identity, source path and symbol

**None.** Port class `RE-DERIVED`: `/Users/sw1/ffn_cellsim` was not opened for this entry and no
source path or symbol was consulted. Every connector class instantiated here is Aleph's own, written
under its own earlier ledger; this entry supplies the `ConnectorSlot`s that put them in a world. A
source commit would be a citation to a repository that played no part.

## 11. Units, domains, singular cases, invariants

Units: length µm, force pN, stiffness pN/µm, energy pN·µm. Evidence class `UNVERIFIED`.

**The singular case is the whole subject of §2.1.** A central coupling is defined for separation
`d > 0`; at `d = 0` the direction `(a−b)/|a−b|` is `0/0`, and the connectors refuse it rather than
regularising. `COINCIDENT_UM = 1.0e-9` is the threshold at which this module declines to *offer* such
a pair, so the refusal arrives as a reported absence at build time instead of an exception at step
time. It is **not** a tolerance on the physics.

**Domain restriction that is easy to miss:** `site_pairs` requires distinct partners. On an 11–16
node owner, nearest-point collapses several sites onto one partner, and the survivors are reported
rather than deduplicated in silence — a coupling carrying three sites where its card assumes eight
has a stiffness that no longer means what it says (§2.2).

**Invariants held:**

| invariant | how it is held |
|---|---|
| `ΣF_a + ΣF_b = 0` on every built coupling | `AdjointPair`, with the two blocks independently normalised |
| net moment closes on the motors | `propose_kinetics` chooses the binding; §4 mistake 5 is what happens when it does not |
| the vertical trio alone reproduces `build_vertical`'s energy | `test_the_vertical_trio_alone_reproduces_the_verified_verticals_energy` |
| one owner, one coupling per pair | `already_built=` (§4 mistake 2) |

## 12. Positive control, negative control

**Positive control** — `test_the_vertical_trio_alone_reproduces_the_verified_verticals_energy`. With
only the membrane, cortex and their two couplings present, the world's total potential energy must
equal `build_vertical`'s to the digit. This is the control that caught mistake 2 in §4: a duplicate
`membrane_cortex_contact` doubled the load while every individual term still looked correct.

**Negative controls** — each of the five refusals in §4 is a negative control that was *received*
rather than designed, and three are pinned as tests:

| what must fail | control |
|---|---|
| a self-coupled owner paired against itself | `ConnectorGeometryError`, asserted in the coupling report's `absent` reason |
| a coupling built twice on one pair | the positive control above turns red |
| a motor with hand-attached heads | moment residual `4.20e-2`, never closing |

The moment residual is the sharpest of the three: **force** closure held throughout while the moment
did not, so a check that stopped at `ΣF = 0` would have passed a motor population applying a net
torque to the cell.

## 13. Numerical and precision envelope

Host `float64` throughout; no `float32` path and no ULP budget, since nothing here is a law kernel
graded against a reference. Force closure holds to round-off rather than bit-exactly — deliberately,
because `AdjointPair` forbids `force_b := −force_a` and that is what makes closure a check instead of
a tautology. The energy comparison against `build_vertical` is the one place a tight tolerance is
used, and it is tight precisely because the two paths should be computing the same sum.

`np.add.at` rather than `+=` on a fancy-indexed array is a **correctness** requirement, not a
precision one: two sites of one coupling may land on one node, and an assignment would keep the last
instead of summing them.

## 14. Production-backend residency and transfer

Host-side NumPy. Owners hold their positions and forces as host arrays and the couplings built here
evaluate on the host; the connectors that scatter through `ctx.backend.scatter_add` do so inside
their own `accumulate`, which this module does not touch. No device residency is claimed and no
transfer is performed by this entry.

`_ScatterBack` is a `np.ndarray` subclass whose `flush()` performs the scatter-add into the owner's
real array. Fancy indexing copies, so a connector writing into `forces[index]` would write into a
temporary and its contribution would vanish — silently, with every energy still finite. That is a
host-memory aliasing hazard rather than a device one, and it is the reason the class exists.

## 15. Comments and docstrings discarded

None to discard — nothing was copied. All prose is new. No provider identifier, path or comment
appears in `whole_cell_couplings.py`.

## 16. Acceptance, reviewer, rollback

**Status `PROPOSED`.** Reviewer: the PI. Acceptance evidence is `tests/scenarios/test_whole_cell.py`,
which at the close of this entry stood at 26 passed including four ablation controls.

Rollback: remove the `build_couplings(...)` and `build_motor_couplings(...)` calls from
`build_whole_cell`. The world returns to two couplings, the twelve rows return to `withdrawn` with
their reasons, and the membrane and cortex become unablatable again — which is the state §6 records
as the reason this entry was written.
