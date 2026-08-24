# PHASE 3 — fourteen connectors, zero families, 82 named questions

**Status:** DECISION QUEUE. Nothing here is a physics result. Produced 2026-08-20 by fourteen parallel
sessions, each owning one module under `aleph/world/families/`, coordinated by Lead session `d380041d`.

**The headline is that nothing was built, and that is the designed outcome.** Every module raises
`ConnectorGapError` and exposes a `FamilySpec` whose `blocked_by` names what it could not answer. A
`BondFamily` cannot be constructed without a `BondCount`, and `BondCount` refuses without *how many,
against what, for which cell, on whose authority*. **Fourteen sessions could not invent a single one.**

    modules 14 · families built 0 · named questions 82

---

## 1. What the 82 are about

| kind | count |
|---|---:|
| how many / density / against-what | **42** |
| kinetics — the rate law bond.py deferred | 8 |
| stiffness / reach / rest length / radius | 9 |
| scope — which cell, which state | 7 |
| authority — no ValidationGate or ModelContract row | 3 |
| other, including two "is this a bond at all" | 13 |

**Half the queue is one question asked fourteen ways: how many.**

---

## 2. Findings that are worth more than the queue

### 2.1 The ERM defect is not one defect. It is a class, and it has five costumes.

The original: one declared connector stood for a POPULATION, nobody had to answer how many, and the
answer turned out to be the icosphere vertex count. Fourteen sessions found the same shape elsewhere,
each in a different disguise:

| where | the number that wore a physiological label |
|---|---|
| ERM (original) | **mesh** vertex count. D1 confirms `bond.py`'s 231.8/µm² exactly — subdiv 7 = 163,842 verts over 4π(7.5)² = 706.9 µm². ⚠ And the three L0 runs of 2026-08-20 ran subdiv 3 = **0.9/µm²**, two and a half orders under even the diagnostic-only proxy |
| `sf_cortex_transient` (D13) | **placement** count — whatever fell inside a 0.06 µm capture radius. Native SF↔cortex median separation is **43×** that radius |
| `if_sf_plectin` (D9) | **twice**: the count is one-per-bead (`span/seg_um`) and the reach is `~l_seg`, so the candidate SET is discretisation too. Demonstrated live: same cell, same biology, `seg 0.50 → 660`, `seg 0.25 → 1260` |
| `actin_cap_linc` (D5) | a **builder signature default** — `n_cap: int = 2` at `engine/sf_population.py:440` |
| `if_nucleus_linc` (D4) | `compartments.py:116` `linc_per_node=1` — the LINC population size IS the nucleus vertex count, 642 at subdiv 3, **the same mesh ERM wore** |

⚠ **And D4 found the incumbent's LINC is not even this bond**: `compartments.py:500` pairs envelope
vertices against `pos_actin`, nucleus↔cortex ACTIN. No intermediate filament on either end.

### 2.2 Three declared connectors are not bonds at all

D2 and D6 reached this independently, from `membrane_cortex_contact` and `nucleus_cortex_contact`.
`bond.py`'s own opening rule: *"field interaction (steric, drag) needs no declaration because it depends
on position alone; a specific persistent pairing does."* Those edges declare `kinetics=False` +
`remap_on_accept=True` — the pair set is regenerated from proximity every accepted step, so **there is
no identity to persist**, which is the only thing a bond adds over a field.

D2's discriminator is measured, not asserted: `chemistry_card is None` **and** `remap_on_accept=True`
selects exactly those three across the whole architecture — `membrane_cortex_contact`,
`nucleus_cortex_contact`, `membrane_ecm_contact`. ⚠ D2 first claimed the pattern held over all 38, then
measured, refuted itself, and narrowed the claim in the module.

**One PI ruling disposes of all three.** And it needs a destination: `aleph/world/` has nowhere for a
non-penetration field — `bond.py` declares pairings, `laws_bind.py` binds per-primitive kernels, neither
is a neighbour query.

### 2.3 Five connectors claim to be "the first family with kinetics", and one of them should not be

`bond.py` defers attach/detach, the free list and snapshot twins to *"the first family that has
kinetics"*. D1, D4, D9, D12/D14 and D7 all reached it. ⚠ D7's warning is the useful one:

> a **processive motor** is the wrong first customer for that machinery. Dynein WALKS, so its bound
> state carries a walk velocity and a stall force as well as on/off rates, and a slip-bond shape
> motivated by a crosslinker family will not cover it.

**D1's ERM is the better first customer** — a passive Bell-slip tether whose kinetics the incumbent
already drives (`erm_bell_kmc_pair_kernel`), and whose shedding IS the mechanistic bleb onset.

### 2.4 The PI has already ruled on the ERM count, and the shortcut is pre-empted

D1 retrieved `VG-ERM-capacity` (PI-authored, `AC_DECISION_CARDS` 2026-07-22 P2 Card 2):

> *"Production `rho_ERM` latch stays FALSE (HOLD). 1/node and zebrafish 600 µm⁻² stay diagnostic-only."*
> *"NECESSARY-ONLY: the back-computed minimum `rho_ERM` is a LOWER BOUND, not the physiological density
> (do NOT adopt the minimum as the value)."*

⚠ **Three of `BondCount`'s four questions already have answers for ERM** — support is the measured mesh
area, cell is MCF7, and `k_erm = 4.6e3 pN/µm` is PI-ratified (KB-3.B1.6, Braunger 2014, `source_audit`
verdict OK). It is a **one-parameter block**, and the PI's own closing order is recorded on the card.

### 2.5 New coupling nobody had declared — D5

If one SUN trimer holds one KASH protein at a time, `actin_cap_linc` / `mt_nucleus_linc` /
`if_nucleus_linc` **compete for one finite SUN-site population** on the inner nuclear membrane. Three
families each sizing themselves against the same envelope over-subscribe it silently. Structurally
identical to PI card B1 (one NMII head bound by two motor connectors double-counts its force).
**Nothing in the arena expresses a shared site budget.** This is a new item, not covered by C1.

### 2.6 Two populations the Round-2 brief said existed, do not

* **D5**: `actin_cap_linc` needs a PERINUCLEAR CAP. `build/stress_fiber.py` builds **ventral** bundles
  only — FA-to-FA, basal plane, and it refuses an axis parallel to the substrate normal. A cap arches
  OVER the nucleus. D5 is blocked ahead of its count, not by it.
* **D8/D9**: the contract endpoint is `sf_arc` — transverse/dorsal **arcs** — and PHASE 1 built
  `stress_fiber`. Different structures, not a rename. Pinned as `xfail(strict=True)` in
  `tests/world/test_families_census.py` together with D4's `nucleus` vs `nuclear_envelope`.

### 2.7 A second KB row that cannot be checked against its own citation

D14: `MogilnerOster2003_BJ` — the paper that decides whether the Brownian ratchet is **elastic** (tips
bend and push, nothing attaches → not a bond) or **tethered** (a subset binds → that subset is the bond
population) — is registered in `source_audit` with **verdict OK and zero chunks**. Registered, audited,
unreadable. Same class as `KB-3.8`/`r_filo` in `PER_CELL_STRUCTURE_COUNTS_2026-08-20.md` §4.

---

## 3. The duplicate map, and the criterion four sessions converged on

    A family's identity is its CHEMISTRY CARD.
    A shared runtime class is not a duplicate. A shared force law is not a duplicate.

Recorded in `aleph/world/families/__init__.py` with each session's route to it. Verdicts:

| folded | into | card |
|---|---|---|
| `lamellipodium_cortex_seam` | `sf_cortex_transient` | `transient_actin_crosslink` |
| `filopodium_membrane_tip` | `lamellipodium_membrane_contact` | `actin_membrane_brownian_ratchet_contact` |
| `sf_cortex_transient` | `dorsal_arc_crosslink` | — |

**Not duplicates**, each argued from physics: the three LINC edges (three cards; and one of them is a
MOTOR, which cannot share a kinetic law with a spring), `if_sf_plectin` vs `mt_sf_spectraplakin`
(different molecules, different binding domains), `filopodium_cortex_root` (different card, and it is
D13's control).

⚠ **A duplicate verdict never answers "how many".** Recorded beside the criterion because the two are
easy to confuse.

---

## 4. What the PI is asked, in the order that unblocks the most

1. **Are the three cardless CONTACT edges bonds at all?** One ruling closes three modules, and it needs
   a destination in `world/` for a non-penetration field. (D2, D6)
2. **ERM density** — one parameter, three already answered, and the PI's own card sets the closing
   order. Unblocks D1 and makes ERM the natural first kinetic family. (D1)
3. **Which connector opens the kinetics machinery**, and not a processive motor. (D1, D4, D7, D9, D12/D14)
4. **The LINC isoform split** — one decision closes D3, D4 and D5 together. (D4)
5. **The SUN-site budget** — a new shared-capacity item. (D5)
6. **`sf_arc` vs `stress_fiber`, `nucleus` vs `nuclear_envelope`** — which name is right. (D4, D8, D9)
7. **A perinuclear cap population**, or D5 moves to the waiting half. (D5)
8. **Leading-edge width** — carried from `PER_CELL_STRUCTURE_COUNTS` §5. (D14)

## 5. What this does not do

It builds no family, selects no count, and ratifies no duplicate beyond recording each module's own
argument. Every number quoted is either read from this repository or from a source the quoting session
named; where a source could not be read, that is said rather than worked around.

---

## 6. ⚠ A defect in one of the two connectors that actually run

Found by D10, **verified by Lead against the source rather than relayed.**

`dorsal_arc_crosslink` is one of only two connectors with a runtime object — it is `DEVICE_RUN` in the
committed census and it is named in `STATE.md` (b)'s row *"The composed native cell bound 2 of 38
connectors"*. Three things are wrong with it at once.

**(a) The contract and the runtime contradict each other.** `ConnectorContract`'s fifth and sixth
positional fields are `kinetics` and `commit_on_accept` (`contracts.py:351-352`), and the declaration
passes `True, True`:

```
contracts.py:905   ConnectorContract("dorsal_arc_crosslink", ConnectorFamily.TRANSIENT_ACTIN,
                                     "sf_arc", "sf_arc", True, True, ...)
```

while its own runtime says the opposite, in its docstring:

```
sf_mechanics.py    "No-op: ``dorsal_arc_crosslink`` declares no kinetics, so it advances no
                    irreversible state."
```

**One of the two is wrong.** If the contract is right, a transient crosslink binds and never unbinds —
a **permanent weld**, which `CLAUDE.md` forbids by name and which `sf_cortex_transient` refuses to
become. If the runtime is right, the contract over-declares and the connector is not kinetic at all.

**(b) It has no chemistry card.** `chemistry_card` defaults to `None` and this declaration does not set
one — the only member of `ConnectorFamily.TRANSIENT_ACTIN` without one. D2 reached the same observation
from the CONTACT side and called it *"a kinetic bond with no identity"*. Under the ratified
`DUPLICATE_CRITERION` its family membership **cannot be evaluated**.

**(c) It reached `DEVICE_RUN` without ever being asked how many.** `sf_population.py:576-581` pairs
every dorsal free end to its nearest arc apex — no radius test, no density, no count:

```python
for fn in dorsal_free_nodes:
    d = np.linalg.norm(arc_xyz - pos_flat[fn], axis=1)
    joints.append((fn, arc_apex_nodes[int(np.argmin(d))]))
```

The population size is `len(dorsal_free_nodes)`, which is a **builder placement**. It is therefore the
ERM class again — and it is in the one connector this project has cited as bound.

⚠ **What this does NOT say.** The committed row is about kernel LAUNCH — `STATE.md` records it as
*"the verdict is the kernel-LAUNCH observation ALONE"* with magnitude explicitly blocked. Nothing here
touches that. What it does say is that the edge's COUNT was never a physiological quantity, so the row
should not later be read as evidence that a dorsal-arc crosslink population was ever sized.

**For the PI:** which of the two declarations is wrong, and does this edge carry a chemistry card.

---

## 7. `BondFamily` cannot express the six `*_cytosol_transfer` edges — a type problem, not a count

Reported by session E (cytosol / FLUID_VOLUME) on 2026-08-20 and **verified here against source, not
relayed.** This is not one of the 82 questions above, because those all ask *how many*. This one says the
declared data structure cannot hold the edge at all.

**What `BondFamily` is.** `bond.py:189-190` declares the endpoints as

```python
node_i: npt.NDArray[np.int64]
node_j: npt.NDArray[np.int64]
```

and `component_pairs` resolves them at `bond.py:208-209`:

```python
a = arena.population_of(Kind.NODE, i) or "<unclaimed>"
b = arena.population_of(Kind.NODE, j) or "<unclaimed>"
```

`Kind.NODE` is **hard-coded**. A bond is a pair of global NODE indices, and the query that derives which
components it joins can only ask about nodes.

**What the cytosol end actually is.** `immersed_transfer.py:141-152` spreads `-α·V_node·grad(p)` onto the
solid force array through a **Peskin 4-point interpolation**. The cytosol-side endpoint is a *stencil
around a position*, not an ID. There is no cytosol node to bond to, and manufacturing one would convert
the interpolation into a lumped surrogate — which is the architectural principle inverted.

**Why this is on the PI's queue rather than a module owner's.** The fix is one `kind` argument once
`Kind.GRID_CELL` exists. But it cannot be written before that decision, and it changes a signature in
`world/bond.py` — which six family modules will then bind against. Sequencing it wrong makes six modules
edit the same file.

**Scope.** Six of the nine cytosol-blocked D connectors are `*_cytosol_transfer`. The remaining three are
unaffected by this item and unblock on the population alone.

**For the PI:** this is a consequence of decision (1) `Kind.GRID_CELL`, not a separate decision — listed
here so that approving `GRID_CELL` is not mistaken for unblocking six connectors that still need a
`bond.py` signature change first.

---

## 8. A cortex population cannot hold an anti-parallel pair — `polarity` is one scalar for all of it

Raised by session H (active channels) on 2026-08-20 while hardening `contraction.py` against the cortex
census, and **verified here against source.** H stated it as *"my check cannot see anti-parallel"*. The
source says something stronger.

`build/cortex.py`'s `StrandPopulation` declares

```python
seed: int
polarity: int          # ← one value, for the whole population
topology_bytes: int
```

with the docstring *"polarity: ``+1`` if the barbed end is the last node of each filament."* **Every
filament in the population.** There is no per-strand polarity array.

**What follows.** A bipolar NMII minifilament binds two ANTI-PARALLEL filaments; a station whose two
filaments run the same direction is a spring wearing a motor's name (session H's own wording, and
`minifilament_topology.py` agrees). In a population where polarity is a single `+1`, **no two filaments
are anti-parallel as recorded** — so an anti-parallel station is not merely unverifiable, it is
unrepresentable.

**Why it surfaces now.** On SF this was hidden behind `STATE.md` (e) 5 (`N_stations: 0`, placement
forbidden). G placed its 442 minifilaments on the **cortical shell**, where nothing blocks first, so the
cortex contraction channel reaches this before (e) 5 is touched.

**This is the ERM defect class in a new costume.** Not a hidden count this time — a hidden *degree of
freedom*. Polarity was recorded per population because that is what a builder needed to place nodes, and
the label reads as a physiological property of each filament. A real cortex is a mixed-polarity network;
that is the premise under which cortical NMII generates tension at all.

**For the PI:** does `StrandPopulation` carry per-strand polarity. It is not a free change — polarity
selects `tip_node` (`strand.py:202`), which is what every barbed-end law addresses, so a per-strand
array propagates into the growth kernels. Owners: A2 (population record) and G (which stations are
legitimate). **Until it is decided, the cortex contraction channel must not be opened** — binding it
would produce tension from same-polarity pairs and report it as cortical.

---

## 9. CORRECTION — `n_heads_per_side` and `F_stall_head` are not open. They are already closed by default.

Raised by session G on 2026-08-20 against **this document's own item #3**, and verified here in both the
code and the KB database. Item #3 asked the PI to choose *10 (AFINES) vs 28 (Billington)*. That is the
wrong question, and asking it would have let a decision that is already live pass as pending.

### 9.1 The runtime already uses 30 heads per side and 2.0 pN per head

`laws/myosin_linear.py:34-35`:

```python
N_HEADS_MINIFIL = 60           # KB-3.18 ~30 myosin-II dimers × 2 heads
F_HEAD_PN = 2.0                # per-head stall force [pN] (KB-3.18 ~2 pN)
```

`resolve_myosin()` (`myosin_linear.py:63-64`) takes **no arguments** and returns `MyosinMinifilament()` —
every field defaulted. So `n_side = 60 // 2 = 30` and `f_stall_pn = 30 × 2.0 = 60 pN`.

`world/active/contraction.py:176` — `motor = motor or resolve_myosin()`. The active-channel binding
inherits both numbers by default.

G declined to accept `F_stall_head` as a builder argument, on the grounds that *a geometry builder must
not hold an unratified magnitude*. That was right, and it does not help: **the value walks in one layer
over.**

### 9.2 Neither candidate traces to a primary measurement — queried, not recalled

| row | status / conf | what the DB actually holds |
|---|---|---|
| `KB-3.18` | **verified / High** | title *"Cortex molecular composition (5-component)"*; citations are **Salbreux2012_TCB, ChughPaluch2018_JCS, Murrell2015_NRMCB, Charras2008_BJ, Ferrer2008_PNAS, Gittes1993_JCB, Isambert1995_JBC, Goldmann2002, PollardBorisy2003_Cell** — **cortex-composition reviews. Neither Billington nor Stam is among them.** |
| `KB-DRAFT-3.B-14` | draft / Medium | citations field ends, in its own text: *"Hill 1938; Kovács 2003 JBC; Stam-Hocky 2015 PNAS; Billington 2013 JBC; Veigel 2002 NCB; Stachowiak 2009 — **ALL absent from SE/KU/corpus**"* |
| `source_audit` | — | `Billington2013` → **`CHECK`** · `StamAlbertsGardelMunro2015` → **`CHECK`** |

So: the 10-head claim rests on a draft row that **declares its own citations absent**, and the 30-head
claim rests on a verified row whose citations are **reviews of cortex composition, not measurements of
minifilament head count**, while both primary papers fail the citation audit. Per CLAUDE.md a KB source
may not be cited in a deliverable unless `source_audit.verdict` is OK; neither primary is.

⚠ **This is the KB-3.21 / KB-3.8 pattern a third time in one day.** `verified` grades the row, not the
scope of the number inside it — and here the row is verified *as a cortex composition claim*, which is
not what the runtime reads out of it.

### 9.3 The question, restated

**Not** *"10 or 28?"* — but:

> The runtime law already runs at 30 heads/side and 2.0 pN/head. Neither this value nor its
> alternative traces to a primary measurement that passes the citation audit. Ratify the standing
> default as a declared axis with a band, or hold the runtime until `Billington2013`'s `CHECK` is
> resolved?

### 9.4 A separate contradiction found in the same query

`KB-DRAFT-7-05` (draft / **High**) is titled *"…ensemble force via Stam minifilament with a **Hill**
force-velocity law"* and, per G, calls the Hill hyperbola the sanctioned oracle. This **contradicts the
2026-07-07 PI ratification** that non-muscle NMII force-velocity is LINEAR, which `myosin_linear.py:7-10`
records and which session H's channel inherits. `myosin_linear.py` notes *"CLAUDE.md to be annotated"*;
**the KB row was never updated.** The engine and the knowledge base now say opposite things about the
same law, and the KB row is the one carrying `confidence: High`.

**For the PI:** the ratification stands in code. Does the KB row get retired, superseded, or re-scoped —
and who is allowed to do it, given that `ValidationGate`/`ModelContract` rows are PI-authored only.

---

## 10. `contraction.py` cannot consume a head-resolved population — the law it binds is two-anchor lumped

Also from G, verified here. `laws/myosin_linear.minifilament_kernel` reads

```python
i = links[t, 0]; j = links[t, 1]
Fmag = f_stall * (1.0 - v_slide[t] / v0)
wp.atomic_add(force, i, f); wp.atomic_add(force, j, -f)
```

The whole minifilament is **one link between two actin anchors**. No head appears. G's population is 34
nodes per minifilament (14 backbone + 20 heads) and head resolution is its entire reason to exist.

Feeding the 442 into this channel demotes them to what `engine/nmii_actuator.py` isolates as
`LegacyNMIIKind.TWO_ANCHOR_LINEAR`, whose own guard at `nmii_actuator.py:492` raises
*"aggregate or two-anchor NMII mechanics is diagnostic-only"*. That is CLAUDE.md's
*"Mechanistic over lumped, at every design decision"* inverted, and it is **a design decision, not a cost
tradeoff** — so no session may take it.

⚠ **Consequence for item #8.** Item #8 held the cortex contraction channel shut pending per-strand
polarity. G's finding is a second, independent reason it must stay shut, and it does not dissolve if
polarity is fixed. G also retracts its earlier *"cortex is open"* — that meant a bond TARGET exists, not
that a station table can be produced: its 442 sit on random tangent axes on the shell and straddle
nothing. **`STATE.md` (e) 5's curved-band argument is SF-specific, but "therefore cortex is open" does
not follow — cortex is the same wall for a different reason.**

**If a cortex contraction slice is wanted today**, G names the honest route and it is not G's population:
place stations from the **cortex builder's own filament pairs** and label the result a **lumped two-anchor
diagnostic**. The head-resolved population then survives undemoted.

---

## 11. The KB failure modes found on 2026-08-20 — three types, three different remedies

Five KB rows failed today, found by four sessions working independently. They are not one problem and
must not be filed as one: each type is created by a different mechanism and closed by a different fix.

### Type 1 — `verified` grades the ROW, not the SCOPE of the number inside it

| row | grade | what broke |
|---|---|---|
| `KB-3.8` | verified / High | `r_filo ~5 µm` — of its 3 citations only 1 is in the corpus, and that one contains no `r_filo` |
| `KB-3.21` | verified / High | `V_cyto ~4200 µm³` is a **R = 10.0 µm** sphere; the arena's sourced MCF7 `R_cell = 7.50` gives 1,211 µm³ — **3.47×** |
| `KB-3.18` | verified / High | read for minifilament head count, but the row is titled *"Cortex molecular composition"* and its 9 citations are composition **reviews** — neither Billington nor Stam is among them |

**Mechanism:** the schema has no field for the scope a value is true for, so a row verified for one
purpose reads as verified for any purpose that quotes it. `bond.SourceClass.UNRATIFIED_PROXY` already
names this exactly — *"perfectly SOURCED and still wrong for the run … the SCOPE is what was never
recorded."*

**Remedy:** a scope field, and a citation-audit rule that a row may only be quoted for the claim its own
title makes.

### Type 2 — the row carries its own warning and the citing session reads only the value

`KB-3.31` (nucleoplasm). Sessions E and this one both cited its bulk η ≈ 52 Pa·s as an argument for a
separate nucleoplasm grid resolution. The row's own text says *"Wire **pore** viscosity into the fluid
channel … poroelastic, **not one lumped viscosity**."* Darcy mobility is `k/µ_pore`; bulk rheology does
not enter CFL at all. **The row we cited refuted the inference we drew from it.** Caught by session F.

**Mechanism:** a value is machine-readable and a caveat is prose, so the value travels and the sentence
does not.

**Remedy:** citation discipline, not schema — read the row, not the field.

### Type 3 — the correction is written where it cannot take effect ⚠ the quietest one

`params_i0b2.yaml:165-179` already declares `eta_nucleoplasm: value: null`, band 1–100 Pa·s, with
`surface_to_pi: "Keep the nucleoplasm viscosity SEPARATE."` **It is already fixed — in a file CLAUDE.md
states plainly that "editing one changes nothing at runtime."** The runtime keeps
`eta = ETA_CYTOPLASM = 65.9` (`laws/motility_warp.py:407`, applied to nucleus beads at `:429`).

**A ledger says one thing, the runtime does another, and only one of them is authoritative — while the
ledger entry reads exactly like a fix that landed.**

**Mechanism:** two stores, one authority, and no check that they agree. This repo has several YAML
ledgers, so the reach of this type is the widest of the three.

**Remedy:** a test that fails when a ledger row's declared value disagrees with the constant its
`symbol` names in code. Not proposed here — it is a gate contract and PI-authored.

### 11.1 What is NOT wrong here, so it is not quoted as such

⚠ The `eta_nucleoplasm` ledger text says the code was *"MIS-SET to 65.9"*, which reads as an
order-of-magnitude error. It is not. Against `KB-3.31`'s bulk nucleoplasm 52 Pa·s (Tseng 2004) it is
**+27%, inside the row's own 1–100 Pa·s band.** A **provenance** defect, not a magnitude defect. Session
F insisted on this and was right to; quoting it as a magnitude error would have been today's fourth
KB-scope accident, committed while cataloguing the first three.

⚠ **The real blocker is neither 65.9 nor 52** — it is that the band spans 100× and no MCF7 point value
exists. The ledger says so itself.

### 11.2 ⚠ RETRACTED as a finding — the July audit already held all of it

This section first reported the `physical_node_gammas` live-caller status as something four sessions
discovered today. **It was already written down on 2026-07-28.**
`T10_GAMMA_NODE_AUDIT_2026-07-28.md:44-47`:

> *Every caller is a legacy `ff_*` driver under `scripts/` (the `scripts-legacy` lane):*
> *`physical_node_gammas` — `ff_sc0_profile.py`, `ff_membrane_cell.py`, `ff_stiffness_sensing.py` (×2),*
> ***`ff_crawl_on_substrate.py`. Five call sites, four files.***

The same audit's table names `motility_warp.py:407-431` and its `6πηR` nucleus Stokes beads at
`ETA_CYTOPLASM = 65.9 Pa·s` explicitly. The tension this section recorded — *"the doc calls it legacy,
the directory says production"* — is resolved by that audit's own lane label, `scripts-legacy`.

**What survives as a finding is one judgement, and it is F's:** 65.9 against `KB-3.31`'s bulk
nucleoplasm 52 Pa·s (Tseng 2004) is **+27%, inside the row's own 1–100 Pa·s band. A provenance defect,
not a magnitude defect.** The `params_i0b2.yaml` phrasing *"MIS-SET to 65.9"* reads as an
order-of-magnitude error and must not be quoted as one. Everything else here is a re-derivation.

⚠ **Two of this section's own claims were also wrong and are corrected here.** (a) The defect is not the
fallback value: `S.get("eta", ETA_CYTOPLASM)` also admits `mda_mb_231 → 12.0`, which is *also* a
cytoplasm viscosity — the defect is the design statement *"the single physiological viscosity η"*, and
filling the scenario dict does not fix it. (b) *"contamination zero"* was too strong: `STATE.md` (c) 2
and (c) 8 retire the **published** outputs, so re-running that driver yields new numbers not yet on any
blocklist. The correct statement is **"nothing currently quotable rests on this value"** — the upstream
is blocked, the path is not clean.

⚠ **And the guard is weaker than either session said.** `ac_gate_t10_medium_native.py:514` reads
*"medium_exterior.assert_single_dissipation_owner keeps it there"* — that is **prose inside a run
record's `notes` string literal, not a call.** The only real call is `test_medium_exterior.py:352`, over
a **hard-coded 2-tuple** (`medium_exterior.py`, `medium_stokes_analytic.py`). The guard is not broken —
its docstring scopes it to *"no medium-bearing module"* — but **the ratchet is maintained by hand**: a
new medium-bearing module is uncovered until somebody adds it to that tuple. T10 lane's call.

### 11.3 Type 4, and it is not a KB defect — it is ours

**An answer that already existed in `aleph/docs/v2_audit/` was re-derived from source by three sessions
in sequence.** F found the site, E corrected F twice, F corrected E once, this session corrected its own
truncated sweep — and the July audit had named the file, the line, the value and the lane a month ago.
Everyone grepped the code first and opened the audit last.

**Mechanism:** with 18 sessions running in parallel, the audit corpus is large enough that reaching for
`grep` is faster than knowing which memo already answers the question. The failure produces *correct*
results, which is why nothing catches it — it costs time and it manufactures false novelty.

**Remedy:** a reading-order rule, not a schema change and not a KB fix. Worth measuring recurrence
across today's sessions before writing one.

⚠ **Method note, recorded because it is the same defect class this document catalogues:** this session
reported the caller sweep as exhaustive from a `grep` piped through `head -12`; there are 22 matches.
Session F made the identical mistake in the same hour (`head`, 10 lines shown, 16 real). **A check
truncated before it reaches its counterexample reads as complete.** Two independent sessions hit it on
the same command in the same hour, which makes it a tooling habit rather than an individual slip.

---
## 12. Half the audit corpus cites dead directories — but only ONE thing needs a person

Found 2026-08-20 by session F, extended by E, counted three ways. `1b895ec0` (2026-08-09,
*"rename ffn_sim to aleph and split the package into layers"*) moved every module; **neither
`aleph/ff/` nor `aleph/ac/` exists.** The memos still cite code by those paths.

⚠ **This section was first filed claiming 13 of the paths need a human to choose between the canonical
and the frozen file. That was wrong, and §12.2 replaces it.** The correction is recorded rather than
edited away because the mistake is instructive: basename collision was read as authorial ambiguity.

### 12.1 The count, and which regex to use

| rule | memos |
|---|---:|
| ① single path segment only — `ac/x.py` | 25 |
| ② `.py` required, subdirectories allowed | 37 |
| **④ ② plus fully-qualified `aleph/ff\|ac/`** | **39 ← use this** |
| ③ ② plus glob / directory forms | 43 — **over-counts** |

⚠ **Numbers move.** Sessions are writing audit memos right now; rule ② went 36 → 37 within half an hour.
Treat any figure here as a snapshot. This session reproduced rule ④ at **37** shortly after F measured
39 — the corpus changed between the two runs, which is the point rather than a discrepancy.

⚠ **Do not use rule ③.** F opened the four memos it uniquely catches, and one is
`CORTEX_STRUCTURE_AUDIT_2026-07-23.md` — **`STATE.md` (d) 7, in the boot read-set.** Its two hits are
`ac/cell_assembled/` (a live OUTPUT path, `aleph/outputs/ac/cell_assembled/`) and `ac/weave` (a lane
label; `aleph/components/weave/` exists). **Neither is a stale code path. (d) 7 is NOT affected by this
item** — attaching a false positive to the memo that decides which cortex results are quotable is
expensive.

The regex for whoever fixes this:

```
grep -rlE '(^|[^A-Za-z0-9_/.])(aleph/)?(ff|ac)/[A-Za-z0-9_/]+\.py' aleph/docs/v2_audit/*.md
```

⚠ E's first sweep undercounted at 28 because its pattern was `[a-z_]+\.py` — no `/` inside the path. So
it dropped every path with a subdirectory, **which is exactly the set that looked most dangerous**
(`aleph/components/nucleus/envelope.py`, `aleph/engine/stress_fiber.py`, `aleph/components/solid/microtubule.py`). A narrow regex
misses precisely the entries that matter. This session's own first count was 23, for a different reason:
it matched `ff|ffn_sim|hoomd` and omitted `ac/` entirely.

### 12.2 ⚠ RETRACTED — the 13 "ambiguous" paths are not ambiguous. Chronology settles them.

The retracted claim: 13 cited paths resolve to two or more current files, one canonical and one frozen,
so *"a person must decide which the memo meant."* Verified basename collisions, e.g.

```
envelope.py    →  world/build/envelope.py     |  components/nucleus/envelope.py
protrusion.py  →  world/active/protrusion.py  |  engine/protrusion.py
```

**But the memos were written before the split.** Measured over the corpus:

| | memos | citations |
|---|---:|---:|
| written BEFORE the 2026-08-09 rename | **39** | **330** |
| written after | 3 | 3 |

When `ARCHITECTURE_COMPLETENESS_2026-07-25.md` wrote `aleph/components/nucleus/envelope.py`, **there was exactly one
file with that name.** `world/` did not exist — it was created 2026-08-15 (`f77a4daf`), six days after
the rename and three weeks after the memo. **The author could not have meant the canonical file, because
it did not exist.**

And git holds the mapping already:

```
R099  ffn_sim/ff/motility_warp.py     →  aleph/laws/motility_warp.py
R100  ffn_sim/ac/nucleus/envelope.py  →  aleph/components/nucleus/envelope.py
R100  ffn_sim/common/compartments.py  →  aleph/laws/compartments.py
```

`git log --follow` answers every one. **The collision was in our matching method — basename — not in the
memos.** Path re-keying is mechanical, 330 of 333 citations, no decision required.

### 12.3 What actually needs a decision — one policy, not thirteen items

Re-keying the path does not re-key the **claim**. `T10_GAMMA_NODE_AUDIT_2026-07-28.md:51` argues safety
from *"`aleph/ac/**` imports none of them"*. Mechanically re-keyed that reads *"`components/` and
`engine/` import none of them"* — still true. But the 2026-08-20 reader wants to know whether it is true
of `world/`, and **`world/` postdates the memo by 18 days. A memo cannot answer a question about a tree
that did not exist when it was written.**

> **THE DECISION: may a conclusion written before `world/` existed be read as a conclusion about
> `world/`?**

* **No** → all 330 re-label as facts about the frozen port source, and anything needed about the
  canonical engine must be re-established.
* **Yes** → a transfer rule is required, stating under what conditions it carries.

**Recommendation: NO by default, with symbol-level re-confirmation as the promotion path.** That is how
T10 actually survived today — its path argument is dead, but its own criterion is symbol-level (*"a
disjoint set of force kernels … none of which carries a drag coefficient"*), and three sessions
independently re-confirmed on the current tree that `physical_node_gammas` is imported by nothing in
`world/`, `laws/` or `engine/`. **The re-confirmation is mechanical — it is an import-graph query — so
the rule is cheap to apply and does not need thirteen judgement calls.**

### 12.4 Boot-path memos: 2

`AC_EXECUTION_PLAN_2026-07-25.md` (17 citations) is `STATE.md` (d)'s minimum read-set;
`AUDIT_WHOLE_REPO_2026-07-25.md` (8) is in the source-of-record list. Every session is told at boot to
read two memos that point at directories which do not exist — **paid once per session, 18 times today.**

### 12.5 The 6 unresolved paths are two different failures, not one

`ac/nucleus/chromatin_net.py` has no deletion record. It is an item on
`COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md:165`'s **"Write:"** list — **never written.** Session F's
`world/build/chromatin.py` (`13dfcff7`) **is that item**, in the new layout.

⚠ **Partial, and F said so first.** The same paragraph also requires nucleus tangents (Helfrich
Gauss-Newton, knee-aware areal, rank-1 volume, WLC pair), a rigid indenter and a nuclear socket
factory — none built — and states *"without them no indentation equilibrium is converged and nothing
quantitative in N2–N5 is admissible."* **N2–N5 are not open.** The other five
(`arclength_multigrid.py`, `fa_clutch_connector.py`, `membrane_shell.py`, `probe_load.py`,
`ff/ecm_modulus_probe.py`) are unsorted and need the same deleted-vs-never-written split.

### 12.6 ⚠ DONE — 290 of 299 re-keyed mechanically, and the residue is a different thing entirely

Executed 2026-08-21. The resolver walks git's full rename chain forward from each cited path until it
stops moving, then checks the destination exists:

```
cited unique            113        citations   299
resolved by rename chain 106        AMBIGUOUS     0
unresolved                 7
```

⚠ **Zero ambiguous.** §12.2's retraction is confirmed by execution, not only by argument: the "13
paths a person must adjudicate" were basename collisions in our matching, and git resolved every one
without a judgement call. **290 citations rewritten across 37 memos.**

**Paths only.** Decision 12 is open and nothing here promotes a conclusion: a memo that said
*"`ac/cell/driver.py` imports none of them"* now says *"`components/incumbent/driver.py` imports none
of them"* — the same claim about the same file under the name it has.

Where the 106 landed, which is the canonical/frozen split made visible:

| destination | paths |
|---|---:|
| `aleph/laws/` — the kernel library `world/` binds | 35 |
| `aleph/engine/` — frozen port source | 30 |
| `aleph/components/incumbent/` — frozen | 17 |
| `components/{weave,fluid,solid,motor,nucleus}` | 20 |
| other | 4 |

### 12.7 The residue is SIX FILES THAT WERE NEVER WRITTEN

Not stale paths. `git log --all --diff-filter=A` finds **no commit that ever added any of them**:

| path | where it was planned |
|---|---|
| `ac/nucleus/chromatin_net.py` | `COMPARTMENT_VALIDATION_TRACKS:165` **"Write:"** — session F's `world/build/chromatin.py` IS this item, partially |
| `ac/engine/probe_load.py` | `:123` **"Write:"** — point / dipole / patch / **indenter**, with adjoint reaction accumulation |
| `ac/engine/membrane_shell.py` | `:140` T4.0, *"land the membrane as a genuine state owner"* |
| `ac/engine/fa_clutch_connector.py` | `:255` M10, *"**one** `fa_clutch_connector.py`, not two"* |
| `ac/cell/arclength_multigrid.py` | `FINE_MESH_MULTIGRID_DESIGN:238` |
| `ff/ecm_modulus_probe.py` | `COMPARTMENT_VALIDATION_TRACKS:446` |

⚠ **Session F guessed this for `chromatin_net.py` from one deletion record. It holds for all six.**

**And one of them is load-bearing for a project goal.** `probe_load.py` was to carry the INDENTER, and
the throughput target this project is measured against is *AFM indentation, ~3 min per cell*. **There
is no indentation forward operator anywhere in this repository** — session C established that
independently on 2026-08-20 while auditing the detectors, finding that `tether` is a bead-pulling
operator and the AFM apparatus is spherical indentation. Two sessions reached the same absence from
opposite directions, a day apart. It is not a stale reference; **it is a hole where a deliverable was
planned.**

### 12.8 What remains

| | work | who |
|---|---|---|
| (a) | the transfer policy in §12.3 | **PI — still open (decision 12)** |
| (b) | ~~re-key~~ | **done, 290 citations** |
| (c) | ~~sort the unresolved~~ | **done: all six were never written** |
| (d) | decide whether the six unwritten items are still wanted, `probe_load.py` first | **PI** |
