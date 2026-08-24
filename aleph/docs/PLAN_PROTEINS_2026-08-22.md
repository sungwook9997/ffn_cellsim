# Plan — myosin motors, integrin, and the rest of the missing molecules

**2026-08-22, Lead session, `engine/main`.** Written on the PI's instruction after the question
*"마이오신 모터랑 인테그린 등 부족한 것들 어떻게 할지"*.

**What this document is.** An ordering and a set of preconditions. It proposes no new representation,
because the representation question is already answered, and it proposes no new module, because seven
of the ones needed already exist and are unbound.

⚠ **What it is NOT.** It is not a result, it establishes no claim, and it authorises nothing. Three of
its steps are PI decisions and are marked as such. Every measurement in it was taken on the built cell
or read from a run record; where a number comes from a document rather than a measurement, the document
is named.

⚠ **It also corrects one step of the audit it builds on** — see §7.

Source it extends: [`ADVERSARIAL_AUDIT_HOW_TO_REPRESENT_A_PROTEIN_2026-08-21.md`](ADVERSARIAL_AUDIT_HOW_TO_REPRESENT_A_PROTEIN_2026-08-21.md).

---

## 1. There is no representation decision left to take

`aleph/world/arena.py:76-78`, in the paragraph that closed the primitive set, already rejected all
three candidate primitives *because each decomposes into existing ones* — and those three are exactly
the three options the question keeps being asked in:

| the question's option | the arena's answer | live today |
|---|---|---|
| add nodes | a molecular attachment is a **`BOND`**, far address on a `SEGMENT` | **0 / 1,000,000** |
| smaller nodes | a binding **site is not an allocation** — resolved from live geometry | never exercised |
| a probability cloud | a species pool is a **`GRID_CELL` channel** + one conserved scalar per consuming range | **0 / 4,000,000** |

Measured in today's native run (`world_phase4/tau3_seed1.json`, `census.live`): `bond` **0**,
`grid_cell` **0**, against `node` 4,558,554 and `segment` 4,280,891.

**So the design is right and both primitives it reserves for exactly these molecules are unused.**

**A node is 50 nm** — the measured median bond length on cortex, microtubule and stress fibre — which
is about **nine actin monomers**. Against that ruler integrin (10–20 nm) and α-actinin (35 nm) are
**sub-node**, and giving either its own node would resolve a receptor more finely than the actin it
binds, so the finer object's detail could not be felt by the coarser one it acts on. **These molecules
are links, not extended things, and a link is a `BOND`.** The one that genuinely wants nodes is
**laminin**, because a basement membrane is a network with its own mechanics — which is what a
`STRAND` population is.

## 2. Cost is not the obstacle, so no argument may reach for it

```
node  headroom   12,000,000 − 4,558,554 = 7,441,446 free (62%)
bond  headroom    1,000,000 −         0 = 1,000,000 free (100%)

whole-cell basement membrane at 50 nm      565,520 nodes   (7.6% of free)
every integrin in the cell, as BONDs      < 200,000 bonds  (< 20% of capacity)
```

## 3. What is actually blocking, measured

### ① No acceptance criterion — `STATE.md` (e) 1, PI queue 18

Today's 30,000-step native run: **0 accepted, 0 rejected, predicate `undefined`.** The charter says a
kinetic connector *"commits ONLY on an accepted physical step"*. **Every chemistry in this design is
therefore unreachable by construction.** An integrin clutch that binds and unbinds is a kinetic
connector. Adding one today produces a population that cannot transition, and whose correctness
**cannot be falsified**, because nothing it would do is permitted to commit.

Sharpened 2026-08-21/22: the equilibration search's candidate grid is `[0 … 14,625]` and the chosen cut
is an **endpoint on all three seeds** — the search is censored by its own range bound. Widening it lets
one seed of three produce a τ; a second never does. **One seed of three is a D-3 input, not D-2.**

### ② The cell is attached to nothing — PI queue 22

`world/` binds 11 of 47 law modules. Verified by grep today, all eight adhesion/ECM modules bound in
**zero** places: `fa_anchor`, `fa_clutch_warp`, `fa_ecm`, `fa_maturation`, `ecm_library`,
`ecm_mechanics`, `ecm_mikado`, `substrate`. The "FA-to-FA" stress fibres have nothing at either end.

### ③ 🆕 The motors do not touch the actin — `STATE.md` (c) 20, measured 2026-08-21

**This blocker did not exist when the audit was written; it was measured the same evening.**

```
442 minifilaments · 32,708 NMII nodes · 26,520 crossbridge capacity if all stationed
nearest cortex node, minimum over the whole population   0.2501 µm
stations at reach 0.20 µm                                0 / 442
```

The **closest** minifilament in the population is farther from actin than the capture proxy. Not "most
cannot bind" — none can, and no rate constant reaches across a gap.

⚠ **③ is independent of ①.** Ruling the acceptance criterion does not give myosin anything to bind to.

⚠ **And upstream of ③ is a guard that cannot fail for what it appears to protect.** `build/nmii.py`
raises when `excursion > thickness/2` — a FIT test — while `nmii.py:156` draws the per-filament radius
over the FULL shell thickness with no inset, so a filament at the outer edge sits that excursion past
the shell whatever the guard returns. Until that is fixed, **any radius chosen from the measured window
is chosen against a shell whose outer edge is not where the build puts its nodes.**

---

## 4. Per-molecule plan

### 4.1 Myosin II motor — the population exists; the coupling does not

**Already standing:** 442 minifilaments, 32,708 nodes, bipolar geometry with head arms built
unstrained; `families/nmii_cortex_crossbridge.py` implements the `BOND`-with-far-address-on-a-`SEGMENT`
pattern; TIER 0 `station_census` needs none of the eight PI-GAPs and runs.

| step | what it is | who |
|---|---|---|
| M1 | Fix the excursion guard: inset the draw, or make the guard test containment instead of fit | **PI** (it changes the cell) |
| M2 | Rule the shell radius. Measured window ≈ 7.20–7.35 µm; resolved only to the ladder's 0.05 step | **PI**, queue 14 |
| M3 | Re-run TIER 0. Expect `nearest_cortex_node_um.min` well under the capture proxy and stations at every reach | Lead, minutes |
| M4 | Rule the eight crossbridge PI-GAPs — **minimum set eight, no runnable subset** | **PI**, queue 11 |
| M5 | Bind the crossbridge connector — `BOND`s, partner as `STATE`. **Blocked on ①** | Lead |

⚠ **M2 without M1 is choosing a constant against a wrong edge.** ⚠ **M5 without ① is a population that
cannot transition.**

### 4.2 Integrin and the adhesion chain — the modules exist; the archetype does not

⚠ **This may be a wrong value rather than a missing one.** `INTEGRIN_A5B1` is fibronectin/α5β1 —
fibroblast/mesenchymal. MCF7 is epithelial: basement membrane, laminin, α6β1/α6β4. Under *"start at the
physiological operating point"* the present constant is a candidate defect, not a gap.

| step | what it is | who |
|---|---|---|
| A1 | Rule the archetype — sphere-in-contact or spread. Decides whether adhesion is required at all and which ligand | **PI**, queue 3 |
| A2 | Choose the ligand/receptor pair the archetype implies | **PI** |
| A3 | Bind — not write — the existing modules. Seven of eight already exist | Lead |
| A4 | Substrate boundary condition + FA count | **blocked on citation, see §7** |
| A5 | The clutch itself. **Blocked on ①** | Lead |

⚠ **The adhesome is one spring.** The real chain is integrin–talin–vinculin–actin and **talin is the
mechanosensor** — force unfolds its rod domains and exposes vinculin sites, which is *the* mechanism of
adhesion reinforcement. `fa_maturation.py` exists and is unbound. **A one-spring clutch cannot show
reinforcement, so any stiffness-sensing result taken on it is a result about a spring.** Represent
reinforcement as **state on the `BOND`**, never as new nodes.

### 4.3 G-actin — the one item that is a REPAIR, not an addition

`world/active/protrusion.py:39,61` passes **`G_actin_uM` as a declared scalar** to a law resolver, with
the module itself flagging *"whether that concentration is this cell's"*. The arena's design says a
species pool is a `GRID_CELL` channel; `grid_cell` live is **0**.

**A spatially uniform constant standing in for a diffusing pool is the lumping the charter forbids, and
it is wrong today rather than merely absent.** It is also the one place where a probability cloud is
not an option but the *correct* answer, because a monomer pool has no identity worth tracking.

⚠ **This is the only item on this plan that does not wait on ①**, because a field update belongs to the
outer physical-time loop rather than to a kinetic connector's commit. **Recommended first.**

### 4.4 Laminin / basement membrane — a `STRAND` population

Absent entirely. It is the one molecule that genuinely wants nodes: 565,520 for a whole-cell shell at
50 nm, 7.6% of free headroom. ⚠ **The mechanical argument for it (an intact BM as an insulator against
collagen remodelling) rests on p13, which may not be cited — §7.** The absence is real regardless of
who established its consequence; what is withdrawn is the citation, not the gap.

### 4.5 Recorded as assumptions rather than built

* **No nucleotide.** `polymerization_warp.py` uses ATP-actin on-rates and `crossbridge_kmc.py` models
  load-slowed ADP release; there is no ATP pool, no hydrolysis, no ADP. **Motors run on infinite fuel.**
  For resting mechanics ATP is saturating and a constant is defensible — it becomes a defect the moment
  anything metabolic, hypoxic or blebbistatin-like is claimed, because those act *through* the
  nucleotide cycle. **Today it is recorded nowhere. Record it.**
* **The membrane is a surface, not a bilayer.** No leaflets, no spontaneous curvature, no
  curvature-sensing protein. Only a defect for curvature/budding/scission questions — **but the tree
  states that limit nowhere**, so a future curvature claim would meet nothing standing in its way.
* **Only the bound side of every reaction exists.** Free integrin, free NMII and free crosslinker have
  no representation, so **the binding reaction has one side.** Ca²⁺, septin, spectrin and the
  glycocalyx have none at all.

### 4.6 NMIIB — last, and for a narrower reason than was reported

⚠ **The earlier report that the corpus has no NMIIA-vs-NMIIB source is RETRACTED.** Nagy 2013 is
registered — `SE549`, `source_audit.verdict = OK`, and `crossbridge_kmc.py:209` was citing it.

⚠ **The correction has its own caveat and it is not a softening: SE549's edges = 0.** It is attached to
no `KnowledgeClaim`, `ModelContract` or `Parameter`. **A passing verdict on a row that links to nothing
reads as coverage and is not coverage.** What remains true is that `hand_kmc.py` has no NMIIB preset,
and a duty band is a cross-check quantity, not the `v0` / `F_stall` / `k_off` set a second isoform is
built from.

---

## 5. Order, and the dependency that decides it

```
                    ① acceptance criterion (queue 18)  ─── blocks ALL kinetics
                                  │
        ┌─────────────────────────┼──────────────────────────┐
        │                         │                          │
   M4 crossbridge            A1 archetype (queue 3)     4.4 basement membrane
   PI-GAPs (queue 11)             │                          │
        │                    A2 ligand                       │
   M1 guard ─ M2 radius           │                          │
   (queue 14)                A3 bind 8 modules (queue 22)    │
        │                         │                          │
   M3 re-measure             A4 substrate + FA count ── blocked on citation (queue 23)
        │                         │
   M5 crossbridge bind       A5 clutch bind
        └────────────┬────────────┘
                     ▼
            first contractile cell

   4.3 G-actin field ─── independent of ① ─── DO THIS FIRST
```

**Read the graph as: nothing on it can be falsified before ①.**

## 6. What each step must PROVE

Written before the work, per the charter, so none of them can be scored on what the run happens to keep.

| step | passes when |
|---|---|
| 4.3 G-actin | `grid_cell` live > 0; a declared projection of the field reproduces the old scalar as its spatial mean in a uniform initial condition; the protrusion law reads the field, not a constant |
| M3 | `nearest_cortex_node_um.min` below the capture proxy **and** `p95` below it — a population, not one filament |
| M5 | `bond` live > 0 **and** a positive control that must be REFUSED is refused |
| A3 | the eight modules bound; a Dirichlet/substrate condition present; the cell's basal disc carries a force, not just geometry |
| A5 | a clutch that binds **and unbinds**, committed on an accepted step, with the un-commit path exercised |

⚠ **"It ran" is not on this list, and neither is a green suite.**

## 7. ⚠ Two steps of the audit's own recommendation cannot be executed as written

1. **"taking the FA count and substrate moduli from p50"** — **p50 cannot be cited.** Queue 23 names p13
   and p50 specifically as the two papers fully present in the content layer and **absent from the
   Contract-Graph entirely** (`se_uid = None`); they are *"the two papers this session came within one
   commit of citing in a deliverable"*. **A4 is therefore blocked on queue 23, not merely on the
   archetype.** Either register them, or take the numbers from a source that has a row.
2. **The PI's ratification that a physiological range may be used where no source exists** changes 4.6.
   It does **not** change ①: *a range with no acceptance criterion under it is still a number nothing
   can reject.*

## 8. What may not be claimed from any of this

* Nothing here supports a magnitude. `STATE.md` (c) 3 and (c) 17 stand.
* No contractility result may be taken on the current cell — (c) 20.
* No MT tensegrity or MT↔nucleus result may be taken on the current cell — (c) 19.
* A stiffness-sensing result on a one-spring adhesome is a result about a spring, not about talin.
* Binding a module is not a physics result; it is a precondition for one.

## 9. Decisions this plan is waiting on

| | | blocks |
|---|---|---|
| **queue 18** | the acceptance criterion | everything kinetic — M5, A5, and all of 4.5 |
| **queue 14 + M1** | the NMII shell radius, and the guard upstream of it | M2, M3, M5 |
| **queue 3** | the archetype | A1–A5, and whether adhesion is required at all |
| **queue 11** | the eight crossbridge PI-GAPs | M4, M5 |
| **queue 23** | citability | A4, and the BM argument in 4.4 |

**Recommended immediate action: 4.3 alone.** It is the only item that is a repair of something wrong
today rather than an addition, and the only one whose correctness does not wait on ①.
