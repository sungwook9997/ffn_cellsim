# Adversarial audit — nodes, smaller nodes, or a probability cloud?

**2026-08-21, session `d380041d`.** PI question: *"노드를 추가해서 이런 단백질들에 의미를 부여하는 것이
맞는지, 아니면 더 작은 노드를 만들어야 하는지, 아니면 확률 구름이나 다른 방식으로 구현해야 되는지 빠진
부분 모두 적대적 감사해서 찾아줘."*

**Audited against the built cell and the arena's own source, not answered from memory.** Every number
below is measured.

---

## 0. The arena already answered this, and the answer is written in its rejection argument

`aleph/world/arena.py:76-78`, in the paragraph that closed the primitive set:

> `HEAD` is a **NODE plus a BOND** whose far address is a point on a SEGMENT; `SITE` is an **address
> resolved from live geometry, not an allocation**; and a species pool is a **named channel on a
> `GRID_CELL` field** plus one conserved scalar per consuming range.

Three candidate primitives were proposed and **rejected because each decomposes into existing ones** —
and those three are exactly the PI's three options:

| the PI's option | the arena's answer | status |
|---|---|---|
| 노드를 추가 (add nodes) | for a **molecular attachment**: `BOND`, far address on a `SEGMENT`. Not a node. | **0 of 1,000,000 claimed** |
| 더 작은 노드 (smaller nodes) | a binding **site is not an allocation at all** — resolved from live geometry | never exercised |
| 확률 구름 (probability cloud) | a **species pool is a `GRID_CELL` channel** + one conserved scalar per consuming range | **0 of 4,000,000 claimed** |

⚠ **So the design is right and it is entirely unused.** The two primitives these proteins need are at
**zero** in the cell that produced every number this session measured.

---

## 1. A node is 50 nm, which changes which question is the right one

Measured on the built cell — median bond length per population:

```
cortex           50.00 nm    (n = 4,128,840)
microtubule      50.00 nm    (n =    73,500)
stress_fiber     49.80 nm    (n =     9,100)
nmii            200.00 nm    (head arm)
```

An actin monomer is **5.4 nm**, so **one node is ~9 monomers**. This is not a monomer-resolved model
and never was. Against that ruler:

| protein | size | vs a 50 nm node |
|---|---|---|
| integrin (ectodomain) | 10–20 nm | **sub-node** |
| α-actinin | ~35 nm | **sub-node** |
| filamin | ~150 nm | 3 nodes |
| laminin (cross arms) | 35–77 nm | ~1 node |
| NMII minifilament | ~300 nm | 6 nodes |
| focal adhesion | 1–2 µm | 20–40 nodes |

⚠ **"Smaller nodes" is the wrong axis, and giving integrin its own node is actively wrong.** It would
resolve a 15 nm receptor more finely than the 50 nm actin it binds to — an asymmetric resolution where
the finer object's detail cannot be felt by the coarser one it acts on. **These molecules are links
between things, not extended things.** A link is a `BOND`.

The one that genuinely wants nodes is **laminin**, because a basement membrane is a *network with its
own mechanics* — the thing a `STRAND` population is for.

---

## 2. Cost, so the argument cannot hide behind expense

```
node headroom   12,000,000 capacity − 4,558,554 live = 7,441,446 free (62%)
bond headroom    1,000,000 capacity −         0 live = 1,000,000 free (100%)

BASEMENT MEMBRANE (laminin-332) at 50 nm, 100 nm thick
  basal contact disc as built (22.8 µm²)          18,240 nodes
  max contact disc for R = 7.5   (177 µm²)       141,600 nodes
  whole cell surface             (706.9 µm²)     565,520 nodes

INTEGRIN CLUTCHES as BONDs  (p50: ~37 nascent FAs, min 6, max 89)
  37 FAs × 500                                    18,500 bonds
  37 FAs × 1,000                                  37,000 bonds
  89 FAs × 2,000                                 178,000 bonds
```

⚠ **Neither is expensive.** A whole-cell basement membrane is 0.6 M nodes against 7.4 M free, and every
integrin in the cell is under 20% of the bond capacity. **The obstacle is not memory, not count, and
not the node budget.** Any argument that reaches for cost here is reaching for the wrong reason.

---

## 3. ⚠ The finding that outranks all of the above

`CLAUDE.md`: *"a kinetic connector commits **ONLY on an accepted physical step** under one
device-resident transaction."*

`aleph/world/bond.py:27`: *"…detaches yet, so a free list would be machinery for a case that does not
exist; **no transaction runs yet**…"*

`aleph/world/step.py:35`: *"the whole transaction sequence executes — and the verdict recorded is
`ACCEPTANCE_UNDEFINED`."* **Zero steps accepted, zero rejected, in every run this session measured.**

**So every chemistry in this design is unreachable by construction.** An integrin clutch that binds and
unbinds is a kinetic connector; a kinetic connector may only commit on an accepted step; no step has
ever been accepted. Adding integrin today produces a population that **cannot transition**, and its
correctness cannot be falsified because nothing it would do is permitted to commit.

⚠ **This is the adversarial answer to the question as asked.** The representation question has a
correct answer already written in `arena.py`, the cost is negligible, the law modules are written —
and **none of that is the blocker.** PI queue **18** is.

---

## 4. What is missing — everything found, by severity

### (a) BLOCKING — no acceptance criterion
Above. Until item 18, any binding kinetics is machinery nothing can accept, reject or falsify.

### (b) BLOCKING — the cell is attached to nothing
`world/` binds **11 of 46** law modules and not one is adhesion or ECM (item 22). No integrin, no FA,
no substrate, no Dirichlet condition. The contact disc is geometry with no physics under it.

### (c) ⚠ The free pool is a **declared scalar**, not a field — this is a live lumping
`world/active/protrusion.py:39,61` takes **`G_actin_uM` as a declared concentration** passed to a law's
resolver, with the module itself flagging *"whether that concentration is this [cell's]"*. The arena's
own design says a species pool is a `GRID_CELL` channel; **`grid_cell` live is 0.**

**A spatially uniform constant standing in for a diffusing pool is exactly the lumping the charter
forbids** — and it is the one place where "확률 구름" is not an option but the *correct* answer, because
a monomer pool genuinely has no identity worth tracking. **Filed as the strongest case for the field
representation, and it is already wrong today rather than merely absent.**

### (d) No nucleotide accounting
`polymerization_warp.py:15` uses ATP-actin on-rates (Pollard 1986) and `crossbridge_kmc.py:5` models
load-slowed **ADP release** (Thirumurugan 2007). Both consume a nucleotide state that **is not
represented**: there is no ATP pool, no hydrolysis, no ADP. Motors and polymerisation therefore run on
infinite fuel.
⚠ **Whether that matters is a question, not a defect** — for resting mechanics ATP is saturating and a
constant is defensible. It becomes a defect the moment anything metabolic, hypoxic or
blebbistatin-like is claimed, because those act *through* the nucleotide cycle. **It should be recorded
as an assumption, and today it is recorded nowhere.**

### (e) No basement membrane object
Absent entirely. p13 makes it **mechanical rather than decorative** — intact laminin-332 BM is a
*mechanical insulator* that prevents the cell from remodelling collagen, and breached BM lets invasion
follow — so **that switch cannot be expressed at all in this tree.**

🔴 ⚠ **But p13 may not be cited, on two independent grounds, and I quoted it above before checking
either.** Session 75's DOI mapping: p13 is `10.1101/2025.01.31.635980`, **unregistered** — there is no
`SourceEvidence` row to carry a verdict. And it is a **bioRxiv preprint, NOT peer reviewed**, which
registration would not repair. **The dossier told me so at `p13.md:219`** — *"flag verdict accordingly
before citing in a deliverable"* — in the same file I read the claim out of. I grepped for `laminin`
and read the lines that answered me, not the line that said whether I could use the answer.

**The gap is real either way**: a tree with no basement membrane cannot express a BM switch regardless
of who established it. What is retracted is the citation, not the absence.

### (f) The membrane is a surface, not a bilayer
163,842 nodes / 409,600 faces of triangulated `FACE`. No leaflets, no lipid, no spontaneous curvature,
no BAR-domain or curvature-sensing protein. ⚠ A defect only for questions about curvature generation,
budding or scission — **but the tree nowhere states that limit**, so a future curvature claim would
have nothing standing in its way.

### (g) The adhesome chain is a single spring
`INTEGRIN_A5B1` is one catch-slip bond. The real chain is integrin–talin–vinculin–actin, and **talin is
the mechanosensor**: force unfolds its rod domains and exposes vinculin sites, which is *the* mechanism
of adhesion reinforcement. `fa_maturation.py` exists and is unbound. ⚠ **A one-spring clutch cannot
show reinforcement**, so any stiffness-sensing result on it would be a result about a spring.

### (h) Wrong ligand for the archetype
`INTEGRIN_A5B1` is fibronectin/α5β1 — fibroblast/mesenchymal. MCF7 is epithelial: basement membrane,
laminin, α6β1/α6β4. Under *"start at the physiological operating point"* this is not a missing addition
but a **possibly incorrect present value** (item 22).

### (i) ⚠ NMIIB — I TOLD THE PI THIS HAD NO SOURCE AND I WAS WRONG

**Nagy 2013 is a registered source with a passing audit.** Session 75 built `kb.duckdb` and queried it:

```
SE549   Nagy2013_JBC   source_type "Direct measurement"   anchor_status OK
        source_audit.verdict = OK
        DOI 10.1074/jbc.M112.424671   -- the same DOI params_i0b3.yaml's comment carries
```

**My report to the PI that the corpus has no NMIIA-vs-NMIIB source was incorrect.** What I searched was
the `docs/corpus/` dossiers, and I concluded from three near-miss hits there that nothing existed —
**a conclusion about the corpus drawn from a search of one layer of it.** The Contract-Graph had the
row all along, and `crossbridge_kmc.py:209` was citing it.

⚠ **And the correction has its own caveat, which is not a softening.** Session 75 also reports
**SE549's edges = 0**: it is attached to no `KnowledgeClaim`, `ModelContract` or `Parameter`. It is an
orphan row with a passing verdict. **A verdict of OK on a row that links to nothing reads as coverage
and is not coverage** — the same shape as everything else logged this session. So NMIIB has a source,
and that source supports nothing yet.

⚠ **What is still true**: `hand_kmc.py`'s roster has no NMIIB *preset*, and a duty band (0.2–0.3) is a
cross-check quantity, not the v0 / F_stall / k_off set a second isoform needs. **The gap is narrower
than "no source" and wider than "one row away".**

### (j) Species with no representation at all
Ca²⁺ and signalling; septin and spectrin; the glycocalyx; free pools of integrin, NMII and crosslinker
(the bound ones exist, the unbound ones do not, so **the binding reaction has only one side**).

---

## 5. Recommendation

**The representation question does not need a decision — `arena.py` already made it.**

```
laminin / basement membrane   ->  STRAND population (a network with its own mechanics)
integrin clutch               ->  BOND, far address on a SEGMENT, partner as STATE
                                  (the pattern nmii_cortex_crossbridge.py:31 already implements)
talin/vinculin reinforcement  ->  state on that BOND, not new nodes
free G-actin / free integrin  ->  GRID_CELL channel + one conserved scalar per consuming range
binding site                  ->  not an allocation; resolved from live geometry
```

**What needs a decision is the order**, and the audit says it is not what it looks like:

1. **Item 18 — an acceptance criterion.** Without it every kinetic connector is un-committable and
   nothing added can be falsified. **This is the whole blocker.**
2. **Item 3 — the archetype.** Decides whether adhesion is required and which ligand.
3. **Then bind**, laminin/α6 rather than fibronectin/α5β1, taking the FA count and substrate moduli
   from p50 — **binding, not writing**: seven modules already exist.
4. **G-actin as a field** can go earlier than the rest, because it is not an addition but the
   **repair of a lumping that is live today**.
5. **NMIIB last** — session 75 has settled it: Nagy 2013 IS registered (SE549, verdict OK), so my
   "no source" claim is retracted. It remains last because SE549 has **zero edges** and a duty band is
   a cross-check, not the v0 / F_stall / k_off set a second isoform is built from.

⚠ **On the PI's instruction that a physiological range may be used where no source exists** — that is
ratified and it changes step 5, but it does not change step 1. A range with no acceptance criterion
under it is still a number nothing can reject.
