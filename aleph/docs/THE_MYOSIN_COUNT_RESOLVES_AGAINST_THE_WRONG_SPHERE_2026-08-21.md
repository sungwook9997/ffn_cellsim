# The myosin count resolves against the wrong sphere — and a self-check pins it there

**2026-08-21, session `d380041d`. NOT FIXED. The number is not changed by this document.**

`build_nmii` takes an areal density and an area, and multiplies them. Four call sites pass an area
that is **not the one the parameter's own docstring asks for**, and the module's self-check asserts
the product that follows from it.

---

## What the parameter says it wants

`aleph/world/build/nmii.py:211`, in `plan_nmii`'s own Args block:

> `support`: the **MEASURED** support `count` resolves against — the **cortical shell area** [µm²] for
> an `areal` count, 1.0 for an `explicit` one. **Measured, never analytic.**

Two requirements in two sentences: *the cortical shell*, and *measured, not analytic*.

## What the four call sites pass

```
aleph/scripts/world_phase1_native.py:178      support=706.86
aleph/scripts/world_render_native.py:110      support=706.86
aleph/scripts/world_export_cell.py:217        support=706.86
aleph/scripts/world_tier0_stations.py:105     support=706.86
```

`4 * pi * 7.5**2 = 706.8583`. **706.86 is the membrane sphere, computed analytically from the nominal
cell radius.** It is the wrong sphere AND it is analytic — it fails both halves of the one sentence
that defines it.

## The two spheres are not the same sphere, and the same call names the other one

Every one of those four calls passes, on a neighbouring line:

```
radius_um=fp.nmii_radius_um          # 6.95 um
```

and `plan_nmii:218` defines that as *"mid-shell radius the minifilaments are placed on"*, while
`geometry.py:68` defines the field as *"the cortical shell radius NMII sits in"*. The module, the
geometry and the call site all agree that the cortical shell is at **6.95 µm**.

So a single call resolves the count against a sphere of 7.5 µm and lays the result on a sphere of
6.95 µm.

```
4*pi*r_cell^2   at 7.50 um  =  706.858 um^2     <- what support is given
4*pi*r_nmii^2   at 6.95 um  =  606.987 um^2     <- the shell the beads are placed on
ratio                        =    1.1645        (+16.5%)

0.625 um^-2  x  706.858  =  441.8  ->  442 minifilaments   (what is built)
0.625 um^-2  x  606.987  =  379.4  ->  379 minifilaments   (the shell reading)
                                        +63 minifilaments, +16.5%
```

⚠ **Which reading is correct is a literature question and is NOT decided here.** If Nie 2015's
0.625 µm⁻² is measured per unit *cell surface*, then 442 is the right total and the shell simply
carries a higher density (442/607 = 0.728 µm⁻²) because it is the smaller sphere. If it is measured
per unit *cortical shell*, then 379 is right and the cell has 63 minifilaments it should not have.
**The defect is not that one of these is wrong. It is that the code does not say which was meant,
and its own docstring asks for the reading it does not use.** That belongs to the PI.

## Why nothing caught it

* `assert_inside_membrane` passes either way — 6.95 + 0.10 half-thickness = 7.05 < 7.5. It asks
  whether nodes are inside the membrane, which they are. It cannot separate two spheres that are
  both inside.
* `assert_partitioned` is a property of ID ranges and says nothing about area.
* The number is **typed**, so there is no derivation to disagree with anything.
* And `nmii.py:445` **asserts it**:

  ```python
  assert plan_nmii(count=nie, support=706.86, **native)["n_minifilaments"] == 442
  ```

  The self-check hard-codes the support the docstring forbids, twelve lines of source below the
  sentence forbidding it, and passes. **A check that ratifies its own subject's error is not a weak
  check; it is the mechanism by which the error became load-bearing.** This is the same shape as the
  other instances logged tonight, and it is the fifth.

## The other half: `support` is typed four times because there is nowhere to read it from

`CellFootprint.__slots__` publishes `r_cell_um` and `nmii_radius_um` and **no area field at all**.
A caller who wanted to obey the docstring has nothing to call. So four call sites each typed a
number, and they agree today — 0.625 / 706.86 / n_bb 14 / 0.301 / 0.200 — by having been copied,
not by having a source.

⚠ This is precisely the defect `aleph/world/populations.py` was written to end. Its own docstring:

> *"The PHASE 1 driver and the renderer each carried their own copy of these calls, with the
> placement arguments written out twice … put **95.5% of the stress fibres and 92.7% of the
> lamellipodium outside the membrane**, while `assert_partitioned` reported success."*

That fix pulled **eight** populations into one builder whose rule is *"nothing is typed here, every
position comes from `fp`"* — and **left NMII outside it, in four copies.** The one population the
consolidation missed is the cell's only active element.

## And PHASE 4 has no copy at all

`world_phase4_native.py` calls `build_all` (3 populations) + `build_remaining_populations` (8) = **11**,
and never calls `build_nmii`. Both landed τ records say so:

```
seed 1: 11 -> [chromatin, cortex, filopodium, intermediate_filament, lamellipodium, lamina,
               membrane, microtubule, nuclear_envelope, sf_arc, stress_fiber]
seed 2: 11 -> (identical)
```

Its own `--core-only` help text says the default is *"the same **twelve** populations the PHASE 1
driver measures"*, and `populations.py`'s title line says *"The **nine** populations that are not
cortex, membrane or envelope"* against a dict with eight keys. **Two documents, each off by one, in
the same direction, and the missing entry in both is NMII.**

⚠ **Every gamma number this session measured was measured on a cell with no myosin in it.** That is
not a correction to those numbers — a passive cortex is a legitimate thing to measure and the records
list what they contain. It is a statement about what they may be compared to.

⚠ And a record that lists **eleven names it has** does not tell a reader it lacks a twelfth. You have
to already know the twelve. That is the same shape as the exporter's `populations_skipped_empty`,
fixed earlier tonight by adding the list of what is *absent* beside the list of what is present.

## What is not done here

The count is **not changed**, the docstrings are **not corrected**, and the four call sites are
**not consolidated**. Changing a population count is a change to the cell, and this session is
holding a running τ series it commissioned. The relevant rule is the one that has held all night:
the person holding the result does not move the quantity the result rests on.

**Filed to the PI.** The literature reading (per cell surface, or per cortical shell) is the decision;
the consolidation into `populations.py` and the two off-by-one docstrings follow from it.
