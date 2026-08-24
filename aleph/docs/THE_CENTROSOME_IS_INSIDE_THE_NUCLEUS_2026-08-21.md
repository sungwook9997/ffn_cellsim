# The centrosome is inside the nucleus, and two thirds of the microtubules are in there with it

**2026-08-21, session `d380041d`. NOT FIXED. No position is changed by this document.**

Found by rendering the cell and then measuring what the picture suggested, rather than by reading
code. The picture showed the nuclear envelope as a thin crescent beside a microtubule aster instead
of the hemisphere a centred sphere should cut to. **Eyes are not evidence, so it was measured.**

## The measurement

`aleph/outputs/ac/world_phase4/tau3_seed2.alephcell`, frame 0, radii about the world origin [µm]:

```
population                        n    r_min    r_med    r_max   centroid (x,y,z)
microtubule                  74,001    0.000    3.700    7.400   (-0.000, -0.000, +0.000)
nuclear_envelope             40,962    5.100    5.100    5.100   (+0.000, +0.000, +0.000)
lamina                        2,043    5.100    5.100    5.100   (+0.000, -0.000, +0.000)
intermediate_filament         2,460    5.250    6.250    7.250   (-0.007, -0.002, +0.000)
cortex                    4,197,654    7.300    7.399    7.500   (+0.011, +0.015, +0.035)
membrane                    163,842    7.500    7.500    7.500   (+0.000, +0.000, +0.000)
```

The nuclear envelope is an exact sphere: **R = 5.100 µm, min = max = median, centred on the origin.**
The microtubule aster starts at **r = 0.000** — the origin — and runs out to 7.400.

```
microtubule nodes inside the nuclear envelope:  50,735 / 74,001  =  68.6%
```

## What the code says the origin is

`aleph/world/build/microtubule.py:203`, in `build_microtubules`' own Args block:

> `centre`: **the MTOC position** [µm].

`aleph/world/populations.py:104`:

```python
out["microtubule"] = run("microtubule", lambda: build_microtubules(
    arena, count=gap(500, "microtubules"), centre=(0.0, 0.0, 0.0), reach_R_um=7.4, seg_um=seg_um))
```

**The MTOC — the centrosome — is placed at (0, 0, 0), which is the geometric centre of the nucleus.**
A centrosome is a cytoplasmic organelle. It is not inside the nucleus in any cell, in any state, at
any point in interphase.

## And the law module states the constraint that is being violated

`aleph/laws/microtubule.py:38-40`, as a statement of physical fact, not a to-do:

> A microtubule is NOT a zero-thickness line: the tube radius is its steric footprint (**it cannot
> interpenetrate the nucleus** or other filaments) and its rendered thickness.

The same module, line 6: *"`n_mt` tubes radiate from ONE MTOC (centrosome)"*, and line 10: the aster
*"positions the nucleus"*. **A centrosome at the nucleus's centre cannot position it** — there is no
direction for it to push, and every arm's first 5.1 µm is inside the thing it is supposed to move.

⚠ This is the same family as the myosin `support` filed as item 19, but a step worse. There, a
docstring asked for one area and the caller passed another — a bookkeeping disagreement. **Here the
law module names the constraint, the builder's own Args block names the organelle, and the built
geometry violates both.** The claim and the code are not in tension about a number; they are in
tension about whether two objects can occupy the same space.

## The line above it got this right

`populations.py:107`, the very next call:

```python
out["intermediate_filament"] = run(..., build_intermediate_filaments(
    arena, count=gap(60, "IF spokes"), centre=(0.0, 0.0, 0.0), R_nuc_um=5.1, R_cortex_um=7.4, ...))
```

`build_intermediate_filaments` **takes `R_nuc_um`**, and its spokes measure out at r_min = 5.250 —
outside the envelope, correctly. `build_microtubules` has **no nucleus argument at all**: `centre`,
`reach_R_um`, `seg_um`, and nothing that could know a nucleus exists.

**Two adjacent lines in one function. One builder is nucleus-aware and the other cannot be.**

## Why nothing caught it

* `assert_inside_membrane` asks whether every node is inside R = 7.5 µm. Every MT node is: 7.400 < 7.500.
  **It passes, correctly, on the question it asks.** It has no concept of an interior exclusion.
* `assert_partitioned` is about ID ranges.
* No test compares `microtubule` to `nuclear_envelope`. The only test naming both
  (`test_families_census.py:119`) checks that a set of population names was built.
* Nothing in the engine evaluates MT–nucleus sterics, so the interpenetration costs no energy and
  produces no force. **It is silent by construction.**

⚠ **And the same aster was audited from the other end and passed.**
`aleph/docs/v2_audit/INTERNAL_DISPLACEMENT_DIAGNOSIS_2026-07-16.md:35`:

> **MT tips don't reach the cortex** — 6 µm MTs from a **central MTOC** end 0.9 µm short of the 7.5 µm
> cortex, so the tip↔cortex contact never engages … → MTs can't act as compression struts.

That audit looked at the aster, found the **outer** end wrong, fixed it — `reach_R_um` is 7.4 today,
and the tips do reach — and wrote *"a central MTOC"* in passing without the root end ever being
asked about. **The audit that touched this geometry recorded the premise and examined only one end
of it.**

## What this does and does not affect

⚠ **It does not invalidate tonight's γ numbers.** γ is measured on the cortex by method of planes;
microtubules carry no bound law in these runs and contribute nothing to it. The runs are what they
are.

⚠ **It does bear directly on tensegrity.** The whole point of the MT compartment
(`laws/microtubule.py:10`) is to be the compression element that balances actomyosin tension and
positions the nucleus. **That mechanism cannot be studied on this geometry** — not because the
numbers would be off, but because the strut's base is inside the object it is supposed to brace.
Any future MT–nucleus or MT–cortex coupling result on this cell would be a result about a
configuration that does not occur.

## What is not done here

No position is changed. Moving the MTOC changes the built cell, and this session is holding a τ
series that rests on that cell — the same reason item 19 was filed rather than fixed. It is also
not obvious what the correct placement is: a centrosome sits adjacent to the nuclear envelope, which
means an **off-centre** MTOC and therefore an anisotropic aster, and that is a modelling decision
about the archetype, not a coordinate to type in. **It connects to PI queue item 3** — sphere-in-contact
or spread cell — because where the centrosome sits is defined relative to a polarity the archetype
has not yet fixed.

**Filed to the PI as item 20.**

---

## What else was measured the same way, and holds

The same radii pass was run over every other population against the placement claims the record
makes, so that nobody repeats it. **These are negative results and are recorded as such.**

Contact disc: radius **2.6926 µm** at z = −7.0, `r_footprint_um` DERIVED as `sqrt(R_cell² − z_basal²)`.

```
population         n       rho_min  rho_med  rho_max    z_min    z_max   outside disc
stress_fiber       9,240     0.025    1.367    2.227   -7.016   -6.984   0  (0.0%)
lamellipodium      4,200     0.018    0.832    2.038   -7.000   -7.000   0  (0.0%)
sf_arc             2,600     0.574    1.801    3.028   -6.616   -6.584   520 (20.0%)
filopodium        61,000     0.411    3.629    6.269   +2.127   +7.215   45,521 (74.6%)
```

* **`stress_fiber` — the claim holds.** `populations.py` says `sf_length_um` is *"a side of the
  inscribed square, so both FA ends land on the contact disc"*. Every one of 9,240 nodes is inside
  the disc, and z spans −7.016 to −6.984 — the basal plane, ±16 nm.
* **`lamellipodium` — holds.** Every node inside the disc, exactly on z = −7.000.
* **`sf_arc` — the 20% is MY comparison being wrong, not a defect.** Arcs are lifted `lift_um = 0.4`
  and sit at z ≈ −6.6, where the membrane's local radius is `sqrt(7.5² − 6.6²) = 3.562 µm`, not 2.693.
  Their rho_max of 3.028 is inside that. **The contact disc is the wrong ruler for a population that
  is deliberately off the contact plane.** Recorded because the wrong ruler looked like a finding.
* **`filopodium` — 74.6% "outside" is also the wrong ruler, and the record already said so.** They are
  apical (z > 0), so the basal disc does not apply. And `placement_envelope` carries, as a field:

  > `filopodium_note`: *"rooted so the tip REACHES the membrane and stops. A real filopodium protrudes
  > and the membrane follows it; this membrane is a rigid icosphere, so protrusion and 'inside the
  > membrane' cannot both hold. **PHASE 2 property, recorded rather than silently traded away.**"*

⚠ **That last one is worth stating on its own.** This document and items 19 and 20 are a catalogue of
places where a record listed what it had and not what it lacked. `filopodium_note` is the same
codebase doing the opposite: an approximation that could not be avoided, named **inside the record**,
in a field that travels with the data, in the record's own words. It is the pattern the rest of
tonight's findings are missing, and it was already here.
