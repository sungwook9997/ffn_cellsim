# A sphere of R = 7.5 µm cannot be an adherent cell, and every basal population inherits that

**Status:** FINDING. Arithmetic, not a run — but it was found by LOOKING at a render, which is the
third time in two days that a picture has caught something the numbers were not being asked.

**What is claimed:** that the current cell geometry cannot produce a physiological basal contact area,
and that every basal population is therefore sized by a number nobody sourced. **What is NOT claimed:**
any spread-area value for MCF7, or that the sphere is the wrong archetype. Both of those are the PI's.

---

## 1. The arithmetic

The cell is a sphere of `R_cell = 7.5 µm` and the substrate plane sits at `z_basal = -7.0 µm`, so the
basal contact disc has radius `sqrt(R² − z²) = 2.693 µm`:

| | value |
|---|---:|
| footprint radius | 2.693 µm |
| **footprint AREA** | **22.8 µm²** |
| nucleus radius | 5.100 µm |
| nucleus projected area | 81.7 µm² |
| **footprint ÷ nucleus** | **0.28** |

**The entire adherent contact region is smaller than the cell's own nucleus** — about a quarter of its
projected area. That is not what an adherent cell looks like from below, and it is visible immediately
in the `basal` view of the native viewer: seven ventral fibres and a few transverse arcs, dwarfed by
the envelope behind them.

## 2. ⚠ And no `z_basal` fixes it

Spread areas for adherent cells are in the hundreds of µm². Solving `π(R² − z²) = A` for the plane:

| target contact area | required footprint radius | possible on R = 7.5? |
|---:|---:|---|
| 300 µm² | 9.77 µm | **no** |
| 700 µm² | 14.93 µm | **no** |
| 1200 µm² | 19.54 µm | **no** |

**A sphere of radius 7.5 µm has a maximum possible contact disc of π·7.5² = 177 µm², and reaching even
that means `z_basal = 0` — the cell cut exactly in half.** There is no plane that makes this geometry
adherent. The constraint is the RADIUS, not the plane.

## 3. What that makes the basal populations

`world/geometry.py` derives the ventral-fibre count, the arc count and the lamellipodium sheet from
this footprint. All three are therefore correct arithmetic on a premise that does not hold:

* `n_stress_fibres = 7`, derived at a held physiological pitch — from a 5.4 µm chord.
* `sf_arc` = 5 arcs (session D), from the same annulus.
* the lamellipodium sheet is 3.2 µm wide, which is **already** flagged as 2.5–10× under the sourced
  leading-edge band (PI decision 11), and this is why: the band cannot be met inside this footprint.

⚠ **So decision 11 is not independently satisfiable.** The PI set the leading edge to a 5 µm test
point; the footprint admits at most 5.4 µm of chord, and only through its centre. The two constraints
touch, and nobody had put them on the same page — which is the same shape as `dx_um` and `k_xl`.

## 4. What is actually undeclared

`z_basal = -7.0` entered in `world/geometry.py` on 2026-08-21 as a default argument, written by this
session while fixing the placement defect. **It has no source and no axis.** It was chosen so a
footprint would exist at all, which is exactly the class of value this project spends its time
catching — and it was written by the session that had just finished writing the module docstring about
that class.

It is not a magnitude error: any `z_basal` gives a footprint too small, per §2. **It is an
undeclared parameter standing in for an unasked question**, and the question is §5.

## 5. For the PI

**Is the adherent archetype a sphere resting on a plane, or a SPREAD cell?**

They are different geometries, not different parameters:

* **a sphere on a plane** is a cell in suspension that happens to touch. Its basal populations are
  small because there is nowhere for them to be, and every count derived from that footprint is
  honest arithmetic about a cell that is not adherent.
* **a spread cell** has a contact area in the hundreds of µm², a flattened dome, and a lamellipodium
  at a leading edge that is a real distance from the nucleus. Its membrane is not an icosphere, and
  `build/membrane.py` builds an icosphere.

⚠ **This does not block tonight's results.** Every measured row stands: the placement invariant, the
step timing and the balance-gate finding are all properties of the machinery and none depends on the
cell being adherent. What it blocks is any basal-population COUNT being read as physiological, and
`world/geometry.py` now says so.

**Recommended, not decided:** declare `z_basal` an axis with its scope stated as *"a sphere in
contact, NOT a spread cell"*, so the current geometry stops implying something it cannot deliver —
and treat the spread archetype as a separate build, because it is one.


---

## 6. ⚠ The sphere assumption survives in exactly FOUR places, and one decision touches all of them

Reported by session D on 2026-08-21 while closing its lane, and verified here. D's point is the one
that makes this item actionable: **"the geometry follows automatically" is not true.** Most of it
does — chord length, arc count, the pitch axis — but four lines encode *sphere* rather than *cell*,
and they are the kind that goes silently wrong the moment `r_cell_um` stops meaning a radius.

| site | line | what it assumes |
|---|---|---|
| `world/build/sf_arc.py` | `r_span = sqrt(fp.r_cell_um² − z_arc²)` | the cell's cross-section at the arc plane is a circle |
| `world/build/sf_arc.py` | `if worst_r >= fp.r_cell_um` | inside means within a radius |
| `world/geometry.py` | `n_out = (norm(pos) > r_cell_um).sum()` | **the PI rule's own checker** |
| `world/build/__init__.py` | `r = norm(pos − centre)` | same shape, in the shell placement |

⚠ **The third is the checker for the rule this document exists to serve.** `assert_inside_membrane`
enforces *"every node lies inside the built membrane"* by comparing a radius, which is exactly right
for an icosphere and means nothing for a spread cell. **A decision that changes the archetype without
changing that line leaves the invariant reading true while measuring nothing** — which is the same
failure as `assert_partitioned` passing on ID ranges while 95.5% of a population sat outside the cell.

**D did not introduce this and neither did the placement fix.** It is a shared assumption that was
correct while the only cell was a sphere, and it is the same shape as `bond.py`'s `node_j` being
topology while every consumer was a static family: right at the time, wrong when a new case arrives.

⚠ **The arc lift is not incidental to it.** Raising a transverse arc off the basal plane is one of the
four properties that distinguish it from a ventral fibre, so `sf_arc` *needs* a cross-section radius at
its own plane. Whatever replaces the sphere has to answer that question, not just remove it.

**So the PI question in §5 carries a fourth line:** a ruling toward the spread archetype is a change to
four call sites, one of which is the invariant checker, and none of them will fail loudly when the
premise moves.
