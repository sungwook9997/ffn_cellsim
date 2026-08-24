# ALEPH-PORT-3643 — the filament is a dotted line, and refining it would make that worse

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3643` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-06` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **INTERNAL FIX.** Segment–segment geometry, re-derived. Nothing is read from a provider for it. |
| Aleph target | `aleph/vertical/cortex_filaments.py`, `aleph/vertical/cortex_surface_coupling.py` |
| Depends on | `ALEPH-PORT-3640` (which needs the same primitive) |
| Blocks | the filament node count becoming a parameter — see §3 |
| Exists because | Two cortical filaments can pass **through** each other and feel exactly zero force. |

---

## 0. The measurement

Two three-node filaments, 300 nm long, crossing at 90°, one translated so the crossing falls
**between** the other's nodes. `σ = 7 nm`, WCA cutoff `r_c = 7.86 nm`:

| gap at the crossing | min node–node distance | pairs inside `r_c` | \|F_steric\| |
|---:|---:|---:|---:|
| 20 nm | 77.62 nm | 0 | **0 pN** |
| 5 nm | 75.17 nm | 0 | **0 pN** |
| **0 nm — the rods intersect** | 75.00 nm | **0** | **0 pN** |

**They occupy the same point in space and nothing resists it.**

## 1. Why: a bead law on a mesh 21× its own bead

`steric_energy_and_forces(positions, pairs, sigma_um, epsilon_pn_um)` is a **node-pair** law:
`U(r) = 4ε[(σ/r)¹² − (σ/r)⁶] + ε` between node *centres*, zero beyond `r_c`.

That is a bead model, and a bead model is a rod only when the beads are spaced about their own
diameter apart. Here:

```
node spacing   150 nm          (a 300 nm filament with three nodes)
sigma            7 nm          (F-actin diameter, sourced)
                 -> the beads are 21x further apart than they are wide
```

So the filament is not a cylinder. **It is a dotted line**, and every crossing that lands between
nodes — which at this spacing is nearly all of them — passes through unopposed.

### This invalidates a gate that has already been used

`ALEPH-PORT-3636` resolves the shell to "zero cross-filament pairs inside `r_c`" and grades on it.
That gate cannot distinguish *no overlap* from *no node near an overlap*, and this measurement shows
the two are different at a 150 nm mesh. The overlap resolution is still worth having — it removes the
interpenetrations it can see, and the `K_A` blow-ups it removed were real — but **"0 pairs" is a
weaker statement than it has been read as**, and every result citing it inherits that.

That is recorded here rather than filed against 3636, because the fix is this entry's and 3636's
gate was correct about the quantity it named.

## 2. The fix: segment–segment, which is also what the crosslink search wants

The excluded-volume interaction of two rods is a **line–line** interaction. The right discrete form
computes the closest approach between two *segments* and applies the potential there — once per
segment pair, at the point of closest approach, with the force split between each segment's two
endpoints by the interpolation weights.

Two properties follow, and the second is the reason this blocks §3:

* **It is a rod at any mesh.** Two crossing filaments interact whether or not a node happens to sit
  near the crossing.
* **It is invariant under refinement.** A node-pair law on a filament refined from 3 nodes to `N`
  gains `O(N²)` pairs between any two filaments and its steric force grows with the discretisation —
  a stiffening with no physics in it. A segment-pair law converges instead of diverging.

`_closest_approach_crosslinks` already needs exactly this primitive and currently approximates it on
two straight sticks. **One geometry serves both**, which is why they are one entry.

## 3. Why this blocks the node count, and not the other way round

`ALEPH-PORT-3638` §C wants `Filament(..., 3, ...)` to become a parameter. **Refining the mesh before
this lands makes the steric law worse, not better**: the dotted line gets denser dots, the `O(N²)`
pair count inflates the force, and the two errors do not cancel.

So the order is fixed by the physics, not by preference:

```
1. segment-segment steric      <- this entry
2. the node count becomes a parameter
3. re-measure everything that used "0 steric pairs" as evidence
```

## 4. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **F-1** | **Two crossing rods repel.** The configuration in §0 at gap 0, with the crossing between nodes. | \|F\| **> 0**, and rising as the gap closes |
| **F-2** | **Refinement invariance.** The steric energy of one fixed crossing configuration at `n` = 3, 5, 9, 17, 33. | varies `< 5%` across `n` |
| **F-3** | **The node-pair law is exactly recovered when the mesh is fine.** At a node spacing below `σ`, segment and node forms agree. | within 1% |
| **F-4** | **Still exactly zero beyond the cutoff.** Two rods at closest approach `> r_c`. | exactly `0.0`, testable with `==` |
| **F-5** | **The gradient is one.** `F = −∂E/∂x` against a central difference. | `≤ 1e-6` relative |
| **F-6** | **The retained control survives.** Default arguments reproduce `K_A = 1081.819088` at ρ = 0.5, level 2. | bitwise |
| **F-7** | **How much did the cortex actually have?** Cross-filament segment pairs inside `r_c` on a shell 3636 declared clean, at ρ = 8, 24, 100. | **measured, no threshold** |

**F-7 is the one that can hurt and has no threshold on purpose.** If a shell the overlap resolution
called clean turns out to carry many *segment* overlaps, then every modulus measured on it was
measured on a structure with interpenetration, and the size of that is a result rather than a thing
to predict. This lane has written a threshold against a guess five times today and been wrong.

**F-2 can fail honestly too.** If the segment form is not refinement-invariant either, the potential
itself is wrong for rods — a Kihara or a proper cylinder excluded volume would be needed — and that
is a different entry.

## 5. What is NOT claimed

- **Not that the cortex's published moduli are wrong by a known amount.** F-7 measures how much
  interpenetration was invisible; until it runs, "some" is the honest word.
- **Not that this makes the filament a cylinder.** A segment–segment WCA is a line with a soft core,
  not a solid rod; it has no torque about its own axis and no finite thickness in the bending law.
- **Not that the node count is addressed here.** It is unblocked, not changed.
- **Not that `ε = 50 pN/µm` is right.** `CortexCard.steric_contact_stiffness_pn_per_um` is row 23 of
  the axis registry, class `axis`, and is marked **UNSOURCED** there. This entry does not source it.

---

## Amendment, 2026-08-06 — F-7, and it is the worst thing measured today

### The gates

| # | outcome |
|---|---|
| F-1 crossing rods repel | **PASS** — 2.72e20 pN at a 0.5 nm gap where the node-pair form gives exactly 0 at every gap |
| F-2 refinement invariance `< 5%` | **PASS** — **0.000%**: 159,306.42 pN·µm at `n` = 3, 5, 9, 17, 33 while the segment-pair count rises 4 → 1,024 |
| F-4 exactly zero beyond the cutoff | **PASS** — `0.0`, testable with `==` |
| F-5 the force is `−∂E/∂x` | **PASS** — 4.158e-10 relative |
| F-6 the retained control | **PASS** — `K_A` 1081.819088, bitwise; the default `steric_law` is `NODE_PAIR` |
| **F-7 what the cortex actually had** | **measured, and it is bad** |

### F-7

Shells built with the overlap resolution of `ALEPH-PORT-3636` run to completion — the ones that
report **zero** cross-filament pairs and were graded clean:

| ρ | filaments | node pairs inside `r_c` | **segment pairs inside `r_c`** | \|F\| node-pair | \|F\| segment |
|---:|---:|---:|---:|---:|---:|
| 8 | 2,493 | **0** | **17** | 0 | 3.56e+17 |
| 24 | 7,480 | **0** | **286** | 0 | 7.71e+32 |

**A shell certified clean carries 286 interpenetrations at ρ = 24, and the node-pair law reports
none of them.** The force magnitudes are the tell: `r⁻¹²` reaching 1e32 means the separation is a
small fraction of σ, so these filaments are not brushing past each other, they are **through** each
other.

### What this does and does not invalidate

**Every cortex modulus this lane measured was measured on a shell with interpenetration.** The
size of the effect is not established by this measurement — a 1e32 pN force is a statement about how
close the rods are, not about how much stress the structure would carry if they were separated —
and quantifying it needs the next paragraph's fix, not more argument.

`ALEPH-PORT-3636`'s gate is not withdrawn. It measured what it named, and the `K_A` blow-ups it
removed were real. What is now on the record is that **"zero steric pairs" was a weaker statement
than every document citing it has read it as**, including this lane's.

### The same defect has a second site

`resolve_cross_filament_overlaps` is **also** node-pair: it separates nodes that are within `r_c` of
each other and cannot see a crossing between them. So the resolution and the law shared one blind
spot, which is why the shell it produces passes its own gate. Fixing the law without fixing the
resolution leaves a shell that is now *measured* as interpenetrating and still built that way.

That is the next entry, not this one. It is named here so the gap is on the record rather than
discovered again.

### Status

`StericLaw.SEGMENT_PAIR` is opt-in and `NODE_PAIR` remains the default, so no committed number
moves. **Nothing in `aleph/` selects the segment law yet**, and until the resolution is fixed too,
selecting it would give a shell whose steric force is 1e32 pN at construction — correct, and
unusable.
