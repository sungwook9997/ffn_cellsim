# ALEPH-PORT-3662 — a blob can be spread by a rotation, and a rotation costs nothing

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3662` |
| Lane | `b4fad06b` |
| Status | `AUDITED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Answers | defect **D1** — *"six compartments cover 0.0–7.2 % of the sky; `_place_in_band` builds each assembly in the first octant and translates it, and a translation moves a blob without making it not a blob"* |
| Follows | `ALEPH-PORT-3661`, which measured that raising an element count 161× moves the coverage from 6.04 % to 5.61 % |

---

## 1. Aleph API

`aleph/scenarios/whole_cell.py`:

- `COMPARTMENT_SPREAD_DEG` — a table beside `COMPARTMENT_BEARINGS`. **Every entry is `0.0`**, which
  is today's behaviour bit-identically. What each owner's value *should* be is a physiological
  question with a different answer per owner, and this entry does not answer six of them at once.
- `_spread_in_band(owner, spread_deg, bearing)` — called after `_place_in_band`, it rotates each
  **declared sub-assembly** about the **cell centre** onto its own bearing, sampled on a Fibonacci
  cap of half-angle `spread_deg` about the compartment's bearing.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** — no provider commit is behind this |
| Source path and symbol | **not applicable** — no `ffn_sim` path and no provider symbol was read |
| Read from | this tree's `_place_in_band`, and `-3661`'s measurement of what a count does to a blob |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. **RE-DERIVED** from this repository's own geometry.

## 4. Physical or mathematical law represented

**A rotation about the cell centre is the one map that spreads an owner and costs nothing**, and
that is the whole design:

| | preserves internal distances | preserves each node's radius | spreads over the sky |
|---|---|---|---|
| translation (`_place_in_band` today) | **yes** | no | **no** |
| scaling (what `_place_in_band` replaced) | no — **4,323 pN·µm injected** | no | no |
| **rotation about the centre** | **yes** | **yes** | **yes** |

Preserving internal distances is what makes it inject **exactly zero** energy — the same argument
`_place_in_band`'s docstring already makes for translation. Preserving each node's radius is
stronger than translation manages and it is what makes the change safe: **the band test
`inner ≤ |r| ≤ outer` is invariant under it, so no owner that fitted before can stop fitting.**

## 5. Units, domains, singular cases, invariants

- **I1 — `spread_deg = 0.0` is bit-identical to no call at all.** `np.array_equal` on every owner's
  node array. This is what makes the table of zeros a true no-op rather than a near-one.
- **I2 — every node's radius is preserved to round-off.** Asserted directly, because it is what
  guarantees I3.
- **I3 — no owner that fitted its band before can fail it after.** A consequence of I2, asserted
  separately anyway, because the band check is the thing a caller notices.
- **I4 — zero energy is injected.** The owner's `potential_energy_pn_um()` is unchanged to
  round-off on a quadratic form evaluated 4 µm from the origin.
- **I5 — an owner with a SATELLITE POINT is REFUSED, by name.** `MicrotubuleNetwork` carries
  `hub_position_um`, a centrosome every rod's minus end is anchored to by a 500 pN/µm spring of rest
  0.1 µm. Rotating rods about the cell centre while one shared hub stays put would stretch every one
  of those anchors — **the `_place_in_band` docstring records that this owner's satellite point has
  already been forgotten once**, and a silent repeat here would inject energy through the one
  attribute the node loop cannot see. Refused rather than special-cased, because the physiological
  answer for microtubules is *"the network is not an aster"* (`-3661` §14, 0.7–3.9 % centriole-
  connected) and that is a topology change, not a placement one.
- **I6 — sub-assemblies are read from the owner's own declarations.** `Rod`, `Bundle` and the
  intermediate filament's cable all publish `node_start` / `node_count`, and each owner asserts its
  ranges tile `[0, N)`. An owner that declares none is left whole, which is a spread of one.

## 6. Source evidence class and known retractions

Not a literature claim. Evidence class: **geometry and measurement on this tree.** The *values* in
`COMPARTMENT_SPREAD_DEG` are all zero precisely so that this entry makes no physiological claim; the
one owner with a sourced answer is `intermediate_filament` (**V-6**, an interpenetrating cell-wide
network) and even there the extent is a named blank.

## 7. Independent oracle or derivation

1. **Rigidity, checkable without running the cell.** A rotation matrix `R` satisfies
   `‖R(x−c) − R(y−c)‖ = ‖x−y‖` for any `x`, `y` in the same sub-assembly, and with `c = 0` also
   `‖Rx‖ = ‖x‖`. So I2 and I4 are theorems and the controls confirm the implementation, not the
   claim.
2. **The baseline this has to beat is measured, not assumed.** `-3661`, whole cell, the coverage
   each owner spans about the cell centre: `filopodium` **0.01 %**, `lamellipodium` 0.64 %,
   `intermediate_filament` **2.16 %**, `sf_arc` 4.33 %, `microtubule` 6.04 %, `ecm` 7.18 %, against
   `membrane`/`nucleus`/`cortex` at **100 / 100 / 99.6 %**.
3. **Prediction, stated before the code and deliberately arithmetic.** A cap of half-angle `θ` covers
   `(1 − cos θ)/2` of the sky. An owner spread over a cap of `θ = 60°` should therefore reach about
   **25 %**, and one at `θ = 180°` about **100 %** — up to the width of the sub-assemblies
   themselves, which adds a little. So the control asserts `intermediate_filament` at
   `spread_deg = 60` lands in **20–35 %** against its 2.16 % today. **A measured coverage that does
   not move refutes the mechanism**; one far above 35 % means the rotation is not confined to the
   cap it was given.

### §7(3) measured, and the prediction was low. Why is the interesting part.

| `spread_deg` | ideal cap `(1−cos θ)/2` | **measured coverage** |
|---:|---:|---:|
| 0 (today) | — | **2.16 %** |
| 60 | 25.0 % | **16.09 %** |
| 90 | 50.0 % | **32.28 %** |
| 120 | 75.0 % | **49.13 %** |
| 180 | 100.0 % | **93.54 %** |

**Predicted 20–40 % at 60°, measured 16.1 %.** Two things the prediction did not account for, and
both are consequences of this owner having only **three** sub-assemblies:

1. the widest Fibonacci sample on a cap of half-angle `t` sits at `cos = cos t + (1−cos t)/(2N)`,
   which at `N = 3`, `t = 60°` is **54.3°**, not 60°;
2. the coverage is measured about the **empirical mean** of the node cloud, and the mean of three
   unbalanced samples is not the axis — taking the measured extent to **47.3°**.

**So the shortfall is a property of how many pieces there are to spread, not of the rotation**, and
the ladder confirms it: the same three pieces reach **93.5 %** at 180°.

**This is the honest statement of how D1 and D3 relate, and it is better than the plan's.**
`§D` of the defect plan says *"D1 gates the others."* `-3661` measured that it does not — it absorbs
them. This entry measures the other direction: **D3 supplies the pieces and D1 spreads them, and
neither is sufficient alone.** Three cables spread over 180° cover 93.5 % of the sky and are still
three cables; 322 microtubules in a bounded fan cover 5.61 % and are still one blob.

### The whole-cell consequence, measured — and it splits the owners in two

§13 made this a precondition of `ACCEPTED` and it is now taken. Whole cell, one build per row:

| | total energy | **owner energy sum** | max residual | audit tally |
|---|---:|---:|---:|---|
| today (all zeros) | 12,421.05 | **9,243.088420906988** | 1,282.172 | 18 / 15 / 1 |
| `intermediate_filament` @ 90° | 12,448.58 (**+0.22 %**) | **identical** | **identical** | unchanged |
| `intermediate_filament` @ 180° | 12,457.45 (**+0.29 %**) | **identical** | **identical** | unchanged |
| `sf_arc` @ 45° | 15,430.91 (**+24.2 %**) | **identical** | 2,287.879 (**+78.5 %**) | unchanged |

**Three readings, and the first is the one this entry was built to earn.**

1. **The owner energy sum is bit-identical in every row** — `9243.088420906988`, to the last digit,
   whether nothing is spread or an owner is spread over the whole sphere. **I4 confirmed at
   whole-cell scale**, and it is the difference between this and the scaling placement that injected
   4,323 pN·µm.
2. **Spreading the uncoupled owner is nearly free.** `intermediate_filament` over the *entire
   sphere* costs **+0.29 %** of the cell's connector energy and moves the maximum residual by
   **nothing at all** — it is not on the residual's critical path either way. Against a source that
   calls this network cell-wide and interpenetrating, and a construction that gives it **2.16 %** of
   the sky, that is a repair available for essentially no price.
3. **Spreading the coupled owner is not.** `sf_arc` at a mere **45°** costs **+24 % of total energy
   and +78.5 % of the residual**, because its focal-adhesion and dorsal-arc rest lengths were taken
   from where it sat. **The cost of D1 is entirely in the connectors and never in the owners**, and
   it scales with how many load paths an owner is on — not with how far it is moved.

**No audit row changes verdict in any configuration**, so nothing is silently unwired by this.

**The default is still not moved, and the reason is not the cost.** For
`intermediate_filament` the cost is 0.29 %. What is missing is the *number*: V-6's cell-wide extent
is a named blank, "180°" is not sourced, and `-3652` §13's rule — the entry that builds a mechanism
does not move a default — is exactly the rule this lane has now watched pay for itself three times.
**The decision now rests on a sourced spread angle, not on affordability.**

**And one thing this entry deliberately does NOT predict:** whether spreading changes the cell's
energy or its convergence. It cannot change the *owner's* energy (I4), but it moves owners relative
to each other, so every **connector** rest length and every nearest-point correspondence is
re-formed somewhere else. That is a whole-cell measurement and it belongs to whoever turns a value
in the table non-zero.

## 8. Positive control

`tests/scenarios/test_compartment_spread.py`:

- `test_zero_spread_is_bit_identical` — I1, `np.array_equal` on every owner in the built cell.
- `test_every_radius_is_preserved` — I2.
- `test_the_owner_energy_does_not_move` — I4.
- `test_the_band_still_holds` — I3, driven through the real placement so the refusal path is live.
- `test_spreading_the_intermediate_filament_covers_the_cap` — §7(3): 2.16 % → 20–35 % at 60°.

## 9. Deliberately failing negative control

- `test_an_owner_with_a_satellite_point_is_refused` — I5, `microtubule`, and the message must name
  the attribute so the reader learns *why* rather than *that*.
- `test_spread_is_not_vacuous` — the vacuity control: a non-zero spread must actually move node
  positions. Without it every invariant above would pass against a function that returned early.
- `test_a_negative_or_oversized_spread_is_refused` — the domain is `[0, 180]`.

## 10. Numerical and precision envelope

float64. I1 is `np.array_equal` — exact, because the claim is bit-identity of a no-op. I2 and I4 are
`np.allclose` at `rtol=1e-12`: a rotation is exact in exact arithmetic and round-off on coordinates
4 µm from the origin is the only difference there can be, so a *tolerance* is right there and an
equality would be wrong.

## 11. Production-backend residency and transfer

**None.** Placement runs at build, on the host, before any backend array is formed. Node positions
change; shapes, dtypes and orderings do not.

## 12. Comments and docstrings to discard

Nothing copied. `_place_in_band`'s docstring gains a pointer to this entry, and keeps its own
account of the 4,323 pN·µm the scaling placement injected — that account is *why* this entry uses a
rotation rather than anything cheaper.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass. **`ACCEPTED` requires a non-zero value in
`COMPARTMENT_SPREAD_DEG` for at least one owner, with the whole-cell consequence measured** —
connector energies, the load-path audit, and the residual — because a mechanism nobody turns on
repairs nothing. **This entry deliberately stops one step short of that**, for the reason `-3652`
§13 gives: the entry that builds a mechanism is not the entry that moves a default.

## 14. Honest limits

- **This does not repair D1. It makes D1 repairable.** Every value is zero, so the cell built today
  is the cell built yesterday, to the bit. What changes is that spreading an owner is now a number
  rather than a rewrite.
- **`microtubule` is refused, and it is the owner D1 was most visible on.** Its single-hub topology
  is what blocks it, the source says that topology is wrong (0.7–3.9 % centriole-connected), and
  neither is fixed here.
- **Three owners should probably never be spread**: both protrusions are localised structures by
  definition, and `ecm` is a substrate placed by contact. A table of six numbers invites someone to
  fill all six; the docstring says which three are already right.
- **`sf_arc` is ventral, not spherical.** Its correct spread is a cap of some tens of degrees, not
  180°, and nothing here sources that number.
- **The sub-assembly spread is deterministic and uncorrelated with anything.** A real cytoskeleton's
  orientation correlates with the cell's polarity, the substrate and the nucleus. This is a
  Fibonacci cap, and `ALEPH-PORT-3638` already records what a Fibonacci lattice costs in realism:
  it is a crystal, 45× more ordered than a random arrangement on the nearest-neighbour statistic.
