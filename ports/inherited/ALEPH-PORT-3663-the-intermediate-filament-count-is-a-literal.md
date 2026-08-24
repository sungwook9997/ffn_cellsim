# ALEPH-PORT-3663 — the intermediate filament count is a literal, and it is what caps D1

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3663` |
| Lane | `b4fad06b` |
| Status | `AUDITED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Answers | defect **D2** for the `intermediate_filament` owner; **V-6** is addressed and stays a named blank, for a reason §6 gives |
| Follows | `ALEPH-PORT-3662`, which measured that three sub-assemblies reach only 16.1 % of a 25 % cap |

---

## 1. Aleph API

`aleph/scenarios/whole_cell.py::_intermediate_filament` gains the `_fx.repeats()` pattern that
`_sf_arc` and `_microtubule` already have: `cytoskeleton_density` copies of the three-cable seed,
fanned across a bounded envelope. **At `cytoskeleton_density = 1` — the default — it is the same
three cables, coordinate for coordinate.**

**No new sourced constant.** V-6 is not filled; see §6.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** — no provider commit is behind this |
| Source path and symbol | **not applicable** — no `ffn_sim` path and no provider symbol was read |
| Read from | this tree's `_sf_arc` and `_microtubule`, which already do this, and `-3662`'s coverage ladder |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. **RE-DERIVED** — the same pattern two sibling owners have.

## 4. Physical or mathematical law represented

**None, and this entry is careful to claim none.** What it changes is that an element count stops
being a literal. `_sf_arc` scales with `cytoskeleton_density` and so does `_microtubule`;
`_intermediate_filament` does not, so **the one owner with a source calling its network cell-wide is
the one owner whose count cannot move at all.**

**Why it matters beyond tidiness, and this is measured rather than argued.** `ALEPH-PORT-3662`
spreads an owner by rotating its **sub-assemblies** onto separate bearings, and its coverage is
therefore bounded by how many there are. Measured there: three cables on a 60° cap reach **16.1 %**
against the cap's own **25 %**, because the widest of three Fibonacci samples sits at 54.3° and the
coverage is read about the empirical mean of three unbalanced points. **The shortfall is the count.**

## 5. Units, domains, singular cases, invariants

- **I1 — `cytoskeleton_density = 1` is bit-identical.** `np.array_equal` on the built node array
  against the three literal cables. Every committed whole-cell number was taken there.
- **I2 — the cable count is `3 × cytoskeleton_density`**, and each cable keeps its population
  (`KERATIN`, `VIMENTIN`, `KERATIN`) so the two-population structure is replicated rather than
  diluted.
- **I3 — the fan is BOUNDED.** `_sf_arc` records what an unbounded offset costs: `1.15·copy − 0.6`
  put copy 7 **7.45 µm** from copy 0 inside a 10 µm cell, and only the old *scaling* placement hid
  it. The band here is `(3.1, 4.2)` — **0.55 µm of half-width**, the tightest in the table — so the
  envelope has to fit inside that at any density.
- **I4 — the owner still builds at zero prestrain.** The cables are placed at cell scale about their
  own centroid and `_place_in_band` translates rigidly; replication must not change that.

## 6. Source evidence class and known retractions

**V-6 is addressed and deliberately NOT filled, and the reason is a units mismatch rather than a
missing paper.**

`docs/design/SOURCED_CONSTRUCTION_VALUES.md` §V-6 carries **≈ 11.7 vimentin filaments per µm³**,
from Renganathan et al. 2025 — 583 filaments reconstructed by 3D FIB-SEM in a ≈50 µm³ window. Two
things stop it becoming a count here:

1. **The window was chosen in the densest region.** It is perinuclear, and the paper says so. The
   density is an **upper bound**, not a cell average.
2. **A `Cable` in this tree is not a filament in that paper's sense.** Aleph's cable is a
   coarse-grained bundle of ~5 nodes spanning ~0.9 µm; the source counts **individual** vimentin
   filaments, and separately measures that they *"form only loosely organized, semi-coherent
   structures"* — i.e. the bundle is exactly the object whose internal structure the source resolves
   and this model does not. **There is no sourced conversion between the two**, and inventing one is
   the defect `PROPOSAL-a-bond-has-no-formation-criterion…` §7 names.

**So the count becomes a parameter and its value stays a named blank.** That is the whole of what
this entry claims, and it is less than `-3661` claimed for the microtubule, on purpose: the
microtubule source counts the same object Aleph builds and this one does not.

## 7. Independent oracle or derivation

1. **I1 is checkable without any physics** — build at density 1, compare arrays.
2. **The prediction, and it is `-3662`'s ladder run again with more pieces.** More sub-assemblies
   should close the gap between the measured coverage and the cap's own `(1 − cos θ)/2`, because
   both causes of the shortfall shrink with `N`: the widest sample sits at
   `cos = cos θ + (1 − cos θ)/(2N)`, and the empirical mean of `N` samples approaches the axis.

   At `θ = 60°` the cap is **25.0 %** and three cables reached **16.1 %**. With
   `cytoskeleton_density = 8` — **24 cables** — the prediction is **22–25 %**, i.e. most of the gap
   closed. **A coverage that does not improve refutes §4's explanation of the shortfall**, which
   would mean the 16.1 % has a different cause and `-3662`'s account of it is wrong.
3. **The band must still hold at density 8**, and if it does not, I3's envelope is too wide and the
   refusal names it — which is `_place_in_band` working, not a failure of this entry.

### §7(2) measured — and it is the one prediction in this lane's run that HELD

| `cytoskeleton_density` | cables | coverage, **no spread** | coverage, **60° spread** | ideal cap |
|---:|---:|---:|---:|---:|
| 1 | 3 | 2.16 % | **16.09 %** | 25.00 % |
| 2 | 6 | 3.13 % | 20.46 % | 25.00 % |
| 4 | 12 | 3.13 % | 23.04 % | 25.00 % |
| **8** | **24** | 3.13 % | **24.41 %** | 25.00 % |
| 16 | 48 | 3.13 % | 25.60 % | 25.00 % |

**Predicted 22–25 % at 24 cables. Measured 24.41 %.** The gap closes monotonically in the count and
overshoots slightly at 48, by the sub-assemblies' own angular width — which is the term §7(2) said
would add a little. Energy is `≈ 2 × 10⁻²⁷ pN·µm` at every density, i.e. zero: I4 holds.

**Three numbers that close the D1/D3 question between them, all on this one owner:**

```text
D3 alone   more pieces, no spread     2.16 %  ->   3.13 %, then FLAT      -- a bounded fan
D1 alone   spread, three pieces       2.16 %  ->  16.09 %                 -- capped by the count
both       24 cables at 60 deg                ->  24.41 % of a 25.00 % cap
```

**Neither is sufficient and together they reach the cap.** That is the honest relation between the
two defects, and it is not `§D`'s *"D1 gates the others"*: `-3661` measured that D1 does not gate
D3, and this measures that D3 does not substitute for D1 either.

**And the first envelope written for this entry was REFUSED, correctly.** `_place_in_band` rejected
it at `r = [3.1432, 4.2497]` against `(3.100, 4.200)` — over by 0.0497, which is exactly the
`0.10 × fraction` this entry had put in **x**, and the bearing is `(−1, 0, 0)` so x is the radial
direction. The offset is tangential only for that reason. **The guard caught it, and it was not
argued with.**

## 8. Positive control

`tests/scenarios/test_intermediate_filament_count.py`:

- `test_density_one_is_bit_identical` — I1.
- `test_the_cable_count_scales_with_the_density` — I2, and that populations repeat rather than drift.
- `test_the_owner_still_builds_unstrained` — I4, at densities 1 and 8.
- `test_more_cables_close_the_gap_to_the_cap` — §7(2), the claim this entry exists for.

## 9. Deliberately failing negative control

- `test_a_denser_fan_still_fits_its_band` — I3 at density 8, through the real placement.
- `test_the_replication_is_not_vacuous` — a vacuity control: density 8 must produce different node
  positions and more cables, not the same three with a bigger array.

## 10. Numerical and precision envelope

float64. I1 is `np.array_equal`. The coverage comparison in §7(2) is a band, because it is a
prediction about a Fibonacci sample's extreme value and not a closed form.

## 11. Production-backend residency and transfer

**None.** Fixture construction on the host; array sizes grow with the count, which every existing
whole-cell control already covers at other densities.

## 12. Comments and docstrings to discard

Nothing copied. The fixture's docstring gains the reason the count was a literal and what the
sourced density cannot say.

## 13. Acceptance

`AUDITED`. §8 and §9 pass — `tests/scenarios/test_intermediate_filament_count.py`, **11 passed**,
alongside `-3662`'s 14. **`ACCEPTED` requires a sourced count**, and §6 says why
there is not one: the source resolves individual filaments and this model builds bundles. **A
conversion between them is the missing measurement**, not a missing paper.

## 14. Honest limits

- **The count is a parameter with no sourced value.** That is strictly better than a literal with no
  sourced value, and strictly worse than `-3661`'s microtubule, where the source counts the object
  the model builds.
- **A denser fan is still a fan.** Without `COMPARTMENT_SPREAD_DEG` non-zero for this owner — which
  is `-3662`'s decision and is still zero — 24 cables occupy the same 2.16 % of the sky that three
  did. **The two entries are only useful together**, and that is the point both of them make.
- **The two-population split is replicated, not modelled.** Keratin and vimentin are different
  proteins in different cells; three cables repeated eight times is eight copies of one arbitrary
  mixture.
- **No length distribution and no interpenetration.** The source's headline structural finding is
  that vimentin and F-actin form an *interpenetrating* network. Nothing here represents that.
