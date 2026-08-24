# ALEPH-PORT-3661 — the microtubule count is a length density, and Aleph declares two rods

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3661` |
| Lane | `b4fad06b` |
| Status | `AUDITED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Answers | **V-5**, and defects **D2** and **D3** for the `microtubule` owner |
| Follows | `filament_count_for_density`, which did exactly this for the cortex |

---

## 1. Aleph API

`aleph/vertical/microtubule_sourcing.py` — a new module, beside `cortex_sourcing.py` and shaped
like it:

- `MICROTUBULE_LENGTH_DENSITY_UM_PER_UM3` — a `LiteratureRange`, the sourced band.
- `microtubule_count_for_length_density(cytoplasm_volume_um3, rod_length_um, *, length_density_um_per_um3) -> int`.

**No default moves.** `_microtubule` in `aleph/scenarios/whole_cell.py` keeps `_fx.repeats()`; what
this entry adds is the function that says what the count **should** be, and a measurement of what
happens when a caller asks for it. `-3652` §13's rule — a default may not move inside the entry that
sources its value — applies here too.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** — no provider commit is behind this |
| Source path and symbol | **not applicable** — no `ffn_sim` path and no provider symbol was read |
| Read from | Müller A. et al. (2020) *3D FIB-SEM reconstruction of microtubule–organelle interaction in whole primary mouse β cells*, **J Cell Biol 219(12):e202010039**, retrieved via PubMed `PMC7748794` and read in full |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. **RE-DERIVED** from a paper.

## 4. Physical or mathematical law represented

**The sourced object is a length density, not a count, and the source is what says so.**

Müller et al. reconstruct **every** microtubule in seven whole primary cells at 4 nm near-isotropic
FIB-SEM, with ≈50 h of manual tracing per cell. Across a glucose stimulation that changes the
**count** by 3× and the **mean length** by 4×, the cumulative length per cell barely moves:

| | low glucose (3 cells) | high glucose (4 cells) |
|---|---|---|
| microtubules per cell | 291 / 342 / 440 | 862 / 932 / 1,007 / 1,101 |
| mean length | 8.54 / 8.82 / 14.63 µm | 2.66 / 3.10 / 3.14 / 3.25 µm |
| **cumulative length** | **3,018 / 3,756 / 4,260 µm** | **2,479 / 2,799 / 3,115 / 3,452 µm** |
| cell volume | 799 / 927 / 997 µm³ | 793 / 834 / 898 / 1,010 µm³ |

Verbatim: *"the cumulative microtubule length per cell showed only modest variation between glucose
conditions"* and *"we did not find major differences in polymerized tubulin density between glucose
conditions."*

```text
rho_L  =  cumulative length / cell volume  =  2.75 … 4.73 um per um^3
count  =  rho_L * cytoplasm_volume / rod_length
```

**This is the same shape as the cortex's areal density and for the same reason.** A count is a
property of the construction; a density is a property of the material. `filament_count_for_density`
already encodes that for the cortex and this is its volumetric twin.

## 5. Units, domains, singular cases, invariants

- **I1** — `count >= 1`; a volume or a rod length that is not finite and positive is refused.
- **I2** — the count scales **linearly** in the volume and **inversely** in the rod length. Both are
  checkable without any measurement and both fail loudly if the formula is inverted.
- **I3** — **the returned count reproduces the source's own cells.** Feeding the paper's cell volumes
  and mean rod lengths back in must return counts inside the paper's own measured range. A sourcing
  function that cannot reproduce its own source is not sourced.
- **I4** — nothing in this module reads the construction. No mesh, no band, no spacing.

## 6. Source evidence class and known retractions

Peer-reviewed, open access, 102 citations, whole-cell volume EM with manual tracing — the strongest
evidence class in this file's V-table. No retraction known.

**The limitation is the cell type, and it is stated in the value's own docstring.** These are mouse
pancreatic β cells of 793–1,010 µm³. Aleph's cell is a 5 µm-radius sphere, 523 µm³, of which
≈410 µm³ is cytoplasm. **The count is not transferable and the density is what carries** — which is
the whole reason the sourced object is a density.

## 7. Independent oracle or derivation

**Three, and the third is the one that matters.**

1. **Round-trip on the source's own cells** (I3): `rho_L × 900 µm³ / 8.8 µm ≈ 281 … 484` against a
   measured 291–440 at low glucose, and `rho_L × 900 / 3.0 ≈ 825 … 1,419` against a measured
   862–1,101 at high. **Both bands overlap the measurement**, which is the least a density owes its
   source.
2. **Applied to Aleph's own geometry**, with the rod length the card already fixes —
   `8 nodes × card.rest_spacing_um = 3.5 µm`:

   ```text
   2.75 … 4.73 um/um^3  x  410 um^3  /  3.5 um   =   322 … 554 rods
   ```

   **Aleph declares TWO.** That is **160–277× short**, and the sourced figure lands in the same
   place as the source's own direct count of 291–440 — from a density and a rod length neither of
   which was chosen to make that happen.

3. **The prediction that this entry expects to be REFUTED by the construction, and says so first.**
   A rod is 3.5 µm; the cytoplasmic band `microtubule` is placed into is `(3.1, 4.9)`, i.e.
   **1.8 µm thick**. `_place_in_band` is a rigid translation that refuses an owner reaching further
   from its centroid than the band's half-width. Two rods fit because they are laid nearly
   tangentially and splayed 37°. **322 rods cannot be**: the fan construction confines every rod to
   `direction = (1.0, 0.35·s_y, 0.15·s_z)`, so raising the count makes a denser fan in the same
   place, not a cell-filling network.

   **So the prediction is that asking for the sourced count REFUSES, and the refusal is the
   measurement of D1.** `docs/superpowers/plans/.../every-defect-found...` §D says *"D1 gates the
   others — raising the element count of a blob gives a denser blob."* This entry is the arithmetic
   that turns that sentence into a number. **If the sourced count builds without complaint, D1 is
   less binding than the plan says and that is a finding too.**

   ### §7(3) IS REFUTED, AND WHAT REPLACES IT IS WORSE

   **It builds.** 322 rods, 2,576 nodes, **19.1 s**, and `potential_energy_pn_um() == 0.0` exactly —
   no prestrain at all. `_place_in_band` never gets near its refusal, because the rods are laid
   nearly tangentially and the fan is a *bounded* envelope.

   | `cytoskeleton_density` | rods | nodes | **sky the owner spans** | θ_max | energy |
   |---:|---:|---:|---:|---:|---:|
   | 1 | 2 | 16 | **6.04 %** | 28.5° | 0.0000 |
   | 8 | 16 | 128 | 5.49 % | 27.1° | 0.0000 |
   | 32 | 64 | 512 | 5.60 % | 27.4° | 0.0000 |
   | **161 (sourced)** | **322** | **2,576** | **5.61 %** | 27.4° | 0.0000 |

   **161× the microtubules and the angular extent does not move.** It falls slightly. Every rod is
   confined to `direction = (1.0, 0.35·s_y, 0.15·s_z)`, so raising the count packs the same wedge
   tighter instead of filling the cell.

   **This is a worse outcome than the predicted refusal, because a refusal tells you.** The plan's
   sentence *"raising the element count of a blob gives a denser blob"* is confirmed, with a number:
   the count is free to move by **two orders of magnitude** while nothing about the geometry does,
   and nothing in the tree complains. **D1 does not *gate* D3 — it silently absorbs it**, which is
   the harder failure to notice and the reason this measurement is the entry's main result rather
   than its footnote.

## 8. Positive control

`tests/vertical/test_microtubule_sourcing.py`:

- `test_the_density_reproduces_the_sources_own_counts` — I3, both glucose conditions.
- `test_the_count_scales_with_volume_and_inversely_with_length` — I2.
- `test_alephs_geometry_asks_for_hundreds_of_rods` — the 322–554 band against the declared 2, pinned
  so the 160–277× gap cannot quietly close by someone editing a literal.
- `test_the_band_carries_its_citation` — the `LiteratureRange` refuses to exist without one, and
  this asserts the citation names the paper rather than being a placeholder.

## 9. Deliberately failing negative control

- `test_a_non_positive_volume_or_rod_length_is_refused` — I1.
- `test_the_count_is_not_read_off_the_construction` — a **vacuity control**: the function is called
  with a volume and a rod length that have nothing to do with this repository's cell and must still
  return the density's answer. A function that had a mesh or a band in it would not.
- `test_the_sourced_count_builds_and_is_still_a_blob` — §7(3) **as corrected by its own
  measurement.** It asserts the build succeeds, at zero prestrain, **and that the angular extent
  does not grow.** If the coverage assertion ever fails, D1 has been repaired and this control is
  what notices.

## 10. Numerical and precision envelope

float64; the count is `int(round(...))`. Band comparisons use the `LiteratureRange` endpoints
directly rather than a midpoint, so nothing in this entry averages two measurements into a third
number nobody made.

## 11. Production-backend residency and transfer

**None.** This module returns an integer. Nothing is allocated, no array crosses a device.

## 12. Comments and docstrings to discard

Nothing copied. The paper is quoted verbatim in §4 and cited.

## 13. Acceptance

`AUDITED`. §8 and §9 pass — `tests/vertical/test_microtubule_sourcing.py`, **16 passed**.

**`ACCEPTED` is withheld, and the reason changed while the entry was being written.** It was going
to be *"the construction cannot yet accept the sourced count"*. It can — it accepts it silently and
puts 322 microtubules in 5.61 % of the cell. So the blocker is not acceptance, it is that
**accepting it means nothing until the placement distributes**, which is defect D1 and is not in
this entry. A sourced count on a blob is a correct number in the wrong place.

## 14. Honest limits

- **One cell type, one tissue, one species.** Mouse pancreatic β cells. Aleph has no cell type at
  all, so this is the density it borrows, and a different tissue would give a different one.
- **The rod length is the card's, not the source's.** `8 × rest_spacing_um = 3.5 µm` against a
  measured 2.66–14.63 µm. It is inside the range but it was not chosen from it, and a rod length
  that moves moves the count inversely.
- **The topology is still wrong and this entry does not touch it.** Every rod in this tree anchors to
  one `hub_position_um`; the source measures **0.7–3.9 %** of microtubules centriole-connected and
  states *"over 80 % of all microtubules are disconnected from these compartments"*. **A count fixed
  on an aster is a correct number of rods in a network shape the source refutes**, and that is a
  larger change than a count.
- **No length distribution.** The cortex got one in `-3656`; this is a single rod length. The source
  plots a length distribution per cell (its Fig 2B) and it is not drawn here.
