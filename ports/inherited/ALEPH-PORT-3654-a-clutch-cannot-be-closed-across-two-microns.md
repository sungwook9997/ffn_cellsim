# ALEPH-PORT-3654 — a clutch cannot be closed across two micrometres

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3654` |
| Lane | `a40e55a2 S-OBSERVE handover, scope extended by PI reassignment` |
| Status | `PROPOSED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Answers | `PROPOSAL-a-bond-has-no-formation-criterion-and-no-rupture-criterion.md` clause **B-1**, and the PI's first observation of 2026-08-07: *"adhesion 거리가 전에 크로스링크 문제처럼 말이 안되게 지정이 되어있는데"* |

---

## 1. Aleph API

`aleph/vertical/adhesion_sourcing.py`, a new module:

- `FA_CORE_HEIGHT_UM`, `TALIN_LENGTH_UM`, `FA_PLAQUE_UM` — `LiteratureRange`, each with its paper.
- `ClutchOutOfReachError(ValueError)`.
- `clutch_reach_um() -> float` — the furthest a closed integrin–talin–actin linkage can span.
- `assert_clutch_can_be_closed(name, separation_um, rest_gap_um) -> None` — refuses, by name and
  with the number, a clutch declared engaged across a distance no linkage reaches.

Called from a probe over `slot.owners` rather than from inside the connector, so `aleph/observe/`
and the scenarios can both ask the question without a step running.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none — no provider file was read.** |
| Source commit | **not applicable**; nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source path | **not applicable** |
| Source symbol(s) | **not applicable**; the Aleph symbols are named in §1 |
| Read from | **primary literature**, recorded in `docs/design/SOURCED_CONSTRUCTION_VALUES.md` §V-1/V-2 |
| Working tree == commit? | **not applicable** — no provider revision cited |

- Kanchanawong P. et al. (2010) *Nature* — *"integrins and actin are vertically separated by a
  ∼40-nm focal adhesion core region"*; the plaque is *"a <200-nm plaque linking the ECM to the
  actin cytoskeleton"*; the companion gives *"a 30-50 nm FA core which is spanned by talin tethers"*.
- Liu J. et al. (2015) *PNAS* — *"Talin was found to be ∼97 nm in length and oriented at ∼15°
  relative to the plasma membrane"*, and length-modified talins moved the FA/stress-fibre interface
  *"in a linear manner"*, which is what makes it a **ruler** rather than a correlate.
- Dedden D. et al. (2019) *Cell* — full-length talin1 reversibly unfolds to *"an ∼60-nm string-like
  conformation"*.

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. Stated rather than left blank.

## 4. Physical or mathematical law represented

Not a force law — the **domain on which the existing force law is defined**.

A clutch is a series linkage: integrin → talin → actin. At rest it spans the focal-adhesion core,
**30–50 nm**, and `ClutchCard.rest_gap_um = 0.05` sits at the top of that — **the rest gap is right
and is not what this entry changes**. What has no bound at all is how far the two ends may be apart
while the clutch is still declared **closed**.

The bound is the linkage's own contour. Talin is the ruler: **97 nm** end to end. A closed
integrin–talin–actin link cannot span further than the molecule that constitutes it, whatever the
builder placed where.

`SeriesJointConnector` computes `extension = distance − rest_gap_um()` and gates delivery on
`is_engaged()`, which is a `BindingState` set at construction and which **nothing in the tree can
change** — `focal_adhesion.py` carries **0** `rate_per_s` against `nmii.py`'s **20**. So an
unbounded distance times a stiffness is delivered as force.

## 5. Units, domains, singular cases, invariants

- µm. `0 < rest_gap ≤ reach`.
- **I1.** `reach = TALIN_LENGTH_UM.hi = 0.097` µm. It does not depend on the mesh, the density, or
  any separation the construction produced.
- **I2.** A clutch with `separation > reach` is **not closable**, independently of its stiffness.
- Singular case: `separation ≤ rest_gap` — a compressed clutch. **Not this entry's business**; the
  spring's own sign handles it and refusing there would forbid a legitimate configuration.

## 6. Source evidence class and known retractions

Peer-reviewed 3D super-resolution (iPALM) and cryo-EM. Kanchanawong 2010 has 1,536 citations and
Liu 2015 has 190. No retraction known for either. Dedden 2019 is used only to record that talin's
*relaxed* form is shorter than its contour, i.e. that 97 nm is a ceiling and not a rest length.

## 7. Independent oracle or derivation

**The arithmetic is the oracle and it is stated before the code runs.** Measured spans at build,
`subdivision_level=1`:

| row | span | reach `0.097 µm` | prediction |
|---|---:|---:|---|
| `filopodium_nascent_fa` | 2.8398 µm | | **refuse, 29.3×** |
| `lamellipodium_nascent_fa` | 1.9892 µm | | **refuse, 20.5×** |
| `alpha2beta1_collagen_series` #0 | 2.2004 µm | | **refuse, 22.7×** |
| `alpha2beta1_collagen_series` #1 | 2.1927 µm | | **refuse, 22.6×** |

**All four refuse. That is the predicted result and it is not a regression** — it is step 2 of the
order pre-registered in `docs/superpowers/plans/2026-08-08-the-construction-layer.md` §2, and the
count is the map for steps 3 and 4.

A second, independent check with no molecule in it: the spans are **11–14× the entire <200 nm
plaque**. A bond longer than the organelle that contains it fails on containment alone.

## 8. Positive control

`tests/vertical/test_adhesion_reach.py`:

- `test_a_clutch_at_its_rest_gap_is_closable` — 0.05 µm passes.
- `test_a_clutch_at_the_talin_contour_is_closable` — 0.097 µm passes; the boundary is inclusive
  because a molecule at full extension is still a molecule.
- `test_the_reach_is_the_talin_ruler_and_not_a_literal` — `clutch_reach_um()` is read from
  `TALIN_LENGTH_UM`, so number and citation cannot drift apart.
- `test_the_rest_gap_sits_inside_the_sourced_core` — 0.05 µm is inside 30–50 nm, i.e. the shipped
  card is **defended, not changed**.

## 9. Deliberately failing negative control

- `test_the_four_shipped_clutches_are_all_out_of_reach` — drives the four measured spans and
  asserts each refuses, with the ratio in the message.
- `test_a_clutch_beyond_the_plaque_is_refused_even_without_the_molecule` — 0.25 µm exceeds the
  <200 nm plaque and refuses, so the finding does not rest on talin alone.
- `test_the_check_can_pass_and_therefore_can_fail` — vacuity control.
- `test_a_compressed_clutch_is_NOT_refused` — a separation below the rest gap must pass; refusing
  there would be a guard that fires on a legitimate configuration.

## 10. Numerical and precision envelope

Exact float64 comparisons against declared constants; the boundary is inclusive at `reach`. No
tolerance — the measured spans are 20–29× outside, nowhere near it.

## 11. Production-backend residency and transfer

**None.** Host-side scalars over an existing owner block; nothing is launched and nothing is
transferred. It reads `slot.owners[*].positions`, which the out-of-step probes already read.

## 12. Comments and docstrings to discard

Nothing copied. Two docstrings are **quoted for refutation** and cited: the registry's *"a stochastic
clutch state, not a spring anchored to a place"* against the implementation's
`extension = distance − rest_gap`, and `focal_adhesion.py`'s *"`is_engaged` … is both sides bound"*.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass and the four refusals are published with their ratios.
**It may not reach `ACCEPTED` while the refusal is not wired into a build** — this entry supplies
the criterion and the count; making `build_whole_cell` act on it moves every adhesion number and
belongs to its own entry, for the same reason `-3652` did not move a default.

## 14. Honest limits

- **This is B-1 only.** B-3 and B-4 — an on-rate and a force-dependent off-rate — are not here. A
  reach makes the bond's *existence* honest; it does not make it *dynamic*, and `focal_adhesion.py`
  still carries no kinetics at all.
- **It refuses; it does not place.** Nothing here moves an integrin closer to its ligand. The
  refusals are a map of where the construction is wrong, and the construction is `_place_in_band`.
- **Talin is the ruler for these three rows.** A clutch built on a different linkage would need its
  own contour, and this entry does not invent one — an unsourced linkage refuses.
- **The reach is a ceiling, not a rest length.** Dedden's 60 nm relaxed form says a real talin at
  97 nm is at full extension and under load. Modelling that load is B-4's problem.
