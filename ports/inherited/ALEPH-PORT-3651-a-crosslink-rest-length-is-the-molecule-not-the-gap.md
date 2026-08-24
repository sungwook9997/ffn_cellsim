# ALEPH-PORT-3651 — a crosslink's rest length is the molecule, not the gap it was born in

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3651` |
| Lane | `a40e55a2 S-OBSERVE handover, scope extended by PI reassignment` |
| Status | `PROPOSED` |
| Written | `2026-08-08` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Follows | `ALEPH-PORT-3650`, deliberately — that entry's refusal is measured **before** this one moves the same numbers |

---

## 1. Aleph API

`aleph/vertical/cortex_sourcing.py` gains:

- `CrosslinkerRestUnsourcedError(UnsourcedCortexError)`
- `crosslinker_rest_length_um(molecule) -> float` — the molecule's **unstrained** end-to-end
  length, which is what a spring's rest length means. Refuses a molecule this entry has not
  sourced, rather than returning its contour length as a stand-in.

`aleph/vertical/cortex_surface_coupling.py` — the two `Crosslink(...)` constructions at the
`FILAMENT_MIDPOINT` and `closest_approach` branches take `rest_length_um` from that function when
a crosslinker is supplied, instead of from `float(np.linalg.norm(point_b - point_a))`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none — no provider file was read for this entry.** |
| Source commit | **not applicable**; nothing was read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source path | **not applicable** |
| Source symbol(s) | **not applicable**; the Aleph symbols approved here are named in §1 |
| Read from | **primary literature** — Meyer & Aebi 1990 *J Cell Biol*; Gorlin et al. 1990 *J Cell Biol*; Suphamungmee et al. 2012 *J Mol Biol* |
| Working tree == commit? | **not applicable** — no provider revision is cited |

## 3. Why source-derived porting beats clean-room

It does not, and nothing here is source-derived. Recorded explicitly rather than left blank.

## 4. Physical or mathematical law represented

A crosslinker is a spring between two filaments. A spring has three separable lengths and this
repository currently conflates two of them:

- **contour length** — the furthest the molecule can reach. Beyond it there is no bond. This is
  `Crosslinker.capture_radius_um`, supplied by `ALEPH-PORT-3640`.
- **rest length** — the length at which the molecule exerts zero force. **This is what
  `Crosslink.rest_length_um` means**, and it is the field taken from the construction today.
- **the actual separation** — a state variable, not a parameter.

Today the third is assigned to the second. The consequence is that **every crosslink is born
unstrained wherever it happens to be**, so the crosslink population can never be wrong: it reports
zero force at construction by definition, whatever geometry produced it.

**And rest ≠ contour, which is why this entry cannot simply reuse the capture radius.** α-actinin is
a rigid dumbbell rod — Meyer & Aebi measure it at 36 nm bound to F-actin, and that is both its
reach and its working length. Filamin A is explicitly *"a molecular leaf spring"* (Gorlin 1990) with
two hinges, whose 160 nm is a **contour** the relaxed molecule does not span. One rule cannot serve
both, so this entry sources them separately and **refuses any molecule it has not sourced**.

## 5. Units, domains, singular cases, invariants

- All lengths **µm**. `rest_length_um > 0`.
- Domain: `rest_length_um <= capture_radius_um` — a molecule cannot rest longer than it can reach.
- Singular case: a molecule with no sourced rest length. **Refused**, not defaulted to the contour.
- **I1.** `rest_length_um` is independent of ρ. Today the value it replaces is not: median
  `0.53469 µm` at ρ=8 and `0.15376 µm` at ρ=100, i.e. `1.512×` and `1.538×` the mean spacing.
- **I2.** After this entry, two cortices built at different densities with the same molecule carry
  **bit-identical** crosslink rest lengths.
- **I3.** A crosslink is no longer born unstrained: `extension = separation − rest` is generally
  non-zero at construction, and its distribution is a measurement rather than a definition.

## 6. Source evidence class and known retractions

Peer-reviewed EM / cryo-EM structural measurements. No retraction known for any of the three.
Gorlin 1990's "leaf spring" is a structural characterisation from rotary-shadowed micrographs and
sequence analysis, used here only to establish that filamin's rest is **not** its contour — not to
supply a number, which is why filamin is refused rather than given one.

## 7. Independent oracle or derivation

Same oracle as `-3650` §7 and it needs nothing from this tree: **a molecule's rest length may not
move when the density moves.** Build at ρ=8 and ρ=100 and compare. The present construction fails
by `1.512` vs `1.538` × spacing; a sourced constant differs by exactly zero. Separately, `I3` gives
a second oracle: after the change the construction-time extension must be **non-zero for some
crosslink**, because a population that is still born unstrained has not actually changed.

## 7a. What the oracle actually returned, and it is the entry's main result

Measured 2026-08-08, **reach frozen at the molecule and the density moved** — the direction the plan
pre-registers, because moving the reach until the network percolates is how the defect arose.

| ρ /µm² | mean spacing | filaments | crosslinks | links/filament | rest length |
|---:|---:|---:|---:|---:|---|
| 8 | 0.3536 µm | — | **REFUSED** | — | — |
| **100** | 0.1000 µm | — | **REFUSED** | — | — |
| 200 | 0.0707 µm | 62,330 | 15,842 | 0.25 | `0.0360` |
| 300 | 0.0577 µm | 93,495 | 175,868 | 1.88 | `0.0360` |
| 500 | 0.0447 µm | 155,826 | 504,908 | 3.24 | `0.0360` |

**Three findings, and the second is the one to act on.**

1. **`I2` holds.** `rest all == 0.0360` at every density that builds — bit-identical across a 2.5×
   change in ρ, where the rule this replaces would have moved by `√2.5 = 1.58×`.

2. **ρ = 100, the density this repository calls the sourced native cortical density, REFUSES.** With
   a 10 nm filament and α-actinin's 60 nm reach, two centres can be at most 70 nm apart and the mean
   spacing is 100 nm — **1.4× too far**. The refusal is correct arithmetic, not a bug, and it
   localises the defect: turn the inequality around and at ρ = 100 the filament must satisfy
   `L ≥ 0.100 − 0.060 = 0.040 µm`. **That is within 10% of the 36 nm floor `ALEPH-PORT-3650` derived
   from an entirely independent argument** — the F-actin crossover repeat. Two unrelated derivations
   landing together is the reason to believe the **filament**, not the density, is what is wrong.

3. **The threshold matches the arithmetic to 2%.** `spacing ≤ reach + filament` predicts
   `ρ ≥ 1/(0.060+0.010)² = 204`; the builder connects at **200**.

**And one convergence worth recording.** `ALEPH-PORT-3640`'s module docstring notes the old
geometric reach gave **20.19 partners per filament against a physiological 3–6**. With the reach
frozen at the molecule, ρ = 500 gives **3.24** — inside that band, reached by fixing the molecule
and letting the density move rather than by tuning the connectivity directly.

## 8. Positive control

`tests/vertical/test_crosslink_rest_length.py`:

- `test_alpha_actinin_rest_is_the_bound_length` — `crosslinker_rest_length_um` returns `0.036`, the
  value `ACTIN_CROSSOVER_REPEAT_UM` carries, and the two agree by construction rather than by
  literal.
- `test_rest_never_exceeds_reach` — for every sourced molecule, `rest <= capture_radius_um`.
- `test_the_rest_length_is_the_same_at_two_densities` — the function is a property of the molecule
  and takes no density, checked by calling it with the same molecule twice and comparing bits.

## 9. Deliberately failing negative control

- `test_filamin_is_refused_because_a_leaf_spring_is_not_its_contour` — filamin raises
  `CrosslinkerRestUnsourcedError`, and the message says why a 160 nm contour is not a rest length.
  **This control fails if someone later "fixes" filamin by returning 0.160.**
- `test_an_unknown_molecule_is_refused` — a molecule not in the sourced table raises rather than
  falling back.
- `test_the_lookup_can_succeed_and_therefore_can_fail` — vacuity control.
- `test_the_old_rule_would_be_density_dependent` — pins the defect: `1.512` and `1.538` are
  asserted to differ from each other by less than 2% *and* to both be ≈1.5, which is the signature
  of a length read off the spacing. If a future change makes the rest length constant, this control
  is what shows the change did something.

## 10. Numerical and precision envelope

Exact float64 equality where two declared constants are compared (`I2`, and the rest-vs-band
agreement in §8), because both sides are literals and a tolerance would be a place for drift to
hide. `I3` uses a strict `> 0` on the maximum absolute extension, not a tolerance.

## 11. Production-backend residency and transfer

`Crosslink.rest_length_um` already reaches the device through `crosslink_arrays()`, which packs it
into a float64 array consumed by `crosslink_energy_and_forces`. **This entry changes the values in
that array and not its shape, dtype, or the transfer**, so residency is unaffected and the existing
host/device parity controls remain the check on it.

## 12. Comments and docstrings to discard

Nothing is copied, so nothing is discarded. One docstring is **quoted for refutation** and must be
corrected as part of this entry: `crosslink_capture_from_spacing_um` states its reach is *"a
property of the population's density rather than of the mesh that happened to seed it"* — which
correctly rejects mesh-dependence and then substitutes density-dependence, and a molecule is
neither.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass and the I3 distribution is measured and published.
It may not reach `ACCEPTED` while filamin is refused — a cortex model that can only use one
crosslinker is a limitation this entry creates knowingly and must record until a source is read.

## 14. Honest limits

- **Only α-actinin is sourced.** Filamin is refused. A cortex asking for it will not build, and
  that is a deliberate cost of not inventing its rest length.
- **This does not change the reach.** `-3640`'s `capture_radius_um = 0.060` for both presets is
  still one lumped number for two molecules that differ by 4.4×, still with no primary citation for
  that field. Reported in `-3650` §2 and routed to `-3640`'s pin; untouched here.
- **This does not choose the filament length.** V-4 is still a named blank; `-3650` supplies only a
  floor.
- **The crosslink population will carry non-zero energy at construction after this lands.** That is
  `I3` working, and it will move every cortex energy number in the repository. Those numbers were
  measurements of a definition; they are now measurements of a configuration.
