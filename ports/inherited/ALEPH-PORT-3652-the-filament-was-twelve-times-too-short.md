# ALEPH-PORT-3652 — the cortical filament was twelve times too short, and that was the whole refusal

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3652` |
| Lane | `a40e55a2 S-OBSERVE handover, scope extended by PI reassignment` |
| Status | `PROPOSED` |
| Written | `2026-08-08` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Closes | the named blank **V-4** in `docs/design/SOURCED_CONSTRUCTION_VALUES.md`, and the refusal `ALEPH-PORT-3651` §7a measured |

---

## 1. Aleph API

`aleph/vertical/cortex_sourcing.py`:

- `CORTICAL_ACTIN_FILAMENT_LENGTH_UM: LiteratureRange` — `[0.060, 0.120]` µm, the Arp2/3-nucleated
  population across the two cell lines the source measures.
- `sourced_filament_half_length_um(cell_line: str = "hela") -> float` — half the sourced length,
  because `build_vertical` parameterises by half-length. Refuses an unsourced cell line.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none — no provider file was read.** |
| Source commit | **not applicable**; nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source path | **not applicable** |
| Source symbol(s) | **not applicable**; the Aleph symbols are named in §1 |
| Read from | **the primary paper, in full**, `https://discovery.ucl.ac.uk/1492897/10/Actin%20kinetics%20shapes%20cortical.pdf` — not an abstract, not a citing paper |
| Working tree == commit? | **not applicable** — no provider revision is cited |

Fritzsche M., Erlenkämper C., Moeendarbary E., Charras G., Kruse K. (2016) *Actin kinetics shapes
cortical network structure and mechanics.* **Science Advances 2:e1501337**, 154 citations.

## 3. Why source-derived porting beats clean-room

It does not; nothing here is source-derived from the provider. Stated explicitly.

## 4. Physical or mathematical law represented

Not a law — a **measured length**, and the one every other cortical length has to be commensurate
with. Verbatim from §Results:

> *"In HeLa cells, we calculated that formin-nucleated filaments had an average length of 1200 nm,
> whereas the average length of Arp2/3-nucleated filaments was 120 nm (M2 cells: 600 and 60 nm)...
> The latter was consistent with results from electron tomography and numerical simulations of the
> lamellipodium."*

and

> *"less than 10% of the filaments, containing 20 to 25% of actin protomers, were nucleated by
> formins."*

## 5. Units, domains, singular cases, invariants

- **µm** throughout. `build_vertical` takes a **half**-length; the source reports full lengths.
- **I1.** `2 × half_length ≥ ACTIN_CROSSOVER_REPEAT_UM` — `-3650`'s floor. Measured: `0.120` against
  `0.036`, clear by **3.3×**.
- **I2.** `reach + filament ≥ mean spacing` at the density being built, or no crosslink exists.
  With `0.060 + 0.120 = 0.180 µm` the percolation threshold is `ρ ≥ 1/0.180² = 31 /µm²`.
- Singular case: an unsourced cell line. Refused rather than defaulted.

## 6. Source evidence class and known retractions

Peer-reviewed, open access, single-molecule fluorescence imaging combined with stochastic
simulation, and the Arp2/3 figure is stated in the paper to be **independently consistent with
electron tomography**. No retraction known. The 154 citations are recorded as an indication of
scrutiny, not as evidence.

## 7. Independent oracle or derivation

**The engine itself, asked a question it could previously only answer with a refusal.**
`ALEPH-PORT-3651` §7a measured that ρ = 100 refuses with the shipped 10 nm filament. If 120 nm is
the right length, that same build must now succeed — and ρ = 8 must still refuse, because
`0.180 µm` reachable against a `0.354 µm` spacing is genuinely too sparse. Both halves are
predictions made before the build was run.

**Measured 2026-08-08:**

| ρ | filament | filaments | crosslinks | links/filament | rest |
|---:|---:|---:|---:|---:|---|
| 100 | 0.010 (shipped) | — | **REFUSED** | — | — |
| **100** | **0.120 (sourced)** | **31,165** | **42,882** | **1.38** | `0.0360` |
| 8 | 0.120 (sourced) | — | **REFUSED**, as predicted | — | — |

The threshold arithmetic moves `204 → 31`, and the native density clears it by 3.2×.

## 8. Positive control

`tests/vertical/test_cortex_sourcing.py`:

- `test_the_sourced_filament_clears_the_crossover_repeat_floor` — `0.120 > 0.036`, and it is read
  from the two bands rather than from literals, so the numbers cannot drift apart.
- `test_the_sourced_filament_percolates_at_the_native_density` — `reach + filament ≥ 1/√100`.
- `test_hela_and_m2_are_both_in_the_band`.

## 9. Deliberately failing negative control

- `test_the_shipped_default_does_not_percolate_at_the_native_density` — pins the defect: with
  `2 × 0.005` the same inequality **fails**, and it fails by 1.43×.
- `test_an_unsourced_cell_line_is_refused` — `sourced_filament_half_length_um("hamster")` raises
  rather than returning the HeLa value.
- `test_the_accessor_can_succeed_and_therefore_can_fail` — vacuity control.

## 10. Numerical and precision envelope

Exact float64 comparisons between declared constants. The percolation inequality is evaluated on
`1/√ρ` with no tolerance; it is not near its boundary at either density tested (3.2× clear at
ρ = 100, 2.0× short at ρ = 8), so no envelope is needed.

## 11. Production-backend residency and transfer

**None.** A build-time scalar. It changes how many filaments and crosslinks are created, which
changes array **sizes** on the device — but not their dtype, layout, or transfer, and the existing
host/device parity controls remain the check.

## 12. Comments and docstrings to discard

Nothing copied, nothing discarded. The paper is quoted verbatim in §4 for refutation and cited.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass. **It may not reach `ACCEPTED` while
`cortex_filament_half_length_um` still defaults to `0.005`** — this entry makes the right value
available and sourced; changing a default that every cortex result in the repository was taken at
is a separate, measured step and is named in §14.

## 14. Honest limits

- **The default is NOT changed here.** `build_vertical` still ships `0.005`. Every cortex number in
  this repository was taken at a filament 12× shorter than the one now sourced, and moving the
  default moves all of them at once. That deserves its own before/after measurement rather than
  being folded into the entry that sourced the value.
- **One population, not two.** The cortex is a mixture: Arp2/3 at 120 nm carrying >90% of filaments
  by number, formin at 1,200 nm carrying <10% by number but 20–25% of protomers, and the source
  calls the formin population *"important determinants of cortical elasticity"*. Aleph has one
  filament-length parameter. Using the number-dominant population is a modelling choice and the
  omission is this entry's largest.
- **HeLa, not "a cell".** 120 nm is HeLa; melanoma M2 is 60 nm, a 2× spread between two lines. The
  band carries both; the accessor's default carries HeLa and says so.
- **It does not fix the placement.** Six compartments are still first-octant blobs
  (`PROPOSAL-only-the-cortex-has-a-resolution.md` §4). A correctly-sized filament in a blob is a
  correctly-sized filament in a blob.
