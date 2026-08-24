# ALEPH-PORT-3658 — the sourced length crossed the overlap threshold, and the default did not follow

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3658` |
| Lane | `b4fad06b` |
| Status | `PROPOSED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Follows | `ALEPH-PORT-3652` (which sourced the length) and `ALEPH-PORT-3636` (which built the sweep) |
| Found by | timing the sourced whole cell for `docs/design/ENGINE_PERFORMANCE_TARGET.md`. The first step's residual read **2.3 × 10¹⁴ pN**, which is not a number a cell has |

---

## 1. Aleph API

No new function and no new argument. What changes is **one default's derivation**, and it is a
derivation that already exists in the tree and was simply never re-run when its input moved:

`aleph/vertical/cortex_surface_coupling.build_radial_cortex_network`, argument
`overlap_relaxation_sweeps`. Its own docstring states the condition:

> *"Any positive value settles the population to steric contact before it is handed to the network,
> which is what a density above `ρ_crit = 1/(2·filament_half_length_um)²` requires; below that
> density it changes nothing, because there are no overlaps to clear."*

**The condition is correct. The default was chosen when the input satisfied it, and the input
moved.**

| | `filament_half_length_um` | `ρ_crit = 1/(2·half)²` | native ρ = 100 |
|---|---:|---:|---|
| shipped, before `-3652` | 0.005 | **10,000 /µm²** | far **below** — no overlaps, `sweeps = 0` correct |
| sourced, after `-3652` | 0.060 | **69.4 /µm²** | **above** — overlaps exist, `sweeps = 0` wrong |

**A 12× change in the filament length is a 145× change in the overlap threshold**, because the
threshold goes as the length squared, and it crossed the density this repository calls native.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** — no provider commit is behind this |
| Source path and symbol | **not applicable** — no `ffn_sim` path and no provider symbol was read |
| Read from | this tree's own `resolve_cross_filament_overlaps` docstring and `ALEPH-PORT-3636`; the length is `ALEPH-PORT-3652`'s |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. **RE-DERIVED**: this entry re-runs an existing derivation
against an input that changed.

## 4. Physical or mathematical law represented

Hard-core exclusion between tangent rods on a shell. Filament centres at areal density `ρ` sit a
mean `1/√ρ` apart; a rod of full length `2h` laid tangentially overlaps its neighbour once
`1/√ρ < 2h`, i.e. once

```text
ρ  >  ρ_crit  =  1 / (2h)²
```

Below the threshold the relaxation is a no-op by construction. Above it, the builder hands the
network a population that is **interpenetrating**, and the steric law is then evaluated at
separations it was never meant to see.

**No new physics and no new constant.** The law, the threshold and the sweep all already exist; this
entry is the observation that one default no longer satisfies the condition its own docstring states.

## 5. Units, domains, singular cases, invariants

- **I1 — below `ρ_crit`, the sweep changes nothing, bit for bit.** This is the property that makes
  the change safe for every existing caller: at the shipped 10 nm filament `ρ_crit` is 10,000 /µm²
  and no density in this repository approaches it. **Every committed cortex number is therefore
  untouched**, and that must be asserted rather than assumed.
- **I2 — the sweep terminates, and the count is not a tuning knob.** It runs out of overlaps, not
  out of iterations.
- **I3 — the relaxation is applied BEFORE the network is built.** Crosslink rest lengths and the
  membrane's material sites are read off the node coordinates, so a shell settled afterwards would
  carry crosslinks that remember where the filaments used to be. The existing code already does
  this and the ordering is load-bearing.
- **I4 — the sweep moves nodes and therefore changes the crosslink set.** It is not energy-neutral
  and does not claim to be; it replaces a configuration the steric law refuses with one it accepts.

## 6. Source evidence class and known retractions

No literature is added. The evidence class is **measurement on this tree**. The filament length it
depends on is `ALEPH-PORT-3652`'s, peer-reviewed and read in full; the threshold is arithmetic.

## 7. Independent oracle or derivation

**Written before the code, three of them, and two are already measured.**

1. **The threshold is crossed and it is arithmetic, not opinion.** `1/(2 × 0.060)² = 69.44 /µm²`
   against a native `ρ = 100`. **44 % above.** At the shipped length, `1/(2 × 0.005)² = 10,000`, and
   100 is **100× below**. No measurement needed to know which side of the line each sits on.
2. **A population above the threshold produces a force no cell has.** *Measured before this entry
   was written*, whole cell at ρ = 100 with the sourced 120 nm filament:

   | `cortex_overlap_relaxation_sweeps` | build | `|F|max` at build | energy |
   |---:|---:|---:|---:|
   | **0** | 32.5 s | **5.858 × 10¹⁶ pN** | 2.925 × 10¹² pN·µm |
   | **2** | 32.4 s | **1.912 × 10⁴ pN** | 2.681 × 10⁶ pN·µm |
   | 8 | 32.2 s | 1.912 × 10⁴ pN | 2.681 × 10⁶ pN·µm |
   | 32 | 32.1 s | 1.912 × 10⁴ pN | 2.681 × 10⁶ pN·µm |

   **A factor of 3.1 × 10¹² in force and 1.1 × 10⁶ in energy, for no measurable build time.**
3. **The closest approach crosses the steric diameter, and that is the SCALE-FREE form of (2).**
   The energy in (2) is set by the single closest pair against a steep law, so it is an
   extreme-value quantity: it explodes with the filament **count**, not with the density.
   Measured on the same construction at two shell sizes —

   | shell | filaments | `E` unrelaxed | `E` 2 sweeps | ratio |
   |---|---:|---:|---:|---:|
   | 1 µm | 1,207 | 1.428e5 | 1.094e5 | **1.31×** |
   | 2 µm | 4,927 | 6.485e5 | 4.328e5 | 1.50× |
   | 5 µm (the cell) | 31,165 | 2.925e12 | 2.681e6 | **1.1e6×** |

   **So the energy ratio is not a portable oracle and this entry does not use it as one.** The
   closest cross-filament approach, in units of the card's steric diameter, is: measured at the 1 µm
   shell, **0.2183× unrelaxed → 1.1337× after two sweeps.** Below 1 the population is
   interpenetrating at any scale; above 1 it is not.

4. **It terminates in two sweeps, which the docstring predicted before it was run.**
   `resolve_cross_filament_overlaps`: *"it terminates in **2 sweeps at every density**: it runs out
   of things it can see."* **8 and 32 sweeps are identical to 2 in every printed digit** — so the
   sweep count is not a knob that buys more, and a caller cannot tune the cortex through it. **If
   32 had differed from 2, this entry's premise would be wrong**, because a "relaxation" that never
   finishes is a parameter and not a repair.

## 8. Positive control

`tests/vertical/test_overlap_threshold.py`:

- `test_the_threshold_is_the_documented_one` — `ρ_crit = 1/(2·half)²` at both the shipped and the
  sourced half-length, so the 145× move is pinned as arithmetic rather than as prose.
- `test_two_sweeps_clear_the_interpenetration_above_the_threshold` — §7(3), on the closest
  approach rather than on the energy, so the assertion is portable across shell sizes. The energy
  form was tried first and **failed on the small shell at 1.31× against a cell-scale 1.1e6×**,
  which is how §7(3)'s extreme-value table came to be measured.
- `test_the_sweep_terminates` — 2, 8 and 32 sweeps agree bit for bit above the threshold.
- `test_below_the_threshold_the_sweep_changes_nothing` — **I1**, the compatibility claim, at the
  shipped 10 nm filament where every committed number was measured. Asserted with
  `np.array_equal` on node positions, because the claim is bit-identity.

## 9. Deliberately failing negative control

- `test_the_unrelaxed_population_really_is_interpenetrating` — a **vacuity control**. Without it
  every control above would pass against a builder that had no steric law at all: it asserts that at
  ρ = 100 with the sourced length there exist node pairs closer than the card's steric diameter, so
  there is genuinely something for the sweep to remove.
- `test_a_negative_sweep_count_is_refused`.

## 10. Numerical and precision envelope

float64. The compatibility claim (I1) is `np.array_equal`, exact. The force comparisons above the
threshold are order-of-magnitude assertions (`< 1e6 pN` against `> 1e10 pN`) rather than tolerances,
because a 10¹²-fold separation does not need one and a tight tolerance there would be a pin on a
number nobody derived.

## 11. Production-backend residency and transfer

**None.** `resolve_cross_filament_overlaps` runs at build, on the host, before any backend array is
formed. Node **positions** change, array shapes do not, and the crosslink count moves with the
settled geometry — which the existing host/device parity controls already cover.

## 12. Comments and docstrings to discard

Nothing copied. `build_radial_cortex_network`'s `overlap_relaxation_sweeps` docstring gains the
worked threshold at both filament lengths, because a condition stated symbolically was not enough to
stop this: the symbol was right and nobody evaluated it when its input moved.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass. **`ACCEPTED` when the whole-cell timing run at the
sourced construction is repeated with the sweep on and the before/after is recorded** in
`docs/design/ENGINE_PERFORMANCE_TARGET.md` — because the first ~25 steps of every sourced run so far
were spent undoing this, at a step size of `dt ≈ 3.4 × 10⁻²⁰ s`, and that is time the target's
arithmetic charged to the solver.

**The library default is NOT moved by this entry.** `overlap_relaxation_sweeps = 0` remains the
retained control on `build_radial_cortex_network`, for `-3652` §13's rule and because I1 says it is
harmless below the threshold — which is every historical caller. What changes is that **a caller
above the threshold is told**, and the `sourced` construction in `scripts/engine_clock_run.py`
turns the sweep on because it is above it.

## 14. Honest limits

- **A refusal would be stronger than a warning, and is not built here.** The honest object is
  `build_radial_cortex_network` refusing a density above `ρ_crit` with `sweeps = 0`, the way
  `-3650` refuses a cortex with no crosslinker. That is a behaviour change with a blast radius
  across `scripts/_ws_*`, and it belongs to its own entry.
- **The settled population is not a measured cortex.** `resolve_cross_filament_overlaps` pushes
  node pairs apart along their separation; `OverlapResolution.NODE_PAIR`, the default, is recorded
  in `ALEPH-PORT-3644` as **known to leave interpenetrations it cannot detect**. So this removes the
  10¹² pN artefact; it does not certify the result is overlap-free.
- **The 2-sweep termination is a property of what the resolver can see**, not proof that no overlap
  remains — same caveat, and it is why §7(3) is framed as "the count is not a knob" rather than as
  "the population is clean".
- **This changes crosslink counts at the sourced density**, because the centres move before the
  network is built. The `-3657` connectivity numbers were taken with `sweeps = 0` and are not
  re-measured here.
