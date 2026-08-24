# ALEPH-PORT-404 — how the (κ, σ) measurement can be defeated

| Field | Value |
|---|---|
| Lane | L4 (analytic oracles) |
| Target | `validation/analytic/degeneracy.py` |
| Reference consulted | `sandbox_fluctuation.py:756-797` (power-law / roll-off generators), `1637-1926` (trusted range, cutoff ladder) |
| Port class | **RE-DERIVED**, scope enlarged (§3) |
| Written | 2026-07-30 |
| Status | LANDED |
| Decision it serves | `ALEPH-DQ-103` — this file is the identifiability statement attached to that decision |

## 1. The three defeats, derived

### D1 — a pure `q⁻²` spectrum is a Helfrich member, so `χ²` can never reject it

Set `κ = 0`. Then `S(q) = kB T/(A σ q²)` **exactly**. A pure `q⁻²` power law with amplitude `C` is
therefore reproduced with zero residual by `(κ, σ) = (0, kB T/(A C))`. It is not "approximately
Helfrich" and it is not "Helfrich in a limit": it is an interior… in fact a *boundary* point of the
model's own parameter set. Consequently

```
min over (κ,σ) of  χ²(q⁻² data)  =  χ²(best fit)  with residual identical to a correct fit
```

and no goodness-of-fit statistic can ever flag it, at any sample size, in any band. Symmetrically,
pure `q⁻⁴` is the member `σ = 0`.

**The only thing that discriminates is the exponent**, because the exponent question is "what is the
slope?" and not "does some `(κ,σ)` reproduce this?". A `χ²` test asks whether the data are *in* the
model; an exponent test asks *where* in the model. Only the second one has an answer here. This is
the single most important sentence in this lane, and (degeneracy suite) `test_chi_square_cannot_reject_a_pure_q2_spectrum`
is its executable form.

Corollary that the tests also assert: `q⁻³` is **not** a member (no `(κ,σ) ≥ 0` gives slope `−3`
except at the single point `q = q*`), so `χ²` *can* reject that one. The model is falsifiable — just
not against its own boundary.

### D2 — amplitude rescaling is exactly absorbed by `(κ,σ) → (κ/f, σ/f)`

```
S(q; κ/f, σ/f) = kB T / ( A ( (κ/f) q⁴ + (σ/f) q² ) )
               = f · kB T / ( A ( κ q⁴ + σ q² ) )
               = f · S(q; κ, σ)                                   for every q, exactly
```

The map is an exact symmetry of the model for all `q` simultaneously — not a near-degeneracy, not a
weakly-constrained direction. It is *one-parameter*: multiplying the measured spectrum by any `f > 0`
lands on another exact member. Therefore:

- `q*= sqrt(σ/κ)` is **invariant** under it (`sqrt((σ/f)/(κ/f)) = sqrt(σ/κ)`), so the crossover is the
  part of the inference that survives an amplitude calibration error.
- The **ratio `σ/κ` is identifiable; the absolute scale of either is not**, unless `kB T / A` is
  independently known. Any claim about an absolute `κ` in pN·μm is a claim about the area calibration
  and the thermostat temperature, silently.
- In the fit's own geometry this is the `(1,1)` direction of `θ = (ln κ, ln σ)`, along which the
  Jacobian identity `f_κ + f_σ = 1` (`ALEPH-PORT-403 §1.2`) makes the model shift rigidly. The
  algebraic identity and the design-matrix degeneracy are the same fact seen twice.

This is a **real limit on `ALEPH-DQ-103`** and is recorded as such, not as a caveat.

### D3 — an unresolved short-wavelength cutoff

Two physically distinct mechanisms, both implemented, because they bias in *opposite* directions and
conflating them is how a fit gets quietly rescued by the wrong story:

**(a) Observation-side box averaging (suppression).** If the height field is reported as an average
over pixels/facets of extent `a`, the observation convolves `h` with a box of width `a`. In Fourier
space that multiplies the *amplitude* by `sinc(qa/2) = sin(qa/2)/(qa/2)` and therefore the
*variance* by

```
G(q) = sinc²(q a / 2)                                            [≤ 1, → 0 at q = 2π/a]
```

`G` decreases with `q`, so the observed spectrum falls *faster* than the truth at high `q`. Fitted
with the continuum model, that reads as extra bending: **`κ` is biased high.**

**(b) Lattice-side finite-difference dispersion (softening).** A nearest-neighbour discrete Laplacian
on spacing `a` has eigenvalue `−(2/a)²sin²(qa/2)`, i.e. the lattice behaves as if `q → q̃`,

```
q̃(q) = (2/a) sin(q a / 2)   ≤ q
```

so the restoring density `κ q̃⁴ + σ q̃²` is *smaller* than the continuum one and the lattice carries
**more** power at high `q`. Fitted with the continuum model, **`κ` is biased low.**

Both are unresolved-cutoff artefacts; both are invisible at low `q`; they have opposite signs. A
future GPU run will have (b) built in by its discretisation and may have (a) from its observation
operator. Aleph therefore records the sign of each and refuses to speak of "the" cutoff damage as a
single multiple — the damage is a *curve* in `q_cutoff/q_max`, which is the form both are reported in.

## 2. What was found in the reference

| Item | Verdict |
|---|---|
| `power_law_variances` pinned through a reference point so the wrong exponent is indistinguishable *at one q* | **Good design.** Carried as an idea, re-implemented. |
| `mesh_truncated_variances` using `exp(−(q/q_c)^p)` | **Defensible as a smooth stand-in, but it is not derived from anything.** It is an ad-hoc roll-off with a free exponent `p`, whereas `sinc²` (a) and the lattice dispersion (b) are *derivable* from a stated observation model. Aleph implements the two derivable ones and keeps the ad-hoc form only as an explicitly-labelled generic stress test. |
| Reference treats the cutoff only as *suppression* (`κ` inflation) | **INCOMPLETE.** The finite-difference case biases `κ` the *other* way. The reference's `measure_cutoff_damage_ladder` refuses to compute a ratio when `κ_true = 0` and reports a ladder rather than a single number — both good — but the sign asymmetry is never raised. |
| `CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE` scope note | **Correct principle.** Adopted. |
| The `q⁻²`-is-a-member degeneracy | **NOT PRESENT as a stated degeneracy** in what I read. The reference has a `WRONG_EXPONENT` failure route and an exponent test, so it *behaves* correctly, but it does not state that `χ²` is structurally incapable here. Aleph states and tests it. |
| The amplitude-rescaling degeneracy | **NOT PRESENT.** The reference has a related note (`EQUIPARTITION_SCALE_DEGENERACY`) about `(κ,σ)` needing to be independent of the amplitudes used for the `kB T` check, which is the *same fact in a different costume*, but the exact `(κ/f, σ/f)` symmetry is not derived. Aleph derives and tests it. |

## 3. Controls shipped
> All controls live under `tests/validation/`, one file per target module, named
> `test_<module>.py`. The rows below name the suite in words and then the test function,
> rather than giving a path: the ledger index scanner reads any bare `validation/...py`
> substring as a claim about a source file at the repository root, and reads any
> `test_<something>` token as a test function that must exist. A literal file path would
> trip both. The function names are exact and are what the discipline test resolves.


| Control | Location | Asserts |
|---|---|---|
| Positive | (degeneracy suite) `test_rescaling_is_an_exact_symmetry` | max rel. deviation `< 1e-14` over 6 decades, for six factors |
| Positive | (degeneracy suite) `test_crossover_is_invariant_under_rescaling` | `q*` and `σ/κ` unchanged to machine precision |
| Positive | (degeneracy suite) `test_pure_power_laws_at_minus_two_and_minus_four_are_exact_members` | exact member `(κ,σ)` recovered |
| Positive | (degeneracy suite) `test_other_exponents_are_not_members` | six other exponents return `None` |
| Positive | (degeneracy suite) `test_box_average_transfer_has_the_derived_expansion` | the `1 − (qa)²/12` expansion, and `G(2π/a) = 0` |
| Positive | (degeneracy suite) `test_lattice_dispersion_has_the_derived_expansion` | the `1 − (qa)²/24` expansion, and `q̃ ≤ q` |
| Positive | (degeneracy suite) `test_a_far_away_cutoff_does_no_measurable_damage` | the low end of the damage curve is nothing |
| Negative (must fail) | (degeneracy suite) `test_chi_square_cannot_reject_a_pure_q2_spectrum` | exact fit on noiseless data; over 120 sampled repeats the α=0.01 gate fires at background rate |
| Negative (must fail) | (degeneracy suite) `test_only_the_exponent_test_rejects_it` | exponent CI excludes `−4` at >20σ on the same data the χ² accepted |
| Negative (must fail) | (degeneracy suite) `test_rescaling_is_invisible_to_chi_square` | rescaled data fit with identical `χ²` and `κ` wrong by exactly `f` |
| Negative (must fail) | (degeneracy suite) `test_an_area_error_is_exactly_an_amplitude_error` | an area miscalibration is the same defeat, wearing the costume it will actually wear |
| Negative (must fail) | (degeneracy suite) `test_box_average_cutoff_biases_kappa_high` | sign and monotonicity of D3(a) |
| Negative (must fail) | (degeneracy suite) `test_cutoff_damage_is_detectable_by_chi_square_only_with_enough_samples` | the bias is always present; detectability is not |
| Negative (must fail) | (degeneracy suite) `test_lattice_dispersion_biases_kappa_low` | opposite sign of D3(b) |
| Negative (must fail) | (degeneracy suite) `test_the_two_mechanisms_transfer_power_in_opposite_directions` | the two transfers bracket 1 from opposite sides |
| Negative (must fail) | (degeneracy suite) `test_a_minus_three_spectrum_is_rejected_by_chi_square` | proves the χ² gate *can* fire, so D1 is a statement about the boundary and not about a dead test |

## 4. Compliance with PLAN §0.2

1. Stripped. ✔ 2. Re-derived (§1). ✔ 3. New prose. ✔ 4. Controls. ✔ 5. Ledger first. ✔
6. No `ffn_cellsim` import. ✔
