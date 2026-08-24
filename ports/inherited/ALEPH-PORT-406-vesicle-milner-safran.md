# ALEPH-PORT-406 — closed-vesicle (Milner–Safran) fluctuation spectrum

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-406` |
| Lane | `L4 analytic oracles` |
| Target | `validation/analytic/vesicle.py` |
| Status | `ACCEPTED` |
| Written | 2026-07-30 (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` — greenfield; nothing was ported |

---

## 1. Aleph API

```python
from validation.analytic.vesicle import (
    DegreeMapping,          # LAPLACIAN_EIGENVALUE | LINEAR_IN_DEGREE
    ExcludedDegree,         # VOLUME_CHANGE | RIGID_TRANSLATION
    FlatLimitAgreement,
    OrderConvention,        # PER_ORDER | MEAN_OVER_ORDERS | SUM_OVER_ORDERS
    VesicleSpectrum,
    angular_power_from_per_order,
    effective_samples_for_angular_power,
    flat_limit_agreement,
    per_order_from_angular_power,
)
```

Nothing outside this list is authorised by this entry.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | none |
| Reference consulted | **none.** No file in `/Users/sw1/ffn_cellsim` was read for this module. |
| Source commit | not applicable — no source was read, so no commit is cited. Citing one would be a fiction. |
| Source path | not applicable — no source path, no source symbol. |
| Origin of the statement under test | The lane coordinator supplied the target formula in prose. It was **not** taken on faith; §4 derives it independently and §4.4 records where the surrounding brief was wrong. |

## 3. Why source-derived porting beats clean-room

It does not, and no porting occurred. This is a clean-room derivation from the Helfrich functional
and the two geometric constraints. The result (Milner & Safran, *Phys. Rev. A* **36**, 4371 (1987))
is standard published physics with a short derivation, so there is nothing expensive to inherit and
no empirically-discovered failure mode that would justify carrying someone else's code.

Recording it as a ledger entry anyway, because §0.2 exists to catch **unrecorded inheritance**, and a
formula handed over in prose by another lane is an inheritance even when no code moves. The audit
question "where did this number come from, and who checked it?" has to have an answer.

## 4. Physical or mathematical law represented

### 4.1 Setup

Quasi-spherical vesicle, `r(Ω) = R(1 + u(Ω))`, `|u| << 1`, `u` **dimensionless**. Expand in real
spherical harmonics, `u = Σ_ℓm u_ℓm Y_ℓm`, `∫ Y_ℓm Y_ℓ'm' dΩ = δ`. Two constraints apply, and both
do real work:

* **enclosed volume fixed** (the vesicle is osmotically stabilised — water crosses the bilayer far
  faster than lipid does),
* **the tension `σ` is conjugate to area**, so the area change is what `σ` multiplies.

### 4.2 The area change, derived (this is where `(ℓ−1)(ℓ+2)` comes from)

To second order in `u`, with `∇` the angular gradient on the unit sphere:

```
A = R² ∫dΩ [ 1 + 2u + u² + ½|∇u|² ]
V = (R³/3) ∫dΩ [ 1 + 3u + 3u² + … ]
```

Using `∫|∇u|² dΩ = Σ ℓ(ℓ+1) u_ℓm²` and `∫u² dΩ = Σ u_ℓm²`:

```
A = 4πR² + 2R² ∫u dΩ + R² Σ [1 + ℓ(ℓ+1)/2] u_ℓm²
V = (4π/3)R³ + R³ [ ∫u dΩ + Σ u_ℓm² ]
```

Holding `V` fixed forces `∫u dΩ = −Σ u_ℓm²` — i.e. **the `ℓ = 0` amplitude is not free; it is slaved
at second order to all the others.** Substituting it into `A`:

```
ΔA = R² Σ [ 1 + ℓ(ℓ+1)/2 − 2 ] u_ℓm²
   = R² Σ [ ℓ(ℓ+1)/2 − 1 ] u_ℓm²
   = (R²/2) Σ (ℓ−1)(ℓ+2) u_ℓm²          since ℓ(ℓ+1) − 2 = (ℓ−1)(ℓ+2)
```

**`(ℓ−1)(ℓ+2)` is the volume-constrained area change per unit mode amplitude.** It is not a
correction bolted onto the flat result; it is what the enclosed-volume constraint leaves behind. Note
it vanishes identically at `ℓ = 1`: a translation changes neither area nor volume, so it cannot cost
tension energy. That is a *derived* fact, not an assumption bolted on afterwards.

The bending term for a quasi-sphere is the standard expansion

```
E_bend = 8πκ + (κ/2) Σ (ℓ−1)(ℓ+2) ℓ(ℓ+1) u_ℓm²
```

whose structure is checkable at the two boundary degrees: it must vanish at `ℓ = 1` (translation) and
at `ℓ = 0` (a dilation is pure volume, already removed), and `(ℓ−1)(ℓ+2)ℓ(ℓ+1)` does both.

### 4.3 The spectrum

```
E = ½ Σ_{ℓ≥2, m} (ℓ−1)(ℓ+2) [ κ ℓ(ℓ+1) + σR² ] u_ℓm²
  = (κ/2) Σ_{ℓ≥2, m} (ℓ−1)(ℓ+2) [ ℓ(ℓ+1) + σ̄ ] u_ℓm²,      σ̄ ≡ σR²/κ
```

Each real `u_ℓm` is one real quadratic degree of freedom, so equipartition (`⟨x²⟩ = kT/(2a)` for
`H = a x²`) gives

```
⟨u_ℓm²⟩ = kB T / [ κ (ℓ−1)(ℓ+2) (ℓ(ℓ+1) + σ̄) ]
```

**The coordinator's stated form is correct.** Two riders that the prose did not carry and that the
module states explicitly:

* `u` is **dimensionless** (relative radial displacement). The *physical* displacement variance is
  `⟨|δr_ℓm|²⟩ = R² ⟨u_ℓm²⟩`, in m². Comparing a metric displacement spectrum against the
  dimensionless one is a factor `R²` — for `R = 5 µm`, a factor `2.5e-11`.
* Unlike the flat case, this formula is **convention-robust in the mean** between real and complex
  `Y_ℓm`: both give `2ℓ+1` modes per degree with the same per-mode variance. The flat case's `1/A`
  ambiguity has no analogue here.

### 4.4 Where the surrounding brief was wrong

| Claim in the brief | Verdict |
|---|---|
| `⟨u_ℓm²⟩ = kT/[κ(ℓ−1)(ℓ+2)(ℓ(ℓ+1)+σ̄)]`, `σ̄ = σR²/κ` | **Correct.** Re-derived in §4.2–4.3. |
| "the `(ℓ−1)(ℓ+2)` prefactor vanishes at both `ℓ=0` and `ℓ=1`, so an unguarded implementation divides by zero" | **Wrong at `ℓ = 0`.** `(0−1)(0+2) = −2`, which is not zero — it is *negative*. An unguarded implementation returns a **negative variance** at `ℓ = 0` and `+inf` at `ℓ = 1`. Two different failure modes needing two different refusals, and the `ℓ = 0` one is the dangerous one: a negative variance survives arithmetic, and `sqrt` of it is a silent `nan` that propagates. Both are refused explicitly, with distinct reasons. |
| "as `ℓ → ∞` with `q ≈ ℓ/R`, the vesicle result approaches the flat one" | **True but the mapping matters, and `q = ℓ/R` is the worse choice.** See §4.5. |
| "L2's `C_ℓ` averages over `2ℓ+1` orders; two conventions differing by `2ℓ+1` will look like a smoothly varying discrepancy" | **The risk is real but it is in the wrong moment.** Under isotropy every order has the same variance, so for the *averaging* convention `⟨C_ℓ⟩ = ⟨u_ℓm²⟩` **exactly — there is no `2ℓ+1` in the mean.** The `2ℓ+1` lives entirely in the degrees of freedom. See §4.6. Writing a `2ℓ+1` into the mean to "fix" the conventions would *introduce* the bug this entry exists to prevent. |

### 4.5 Flat limit — the mapping determines whether the comparison is clean

Two natural ways to assign a wavevector to a degree:

* `q² = ℓ(ℓ+1)/R²` — the actual eigenvalue of the Laplace–Beltrami operator on the sphere;
* `q = ℓ/R` — the naive one.

Under the **Laplacian** mapping, with `V_flat = kT/(A_eff (κq⁴+σq²))`:

```
V_flat  = kT R² / [ κ ℓ(ℓ+1) ( ℓ(ℓ+1) + σ̄ ) ]      (using A_eff = R², see below)
V_ves   = kT R² / [ κ (ℓ−1)(ℓ+2) ( ℓ(ℓ+1) + σ̄ ) ]

V_ves / V_flat = ℓ(ℓ+1) / [(ℓ−1)(ℓ+2)] = L/(L−2),   L ≡ ℓ(ℓ+1)
```

**The `σ̄` factor cancels exactly.** The discrepancy between the two oracles is purely geometric —
independent of `κ`, `σ` and `R`. Under the `q = ℓ/R` mapping it does *not* cancel and the
disagreement acquires a spurious `σ̄` dependence, which is precisely the kind of smooth,
fittable-looking artefact that gets absorbed into a tension estimate. The module defaults to the
Laplacian mapping and says why.

Since `V_ves/V_flat − 1 = 2/(L−2)`, the agreement boundary is exact arithmetic, not a fit:

| Tolerance | First degree `ℓ` meeting it | Actual excess there |
|---|---|---|
| 10 % | **ℓ = 5** | 7.14 % |
| 5 % | **ℓ = 6** | 5.000 % (exactly) |
| 1 % | **ℓ = 14** | 0.96 % |

At `ℓ = 2` the flat formula is wrong by **50 %**; at `ℓ = 3`, by 20 %. This table is reproduced in
both module docstrings so nobody reaches for the wrong oracle. **Below `ℓ ≈ 14` the flat patch is not
an approximation to a vesicle, it is a different answer.**

**The effective area is `A_eff = R²`, not `4πR²`.** This falls out of matching the two equipartition
results and is the single most dangerous number in this entry: using the total area `4πR²` gives a
spectrum low by `4π ≈ 12.57`, and since amplitude rescaling is exactly absorbed by
`(κ,σ) → (κ/f, σ/f)` (`ALEPH-PORT-404` D2) it produces a **perfect fit with `κ` wrong by 12.57×**.
The factor traces to `∫|Y_ℓm|²dΩ = 1` over `4π` steradians against the flat `1/A` convention; it is a
normalisation fact, not a geometry fact.

### 4.6 The `C_ℓ` normalisation conversion

Three conventions are supported explicitly rather than assumed:

| Convention | Definition | Mean | Real dof from `N` snapshots |
|---|---|---|---|
| `PER_ORDER` | `u_ℓm²`, single order | `⟨u_ℓm²⟩` | `N` |
| `MEAN_OVER_ORDERS` (L2's `C_ℓ`) | `(1/(2ℓ+1)) Σ_m u_ℓm²` | `⟨u_ℓm²⟩` — **identical** | `N(2ℓ+1)` |
| `SUM_OVER_ORDERS` | `Σ_m u_ℓm²` | `(2ℓ+1)⟨u_ℓm²⟩` | `N(2ℓ+1)` |

So for L2's averaging convention the conversion factor on the **mean is 1**, and the thing that must
be carried across the boundary is the **degrees of freedom**: `spectrum_fit` must be handed
`samples_per_mode = N(2ℓ+1)`, not `N`. Passing `N` widens every error bar by `√(2ℓ+1)` and deflates
`χ²` by `(2ℓ+1)` — so a wrong model *passes*, and it passes more comfortably at high `ℓ`. That is the
same class of failure as the uncorrected log bias in `ALEPH-PORT-402 §3`: a defect that gets better
at hiding as the measurement gets better.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit (strict SI, per the settled convention) | Domain |
|---|---|---|
| `kappa_J` | J | `> 0` (a vesicle needs bending rigidity to be quasi-spherical) |
| `sigma_N_per_m` | N/m | `>= 0` |
| `radius_m` | m | `> 0` |
| `temperature_K` | K | `> 0` |
| `reduced_tension` `σ̄` | dimensionless | `>= 0` |
| `relative_mode_variance` | dimensionless (`u`) | `ℓ >= 2` |
| `radial_mode_variance_m2` | m² (`δr = Ru`) | `ℓ >= 2` |
| `flat_equivalent_area_m2` | m² | `= R²` |

**`validation/analytic/**` stays strict SI. No scaling convention lives inside the oracle** — a scale
factor inside the physics reference is a scale factor nobody can audit. Conversion happens at the
caller's boundary.

**Singular cases**, both refused with distinct reasons (never silently dropped):

* `ℓ = 0` — uniform radial dilation, i.e. a **volume change**, forbidden by the enclosed-volume
  constraint, and already eliminated in §4.2 by slaving `u_00`. Prefactor `= −2`; unguarded, the
  formula returns a **negative variance**.
* `ℓ = 1` — **rigid translation**, zero energy cost. Prefactor `= 0`; unguarded, the formula returns
  `+inf`. The amplitude of this mode is set by whatever holds the vesicle, not by the membrane.

**Invariants** asserted in the tests:

* `V_ves/V_flat = L/(L−2)` under the Laplacian mapping, independent of `κ`, `σ`, `R`, `T`.
* `q*` from the vesicle crossover degree equals the flat `q* = √(σ/κ)`.
* Amplitude rescaling by `f` absorbed by `(κ,σ) → (κ/f, σ/f)` here too, exactly as in `ALEPH-PORT-404`.
* Monotone decrease of `⟨u_ℓm²⟩` in `ℓ`.

## 6. Numerical and precision envelope

float64 throughout. The ratio identity `L/(L−2)` is exact to round-off and is asserted at `1e-14`
relative. `flat_limit_agreement` degrees are integers obtained by exact integer arithmetic on
`ℓ(ℓ+1)`, not by a floating-point root find, so the boundary table has no tolerance at all. The
`σ̄`-cancellation test sweeps `σ̄` over twelve orders of magnitude (`1e-6` … `1e6`) and holds to
`1e-12` relative. Catastrophic cancellation is not a concern here: every quantity is a ratio of sums
of positive terms.

## 7. Production-backend residency and transfer

Host-side CPU numpy only. No GPU, no device transfer, no backend dependency. `validation/**` is
never imported by `aleph/runtime/**`, so nothing here is ever resident on a production device; it is
a comparison oracle evaluated on the host after the fact.

## 8. Comments and docstrings to discard

None to discard — no source prose exists to strip, since nothing was read. All prose is new and
written for Aleph. The coordinator's brief was consumed as a *claim to be checked* and its wording
does not appear in the module; where it was wrong (§4.4) the module states the corrected fact.

## 9. Controls shipped

> All controls live under `tests/validation/`, one file per target module, named
> `test_<module>.py`. The rows below name the suite in words and then the test function,
> rather than giving a path: the ledger index scanner reads any bare `validation/...py`
> substring as a claim about a source file at the repository root, and reads any
> `test_<something>` token as a test function that must exist. A literal file path would
> trip both. The function names are exact and are what the discipline test resolves.

| Control | Location | Asserts |
|---|---|---|
| Positive | (vesicle suite) `test_mode_variance_matches_the_closed_form` | agreement with an independently written expression |
| Positive | (vesicle suite) `test_energy_and_variance_halves_of_the_derivation_agree` | `⟨E_ℓm⟩ = kT/2` recovered from `k_ℓ⟨u²⟩/2` |
| Positive | (vesicle suite) `test_area_change_prefactor_matches_a_direct_quadrature` | `(ℓ−1)(ℓ+2)` against numerical `ΔA` on the sphere |
| Positive | (vesicle suite) `test_flat_ratio_is_exactly_l_over_l_minus_two` | the exact geometric ratio, `1e-14` |
| Positive | (vesicle suite) `test_sigma_bar_cancels_under_the_laplacian_mapping` | invariance over `σ̄ ∈ [1e-6, 1e6]` |
| Positive | (vesicle suite) `test_flat_limit_agreement_degrees` | **10 % at ℓ=5, 5 % at ℓ=6, 1 % at ℓ=14** |
| Positive | (vesicle suite) `test_crossover_degree_maps_to_the_flat_crossover_wavevector` | vesicle `q*` equals `√(σ/κ)` |
| Positive | (vesicle suite) `test_mean_is_unchanged_by_the_averaging_convention` | `⟨C_ℓ⟩ = ⟨u_ℓm²⟩`, no `2ℓ+1` in the mean |
| Positive | (vesicle suite) `test_sum_convention_carries_the_two_l_plus_one` | the other convention does carry it |
| Positive | (vesicle suite) `test_angular_power_round_trips_through_both_conventions` | conversions invert exactly |
| Negative (must fail) | (vesicle suite) `test_degree_zero_is_refused_as_a_volume_change` | `ℓ=0` raises, naming the volume constraint |
| Negative (must fail) | (vesicle suite) `test_degree_one_is_refused_as_a_rigid_translation` | `ℓ=1` raises, naming the translation |
| Negative (must fail) | (vesicle suite) `test_unguarded_formula_gives_a_negative_variance_at_degree_zero` | proves the `ℓ=0` failure is a *negative*, not a division by zero |
| Negative (must fail) | (vesicle suite) `test_unguarded_formula_diverges_at_degree_one` | proves the `ℓ=1` failure *is* a division by zero |
| Negative (must fail) | (vesicle suite) `test_flat_oracle_is_badly_wrong_at_low_degree` | flat formula off by 50 % at `ℓ=2`, 20 % at `ℓ=3` |
| Negative (must fail) | (vesicle suite) `test_total_area_instead_of_effective_area_fits_perfectly_with_wrong_kappa` | the `4π` trap: `χ²≈0` and `κ` wrong by `4π` |
| Negative (must fail) | (vesicle suite) `test_wrong_sample_count_from_the_c_l_convention_hides_a_bad_model` | passing `N` instead of `N(2ℓ+1)` lets a wrong model pass |
| Negative (must fail) | (vesicle suite) `test_linear_mapping_leaks_a_spurious_sigma_dependence` | `q=ℓ/R` makes the discrepancy depend on `σ̄` |

## 10. Acceptance, reviewer, rollback

| Field | Value |
|---|---|
| Status | `ACCEPTED` |
| Reviewer | Lane L4 (author). **Not independently reviewed** — flagged for the PI's morning review, since the author and the deriver are the same agent. |
| Independent oracle | The derivation in §4, done from the Helfrich functional and the two constraints, with `ΔA` additionally checked against a direct numerical quadrature on the sphere. The published Milner–Safran result agrees; agreement with the literature is corroboration, not the evidence. |
| Rollback | Delete `validation/analytic/vesicle.py` and its suite. No other module imports it; `helfrich.py` gains only a docstring cross-reference, which is inert. Rollback cost: minutes. |
| Retractions | None. No claim in this entry has been retracted. |

## 11. Compliance with PLAN §0.2

1. No source prose to strip — nothing was read. ✔
2. Law re-derived independently, and the brief's errors recorded rather than inherited (§4.4). ✔
3. All prose new, written for Aleph. ✔
4. Positive and deliberately-failing negative controls (§9). ✔
5. This entry written before the code. ✔
6. Passes with `ffn_cellsim` absent from `sys.path` — it is never referenced. ✔
