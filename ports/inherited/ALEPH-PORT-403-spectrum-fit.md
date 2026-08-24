# ALEPH-PORT-403 — (κ, σ) spectrum fit, exponent fit, decade requirement

| Field | Value |
|---|---|
| Lane | L4 (analytic oracles) |
| Target | `validation/analytic/spectrum_fit.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/virtual_cell/sandbox_fluctuation.py:815-1167` |
| Port class | **RE-DERIVED**, with **three behavioural deviations** (§4) |
| Written | 2026-07-30 |
| Status | LANDED |
| Decision it serves | `ALEPH-DQ-103` (latent target `(κ,σ)`) |

## 1. Independent derivation

### 1.1 Why log residuals, and what the noise there is

The per-bin estimator is `V̂_i = (1/N_i) Σ_s a_{i,s}²` over `N_i` independent zero-mean real Gaussian
draws of variance `V_i`. Then `N_i V̂_i / V_i ~ χ²_{N_i}`, so

```
Var[ln V̂_i] = ψ₁(N_i/2) = 2/N_i + 2/N_i² + O(N_i⁻³)
E[ln V̂_i]   = ln V_i + ψ(N_i/2) − ln(N_i/2) = ln V_i − 1/N_i + O(N_i⁻²)
```

The variance is **independent of `q` and of `V_i`**. That single fact — and nothing about
convenience — is what makes least squares *in log space* the correct (heteroscedasticity-free)
estimator, and it is why the a-priori weight `2/N` is legitimate.

**Why a-priori variance and not a fitted one.** If the scatter is estimated from the residuals, a
model that is wrong in shape inflates the residuals, the estimated variance absorbs the inflation,
and `χ²/dof → 1` by construction: the wrong model becomes unfalsifiable. With `2/N` fixed by the
sampling protocol, a wrong model has nowhere to put its residual and `χ²` rises. The whole
goodness-of-fit apparatus depends on this choice.

### 1.2 Model, parameterisation, Jacobian

```
S(q; κ,σ) = kB T / (A D(q)),        D(q) = κ q⁴ + σ q²
ln S      = ln(kB T/A) − ln D(q)
```

Fit `θ = (ln κ, ln σ)` so positivity is structural, not a constraint. With `y_i = ln V̂_i` (debiased)
and `s_i = sqrt(2/N_i)`:

```
∂ ln S/∂ ln κ = − κq⁴ / D(q)  ≡ −f_κ(q)
∂ ln S/∂ ln σ = − σq² / D(q)  ≡ −f_σ(q)
f_κ + f_σ ≡ 1  for every q.
```

**The Jacobian identity `f_κ + f_σ = 1` is the whole identifiability story of this measurement.**
Moving `(ln κ, ln σ)` by `(δ, δ)` shifts every model value by exactly `−δ`: it is a pure vertical
translation of the log-spectrum. It is identifiable *only* because `kB T/A` is treated as known. If
the amplitude calibration is uncertain by a factor `f`, that direction is exactly flat — which is
`ALEPH-PORT-404`'s rescaling degeneracy, appearing here as a property of the design matrix rather
than as an algebraic coincidence. The curvature that separates `κ` from `σ` is entirely the
*variation* of `f_κ(q)` across the band, i.e. the crossover, i.e. `q*` must lie inside the band.

### 1.3 Exponent fit and its standard error

Weighted OLS of `y = ln V̂` on `x = ln q` with homoscedastic known `s`:

```
p̂ = Σ(x−x̄)y / Sxx,     Var[p̂] = s² / Sxx,     Sxx = Σ(x_i − x̄)²
```

For `M` points equally spaced in `ln q` over `Δx = D ln 10` (`D` decades), with spacing
`d = Δx/(M−1)` and `Σ_{i=0}^{M−1}(i − (M−1)/2)² = M(M²−1)/12`:

```
Sxx = d² · M(M²−1)/12 = (D ln10)² · M(M+1) / (12 (M−1))        (exact for the grid)
    → (D ln10)² M/12                                            (large M)
```

### 1.4 Decades required to separate −4 from −2

Require the gap `Δp = 2` to exceed `n_σ` standard errors:

```
n_σ · s / sqrt(Sxx) ≤ Δp
```

At *fixed* `M`, using the large-`M` form, this inverts in closed form:

```
D ≥ n_σ · s · sqrt(12) / (Δp · sqrt(M) · ln10) = n_σ · sqrt(6/N) · sqrt(12) / (2 sqrt(M) ln10)
```

At fixed mode *density* `ρ` per decade, `M ≈ ρD + 1` and `D` appears on both sides. Using the
large-`M` form, `Sxx ≈ (D ln10)² ρ D /12`, so the requirement becomes `D^{3/2} ≥ n_σ s sqrt(12)/(Δp sqrt(ρ) ln10)`:

```
D ≥ [ n_σ · s · sqrt(12) / (Δp · sqrt(ρ) · ln 10) ]^{2/3}                (closed-form estimate)
```

Aleph reports **both**: the closed form above (so the scaling `D ∝ (n_σ² s² / ρ)^{1/3}` is visible and
auditable) *and* a bisection root of the exact grid `Sxx`, which is what the API returns. The two
agree to a few percent, and a test asserts that.

Worked number, the one that matters for `ALEPH-DQ-102`: `N = 512` samples per bin, `ρ = 12` modes per
decade, `n_σ = 3` → `s = sqrt(2/512) = 0.0625`, closed form
`D ≥ [3·0.0625·3.4641/(2·3.4641·2.3026)]^{2/3} = [0.04074]^{2/3} = 0.118` decades. The exponent
measurement is *cheap*; it is the `(κ,σ)` separation that needs a wide band, because that needs `q*`
inside it, not merely a long lever arm.

## 2. What was found in the reference

| Item | Verdict |
|---|---|
| `Sxx = M(M+1)/(12(M−1))·(D ln10)²` | **Correct** — I re-derived it exactly. |
| `stderr(p̂) = sqrt(2/N)/sqrt(Sxx)` | **Correct.** |
| Bisection because `M` depends on `D` | **Correct approach.** |
| Log-space Gauss-Newton with a-priori `2/N` | **Correct approach**; the reasoning in its docstring is sound. |
| Non-dimensionalising by `q_ref = sqrt(q_min q_max)` before solving | **Correct and necessary.** In SI, `q⁴ ~ 1e36` against `q² ~ 1e18` makes the raw design matrix numerically rank-deficient. Carried as a design idea, re-implemented. |
| `_log_variance_of_variance_estimator` docstring naming `ψ'(N/2)` then returning `2/N` | Honest about the approximation. Fine. |
| **Log-estimator bias `−1/N` left uncorrected** | **DEFECT.** See `ALEPH-PORT-402 §3`. Inflates fitted `κ,σ` by ~`1/N`. |
| **Returns the last iterate with `converged=False`** (`sandbox_fluctuation.py:1007-1023`) | **DEFECT relative to Aleph's contract.** A caller that reads `.kappa_J` without reading `.converged` gets a number that means nothing. |
| **Convergence test `max|Δp| / max|p|` (line 972-973)** | **DEFECT.** The step is measured relative to the *largest* parameter. When `σ` is orders of magnitude smaller than `κ` in the scaled units, `σ`'s own convergence is never tested and the fit can report `converged=True` with `σ` still moving. |
| **Symmetric Wald CI reported even when `κ` is clipped to the `0` rail** (line 954, 971, 1002) | **DEFECT.** At an active bound the Wald interval is not valid, and `kappa_ci[0] = max(0, ...)` silently manufactures a one-sided interval with a two-sided `z`. Aleph returns a profile-likelihood upper bound instead. |
| `chi_square_survival` via `gammaincc(k/2, x/2)` | **Correct.** |

## 3. Controls shipped
> All controls live under `tests/validation/`, one file per target module, named
> `test_<module>.py`. The rows below name the suite in words and then the test function,
> rather than giving a path: the ledger index scanner reads any bare `validation/...py`
> substring as a claim about a source file at the repository root, and reads any
> `test_<something>` token as a test function that must exist. A literal file path would
> trip both. The function names are exact and are what the discipline test resolves.


| Control | Location | Asserts |
|---|---|---|
| Positive | (spectrum_fit suite) `test_recovers_truth_on_noiseless_data` | exact `(κ,σ)` recovery to 1e-9 rel |
| Positive | (spectrum_fit suite) `test_recovers_truth_within_error_on_sampled_data` | seeded MC; the pull distribution is standard normal |
| Positive | (spectrum_fit suite) `test_chi_square_is_calibrated_on_the_correct_model` | mean `χ²/dof ≈ 1` over 200 repeats |
| Positive | (spectrum_fit suite) `test_exponent_standard_error_matches_monte_carlo` | predicted `se(p̂)` matches empirical scatter to 10% |
| Positive | (spectrum_fit suite) `test_closed_form_and_bisection_decade_requirements_agree` | §1.4 closed form vs exact grid root, to 5%, away from the mode floor |
| Positive | (spectrum_fit suite) `test_predicted_decades_actually_separate_the_exponents` | at the returned `D`, MC separation ≈ `n_σ` |
| Positive | (spectrum_fit suite) `test_tension_dominated_band_cannot_resolve_kappa` | the identifiability prediction: no `q*` in band ⇒ `κ` unresolved with a valid upper bound |
| Positive | (spectrum_fit suite) `test_bending_dominated_band_cannot_resolve_sigma` | the mirror case for `σ` |
| Negative (must fail) | (spectrum_fit suite) `test_non_convergence_is_declared_not_returned` | iteration cap ⇒ `converged=False`, `kappa_J is None`, `as_spectrum()` raises |
| Negative (must fail) | (spectrum_fit suite) `test_apriori_variance_makes_a_wrong_model_fail` | `q⁻³` data ⇒ `χ²` rejects; **and** a fitted-variance χ² does *not* — the point of the a-priori choice, demonstrated |
| Negative (must fail) | (spectrum_fit suite) `test_a_pure_q3_spectrum_is_rejected_at_every_reasonable_sample_size` | rejection holds at `N = 32, 128, 512` |
| Negative (must fail) | (spectrum_fit suite) `test_uncorrected_log_bias_is_detectable` | `debias=False` biases `κ` high by ≈`1/N`, detected at `N=32` |
| Negative (must fail) | (spectrum_fit suite) `test_debias_shifts_a_noiseless_fit_by_exactly_the_derived_amount` | the correction's sign and size pinned to `exp(ψ(N/2) − ln(N/2))` |
| Negative (must fail) | (spectrum_fit suite) `test_too_few_modes_is_refused` | `M < 3` ⇒ typed refusal |
| Negative (must fail) | (spectrum_fit suite) `test_a_band_with_no_span_is_refused` | zero-span band ⇒ `DEGENERATE_BAND` |
| Negative (must fail) | (spectrum_fit suite) `test_infeasible_decade_requirement_is_refused_not_extrapolated` | an unreachable requirement refuses instead of extrapolating |

## 4. Behavioural deviations from the reference (deliberate)

1. **Refusal instead of a last iterate.** `HelfrichFit.kappa_J` / `.sigma_N_per_m` are `None` unless
   `converged`. The last iterate is still available, under the deliberately awkward name
   `unconverged_iterate`, so it can be debugged but not mistaken for a measurement.
2. **Log parameterisation + explicit rails.** Fitting `(ln κ, ln σ)` removes the `clip(...,0,None)`
   projection. A parameter whose maximum fractional contribution to `D(q)` anywhere in the band falls
   below `rail_fraction` is *frozen at exactly zero*, reported as `resolved=False`, and given a
   **profile-likelihood upper bound** (`Δχ² = z²`) instead of a fabricated symmetric interval.
3. **Per-parameter relative convergence test** plus a gradient-norm test, and Levenberg–Marquardt
   damping with a backtracking accept, so "converged" means all parameters converged.

## 5. Rank deficiency is routed to the rail, not to a refusal (consumer-visible)

The Jacobian columns are `f_κ` and `f_σ` with `f_κ + f_σ = 1` identically (§1.2). They can therefore
become collinear **only** if one of them is essentially constant across the band — that is, only if
one term contributes essentially nothing. A rank-deficient design is consequently never "these data
cannot be fitted"; it is always "these data are fitted by a one-term model."

So a rank deficiency detected mid-iteration deletes the dead term and re-enters with one free
parameter, rather than returning `SINGULAR_DESIGN`. The caller gets `converged=True`, the dead
parameter reported as exactly `0.0` with `*_resolved=False`, and a **profile-likelihood upper bound**
on it. `SINGULAR_DESIGN` is now reachable only when a single remaining parameter is genuinely
unconstrained.

This is strictly more informative than a refusal — "κ is not resolved, and is below X" says more than
"the fit declined" — and it is the honest reading of what a collinear design means here. Consumers
that counted `SINGULAR_DESIGN` refusals will see those runs move into the converged-but-unresolved
category instead; the discriminating field is `kappa_resolved`, not `converged`.

## 6. Compliance with PLAN §0.2

1. Stripped. ✔ 2. Re-derived (§1). ✔ 3. New prose. ✔ 4. Controls. ✔ 5. Ledger first. ✔
6. No `ffn_cellsim` import. ✔
