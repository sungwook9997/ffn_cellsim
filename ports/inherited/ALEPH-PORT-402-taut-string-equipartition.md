# ALEPH-PORT-402 — 1-D taut string, equipartition as an identity

| Field | Value |
|---|---|
| Lane | L4 (analytic oracles) |
| Target | `validation/analytic/taut_string.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/virtual_cell/sandbox_fluctuation.py:594-672, 1490-1564` |
| Port class | **RE-DERIVED**, and **scope changed** (see §4) |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. Independent derivation

String of rest length `L`, held at tension `τ` (units of force, `[τ] = M T⁻²·L = F`), transverse
displacement `h(x)`, fixed ends `h(0) = h(L) = 0`, temperature `T`. Quadratic (small-slope) energy:

```
E[h] = (τ/2) ∫₀^L (∂h/∂x)² dx
```

Fixed ends make the Dirichlet sine basis complete and *exact* — not an approximation, not a truncation:

```
h(x) = Σ_{n≥1} a_n sin(q_n x),     q_n = nπ/L
```

Using `∫₀^L cos(q_n x) cos(q_m x) dx = (L/2) δ_{nm}` (for `n,m ≥ 1`):

```
E = (τ/2) Σ_{n,m} a_n a_m q_n q_m ∫₀^L cos(q_n x)cos(q_m x) dx
  = (τ L /4) Σ_n q_n² a_n²
  = Σ_n (1/2) k_n a_n²,        k_n ≡ τ L q_n² / 2
```

The energy is **exactly diagonal**. Each `a_n` is one real quadratic degree of freedom, so the
Boltzmann distribution factorises exactly and

```
⟨E_n⟩ = kB T / 2          (an identity — holds for every n, every L, every τ, at every temperature)
⟨a_n²⟩ = kB T / k_n = 2 kB T / (τ L q_n²)
```

There is no limit, no continuum approximation, no high-`n` asymptotics anywhere in that chain. This
is precisely why it is the estimator's calibration case: **any fitting routine that cannot return
`⟨E_n⟩/(kB T/2) = 1` here is broken as software, independent of any physics question.**

**Exact sampling distribution of the verdict.** With `N` independent configurations and `M` modes,
`2E_n/(kB T)` per sample is `χ²₁`, and the `N·M` values are independent, so

```
X ≡ (2/(kB T)) Σ_{n,s} E_{n,s} ~ χ²_{N M}          exactly
```

That gives a goodness-of-fit test with no asymptotic content at all, which is what
`check_equipartition` reports. Per-mode, `N·⟨E_n⟩_meas/(kB T/2) ~ χ²_N`, so the relative standard
error on a single mode's `⟨E_n⟩` is exactly `sqrt(2/N)`.

## 2. What was found in the reference

| Item | Verdict |
|---|---|
| `E = Σ_n (τL/4) q_n² a_n²`, `k_n = τ L q_n²/2`, `⟨a_n²⟩ = 2kBT/(τ L q_n²)` | **Correct** — I reproduced all three. |
| `mean_energy_per_quadratic_dof` reporting in units of `kB T/2` rather than joules | **Correct and a good idea.** Carried as a design principle, re-implemented: it makes a factor-of-2 error visible instead of plausible. |
| `measure_mode_variance` using `mean(a²)` and not `np.var` | **Correct.** The mode mean is zero by symmetry; estimating it would burn a degree of freedom and bias the estimate low by `1/N`. |
| The reference's `check_equipartition` (line 1510) | **Different function entirely** — it tests the *Helfrich* patch for a wrong `kB T` scale given externally supplied `(κ,σ)`. It is *not* the taut-string identity check this lane was asked for. |
| Reference `check_equipartition` uses `stderr = sqrt((2/N)/M)` on a mean of logs | **DEFECT FOUND — see §3.** |

## 3. Defect found in the reference (log-estimator bias)

`V̂ = V · χ²_N / N`, so `E[ln V̂] = ln V + ψ(N/2) − ln(N/2)`, and the asymptotic expansion of the
digamma gives

```
bias(N) ≡ ψ(N/2) − ln(N/2) = −1/N − 1/(3N²) + O(N⁻⁴)
```

The reference averages `M` log-estimates and quotes a standard error `sqrt(2/(N M))` while leaving
this `−1/N` bias in. The bias-to-noise ratio is

```
|bias| / stderr = (1/N) / sqrt(2/(N M)) = sqrt(M / (2N))
```

which **grows with the number of modes**. At `N = 100, M = 50` it is `0.5σ`; at `N = 100, M = 400`
it is `1.41σ`; at `N = 50, M = 400` it is `2σ`. So the reference's `consistent` flag drifts towards
*false* on perfectly correct data as the measurement gets better in the mode direction. That is a
false-alarm mechanism that gets worse with more data — the worst kind.

The same bias propagates into the reference's log-space `(κ,σ)` fit: the fitted restoring density is
inflated by roughly `exp(1/N)`, i.e. `κ` and `σ` both high by about `1/N` relative
(1% at `N = 100`), which is a systematic that the reported statistical error bar does not cover.

**Aleph fixes this**: `validation/analytic/spectrum_fit.py` exposes `log_estimator_bias(N)` and
applies it by default (`debias=True`), and `test_spectrum_fit.py` contains a
Monte-Carlo test that *fails* with `debias=False` at the `N`/`M` where the reference would fail.

A second, smaller point: `Var[ln V̂] = ψ₁(N/2) = 2/N + 2/N² + O(N⁻³)`, so the mandated a-priori `2/N`
understates the true log-variance by a relative `1/N`. Aleph uses `2/N` as the contract default
(as instructed) and offers `log_variance_mode="exact"` for `ψ₁(N/2)`; the difference is recorded, not
hidden.

## 4. Scope change from the reference

The reference has no taut-string equipartition verdict; it only has the spectrum classes plus a
Helfrich-scale check that borrows the name. Aleph's `check_equipartition` is a new function with an
*exact* `χ²_{NM}` sampling distribution rather than a log-space normal approximation.

## 5. Controls shipped
> All controls live under `tests/validation/`, one file per target module, named
> `test_<module>.py`. The rows below name the suite in words and then the test function,
> rather than giving a path: the ledger index scanner reads any bare `validation/...py`
> substring as a claim about a source file at the repository root, and reads any
> `test_<something>` token as a test function that must exist. A literal file path would
> trip both. The function names are exact and are what the discipline test resolves.


| Control | Location | Asserts |
|---|---|---|
| Positive | (taut_string suite) `test_equipartition_is_an_identity_not_a_limit` | `⟨E_n⟩ = kBT/2` exactly, for every `n`, analytically |
| Positive | (taut_string suite) `test_energy_diagonalises_against_a_direct_real_space_quadrature` | sine-basis energy matches numerical `∫(dh/dx)²` |
| Positive | (taut_string suite) `test_sampled_amplitudes_recover_equipartition` | seeded draw passes the verdict |
| Positive | (taut_string suite) `test_chi_square_statistic_has_the_stated_distribution` | KS test of the `χ²_{NM}` statistic over 600 repeats |
| Positive | (stats suite) `test_monte_carlo_confirms_the_bias_and_the_variance` | the `-1/N` bias and `ψ₁(N/2)` variance confirmed by 200 000 draws |
| Negative (must fail) | (taut_string suite) `test_temperature_scaled_amplitudes_are_rejected` | amplitudes drawn at 1.2× and 0.8× the temperature are rejected |
| Negative (must fail) | (taut_string suite) `test_mode_dependent_error_with_a_correct_total_is_still_caught` | correct total energy, wrong distribution across modes |
| Negative (must fail) | (stats suite) `test_the_contract_two_over_n_is_demonstrably_not_exact` | `2/N` shown to be an approximation, not a fact |
| Negative (must fail) | (taut_string suite) `test_verdict_can_fail_on_a_wrong_stiffness` | a factor-2 stiffness error is caught |
| Negative (must fail) | (spectrum_fit suite) `test_uncorrected_log_bias_is_detectable` | the reference's uncorrected `-1/N` bias, measured at `N = 32` |

## 6. Compliance with PLAN §0.2

1. Stripped. ✔ 2. Re-derived (§1). ✔ 3. New prose. ✔ 4. Controls. ✔ 5. Ledger first. ✔
6. No `ffn_cellsim` import. ✔
