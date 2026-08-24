# ALEPH-PORT-405 — exact Ornstein–Uhlenbeck propagation and the Boltzmann control

| Field | Value |
|---|---|
| Lane | L4 (analytic oracles) |
| Target | `validation/analytic/ou_process.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/virtual_cell/stochastic_node.py:1870-1910, 2034-2060` |
| Port class | **RE-DERIVED**, scope enlarged (§3) |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. Independent derivation

### 1.1 The physical starting point

Overdamped Langevin for a particle in a quadratic well, friction `γ`, stiffness matrix `K` (symmetric
positive definite), temperature `T`:

```
γ dx = −K x dt + sqrt(2 γ kB T) dW
dx   = −(K/γ) x dt + sqrt(2 kB T/γ) dW
```

The noise amplitude is *not* free: it is fixed by the fluctuation–dissipation theorem, which is the
only reason a Boltzmann stationary state exists. Writing the general linear SDE

```
dx = −B x dt + L dW,        Σ ≡ L Lᵀ                (B = K/γ, Σ = (2 kB T/γ) I)
```

### 1.2 Stationary covariance — Lyapunov

Let `P = lim ⟨x xᵀ⟩`. Itô on `x xᵀ` in the stationary state gives `d⟨xxᵀ⟩/dt = 0 = −BP − PBᵀ + Σ`:

```
B P + P Bᵀ = Σ                                                  (continuous Lyapunov)
```

For the gradient case, substitute the Boltzmann guess `P = kB T K⁻¹`:

```
(K/γ)(kB T K⁻¹) + (kB T K⁻¹)(K/γ) = 2 kB T/γ · I = Σ   ✔
```

so `P = kB T K⁻¹` **is** the stationary covariance, which is exactly the covariance of
`p(x) ∝ exp(−xᵀKx / (2 kB T))`. Boltzmann is recovered, not assumed. In 1-D: `Var = kB T/k`.

### 1.3 Exact finite-time propagation

The solution of the linear SDE over `dt` is

```
x(t+dt) = e^{−B dt} x(t) + ∫₀^{dt} e^{−B(dt−s)} L dW(s)
```

so the update is exactly Gaussian with mean `e^{−B dt} x(t)` and covariance
`C(dt) = ∫₀^{dt} e^{−Bs} Σ e^{−Bᵀ s} ds`. Rather than quadrature, use the Lyapunov solution:

```
d/ds [ e^{−Bs} P e^{−Bᵀ s} ] = −e^{−Bs}(B P + P Bᵀ)e^{−Bᵀ s} = −e^{−Bs} Σ e^{−Bᵀ s}
```

Integrating `0 → dt`:

```
C(dt) = P − e^{−B dt} P e^{−Bᵀ dt}
```

This costs one `expm` and one Lyapunov solve and is **exact for every `dt`, including `dt = ∞`**
(`C → P`). Moment propagation from an arbitrary Gaussian `(m₀, C₀)`:

```
m(dt) = e^{−B dt} m₀
C(dt) = P + e^{−B dt}(C₀ − P) e^{−Bᵀ dt}
```

Stationary autocovariance at lag `τ ≥ 0`: `⟨x(t+τ)x(t)ᵀ⟩ = e^{−Bτ} P`.

### 1.4 The negative control, in closed form

Euler–Maruyama in 1-D on `dx = −θ x dt + sqrt(2D) dW`:

```
x_{n+1} = (1 − θ dt) x_n + sqrt(2 D dt) ξ_n,      ξ ~ N(0,1) i.i.d.
```

Its stationary variance is a geometric series, available in closed form:

```
Var_EM = 2 D dt / (1 − (1 − θ dt)²) = 2D dt / (2θ dt − θ² dt²) = (D/θ) · 1/(1 − θ dt/2)
```

so with `Var_exact = D/θ`:

```
Var_EM / Var_exact = 1 / (1 − θ dt / 2)      = 1 + θdt/2 + O((θdt)²)
```

**Euler–Maruyama systematically over-heats the particle, by a first-order-in-`dt` amount.** Concrete
values: `θdt = 0.01 → +0.50%`; `θdt = 0.1 → +5.3%`; `θdt = 1 → ×2 exactly`; `θdt → 2⁻` → divergent;
`θdt ≥ 2` → unstable, no stationary variance exists. `θdt = 1` giving *exactly* twice the correct
variance is a clean, memorable, exactly-assertable failure — this is the mandated
**positive control that must fail**: at `dt → 0` the integrator passes a 1% tolerance, at `θdt = 1`
it fails the same tolerance by a factor of two, and the test asserts *both*, so the test is proven
capable of failing.

Multivariate: `x_{n+1} = (I − B dt)x_n + sqrt(dt) L ξ` has stationary covariance solving the
**discrete** Lyapunov equation `P_EM = M P_EM Mᵀ + Σ dt` with `M = I − B dt`, and is stable iff the
spectral radius of `M` is `< 1`.

## 2. What was found in the reference

| Item | Verdict |
|---|---|
| `stochastic_node.py:1900` `solve_continuous_lyapunov(A, −2 D)` with `A` the (negative-definite) drift | **Correct.** SciPy solves `A X + X Aᵀ = Q`; with `A = −B` and `Q = −Σ` this is `B P + P Bᵀ = Σ`. Their `noise` is `D` in the `sqrt(2D)dW` convention, so `Σ = 2D`, and the factor of 2 is right. |
| `stochastic_node.py:1902-1905` `C = P + e^{At}(C₀ − P)e^{Aᵀt}` | **Correct** — identical to my §1.3 result. |
| `stationary_mean = solve(A, −offset)` for an affine drift | **Correct.** |
| Guard that refuses an unstable drift matrix before quoting a stationary covariance | **Correct and important.** Adopted. |
| Framing: "comparing Euler–Maruyama against this separates DISCRETISATION error from CLOSURE error" | **Correct and a genuinely good idea.** Adopted as the reason this module exists. |
| `boltzmann_density` for a 1-D gradient system | **Correct in form** (`p ∝ exp(−U/D)` for `dX = −U′dt + sqrt(2D)dW`). Note this is only closed-form-normalisable for the harmonic case; for a general `U` the normalisation is a quadrature and must be labelled semi-analytic, which Aleph does. |
| Closed-form Euler–Maruyama stationary variance `1/(1 − θdt/2)` | **NOT PRESENT.** The reference compares EM to exact numerically. Having the error in closed form turns the negative control from "it drifted" into "it is wrong by exactly ×2, assert it". Aleph adds it. |

**No error found in the reference's OU algebra.** This module is the one where the reference is
cleanly right; the additions are in the controls, not in corrections.

## 3. Controls shipped
> All controls live under `tests/validation/`, one file per target module, named
> `test_<module>.py`. The rows below name the suite in words and then the test function,
> rather than giving a path: the ledger index scanner reads any bare `validation/...py`
> substring as a claim about a source file at the repository root, and reads any
> `test_<something>` token as a test function that must exist. A literal file path would
> trip both. The function names are exact and are what the discipline test resolves.


| Control | Location | Asserts |
|---|---|---|
| Positive | (ou_process suite) `test_lyapunov_solution_satisfies_its_own_equation` | residual `‖BP+PBᵀ−Σ‖ < 1e-12` rel |
| Positive | (ou_process suite) `test_gradient_system_stationary_covariance_is_boltzmann` | `P = kB T K⁻¹` to 1e-12, derived two independent ways |
| Positive | (ou_process suite) `test_step_covariance_matches_brute_force_quadrature` | `C(dt)` vs a 20 001-point Simpson integral |
| Positive | (ou_process suite) `test_exact_propagator_is_stationary_at_any_dt` | exact at `θdt` from 1e-6 to 1e6 |
| Positive | (ou_process suite) `test_boltzmann_density_is_correctly_normalised` | the closed-form normalisation integrates to 1 |
| Positive | (ou_process suite) `test_exact_sampler_reproduces_stationary_variance` | seeded MC within an effective-sample-size error bar |
| Positive | (ou_process suite) `test_autocovariance_decays_as_the_propagator` | `⟨x(τ)x(0)⟩ = e^{−Bτ}P` |
| **Negative (must fail)** | (ou_process suite) `test_euler_maruyama_fails_at_large_dt` | at `θdt = 1`, EM variance is `2×` truth and **fails** the 1% gate the same code passes at `θdt = 1e-4` |
| Negative (must fail) | (ou_process suite) `test_euler_maruyama_closed_form_matches_simulation` | the `1/(1−θdt/2)` law confirmed by seeded MC |
| Negative (must fail) | (ou_process suite) `test_euler_maruyama_converges_to_the_exact_answer_as_dt_goes_to_zero` | first order in `dt` with the derived coefficient `1/2` |
| Negative (must fail) | (ou_process suite) `test_euler_maruyama_is_unstable_beyond_theta_dt_two` | `None`/`inf`, not a silent large number |
| Negative (must fail) | (ou_process suite) `test_unstable_drift_is_refused` | a non-decaying eigenvalue ⇒ refusal |
| Negative (must fail) | (ou_process suite) `test_non_confining_well_is_refused` | an indefinite stiffness has no normalisable Boltzmann state |

## 4. Compliance with PLAN §0.2

1. Stripped. ✔ 2. Re-derived (§1). ✔ 3. New prose. ✔ 4. Controls. ✔ 5. Ledger first. ✔
6. No `ffn_cellsim` import. ✔
