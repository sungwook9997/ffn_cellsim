# ALEPH-PORT-501 — stationarity verdict and the error bar on a correlated mean

| Field | Value |
|---|---|
| Lane | L5 (observation layer) |
| Target | `aleph/observe/stationarity.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/ac/engine/observe/stationarity.py` (472 lines) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` |
| Source file digest | `sha256:10dccab8f7b390471ddb82b27c4ea74fc45dea543b05a1b17b860bd5fbd521d2` |
| Source test status | Reference has live tests; not re-run — the statistics were re-derived and re-verified numerically instead. |
| Port class | **RE-DERIVED**, with the `TOO_SHORT` criterion changed from an externally supplied driving time to the series' own measured correlation time |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. What the reference does

Four things: an FFT autocorrelation with a Sokal automatic window for `tau_int`; a Chodera-style
equilibration search that maximises the effective sample size of the retained tail; a signed
least-squares drift test; and `sem = sigma * sqrt(2 * tau_int / T)`. It returns a `TOO_SHORT` verdict
and withholds the mean when the series is too short.

The audit rated it the best-written module in that repository, and that assessment holds: the
statistics are correct and the module knows why it exists. This port is therefore about **what
changes**, not about whether the reference is wrong.

## 2. Independent re-derivation

### 2.1 The error on a correlated mean

For a stationary series `x_1..x_N` with sample spacing `dt`, variance `sigma^2`, and normalised
autocorrelation `rho_k`:

    Var(mean) = (sigma^2 / N) * [ 1 + 2 * sum_{k=1..N-1} (1 - k/N) * rho_k ]

The bracket is the *variance inflation factor*. As `N -> inf` with summable `rho`, it tends to
`1 + 2*sum_k rho_k`, and defining the integrated autocorrelation time in samples as

    tau_samples = 1/2 + sum_{k>=1} rho_k

the bracket is exactly `2 * tau_samples`. Hence

    Var(mean) = sigma^2 * 2 * tau_samples / N
    sem       = sigma * sqrt(2 * tau_samples / N) = sigma * sqrt(2 * tau_int / T)

with `tau_int = tau_samples * dt` and `T = N * dt`. The naive `sigma / sqrt(N)` is the special case
`tau_samples = 1/2`, i.e. an uncorrelated series — which is why `1/2` is the right floor for
`tau_samples` and not `0`.

The ratio of the two forms is `sqrt(2 * tau_samples)`. **This is not a small correction.** For AR(1)
with coefficient `phi`, `rho_k = phi^k`, so `tau_samples = (1 + phi) / (2 * (1 - phi))` and

    sem_correct / sem_naive = sqrt( (1 + phi) / (1 - phi) )

which at `phi = 0.95` is `6.24`. A published error bar six times too small is not a rounding issue;
it is a false claim, and it is false in the direction that makes a result look significant.

### 2.2 The effective sample size

`n_eff = T / (2 * tau_int) = N / (2 * tau_samples)`. Equivalently, `n_eff` is the `N` an uncorrelated
series would need to reach the same `sem`. `n_eff <= N` always, with equality iff uncorrelated.

### 2.3 Why the sum must be windowed

The estimator `rho_hat_k` has noise of order `1/sqrt(N)` per lag and the true `rho_k` decays. Summing
all `N-1` lags therefore adds a random walk of `N` steps each of size `1/sqrt(N)` — an `O(1)`
contribution with no sign preference, which in practice inflates `tau` because the sum is truncated
at whichever lag happened to be running high. The Sokal automatic window takes the smallest `M` with
`M >= c * tau(M)`: the sum is cut a fixed number of correlation times out, so the discarded tail is
exponentially small while the accumulated noise is bounded. `c = 5` is used here rather than the
usual 6 because Aleph's series are short and a larger `c` frequently fails to close the window at
all; when it does not close, that is **reported**, and `tau_int` is labelled a lower bound rather
than quietly returned as if converged.

### 2.4 What changed from the reference — the `TOO_SHORT` criterion

The reference asks whether the series covers `min_windows` multiples of an *externally supplied*
`driving_correlation_time_s` — for it, the myosin bound-head lifetime. That coupling is
domain-specific and it has a failure mode: the driving time is a property of the forcing, and the
collective relaxation of the observable can be orders of magnitude slower. A series can clear the
driving-time bar and still be far too short for its own correlation time.

Aleph's criterion is **self-referential**: a series is `TOO_SHORT` when

    T < min_tau_windows * tau_int          [default min_tau_windows = 20]

i.e. when it does not contain enough independent samples of *the thing being measured*. This needs
no external physics, so `aleph/observe/stationarity.py` is domain-independent numerical-experiment
hygiene and nothing in it knows what a membrane is. It also composes: `n_eff = T / (2 * tau_int)`, so
the rule is exactly `n_eff < min_tau_windows / 2`, i.e. at least 10 independent samples by default.

A caller who *does* have an external time scale can still impose it, via the optional
`min_span_s` argument, and both criteria are recorded in the returned contract. Neither has a
default that was chosen after looking at data.

### 2.5 The equilibration search, and its honest limitation

The discard point is the truncation whose retained tail maximises `n_eff`. Discard too little and
the transient inflates both `sigma` and `tau_int`; discard too much and `T` shrinks. The maximum
trades these off by a rule rather than by eye.

The reference documents a real failure of this rule, found on a live run: a series that rises and
then flattens can score its maximum `n_eff` at `start = 0`, because `T` is maximal there and the ACF
estimator does not penalise the trend enough to compensate — and the retained transient *inflates*
`sigma`, which makes the drift test easier to pass. The verdict is then self-reinforcing.

That failure is real and Aleph reproduces the guard against it. The re-derivation changed the
statistic, for reasons found by measurement rather than by reading. Aleph compares the mean of the
retained window's leading tenth against the mean of its trailing tenth, normalised by the **standard
error of that difference under the correlated null**:

    sd(lead - trail) = sigma_trail * sqrt(4 * tau_trail / m)

with both `sigma_trail` and `tau_trail` measured on the trailing tenth alone. Two requirements have
to hold at once, and each is violated by an obvious construction:

- **Correlation-aware.** The reference normalises by the raw *scatter* of the trailing tenth. On a
  correlated series the scatter of a block says almost nothing about the uncertainty in that block's
  mean, so the test becomes wildly over-sensitive: a strongly correlated stationary wander gets
  flagged as a transient. Since the right verdict for such a series is `TOO_SHORT`, that converts an
  honest "we could not tell" into a confident and wrong "it is drifting". Measured with the
  raw-scatter denominator: AR(1) at `phi = 0.99` over 300 samples was called `DRIFTING` in 37 of 60
  realisations.
- **Transient-proof.** Estimating `sigma` and `tau` on the *whole* retained window restores the
  self-reinforcing loop in a new place — a transient inflates `tau_int` enormously, since a monotone
  segment is correlated at every lag, so a *larger* transient becomes *harder* to detect. Measured:
  on a 300-sample rise followed by 12,000 flat samples the whole-window version failed to fire at
  all; the trailing-tenth version fires at 580 sigma.

The guard **declines** rather than guessing when the trailing tenth spans fewer than 10 correlation
times, because below that its own `tau` estimate is unusable. Such series are `TOO_SHORT`, and the
length gate is the right reporter for them.

### 2.6 The drift test needed the same correction, and one more

The reference's drift criterion is `|drift across the window| < drift_sigma * sigma`. That ratio has
**no null distribution**: it says how large a trend is relative to the scatter, which is a
description rather than a test. Worse, it is bounded — a *perfect* ramp scores at most `sqrt(12)`
≈ 3.46, because the ramp itself sets the `sigma` it is divided by. Measured: a clean ramp scores
3.44 while stationary AR(1) wanders reach 2.39, so no threshold on that ratio separates them.

Aleph tests a significance instead:

    drift_sigma = |b * T| / ( s_residual * sqrt(12 / n_eff_residual) )

from `Var(b) = s^2 / (N * Var(t))` with `Var(t) = T^2/12` for evenly spaced samples. Two choices
carry the weight:

- **The scale is the residual scatter about the fitted line, not the series scatter.** Using the
  series' own `sigma` and `tau_int` conflates the trend being tested with the noise it is tested
  against; for a strong trend both are dominated by the trend, and the statistic goes blind to
  exactly what it is looking for. Detrending fixes this at the root — the residuals of a ramp with
  white noise *are* white noise. Measured: the same clean ramp scores **1238 sigma** by this
  statistic against a loudest-wander 5.6, a separation of 220x where the reference's ratio had none.
- **The count is `n_eff`, not `N`** — the same correction, and the same mistake avoided, as for the
  mean.

The threshold is 3 sigma rather than 1, because on a correctly normalised statistic a one-sigma
threshold rejects a third of perfectly stationary series. Measured false-`DRIFTING` rate on
genuinely stationary AR(1): 1–2%.

### 2.7 Verdict order is a claim, not an implementation detail

1. an explicit `min_span_s`, if declared;
2. the transient guard;
3. the drift test;
4. the length gate;
5. otherwise `STATIONARY`.

The drift test runs **before** the length gate deliberately. A significant slope is a positive
finding at any length, because the statistic is already error-barred against the correlated null and
therefore cannot mistake a short wander for a trend. Running the gate first hides real results: a
clean ramp has an enormous `tau_int` *precisely because* it is a ramp, so the gate calls it
`TOO_SHORT` and the caller goes away and runs longer to observe the same ramp again. This was not
hypothetical — it is what the first implementation did, and the ramp control caught it.

### 2.8 The gate's own limitation, measured and stated

The length gate tests against the *measured* `tau_int`, and on a short series that measurement is
**biased low**. Measured on AR(1) at `phi = 0.99` (true `tau_samples` = 99.5):

| N | `tau` measured | ratio to truth |
|---|---|---|
| 1,000 | 64 | 0.64x |
| 2,000 | 56 | 0.56x |
| 5,000 | 98 | 0.98x |
| 100,000 | 109 | 1.10x |

The bias runs toward *fewer* correlation times, hence *more* apparent independent samples, hence an
error bar that is **still too small** — the same direction as the naive form, far less extreme. So
the gate is necessary and not sufficient. Aleph responds by setting the default to 50 correlation
times (25 independent samples) rather than a handful, by reporting `TauEstimate.reliability` =
`T / tau_int` on every result, and by annotating any report below 50 that its `sem` is a lower bound.

Aleph states §2.5–§2.8 as **known limitations with guards**, not as solved problems.

## 3. Numerical verification performed before the code landed

AR(1), `phi = 0.95`, `sigma_x = 1`, `N = 4000`, 400 independent realisations:

| Quantity | Value |
|---|---|
| Asymptotic inflation `sqrt((1+phi)/(1-phi))` | 6.245 |
| Exact finite-`N` sd of the mean (closed form, §2.1) | 0.09850 |
| Empirical sd of the mean over 400 runs | 0.09656 |
| Mean naive `sigma/sqrt(N)` | 0.01569 |
| Mean `tau_int`-corrected sem | 0.09614 |
| Empirical / naive | **6.15** |
| Corrected / empirical | **0.996** |
| Measured `tau_samples` (mean) vs true 19.5 | 18.96 |
| 95% CI coverage using naive sem | **0.225** |
| 95% CI coverage using corrected sem | **0.945** |

The coverage row is the one that matters. A nominal 95% interval built the naive way contains the
truth 22.5% of the time. That is the counterexample test, and it is the reason this module exists.

## 4. Controls shipped

All under `tests/observe/test_stationarity.py`.

| Control | Test | Asserts |
|---|---|---|
| **Counterexample** | `test_naive_error_bar_understates_the_truth_on_a_correlated_series` | naive is 6.15x too small; corrected is right to 0.5% |
| **Counterexample** | `test_naive_confidence_interval_fails_to_cover` | nominal 95% covers 22.5%; corrected covers 94.5% |
| Positive | `test_white_noise_has_tau_int_of_half_a_sample` | uncorrelated limit `tau_samples -> 1/2` |
| Positive | `test_white_noise_effective_sample_size_is_the_sample_count` | `n_eff ~= N` |
| Positive | `test_white_noise_makes_the_two_error_bars_agree` | corrected equals naive when uncorrelated |
| Positive | `test_ar1_tau_int_matches_the_closed_form` | `tau_samples = (1+phi)/(2(1-phi))` at four `phi` |
| Positive | `test_ar1_autocorrelation_is_phi_to_the_k` | the ACF estimator itself |
| Positive | `test_tau_int_is_reported_in_the_time_unit_it_was_given` | `tau_int` scales with `dt`, `tau_samples` does not |
| Positive | `test_variance_inflation_and_effective_sample_size_are_consistent` | `2*tau` and `N/(2*tau)` agree |
| Positive | `test_effective_sample_size_never_exceeds_the_sample_count` | the floor at `tau_samples = 1/2` bites |
| Positive | `test_constant_series_is_handled_rather_than_dividing_by_zero` | zero-variance path |
| Positive | `test_equilibration_discards_an_injected_transient` | recovers a planted burn-in |
| Positive | `test_equilibration_on_a_clean_series_discards_little` | no spurious discard |
| Positive | `test_stationary_series_returns_a_mean_and_a_sem` | happy path, mean recovered |
| Positive | `test_both_z_scores_are_reported_so_the_difference_is_visible` | `z_naive/z_ess = sqrt(2*tau_samples)` exactly |
| Positive | `test_drift_significance_is_enormous_for_a_clean_ramp` | §2.6: a ramp scores >100 sigma |
| Positive | `test_stationary_series_are_rarely_called_drifting` | false-positive rate at the 3-sigma threshold |
| Positive | `test_optional_absolute_span_bound_is_honoured_and_recorded` | `min_span_s` |
| Positive | `test_the_contract_is_copied_into_the_report` | thresholds travel with the verdict |
| Positive | `test_tau_reliability_is_the_number_of_correlation_times_spanned` | §2.8 reporting |
| **Negative** | `test_tau_is_biased_low_on_a_short_series` | §2.8: the estimator's own limitation, measured |
| **Negative** | `test_a_window_opening_on_a_transient_is_detected` | §2.5 guard fires on a planted rise |
| **Negative** | `test_no_mean_is_ever_reported_across_a_rise` | the invariant, whichever branch is taken |
| **Negative** | `test_the_transient_guard_does_not_fire_on_a_stationary_correlated_series` | §2.5: the over-sensitivity the raw-scatter denominator would cause |
| **Negative** | `test_a_window_too_short_to_split_into_tenths_declines_rather_than_guessing` | declines instead of inventing a verdict |
| **Negative** | `test_normalising_the_drift_by_the_series_variance_hides_a_clean_ramp` | §2.6: the reference's ratio cannot separate a ramp from a wander; the significance separates them 220x |
| **Negative** | `test_series_short_relative_to_its_own_tau_is_refused` | `TOO_SHORT`, `mean is None` |
| **Negative** | `test_a_strongly_correlated_short_series_is_refused_in_practice` | 0/40 short AR(1) series certified `STATIONARY` |
| **Negative** | `test_too_short_verdict_withholds_the_sem_as_well` | no error bar either |
| **Negative** | `test_naive_sem_is_reported_even_on_a_refusal` | the correction's size is in every record |
| **Negative** | `test_drifting_series_returns_no_mean` | `DRIFTING`, `mean is None` |
| **Negative** | `test_drift_is_reported_signed` | rising and falling are different findings |
| **Negative** | `test_non_finite_sample_is_refused` | no silent NaN |
| **Negative** | `test_iteration_index_as_dt_is_rejected` | `dt` is required and unitful |
| **Negative** | `test_too_few_samples_is_rejected` | fewer than 4 samples |

## 4a. Units, domains, singular cases, invariants

**Units.** `dt` is a sample interval in seconds and is **required with no default**, because passing
an iteration index in its place silently returns a correlation time measured in iterations — the
exact confusion the module exists to prevent. Consequently `tau_int_s` carries whatever time unit
`dt` carried, while `tau_samples`, `n_effective`, `sem_inflation`, `z_ess`, `z_naive`,
`drift_sigma_measured` and `reliability` are dimensionless. `mean`, `sem`, `sem_naive` and `std`
carry the observable's own units; `drift_per_s` carries them per unit time.

`dt` must be the **sampling** interval, not the integrator step. With a recording stride the two
differ, and handing the integrator step to the estimator understates `tau_int` by exactly that
stride — an error in the flattering direction.

**Domain.** At least 4 finite samples, `dt` positive and finite. Everything else is a *verdict*
rather than a rejection, which is the point of the module.

**Singular cases.**

- Zero variance: `tau_samples` is floored at the uncorrelated value 1/2 and `zero_variance` is set,
  rather than dividing by an autocorrelation of zero at lag 0.
- Sokal window never closes: reported via `window_closed=False`, and `tau_int` is labelled a
  **lower** bound with `n_effective` an upper one, rather than being returned as if converged.
- Trailing tenth too short to carry its own `tau`: the transient guard **declines** (returns
  `(False, 0.0)`) instead of inventing a verdict; the length gate then reports `TOO_SHORT`.
- Fewer than 4 samples, non-finite samples, or a non-positive `dt`: `ValueError`. These are
  programming errors in the caller, not findings about the data, which is why they raise where a
  verdict would otherwise be returned.

**Invariants.**

- `tau_samples >= 1/2` always, so `n_effective <= n_samples` always, with equality only for an
  uncorrelated series.
- `sem / sem_naive == sem_inflation == sqrt(2 * tau_samples)` exactly — asserted to 1e-9.
- `z_naive / z_ess == sqrt(2 * tau_samples)` exactly, the same factor, asserted to 1e-9.
- `mean` and `sem` are `None` together, and both only when the verdict is `STATIONARY`.
- `sem_naive` is reported on **every** report including refusals, so the size of the correction is
  in the record rather than only in a docstring.
- `drift_per_s` is signed, always: rising toward a steady state and falling toward one are different
  physics.

## 4b. Numerical and precision envelope, and production-backend residency

**Numerical.** All float64. The autocorrelation is computed by FFT on the mean-subtracted series,
zero-padded to at least `2N-1` so the circular correlation the FFT computes equals the linear one;
cost is `O(N log N)`. The **biased** normalisation (dividing by `N`, not `N-k`) is deliberate: the
unbiased form's variance grows without bound as `k` approaches `N`, where a handful of pairs
contribute, and feeding those lags into a sum is how a correlation time acquires a long random tail.
Tolerances in the controls are statistical, not numerical — 10–15% on estimator agreement, driven by
the sampling error at the tested `N`, and stated per test.

**Residency.** Host-side, CPU only, `numpy` plus `numpy.fft`. This is deliberate and should stay
that way. The module consumes a *recorded telemetry series* — one scalar per observable per recorded
frame — which is four to six orders of magnitude smaller than the state it was reduced from, so the
transfer has already happened by the time this code runs and there is nothing to gain from a device
implementation. The transfer boundary is: the runtime reduces on device, copies scalars to host, and
this module reads the host array. No device memory is touched here.

## 5. Compliance with PLAN §0.2

1. All reference comments and docstrings discarded; nothing carried. ✔
2. Every statistic re-derived (§2) and re-verified numerically (§3) before the code landed. ✔
3. New prose written for Aleph's situation; the reference's cell-biology framing is gone and the
   module is domain-independent. ✔
4. Positive, negative, and counterexample controls (§4). ✔
5. Entry written before the code. ✔
6. No `ffn_cellsim` import; no `validation` import. ✔

## 6. Verdict on the reference

**Statistics correct; coupling too tight.** The estimators are right and the module's self-awareness
about its own failure modes is unusual and worth preserving. What Aleph changes is the dependency:
the reference's refusal criterion needs a domain constant supplied from outside, which makes a piece
of general numerical hygiene into a cell-biology module. Aleph's criterion is measured from the
series itself, and the external bound is optional rather than required.
