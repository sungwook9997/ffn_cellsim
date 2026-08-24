# stage 0.5 — the timescale separation, measured

**Run:** `ac_gate_b_cortex_motor_native.py`, build `4095882c`, A5000 `cuda:0`, 2026-07-28 (night).
**Population:** 70,686 cortical filaments / 494,802 actin nodes / 674,314 total nodes / 8,840 NMII heads —
**FULL NATIVE**. 1,500 accepted steps at `dt_phys` 0.01 = **15.0 s of physical time**, sampled every 5th
step, `--catch-slip`, `k_xb` = 1000 pN/µm, membrane subdiv 7.

## What this run was for

Nothing in the repository had ever measured how fast the cortex's tension decorrelates, so the
stationarity contract — "the analysed window must be at least N driving correlation times" — was being
argued over a denominator nobody had measured. `observe/stationarity.py`'s own docstring assumes the
collective relaxation is *slower* than the process driving it. This run tests that.

## Result — all four observables STATIONARY

| Observable | Verdict | Mean | sem | τ_int [s] | n_eff |
|---|---|---|---|---|---|
| γ_total | STATIONARY | 4.3256 | 0.00047 | 0.0271 | 269 |
| γ_network | STATIONARY | 4.3030 | 0.00046 | 0.0268 | 272 |
| γ_source (crossbridge) | STATIONARY | 0.02259 | 1.4e-05 | 0.0250 | 300 |
| bound fraction | STATIONARY | 0.98261 | 0.00032 | 0.0295 | 255 |

Contract, declared before the run: window ≥ **5 × driving τ**, `|drift| < 1σ`. Driving τ derived from the
run's own kinetics rather than assumed — Pereverzev at zero load is `k_catch0 + k_slip0` = 0.70/s, and
the mean over the bound heads' actual loads gives **1.456 s**, the larger and therefore the one used.
Transient discarded: **0.45 s**, chosen by maximising n_eff. Drift over the analysed window: **0.24σ**.

## The finding — the module's own assumption is refuted

**γ decorrelates at least 54× FASTER than the process driving it.** With 8,840 heads, individual binding
and unbinding averages out almost immediately at whole-cell level.

**τ_int is protocol-floored and must be read as an upper bound.** The estimator's floor is
`sample_dt/2` = 0.025 s and the measurements sit at 0.025–0.030 — the autocorrelation has already fallen
to ~0.09 by one sample interval. What was measured is *"uncorrelated at 0.05 s sampling"*, not
*"τ_int = 0.027 s"*. Resolving the true value needs per-step sampling. Either way the conclusion holds:
the observable's correlation time is far below the driving time, so the driving-time floor is the binding
constraint on the contract and is ~2.4× more conservative than the statistics require.

## What may NOT be quoted from this run

- **The magnitudes.** Every NMII parameter here is a PI-GAP (`k_xb`, `f_stall`, `kappa`, `N_side`,
  `L_bb`, `k_on`); `QuantitativeClaim` is BLOCKED and the record says so.
- **γ as a converged cortical tension.** Every step is **force-accepted** (`accepted_d = ones`); the
  inner solve ran a fixed 40 iterations and was never asked whether it had converged. What settled is the
  trajectory of a force-accepted integration, which is a statement about the observable's time series and
  not about mechanical balance. `timing.comparable` is correspondingly `false`.
- **A stationarity PASS as a gate result.** The amendment that would make settling an acceptance
  criterion is PROPOSED and unratified; these verdicts are recorded, not scored.

## Figures

- `figs/stationarity_trajectory.png` — γ_total with the discarded transient shaded, the steady-state mean
  and its ±sem band; the verdict and the contract that produced it printed together.
- `figs/stationarity_source_network_split.png` — source and network on separate axes (they differ by more
  than two decades) plus the closure residual, so a split that does not sum would be visible.
- `figs/stationarity_autocorrelation.png` — the ACF against the estimator floor, which is where a reader
  sees that τ_int is set by the sampling protocol rather than by the physics.

Regenerate: `python ffn_sim/scripts/ac_stationarity_vis.py --record <this dir>/record.json`
