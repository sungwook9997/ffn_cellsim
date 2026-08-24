# τ grows with the window — and my τ_relax estimate was the FASTEST mode, not the slowest

**Status:** MEASUREMENT + a correction to this session's own arithmetic. Native, `cuda:0`,
4,361,496 nodes, cortex bound, thermostat at mobility 1e-6. 1,600 steps, 109 minutes.
Artifact `world_phase4/gthermo_eq.json`.

⚠ **No γ magnitude is claimed.** The mobility is unsourced, the coefficients are test points, and the
γ values here reach 252 pN/µm — orders above any physiological band. Nothing in this run is quotable
as a cortical tension. What is measured is a property of the SERIES.

---

## 1. The question this run was sized to answer

`THERMOSTAT_RUN_LENGTH_IS_DERIVED_2026-08-21.md` §5.1: `τ_relax` computes to 227 steps but the
measured `tau_samples` was 35, a factor of six apart, and the next run's cost depends on which is
right. Two readings:

* **(a)** the series has not reached its slow modes, so τ will GROW.
* **(b)** γ's correlation time is genuinely shorter than a single node's, because γ is a global sum
  over ~8e6 elements averaged over 64 orientations and individual slow modes partly cancel.

One measurement decides it: **does `tau_int` stop growing as the series lengthens?**

## 2. It is (a), and unambiguously

| window | `tau_int` | in samples | `n_eff` | `tau_reliability` |
|---:|---:|---:|---:|---:|
| 300 | 1.139 s | 22.8 | 4.28 | 8.56 |
| 800 | 2.575 s | 51.5 | 4.17 | 8.35 |
| **1,600** | **9.384 s** | **187.7** | **4.16** | **8.31** |

**τ grew 8.2× while the window grew 5.3×, and `n_eff` did not move.** That is the signature, and it
is sharper than "τ is still growing":

> `n_eff` pinned at ~4.16 across a 5× change in window length means **τ ≈ T/4 at every length.** The
> autocorrelation is not decaying at all — it is tracking the window, because the **trend** dominates
> it. **You cannot measure a correlation time on a drifting series**; the number you get back is a
> restatement of how long you looked.

γ over the run confirms it: 0.189 → 251, still climbing at the end, with the up/down ratio moving
toward balance (136/63 → 109/90) but not reaching it.

## 3. ⚠ CORRECTION — `τ_relax = 1/(µk)` is the fastest mode, and I sized the run with it

The 1,600-step length came from `τ_relax = 1/(mobility · k) = 227 steps`, with `k = 88,000 pN/µm` the
**axial stiffness of one bond**. That is the relaxation time of a SINGLE SPRING against the medium.

**A network's slowest mode is not its stiffest bond.** The collective modes — the ones that carry a
cortex from a relaxed build to a thermal steady state — relax far more slowly than any individual
bond, and they are what equilibration waits on. **Seven "τ" by the single-bond measure is not seven
relaxation times of the thing that has to equilibrate.**

That is my error, and it is a specific one: I took a per-element constant as a system-level timescale
because the arithmetic was clean and the units were right. The same shape as the near-misses beside
`d0` — dimensionally correct, well-founded, and about the wrong object.

## 4. What this leaves, honestly

**Not measurable from this run:** the length a stationary run needs. §5.2's table priced the options
from `τ = 35` and `τ = 227` and **both were readings of a series that has no correlation time yet**.
The table's structure stands — the length bar is still a factor of ten in GPU hours — but the τ column
it multiplies is not yet a number.

**What would give it:** equilibration first, judged by γ's own trend flattening rather than by a
timescale computed from a bond. `min_tau_reliability` is the field that says when a τ estimate is
trustworthy — `observe/stationarity.py`'s own docstring puts that threshold at 50 — and it sat at
**8.3 for all three windows**, so every τ above is biased low and every `n_eff` is an upper bound.

⚠ **And a cheap way to get there is available and should be refused for the same reason as before.**
Raising the mobility shortens every timescale in the problem. It would produce an equilibrated series
in minutes, and the number that made it possible would be an unsourced constant in a driver flag.
**The run is expensive because the physics is slow, not because the code is.**

## 5. What the run DID establish

* **The thermostat produces a genuinely fluctuating series** — 136 up / 63 down in the first eighth,
  109 / 90 in the last. The frozen arm is degenerate and the relaxed arm is monotone; this is neither.
* **The within-eighth spread falls as it climbs** — 13.1 → 4.2 — which is what approaching a steady
  state looks like, even though it has not arrived.
* **The drift test is not fooled by any of it.** Every window returns `DRIFTING`, correctly, and the
  degeneracy check added earlier the same night would catch the opposite failure. **Both ends of the
  gate work; what is missing is the run, not the criterion.**
