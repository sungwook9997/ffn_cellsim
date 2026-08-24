# The equilibration run length is DERIVED, and the tempting shortcut is a charter violation

**Status:** ARITHMETIC + a recorded refusal. No result yet; the run it sizes is in flight.

---

## 1. What the 300-step run showed

γ under the thermostat, mobility 1e-6 µm/(pN·s), cortex bound at k = 88,000 pN/µm:

| quarter | mean γ | std | up / down |
|---|---:|---:|---:|
| 1 | 18.52 | 7.49 | 51 / 23 |
| 2 | 35.76 | 4.93 | 51 / 23 |
| 3 | 50.51 | 4.11 | 47 / 27 |
| 4 | 60.38 | 5.84 | 47 / 27 |

`DRIFTING`, 14.2 σ, slope +3.26 /s. **It fluctuates and it has not settled.** Both halves matter: the
fluctuation is what the frozen and relaxed arms could not produce, and the drift is why it is not yet
a stationary series.

## 2. Why, in one line of arithmetic

```
tau_relax = 1 / (mobility * k)  =  1 / (1e-6 * 88000)  =  11.36 s  =  227 steps at dt = 0.05
300 steps = 1.32 tau
```

**The run was one relaxation time long.** Each step relaxes 0.44% of a displacement
(`mobility*k*dt = 4.4e-3`), so a network released from a relaxed configuration is still climbing
toward its thermal amplitude when the run ends. Equilibration wants roughly 5 τ, and the stationary
window sits on top of that.

Sanity, from the other direction: the equilibrium RMS displacement per degree of freedom is
`sqrt(kT/k) = 2.21e-04 µm` and the per-step kick is `2.07e-05 µm`, so the walk needs ~100 steps to
reach amplitude and several τ for the FORCE field to follow. The two estimates agree.

## 3. ⚠ The shortcut, and why it is refused

The obvious way to finish sooner is to raise the mobility, because τ is `1/(µk)`:

| mobility | τ | 5τ | dt/τ |
|---:|---:|---:|---:|
| 1e-6 | 11.36 s | **1,136 steps** | 0.004 |
| 1e-5 | 1.14 s | 114 steps | 0.044 |
| 1e-4 | 0.11 s | **11 steps** | 0.440 |

A hundred-fold shorter run is one keystroke away, and the step stays inside its stability bound.

⚠ **And mobility is precisely the UNSOURCED quantity in this run.** Turning it up so a run finishes is
choosing an unsourced constant to make a gate pass — the thing the charter names first and the thing
this engine spent 2026-08-20 cataloguing in five other disguises. **It would also be invisible in the
result**: the series would equilibrate, pass the stationarity contract honestly, and the number that
made it possible would sit in a driver flag.

**So the run length is derived from the constants and the run is 1,600 steps = 7 τ ≈ 110 min.**
The arithmetic is here so the choice is checkable, and so that a future session reading a long
runtime does not "fix" it.

## 4. What this does not claim

Nothing about γ's magnitude — the mobility is unsourced, so the TIMESCALE is unsourced, and τ above
is `1/(µk)` in whatever units µ turns out to have. **What is sourced is the RATIO**: a run must be
several relaxation times long, and that statement survives whatever the drag law turns out to be. The
number of steps changes; the requirement does not.

⚠ And if the series does equilibrate, it will be the first in this engine to pass a stationarity
contract for a true reason — the frozen arm is degenerate (`ecfd3f72`: a constant certifies as a
perfect steady state with `sem = 0.0` if run long enough) and the relaxed arm drifts. Passing is not
the same as being right, and the report's `min_n_effective` is the number to read first.


---

## 5. ⚠ CORRECTION — 7 τ buys EQUILIBRATION, not a MEASUREMENT

Session C, on reading §3. The arithmetic in §2 sized the run for the transient and stopped there;
**the transient is then DISCARDED and the contract is evaluated on what remains.**

```
1,600 steps  −  7 τ transient (1,589)  =  an 11-step analysed window
```

The contract wants `n_eff >= min_tau_windows/2 = 25`, i.e. an analysed window of `50 × tau_samples`
**on top of** the discarded transient:

| length bar | τ = 35 (measured) | τ = 227 (= τ_relax) |
|---|---|---|
| **50 τ** — `observe/stationarity`'s default | 2,885 steps = 3.3 h/seed · **9.8 h × 3** | 12,485 = 14.1 h/seed · **42.2 h × 3** |
| **5 τ** — engine's `min_windows`, PI-open | 1,310 steps = 1.5 h/seed · **4.4 h × 3** | 2,270 = 2.6 h/seed · **7.7 h × 3** |

**None of these fits the 2.2 h remaining on this allocation**, and the 1,600-step run in flight is
1.8 h.

⚠ **And the 300-step run was worse than its refusal showed.** `min_tau_reliability = T/tau_int` was
`300/35 = 8.6`, and `observe/stationarity.py`'s own docstring says that below 50 the `tau_int`
estimate is biased LOW — so `n_effective` is OVERESTIMATED. The reported `n_eff = 4.28` is an upper
bound, not a value. C surfaced both fields (`fca6067c`) because the module had computed them and was
not saying them.

### 5.1 The run in flight stays — its PURPOSE changes

It cannot pass the contract and it was launched expecting to. It is kept because it is **the
measurement that decides the next length**, and that is worth more than a pass:

⚠ **τ_relax is 227 steps but the measured `tau_samples` is 35.** Two readings, and they differ by a
factor of six in what the next run costs:

* **(a)** the series has not reached its slow modes yet, so τ will GROW — the target recedes.
* **(b)** **γ's correlation time is genuinely shorter than a single node's τ_relax.** γ is a GLOBAL
  quantity — ~8e6 elements summed across 64 planes — so individual nodes' slow modes partly cancel.

**Neither can be settled by argument, and one measurement decides it: does `tau_int` stop growing as
the series lengthens?** `tau_int` at n = 1,600 against n = 300 is exactly that comparison, and
`min_tau_reliability` says when the estimate becomes trustworthy.

**So the next run's length will be derived from this one rather than estimated** — which is the same
discipline §3 applied to the mobility, one level up. Running three seeds now, without knowing τ,
risks all three being short for the same reason.

### 5.2 ⚠ The (e) 1 threshold is not a preference — it is the price of one series

`min_tau_windows = 50` here against the engine's `min_windows = 5.0` is **open on the PI queue**, and
the table above is what it costs: **a factor of 10 in GPU hours for one stationary γ series.** That
table belongs on the decision card. The knob is not a taste in statistical rigour; it is the budget.

⚠ **And the two thresholds are coupled.** `2cbd2c4d` established `min_tau_windows >= (4/3)·sokal_c`,
so at `sokal_c = 5` the floor is 6.67 — **`min_windows = 5.0` is below it and the contract now
raises on it.** A PI ruling of 5 requires lowering `sokal_c` in the same decision. Also for the card.
