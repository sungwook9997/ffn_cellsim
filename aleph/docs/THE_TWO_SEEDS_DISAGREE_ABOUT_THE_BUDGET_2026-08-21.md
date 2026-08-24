# The two seeds disagree about whether 30,000 was enough, and the disagreement is not quantifiable

**2026-08-21, session `d380041d`. NOT a recommendation. NO run length is proposed here.**

Written before seed 3 landed, because the question the PI meets immediately afterwards is *"how long
does the next run need to be"*, and the answer is derivable **only as far as showing it is not
derivable.** Putting that arithmetic in front of the decision is different from proposing a number.

## What each seed says

**seed 1 — the only seed with any window above the `n_eff` bar of 25:**

```
transient ends            ~18,900 samples
τ over its VALID cuts      102 … 116        (n_eff 54.4 at 18,900 and 34.5 at 22,000)
transient + 50τ           ~24,700 samples   = 1.7 h at 0.2473 s/step
```

**Seed 1 says 30,000 was already enough, with ~5,300 samples of margin.**

**seed 2 — no cut anywhere clears the bar (`n_eff` 5.7 / 7.2 / 14.5):**

```
transient ends            NOT REACHED inside 30,000 samples
τ                         971 → 556 → 138, monotone, not converged
```

**Seed 2 gives no length at all.** Its τ is not a lower bound on anything: a τ measured on a window
that is still rising is not that window's τ. The tail slope is uniformly positive and *increases* as
the window shortens (+1.0e-03 → +1.4e-03), which is the same fact from the other side.

## The joint number the module prints — and why it may not be used as a budget

The three-gate read prints:

```
required_window_samples   136,371
window_samples available   15,375        ratio 8.9×
```

Taken at face value: `136,371 + a ~19,000 transient ≈ 155,000 steps ≈ 10.7 h per seed`, and
**32.0 h for three seeds** at the measured 0.2473 s/step.

⚠ **It may not be taken at face value.** `required_window_samples` is `50 × max_tau_samples`, and
`max_tau_samples` is **2,727 — seed 2's τ**, measured on a window that never settled. Sizing the next
run from it means sizing it from the quantity the next run exists to measure. **That is the
circularity two sessions have refused all night, and it does not stop being circular by arriving
inside a field instead of inside a sentence.**

The 32 h is written here only so that nobody arrives at it independently and treats it as clean.

## What is actually established

* **seed 1: 30,000 sufficient, with margin.** A measurement.
* **seed 2: 30,000 insufficient, by an amount this data does not bound above.** The absence of one.
* **The two do not average.** One is a measurement and the other is a failure to obtain one, and the
  mean of a number and a non-number is not a number.

⚠ **And the seed-to-seed difference is itself the finding, not noise around a budget.** `STATE.md`
(e) 1 D-3 asks for *"initial-condition independence: without this, stationary is compatible with
stuck where it started"*. Two seeds of the same configuration reaching stationarity at different
times — or one of them not reaching it — **is the thing D-3 asks about**, and it is a different
question from D-2's seed scatter. Which of the two this belongs to is the PI's call.

## What seed 3 can and cannot change

* **If seed 3 has a cut with `n_eff ≥ 25`** that also passes `opens` and `drift`: there are then two
  measurements and one absence, and the first genuine seed-to-seed τ comparison exists.
* **If it does not**: two absences and one measurement, and the case that 30,000 is inadequate for
  this configuration gets stronger — still without a number attached.

⚠ **In neither case does the three-seed report pass.** It requires all three replicates STATIONARY,
and seed 2 is not. That is by design (session 21's): dropping a replicate that did not settle and
averaging the rest is selection on the outcome.

**No length is proposed, in either branch.** The only seed that could support one is a single seed,
and sizing a scatter measurement from one member assumes the scatter it exists to measure.
