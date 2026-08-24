# ALEPH-PORT-3649 — the line search pays a retraction it never charged the reference

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3649` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-07` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **INTERNAL FIX.** Four lines in an acceptance test. No new physics, no new law. |
| Aleph target | `aleph/runtime/resident_cortex.py` |
| Depends on | `ALEPH-PORT-3641` (the CG driver this fixes) |
| Exists because | CG stops at the same iteration whatever the budget and whatever the backtrack count, and the reason is a rounding the acceptance test does not charge symmetrically. |

---

## 0. The measurement, and it is exact

ρ = 8, `n` = 5, filamin, END_CORRECTED, device-side filament-pair steric. CG driven to iteration
4,063, then **one more iteration by hand** with the direction forced to steepest descent — the safest
direction that exists — and the line search's own trial sequence printed:

```
E0 (the stored reference)              = 0.15230932101996544
<d,r> = 1.096444e-01                   (positive: it IS a descent direction)

  trial step        trial energy                       dE      accept?
  1.0000e-09        0.15230932091088401       -1.0908e-10        YES
  1.0000e-11        0.15230932101929801       -6.6744e-13        YES
  1.0000e-12        0.15230932102036207       +3.9663e-13         no
  1.0000e-16        0.15230932102048861       +5.2316e-13         no
  1.0000e-26        0.15230932102048861       +5.2316e-13         no    <- the limit
```

**As the step goes to zero the trial energy does not converge to `E0`.** It converges to
`0.15230932102048861`, which is **`5.2316e-13` above it**, and stays there for fourteen further
decades of step size.

## 1. What that gap is

`k_step_along` computes `base + trial · direction` and then **re-pins every node to its stored
radius** — the retraction that keeps the relaxation on the constraint manifold. At `trial = 0` that
is arithmetically the identity, and in float64 it is **not**: re-normalising a vector and
re-scaling it by a radius does not reproduce the bits of a point that was itself produced by an
earlier re-pin along a different rounding path.

So the acceptance test

```python
if trial_energy <= energy:   # energy: measured WITHOUT a retraction
                             # trial_energy: measured WITH one
```

**compares two different operators**, and charges the trial side a fixed penalty of about `5e-13`
that the reference side never pays. The penalty does not shrink when the step does, because it is
not a function of the step.

**So the line search can accept only while the true decrease exceeds the retraction rounding.** Once
per-step progress falls below it — which it must, as the iterate approaches a minimum — no step is
acceptable, at any size, and CG stops.

## 2. This explains three things that were previously unexplained

* **`ALEPH-PORT-3641`'s amendment**: *"CG stops at 15,300 however large the budget."* Measured again
  here: identical energy, to the last bit, at budgets 4,064 and 20,000.
* **Backtracking does not help.** Measured at 30, 50, 80 and 120 backtracks:
  **`0.1523093210199654` at every one, bitwise, and 4,064 iterations at every one.** A fixed penalty
  is not something a smaller step can escape.
* **`ALEPH-PORT-3648` K-5**: the iteration count *falls* as the mesh is refined — 34,583 → 11,359 →
  **4,357** at ρ = 24 for 22,440 → 37,400 → 67,320 nodes. More nodes means more retraction rounding
  in the same energy sum, so the penalty grows while per-step progress does not. **The finer mesh
  stops sooner because it is noisier, not because it is easier**, and K-5's whole refinement drift
  is therefore suspect.

## 3. The fix

Measure the reference **through the same retraction**. Before the trial loop, apply `k_step_along`
with `trial = 0`, re-assemble, and use *that* as the value to beat.

```
base <- positions
step_along(base, direction, 0.0)  ->  positions      # pay the retraction once, on the reference
assemble_forces();  energy <- E                       # now both sides have paid it
... trial loop unchanged ...
```

Both sides then carry the same rounding and the comparison is between two configurations rather than
between two operators. It costs **one extra force assembly per iteration** — about 1 in 30 of the
line search's own cost — and no host traffic, so `ALEPH-PORT-3637`'s residency invariant is untouched.

**What is deliberately NOT done.** Switching the acceptance test to an Armijo condition on
`⟨d, r⟩`, or to a residual-norm decrease. Both are defensible and both change *what the method is*;
this entry changes only what it compares, and if the stall survives that then the problem is the
method and not the arithmetic, which is a different entry.

## 4. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **L-1** | **The reference is the limit.** At the configuration in §0, the trial energy as `trial → 0` equals the stored reference. | agree to `< 1e-15` relative, where they now differ by `5.2e-13` absolute |
| **L-2** | **It gets past 4,064.** Same configuration, same budget of 60,000. | iterations `> 4,064`, and the RMS residual **lower** than 0.002966 |
| **L-3** | **Budget starts to matter again.** Energy at budgets 20,000 and 60,000. | they differ, where today they are bitwise identical |
| **L-4** | **The mesh ordering inverts.** Iterations at `n` = 3, 5, 9, ρ = 24. | **rises** with node count, where today it falls 34,583 → 4,357 |
| **L-5** | **The energy is still monotone across accepted steps.** | asserted, no step raises it |
| **L-6** | **Residency.** Mesh-sized host crossings between `push` and `pull`. | exactly `0` |
| **L-7** | **What it costs.** Wall clock per iteration, before and after. | **measured, no threshold** |
| **L-8** | **What it moves.** Relaxed `K_A` at ρ = 8, 24, 100, before and after. | **measured, no threshold** |

**L-2 and L-4 are the port's reason and can fail honestly.** If CG still stops early with the
reference corrected, then the retraction rounding was real but not binding, and the stall is the
method — an ill-conditioned nonlinear CG that needs a restart strategy or a Newton–CG, which
`-3641` §2 already named as the follow-on.

**L-8 carries no threshold.** Every relaxed modulus in this repository was measured with the solver
stopping early. If they move, the ones already published move with them, and that is a result to
report rather than a size to predict — this lane has predicted the size of an effect and been wrong
eight times.

## 5. What is NOT claimed

- **Not that the relaxed moduli are wrong.** They were measured at a stopping point that was not a
  minimum. Whether the minimum is elsewhere is L-8.
- **Not that the retraction is wrong.** Re-pinning to a stored radius is the constraint, and it is
  correct. What is wrong is comparing an energy that has been through it against one that has not.
- **Not that this is the last solver defect.** `-3641` §2 named truncated Newton–CG as the follow-on
  and this does not deliver it.
- **Not that `-3641`'s existing findings are withdrawn.** The descent guard `⟨d,r⟩ > 0` and the
  vector transport are both still necessary; §0's run has them and still stalls.

---

## Amendment, 2026-08-07 — the fix works, and two of its own gates were written against a guess

CPU-Warp, ρ = 8, `n` = 5, filamin, END_CORRECTED, device-side filament-pair steric.

| | before | **after** |
|---|---:|---:|
| iterations | 4,064 | **13,645** |
| energy | 0.1523093210199654 | **0.1517941752616385** |
| RMS residual | 0.002966 | **0.0002349** |
| L∞ residual | 0.2301 | **0.01639** |

**3.36× the iterations, 12.6× the RMS, 14× the L∞** — from measuring the reference through the same
retraction the trials go through, and nothing else.

### The gates

| # | outcome |
|---|---|
| **L-2** it gets past 4,064 with a lower residual | **PASS** — 13,645 and 0.000235 against 0.002966 |
| **L-5** the energy is monotone | **PASS** |
| **L-6** residency | **PASS** — 0 mesh-sized crossings |
| **L-7** what it costs | one extra force assembly per iteration, ~1/30 of the line search's own |
| **L-8** what it moves | `K_A` at ρ = 8, `n` = 3: 175.2095 → **175.1137**, **−0.055%** |
| **L-3** budget starts to matter | **FAILED as written** |
| **L-4** the mesh ordering inverts | **FAILED as written** |

### L-3 and L-4 failed because they predicted the shape of the fix, not just its direction

**L-3** asked that the energy differ between budgets of 20,000 and 60,000. It does not — 13,645
iterations and the same energy at 20,000, 60,000 and 200,000. **That is the fix working better than
the gate expected**: CG now stops at a genuine stationary point, and more budget correctly changes
nothing. The gate assumed the stall would become budget-limited; it became *converged*.

**L-4** asked that iterations *rise* with node count. They do not — 16,531 / 13,645 / 13,741 at
`n` = 3 / 5 / 9. But the pathology it was written to catch is gone: before the fix the same three
were **11,465 / 3,495 / 11,955**, with the middle case stopping at a residual an order of magnitude
worse. Flat is not rising, and flat is what makes a comparison across `n` valid, which is what the
gate was for.

**Two more thresholds written against a guess, which makes ten.** They are recorded FAILED rather
than reworded. The pattern is now so consistent that it is worth stating as a rule this lane keeps
relearning: **a gate should say what must be true, not what the fixed version will look like.**

### And K-5's refinement drift SURVIVES the fix, which was not the expected outcome

`ALEPH-PORT-3648`'s K-5 amendment suspected the +0.740% drift at ρ = 8 was contaminated by the finer
meshes converging worse. With the line search corrected and all three residuals now comparable:

| `n` | before the fix | **after** |
|---:|---:|---:|
| 3 | 175.2095 | 175.1137 |
| 5 | 176.2494 | 175.4424 |
| 9 | 176.5053 | 176.3999 |
| **drift 3 → 9** | **+0.740%** | **+0.735%** |

**The drift is real.** It is not the solver. `-3648` K-5's suspicion is discharged and its question
now has an answer at ρ = 8: the relaxed areal modulus rises by about 0.74% between `n` = 3 and
`n` = 9, above the ±0.5% run-to-run spread, so the cortex's relaxed modulus is **not yet
mesh-independent**.

**One caveat and it is live.** A best-of-nine run under the *old* solver, at the same density, gave
`+0.229%` — below the spread and so not resolvable. Single-run-fixed-solver and best-of-nine-old-
solver disagree, so the definitive number needs both together, and that run is queued rather than
claimed.

### What every relaxed modulus in this repository was measured at

All of them were measured with the line search stopping early. At ρ = 8 that cost 0.055% in `K_A`
and a factor of 12.6 in the residual it was reported at — so the published *numbers* barely move,
and the *convergence claims* attached to them were about a stopping point that was not a minimum.
`ALEPH-PORT-3641`'s amendment already said the L∞ criterion had not noticed convergence; this says
the line search had not reached it.

---

## Amendment 2, 2026-08-07 — **RETRACTION.** "The drift is real. It is not the solver."

The amendment above concluded, from single runs before and after the fix:

> *"The drift is real. It is not the solver. `-3648` K-5's suspicion is discharged and its question
> now has an answer at ρ = 8: the relaxed areal modulus rises by about 0.74% between `n` = 3 and
> `n` = 9, above the ±0.5% run-to-run spread, so the cortex's relaxed modulus is **not yet
> mesh-independent**."*

**That is retracted.** It named its own caveat — *"a best-of-nine run under the old solver gave
+0.229%, so the definitive number needs both together"* — and then stated the conclusion anyway.

### Both corrections together, which is what the caveat asked for

Fixed line search **and** best-of-nine matched convergence, ρ = 8:

| `n` | nodes | `K_A` relaxed | iterations | rms |
|---:|---:|---:|---:|---:|
| 3 | 7,479 | 175.77 | 6,346 | 0.000357 |
| 5 | 12,465 | **175.52** | 13,552 | 0.000210 |
| 9 | 22,437 | 176.528 | 13,998 | 0.000271 |

**Drift `n` = 3 → 9: `+0.431%` — below the ±0.5% run-to-run spread, and therefore not resolvable.**

And it is **not monotonic**: 175.77 → 175.52 → 176.53 goes down and then up. A real mesh dependence
would not.

### What is actually known

| measurement | drift at ρ = 8 | verdict |
|---|---:|---|
| one run, old solver | +0.740% | contaminated — finer meshes stopped 3× earlier |
| best-of-9, old solver | +0.229% | not resolvable |
| one run, fixed solver | +0.735% | **this is what the retracted claim rested on** |
| **best-of-9, fixed solver** | **+0.431%** | **not resolvable, and non-monotonic** |

Three of the four say the drift is at or below the noise. The one that said otherwise was the one
with the smallest sample, and it is the one that got written up.

**`ALEPH-PORT-3648` K-5's question is still open at ρ = 8**, and the honest answer is *"the
measurement cannot tell"* rather than either "it converges" or "it does not". Resolving it needs
either a tighter solver spread than ±0.5% or a wider refinement range than 3 → 9.

### Nothing else in this entry moves

The line-search fix stands on its own evidence and none of it came from K-5: 4,064 → 13,645
iterations, RMS 0.002966 → 0.000235, L∞ 0.2301 → 0.01639, and the `5.23e-13` limit measurement in
§0 that diagnosed it. **What is retracted is a downstream inference, not the fix.**

**This is the eleventh time on this project that a conclusion has been drawn from an under-powered
measurement, and the third time today.** The pattern is specific enough to name: *this lane
generalises from the smallest sample that produces a clean-looking number.* The guard that keeps
catching it is running the same thing again with more of it.
