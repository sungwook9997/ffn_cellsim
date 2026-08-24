# The seed-scatter report refused twice, and the two refusals are the finding

**Status:** MEASUREMENT + a blocker that only became visible once γ existed. Native, `cuda:0`,
4,361,496 nodes, cortex bound. Six runs: three seeds × two arms. Artifacts
`world_phase4/gamma_seed{11,12,13}.json` and `gfrozen{11,12,13}.json`.

⚠ **NO γ MAGNITUDE IS CLAIMED.** `STATE.md` (c) 3 and (c) 17 stand. What is reported is the SPREAD
across independent builds and the two reasons a stationarity report could not be produced from it.

---

## 1. What was attempted, and why

`PI_DECISION_e1_STATIONARITY_2026-08-16.md` D-2 requires γ to be judged against **seed scatter**, not
within-run `sem` — three replicates give S = 1.95% of γ, **67× the sem**, and building a criterion on
within-run variance is called *"the single most likely way to write the amendment and still be
measuring nothing"*. Session C's module enforces that: it accepts only seed-keyed replicates and has
**no entry point that returns a single run's error bar**.

`build_all` did not forward a seed, so three runs at three `--seed` values produced the SAME CELL.
Forwarded 2026-08-21 — **without it, three identical runs offered as a scatter would have been that
exact mistake with extra steps.**

## 2. Refusal one — DRIFTING at 418,065 σ

Arm A: relaxation on at 1e-9 µm/pN, 24 steps.

```
seed 11  DRIFTING   slope -1.722e-05 /s   drift significance 418,065 σ   (declared threshold: 3)
seed 12  DRIFTING                          372,623 σ
```

The report is refused **whole**, not with the drifting replicates dropped — dropping the ones that did
not settle and averaging the rest is selection on the outcome, and session C's module says so.

⚠ **The interesting part is the size of that σ.** A drift of 1.7e-05 per second, five decimal places
below the value, is detected at **four hundred thousand sigma**. That is not a mis-scaled test: γ's
own noise is essentially zero (§3), so the denominator of the significance is tiny and any real trend
is astronomically significant. **The instrument is quiet enough that its drift test is almost
unfoolable** — which is the property a stationarity gate most needs and the residual cannot offer.

## 3. Refusal two — TOO_SHORT, and this one is structural

Arm B: relaxation OFF, 20 steps, geometry frozen.

| seed | mean γ [pN/µm] | within-run spread |
|---:|---:|---:|
| 11 | 0.18858987 | 2.776e-17 |
| 12 | 0.18860334 | **0.000e+00** |
| 13 | 0.18862015 | **0.000e+00** |

**Across seeds: std(ddof=1) = 1.5169e-05 pN/µm = 0.008% of γ.**

The report still refused: *"0.8 s of analysable window against a required 17.5 s = 50 × the measured
correlation time"*, `n_eff = 1.14` against a contract minimum of 25.

⚠ **And more steps would not fix it.** With the geometry frozen, γ is a CONSTANT — the within-run
spread is zero or one ULP. A constant series has no correlation time; the 0.35 s the assessor reports
is fitted to numerical dust. Running 350 steps produces 350 copies of the same number.

## 4. ⚠ The blocker this exposes: the engine cannot produce a STATIONARY γ series

The two arms are the only two this machinery has, and neither is stationary:

| arm | γ behaviour | verdict |
|---|---|---|
| `relax = 0` | exactly constant | **degenerate** — no fluctuation to be stationary about |
| `relax > 0` | monotone descent | **DRIFTING** at 4e+05 σ |

A stationary series is one that FLUCTUATES about a steady state. That requires the system to be
driven and damped — thermal forcing and a drag law — and **neither is bound**: the arena's per-node
`mobility` is still zero because a mobility belongs to a drag law nobody has written, and `relax` is a
declared numerical relaxation that says so in its own help text.

> **Emitting γ removed the blocker that made (e) 1 unwritable. It revealed the next one: (e) 1's
> criterion needs dynamics the engine does not yet have.**

That is progress and should be read as such — the previous blocker was *"we cannot write the gate"*
and this one is *"we can write it and have nothing stationary to run it on"*. It is also not a
surprise in hindsight: `STATE.md` (e) 1 itself records that a living cell is not in static mechanical
equilibrium and that force balance holds *statistically*. **A statistical balance needs statistics,
and statistics need noise.**

## 4.1 ⚠ CORRECTION, and it is worse than §3 said: MORE STEPS WOULD HAVE MADE IT PASS

§3 concluded *"more steps would not fix it"*. Session C went to check that and found the reverse.
**Measured, on the shared assessor:**

```
exactly constant, 6000 samples   ->  STATIONARY   mean 0.18860334   sem 0.0     n_eff 6000
```

**A frozen cell certifies as a perfect steady state with a zero error bar, if the series is long
enough.** The 20-step arm was refused for being SHORT (`n_eff` 1.14 against 25), not for being
constant — and 300 samples clears that bar.

⚠ **And `TOO_SHORT` tells the reader to do the thing that breaks it.** The refusal's own remedy is
*run longer*, and running longer converts a correct refusal into a false certificate. A gate whose
error message points at the failure mode is worse than one that is merely silent.

**This is `PI_DECISION_e1_STATIONARITY_2026-08-16.md` §1 objection (b) arriving by a route nobody was
watching.** Objection (b) is that *"the solver stopped moving"* must not be promoted to *"this is the
resting cell"*, and it was written about the SOLVER. It came in through the STATISTICS instead: the
solver never claimed anything, and the stationarity test certified the stop.

### The fix, and why it is not a threshold

`ecfd3f72`. Degeneracy is judged BEFORE verdict and length, against machine epsilon:

```
std <= eps * |mean|      ->  DEGENERATE, seeds named, and the refusal says
                             "MORE SAMPLES WILL NOT HELP: they make it pass"
```

**That is a machine property, not a tuning constant.** A threshold here would be a number deciding
which runs pass, which the charter forbids; `eps` is not chosen and is grid-invariant.

⚠ **A `zero_variance` flag would not have caught it.** A round-off constant has variance > 0, so the
flag reads False — and seed 11's `2.776e-17` spread is exactly that case. The comparison has to be
against `eps * |mean|`, not against zero.

**What it does to the three results above:** refusal ① (DRIFTING) is unchanged — that was a real
drift, correctly caught. Refusal ② now returns **DEGENERATE**: the same refusal for the true reason.
And the 0.008% scatter no longer passes at all — **what §5 says by hand, the module now refuses by
construction.**

## 4.2 A second correction: seconds were laundering an unsourced timescale

Also `ecfd3f72`, and raised by the thermostat's own gap. The report emitted `window_s`, `tau_int_s`
and `equilibration_s` in SECONDS — but seconds are `samples × dt`, and with an unsourced mobility the
run's `dt` is a step index wearing a unit. **The record was quietly converting a count into physical
time**, which is `CLAUDE.md`'s *"a magnitude quoted without the axis that says whether it may be
quoted at all"*.

`window_samples` now travels beside `window_s`, and the contract says which is which: *"seconds here
are samples × sample_dt_s and are only as sourced as the caller's dt_s. window_samples and
n_effective are counts and are always sourced — read those when the run's time axis is not."*

## 5. What the 0.008% is, and what it is not

**It is:** the build-to-build variation of γ across three independent cortex seeds at frozen
geometry — filaments placed differently on the same shell, same density, same discretisation.

**It is NOT the seed scatter D-2 asks for.** D-2's S = 1.95% was measured on runs that had dynamics;
this 0.008% is **240× tighter** precisely because nothing here fluctuates. Quoting it as the error bar
for a γ criterion would understate the real one by more than two orders of magnitude — the same
direction of error, and a worse magnitude, than the within-run `sem` mistake D-2 warns against.

**Recorded so it cannot be picked up as an error bar.**

## 6. For (e) 1, three things to carry into the wording

1. **γ, not the residual.** γ's noise floor is zero and it responds monotonically; the residual has an
   irreducible ~1e-13 pN atomic-add floor and is blind to 8.2 M added force terms.
2. ⚠ **A residual-based tolerance is loosened by the events it judges** (session A). A binding event
   ADDS a force term, a Higham-style tolerance grows with the term count, so more binding means an
   easier gate — and binding is what is being judged. If (e) 1 stays residual-based, the tolerance
   must be fixed against the declared **CAPACITY**, never the live term count. Occupancy is emergent
   and must not be a gate parameter.
3. **A stationarity criterion needs a driven-damped run.** Whatever statistic and window are chosen,
   there is nothing to run them on until a drag law and a thermal forcing are bound. That is a PHASE 2
   dependency of (e) 1 which was invisible while γ was missing.
