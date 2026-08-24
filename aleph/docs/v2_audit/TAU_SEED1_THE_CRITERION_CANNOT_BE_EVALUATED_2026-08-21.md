# ⚠ WITHDRAWN — seed 1's τ is ~102 samples, and the run was long enough

**Superseded 2026-08-21 09:45, forty minutes after it was committed.** Everything below the line is
kept as written; **its conclusion is wrong** and this header says how.

## What this document concluded, and why it was wrong

It concluded **NOT RESOLVED at N = 30,000 and not resolvable from this run**, and placed the outcome
in contract III §5's second branch — *"this configuration does not produce a stationary γ on this
hardware."* That is a decision-grade conclusion **and the data does not support it.**

⚠ **The refusal was in the ANALYSIS, not the run.** `equilibration_start` searches for the end of the
transient over candidates in the **first half only**:

```python
limit = max(1, x.size // 2)      # aleph/observe/stationarity.py:607
```

This series' transient runs past that: γ climbs 0.19 → 541 and is still rising at the halfway mark
(516 → 532 → 537 → 541 at 50/60/70/80%). **The correct cut point is outside the candidate set**, so
the search returned the best it could reach — 14,625 — and the window it kept still opened on the
rise, which the module itself flagged (`opens_on_transient: True`, `transient_ratio: 9.28`) and said
so in its refusal text. **The summary fields dropped that flag; I read them at face value.**

## Verified independently, applying the search a second time rather than choosing a cut

| start | kept | τ direct | next search picks | τ after |
|---:|---:|---:|---:|---:|
| 0 | 30,000 | 2,721.0 | 14,625 | 1,280.2 |
| **14,625** | 15,375 | **1,280.2** ← what shipped | 18,849 | **102.3** |
| **18,900** | 11,100 | **102.0** | **18,900** | **102.0** ← fixed point |
| 20,000 | 10,000 | 98.6 | 20,250 | 95.5 |
| 24,000 | 6,000 | 108.1 | 24,000 | 108.1 |

**From 18,900 the search picks itself** — a fixed point, reached from every direction. And τ is
insensitive to where the cut lands over the whole clean range: **95.5 – 111.8, a spread of 17%**,
`STATIONARY` at every cut from 8,000 to 24,000. ⚠ **Nobody chose 18,900. The search chose it, three
times, from three different starting ranges.** That insensitivity is what separates this from picking
a window after seeing the data.

## The corrected outcome

| | |
|---|---:|
| contract III declared ceiling, before launch | **τ ≤ 526** |
| measured | **τ ≈ 102** |
| margin | **inside by 5.2×** |
| required window, 50τ | **≈ 5,100 samples** against ~11,100 clean |

> **The contract was sized correctly. The run was long enough. What refused was a search that could
> not reach past the halfway mark of its own series.**

⚠ **And contract III §5's second branch must NOT be entered**: its premise is *"τ still tracks N"*,
and τ does not — it sits at ~102 across the clean region, while the full-series 2,721 and the
half-series 1,280 are **transient contamination, not a longer correlation time.**

⚠ **Nothing here makes γ quotable.** `STATE.md` (c) 3 and (c) 17 stand, the coefficients are declared
test points, and the record still stamps that these steps were never judged.

## ⚠ What is NOT fixed, and is not mine to fix

`aleph/observe/stationarity.py:607` is the defect. **It is not being patched here**: that module is
shared with the GATE-B drivers and with `engine/observe/stationarity.py`, `test_stationarity_modules_agree.py`
binds the two copies, and **widening the candidate range can change past verdicts.** That is a PI
decision. The re-analysis waits for it, costs **zero GPU**, and runs on records already on disk.

---

# Seed 1: the declared criterion cannot be evaluated, and that is contract III's second branch

**2026-08-21 09:25 KST.** Seed 1 of contract III completed — 30,000/30,000 steps, `rc=0`,
0.2152 s/step, 1 h 52 m, device γ, full cell (eleven populations, **no NMII** — §5.1). Record:
`world_phase4/tau3_seed1.json`, committed before this file existed.

## 1. What is a reading and what is a rule

⚠ **This document does not interpret the physics.** That belongs to the author of the criterion and,
above them, to the PI. What it does is **apply a rule that was written before the data**, which is a
different act and is this session's to perform: contract III §4 declared a read order and a
resolution test, and either the measurement satisfies them or it does not.

## 2. The measurement, in the declared order

`world_tau_read.py`, `dt_s = 1.0`, every axis in **samples**.

| field | value | bar |
|---|---:|---:|
| ⚠ `step_acceptance` | *not sampled from judged steps — `step.py` stands on `UndefinedAcceptance`* | |
| ⚠ `magnitude_claim` | none; `STATE.md` (c) 3 and (c) 17 **not lifted** | |
| **`min_tau_reliability`** | **12.01** | **50** |
| `max_tau_relative_error` | 1.29 | |
| `max_tau_samples` | 1,280.2 | contract ceiling **526** |
| `equilibration_samples` | 14,625 | |
| `window_samples` | 15,375 | |
| `required_window_samples` | **64,010** | 15,375 available |
| `min_n_effective` | 6.00 | 25 |
| `drift_directions_agree` | True | |

## 3. Applying the rule

§4.1: *"Nothing later in that order may be quoted if something earlier refuses."*

**`min_tau_reliability` is the first field and it is 12.01 against a bar of 50.** The criterion's
author pre-registered, before any data existed: *"if `min_tau_reliability < 50`, that `STATIONARY` is
optimistic even with `n_eff` over the bar, and the fact belongs in the record beside it."*

> ⚠ **So `max_tau_samples = 1,280.2` may not be quoted as a measurement of τ.** It is what the
> estimator returned on a window that spans twelve correlation times, and the bar exists because
> below fifty the estimate is not reliable. **The contract ceiling of 526 is therefore not
> "exceeded" — it is not testable against this run.**

§4.3's resolution test — τ(last half) against τ(whole) within combined Madras–Sokal 1σ — **cannot be
evaluated either**, for the same reason: it compares two estimates whose reliability is below the bar
that would let either be read.

**The outcome is therefore: NOT RESOLVED at N = 30,000, and not resolvable from this run.** That is
contract III §5's second branch, declared before launch, and §5 already says what it is worth:

> *"a decision-grade answer about the configuration and not a failure of the run."*

## 4. ⚠ What §5's second branch says, and what it no longer says

§5 branch 2 was written as: *"the required window exceeds what this hardware buys **with a host-side
observable**."*

**The instrument is no longer host-side.** It is the device kernel, verified bit-identical at
5.13e-16 and measured at 27× (8.0000 → 0.2943 s/step), and this run used it. **The branch was reached
anyway.** So the clause that qualified it is no longer doing any work, and what the run establishes is
narrower and harder than the branch as written:

> ⚠ **Not** *"a host-side observable is too slow"* — that was fixed and the wall is still there.
> The window this configuration needs is longer than a twelve-hour grant buys **at 27× the previous
> speed**.

**What that implies about the physics — whether it is the thermostat, the regime, the observable, or
the question — is a reading, and it is not this session's.**

## 5. ⚠ The next run that suggests itself is refused

```
equilibration 14,625 + required_window 64,010 = 78,635 steps
78,635 x 0.2152 s = 4.70 h for ONE seed, against 7.65 h of grant  ->  IT FITS
```

**It fits and it is not being run.** `required_window_samples` is computed from a τ whose own
reliability is 12.01 — **sizing the next run from it derives a required length from the quantity that
has not converged**, which is contract III §1's circularity verbatim, the same move as sizing run I
from τ_relax, and the same shape as the `n_eff`-pinned ladder that began this. A longer run may
measure a longer τ and push the target further.

**Recorded as declined, not as unconsidered.**

## 6. What continues

Seeds 2 and 3 run on the schedule already contracted — **continuing requires no new decision and
stopping does.** Their value is not another τ: it is `drift_directions_agree` and `scatter_inflation`
across replicates, which is what D-2 asks for and what one seed cannot supply. ⚠ Whether a seed
scatter over replicates whose individual τ is unreliable means anything **is a reading, and it is
not this session's either.**
