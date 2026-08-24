# Contract for the long τ run — written BEFORE the run, on a specification I did not author

**Written 2026-08-21 05:50 KST, before any step executes.** The charter says a gate is a contract
written before the run and may not be edited inline or re-thresholded after seeing the data. This
file is that contract. The specification is Lane C's (`0f47656e`); the execution is mine; **the
statistics choices are neither of ours and are not being made here.**

## 1. Why a long run at all — the previous τ was not a measurement

Three windows on the 1,600-step run returned τ = 1.139 → 2.575 → 9.384 s with `n_eff` pinned at 4.16.
Lane C read the pairing and it is decisive:

    n_eff = N / (2·tau_samples)     n_eff fixed at 4.16  ⟺  N/tau_samples = 8.32, constant

**τ was tracking the run length.** The estimator returns ≈ N/8.32 whatever the truth is, so those are
not three measurements of τ — they are three values of the *budget*. Deriving a required run length
from such an estimate is circular: the quantity has not converged, and it is the quantity being used
to decide how long to run.

⚠ **And my own sizing was wrong in the same direction.** I derived 1,600 steps from
τ_relax = 1/(µk), which is a **single bond** relaxing — the *fastest* mode. Equilibration waits on the
network's slowest collective mode. Fourth correction of the night with the same shape: the quantity I
sized on was not the quantity the question is about.

## 2. What breaks the circle — τ's own error bar

Lane C added it (`0f47656e`): Madras & Sokal 1988, `σ(τ)/τ ≈ sqrt(2(2M+1)/N)` with `M` the Sokal
window lag the estimator already carried. Validated against 200 AR(1) realisations — predicted
0.217/0.349 vs realised 0.171/0.257 at ρ = 0.9/0.98, so **~30% conservative**, i.e. reluctant to call
a change significant. For something a gate rests on, that is the correct direction to be wrong in.

Against that bar the observed 2.26× and 3.64× growths sit far outside ~1.5 relative 1σ. **The growth
is real and τ is unresolved** — now a field (`max_tau_relative_error`), not an argument.

## 3. The run

**One seed. γ recorded every step. Every step the grant will buy.** Job 73, RTX 4090-1, thermostat
configuration as `gthermo_long.json` (cortex bound, k_axial 88000 pN/µm, κ 0.07 pN·µm², mobility
1e-06 µm/pN·s ⚠ UNSOURCED, T = 310 K).

At the measured 4.06 s/step, 11 hours is **≈ 9,754 steps**.

## 4. Declared before the run — read order, windows, and the resolution criterion

**4.1 Read order.** `min_tau_reliability` → `max_tau_relative_error` → `required_window_samples`.
Nothing later in that order may be quoted if something earlier refuses.

**4.2 Windows are TRAILING, not leading**: the last 1/8, 1/4, 1/2, and the whole series. Trailing,
because the transient leaves the window as it grows and τ converges from below. A leading window
requires knowing the transient's length, which is unknown, which is the circularity again.

**4.3 Resolution criterion — Lane C's, not chosen by me and not chosen after seeing data:**

> **RESOLVED** when τ(last half) and τ(whole) agree **within their combined Madras–Sokal 1σ**.

A value inside its own error bar is the same value.

## 5. Both outcomes are decisive — that is the point of the specification

    4.06 s/step × 11 h = 9,754 steps.    last measured τ = 188 steps (a lower bound).
    the 50-τ bar passes iff  50·T* ≤ 9,754  ⟺  T* ≤ 195 steps

⚠ **188 against 195 — this budget buys exactly one window if τ has just flattened.** Not a
coincidence to celebrate; a knife-edge to state in advance.

* **τ flattens at T\* ≤ 195 steps** → `n_eff` ≥ 25 → **the first stationary γ series this engine has
  produced.** ⚠ One seed, so **D-2 remains unanswered** — seed scatter is a separate question and this
  run cannot speak to it.
* **τ keeps tracking N** → τ(9,754) ≈ 1,172 steps → required window ≈ **58,600 steps ≈ 66 h/seed, 198 h
  for three.** ⚠ **That is not a failure, it is a decision-grade answer**: this configuration does not
  produce a stationary γ on this hardware, and one of three things must change — the γ kernel moves
  device-side, the length bar moves (50τ vs 5τ, PI), or the physical regime does.

## 6. What this run may NOT be used for, stated now

* ⚠ **No magnitude.** The mobility is UNSOURCED (the arena's per-node mobility is zero because a
  mobility is a drag law's parameter and the builders decline to invent one). **The timescale of any
  series here is unsourced**, so τ in seconds is not a physical time. τ in *steps* is what the
  criterion uses, and that is deliberate.
* ⚠ **No seed-scatter claim** (D-2). One seed.
* ⚠ **No stationarity verdict by me.** The criterion above decides; if it is not met, the answer is
  "unresolved at N = ⟨actual⟩", not a judgement call.
* ⚠ **The step count is whatever the grant buys.** If the job ends early the analysis uses the actual
  N and says so. **The criterion does not depend on reaching 9,754** — it is an agreement test between
  two trailing windows of the series that exists.

## 7. Failure handling, declared in advance

Snapshots are written periodically so a truncated run is still analysable. **A truncated run is
analysed as a shorter run, never padded, never restarted-and-concatenated** — a concatenation across a
restart is not one trajectory and its autocorrelation is not the series' autocorrelation.

---

## APPENDED 2026-08-21 05:55, MID-RUN — an annotation, NOT an amendment

⚠ **Nothing above this line is changed.** The read order, the trailing windows and the resolution
criterion in §4 are exactly as committed before the first step, and they stay that way whatever this
note says. **A contract that can be edited once the run has started is not a contract**, and that
holds even when the edit would only fix my own transcription.

**What this note records**: §5's branch-2 wording — *"this configuration does not produce a stationary
γ on this hardware"* — **is false as written, and was false when I wrote it.**

96% of the step is a host-side observable
([`GAMMA_IS_96_PERCENT_OF_THE_STEP_2026-08-21.md`](GAMMA_IS_96_PERCENT_OF_THE_STEP_2026-08-21.md),
measured at 05:45 while this run was executing). The true statement is **"does not produce one with a
host-side observable"**: the same 58,600-step series at the no-γ step cost is **2.73 h for one seed**,
not 67 h.

The author of the specification accepted the correction in the same terms — *"I transcribed a budget
limit as a physical conclusion"* — and it is theirs and mine jointly, since I committed the wording.

⚠ **This changes nothing about the run.** The trajectory is identical, τ in steps is unaffected, and
the §4 criterion decides exactly as declared. It changes **only what may be concluded** if branch 2 is
the one that lands. Recorded here so a reader of this file alone is not misled, and dated so the order
is checkable: **the correction predates the data.**


---

## APPENDED 06:45, MID-RUN — §7 is FALSE about the run it governs

⚠ **§7 says "Snapshots are written periodically so a truncated run is still analysable." That is not
what the driver does, and I wrote the sentence without checking.**

`--snapshot-every 200` appends frames to a list **in host memory**; the record — including every
frame — is written **after the step loop returns**. So a run killed at the wall clock does not leave a
shorter run. **It leaves nothing.**

Which makes the second half of §7 wrong in the same breath: *"A truncated run is analysed as a shorter
run, never padded, never restarted-and-concatenated."* The prohibition still stands and is still
right; what is false is the premise that there would be anything to analyse.

### What this costs the run in flight

| rate | source | total | still needed at 06:40 | grant left | verdict |
|---|---|---:|---:|---:|---|
| 4.06 s/step | measured in a 300-step run of **this exact configuration** | 10.15 h | 9.14 h | 10.37 h | **fits, 1.23 h margin** |
| 6.14 s/step | this morning's isolated timing of the same two calls | 15.35 h | 14.34 h | 10.37 h | **killed at ~6,672 steps, nothing written** |

⚠ **The discrepancy I recorded as open now has a consequence.** Under the run-based number it
finishes; under the isolated number it produces nothing at all. I did not know which, and **had no way
to find out**, because the job emits no progress of any kind — no file, no line, only a process.

### Why it was not restarted

Restarting with a fixed driver costs the hour already spent, leaving **9.36 h for a 10.15 h run** —
which converts a run that probably fits into one that certainly does not. **The current run continues,
blind, at the risk stated above.** That is the choice, made explicitly rather than by default.

### What changed instead

`_snapshot` now writes an atomic `<out>.progress.json` on every snapshot — steps done, measured
s/step, ETA, frames held — and says **in the file** that it is a heartbeat and not a checkpoint, so no
one reads it as the durability §7 promised. ⚠ **The job now running predates that change and remains
unobservable.**

⚠ **§1–§5 are untouched.** The read order, the trailing windows and the resolution criterion are
exactly as committed before the first step. This appendix corrects a statement of FACT about the
driver; it does not move a threshold, and it could not, which is the whole reason the two are written
in separate sections.
