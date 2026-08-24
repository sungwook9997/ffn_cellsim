# Contract II for the τ run — written before this run, and it exists because contract I was sized wrong

**2026-08-21 07:15 KST, before any step of run II executes.**
[Contract I](TAU_RUN_CONTRACT_2026-08-21.md) stands as the record of what was attempted and why it
failed. **Everything in its §4 — the read order, the trailing windows, the resolution criterion — is
carried here unchanged.** Only the size changes, and the reason the size changes is a defect.

## 1. Why run I is being stopped

⚠ **It was sized from a rate measured on a different cell.**

| record | nodes | populations | s/step |
|---|---:|---:|---:|
| `gthermo_long.json` — what I sized from | 4,361,496 | **0** (core only) | **4.1212** |
| `ab_host.json` — the τ configuration, measured today | 4,558,554 | **11** | **6.4204** |

τ runs the full cell. At 6.4204 s/step, 9,000 steps is **16.05 h** against a ~10 h grant: it is killed
at **~6,392 steps**, and because snapshots accumulate in RAM and the record is written after the loop,
it saves **nothing**.

⚠ **This is the same error as the balance-gate table that compared two different cells** — which I
committed and corrected five hours earlier, in a document about exactly this. The correction did not
transfer, because I had stopped thinking of `s/step` as a quantity that has a cell attached to it.

## 2. What is NOT changing

**The instrument stays the host estimator.** `--gamma-device` was wired, measured at **0.1985 s/step
against 6.4204 — 32.3×** — and then **rejected on an A/B**:

```
host seed 7, run 1  vs  host seed 7, run 2   worst rel  2.17e-16     ← the engine is deterministic
host seed 7         vs  device seed 7        worst rel  1.54e-01     ← 15%
```

Two host runs of the same seed agree to machine precision, so the engine is deterministic and **the
device path is producing a different quantity.** Cause identified, in the estimator's own docstring:
`centre` *"defaults to the ACTIN centroid — see the module docstring on why the two differ"*, and I
passed the **origin**. ⚠ **The plane-sum equivalence that passed at native (3.87e-12) was a test of
the kernel, on inputs I supplied; it could not catch a caller handing it the wrong geometry.**

⚠ **A 32× speedup that changes the answer by 15% is not a speedup, and it is not being used.**

## 3. Run II

**One seed, γ recorded every step, host estimator, the same thermostat configuration as run I** —
cortex bound, `k_axial` 88000 pN/µm, κ 0.07 pN·µm², mobility 1e-06 µm/pN·s ⚠ UNSOURCED, T = 310 K,
full cell, no `--core-only`.

⚠ **`--core-only` would run at 4.12 s/step and buy 8,200 steps instead of 5,200. It is not being
used.** Shrinking the cell to fit the budget is choosing the model to make the arithmetic work, and
`CLAUDE.md`'s conclusion rule forbids concluding below native regardless.

**5,200 steps ≈ 9.27 h**, against 9.85 h of grant — a declared 0.5 h margin plus slack.

⚠ **The 50-τ bar now passes only if T\* ≤ 104 steps**, against run I's 195 and a last measured lower
bound of 188. **The budget got tighter, not looser, and it is stated before the run rather than
discovered after.** If τ has not flattened below 104 steps, run II answers *"unresolved at N = 5,200"*
— which is the second branch of §5 below and is still decision-grade.

## 4. Carried unchanged from contract I §4

**4.1 Read order** — `min_tau_reliability` → `max_tau_relative_error` → `required_window_samples`.
Nothing later may be quoted if something earlier refuses.

**4.2 Windows are TRAILING** — the last 1/8, 1/4, 1/2, and the whole series. Trailing because the
transient leaves the window as it grows; a leading window needs the transient's length, which is
unknown, which is the circularity again.

**4.3 RESOLVED** when τ(last half) and τ(whole) agree **within their combined Madras–Sokal 1σ**.

## 5. Both outcomes remain decisive

* **τ flattens at T\* ≤ 104 steps** → `n_eff` ≥ 50 → the first stationary γ series this engine has
  produced. ⚠ One seed, so **D-2 stays unanswered**.
* **τ still tracks N** → the required window exceeds what this hardware buys with a host-side
  observable, which is a decision-grade answer about the configuration and not a failure of the run.

## 6. What is fixed for run II that was broken in run I

`_snapshot` now writes an atomic `<out>.progress.json` — steps done, measured s/step, ETA, frames
held — so **the rate is observable while it runs.** ⚠ It says in the file that it is a **heartbeat and
not a checkpoint**: the frames still accumulate in RAM and are still written after the loop, so a
wall-clock kill still loses them. Contract I §7 claimed otherwise and was false; this one does not.

## 7. What may not be done with run II

⚠ No magnitude — the mobility is unsourced, so τ in seconds is not a physical time and the criterion
uses τ in **steps**. ⚠ No seed-scatter claim. ⚠ No stationarity verdict except by §4.3. ⚠ **A truncated
run is analysed as a shorter run, never padded and never restarted-and-concatenated** — and run I's
frames are discarded entirely rather than joined to run II's, which is the same rule applied to the
run this one replaces.


---

## APPENDED 07:30 — run II FITS THE GRANT AND CANNOT ANSWER THE QUESTION

The heartbeat landed at **250 steps, 4.8619 s/step**, inside the window declared before the data
(4.80 – 6.82) and six minutes ahead of the declared deadline. **So run II fits**: 5,250 steps is
7.09 h against ~9.5 h of grant, with 2.4 h of margin. The rate is faster than both prior host
measurements (6.4204 beside run I, 8.0000 beside run II's own probes) because run II has the machine
to itself at 4441% CPU — which settles the contention question quantitatively.

⚠ **And it is being stopped anyway, for a reason the wall clock never mentions.**

§4.3 needs a discarded transient plus an analysed window — `7τ + 50τ = 57τ ≤ N`:

| | steps | **covers τ up to** |
|---|---:|---:|
| **run II** | 5,250 | **92 steps** |
| run I's last measured τ | | **188 steps, and that was a LOWER BOUND** |

> **Run II covers a τ ceiling less than half the smallest value τ has ever been measured at. It fits
> the budget and is guaranteed to return REFUSED.**

⚠ **This is a defect in contract II and it is mine — but not the version of it I first wrote.**

I recorded it as *"I should have computed the ceiling and did not"*. The criterion's author corrected
that: **`57τ ≤ N` did not exist in computable form when contract II was written.**
`required_window_samples` was only surfaced at 02:40, and the ceiling was first written down for
contract III. **So the specific accusation is false and is withdrawn.**

⚠ **What stands is the habit underneath it, and it needs no formula.** I sized by wall clock because
wall clock was the number I had, and I never stopped to ask whether it was the number that decides.
*"Fits the grant"* and *"can answer the question"* are different predicates, and the read order I had
already been given says in its own terms that the criterion consumes **τ**. **Noticing that they are
different does not require knowing the constant that relates them.**

**The same shape as sizing run I from a rate measured on a different cell**: right arithmetic, wrong
quantity, third time tonight.

It was found because the author of the criterion computed the ceiling for **contract III** before its
launch, and applying the same arithmetic backwards to II answered a question nobody had asked of II.
**Nobody did it for II because II was written in a hurry to replace I, and the hurry is the whole
explanation.**

**Run III supersedes: 3 seeds × 30,000 steps, covering τ ≤ 526.**
