# Pre-registration — `n_inner` sensitivity, DOWNWARD

**Written and committed BEFORE the run.** A gate is a contract written before the data; this file exists
so the reading rule cannot be chosen after seeing it.

## The question, and why it is not already answered

Measured at full native, the outer step costs **1.51 s** and carries **43 inner iterations**
(`outputs/ac/profile/outer_step_profile.json`). 43 is a **budget**, not a convergence count:

* `STATE.md` tier-(a), negative result — at the resting configuration **66× the budget ends HIGHER than
  the start** (60 iterations −0.52%, 4,000 **+0.14%**). Raising it buys nothing.
* `ROADMAP.md`:93-95 — doubling the inner relaxation budget moves γ by **6e-08**. The observable does not
  depend on the inner budget in the UP direction.

**Nobody has measured the DOWN direction.** If the observable is as insensitive below 40 as it is above,
the budget is a straight multiplier available today at no modelling cost. That is the whole question.

## What is run

`aleph/scripts/ac_gate_b_cortex_motor_native.py` — an existing driver, an existing knob. No new physics
code is authored for this measurement.

| | |
|---|---|
| swept | `--outer` (the driver's name for `n_inner`) ∈ **{80, 40, 20, 10, 5}** |
| `--outer 80` | UP-control. Must reproduce the known insensitivity, or the sweep is not measuring what it thinks. |
| `--outer 40` | the incumbent default — the reference point |
| replicates | `--seed` ∈ **{0, 1, 2}** at every point |
| held fixed | `--steps 30 --dt 0.01 --filaments 70686`, all else default |
| population | full native (70,686 cortical filaments) |

## The projection, declared now

**Primary observable: `max_pf_pn`** — the whole-cell projected force residual max\|PF\| at the **final
accepted step**, read from the stored per-step ledger the driver already writes. Chosen before the run,
read from the field, not selected afterwards from what the run happened to keep.

**Secondary: `wall_seconds / accepted step`.** This is the payoff side and is expected to fall roughly
linearly in `n_inner`; it is reported, never argued from.

**γ is NOT an observable here.** `STATE.md` (c) 17 blocks every cortex-motor γ at `dt`=0.01 —
**magnitudes AND ratios** — so no γ appears in the result, in either form.

## The reading rule, declared now

Judgement is against the **across-seed scatter at `--outer 40`**, never a within-run `sem`
(`STATE.md` MEASUREMENT RULE; one seed cannot resolve a difference below ~5%).

* **INSENSITIVE** at a point ⇔ its 3-seed mean max\|PF\| lies within the `--outer 40` seed scatter.
* **SENSITIVE** ⇔ it lies outside.
* The lowest `--outer` that is still INSENSITIVE is the reportable result, *and nothing more than that.*

## What this run may NOT be used to claim

* **NOT that a lower `n_inner` converges.** Nothing here converges — see `STATE.md` (c) 4. Insensitivity
  means 40 is not buying what it costs; it does not mean 5 is enough for a physics conclusion.
* **NOT a speed number for the engine.** A timing is a property of an engine *at a stated convergence*,
  and this run states none.
* **NOT a new default.** Changing the production `n_inner` is a PI call; this measurement is its input.
* No physics magnitude of any kind, and no γ in any form.

`kind: diagnostic`.
