# `n_inner` sensitivity, DOWNWARD — result

Pre-registered before the run: [`PREREGISTRATION.md`](PREREGISTRATION.md) (`a3abe347`).
15 runs = 5 budgets × 3 seeds, 30 accepted steps at `dt` = 0.01 s, 70,686 cortical filaments
(full native), RTX 4090-1, Slurm job 81. Zero failed runs. `kind: diagnostic`.

## 1. The pre-registered verdict: **do not cut**, and the rule cannot be trusted either way

| `--outer` (`n_inner`) | max\|PF\| final step [pN] | shift from 40 | pre-registered verdict |
|---:|---:|---:|---|
| 80 (up-control) | 0.767527 | 0.004508 | SENSITIVE |
| **40** (reference) | **0.772035** | — | INSENSITIVE |
| 20 | 0.774538 | 0.002503 | SENSITIVE |
| 10 | 0.813476 | 0.041441 | SENSITIVE |
| 5 | 0.887813 | 0.115778 | SENSITIVE |

Lowest INSENSITIVE budget = **40**, i.e. the declared rule says cut nothing.

**The up-control also came back SENSITIVE, so the rule failed its own control.** The seed scatter at
`--outer 40` is `S₄₀ = 1.09e-4 pN` — **0.014% of the value**, two orders tighter than the γ scatter
`STATE.md`'s MEASUREMENT RULE was written from (1.95%). A band that narrow separates changes with no
physical meaning, so it returns SENSITIVE for everything and its INSENSITIVE is not informative either.

**The rule was NOT changed after seeing this.** Re-thresholding a gate against its own data is what the
charter forbids, and the pre-registration exists precisely so that this temptation is refusable. The
verdict above stands as declared. What this run additionally produced is the finding that
**an across-seed scatter is the wrong band for this observable** — the next pre-registration that reaches
for one must justify it against a physical scale, not inherit it from the γ rule.

## 2. The magnitudes — declared as the primary observable, and they show a knee

Relative to the `--outer 40` reference:

| `--outer` | change vs 40 |
|---:|---:|
| 80 | −0.59% |
| 20 | **+0.32%** |
| 10 | +5.37% |
| 5 | +15.00% |

The residual is **flat from 80 down to 20** (0.9% across a 4× budget change) and then **rises steeply
below 20**. So the budget is buying nothing between 80 and 20, and buying something real below it.

**This refutes the hypothesis the sweep was built to test.** The session's working claim was that
`n_inner` could go 40 → 5 nearly free, worth 4–8×. It cannot: 5 costs **+15%** in the observable.
What survives is the smaller, measured statement: **40 → 20 costs 0.32%**.

## 3. Cost — NOT from this sweep

This sweep's wall-clock cannot measure the lever and is reported only for completeness (6.889 → 6.426
s/step from `--outer` 80 → 5, i.e. **6.7% for a 16× budget change**). The driver's own `device_note`
says why: telemetry (γ readback + host force norm) runs between steps, here at the default
`--gamma-every 1`, and it dominates the per-step wall-clock. Leaving that default was a design fault in
the sweep, not a property of the engine.

The cost side is already measured elsewhere, on a telemetry-light configuration:
`outputs/ac/profile/` gives 29.66 ms per inner iteration and 1513.76 ms per outer step at 43 inner
iterations — inner work is **84.2%** of the step there. Halving the budget therefore removes ~638 ms of
1514 ms: **≈1.7×**.

**So the honest size of this lever is ~1.7× at +0.32% in the observable** — real, and an order smaller
than the 4–8× claimed before it was measured.

## 4. What may not be taken from this

* **Not that `--outer 20` converges.** Nothing here converges (`STATE.md` (c) 4). "Flat" means 40 is not
  buying what it costs, and says nothing about whether either budget is sufficient.
* **Not a speed number for the engine.** Every record carries `timing.comparable: false`: a timing is a
  property of an engine at a stated convergence, and no run here states one.
* **Not a new production default.** That is a PI call; this is its input.
* **No γ**, in magnitude or ratio form — `STATE.md` (c) 17.

## 5. Provenance — two gaps, both recorded rather than repaired

1. **The 15 records carry no build commit** (`build.commit: "unknown"`, `source: "unavailable"`):
   `ac_gate_b_cortex_motor_native.py` has no `--commit` flag and the run host is not a git checkout.
   What *was* verified: the driver's whole first-party import closure, **120/120 files hash-matched**
   against the working tree immediately before launch.
2. **A file in that closure changed mid-sweep.** `aleph/engine/composed_native.py` was re-synced after
   run 6 of 15 — so `--outer` 80 and 40 ran against one working tree and 20/10/5 against another, with
   the boundary falling exactly on the reference point. **No executed code differs**: both changed hunks
   are inside `compose_native_cell_world` and `build_native_composed_cell_world`, and the only callers of
   either are `ac_connector_devicerun_census.py` and `tests/ac/engine/test_composed_native_smoke.py`.
   This driver calls neither, and no module-level code changed. Stated here rather than left for a reader
   to discover, because a closure hash taken before launch does not cover a file replaced during it.
