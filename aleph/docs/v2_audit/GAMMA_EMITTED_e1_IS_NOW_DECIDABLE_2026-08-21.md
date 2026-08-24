# γ is emitted on the resting path, and (e) 1 is no longer undecidable

**Status:** MEASUREMENT + the removal of a blocker. Native, `cuda:0`, 4,361,496 nodes with the cortex
bound. Artifacts `world_phase4/phase4_gamma.json` and `phase4_gamma_relax.json`.

⚠ **NO γ VALUE HERE IS QUOTABLE, and that is not a hedge — it is the point.** `STATE.md` (c) 3 and
(c) 17 retire every γ this engine has produced; the coefficients in these runs are declared test
points; and one run cannot supply the **seed scatter** that (e) 1's D-2 requires. What is claimed is
about the OBSERVABLE's behaviour, not its magnitude, and that distinction is what makes the claim
possible at all.

---

## 1. What was blocking (e) 1

`PI_DECISION_e1_STATIONARITY_2026-08-16.md` §2.2 states it exactly: *"A stationarity gate cannot be
written against an observable the run does not produce."* Cortical tension γ is that observable — it
is what a resting cell IS mechanically and what every external measurement this project compares
against reports — and it existed only in the GATE-B slice drivers. **The amendment was not pending.
It was unwritable.**

Session C landed `world/observe_gamma.py` on 2026-08-21; this run is it, on the resting path, at
native, feeding `StepContext.observables`.

## 2. What was measured

Six steps, two arms, same build and same coefficients. The second differs only by a declared
numerical relaxation of 1e-8 µm/pN.

| step | frozen geometry | relaxation on |
|---:|---:|---:|
| 0 | 0.18859385 | 0.18859378 |
| 1 | 0.18859385 | 0.18858507 |
| 3 | 0.18859385 | 0.18856764 |
| 5 | 0.18859385 | 0.18855022 |
| **spread** | **0.000e+00** | **4.356e-05** |

## 3. ⚠ The result, and it is about noise floors rather than magnitudes

**Frozen, γ's spread is EXACTLY ZERO.** Not small — zero, to every digit float64 holds, across six
steps.

Compare the residual, from the same two runs: its spread is **2.111e-13 pN even when nothing moves.**
That floor is atomic-add ordering in the device reduction, and it is irreducible.

> **γ has no noise floor of its own; the residual does. A criterion written on γ therefore does not
> inherit the ~1e-13 pN floor that a residual-difference criterion must clear.**

And γ **responds**: under relaxation it falls monotonically over every step, a signal of 4.4e-05
against a noise floor of zero. Monotone, not scattered — the estimator is tracking the configuration
rather than sampling it.

**Those two facts together are what a stationarity criterion needs**, and they are exactly what the
residual cannot offer: the residual is noisy when nothing happens and blind when something does
(`BALANCE_GATE_MEASURED_2026-08-21.md` — 8.2 M added force terms left it unchanged while its own
tolerance grew 86×).

## 4. What is now decidable, and what is still the PI's

**Decidable:** (e) 1 can now be WRITTEN. There is an observable on the resting path, it is
deterministic, and it responds to the thing a stationarity gate is meant to detect.

**Still the PI's, and none of it is implied by the above:**

1. **Which observables define the resting state, and in what band** — (e) 1's D-1. γ is a candidate,
   not the answer; radius, volume, bound-myosin fraction and cortex thickness are the others.
2. **The statistic and its window** — D-2.
3. ⚠ **The variance it is judged against.** `STATE.md` (f) carries a standing rule that γ is compared
   against **SEED SCATTER, not within-run `sem`**: three replicates give S = 1.95% of γ, **67× the
   sem**. Session C's module has **no entry point that returns a single run's error bar** — it accepts
   only seed-keyed replicates — because D-2 calls within-run variance *"the single most likely way to
   write the amendment and still be measuring nothing"*. **This run is one seed. It cannot answer D-2
   and does not try to.**
4. **Stable versus conditional** — D-3, perturbation return and initial-condition independence. No
   single run answers it.

## 5. Two things the run also settled, recorded so they are not re-derived

**The estimator reads the force field that is actually integrated.** `force_mask` is a required
argument with no default in session C's module, because 2026-07-29 measured that masking the SOLVE
and not the MEASUREMENT left the dynamics bit-identical while γ differed by a constant
**0.41526 pN/µm** — the estimator was reading a different force field from the one being integrated.
This run passes `force_mask=()` explicitly.

**γ costs a host readback.** 8.34 s/step with γ against 0.16 s/step without: **a 50× step cost**,
because the estimator takes host arrays by design (`world/` may not import `engine/`, and that
constraint is what lets the whole path self-test on a machine with no CUDA). For a stationarity
window of thousands of steps that is the wrong side of the trade, and moving the method-of-planes cut
onto the device is the obvious follow-on. **Not done, and not needed to make (e) 1 decidable** — it is
needed to make it cheap.

## 6. Provenance

* build: `world/observe_gamma.py` (session C), `world/strand_bind.py`, `world/step.py`, 2026-08-21
* device: `cuda:0`, RTX 4090, Slurm job 72
* cortex bound: 4,128,840 axial segments + 4,060,026 bending triples at the BUILT chord
* every step: `ACCEPTANCE_UNDEFINED`. **γ being emitted does not accept anything** — the predicate is
  still `UndefinedAcceptance`, and it will stay that way until the PI answers §4.
