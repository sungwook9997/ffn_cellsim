# `--engine-cortex` A/B — `ac/engine` computes the cortex in production, and the answer does not move

**Lane** `lead` · **driver** `ffn_sim/scripts/ac_gate_b_cortex_motor_native.py`
**Population (both arms, FULL NATIVE)** 70,686 cortical filaments / 494,802 actin nodes / 674,314 total
nodes / 8,840 NMII heads, membrane subdiv 7, `k_xb` = 1000 pN/µm, catch-slip, seed 0, 30 steps, A5000.
**Verdict** PASS. **`quantitative_claim_status`** BLOCKED — this is a refactor-equivalence check, not a
magnitude.

## The criterion, declared before the run

> With the flag ON, `ac/engine` computes the cortex channels and the incumbent omits them via
> `_accumulate_all(omit=)`. **ON and OFF must agree node-by-node, or the swap moved force instead of
> relocating it.**

## Result

| step | γ_total OFF | γ_total ON | \|diff\| |
|---|---|---|---|
| 4 | 3.916381 | 3.916381 | 0 |
| 9 | 4.287321 | 4.287321 | 0 |
| 14 | 4.321879 | 4.321879 | 0 |
| 19 | 4.321579 | 4.321579 | 0 |
| 24 | 4.324925 | 4.324925 | 8.9e-16 |
| 29 | 4.320347 | 4.320347 | 0 |

| observable | max relative difference |
|---|---|
| γ_total | 2.05e-16 |
| γ_network | 1.14e-16 |
| γ_source | 0 (bit-identical) |
| `n_bound` | 0 (bit-identical) |
| `max_pf_pn` | 3.79e-16 |
| `max_f_cortex_pn` | 0 (bit-identical) |

Everything sits at float64 atomic-arrival order. **The engine component and the incumbent produce the same
cortex, in a production run, on the full native population.**

## This gate FAILED first, and the failure is the more useful half of the report

The first ON arm disagreed by a **constant 0.41526 pN/µm** at every sample — while `n_bound`, `γ_source`
and `max_f_cortex_pn` were bit-identical throughout. A constant offset in one observable with identical
dynamics is not a physics disagreement: the trajectories never diverged, so the relocation was already
doing what it claimed and the **measurement** was reading a different force field from the one being
integrated.

**Cause.** `make_inner_solve(..., omit=omit_channels)` masked the incumbent's cortex channels for the
SOLVE. `_measure_gamma` called `_accumulate_all` with no mask. The engine's cortex arrives through the
`cell.myosin` hook, which is not an omittable channel, so an unmasked assembly summed the incumbent's
bending/crosslink **and** the engine's — double-counting in the estimator only. Fixed by threading `omit`
through `_measure_gamma`, `_dump_viz` and all three `_accumulate_all` sites.

**Why nothing else could have caught it.** No exception, no NaN, a plausible γ, and a bit-identical
trajectory that positively *argues* the swap was clean. It surfaced only because an A/B ran with its
criterion declared in advance. The invariant is now a static test
(`tests/ac/engine/test_measurement_uses_the_solve_mask.py`), which immediately found a third unmasked call
site that had been missed by hand — and whose own first assertion, "the driver exists where this test
looks for it", caught a wrong `parents[4]` that would have made the other four assertions pass vacuously
on an empty file list.

**A second, cheaper failure on the way here.** Re-run #1 was queued from memory rather than from the arm's
recorded command and dropped `--out`, `--membrane-subdiv 7` and `--k-xb 1000`. It burned a device slot,
exited 0, and wrote nothing — while the run index recorded an artifact that did not exist. `ffn_gpu.py run`
now probes the declared artifact afterwards and exits 5 if it is absent.

## ⚠ The scope of this PASS: 30 steps is 0.3 s, and this repo has retired a conclusion at that duration

Both arms ran **30 steps × dt 0.01 = 0.3 s of physical time**, against a bound-head lifetime of 1.456 s —
**21% of one turnover.** That is the same duration class as the 2026-07-28 GATE-B reading that was retired
for being 12% of a turnover, and as the `sf_motor` 0.20 s runs that turned out to sit inside an induction
period. For a *refactor-equivalence* check the transient is not obviously disqualifying — both arms see the
same transient — but two things are untested by 30 steps:

- whether round-off differences **amplify** over a settled trajectory rather than staying at 1e-16;
- whether the equivalence holds once the bound population has saturated, which it has not at 0.3 s.

**Both arms are therefore re-queued `--until-stationary`** (2026-07-29), and the outcome is declared here
first: *if ON and OFF still agree at atomic-arrival order once settled, the relocation claim holds
unqualified; if they diverge, this REPORT must say the equivalence is **transient-only** and the headline
above narrows accordingly.* Until that lands, the PASS stands **as measured at 0.3 s** and not beyond.

### ANSWERED — the equivalence holds at stationarity, and the differences do not amplify

The ON arm ran to STATIONARY at **800 steps / 8.0 s physical time**, 27× longer than the original A/B, and
settles at γ_total = **4.32563718547168** — the same 15 digits as the incumbent `k_xb`=1000 point.
Compared step by step over all 160 samples of the settled trajectory:

| observable | at 30 steps | **at 800 steps (settled)** |
|---|---|---|
| γ_total | 2.05e-16 | **2.06e-16** |
| γ_network | 1.14e-16 | 2.07e-16 |
| `max_pf_pn` | 3.79e-16 | 1.40e-15 |
| γ_source, `n_bound` | bit-identical | bit-identical |

Everything stays at float64 atomic-arrival order; nothing grows with trajectory length. `max_pf_pn` is the
largest at 1.4e-15 and is expected to be — it is a max over 494,802 nodes of a residual, the most
sensitive quantity available. **So the relocation claim holds unqualified**, which is the first of the two
branches declared above.

**The comparison is against three runs on three builds, not two.** ON (`8615b917`), the incumbent
`k_xb`=1000 point (`6a5ef829`) and its cross-build rebuild (`8941bcb5`) all settle at the same 15 digits.
The only configuration difference between ON and OFF is `--engine-cortex` and a `--max-wall-minutes` cap
that neither run reached.

### The OFF arm did not run, and the reason was my own probe

`settled_off` was refused with exit 3 — duplicate `(build, config)` — because I had earlier run that exact
command through `ffn_gpu.py run` as a **diagnostic probe** and then killed it. The probe registered the
`(build, config)` pair in the run index and wrote **nothing**, so the real run was refused by a guard doing
precisely its job. That index row is the very defect the post-run artifact probe was added for this
morning: it names `outputs/ac/engine_cortex_ab/settled_off/record.json`, and no such file exists.

**Using `ffn_gpu.py run` as a diagnostic is the mistake**, not the guard's behaviour. A probe should not
consume the duplicate-run slot of the experiment it is probing. The comparison above did not need the
fresh OFF arm — the incumbent point is one, on a build already shown bit-identical — so the run was not
re-forced; forcing past a duplicate guard to repair a self-inflicted collision would have been the worse
of the two options.

## What this does and does not establish

- **Does**: `ac/engine` fires production kernels on the full native population and the composed answer is
  unchanged to round-off. The one-channel handoff mechanism (`_accumulate_all(omit=)`, PI
  `COMPARTMENT_VALIDATION_TRACKS` §6 D2) works end to end.
- **Does not**: end co-location. The arrays are still shared — the component computes the cortex channels
  but writes into the incumbent's `f`, because a genuinely private force array needs the inner solve to
  relax two arrays together, and that is `SOLVE_COUPLED`, still an empty phase (`STATE.md` (c) 4).
- **Does not**: make any γ magnitude quotable. Every NMII parameter is a PI-GAP and every step is
  force-accepted.

## Figures

Neither arm writes figures: the comparison is a table of two trajectories, and a plot of two curves that
coincide to 1e-16 conveys less than the numbers above. The `--dump-viz` path now carries the same omit
mask, so a rendered scene's γ title matches the solve it depicts.
