# Blebbistatin dose response, run to STATIONARITY — three points, no ratio yet

**Lane** `lead` · **driver** `ffn_sim/scripts/ac_gate_b_cortex_motor_native.py --until-stationary`
**Population (every point, FULL NATIVE)** 70,686 cortical filaments / 494,802 actin nodes, membrane
subdiv 7, 8,840 NMII heads, A5000. **Verdict** `quantitative_claim_status: BLOCKED`, `comparable: false`.

## What this replaces, and why the previous attempt was void

The 2026-07-28 sweep ran six points for a **fixed 60 steps** and reported that `k_xb=3` gave higher tension
than `k_xb=10` — a non-physical ordering. The cause was not the physics and not the engine: only `k_xb=1000`
had settled in 60 steps; the others were still rising. **A softer crossbridge relaxes more slowly, which is
correct physics — so when the swept parameter sets the relaxation time, a fixed step count is the wrong
protocol by construction.** Every point here therefore runs until its own settling criterion is met, and the
criterion is declared before the run and recorded in the artifact next to the verdict it produced.

## Result

| point | `k_xb` (pN/µm) | steps | wall | transient | τ_int | γ_total (pN/µm) | sem | drift/σ | bound |
|---|---|---|---|---|---|---|---|---|---|
| `kxb_1000` | 1000 | 800 | 38.7 min | 0.40 s | 0.025 s | **4.3256** | 0.00060 | 0.13 | 0.9823 |
| `kxb_100` | 100 | 800 | 38.7 min | 0.40 s | 0.073 s | **3.4807** | 0.00172 | 0.19 | 0.9823 |
| `kxb_10` | 10 | 1400 | 67.2 min | 3.90 s | 0.237 s | **2.8839** | 0.00312 | 0.98 | 0.9829 |

Contract, identical for all three and declared before each run: window ≥ **5.0 ×** the driving correlation
time (1.456 s, the catch-slip zero-load bound-head lifetime), `|drift| < 1.0 σ` over that window. `sem` is
`σ·√(2·τ_int/T)` — the error on a *correlated* mean, not `σ/√N`.

**Bound fraction is flat at 98.2–98.3% across two decades of `k_xb`.** This is the negative control and it
passes: the sweep changes how much force a bound head transmits, not how many heads are bound, so the γ
change is force transmission and not a kinetic side effect.

## Two disclosures the STATIONARY label hides

1. **`kxb_10` passed by 2%** — 0.98 σ against a 1.0 threshold, with `n_eff = 21`. It is a pass, not a
   comfortable one.
2. **At `kxb_10` the `γ_source` component did NOT settle** (1.18 σ) and the detector returned no mean for
   it. R is a `γ_total` ratio, so this does not block it — but the point is not "settled" without
   qualification.

Read together these say something about the criterion rather than the physics: at `k_xb=10` the three
observables sat at **0.92 / 0.98 / 1.18 σ**, so a 1.0 threshold sorted a set of physically similar states
into pass and fail. The analyser now prints drift/σ for every point and names both cases, because a
refusal the reader never sees is the defect this sweep exists to avoid, one level down.

## The provenance block was lifted by measurement — and the parity is stronger than the criterion asked

The three points carry three different build commits with no closure digest, so the analyser refused: not
because the physics was in doubt but because *whether the code differed could not be established*. A
cross-build parity run settled it. `k_xb=1000` was re-measured on build `8941bcb5` against the original
`6a5ef829`:

| observable | `6a5ef829` | `8941bcb5` | identical |
|---|---|---|---|
| γ_total | 4.32563718547168 | 4.32563718547168 | **yes, to the last bit** |
| γ_network | 4.303042713651803 | 4.303042713651803 | yes |
| γ_source | 0.022581453932423 | 0.022581453932423 | yes |
| bound fraction | 0.9823069852941175 | 0.9823069852941175 | yes |

The declared criterion (before the run) was agreement within **3 combined sems**. The measured agreement is
**0.00 σ — bit-identical on all four observables**, i.e. the two builds produced the same float64
trajectory, not merely a compatible mean. The commit difference was a false difference, exactly as the
refusal's own wording allowed for.

**The limit of that evidence, stated because it is easy to overclaim.** The pair covers `6a5ef829` and
`8941bcb5`. Builds `efe74e52` and `8dda13a4` carry the `k_xb`=100 and 10 points and were **not** re-run;
their equivalence is INFERRED from one pair spanning the range, not measured. The analyser prints this
every time it lifts the block.

## R — superseded twice in one night; read §"RESULT" and §"CONCLUSION" below instead

This section first read **R = 0.6667 ± 0.0007, band INCONCLUSIVE (an upper bound)**. Both the number and
the framing are now void: `k_xb` = 3 and 1 landed afterwards, γ turned around, and the sweep turned out to
run outside the model's own validity window. Kept only so the sequence of readings is visible — three
successive statements of "what R is", each retired by the next point that landed.

## Two points were dropped tonight, by arithmetic, not by preference

`k_xb = 1` and `k_xb = 0.1` were queued and removed before running. The transient grew **0.40 → 0.40 →
3.90 s** across the three decades measured, i.e. roughly ∝ 1/`k_xb` once it leaves the detector's floor.
Extrapolated, `k_xb=1` needs ≈ 4,600 steps ≈ 220 min and `k_xb=0.1` an order of magnitude more — against
`max_steps` 2,500 and a 75-minute wall. **Neither could return a verdict at any affordable setting**, and a
point that cannot settle contributes nothing to a dose response. `k_xb=3` was stopped 3 minutes in for the
same reason (≈ 2,300 steps ≈ 111 min > 75 min wall).

This is recorded rather than silently omitted: the sweep spans 2 decades, not 4, and the reason is a
measured cost, not a choice about which points would look better.

### That justification rests on two data points, so it is being tested — prediction declared first

"Transient ∝ 1/`k_xb`" is an extrapolation from **0.40 → 3.90 s**, i.e. two points, one of which sits at
the detector's resolution floor. It is load-bearing: it is the whole argument for dropping three doses. So
`k_xb = 3` is queued at a wall sized from that scaling rather than from the 75-minute default that made it
fail before, and the prediction is written here **before the run**:

| | predicted |
|---|---|
| transient | ≈ 13 s |
| steps to STATIONARY | ≈ 2,030 (cap 4,000) |
| wall | ≈ 105 min (cap 135) |
| verdict | STATIONARY |

**What each outcome means, also declared in advance.** *Settles near 2,030 steps* → the scaling holds and
dropping `k_xb` = 1 and 0.1 was justified. *Settles far sooner* → the scaling is wrong, the drop was NOT
justified, and `k_xb=1` goes back on the queue. *Does not settle by 4,000 steps* → the growth is worse
than 1/`k_xb` and this REPORT **understates** the cost of the low-`k_xb` end.

### RESULT — the prediction was REFUTED, and it took the sweep's framing with it

| | predicted | measured |
|---|---|---|
| transient | ≈ 13 s | **0.75 s** |
| steps | ≈ 2,030 | 2,000 |
| wall | ≈ 105 min | 96 min |
| verdict | STATIONARY | STATIONARY |

The step count matched by coincidence — the run was long because τ_int grew, not because the transient
did. **Transient measured across the sweep: 0.40 → 0.40 → 3.90 → 0.75 s. It is not monotonic in `k_xb` and
it does not scale as 1/`k_xb`.** So the argument that dropped `k_xb` = 1 and 0.1 was wrong, and `k_xb`=1
is back on the queue. Per the middle branch above, that is the outcome this run was built to force.

**And the larger finding, which is not about cost at all: γ(`k_xb`=3) = 12.2425 ± 0.0674 — 4.25× the value
at `k_xb`=10.** The dose response has a **turning point** inside the swept range.

| `k_xb` | γ_total | max\|PF\| | maxF_cortex | γ_source | signal/residual |
|---|---|---|---|---|---|
| 1000 | 4.3256 | 0.5730 | 6.5203 | 0.02258 | 7.5 |
| 100 | 3.4807 | 0.5291 | 6.2784 | 0.02694 | 6.6 |
| 10 | 2.8839 | 0.4999 | 5.2223 | — (drifting) | 5.8 |
| **3** | **12.2425** | 0.4914 | **2.3909** | 0.10911 | **24.9** |

**The residual explanation is disfavoured** — the one that retired (c) 3, so it was checked first.
`max|PF|` is flat-to-falling across the sweep (0.573 → 0.491) while γ rises 4.25×; if γ were reading the
unrelaxed residual the two would move together. At `k_xb`=3 the signal sits **24.9×** above the residual,
the best separation in the sweep, so this point is the *least* residual-contaminated, not the most.

What does move is the load distribution: **`maxF_cortex` more than halves (5.22 → 2.39) while network
tension quadruples.** Peak node force down, integrated tension up — a softer crossbridge spreads the load
instead of concentrating it, which is the same axis as the existing GATE-A finding that the resting
residual is "concentrated on ~1% of sites, intrinsic to a discrete point load". Stated as a candidate, not
a conclusion: one alternative has been excluded, which is not the same as establishing the mechanism.

**Consequence for R: there is none to quote.** `R = γ(k_xb→0)/γ(reference)` is a LIMIT, and every reading
of it presumes γ descends toward that limit. With a turning point in range the softest point is not an
extreme and dividing by it answers no question. The analyser now refuses on non-monotonicity
(`MONOTONIC_SIGMA = 3.0`, declared with the other two thresholds) — before this it would have reported
**R = 2.83** with a plateau note and no objection, because every other guard it had asked only whether the
points had SETTLED.

**A correction I owe the 2026-07-28 sweep.** That sweep also put `k_xb`=3 above `k_xb`=10, and I attributed
it entirely to non-stationarity. Retiring it was right — those points genuinely had not settled — but the
*explanation* was wrong: the ordering survives when every point is settled. The anomaly was real and I
explained it away.

## `k_xb` = 1 landed: RUNAWAY, not a turning point — and the repo already knew why

γ(`k_xb`=1) = **49.8001 ± 0.0132**, another **4.07×**. The sweep in full:

| `k_xb` | γ_total | γ_network | γ_source | max\|PF\| | maxF_cortex | bound |
|---|---|---|---|---|---|---|
| 1000 | 4.3256 | 4.3030 | 0.02258 | 0.5730 | 6.5203 | 0.9823 |
| 100 | 3.4807 | 3.4538 | 0.02694 | 0.5291 | 6.2784 | 0.9823 |
| 10 | 2.8839 | 2.8069 | — | 0.4999 | 5.2223 | 0.9829 |
| 3 | 12.2425 | 12.1376 | 0.10911 | 0.4914 | 2.3909 | 0.9825 |
| 1 | **49.8001** | 49.7017 | — | 0.3696 | 1.9068 | 0.9823 |

**Two regimes, and the exponent flips sign.** Fitting γ ∝ `k_xb`⁻ⁿ between adjacent points: n = −0.09,
−0.08 across 1000→100→10, then **+1.20, +1.28** across 10→3→1. Every force falls monotonically
(`max|PF|` 0.573→0.370, `maxF_cortex` 6.52→1.91) while γ rises 17×. γ_network is >99% of γ throughout.

### The cause was already written down, in a function whose job is exactly this

`ac/motor/minifilament_topology.working_stroke_strain` calls itself **"the k_xb MASTER-knob arbiter"** and
its docstring says: *the physiological working stroke of 5–20 nm needs `k_xb` ~100–1000 pN/µm; `k_xb` = 1
pN/µm gives a strain larger than the whole minifilament, so the stall geometry collapses.* Evaluated on the
swept doses at `f_stall` = 0.5 pN/head:

| `k_xb` | working stroke | vs 5–20 nm |
|---|---|---|
| 1000 | 0.5 nm | too STIFF |
| **100** | **5.0 nm** | **the only dose inside** |
| 10 | 50 nm | 2× too soft |
| 3 | 167 nm | 8× too soft |
| 1 | **500 nm** | 25× too soft — and larger than the head offset itself (200 nm) |

**The exponent flips exactly where the stroke leaves the window**, and the runaway begins where the strain
exceeds the head's own geometry. So the shape is not a prediction the model is making; it is what running
the model outside its declared domain looks like. **That function existed while this sweep was designed,
queued and run, and nothing — including me — consulted it.** The analyser now does, and refuses.

### A second confound I did not declare

`ac/motor/backbone_warp` sets `k_θ,arm = k_xb · r0_head²` — the lever arm's transverse stiffness is TIED to
`k_xb` (documented, "rigid-lever limit"). So sweeping `k_xb` sweeps the crossbridge spring **and** the
lever arm together. `k_xb` is not a one-parameter knob, and the sweep never said so.

## CONCLUSION — this is not a blebbistatin dose response, and the bound fraction says so

**The bound fraction is 98.2–98.3% at every dose.** I reported that as the negative control passing, and
for *this* sweep it does — the γ change is force transmission, not binding. But blebbistatin's actual
mechanism is to hold myosin in a weakly-bound, low-force state: it lowers the **duty ratio and stall
force**. A sweep that holds the bound fraction fixed to four significant figures has not touched the
mechanism it is named after.

So there are two independent reasons no R may be quoted from this sweep, and neither is about statistics:

1. **Only one dose (`k_xb`=100) is inside the model's validity window**, so there is no valid dose *range*.
2. **`k_xb` is the wrong knob** for the perturbation being modelled — and it is a compound knob at that.

**What survives.** The protocol (run each point to its own stationarity criterion), the instrument (four
guards now declared and tested: parity, plateau, monotonicity, working-stroke validity), the cross-build
parity result, and the measured fact that γ falls only weakly with `k_xb` in the stiff regime (n ≈ −0.08)
while the load redistributes sharply in the soft one. **What does not survive is the experiment's design.**
A blebbistatin analogue needs `f_stall` and the duty ratio swept, at `k_xb` held inside 100–1000 pN/µm —
**a PI question**, since `f_stall` is itself a PI-GAP.

## Cost of reaching stationarity — computed from these runs, and it is flat in `k_xb`

`ROADMAP.md` calls the cost of a native step "the one measurement that unlocks the budget". These five runs
already contain it, so no additional device time was spent:

| `k_xb` | steps | wall | s/step | physical time | real-time factor |
|---|---|---|---|---|---|
| 1000 | 800 | 38.7 min | 2.903 | 8.0 s | 290× |
| **100** | **800** | **38.7 min** | **2.900** | **8.0 s** | **290×** |
| 10 | 1400 | 67.2 min | 2.882 | 14.0 s | 288× |
| 3 | 2000 | 96.0 min | 2.879 | 20.0 s | 288× |
| 1 | 1160 | 55.5 min | 2.872 | 11.6 s | 287× |

**s/step is 2.87–2.90 across a 1000× range of `k_xb`** — the per-step cost does not depend on the swept
parameter at all, so the cost of reaching stationarity is set *entirely* by how many steps settling takes.
That is the useful form for a sweep budget: cost = (settling steps) × 2.9 s, and only the first factor
varies. At the one physiological dose it is **38.7 min of A5000 for 8.0 s of physical time.**

⚠ **The caveat travels with the number**: every step here is FORCE-ACCEPTED. This is the cost of reaching
a stationary *force-accepted* trajectory, which is **not** the "converged native step" stage 2 asks for.
Quoting it as stage 2's measurement would be exactly the promotion this repo forbids.

## Still unchecked: are these γ values converged in `dt`?

Every run above used `dt = 0.01`. Nothing has tested whether the settled γ depends on it. A `dt = 0.005`
run at `k_xb`=100 — the only physiological dose — is queued, with the prediction declared first:
**γ(dt=0.005) must equal 3.4807 within 3 combined sems.** If it does not, `dt=0.01` is not time-converged
and **every settled γ in this sweep inherits that error.**

`dt` is swept *downward* deliberately: `k_on·dt` = 0.5 at `dt`=0.01 is already coarse, and `dt`=0.02 is
recorded as breaking (2026-07-28g), so the safe direction is the only informative one.

### RESULT — REFUTED. `dt = 0.01` is NOT time-converged, and this is the largest finding of the sweep

**γ(dt=0.005) = 3.1998 ± 0.0011 against γ(dt=0.01) = 3.4807 ± 0.0017 — a 137 σ, 8.07% difference.** Both
points are STATIONARY over comparable physical time (8.8 s vs 8.0 s), so this is not a duration artifact.

**The shift is not in the kinetics.** Decomposed:

| observable | dt=0.010 | dt=0.005 | relative | σ |
|---|---|---|---|---|
| **γ_network** | 3.45379 | 3.17296 | **8.13%** | **135.7** |
| `max_f_cortex` | 6.2784 | 5.9263 | 5.61% | — |
| γ_source | 0.02694 | 0.02663 | 1.16% | 1.2 |
| `bound_fraction` | 0.98229 | 0.98472 | **0.25%** | 4.1 |

The binding kinetics are converged to 0.25% and the crossbridge stress term to ~1%. **The entire dt
dependence sits in the mechanical relaxation term.**

**What this costs.** Every γ in this REPORT was measured at `dt = 0.01` and therefore carries an ~8% dt
error at the one dose where it has been measured — direction and size at other `k_xb` unknown. Crucially
this is **not** covered by the "R is a ratio, so it survives the per-parameter PI-GAPs" argument used
earlier: a ratio only survives an error that CANCELS, and nothing shows the dt error is the same at two
different `k_xb`. **The lane cannot quote γ ratios either until dt convergence is established.**

**What it does NOT touch.** The A/B equivalence and the cross-build parity are comparisons at the *same*
`dt`, so both stand exactly as reported. This is a statement about absolute values, not about whether two
code paths agree.

**And it moves stage 2 further away, not closer.** `ROADMAP.md` stage 2 asks for `dt` to be **RAISED** and
the answer shown `dt`-independent. `dt` cannot be raised from a point that is not converged.

### Two readings remain, and the run that separates them is queued

**(A) Relaxation depth** — halving `dt` doubles inner iterations per physical second, so less unrelaxed
stress survives and γ falls. Under this reading the tension this lane reports contains unrelaxed residual,
which is the same failure mode as `STATE.md` (c) 3.
**(B) Integrator order** — ordinary first-order discretisation error, with no residual interpretation.

Both predict "γ moves with `dt`", so the sweep cannot distinguish them. A run at **`dt`=0.01 with `--outer`
40 → 80** matches (A)'s relaxation budget *without* changing `dt`, and is queued with the reading declared
first: **γ falling toward 3.20 ⇒ (A); γ staying near 3.4807 ⇒ (B).**

### SEPARATED — it is (B), the outer time integrator. The inner solve is not the problem

| run | `dt` | `--outer` | inner iters | wall | γ_total |
|---|---|---|---|---|---|
| baseline | 0.01 | 40 | 32,000 | 38.7 min | 3.4807438591 |
| `dt` halved | 0.005 | 40 | 70,400 | 68.9 min | **3.1998352984** |
| **outer doubled** | 0.01 | **80** | **64,000** | 61.5 min | **3.4807439934** |

Doubling the relaxation budget moved γ by **6.0e-08 relative**. Halving `dt` moved it by **8.07%**. The
flag is plumbed and the work was done — inner iterations 32,000 → 64,000, wall ×1.59 — and `max|PF|` did
improve, 0.5291 → 0.4998 (5.5%). **So the residual falls with more iterations while γ does not move at
all.** `bound_fraction` and `n_bound` are bit-identical across the pair, as expected.

Two things follow, and the second is the useful one:

1. **γ is not reading unrelaxed residual.** A 2× relaxation budget that visibly lowers `max|PF|` leaves γ
   unchanged to 7 decimal places. This is direct evidence for the conclusion I had softened above — the
   `max|PF|`-based argument was weak, but the answer it reached was right.
2. **The `dt` dependence is in the OUTER physical-time integration**, not the inner mechanical solve. That
   matters for where effort goes: this lane's convergence work has been inner-solver conditioning
   (multigrid, preconditioners, IMEX step scale), and none of it addresses this.

### Order still unknown — the third point is queued with both predictions declared

Two `dt` values cannot give a convergence order. Extrapolating from the two:

| assumed order | γ(`dt`→0) | predicts γ(`dt`=0.0025) |
|---|---|---|
| first | 2.9189 | **3.0594** |
| second | 3.1062 | **3.1296** |

Separation **0.070 pN/µm** against sems ~0.0015 — a ~47 σ discriminator. `dt`=0.0025 is queued.

### RESULT — both predictions refuted, and the detector was wrong about its own window

Measured γ(`dt`=0.0025) = **3.4148**, against 3.0594 (first order) and 3.1296 (second). Neither fits.

**But the run first had to be re-measured, because the detector mis-set its own transient.**
`equilibration_point` chose `start = 0` on a series that rose 0.8 → 3.4 in 0.5 s and was then flat for 7 s,
so the "steady-state mean" was taken ACROSS the rise: **3.373 with std 0.2538**, against the plateau's
**3.415 with std 0.0106**. The figure shows a clean plateau; the number did not.

It also hid itself. The drift verdict is `|drift| < σ·std`, and keeping the transient inflated `std`
twentyfold — **a larger transient makes the drift test easier to pass.** The run was certified STATIONARY
at 0.94 σ. Fixed: `window_opens_on_a_transient` compares the retained window's leading tenth to its
trailing tenth, normalised by the **trailing** tenth's scatter (transient-free, so it cannot inherit the
same loop), and returns DRIFTING with no mean when they differ by more than one standard deviation. It
reads **0.0 / 0.8 / 33.5** on the three runs — the two good cuts keep, the contaminated one rejects.

**Re-measured on an identical window** (`t ≥ 2.0 s` for all three, since each auto-selected a different
transient — 0.40 / 1.30 / 0.00 s — so the comparison had not been like-for-like):

| `dt` | auto-cut mean | **identical-window mean** | std |
|---|---|---|---|
| 0.01 | 3.48074 | **3.48131** | 0.0113 |
| 0.005 | 3.19984 | **3.20026** | 0.0104 |
| 0.0025 | 3.37312 | **3.41477** | 0.0106 |

**The sequence is NOT MONOTONE: 3.4813 → 3.2003 → 3.4148.** It falls, then rises. No convergence order
fits, and Richardson extrapolation does not apply.

⚠ **RETRACTED — "either way the `dt`=0.01 value is 10–19% high."** That rested on assuming first or second
order; both are refuted. The true `dt`→0 value is **unknown**. This does not weaken (c) 17, it strengthens
it: an error of unknown size and unknown sign cannot be corrected for.

### The control that should have come first, and has not been run: SEED

**Every run in this REPORT used `seed = 0`.** The seed-to-seed scatter of the settled γ has never been
measured — yet every comparison here presumes it is small next to the differences being read. The `dt`
sequence spans **0.28** pN/µm; if replicates at fixed configuration scatter by that much, the sequence is
noise and there is no `dt` story to explain.

Two replicates at `k_xb`=100, `dt`=0.01 (`seed` 1 and 2) are queued, with the reading declared first:

- **agree within a few combined sems** → seed scatter is negligible, the `dt` differences are real, and
  the non-monotonicity above is a genuine property to explain;
- **scatter by ~0.2** → the `dt` sequence is seed noise, and (c) 17's premise must be re-examined.

Either way the `k_xb` 10 → 3 jump (2.88 → 12.24, a factor 4.25) survives, being ~40× larger than the
scatter in question. This is a control on the `dt` claim, not on the sweep's headline.

#### First replicate: the `dt` differences ARE real — and the sem was the wrong yardstick all night

| | identical-window mean (`t ≥ 2.0 s`) | std |
|---|---|---|
| `seed 0` | 3.48131 | 0.01132 |
| `seed 1` | **3.46411** | 0.01052 |

Difference **0.0172 pN/µm** — **0.49%** of γ, and **6.1%** of the `dt` spread. That is the first branch
declared above: the seed scatter is an order of magnitude below the `dt` effect, so the `dt` dependence is
real and (c) 17's premise holds. The `k_xb` 10 → 3 jump is ~550× the scatter and was never in doubt.

**But a third fact falls out, and it corrects how tonight's numbers were expressed.** 0.0172 is **13 naive
combined sems**. The within-run `sem` (~0.0011) therefore *overstates* run-to-run reproducibility by an
order of magnitude: re-running the identical configuration moves γ by 13× the error bar the run reports
for itself. **The right yardstick for comparing two γ values is the seed scatter, not the sem** — so every
σ figure quoted tonight ("137 σ", "168 σ", "135.7 σ") is inflated, because all of them are sem-based. The
orderings they support are unchanged (the `dt` effect is still ~16× the scatter); the *significances* are
not what was written. The second replicate is running so the scatter can be estimated rather than inferred
from a single pair, after which those figures will be restated against it.

#### Three replicates: the scatter is 1.95% of γ, and it retires one of tonight's conclusions

| seed | identical-window mean | within-run sem |
|---|---|---|
| 0 | 3.48131 | 0.00103 |
| 1 | 3.46411 | 0.00087 |
| 2 | **3.35797** | 0.00103 |

**S = 0.0668 pN/µm = 1.95% of γ — 67× the within-run sem.** Two single-seed values differ with uncertainty
`√2·S = 0.094`. Every comparison made tonight, restated on that basis:

| comparison | \|Δ\| | proper σ | what I said | verdict |
|---|---|---|---|---|
| `k_xb` 10 vs 3 | 9.3586 | **99.1** | — | established |
| `k_xb` 1000 vs 100 | 0.8411 | **8.9** | — | established |
| `k_xb` 100 vs 10 | 0.5968 | **6.3** | ~~168 σ~~ | established |
| `dt` 0.01 vs 0.005 | 0.2811 | **3.0** | ~~137 σ~~ | marginal |
| `dt` 0.005 vs 0.0025 | 0.2145 | **2.3** | — | **NOT established** |

**Survives, and strengthens**: the `k_xb` runaway at 99 σ, and the working-stroke argument, which is a
statement about the model's domain and not a statistical one at all.

**Weakens**: `dt` non-convergence is real but **3.0 σ, not 137**.

⚠ **RETRACTED — the `dt` non-monotonicity.** At **2.3 σ** it is not separable from seed noise. With one seed
per `dt` point, "non-monotone `dt` dependence" and "monotone `dt` dependence plus scatter" cannot be told
apart by these data, so "no convergence order fits" and "Richardson does not apply" are withdrawn as
unsupported. I read a 2.3 σ wiggle as a structural finding and built two paragraphs on it.

**The methodological consequence outlives the sweep**: a single-seed run cannot resolve a γ difference below
~5%. That is a standing limit on every sweep this project runs — and stage 4's whole business is separating
parameters by small changes in an observable, so **1.95% sets a floor on identifiability** before any Fisher
argument is made.


⚠ This also reopens a judgement I made earlier tonight. I ruled the residual explanation "DISFAVOURED" for
the `k_xb` runaway because `max|PF|` stayed flat while γ rose. But `max|PF|` is a **max over nodes** and
γ_network is an **integral** — the max can hold steady while integrated unrelaxed stress moves. That
argument was weaker than I presented it.

## Figures

- `kxb_*/figs/stationarity_trajectory.png` — γ_total vs physical time; shaded transient, mean ± sem band,
  and the acceptance contract printed on the figure so the verdict is judged against its own criterion.
- `kxb_*/figs/stationarity_source_network_split.png` — source and network terms on separate axes with the
  closure residual; this is where the `γ_source` drift at `k_xb=10` is visible.
- `kxb_*/figs/stationarity_autocorrelation.png` — ACF against the estimator's sampling floor, which is what
  makes τ_int an upper bound rather than a value at the two stiff points.
