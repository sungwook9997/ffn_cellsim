# ROADMAP — the destination and its stage gates

**This is a STATUS file, not a charter.** `CLAUDE.md` holds what does not change; everything here does.
It is moved out of the boot context on purpose: it was 46% of `CLAUDE.md` and most of it was current
state wearing the clothes of a rule.

## ⚠ What in this file changes — check before quoting any of it

| What | Changes when | Verify against |
|---|---|---|
| **Which stage is current, and which gates have passed** | every session that lands a run | `STATE.md` (b) and (c) — this file states the gate, never the score |
| **The gate WORDING itself** | **a PROPOSED amendment is pending PI ratification** — see the ⚠ block below the table. Until it is signed, the ratified table stands | the amendment block's own header |
| **K and M — the cell types and states** | a PI decision | `CellState` in the Notion Contract-Graph (the 9th DB), which is the machine-readable copy |
| **The stage-0 observable set** | when a measurement contract is signed | the `ValidationGate` rows named in each entry |
| **Every number in the budget arithmetic** (cost per step, sweep sizes, GPU-days) | any engine or solver change; several were re-measured the day after they were written | `STATE.md` (f) and the run records under `outputs/ac/` |
| **The five corrections** | closed — kept because each was a real reversal | — |

**Rule of thumb:** this file says *what a stage must prove*. It never says *whether it has been proved* —
that is `STATE.md` (b), and if the two ever disagree, `STATE.md` wins and this file is the defect.

---

> ✅ **RATIFIED PI 2026-07-28.** Drafted from the PI's own six-stage outline plus five corrections argued
> from measured repo facts, then closed out with the three decisions that were blocking it: the cell types
> and states (K, M), the stage-0 observable set, and the finding that a brute-force sweep is impossible by
> three to five orders of magnitude. All three are recorded below with their evidence.

**Destination.** A platform that infers **per-cell-type, per-cell-state molecular parameters with their
uncertainty** from fine-grained mechanistic simulation (PI reframe 2026-07-25). Parameters are OUTPUTS;
literature supplies **priors**, never fixed point values. The deliverable is therefore stage 5 — a cell
library — and every earlier stage exists to make that stage's comparison meaningful. **Stage 5 also defines
the denominator for "how far along are we"**: `cell types × states × parameters recovered`. Without it that
question has no answer, which is why it had none before this section existed.

**K = 2, M = 2 (PI 2026-07-28).** The library is **MCF7** and **MDA-MB-231**, each in **suspended** and
**adherent on 2D collagen**. Not a preference — the KB is already deepest on exactly these two (32 vs 22
SourceEvidence rows, and `VG-231-G0…G6` drafted), and they are *the* contrast pair in breast-cancer
mechanobiology: same tissue origin, opposite motility phenotype. That matters because **one cell type cannot
validate an inference platform** — you can always fit one. Two is the minimum that makes the claim
falsifiable: the platform must recover the KNOWN difference between them *without being told it*. A third
line triples cost and adds no falsifiability; NIH/3T3 already appears inside `VG-231-G5` as a cross-check,
not as a library entry. **Migration in 3D collagen is deliberately NOT a state** — it sits behind T10 and
`BLOCKS_CRAWL`, so including it now would keep the metric red for a reason unrelated to inference. It is a
stage-6 branch.

**Progress denominator: `K × M × (identifiable parameters, MEASURED)`.** The parameter count is deliberately
not fixed. 33 PI-GAP cards exist, but an observable set does not constrain 33 parameters — and measuring how
many it does constrain is exactly what the Fisher information in `ac/engine/observe/sensitivity.py` is for.
Leaving the third factor measured rather than declared makes the metric self-correcting: claiming 33 when
the data constrains 5 is then arithmetically visible rather than a matter of opinion.

**Stage-0 observable set (6).** Each is reported by experiment, computable by this engine, and already has a
KB band. *Suspended*: cortical tension γ · intracellular pressure ΔP · AFM apparent modulus E + cytoplasm
relaxation τ (`VG-γ-MCF7-suspended`, `VG-231-G0`, `VG-231-G1`). *Adherent-2D*: total traction |F| vs
substrate stiffness · spreading area A (top-down silhouette, the PI-ratified metric) · 2D migration speed +
persistence + aspect ratio (`VG-231-G2`, `VG-231-G3`). Several are currently BLOCKED — γ sits in `STATE.md`
(c) and native crawl is 0.12 nm/s — which is fine: **stage 0 is about fixing protocol, prior and pass
condition BEFORE a run, not about passing today.**

**A brute-force sweep cannot build the library, and the arithmetic is not close.** Measured: the only native
parameter sweep that exists (`outputs/ac/resting_native/kxb_sensitivity.json`) is **2 points at ~20 min each
on 70,686 filaments — and its residual did not move** (5.789 → 5.789), so it measured nothing. At that cost a
5-parameter × 5-level grid is 43 days, an 8 × 5 grid is **14.9 years**, and both are ×4 for the cell states
and ×~150 again for a step that actually converges (172 s vs 1.1 s measured). **The bottleneck is one run,
not the sweep.** Four multipliers make it affordable, and every one is already a stage: raising `dt` under a
physical predicate (stage 2, largest and still unmeasured), cutting dimension by Fisher (~125×), replacing
the grid with the ALREADY-VALIDATED analytic gradient (~20–60×), and a surrogate with sequential design.
Together: ~390,000 runs → ~100–300, i.e. ~17 GPU-days for all four cell states. **This is why stage 4 is
"sensitivity-directed" and why stage 2 precedes it.**

**The one measurement that unlocks the budget, and it is cheap.** Nobody knows what a CONVERGED native run
costs, because no converged native run exists — 172 s/step is a converging *attempt*. Every native run must
now fill `observe/artifact.py::timing_block`, whose own rule is the point: a timing is a property of an
engine **at a stated convergence**, since any solver is arbitrarily fast if allowed to stop early. One
GPU-day turns the sweep budget from guesswork into arithmetic.

**Measured 2026-07-29 — the cost half is done, and it is simpler than the plan assumed.** Five settled
full-native cortex-motor runs give **2.87–2.90 s/step across a 1000× range of `k_xb`**: the per-step cost
does not depend on the swept parameter at all. So the cost of reaching stationarity is
`(settling steps) × 2.9 s` and **only the first factor varies** — which is the form a sweep budget needs.
At the one physiological dose it is **38.7 min of A5000 for 8.0 s of physical time** (RTF ≈ 290×).
⚠ Every step in those runs is FORCE-ACCEPTED, so this is the cost of reaching a stationary *force-accepted*
trajectory and **not** the "converged native step" this stage asks for. Quoting it as stage 2's measurement
would be the promotion the charter forbids.

**And the `dt` half moved AWAY, not closer.** Stage 2 asks for `dt` to be **RAISED** and the answer shown
`dt`-independent. Measured at the physiological dose: γ(`dt`=0.005) is **8.07% below** γ(`dt`=0.01) at
137 σ, so `dt` = 0.01 is not converged — and `dt` cannot be raised from a point that is not converged. A
third point (`dt`=0.0025) refuted both first- and second-order predictions and left the sequence
**non-monotone** (3.4813 → 3.2003 → 3.4148), so no convergence order fits and the error's size is unknown.
Doubling the inner relaxation budget moves γ by 6e-08, so the dependence is in the **outer** time
integration, not the inner solve — which matters because this lane's convergence work to date has all been
inner-solver conditioning. Detail and the retraction that came with it: `STATE.md` (c) 17.
⚠ **One control is still missing and it undercuts the above if it fails**: every run used `seed = 0`, so
the seed-to-seed scatter has never been measured. Replicates are running with the reading declared first.

| # | Stage | Exit gate — falsifiable, declared before the run | Prereq |
|---|---|---|---|
| **0** | Measurement contract | every observable's protocol, literature prior and pass condition exist as a signed gate contract BEFORE any run | — |
| **1** | Assembly = transmission | one lane's coupled solve **converges** and acceptance gates on `balance_ok_d`; connector dispatch count is recorded and each dispatched connector's two-sided resultant closes as an adjoint pair | — |
| **2** | Time-scale **and cost** | under that physical predicate `dt` is **raised** and the answer shown `dt`-independent, AND the cost of a converged native step is MEASURED via `timing_block` — this is what opens stage 4's budget | runs with 1 |
| **T10** | Exterior solve | Stokes/Brinkman lands **and `γ_node = 6πηR/Nc` is removed in the same change**; the six whole-cell rigid modes are set by physics, not by numerical regularisation. Releases `BLOCKS_CRAWL` | 1 |
| **3-1** | State: **suspended** | γ, ΔP, E, τ at native FULL population with a converged solve, inside stage-0's band | 0,1,2 |
| **3-2** | State: **adherent on 2D collagen** | traction \|F\|, spreading area A, and 2D migration speed/persistence/aspect inside their bands. Monomer/G-actin flux is a PREREQ, not a stage-6 branch — without treadmilling a protrusion does not persist | T10, 3-1 |
| **4** | Sensitivity-directed inference | only the axes native Fisher says are identifiable, driven by the analytic gradient + a surrogate rather than a grid; every point passes the stage-1 convergence gate | 2, 3-1 |
| **5** | **Cell library (the product)** | for **2 cell types × 2 states**, parameters recovered **with uncertainty** and checked against held-out experiments not used to fit | 4 |
| **6** | Branches | migration/invasion in 3D collagen (needs T10 + MMP) · 3D spreading · intravasation · a third cell line if ever — each with its own stage-0 contract | 5 |

### ⚠️ PROPOSED AMENDMENT 2026-07-28 (night) — the gates above score a STATIC object that does not exist

**Not yet ratified. The table above stands until the PI signs this.** Raised because the premise it rests
on was measured false, and the measurement is a year-old habit rather than a new result.

**What was measured.** `GATE_A_RESTING_CONVERGENCE_2026-07-23.md` §10: the resting cell with myosin **OFF**
converges seed-robustly (`max|PF|` 0.1398 / 0.1416 / 0.1586, ≥24% under the floor). Seed 4,420 discrete
isometric heads and the residual is ~1.5 and **no solver removes it** — `PF ≈ f_head`, ∝ `f_head`, inherent
to discrete point loads, needing ~7× more heads (unphysical) to fall under the gate. Today's native
`resting_cost_curve` re-measured the same wall: STALLED, 6.30× the derived floor, 6.12% of the residual
removed per ~2,000 iterations.

**Why.** A living cortex's tension is maintained **dynamically** by myosin turnover — a non-equilibrium
steady state. The force field is measurably **non-conservative** (`∮F·dx ≠ 0` with one bound head; (b)), so
no potential exists to minimise, and the solver's **symmetric tangent silently assumes one** — the
assumption lived in the linear algebra, not in anyone's head, which is why it was never checked. Asking
that system for a static equilibrium is a category error, and the irreducible residual is the error message.

**What this changes — gates, not stages.** The stage ladder survives; what each gate SCORES does not:

| Stage | Ratified gate | Proposed |
|---|---|---|
| **1** | "one lane's coupled solve **converges**" | inner (chemistry frozen): ~~accept on `balance_ok_d`~~ — **⚠ this line is DEFECTIVE and must be rewritten before signature.** `sf_motor_slice.py:65` states that predicate "tests the adjoint wiring, **not the convergence**", so as written the amendment rests a convergence gate on something that cannot fail for non-convergence. The two must be named separately. Outer: **stationarity** — observable drift below its own fluctuation over ≥N bound-head lifetimes |
| **2** | cost of a **converged** native step | cost of **reaching stationarity**; needs the turnover time, so it is defined only after 0.5 |
| **3-1/3-2** | γ, ΔP, E, τ inside the band | the same observables as **time averages with fluctuation bands** — which is what the experiments do too (Wang 2021 moves γ 32→50 mN/m on suction pressure alone) |
| **NEW 0.5** | — | **measure the timescale separation**: mechanical relaxation vs bound-head lifetime (`k_off0` = 0.4/s → 2.5 s, Nagy 2013). If they are comparable, even quasi-static-at-fixed-chemistry fails |

**The static turgor baseline is not discarded — it is relocated to the INITIAL CONDITION.** Start from the
converged passive equilibrium, switch myosin on, integrate to stationarity. That is what this file's own
physiological-baseline rule already says ("measure on the turgor-pressurised cell, then add myosin as a
modulator"); the runtime was doing the opposite.

**Status 2026-07-29 (overnight), so the amendment is judged against what now exists rather than what it
assumed.** Three of its premises have moved:

- **The stationarity detector is built, wired and tested.** It runs in the cortex-motor driver, judges each
  observable, and its verdicts are RECORDED rather than scored — exactly the state a pending amendment should
  leave them in. The "nothing in the repo has one" line below is now historical.
- **Stage 0.5 has its first measurement.** At full native, γ decorrelates at least 54× FASTER than the
  bound-head lifetime — the OPPOSITE of the detector's own assumption that collective relaxation is slower than
  its driver. And τ came out at the estimator's sampling floor, so it is an upper bound: the honest statement is
  "uncorrelated at 0.05 s", not a value. That makes the driving-time floor ~2.4× more conservative than the
  statistics require, which is now a known margin rather than an unexamined one.
- **The driving time itself was wrong in the code by 1.75×** (a Bell prefactor where the run used catch-slip,
  whose zero-load rate is the SUM of two pathways). Stage 2 is denominated in that number, so it could not have
  been evaluated correctly before this was fixed.
- **Stage 1 gained a mechanism it did not have: a component can now be handed one channel of the force
  assembly** (PI decision `COMPARTMENT_VALIDATION_TRACKS` §6 D2, approved 2026-07-29 — name the
document; `AC_EXECUTION_PLAN` §9 has a different D2). Before it, binding an engine component into a production
  run would have added its force TWICE, so "assembly = transmission" could only be tested on a slice rig. It is
  a mechanism, not evidence — the native A/B and the double-count gate are queued, and stage 1's gate is
  unchanged.

**What must still be built — narrowed 2026-07-29, because the original sentence has been half-answered.** The
amendment said "a stationarity detector must be built; nothing in the repo has one". The detector now exists
(`ac/engine/observe/stationarity.py`), runs in the production driver, and every sweep point since carries its
own verdict. What is NOT done is the part the amendment actually turns on:

- **The verdict is recorded, not scored.** No gate reads it. That is deliberate while the amendment is
  unsigned — a detector that gated runs before the PI ratified the gate wording would be a lane scoring itself
  — but it means the amendment is not implemented, only *implementable*.
- **`SOLVE_COUPLED` is still empty** ((c) 4). Stationarity is measured on an observable trajectory the
  incumbent driver produces; nothing in `ac/engine` drives a coupled solve to stationarity yet.
- **The outer gate has no declared N.** "≥N bound-head lifetimes" — the sweep runs `min_windows=5` and the
  choice is declared before each run, which is honest but is a driver flag, not a contract.

**Scale of the miss, for calibration:** GATE B's headline tension emergence ((c) 3) ran 30 steps at
`dt_phys` 0.01 = **0.3 s of physical time against a 2.5 s bound-head lifetime — 12% of one turnover.** It
observed a transient and no stationarity test existed to say so.

**Five corrections to the outline, with their evidence** — kept because each was a real reversal:

1. "It runs" is **already true and means nothing**: the full-native cell builds and steps while recording
   `inner_converged: false` and `overall_status: INCOMPLETE_OR_FAIL`. Stage 1's gate is convergence.
2. `dt` goes **up, not down**. Shrinking `dt` is the explicit path's remedy, and the solver-side optimisation
   ceiling was already measured at ~9% ("do not retry"). The lever is step size under a physical predicate.
3. Crawl is **not simply after** the static test — it is after T10, and trap #4 (`γ_node`) makes a partial
   T10 exactly a 2× error.
4. Stage 0 did not exist in the outline. The AFM "too stiff" result was **three protocol artifacts**, not
   physics; without a protocol fixed in advance, stage 3-1 measures artifacts again.
5. A naive sweep is both unaffordable and meaningless on an unconverged observable. Sensitivity ranking is
   what makes stage 4 finite.

**What is NOT a gate, ever** — each of these has passed while the physics was wrong: a run completing; a green
test suite (`structural tests passing ≠ production`); a conclusion drawn on a slice (**develop on a slice,
conclude at native** — 2026-07-28: four headline numbers measured at 0.18% of native were retired one day
after landing, see `STATE.md` (c) 15); and a magnitude quoted without its `QuantitativeClaim` axis.

**Where the detail lives.** Per-stage work: `aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md` (read its
§10 for the current ordering) and `COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md`. Current status against these
gates: **`STATE.md`** — that file, not this one, says where we actually are.

**Historical, closed, not current state** (all four pre-date the Warp-only contract and the ×40 retirement —
`STATE.md` (d) lists them as the documents most likely to mislead a fresh reader): `PHASE_0_CLOSEOUT.md`,
`PHASE_0_3_DECISIONS.md`, `AFINES_ALGORITHM_NOTES.md`, `docs/briefs/H{1,2,3,4}_*.md`. The H.1→H.7 chain they
describe has had no work since 2026-06-11 and h7/h8/h9 have no closeout.

