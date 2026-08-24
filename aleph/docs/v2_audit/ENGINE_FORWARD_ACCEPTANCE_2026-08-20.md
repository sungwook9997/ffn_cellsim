# Convergence was never the hard question. It was the malformed one.

**Status:** ANALYSIS + PLAN. Decides nothing; selects no observable, band, window or constant.
Written 2026-08-20 by session `d380041d` at PI instruction: *"수렴의 의미는 말이 안된다고 생각함. 수렴을
애초에 그리고 보는데 기준들이 너무 많고 적음."*

Supersedes nothing. It **extends** `PI_DECISION_e1_STATIONARITY_2026-08-16.md`, which established that the
stationarity amendment is undecidable as written. That document asked the PI for three decisions. This one
argues that a decision was not the missing thing.

---

## 1. The two sentences, made exact

**"수렴의 의미는 말이 안된다."** Convergence is a property of *a solver applied to a problem that has a
solution*. Three independent lines of evidence say this problem has none. Then "which solver converges?"
is not a hard question — it is a question with no referent, and we asked it five times.

**"기준들이 너무 많고 적음."** Measured on this tree today: **17 distinct verdict strings** in engine code,
an **8-rung evidence ladder**, **3 evidence tiers** in `STATE.md`, and at least **7 ratchets**
(`DOCS_DEBT`, `REPORT_DEBT`, `UNFINGERPRINTED_ENTRIES`, `CUDA_ORDINAL_DEBT`, run coverage, params coverage, the oracle-breach counter). And
**exactly one predicate decides whether a physical step is accepted** — `ledger.py:250` — and it tests
floating-point roundoff.

Both sentences are the same observation from opposite ends.

---

## 2. Why there is no target — three lines, none of them a solver property

**(a) The driver says so, in writing.** `driver.py:260-261`, in `_preload_erm_resting_balance` (`driver.py:254`; the e1 document cites `233-241`, which the 2026-08-19 merges shifted):

> *"a relaxed cortex plus a positive turgor is not a valid resting baseline — there is no force-balanced
> equilibrium there, which is why no inner solver converges from it."*

Of 23 recorded runs carrying a resting residual, **none applied either preload**.

**(b) The active tension source cannot run, by construction.** Two constants on the resting path are
`PI-GAP` with **no default** (`driver.py:1592, 1597`): the fraction of NMII heads bound at rest, and the
resting isometric per-head tension. Without them there is no active cortical hoop tension. Turgor then has
nothing to balance against, and `myosin_bound_after_step` reads **0** in the committed full-native record.

This is the causal chain, and every link is on disk:

    two PI-GAP constants have no value
      -> no active cortical tension is generated at rest
        -> the turgor load has no opposing hoop stress
          -> no force-balanced configuration exists
            -> no solver descends, and five did not
              -> and the predicate that accepts steps tests roundoff, so nothing ever said so

**(c) The empirical signature is not "under-resolved".** Explicit relaxation is **non-monotone at two
population scales five orders of magnitude apart** (156 DOF and 1.65 M DOF). A 66× inner budget ends
**higher** than it started (60 iterations −0.52%, 4,000 **+0.14%**). More iterations make it worse. A solver
short of budget gets slowly better; a solver descending toward nothing does this.

⚠ Five candidate causes were eliminated on 2026-08-15 — non-conservativity, cross-link pre-stress, the
iteration budget, solid–fluid coupling, the infeasible shell. **Those eliminations stand.** What does not
stand is the conclusion drawn from them, *"the operator is the only cause left"*: it presupposes a target.
The e1 document already retracted it to PROVISIONAL. This document is the second reason.

---

## 3. Why 23 runs did not notice

The acceptance predicate. `aleph/engine/ledger.py:239` (`_assemble_force_balance_kernel`, the predicate at `:250-259`):

```python
r = reaction_resultant[0] + traction_resultant[0]
n2 = wp.dot(r, r)
balance_ok[0] = 1 if n2 <= tol_sq[0] else 0
```

and `tol_sq` comes from `assemble_balance_tolerance_ratio(n_terms)` — the **Higham summation bound**, i.e.
`n_terms × eps64`. So the gate asks: *do the two independently sourced force channels cancel to within
float64 roundoff?* That is Newton's third law holding **in the wiring**. It is true at any configuration
where the connectors are correctly attached — resting, mid-transient, or absurd.

`STATE.md` already records this of the SF-motor lane: *"the predicate tests ADJOINT CLOSURE, not
convergence... so a non-converged step is not what it rejects."* The same predicate is what
`composed_native.py`, `cortex_motor_slice.py`, `erm_cortex_slice.py` and `interior_column_slice.py` pass as
`accepted_d`.

`inner_converged` **is** computed. It is recorded in the run record. **No gate reads it.**

So the engine has held, simultaneously, a gate that cannot fail and a quantity that cannot succeed. Neither
is a bug in isolation. Together they are a control loop with the feedback wire cut.

---

## 4. The census — what "too many and too few" measures

| kind of criterion | count | what it decides |
|---|---:|---|
| verdict strings in engine code | 17 | how to *label* a result |
| evidence-ladder rungs | 8 | how far a *component* has come |
| `STATE.md` evidence tiers | 3 | whether a *claim* may be quoted |
| ratchets | ≥7 | whether *debt* grew |
| `run-record@2` stamps | — | whether a *record* is self-describing |
| **predicates gating a physical step** | **1** | floating-point roundoff |
| **criteria stating what the resting cell IS** | **0** | — |

The asymmetry is the finding. **This project has a rigorous epistemology of claims and no definition of the
state.** Every one of those mechanisms answers *"may we say this?"* Not one answers *"is the cell right?"*
That is not an oversight of effort — it is where all the effort went, and it is why the machinery is
excellent at refusing to overclaim and silent about being wrong.

---

## 5. The surprise: most of what is needed is built, and unwired

This changes the cost of everything below, so it is stated before the plan.

| piece | state | gap |
|---|---|---|
| stationarity detector | **`aleph/observe/stationarity.py`, 1,030 lines** — Sokal automatic windowing, `tau_int`, `n_eff >= 25`, an explicit transient guard, 5 verdicts. Runs in the production driver; every sweep point carries one. | **no gate reads the verdict** |
| cortical tension γ | **`cortical_tension.py::measure_cortical_stress`** — method of planes, `gamma_source` / `gamma_network` separated, Newton-closure diagnostics kept distinct from γ | **not on the resting path** |
| the balance MECHANISM | **`resting_balance_oracle.py`** already composes turgor + membrane area tension + unilateral ERM + discrete contractile cortex dipoles, and reads γ back by method-of-planes *"never imposed"* | proven on a **coarse** mesh; §4 of the diagnosis says the native solver cannot equilibrate it |
| seed-scatter rule | measured: **S = 1.95% of γ, 67× the `sem`** | recorded in `STATE.md` (f); nothing enforces it |

⚠ **A defect found while checking this, then MEASURED rather than left as a reading of the source.**
There are two different `stationarity.py`: `aleph/observe/stationarity.py` (1,030 lines, imported by the
`observe/` stack — `fluctuation`, `kinematics`, `tether`) and `aleph/engine/observe/stationarity.py`
(472 lines, imported by the production GATE-B driver `ac_gate_b_cortex_motor_native.py:114`).

**The statistic is identical.** `integrated_autocorrelation_time` agrees to `ratio = 1.00000` on white
noise, AR(1) ρ=0.9, AR(1) ρ=0.99 and a drifting series. Both use Sokal `c = 5.0`, both carry the same
three verdicts. So my first reading — *"which one is in force decides which criterion"* — was **wrong about
the statistic** and is corrected here.

**The VERDICT contracts differ, and they disagree at a measured rate:**

| default | engine/observe (production driver) | observe (observation stack) |
|---|---:|---:|
| `min_windows` / `min_tau_windows` | **5.0** | **50.0** — 10× longer window required |
| `drift_sigma` | **1.0** | **3.0** — 3× more drift tolerated |
| `driving_correlation_time_s` | **required argument** | absent |

Swept over 420 AR(1) series (ρ ∈ {0.5, 0.9, 0.99} × 7 drift levels × 20 seeds, N = 4,000, dt = 1e-3):

    disagreement                                        129 / 420 = 30.7%
      engine STATIONARY  while  observe DRIFTING         74
      engine DRIFTING    while  observe TOO_SHORT        31
      engine STATIONARY  while  observe TOO_SHORT        16
      engine DRIFTING    while  observe STATIONARY        8

⚠ **Both directions occur**, so neither module is the conservative one and this cannot be settled by
picking it. And the dominant mode — 74 of the 129 — is the **production driver calling STATIONARY what
the observation stack calls DRIFTING**. That is the direction that matters: it is the permissive one, and
it is the one a stationarity gate would be built on today.

⚠ **This paragraph was wrong twice before it was measured.** The first version claimed the two modules
implement different statistics (they do not — `tau_int` agrees to 1.00000). The second quoted a single
series and got the disagreement direction backwards. Both errors came from reading the source instead of
running it, on a question where running it costs seconds. Pinned by
`aleph/tests/architecture/test_stationarity_modules_agree.py`.

**PARTIALLY SETTLED the same day by PI decision** ("일단 B로 진행"), and only partially — on purpose:

* **ADOPTED — the self-referential refusal.** `min_windows` and `min_tau_windows` were **never the same
  knob**, which this document's earlier draft also got wrong: the engine's counts multiples of an
  EXTERNALLY supplied driving time (the 2.5 s NMII bound-head lifetime → 12.5 s required), the observation
  stack's counts multiples of the MEASURED `tau_int` (0.08 s to 29 s over AR(1) ρ 0.5–0.999). Neither
  dominates. The observation stack names the external-bound-only design as a failure mode — *"A series
  can clear that bar and still be hopeless for its own correlation time"* — so the engine module gained
  `min_tau_windows = 50.0` (`n_eff ≥ 25`) **in addition to** its driving-time bound, which is exactly the
  composition that module recommends. Both are recorded in the returned contract. Verified with a positive
  control (N = 2,500, `n_eff` 21.1 → TOO_SHORT on the new reason) and a negative one (N = 6,000,
  `n_eff` 39.8 → STATIONARY). Disagreement **30.7% → 24.8%**; 84 existing tests still pass.
* **HELD — `drift_sigma`.** The adopted module is **looser** here (3.0 vs 1.0), and the decision to adopt
  it rested on this session presenting the two window knobs as comparable. Applying it flipped cases to
  STATIONARY — at ρ=0.9 with a linear drift a window scoring `|drift|/std = 1.00` is admitted by 3.0 and
  refused by 1.0 — and moved the dominant disagreement from 74 cases to **89**: loosening made agreement
  *worse*. Moving a gate threshold on a mistaken premise is what the charter forbids, so it was reverted
  and is surfaced here rather than applied.

  **RESOLVED 2026-08-20 when the PI asked whether 1.0 is too strict. It is not, and 3.0 is not merely
  looser — it is above what the statistic can reach.** `drift_over_std` is computed on the window the
  equilibration search KEEPS, and that search cuts the drifting head off, so the ratio saturates. Measured
  over AR(1) ρ=0.9 series carrying a real linear drift of 3× and 6× their own standard deviation:

      real drift       caught at 1.0   at 2.0   at 3.0        peak drift_over_std observed
      1× std                 10%          0%       0%
      2× std                 52%          0%       0%
      3× std                 94%          1%       0%              2.10
      6× std                100%         97%       0%              2.73

  and under the null — genuinely stationary, burned-in AR(1) — the **false-DRIFTING rate at 1.0 is
  0–2.4%** across n_eff from 41 to 1,361.

  So 1.0 sits inside the statistic's reachable range at a ~2% significance level, and **3.0 sits outside
  it: a threshold no series can attain is a gate that cannot fail** — the same defect §3 of this document
  diagnoses in `balance_ok`. Adopting it would have reproduced, inside the stationarity criterion, the
  exact failure this whole analysis exists to name.

  ⚠ Scope, so this is not over-read: a 3.0 run does not call everything stationary — the transient guard
  still fires on an obvious ramp. What dies is the DRIFT TEST, so the cases that flip are the intermediate
  ones: real drift with no sharp transient for the guard to catch. Pinned by
  `test_drift_sigma_3_would_be_a_gate_that_cannot_fire`.

⚠ Still unsettled, and now visible: 24.8% of series still disagree with **identical** window bounds, so
the residual is in the IMPLEMENTATIONS — equilibration search, transient guard, and which window the
drift is measured over — not in the constants. Matching numbers was necessary and is not sufficient;
consolidation needs one implementation.

---

## 6. The plan — four levels, in order. Each has one thing it must prove.

The ordering is not preference. Each level is meaningless before the one above it, and that is the whole
argument for not starting at L1 as the e1 document's D-1..D-3 implies.

### L0 — make a target exist. *Nothing about solvers happens here.*

1. **The two resting bound-myosin constants** (bound fraction, per-head isometric tension). PI-GAP, no
   default, and they are what generates the only stress that can oppose turgor.

   ⚠ **AMENDED 2026-08-20, same day, by the PI: this is not a request for a sourced point value.**
   *"우리는 어떤 것에 대한 기준을 잡지 않고 생리 범위만을 일단 테스트용으로 잡는다. 그리고 모든 것을
   믿지 않고 오직 테스트로 생각해서 돌린다."* That is **standing ruling #1 applied here** — *every
   physiological value is a declared axis with a band and the scope that band is true for* — and it is
   already what this module was built for: `resting_setpoint.py:29` reads *"The physiological values do
   NOT live here; a TEST fraction proves the MECHANISM, never the native gate."* The ask is therefore a
   **band + its scope**, not a citation, and the first version of this section asked wrongly.

   Three things this does NOT relax, each of which is how a band quietly becomes a value:

   * **The whole band is the result.** If only part of the band admits a balance point, the finding is
     *"the model balances over this sub-range"*. Reporting the sub-range as the answer is choosing a
     constant after seeing the data, which is the one thing the charter names. A band protects against
     that only when the band is what gets reported.
   * **The scope clause is the load-bearing half.** `k_xl` went wrong with a correctly-cited number:
     820 pN/nm IS in Ferrer 2008, and it is true of a 2 Å binding barrier, not of a 160 nm dimer under
     network load. A band without the scope it holds over reproduces that error with more digits.
   * **Standing ruling #1 caps its own yield** — at most the SVD rank of the log-coordinate Jacobian,
     computed BEFORE the sweep, may be reported as inferred.

   And one mechanical consequence: `RestingBoundMyosinSetpoint.__post_init__` **raises** on an empty
   `source` (`resting_setpoint.py:80`). A test axis does not evade that guard — it satisfies it with a
   string that says what it is, so the run record carries *"test axis, not a sourced value"* into every
   artifact instead of leaving the field to be misread later as provenance.

   ⚠ **And the FRACTION may not need a band at all.** `hand.py:187` and `:230` both state that the bound
   fraction *"EMERGES from k_on/k_off (P2), never an imposed duty ratio"*, and
   `engaged_fraction_steady(k_on, f_stall_head, k_off0, f0)` already computes it. Declaring a band on a
   quantity this engine generates is the lumped choice; the mechanistic one puts the band on the KINETICS
   and lets the duty ratio come out. The obstacle is that the resting inner solve **freezes kinetics**, so
   today the seeded fraction is an initial condition that persists — which is a solver-architecture
   question, not a constant. Recorded as the fork it is: **band the fraction, or band `k_on`/`k_off0`/`f0`
   and let the fraction emerge.**

   A bracket already exists on the other one and was never surfaced as a band:
   `ensemble_stall_analytic.py:67` — **per-head isometric stall, I0-B3 GAP: 0.5 vs 2.0 pN** (and
   `:66`, heads per half-filament: 10 AFINES vs 28-30 Billington).
2. **The `k_xl` axis.** `CROSSLINK_STIFFNESS_AXIS_2026-08-16.md` shows production runs
   `k_xl = 8.2e5 pN/µm` — a **binding-barrier curvature over ~2 Å**, used as the spring constant of a
   35–160 nm protein. It is **9.3× stiffer than the filament it connects**, while `kim_network.py`'s own
   docstring states the model's valid regime as `k_xl << EA/L_seg`. We are an order of magnitude into the
   opposite regime, and `dt_mu = 0.1/kmax` is set by it.

   ⚠ These are separate: (1) decides whether a balance point exists, (2) decides whether the cortex
   mechanics are in the right regime at all. **Neither is a solver fix and both are PI decisions.**

**PROVES:** a configuration exists whose residual can decrease.
**NEGATIVE CONTROL, declared now:** with the tension source off, the residual must **not** descend. If it
descends in both arms, the descent is not the mechanism and L0 has failed, not passed.

### L1 — define the state, not the solver

D-1 (observables + physiological bands), D-2 (stationarity statistic, window, **judged against seed scatter
— never `sem`, which is wrong by 67×**), D-3 (perturbation return + initial-condition independence).
Engineering prerequisite: **γ onto the resting path** — the estimator exists; wiring it is not a decision.

**PROVES:** acceptance is a statement about the cell that can be checked **without knowing which solver
ran**. That property is the entire point, and it is what `balance_ok` lacks.

⚠ γ sits on the blocklist ((c) 3 and (c) 17). Using it to *define* a state is not the same as *quoting* a
magnitude, and the distinction has to be written into D-1 explicitly or the first run will be misread.

### L2 — a gate that cannot fail is not a gate

* `balance_ok` **stays**, renamed to what it tests — adjoint/wiring closure. It is a good check. It is not
  acceptance.
* every criterion promoted to gating status ships with a **demonstrated negative control first**. This
  project already does this well (GATE-B's positive control injects one head's stall force one-sidedly and
  the same path rejects it). The rule is only that it becomes mandatory rather than admirable.

### L3 — collapse the thicket

**One** definition of an accepted physical step. Everything in the §4 census that does not gate something
is demoted to *diagnostic* in name as well as in fact. 17 verdict strings is not rigour; it is the absence
of one decision, distributed.

---

## 7. What changes about how the engine is run

**Stop asking "which solver".** Five have been tried; the answer will not be the sixth. The 2026-08-15
eliminations are worth keeping precisely because they close that direction: budget, construction,
conditioning, coupling and shell feasibility are all excluded.

**Stop treating the arena as blocked on this.** Its ID-range addressing is force-identical by measurement
(`78942fc4`: relative residual 1.108e-16, leak exactly 0.0 pN). It changes no physics number by
construction, so it neither helps nor hurts L0–L3 and should not wait on them.

**Expect L0 to be cheap and L1 to be the work.** L0 is two constants and one axis decision. L1 is wiring
γ, choosing bands, and running the two protocols D-3 needs — and it is where a real GPU budget goes.

---

## 8. What this document does not do

It selects no observable, sets no band, writes no threshold, picks no window, and assigns no value to
either PI-GAP constant or to `k_xl`. Doing any of those **after** seeing 2026-08-15's numbers is what the
charter forbids. It also does not retract any of the five eliminations; it retracts only the inference
drawn from them, which the e1 document had already marked PROVISIONAL.

**Open, and in this order:** the `stationarity.py` duplication (a defect, not a decision); the two
resting-myosin constants; the `k_xl` axis; then D-1, D-2, D-3.
