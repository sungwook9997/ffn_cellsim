# (e) 1 — the gate amendment, put into a form the PI can actually decide

**Status:** PROPOSED → this document does not decide anything. Written 2026-08-16 by the `arena/solver`
session at PI request, after the PI's own objection made it clear the amendment is not merely *pending*
but currently **undecidable as stated**.

**What it is about:** `STATE.md` (e) 1 — replace the stage-1 acceptance criterion *convergence* with
*stationarity*. PI HELD it 2026-07-28. It has been held for 19 days, and every native run in that window,
including all nine of 2026-08-15, ran under the un-amended criterion.

---

## 1. The PI's objection, stated precisely

> *The initial state is not an equilibrium to begin with, so assuming it converges makes no sense. Without
> knowing whether that state is a stable convergence or a state-limited one — without specifying the state —
> how would you know it converges at all?*

This is two objections and both stand.

**(a) The target may not exist.** The current criterion asks the residual to go to zero. A living cell is
not in static mechanical equilibrium: myosin binds and unbinds, cross-links turn over, actin treadmills.
Force balance holds *statistically*, not exactly. `residual → 0` is a criterion imported from passive solid
mechanics, and applied to an active material it can be unreachable in principle — no solver fixes that.
This was already found on 2026-07-28: *"the runtime was solving a static equilibrium that does not exist —
the solver's symmetric tangent quietly asserted a potential existed, and it had never been checked."* That
finding is what produced (e) 1 in the first place.

**(b) Even reaching it proves nothing.** Suppose the residual did go to zero. What is then established is
*"the solver stopped moving"*. It does not distinguish

* the physiological resting state,
* a local basin the initial condition happened to fall into,
* a configuration that a different initial condition would not reproduce.

Nothing in the current criterion separates these, **because the target state is never specified
independently of the solver's own stopping behaviour.** A criterion that is defined by the process cannot
certify the process.

⚠ Objection (b) is the one that survives every fix to (a). Giving the cell cortical tension may make a
balance point exist; it does not make "we stopped here" mean "this is the resting cell".

---

## 2. Why the amendment cannot be ratified as currently written

The amendment says *stationarity instead of convergence*. Stationarity **of what** is not specified, and
that is not a wording gap — the quantities it would have to be stationary in are largely not measured.

### 2.1 What the resting driver emits today

Read off `StabilityReport` (`driver.py:986-1031`) and the build ledger of a committed full-native run
(`arena_prestudy/s1_pressure_on.json`):

| Candidate state observable | Emitted by the resting run? | Field |
|---|---|---|
| Cell radius (mean / min / max) | **YES** | `r_mean_start`, `r_mean_end`, `r_min_end`, `r_max_end` |
| Fluid content (volume proxy) | **YES** | `content_start`, `content_end`, `content_delta` |
| Centre-of-mass drift | **YES** | `com_drift_um` |
| Bound myosin head count | **YES** | `myosin_bound_after_step` — currently `0` |
| Bound ERM tether count | **YES** | `erm_bound_after_step` — currently `642` |
| Control volume | **YES** | `pressure_control_volume_um3` |
| **Cortical tension γ** | **NO** | not in the report, not in the ledger |
| **Cortex thickness** | **NO** | not emitted |
| Per-population residual | **NO** | only the scalar `max\|PF\|` |

### 2.2 The observable the state most needs is the one that is missing

γ is what a resting cell *is*, mechanically, and it is what every external measurement this project intends
to compare against reports. The resting driver does not compute it. γ exists only in the GATE-B slice
drivers (`ac_gate_b_sf_motor_native.py`, the cortex-motor gates).

**And γ is simultaneously on the blocklist.** `STATE.md` (c) 3 retires *"γ = 3.72 pN/µm and every SF/cortex
motor tension"*; (c) 17 retires *"every cortex-motor γ at dt=0.01 — magnitudes AND ratios"*. So the
quantity that would define the state is (i) not measured on the resting path and (ii) not quotable where it
is measured.

⚠ This is the actual blocker on (e) 1, and it is not a wording problem. **A stationarity gate cannot be
written against an observable the run does not produce.**

### 2.3 And the current gate cannot fail for non-convergence anyway

Already recorded on (e) 1: stage 1's inner gate accepts on `balance_ok_d`, but `sf_motor_slice.py:65` says
that predicate tests **adjoint closure**, not convergence. So a convergence gate is resting on a predicate
that cannot fail for non-convergence. Amending *what* the gate tests without fixing *what it reads* moves
the problem rather than closing it.

---

## 3. What "specifying the state" would concretely require

Three things, and each is a separate decision.

### D-1 — Which observables define the resting state, and in what band

The state must be declared by quantities that are measurable **both in the engine and in a real cell**,
so the comparison the project exists for remains possible. Candidates, with current feasibility:

| Observable | In engine now | In literature | Note |
|---|---|---|---|
| cell radius R | yes | yes | already emitted; weakest discriminator (geometry is imposed at build) |
| cell volume V | yes (fluid content) | yes | responds to turgor; currently pinned by Π₀ |
| cortical tension γ | **needs work** | yes (AFM, aspiration, dynamic AFM) | the strongest discriminator; **not emitted on the resting path** |
| bound myosin fraction | yes (count) | yes | currently 0 by construction; becomes live only with the resting setpoint |
| cortex thickness | **needs work** | yes | not emitted |

**Decision:** which of these constitute the state, and what band each must lie in for the configuration to
be called *resting*. The band is a physiological claim, not a numerical tolerance, and per the charter it
must be written before any run.

### D-2 — What stationarity means, quantitatively

A criterion of the form: *observable O is stationary over window W if its mean is stable and its
fluctuations are bounded*. Two parameters have to be fixed before a run:

* **the window W.** Already partly held under (e) 1: `min_windows` — 5 driving lifetimes = 7.28 s against
  statistics needing ≤ 2.5 s, a 2.4× margin **measured, not assumed**.
* **the variance the mean is judged against.** ⚠ `STATE.md` (f) carries a standing measurement rule that
  binds this directly: compare γ against the **seed scatter**, not the within-run `sem`. Three replicates
  give S = 0.0668 pN/µm = 1.95% of γ, **67× the sem**, so one seed cannot resolve a γ difference below
  ~5%. **A stationarity test built on within-run variance would be wrong by that factor**, and this is the
  single most likely way to write the amendment and still be measuring nothing.

### D-3 — How stable-versus-conditional is distinguished

This is the PI's objection (b), and it cannot be answered by any single run. It needs a protocol:

* **Perturbation return.** Displace the converged state and check it comes back. A basin test.
* **Initial-condition independence.** Start from different constructions and check the same observables
  land in the same band. Without this, "stationary" is compatible with "stuck where it started".

⚠ Neither exists today. Both are cheap once D-1 and D-2 are fixed, and both are meaningless before.

---

## 4. What the 2026-08-15 results contribute to this decision

Five candidate causes for the standing non-descent were eliminated. **These stand under either criterion**
— they would have been confounds for a stationarity gate exactly as they were for a convergence gate:

* non-conservativity — excluded (no head bound at rest; the inner solve freezes kinetics)
* cross-link pre-stress — excluded (extension exactly 0.0 across all 1,413,720)
* iteration budget — excluded (20× and 66× both end **higher**)
* solid–fluid coupling — excluded (Π₀ = 0 with Biot fully composed, 67 sub-cycles, reproduces the
  no-fluid arm to the digit)
* infeasible shell — excluded (capacity ratio 0.33 → 16.6, still diverges, 900× faster)

Two results bear on the amendment directly:

* **Explicit is non-monotone at two population scales five orders apart** (156 DOF and 1.65 M DOF). More
  iterations make it worse. That is not the signature of a solver that needs a bigger budget; it is
  consistent with a target that is not there.
* **The driver documents its own resting configuration as having no solution.**
  `_preload_erm_resting_balance` (driver.py:233-241): *"a relaxed cortex plus a positive turgor is not a
  valid resting baseline — there is no force-balanced equilibrium there, which is why no inner solver
  converges from it."* Of 23 recorded runs carrying a resting residual, **none** ever applied either
  preload.

⚠ **RETRACTED HERE:** the 2026-08-15 conclusion *"the operator is the only cause left for the resting
non-descent"* presupposes that a target exists and is therefore not established. It is withdrawn to
PROVISIONAL pending D-1..D-3. `STATE.md` already marks it provisional; this document is the reason.

---

## 5. What the PI is actually being asked

1. **D-1** — which observables define the resting state, and in what bands.
2. **D-2** — the stationarity statistic and its window, judged against seed scatter rather than `sem`.
3. **D-3** — whether the stability protocol (perturbation return + initial-condition independence) is part
   of the gate or a separate follow-on.

And one prerequisite that is engineering rather than decision, listed so it is not mistaken for one:

4. **γ must be emitted on the resting path.** Whatever D-1 selects, if γ is in it, the resting driver has
   to compute and record it. It currently does not.

## 6. What this document does NOT decide

It selects no observable, sets no band, writes no threshold, and picks no window. Doing any of those after
seeing 2026-08-15's numbers is precisely what the charter forbids. It also does not touch the three
resting bound-myosin PI-GAP constants — those are a separate open item, and they gate whether a balance
point exists at all (level 1), not what it would mean to have found one (level 2).
