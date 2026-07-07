# DCM hard-constraint reformulation — large-dt implementation plan (2026-07-07)

**Motivation.** DCM's dt is CFL-capped by STIFF penalties (turgor K_vol=7.73e5, contact rep=2e8,
edges k_edge). baoab dt=8e-6; implicit IMEX dt=8e-4 (100×). Reaching biological minutes needs ~1e6
steps → infeasible. This is fidelity red-flag #2 (timescale). FF/Cytosim avoids it by keeping stiff
crosslinks out of the explicit step (discrete KMC events) — but DCM's contact/turgor are CONTINUOUS
constraints, so the analog is not "exclude" but "solve as hard constraints." PI 2026-07-07: build this
(#1), on the user's go.

Design workflow: XPBD vs proper-IPC vs hybrid, adversarially reviewed, synthesized.

## RECOMMENDATION: finish the IPC (projected-Newton incremental potential), NOT XPBD

The codebase is already ~80% toward IPC (`dcm_contact_implicit_warp.py`: `nearest_face_ipc_kernel`,
`_ipc_b/_bp/_bpp` log-barrier, `ccd_alpha`, `make_contact_hess_apply`). Complete it and route turgor
through the same Newton loop. Reasons, each tied to a hard rule:

1. **Cells NEVER interpenetrate** → log-barrier guarantees gap d>0 *by construction* (barrier diverges
   as d→0) + CCD-filtered Armijo line-search keeps every Newton iterate strictly feasible at ANY dt.
   XPBD-compliant contact still permits load-proportional penetration (the penalty failure mode);
   XPBD-hard is a one-shot projection needing CCD substepping anyway. IPC = the STRONGER guarantee,
   lower marginal effort.
2. **No magic numbers** → IPC's κ=rep, d̂=c_rep are already derived; every other IPC constant (tol,
   max_newton, Armijo c1, backtrack β, CCD η) is a grid-invariant numerical control, not physics.
   XPBD compliance risks a dt-dependent *hidden* effective stiffness (invisible at a single dt) unless
   the overdamped α̃=α (no /dt²) mapping is exactly right. **Reject Li et al.'s adaptive-κ** (= tuning
   a constant to converge = rule violation); keep κ=rep fixed, add iterations if Newton struggles.
3. **Physiological dP0 + finite K_vol** → turgor must stay COMPRESSIBLE; a hard V=V0 equality is wrong
   physics (incompressible, drops dP0). Turgor needs NO reformulation: it's a smooth C¹ implicit force
   with an analytic K_vol diagonal (`build_diagA`), converges in 3-5 CG iters. **The same Newton loop
   that finishes contact lifts the turgor large-dt accuracy cap for free, with exact parity.**

**Why the current `ipc=True` silently stalls** (confirmed in code, `run_decohesion:1189-1220`): it
does a SINGLE linearized `device_cg` solve then scales by `ccd_alpha`. Under a strong cadherin bundle
CCD clamps α small → the step stalls at a CCD-throttled sliver (runs + looks stable but never reaches
the full-dt feasible minimizer). Missing piece = Newton re-linearization loop + energy line-search
around that existing single step. The `a=γ/dt` effective-mass regularizer is already present.

## Build order (parity gate at each step; NO dt push until the prior gate passes)

All gates small-N (N=2..12), on-device, at dt=8e-6 unless noted.

- **Step 0 — Telemetry.** Log realized-α, Newton residual, iter count on the current `ipc=True` branch.
  *Gate 0:* reproduce the silent stall (α clamps small, V/V0·pen drift vs penalty) — the baseline to beat.
- **Step 1 — Energy kernels** (`dcm_contact_implicit_warp.py`): `ipc_barrier_energy_kernel` (mirror the
  force kernel's IDENTICAL nearest-face pick), + edge/turgor/inertial energy reductions + host
  `total_energy(x)`. *Gate 1:* `FD(total_energy) == −F_total` to 1e-6 (energy–gradient consistency; the
  subtle line-search trap). Blocks Step 2.
- **Step 2 — Device Newton driver** (`dcm_warp_implicit.py`): `ipc_newton_step` wrapping the existing
  `device_cg` as inner solve; residual G=a(x−xn)−F_total(x); re-linearize barrier+turgor each iterate;
  line-search α0=ccd_alpha then Armijo backtrack. *Gate 2:* energy monotone-decrease, residual clears
  tol in bounded iters, pen==0 every iterate.
- **Step 3 — Driver wiring** (`run_decohesion:1189-1220`): replace single-step+one-shot-CCD with the
  Newton loop. Broad-phase once/step on xn; soft drivers (cadherin, wetting, ECM, gravity) explicit,
  lagged at xn; turgor re-evaluated each iter. Keep de-cohesion cadherin-rupture-EMERGENT. **Gate 3
  (THE hard penalty stiff-limit parity gate, PI-visible):** at dt=8e-6, finished-IPC vs penalty on
  V/V0, max node force, pen_frac, cadherin churn (N=2 and N=12). Fail → HALT to PI. No dt push before.
- **Step 4 — dt ramp** 8e-4→8e-3→8e-2, assert pen==0 every step + V/V0 in the penalty-parity band.
  *Gate 4:* physics quantities dt-INVARIANT across 1×…1000×; report Newton-iters + wall/step per dt.
- **Step 5 (optional, flagged)** — XPBD-compliant turgor warm-start (α_vol=V0/K_vol, overdamped α̃=α,
  dP0 explicit) ONLY if Step-4 profiling shows the dense per-cell volume coupling dominates Newton cost.
  *Gate 5:* dimensional sanity + cross-dt parity that the warm-start does not shift converged V/V0.

Test: extend `_ipc_full_test` into a Newton-loop gate covering Gates 1-4 (measure realized-α + residual
explicitly, not just final pen, or the silent-throttle failure stays hidden).

## Expected gain + honest cost

- **dt:** 8e-4 → **8e-3…8e-2 (10-100×)**, minutes-regime in O(1e2-1e3) steps (vs O(1e6) baoab).
- **Per-step cost:** ~2-6 Newton iters/step × (one device_cg 3-5 iters + 1 barrier + 1-3 line-search
  evals) ≈ 4-8× the current single-step cost → **net ~2-20× wall-clock, WITH guaranteed non-penetration**
  replacing the pen 1.5-2.6 tunnelling. Ceiling: at extreme dt the barrier ≫ a/dt, Newton needs more
  iters → cap `max_newton` + adaptive-dt fallback (shrink dt for that step, report iters).

## No-magic-number check
Every new constant is a numerical control, not physics: κ=rep, d̂=c_rep (derived); tol/max_newton/c1/β/η
(grid-invariant convergence controls). Optional Step-5 α_vol=V0/K_vol, α_edge=1/k_edge, w=dt/γ — all
derived from existing sourced constants. dP0 stays an explicit physiological preload, folded into nothing.

## Fallback if Gate 3 (parity) fails
1. Suspect energy–gradient inconsistency (nearest-face mismatch energy↔force) → re-run Gate 1.
2. If parity holds at 8e-6 but breaks on the dt ramp → it's ACCURACY not method: add Newton iters /
   tighten tol; do NOT touch κ or K_vol.

---

## PROGRESS LOG

- **Step 0 ✅ (d6fa1d6)** — IPC_TELEM telemetry. Found: single-step stable to dt=8e-2 but cadherin
  kinetics (36ms) vs step (80ms) is the real multiscale bottleneck; σ_vm 2× error at 8e-2 → Newton needed.
- **Step 1 ✅ (8053b33 energy kernels + Gate-1 flag; 5e30c0d COMPLETE)** — barrier energy↔force match.
  **PI decision A** (add the exact ∇area term, not freeze-area — full variational IPC gradient per the
  fine-grained rule). `nearest_face_ipc_kernel` gained `F=−dE/darea·∇area` on the 3 neighbour face
  vertices (∂A/∂a=½ n̂_f×(c−b) cyclic, Σ=0 momentum-conserving), both barrier + inside-feasibilization
  branches. **Gate 1 FD ∇E==−F: barrier 12%→1.97e-9 PASS** (edge/turgor machine-precision). Node normal
  force unchanged (2.2e-16), momentum conserved (2.18e-16), IPC-vs-penalty parity harness un-regressed.
- **Step 2 ✅ (this commit)** — device projected-Newton driver `ipc_newton_step` in `dcm_warp_implicit.py`
  (wraps `device_cg`; residual G=a(x−xn)−F_total; CCD-filtered Armijo energy line-search on
  Φ=½a|x−xn|²+U). Energy kernels for the merit: barrier+edge+turgor (Step 1) + new inertial +
  soft-linear (lagged drivers) in `dcm_ipc_energy.py`. **Gate 2 (`tests/test_ipc_newton_gate2.py`, CPU):
  PASS** — feasible-start + 4·Fcap pull, 30-step trajectory: worst_pen=0.00 at dt 1×/100×/1000×,
  energy monotone, residual→tol; Newton iters 3.0→4.6→7.4 as dt scales (loop earns its keep at big dt).
  Note: deep penetrating-START recovery is slow (only OUTSIDE barrier diverges) → **production must start
  feasible/watertight** (the Voronoi-confluent/dispersed init — matches the init strategy).
- **Step 3 ✅ (baee34f)** — wired the Newton loop into `run_decohesion` additively (`--ipc-newton`,
  needs `--ipc`); working `ipc=True` single-step path untouched (else-branch, no regression). RHS split
  once/step: F_soft = F_all(xₙ) − stiff{turgor,edges,barrier}(xₙ), lagged; the loop re-linearises the
  stiff set. Added `turgor_energy_pc_kernel` (per-cell dP0, matches pc + osmotic) + a "stalled → converged"
  early-exit (the relative |G|/|G0| test is noisy near equilibrium). **Gate 3 (dt=8e-6, N=2/12): PASS** —
  Newton-vs-single V/V0 Δ=3.8e-5 (N2) / 5.6e-5 (N12), pen=0, all gates PASS; no regression on single/penalty.
- **Step 4 ✅ (dt ramp, `tests/test_ipc_newton_gate4.py`) — split into two contracts:**
  - **4a STABILITY + NON-PENETRATION invariance = PASS (the IPC payoff).** Newton stays finite +
    penetration-free at 1×…10000× dt; at 8e-2 (10000×) Newton pen_peak=0.000 while the single linearised
    step **DIVERGES** (pen_peak=3.45, V/V0 collapses 0.917). This is exactly what finishing IPC buys.
  - **4b FULL-PHYSICS V/V0 invariance = FAILS by ~2% at ≥100× → SOFT-FORCE LAG ceiling (SURFACE TO PI).**
    V/V0 drifts 1.0204(base)→1.0000(≥100×). Verified NOT fixable by more Newton iters / tighter tol
    (stiff solve already converged) and NOT a wiring bug (the single-step drifts identically). It's the
    LAGGED soft/cadherin forces being under-counted over a huge step — the compaction driver. **Did NOT
    loosen the V/V0 contract.** The honest conclusion: mechanics takes large dt STABLY, but accurate
    large-dt physics for the 24–48h aggregate needs the soft/cadherin drivers SUB-CYCLED at their own
    (~36 ms) timescale — the multiscale layer the memory + this plan's honest-cost note anticipated.
- **Step 5 🔜 (PI decision)** — cadherin/soft sub-cycling (mechanics large-dt + soft micro-steps), OR
  accept mechanics-only large-dt with the documented soft-lag caveat, THEN the Voronoi-confluent aggregate
  run on gbook GPU. PI picks the path.

---
*Full workflow transcript (XPBD + IPC designs + synthesis): run wf_056df109-3dd. The adversarial-reviewer
agent hit the StructuredOutput retry cap; its concerns were internalized by the synthesis (energy–gradient
trap, dt-dependent-compliance magic-number risk, incompressibility error of a hard V=V0 constraint).*
