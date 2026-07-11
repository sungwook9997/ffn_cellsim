# DCM timescale-gap attack — 2026-07-11 (PI: "attack the timescale gap directly")

## Goal
The DCM (and FF) reproduce the *mechanism* + *kinetics* of active unjamming/compaction but cannot reach the
*biology-time equilibrium*: physiological unjamming/rearrangement needs minutes–hours, and at the feasible integration
step accel_dt=8e-4 that is ~10⁶ steps ≈ 15 h wall/run (infeasible). PI directed attacking this gap directly. The
cleanest first lever — since the runtime already has an **implicit integrator + IPC** (unconditionally stable for stiff
forces) — is to check whether cranking `accel_dt` up 10–100× stays *accurate*, which would cover 10–100× more physical
time at the same wall-clock. (Kim's τ_p=2^N·τ is compaction-rate-specific; multi-scale subcycling is more complex — so
large-dt was tried first.)

## Step 1 — naive large-dt push in v0-mode — **FAILS (non-convergent)**
v0-mode chosen (per-node force = γ·v0 → slow drift, v0·dt tiny at any dt, so CFL is bounded — the honest large-dt-safe
mode). Held physical time = 24 s constant (v0=20 µm/min, τ_p=600 s, N=100), varied accel_dt; settle held at 0.8 s
physical. If implicit+IPC were accurate at large dt, A/A0 + shape index would converge across dt.

| accel_dt | steps | wall | A/A0 | V/V0 | shape index s (→) | per-cell disp mean/max | verdict |
|---|---|---|---|---|---|---|---|
| 8e-4 | 30000 | 40 min | 1.154 | 1.000 | 4.855→**5.003** | **3.07** / 21.2 µm | **accurate** (reference; mild unjam) |
| 2e-3 | 10000 | 4.5 min | 1.002 | 1.000 | 4.855→**4.892** | **0.153** / 1.26 µm | **FROZEN already at 2.5× dt** (1/20 the drift) |
| 8e-3 | 3000 | 77 s | **31.7** | 1.058 | 4.855→**13.72** (unphysical) | 0.37 / 29.4 µm | **UNSTABLE** — cell mesh inversion → blow-up |
| 2e-2 | 1200 | 68 s | 1.002 | 1.000 | 4.855→**4.916** | 0.287 / 3.5 µm | **FROZEN** — drift over-damped |

The accurate window is essentially **dt=8e-4 only** — even 2.5× (dt=2e-3) loses 95 % of the creep. So the force-based
path gives ~**1× acceleration** = the gap is NOT closable by dt on this integrator.

The three results are wildly different → **not converged**. Visual-verified: dt=8e-3 collapses to a tiny mangled blob
with an inverted-cell red spike (`s14_dt8e-3_last.png`) = numerical blow-up, not physics; dt=2e-2 aggregate is unchanged
(frozen).

## Root cause — backward-Euler over-damping of SLOW active forcing (+ intermediate-dt contact instability)
The "implicit handles large dt" assumption holds for **stiff, fast** forces (elastic, contact) — backward-Euler is
L-stable, so it damps the fast elastic relaxation and stays stable. But that same L-stability **over-damps SLOW
forcing**: the active motility drift (a near-constant force) is under-resolved at large dt — backward-Euler relaxes the
cell back elastically faster than the slow drift can move it, so the net motion is suppressed (dt=2e-2 froze at 1/10 the
reference displacement). And at intermediate dt (8e-3) the contact/elastic solve loses robustness and a cell inverts →
blow-up. So there is **no usable large-dt window** on the current backward-Euler + explicit-safeguard integrator: too
small = infeasible, intermediate = unstable, large = frozen. **Cranking dt does NOT close the gap.**

## The real fix (integrator-level, for PI decision)
The gap is not a parameter; it is the time integrator's response to slow active forcing. Options:
- **(a) A-stable, non-over-damped integrator** — implicit midpoint / trapezoidal / BDF2 instead of backward-Euler.
  A-stable-but-not-L-stable schemes resolve slow forcing accurately at large dt (no artificial over-damping) while
  keeping stability. This is the principled fix; it is a real integrator change in the Warp runtime.
- **(b) Operator-split exact drift** — advance the motility as a *prescribed* slow translation of cell centers
  (exact at any dt: Δx = v0·dt·p̂), then relax elastic/contact implicitly. Decouples the drift from the over-damping;
  the SPV picture (centers follow active velocity, tissue relaxes). Smaller change than (a); directly targets the freeze.
- **(c) Robustify contact + moderate dt** — fix the intermediate-dt cell-inversion (continuous collision / tighter IPC
  barrier) to open a 5–10× window; less ambitious, and 5–10× likely still short of the minutes–hours gap.

Recommendation: **(b) operator-split exact drift** is the smallest change that directly addresses the observed freeze
and is SPV-faithful; **(a)** is the more general/correct fix if (b) is insufficient. **(c)** alone won't span the gap.

## Step 2 — OPERATOR-SPLIT drift (option b) — IMPLEMENTED, under test
Prototyped (b): `active_drift_translate_kernel` in `dcm_active_motility_warp.py` + `--motility-split` in
`run_spread_overnight` (wired through `run_decohesion(motility_split=…)`). At the top of `step_once`, before grids /
forces / the xₙ freeze, each node is translated by `v0·dt·p̂_c`; the force-based self-propulsion is skipped (v0-mode).
So the backward-Euler anchor is x* = xₙ + drift and the implicit IPC-Newton relax runs from there → Lie-Trotter
drift-then-relax, drift exact at any dt. Local CPU smoke PASS (compiles, runs, V/V0 conserved).

**Convergence test running** (s15, same config as step 1 + `--motility-split`, dt ∈ {8e-4, 8e-3, 2e-2}):
- **Falsifiable predictions:** split @ dt=8e-4 must match the force-based reference (per-cell disp ≈ 3.07 µm, s→5.0)
  — small-dt consistency; and if the split works, split @ dt=2e-2 must ALSO give ≈3 µm (not freeze) = **gap closed,
  25× acceleration**. If split @ 2e-2 also freezes, the jammed-creep is intrinsically small-dt (T1-event-limited) and
  the gap is deeper than the integrator — escalate to (a) an A-stable higher-order scheme or accept the honest limit.

## Fix options (for PI if step 2 is insufficient)
- **(a) A-stable, non-over-damped integrator** — implicit midpoint / trapezoidal / BDF2 instead of backward-Euler;
  resolves slow forcing at large dt without artificial over-damping. The general/correct fix; a real runtime change.
- **(b) Operator-split exact drift** — IMPLEMENTED above (under test).
- **(c) Robustify contact + moderate dt** — fix the intermediate-dt cell-inversion to open a 5–10× window; likely
  still short of the minutes–hours gap on its own.

## Status
- Naive large-dt push (force-based) characterized + REFUTED (accurate window ≈ dt=8e-4 only), visual-verified.
- Operator-split (b) implemented behind `--motility-split`; smoke PASS; convergence test in flight (s15).
- Next: read the s15 convergence — if it closes, a physiological-v0 + large-dt hero run to τ_p=600 s; else escalate to (a)/PI.
