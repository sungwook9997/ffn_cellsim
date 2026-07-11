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

**Convergence test (s15) — RESULT: the split FIXES the freeze + instability (partial success).**

| dt | method | A/A0 | V/V0 | shape index (→) | per-cell disp | verdict |
|---|---|---|---|---|---|---|
| 8e-4 | force (ref) | 1.154 | 1.000 | →5.003 | 3.07 µm | accurate reference |
| 8e-4 | **split** | 1.174 | 1.000 | →4.998 | **3.13 µm** | **matches ref** ✓ (small-dt consistency) |
| 8e-3 | force | 31.7 | 1.058 | →13.72 | 0.37 µm | blow-up |
| 8e-3 | **split** | 1.277 | 1.000 | →5.185 | **7.24 µm** | **stable** ✓ (blow-up fixed) |
| 2e-2 | force | 1.002 | 1.000 | →4.916 | 0.287 µm | FROZEN |
| 2e-2 | **split** | 1.277 | 1.000 | →**5.307** (frac 0.31) | **6.69 µm** | **NOT frozen** ✓ genuine unjam |

**What the split achieves:** (1) at dt=8e-4 it matches the force reference (disp 3.13 vs 3.07) → the split is correct at
small dt; (2) it removes the dt=8e-3 blow-up (cells stay intact, V/V0=1.0); (3) at dt=2e-2 (25×) it PRESERVES the drift
(6.69 µm vs the force path's frozen 0.287 µm) → the aggregate genuinely unjams (s→5.307, 31 % of cells past s0*).
**Visual-verified** (`s15_split_2e-2_last.png`): 100 intact cells, loosened + rearranged, cadherin bonds intact,
moderate stress — real collective unjamming at 25× dt, NOT the force path's mangled blob.

**The caveat (honest):** the large-dt split OVER-drifts ~2× vs the small-dt reference (6.7 µm vs 3.1 µm) and is not
cleanly monotone (8e-3: 7.24, 2e-2: 6.69; the 8e-3 run also carried pen 0.058). This is a Lie-Trotter O(dt) splitting
error (fewer relax steps between drift kicks → less elastic push-back → over-drift) plus some large-dt contact-resolution
slack. So the split **qualitatively closes the gap** (stable, physical, unjamming large-dt runs → 25× ≈ τ_p=600 s in
~40 min instead of ~15 h) with a **~2× quantitative accuracy caveat**.

**Next (fork):** tighten the accuracy — **Strang splitting** (symmetric half-kick / relax / half-kick → O(dt²)) and/or
more IPC-Newton iters at large dt — before quantitative claims; OR proceed to a physiological-v0 large-dt hero run now
as the qualitative gap-closing demonstration (unjamming at physiological v0 over ~physiological time, ±2× rate). Both
are legitimate; the split is the enabling result either way.

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
