# motor-nmii track — REPORT (I3 head-resolved NMII, analytic foundation)

Track: **motor-nmii** (I3 head-resolved Stam-Hocky NMII) · branch `ac/motor-nmii` · base `ac/new-engine`.
Status: I3 analytic oracle suite committed + visualized (41 host-numpy tests, all green) + Warp kernel SOURCE
authored (`hand.py` device KMC API + `minifilament_warp.py` `MyosinForce.accumulate`; CPU-codegen type-checked,
NOT launched — native gates are gbook-owned, serial). Follows the fluid-spine visualization pattern
(AC_PARALLEL_SESSIONS §1.9).

## What landed

- **`ac/motor/params_i0b3.yaml`** — the I0-B3 ledger. 11 GAP (`value: null`, surface to PI): `F_stall_head`,
  `N_side`, `v0`, `kappa_hill`, `k_xb`, `r0_head`, `L_bb`, `n_bb`, `k_on`, `d_step`, `dG_atp`; 3 provisional
  (`k_off0`, `x_beta`, `duty`). ⭐ **FV-form reconciliation** recorded as a THREE-way conflict on `kappa_hill`
  (Hill-1938 muscle 0.25 / Kovács-2003 archived NMII 0.5 / PI-2026-07-07 linear κ→∞) — the Hill machinery is
  built, the curvature is left to PI. ⭐ **`k_xb` MASTER-knob** flagged: physical 100–1000 pN/µm; the archived
  1 pN/µm is the broken-soft surrogate that breaks stall.
- **4 analytic oracle modules** (host numpy, ZERO Warp): `hill_fv_analytic`, `bell_kinetics_analytic`,
  `ensemble_stall_analytic`, `minifilament_topology`. All 7 I3 analytic gates green.
- **Warp SOURCE** (authored, not run): `hand.py` (per-head hand/KMC device API — I3 OWNS the §1.4 frozen
  interface) + `minifilament_warp.py` (`MyosinForce.accumulate` + topology build + head-load loop).
- **`INTEGRATION.md`** — the `ff/network_warp.py` lumped-`f_myo` / `myosin_kernel` deletion patch-notes
  (same-commit-as-fine-motor guard, double-count #6) + the full native-gate spec.

## Figures

Regenerate all: `PYTHONPATH=<repo> python ffn_sim/scripts/ac_motor_vis.py` (pure numpy/matplotlib on the
committed oracles). Every numeric constant is an I0-B3 GAP/provisional value **swept as an oracle variable**,
labeled as such on every figure — NOT a chosen value.

- **`figs/i3_hill_fv.png`** — single-head Hill force-velocity vs Hill-1938. Hyperbolas for κ=0.25 (muscle),
  0.5 (Kovács/archived NMII), 1.0, with the κ→∞ linear law (PI 2026-07-07) overlaid dashed; the Hill-1938
  reference point v/v₀=1/6 @ F=F_s/2 sits exactly on the κ=0.25 curve. Right panel: mechanical power P=F·v ≥ 0
  on [0,F_s], zero at v₀ and at stall, single interior max (the work-sign gate).
- **`figs/i3_bell_kinetics.png`** — Bell slip off-rate p_off=k_off0·exp(|f|/f₀) monotone in load (log axis,
  noted), p_off(0)=k_off0; and the EMERGENT engaged fraction φ_b=k_on/(k_on+p_off) that is self-limiting under
  load. ⚠ honest finding shown: with f₀≈7.1 pN ≫ F_head, φ_b≈0.99 even at per-head stall — so the ensemble
  stall is only slightly below the naive product HERE, but it is COMPUTED from kinetics, not imposed.
- **`figs/i3_ensemble_stall.png`** — ensemble stall EMERGES. Left: the stochastic per-realisation
  distribution (binomial heads) with the mean-field E[F]=N_side·φ_b·F_head < the naive N_side·F_head (imposed)
  line, the KB-3.18 50–100 pN band as a REFERENCE overlay (not a target). Right: emergent E[F] falls as the
  off-rate rises (self-limiting) — a hard-coded N_side·F_head would be flat.
- **`figs/i3_newton_kxb.png`** — per-head Newton closure residual → below 1e-12 in one step (linear system,
  log axis); and the k_xb MASTER-knob arbiter (loglog): working-stroke strain = F_head/k_xb lands in the 5–20
  nm physical band only for k_xb 100–1000 pN/µm, while the archived 1 pN/µm gives 2000 nm ≫ the 301 nm
  minifilament → "1 pN/µm breaks stall".

All figures overlay the closed-form oracle/band on the numeric result, annotate engine units (pN, µm/s, 1/s),
do not truncate axes, and (ensemble) show per-realisation spread + the mean overlay — the visualization-
integrity rules.

## I4-weave — the actin↔motor JOINT (production motor bound to the woven network)

Branch `ac/i4-weave-motor` (worktree). **What this increment delivers:** the production point-to-segment NMII
motor (`ac/motor/segment_motor.py` + `segment_query.py`, device-parity 5/5) is now bound to the ONE woven
actin network through a **first-class, reusable joint on `WovenCell`** — so the head-resolved minifilaments
grab points on the REAL filaments (cortex + SF + arc + …) and each bound head's `walk_dir` is overwritten from
the segment's actin barbed-end polarity. No lumped `f_myo`, no imposed `N_side·F_head` single link.

**Where the polarity overwrite happens (already device-resident, confirmed):** `attach_segment_gated_kernel`
sets `walk_dir[h] = q_barbed[h]` (the queried segment's barbed direction) on (re)binding, and
`_refresh_bound_walk_dir_kernel` keeps a bound head following its live segment; on detach `walk_dir` resets to
zero (passive). The bipolar-geometry default from `build_minifilament_nodes` is the OFF/regression path only
(used when no segment runtime is enabled). So production heads walk the actin they grabbed, not the
minifilament axis (INTEGRATION.md §1, physiological-baseline rule).

**What landed here (in scope: `segment_motor.py` · `segment_query.py` · `walk_dir` · `woven_cell.py` joint):**

- **`ac/weave/woven_cell.py` — the canonical joint.** New module function `actin_segment_topology(fiber_offsets,
  polarity)` and the `WovenCell.actin_segment_topology()` method flatten the woven fibers to per-bond segments
  `(seg_node_a, seg_node_b, seg_polarity)` — the exact three device arrays
  `MyosinForce.enable_segment_runtime` / `SegmentMotorRuntime` consume. `WovenCell.segment_barbed_directions()`
  is the host twin of the device `refresh_segment_barbed_kernel` (sign·(b−a) normalised), so the polarity gate
  runs on the CPU-only Mac.
- **`ac/cell/assemble.py` now delegates** its private `_actin_segment_topology` to the woven joint — single
  source of truth: the production cell binds the motor to the SAME topology the weave network exposes (the
  existing `build_cell` already calls `myo.enable_segment_runtime(...)`; this increment makes the topology it
  passes the canonical woven joint).
- **Gates — `tests/ac/weave/test_actin_motor_joint.py` (10 host + 1 CUDA-gated):**
  1. *topology* — flattened adjacent-node segments, one per bond, polarity broadcast per fiber; segment count
     = Σ(len_fiber − 1); no segment spans a fiber boundary; singleton-fiber / bad-polarity rejected.
  2. *barbed-end DIRECTION correctness* (the polarity gate) — `segment_barbed_directions` reproduces the device
     formula and is unit; on a woven cortex every segment's barbed dir projects **positively** on the
     node→barbed-end direction (correct sign, `min = 0.986`); on a straight fiber it equals the node-anchored
     `walk_dir` exactly (the segment-anchored and node-anchored polarity paths agree); a polarity flip flips
     every one of a fiber's segment directions (the reattach hand-off).
  3. *reaction symmetry + transmission* — a head bound to a woven segment closes Newton's 3rd law across the
     two nodes (`f_head + f_seg_a + f_seg_b = 0`, `< 1e-9`); a collinear isometric hold transmits the
     crossbridge tension 1:1 (Hill load == Bell load == `k_xb·|Δ|`), and the head is pulled toward the barbed
     end (directed contraction).
  4. *(CUDA-gated, native/gbook)* end-to-end — the woven joint feeds a `SegmentMotorRuntime`; an accepted
     commit binds a head whose device `walk_dir` equals the woven segment's barbed direction (the overwrite
     proven on device).
- **Scope honesty.** Host gates certify the JOINT geometry + the crossbridge oracle formulas (dev-Mac I0-A:
  CUDA kernels not launched here). The end-to-end device bind, GPU-residency, and the *emergent* contraction on
  the full native population are gbook-owned native gates (see §3 + below). No I0-B3 magnitude was chosen or
  tuned — `walk_dir` is seeded geometry, and the joint touches no `k_xb`/`f_stall`/`v0`/`kappa` GAP.

Full `ffn_sim/tests/ac/` suite: **412 passed, 43 skipped** (skips = CUDA-only native gates on the CPU Mac).

### Figures (I4-weave)

Regenerate all: `PYTHONPATH=<repo> python ffn_sim/scripts/ac_motor_vis.py` (one entry-point, extended).

- **`figs/i4_weave_joint.png`** — LEFT: a woven cortex patch (`weave_cell([CORTEX])`, 90 fibers → 540 segments)
  with each segment's barbed-end unit direction (the `walk_dir` a bound head is overwritten with) quivered at
  the segment midpoint — the arrows lie tangent to the shell, oriented toward each fiber's barbed end. RIGHT:
  the polarity gate — all 540 segments project positively (min 0.986) on the node→barbed-end direction, so the
  joint orients every segment correctly across the curved network. Geometry only (no magnitude); the on-device
  overwrite + emergent contraction are the native gates.

## Native-gate viz spec (lead runs on gbook)

Interactive 3-D cell-morphology **HTML** rendering the **engaged-head map** (per-head bound/free state colored
on the cortex) + **minifilament geometry** (backbone beads + the two anti-parallel head arms) on the LIVE
native cell; full-res / no downsampling; browser-verified (`browser_check.py`). See `INTEGRATION.md` for the
named fields and the full native-gate list (ensemble-stall emergence, GPU-residency zero-roundtrip, native
γ-floor report-not-tune).
