---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# FF crawl (C1+C2) — native protrusion works, crawl mechanism validated, native SPEED diagnosed (2026-07-09)

Overnight autonomous run (PI /goal 8h). Goal: drive the validated adherent cell through its motility program and
VISUALIZE the development. C1 (protrusion) and the C2 crawl mechanism are validated; the native crawl SPEED is
grid-drag-limited — a diagnosed open item. This records what works and the four diagnostic steps, so it is not
re-derived.

## What works (validated)

- **C1 — leading-edge protrusion in the native path.** `gpu_force_fn` (the ONLY native large-dt path) never applied
  the protrusion → native disp≡0. Added the `elif protrude` branch (mirrors `full_force`/the explicit loop). NATIVE
  run (A5000, Nc=266,000, `--from-resting --microtubules --implicit`): reproduces the resting checkpoint (ΔP=40 Pa →
  γ=0.151 mN/m, V/V0=1.000), adheres flat (basal-gap 0.007 µm, 12475 contacts), and the native displacement is now
  non-zero — the disp=0 blocker is fixed. Committed `14e2cc5`.
- **C2 — the crawl is EMERGENT from protrusion + clutch turnover, and traction-driven.** No explicit front-rear
  "treadmill" bias is needed: with the clutch catch-slip turnover firing (`--kmc-every` small enough), the
  protrusion-directed (phat) motion emerges. **Traction-driven, decisively**: the clutches-OFF audit gives net
  displacement ≈ 0 (protrusion alone is internal, a Newton pair) while clutches-ON crawls — audit ratio 1.1e4–8.1e4×
  PASS at both native and coarse scale. At **coarse resolution the crawl is PHYSIOLOGICAL: ~60 nm/s** (1.8 µm in 30 s;
  reality band 10–100 nm/s, KU-3.12). **SEED-ROBUST** (not a single-seed artifact): seeds 7/11/17/23 all crawl
  FORWARD (+phat) at 28–45 nm/s (28.6 / 36.3 / 32.0 / 44.9), tight cluster, mean ≈ 35 nm/s.

## The native crawl SPEED — four diagnostic steps (root cause: grid-dependent drag)

The native crawl runs but is ~500× too slow (0.12 nm/s vs the coarse 60 nm/s). Diagnosed by elimination:

1. **kmc_every ≥ steps froze the clutches.** The default `kmc_every=2000` (tuned for the explicit tiny-dt path) ≥ the
   implicit-run step count → the catch-slip turnover + nascent rebind block NEVER fired → clutches stayed all-bound
   (bound=1.0) → no crawl. Fixed: `--kmc-every` (set so dt·kmc_every ≈ 1–2 s resolves the ~1 s clutch lifetime).
2. **Explicit front-bias "treadmill" — REJECTED.** Making nascent rebind fire only under the leading edge OVER-de-adheres:
   as the COM advances, formerly-front clutches fall behind it and stop re-forming → bound→0 → traction collapses →
   LESS crawl (46 vs 61 nm/s, bound 0.00). The emergent uniform-turnover path is both simpler and better. Kept as a
   documented dead-end (`--treadmill`, off by default).
3. **Clutch DENSITY — REFUTED.** Hypothesis: 10953 per-node clutches at native dilute per-clutch load (0.088 pN ≪ F*≈7 pN)
   → no catch-slip rupture → symmetric pinning. Tested with `--n-fa` (cap to N discrete FA sites via farthest-point
   sampling): capping made the crawl SLOWER, not faster (dense 40 sites → 7 nm/s vs sparse 44 → 60 nm/s at equal count),
   and sparse-vs-dense cortex differ 5× at equal clutch count → density is not the driver. (`--n-fa` kept as a
   physically-reasonable discrete-FA option, but it is not the speed lever.)
4. **ROOT CAUSE — grid-dependent drag Σγ ∝ Nc.** `physical_node_gammas` assigns the NF2007 **single-fiber log-drag per
   node**, so the whole-cell COM drag Σγ scales with the node count. At native Nc it is ~100–300× the physical
   whole-cell Stokes drag 6πηR → v = F_net/Σγ ∝ 1/Nc. Confirmed: v = 60 → 13 → 0.12 nm/s as Nc = 840 → 3500 → 266,000.
   The **physically-correct drag is Σγ_cortex = 6πηR** (η=65.9 Pa·s — the SAME bulk-η calibration validated for the AFM
   indenter, `FF_RESULTS_LOG` 7521647). **⚠️ But imposing it destabilizes the implicit solver**: the per-node
   γ=6πηR/Nc ≈ 0.035 makes the solver diagonal γ/dt ≪ the crosslink stiffness K (~1e6) → the cortex rigid-body modes go
   unregularized → NaN. So grid-consistent crawl drag is an **OPEN solver-side item** (regularize the rigid modes, or
   add an inertial/mass term so the COM drag is decoupled from the per-node stiffness scaling) — not a safe one-liner.

   **Two fixes attempted, BOTH fail** (⇒ the drag must be fixed INSIDE the solver): **`--bulk-drag`** (in-solve
   Σγ=6πηR) → per-node γ/dt ≪ K → rigid modes unregularized → **NaN** (coarse fil-120 + fil-500 diverge at ~T=2.5 s).
   **`--com-drag`** (post-hoc: keep per-node γ, override the COM translation to 6πηR after each step) → **RUNAWAY**
   (767 µm/30 s, clutches rip off): overriding the COM post-solve desynchronises it from the implicit clutch springs
   → the shift stretches the clutches → the measured external force grows → the next shift grows = positive feedback.
   The grid-consistent drag cannot be patched after the solve while position-dependent external forces (clutches) are
   in the loop — it must be solved WITH the clutch coupling (an in-solver rigid-mode-regularized whole-cell drag).
   Both flags are left in as documented dead-ends (off by default).

## Deliverable + recommendation

- **Visualised**: the native adhered+protruding cell (`ff_dev_S1_native_morph.html`, browser-verified) + a coarse-scale
  VISIBLE crawl demo (`ff_dev_crawl_demo_*`) showing the validated mechanism producing physiological directed motility
  (frame-animated, COM trajectory trace) + dev-curves (COM-disp / polarization / bound-FA / σ_vm vs time, OFF-audit
  overlay). The native full-res interactive HTMLs are kept local (>100 MB, GitHub limit); screenshots + curves are the
  committed record.
- **For PI**: the motility machinery (protrusion + clutch turnover) is mechanistically complete and validated
  (traction-driven, physiological at coarse scale). The single blocker to a native-scale physiological crawl is the
  **grid-consistent whole-cell drag under the implicit solver** — a focused numerics task (rigid-mode regularization /
  inertial term), the natural next step. No magic numbers were introduced; F* stays 7 pN as-recorded; no gate loosened.

New crawl-driver knobs (all additive, default = prior behaviour): `--kmc-every` (turnover cadence — REQUIRED < steps for
the implicit path), `--seed` (cross-seed ensemble), `--n-fa` (discrete FA sites), `--treadmill` (rejected front-bias,
documented), `--bulk-drag` (grid-invariant drag diagnosis, destabilizing — off).
