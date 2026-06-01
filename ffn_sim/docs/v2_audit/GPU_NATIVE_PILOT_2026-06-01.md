# GPU native myosin-ON pilot — pipeline validated + blow-up characterization (2026-06-01)

> Autonomous follow-on to the cupy port (#1/#1b). PI chose "pilot native run
> first". Goal: run a native (n_fil=1000) myosin-ON grip-walk run on GPU and
> see whether **adv** (bead re-targeting = KU-3.5 floor layer-3 transport
> signal) climbs. Outcome: **the GPU native pipeline is validated end-to-end**;
> the grip-walk *production config* (v0/force-scaling/dt/run-length to reach
> adv>0 *stably* at native scale) is iterative KU-3.5 science tuning, **PI-gated**.

## 1. Pipeline VALIDATED end-to-end (binned_r0 arm, GPU, n_fil=1000)

The `binned_r0` arm (proxy floor, no force-scaling) ran **8/8 samples clean on
GPU** at n_fil=1000 (`h3_ku35_gripwalk_tier1.py --device gpu`):
- constraint drift ~1e-14 every sample (the GPU constrained Action is correct
  at native scale),
- rigid + soft method-of-planes tension measured (so the `chains_tag_stacked`
  consumer fix works),
- adv (step_advances) tracked 1632→8551 (monotonic),
- soft-start equilibration prevented the int32 blow-up.

So **the whole GPU-resident stack — constrained BAOAB Action + M-SHAKE + Fixman
+ method-of-planes tension + soft-start equilibration — works end-to-end in the
real KU-3.5 driver on GPU at native scale.** This is the load-bearing result.

### Port integration bug found + fixed (commit 5147312)
The first pilot crashed in the rigid-bond tension M-OP:
`pos_byTag[chains]` indexed host numpy positions with the Action's
**device-resident cupy** `_chains_tag_stacked` → "Implicit conversion to a NumPy
array is not allowed". Fix: `chains_tag_stacked` property (host-numpy view,
mirrors `lambda_buf`); `prv_rnds` hardened likewise; consumer audit confirmed
these were the only leaks. CPU bit-invariance preserved (regression green).

## 2. grip-walk arm blow-up — characterized (NOT a port bug)

The `grip_walk` arm blew up at the FIRST production step in every config tried.
The int32 wrap-guard / M-SHAKE non-convergence **correctly caught** an
unphysical over-contraction — the guards work; the *config* is too hot at
native scale:

| config (n_fil=1000, grip_walk) | failure | worst |
|---|---|---|
| v0_accel=300, force-scaling ON, dt 1e-3 | int32 wrap | 2.342e8 |
| v0_accel=30,  force-scaling ON, dt 1e-3 | int32 wrap | **2.342e8** (identical → v0_accel-independent) |
| v0_accel=300, force-scaling OFF, dt 1e-3 | M-SHAKE drift 1.13e3 | — |
| v0_accel=10,  force-scaling ON, dt 2e-4 (Route-B CFL) | **STABLE — 4/4 clean, no blow-up** | drift 3e-15 |

**Two compounding "too hot" causes at native scale:**
1. **Force-scaling ×37.7 (Route B) at default dt** — the identical int32 worst
   value across v0_accel=300/30 proves the blow-up is the **stiffened head
   springs (k_head ×37.7) violating CFL at dt_factor=1e-3**. This is exactly the
   Route-B dt-tax the KU-3.5 record flagged (needs dt_factor≈2e-4, ~5× more
   steps).
2. **v0_accel=300 accelerated walk over-drives the larger network** — grip_walk
   with force-scaling OFF still blew (M-SHAKE) because the ×300 diagnostic
   accelerant over-contracts n_fil=1000 (it was tuned for n_fil~60–150). The
   `binned_r0` arm has no walk (s_grip≡0) → ran fine, isolating the walk as the
   second hot factor.

## 3. The real implication — the GPU unlock obsoletes the ×300 accelerant

The `v0_accel` ×300 was a *diagnostic* to make grip-walk observable in the SHORT
runs CPU allowed. **The GPU port's whole purpose is to make LONG runs feasible**
— so the principled native config is **literal (or near-literal) v0 + Route-B
force-scaling at dt_factor≈2e-4 + a long run** (run-length sized to ℓ₀/v0 ≈ 2.5 s
per bead, the KU-3.5 run-length finding), NOT the ×300 accelerant artifact. That
is a genuine multi-hour/day GPU **production** run and a multi-parameter science
choice (v0, force factor, dt, run-length, n_fil) the project reserves for PI
(KU-3.5 Tier-3 is gated on the §4b force-budget resolution).

## 4. Verdict + hand-off

- ✅ **GPU native myosin-ON pipeline is validated and stable** — both the
  binned_r0 floor (8/8) AND the **principled grip-walk production config**
  (force-scaling ON + dt_factor=2e-4 + v0_accel=10): 4/4 clean at n_fil=1000,
  drift 3e-15, **s_grip grows monotonically 1.2→2.7→4.4→6.0 nm** (the grip-walk
  transport half engages), engaged 51→160, no blow-up. The constrained-Action
  port + consumer fix + soft-start all work in the real KU-3.5 driver on GPU.
- 🎯 **adv=0 only because the run is short** — s_grip 6 nm ≪ ℓ₀, so heads have
  not yet walked a full bead (matches the KU-3.5 run-length finding exactly).
  At this rate (~0.15 nm / 1000 steps, v0_accel=10), reaching s_grip≈ℓ₀ to fire
  adv>0 needs ~10⁶ steps (≈ hours on GPU at ~90 steps/s) — **feasible on GPU,
  was infeasible on CPU**. That long run = the native KU-3.5 transport-regime
  measurement, the literal floor-layer-3 test. **This is the PI-gated production
  launch the whole port was built to enable.**
- ⚙️ **grip-walk production config is PI's call**: pick (force factor, v0,
  dt_factor, run-length, n_fil) for the native KU-3.5 transport-regime run. Data
  says: use dt_factor≈2e-4 with force-scaling (CFL), drop the ×300 accelerant in
  favour of a long run at realistic v0, and scale warmup with N.
- 🔧 Driver: `h3_ku35_gripwalk_tier1.py` gained `--device {cpu,gpu}`, soft-start
  warm-up, and `--arm {both,binned_r0,grip_walk}` (fast single-arm iteration).
- 🔴 PI-gated (unchanged): the long native production run itself; full-length
  GPU L_p/equipartition contract gates; default-device flip.

## 4b. Native production run attempt — Route-B stability ceiling at s_grip≈74 nm

PI picked: Route B ON, v0 accel-first, n_fil=3000, dt_factor=2e-4. Sequence:
- v0×50 short probe (n_fil=3000, dt 2e-4): **STABLE 4/4, drift 3e-15, s_grip→27 nm.**
- v0×100 production launch (20×60k=1.2M steps): **sample 1 reached s_grip=74 nm
  (drift 3e-15, clean), then M-SHAKE non-convergence (drift 1.87e3) in sample 2.**

So the native Route-B grip-walk run has a **stability ceiling around s_grip≈50–74 nm**
at dt_factor=2e-4, n_fil=3000 — well below ℓ₀=500 nm, so adv=0 was never crossed.
Notes:
- NOT a port bug — the constrained Action ran clean to 74 nm (drift 3e-15); the
  M-SHAKE guard correctly caught an over-driven backbone bond.
- v0_accel changes only how fast s_grip reaches the ceiling, not the ceiling
  itself (the scaled head force at a given s_grip is v0-independent) — so lowering
  v0 only delays the blow-up.
- The doc's earlier n_fil=150 run reached s_grip 160 nm → the ceiling is
  **density/scale-dependent** (n_fil=3000 ≫ 150), pointing at contraction-induced
  inter-filament LJ overlap at native density (not a simple spring-CFL, since the
  scaled-spring overdamped CFL is s_grip-independent).
- **Characterization probe running**: dt_factor=1e-4 (2× smaller), v0×100,
  n_fil=3000 — if s_grip clears 74 nm, dt is the lever (PI: smaller dt, ~2× more
  steps); if it blows at ≈74 nm again, it's contraction-induced overlap (needs a
  soft-core LJ / lower force-factor / lower density, not dt). **[result pending]**

**PI decision needed** for the native transport-regime run: how to lift the
Route-B native-scale stability ceiling — smaller dt vs soft-core/excluded-volume
handling during contraction vs lower force-factor vs lower n_fil. This is the
Route-B-at-native numerical-stability tuning the KU-3.5 record reserves for PI.

## 4c. Native 38k transport run (v0×1500, 2026-06-02) — TRANSPORT ENGAGED

PI ran the native cortex (n_fil=38,000 = 266k particles, force-scaling OFF,
dt 1e-4, v0×1500, dynamic-topology GSD) on the A5000. 9 samples then stopped.

| signal | result |
|---|---|
| **adv (bead re-targets)** | **0→0→0→0→17→53→81→100→131** — onset at sample 5; the myosin TRANSPORT half ENGAGES at native scale ✅ |
| s_grip | 70→265 nm then plateaus ~230–280 nm (walk↔re-target steady state) |
| drift | 3.4–3.6e-15 throughout — **stable** (force-scaling OFF passed s_grip 280nm; force-scaling ON had blown at ~74nm) |
| γ_total | ~2–5 ×10⁻³ mN/m, **flat, ~100× below the KU-3.5 band** — no rise |
| network forming | xlink_attach 8→31, myosin_attach 40→291 (captured now that GSD uses `dynamic=[property,topology]`) |
| GPU | 97% in production, ~25 steps/s, oscillating 0↔97% (per-sample O(N) host M-OP/snapshot) |

**Interpretation — exactly the predicted split:**
- ✅ **Transport regime reached at native scale on GPU** (adv climbs, stable) — the
  literal floor-layer-3 *transport* signal, infeasible on CPU, now demonstrated.
- ❌ **γ stays ~100× short** because v0×1500 accelerates the myosin WALK but NOT
  crosslink binding (literal kinetics) → the two timescales **desync** → the
  network never percolates (only 31 of ~2000 xlink heads bound) → no long-range
  stress transmission → no cortical tension. This **confirms** the coherence-half
  requirement: the true γ needs **literal (near-1×) timescales so crosslinking
  keeps pace with contraction** — a long literal run, not the accelerant.

**Bottom line:** v0_accel is a valid *transport diagnostic* (does myosin walk at
native scale — YES) but NOT a γ measurement (the coherence/percolation half is
broken by the walk/crosslink desync). The gold-standard γ run is literal-v0 +
long, on GPU.

Figures: `outputs/h3/figs/native38k_cell_structure.png`,
`native38k_bonds_by_type.png`, `native38k_v1500_trajectory.png`,
`native38k_v1500_structure.png`. Viz: `scripts/viz_gsd_structure.py` (3D GSD
render), `scripts/viz_ku35_trajectory.py` (per-sample s_grip/adv/γ).

## 5. Files

- `scripts/h3_ku35_gripwalk_tier1.py` — `--device`, `--arm`, soft-start warm-up.
- (port consumer fix + driver `--device`/soft-start already in commit 5147312.)
