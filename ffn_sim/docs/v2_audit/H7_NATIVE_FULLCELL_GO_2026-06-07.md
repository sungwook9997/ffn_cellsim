# H.7 native constrained integrator — FULL-CELL GO measurement (Lead, 2026-06-07)

**Task:** PI directed (A)+(B): wire the Stage-1b native constrained BAOAB plugin
(`cacb27b`+`0e2d490`+`24f12ad`+`398e96e`) into the production driver and measure
the **full-cell** end-to-end speedup + Gate-A/B feasibility on gbook (RTX A5000).
Driver: `scripts/h7_native_fullcell_go.py` (runtime integrator swap, additive).

## Foundation reproduced in-hand (gbook)
- Plugin builds; `import hoomd, _ffn_native` OK (hoomd must load first).
- `test_constrained_parity.py`: kT=0 max|nat−cupy| 3.7e-9 PASS; kT>0 500-step
  drift 9.8e-13 PASS.
- `bench_constrained.py`: integrator-only **39.47×** (native 5598 vs cupy 142
  steps/s @ 17,500-bead synthetic chains). Matches the sub-session's ~40×.

## Full-cell GO result — three criteria
Two-phase stable build (unconstrained warm-up with cytoplasm η → fresh
constrained production seeded from warm positions), n_fil=1000, all compartments
ON, pure-stepping timing (γ sampled outside the timed region):

| | cupy | native | speedup |
|---|---|---|---|
| steps/s | 82.8 | **201.4** | **2.43×** |
| 2×10⁸-step ETA | 671 h | **276 h (≈11.5 d)** | — |
| nonconverged_count | — | **0** | — |
| γ_soft (mN/m) | 2.96e-7 | 2.88e-7 | parity PASS |
| γ_rigid (mN/m) | 0.0835 | 0.0831 | parity PASS |

**GO criteria all PASS** (native available, nonconverged=0, γ parity) → the
native constrained integrator is **correct + stable end-to-end**.

## ⚠️ Reconciliation: 2.43× (this) vs 6.32× (sub-session) — a CONFIG difference, not error
- Pure-stepping timing did **not** change my number (2.47×→2.43×) → it is NOT
  γ-measurement overhead.
- The full physiological cell wires **3 `md.force.Custom`** every step —
  `EnclosedVolumePressure` (turgor), `MembraneSurfaceTension`,
  `NucleusConfinement` — plus the Xlink/Myosin binder updaters, and is
  compartment-particle-dominated (n_fil=120 → 12,840 particles, only 840 cortex).
- The sub-session's profiled build (398e96e) has its **own note: "0
  `md.force.Custom`"** — i.e. cortex+binders WITHOUT the compartment stack. That
  is why it reaches 734 steps/s / 6.32× / 3.2 days.
- The **physiological-baseline HARD rule** requires every compartment ON at its
  setpoint for production. So Gate-A/B MUST run on the full cell, where the real
  number is **2.43× / 201 steps/s / ≈11.5 days per 2×10⁸-step run** — NOT 3.2 days.

## Honest conclusion + remaining levers
- The integrator lane is solved+validated (39× microbench, 2.43× full-cell,
  parity, nonconverged 0). The **integrator is GO**.
- At the real operating point the step is dominated by, in rough order:
  the **3 compartment Custom forces** (turgor/membrane/nucleus — excluded from the
  sub-session's profile) + **binder host-sync** (myosin/xlink/integrin
  `cpu_local_snapshot`, ~584 µs in the no-compartment profile) + **LJ** (~606 µs).
- Next throughput levers, in impact order at the production operating point:
  1. GPU-resident port of the **compartment `md.force.Custom`** (turgor/membrane/
     nucleus) — the sub-session called Stage-2 "moot" only because their build had
     none; at the real operating point they are a dominant cost.
  2. **Binder host-sync** port (myosin/xlink/integrin off per-step
     `cpu_local_snapshot`) — the GPU-main binder lane.
  3. LJ force-stack (HOOMD nlist already optimal — research, not a quick win).
- **Gate-A/B feasibility:** ~11.5 days/run native (vs ~28 days cupy). Runnable but
  week-scale; the compartment-force + binder ports would cut it materially.

Artifacts: `outputs/h7/production/h7_native_fullcell_go.json` (gbook, Syncthing).
Driver committed; no core modules modified (runtime swap). Gate-B still needs the
separate relaxed-M-SHAKE mode (the rigid native updater cannot relax constraints).

---

## Compartment-force GPU port — MEASURED (corrects the attribution above)

Wired the existing device-resident compartment twins (`enclosed_volume_gpu` /
`membrane_surface_gpu` / `nucleus_confinement_gpu`) into the build opt-in
(`FFN_GPU_DEVICE_COMPARTMENTS=1`, post-attach swap in `cell.py`, default-off zero
regression) and measured the full-cell throughput+γ matrix
(`scripts/h7_compartment_gpu_port.py`, gbook, n_fil=1000):

| arm | steps/s | ×base | γ_soft (mN/m) | γ_rigid | nonconv |
|---|---|---|---|---|---|
| cupy + CPU-comp (baseline) | 84.0 | 1.00× | 2.957e-7 | 8.353e-2 | — |
| cupy + GPU-comp | 89.8 | **1.07×** | 2.957e-7 (EXACT) | 8.353e-2 (EXACT) | — |
| native + GPU-comp (full GPU-main) | 243.0 | **2.89×** | 2.879e-7 | 8.311e-2 | 0 |

**The compartment GPU port is CORRECT but MODEST (1.07× alone).** γ is bit-exact
vs the CPU forces → physics validated. **It corrects this doc's earlier ≈72%
attribution**: that was an inference from comparing to the sub-session's
no-compartment build (which differs in particle count + binder config, not just
the 3 Custom forces). The DIRECT measurement (same cell, CPU-comp vs GPU-comp)
shows the compartment Custom-force host-sync is only ~17% of the step (~860 µs).
The **native integrator is the dominant lever (2.7×)**; combined the full GPU-main
stack is **2.89× / 243 steps/s / ~9.5 days per 2×10⁸**.

⇒ The remaining wall is **binder host-sync (myosin/xlink/integrin) + LJ**, not the
compartment forces. Next ~2× lever = the **binder GPU-main port** (the GPU-main
binder lane). Gate-A/B is ~9.5 days/run on the full GPU-main stack today.

### Native ForceCompute upgrade (NativeRadialShellForce) — full GPU-main stack

Replaced the cupy compartment twins with the sub-session's native C++
`FFNRadialShellForce` (b659674/2a18203; `NativeRadialShellForce.from_{turgor,
membrane,nucleus}`, one ForceCompute + shared-mem block-reduce). gbook n_fil=1000:

| arm | steps/s | ×base | 2e8 ETA | γ_soft (mN/m) |
|---|---|---|---|---|
| cupy-int + CPU-comp (baseline) | 83.7 | 1.00× | 664h | 2.957e-7 |
| cupy-int + native-comp | 111.8 | 1.34× | 497h | 2.957e-7 (EXACT) |
| native-int + native-comp (full GPU-main) | **490.7** | **5.86×** | **113h ≈ 4.7 d** | 2.879e-7 |

Going cupy→native on the compartment forces lifted the native arm 243→491 steps/s
(**2×** — confirms the cupy gpu_local_snapshot launch-bound floor). γ bit-exact.
**Full GPU-main stack (native integrator + native compartment ForceCompute) =
5.86× / 491 steps/s / ~4.7 days per 2×10⁸** at the real production operating point
(all binders + compartments ON). The sub-session's 8.3×/3.3d is a lighter-binder
config; the residual gap to it is the **binder host-sync** (myosin/xlink/integrin
`cpu_local_snapshot`) — the next ~2× lever (the GPU-main binder lane), which would
bring Gate-A/B toward ~2.7 d. Below that, LJ is HOOMD-optimal (Tree nlist tested);
sub-µs would need mesh-restricted excluded-volume (manifold lane, research).

**Gate-A/B feasibility TODAY: ~4.7 days/run on the full GPU-main stack** (vs ~28 d
cupy baseline). Gate-B still needs the separate relaxed-M-SHAKE mode.
