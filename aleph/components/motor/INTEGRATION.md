# ac/motor — INTEGRATION notes (I3 head-resolved NMII → lead-owned serial native integration)

Track **motor-nmii** (`ac/motor-nmii`), increment **I3**. This file is the contract between this worktree's
NET-NEW `ac/motor/` package and the shared `ff/` runtime the **lead** edits in spine order (§1.3). Nothing
here is edited into `ff/` from this worktree.

Contents: (1) the frozen `hand.py` device API downstream tracks consume; (2) the `ff/` patch-notes (exact
file · symbol · before→after · double-count guard) the lead applies WHEN I3 lands native; (3) the native-gate
spec (gbook, lead-run); (4) honest-scope + forbidden-shortcut notes.

---

## 1. Frozen interface this track OWNS (§1.4)

- **`ac/motor/hand.py` — the per-head hand/KMC device API. I3 OWNS it.** `NMIIHandParams` (device struct) +
  `bell_off_rate` / `hill_velocity` / `attach_prob` / `detach_prob` (`@wp.func`) + `attach_kernel` /
  `step_detach_kernel` (device KMC sub-steps, device RNG, no host draw) + `allocate_hand_state`.
  Downstream **I4 (weave) / I6 (adhesion) / I8 (mt-coupling)** CONSTRUCT hands via this API in their OWN
  `ac/<pkg>/presets.py`; they **never edit this KMC core**. (I8 dynein is a POPULATION of these hands with a
  dynein `NMIIHandParams` set — engaged count emerges; never a lumped `N×f`.)
  - ⭐ **`walk_dir` hand-off (NEW, HARD for downstream).** `allocate_hand_state` now also carries a per-head
    `walk_dir` (vec3d) = the bound actin filament's **barbed-end polarity**, along which the power stroke
    advances the crossbridge attachment point. `MyosinForce` seeds a bipolar-geometry default (outward along
    the backbone axis, from `build_minifilament_nodes`). **I4-weave MUST overwrite `state["walk_dir"][h]` from
    the true actin polarity when a head (re)binds** (in its `attach` preset) — otherwise heads walk along the
    minifilament axis rather than along the actin they grabbed. A ZERO `walk_dir` = a passive head (no directed
    force) — the OFF/regression default only; production must set it (physiological-baseline rule).
  - ⭐ **rest-length split (NEW).** `NMIIHandParams` now has TWO rests: `r0_head` = the head↔BACKBONE **arm**
    rest (= topology `head_offset_um`), and `r0_xb` = the head↔ACTIN **crossbridge** rest (~0). They were
    previously one overloaded field; downstream presets that build `NMIIHandParams` must set BOTH (see
    params_i0b3.yaml — `r0_xb` is grounded at 0, not a GAP).
- **`ac/motor/minifilament_warp.py :: MyosinForce`** — the inner-force primitive per §1.4:
  `MyosinForce.accumulate(pos, out_force)` adds the minifilament's backbone + head↔backbone-arm + **power-stroke
  crossbridge** forces to the global node force. The **lead-owned integrator** sums `MyosinForce` +
  `StericForce` (I2b) + `PressureCoupling` (I1a). Loop order per outer tick: mechanical solve →
  `MyosinForce.compute_loads(pos)` → `MyosinForce.step_kinetics(...)` (KMC) — so the per-head load fed to
  Hill/Bell is the LIVE **tangential** crossbridge tension. ⭐ **The crossbridge is now the ACTIVE term**: its
  attachment point on actin ADVANCES along `walk_dir` by the walked `abscissa` (the myosin power stroke), so a
  stepping head makes a **directed contractile force** and force-velocity self-limits at `f_stall` (the load
  rises with abscissa until `hill_velocity → 0`). The engaged fraction / ensemble stall EMERGE (never imposed).
  Constructor signature gained a `walk_dir` arg (default geometry from `build_minifilament_nodes`, which now
  also returns `walk_dir`).

The host-numpy analytic oracles (`hill_fv_analytic`, `bell_kinetics_analytic`, `ensemble_stall_analytic`,
`minifilament_topology`) carry the SAME closed forms as the device `@wp.func`s, so the native gate reproduces
the CPU-green analytic gates bit-for-formula.

---

## 2. `ff/` patch-notes — lead applies these WHEN I3 lands native (do NOT apply from this worktree)

⚠ **Double-count guard #6 (`f_myo` + fine-motor), HARD:** the lumped myosin below and the head-resolved
`MyosinForce` are the SAME physical force at two fidelities. They must **never both be active**. The deletions
here land in the **SAME commit** as `MyosinForce` is wired into the native driver (the "same-commit-as-fine-
motor guard", build-plan §6 line 545). Running both = counting myosin contractility twice.

### 2a. `ff/network_warp.py` — DELETE the lumped constant-`f_myo` myosin

| Symbol | File:line (current) | Before | After (I3) |
|---|---|---|---|
| `myosin_kernel` | `ff/network_warp.py:105-118` | constant-`f_myo` dipole `f=(f_myo/L)·d` pulling `links[t]=(i,j)` together | **DELETE as production.** Replaced by `MyosinForce.accumulate` (explicit backbone beads + heads + per-head crossbridge). May survive ONLY as a Warp diagnostic clearly named non-production. |
| `myosin_bound_kernel` | `ff/network_warp.py:170-197` | `f_myo` with a bound-mask + `p0` KMC | **DELETE.** The bound-mask/KMC is now per-head device state in `ac/motor/hand.py`. |
| `myosin_force_np` | `ff/network_warp.py:502-510` | host wrapper launching `myosin_kernel` | **DELETE** (host wrapper of a deleted kernel). |
| `f_myo` scalar arg | threaded through `simulate_*` drivers (`:672, 950, 1234, 1327, 1385, 1492/1547, 1561/1652, 1665/1768, 1803/1874/1885`) + `relax_on_device(myo_links, f_myo)` | a single swept contractile magnitude | **REMOVE** the `f_myo` scalar + `myo_links`/`cortex.myo_i,myo_j` two-anchor path from the native driver; the motor is now `MyosinForce` (its own device state). Diagnostic drivers may keep a `f_myo` path if renamed non-production. |
| `cortex.myo_i` / `cortex.myo_j` | two-anchor myosin link index arrays | one contractile link per placed minifilament (aggregate) | **REPLACE** with explicit minifilament topology (`build_minifilament_nodes`: backbone beads + 2·N_side heads + dynamic crossbridges). |

### 2b. `ff/gamma_floor.py` — DELETE the lumped `_myosin_force`

| Symbol | File:line | Before | After (I3) |
|---|---|---|---|
| `_myosin_force` | `ff/gamma_floor.py:251-260` | host accumulate `f = f_myo·u` (constant per link) | **DELETE as production.** The γ-floor harness reads the head-resolved contractile force from `MyosinForce` instead. |
| `external_force` | `ff/gamma_floor.py:302-313` | calls `_myosin_force(...)` | rewire to sum `MyosinForce.accumulate` (lead-owned integrator). |
| `equilibrate` / `measure_gamma` `f_myo` arg | `ff/gamma_floor.py:320, 419` | sweep `f_myo` as the controlled prestress | the controlled variable is now the head-resolved minifilament (its I0-B3 params); the emergent γ is MEASURED from it. |
| `PROD_N_MYO=442` + Nie-2015 density | `ff/gamma_floor.py:116-118` | 0.625 minifil/µm² (the only direct cortical NMII datum) | KEEP as the native minifilament COUNT; the native γ-floor is **DENSITY-floored** — see gate 3d. |

### 2c. `ff/myosin_linear.py` — DEMOTE to diagnostic

`myosin_linear.minifilament_kernel` (aggregate two-anchor inverted-linear-FV) is **NOT the production motor**
(build-plan §2 REBUILD line 159, I3 line 243). It may remain a Warp **diagnostic/control** only — e.g. the
κ→∞ linear-limit control for the Hill FV native gate — never wired as the runtime motor. Its `V0/F_HEAD/
N_HEADS/DUTY` module constants are the archived defaults, superseded by the I0-B3 ledger (do NOT re-import
them as chosen values).

### 2d. `ff/hand_kmc.py` — REUSE read-only, port to device

`ac/motor/hand.py` ports the NF2007 attach/step/detach semantics + Bell slip law from `ff/hand_kmc.py`
(host) to device-resident per-head state. `ff/hand_kmc.py` stays untouched; the `NMIIA_MYOSIN` preset there is
the parameter LINEAGE (k_off0=0.35, x_β=0.6 nm, k_on=50) — but I0-B3 re-audits those for head-level use
(several are GAP; see params_i0b3.yaml).

---

## 3. Native-gate spec (lead runs on the gbook A5000, in spine order I3 after I2b)

The dev Mac cannot run Warp-CUDA (CPU-only build) — these gates are the lead's. Each maps to a CPU-green
analytic oracle in this package.

- **3a. Topology/count invariants (native).** On the full-population cell, assert per minifilament
  `n_particles = n_bb + 2·N_side`, `n_static_bonds = (n_bb-1) + 2·N_side`, exactly N_side heads per side,
  backbone a connected chain. Log unique-active minifilament IDs + head/backbone/bond totals + dormant
  allocation + peak GPU bytes (native-population ledger, HARD). Oracle: `minifilament_topology`.
- **3b. Single-head Hill FV vs Hill-1938 (native).** A single device minifilament under a controlled load
  ramp reproduces the host `hill_velocity` (chosen κ) to tolerance; the κ→∞ run matches the linear control.
  Oracle: `hill_fv_analytic`. Force-free at v0; zero velocity at F_stall.
- **3c. Bell slip + per-head Newton closure (native).** Device `bell_off_rate` monotone in load == host;
  the inner mechanical solve's per-head force residual → 0 (Newton's-3rd-law crossbridge). Oracle:
  `bell_kinetics_analytic`, `minifilament_topology.head_newton_residual`.
- **3d. Ensemble stall EMERGES + native γ-floor (report-not-tune, HARD).** On the native cortex the
  per-minifilament stall must EMERGE from the Bell-governed bound-head population (matches
  `ensemble_stall_meanfield` within stochastic spread), NOT be imposed as `N_side·F_head`. The resulting
  native γ is **density-floored → a FINDING** (§6.2): report the deviation from the target band honestly;
  **NEVER add heads / raise density / stiffen k_xb to close the floor.** (The historical ~530× γ deficit is
  force-magnitude/density-bound; the head-resolved motor does not by itself close it — that is the expected,
  reported outcome.)
- **3e. ATP/step work sign (native).** Per-head mechanical work `F·d_step ≥ 0` below stall, efficiency
  `η=F·d/ΔG_ATP ∈ [0,1)`. Oracle: `hill_fv_analytic.step_work` (needs `d_step`, `dG_atp` — I0-B3 GAP).
- **3f. GPU-residency zero-roundtrip (HARD, I0-A).** A profiler gate shows ZERO authoritative GPU→CPU
  roundtrips for myosin state inside the outer physical-time loop: `bound`/`anchor`/`abscissa`/`walk_dir`/
  `loads` stay device-resident; only the scalar RNG seed + tick index cross from host. No per-step host
  snapshot/mutation. (`walk_dir` is refreshed on-device by the I4-weave attach preset from actin polarity.)
- **3h. Power-stroke → force coupling (NEW, native).** The whole point of I3: a single engaged device
  minifilament under an isometric hold must build a DIRECTED contractile crossbridge force that self-limits at
  `f_stall` as the walked abscissa lifts the tangential load (`hill_velocity → 0`); zero abscissa ⇒ zero active
  force. Oracle: `powerstroke_analytic.{crossbridge_force,step_engaged_head}` — the device kernels are already
  verified bit-identical to it on CPU codegen. The settled strain `f_stall/k_xb` must equal
  `minifilament_topology.working_stroke_strain` (5–20 nm physical band).
- **3g. OFF regression (bit-identical).** With no bound heads (or `MyosinForce` off), the native step is
  bit-identical to the passive cell — the double-count guard's regression anchor.

### Native viz (lead, HTML — name the fields):
Interactive 3-D cell-morphology **HTML**: (i) **engaged-head map** — per-head bound/free state colored on the
cortex (shows the emergent engaged fraction spatially); (ii) **minifilament geometry** — backbone beads + the
two anti-parallel head arms + live crossbridges to actin. Full-res, no downsampling (PI 2026-07-07);
browser-verified via `browser_check.py` (PI 2026-07-07/08). ANOMALY → STOP, render, surface to PI (never
self-correct / never "explosion").

---

## 4. Honest-scope + forbidden-shortcut notes (state these; do not overclaim)

- **⭐ DEFECT FIXED 2026-07-16 (power-stroke ↔ force decoupling).** The first I3 source advanced the head
  `abscissa` in `step_detach_kernel` but the crossbridge/load kernels pulled the head toward a FIXED anchor and
  never READ the abscissa — so the Warp source realized only a passive spring network + KMC bookkeeping and the
  ACTIVE directed contraction (the whole point of P2) was NOT generated; the host oracles missed it because the
  kernels are not launched on the Mac. FIX: the crossbridge attachment point now advances along `walk_dir` by
  the walked `abscissa` (`x_att = anchor + abscissa*walk_dir`), so a stepping head makes a directed contractile
  force and force-velocity emerges (self-limits at `f_stall`). Also split the overloaded `r0_head` into the
  arm rest (`r0_head`) and the crossbridge rest (`r0_xb`~0) — the old overload pre-stressed every freshly-bound
  head by `k_xb·0.2 µm`. New LOCAL gates: `tests/ac/motor/test_powerstroke_coupling_oracle.py` (numpy
  behavioural) + `test_powerstroke_kernel_wiring.py` (AST gate that fails if the kernel is re-decoupled + a
  Warp codegen check). No `ff/` file touched; no I0-B magnitude tuned (`k_xb`/`f_stall`/`v0`/`kappa` stay GAP,
  `r0_xb`=0 is grounded, `walk_dir` is geometry data not a magnitude).
- **P2 honest scope:** at the quasi-static inner equilibrium `v_slide → 0`, so per-head force → F_stall (a
  constant) and **force-velocity bites only the transient**. The fine-grained content that SURVIVES at the
  settled state is (1) the DERIVED contractile magnitude, (2) the **Bell-emergent engaged fraction** (fewer
  heads bound under load — self-limiting), (3) **topology reformation**. Do NOT claim FV changes the settled
  state. (Build-plan §1 P2 line 113-115.)
- **⚠ `mesoscale_force_scaling` is FORBIDDEN in native.** The archived `cortex/myosin.py:305-347` scaled
  `k_head`, `k_xb`, `F_stall` by `(native_motor_count / config_motor_count)` to run few coarse motors as
  many. The new engine runs **explicit native minifilaments at full physiological population** (PI HARD:
  validate at full native population; never lower density / coarse-grain to fit memory). If native motor+head
  state exceeds 16 GB, optimize device layout or use larger/multi-GPU hardware — do NOT re-enable mesoscale
  scaling.
- **⚠ `k_xb` MASTER force knob (params_i0b3.yaml).** Physical crossbridge stiffness 100–1000 pN/µm (working
  stroke 5–20 nm). The archived cortex `k_head_actin = 1 pN/µm` is the broken-soft AFINES surrogate (strain
  ~2 µm ≫ the ~0.3 µm minifilament → stall geometry collapses; `minifilament_topology.working_stroke_strain`
  demonstrates it). GAP to PI; NEVER tuned to a γ/force band.
- **⚠ FV-form conflict (params_i0b3.yaml::kappa_hill).** THREE positions: Hill-1938 muscle κ=0.25 (I0-A
  mandate) / Kovács-2003 archived NMII κ=0.5 (the value the archived head-resolved code actually used, so a
  non-muscle Hill curvature IS sourced) / PI-2026-07-07 non-muscle-LINEAR (κ→∞). The Hill machinery is built
  (κ is a device param); the CURVATURE is a GAP for PI to reconcile — NOT chosen here.
- **⚠ engaged-fraction vs duty (params_i0b3.yaml::k_on/duty).** The k_on/k_off equilibrium engaged fraction
  (~0.99 at the archived k_on=50/s) differs from the measured Kovács unloaded duty ~0.1. Flag for PI: is the
  head-level `k_on` right, or is the duty a different (ATPase-cycle) quantity? The engaged fraction EMERGES
  from k_on/k_off; `duty` is a cross-check only, never a runtime input.
