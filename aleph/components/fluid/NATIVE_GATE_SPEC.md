# ac/fluid — NATIVE-GATE SPEC (lead runs on the gbook A5000, in spine order)

The fluid-spine track hands the lead a **CPU-green** module set (the NumPy oracles below all pass on the
dev Mac) + this spec. The dev Mac cannot run Warp-CUDA (I0-A), so every gate that needs the *device* is
listed here for the lead to run on the A5000, **in spine order** (I1a → I1b → I1c; each `--from-resting`
gate literally contains the previous increment's live state). Precondition for opening a device slot:
`python -m pytest aleph/tests/ac/fluid/ -q` is green (65 tests, no CUDA).

Legend: **[STRUCTURAL]** = regression/conservation (no physiological magnitude needed, runnable now);
**[MAGNITUDE — blocked]** = a verdict HELD until its I0-B parameter closes (do NOT tune to pass).

---

## NG-0 — Warp↔NumPy source faithfulness (I1a, run first)
On a small case (e.g. `16³`, random p, no-flux box), launch `biot_pmass_update_kernel` on CUDA and compare
to `fv_reference.BiotFVReference.step` to **round-off**. Same for `rad_transport_kernel` vs
`transport_reference.RADTransportReference.step` and `darcy_discharge_kernel` vs `darcy_analytic.darcy_flux`.
This certifies the device kernel IS the CPU-verified stencil. **[STRUCTURAL]**

## NG-1 — FSI-OFF regression == Warp-FF bit-identical (I1a)
With the fluid channel forced OFF, the `ac/` cell reproduces the frozen Warp-FF result **bit-identically**
(the `biot_fsi=None` regression oracle, PN-5). Any drift = an unintended coupling. **[STRUCTURAL]**

## NG-2 — device conservation suite (I1a)
On the live membrane-minus-nucleus `Domain`:
- **impermeable limit** (all membrane flux 0): total content `S∫p dV` conserved to round-off over N steps;
- **content == flux**: `d/dt(total content)` equals the integrated membrane hydraulic flux (`boundary`) +
  `∫(s_water − α∇·v_s)dV` to machine precision (discrete Gauss);
- **constant-state**: a uniform `p=Π₀` with zero source + frozen solid does not drift (incl. under the
  moving-domain remap — the `moving_face` accounting must close);
- **Terzaghi / Green** on a device box match the analytic oracles (the `fv_reference` CPU gates, on CUDA). **[STRUCTURAL]**

## NG-3 — resting Π₀ = 40 Pa pre-tension at t0 (I1a)
Seed `p_bar = Π₀ = 40 Pa` (I0-B1). At `t0` the shell is pre-tensioned to `γ = ΔP·R/2` (the FF baseline is
preserved), NOT a floppy unpressurised bag. Physiological-baseline HARD rule. **[MAGNITUDE — Π₀ is
PI-ratified-as-proxy; MCF7 datum is the expiry trigger]**

## NG-4 — drained↔undrained rate dependence + spatial τ_p ≈ 1 s (I1a, the physics payload)
Ramp a load fast vs slow: the undrained (fast) response is STIFFER than the drained (slow) one (Skempton).
The spatial poroelastic relaxation time reads `τ_p = R²/c_v ≈ (7.5 µm)²/(50 µm²/s) ≈ 1.1 s` at `R=7.5 µm`
(I0-B1 `c_v` anchor). A single scalar-turgor balloon CANNOT show either — this is the fluid-first payload.
**[MAGNITUDE — c_v is the KB-3.B3.2 anchor (cross-cell); MCF7-specific D_p is the expiry trigger]**

## NG-5 — adjoint transfer + closed-system net-force + nucleus no-flux (I1a)
On the native cell: the fluid→solid `PressureCoupling` transfer does non-negative dissipation (work sign),
the closed system shows **no COM translation** (`Σ internal force ≈ 0`, `test_i1a_ibm_gates` on CUDA), and
`NucleusNoFluxBC.leaked_flux()` == 0 (relative no-flux across the moving envelope). **[STRUCTURAL]**

## NG-6 — device-residency zero-roundtrip (I1a)
A profiler trace shows **zero authoritative GPU→CPU roundtrips inside the physical-time loop** (I0-A). The
only host reads are out-of-hot-loop diagnostics (`total_content`, `leaked_flux`, reports). **[STRUCTURAL]**

## NG-7 — Darcy kinematics on device (I1b)
`v_f − v_s == q/φ` cell-by-cell; a slab's discharge is linear in Δp and in `k/µ`; a passive tracer follows
`v_f`, NOT the discharge `q` (advection-velocity gate). **[MAGNITUDE — blocked on φ (I0-B1b GAP)]**

## NG-8 — conservative monomer transport on device (I1c)
`∫φc dV + N_polymer` conserved to **machine precision** through barbed consumption + pointed release (the
non-vacuous invariant); FRAP recovery time reads back `D_c`; the advection front follows `v_f`; `c ≥ 0`
(positivity). **[MAGNITUDE — blocked on c_0 (I0-B1c GAP, MCF7 free pool) + φ]**

## NG-9 — full-native wall-time benchmark
Benchmark ONE full-config `--from-resting` FSI-ON outer step at the native population (70,686 cortical
F-actin + compartments) BEFORE any conclusion — the single biggest unmeasured number (hard-truth #5).
Publish the §A.1 ledger: `N_unique_active`, `N_allocated`, `N_nodes`, `N_state`, **peak `GPU_bytes`**, and
`Δt_outer`/inner-iteration counts. Never lower a biological density to fit memory. **[STRUCTURAL — measurement]**

## NG-10 — CFD-ON cell-morphology visualization (viz obligation, §1.9)
Interactive **3-D CELL-MORPHOLOGY HTML** on the REAL geometry (peel/slab/cut views, **full-res, no
downsampling**, browser-verified via `browser_check.py` — a WebGL grep is meaningless). Render, on the live
membrane-minus-nucleus domain:
- the spatial **pore-pressure** field `p_excess` [Pa] (color on a slab cut through the cell);
- the reconstructed **pore-fluid velocity** `|v_f|` [µm/s] + `q` **discharge streamlines** (the fluid FLOWS);
- the membrane hydraulic-flux hotspots + the nucleus relative-no-flux shell.
An anomalous field → STOP, render, surface to PI; never self-correct or wave it off (§1.9 ANOMALY).

---

## I0-B GAPS surfaced to PI (do NOT choose — magnitude gates INVALID until closed)

| Param | Increment | Status | Blocks | File |
|---|---|---|---|---|
| **φ** (porosity) | I1b/I1c | **GAP — PI** (water fraction ~0.7 is an UPPER proxy, not the Darcy pore fraction) | NG-7, NG-8 native | `params_i0b1b.yaml`, `params_i0b1c.yaml` |
| **c₀** (initial free G-actin, MCF7) | I1c | **GAP — PI** (total ~100 µM; free ~tens µM, MCF7-specific) | NG-8 native | `params_i0b1c.yaml` |
| **D_c** (monomer diffusivity) | I1c | **draft/Medium** (~2-6 µm²/s crowded; ratification + MCF7 applicability OPEN) | NG-8 magnitude | `params_i0b1c.yaml` |
| k (permeability), µ_pore | I1b | draft/Medium (I0-B1, consistent with c_v anchor) | NG-4/NG-7 tighten | `params_i0b1.yaml` |
| L_p (membrane hydraulic cond.) | I1a | draft/Medium (MCF7/AQP5; unit reconciliation open) | NG-4 tighten | `params_i0b1.yaml` |

Already CLOSED for I1a native: `c_v` (anchor), `α=1.0` (ratified), `M/S` (derived), `Π₀=40 Pa` (proxy).
The **analytic gates are parameter-agnostic** where a GAP exists (φ, c₀ swept as oracle parameters), so the
CPU suite is fully green now; only the NATIVE magnitude verdicts wait on PI.
