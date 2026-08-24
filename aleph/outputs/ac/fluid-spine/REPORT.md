# fluid-spine track — REPORT (I1a + I1b + I1c authored, CPU-green)

Track: **fluid-spine** (I1a → I1b → I1c) · branch `ac/fluid-spine` · base `ac/new-engine`.
Status: the conservative fluid-first Biot substrate (I1a), the mixed q-p pore-fluid-velocity form (I1b), and
the conservative G-actin monomer transport (I1c) are all **authored and CPU-green**: 65 host-numpy tests
pass (no Warp/CUDA), the Warp-CUDA kernel source type-checks (codegen-clean, device-guarded to CUDA), and the
native gates + INTEGRATION patch-notes + I0-B gaps are handed to the lead. This REPORT + `figs/` is the **ac/
visualization pattern** the other tracks copy (AC_PARALLEL_SESSIONS §1.9).

## What landed

| Piece | Files | Gate (CPU, local) |
|---|---|---|
| **I1a** conservative p/mass Biot | Warp: `field_grid.py` `biot_substrate.py` `boundary.py` `domain.py` (OWNS `Domain`+`set_nucleus_boundary`) `scheduler.py`; NumPy oracle: `fv_reference.py` `ibm_reference.py` | 21 tests: constant-state · impermeable mass conservation · content==flux · pressure-work sign (SPD) · CFL · Terzaghi/Green 2nd-order convergence · Peskin adjoint `Interp=h³Spreadᵀ` · net-force projection |
| **I1b** mixed q-p / `v_f` | Warp: `velocity.py`; oracle: `darcy_analytic.py` | 12 tests: flux linear in Δp + mobility · `v_f−v_s=q/φ` (φ swept) · hydrostatic zero-discharge · radial conservation · no relative flow at equilibrium |
| **I1c** conservative RAD transport | Warp: `transport.py`; oracle: `transport_analytic.py` `transport_reference.py` | 7 tests: total-actin `∫φc+N_polymer` conserved to **machine precision** · positivity · FRAP↔D_c · advection follows `v_f` not `q` · uniform-state · treadmill balance |
| Params | `params_i0b1.yaml` (I1a, closed) `params_i0b1b.yaml` (φ GAP) `params_i0b1c.yaml` (D_c draft, c₀ GAP) | I0-B ledger |
| Hand-off | `INTEGRATION.md` (6 `ff/` retirement patch-notes) · `NATIVE_GATE_SPEC.md` (NG-0…NG-10, spine order) | lead applies |

**Why analytic-first + a NumPy discrete reference:** the dev Mac cannot run Warp-CUDA (I0-A), so every physics
kernel is authored as SOURCE and gated on the A5000 by the lead. Locally, the closed-form oracles
(Terzaghi/Green/manufactured/Darcy/FRAP) are the CONTINUUM ground truth, and `fv_reference` /
`transport_reference` implement the EXACT conservative stencil the Warp kernel ports — so the discrete
conservation identities a closed form cannot express (content==flux, `∫φc+N_polymer` to round-off, adjoint
IBM) are green before the GPU slot opens. The Warp↔NumPy round-off parity is native gate NG-0.

## Figures

Regenerate all: `python ffn_sim/scripts/ac_fluid_vis.py` (pure numpy/matplotlib on the committed oracles;
run with the worktree root on `PYTHONPATH`). Every figure overlays the closed-form oracle/band on the
numeric result, annotates SI units, and does not truncate axes (the visualization-integrity rules).

- **`figs/i1a_terzaghi.png`** — Terzaghi 1-D consolidation. U(T_v) oracle series vs the Terzaghi-1943
  tabulated band (U=0.5 @ T_v≈0.197); excess-pressure profiles dissipating symmetrically (drained z/H=0,2).
- **`figs/i1a_greens.png`** — diffusion heat kernel. Impulse spread (σ=√(2 c_v t)); numeric ⟨r²⟩=∫r²G dV on
  the analytic 2·d·c_v·t line (conservation + the diffusion law reading c_v back).
- **`figs/i1a_manufactured.png`** — cosine eigenmode decays at λ=c_v k² (log axis); the moving-boundary
  conservation identity dΦ/dt = moving-face + Darcy flux + source matches a finite-difference of Φ(t).
- **`figs/i1a_convergence.png`** — the DISCRETE FV reference (the Warp kernel's CPU twin) sits on the
  Terzaghi oracle for n=40/80/160 and its L2 error follows the 2nd-order slope ∝ dx² — the discretisation
  the device kernel ports is verified to converge.
- **`figs/i1b_darcy.png`** — Darcy discharge linear in Δp and in mobility k/µ; `v_f−v_s=q/φ` for a sweep of
  the (GAP) porosity φ (slope 1/φ). The fluid FLOWS.
- **`figs/i1c_transport.png`** — total-actin conservation residual at ~2e-16 (machine precision, far under
  the 1e-11 gate); FRAP cosine bleach recovery on exp(−D_c k² t) (reads D_c); the advection blob centroid
  rides `x₀+v_f·t` (correct) and is clearly separated from `x₀+q·t` (q=φv_f, the wrong solute velocity).

## Native-render viz spec (lead runs on gbook — NG-10)

Interactive 3-D cell-morphology **HTML** on the LIVE membrane-minus-nucleus geometry (peel/slab/cut,
**full-res, no downsampling**, browser-verified via `browser_check.py` — a WebGL grep is meaningless):
the spatial **pore-pressure** field `p_excess` [Pa], the reconstructed **pore-fluid velocity** `|v_f|`
[µm/s] + `q` discharge streamlines, and the membrane hydraulic-flux / nucleus no-flux shell. An anomalous
field → STOP + render + surface to PI (§1.9). See `NATIVE_GATE_SPEC.md` for NG-0…NG-9.

## I0-B GAPs surfaced to PI (do NOT choose)

- **φ** (porosity) — I1b/I1c; the ~0.7 water fraction is an UPPER proxy, not the Darcy pore fraction (GAP).
- **c₀** (initial free G-actin, MCF7-specific) — I1c; total ~100 µM, free ~tens µM (GAP).
- **D_c** (monomer diffusivity) — I1c; draft ~2-6 µm²/s (KB-DRAFT-7-03), ratification + MCF7 applicability OPEN.

The analytic gates are parameter-agnostic where a GAP exists (φ, c₀ swept as oracle parameters), so the CPU
suite is fully green now; only the NATIVE magnitude verdicts (NG-4/7/8) wait on PI.
