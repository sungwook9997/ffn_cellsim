# Phase 1 — Unit 1.2 Progress Report

**Date**: 2026-05-18 (Day 1, third pass — Unit 1.2)
**Scope**: ECM Unit 1.2 — overdamped Langevin time integration, CFL Sanity Gate, equipartition validation, **preliminary** strain-stiffening signal check.
**Repo HEAD**: `4ffcbbe` (`acs_kb/` untracked)
**KU sources**: KU-1.1, KU-1.2, KU-1.4, KU-1.7, KU-1.9, KU-1.22, KU-1.24, KU-1.26, KU-1.27, KU-1.28, KU-1.29, KU-1.30
**Mode**: production (`demo_mode = false`)

---

## What was built (delta vs Unit 1.1)

```
acs_kb/
├── common/
│   ├── derived_params.py   ← + Stokes γ_b, τ_xl/τ_stretch/τ_bend, dt_cfl
│   └── sanity_gate.py      ← + gate_unit1_2_dynamics, gate_equipartition
├── configs/
│   └── phase1_unit1.yaml   ← + water_viscosity, dynamics block, derived.γ_b/τ/dt
├── ecm/
│   ├── integrator.py       ← NEW: Integrator protocol, EulerMaruyama, BAOAB stub, run()
│   ├── diagnostics.py      ← NEW: link_extension_energies
│   ├── shear_protocol.py   ← NEW: apply_affine_shear, compute_virial_shear_stress_2d,
│   │                              preliminary_affine_G_curve
│   └── (fiber_*, cross_links unchanged from Unit 1.1)
├── tests/
│   ├── test_dynamics.py    ← NEW: 7 tests (CFL, free diffusion, equipartition,
│   │                              shear sanity, virial-zero-at-rest, KU-1.30 OOM)
│   └── test_ecm.py         (unchanged, still 7)
├── notebooks/
│   └── 02_unit1_2_dynamics.py    ← NEW: free diffusion + equipartition + G(γ)
└── outputs/phase1/
    ├── REPORT_unit1_2.md         ← this file
    ├── 04_free_diffusion.png
    ├── 05_equipartition.png
    ├── 06_preliminary_strain_stiffening.png
    └── 02_summary.json
```

## KU validation numbers

| Quantity | Measured | Theory / reference | KU | Status |
|---|---|---|---|---|
| Stokes drag γ_b | **6.515e-10 N·s/m** | 6 π η_water bead_radius = 6.5e-10 | KU-1.26 | derived |
| τ_xl | 6.5e-7 s | γ_b / k_xl | KU-1.26 | derived |
| τ_stretch | **1.9e-7 s** (smallest) | γ_b ℓ₀ / μ | KU-1.26 | sets CFL |
| τ_bend | 1.5e-1 s | γ_b ℓ₀³ / κ | KU-1.26 | non-binding |
| dt (default) | **1.88e-8 s** | 0.99 · 0.1 · τ_min | KU-1.26 CFL | safe by 1 % margin |
| Free-bead MSD vs 4 D t | rel err **1.6 %** at t = 200 ns | thermal floor ~ 1/√n_walkers | KU-1.26 | PASS |
| Equipartition ⟨½k(|d|−r₀)²⟩ / ½ k_BT | **1.031** (3 % above) | 1.0 | KU-1.26 | PASS (tol 12 %) |
| G_affine(γ→0) using h=ξ=2 μm | **32 Pa** | 36 Pa (KU-1.30 #1) | KU-1.30 | OOM PASS |

## Test suite

```
acs_kb/tests/ — 14 passed in 4.93 s   (7 Unit 1.1 + 7 Unit 1.2)

test_dynamics.py:
  ✓ test_cfl_gate_passes_on_default
  ✓ test_cfl_gate_rejects_excess_dt              # 100× over CFL must raise
  ✓ test_free_bead_diffusion                     # ⟨|Δr|²⟩ = 4 D t to < 5 %
  ✓ test_equipartition_on_small_network          # < 12 % from ½ k_BT
  ✓ test_affine_shear_is_pure_x_shift
  ✓ test_virial_stress_zero_on_rest_network      # at rest, σ_xy = 0 exact
  ✓ test_preliminary_G_in_KU130_order_of_magnitude
```

## Sanity Gate outputs

```
unit1_2_overdamped_langevin
  PASS cfl_dt_below_alpha_tau_min        dt=1.88e-08, α·τ_min=1.89e-08 (α=0.1)
  PASS integrator_supported_in_phase_1   euler_maruyama

unit1_2_equipartition
  PASS link_energy_in_canonical_ensemble ⟨½k|Δr|²⟩=2.21e-21 J vs ½k_BT=2.14e-21 J,
                                          rel err 0.031, tol 0.12, n_samples=22200
```

## Performance

For default `n_fibers = 3142` network only used for assertions; the Unit 1.2 walkthrough uses smaller networks tuned to keep wall-clock manageable:

| Step | Wall time |
|---|---|
| Free diffusion (4000 walkers × 400 steps) | 87 ms |
| Equipartition (400-fiber net, 12 k equilibrate + 200 × 50 sample steps) | 6.0 s |
| Affine G(γ) sweep (13 γ values, no relax) | 0.22 s |

The equipartition test dominates wall time; further reduction would need a tinier network or Numba'd force kernel — not warranted at Phase 1.

## Why the strain-stiffening reporting differs from the brief

The initial brief asked for a Lees-Edwards quasi-static shear + Storm-MacKintosh `G(γ) ∼ G_0/(1−γ/γ_c)²` fit. Per PI direction:

> "Before quantitative shear/strain-stiffening publication-grade validation, revisit exact intersection-node insertion or finer bead discretization."

Phase 1 therefore ships only a preliminary affine measurement:

1. **No Lees-Edwards.** Without it, a Langevin run after affine shear relaxes back to the unsheared equilibrium (PBC do not pin the deformation). σ_xy → 0 was the symptom of this in the first try; the second cut measures the *instantaneous affine virial* with no relaxation, which is the affine upper bound.
2. **No quantitative G_0 / γ_c fit.** With our nearest-bead XL anchoring and only 5 beads/fiber, the bond-virial undercounts non-affine reorganisation and overcounts compressive bonds at large γ. The 2D virial is reported in N/m (its natural units); only when divided by an explicit "ECM slab thickness" `h = ξ = 2 μm` does it convert to 3D Pa for the optical comparison.
3. **What we observe.** With h = ξ = 2 μm convention, |G_affine(γ→0)| ≈ 32 Pa — in the KU-1.30 #1 acceptance band [15, 200] Pa around the experimental 36 Pa. Beyond γ ≈ 0.04 the affine assumption breaks (some bonds compress, the virial flips sign). This is **order-of-magnitude agreement, not validation**.
4. **What it tells us.** The microscopic ingredients (μ, κ, ρ_L, k_xl, ℓ₀) plus the 2D→3D thickness convention give the right ballpark for collagen at 1.5 mg/mL. A publication-grade comparison needs (i) exact intersection-node insertion (so XL forces sit at the true crossing point, not the nearest bead), (ii) finer bead discretization (so bond reorientation is resolved), (iii) Lees-Edwards PBC (so non-affine relaxation can actually develop).

## Out of scope (Phase 1)

- Lees-Edwards periodic BC (deferred Phase 2+, KU-1.29).
- BAOAB integrator (stub only, Phase 2+).
- Exact intersection-node insertion (replaces nearest-bead anchoring at quantitative stage).
- Bell-Evans XL kinetics (KU-1.15, Phase 2+).
- 3D promotion (collagen biology only fully matched in 3D).

## Suggested next step

Two reasonable forks; PI decides:

**(a) Tighten Unit 1.2 to publication grade first.** Implement Lees-Edwards minimum-image, add exact intersection-node insertion or 10×–20× bead refinement, sweep G(γ) with seed averaging, fit Storm-MacKintosh, publish G_0 ≈ 36 Pa / γ_c ≈ 0.16 (KU-1.30 #1, #2). Estimated effort: 1–2 days.

**(b) Advance to Unit 2.1 (ECM-Cell Bridge) first.** Bring the focal-adhesion / motor-clutch stack online with the current preliminary-grade ECM; the FA dynamics validation does not need quantitative G_0. Strain-stiffening tightening becomes a parallel polish track. Estimated effort: 3–5 days. (The freshly added `resolve_bridge` / `_pereverzev_peak` helpers in `derived_params.py` suggest groundwork for this path is already starting.)

**Awaiting PI direction.**
