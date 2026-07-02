# FF active movement — piece 2/5: cell polarization (actomyosin contractile-flow instability) 2026-07-02

**Date:** 2026-07-02  **Engine:** FF (Warp/numpy, µm·pN·s)  **Branch:** dcm/main  **Scope: 2 of 5.**
Piece-1 = protrusion engine (polymerization). Piece-2 = the front-rear axis a cell must break BEFORE it can
migrate directionally. Mechanism chosen (FF_POLARIZATION_LITERATURE_2026-07-02): **A — actomyosin
contractile-flow instability** (Bois-Jülicher-Grill 2011), over Rho wave-pinning, per the mechanistic-not-lumped
rule.

## What this is

`ff/polarization_activegel.py` — the active-gel reduction of the actomyosin cortex on a 1-D periodic ring
(arclength x∈[0,L=2πR)), the mechanistic hydrodynamic theory of the explicit myosin the engine carries:

```
myosin :  ∂_t c = -∂_x(c v) + D ∂²_x c - k_off (c - c0)          (advection + diffusion + turnover)
force  :  γ v - η ∂²_x v = ∂_x σ_a ,  σ_a = ζ f(c) , f(c)=c/(1+c/c*) ;  ℓ=√(η/γ)
```
Positive feedback: myosin → active stress → flow → flow advects myosin → density ↑. Above ζ_c a single
high-myosin cap (the rear) + steady cortical flow emerge = polarity. Solved spectrally (periodic FFT).

## Validation (analytic ground truth FIRST, per oracle-is-crosscheck)

`tests/ff/test_polarization_activegel.py` (4/4 PASS):
- **Dispersion:** the spectral solver's per-mode numerical growth rate == the analytic λ(k) =
  `[c0 ζ f'(c0)/γ]·k²/(1+ℓ²k²) − D k² − k_off` to **rel_err <0.1%** across modes 1–10 (unstable low-k,
  stable high-k band). *(This verifies the solver integrates the model PDEs correctly — a manufactured-solution
  check against the model's own linear theory, not an experimental validation.)*
- **Threshold:** homogeneous state stable (λ<0) at ζ=0.9·ζ_c, unstable (λ>0) at ζ=1.1·ζ_c — onset brackets ζ_c.
- **Nonlinear:** from random noise (ζ=4·ζ_c) a **single** high-myosin cap forms (dominant Fourier mode = 1,
  density contrast 12.3×) with cortical flow converging into it — spontaneous symmetry-breaking.

## Constants (µm·pN·s; DOIs in FF_POLARIZATION_LITERATURE; flagged for PI KB-registration)

- D=1 µm²/s, k_off=0.1 s⁻¹ — Bois 2011 (author order-of-magnitude estimates).
- ℓ=14 µm — Mayer 2010 MEASURED hydrodynamic length (⚠️ C. elegans; **no MCF7 datum → PI-gated**).
- ζ (contractility) — the SWEPT control variable (like the γ-floor myosin sweep; NOT tuned to an outcome).
- c0=c*=1 normalized; γ=1 normalization.

## Honest scope + caveats

- **1-D active-gel REDUCTION, not the full 3D fine-grained cortex.** This is the validated theoretical core +
  the analytic oracle for polarization. The production increment is the *fine-grained emergence*: explicit
  myosin minifilaments advecting + turning over on the 3D cortex mesh (network_warp), where the cap forms from
  the particle dynamics. Mirrors piece-1 (kernel validated vs analytic v(f) → then coupled into the cell).
- **Absolute flow speed is NOT claimed.** The nonlinear cap flow (~4.6 µm/s here) is ζ/γ-scaled; γ (cortical
  friction) is a normalization (Mayer gives only η/γ=ℓ²). Only the **length scale ℓ** and the **instability
  onset** are physical. An absolute cortical friction/viscosity for MCF7 → surface to PI.
- **Dispersion match is solver-verification** (against the model's own linear theory), not experiment. The
  MODEL (Bois active-gel) is the lit-anchored mechanism; its predictions (threshold, single-cap coarsening) are
  genuine.

## Files
- `ff/polarization_activegel.py` — resolve/dispersion/zeta_critical/step/measure_growth_rate.
- `tests/ff/test_polarization_activegel.py` (4/4).
- `outputs/ff/figs/polarization_activegel.png` — dispersion + instability window + kymograph + polarized state.
- `outputs/ff/figs/polarization_activegel.html` — interactive: myosin density on the cortex ring (radial bulge =
  cap) animated as the cell polarizes (▶ play).

Related: [[project-ff-active-movement-pieces]], FF_POLARIZATION_LITERATURE_2026-07-02, FF_POLYMERIZATION_PIECE1.
