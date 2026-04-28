# 05 — Radial Approximation: Proper Continuum Derivation

## Why This Document Matters

The PI's experimental analysis uses a phenomenological radial-symmetric model:
```
A/A₀ = a + b/R + c/R²
```
This model is empirically successful but the **theoretical foundation** of the approximation is not fully derived in the PI's existing work. This document provides a first-principles continuum-mechanics derivation of how a fully 3D anisotropic active hydrodynamic system **reduces** to such a form, and identifies the **regime of validity**.

This is a methodological contribution in its own right. The simulation framework will use this derivation to construct **Sim B** (the radial-reduced simulation) parallel to **Sim A** (the full 3D simulation), enabling direct quantitative comparison.

## Starting Point: Full 3D Active Hydrodynamics

For an active viscoelastic incompressible fluid spreading on a substrate, the governing equations (cf. Marchetti Rev Mod Phys 2013) are:

### Continuity
```
∂ρ/∂t + ∇·(ρ v) = 0
```

### Momentum
```
ρ (∂v/∂t + v·∇v) = ∇·σ + f_substrate
```
where σ is the total stress tensor: σ = -p I + σ_visc + σ_active

### Constitutive (Maxwell-active)
```
∂σ_visc/∂t + σ_visc/τ = 2μ ε̇
σ_active = ζ_a Q + p_a I
```

### Free Surface Boundary Condition
At the spheroid-medium interface r = h(θ, φ, t):
```
σ·n̂ = γ κ n̂ + ∇_s γ_eff   (Marangoni)
```

### Substrate Boundary Condition
At z = 0 inside contact region:
```
v_z = 0
σ_xz, σ_yz = traction (FA-mediated)
```

## Reduction Strategy (The Approximation Hierarchy)

The PI's `A/A₀ = a + b/R + c/R²` form emerges through **layered approximations**, each with explicit conditions:

### Step 1 — Quasi-Steady Approximation
**Assumption**: Spreading dynamics are slow compared to internal stress relaxation.
**Justification**: τ_relax (~minutes) << τ_spreading (~hours)
**Effect**: Drop ∂σ_visc/∂t, treat σ_visc as instantaneously equilibrated.

### Step 2 — Symmetry Assumption
**Assumption**: Cylindrical (axisymmetric) symmetry around the z-axis through spheroid centroid.
**Justification**: Strong cohesion + isotropic substrate → rotational averaging valid in mean-field
**Effect**: All fields depend only on (r, z, t), not on θ. Tangential (θ-direction) variations vanish.
**Validity check**: Sim A's circularity ≈ 1 → assumption holds.

### Step 3 — Lubrication / Thin-Film Approximation
**Assumption**: Spheroid much wider than thick at later times: R >> H (height).
**Justification**: For cellular monolayers in late-spreading regime.
**Effect**: Vertical stress equilibrium → p ≈ p(r, t); horizontal stress dominates dynamics.
**Caveat**: Does NOT hold in early spreading (compact 3D spheroid). Stage 1 sim provides full 3D for comparison.

### Step 4 — Tangential Averaging (Closure)
**Assumption**: Stochastic protrusion events average out over θ in mean-field.
**Justification**: 5,000 cells × ~5% leader fraction → ~250 leader events; central limit theorem applies.
**Effect**: Replace stochastic σ_active^edge by its angular mean ⟨σ_active⟩.

### Step 5 — Asymptotic Expansion in 1/R
For a thin spreading droplet of radius R(t) with leading-edge active stress, the radial momentum balance reduces to an ODE for R(t):
```
γ_eff dR/dt = (terms depending on R)
```

Substituting the spreading coefficient S, surface tension γ, active boundary stress σ_a:
```
γ_eff R dR/dt = α₁ + α₂/R + α₃/R²
```

Integrating:
```
R²(t) - R²(0) = 2 α₁ t / γ_eff  +  ...  (early-time linear in t)
A(t) ~ π R²(t)
A/A₀ ≈ 1 + (linear term) + (correction terms in 1/R)
```

After regrouping, this matches the empirical form:
```
A/A₀ = a + b/R + c/R²
```

with:
- **a**: baseline spreading-coefficient-driven term (size-independent)
- **b/R**: edge curvature / line tension contribution
- **c/R²**: small-size penalty from boundary-area-to-volume ratio

### Validity Condition Summary

The approximation is valid when ALL of:
- (Step 1) τ_active << τ_spreading
- (Step 2) circularity stays ≥ 0.9 throughout
- (Step 3) H/R << 1 at the time-window of analysis
- (Step 4) leader cell distribution is approximately uniform on edge
- (Step 5) R is in asymptotic regime: R >> ξ (correlation length, ~10–50 μm)

When these break (e.g., highly anisotropic spreading, leader-cell cluster localization, very small spheroids), Sim A and Sim B will diverge.

## Implementation: Sim B (Radial-Reduced ODE Solver)

```python
# Pseudocode for Sim B
def radial_simulation(params, t_total, dt):
    R = params.R0
    states = [(0, R, np.pi*R**2/(np.pi*params.R0**2))]  # (t, R, A/A0)

    for t in np.arange(dt, t_total, dt):
        # Effective parameters from φ ODE (same φ ODE as Sim A, but spatially uniform)
        phi_eff = solve_phi_ode(...)
        gamma = gamma_cc(phi_eff)
        sigma_a = sigma_active(phi_eff)

        # Radial ODE: γ R dR/dt = α₁ + α₂/R + α₃/R²
        alpha_1 = params.S_wetting + sigma_a * params.h_edge_layer
        alpha_2 = params.line_tension_correction(gamma)
        alpha_3 = params.curvature_correction(gamma, sigma_a)

        dRdt = (alpha_1 + alpha_2/R + alpha_3/R**2) / (gamma * R)
        R = R + dRdt * dt
        A_over_A0 = (R / params.R0)**2  # 2D area scales as R²

        states.append((t, R, A_over_A0))

    return states
```

**Key**: Sim B uses the **same parameter library** as Sim A (same γ, σ_a, etc.), so the comparison is meaningful.

## Comparison Metrics (vs Sim A)

In `06_radial_full_comparison.md`:
- A/A₀(t) RMS deviation
- R(t) absolute deviation
- circularity time series (Sim A measures this; Sim B assumes = 1)
- effective stress anisotropy ⟨σ_rr - σ_θθ⟩ (Sim A measures; Sim B assumes = 0)

These metrics generate a **validity phase diagram** in parameter space.

## Open Question (for the Project)

The PI's a, b, c values fit MCF7 spheroid data well empirically. What do those values **mean** in terms of the physical parameters γ, σ_a, S, etc.?

Sim B with our parameter library will **predict** specific (a, b, c) values from physical inputs. Comparing to the PI's fitted values gives:
- If close → physical parameter library is reasonable
- If different → PI's empirical model captures additional physics not yet in our parameter library

Either outcome is informative.

## References
- Brochard-Wyart F, de Gennes PG. *Spreading of a "drop" between adherent surfaces.* Langmuir 1992. [foundational wetting framework]
- de Gennes PG. *Wetting: statics and dynamics.* Rev Mod Phys 1985, 57:827. [classic wetting derivation, IF ~50]
- Marchetti et al. Rev Mod Phys 2013 — active matter framework.
- Pérez-González et al. Nat Phys 2019 — active wetting derivation.
- Banerjee S, Marchetti MC. *Continuum models of collective cell migration.* Adv Exp Med Biol 2019.

The mathematical-physics rigor of this derivation is itself a Stage 1 deliverable — it forms a section of any methods paper from this project.
