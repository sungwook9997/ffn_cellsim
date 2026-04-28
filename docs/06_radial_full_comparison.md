# 06 — Sim A vs Sim B Comparison Protocol

## Purpose
Independent comparison between **Sim A (full 3D anisotropic)** and **Sim B (radial-reduced ODE)** to quantify the validity of the radial approximation in different parameter regimes.

See `05_radial_approximation.md` for the derivation of Sim B from first principles.

## Protocol

### Setup
1. Choose a parameter set P (one of the 5 sweep points from `03_adhesion_dynamics.md`)
2. Generate identical initial conditions:
   - Same initial radius R₀
   - Same initial φ field (uniform)
   - Same physical parameter library
3. Run **Sim A** with full 3D MPM (~5 hours wall-clock)
4. Run **Sim B** with 1D radial ODE solver (~seconds)

Note: Sim B is computationally trivial; the cost is essentially zero. The expensive comparison is running Sim A across the parameter space.

### Metrics

#### Bulk-Level (1D summaries)
- A/A₀(t) RMS deviation: `√⟨(A_A - A_B)² ⟩_t / ⟨A_A⟩`
- R(t) RMS deviation: similarly
- Final A/A₀ deviation: percent error at t = 80 hr

#### Anisotropy Diagnostics (Sim A only)
- Circularity time series: c(t) = 4π·A(t) / P(t)²
- Aspect ratio of fitted ellipse
- Stress anisotropy: ⟨σ_rr - σ_θθ⟩ / ⟨σ⟩

#### Spatial Heterogeneity (Sim A only, missing in Sim B)
- Coffee-ring index CRI(t)
- Marangoni driving force ⟨|∇_s γ_eff|⟩
- Vorticity magnitude ⟨|ω|⟩
- Density gradient (Layer 5 active): max(ρ) - min(ρ)

### Validity Region Definition

The radial approximation is valid where:
```
A/A₀ deviation < 5% AND circularity > 0.9 AND stress anisotropy < 0.1
```

Mark this region in parameter space.

### Phase Diagram Output
For each sweep point, classify as:
- **Radial-valid**: all metrics within tolerance
- **Marginal**: some metrics borderline
- **Radial-invalid**: significant deviations

Plot in 2D parameter slices (e.g., φ_steady vs σ_a) and identify boundary curves.

## Deliverable
A figure (and quantitative table) showing where in the cellular-mechanical parameter space the PI's `A/A₀ = a + b/R + c/R²` form is valid. This is a methodological contribution.
