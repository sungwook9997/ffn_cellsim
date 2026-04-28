# 07 — Internal Flow Dynamics (Drying-Thin-Film Inspired)

## Conceptual Origin

This layer is the project's **unique angle**. The PI's prior research includes experimental work on **thin-film drying** with self-alignment phenomena driven by Marangoni and capillary effects. The hypothesis explored in this simulation: **analogous mechanisms exist in cellular spreading**, and capturing them computationally generates new mechanobiological insight.

Three correspondences are exploited:

1. **Marangoni**: γ-gradient → tangential flow at interface → in cells: φ-gradient creates effective γ-gradient → cellular Marangoni flow
2. **Coffee-ring (inverse)**: drying creates edge accumulation → in cells: active edge protrusion creates inverse edge density profile (cells thin out at center, accumulate at edge)
3. **Drying-induced concentration**: water evaporation increases solute density → in cells: spreading-driven water efflux increases cytoplasm density (handled in Layer 5, complementary)

These correspondences are **academically grounded** (not invented):
- Pajic-Lijakovic & Milivojevic, *Marangoni effect and cell spreading*, Eur Biophys J 2022.
- Fütterer et al., *Gradients in solid surface tension drive Marangoni-like motions in cell aggregates*, Phys Rev Fluids 2022 (laser-ablation experiment).

But the **systematic incorporation** of all three drying-thin-film correspondences into a single cellular simulation framework is novel.

## §1 — Cellular Marangoni Flow

### Effective Surface Tension Field
At every material point, compute:
```
γ_eff_i = γ_cc(φ_i) · w_local_density_i
```
where w_local_density accounts for cell density (sparse regions have weaker effective surface tension).

### Surface Gradient
On the spheroid free surface, compute tangential gradient:
```
∇_s γ_eff = (I - n̂ ⊗ n̂) · ∇γ_eff
```
n̂ is the local outward normal at the free surface.

### Marangoni Force
Tangential body force per unit interfacial area:
```
f_Marangoni = ∇_s γ_eff
```
This drives flow from low-γ to high-γ regions along the interface.

In our simulation:
- Edge regions enriched in laminin/integrin (high φ) → low γ_cc → low γ_eff
- Spheroid core (low φ) → high γ_cc → high γ_eff
- Result: Marangoni flow pulls material from edge toward center — but combined with active spreading at edge, complex internal recirculation emerges.

### Implementation Notes
- Compute γ_eff on background grid via P2G transfer
- Take gradient on grid
- Project to surface (using level-set or implicit surface representation)
- Apply tangential force to surface material points

## §2 — Active Nematic Order

### Order Parameter Q
For each material point with polarization vector **p** (3D), construct local nematic tensor:
```
Q_i = p_i ⊗ p_i - (1/3) I
```
Coarse-grain to grid:
```
Q_grid(x) = ⟨Q_i⟩_local
```

### Active Nematic Stress
```
σ_active_nematic = -ζ_Q · Q_grid
```
- ζ_Q > 0: extensile (rod-like cells push outward along their long axis)
- ζ_Q < 0: contractile

### Topological Defects
For 2D projection at substrate (z = 0+) of the spreading monolayer:
- ±1/2 defects can form
- Detection: contour integral of nematic director angle around closed loops
- These defects are biologically meaningful (cell extrusion sites, leader cell positions)

### Reference
- Saw TB et al. *Topological defects in epithelia govern cell death and extrusion.* Nature 2017, 544:212. [IF 65]
- Duclos G et al. *Topological defects in confined populations.* Nat Phys 2017, 13:58.

## §3 — Internal Flow Diagnostics

### Vorticity Field
On background grid:
```
ω = ∇ × v
```

For 3D simulation, ω is a vector field. For axisymmetric reduction, only ω_θ matters (toroidal recirculation).

### Streamlines
Numerical particle tracing through the velocity field. Visualization shows internal recirculation patterns predicted by Phys Rev Fluids 2022.

### Recirculation Identification
- Detect closed streamlines via Q-criterion or λ₂-criterion (standard CFD techniques)
- Quantify recirculation strength as ⟨|ω|²⟩ over recirculating volume

## §4 — Coffee-Ring Analog

### Radial Density Profile
At any time t, compute the radial cell density profile:
```
n_cells(r) = ∫∫ δ(r - r_i) dA  for material points within thickness h
```

### Coffee-Ring Index
```
CRI(t) = max[n_cells(r) over edge region] / mean[n_cells(r) over interior region]
```
- CRI > 1: edge accumulation (positive coffee-ring)
- CRI ≈ 1: uniform
- CRI < 1: edge depletion (inverse coffee-ring; possibly seen in necrosis-driven systems, Stage 2)

### Drying-Inspired Hypothesis
For Stage 1 active spreading on Col1: leader cells + lamellipodia at edges should produce CRI > 1 → analog of standard coffee-ring. Simulation will quantify this.

## §5 — New Dimensionless Number: Cellular Marangoni Number

Define:
```
Ma_cell = (Δγ_cc · L) / (η_eff · v_spreading)
```
where:
- Δγ_cc = γ_cc(φ=0) - γ_cc(φ=1) ≈ 1.7 mJ/m²
- L = spheroid radius scale ~ 200 μm
- η_eff = effective tissue viscosity ~ 10⁵ Pa·s
- v_spreading ~ 0.1 μm/s

Order-of-magnitude: Ma_cell ~ (10⁻³ × 10⁻⁴) / (10⁵ × 10⁻⁷) = 10⁻⁵

Compare to thin-film: Ma_film ~ 1–100. So cellular Marangoni is **much weaker** than thin-film Marangoni in absolute terms — but acts over **much longer timescales** (hours vs seconds).

This dimensional analysis is a Stage 1 deliverable. It defines the regime where Marangoni effects are observable in cellular spreading.

## §6 — Visualizations Specific to This Layer

In `09_visualization.md`, this layer adds:
- γ_eff field heatmap on free surface
- Tangential Marangoni force vector field
- Q tensor ellipsoids on background grid
- Streamline animation through interior
- Vorticity isosurfaces
- Defect tracking in 2D substrate slice
- Radial density profile time series (coffee-ring evolution)

## Output Metrics (for Analysis Dashboard)
Stored per frame in metrics CSV:
- mean(|∇_s γ_eff|) — Marangoni driving force magnitude
- mean(|ω|) — vorticity magnitude
- nematic order parameter S = ⟨3 cos²θ - 1⟩/2
- defect count
- CRI(t) — coffee-ring index

## Limitations
- 3D nematic order is more complex than scalar S; we use coarse approximation
- Defect detection in 3D is non-trivial (defect lines, not points)
- Marangoni force can interact with substrate friction in non-trivial ways — sanity-check with momentum conservation
