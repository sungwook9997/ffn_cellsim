# 01 — Model Architecture (5+1 Layers)

The simulation is structured as **independent toggleable layers**, each with clear physical meaning and academic provenance. Layers stack additively; each adds new physics. Staged activation (see `10_dev_roadmap.md`) starts with Layers 1+2 only.

```
┌─────────────────────────────────────────────────────────────┐
│  Layer 6 — Chemistry  (Stage 2 hook, NOT active in Stage 1) │
│  • O₂/nutrient diffusion fields                             │
│  • Cell viability state v                                   │
│  • Necrosis dynamics                                        │
└─────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│  Layer 5 — Mechano-Osmotic Coupling (Tier 2, phenomenol.)   │
│  • Material-point density ρ(t)                              │
│  • Spreading-rate-driven water efflux                       │
│  • Density-modulated viscosity & cortical stiffness         │
└─────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│  Layer 4 — Internal Flow Dynamics (drying-thin-film inspir.)│
│  • Cellular Marangoni (γ-gradient tangential flow)          │
│  • Nematic order tensor Q                                   │
│  • Vorticity, streamlines, topological defects              │
│  • Coffee-ring radial-density profile metric                │
└─────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│  Layer 3 — Adhesion Network Dynamics                        │
│  • Internal state φ ∈ [0,1]: E-cad ↔ Int-β1 ratio           │
│  • Substrate-driven ODE: dφ/dt = k₊·S(t)·(1-φ) - k₋·φ       │
│  • φ modulates Layer 1+2 parameters (γ_cc, σ_a, ρ_int, ...) │
└─────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│  Layer 2 — Boundary Cell Biology                            │
│  • Edge-layer identification (free-surface tracking)        │
│  • Lamellipodia: stochastic protrusion events               │
│  • Filopodia: directional ECM probing                       │
│  • Discrete focal adhesion patches (catch/slip-bond)        │
│  • Leader cell heterogeneity (stochastic polarization)      │
└─────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│  Layer 1 — Bulk Active Hydrodynamics (foundational)         │
│  • Hybrid Eulerian-Lagrangian MPM (MLS-MPM)                 │
│  • Cell-equivalent material points (5,000)                  │
│  • Background grid for continuum fields                     │
│  • Active viscoelastic stress tensor (Maxwell + active)     │
│  • Cohesion (surface tension γ), traction, gravity, buoyancy│
│  • Substrate contact (Hertz + adhesion)                     │
└─────────────────────────────────────────────────────────────┘
```

## Key Architectural Decision: Hybrid Eulerian-Lagrangian

Material points are treated as **cell-equivalent carriers** (each carries mass, polarization, internal state φ, density ρ — i.e., one point ≈ one cell of "stuff"), while continuum fields (stress, velocity, vorticity) live on a background grid via standard MLS-MPM particle-to-grid (P2G) and grid-to-particle (G2P) transfers.

This choice gives:
- **Cell-level interpretability**: trace individual material points, plot their trajectories, attach state variables
- **Continuum mechanics rigor**: bulk stress, vorticity, divergence computed properly on grid
- **Both visualizations possible**: cell-tracking video AND continuum field heatmap

References for this hybrid approach:
- Hu et al. SIGGRAPH 2018 (MLS-MPM): the algorithmic foundation
- Cellular Potts hybrid extensions (Glazier group): biological precedent
- Active hydrodynamics with discrete carriers (Marchetti review): theoretical grounding

## Independence vs Coupling Between Layers

| Layer | Inputs from below | Outputs to above |
|---|---|---|
| 1 | (foundational) | Stress σ, velocity v, density n_cells |
| 2 | Free-surface ID, local stress | Active boundary force, FA forces |
| 3 | Substrate contact area, integrin engagement signal | φ field, modulated parameters |
| 4 | φ field, velocity v, polarization | γ-gradient force, Q tensor diagnostics |
| 5 | Local strain rate ε̇ | ρ field, modulated viscosity & stiffness |
| 6 | (Stage 2 only) ρ, position | Viability state v, source/sink terms |

Each layer has a clean interface; layers can be disabled via config without breaking the model below.

## Configuration Schema (config YAML excerpt)
```yaml
model:
  layers:
    bulk_hydrodynamics: true       # Layer 1: always on
    boundary_biology: true         # Layer 2
    adhesion_dynamics: false       # Layer 3 (Stage 1b activates)
    internal_flow: false           # Layer 4 (Stage 1d activates)
    mechano_osmotic: false         # Layer 5 (Stage 1c activates)
    chemistry: false               # Layer 6 (Stage 2 only)

  parameters:
    n_material_points: 5000
    domain_size_um: 1000
    grid_resolution: 128
    # ... see other docs for layer-specific parameters
```

## Detail References
- Layer 1 force terms: `02_force_models.md` §1
- Layer 2 boundary biology: `02_force_models.md` §2
- Layer 3 φ ODE: `03_adhesion_dynamics.md`
- Layer 4 internal flow: `07_internal_flow_dynamics.md`
- Layer 5 mechano-osmotic: `08_mechano_osmotic.md`
- Layer 6 (deferred): `10_dev_roadmap.md` Stage 2 section
