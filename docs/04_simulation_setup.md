# 04 — Simulation Setup (Initial & Boundary Conditions, Time Stepping)

## Domain & Grid

### Spatial Domain
- 3D box, axis-aligned
- Pilot: 600 × 600 × 600 μm
- Production: 1000 × 1000 × 600 μm (anisotropic; spreading is mostly horizontal)
- Substrate plane at z = 0; medium fills z > 0

### Background Grid (MLS-MPM)
- Pilot: 64³ cells (~9 μm grid spacing)
- Production: 128 × 128 × 64 (~8 μm grid spacing)
- Boundary conditions on grid:
  - z = 0: rigid wall + adhesion energy (substrate)
  - x, y, z_max: free outflow (Neumann zero-flux)
  - Top of medium (z = z_max): traction-free

## Material Points

### Initialization
- N points = 1000 (pilot) or 5000 (production)
- Distributed inside an initial sphere of radius R₀
  - Pilot R₀ ~ 100 μm
  - Production R₀ chosen from {100, 150, 200, 250, 300} μm to span PI's experimental range
- Random uniform packing inside sphere (rejection sampling)
- Center of sphere at (Lx/2, Ly/2, R₀ + small gap) — just above substrate
- Initial velocity = 0
- Initial φ = 0.05 (E-cad dominant default)
- Initial ρ = ρ₀ (equilibrium density, scaled to 1.0)
- Initial Q (nematic) = 0 (no preferred orientation)

### Equilibration Phase
**Before** the "spreading clock" (t = 0) starts:
- Run sim for τ_eq ~ 5–10 sim minutes with substrate adhesion OFF
- Spheroid relaxes to round shape under cohesion + cortical tension only
- This eliminates artifacts from initialization
- After equilibration: enable substrate adhesion → t = 0 begins

## Substrate Model (Stage 1)
- Plane z = 0
- Effective adhesion energy density γ_substrate(x,y) = constant for Stage 1 (uniform Col1 coating)
- Stiffness: rigid (deformation negligible compared to cells)
- Friction: discrete focal adhesion patches (Layer 2)

## Gravity & Buoyancy
- g = 9.81 m/s² in -z direction
- Net force on each material point: F_gz = (ρ_cell - ρ_medium) · V_point · g
- Net density difference is small (~5%) → weak driving force, mostly relevant for initial settling

## Time Stepping

### Adaptive Δt (CFL-based)
For MLS-MPM with explicit time integration:
```
Δt = CFL · min(grid_spacing / v_max, grid_spacing² / (2 · ν))
```
where ν is the effective kinematic viscosity. CFL ~ 0.4 (conservative).

Typical Δt:
- Quiescent phase: Δt ~ 10⁻³ s (1 ms)
- Active spreading: Δt ~ 10⁻⁴ s (0.1 ms)

### Frame Output Interval
- Pilot: every 60 sim min (= 1 frame/hour)
- Production: every 15 sim min (matches PI's experimental imaging cadence)
- Per frame: HDF5 snapshot with all material point states + grid stress fields

### Simulation Duration
- Pilot: 4–8 sim hours
- Production target: 80 sim hours (matches PI's experimental duration)
- Production extension: up to 100+ sim hours when validating long-term behavior

## Stochastic Events (Layer 2)
Random number generation:
- Use **independent RNG seeds per material point** (deterministic given seed)
- Seed = `f(global_seed, material_point_id)` for reproducibility
- Document seed in run config

Events handled:
- Lamellipodia formation: Poisson process, rate λ_lam(stress, φ, ρ_FA)
- Filopodia probing: similar
- FA formation/dissolution: rate-based with force feedback
- Leader cell stochastic transition: bistable Markov chain

## Numerical Schemes

### Spatial Discretization
- MLS-MPM (Hu et al. SIGGRAPH 2018) — affine velocity reconstruction
- B-spline (cubic) interpolation between grid and particles

### Time Integration (Stage 1)
- **Pilot**: Symplectic Euler (explicit, simple, fast prototyping)
- **Production**: Semi-implicit / IMEX scheme for stiff terms (cortical pressure)
- **Future**: Implicit MPM (Stomakhin et al.) for ~10× speedup if needed

### Conservation Properties (Must Verify)
- Mass conservation: exact (MLS-MPM is conservative by construction)
- Linear momentum: conserved in absence of external force (gravity)
- Angular momentum: conserved (verify with no-substrate test)
- Energy: monotonically dissipative under viscoelastic relaxation

## Convergence & Stability Checks
For each new physical mechanism added:
1. Run with mechanism in isolation
2. Verify conservation laws hold
3. Check that Δt → Δt/2 doesn't change result by more than 5%
4. Verify resolution convergence (grid → 2× grid changes result < 5%)

## Configuration Schema (configs/*.yaml)
```yaml
domain:
  size_um: [600, 600, 600]
  grid_resolution: [64, 64, 64]
  substrate_z: 0

initialization:
  n_material_points: 1000
  initial_radius_um: 100
  initial_z_offset_um: 105    # just above substrate
  initial_phi: 0.05
  initial_rho_relative: 1.0
  equilibration_sim_min: 5

time_stepping:
  cfl: 0.4
  adaptive: true
  max_dt_s: 1.0e-3
  total_sim_hours: 8           # pilot
  frame_interval_sim_min: 60   # pilot

random_seeds:
  global_seed: 42

output:
  hdf5_path: "results/run_{run_id}/snapshots.h5"
  metrics_csv: "results/run_{run_id}/metrics.csv"
  config_dump: "results/run_{run_id}/config.yaml"
  git_commit: auto             # auto-fill from git rev-parse HEAD
```
