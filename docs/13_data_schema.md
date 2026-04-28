# 13 — Data Schema (Output Storage)

## Directory Structure (per run)

```
results/
└── {run_id}/                          # e.g., "20260501_pilot_p1_seed42"
    ├── config.yaml                    # exact config used (auto-saved)
    ├── git_commit.txt                 # commit hash
    ├── sim.log                        # log file
    ├── snapshots.h5                   # HDF5 with all material-point states
    ├── grid_fields.h5                 # HDF5 with background grid fields
    ├── metrics.csv                    # time series of summary metrics
    ├── frames/                        # PNG frame sequence (Vispy)
    │   └── frame_{t:05d}.png
    ├── render/
    │   └── cinematic.mp4              # Blender render
    ├── figures/                       # static matplotlib outputs
    └── dashboard.html                 # interactive plotly dashboard
```

## HDF5 Schema (`snapshots.h5`)

```
/
├── metadata/
│   ├── n_material_points         (int, attrs: scalar)
│   ├── domain_size_um            (float[3])
│   ├── total_sim_seconds         (float)
│   ├── frame_interval_seconds    (float)
│   ├── n_frames                  (int)
│   └── created_utc               (string)
└── frames/
    ├── 00000/
    │   ├── time_seconds          (float)
    │   ├── position              (float[N, 3], μm)
    │   ├── velocity              (float[N, 3], μm/s)
    │   ├── deformation_gradient  (float[N, 3, 3])
    │   ├── stress_tensor         (float[N, 3, 3], Pa)
    │   ├── phi                   (float[N], dimensionless)
    │   ├── rho_relative          (float[N])
    │   ├── polarization          (float[N, 3])
    │   ├── is_boundary           (bool[N])
    │   ├── is_leader_cell        (bool[N])
    │   └── lamellipodia_active   (bool[N])
    ├── 00001/
    │   └── ...
    └── ...
```

## HDF5 Schema (`grid_fields.h5`)

Background grid fields, written less frequently than material points (each k-th frame):

```
/
├── metadata/...
└── frames/
    ├── 00000/
    │   ├── time_seconds
    │   ├── density               (float[Nx, Ny, Nz])
    │   ├── velocity              (float[Nx, Ny, Nz, 3])
    │   ├── stress                (float[Nx, Ny, Nz, 3, 3])
    │   ├── vorticity             (float[Nx, Ny, Nz, 3])
    │   ├── nematic_Q             (float[Nx, Ny, Nz, 3, 3])  # if Layer 4 active
    │   └── gamma_eff             (float[Nx, Ny, Nz])         # if Layer 4 active
    └── ...
```

## Metrics CSV Schema (`metrics.csv`)

One row per output frame:

| column | unit | description |
|---|---|---|
| time_seconds | s | simulation time |
| frame_index | - | frame counter |
| area_um2 | μm² | projected 2D area at z = 0+ slice |
| A_over_A0 | - | normalized area ratio |
| effective_radius_um | μm | √(A/π) |
| circularity | - | 4πA/P² |
| height_max_um | μm | max z extent of spheroid |
| centroid_z_um | μm | z-coordinate of mass centroid |
| mean_phi | - | spatial mean of φ |
| max_phi | - | max φ (typically at edge) |
| mean_rho_relative | - | spatial mean of ρ |
| total_mass_g | g | total mass (conservation check) |
| total_kinetic_energy_J | J | sum of (1/2)mv² over points |
| total_strain_energy_J | J | elastic strain energy |
| stress_radial_mean_Pa | Pa | ⟨σ_rr⟩ (for radial-aware analysis) |
| stress_tangential_mean_Pa | Pa | ⟨σ_θθ⟩ |
| n_lamellipodia_active | - | count of active lamellipodia events |
| n_leader_cells | - | count of leader cells |
| coffee_ring_index | - | edge-to-core density ratio |
| marangoni_force_mean | N/m² | mean ∇_s γ_eff magnitude |
| vorticity_mean | 1/s | mean |ω| |
| nematic_order_S | - | scalar nematic order parameter |

If a metric is not applicable (e.g., Marangoni metric in a Stage 1a run with Layer 4 disabled), record NaN.

## Run Naming Convention

```
{date}_{stage}_{point}_{seed}
```
- date: YYYYMMDD (start date)
- stage: pilot, prod, sweep
- point: condition identifier (p1, p2, p3, p4, p5 for sweep; baseline for pilot)
- seed: RNG seed (for reproducibility)

Examples:
- `20260501_pilot_baseline_seed42`
- `20260515_prod_p1_seed1`
- `20260515_prod_p1_seed2` (replicate)

## Cross-Run Aggregation

A separate database (`results/index.csv`) lists every run with metadata for sweep aggregation:

```
run_id, stage, point_id, seed, status, start_utc, end_utc, wall_clock_s, hardware, git_commit, config_path
```

This enables `pandas` queries like "all production runs for point p1" → list, load each metric.csv, compute statistics.

## Compression
HDF5 with gzip level 4 (default in h5py). Material point arrays are compressible to ~30% original size.

## Reading the Data

```python
import h5py
import pandas as pd

# Material points at frame 50
with h5py.File('results/run_id/snapshots.h5', 'r') as f:
    pos = f['frames/00050/position'][:]  # [N, 3]
    phi = f['frames/00050/phi'][:]       # [N]

# Time series
metrics = pd.read_csv('results/run_id/metrics.csv')
metrics.plot(x='time_seconds', y='A_over_A0')
```

## Hooks for Future Extensions
Layer 6 (chemistry, future):
- Add `chemistry/oxygen_field`, `chemistry/glucose_field`, `chemistry/viability_v` to schema

Cell heterogeneity (future):
- Add `material_properties/cortex_stiffness_Pa` per material point (currently uniform)
