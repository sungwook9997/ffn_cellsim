# 09 — Visualization Pipeline (3-Tier Auto-Generation)

## Strategy
Three independent visualization tiers run **automatically** at the end of each simulation:

| Tier | Tool | Purpose | Cost |
|---|---|---|---|
| **1. Real-time / dev viewer** | Vispy (offline mode for SSH) | Quick monitoring, debugging | Low |
| **2. Production render** | Blender Cycles (headless) | Publication-quality MP4 | Medium-High |
| **3. Analysis dashboard** | matplotlib + plotly | Quantitative metrics, time series | Low |

All run **headless** (no GUI required) — compatible with SSH workflow.

## Tier 1 — Vispy Frame Sequence

**Output**: `results/{run_id}/frames/frame_{t:05d}.png` sequence

Per frame, render:
- Material points (3D scatter, color = φ or stress or ρ depending on config)
- Spheroid free surface (reconstructed via marching cubes from density grid)
- Substrate plane
- Lamellipodia/filopodia event overlays (small icons)
- Time stamp + key metrics text overlay

Compose to MP4 with ffmpeg automatically post-run.

**Why Vispy**: Lightweight, fast, OpenGL-based, headless via osmesa or virtual framebuffer (`xvfb-run`).

**Why offline mode (no live GUI)**: SSH X-forwarding is too slow for 60+ fps interactive 3D. Frame export → MP4 is faster end-to-end and works headless.

## Tier 2 — Blender Cycles Production Render

**Output**: `results/{run_id}/render/cinematic.mp4`

Workflow:
1. Python script in repo: `viz/blender_render.py`
2. Calls Blender headless: `blender --background --python viz/blender_render.py -- --run_id {id}`
3. Inside Blender: imports HDF5 frames, builds scene per frame, renders with Cycles

Scene composition:
- Spheroid as deformable mesh (subdivided icosphere fitted to free surface)
- Subsurface scattering material (cells look "alive", semi-translucent)
- Stress heatmap as material color (vertex colors driven by stress field)
- Substrate plane with Col1-coating texture
- Lamellipodia visualized as Niagara-style particle bursts
- Camera movement: orbit + zoom-in on key events
- Volumetric lighting for depth

GPU-accelerated rendering on RTX A5000 via Cycles CUDA backend.

**Why Blender**: Pixar-level rendering; subsurface scattering for cells; can be fully scripted; free; great GPU usage of A5000.

**Cost**: ~30 min – 2 hr per simulation render (post-sim, runs after sim completes).

## Tier 3 — Analysis Dashboard

**Output**: `results/{run_id}/dashboard.html` (Plotly interactive) + `results/{run_id}/figures/*.png` (static matplotlib)

Panels per run:
1. **A/A₀ vs t**: simulation curve, optionally overlay PI's experimental curves (read from `data/experimental/`)
2. **R(t) and circularity**
3. **φ field evolution**: heatmap or radial profile vs time
4. **ρ density field** (if Layer 5 active): radial profile vs time
5. **Stress tensor diagnostics**: ⟨σ_rr⟩, ⟨σ_θθ⟩, anisotropy
6. **Marangoni driving force** (if Layer 4 active)
7. **Vorticity magnitude** (if Layer 4 active)
8. **Coffee-ring index** (if Layer 4 active)
9. **Sim B vs Sim A overlay** (if comparison run)
10. **Conservation diagnostics**: total mass, total momentum, total energy time series

Interactive Plotly version allows hovering, zoom, time-scrubbing.

## Cross-Run Comparison Dashboard

**Output**: `results/sweep_{date}/comparison_dashboard.html`

For each parameter sweep:
- All 5 condition spreading curves on one plot
- Phase diagram (radial-valid regions)
- Final A/A₀ scatter vs parameter
- φ trajectories overlaid

## Configuration

```yaml
visualization:
  tier1_vispy:
    enabled: true
    color_mode: "phi"        # "phi" | "stress" | "rho" | "velocity"
    fps: 30
    resolution: [1280, 720]
  tier2_blender:
    enabled: true
    quality: "production"    # "draft" | "production"
    samples: 256
    color_mode: "stress"
  tier3_dashboard:
    enabled: true
    overlay_experimental_data: true   # for the conditions matching Bare/Pre/Lam4
    interactive_html: true
```

## Implementation Hooks
- `viz/vispy_renderer.py` — frame export
- `viz/blender_render.py` — Blender script
- `viz/dashboard.py` — plotly + matplotlib generation
- `viz/cross_run.py` — sweep aggregation

All called automatically by the post-simulation pipeline (`run.py`).
