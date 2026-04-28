# 11 — Performance Protocol & GPU Portability

## Performance Targets

| Mode | N points | Sim time | Wall-clock target | Hardware |
|---|---|---|---|---|
| Pilot | 1,000 | 4–8 hr | **≤ 30 min** | RTX A5000 |
| Production | 5,000 | 80 hr | 5–15 hr | RTX A5000 |
| Production large | 10,000 | 80 hr | 20–40 hr | RTX A5000 (or 4090×2 split) |

Pilot must complete in 30 min for fast iteration during development. If it does not, reduce N or simplify before optimizing.

## Optimization Priorities (in order)

1. **Use Taichi's official MLS-MPM example as starting point** — already optimized.
2. **Adaptive time step** (CFL-based) — reduces step count by 3–5×.
3. **Spatial sorting / cell list** — keeps locality for GPU caching.
4. **Avoid Python-level loops** — all hot paths in `@ti.kernel`.
5. **Implicit time integration** for stiff terms (cortex pressure) — 10× speedup but complex; consider for production only.
6. **Multi-GPU split** (4090×2) — for large production runs.

Avoid premature optimization. Profile first (Stage 1a).

## GPU Portability Mandate

The same code must run on:
- RTX A5000 (primary, CUDA, 24 GB)
- RTX 4090 ×2 (secondary, CUDA, 24 GB each — single GPU per process initially)
- Google Colab (T4 / A100, fallback)
- CPU (debugging fallback)

### Implementation

```python
# In every entry point
import os
import taichi as ti

backend_str = os.environ.get("ACS_GPU_BACKEND", "auto")
if backend_str == "cuda":
    ti.init(arch=ti.cuda)
elif backend_str == "vulkan":
    ti.init(arch=ti.vulkan)
elif backend_str == "cpu":
    ti.init(arch=ti.cpu)
elif backend_str == "auto":
    ti.init(arch=ti.gpu)  # Taichi auto-selects
else:
    raise ValueError(f"Unknown ACS_GPU_BACKEND: {backend_str}")
```

### Environment Detection
Add `scripts/check_env.py` that prints:
- Detected GPU
- VRAM available
- CUDA / Taichi version
- Suggested config (pilot or production)

## Benchmarking Protocol (Stage 1a Deliverable)

Run on each available hardware:
- Pilot config (1,000 points, 4 sim hr)
- Record: total wall-clock, mean step time, memory usage, GPU utilization

Record in `benchmarks/{hardware}_{date}.json`:
```json
{
  "hardware": "RTX A5000",
  "n_points": 1000,
  "sim_hours": 4,
  "wall_clock_seconds": 1620,
  "mean_step_time_ms": 0.45,
  "peak_vram_gb": 4.2,
  "mean_gpu_util_pct": 78
}
```

This data informs Production-config tuning.

## Output Storage

Each run produces:
- HDF5 snapshots: ~100 MB – 1 GB per run
- Frame PNGs: ~50–200 MB
- Blender MP4: ~100–500 MB
- Dashboard HTML + figures: ~10 MB

Total per run: ~250 MB – 2 GB

Workflow on Windows workstation:
- Active runs in `D:/ActiveCellSim/results/`
- After analysis, move to NAS or Dropbox sync
- Mac side: only download dashboard.html + summary metrics for review

## Memory Budget (RTX A5000, 24 GB VRAM)

For 5,000 material points:
- Material point arrays (position, velocity, F, σ, φ, ρ, ...): ~5 KB/point × 5,000 = 25 MB
- Background grid (128³): 128³ × ~10 fields × 4 bytes = 80 MB
- Particle-grid transfers, intermediates: ~200 MB
- Taichi runtime overhead: ~500 MB

Total: < 1 GB — well within budget. Even 100,000 points fits comfortably; we are CPU-bound, not memory-bound.

## SSH Workflow (Headless)

### Long-running simulation invocation
```bash
# On Windows workstation (via SSH from Mac)
cd D:/ActiveCellSim
python run.py --config configs/production.yaml --run-id sweep_p1_seed42 &
disown
```

Or use `tmux` / `screen` for session persistence. The PI's existing `train_win` alias pattern can be adapted.

### Monitoring
```bash
# From Mac
ssh win "tail -f D:/ActiveCellSim/results/sweep_p1_seed42/sim.log"
ssh win "nvidia-smi"
```

### Result retrieval
```bash
# From Mac, after completion
scp -r win:D:/ActiveCellSim/results/sweep_p1_seed42/dashboard.html ~/Desktop/
# Don't pull large HDF5 unless needed for offline reanalysis
```

## Failure Modes & Recovery

- **OOM on GPU**: reduce N or grid resolution; run on 4090
- **Numerical instability**: check Δt, reduce CFL, verify initial conditions sane
- **Blender render fails**: regenerate from saved HDF5; render is a separate process
- **Long-tail simulation**: kill if stuck; debug with profiler dump

Add health check pings: every 1 sim-hour, dump a "still alive" log line. If silent for > 30 min wall-clock, suspect hang.
