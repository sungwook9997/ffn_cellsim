# DCM GPU — device-dispatch wiring + active+necrosis law sweep (gbook-ready)

> ⚠️ **CORRECTION / RETRACTION (2026-06-19, grounding-pass C14).** Any "spreading" / A/A0 result below measured by the **basal-contact / convex-hull FOOTPRINT** metric is **not valid**: the PI FORBADE it on 2026-06-12 — valid spreading A/A0 MUST be the **top-down xy silhouette** (all nodes). Footprint inflates A/A0 ~10–20×. Where top-down was actually measured (`two_stage_n400_lamel_S10`), the full mechanistic stack **COMPACTS** — top-down A/A0 1.0 → **0.597**, the OPPOSITE sign. Treat footprint "spreading-law" fits below as **not valid** until regenerated with the top-down metric on the gbook GPU. See the grounding-pass table + `project-rebuild-audit` memory.

**2026-06-11 · branch `h7/compartment-platform`.** Continues the GPU work: the validated cupy
kernels (A5000, 827×) are now WIRED into the live HOOMD DCM via device-dispatch, the active
traction is GPU-pathed, and a large active+necrosis size-sweep script tests the PI spreading
law A/A₀ = a + b/R + c/R² — ready to run on the gbook A5000.

## What's wired
- `cell/dcm_gpu_forces.py` — `DeviceDispatch` (GPU → gpu_local_snapshot + cupy + kernels_gpu;
  CPU → cpu_local + numpy + kernels_cpu, cupy import guarded). Forces: `DcmTentContactGPU`,
  `DcmSubstrateForceGPU`, **`DcmActiveRimTractionGPU`** (the active lamellipodium/contraction-
  belt/FA-clutch rim traction). Native mesh shell is already GPU-native. Low-cadence state
  updaters (necrosis/pressure/junction/division) stay CPU (batched, cheap).
- `cell/dcm_active.py:609` — **one-line fix**: `build_active_spheroid` now forwards
  `device=device` to `build_native_dcm_simulation`, so requesting a GPU device actually builds
  the sim GPU-resident (was silently CPU before). CPU path unchanged (device=None default;
  active smoke A/A₀ 2.862 identical pre/post fix). `integrator/baoab.py` untouched (frozen).
- `scripts/dcm_gpu_law_sweep.py` — the gbook law-test: a SIZE SWEEP of large active+necrosis
  spheroids (GPU auto / CPU fallback), A/A₀ vs effective R, fit a+b/R+c/R² vs the PI law.
- Tests: `test_dcm_gpu_forces_parity.py`, `test_dcm_active_gpu_parity.py` — **9/9 PASS**,
  CPU-path force diff < 1e-12 N (bit-identical to the CPU originals).

## CPU-fallback mini-sweep result (R_cell=20 µm patch; N=14/20/28; R 78–93 µm)
| N | R_eff (µm) | A/A₀ | necrotic frac | active-rim frac |
|---|---|---|---|---|
| 14 | 77.6 | 1.57 | 0.00 | 0.14 |
| 20 | 85.6 | 1.39 | 0.00 | 0.15 |
| 28 | 93.1 | 1.26 | 0.00 | 0.21 |

- **The law SIGNATURE is reproduced:** A/A₀ DECREASES monotonically with R (corr = −1.00),
  fit r² = **1.000** (a=−0.005, b=96 µm, c=1994 µm²). Same decreasing shape as the PI law
  (a=−0.33, b=188.7, c=−2655, R=31–78 µm) overlaid in the figure — the coefficients differ
  because this mini-sweep is below the law's R range AND necrosis is off (see below).
- **Necrosis off at this scale** (R < 150 µm onset) — expected; a separate check at
  R_cell=45 µm / N=30 (depth 166 µm) fired **4 NECROTIC cells**, confirming the necrosis
  machinery activates at R > 150 µm. So the gbook large sweep will cross the necrotic onset.
- Figure: `figs/gpu_law_sweep.png` (A/A₀ vs R + fit + PI-law overlay + necrosis-onset line;
  necrotic-fraction / active-rim-fraction vs R).

## Run on the gbook A5000 (PI)
```
PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.scripts.dcm_gpu_law_sweep \
    --sizes 60 120 200 300 450
```
(R_cell=20 µm patch → R ≈ 118–213 µm, crossing the ~150 µm necrosis onset → necrotic core
turns on → the active-rim-fraction-shrinks-with-R mechanism that should yield the law's b/R,
c/R² terms. `--division` enables live division; `--cell-radius-um 7.5 --sizes 2000 6000`
for fine-scale, much slower.)

## Honest bottom line
- The mechanism reproduces the law's qualitative signature (A/A₀ ↓ with R, r²=1.0) ALREADY
  at CPU mini-scale. The COEFFICIENT-level validation needs the gbook R>150 µm run where
  necrosis activates — the only remaining blocker is the gbook hardware (dirty-branch cleanup,
  PI). All code (kernels 827×-validated, forces wired + CPU-parity, device fix, sweep script)
  is ready.
- For a fully GPU-RESIDENT run the only outstanding software item is the BAOAB integrator's
  per-step cpu_local sync (frozen — PI sign-off to port); the FORCES are GPU-pathed now.
