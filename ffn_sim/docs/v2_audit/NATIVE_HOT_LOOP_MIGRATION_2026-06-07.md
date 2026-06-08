# Native Hot-Loop Migration Plan — 2026-06-07

## Decision

Production H.7 cannot reach the required wall time by removing snapshots alone.
The next architecture move is a **VRAM-resident hot loop**:

- Python remains the orchestration layer: configuration, construction, launch,
  sparse checkpointing, plotting, and validation.
- Per-step physics moves out of Python: BAOAB, fixed-pool binder forces,
  custom external forces, and manifold contact must run in HOOMD native
  C++/CUDA operations or equivalent compiled kernels.
- Global `get_snapshot()` / `set_snapshot()` stays out of the production
  timestep loop. It is allowed for construction, debug, checkpoint, and final
  measurement.

This does **not** mean rewriting the whole simulator in C++. It means removing
Python re-entry and host-memory traffic from the per-step path.

## Evidence From GBook

Hardware: RTX A5000 Laptop GPU, current driver power limit 80 W (max reported
110 W). GEMM reaches 79-80 W and 100% SM utilization, so the GPU can draw its
current cap when given a dense GPU workload.

H.7 Gate-A speed probes use:

- `FFN_GPU_DEVICE_BAOAB=1`
- `MEASURE_MODE=final`
- `FREEZE_XLINKS=1`
- `SNAPSHOT_FREE_MYOSIN=1`
- `MYOSIN_WARM_BIND_TICKS=0`

These are speed diagnostics, not Gate-A correctness verdicts, because dynamic
myosin topology is frozen.

| n_fil | particles (approx) | BAOAB fast steps/s | GPU power / SM | Native dt=0 ceiling | Notes |
|---:|---:|---:|---|---:|---|
| 150 | 5,498 | 2,123 | 38-41 W / 56-65% | 2,861 | Small workload, Python BAOAB visible. |
| 1,000 | 17,660 | 1,543 | 45-48 W / 68-75% | 1,864 | Ratified mesoscopic production scale. |
| 3,000 | ~45k | 1,036 | 56-59 W / 80-86% | not rerun here | GPU occupancy improves but wall time worsens. |
| 5,000 | 75,068 | fails first batch | out-of-bounds | 719 | Stability wall before speed wall. |
| 38,000 | 550,016 | fails first 5k steps | out-of-bounds | 122 | Native biological count is not production-feasible as explicit topology. |

The native dt=0 ceiling is measured with built-in HOOMD bond/angle/LJ force
scheduling and no movement (`dt=0`, `kT=0`). It is a timing ceiling for the
current explicit topology, not a physical run.

## Native Plugin Probe — Same Day

A native HOOMD plugin scaffold now exists under
`native/ffn_hoomd_plugin/`. It is intentionally outside the normal Python
package build so experimental compiled operations cannot regress ordinary test
runs.

Validated on gbook / HOOMD-blue 7.0.1:

- `_ffn_native.build_info()` compile/load probe imports from Python.
- `FFNNoOpUpdater` attaches as a true HOOMD C++ updater and is called exactly
  once per timestep on CPU and GPU.
- `FFNPositionKickUpdater` mutates `ParticleData::getPositions()` directly.
  CPU runs use a host `ArrayHandle`; GPU runs use a device `ArrayHandle` and a
  CUDA/HIP kernel.

Smoke results:

| Probe | Device | Steps | Result |
|---|---|---:|---|
| native no-op updater | GPU | 16 | count = 16, last timestep = 15 |
| native position kick | GPU | 8 | x_final = 1.0, count = 8 |
| native position kick | CPU | 8 | x_final = 1.0, count = 8 |

Native updater overhead probe (`benchmark_updater_overhead.py`, GPU, 200k
steps, one particle):

| Mode | Steps/s | Meaning |
|---|---:|---|
| no updater | 859,979 | HOOMD empty run-loop ceiling |
| Python no-op `CustomUpdater` | 763,130 | Empty Python callback alone is not the whole H.7 wall |
| native no-op updater | 829,522 | C++ scheduler attachment cost is near empty-loop ceiling |
| native position-kick CUDA kernel | 330,985 | One native device kernel per step is still ~150x faster than H.7 Gate-A |

Native position-kick scaling (`dx=0`, GPU, 100k steps):

| Particles | Steps/s |
|---:|---:|
| 17,660 | 189,588 |
| 75,068 | 99,811 |
| 550,016 | 9,093 |

Interpretation: the present H.7 speed wall is not simply that Python exists at
launch time. It is the per-step combination of Python re-entry, custom force
bookkeeping, local snapshot/context paths, and multiple Python-launched device
kernels. A compiled HOOMD updater can attach, run, and mutate GPU-resident
state at the right order of magnitude for the `~1000` effective-filament
production scale.

## Interpretation

1. Snapshot removal was necessary but not sufficient.
2. At `n_fil=150`, compiled BAOAB could recover some overhead, but native
   HOOMD ceiling is only ~1.35x above the current fast path.
3. At `n_fil=1000`, Python BAOAB overhead is smaller; the force/topology stack
   dominates.
4. Increasing filament count improves GPU occupancy but reduces steps/s.
5. Native biological filament count (`~38k`) as explicit bead/bond topology is
   outside the required production time budget.

The ratified `~1000` effective-filament scale is therefore not a convenience. It
is the workable production scale for current hardware unless the topology and
interaction representation change.

## Migration Architecture

### Stage 1 — Native BAOAB Method

Target:

- Replace `hoomd.custom.Action` BAOAB per-step callback with a compiled HOOMD
  operation.
- Preserve the Leimkuhler-Matthews update:
  `r += (F/gamma) dt + sqrt(kT/gamma) (W_n + W_{n-1}) dt`.
- Keep per-tag previous random vectors in device memory.
- Keep CPU fallback and the current Python implementation for validation.

Validation gates:

- CPU vs existing Python CPU bit/path parity where deterministic.
- GPU vs CPU distribution tests: diffusion, drift, harmonic equilibrium.
- Full H.7 smoke: no regression in `n_fil=150` and `n_fil=1000`.

Expected effect:

- Removes per-step Python callback from BAOAB.
- Helps, but cannot by itself reach 28k-56k steps/s because native force-stack
  ceilings are lower at realistic topology sizes.

### Stage 2 — Fixed-Pool Binder Forces

Current binder path:

```text
get_snapshot()
host cKDTree / Python loops
rewrite HOOMD bonds.group/typeid
set_snapshot()
```

Target:

```text
fixed attachment pool in VRAM:
  active flag
  head tag
  actin tag or local manifold/mesh coordinate
  grip_s
  bound filament / bead position

CUDA batch kernels:
  Bell-Evans/Pereverzev break
  bind candidate selection
  grip-walk update
  force accumulation
```

HOOMD bond count should not change during production. Binding state changes
toggle device-side records, not HOOMD topology.

Validation gates:

- Single-head Bell-Evans off-rate oracle.
- Static attachment pool reproduces HOOMD harmonic force for a known pair.
- Grip-walk dipole sign tests from the existing Stage-1 suite.
- H.7 active-tension trend compared against current exact topology path at
  small `n_fil`.

### Stage 3 — Native Custom Forces

Port production `md.force.Custom` forces that run every step:

- ERM tether
- enclosed-volume pressure / turgor
- membrane surface tension / reaction forces
- future manifold confinement/contact

Prefer built-in HOOMD primitives when they exactly match the physics. Otherwise
use compiled force kernels.

Validation gates:

- Force sign and units against existing Python custom force tests.
- Energy/pressure sanity where applicable.
- No contribution to active cortical-tension estimator unless physically
  intended.

### Stage 4 — Manifold / Mesh Contact

The explicit bead network remains the mechanistic state, but search and contact
are restricted by a mesh/manifold layer:

- patch / triangle / shell coordinate stored in VRAM
- local candidate lists for binder/contact
- normal-only confinement and contact as compiled kernels

This is the path that makes higher effective complexity possible without
falling back to all-to-all 3D host searches.

## Immediate Engineering Tasks

1. Keep `h7_hotloop_scaling.py` as the reproducible benchmark harness.
2. Open a native plugin directory outside the Python package build so it cannot
   break normal tests. **Done:** `native/ffn_hoomd_plugin/`.
3. Implement only a compile/load probe first. **Done:** `_ffn_native`.
4. Prove native HOOMD updater attach and GPU state mutation. **Done:**
   no-op updater + position-kick CUDA kernel.
5. Implement native BAOAB as the first real operation.
6. Add fixed-pool binder data structures in Python first, then move kernels.
7. Only after parity gates pass, wire native mode behind an opt-in env var.

## Non-Negotiables

- No gate loosening.
- No empirical speed hacks that change physics silently.
- Dynamic topology replacement must be validated against the exact topology
  path at small scale before becoming production.
- `n_fil=38,000` explicit bead/bond topology is a diagnostic stress test, not
  the production baseline on the A5000.
