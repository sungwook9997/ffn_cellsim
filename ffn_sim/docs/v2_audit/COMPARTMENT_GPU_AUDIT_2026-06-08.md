# Compartment GPU / Native Readiness Audit (2026-06-08)

Per-compartment GPU-main readiness for the explicit single-cell stack, generated
against the registry's declared *Compartment Performance Contracts*
(`cell/compartment_registry.py`) and cross-checked with the build-path survey and
the GPU-main port docs (`GPU_MAIN_PORT_2026-05-31.md`,
`NATIVE_HOT_LOOP_MIGRATION_2026-06-07.md`, `H7_NATIVE_FULLCELL_GO_2026-06-07.md`,
`HOOMD_GPU_BESTPRACTICE_GAP_2026-06-07.md`).

Reproduce the live table:
`PYTHONPATH=. python ffn_sim/scripts/compartment_force_profile.py --table`
Microbench a built cell:
`PYTHONPATH=. python ffn_sim/scripts/compartment_force_profile.py --build --n-fil 200`

## 1. Env-flag conventions (default-OFF, opt-in, parity-validated)

| Flag | Effect | Parity |
|---|---|---|
| `FFN_GPU_DEVICE_BAOAB=1` | BAOAB integrator runs device-resident (`gpu_local_snapshot` + cupy), removing the per-step host round-trip | `tests/test_baoab_device.py` |
| `FFN_GPU_DEVICE_COMPARTMENTS=1` | turgor/membrane/nucleus Custom forces → GPU/native ForceCompute | bit-exact vs CPU (`H7_NATIVE_FULLCELL_GO`) |

Native constrained-BAOAB plugin (`native/`) is a separate, additive opt-in
(`constrained=True`), already validated (γ-parity, bit-parity, nonconv=0).

## 2. P0 — always-on per-step (production hot path)

All five P0 compartments have a planned/landed GPU-resident or HOOMD-builtin path
— **no P0 GPU debt**:

| Compartment | Per-step | GPU path now | Note |
|---|---|---|---|
| cortex | bond/angle | HOOMD-builtin | harmonic bond/angle already GPU-native |
| cytoplasm | gamma_map scale | folded into BAOAB | `FFN_GPU_DEVICE_BAOAB` |
| enclosed_volume (turgor) | `EnclosedVolumePressure` Custom | CuPy / native ForceCompute | `FFN_GPU_DEVICE_COMPARTMENTS`; `enclosed_volume_gpu.py` |
| nucleus | `NucleusConfinement` Custom | CuPy / native | `nucleus_confinement_gpu.py`; stiff lamin sets CFL |
| membrane_surface | `MembraneSurfaceTension` Custom | CuPy / native | `membrane_surface_gpu.py` |

Measured stack (gbook RTX A5000, `H7_NATIVE_FULLCELL_GO`): native integrator
2.43× over cupy; full GPU-main (native int + native compartment ForceCompute)
5.86× / 491 steps/s / ~4.7 d per 2×10⁸ steps. Compartment-force GPU swap alone
1.07×; full GPU-main 2.89×, γ bit-exact.

## 3. P1 — recipe-on per-step (optimization DEBT, CPU today)

Every P1 compartment is a per-step force or updater that reads
`cpu_local_snapshot` — i.e. a device→host→device sync per call. On a GPU-main run
this throttles the C++ force eval (Amdahl; `GPU_MAIN_PORT`,
`HOOMD_GPU_BESTPRACTICE_GAP`). These are **explicit optimization debts** to clear
before any P1 compartment runs at production scale:

| Compartment | Op | Debt |
|---|---|---|
| fa | `IntegrinBondUpdater` + `SubstrateLigandPin` | binder host-sync; P1 microbench required pre-production |
| rigid_ligand_coating | `SubstrateLigandPin` (per-step ligand pin) | host-sync |
| substrate | `SubstrateAnchorSpring` Custom | host-sync; **native ForceCompute candidate** |
| lamellipodium | `WaveMembranePin` + 3 updaters | host-sync + runtime topology growth |
| membrane_load | `MembraneLoadUpdater` | host-sync |
| erm | `ERMHarmonic` Custom | host-sync; **native ForceCompute candidate** (simple radial spring) |

The two `md.force.Custom` P1 compartments with no broad-phase (substrate, erm)
are the cheapest native-ForceCompute wins — they mirror the already-ported
turgor/membrane/nucleus forces (radial/normal springs over a tag range).

## 4. P2 — batch updaters / binding / broad-phase

The binder updaters (crosslinkers, myosin, turnover, fa) and every experimental
compartment are P2. The cortex binders (myosin/xlink) are the standing **P3 binder
GPU-port** in `GPU_MAIN_PORT` (cKDTree → on-device nlist reuse; measured −6…−7% at
native only, so deprioritized behind the integrator/compartment ports). New
experimental compartments inherit the same posture: CPU batch updater with an
explicit debt note, microbench at authoring, GPU port deferred until LIVE.

CFL watch: **microtubules** declare a very high bending stiffness; their bending
CFL may dominate the cell `dt` — gate `dt` before enabling on the full cell
(flagged in the module + registry).

## 5. P3 — geometry / diagnostics

`surface_manifold` (NumPy/scipy geometry, cKDTree broad-phase) and this profiler
are P3 — not in the force budget, CPU acceptable. The manifold must remain
force-free (one optional soft normal confinement, PI-gated) — see
`H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX`.

## 6. Native ForceCompute candidates (ranked)

1. **erm** — single radial harmonic over a tag range; trivial native kernel,
   mirrors `EnclosedVolumePressure`.
2. **substrate** — single anchor spring over ligand tags.
3. **membrane_reservoir** (when LIVE) — breakable tethers; native after the
   kinetics stabilize.
The binder updaters (myosin/xlink/fa/cadherin) are NOT good native-ForceCompute
candidates (topology-mutating, broad-phase); they belong on the cupy/on-device
nlist port track.

## 7. Audit verdict

- **P0 GPU-main: GREEN.** Every always-on compartment has a device-resident or
  builtin path; the full GPU-main stack is measured and γ bit-exact.
- **P1 GPU-main: DEBT, documented.** Six recipe-on compartments are CPU/host-sync;
  each carries an explicit `optimization_debt` in its `GpuReadiness`. None block
  the suspended (Gate-A/B) production path, which uses only P0 compartments.
- **New experimental compartments:** default-OFF, CPU, with debt notes; no GPU
  work required until each clears its pairwise gate and is promoted to LIVE.
