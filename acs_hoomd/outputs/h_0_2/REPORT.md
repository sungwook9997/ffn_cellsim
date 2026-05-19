# Phase 0.2 — HOOMD-blue environment verification

**Date**: 2026-05-19
**Branch**: `v2/foundation`
**Plan reference**: Notion Build Plan v2 §2 Unit H.0.2 — <https://www.notion.so/365120daec5d81799efefcf078f2039e>

## Environment

| Layer | Detail |
| --- | --- |
| Hardware | M1 Max MacBook (10-core), macOS 26.4.1, `arm64` |
| Conda | `conda 26.1.1`, env name `acs_hoomd` |
| Python | 3.13.13 |
| HOOMD-blue | **7.0.1** (cpu build, MD module enabled, MPI off, GPU off, float64/float32) |
| Companion deps | gsd 5.0.1, freud 3.5.0, numpy 2.4.6, scipy 1.17.1, matplotlib 3.10.9, pytest 9.0.3, pyyaml 6.0.3 |

Install one-liner (executed):

```bash
conda create -n acs_hoomd python=3.13 -y
conda install -n acs_hoomd -c conda-forge "hoomd=7.0.1=cpu*" gsd freud numpy scipy matplotlib pytest pyyaml -y
```

**Note on Plan vs reality**: Plan §2 H.0.2 specified HOOMD ≥ 4.4. Current conda-forge stable on `osx-arm64` is 7.0.1; we are jumping straight to the current series. The 4 → 7 API change retired `hoomd.version.cuda_built` (now `gpu_enabled`), the rest of the surface we use (`bond.Harmonic`, `angle.Harmonic`, `methods.Brownian`, `methods.Langevin`, `Integrator`) is intact.

`fresnel` was deferred — not required by the polymer smoke test, and the v2 visualization pipeline (Plan §8 cross-cutting) will be specified separately.

## Polymer smoke test

Script: [acs_hoomd/scripts/hoomd_polymer_sanity.py](../../../acs_hoomd/scripts/hoomd_polymer_sanity.py)

| Setting | Value |
| --- | --- |
| Beads | 100 |
| Bonds | 99 (harmonic, k=100, r0=1.0) |
| Angles | 98 (harmonic, k=5, t0=π) |
| Thermostat | Langevin, kT=1.0, γ=1.0 |
| dt | 0.005 |
| Box | cubic, side = max(2·N·r0, 50) |

### Results

| Run | Wall time | Steps/s | Final KE | Final PE | T_eff |
| --- | --- | --- | --- | --- | --- |
| 5 000 steps | 0.15 s | 33 194 | 160.064 | 122.037 | 1.067 (target 1.0) |
| 50 000 steps | 0.56 s | **89 297** | 165.846 | 139.835 | 1.106 (target 1.0) |

`T_eff` is computed as `(2/3) · KE / N` — instantaneous, single-frame. Equilibrium fluctuation ~1/√N ≈ 10 % is consistent with the ~7–11 % residual; a time-averaged read over the latter half of the trajectory would tighten to a few %. For sanity-test purposes the integrator, bonds, angles, and thermostat are all working.

### Extrapolation

89 297 steps/s gives **~11 s wall-time for 1 M timesteps** on this 100-bead chain — well inside the Plan's "<30 s on A100 reference" target. (For small systems CPU often beats GPU due to launch overhead; the GPU win arrives at ≳ 10⁴ particles.)

### Artifacts (in this directory)

- `initial.gsd` — initial straight-chain configuration.
- `trajectory.gsd` — frames every 1 000 steps from the 50 k-step run.

## Phase 0.2 verdict

**PASS.** Environment is reproducible (single `conda install`), HOOMD imports cleanly, the canonical HOOMD primitives we need for Phase 1 (`bond.Harmonic`, `angle.Harmonic`, `Langevin`, `Brownian`) all instantiate, and a 100-bead Langevin polymer reaches a sensible equipartition temperature. Ready for Phase 0.3 (AFINES algorithm extraction) and Phase 0.4 (`acs_hoomd/` skeleton commit — already laid down in this PR).
