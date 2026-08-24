# FF Stage 6R — comprehensive GPU↔CPU bit-parity (the A5000 production path)

**Date:** 2026-07-01  **Engine:** FF (Warp, GPU-native)  **Hardware:** gbook RTX A5000  **Branch:** dcm/main

FF is GPU-native (the reason it chose Warp over Cytosim C++): the production path is `device="cuda:0"`.
This stage validates that path as **bit-faithful to the `device="cpu"` reference across EVERY on-device
kernel and integrator at native scale** — the standing guarantee behind every A5000 γ-floor number.

## Sweep result (identical input, CPU vs cuda:0, max|Δ|)

| kernel / integrator | scale | max \|Δ\| | note |
|---|---|---|---|
| `reshape` (inextensibility projection) | N=494 802 nodes | **2.2e-16** | exact-arithmetic identical (machine ε) |
| `link_spring_force` | 70 686 crosslinks | **2.2e-11** | float64 atomic-add ordering only |
| `myosin_force` | 442 dipoles | **1.8e-15** | — |
| `branch_angle_force` (Arp2/3) | 100 triples | **1.8e-16** | — |
| `relax_on_device` (bending+xl+reshape, 400 steps) | 140 000 nodes | **9.8e-15** | per-step atomic Δ does NOT amplify |
| `simulate_loaded_shell_on_device` (turgor, 4000 steps) | 35 000 nodes | **4.5e-14** | V/V0 & ΔP identical to 4 dp |
| `simulate_turnover_on_device` (Hand KMC) | 500 motors | engaged frac Δ = **0.0000** | stochastic → statistical parity (same Bell steady state) |

**Reading:** the deterministic kernels are bit-identical up to float64 atomic-add reduction ordering
(the `link_spring` 2.2e-11 is the worst case, from parallel `wp.atomic_add` summation order — physically
zero). The two iterated integrators stay at ~1e-14 even after thousands of steps (no error amplification).
The KMC turnover is stochastic (independent `wp.rand` streams) so it is validated STATISTICALLY — the
engaged fraction converges to the same analytic Bell steady state (k_on/(k_on+p_off)) on both devices.

## Codified

`tests/ff/test_gpu_parity.py` — runs the force-kernel + relax-integrator parity on the A5000, **skips
offline** (`wp.is_cuda_available()` is False on the dev Mac) — the same skip-if-absent pattern as the
Cytosim parity oracle. On the A5000: **4 passed in 1.5 s**. So the GPU production path is now an
always-on (GPU-present) regression guard, not a one-off check.

## Sanity gate

- **Boundary:** CPU and cuda:0 are the SAME numerical algorithm; any Δ above float64 atomic ordering would be a kernel bug. ✓ (all ≤ 2.2e-11)
- **No amplification:** 400- and 4000-step integrators stay ~1e-14 — the per-step Δ is not a growing mode. ✓
- **Stochastic handled honestly:** KMC is statistical (not bit) parity — engaged fraction = analytic Bell steady state on both. ✓
- **Native scale:** validated at the lit-faithful native cortex (N=70686 / 494 802 nodes), not a toy. ✓

Related: [[project-ff-gpu-native]] (A5000 ~3e-14 core-kernel parity, now extended to ALL kernels + integrators),
FF_STAGE6J_GPU_NATIVE, FF_STAGE6Q (the lit-faithful densities this was run at).
