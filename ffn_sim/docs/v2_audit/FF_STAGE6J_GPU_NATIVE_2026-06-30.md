---
kb_record:
  doc_id: FF_STAGE6J_GPU_NATIVE_2026-06-30
  title: FF engine is GPU-native (Warp) — forces + reshape + relaxation on-device, A5000-validated
  authoritative_as_of: 2026-06-30
  supersedes: []
  status: current
---

# FF Stage 6j — the FF cortex runs fully GPU-resident on the gbook A5000

**Date:** 2026-06-30 · **Engine:** `ffn_sim/ff/` · branch `dcm/main`

## Why (PI correction)

The whole point of building FF on **NVIDIA Warp** (instead of reusing Cytosim's C++ core) is that
Warp is **GPU-native + differentiable**: the same kernel runs on CPU (Mac dev) or CUDA (gbook A5000)
by passing `device=`. An earlier FF relaxation/force pass was written in **numpy**, which defeats that
— it pins the work to the CPU and throws away the Warp advantage. PI flagged it
("gpu 포팅 처음에 안해도 된다면서 c++이라서" — Warp should be GPU-native from the start, that's the
reason we chose it over the C++ engine). This stage closes that gap: every FF cortex force and the
relaxation loop are Warp kernels, and the production γ-floor was run on the real A5000.

## What was built (`ff/network_warp.py`)

All one-thread-per-element Warp kernels, float64, device-agnostic (the DCM `network_warp` pattern):

| kernel | physics |
|---|---|
| `link_spring_kernel` | crosslinker / actin axial Hookean spring `k·(L−r0)·û`, `atomic_add ±` |
| `myosin_kernel` | myosin contractile `f_myo·û` pulling a node pair together |
| `turgor_kernel` | radial osmotic `ΔP·(area/node)·r̂` per node (Young–Laplace shell) |
| `axpy_kernel` | explicit overdamped step `x += (dt/γ)·F` |
| `reshape_kernel` | **NF2007 §5.3 inextensibility reshape**, one thread per fiber (disjoint node ranges → no races), restores segment rest-lengths conserving COG |
| `relax_on_device` | the whole relaxation: zero→bending→[link]→[myosin]→axpy→periodic reshape, **all on device buffers, no per-step numpy round-trip** |

The bending kernel (`forces_warp.cytosim_bending_kernel`) and the implicit CG
(`dcm.dcm_warp_implicit.device_cg`) were already Warp; this stage adds the remaining network forces +
the reshape + the driving loop, so nothing in the FF cortex relaxation leaves the GPU.

## Validation

**Bit-parity (CPU).** Each kernel is parity-checked against the original numpy force
(`tests/ff/test_network_warp.py`, 5 tests):
- `link_spring` / `myosin` vs the numpy `_link_spring_force` / `_myosin_force`: ≤1e-12.
- `reshape_kernel` vs numpy `constraints.reshape`: **1.8e-15** (bit-level), segments held.
- `relax_on_device` vs the numpy bending+reshape relax: **2.7e-14 µm** (N=1000), E→0, segments held.

**Real A5000 (gbook, Warp 1.14, CUDA 12.9, sm_86).** FF source rsync'd to `~/ff_scratch` (NOT the
gbook repo — the DCM `/loop` was actively running `dcm_warp_decohesion` on `cuda:0`; we shared the GPU
and never touched its working tree):
- `relax_on_device` GPU vs CPU bit-parity **3e-14 µm**.
- N=4000 cortex: **A5000 0.06 s vs CPU 0.28 s (4.5×)** after kernel warmup; speedup grows with N.

**Production γ-floor, GPU-resident.** `equilibrate(method="device", device="cuda:0")` →
`gamma_floor_run(..., method="device", device="cuda:0")` runs the resting-shell settle on the A5000.
Production point (N=1000 / n_xl=1000 / n_myo=100, f_myo = NMIIA minifilament stall 5.0 pN, 4 real, 0.8 s):

> **γ_active = 1.413e-4 mN/m** (± 2.4e-2 pN/µm) — this **exactly matches** the archived BAOAB-MD
> `g_soft` (~1.4e-4 mN/m). The GPU-resident, MD-free FF reproduces the established actomyosin
> cortical-tension floor (~2500–4600× under the Salbreux 0.35–0.65 mN/m band).

So the GPU path is not just fast — it lands on the same physics number the CPU/MD path established
(γ-floor RESOLVED, `FF_STAGE6D_GAMMA_FLOOR`). The "numpy-not-GPU" gap is closed **end-to-end**.

## Notes / scope

- `relax_on_device` uses **periodic reshape** for inextensibility (the robust path), not the
  per-fiber projector — same resting shell, no singular-prone linear solve in the hot loop.
- Turgor is **not** in the device loop yet: the resting baseline is turgor-free per the
  physiological-baseline rule; the loaded/turgor-pressurised shell (which has no static equilibrium —
  itself the floor signature) uses the implicit host path. A state-dependent on-device turgor
  (centroid + R_mean reduction → Guo closure → `turgor_kernel`) is the next on-device extension.
- gbook sync was **non-destructive** (rsync to scratch). Promoting the gbook *repo* to the current
  commit (Syncthing-share vs git-update) is still a PI decision — the DCM `/loop` runs from that repo.

## Files

`ff/network_warp.py` (kernels + `relax_on_device`), `ff/gamma_floor.py` (`equilibrate`/`gamma_floor_run`
`method="device"`/`device=`), `tests/ff/test_network_warp.py` (5). Commits `d169a70`, `b588c7f`, `5df5c1a`.
