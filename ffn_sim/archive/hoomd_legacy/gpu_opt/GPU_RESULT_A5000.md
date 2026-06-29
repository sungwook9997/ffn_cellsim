# GPU port result — RTX A5000 (gbook), cupy 14.1.0 — ALL KERNELS PASS

The 4 DCM hot-loop kernels (`kernels_cpu.py`) were ported to cupy (`kernels_gpu.py`) by the
Fable model + a multi-agent review (2 defects fixed; numpy-shim parity 11/11) and validated
GPU-resident on the lab RTX A5000. **All kernels bit-parity (max abs ≤ 2.5e-21 N, far inside
the 1e-12 gate); the dominant cell-cell kernel is 827× faster.**

## Bench (measured on the A5000)

`--n-groups 40 --n-per 42` (N=1,680):
```
mesh_pressure  CPU 0.30 ms | GPU 0.64 ms (0.5x) | parity OK (1.06e-22)
plane_well     CPU 0.02 ms | GPU 0.13 ms (0.1x) | parity OK (0)
bond_spring    CPU 0.02 ms | GPU 0.31 ms (0.1x) | parity OK (0)
group_pair     CPU 44.78 ms| GPU 13.04 ms (3.4x)| parity OK (1.03e-25)
ALL KERNELS PASS
```
`--n-groups 200 --n-per 162` (N=32,400):
```
mesh_pressure  CPU 6.73 ms   | GPU 0.71 ms  (9.5x)  | parity OK (2.54e-21)
plane_well     CPU 0.36 ms   | GPU 0.18 ms  (1.9x)  | parity OK (0)
bond_spring    CPU 0.05 ms   | GPU 0.35 ms  (0.1x)  | parity OK (0)
group_pair     CPU 11706 ms  | GPU 14.14 ms (827.7x)| parity OK (4.24e-22)
ALL KERNELS PASS
```

## What this unblocks
`group_pair` (the cell-cell soft-core + adhesion neighbour search) was the #1 bottleneck
(O(N²); 11.7 s/call at N=32k on CPU). The cupy uniform-grid neighbour list makes it
**14 ms (827×)** — so a LARGE many-cell spheroid (the scale needed to show the necrotic
core >500 µm and fit A/A₀=a+b/R+c/R² across R 31-78 µm) is now computationally reachable.
The small-N micro-kernels (plane_well/bond_spring) stay launch-overhead-bound — expected,
not a port issue.

## Port design (kernels_gpu.py)
K1-K3 direct cupy with `cupyx.scatter_add` + a manual component-wise cross product
(FP-identical to np.cross). K4 = fully-vectorized uniform-grid neighbour list: bin on a
grid of edge r_cut, gather 27-cell-stencil candidates via sorted-cell-id `searchsorted` +
cumsum segment expansion (deliberately avoiding `cp.repeat` with ndarray repeats, which
the installed cupy rejects), `i<j` dedup, then the exact reference per-pair math with
Newton-3 ± scatters.

## NEXT — Phase 2 (Opus): wire into the HOOMD DCM Custom forces
The validated cupy kernels now back the HOOMD `md.force.Custom` classes via
`gpu_local_snapshot` + `gpu_local_force_arrays` (GPU-resident), CPU path bit-identical:
`DcmTurgorForce`/`DcmSubstrateForce` (dcm.py), `FaClutchForce` (dcm_ecm.py),
`DcmCellCellAdhesion` (dcm_prolif.py). Then the large-spheroid production run on the A5000.
Files in `ffn_sim/gpu_opt/` + on the gbook `~/kernel_port/` (HANDOFF.md, test_shim_parity.py).
