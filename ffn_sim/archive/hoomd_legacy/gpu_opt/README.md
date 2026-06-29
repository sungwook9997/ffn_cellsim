# gpu_opt/ — DCM hot-loop GPU-optimization bundle (for the Fable optimizer)

Self-contained, **biology-free** extraction of the DCM many-cell spheroid's per-step force
kernels, so they can be ported to GPU (cupy / `gpu_local_snapshot`) for RTX-A5000 production
WITHOUT needing any cell-biology knowledge. The DCM physics currently runs through HOOMD
`cpu_local_snapshot` (a GPU→CPU sync every step) — this bundle is the port surface.

## Files
- **`FABLE_GPU_OPT_PROMPT.md`** — the copy-paste task for the Fable optimizer (Phase 1: port
  4 kernels to cupy + parity; Phase 2: wire them into the HOOMD Custom forces via
  `gpu_local_snapshot`, CPU path bit-identical).
- **`dcm_kernels_cpu.py`** — pure-numpy REFERENCE kernels (the math, extracted verbatim from
  the in-tree Custom forces) + a synthetic many-cell data generator. No HOOMD, no biology.
- **`parity_bench.py`** — compares the GPU kernels (to be written as `dcm_kernels_gpu.py`)
  against the numpy reference on synthetic arrays + times both. Acceptance = parity rtol 1e-9
  + GPU faster on large N.
- `dcm_kernels_gpu.py` — **written by the optimizer** (cupy).

## The 4 kernels (array contracts; all SI, float64, return per-node force (N,3))
| kernel | inputs | source class |
|---|---|---|
| `turgor_forces` | pos(N,3), faces(F,3), face_cell(F,), n_cells, V0, turgor_dP0, K_vol | DcmTurgorForce (dcm.py / dcm_ecm.py / dcm_prolif.py) |
| `substrate_forces` | pos, z0, k_well, rng, actin_mask | DcmSubstrateForce (dcm.py) |
| `fa_clutch_forces` | pos, pairs(M,2), r0(M,), k_fa | FaClutchForce (dcm_ecm.py) |
| `cell_cell_adhesion_forces` | pos, cell_of_node(N,), sigma, r_cut, k_core, F_adh0 | DcmCellCellAdhesion (dcm_prolif.py) |

## Measured bottleneck (CPU reference, 1260 nodes)
`cell_cell_adhesion ≈ 21 ms/call` (O(N²)) ≫ turgor 0.2 / substrate 0.02 / fa_clutch 0.01 ms.
→ **the cell-cell adhesion neighbour search is the #1 GPU target** (use a cupy cell list).

## Run
```
source ~/miniconda3/etc/profile.d/conda.sh && conda activate ffn_sim
cd /Users/sw1/ffn_cellsim-platform
python -m ffn_sim.gpu_opt.parity_bench --n-cells 40 --n-per-cell 42      # parity
python -m ffn_sim.gpu_opt.parity_bench --n-cells 200 --n-per-cell 162    # large-N speedup
```

## Out of scope for Fable
`integrator/baoab.py` (the BAOAB GPU-main port is a separate PI-frozen effort) and all
biology/config files. Native `md.bond/angle/pair` are already GPU-native — leave as-is.
