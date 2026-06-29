# CUDA/cupy port task — 4 array kernels (copy-paste)

Pure numerical-computing task: port four small NumPy array kernels to **cupy** for an
NVIDIA RTX A5000, keeping the math bit-identical, and pass the included parity + benchmark.

**Environment.** Work from the project root — the directory that contains the `ffn_sim/`
package (`cd` there; run python from it so `import ffn_sim` resolves). Activate the env:
`source ~/miniconda3/etc/profile.d/conda.sh && conda activate ffn_sim`. `cupy`
is available on the A5000 box. (Off the GPU box, still write the cupy code; the harness is
import-guarded so it runs CPU-only locally.)

**Reference.** `ffn_sim/gpu_opt/kernels_cpu.py` has 4 pure-NumPy kernels + a synthetic
data generator. They are plain vectorized linear algebra over float64 arrays.

## Task — write `ffn_sim/gpu_opt/kernels_gpu.py`

Expose the SAME 4 function names + signatures as `kernels_cpu.py`, operating on **cupy**
arrays (accept numpy or cupy in, return cupy out):

1. `mesh_pressure_forces(pos, tris, tri_group, n_groups, V0, p0, K)` — per-facet normal
   force from a per-group divergence-theorem volume. Reference uses `np.add.at` for the
   per-group volume reduction and the per-vertex force scatter; on cupy use
   `cupyx.scatter_add` (or a segment-sum). Math identical.
2. `plane_well_forces(pos, z0, k, w, mask=None)` — z-only capped-harmonic; elementwise.
3. `bond_spring_forces(pos, bonds, r0, k)` — per-bond Hookean; scatter ±force to the two
   endpoints with `cupyx.scatter_add`. Handle `bonds.shape[0]==0`.
4. `group_pair_forces(pos, group_id, sigma, r_cut, k_core, f_well0)` — **the dominant cost**
   (the NumPy reference is brute-force O(N²), measured ~21 ms/call at N≈1.3k). **Replace it
   with a cupy grid/neighbour list**: bin particles on a uniform grid of bin edge `r_cut`,
   gather candidate pairs within `r_cut`, drop pairs with equal `group_id` and any
   `group_id < 0`. The surviving-pair math (soft-core repulsion `r<sigma`, linear attractive
   well `sigma<=r<r_cut`, Newton-3 ± scatter) must stay bit-identical to the reference.

**Constraints:** pure NumPy→cupy port — do NOT change the math. Deterministic, float64,
exact Newton-3 symmetry on pair/bond/facet scatters.

## Verify — must pass

`ffn_sim/gpu_opt/bench.py` compares your kernels to the reference + times both:
```
python -m ffn_sim.gpu_opt.bench --n-groups 40 --n-per 42       # parity (rtol 1e-9)
python -m ffn_sim.gpu_opt.bench --n-groups 200 --n-per 162     # large-N speedup
```
- Every kernel: **parity OK** (max abs error ≲ 1e-12; rtol 1e-9).
- Large system: cupy faster than NumPy, especially `group_pair`.

## Deliverables
1. `ffn_sim/gpu_opt/kernels_gpu.py` (cupy, parity-passing).
2. The `bench.py` output (parity table + speedups).
3. One paragraph: what you did + the measured A5000 speedup at `--n-groups 200 --n-per 162`.

Do not modify any other files.
