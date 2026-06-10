# Fable GPU-optimization prompt — DCM hot-loop kernels (copy-paste to Fable)

> Self-contained GPU / numerics task. **No biology knowledge required** — these are pure
> array-math kernels. Goal: make the DCM many-cell spheroid runnable as GPU-resident
> production on an NVIDIA RTX A5000 (the lab "gbook").

---

You are optimizing the per-step hot-loop force kernels of a HOOMD-blue particle simulator
for NVIDIA GPU. The reference implementations are pure NumPy in
`ffn_sim/gpu_opt/dcm_kernels_cpu.py`; right now the real simulation computes them via
HOOMD `cpu_local_snapshot`, which forces a **GPU→CPU sync every step** and kills GPU
throughput. Your job is to make them **GPU-resident (cupy)** and verify bit-parity.

**Environment.** Repo root `/Users/sw1/ffn_cellsim-platform`. Run python FROM THE ROOT so
`import ffn_sim` works. Activate: `source ~/miniconda3/etc/profile.d/conda.sh && conda
activate ffn_sim`. Target hardware: RTX A5000 (CUDA), `cupy` available there. (If you are
NOT on the GPU box, still write the cupy code + make the parity harness import-guard so it
runs CPU-only locally — the lab runs the GPU parity on the A5000.)

## Phase 1 — port the 4 kernels to cupy (REQUIRED)

Create `ffn_sim/gpu_opt/dcm_kernels_gpu.py` exposing the SAME 4 function names + signatures
as `dcm_kernels_cpu.py`, but operating on **cupy** arrays (accept numpy or cupy in, return
cupy out — the harness moves data to device and back):

1. `turgor_forces(pos, faces, face_cell, n_cells, V0, turgor_dP0, K_vol)` — per-cell exact
   triangulated enclosed-volume → outward face-normal force. The numpy uses `np.add.at` for
   the per-cell volume reduction and the per-node force scatter; on cupy use
   `cupyx.scatter_add` (or segment-sum) for both. Math identical.
2. `substrate_forces(pos, z0, k_well, rng, actin_mask=None)` — z-only capped-harmonic well.
   Trivially elementwise on cupy.
3. `fa_clutch_forces(pos, pairs, r0, k_fa)` — per-bond harmonic; scatter ±force to bond
   ends with `cupyx.scatter_add`. Handle `pairs.shape[0]==0`.
4. `cell_cell_adhesion_forces(pos, cell_of_node, sigma, r_cut, k_core, F_adh0)` — the
   DOMINANT cost for many cells. The reference is brute-force O(N²); **replace it with a
   cupy cell/neighbour list** (uniform grid bin by `r_cut`, gather candidate pairs within
   `r_cut`, drop same-cell pairs `cell_of_node[i]==cell_of_node[j]` and dormant `<0`) so it
   scales to ~10⁴–10⁵ nodes. The PER-PAIR math (soft-core repulsion r<sigma, linear
   adhesive well sigma≤r<r_cut, Newton-3 ± scatter) is exactly the reference — keep it
   bit-identical for the pairs that survive the cutoff.

**Constraints:** do NOT change any physics/math — this is a pure numpy→cupy port.
Deterministic (no RNG in these kernels). Float64. Match Newton's-3rd-law symmetry exactly.

## Verify — parity + benchmark (REQUIRED, must pass)

`ffn_sim/gpu_opt/parity_bench.py` already compares your GPU kernels to the CPU reference on
synthetic many-cell arrays and times both. Acceptance:
```
python gpu_opt/parity_bench.py --n-cells 40 --n-per-cell 42      # parity (rtol 1e-9)
python gpu_opt/parity_bench.py --n-cells 200 --n-per-cell 162    # large-N speedup
```
- Every kernel must report **parity OK** (max abs error ≲ 1e-12 N; rtol 1e-9).
- On the large system the cupy kernels must be **faster** than CPU (esp. cell_cell_adhesion).
Report the parity table + the speedups.

## Phase 2 — wire cupy kernels into the HOOMD Custom forces (do if Phase 1 passes)

Re-wire these in-tree `md.force.Custom` / `hoomd.custom.Action` classes to compute on the
GPU via `gpu_local_snapshot` + `gpu_local_force_arrays` (calling your cupy kernels) when the
device is a GPU, with a **device flag that keeps the CPU path bit-identical** (so CPU runs
are unchanged). Follow the project's established GPU-main pattern — see
`docs/v2_audit/GPU_MAIN_PORT_PHASE1_2026-06-01.md`,
`GPU_MAIN_PORT_PHASE1b_BINDER_READINESS_2026-06-01.md`,
`GPU_OPTIMIZATION_ROADMAP_2026-06-01.md`, and the validated cupy path in `native/`.

Classes to re-wire (replace the `with ...cpu_local_snapshot` / `cpu_local_force_arrays`
blocks with the gpu_local equivalents calling the cupy kernels; keep the CPU branch):
- `ffn_sim/cell/dcm.py`: `DcmTurgorForce.set_forces`, `DcmSubstrateForce.set_forces`
- `ffn_sim/cell/dcm_ecm.py`: `FaClutchForce.set_forces` (+ the `FaCatchSlipUpdater.act`
  topology rebind is inherently a host op — leave it batched/low-cadence, do NOT block on it)
- `ffn_sim/cell/dcm_prolif.py`: `DcmCellCellAdhesion.set_forces`, the turgor force, and the
  `ProliferationUpdater.act` (division `set_snapshot` is host — keep low cadence)

Do NOT touch `integrator/baoab.py` — the BAOAB integrator GPU-main port is a SEPARATE,
PI-frozen effort (coordinate, don't duplicate). Do NOT change any other files, configs, or
the physics. Keep `md.bond.Harmonic` / `md.angle.Harmonic` / `md.pair.LJ` (already
GPU-native) as-is.

## Deliverables
1. `ffn_sim/gpu_opt/dcm_kernels_gpu.py` (cupy kernels, parity-passing).
2. The parity+bench output (table + speedups).
3. (Phase 2) the gpu_local re-wiring of the listed Custom classes, CPU path bit-identical.
4. A one-paragraph note of what you changed + the measured A5000 speedup on the large system.
