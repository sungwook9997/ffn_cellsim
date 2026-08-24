r"""Native-scale equivalence: device plane sums vs the host estimator, on the real cortex.

Answers a question the module's author pre-registered BEFORE any native run: the per-block error grew
14.7x from 1 to 31 elements/block, and at native (~4,000 per block) a linear extrapolation crosses the
1e-11 tolerance while a sqrt one does not. This differences the two paths on the same arrays.
"""
import sys, time
import numpy as np, warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.build import build_all
from aleph.world.observe_gamma import (gamma_planes_device, segment_offsets, plane_force_sums,
                                       actin_axial_tension, _PLANE_BLOCKS)

wp.init()
dev = wp.get_device("cuda:0")
assert dev.is_cuda, dev
print(f"device {dev}  warp {wp.config.version}  _PLANE_BLOCKS={_PLANE_BLOCKS}", flush=True)

arena = WorldArena(capacity={
    Kind.NODE: 12_000_000, Kind.SEGMENT: 12_000_000, Kind.ANGLE3: 12_000_000,
    Kind.ANGLE4: 2_000_000, Kind.FACE: 1_500_000, Kind.STRAND: 200_000,
    Kind.BOND: 1_000_000, Kind.GRID_CELL: 4_000_000}, device="cuda:0")
cell = build_all("cuda:0", arena=arena)
wp.synchronize_device(dev)

n_str, n_per = int(cell.cortex.n_strands), int(cell.cortex.nodes_per_strand)
lo = int(cell.cortex.nodes.lo)
offsets = lo + np.arange(n_str + 1, dtype=np.int64) * n_per
seg_off = segment_offsets(offsets)
n_actin = lo + n_str * n_per
n_seg = int(seg_off[-1])
print(f"cortex {n_str:,} filaments, {n_seg:,} segments -> "
      f"{n_seg / _PLANE_BLOCKS:.0f} elements per block", flush=True)

pos_d, f_d = arena.node_arrays["position"], arena.node_arrays["force"]
R = float(cell.cortex.radius_um)
CENTRE = (0.0, 0.0, 0.0)
N_PLANES = 64

# ⚠ A zero force field would make every plane sum zero and the comparison vacuous — the ratio would
# be 0/0. Give the nodes a deterministic, non-physical force so both paths have something to sum.
# NOTHING here is a physics claim; this is an arithmetic equivalence check on identical inputs.
rng = np.random.default_rng(20260821)
f_host = np.zeros((pos_d.shape[0], 3), np.float64)
f_host[lo:n_actin] = rng.normal(size=(n_actin - lo, 3)) * 0.5
f_d.assign(wp.array(f_host, dtype=wp.vec3d, device=dev))
wp.synchronize_device(dev)

t = time.perf_counter()
devres = gamma_planes_device(
    pos_d=pos_d, f_ext_d=f_d,
    fiber_offsets_d=wp.array(offsets.astype(np.int32), dtype=wp.int32, device=dev),
    seg_offsets_d=wp.array(seg_off.astype(np.int32), dtype=wp.int32, device=dev),
    n_segments=n_seg, n_actin=n_actin, R_um=R, centre=CENTRE, n_planes=N_PLANES, device=dev)
wp.synchronize_device(dev)
print(f"device path  {1000*(time.perf_counter()-t):8.1f} ms", flush=True)

pos_h = pos_d.numpy()
t = time.perf_counter()
sa, sb, t_actin = actin_axial_tension(f_host, pos_h, offsets)
host_actin = plane_force_sums(pos_h[sa], pos_h[sb], t_actin, R_um=R, centre=CENTRE,
                              n_planes=N_PLANES)
print(f"host path    {1000*(time.perf_counter()-t):8.1f} ms", flush=True)

d, h = np.asarray(devres["actin"], np.float64), np.asarray(host_actin, np.float64)
assert d.shape == h.shape, (d.shape, h.shape)
worst = float(np.max(np.abs(d - h) / (np.abs(h) + 1e-300)))
print(f"n_planes {d.size}  |host| median {np.median(np.abs(h)):.6e}", flush=True)
print(f"WORST RELATIVE ERROR  {worst:.6e}", flush=True)
print(f"tolerance             1.000000e-11   -> {'WITHIN' if worst <= 1e-11 else 'EXCEEDS'}", flush=True)


# ⚠ WHY THIS IS A SCRIPT AND NOT A TEST. It needs the native cortex on a real card, and the production
# environment on the GPU host has no test runner (PI queue item 15). Run it as:
#
#   srun --jobid=<N> --overlap env CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=<n> \
#       python aleph/scripts/world_gamma_native_equiv.py
#
# `srun --jobid --overlap` alone leaves CUDA_VISIBLE_DEVICES unset and all cards visible, and a bare
# process without SLURM_JOB_ID is SIGTERMed by gpu-bypass-watch after ~20 s. Both halves are required.
