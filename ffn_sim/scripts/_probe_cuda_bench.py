"""Micro-benchmark the GPU filopodia probe on cuda — isolate grid-build vs kernel cost,
and the effect of the query radius (tip_capture+R vs tip_capture+mean_edge)."""
import sys, time
import numpy as np
import warp as wp
from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
from ffn_sim.warp_port.dcm_filopodia_probe_warp import probe_faces_gpu, filopodia_probe_kernel
from ffn_sim.warp_port.dcm_neighbor_warp import pos_to_f32, face_centroids_f32

dev = sys.argv[1] if len(sys.argv) > 1 else "cuda:0"
wp.init()
p = ResolvedDCM(subdivisions=2); R = p.R_cell
v1, e1, f1 = icosphere_mesh(R, 2); npc = v1.shape[0]
me = float(np.linalg.norm(v1[e1[:, 0]] - v1[e1[:, 1]], axis=1).mean())

# build an N-cell fcc-ish cluster
N = 100
import itertools
side = int(np.ceil(N ** (1 / 3)))
centers = np.array(list(itertools.product(range(side), repeat=3)), float)[:N] * (1.8 * R)
verts = np.concatenate([v1 + c for c in centers])
faces = np.concatenate([f1 + k * npc for k in range(N)]).astype(np.int64)
fcell = np.repeat(np.arange(N), f1.shape[0]).astype(np.int64)
nf = faces.shape[0]
print(f"N={N} cells, {verts.shape[0]} nodes, {nf} faces, mean_edge={me*1e6:.2f}um, R={R*1e6:.2f}um")

# free tips: ~3000 near random surface nodes
rng = np.random.default_rng(0)
ntip = 3000
base = rng.integers(0, verts.shape[0], ntip)
tips = verts[base] + rng.normal(0, 0.3e-6, (ntip, 3))
owncell = (base // npc).astype(np.int32)
tip_capture = 0.8e-6

def timeit(radius_R, label, reps=3):
    # patch the radius by passing R or a small value through the R arg
    t = []
    for _ in range(reps):
        wp.synchronize_device(dev)
        t0 = time.perf_counter()
        bf, bb, bd = probe_faces_gpu(free_tips_xyz=tips, free_tip_owncell=owncell, pos=verts,
                                     faces=faces, fcell=fcell, tip_capture=tip_capture,
                                     R=radius_R, device=dev)
        wp.synchronize_device(dev)
        t.append(time.perf_counter() - t0)
    print(f"  {label:28} radius={tip_capture+radius_R+tip_capture:.1e}m  "
          f"time={min(t)*1e3:.1f}ms  adhered={(bf>=0).sum()}")

print("=== probe_faces_gpu timing (min of 3) ===")
timeit(R, "tip_capture+R (current 9.1um)")
timeit(2.0 * me, "tip_capture+2*mean_edge")
timeit(me, "tip_capture+mean_edge")
