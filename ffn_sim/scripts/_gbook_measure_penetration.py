"""Measure the REAL interpenetration depth (µm, R-units) of a saved DCM frames npz — decoupled
from the gate's pen_frac (= max_pen / STALE initial mean_edge, ~30× inflated per N2000_CONFLUENT_FINDINGS).

Reuses the engine's own diagnostic penetration_depth_kernel (Ericson closest-point, node-on-inner-side-of-
neighbour-face sign test) but reports the ABSOLUTE depth distribution + normalization by the CURRENT-frame
mean_edge, so we can tell whether pen_frac≈75 is real overlap or a metric artifact.

Env: NPZ (path), RADIUS_R (query radius in R units, default 1.0), RCELL_UM (default 7.5).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser("~/ff_scratch"))
import warp as wp

wp.init()
from ffn_sim.dcm.dcm_neighbor_warp import penetration_depth_kernel

NPZ = os.environ["NPZ"]
RCELL = float(os.environ.get("RCELL_UM", "7.5")) * 1e-6      # m
RADIUS = float(os.environ.get("RADIUS_R", "1.0")) * RCELL    # query radius (m)
DEV = os.environ.get("DEV", "cuda:0")

d = np.load(NPZ)
frames = d["frames"]
faces = d["faces"].astype(np.int32)
cof = d["cof"].astype(np.int32)
fcell = cof[faces[:, 0]].astype(np.int32)
N = frames.shape[1]
nf = faces.shape[0]
print(f"[pen] {os.path.basename(NPZ)}  N_nodes={N} faces={nf} cells={int(cof.max())+1}  "
      f"radius={RADIUS*1e6:.2f}um", flush=True)

faces_d = wp.array(faces, dtype=wp.int32, device=DEV)
cof_d = wp.array(cof, dtype=wp.int32, device=DEV)
fcell_d = wp.array(fcell, dtype=wp.int32, device=DEV)
grid = wp.HashGrid(128, 128, 128, device=DEV)
pen_d = wp.zeros(N, dtype=wp.float64, device=DEV)


def measure(pos):
    # current-frame mean edge (from faces, deduped by just averaging the 3 sides)
    e = np.concatenate([
        np.linalg.norm(pos[faces[:, 0]] - pos[faces[:, 1]], axis=1),
        np.linalg.norm(pos[faces[:, 1]] - pos[faces[:, 2]], axis=1),
        np.linalg.norm(pos[faces[:, 2]] - pos[faces[:, 0]], axis=1)])
    mean_edge = float(e.mean())
    cent = pos[faces].mean(axis=1).astype(np.float32)
    pos_d = wp.array(pos, dtype=wp.vec3d, device=DEV)
    node_f32 = wp.array(pos.astype(np.float32), dtype=wp.vec3, device=DEV)
    cent_d = wp.array(cent, dtype=wp.vec3, device=DEV)
    grid.build(points=cent_d, radius=RADIUS)
    pen_d.zero_()
    wp.launch(penetration_depth_kernel, dim=N,
              inputs=[grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d, wp.float32(RADIUS), pen_d],
              device=DEV)
    wp.synchronize_device(DEV)
    pen = pen_d.numpy() * 1e6                        # µm (DCM pos in metres)
    return mean_edge * 1e6, pen


RCELL_UM = RCELL * 1e6
ALL = os.environ.get("ALL_FRAMES", "0") == "1"       # 1 → per-frame time series (interpenetration stability)
if ALL:
    steps = d["step"] if "step" in d.files else np.arange(len(frames))
    frame_list = [(f"f{t:02d}s{int(steps[t]):>6}", t) for t in range(len(frames))]
else:
    frame_list = [("INIT ", 0), ("FINAL", len(frames) - 1)]
for label, idx in frame_list:
    me, pen = measure(frames[idx].astype(np.float64))
    inside = pen > 0
    deep = pen > 0.2 * RCELL_UM
    p = pen[inside]
    med = float(np.median(p)) if p.size else 0.0
    mx = float(pen.max())
    # which cells own the deep-penetrating nodes?
    deep_cells = np.unique(cof[deep])
    print(f"[{label}] mean_edge={me:.2f}um  n_inside={int(inside.sum())} ({inside.mean()*100:.3f}%)  "
          f"n_deep>0.2R={int(deep.sum())} in {deep_cells.size} cells  "
          f"max_pen={mx:.2f}um ({mx/RCELL_UM:.2f}R)  median_inside={med:.2f}um ({med/RCELL_UM:.2f}R)  "
          f"|| gate max/init_mean_edge would be ~{mx/me:.1f} (vs gate's stale ~75)", flush=True)
print("PEN DONE", flush=True)
