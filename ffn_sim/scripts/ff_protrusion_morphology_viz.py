"""Render a ff_protrusion_into_ecm result as CELL MORPHOLOGY — membrane surface + filopodia FINGERS + nucleus
+ ECM — instead of bare actin centerlines (PI 2026-07-02 "이건 필로포디움이 아니라 그냥 작대기").

Per frame it reconstructs each filopodium's base→tip from the grown actin bundle and wraps it in a membrane
finger-tube (ff_cell_morphology.finger_tube), so the fingers GROW into the ECM as proper 3D fingers. Membrane +
nucleus are surfaces; the ECM Mikado is grey lines; the actin cores are faint lines inside the fingers.
"""
from __future__ import annotations

import argparse

import numpy as np

from ffn_sim.scripts.ff_viewer_html import build_viewer
from ffn_sim.scripts.ff_cell_morphology import uv_sphere, finger_tube, merge


def _finger_meshes(filo_nodes, finger_of_fiber, nb_filo, c, R, radius=0.18):
    """From this frame's filo nodes (F,nb,3) + finger grouping → one merged finger-tube mesh (verts, faces)."""
    F = filo_nodes.shape[0]
    meshes = []
    for fid in range(int(finger_of_fiber.max()) + 1):
        fibs = filo_nodes[finger_of_fiber == fid]                 # (nf, nb, 3)
        if fibs.shape[0] == 0:
            continue
        tip = fibs[:, -1, :].mean(0)                              # grown barbed tips
        graft = fibs[:, 0, :].mean(0)                             # bundle base (inside cortex)
        axis = tip - graft; L = np.linalg.norm(axis)
        if L < 1e-3:
            continue
        axis /= L
        base_surf = c + R * axis                                  # finger emerges at the membrane surface
        length = max(float(np.linalg.norm(tip - base_surf)), 0.2)
        meshes.append(finger_tube(base_surf, axis, length, radius))
    return merge(meshes) if meshes else (np.zeros((0, 3), np.float32), np.zeros((0, 3), np.uint32))


def _ecm_lines(ecm_nodes, n_ecm_fib, nb_ecm):
    """ECM node array → line-segment endpoint pairs (S,2,3)."""
    p = ecm_nodes.reshape(n_ecm_fib, nb_ecm, 3)
    segs = np.stack([p[:, :-1, :], p[:, 1:, :]], 2).reshape(-1, 2, 3)
    return segs.astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="ffn_sim/outputs/ff/figs/protrusion_ecm_native7.npz")
    ap.add_argument("--out", default="ffn_sim/outputs/ff/figs/protrusion_ecm_morphology")
    ap.add_argument("--r-nuc-frac", type=float, default=0.25)     # sim value; MCF7 nucleus is larger (note)
    args = ap.parse_args()
    d = np.load(args.npz, allow_pickle=True)
    ncs = int(d["n_cortex_sub"]); nff = int(d["n_filo_fib"]); nbf = int(d["nb_filo"])
    nef = int(d["n_ecm_fib"]); nbe = int(d["nb_ecm"])
    fof = d["finger_of_fiber"] if "finger_of_fiber" in d else np.zeros(nff, np.int32)
    R = float(d["R"]); c = np.asarray(d["c"], float)
    fe = d["frames_ecm"]                                          # (macro, Nviz, 3)
    filo0, ecm0, nuc0 = ncs, ncs + nff * nbf, ncs + nff * nbf + nef * nbe

    # static compartments: membrane + nucleus spheres
    mem_v, mem_f = uv_sphere(c, R)
    nuc_v, nuc_f = uv_sphere(c, args.r_nuc_frac * R, nu=22, nv=14)
    # ECM lines from the final frame
    ecm_lines = _ecm_lines(fe[-1][ecm0:nuc0], nef, nbe)

    # animated finger meshes + actin cores, per frame
    finger_frames, actin_frames = [], []
    f0_verts = None
    for fr in fe:
        filo = fr[filo0:ecm0].reshape(nff, nbf, 3)
        fv, ff = _finger_meshes(filo, fof, nbf, c, R)
        if f0_verts is None:
            f0_verts, f0_faces = fv, ff
        # pad/truncate to constant vertex count (fingers can vanish if L<thresh early; keep count stable)
        if fv.shape[0] != f0_verts.shape[0]:
            fv2 = np.zeros_like(f0_verts); fv2[:min(len(fv), len(f0_verts))] = fv[:len(f0_verts)]; fv = fv2
        finger_frames.append(fv)
        # actin cores as line segments
        segs = np.stack([filo[:, :-1, :], filo[:, 1:, :]], 2).reshape(-1, 2, 3).astype(np.float32)
        actin_frames.append(segs)

    scenes = {"MCF7 cell — filopodia into ECM (morphology)": [
        {"name": "plasma membrane", "kind": "mesh", "color": "#7db8e8", "opacity": 0.22,
         "verts": mem_v, "faces": mem_f},
        {"name": f"nucleus ({args.r_nuc_frac:.2f}R)", "kind": "mesh", "color": "#c65b7c", "opacity": 0.95,
         "verts": nuc_v, "faces": nuc_f},
        {"name": f"filopodia ({int(fof.max())+1} fingers, growing)", "kind": "mesh", "color": "#ffb000",
         "opacity": 1.0, "verts": f0_verts, "faces": f0_faces, "frames": finger_frames},
        {"name": "actin bundle cores", "kind": "lines", "color": "#e07a00",
         "verts": actin_frames[0], "frames": actin_frames},
        {"name": "collagen ECM (Mikado)", "kind": "lines", "color": "#8a8a8a", "verts": ecm_lines},
    ]}
    build_viewer(scenes, out=f"{args.out}.html",
                 title="FF MCF7 cell — filopodia protruding into a collagen ECM (membrane + fingers + nucleus)")
    print("wrote", f"{args.out}.html", f"({int(fof.max())+1} fingers, {len(fe)} frames)")


if __name__ == "__main__":
    main()
