"""Animated 3D viewer — the FF cortex DEFORMING under parallel-plate compression, with the
per-filament von-Mises STRESS field colour-mapped as it deforms.

Runs the compression driver at a sequence of strains (rigid plate) and, at each, captures both
the deformed cortex geometry AND its per-node von-Mises stress (ff_virial_stress). Two scenes:
a plain-filament deformation and a turbo-coloured stress field, both animated (press ▶). Full
resolution (no downsampling). Self-contained double-clickable HTML (three.js), not a CSP artifact.

Usage: python -m aleph.scripts.ff_cortex_deform_viewer --nfil 70686 --device cuda:0
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy.spatial import ConvexHull

from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device
from aleph.laws.ff_virial_stress import cortex_node_stress, cortex_node_areal_strain, face_areas
from aleph.scripts.ff_s1_sphere import LOAD_PHYSIO  # rigid plate (correct native measurement)
from aleph.scripts.ff_resting_viewer import _turbo, _turbo_gradient
from aleph.scripts.ff_viewer_html import build_viewer


def compression_frames(nfil, strains, n_steps, device, seed=1):
    """Return (geom_frames, svm_frames, seg, R0, forces): deformed cortex + per-node stress per strain."""
    geom, svm_frames, strain_frames, forces = [], [], [], []
    seg = faces = area0 = None
    R0 = float(np.linalg.norm(
        (b := build_crosslinked_cortex(CortexParams(), n_filaments=nfil, n_xl=nfil,
                                       n_myo=max(1, nfil // 10), rng=np.random.default_rng(seed))).net.pos
        - b.net.pos.mean(0), axis=1).mean())
    V0 = (4.0 / 3.0) * np.pi * R0 ** 3
    for s in strains:
        cx = build_crosslinked_cortex(CortexParams(), n_filaments=nfil, n_xl=nfil,
                                      n_myo=max(1, nfil // 10), rng=np.random.default_rng(seed))
        cx.R0_mean = R0
        pos_all, m = simulate_whole_cell_compression_on_device(
            cx, NMIIA_MINIFIL_STALL_PN, strain=float(s), nucleus=None, membrane=None,
            microtubule=None, n_steps=n_steps,
            turgor_every=(50 if device.startswith("cuda") else 20), device=device, **LOAD_PHYSIO)
        Nc = cx.net.n_nodes
        pcx = np.asarray(pos_all)[:Nc]
        if seg is None:
            seg = cx.net.segments
            faces = ConvexHull(pcx).simplices.astype(np.int64)   # fixed triangulation from the rest frame
            area0 = face_areas(pcx, faces)
        # remap sim-z (plate-compression axis) → viewer-y (up) so the flattening is visible.
        geom.append(pcx[seg].reshape(-1, 3)[:, [0, 2, 1]])
        strain_frames.append(cortex_node_areal_strain(pcx, faces, area0))
        svm, _ = cortex_node_stress(
            pcx, np.full(Nc, V0 / Nc),
            xl_ij=np.stack([cx.xl_i, cx.xl_j], 1).astype(np.int64), k_xl=cx.xl_k, r0_xl=cx.xl_rest,
            myo_ij=np.stack([cx.myo_i, cx.myo_j], 1).astype(np.int64),
            f_myo=float(NMIIA_MINIFIL_STALL_PN), dP=float(m["dP_turgor_Pa"]))
        svm_frames.append(svm)
        forces.append(float(m["F_plate_pN"]))
        print(f"strain {s*100:4.1f}%  F={m['F_plate_pN']:9.0f} pN  svm_p95={np.percentile(svm,95):.0f} Pa", flush=True)
    return geom, svm_frames, strain_frames, seg, R0, forces


def main() -> None:
    ap = argparse.ArgumentParser(description="FF cortex deformation + stress-field 3D viewer.")
    ap.add_argument("--nfil", type=int, default=70686)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--strains", default="0.0,0.05,0.10,0.15,0.20,0.25,0.30")
    ap.add_argument("--out", default="aleph/outputs/mech_hier/figs/cortex_deform_3d.html")
    a = ap.parse_args()

    strains = [float(x) for x in a.strains.split(",")]
    geom, svm_frames, strain_frames, seg, R0, forces = compression_frames(a.nfil, strains, a.steps, a.device)
    nseg = geom[0].shape[0] // 2
    tag = f"NF={a.nfil}, {nseg} segments"
    smax = int(strains[-1] * 100)

    plain = {"name": f"cortex ({tag})", "kind": "lines", "verts": geom[0], "frames": geom,
             "color": "#8fbff0", "size": 1.3, "opacity": 0.5}
    hi = float(max(np.percentile(np.concatenate(svm_frames), 98), 1e-6))
    cframes = [_turbo(svm[seg].reshape(-1), 0.0, hi) for svm in svm_frames]
    stress = {"name": f"von-Mises stress ({tag})", "kind": "lines", "verts": geom[0], "frames": geom,
              "color": "#8fbff0", "size": 1.6, "opacity": 0.75, "color_frames": cframes, "cbar": "σ_vm [Pa]"}
    smag = float(max(np.percentile(np.abs(np.concatenate(strain_frames)), 98), 1e-6))
    sframes = [_turbo(st[seg].reshape(-1), -smag, smag) for st in strain_frames]
    strain = {"name": f"areal strain ({tag})", "kind": "lines", "verts": geom[0], "frames": geom,
              "color": "#8fbff0", "size": 1.6, "opacity": 0.75, "color_frames": sframes, "cbar": "areal strain"}

    s_plain = f"plate compression 0→{smax}% (plain, press ▶)"
    s_stress = f"von-Mises STRESS 0→{smax}% (press ▶)"
    s_strain = f"areal STRAIN 0→{smax}% (press ▶)"
    scenes = {s_stress: [stress], s_strain: [strain], s_plain: [plain]}
    cbars = {s_stress: {"lo": 0.0, "hi": hi, "label": "cortex σ_vm", "unit": "Pa", "grad": _turbo_gradient()},
             s_strain: {"lo": -smag, "hi": smag, "label": "areal strain", "unit": "", "grad": _turbo_gradient()}}
    build_viewer(scenes, out=a.out, cbars=cbars,
                 title=f"FF cortex under compression — deformation + stress (NF={a.nfil}, R0={R0:.2f}µm, "
                       f"F {forces[0]:.0f}→{forces[-1]:.0f} pN)")
    print(f"wrote {a.out}  (frames {[f'{s*100:.0f}%' for s in strains]}, σ_vm hi={hi:.0f} Pa)", flush=True)


if __name__ == "__main__":
    main()
