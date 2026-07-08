"""Resting full-compartment cell viewer — the persistent interactive record of the validated FF cortical-mechanics
state (PI viz rule: the actual 3D cell MORPHOLOGY, explorable, compartments visible; NO fiber downsampling).

Renders, from ff_resting_full_compartment.py's npz: the woven actin CORTEX (all filaments, full-res), the plasma
MEMBRANE envelope (translucent), and the NUCLEUS — with a cut plane (clip) so the interior compartments show, and a
σ_vm stress scene + turbo colorbar. The cytoplasm is the pressure-borne turgor interior (ΔP annotated in the title).

Usage: python ff_resting_viewer.py <resting_cell.npz> <out.html>
"""
import sys
import numpy as np
import matplotlib.cm as cm
from scipy.spatial import ConvexHull

from ffn_sim.scripts.ff_viewer_html import build_viewer


def _turbo(scalar, lo, hi):
    rng = max(float(hi) - float(lo), 1e-12)
    return (cm.get_cmap("turbo")(np.clip((np.asarray(scalar) - lo) / rng, 0.0, 1.0))[:, :3] * 255).astype(np.uint8)


def _turbo_gradient(n=16):
    t = cm.get_cmap("turbo")
    return ["#%02x%02x%02x" % tuple((np.asarray(t(x)[:3]) * 255).astype(int)) for x in np.linspace(0, 1, n)]


def build(npz_path: str, out: str, title: str | None = None) -> str:
    d = np.load(npz_path)
    pos = np.asarray(d["frames"], np.float32)[0]                 # (N,3) resting frame
    Nc, Ne, n_nuc = int(d["Nc"]), int(d["Ne"]), int(d["n_nuc"])
    foff = np.asarray(d["foff"], np.int64)
    faces = np.asarray(d["faces"], np.int64)
    R, R_nuc = float(d["R"]), float(d["R_nuc"])
    dP = float(d["dP"]) if "dP" in d.files else 0.0
    gamma = float(d["gamma_mN_m"]) if "gamma_mN_m" in d.files else 0.0
    pcx = pos[:Nc]
    nuc = pos[Ne:] if n_nuc > 0 else None

    # cortex filament segments — ALL fibers, NO downsampling (PI rule)
    cort_fibers = [f for f in range(len(foff) - 1) if int(foff[f]) < Nc]
    seg = np.array([(n, n + 1) for f in cort_fibers
                    for n in range(int(foff[f]), int(foff[f + 1]) - 1)], dtype=np.int64)
    fil_flat = pos[seg].reshape(-1, 3)                           # (2S,3) segment endpoints
    fib_name = f"cortex filaments (woven actin — {len(cort_fibers)} fibers, {seg.shape[0]} segments, FULL-RES)"

    membrane = {"name": "plasma membrane (reservoir envelope)", "kind": "mesh", "verts": pcx, "faces": faces,
                "color": "#7db8e8", "opacity": 0.12, "clip": True}

    def nucleus_layer():
        nfaces = ConvexHull(nuc).simplices.astype(np.int64)
        return {"name": f"nucleus ({n_nuc} beads, R_nuc=0.70R={R_nuc:.2f}µm)", "kind": "mesh", "verts": nuc,
                "faces": nfaces, "color": "#d17fe0", "opacity": 0.97}

    # pre-cut "cut-away" — drop the front cortex hemisphere (segment midpoint on the +x side of the centroid) so
    # the nucleus + interior compartments are visible in the DEFAULT view (native cortex is opaque from outside).
    cx = float(pcx[:, 0].mean())
    seg_mid_x = 0.5 * (pos[seg[:, 0], 0] + pos[seg[:, 1], 0])
    back = seg_mid_x <= cx
    fil_back = pos[seg[back]].reshape(-1, 3)
    cutaway = [{"name": f"cortex (back hemisphere cut-away — {int(back.sum())} of {seg.shape[0]} segments)",
                "kind": "lines", "verts": fil_back, "color": "#8fbff0", "size": 1.3, "opacity": 0.7}]
    if nuc is not None:
        cutaway.append(nucleus_layer())
    cutaway.append({"name": "plasma membrane (envelope)", "kind": "mesh", "verts": pcx, "faces": faces,
                    "color": "#7db8e8", "opacity": 0.06})

    morph = [{"name": fib_name, "kind": "lines", "verts": fil_flat, "color": "#8fbff0", "size": 1.2, "opacity": 0.45, "clip": True},
             membrane]
    if nuc is not None:
        morph.append(nucleus_layer())
    scenes = {"cut-away (compartments visible)": cutaway, "full morphology (all filaments)": morph}
    cbars = {}

    svm = np.asarray(d["svm"], np.float32)[0] if "svm" in d.files else None
    if svm is not None:
        lo, hi = 0.0, float(max(np.percentile(svm, 95), 1e-6))
        cf = _turbo(svm[seg].reshape(-1), lo, hi)               # per-segment-endpoint colour
        sname = "cortex von-Mises stress"
        svm_layers = [{"name": fib_name, "kind": "lines", "verts": fil_flat, "color": "#8fbff0", "size": 1.4,
                       "opacity": 0.55, "clip": True, "color_frames": [cf]}, dict(membrane)]
        if nuc is not None:
            svm_layers.append(nucleus_layer())
        scenes[sname] = svm_layers
        cbars[sname] = {"lo": lo, "hi": hi, "label": "cortex σ_vm [Pa]", "stops": _turbo_gradient()}

    ttl = title or (f"FF MCF7 cell — VALIDATED resting full-compartment · R={R:.1f}µm · "
                    f"pressure-borne turgor ΔP={dP:.0f} Pa → γ={gamma:.2f} mN/m · cortex+membrane+nucleus+cytoplasm")
    return build_viewer(scenes, out=out, title=ttl, cbars=cbars or None)


if __name__ == "__main__":
    npz = sys.argv[1] if len(sys.argv) > 1 else "resting_cell.npz"
    out = sys.argv[2] if len(sys.argv) > 2 else "ffn_sim/outputs/ff/figs/resting_full_compartment.html"
    print(build(npz, out))
