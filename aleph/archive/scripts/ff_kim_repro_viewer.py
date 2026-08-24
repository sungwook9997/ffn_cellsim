"""Interactive 3D HTML morphology viewer for the FF ⟷ Slater/Kim 2021 reproduction — matching Kim Fig 1B/C:
the contracted cortex ring + collagen fibers colored by TENSILE force (the radial tension chains that carry
the cell's contraction to the boundary) + crosslinks, in the emergent remodeled state.

Reads the geometry npz saved by ``ff_kim_repro --save-geom`` and reuses ``ff_viewer_html.build_viewer`` with a
per-vertex diverging tension colormap (blue=compression → red=tension) + colorbar legend. Two scenes: the
relaxed initial matrix (Kim Fig 1B) and the contracted, tension-colored state (Kim Fig 1C).

Run:  python -m aleph.scripts.ff_kim_repro_viewer [--npz PATH]
Out:  aleph/outputs/ff/kim_repro/figs/kim_morphology.html
"""

from __future__ import annotations

import argparse

import numpy as np

from aleph.scripts.ff_viewer_html import build_viewer

OUT = "aleph/outputs/ff/kim_repro"

# Diverging tension bands (Kim Fig 1C colour scaling) → one single-colour ``lines`` layer per band (robust;
# avoids the per-vertex color_frames path, which needs a frames array and breaks for a static field).
_BANDS = [(-1e9, -0.5, "#2166ac", "compression (strong)"), (-0.5, -0.15, "#67a9cf", "compression"),
          (-0.15, 0.15, "#999999", "~neutral"), (0.15, 0.5, "#ef8a62", "tension"),
          (0.5, 1e9, "#b2182b", "tension (strong)")]


def _tension_band_layers(seg_pairs, seg_tension):
    """Split fibers into diverging tension bands (normalised by the 90th-pct |tension|), one layer each."""
    t = np.asarray(seg_tension, float)
    scale = float(np.percentile(np.abs(t), 90)) or 1.0
    tn = t / scale
    layers = []
    for lo, hi, col, name in _BANDS:
        m = (tn >= lo) & (tn < hi)
        if m.any():
            is_tension = "tension" in name
            layers.append({"name": f"{name}", "kind": "lines", "verts": seg_pairs[m], "color": col,
                           "opacity": 0.35 if "neutral" in name else (0.9 if is_tension else 0.6),
                           "size": 2.0 if "strong" in name else (1.5 if is_tension else 1.0),
                           "on_top": is_tension})               # draw the radial tension chains over the compressed bulk
    return layers, scale


def build(npz, out=None):
    d = np.load(npz)
    mi, mf = d["matrix_init"], d["matrix_final"]
    seg_i, seg_j = d["seg_i"], d["seg_j"]
    seg_f = np.stack([mf[seg_i], mf[seg_j]], axis=1)                 # (S,2,3) contracted
    seg_0 = np.stack([mi[seg_i], mi[seg_j]], axis=1)                 # (S,2,3) relaxed

    band_layers, scale = _tension_band_layers(seg_f, d["seg_tension"])
    contracted = band_layers + [
        {"name": f"cortex (contracted, tension ±{scale:.0f} pN scale)", "kind": "points",
         "verts": d["cortex_final"], "color": "#39d353", "size": 6.0, "on_top": True}]
    if d["xl_i"].shape[0]:
        xl = np.stack([mf[d["xl_i"]], mf[d["xl_j"]]], axis=1)
        contracted.append({"name": "crosslinks", "kind": "lines", "verts": xl,
                           "color": "#ffffff", "opacity": 0.1, "size": 0.6})
    initial = [{"name": "fibers (relaxed)", "kind": "lines", "verts": seg_0,
                "color": "#e8a33d", "opacity": 0.4, "size": 1.0}]

    scenes = {"Contracted — fibers colored by tension (≈ Kim Fig 1C)": contracted,
              "Initial — relaxed collagen matrix (≈ Kim Fig 1B)": initial}
    out = out or f"{OUT}/figs/kim_morphology.html"
    build_viewer(scenes, out,
                 title="FF ⟷ Slater/Kim 2021 — a contractile cell remodels a viscoelastic collagen matrix (tension field)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=f"{OUT}/kim_geom_geom.npz")
    a = ap.parse_args()
    print("wrote", build(a.npz), flush=True)


if __name__ == "__main__":
    main()
