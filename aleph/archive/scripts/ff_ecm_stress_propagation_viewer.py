"""FF ECM — interactive HTML viewer for stress propagation from a contractile inclusion (ROADMAP MID, KB-1.10).

The morphology counterpart to `figs/stress_propagation_sp_native.png`: SHOWS the displacement field a contractile
inclusion drives through collagen-I, ISOTROPIC vs ALIGNED, so the reader can SEE the directional channeling that
makes the aligned exponent smaller (stress reaches farther along the fiber director). Each scene renders the
(deformed) fiber network faint, the frozen contractile inclusion (red), and every matrix node coloured by its
displacement magnitude in discrete bins (blue = still → red = strongly pulled). For the aligned scene the pulled
region is elongated along the green director; for the isotropic scene it is round.

A representative native-DENSITY REV (default 40 µm) — the mechanism view; the quantitative n(S) exponent is the
native box=90 µm result in `ff_ecm_stress_propagation.py`. Output (standard ecm_lib location):
aleph/outputs/ff/ecm_lib/stress_propagation_viewer.html

Run:  python -m aleph.scripts.ff_ecm_stress_propagation_viewer            # iso vs aligned, ~40 µm REV
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import warp as wp

from aleph.laws import ecm_library as L
from aleph.scripts.ff_ecm_stress_propagation import _relax_inclusion
from aleph.scripts.ff_ecm_viewer import ecm_layers
from aleph.scripts.ff_contact_guidance_viewer import director_line
from aleph.scripts.ff_viewer_html import build_viewer

OUT = "aleph/outputs/ff/ecm_lib"
# blue → cyan → green → orange → red, by displacement bin
_BINS = ["#2b5fff", "#00c2c2", "#39c24a", "#ff9f1c", "#ff3333"]


def _disp_layers(pos, pos0):
    """Displacement-magnitude layers that ISOLATE the pulled region: the still bulk (below the 80th displacement
    percentile) is a faint gray backdrop; the moving top-20% is split into blue→red bins by percentile, so the
    SHAPE of the pulled region reads (round for isotropic, elongated along the director for aligned)."""
    disp = np.linalg.norm(pos - pos0, axis=1)
    thr = float(np.percentile(disp, 80))                          # below this = "still" backdrop
    still = disp < thr
    layers = [{"name": "still matrix (<80th pct)", "kind": "points", "verts": pos[still],
               "color": "#334155", "size": 1.0, "opacity": 0.10}]
    move = disp >= thr
    if move.sum() == 0:
        return layers, disp
    dm = disp[move]
    edges = np.percentile(dm, np.linspace(0, 100, len(_BINS) + 1))
    for b in range(len(_BINS)):
        sel = move & (disp >= edges[b]) & ((disp < edges[b + 1]) if b < len(_BINS) - 1 else (disp <= edges[b + 1]))
        if sel.sum() == 0:
            continue
        layers.append({"name": f"|Δ|∈[{edges[b]:.2f},{edges[b+1]:.2f})µm", "kind": "points",
                       "verts": pos[sel], "color": _BINS[b], "size": 2.4 + 0.9 * b, "on_top": b >= len(_BINS) - 2})
    return layers, disp


def build(box=40.0, conc=1.5, R_incl=None, eps=0.2, steps=6000, S_list=(0.0, 0.83),
          device="cpu", out=None):
    os.makedirs(OUT, exist_ok=True)
    out = out or f"{OUT}/stress_propagation_viewer.html"
    R_incl = R_incl if R_incl is not None else box / 6.0
    spec = L.get_spec("collagen_I")
    lo, hi = [0.0, 0.0, 0.0], [box, box, box]
    center = 0.5 * (np.asarray(lo) + np.asarray(hi))
    wp.init()
    scenes = {}
    for S in S_list:
        rng = np.random.default_rng(7)
        ecm = L.build_fibrillar_ecm(spec, lo, hi, concentration=conc, dim=3, alignment_S=float(S),
                                    director=(1.0, 0.0, 0.0), pin_faces=(), target_z=3.2, rng=rng)
        pos0 = np.ascontiguousarray(ecm.net.pos, np.float64).copy()
        pos, meta = _relax_inclusion(ecm, center, R_incl, eps, steps=steps, device=device)
        dlayers, disp = _disp_layers(pos, pos0)
        # deformed fibers (faint) — rebuild layers from the deformed positions
        ecm.net.pos[:] = pos
        flayers = [L_ for L_ in ecm_layers(ecm, name=f"S={S:g}") if L_["name"] in ("fibers", "crosslinks")]
        for L_ in flayers:
            L_["opacity"] = 0.12                                   # faint so the displacement points read
        incl_mask = np.linalg.norm(pos0 - center, axis=1) < R_incl
        incl_layer = {"name": "contractile inclusion", "kind": "points", "verts": pos[incl_mask],
                      "color": "#ff2266", "size": 3.0, "on_top": True}
        tag = (f"collagen-I  S={S:g} ({'ISOTROPIC' if S < 0.1 else 'ALIGNED'})  ε={eps:.0%} inclusion  "
               f"max|Δ|={disp.max():.2f}µm")
        scenes[tag] = flayers + [incl_layer] + dlayers + [director_line(lo, hi, (1.0, 0.0, 0.0))]
    title = ("FF ECM MID — stress propagation: displacement field a contractile inclusion drives through "
             "collagen-I (ISOTROPIC vs ALIGNED; blue=still→red=pulled; green=director) — the morphology behind "
             "the r⁻ⁿ channeling")
    build_viewer(scenes, out, title=title)
    print(f"wrote {out}  ({len(scenes)} scenes)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--box", type=float, default=40.0)
    ap.add_argument("--conc", type=float, default=1.5)
    ap.add_argument("--eps", type=float, default=0.2)
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--S", default="0.0,0.83")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    S_list = tuple(float(x) for x in a.S.split(",") if x.strip())
    build(box=a.box, conc=a.conc, eps=a.eps, steps=a.steps, S_list=S_list, device=a.device, out=a.out)


if __name__ == "__main__":
    main()
