"""Standalone interactive HTML viewer for FF ECM networks (fibers + crosslinks + BC + indenter) — no cell.

The existing FF viewers (``ff_crawl_viewer``, ``ff_protrusion_ecm_viz``) render the ECM only as an overlay
on a cell. This renders the ECM matrix ON ITS OWN so each pre-built environment (material × alignment × 2D/3D)
can be visually inspected — the PI's "여러 상태로 미리 구현해두면 … 미리 상태 파악". Reuses ``ff_viewer_html.build_viewer``:
fibers as a ``lines`` layer (from ``net.segments``), crosslinks as a second ``lines`` layer, pinned nodes as
``points``, and an optional spherical indenter as a ``points`` shell. Many states go into ONE HTML with a
scene dropdown, so a whole gallery is one double-clickable file.

Run:  python -m aleph.scripts.ff_ecm_viewer            # build the default gallery
Out:  aleph/outputs/ff/ecm_lib/figs/ecm_gallery.html
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from aleph.laws import ecm_library as L
from aleph.scripts.ff_viewer_html import build_viewer

OUT = "aleph/outputs/ff/ecm_lib/figs"

_FIBER_COLORS = {"collagen_I": "#e8a33d", "fibrin": "#d94f4f", "agarose": "#7bc86c"}
_GEL_COLOR = "#5aa9e6"


def ecm_layers(ecm, *, name="ecm", indenter=None) -> list:
    """Layers for one ECM network: fibers (or gel bonds) + crosslinks + pinned nodes [+ indenter shell]."""
    pos = ecm.net.pos
    layers = []
    is_gel = ecm.net.bend_triples.shape[0] == 0
    # fiber / lattice segments as line pairs
    if is_gel and ecm.seg_i.size:
        seg = np.stack([pos[ecm.seg_i], pos[ecm.seg_j]], axis=1)      # (M,2,3)
        layers.append({"name": "gel bonds", "kind": "lines", "verts": seg, "color": _GEL_COLOR,
                       "opacity": 0.25, "size": 1.0})
    elif ecm.net.segments.shape[0]:
        seg = pos[ecm.net.segments]                                    # (S,2,3)
        col = _FIBER_COLORS.get(ecm.material.split("+")[0], "#e8a33d")
        layers.append({"name": "fibers", "kind": "lines", "verts": seg, "color": col,
                       "opacity": 0.55, "size": 1.2})
    # crosslinks
    if ecm.xl_i.size:
        xl = np.stack([pos[ecm.xl_i], pos[ecm.xl_j]], axis=1)
        layers.append({"name": "crosslinks", "kind": "lines", "verts": xl, "color": "#ffffff",
                       "opacity": 0.5, "size": 2.0, "on_top": True})
    # pinned BC nodes
    if ecm.pinned.any():
        layers.append({"name": "pinned BC", "kind": "points", "verts": pos[ecm.pinned],
                       "color": "#3355ff", "size": 3.0})
    # indenter shell (visual)
    if indenter is not None:
        c, R = indenter
        u = np.random.default_rng(0).normal(size=(1200, 3))
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        layers.append({"name": "indenter", "kind": "points", "verts": c + R * u,
                       "color": "#ff3333", "size": 2.0, "on_top": True})
    return layers


def layers_from_npz(path, name):
    """Build viewer layers from a saved native ECM npz (pos/segments/xl/pinned)."""
    d = np.load(path, allow_pickle=True)
    pos = d["pos"].astype(np.float32)
    layers = []
    seg = d["segments"]
    if seg.shape[0]:
        col = _FIBER_COLORS.get(str(d["material"]).split("+")[0], "#e8a33d")
        layers.append({"name": "fibers", "kind": "lines", "verts": pos[seg], "color": col,
                       "opacity": 0.5, "size": 1.0})
    if d["xl_i"].shape[0]:
        xl = np.stack([pos[d["xl_i"]], pos[d["xl_j"]]], axis=1)
        layers.append({"name": "crosslinks", "kind": "lines", "verts": xl, "color": "#ffffff",
                       "opacity": 0.4, "size": 1.5, "on_top": True})
    if d["pinned"].any():
        layers.append({"name": "pinned BC", "kind": "points", "verts": pos[d["pinned"]],
                       "color": "#3355ff", "size": 2.0})
    return layers


def native_viewer(npz_dir="aleph/outputs/ff/ecm_lib/native", out=None):
    """Render the native full-extent (5R×5R) ECM builds into one HTML."""
    out = out or f"{OUT}/ecm_native_gallery.html"
    scenes = {}
    for tag, fn in [("collagen 5R×5R isotropic (native 94.5k nodes)", "collagen_iso_native.npz"),
                    ("collagen 5R×5R aligned S=0.6 (native)", "collagen_aln_native.npz")]:
        p = os.path.join(npz_dir, fn)
        if os.path.exists(p):
            scenes[tag] = layers_from_npz(p, tag)
    if scenes:
        build_viewer(scenes, out, title="FF ECM library — native 5R×5R full extent")
        print(f"wrote {out}  ({len(scenes)} scenes)")
    return out


def default_gallery(out=None, box=22.0, device="cpu"):
    """Build a gallery of representative ECM states into ONE multi-scene HTML."""
    os.makedirs(OUT, exist_ok=True)
    out = out or f"{OUT}/ecm_gallery.html"
    lo, hi = [0, 0, 0], [box, box, box]
    scenes = {}
    rng = lambda s: np.random.default_rng(s)
    # collagen: isotropic / TACS-3 aligned / tendon-aligned / 2D sheet
    scenes["collagen 3D isotropic (dermis)"] = ecm_layers(
        L.build_fibrillar_ecm(L.get_spec("collagen_I"), lo, hi, concentration=2.0, dim=3,
                              alignment_S=0.0, pin_faces=("z_lo",), rng=rng(1)))
    scenes["collagen 3D aligned S=0.6 (TACS-3)"] = ecm_layers(
        L.build_fibrillar_ecm(L.get_spec("collagen_I"), lo, hi, concentration=2.0, dim=3,
                              alignment_S=0.6, director=(1, 0, 0), pin_faces=("z_lo",), rng=rng(2)))
    scenes["collagen 3D aligned S=0.85 (tendon)"] = ecm_layers(
        L.build_fibrillar_ecm(L.get_spec("collagen_I"), lo, hi, concentration=2.0, dim=3,
                              alignment_S=0.85, director=(1, 0, 0), pin_faces=("z_lo",), rng=rng(3)))
    scenes["collagen 2D sheet (aligned S=0.5)"] = ecm_layers(
        L.build_fibrillar_ecm(L.get_spec("collagen_I"), [0, 0, 0], [box, box, 3.0], concentration=2.0,
                              dim=2, alignment_S=0.5, director=(1, 0, 0), pin_faces=(), rng=rng(4)))
    scenes["fibrin 3D isotropic (clot)"] = ecm_layers(
        L.build_fibrillar_ecm(L.get_spec("fibrin"), lo, hi, concentration=2.0, dim=3,
                              alignment_S=0.0, pin_faces=("z_lo",), rng=rng(5)))
    # continuum gels
    scenes["PA gel (continuum lattice)"] = ecm_layers(
        L.build_continuum_ecm(L.get_spec("pa_gel"), lo, hi, node_spacing_um=2.0, pin_faces=("z_lo",), rng=rng(6)))
    scenes["Matrigel (basement membrane lattice)"] = ecm_layers(
        L.build_continuum_ecm(L.get_spec("matrigel"), lo, hi, node_spacing_um=2.0, pin_faces=("z_lo",), rng=rng(7)))
    # composite (interpenetrating collagen + matrigel)
    comp = L.build_composite([{"material": "collagen_I", "concentration": 1.5, "alignment_S": 0.3,
                               "director": (1, 0, 0)},
                              {"material": "matrigel"}], lo, hi, dim=3, interlink_um=1.0,
                             pin_faces=("z_lo",), rng=rng(8))
    scenes["composite collagen+Matrigel"] = ecm_layers(comp)
    build_viewer(scenes, out, title="FF ECM library — network gallery")
    print(f"wrote {out}  ({len(scenes)} scenes)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--box", type=float, default=22.0)
    a = ap.parse_args()
    default_gallery(out=a.out, box=a.box)


if __name__ == "__main__":
    main()
