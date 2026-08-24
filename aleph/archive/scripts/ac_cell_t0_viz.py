"""FIRST interactive 3-D view of the new Active Cell (ac) engine — the WHOLE cell at physiological t0.

The professor's plan ("initialize at the physiological operating point") wants to SEE the composed cell,
not just the analytic gate-plots (Terzaghi curves, WCA potentials) the ac tracks have emitted so far. This
net-new script assembles the STATIC resting-baseline GEOMETRY that each ac track / reused ff engine builds
and emits ONE self-contained interactive WebGL HTML (rotate / zoom / pan + a live cut plane + a scene
dropdown), reusing the existing ``ff_viewer_html.build_viewer`` three.js infrastructure (the same viewer
behind ff_resting_viewer / dcm_mesh_viewer).

Composed compartments (all at the MCF7 R≈7.5 µm suspended resting setpoint):
  * plasma MEMBRANE — the Helfrich icosphere lipid sheet (ff/membrane_surface.build_membrane_mesh).
  * cortical F-ACTIN — the crosslinked cortical mesh (ff/gamma_floor.build_crosslinked_cortex) at the FULL
    physiological population (100 µm⁻² × 4π·7.5² = 70,686 filaments), rendered FULL-RES / NO downsampling
    (feedback-viewer-no-downsample: ~483 k nodes is fine, big HTML is fine).
  * deformable-mesh NUCLEUS — the lamina icosphere (ac/nucleus/geometry.build_oblate_mesh), placed inside,
    off-centre (illustrative eccentricity; the MCF7 oblate aspect is a ratified GAP → rendered ROUND, the
    honest suspended-resting shape; the adherent flatten is an emergent I7 output, NOT imposed here).
  * NMII MINIFILAMENTS — the head-resolved Stam-Hocky backbone+heads topology (ac/motor
    minifilament_topology + minifilament_warp.build_minifilament_nodes) seeded at the cortex's OWN myosin
    link sites (cortex.myo_i/myo_j), coloured backbone (gold) vs heads (red).

STATIC GEOMETRY ONLY. This is numpy on the dev Mac (macOS/CPU) — NO Warp/CUDA kernel is launched; the t0
resting shell needs no solve (fibers are constructed on great circles at R, the resting geometry). The
DYNAMICS / physics-field views (pore-pressure NG-10, myosin power-stroke, flatten) need the native Warp run
on the gbook A5000 and are out of scope here.

    PYTHONPATH=/Users/sw1/ffn_cellsim python aleph/scripts/ac_cell_t0_viz.py
        -> aleph/outputs/ac/cell_t0/ac_cell_t0.html  (+ browser-verified screenshots)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from aleph.components.motor.minifilament_topology import MinifilamentTopology
from aleph.components.motor.minifilament_warp import build_minifilament_nodes
from aleph.components.nucleus.geometry import build_oblate_mesh
from aleph.laws.gamma_floor import PROD_N_FIL, PROD_N_MYO, PROD_N_XL, build_crosslinked_cortex
from aleph.laws.membrane_surface import build_membrane_mesh
from aleph.scripts.ff_cell_morphology import uv_sphere
from aleph.scripts.ff_viewer_html import build_viewer

R_CELL = 7.5                     # µm — MCF7 radius (Wagner 2011; CortexParams / membrane default)
R_NUC = 0.70 * R_CELL            # µm — nucleus radius, FF resting-viewer convention (~0.7·R_cell)
NUC_OFFSET = np.array([1.1, 0.0, -0.9])   # µm — ILLUSTRATIVE off-centre placement (see module docstring)

# NMII bipolar minifilament reference topology (I0-B3 GAP reference values, minifilament_topology docstring):
# n_bb=14, H=10 (Stam-Hocky/AFINES coarse anchor, = NMIIA_HEADS_PER_SIDE), L_bb=0.301 µm (Billington 2013 EM),
# head offset 0.200 µm (archived). True-scale (~0.3 µm) — small dashes on the 7.5 µm shell, as they should be.
NMII_TOPO = MinifilamentTopology(n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.200)

COL = {
    "cortex": "#8fbff0",     # light blue — cortical F-actin
    "xlink": "#5fd8a8",      # teal — α-actinin/filamin crosslinkers
    "membrane": "#7db8e8",   # blue translucent — plasma membrane
    "nucleus": "#d17fe0",    # violet — nuclear lamina envelope
    "cyto": "#8fd0c8",       # faint teal — implicit turgor interior
    "nmii_bb": "#f4a742",    # gold — NMII backbone
    "nmii_head": "#ff3b30",  # red — NMII motor heads
}

OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "cell_t0"


def _fiber_segments(pos: np.ndarray, offsets: np.ndarray) -> np.ndarray:
    """Vectorized internal segment endpoint pairs of every fiber -> (2·Nseg, 3) for THREE.LineSegments.

    A node is a segment START unless it is the last node of its fiber. Fully vectorized (no per-fiber
    Python loop) so the full 70,686-filament / 424 k-segment cortex builds instantly.
    """
    n = pos.shape[0]
    is_last = np.zeros(n, bool)
    is_last[offsets[1:] - 1] = True
    starts = np.arange(n)[~is_last]
    seg = np.stack([starts, starts + 1], axis=1)            # (Nseg, 2) node-index pairs
    return pos[seg].reshape(-1, 3).astype(np.float32)       # (2·Nseg, 3) endpoint pairs


def _xlink_segments(pos: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Crosslinker link endpoint pairs -> (2·Nxl, 3)."""
    seg = np.stack([i, j], axis=1)
    return pos[seg].reshape(-1, 3).astype(np.float32)


def _build_nmii(pos: np.ndarray, myo_i: np.ndarray, myo_j: np.ndarray):
    """Build head-resolved NMII minifilaments seeded at the cortex's OWN myosin link sites.

    Each minifilament centre = the myosin link midpoint; backbone axis = the link direction (tangent to the
    shell). Returns (backbone_segments (2·Sbb,3), head_bond_segments (2·Nheads,3), n_minifil, n_heads_total).

    Both are rendered as LINE segments (not points): the viewer forces all point layers to a single world-unit
    ``sizeAttenuation`` size, which would make ~9 k heads giant occluding blobs. The head↔backbone BOND stubs
    read as thin explicit heads AND faithfully show the head-resolved attachment topology (each of the 2H
    heads tethered to one backbone bead).
    """
    bb_segs, head_bond_segs = [], []
    n_heads_total = 0
    for k in range(myo_i.shape[0]):
        centre = 0.5 * (pos[myo_i[k]] + pos[myo_j[k]])
        axis = pos[myo_j[k]] - pos[myo_i[k]]
        d = build_minifilament_nodes(NMII_TOPO, centre.astype(np.float64), axis.astype(np.float64))
        p = d["positions"]
        bb_segs.append(p[d["backbone_bonds"]].reshape(-1, 3))    # (n_bb-1, 2) backbone chain
        head_bond_segs.append(p[d["head_bonds"]].reshape(-1, 3))  # (2H, 2) head<->backbone stubs
        n_heads_total += d["head_node"].shape[0]
    bb = (np.concatenate(bb_segs, 0).astype(np.float32) if bb_segs else np.zeros((0, 3), np.float32))
    hb = (np.concatenate(head_bond_segs, 0).astype(np.float32) if head_bond_segs else np.zeros((0, 3), np.float32))
    return bb, hb, myo_i.shape[0], n_heads_total


def build(out_html: Path, n_filaments: int = PROD_N_FIL) -> tuple[str, dict]:
    """Assemble the composed t0 cell geometry and emit the interactive HTML. Returns (path, counts)."""
    # ---- cortical F-actin (crosslinked cortex; static resting shell, no Warp solve) ----------------
    n_xl = n_filaments if n_filaments == PROD_N_FIL else n_filaments
    cortex = build_crosslinked_cortex(n_filaments=n_filaments, n_xl=n_xl, n_myo=PROD_N_MYO)
    pos = np.asarray(cortex.net.pos, np.float64)
    centroid = pos.mean(0)
    off = np.asarray(cortex.net.fiber_offsets, np.int64)
    cortex_seg = _fiber_segments(pos, off)
    xlink_seg = _xlink_segments(pos, cortex.xl_i, cortex.xl_j)
    n_seg = cortex_seg.shape[0] // 2
    n_xl_rendered = xlink_seg.shape[0] // 2

    # ---- plasma membrane (Helfrich icosphere) -------------------------------------------------------
    mm = build_membrane_mesh(R_CELL, subdivisions=4, centre=tuple(centroid))
    mem_v = np.asarray(mm.verts, np.float32)
    mem_f = np.asarray(mm.faces, np.int64)

    # ---- deformable-mesh nucleus (lamina envelope; ROUND suspended resting, off-centre) -------------
    nuc_centre = centroid + NUC_OFFSET
    nuc_v, nuc_f = build_oblate_mesh(R_NUC, aspect=1.0, subdivisions=3, centre=tuple(nuc_centre))
    nuc_v = nuc_v.astype(np.float32)

    # ---- implicit turgor cytoplasm interior (faint ghost sphere just inside the cortex) -------------
    cyto_v, cyto_f = uv_sphere(centroid, 0.90 * R_CELL, nu=32, nv=20)

    # ---- head-resolved NMII minifilaments (seeded at the cortex's own myosin sites) -----------------
    nmii_bb, nmii_hb, n_minifil, n_heads = _build_nmii(pos, cortex.myo_i, cortex.myo_j)

    counts = {
        "membrane_verts": int(mem_v.shape[0]), "membrane_faces": int(mem_f.shape[0]),
        "cortex_filaments": int(cortex.net.n_fibers), "cortex_nodes": int(cortex.net.n_nodes),
        "cortex_segments": int(n_seg), "crosslinkers": int(n_xl_rendered),
        "nucleus_verts": int(nuc_v.shape[0]), "nucleus_faces": int(nuc_f.shape[0]),
        "nmii_minifilaments": int(n_minifil), "nmii_heads": int(n_heads),
        "nmii_backbone_beads": int(n_minifil * NMII_TOPO.n_bb),
    }

    # ---- compose scenes (reuse ff_viewer_html.build_viewer) -----------------------------------------
    cortex_full = {"name": f"cortical F-actin ({counts['cortex_filaments']:,} filaments · "
                           f"{counts['cortex_segments']:,} segments · FULL-RES)", "kind": "lines",
                   "verts": cortex_seg, "color": COL["cortex"], "size": 1.1, "opacity": 0.42, "clip": True}
    membrane = {"name": f"plasma membrane (Helfrich icosphere · {counts['membrane_verts']} nodes)",
                "kind": "mesh", "verts": mem_v, "faces": mem_f, "color": COL["membrane"],
                "opacity": 0.11, "clip": True}
    membrane_ghost = dict(membrane, opacity=0.05, name="plasma membrane (envelope, ghost)")
    nucleus = {"name": f"nuclear lamina envelope ({counts['nucleus_verts']} nodes · R≈{R_NUC:.1f}µm · "
                       f"off-centre)", "kind": "mesh", "verts": nuc_v, "faces": nuc_f,
               "color": COL["nucleus"], "opacity": 0.92}
    cyto = {"name": "cytoplasm (implicit turgor interior ΔP≈40 Pa)", "kind": "mesh",
            "verts": cyto_v, "faces": cyto_f, "color": COL["cyto"], "opacity": 0.05, "clip": True}
    nmii_bb_layer = {"name": f"NMII backbones ({n_minifil} minifilaments · Stam-Hocky, Nie-2015 density)",
                     "kind": "lines", "verts": nmii_bb, "color": COL["nmii_bb"], "size": 2.6,
                     "opacity": 1.0, "on_top": True}
    nmii_head_layer = {"name": f"NMII motor heads ({n_heads} explicit heads · Hill FV / Bell kinetics)",
                       "kind": "lines", "verts": nmii_hb, "color": COL["nmii_head"], "size": 1.6,
                       "opacity": 1.0, "on_top": True}

    # cut-away (interior visible by default): keep the BACK hemisphere of the cortex/NMII so nucleus shows
    cx = float(centroid[0])
    seg_mid_back = 0.5 * (cortex_seg.reshape(-1, 2, 3)[:, 0, 0] + cortex_seg.reshape(-1, 2, 3)[:, 1, 0]) <= cx
    cortex_back = cortex_seg.reshape(-1, 2, 3)[seg_mid_back].reshape(-1, 3)
    bb_back_mask = 0.5 * (nmii_bb.reshape(-1, 2, 3)[:, 0, 0] + nmii_bb.reshape(-1, 2, 3)[:, 1, 0]) <= cx
    nmii_bb_back = nmii_bb.reshape(-1, 2, 3)[bb_back_mask].reshape(-1, 3)
    hb_back_mask = 0.5 * (nmii_hb.reshape(-1, 2, 3)[:, 0, 0] + nmii_hb.reshape(-1, 2, 3)[:, 1, 0]) <= cx
    nmii_hb_back = nmii_hb.reshape(-1, 2, 3)[hb_back_mask].reshape(-1, 3)
    cyto_f_back = cyto_f[cyto_v[cyto_f].mean(1)[:, 0] <= cx]

    cutaway = [
        {"name": f"cortex (back hemisphere cut-away · {int(seg_mid_back.sum()):,} segments)", "kind": "lines",
         "verts": cortex_back, "color": COL["cortex"], "size": 1.2, "opacity": 0.62},
        {"name": "cytoplasm (turgor interior, back half)", "kind": "mesh", "verts": cyto_v,
         "faces": cyto_f_back, "color": COL["cyto"], "opacity": 0.09},
        nucleus,
        {"name": nmii_bb_layer["name"] + " (back half)", "kind": "lines", "verts": nmii_bb_back,
         "color": COL["nmii_bb"], "size": 2.6, "opacity": 1.0, "on_top": True},
        {"name": nmii_head_layer["name"] + " (back half)", "kind": "lines", "verts": nmii_hb_back,
         "color": COL["nmii_head"], "size": 1.6, "opacity": 1.0, "on_top": True},
        membrane_ghost,
    ]

    whole = [cortex_full, cyto, membrane, nucleus, nmii_bb_layer, nmii_head_layer]

    motors = [
        {"name": f"cortical F-actin (context, faint · {counts['cortex_filaments']:,} filaments)",
         "kind": "lines", "verts": cortex_seg, "color": COL["cortex"], "size": 0.9, "opacity": 0.16,
         "clip": True},
        {"name": f"crosslinkers ({counts['crosslinkers']:,} α-actinin/filamin links)", "kind": "lines",
         "verts": xlink_seg, "color": COL["xlink"], "size": 1.4, "opacity": 0.5, "clip": True},
        nmii_bb_layer, nmii_head_layer, membrane_ghost,
    ]

    scenes = {
        "cut-away · interior compartments": cutaway,
        f"whole cell · {counts['cortex_filaments']:,} filaments (use CUT slider)": whole,
        "crosslinked cortex + NMII motors": motors,
    }

    title = (f"AC engine — MCF7 cell at physiological t0 (SUSPENDED resting baseline) · R={R_CELL:.1f} µm · "
             f"membrane(Helfrich) + cortex({counts['cortex_filaments']:,} F-actin, FULL-RES) + "
             f"nucleus(lamina) + NMII({n_minifil} head-resolved minifilaments) · units µm · STATIC geometry")
    build_viewer(scenes, out=str(out_html), title=title)
    return str(out_html), counts


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    out_html = OUT / "ac_cell_t0.html"
    path, counts = build(out_html)
    print(f"wrote {path}")
    for k, v in counts.items():
        print(f"  {k:22s} {v:,}")


if __name__ == "__main__":
    main()
