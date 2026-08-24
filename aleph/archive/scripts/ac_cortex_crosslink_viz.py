#!/usr/bin/env python3
r"""Visualise the cortex α-actinin CROSSLINK NETWORK — coarse vs fine mesh — so the density is legible.

WHY (PI catch, 2026-07-24). The PI flagged that the cortex crosslinks "look too few". The existing dynamic-
tension viewer (``gate_b_dynamic.html``) draws only actin nodes + myosin — NOT the ~1.4 M α-actinin crosslink
bonds — so the network density is INVISIBLE and the concern is unverifiable. This script draws the crosslinks as
EXPLICIT EDGES (a line segment between each ``xl_i`` and ``xl_j`` node), for the coarse baseline and the derived
fine config, so the mesh-density difference is directly checkable. See
``docs/v2_audit/CORTEX_MESH_FIDELITY_2026-07-24.md`` §7.

Configs (memo §7, DERIVED not tuned):
  * COARSE baseline : ``cortex_seg_um=0.5``  ``cortex_density_per_fil=20`` → ~500 nm mesh.
    Full cell (70,686 fil): 494,802 actin nodes / 1,413,720 crosslinks.
  * FINE (physiol.) : ``cortex_seg_um=0.075`` ``cortex_density_per_fil=40`` → ~75 nm mesh (Morone/Bovellan/Chugh).
    Full cell (70,686 fil): 2,898,126 actin nodes / 2,827,440 crosslinks.

CPU-ONLY FEASIBILITY (I0-A dev-Mac path). The cortex WEAVE (node positions + crosslink index arrays) is pure
host NumPy — ``ac.weave.woven_cell.weave_cell`` → ``WovenCell.to_crosslinked_cortex`` — with NO CUDA mechanics
(``build_cell`` is the CUDA-gated composition, which this script bypasses). So the crosslinked cortex is computed
directly from the weave on the Mac. Verified: full coarse = 7.5 s / ~1 GB; full fine = 135 s / ~9.7 GB.

PATCH (default). The FULL fine cortex (2.9 M nodes, 2.8 M crosslink edges → ~5.6 M line verts) is heavy to weave
AND to render legibly. To keep the coarse-vs-fine comparison apples-to-apples and legible, this renders a
REPRESENTATIVE PATCH: a moderate full-areal-density (100 filaments/µm²) cortex sphere clipped to an identical
small near-cap window at BOTH resolutions. The mesh texture (pore size ≈ ``seg_um``; crosslink spacing ≈
``length/density``) is a LOCAL, per-area quantity — statistically identical to a patch of the real 70,686-filament
R=7.4 µm cell (only the global curvature differs across the small window) — so the patch is a faithful proxy for
the density question. The full-cell counts are printed in the title. ``--full-cell`` builds the real R=7.4 µm,
70,686-filament weave and clips a real patch from it instead (slow/heavy for the fine config).

Runtime: host NumPy + matplotlib only (no Warp CUDA); Warp imports as its CPU build purely as a transitive dep of
the weave modules and is never launched. Render reuses ``ff_viewer_html.build_viewer`` (self-contained WebGL);
browser-verify each with ``scripts/browser_check.py``.
"""
from __future__ import annotations

import argparse
import dataclasses
import time
from pathlib import Path

import numpy as np

from aleph.components.weave.regions import CORTEX_REGION
from aleph.components.weave.woven_cell import weave_cell
from aleph.scripts.ac_viz_common import seg_pairs_from_offsets
from aleph.scripts.ff_viewer_html import build_viewer

# ── physiological setpoints (mirror assemble.py; kept local so this script needs NO CUDA import) ──────────
R_CELL_UM = 7.5
CORTEX_MEMBRANE_GAP_UM = 0.10
R_CORTEX_UM = R_CELL_UM - CORTEX_MEMBRANE_GAP_UM          # 7.40 µm — the real production cortex-shell radius
AREAL_DENSITY_PER_UM2 = 100.0                            # 100 filaments/µm² (= 70,686 on the 4π·7.4² shell)
LENGTH_UM = 3.0                                          # representative cortical filament contour length

# memo §7 configs
COARSE = dict(seg_um=0.5, density_per_fil=20.0, mesh_nm=500)
FINE = dict(seg_um=0.075, density_per_fil=40.0, mesh_nm=75)

# full-cell reference counts (memo §7.2; VERIFIED by this weave on the Mac) — printed so the patch is honest
FULL_CELL = {
    "coarse": {"nodes": 494_802, "xl": 1_413_720},
    "fine": {"nodes": 2_898_126, "xl": 2_827_440},
}

# α-actinin crosslink edge colour: bright cyan, held clear of the grey actin so the network reads at a glance.
XL_COLOR = "#00e5ff"
ACTIN_COLOR = "#5a6b8c"


def cortex_region(n_filaments: int, seg_um: float, density_per_fil: float, R_um: float):
    """A CORTEX region spec at a chosen filament count + geometry (mirrors ``assemble._cortex_region``).

    Threads ``seg_um`` / ``length_um`` / ``density_per_fil`` / ``R_um`` onto a COPY of ``CORTEX_REGION.arch`` via
    :func:`dataclasses.replace` (the module constant is never mutated), so the default coarse call reproduces the
    production cortex arch field-for-field.
    """
    base = CORTEX_REGION.arch
    fil = dataclasses.replace(base.filament, n_filaments=n_filaments, seg_um=seg_um, length_um=LENGTH_UM)
    xl = dataclasses.replace(base.crosslinker, density_per_fil=density_per_fil)
    arch = dataclasses.replace(base, filament=fil, crosslinker=xl, R_um=R_um)
    return dataclasses.replace(CORTEX_REGION, arch=arch)


def build_weave(seg_um: float, density_per_fil: float, n_filaments: int, R_um: float, seed: int = 0):
    """Compute the crosslinked cortex on CPU → ``(pos, fiber_offsets, xl_i, xl_j, meta)`` (no CUDA).

    Args:
        seg_um: segment rest length ℓ₀ [µm]; caps the mesh/pore size.
        density_per_fil: crosslinks per filament (n_xl = round(n_filaments·density)).
        n_filaments: number of representative cortical filaments woven on the sphere.
        R_um: cortex-shell radius [µm].
        seed: RNG seed (weave is RNG-deterministic; seed 0 matches the production build).

    Returns:
        pos ``(N,3)`` node positions [µm]; fiber_offsets ``(F+1,)``; xl_i/xl_j ``(X,)`` crosslink node index
        pairs; meta dict with counts + wall time.
    """
    t0 = time.time()
    wc = weave_cell([cortex_region(n_filaments, seg_um, density_per_fil, R_um)],
                    rng=np.random.default_rng(seed), overlap_free=False)
    cx = wc.to_crosslinked_cortex()
    meta = {"n_filaments": n_filaments, "n_nodes": int(wc.n_nodes), "n_xl": int(cx.xl_i.size),
            "R_um": R_um, "seg_um": seg_um, "density_per_fil": density_per_fil,
            "wall_s": time.time() - t0}
    return (np.ascontiguousarray(wc.pos, np.float64), np.ascontiguousarray(wc.fiber_offsets, np.int64),
            np.ascontiguousarray(cx.xl_i, np.int64), np.ascontiguousarray(cx.xl_j, np.int64), meta)


def clip_patch(pos, fiber_offsets, xl_i, xl_j, half_um: float, z_min_frac: float = 0.0):
    """Clip to a near-cap window ``|x|<half, |y|<half, z>z_min`` and REINDEX nodes/segments/crosslinks.

    A segment or crosslink is kept only when BOTH endpoints fall inside the window (so no edge dangles outside the
    patch). Returns ``(pos_p, seg_p, xl_p, kept_mask)`` in patch-local node indices.
    """
    z_min = z_min_frac * float(np.abs(pos[:, 2]).max())
    inside = (np.abs(pos[:, 0]) < half_um) & (np.abs(pos[:, 1]) < half_um) & (pos[:, 2] > z_min)
    # reindex map: global node id → patch-local id (or -1)
    remap = np.full(pos.shape[0], -1, np.int64)
    remap[inside] = np.arange(int(inside.sum()))
    pos_p = pos[inside]
    seg = seg_pairs_from_offsets(fiber_offsets, pos.shape[0])
    seg_keep = inside[seg[:, 0]] & inside[seg[:, 1]]
    seg_p = remap[seg[seg_keep]]
    xl = np.stack([xl_i, xl_j], axis=1)
    xl_keep = inside[xl[:, 0]] & inside[xl[:, 1]]
    xl_p = remap[xl[xl_keep]]
    return pos_p, seg_p, xl_p, inside


def scene_layers(pos, seg, xl, *, label: str, xl_opacity: float, actin_opacity: float) -> list[dict]:
    """Faint actin filament segments + α-actinin crosslink EDGES for one resolution.

    The viewer draws every line at 1 px (three.js r128 ignores WebGL line width), and the auto-fit camera scales
    the patch to fill the frame — so with thousands of edges tiled over a thin curved sheet a HIGH opacity
    saturates to a solid fill and the coarse-vs-fine difference vanishes. LOW opacity instead lets the accumulated
    brightness encode LOCAL EDGE DENSITY: the coarse ~500 nm mesh reads as sparse strands with dark pores between
    them, the fine ~75 nm mesh as a brighter, denser web — which is exactly the mesh-density difference to show.
    """
    layers: list[dict] = []
    if seg.size:
        layers.append({"name": f"{label} · cortex F-actin ({seg.shape[0]:,} seg, faint)", "kind": "lines",
                       "verts": pos[seg].reshape(-1, 3), "color": ACTIN_COLOR, "swatch": ACTIN_COLOR,
                       "group": "compartment", "size": 1.0, "opacity": actin_opacity, "clip": False})
    if xl.size:
        layers.append({"name": f"{label} · α-actinin CROSSLINKS ({xl.shape[0]:,} bonds)", "kind": "lines",
                       "verts": pos[xl].reshape(-1, 3), "color": XL_COLOR, "swatch": XL_COLOR,
                       "group": "connector", "size": 1.0, "opacity": xl_opacity, "on_top": True})
    return layers


def _fmt_full(kind: str) -> str:
    f = FULL_CELL[kind]
    return f"{f['nodes']:,} nodes / {f['xl']:,} xl"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=Path("aleph/outputs/cortex_mesh/cortex_crosslink_coarse_vs_fine.html"),
                    help="output HTML path")
    ap.add_argument("--patch-half-um", type=float, default=1.0, help="half-width of the square near-cap window [µm]")
    ap.add_argument("--patch-R-um", type=float, default=3.0,
                    help="radius of the moderate full-density proxy sphere the patch is cut from [µm]")
    ap.add_argument("--full-cell", action="store_true",
                    help="build the REAL 70,686-filament R=7.4µm cortex and clip a real patch (slow: fine ~135s/~10GB)")
    ap.add_argument("--xl-opacity", type=float, default=0.16,
                    help="crosslink-edge opacity (LOW so accumulated brightness encodes local edge density)")
    ap.add_argument("--actin-opacity", type=float, default=0.05, help="faint actin filament opacity")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.full_cell:
        R_um = R_CORTEX_UM
        n_fil = int(round(AREAL_DENSITY_PER_UM2 * 4.0 * np.pi * R_um * R_um))  # = 70,686
        mode = f"REAL full cell (R={R_um:g}µm, {n_fil:,} filaments) → clipped patch"
    else:
        R_um = args.patch_R_um
        n_fil = int(round(AREAL_DENSITY_PER_UM2 * 4.0 * np.pi * R_um * R_um))
        mode = (f"representative PATCH (full-density R={R_um:g}µm proxy sphere, {n_fil:,} filaments — "
                f"local mesh identical to the real cell)")

    half = args.patch_half_um
    scenes: dict[str, list] = {}
    for kind, cfg in (("coarse", COARSE), ("fine", FINE)):
        pos, off, xi, xj, meta = build_weave(cfg["seg_um"], cfg["density_per_fil"], n_fil, R_um, args.seed)
        pos_p, seg_p, xl_p, inside = clip_patch(pos, off, xi, xj, half)
        label = f"{kind.upper()} seg={cfg['seg_um']}µm ρ={cfg['density_per_fil']:.0f}/fil (~{cfg['mesh_nm']}nm mesh)"
        key = (f"◆ {kind.upper()} ~{cfg['mesh_nm']}nm mesh · {xl_p.shape[0]:,} xl in {2*half:g}×{2*half:g}µm patch "
               f"· FULL CELL {_fmt_full(kind)}")
        scenes[key] = scene_layers(pos_p, seg_p, xl_p, label=label,
                                   xl_opacity=args.xl_opacity, actin_opacity=args.actin_opacity)
        print(f"[{kind:6s}] weave {meta['n_nodes']:,} nodes / {meta['n_xl']:,} xl in {meta['wall_s']:.1f}s "
              f"→ PATCH {int(inside.sum()):,} nodes / {xl_p.shape[0]:,} xl / {seg_p.shape[0]:,} actin seg")

    title = (f"Cortex α-actinin CROSSLINK NETWORK — COARSE (~500nm) vs FINE (~75nm) mesh · {mode} · "
             f"{2*half:g}×{2*half:g}µm cap window · FULL-CELL coarse {_fmt_full('coarse')} · fine {_fmt_full('fine')} "
             f"· units µm · FULL-RES (no downsample) · select COARSE/FINE in the preset dropdown")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    build_viewer(scenes, str(args.out), title=title, cbars=None)
    print(f"[OK] wrote {args.out}")
    print("     scenes:", *scenes.keys(), sep="\n       ")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
