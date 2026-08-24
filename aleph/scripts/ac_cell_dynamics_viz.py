"""Interactive 3-D render of the FLOWING physics of the fluid-coupled Active Cell (NG-10-style dynamics).

This is the poroelastic-FLOW counterpart of ``ac_cell_assembled_viz.py``. There, the resting cell ran at a
UNIFORM Π₀ ⇒ grad p = 0 (a balloon: one pressure everywhere, no flow). HERE the geometry AND the fields come
from the FSI-COUPLED deformation run dumped on the gbook A5000
(``aleph/outputs/ac/fsi/coupled_state.npz``): a controlled compression squeezed the poroelastic solid, a
real pore-pressure GRADIENT developed, and the pore fluid genuinely FLOWS down that gradient. This script
(numpy on the dev Mac — NO Warp/CUDA) reads the dumped conservative-field GRID (43³, dx=0.5 µm) + the actin
node cloud and emits ONE self-contained interactive WebGL HTML via the reused
``ff_viewer_html.build_viewer`` (OrbitControls rotate/zoom/pan + a live CUT plane + a scene dropdown + a
labelled colour scale), so the Lead/PI can SEE the difference between a balloon and a poroelastic cell.

What the milestone-2 coupled dump actually contains (masked to the cell interior, R=7.5 µm sphere):
  * pore pressure  p [Pa]      — excess over the resting Π₀=40 Pa spans ~17 Pa (centre) → ~427 Pa (shell):
                                 a genuine radial GRADIENT, NOT the resting uniform ΔP.
  * pore-fluid velocity v_f    — |v_f| up to ~1.79 µm/s, pointing radially INWARD (100 % inward): the cytosol
    [µm/s]                       flows down the pressure gradient as the compressed solid squeezes it inward.
  * Darcy discharge q [µm/s]   — the volumetric flux |q| up to ~0.89 µm/s (q = φ·v_f up to the φ GAP).
  * ∇·v_s (div_vs)             — the solid dilatation-rate SINK (≤0) that DRIVES the flow (the compression).
  * actin nodes + composed |F| — the 70 686-filament cortex, coloured by the composed per-node force [pN].

The pore pressure/velocity live on the conservative field GRID (a 3-D array over the cell domain, masked to
the sphere); the actin is a 494 802-node cloud. The grid is slab-cut / clip-tagged so the INTERIOR gradient
is visible — the whole point of the figure is that the interior is NOT one flat colour.

Colour ramps (perceptually-ordered **turbo**, §1.9 integrity — no axis truncation, units annotated on every
colour bar): pressure & velocity & q are LINEAR (smooth fields, a linear ramp reads the gradient honestly);
the composed actin |F| is LOG₁₀ (it spans >4 decades — resting floor vs excluded-volume hotspots).

CAVEAT (honest, printed in the title): this is the ISOLATED fluid-payload run — STERIC OFF and a CONTROLLED
ε=1 % compression, NOT the full active adherent cell. The porosity φ=0.5 is a provisional GAP that linearly
scales v_f = q/φ. The dump carries no explicit NMII head positions, so motors are reported as a count, not
drawn (the composed |F| already folds in the myosin contribution).

    PYTHONPATH=/Users/sw1/ffn_cellsim python aleph/scripts/ac_cell_dynamics_viz.py
        -> aleph/outputs/ac/cell_dynamics/ac_cell_dynamics.html  (+ browser-verified screenshots)
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np

from aleph.scripts.ff_cell_morphology import uv_sphere
from aleph.scripts.ff_viewer_html import build_viewer

R_CELL = 7.5                          # µm — MCF7 radius (assemble.R_CELL_UM)
PI0_PA = 40.0                         # Pa — resting osmotic setpoint Π₀ (the balloon pressure); p_excess = p - Π₀
NPZ = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "fsi" / "coupled_state.npz"
OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "cell_dynamics"

_TURBO = matplotlib.colormaps["turbo"]

COL = {
    "envelope": "#5a6b8c",   # faint blue — R=7.5 µm cell boundary reference
    "actin": "#39435a",      # muted slate — faint cortex context
}


# --------------------------------------------------------------------------------------------------------- #
# colour helpers (turbo; linear for the smooth fields, log₁₀ for the >4-decade force)
# --------------------------------------------------------------------------------------------------------- #
def _lin_rgb(val: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    """Map a scalar field → (N,3) uint8 turbo colours on a clamped LINEAR axis."""
    t = (np.clip(val.astype(np.float64), vmin, vmax) - vmin) / (vmax - vmin + 1e-30)
    return (_TURBO(t)[:, :3] * 255.0 + 0.5).astype(np.uint8)


def _log_rgb(fmag: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    """Map |F| [pN] → (N,3) uint8 turbo colours on a clamped log₁₀ axis."""
    f = np.clip(fmag.astype(np.float64), vmin, vmax)
    t = (np.log10(f) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))
    return (_TURBO(t)[:, :3] * 255.0 + 0.5).astype(np.uint8)


def _cbar(label: str, unit: str, vmin: float, vmax: float, log: bool = False, n: int = 7) -> dict:
    """Turbo colour-bar spec for build_viewer (grad bottom→top; ticks span [vmin, vmax])."""
    grad = [matplotlib.colors.to_hex(_TURBO(i / (n - 1))) for i in range(n)]
    lo, hi = (np.log10(vmin), np.log10(vmax)) if log else (vmin, vmax)
    return {"label": label, "unit": (("log₁₀ " + unit) if log else unit), "grad": grad,
            "lo": float(lo), "hi": float(hi)}


def _fiber_seg_idx(offsets: np.ndarray, n_nodes: int) -> np.ndarray:
    """Internal segment node-index pairs of every fiber → (Nseg, 2) (vectorized; last node of a fiber skipped)."""
    is_last = np.zeros(n_nodes, bool)
    is_last[offsets[1:] - 1] = True
    starts = np.arange(n_nodes)[~is_last]
    return np.stack([starts, starts + 1], axis=1)


def _grid_coords(shape, origin, dx) -> np.ndarray:
    """Physical (µm) coordinate of every grid node → (nx,ny,nz,3)."""
    nx, ny, nz = shape
    ax = origin[0] + np.arange(nx) * dx
    ay = origin[1] + np.arange(ny) * dx
    az = origin[2] + np.arange(nz) * dx
    X, Y, Z = np.meshgrid(ax, ay, az, indexing="ij")
    return np.stack([X, Y, Z], axis=-1).astype(np.float32)


def _arrows(pts: np.ndarray, vecs: np.ndarray, scale: float) -> np.ndarray:
    """Vector field → line-segment endpoint pairs P→P+scale·V; returns (2·M, 3) [p0,e0,p1,e1,…]."""
    ends = pts + scale * vecs
    seg = np.stack([pts, ends], axis=1)         # (M, 2, 3)
    return seg.reshape(-1, 3).astype(np.float32)


def build(out_html: Path) -> tuple[str, dict]:
    """Read the FSI-coupled field dump and emit the flowing-physics interactive HTML."""
    d = np.load(NPZ, allow_pickle=True)
    report = json.loads(str(d["report_json"]))
    ledger = json.loads(str(d["ledger_json"]))

    # ---- fields on the conservative grid, masked to the cell interior --------------------------------------
    shape = tuple(int(v) for v in d["grid_shape"])
    dx = float(d["grid_dx"]); origin = np.asarray(d["grid_origin"], float)
    XYZ = _grid_coords(shape, origin, dx)                       # (nx,ny,nz,3)
    mask = d["mask"] > 0                                        # (nx,ny,nz) cell interior
    p = d["p_field"].astype(np.float64)                         # (nx,ny,nz) pore pressure [Pa]
    p_ex = p - PI0_PA                                           # excess over resting Π₀ [Pa]
    vf = d["v_f_field"].astype(np.float64)                      # (nx,ny,nz,3) pore-fluid velocity [µm/s]
    q = d["q_field"].astype(np.float64)                         # (nx,ny,nz,3) Darcy discharge [µm/s]
    dvs = d["div_vs"].astype(np.float64)                        # (nx,ny,nz) solid dilatation rate [1/s]
    vf_mag = np.linalg.norm(vf, axis=-1)
    q_mag = np.linalg.norm(q, axis=-1)

    # masked stats (drive the colour ranges + the report)
    pex_m = p_ex[mask]; vfm_m = vf_mag[mask]; qm_m = q_mag[mask]; dvs_m = dvs[mask]
    PEX_MAX = float(np.ceil(pex_m.max() / 10) * 10)            # ~430 Pa
    VF_MAX = float(np.ceil(vfm_m.max() * 10) / 10)            # ~1.8 µm/s
    Q_MAX = float(np.ceil(qm_m.max() * 100) / 100)           # ~0.89 µm/s
    DVS_MIN = float(np.floor(dvs_m.min() * 10) / 10)         # ~-2.8 (sink)

    # radial inflow diagnostic (for the caption): fraction of masked cells whose v_f points inward
    r = np.linalg.norm(XYZ, axis=-1)
    rhat = XYZ / (r[..., None] + 1e-9)
    inward_frac = float((np.sum(vf * rhat, axis=-1)[mask] < 0).mean())

    idx = np.argwhere(mask)                                     # (M,3) grid indices of interior cells
    gp = XYZ[mask]                                              # (M,3) interior grid points [µm]
    zc = gp[:, 2]
    slab = np.abs(zc) < 0.75                                    # central z-slab (3 planes, ~1.5 µm thick)
    M = int(mask.sum())

    # ---- reference envelope (R=7.5 µm cell boundary, faint) ----------------------------------------------
    ref_v, ref_f = uv_sphere(np.zeros(3), R_CELL, nu=48, nv=28)
    ref_layer = {"name": "R=7.5 µm cell boundary (reference envelope)", "kind": "mesh",
                 "verts": ref_v, "faces": ref_f, "color": COL["envelope"], "opacity": 0.05, "clip": True}

    # ---- arrow scales (max arrow ≈ 1.2 µm so the flow reads without over-cluttering the slab) -------------
    SV = 1.2 / max(VF_MAX, 1e-9)
    SQ = 1.2 / max(Q_MAX, 1e-9)

    # ---- PRESSURE colours (linear, 0 → PEX_MAX Pa; blue≈resting, red=high excess) ------------------------
    p_rgb_all = _lin_rgb(p_ex[mask], 0.0, PEX_MAX)             # (M,3)
    p_rgb_slab = _lin_rgb(p_ex[mask][slab], 0.0, PEX_MAX)

    # ================= SCENE 1: pressure central slab (the money-shot cross-section) ======================
    n_slab = int(slab.sum())
    s1 = [
        {"name": f"pore-pressure p_excess · central z-slab ({n_slab:,} field cells · linear turbo)",
         "kind": "points", "verts": gp[slab], "color": "#ffffff", "size": 4.0,
         "color_frames": [p_rgb_slab]},
        ref_layer,
    ]

    # ================= SCENE 2: pressure whole interior (clip to reveal the 3-D gradient) =================
    s2 = [
        {"name": f"pore-pressure p_excess · interior ({M:,} field cells · CLIP to cut · linear turbo)",
         "kind": "points", "verts": gp, "color": "#ffffff", "size": 3.0, "clip": True,
         "color_frames": [p_rgb_all]},
        ref_layer,
    ]

    # ================= SCENE 3: pore-fluid velocity v_f — slab arrows (the FLOW) ===========================
    vf_rgb_slab = _lin_rgb(vf_mag[mask][slab], 0.0, VF_MAX)
    vf_arr_slab = _arrows(gp[slab], vf[mask][slab], SV)
    s3 = [
        {"name": f"pore-fluid velocity v_f · central z-slab ({n_slab:,} arrows · P→P+{SV:.2f}·v_f · "
                 f"linear turbo)", "kind": "lines", "verts": vf_arr_slab, "color": "#ffffff",
         "size": 1.5, "opacity": 0.95, "color_frames": [np.repeat(vf_rgb_slab, 2, axis=0)]},
        ref_layer,
    ]

    # ================= SCENE 4: pore-fluid velocity v_f — whole interior arrows (clip) =====================
    vf_rgb_all = _lin_rgb(vf_mag[mask], 0.0, VF_MAX)
    vf_arr_all = _arrows(gp, vf[mask], SV)
    s4 = [
        {"name": f"pore-fluid velocity v_f · interior ({M:,} arrows · CLIP to cut · linear turbo)",
         "kind": "lines", "verts": vf_arr_all, "color": "#ffffff", "size": 1.2, "opacity": 0.9,
         "clip": True, "color_frames": [np.repeat(vf_rgb_all, 2, axis=0)]},
        ref_layer,
    ]

    # ================= SCENE 5: Darcy discharge q — slab arrows ============================================
    q_rgb_slab = _lin_rgb(q_mag[mask][slab], 0.0, Q_MAX)
    q_arr_slab = _arrows(gp[slab], q[mask][slab], SQ)
    s5 = [
        {"name": f"Darcy discharge q · central z-slab ({n_slab:,} arrows · P→P+{SQ:.2f}·q · linear turbo)",
         "kind": "lines", "verts": q_arr_slab, "color": "#ffffff", "size": 1.5, "opacity": 0.95,
         "color_frames": [np.repeat(q_rgb_slab, 2, axis=0)]},
        ref_layer,
    ]

    # ================= SCENE 6: solid dilatation-rate ∇·v_s (the sink that DRIVES the flow) ================
    #   div_vs ≤ 0 everywhere (compression). Colour on |div_vs| 0→|DVS_MIN| so the strongest sink is red.
    dvs_rgb_slab = _lin_rgb(-dvs[mask][slab], 0.0, -DVS_MIN)
    s6 = [
        {"name": f"solid dilatation-rate −∇·v_s (compression SINK, drives inflow) · z-slab ({n_slab:,} "
                 f"cells · linear turbo)", "kind": "points", "verts": gp[slab], "color": "#ffffff",
         "size": 4.0, "color_frames": [dvs_rgb_slab]},
        ref_layer,
    ]

    # ================= SCENE 7: cortex actin coloured by composed per-node |F| =============================
    n_actin = int(report["n_actin"])
    pos_actin = d["actin_pos"].astype(np.float32)              # (n_actin, 3)
    f_mag = d["f_node_mag"].astype(np.float64)                 # (n_actin,) composed |F| [pN]
    offsets = d["fiber_offsets"].astype(np.int64)
    seg = _fiber_seg_idx(offsets, n_actin)                     # (Nseg, 2)
    seg_verts = pos_actin[seg].reshape(-1, 3)
    F_VMIN, F_VMAX = 1.0e-2, 2.5e3
    f_rgb = _log_rgb(f_mag, F_VMIN, F_VMAX)
    seg_col = f_rgb[seg].reshape(-1, 3)
    n_seg = int(seg.shape[0])
    n_bound = int(report.get("n_bound_heads", 0))
    n_heads = int(ledger.get("myosin_n_heads", 0))
    n_minifil = int(ledger.get("myosin_n_minifilaments", 0))
    s7 = [
        {"name": f"cortex F-actin · composed |F| ({report['n_filaments']:,} fil · {n_seg:,} seg · "
                 f"FULL-RES · log turbo · folds in myosin)", "kind": "lines", "verts": seg_verts,
         "color": "#ffffff", "size": 1.1, "opacity": 0.55, "clip": True, "color_frames": [seg_col]},
        ref_layer,
    ]

    # ---- scene registry + colour bars --------------------------------------------------------------------
    sc1 = "① pore-pressure p_excess · central slab [Pa]"
    sc2 = "② pore-pressure p_excess · interior (CLIP) [Pa]"
    sc3 = "③ pore-fluid velocity |v_f| · slab arrows [µm/s]"
    sc4 = "④ pore-fluid velocity |v_f| · interior (CLIP) [µm/s]"
    sc5 = "⑤ Darcy discharge |q| · slab arrows [µm/s]"
    sc6 = "⑥ solid dilatation −∇·v_s · slab (flow driver) [1/s]"
    sc7 = "⑦ cortex actin · composed |F| [pN, log]"
    scenes = {sc1: s1, sc2: s2, sc3: s3, sc4: s4, sc5: s5, sc6: s6, sc7: s7}
    cbars = {
        sc1: _cbar("pore-pressure excess p−Π₀", "Pa", 0.0, PEX_MAX),
        sc2: _cbar("pore-pressure excess p−Π₀", "Pa", 0.0, PEX_MAX),
        sc3: _cbar("pore-fluid speed |v_f|", "µm/s", 0.0, VF_MAX),
        sc4: _cbar("pore-fluid speed |v_f|", "µm/s", 0.0, VF_MAX),
        sc5: _cbar("Darcy discharge |q|", "µm/s", 0.0, Q_MAX),
        sc6: _cbar("compression sink −∇·v_s", "1/s", 0.0, -DVS_MIN),
        sc7: _cbar("composed per-node |F|", "pN", F_VMIN, F_VMAX, log=True),
    }

    title = (
        f"AC engine — fluid-coupled MCF7 cell, FLOWING physics (FSI milestone-2, coupled_state.npz) · "
        f"pore-pressure p_excess ∈ [{pex_m.min():.0f}, {pex_m.max():.0f}] Pa "
        f"(a real radial GRADIENT — NOT the resting uniform Π₀={PI0_PA:.0f} Pa balloon) · "
        f"pore-fluid |v_f| up to {vfm_m.max():.2f} µm/s, {inward_frac*100:.0f}% INWARD (cytosol flows down "
        f"grad p) · Darcy |q| up to {qm_m.max():.2f} µm/s · {report['n_filaments']:,} cortex F-actin / "
        f"{n_actin:,} nodes · {n_minifil} NMII minifil / {n_heads} heads ({n_bound} bound; positions not in "
        f"dump) · units µm · CAVEAT: isolated fluid-payload run, STERIC OFF, controlled ε="
        f"{report['eps_compression']*100:.0f}% compression, φ={report['phi_provisional_GAP']} GAP scales v_f "
        f"(NOT the full active adherent cell)"
    )
    build_viewer(scenes, out=str(out_html), title=title, cbars=cbars)

    counts = {
        "field_cells_interior": M, "field_cells_slab": n_slab,
        "p_excess_min_Pa": float(pex_m.min()), "p_excess_max_Pa": float(pex_m.max()),
        "vf_max_um_s": float(vfm_m.max()), "vf_mean_um_s": float(vfm_m.mean()),
        "q_max_um_s": float(qm_m.max()), "inward_fraction": inward_frac,
        "div_vs_min": float(dvs_m.min()), "cortex_filaments": int(report["n_filaments"]),
        "actin_nodes": n_actin, "actin_segments": n_seg, "nmii_minifilaments": n_minifil,
        "nmii_heads": n_heads, "nmii_bound_heads": n_bound,
        "f_composed_max_pN": float(f_mag.max()), "f_composed_mean_pN": float(f_mag.mean()),
    }
    return str(out_html), counts


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    out_html = OUT / "ac_cell_dynamics.html"
    path, counts = build(out_html)
    print(f"wrote {path}")
    for k, v in counts.items():
        print(f"  {k:24s} {v:,}" if isinstance(v, int) else f"  {k:24s} {v:.4g}")


if __name__ == "__main__":
    main()
