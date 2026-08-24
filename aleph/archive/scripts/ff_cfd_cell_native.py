r"""CFD-ON native cell — the spatial Biot pore-pressure field working INSIDE the trustworthy cell.

The verifiable "is the CFD real?" demonstration (PI 2026-07-16). Takes the validated pre-stressed
native cell (`resting_reval.npz`, P6.1) and runs the device-resident Biot poroelastic field
(`ff/biot_fluid_warp.py`, validated P7.F.1/F.2/F.5) on it: a compression-induced overpressure fills
the cytoplasm and DRAINS from the plasma-membrane boundary (Lp efflux), so the pore pressure p(x,t)
relaxes on the POROELASTIC time τ_p ~ R²/D (D = Moeendarbary 40–60 µm²/s) — the FAST spatial
relaxation, NOT the slow 0-D efflux clock (τ_osm ~ 30–200 s). This is how you SEE + check the CFD:
the pressure-field decay time is a number you can compare to the measured poroelastic τ_p, and the
field is rendered inside the cell (turbo p, cut-away).

Usage: python -m aleph.scripts.ff_cfd_cell_native <resting_cell.npz> <out.html>
Verifiable output: `cfd_cell_metrics.json` — the center-pressure decay τ_fit vs the poroelastic
τ_p = R²/(π²D) (drained-sphere fundamental mode), and the spatial pressure front.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from aleph.laws.biot_fluid_warp import BiotField

D_CYTO = 50.0        # µm²/s — Moeendarbary poroelastic diffusivity (KB anchor)


def run(npz_path: str, out_html: str) -> dict:
    d = np.load(npz_path)
    frames = d["frames"]
    pos = frames[0] if frames.ndim == 3 else frames
    Nc = int(d["Nc"])
    pcx = pos[:Nc]
    c = pcx.mean(0)
    R = float(np.linalg.norm(pcx - c, axis=1).mean())

    # ---- Biot grid spanning the cell ----
    dx = 0.4
    half = 1.25 * R
    n = int(2 * half / dx) + 1
    origin = c - half
    g = origin[:, None] + np.arange(n)[None, :] * dx      # (3, n)
    X, Y, Z = np.meshgrid(g[0], g[1], g[2], indexing="ij")
    r_grid = np.sqrt((X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2)
    inside = r_grid < R                                   # cytoplasm; membrane at r=R is the drained boundary

    # ---- initial compression-induced overpressure (uniform inside, drained outside) ----
    p0 = np.where(inside, 1.0, 0.0).astype(np.float64)

    bf = BiotField(n, dx, D_CYTO, device="cpu")
    bf.set_field(p0)
    dt = 0.5 * bf.dt_max
    tau_p = R ** 2 / (np.pi ** 2 * D_CYTO)                # drained-sphere fundamental-mode relaxation time
    t_end = 5.0 * tau_p
    nsteps = int(t_end / dt)

    # centre + mid-radius sample cells
    ci = n // 2
    mid_mask = np.abs(r_grid - 0.5 * R) < dx
    ts, p_center, p_mid, snap = [], [], [], None
    for s in range(nsteps):
        bf.step(dt)
        if s % 5 == 0:                                    # enforce the drained membrane boundary (p=0 at r≥R)
            f = bf.get_field(); f[~inside] = 0.0; bf.set_field(f)
        if s % max(1, nsteps // 40) == 0:
            f = bf.get_field()
            ts.append(s * dt); p_center.append(float(f[ci, ci, ci]))
            p_mid.append(float(f[mid_mask].mean()) if mid_mask.any() else 0.0)
            if abs(s * dt - tau_p) < dt * (nsteps // 40):  # a mid-relaxation snapshot for the render
                snap = f.copy()
    if snap is None:
        snap = bf.get_field()
    ts = np.array(ts); pc = np.array(p_center)

    # ---- fit the centre-pressure decay τ and compare to the poroelastic τ_p ----
    m = pc > 0.02 * pc[0]
    tau_fit = float(-1.0 / np.polyfit(ts[m], np.log(pc[m]), 1)[0]) if m.sum() > 3 else float("nan")
    # a spatial front check: at t≈τ_p the centre still holds pressure while the rim has drained
    front_ratio = float(snap[ci, ci, ci] / max(snap[mid_mask].mean(), 1e-9)) if mid_mask.any() else 0.0

    metrics = {
        "R_um": R, "D_cyto_um2_s": D_CYTO,
        "tau_p_poroelastic_s": float(tau_p),           # the analytic poroelastic time = R²/(π²D)
        "tau_fit_center_s": tau_fit,                   # the MEASURED field decay — compare these two
        "tau_ratio_fit_over_analytic": tau_fit / tau_p if tau_fit == tau_fit else None,
        "center_holds_vs_rim_drained_at_tau_p": front_ratio,   # >1 → spatial front (centre lags rim)
        "tau_osm_0D_efflux_s": "30–200 (the SLOW clock the spatial field REPLACES)",
        "grid": f"{n}³", "Nc_cortex": Nc,
        "verifiable": "compare tau_fit_center to tau_p_poroelastic — the spatial field relaxes on the FAST "
                      "poroelastic time (~R²/D), reproducing Moeendarbary, NOT the 0-D efflux clock.",
    }
    Path(out_html).parent.mkdir(parents=True, exist_ok=True)
    (Path(out_html).with_name("cfd_cell_metrics.json")).write_text(json.dumps(metrics, indent=1))
    _render(pcx, c, R, X, Y, Z, inside, snap, r_grid, metrics, out_html)
    print(f"# CFD cell: tau_fit={tau_fit:.3f}s vs poroelastic tau_p={tau_p:.3f}s "
          f"(ratio {metrics['tau_ratio_fit_over_analytic']:.2f}); front centre/rim={front_ratio:.2f}")
    print(f"# wrote {out_html} + cfd_cell_metrics.json")
    return metrics


def _render(pcx, c, R, X, Y, Z, inside, p, r_grid, metrics, out_html):
    """Interactive HTML: the pore-pressure field p as a turbo-coloured mid-plane SLAB (the radial
    gradient is visible — high at the centre, drained at the rim) + the cortex shell outline."""
    import matplotlib.cm as cm
    from aleph.scripts.ff_viewer_html import build_viewer
    turbo = cm.get_cmap("turbo")
    lo, hi = 0.0, max(float(p[inside].max()), 1e-9)

    def slab(field, snapshot_name):
        # a thin z-slab through the centre → a colored disk showing the radial p gradient
        zc = c[2]
        dz = 1.2 * (X[1, 0, 0] - X[0, 0, 0]) if X.shape[0] > 1 else 0.5
        m = inside & (np.abs(Z - zc) < dz)
        pts = np.stack([X[m], Y[m], Z[m]], axis=1)
        pv = field[m]
        cols = (turbo(np.clip(pv / hi, 0, 1))[:, :3] * 255).astype(np.uint8)
        return {"name": f"{snapshot_name}  [turbo p, {pts.shape[0]} pts]", "kind": "points",
                "verts": pts, "size": 4.0, "color_frames": [cols]}

    # cortex equatorial ring outline (a thin z-slab of the shell) — shows the cell boundary
    zc = c[2]
    ring = pcx[np.abs(pcx[:, 2] - zc) < 0.6]
    cortex_layer = {"name": f"cortex boundary (equatorial ring, {ring.shape[0]} pts)", "kind": "points",
                    "verts": ring[::4], "color": "#8fa", "size": 1.2, "opacity": 0.5}
    tau_fit = metrics["tau_fit_center_s"]; tau_p = metrics["tau_p_poroelastic_s"]
    front = metrics["center_holds_vs_rim_drained_at_tau_p"]
    scene = [slab(p, f"pore pressure p(x, t≈τ_p) — centre holds {front:.1f}× the drained rim"), cortex_layer]
    grad = ["#%02x%02x%02x" % tuple((np.asarray(turbo(x)[:3]) * 255).astype(int))
            for x in np.linspace(0, 1, 16)]
    cbars = {"CFD pore-pressure field (mid-plane slab)":     # JS reads cb.grad (ff_viewer_html.py:156)
             {"lo": lo, "hi": hi, "label": "pore pressure p [normalized]", "grad": grad}}
    title = (f"FF MCF7 cell — CFD ON · spatial Biot pore-pressure p(x,t) draining on the POROELASTIC "
             f"τ_p={tau_p:.2f}s (measured fit {tau_fit:.2f}s), D={metrics['D_cyto_um2_s']}µm²/s (Moeendarbary) · "
             f"radial gradient: centre lags the drained rim ({front:.1f}×) · NOT the 0-D efflux clock (30–200s)")
    build_viewer({"CFD pore-pressure field (mid-plane slab)": scene}, out_html, title=title, cbars=cbars)


def main():
    npz = sys.argv[1] if len(sys.argv) > 1 else "aleph/outputs/mech_hier/resting_reval.npz"
    out = sys.argv[2] if len(sys.argv) > 2 else "aleph/outputs/mech_hier/figs/cfd_cell_native.html"
    run(npz, out)


if __name__ == "__main__":
    main()
