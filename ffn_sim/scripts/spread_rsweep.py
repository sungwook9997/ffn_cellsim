"""R-sweep experiment: the legit test of the PI law A/A0 = a + b/R + c/R².

The PI's empirical MCF7-spheroid spreading law states the normalised top-down spread
area A/A0 of a spheroid on an adhesive substrate is a quadratic in 1/R, where R is the
spheroid radius. The ONLY honest way to probe this is to *vary R* (the controlled
independent variable) while holding EVERY derived mechanistic parameter fixed — the
cadherin bundle multiplicity, the ECM-clutch bundle multiplicity, the integrator, the
contact stiffness, etc. are NEVER tuned per-N (tuning params to the outcome is forbidden,
CLAUDE.md "no fitting to outcome").

R is controlled here through the cell count N: more cells → a larger rested aggregate →
larger R. For each N we run the full mechanistic spreading config (cadherin catch-bonds +
Pereverzev ECM clutch + lamellipodium crawl + nucleus + surface tension + bending +
edge–edge contact), settle the cube into a ball, then spread it on the dish, and read the
driver-reported A/A0 (top-down silhouette, PI rule — never basal contact). R is the
effective spheroid radius of the *rested* aggregate, R = (3·V_total/4π)^(1/3).

Then we fit A/A0 = a + b/R + c/R² (a least-squares polynomial in 1/R) and report a,b,c,r²
with a figure overlaying the fit on the per-N points.

This is a SMOKE-able experiment harness, not a production launcher: the default is a tiny
2-point quick test so the pipeline can be validated fast. A real sweep is a long GPU job
the operator launches deliberately (more N values, more steps), NOT something this script
kicks off on its own.

Run (tiny smoke):
    PYTHONPATH=. python -m ffn_sim.scripts.spread_rsweep --device cpu --ns "7,13" --steps 3000

Run (a real sweep — operator-launched, long):
    PYTHONPATH=. python -m ffn_sim.scripts.spread_rsweep --device cuda:0 \
        --ns "13,24,48,100" --steps 40000
"""

from __future__ import annotations

import argparse
import os
import tempfile

import numpy as np

from ffn_sim.warp_port.dcm_warp_decohesion import run_decohesion

# Derived mechanistic parameters — FIXED across every N (never tuned per-N).
CAD_BUNDLE = 40.0          # cadherin catch-bond bundle multiplicity (derived)
ECM_BUNDLE = 167.0         # Pereverzev ECM-clutch bundle multiplicity (derived)
ACCEL_DT = 8.0e-4          # implicit IMEX dt
GAMMA_SURF = 1.0e-4        # cortical surface tension [N/m]

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "warp_decohesion", "figs")


def effective_radius_um(vol_sum_m3: float) -> float:
    """Effective spheroid radius R [µm] from the total aggregate volume.

    R = (3·V_total / 4π)^(1/3): the radius of the equivalent solid sphere of the same
    volume as the summed per-cell volumes of the rested aggregate.

    Args:
        vol_sum_m3: total aggregate volume ∑ V_cell [m³].

    Returns:
        Effective spheroid radius R [µm].
    """
    R_m = (3.0 * max(vol_sum_m3, 0.0) / (4.0 * np.pi)) ** (1.0 / 3.0)
    return float(R_m * 1e6)


def gyration_radius_um(P: np.ndarray, cof: np.ndarray) -> float:
    """Gyration radius R_g [µm] of the live nodes (alternative R measure / cross-check).

    Args:
        P: (N,3) node positions [m].
        cof: (N,) per-node cell index; live nodes are cof>=0.

    Returns:
        Radius of gyration R_g [µm] of the live aggregate.
    """
    live = P[cof >= 0]
    if live.shape[0] < 2:
        return 0.0
    rel = live - live.mean(0)
    return float(np.sqrt((rel ** 2).sum(axis=1).mean()) * 1e6)


def run_one(n_cells: int, steps: int, device: str) -> dict:
    """Run one spreading point at cell count ``n_cells`` and extract (R, A/A0).

    The full mechanistic spreading config is used with the derived parameters held FIXED
    (see module constants); only ``n_cells`` and ``steps`` vary. A settle phase rounds the
    cube into a ball, then an equal-length spread phase wets it on the dish.

    Args:
        n_cells: number of cells (the controlled R lever).
        steps: spread-phase step count (settle uses the same count).
        device: Warp device string, e.g. "cpu" or "cuda:0".

    Returns:
        Dict with n_cells, R_um (volume-equiv), Rg_um (gyration), aa0_final, vv0_final.
    """
    tmp = os.path.join(tempfile.gettempdir(), f"_rsweep_n{n_cells}.npz")
    out = run_decohesion(
        n_cells=n_cells, subdiv=2, steps=steps, frames=4, device=device,
        dt=8e-6, warmup=1000, settle_steps=steps, settle_frames=2, gap=2.05,
        cadherin=True, cad_bundle=CAD_BUNDLE,
        ecm_clutch=True, ecm_bundle=ECM_BUNDLE,
        lamellipodium=True, nucleus=True,
        surface_tension=True, gamma_surf=GAMMA_SURF,
        bending=True, edge_edge=True,
        substrate_wetting=True, use_substrate_well=True,
        builder="fcc", integrator="implicit", accel_dt=ACCEL_DT,
        save_frames=tmp)
    if "error" in out:
        raise RuntimeError(f"N={n_cells} run failed: {out['error']}")
    d = np.load(tmp, allow_pickle=True)
    P = d["frames"][-1].astype(np.float64)
    cof = d["cof"]
    # rested-aggregate volume: vv0_final · V0_total. We recover V_total from vv0 and the
    # per-cell rested volume is not directly returned, so use the gyration radius of the
    # settled ball as the volume proxy when needed; but the volume-equivalent R is the
    # primary measure, reconstructed from the summed live-cell volumes of the final frame.
    R_vol = _volume_radius_um(P, d["faces"].astype(np.int64), cof)
    return {
        "n_cells": n_cells,
        "R_um": R_vol,
        "Rg_um": gyration_radius_um(P, cof),
        "aa0_final": float(out["aa0_final"]),
        "aa0_peak": float(out["aa0_peak"]),
        "vv0_final": float(out["vv0_final"]),
    }


def _volume_radius_um(P: np.ndarray, faces: np.ndarray, cof: np.ndarray) -> float:
    """Effective radius R [µm] from the summed signed per-face tetra volume of live cells.

    R = (3·V_total/4π)^(1/3). Faces index into the FULL node array; dormant/parked nodes
    (cof<0) are far away but their faces still close on themselves, so we mask to faces
    whose vertices are all live.

    Args:
        P: (N,3) node positions [m].
        faces: (F,3) triangle vertex indices.
        cof: (N,) per-node cell index; live nodes are cof>=0.

    Returns:
        Volume-equivalent spheroid radius R [µm].
    """
    live = cof >= 0
    keep = live[faces[:, 0]] & live[faces[:, 1]] & live[faces[:, 2]]
    f = faces[keep]
    if f.shape[0] == 0:
        return 0.0
    v0, v1, v2 = P[f[:, 0]], P[f[:, 1]], P[f[:, 2]]
    vol = abs(float(np.einsum("ij,ij->i", v0, np.cross(v1 - v0, v2 - v0)).sum() / 6.0))
    return effective_radius_um(vol)


def fit_law(R: np.ndarray, aa0: np.ndarray) -> tuple[float, float, float, float]:
    """Fit A/A0 = a + b/R + c/R² (least squares in 1/R) and return (a, b, c, r²).

    Args:
        R: (k,) effective radii [µm].
        aa0: (k,) measured A/A0 at each radius.

    Returns:
        (a, b, c, r²). With fewer than 3 points the fit degrades to the available order
        (c=0 for 2 points, b=c=0 for 1 point) and r² is set to nan when undefined.
    """
    x = 1.0 / np.asarray(R, float)
    y = np.asarray(aa0, float)
    order = min(2, x.size - 1)                       # 2 pts → linear, 1 pt → constant
    coeffs = np.polyfit(x, y, order)                 # highest power first
    # pad to [c, b, a] (c is the 1/R² coefficient)
    full = np.concatenate([np.zeros(3 - coeffs.size), coeffs])
    c, b, a = full[0], full[1], full[2]
    yhat = c * x ** 2 + b * x + a
    ss_res = float(((y - yhat) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = float("nan") if ss_tot <= 0 else 1.0 - ss_res / ss_tot
    return float(a), float(b), float(c), r2


def save_figure(rows: list[dict], a: float, b: float, c: float, r2: float, path: str) -> None:
    """Save the A/A0-vs-R figure: per-N points + the fitted a + b/R + c/R² curve.

    Visualization-integrity rules (CLAUDE.md): no axis truncation, units annotated, the
    fit overlaid on the measurements.

    Args:
        rows: per-N result dicts (need R_um, aa0_final, n_cells).
        a, b, c, r2: fitted law coefficients and r².
        path: output PNG path.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    R = np.array([r["R_um"] for r in rows])
    aa0 = np.array([r["aa0_final"] for r in rows])
    ns = [r["n_cells"] for r in rows]

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.scatter(R, aa0, s=60, color="C0", zorder=3, label="per-N runs")
    for ri, ai, ni in zip(R, aa0, ns):
        ax.annotate(f"N={ni}", (ri, ai), textcoords="offset points", xytext=(6, 4), fontsize=8)
    if R.size >= 2:
        rr = np.linspace(R.min() * 0.95, R.max() * 1.05, 200)
        yy = a + b / rr + c / rr ** 2
        ax.plot(rr, yy, "C3-", lw=1.6,
                label=f"fit a+b/R+c/R²  (a={a:.3g}, b={b:.3g}, c={c:.3g}, r²={r2:.3f})")
    ax.set_xlabel("effective spheroid radius R  [µm]")
    ax.set_ylabel("spread area  A/A0  (top-down silhouette)")
    ax.set_title("R-sweep: PI law A/A0 = a + b/R + c/R²")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(f"  saved figure -> {os.path.abspath(path)}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--device", default="cuda:0",
                    help="Warp device (cuda:0 default; use cpu for a dev smoke)")
    ap.add_argument("--ns", default="7,13",
                    help="comma-separated cell counts N (the R lever); SMALL default for smoke")
    ap.add_argument("--steps", type=int, default=3000,
                    help="spread-phase steps (settle uses the same); SMALL default for smoke")
    ap.add_argument("--fig", default=os.path.join(FIG_DIR, "rsweep_aa0_vs_R.png"),
                    help="output figure path")
    args = ap.parse_args()

    ns = [int(x) for x in args.ns.split(",") if x.strip()]
    print(f"R-sweep — N={ns}, steps={args.steps}, device={args.device}")
    print(f"  FIXED derived params: cad_bundle={CAD_BUNDLE:.0f}, ecm_bundle={ECM_BUNDLE:.0f}, "
          f"accel_dt={ACCEL_DT:.0e}, integrator=implicit  (NEVER tuned per-N)\n")
    print(f"{'N':>5} {'R[µm]':>9} {'Rg[µm]':>9} {'A/A0_final':>11} {'A/A0_peak':>10} {'V/V0':>7}")
    rows = []
    for n in ns:
        r = run_one(n, args.steps, args.device)
        rows.append(r)
        print(f"{r['n_cells']:5d} {r['R_um']:9.3f} {r['Rg_um']:9.3f} "
              f"{r['aa0_final']:11.4f} {r['aa0_peak']:10.4f} {r['vv0_final']:7.3f}", flush=True)

    R = np.array([r["R_um"] for r in rows])
    aa0 = np.array([r["aa0_final"] for r in rows])
    a, b, c, r2 = fit_law(R, aa0)
    print(f"\n[FIT] A/A0 = a + b/R + c/R²")
    print(f"      a = {a:.5g}")
    print(f"      b = {b:.5g}  [µm]")
    print(f"      c = {c:.5g}  [µm²]")
    print(f"      r² = {r2:.4f}  ({R.size} points)")
    if R.size < 3:
        print("      (note: <3 points → reduced-order fit, c may be 0; this is a smoke run)")
    save_figure(rows, a, b, c, r2, args.fig)


if __name__ == "__main__":
    main()
