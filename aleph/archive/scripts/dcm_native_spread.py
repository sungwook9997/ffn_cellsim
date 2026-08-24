"""NATIVE full-compartment single-cell DCM spreading run — authoritative G3/G4/G5.

The authoritative test for the DCM spreading mechanism (T-2 active area generation + the
energy-based conservative substrate adhesion), at NATIVE resolution with EVERY available
Warp-DCM compartment ON, from a warmed resting physiological baseline (PI 2026-07-13 rules).

Compartments (all ON): cytoplasm (turgor + K_vol), cortex (edge springs + bending), membrane
(surface tension gamma_surf), nucleus (E_nuc + lamin), substrate (rigid floor + energy adhesion).
NB the Warp DCM has NO microtubules / no `--from-resting` checkpoint (FF-engine only) — "native
full" for the DCM = subdiv 2 (162 nodes, the validated resolution) + all the above compartments,
warmed to the resting turgor/tension baseline before the spread phase.

GROUNDED params (no magic numbers; see docs/v2_audit/DCM_SPREADING_ACTIVE_AREA_PLAN + KB):
  * gamma_surf = 5e-4 N/m   — slow-band MCF7 cortical tension (KB-3.5 Moazzeni; PI-confirmed 2026-07-13)
  * w_cs       = 1e-3 J/m²  — cell-substrate adhesion (KB-2.19, grounded top of 0.1-1 mJ/m²; NOT the
                              misattributed 2.85e-3 = Gil-Redondo strain-energy x1000 slip)
  * k_area     = K_A/A0     — Brückner 2015 membrane area modulus K_A=0.13 N/m / A0 (kernel N/m³)
  * reservoir  = 1.4        — KB membrane-area reservoir (Brückner/Figard); A/A0 target ~2*res ~2.8
  * v_protrusion            — lamellipodial edge advance (KB-3.6 Mogilner-Oster); physical 1e-7 m/s
                              maps to ~min real-time; accelerated here for tractable steps (REPORTED)
  * adh_energy_range        — adhesion engagement reach; passive ~0.3-0.5 µm (active-protrusion larger)

Run on gbook A5000:
  ~/miniconda3/envs/ffn_sim/bin/python -m aleph.scripts.dcm_native_spread --device cuda:0 [flags]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from aleph.dcm.dcm_warp_decohesion import run_decohesion, _cell_volumes

R_CELL = 7.5e-6
K_A_MEMBRANE = 0.13                       # Brückner 2015 area modulus [N/m]
A0_SPHERE = 4.0 * np.pi * R_CELL ** 2
K_AREA = K_A_MEMBRANE / A0_SPHERE         # kernel N/m³ (= K_A / A0)

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "h_dcm_spreading")
os.makedirs(os.path.join(OUT, "figs"), exist_ok=True)


def _silhouette(pos):
    """Top-down xy silhouette area (all nodes) — the A/A0 convention (feedback-aa0-topdown-area)."""
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(pos[:, :2]).volume)
    except Exception:  # noqa: BLE001
        return float(np.pi * 0.25 * np.ptp(pos[:, 0]) * np.ptp(pos[:, 1]))


def _metrics(fr):
    c = fr.mean(0)
    z = fr[:, 2]
    rxy = np.linalg.norm(fr[:, :2] - c[:2], axis=1)
    p90 = np.percentile(rxy, 90)
    return dict(
        height_um=(z.max() - z.min()) * 1e6,
        sil_um2=_silhouette(fr) * 1e12,
        splay=float(_silhouette(fr) / (np.pi * p90 ** 2)) if p90 > 0 else 0.0,
        basal_frac=float((z < z.min() + 2e-6).mean()),
        rxy_max_um=float(rxy.max()) * 1e6,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--subdiv", type=int, default=2)          # native DCM resolution
    ap.add_argument("--steps", type=int, default=30000)
    ap.add_argument("--warmup", type=int, default=1500)       # settle to resting baseline (do_spread=False)
    ap.add_argument("--frames", type=int, default=30)
    ap.add_argument("--dt", type=float, default=8.0e-6)
    ap.add_argument("--gamma-surf", type=float, default=5.0e-4)
    ap.add_argument("--w-cs", type=float, default=1.0e-3)     # KB-2.19 grounded top
    ap.add_argument("--reservoir", type=float, default=1.4)   # KB Brückner
    ap.add_argument("--v-protrusion", type=float, default=1.0e-4)  # accelerated (report mapping)
    ap.add_argument("--adh-range", type=float, default=0.5e-6)
    ap.add_argument("--e-nuc", type=float, default=399.0)
    ap.add_argument("--no-nucleus", action="store_true")
    ap.add_argument("--no-bending", action="store_true")
    ap.add_argument("--no-gravity", action="store_true")
    ap.add_argument("--no-area-growth", action="store_true")
    ap.add_argument("--no-adhesion-energy", action="store_true")
    ap.add_argument("--edge-traction", action="store_true", help="T-3 edge-advancing clutch traction (the pancake flattening driver)")
    ap.add_argument("--edge-cap", type=float, default=5.0e-9, help="T-3 per-leading-node traction cap [N] (KB-2.12 per-FA ~5nN)")
    ap.add_argument("--edge-total", type=float, default=0.0, help="T-3 GRID-INVARIANT total per-cell traction [N] (0=per-node cap; Gil-Redondo MCF7 ~102nN); per-node=total/n_lead")
    ap.add_argument("--edge-advance-gap", type=float, default=1.5e-6)
    ap.add_argument("--membrane-add", action="store_true", help="option (i): reservoir cap grows over time (exocytosis) r_fold->r_deep")
    ap.add_argument("--r-deep", type=float, default=4.0, help="deep membrane reservoir ceiling (Gauthier-Masters-Sheetz 2-4x)")
    ap.add_argument("--remesh-period", type=int, default=0, help="remesh cadence [steps] (0=off; ON avoids mesh degeneracy when spreading)")
    ap.add_argument("--tag", default="native")
    args = ap.parse_args()

    npz = os.path.join(OUT, f"native_spread_{args.tag}.npz")
    print("=== NATIVE full-compartment single-cell DCM spread ===", flush=True)
    print(f"    device={args.device} subdiv={args.subdiv} ({'162' if args.subdiv==2 else '?'} nodes) "
          f"warmup={args.warmup} steps={args.steps}", flush=True)
    print(f"    gamma_surf={args.gamma_surf:.1e} w_cs={args.w_cs:.1e} k_area={K_AREA:.2e} "
          f"reservoir={args.reservoir} v_prot={args.v_protrusion:.1e} adh_range={args.adh_range:.1e}", flush=True)
    print(f"    nucleus={not args.no_nucleus} bending={not args.no_bending} gravity={not args.no_gravity} "
          f"area_growth={not args.no_area_growth} adhesion_energy={not args.no_adhesion_energy}", flush=True)

    r = run_decohesion(
        n_cells=1, subdiv=args.subdiv, device=args.device, dt=args.dt,
        steps=args.steps, warmup=args.warmup, frames=args.frames,
        k_vol=7.73e5,
        # membrane surface tension + area elasticity (T-2 setpoint lives here)
        surface_tension=True, gamma_surf=args.gamma_surf, k_area=K_AREA,
        # T-2 active area generation
        area_growth=not args.no_area_growth, area_reservoir=args.reservoir,
        v_protrusion=args.v_protrusion,
        # substrate: energy-based conservative adhesion (descent+spread) + rigid floor; split wetting OFF
        substrate_wetting=False, use_substrate_well=True,
        adhesion_energy=not args.no_adhesion_energy, adh_energy_range=args.adh_range,
        w_cs_jm2=args.w_cs,
        # T-3 edge-advancing clutch traction (the pancake flattening driver)
        edge_traction=args.edge_traction, edge_traction_cap=args.edge_cap,
        edge_total_traction=args.edge_total,
        edge_advance_gap=args.edge_advance_gap,
        # option (i): membrane addition over time + remesh (avoid mesh degeneracy)
        membrane_add=args.membrane_add, r_deep=args.r_deep,
        remesh_period=args.remesh_period,
        # nucleus + bending + gravity (physiological baseline)
        nucleus=not args.no_nucleus, E_nuc=args.e_nuc, R_nuc_factor=0.7, ratio_lamin=1.4,
        bending=not args.no_bending, k_bend=1.0e-5,
        gravity=not args.no_gravity, delta_rho=55.0,
        settle_steps=0, settle_frames=0,
        save_frames=npz,
    )

    d = np.load(npz, allow_pickle=True)
    frames = d["frames"]; faces = d["faces"]; cof = d["cof"]
    fcell = cof[faces[:, 0]]; nc = int(cof.max() + 1)
    A0 = _silhouette(frames[0]); V0 = float(np.abs(_cell_volumes(frames[0], faces, fcell, nc)).sum())
    traj = []
    for fr in frames:
        m = _metrics(fr)
        V = float(np.abs(_cell_volumes(fr, faces, fcell, nc)).sum())
        traj.append(dict(AoverA0=m["sil_um2"] * 1e-12 / A0, VoverV0=V / V0, **m))
    aa = [t["AoverA0"] for t in traj]
    tf = traj[-1]
    # real-time mapping of the accelerated protrusion (physical v0=1e-7 m/s)
    accel = args.v_protrusion / 1.0e-7
    summary = dict(
        tag=args.tag, finite=bool(r.get("finite", True)), device=args.device, subdiv=args.subdiv,
        A0_um2=A0 * 1e12, V0_um3=V0 * 1e18,
        AoverA0_final=round(tf["AoverA0"], 3), AoverA0_max=round(max(aa), 3),
        VoverV0_final=round(tf["VoverV0"], 3),
        height_final_um=round(tf["height_um"], 2), basal_frac_final=round(tf["basal_frac"], 3),
        splay_final=round(tf["splay"], 2),
        v_protrusion=args.v_protrusion, accel_factor=accel,
        params=dict(gamma_surf=args.gamma_surf, w_cs=args.w_cs, k_area=K_AREA,
                    reservoir=args.reservoir, adh_range=args.adh_range),
        # G-gate readouts
        G3_volume_conserved=bool(abs(tf["VoverV0"] - 1.0) < 0.05),
        G4_spread_2to8=bool(2.0 <= max(aa) <= 8.0),
        G5_real_flatten=bool(tf["height_um"] < 0.85 * traj[0]["height_um"] and tf["splay"] < 1.5
                             and tf["basal_frac"] > 0.25),
        trajectory=traj,
    )
    def _jsonable(o):
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))

    jpath = os.path.join(OUT, f"native_spread_{args.tag}.json")
    with open(jpath, "w") as f:
        json.dump(summary, f, indent=2, default=_jsonable)
    print(f"\n  A0={A0*1e12:.1f} µm²  V0={V0*1e18:.1f} µm³", flush=True)
    print(f"  A/A0 {aa[0]:.2f} -> {tf['AoverA0']:.2f} (max {max(aa):.2f})  V/V0={tf['VoverV0']:.3f}  "
          f"height {traj[0]['height_um']:.1f}->{tf['height_um']:.1f}µm  basal {tf['basal_frac']:.2f}  "
          f"splay {tf['splay']:.2f}", flush=True)
    print(f"  GATES: G3(V/V0=1)={summary['G3_volume_conserved']}  G4(A/A0 2-8)={summary['G4_spread_2to8']}  "
          f"G5(real flatten)={summary['G5_real_flatten']}", flush=True)
    print(f"  json -> {os.path.relpath(jpath)}  frames -> {os.path.relpath(npz)}", flush=True)
    try:
        _figure(traj, summary, args.tag)
    except Exception as e:  # noqa: BLE001
        print(f"  (figure skipped: {e})", flush=True)
    return 0


def _figure(traj, summary, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    aa = [t["AoverA0"] for t in traj]; vv = [t["VoverV0"] for t in traj]
    h = [t["height_um"] for t in traj]; bf = [t["basal_frac"] for t in traj]
    x = list(range(len(traj)))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    ax[0].axhspan(2, 8, color="#c8e6c9", alpha=0.5, label="G4 target 2-8×")
    ax[0].plot(x, aa, "o-", color="C3"); ax[0].set_ylabel("A/A0 (top-down silhouette)")
    ax[0].set_xlabel("frame"); ax[0].set_title(f"(a) spread A/A0 (final {summary['AoverA0_final']})")
    ax[0].legend(fontsize=8)
    ax[1].plot(x, vv, "s-", color="C0"); ax[1].axhline(1.0, ls=":", color="k")
    ax[1].set_ylim(0.9, 1.1); ax[1].set_ylabel("V/V0"); ax[1].set_xlabel("frame")
    ax[1].set_title(f"(b) volume conservation G3 (final {summary['VoverV0_final']})")
    ax[2].plot(x, h, "o-", color="C4", label="height µm"); ax[2].set_ylabel("height [µm]", color="C4")
    a2 = ax[2].twinx(); a2.plot(x, bf, "^--", color="C2", label="basal frac"); a2.set_ylabel("basal frac", color="C2")
    ax[2].set_xlabel("frame"); ax[2].set_title(f"(c) flatten G5 (h {h[0]:.1f}->{h[-1]:.1f}, splay {summary['splay_final']})")
    fig.suptitle(f"NATIVE full-compartment DCM single-cell spread [{tag}] — "
                 f"gamma={summary['params']['gamma_surf']:.0e} w_cs={summary['params']['w_cs']:.0e} "
                 f"reservoir={summary['params']['reservoir']} (grounded)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUT, "figs", f"native_spread_{tag}.png")
    fig.savefig(p, dpi=120); print(f"  figure -> {os.path.relpath(p)}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
