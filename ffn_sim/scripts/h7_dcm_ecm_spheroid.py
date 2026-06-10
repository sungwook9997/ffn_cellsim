"""H.7 DCM⊗ECM SPHEROID — a CLUSTER of deformable cells collectively remodeling
an explicit cross-linked fiber ECM (the PI's ultimate goal: spheroid spreading
on ECM).

Multicell generalisation of ``h7_dcm_ecm_remodel.py``. Places a small cluster of
DCM cells (default 7, 2D hex rosette) on the explicit cross-linked Mikado fiber
bed (``cell/dcm_ecm.build_dcm_ecm_spheroid_simulation``), each cell:
  * exact-turgor deformable shell with contractile edge springs,
  * Pereverzev catch-slip FA clutch gripping the bed (its basal nodes → nearest
    fiber beads, force-free at build),
  * cadherin-scale cell-cell adhesion to neighbours + WCA excluded volume,
runs it on BAOAB, and measures COLLECTIVE matrix remodeling (inward fiber
displacement toward the cluster centroid) and the CLUSTER FOOTPRINT over time.

Uses the best single-cell sweep config (contractility 0.45, turgor 133 Pa,
dt 5e-10) and widens the bed footprint so it spans the whole cluster. Conservative
stability: small dt, gapped start, per-bond force-free bonds.

Usage:
    python -m ffn_sim.scripts.h7_dcm_ecm_spheroid --n-cells 7 --steps 30000 \
        --sample-every 5000 --contractility 0.45 --turgor 133 --dt 5e-10 \
        --device cpu --allow-cpu-dev \
        --out ffn_sim/outputs/h7/dcm_ecm/spheroid.json \
        --fig ffn_sim/outputs/h7/figs/dcm_ecm/dcm_ecm_spheroid.png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _cluster_footprint(mem_xy: np.ndarray, centroid_xy: np.ndarray) -> dict:
    """Cluster footprint observables from all membrane-node xy positions.

    Returns the convex-hull area (the spread footprint), the radius of gyration
    of the membrane nodes about the cluster centroid, and the mean node–centroid
    distance — three complementary footprint measures.
    """
    d = mem_xy - centroid_xy[None, :]
    rg = float(np.sqrt(np.mean(np.einsum("ij,ij->i", d, d))))
    mean_r = float(np.mean(np.linalg.norm(d, axis=1)))
    # convex-hull area (footprint); fall back to bbox if scipy missing / degenerate
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(mem_xy)
        area = float(hull.volume)  # 2D: .volume is the area
    except Exception:
        area = float((mem_xy[:, 0].max() - mem_xy[:, 0].min())
                     * (mem_xy[:, 1].max() - mem_xy[:, 1].min()))
    return {"hull_area_um2": area * 1e12, "rg_um": rg * 1e6,
            "mean_node_r_um": mean_r * 1e6}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=7,
                    help="cells in the cluster (7 = hex rosette)")
    ap.add_argument("--cluster-mode", default="2d", choices=["2d", "3d"],
                    help="2d hex monolayer (recommended on a flat bed) | 3d FCC ball")
    ap.add_argument("--spacing-factor", type=float, default=2.3,
                    help="cell-centre spacing = factor·R_cell (gapped start)")
    ap.add_argument("--steps", type=int, default=30000)
    ap.add_argument("--sample-every", type=int, default=5000)
    ap.add_argument("--contractility", type=float, default=0.45,
                    help="edge r0 = contractility·mean_edge (<1 contracts); "
                         "0.45 = best single-cell sweep config")
    ap.add_argument("--turgor", type=float, default=133.0,
                    help="osmotic turgor dP0 (Pa); 133 = physiological MCF7")
    ap.add_argument("--footprint-factor", type=float, default=None,
                    help="bed half-window = factor·R_cell; default auto-fit to "
                         "span the cluster + a cell radius")
    ap.add_argument("--W-cc", type=float, default=0.2e-3,
                    help="cell-cell adhesion energy density (J/m², cadherin scale)")
    ap.add_argument("--k-fa", type=float, default=None)
    ap.add_argument("--fa-capture", type=float, default=None)
    ap.add_argument("--dt", type=float, default=5.0e-10)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--allow-cpu-dev", action="store_true")
    ap.add_argument("--out", default="ffn_sim/outputs/h7/dcm_ecm/spheroid.json")
    ap.add_argument("--fig", default="ffn_sim/outputs/h7/figs/dcm_ecm/dcm_ecm_spheroid.png")
    args = ap.parse_args()

    import hoomd
    from ffn_sim.cell.manifest import load_manifest
    from ffn_sim.cell.dcm import _cluster_centers
    from ffn_sim.cell.dcm_ecm import (
        ResolvedDcmEcm, build_dcm_ecm_spheroid_simulation,
        read_positions, ecm_radial_displacement, fiber_tension_vs_r,
    )

    h1 = load_manifest("phase1_h1.yaml")

    # Auto-fit the bed footprint to span the cluster + a cell radius (+margin) so
    # every cell sits over fibers (else the outer cells overhang the bed edge).
    # (slots=True dataclass: read the default off a throwaway instance, not the class.)
    R_cell = ResolvedDcmEcm().R_cell
    if args.footprint_factor is None:
        centres = _cluster_centers(args.n_cells, args.spacing_factor * R_cell,
                                   0.0, R_cell, mode=args.cluster_mode)
        cluster_halfspan = float(np.abs(centres[:, :2]).max())
        # window must reach the outer cell rim + a little margin
        need_W = cluster_halfspan + R_cell * 1.15
        footprint_factor = need_W / R_cell
    else:
        footprint_factor = args.footprint_factor

    kw = dict(contractility=args.contractility, turgor_dP0=args.turgor,
              dt=args.dt, footprint_factor=footprint_factor)
    if args.k_fa is not None:
        kw["k_fa"] = args.k_fa
    if args.fa_capture is not None:
        kw["fa_capture"] = args.fa_capture
    p = ResolvedDcmEcm(**kw)

    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    h = build_dcm_ecm_spheroid_simulation(
        p, h1, args.n_cells, cluster_mode=args.cluster_mode,
        spacing_factor=args.spacing_factor, W_cc_Jm2=args.W_cc, device=dev)
    sim = h["sim"]
    ecm_tags = h["ecm_bead_tags"]
    ecm0 = h["ecm_pos0"]
    ecm_bonds = h["ecm_bonds"]
    ecm_bond_r0 = h["ecm_bond_r0"]
    ecm_bond_k = h["ecm_bond_k"]
    mem_tags = h["mem_tags"]
    cell_centres0 = h["cell_centres0"]
    centroid_xy = h["cluster_centroid0_xy"]   # cluster centroid (inward reference)

    print(f"[spheroid] BUILD: {sim.state.N_particles} particles "
          f"({h['n_ecm']} ECM beads + {h['n_cells']}×{h['nv']} cell nodes), "
          f"footprint_factor={footprint_factor:.2f} "
          f"(W={footprint_factor*R_cell*1e6:.1f} µm), "
          f"FA F*={h['F_star']:.2e} N, FA engaged={h['fa_force'].pairs.shape[0]}, "
          f"contractility={args.contractility} turgor={args.turgor} Pa", flush=True)

    def cluster_mean_radius():
        """Mean per-cell membrane radius (cell-shell size; contraction proxy)."""
        pos, _ = read_positions(sim)
        mp = pos[mem_tags]
        cot = h["cell_of_memtag"]
        rr = []
        for c in range(h["n_cells"]):
            m = cot == c
            rr.append(np.linalg.norm(mp[m] - cell_centres0[c][None, :], axis=1).mean())
        return float(np.mean(rr))

    r_cell0 = cluster_mean_radius()
    series = []
    stable = True

    def sample(step):
        nonlocal stable
        pos, _ = read_positions(sim)
        if not np.all(np.isfinite(pos)) or np.abs(pos).max() > 10.0 * sim.state.box.Lx:
            stable = False
        # collective inward fiber displacement toward the CLUSTER CENTROID
        r0, inward = ecm_radial_displacement(pos, ecm0, ecm_tags, centroid_xy)
        mem_xy = pos[mem_tags][:, :2]
        fp = _cluster_footprint(mem_xy, centroid_xy)
        rec = {"step": int(step),
               "mean_inward_nm": float(np.mean(inward) * 1e9),
               "max_inward_nm": float(np.max(inward) * 1e9),
               "cluster_mean_cell_radius_nm": cluster_mean_radius() * 1e9,
               "n_fa_engaged": int(h["fa_force"].pairs.shape[0]),
               **fp}
        series.append(rec)
        return rec

    r = sample(0)
    print(f"[spheroid] t=0  mean_inward={r['mean_inward_nm']:.2f} nm  "
          f"hull_area={r['hull_area_um2']:.1f} µm²  Rg={r['rg_um']:.2f} µm  "
          f"FA={r['n_fa_engaged']}", flush=True)

    n_chunks = max(1, args.steps // args.sample_every)
    for c in range(1, n_chunks + 1):
        sim.run(args.sample_every)
        rec = sample(c * args.sample_every)
        print(f"[spheroid] step {rec['step']:>7}  mean_inward={rec['mean_inward_nm']:6.2f} nm  "
              f"max_inward={rec['max_inward_nm']:7.1f} nm  "
              f"hull_area={rec['hull_area_um2']:7.1f} µm²  Rg={rec['rg_um']:.2f} µm  "
              f"FA={rec['n_fa_engaged']}  stable={stable}", flush=True)
        if not stable:
            print("[spheroid] BLOWUP detected — stopping early", flush=True)
            break

    # Final remodeling fields for the figure + the 1/r tension check (raw bed
    # tension toward the cluster centroid; no per-cell baseline twin for the
    # multicell run — the headline observable here is COLLECTIVE remodeling +
    # footprint, the single-cell driver owns the cell-induced-tension validation).
    pos, _ = read_positions(sim)
    r0_b, inward = ecm_radial_displacement(pos, ecm0, ecm_tags, centroid_xy)
    pos_cell_ecm = pos[ecm_tags]
    r_mid, tension = fiber_tension_vs_r(
        pos_cell_ecm, ecm_bonds, ecm_bond_r0, ecm_bond_k, centroid_xy)

    tens = tension > 0
    rr = r_mid[tens] * 1e6
    tt = tension[tens]
    cluster_R = float(np.linalg.norm(cell_centres0[:, :2] - centroid_xy, axis=1).max()
                      + p.R_cell)
    band = (rr > cluster_R * 1e6) & (rr < 0.95 * footprint_factor * p.R_cell * 1e6)
    slope = float("nan")
    if band.sum() > 5:
        lr, lt = np.log(rr[band]), np.log(tt[band] + 1e-30)
        slope = float(np.polyfit(lr, lt, 1)[0])

    r_cell_final = cluster_mean_radius()
    cell_radius_change_pct = 100.0 * (r_cell_final - r_cell0) / r_cell0
    area0 = series[0]["hull_area_um2"]
    area_final = series[-1]["hull_area_um2"]
    footprint_change_pct = 100.0 * (area_final - area0) / area0

    meta = {
        "n_particles": int(sim.state.N_particles), "n_ecm": int(h["n_ecm"]),
        "nv": int(h["nv"]), "n_cells": int(h["n_cells"]),
        "cluster_mode": args.cluster_mode, "spacing_factor": args.spacing_factor,
        "footprint_factor": float(footprint_factor),
        "contractility": args.contractility, "turgor_dP0": args.turgor,
        "W_cc_Jm2": args.W_cc, "k_fa": float(p.k_fa),
        "fa_capture": float(p.fa_capture), "dt": args.dt, "steps": args.steps,
        "stable": bool(stable),
        "F_star_N": float(h["F_star"]),
        "tension_decay_loglog_slope": slope,
        "cluster_mean_cell_radius_nm_0": r_cell0 * 1e9,
        "cluster_mean_cell_radius_nm_final": r_cell_final * 1e9,
        "cell_radius_change_pct": cell_radius_change_pct,
        "hull_area_um2_0": area0, "hull_area_um2_final": area_final,
        "footprint_change_pct": footprint_change_pct,
        "n_fa_engaged_0": series[0]["n_fa_engaged"],
        "n_fa_engaged_final": int(h["fa_force"].pairs.shape[0]),
        "final_mean_inward_nm": series[-1]["mean_inward_nm"],
        "final_max_inward_nm": series[-1]["max_inward_nm"],
    }
    print(f"\n[spheroid] DONE  stable={stable}  "
          f"mean_inward {series[0]['mean_inward_nm']:.2f}->"
          f"{series[-1]['mean_inward_nm']:.2f} nm  "
          f"hull_area {area0:.1f}->{area_final:.1f} µm² ({footprint_change_pct:+.2f}%)  "
          f"cell_r {r_cell0*1e9:.1f}->{r_cell_final*1e9:.1f} nm "
          f"({cell_radius_change_pct:+.2f}%)  "
          f"FA={h['fa_force'].pairs.shape[0]}  tension~r^{slope:.2f}", flush=True)

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump({"meta": meta, "series": series}, fh, indent=2)
        print(f"[spheroid] wrote {args.out}")

    if args.fig:
        _figure(p, h, ecm0, pos, ecm_tags, mem_tags, centroid_xy, cell_centres0,
                r0_b, inward, r_mid, tension, slope, series, footprint_factor,
                args.fig)
        print(f"[spheroid] wrote {args.fig}")
    return 0


def _figure(p, h, ecm0, pos, ecm_tags, mem_tags, centroid_xy, cell_centres0,
            r0_b, inward, r_mid, tension, slope, series, footprint_factor,
            out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(13.5, 10), constrained_layout=True)
    um = 1e6

    # Panel 1: collective remodeling map — fibers colored by inward displacement,
    # cell membrane nodes overlaid (the cluster footprint on the bed).
    a = ax[0, 0]
    now = pos[ecm_tags]
    vmax = np.abs(inward * 1e9).max() or 1.0
    sc = a.scatter(now[:, 0] * um, now[:, 1] * um, c=inward * 1e9, s=5,
                   cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    mem_now = pos[mem_tags]
    a.scatter(mem_now[:, 0] * um, mem_now[:, 1] * um, s=3, c="0.15",
              alpha=0.55, label="cell nodes")
    th = np.linspace(0, 2 * np.pi, 48)
    for c in range(h["n_cells"]):
        a.plot(cell_centres0[c, 0] * um + p.R_cell * um * np.cos(th),
               cell_centres0[c, 1] * um + p.R_cell * um * np.sin(th),
               "k--", lw=0.7, alpha=0.6)
    a.set_aspect("equal"); a.set_xlabel("x (µm)"); a.set_ylabel("y (µm)")
    a.set_title(f"collective ECM remodeling ({h['n_cells']}-cell cluster)")
    plt.colorbar(sc, ax=a, label="inward toward cluster (nm)")
    a.legend(fontsize=8, loc="upper right")

    # Panel 2: tension vs r (log-log) toward the cluster centroid + 1/r ref.
    a = ax[0, 1]
    tens = tension > 0
    rr = r_mid[tens] * um
    tt = tension[tens] * 1e12  # pN
    a.loglog(rr, tt, ".", ms=3, color="#c0392b", alpha=0.45)
    if np.isfinite(slope) and rr.size and (rr > 0).any():
        rx = np.array([rr.min(), rr.max()])
        med_r = np.median(rr); med_t = np.median(tt)
        a.loglog(rx, med_t * (med_r / rx), "k--", lw=1.5, label="~1/r (d0sm01911a 2D)")
    a.set_xlabel("distance from cluster centroid r (µm)")
    a.set_ylabel("fiber tension (pN)")
    a.set_title(f"radial stress decay: tension ~ r^{slope:.2f}  (2D predicts ~ -1)")
    a.legend(fontsize=8); a.grid(alpha=0.3, which="both")

    # Panel 3: cluster footprint (convex-hull area) + Rg over time.
    a = ax[1, 0]
    t = [r["step"] for r in series]
    area = [r["hull_area_um2"] for r in series]
    a.plot(t, area, "-o", color="#8e44ad", lw=2, label="convex-hull footprint")
    a.set_xlabel("BAOAB step"); a.set_ylabel("cluster footprint area (µm²)",
                                              color="#8e44ad")
    a.tick_params(axis="y", labelcolor="#8e44ad")
    a2 = a.twinx()
    rg = [r["rg_um"] for r in series]
    a2.plot(t, rg, "-s", color="#16a085", lw=1.6, label="R_g (membrane)")
    a2.set_ylabel("membrane radius of gyration (µm)", color="#16a085")
    a2.tick_params(axis="y", labelcolor="#16a085")
    a.set_title("cluster footprint over time (collective spread/compaction)")
    a.grid(alpha=0.3)

    # Panel 4: collective inward fiber displacement + FA engagement over time.
    a = ax[1, 1]
    mi = [r["mean_inward_nm"] for r in series]
    mx = [r["max_inward_nm"] for r in series]
    a.plot(t, mi, "-o", color="#2f6fb0", lw=2, label="mean inward")
    a.plot(t, mx, "-s", color="#e67e22", lw=1.6, label="max inward")
    a.set_xlabel("BAOAB step"); a.set_ylabel("ECM inward displacement (nm)")
    a3 = a.twinx()
    fa = [r["n_fa_engaged"] for r in series]
    a3.plot(t, fa, "-^", color="#7f8c8d", lw=1.3, label="FA engaged")
    a3.set_ylabel("FA clutches engaged", color="#7f8c8d")
    a3.tick_params(axis="y", labelcolor="#7f8c8d")
    a.set_title("collective matrix remodeling (catch-slip FA + contraction)")
    a.legend(fontsize=8, loc="upper left"); a.grid(alpha=0.3)

    fig.suptitle(
        "H.7 SPHEROID: a CLUSTER of DCM deformable cells collectively remodeling "
        "an explicit cross-linked fiber ECM via catch-slip FA clutches\n"
        f"({h['n_ecm']} ECM beads + {h['n_cells']}×{h['nv']} cell nodes; multicell "
        "extension of Slater…Kim Soft Matter 2021 — the PI's spheroid-on-ECM goal)",
        fontsize=11, fontweight="bold")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
