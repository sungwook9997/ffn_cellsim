"""H.7 DCM⊗ECM — a deformable cell remodeling an explicit cross-linked fiber ECM
through a catch-slip focal-adhesion clutch (the "beyond-the-papers" headline).

Builds the integration (``cell/dcm_ecm.py``): a DCM exact-turgor deformable shell
resting on an explicit cross-linked Mikado fiber ECM, gripping it via a Pereverzev
catch-slip FA clutch, and CONTRACTING (membrane edge-spring prestress). Runs it on
BAOAB and measures the matrix remodeling, validating against Slater…T.Kim Soft
Matter 2021 (d0sm01911a): inward fiber displacement, radial fiber tension that
decays ~1/r (2D), and radial alignment of tensed fibers near the cell.

Usage:
    python -m ffn_sim.scripts.h7_dcm_ecm_remodel --steps 40000 --sample-every 4000 \
        --device cpu --allow-cpu-dev \
        --out outputs/h7/dcm_ecm/remodel.json --fig outputs/h7/figs/dcm_ecm/dcm_ecm_remodel.png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--sample-every", type=int, default=4000)
    ap.add_argument("--contractility", type=float, default=0.7)
    ap.add_argument("--dt", type=float, default=5.0e-10)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--allow-cpu-dev", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--fig", default=None)
    args = ap.parse_args()

    import hoomd
    from ffn_sim.cell.manifest import load_manifest
    from ffn_sim.cell.dcm_ecm import (
        ResolvedDcmEcm, build_dcm_ecm_simulation, read_positions,
        ecm_radial_displacement, fiber_tension_vs_r, fiber_radial_alignment,
    )

    h1 = load_manifest("phase1_h1.yaml")
    p = ResolvedDcmEcm(contractility=args.contractility, dt=args.dt)
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    h = build_dcm_ecm_simulation(p, h1, device=dev)
    sim = h["sim"]
    ecm_tags = h["ecm_bead_tags"]
    ecm0 = h["ecm_pos0"]
    ecm_bonds = h["ecm_bonds"]
    ecm_bond_r0 = h["ecm_bond_r0"]
    p_h1 = h["p_h1"]
    bond_k = float(getattr(p_h1, "bond_k", getattr(p_h1, "k_bond", 1.0e-3)))
    cell_xy = h["cell_centre0"][:2]

    print(f"[dcm-ecm] BUILD: {sim.state.N_particles} particles "
          f"({h['n_ecm']} ECM beads + {h['nv']} cell nodes), "
          f"FA F*={h.get('F_star', float('nan')):.2e} N, contractility={args.contractility}",
          flush=True)

    series = []

    def sample(step):
        pos, _ = read_positions(sim)
        r0, inward = ecm_radial_displacement(pos, ecm0, ecm_tags, cell_xy)
        n_eng = int(getattr(h["fa_updater"], "n_engaged", lambda: 0)()) if callable(
            getattr(h["fa_updater"], "n_engaged", None)) else int(
            getattr(h["fa_updater"], "n_engaged", 0))
        rec = {"step": int(step),
               "mean_inward_nm": float(np.mean(inward) * 1e9),
               "max_inward_nm": float(np.max(inward) * 1e9),
               "n_fa_engaged": n_eng}
        series.append(rec)
        return rec

    r = sample(0)
    print(f"[dcm-ecm] t=0  mean_inward={r['mean_inward_nm']:.1f} nm  "
          f"FA_engaged={r['n_fa_engaged']}", flush=True)
    n_chunks = max(1, args.steps // args.sample_every)
    for c in range(1, n_chunks + 1):
        sim.run(args.sample_every)
        rec = sample(c * args.sample_every)
        print(f"[dcm-ecm] step {rec['step']:>7}  mean_inward={rec['mean_inward_nm']:.1f} nm  "
              f"max_inward={rec['max_inward_nm']:.1f} nm  FA_engaged={rec['n_fa_engaged']}",
              flush=True)

    # Final remodeling fields for the figure + the 1/r stress check.
    pos, _ = read_positions(sim)
    r0_b, inward = ecm_radial_displacement(pos, ecm0, ecm_tags, cell_xy)
    r_mid, tension = fiber_tension_vs_r(pos, ecm_bonds, ecm_bond_r0, bond_k, cell_xy)
    r_al, align = fiber_radial_alignment(pos, ecm_bonds, cell_xy)

    # 1/r fit on the tensile (stretched) fibers, in the cell-to-mid-field band.
    tens = tension > 0
    rr = r_mid[tens] * 1e6
    tt = tension[tens]
    band = (rr > p.R_cell * 1e6) & (rr < 0.9 * p.footprint_factor * p.R_cell * 1e6)
    slope = float("nan")
    if band.sum() > 5:
        lr, lt = np.log(rr[band]), np.log(tt[band] + 1e-30)
        slope = float(np.polyfit(lr, lt, 1)[0])  # ~ -1 expected (2D ~1/r)

    meta = {"n_particles": int(sim.state.N_particles), "n_ecm": int(h["n_ecm"]),
            "nv": int(h["nv"]), "contractility": args.contractility, "dt": args.dt,
            "steps": args.steps, "bond_k": bond_k,
            "tension_decay_loglog_slope": slope,
            "final_mean_inward_nm": series[-1]["mean_inward_nm"],
            "final_max_inward_nm": series[-1]["max_inward_nm"]}
    print(f"\n[dcm-ecm] DONE  mean_inward {series[0]['mean_inward_nm']:.1f}->"
          f"{series[-1]['mean_inward_nm']:.1f} nm  "
          f"tension~r^{slope:.2f} (d0sm01911a 2D predicts ~ -1)", flush=True)

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump({"meta": meta, "series": series}, fh, indent=2)
        print(f"[dcm-ecm] wrote {args.out}")

    if args.fig:
        _figure(p, h, ecm0, pos, ecm_tags, cell_xy, r0_b, inward,
                r_mid, tension, r_al, align, slope, series, args.fig)
        print(f"[dcm-ecm] wrote {args.fig}")
    return 0


def _figure(p, h, ecm0, pos, ecm_tags, cell_xy, r0_b, inward, r_mid, tension,
            r_al, align, slope, series, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5), constrained_layout=True)
    um = 1e6

    # Panel 1: ECM remodeling map — fiber beads, colored by inward displacement.
    a = ax[0, 0]
    now = pos[ecm_tags]
    sc = a.scatter(now[:, 0] * um, now[:, 1] * um, c=inward * 1e9, s=6,
                   cmap="RdBu_r", vmin=-np.abs(inward * 1e9).max(),
                   vmax=np.abs(inward * 1e9).max())
    th = np.linspace(0, 2 * np.pi, 64)
    a.plot(cell_xy[0] * um + p.R_cell * um * np.cos(th),
           cell_xy[1] * um + p.R_cell * um * np.sin(th), "k--", lw=1, label="cell")
    a.set_aspect("equal"); a.set_xlabel("x (µm)"); a.set_ylabel("y (µm)")
    a.set_title("ECM fiber bead inward displacement (remodeling)")
    plt.colorbar(sc, ax=a, label="inward (nm)"); a.legend(fontsize=8)

    # Panel 2: tension vs r (log-log) + 1/r reference.
    a = ax[0, 1]
    tens = tension > 0
    rr = r_mid[tens] * um
    tt = tension[tens] * 1e12  # pN
    a.loglog(rr, tt, ".", ms=3, color="#c0392b", alpha=0.5)
    if np.isfinite(slope) and (rr > 0).any():
        rx = np.array([rr.min(), rr.max()])
        # 1/r reference anchored at the median
        med_r = np.median(rr); med_t = np.median(tt)
        a.loglog(rx, med_t * (med_r / rx), "k--", lw=1.5, label="~1/r (d0sm01911a 2D)")
    a.set_xlabel("distance from cell r (µm)"); a.set_ylabel("fiber tension (pN)")
    a.set_title(f"radial stress decay: tension ~ r^{slope:.2f}  (2D predicts ~ -1)")
    a.legend(fontsize=8); a.grid(alpha=0.3, which="both")

    # Panel 3: radial alignment vs r (binned mean).
    a = ax[1, 0]
    rb = r_al * um
    bins = np.linspace(0, p.footprint_factor * p.R_cell * um, 12)
    idx = np.digitize(rb, bins)
    bx, by = [], []
    for i in range(1, len(bins)):
        m = idx == i
        if m.sum() > 3:
            bx.append(0.5 * (bins[i - 1] + bins[i])); by.append(np.mean(align[m]))
    a.plot(bx, by, "-o", color="#16a085", lw=2)
    a.axhline(0.5, color="0.6", ls=":", label="isotropic (0.5)")
    a.set_xlabel("distance from cell r (µm)"); a.set_ylabel("|cos(bond, radial)|")
    a.set_title("tensed fibers align radially near the cell (Kim 2021)")
    a.set_ylim(0, 1); a.legend(fontsize=8); a.grid(alpha=0.3)

    # Panel 4: mean inward displacement over time.
    a = ax[1, 1]
    t = [r["step"] for r in series]
    mi = [r["mean_inward_nm"] for r in series]
    mx = [r["max_inward_nm"] for r in series]
    a.plot(t, mi, "-o", color="#2f6fb0", lw=2, label="mean inward")
    a.plot(t, mx, "-s", color="#e67e22", lw=1.6, label="max inward")
    a.set_xlabel("BAOAB step"); a.set_ylabel("ECM inward displacement (nm)")
    a.set_title("matrix remodeling builds over time (catch-slip FA + contraction)")
    a.legend(fontsize=8); a.grid(alpha=0.3)

    fig.suptitle(
        "H.7 BEYOND-papers: DCM deformable cell remodeling an explicit cross-linked "
        "fiber ECM via catch-slip FA clutch\n"
        f"({h['n_ecm']} ECM beads + {h['nv']} cell nodes; reproduces+exceeds Slater…"
        "Kim Soft Matter 2021: catch-slip FA + turgor shell + cross-linked fibers)",
        fontsize=11, fontweight="bold")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
