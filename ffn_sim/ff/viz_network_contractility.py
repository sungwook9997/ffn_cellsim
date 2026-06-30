"""Figure: contractile-network σ vs per-link force and vs connectivity z (Stage 6h decisive test).

Shows the active contractile stress σ stays deep under the cortical band across both the contractile
force f_act AND the network connectivity z (incl. the sub-isostatic regime) — the floor survives in
the buckling-capable inextensible network. Writes outputs/ff/figs/network_contractility.png.
Run: python -m ffn_sim.ff.viz_network_contractility
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.ff.gamma_estimator import SALBREUX_BAND_PN_UM
from ffn_sim.ff.network_contractility import build_patch, measure_sigma

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "ff")
NMIIA_STALL = 5.0  # pN, per-side NMIIA minifilament stall (the physiological anchor)


def run_sweeps(n_steps=6000):
    f_vals = [1.0, 3.0, 8.0, 20.0, 50.0, 100.0]
    sig_f = []
    for f in f_vals:
        p = build_patch(nx=9, ny=9, keep_frac=1.0, rng=np.random.default_rng(0))
        o = measure_sigma(p, f_act=f, n_steps=n_steps)
        sig_f.append(o["sigma_pN_um"] if o else np.nan)
        print(f"  f_act={f:5.0f} → σ={sig_f[-1]:.2f} pN/µm", flush=True)
    keep = [0.55, 0.65, 0.75, 0.85, 1.0]
    z_vals, sig_z = [], []
    for kf in keep:
        rs = []
        zs = []
        for s in range(3):
            p = build_patch(nx=11, ny=11, keep_frac=kf, rng=np.random.default_rng(s))
            if p is None:
                continue
            o = measure_sigma(p, f_act=20.0, n_steps=n_steps)
            if o:
                rs.append(o["sigma_pN_um"]); zs.append(o["z_mean"])
        if rs:
            z_vals.append(np.mean(zs)); sig_z.append(np.mean(rs))
            print(f"  keep={kf} z={z_vals[-1]:.2f} → σ={sig_z[-1]:.2f} pN/µm", flush=True)
    return f_vals, sig_f, z_vals, sig_z


def render(f_vals, sig_f, z_vals, sig_z, outdir=OUTDIR):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.join(outdir, "figs"), exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax in (ax1, ax2):
        ax.axhspan(*SALBREUX_BAND_PN_UM, color="tab:green", alpha=0.15,
                   label="Salbreux band 350–650 pN/µm")
        ax.set_yscale("log"); ax.set_ylabel("contractile σ [pN/µm]  (band = 0.35–0.65 mN/m)")
        ax.grid(True, which="both", alpha=0.25)
    ax1.plot(f_vals, sig_f, "-o", color="tab:red", label="σ (connected network)")
    ax1.axvline(NMIIA_STALL, color="k", ls=":", label=f"NMIIA stall {NMIIA_STALL:.0f} pN")
    ax1.set_xlabel("per-link contractile force f_act [pN]"); ax1.set_xscale("log")
    ax1.set_title("σ vs contractility (z≈6 full lattice)")
    ax1.legend(fontsize=8)
    ax2.plot(z_vals, sig_z, "-s", color="tab:purple", label="σ at f_act=20 pN")
    ax2.axvline(4.0, color="gray", ls="--", label="2D isostatic z=4")
    ax2.set_xlabel("network connectivity z (diluted lattice)")
    ax2.set_title("σ vs connectivity (sub-isostatic → over-isostatic)")
    ax2.legend(fontsize=8)
    fig.suptitle("Contractile inextensible network: σ stays ~30–100× UNDER band across f_act AND z\n"
                 "→ buckling/connectivity does not lift the active-γ floor (FF/Cytosim inextensible regime)")
    fig.tight_layout()
    path = os.path.join(outdir, "figs", "network_contractility.png")
    fig.savefig(path, dpi=130); plt.close(fig)
    print(f"wrote {path}")
    return path


def main():
    print("contractile-network σ sweeps (f_act, z)...", flush=True)
    f_vals, sig_f, z_vals, sig_z = run_sweeps()
    np.savez(os.path.join(OUTDIR, "network_contractility.npz"),
             f_vals=f_vals, sig_f=sig_f, z_vals=z_vals, sig_z=sig_z, band=SALBREUX_BAND_PN_UM)
    render(f_vals, sig_f, z_vals, sig_z)


if __name__ == "__main__":
    main()
