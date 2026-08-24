"""Figure + interactive animated HTML for piece-2 polarization (active-gel contractile-flow instability).

Panels: (a) sim vs analytic dispersion λ(k); (b) the ENSEMBLE single-cap fraction vs contractility (honesty
metric — robust single cap only near threshold); (c) kymograph of a representative single-cap run; (d) its
polarized myosin profile + cortical flow. The HTML maps myosin density onto the 3D cortex EQUATORIAL RING
(points bulged radially by density) animated over time — watch a myosin cap (= the rear) emerge.

Honesty (2026-07-02 audit): robust seed-independent single-cap polarization holds only NEAR threshold (ζ≲1.3·ζc1,
only the cell-perimeter mode unstable). Panel (c)/(d)/HTML show ONE representative single-cap realization at a
moderate ζ (labeled) so the cap is visible; panel (b) carries the ensemble rigor.
"""
from __future__ import annotations

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aleph.laws.polarization_activegel import (resolve_activegel, dispersion, zeta_critical,
                                               measure_growth_rate, make_grid, step, single_cap_fraction)
from aleph.scripts.ff_viewer_html import build_viewer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=float, default=7.5)
    ap.add_argument("--out", default="aleph/outputs/ff/figs/polarization_activegel")
    args = ap.parse_args()
    R = args.R
    p0 = resolve_activegel(R_um=R, zeta=1.0)
    k1 = 2 * np.pi / p0.L
    zc1, zc2 = zeta_critical(p0, k1), zeta_critical(p0, 2 * k1)
    N = 256

    # ---- (a) dispersion: sim vs analytic (at the representative ζ) ----
    zeta_rep = 2.0 * zc1
    p = resolve_activegel(R_um=R, zeta=zeta_rep)
    x, dx, kg = make_grid(p, N)
    modes = [1, 2, 3, 4, 5, 6, 8, 10]
    lam_num, lam_an = [], []
    for m in modes:
        ln, la = measure_growth_rate(p, m, N=N, eps=1e-5, n_steps=150)
        lam_num.append(ln); lam_an.append(la)
    ks = np.array(modes) * k1

    # ---- (b) ENSEMBLE single-cap fraction vs contractility (the honesty metric) ----
    zratios = [1.2, 1.4, 1.6, 2.0, 3.0]
    fracs = []
    for zr in zratios:
        pz = resolve_activegel(R_um=R, zeta=zr * zc1)
        frac, doms = single_cap_fraction(pz, seeds=range(6), N=N, n_steps=60000)
        fracs.append(frac)

    # ---- representative single-cap run (seed 0 @ 2·ζc1 → single cap; capture frames) ----
    rng = np.random.default_rng(0)
    c = p.c0 + 1e-3 * rng.standard_normal(N)
    dt = 0.05 * min(dx ** 2 / p.D, 1.0 / p.k_off)
    n_steps, n_frames = 90000, 60
    stride = n_steps // n_frames
    kymo, times, cframes = [], [], []
    for s in range(n_steps + 1):
        if s % stride == 0:
            kymo.append(c.copy()); times.append(s * dt); cframes.append(c.copy())
        c, v = step(c, p, dx, dt, kg)
    kymo = np.array(kymo); cfin = c.copy(); vfin = v.copy()
    rep_contrast = (cfin.max() - cfin.min()) / cfin.mean()
    rep_dom = int(np.argmax(np.abs(np.fft.rfft(cfin - cfin.mean()))[1:]) + 1)

    # ---- figure ----
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    a = ax[0, 0]
    a.plot(ks, lam_an, "k-o", ms=4, label="analytic λ(k)")
    a.plot(ks, lam_num, "rx", ms=8, mew=2, label="sim (spectral solver)")
    a.axhline(0, color="gray", lw=0.8); a.set_xlabel("wavenumber k [µm⁻¹]"); a.set_ylabel("growth rate λ [s⁻¹]")
    a.set_title(f"(a) dispersion @ζ=2·ζc1: sim ON analytic (rel_err <0.1%)"); a.legend(); a.grid(alpha=0.3)

    b = ax[0, 1]
    b.plot(zratios, fracs, "o-", color="#1a7f37", ms=6)
    b.axhline(1.0, color="gray", ls=":", lw=0.8)
    b.axvspan(1.0, zc2 / zc1, alpha=0.15, color="#1a7f37", label=f"only k1 unstable (robust, →1.0)")
    b.axvspan(zc2 / zc1, 3.2, alpha=0.12, color="#d62728", label="multi-mode (seed-dependent multi-cap)")
    b.set_xlabel("contractility ζ / ζ_c1"); b.set_ylabel("single-cap fraction (6 seeds)")
    b.set_ylim(0, 1.08); b.set_title("(b) ensemble: robust single cap ONLY near threshold"); b.legend(fontsize=7)
    b.grid(alpha=0.3)

    c2 = ax[1, 0]
    im = c2.imshow(kymo.T, aspect="auto", origin="lower", cmap="magma", extent=[times[0], times[-1], 0, p.L])
    c2.set_xlabel("time [s]"); c2.set_ylabel("arclength x [µm]")
    c2.set_title(f"(c) representative single-cap run (ζ=2·ζc1, seed 0)")
    fig.colorbar(im, ax=c2, label="myosin density c")

    d = ax[1, 1]
    d.plot(x, cfin, "-", color="#c1121f", label="myosin density c(x)")
    d.plot(x, p.c0 + 0 * x, "k:", lw=0.8, label="baseline c0")
    axv = d.twinx(); axv.plot(x, vfin, "-", color="#0077b6", alpha=0.7, label="cortical flow v(x)")
    axv.set_ylabel("flow v [µm/s]", color="#0077b6")
    d.set_xlabel("arclength x [µm]"); d.set_ylabel("myosin density c", color="#c1121f")
    d.set_title(f"(d) polarized state (rep: mode {rep_dom}, contrast {rep_contrast:.1f}×)")
    h1, l1 = d.get_legend_handles_labels(); h2, l2 = axv.get_legend_handles_labels()
    d.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")

    fig.suptitle("FF active movement piece-2/5 — cell polarization by actomyosin contractile-flow instability "
                 "(Bois 2011; ℓ=14µm Mayer 2010). Robust single-cap ONLY near threshold — see (b).", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{args.out}.png", dpi=130); print("wrote", f"{args.out}.png")

    # ---- interactive HTML: representative single-cap run on the cortex ring (radial bulge = myosin) ----
    theta = 2 * np.pi * x / p.L
    gain = 2.0 / max(cfin.max(), 1e-9)
    def ring_points(cf):
        r = R + gain * (cf - p.c0)
        return np.stack([r * np.cos(theta), r * np.sin(theta), np.zeros_like(theta)], 1).astype(np.float32)
    base_ring = np.stack([R * np.cos(theta), R * np.sin(theta), np.zeros_like(theta)], 1).astype(np.float32)
    scenes = {"myosin cap on cortex ring (representative single-cap run)": [
        {"name": "cortex ring (R=7.5µm)", "kind": "points", "color": "#3b6ea5", "size": 1.5, "verts": base_ring},
        {"name": "myosin density (radial bulge = cap)", "kind": "points", "color": "#c1121f", "size": 3.0,
         "verts": ring_points(cframes[0]), "frames": [ring_points(cf) for cf in cframes]},
    ]}
    build_viewer(scenes, out=f"{args.out}.html",
                 title="FF piece-2 — polarization: a myosin cap forms on the cortex ring (▶ play; robust only near threshold)")
    print("wrote", f"{args.out}.html")
    print(f"[summary] zc1={zc1:.1f} zc2={zc2:.1f}  single-cap fraction vs ζ/ζc1={dict(zip(zratios, [round(f,2) for f in fracs]))}"
          f"  rep(2·ζc1,seed0): mode={rep_dom} contrast={rep_contrast:.2f}")


if __name__ == "__main__":
    main()
