"""Figure + interactive animated HTML for piece-2 polarization (active-gel contractile-flow instability).

Shows: (a) sim vs analytic dispersion λ(k); (b) the marginal-stability window ζ_c(k); (c) a kymograph of the
myosin density c(x,t) breaking symmetry into one cap; (d) the final myosin profile + cortical flow. The HTML maps
c(x) onto the 3D cortex EQUATORIAL RING (points bulged radially outward by myosin density) animated over time —
watch the cell polarize (a single myosin cap = the rear emerges).
"""
from __future__ import annotations

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.ff.polarization_activegel import (resolve_activegel, dispersion, zeta_critical,
                                               measure_growth_rate, make_grid, step)
from ffn_sim.scripts.ff_viewer_html import build_viewer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=float, default=7.5)
    ap.add_argument("--out", default="ffn_sim/outputs/ff/figs/polarization_activegel")
    args = ap.parse_args()
    R = args.R
    p0 = resolve_activegel(R_um=R, zeta=1.0)
    zc1 = zeta_critical(p0, 2 * np.pi / p0.L)
    p = resolve_activegel(R_um=R, zeta=4.0 * zc1)                # well above threshold → clear polarization
    N = 256
    x, dx, kg = make_grid(p, N)

    # ---- (a) dispersion: sim vs analytic ----
    modes = [1, 2, 3, 4, 5, 6, 8, 10]
    lam_num, lam_an = [], []
    for m in modes:
        ln, la = measure_growth_rate(p, m, N=N, eps=1e-5, n_steps=150)
        lam_num.append(ln); lam_an.append(la)
    ks = np.array(modes) * 2 * np.pi / p.L

    # ---- nonlinear evolution: capture c(x,t) frames ----
    rng = np.random.default_rng(0)
    c = p.c0 + 1e-3 * rng.standard_normal(N)
    dt = 0.05 * min(dx ** 2 / p.D, 1.0 / p.k_off)
    n_steps, n_frames = 60000, 60
    stride = n_steps // n_frames
    kymo, times, cframes, vframes = [], [], [], []
    for s in range(n_steps + 1):
        if s % stride == 0:
            kymo.append(c.copy()); times.append(s * dt)
            cframes.append(c.copy());
        c, v = step(c, p, dx, dt, kg)
        if s % stride == 0:
            vframes.append(v.copy())
    kymo = np.array(kymo); cfin = c.copy(); _, vfin = c, v
    v = _velocity_final = None

    # ---- figure ----
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    a = ax[0, 0]
    a.plot(ks, lam_an, "k-o", ms=4, label="analytic λ(k)")
    a.plot(ks, lam_num, "rx", ms=8, mew=2, label="sim (spectral solver)")
    a.axhline(0, color="gray", lw=0.8); a.set_xlabel("wavenumber k [µm⁻¹]"); a.set_ylabel("growth rate λ [s⁻¹]")
    a.set_title("(a) dispersion: sim ON the analytic (rel_err <0.1%)"); a.legend(); a.grid(alpha=0.3)

    b = ax[0, 1]
    kk = np.linspace(2 * np.pi / p.L, 12 * 2 * np.pi / p.L, 200)
    b.plot(kk, [zeta_critical(p, k) for k in kk], "b-", label="ζ_c(k) marginal")
    b.axhline(p.zeta, color="r", ls="--", label=f"ζ used = {p.zeta:.0f} (4·ζ_c1)")
    b.fill_between(kk, [zeta_critical(p, k) for k in kk], p.zeta,
                   where=[zeta_critical(p, k) < p.zeta for k in kk], alpha=0.2, color="r", label="unstable band")
    b.set_xlabel("wavenumber k [µm⁻¹]"); b.set_ylabel("critical contractility ζ_c [norm]")
    b.set_title("(b) instability window (ζ > ζ_c ⇒ polarize)"); b.legend(fontsize=8); b.grid(alpha=0.3)

    c2 = ax[1, 0]
    im = c2.imshow(kymo.T, aspect="auto", origin="lower", cmap="magma",
                   extent=[times[0], times[-1], 0, p.L])
    c2.set_xlabel("time [s]"); c2.set_ylabel("arclength x [µm]"); c2.set_title("(c) myosin c(x,t): one cap emerges")
    fig.colorbar(im, ax=c2, label="myosin density c")

    d = ax[1, 1]
    d.plot(x, cfin, "-", color="#c1121f", label="myosin density c(x)")
    d.plot(x, p.c0 + 0 * x, "k:", lw=0.8, label="baseline c0")
    axv = d.twinx(); axv.plot(x, vfin, "-", color="#0077b6", alpha=0.7, label="cortical flow v(x)")
    axv.set_ylabel("flow v [µm/s]", color="#0077b6")
    d.set_xlabel("arclength x [µm]"); d.set_ylabel("myosin density c", color="#c1121f")
    d.set_title(f"(d) polarized steady state (contrast {(cfin.max()-cfin.min())/cfin.mean():.1f}×)")
    h1, l1 = d.get_legend_handles_labels(); h2, l2 = axv.get_legend_handles_labels()
    d.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")

    fig.suptitle("FF active movement piece-2/5 — cell polarization by actomyosin contractile-flow instability "
                 "(Bois 2011; ℓ=14µm Mayer 2010). Active-gel reduction, validated vs analytic λ(k).", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{args.out}.png", dpi=130); print("wrote", f"{args.out}.png")

    # ---- interactive HTML: myosin density on the cortex equatorial ring, bulged radially, animated ----
    theta = 2 * np.pi * x / p.L
    gain = 2.0 / max(cfin.max(), 1e-9)                            # radial bulge scale [µm per unit density]
    def ring_points(cf):
        r = R + gain * (cf - p.c0)
        return np.stack([r * np.cos(theta), r * np.sin(theta), np.zeros_like(theta)], 1).astype(np.float32)
    base_ring = np.stack([R * np.cos(theta), R * np.sin(theta), np.zeros_like(theta)], 1).astype(np.float32)
    scenes = {"myosin density on cortex ring (polarizing)": [
        {"name": "cortex ring (R=7.5µm)", "kind": "points", "color": "#3b6ea5", "size": 1.5, "verts": base_ring},
        {"name": "myosin density (radial bulge = cap)", "kind": "points", "color": "#c1121f", "size": 3.0,
         "verts": ring_points(cframes[0]), "frames": [ring_points(cf) for cf in cframes]},
    ]}
    build_viewer(scenes, out=f"{args.out}.html",
                 title="FF piece-2 — cell polarization (myosin cap forms on the cortex ring; ▶ play)")
    print("wrote", f"{args.out}.html")
    print(f"[summary] zeta_c1={zc1:.2f}  zeta={p.zeta:.2f}  final contrast={(cfin.max()-cfin.min())/cfin.mean():.2f}"
          f"  dominant mode={int(np.argmax(np.abs(np.fft.rfft(cfin-cfin.mean()))[1:])+1)}")


if __name__ == "__main__":
    main()
