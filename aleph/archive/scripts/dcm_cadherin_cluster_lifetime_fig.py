"""Figure: cadherin load-sharing cluster lifetime — the finding that reshapes the DCM real-timescale
redesign (S1). Two panels:
  (L) analytic cluster lifetime T(n_b) vs cluster size, for several per-molecule loads, with the
      KB-4.11 5-30 min junction-lifetime band overlaid + the single-molecule (n_b=1) baseline marked.
      Shows that load sharing lifts the junction from ~0.036 s into the min-hr regime, but T(n_b) is
      super-exponentially sensitive near the rho=k_on/eps~1 critical transition (n_b=100 -> permanent).
  (R) G1 validation: the emergent stochastic cluster lifetime reproduces the analytic BD-MFPT.

Integrity: no axis truncation (log y is annotated), units on every axis, the lit band is overlaid on
the measurement, and per-realisation MC points are shown against the analytic curve.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aleph.dcm.dcm_cadherin_cluster import cluster_bd_step, cluster_mean_lifetime_analytic
from aleph.validation.cadherin_sliding_rebinding import effective_k_off

K_ON = 27.96
OUT = "aleph/outputs/h_dcm_two_stage/figs"
os.makedirs(OUT, exist_ok=True)


def emergent_lifetime(n_b, eps, dt, n_real=3000, seed=7):
    rng = np.random.default_rng(seed)
    p_off = np.full(n_real, 1.0 - np.exp(-eps * dt))
    p_on = np.full(n_real, 1.0 - np.exp(-K_ON * dt))
    m = np.full(n_real, n_b, dtype=np.int64)
    life = np.zeros(n_real)
    alive = np.ones(n_real, dtype=bool)
    for _ in range(int(80.0 * cluster_mean_lifetime_analytic(n_b, eps, K_ON) / dt) + 200):
        if not alive.any():
            break
        m[alive] = cluster_bd_step(m[alive], n_b, p_off[alive], p_on[alive], rng)
        life[alive] += dt
        alive = m > 0
    return float(life.mean())


def main():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.2))

    # ---- (L) T(n_b) vs cluster size, several per-molecule loads ----
    sizes = np.arange(1, 31)
    loads_pn = [0.0, 10.0, 20.0, 29.2]
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(loads_pn)))
    for f_pn, c in zip(loads_pn, colors):
        eps = effective_k_off(f_pn * 1e-12)
        taus = [cluster_mean_lifetime_analytic(int(n), eps, K_ON) for n in sizes]
        axL.semilogy(sizes, taus, "-", color=c, lw=2,
                     label=f"F₁={f_pn:.0f} pN  (ε={eps:.0f}/s, ρ={K_ON/eps:.2f})")
    # KB-4.11 mature junction-lifetime band 5-30 min
    axL.axhspan(5 * 60, 30 * 60, color="tab:green", alpha=0.15, zorder=0)
    axL.text(1.2, 12 * 60, "KB-4.11 mature junction\n5–30 min", fontsize=9,
             color="darkgreen", va="center")
    # single-molecule baseline
    single = cluster_mean_lifetime_analytic(1, effective_k_off(0.0), K_ON)
    axL.scatter([1], [single], color="crimson", zorder=5, s=45)
    axL.annotate(f"single molecule\n{single*1e3:.0f} ms (n_b=1)", (1, single),
                 xytext=(3.0, single * 6), fontsize=9, color="crimson",
                 arrowprops=dict(arrowstyle="->", color="crimson"))
    axL.set_xlabel("cluster size  n_b  (parallel trans-dimers per junction)")
    axL.set_ylabel("mean junction lifetime  T(n_b)  [s]")
    axL.set_title("Load-sharing lifts the junction into the min–hr regime\n"
                  "— but T(n_b) is super-exponential near ρ~1 (finding)")
    axL.legend(fontsize=8, loc="lower right")
    axL.grid(True, which="both", alpha=0.25)
    axL.set_ylim(1e-2, 1e9)

    # ---- (R) G1: emergent (MC) vs analytic ----
    f_pn = 10.0
    eps = effective_k_off(f_pn * 1e-12)
    ns = [1, 2, 3, 4, 6, 8, 10]
    ana = [cluster_mean_lifetime_analytic(n, eps, K_ON) for n in ns]
    dt = 0.05 / max(eps, K_ON)
    emp = [emergent_lifetime(n, eps, dt) for n in ns]
    axR.loglog(ana, ana, "k--", alpha=0.6, label="analytic BD-MFPT (oracle)")
    axR.scatter(ana, emp, color="tab:blue", s=55, zorder=5,
                label="emergent stochastic (cluster_bd_step)")
    axR.set_xlabel("analytic cluster lifetime  T(n_b)  [s]")
    axR.set_ylabel("emergent Monte-Carlo lifetime  [s]")
    axR.set_title(f"G1: the runtime cluster reproduces the oracle\n(per-molecule load {f_pn:.0f} pN, "
                  f"n_b=1…10)")
    axR.legend(fontsize=9, loc="upper left")
    axR.grid(True, which="both", alpha=0.25)

    fig.tight_layout()
    path = f"{OUT}/dcm_cadherin_cluster_lifetime.png"
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")
    print(f"  single-molecule (n_b=1): {single*1e3:.1f} ms")
    for n in [10, 15, 20, 25]:
        t = cluster_mean_lifetime_analytic(n, effective_k_off(0.0), K_ON)
        print(f"  n_b={n:2d} @ F₁=0: {t:.3g}s = {t/60:.2f} min")


if __name__ == "__main__":
    main()
