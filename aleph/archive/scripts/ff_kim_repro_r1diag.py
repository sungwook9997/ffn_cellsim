"""Diagnose the r^-1 stress-transmission exponent for the FF ⟷ Slater/Kim 2021 reproduction.

The per-seed whole-band cylindrical exponent is noisy (−0.8..1.2, even unphysical negatives). Hypothesis (from
the DCM 1/r case, project-dcm-ecm-kim-1r-gap): the finite-cell NEAR field and the pinned-boundary FAR field
contaminate the whole-band fit, hiding a clean r^-1 in the intermediate field. This ensemble-averages the
cylindrical σ(r) profiles across the WLC ⟨z⟩≈3.2 seeds (same r grid) and fits NEAR / MID / FAR sub-bands
separately — the honest test of whether the r^-1 is there in the valid regime.

Run:  python -m aleph.scripts.ff_kim_repro_r1diag
Out:  aleph/outputs/ff/kim_repro/figs/kim_r1_diag.png  (+ printed sub-band exponents)
"""

from __future__ import annotations

import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAT = "aleph/outputs/ff/kim_repro/native"
FIGS = "aleph/outputs/ff/kim_repro/figs"


def _fit(r, s):
    ok = (r > 0) & (s > 0) & np.isfinite(s)
    if ok.sum() < 2:
        return float("nan")
    return float(-np.polyfit(np.log(r[ok]), np.log(s[ok]), 1)[0])


def main():
    os.makedirs(FIGS, exist_ok=True)
    paths = sorted(glob.glob(f"{NAT}/kim_repro_ens*.json")) + [f"{NAT}/kim_native_wlc.json", f"{NAT}/kim_repro_geom.json"]
    profs = []
    r_ref = None
    for p in paths:
        if not os.path.exists(p):
            continue
        pf = json.load(open(p))["on"]["profiles"][-1]
        r = np.asarray(pf["r_cyl"]); s = np.asarray(pf["sigma_cyl"])
        if r_ref is None:
            r_ref = r
        if r.shape == r_ref.shape and np.allclose(r, r_ref):        # same grid → stack directly
            profs.append(s)
    profs = np.array(profs)                                          # (N_seed, N_shell)
    n_seed = profs.shape[0]

    # ensemble mean profile (mean of σ per shell), normalized at the first MID shell for display
    mean_s = profs.mean(0)
    r = r_ref
    # sub-bands by radius: NEAR = inner third (finite cell), MID = middle third (clean), FAR = outer third (pinned bed)
    q1, q2 = np.percentile(r, [33, 66])
    near = r <= q1; mid = (r > q1) & (r <= q2); far = r > q2
    n_near, n_mid, n_far = _fit(r[near], mean_s[near]), _fit(r[mid], mean_s[mid]), _fit(r[far], mean_s[far])
    n_whole = _fit(r, mean_s)
    # per-seed MID exponents (robustness of the MID field across seeds)
    mid_per_seed = np.array([_fit(r[mid], profs[i][mid]) for i in range(n_seed)])

    print(f"[r1diag] N_seed={n_seed}  ensemble-mean σ(r) sub-band exponents:", flush=True)
    print(f"         whole = {n_whole:+.2f}   NEAR(r≤{q1:.0f}) = {n_near:+.2f}   "
          f"MID({q1:.0f}<r≤{q2:.0f}) = {n_mid:+.2f}   FAR(r>{q2:.0f}) = {n_far:+.2f}   (Kim r^-1 = 1.00)", flush=True)
    print(f"         per-seed MID exponent = {np.nanmean(mid_per_seed):+.2f} ± {np.nanstd(mid_per_seed):.2f} "
          f"(N={np.isfinite(mid_per_seed).sum()})", flush=True)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    for i in range(n_seed):
        a1.loglog(r, profs[i], color="crimson", alpha=0.2, lw=1)
    a1.loglog(r, mean_s, "o-", color="navy", lw=2.2, label="ensemble mean σ(r)")
    a1.loglog(r[mid], mean_s[mid][0] * (r[mid] / r[mid][0]) ** -1.0, "k--", lw=1.6, label="Kim r$^{-1}$ (over MID)")
    for m, c, lbl in [(near, "#888", "NEAR"), (mid, "#2ca02c", "MID"), (far, "#888", "FAR")]:
        a1.axvspan(r[m].min(), r[m].max(), color=c, alpha=0.06)
    a1.set_xlabel("r [µm]"); a1.set_ylabel("σ(r) cylindrical flux [Pa]")
    a1.set_title(f"Ensemble σ(r): MID-field exponent = {n_mid:.2f} vs Kim 1\n(whole {n_whole:.2f}, NEAR {n_near:.2f}, FAR {n_far:.2f})")
    a1.legend(fontsize=8); a1.grid(alpha=0.3, which="both")

    a2.hist(mid_per_seed[np.isfinite(mid_per_seed)], bins=np.linspace(-1, 3, 13), color="#2ca02c", alpha=0.7, edgecolor="k")
    a2.axvline(1.0, color="k", ls="--", lw=1.6, label="Kim r$^{-1}$")
    a2.axvline(np.nanmean(mid_per_seed), color="navy", lw=1.6, label=f"mean {np.nanmean(mid_per_seed):.2f}±{np.nanstd(mid_per_seed):.2f}")
    a2.set_xlabel("per-seed MID-field exponent"); a2.set_ylabel("count")
    a2.set_title(f"MID-field exponent per seed (N={n_seed})"); a2.legend(fontsize=8); a2.grid(alpha=0.3)

    fig.suptitle("FF ⟷ Slater/Kim 2021 — r$^{-1}$ diagnosis: is a clean r$^{-1}$ hiding in the intermediate field?", fontsize=12)
    fig.tight_layout()
    out = f"{FIGS}/kim_r1_diag.png"
    fig.savefig(out, dpi=130)
    print(f"[out] {out}", flush=True)


if __name__ == "__main__":
    main()
