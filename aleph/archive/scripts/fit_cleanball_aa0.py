"""Fit the DCM clean-ball wetting sweep to the experimental law A/A0 = a + b/R + c/R^2.

Reads the per-N spread pickles (cleanball_n{N}.pkl), extracts A/A0(t) (peak and
final/equilibrium) and the spheroid radius R0 = N^(1/3)·R_cell (equivalent-sphere
radius of the total cell volume), least-squares fits a+b/R+c/R^2 to BOTH the peak and
the final A/A0, and plots A/A0 vs R with the fit + per-N A/A0(step) trajectories.
Output → figs/. The magnitude is expected << the experimental 7-10 (structural limit
of a fine-grained 12-200 cell model); the deliverable is the FORM/trend + the fit.
"""
from __future__ import annotations
import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("aleph/outputs/h_dcm_two_stage")
R_CELL_UM = 7.5


def _fit(R, A):
    """Least squares A = a + b/R + c/R^2."""
    M = np.column_stack([np.ones_like(R), 1.0 / R, 1.0 / R ** 2])
    coef, *_ = np.linalg.lstsq(M, A, rcond=None)
    pred = M @ coef
    ss_res = float(((A - pred) ** 2).sum())
    ss_tot = float(((A - A.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return coef, r2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=int, nargs="+", default=[12, 30, 60, 100, 200])
    ap.add_argument("--prefix", default="cleanball_n")
    ap.add_argument("--out", default="figs/cleanball_aa0_fit.png")
    args = ap.parse_args()

    rows = []
    trajectories = []
    for N in args.ns:
        p = OUT / f"{args.prefix}{N}.pkl"
        if not p.exists():
            print(f"  skip N={N}: {p} missing")
            continue
        with open(p, "rb") as fh:
            d = pickle.load(fh)
        s2 = d["spread"]
        A0 = s2.get("A0_topdown_um2") or s2["diags"][0]["topdown_um2"]
        aa = np.array([g["topdown_um2"] / A0 for g in s2["diags"]])
        steps = np.array(s2["steps"])
        R0 = (N ** (1.0 / 3.0)) * R_CELL_UM
        rows.append((N, R0, float(aa.max()), float(aa[-1])))
        trajectories.append((N, steps, aa))
        print(f"  N={N:>3} R0={R0:5.1f}um  A/A0 peak={aa.max():.3f} final={aa[-1]:.3f}")

    if len(rows) < 3:
        print(f"only {len(rows)} points — need >=3 for the 3-param fit; plotting raw")

    rows = np.array(rows, dtype=float)
    N_, R, Apk, Afi = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]

    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
    Rg = np.linspace(R.min() * 0.95, R.max() * 1.05, 200)
    for A, lab, col in [(Apk, "peak", "C3"), (Afi, "final/equil", "C0")]:
        ax[0].plot(R, A, "o", color=col, ms=8, label=f"{lab} A/A0")
        if len(rows) >= 3:
            coef, r2 = _fit(R, A)
            ax[0].plot(Rg, coef[0] + coef[1] / Rg + coef[2] / Rg ** 2, "-", color=col,
                       label=f"  fit a={coef[0]:.2f} b={coef[1]:.1f} c={coef[2]:.0f} (r²={r2:.3f})")
    ax[0].axhline(1.0, ls=":", color="0.6")
    ax[0].set_xlabel("spheroid radius R0 = N^(1/3)·R_cell  [µm]")
    ax[0].set_ylabel("A/A0 (top-down silhouette)")
    ax[0].set_title("A/A0 = a + b/R + c/R²  (clean ball + wetting, proxy OFF)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    for N, steps, aa in trajectories:
        ax[1].plot(steps, aa, "-o", ms=3, label=f"N={int(N)}")
    ax[1].axhline(1.0, ls=":", color="0.6")
    ax[1].set_xlabel("spread step"); ax[1].set_ylabel("A/A0")
    ax[1].set_title("A/A0(step) per N — overshoot then relax")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    fig.suptitle("DCM spheroid spreading sweep — clean FCC ball + substrate wetting "
                 "(MCF7 W_cs=2.85e-3, proxy OFF)\nmagnitude bounded << experimental 7-10 "
                 "(fine-grained-model structural limit); FORM is the deliverable",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    pth = OUT / args.out
    pth.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pth, dpi=140)
    print(f"wrote {pth}")


if __name__ == "__main__":
    main()
