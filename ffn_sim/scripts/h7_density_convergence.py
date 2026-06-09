"""H.7 per-SF traction sign-convergence vs engaged-head count (PI-directed high-density test).

The loop23 "decisive +131 pN" was refuted as a single-seed fluke: at ~50 engaged heads the per-SF
coherent traction differential is sign-unstable across realizations (mean ≈ 0, std ≈ 350 pN). The
DENSITY hypothesis (H7_SF_ARRAY_TRACTION_SCALEUP §3): as the number of PARALLEL engaged heads grows,
a geometric contractile bias from the sarcomeric organization should dominate the small-N placement
noise → the per-SF differential converges to a stable, sign-definite contractile value.

This script ingests the parallel-density sweep (wider bundle → more parallel minifilaments) +
the ~50-head baseline, groups by density level, and tests convergence: per level it reports the
across-seed mean ± std and the coefficient of variation |std/mean|, and plots mean±std vs mean
engaged heads. VERDICT: convergence = CV shrinks AND the sign becomes definite as N_heads grows
(→ rectification real, magnitude is the density lever); persistence = sign stays unstable (→ the
sarcomeric construction does not rectify; halt the traction line).

Usage: python ffn_sim/scripts/h7_density_convergence.py
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_PROD = Path("ffn_sim/outputs/h7/production")
_FIG = Path("ffn_sim/outputs/h7/figs/h7_density_convergence.png")


def _load(pattern):
    """Return list of (engaged, differential_pN) for a glob."""
    out = []
    for f in sorted(glob.glob(str(_PROD / pattern))):
        try:
            d = json.loads(Path(f).read_text())
            out.append((d.get("engaged_heads_mean", float("nan")),
                        d["coherent_differential_pN"]))
        except Exception:
            pass
    return out


def main() -> int:
    # level label → (glob); baseline = the refutation ensemble (~50 heads, br=400, 16 motors)
    levels = [
        ("P1 ~50 heads\n(br400, 16mf)", "h7_sarc_cpu_s*.json"),
        ("P2 (br800, 28mf)", "h7_dens_P2_*.json"),
        ("P3 (br1200, 55mf)", "h7_dens_P3_*.json"),
        ("P4 (br1600, 100mf)", "h7_dens_P4_*.json"),
    ]
    rows = []
    for label, pat in levels:
        data = _load(pat)
        if not data:
            print(f"  [skip] {label}: no jsons ({pat})")
            continue
        eng = np.array([d[0] for d in data])
        vals = np.array([d[1] for d in data])
        mean = vals.mean()
        std = vals.std(ddof=1) if len(vals) > 1 else float("nan")
        cv = abs(std / mean) if mean != 0 else float("inf")
        npos = int((vals > 0).sum()); nneg = int((vals < 0).sum())
        rows.append((label, np.nanmean(eng), mean, std, cv, len(vals), npos, nneg))
        print(f"  {label.splitlines()[0]:22s} | N={len(vals)} | engaged≈{np.nanmean(eng):5.0f} | "
              f"mean={mean:+8.1f} ± {std:6.1f} pN | CV={cv:5.2f} | {npos}+/{nneg}-")

    if not rows:
        print("  no data yet"); return 0

    # verdict
    cvs = [r[4] for r in rows if np.isfinite(r[4])]
    converging = len(cvs) >= 2 and cvs[-1] < cvs[0] and rows[-1][2] > 0 and rows[-1][6] >= rows[-1][7]
    verdict = ("CONVERGING → rectification real, magnitude = density lever" if converging
               else "SIGN STILL UNSTABLE → sarcomeric construction does not robustly rectify")
    print(f"\n  VERDICT: {verdict}")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    engs = [r[1] for r in rows]
    means = [r[2] for r in rows]
    stds = [r[3] for r in rows]
    ax[0].errorbar(engs, means, yerr=stds, fmt="o-", capsize=5, color="#1f77b4", lw=1.5,
                   label="across-seed mean ± std")
    ax[0].axhline(0, color="k", lw=1, ls=":")
    ax[0].set_xlabel("mean engaged heads per SF")
    ax[0].set_ylabel("coherent traction differential (ON−OFF) [pN]")
    ax[0].set_title("(a) per-SF traction vs parallel engaged-head count")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    cv_vals = [r[4] for r in rows]
    ax[1].plot(engs, cv_vals, "s-", color="#d62728", lw=1.5)
    ax[1].axhline(1.0, color="k", lw=0.8, ls=":", label="CV=1 (sign indeterminate)")
    ax[1].set_xlabel("mean engaged heads per SF")
    ax[1].set_ylabel("coefficient of variation |std/mean|")
    ax[1].set_title("(b) sign-stability vs density\n(↓ = converging to definite sign)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    fig.suptitle(f"H.7 per-SF traction sign-convergence — {verdict}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"  fig → {_FIG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
