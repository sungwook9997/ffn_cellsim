"""Plot KU-3.5 grip-walk per-sample trajectories from a run log.

Parses the `[grip_walk] sample=.. g_tot=.. s_grip=..nm engaged=.. adv=..` lines
and plots s_grip, adv (bead re-targets), gamma_total, and engaged-head count vs
sample index. The at-a-glance "did transport engage / did tension build" figure.

Usage:
    python ffn_sim/scripts/viz_ku35_trajectory.py <run.log> <out.png>
"""
from __future__ import annotations

import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAT = re.compile(
    r"sample=(\d+)/\d+.*?g_tot=([\d.eE+-]+).*?s_grip=([\d.eE+-]+)nm"
    r"\s+engaged=(\d+)\s+adv=(\d+)"
)


def main():
    log, out = sys.argv[1], sys.argv[2]
    rows = []
    for line in open(log):
        m = PAT.search(line)
        if m:
            rows.append([int(m.group(1)), float(m.group(2)), float(m.group(3)),
                         int(m.group(4)), int(m.group(5))])
    if not rows:
        print("no sample lines parsed"); return
    a = np.array(rows, dtype=float)
    s, g_tot, s_grip, engaged, adv = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))

    ax[0, 0].plot(s, s_grip, "o-", color="purple")
    ax[0, 0].set_title("commanded grip stretch s_grip"); ax[0, 0].set_ylabel("s_grip (nm)")
    ax[0, 0].axhline(500, ls="--", c="gray", lw=0.8, label="ℓ₀=500 nm")
    ax[0, 0].legend(fontsize=8)

    ax[0, 1].plot(s, adv, "o-", color="crimson")
    ax[0, 1].set_title("bead re-targets  adv  (TRANSPORT signal)")
    ax[0, 1].set_ylabel("step_advances (cumulative)")
    onset = s[adv > 0]
    if onset.size:
        ax[0, 1].axvline(onset[0], ls=":", c="crimson", lw=1,
                         label=f"adv>0 onset @ sample {int(onset[0])}")
        ax[0, 1].legend(fontsize=8)

    ax[1, 0].plot(s, g_tot * 1e3, "o-", color="teal")
    ax[1, 0].axhspan(0.35, 0.65, color="green", alpha=0.15, label="KU-3.5 band")
    ax[1, 0].set_yscale("log")
    ax[1, 0].set_title("γ_total (DISTORTED at high v0_accel — no percolation)")
    ax[1, 0].set_ylabel("γ_total (×10⁻³ mN/m)"); ax[1, 0].set_xlabel("sample")
    ax[1, 0].legend(fontsize=8)

    ax[1, 1].plot(s, engaged, "o-", color="darkorange")
    ax[1, 1].set_title("engaged myosin heads"); ax[1, 1].set_ylabel("n_engaged")
    ax[1, 1].set_xlabel("sample")

    fig.suptitle(f"{log.split('/')[-1]} — native n_fil=38000, v0×1500 (transport diagnostic)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=130)
    print(f"wrote {out}  ({len(rows)} samples; adv {int(adv.min())}→{int(adv.max())}, "
          f"s_grip max {s_grip.max():.0f}nm, γ_total {g_tot.min()*1e3:.2e}-{g_tot.max()*1e3:.2e} mN/m)")


if __name__ == "__main__":
    main()
