"""Visualize the cortical-tension γ A/B: CONNECTED mesh vs FRAGMENTED baseline.

The CORTEX REBUILD payoff test (Kadzik-Munro 2026 ⭐: connectivity is the
prerequisite for force transmission). Parses the two production logs from
``mcf7_fullcell_stage1`` (``--connected-mesh`` vs the fragmented baseline, same
params) and renders γ_soft / γ_rigid / γ_total per sample vs the KU-3.5 band
[0.35, 0.65] mN/m — does building a CONNECTED percolated cortex lift the
γ-floor?

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.viz_connected_gamma_ab \
          --connected outputs/h3/production/connected_mesh_ab/CONNECTED_gpu.log \
          --fragmented outputs/h3/production/connected_mesh_ab/FRAGMENTED_cpu.log
Output: ffn_sim/outputs/h3/figs/connected_gamma_ab.png
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parents[1]
FIG = PKG / "outputs" / "h3" / "figs"

_LINE = re.compile(
    r"r/r0=([0-9.eE+-]+)\s+g_soft=([0-9.eE+-]+)\s+g_rigid=([0-9.eE+-]+)\s+"
    r"g_tot=([0-9.eE+-]+)"
)


def parse_log(path: Path):
    """Return (r_over_r0, g_soft, g_rigid, g_tot) arrays in mN/m per sample."""
    if not path.exists():
        return None
    rows = []
    for m in _LINE.finditer(path.read_text()):
        rows.append([float(m.group(i)) for i in range(1, 5)])
    if not rows:
        return None
    a = np.array(rows)
    return dict(r=a[:, 0], g_soft=a[:, 1], g_rigid=a[:, 2], g_tot=a[:, 3])


def make_figure(conn, frag, out_png, band=(0.35, 0.65)):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    comps = [("g_soft", "γ_soft (compliant)"),
             ("g_rigid", "γ_rigid (constrained)"),
             ("g_tot", "γ_total")]
    for ax, (key, label) in zip(axes, comps):
        ax.axhspan(band[0], band[1], color="tab:green", alpha=0.15,
                   label=f"KU-3.5 band [{band[0]},{band[1]}]")
        for name, d, c in (("FRAGMENTED (z≈1.3)", frag, "tab:red"),
                           ("CONNECTED (z≈3.3)", conn, "tab:blue")):
            if d is None:
                continue
            y = d[key]
            x = np.arange(1, len(y) + 1)
            ax.plot(x, y, "o-", color=c, label=f"{name}: ⟨{y.mean():.3e}⟩")
        ax.set_title(label)
        ax.set_xlabel("sample"); ax.set_ylabel("γ (mN/m)")
        ax.set_yscale("log")
        ax.legend(fontsize=8)

    cg = conn["g_tot"].mean() if conn else float("nan")
    fg = frag["g_tot"].mean() if frag else float("nan")
    lift = cg / fg if (frag and fg > 0) else float("nan")
    verdict = (f"CONNECTED ⟨γ_tot⟩={cg:.3e} vs FRAGMENTED ⟨γ_tot⟩={fg:.3e} mN/m  "
               f"→ connectivity ×{lift:.1f} on γ_total" if frag and conn else
               f"CONNECTED ⟨γ_tot⟩={cg:.3e} mN/m")
    fig.suptitle(
        f"Cortex REBUILD payoff — does connectivity lift the γ-floor?  "
        f"(Kadzik-Munro 2026: connectivity ⇒ force transmission)\n{verdict}",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return verdict


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--connected", type=Path,
                    default=PKG / "outputs/h3/production/connected_mesh_ab/CONNECTED_gpu.log")
    ap.add_argument("--fragmented", type=Path,
                    default=PKG / "outputs/h3/production/connected_mesh_ab/FRAGMENTED_cpu.log")
    args = ap.parse_args()
    conn = parse_log(args.connected)
    frag = parse_log(args.fragmented)
    if conn is None and frag is None:
        print("[viz] no γ samples parsed yet — runs may still be warming up.")
        return 1
    out = FIG / "connected_gamma_ab.png"
    verdict = make_figure(conn, frag, out)
    print(f"[viz] {verdict}")
    print(f"[viz] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
