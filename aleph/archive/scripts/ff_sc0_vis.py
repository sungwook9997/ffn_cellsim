"""SC0 solver-profile figure — CG-domination (robust across resolution + full diag) + the CG-iteration ladder.

Renders the reinforced SC0 finding (native NF=38000 & 70686, A5000):
  - the implicit step is CG-DOMINATED (assembly 8.6-22.7%, CG 69-89%) at BOTH native definitions and WITH the
    full clutch/substrate diagonal → §4.2 matrix-free is the wrong lever;
  - the CG iteration count is set by conditioning targets a preconditioner/deflation can remove: the a_com modal
    rigid-COM mode (71→135) AND the stiff heterogeneous basal clutch/substrate diagonal (+diag RAISES iters,
    135→210), and it grows with N — so §4.4 preconditioner is the win, worsening (hence more valuable) at native.

    python ff_sc0_vis.py   # reads outputs/ff_single_opt/sc0_profiles/*.json
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROF = os.path.join(HERE, "outputs", "ff_single_opt", "sc0_profiles")


def _load(name):
    with open(os.path.join(PROF, name)) as fh:
        return json.load(fh)


def main() -> None:
    nocom = _load("sc0_native_nocom.json")          # 70k, no a_com
    com70 = _load("sc0_native_comdrag.json")        # 70k, a_com
    full70 = _load("sc0_native70k_full.json")       # 70k, a_com + clutch/substrate diag
    com38 = _load("sc0_native38k_comdrag.json")     # 38k, a_com
    full38 = _load("sc0_native38k_full.json")       # 38k, a_com + diag

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4))

    # ---- Panel 1: CG-domination is robust across resolution + diag ----
    configs = [("70k\nno-a_com", nocom), ("70k\n+a_com", com70), ("70k\n+a_com\n+diag", full70),
               ("38k\n+a_com", com38), ("38k\n+a_com\n+diag", full38)]
    labels = [c[0] for c in configs]
    asm = [c[1]["assemble_frac"] * 100 for c in configs]
    cgp = [c[1]["cg_frac"] * 100 for c in configs]
    frc = [c[1]["force_frac"] * 100 for c in configs]
    x = range(len(configs))
    ax1.bar(x, asm, label="assemble (COO→CSR)", color="#4C78A8")
    ax1.bar(x, cgp, bottom=asm, label="CG solve", color="#E45756")
    ax1.bar(x, frc, bottom=[a + c for a, c in zip(asm, cgp)], label="force eval", color="#59A14F")
    for i, c in enumerate(configs):
        ax1.text(i, 101, f"{cgp[i]:.0f}%", ha="center", va="bottom", fontsize=10, weight="bold", color="#B0281A")
        ax1.text(i, 4, f"{c[1]['mean_t_full_ms']:.0f}ms", ha="center", va="bottom", fontsize=8, color="white")
    ax1.set_xticks(list(x)); ax1.set_xticklabels(labels, fontsize=8.5)
    ax1.set_ylabel("% of implicit step wall-time"); ax1.set_ylim(0, 112)
    ax1.axhline(50, color="gray", ls=":", lw=0.8)
    ax1.set_title("CG-DOMINATED at every native config (assembly 9-23%, CG 70-89%)\n→ §4.2 matrix-free is the wrong lever", fontsize=10.5)
    ax1.legend(loc="lower center", fontsize=8.5, ncol=3)

    # ---- Panel 2: the CG-iteration ladder (conditioning targets a preconditioner removes) ----
    ladder = [("70k\nno-a_com", nocom["mean_cg_iters"], "#54A24B"),
              ("70k\n+a_com", com70["mean_cg_iters"], "#F58518"),
              ("70k\n+a_com+diag", full70["mean_cg_iters"], "#E45756"),
              ("38k\n+a_com", com38["mean_cg_iters"], "#F58518"),
              ("38k\n+a_com+diag", full38["mean_cg_iters"], "#E45756")]
    lx = range(len(ladder))
    ax2.bar(lx, [v for _, v, _ in ladder], color=[c for _, _, c in ladder])
    for i, (_, v, _) in enumerate(ladder):
        ax2.text(i, v + 3, f"{v:.0f}", ha="center", va="bottom", fontsize=11, weight="bold")
    ax2.set_xticks(list(lx)); ax2.set_xticklabels([l for l, _, _ in ladder], fontsize=8.5)
    ax2.set_ylabel("unpreconditioned CG iterations / step"); ax2.set_ylim(0, 245)
    ax2.annotate("a_com rigid mode\n(§4.4 deflation, closed-form)", xy=(1, com70["mean_cg_iters"]), xytext=(0.15, 175),
                 arrowprops=dict(arrowstyle="->", color="#F58518", lw=1.5), fontsize=8.5, color="#B5650E")
    ax2.annotate("stiff basal diag\n(diagonally-dominant → Jacobi-amenable)", xy=(2, full70["mean_cg_iters"]),
                 xytext=(2.05, 228), arrowprops=dict(arrowstyle="->", color="#E45756", lw=1.5), fontsize=8.5, color="#B0281A")
    ax2.set_title("CG iters set by TWO conditioning targets a preconditioner removes\n(a_com 71→135; +stiff basal diag →210; grows with N)", fontsize=10.5)

    fig.suptitle("SC0 native profile (A5000) — CG-domination is robust; the lever is CG iteration count (deflation+block preconditioner), not assembly (matrix-free)",
                 fontsize=10.5, y=1.01)
    fig.tight_layout()
    out = os.path.join(PROF, "figs", "sc0_solver_decomposition.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
