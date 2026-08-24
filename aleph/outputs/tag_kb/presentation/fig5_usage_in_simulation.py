"""Figure 5 — How the KB is used continuously in the simulation workflow.

Clean vertical provenance loop (left/center) around a single simulation constant
+ a clearly separated live grounding-pass panel (right). No overlapping boxes.
Renders to figs/usage_in_simulation.png.
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ---- shared style ---------------------------------------------------------
C_TAG = "#3C8D7A"        # teal
C_SYS = "#5B4B8A"        # indigo (THIS system)
C_EDGE = "#444"          # neutral box edge / arrows
C_FAIL = "#C0392B"       # forbidden / failure accent (red)
C_PASS = "#2E8B57"       # verified / pass (green)
C_FILL = "#F5F6FA"       # light fill

OUT = "/Users/sw1/ffn_cellsim-platform/aleph/outputs/tag_kb/presentation/figs/usage_in_simulation.png"


def box(ax, x, y, w, h, text, edge=C_SYS, fc=C_FILL, fs=11, weight="normal", tc="#222"):
    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.3,rounding_size=0.08",
            facecolor=fc, edgecolor=edge, linewidth=2, zorder=2,
        )
    )
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, weight=weight,
            color=tc, zorder=3)


def arrow(ax, p0, p1, color=C_EDGE, rad=0.0, label=None, lab_pos=None,
          lab_color="#555", lw=1.9):
    ax.add_patch(
        FancyArrowPatch(
            p0, p1, arrowstyle="-|>", mutation_scale=18, color=color,
            lw=lw, connectionstyle=f"arc3,rad={rad}", zorder=1,
            shrinkA=2, shrinkB=2,
        )
    )
    if label:
        if lab_pos is None:
            lab_pos = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
        ax.text(*lab_pos, label, ha="center", va="center", fontsize=9,
                style="italic", color=lab_color, zorder=4)


fig, ax = plt.subplots(figsize=(13.5, 8.2))
ax.set_xlim(0, 13.5)
ax.set_ylim(0, 8.2)
ax.set_axis_off()

# Title + subtitle
fig.text(0.5, 0.965, "How the KB is used continuously in the simulation workflow",
         ha="center", fontsize=15, weight="bold", color="#222")
fig.text(0.5, 0.923,
         "every simulation constant is re-bound to a verified citation on every push — drift blocks the merge",
         ha="center", fontsize=10, color="#777")

# ---- vertical provenance loop (single column, top -> bottom) --------------
CX = 4.3            # column centre x
W = 4.2             # box width
ys = [7.05, 5.25, 3.45, 1.55]   # config, param_audit, source_audit, CI gate

box(ax, CX, ys[0], W, 1.05,
    "config YAML   # KU-3.5\n$\\gamma_{cortex}=0.5\\times10^{-3}$ N/m   (the constant)",
    edge=C_SYS, fc="#ECE8F5", fs=11, weight="bold", tc=C_SYS)
box(ax, CX, ys[1], W, 1.05,
    "param_audit\nconstant  →  KU  →  citation_key", edge=C_TAG, fs=11)
box(ax, CX, ys[2], W, 1.05,
    "source_audit verdict\n(CrossRef:  OK / CHECK / DEAD)", edge=C_TAG, fs=11)
box(ax, CX, ys[3], W, 1.05,
    "CI params-gate\nexit 1 on DRIFT  →  blocks the merge",
    edge=C_PASS, fs=11, weight="bold", tc=C_PASS)

half = 0.525
# straight down-arrows between consecutive boxes, with the step + frequency
steps = [
    ("1.  bind:  KU tag → citation", "per refresh / push", C_EDGE),
    ("2.  is that citation verified?", "per result claim", C_TAG),
    ("3.  enforce on every push", "per push (CI)", C_PASS),
]
for i, (lab, freq, col) in enumerate(steps):
    y0, y1 = ys[i] - half, ys[i + 1] + half
    arrow(ax, (CX, y0), (CX, y1), color=col)
    ax.text(CX + 0.25, (y0 + y1) / 2 + 0.10, lab, ha="left", va="center",
            fontsize=9.5, color="#333", zorder=4)
    ax.text(CX + 0.25, (y0 + y1) / 2 - 0.20, freq, ha="left", va="center",
            fontsize=8.5, style="italic", color=col, zorder=4)

# return arrow (gate passes -> the constant stays), bowing out to the LEFT
arrow(ax, (CX - W / 2, ys[3]), (CX - W / 2, ys[0]), color=C_PASS, rad=-0.32)
ax.text(0.85, (ys[0] + ys[3]) / 2, "gate passes\n→ constant stays\n\n(DRIFT → blocked)",
        ha="center", va="center", fontsize=9, style="italic", color=C_PASS, zorder=4)

# ---- right panel: live grounding-pass catches (clearly separated) ---------
px0, px1 = 7.7, 12.95
py0, py1 = 0.85, 6.95
ax.add_patch(
    FancyBboxPatch(
        (px0, py0), px1 - px0, py1 - py0,
        boxstyle="round,pad=0.3,rounding_size=0.05",
        facecolor="#FBEEEC", edgecolor=C_FAIL, linewidth=2, zorder=2,
    )
)
pcx = (px0 + px1) / 2
ax.text(pcx, py1 - 0.55, "live grounding pass — 2026-06-19/20",
        ha="center", va="center", fontsize=11.5, weight="bold", color=C_FAIL, zorder=3)
ax.text(pcx, py1 - 1.05, "the gates' ongoing catches (not the old 3-citation story)",
        ha="center", va="center", fontsize=9, style="italic", color="#9A5B52", zorder=3)
catches = [
    ("✗  forbidden basal-footprint A/A0",
     "headlined in 5 DCM REPORTs — the sanctioned\ntop-down metric shows COMPACTION (0.597)"),
    ("✗  PNG-only \"breakthrough\" numbers",
     "1.87 / r²=0.96 lived only in figure titles,\nno committed data  →  NEEDS_REGEN"),
    ("✗  prose-inflated \"44×\" speedup",
     "corrected to the project's own record  →  39.47×"),
]
yy = py1 - 1.95
for head, body in catches:
    ax.text(px0 + 0.35, yy, head, ha="left", va="center",
            fontsize=10, weight="bold", color=C_FAIL, zorder=3)
    ax.text(px0 + 0.55, yy - 0.50, body, ha="left", va="top",
            fontsize=8.8, color="#444", zorder=3)
    yy -= 1.62

# one arrow from the CI gate to the panel
arrow(ax, (CX + W / 2, ys[3]), (px0 - 0.05, py0 + 0.9), color=C_FAIL,
      rad=-0.15, label="gates fire", lab_pos=(6.55, 1.05), lab_color=C_FAIL)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print(OUT)
