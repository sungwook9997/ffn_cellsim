"""fig4_three_tier_comparison.py

Renders a clean comparison matrix (NOT a pipeline) contrasting three tiers of
knowledge-system capability for a wet-lab / biophysics audience:

    RAG  ->  TAG  ->  Contract-Graph (this system)

Self-contained: writes the PNG to a fixed absolute path and prints it.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---- shared style constants -------------------------------------------------
RAG = "#6B9AC4"      # muted blue
TAG = "#3C8D7A"      # teal
SYS = "#5B4B8A"      # indigo
EDGE = "#444"        # neutral box edge
RED = "#C0392B"      # failure / forbidden
GREEN = "#2E8B57"    # verified / pass
LIGHT = "#F5F6FA"    # light fill

OUT = (
    "/Users/sw1/ffn_cellsim-platform/ffn_sim/outputs/tag_kb/"
    "presentation/figs/three_tier_comparison.png"
)

# ---- table content ----------------------------------------------------------
COL_HEADERS = [
    ("RAG", RAG),
    ("TAG", TAG),
    ("Contract-Graph\n(this system)", SYS),
]

# Each row: (row label, [cell_RAG, cell_TAG, cell_SYS])
# A cell is (kind, text) where kind in {"x", "check", "text"}.
ROWS = [
    (
        "Knowledge form",
        [
            ("text", "text chunks"),
            ("text", "relational DB"),
            ("text", "curated\ncontract-graph"),
        ],
    ),
    (
        "Relational +\nmulti-hop joins",
        [("x", "No"), ("check", "Yes"), ("check", "Yes")],
    ),
    (
        "Provenance\nenforced",
        [("x", "No"), ("x", "No"), ("check", "Yes\n3 disk gates")],
    ),
    (
        "Stops a wrong number\nreaching the sim",
        [("x", "No"), ("x", "No"), ("check", "Yes")],
    ),
]

# ---- layout -----------------------------------------------------------------
LABEL_W = 3.0          # width of the row-label column
COL_W = 2.7            # width of each tier column
HEADER_H = 1.15        # height of the header row
ROW_H = 1.15           # height of each data row
GAP = 0.12             # gap between cells

n_rows = len(ROWS)
n_cols = len(COL_HEADERS)

fig_w = LABEL_W + n_cols * COL_W + 0.6
fig_h = HEADER_H + n_rows * ROW_H + 1.6

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
# pad the x-range so the outer columns aren't cropped flush by bbox_inches='tight'
ax.set_xlim(-0.35, LABEL_W + n_cols * COL_W + 0.35)
ax.set_ylim(-0.2, HEADER_H + n_rows * ROW_H + 0.1)
ax.set_axis_off()

top = HEADER_H + n_rows * ROW_H


def cell_xy(col_idx: int, row_idx: int):
    """Return lower-left (x, y) of the data cell at (col, row)."""
    x = LABEL_W + col_idx * COL_W
    y = top - HEADER_H - (row_idx + 1) * ROW_H
    return x, y


def draw_box(x, y, w, h, facecolor, edgecolor, lw=2):
    box = FancyBboxPatch(
        (x + GAP / 2, y + GAP / 2),
        w - GAP,
        h - GAP,
        boxstyle="round,pad=0.3,rounding_size=0.08",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        mutation_aspect=h / w if w else 1.0,
    )
    ax.add_patch(box)


# ---- column headers ---------------------------------------------------------
for c, (label, color) in enumerate(COL_HEADERS):
    x = LABEL_W + c * COL_W
    y = top - HEADER_H
    draw_box(x, y, COL_W, HEADER_H, facecolor=color, edgecolor=color, lw=2)
    ax.text(
        x + COL_W / 2,
        y + HEADER_H / 2,
        label,
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color="white",
    )

# little arrows between the three headers to convey RAG -> TAG -> system
hy = top - HEADER_H / 2
for c in range(n_cols - 1):
    x0 = LABEL_W + (c + 1) * COL_W - GAP / 2
    x1 = LABEL_W + (c + 1) * COL_W + GAP / 2
    ax.annotate(
        "",
        xy=(x1 + 0.02, hy),
        xytext=(x0 - 0.02, hy),
        annotation_clip=False,
        arrowprops=dict(arrowstyle="-|>", color=EDGE, lw=1.8),
    )

# ---- row labels + data cells ------------------------------------------------
for r, (row_label, cells) in enumerate(ROWS):
    y = top - HEADER_H - (r + 1) * ROW_H
    # row label (left column)
    draw_box(0, y, LABEL_W, ROW_H, facecolor=LIGHT, edgecolor=EDGE, lw=1.6)
    ax.text(
        LABEL_W / 2,
        y + ROW_H / 2,
        row_label,
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color="#222",
    )

    for c, (kind, text) in enumerate(cells):
        cx, cy = cell_xy(c, r)
        col_color = COL_HEADERS[c][1]
        draw_box(cx, cy, COL_W, ROW_H, facecolor=LIGHT, edgecolor=col_color, lw=2)
        ccx = cx + COL_W / 2
        ccy = cy + ROW_H / 2

        if kind == "text":
            ax.text(
                ccx,
                ccy,
                text,
                ha="center",
                va="center",
                fontsize=11,
                color="#222",
            )
        elif kind == "x":
            ax.text(
                ccx,
                ccy + 0.14,
                "✗",  # heavy ballot X
                ha="center",
                va="center",
                fontsize=24,
                fontweight="bold",
                color=RED,
            )
            ax.text(
                ccx,
                ccy - 0.34,
                text,
                ha="center",
                va="center",
                fontsize=10,
                color=RED,
            )
        elif kind == "check":
            ax.text(
                ccx,
                ccy + 0.14,
                "✓",  # check mark
                ha="center",
                va="center",
                fontsize=24,
                fontweight="bold",
                color=GREEN,
            )
            ax.text(
                ccx,
                ccy - 0.34,
                text,
                ha="center",
                va="center",
                fontsize=9.5,
                color=GREEN,
            )

# ---- title ------------------------------------------------------------------
fig.suptitle(
    "Why each step is needed: RAG → TAG → Contract-Graph",
    fontsize=15,
    fontweight="bold",
    y=0.99,
)
fig.text(
    0.5,
    0.935,
    "Each tier keeps the previous one's power and adds a capability the previous "
    "one structurally cannot provide",
    ha="center",
    fontsize=10,
    color="grey",
)

plt.tight_layout(rect=(0, 0, 1, 0.92))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print(OUT)
