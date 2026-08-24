"""Tier 1 — Conventional RAG mechanism figure (presentation set).

Renders a left-to-right RAG pipeline with a red callout box listing the three
limits of conventional RAG for a literature-anchored simulator. Self-contained:
writes the PNG to the absolute path below and prints it.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# Shared presentation palette
RAG = "#6B9AC4"
EDGE = "#444"
RED = "#C0392B"
LIGHT = "#F5F6FA"

OUT = (
    "/Users/sw1/ffn_cellsim-platform/aleph/outputs/tag_kb/"
    "presentation/figs/rag_mechanism.png"
)


def _box(ax, x, y, w, h, text, edgecolor):
    """Draw a rounded stage box centered text."""
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.3,rounding_size=0.08",
        facecolor=LIGHT,
        edgecolor=edgecolor,
        linewidth=2,
    )
    ax.add_patch(box)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=11,
        wrap=True,
    )


def _arrow(ax, x0, y0, x1, y1, label=None):
    arr = FancyArrowPatch(
        (x0, y0),
        (x1, y1),
        arrowstyle="-|>",
        mutation_scale=18,
        color=EDGE,
        lw=1.8,
    )
    ax.add_patch(arr)
    if label:
        ax.text(
            (x0 + x1) / 2,
            (y0 + y1) / 2 + 0.18,
            label,
            ha="center",
            va="bottom",
            fontsize=9,
            style="italic",
            color=EDGE,
        )


def main() -> str:
    fig, ax = plt.subplots(figsize=(13.5, 6.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.set_axis_off()

    # Pipeline stages: (text, width)
    stages = [
        "Documents\n(papers, notes)",
        "chunk\n+ embed",
        "vector\nstore",
        "top-k similarity\nretrieve",
        "stuff into\nLLM prompt",
        "answer",
    ]
    n = len(stages)
    bw = 12.5
    bh = 9.0
    gap = (100 - n * bw) / (n + 1)
    y = 30.0

    centers = []
    for i, text in enumerate(stages):
        x = gap + i * (bw + gap)
        _box(ax, x, y, bw, bh, text, RAG)
        centers.append((x, x + bw))

    arrow_labels = [None, None, "vectors", "query", "context", None]
    for i in range(n - 1):
        x0 = centers[i][1]
        x1 = centers[i + 1][0]
        _arrow(ax, x0, y + bh / 2, x1, y + bh / 2, arrow_labels[i + 1])

    # Red callout: the three limits
    cx, cy, cw, ch = 12, 4, 76, 16
    callout = FancyBboxPatch(
        (cx, cy),
        cw,
        ch,
        boxstyle="round,pad=0.3,rounding_size=0.08",
        facecolor="#FBEDEC",
        edgecolor=RED,
        linewidth=2,
    )
    ax.add_patch(callout)
    ax.text(
        cx + cw / 2,
        cy + ch - 3.2,
        "Limits for a literature-anchored simulator",
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="bold",
        color=RED,
    )
    limits = [
        "no joins / multi-hop",
        "no provenance guarantee",
        "can fabricate a citation",
    ]
    for j, lim in enumerate(limits):
        lx = cx + cw * (j + 0.5) / 3
        ax.text(
            lx,
            cy + 4.5,
            "✖  " + lim,
            ha="center",
            va="center",
            fontsize=10.5,
            color=RED,
        )

    # connector from pipeline down to callout
    _arrow(ax, 50, y - 0.5, 50, cy + ch + 0.3)

    fig.suptitle(
        "Tier 1 — Conventional RAG: similarity retrieval over text",
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )
    ax.text(
        50,
        47,
        "embed text → vector top-k → stuff context → generate",
        ha="center",
        va="center",
        fontsize=10,
        color="grey",
    )

    plt.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(OUT)
    return OUT


if __name__ == "__main__":
    main()
