"""Tier 2 — TAG mechanism figure for the KB-architecture talk.

Renders a self-contained schematic of Table-Augmented Generation:
    [NL question] -> [syn: LLM writes SQL] -> [exec: DuckDB over relational KB]
    -> [gen: answer + citations]
with a small relational-table icon feeding the exec stage and a green callout
naming what TAG adds over conventional RAG (exact joins, aggregation, multi-hop).
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# --- shared talk style -----------------------------------------------------
TAG = "#3C8D7A"        # TAG tier teal
EDGE = "#444"          # neutral box edge / arrows
GREEN = "#2E8B57"      # verified / new-capability accent
LIGHT = "#F5F6FA"      # light fill

OUT = "/Users/sw1/ffn_cellsim-platform/ffn_sim/outputs/tag_kb/presentation/figs/tag_mechanism.png"


def box(ax, x, y, w, h, text, *, edgecolor=TAG, fc=LIGHT, fs=11):
    """Rounded stage box with centered, wrapped label."""
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.3,rounding_size=0.08",
            facecolor=fc,
            edgecolor=edgecolor,
            linewidth=2,
            zorder=2,
        )
    )
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fs,
        wrap=True,
        zorder=3,
    )


def arrow(ax, x0, y0, x1, y1, label=None):
    ax.add_patch(
        FancyArrowPatch(
            (x0, y0),
            (x1, y1),
            arrowstyle="-|>",
            mutation_scale=18,
            color=EDGE,
            lw=1.8,
            zorder=1,
        )
    )
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


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    fig, ax = plt.subplots(figsize=(13, 6))
    ax.set_xlim(0, 26)
    ax.set_ylim(0, 12)
    ax.set_axis_off()

    # title + subtitle
    fig.suptitle(
        "Tier 2 — TAG: ask the database in natural language (Biswal 2024)",
        fontsize=15,
        fontweight="bold",
        y=0.97,
    )
    ax.text(
        13,
        11.0,
        "Table-Augmented Generation: syn (NL→SQL) → exec (DuckDB) → gen (answer + citations)",
        ha="center",
        va="center",
        fontsize=10,
        color="grey",
    )

    # --- main pipeline row (y centered ~6.5) ---
    yb, h = 5.6, 1.9
    box(ax, 0.4, yb, 4.0, h, "NL question\n(natural language)")
    box(ax, 6.0, yb, 4.2, h, "syn\nLLM writes SQL")
    box(ax, 11.8, yb, 4.6, h, "exec\nDuckDB over the\nrelational KB")
    box(ax, 18.0, yb, 4.4, h, "gen\nanswer + citations", edgecolor=GREEN)

    yc = yb + h / 2
    arrow(ax, 4.4, yc, 6.0, yc)
    arrow(ax, 10.2, yc, 11.8, yc, "SQL")
    arrow(ax, 16.4, yc, 18.0, yc, "rows")

    # --- relational table icon feeding exec ---
    tx, tw = 11.8, 4.6
    ty, th = 1.3, 2.6
    n_rows = 4
    rh = th / n_rows
    # outer frame
    ax.add_patch(
        FancyBboxPatch(
            (tx, ty),
            tw,
            th,
            boxstyle="round,pad=0.15,rounding_size=0.04",
            facecolor="white",
            edgecolor=TAG,
            linewidth=2,
            zorder=2,
        )
    )
    # header row (teal fill)
    ax.add_patch(
        plt.Rectangle((tx, ty + th - rh), tw, rh, facecolor=TAG, edgecolor=TAG, zorder=3)
    )
    ax.text(
        tx + tw / 2,
        ty + th - rh / 2,
        "Parameter  |  KU  |  source",
        ha="center",
        va="center",
        fontsize=8.5,
        color="white",
        fontweight="bold",
        zorder=4,
    )
    rows = [
        "k_off0     KU-4.2   Ferrer08",
        "γ_cortex   KU-3.5   Chugh17",
        "W_cs       KU-6.2   Nagle22",
    ]
    for i, r in enumerate(rows):
        ry = ty + th - rh * (i + 2)
        if i % 2 == 0:
            ax.add_patch(
                plt.Rectangle((tx, ry), tw, rh, facecolor=LIGHT, edgecolor="none", zorder=2)
            )
        ax.text(
            tx + 0.25,
            ry + rh / 2,
            r,
            ha="left",
            va="center",
            fontsize=8,
            family="monospace",
            zorder=4,
        )
        ax.plot([tx, tx + tw], [ry, ry], color="#ccc", lw=0.6, zorder=3)

    ax.text(
        tx + tw / 2,
        ty - 0.35,
        "relational Contract-Graph (8 DBs in DuckDB)",
        ha="center",
        va="top",
        fontsize=8.5,
        style="italic",
        color=EDGE,
    )
    # feed arrow from table up into exec
    arrow(ax, tx + tw / 2, ty + th, tx + tw / 2, yb)

    # --- green callout: NEW vs RAG ---
    ax.add_patch(
        FancyBboxPatch(
            (0.6, 0.6),
            10.4,
            3.2,
            boxstyle="round,pad=0.3,rounding_size=0.08",
            facecolor="#EAF4EE",
            edgecolor=GREEN,
            linewidth=2,
            zorder=2,
        )
    )
    ax.text(
        5.8,
        3.2,
        "NEW vs RAG",
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="bold",
        color=GREEN,
    )
    ax.text(
        5.8,
        1.85,
        "exact joins, aggregation, multi-hop\n— answers RAG structurally cannot",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#1d5b39",
    )

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(OUT)


if __name__ == "__main__":
    main()
