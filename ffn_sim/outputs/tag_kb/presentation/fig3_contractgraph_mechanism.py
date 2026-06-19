"""Tier 3 — Contract-Graph + provenance gates (presentation figure).

Self-contained. Renders the THIS-SYSTEM tier: a Notion 8-DB Contract-Graph as
single source of truth, materialized to an Obsidian graph and a DuckDB+TAG
layer, with a lower band of three disk-grounded integrity gates wired into CI.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Shared palette
C_THIS = "#5B4B8A"   # indigo (this-system tier)
EDGE = "#444"
RED = "#C0392B"
GREEN = "#2E8B57"
LIGHT = "#F5F6FA"

BOX = dict(boxstyle="round,pad=0.3,rounding_size=0.08", facecolor=LIGHT,
           edgecolor=C_THIS, linewidth=2)

OUT = ("/Users/sw1/ffn_cellsim-platform/ffn_sim/outputs/tag_kb/"
       "presentation/figs/contractgraph_mechanism.png")


def box(ax, x, y, w, h, text, *, edgecolor=C_THIS, fontsize=11,
        fontweight="normal", textcolor="black"):
    """Draw a rounded box centered at (x, y) with wrapped label."""
    bb = dict(BOX)
    bb["edgecolor"] = edgecolor
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, **bb))
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            fontweight=fontweight, color=textcolor, wrap=True, zorder=5)


def arrow(ax, xy_from, xy_to, label=None, color=EDGE, lw=1.8):
    ax.add_patch(FancyArrowPatch(xy_from, xy_to, arrowstyle="-|>",
                                 mutation_scale=18, color=color, lw=lw,
                                 zorder=4))
    if label:
        mx = (xy_from[0] + xy_to[0]) / 2
        my = (xy_from[1] + xy_to[1]) / 2
        ax.text(mx + 0.4, my, label, fontsize=9, style="italic",
                color=color, ha="left", va="center")


def main() -> None:
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_axis_off()

    fig.suptitle("Tier 3 — Contract-Graph + provenance gates",
                 fontsize=15, fontweight="bold", y=0.97)
    ax.text(50, 95.5, "domain ontology (8-DB SoT) + disk-grounded integrity "
            "governance — beyond RAG and vanilla TAG",
            ha="center", va="center", fontsize=10, color="grey")

    # --- Top: Notion 8-DB Contract-Graph = SoT ---
    nodes = ("SourceEvidence  ·  KnowledgeClaim/KU  ·  ModelContract  ·  "
             "Parameter\n"
             "ValidationGate  ·  CodeMapping  ·  RunResult  ·  DecisionLedger")
    box(ax, 50, 86,
        78, 13,
        "Notion 8-DB Contract-Graph  =  Single Source of Truth\n\n" + nodes,
        edgecolor=C_THIS, fontsize=10.5, fontweight="bold")

    # --- materialize arrows down to two boxes ---
    arrow(ax, (38, 79), (28, 70), label=None)
    arrow(ax, (62, 79), (72, 70), label=None)
    ax.text(50, 75.5, "materialize  (refresh.sh)", ha="center", va="center",
            fontsize=9.5, style="italic", color=EDGE)

    box(ax, 26, 64, 38, 8,
        "Obsidian graph\nvisual neighbourhood browsing",
        edgecolor=C_THIS, fontsize=10)
    box(ax, 74, 64, 38, 8,
        "DuckDB + TAG\nrelational / multi-hop queries",
        edgecolor=C_THIS, fontsize=10)

    # --- arrow down into the integrity-gate band ---
    arrow(ax, (50, 60), (50, 50.5))

    # --- Integrity gate band ---
    # band backdrop
    ax.add_patch(FancyBboxPatch((6, 16), 88, 33,
                                boxstyle="round,pad=0.4,rounding_size=0.06",
                                facecolor="white", edgecolor=C_THIS,
                                linewidth=2.4, linestyle="--", zorder=1))
    ax.text(50, 46, "Integrity gates  (disk-grounded, run in CI on every push)",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color=C_THIS)

    gate_y = 35
    gw, gh = 27, 13
    box(ax, 22, gate_y, gw, gh,
        "verify_sources\n→ source_audit\n\ncitations real?\n(CrossRef DOI match)",
        edgecolor=GREEN, fontsize=9.5)
    box(ax, 50, gate_y, gw, gh,
        "verify_runs\n→ run_audit\n\nresults hold on disk?\n(sanctioned metric)",
        edgecolor=GREEN, fontsize=9.5)
    box(ax, 78, gate_y, gw, gh,
        "verify_params\n→ param_audit\n\nconstants traceable\nto an OK source?",
        edgecolor=GREEN, fontsize=9.5)

    ax.text(50, 27, "Gate = exit 1 on DRIFT  →  blocks the merge",
            ha="center", va="center", fontsize=9, style="italic", color=RED)

    # --- 5-hop provenance arrow spanning the three gates ---
    arrow(ax, (10, 21), (90, 21), color=C_THIS, lw=2.2)
    ax.text(50, 19, "5-hop provenance:  run_result → validation_gate → "
            "model_contract → parameter → source_evidence",
            ha="center", va="center", fontsize=9.5, style="italic",
            color=C_THIS)

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(OUT)


if __name__ == "__main__":
    main()
