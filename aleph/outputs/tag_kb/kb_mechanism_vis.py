#!/usr/bin/env python3
"""Mechanism / architecture figures for the ffn_cellsim KB (paper 'method' figures).

These are schematic diagrams of HOW the system works (not benchmark results):
  fig_mechanism_architecture.png — Notion SoT -> materialised Obsidian graph +
        DuckDB table/content layer -> TAG engine (the data-flow).
  fig_mechanism_tag_loop.png     — the canonical TAG loop syn->exec->gen with the
        BM25 PDF-content branch (the query mechanism).
  fig_mechanism_ablation.png     — the 4 benchmark conditions as strictly-additive
        context layers feeding one fixed answerer (C1 ⊂ C2 ⊂ C3 ⊂ C4).

Pure schematic — no data dependency.
RUN:  conda activate ffn_sim && python kb_mechanism_vis.py
"""
from __future__ import annotations

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = pathlib.Path(__file__).parent
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)

C_SOT = "#b30000"      # source of truth
C_TBL = "#045a8d"      # duckdb / tables
C_GRAPH = "#2b8cbe"    # obsidian graph
C_PDF = "#807dba"      # pdf corpus
C_ENGINE = "#238b45"   # tag engine
C_BOX = "#f7f7f7"


def box(ax, x, y, w, h, text, fc=C_BOX, ec="#333", fs=9, tc="#111", bold=False, round=0.02):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.004,rounding_size={round}",
                       linewidth=1.4, edgecolor=ec, facecolor=fc, zorder=2)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, zorder=3, fontweight="bold" if bold else "normal", wrap=True)


def arrow(ax, x1, y1, x2, y2, text="", color="#444", style="-|>", lw=1.8, ls="-",
          rad=0.0, fs=8):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=16,
                        linewidth=lw, color=color, linestyle=ls,
                        connectionstyle=f"arc3,rad={rad}", zorder=1)
    ax.add_patch(a)
    if text:
        ax.text((x1 + x2) / 2, (y1 + y2) / 2, text, ha="center", va="center",
                fontsize=fs, color=color, style="italic",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))


def _ax(figsize):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    return fig, ax


# --------------------------------------------------------------------------- #
def fig_architecture():
    fig, ax = _ax((12, 7))
    ax.set_title("ffn_cellsim knowledge-base architecture — one source of truth, two read layers, one query engine",
                 fontsize=12, fontweight="bold")
    # SoT
    box(ax, 3, 60, 30, 30,
        "Notion Contract-Graph\n(SOURCE OF TRUTH)\n\n8 atomic, bidirectionally-\nrelated databases:\nSourceEvidence · KnowledgeClaim\nModelContract · Parameter\nValidationGate · CodeMapping\nRunResult · DecisionLedger",
        fc="#fde0dd", ec=C_SOT, fs=8.5, bold=False)
    ax.text(18, 91.5, "edit here ONLY", ha="center", fontsize=8, color=C_SOT, style="italic")
    # reference PDFs
    box(ax, 3, 30, 30, 18,
        "references/  (≈128 PDFs)\nspheroid · cancer · cortex ·\nadhesion · ECM literature",
        fc="#efedf5", ec=C_PDF, fs=8.5)

    # materialize arrows
    arrow(ax, 33, 78, 50, 84, "materialize\n(Notion API)", color=C_GRAPH, rad=0.12)
    arrow(ax, 33, 70, 50, 60, "materialize\n(notion_to_duckdb)", color=C_TBL, rad=-0.12)
    arrow(ax, 33, 38, 50, 47, "fitz chunk\n(references_ingest)", color=C_PDF, rad=-0.1)

    # Obsidian
    box(ax, 50, 78, 30, 14,
        "Obsidian vault  [GRAPH layer]\n~640 nodes · typed wikilinks ·\ncolor groups · Dev-Log nodes\n→ visual graph exploration",
        fc="#d0e1f2", ec=C_GRAPH, fs=8.2)
    # DuckDB
    box(ax, 50, 44, 34, 26,
        "kb.duckdb  [TABLE + CONTENT layer]\n\n8 node tables (rows = Notion rows)\nedges(808)  ← relation triples\nreferences(128) · paper_chunks(5327)\n+ DuckDB BM25 full-text index",
        fc="#cfe8f3", ec=C_TBL, fs=8.2)

    # TAG engine
    box(ax, 64, 12, 30, 22,
        "TAG engine  (tag_query.py)\n\nsyn → exec → gen\nNL question → DuckDB SQL\n→ rows → grounded answer\n+ citations",
        fc="#c7e9c0", ec=C_ENGINE, fs=8.5, bold=True)
    arrow(ax, 67, 44, 75, 34, "schema + edges\nvocab + chunks", color=C_ENGINE, rad=-0.1)

    # refresh note
    box(ax, 50, 30, 34, 9,
        "refresh.sh  — both layers are REGENERATED from Notion\n(never hand-edited; Notion stays canonical)",
        fc="#ffffe0", ec="#999", fs=8)
    ax.text(50, 6, "Which layer? See structure / follow links → Obsidian   ·   exact relational/aggregate answer + citations → TAG",
            ha="center", fontsize=8.5, color="#333", style="italic")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_mechanism_architecture.png", dpi=150)
    plt.close(fig)


def fig_tag_loop():
    fig, ax = _ax((12, 5.2))
    ax.set_title("The TAG loop (Biswal et al. 2024): syn → exec → gen — exact SQL fused with PDF-content semantics",
                 fontsize=11.5, fontweight="bold")
    y = 60
    box(ax, 2, y, 13, 16, "NL question R\n\n“which Params feed\na failing gate?”", fc="#f0f0f0", fs=8.5)
    box(ax, 20, y, 17, 16, "syn  (Claude)\n\nlive schema +\nedges relation vocab\n→ one DuckDB SELECT", fc="#c7e9c0", ec=C_ENGINE, fs=8.2)
    box(ax, 42, y, 17, 16, "exec  (DuckDB)\n\nread-only · SELECT-guard\n1× self-repair on error\n→ result rows T", fc="#cfe8f3", ec=C_TBL, fs=8.2)
    box(ax, 64, y, 16, 16, "gen  (Claude)\n\n(R, T, excerpts)\n→ answer A\n+ citations", fc="#c7e9c0", ec=C_ENGINE, fs=8.2)
    box(ax, 84, y, 14, 16, "Answer A\n\ngrounded in rows;\ncites KB-id / DOI /\ncitation_key", fc="#f0f0f0", fs=8.2)
    arrow(ax, 15, y + 8, 20, y + 8)
    arrow(ax, 37, y + 8, 42, y + 8, "SQL Q")
    arrow(ax, 59, y + 8, 64, y + 8, "rows T")
    arrow(ax, 80, y + 8, 84, y + 8)
    arrow(ax, 50, y, 30, y, "DuckDB error → repair", color="#d95f0e", ls="--", rad=-0.4, fs=7.5)
    # content branch
    box(ax, 30, 24, 30, 14, "paper_chunks  (5327)\nBM25 / full-text index\n→ top-k evidence excerpts",
        fc="#efedf5", ec=C_PDF, fs=8.2)
    arrow(ax, 45, 38, 70, y, "excerpts (semantic step)", color=C_PDF, rad=0.18)
    ax.text(50, 12, "Text2SQL alone can't do the semantic step; RAG alone can't do the exact join/aggregate — TAG unifies both.",
            ha="center", fontsize=8.5, color="#333", style="italic")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_mechanism_tag_loop.png", dpi=150)
    plt.close(fig)


def fig_ablation():
    fig, ax = _ax((12, 6.6))
    ax.set_title("Benchmark conditions = strictly-additive context layers (C1 ⊂ C2 ⊂ C3 ⊂ C4), one fixed answerer",
                 fontsize=11.5, fontweight="bold")
    conds = [
        ("C1  no-KB", ["(closed-book —", "parametric memory only)"], "#bdbdbd",
         "ActiveCellSim era"),
        ("C2  RAG", ["+ BM25 passages over", "the reference-PDF corpus"], "#74a9cf",
         "unstructured text"),
        ("C3  +Obsidian", ["+ Contract-graph nodes", "& their typed edges"], "#3690c0",
         "graph structure"),
        ("C4  +TAG", ["+ synthesised SQL result", "+ source_audit table"], "#045a8d",
         "exact relational query"),
    ]
    x0, w, gap = 4, 21, 2
    for i, (name, lines, col, tag) in enumerate(conds):
        x = x0 + i * (w + gap)
        # cumulative stack: each condition shows all blocks up to its level
        blocks = []
        if i >= 0: blocks.append(("closed-book LLM", "#eeeeee"))
        if i >= 1: blocks.append(("RAG passages", "#deebf7"))
        if i >= 2: blocks.append(("graph nodes+edges", "#c6dbef"))
        if i >= 3: blocks.append(("TAG SQL + audit", "#c7e9c0"))
        by = 20
        bh = 9
        for label, bc in blocks:
            box(ax, x, by, w, bh, label, fc=bc, fs=8, ec="#888")
            by += bh + 1.5
        # header
        box(ax, x, by + 1, w, 9, name, fc=col, ec="#222", fs=10, bold=True,
            tc="white" if i >= 2 else "#111")
        ax.text(x + w / 2, by + 0.0, tag, ha="center", va="top", fontsize=7.5,
                color="#444", style="italic")
        ax.text(x + w / 2, 15, "↓ same answerer\n(tools OFF, scratch cwd)", ha="center",
                va="top", fontsize=7.5, color="#333")
    # the shared answerer bar
    box(ax, 4, 4, 4 * w + 3 * gap, 8,
        "FIXED ANSWERER  (claude-sonnet-4-6 / -haiku-4-5)  ·  identical prompt  ·  only the context block above differs  →  judged by claude-opus-4-8 (G-Eval)",
        fc="#fff7bc", ec="#cc9", fs=8.5, bold=True)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_mechanism_ablation.png", dpi=150)
    plt.close(fig)


def main():
    fig_architecture()
    fig_tag_loop()
    fig_ablation()
    print(f"wrote 3 mechanism figures to {FIGS}/")


if __name__ == "__main__":
    main()
