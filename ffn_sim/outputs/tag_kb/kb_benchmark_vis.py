#!/usr/bin/env python3
"""Figures for the KB-architecture ablation benchmark (kb_benchmark.py).

Reads benchmark_results.json and renders paper-ready figures into figs/:
  fig_bench_headline.png   — accuracy vs hallucination across C1..C4 (the story)
  fig_bench_ragas.png      — RAGAS triad (faithfulness / relevancy / context recall)
  fig_bench_by_class.png   — accuracy by question class × condition (where each
                             capability switches on)
  fig_bench_heatmap.png    — per-question correctness & hallucination grid

Visualization-integrity rules (CLAUDE.md): no axis truncation (0–100 / 0–1
fixed), units/percent annotated, every bar value labelled, condition order is
the KB-evolution order C1→C4.

RUN:  conda activate ffn_sim && python kb_benchmark_vis.py
"""
from __future__ import annotations

import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).parent
RES = HERE / "benchmark_results.json"
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)

CONDS = ["C1", "C2", "C3", "C4"]
LABELS = {"C1": "C1\nno-KB\n(closed-book)", "C2": "C2\nRAG",
          "C3": "C3\nRAG+Obsidian", "C4": "C4\nRAG+TAG\n+Obsidian"}
# colourblind-safe sequential (light→dark = more KB structure)
COLORS = {"C1": "#bdbdbd", "C2": "#74a9cf", "C3": "#2b8cbe", "C4": "#045a8d"}


def _mean(xs):
    xs = [x for x in xs if isinstance(x, (int, float))]
    return sum(xs) / len(xs) if xs else float("nan")


def load():
    d = json.loads(RES.read_text())
    return d, d["results"]


def by_cond(rows, cond):
    return [r for r in rows if r["cond"] == cond]


# --------------------------------------------------------------------------- #
def fig_headline(d, rows):
    accs, halls, fact_halls = [], [], []
    for c in CONDS:
        rs = by_cond(rows, c)
        accs.append(100 * _mean([r["accuracy"] for r in rs if r["accuracy"] is not None]))
        halls.append(100 * _mean([r["hallucination"] for r in rs]))
        js = [r.get("judge", {}) for r in rs if r.get("judge")]
        tc = sum((j.get("n_atomic_claims") or 0) for j in js)
        th = sum((j.get("n_hallucinated") or 0) for j in js)
        fact_halls.append(100 * th / tc if tc else float("nan"))

    x = np.arange(len(CONDS))
    w = 0.26
    fig, ax = plt.subplots(figsize=(9, 5.5))
    b1 = ax.bar(x - w, accs, w, label="Accuracy (correct answers)", color="#1a9850")
    b2 = ax.bar(x, halls, w, label="Hallucination rate (programmatic)", color="#d73027")
    b3 = ax.bar(x + w, fact_halls, w, label="FActScore hallucinated-fact rate", color="#fc8d59")
    for bars in (b1, b2, b3):
        for b in bars:
            h = b.get_height()
            if not np.isnan(h):
                ax.annotate(f"{h:.0f}", (b.get_x() + b.get_width() / 2, h),
                            ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[c] for c in CONDS])
    ax.set_ylabel("percent (%)")
    ax.set_ylim(0, 105)
    ax.set_title("KB architecture vs answer quality — ffn_cellsim knowledge base\n"
                 f"answerer={d['model']}, judge=G-Eval, {d['n_questions']} gold questions",
                 fontsize=11)
    ax.legend(loc="upper center", fontsize=9, framealpha=0.9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_bench_headline.png", dpi=150)
    plt.close(fig)


def fig_ragas(d, rows):
    metrics = [("answer_relevancy", "Answer relevancy"),
               ("faithfulness", "Faithfulness (grounded-in-context)"),
               ("context_recall", "Context recall")]
    x = np.arange(len(CONDS))
    w = 0.25
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for i, (key, lab) in enumerate(metrics):
        vals = []
        for c in CONDS:
            js = [r.get("judge", {}) for r in by_cond(rows, c) if r.get("judge")]
            vals.append(_mean([j.get(key) for j in js]))
        bars = ax.bar(x + (i - 1) * w, vals, w, label=lab)
        for b in bars:
            h = b.get_height()
            if not np.isnan(h):
                ax.annotate(f"{h:.2f}", (b.get_x() + b.get_width() / 2, h),
                            ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[c] for c in CONDS])
    ax.set_ylabel("RAGAS score (0–1)")
    ax.set_ylim(0, 1.08)
    ax.set_title("RAGAS triad across KB architectures (LLM-judge)\n"
                 "context recall is 0 without retrieval (C1/C2); faithfulness stays "
                 "high as honest refusals count as faithful",
                 fontsize=11)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_bench_ragas.png", dpi=150)
    plt.close(fig)


def fig_by_class(d, rows):
    classes = sorted({r["cls"] for r in rows})
    data = np.full((len(classes), len(CONDS)), np.nan)
    for i, cls in enumerate(classes):
        for j, c in enumerate(CONDS):
            rs = [r for r in rows if r["cls"] == cls and r["cond"] == c]
            accs = [r["accuracy"] for r in rs if r["accuracy"] is not None]
            # citation_trap has no `accuracy`; use judge `correct` (caught fabrication)
            if not accs:
                accs = [1 if (r.get("judge") or {}).get("correct") else 0 for r in rs]
            data[i, j] = 100 * _mean(accs)

    fig, ax = plt.subplots(figsize=(8.5, 5))
    im = ax.imshow(data, cmap="RdYlGn", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(CONDS)))
    ax.set_xticklabels([c for c in CONDS])
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes)
    for i in range(len(classes)):
        for j in range(len(CONDS)):
            v = data[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        color="black", fontsize=9, fontweight="bold")
    ax.set_title("Correctness (%) by question class × KB architecture\n"
                 "relational/aggregation/citation-trap switch on only with TAG+graph",
                 fontsize=11)
    fig.colorbar(im, ax=ax, label="correct (%)")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_bench_by_class.png", dpi=150)
    plt.close(fig)


def fig_heatmap(d, rows):
    qids = sorted({r["qid"] for r in rows}, key=lambda s: int(s[1:]))
    # value: +1 correct, 0 wrong-but-honest, -1 hallucinated
    grid = np.zeros((len(qids), len(CONDS)))
    for i, q in enumerate(qids):
        for j, c in enumerate(CONDS):
            r = next((x for x in rows if x["qid"] == q and x["cond"] == c), None)
            if r is None:
                grid[i, j] = np.nan
                continue
            if r["hallucination"]:
                grid[i, j] = -1
            elif r["accuracy"] == 1 or (r.get("judge") or {}).get("correct"):
                grid[i, j] = 1
            else:
                grid[i, j] = 0

    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(["#d73027", "#ffffbf", "#1a9850"])  # halluc / honest-miss / correct
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(CONDS)))
    ax.set_xticklabels(CONDS)
    ax.set_yticks(range(len(qids)))
    cls_by_q = {r["qid"]: r["cls"] for r in rows}
    ax.set_yticklabels([f"{q} ({cls_by_q[q]})" for q in qids], fontsize=8)
    sym = {1: "✓", 0: "·", -1: "✗"}
    for i in range(len(qids)):
        for j in range(len(CONDS)):
            v = grid[i, j]
            if not np.isnan(v):
                ax.text(j, i, sym[int(v)], ha="center", va="center", fontsize=11)
    ax.set_title("Per-question outcome × KB architecture\n"
                 "green ✓ correct · yellow · honest miss · red ✗ hallucination",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_bench_heatmap.png", dpi=150)
    plt.close(fig)


def fig_model_compare():
    """Answerer-robustness panel: accuracy & FActScore-hallucination across
    conditions for two answerer models (primary vs weaker), if both exist."""
    haiku = HERE / "benchmark_results_haiku.json"
    if not haiku.exists():
        return False
    da = json.loads(RES.read_text())
    dh = json.loads(haiku.read_text())

    def agg(rows):
        acc, fh = [], []
        for c in CONDS:
            rs = by_cond(rows, c)
            acc.append(100 * _mean([r["accuracy"] for r in rs if r["accuracy"] is not None]))
            js = [r.get("judge", {}) for r in rs if r.get("judge")]
            tc = sum((j.get("n_atomic_claims") or 0) for j in js)
            th = sum((j.get("n_hallucinated") or 0) for j in js)
            fh.append(100 * th / tc if tc else float("nan"))
        return acc, fh

    a_acc, a_fh = agg(da["results"])
    h_acc, h_fh = agg(dh["results"])
    x = np.arange(len(CONDS))
    w = 0.2
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (sa, ha, ylab, ttl) in zip(
            (ax1, ax2),
            [((a_acc, h_acc), None, "accuracy (%)", "Accuracy"),
             ((a_fh, h_fh), None, "FActScore hallucinated-fact (%)", "Atomic hallucination")]):
        svals, hvals = sa
        b1 = ax.bar(x - w / 2, svals, w, label=f"answerer: {da['model']}", color="#2b8cbe")
        b2 = ax.bar(x + w / 2, hvals, w, label=f"answerer: {dh['model']}", color="#fdae61")
        for bars in (b1, b2):
            for b in bars:
                h = b.get_height()
                if not np.isnan(h):
                    ax.annotate(f"{h:.0f}", (b.get_x() + b.get_width() / 2, h),
                                ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels(CONDS)
        ax.set_ylabel(ylab); ax.set_title(ttl); ax.set_ylim(0, 105)
        ax.grid(axis="y", alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("Answerer-model robustness: KB-architecture effect holds across a "
                 "stronger and a weaker model (judge fixed = opus)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_bench_model_compare.png", dpi=150)
    plt.close(fig)
    return True


def main():
    if not RES.exists():
        raise SystemExit(f"{RES} not found — run kb_benchmark.py first.")
    d, rows = load()
    fig_headline(d, rows)
    fig_ragas(d, rows)
    fig_by_class(d, rows)
    fig_heatmap(d, rows)
    n = 4 + (1 if fig_model_compare() else 0)
    print(f"wrote {n} figures to {FIGS}/")


if __name__ == "__main__":
    main()
