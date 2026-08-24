"""Developmental storyboard — compose the per-stage cell-morphology screenshots into ONE figure so the whole
motility DEVELOPMENT (resting → adhered → protruding → translocating → …) reads at a glance. Each panel is a
browser-verified screenshot of the interactive full-compartment viewer at that stage.

Usage:
  python -m aleph.scripts.ff_dev_storyboard --panel "S0 resting:.../S0.png" --panel "S1 adhered:.../S1.png" \
      --out .../ff_development_storyboard.png --title "FF single-cell motility development (native, A5000)"
"""
from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt


def build(panels: list[tuple[str, str]], out: str, title: str | None = None, caption: str | None = None):
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.0))
    if n == 1:
        axes = [axes]
    for ax, (label, path) in zip(axes, panels):
        try:
            ax.imshow(mpimg.imread(path))
        except Exception as e:                                 # missing panel → labelled placeholder, never crash
            ax.text(0.5, 0.5, f"(missing)\n{e}", ha="center", va="center", fontsize=8, wrap=True)
        ax.set_title(label, fontsize=12, fontweight="bold")
        ax.axis("off")
    if title:
        fig.suptitle(title, fontsize=15, y=0.99)
    if caption:
        fig.text(0.5, 0.01, caption, ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(out, dpi=110, bbox_inches="tight")
    return {"out": out, "panels": n}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", action="append", default=[], help="'label:path' — repeatable, left→right in dev order")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default=None)
    ap.add_argument("--caption", default=None)
    args = ap.parse_args()
    panels = [(p.split(":", 1)[0], p.split(":", 1)[1]) for p in args.panel]
    info = build(panels, args.out, title=args.title, caption=args.caption)
    print(f"wrote {info['out']}  ({info['panels']} panels)")


if __name__ == "__main__":
    main()
