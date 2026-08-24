#!/usr/bin/env python3
"""Render accepted AFM F-delta ensembles in the established static-figure style."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def render(record_path: Path, output: Path) -> None:
    """Plot mean reaction with between-seed SD; refuse anything but an accepted ensemble."""
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if (
        record.get("schema") != "afm-force-indentation-ensemble@1"
        or record.get("status") != "ACCEPTED_ENSEMBLE"
    ):
        raise ValueError("F-delta visualization requires an accepted AFM ensemble")
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    speeds = sorted({float(row["speed_um_s"]) for row in record["points"]})
    for speed in speeds:
        rows = [row for row in record["points"] if float(row["speed_um_s"]) == speed]
        depth = np.asarray([row["depth_um"] for row in rows], dtype=np.float64)
        mean = np.asarray([row["reaction_mean_pn"] for row in rows], dtype=np.float64)
        sd = np.asarray([row["reaction_sample_sd_pn"] for row in rows], dtype=np.float64)
        ax.plot(depth, mean, marker="o", linewidth=1.6, label=f"{speed:g} µm/s")
        ax.fill_between(depth, mean - sd, mean + sd, alpha=0.18)
    ax.set_xlabel("indentation depth δ [µm]")
    ax.set_ylabel("vertical indenter reaction F [pN]")
    ax.set_title("Accepted suspended-cell AFM force–indentation ensemble")
    ax.grid(alpha=0.22)
    ax.legend(title="loading speed", frameon=False)
    fig.text(
        0.5, 0.012,
        "Bands are sample SD across distinct seeds (≥3), not within-run SEM. Every plotted point passed "
        "the unchanged physical-step predicate.",
        ha="center", fontsize=7.5,
    )
    fig.tight_layout(rect=(0.0, 0.05, 1.0, 1.0))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    render(args.record, args.output)


if __name__ == "__main__":
    main()
