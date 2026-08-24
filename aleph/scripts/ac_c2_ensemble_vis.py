#!/usr/bin/env python3
"""Render the accepted full-native C-2 seed ensemble as a numerical gate figure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def render(record_path: Path, output: Path) -> None:
    """Plot each accepted terminal residual against its predeclared threshold."""
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if record.get("schema") != "afm-c2-accepted-ensemble@1":
        raise ValueError("unsupported C-2 ensemble schema")
    rows = record["records"]
    if len(rows) < 3:
        raise ValueError("accepted C-2 visualization requires at least three seed records")
    seeds = [int(row["seed"]) for row in rows]
    forces = [float(row["max_projected_force_pn"]) for row in rows]
    tolerances = [float(row["projected_force_tolerance_pn"]) for row in rows]
    if any(force > tolerance for force, tolerance in zip(forces, tolerances, strict=True)):
        raise ValueError("an accepted C-2 figure cannot contain a force-rejected record")

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    x = range(len(seeds))
    ax.bar(x, forces, width=0.58, color="#2166ac", alpha=0.88, label="accepted terminal maxPF")
    ax.scatter(x, tolerances, marker="_", s=430, linewidths=2.0, color="#b2182b",
               label="unchanged per-seed C-2 threshold")
    for index, (force, row) in enumerate(zip(forces, rows, strict=True)):
        ax.text(index, force - 0.00055, f"{force:.6f}\n{int(row['reported_iterations'])} iterations",
                ha="center", va="top", color="white", fontsize=8)
    ax.set_xticks(list(x), [f"seed {seed}" for seed in seeds])
    ax.set_ylabel("maximum projected nodal force [pN]")
    ax.set_ylim(0.0, max(tolerances) * 1.12)
    ax.set_title("Full-native C-2 accepted ensemble — membrane subdivision 8")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    fig.text(
        0.5,
        0.012,
        "RTX 3090 / Slurm; 1,150,806 nodes and all compartments. Numerical readiness only — not an "
        "F–δ curve or cortical-tension magnitude.",
        ha="center",
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0.0, 0.055, 1.0, 1.0))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.record, args.output)


if __name__ == "__main__":
    main()
