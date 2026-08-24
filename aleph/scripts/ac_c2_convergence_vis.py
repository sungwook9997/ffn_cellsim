#!/usr/bin/env python3
"""Render full-native C-2 projected-force histories without turning them into a material curve."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_trace(path: Path) -> tuple[str, list[int], list[float], float, str]:
    """Read one driver diagnostic and return only its numerical-solver trace."""
    record = json.loads(path.read_text(encoding="utf-8"))
    solver = record["solver"]
    history = solver["convergence_history"]
    label = f"{solver['name']} ({solver['iteration_budget']:,}-iteration budget)"
    objective = str(solver.get("line_search_objective", "max_force"))
    if objective != "max_force":
        label += f" [{objective} candidate selection]"
    return (
        label,
        [int(value) for value in history["iteration"]],
        [float(value) for value in history["projected_force_pN"]],
        float(solver["projected_force_tolerance_pN"]),
        str(record["status"]),
    )


def render(inputs: list[Path], output: Path) -> None:
    """Plot solver residual histories and the pre-existing acceptance threshold."""
    traces = [load_trace(path) for path in inputs]
    if not traces or any(not iterations for _, iterations, _, _, _ in traces):
        raise ValueError("each C-2 diagnostic must contain at least one convergence-history sample")
    tolerances = {trace[3] for trace in traces}
    if len(tolerances) != 1:
        raise ValueError("C-2 traces with different projected-force thresholds require separate figures")

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    for label, iterations, force, _, status in traces:
        verdict = "accepted" if status == "ACCEPTED" else "rejected, rolled back"
        ax.plot(iterations, force, marker="o", markersize=2.5, linewidth=1.3,
                label=f"{label} - {verdict}")
    tolerance = next(iter(tolerances))
    ax.axhline(tolerance, color="#b2182b", linestyle="--", linewidth=1.3,
               label=f"unchanged C-2 threshold ({tolerance:.3g} pN)")
    ax.set_yscale("log")
    ax.set_xlabel("inner mechanical iteration [-]")
    ax.set_ylabel("maximum projected nodal force [pN]")
    ax.set_title("Full-native C-2 numerical convergence diagnostic")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.text(
        0.5, 0.012,
        "Rejected candidates are rolled back. Lines connect solver diagnostics; no material response or "
        "convergence-rate fit.",
        ha="center", fontsize=7.5,
    )
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 1.0))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.inputs, args.output)


if __name__ == "__main__":
    main()
