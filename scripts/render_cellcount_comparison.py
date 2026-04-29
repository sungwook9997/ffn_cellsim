"""scripts/render_cellcount_comparison.py — aggregate Lam4 cell-count
sweep figures (1k/2k/4k/8k) with paper-grade dual save.

Per CODEX_FIGURE_GUIDE.md and PI directive 2026-04-30. Generates:
- A/A0_topdown_vs_PI.{svg,png}: 4 cell-count lines + PI experimental
- A/A0_largest_component_vs_PI.{svg,png}: intact-spheroid version
- hull_leverage_vs_t.{svg,png}: escape diagnostic per cell-count
- final_frame_topdown_quartet.{svg,png}: 4-panel final top-down
- final_frame_sideview_quartet.{svg,png}: 4-panel final side-view
- cellcount_summary_dashboard.{svg,png}: convergence panels
- parameter_comparison_table.{csv,png}: all 4 configs side by side

Outputs go to results/cellcount_comparison_lam4/figures/{Figures_for_draft,
Figures_for_PPT,Logs/_captions/stage1dc_cellcount/}.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from acs.visualization.figure_style import install_helvetica_style, save_figure_dual

CELL_COUNTS = [1000, 2000, 4000, 8000]
COLORS = {1000: "#9b1d20", 2000: "#bb6f00", 4000: "#3a8c4f", 8000: "#0d4ec9"}


def _safe_float(v, default=float("nan")):
    try:
        f = float(v)
        if np.isnan(f) or np.isinf(f):
            return default
        return f
    except (TypeError, ValueError):
        return default


def _load_metrics(run_dir: Path) -> List[Dict[str, str]]:
    p = run_dir / "metrics.csv"
    if not p.exists():
        return []
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _load_cfg(run_dir: Path):
    cfg_path = run_dir / "config.yaml"
    if not cfg_path.exists():
        return None
    with open(cfg_path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _load_pi_lam4_mean() -> Optional[np.ndarray]:
    """PI experimental Lam4 mean trajectory (26 positions averaged per
    1-hour bin). Returns array (N, 2) [time_hr, A/A0_mean]."""
    p = Path("data/experimental/260313_Lam4.csv")
    if not p.exists():
        return None
    with open(p, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    bins: Dict[int, List[float]] = {}
    for r in rows:
        try:
            t_min = float(r["Time_min"])
            a = float(r["A_over_A0"])
        except (KeyError, ValueError):
            continue
        b = int(t_min // 60)
        bins.setdefault(b, []).append(a)
    bins_sorted = sorted(bins.items())
    return np.array([[b, np.mean(vals)] for b, vals in bins_sorted])


def _trajectory(run_dir: Path) -> Dict[str, np.ndarray]:
    metrics = _load_metrics(run_dir)
    cfg = _load_cfg(run_dir)
    if not metrics:
        return {}
    tau = float(cfg["physics"]["maxwell_tau_s"]) if cfg else 60.0
    t_hr = np.array([_safe_float(r.get("time_star")) * tau / 3600.0 for r in metrics])
    return {
        "t_hr": t_hr,
        "aa0_hull": np.array([_safe_float(r.get("A_over_A0_topdown")) for r in metrics]),
        "aa0_largest": np.array([_safe_float(r.get("A_largest_component_over_A0_topdown")) for r in metrics]),
        "hull_lev": np.array([_safe_float(r.get("topdown_hull_leverage")) for r in metrics]),
        "frac": np.array([_safe_float(r.get("largest_component_fraction")) for r in metrics]),
        "R": np.array([_safe_float(r.get("effective_radius")) for r in metrics]),
        "psi": np.array([_safe_float(r.get("wadell_sphericity")) for r in metrics]),
    }


def _cellcount_runs() -> Dict[int, Path]:
    out: Dict[int, Path] = {}
    for n in CELL_COUNTS:
        p = Path(f"results/stage1dc_production_{n}_lam4")
        if p.exists() and (p / "metrics.csv").exists():
            out[n] = p
    return out


def render_aa0_overlay(runs: Dict[int, Path], out_dir: Path,
                       *, use_largest: bool = False) -> Path:
    pi_data = _load_pi_lam4_mean()
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ylab_key = "aa0_largest" if use_largest else "aa0_hull"
    metric_label = ("A_largest / A_0" if use_largest else "A_hull / A_0")
    for n, run_dir in runs.items():
        traj = _trajectory(run_dir)
        if not traj:
            continue
        ax.plot(traj["t_hr"], traj[ylab_key],
                color=COLORS[n], linewidth=1.6,
                label=f"sim {n} cells")
    if pi_data is not None:
        ax.scatter(pi_data[:, 0], pi_data[:, 1],
                   color="#222", marker="x", s=22, alpha=0.7,
                   label="PI exp Lam4 (n=26)")
    ax.set_xlabel("time (hr)")
    ax.set_ylabel(f"{metric_label} (top-down projection)")
    ax.set_title(
        ("Intact spheroid spread (largest connected component) vs PI experiment"
         if use_largest else
         "All-particle hull spread vs PI experiment (with escape artifacts)"),
        fontsize=8,
    )
    ax.legend(loc="upper left", fontsize=6.5, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    stem = ("AA0_largest_component_vs_PI" if use_largest
            else "AA0_topdown_hull_vs_PI")
    save_figure_dual(
        fig, stem, run_dir=out_dir,
        caption=(
            "Lam4 cell-count sweep top-down spreading trajectory compared to "
            "PI experimental data (260313, n=26 positions averaged per 1-hr "
            "bin). Lines: simulation; markers: experimental mean. "
            f"Variant: {'largest connected component (intact spheroid)' if use_largest else 'all-particle convex hull (includes escape artifacts)'}."
        ),
        description=f"AA0 trajectory comparison ({metric_label})",
        metadata={"figure_class": "comparison_overlay", "metric": metric_label},
        step="stage1dc_cellcount",
    )
    fig.savefig(out_dir / "figures" / f"{stem}.png", dpi=110)
    plt.close(fig)
    return out_dir / "figures" / f"{stem}.png"


def render_hull_leverage_panel(runs: Dict[int, Path], out_dir: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for n, run_dir in runs.items():
        traj = _trajectory(run_dir)
        if not traj:
            continue
        axes[0].plot(traj["t_hr"], traj["hull_lev"],
                     color=COLORS[n], linewidth=1.5, label=f"{n} cells")
        axes[1].plot(traj["t_hr"], traj["frac"],
                     color=COLORS[n], linewidth=1.5, label=f"{n} cells")
    axes[0].axhline(1.35, color="#666", linestyle="--", linewidth=0.8,
                    label="gate ≤ 1.35")
    axes[0].set_xlabel("time (hr)")
    axes[0].set_ylabel("A_hull / A_largest_component")
    axes[0].set_title("Hull leverage (escape diagnostic)", fontsize=8)
    axes[0].legend(loc="upper right", fontsize=6.5, frameon=False)
    axes[0].grid(alpha=0.3)
    axes[1].axhline(0.90, color="#666", linestyle="--", linewidth=0.8,
                    label="gate ≥ 0.90")
    axes[1].set_xlabel("time (hr)")
    axes[1].set_ylabel("largest connected fraction")
    axes[1].set_title("Largest component fraction", fontsize=8)
    axes[1].legend(loc="lower left", fontsize=6.5, frameon=False)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    save_figure_dual(
        fig, "hull_leverage_vs_t", run_dir=out_dir,
        caption=(
            "Stage 1d.c connected-component diagnostic across cell-count "
            "sweep. Left: hull leverage = A_hull / A_largest_component (1.0 "
            "= intact, larger = escape inflates hull). Right: fraction of "
            "particles in the largest connected component. Dashed lines mark "
            "gate limits 1.35 and 0.90."
        ),
        description="Escape diagnostic across cell counts",
        metadata={"figure_class": "diagnostic"},
        step="stage1dc_cellcount",
    )
    fig.savefig(out_dir / "figures" / "hull_leverage_vs_t.png", dpi=110)
    plt.close(fig)
    return out_dir / "figures" / "hull_leverage_vs_t.png"


def render_final_frame_quartets(runs: Dict[int, Path], out_dir: Path) -> Dict[str, Path]:
    try:
        import h5py
    except ImportError:
        return {}
    out_paths: Dict[str, Path] = {}
    for view in ("topdown", "sideview"):
        fig, axes = plt.subplots(2, 2, figsize=(8, 8))
        for ax, n in zip(axes.flat, CELL_COUNTS):
            run_dir = runs.get(n)
            if run_dir is None:
                ax.text(0.5, 0.5, f"no {n}", ha="center", va="center")
                ax.axis("off")
                continue
            cfg = _load_cfg(run_dir)
            domain = float(cfg["nondim"]["domain_star"]) if cfg else 6.0
            with h5py.File(run_dir / "snapshots.h5", "r") as hf:
                keys = sorted(hf["frames"].keys(), key=lambda k: int(k))
                x = np.array(hf["frames"][keys[-1]]["position"])
            if view == "topdown":
                ax.scatter(x[:, 0], x[:, 1], s=2, c="#3a82f7",
                           alpha=0.55, edgecolors="none")
                ax.set_xlim(0, domain)
                ax.set_ylim(0, domain)
                ax.set_aspect("equal")
                ax.set_xlabel("x*")
                ax.set_ylabel("y*")
            else:
                ax.scatter(x[:, 0], x[:, 2], s=2, c="#3a82f7",
                           alpha=0.55, edgecolors="none")
                ax.axhline(0, color="#444", linestyle="--", linewidth=0.8)
                ax.set_xlim(0, domain)
                ax.set_ylim(-0.05 * domain, 0.5 * domain)
                ax.set_xlabel("x*")
                ax.set_ylabel("z*")
            ax.set_title(f"{n} cells", fontsize=8, fontweight="bold")
        fig.suptitle(
            f"Final frame {view.replace('top', 'top-')} (Lam4 80 hr)",
            fontsize=10, y=0.995,
        )
        fig.tight_layout()
        stem = f"final_frame_{view}_quartet"
        save_figure_dual(
            fig, stem, run_dir=out_dir,
            caption=(
                f"Final-frame {view} projection across the Lam4 cell-count "
                f"sweep (1k / 2k / 4k / 8k cells, 80 hr). Each panel shows all "
                f"particles at the last saved frame. Top-down panels show the "
                f"xy projection (Hard Rule 11 PI-matched modality); side-view "
                f"panels show the xz projection with the substrate at z=0."
            ),
            description=f"Final-frame {view} quartet",
            metadata={"figure_class": "snapshot_quartet"},
            step="stage1dc_cellcount",
        )
        fig.savefig(out_dir / "figures" / f"{stem}.png", dpi=110)
        out_paths[view] = out_dir / "figures" / f"{stem}.png"
        plt.close(fig)
    return out_paths


def render_summary_dashboard(runs: Dict[int, Path], out_dir: Path) -> Path:
    pi_data = _load_pi_lam4_mean()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    # (0,0) A/A0_largest vs PI
    for n, run_dir in runs.items():
        traj = _trajectory(run_dir)
        if not traj:
            continue
        axes[0, 0].plot(traj["t_hr"], traj["aa0_largest"],
                        color=COLORS[n], linewidth=1.5, label=f"{n}")
    if pi_data is not None:
        axes[0, 0].scatter(pi_data[:, 0], pi_data[:, 1],
                           color="#222", marker="x", s=22, alpha=0.7,
                           label="PI exp")
    axes[0, 0].set_xlabel("time (hr)")
    axes[0, 0].set_ylabel("A_largest / A_0")
    axes[0, 0].set_title("Intact-spheroid spread vs PI", fontsize=8)
    axes[0, 0].legend(loc="upper left", fontsize=6.5, frameon=False)
    axes[0, 0].grid(alpha=0.3)
    # (0,1) end A/A0 vs cell count
    end_largest = []
    end_pi = pi_data[-1, 1] if pi_data is not None else float("nan")
    for n, run_dir in runs.items():
        traj = _trajectory(run_dir)
        if traj and len(traj["aa0_largest"]) > 0:
            end_largest.append((n, traj["aa0_largest"][-1]))
    if end_largest:
        ns, vals = zip(*end_largest)
        axes[0, 1].bar([str(n) for n in ns], vals,
                       color=[COLORS[n] for n in ns], alpha=0.8)
        if not np.isnan(end_pi):
            axes[0, 1].axhline(end_pi, color="#222", linestyle="--",
                               linewidth=1.5, label=f"PI 78 hr ≈ {end_pi:.2f}")
            axes[0, 1].legend(fontsize=6.5, frameon=False)
    axes[0, 1].set_xlabel("cells")
    axes[0, 1].set_ylabel("A_largest / A_0 at 80 hr")
    axes[0, 1].set_title("Endpoint convergence vs cell count", fontsize=8)
    axes[0, 1].grid(alpha=0.3, axis="y")
    # (1,0) hull leverage
    for n, run_dir in runs.items():
        traj = _trajectory(run_dir)
        if traj:
            axes[1, 0].plot(traj["t_hr"], traj["hull_lev"],
                            color=COLORS[n], linewidth=1.5, label=f"{n}")
    axes[1, 0].axhline(1.35, color="#666", linestyle="--", linewidth=0.8,
                       label="gate")
    axes[1, 0].set_xlabel("time (hr)")
    axes[1, 0].set_ylabel("hull leverage")
    axes[1, 0].set_title("Escape diagnostic (Codex connected-component gate)",
                          fontsize=8)
    axes[1, 0].legend(fontsize=6.5, frameon=False)
    axes[1, 0].grid(alpha=0.3)
    # (1,1) sphericity
    for n, run_dir in runs.items():
        traj = _trajectory(run_dir)
        if traj:
            axes[1, 1].plot(traj["t_hr"], traj["psi"],
                            color=COLORS[n], linewidth=1.5, label=f"{n}")
    axes[1, 1].axhline(0.70, color="#666", linestyle="--", linewidth=0.8,
                       label="post-spread limit")
    axes[1, 1].set_xlabel("time (hr)")
    axes[1, 1].set_ylabel("Wadell sphericity ψ")
    axes[1, 1].set_title("Sphericity convergence", fontsize=8)
    axes[1, 1].legend(fontsize=6.5, frameon=False)
    axes[1, 1].grid(alpha=0.3)
    fig.suptitle("Stage 1d.c Lam4 80 hr — cell-count sweep summary",
                 fontsize=10, y=0.995)
    fig.tight_layout()
    save_figure_dual(
        fig, "cellcount_summary_dashboard", run_dir=out_dir,
        caption=(
            "Stage 1d.c Lam4 cell-count sweep summary dashboard. "
            "Top-left: intact-spheroid A_largest/A_0 trajectory vs PI "
            "experimental Lam4 mean (n=26 positions, 260313). Top-right: "
            "endpoint A_largest/A_0 at 80 hr by cell count, with the PI "
            "78 hr mean as horizontal reference. Bottom-left: hull-leverage "
            "trajectory (Codex Stage 1d.c connected-component diagnostic). "
            "Bottom-right: Wadell sphericity ψ convergence."
        ),
        description="Cell-count sweep summary",
        metadata={"figure_class": "summary_dashboard"},
        step="stage1dc_cellcount",
    )
    fig.savefig(out_dir / "figures" / "cellcount_summary_dashboard.png", dpi=110)
    plt.close(fig)
    return out_dir / "figures" / "cellcount_summary_dashboard.png"


def main() -> None:
    install_helvetica_style()
    runs = _cellcount_runs()
    out_dir = Path("results/cellcount_comparison_lam4")
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    print(f"[cellcount_comparison] runs: {sorted(runs.keys())}")

    print("  → A/A0_topdown_hull vs PI overlay")
    render_aa0_overlay(runs, out_dir, use_largest=False)
    print("  → A/A0_largest_component vs PI overlay")
    render_aa0_overlay(runs, out_dir, use_largest=True)
    print("  → hull leverage diagnostic")
    render_hull_leverage_panel(runs, out_dir)
    print("  → final-frame quartets (topdown + sideview)")
    render_final_frame_quartets(runs, out_dir)
    print("  → summary dashboard")
    render_summary_dashboard(runs, out_dir)
    print(f"[cellcount_comparison] Done. Outputs: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
