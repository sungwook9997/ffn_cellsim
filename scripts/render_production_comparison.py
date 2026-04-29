"""scripts/render_production_comparison.py — aggregate Bare/Pre/Lam4
production comparison figures.

Per PI visualization directive 2026-04-29.

Reads multiple run directories and produces the cross-condition
package under results/production_comparison/figures/:
- phenotype_triptych.png         (final top-down 3-panel: Bare/Pre/Lam4)
- final_sideview_triptych.png    (final side-view 3-panel)
- AA0_experiment_overlay.png     (sim A/A₀_topdown(t) + PI experimental)
- production_summary_dashboard.png
- parameter_comparison_table.csv + .png

Usage:
    python -m scripts.render_production_comparison \
        --bare results/<bare_run> \
        --pre results/<pre_run> \
        --lam4 results/<lam4_run> \
        [--out results/production_comparison]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml


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


def _load_pi_experimental(condition: str) -> Optional[np.ndarray]:
    """Load PI experimental Area_um2 column for a given condition.

    Looks under data/experimental/260313_<condition>.csv (file naming
    matches PI 2026-03-13 experimental dataset). Returns array of
    shape (N, 2) with columns [time_hr, A/A0] or None.
    """
    candidates = [
        Path(f"data/experimental/260313_{condition}.csv"),
        Path(f"data/experimental/{condition}.csv"),
    ]
    for p in candidates:
        if p.exists():
            try:
                with open(p, encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                if not rows:
                    continue
                # Try common column patterns.
                t_keys = ["Time_hr", "time_hr", "time", "Hr", "T"]
                a_keys = ["Area_um2", "A_over_A0", "Area", "Area_norm"]
                t_col = next((k for k in t_keys if k in rows[0]), None)
                a_col = next((k for k in a_keys if k in rows[0]), None)
                if not t_col or not a_col:
                    continue
                t_arr = np.array([_safe_float(r[t_col]) for r in rows])
                a_arr = np.array([_safe_float(r[a_col]) for r in rows])
                # Normalize to A/A0 if raw area.
                a0 = a_arr[0] if a_arr.size > 0 and not np.isnan(a_arr[0]) else 1.0
                if a0 > 1.5:  # raw μm² → normalize
                    a_arr = a_arr / a0
                return np.column_stack([t_arr, a_arr])
            except Exception:
                continue
    return None


def _condition_label(cfg) -> str:
    if cfg is None:
        return "?"
    phi = cfg.get("layer3", {}).get("phi_initial", None)
    if phi is None:
        return "?"
    if phi <= 0.40:
        return "Bare"
    if phi <= 0.65:
        return "Pre"
    return "Lam4"


def _t_hr(metrics, cfg):
    tau = float(cfg["physics"]["maxwell_tau_s"]) if cfg else 60.0
    t_star = np.array([_safe_float(r.get("time_star", 0)) for r in metrics])
    return t_star * tau / 3600.0


def render_triptych_topdown(runs: Dict[str, Path], out_path: Path) -> None:
    """3-panel final top-down comparison."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import h5py
    except ImportError:
        return
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (label, run_dir) in zip(axes, runs.items()):
        snaps_path = run_dir / "snapshots.h5"
        cfg = _load_cfg(run_dir)
        domain = float(cfg["nondim"]["domain_star"]) if cfg else 6.0
        if not snaps_path.exists():
            ax.text(0.5, 0.5, f"no snapshots ({label})",
                    ha="center", va="center")
            ax.axis("off")
            continue
        with h5py.File(snaps_path, "r") as hf:
            keys = sorted([k for k in hf.keys() if k.startswith("frame_")],
                          key=lambda s: int(s.split("_")[1]))
            if not keys:
                ax.text(0.5, 0.5, "empty", ha="center", va="center")
                ax.axis("off")
                continue
            x = np.array(hf[keys[-1]]["x"])
        ax.scatter(x[:, 0], x[:, 1], s=4, c="#3a82f7", alpha=0.6, edgecolors="none")
        try:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(x[:, :2])
            poly = np.vstack([x[:, :2][hull.vertices], x[:, :2][hull.vertices][:1]])
            ax.plot(poly[:, 0], poly[:, 1], color="#0d4ec9", linewidth=1.5)
        except Exception:
            pass
        ax.set_xlim(0, domain)
        ax.set_ylim(0, domain)
        ax.set_aspect("equal")
        ax.set_title(label, fontsize=14, fontweight="bold")
        ax.set_xlabel("x*")
        ax.set_ylabel("y*")
    fig.suptitle("Phenotype final top-down comparison (Bare / Pre / Lam4)",
                 fontsize=14, y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def render_triptych_sideview(runs: Dict[str, Path], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import h5py
    except ImportError:
        return
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (label, run_dir) in zip(axes, runs.items()):
        snaps_path = run_dir / "snapshots.h5"
        cfg = _load_cfg(run_dir)
        domain = float(cfg["nondim"]["domain_star"]) if cfg else 6.0
        if not snaps_path.exists():
            ax.text(0.5, 0.5, f"no snapshots ({label})",
                    ha="center", va="center")
            ax.axis("off")
            continue
        with h5py.File(snaps_path, "r") as hf:
            keys = sorted([k for k in hf.keys() if k.startswith("frame_")],
                          key=lambda s: int(s.split("_")[1]))
            x = np.array(hf[keys[-1]]["x"])
        ax.scatter(x[:, 0], x[:, 2], s=4, c="#3a82f7", alpha=0.6, edgecolors="none")
        ax.axhline(0.0, color="#444", linestyle="--", linewidth=1.0)
        ax.set_xlim(0, domain)
        ax.set_ylim(-0.05 * domain, 0.5 * domain)
        ax.set_title(label, fontsize=14, fontweight="bold")
        ax.set_xlabel("x*")
        ax.set_ylabel("z*")
    fig.suptitle("Phenotype final side-view comparison",
                 fontsize=14, y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def render_aa0_experiment_overlay(runs: Dict[str, Path], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = {"Bare": "#9b1d20", "Pre": "#bb6f00", "Lam4": "#0d4ec9"}
    for label, run_dir in runs.items():
        cfg = _load_cfg(run_dir)
        metrics = _load_metrics(run_dir)
        if not metrics:
            continue
        t = _t_hr(metrics, cfg)
        a = np.array([_safe_float(r.get("A_over_A0_topdown", 1.0)) for r in metrics])
        col = colors.get(label, "#444")
        ax.plot(t, a, color=col, linewidth=2.0, label=f"sim {label}")
        # Experimental overlay
        exp = _load_pi_experimental(label)
        if exp is not None:
            ax.scatter(exp[:, 0], exp[:, 1], color=col, marker="x", s=36,
                       label=f"exp {label} (PI 2026-03-13)")
    ax.set_xlabel("t (hr)")
    ax.set_ylabel("A/A₀_topdown")
    ax.set_title("Sim vs PI experimental — top-down projection (PI is comparison only, NOT fit)")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


def render_summary_dashboard(runs: Dict[str, Path], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    colors = {"Bare": "#9b1d20", "Pre": "#bb6f00", "Lam4": "#0d4ec9"}
    for label, run_dir in runs.items():
        cfg = _load_cfg(run_dir)
        metrics = _load_metrics(run_dir)
        if not metrics:
            continue
        t_h = _t_hr(metrics, cfg)
        a = np.array([_safe_float(r.get("A_over_A0_topdown", 1.0)) for r in metrics])
        R = np.array([_safe_float(r.get("effective_radius", 1.0)) for r in metrics])
        psi = np.array([_safe_float(r.get("wadell_sphericity", 1.0)) for r in metrics])
        rho = np.array([_safe_float(r.get("rho_osm_max", 1.0)) for r in metrics])
        col = colors.get(label, "#444")
        axes[0, 0].plot(t_h, a, color=col, label=label)
        axes[0, 1].plot(t_h, R, color=col, label=label)
        axes[1, 0].plot(t_h, psi, color=col, label=label)
        axes[1, 1].plot(t_h, rho, color=col, label=label)
    axes[0, 0].set_title("A/A₀_topdown")
    axes[0, 0].set_xlabel("t (hr)")
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].legend(fontsize=9)
    axes[0, 1].set_title("R/R₀ (effective radius)")
    axes[0, 1].set_xlabel("t (hr)")
    axes[0, 1].grid(alpha=0.3)
    axes[1, 0].set_title("ψ Wadell sphericity")
    axes[1, 0].set_xlabel("t (hr)")
    axes[1, 0].grid(alpha=0.3)
    axes[1, 1].set_title("ρ_osm max")
    axes[1, 1].set_xlabel("t (hr)")
    axes[1, 1].grid(alpha=0.3)
    fig.suptitle("Production summary dashboard (Bare / Pre / Lam4)",
                 fontsize=14, y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def render_parameter_comparison_table(runs: Dict[str, Path],
                                      out_csv: Path,
                                      out_png: Path) -> None:
    """3-condition parameter comparison: rows = parameters, columns =
    Bare/Pre/Lam4."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    from acs.visualization.parameter_tables import PARAM_TABLE_SPEC, _resolve
    cfgs = {label: _load_cfg(rd) or {} for label, rd in runs.items()}
    rows = []
    rows.append({
        "poster_label": "Run name",
        "code_name": "run.name",
        **{label: cfgs[label].get("run", {}).get("name", "?") for label in runs},
        "unit": "—", "source": "config",
    })
    for poster_label, code_name, unit, source in PARAM_TABLE_SPEC:
        row = {"poster_label": poster_label, "code_name": code_name,
               "unit": unit, "source": source}
        for label in runs:
            row[label] = str(_resolve(cfgs[label], code_name))
        rows.append(row)
    headers = ["poster_label", "code_name"] + list(runs.keys()) + ["unit", "source"]
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        w.writerows(rows)
    # Render figure.
    fig_h = max(4.0, 0.30 * len(rows) + 1.0)
    fig, ax = plt.subplots(figsize=(13, fig_h))
    ax.axis("off")
    cell_data = [[r[h] for h in headers] for r in rows]
    table = ax.table(cellText=cell_data, colLabels=headers, cellLoc="left",
                     loc="upper left")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.2)
    for i in range(len(headers)):
        table[(0, i)].set_text_props(fontweight="bold", color="white")
        table[(0, i)].set_facecolor("#0d4ec9")
    fig.suptitle("Parameter comparison — Bare / Pre / Lam4", fontsize=12, y=0.99)
    fig.savefig(out_png, dpi=110, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bare", type=Path, required=False, default=None)
    ap.add_argument("--pre", type=Path, required=False, default=None)
    ap.add_argument("--lam4", type=Path, required=False, default=None)
    ap.add_argument("--out", type=Path, default=Path("results/production_comparison"))
    args = ap.parse_args()
    runs: Dict[str, Path] = {}
    if args.bare and args.bare.exists():
        runs["Bare"] = args.bare
    if args.pre and args.pre.exists():
        runs["Pre"] = args.pre
    if args.lam4 and args.lam4.exists():
        runs["Lam4"] = args.lam4
    if not runs:
        raise SystemExit("No valid run directories provided.")
    fig_dir = args.out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    print(f"[production_comparison] Conditions: {list(runs.keys())}")
    print("  → phenotype_triptych_topdown.png")
    render_triptych_topdown(runs, fig_dir / "phenotype_triptych.png")
    print("  → final_sideview_triptych.png")
    render_triptych_sideview(runs, fig_dir / "final_sideview_triptych.png")
    print("  → AA0_experiment_overlay.png")
    render_aa0_experiment_overlay(runs, fig_dir / "AA0_experiment_overlay.png")
    print("  → production_summary_dashboard.png")
    render_summary_dashboard(runs, fig_dir / "production_summary_dashboard.png")
    print("  → parameter_comparison_table.{csv,png}")
    render_parameter_comparison_table(
        runs,
        fig_dir / "parameter_comparison_table.csv",
        fig_dir / "parameter_comparison_table.png",
    )
    print("[production_comparison] Done.")


if __name__ == "__main__":
    main()
