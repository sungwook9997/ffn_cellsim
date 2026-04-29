"""acs.visualization.dashboard — production run dashboard figure.

Per PI visualization directive 2026-04-29: each production run gets a
single dashboard PNG containing:
- A/A₀_topdown(t)
- R/R₀(t)
- Wadell sphericity ψ(t)
- φ_memory + contact_activation summaries (boundary mean, contact mean,
  global mean as diagnostic only)
- ρ_osm mean / min / max
- ecm_strength + mmp_total (Layer 6)
- contact area diagnostic
- horizontal momentum drift absolute
- Gate PASS/FAIL summary text
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml


def _safe_float(v: Any, default: float = float("nan")) -> float:
    try:
        f = float(v)
        if np.isnan(f) or np.isinf(f):
            return default
        return f
    except (TypeError, ValueError):
        return default


def _load_metrics_csv(path: Path) -> List[Dict[str, str]]:
    import csv
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _load_gate_report(path: Path) -> Dict[str, Any]:
    """Parse gate_report.md to extract Overall + per-check PASS/FAIL."""
    if not path.exists():
        return {"overall": "?", "checks": []}
    overall = "?"
    checks: List[Dict[str, str]] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("- **Overall**"):
                overall = line.split(":", 1)[1].strip().strip("*").strip()
            elif line.startswith("- **[PASS]**") or line.startswith("- **[FAIL]**"):
                status = "PASS" if "[PASS]" in line else "FAIL"
                # Strip the bullet prefix and bold marker.
                rest = line.split("**", 2)[2].lstrip(":").strip()
                # Take the colon-prefixed gate name.
                if ":" in rest:
                    name = rest.split(":")[0].strip()
                else:
                    name = rest
                checks.append({"status": status, "name": name})
    return {"overall": overall, "checks": checks}


def render_dashboard(run_dir: Path, out_path: Optional[Path] = None) -> Path:
    """Render the run dashboard PNG.

    Auto-saves to `run_dir/figures/dashboard/run_dashboard.png` if
    `out_path` is None.
    """
    run_dir = Path(run_dir)
    metrics = _load_metrics_csv(run_dir / "metrics.csv")
    gate = _load_gate_report(run_dir / "gate_report.md")
    cfg_path = run_dir / "config.yaml"
    cfg = None
    if cfg_path.exists():
        with open(cfg_path, encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
    if out_path is None:
        out_path = run_dir / "figures" / "dashboard" / "run_dashboard.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not metrics:
        # Empty-data placeholder.
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, f"No metrics for {run_dir.name}", ha="center",
                va="center", fontsize=14)
        ax.set_axis_off()
        fig.savefig(out_path, dpi=110)
        plt.close(fig)
        return out_path

    # Convert relevant columns to arrays (numeric-safe).
    def col(name: str) -> np.ndarray:
        return np.array([_safe_float(r.get(name, float("nan"))) for r in metrics])

    t = col("time_star")
    aa0 = col("A_over_A0_topdown")
    R = col("effective_radius")
    psi = col("wadell_sphericity")
    n_part = float(cfg.get("simulation", {}).get("n_material_points", 1))
    if n_part <= 0:
        n_part = 1.0
    phi_mem = col("phi_memory_sum") / n_part
    c_act = col("c_act_sum") / n_part
    c_act_band_mean = np.array([
        _safe_float(r.get("c_act_band_sum", 0.0)) / max(_safe_float(r.get("c_act_band_count", 0.0), 1.0), 1.0)
        for r in metrics
    ])
    rho_osm = col("rho_osm_sum") / n_part
    rho_min = col("rho_osm_min")
    rho_max = col("rho_osm_max")
    ecm = col("ecm_strength")
    mmp = col("mmp_total")
    contact_area = col("contact_area_xy_hull")
    px = col("momentum_x")
    py = col("momentum_y")
    px0 = px[0] if px.size else 0.0
    py0 = py[0] if py.size else 0.0
    dpxy_abs = np.sqrt((px - px0) ** 2 + (py - py0) ** 2)

    # 4×3 grid layout.
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 4, hspace=0.4, wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(t, aa0, color="#0d4ec9")
    ax.set_xlabel("t* (τ_relax)")
    ax.set_ylabel("A/A₀_topdown")
    ax.set_title("A/A₀ (top-down projection)")
    ax.axhline(1.0, color="#aaa", linestyle="--", linewidth=0.8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(t, R, color="#9b1d20")
    ax.set_xlabel("t*")
    ax.set_ylabel("R/R₀ (effective radius)")
    ax.set_title("Effective radius")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 2])
    ax.plot(t, psi, color="#3a8c4f")
    ax.axhline(0.7, color="#999", linestyle="--", linewidth=0.8,
               label="post-spread limit")
    ax.set_xlabel("t*")
    ax.set_ylabel("ψ Wadell sphericity")
    ax.set_title("Wadell sphericity")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 3])
    ax.plot(t, dpxy_abs, color="#bb6f00")
    ax.axhline(2.0e-3, color="#999", linestyle="--", linewidth=0.8,
               label="abs limit 2e-3")
    ax.set_xlabel("t*")
    ax.set_ylabel("|Δp_xy|")
    ax.set_title("Horizontal momentum drift (abs)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(t, phi_mem, color="#0d4ec9", label="<φ_memory>")
    ax.plot(t, c_act, color="#bb6f00", label="<c_act> global", linestyle="--")
    ax.plot(t, c_act_band_mean, color="#9b1d20", label="<c_act>_band")
    ax.set_xlabel("t*")
    ax.set_ylabel("φ / c_act")
    ax.set_title("Layer 3 split (Stage 1b.b)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 1])
    ax.plot(t, rho_osm, color="#0d4ec9", label="<ρ_osm>")
    ax.plot(t, rho_min, color="#888", linestyle=":", label="min")
    ax.plot(t, rho_max, color="#888", linestyle="--", label="max")
    ax.set_xlabel("t*")
    ax.set_ylabel("ρ_osm")
    ax.set_title("Layer 5 osmotic state")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 2])
    ax.plot(t, ecm, color="#3a8c4f", label="ecm_strength")
    ax2 = ax.twinx()
    ax2.plot(t, mmp, color="#9b1d20", label="mmp_total")
    ax.set_xlabel("t*")
    ax.set_ylabel("ecm_strength", color="#3a8c4f")
    ax2.set_ylabel("mmp_total", color="#9b1d20")
    ax.set_title("Layer 6 chemistry / ECM")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 3])
    ax.plot(t, contact_area, color="#0d4ec9")
    ax.set_xlabel("t*")
    ax.set_ylabel("contact_area_xy_hull")
    ax.set_title("Contact patch area (diagnostic)")
    ax.grid(alpha=0.3)

    # Gate summary panel (spans two columns).
    ax = fig.add_subplot(gs[2, :2])
    ax.axis("off")
    overall = gate.get("overall", "?")
    color = {"PASS": "#2a8a4d", "FAIL": "#a02828"}.get(overall, "#666")
    ax.text(0.0, 1.0, f"Overall: {overall}", fontsize=14,
            fontweight="bold", color=color, va="top")
    y = 0.92
    n_pass = sum(1 for c in gate.get("checks", []) if c["status"] == "PASS")
    n_fail = sum(1 for c in gate.get("checks", []) if c["status"] == "FAIL")
    ax.text(0.0, y, f"PASS: {n_pass} / FAIL: {n_fail}", fontsize=11,
            va="top")
    y -= 0.07
    for c in gate.get("checks", []):
        col_name = "#2a8a4d" if c["status"] == "PASS" else "#a02828"
        ax.text(0.0, y, f"[{c['status']}] {c['name']}",
                fontsize=8, color=col_name, va="top")
        y -= 0.045
        if y < 0.0:
            break

    # Run identity panel.
    ax = fig.add_subplot(gs[2, 2:])
    ax.axis("off")
    label = "?"
    if cfg:
        phi_init = cfg.get("layer3", {}).get("phi_initial", None)
        if phi_init is not None:
            if phi_init <= 0.40:
                label = "Bare"
            elif phi_init <= 0.65:
                label = "Pre"
            else:
                label = "Lam4"
    n_p = int(cfg.get("simulation", {}).get("n_material_points", 0)) if cfg else 0
    total_t = float(cfg.get("simulation", {}).get("total_sim_time_s", 0)) if cfg else 0
    info_lines = [
        f"Run: {run_dir.name}",
        f"Condition: {label}   n_particles: {n_p}",
        f"Total sim time: {total_t/3600:.1f} hr",
    ]
    if cfg:
        seed = cfg.get("run", {}).get("seed", "?")
        info_lines.append(f"Seed: {seed}")
        layer3_split = cfg.get("layer3", {}).get("split", False)
        if layer3_split:
            info_lines.append("Layer 3: split (φ_memory + c_act)")
        layer4_dyn = cfg.get("layer4", {}).get("dynamic_gamma", False)
        if layer4_dyn:
            aA = cfg["layer4"].get("alpha_A_star", 0.0)
            aF = cfg["layer4"].get("alpha_F_star", 0.0)
            info_lines.append(f"Layer 4: dynamic γ (α_A={aA}, α_F={aF})")
        layer2_b = cfg.get("layer2_b", {})
        if layer2_b.get("enabled", False):
            info_lines.append(
                f"Layer 2.b: stochastic events (λ={layer2_b.get('lambda_lam_star', 0)},"
                f" impulse={layer2_b.get('impulse_lam_star', 0)})"
            )
    y = 1.0
    for ln in info_lines:
        ax.text(0.0, y, ln, fontsize=10, va="top", family="monospace")
        y -= 0.08

    fig.suptitle(f"Run dashboard — {run_dir.name}", fontsize=14, y=0.995)
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return out_path
