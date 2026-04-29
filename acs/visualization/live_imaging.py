"""acs.visualization.live_imaging — microscope-style live imaging frames.

Reads HDF5 snapshots + metrics.csv from a results/{run.name}/ directory
and renders microscope-style top-down (xy projection) and side-view
(xz projection) PNG frames + MP4 movies + GIF previews per the PI
visualization directive 2026-04-29.

PI rules followed:
- Top-down = xy projection of ALL particles (CLAUDE.md Hard Rule 11
  matching PI experimental Area_um2 imaging modality).
- Per-frame text overlay: time stamp in hr, scale bar in μm, condition
  label, A/A₀_topdown value, R/R₀.
- Side-view shows substrate plane and contact band.
- Headless SSH compatible (Agg backend).
- Outputs: figures/live_topdown/*.png, figures/live_sideview/*.png,
  figures/movies/topdown_live.mp4, figures/movies/sideview_live.mp4.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

try:
    import h5py
except ImportError:
    h5py = None

try:
    import imageio.v2 as imageio
except ImportError:
    try:
        import imageio
    except ImportError:
        imageio = None


def _load_run(run_dir: Path) -> Dict[str, Any]:
    """Load snapshots, metrics, config from a run directory."""
    run_dir = Path(run_dir)
    # Config
    config_paths = [run_dir / "config.yaml", run_dir / ".." / "configs" /
                    f"{run_dir.name}.yaml"]
    cfg = None
    for cp in config_paths:
        if cp.exists():
            with open(cp, encoding="utf-8") as f:
                cfg = yaml.safe_load(f)
            break
    # Snapshots (.h5) — actual schema: f['frames']/<5-digit-id>/{position,
    # velocity, deformation_gradient, is_boundary, stress_tensor}.
    snaps_path = run_dir / "snapshots.h5"
    snaps: Dict[str, Any] = {"frames": [], "x": [], "phi": [], "n_frames": 0}
    if h5py is not None and snaps_path.exists():
        with h5py.File(snaps_path, "r") as hf:
            grp_root = hf["frames"] if "frames" in hf else hf
            keys = sorted(grp_root.keys(), key=lambda s: int(s))
            for k in keys:
                grp = grp_root[k]
                snaps["frames"].append(int(k))
                # Position field naming.
                if "position" in grp:
                    snaps["x"].append(np.array(grp["position"]))
                elif "x" in grp:
                    snaps["x"].append(np.array(grp["x"]))
                else:
                    continue
                # Optional φ field.
                if "phi" in grp:
                    snaps["phi"].append(np.array(grp["phi"]))
            snaps["n_frames"] = len(snaps["x"])
    # metrics.csv
    metrics_path = run_dir / "metrics.csv"
    metrics: List[Dict[str, Any]] = []
    if metrics_path.exists():
        import csv
        with open(metrics_path, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                metrics.append(row)
    return {"snaps": snaps, "metrics": metrics, "cfg": cfg, "run_dir": run_dir}


def _condition_label(cfg: Optional[Dict[str, Any]]) -> str:
    """Infer Bare/Pre/Lam4 from config phi_initial."""
    if not cfg:
        return "?"
    phi_init = cfg.get("layer3", {}).get("phi_initial", None)
    if phi_init is None:
        return "?"
    if phi_init <= 0.40:
        return "Bare"
    elif phi_init <= 0.65:
        return "Pre"
    else:
        return "Lam4"


def _t_hr_from_metrics(row: Dict[str, Any], cfg: Optional[Dict[str, Any]]) -> float:
    """Convert frame time_star to physical hours via maxwell_tau_s."""
    try:
        ts = float(row.get("time_star", 0.0))
    except (TypeError, ValueError):
        return 0.0
    tau = 60.0
    if cfg:
        try:
            tau = float(cfg["physics"]["maxwell_tau_s"])
        except (KeyError, TypeError, ValueError):
            pass
    return ts * tau / 3600.0  # seconds → hours


def _scale_bar_um_from_cfg(cfg: Optional[Dict[str, Any]]) -> float:
    """Return scale bar length in star units corresponding to ~50 μm."""
    if not cfg:
        return 0.5  # in star units
    try:
        R0_um = float(cfg["initialization"]["radius_um"])
    except (KeyError, TypeError, ValueError):
        R0_um = 100.0
    # Scale bar = 50 μm in physical, convert to star (R₀ = 1 star)
    return 50.0 / R0_um


def render_topdown_frame(
    x: np.ndarray,
    metrics_row: Dict[str, Any],
    cfg: Optional[Dict[str, Any]],
    out_path: Path,
    fig_size: Tuple[float, float] = (6, 6),
) -> None:
    """Render a single top-down (xy projection) frame.

    Microscope-style logic per PI directive:
    - All particles plotted as small filled circles.
    - Convex-hull outline of xy projection drawn over.
    - Time stamp (hr), scale bar (μm), condition label (Bare/Pre/Lam4),
      A/A₀_topdown, R/R₀ overlaid as text.
    """
    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_aspect("equal")
    ax.scatter(x[:, 0], x[:, 1], s=4, c="#3a82f7", alpha=0.55,
               edgecolors="none")
    # Convex hull outline
    try:
        from scipy.spatial import ConvexHull
        if x.shape[0] >= 3:
            hull = ConvexHull(x[:, :2])
            poly = x[:, :2][hull.vertices]
            poly = np.vstack([poly, poly[:1]])
            ax.plot(poly[:, 0], poly[:, 1], color="#0d4ec9", linewidth=1.5)
    except Exception:
        pass
    # Frame box: domain
    if cfg:
        try:
            dom = float(cfg["nondim"]["domain_star"])
            ax.set_xlim(0, dom)
            ax.set_ylim(0, dom)
        except (KeyError, TypeError, ValueError):
            pass
    # Text overlays
    label = _condition_label(cfg)
    t_hr = _t_hr_from_metrics(metrics_row, cfg)
    aa0 = metrics_row.get("A_over_A0_topdown", "?")
    try:
        aa0_str = f"{float(aa0):.3f}"
    except (TypeError, ValueError):
        aa0_str = str(aa0)
    R = metrics_row.get("effective_radius", "?")
    try:
        R_str = f"{float(R):.3f}"
    except (TypeError, ValueError):
        R_str = str(R)
    ax.text(0.02, 0.98, f"{label}", transform=ax.transAxes,
            ha="left", va="top", fontsize=14, fontweight="bold",
            color="#0d4ec9")
    ax.text(0.02, 0.92, f"t = {t_hr:.2f} hr", transform=ax.transAxes,
            ha="left", va="top", fontsize=11, color="#222")
    ax.text(0.02, 0.87, f"A/A₀ = {aa0_str}", transform=ax.transAxes,
            ha="left", va="top", fontsize=11, color="#222")
    ax.text(0.02, 0.82, f"R/R₀ = {R_str}", transform=ax.transAxes,
            ha="left", va="top", fontsize=11, color="#222")
    # Scale bar
    sb_star = _scale_bar_um_from_cfg(cfg)
    if cfg:
        try:
            dom = float(cfg["nondim"]["domain_star"])
            x0 = 0.07 * dom
            y0 = 0.07 * dom
            ax.plot([x0, x0 + sb_star], [y0, y0], "k-", linewidth=3)
            ax.text(x0 + sb_star / 2, y0 + 0.01 * dom, "50 μm",
                    ha="center", va="bottom", fontsize=10, color="#000")
        except (KeyError, TypeError, ValueError):
            pass
    ax.set_xlabel("x*")
    ax.set_ylabel("y*")
    ax.set_title("Top-down projection (live-imaging style)", fontsize=11)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def render_sideview_frame(
    x: np.ndarray,
    metrics_row: Dict[str, Any],
    cfg: Optional[Dict[str, Any]],
    out_path: Path,
    fig_size: Tuple[float, float] = (8, 5),
) -> None:
    """Render a single side-view (xz) frame.

    Substrate plane drawn at z=0; contact band shaded; spheroid COM-z
    annotated. Time stamp + condition + R/R₀ overlaid.
    """
    fig, ax = plt.subplots(figsize=fig_size)
    ax.scatter(x[:, 0], x[:, 2], s=4, c="#3a82f7", alpha=0.55,
               edgecolors="none")
    # Substrate plane
    if cfg:
        try:
            dom = float(cfg["nondim"]["domain_star"])
            ax.axhline(0.0, color="#444", linewidth=1.5, linestyle="--")
            # Contact band
            n_band = int(cfg.get("substrate", {}).get("n_contact_band", 3))
            dx_star = float(cfg["nondim"]["domain_star"]) / int(
                cfg["numerics"].get("background_grid_resolution", 64))
            h_band = n_band * dx_star
            ax.axhspan(0.0, h_band, color="#fde2c8", alpha=0.5,
                       label="Contact band")
            ax.set_xlim(0, dom)
            ax.set_ylim(-0.05 * dom, 0.5 * dom)
        except (KeyError, TypeError, ValueError):
            pass
    # COM-z line
    com_z = x[:, 2].mean()
    ax.axhline(com_z, color="#9b1d20", linewidth=1.0, linestyle=":",
               label=f"COM_z* = {com_z:.3f}")
    label = _condition_label(cfg)
    t_hr = _t_hr_from_metrics(metrics_row, cfg)
    R = metrics_row.get("effective_radius", "?")
    try:
        R_str = f"{float(R):.3f}"
    except (TypeError, ValueError):
        R_str = str(R)
    ax.text(0.02, 0.95, f"{label}  t = {t_hr:.2f} hr  R/R₀ = {R_str}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=12, fontweight="bold", color="#0d4ec9")
    ax.set_xlabel("x*")
    ax.set_ylabel("z* (substrate at z=0)")
    ax.set_title("Side-view (xz projection)", fontsize=11)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def render_topdown_sequence(run_dir: Path, max_frames: Optional[int] = None) -> List[Path]:
    """Render every snapshot frame as a top-down PNG.

    Returns list of PNG paths (sorted by frame index).
    """
    data = _load_run(run_dir)
    snaps = data["snaps"]
    metrics = data["metrics"]
    cfg = data["cfg"]
    out_dir = Path(run_dir) / "figures" / "live_topdown"
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    n = snaps["n_frames"]
    if max_frames is not None:
        n = min(n, max_frames)
    for i in range(n):
        frame_idx = snaps["frames"][i]
        x = snaps["x"][i]
        if i < len(metrics):
            mrow = metrics[i]
        else:
            mrow = {}
        out_path = out_dir / f"frame_{frame_idx:05d}.png"
        render_topdown_frame(x, mrow, cfg, out_path)
        paths.append(out_path)
    return paths


def render_sideview_sequence(run_dir: Path, max_frames: Optional[int] = None) -> List[Path]:
    """Render every snapshot frame as a side-view (xz) PNG."""
    data = _load_run(run_dir)
    snaps = data["snaps"]
    metrics = data["metrics"]
    cfg = data["cfg"]
    out_dir = Path(run_dir) / "figures" / "live_sideview"
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    n = snaps["n_frames"]
    if max_frames is not None:
        n = min(n, max_frames)
    for i in range(n):
        frame_idx = snaps["frames"][i]
        x = snaps["x"][i]
        if i < len(metrics):
            mrow = metrics[i]
        else:
            mrow = {}
        out_path = out_dir / f"frame_{frame_idx:05d}.png"
        render_sideview_frame(x, mrow, cfg, out_path)
        paths.append(out_path)
    return paths


def write_movie(png_paths: List[Path], out_path: Path,
                fps: int = 8, gif: bool = False) -> Optional[Path]:
    """Write an MP4 (or GIF) movie from a sequence of PNG frames.

    Returns out_path on success, None if imageio is missing or there
    are no frames.
    """
    if imageio is None or not png_paths:
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames = []
    for p in png_paths:
        try:
            frames.append(imageio.imread(str(p)))
        except Exception:
            continue
    if not frames:
        return None
    try:
        if gif:
            imageio.mimsave(str(out_path), frames, fps=fps, loop=0)
        else:
            imageio.mimsave(str(out_path), frames, fps=fps,
                            codec="libx264", quality=8)
    except Exception as e:
        # Common fallback: GIF if MP4 codec missing.
        if not gif:
            try:
                gif_path = out_path.with_suffix(".gif")
                imageio.mimsave(str(gif_path), frames, fps=fps, loop=0)
                return gif_path
            except Exception:
                pass
        return None
    return out_path


def render_final_frame_pair(run_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Render only the FINAL frame in both modalities.

    Useful for quick poster/manuscript snapshots without rendering the
    full sequence.
    """
    data = _load_run(run_dir)
    snaps = data["snaps"]
    metrics = data["metrics"]
    cfg = data["cfg"]
    if snaps["n_frames"] == 0:
        return None, None
    i = snaps["n_frames"] - 1
    x = snaps["x"][i]
    mrow = metrics[i] if i < len(metrics) else {}
    fig_dir = Path(run_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    top_path = fig_dir / "final_frame_topdown.png"
    side_path = fig_dir / "final_frame_sideview.png"
    render_topdown_frame(x, mrow, cfg, top_path)
    render_sideview_frame(x, mrow, cfg, side_path)
    return top_path, side_path
