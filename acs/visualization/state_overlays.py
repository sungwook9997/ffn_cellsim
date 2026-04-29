"""acs.visualization.state_overlays — multi-channel particle-state overlays.

Per PI visualization directive 2026-04-29: optional overlay panels
showing per-particle state fields colorized over the xy projection.
Generates a 2×3 dashboard PNG per saved frame for a quick overview of
all state channels in one image.

Channels rendered (when present in HDF5 snapshots or metrics row):
- φ_memory map (formation phenotype)
- contact_activation map
- φ_eff map (legacy phi_p, post Stage 1b.b consumed by Layer 4)
- ρ_osm map
- ecm_strength (scalar; rendered as title only)
- γ_p (Stage 1d.b dynamic γ)

Headless / SSH compatible.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import h5py
except ImportError:
    h5py = None

import yaml


def _load_snapshot(run_dir: Path, frame_index: int) -> Optional[Dict[str, np.ndarray]]:
    """Load a single snapshot from the HDF5 file.

    Schema: f['frames']/<5-digit-id>/{position, velocity, ...}.
    Returns dict with renamed convenience aliases: x = position, v = velocity.
    """
    snaps_path = Path(run_dir) / "snapshots.h5"
    if h5py is None or not snaps_path.exists():
        return None
    with h5py.File(snaps_path, "r") as hf:
        grp_root = hf["frames"] if "frames" in hf else hf
        key = f"{frame_index:05d}"
        if key not in grp_root:
            for k in grp_root.keys():
                if k.startswith("frame_") and int(k.split("_")[1]) == frame_index:
                    key = k
                    break
                if int(k) == frame_index:
                    key = k
                    break
            else:
                return None
        grp = grp_root[key]
        out: Dict[str, np.ndarray] = {}
        for k in grp.keys():
            out[k] = np.array(grp[k])
        # Aliases for convenience
        if "position" in out and "x" not in out:
            out["x"] = out["position"]
        if "velocity" in out and "v" not in out:
            out["v"] = out["velocity"]
        return out


def render_state_overlay(
    run_dir: Path,
    frame_index: int,
    metrics_row: Optional[Dict[str, Any]] = None,
    out_path: Optional[Path] = None,
) -> Optional[Path]:
    """Render a 2×3 state-overlay dashboard for one frame.

    If the snapshot does not have the requested channel (e.g., older
    runs without phi_memory), that panel shows a placeholder.
    """
    snap = _load_snapshot(run_dir, frame_index)
    if snap is None:
        return None
    if out_path is None:
        out_path = (Path(run_dir) / "figures" / "overlays" /
                    f"overlay_{frame_index:05d}.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    x = snap.get("x")
    if x is None:
        return None

    # Stage 1d.c (PI directive 2026-04-30): expanded to 3×4 grid to
    # cover all per-particle channels including ECM/protrusion fields.
    fig, axes = plt.subplots(3, 4, figsize=(18, 12))
    # Channel inventory: (axis, key, title, cmap, vmin, vmax, special)
    channels = [
        # Row 1: Layer 3/4/5 state
        ((0, 0), "phi_memory", "φ_memory (formation)", "Blues", 0.0, 1.0, None),
        ((0, 1), "c_act", "contact_activation", "Oranges", 0.0, 1.0, None),
        ((0, 2), "phi", "φ_eff (Layer 4 input)", "viridis", 0.0, 1.0, None),
        ((0, 3), "rho_osm", "ρ_osm", "RdBu_r", 0.5, 1.6, None),
        # Row 2: Marangoni / kinematic
        ((1, 0), "gamma_p", "γ_p (dynamic)", "magma", None, None, None),
        ((1, 1), "v_norm", "‖v‖ speed", "plasma", 0.0, 0.05, None),
        ((1, 2), "pressure", "pressure (-tr/3)", "RdBu_r", None, None, None),
        ((1, 3), "dev_norm", "‖τ_dev‖ (deviatoric)", "viridis", None, None, None),
        # Row 3: Stage 1d.c ECM/protrusion
        ((2, 0), "protrusion_state", "protrusion state (0-4)", "Set1", -0.5, 4.5, "discrete"),
        ((2, 1), "fa_strength", "FA strength (0-1)", "Greens", 0.0, 1.0, None),
        ((2, 2), "ecm_signal", "|u_ecm| at particle", "magma", None, None, None),
        ((2, 3), "polarity", "polarity (xy arrows)", None, None, None, "vector"),
    ]
    cfg_path = Path(run_dir) / "config.yaml"
    cfg = None
    if cfg_path.exists():
        with open(cfg_path, encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
    domain = float(cfg["nondim"]["domain_star"]) if cfg else 6.0
    for entry in channels:
        (i, j), key, title, cmap, vmin, vmax = entry[:6]
        special = entry[6] if len(entry) > 6 else None
        ax = axes[i, j]
        # Per-channel field extraction.
        if key == "v_norm":
            v = snap.get("v")
            if v is None:
                ax.text(0.5, 0.5, "no v field", ha="center", va="center")
                ax.axis("off")
                continue
            field = np.linalg.norm(v, axis=1)
        elif special == "vector":
            # Polarity drawn as arrow quiver, colored by ‖p‖.
            pol = snap.get(key)
            if pol is None or pol.shape[0] != x.shape[0]:
                ax.text(0.5, 0.5, f"no {key} channel", ha="center", va="center")
                ax.axis("off")
                continue
            mag = np.linalg.norm(pol[:, :2], axis=1)
            # Subsample for clarity if too many particles.
            stride = max(1, x.shape[0] // 600)
            ax.quiver(
                x[::stride, 0], x[::stride, 1],
                pol[::stride, 0], pol[::stride, 1],
                mag[::stride],
                cmap="plasma", scale=20, width=0.003, alpha=0.85,
            )
            ax.set_xlim(0, domain)
            ax.set_ylim(0, domain)
            ax.set_aspect("equal")
            ax.set_title(title, fontsize=10)
            continue
        else:
            field = snap.get(key)
        if field is None or field.shape[0] != x.shape[0]:
            ax.text(0.5, 0.5, f"no {key} channel", ha="center", va="center")
            ax.axis("off")
            continue
        sc = ax.scatter(x[:, 0], x[:, 1], c=field, s=4, cmap=cmap,
                        vmin=vmin, vmax=vmax, alpha=0.85,
                        edgecolors="none")
        ax.set_xlim(0, domain)
        ax.set_ylim(0, domain)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    fig.suptitle(
        f"State overlays — {Path(run_dir).name} frame {frame_index}",
        fontsize=12, y=0.98,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=100, bbox_inches="tight")
    plt.close(fig)
    return out_path
