"""Cell / cortex matplotlib visualisation (no KU; supports Unit 3.1 figs)."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from acs_kb.cell.cortex import Cortex


def _orientation_colors(angles: np.ndarray) -> np.ndarray:
    hue = (np.mod(angles, np.pi) / np.pi).astype(float)
    hsv = np.stack(
        [hue, np.full_like(hue, 0.85), np.full_like(hue, 0.85)], axis=1
    )
    return mcolors.hsv_to_rgb(hsv)


def plot_cell_static(
    cortex: Cortex,
    *,
    ax: plt.Axes | None = None,
    save_path: str | Path | None = None,
    show: bool = False,
    title: str | None = None,
    show_xls: bool = True,
    show_center: bool = True,
    pad_factor: float = 1.4,
) -> plt.Figure:
    """Plot a single cortex: fibers (orientation-colored) + XLs + centre."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure

    colors = _orientation_colors(cortex.fiber_orientations)
    for f in range(cortex.n_fibers):
        pts = cortex.bead_positions[f]
        ax.plot(pts[:, 0], pts[:, 1], color=colors[f], lw=0.8)

    if show_xls and cortex.cross_links:
        flat = cortex.bead_positions.reshape(-1, 2)
        N = cortex.n_beads_per_fiber
        xs, ys = [], []
        for xl in cortex.cross_links:
            ia = xl.fiber_a * N + xl.bead_a
            ib = xl.fiber_b * N + xl.bead_b
            mid = 0.5 * (flat[ia] + flat[ib])
            xs.append(mid[0])
            ys.append(mid[1])
        ax.scatter(
            xs, ys, s=6, c="red", alpha=0.7, zorder=3,
            label=f"XL ({len(cortex.cross_links)})",
        )

    if show_center:
        pts = cortex.bead_positions.reshape(-1, 2)
        centroid = pts.mean(axis=0)
        ax.scatter(
            [centroid[0]], [centroid[1]], marker="+", c="black", s=80,
            zorder=4, label="centroid",
        )

    pts = cortex.bead_positions.reshape(-1, 2)
    centroid = pts.mean(axis=0)
    span = pad_factor * np.linalg.norm(pts - centroid, axis=1).max()
    ax.set_xlim(centroid[0] - span, centroid[0] + span)
    ax.set_ylim(centroid[1] - span, centroid[1] + span)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    if title:
        ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()

    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
    if show:
        plt.show()
    return fig


def plot_aspect_ratio_timeseries(
    times: np.ndarray,
    aspect_ratios: np.ndarray,
    *,
    acceptance: float | None = 1.2,
    acceptance_window_s: float | None = 60.0,
    save_path: str | Path | None = None,
    label: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """Plot aspect-ratio vs time with the KU-3.1 acceptance band overlay."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    else:
        fig = ax.figure
    ax.plot(times, aspect_ratios, lw=1.5, label=label or "aspect ratio")
    if acceptance is not None:
        ax.axhline(acceptance, color="red", ls="--", lw=1.0,
                   label=f"KU-3.1 cap = {acceptance:.2f}")
    if acceptance_window_s is not None:
        ax.axvline(acceptance_window_s, color="grey", ls=":", lw=1.0,
                   label=f"window = {acceptance_window_s:.0f} s")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("aspect ratio  √(λ_max / λ_min)")
    ax.set_title("Cortex rounding (KU-3.1 VALIDATION)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
    return fig


def plot_rounding_before_after(
    cortex_initial: Cortex,
    cortex_final: Cortex,
    *,
    save_path: str | Path | None = None,
) -> plt.Figure:
    """Side-by-side static plot: initial vs final cortex."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    plot_cell_static(cortex_initial, ax=axes[0], title="t = 0 s (initial)")
    plot_cell_static(cortex_final, ax=axes[1], title="t = final (relaxed)")
    fig.tight_layout()
    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
    return fig


def animate_cell_rounding(
    frames: Sequence[np.ndarray],
    times: Sequence[float],
    *,
    fiber_orientations: np.ndarray,
    save_path: str | Path,
    fps: int = 10,
    pad_factor: float = 1.4,
) -> Path:
    """Save an animation of cortex bead positions across frames.

    Uses matplotlib's ``FuncAnimation``. Output format inferred from
    the suffix of ``save_path`` (``.gif`` uses Pillow writer, ``.mp4``
    uses ffmpeg if available; fall back to GIF on ImportError).
    """
    from matplotlib.animation import FuncAnimation, PillowWriter

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)

    colors = _orientation_colors(fiber_orientations)
    fig, ax = plt.subplots(figsize=(6, 6))

    pts0 = np.asarray(frames[0])
    n_fibers, n_beads, _ = pts0.shape
    lines = [
        ax.plot(pts0[f, :, 0], pts0[f, :, 1], color=colors[f], lw=0.8)[0]
        for f in range(n_fibers)
    ]
    centroid0 = pts0.reshape(-1, 2).mean(axis=0)
    spans = []
    for fr in frames:
        flat = np.asarray(fr).reshape(-1, 2)
        cen = flat.mean(axis=0)
        spans.append(np.linalg.norm(flat - cen, axis=1).max())
    span = pad_factor * max(spans)
    ax.set_xlim(centroid0[0] - span, centroid0[0] + span)
    ax.set_ylim(centroid0[1] - span, centroid0[1] + span)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    time_text = ax.text(
        0.02, 0.97, "", transform=ax.transAxes, va="top",
        fontsize=10, family="monospace",
    )

    def update(frame_idx: int):
        pts = np.asarray(frames[frame_idx])
        for f in range(n_fibers):
            lines[f].set_data(pts[f, :, 0], pts[f, :, 1])
        time_text.set_text(f"t = {times[frame_idx]:.1f} s")
        return lines + [time_text]

    anim = FuncAnimation(
        fig, update, frames=len(frames), blit=False, interval=1000 / fps
    )
    writer = PillowWriter(fps=fps)
    anim.save(save_path, writer=writer)
    plt.close(fig)
    return save_path
