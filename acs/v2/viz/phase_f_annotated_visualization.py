"""B-tier annotated visualization for Phase F minimal motility pilot.

Codex ``id=2103`` PI-readability fix: the bare frame replay / movie
showed ~0.016 um centroid displacement that is invisible at the
default 4 um polygon scale. This module produces annotated artifacts
that overlay the initial and current contour, mark the attached FA in
a distinct color, draw displacement vectors and centroid trace, and
print on-figure metrics so the PI can immediately read what moved
and why it is local pulling rather than rigid translation.

Per Codex ``id=2103`` honest-scope guardrail: when displacement is
amplified for readability, the figure title and console output state
``visual displacement scale = N×`` and the metadata records
``displacement_amplification`` so downstream consumers cannot mistake
the amplified visual for a physical displacement.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import re
import tempfile
from dataclasses import dataclass
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from acs.v2.output.frame_dump import read_frame


_FRAME_PATTERN = re.compile(r"^frame_(\d+)\.h5$")


@dataclass(frozen=True, slots=True)
class AnnotatedPhaseFArtifacts:
    """File paths emitted by Phase F annotated visualization helpers."""

    run_dir: str
    dashboard_png_path: Optional[str]
    dashboard_html_path: Optional[str]
    movie_mp4_path: Optional[str]
    movie_gif_path: Optional[str]
    n_frames_total: int
    displacement_amplification: float
    backend_available: bool


def _read_metadata(run_dir: str) -> dict[str, object]:
    metadata_path = os.path.join(run_dir, "metadata.json")
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(
            f"metadata.json missing at {metadata_path!r}; expected runner artifact"
        )
    with open(metadata_path, encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"{metadata_path!r} top-level is not a JSON object")
    runner = data.get("runner")
    if runner != "phase_f_minimal_motility_pilot":
        raise ValueError(
            f"{metadata_path!r} runner field {runner!r} is not "
            f"'phase_f_minimal_motility_pilot'; refusing to render"
        )
    return data


def _read_diagnostics(run_dir: str) -> list[dict[str, float]]:
    csv_path = os.path.join(run_dir, "diagnostics.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"diagnostics.csv missing at {csv_path!r}; expected runner artifact"
        )
    rows: list[dict[str, float]] = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            row: dict[str, float] = {}
            for k, v in raw.items():
                try:
                    row[k] = float(v)
                except (TypeError, ValueError):
                    row[k] = float("nan")
            rows.append(row)
    return rows


def _discover_frames(run_dir: str) -> list[str]:
    abs_run_dir = os.path.abspath(run_dir)
    candidates = sorted(glob.glob(os.path.join(abs_run_dir, "frame_*.h5")))
    return [
        path for path in candidates if _FRAME_PATTERN.match(os.path.basename(path))
    ]


def _read_frame_vertices_and_fa(path: str):
    """Return (vertices [N,2], fa_positions [(x,y),...], time_s)."""
    fr = read_frame(path)
    cells = list(fr.cluster.cells.values())
    if not cells:
        return None, [], fr.time_s
    cell = cells[0]
    import numpy as np
    verts = np.asarray(cell.measurement_boundary.vertices_xy_um, dtype=float)
    fa_positions = [
        (float(fa.position_um_xy[0]), float(fa.position_um_xy[1]))
        for fa in (cell.adhesions or [])
    ]
    return verts, fa_positions, float(fr.time_s)


def _attached_vertex_index(initial_vertices, fa_position) -> int:
    """Find the polygon vertex closest to the FA position (used as
    attached vertex marker; FA position follows attached vertex by
    Phase F contract, so initial closest vertex == attached)."""
    import numpy as np
    fp = np.asarray(fa_position, dtype=float)
    dists = np.linalg.norm(initial_vertices - fp, axis=1)
    return int(dists.argmin())


def _amplify_vertices(initial, current, scale: float):
    """Return amplified current vertices = initial + scale * (current - initial)."""
    return initial + scale * (current - initial)


def _draw_annotated_axes(
    ax,
    initial_verts,
    current_verts,
    fa_positions,
    attached_idx: int,
    centroid_trace,
    metrics: dict[str, float],
    *,
    displacement_amplification: float,
    title_prefix: str,
):
    """Draw one annotated panel (used by dashboard PNG and per-frame movie)."""
    import numpy as np

    # Closed polygons
    init_closed = np.vstack([initial_verts, initial_verts[:1]])
    curr_amp = _amplify_vertices(initial_verts, current_verts, displacement_amplification)
    curr_closed = np.vstack([curr_amp, curr_amp[:1]])

    ax.plot(init_closed[:, 0], init_closed[:, 1], color="#9ca0a6", linestyle="--",
            linewidth=1.0, label="initial contour")
    ax.plot(curr_closed[:, 0], curr_closed[:, 1], color="#1f6fbf", linestyle="-",
            linewidth=2.0, label="current contour (amp)" if displacement_amplification > 1.0 else "current contour")

    # Per-vertex displacement arrows (amplified)
    for i in range(initial_verts.shape[0]):
        dx = curr_amp[i, 0] - initial_verts[i, 0]
        dy = curr_amp[i, 1] - initial_verts[i, 1]
        if abs(dx) + abs(dy) < 1e-15:
            continue
        color = "#d62728" if i == attached_idx else "#bbbbbb"
        ax.arrow(
            initial_verts[i, 0],
            initial_verts[i, 1],
            dx,
            dy,
            head_width=0.02,
            head_length=0.02,
            fc=color,
            ec=color,
            alpha=0.85 if i == attached_idx else 0.45,
            length_includes_head=True,
        )

    # Attached vertex marker (initial + current)
    ax.plot(initial_verts[attached_idx, 0], initial_verts[attached_idx, 1],
            "o", color="#d62728", markersize=8, label="attached FA vertex (init)")
    ax.plot(curr_amp[attached_idx, 0], curr_amp[attached_idx, 1],
            "*", color="#d62728", markersize=12, label="attached FA vertex (now)")

    # FA positions (current)
    for fp in fa_positions:
        ax.plot(fp[0], fp[1], "x", color="#9467bd", markersize=8, alpha=0.8)

    # Centroid trace
    if len(centroid_trace) > 1:
        amp_trace = []
        c0 = np.asarray(centroid_trace[0], dtype=float)
        for c in centroid_trace:
            cv = np.asarray(c, dtype=float)
            amp_trace.append(c0 + displacement_amplification * (cv - c0))
        amp_trace = np.asarray(amp_trace)
        ax.plot(amp_trace[:, 0], amp_trace[:, 1], color="#2ca02c",
                linestyle="-", linewidth=1.0, alpha=0.8, label="centroid trace (amp)")
        ax.plot(amp_trace[-1, 0], amp_trace[-1, 1], "P", color="#2ca02c",
                markersize=9, label="centroid (now)")

    # Aspect + zoom
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, alpha=0.3)

    # On-figure metrics
    metrics_text = (
        f"t = {metrics.get('time_s', 0.0):.3f} s\n"
        f"centroid dx = {metrics.get('centroid_dx_um', 0.0):.5f} um\n"
        f"attached/non-attached peak ratio = {metrics.get('ratio', 0.0):.2f}\n"
        f"area delta = {metrics.get('area_delta_um2', 0.0):+.5f} um^2\n"
        f"perimeter delta = {metrics.get('perimeter_delta_um', 0.0):+.5f} um"
    )
    ax.text(
        0.02,
        0.98,
        metrics_text,
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#bbbbbb", alpha=0.9),
    )

    if displacement_amplification > 1.0:
        scale_label = f"visual displacement scale = {displacement_amplification:g}x (NOT physical)"
    else:
        scale_label = "visual displacement scale = 1x (physical)"
    ax.set_title(f"{title_prefix} · {scale_label}", fontsize=9)
    ax.legend(loc="lower right", fontsize=7, framealpha=0.85)


def _build_centroid_trace(rows: list[dict[str, float]]):
    return [
        (row.get("centroid_x_um", 0.0), row.get("centroid_y_um", 0.0))
        for row in rows
    ]


def _metrics_for_step(initial_row, current_row, ratio_fallback=None):
    """Compute on-figure metrics dict for a given step row."""
    centroid_dx = current_row.get("centroid_x_um", 0.0) - initial_row.get(
        "centroid_x_um", 0.0
    )
    attached = current_row.get("attached_vertex_displacement_um", 0.0)
    non_attached = current_row.get("max_non_attached_vertex_displacement_um", 0.0)
    if non_attached > 0.0:
        ratio = attached / non_attached
    else:
        ratio = float("inf") if attached > 0.0 else 0.0
    if ratio_fallback is not None and (ratio == float("inf") or ratio == 0.0):
        ratio = ratio_fallback
    area_delta = current_row.get("area_um2", 0.0) - initial_row.get("area_um2", 0.0)
    perimeter_delta = current_row.get("perimeter_um", 0.0) - initial_row.get(
        "perimeter_um", 0.0
    )
    return {
        "time_s": current_row.get("time_s", 0.0),
        "centroid_dx_um": centroid_dx,
        "ratio": ratio,
        "area_delta_um2": area_delta,
        "perimeter_delta_um": perimeter_delta,
    }


def render_annotated_phase_f_dashboard(
    run_dir: str,
    *,
    displacement_amplification: float = 1.0,
    write_html: bool = True,
) -> AnnotatedPhaseFArtifacts:
    """Render annotated Phase F dashboard PNG (and HTML rollup).

    Renders one annotated panel showing initial vs final contour with
    displacement arrows, attached FA marker, centroid trace, and
    on-figure metrics. Caller may set ``displacement_amplification > 1``
    to amplify the visual displacement; the title and metadata
    explicitly state the scale and that the amplification is NOT
    physical.
    """
    if displacement_amplification <= 0.0:
        raise ValueError(
            f"displacement_amplification must be positive, "
            f"got {displacement_amplification!r}"
        )
    abs_run_dir = os.path.abspath(run_dir)
    metadata = _read_metadata(abs_run_dir)
    rows = _read_diagnostics(abs_run_dir)
    frames = _discover_frames(abs_run_dir)
    if not frames:
        raise FileNotFoundError(
            f"no frame_*.h5 files in run_dir {abs_run_dir!r}"
        )

    init_verts, init_fa, t0 = _read_frame_vertices_and_fa(frames[0])
    final_verts, final_fa, tF = _read_frame_vertices_and_fa(frames[-1])
    attached_idx = _attached_vertex_index(init_verts, init_fa[0]) if init_fa else 0

    # Build centroid trace from diagnostics CSV
    centroid_trace = _build_centroid_trace(rows)
    # Initial row = pseudo at t=0 from initial vertices' centroid
    if rows:
        first_row = rows[0]
        last_row = rows[-1]
    else:
        import numpy as np
        first_row = {
            "centroid_x_um": float(init_verts[:, 0].mean()),
            "centroid_y_um": float(init_verts[:, 1].mean()),
            "time_s": 0.0,
            "perimeter_um": 0.0,
            "area_um2": 0.0,
            "attached_vertex_displacement_um": 0.0,
            "max_non_attached_vertex_displacement_um": 0.0,
        }
        last_row = first_row

    # Use peak attached / peak non-attached over the entire run for ratio
    if rows:
        peak_attached = max(r.get("attached_vertex_displacement_um", 0.0) for r in rows)
        peak_non_attached = max(
            r.get("max_non_attached_vertex_displacement_um", 0.0) for r in rows
        )
        if peak_non_attached > 0.0:
            run_ratio = peak_attached / peak_non_attached
        else:
            run_ratio = float("inf") if peak_attached > 0.0 else 0.0
    else:
        run_ratio = 0.0

    metrics = _metrics_for_step(first_row, last_row, ratio_fallback=run_ratio)
    metrics["ratio"] = run_ratio

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    _draw_annotated_axes(
        ax,
        init_verts,
        final_verts,
        final_fa,
        attached_idx,
        centroid_trace,
        metrics,
        displacement_amplification=displacement_amplification,
        title_prefix=f"Phase F annotated · {os.path.basename(abs_run_dir)}",
    )
    fig.tight_layout()
    png_path = os.path.join(abs_run_dir, "annotated_phase_f_dashboard.png")
    fig.savefig(png_path, dpi=140)
    plt.close(fig)

    html_path: Optional[str] = None
    if write_html:
        html_path = os.path.join(abs_run_dir, "annotated_phase_f_dashboard.html")
        amp_text = (
            f"visual displacement scale = {displacement_amplification:g}x (NOT physical)"
            if displacement_amplification > 1.0
            else "visual displacement scale = 1x (physical)"
        )
        run_label = os.path.basename(abs_run_dir)
        html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Phase F annotated dashboard - {run_label}</title>
<style>body {{font-family:-apple-system,sans-serif;margin:24px;color:#222}}
h1 {{font-size:18px}} h2 {{font-size:14px;margin-top:18px}}
table {{border-collapse:collapse;font-size:12px}}
th,td {{border:1px solid #bbb;padding:4px 8px;text-align:left}}
img {{max-width:100%;border:1px solid #ddd}}
.note {{color:#777;font-size:11px}}</style></head><body>
<h1>Phase F annotated dashboard - {run_label}</h1>
<p class="note">{amp_text}. Display-only annotated visualization over the locked
Phase F minimal motility pilot output. NOT a validation metric.</p>
<h2>Annotated panel</h2>
<img src="annotated_phase_f_dashboard.png" alt="annotated phase f panel">
<h2>Run metrics summary</h2>
<table>
<tr><th>centroid dx (final)</th><td>{metrics['centroid_dx_um']:+.5f} um</td></tr>
<tr><th>attached/non-attached peak ratio</th><td>{metrics['ratio']:.2f}</td></tr>
<tr><th>area delta (final)</th><td>{metrics['area_delta_um2']:+.5f} um^2</td></tr>
<tr><th>perimeter delta (final)</th><td>{metrics['perimeter_delta_um']:+.5f} um</td></tr>
<tr><th>final time</th><td>{metrics['time_s']:.3f} s</td></tr>
<tr><th>n frames</th><td>{len(frames)}</td></tr>
<tr><th>displacement_amplification</th><td>{displacement_amplification:g}x</td></tr>
</table>
</body></html>
"""
        with open(html_path, "w", encoding="utf-8") as fh:
            fh.write(html)

    return AnnotatedPhaseFArtifacts(
        run_dir=abs_run_dir,
        dashboard_png_path=png_path,
        dashboard_html_path=html_path,
        movie_mp4_path=None,
        movie_gif_path=None,
        n_frames_total=len(frames),
        displacement_amplification=float(displacement_amplification),
        backend_available=True,
    )


def render_annotated_phase_f_movie(
    run_dir: str,
    *,
    displacement_amplification: float = 1.0,
    fps: int = 5,
    gif_fps: int = 5,
    write_gif: bool = True,
    keep_frame_pngs: bool = False,
) -> AnnotatedPhaseFArtifacts:
    """Render annotated Phase F MP4 (+ optional GIF). One annotated
    panel per frame; explicit visual amplification label."""

    if displacement_amplification <= 0.0:
        raise ValueError(
            f"displacement_amplification must be positive, "
            f"got {displacement_amplification!r}"
        )
    if fps <= 0 or gif_fps <= 0:
        raise ValueError("fps and gif_fps must be positive")

    abs_run_dir = os.path.abspath(run_dir)
    _read_metadata(abs_run_dir)
    rows = _read_diagnostics(abs_run_dir)
    frames = _discover_frames(abs_run_dir)
    if not frames:
        raise FileNotFoundError(
            f"no frame_*.h5 files in run_dir {abs_run_dir!r}"
        )

    init_verts, init_fa, _ = _read_frame_vertices_and_fa(frames[0])
    attached_idx = _attached_vertex_index(init_verts, init_fa[0]) if init_fa else 0
    centroid_trace = _build_centroid_trace(rows)

    # Compute peak ratio for the whole run (consistent across frames)
    if rows:
        peak_attached = max(r.get("attached_vertex_displacement_um", 0.0) for r in rows)
        peak_non_attached = max(
            r.get("max_non_attached_vertex_displacement_um", 0.0) for r in rows
        )
        if peak_non_attached > 0.0:
            run_ratio = peak_attached / peak_non_attached
        else:
            run_ratio = float("inf") if peak_attached > 0.0 else 0.0
    else:
        run_ratio = 0.0

    if rows:
        first_row = rows[0]
    else:
        first_row = {
            "centroid_x_um": float(init_verts[:, 0].mean()),
            "centroid_y_um": float(init_verts[:, 1].mean()),
            "time_s": 0.0,
            "perimeter_um": 0.0,
            "area_um2": 0.0,
        }

    backend_available = False
    write_mp4_and_gif = None
    try:
        from acs.visualization import live_imaging as _live_imaging
        write_mp4_and_gif = _live_imaging.write_mp4_and_gif
        backend_available = _live_imaging.imageio is not None
    except Exception:  # pragma: no cover
        pass

    if keep_frame_pngs:
        png_dir = os.path.join(abs_run_dir, "annotated_movie_frames")
        os.makedirs(png_dir, exist_ok=True)
        cleanup_dir: Optional[str] = None
    else:
        png_dir = tempfile.mkdtemp(prefix="annotated_phase_f_")
        cleanup_dir = png_dir

    mp4_path: Optional[str] = None
    gif_path: Optional[str] = None

    try:
        from pathlib import Path as _Path
        png_paths: list[_Path] = []
        for i, frame_path in enumerate(frames):
            verts, fa_positions, t_s = _read_frame_vertices_and_fa(frame_path)
            # Per-frame metrics
            if i < len(rows):
                row_now = rows[i]
            else:
                row_now = first_row
            metrics = _metrics_for_step(
                first_row, row_now, ratio_fallback=run_ratio
            )
            metrics["ratio"] = run_ratio
            fig, ax = plt.subplots(figsize=(8.0, 6.0))
            _draw_annotated_axes(
                ax,
                init_verts,
                verts,
                fa_positions,
                attached_idx,
                centroid_trace[: i + 1] if rows else [],
                metrics,
                displacement_amplification=displacement_amplification,
                title_prefix=f"Phase F annotated · frame {i}",
            )
            fig.tight_layout()
            base = os.path.splitext(os.path.basename(frame_path))[0]
            out_png = os.path.join(png_dir, f"{base}_annotated.png")
            fig.savefig(out_png, dpi=120)
            plt.close(fig)
            png_paths.append(_Path(out_png))

        if write_mp4_and_gif is not None and backend_available:
            mp4_target = _Path(os.path.join(abs_run_dir, "annotated_phase_f_movie.mp4"))
            result = write_mp4_and_gif(png_paths, mp4_target, fps=fps, gif_fps=gif_fps)
            mp4_returned = result.get("mp4")
            gif_returned = result.get("gif")
            if mp4_returned is not None and _Path(mp4_returned).exists():
                if _Path(mp4_returned).suffix.lower() == ".mp4":
                    mp4_path = str(_Path(mp4_returned).resolve())
                elif _Path(mp4_returned).suffix.lower() == ".gif" and write_gif:
                    gif_path = str(_Path(mp4_returned).resolve())
            if write_gif and gif_returned is not None and _Path(gif_returned).exists():
                if _Path(gif_returned).suffix.lower() == ".gif":
                    gif_path = str(_Path(gif_returned).resolve())
            elif not write_gif and gif_returned is not None:
                try:
                    os.remove(str(gif_returned))
                except OSError:
                    pass

        return AnnotatedPhaseFArtifacts(
            run_dir=abs_run_dir,
            dashboard_png_path=None,
            dashboard_html_path=None,
            movie_mp4_path=mp4_path,
            movie_gif_path=gif_path,
            n_frames_total=len(frames),
            displacement_amplification=float(displacement_amplification),
            backend_available=backend_available,
        )
    finally:
        if cleanup_dir is not None and os.path.isdir(cleanup_dir):
            for fname in os.listdir(cleanup_dir):
                try:
                    os.remove(os.path.join(cleanup_dir, fname))
                except OSError:
                    pass
            try:
                os.rmdir(cleanup_dir)
            except OSError:
                pass
