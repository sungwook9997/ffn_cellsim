"""B-tier visual evidence frame replay viewer for the Phase E v2 pilot runner.

Companion to ``acs.v2.viz.phase_e_v2_pilot_dashboard``. Where the dashboard
renders ``diagnostics.csv`` trajectories, this module renders the per-step
HDF5 frame artifacts produced by ``acs.v2.phase_e_v2_pilot_runner`` as a
multi-panel composite PNG plus an optional HTML rollup.

Per Codex ``id=1728`` 4 가드레일 (B-tier scope discipline):

1. 기존 ``acs.v2.viz.stub3d`` + frame read path 재사용. Replay panels
   reuse :func:`acs.v2.viz.stub3d._render_axes` so the per-frame
   visualization stays consistent with single-frame ``render_frame_png``.
2. Dependency-light: matplotlib only (Agg backend, headless). No GIF /
   MP4 / external video tooling.
3. Tests: frame-count / schema / path / non-empty output. NOT pixel-perfect.
4. Errors: missing ``frame_*.h5`` files → :class:`FileNotFoundError`;
   wrong runner / malformed metadata → :class:`ValueError`. No silent
   degradation past obvious schema mismatches.

Display-only artifact: replay panels are visualization, NOT validation
metrics. No new measurement extraction from frames; titles/HTML labels
explicitly say "display-only" to preserve the wording boundary.
"""

from __future__ import annotations

import glob
import json
import os
import re
from dataclasses import dataclass
from typing import Optional

import matplotlib

matplotlib.use("Agg")  # headless / SSH-friendly
import matplotlib.pyplot as plt  # noqa: E402

from acs.v2.output.frame_dump import read_frame
from acs.v2.viz.stub3d import _render_axes


_REPLAY_PNG_NAME = "frame_replay.png"
_REPLAY_HTML_NAME = "frame_replay.html"
_FRAME_PATTERN = re.compile(r"^frame_(\d+)\.h5$")
_DEFAULT_MAX_FRAMES = 12


@dataclass(frozen=True, slots=True)
class PilotFrameReplayArtifacts:
    """File paths emitted by :func:`render_phase_e_v2_pilot_frame_replay`."""

    run_dir: str
    png_path: str
    html_path: Optional[str]
    frame_paths_rendered: tuple[str, ...]
    n_frames_total: int
    n_frames_rendered: int


def _discover_frames(run_dir: str) -> list[str]:
    abs_run_dir = os.path.abspath(run_dir)
    candidates = sorted(glob.glob(os.path.join(abs_run_dir, "frame_*.h5")))
    valid = [
        path for path in candidates if _FRAME_PATTERN.match(os.path.basename(path))
    ]
    return valid


def _read_metadata_json(metadata_path: str) -> dict[str, object]:
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(
            f"metadata.json missing at {metadata_path!r}; expected runner artifact"
        )
    with open(metadata_path, encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"{metadata_path!r} top-level is not a JSON object")
    runner = data.get("runner")
    if runner != "phase_e_v2_pilot":
        raise ValueError(
            f"{metadata_path!r} runner field {runner!r} is not "
            f"'phase_e_v2_pilot'; refusing to render unrelated artifacts"
        )
    return data


def _select_frames(frame_paths: list[str], *, max_frames: int) -> list[str]:
    if max_frames <= 0:
        raise ValueError(f"max_frames must be positive, got {max_frames!r}")
    if len(frame_paths) <= max_frames:
        return list(frame_paths)
    # max_frames == 1: return only the first frame (Codex id=1756 fix —
    # the evenly-spaced formula divides by (max_frames - 1) which would
    # zero-divide; the natural viewer-cap semantics for "at most one panel"
    # is the initial frame).
    if max_frames == 1:
        return [frame_paths[0]]
    # Evenly-spaced subsample including first + last to preserve trajectory ends.
    n = len(frame_paths)
    indices = [round(i * (n - 1) / (max_frames - 1)) for i in range(max_frames)]
    # Deduplicate while preserving order (rounding can collide for small spans).
    seen: set[int] = set()
    selected: list[str] = []
    for idx in indices:
        if idx in seen:
            continue
        seen.add(idx)
        selected.append(frame_paths[idx])
    return selected


def _grid_shape(n_panels: int) -> tuple[int, int]:
    if n_panels <= 1:
        return (1, 1)
    if n_panels <= 2:
        return (1, 2)
    if n_panels <= 4:
        return (2, 2)
    if n_panels <= 6:
        return (2, 3)
    if n_panels <= 9:
        return (3, 3)
    if n_panels <= 12:
        return (3, 4)
    return (4, 4)


def _render_replay_png(
    frame_paths: list[str],
    png_path: str,
    *,
    run_label: str,
) -> tuple[str, ...]:
    if not frame_paths:
        # Edge case: zero-step run produces only the initial frame; if that
        # somehow goes missing, render an explicit placeholder PNG.
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        ax.text(
            0.5,
            0.5,
            "no frame_*.h5 files in run directory",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        fig.suptitle(
            f"Phase E v2 pilot replay · {run_label} (display-only viewer)",
            fontsize=10,
        )
        fig.tight_layout()
        fig.savefig(png_path, dpi=120)
        plt.close(fig)
        return ()

    rows, cols = _grid_shape(len(frame_paths))
    fig, axes = plt.subplots(rows, cols, figsize=(3.5 * cols, 3.5 * rows))
    if rows * cols == 1:
        flat_axes = [axes]
    else:
        flat_axes = list(axes.flat)

    rendered: list[str] = []
    for ax, path in zip(flat_axes, frame_paths):
        result = read_frame(path)
        _render_axes(ax, result.cluster, result.time_s, height_tint_scale_um=10.0)
        ax.set_title(
            f"{os.path.basename(path)} · t={result.time_s:.3f}s",
            fontsize=8,
        )
        rendered.append(path)

    # Hide any unused axes in the grid.
    for ax in flat_axes[len(frame_paths):]:
        ax.set_axis_off()

    fig.suptitle(
        f"Phase E v2 pilot replay · {run_label} "
        f"(display-only viewer; NOT a validation metric)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=120)
    plt.close(fig)
    return tuple(rendered)


def _render_replay_html(
    frame_paths: list[str],
    metadata: dict[str, object],
    html_path: str,
    *,
    png_basename: str,
    run_label: str,
    n_total: int,
) -> None:
    notes = metadata.get("notes", [])
    git_commit_hash = metadata.get("git_commit_hash", "unrecorded")

    rows: list[str] = []
    for path in frame_paths:
        try:
            result = read_frame(path)
            t_s = f"{result.time_s:.6f}"
            n_cells = len(result.cluster.cells)
        except Exception as exc:  # pragma: no cover - defensive
            t_s = f"read-error: {exc!s}"
            n_cells = -1
        rows.append(
            f"<tr><td>{os.path.basename(path)}</td><td>{t_s}</td>"
            f"<td>{n_cells}</td></tr>"
        )
    table_rows = "".join(rows)
    notes_items = "".join(
        f"<li>{n}</li>" for n in (notes if isinstance(notes, list) else [])
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Phase E v2 pilot frame replay - {run_label}</title>
<style>
body {{ font-family: -apple-system, sans-serif; margin: 24px; color: #222; }}
h1 {{ font-size: 18px; }}
h2 {{ font-size: 14px; margin-top: 18px; }}
table {{ border-collapse: collapse; font-size: 12px; }}
th, td {{ border: 1px solid #bbb; padding: 4px 8px; text-align: left; }}
.note {{ color: #777; font-size: 11px; }}
img {{ max-width: 100%; border: 1px solid #ddd; }}
</style>
</head>
<body>
<h1>Phase E v2 pilot frame replay - {run_label}</h1>
<p class="note">Display-only viewer over runner frame artifacts.
NOT a validation metric. No new measurement extraction; this page
visualizes existing frame contents only.</p>

<h2>Run metadata</h2>
<table>
<tr><th>git_commit_hash</th><td>{git_commit_hash}</td></tr>
<tr><th>n_frames_total_in_run</th><td>{n_total}</td></tr>
<tr><th>n_frames_rendered</th><td>{len(frame_paths)}</td></tr>
</table>

<h2>Replay panel (grid)</h2>
<img src="{png_basename}" alt="Phase E v2 pilot frame replay grid">

<h2>Frame index (display-only)</h2>
<table><tr><th>frame</th><th>time_s</th><th>n_cells</th></tr>
{table_rows}</table>

<h2>Scope notes (inherited from runner metadata)</h2>
<ul>{notes_items}</ul>
</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)


def render_phase_e_v2_pilot_frame_replay(
    run_dir: str,
    *,
    write_html: bool = True,
    max_frames: int = _DEFAULT_MAX_FRAMES,
) -> PilotFrameReplayArtifacts:
    """Render an HDF5 frame replay grid PNG (and optional HTML rollup).

    Args:
        run_dir: Directory previously produced by
            :func:`acs.v2.phase_e_v2_pilot_runner.run_phase_e_v2_pilot`.
        write_html: If True, also emit a self-contained HTML rollup that
            embeds the replay grid PNG and a per-frame index table.
        max_frames: Cap on grid panel count. If the run has more frames,
            evenly-spaced subsampling preserves first + last + intermediate.

    Returns:
        :class:`PilotFrameReplayArtifacts` with absolute paths of the
        generated artifacts and the list of frame paths actually rendered.

    Raises:
        FileNotFoundError: if ``run_dir`` does not exist or contains no
            ``frame_*.h5`` files (subject to ``metadata.json`` validation
            having passed first).
        ValueError: if ``metadata.json`` is missing the
            ``runner == "phase_e_v2_pilot"`` field, or if ``max_frames``
            is non-positive.
    """

    if not os.path.isdir(run_dir):
        raise FileNotFoundError(
            f"run_dir {run_dir!r} is not an existing directory"
        )
    if max_frames <= 0:
        raise ValueError(f"max_frames must be positive, got {max_frames!r}")

    abs_run_dir = os.path.abspath(run_dir)
    metadata_path = os.path.join(abs_run_dir, "metadata.json")
    metadata = _read_metadata_json(metadata_path)

    frame_paths = _discover_frames(abs_run_dir)
    if not frame_paths:
        raise FileNotFoundError(
            f"no frame_*.h5 files in run_dir {abs_run_dir!r}; "
            f"runner must have emitted at least one frame"
        )

    selected = _select_frames(frame_paths, max_frames=max_frames)
    run_label = os.path.basename(abs_run_dir.rstrip(os.sep)) or abs_run_dir
    png_path = os.path.join(abs_run_dir, _REPLAY_PNG_NAME)
    rendered = _render_replay_png(selected, png_path, run_label=run_label)

    html_path: Optional[str] = None
    if write_html:
        html_path = os.path.join(abs_run_dir, _REPLAY_HTML_NAME)
        _render_replay_html(
            list(rendered),
            metadata,
            html_path,
            png_basename=_REPLAY_PNG_NAME,
            run_label=run_label,
            n_total=len(frame_paths),
        )

    return PilotFrameReplayArtifacts(
        run_dir=abs_run_dir,
        png_path=png_path,
        html_path=html_path,
        frame_paths_rendered=tuple(rendered),
        n_frames_total=len(frame_paths),
        n_frames_rendered=len(rendered),
    )
