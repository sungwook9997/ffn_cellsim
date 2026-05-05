"""B-tier MP4/GIF movie generator for the Phase E v2 pilot runner.

Companion to ``acs.v2.viz.phase_e_v2_pilot_dashboard`` (commit ``107c4fb``)
and ``acs.v2.viz.phase_e_v2_pilot_frame_replay`` (commit
``40b2edf``+``21495a2``). Where the dashboard renders ``diagnostics.csv``
and the frame replay renders a static multi-panel grid, this module
renders the per-step HDF5 frame artifacts into an animated MP4 (with
optional GIF preview).

Per Codex ``id=1783`` Cycle 15 mp4 routing scope:

1. Reuse v1 :func:`acs.visualization.live_imaging.write_mp4_and_gif`
   (PNG sequence -> MP4 + downsampled GIF preview via ``imageio``).
2. Render single-panel per-frame PNGs into a temporary directory using
   :func:`acs.v2.viz.stub3d.render_frame_png`, then composite via the
   v1 helper.
3. Graceful fallback when ``imageio`` or the MP4 codec is unavailable:
   :func:`render_phase_e_v2_pilot_movie` returns artifacts whose
   ``mp4_path`` and/or ``gif_path`` are ``None`` instead of raising.
4. Display-only artifact: the movie is visualization, NOT a validation
   metric. No new measurement extraction; the per-frame PNGs are
   discarded after compositing unless ``keep_frame_pngs=True``.

Per Codex ``id=1728`` (sister-cycle) gardrails inherited:
- Runner artifact consumer (frame_*.h5 + metadata.json) only
- No new physics, no new measurement semantics, no PI data use
- Schema-correct error surfacing: missing frames -> FileNotFoundError,
  wrong runner -> ValueError
"""

from __future__ import annotations

import glob
import json
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from acs.v2.viz.stub3d import render_frame_png


_FRAME_PATTERN = re.compile(r"^frame_(\d+)\.h5$")
_DEFAULT_MP4_NAME = "pilot_movie.mp4"
_DEFAULT_GIF_NAME = "pilot_movie.gif"
_DEFAULT_FPS = 2
_DEFAULT_GIF_FPS = 2
_DEFAULT_GIF_MAX_FRAMES = 80


@dataclass(frozen=True, slots=True)
class PilotMovieArtifacts:
    """File paths emitted by :func:`render_phase_e_v2_pilot_movie`."""

    run_dir: str
    mp4_path: Optional[str]
    gif_path: Optional[str]
    n_frames_total: int
    n_frames_rendered: int
    frame_png_paths: tuple[str, ...]
    fps: int
    gif_fps: int
    backend_available: bool


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
    accepted_runners = ("phase_e_v2_pilot", "phase_f_minimal_motility_pilot")
    if runner not in accepted_runners:
        raise ValueError(
            f"{metadata_path!r} runner field {runner!r} is not in "
            f"{accepted_runners}; refusing to render unrelated artifacts"
        )
    return data


def render_phase_e_v2_pilot_movie(
    run_dir: str,
    *,
    fps: int = _DEFAULT_FPS,
    gif_fps: int = _DEFAULT_GIF_FPS,
    write_gif: bool = True,
    keep_frame_pngs: bool = False,
) -> PilotMovieArtifacts:
    """Render an MP4 (and optional GIF preview) from the pilot runner's frames.

    Workflow:
      1. Validate ``metadata.json`` (must have ``runner == "phase_e_v2_pilot"``).
      2. Discover ``frame_*.h5`` files in ``run_dir`` (sorted, regex-filtered).
      3. Render each frame to a temporary single-panel PNG via
         :func:`acs.v2.viz.stub3d.render_frame_png` (Agg backend, no display).
      4. Composite the PNG sequence into MP4 (and optional GIF) via the
         v1 helper :func:`acs.visualization.live_imaging.write_mp4_and_gif`.
      5. Save MP4 / GIF to ``<run_dir>/pilot_movie.mp4`` /
         ``<run_dir>/pilot_movie.gif``.
      6. Clean up the temporary PNGs unless ``keep_frame_pngs=True``.

    Args:
        run_dir: Pilot runner output directory (must contain
            ``frame_*.h5`` and ``metadata.json``).
        fps: MP4 frames-per-second (default 2 — slow enough to read each
            step at the small ``n_steps`` typical of the pilot runner).
        gif_fps: GIF preview frames-per-second (default 2).
        write_gif: If True (default), also emit a GIF preview alongside
            the MP4. The v1 helper downsamples the GIF if frame count
            exceeds 80 (``gif_max_frames``).
        keep_frame_pngs: If True, keep the per-frame PNGs in a
            ``movie_frames/`` subdirectory of ``run_dir``. Default False
            (PNGs are written to a temp dir and discarded after
            compositing).

    Returns:
        :class:`PilotMovieArtifacts` with absolute paths to the MP4 /
        GIF artifacts (each ``None`` if the corresponding artifact could
        not be produced; see ``backend_available`` for whether the
        ``imageio`` backend was importable at all).

    Raises:
        FileNotFoundError: ``run_dir`` doesn't exist, or ``metadata.json``
            is missing, or no ``frame_*.h5`` files are present.
        ValueError: ``metadata.runner != "phase_e_v2_pilot"``, or
            ``fps`` / ``gif_fps`` is non-positive.
    """

    if not os.path.isdir(run_dir):
        raise FileNotFoundError(
            f"run_dir {run_dir!r} is not an existing directory"
        )
    if fps <= 0:
        raise ValueError(f"fps must be positive, got {fps!r}")
    if gif_fps <= 0:
        raise ValueError(f"gif_fps must be positive, got {gif_fps!r}")

    abs_run_dir = os.path.abspath(run_dir)
    metadata_path = os.path.join(abs_run_dir, "metadata.json")
    _read_metadata_json(metadata_path)  # Validate runner

    frame_paths = _discover_frames(abs_run_dir)
    if not frame_paths:
        raise FileNotFoundError(
            f"no frame_*.h5 files in run_dir {abs_run_dir!r}; "
            f"runner must have emitted at least one frame"
        )

    # Lazy import the v1 helper module. Per Codex `id=1814` BLOCKER 1:
    # the v1 module imports `imageio` defensively (sets it to None on
    # ImportError) so a successful `from ... import write_mp4_and_gif`
    # does NOT imply the imageio backend is actually available.
    # Confirm `live_imaging.imageio is not None` before treating the
    # backend as available.
    try:
        from acs.visualization import live_imaging as _live_imaging
        write_mp4_and_gif = _live_imaging.write_mp4_and_gif
        backend_available = _live_imaging.imageio is not None
    except Exception:  # pragma: no cover - defensive
        write_mp4_and_gif = None  # type: ignore[assignment]
        backend_available = False

    if keep_frame_pngs:
        png_dir = os.path.join(abs_run_dir, "movie_frames")
        os.makedirs(png_dir, exist_ok=True)
        cleanup_dir: Optional[str] = None
    else:
        png_dir = tempfile.mkdtemp(prefix="phase_e_v2_pilot_movie_")
        cleanup_dir = png_dir

    try:
        png_paths: list[Path] = []
        for h5_path in frame_paths:
            base = os.path.splitext(os.path.basename(h5_path))[0]
            png_path = os.path.join(png_dir, f"{base}.png")
            render_frame_png(h5_path, png_path)
            png_paths.append(Path(png_path))

        mp4_target = os.path.join(abs_run_dir, _DEFAULT_MP4_NAME)
        mp4_path: Optional[str] = None
        gif_path: Optional[str] = None

        if write_mp4_and_gif is not None and backend_available:
            result = write_mp4_and_gif(
                png_paths,
                Path(mp4_target),
                fps=fps,
                gif_fps=gif_fps,
                gif_max_frames=_DEFAULT_GIF_MAX_FRAMES,
            )
            mp4_returned = result.get("mp4")
            gif_returned = result.get("gif")

            # Per Codex `id=1814` BLOCKER 2: v1 `write_movie(..., gif=False)`
            # falls back to writing a `.gif` if the MP4 codec is unavailable
            # and returns that GIF path under `result['mp4']`. Classify by
            # suffix + existence so a `.gif` returned in the mp4 slot is
            # promoted to gif_path (when write_gif=True) or removed (when
            # write_gif=False), preventing mp4_path from pointing at a .gif.
            mp4_resolved: Optional[Path] = (
                Path(mp4_returned) if mp4_returned is not None else None
            )
            gif_resolved: Optional[Path] = (
                Path(gif_returned) if gif_returned is not None else None
            )

            if mp4_resolved is not None and mp4_resolved.exists():
                if mp4_resolved.suffix.lower() == ".mp4":
                    mp4_path = str(mp4_resolved.resolve())
                elif mp4_resolved.suffix.lower() == ".gif":
                    # Codec-fallback GIF surfaced in the mp4 slot.
                    if write_gif and gif_resolved is None:
                        # Promote to gif_path when no separate GIF returned.
                        gif_resolved = mp4_resolved
                    elif not write_gif:
                        try:
                            os.remove(str(mp4_resolved))
                        except OSError:
                            pass

            if gif_resolved is not None and gif_resolved.exists():
                if gif_resolved.suffix.lower() == ".gif":
                    if write_gif:
                        gif_path = str(gif_resolved.resolve())
                    else:
                        # write_gif=False: caller asked to suppress GIF.
                        try:
                            os.remove(str(gif_resolved))
                        except OSError:
                            pass

        return PilotMovieArtifacts(
            run_dir=abs_run_dir,
            mp4_path=mp4_path,
            gif_path=gif_path,
            n_frames_total=len(frame_paths),
            n_frames_rendered=len(png_paths),
            frame_png_paths=tuple(str(p) for p in png_paths) if keep_frame_pngs else (),
            fps=fps,
            gif_fps=gif_fps,
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
