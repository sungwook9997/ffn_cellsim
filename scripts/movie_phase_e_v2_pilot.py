"""Render the Phase E v2 pilot MP4 (+ optional GIF) from a runner output dir.

B-tier movie generator over the ``frame_*.h5`` artifacts emitted by
``scripts/run_phase_e_v2_pilot.py``. Reuses v1 ffmpeg/imageio backend
via :func:`acs.visualization.live_imaging.write_mp4_and_gif`.

Plumbing only: no new physics, no new measurement semantics, no PI data use.
If the ``imageio`` backend is unavailable, the script reports the missing
backend and exits cleanly without raising (display-only artifact;
graceful-skip semantics).
"""

from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        help="Repository root used for imports.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Pilot runner output directory (must contain frame_*.h5 + "
        "metadata.json).",
    )
    parser.add_argument("--fps", type=int, default=2)
    parser.add_argument("--gif-fps", type=int, default=2)
    parser.add_argument(
        "--no-gif",
        action="store_true",
        help="Skip the GIF preview; emit only the MP4.",
    )
    parser.add_argument(
        "--keep-frame-pngs",
        action="store_true",
        help="Keep the per-frame PNGs in <run_dir>/movie_frames/.",
    )
    args = parser.parse_args()

    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    from acs.v2.viz.phase_e_v2_pilot_movie import (  # noqa: E402
        render_phase_e_v2_pilot_movie,
    )

    artifacts = render_phase_e_v2_pilot_movie(
        args.run_dir,
        fps=args.fps,
        gif_fps=args.gif_fps,
        write_gif=not args.no_gif,
        keep_frame_pngs=args.keep_frame_pngs,
    )
    print("Phase E v2 pilot movie PASS")
    print(f"run_dir: {artifacts.run_dir}")
    print(f"backend_available: {artifacts.backend_available}")
    print(f"mp4: {artifacts.mp4_path}")
    print(f"gif: {artifacts.gif_path}")
    print(f"n_frames_total: {artifacts.n_frames_total}")
    print(f"n_frames_rendered: {artifacts.n_frames_rendered}")
    if not artifacts.backend_available:
        print(
            "NOTE: imageio backend unavailable; mp4/gif not produced. "
            "Install imageio + ffmpeg to enable this output. "
            "Display-only artifact; graceful skip per B-tier scope."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
