"""Render annotated Phase F dashboard + movie from a runner output dir.

B-tier visualization-packaging per Codex `id=2103` (PI urgent
readability fix). Overlays initial vs current contour, marks attached
FA vertex, draws displacement arrows + centroid trace, and prints
on-figure metrics. Optional ``--displacement-amplification N`` makes
small displacements visible; the title and metadata explicitly state
``visual displacement scale = N×ㅡ NOT physical``.
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
    )
    parser.add_argument("--run-dir", required=True)
    parser.add_argument(
        "--displacement-amplification",
        type=float,
        default=1.0,
        help="Visual amplification factor for displacement vectors / current contour. "
        "Title and metadata explicitly state if > 1 (NOT physical).",
    )
    parser.add_argument("--fps", type=int, default=5)
    parser.add_argument("--gif-fps", type=int, default=5)
    parser.add_argument("--no-gif", action="store_true")
    parser.add_argument("--no-html", action="store_true")
    parser.add_argument("--no-movie", action="store_true")
    parser.add_argument("--keep-frame-pngs", action="store_true")
    args = parser.parse_args()

    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    from acs.v2.viz.phase_f_annotated_visualization import (  # noqa: E402
        render_annotated_phase_f_dashboard,
        render_annotated_phase_f_movie,
    )

    dashboard = render_annotated_phase_f_dashboard(
        args.run_dir,
        displacement_amplification=args.displacement_amplification,
        write_html=not args.no_html,
    )
    print("Phase F annotated dashboard PASS")
    print(f"png: {dashboard.dashboard_png_path}")
    if dashboard.dashboard_html_path:
        print(f"html: {dashboard.dashboard_html_path}")
    print(f"n_frames_total: {dashboard.n_frames_total}")
    print(f"displacement_amplification: {dashboard.displacement_amplification}")

    if not args.no_movie:
        movie = render_annotated_phase_f_movie(
            args.run_dir,
            displacement_amplification=args.displacement_amplification,
            fps=args.fps,
            gif_fps=args.gif_fps,
            write_gif=not args.no_gif,
            keep_frame_pngs=args.keep_frame_pngs,
        )
        print("Phase F annotated movie PASS")
        print(f"backend_available: {movie.backend_available}")
        print(f"mp4: {movie.movie_mp4_path}")
        print(f"gif: {movie.movie_gif_path}")
        if not movie.backend_available:
            print(
                "NOTE: imageio backend unavailable; mp4/gif not produced. "
                "Install imageio + ffmpeg to enable."
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
