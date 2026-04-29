"""scripts/render_run_visuals.py — generate full visualization package
for a single run directory.

Per PI visualization directive 2026-04-29.

Usage:
    python -m scripts.render_run_visuals results/<run_name> \
        [--no-movies] [--max-frames 80] [--no-overlays]

Outputs (under results/<run_name>/figures/):
- live_topdown/frame_*.png
- live_sideview/frame_*.png
- movies/topdown_live.mp4 (or .gif fallback)
- movies/sideview_live.mp4 (or .gif fallback)
- overlays/overlay_*.png (1 per snapshot)
- dashboard/run_dashboard.png
- parameter_tables/run_parameters.csv
- parameter_tables/run_parameters.png
- final_frame_topdown.png
- final_frame_sideview.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

from acs.visualization import (
    install_helvetica_style,
    render_dashboard,
    render_final_frame_pair,
    render_parameter_table_figure,
    render_sideview_sequence,
    render_state_overlay,
    render_topdown_sequence,
    save_figure_dual,
    write_movie,
    write_mp4_and_gif,
    write_parameter_table_csv,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--no-movies", action="store_true",
                    help="Skip MP4/GIF rendering (PNG sequences still produced)")
    ap.add_argument("--no-overlays", action="store_true",
                    help="Skip per-frame state overlays (heaviest step)")
    ap.add_argument("--max-frames", type=int, default=None,
                    help="Cap number of frames rendered for testing")
    ap.add_argument("--fps", type=int, default=8)
    args = ap.parse_args()

    run_dir = args.run_dir.resolve()
    if not run_dir.exists():
        raise SystemExit(f"Run directory not found: {run_dir}")

    # CODEX_FIGURE_GUIDE.md: Helvetica typography, paper-grade rcParams.
    install_helvetica_style()

    print(f"[render_run_visuals] Run: {run_dir.name}")

    print("  → top-down PNG sequence")
    top_paths = render_topdown_sequence(run_dir, max_frames=args.max_frames)
    print(f"     ({len(top_paths)} frames)")

    print("  → side-view PNG sequence")
    side_paths = render_sideview_sequence(run_dir, max_frames=args.max_frames)
    print(f"     ({len(side_paths)} frames)")

    print("  → final-frame snapshots")
    render_final_frame_pair(run_dir)

    if not args.no_movies:
        print("  → top-down movie + GIF preview")
        out = write_mp4_and_gif(
            top_paths,
            run_dir / "figures" / "movies" / "topdown_live.mp4",
            fps=args.fps,
        )
        for k, p in out.items():
            print(f"    {k}: {p.name if p else '(none)'}")

        print("  → side-view movie + GIF preview")
        out = write_mp4_and_gif(
            side_paths,
            run_dir / "figures" / "movies" / "sideview_live.mp4",
            fps=args.fps,
        )
        for k, p in out.items():
            print(f"    {k}: {p.name if p else '(none)'}")

    if not args.no_overlays:
        print("  → state overlays")
        # Use the same frame indices as the top-down sequence.
        from acs.visualization.live_imaging import _load_run
        data = _load_run(run_dir)
        n = data["snaps"]["n_frames"]
        if args.max_frames is not None:
            n = min(n, args.max_frames)
        for i in range(n):
            frame_idx = data["snaps"]["frames"][i]
            mrow = data["metrics"][i] if i < len(data["metrics"]) else None
            render_state_overlay(run_dir, frame_idx, mrow)
        print(f"     ({n} overlay frames)")

    print("  → dashboard")
    render_dashboard(run_dir)

    print("  → parameter table")
    write_parameter_table_csv(run_dir)
    render_parameter_table_figure(run_dir)

    print("[render_run_visuals] Done.")


if __name__ == "__main__":
    main()
