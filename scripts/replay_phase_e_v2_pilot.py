"""Render the Phase E v2 pilot frame replay grid from a runner output dir.

B-tier viewer over the ``frame_*.h5`` artifacts emitted by
``scripts/run_phase_e_v2_pilot.py``. Produces ``frame_replay.png`` and
optionally ``frame_replay.html`` inside the same run directory.

Plumbing only: no new physics, no new measurement semantics, no PI data use.
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
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Skip the HTML rollup; emit only the replay grid PNG.",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=12,
        help="Cap on replay grid panel count (default: 12). Evenly-spaced "
        "subsampling preserves first + last frames.",
    )
    args = parser.parse_args()

    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    from acs.v2.viz.phase_e_v2_pilot_frame_replay import (  # noqa: E402
        render_phase_e_v2_pilot_frame_replay,
    )

    artifacts = render_phase_e_v2_pilot_frame_replay(
        args.run_dir,
        write_html=not args.no_html,
        max_frames=args.max_frames,
    )
    print("Phase E v2 pilot frame replay PASS")
    print(f"run_dir: {artifacts.run_dir}")
    print(f"png: {artifacts.png_path}")
    if artifacts.html_path:
        print(f"html: {artifacts.html_path}")
    print(f"n_frames_total: {artifacts.n_frames_total}")
    print(f"n_frames_rendered: {artifacts.n_frames_rendered}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
