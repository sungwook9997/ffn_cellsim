"""Render the Phase E v2 pilot dashboard from a runner output directory.

B-tier viewer over the artifacts emitted by ``scripts/run_phase_e_v2_pilot.py``
(``diagnostics.csv``, ``metadata.json``). Produces ``dashboard.png`` and
optionally ``dashboard.html`` inside the same run directory. Plumbing only:
no new physics, no new measurement semantics, no PI data use.
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
        help="Pilot runner output directory (must contain diagnostics.csv "
        "and metadata.json).",
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Skip the HTML rollup; emit only the PNG dashboard.",
    )
    args = parser.parse_args()

    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    from acs.v2.viz.phase_e_v2_pilot_dashboard import (  # noqa: E402
        render_phase_e_v2_pilot_dashboard,
    )

    artifacts = render_phase_e_v2_pilot_dashboard(
        args.run_dir,
        write_html=not args.no_html,
    )
    print("Phase E v2 pilot dashboard PASS")
    print(f"run_dir: {artifacts.run_dir}")
    print(f"png: {artifacts.png_path}")
    if artifacts.html_path:
        print(f"html: {artifacts.html_path}")
    print(f"n_diagnostic_rows: {artifacts.n_diagnostic_rows}")
    print(f"n_frames: {artifacts.n_frames}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
