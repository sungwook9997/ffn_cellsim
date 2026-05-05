"""Run the B-tier Phase E v2 pilot runner.

The runner uses a synthetic one-cell/two-FA fixture and writes HDF5 frames,
``diagnostics.csv``, and ``metadata.json``. It is plumbing only: no new
physics, no new measurement semantics, and no PI-data fitting path.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import subprocess
import sys


def _git_commit_hash(repo: str) -> str:
    try:
        out = subprocess.check_output(
            ["git", "-C", repo, "rev-parse", "HEAD"], text=True
        )
        return out.strip()
    except subprocess.CalledProcessError:
        return ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        help="Repository root used for imports, git hash, and default run output.",
    )
    parser.add_argument("--output-dir", default=None, help="Output directory.")
    parser.add_argument("--n-steps", type=int, default=4)
    parser.add_argument("--dt-s", type=float, default=0.0)
    parser.add_argument("--frame-interval", type=int, default=1)
    parser.add_argument("--grid", type=int, default=4)
    parser.add_argument("--spacing-um", type=float, default=1.0)
    parser.add_argument("--k-active", type=float, default=1.0)
    args = parser.parse_args()

    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    from acs.v2.phase_e_v2_pilot_runner import (  # noqa: E402
        PhaseEV2PilotConfig,
        run_phase_e_v2_pilot,
    )

    output_dir = args.output_dir
    if output_dir is None:
        timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        output_dir = os.path.join(repo, "runs", f"{timestamp}_phase_e_v2_pilot")
    output_dir = os.path.abspath(output_dir)

    config = PhaseEV2PilotConfig(
        n_steps=args.n_steps,
        dt_s=args.dt_s,
        frame_interval=args.frame_interval,
        grid_n=args.grid,
        spacing_um=args.spacing_um,
        k_active=args.k_active,
    )
    result = run_phase_e_v2_pilot(
        output_dir,
        config,
        git_commit_hash=_git_commit_hash(repo),
    )
    print("Phase E v2 pilot PASS")
    print(f"output_dir: {result.output_dir}")
    print(f"metadata: {result.metadata_path}")
    print(f"diagnostics: {result.csv_path}")
    print(f"frames: {len(result.frame_paths)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
