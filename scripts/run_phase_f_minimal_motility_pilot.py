"""Run the B-tier Phase F minimal cell motility pilot runner.

The runner uses a synthetic single-cell + multi-FA fixture and writes
HDF5 frames + ``diagnostics.csv`` + ``metadata.json``. This is plumbing
only over locked ``step_phase_f_minimal_motility`` (commit ``bf03e83``);
no new physics, no new measurement semantics, no PI data use.

This pilot shows FA-traction-driven contour motion with Phase E v2
ECM-side diagnostics running in the loop. It does not yet implement
ECM->cell motility feedback, because HB#4 rate multipliers are
diagnostic-only here.
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
    parser.add_argument("--n-steps", type=int, default=200)
    parser.add_argument("--dt-cell-s", type=float, default=1e-3)
    parser.add_argument("--frame-interval", type=int, default=20)
    parser.add_argument("--grid", type=int, default=8)
    parser.add_argument("--spacing-um", type=float, default=0.5)
    parser.add_argument("--n-vertices", type=int, default=16)
    parser.add_argument("--radius-um", type=float, default=1.0)
    parser.add_argument("--lambda-c-nN", type=float, default=0.05)
    parser.add_argument("--sigma-c", type=float, default=0.05)
    parser.add_argument("--k-active", type=float, default=1.0)
    parser.add_argument("--fa-traction-nN", type=float, default=0.5)
    parser.add_argument("--fa-count", type=int, default=1)
    args = parser.parse_args()

    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    from acs.v2.phase_f_minimal_motility_pilot_runner import (  # noqa: E402
        PhaseFPilotConfig,
        run_phase_f_minimal_motility_pilot,
    )

    output_dir = args.output_dir
    if output_dir is None:
        timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        output_dir = os.path.join(
            repo, "runs", f"{timestamp}_phase_f_minimal_motility_pilot"
        )
    output_dir = os.path.abspath(output_dir)

    config = PhaseFPilotConfig(
        n_steps=args.n_steps,
        dt_cell_s=args.dt_cell_s,
        frame_interval=args.frame_interval,
        grid_n=args.grid,
        spacing_um=args.spacing_um,
        n_vertices=args.n_vertices,
        radius_um=args.radius_um,
        lambda_c_nN=args.lambda_c_nN,
        sigma_c_nN_per_um=args.sigma_c,
        k_active=args.k_active,
        fa_traction_nN=args.fa_traction_nN,
        fa_count=args.fa_count,
    )
    result = run_phase_f_minimal_motility_pilot(
        output_dir, config, git_commit_hash=_git_commit_hash(repo)
    )
    print("Phase F minimal motility pilot PASS")
    print(f"output_dir: {result.output_dir}")
    print(f"metadata: {result.metadata_path}")
    print(f"diagnostics: {result.csv_path}")
    print(f"frames: {len(result.frame_paths)}")
    if result.diagnostics:
        first = result.diagnostics[0]
        last = result.diagnostics[-1]
        print(
            f"centroid x: {first.centroid_x_um:.5f} -> {last.centroid_x_um:.5f} um"
        )
        print(
            f"attached vertex peak displacement: "
            f"{max(d.attached_vertex_displacement_um for d in result.diagnostics):.5f} um"
        )
        print(
            f"non-attached vertex peak displacement: "
            f"{max(d.max_non_attached_vertex_displacement_um for d in result.diagnostics):.5f} um"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
