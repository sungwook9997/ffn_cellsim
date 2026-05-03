"""Run the P1 alpha active-contour gate harness end-to-end and persist
artifacts to a UTC-timestamped directory under ``runs/``.

Produces, per gate-test, the visual deliverable bundle defined in
``docs/v2_p1_active_contour_sanity_gate.md`` §8: HDF5 frame_dump
sequence at every K steps, per-frame PNG/HTML, diagnostic 4-panel
plot, harness-generated ``summary.html`` with the §8 status table,
``metadata.json`` mirroring the same status payload, and
(Test 4 only) ``reference_R_star_*.png`` overlay against the
analytic equilibrium circle.

Usage (from repo root, .venv-collab Python):

    /Users/sw1/ActiveCellSim/.venv-collab/bin/python3 \
        scripts/run_p1_alpha_gate.py

The script does not commit its outputs. ``runs/`` is intentionally
excluded from version control (artifacts are large and reproducible).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys

import numpy as np


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
        help="Repository root used for git commit lookup and runs/ output.",
    )
    parser.add_argument(
        "--frame-interval",
        type=int,
        default=50,
        help="Per-test K (write a frame every K steps). Default 50 to keep "
        "the example artifact bundle compact while still producing a "
        "multi-frame sequence.",
    )
    args = parser.parse_args()
    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    from acs.v2.active_contour import (  # noqa: E402
        ActiveContourParameters,
        ActiveContourState,
        ellipse_polygon_vertices,
        regular_polygon_vertices,
    )
    from acs.v2.active_contour_harness import (  # noqa: E402
        compute_finite_n_residual,
        equilibrium_radius_from_cubic,
        run_active_contour_test,
    )

    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%MZ")
    run_root = os.path.join(repo, "runs", f"{timestamp}_p1_alpha")
    os.makedirs(run_root, exist_ok=True)
    git_hash = _git_commit_hash(repo)

    summary_index = []

    # Test 1: zero-force baseline
    n1, r1 = 8, 1.5
    params1 = ActiveContourParameters(
        k_a_nN_per_um=0.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=n1,
        lambda_c_nN=0.0,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    state1 = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n1, radius_um=r1),
        params=params1,
    )
    run1 = run_active_contour_test(
        "test_1_zero_force",
        state1,
        os.path.join(run_root, "test_1_zero_force"),
        n_steps=50,
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
    )
    summary_index.append(
        {"test": "test_1_zero_force", "summary": run1.summary_path}
    )

    # Test 2: area-only — start above target, drive inward
    n2, r2 = 16, 2.0
    params2 = ActiveContourParameters(
        k_a_nN_per_um=1.0,
        target_area_um2=4.0,
        dt_cell_s=1e-3,
        n_vertices=n2,
        lambda_c_nN=0.0,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    state2 = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n2, radius_um=r2),
        params=params2,
    )
    run2 = run_active_contour_test(
        "test_2_area_only",
        state2,
        os.path.join(run_root, "test_2_area_only"),
        n_steps=200,
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
    )
    summary_index.append(
        {"test": "test_2_area_only", "summary": run2.summary_path}
    )

    # Test 3: cortex-only — regular polygon shrinks self-similarly
    n3, r3 = 16, 2.0
    params3 = ActiveContourParameters(
        k_a_nN_per_um=0.0,
        target_area_um2=10.0,
        dt_cell_s=1e-3,
        n_vertices=n3,
        lambda_c_nN=1.0,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    state3 = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n3, radius_um=r3),
        params=params3,
    )
    run3 = run_active_contour_test(
        "test_3_cortex_only",
        state3,
        os.path.join(run_root, "test_3_cortex_only"),
        n_steps=100,
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
    )
    summary_index.append(
        {"test": "test_3_cortex_only", "summary": run3.summary_path}
    )

    # Test 4: coupled ellipse — diagnostic only
    n4 = 64
    semi_a, semi_b = 1.5, 1.0
    target_area = float(np.pi * semi_a * semi_b)
    params4 = ActiveContourParameters(
        k_a_nN_per_um=1.0,
        target_area_um2=target_area,
        dt_cell_s=1e-3,
        n_vertices=n4,
        lambda_c_nN=0.05,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    state4 = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=ellipse_polygon_vertices(
            n4, semi_a_um=semi_a, semi_b_um=semi_b
        ),
        params=params4,
    )
    r_star = equilibrium_radius_from_cubic(
        lambda_c_nN=params4.lambda_c_resolved_nN,
        k_a_nN_per_um=params4.k_a_nN_per_um,
        target_area_um2=target_area,
    )
    reference_state = ActiveContourState(
        cell_id="cell-A-ref",
        vertices_xy_um=regular_polygon_vertices(n4, radius_um=r_star),
        params=params4,
    )
    finite_n_residual = compute_finite_n_residual(reference_state)
    run4 = run_active_contour_test(
        "test_4_coupled_ellipse",
        state4,
        os.path.join(run_root, "test_4_coupled_ellipse"),
        n_steps=200,
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
        extra_metrics={
            "r_star_um": r_star,
            "finite_n_residual_nN": finite_n_residual,
        },
    )
    summary_index.append(
        {"test": "test_4_coupled_ellipse", "summary": run4.summary_path}
    )

    index_path = os.path.join(run_root, "index.json")
    with open(index_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "git_commit_hash": git_hash,
                "frame_interval": args.frame_interval,
                "tests": summary_index,
                "run_root": run_root,
            },
            fh,
            indent=2,
        )
    print(f"P1 alpha gate run complete. artifacts: {run_root}")
    print(f"index: {index_path}")
    for entry in summary_index:
        print(f"  {entry['test']} -> {entry['summary']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
