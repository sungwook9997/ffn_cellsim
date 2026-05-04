"""Run the open-loop ECM preflight harness across the four open-loop
channels and persist artifacts to a UTC-timestamped directory under
``runs/``.

Produces, per scenario, the visual deliverable bundle defined in
``acs.v2.ecm_open_loop_harness``: HDF5 frame_dump sequence at every K
steps, per-frame PNG/HTML, multi-row diagnostic plot, harness-generated
``summary.html`` with the active-channels + per-channel field extrema
status table, and ``metadata.json`` mirroring the same status payload.
A top-level ``index.json`` lists every scenario summary so PI can open
one HTML index and walk the bundle.

Scenarios are caller-supplied test inputs only — no model defaults, no
biological rates, no closed-loop feedback. The script exists to make
the harness's plumbing visible to PI in browser, not to claim a
physical interpretation.

Usage (from repo root, .venv-collab Python):

    /Users/sw1/ActiveCellSim/.venv-collab/bin/python3 \
        scripts/run_ecm_ol_harness.py

The script does not commit its outputs. ``runs/`` is intentionally
excluded from version control via ``.gitignore`` (``*.h5`` / large
artefacts are reproducible from this script + the harness module).
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


def _uniform_grid_factory(shape, value: float):
    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        return np.full(shape, float(value), dtype=np.float64)
    return _factory


def _ramp_grid_factory(shape, low: float, high: float):
    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        return np.linspace(
            float(low), float(high), int(shape[0] * shape[1]), dtype=np.float64
        ).reshape(shape)
    return _factory


def _uniform_diagonal_orientation_rate_factory(shape, diag_value: float):
    """Symmetric diagonal-only orientation rate; choose ``diag_value``
    negative when the initial ECM has identity orientation
    (T_00=T_11=1.0) so post-state stays inside ``|T_ij| <= 1.0``."""

    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        rate = np.zeros((*shape, 2, 2), dtype=np.float64)
        rate[..., 0, 0] = float(diag_value)
        rate[..., 1, 1] = float(diag_value)
        return rate
    return _factory


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
        default=10,
        help="Per-scenario K (write a frame every K steps). Default 10 to "
        "keep the example artifact bundle compact while still producing a "
        "multi-frame sequence.",
    )
    parser.add_argument(
        "--n-steps",
        type=int,
        default=20,
        help="Number of steps per scenario. Default 20 with frame_interval=10 "
        "yields step 0, 10, 20 = 3 frame artifacts per scenario.",
    )
    parser.add_argument(
        "--grid",
        type=int,
        default=4,
        help="ECM grid edge length (nx == ny). Default 4 → 16 cells; tiny "
        "by design so every frame_dump round-trip stays cheap.",
    )
    parser.add_argument(
        "--dt-s",
        type=float,
        default=0.1,
        help="Per-step dt_s [s]. Default 0.1.",
    )
    args = parser.parse_args()
    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    from acs.v2.ecm_open_loop_harness import (  # noqa: E402
        EcmOlScenario,
        make_default_ecm,
        run_ecm_ol_scenario,
    )

    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = os.path.join(repo, "runs", f"{timestamp}_ecm_ol")
    os.makedirs(run_root, exist_ok=False)
    git_hash = _git_commit_hash(repo)

    nx = ny = int(args.grid)
    grid_shape = (nx, ny)
    summary_index: list[dict] = []

    def _record(scenario_obj, run_obj) -> None:
        status = "FAIL" if run_obj.failure is not None else "PASS"
        summary_index.append(
            {
                "scenario": scenario_obj.name,
                "status": status,
                "summary": run_obj.summary_path,
                "metadata": run_obj.metadata_path,
                "failure": run_obj.failure,
            }
        )

    # Scenario 1: traction-only (uniform prescribed traction)
    ecm1 = make_default_ecm(nx=nx, ny=ny)
    scenario1 = EcmOlScenario(
        name="traction_only_uniform",
        initial_ecm=ecm1,
        n_steps=args.n_steps,
        dt_s=args.dt_s,
        traction_factory=_uniform_grid_factory(grid_shape, 0.5),
    )
    run1 = run_ecm_ol_scenario(
        scenario1,
        os.path.join(run_root, scenario1.name),
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
    )
    _record(scenario1, run1)

    # Scenario 2: stiffness-rate-only (positive ramp)
    ecm2 = make_default_ecm(nx=nx, ny=ny)
    scenario2 = EcmOlScenario(
        name="stiffness_rate_only_ramp",
        initial_ecm=ecm2,
        n_steps=args.n_steps,
        dt_s=args.dt_s,
        stiffness_rate_factory=_ramp_grid_factory(grid_shape, 0.1, 1.0),
    )
    run2 = run_ecm_ol_scenario(
        scenario2,
        os.path.join(run_root, scenario2.name),
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
    )
    _record(scenario2, run2)

    # Scenario 3: density-rate-only (ligand grows, fiber shallower growth)
    ecm3 = make_default_ecm(nx=nx, ny=ny)
    scenario3 = EcmOlScenario(
        name="density_rate_only_mixed",
        initial_ecm=ecm3,
        n_steps=args.n_steps,
        dt_s=args.dt_s,
        ligand_density_rate_factory=_uniform_grid_factory(grid_shape, 0.2),
        fiber_density_rate_factory=_uniform_grid_factory(grid_shape, 0.05),
    )
    run3 = run_ecm_ol_scenario(
        scenario3,
        os.path.join(run_root, scenario3.name),
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
    )
    _record(scenario3, run3)

    # Scenario 4: orientation-rate-only (negative diagonal from identity)
    ecm4 = make_default_ecm(nx=nx, ny=ny)
    scenario4 = EcmOlScenario(
        name="orientation_rate_only_negative_diagonal",
        initial_ecm=ecm4,
        n_steps=args.n_steps,
        dt_s=args.dt_s,
        orientation_rate_factory=_uniform_diagonal_orientation_rate_factory(
            grid_shape, -0.05
        ),
    )
    run4 = run_ecm_ol_scenario(
        scenario4,
        os.path.join(run_root, scenario4.name),
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
    )
    _record(scenario4, run4)

    # Scenario 5: all-channel smoke (exercise the sequential application order)
    ecm5 = make_default_ecm(nx=nx, ny=ny)
    scenario5 = EcmOlScenario(
        name="all_channel_smoke",
        initial_ecm=ecm5,
        n_steps=args.n_steps,
        dt_s=args.dt_s,
        traction_factory=_uniform_grid_factory(grid_shape, 0.3),
        stiffness_rate_factory=_uniform_grid_factory(grid_shape, 0.4),
        ligand_density_rate_factory=_uniform_grid_factory(grid_shape, 0.1),
        fiber_density_rate_factory=_uniform_grid_factory(grid_shape, 0.05),
        orientation_rate_factory=_uniform_diagonal_orientation_rate_factory(
            grid_shape, -0.05
        ),
    )
    run5 = run_ecm_ol_scenario(
        scenario5,
        os.path.join(run_root, scenario5.name),
        frame_interval=args.frame_interval,
        git_commit_hash=git_hash,
    )
    _record(scenario5, run5)

    failed_scenarios = [
        entry["scenario"] for entry in summary_index if entry["status"] == "FAIL"
    ]
    aggregate_status = "FAIL" if failed_scenarios else "PASS"
    index_path = os.path.join(run_root, "index.json")
    with open(index_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "git_commit_hash": git_hash,
                "frame_interval": args.frame_interval,
                "n_steps": args.n_steps,
                "grid_shape": list(grid_shape),
                "dt_s": args.dt_s,
                "scenarios": summary_index,
                "aggregate_status": aggregate_status,
                "failed_scenarios": failed_scenarios,
                "run_root": run_root,
            },
            fh,
            indent=2,
        )
    print(f"ECM-OL harness run {aggregate_status}. artifacts: {run_root}")
    print(f"index: {index_path}")
    for entry in summary_index:
        print(f"  [{entry['status']}] {entry['scenario']} -> {entry['summary']}")
    if failed_scenarios:
        print(f"FAILED scenarios: {', '.join(failed_scenarios)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
