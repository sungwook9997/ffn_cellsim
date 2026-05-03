"""Active contour test/gate harness: simulate + diagnostics + visual deliverables.

Glues together the schema (`acs.v2.active_contour`), dynamics
(`acs.v2.dynamics.active_contour`), HDF5 frame writer
(`acs.v2.output.frame_dump`), and stub3d renderer
(`acs.v2.viz.stub3d`) into the four gate tests defined in
``docs/v2_p1_active_contour_sanity_gate.md`` §4 / §8 and lock §6.

The harness produces every artifact §8 lists for each test, on PASS
or FAIL. Quantitative gate decisions live in the tests, not the
harness; the harness only records the diagnostics they consume.

This module is plumbing only — no new physics or magic numbers. The
two integer knobs (``frame_interval`` for HDF5/PNG sampling and
``n_steps`` for run length) are caller-supplied, not module
defaults.
"""

from __future__ import annotations

import json
import os
import time
import traceback
from dataclasses import dataclass, field
from typing import Optional

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from acs.v2.active_contour import ActiveContourState, compute_rate_max  # noqa: E402
from acs.v2.cell_cluster import CellClusterState  # noqa: E402
from acs.v2.dynamics.active_contour import (  # noqa: E402
    ActiveContourStepError,
    compute_area_forces,
    compute_cortex_forces,
    energy_area,
    energy_cortex,
    step,
)
from acs.v2.ecm_substrate import ECMSubstrateState  # noqa: E402
from acs.v2.measurement_boundary import MeasurementBoundaryError  # noqa: E402
from acs.v2.output.frame_dump import write_frame  # noqa: E402
from acs.v2.single_cell import SingleCellState  # noqa: E402
from acs.v2.viz.stub3d import render_frame_html, render_frame_png  # noqa: E402


@dataclass
class StepDiagnostics:
    step_index: int
    time_s: float
    energy_cortex_nNum: float
    energy_area_nNum: float
    force_norm_nN: float
    margin_ratio: float
    area_um2: float
    perimeter_um: float
    anisotropy: float


@dataclass
class HarnessRun:
    test_name: str
    output_dir: str
    history: list[StepDiagnostics] = field(default_factory=list)
    final_state: Optional[ActiveContourState] = None
    final_step_index: int = 0
    frame_interval: int = 100
    n_steps: int = 0
    wall_clock_s: float = 0.0
    failure: Optional[dict] = None
    frame_paths: list[str] = field(default_factory=list)
    png_paths: list[str] = field(default_factory=list)
    html_paths: list[str] = field(default_factory=list)
    diagnostic_plot_path: Optional[str] = None
    summary_path: Optional[str] = None
    failure_report_path: Optional[str] = None
    metadata_path: Optional[str] = None
    reference_plot_path: Optional[str] = None
    extra_metrics: dict = field(default_factory=dict)


def _placeholder_ecm() -> ECMSubstrateState:
    """Tiny ECM grid the frame_dump writer needs for cluster validation."""

    nx, ny = 2, 2
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    return ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=np.zeros((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
    )


def _record_diagnostics(state: ActiveContourState, step_index: int) -> StepDiagnostics:
    rate_info = compute_rate_max(state)
    forces = compute_cortex_forces(state) + compute_area_forces(state)
    f_norm = float(np.linalg.norm(forces))
    centroid = state.vertices_xy_um.mean(axis=0)
    radii = np.linalg.norm(state.vertices_xy_um - centroid, axis=1)
    r_mean = float(radii.mean()) if radii.size else 0.0
    anisotropy = (
        float(radii.max() - radii.min()) / r_mean if r_mean > 0.0 else 0.0
    )
    return StepDiagnostics(
        step_index=step_index,
        time_s=float(step_index) * state.params.dt_cell_s,
        energy_cortex_nNum=energy_cortex(state),
        energy_area_nNum=energy_area(state),
        force_norm_nN=f_norm,
        margin_ratio=float(rate_info["margin_ratio"]),
        area_um2=state.area_um2(),
        perimeter_um=state.perimeter_um(),
        anisotropy=anisotropy,
    )


def _write_frame_artifacts(
    state: ActiveContourState, run: HarnessRun, step_index: int
) -> None:
    boundary = state.to_measurement_boundary(source_modality="active_contour_harness")
    cell = SingleCellState(
        cell_id=state.cell_id,
        time_s=float(step_index) * state.params.dt_cell_s,
        measurement_boundary=boundary,
        cell_state="alive",
    )
    cluster = CellClusterState(
        cells={state.cell_id: cell},
        ecm=_placeholder_ecm(),
    )
    frame_path = os.path.join(run.output_dir, f"frame_{step_index:06d}.h5")
    png_path = os.path.join(run.output_dir, f"frame_{step_index:06d}.png")
    html_path = os.path.join(run.output_dir, f"frame_{step_index:06d}.html")
    write_frame(frame_path, cluster, time_s=float(step_index) * state.params.dt_cell_s)
    render_frame_png(frame_path, png_path)
    render_frame_html(frame_path, html_path)
    run.frame_paths.append(frame_path)
    run.png_paths.append(png_path)
    run.html_paths.append(html_path)


def _make_diagnostic_plot(run: HarnessRun) -> Optional[str]:
    if not run.history:
        return None
    steps = np.array([d.step_index for d in run.history])
    e_c = np.array([d.energy_cortex_nNum for d in run.history])
    e_a = np.array([d.energy_area_nNum for d in run.history])
    f_norm = np.array([d.force_norm_nN for d in run.history])
    margin = np.array([d.margin_ratio for d in run.history])
    aniso = np.array([d.anisotropy for d in run.history])

    fig, axes = plt.subplots(4, 1, figsize=(6.0, 8.0), sharex=True)
    axes[0].plot(steps, e_c + e_a, label="E_c + E_A")
    axes[0].plot(steps, e_c, label="E_c", linestyle="--")
    axes[0].plot(steps, e_a, label="E_A", linestyle=":")
    axes[0].set_ylabel("energy (nN·μm)")
    axes[0].legend(loc="best", fontsize="x-small")
    axes[1].plot(steps, f_norm)
    axes[1].set_ylabel("||F||_2 (nN)")
    axes[2].plot(steps, margin)
    axes[2].axhline(1.0, color="red", linestyle="--", linewidth=0.5)
    axes[2].set_ylabel("margin_ratio")
    axes[3].plot(steps, aniso)
    axes[3].set_ylabel("anisotropy")
    axes[3].set_xlabel("step index")
    fig.suptitle(f"{run.test_name} diagnostics")
    fig.tight_layout()
    plot_path = os.path.join(run.output_dir, f"diagnostic_{run.test_name}.png")
    fig.savefig(plot_path, dpi=120)
    plt.close(fig)
    return plot_path


def _build_status_payload(
    run: HarnessRun, status: str, git_commit_hash: str
) -> dict:
    """Build the gate §8 status payload (also written to metadata.json)."""

    final_state = run.final_state
    final_area = final_state.area_um2() if final_state is not None else 0.0
    target_area = (
        final_state.params.target_area_um2 if final_state is not None else 0.0
    )
    final_residual = run.history[-1].force_norm_nN if run.history else 0.0
    finite_n_residual = run.extra_metrics.get("finite_n_residual_nN", None)
    return {
        "test_name": run.test_name,
        "status": status,
        "n_steps": run.n_steps,
        "final_step_index": run.final_step_index,
        "frame_interval": run.frame_interval,
        "wall_clock_s": run.wall_clock_s,
        "git_commit_hash": git_commit_hash or "unrecorded",
        "final_residual_nN": final_residual,
        "finite_N_residual_nN": finite_n_residual,
        "final_area_um2": final_area,
        "target_area_um2": target_area,
        "extra_metrics": run.extra_metrics,
        "failure": run.failure,
    }


def _write_metadata_json(run: HarnessRun, payload: dict) -> str:
    metadata_path = os.path.join(run.output_dir, "metadata.json")
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
    return metadata_path


def _make_reference_r_star_plot(run: HarnessRun) -> Optional[str]:
    """Test 4 reference plot: final polygon overlaid with the analytic
    equilibrium circle of radius `r_star_um` (recorded in extra_metrics)."""

    if run.final_state is None:
        return None
    r_star = run.extra_metrics.get("r_star_um")
    if r_star is None:
        return None
    centroid = run.final_state.vertices_xy_um.mean(axis=0)
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    verts = run.final_state.vertices_xy_um
    closed = np.vstack([verts, verts[:1]])
    ax.plot(closed[:, 0], closed[:, 1], color="#0d3d66", linewidth=1.4, label="final polygon")
    theta = np.linspace(0.0, 2.0 * np.pi, 256)
    ax.plot(
        centroid[0] + r_star * np.cos(theta),
        centroid[1] + r_star * np.sin(theta),
        color="red",
        linestyle="--",
        linewidth=1.0,
        label=f"R* = {r_star:.4f} μm",
    )
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="best", fontsize="x-small")
    ax.set_title(f"{run.test_name} reference equilibrium")
    fig.tight_layout()
    plot_path = os.path.join(run.output_dir, f"reference_R_star_{run.test_name}.png")
    fig.savefig(plot_path, dpi=120)
    plt.close(fig)
    return plot_path


def _make_summary_html(
    run: HarnessRun, *, status: str, git_commit_hash: str = "", payload: dict
) -> str:
    rows = "".join(
        f'<tr><td>{d.step_index}</td><td>{d.time_s:.4e}</td>'
        f'<td>{d.energy_cortex_nNum + d.energy_area_nNum:.4e}</td>'
        f'<td>{d.force_norm_nN:.4e}</td><td>{d.margin_ratio:.4e}</td>'
        f'<td>{d.area_um2:.4e}</td><td>{d.anisotropy:.4e}</td></tr>'
        for d in run.history
    )
    thumbs = "".join(
        f'<a href="{os.path.basename(p)}"><img src="{os.path.basename(p)}" '
        f'width="120" alt="frame"/></a>'
        for p in run.png_paths
    )
    diag_block = ""
    if run.diagnostic_plot_path:
        diag_block = (
            f'<h2>Diagnostic plot</h2>'
            f'<img src="{os.path.basename(run.diagnostic_plot_path)}" />'
        )
    reference_block = ""
    if run.reference_plot_path:
        reference_block = (
            f'<h2>Reference equilibrium</h2>'
            f'<img src="{os.path.basename(run.reference_plot_path)}" />'
        )
    failure_block = ""
    if run.failure is not None:
        failure_block = (
            f'<h2>Failure</h2><pre>{json.dumps(run.failure, indent=2)}</pre>'
        )
    status_table_rows = []
    for key in (
        "test_name",
        "status",
        "n_steps",
        "final_step_index",
        "frame_interval",
        "wall_clock_s",
        "git_commit_hash",
        "final_residual_nN",
        "finite_N_residual_nN",
        "final_area_um2",
        "target_area_um2",
    ):
        value = payload.get(key)
        status_table_rows.append(f"<tr><th>{key}</th><td>{value}</td></tr>")
    if payload.get("extra_metrics"):
        status_table_rows.append(
            f'<tr><th>extra_metrics</th><td><pre>'
            f'{json.dumps(payload["extra_metrics"], indent=2, default=str)}'
            f'</pre></td></tr>'
        )
    status_table = "".join(status_table_rows)
    body = (
        f"<!doctype html><html><head><meta charset=\"utf-8\">"
        f"<title>{run.test_name} summary</title>"
        f"<style>body{{font-family:sans-serif;margin:1em;}}"
        f"table{{border-collapse:collapse;margin-top:1em;}}"
        f"td,th{{border:1px solid #999;padding:0.25em 0.5em;text-align:left;}}"
        f"img{{margin:0.25em;}}"
        f"pre{{margin:0;}}"
        f"</style></head><body>"
        f"<h1>{run.test_name} summary</h1>"
        f"<h2>Status</h2><table><tbody>{status_table}</tbody></table>"
        f"{diag_block}"
        f"{reference_block}"
        f"{failure_block}"
        f"<h2>Frame thumbnails</h2><div>{thumbs}</div>"
        f"<h2>Per-step metrics</h2>"
        f"<table><thead><tr><th>step</th><th>time_s</th><th>E_total</th>"
        f"<th>||F||_2</th><th>margin_ratio</th><th>area_um2</th>"
        f"<th>anisotropy</th></tr></thead>"
        f"<tbody>{rows}</tbody></table></body></html>"
    )
    summary_path = os.path.join(run.output_dir, "summary.html")
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(body)
    return summary_path


def _write_failure_report(run: HarnessRun) -> str:
    last_state = run.history[-1] if run.history else None
    payload = {
        "test_name": run.test_name,
        "final_step_index": run.final_step_index,
        "failure": run.failure,
        "last_metrics": last_state.__dict__ if last_state is not None else None,
        "frame_paths": [os.path.basename(p) for p in run.frame_paths],
        "diagnostic_plot": (
            os.path.basename(run.diagnostic_plot_path)
            if run.diagnostic_plot_path
            else None
        ),
    }
    body = "# Active contour harness — failure report\n\n```json\n"
    body += json.dumps(payload, indent=2, default=str)
    body += "\n```\n"
    path = os.path.join(run.output_dir, "failure_report.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(body)
    return path


def _require_non_negative_int(value, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be int, got {type(value).__name__}")
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value!r}")
    return value


def _require_positive_int(value, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be int, got {type(value).__name__}")
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer, got {value!r}")
    return value


def run_active_contour_test(
    test_name: str,
    initial_state: ActiveContourState,
    output_dir: str,
    *,
    n_steps: int,
    frame_interval: int = 100,
    git_commit_hash: str = "",
    extra_metrics: Optional[dict] = None,
) -> HarnessRun:
    """Run one Sanity Gate test and produce all §8 visual artifacts.

    Parameters:
        test_name: Used for filenames and the summary heading.
        initial_state: Validated :class:`ActiveContourState`.
        output_dir: Per-test output directory; created if missing.
        n_steps: Non-negative integer count of overdamped Euler steps.
            Bools and floats are rejected.
        frame_interval: Positive integer K. Write a frame_dump + PNG + HTML
            every K steps (plus step 0 and the final step). Default 100
            per gate §8 (Tests 2/3/4 multi-frame sequence). Bools and
            floats are rejected.
        git_commit_hash: Optional value embedded in summary.html and
            metadata.json for reproducibility tracking.
        extra_metrics: Optional pre-computed metrics (e.g.
            ``finite_n_residual_nN`` and ``r_star_um`` for Test 4)
            recorded in the run, persisted in metadata.json, and shown
            in the summary status table.

    Returns the populated :class:`HarnessRun`. On failure the run's
    ``failure`` field is set, ``failure_report.md`` is written, and
    the harness returns (does not raise) so the caller can decide
    PASS / DIAGNOSTIC / BLOCKER from the recorded state.
    """

    _require_non_negative_int(n_steps, "n_steps")
    _require_positive_int(frame_interval, "frame_interval")
    os.makedirs(output_dir, exist_ok=True)
    run = HarnessRun(
        test_name=test_name,
        output_dir=output_dir,
        frame_interval=frame_interval,
        n_steps=n_steps,
        extra_metrics=dict(extra_metrics) if extra_metrics else {},
    )
    wall_clock_start = time.perf_counter()
    state = initial_state
    run.history.append(_record_diagnostics(state, step_index=0))
    try:
        _write_frame_artifacts(state, run, step_index=0)
    except Exception as exc:  # pragma: no cover - I/O failure is contract-level
        run.failure = {
            "kind": "frame_write_failure",
            "step_index": 0,
            "exception_type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        run.failure_report_path = _write_failure_report(run)
        return run

    last_step_emitted = 0
    for i in range(1, n_steps + 1):
        try:
            state = step(state)
        except ActiveContourStepError as exc:
            run.failure = {
                "kind": "dt_violation",
                "step_index": i,
                "failure_kind": exc.failure_kind,
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
            break
        except MeasurementBoundaryError as exc:
            run.failure = {
                "kind": "measurement_boundary_violation",
                "step_index": i,
                "failure_kind": exc.failure_kind,
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
            break
        run.history.append(_record_diagnostics(state, step_index=i))
        if i % frame_interval == 0 or i == n_steps:
            try:
                _write_frame_artifacts(state, run, step_index=i)
                last_step_emitted = i
            except Exception as exc:  # pragma: no cover
                run.failure = {
                    "kind": "frame_write_failure",
                    "step_index": i,
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                    "traceback": traceback.format_exc(),
                }
                break
    else:
        last_step_emitted = n_steps

    run.final_state = state
    run.final_step_index = state.step_count
    run.wall_clock_s = float(time.perf_counter() - wall_clock_start)
    run.diagnostic_plot_path = _make_diagnostic_plot(run)
    run.reference_plot_path = _make_reference_r_star_plot(run)
    status = "FAIL" if run.failure is not None else "PASS"
    payload = _build_status_payload(run, status, git_commit_hash)
    run.metadata_path = _write_metadata_json(run, payload)
    run.summary_path = _make_summary_html(
        run, status=status, git_commit_hash=git_commit_hash, payload=payload
    )
    if run.failure is not None:
        run.failure_report_path = _write_failure_report(run)
    return run


def compute_finite_n_residual(state: ActiveContourState) -> float:
    """Helper for Test 4 diagnostic: ``||F_cortex + F_area||`` evaluated
    on the supplied state."""

    forces = compute_cortex_forces(state) + compute_area_forces(state)
    return float(np.linalg.norm(forces))


def equilibrium_radius_from_cubic(
    *, lambda_c_nN: float, k_a_nN_per_um: float, target_area_um2: float
) -> float:
    """Solve the force-balance cubic in radius via :func:`numpy.roots`.

    Starting from ``K_A · (πR² − A_0)/A_0 + λ_c/R = 0`` and multiplying
    by ``R`` gives the cubic
    ``(K_A · π / A_0) · R³ − K_A · R + λ_c = 0``. (The original
    fraction is degree-3 in ``R`` after clearing the ``1/R`` term, not
    quartic.) Returns the positive real root closest to
    ``sqrt(A_0 / π)``. Raises ``ValueError`` when no positive real
    root exists in the expected range.
    """

    if lambda_c_nN <= 0.0 or k_a_nN_per_um <= 0.0 or target_area_um2 <= 0.0:
        raise ValueError(
            "equilibrium radius is only defined for positive lambda_c, K_A, "
            "and target_area"
        )
    a = k_a_nN_per_um * np.pi / target_area_um2  # coefficient on R^3
    b = -k_a_nN_per_um  # coefficient on R
    c = lambda_c_nN  # constant
    coeffs = [a, 0.0, b, c]  # cubic: a*R^3 + 0*R^2 + b*R + c
    roots = np.roots(coeffs)
    real_positive = [
        float(r.real) for r in roots if abs(r.imag) < 1e-9 and r.real > 0.0
    ]
    if not real_positive:
        raise ValueError(f"no positive real root for coefficients {coeffs!r}")
    r_init_guess = float(np.sqrt(target_area_um2 / np.pi))
    return min(real_positive, key=lambda r: abs(r - r_init_guess))


# Backward-compat alias (deprecated): the public name was briefly
# `equilibrium_radius_from_quartic` in the immediate-prior commit; the
# math is cubic, not quartic. The alias keeps existing callers working
# while the renamed function is the canonical export.
equilibrium_radius_from_quartic = equilibrium_radius_from_cubic
