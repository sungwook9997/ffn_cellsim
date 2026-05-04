"""Open-loop ECM preflight harness: drive the four open-loop preflight
paths over scenarios + persist artifacts.

Composes the already-tested open-loop preflight functions in
:mod:`acs.v2.dynamics.ecm_open_loop`:

- :func:`accumulate_prescribed_traction`
- :func:`apply_prescribed_stiffness_rate`
- :func:`apply_prescribed_density_rate`
- :func:`apply_prescribed_orientation_rate`

into a multi-step / multi-scenario sanity-test driver. A scenario
defines, per channel, an optional caller-supplied factory that
returns the *prescribed* field at a given integer step (no closed-loop
feedback). Per step the harness applies the active channels in a
fixed order ``traction → stiffness_rate → density_rate →
orientation_rate``, records HDF5 frame dumps, per-step diagnostics,
a multi-row diagnostic plot, a status table in ``summary.html``, and
a JSON metadata payload that mirrors the same status table — same
artifact contract as the P1 alpha gate harness.

This module is plumbing only — no new physics, no new tunable
constants, no new biological default rates. Caller-supplied
prescribed factories are *test inputs*, not model defaults; the
harness records which channels were active under
``metadata['active_channels']`` so a downstream comparison can
re-run with a different prescribed pattern without changing this
file. Sequential application is exactly four library calls per
step; there is no inter-channel coupling beyond the order in which
they execute.

Sanity Gate scope (ecm_open_loop_harness):

- §1 dimensional: every channel application uses an
  already-verified open-loop preflight whose Rule 10 unit chain was
  established in ``acs.v2.dynamics.ecm_open_loop`` (traction
  ``[nN/μm²]·[s]=[nN·s/μm²]``; stiffness ``[kPa/s]·[s]=[kPa]``;
  density ``[1/s]·[s]=[dimensionless]``; orientation
  ``[1/s]·[s]=[dimensionless]``). The harness only records per-step
  scalars (max / min / mean of the relevant field) — it introduces
  no new unit reduction.
- §2 boundary cases: at least one channel factory is required (a
  scenario with no factory raises ``ValueError`` at construction);
  ECM ``output_dir`` is required; ``n_steps`` is non-negative
  integer; ``frame_interval`` is positive integer. Bool/float are
  rejected for both knobs.
- §3 conservation/provenance: the harness threads a single ECM
  state through the four channel calls per step without ever
  mutating an earlier step's state. Each ``frame_dump`` persists a
  ``CellClusterState`` whose ``ECMSubstrateState`` field is the
  canonical evolved state, so a downstream ``read_frame``
  round-trip materialises the same ECM that drove the next step.
- §4 numerical sanity: float64 throughout via the underlying
  preflight functions.
- §5 sign/sense: prescribed traction is non-negative by the
  open-loop preflight contract; stiffness / density rates are
  signed and may push fields into invalid bounds — the harness
  records the resulting ``ECMOpenLoopError.failure_kind`` rather
  than clamping. Orientation rates are signed and must be
  symmetric; the schema bound check is component-wise.
- §6 measurement-protocol consistency: per Hard Rule 11, the
  harness's measurement modality is the per-step post-state ECM
  field grids persisted by ``frame_dump`` and re-rendered by
  ``stub3d``; no new measurement is invented. The diagnostic plot
  panels read the same per-step scalars that ``metadata.json``
  reports.

Magic-Number Block: this module declares no tunable numeric. All
scenario parameters (traction amplitudes, rates, n_steps, dt_s,
frame_interval) are caller-supplied test inputs.
"""

from __future__ import annotations

import json
import os
import time
import traceback
from dataclasses import dataclass, field
from typing import Callable, Optional

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from acs.v2.cell_cluster import CellClusterState  # noqa: E402
from acs.v2.dynamics.ecm_open_loop import (  # noqa: E402
    ECMOpenLoopError,
    accumulate_prescribed_traction,
    apply_prescribed_density_rate,
    apply_prescribed_orientation_rate,
    apply_prescribed_stiffness_rate,
)
from acs.v2.ecm_substrate import ECMSubstrateState  # noqa: E402
from acs.v2.measurement_boundary import MeasurementBoundary  # noqa: E402
from acs.v2.output.frame_dump import write_frame  # noqa: E402
from acs.v2.single_cell import SingleCellState  # noqa: E402
from acs.v2.viz.stub3d import render_frame_html, render_frame_png  # noqa: E402


@dataclass
class EcmOlStepDiagnostics:
    step_index: int
    time_s: float
    accumulated_max_nNs_per_um2: float
    accumulated_mean_nNs_per_um2: float
    accumulated_sum_nNs_per_um2: float
    nonzero_fraction: float
    traction_max_nN_per_um2: float
    traction_mean_nN_per_um2: float
    stiffness_kpa_max: float
    stiffness_kpa_min: float
    ligand_density_max: float
    ligand_density_min: float
    fiber_density_max: float
    fiber_density_min: float
    orientation_tensor_abs_max: float


@dataclass
class EcmOlScenario:
    """One prescribed open-loop scenario for the ECM-OL preflight.

    At least one of the channel factories must be supplied. Each
    factory takes a 1-based integer step index and returns the
    caller-supplied prescribed field for that step in the units
    contracted by the underlying preflight function:

    - ``traction_factory(i) -> (nx, ny) float64`` ``[nN/μm²]``,
      non-negative.
    - ``stiffness_rate_factory(i) -> (nx, ny) float64`` ``[kPa/s]``,
      signed.
    - ``ligand_density_rate_factory(i) -> (nx, ny) float64``
      ``[1/s]``, signed.
    - ``fiber_density_rate_factory(i) -> (nx, ny) float64``
      ``[1/s]``, signed.
    - ``orientation_rate_factory(i) -> (nx, ny, 2, 2) float64``
      ``[1/s]``, symmetric component-wise.

    Step 0 records only the initial state; channel factories are
    called for steps ``1..n_steps``.
    """

    name: str
    initial_ecm: ECMSubstrateState
    n_steps: int
    dt_s: float
    traction_factory: Optional[Callable[[int], np.ndarray]] = None
    stiffness_rate_factory: Optional[Callable[[int], np.ndarray]] = None
    ligand_density_rate_factory: Optional[Callable[[int], np.ndarray]] = None
    fiber_density_rate_factory: Optional[Callable[[int], np.ndarray]] = None
    orientation_rate_factory: Optional[Callable[[int], np.ndarray]] = None

    def __post_init__(self) -> None:
        if all(
            f is None
            for f in (
                self.traction_factory,
                self.stiffness_rate_factory,
                self.ligand_density_rate_factory,
                self.fiber_density_rate_factory,
                self.orientation_rate_factory,
            )
        ):
            raise ValueError(
                "EcmOlScenario requires at least one channel factory "
                "(traction_factory, stiffness_rate_factory, "
                "ligand_density_rate_factory, fiber_density_rate_factory, "
                "or orientation_rate_factory)"
            )

    def active_channels(self) -> list[str]:
        active: list[str] = []
        if self.traction_factory is not None:
            active.append("traction")
        if self.stiffness_rate_factory is not None:
            active.append("stiffness_rate")
        if self.ligand_density_rate_factory is not None:
            active.append("ligand_density_rate")
        if self.fiber_density_rate_factory is not None:
            active.append("fiber_density_rate")
        if self.orientation_rate_factory is not None:
            active.append("orientation_rate")
        return active


@dataclass
class EcmOlRun:
    scenario_name: str
    output_dir: str
    history: list[EcmOlStepDiagnostics] = field(default_factory=list)
    final_ecm: Optional[ECMSubstrateState] = None
    final_step_index: int = 0
    n_steps: int = 0
    dt_s: float = 0.0
    frame_interval: int = 1
    wall_clock_s: float = 0.0
    failure: Optional[dict] = None
    frame_paths: list[str] = field(default_factory=list)
    png_paths: list[str] = field(default_factory=list)
    html_paths: list[str] = field(default_factory=list)
    diagnostic_plot_path: Optional[str] = None
    summary_path: Optional[str] = None
    metadata_path: Optional[str] = None
    failure_report_path: Optional[str] = None
    active_channels: list[str] = field(default_factory=list)


def _placeholder_cell_state(ecm: ECMSubstrateState, time_s: float) -> SingleCellState:
    """Tiny CCW square at the ECM grid origin so frame_dump's cluster
    contract has a cell to wrap around the evolving ECM."""

    nx, ny = ecm.grid_shape
    extent_x = ecm.spacing_um * nx
    extent_y = ecm.spacing_um * ny
    side = max(0.1 * min(extent_x, extent_y), ecm.spacing_um * 0.5)
    ox, oy = ecm.origin_um_xy
    cx = ox + 0.5 * extent_x
    cy = oy + 0.5 * extent_y
    half = 0.5 * side
    boundary = MeasurementBoundary.from_array(
        np.array(
            [
                [cx - half, cy - half],
                [cx + half, cy - half],
                [cx + half, cy + half],
                [cx - half, cy + half],
            ],
            dtype=np.float64,
        ),
        coordinate_convention="world_um_y_up",
        source_modality="ecm_ol_harness_placeholder",
    )
    return SingleCellState(
        cell_id="ecm-ol-placeholder",
        time_s=time_s,
        measurement_boundary=boundary,
        cell_state="alive",
    )


def _record_diagnostics(
    ecm: ECMSubstrateState,
    traction: Optional[np.ndarray],
    step_index: int,
    dt_s: float,
) -> EcmOlStepDiagnostics:
    accumulated = np.asarray(ecm.accumulated_traction_nNs_per_um2, dtype=np.float64)
    nonzero = float((accumulated > 0.0).mean()) if accumulated.size else 0.0
    stiffness = np.asarray(ecm.stiffness_kpa, dtype=np.float64)
    ligand = np.asarray(ecm.ligand_density, dtype=np.float64)
    fiber = np.asarray(ecm.fiber_density, dtype=np.float64)
    orientation = np.asarray(ecm.orientation_tensor, dtype=np.float64)
    traction_arr = (
        np.asarray(traction, dtype=np.float64)
        if traction is not None
        else np.zeros(ecm.grid_shape, dtype=np.float64)
    )
    return EcmOlStepDiagnostics(
        step_index=step_index,
        time_s=float(step_index) * float(dt_s),
        accumulated_max_nNs_per_um2=float(accumulated.max()) if accumulated.size else 0.0,
        accumulated_mean_nNs_per_um2=float(accumulated.mean()) if accumulated.size else 0.0,
        accumulated_sum_nNs_per_um2=float(accumulated.sum()),
        nonzero_fraction=nonzero,
        traction_max_nN_per_um2=float(traction_arr.max()) if traction_arr.size else 0.0,
        traction_mean_nN_per_um2=float(traction_arr.mean()) if traction_arr.size else 0.0,
        stiffness_kpa_max=float(stiffness.max()) if stiffness.size else 0.0,
        stiffness_kpa_min=float(stiffness.min()) if stiffness.size else 0.0,
        ligand_density_max=float(ligand.max()) if ligand.size else 0.0,
        ligand_density_min=float(ligand.min()) if ligand.size else 0.0,
        fiber_density_max=float(fiber.max()) if fiber.size else 0.0,
        fiber_density_min=float(fiber.min()) if fiber.size else 0.0,
        orientation_tensor_abs_max=(
            float(np.abs(orientation).max()) if orientation.size else 0.0
        ),
    )


def _write_frame_artifacts(
    ecm: ECMSubstrateState, run: EcmOlRun, step_index: int, time_s: float
) -> None:
    cell = _placeholder_cell_state(ecm, time_s)
    cluster = CellClusterState(cells={cell.cell_id: cell}, ecm=ecm)
    frame_path = os.path.join(run.output_dir, f"frame_{step_index:06d}.h5")
    png_path = os.path.join(run.output_dir, f"frame_{step_index:06d}.png")
    html_path = os.path.join(run.output_dir, f"frame_{step_index:06d}.html")
    write_frame(frame_path, cluster, time_s=time_s)
    render_frame_png(frame_path, png_path)
    render_frame_html(frame_path, html_path)
    run.frame_paths.append(frame_path)
    run.png_paths.append(png_path)
    run.html_paths.append(html_path)


def _make_diagnostic_plot(run: EcmOlRun) -> Optional[str]:
    if not run.history:
        return None
    steps = np.array([d.step_index for d in run.history])
    rows = [
        ("traction max (nN/μm²)", [d.traction_max_nN_per_um2 for d in run.history]),
        ("accum max (nN·s/μm²)", [d.accumulated_max_nNs_per_um2 for d in run.history]),
        ("stiffness max (kPa)", [d.stiffness_kpa_max for d in run.history]),
        ("ligand max", [d.ligand_density_max for d in run.history]),
        ("fiber max", [d.fiber_density_max for d in run.history]),
        ("orient |max|", [d.orientation_tensor_abs_max for d in run.history]),
    ]
    fig, axes = plt.subplots(len(rows), 1, figsize=(6.0, 1.6 * len(rows)), sharex=True)
    if len(rows) == 1:
        axes = [axes]
    for ax, (label, series) in zip(axes, rows):
        ax.plot(steps, np.asarray(series, dtype=np.float64))
        ax.set_ylabel(label)
    axes[-1].set_xlabel("step index")
    fig.suptitle(f"{run.scenario_name} ECM-OL diagnostics")
    fig.tight_layout()
    plot_path = os.path.join(run.output_dir, f"diagnostic_{run.scenario_name}.png")
    fig.savefig(plot_path, dpi=120)
    plt.close(fig)
    return plot_path


def _build_status_payload(
    run: EcmOlRun, status: str, git_commit_hash: str
) -> dict:
    final_state = run.final_ecm

    def _max(field_name: str) -> float:
        if final_state is None:
            return 0.0
        arr = np.asarray(getattr(final_state, field_name))
        return float(arr.max()) if arr.size else 0.0

    def _min(field_name: str) -> float:
        if final_state is None:
            return 0.0
        arr = np.asarray(getattr(final_state, field_name))
        return float(arr.min()) if arr.size else 0.0

    def _mean(field_name: str) -> float:
        if final_state is None:
            return 0.0
        arr = np.asarray(getattr(final_state, field_name))
        return float(arr.mean()) if arr.size else 0.0

    def _abs_max(field_name: str) -> float:
        if final_state is None:
            return 0.0
        arr = np.asarray(getattr(final_state, field_name))
        return float(np.abs(arr).max()) if arr.size else 0.0

    return {
        "scenario_name": run.scenario_name,
        "status": status,
        "n_steps": run.n_steps,
        "final_step_index": run.final_step_index,
        "frame_interval": run.frame_interval,
        "dt_s": run.dt_s,
        "wall_clock_s": run.wall_clock_s,
        "git_commit_hash": git_commit_hash or "unrecorded",
        "active_channels": list(run.active_channels),
        "final_accumulated_max_nNs_per_um2": _max("accumulated_traction_nNs_per_um2"),
        "final_accumulated_mean_nNs_per_um2": _mean("accumulated_traction_nNs_per_um2"),
        "final_stiffness_kpa_max": _max("stiffness_kpa"),
        "final_stiffness_kpa_min": _min("stiffness_kpa"),
        "final_ligand_density_max": _max("ligand_density"),
        "final_ligand_density_min": _min("ligand_density"),
        "final_fiber_density_max": _max("fiber_density"),
        "final_fiber_density_min": _min("fiber_density"),
        "final_orientation_tensor_abs_max": _abs_max("orientation_tensor"),
        "failure": run.failure,
    }


def _write_metadata_json(run: EcmOlRun, payload: dict) -> str:
    metadata_path = os.path.join(run.output_dir, "metadata.json")
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
    return metadata_path


def _make_summary_html(run: EcmOlRun, *, payload: dict) -> str:
    rows = "".join(
        f'<tr><td>{d.step_index}</td><td>{d.time_s:.4e}</td>'
        f'<td>{d.accumulated_max_nNs_per_um2:.4e}</td>'
        f'<td>{d.accumulated_mean_nNs_per_um2:.4e}</td>'
        f'<td>{d.accumulated_sum_nNs_per_um2:.4e}</td>'
        f'<td>{d.nonzero_fraction:.4e}</td>'
        f'<td>{d.traction_max_nN_per_um2:.4e}</td>'
        f'<td>{d.stiffness_kpa_max:.4e}</td>'
        f'<td>{d.ligand_density_max:.4e}</td>'
        f'<td>{d.fiber_density_max:.4e}</td>'
        f'<td>{d.orientation_tensor_abs_max:.4e}</td></tr>'
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
    failure_block = ""
    if run.failure is not None:
        failure_block = (
            f'<h2>Failure</h2><pre>{json.dumps(run.failure, indent=2)}</pre>'
        )
    status_table_rows = []
    for key in (
        "scenario_name",
        "status",
        "n_steps",
        "final_step_index",
        "frame_interval",
        "dt_s",
        "wall_clock_s",
        "git_commit_hash",
        "active_channels",
        "final_accumulated_max_nNs_per_um2",
        "final_accumulated_mean_nNs_per_um2",
        "final_stiffness_kpa_max",
        "final_stiffness_kpa_min",
        "final_ligand_density_max",
        "final_ligand_density_min",
        "final_fiber_density_max",
        "final_fiber_density_min",
        "final_orientation_tensor_abs_max",
    ):
        status_table_rows.append(
            f"<tr><th>{key}</th><td>{payload.get(key)}</td></tr>"
        )
    status_table = "".join(status_table_rows)
    body = (
        f"<!doctype html><html><head><meta charset=\"utf-8\">"
        f"<title>{run.scenario_name} ECM-OL summary</title>"
        f"<style>body{{font-family:sans-serif;margin:1em;}}"
        f"table{{border-collapse:collapse;margin-top:1em;}}"
        f"td,th{{border:1px solid #999;padding:0.25em 0.5em;text-align:left;}}"
        f"img{{margin:0.25em;}}"
        f"</style></head><body>"
        f"<h1>{run.scenario_name} ECM-OL summary</h1>"
        f"<h2>Status</h2><table><tbody>{status_table}</tbody></table>"
        f"{diag_block}"
        f"{failure_block}"
        f"<h2>Frame thumbnails</h2><div>{thumbs}</div>"
        f"<h2>Per-step metrics</h2>"
        f"<table><thead><tr><th>step</th><th>time_s</th><th>max accum</th>"
        f"<th>mean accum</th><th>sum accum</th><th>nonzero frac</th>"
        f"<th>traction max</th><th>stiffness max</th><th>ligand max</th>"
        f"<th>fiber max</th><th>orient |max|</th></tr></thead>"
        f"<tbody>{rows}</tbody></table></body></html>"
    )
    summary_path = os.path.join(run.output_dir, "summary.html")
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(body)
    return summary_path


def _write_failure_report(run: EcmOlRun) -> str:
    last_state = run.history[-1] if run.history else None
    payload = {
        "scenario_name": run.scenario_name,
        "final_step_index": run.final_step_index,
        "active_channels": list(run.active_channels),
        "failure": run.failure,
        "last_metrics": last_state.__dict__ if last_state is not None else None,
        "frame_paths": [os.path.basename(p) for p in run.frame_paths],
        "diagnostic_plot": (
            os.path.basename(run.diagnostic_plot_path)
            if run.diagnostic_plot_path
            else None
        ),
    }
    body = "# ECM-OL preflight harness — failure report\n\n```json\n"
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


def _apply_step_channels(
    ecm: ECMSubstrateState,
    scenario: EcmOlScenario,
    step_index: int,
) -> tuple[ECMSubstrateState, Optional[np.ndarray]]:
    """Apply the four open-loop preflight channels for one step in
    fixed order: traction → stiffness → density → orientation. Each
    channel is skipped if its factory is ``None``. Returns the
    evolved ``ECMSubstrateState`` and the prescribed traction array
    (or ``None``) for diagnostic recording."""

    traction: Optional[np.ndarray] = None
    if scenario.traction_factory is not None:
        traction = scenario.traction_factory(step_index)
        ecm = accumulate_prescribed_traction(ecm, traction, dt_s=scenario.dt_s)
    if scenario.stiffness_rate_factory is not None:
        stiffness_rate = scenario.stiffness_rate_factory(step_index)
        ecm = apply_prescribed_stiffness_rate(
            ecm, stiffness_rate, dt_s=scenario.dt_s
        )
    if (
        scenario.ligand_density_rate_factory is not None
        or scenario.fiber_density_rate_factory is not None
    ):
        ligand_rate = (
            scenario.ligand_density_rate_factory(step_index)
            if scenario.ligand_density_rate_factory is not None
            else np.zeros(ecm.grid_shape, dtype=np.float64)
        )
        fiber_rate = (
            scenario.fiber_density_rate_factory(step_index)
            if scenario.fiber_density_rate_factory is not None
            else np.zeros(ecm.grid_shape, dtype=np.float64)
        )
        ecm = apply_prescribed_density_rate(
            ecm, ligand_rate, fiber_rate, dt_s=scenario.dt_s
        )
    if scenario.orientation_rate_factory is not None:
        orientation_rate = scenario.orientation_rate_factory(step_index)
        ecm = apply_prescribed_orientation_rate(
            ecm, orientation_rate, dt_s=scenario.dt_s
        )
    return ecm, traction


def run_ecm_ol_scenario(
    scenario: EcmOlScenario,
    output_dir: str,
    *,
    frame_interval: int = 100,
    git_commit_hash: str = "",
) -> EcmOlRun:
    """Run one ECM-OL scenario and produce frame/PNG/HTML/summary/metadata
    artifacts. Per step the harness applies whichever channels the
    scenario provides, in the fixed order
    ``traction → stiffness_rate → density_rate → orientation_rate``.
    On failure (e.g. a factory returning a malformed array, or any
    underlying preflight raising ``ECMOpenLoopError``), the harness
    records the failure with ``failure_kind`` and writes
    ``failure_report.md`` — same pattern as
    :mod:`acs.v2.active_contour_harness`.
    """

    _require_non_negative_int(scenario.n_steps, "scenario.n_steps")
    _require_positive_int(frame_interval, "frame_interval")
    if not np.isfinite(scenario.dt_s) or scenario.dt_s < 0.0:
        raise ValueError(f"scenario.dt_s must be finite non-negative, got {scenario.dt_s!r}")
    os.makedirs(output_dir, exist_ok=True)
    run = EcmOlRun(
        scenario_name=scenario.name,
        output_dir=output_dir,
        n_steps=scenario.n_steps,
        dt_s=float(scenario.dt_s),
        frame_interval=frame_interval,
        active_channels=scenario.active_channels(),
    )
    wall_clock_start = time.perf_counter()
    ecm = scenario.initial_ecm
    run.history.append(_record_diagnostics(ecm, None, 0, scenario.dt_s))
    try:
        _write_frame_artifacts(ecm, run, step_index=0, time_s=0.0)
    except Exception as exc:  # pragma: no cover
        run.failure = {
            "kind": "frame_write_failure",
            "step_index": 0,
            "exception_type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        run.failure_report_path = _write_failure_report(run)
        return run

    for i in range(1, scenario.n_steps + 1):
        try:
            ecm, traction = _apply_step_channels(ecm, scenario, i)
        except ECMOpenLoopError as exc:
            run.failure = {
                "kind": "ecm_open_loop_violation",
                "step_index": i,
                "failure_kind": exc.failure_kind,
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
            break
        except Exception as exc:  # pragma: no cover
            run.failure = {
                "kind": "scenario_factory_failure",
                "step_index": i,
                "exception_type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            }
            break
        run.history.append(_record_diagnostics(ecm, traction, i, scenario.dt_s))
        if i % frame_interval == 0 or i == scenario.n_steps:
            try:
                _write_frame_artifacts(
                    ecm, run, step_index=i, time_s=float(i) * scenario.dt_s
                )
            except Exception as exc:  # pragma: no cover
                run.failure = {
                    "kind": "frame_write_failure",
                    "step_index": i,
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                    "traceback": traceback.format_exc(),
                }
                break

    run.final_ecm = ecm
    run.final_step_index = run.history[-1].step_index if run.history else 0
    run.wall_clock_s = float(time.perf_counter() - wall_clock_start)
    run.diagnostic_plot_path = _make_diagnostic_plot(run)
    status = "FAIL" if run.failure is not None else "PASS"
    payload = _build_status_payload(run, status, git_commit_hash)
    run.metadata_path = _write_metadata_json(run, payload)
    run.summary_path = _make_summary_html(run, payload=payload)
    if run.failure is not None:
        run.failure_report_path = _write_failure_report(run)
    return run


def make_default_ecm(
    *, nx: int = 8, ny: int = 8, spacing_um: float = 1.0, source: str = "synthetic"
) -> ECMSubstrateState:
    """Convenience helper for tests / scripts: a uniform-zero ECM grid."""

    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    return ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=spacing_um,
        stiffness_kpa=np.zeros((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source=source,  # type: ignore[arg-type]
    )
