"""Run the open-loop ECM sensitivity sweep harness (Phase C of the
closed-loop ECM gate phased plan) and persist artifacts to a
UTC-timestamped directory under ``runs/``.

Per ``docs/v2_ecm_ol_sweep_sanity_gate.md`` and the locked phased
plan (``docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`` §1
Phase C), this is **open-loop sweep baseline** evidence — it
varies ``(grid_n, spacing_um, dt_s)`` over the four ECM-OL
preflight functions and records per-tuple per-channel summary
metrics. It does **NOT** satisfy closed-loop ECM gate Item 5
("closed-loop grid spacing / dt sensitivity"), which requires
scattering geometry + response law (Phase D / E).

Wording lock (Hard Rule 11 protection):
- Every output (summary.html, metadata.json, index.json, console)
  uses "open-loop sweep baseline".
- No output may say "Item 5 satisfied" or "closed-loop sensitivity".

Caller-supplied contract (Codex impl review id=1248 strict mode):
- The harness function ``run_ecm_ol_sensitivity_sweep`` takes
  every input as a mandatory argument; no ``=...`` defaults in
  the signature for the locked inputs (tuple list, per-channel
  scenarios, reduction choice, ``max_grid_cells_total``).
- The driver script's ``_build_scenarios`` helper carries an
  **explicit non-production fixture** labeled
  ``fixture_kind="non_production_smoke"`` in metadata.json +
  index.json. CLI may override every fixture entry.
- Memory cap is mandatory: ``max_grid_cells_total`` rejects any
  tuple whose ``grid_n[0] * grid_n[1]`` exceeds the cap, BEFORE
  any ECM is constructed.

Per-tuple per-channel isolation (gate §3):
- Each (tuple, channel) builds a fresh ``ECMSubstrateState`` via
  ``make_default_ecm``. No shared mutable state.

f8cdff3 hygiene:
- Seconds-resolution timestamp ``%Y%m%dT%H%M%SZ``.
- ``os.makedirs(exist_ok=False)`` so a colliding run-root fails
  loudly (FileExistsError) instead of silently merging artifacts.
- Aggregate status computed against ``expected_status`` per
  scenario; ``index.json`` exposes
  ``intentional_failures_observed``, ``unexpected_failures``,
  ``unexpected_passes`` separately.

Usage (from repo root, .venv-collab Python):

    /Users/sw1/ActiveCellSim/.venv-collab/bin/python3 \
        scripts/run_ecm_ol_sensitivity_sweep.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import traceback
from dataclasses import dataclass, field
from typing import Callable, Mapping, Optional

# Ensure the repo root is importable before pulling in acs.v2 modules.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import matplotlib  # noqa: E402

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from acs.v2.dynamics.ecm_open_loop import (  # noqa: E402
    ECMOpenLoopError,
    accumulate_prescribed_traction,
    apply_prescribed_density_rate,
    apply_prescribed_orientation_rate,
    apply_prescribed_stiffness_rate,
)
from acs.v2.ecm_open_loop_harness import make_default_ecm  # noqa: E402
from acs.v2.ecm_substrate import ECMSubstrateState  # noqa: E402

_ALLOWED_REDUCTIONS = ("max", "mean", "sum")
_CHANNELS = (
    "traction",
    "stiffness_rate",
    "density_rate",
    "orientation_rate",
)


class SweepValidationError(Exception):
    """Caller-supplied sweep input failed validation. Carries a
    machine-readable ``failure_kind`` for the index.json payload."""

    def __init__(self, failure_kind: str, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


@dataclass(frozen=True)
class SweepTuple:
    """One sweep point: caller-supplied grid + dt + step parameters."""

    grid_n: tuple[int, int]
    spacing_um: float
    dt_s: float
    n_steps: int
    frame_interval: int
    expected_status: str  # "PASS" | "FAIL"
    label: str  # short human-readable identifier, used in directory names


@dataclass(frozen=True)
class ChannelScenario:
    """Per-channel caller-supplied prescribed-input factory.

    The factory returns the prescribed input for the given
    ``(grid_shape, step_index)``. Shape contract per channel:

    - traction: ``(grid_shape) -> (nx, ny) float64``, non-negative
    - stiffness_rate: ``(grid_shape) -> (nx, ny) float64``, signed
    - density_rate: ``(grid_shape) -> tuple[(nx,ny), (nx,ny)]``
      (ligand, fiber) signed
    - orientation_rate: ``(grid_shape) -> (nx, ny, 2, 2) float64``,
      symmetric
    """

    name: str  # one of _CHANNELS
    factory: Callable


@dataclass
class TupleChannelRecord:
    tuple_label: str
    channel: str
    status: str = "PASS"
    expected_status: str = "PASS"
    failure: Optional[dict] = None
    final_field_max: float = 0.0
    final_field_min: float = 0.0
    final_field_mean: float = 0.0
    final_field_sum: float = 0.0
    selected_reduction_value: float = 0.0
    reduction_used: str = ""
    n_steps_completed: int = 0
    per_step_max: list[float] = field(default_factory=list)
    per_step_min: list[float] = field(default_factory=list)
    per_step_mean: list[float] = field(default_factory=list)


@dataclass
class SweepRun:
    output_dir: str
    fixture_kind: str
    reduction_choice: str
    max_grid_cells_total: int
    records: list[TupleChannelRecord] = field(default_factory=list)
    aggregate_status: str = "PASS"
    intentional_failures_observed: list[str] = field(default_factory=list)
    unexpected_failures: list[str] = field(default_factory=list)
    unexpected_passes: list[str] = field(default_factory=list)


def _git_commit_hash(repo: str) -> str:
    try:
        out = subprocess.check_output(
            ["git", "-C", repo, "rev-parse", "HEAD"], text=True
        )
        return out.strip()
    except subprocess.CalledProcessError:
        return ""


def _validate_int(value, name: str, *, positive: bool = True) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise SweepValidationError(
            f"sweep_{name}_invalid",
            f"{name} must be int, got {type(value).__name__}",
        )
    if positive and value <= 0:
        raise SweepValidationError(
            f"sweep_{name}_invalid",
            f"{name} must be positive, got {value!r}",
        )
    if not positive and value < 0:
        raise SweepValidationError(
            f"sweep_{name}_invalid",
            f"{name} must be non-negative, got {value!r}",
        )
    return value


def _validate_float(value, name: str, *, positive: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SweepValidationError(
            f"sweep_{name}_invalid",
            f"{name} must be a finite scalar, got {type(value).__name__}",
        )
    v = float(value)
    if not np.isfinite(v):
        raise SweepValidationError(
            f"sweep_{name}_invalid",
            f"{name} must be finite, got {value!r}",
        )
    if positive and v <= 0.0:
        raise SweepValidationError(
            f"sweep_{name}_invalid",
            f"{name} must be > 0, got {value!r}",
        )
    if not positive and v < 0.0:
        raise SweepValidationError(
            f"sweep_{name}_invalid",
            f"{name} must be >= 0, got {value!r}",
        )
    return v


def _validate_sweep_tuple(t: SweepTuple, max_grid_cells_total: int) -> None:
    nx, ny = t.grid_n
    _validate_int(nx, "grid_n", positive=True)
    _validate_int(ny, "grid_n", positive=True)
    _validate_float(t.spacing_um, "spacing", positive=True)
    _validate_float(t.dt_s, "dt", positive=False)
    _validate_int(t.n_steps, "n_steps", positive=True)
    _validate_int(t.frame_interval, "n_steps", positive=True)
    if t.expected_status not in {"PASS", "FAIL"}:
        raise SweepValidationError(
            "sweep_n_steps_invalid",
            f"expected_status must be PASS or FAIL, got {t.expected_status!r}",
        )
    if nx * ny > max_grid_cells_total:
        raise SweepValidationError(
            "sweep_memory_cap_exceeded",
            f"tuple {t.label!r} grid {nx}×{ny} = {nx*ny} cells exceeds "
            f"max_grid_cells_total = {max_grid_cells_total}",
        )


def _reduce(field: np.ndarray, reduction: str) -> float:
    if field.size == 0:
        return 0.0
    if reduction == "max":
        return float(field.max())
    if reduction == "mean":
        return float(field.mean())
    if reduction == "sum":
        return float(field.sum())
    raise SweepValidationError(
        "sweep_n_steps_invalid",
        f"unknown reduction_choice {reduction!r}; allowed: {_ALLOWED_REDUCTIONS}",
    )


def _post_state_for_record(
    final_ecm: ECMSubstrateState, channel: str
) -> np.ndarray:
    if channel == "traction":
        return np.asarray(
            final_ecm.accumulated_traction_nNs_per_um2, dtype=np.float64
        )
    if channel == "stiffness_rate":
        return np.asarray(final_ecm.stiffness_kpa, dtype=np.float64)
    if channel == "density_rate":
        # Combined: report ligand. (Detailed per-field metrics in metadata.)
        return np.asarray(final_ecm.ligand_density, dtype=np.float64)
    if channel == "orientation_rate":
        # |T_ij| max -> use abs
        return np.abs(np.asarray(final_ecm.orientation_tensor, dtype=np.float64))
    raise SweepValidationError(
        "sweep_n_steps_invalid",
        f"unknown channel {channel!r}; allowed: {_CHANNELS}",
    )


def _record_per_step(rec: TupleChannelRecord, ecm: ECMSubstrateState) -> None:
    """Record the post-step max/min/mean of the channel's field
    on this ECM state, so the per-(tuple, channel) diagnostic plot
    can show the trajectory."""

    field = _post_state_for_record(ecm, rec.channel)
    if field.size == 0:
        rec.per_step_max.append(0.0)
        rec.per_step_min.append(0.0)
        rec.per_step_mean.append(0.0)
        return
    rec.per_step_max.append(float(field.max()))
    rec.per_step_min.append(float(field.min()))
    rec.per_step_mean.append(float(field.mean()))


def _run_one_tuple_channel(
    t: SweepTuple,
    scenario: ChannelScenario,
    reduction: str,
) -> TupleChannelRecord:
    rec = TupleChannelRecord(
        tuple_label=t.label,
        channel=scenario.name,
        expected_status=t.expected_status,
        reduction_used=reduction,
    )
    nx, ny = t.grid_n
    grid_shape = (nx, ny)
    # Per-(tuple, channel) ECM isolation: fresh state per call.
    ecm = make_default_ecm(nx=nx, ny=ny, spacing_um=t.spacing_um)
    # Step 0 baseline (initial state) for the per-step series.
    _record_per_step(rec, ecm)
    try:
        for i in range(1, t.n_steps + 1):
            if scenario.name == "traction":
                traction = scenario.factory(grid_shape, i)
                ecm = accumulate_prescribed_traction(ecm, traction, dt_s=t.dt_s)
            elif scenario.name == "stiffness_rate":
                rate = scenario.factory(grid_shape, i)
                ecm = apply_prescribed_stiffness_rate(ecm, rate, dt_s=t.dt_s)
            elif scenario.name == "density_rate":
                ligand_rate, fiber_rate = scenario.factory(grid_shape, i)
                ecm = apply_prescribed_density_rate(
                    ecm, ligand_rate, fiber_rate, dt_s=t.dt_s
                )
            elif scenario.name == "orientation_rate":
                rate = scenario.factory(grid_shape, i)
                ecm = apply_prescribed_orientation_rate(ecm, rate, dt_s=t.dt_s)
            else:
                raise SweepValidationError(
                    "sweep_n_steps_invalid",
                    f"unknown channel {scenario.name!r}",
                )
            rec.n_steps_completed = i
            _record_per_step(rec, ecm)
    except ECMOpenLoopError as exc:
        rec.status = "FAIL"
        rec.failure = {
            "kind": "ecm_open_loop_violation",
            "step_index": rec.n_steps_completed + 1,
            "failure_kind": exc.failure_kind,
            "exception_type": type(exc).__name__,
            "message": str(exc),
        }
        return rec
    except SweepValidationError as exc:
        rec.status = "FAIL"
        rec.failure = {
            "kind": "sweep_validation_error",
            "step_index": rec.n_steps_completed + 1,
            "failure_kind": exc.failure_kind,
            "exception_type": type(exc).__name__,
            "message": str(exc),
        }
        return rec
    except Exception as exc:  # pragma: no cover
        rec.status = "FAIL"
        rec.failure = {
            "kind": "scenario_factory_failure",
            "step_index": rec.n_steps_completed + 1,
            "exception_type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        return rec

    field_post = _post_state_for_record(ecm, scenario.name)
    rec.final_field_max = float(field_post.max()) if field_post.size else 0.0
    rec.final_field_min = float(field_post.min()) if field_post.size else 0.0
    rec.final_field_mean = float(field_post.mean()) if field_post.size else 0.0
    rec.final_field_sum = float(field_post.sum()) if field_post.size else 0.0
    rec.selected_reduction_value = _reduce(field_post, reduction)
    return rec


def run_ecm_ol_sensitivity_sweep(
    sweep_tuples,
    channel_scenarios,
    reduction_choice,
    max_grid_cells_total,
    output_dir,
    *,
    git_commit_hash,
    fixture_kind,
):
    """Run the open-loop ECM sensitivity sweep over caller-supplied
    tuples and per-channel scenarios. Every input is mandatory — no
    ``=...`` defaults for the locked inputs (per Codex impl review
    id=1248).

    Open-loop sweep baseline only. Does NOT satisfy closed-loop ECM
    gate Item 5.

    Returns a :class:`SweepRun` with per-(tuple, channel) records,
    aggregate status against ``expected_status`` per tuple, and the
    paths of the per-tuple summary HTML / metadata JSON / diagnostic
    plot. The cross-tuple summary plot lives at
    ``output_dir/sweep_summary.png``.
    """

    # Validate reduction choice up front so an invalid value fails
    # before any per-tuple work runs.
    if reduction_choice not in _ALLOWED_REDUCTIONS:
        raise SweepValidationError(
            "sweep_n_steps_invalid",
            f"reduction_choice must be one of {_ALLOWED_REDUCTIONS}, "
            f"got {reduction_choice!r}",
        )
    _validate_int(max_grid_cells_total, "memory_cap", positive=True)

    # Validate all sweep tuples up front so an invalid tuple cannot
    # leave the run-root half-populated.
    for t in sweep_tuples:
        _validate_sweep_tuple(t, max_grid_cells_total)

    # Validate channel scenarios: must be EXACTLY one ChannelScenario per
    # known channel in _CHANNELS, no missing, no duplicates. Per Codex
    # impl review id=1253, accepting subsets or duplicates would let a
    # caller publish a "PASS" bundle that does not exercise all four
    # ECM-OL preflight functions per Phase C scope.
    seen: dict[str, int] = {}
    for scen in channel_scenarios:
        if scen.name not in _CHANNELS:
            raise SweepValidationError(
                "sweep_channel_invalid",
                f"channel name {scen.name!r} not in {_CHANNELS}",
            )
        seen[scen.name] = seen.get(scen.name, 0) + 1
    for ch in _CHANNELS:
        if seen.get(ch, 0) == 0:
            raise SweepValidationError(
                "sweep_channel_invalid",
                f"channel scenarios missing required channel {ch!r}; "
                f"Phase C requires exactly one ChannelScenario per "
                f"channel in {list(_CHANNELS)}",
            )
    duplicates = [ch for ch, count in seen.items() if count > 1]
    if duplicates:
        raise SweepValidationError(
            "sweep_channel_invalid",
            f"channel scenarios contain duplicate channels {duplicates}; "
            f"each channel must have exactly one ChannelScenario",
        )

    os.makedirs(output_dir, exist_ok=False)
    run = SweepRun(
        output_dir=output_dir,
        fixture_kind=fixture_kind,
        reduction_choice=reduction_choice,
        max_grid_cells_total=max_grid_cells_total,
    )

    for t in sweep_tuples:
        for scen in channel_scenarios:
            rec = _run_one_tuple_channel(t, scen, reduction_choice)
            run.records.append(rec)

    # Aggregate against expected_status per tuple.
    # A tuple is treated as "PASS" overall if every channel's status
    # equals the tuple's expected_status; otherwise the tuple counts
    # as a mismatch.
    tuple_actual_status: dict[str, str] = {}
    tuple_expected: dict[str, str] = {t.label: t.expected_status for t in sweep_tuples}
    for rec in run.records:
        prev = tuple_actual_status.get(rec.tuple_label)
        if prev is None:
            tuple_actual_status[rec.tuple_label] = rec.status
        elif rec.status == "FAIL":
            tuple_actual_status[rec.tuple_label] = "FAIL"

    for label, actual in tuple_actual_status.items():
        expected = tuple_expected[label]
        if actual == "FAIL" and expected == "FAIL":
            run.intentional_failures_observed.append(label)
        elif actual == "FAIL" and expected != "FAIL":
            run.unexpected_failures.append(label)
        elif actual != "FAIL" and expected == "FAIL":
            run.unexpected_passes.append(label)

    if run.unexpected_failures or run.unexpected_passes:
        run.aggregate_status = "FAIL"

    # Per-tuple summary.html + metadata.json + diagnostic per channel.
    for t in sweep_tuples:
        tuple_dir = os.path.join(output_dir, t.label)
        os.makedirs(tuple_dir, exist_ok=False)
        tuple_records = [r for r in run.records if r.tuple_label == t.label]
        per_channel = {r.channel: r for r in tuple_records}
        payload = {
            "tuple_label": t.label,
            "fixture_kind": fixture_kind,
            "reduction_choice": reduction_choice,
            "max_grid_cells_total": max_grid_cells_total,
            "grid_n": list(t.grid_n),
            "spacing_um": t.spacing_um,
            "dt_s": t.dt_s,
            "n_steps": t.n_steps,
            "frame_interval": t.frame_interval,
            "expected_status": t.expected_status,
            "actual_tuple_status": tuple_actual_status[t.label],
            "git_commit_hash": git_commit_hash or "unrecorded",
            "evidence_kind": "open-loop sweep baseline",
            "channels": {
                name: {
                    "status": rec.status,
                    "expected_status": rec.expected_status,
                    "failure": rec.failure,
                    "final_field_max": rec.final_field_max,
                    "final_field_min": rec.final_field_min,
                    "final_field_mean": rec.final_field_mean,
                    "final_field_sum": rec.final_field_sum,
                    "selected_reduction_value": rec.selected_reduction_value,
                    "reduction_used": rec.reduction_used,
                    "n_steps_completed": rec.n_steps_completed,
                    "per_step_max": rec.per_step_max,
                    "per_step_min": rec.per_step_min,
                    "per_step_mean": rec.per_step_mean,
                }
                for name, rec in per_channel.items()
            },
        }
        meta_path = os.path.join(tuple_dir, "metadata.json")
        with open(meta_path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, default=str)
        summary_path = os.path.join(tuple_dir, "summary.html")
        rows = "".join(
            f"<tr><th>{name}</th><td>{rec.status}</td>"
            f"<td>{rec.expected_status}</td>"
            f"<td>{rec.final_field_max:.4e}</td>"
            f"<td>{rec.final_field_min:.4e}</td>"
            f"<td>{rec.final_field_mean:.4e}</td>"
            f"<td>{rec.reduction_used}</td>"
            f"<td>{rec.n_steps_completed}</td>"
            f"<td>{rec.failure}</td></tr>"
            for name, rec in per_channel.items()
        )
        body = (
            f"<!doctype html><html><head><meta charset=\"utf-8\">"
            f"<title>{t.label} open-loop sweep baseline</title>"
            f"<style>body{{font-family:sans-serif;margin:1em;}}"
            f"table{{border-collapse:collapse;margin-top:1em;}}"
            f"td,th{{border:1px solid #999;padding:0.25em 0.5em;text-align:left;}}"
            f"</style></head><body>"
            f"<h1>{t.label}: open-loop sweep baseline</h1>"
            f"<p>Evidence kind: <strong>open-loop sweep baseline</strong>. "
            f"This artifact does NOT claim closed-loop ECM gate Item 5 "
            f"satisfaction.</p>"
            f"<p>fixture_kind: {fixture_kind}; expected_status: "
            f"{t.expected_status}; actual_tuple_status: "
            f"{tuple_actual_status[t.label]}</p>"
            f"<table><thead><tr><th>channel</th><th>status</th>"
            f"<th>expected</th><th>final_max</th><th>final_min</th>"
            f"<th>final_mean</th><th>reduction</th>"
            f"<th>steps_completed</th><th>failure</th></tr></thead>"
            f"<tbody>{rows}</tbody></table></body></html>"
        )
        with open(summary_path, "w", encoding="utf-8") as fh:
            fh.write(body)

    # Cross-tuple summary plot per channel: selected_reduction_value
    # (from _reduce(field, reduction)) vs tuple ordering, one trace
    # per channel.
    if run.records:
        fig, axes = plt.subplots(
            len(_CHANNELS), 1, figsize=(7.0, 1.8 * len(_CHANNELS)), sharex=True
        )
        if len(_CHANNELS) == 1:
            axes = [axes]
        labels = [t.label for t in sweep_tuples]
        for ax, channel in zip(axes, _CHANNELS):
            xs = list(range(len(labels)))
            ys = []
            for label in labels:
                matching = [
                    r
                    for r in run.records
                    if r.tuple_label == label and r.channel == channel
                ]
                if not matching or matching[0].status == "FAIL":
                    ys.append(np.nan)
                    continue
                ys.append(matching[0].selected_reduction_value)
            ax.plot(xs, ys, marker="o")
            ax.set_ylabel(f"{channel}\n({reduction_choice})", fontsize=8)
        axes[-1].set_xticks(list(range(len(labels))))
        axes[-1].set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        axes[-1].set_xlabel("tuple")
        fig.suptitle("open-loop sweep baseline (NOT Item 5 satisfaction)")
        fig.tight_layout()
        plot_path = os.path.join(output_dir, "sweep_summary.png")
        fig.savefig(plot_path, dpi=120)
        plt.close(fig)

    # Per-(tuple, channel) diagnostic plot: per-step max/min/mean
    # series of the post-step field over the n_steps run, so a
    # downstream reader can see the trajectory at each (grid, dt)
    # point. Sanity Gate §0/§8 contract.
    for t in sweep_tuples:
        tuple_dir = os.path.join(output_dir, t.label)
        for scen in channel_scenarios:
            rec = next(
                r
                for r in run.records
                if r.tuple_label == t.label and r.channel == scen.name
            )
            if not rec.per_step_max:
                continue
            fig, ax = plt.subplots(figsize=(6.0, 3.0))
            steps = list(range(len(rec.per_step_max)))
            ax.plot(steps, rec.per_step_max, label="max", marker="o")
            ax.plot(steps, rec.per_step_min, label="min", marker="s")
            ax.plot(steps, rec.per_step_mean, label="mean", marker="^")
            ax.set_xlabel("step index")
            ax.set_ylabel(f"{scen.name} field reduction")
            ax.set_title(
                f"{t.label} {scen.name} per-step (open-loop sweep baseline)"
            )
            ax.legend(loc="best", fontsize=8)
            fig.tight_layout()
            chan_plot_path = os.path.join(
                tuple_dir, f"diagnostic_{scen.name}.png"
            )
            fig.savefig(chan_plot_path, dpi=120)
            plt.close(fig)

    # index.json
    index_path = os.path.join(output_dir, "index.json")
    with open(index_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "git_commit_hash": git_commit_hash or "unrecorded",
                "fixture_kind": fixture_kind,
                "evidence_kind": "open-loop sweep baseline",
                "reduction_choice": reduction_choice,
                "max_grid_cells_total": max_grid_cells_total,
                "tuples": [
                    {
                        "label": t.label,
                        "grid_n": list(t.grid_n),
                        "spacing_um": t.spacing_um,
                        "dt_s": t.dt_s,
                        "n_steps": t.n_steps,
                        "expected_status": t.expected_status,
                        "actual_status": tuple_actual_status[t.label],
                    }
                    for t in sweep_tuples
                ],
                "aggregate_status": run.aggregate_status,
                "intentional_failures_observed": run.intentional_failures_observed,
                "unexpected_failures": run.unexpected_failures,
                "unexpected_passes": run.unexpected_passes,
                "run_root": output_dir,
            },
            fh,
            indent=2,
            default=str,
        )
    return run


def _build_scenarios(grid_shape: tuple[int, int]):
    """Explicit non-production fixture for the script-driven smoke run.

    Per Codex impl review id=1248, this fixture is labeled
    `non_production_smoke` in metadata.json + index.json. The CLI
    can override every value. The harness function does NOT take
    these values from a default; the script passes them explicitly.
    """

    def _traction_uniform(grid_shape, step_index):  # noqa: ARG001
        return np.full(grid_shape, 0.3, dtype=np.float64)

    def _stiffness_rate_uniform(grid_shape, step_index):  # noqa: ARG001
        return np.full(grid_shape, 0.2, dtype=np.float64)

    def _density_rate_uniform(grid_shape, step_index):  # noqa: ARG001
        ligand = np.full(grid_shape, 0.05, dtype=np.float64)
        fiber = np.full(grid_shape, 0.02, dtype=np.float64)
        return ligand, fiber

    def _orientation_rate_diagonal(grid_shape, step_index):  # noqa: ARG001
        rate = np.zeros((*grid_shape, 2, 2), dtype=np.float64)
        rate[..., 0, 0] = -0.05
        rate[..., 1, 1] = -0.05
        return rate

    return [
        ChannelScenario(name="traction", factory=_traction_uniform),
        ChannelScenario(name="stiffness_rate", factory=_stiffness_rate_uniform),
        ChannelScenario(name="density_rate", factory=_density_rate_uniform),
        ChannelScenario(name="orientation_rate", factory=_orientation_rate_diagonal),
    ]


def _build_tuples_smoke_fixture():
    """Explicit non-production smoke fixture: 3 grid sizes × 3 dt
    values = 9 tuples. Caller must explicitly pass this fixture to
    the harness; it is NOT a runtime default."""

    grid_sizes = [(2, 2), (4, 4), (6, 6)]
    dt_values = [0.1, 0.05, 0.01]
    tuples = []
    for nx, ny in grid_sizes:
        for dt_s in dt_values:
            tuples.append(
                SweepTuple(
                    grid_n=(nx, ny),
                    spacing_um=1.0,
                    dt_s=dt_s,
                    n_steps=8,
                    frame_interval=4,
                    expected_status="PASS",
                    label=f"g{nx}x{ny}_dt{dt_s:g}",
                )
            )
    return tuples


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        default=_REPO_ROOT,
        help="Repository root used for git commit lookup and runs/ output.",
    )
    parser.add_argument(
        "--max-grid-cells-total",
        type=int,
        required=True,
        help=(
            "Mandatory memory cap (cells per grid). The harness function has "
            "no default; the CLI surfaces this requirement. The Sanity Gate "
            "doc suggests 4096 for Phase 1 baseline as a worked example."
        ),
    )
    parser.add_argument(
        "--reduction",
        choices=_ALLOWED_REDUCTIONS,
        required=True,
        help="Cross-tuple reduction choice (mandatory; no default).",
    )
    args = parser.parse_args()
    repo = os.path.abspath(args.repo)

    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = os.path.join(repo, "runs", f"{timestamp}_ecm_ol_sensitivity")
    git_hash = _git_commit_hash(repo)

    sweep_tuples = _build_tuples_smoke_fixture()
    grid_shape_for_scenarios = sweep_tuples[0].grid_n if sweep_tuples else (2, 2)
    channel_scenarios = _build_scenarios(grid_shape_for_scenarios)

    run = run_ecm_ol_sensitivity_sweep(
        sweep_tuples=sweep_tuples,
        channel_scenarios=channel_scenarios,
        reduction_choice=args.reduction,
        max_grid_cells_total=args.max_grid_cells_total,
        output_dir=run_root,
        git_commit_hash=git_hash,
        fixture_kind="non_production_smoke",
    )

    print(
        f"open-loop sweep baseline run {run.aggregate_status}. "
        f"artifacts: {run.output_dir}"
    )
    print(f"index: {os.path.join(run.output_dir, 'index.json')}")
    if run.intentional_failures_observed:
        print(
            f"intentional failure tuples (as expected): "
            f"{', '.join(run.intentional_failures_observed)}"
        )
    if run.unexpected_failures:
        print(
            f"UNEXPECTED FAILURES: {', '.join(run.unexpected_failures)}",
            file=sys.stderr,
        )
    if run.unexpected_passes:
        print(
            f"UNEXPECTED PASSES: {', '.join(run.unexpected_passes)}",
            file=sys.stderr,
        )
    return 0 if run.aggregate_status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
