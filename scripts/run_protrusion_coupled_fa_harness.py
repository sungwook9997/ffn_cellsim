"""Run the 6.3b protrusion-coupled FA dynamics over a small set of
caller-supplied scenarios and persist artifacts to a UTC-timestamped
directory under ``runs/``.

Produces, per scenario, a status table HTML
(``summary.html``), a JSON metadata payload (``metadata.json``), and
a matplotlib diagnostic plot (``diagnostic_<scenario>.png``) showing
per-FA maturity / bound-fraction trajectories plus the per-step
maximum effective rate per rate name. A top-level ``index.json``
lists every scenario summary path with an ``aggregate_status``
(PASS / FAIL) and three categorised lists —
``intentional_failures_observed``, ``unexpected_failures``, and
``unexpected_passes`` — so a downstream consumer can tell at a
glance whether any scenario broke without mistaking the intentional
demonstration for a regression.

Scenarios are caller-supplied test inputs only — no model defaults,
no biology decisions, no RNG. The script exists to make the 6.3b
wrapper visible to PI in a browser, not to claim a physical
interpretation. Reference biology multipliers from
``docs/v2/v2_63b_protrusion_coupling_locked.md`` §3 are used as test
inputs for the boost scenario; the script keeps them inline as
explicit test values, NOT as production defaults.

Usage (from repo root, .venv-collab Python):

    /Users/sw1/ActiveCellSim/.venv-collab/bin/python3 \
        scripts/run_protrusion_coupled_fa_harness.py

The script does not commit its outputs. ``runs/`` is intentionally
excluded from version control via ``.gitignore``.

Aggregate status / exit code: each scenario carries an
``expected_status`` (PASS or FAIL). The aggregate is computed
against expectation, not raw outcome — a scenario whose intentional
FAIL exercises the wrapper's failure path is **not** a regression.
Concretely:

- ``aggregate_status == "PASS"`` ⇔ every scenario's actual status
  matches its ``expected_status``. The intentional FAIL scenario
  ``dt_rate_violation_failure`` lands here when it correctly raises
  ``dt_rate_violation`` — that is the "as expected" outcome.
- ``aggregate_status == "FAIL"`` ⇔ at least one scenario's status
  does not match expectation. Reported via ``unexpected_failures``
  (expected PASS but got FAIL) and ``unexpected_passes`` (expected
  FAIL but got PASS). Exit 1 only in this case.

This refinement of the f8cdff3 ECM-OL contract preserves the
"silent unintended FAIL is impossible" guarantee while letting an
intentional demonstration of the failure path coexist with the
healthy scenarios in the same run. Per Codex review id=1216, the
intentional FAIL is also visibly tagged in ``index.json`` and in
the per-scenario ``summary.html`` status table so readers do not
mistake the whole deliverable for a regression.
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

from acs.v2.dynamics.focal_adhesion import (  # noqa: E402
    FocalAdhesionDynamicsError,
    FocalAdhesionDynamicsParameters,
)
from acs.v2.dynamics.protrusion_coupled_focal_adhesion import (  # noqa: E402
    ProtrusionStateMultipliers,
    step_protrusion_coupled_focal_adhesions,
)
from acs.v2.focal_adhesion import FocalAdhesionState  # noqa: E402
from acs.v2.protrusion import ProtrusionEvent  # noqa: E402


def _git_commit_hash(repo: str) -> str:
    try:
        out = subprocess.check_output(
            ["git", "-C", repo, "rev-parse", "HEAD"], text=True
        )
        return out.strip()
    except subprocess.CalledProcessError:
        return ""


@dataclass
class Scenario:
    name: str
    description: str
    centroid_um_xy: tuple[float, float]
    params: FocalAdhesionDynamicsParameters
    multipliers: ProtrusionStateMultipliers
    initial_adhesions: tuple[FocalAdhesionState, ...]
    protrusion_factory: Callable[[int], Mapping[str, ProtrusionEvent]]
    n_steps: int
    expected_status: str  # "PASS" | "FAIL"


@dataclass
class StepRecord:
    step_index: int
    time_s: float
    per_fa_maturity: list[float]
    per_fa_bound_fraction: list[float]
    max_effective_rate_per_name: dict
    reciprocal_missing: int


@dataclass
class ScenarioRun:
    scenario_name: str
    output_dir: str
    history: list[StepRecord] = field(default_factory=list)
    status: str = "PASS"
    expected_status: str = "PASS"
    failure: Optional[dict] = None
    final_step_index: int = 0
    n_steps: int = 0
    summary_path: Optional[str] = None
    metadata_path: Optional[str] = None
    diagnostic_plot_path: Optional[str] = None


def _record_step(
    step_index: int,
    time_s: float,
    adhesions: tuple[FocalAdhesionState, ...],
    diagnostics: dict,
) -> StepRecord:
    return StepRecord(
        step_index=step_index,
        time_s=time_s,
        per_fa_maturity=[fa.maturity for fa in adhesions],
        per_fa_bound_fraction=[fa.bound_fraction for fa in adhesions],
        max_effective_rate_per_name=dict(diagnostics["max_effective_rate_per_name"]),
        reciprocal_missing=int(diagnostics["reciprocal_missing"]),
    )


def _make_diagnostic_plot(run: ScenarioRun, fa_ids: list[str]) -> Optional[str]:
    if not run.history:
        return None
    steps = [h.step_index for h in run.history]
    n_fa = len(fa_ids)
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 9.0), sharex=True)
    for j, fa_id in enumerate(fa_ids):
        axes[0].plot(
            steps, [h.per_fa_maturity[j] for h in run.history], label=fa_id
        )
        axes[1].plot(
            steps, [h.per_fa_bound_fraction[j] for h in run.history], label=fa_id
        )
    axes[0].set_ylabel("maturity")
    axes[0].set_ylim(-0.05, 1.05)
    if n_fa > 0:
        axes[0].legend(loc="best", fontsize=8)
    axes[1].set_ylabel("bound_fraction")
    axes[1].set_ylim(-0.05, 1.05)
    if n_fa > 0:
        axes[1].legend(loc="best", fontsize=8)
    rate_names = ["k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s"]
    for name in rate_names:
        axes[2].plot(
            steps,
            [h.max_effective_rate_per_name.get(name, 0.0) for h in run.history],
            label=name,
        )
    axes[2].set_ylabel("max effective rate (1/s)")
    axes[2].set_xlabel("step index")
    axes[2].legend(loc="best", fontsize=8)
    fig.suptitle(f"{run.scenario_name} 6.3b diagnostics")
    fig.tight_layout()
    plot_path = os.path.join(run.output_dir, f"diagnostic_{run.scenario_name}.png")
    fig.savefig(plot_path, dpi=120)
    plt.close(fig)
    return plot_path


def _build_status_payload(
    run: ScenarioRun,
    scenario: Scenario,
    git_commit_hash: str,
) -> dict:
    return {
        "scenario_name": run.scenario_name,
        "description": scenario.description,
        "status": run.status,
        "expected_status": run.expected_status,
        "n_steps": run.n_steps,
        "final_step_index": run.final_step_index,
        "git_commit_hash": git_commit_hash or "unrecorded",
        "centroid_um_xy": list(scenario.centroid_um_xy),
        "params": {
            "dt_fa_s": scenario.params.dt_fa_s,
            "traction_scale_nN": scenario.params.traction_scale_nN,
            "k_maturity_per_s": scenario.params.k_maturity_per_s,
            "k_bind_per_s": scenario.params.k_bind_per_s,
            "k_unbind_per_s": scenario.params.k_unbind_per_s,
            "max_traction_nN": scenario.params.max_traction_nN,
        },
        "multipliers_table": {
            state: dict(rates)
            for state, rates in scenario.multipliers.table.items()
        },
        "fa_ids": [fa.adhesion_id for fa in scenario.initial_adhesions],
        "protrusion_state_at_step_1": {
            pid: prot.state
            for pid, prot in scenario.protrusion_factory(1).items()
        },
        "failure": run.failure,
    }


def _make_summary_html(
    run: ScenarioRun,
    scenario: Scenario,
    payload: dict,
) -> str:
    fa_ids = [fa.adhesion_id for fa in scenario.initial_adhesions]
    rows = "".join(
        f"<tr><td>{h.step_index}</td><td>{h.time_s:.4e}</td>"
        + "".join(f"<td>{m:.4e}</td>" for m in h.per_fa_maturity)
        + "".join(f"<td>{b:.4e}</td>" for b in h.per_fa_bound_fraction)
        + "".join(
            f"<td>{h.max_effective_rate_per_name.get(name, 0.0):.4e}</td>"
            for name in ("k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s")
        )
        + f"<td>{h.reciprocal_missing}</td></tr>"
        for h in run.history
    )
    fa_maturity_headers = "".join(f"<th>{fa_id} maturity</th>" for fa_id in fa_ids)
    fa_bound_headers = "".join(
        f"<th>{fa_id} bound_fraction</th>" for fa_id in fa_ids
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
        "description",
        "status",
        "expected_status",
        "n_steps",
        "final_step_index",
        "git_commit_hash",
        "centroid_um_xy",
        "params",
        "multipliers_table",
        "fa_ids",
        "protrusion_state_at_step_1",
    ):
        status_table_rows.append(
            f"<tr><th>{key}</th><td>{payload.get(key)}</td></tr>"
        )
    status_table = "".join(status_table_rows)
    body = (
        f"<!doctype html><html><head><meta charset=\"utf-8\">"
        f"<title>{run.scenario_name} 6.3b summary</title>"
        f"<style>body{{font-family:sans-serif;margin:1em;}}"
        f"table{{border-collapse:collapse;margin-top:1em;}}"
        f"td,th{{border:1px solid #999;padding:0.25em 0.5em;text-align:left;}}"
        f"img{{margin:0.25em;}}</style></head><body>"
        f"<h1>{run.scenario_name} 6.3b summary</h1>"
        f"<h2>Status</h2><table><tbody>{status_table}</tbody></table>"
        f"{diag_block}{failure_block}"
        f"<h2>Per-step trajectories</h2>"
        f"<table><thead><tr><th>step</th><th>time_s</th>"
        f"{fa_maturity_headers}{fa_bound_headers}"
        f"<th>max k_maturity_per_s</th><th>max k_bind_per_s</th>"
        f"<th>max k_unbind_per_s</th><th>reciprocal_missing</th></tr></thead>"
        f"<tbody>{rows}</tbody></table></body></html>"
    )
    summary_path = os.path.join(run.output_dir, "summary.html")
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(body)
    return summary_path


def _run_scenario(
    scenario: Scenario, output_dir: str, git_commit_hash: str
) -> ScenarioRun:
    os.makedirs(output_dir, exist_ok=False)
    run = ScenarioRun(
        scenario_name=scenario.name,
        output_dir=output_dir,
        n_steps=scenario.n_steps,
        expected_status=scenario.expected_status,
    )

    adhesions = scenario.initial_adhesions
    initial_diag = {
        "max_effective_rate_per_name": {
            "k_maturity_per_s": 0.0,
            "k_bind_per_s": 0.0,
            "k_unbind_per_s": 0.0,
        },
        "reciprocal_missing": 0,
    }
    run.history.append(_record_step(0, 0.0, adhesions, initial_diag))

    for i in range(1, scenario.n_steps + 1):
        registry = scenario.protrusion_factory(i)
        try:
            result = step_protrusion_coupled_focal_adhesions(
                adhesions,
                scenario.centroid_um_xy,
                scenario.params,
                registry,
                scenario.multipliers,
            )
        except FocalAdhesionDynamicsError as exc:
            run.status = "FAIL"
            run.failure = {
                "kind": "fa_dynamics_error",
                "step_index": i,
                "failure_kind": exc.failure_kind,
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
            break
        except Exception as exc:  # pragma: no cover
            run.status = "FAIL"
            run.failure = {
                "kind": "unexpected_exception",
                "step_index": i,
                "exception_type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            }
            break
        adhesions = result.updated_adhesions
        run.history.append(
            _record_step(
                i, float(i) * scenario.params.dt_fa_s, adhesions, result.diagnostics
            )
        )

    run.final_step_index = run.history[-1].step_index if run.history else 0
    payload = _build_status_payload(run, scenario, git_commit_hash)
    metadata_path = os.path.join(output_dir, "metadata.json")
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
    run.metadata_path = metadata_path
    fa_ids = [fa.adhesion_id for fa in scenario.initial_adhesions]
    run.diagnostic_plot_path = _make_diagnostic_plot(run, fa_ids)
    run.summary_path = _make_summary_html(run, scenario, payload)
    return run


def _build_scenarios() -> list[Scenario]:
    """Build the caller-supplied scenarios. Reference biology multipliers
    from the lock §3 doc are used here as **test inputs**, not as
    production defaults."""

    fa_centered = FocalAdhesionState(
        adhesion_id="fa-east",
        cell_id="cell-A",
        position_um_xy=(1.0, 0.0),
        age_s=0.0,
        maturity=0.4,
        bound_fraction=0.5,
        state="mature",
        linked_protrusion_id="p-1",
    )
    fa_north = FocalAdhesionState(
        adhesion_id="fa-north",
        cell_id="cell-A",
        position_um_xy=(0.0, 1.0),
        age_s=0.0,
        maturity=0.3,
        bound_fraction=0.4,
        state="mature",
        linked_protrusion_id="p-2",
    )
    centroid = (0.0, 0.0)

    base_params = FocalAdhesionDynamicsParameters(
        dt_fa_s=0.1,
        traction_scale_nN=1.0,
        k_maturity_per_s=0.5,
        k_bind_per_s=0.5,
    )

    def _registry_growing(step_index: int) -> Mapping[str, ProtrusionEvent]:  # noqa: ARG001
        return {
            "p-1": ProtrusionEvent(
                cell_id="cell-A",
                start_time_s=0.0,
                boundary_angle_rad=0.0,
                length_um=1.0,
                event_id="p-1",
                state="growing",
                associated_adhesion_ids=("fa-east",),
            )
        }

    def _registry_retracting(step_index: int) -> Mapping[str, ProtrusionEvent]:  # noqa: ARG001
        return {
            "p-1": ProtrusionEvent(
                cell_id="cell-A",
                start_time_s=0.0,
                boundary_angle_rad=0.0,
                length_um=1.0,
                event_id="p-1",
                state="retracting",
                associated_adhesion_ids=("fa-east",),
            )
        }

    def _registry_mixed(step_index: int) -> Mapping[str, ProtrusionEvent]:  # noqa: ARG001
        return {
            "p-1": ProtrusionEvent(
                cell_id="cell-A",
                start_time_s=0.0,
                boundary_angle_rad=0.0,
                length_um=1.0,
                event_id="p-1",
                state="growing",
                associated_adhesion_ids=("fa-east",),
            ),
            "p-2": ProtrusionEvent(
                cell_id="cell-A",
                start_time_s=0.0,
                boundary_angle_rad=1.5708,
                length_um=1.0,
                event_id="p-2",
                state="stalled",
                associated_adhesion_ids=("fa-north",),
            ),
        }

    neutral = ProtrusionStateMultipliers(
        table={
            "growing": {"k_maturity_per_s": 1.0, "k_bind_per_s": 1.0},
            "retracting": {"k_maturity_per_s": 1.0, "k_bind_per_s": 1.0},
        }
    )
    growing_boost = ProtrusionStateMultipliers(
        table={
            "growing": {"k_maturity_per_s": 2.0, "k_bind_per_s": 1.5},
        }
    )
    mixed_boost = ProtrusionStateMultipliers(
        table={
            "growing": {"k_maturity_per_s": 2.0, "k_bind_per_s": 1.5},
            "stalled": {"k_maturity_per_s": 1.5, "k_bind_per_s": 1.0},
        }
    )
    # dt_fa = 0.1, base k_maturity = 0.5, multiplier 12 → effective 6/s,
    # dt·effective = 0.6 > 0.5 → dt_rate_violation expected.
    failure_table = ProtrusionStateMultipliers(
        table={"growing": {"k_maturity_per_s": 12.0}}
    )

    return [
        Scenario(
            name="neutral_multipliers_baseline",
            description=(
                "All multipliers 1.0 — 6.3b reduces to 6.3a baseline. "
                "Demonstrates the deterministic-wrapper neutral path."
            ),
            centroid_um_xy=centroid,
            params=base_params,
            multipliers=neutral,
            initial_adhesions=(fa_centered,),
            protrusion_factory=_registry_growing,
            n_steps=15,
            expected_status="PASS",
        ),
        Scenario(
            name="growing_boost_2x",
            description=(
                "Single FA linked to a growing protrusion; multiplier "
                "2.0 / 1.5 boosts maturity / bind rate vs the neutral "
                "baseline."
            ),
            centroid_um_xy=centroid,
            params=base_params,
            multipliers=growing_boost,
            initial_adhesions=(fa_centered,),
            protrusion_factory=_registry_growing,
            n_steps=15,
            expected_status="PASS",
        ),
        Scenario(
            name="retracting_neutral",
            description=(
                "Single FA linked to a retracting protrusion; neutral "
                "multiplier table — no boost / no penalty (base 6.3a "
                "behavior)."
            ),
            centroid_um_xy=centroid,
            params=base_params,
            multipliers=neutral,
            initial_adhesions=(fa_centered,),
            protrusion_factory=_registry_retracting,
            n_steps=15,
            expected_status="PASS",
        ),
        Scenario(
            name="mixed_states_demo",
            description=(
                "Two FAs linked to different protrusions (growing + "
                "stalled). Per-FA effective rate differs; the dt-rate "
                "gate uses the global maximum across both FAs."
            ),
            centroid_um_xy=centroid,
            params=base_params,
            multipliers=mixed_boost,
            initial_adhesions=(fa_centered, fa_north),
            protrusion_factory=_registry_mixed,
            n_steps=15,
            expected_status="PASS",
        ),
        Scenario(
            name="dt_rate_violation_failure",
            description=(
                "Intentional FAIL: multiplier 12.0 inflates the "
                "effective rate past the dt-rate safety margin. "
                "Verifies the wrapper raises dt_rate_violation BEFORE "
                "stepping, the runner records the failure_kind, and "
                "the index aggregator flags this scenario as FAIL."
            ),
            centroid_um_xy=centroid,
            params=base_params,
            multipliers=failure_table,
            initial_adhesions=(fa_centered,),
            protrusion_factory=_registry_growing,
            n_steps=5,
            expected_status="FAIL",
        ),
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        help="Repository root used for git commit lookup and runs/ output.",
    )
    args = parser.parse_args()
    repo = os.path.abspath(args.repo)
    sys.path.insert(0, repo)

    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = os.path.join(repo, "runs", f"{timestamp}_protrusion_coupled_fa")
    os.makedirs(run_root, exist_ok=False)
    git_hash = _git_commit_hash(repo)

    scenarios = _build_scenarios()
    summary_index: list[dict] = []

    for scenario in scenarios:
        scenario_dir = os.path.join(run_root, scenario.name)
        run = _run_scenario(scenario, scenario_dir, git_hash)
        summary_index.append(
            {
                "scenario": scenario.name,
                "status": run.status,
                "expected_status": run.expected_status,
                "summary": run.summary_path,
                "metadata": run.metadata_path,
                "failure": run.failure,
            }
        )

    # Aggregate against expectation. An intentional FAIL scenario whose
    # recorded status matches expected_status is "as expected" and
    # contributes to PASS aggregate; only unexpected mismatches drive
    # exit 1.
    intentional_failures_observed = [
        entry["scenario"]
        for entry in summary_index
        if entry["status"] == "FAIL" and entry["expected_status"] == "FAIL"
    ]
    unexpected_failures = [
        entry["scenario"]
        for entry in summary_index
        if entry["status"] == "FAIL" and entry["expected_status"] != "FAIL"
    ]
    unexpected_passes = [
        entry["scenario"]
        for entry in summary_index
        if entry["status"] == "PASS" and entry["expected_status"] != "PASS"
    ]
    aggregate_status = (
        "FAIL" if (unexpected_failures or unexpected_passes) else "PASS"
    )
    index_path = os.path.join(run_root, "index.json")
    with open(index_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "git_commit_hash": git_hash,
                "scenarios": summary_index,
                "aggregate_status": aggregate_status,
                "intentional_failures_observed": intentional_failures_observed,
                "unexpected_failures": unexpected_failures,
                "unexpected_passes": unexpected_passes,
                "run_root": run_root,
            },
            fh,
            indent=2,
            default=str,
        )
    print(
        f"6.3b protrusion-coupled FA harness run {aggregate_status}. "
        f"artifacts: {run_root}"
    )
    print(f"index: {index_path}")
    for entry in summary_index:
        if entry["status"] == entry["expected_status"]:
            tag = (
                "[PASS]"
                if entry["status"] == "PASS"
                else "[FAIL-as-expected]"
            )
            note = ""
        else:
            tag = "[UNEXPECTED]"
            note = (
                f" (got {entry['status']}, expected "
                f"{entry['expected_status']})"
            )
        print(f"  {tag} {entry['scenario']}{note} -> {entry['summary']}")
    if intentional_failures_observed:
        print(
            f"intentional failure scenarios (as expected): "
            f"{', '.join(intentional_failures_observed)}"
        )
    if unexpected_failures:
        print(
            f"UNEXPECTED FAILURES: {', '.join(unexpected_failures)}",
            file=sys.stderr,
        )
    if unexpected_passes:
        print(
            f"UNEXPECTED PASSES: {', '.join(unexpected_passes)}",
            file=sys.stderr,
        )
    return 0 if aggregate_status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
