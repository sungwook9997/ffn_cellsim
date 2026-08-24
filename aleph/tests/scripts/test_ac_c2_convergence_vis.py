"""C-2 diagnostic visualization must remain a solver plot, never an F-delta curve."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from aleph.scripts.ac_c2_convergence_vis import load_trace, render


def _record(path: Path, *, tolerance: float = 0.2, objective: str | None = None) -> Path:
    payload = {
        "status": "C2_REJECTED_ROLLED_BACK",
        "solver": {
            "name": "erm_gauss_seidel",
            "iteration_budget": 40,
            "projected_force_tolerance_pN": tolerance,
            "convergence_history": {
                "iteration": [20, 40],
                "projected_force_pN": [0.8, 0.7],
            },
        },
    }
    if objective is not None:
        payload["solver"]["line_search_objective"] = objective
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_load_trace_preserves_rejected_status_and_solver_axes(tmp_path: Path) -> None:
    label, iteration, force, tolerance, status = load_trace(_record(tmp_path / "trace.json"))
    assert label == "erm_gauss_seidel (40-iteration budget)"
    assert iteration == [20, 40]
    assert force == [0.8, 0.7]
    assert tolerance == 0.2
    assert status == "C2_REJECTED_ROLLED_BACK"


def test_load_trace_labels_nondefault_candidate_objective(tmp_path: Path) -> None:
    label, *_ = load_trace(_record(tmp_path / "trace.json", objective="l2_squared"))
    assert label.endswith("[l2_squared candidate selection]")


def test_render_writes_a_solver_figure(tmp_path: Path) -> None:
    output = tmp_path / "c2.png"
    render([_record(tmp_path / "trace.json")], output)
    assert output.read_bytes().startswith(b"\x89PNG")


def test_traces_with_different_gate_thresholds_are_not_overlaid(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="different projected-force thresholds"):
        render(
            [_record(tmp_path / "a.json", tolerance=0.2),
             _record(tmp_path / "b.json", tolerance=0.3)],
            tmp_path / "bad.png",
        )
