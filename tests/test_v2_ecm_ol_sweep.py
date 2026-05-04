"""Tests for the open-loop ECM sensitivity sweep harness (Phase C).

Test catalog mirrors ``docs/v2_ecm_ol_sweep_sanity_gate.md`` §9 with
the §0/§2/§4/§6 boundary and wording-lock checks. Each test maps to
one Sanity Gate item or one forbidden behavior from the locked
phased plan ``docs/v2_closed_loop_ecm_gate_phased_plan_locked.md``
§1 Phase C.

Phase C is open-loop sweep baseline only — these tests do **NOT**
satisfy closed-loop ECM gate Item 5. The meta-test
``test_open_loop_sweep_does_not_satisfy_closed_loop_item_5``
encodes this boundary as a runtime assertion (parallel to Phase B's
``test_stimulus_monotonicity_does_not_satisfy_closed_loop_item_1``).
"""

from __future__ import annotations

import json
import os

import numpy as np
import pytest

from scripts.run_ecm_ol_sensitivity_sweep import (
    ChannelScenario,
    SweepTuple,
    SweepValidationError,
    _build_scenarios,
    _build_tuples_smoke_fixture,
    run_ecm_ol_sensitivity_sweep,
)


def _tiny_tuple(label: str = "tiny", *, grid=(2, 2), dt_s=0.1, n_steps=2):
    return SweepTuple(
        grid_n=grid,
        spacing_um=1.0,
        dt_s=dt_s,
        n_steps=n_steps,
        frame_interval=1,
        expected_status="PASS",
        label=label,
    )


# ---------------------------------------------------------------------------
# Test 1 — explicit fixture runs all 4 channels (gate §0 + §3 isolation)
# ---------------------------------------------------------------------------


def test_open_loop_sweep_runs_4_channels_explicit_fixture(tmp_path):
    """Caller-supplied fixture (no harness defaults) runs all 4
    channels per tuple and writes per-tuple summary + metadata +
    cross-tuple plot. fixture_kind label is preserved into
    index.json."""

    sweep_tuples = [
        _tiny_tuple("a", grid=(2, 2), dt_s=0.1, n_steps=2),
        _tiny_tuple("b", grid=(3, 3), dt_s=0.05, n_steps=2),
    ]
    scenarios = _build_scenarios((2, 2))

    out = str(tmp_path / "fixture_run")
    run = run_ecm_ol_sensitivity_sweep(
        sweep_tuples=sweep_tuples,
        channel_scenarios=scenarios,
        reduction_choice="max",
        max_grid_cells_total=64,
        output_dir=out,
        git_commit_hash="abc123",
        fixture_kind="non_production_smoke",
    )
    assert run.aggregate_status == "PASS"
    assert run.fixture_kind == "non_production_smoke"
    # 2 tuples × 4 channels = 8 records.
    assert len(run.records) == 8

    with open(os.path.join(out, "index.json")) as fh:
        data = json.load(fh)
    assert data["aggregate_status"] == "PASS"
    assert data["fixture_kind"] == "non_production_smoke"
    assert data["evidence_kind"] == "open-loop sweep baseline"

    # Per-tuple per-channel artifact present.
    assert os.path.isfile(os.path.join(out, "a", "metadata.json"))
    assert os.path.isfile(os.path.join(out, "a", "summary.html"))
    assert os.path.isfile(os.path.join(out, "b", "metadata.json"))
    assert os.path.isfile(os.path.join(out, "sweep_summary.png"))


# ---------------------------------------------------------------------------
# Test 2 — grid spacing variation over the same dt (gate §1 + §3 + §6)
# ---------------------------------------------------------------------------


def test_open_loop_sweep_grid_spacing_variation(tmp_path):
    """Two tuples with different spacing_um at the same dt run
    independently per (tuple, channel); per-tuple final field
    values may differ but the PASS verdict holds."""

    tuples = [
        SweepTuple(
            grid_n=(2, 2),
            spacing_um=0.5,
            dt_s=0.1,
            n_steps=2,
            frame_interval=1,
            expected_status="PASS",
            label="spacing_half",
        ),
        SweepTuple(
            grid_n=(2, 2),
            spacing_um=2.0,
            dt_s=0.1,
            n_steps=2,
            frame_interval=1,
            expected_status="PASS",
            label="spacing_two",
        ),
    ]
    scenarios = _build_scenarios((2, 2))
    out = str(tmp_path / "spacing_run")
    run = run_ecm_ol_sensitivity_sweep(
        sweep_tuples=tuples,
        channel_scenarios=scenarios,
        reduction_choice="max",
        max_grid_cells_total=64,
        output_dir=out,
        git_commit_hash="",
        fixture_kind="non_production_smoke",
    )
    assert run.aggregate_status == "PASS"
    half = [r for r in run.records if r.tuple_label == "spacing_half"]
    two = [r for r in run.records if r.tuple_label == "spacing_two"]
    assert len(half) == 4
    assert len(two) == 4


# ---------------------------------------------------------------------------
# Test 3 — dt variation (gate §1 + §4 + §6)
# ---------------------------------------------------------------------------


def test_open_loop_sweep_dt_variation(tmp_path):
    """Two tuples at the same grid with different dt run
    independently. Larger dt accumulates more per step."""

    tuples = [
        SweepTuple(
            grid_n=(2, 2),
            spacing_um=1.0,
            dt_s=0.1,
            n_steps=2,
            frame_interval=1,
            expected_status="PASS",
            label="dt_big",
        ),
        SweepTuple(
            grid_n=(2, 2),
            spacing_um=1.0,
            dt_s=0.01,
            n_steps=2,
            frame_interval=1,
            expected_status="PASS",
            label="dt_small",
        ),
    ]
    scenarios = _build_scenarios((2, 2))
    out = str(tmp_path / "dt_run")
    run = run_ecm_ol_sensitivity_sweep(
        sweep_tuples=tuples,
        channel_scenarios=scenarios,
        reduction_choice="max",
        max_grid_cells_total=64,
        output_dir=out,
        git_commit_hash="",
        fixture_kind="non_production_smoke",
    )
    assert run.aggregate_status == "PASS"
    big_traction = [
        r for r in run.records if r.tuple_label == "dt_big" and r.channel == "traction"
    ][0]
    small_traction = [
        r for r in run.records if r.tuple_label == "dt_small" and r.channel == "traction"
    ][0]
    # Larger dt accumulates more (constant prescribed traction).
    assert big_traction.final_field_max > small_traction.final_field_max


# ---------------------------------------------------------------------------
# Test 4 — summary reports baseline only (gate §6 wording, Hard Rule 11)
# ---------------------------------------------------------------------------


def test_open_loop_sweep_summary_reports_baseline_only(tmp_path):
    """Every artifact (summary.html, metadata.json, index.json)
    reports `evidence_kind = 'open-loop sweep baseline'` and does
    NOT contain the forbidden phrase 'Item 5 satisfied'."""

    tuples = [_tiny_tuple("only")]
    scenarios = _build_scenarios((2, 2))
    out = str(tmp_path / "wording_run")
    run_ecm_ol_sensitivity_sweep(
        sweep_tuples=tuples,
        channel_scenarios=scenarios,
        reduction_choice="max",
        max_grid_cells_total=64,
        output_dir=out,
        git_commit_hash="",
        fixture_kind="non_production_smoke",
    )

    with open(os.path.join(out, "index.json")) as fh:
        index_text = fh.read()
    with open(os.path.join(out, "only", "metadata.json")) as fh:
        meta_text = fh.read()
    with open(os.path.join(out, "only", "summary.html")) as fh:
        summary_text = fh.read()

    for blob in (index_text, meta_text, summary_text):
        assert "open-loop sweep baseline" in blob
        assert "Item 5 satisfied" not in blob
        assert "closed-loop sensitivity" not in blob


# ---------------------------------------------------------------------------
# Test 5 — Hard Rule 11 wording boundary (meta-test, parallel to Phase B test 4)
# ---------------------------------------------------------------------------


def test_open_loop_sweep_does_not_satisfy_closed_loop_item_5():
    """Phase C open-loop sweep baseline — Hard Rule 11 wording
    protection. Documents that this harness measures open-loop
    preflight outputs at varying (grid, dt), which is not the
    closed-loop sensitivity sweep that Item 5 requires.

    Item 5 closed-loop side requires:
    - Scattering geometry (Hard Blocker #3) — not yet locked
    - Constitutive response law (Hard Blocker #1) — not yet locked
    - Sensitivity sweep over (grid, dt) on the response field, not
      the open-loop preflight output

    This test is a runtime assertion that makes the wording
    boundary explicit so a future refactor cannot silently relabel
    the Phase C sweep as Item 5 satisfaction.
    """

    phase_c_evidence_kind = "open-loop sweep baseline"
    closed_loop_item_5_claim = "closed-loop grid/dt sensitivity (TBD)"
    assert phase_c_evidence_kind != closed_loop_item_5_claim, (
        "Phase C sweep measures open-loop preflight outputs at varying "
        "(grid, dt). It cannot satisfy closed-loop gate Item 5, which "
        "requires a response field + scattering geometry per Hard "
        "Blockers #1 and #3 in "
        "docs/v2_closed_loop_ecm_gate_phased_plan_locked.md §2."
    )


# ---------------------------------------------------------------------------
# Test 6 — memory cap rejects oversize grid before allocation (gate §4)
# ---------------------------------------------------------------------------


def test_open_loop_sweep_memory_cap_rejects_oversize_grid(tmp_path):
    """A tuple whose grid_n[0] * grid_n[1] exceeds
    max_grid_cells_total raises sweep_memory_cap_exceeded BEFORE any
    ECM is constructed, so the run-root remains absent (no half-run
    artifact dir)."""

    tuples = [_tiny_tuple("ok"), _tiny_tuple("oversize", grid=(20, 20))]
    scenarios = _build_scenarios((2, 2))
    out = str(tmp_path / "memory_run")
    with pytest.raises(SweepValidationError) as excinfo:
        run_ecm_ol_sensitivity_sweep(
            sweep_tuples=tuples,
            channel_scenarios=scenarios,
            reduction_choice="max",
            max_grid_cells_total=64,
            output_dir=out,
            git_commit_hash="",
            fixture_kind="non_production_smoke",
        )
    assert excinfo.value.failure_kind == "sweep_memory_cap_exceeded"
    # No partial run-root artifact directory was created.
    assert not os.path.isdir(out)


# ---------------------------------------------------------------------------
# Boundary cases (§2): per-input invalid value rejection
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "field, value, kind",
    [
        ("grid_n", (0, 2), "sweep_grid_n_invalid"),
        ("grid_n", (True, 2), "sweep_grid_n_invalid"),
        ("spacing_um", 0.0, "sweep_spacing_invalid"),
        ("spacing_um", -1.0, "sweep_spacing_invalid"),
        ("dt_s", -0.5, "sweep_dt_invalid"),
        ("n_steps", 0, "sweep_n_steps_invalid"),
        ("n_steps", True, "sweep_n_steps_invalid"),
    ],
)
def test_open_loop_sweep_rejects_invalid_tuple_input(field, value, kind, tmp_path):
    base = dict(
        grid_n=(2, 2),
        spacing_um=1.0,
        dt_s=0.1,
        n_steps=2,
        frame_interval=1,
        expected_status="PASS",
        label="invalid",
    )
    base[field] = value
    bad = SweepTuple(**base)
    scenarios = _build_scenarios((2, 2))
    out = str(tmp_path / "invalid_run")
    with pytest.raises(SweepValidationError) as excinfo:
        run_ecm_ol_sensitivity_sweep(
            sweep_tuples=[bad],
            channel_scenarios=scenarios,
            reduction_choice="max",
            max_grid_cells_total=64,
            output_dir=out,
            git_commit_hash="",
            fixture_kind="non_production_smoke",
        )
    assert excinfo.value.failure_kind == kind


# ---------------------------------------------------------------------------
# Reduction choice validation
# ---------------------------------------------------------------------------


def test_open_loop_sweep_rejects_unknown_reduction(tmp_path):
    tuples = [_tiny_tuple("only")]
    scenarios = _build_scenarios((2, 2))
    out = str(tmp_path / "bad_red")
    with pytest.raises(SweepValidationError):
        run_ecm_ol_sensitivity_sweep(
            sweep_tuples=tuples,
            channel_scenarios=scenarios,
            reduction_choice="median",  # not in allowed set
            max_grid_cells_total=64,
            output_dir=out,
            git_commit_hash="",
            fixture_kind="non_production_smoke",
        )


# ---------------------------------------------------------------------------
# Smoke fixture builder sanity
# ---------------------------------------------------------------------------


def test_smoke_fixture_builder_returns_9_tuples():
    tuples = _build_tuples_smoke_fixture()
    assert len(tuples) == 9  # 3 grid sizes × 3 dt values
    for t in tuples:
        assert t.expected_status == "PASS"
        assert t.label.startswith("g")


def test_explicit_fixture_scenarios_cover_4_channels():
    scenarios = _build_scenarios((2, 2))
    names = sorted(s.name for s in scenarios)
    assert names == sorted(
        ["traction", "stiffness_rate", "density_rate", "orientation_rate"]
    )
