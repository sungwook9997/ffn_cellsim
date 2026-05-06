"""B-tier 테스트: Phase F minimal motility pilot runner.

5 focused tests:
1. fixture: single-FA + leftmost vertex + traction (-x) -> visible +x motion
2. run produces frames/csv/metadata + centroid moves +x
3. config validation: invalid values -> ValueError
4. metadata labels (visual_smoke_only=False, not_mechanistic=True,
   ecm_to_cell_feedback=False) + honest scope notes
5. CLI smoke
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys

import numpy as np
import pytest

from acs.v2.phase_f_minimal_motility_pilot_runner import (
    PhaseFPilotConfig,
    PhaseFPilotRunResult,
    make_phase_f_pilot_fixture,
    run_phase_f_minimal_motility_pilot,
)


def test_phase_f_pilot_fixture_single_fa_left_vertex_traction_negative_x():
    """Test 1: single-FA fixture default — FA at leftmost vertex,
    traction along -x (Newton 3rd reaction on cell vertex = +x)."""
    config = PhaseFPilotConfig(n_vertices=12, radius_um=1.0, fa_count=1)
    contour, ecm, adhesions, fa_to_vertex = make_phase_f_pilot_fixture(config)

    assert len(adhesions) == 1
    assert len(fa_to_vertex) == 1
    fa = adhesions[0]
    v_idx = fa_to_vertex[fa.adhesion_id]
    leftmost_x = float(contour.vertices_xy_um[:, 0].min())
    assert float(contour.vertices_xy_um[v_idx, 0]) == pytest.approx(leftmost_x, abs=1e-12)
    assert fa.traction_force_nN_xy == pytest.approx((-config.fa_traction_nN, 0.0))


def test_phase_f_pilot_run_writes_frames_csv_metadata_with_visible_x_motion(tmp_path):
    """Test 2: end-to-end pilot run produces all artifacts + centroid
    moves visibly +x (Phase F local pulling/deformation evidence)."""
    config = PhaseFPilotConfig(
        n_steps=50,
        dt_cell_s=1e-3,
        frame_interval=10,
        grid_n=8,
        spacing_um=0.5,
        n_vertices=16,
        radius_um=1.0,
        lambda_c_nN=0.05,
        sigma_c_nN_per_um=0.05,
        k_active=1.0,
        fa_traction_nN=0.5,
        fa_count=1,
    )
    result = run_phase_f_minimal_motility_pilot(
        str(tmp_path / "phase_f_run"), config, git_commit_hash="test-pf"
    )

    assert isinstance(result, PhaseFPilotRunResult)
    assert os.path.isfile(result.csv_path)
    assert os.path.isfile(result.metadata_path)
    assert len(result.diagnostics) == 50
    # Initial frame + 5 step frames (50 / 10) = 6 frames
    assert len(result.frame_paths) == 6

    # Centroid moves +x (FA pulls cell along -x; Newton 3rd reaction +x)
    first = result.diagnostics[0]
    last = result.diagnostics[-1]
    assert last.centroid_x_um > first.centroid_x_um

    # Attached vertex peak displacement > non-attached peak (local pulling
    # NOT rigid translation — PI id=1949 acceptance)
    attached_peak = max(
        d.attached_vertex_displacement_um for d in result.diagnostics
    )
    non_attached_peak = max(
        d.max_non_attached_vertex_displacement_um for d in result.diagnostics
    )
    assert attached_peak > non_attached_peak


def test_phase_f_pilot_config_validation_rejects_invalid():
    """Test 3: config validation — invalid values -> ValueError."""
    with pytest.raises(ValueError, match="dt_cell_s"):
        PhaseFPilotConfig(dt_cell_s=0.0).validate()
    with pytest.raises(ValueError, match="k_active"):
        PhaseFPilotConfig(k_active=-1.0).validate()
    with pytest.raises(ValueError, match="fa_count"):
        PhaseFPilotConfig(n_vertices=4, fa_count=10).validate()
    with pytest.raises(ValueError, match="fa_traction_nN"):
        PhaseFPilotConfig(fa_traction_nN=float("nan")).validate()


def test_phase_f_pilot_metadata_labels_honest_scope(tmp_path):
    """Test 4: metadata labels visual_smoke_only=False, not_mechanistic=True,
    ecm_to_cell_feedback=False; notes include honest-scope wording."""
    config = PhaseFPilotConfig(n_steps=10, frame_interval=5)
    result = run_phase_f_minimal_motility_pilot(
        str(tmp_path / "metadata"), config, git_commit_hash="test-md"
    )
    with open(result.metadata_path, encoding="utf-8") as fh:
        metadata = json.load(fh)
    assert metadata["runner"] == "phase_f_minimal_motility_pilot"
    assert metadata["visual_smoke_only"] is False
    assert metadata["not_mechanistic"] is True
    assert metadata["ecm_to_cell_feedback"] is False
    assert metadata["physics_layer"] == "FA-vertex Newton 3rd law (locked)"
    notes = metadata["notes"]
    assert any("FA traction drives attached cell vertex" in n for n in notes)
    assert any("HB#4 rate multipliers are diagnostic-only" in n for n in notes)
    assert any(
        "does not yet implement ECM->cell motility feedback" in n for n in notes
    )
    assert any("no PI experimental data use" in n for n in notes)


def test_phase_f_pilot_cli_smoke(tmp_path):
    """Test 5: CLI smoke — runner CLI generates all artifacts."""
    out_dir = tmp_path / "cli_pf"
    proc = subprocess.run(
        [
            sys.executable,
            "scripts/run_phase_f_minimal_motility_pilot.py",
            "--output-dir",
            str(out_dir),
            "--n-steps",
            "20",
            "--frame-interval",
            "5",
            "--n-vertices",
            "12",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "Phase F minimal motility pilot PASS" in proc.stdout
    assert (out_dir / "metadata.json").exists()
    assert (out_dir / "diagnostics.csv").exists()
    # Initial + 4 step frames (20/5) = 5 frames
    assert (out_dir / "frame_000000.h5").exists()
    assert (out_dir / "frame_000020.h5").exists()
