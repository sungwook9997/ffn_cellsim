from __future__ import annotations

import csv
import json
import subprocess
import sys

import numpy as np
import pytest

from acs.v2.output.frame_dump import read_frame
from acs.v2.phase_e_v2_pilot_runner import (
    PhaseEV2PilotConfig,
    make_phase_e_v2_pilot_fixture,
    run_phase_e_v2_pilot,
)


def test_phase_e_v2_pilot_fixture_matches_anti_collapse_geometry():
    ecm, cell, adhesions = make_phase_e_v2_pilot_fixture(grid_n=4, spacing_um=1.0)

    assert ecm.grid_shape == (4, 4)
    np.testing.assert_allclose(ecm.orientation_tensor[..., 0, 0], 1.0)
    np.testing.assert_allclose(ecm.orientation_tensor[..., 1, 1], -1.0)
    assert cell.cell_id == "pilot-cell"
    assert len(adhesions) == 2
    assert adhesions[0].position_um_xy == adhesions[1].position_um_xy
    assert adhesions[0].traction_force_nN_xy == pytest.approx((1.0, 0.0))
    assert adhesions[1].traction_force_nN_xy == pytest.approx((0.0, 1.0))


def test_phase_e_v2_pilot_writes_frames_csv_and_metadata(tmp_path):
    config = PhaseEV2PilotConfig(n_steps=3, dt_s=0.0, frame_interval=2, grid_n=4)
    result = run_phase_e_v2_pilot(
        str(tmp_path / "pilot"),
        config,
        git_commit_hash="abc123",
    )

    assert len(result.diagnostics) == 3
    assert [d.step_index for d in result.diagnostics] == [1, 2, 3]
    # Initial frame plus step 2 and final step 3.
    assert len(result.frame_paths) == 3

    restored = read_frame(result.frame_paths[-1])
    assert restored.cluster.ecm.grid_shape == (4, 4)
    assert tuple(restored.cluster.cells.keys()) == ("pilot-cell",)

    with open(result.csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 3
    assert rows[0]["step_index"] == "1"
    assert float(rows[0]["multiplier_max"]) > 1.0
    assert float(rows[0]["multiplier_min"]) < 1.0

    with open(result.metadata_path, encoding="utf-8") as fh:
        metadata = json.load(fh)
    assert metadata["runner"] == "phase_e_v2_pilot"
    assert metadata["tier"] == "B"
    assert metadata["git_commit_hash"] == "abc123"
    assert metadata["config"]["n_steps"] == 3
    assert metadata["n_diagnostic_rows"] == 3
    assert metadata["n_frames"] == 3
    assert "no new physics law" in metadata["notes"]


def test_phase_e_v2_pilot_preserves_exact_dt_zero_anti_collapse_values(tmp_path):
    config = PhaseEV2PilotConfig(n_steps=1, dt_s=0.0, frame_interval=1, k_active=1.0)
    result = run_phase_e_v2_pilot(str(tmp_path / "pilot_exact"), config)

    phase_e = result.final_phase_e_result
    assert phase_e is not None
    multipliers = phase_e.ecm_to_fa_bias.multipliers_per_fa
    np.testing.assert_allclose(multipliers[0, :], np.exp(1.0), rtol=1e-12)
    np.testing.assert_allclose(multipliers[1, :], np.exp(-1.0), rtol=1e-12)
    assert result.diagnostics[0].multiplier_max == pytest.approx(np.exp(1.0))
    assert result.diagnostics[0].multiplier_min == pytest.approx(np.exp(-1.0))


def test_phase_e_v2_pilot_config_rejects_bad_values():
    with pytest.raises(TypeError, match="n_steps"):
        PhaseEV2PilotConfig(n_steps=True).validate()  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="frame_interval"):
        PhaseEV2PilotConfig(frame_interval=0).validate()
    with pytest.raises(ValueError, match="k_active"):
        PhaseEV2PilotConfig(k_active=-1.0).validate()


def test_phase_e_v2_pilot_cli_smoke(tmp_path):
    out_dir = tmp_path / "cli_pilot"
    proc = subprocess.run(
        [
            sys.executable,
            "scripts/run_phase_e_v2_pilot.py",
            "--output-dir",
            str(out_dir),
            "--n-steps",
            "1",
            "--dt-s",
            "0",
            "--frame-interval",
            "1",
            "--grid",
            "3",
        ],
        check=True,
        text=True,
        capture_output=True,
    )

    assert "Phase E v2 pilot PASS" in proc.stdout
    assert (out_dir / "metadata.json").exists()
    assert (out_dir / "diagnostics.csv").exists()
    assert (out_dir / "frame_000000.h5").exists()
    assert (out_dir / "frame_000001.h5").exists()
