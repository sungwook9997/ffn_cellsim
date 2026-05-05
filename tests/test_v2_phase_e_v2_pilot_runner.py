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


def test_phase_e_v2_pilot_motility_off_default_preserves_existing_behavior(tmp_path):
    """Cycle 17 motility hook: default off → 기존 38b679f 동작 회귀 0개.
    2-FA anti-collapse fixture + multipliers exp(±1) 그대로 유지."""
    config = PhaseEV2PilotConfig(n_steps=2, dt_s=0.0, frame_interval=1, k_active=1.0)
    assert config.cell_motility == "off"
    assert config.smoke_displacement_um_per_step == 0.0

    result = run_phase_e_v2_pilot(str(tmp_path / "default_off"), config)

    # 2-FA anti-collapse fixture가 그대로 사용되었는지 확인
    multipliers = result.final_phase_e_result.ecm_to_fa_bias.multipliers_per_fa
    assert multipliers.shape == (2, 3)
    np.testing.assert_allclose(multipliers[0, :], np.exp(1.0), rtol=1e-12)
    np.testing.assert_allclose(multipliers[1, :], np.exp(-1.0), rtol=1e-12)


def test_phase_e_v2_pilot_motility_smoke_on_translates_cell_boundary(tmp_path):
    """Cycle 17 motility smoke: cell boundary가 traction 방향으로 translate."""
    config = PhaseEV2PilotConfig(
        n_steps=3,
        dt_s=0.0,
        frame_interval=1,
        k_active=1.0,
        cell_motility="translation_smoke",
        smoke_displacement_um_per_step=0.2,
    )
    result = run_phase_e_v2_pilot(str(tmp_path / "smoke_on"), config)

    # Single-FA fixture (traction (1, 0)) → smoke가 매 step boundary를 +x로 0.2um shift
    # 첫 frame (step 0)과 마지막 frame을 read해서 centroid 비교
    from acs.v2.output.frame_dump import read_frame
    f0 = read_frame(result.frame_paths[0])
    f_last = read_frame(result.frame_paths[-1])
    cell_0 = next(iter(f0.cluster.cells.values()))
    cell_last = next(iter(f_last.cluster.cells.values()))
    verts_0 = np.asarray(cell_0.measurement_boundary.vertices_xy_um, dtype=np.float64)
    verts_last = np.asarray(cell_last.measurement_boundary.vertices_xy_um, dtype=np.float64)
    centroid_0 = verts_0.mean(axis=0)
    centroid_last = verts_last.mean(axis=0)
    delta = centroid_last - centroid_0
    # n_steps=3 step 후 +x로 0.2*3 = 0.6um shift 예상 (frame_interval=1, 매 step보임)
    np.testing.assert_allclose(delta, [0.6, 0.0], atol=1e-9)


def test_phase_e_v2_pilot_motility_config_validation(tmp_path):
    """Cycle 17 motility config validation: 잘못된 값 → ValueError."""
    # cell_motility 값 검증
    with pytest.raises(ValueError, match="cell_motility"):
        PhaseEV2PilotConfig(cell_motility="invalid_mode").validate()  # type: ignore[arg-type]
    # smoke_displacement 음수 거부
    with pytest.raises(ValueError, match="smoke_displacement"):
        PhaseEV2PilotConfig(smoke_displacement_um_per_step=-0.1).validate()
    # smoke_displacement 비유한 거부
    with pytest.raises(ValueError, match="smoke_displacement"):
        PhaseEV2PilotConfig(smoke_displacement_um_per_step=float("nan")).validate()
    # translation_smoke + displacement=0 → ValueError (visible motion 안 보임)
    with pytest.raises(ValueError, match="smoke_displacement_um_per_step > 0"):
        PhaseEV2PilotConfig(
            cell_motility="translation_smoke",
            smoke_displacement_um_per_step=0.0,
        ).validate()


def test_phase_e_v2_pilot_motility_metadata_labels_visual_smoke(tmp_path):
    """Cycle 17 motility smoke: metadata에 visual_smoke_only/not_mechanistic 라벨."""
    config = PhaseEV2PilotConfig(
        n_steps=1,
        dt_s=0.0,
        frame_interval=1,
        cell_motility="translation_smoke",
        smoke_displacement_um_per_step=0.1,
    )
    result = run_phase_e_v2_pilot(str(tmp_path / "metadata_smoke"), config)

    with open(result.metadata_path, encoding="utf-8") as fh:
        metadata = json.load(fh)
    assert metadata["visual_smoke_only"] is True
    assert metadata["not_mechanistic"] is True
    assert metadata["cell_motility_mode"] == "translation_smoke"
    notes = metadata["notes"]
    assert any("visual smoke only" in n for n in notes)
    assert any("display parameter" in n.lower() or "display param" in n.lower() for n in notes)
    assert any("no FA-to-cell force coupling" in n for n in notes)


def test_phase_e_v2_pilot_motility_off_metadata_labels_correct(tmp_path):
    """Cycle 17 motility off: metadata에 visual_smoke_only=false, mode=off."""
    config = PhaseEV2PilotConfig(n_steps=1, dt_s=0.0, frame_interval=1)
    result = run_phase_e_v2_pilot(str(tmp_path / "metadata_off"), config)

    with open(result.metadata_path, encoding="utf-8") as fh:
        metadata = json.load(fh)
    assert metadata["visual_smoke_only"] is False
    assert metadata["not_mechanistic"] is True
    assert metadata["cell_motility_mode"] == "off"


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
