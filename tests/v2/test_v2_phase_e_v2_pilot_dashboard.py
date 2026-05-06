"""B-tier 테스트: Phase E v2 pilot dashboard renderer.

5 focused tests per Codex `id=1694` 가드레일 (artifact-schema 검증
위주, pixel-perfect 이미지 비교 안 함):

1. 정상 dashboard 생성 (PNG + HTML 비어있지 않음)
2. CSV 컬럼 schema 검증 (잘못된 schema는 ValueError)
3. metadata.json runner 필드 검증 (다른 runner는 거부)
4. CLI smoke (subprocess + 빈 PNG 아님)
5. 빈 diagnostic rows (zero-step 케이스) — 명시적 placeholder PNG 생성
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys

import pytest

from acs.v2.phase_e_v2_pilot_runner import (
    PhaseEV2PilotConfig,
    run_phase_e_v2_pilot,
)
from acs.v2.viz.phase_e_v2_pilot_dashboard import (
    PilotDashboardArtifacts,
    render_phase_e_v2_pilot_dashboard,
)


def test_dashboard_renders_png_and_html_from_pilot_run(tmp_path):
    """Test 1: 정상 케이스 — pilot runner 결과로부터 PNG + HTML 생성."""
    config = PhaseEV2PilotConfig(n_steps=3, dt_s=0.0, frame_interval=2, grid_n=4)
    pilot = run_phase_e_v2_pilot(str(tmp_path / "pilot"), config, git_commit_hash="dash1")

    artifacts = render_phase_e_v2_pilot_dashboard(pilot.output_dir, write_html=True)

    assert isinstance(artifacts, PilotDashboardArtifacts)
    assert os.path.isfile(artifacts.png_path)
    assert os.path.getsize(artifacts.png_path) > 0
    assert artifacts.html_path is not None
    assert os.path.isfile(artifacts.html_path)
    assert os.path.getsize(artifacts.html_path) > 0
    assert artifacts.n_diagnostic_rows == 3
    # n_frames metadata는 runner가 기록한 값 그대로 전달
    assert artifacts.n_frames >= 1

    with open(artifacts.html_path, encoding="utf-8") as fh:
        html = fh.read()
    # HTML rollup이 핵심 섹션을 포함
    assert "Phase E v2 pilot dashboard" in html
    assert "display-only" in html.lower()
    # PNG는 HTML과 같은 디렉토리에 있어 상대 경로로 참조됨
    assert "dashboard.png" in html
    # diagnostic 컬럼이 표 헤더로 들어감
    assert "v_active_um2" in html
    assert "multiplier_max" in html


def test_dashboard_rejects_csv_with_missing_columns(tmp_path):
    """Test 2: CSV schema 검증 — 핵심 컬럼이 빠지면 ValueError."""
    run_dir = tmp_path / "broken_csv"
    run_dir.mkdir()
    csv_path = run_dir / "diagnostics.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["step_index", "time_s"])
        writer.writeheader()
        writer.writerow({"step_index": 1, "time_s": 0.0})
    metadata_path = run_dir / "metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump({"runner": "phase_e_v2_pilot", "config": {}, "notes": []}, fh)

    with pytest.raises(ValueError, match="missing expected diagnostic columns"):
        render_phase_e_v2_pilot_dashboard(str(run_dir))


def test_dashboard_rejects_unrelated_runner_metadata(tmp_path):
    """Test 3: metadata.json runner 필드가 다르면 거부 — 다른 runner 결과는
    이 dashboard로 렌더하지 않음."""
    run_dir = tmp_path / "wrong_runner"
    run_dir.mkdir()
    csv_path = run_dir / "diagnostics.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "step_index",
                "time_s",
                "v_active_um2",
                "multiplier_min",
                "multiplier_max",
                "multiplier_mean",
                "max_traction_norm_nN_per_um2",
                "orientation_tensor_abs_max",
            ],
        )
        writer.writeheader()
    metadata_path = run_dir / "metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump({"runner": "some_other_runner", "config": {}}, fh)

    with pytest.raises(ValueError, match="phase_e_v2_pilot"):
        render_phase_e_v2_pilot_dashboard(str(run_dir))


def test_dashboard_cli_smoke(tmp_path):
    """Test 4: CLI smoke — pilot runner CLI + dashboard CLI 연쇄 호출."""
    run_dir = tmp_path / "cli_dash"
    # 1) Pilot runner CLI로 작은 출력 생성
    proc_run = subprocess.run(
        [
            sys.executable,
            "scripts/run_phase_e_v2_pilot.py",
            "--output-dir",
            str(run_dir),
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
    assert "Phase E v2 pilot PASS" in proc_run.stdout

    # 2) Dashboard CLI로 PNG + HTML 생성
    proc_dash = subprocess.run(
        [
            sys.executable,
            "scripts/render_phase_e_v2_pilot.py",
            "--run-dir",
            str(run_dir),
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "Phase E v2 pilot dashboard PASS" in proc_dash.stdout
    assert (run_dir / "dashboard.png").exists()
    assert (run_dir / "dashboard.png").stat().st_size > 0
    assert (run_dir / "dashboard.html").exists()
    assert (run_dir / "dashboard.html").stat().st_size > 0


def test_dashboard_handles_zero_step_pilot_with_placeholder_png(tmp_path):
    """Test 5: zero-step pilot run (n_steps=0) → 빈 diagnostic rows이지만
    PNG는 placeholder로 생성. CSV는 헤더만 있고 데이터 행 없음."""
    config = PhaseEV2PilotConfig(n_steps=0, dt_s=0.0, frame_interval=1, grid_n=3)
    pilot = run_phase_e_v2_pilot(str(tmp_path / "zero_step"), config, git_commit_hash="dash5")

    # Zero-step run은 0 diagnostic rows + 1 initial frame을 생성
    assert len(pilot.diagnostics) == 0
    assert len(pilot.frame_paths) == 1

    artifacts = render_phase_e_v2_pilot_dashboard(pilot.output_dir, write_html=False)

    assert artifacts.n_diagnostic_rows == 0
    assert os.path.isfile(artifacts.png_path)
    assert os.path.getsize(artifacts.png_path) > 0
    # write_html=False 시 HTML 미생성
    assert artifacts.html_path is None
