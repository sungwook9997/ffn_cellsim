"""B-tier 테스트: Phase E v2 pilot frame replay viewer.

5 focused tests per Codex `id=1728` 가드레일 (frame-count / schema /
path / non-empty 검증; pixel-perfect 비교 안 함):

1. 정상 replay 생성 (PNG + HTML 비어있지 않음, frame 수 전달)
2. metadata.json runner 필드 mismatch → ValueError
3. frame_*.h5 없는 run_dir → FileNotFoundError
4. max_frames subsampling (frame 수 > max → 균등-간격 subsample)
5. CLI smoke (subprocess + non-empty output)
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

from acs.v2.phase_e_v2_pilot_runner import (
    PhaseEV2PilotConfig,
    run_phase_e_v2_pilot,
)
from acs.v2.viz.phase_e_v2_pilot_frame_replay import (
    PilotFrameReplayArtifacts,
    render_phase_e_v2_pilot_frame_replay,
)


def test_frame_replay_renders_grid_png_and_html_from_pilot_run(tmp_path):
    """Test 1: 정상 케이스 — pilot runner 결과로부터 multi-panel PNG + HTML."""
    config = PhaseEV2PilotConfig(n_steps=3, dt_s=0.0, frame_interval=1, grid_n=4)
    pilot = run_phase_e_v2_pilot(str(tmp_path / "replay1"), config, git_commit_hash="r1")

    artifacts = render_phase_e_v2_pilot_frame_replay(pilot.output_dir, write_html=True)

    assert isinstance(artifacts, PilotFrameReplayArtifacts)
    assert os.path.isfile(artifacts.png_path)
    assert os.path.getsize(artifacts.png_path) > 0
    assert artifacts.html_path is not None
    assert os.path.isfile(artifacts.html_path)
    assert os.path.getsize(artifacts.html_path) > 0
    # 4 frames total (initial + 3 steps with frame_interval=1)
    assert artifacts.n_frames_total == 4
    assert artifacts.n_frames_rendered == 4
    assert len(artifacts.frame_paths_rendered) == 4

    with open(artifacts.html_path, encoding="utf-8") as fh:
        html = fh.read()
    assert "Phase E v2 pilot frame replay" in html
    assert "display-only" in html.lower()
    assert "frame_replay.png" in html
    # per-frame index 표 헤더 + 각 frame basename
    assert "<th>frame</th>" in html
    assert "<th>time_s</th>" in html
    for path in artifacts.frame_paths_rendered:
        assert os.path.basename(path) in html


def test_frame_replay_rejects_unrelated_runner_metadata(tmp_path):
    """Test 2: metadata.json runner != phase_e_v2_pilot → ValueError."""
    run_dir = tmp_path / "wrong_runner"
    run_dir.mkdir()
    # 빈 frame 파일 하나라도 있어야 metadata 단계까지 도달; 하지만 metadata
    # 검사가 먼저 일어나야 함 (lock 안 된 schema 보호).
    metadata_path = run_dir / "metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump({"runner": "some_other_runner", "config": {}, "notes": []}, fh)

    with pytest.raises(ValueError, match="phase_e_v2_pilot"):
        render_phase_e_v2_pilot_frame_replay(str(run_dir))


def test_frame_replay_raises_when_no_frames(tmp_path):
    """Test 3: metadata 정상이지만 frame_*.h5 없으면 FileNotFoundError."""
    run_dir = tmp_path / "no_frames"
    run_dir.mkdir()
    metadata_path = run_dir / "metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump({"runner": "phase_e_v2_pilot", "config": {}, "notes": []}, fh)

    with pytest.raises(FileNotFoundError, match="no frame_.*\\.h5 files"):
        render_phase_e_v2_pilot_frame_replay(str(run_dir))


def test_frame_replay_subsamples_when_exceeding_max_frames(tmp_path):
    """Test 4: frame 수 > max_frames → 균등-간격 subsample, 첫/마지막 보존."""
    config = PhaseEV2PilotConfig(n_steps=10, dt_s=0.0, frame_interval=1, grid_n=3)
    pilot = run_phase_e_v2_pilot(str(tmp_path / "replay4"), config, git_commit_hash="r4")
    # initial frame + 10 step frames = 11 total
    assert len(pilot.frame_paths) == 11

    artifacts = render_phase_e_v2_pilot_frame_replay(
        pilot.output_dir, write_html=False, max_frames=4
    )
    assert artifacts.n_frames_total == 11
    assert artifacts.n_frames_rendered == 4
    # 첫번째와 마지막 frame이 보존됨 (균등-간격 subsample)
    rendered_basenames = [os.path.basename(p) for p in artifacts.frame_paths_rendered]
    all_basenames = sorted(os.path.basename(p) for p in pilot.frame_paths)
    assert rendered_basenames[0] == all_basenames[0]  # frame_000000.h5
    assert rendered_basenames[-1] == all_basenames[-1]  # frame_000010.h5


def test_frame_replay_cli_smoke(tmp_path):
    """Test 5: CLI smoke — pilot runner CLI + replay CLI 연쇄."""
    run_dir = tmp_path / "cli_replay"
    proc_run = subprocess.run(
        [
            sys.executable,
            "scripts/run_phase_e_v2_pilot.py",
            "--output-dir",
            str(run_dir),
            "--n-steps",
            "2",
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

    proc_replay = subprocess.run(
        [
            sys.executable,
            "scripts/replay_phase_e_v2_pilot.py",
            "--run-dir",
            str(run_dir),
            "--max-frames",
            "6",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "Phase E v2 pilot frame replay PASS" in proc_replay.stdout
    assert (run_dir / "frame_replay.png").exists()
    assert (run_dir / "frame_replay.png").stat().st_size > 0
    assert (run_dir / "frame_replay.html").exists()
    assert (run_dir / "frame_replay.html").stat().st_size > 0
