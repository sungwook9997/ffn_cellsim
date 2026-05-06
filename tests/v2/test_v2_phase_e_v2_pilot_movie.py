"""B-tier 테스트: Phase E v2 pilot movie generator.

5 focused tests per Codex `id=1783` Cycle 15 scope:

1. imageio 백엔드 가용 여부와 무관하게 graceful 결과 (PilotMovieArtifacts
   반환; mp4/gif는 None 가능)
2. metadata.json runner 필드 mismatch → ValueError
3. frame_*.h5 없는 run_dir → FileNotFoundError
4. fps/gif_fps 비양수 → ValueError
5. CLI smoke + graceful "backend unavailable" 메시지 처리

테스트는 pixel-perfect 비교 안 하며, 의도적으로 imageio가 없을 수 있는
venv에서도 PASS하도록 graceful semantics 검증 위주.
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
from acs.v2.viz.phase_e_v2_pilot_movie import (
    PilotMovieArtifacts,
    render_phase_e_v2_pilot_movie,
)


def test_movie_renders_artifacts_dataclass_from_pilot_run(tmp_path):
    """Test 1: 정상 케이스 — PilotMovieArtifacts 반환, n_frames 일치.
    imageio 가용 시 mp4_path 비어있지 않음; 미가용 시 mp4_path=None.
    어느 경우든 graceful (raise 없음)."""
    config = PhaseEV2PilotConfig(n_steps=2, dt_s=0.0, frame_interval=1, grid_n=3)
    pilot = run_phase_e_v2_pilot(str(tmp_path / "movie1"), config, git_commit_hash="m1")

    artifacts = render_phase_e_v2_pilot_movie(pilot.output_dir)
    assert isinstance(artifacts, PilotMovieArtifacts)
    assert artifacts.n_frames_total == 3
    assert artifacts.n_frames_rendered == 3
    assert artifacts.fps == 2
    assert artifacts.gif_fps == 2

    if artifacts.backend_available and artifacts.mp4_path is not None:
        # imageio + ffmpeg 환경에서 mp4 산출됨
        assert os.path.isfile(artifacts.mp4_path)
        assert os.path.getsize(artifacts.mp4_path) > 0
    else:
        # 백엔드 미가용 → graceful None
        assert artifacts.mp4_path is None
        assert artifacts.gif_path is None

    # PNG temporary cleanup 검증 (keep_frame_pngs=False default)
    movie_frames_dir = os.path.join(pilot.output_dir, "movie_frames")
    assert not os.path.isdir(movie_frames_dir)


def test_movie_rejects_unrelated_runner_metadata(tmp_path):
    """Test 2: metadata.json runner != phase_e_v2_pilot → ValueError."""
    run_dir = tmp_path / "wrong_runner"
    run_dir.mkdir()
    metadata_path = run_dir / "metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump({"runner": "some_other_runner", "config": {}, "notes": []}, fh)

    with pytest.raises(ValueError, match="phase_e_v2_pilot"):
        render_phase_e_v2_pilot_movie(str(run_dir))


def test_movie_raises_when_no_frames(tmp_path):
    """Test 3: metadata 정상이지만 frame_*.h5 없으면 FileNotFoundError."""
    run_dir = tmp_path / "no_frames"
    run_dir.mkdir()
    metadata_path = run_dir / "metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump({"runner": "phase_e_v2_pilot", "config": {}, "notes": []}, fh)

    with pytest.raises(FileNotFoundError, match="no frame_.*\\.h5 files"):
        render_phase_e_v2_pilot_movie(str(run_dir))


def test_movie_rejects_non_positive_fps(tmp_path):
    """Test 4: fps/gif_fps ≤ 0 → ValueError (CLI 잘못된 인자 방지)."""
    config = PhaseEV2PilotConfig(n_steps=1, dt_s=0.0, frame_interval=1, grid_n=3)
    pilot = run_phase_e_v2_pilot(str(tmp_path / "fps"), config, git_commit_hash="m4")

    with pytest.raises(ValueError, match="fps"):
        render_phase_e_v2_pilot_movie(pilot.output_dir, fps=0)
    with pytest.raises(ValueError, match="fps"):
        render_phase_e_v2_pilot_movie(pilot.output_dir, fps=-1)
    with pytest.raises(ValueError, match="gif_fps"):
        render_phase_e_v2_pilot_movie(pilot.output_dir, gif_fps=0)


def test_movie_classifies_mp4_codec_fallback_gif_correctly(tmp_path, monkeypatch):
    """Test 6 (Codex `id=1814` BLOCKER 2 fix regression): v1
    `write_movie(..., gif=False)`가 MP4 codec 실패 시 `.gif` path를
    `result['mp4']`에 반환하는 케이스 — wrapper가 suffix로 분류해서
    `mp4_path`가 `.gif`를 가리키지 않도록 보장."""
    config = PhaseEV2PilotConfig(n_steps=1, dt_s=0.0, frame_interval=1, grid_n=3)
    pilot = run_phase_e_v2_pilot(str(tmp_path / "fallback"), config, git_commit_hash="m6")

    from acs.v2.viz import phase_e_v2_pilot_movie as movie_mod

    # imageio backend가 가용하다고 강제 표시 + write_mp4_and_gif가
    # codec fallback 시나리오를 시뮬레이션 (mp4 슬롯에 .gif 반환).
    fake_gif_path = str(tmp_path / "fallback" / "pilot_movie.gif")

    def _fake_write_mp4_and_gif(png_paths, mp4_path, *, fps, gif_fps, gif_max_frames):
        # Codec fallback: write a .gif but place it in the mp4 slot
        # (mirrors v1 write_movie's gif fallback behavior).
        with open(fake_gif_path, "wb") as fh:
            fh.write(b"\x00\x01gif-bytes")
        from pathlib import Path as _Path
        return {"mp4": _Path(fake_gif_path), "gif": _Path(fake_gif_path)}

    # Monkeypatch the live module attributes directly so the wrapper's
    # lazy import sees backend as available + the fake helper.
    import acs.visualization.live_imaging as _live_imaging_real
    monkeypatch.setattr(_live_imaging_real, "imageio", object())
    monkeypatch.setattr(_live_imaging_real, "write_mp4_and_gif", _fake_write_mp4_and_gif)

    artifacts = render_phase_e_v2_pilot_movie(pilot.output_dir)
    # mp4_path는 .gif suffix를 가리키지 않아야 함
    assert artifacts.mp4_path is None or artifacts.mp4_path.endswith(".mp4")
    # write_gif=True (default) → gif_path가 .gif를 가리킴
    assert artifacts.gif_path is not None
    assert artifacts.gif_path.endswith(".gif")
    assert os.path.isfile(artifacts.gif_path)


def test_movie_cli_smoke_with_graceful_backend_handling(tmp_path):
    """Test 5: CLI smoke — pilot runner CLI + movie CLI 연쇄. backend
    가용 여부와 무관하게 exit code 0; 백엔드 미가용 시 NOTE 메시지 출력."""
    run_dir = tmp_path / "cli_movie"
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

    proc_movie = subprocess.run(
        [
            sys.executable,
            "scripts/movie_phase_e_v2_pilot.py",
            "--run-dir",
            str(run_dir),
            "--fps",
            "2",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "Phase E v2 pilot movie PASS" in proc_movie.stdout
    assert "backend_available:" in proc_movie.stdout
    assert "n_frames_total: 2" in proc_movie.stdout

    if "backend_available: False" in proc_movie.stdout:
        # 백엔드 미가용 시 graceful skip 메시지
        assert "NOTE: imageio backend unavailable" in proc_movie.stdout
    else:
        # 백엔드 가용 시 mp4 파일 존재 (imageio + ffmpeg 환경)
        assert (run_dir / "pilot_movie.mp4").exists() or "mp4: None" in proc_movie.stdout
