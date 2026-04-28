"""Stage 0 smoke tests — verify the harness is wired correctly.

These do not validate any physics. They check:
- The `acs` package imports cleanly.
- The GPU dispatcher resolves backends without touching Taichi.
- Config + provenance helpers produce the expected shapes.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import acs
from acs import config, gpu, provenance


def test_package_importable():
    assert hasattr(acs, "__version__")


@pytest.mark.parametrize(
    "value, expected",
    [
        ("cuda", "cuda"),
        ("CUDA", "cuda"),
        ("vulkan", "vulkan"),
        ("cpu", "cpu"),
        ("auto", "auto"),
        (None, "auto"),
    ],
)
def test_gpu_backend_resolution(monkeypatch, value, expected):
    monkeypatch.delenv("ACS_GPU_BACKEND", raising=False)
    if value is None:
        assert gpu._resolve_backend(None) == expected
    else:
        assert gpu._resolve_backend(value) == expected


def test_gpu_backend_env_var(monkeypatch):
    monkeypatch.setenv("ACS_GPU_BACKEND", "vulkan")
    assert gpu._resolve_backend(None) == "vulkan"


def test_gpu_backend_rejects_unknown():
    with pytest.raises(ValueError):
        gpu._resolve_backend("rocm")  # not yet supported


def test_load_config_dev():
    cfg = config.load_config(Path("configs/dev.yaml"))
    assert cfg["run"]["name"] == "dev"
    # Stage 0: only L1 should be on; everything else off.
    assert cfg["layers"]["L1_bulk_hydrodynamics"] is True
    assert cfg["layers"]["L2_boundary_biology"] is False


def test_provenance_manifest(tmp_path):
    fake_cfg = {"run": {"name": "smoke"}, "layers": {"L1_bulk_hydrodynamics": True}}
    manifest = provenance.build_manifest(fake_cfg, repo_root=".")
    assert "timestamp_utc" in manifest
    assert manifest["config"]["run"]["name"] == "smoke"
    # git hash is None when not in a repo; that's acceptable for the smoke test.
    assert "git_commit" in manifest

    out = tmp_path / "manifest.json"
    provenance.write_manifest(manifest, out)
    reloaded = json.loads(out.read_text(encoding="utf-8"))
    assert reloaded["config"]["run"]["name"] == "smoke"
