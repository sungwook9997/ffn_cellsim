"""Taichi GPU initialization test — separate file so the rest of the smoke suite
can run on machines without a Taichi install.

Skipped if Taichi is not importable.
"""

from __future__ import annotations

import os

import pytest

ti = pytest.importorskip("taichi")


def test_taichi_init_via_acs(monkeypatch):
    """Importing the dispatcher and calling init_taichi should select an arch."""
    from acs.gpu import init_taichi

    # Honor env var if already set (real GPU runs); otherwise let Taichi pick.
    backend = os.environ.get("ACS_GPU_BACKEND") or "auto"
    arch = init_taichi(backend=backend)
    assert arch is not None
    # Ensure the chosen arch is one of the legal Taichi archs.
    assert str(arch) in {"Arch.cuda", "Arch.vulkan", "Arch.opengl",
                         "Arch.metal", "Arch.x64", "Arch.arm64", "Arch.cpu"}
    ti.reset()  # clean slate so other Taichi tests don't inherit state
