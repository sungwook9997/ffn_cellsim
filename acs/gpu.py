"""GPU backend portability layer.

Reads `ACS_GPU_BACKEND` to pick the Taichi arch at runtime so the same code runs
on the primary RTX A5000 (cuda), the lab's RTX 4090×2 (cuda), Vulkan-only hosts,
or CPU-only fallback (Colab T4-less, CI).

Allowed values (case-insensitive):
    cuda   — NVIDIA via CUDA (default; required for production)
    vulkan — cross-vendor GPU
    opengl — last-resort GPU
    cpu    — fallback / CI
    auto   — let Taichi pick (`ti.gpu` then `ti.cpu`)

Usage:
    from acs.gpu import init_taichi
    arch = init_taichi()
    print("Taichi running on", arch)
"""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)

_BACKEND_MAP = {
    "cuda": "cuda",
    "vulkan": "vulkan",
    "opengl": "opengl",
    "metal": "metal",
    "cpu": "cpu",
    "auto": "auto",
}


def _resolve_backend(explicit: str | None) -> str:
    raw = (explicit or os.environ.get("ACS_GPU_BACKEND") or "auto").strip().lower()
    if raw not in _BACKEND_MAP:
        raise ValueError(
            f"ACS_GPU_BACKEND={raw!r} not recognized. "
            f"Choose from: {sorted(_BACKEND_MAP)}"
        )
    return _BACKEND_MAP[raw]


def init_taichi(backend: str | None = None, **ti_kwargs) -> str:
    """Initialize Taichi with the requested arch and return the chosen backend name.

    Any extra kwargs are forwarded to ``ti.init`` so callers can pass
    ``device_memory_GB``, ``random_seed``, ``debug``, etc.
    """
    import taichi as ti  # local import to keep module importable without Taichi

    name = _resolve_backend(backend)
    arch_lookup = {
        "cuda": ti.cuda,
        "vulkan": ti.vulkan,
        "opengl": ti.opengl,
        "metal": ti.metal,
        "cpu": ti.cpu,
        "auto": ti.gpu,
    }
    arch = arch_lookup[name]
    ti.init(arch=arch, **ti_kwargs)
    chosen = ti.cfg.arch
    logger.info("Taichi initialized: requested=%s resolved=%s", name, chosen)
    return str(chosen)
