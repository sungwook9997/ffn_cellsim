"""Production-run guardrails for GPU-main physiological simulations.

These helpers turn project policy into executable checks.  The point is to make
invalid production states fail at launch instead of silently falling back to
CPU, water drag, or an unpressurised shell.
"""

from __future__ import annotations

import argparse
import importlib
import math
from typing import Any


DEFAULT_MCF7_TURGOR_PA: float = 40.0


def hoomd_has_gpu_build(hoomd_module: Any) -> bool:
    """True iff the imported HOOMD module is a GPU-enabled build."""
    return bool(getattr(hoomd_module.version, "gpu_enabled", False))


def require_production_device(
    device: str,
    *,
    allow_cpu_dev: bool,
    hoomd_module: Any,
) -> None:
    """Enforce GPU-main production device policy.

    CPU is still useful for local smoke/dev work, but only when explicitly
    requested with a dev escape hatch.  GPU requests also require a GPU-enabled
    HOOMD build so a CPU-only conda env cannot masquerade as production.
    """
    if device == "cpu" and not allow_cpu_dev:
        raise ValueError(
            "--device cpu is dev-only. Use --device gpu for production, or pass "
            "--allow-cpu-dev for an explicit local smoke/dev run."
        )
    if device == "gpu" and not hoomd_has_gpu_build(hoomd_module):
        raise ValueError(
            "Requested --device gpu, but this HOOMD build reports gpu_enabled=False. "
            "Install/activate the CUDA GPU environment before starting production."
        )


def add_production_device_args(
    parser: argparse.ArgumentParser,
    *,
    default: str = "gpu",
) -> None:
    """Add the standard GPU-main device flags to an argparse parser."""
    parser.add_argument("--device", choices=["cpu", "gpu"], default=default)
    parser.add_argument(
        "--allow-cpu-dev",
        action="store_true",
        help="explicitly allow --device cpu for local smoke/dev runs; production is GPU",
    )


def validate_production_device_args(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    *,
    hoomd_module: Any | None = None,
) -> None:
    """Validate standard device args and report failures as argparse errors."""
    if hoomd_module is None:
        hoomd_module = importlib.import_module("hoomd")
    try:
        require_production_device(
            args.device,
            allow_cpu_dev=args.allow_cpu_dev,
            hoomd_module=hoomd_module,
        )
    except ValueError as exc:
        parser.error(str(exc))


def require_physiological_cytoplasm(
    p_cytoplasm: Any,
    *,
    allow_water_dev: bool = False,
) -> None:
    """Require a resolved non-water cytoplasm viscosity for production."""
    if p_cytoplasm is None:
        if allow_water_dev:
            return
        raise ValueError(
            "Production full-cell runs require p_cytoplasm with measured "
            "cytoplasm viscosity; p_cytoplasm=None would leave water/no-op drag."
        )
    eta_eff = float(getattr(p_cytoplasm, "eta_eff"))
    eta_water = float(getattr(p_cytoplasm, "eta_water"))
    cell_type = getattr(p_cytoplasm, "cell_type", None)
    if not (math.isfinite(eta_eff) and eta_eff > 0.0):
        raise ValueError(f"p_cytoplasm.eta_eff must be finite and > 0; got {eta_eff!r}.")
    if (
        eta_eff == eta_water
        or bool(getattr(p_cytoplasm, "is_noop", False))
        or cell_type == "water"
    ):
        if allow_water_dev:
            return
        raise ValueError(
            "Production full-cell runs require measured cytoplasm viscosity, "
            f"not water/no-op drag (eta_eff={eta_eff:g}, cell_type={cell_type!r})."
        )


def require_physiological_turgor(
    p_enclosed_volume: Any,
    *,
    allow_unpressurized_dev: bool = False,
) -> None:
    """Require positive baseline osmotic turgor for production."""
    if p_enclosed_volume is None:
        if allow_unpressurized_dev:
            return
        raise ValueError(
            "Production full-cell runs require p_enclosed_volume with positive "
            "baseline osmotic turgor."
        )
    raw_turgor = getattr(p_enclosed_volume, "turgor_dP0")
    if raw_turgor is not None:
        turgor = float(raw_turgor)
        if math.isfinite(turgor) and turgor > 0.0:
            return
    if allow_unpressurized_dev:
        return
    raise ValueError(
        "Production full-cell runs require baseline osmotic turgor "
        f"(turgor_dP0 > 0; default {DEFAULT_MCF7_TURGOR_PA:g} Pa). "
        "Zero/None pressure is dev/attribution only."
    )


def require_full_cell_physiological_baseline(
    compartments: dict[str, Any],
    *,
    allow_water_dev: bool = False,
    allow_unpressurized_dev: bool = False,
) -> None:
    """Require the two recurring full-cell physiological baseline contracts."""
    require_physiological_cytoplasm(
        compartments.get("p_cytoplasm"),
        allow_water_dev=allow_water_dev,
    )
    require_physiological_turgor(
        compartments.get("p_enclosed_volume"),
        allow_unpressurized_dev=allow_unpressurized_dev,
    )
