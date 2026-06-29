#!/usr/bin/env python
"""Production preflight for GPU-main physiological ffn_cellsim runs.

This is intentionally a launch-time guard, not documentation.  A production
run should fail here before it can burn hours on the wrong machine or with a
non-physiological full-cell baseline.
"""

from __future__ import annotations

import argparse
import importlib
import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from typing import Any

import hoomd

from ffn_sim.common.production_policy import (
    DEFAULT_MCF7_TURGOR_PA,
    add_production_device_args,
    require_full_cell_physiological_baseline,
    require_production_device,
)
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.scripts.mcf7_fullcell_stage1 import _cfg_for, _resolve_compartments


@dataclass(slots=True)
class PreflightCheck:
    name: str
    passed: bool
    detail: str


def _check(name: str, passed: bool, detail: str) -> PreflightCheck:
    return PreflightCheck(name=name, passed=bool(passed), detail=detail)


def _check_device(device: str, allow_cpu_dev: bool) -> PreflightCheck:
    if device != "gpu":
        return _check(
            "device",
            False,
            f"--device {device!r} is not a production device; production must use gpu.",
        )
    try:
        require_production_device(
            device, allow_cpu_dev=allow_cpu_dev, hoomd_module=hoomd
        )
    except ValueError as exc:
        return _check("device", False, str(exc))
    return _check("device", True, "--device gpu requested and HOOMD accepts it")


def _check_hoomd_gpu_build() -> PreflightCheck:
    gpu_enabled = bool(getattr(hoomd.version, "gpu_enabled", False))
    detail = f"HOOMD {hoomd.version.version}, gpu_enabled={gpu_enabled}"
    return _check("hoomd_gpu_build", gpu_enabled, detail)


def _check_cupy() -> PreflightCheck:
    try:
        cupy = importlib.import_module("cupy")
    except Exception as exc:  # noqa: BLE001
        return _check("cupy", False, f"cupy import failed: {exc}")
    return _check("cupy", True, f"cupy {getattr(cupy, '__version__', 'unknown')}")


def _check_nvidia_smi() -> PreflightCheck:
    exe = shutil.which("nvidia-smi")
    if exe is None:
        return _check("nvidia_smi", False, "nvidia-smi not found on PATH")
    try:
        proc = subprocess.run(
            [
                exe,
                "--query-gpu=name,driver_version,memory.total",
                "--format=csv,noheader",
            ],
            check=True,
            text=True,
            capture_output=True,
            timeout=10,
        )
    except Exception as exc:  # noqa: BLE001
        return _check("nvidia_smi", False, f"nvidia-smi failed: {exc}")
    gpu_line = proc.stdout.strip().splitlines()[0] if proc.stdout.strip() else exe
    return _check("nvidia_smi", True, gpu_line)


def _check_mcf7_baseline(turgor_pa: float | None) -> PreflightCheck:
    try:
        cfg = _cfg_for(
            n_fil=200,
            n_motors=20,
            n_xl=200,
            stepping_mode="grip_walk",
            force_scaling=True,
            backbone_nm=700.0,
        )
        p_cortex = resolve_h3_derived(cfg)
        compartments = _resolve_compartments(p_cortex, turgor_pa=turgor_pa)
        require_full_cell_physiological_baseline(compartments)
    except Exception as exc:  # noqa: BLE001
        return _check("mcf7_fullcell_baseline", False, str(exc))

    cyto = compartments["p_cytoplasm"]
    ev = compartments["p_enclosed_volume"]
    return _check(
        "mcf7_fullcell_baseline",
        True,
        (
            f"cytoplasm eta={cyto.eta_eff:g} Pa*s, "
            f"cell_type={cyto.cell_type}; turgor_dP0={ev.turgor_dP0:g} Pa"
        ),
    )


def run_preflight(
    *,
    device: str = "gpu",
    allow_cpu_dev: bool = False,
    turgor_pa: float | None = DEFAULT_MCF7_TURGOR_PA,
) -> tuple[bool, list[PreflightCheck]]:
    """Run production preflight checks.

    Returns:
        ``(passed, checks)``.  ``passed`` is true only when every check passes.
    """
    checks = [
        _check_device(device, allow_cpu_dev),
        _check_hoomd_gpu_build(),
        _check_cupy(),
        _check_nvidia_smi(),
        _check_mcf7_baseline(turgor_pa),
    ]
    return all(c.passed for c in checks), checks


def _print_human(checks: list[PreflightCheck]) -> None:
    print("=== ffn_cellsim production preflight ===", flush=True)
    for check in checks:
        mark = "PASS" if check.passed else "FAIL"
        print(f"[{mark:4s}] {check.name}: {check.detail}", flush=True)
    passed = all(c.passed for c in checks)
    print(
        "\nRESULT " + ("production-ready" if passed else "NOT production-ready"),
        flush=True,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    add_production_device_args(ap, default="gpu")
    ap.add_argument(
        "--turgor-pa",
        type=float,
        default=DEFAULT_MCF7_TURGOR_PA,
        help="MCF7 enclosed-volume baseline osmotic turgor to validate [Pa]",
    )
    ap.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    args = ap.parse_args()

    passed, checks = run_preflight(
        device=args.device,
        allow_cpu_dev=args.allow_cpu_dev,
        turgor_pa=args.turgor_pa,
    )
    if args.json:
        payload: dict[str, Any] = {
            "passed": passed,
            "checks": [asdict(c) for c in checks],
        }
        print(json.dumps(payload, indent=2), flush=True)
    else:
        _print_human(checks)
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
