#!/usr/bin/env python
"""Single entrypoint for production ffn_cellsim runs.

The target simulations still live in their existing modules.  This launcher is
the narrow production gate: run preflight first, then call the selected driver
with standard GPU-main flags.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass

from ffn_sim.common.production_policy import (
    DEFAULT_MCF7_TURGOR_PA,
    add_production_device_args,
)
from ffn_sim.scripts.production_preflight import run_preflight


@dataclass(frozen=True, slots=True)
class ProductionTarget:
    module: str
    description: str


TARGETS = {
    "mcf7-fullcell": ProductionTarget(
        module="ffn_sim.scripts.mcf7_fullcell_stage1",
        description="MCF7 full-cell Stage-1 driver",
    ),
}


def _has_option(argv: list[str], *names: str) -> bool:
    return any(
        any(arg == name or arg.startswith(name + "=") for name in names)
        for arg in argv
    )


def _inject_standard_args(
    driver_args: list[str],
    *,
    device: str,
    allow_cpu_dev: bool,
    turgor_pa: float,
) -> list[str]:
    out = list(driver_args)
    if not _has_option(out, "--device"):
        out.extend(["--device", device])
    if allow_cpu_dev and not _has_option(out, "--allow-cpu-dev"):
        out.append("--allow-cpu-dev")
    if not _has_option(out, "--turgor-pa"):
        out.extend(["--turgor-pa", str(turgor_pa)])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    add_production_device_args(ap, default="gpu")
    ap.add_argument(
        "--turgor-pa",
        type=float,
        default=DEFAULT_MCF7_TURGOR_PA,
        help="MCF7 enclosed-volume baseline osmotic turgor [Pa]",
    )
    ap.add_argument(
        "--skip-preflight",
        action="store_true",
        help="skip production preflight; intended only for tightly controlled dev runs",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print the target command after preflight without executing it",
    )
    ap.add_argument("target", choices=sorted(TARGETS))
    args, driver_args = ap.parse_known_args()

    if not args.skip_preflight:
        passed, checks = run_preflight(
            device=args.device,
            allow_cpu_dev=args.allow_cpu_dev,
            turgor_pa=args.turgor_pa,
        )
        if not passed:
            for check in checks:
                mark = "PASS" if check.passed else "FAIL"
                print(f"[{mark:4s}] {check.name}: {check.detail}", flush=True)
            return 2

    target = TARGETS[args.target]
    final_driver_args = _inject_standard_args(
        driver_args,
        device=args.device,
        allow_cpu_dev=args.allow_cpu_dev,
        turgor_pa=args.turgor_pa,
    )
    cmd = [sys.executable, "-m", target.module, *final_driver_args]
    print(f"[production] {target.description}", flush=True)
    print("[production] command: " + " ".join(cmd), flush=True)
    if args.dry_run:
        return 0
    return subprocess.run(cmd, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
