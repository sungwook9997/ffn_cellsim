"""Stage 0 verification — corresponds to docs/10_dev_roadmap.md "Validation".

Run after env install completes:
    python scripts/verify_env.py
or with a forced backend:
    ACS_GPU_BACKEND=vulkan python scripts/verify_env.py

Reports Taichi arch, CUDA visibility (when applicable), and whether the
expected core libraries import cleanly.
"""

from __future__ import annotations

import importlib
import sys

REQUIRED = ["taichi", "numpy", "scipy", "h5py", "yaml", "matplotlib", "plotly"]
OPTIONAL = ["vispy"]   # vispy needs OpenGL libs that aren't always present headlessly


def _try_import(name: str) -> tuple[bool, str]:
    try:
        mod = importlib.import_module(name)
        return True, getattr(mod, "__version__", "unknown")
    except Exception as exc:  # noqa: BLE001 — surface anything that breaks the import
        return False, f"{type(exc).__name__}: {exc}"


def main() -> int:
    print("=== ActiveCellSim Stage 0 environment check ===\n")
    print(f"Python:  {sys.version.split()[0]}")
    print(f"Sys exe: {sys.executable}\n")

    failures: list[str] = []

    print("Required imports:")
    for name in REQUIRED:
        ok, info = _try_import(name)
        marker = "OK " if ok else "FAIL"
        print(f"  [{marker}] {name:<10} {info}")
        if not ok:
            failures.append(name)

    print("\nOptional imports:")
    for name in OPTIONAL:
        ok, info = _try_import(name)
        marker = "OK " if ok else "WARN"
        print(f"  [{marker}] {name:<10} {info}")

    print("\nTaichi GPU dispatch:")
    try:
        from acs.gpu import init_taichi

        arch = init_taichi()
        print(f"  arch = {arch}")
    except Exception as exc:  # noqa: BLE001
        print(f"  FAIL: {type(exc).__name__}: {exc}")
        failures.append("taichi-init")

    if failures:
        print(f"\nStage 0 verification FAILED on: {failures}")
        return 1
    print("\nStage 0 verification PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
