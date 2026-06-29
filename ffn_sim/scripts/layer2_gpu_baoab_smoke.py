"""gbook-ONLY GPU-vs-CPU physics smoke for the device BAOAB on a small catch-bond growth.

Companion to ``layer2_gpu_scaleup.py``: where that driver runs a production native-N sweep,
this is a fast, small SAFETY check that the GPU device path (``integrator/baoab_device.py`` →
``gpu_local_snapshot`` + cupy) produces PHYSICALLY SOUND, finite, CPU-consistent output before
a long sweep is launched. It is deliberately tiny (N0 ~ 100, a few hours of biological time)
so it runs in seconds-to-minutes and is cheap to re-run after any GPU-path change.

What it asserts (the device-branch contract):
  1. The GPU run produces FINITE positions (no NaN/Inf escaping the device path).
  2. The GPU connected-core A/A0 agrees with the CPU run within a documented seed-noise band.

Why a BAND, not bit-exact parity: the GPU BAOAB Action uses cupy's own RNG stream (see the
``baoab_device`` correctness contract — CPU path is bit-identical to the frozen Action, GPU
path is *statistical*). Different Brownian noise → different microscopic trajectory → A/A0
differs by seed noise even at the same ``seed``. The band below is the policy tolerance.

SKIP semantics: if cupy / a CUDA GPU is unavailable (the CPU dev box), this prints a clear
SKIP message and exits 0 — so it is safe to run anywhere, and only actually exercises the GPU
branch on gbook.

Usage
-----
    python -m ffn_sim.scripts.layer2_gpu_baoab_smoke           # N0=100, ~4 h biological time
    python -m ffn_sim.scripts.layer2_gpu_baoab_smoke <N0> <total_time_h> <seed>
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.archive.hoomd_legacy.spheroid.observables import effective_radius
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.archive.hoomd_legacy.spheroid.proliferation import run_growth_pooled

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"

# Seed-noise band for |dA/A0_core| between the CPU (frozen numpy BAOAB) and GPU (cupy-RNG
# BAOAB) runs of the SAME (N0, seed). NOT a physics constant — a numerical-policy tolerance:
# the GPU path draws an independent Brownian stream (baoab_device correctness contract), so
# A/A0 parity is statistical. 0.2 is set conservatively wider than the ~0.05-0.15 seed-noise
# scatter the scaleup `parity` sub-command documents, so this fast small-N smoke never flags a
# benign RNG-stream difference as a regression — it catches a BROKEN device path (wrong update,
# lost force, non-finite drift), not microscopic noise. A real device-path bug shifts A/A0 far
# beyond 0.2 (or makes it non-finite), which this still catches.
_CORE_PARITY_BAND = 0.2


def _gpu_available() -> bool:
    """True iff cupy imports AND a CUDA GPU is actually visible (else SKIP cleanly).

    Mirrors the deferred-import policy of ``constrained_baoab.array_backend``: cupy is a
    GPU-host-only dep, so a bare ``import cupy`` on the CPU dev box must not hard-fail. We also
    probe ``runtime.getDeviceCount`` because a cupy install with no visible device would import
    yet still have no GPU to run on.
    """
    try:
        import cupy as cp  # deferred: GPU hosts only
    except Exception as exc:  # ImportError, or a CUDA-less build raising on import
        print(f"[smoke] SKIP — cupy unavailable ({type(exc).__name__}: {exc}).")
        return False
    try:
        if cp.cuda.runtime.getDeviceCount() < 1:
            print("[smoke] SKIP — cupy present but no CUDA device visible.")
            return False
    except Exception as exc:
        print(f"[smoke] SKIP — CUDA runtime query failed ({type(exc).__name__}: {exc}).")
        return False
    return True


def _run(kind: str, n0: int, total_time: float, seed: int) -> dict:
    """One small catch-bond growth on the requested device → the run record dict.

    Mirrors ``layer2_gpu_scaleup.run_one`` but with smoke-sized arguments and the device
    constructed/asserted here (a ``hoomd.device.GPU`` that silently fell back to CPU would mis-
    label the run; we raise instead).
    """
    import hoomd
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved)
    if kind == "gpu":
        device = hoomd.device.GPU(notice_level=0)
        if not isinstance(device, hoomd.device.GPU):
            raise RuntimeError("requested GPU but no CUDA device was constructed.")
    else:
        device = hoomd.device.CPU(notice_level=0)
    max_cells = int(5 * n0)  # x4 growth headroom + margin (same rule as the scaleup pool)
    res = run_growth_pooled(
        resolved, prolif, n_cells_init=n0, total_time=total_time, seed=seed,
        cohesion="catch", cad=cad, max_cells=max_cells, device=device,
    )
    device_actual = "gpu" if isinstance(device, hoomd.device.GPU) else "cpu"
    return {
        "device": device_actual,
        "R0_um": float(effective_radius(res["a0"]) * 1e6),
        "aa0_core": float(res["area_core_over_a0"][-1]),
        "aa0_hull": float(res["area_over_a0"][-1]),
        "n_final": int(res["n_cells"][-1]),
        "pos_final": np.asarray(res["pos_final"], dtype=np.float64),
    }


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    n0 = int(args[0]) if len(args) > 0 else 100
    total_time_h = float(args[1]) if len(args) > 1 else 4.0
    seed = int(args[2]) if len(args) > 2 else 0
    total_time = total_time_h * 3600.0  # h -> s (biological time)

    print(f"[smoke] device BAOAB GPU-vs-CPU catch-bond growth — "
          f"N0={n0} total_time={total_time_h} h seed={seed}")

    if not _gpu_available():
        # Safe on the CPU dev box: no GPU branch to exercise, so report SKIP and pass.
        print("[smoke] PASS (SKIP) — GPU branch not exercised on this host.")
        return 0

    cpu = _run("cpu", n0, total_time, seed)
    gpu = _run("gpu", n0, total_time, seed)
    print(f"  CPU: A/A0_core={cpu['aa0_core']:.3f} hull={cpu['aa0_hull']:.3f} "
          f"R0={cpu['R0_um']:.1f} um n_final={cpu['n_final']}")
    print(f"  GPU: A/A0_core={gpu['aa0_core']:.3f} hull={gpu['aa0_hull']:.3f} "
          f"R0={gpu['R0_um']:.1f} um n_final={gpu['n_final']}")

    # check 1: GPU positions finite (the device path did not leak a NaN/Inf)
    finite_ok = bool(np.all(np.isfinite(gpu["pos_final"])))
    # check 2: core A/A0 agrees within the documented seed-noise band
    d_core = abs(cpu["aa0_core"] - gpu["aa0_core"])
    parity_ok = bool(np.isfinite(d_core) and d_core <= _CORE_PARITY_BAND)

    print(f"  GPU positions finite: {'PASS' if finite_ok else 'FAIL'}")
    print(f"  |dA/A0_core| = {d_core:.3f}  (band <= {_CORE_PARITY_BAND}): "
          f"{'PASS' if parity_ok else 'FAIL'}")

    if finite_ok and parity_ok:
        print("[smoke] PASS — device BAOAB GPU branch finite + CPU-consistent.")
        return 0
    print("[smoke] FAIL — device BAOAB GPU branch diverged from CPU "
          "(non-finite or A/A0 beyond the seed-noise band).")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
