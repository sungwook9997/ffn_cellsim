"""Per-compartment force/updater profiler + GPU-readiness report (H.7 platform).

A P3 diagnostic tool that answers, for the explicit compartment stack:
  * what is each compartment's *declared* Compartment Performance Contract
    (hot-path priority, per-step force?, cpu_local_snapshot?, GPU path now)?
  * for a *built* cell, how expensive is each attached ``md.force.Custom`` and
    each per-batch ``Updater`` in wall-clock per call (microbench of
    ``set_forces`` / ``act``)?
  * which GPU env flags are active and which P0 compartments still lack a
    device-resident path (the GPU-main debt)?

It is intentionally cheap by default (``--table`` reads the registry only, no
HOOMD). ``--build`` assembles a baseline cell and microbenches the live
operations. Use it as the microbench hook the registry's ``GpuReadiness`` points
P0/P1 compartments at.

Examples
--------
    # registry table only (no HOOMD, instant):
    python ffn_sim/scripts/compartment_force_profile.py --table

    # microbench the live force/updater stack on a small suspended cell:
    python ffn_sim/scripts/compartment_force_profile.py --build --n-fil 200 --iters 50

    # GPU paths (set the env flags first):
    FFN_GPU_DEVICE_COMPARTMENTS=1 python ffn_sim/scripts/compartment_force_profile.py --build
"""

from __future__ import annotations

import argparse
import os
import time
from typing import Any

from ffn_sim.cell.compartment_registry import REGISTRY, GpuPath, HotPathPriority

_GPU_ENV_FLAGS = ("FFN_GPU_DEVICE_BAOAB", "FFN_GPU_DEVICE_COMPARTMENTS")


# ---------------------------------------------------------------------------
# Registry table (no HOOMD)
# ---------------------------------------------------------------------------
def print_contract_table() -> None:
    """Print the declared Compartment Performance Contract for every compartment."""
    rows = []
    for spec in REGISTRY.all():
        pc = spec.performance_contract
        rows.append(
            (
                spec.name,
                spec.status.value,
                pc.hot_path_priority.value,
                "force" if pc.per_step_force else ("upd" if pc.per_batch_updater else "-"),
                "cpu_snap" if pc.uses_cpu_local_snapshot else "-",
                "kdtree" if pc.uses_broad_phase else "-",
                pc.gpu_path_now.value,
                "Y" if spec.gpu_readiness.device_resident_now else "N",
            )
        )
    hdr = ("compartment", "status", "P", "hot", "snapshot", "search", "gpu_now", "dev")
    widths = [max(len(str(r[i])) for r in (rows + [hdr])) for i in range(len(hdr))]

    def fmt(r: tuple) -> str:
        return "  ".join(str(c).ljust(widths[i]) for i, c in enumerate(r))

    print("\n=== Compartment Performance Contracts (declared) ===")
    print(fmt(hdr))
    print("  ".join("-" * w for w in widths))
    for r in rows:
        print(fmt(r))


def print_gpu_readiness() -> None:
    """Print GPU env-flag status + the P0/P1 device-residency gap."""
    print("\n=== GPU env flags ===")
    for flag in _GPU_ENV_FLAGS:
        print(f"  {flag} = {os.environ.get(flag, '(unset)')}")

    print("\n=== P0/P1 GPU-readiness gap ===")
    for spec in REGISTRY.all():
        pc = spec.performance_contract
        if pc.hot_path_priority not in (HotPathPriority.P0, HotPathPriority.P1):
            continue
        dr = spec.gpu_readiness.device_resident_now
        builtin = pc.gpu_path_now in (GpuPath.HOOMD_BUILTIN, GpuPath.CUPY, GpuPath.NATIVE_PLUGIN)
        status = "OK" if (dr or builtin) else "DEBT"
        debt = spec.gpu_readiness.optimization_debt or ""
        print(f"  [{pc.hot_path_priority.value}] {spec.name:<22} {status:<5} {debt}")


# ---------------------------------------------------------------------------
# Live microbench (requires HOOMD + a built cell)
# ---------------------------------------------------------------------------
def _bench_callable(fn, iters: int) -> float:
    """Mean wall-time per call (seconds) of ``fn(timestep)`` over ``iters`` calls."""
    # warmup
    try:
        fn(0)
    except Exception as exc:  # noqa: BLE001 — report, don't crash the whole profile
        return float("nan")
    t0 = time.perf_counter()
    for k in range(iters):
        fn(k)
    return (time.perf_counter() - t0) / max(iters, 1)


def microbench_live_cell(
    *,
    n_fil: int,
    iters: int,
    allow_no_nucleus: bool,
) -> None:
    """Build a small baseline cell and microbench each Custom force / updater.

    Times each attached ``md.force.Custom.set_forces`` and each Updater's
    ``act``/``set_forces`` directly (not a full integrated run), so the per-op
    cost is isolated. This is the per-compartment microbench the registry's
    GpuReadiness.microbench_hook points at.
    """
    import hoomd  # noqa: F401  (ensure HOOMD present)

    from ffn_sim.cell.manifest import build_baseline_cell, load_manifest

    manifest = load_manifest("mcf7_baseline.yaml")
    # Shrink the cortex for a cheap profile (physics-irrelevant for *timing*).
    cov = manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
    cov["n_filaments"] = int(n_fil)

    print(f"\n=== Building baseline cell (n_fil={n_fil}, allow_no_nucleus={allow_no_nucleus}) ===")
    t0 = time.perf_counter()
    cell = build_baseline_cell(
        manifest=manifest,
        allow_no_nucleus=allow_no_nucleus,
        with_baoab=True,
    )
    print(f"  build: {time.perf_counter() - t0:.2f} s")

    sim = getattr(cell, "sim", None) or getattr(cell, "simulation", None)
    if sim is None:
        # Cell may expose the integrator/operations differently; try common attrs.
        sim = getattr(cell, "_sim", None)
    if sim is None:
        print("  [warn] could not locate the hoomd.Simulation on the Cell object; "
              "microbench skipped. (Use --table for the declared contracts.)")
        return

    ig = sim.operations.integrator
    print(f"\n=== Per-Custom-force microbench ({iters} calls each) ===")
    forces = list(getattr(ig, "forces", []) or [])
    for f in forces:
        name = type(f).__name__
        fn = getattr(f, "set_forces", None)
        if fn is None:
            print(f"  {name:<34} (no set_forces)")
            continue
        dt = _bench_callable(fn, iters)
        print(f"  {name:<34} {dt * 1e6:9.1f} us/call")

    print(f"\n=== Per-Updater microbench ({iters} calls each) ===")
    updaters = list(getattr(sim.operations, "updaters", []) or [])
    for u in updaters:
        action = getattr(u, "action", u)
        name = type(action).__name__
        fn = getattr(action, "act", None) or getattr(action, "set_forces", None)
        if fn is None:
            print(f"  {name:<34} (no act)")
            continue
        dt = _bench_callable(fn, iters)
        print(f"  {name:<34} {dt * 1e6:9.1f} us/call")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", action="store_true", help="print declared contracts (no HOOMD)")
    p.add_argument("--build", action="store_true", help="build a small cell + microbench live ops")
    p.add_argument("--n-fil", type=int, default=200, help="cortex filament count for the microbench cell")
    p.add_argument("--iters", type=int, default=50, help="calls per op in the microbench")
    p.add_argument("--with-nucleus", action="store_true",
                   help="keep the nucleus (default uses the sanctioned no-nucleus build for speed)")
    args = p.parse_args(argv)

    if not (args.table or args.build):
        args.table = True  # default: cheap table

    if args.table:
        print_contract_table()
        print_gpu_readiness()

    if args.build:
        try:
            microbench_live_cell(
                n_fil=args.n_fil,
                iters=args.iters,
                allow_no_nucleus=not args.with_nucleus,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"\n[error] live microbench failed: {type(exc).__name__}: {exc}")
            print("  (The declared --table contracts are still valid. The live "
                  "microbench needs a working ffn_sim env + buildable baseline cell.)")
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
