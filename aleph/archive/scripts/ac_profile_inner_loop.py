"""Where do the 192 ms of one native inner iteration actually go?

WHY THIS DID NOT EXIST.  The repository profiles exactly one thing — ``wp.TIMING_MEMCPY``, to prove the
GPU-residency gate is non-vacuous.  Nothing has ever profiled ``wp.TIMING_KERNEL``, so every statement
about where native time goes has been an inference from wall-clock divided by a loop count.  Measured
tonight: the resting solve costs 181.8 ms per inner iteration and the cortex-motor step costs 7.7 s for
40 inner iterations = 192 ms each, which agrees — but agreement on the TOTAL says nothing about which
kernel owns it, and optimising the wrong one is free of any signal that it was wrong.

WHAT IT MEASURES.  One native build, then N inner force-assembly + descent iterations with kernel timing
active, aggregated per kernel: total ms, launch count, ms per launch, and share of the loop.  The build
is excluded — it is paid once per run and is not what a longer trajectory spends its time on.

WHY THIS IS NOT A BENCHMARK, AND SAYS SO.  Kernel timing serialises the stream: Warp inserts events
around every launch, so the SUM here exceeds the un-instrumented wall-clock and the per-kernel SHARES are
what transfer, not the absolute totals.  The record carries both, and the uninstrumented wall-clock of
the same loop is measured separately in the same run so the instrumentation overhead is visible rather
than assumed.

Sanity Gate:
    * dimensional: every reported time is milliseconds; shares are dimensionless and sum to 1 over the
      profiled window.
    * boundary: a kernel that never launches does not appear, rather than appearing at zero — an absent
      row means "not exercised in this configuration", which is a different statement from "free".
    * conservation/invariant: the per-kernel totals sum to the instrumented loop total by construction;
      that sum is reported next to the uninstrumented wall-clock so the overhead is explicit.
    * numerical: no physics is judged here and no tolerance is applied; this module computes a cost
      profile, not a result.
    * measurement-protocol: profiling runs the SAME inner path the production driver runs, at the SAME
      native population — a profile taken on a slice would rank kernels by a different work distribution
      and is the specific mistake this file exists to avoid.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import warp as wp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from aleph.components.incumbent.assemble import CellConfig, build_cell  # noqa: E402
from aleph.components.incumbent.driver import _accumulate_all, make_inner_solve  # noqa: E402

NATIVE_FILAMENTS = 70_686


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--filaments", type=int, default=NATIVE_FILAMENTS)
    ap.add_argument("--iterations", type=int, default=20,
                    help="inner iterations to profile; the build is excluded from the profile")
    ap.add_argument("--myosin-fraction", type=float, default=0.5,
                    help="resting bound-myosin fraction — the crossbridge family is a live cost and "
                         "profiling without it would rank the passive kernels against a cell nobody runs")
    ap.add_argument("--inner-solver", default="explicit",
                    help="the inner path to profile; 'explicit' is the production default and "
                         "'erm_jacobi_pure' is what the resting driver runs")
    ap.add_argument("--dt-phys", type=float, default=0.01)
    ap.add_argument("--stage", choices=("inner", "outer"), default="inner",
                    help="inner: the mechanical chunk only. outer: ONE full physical step through the "
                         "PhysicalScheduler — fluid subcycles, kinetics, reshape and commit included. "
                         "Measured 2026-07-28: 40 inner iterations account for 1.19 s of a 7.7 s "
                         # `%%`: argparse runs help through `% params`, so a bare `% o` is a conversion
                         # specifier wanting an int. See tests/scripts/test_argparse_help_is_formattable.
                         "production step, so 85%% of the cost is in `outer` and had never been profiled")
    ap.add_argument("--out", type=Path,
                    default=Path("outputs/ac/profile/inner_loop_profile.json"))
    args = ap.parse_args()

    cfg = CellConfig(
        n_filaments=int(args.filaments), overlap_free_cortex=True, erm_radial_pairing=True,
        membrane_subdivisions=6,
        resting_bound_myosin_fraction=float(args.myosin_fraction),
        resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="profile only — cost measurement, no physics claim",
        resting_bound_myosin_capture_um=0.6,
    )
    print(f"[profile] building native cell: {args.filaments:,} filaments ...", flush=True)
    t0 = time.perf_counter()
    cell = build_cell(cfg)
    build_s = time.perf_counter() - t0
    print(f"[profile] built in {build_s:.1f} s — {cell.n_actin:,} actin / {cell.n_total:,} total nodes",
          flush=True)

    n_iter = int(args.iterations)

    # PROFILE THE WHOLE INNER ITERATION, not just the force assembly. The first version of this file
    # profiled `_accumulate_all` alone and measured 24.6 ms against a production iteration of 192 ms —
    # so assembly is 13% and 87% was left unattributed, which is exactly the guess this module was
    # written to stop. `inner_solve(dt_phys)` runs the real chunk: position copy, assembly, force copy,
    # the inextensibility projection over every fiber, and the descent step.
    inner_solve = make_inner_solve(
        cell, n_iter, reshape_every=max(2, n_iter), max_inner_retries=0,
        inner_solver=str(args.inner_solver),
    )

    if args.stage == "outer":
        from aleph.components.fluid.scheduler import PhysicalScheduler
        from aleph.components.incumbent.driver import PI_0_PA
        sched = PhysicalScheduler(
            substrate=cell.substrate, domain=cell.domain, membrane_bc=cell.membrane_bc,
            inner_solve=inner_solve, osmotic_difference=lambda _t: PI_0_PA,
        )
        divisor = 1  # per STEP, not per inner iteration
        run = lambda: sched.outer_step(float(args.dt_phys), biot_cfl_safety=0.9)  # noqa: E731
    else:
        divisor = n_iter
        run = lambda: inner_solve(float(args.dt_phys))  # noqa: E731

    # Uninstrumented wall-clock FIRST, so instrumentation overhead is measured rather than asserted.
    # One warm pass so module load/JIT is not counted.
    run()
    wp.synchronize_device(cell.device)
    t0 = time.perf_counter()
    run()
    wp.synchronize_device(cell.device)
    bare_ms = (time.perf_counter() - t0) * 1e3

    wp.timing_begin(wp.TIMING_KERNEL)
    run()
    records = wp.timing_end(synchronize=True)
    n_iter = divisor

    totals: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    for record in records:
        totals[record.name] += float(record.elapsed)
        counts[record.name] += 1
    instrumented_ms = sum(totals.values())

    rows = sorted(totals.items(), key=lambda kv: kv[1], reverse=True)
    print(f"\n[profile] ONE inner_solve chunk of {n_iter} iterations, solver={args.inner_solver!r}, at {cell.n_total:,} nodes")
    print(f"[profile] uninstrumented {bare_ms / n_iter:8.2f} ms/iteration   "
          f"instrumented {instrumented_ms / n_iter:8.2f} ms/iteration "
          f"(overhead {instrumented_ms / bare_ms:.2f}x — shares transfer, absolutes do not)\n")
    print(f"{'kernel':52s} {'ms/iter':>9s} {'share':>7s} {'launches':>9s} {'ms/launch':>10s}")
    for name, total in rows[:18]:
        print(f"{name[:52]:52s} {total / n_iter:9.3f} {100 * total / instrumented_ms:6.1f}% "
              f"{counts[name] // n_iter:9d} {total / counts[name]:10.4f}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "schema": "ac.profile/inner-loop@1",
        "census": {"n_filaments": int(args.filaments), "n_actin": int(cell.n_actin),
                   "n_total": int(cell.n_total),
                   "n_crosslinks": int(getattr(cell, "n_xl", 0)),
                   "fraction_of_native": float(args.filaments) / NATIVE_FILAMENTS},
        "device": str(cell.device),
        "build_seconds": build_s,
        "iterations": n_iter,
        "inner_solver": str(args.inner_solver),
        "uninstrumented_ms_per_iteration": bare_ms / n_iter,
        "instrumented_ms_per_iteration": instrumented_ms / n_iter,
        "instrumentation_overhead_factor": instrumented_ms / bare_ms if bare_ms else None,
        "note": ("kernel timing serialises the stream, so per-kernel SHARES are the transferable "
                 "quantity and the absolute totals are inflated; the uninstrumented figure above is the "
                 "one to compare against production wall-clock"),
        "kernels": [{"name": n, "ms_per_iteration": t / n_iter,
                     "share": t / instrumented_ms,
                     "launches_per_iteration": counts[n] // n_iter,
                     "ms_per_launch": t / counts[n]} for n, t in rows],
    }, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
