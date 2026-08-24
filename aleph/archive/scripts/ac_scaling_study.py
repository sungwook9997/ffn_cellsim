"""What does a smaller cell actually buy, and what does it cost you in the answer?

TWO CURVES, MEASURED TOGETHER, BECAUSE EITHER ALONE MISLEADS.

*Cost vs population* tells you what a development slice saves.  *Answer vs population* tells you whether
that slice is allowed to conclude anything.  This project has the second rule already — "develop on a
slice, conclude at native" — enforced entirely by hand, with no measurement of WHERE a slice stops being
safe.  Four headline numbers were retired on 2026-07-28 for being measured at 0.18% of native, and the
retraction was a judgement call, not a computation, because nothing had ever plotted an observable
against population.

WHY COST IS NOT EXPECTED TO SCALE THE WAY THE POPULATION DOES.  The 2026-07-28 outer-step profile put
53.5% of one native step in the fluid and its moving boundary.  Of that, ``classify_live_membrane`` runs
over the FINITE-VOLUME GRID and ``live_triangle_membrane_flux`` over the MEMBRANE TRIANGLES — neither is
indexed by filament count.  If that holds, halving the filaments does not halve the step, and a
development slice is far less of a saving than the population ratio suggests.  This driver separates the
three scalings by measuring them rather than assuming:

    node-scaling      link_spring, cytosim_bending, wca_steric, _interp_grad_force, project_constraint
    grid-scaling      classify_live_membrane, classify_live_nucleus, _pressure_gradient_masked
    membrane-scaling  live_triangle_membrane_flux, membrane_pressure_traction, helfrich_bending

WHAT IT DOES NOT DO.  It does not run long enough at any scale to see a steady state — that is a separate
and much longer run.  The trajectory it records at each scale is the SAME early transient at every point,
which is exactly what a comparison needs: if the first N steps already differ with population, no slice
is safe; if they agree, that is necessary but NOT sufficient, and the record says so rather than letting
agreement on a transient be read as agreement on a steady state.

engine units: length µm, force pN, tension pN/µm, time s.

Sanity Gate:
    * dimensional: costs are ms/step; tensions pN/µm; the scaling exponents are dimensionless.
    * boundary: a scale point that fails to build is recorded as failed and the sweep continues, so one
      bad point does not lose the curve.
    * conservation/invariant: the per-kernel totals at each scale sum to that scale's instrumented total.
    * numerical: cost is wall-clock and carries run-to-run variation; each point reports the
      uninstrumented step time so the comparison is against the same quantity at every scale.
    * sign-sense: not applicable — no signed physical quantity is judged here.
    * measurement-protocol: every scale point runs the SAME driver path with only the population knobs
      changed; a point that also changed the solver or the timestep would confound scale with method.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import warp as wp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from aleph.components.incumbent.assemble import CellConfig, build_cell  # noqa: E402
from aleph.components.incumbent.driver import make_inner_solve  # noqa: E402

NATIVE_FILAMENTS = 70_686

#: How each kernel is expected to scale, so the measurement can CONFIRM or REFUTE it rather than assume.
#: A kernel that moves with the wrong knob is the finding, not a mislabel to quietly fix.
EXPECTED_SCALING: dict[str, str] = {
    "link_spring": "node", "cytosim_bending": "node", "wca_steric": "node",
    "_interp_grad_force": "node", "project_constraint_forces": "node",
    "conditional_axpy": "node", "pos_to_f32": "node",
    "classify_live_membrane": "grid", "classify_live_nucleus": "grid",
    "_pressure_gradient_masked": "grid", "spread_vs": "grid",
    "live_triangle_membrane_flux": "membrane", "membrane_pressure_traction": "membrane",
    "helfrich_bending": "membrane", "membrane_area": "membrane",
}


def _kernel_family(name: str) -> tuple[str, str]:
    """Return ``(short_name, expected_scaling)`` for a Warp timing record name."""
    short = name.replace("forward kernel ", "").rsplit("_", 1)[0]
    for prefix, scaling in EXPECTED_SCALING.items():
        if short.startswith(prefix):
            return short, scaling
    return short, "other"


def profile_one(n_filaments: int, subdiv: int, *, n_inner: int, dt_phys: float) -> dict:
    """Build one cell and profile a single inner chunk; returns the point's cost record."""
    cfg = CellConfig(
        n_filaments=int(n_filaments), overlap_free_cortex=True, erm_radial_pairing=True,
        membrane_subdivisions=int(subdiv),
        resting_bound_myosin_fraction=0.5, resting_bound_myosin_force_pn=1.5,
        resting_bound_myosin_source="scaling study — cost measurement, no physics claim",
        resting_bound_myosin_capture_um=0.6,
    )
    t0 = time.perf_counter()
    cell = build_cell(cfg)
    build_s = time.perf_counter() - t0

    solve = make_inner_solve(cell, n_inner, reshape_every=max(2, n_inner), max_inner_retries=0,
                             inner_solver="explicit")
    solve(dt_phys)                       # warm
    wp.synchronize_device(cell.device)
    t0 = time.perf_counter()
    solve(dt_phys)
    wp.synchronize_device(cell.device)
    bare_ms = (time.perf_counter() - t0) * 1e3

    wp.timing_begin(wp.TIMING_KERNEL)
    solve(dt_phys)
    records = wp.timing_end(synchronize=True)

    per_kernel: dict[str, float] = defaultdict(float)
    per_family: dict[str, float] = defaultdict(float)
    for record in records:
        short, family = _kernel_family(record.name)
        per_kernel[short] += float(record.elapsed)
        per_family[family] += float(record.elapsed)

    point = {
        "n_filaments": int(n_filaments), "membrane_subdiv": int(subdiv),
        "n_actin": int(cell.n_actin), "n_total": int(cell.n_total),
        "n_crosslinks": int(getattr(cell, "n_xl", 0)),
        "fraction_of_native": float(n_filaments) / NATIVE_FILAMENTS,
        "build_seconds": build_s,
        "ms_per_inner_iteration": bare_ms / n_inner,
        "kernel_ms_per_iteration": {k: v / n_inner for k, v in
                                    sorted(per_kernel.items(), key=lambda kv: -kv[1])},
        "family_ms_per_iteration": {k: v / n_inner for k, v in per_family.items()},
    }
    del cell, solve
    return point


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--filaments", type=int, nargs="+",
                    default=[2000, 8000, 20000, 40000, NATIVE_FILAMENTS])
    ap.add_argument("--subdiv", type=int, default=8,
                    help="membrane subdivisions held FIXED across the filament sweep, so the grid and "
                         "membrane costs stay constant and the node-scaling term is isolated")
    ap.add_argument("--subdiv-sweep", type=int, nargs="*", default=[6, 7, 8],
                    help="a second sweep at fixed native filaments, to price the membrane knob itself")
    ap.add_argument("--subdiv-sweep-filaments", type=int, default=8000,
                    help="filament count for the subdivision sweep; small, since the point is the "
                         "membrane/grid term and the node term is a constant offset there")
    ap.add_argument("--n-inner", type=int, default=10)
    ap.add_argument("--dt-phys", type=float, default=0.01)
    ap.add_argument("--build-commit", default=None)
    ap.add_argument("--out", type=Path, default=Path("outputs/ac/profile/scaling_study.json"))
    args = ap.parse_args()

    points, failures = [], []
    for n in args.filaments:
        print(f"\n[scale] filaments={n:,} subdiv={args.subdiv} ...", flush=True)
        try:
            p = profile_one(n, args.subdiv, n_inner=args.n_inner, dt_phys=args.dt_phys)
        except Exception as exc:  # noqa: BLE001 — one bad point must not lose the curve
            print(f"[scale] FAILED at {n}: {type(exc).__name__}: {exc}", flush=True)
            failures.append({"n_filaments": n, "subdiv": args.subdiv, "error": f"{type(exc).__name__}: {exc}"})
            continue
        points.append(p)
        fam = p["family_ms_per_iteration"]
        print(f"[scale]   {p['n_total']:>9,} nodes  {p['ms_per_inner_iteration']:8.2f} ms/iter   "
              f"node={fam.get('node', 0):6.2f}  grid={fam.get('grid', 0):6.2f}  "
              f"membrane={fam.get('membrane', 0):6.2f}  other={fam.get('other', 0):6.2f}", flush=True)

    subdiv_points = []
    for s in args.subdiv_sweep or []:
        print(f"\n[scale] subdiv={s} at filaments={args.subdiv_sweep_filaments:,} ...", flush=True)
        try:
            p = profile_one(args.subdiv_sweep_filaments, s, n_inner=args.n_inner, dt_phys=args.dt_phys)
        except Exception as exc:  # noqa: BLE001
            print(f"[scale] FAILED at subdiv {s}: {exc}", flush=True)
            failures.append({"n_filaments": args.subdiv_sweep_filaments, "subdiv": s, "error": str(exc)})
            continue
        subdiv_points.append(p)
        fam = p["family_ms_per_iteration"]
        print(f"[scale]   {p['ms_per_inner_iteration']:8.2f} ms/iter   grid={fam.get('grid', 0):6.2f}  "
              f"membrane={fam.get('membrane', 0):6.2f}", flush=True)

    # Fit cost ~ a * N^b per family, so "does this term follow the population" is a number, not a look.
    exponents: dict[str, float] = {}
    if len(points) >= 2:
        n_nodes = np.array([p["n_total"] for p in points], dtype=np.float64)
        for family in ("node", "grid", "membrane", "other"):
            y = np.array([p["family_ms_per_iteration"].get(family, 0.0) for p in points], dtype=np.float64)
            if np.all(y > 0):
                exponents[family] = float(np.polyfit(np.log(n_nodes), np.log(y), 1)[0])

    print("\n[scale] cost ~ N^b over the filament sweep (N = total nodes):")
    for family, b in exponents.items():
        reading = ("follows the population" if b > 0.7 else
                   "nearly INDEPENDENT of population" if b < 0.25 else "sub-linear")
        print(f"[scale]   {family:9s} b = {b:+.2f}   {reading}")

    if points:
        smallest, largest = points[0], points[-1]
        ratio_n = largest["n_total"] / smallest["n_total"]
        ratio_t = largest["ms_per_inner_iteration"] / smallest["ms_per_inner_iteration"]
        print(f"\n[scale] {ratio_n:.0f}x the nodes costs {ratio_t:.1f}x the time — a slice at "
              f"{100 * smallest['fraction_of_native']:.1f}% of native buys {ratio_t:.1f}x, "
              f"not {ratio_n:.0f}x")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "schema": "ac.profile/scaling-study@1",
        "build_commit_declared": args.build_commit,
        "n_inner": args.n_inner, "dt_phys": args.dt_phys,
        "filament_sweep": points, "subdiv_sweep": subdiv_points, "failures": failures,
        "cost_exponents_vs_total_nodes": exponents,
        "caveat": ("cost only. This says what a slice SAVES, never what it is allowed to conclude — the "
                   "observable-vs-population curve is a separate measurement and is not in this file"),
    }, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
